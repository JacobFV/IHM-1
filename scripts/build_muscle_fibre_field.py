"""Per-element fibre direction field for the repaired canonical muscular volumes.

Run with the isolated libigl environment, not the physiological runtime:
data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/build_muscle_fibre_field.py --self-test
data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/build_muscle_fibre_field.py \
  --output data/derived/muscle-fibre-field-v1

Pipeline, each stage recorded per entity:
  1 register        one whole-body rigid fit of the 22 source segment mass centres onto their canonical
                    bone-envelope centres seeds a per-body trimmed ICP of the source bone display meshes
                    onto the canonical bone group named by mechanics.json, at one clamped uniform PCA
                    scale. Rigid plus uniform scale only; no anatomical landmark correspondence is
                    claimed, and the seed is what stops nearly symmetric clouds registering mirrored.
  2 correspond      explicit OpenSim-muscle to canonical-entity name table, then a geometric receipt
                    (path-point to mesh-surface distance). Muscles with no counterpart are declared.
  3 tetrahedralise  TetGen pYq1.414 on the repaired surface, the same flags the tet-readiness proof used
  4 ends            surface nodes ranked by arc-length position along the registered OpenSim path
                    (terminal segments extended); the extreme deciles become the two Dirichlet patches
  5 harmonic        linear-tet FEM Laplace solve, phi=0 on the origin patch, phi=1 on the insertion
                    patch, natural (zero-flux) elsewhere; per-tet fibre = normalised constant grad phi
  6 pennate         where OpenSim declares pennation_angle_at_optimal>0 every fibre is rotated by that
                    angle in the plane spanned by the local fibre and the muscle volume's thinnest
                    principal axis. A uniform unipennate prior, not a measured per-element pennation.
Entities with no OpenSim counterpart get the same harmonic solve with the ends taken from the canonical
principal_axis and are marked inferred. Nothing in this file measures fibre architecture; the derived
field is constrained by a line-actuator path and the inferred field by a bounding axis.
"""
from pathlib import Path
import argparse
import ctypes
import gzip
import hashlib
import importlib.metadata
import json
import os
import shutil
import sys
import tempfile
import time
import xml.etree.ElementTree as ET

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from igl.copyleft import tetgen

LIBC=ctypes.CDLL(None)
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.spatial.vtk import surface as vtk_surface

MILLARD_MODEL='data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/example3DWalking/subject_walk_scaled.osim'
THELEN_MODEL='data/derived/support-physical-root-nayktxy6/assembled_model.osim'
SNAPSHOT='data/derived/native-stream-smoke-n6e0pvqi/supine/smoke.json'
TET_BUILD='data/derived/muscle-tet-ready-v1'
ANATOMY='data/derived/canonical/anatomy.json'
MECHANICS='data/derived/canonical/mechanics.json'
GEOMETRY_DIRS=('data/raw/anatomy/opensim-models/source/Models/Rajagopal/Geometry',
               'data/raw/anatomy/opensim-models/source/Geometry')
FLAGS='pYq1.414'
END_QUANTILE=.10
SCALE_CLAMP=(.75,1.35)

# OpenSim actuator -> canonical entity name template. Every row is a name identity assertion made here,
# not a datum carried by either source. Compartmented actuators map many-to-one onto the single canonical
# belly; that multiplicity is recorded and the paths are averaged.
NAME_MAP={
 'addbrev':'{s} adductor brevis','addlong':'{s} adductor longus',
 'addmagDist':'{s} adductor magnus','addmagIsch':'{s} adductor magnus',
 'addmagMid':'{s} adductor magnus','addmagProx':'{s} adductor magnus',
 'bflh':'long head of {s} biceps femoris','bfsh':'short head of {s} biceps femoris',
 'edl':'{s} extensor digitorum longus','ehl':'{s} extensor hallucis longus',
 'fdl':'{s} flexor digitorum longus','fhl':'{s} flexor hallucis longus',
 'gaslat':'lateral head of {s} gastrocnemius','gasmed':'medial head of {s} gastrocnemius',
 'glmax1':'{s} gluteus maximus','glmax2':'{s} gluteus maximus','glmax3':'{s} gluteus maximus',
 'glmed1':'{s} gluteus medius','glmed2':'{s} gluteus medius','glmed3':'{s} gluteus medius',
 'glmin1':'{s} gluteus minimus','glmin2':'{s} gluteus minimus','glmin3':'{s} gluteus minimus',
 'grac':'{s} gracilis','iliacus':'{s} iliacus',
 'perbrev':'{s} fibularis brevis','perlong':'{s} fibularis longus',
 'piri':'{s} piriformis','psoas':'{s} psoas major',
 'recfem':'{s} rectus femoris','sart':'{s} sartorius',
 'semimem':'{s} semimembranosus','semiten':'{s} semitendinosus','soleus':'{s} soleus',
 'tfl':'{s} tensor fasciae latae','tibant':'{s} tibialis anterior','tibpost':'{s} tibialis posterior',
 'vasint':'{s} vastus intermedius','vaslat':'{s} vastus lateralis','vasmed':'{s} vastus medialis',
 'arm26_TRIlong':'long head of {s} triceps brachii','arm26_TRIlat':'lateral head of {s} triceps brachii',
 'arm26_TRImed':'medial head of {s} triceps brachii','arm26_BIClong':'long head of {s} biceps brachii',
 'arm26_BICshort':'short head of {s} biceps brachii','arm26_BRA':'{s} brachialis',
 'gait2392_extobl':'{s} external oblique'}
UNMAPPED={'gait2392_intobl':'no canonical internal oblique surface exists in the atlas',
          'gait2392_ercspn':'lumped erector spinae actuator; the atlas resolves iliocostalis, longissimus '
                            'and spinalis separately and no unique counterpart exists'}


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write(path,value):Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
def unit(v,axis=-1):
    n=np.linalg.norm(v,axis=axis,keepdims=True)
    return np.divide(v,n,out=np.zeros_like(v),where=n>0)
def deg(cosine):return float(np.degrees(np.arccos(np.clip(cosine,-1.,1.))))


def euler_xyz(angles):
    x,y,z=angles
    cx,sx,cy,sy,cz,sz=np.cos(x),np.sin(x),np.cos(y),np.sin(y),np.cos(z),np.sin(z)
    rx=np.array([[1,0,0],[0,cx,-sx],[0,sx,cx]]);ry=np.array([[cy,0,sy],[0,1,0],[-sy,0,cy]])
    rz=np.array([[cz,-sz,0],[sz,cz,0],[0,0,1]])
    return rx@ry@rz


# ---------------------------------------------------------------- geometry io

def load_surface(path):
    payload=json.loads(gzip.decompress(Path(path).read_bytes()))
    v=np.ascontiguousarray(np.asarray(payload['positions'],float).reshape(-1,3))
    f=np.ascontiguousarray(np.asarray(payload['indices'],np.int64).reshape(-1,3))
    return v,f


def tetrahedralize(v,f,flags=FLAGS):
    """TetGen in a forked child, mesh returned through a temporary npz.

    Same isolation as scripts/verify_muscle_tet_ready_surfaces.py: TetGen aborts the whole process on
    some inputs instead of raising, so a signalled child is a recorded failure rather than a lost run."""
    handle,log=tempfile.mkstemp(suffix='.tetgen');os.close(handle)
    handle,payload=tempfile.mkstemp(suffix='.npz');os.close(handle)
    reader,writer=os.pipe();started=time.monotonic();outcome={'flags':flags}
    pid=os.fork()
    if pid==0:
        try:
            os.close(reader);sink=os.open(log,os.O_WRONLY)
            sys.stdout.flush();sys.stderr.flush();LIBC.fflush(None);os.dup2(sink,1);os.dup2(sink,2)
            child={}
            try:
                result=tetgen.tetrahedralize(v,f,flags=flags)
                child.update(status=int(result[-1]),tet_vertices=int(len(result[0])),tets=int(len(result[1])))
                if child['status']==0 and len(result[1]):
                    np.savez(payload,tv=np.ascontiguousarray(np.asarray(result[0],float)),
                             tt=np.ascontiguousarray(np.asarray(result[1],np.int64)))
            except BaseException as error:
                child.update(status=None,exception=type(error).__name__,message=str(error))
            sys.stdout.flush();sys.stderr.flush();LIBC.fflush(None)
            os.write(writer,json.dumps(child,allow_nan=False).encode());os.close(writer)
        finally:os._exit(0)
    os.close(writer);chunks=[]
    while True:
        chunk=os.read(reader,65536)
        if not chunk:break
        chunks.append(chunk)
    os.close(reader);_,status=os.waitpid(pid,0)
    if chunks:outcome.update(json.loads(b''.join(chunks).decode()))
    if os.WIFSIGNALED(status):
        outcome.update(status=None,exception='ProcessAborted',
                       message='TetGen terminated the process with signal %d'%os.WTERMSIG(status))
    elif not chunks:outcome.update(status=None,exception='NoResult',message='TetGen child exited without a result')
    outcome['tetgen_stdout']=Path(log).read_text(errors='replace').strip().splitlines();Path(log).unlink()
    tv=tt=None
    if outcome.get('status')==0 and outcome.get('tets',0)>0 and Path(payload).stat().st_size:
        with np.load(payload) as data:tv=data['tv'];tt=data['tt']
    Path(payload).unlink()
    outcome['succeeded']=tv is not None
    outcome['seconds']=time.monotonic()-started
    return outcome,tv,tt


def boundary_faces(tt):
    faces=np.concatenate([tt[:,[0,2,1]],tt[:,[0,1,3]],tt[:,[1,2,3]],tt[:,[0,3,2]]])
    key=np.sort(faces,axis=1)
    _,index,count=np.unique(key,axis=0,return_index=True,return_counts=True)
    return faces[index[count==1]]


def tet_gradients(tv,tt):
    """Constant per-tet shape-function gradients (m^-1) and signed volumes (m^3)."""
    d=np.stack([tv[tt[:,1]]-tv[tt[:,0]],tv[tt[:,2]]-tv[tt[:,0]],tv[tt[:,3]]-tv[tt[:,0]]],axis=1)
    volume=np.linalg.det(d)/6.
    inverse=np.linalg.inv(d)
    g=np.empty((len(tt),4,3))
    g[:,1,:]=inverse[:,:,0];g[:,2,:]=inverse[:,:,1];g[:,3,:]=inverse[:,:,2]
    g[:,0,:]=-g[:,1:,:].sum(axis=1)
    return g,volume


def harmonic(tv,tt,gradients,volume,fixed,values):
    n=len(tv);weight=np.abs(volume)
    rows=np.repeat(tt,4,axis=1).ravel();cols=np.tile(tt,(1,4)).ravel()
    data=(weight[:,None,None]*np.einsum('nid,njd->nij',gradients,gradients)).ravel()
    k=sp.coo_matrix((data,(rows,cols)),shape=(n,n)).tocsr()
    free=np.ones(n,bool);free[fixed]=False
    phi=np.zeros(n);phi[fixed]=values
    idx=np.flatnonzero(free)
    if len(idx):
        a=k[idx][:,idx].tocsc();b=-(k[idx][:,fixed]@values)
        phi[idx]=spl.spsolve(a,b)
    if not np.isfinite(phi).all():raise ValueError('nonfinite harmonic potential')
    return phi


def polyline_parameter(points,poly):
    """Arc length along `poly` of the closest point, terminal segments extended to infinity."""
    d=poly[1:]-poly[:-1];length=np.linalg.norm(d,axis=1)
    if not np.all(length>0):raise ValueError('degenerate polyline segment')
    cumulative=np.r_[0.,np.cumsum(length)]
    best_distance=np.full(len(points),np.inf);best_tau=np.zeros(len(points))
    for i in range(len(d)):
        t=((points-poly[i])@d[i])/length[i]**2
        lo=-np.inf if i==0 else 0.;hi=np.inf if i==len(d)-1 else 1.
        t=np.clip(t,lo,hi)
        distance=np.linalg.norm(points-(poly[i]+t[:,None]*d[i]),axis=1)
        take=distance<best_distance
        best_distance[take]=distance[take];best_tau[take]=cumulative[i]+t[take]*length[i]
    return best_tau,best_distance


def end_patches(tau,surface_nodes,quantile=END_QUANTILE):
    value=tau[surface_nodes]
    low=np.quantile(value,quantile);high=np.quantile(value,1.-quantile)
    if not high>low:raise ValueError('degenerate end separation')
    a=surface_nodes[value<=low];b=surface_nodes[value>=high]
    if not len(a) or not len(b) or len(np.intersect1d(a,b)):raise ValueError('degenerate end patches')
    return a,b


def pennate(fibre,centroids,volume,alpha):
    """Tilt every fibre by alpha in the plane spanned by the local fibre and the muscle's thinnest axis.

    One global tilt axis, so the result stays as smooth as the harmonic field. This is a uniform
    unipennate prior standing in for architecture the atlas does not carry: no aponeurosis is located,
    no bipennate reversal is represented and the sense of the tilt is arbitrary."""
    if alpha<=0:return fibre.copy(),0,None
    weighted=centroids-(centroids*volume[:,None]).sum(0)/volume.sum()
    covariance=(weighted*volume[:,None]).T@weighted/volume.sum()
    minor=np.linalg.eigh(covariance)[1][:,0]
    w=minor-(fibre@minor)[:,None]*fibre
    norm=np.linalg.norm(w,axis=1);usable=norm>1e-9
    w=np.divide(w,norm[:,None],out=np.zeros_like(w),where=norm[:,None]>1e-9)
    out=fibre.copy()
    out[usable]=np.cos(alpha)*fibre[usable]+np.sin(alpha)*w[usable]
    return unit(out),int((~usable).sum()),minor


# ---------------------------------------------------------------- registration

def rigid(a,b,weight=None):
    weight=np.ones(len(a)) if weight is None else weight
    weight=weight/weight.sum();ac=weight@a;bc=weight@b
    u,_,vh=np.linalg.svd(((a-ac)*weight[:,None]).T@(b-bc))
    r=vh.T@np.diag([1.,1.,np.sign(np.linalg.det(vh.T@u.T))])@u.T
    return r,bc-r@ac


def icp(src,dst,r,t,iterations=120,trim=.9):
    tree=cKDTree(dst);previous=None
    for _ in range(iterations):
        distance,index=tree.query((r@src.T).T+t)
        keep=distance<=max(np.quantile(distance,trim),1e-12)
        r,t=rigid(src[keep],dst[index][keep])
        value=float(np.sqrt((distance[keep]**2).mean()))
        if previous is not None and abs(previous-value)<1e-12:break
        previous=value
    distance,_=tree.query((r@src.T).T+t)
    return r,t,float(np.sqrt((distance[keep]**2).mean())),float(np.median(distance)),float(distance.max())


def seeded_icp(src,dst,seed,trim=.9):
    """Clamped uniform PCA scale, then trimmed ICP seeded from the whole-body rigid fit.

    Seeding matters: a free principal-axis initialisation lets nearly bilaterally symmetric clouds
    (torso, foot) converge onto their own mirror image at an almost identical residual."""
    sc=src-src.mean(0);dc=dst-dst.mean(0)
    scale=float(np.clip(np.sqrt((dc**2).sum()/len(dc)/((sc**2).sum()/len(sc))),*SCALE_CLAMP))
    origin=src.mean(0);scaled=origin+scale*sc
    return scale,origin,icp(scaled,dst,seed,dst.mean(0)-seed@scaled.mean(0),trim=trim)


def opensim_bodies(model):
    return {b.get('name'):b for b in ET.parse(model).getroot().findall('.//BodySet/objects/Body')}


def body_cloud(body,resolve):
    out=[]
    for mesh in body.findall('./attached_geometry/Mesh'):
        if mesh.findtext('socket_frame').strip()!='..':raise ValueError('offset-framed display mesh')
        scale=np.fromstring(mesh.findtext('scale_factors'),sep=' ')
        vertices,_=vtk_surface(resolve(mesh.findtext('mesh_file')))
        out.append(np.asarray(vertices,float)*scale)
    return np.concatenate(out)


def register(root,spec,trim=.9):
    mechanics=json.loads((root/MECHANICS).read_bytes())
    snapshot=json.loads((root/SNAPSHOT).read_bytes())['initial']
    bodies=opensim_bodies(root/MILLARD_MODEL)
    def resolve(name):
        for directory in GEOMETRY_DIRS:
            path=root/directory/name
            if path.exists():return path
        raise FileNotFoundError(name)
    source_landmarks=[];canonical_landmarks=[]
    for name,entry in mechanics['registration'].items():
        ground=np.asarray(snapshot['bodies'][name]['transform_ground'],float)
        source_landmarks.append(ground[:3,:3]@np.asarray(snapshot['bodies'][name]['mass_center_local_m'],float)+ground[:3,3])
        canonical_landmarks.append((np.min([spec[i]['bounds_m']['min'] for i in entry['canonical_bones']],axis=0)+
                                    np.max([spec[i]['bounds_m']['max'] for i in entry['canonical_bones']],axis=0))/2)
    source_landmarks=np.array(source_landmarks);canonical_landmarks=np.array(canonical_landmarks)
    seed,offset=rigid(source_landmarks,canonical_landmarks)
    global_residual=(source_landmarks@seed.T+offset)-canonical_landmarks
    out={'__global__':{'rotation':seed,'translation':offset,
        'basis':'Proper-rigid fit of the 22 source segment mass centres onto their canonical bone-envelope '
                'centres. Neither side is a measured landmark; this only fixes the gross frame orientation '
                'that seeds every per-body ICP.',
        'rms_landmark_residual_m':float(np.sqrt((global_residual**2).sum(1).mean())),
        'max_landmark_residual_m':float(np.linalg.norm(global_residual,axis=1).max())}}
    for name,entry in mechanics['registration'].items():
        source=body_cloud(bodies[name],resolve)
        target=np.concatenate([load_surface(root/spec[i]['reference_geometry']['path'])[0]
                               for i in entry['canonical_bones']])
        ground=np.asarray(snapshot['bodies'][name]['transform_ground'],float)
        posed=(ground[:3,:3]@source.T).T+ground[:3,3]
        scale,scale_origin,(r,t,rms,median,worst)=seeded_icp(posed,target,seed,trim=trim)
        out[name]={'rotation':r,'translation':t,'scale':scale,'scale_origin':scale_origin,'ground':ground,
                   'source_points':int(len(source)),'canonical_points':int(len(target)),
                   'canonical_bones':entry['canonical_bones'],
                   'rotation_vs_global_seed_deg':deg((np.trace(r@seed.T)-1)/2),
                   'trimmed_rms_m':rms,'median_nearest_m':median,'max_nearest_m':worst}
    return out,mechanics,snapshot


def to_canonical(registration,body,local):
    """Body-local metres -> canonical world metres. ICP was fitted on scale*(posed - mean) + mean."""
    entry=registration[body];ground=entry['ground']
    posed=(ground[:3,:3]@np.atleast_2d(np.asarray(local,float)).T).T+ground[:3,3]
    scaled=entry['scale_origin']+entry['scale']*(posed-entry['scale_origin'])
    return (entry['rotation']@scaled.T).T+entry['translation']


def direction_to_canonical(registration,body,vector):
    return registration[body]['rotation']@np.asarray(vector,float)


# ---------------------------------------------------------------- opensim paths

def joint_table(model):
    joints={}
    for joint in ET.parse(model).getroot().findall('.//JointSet/objects/*'):
        frames={f.get('name'):f for f in joint.findall('frames/PhysicalOffsetFrame')}
        def side(socket):
            frame=frames[joint.findtext(socket).rsplit('/',1)[-1]]
            return (frame.findtext('socket_parent').rsplit('/',1)[-1],
                    np.fromstring(frame.findtext('translation'),sep=' '),
                    euler_xyz(np.fromstring(frame.findtext('orientation'),sep=' ')))
        parent,parent_offset,_=side('socket_parent_frame')
        child,child_offset,child_rotation=side('socket_child_frame')
        axes=[]
        if joint.tag=='PinJoint':
            axes.append({'coordinate':joint.find('coordinates/Coordinate').get('name'),'axis':np.array([0.,0,1])})
        coupled=[]
        for ta in joint.findall('SpatialTransform/TransformAxis'):
            coordinate=(ta.findtext('coordinates') or '').strip()
            if not coordinate:continue
            function=[c.tag for c in ta if c.tag not in ('coordinates','axis')]
            # Only a unit LinearFunction rotation axis carries a moment arm this file can state without
            # differentiating the coupler; polynomial and multiplier axes are recorded and skipped.
            if ta.get('name','').startswith('rotation') and function==['LinearFunction'] and \
               float(ta.findtext('LinearFunction/coefficients').split()[0])==1.:
                axes.append({'coordinate':coordinate,'axis':np.fromstring(ta.findtext('axis'),sep=' ')})
            else:coupled.append({'axis':ta.get('name'),'coordinate':coordinate,'function':function})
        joints[joint.get('name')]={'type':joint.tag,'parent':parent,'child':child,
                                   'child_offset_m':child_offset,'child_rotation':child_rotation,
                                   'parent_offset_m':parent_offset,'axes':axes,'coupled_axes_ignored':coupled}
    return joints


def subtree(joints,body):
    seen={body};frontier=[body]
    while frontier:
        current=frontier.pop()
        for j in joints.values():
            if j['parent']==current and j['child'] not in seen:
                seen.add(j['child']);frontier.append(j['child'])
    return seen


def muscle_paths(root):
    rows={}
    for label,relative,tag in (('millard',MILLARD_MODEL,'Millard2012EquilibriumMuscle'),
                               ('thelen',THELEN_MODEL,'Thelen2003Muscle')):
        tree=ET.parse(root/relative).getroot()
        for muscle in tree.iter(tag):
            path=muscle.find('GeometryPath')
            if path is None:continue
            points=[(p.findtext('socket_parent_frame').rsplit('/',1)[-1],
                     np.fromstring(p.findtext('location'),sep=' ')) for p in path.findall('PathPointSet/objects/*')]
            kinds={p.tag for p in path.findall('PathPointSet/objects/*')}
            if kinds!={'PathPoint'}:raise ValueError('unsupported path point type: '+str(kinds))
            rows[muscle.get('name')]={'model':label,'source_path':relative,'points':points,
                'wrap_objects':[w.findtext('wrap_object') for w in path.findall('PathWrapSet/objects/PathWrap')],
                'pennation_angle_at_optimal_rad':float(muscle.findtext('pennation_angle_at_optimal')),
                'optimal_fiber_length_m':float(muscle.findtext('optimal_fiber_length')),
                'tendon_slack_length_m':float(muscle.findtext('tendon_slack_length')),
                'max_isometric_force_n':float(muscle.findtext('max_isometric_force'))}
    return rows


def moment_arms(joints,registration,path_points):
    """Signed moment arm of unit tension about every rotational axis of every crossed joint, in the
    OpenSim ground frame so no registration enters the model side.

    For a path with no active wrapping this is the model's own moment arm by virtual work: the segment
    straddling the joint carries the whole force onto the distal subtree. Wrapping is ignored, so a
    wrapped path's value is a straight-line proxy."""
    bodies=[b for b,_ in path_points]
    ground=np.array([registration[b]['ground'][:3,:3]@np.asarray(p,float)+registration[b]['ground'][:3,3]
                     for b,p in path_points])
    out=[]
    for name,joint in joints.items():
        if not joint['axes']:continue
        distal=subtree(joints,joint['child'])
        inside=np.array([b in distal for b in bodies])
        if inside.all() or not inside.any():continue
        crossing=[i for i in range(len(bodies)-1) if inside[i]!=inside[i+1]]
        if not crossing:continue
        i=crossing[0]
        far,near=(i,i+1) if inside[i] else (i+1,i)
        tension=unit(ground[near]-ground[far])
        child=registration[joint['child']]['ground']
        centre=child[:3,:3]@joint['child_offset_m']+child[:3,3]
        lever=ground[far]-centre
        row={'joint':name,'type':joint['type'],'child':joint['child'],
             'crossing_segment':[int(far),int(near)],
             'joint_centre_source_m':centre.tolist(),
             'joint_centre_canonical_m':to_canonical(registration,joint['child'],joint['child_offset_m'])[0].tolist(),
             'child_body_scale':registration[joint['child']]['scale'],
             'coupled_axes_ignored':joint['coupled_axes_ignored'],'axes':{}}
        for axis in joint['axes']:
            direction=unit(child[:3,:3]@joint['child_rotation']@axis['axis'])
            row['axes'][axis['coordinate']]={'axis_source':direction.tolist(),
                'axis_canonical':unit(direction_to_canonical(registration,joint['child'],
                                                             joint['child_rotation']@axis['axis'])).tolist(),
                'moment_arm_m':float(direction@np.cross(lever,tension))}
        out.append(row)
    return out


# ---------------------------------------------------------------- field per entity

def solve_entity(tv,tt,tau,mode):
    gradients,volume=tet_gradients(tv,tt)
    if (volume<=0).any():
        flip=volume<0;tt=tt.copy();tt[flip]=tt[flip][:,[0,1,3,2]]
        gradients,volume=tet_gradients(tv,tt)
    if (volume<=0).any():raise ValueError('degenerate tet volume')
    faces=boundary_faces(tt);surface_nodes=np.unique(faces)
    rows=tt[:,[0,0,0,1,1,2]].ravel();cols=tt[:,[1,2,3,2,3,3]].ravel()
    label=connected_components(sp.coo_matrix((np.ones(len(rows)),(rows,cols)),shape=(len(tv),len(tv))),
                               directed=False)[1]
    phi=np.zeros(len(tv));components=[];origin_nodes=[];insertion_nodes=[]
    for c in np.unique(label[tt[:,0]]):
        nodes=np.flatnonzero(label==c)
        local_surface=np.intersect1d(nodes,surface_nodes)
        try:
            a,b=end_patches(tau,local_surface);resolved='harmonic'
        except ValueError:
            a=b=None;resolved='axis_fallback'
        if resolved=='harmonic':
            sub=tt[label[tt[:,0]]==c];remap=-np.ones(len(tv),np.int64);remap[nodes]=np.arange(len(nodes))
            local=remap[sub];g,v=tet_gradients(tv[nodes],local)
            fixed=np.r_[remap[a],remap[b]];values=np.r_[np.zeros(len(a)),np.ones(len(b))]
            phi[nodes]=harmonic(tv[nodes],local,g,v,fixed,values)
            origin_nodes.append(a);insertion_nodes.append(b)
        components.append({'component':int(c),'nodes':int(len(nodes)),'resolution':resolved,
                           'origin_patch_nodes':0 if a is None else int(len(a)),
                           'insertion_patch_nodes':0 if b is None else int(len(b))})
    fibre=np.einsum('nid,ni->nd',gradients,phi[tt])
    magnitude=np.linalg.norm(fibre,axis=1)
    fallback=np.zeros(len(tt),bool)
    for record in components:
        if record['resolution']=='axis_fallback':
            fallback|=label[tt[:,0]]==record['component']
    axis=unit(np.asarray(mode,float))
    if fallback.any():fibre[fallback]=axis;magnitude[fallback]=np.nan
    # A tet whose four nodes all sit in one Dirichlet patch carries no potential difference and so no
    # gradient. That is a boundary-condition artefact at the two end caps, not an interior singularity:
    # such tets inherit the nearest resolved tet's direction and are counted separately.
    solved=~fallback
    resolved=solved&(magnitude>1e-6*(np.nanmedian(magnitude[solved]) if solved.any() else 1.))
    flat=solved&~resolved
    if flat.any():
        centroids=tv[tt].mean(axis=1)
        if resolved.any():
            _,index=cKDTree(centroids[resolved]).query(centroids[flat])
            fibre[flat]=fibre[resolved][index]
        else:fibre[flat]=axis
    ends=(np.concatenate(origin_nodes) if origin_nodes else np.zeros(0,np.int64),
          np.concatenate(insertion_nodes) if insertion_nodes else np.zeros(0,np.int64))
    return (tt,gradients,volume,phi,unit(fibre),magnitude,surface_nodes,faces,components,fallback,flat,ends)


def field_statistics(tt,volume,fibre,magnitude,faces_shared,flat=None,phi=None):
    total=float(volume.sum())
    norms=np.linalg.norm(fibre,axis=1)
    flat=np.zeros(len(tt),bool) if flat is None else flat
    interior=np.isfinite(magnitude)&~flat
    stats={'tets':int(len(tt)),'volume_m3':total,
           'unit_norm_max_deviation':float(np.abs(norms-1).max()),
           'constant_potential_tets':int(flat.sum()),
           'constant_potential_volume_fraction':float(volume[flat].sum()/total),
           'gradient_min_per_m':float(magnitude[interior].min()) if interior.any() else None,
           'gradient_median_per_m':float(np.median(magnitude[interior])) if interior.any() else None,
           'near_singular_tets':None,'volume_fraction_below_decile_gradient':None}
    if phi is not None:
        # A linear-tet stiffness matrix is not an M-matrix on obtuse elements, so the discrete maximum
        # principle can fail. The overshoot is recorded rather than clamped.
        stats['potential_min']=float(phi.min());stats['potential_max']=float(phi.max())
        stats['potential_overshoot']=float(max(-phi.min(),phi.max()-1.,0.))
        outside=np.any((phi[tt]<-1e-9)|(phi[tt]>1+1e-9),axis=1)
        stats['overshoot_volume_fraction']=float(volume[outside].sum()/total)
    if interior.any():
        median=float(np.median(magnitude[interior]))
        stats['near_singular_tets']=int(np.count_nonzero(magnitude[interior]<1e-3*median))
        stats['volume_fraction_below_decile_gradient']=float(volume[interior][magnitude[interior]<.1*median].sum()/total)
    if len(faces_shared):
        a,b=faces_shared[:,0],faces_shared[:,1]
        angle=np.degrees(np.arccos(np.clip(np.einsum('ij,ij->i',fibre[a],fibre[b]),-1,1)))
        w=.5*(volume[a]+volume[b])
        stats.update(interface_count=int(len(a)),
                     interface_angle_median_deg=float(np.median(angle)),
                     interface_angle_p95_deg=float(np.percentile(angle,95)),
                     interface_angle_max_deg=float(angle.max()),
                     interface_angle_volume_weighted_mean_deg=float((angle*w).sum()/w.sum()),
                     interface_volume_fraction_over_30deg=float(w[angle>30].sum()/w.sum()))
        # Where the potential is flat the gradient direction is undefined, so a discontinuity there is a
        # statement about the domain, not about the field. Report the conditioned subset separately.
        if interior.any():
            good=interior[a]&interior[b]&(magnitude[a]>.1*median)&(magnitude[b]>.1*median)
            if good.any():
                stats.update(conditioned_interface_count=int(good.sum()),
                             conditioned_interface_angle_median_deg=float(np.median(angle[good])),
                             conditioned_interface_angle_p95_deg=float(np.percentile(angle[good],95)),
                             conditioned_interface_angle_max_deg=float(angle[good].max()))
    return stats


def local_direction(centroids,volume,fibre,point,fraction=.1):
    """Volume-weighted mean fibre over the tets nearest an attachment patch.

    A muscle whose canonical surface includes its long tendon has a whole-muscle mean direction that is
    the belly's, not the tendon's; the force it hands to the distal bone leaves along the local one."""
    distance=np.linalg.norm(centroids-point,axis=1)
    near=distance<=max(np.quantile(distance,fraction),0.)
    if not near.any():near=distance<=distance.min()
    return unit((fibre[near]*volume[near,None]).sum(0)),int(near.sum())


def shared_faces(tt):
    faces=np.concatenate([tt[:,[0,2,1]],tt[:,[0,1,3]],tt[:,[1,2,3]],tt[:,[0,3,2]]])
    owner=np.tile(np.arange(len(tt)),4)
    key=np.sort(faces,axis=1)
    order=np.lexsort(key.T[::-1]);ordered=key[order];owned=owner[order]
    same=np.all(ordered[1:]==ordered[:-1],axis=1)
    index=np.flatnonzero(same)
    return np.stack([owned[index],owned[index+1]],axis=1)


# ---------------------------------------------------------------- build

def build(output,limit=0,only=None,trim=.9):
    output=Path(output).resolve()
    if output.exists() or not output.is_relative_to(ROOT):raise ValueError('Fresh owned output directory required')
    started=time.monotonic();output.mkdir(parents=True)
    (output/'fibre').mkdir();(output/'inputs').mkdir()
    shutil.copyfile(__file__,output/'inputs/build_muscle_fibre_field.py')
    anatomy=json.loads((ROOT/ANATOMY).read_bytes())
    spec={e['id']:e for e in anatomy['entities']}
    tet_build=ROOT/TET_BUILD
    records=[json.loads(l) for l in (tet_build/'entities.jsonl').read_text().splitlines()]
    by_name={}
    for record in records:by_name.setdefault(spec[record['entity_id']]['name'],[]).append(record['entity_id'])
    registration,mechanics,snapshot=register(ROOT,spec,trim=trim)
    joints=joint_table(ROOT/MILLARD_MODEL)
    paths=muscle_paths(ROOT)
    write(output/'registration.json',{'schema':'ihm.muscle-fibre-registration.v1',
        'method':'per-OpenSim-body trimmed point-to-point ICP of the scaled source bone display meshes onto the '
                 'union of the canonical bone surfaces named by mechanics.json registration, at one clamped '
                 'uniform PCA scale, seeded by the whole-body rigid mass-centre fit recorded in global_seed',
        'trim_fraction':trim,'scale_clamp':list(SCALE_CLAMP),'snapshot':SNAPSHOT,
        'basis':'Bone-surface shape agreement only. No homologous landmark is measured and the residual is not '
                'an anatomical registration accuracy claim.',
        'global_seed':{k:(v.tolist() if isinstance(v,np.ndarray) else v)
                       for k,v in registration['__global__'].items()},
        'bodies':{name:{k:(v.tolist() if isinstance(v,np.ndarray) else v)
                        for k,v in entry.items() if k!='ground'}
                  for name,entry in registration.items() if name!='__global__'}})

    # correspondence
    assignments={};correspondence=[]
    for muscle,row in sorted(paths.items()):
        base=muscle[:-2] if muscle[-2:] in ('_r','_l') else muscle
        side='right' if muscle.endswith('_r') else 'left' if muscle.endswith('_l') else None
        if base in UNMAPPED:
            correspondence.append({'opensim_muscle':muscle,'canonical_entity_id':None,
                                   'reason':UNMAPPED[base]});continue
        if base not in NAME_MAP or side is None:raise ValueError('Unmapped source actuator: '+muscle)
        name=NAME_MAP[base].format(s=side)
        candidates=by_name.get(name,[])
        if len(candidates)!=1:raise ValueError('Canonical target is not unique: '+name)
        ident=candidates[0]
        polyline=np.concatenate([to_canonical(registration,body,local) for body,local in row['points']])
        assignments.setdefault(ident,[]).append({'opensim_muscle':muscle,'polyline_m':polyline,**row})
        correspondence.append({'opensim_muscle':muscle,'canonical_entity_id':ident,'canonical_name':name,
                               'model':row['model'],'source_path':row['source_path'],
                               'path_points':len(row['points']),'wrap_objects':row['wrap_objects'],
                               'pennation_angle_at_optimal_rad':row['pennation_angle_at_optimal_rad']})
    write(output/'correspondence.json',{'schema':'ihm.muscle-fibre-correspondence.v1',
        'basis':'Explicit anatomical name identity asserted in this script between the OpenSim actuator and the '
                'canonical surface. Neither source carries the mapping. Compartmented actuators are many-to-one '
                'onto one canonical belly and their registered paths are averaged.',
        'name_map':NAME_MAP,'declared_unmapped':UNMAPPED,
        'opensim_actuators':len(paths),'mapped_actuators':sum(r['canonical_entity_id'] is not None for r in correspondence),
        'canonical_entities_with_path':len(assignments),'rows':correspondence})

    entities=sorted(records,key=lambda r:spec[r['entity_id']]['name'])
    if only:entities=[r for r in entities if r['entity_id'] in only]
    if limit:entities=entities[:limit]
    emitted=[];failures=[]
    for count,record in enumerate(entities,1):
        ident=record['entity_id'];entity=spec[ident]
        geometry=tet_build/record['output_path']
        v,f=load_surface(geometry)
        outcome,tv,tt=tetrahedralize(v,f)
        if not outcome['succeeded']:
            failures.append({'entity_id':ident,'name':entity['name'],'stage':'tetrahedralize',
                             'tetgen':{k:outcome.get(k) for k in ('status','tets','exception','message')},
                             'tetgen_stdout':outcome['tetgen_stdout']})
            continue
        tt=tt.astype(np.int64)
        assigned=assignments.get(ident,[])
        if assigned:
            taus=[]
            for row in assigned:
                tau,_=polyline_parameter(tv,row['polyline_m'])
                span=tau.max()-tau.min()
                if span<=0:raise ValueError('degenerate path span')
                taus.append((tau-tau.min())/span)
            tau=np.mean(taus,axis=0)
            provenance='derived'
            pennation=float(np.mean([r['pennation_angle_at_optimal_rad'] for r in assigned]))
            axis_prior=unit(assigned[0]['polyline_m'][-1]-assigned[0]['polyline_m'][0])
        else:
            axis=np.asarray(entity['principal_axis'],float);axis=unit(axis)
            if axis[1]<0:axis=-axis
            tau=(tv-tv.mean(0))@axis
            tau=(tau-tau.min())/max(tau.max()-tau.min(),1e-12)
            provenance='inferred';pennation=0.;axis_prior=axis
        try:
            tt,gradients,volume,phi,harmonic_fibre,magnitude,surface_nodes,faces,components,fallback,flat,ends=\
                solve_entity(tv,tt,tau,axis_prior)
        except Exception as error:
            failures.append({'entity_id':ident,'name':entity['name'],'stage':'harmonic',
                             'exception':type(error).__name__,'message':str(error)});continue
        centroids=tv[tt].mean(axis=1)
        boundary_points=tv[surface_nodes]
        fibre,unpennated,tilt_axis=(harmonic_fibre,0,None) if pennation<=0 else \
            pennate(harmonic_fibre,centroids,volume,pennation)
        interfaces=shared_faces(tt)
        stats=field_statistics(tt,volume,fibre,magnitude,interfaces,flat,phi)
        stats_harmonic=field_statistics(tt,volume,harmonic_fibre,magnitude,interfaces,flat)
        weight=volume[:,None]
        integrated=unit((harmonic_fibre*weight).sum(0))
        integrated_pennate=unit((fibre*weight).sum(0))
        coherence=float(np.linalg.norm((harmonic_fibre*weight).sum(0))/volume.sum())
        origin_centre=tv[ends[0]].mean(0) if len(ends[0]) else centroids.mean(0)
        insertion_centre=tv[ends[1]].mean(0) if len(ends[1]) else centroids.mean(0)
        mesh_chord=unit(insertion_centre-origin_centre)
        volume_centre=(centroids*weight).sum(0)/volume.sum()
        local_origin,origin_local_tets=local_direction(centroids,volume,harmonic_fibre,origin_centre)
        local_insertion,insertion_local_tets=local_direction(centroids,volume,harmonic_fibre,insertion_centre)
        row={'entity_id':ident,'name':entity['name'],'provenance':provenance,
             'geometry_path':record['output_path'],'geometry_sha256':record['output_sha256'],
             'tetgen':{'flags':FLAGS,'tets':int(len(tt)),'tet_vertices':int(len(tv)),'seconds':outcome['seconds']},
             'components':components,'axis_fallback_tets':int(fallback.sum()),
             'pennation_angle_applied_rad':pennation,'pennation_degenerate_tets':unpennated,
             'pennation_tilt_axis':None if tilt_axis is None else tilt_axis.tolist(),
             'pennation_basis':None if pennation<=0 else
                'OpenSim pennation_angle_at_optimal applied uniformly, rotating every fibre by that angle in the '
                'plane spanned by the local fibre and the muscle volume\'s thinnest principal axis; a whole-muscle '
                'unipennate prior, not a measured per-element pennation',
             'fibre_orientation_convention':'origin patch (phi=0) toward insertion patch (phi=1)',
             'integrated_fibre_direction':integrated.tolist(),
             'integrated_fibre_direction_pennated':integrated_pennate.tolist(),
             'fibre_coherence':coherence,'field':stats,'field_harmonic':stats_harmonic,
             'surface_nodes':int(len(surface_nodes)),
             'origin_patch_nodes':int(len(ends[0])),'insertion_patch_nodes':int(len(ends[1])),
             'origin_patch_centroid_m':origin_centre.tolist(),
             'insertion_patch_centroid_m':insertion_centre.tolist(),
             'centre_of_volume_m':volume_centre.tolist(),
             'end_patch_chord_direction':mesh_chord.tolist(),
             'local_origin_fibre_direction':local_origin.tolist(),
             'local_insertion_fibre_direction':local_insertion.tolist(),
             'local_patch_tets':[origin_local_tets,insertion_local_tets],
             'integrated_fibre_vs_end_patch_chord_deg':deg(abs(float(integrated@mesh_chord)))}
        if assigned:
            row['opensim']=[]
            surface_tree=cKDTree(boundary_points)
            for entry in assigned:
                poly=entry['polyline_m']
                chord=unit(poly[-1]-poly[0])
                distance,_=surface_tree.query(poly)
                segments=poly[1:]-poly[:-1];length=np.linalg.norm(segments,axis=1)
                tangent=unit((segments*length[:,None]).sum(0))
                bodies=[b for b,_ in entry['points']]
                row['opensim'].append({'opensim_muscle':entry['opensim_muscle'],'model':entry['model'],
                    'source_path':entry['source_path'],'bodies':bodies,'wrap_objects':entry['wrap_objects'],
                    'path_length_m':float(length.sum()),
                    'optimal_fiber_length_m':entry['optimal_fiber_length_m'],
                    'tendon_slack_length_m':entry['tendon_slack_length_m'],
                    'max_isometric_force_n':entry['max_isometric_force_n'],
                    'pennation_angle_at_optimal_rad':entry['pennation_angle_at_optimal_rad'],
                    'chord_direction':chord.tolist(),
                    'integrated_fibre_vs_chord_deg':deg(abs(float(integrated@chord))),
                    'end_patch_chord_vs_path_chord_deg':deg(abs(float(mesh_chord@chord))),
                    'local_insertion_fibre_vs_terminal_segment_deg':deg(abs(float(local_insertion@unit(poly[-1]-poly[-2])))),
                    'local_origin_fibre_vs_initial_segment_deg':deg(abs(float(local_origin@unit(poly[1]-poly[0])))),
                    'integrated_fibre_vs_length_weighted_tangent_deg':deg(abs(float(integrated@tangent))),
                    'pennated_fibre_vs_chord_deg':deg(abs(float(integrated_pennate@chord))),
                    'path_point_to_surface_min_m':float(distance.min()),
                    'path_point_to_surface_median_m':float(np.median(distance)),
                    'path_point_to_surface_max_m':float(distance.max()),
                    'moment_arms':moment_arms(joints,registration,entry['points'])})
            row['mesh_moment_arms']=mesh_moment_arms(joints,registration,assigned,integrated,
                (local_origin,local_insertion),(origin_centre,insertion_centre),volume_centre)
        np.savez_compressed(output/'fibre'/(ident+'.npz'),tet_vertices_m=tv,tets=tt.astype(np.int32),
                            potential=phi,fibre=fibre.astype(np.float32),
                            fibre_harmonic=harmonic_fibre.astype(np.float32))
        row['fibre_path']='fibre/'+ident+'.npz';row['fibre_sha256']=sha(output/'fibre'/(ident+'.npz'))
        emitted.append(row)
        with (output/'entities.jsonl').open('a') as handle:handle.write(json.dumps(row,allow_nan=False)+'\n')
        if count%25==0:print('Field for',count,'of',len(entities),flush=True)
    summarise(output,emitted,failures,entities,assignments,spec,records,paths,registration,started)


def mesh_moment_arms(joints,registration,assigned,integrated,local,ends,volume_centre):
    """Geometric moment arm of the mesh-derived line of action about the same registered joint centres.

    Primary line: the volume-weighted integrated fibre direction carried through the distal attachment
    patch centroid. Three alternatives sit beside it - the fibre direction local to that patch, the
    integrated direction through the belly centre of volume, and the straight end-patch chord - because
    the choice of line is a modelling decision the geometry does not settle, and their spread is reported
    as the sensitivity it is. None of the four is a continuum moment arm: a muscle that presses along a
    bone or spans two joints transmits force through its whole course, not through one attachment, and
    where the canonical surface stops short of the modelled bony attachment none can agree with the path."""
    origin,insertion=ends;local_origin,local_insertion=local
    out=[]
    for entry in assigned:
        bodies=[b for b,_ in entry['points']]
        for name,joint in joints.items():
            if not joint['axes']:continue
            distal=subtree(joints,joint['child'])
            inside=np.array([b in distal for b in bodies])
            if inside.all() or not inside.any():continue
            if not any(inside[i]!=inside[i+1] for i in range(len(bodies)-1)):continue
            centre=to_canonical(registration,joint['child'],joint['child_offset_m'])[0]
            # tau runs path start -> path end, so the mesh insertion patch is the distal end exactly when
            # the last path point sits on the distal side of this joint.
            distal_at_insertion=bool(inside[-1])
            attachment=insertion if distal_at_insertion else origin
            opposite=origin if distal_at_insertion else insertion
            sign=-1. if distal_at_insertion else 1.
            tension=sign*(local_insertion if distal_at_insertion else local_origin)
            chord=unit(opposite-attachment)
            row={'opensim_muscle':entry['opensim_muscle'],'joint':name,
                 'distal_end':'insertion_patch' if distal_at_insertion else 'origin_patch',
                 'attachment_m':attachment.tolist(),'tension_direction':(sign*integrated).tolist(),
                 'axes':{},'axes_local_attachment_direction':{},'axes_at_centre_of_volume':{},
                 'axes_end_patch_chord':{}}
            for axis in joint['axes']:
                direction=unit(direction_to_canonical(registration,joint['child'],
                                                      joint['child_rotation']@axis['axis']))
                row['axes'][axis['coordinate']]=float(direction@np.cross(attachment-centre,sign*integrated))
                row['axes_local_attachment_direction'][axis['coordinate']]=\
                    float(direction@np.cross(attachment-centre,tension))
                row['axes_at_centre_of_volume'][axis['coordinate']]=\
                    float(direction@np.cross(volume_centre-centre,sign*integrated))
                row['axes_end_patch_chord'][axis['coordinate']]=\
                    float(direction@np.cross(attachment-centre,chord))
            out.append(row)
    return out


def subset_summary(rows,basis):
    if not rows:return {'count':0,'basis':basis}
    difference=np.abs([c['difference_m'] for c in rows])
    ratio=[c['mesh_moment_arm_m']/c['opensim_moment_arm_canonical_m'] for c in rows
           if c['opensim_moment_arm_canonical_m']]
    return {'count':len(rows),'basis':basis,
            'sign_agreement':float(np.mean([c['sign_agrees'] for c in rows])),
            'absolute_difference_median_m':float(np.median(difference)),
            'absolute_difference_p90_m':float(np.percentile(difference,90)),
            'absolute_difference_max_m':float(difference.max()),
            'mesh_over_opensim_ratio_median':float(np.median(ratio)) if ratio else None}


def summarise(output,emitted,failures,entities,assignments,spec,records,paths,registration,started):
    derived=[r for r in emitted if r['provenance']=='derived']
    inferred=[r for r in emitted if r['provenance']=='inferred']
    def collect(rows,key,sub='field'):
        return [r[sub][key] for r in rows if r[sub].get(key) is not None]
    comparisons=[]
    for row in derived:
        mesh={(m['opensim_muscle'],m['joint']):m for m in row.get('mesh_moment_arms',[])}
        for entry in row['opensim']:
            for joint in entry['moment_arms']:
                block=mesh.get((entry['opensim_muscle'],joint['joint']))
                if block is None:continue
                for coordinate,axis in joint['axes'].items():
                    if coordinate not in block['axes']:continue
                    scaled=axis['moment_arm_m']*joint['child_body_scale']
                    value=block['axes'][coordinate]
                    comparisons.append({'entity_id':row['entity_id'],'name':row['name'],
                        'opensim_muscle':entry['opensim_muscle'],'joint':joint['joint'],
                        'coordinate':coordinate,'wrapped':bool(entry['wrap_objects']),
                        'opensim_moment_arm_m':axis['moment_arm_m'],
                        'opensim_moment_arm_canonical_m':scaled,'mesh_moment_arm_m':value,
                        'mesh_moment_arm_local_attachment_direction_m':block['axes_local_attachment_direction'][coordinate],
                        'mesh_moment_arm_at_centre_of_volume_m':block['axes_at_centre_of_volume'][coordinate],
                        'mesh_moment_arm_end_patch_chord_m':block['axes_end_patch_chord'][coordinate],
                        'difference_m':value-scaled,
                        'sign_agrees':bool(np.sign(value)==np.sign(scaled)),
                        'free_tendon_fraction':entry['tendon_slack_length_m']/entry['path_length_m'],
                        'tendon_slack_length_m':entry['tendon_slack_length_m'],
                        'path_length_m':entry['path_length_m']})
    write(output/'moment_arms.json',{'schema':'ihm.muscle-fibre-moment-arm.v1',
        'basis':'Both sides are geometric moment arms of unit tension about the same registered joint centre and '
                'axis, in metres. The OpenSim value equals the model -dL/dq only where the path carries no active '
                'wrap object; its canonical counterpart is the same length multiplied by the distal body ICP '
                'scale. The mesh value carries the fibre direction local to the distal attachment patch through '
                'that patch centroid, with whole-muscle-direction, centre-of-volume and end-patch-chord '
                'alternatives beside it. Where the canonical surface stops short of the modelled bony attachment '
                'the two cannot agree; free_tendon_fraction and the path-point-to-surface distances say where.',
        'comparisons':comparisons})
    angles=[e['integrated_fibre_vs_chord_deg'] for r in derived for e in r['opensim']]
    chords=[e['end_patch_chord_vs_path_chord_deg'] for r in derived for e in r['opensim']]
    terminal=[a for r in derived for e in r['opensim']
              for a in (e['local_insertion_fibre_vs_terminal_segment_deg'],
                        e['local_origin_fibre_vs_initial_segment_deg'])]
    unwrapped=[c for c in comparisons if not c['wrapped']]
    summary={'schema':'ihm.muscle-fibre-field.v1',
        'canonical_muscle_entities':len(records),'entities_attempted':len(entities),
        'entities_emitted':len(emitted),'entities_failed':len(failures),
        'derived_entities':len(derived),'inferred_entities':len(inferred),
        'canonical_entities_with_opensim_path':len(assignments),
        'opensim_actuators':len(paths),
        'opensim_actuators_mapped':sum(len(v) for v in assignments.values()),
        'total_tets':int(sum(r['tetgen']['tets'] for r in emitted)),
        'derived_tets':int(sum(r['tetgen']['tets'] for r in derived)),
        'inferred_tets':int(sum(r['tetgen']['tets'] for r in inferred)),
        'total_volume_m3':float(sum(r['field']['volume_m3'] for r in emitted)),
        'derived_volume_m3':float(sum(r['field']['volume_m3'] for r in derived)),
        'unit_norm_max_deviation':max(r['field']['unit_norm_max_deviation'] for r in emitted) if emitted else None,
        'entities_with_axis_fallback_tets':sum(r['axis_fallback_tets']>0 for r in emitted),
        'axis_fallback_tets':int(sum(r['axis_fallback_tets'] for r in emitted)),
        'near_singular_tets':int(sum(r['field']['near_singular_tets'] or 0 for r in emitted)),
        'constant_potential_volume_fraction':float(sum(r['field']['constant_potential_volume_fraction']*
            r['field']['volume_m3'] for r in emitted)/sum(r['field']['volume_m3'] for r in emitted)) if emitted else None,
        'max_potential_overshoot':max(r['field']['potential_overshoot'] for r in emitted) if emitted else None,
        'entities_with_potential_overshoot_over_1pct':[r['entity_id'] for r in emitted
                                                       if r['field']['potential_overshoot']>.01],
        'overshoot_volume_fraction':float(sum(r['field']['overshoot_volume_fraction']*r['field']['volume_m3']
            for r in emitted)/sum(r['field']['volume_m3'] for r in emitted)) if emitted else None,
        'entities_with_near_singular_tets':sum(bool(r['field']['near_singular_tets']) for r in emitted),
        'interface_angle_median_deg_median':float(np.median(collect(emitted,'interface_angle_median_deg'))) if emitted else None,
        'interface_angle_p95_deg_median':float(np.median(collect(emitted,'interface_angle_p95_deg'))) if emitted else None,
        'interface_angle_max_deg_max':float(max(collect(emitted,'interface_angle_max_deg'))) if emitted else None,
        'fibre_coherence_median':float(np.median([r['fibre_coherence'] for r in emitted])) if emitted else None,
        'line_of_action_vs_path_chord_deg':{'count':len(angles),
            'median':float(np.median(angles)) if angles else None,
            'p90':float(np.percentile(angles,90)) if angles else None,
            'max':float(max(angles)) if angles else None} ,
        'moment_arm_comparisons':len(comparisons),
        'moment_arm_unwrapped_comparisons':len(unwrapped),
        'moment_arm_sign_agreement':float(np.mean([c['sign_agrees'] for c in comparisons])) if comparisons else None,
        'moment_arm_unwrapped_sign_agreement':float(np.mean([c['sign_agrees'] for c in unwrapped])) if unwrapped else None,
        'moment_arm_absolute_difference_m':{'median':float(np.median([abs(c['difference_m']) for c in comparisons])) if comparisons else None,
            'p90':float(np.percentile([abs(c['difference_m']) for c in comparisons],90)) if comparisons else None,
            'max':float(max(abs(c['difference_m']) for c in comparisons)) if comparisons else None},
        'moment_arm_line_variants':{key:subset_summary(
                [dict(c,mesh_moment_arm_m=c[key],difference_m=c[key]-c['opensim_moment_arm_canonical_m'],
                      sign_agrees=bool(np.sign(c[key])==np.sign(c['opensim_moment_arm_canonical_m'])))
                 for c in comparisons if abs(c['opensim_moment_arm_canonical_m'])>=.01],
                'axes where the model moment arm is at least 10 mm, using this line')
            for key in ('mesh_moment_arm_m','mesh_moment_arm_local_attachment_direction_m',
                        'mesh_moment_arm_at_centre_of_volume_m','mesh_moment_arm_end_patch_chord_m')},
        'moment_arm_subsets':[
            subset_summary(comparisons,'every crossed rotational axis of every mapped actuator'),
            subset_summary([c for c in comparisons if abs(c['opensim_moment_arm_canonical_m'])>=.01],
                           'axes where the model moment arm is at least 10 mm'),
            subset_summary([c for c in comparisons if not c['wrapped']],
                           'actuators whose path carries no wrap object'),
            subset_summary([c for c in comparisons if c['free_tendon_fraction']<.5],
                           'actuators whose tendon slack length is under half the path length'),
            subset_summary([c for c in comparisons if c['free_tendon_fraction']<.5
                            and abs(c['opensim_moment_arm_canonical_m'])>=.01],
                           'short-tendon actuators on axes where the model moment arm is at least 10 mm')],
        'local_attachment_fibre_vs_terminal_path_segment_deg':{'count':len(terminal),
            'median':float(np.median(terminal)) if terminal else None,
            'p90':float(np.percentile(terminal,90)) if terminal else None,
            'max':float(max(terminal)) if terminal else None},
        'end_patch_chord_vs_path_chord_deg':{'count':len(chords),
            'median':float(np.median(chords)) if chords else None,
            'p90':float(np.percentile(chords,90)) if chords else None,
            'max':float(max(chords)) if chords else None},
        'registration_trimmed_rms_m':{name:entry['trimmed_rms_m'] for name,entry in registration.items()
                                      if name!='__global__'},
        'registration_worst_trimmed_rms_m':max(e['trimmed_rms_m'] for n,e in registration.items() if n!='__global__'),
        'registration_worst_rotation_vs_global_seed_deg':max(e['rotation_vs_global_seed_deg'] for n,e in
                                                             registration.items() if n!='__global__'),
        'registration_global_seed_rms_m':registration['__global__']['rms_landmark_residual_m'],
        'failures':failures,'elapsed_s':time.monotonic()-started}
    write(output/'summary.json',summary)
    inputs={p:sha(ROOT/p) for p in (MILLARD_MODEL,THELEN_MODEL,SNAPSHOT,ANATOMY,MECHANICS,
                                    TET_BUILD+'/manifest.json',TET_BUILD+'/entities.jsonl')}
    artifacts={str(p.relative_to(output)):sha(p) for p in sorted(output.rglob('*')) if p.is_file()}
    write(output/'manifest.json',{'schema':'ihm.muscle-fibre-field.v1',
        'packages':{n:importlib.metadata.version(n) for n in ('libigl','numpy','scipy')},'python':sys.version,
        'inputs_sha256':inputs,'artifacts_sha256':artifacts,'tetgen_flags':FLAGS,'end_quantile':END_QUANTILE,
        'method':'whole-body-seeded per-body trimmed ICP registration of OpenSim bone display meshes onto '
                 'canonical bone surfaces; TetGen pYq1.414; Dirichlet end patches from the extreme deciles of '
                 'arc-length position along the registered path (or the canonical principal axis when no path '
                 'exists); linear-tet FEM Laplace solve; per-tet fibre = normalised grad phi; declared pennation '
                 'applied as a uniform rotation in the plane of the fibre and the volume thinnest axis',
        'canonical_geometry_modified':False,
        'limitations':[
          'No fibre architecture is measured anywhere in this build. A derived field is a harmonic interpolation '
          'between two end patches chosen by a line-actuator path; an inferred field is the same interpolation '
          'between the extremes of a bounding principal axis.',
          'The OpenSim gait model covers a minority of the canonical muscular atlas. Every entity outside that '
          'cover is an engineering prior and is marked inferred.',
          'Registration is rigid plus one uniform scale per body between two different specimens. Its residual '
          'is bone-surface shape agreement, not homologous landmark accuracy, and it bounds every reported angle.',
          'Path wrapping is ignored. A wrapped path is a straight polyline here, so its moment arm and its '
          'direction between via points are proxies.',
          'The canonical surface and the modelled bony attachments do not coincide and do so inconsistently: most '
          'bellies stop short of the OpenSim attachment while some, the long toe flexors among them, carry their '
          'tendon the whole way. Every end patch therefore sits somewhere the path model does not, and the four '
          'reported moment-arm lines disagree by that much. Distances from each path point to the surface are '
          'recorded per actuator so the reader can see which case a muscle is.',
          'Pennation is one declared whole-muscle angle rotated in the plane of the local fibre and the volume '
          'thinnest axis, with an arbitrary sense. Bipennate and multipennate architecture, aponeurosis geometry '
          'and intramuscular pennation gradients are absent.',
          'Compartmented actuators map many-to-one onto one canonical belly; their paths are averaged, which '
          'discards the compartment structure the source model resolves.',
          'The linear-tet Laplace operator is not an M-matrix on obtuse elements, so the discrete maximum '
          'principle fails locally on sliver tets. The overshoot is recorded per entity and not clamped; where '
          'it is large the local gradient direction there is not trustworthy.']})
    print(json.dumps({k:v for k,v in summary.items() if k not in ('failures','registration_trimmed_rms_m')},indent=2))


# ---------------------------------------------------------------- self test

def self_test():
    a=np.random.default_rng(0).normal(size=(50,3))
    angle=.7;r=euler_xyz([.3,-.2,angle]);t=np.array([1.,-2,3])
    fitted,shift=rigid(a,(r@a.T).T+t)
    assert np.allclose(fitted,r,atol=1e-12) and np.allclose(shift,t,atol=1e-12)
    scale,_,(ricp,ticp,rms,median,worst)=seeded_icp(a,(r@a.T).T+t,np.eye(3))
    assert rms<1e-6 and abs(scale-1)<1e-6 and np.allclose(ricp,r,atol=1e-6),(rms,scale)
    mirror=a@np.diag([1.,1,-1])
    _,_,(bad,_,_,_,_)=seeded_icp(a,(r@mirror.T).T+t,np.eye(3))
    assert deg((np.trace(bad@r.T)-1)/2)>1.,'a mirrored target must not be reported as the same rigid fit'

    box=np.array([[0.,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]])
    quads=[(0,3,2,1),(4,5,6,7),(0,1,5,4),(1,2,6,5),(2,3,7,6),(3,0,4,7)]
    faces=np.array([[q[0],q[i],q[i+1]] for q in quads for i in (1,2)],np.int64)
    outcome,tv,tt=tetrahedralize(box,faces,'pq1.414a0.01')
    assert outcome['succeeded'] and len(tt)>50,outcome
    tau=tv[:,0]
    tt2,gradients,volume,phi,fibre,magnitude,surface_nodes,_,components,fallback,flat,ends=solve_entity(tv,tt.astype(np.int64),tau,[1.,0,0])
    assert not fallback.any() and not flat.any() and len(components)==1
    assert tv[ends[0]][:,0].max()<tv[ends[1]][:,0].min(),'end patches must sit at opposite ends of the axis'
    assert abs(volume.sum()-1)<1e-9,volume.sum()
    assert np.abs(phi-tv[:,0]).max()<1e-9,np.abs(phi-tv[:,0]).max()
    assert np.abs(fibre-np.array([1.,0,0])).max()<1e-8
    assert np.abs(np.linalg.norm(fibre,axis=1)-1).max()<1e-12
    assert abs(np.nanmedian(magnitude)-1)<1e-8

    centroids=tv[tt2].mean(axis=1)
    tilted,degenerate,tilt_axis=pennate(fibre,centroids,volume,.3)
    assert degenerate==0 and abs(abs(tilt_axis@np.array([1.,0,0]))-0)<.9
    dot=np.clip(np.einsum('ij,ij->i',tilted,fibre),-1,1)
    assert np.abs(np.arccos(dot)-.3).max()<1e-9
    assert np.abs(np.linalg.norm(tilted,axis=1)-1).max()<1e-12

    poly=np.array([[0.,0,0],[1,0,0],[1,1,0]])
    tau,distance=polyline_parameter(np.array([[.5,0,0],[1,.5,0],[-1,0,0],[1,2,0]]),poly)
    assert np.allclose(tau,[.5,1.5,-1.,3.]) and np.allclose(distance,0)
    low,high=end_patches(np.array([0.,.1,.5,.9,1.]),np.arange(5),quantile=.25)
    assert set(low.tolist())=={0,1} and set(high.tolist())=={3,4}

    joints={'pin':{'type':'PinJoint','parent':'a','child':'b','child_offset_m':np.zeros(3),
                   'child_rotation':np.eye(3),'parent_offset_m':np.zeros(3),'coupled_axes_ignored':[],
                   'axes':[{'coordinate':'q','axis':np.array([0.,0,1])}]}}
    registration={k:{'rotation':np.eye(3),'translation':np.zeros(3),'scale':1.,
                     'scale_origin':np.zeros(3),'ground':np.eye(4)} for k in ('a','b')}
    arms=moment_arms(joints,registration,[('a',np.array([.05,1,0])),('b',np.array([.05,0,0]))])
    assert len(arms)==1 and abs(arms[0]['axes']['q']['moment_arm_m']-.05)<1e-12,arms
    assert not set(NAME_MAP)&set(UNMAPPED)
    assert to_canonical(registration,'b',[1.,2,3]).shape==(1,3)
    print('PASS rigid/ICP recovery, exact linear harmonic field, pennation angle, polyline parameter, '
          'end patches and unit moment arm')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--self-test',action='store_true');p.add_argument('--output',type=Path)
    p.add_argument('--limit',type=int,default=0);p.add_argument('--only',nargs='*')
    p.add_argument('--trim',type=float,default=.9)
    a=p.parse_args()
    if a.self_test:self_test()
    if a.output:build(a.output,a.limit,a.only,a.trim)
    if not a.self_test and not a.output:p.error('Select --self-test or --output')
