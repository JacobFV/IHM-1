#!/usr/bin/env python3
"""Cut the canonical skin into one contact surface per segment.

Contact with the world is never bone against world.  It is skin, over fat and
muscle, over bone, and it is the skin that meets the floor.  This produces the
outer end of that chain in the form the engine can already carry: one closed
triangle mesh per driven segment, in that segment's own frame, admissible as an
OpenSim `ContactMesh`.  The soft tissue between skin and bone is then the elastic
foundation's own layer -- its documented stiffness law is
`k = (1-p)E/((1+p)(1-2p)h)` for a uniform elastic layer of thickness h over a
rigid substrate, and E, p and h are taken from the canonical skin-layer entities
rather than typed in.

What is real here and what is not:

* The SURFACE is the measured canonical exterior skin -- 109,183 triangles,
  1.78 m² -- not a projection, an envelope or a sphere.  Every segment gets the
  skin that rides it, so contact works in any pose, not only the one the
  quadrature was rasterized for.
* The partition is the repo's own `continuous_surface_binding`, a graph-diffused
  skinning weight per skin vertex per segment.  A triangle goes to the argmax of
  its three vertices' mean weight.  It is a HARD partition of a surface that is
  really continuous, so every segment boundary is a seam that in the real body
  does not exist.
* The skin is carried RIGIDLY by its segment.  It does not stretch, slide or
  deform in-plane; the only compliance is the foundation's normal layer.  That is
  the honest limit of this stack: Simbody is a rigid multibody engine and there
  is no deformable continuum anywhere in it.
* Cutting an open surface into pieces leaves each piece open, and SimTK will not
  accept an open mesh.  Each piece is therefore CAPPED, and the area and volume
  the capping invents are reported per segment rather than absorbed.
"""
from pathlib import Path
import argparse,gzip,hashlib,importlib.util,json,sys
import numpy as np
import xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(ROOT/'scripts'))
from ihm.assembly.articulated import CanonicalRegistration
from build_supine_surface_contact import layer_thickness
from build_segment_contact_meshes import measure,controls,simtk_precondition,repair,cap_boundaries

DEFAULT_REFERENCE='data/derived/supine-support-5ma720yd/initial_native.json'
BINDING='data/derived/canonical/continuous_surface_binding.json.gz'
EVIDENCE='data/research/engineered_skin_territories/materialization.json'

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def binding_registration():
    """The OTHER canonical->source map this repo carries, and the pose it was fitted at.

    `CanonicalRegistration.global_map` is an unweighted proper-rigid fit of 22
    approximate COM / bone-envelope-centre correspondences with NO SCALE, and it
    reports its own residual as 123.4 mm RMS.  `bind_anatomy_to_segments.py`
    solved a different problem -- the 33 model coordinates AND one similarity,
    jointly, against 22 bone-group centroids and their principal axes -- and got
    scale 0.96303 with a much smaller residual.  Because the POSE was free in
    that fit, the atlas matches the model only at `reference_pose_rad`, so the
    segment frames have to be taken there and not at the zero pose.

    These two maps are not the same map.  Choosing between them is the whole
    difference between skin that touches the floor and skin that hovers 96 mm
    above it, so the choice is explicit and the gate below decides it.
    """
    binding=json.loads((ROOT/'data/derived/anatomy-segment-binding/binding.json').read_text())
    similarity=np.asarray(binding['similarity_atlas_from_opensim_ground'],dtype=float)
    spec=importlib.util.spec_from_file_location('_render_body_3d',ROOT/'scripts/render_body_3d.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    model=module.OsimModel(ROOT/'data/models/engineering_stance_v1/model.osim')
    rest=model.forward(binding['reference_pose_rad'])
    return np.linalg.inv(similarity),{name:np.asarray(value,dtype=float) for name,value in rest.items()},binding

def bone_clouds():
    """Every body's own bone-mesh vertices, in the body frame, scale baked in.

    The gate a registration has to pass is the one whose answer everybody knows:
    **a body's bones are inside its skin.**  These are the points to test.
    """
    from ihm.spatial.vtk import surface as read_surface
    geometry=ROOT/'data/raw/anatomy/opensim-models/source/Geometry'
    root=ET.parse(ROOT/'data/models/engineering_stance_v1/model.osim').getroot().find('Model')
    out={}
    for body in root.iter('Body'):
        parts=[]
        for mesh in body.iter('Mesh'):
            factors=np.fromstring(mesh.findtext('scale_factors'),sep=' ')
            points,_=read_surface(geometry/mesh.findtext('mesh_file'))
            parts.append(points*factors)
        if parts:out[body.get('name')]=np.concatenate(parts)
    return out

def enclosure(skin_vertices,skin_faces,bone_points,samples=1500,seed=0):
    """Share of a segment's bone vertices that lie inside its skin surface.

    Moller-Trumbore ray parity against the skin's own triangles, written out
    rather than called, because trimesh's own `contains` needs an rtree that is
    not installed here and silently raises.  Two rays per point in opposite
    directions; a point counts as inside only if BOTH give an odd crossing
    count, which throws away the grazing hits that a single ray gets wrong.

    A correct registration prints ~1.0 and a wrong one prints ~0.  That is what
    makes this a gate rather than a diagnostic: a body's bones are inside its
    skin, and no argument about registration quality survives a 0.
    """
    generator=np.random.default_rng(seed)
    points=np.asarray(bone_points,dtype=float)
    if len(points)>samples:points=points[generator.choice(len(points),samples,replace=False)]
    triangles=np.asarray(skin_vertices,dtype=float)[np.asarray(skin_faces)]
    a=triangles[:,0];edge1=triangles[:,1]-a;edge2=triangles[:,2]-a
    direction=generator.normal(size=3);direction/=np.linalg.norm(direction)
    counts=[]
    for sign in (1.,-1.):
        d=sign*direction
        pvec=np.cross(d,edge2);det=(edge1*pvec).sum(axis=1)
        parallel=np.abs(det)<1e-14;inverse=np.where(parallel,0.,1./np.where(parallel,1.,det))
        hits=np.zeros(len(points),dtype=np.int64)
        for start in range(0,len(points),256):
            block=points[start:start+256]
            tvec=block[:,None,:]-a[None,:,:]
            u=(tvec*pvec[None,:,:]).sum(axis=2)*inverse[None,:]
            qvec=np.cross(tvec,edge1[None,:,:])
            v=(qvec*d).sum(axis=2)*inverse[None,:]
            t=(qvec*edge2[None,:,:]).sum(axis=2)*inverse[None,:]
            ok=(~parallel[None,:])&(u>=0)&(v>=0)&(u+v<=1)&(t>1e-9)
            hits[start:start+256]=ok.sum(axis=1)
        counts.append(hits%2==1)
    return float(np.mean(counts[0]&counts[1]))

def skin_layers(mechanics,surface_area):
    """The declared skin layers, so the foundation's E, p and h are measured.

    Thickness comes from the same `layer_thickness` the supine foundation uses,
    which is the explicit shell thickness the canonical entity declares.  The
    elastic foundation's own documented reading of its stiffness is a uniform
    elastic layer of thickness h over a rigid substrate, so the layer this body
    declares IS the parameter, and the stiffness follows from it rather than
    being chosen to make a number come out.
    """
    layers=[e for e in mechanics['entities'] if e['role']=='skin_layer']
    if not layers:raise ValueError('no declared skin layers')
    parameters={(e['material']['young_modulus']['value'],e['material']['poisson_ratio']['value']) for e in layers}
    if len(parameters)!=1:raise ValueError('skin layers disagree on material; no silent averaging')
    young,poisson=parameters.pop()
    thicknesses=[]
    for entity in layers:
        thickness,basis=layer_thickness(entity,surface_area)
        thicknesses.append(dict(id=entity['id'],thickness_m=thickness,thickness_basis=basis))
    total=sum(t['thickness_m'] for t in thicknesses)
    stiffness=(1-poisson)*young/((1+poisson)*(1-2*poisson)*total)
    return dict(layers=thicknesses,youngs_modulus_pa=float(young),poissons_ratio=float(poisson),
                layer_thickness_m=total,stiffness_pa_per_m=stiffness,
                stiffness_basis='k=(1-p)E/((1+p)(1-2p)h) -- the elastic foundation\'s own law for a uniform elastic layer of thickness h over a rigid substrate, with E, p and h taken from the declared skin layers.')

def build(out_dir,reference_path,minimum_faces,registration_choice,warp=None):
    out_dir=Path(out_dir).resolve();meshes=out_dir/'meshes';meshes.mkdir(parents=True,exist_ok=True)
    mechanics=json.loads((ROOT/'data/derived/canonical/mechanics.json').read_text())
    skin=next(e for e in mechanics['entities'] if e['role']=='skin')
    mesh_path=ROOT/skin['reference_geometry']['path']
    if sha(mesh_path)!=skin['reference_geometry']['sha256']:raise ValueError('Skin identity mismatch')
    geometry=json.loads(gzip.decompress(mesh_path.read_bytes()))
    canonical=np.asarray(geometry['positions'],dtype=float).reshape(-1,3)
    faces=np.asarray(geometry['indices'],dtype=np.int64).reshape(-1,3)
    # Only the exterior connected component is skin that can touch anything; the
    # other 99 components are interior surfaces of the same acquired body.
    evidence=json.loads((ROOT/EVIDENCE).read_text())
    exterior=np.asarray(evidence['contact_eligible_triangle_ids'],dtype=np.int64)
    if registration_choice=='binding':
        transform,frames,binding=binding_registration()
        registration_report=dict(choice='binding',
            map='inverse of binding.json similarity_atlas_from_opensim_ground',
            scale=binding['registration']['scale'],method=binding['registration']['method'],
            pose='binding.json reference_pose_rad (the pose the similarity was jointly fitted at)')
    else:
        reference=json.loads((ROOT/reference_path).read_text())
        registration=CanonicalRegistration(mechanics,reference)
        transform=np.linalg.inv(registration.global_map)
        frames={name:np.asarray(value['transform_ground'],dtype=float) for name,value in reference['bodies'].items()}
        registration_report=dict(choice='canonical',map='inverse of CanonicalRegistration.global_map',
            scale=1.0,method=registration.global_fit['basis'],
            rms_landmark_residual_m=registration.global_fit['rms_landmark_residual_m'],
            maximum_landmark_residual_m=registration.global_fit['maximum_landmark_residual_m'],
            pose=str(reference_path)+' t=0 body transforms')
    if warp is None:
        source=canonical@transform[:3,:3].T+transform[:3,3]
    else:
        # One smooth space warp (scripts/skin_warp.py) on top of the binding map, applied to the
        # WHOLE skin before it is cut.  Its base must be exactly the map chosen above, so a zero
        # displacement reproduces the unwarped bundle bit for bit.
        from skin_warp import load_warp
        if registration_choice!='binding':raise ValueError('a skin warp is defined on the binding map only')
        field=load_warp(ROOT/warp)
        if not np.array_equal(field.base,transform):raise ValueError('the warp\'s base is not the binding map')
        source=field.apply(canonical)
        registration_report['warp']=dict(path=str(warp),sha256=sha(ROOT/warp),centres=int(len(field.centres)),
            steps=getattr(field,'steps',None),
            form=('W(x) = flow_1(G x), the time-1 flow of a stationary velocity field integrated by scaling and squaring'
                  if hasattr(field,'steps') else
                  'W(x) = G x + d(G x), d a regularised 3D thin-plate spline')+' (scripts/skin_warp.py); G the binding map above',
            meta=field.meta)
    bones=bone_clouds()
    binding=json.loads(gzip.decompress((ROOT/BINDING).read_bytes()))
    segments=[s['id'] for s in binding['segments']]
    weights=np.asarray(binding['weights'],dtype=np.float32)
    if weights.shape!=(len(canonical),len(segments)):raise ValueError('binding width does not match the skin mesh')
    owner=((weights[faces[:,0]]+weights[faces[:,1]]+weights[faces[:,2]])/3).argmax(axis=1)
    triangles=source[faces[exterior]]
    exterior_area=float(np.linalg.norm(np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0]),axis=1).sum()/2)
    records=[];rejected=[]
    for index,segment in enumerate(segments):
        selected=exterior[owner[exterior]==index]
        if len(selected)<minimum_faces:
            rejected.append(dict(body=segment,reason='fewer than %d exterior triangles'%minimum_faces,faces=int(len(selected)),admitted=False))
            continue
        used,inverse=np.unique(faces[selected],return_inverse=True)
        piece_faces=inverse.reshape(-1,3)
        world=np.linalg.inv(frames[segment])
        local=source[used]@world[:3,:3].T+world[:3,3]
        cut=measure(local,piece_faces)
        reason=simtk_precondition(local,piece_faces);capped=False;cap_area=0.
        if reason is not None:
            clean_vertices,clean_faces=repair(local,piece_faces)
            try:fixed_vertices,fixed_faces,cap_area=cap_boundaries(clean_vertices,clean_faces)
            except Exception as error:
                rejected.append(dict(body=segment,reason='capping failed: '+str(error),raw_reason=reason,
                                     faces=int(len(piece_faces)),cut=cut,admitted=False));continue
            fixed_reason=simtk_precondition(fixed_vertices,fixed_faces)
            if fixed_reason is not None:
                rejected.append(dict(body=segment,reason=fixed_reason,raw_reason=reason,
                                     faces=int(len(piece_faces)),cut=cut,admitted=False));continue
            local,piece_faces,capped=fixed_vertices,fixed_faces,True
        closed=measure(local,piece_faces)
        bone=bones.get(segment)
        contained=None if bone is None else enclosure(local,piece_faces,bone)
        bone_offset=None if bone is None else float(local[:,1].min()-bone[:,1].min())
        target=meshes/('skin_'+segment+'.obj')
        with target.open('w') as handle:
            handle.write('# skin_'+segment+' from '+skin['id']+' exterior component, segment-local\n')
            for v in local:handle.write('v %.9g %.9g %.9g\n'%tuple(v))
            for f in piece_faces:handle.write('f %d %d %d\n'%(f[0]+1,f[1]+1,f[2]+1))
        records.append(dict(body=segment,element='skin_'+segment,mesh_file=target.name,
                            written_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
                            source_geometry=str(mesh_path.relative_to(ROOT)),source_sha256=sha(mesh_path),
                            scale_factors=[1.,1.,1.],simtk_precondition_repaired=capped,
                            exterior_triangles=int(len(selected)),
                            bone_vertices_inside_skin=contained,
                            skin_minus_bone_minimum_y_m=bone_offset,
                            cut_surface_area_m2=cut['surface_area_m2'],
                            capped_surface_area_m2=closed['surface_area_m2'],
                            capped_area_added_m2=closed['surface_area_m2']-cut['surface_area_m2'],
                            joint_cap_area_m2=cap_area,
                            capped_area_added_fraction=(closed['surface_area_m2']-cut['surface_area_m2'])/cut['surface_area_m2'],
                            unscaled=cut,scaled=closed))
    if not records:raise ValueError('no admissible skin surfaces')
    material=skin_layers(mechanics,exterior_area)
    report=dict(schema='ihm.segment-contact-meshes.v1',layer='skin',
                model=str((ROOT/reference_path)),scope='skin_exterior_component',
                model_sha256=sha(ROOT/reference_path),
                geometry_dirs=[str(mesh_path.parent)],
                bodies=sorted(r['body'] for r in records),meshes=len(records),
                total_faces=sum(r['scaled']['faces'] for r in records),
                total_vertices=sum(r['scaled']['vertices'] for r in records),
                watertight_meshes=sum(1 for r in records if r['scaled']['watertight']),
                concave_meshes=sum(1 for r in records if r['scaled']['volume_over_hull_volume']<.999),
                repaired_meshes=sum(1 for r in records if r['simtk_precondition_repaired']),
                refused_meshes=len(rejected),refused=rejected,
                faces_by_body={r['body']:r['scaled']['faces'] for r in records},
                exterior_triangles=int(len(exterior)),
                exterior_surface_area_m2=exterior_area,
                admitted_cut_area_m2=sum(r['cut_surface_area_m2'] for r in records),
                admitted_capped_area_m2=sum(r['capped_surface_area_m2'] for r in records),
                capped_area_added_m2=sum(r['capped_area_added_m2'] for r in records),
                skin_material=material,
                partition='continuous_surface_binding graph-diffused skinning weights; a triangle goes to the argmax of its three vertices\' mean weight. Hard partition of a continuous surface: every segment boundary is a seam the real body does not have.',
                registration=registration_report,
                bone_vertices_inside_skin=float(np.mean([r['bone_vertices_inside_skin'] for r in records if r['bone_vertices_inside_skin'] is not None])) if any(r['bone_vertices_inside_skin'] is not None for r in records) else None,
                segments_enclosing_their_bone=sum(1 for r in records if (r['bone_vertices_inside_skin'] or 0)>=.99),
                enclosure_gate='Share of a segment\'s own bone-mesh vertices lying inside ITS OWN skin piece. This is the PARTITION-CONFOUNDED measure and must not be read as registration quality: a bone that crosses into a neighbour\'s piece, or a segment whose piece is a strip or a patch, scores low under ANY registration. Measured on this bundle: radius 0.07 and ulna 0.12-0.29 own strips, patella 0.10-0.15 owns a patch, the calcaneus sits inside the TOES piece (0.13) and the femoral head inside the PELVIS piece (0.56-0.64) -- while the same bones against the WHOLE skin read 1.000, 1.000, 1.000, 0.91 and 1.000. Segments that own a closed region containing their own bone agree between the two measures (hand 0.78 vs 0.79, toes 0.85 vs 0.85). The partition-free measure is scripts/measure_skin_enclosure_whole.py and it is the one a registration is judged by.',
                segments_enclosing_their_bone_basis='COUNT OF SEGMENTS REACHING 0.99 on the partition-confounded measure above -- NOT a count of segments whose bone is inside their skin. It is 0 for every registration this repo has tried, including the body\'s own anatomical bones, because every bone cloud has vertices near a joint boundary that fall outside its own piece. Do not read 0 as "no segment encloses its bone": the same bundle\'s per-segment values run to 0.89, and against the whole skin 12 of 22 segments reach 0.99.',
                reference_pose=str(reference_path),
                reference_pose_basis='Segment-local stations are taken through the reference run\'s t=0 body transforms, which are the model\'s zero-coordinate neutral pose (only pelvis_ty is nonzero). The supine environment rotates GRAVITY, not the body, so the same transforms serve upright.',
                basis='Canonical exterior skin surface cut per segment and capped so SimTK will accept it. The skin is carried RIGIDLY by its segment: no in-plane stretch, no sliding, no deformable continuum anywhere in this engine. The only compliance is the elastic foundation\'s normal layer.',
                simtk_precondition='SimTK::ContactGeometry::TriangleMesh requires a closed, consistently oriented, non-degenerate edge-2-manifold.',
                controls=controls(),records=records)
    (out_dir/'manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    return report

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',required=True)
    parser.add_argument('--reference',default=DEFAULT_REFERENCE)
    parser.add_argument('--minimum-faces',type=int,default=64)
    parser.add_argument('--registration',choices=('canonical','binding'),default='binding')
    parser.add_argument('--warp',default=None,help='a scripts/skin_warp.py warp (.npz, path relative to the repo) applied to the whole skin on top of the binding map before it is cut; default: none')
    args=parser.parse_args()
    report=build(args.out,args.reference,args.minimum_faces,args.registration,args.warp)
    print(json.dumps({k:v for k,v in report.items() if k not in ('records','controls')},indent=2))
