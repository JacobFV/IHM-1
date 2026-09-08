"""Resolve display-scene aliases, triage the residue, promote every promotable structure.

CANDIDATE. Writes only under --output. Nothing in data/derived/app, data/derived/canonical or
data/raw is read-write; every input is hashed and every hash is recorded.

Run with the isolated libigl environment, not the physiological runtime:
  data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/build_display_promotion_candidate.py --self-test
  data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/build_display_promotion_candidate.py \
      --output data/derived/display-promotion-candidate-v1 --workers 6

Stages, each with its own receipt file:
  1 aliases      every display id that is not a canonical id is tested against the canonical set by
                 id ('body-'+id), by normalized name, and by geometry. Geometry agreement is bbox
                 IoU plus symmetric mean surface distance normalised by the bbox diagonal; FMA
                 concept sets are compared where both records carry them. A name match whose
                 geometry disagrees is reported as a naming collision, not a duplicate.
  2 triage       non-aliases become promote / keep_model_projection / delete, with a reason.
  3 containment  winding-number vertex fraction of each candidate inside every canonical repaired
                 surface whose bounding box could hold it. This is the part-of relation the
                 cross-structure ownership rule lacks; it is written as a first-class field.
  4 promote      registration (axis rotation then the stored z_anatomy affine+TPS), the muscular
                 repair order imported verbatim, role, system, material by role, centroid, bounds,
                 mass where a density exists, and full provenance.
  5 tetgen       the volume-gated forked harness from verify_muscle_tet_ready_surfaces, per role.
No hole is filled and no open sheet is closed; those entities are flagged and left open.
"""
from pathlib import Path
import argparse
import gzip
import importlib.metadata
import json
import multiprocessing as mp
import os
import shutil
import signal
import subprocess
import sys
import tempfile
import time

import numpy as np
import igl

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'scripts'))
from ihm.assembly.anatomy import (FRAME,MODEL_ID,ROTATION,SURFACE_REGION_SYSTEMS,SURFACE_REGION_TOKENS,
                                  LandmarkRegistration,mesh_properties,name_key,normalized_name,
                                  physical_role,refine_role)
from build_muscle_tet_ready_surfaces import (sha,write,signed_volume,weld_exact,diagnose,repair)
from build_entity_tet_ready_surfaces import shape_metrics,analyse_patches
from verify_muscle_tet_ready_surfaces import tetrahedralize,VOLUME_GATE

ANNIHILATION=ROOT/'data/derived/conflict-free-atlas-v1/annihilation.json'
ENVELOPE=ROOT/'data/derived/outer-envelope/outer-envelope.npz'
PROVENANCE_SCHEMA=ROOT/'data/derived/structure-provenance-candidate-v1/schema.json'
MANIFEST=ROOT/'data/derived/app/manifest.json'
ANATOMY=ROOT/'data/derived/canonical/anatomy.json'
ZA_INDEX=ROOT/'data/derived/anatomy/extended/source_index.json'
MATERIALS=ROOT/'data/derived/tissue-material-candidate-v1/materials.json'
MATERIAL_MANIFEST=ROOT/'data/derived/tissue-material-candidate-v1/manifest.json'
PARENT_BUILDS=(ROOT/'data/derived/entity-tet-ready-v1',ROOT/'data/derived/muscle-tet-ready-v1')
RECIPE_SCRIPTS=(ROOT/'scripts/build_muscle_tet_ready_surfaces.py',
                ROOT/'scripts/verify_muscle_tet_ready_surfaces.py',
                ROOT/'scripts/build_entity_tet_ready_surfaces.py')

# role -> (key in materials.json linear block, why that key). Roles absent here have no defensible
# density in the candidate table and their mass stays null rather than being silently filled.
DENSITY_KEY={'muscle':'density','rigid_bone':'density_whole_skeleton','vascular':'density',
             'soft_organ':'density_soft_tissue','skin_layer':'density','skin':'density',
             'adipose':'density'}
# The topographic surface-region rule and the cross-source name key live in ihm.assembly.anatomy so
# the canonical assembly and this lane cannot drift apart on what a role or a name match is.
PART_TOKENS=(' head of ',' part of ',' belly of ',' branch of ',' layer of ',' portion of ',
             ' fibres of ',' fibers of ',' bundle of ',' root of ',' segment of ')
FOREIGN_SPECIMEN_MODELS={'opensim-rajagopal':'Rajagopal2016 published model; the source record says '
                         '"not the BodyParts3D subject". Bone meshes and wrapped muscle path '
                         'polylines from a different specimen in an unregistered frame.',
                         'vascular-aorta':'pediatric Fontan patient aorta; a different human.',
                         'vascular-cerebral':'Vascular Model Repository case; a different human.',
                         'vascular-pulmonary':'Vascular Model Repository case; a different human.'}
NON_ANATOMY_MODELS={'betse-tissue':'BETSE cell-scale bioelectric simulation output; a computational '
                    'lattice, not a structure of this body.'}


def load_json(path):
    return json.loads(Path(path).read_text())


def read_gz(path):
    return json.loads(gzip.decompress(Path(path).read_bytes()))


def write_gz(path,payload):
    Path(path).write_bytes(gzip.compress(json.dumps(payload,allow_nan=False).encode(),mtime=0))


def surface_from_gz(path):
    raw=read_gz(path)
    v=np.asarray(raw['positions'],float).reshape(-1,3)
    f=np.ascontiguousarray(np.asarray(raw['indices'],np.int64).reshape(-1,3))
    return v,f


def bbox(v):
    return v.min(0),v.max(0)


def bbox_iou(a_lo,a_hi,b_lo,b_hi):
    lo=np.maximum(a_lo,b_lo);hi=np.minimum(a_hi,b_hi)
    inter=float(np.prod(np.clip(hi-lo,0,None)))
    if inter<=0:return 0.0
    va=float(np.prod(a_hi-a_lo));vb=float(np.prod(b_hi-b_lo))
    union=va+vb-inter
    return inter/union if union>0 else 0.0


def sample_rows(a,count):
    if len(a)<=count:return a
    return a[np.linspace(0,len(a)-1,count).astype(int)]


def surface_distance(av,af,bv,bf,samples=400):
    """Symmetric mean point-to-surface distance in metres, sampled at vertices both ways."""
    d1,_,_=igl.point_mesh_squared_distance(np.ascontiguousarray(sample_rows(av,samples)),bv,bf)
    d2,_,_=igl.point_mesh_squared_distance(np.ascontiguousarray(sample_rows(bv,samples)),av,af)
    return float((np.sqrt(np.asarray(d1,float)).mean()+np.sqrt(np.asarray(d2,float)).mean())/2)


def geometry_agreement(av,af,bv,bf):
    a_lo,a_hi=bbox(av);b_lo,b_hi=bbox(bv)
    diag=float(np.linalg.norm(a_hi-a_lo))
    msd=surface_distance(av,af,bv,bf)
    va=abs(signed_volume(av,af));vb=abs(signed_volume(bv,bf))
    return dict(bbox_iou=bbox_iou(a_lo,a_hi,b_lo,b_hi),
                centroid_distance_m=float(np.linalg.norm((a_lo+a_hi)/2-(b_lo+b_hi)/2)),
                candidate_bbox_diagonal_m=diag,
                symmetric_mean_surface_distance_m=msd,
                surface_distance_over_diagonal=msd/diag if diag>0 else None,
                signed_volume_ratio=va/vb if vb>0 else None)


# A pair is the same structure when the boxes agree and the surfaces sit within a tenth of the
# candidate's own size of each other. Both bounds are geometric, not anatomical.
ALIAS_IOU=.5
ALIAS_DISTANCE_FRACTION=.10
DUPLICATE_IOU=.8
DUPLICATE_DISTANCE_FRACTION=.02
# A name match whose surfaces disagree is still the same named structure unless the two sit in
# different places. 0.25 m is a quarter of the standing half-height of this body.
COLLISION_CENTROID_DISTANCE_M=.25
TETGEN_TIMEOUT_S=180


def alias_verdict(agreement):
    r=agreement['surface_distance_over_diagonal']
    return bool(agreement['bbox_iou']>=ALIAS_IOU and r is not None and r<=ALIAS_DISTANCE_FRACTION)


def duplicate_verdict(agreement):
    r=agreement['surface_distance_over_diagonal']
    return bool(agreement['bbox_iou']>=DUPLICATE_IOU and r is not None and r<=DUPLICATE_DISTANCE_FRACTION)


def name_part_hint(name):
    low=' '+name.lower().replace('(',' ').replace(')',' ')+' '
    for token in PART_TOKENS:
        if token in low:
            return dict(token=token.strip(),parent_phrase=low.split(token,1)[1].strip().rstrip('.lr').strip())
    return None


def material_for(role,materials):
    entry=materials['materials'].get(role)
    if entry is None:
        return dict(role=role,tissue=None,material_row=None,density_key=None,density_kg_m3=None,
                    density_tier='absent',density_source=None,
                    density_note='no candidate material row exists for this role in '
                                 'data/derived/tissue-material-candidate-v1/materials.json',
                    linear={})
    linear=entry.get('linear',{})
    key=DENSITY_KEY.get(role)
    cell=linear.get(key) if key else None
    density=cell.get('value') if cell else None
    return dict(role=role,tissue=entry.get('tissue'),material_row='materials.json#materials.'+role,
                density_key=key,density_kg_m3=density,
                density_tier=cell.get('tier') if cell else 'absent',
                density_source=cell.get('source') if cell else None,
                density_note=cell.get('note') if cell else
                'no density key selected for this role; mass is deliberately null',
                linear=linear)


def git_build_record(script):
    def run(*args):
        try:
            return subprocess.run(('git',)+args,cwd=str(ROOT),capture_output=True,text=True,
                                  timeout=30).stdout.strip()
        except Exception:
            return ''
    commit=run('rev-parse','HEAD') or None
    dirty=bool(run('status','--porcelain'))
    return dict(script=str(Path(script).resolve().relative_to(ROOT)),script_sha256=sha(script),
                commit=commit,commit_covers_working_tree=not dirty,uncommitted_changes=dirty)


def structure_provenance(structure,source,geometry,transforms,assumptions,build_record,
                         source_vertices,source_faces):
    """One record in the shape data/derived/structure-provenance-candidate-v1/schema.json declares,
    so the universal provenance lane can index a promoted structure without a special case."""
    return dict(schema='ihm.structure-provenance.v1',
        structure_id=structure['display_id'],model_id=MODEL_ID,
        canonical_entity_id=structure['entity_id'],name=structure['name'],
        system=structure['system'],role=structure['role'],evidence_kind='registered_geometry',
        dataset=dict(id='z-anatomy',label=source['label'],version=None,revision=source.get('revision'),
                     url=source['url'],specimen=source['specimen'],units=source['units'],
                     frame=source['frame'],attribution=source['attribution'],
                     acquisition_status=source['status'],license=source['license'],
                     license_url=source.get('license_url'),
                     license_evidence=dict(path=str(ZA_INDEX.relative_to(ROOT)),
                                           field='structures[].source.license'),
                     license_absent_reason=None),
        source_file=dict(path=source['source_geometry_path'],sha256=source['source_geometry_sha256'],
                         sha256_verified=source['sha256_verified'],bytes=source['bytes'],
                         absent_reason=None),
        build=build_record,
        geometry=dict(path=geometry['path'],sha256=geometry['sha256'],sha256_verified=True,
                      representation='triangular_surface',frame=FRAME,units='m',
                      source_vertex_count=source_vertices,source_face_count=source_faces),
        transforms=transforms,
        transform_note='axis permutation into the display frame, then the stored z_anatomy '
                       'landmark registration; the repair step applies to the tet-ready surface only',
        frame_relation='registered_into_canonical',
        tier='transferred',
        tier_basis='geometry from a separate source placed into this body by a fitted transform',
        tier_evidence={'source.specimen':source['specimen'],'source.status':source['status'],
                       'registration.held_out_rms_m':None},
        assumptions=assumptions,derived_artifacts=[],display_present=True)


def registration_from_anatomy(anatomy):
    reg=anatomy['registrations']['z_anatomy']
    source=np.array([l['source_rotated_m'] for l in reg['landmarks']],float)
    target=np.array([l['target_m'] for l in reg['landmarks']],float)
    return LandmarkRegistration(source,target),reg


def register(vertices,rotation,registration):
    return registration.transform(np.asarray(vertices,float)@rotation.T)


# ---------------------------------------------------------------- containment

def containment(candidate_v,parents,fraction=.9,partial=.25,box_overlap=.35,samples=64,max_parents=60):
    """Winding-number vertex fraction of the candidate inside each plausible container.

    parents: list of dicts with entity_id, volume, lo, hi and a zero-argument surface loader.
    Containers are tried smallest first so the tightest one is found before the body wall.
    """
    lo,hi=bbox(candidate_v)
    box=float(np.prod(np.maximum(hi-lo,1e-9)))
    tried=[]
    for r in parents:
        overlap=np.clip(np.minimum(hi,r['hi'])-np.maximum(lo,r['lo']),0,None)
        if float(np.prod(overlap))>=box_overlap*box and float(np.prod(r['hi']-r['lo']))>=box:
            tried.append(r)
    tried.sort(key=lambda r:r['volume'])
    points=np.ascontiguousarray(sample_rows(candidate_v,samples))
    out=[]
    for record in tried[:max_parents]:
        pv,pf=record['load']()
        w=np.asarray(igl.winding_number(pv,pf,points),float)
        f=float((np.abs(w)>.5).mean())
        if f>=partial:
            out.append(dict(parent_entity_id=record['entity_id'],parent_name=record['name'],
                            parent_role=record['role'],parent_volume_m3=record['volume'],
                            enclosed_vertex_fraction=f,
                            relation='part_of' if f>=fraction else 'overlaps',
                            method='winding_number at %d sampled candidate vertices against the '
                                   'repaired canonical surface'%len(points)))
    out.sort(key=lambda r:(-r['enclosed_vertex_fraction'],r['parent_volume_m3']))
    return dict(candidates_tested=len(tried[:max_parents]),bbox_containers=len(tried),
                relations=out,
                part_of=[r['parent_entity_id'] for r in out if r['relation']=='part_of'],
                overlaps=[r['parent_entity_id'] for r in out if r['relation']=='overlaps'],
                tightest_container=out[0]['parent_entity_id'] if out and out[0]['relation']=='part_of' else None)


# ---------------------------------------------------------------- workers

_STATE={}


def _init(state):
    _STATE.update(state)
    _STATE['registration']=LandmarkRegistration(np.asarray(state['reg_source']),
                                                np.asarray(state['reg_target']))


def _parents():
    if 'parents' in _STATE:return _STATE['parents']
    parents=[]
    for record in _STATE['parent_index']:
        path=record['path']
        parents.append(dict(entity_id=record['entity_id'],name=record['name'],role=record['role'],
                            volume=record['volume'],lo=np.asarray(record['lo']),hi=np.asarray(record['hi']),
                            load=(lambda p=path:surface_from_gz(p))))
    _STATE['parents']=parents
    return parents


def _promote_one(job):
    ident,source_id,name,system,path,digest,out=job['id'],job['source_id'],job['name'],job['system'],job['path'],job['sha256'],Path(job['out'])
    try:
        with np.load(path) as mesh:
            raw_v=np.asarray(mesh['vertices'],float);faces=np.ascontiguousarray(np.asarray(mesh['faces'],np.int64))
        if not np.isfinite(raw_v).all() or not len(faces) or faces.min()<0 or faces.max()>=len(raw_v):
            raise ValueError('Invalid source geometry: '+source_id)
        v=_STATE['registration'].transform(raw_v@np.asarray(_STATE['rotation']).T)
        reference=out/'geometry'/(ident+'.json.gz')
        write_gz(reference,dict(positions=[float(x) for x in v.ravel()],
                                indices=[int(i) for i in faces.ravel()],units='m',frame=FRAME,
                                registration_id='z_anatomy',source_geometry_sha256=digest))
        raw_diag=diagnose(v,faces,intersections=False)
        wv,wf=weld_exact(v,faces)
        before=diagnose(wv,wf)
        before_patches=analyse_patches(wv,wf)
        rv,rf,log=repair(v,faces)
        after=diagnose(rv,rf)
        after_patches=analyse_patches(rv,rf,seed_points=True)
        tet=out/'tet-ready'/(ident+'.json.gz')
        write_gz(tet,dict(schema='ihm.promoted-tet-ready-surface.v1',entity_id=ident,name=name,
                          system=system,units='m',frame=FRAME,representation='triangular_surface',
                          source_path=str(Path(path).relative_to(ROOT)),source_sha256=digest,
                          positions=[float(x) for x in rv.ravel()],indices=[int(i) for i in rf.ravel()]))
        drift=after['signed_volume_m3']-before['signed_volume_m3']
        contain=containment(rv if after['closed'] else wv,_parents())
        properties=mesh_properties(v,faces)
        return dict(entity_id=ident,source_id=source_id,name=name,system=system,
                    reference_path=str(reference.relative_to(out)),reference_sha256=sha(reference),
                    tet_ready_path=str(tet.relative_to(out)),tet_ready_sha256=sha(tet),
                    mesh_properties=properties,raw=raw_diag,before=before,after=after,operations=log,
                    before_patch_analysis=before_patches,after_patch_analysis=after_patches,
                    shape_after=shape_metrics(rv,rf),containment=contain,
                    signed_volume_drift_m3=drift,
                    signed_volume_relative_drift=drift/before['signed_volume_m3'] if before['signed_volume_m3'] else None,
                    hole_filling_required=after['boundary_edges']>0,
                    open_sheet=before['boundary_edges']>0 and after['boundary_edges']>0,
                    vertex_links_manifold=after['nonmanifold_vertices']==0,
                    cavity_convention_pending=bool(after_patches.get('void_patches') or
                                                   before_patches.get('authored_cavities') or
                                                   not after_patches.get('analysed',True)),
                    tet_ready=after['closed'] and after['self_intersection_free'] and
                              after['orientation_consistent'] and after['signed_volume_m3']>0)
    except BaseException as error:
        return dict(entity_id=ident,source_id=source_id,name=name,system=system,failed=True,
                    error_type=type(error).__name__,error=str(error))


def _agreement_one(job):
    av,af=surface_from_gz(job['a'])
    if job.get('register'):
        av=_STATE['registration'].transform(av@np.asarray(_STATE['rotation']).T)
    if job.get('a_npz'):
        with np.load(job['a']) as mesh:
            av=np.asarray(mesh['vertices'],float);af=np.ascontiguousarray(np.asarray(mesh['faces'],np.int64))
        av=_STATE['registration'].transform(av@np.asarray(_STATE['rotation']).T)
    bv,bf=surface_from_gz(job['b'])
    return dict(key=job['key'],**geometry_agreement(av,af,bv,bf))


def _npz_agreement(job):
    with np.load(job['a']) as mesh:
        av=np.asarray(mesh['vertices'],float);af=np.ascontiguousarray(np.asarray(mesh['faces'],np.int64))
    av=_STATE['registration'].transform(av@np.asarray(_STATE['rotation']).T)
    bv,bf=surface_from_gz(job['b'])
    return dict(key=job['key'],display_id=job['display_id'],canonical_id=job['canonical_id'],
                **geometry_agreement(av,af,bv,bf))


# ---------------------------------------------------------------- build

class _TetgenTimeout(Exception):
    pass


def _own_children():
    mine=str(os.getpid())
    out=[]
    for entry in Path('/proc').iterdir():
        if not entry.name.isdigit():continue
        try:
            fields=(entry/'stat').read_text().rsplit(')',1)[1].split()
        except OSError:
            continue
        if fields[1]==mine:out.append(int(entry.name))
    return out


def tetgen_with_timeout(v,f,flags,seconds):
    """The imported harness forks and blocks on the child's pipe with no deadline. A 207k-face open
    retinal-artery sheet ran TetGen past 25 minutes at full CPU, so a deadline is imposed here: the
    alarm interrupts the parent's read, the forked child is killed and reaped, and the timeout is
    recorded as a failure exactly like an abort."""
    def handler(signum,frame):
        raise _TetgenTimeout()
    before=set(_own_children())
    previous=signal.signal(signal.SIGALRM,handler);signal.alarm(seconds)
    try:
        outcome=tetrahedralize(v,f,flags)
        outcome['timeout_s']=seconds
        return outcome
    except _TetgenTimeout:
        for pid in set(_own_children())-before:
            try:
                os.kill(pid,signal.SIGKILL);os.waitpid(pid,0)
            except (ProcessLookupError,ChildProcessError):
                pass
        return dict(flags=flags,status=None,exception='Timeout',timeout_s=seconds,
                    message='TetGen did not return within %d s; the forked child was killed'%seconds,
                    tetgen_stdout=[],volume_gate=VOLUME_GATE,
                    tet_volume_vs_surface_relative_error=None,succeeded=False,seconds=float(seconds))
    finally:
        signal.alarm(0);signal.signal(signal.SIGALRM,previous)


def parent_index():
    records=[]
    for build in PARENT_BUILDS:
        for line in (build/'entities.jsonl').read_text().splitlines():
            r=json.loads(line)
            if r.get('failed') or not r.get('tet_ready'):continue
            path=build/r['output_path']
            v,_=surface_from_gz(path)
            lo,hi=bbox(v)
            records.append(dict(entity_id=r['entity_id'],name=r['name'],role=r.get('role','muscle'),
                                volume=abs(r['after']['signed_volume_m3']),path=str(path),
                                lo=lo.tolist(),hi=hi.tolist()))
    return records


def build(output,workers,limit,flags,skip_tetgen,timeout):
    out=Path(output).resolve()
    if out.exists():raise ValueError('Choose a fresh output directory')
    packages={n:importlib.metadata.version(n) for n in ('libigl','numpy','scipy')}
    if packages['libigl']!='2.6.2':raise ValueError('This build pins libigl 2.6.2')
    started=time.monotonic()
    manifest=load_json(MANIFEST);anatomy=load_json(ANATOMY);za_index=load_json(ZA_INDEX)
    materials=load_json(MATERIALS)
    entities={e['id']:e for e in anatomy['entities']}
    by_name={}
    for e in anatomy['entities']:by_name.setdefault(normalized_name(e['name']),[]).append(e['id'])
    structures={s['id']:s for s in manifest['structures']}
    za={m['id']:m for m in za_index['meshes']}
    registration,reg_report=registration_from_anatomy(anatomy)
    out.mkdir(parents=True);(out/'geometry').mkdir();(out/'tet-ready').mkdir();(out/'inputs').mkdir()
    shutil.copyfile(__file__,out/'inputs'/Path(__file__).name)
    print('parent index...',flush=True)
    parents=parent_index()
    state=dict(rotation=ROTATION.tolist(),reg_source=[l['source_rotated_m'] for l in reg_report['landmarks']],
               reg_target=[l['target_m'] for l in reg_report['landmarks']],parent_index=parents)
    pool=mp.get_context('fork').Pool(workers,initializer=_init,initargs=(state,))
    _init(state)

    by_key={}
    for e in anatomy['entities']:by_key.setdefault(name_key(e['name']),[]).append(e['id'])

    # ---- stage 1: aliases -------------------------------------------------
    non_canonical=[s for s in manifest['structures'] if s['id'] not in entities]
    # A display id whose canonical form was collapsed into a duplicate survivor is still an alias of
    # that survivor, not a structure this body lacks. Without this the five collapsed BodyParts3D
    # rows would fall through the id test, fail the name test, and be deleted as an unclassified
    # model on the next run of a lane whose whole job is to not lose anatomy.
    collapsed={d['id']:d['survivor'] for d in anatomy.get('duplicate_surface_collapse',{}).get('dropped',[])}
    aliases=[];residue=[];collisions=[];extent=[]
    name_jobs=[]
    for s in non_canonical:
        cid='body-'+s['id']
        if cid in collapsed:
            survivor=entities[collapsed[cid]]
            aliases.append(dict(display_id=s['id'],display_model=s['model_id'],name=s['name'],
                                surviving_id=survivor['id'],tier='collapsed_duplicate_alias',
                                evidence='the canonical row for this display id was collapsed into an '
                                         'identical duplicate-authored surface; the survivor is the '
                                         'same structure, vertex for vertex',
                                canonical_name=survivor['name'],
                                collapse_evidence=anatomy['duplicate_surface_collapse']['evidence'],
                                concept_sets_equal=set(c['concept_id'] for c in s.get('concepts',[]))==
                                                   set(c['concept_id'] for c in survivor.get('concepts',[])),
                                display_concepts=len(s.get('concepts',[])),
                                canonical_concepts=len(survivor.get('concepts',[])),verified=True))
            continue
        if cid in entities:
            e=entities[cid]
            record=dict(display_id=s['id'],display_model=s['model_id'],name=s['name'],
                        surviving_id=cid,tier='id_alias',verified=None,
                        concept_sets_equal=set(c['concept_id'] for c in s.get('concepts',[]))==
                                           set(c['concept_id'] for c in e.get('concepts',[])),
                        display_concepts=len(s.get('concepts',[])),canonical_concepts=len(e.get('concepts',[])))
            if s['model_id']=='bodyparts3d':
                dp=ROOT/'data/derived/app/geometry'/(s['id']+'.json.gz')
                record.update(evidence='byte_identical_geometry',
                              display_geometry_sha256=sha(dp),
                              canonical_geometry_sha256=e['reference_geometry']['sha256'],
                              verified=sha(dp)==e['reference_geometry']['sha256'])
            elif s['model_id']=='z-anatomy':
                z=za[s['id']]
                with np.load(ROOT/z['source_geometry_path']) as mesh:
                    v=np.asarray(mesh['vertices'],float);f=np.asarray(mesh['faces'],np.int64)
                t=register(v,ROTATION,registration)
                cv,cf=surface_from_gz(ROOT/e['reference_geometry']['path'])
                delta=float(np.abs(t-cv).max()) if cv.shape==t.shape else None
                record.update(evidence='registration_reproduces_canonical_vertices',
                              max_abs_vertex_delta_m=delta,
                              faces_identical=bool(cf.shape==f.shape and np.array_equal(cf,f)),
                              verified=bool(delta is not None and delta<1e-12 and np.array_equal(cf,f)))
            else:
                record.update(evidence='identity_by_construction; non-surface structure',verified=True)
            aliases.append(record)
            continue
        key=normalized_name(s['name'])
        token=name_key(s['name'])
        if s['model_id']=='z-anatomy' and (key in by_name or token in by_key):
            tier='name_alias' if key in by_name else 'token_name_alias'
            target=(by_name[key] if key in by_name else by_key[token])[0]
            name_jobs.append(dict(key=key,tier=tier,display_id=s['id'],canonical_id=target,
                                  a=str(ROOT/za[s['id']]['source_geometry_path']),
                                  b=str(ROOT/entities[target]['reference_geometry']['path'])))
        else:
            residue.append(s)
    print('name-alias geometry: %d pairs'%len(name_jobs),flush=True)
    tiers={j['display_id']:j['tier'] for j in name_jobs}
    for agreement in pool.imap_unordered(_npz_agreement,name_jobs,chunksize=8):
        s=structures[agreement['display_id']]
        confirmed=alias_verdict(agreement)
        record=dict(display_id=s['id'],display_model=s['model_id'],name=s['name'],
                    surviving_id=agreement['canonical_id'],tier=tiers[s['id']],
                    evidence='name-token equality plus geometry agreement',
                    canonical_name=entities[agreement['canonical_id']]['name'],
                    agreement={k:v for k,v in agreement.items() if k not in ('key','display_id','canonical_id')},
                    concept_sets_equal=None,verified=confirmed)
        if confirmed:
            aliases.append(record)
        elif agreement['centroid_distance_m']>COLLISION_CENTROID_DISTANCE_M:
            record['tier']='naming_collision'
            record['collided_with']=agreement['canonical_id']
            record['surviving_id']=None
            collisions.append(record);residue.append(s)
        else:
            record['tier']='name_alias_extent_mismatch'
            record['evidence']=('name-token equality; the two authorings of the same named structure '
                                'disagree in extent. Same structure, cross-source geometry conflict.')
            record['verified']=True
            aliases.append(record);extent.append(record)

    # ---- stage 2: triage --------------------------------------------------
    triage=[];promote=[]
    for s in residue:
        model=s['model_id']
        base=dict(display_id=s['id'],display_model=model,name=s['name'],system=s.get('system'),
                  kind=s.get('kind'))
        if model in FOREIGN_SPECIMEN_MODELS:
            triage.append(dict(**base,decision='delete',category='foreign_specimen',
                               reason=FOREIGN_SPECIMEN_MODELS[model]))
        elif model in NON_ANATOMY_MODELS:
            triage.append(dict(**base,decision='delete',category='not_anatomy',
                               reason=NON_ANATOMY_MODELS[model]))
        elif model==MODEL_ID:
            triage.append(dict(**base,decision='keep',category='model_projection',
                               reason='derived from this model by another lane; it is a projection '
                                      'of the implicit human model and is not a display artifact. '
                                      'No triangular source surface, so not promoted here.'))
        elif model=='z-anatomy':
            z=za[s['id']]
            collision=next((c for c in collisions if c['display_id']==s['id']),None)
            role,refinement=refine_role(s['name'],s['system'],physical_role(s['name'],s['system']))
            hint=name_part_hint(s['name'])
            entry=dict(**base,decision='promote',category='surface_region' if role=='surface_region' else 'anatomical_structure',
                       role=role,role_refinement=refinement,name_part_hint=hint,
                       source_vertices=z['source_vertices'],source_triangles=z['source_triangles'],
                       source_collections=z['source_collections'],
                       name_collision_with=collision['collided_with'] if collision else None,
                       reason='absent from the canonical model by id and by name'+
                              ('; its name matches canonical %s but the geometry disagrees, so it is '
                               'a naming collision and a distinct structure'%collision['collided_with']
                               if collision else ''))
            triage.append(entry);promote.append((s,z,entry))
        else:
            triage.append(dict(**base,decision='delete',category='unclassified_model',
                               reason='display model %s is neither the canonical model nor a source '
                                      'the canonical model was built from'%model))
    if limit:promote=promote[:limit]

    # ---- stage 3+4: containment and promotion -----------------------------
    jobs=[dict(id='body-'+s['id'],source_id=s['id'],name=s['name'],system=s['system'],
               path=str(ROOT/z['source_geometry_path']),sha256=z['source_geometry_sha256'],out=str(out))
          for s,z,_ in promote]
    for job in jobs:
        if job['id'] in entities:raise ValueError('Promoted id collides with a canonical id: '+job['id'])
    print('promoting %d structures on %d workers...'%(len(jobs),workers),flush=True)
    built={}
    done=0
    for record in pool.imap_unordered(_promote_one,jobs,chunksize=4):
        built[record['entity_id']]=record;done+=1
        if done%200==0:print('  repaired %d/%d'%(done,len(jobs)),flush=True)
    pool.close();pool.join()

    # ---- stage 4b: geometric alias sweep over the repaired candidates -----
    dup=[]
    parent_lo=np.array([p['lo'] for p in parents]);parent_hi=np.array([p['hi'] for p in parents])
    for ident,record in built.items():
        if record.get('failed'):continue
        b=record['mesh_properties']['bounds_m']
        lo=np.array(b['min']);hi=np.array(b['max'])
        ious=np.array([bbox_iou(lo,hi,parent_lo[i],parent_hi[i]) for i in range(len(parents))])
        order=np.argsort(-ious)[:5]
        av,af=None,None
        for index in order:
            if ious[index]<DUPLICATE_IOU:break
            if av is None:av,af=surface_from_gz(out/record['tet_ready_path'])
            bv,bf=surface_from_gz(parents[int(index)]['path'])
            agreement=geometry_agreement(av,af,bv,bf)
            if duplicate_verdict(agreement):
                dup.append(dict(display_id=record['source_id'],name=record['name'],
                                canonical_id=parents[int(index)]['entity_id'],
                                canonical_name=parents[int(index)]['name'],agreement=agreement,
                                verdict='geometric_duplicate_under_a_different_name'))
                break
    write(out/'geometric-duplicates.json',dup)
    duplicate_ids={d['display_id'] for d in dup}
    for d in dup:
        s=structures[d['display_id']]
        aliases.append(dict(display_id=d['display_id'],display_model=s['model_id'],name=s['name'],
                            surviving_id=d['canonical_id'],tier='geometric_alias',
                            evidence='no name match; the registered surface duplicates a canonical '
                                     'surface within the duplicate thresholds',
                            canonical_name=d['canonical_name'],agreement=d['agreement'],
                            concept_sets_equal=None,verified=True))
        for row in triage:
            if row['display_id']==d['display_id']:
                row.update(decision='collapse',category='geometric_alias',
                           surviving_id=d['canonical_id'],
                           reason='the registered surface duplicates canonical %s (%s) within the '
                                  'duplicate thresholds; promoting it would double the structure'
                                  %(d['canonical_id'],d['canonical_name']))
    with (out/'aliases.jsonl').open('w') as handle:
        for r in sorted(aliases,key=lambda r:r['display_id']):handle.write(json.dumps(r,allow_nan=False)+'\n')
    write(out/'naming-collisions.json',dict(
        definition='a name-token match whose surfaces disagree AND whose bounding-box centres are '
                   'more than %.2f m apart'%COLLISION_CENTROID_DISTANCE_M,
        count=len(collisions),rows=sorted(collisions,key=lambda r:r['display_id'])))
    write(out/'name-alias-extent-mismatch.json',dict(
        definition='the same named structure authored twice at different extent; the canonical id '
                   'survives and the measured disagreement is handed to the cross-structure conflict lane',
        count=len(extent),
        max_centroid_distance_m=max((r['agreement']['centroid_distance_m'] for r in extent),default=0.0),
        rows=sorted(extent,key=lambda r:-r['agreement']['centroid_distance_m'])))
    write(out/'alias-map.json',{r['display_id']:r['surviving_id'] for r in sorted(aliases,key=lambda r:r['display_id'])})

    # ---- stage 5: entity records -----------------------------------------
    za_model=next(m for m in manifest['models'] if m['id']=='z-anatomy')
    recipe_hashes={str(p.relative_to(ROOT)):sha(p) for p in RECIPE_SCRIPTS}
    inputs={str(p.relative_to(ROOT)):sha(p) for p in (MANIFEST,ANATOMY,ZA_INDEX,MATERIALS,MATERIAL_MANIFEST)}
    for build_dir in PARENT_BUILDS:
        inputs[str((build_dir/'manifest.json').relative_to(ROOT))]=sha(build_dir/'manifest.json')
    transform_record=[dict(op='axis_rotation',matrix=ROTATION.tolist(),
                           source_frame='z-anatomy-blender-world',target_frame='bodyparts3d-display-axes'),
                      dict(op='landmark_registration',registration_id='z_anatomy',
                           method=reg_report['method'],affine_4x3=reg_report['affine_4x3'],
                           smoothing_prior=reg_report['smoothing_prior'],
                           landmark_count=reg_report['landmark_count'],
                           fit_rms_m=reg_report['fit_rms_m'],held_out_rms_m=reg_report['held_out_rms_m'],
                           held_out_max_m=reg_report['held_out_max_m'],
                           source='data/derived/canonical/anatomy.json#registrations.z_anatomy',
                           anatomy_sha256=inputs['data/derived/canonical/anatomy.json'])]
    build_record=git_build_record(Path(__file__).resolve())
    provenance_transforms=[
        dict(op='axis_rotation',matrix=ROTATION.tolist(),
             source_frame='z-anatomy-blender-world',target_frame='bodyparts3d-display-axes',
             residual=None,residual_reason='exact orthonormal axis permutation; no fit'),
        dict(op='landmark_registration',registration_id='z_anatomy',method=reg_report['method'],
             affine_4x3=reg_report['affine_4x3'],smoothing_prior=reg_report['smoothing_prior'],
             landmark_count=reg_report['landmark_count'],
             source='data/derived/canonical/anatomy.json#registrations.z_anatomy',
             anatomy_sha256=inputs['data/derived/canonical/anatomy.json'],
             residual=dict(metric='held_out_rms_m',value=reg_report['held_out_rms_m'],
                           method='deterministic five-fold bone-centroid holdout'))]
    promoted=[]
    for s,z,entry in promote:
        ident='body-'+s['id']
        record=built.get(ident)
        if record is None or record.get('failed') or s['id'] in duplicate_ids:continue
        role=entry['role']
        material=material_for(role,materials)
        props=record['mesh_properties']
        closed=record['after']['closed'] and record['after']['signed_volume_m3']>0
        volume=abs(record['after']['signed_volume_m3']) if closed else None
        is_tissue=role!='surface_region'
        density=material['density_kg_m3'] if is_tissue else None
        mass=volume*density if (volume is not None and density) else None
        e=dict(id=ident,model_id=MODEL_ID,source_id=s['id'],name=normalized_name(s['name']),
               display_name=s['name'],system=s['system'],role=role,
               role_basis='ihm.assembly.anatomy.physical_role'+(' plus surface-region refinement'
                          if entry['role_refinement'] else ''),
               role_refinement_note=entry['role_refinement'],
               evidence_kind='registered_geometry',
               reference_geometry=dict(path=str((out/record['reference_path']).relative_to(ROOT)),
                                       sha256=record['reference_sha256'],frame=FRAME,units='m',
                                       representation='triangular_surface'),
               tet_ready_geometry=dict(path=str((out/record['tet_ready_path']).relative_to(ROOT)),
                                       sha256=record['tet_ready_sha256'],frame=FRAME,units='m',
                                       representation='triangular_surface',
                                       closed=record['after']['closed'],
                                       self_intersection_free=record['after']['self_intersection_free'],
                                       hole_filling_required=record['hole_filling_required'],
                                       open_sheet=record['open_sheet']),
               bounds_m=props['bounds_m'],centroid_m=props['centroid_m'],
               centroid_definition=props['centroid_definition'],
               surface_area_m2=record['shape_after']['surface_area_m2'],
               volume_m3=volume,
               volume_method='absolute signed surface integral of the repaired shell; closed, oriented '
                             'and self-intersection free' if closed else
                             'unknown: the repaired surface still carries boundary or nonmanifold edges '
                             'and no hole was filled',
               source_vertex_count=props['source_vertex_count'],source_face_count=props['source_face_count'],
               principal_axis=props['principal_axis'],
               material=dict(**{k:v for k,v in material.items() if k!='linear'},
                             applies_to_mass=is_tissue,
                             material_manifest_sha256=inputs['data/derived/tissue-material-candidate-v1/manifest.json']),
               mass_kg=mass,
               mass_method=('volume times role density from the tissue-material candidate'
                            if mass is not None else
                            'absent: no defensible density for this role, or the shell is not closed, '
                            'or the structure is a surface region carrying no tissue volume'),
               containment=record['containment'],
               connections=[dict(entity_id=r['parent_entity_id'],relation=r['relation'],
                                 enclosed_vertex_fraction=r['enclosed_vertex_fraction'],
                                 evidence=r['method']) for r in record['containment']['relations']],
               name_part_hint=entry['name_part_hint'],
               name_collision_with=entry.get('name_collision_with'),
               repair=dict(operations=record['operations'],before=record['before'],after=record['after'],
                           signed_volume_relative_drift=record['signed_volume_relative_drift'],
                           cavity_convention_pending=record['cavity_convention_pending'],
                           patch_analysis=record['after_patch_analysis']),
               shape=record['shape_after'],
               provenance=dict(
                   source_dataset=dict(id='z-anatomy',label=za_model['source']['label'] if 'source' in za_model else 'Z-Anatomy',
                                       url=s['source']['url'],revision=s['source'].get('revision'),
                                       archive_sha256=s['source']['sha256'],license=s['source']['license'],
                                       attribution=s['source']['attribution'],
                                       specimen=s['source']['specimen']),
                   source_files=[dict(path=z['source_geometry_path'],sha256=z['source_geometry_sha256'],
                                      role='authored source surface, Blender world frame'),
                                 dict(path=str(ZA_INDEX.relative_to(ROOT)),sha256=inputs[str(ZA_INDEX.relative_to(ROOT))],
                                      role='source index'),
                                 dict(path=str(MANIFEST.relative_to(ROOT)),sha256=inputs[str(MANIFEST.relative_to(ROOT))],
                                      role='display scene the structure was carried in')],
                   source_object=z['source_object'],source_collections=z['source_collections'],
                   extraction_scripts=[dict(path='scripts/build_extended_anatomy.py',
                                            sha256=sha(ROOT/'scripts/build_extended_anatomy.py'),
                                            role='acquired and evaluated the source surfaces')],
                   promotion_script=dict(path=str(Path(__file__).resolve().relative_to(ROOT)),
                                         sha256=sha(Path(__file__).resolve())),
                   repair_scripts=recipe_hashes,
                   transforms=transform_record+[dict(op='surface_repair',
                       method='exact weld; drop repeated-index and exactly-collinear faces; collapse '
                              'coincident faces; bfs_orient patches flipped outward; CGAL '
                              'remesh_self_intersections; CGAL self-union where nonmanifold edges remain',
                       script='scripts/build_muscle_tet_ready_surfaces.py',
                       sha256=recipe_hashes['scripts/build_muscle_tet_ready_surfaces.py'],
                       applied_to='tet_ready_geometry only; reference_geometry is the registered '
                                  'source surface with no repair')],
                   outputs=dict(reference_geometry_sha256=record['reference_sha256'],
                                tet_ready_geometry_sha256=record['tet_ready_sha256']),
                   dependency_group='bodyparts3d-derived-reference',
                   source_ids=[s['id']]),
               assumptions=['CANONICAL-GENERIC-REFERENCE','Z-LANDMARK-REGISTRATION']+
                           (['SURFACE-REGION-NO-TISSUE-VOLUME'] if role=='surface_region' else []),
               structure_provenance=structure_provenance(
                   dict(display_id=s['id'],entity_id=ident,name=normalized_name(s['name']),
                        system=s['system'],role=role),
                   dict(s['source'],source_geometry_path=z['source_geometry_path'],
                        source_geometry_sha256=z['source_geometry_sha256'],
                        sha256_verified=sha(ROOT/z['source_geometry_path'])==z['source_geometry_sha256'],
                        bytes=(ROOT/z['source_geometry_path']).stat().st_size),
                   dict(path=str((out/record['reference_path']).relative_to(ROOT)),
                        sha256=record['reference_sha256']),
                   provenance_transforms+[dict(op='surface_repair',
                       script='scripts/build_muscle_tet_ready_surfaces.py',
                       sha256=recipe_hashes['scripts/build_muscle_tet_ready_surfaces.py'],
                       applies_to=str((out/record['tet_ready_path']).relative_to(ROOT)),
                       output_sha256=record['tet_ready_sha256'],
                       residual=dict(metric='signed_volume_relative_drift',
                                     value=record['signed_volume_relative_drift'],
                                     method='repaired minus welded signed surface integral'))],
                   ['CANONICAL-GENERIC-REFERENCE','Z-LANDMARK-REGISTRATION'],build_record,
                   props['source_vertex_count'],props['source_face_count']),
               uncertainty=dict(biological=dict(status='generic authored reference; population '
                                                       'variability not estimated',
                                                confidence_percent=None,independent_subject_count=0,
                                                subject_calibrated=False),
                                registration_or_synthesis=dict(kind='registered_geometry',
                                                               registration_id='z_anatomy',
                                                               held_out_rms_m=reg_report['held_out_rms_m'],
                                                               calibrated_probability=None),
                                display_numerics=dict(surface_only=True,display_reduction_applied=False,
                                                      position_units='m',rounding='JSON float serialization')))
        promoted.append(e)
    with (out/'entities.jsonl').open('w') as handle:
        for e in sorted(promoted,key=lambda e:e['id']):handle.write(json.dumps(e,allow_nan=False)+'\n')
    with (out/'decisions.jsonl').open('w') as handle:
        for r in sorted(triage,key=lambda r:r['display_id']):handle.write(json.dumps(r,allow_nan=False)+'\n')
    failures=[r for r in built.values() if r.get('failed')]
    write(out/'promotion-failures.json',failures)

    # ---- stage 6: tetgen --------------------------------------------------
    trials=[]
    if not skip_tetgen:
        scratch=tempfile.mkdtemp(prefix='promotion-tetgen-')
        here=os.getcwd()
        try:
            os.chdir(scratch)
            for i,e in enumerate(sorted(promoted,key=lambda e:e['id'])):
                rv,rf=surface_from_gz(ROOT/e['tet_ready_geometry']['path'])
                repaired=tetgen_with_timeout(rv,rf,flags,timeout)
                cv,cf=surface_from_gz(ROOT/e['reference_geometry']['path'])
                cv,cf=weld_exact(cv,cf)
                control=tetgen_with_timeout(np.ascontiguousarray(cv),np.ascontiguousarray(cf),flags,timeout)
                trial=dict(entity_id=e['id'],name=e['name'],role=e['role'],system=e['system'],
                           faces_repaired=int(len(rf)),
                           hole_filling_required=e['tet_ready_geometry']['hole_filling_required'],
                           open_sheet=e['tet_ready_geometry']['open_sheet'],
                           welded_control=control,repaired=repaired)
                trials.append(trial)
                with (out/'tetgen-trials.jsonl').open('a') as handle:
                    handle.write(json.dumps(trial,allow_nan=False)+'\n')
                if (i+1)%200==0:print('  tetgen %d/%d'%(i+1,len(promoted)),flush=True)
        finally:
            os.chdir(here);shutil.rmtree(scratch,ignore_errors=True)
        stray=sorted(str(p.relative_to(ROOT)) for p in ROOT.glob('tetgen-tmpfile_skipped.*'))
        write(out/'tetgen-stray-files.json',dict(found_at_repo_root=stray))

    # ---- stage 6b: the part-of relation against the annihilation defect ----
    write(out/'part-of-defect-probe.json',defect_probe(parents))
    envelope_check(out)

    # ---- stage 7: reports -------------------------------------------------
    write(out/'deletions.json',deletion_report(triage,aliases,dup,structures))
    write(out/'provenance-coverage.json',provenance_coverage(anatomy))
    summary=summarise(anatomy,manifest,aliases,collisions,extent,triage,promoted,trials,dup,failures,
                      started,flags,skip_tetgen)
    write(out/'summary.json',summary)
    artifacts={str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file()}
    write(out/'manifest.json',dict(schema='ihm.display-promotion-candidate.v1',
        status='CANDIDATE. Not canonical, not appended to any manifest, not consumed by any runtime.',
        packages=packages,python=sys.version,inputs_sha256=inputs,recipe_sha256=recipe_hashes,
        canonical_assets_modified=False,display_scene_modified=False,artifacts_sha256=artifacts,
        thresholds=dict(alias_bbox_iou=ALIAS_IOU,alias_surface_distance_over_diagonal=ALIAS_DISTANCE_FRACTION,
                        duplicate_bbox_iou=DUPLICATE_IOU,
                        duplicate_surface_distance_over_diagonal=DUPLICATE_DISTANCE_FRACTION,
                        containment_part_of_fraction=.9,containment_overlap_fraction=.25,
                        containment_bbox_overlap_fraction=.35,
                        tetgen_volume_gate=VOLUME_GATE,tetgen_timeout_s=timeout),
        limitations=[
            'Promotion adds Z-Anatomy authored surfaces registered by an affine plus thin-plate-spline '
            'fit to shared bone bounding-box centres. It adds no independent measured subject.',
            'No hole is filled and no open sheet is closed. Open shells are promoted as entities with '
            'volume_m3 null and are reported, not repaired into solids.',
            'Containment is a geometric winding-number relation against repaired canonical surfaces. '
            'It is not an FMA part-of assertion and does not certify tissue ownership.',
            'Mass is computed only for roles carrying a density in the tissue-material candidate. '
            'Roles whose density is absent have mass null by design.',
            'Inter-entity overlap between promoted and canonical surfaces is unresolved here; the '
            'cross-structure conflict lane owns it.']))
    print(json.dumps({k:v for k,v in summary.items() if not isinstance(v,list)},indent=2))
    return summary


def defect_probe(parents):
    """The cross-structure ownership rule deleted 80 entities by handing every disputed tet to a
    priority owner. Each was a part inside its own container. Recompute containment for exactly
    those entities and report how many the relation recovers."""
    if not ANNIHILATION.exists():
        return dict(available=False,reason='data/derived/conflict-free-atlas-v1/annihilation.json absent')
    report=load_json(ANNIHILATION)
    prepared=[dict(entity_id=r['entity_id'],name=r['name'],role=r['role'],volume=r['volume'],
                   lo=np.asarray(r['lo']),hi=np.asarray(r['hi']),path=r['path'],
                   load=(lambda p=r['path']:surface_from_gz(p))) for r in parents]
    by_id={r['entity_id']:r for r in prepared}
    rows=[]
    for victim in report['entities_reduced_to_nothing']:
        parent=by_id.get(victim['entity_id'])
        if parent is None:
            rows.append(dict(entity_id=victim['entity_id'],name=victim['name'],
                             available=False,reason='no repaired surface in the parent index'))
            continue
        v,_=surface_from_gz(parent['path'])
        others=[r for r in prepared if r['entity_id']!=victim['entity_id']]
        found=containment(v,others)
        owners=set(victim['owners'])
        rows.append(dict(entity_id=victim['entity_id'],name=victim['name'],role=victim['role'],
                         system=victim['system'],volume_before_m3=victim['volume_before_m3'],
                         recorded_owners=sorted(owners),
                         containment_part_of=found['part_of'],
                         containment_overlaps=found['overlaps'],
                         max_enclosed_vertex_fraction=max((r['enclosed_vertex_fraction']
                                                           for r in found['relations']),default=0.0),
                         owners_overlapping=sorted(owners&set(found['overlaps'])),
                         tightest_container=found['tightest_container'],
                         owners_recovered=sorted(owners&set(found['part_of'])),
                         recovered=bool(owners&set(found['part_of'])),
                         available=True))
    usable=[r for r in rows if r.get('available')]
    return dict(schema='ihm.part-of-defect-probe.v1',
        source='data/derived/conflict-free-atlas-v1/annihilation.json',
        source_sha256=sha(ANNIHILATION),
        entities_reduced_to_nothing=report['entities_reduced_to_nothing_count'],
        probed=len(usable),
        with_a_recovered_container=sum(r['recovered'] for r in usable),
        with_any_container=sum(bool(r['containment_part_of']) for r in usable),
        with_a_recovered_owner_as_part_or_overlap=sum(
            bool(r['owners_recovered'] or r['owners_overlapping']) for r in usable),
        with_a_partial_container=sum(bool(r['containment_overlaps']) for r in usable),
        finding='the entities the ownership rule annihilated are parts lying inside the structures '
                'that took their volume. The winding-number relation recovers the enclosure, fully '
                'for compact parts and partially for branching parts that leave their container at '
                'a hilum, so a part-of aware rule would concede volume without deleting the part.',
        rows=rows)


def envelope_check(output):
    """Every promoted surface must land inside the watertight body envelope. A registration that
    put a structure outside the body would show here and nowhere else in this build."""
    out=Path(output).resolve()
    with np.load(ENVELOPE) as mesh:
        ev=np.ascontiguousarray(np.asarray(mesh['positions'],float).reshape(-1,3))
        ef=np.ascontiguousarray(np.asarray(mesh['indices'],np.int64).reshape(-1,3))
    rows=[]
    for line in (out/'entities.jsonl').read_text().splitlines():
        e=json.loads(line)
        v,_=surface_from_gz(ROOT/e['tet_ready_geometry']['path'])
        points=np.ascontiguousarray(sample_rows(v,64))
        w=np.asarray(igl.winding_number(ev,ef,points),float)
        f=float((np.abs(w)>.5).mean())
        rows.append(dict(entity_id=e['id'],name=e['name'],system=e['system'],role=e['role'],
                         enclosed_vertex_fraction=f))
    inside=[r for r in rows if r['enclosed_vertex_fraction']>=.9]
    report=dict(schema='ihm.promoted-envelope-containment.v1',
        envelope=str(ENVELOPE.relative_to(ROOT)),envelope_sha256=sha(ENVELOPE),
        envelope_faces=int(len(ef)),promoted=len(rows),
        fully_inside=len(inside),
        partly_outside=[r for r in rows if .0<r['enclosed_vertex_fraction']<.9],
        entirely_outside=[r for r in rows if r['enclosed_vertex_fraction']==0.0],
        median_enclosed_vertex_fraction=float(np.median([r['enclosed_vertex_fraction'] for r in rows])),
        interpretation='a promoted surface not inside the body envelope is a registration failure, '
                       'not an anatomical finding')
    write(out/'envelope-containment.json',report)
    print(json.dumps({k:v for k,v in report.items() if not isinstance(v,list)},indent=2))
    return report


def deletion_report(triage,aliases,duplicates,structures):
    rows=[dict(display_id=r['display_id'],name=r['name'],model=r['display_model'],
               category=r['category'],reason=r['reason'])
          for r in triage if r['decision']=='delete']
    collapse=[dict(display_id=r['display_id'],name=r['name'],model=r['display_model'],
                   category='duplicate_of_canonical_'+r['tier'],
                   surviving_id=r['surviving_id'],
                   reason='the same structure is already carried under the canonical id; the display '
                          'scene lists it twice')
              for r in aliases]
    return dict(schema='ihm.display-deletion-list.v1',
                statement='A list only. Nothing is deleted here and data/derived/app is regenerable '
                          'from its builder.',
                remove_as_not_a_projection_of_the_model=rows,
                collapse_onto_canonical_id=collapse,
                geometric_duplicates_under_a_different_name=duplicates,
                counts=dict(remove=len(rows),collapse=len(collapse),geometric_duplicates=len(duplicates)),
                downstream_warning='opensim-rajagopal structures are inputs to the mechanics lane. '
                                   'Removing them from the anatomical scene does not remove the model; '
                                   'check ihm/ before the reviewed deletion step.')


REQUIRED_PROVENANCE=('source dataset id','source file path','source file sha256','extraction script',
                     'transform record')


def provenance_coverage(anatomy):
    rows=[]
    for e in anatomy['entities']:
        p=e.get('provenance',{})
        rows.append(dict(entity_id=e['id'],
                         has_source_ids=bool(p.get('source_ids')),
                         has_source_files=bool(p.get('files')),
                         has_source_file_hashes=all('sha256' in f for f in p.get('files',[])) and bool(p.get('files')),
                         has_source_block=bool(p.get('source')),
                         has_extraction_script=False,
                         has_transform_record=bool(p.get('source_to_canonical') or p.get('registration_id')),
                         has_geometry_hash=bool(e.get('reference_geometry',{}).get('sha256'))))
    keys=[k for k in rows[0] if k!='entity_id']
    return dict(schema='ihm.provenance-coverage.v1',
                entities=len(rows),
                required_fields=list(REQUIRED_PROVENANCE),
                coverage={k:sum(r[k] for r in rows) for k in keys},
                gap_note='No canonical entity records the extraction script that produced its source '
                         'surface. Every promoted record in this candidate carries '
                         'provenance.extraction_scripts, provenance.promotion_script, '
                         'provenance.repair_scripts and provenance.transforms so the universal '
                         'provenance lane can index newcomers and back-fill incumbents.',
                rows=rows)


def summarise(anatomy,manifest,aliases,collisions,extent,triage,promoted,trials,duplicates,failures,
              started,flags,skip_tetgen):
    from collections import Counter
    before_entities=len(anatomy['entities'])
    before_systems=Counter(e['system'] for e in anatomy['entities'])
    before_roles=Counter(e['role'] for e in anatomy['entities'])
    add_systems=Counter(e['system'] for e in promoted)
    add_roles=Counter(e['role'] for e in promoted)
    volumes=[e['volume_m3'] for e in promoted if e['volume_m3'] is not None and e['role']!='surface_region']
    masses=[e['mass_kg'] for e in promoted if e['mass_kg'] is not None]
    region_volumes=[e['volume_m3'] for e in promoted if e['volume_m3'] is not None and e['role']=='surface_region']
    no_mass=Counter(e['role'] for e in promoted if e['mass_kg'] is None)
    by_role={}
    for t in trials:
        r=by_role.setdefault(t['role'],dict(trials=0,repaired_ok=0,control_ok=0,open_sheet=0,aborts=0))
        r['trials']+=1;r['repaired_ok']+=bool(t['repaired']['succeeded'])
        r['control_ok']+=bool(t['welded_control']['succeeded']);r['open_sheet']+=bool(t['open_sheet'])
        r['aborts']+=t['repaired'].get('exception')=='ProcessAborted'
    closed=[e for e in promoted if e['volume_m3'] is not None]
    return dict(schema='ihm.display-promotion-candidate.v1',
        display_structures=len(manifest['structures']),canonical_entities_before=before_entities,
        display_ids_not_canonical_ids=len(manifest['structures'])-sum(
            1 for s in manifest['structures'] if s['id'] in {e['id'] for e in anatomy['entities']}),
        aliases_resolved=len(aliases),
        aliases_by_tier=dict(Counter(r['tier'] for r in aliases)),
        aliases_verified=sum(bool(r['verified']) for r in aliases),
        aliases_unverified=[r['display_id'] for r in aliases if not r['verified']],
        naming_collisions=len(collisions),
        name_alias_extent_mismatch=len(extent),
        name_alias_extent_mismatch_max_centroid_distance_m=max(
            (r['agreement']['centroid_distance_m'] for r in extent),default=0.0),
        triaged=len(triage),
        triage_decisions=dict(Counter(r['decision'] for r in triage)),
        triage_categories=dict(Counter(r['category'] for r in triage)),
        promotion_attempted=sum(r['decision']=='promote' for r in triage),
        promoted=len(promoted),promotion_failures=len(failures),
        promoted_by_system=dict(add_systems),promoted_by_role=dict(add_roles),
        entities_after=before_entities+len(promoted),
        systems_before=dict(before_systems),
        systems_after={k:before_systems.get(k,0)+add_systems.get(k,0)
                       for k in set(before_systems)|set(add_systems)},
        roles_before=dict(before_roles),
        roles_after={k:before_roles.get(k,0)+add_roles.get(k,0) for k in set(before_roles)|set(add_roles)},
        closed_shells=len(closed),open_shells=len(promoted)-len(closed),
        open_shells_by_system=dict(Counter(e['system'] for e in promoted if e['volume_m3'] is None)),
        added_tissue_volume_m3=float(sum(volumes)),added_tissue_volume_L=float(sum(volumes)*1000),
        surface_region_enclosed_volume_m3=float(sum(region_volumes)),
        surface_region_note='excluded from added tissue volume: topographic regions lie on the skin '
                            'and their enclosed volume is not additional tissue',
        added_mass_kg=float(sum(masses)),
        entities_with_mass=len(masses),
        entities_without_mass_by_role=dict(no_mass),
        mass_gap_note='roles with no defensible density in the tissue-material candidate carry '
                      'mass null; the largest such role is reported above',
        containment_relations=sum(len(e['containment']['relations']) for e in promoted),
        with_part_of=sum(bool(e['containment']['part_of']) for e in promoted),
        with_tightest_container=sum(e['containment']['tightest_container'] is not None for e in promoted),
        with_name_part_hint=sum(e['name_part_hint'] is not None for e in promoted),
        with_overlap_relations=sum(bool(e['containment']['overlaps']) for e in promoted),
        geometric_duplicates=len(duplicates),
        hole_filling_required=sum(e['tet_ready_geometry']['hole_filling_required'] for e in promoted),
        open_sheets_flagged_not_closed=sum(e['tet_ready_geometry']['open_sheet'] for e in promoted),
        cavity_convention_pending=sum(e['repair']['cavity_convention_pending'] for e in promoted),
        tetgen_flags=flags,tetgen_skipped=skip_tetgen,tetgen_trials=len(trials),
        tetgen_timeout_s=TETGEN_TIMEOUT_S,
        tetgen_timeouts=sum(t['repaired'].get('exception')=='Timeout' for t in trials),
        tetgen_control_timeouts=sum(t['welded_control'].get('exception')=='Timeout' for t in trials),
        tetgen_succeeded=sum(t['repaired']['succeeded'] for t in trials),
        tetgen_control_succeeded=sum(t['welded_control']['succeeded'] for t in trials),
        tetgen_success_rate=sum(t['repaired']['succeeded'] for t in trials)/len(trials) if trials else None,
        tetgen_success_rate_closed_only=(sum(t['repaired']['succeeded'] for t in trials if not t['open_sheet'])/
                                         max(1,sum(1 for t in trials if not t['open_sheet']))) if trials else None,
        tetgen_by_role=by_role,
        tetgen_total_tets=sum(t['repaired'].get('tets',0) for t in trials if t['repaired']['succeeded']),
        tetgen_process_aborts=sum(t['repaired'].get('exception')=='ProcessAborted' for t in trials),
        tetgen_failures=[dict(entity_id=t['entity_id'],name=t['name'],role=t['role'],
                              open_sheet=t['open_sheet'],
                              exception=t['repaired'].get('exception'),
                              status=t['repaired'].get('status'),
                              message=t['repaired'].get('message'),
                              tetgen_stdout=t['repaired']['tetgen_stdout'][-8:])
                         for t in trials if not t['repaired']['succeeded']],
        elapsed_s=time.monotonic()-started)


# ---------------------------------------------------------------- self test

def _cube(scale,offset):
    c=np.array([[0.,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1.],[1,0,1],[1,1,1],[0,1,1]])*scale+offset
    q=np.array([[0,2,1],[0,3,2],[4,5,6],[4,6,7],[0,1,5],[0,5,4],[1,2,6],[1,6,5],[2,3,7],[2,7,6],
                [3,0,4],[3,4,7]],np.int64)
    return c,q


def self_test():
    with tempfile.TemporaryDirectory() as scratch:
        scratch=Path(scratch)
        ov,of=_cube(1.,np.zeros(3))
        same=geometry_agreement(ov,of,ov,of)
        assert same['bbox_iou']==1.0 and same['symmetric_mean_surface_distance_m']==0.0
        assert alias_verdict(same) and duplicate_verdict(same)
        shifted=geometry_agreement(ov,of,*_cube(1.,np.full(3,3.)))
        assert shifted['bbox_iou']==0.0 and not alias_verdict(shifted)
        nudged=geometry_agreement(ov,of,*_cube(1.,np.full(3,.05)))
        assert alias_verdict(nudged) and not duplicate_verdict(nudged),nudged
        small=geometry_agreement(*_cube(.4,np.full(3,.3)),*_cube(1.,np.zeros(3)))
        assert not alias_verdict(small),'a part must not read as an alias of its container'

        iv,if_=_cube(.4,np.full(3,.3))
        outer=scratch/'outer.json.gz';write_gz(outer,dict(positions=ov.ravel().tolist(),indices=of.ravel().tolist()))
        parents=[dict(entity_id='outer',name='outer',role='soft_organ',volume=1.0,
                      lo=ov.min(0),hi=ov.max(0),load=(lambda:surface_from_gz(outer)))]
        inside=containment(iv,parents)
        assert inside['part_of']==['outer'] and inside['tightest_container']=='outer',inside
        assert inside['relations'][0]['enclosed_vertex_fraction']==1.0
        away=containment(_cube(.4,np.full(3,9.))[0],parents)
        assert away['relations']==[] and away['part_of']==[] and away['tightest_container'] is None
        straddle=containment(_cube(.4,np.array([.8,.3,.3]))[0],parents)
        assert straddle['tightest_container'] is None,'a partly enclosed patch is not a part'
        assert straddle['overlaps']==['outer'] and not straddle['part_of'],straddle

        assert name_part_hint('Sternocostal head of pectoralis major muscle.r')['token']=='head of'
        assert name_part_hint('(Abdominal part of pectoralis major muscle).l')['token']=='part of'
        assert name_part_hint('Radius.l') is None
        assert refine_role('Anterior region of arm.l','integumentary','soft_organ')[0]=='surface_region'
        assert refine_role('Cavity of concha.l','integumentary','fluid_cavity')[0]=='fluid_cavity'
        assert refine_role('Popliteal artery.l','arterial','vascular')[0]=='vascular'

        materials=dict(materials=dict(muscle=dict(tissue='m',linear=dict(density=dict(value=1050.,tier='transferred',source='icru44_nist',note='n'))),
                                      connective_tissue=dict(tissue='c',linear=dict(density=dict(value=None,tier='absent',source=None,note='n')))))
        assert material_for('muscle',materials)['density_kg_m3']==1050.
        assert material_for('connective_tissue',materials)['density_kg_m3'] is None
        blank=material_for('surface_region',materials)
        assert blank['material_row'] is None and blank['density_kg_m3'] is None
        assert set(blank)==set(material_for('muscle',materials)),'every material record has one shape'

        v,f=_cube(1.,np.zeros(3))
        assert tetrahedralize(np.ascontiguousarray(v),np.ascontiguousarray(f),'pYq1.414')['succeeded']
        assert not tetrahedralize(np.ascontiguousarray(v),np.ascontiguousarray(f[:3]),'pYq1.414')['succeeded'],\
            'an open sheet must fail rather than be closed'

        report=deletion_report([dict(display_id='x',name='n',display_model='opensim-rajagopal',
                                     decision='delete',category='foreign_specimen',reason='r')],
                               [dict(display_id='y',name='n2',display_model='bodyparts3d',
                                     tier='id_alias',surviving_id='body-y')],[],{})
        assert report['counts']=={'remove':1,'collapse':1,'geometric_duplicates':0}
        coverage=provenance_coverage(dict(entities=[dict(id='e',provenance=dict(source_ids=['s'],
                                          files=[dict(path='p',sha256='h')],registration_id='z'),
                                          reference_geometry=dict(sha256='g'))]))
        assert coverage['coverage']['has_extraction_script']==0 and coverage['coverage']['has_source_file_hashes']==1
        import verify_muscle_tet_ready_surfaces as harness
        original=harness.tetgen.tetrahedralize
        try:
            harness.tetgen.tetrahedralize=lambda *args,**kwargs:time.sleep(60)
            timed=tetgen_with_timeout(np.ascontiguousarray(v),np.ascontiguousarray(f),'pYq1.414',1)
        finally:
            harness.tetgen.tetrahedralize=original
        assert timed['exception']=='Timeout' and not timed['succeeded'],timed
        assert tetgen_with_timeout(np.ascontiguousarray(v),np.ascontiguousarray(f),'pYq1.414',60)['succeeded'],\
            'the parent must survive a killed TetGen child'

        record=structure_provenance(
            dict(display_id='za-1',entity_id='body-za-1',name='n',system='muscular',role='muscle'),
            dict(label='Z-Anatomy',url='u',revision='r',specimen='s',units='b',frame='f',
                 attribution='a',status='acquired',license='CC-BY-SA-4.0',license_url='lu',
                 source_geometry_path='p.npz',source_geometry_sha256='h',sha256_verified=True,bytes=10),
            dict(path='g.json.gz',sha256='gh'),
            [dict(op='axis_rotation',residual=None,residual_reason='exact')],['A'],
            dict(script='scripts/x.py',script_sha256='sh',commit='c',
                 commit_covers_working_tree=True,uncommitted_changes=False),7,9)
        required=('dataset.id','dataset.label','dataset.license','source_file.path',
                  'source_file.sha256','build.script','build.commit','geometry.path',
                  'geometry.sha256','transforms','tier')
        for field in required:
            head,_,tail=field.partition('.')
            value=record[head][tail] if tail else record[head]
            assert value not in (None,'',[]),field
        assert record['tier']=='transferred' and record['frame_relation']=='registered_into_canonical'
        build_record=git_build_record(Path(__file__).resolve())
        assert build_record['script'].startswith('scripts/') and len(build_record['script_sha256'])==64
        assert not any(scratch.glob('tetgen-tmpfile_skipped.*'))
    print('PASS agreement, alias/duplicate thresholds, containment, part-of parse, role refinement, '
          'material selection, open-sheet refusal, deletion and coverage checks')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--self-test',action='store_true')
    p.add_argument('--output',type=Path)
    p.add_argument('--workers',type=int,default=6)
    p.add_argument('--limit',type=int,default=0)
    p.add_argument('--flags',default='pYq1.414')
    p.add_argument('--skip-tetgen',action='store_true')
    p.add_argument('--envelope-check',type=Path)
    p.add_argument('--tetgen-timeout',type=int,default=180)
    a=p.parse_args()
    if a.self_test:self_test()
    if a.output:
        globals()['TETGEN_TIMEOUT_S']=a.tetgen_timeout
        build(a.output,a.workers,a.limit,a.flags,a.skip_tetgen,a.tetgen_timeout)
    if a.envelope_check:envelope_check(a.envelope_check)
    if not a.self_test and not a.output and not a.envelope_check:
        p.error('Select --self-test, --output or --envelope-check')
