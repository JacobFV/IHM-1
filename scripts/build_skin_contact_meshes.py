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
import argparse,gzip,hashlib,json,sys
import numpy as np
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

def build(out_dir,reference_path,minimum_faces):
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
    reference=json.loads((ROOT/reference_path).read_text())
    registration=CanonicalRegistration(mechanics,reference)
    transform=np.linalg.inv(registration.global_map)
    source=canonical@transform[:3,:3].T+transform[:3,3]
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
        world=np.linalg.inv(np.asarray(reference['bodies'][segment]['transform_ground']))
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
    args=parser.parse_args()
    report=build(args.out,args.reference,args.minimum_faces)
    print(json.dumps({k:v for k,v in report.items() if k not in ('records','controls')},indent=2))
