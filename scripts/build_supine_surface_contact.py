"""Build a retained-skin posterior envelope and bounded contact-law fixtures."""
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import sys
import tempfile
import time
import xml.etree.ElementTree as ET
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from ihm.assembly.articulated import CanonicalRegistration
from ihm.assembly.supine_contact import posterior_envelope,foundation


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def layer_thickness(entity,surface_area,legacy_basis=None):
    shell=entity.get('shell') or {}
    if 'thickness_m' in shell:
        thickness=float(shell['thickness_m'])
        basis=dict(method='explicit_shell_thickness',prior_source=shell.get('prior_source'),shell=shell,
                   physical_surface_support=entity.get('physical_surface_support'))
    else:
        if entity.get('physical_surface_support') is not None:
            raise ValueError('Physical-support layer requires explicit shell thickness')
        if legacy_basis!='legacy_raw_full_skin_area_volume_ratio':
            raise ValueError('Legacy layer requires explicit legacy_raw_full_skin_area_volume_ratio basis')
        if not np.isfinite(surface_area) or surface_area<=0:raise ValueError('Invalid legacy surface area')
        thickness=float(entity['volume_m3'])/surface_area
        basis=dict(method=legacy_basis,warning='Legacy volume inferred from full raw skin area; not valid for exterior-area-regenerated layers',surface_area_m2=surface_area)
    if not np.isfinite(thickness) or thickness<=0:raise ValueError('Finite positive shell thickness required')
    return thickness,basis


def fixtures():
    vertices=np.array([[0,0,0],[0,.02,0],[0,.02,.02],[0,0,.02],[.001,0,0],[.001,.02,0],[.001,.02,.02],[.001,0,.02]])
    faces=np.array([[0,1,2],[0,2,3],[4,5,6],[4,6,7]])
    envelope=posterior_envelope(vertices,faces,.005)
    assert len(envelope['points_source_m'])==16 and np.all(envelope['face_indices']<2)
    assert np.isclose(envelope['area_m2'].sum(),.0004)
    material=dict(total_layer_thickness_m=.0066,shear_modulus_pa=1034.4827586206898,lame_lambda_pa=9310.34482758621,
                  minimum_thickness_ratio=.5,dissipation_s_m=2.,dynamic_friction=.8,viscous_friction=.5,transition_velocity_m_s=.2)
    p=np.array([[-.001,.3,.2],[-.002,-.1,.1]]);v=np.array([[-.1,.2,.3],[.05,-.2,.1]])
    area=np.array([.001,.002]);owners=np.array([0,1]);origins=np.array([[0,.2,0],[0,-.2,0]])
    result=foundation(p,v,area,owners,origins,0.,material)
    assert np.linalg.norm(result['force_balance_residual_n'])<1e-12 and np.linalg.norm(result['moment_balance_residual_nm'])<1e-12
    assert result['dissipative_power_w']<0
    zero=foundation(p,np.zeros_like(p),area,owners,origins,0.,material)
    eps=1e-8;plus=p.copy();minus=p.copy();plus[:,0]+=eps;minus[:,0]-=eps
    derivative=(foundation(plus,np.zeros_like(p),area,owners,origins,0.,material)['elastic_energy_j']-foundation(minus,np.zeros_like(p),area,owners,origins,0.,material)['elastic_energy_j'])/(2*eps)
    assert np.isclose(-derivative,zero['point_forces_n'][:,0].sum(),rtol=1e-7)
    separated=p.copy();separated[:,0]=.001
    assert foundation(separated,v,area,owners,origins,0.,material)['elastic_energy_j']==0
    try:foundation(p-.01,np.zeros_like(p),area,owners,origins,0.,material)
    except ValueError:pass
    else:raise AssertionError('Excessive compression accepted')
    return dict(passed=True,native_run=False,checks=['exact retained-face envelope','inner surface not double counted',
         'projected area','equal/opposite force and moment','force is negative energy gradient','dissipative work','no tension','compression domain rejection'])


def build(legacy_thickness_basis=None):
    started=time.monotonic()
    mechanics_path=ROOT/'data/derived/canonical/mechanics.json';mechanics=json.loads(mechanics_path.read_text())
    skin=next(e for e in mechanics['entities'] if e['role']=='skin')
    mesh_path=ROOT/skin['reference_geometry']['path']
    if sha(mesh_path)!=skin['reference_geometry']['sha256']:raise ValueError('Skin identity mismatch')
    geometry=json.loads(gzip.decompress(mesh_path.read_bytes()))
    canonical=np.array(geometry['positions']).reshape(-1,3);faces=np.array(geometry['indices']).reshape(-1,3)
    reference_path=ROOT/'data/derived/supine-support-5ma720yd/initial_native.json'
    native=json.loads(reference_path.read_text());registration=CanonicalRegistration(mechanics,native)
    transform=np.linalg.inv(registration.global_map)
    source=canonical@transform[:3,:3].T+transform[:3,3]
    envelope=posterior_envelope(source,faces)
    points=envelope['points_source_m'];canonical_points=points@registration.basis.T+registration.global_map[:3,3]
    bodies=list(registration.groups)
    distances=[]
    for body in bodies:
        group=registration.groups[body]
        distances.append(np.linalg.norm(np.maximum(np.maximum(group['bounds_min_m']-canonical_points,canonical_points-group['bounds_max_m']),0),axis=1))
    ranking=np.argsort(np.array(distances),axis=0,kind='stable');owners=ranking[0]
    local=np.zeros_like(points)
    for index,body in enumerate(bodies):
        inverse=np.linalg.inv(np.array(native['bodies'][body]['transform_ground']))
        selected=owners==index;local[selected]=points[selected]@inverse[:3,:3].T+inverse[:3,3]
    triangles=canonical[faces];surface_area=float(np.linalg.norm(np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0]),axis=1).sum()/2)
    layers=[]
    for entity in mechanics['entities']:
        if entity['role']=='skin_layer':
            thickness,basis=layer_thickness(entity,surface_area,legacy_thickness_basis)
            layers.append(dict(id=entity['id'],thickness_m=thickness,thickness_basis=basis,
                               young_modulus=entity['material']['young_modulus'],poisson_ratio=entity['material']['poisson_ratio']))
    parameters={(layer['young_modulus']['value'],layer['poisson_ratio']['value']) for layer in layers}
    if len(parameters)!=1:raise ValueError('Series law requires identical held layer material; no silent effective averaging')
    young,nu=parameters.pop();mu=young/(2*(1+nu));lam=young*nu/((1+nu)*(1-2*nu))
    contact_path=ROOT/'data/derived/supine-support-5ma720yd/plant/native/inputs/subject_walk_scaled_ContactForceSet.xml'
    contact=next(ET.fromstring(contact_path.read_bytes()).iter('SmoothSphereHalfSpaceForce'))
    static=float(contact.findtext('static_friction'));dynamic=float(contact.findtext('dynamic_friction'))
    if static!=dynamic:raise ValueError('No static/dynamic blending silently introduced')
    material=dict(total_layer_thickness_m=sum(l['thickness_m'] for l in layers),shear_modulus_pa=mu,lame_lambda_pa=lam,
        minimum_thickness_ratio=.5,dissipation_s_m=float(contact.findtext('dissipation')),dynamic_friction=dynamic,
        viscous_friction=float(contact.findtext('viscous_friction')),transition_velocity_m_s=float(contact.findtext('transition_velocity')))
    output=Path(tempfile.mkdtemp(prefix='supine-surface-contact-',dir=ROOT/'data/derived'))
    arrays=output/'quadrature.npz';np.savez_compressed(arrays,stations_local_m=local,reference_points_source_m=points,
         body_indices=owners,face_indices=envelope['face_indices'],area_m2=envelope['area_m2'])
    # Plane is the full retained posterior surface minimum, not an inertia radius.
    plane=float(source[:,0].min())
    native_input=output/'supine_surface_foundation.txt'
    order=('total_layer_thickness_m','shear_modulus_pa','lame_lambda_pa','minimum_thickness_ratio','dissipation_s_m','dynamic_friction','viscous_friction','transition_velocity_m_s')
    header=['IHM_SURFACE_FOUNDATION_V1',str(plane),*[str(material[k]) for k in order],str(len(points))]
    with native_input.open('w') as handle:
        handle.write(' '.join(header)+'\n')
        for index,station,area in zip(owners,local,envelope['area_m2']):
            handle.write(' '.join([bodies[index],*map(str,station),str(area)])+'\n')
    per_body=[]
    for index,body in enumerate(bodies):
        selected=owners==index
        if selected.any():per_body.append(dict(body=body,points=int(selected.sum()),projected_area_m2=float(envelope['area_m2'][selected].sum()),
             minimum_gap_m=float(points[selected,0].min()-plane)))
    manifest=dict(schema='ihm.supine-skin-foundation.v1',accepted_support=False,native_integration=False,
        native_input_path=str(native_input.relative_to(ROOT)),native_input_sha256=sha(native_input),
        bodies=bodies,arrays_path=str(arrays.relative_to(ROOT)),arrays_sha256=sha(arrays),plane_source_x_m=plane,
        source_files={str(p.relative_to(ROOT)):sha(p) for p in (mechanics_path,mesh_path,reference_path,contact_path,Path(__file__).resolve(),ROOT/'ihm/assembly/supine_contact.py')},
        registration=registration.global_fit,attachment_basis='Nearest named bone envelope per ray/triangle-linked point; uncertain anatomical ownership, no added inertia',
        posterior_selection='Minimum source-X triangle intersection per 5 mm YZ raster cell; full source triangle index retained; inner/back surfaces do not duplicate ray area',
        points=len(points),source_faces=len(faces),source_skin_surface_area_m2=surface_area,projected_area_m2=float(envelope['area_m2'].sum()),
        sampled_minimum_gap_m=float(points[:,0].min()-plane),spacing_m=.005,per_body=per_body,
        material=material,layers=layers,material_scope='Retained generic skin priors; no patient or mattress calibration; confined neo-Hookean columns in series',
        compression_domain='Thickness ratio >=0.5 is an explicit engineering validity ceiling, not empirically validated tissue strain',
        limits=['Reference posterior envelope only; large rotations exposing unsampled surfaces require rebuilt quadrature',
                '5 mm reference projected-area quadrature; boundary-cell area and between-ray extrema error unresolved',
                'Global COM-to-bone-envelope registration is approximate; contact geometry is not registration validation',
                'Stationary rigid plane; local skin-column compression, not a calibrated deforming mattress',
                'Alternative contact owner must replace old spheres, never add them'],wall_s=time.monotonic()-started)
    (output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(dict(output=str(output),points=len(points),plane_source_x_m=plane,projected_area_m2=manifest['projected_area_m2'],per_body=per_body,wall_s=manifest['wall_s']),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--fixture-check',action='store_true');parser.add_argument('--legacy-thickness-basis',choices=['legacy_raw_full_skin_area_volume_ratio']);args=parser.parse_args()
    if args.fixture_check:print(json.dumps(fixtures(),indent=2))
    else:build(args.legacy_thickness_basis)
