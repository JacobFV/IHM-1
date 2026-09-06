"""Prepare separately bound 98-muscle, seam-omitted static reference inputs."""
from pathlib import Path
import argparse,copy,json,hashlib,shutil,tempfile
import numpy as np
from frozen_static_stream import validate
ROOT=Path(__file__).resolve().parents[1]


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def retained_indices(face_indices,components,exterior):
    return np.where(np.asarray(components)[np.asarray(face_indices)]==exterior)[0]


def prepare():
    variant=ROOT/'data/derived/lumbar-muscle-native-lb45uirs/variant/registration.json';v=json.loads(variant.read_text())
    acceptance=ROOT/'data/derived/lumbar-muscle-native-lb45uirs/report.json';accepted=json.loads(acceptance.read_text())
    if not accepted['passed'] or accepted['muscles']!=98 or v['muscle_count']!=98:raise ValueError('Accepted98-muscle variant required')
    for name,digest in v['sources'].items():
        if sha(ROOT/name)!=digest:raise ValueError('Variant source changed: '+name)
    for key in ('model','catalog','insert'):
        if sha(ROOT/v[key+'_path'])!=v[key+'_sha256']:raise ValueError('Variant input hash changed')
    frozen=ROOT/'data/derived/frozen-static-stream-oay9phwy/manifest.json';validate(frozen)
    original=ROOT/'data/derived/supine-surface-contact-exmzq9pq/manifest.json';m=json.loads(original.read_text())
    topology=ROOT/'data/research/engineered_skin_territories/materialization.json';t=json.loads(topology.read_text())
    source=ROOT/m['arrays_path'];native=ROOT/m['native_input_path']
    if sha(source)!=m['arrays_sha256'] or sha(native)!=m['native_input_sha256']:raise ValueError('Old contact artifact changed')
    a=dict(np.load(source));keep=retained_indices(a['face_indices'],t['surface_diagnostic']['face_component_ids'],t['surface_diagnostic']['selected_component_id'])
    removed=np.setdiff1d(np.arange(len(a['face_indices'])),keep)
    if removed.tolist()!=[15948] or int(a['face_indices'][15948])!=165251:raise ValueError('Unexpected unresolved contact mask')
    output=Path(tempfile.mkdtemp(prefix='lumbar-supine-reference-',dir=ROOT/'data/derived'));contact=output/'contact';contact.mkdir()
    shutil.copytree(frozen.parent,output/'frozen_engine',symlinks=True);validate(output/'frozen_engine/manifest.json')
    revised={k:value[keep] for k,value in a.items()};revised['original_quadrature_indices']=keep
    arrays=contact/'quadrature.npz';np.savez_compressed(arrays,**revised)
    lines=native.read_text().splitlines();header=lines[0].split()
    if int(header[-1])!=len(a['face_indices']) or len(lines)!=len(a['face_indices'])+1:raise ValueError('Native contact row count mismatch')
    header[-1]=str(len(keep));newnative=contact/'supine_surface_foundation.txt';newnative.write_text(' '.join(header)+'\n'+'\n'.join(lines[k+1] for k in keep)+'\n')
    revised_manifest=copy.deepcopy(m)
    revised_manifest.update(native_input_path=str(newnative.relative_to(ROOT)),native_input_sha256=sha(newnative),arrays_path=str(arrays.relative_to(ROOT)),arrays_sha256=sha(arrays),points=len(keep),projected_area_m2=float(revised['area_m2'].sum()))
    revised_manifest['source_files']={str(p.relative_to(ROOT)):sha(p) for p in (original,source,native,topology,Path(__file__))}
    revised_manifest['omitted_unresolved_contact']=dict(original_quadrature_indices=removed.tolist(),source_face_indices=a['face_indices'][removed].tolist(),projected_area_m2=float(a['area_m2'][removed].sum()),basis='Exclude tiny seam component without substituting119mm deeper opposite ray hit; missing cell and open exterior uncertainty retained')
    revised_manifest['posterior_selection']='Frozen minimum-X quadrature subset: explicit unresolved seam cell omitted; all remaining face/station/area rows byte-equivalent numerically, no rerasterization'
    for row in revised_manifest['per_body']:
        select=revised['body_indices']==m['bodies'].index(row['body']);row['points']=int(select.sum());row['projected_area_m2']=float(revised['area_m2'][select].sum())
    surface=contact/'manifest.json';surface.write_text(json.dumps(revised_manifest,indent=2)+'\n')
    seed_source=ROOT/'data/derived/constrained-supine-f58do2qb/last_optimizer_iterate.json';old=json.loads(seed_source.read_text())
    seed={name:c['value'] for name,c in old['native']['coordinates'].items() if c['independent']}
    seed_path=output/'old_q_seed.json';seed_path.write_text(json.dumps(dict(coordinates=seed,basis='Old92-muscle pose only; unaccepted under98-muscle/seam-omitted physics; zero cached responses transferred',source_path=str(seed_source.relative_to(ROOT)),source_sha256=sha(seed_source)),indent=2)+'\n')
    paths=[variant,acceptance,ROOT/v['model_path'],ROOT/v['catalog_path'],surface,arrays,newnative,seed_path,output/'frozen_engine/manifest.json',ROOT/'ihm/assembly/bed_compression.py',Path(__file__)]
    record=dict(schema='ihm.lumbar-supine-static-reference.v1',frozen_engine=str((output/'frozen_engine/manifest.json').relative_to(ROOT)),variant_registration=str(variant.relative_to(ROOT)),surface_manifest=str(surface.relative_to(ROOT)),old_q_seed=str(seed_path.relative_to(ROOT)),material='MM',mass_kg=77.6122029,default_activation='Unmodified native/model defaults; six added lumbar muscles .05',cached_responses_reused=0,accepted_equilibrium=False,omitted_unresolved_contact=revised_manifest['omitted_unresolved_contact'],files={str(p.relative_to(ROOT)):sha(p) for p in paths})
    path=output/'manifest.json';path.write_text(json.dumps(record,indent=2)+'\n');return path

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--prepare',action='store_true');a=p.parse_args()
    if not a.prepare:raise SystemExit('Explicit source-only --prepare required')
    print(prepare())
