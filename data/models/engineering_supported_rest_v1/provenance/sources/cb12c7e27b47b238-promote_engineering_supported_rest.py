"""Freeze exact supported-rest engineering model/pose/evidence without research runtime paths."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def sha(raw):return hashlib.sha256(raw).hexdigest()
def main():
    out=ROOT/'data/models/engineering_supported_rest_v1'
    if out.exists():raise ValueError('Fresh bundle directory required')
    variant=ROOT/'data/derived/supine-selected-defaults-ih90ww49';reg=json.loads((variant/'registration.json').read_bytes());pose=json.loads((variant/'initial_pose.json').read_bytes());surface=json.loads((ROOT/pose['support_geometry']['manifest_path']).read_bytes());initial=json.loads((ROOT/'data/derived/supine-equilibrium-initial-xfhqpo4e/initialized.json').read_bytes())
    if reg['model_sha256']!='16af7ffdc6d346b553d5ba9ba24a85f69248a493f844197ca15af2de4409769d' or pose['target_mass_kg']!=77.6122029:raise ValueError('Unexpected owned model or mass')
    out.mkdir(parents=True);files={};origins={};by_digest={};canonical={}
    def write(name,raw,origin=None):
        path=out/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(raw);relative=str(path.relative_to(ROOT));digest=sha(raw);files[relative]=digest;by_digest[digest]=relative
        if origin:origins[str(origin)]={'bundle_path':relative,'sha256':digest}
        return relative
    def retain(origin,digest=None):
        raw=(ROOT/origin).read_bytes();actual=sha(raw)
        if digest is not None and actual!=digest:raise ValueError('Changed source '+origin)
        if actual in by_digest:
            origins[origin]={'bundle_path':by_digest[actual],'sha256':actual};return by_digest[actual]
        return write('provenance/sources/'+actual[:16]+'-'+Path(origin).name,raw,origin)
    def source_binding(origin,digest):
        if origin.startswith('data/derived/canonical/'):
            if sha((ROOT/origin).read_bytes())!=digest:raise ValueError('Patient canonical source mismatch')
            canonical[origin]=digest;return origin
        return retain(origin,digest)
    promoted=dict(reg);promoted['model_path']=write('model.osim',(ROOT/reg['model_path']).read_bytes(),reg['model_path'])
    if files[promoted['model_path']]!=reg['model_sha256']:raise ValueError('Model changed')
    promoted['sources']={source_binding(p,d):d for p,d in reg['sources'].items()}
    catalog=json.loads((ROOT/reg['catalog_path']).read_bytes());retain(reg['catalog_path'],reg['catalog_sha256'])
    for row in catalog:
        original=row['source_path'];row['source_path']=source_binding(original,row['source_sha256'])
        if promoted['sources'].get(row['source_path'])!=row['source_sha256']:raise ValueError('Catalog source ownership mismatch')
    promoted['catalog_path']=write('catalog.json',(json.dumps(catalog,indent=2)+'\n').encode());promoted['catalog_sha256']=files[promoted['catalog_path']]
    promoted['insert_path']=retain(reg['insert_path'],reg['insert_sha256'])
    for key in ('base_model_path','rematerialization_receipt'):promoted[key]=retain(reg[key])
    promoted['activation_initialization']={**reg['activation_initialization'],'pose_path':str((out/'initial_pose.json').relative_to(ROOT)),'scope':'Explicit26-muscle inverse-statics resting tone; not learned cortex. All98 excitations retained separately.'}
    promoted['target_mass_kg']=pose['target_mass_kg'];promoted['default_enabled']=False;promoted['native_acceptance_complete']=True;promoted['canonical_source_bindings']=canonical
    write('provenance/original_registration.json',(variant/'registration.json').read_bytes(),str((variant/'registration.json').relative_to(ROOT)))
    support=dict(surface);support['source_files']={source_binding(p,d):d for p,d in surface['source_files'].items()}
    for key,name in [('native_input','surface/supine_surface_foundation.txt'),('arrays','surface/quadrature.npz')]:
        raw=(ROOT/surface[key+'_path']).read_bytes()
        if sha(raw)!=surface[key+'_sha256']:raise ValueError('Surface artifact changed')
        support[key+'_path']=write(name,raw,surface[key+'_path'])
    def remap_known(value):
        if isinstance(value,dict):return {k:remap_known(v) for k,v in value.items()}
        if isinstance(value,list):return [remap_known(v) for v in value]
        if isinstance(value,str) and value in origins:return origins[value]['bundle_path']
        return value
    support=remap_known(support);supportpath=write('surface/manifest.json',(json.dumps(support,indent=2)+'\n').encode())
    write('provenance/original_surface_manifest.json',(ROOT/pose['support_geometry']['manifest_path']).read_bytes(),pose['support_geometry']['manifest_path'])
    beddir=ROOT/'data/research/bed_material';bed=json.loads((beddir/'manifest.json').read_bytes())
    for row in bed['sources']:
        raw=(beddir/row['path']).read_bytes()
        if sha(raw)!=row['sha256']:raise ValueError('Bed evidence changed')
        write('bed/'+row['path'],raw,str((beddir/row['path']).relative_to(ROOT)))
    table=bed['digitization']['path'];raw=(beddir/table).read_bytes()
    if sha(raw)!=bed['derived_files'][Path(table).name]['sha256']:raise ValueError('Bed curve changed')
    write('bed/'+table,raw,str((beddir/table).relative_to(ROOT)));bedpath=write('bed/manifest.json',(beddir/'manifest.json').read_bytes(),str((beddir/'manifest.json').relative_to(ROOT)))
    selected=sorted(pose['activations']);alltones={n:row['activation'] for n,row in sorted(initial['muscles'].items())}
    if len(selected)!=26 or len(alltones)!=98 or any(alltones[n]!=pose['activations'][n] for n in selected):raise ValueError('Exact selected tone receipt differs')
    promotedpose={**pose,'coordinates':dict(sorted(pose['coordinates'].items())),'registration_path':str((out/'registration.json').relative_to(ROOT)),'support_geometry':{'manifest_path':supportpath,'manifest_sha256':files[supportpath],'bed_material':'MM','bed_material_manifest_path':bedpath}}
    promotedpose['activation_selection_path']=retain(pose['activation_selection_path']);promotedpose['dynamic_seed']=remap_known(pose.get('dynamic_seed',{}))
    posepath=write('initial_pose.json',(json.dumps(promotedpose,indent=2)+'\n').encode());write('equilibrium_excitations.json',(json.dumps(alltones,indent=2)+'\n').encode());promoted['activation_initialization']['pose_sha256']=files[posepath]
    for path in (promoted['model_path'],promoted['catalog_path'],posepath):promoted['sources'][path]=files[path]
    registrationpath=write('registration.json',(json.dumps(promoted,indent=2)+'\n').encode())
    receipts={}
    for label,path in [('static_initializer','data/derived/supine-equilibrium-initial-xfhqpo4e/report.json'),('static_state','data/derived/supine-equilibrium-initial-xfhqpo4e/static_actual_initialized_coordinates.json'),('initialized_state','data/derived/supine-equilibrium-initial-xfhqpo4e/initialized.json'),('short_forward','data/derived/supine-equilibrium-initial-xfhqpo4e/advanced_1us.json'),('initializer_equivalence','data/derived/supine-lumbar-default-receipt-5zunf4qh/report.json'),('coupled_world_1s','data/derived/supine-tonic-world-6u1wzfk9/report.json'),('solver','data/derived/supine-supported-tones-_wh_r3gi/solver.py')]:receipts[label]=retain(path)
    retain(str(Path(__file__).relative_to(ROOT)))
    from ihm.native.mechanical_stream import SOURCE_FILES
    explicit_inputs={}
    for name in SOURCE_FILES:
        if name=='subject_walk_scaled.osim':continue
        path='data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/example3DWalking/'+name;explicit_inputs[path]=sha((ROOT/path).read_bytes())
    manifest=dict(schema='ihm.engineering-supported-rest-bundle.v1',target_mass_kg=pose['target_mass_kg'],native_model_sha256=reg['model_sha256'],selected_tonic_muscles=selected,surface_contact_manifest_path=supportpath,bed_material='MM',bed_material_manifest_path=bedpath,bed_material_manifest_sha256=files[bedpath],registration_path=registrationpath,initial_pose_path=posepath,files=files,origins=origins,acceptance=receipts,canonical_source_bindings=canonical,explicit_native_source_inputs=explicit_inputs,initial_metabolic_reference=initial['metabolic_reference'],mass_ownership='Single77.6122029kg patient-matched native body; no instance mass override validated',default_enabled=False,brain_trained=False,walking_demonstrated=False,scope='Numerical unloaded static equilibrium plus retained1s coupled-cloth measurement; engineering baseline is explicit, not cortical output',archival_paths='Historical origin strings inside immutable receipts are archival identifiers; origins maps retained files. Runtime registration/surface sources resolve only within bundle or explicit canonical bindings.')
    total=sum(p.stat().st_size for p in out.rglob('*') if p.is_file());manifest['closure_bytes_excluding_manifest']=total
    if total>=64*1024*1024:raise ValueError('Bundle exceeds64MiB')
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    for path,digest in files.items():
        if sha((ROOT/path).read_bytes())!=digest:raise ValueError('Bundle integrity failure')
    print(json.dumps(dict(output=str(out),files=len(files),closure_bytes=sum(p.stat().st_size for p in out.rglob('*') if p.is_file()),canonical_source_bindings=canonical,explicit_native_source_inputs=explicit_inputs),indent=2))
if __name__=='__main__':main()
