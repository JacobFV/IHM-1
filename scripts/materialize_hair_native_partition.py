"""Bounded source-only inferred hair partition; never writes a native model."""
from pathlib import Path
import io,json,sys,zipfile
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.assembly.source_skin import file_sha256,read_small_member,SourceSkinPatch
from ihm.assembly.hair_source_factory import prepare_hair_factory,guide_mass_properties
from ihm.assembly.hair_source_ownership import classify_root_eligibility,select_guide_geometry,select_source_faces,infer_source_owners,partition_represented_guides
TOPOLOGY='data/research/engineered_skin_territories/materialization.json'
TOPOLOGY_SHA='de669f96f6db5b150dd467f149464677e19f22ec2e7be5a77d895d873b3fef1a'


def materialize(root=ROOT):
    root=Path(root);dependencies={}
    def retained(path,expected=None):
        path=Path(path);path=path if path.is_absolute() else root/path
        with path.open('rb') as handle:raw=handle.read(32*1024*1024+1)
        if len(raw)>32*1024*1024:raise ValueError('Retained input byte budget exceeded')
        digest=file_sha256(raw);key=str(path.relative_to(root))
        if expected is not None and digest!=expected:raise ValueError('Source dependency hash mismatch: '+key)
        if key in dependencies and dependencies[key]!=digest:raise ValueError('Input changed during materialization')
        dependencies[key]=digest;return raw
    base=json.loads(retained('data/research/hair_source_factory.json'))
    for path,digest in base['source_dependencies'].items():retained(path,digest)
    native=base['actual_frozen_native_inertia'];registration=json.loads(retained(native['registration_path'],native['registration_sha256']))
    contact_path=str(Path(native['assembled_path']).parent/'inputs/surface_contact_manifest.json');contact=json.loads(retained(contact_path))
    if contact['bodies']!=list(registration['groups']) or not np.array_equal(contact['registration']['source_to_canonical_ground'],registration['global_rigid_fit']['source_to_canonical_ground']):raise ValueError('Neutral contact/native registration mismatch')
    topology=json.loads(retained(TOPOLOGY,TOPOLOGY_SHA));diagnostic=topology['surface_diagnostic'];labels=np.asarray(diagnostic['face_component_ids']);selected=diagnostic['selected_component_id'];permitted=topology['contact_eligible_triangle_ids']
    groups={};patches={};geometries={};populations={}
    for name,group in base['groups'].items():
        patch=group['source_patch'];source=next(r for r in topology['source_receipts'] if r['entity_id']==patch['source_id'])
        if source['sha256']!=patch['source_sha256'] or len(labels)!=patch['source_face_count'] or contact['source_files'].get(source['path'])!=source['sha256']:raise ValueError('Exact topology/contact source identity mismatch')
        path=f'data/derived/hair/elastic_v3/body-dynamic-{name}-hair.json.gz';raw=retained(path,base['source_dependencies'][path]);geometry={k:read_small_member(raw,k) for k in ('attachment','strands')}
        raw=retained(group['population_path'],group['population_sha256']);keys=('ids','face_index','barycentric','roots_m','normals','radius_m','length_m','source_entity_id','source_geometry_sha256')
        with zipfile.ZipFile(io.BytesIO(raw)) as archive:
            if sum(i.file_size for i in archive.infolist() if i.filename in {k+'.npy' for k in keys})>128*1024*1024:raise ValueError('Population decode budget exceeded')
        with np.load(io.BytesIO(raw),allow_pickle=False) as archive:population={k:archive[k] for k in keys if k in archive}
        eligibility=classify_root_eligibility(group['mapping'],labels,permitted,selected_component_id=selected)
        for rejected in eligibility['ineligible_guides']:
            rejected['represented_mass_kg']=guide_mass_properties(select_guide_geometry(geometry,[rejected['guide_index']])['strands'])['mass_kg']
        indices=eligibility['eligible_guide_indices'];geometries[name]=select_guide_geometry(geometry,indices);populations[name]=population
        patches[name]=select_source_faces(patch,np.asarray(group['mapping']['source_face_indices'])[indices])
        faces=population['face_index']
        if faces.dtype.kind not in 'iu' or np.any(faces<0) or np.any(faces>=len(labels)):raise ValueError('Invalid original population face index')
        components=labels[faces];mass=geometry['strands']['density_kg_m3']*np.pi*population['radius_m']**2*population['length_m'];eligible=(components==selected)&np.isin(faces,permitted)
        groups[name]={'guide_eligibility':eligibility,'population_eligibility':{'original_count':len(faces),'eligible_count':int(eligible.sum()),'ineligible_count':int((~eligible).sum()),'eligible_prior_mass_kg':float(mass[eligible].sum()),'ineligible_prior_mass_kg':float(mass[~eligible].sum()),'components':[{'component_id':int(c),'count':int((components==c).sum()),'prior_mass_kg':float(mass[components==c].sum())} for c in np.unique(components)]},'source_patch':patches[name]}
    ownership=infer_source_owners(patches,registration,registration_sha256=native['registration_sha256']);prepared={}
    for name in groups:
        bound=SourceSkinPatch.bind(patches[name],ownership['bindings'][name],expected_registration_sha256=native['registration_sha256'],owner_names=list(registration['groups']))
        prepared[name]=prepare_hair_factory(geometries[name],populations[name],bound);groups[name]['prepared']=prepared[name]
    partition=partition_represented_guides(prepared,native['bodies'])
    for path in (Path(__file__),root/'ihm/assembly/hair_source_ownership.py'):retained(path)
    return {'schema':'ihm.hair-native-partition-preparation.v1','native_enabled':False,'installed':False,'source_dependencies':dependencies,'groups':groups,'ownership':ownership,'partition':partition,
        'topology_evidence':{'artifact':TOPOLOGY,'sha256':TOPOLOGY_SHA,'selected_component_id':selected,'selection_evidence':diagnostic['selection_evidence'],'physical_surface_exclusivity_validated':diagnostic['physical_surface_exclusivity_validated'],'signed_integral_is_closed_volume':diagnostic['signed_integral_is_closed_volume']},
        'scope':'Current reduced 22-body reference topology only. Scalp maps to torso because no separate native head exists; this is a structural approximation, not anatomical scalp ownership. Revised cervical topology requires regeneration.',
        'mass_basis':'Explicit inferred partition of actual frozen patient-scaled native segment inertia for 94 independent exterior guides. Canonical hair proxy metadata is not a separately integrated mass and is never debited. All unmaterialized hair remains lumped in residual bodies.',
        'limits':['Two retained body roots on inner shell remain unchanged and physically ineligible; population and original guide assets are preserved','Exterior selection is an engineering surface proxy with open boundaries, not closed-volume or anatomical exclusivity proof','Sparse follicle faces are not complete moving body/garment contact coverage','Nearest bone envelopes and connected-face coherence are explicit inferred priors; ambiguity and overrides remain retained','Native installation, matched moving-state initialization, single coupled interval and authoritative viewer synchronization remain absent']}


if __name__=='__main__':
    report=materialize();path=ROOT/'data/research/hair_native_partition.json';path.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'path':str(path),'ownership':report['ownership']['summary'],'partition':{k:report['partition'][k] for k in ('represented_guide_count','represented_mass_kg','modified_owner_count')},'native_enabled':False}))
