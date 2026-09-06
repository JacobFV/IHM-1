"""Prepare retained hair mapping/inertia evidence without loading a whole mesh.

Small source-only run. No native process, nearest-face query, owner inference,
solver advance, population reduction, or live factory enablement.
"""
from pathlib import Path
import io,json,sys,xml.etree.ElementTree as ET,zipfile
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.assembly.source_skin import file_sha256,read_small_member,extract_source_patch
from ihm.assembly.hair_source_factory import map_retained_guides,validate_retained_source,guide_mass_properties,population_mass_properties


def audit(root=ROOT):
    root=Path(root);dependencies={}
    def retained(path):
        with path.open('rb') as handle:raw=handle.read(32*1024*1024+1)
        if len(raw)>32*1024*1024:raise ValueError('Retained input exceeds explicit source audit byte budget')
        key=str(path.relative_to(root));digest=file_sha256(raw)
        if key in dependencies and dependencies[key]!=digest:raise ValueError('Retained input changed during audit')
        dependencies[key]=digest;return raw
    folder=root/'data/derived/hair/elastic_v3';manifest_path=folder/'manifest_fragment.json';manifest=json.loads(retained(manifest_path))
    source_path=root/manifest['source_geometry_path'];source_raw=retained(source_path);source_sha=manifest['source_geometry_sha256'];groups={}
    for name in ('scalp','body'):
        record=next(s for s in manifest['structures'] if s['id']==f'body-dynamic-{name}-hair');geometry_path=root/record['geometry_path']
        geometry_raw=retained(geometry_path)
        if file_sha256(geometry_raw)!=record['geometry_sha256']:raise ValueError('Retained elastic geometry hash mismatch')
        geometry={key:read_small_member(geometry_raw,key) for key in ('attachment','strands')}
        population_path=folder/'scalp_population.npz' if name=='scalp' else root/'data/derived/canonical/hair_samples.npz'
        population_raw=retained(population_path);population_sha=file_sha256(population_raw)
        if name=='body' and population_sha!=manifest['body_population_sha256']:raise ValueError('Original body population hash mismatch')
        keys=('ids','face_index','barycentric','roots_m','normals','radius_m','length_m','region','source_entity_id','source_geometry_sha256')
        with zipfile.ZipFile(io.BytesIO(population_raw)) as archive:
            requested=[i for i in archive.infolist() if i.filename in {k+'.npy' for k in keys}]
            if sum(i.file_size for i in requested)>128*1024*1024:raise ValueError('Population members exceed bounded decode budget')
        with np.load(io.BytesIO(population_raw),allow_pickle=False) as archive:
            population={key:archive[key] for key in keys if key in archive}
        mapping=map_retained_guides(geometry['attachment'],population,source_id=geometry['attachment']['skin_entity_id'],source_sha256=source_sha)
        patch=extract_source_patch(source_raw,source_sha,mapping['source_face_indices'],source_id=geometry['attachment']['skin_entity_id'])
        validated=validate_retained_source(geometry,population,patch)
        groups[name]={'mapping':validated,'source_patch':patch,'represented_guides':guide_mass_properties(geometry['strands']),
                      'full_population_prior':population_mass_properties(population,geometry['strands']['density_kg_m3']),
                      'population_path':str(population_path.relative_to(root)),'population_sha256':population_sha,
                      'population_hash_basis':'Original body population digest retained in elastic manifest' if name=='body' else 'Current retained scalp bytes pinned by this audit; older elastic manifest omitted a population digest',
                      'owner_binding':None,'native_enabled':False}
        if 'region' in population:
            names,counts=np.unique(population['region'],return_counts=True);groups[name]['population_regions']=dict(zip(names.tolist(),map(int,counts)))
    mechanics_path=root/'data/derived/canonical/mechanics.json';mechanics=json.loads(retained(mechanics_path))
    proxies=[{k:e.get(k) for k in ('id','name','mass_kg','mass_role','centroid_m','inertia_diagonal_kg_m2','reference_geometry')} for e in mechanics['entities'] if e['id'] in ('body-bp3d-FJ2813','body-bp3d-FJ2815')]
    registration_path=root/'data/derived/mechanics/whole_body_arm26_v2/registration.json';registration=json.loads(retained(registration_path));model_path=root/registration['model_path'];model_raw=retained(model_path)
    if file_sha256(model_raw)!=registration['model_sha256']:raise ValueError('Native source model provenance mismatch')
    bodies=[]
    for body in ET.fromstring(model_raw).findall('.//BodySet/objects/Body'):
        mass=float(body.findtext('mass'));center=list(map(float,body.findtext('mass_center').split()));inertia=list(map(float,body.findtext('inertia').split()))
        bodies.append({'owner':body.attrib['name'],'mass_kg':mass,'mass_center_local_m':center,'inertia_com_local_xx_yy_zz_xy_xz_yz_kg_m2':inertia})
    frozen=root/'data/derived/audits/cutaneous-factory-ylro2d66/body/mechanics';assembled_path=frozen/'native/assembled_model.osim';identity_path=frozen/'identity.json';reference_path=frozen/'registration.json'
    frozen_identity=json.loads(retained(identity_path));reference_raw=retained(reference_path);reference=json.loads(reference_raw);assembled_raw=retained(assembled_path);assembled=[]
    for body in ET.fromstring(assembled_raw).findall('.//BodySet/objects/Body'):
        owner=body.attrib['name'];mass=float(body.findtext('mass'));center=list(map(float,body.findtext('mass_center').split()));v=list(map(float,body.findtext('inertia').split()))
        tensor=[[v[0],v[3],v[4]],[v[3],v[1],v[5]],[v[4],v[5],v[2]]]
        group=reference['groups'][owner]
        assembled.append({'owner':owner,'frame':'native-body-local:'+owner,'mass_kg':mass,'centroid_m':center,'inertia_com_kg_m2':tensor,
                          'canonical_reference_to_body_local':group['canonical_reference_embedding_in_source_body'],'native_reference_transform':group['native_reference_transform']})
    if abs(sum(b['mass_kg'] for b in assembled)-frozen_identity['effective_native_body_mass_kg'])>1e-9:raise ValueError('Frozen assembled inertia disagrees with actual factory mass identity')
    for path in (root/'data/measurements/hair/strand_reference.json',Path(__file__),root/'ihm/assembly/source_skin.py',root/'ihm/assembly/hair_source_factory.py',root/'ihm/assembly/hair_dynamics.py'):retained(path)
    return {'schema':'ihm.hair-source-factory-preparation.v1','native_enabled':False,'source_dependencies':dependencies,'groups':groups,
            'canonical_proxy_metadata':proxies,
            'actual_frozen_native_inertia':{'assembled_path':str(assembled_path.relative_to(root)),'assembled_sha256':file_sha256(assembled_raw),'identity_path':str(identity_path.relative_to(root)),
                'registration_path':str(reference_path.relative_to(root)),'registration_sha256':file_sha256(reference_raw),'effective_mass_kg':frozen_identity['effective_native_body_mass_kg'],'bodies':assembled,
                'status':'Actual retained cutaneous-factory assembled body inertias and native reference transforms. No running session is started or modified; exact node-owner binding and inferred hair partition remain absent.'},
            'native_source_inertia':{'model_path':str(model_path.relative_to(root)),'model_sha256':registration['model_sha256'],'total_source_mass_kg':sum(b['mass_kg'] for b in bodies),'bodies':bodies,
                'status':'Retained source segment inertias, not a running patient-scaled snapshot. Articulated uses native segment inertia exclusively; canonical hair proxies are not separately integrated.'},
            'debit_recipe':[
                'Freeze the actual active patient-scaled native segment mass, COM, tensor and common reference transforms; source XML values alone are insufficient.',
                'Supply one source-hash/registration-hash-bound fixed owner for each retained source node. Join sample IDs to population face indices exactly; never nearest-face reattach.',
                'Compute full-population physical mass/moments from retained density/radius/length priors as plausibility evidence. They are not calibrated subject partitions or guide weights.',
                'Declare an inferred hair partition of each actual native owner. Debit only the explicitly represented independent guide nodal mass, first moment and second moment in that same body reference frame. Unrepresented population remains in the residual owner budget.',
                'Require positive residual mass and centered second-moment tensor; preserve total mass, first moment and full inertia when residual native owner and explicit guides are recombined. Do not subtract canonical proxy masses or weight guides by rendered population.',
                'Apply no change until a source-bound native mass/COM/inertia installation path and one coupled common-clock interval are verified. Compose hair with existing mechanical feedback without independently advancing the body twice.'
            ],
            'blockers':['No frozen exact source-node owner binding supplied; correspondence is ready but SourceSkinPatch.bind has not run on anatomical data',
                        'Actual frozen patient-scaled native inertias are retained, but no inferred segment partition/debit is selected or installed; proxy metadata is not a native debit budget',
                        'Original scalp population lacks predecessor digest; this receipt pins retained bytes and verifies their exact source correspondence',
                        'Only requested follicle faces are extracted; this sparse patch is not complete nearby body/garment contact coverage',
                        'No native/live owner, coupled interval composition or authoritative viewer synchronization is enabled'],
            'scope':'Source-only exact retained guide preparation. Population counts, fibers, source geometry and native inertia are unchanged.'}


if __name__=='__main__':
    report=audit();destination=ROOT/'data/research/hair_source_factory.json';destination.parent.mkdir(parents=True,exist_ok=True);destination.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'artifact':str(destination.relative_to(ROOT)),'native_source_mass_kg':report['native_source_inertia']['total_source_mass_kg'],
                      'actual_frozen_native_mass_kg':report['actual_frozen_native_inertia']['effective_mass_kg'],
                      'groups':{k:{'guides':v['represented_guides']['guide_count'],'guide_mass_kg':v['represented_guides']['mass_kg'],'population':v['full_population_prior']['strand_count'],'population_prior_mass_kg':v['full_population_prior']['mass_kg'],'retained_faces':len(v['source_patch']['source_face_indices']),'retained_nodes':len(v['source_patch']['source_node_indices'])} for k,v in report['groups'].items()},'native_enabled':False},indent=2))
