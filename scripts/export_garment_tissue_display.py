"""Export verified physical trajectories; boundary extraction does not alter volumes."""
from pathlib import Path
import argparse,gzip,hashlib,json
import numpy as np

BASE=Path(__file__).resolve().parents[1]
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def _resolve(path):
    p=Path(path);return p if p.is_absolute() else BASE/p
def _hash_array(value,dtype):return hashlib.sha256(np.asarray(value,dtype=dtype).tobytes()).hexdigest()

def export_display(experiment,output):
    source=_resolve(experiment);out=_resolve(output)
    if out.exists():raise ValueError('Retain existing display artifacts; choose a fresh output')
    report=json.loads((source/'report.json').read_text());config=json.loads((source/'configuration.json').read_text())
    if report['status']!='completed_exploratory_fixture' or report['minimum_jacobian']<=0:raise ValueError('Only completed positive-volume cases can be exported')
    for filename,digest in report['artifacts_sha256'].items():
        if sha(source/filename)!=digest:raise ValueError('Changed physical experiment artifact')
    for path,digest in config['source_sha256'].items():
        if sha(source/'source-snapshot'/Path(path).name)!=digest:raise ValueError('Changed solver snapshot')
    manifest_path=_resolve(config['tissue_manifest_path'])
    if sha(manifest_path)!=config['tissue_manifest_sha256']:raise ValueError('Changed tissue manifest')
    domain=json.loads(manifest_path.read_text())
    volume_path=manifest_path.parent/'pelvic-domain.npz'
    if sha(volume_path)!=domain['artifacts']['pelvic-domain.npz']:raise ValueError('Changed source volume')
    volume=np.load(volume_path,allow_pickle=False);raw=np.load(source/'trajectory.npz',allow_pickle=False)
    tri=raw['tissue_surface_triangles'];ids,reverse=np.unique(tri,return_inverse=True)
    tets=raw['tissue_tetrahedra'];labels=volume['material_index']
    face_owners={}
    for local in ((1,2,3),(0,3,2),(0,1,3),(0,2,1)):
        for face,label in zip(tets[:,local],labels):face_owners[tuple(sorted(face))]=int(label)
    face_labels=[face_owners[tuple(sorted(face))] for face in tri]
    source_garment=_resolve(config['panel']['source_path'])
    if sha(source_garment)!=config['panel']['source_sha256']:raise ValueError('Changed source garment')
    garment=next(g for g in json.loads(source_garment.read_text())['garments'] if g['id']=='shorts')
    source_nodes=config['panel']['source_node_indices'];source_faces=config['panel']['source_triangle_indices']
    # The frontend can compare these source arrays without reproducing Python JSON hashing.
    source_positions=np.asarray(garment['positions']).reshape(-1,3)[source_nodes]
    source_triangles=np.asarray(garment['indices']).reshape(-1,3)[source_faces]
    positions_a=raw['tissue_positions_m'][:,ids];positions_b=raw['panel_positions_m']
    if not np.array_equal(positions_b[0],source_positions):raise ValueError('Panel no longer corresponds to recorded garment nodes')
    quantization=max(float(np.max(np.abs(p-p.astype(np.float32).astype(float)))) for p in (positions_a,positions_b))
    if quantization>1e-7:raise ValueError('Display quantization exceeds 0.1 micrometre')
    def array(value):return np.asarray(value,dtype=np.float32).tolist()
    runtime_sources=config['source_sha256']
    source_hashes={str(manifest_path.relative_to(BASE)):sha(manifest_path),str(volume_path.relative_to(BASE)):sha(volume_path),
        str(source_garment.relative_to(BASE)):sha(source_garment),**{s['path']:s['source_sha256'] for s in domain['source_surfaces']}}
    source_hashes['app/src/clothing.js']=config['panel']['garment_constructor_sha256']
    source_hashes[config['panel']['source_skin']['path']]=config['panel']['source_skin']['geometry_sha256']
    payload={'schema_version':1,'kind':'computed_local_garment_tissue_contact','frame':domain['frame'],
        'time_s':raw['time_s'].tolist(),'replaced_body_entity_ids':list(config['material_owner_ids']),
        'tissue':{'positions_m':array(positions_a),'triangles':reverse.reshape(-1,3).tolist(),
            'triangle_material_index':face_labels,'material_regions':domain['material_regions'],
            'source_volume_node_indices':ids.tolist(),'contact_force_n':array(raw['tissue_contact_force_n'][:,ids])},
        'panel':{'positions_m':array(positions_b),'triangles':raw['panel_triangles'].tolist(),
            'contact_force_n':array(raw['panel_contact_force_n']),'source_node_indices':source_nodes,
            'source_triangle_indices':source_faces,'source_garment_id':'shorts',
            'source_positions_m':source_positions.tolist(),'source_triangles':source_triangles.tolist(),
            'source_positions_float64_le_sha256':_hash_array(garment['positions'],'<f8'),
            'source_indices_int64_le_sha256':_hash_array(garment['indices'],'<i8')},
        'force_observation':'contact impulse divided by preceding 0.005 s interval; initial frame is zero; N per node',
        'frames':[json.loads(line) for line in (source/'frames.jsonl').read_text().splitlines()],
        'configuration':config,'report':report,'source_hashes':source_hashes,'runtime_sources':runtime_sources,
        'artifacts':{'experiment_path':str(source.relative_to(BASE)),
            'experiment_report_sha256':sha(source/'report.json'),'original_volume_path':str(volume_path.relative_to(BASE)),
            'original_volume_sha256':sha(volume_path),'projection':'all original boundary faces; float32 display positions; full tetrahedra retained independently'}}
    out.mkdir(parents=True)
    encoded=json.dumps(payload,separators=(',',':'),allow_nan=False).encode()
    (out/'display.json.gz').write_bytes(gzip.compress(encoded,mtime=0))
    receipt={'schema_version':1,'status':'exported_computed_exploratory_fixture','display_path':'display.json.gz',
        'display_sha256':sha(out/'display.json.gz'),'bytes':(out/'display.json.gz').stat().st_size,
        'experiment_path':str(source.relative_to(BASE)),'frames':len(raw['time_s']),
        'tissue_boundary_nodes':len(ids),'tissue_boundary_triangles':len(tri),
        'panel_nodes':len(positions_b[0]),'panel_triangles':len(raw['panel_triangles']),
        'position_quantization_max_error_m':quantization,'source_hashes':source_hashes,'runtime_sources':runtime_sources,
        'exporter_sha256':sha(__file__),'whole_garment_containment_validated':False,'biological_calibration':False}
    (out/'manifest.json').write_text(json.dumps(receipt,indent=2)+'\n');return receipt

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('experiment');p.add_argument('output');a=p.parse_args()
    print(json.dumps(export_display(a.experiment,a.output),indent=2))
