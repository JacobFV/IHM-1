"""Verify full-resolution atlas projection against every authoritative source mesh."""
from pathlib import Path
import gzip,json
import numpy as np
import trimesh
from ihm.forge.acquisition import sha256

root=Path(__file__).resolve().parents[1]
manifest=json.loads((root/'data/derived/app/manifest.json').read_text())
model=next(m for m in manifest['models'] if m['id']=='bodyparts3d')
transform=model['display_transform'];rotation=np.array(transform['rotation']);translation=np.array(transform['translation']);scale=transform['scale']
assert np.allclose(rotation.T@rotation,np.eye(3),atol=1e-15) and np.isclose(np.linalg.det(rotation),1)
index=json.loads((root/'data/derived/anatomy/bodyparts3d_index.json').read_text())
structures={s['id']:s for s in manifest['structures']}
counts=dict(meshes=0,source_vertices=0,source_triangles=0,display_triangles=0,maximum_coordinate_error_m=0.,source_degenerate_triangles=0)
receipts=[]
for entry in index['meshes']:
    s=structures['bp3d-'+entry['element_id']];source=root/entry['source_path']
    assert sha256(source)==entry['sha256']==s['source']['sha256']
    geometry=root/'data/derived/app/geometry'/(s['id']+'.json.gz')
    assert sha256(geometry)==s['geometry_sha256']
    with gzip.open(geometry,'rt') as f:g=json.load(f)
    original=trimesh.load(source,force='mesh',process=False)
    vertices=np.array(g['positions']).reshape(-1,3);faces=np.array(g['indices']).reshape(-1,3)
    assert not g['display_decimation'] and len(faces)==entry['faces']
    assert np.array_equal(faces,original.faces),s['id']
    expected=np.array(original.vertices)@rotation.T*scale+translation
    assert vertices.shape==expected.shape and np.isfinite(vertices).all()
    error=float(np.abs(vertices-expected).max())
    assert error<5e-13,(s['id'],error)
    triangles=np.asarray(original.vertices)[np.asarray(original.faces)]
    degenerate=int(np.count_nonzero(np.linalg.norm(np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0]),axis=1)==0))
    counts['meshes']+=1;counts['source_vertices']+=len(original.vertices);counts['source_triangles']+=len(original.faces);counts['display_triangles']+=len(faces)
    counts['source_degenerate_triangles']+=degenerate;counts['maximum_coordinate_error_m']=max(error,counts['maximum_coordinate_error_m'])
    receipts.append(dict(id=s['id'],source_path=entry['source_path'],source_sha256=entry['sha256'],geometry_sha256=s['geometry_sha256'],triangles=len(faces),source_degenerate_triangles=degenerate))
assert counts['source_triangles']==counts['display_triangles'] and counts['meshes']==2234
report=dict(schema_version=1,summary=counts,receipts=receipts,transform=transform,limitations=['Topology and source coordinates reproduced; source anatomy is not independently measured here.','Degenerate source triangles are reported and retained, not silently repaired.','GPU float32 conversion is a display operation; JSON projection and authoritative source retain full precision.'])
out=root/'data/derived/anatomy/fidelity.json';out.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(counts))
