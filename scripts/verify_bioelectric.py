"""Audit exported actual native tissue values and geometric cell identity."""
from pathlib import Path
import gzip,json,hashlib
import numpy as np
root=Path(__file__).resolve().parents[1]
m=json.loads((root/'data/derived/app/manifest.json').read_text())
s=next(s for s in m['structures'] if s['id']=='betse-tissue-cells')
p=root/'data/derived/app/geometry/betse-tissue-cells.json.gz'
assert hashlib.sha256(p.read_bytes()).hexdigest()==s['geometry_sha256']
g=json.loads(gzip.decompress(p.read_bytes()));v=np.array(g['positions']).reshape(-1,3);f=np.array(g['indices']).reshape(-1,3);ids=np.array(g['cell_ids'])
assert len(v)==len(ids) and set(ids)==set(range(212))
assert all(len(set(row))==1 for row in ids[f])
assert np.isfinite(v).all() and np.all(v[:,2]==0)
assert np.all(np.linalg.norm(np.cross(v[f[:,1]]-v[f[:,0]],v[f[:,2]]-v[f[:,0]]),axis=1)>1e-12)
assert len(g['times'])==34 and np.all(np.diff(g['times'])>0)
for field in g['scalar_fields']:
 a=np.array([frame['values'][field['id']] for frame in g['frames']])
 assert a.shape==(34,212) and np.isfinite(a).all()
 assert [float(a.min()),float(a.max())]==field['range']
assert np.isclose(np.mean(g['frames'][-1]['values']['Vmem']),-.04312935619067724)
assert g['scalar_fields'][0]['unit']=='V'
assert s['source']['specimen']=='generic computational tissue; no human subject'
print('verified 212 native planar cells, 34 actual field frames, cell identity per triangle, bounds and Vmem units')
mem=np.array(g['source_membrane_to_cell']);count=np.bincount(mem,minlength=212)
expected=np.array([np.bincount(mem,weights=x,minlength=212)/count for x in g['source_membrane_voltages_V']])
actual=np.array([x['values']['Vmem'] for x in g['frames']])
assert np.allclose(actual,expected,rtol=1e-12,atol=1e-14)
