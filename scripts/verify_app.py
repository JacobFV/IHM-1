"""Exercise local API and its asset/run boundaries over real HTTP."""
from pathlib import Path
import json
import threading
from urllib.request import urlopen,Request
from urllib.error import HTTPError
from ihm.app import create_server

server=create_server(Path(__file__).resolve().parents[1],port=0)
thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
base=f'http://127.0.0.1:{server.server_port}'
def get(path):
 with urlopen(base+path) as r:return json.load(r)
def rejected(path,data=None,headers=None):
 try:urlopen(Request(base+path,data=None if data is None else json.dumps(data).encode(),headers=headers or {'Content-Type':'application/json'}))
 except HTTPError as e:return e.code
 raise AssertionError('invalid request accepted')
try:
 try:create_server(Path(__file__).resolve().parents[1],port=server.server_port)
 except OSError:pass
 else:raise AssertionError('duplicate bind unexpectedly succeeded')
 m=get('/api/manifest');assert len(m['structures'])>=2234
 # A materialization draws a subset of the implicit model, and every tier's membership
 # is derived from the manifest and the colocation audits rather than written down.
 tiers=get('/api/body/materializations');ids={s['id'] for s in m['structures'] if s['model_id']=='ihm-body'}
 assert len(tiers['tiers'])>=3 and any(t['value']==tiers['default'] for t in tiers['tiers'])
 assert all(set(t['omits'])<=ids and t['structures']==len(ids)-len(t['omits']) for t in tiers['tiers'])
 assert tiers['tiers'][-1]['omits']==[] and tiers['tiers'][-1]['structures']==len(ids)
 assert tiers['tiers'][0]['structures']<tiers['tiers'][-1]['structures']
 g=get(m['structures'][0]['geometry_url']);assert len(g['positions'])>10 and len(g['indices'])>3
 coverage=get('/api/anatomy/coverage');assert coverage
 assert get('/api/anatomy/fidelity')['summary']['source_triangles']==6681030
 lymph=get('/api/lymphatic');assert len(lymph['nodes'])==996 and len(lymph['edges'])==1117
 assert lymph['directed'] is False and all(e['flow_ml_s'] is None for e in lymph['edges'])
 thermal=get('/api/thermal/index');assert len(thermal['runs'])==4
 assert get('/api/thermal?run=lying_default')['configuration']['posture']=='lying'
 assert rejected('/api/thermal?run=../../secret')==404
 assert get('/api/native-targets')['runs']
 p=get('/api/physiology');assert len(p['time_s'])==3000 and 'ArterialPressure(mmHg)' in p['values']
 assert rejected('/api/geometry/../../pyproject.toml') in (400,404)
 assert rejected('/%2e%2e/pyproject.toml') in (400,404)
 assert rejected('/api/scenarios',{'seconds':float('nan')})==400
 assert rejected('/api/scenarios',{'seconds':10000})==400
 assert rejected('/api/scenarios',{'state_path':'/etc/passwd'})==400
 assert rejected('/api/scenarios',{'patient':'../../secret'})==400
 assert rejected('/api/scenarios',{'seconds':60},headers={'Content-Type':'application/json','Origin':'https://untrusted.example'})==403
 assert get('/api/scenarios')['available'] is True
 print('verified: local HTTP assets/physiology, traversal/input/origin rejection, backend discovery')
finally:server.shutdown();server.server_close()
