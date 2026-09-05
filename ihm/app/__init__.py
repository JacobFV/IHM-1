"""Local scientific workbench API. Source assets are served only by declared IDs."""
from concurrent.futures import ThreadPoolExecutor
from http.server import BaseHTTPRequestHandler,ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlparse,parse_qs,unquote
import gzip
import hashlib
import json
import mimetypes
import re
import socket
import threading
import time
import uuid

SAFE_ID=re.compile(r'^[A-Za-z0-9_-]+$')

def read_json(path):return json.loads(Path(path).read_text())

class Jobs:
    def __init__(self,root):
        self.root=root;self.pool=ThreadPoolExecutor(max_workers=1,thread_name_prefix='ihm-native')
        self.lock=threading.Lock();self.runs=[];self.variant_cache={}
        self.directory=root/'data/derived/scenarios';self.directory.mkdir(parents=True,exist_ok=True)
        for file in sorted(self.directory.glob('*/job.json')):
            try:
                job=read_json(file)
                if job['status'] in ('running','queued'):job.update(status='interrupted',error='Server stopped before this run completed')
                self.runs.append(job)
            except (ValueError,KeyError):continue
    def list(self):
        from ihm.native import RUNTIME,available_patients
        with self.lock:runs=[dict(r) for r in self.runs]
        variants=[{'id':'upstream','label':'Original upstream source','available':True,'status':'source reference; known female initialization bounds defect'}]
        for name,label in [('saturation_bounds','Saturation bounds correction'),('saturation_bounds_heatflux','Bounds and evaporation telemetry corrections')]:
            directory=RUNTIME/'variants'/name;manifest=directory/'manifest.json';library=directory/'libbiogears.so.8.0.0'
            if manifest.is_file() and library.is_file():
                metadata=read_json(manifest)
                stamp=(library.stat().st_mtime_ns,library.stat().st_size,metadata['library_sha256'])
                cached=self.variant_cache.get(name)
                if cached is None or cached[0]!=stamp:
                    with library.open('rb') as f:matches=hashlib.file_digest(f,'sha256').hexdigest()==metadata['library_sha256']
                    self.variant_cache[name]=(stamp,matches)
                else:matches=cached[1]
                variants.append(dict(id=name,label=label,available=matches,status='local source patch; execution regression checked, not independent clinical validation',scope=metadata['scope'],library_sha256=metadata['library_sha256']))
        preferred='saturation_bounds_heatflux' if any(v['id']=='saturation_bounds_heatflux' and v['available'] for v in variants) else 'upstream'
        return {'available':(RUNTIME/'native_biogears_rest').is_file(),'patients':available_patients(),'runs':runs,'engine_variants':variants,'default_engine_variant':preferred,
                'limits':{'max_seconds':600,'max_pending':4,'parallel_runs':1}}
    def submit(self,data,canonical=False):
        from ihm.native import NativeConfig
        if not isinstance(data,dict) or 'state_path' in data:raise ValueError('State paths are not accepted by the web API')
        if canonical:
            from ihm.assembly.body import CanonicalBody
            body=CanonicalBody.from_workspace(self.root)
            data=dict(data)
            patient=body.payload['profile']['native_patient']
            if data.get('patient',patient)!=patient:raise ValueError('Canonical scenarios require the shared generic body profile')
            data['patient']=patient
            body_interventions=data.pop('body_interventions',[])
        elif 'body_interventions' in data:raise ValueError('Body interventions require the canonical body')
        config=NativeConfig.from_dict(data)
        if canonical:
            from ihm.assembly.body_protocol import validate_interventions
            body_interventions=validate_interventions(body_interventions,body.assets['peripheral'],config.seconds)
        if config.seconds>600:raise ValueError('Web scenarios are bounded to 600 seconds')
        if len(config.interventions)>20:raise ValueError('Web scenarios permit at most 20 actions')
        from dataclasses import asdict
        with self.lock:
            if sum(r['status'] in ('queued','running') for r in self.runs)>=4:raise ValueError('Scenario queue full; wait for a run to finish')
            job={'id':time.strftime('%Y%m%d-%H%M%S')+'-'+uuid.uuid4().hex[:8],'status':'queued','config':asdict(config),'created_unix':time.time(),'canonical_body':canonical}
            if canonical:job['canonical_sources']=body.payload['sources'];job['canonical_runtime_sources']=body.payload['runtime_sources'];job['canonical_patient_sha256']=body.payload['profile']['native_patient_sha256'];job['body_interventions']=body_interventions
            self.runs.append(job);self._save(job)
        self.pool.submit(self._execute,job,config)
        return dict(job)
    def _save(self,job):
        p=self.directory/job['id'];p.mkdir(exist_ok=True)
        temporary=p/'job.partial.json';temporary.write_text(json.dumps(job,allow_nan=False));temporary.replace(p/'job.json')
    def _execute(self,job,config):
        from ihm.native import run_native
        with self.lock:job['status']='running';self._save(job)
        try:
            if job.get('canonical_body'):
                from ihm.assembly.body import CanonicalBody
                body=CanonicalBody.from_workspace(self.root)
                if body.payload['sources']!=job['canonical_sources'] or body.payload['runtime_sources']!=job['canonical_runtime_sources'] or body.payload['profile']['native_patient_sha256']!=job['canonical_patient_sha256']:raise ValueError('Canonical body changed after scenario submission; submit a new run')
            summary=run_native(config,self.directory/job['id']/'output')
            if job.get('canonical_body'):
                from ihm.assembly.body import CanonicalBody
                current=CanonicalBody.from_workspace(self.root)
                if current.payload['sources']!=job['canonical_sources'] or current.payload['runtime_sources']!=job['canonical_runtime_sources']:raise ValueError('Canonical body changed during native execution; original native output retained')
                trajectory=body.simulate(self.directory/job['id']/'output',self.directory/job['id']/'body-trajectory.json',body_interventions=job.get('body_interventions',[]))
                summary['canonical_body']={'frames':len(trajectory['frames']),'clock':trajectory['clock'],'audit':trajectory['audit']}
            with self.lock:job.update(status='completed',summary=summary);self._save(job)
        except Exception as e:
            with self.lock:job.update(status='failed',error=str(e));self._save(job)

class Server(ThreadingHTTPServer):
    daemon_threads=True
    def manifest(self):
        file=self.root/'data/derived/app/manifest.json';stamp=file.stat().st_mtime_ns
        with self.manifest_lock:
            if getattr(self,'manifest_stamp',None)!=stamp:
                self.manifest_data=read_json(file)
                self.structure_ids={s['id'] for s in self.manifest_data['structures']}
                self.manifest_stamp=stamp
            return self.manifest_data,self.structure_ids
    def server_close(self):
        super().server_close()
        if hasattr(self,'jobs'):self.jobs.pool.shutdown(wait=False,cancel_futures=True)

def create_server(root=None,port=8765,host='127.0.0.1'):
    if host not in ('127.0.0.1','localhost','::1'):raise ValueError('Workbench binds loopback only')
    root=Path(root or Path(__file__).resolve().parents[2]).resolve()
    class Handler(BaseHTTPRequestHandler):
        def log_message(self,format,*args):pass
        def _authorized(self,post=False):
            hostname=urlparse('http://'+self.headers.get('Host','')).hostname
            if hostname not in ('127.0.0.1','localhost','::1'):return False
            origin=self.headers.get('Origin')
            if origin:
                parsed=urlparse(origin)
                if parsed.scheme!='http' or parsed.hostname not in ('127.0.0.1','localhost','::1'):return False
            return True
        def _send(self,data,status=200,content_type='application/json',encoding=None):
            if not isinstance(data,bytes):data=json.dumps(data,allow_nan=False,separators=(',',':')).encode()
            if not encoding and len(data)>2000 and 'gzip' in self.headers.get('Accept-Encoding',''):
                data=gzip.compress(data,compresslevel=1,mtime=0);encoding='gzip'
            self.send_response(status);self.send_header('Content-Type',content_type)
            self.send_header('Content-Length',str(len(data)));self.send_header('X-Content-Type-Options','nosniff')
            self.send_header('Cache-Control','no-cache')
            if encoding:self.send_header('Content-Encoding',encoding)
            self.end_headers()
            try:self.wfile.write(data)
            except (BrokenPipeError,ConnectionResetError):pass
        def _error(self,message,status=400):self._send({'error':message},status)
        def do_GET(self):
            if not self._authorized():return self._error('Only local workbench requests are accepted',403)
            parsed=urlparse(self.path);path=unquote(parsed.path);query=parse_qs(parsed.query)
            if '..' in path.split('/') or '\\' in path or '\x00' in path:return self._error('Invalid path',400)
            try:
                derived=root/'data/derived'
                if path=='/api/manifest':return self._send(self.server.manifest()[0])
                if path=='/api/body':
                    from ihm.assembly.body import CanonicalBody
                    return self._send(CanonicalBody.from_workspace(root).describe())
                if path=='/api/body/coverage':
                    return self._send(read_json(derived/'audits/execution-coverage.json'))
                if path.startswith('/api/body/experiments/'):
                    from ihm.app.experiments import read_experiment
                    return self._send(read_experiment(root,path.removeprefix('/api/body/experiments/')))
                if path=='/api/body/peripheral':
                    from ihm.assembly.body import CanonicalBody
                    return self._send(CanonicalBody.from_workspace(root).assets['peripheral'])
                if path=='/api/body/certainty':
                    from ihm.assembly.body import CanonicalBody
                    return self._send(CanonicalBody.from_workspace(root).certainty(query.get('entity',[''])[0]))
                if path=='/api/body/spectra':
                    from ihm.assembly.body import CanonicalBody,read_native
                    spectra=read_json(derived/'canonical/trajectory-spectra.json');body=CanonicalBody.from_workspace(root)
                    if spectra.get('canonical_sources')!=body.payload['sources'] or spectra.get('runtime_sources')!=body.payload['runtime_sources']:return self._error('Canonical spectra are stale; rematerialize the native run',409)
                    native=read_native(spectra['native_directory'])
                    if spectra.get('native_summary_sha256')!=native['input_hashes']['summary.json'] or any(native['input_hashes'][key]!=value for key,value in spectra['source_files'].items()):return self._error('Native spectra inputs changed; rematerialize',409)
                    return self._send(spectra)
                if path=='/api/body/trajectory':
                    run=query.get('run',[None])[0]
                    if run:
                        if not SAFE_ID.fullmatch(run):raise ValueError('Invalid run ID')
                        found=next((j for j in self.server.jobs.list()['runs'] if j['id']==run and j['status']=='completed' and j.get('canonical_body')),None)
                        if not found:return self._error('Completed canonical run not found',404)
                        file=derived/'scenarios'/run/'body-trajectory.json'
                    else:file=derived/'canonical/trajectory.json'
                    raw=file.read_bytes();trajectory=json.loads(raw)
                    from ihm.assembly.body import CanonicalBody
                    body=CanonicalBody.from_workspace(root)
                    if trajectory.get('runtime_sources')!=body.payload['runtime_sources'] or any(trajectory['sources'].get(key)!=value for key,value in body.payload['sources'].items()):return self._error('Canonical trajectory is stale; rematerialize the native run',409)
                    from ihm.assembly.body import read_native
                    native=read_native(trajectory['sources']['native_run']['directory'])
                    if native['input_hashes']!=trajectory['sources']['native_run']['files']:return self._error('Native trajectory inputs changed; rematerialize',409)
                    view=query.get('view',['full'])[0]
                    if view not in ('full','display'):raise ValueError('Unknown trajectory projection')
                    if view=='display':
                        from ihm.app.projection import display_trajectory
                        trajectory=display_trajectory(trajectory,source_sha256=hashlib.sha256(raw).hexdigest())
                    return self._send(trajectory)
                if path.startswith('/api/geometry/'):
                    id=path.removeprefix('/api/geometry/')
                    if not SAFE_ID.fullmatch(id):return self._error('Invalid structure ID')
                    if id not in self.server.manifest()[1]:return self._error('Unknown structure',404)
                    file=(derived/'app/geometry'/f'{id}.json.gz').resolve()
                    if not file.is_relative_to((derived/'app/geometry').resolve()):return self._error('Invalid asset',400)
                    compressed=file.read_bytes()
                    if 'gzip' in self.headers.get('Accept-Encoding',''):return self._send(compressed,encoding='gzip')
                    return self._send(gzip.decompress(compressed))
                if path=='/api/physiology':
                    from ihm.native import load_trajectory
                    run=query.get('run',[None])[0]
                    if run:
                        if not SAFE_ID.fullmatch(run):raise ValueError('Invalid run ID')
                        found=next((j for j in self.server.jobs.list()['runs'] if j['id']==run and j['status']=='completed'),None)
                        if not found:return self._error('Completed run not found',404)
                        directory=derived/'scenarios'/run/'output'
                    else:directory=derived/'physiology/biogears_native_run'
                    data=load_trajectory(directory);data['metadata']={'source_kind':'native_simulation','source':'BioGears','independently_validated':False}
                    summary=directory/'summary.json'
                    if summary.exists():data['metadata']['execution']=read_json(summary)
                    return self._send(data)
                if path=='/api/temporal':return self._send(read_json(derived/'temporal/index.json'))
                if path=='/api/calibration':return self._send(read_json(derived/'calibration/index.json'))
                if path=='/api/human':
                    from ihm.human import ImplicitHuman
                    return self._send(ImplicitHuman.open(root).describe())
                if path=='/api/anatomy/coverage':return self._send(read_json(derived/'anatomy/anatomy_coverage.json'))
                if path=='/api/anatomy/fidelity':return self._send(read_json(derived/'anatomy/fidelity.json'))
                if path=='/api/lymphatic':return self._send(read_json(derived/'lymphatic/graph.json'))
                if path=='/api/native-targets':return self._send(read_json(derived/'calibration/native-target-audit.json'))
                if path=='/api/thermal/index':return self._send(read_json(derived/'thermal/index.json'))
                if path=='/api/thermal':
                    index=read_json(derived/'thermal/index.json');run=query.get('run',[None])[0]
                    selected=next((m for m in index['models'] if m['id']==run),None) if run else index['models'][0]
                    if selected is None:return self._error('Unknown thermal run',404)
                    file=(root/selected['trajectory_path']).resolve()
                    if not file.is_relative_to((derived/'thermal').resolve()):return self._error('Invalid thermal asset',400)
                    return self._send(read_json(file))
                if path=='/api/csf/index':return self._send(read_json(derived/'csf/index.json'))
                if path=='/api/csf':
                    run=query.get('run',['baseline'])[0]
                    if run not in ('baseline','native_map_driven','hypotension'):return self._error('Unknown CSF run',404)
                    return self._send(read_json(derived/'csf'/f'{run}.json'))
                if path=='/api/bioelectric':return self._send(read_json(derived/'bioelectric/tissue.json'))
                if path=='/api/reproductive':return self._send(read_json(derived/'reproductive/trajectory.json'))
                if path=='/api/vascular/audit':return self._send(read_json(derived/'vascular/audit.json'))
                if path=='/api/vascular/cap-flow':return self._send(read_json(derived/'vascular/cap-flow.json'))
                if path=='/api/coupling/skin-lymph':return self._send(read_json(derived/'coupling/native-skin-circuit.json'))
                if path=='/api/coupling/hydrostatic':return self._send(read_json(derived/'coupling/hydrostatic.json'))
                if path=='/api/coverage':return self._send(read_json(derived/'system-coverage.json'))
                if path=='/api/coupling':return self._send(read_json(derived/'coupling/index.json'))
                if path=='/api/evidence':
                    from ihm.forge.catalog import EvidenceCatalog
                    catalog=EvidenceCatalog(derived/'evidence-catalog.sqlite')
                    if 'term' in query:
                        term=query['term'][0]
                        if len(term)>200:raise ValueError('Search term too long')
                        return self._send({'records':catalog.search(term,30,query.get('source',[None])[0])})
                    return self._send(catalog.summary())
                if path=='/api/scenarios':return self._send(self.server.jobs.list())
                if path.startswith('/api/'):return self._error('Endpoint not found',404)
                static=(root/'app/dist').resolve();file=(static/(path.lstrip('/') or 'index.html')).resolve()
                if not file.is_relative_to(static) or not file.is_file():return self._error('Asset not found; build app first',404)
                return self._send(file.read_bytes(),content_type=mimetypes.guess_type(file.name)[0] or 'application/octet-stream')
            except FileNotFoundError:return self._error('Required local artifact unavailable; run the corresponding build script',404)
            except (ValueError,KeyError,TypeError) as e:return self._error(str(e),400)
        def do_POST(self):
            if not self._authorized(post=True):return self._error('Only local workbench requests are accepted',403)
            if self.path not in ('/api/scenarios','/api/body/scenarios'):return self._error('Endpoint not found',404)
            if self.headers.get('Content-Type','').split(';')[0]!='application/json':return self._error('Expected application/json',415)
            try:
                length=int(self.headers.get('Content-Length','0'))
                if not 0<length<=32768:raise ValueError('JSON request must be 1–32768 bytes')
                data=json.loads(self.rfile.read(length),parse_constant=lambda v:(_ for _ in ()).throw(ValueError('Nonfinite JSON number')))
                if not self.server.jobs.list()['available']:return self._error('Native backend unavailable; build it locally',503)
                return self._send(self.server.jobs.submit(data,canonical=self.path=='/api/body/scenarios'),202)
            except FileNotFoundError:return self._error('Required canonical/native artifact unavailable; build the body first',503)
            except (ValueError,TypeError,KeyError) as e:return self._error(str(e),400)
    server_class=type('IPv6Server',(Server,),{'address_family':socket.AF_INET6}) if host=='::1' else Server
    server=server_class((host,port),Handler);server.root=root;server.jobs=Jobs(root);server.manifest_lock=threading.Lock();return server

def serve(root=None,port=8765):
    server=create_server(root,port)
    print(f'IHM workbench: http://127.0.0.1:{server.server_port}',flush=True)
    try:server.serve_forever()
    except KeyboardInterrupt:pass
    finally:server.server_close()
