"""Prepare then explicitly run isolated retained-library wrap observation/candidate."""
from pathlib import Path
import argparse,csv,hashlib,json,os,signal,subprocess,tempfile,time
import xml.etree.ElementTree as ET
from materialize_stationary_wrap_probe import ROOT,materialize,SOURCE
from build_effective_potential_probe import verify
BUILD=ROOT/'data/derived/effective-potential-build-3a9juno_/manifest.json'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def prepare():
    build=verify(BUILD);out=Path(tempfile.mkdtemp(prefix='stationary-wrap-candidate-',dir=ROOT/'data/derived'));materialize(out/'models')
    for n in ('native_stationary_wrap_probe.cpp','arm26_wrap_candidate.h','arm26_stationary_ellipsoid.h'):(out/n).write_bytes((ROOT/'scripts'/n).read_bytes())
    case=json.loads((ROOT/build['reference_manifest']).read_text());seedpath=ROOT/case['old_q_seed'];seed=json.loads(seedpath.read_text())['coordinates']
    raw=(ROOT/SOURCE).read_bytes();bounds={c.get('name'):list(map(float,c.findtext('range').split())) for c in ET.fromstring(raw).iter('Coordinate')}
    coordinates=[stem+s for s in ('r','l') for stem in ('arm_flex_','arm_add_','arm_rot_','elbow_flex_')]
    cases=[]
    for c in coordinates:
        for h in (-1e-3,-1e-4,-1e-5,0,1e-5,1e-4,1e-3):
            q=dict(seed);q[c]+=h
            if bounds[c][0]<=q[c]<=bounds[c][1]:cases.append((f'local_{c}_{h:g}',c,q))
        for t in (0,.25,.5,.75,1):
            q=dict(seed);q[c]=bounds[c][0]+t*(bounds[c][1]-bounds[c][0]);cases.append((f'range_{c}_{t:g}',c,q))
    # Held out combined perturbations; no parameter fitting uses these points.
    for i,c in enumerate(coordinates):
        q=dict(seed)
        for j,n in enumerate(coordinates):q[n]=max(bounds[n][0],min(bounds[n][1],q[n]+(.017 if (i+j)%2 else -.023)))
        cases.append((f'heldout_{i}',c,q))
    (out/'cases.txt').write_text(''.join(label+' '+c+' '+' '.join(n+' '+format(v,'.17g') for n,v in q.items())+'\n' for label,c,q in cases))
    command=[str(out/'native_stationary_wrap_probe.cpp') if s.endswith('/native_probe.cpp') else str(out/'probe') if s==str(ROOT/build['executable']) else s for s in build['command']]
    files={str(p.relative_to(ROOT)):sha(p) for p in [out/'native_stationary_wrap_probe.cpp',out/'arm26_wrap_candidate.h',out/'arm26_stationary_ellipsoid.h',out/'cases.txt',*(out/'models').iterdir(),seedpath,ROOT/'scripts/run_stationary_wrap_probe.py',ROOT/'scripts/materialize_stationary_wrap_probe.py',ROOT/'scripts/materialize_arm26_wrap_probe.py']}
    m={'schema':'ihm.arm26-wrap-candidate-probe.v1','base_build':str(BUILD.relative_to(ROOT)),'base_build_sha256':sha(BUILD),'archive_manifest':build['archive_manifest'],'files':files,'command':command,'cases':len(cases),'coordinate_bounds':{c:bounds[c] for c in coordinates},'native_modified':False,'physical_time_advanced_s':0,'candidate_scope':'Private BRA cylinders and stationary BIClong ellipsoid geodesics; same native-seeded branches or explicit rejection; other94 paths untouched.'}
    (out/'manifest.json').write_text(json.dumps(m,indent=2)+'\n');print(out,flush=True);return out

def run(out):
    out=Path(out).resolve();m=json.loads((out/'manifest.json').read_text());verify(ROOT/m['base_build'])
    assert sha(ROOT/m['base_build'])==m['base_build_sha256']
    for p,h in m['files'].items():assert sha(ROOT/p)==h,p
    started=time.monotonic();archive=ROOT/m['archive_manifest'];ar=json.loads(archive.read_text());libs=archive.parent/'libraries'
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',LD_LIBRARY_PATH=str(libs))
    def execute(command,logname):
        with (out/logname).open('w') as log:
            p=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,env=env,start_new_session=True)
            try:
                result=p.wait(timeout=max(1,min(60,90-(time.monotonic()-started))))
                if result:raise RuntimeError(f'Native step failed: {logname} exit{result}')
            finally:
                if p.poll() is None:os.killpg(p.pid,signal.SIGKILL)
                p.wait()
    execute(m['command'],'compile.log')
    for mode in ('observed','candidate'):
        execute(['prlimit','--as=4294967296','--','nice','-n','10',str(archive.parent/ar['dynamic_loader']),'--library-path',str(libs),str(out/'probe'),str(out/'models'/f'{mode}.osim'),str(BUILD.parent/'inputs/subject_walk_scaled_FunctionBasedPathSet.xml'),str(out/'cases.txt'),str(out/f'{mode}.csv'),str(out/f'{mode}.jsonl')],mode+'.log')
    assert all('cylinder_boundary_rejections=4' in (out/(mode+'.log')).read_text() for mode in ('observed','candidate'))
    rows={mode:list(csv.DictReader((out/f'{mode}.csv').open())) for mode in ('observed','candidate')}
    assert all(len(v)==m['cases']*98 for v in rows.values())
    baseline={(r['case'],r['muscle']):r for r in rows['observed']};errors=[];reject=[];unchanged=0
    for r in rows['candidate']:
        b=baseline[r['case'],r['muscle']]
        if not r['muscle'].startswith(('arm26_BRA_','arm26_BIClong_')):
            assert r['status']==b['status']
            if r['status']=='ok':assert abs(float(r['length'])-float(b['length']))<1e-12
            unchanged+=1
        if r['status']!='ok':reject.append({'case':r['case'],'muscle':r['muscle']});continue
        if r['muscle'].startswith(('arm26_BRA_','arm26_BIClong_')):
            errors.append(abs(float(r['moment_arm'])-float(r['fd'])))
            assert float(r['unit_resultant_force'])<1e-10 and float(r['unit_resultant_moment'])<1e-10
            assert abs(float(r['lengthening_speed'])+.07*float(r['fd']))<1e-8
            assert abs(float(r['unit_body_power'])+float(r['lengthening_speed']))<1e-10
    for p,h in m['files'].items():assert sha(ROOT/p)==h,p
    verify(ROOT/m['base_build'])
    result={'status':'analyzed','output':str(out.relative_to(ROOT)),'cases_per_variant':m['cases'],'unchanged_other_muscle_observations':unchanged,'accepted_corrected_wrap_observations':len(errors),'rejected_observations':reject,'maximum_accepted_wrap_virtual_work_error_m':max(errors,default=None),'candidate_gradient_gate_passed':bool(errors) and max(errors)<1e-6,'whole_source_range_certified':False,'native_time_advanced_s':0,'wall_s':time.monotonic()-started,'executable_sha256':sha(out/'probe')}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--prepare',action='store_true');p.add_argument('--run',type=Path);a=p.parse_args()
    if a.prepare:prepare()
    elif a.run:run(a.run)
    else:raise SystemExit('Select source preparation or coordinated native run')
