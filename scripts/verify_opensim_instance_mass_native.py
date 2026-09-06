"""Isolated actual92-muscle baseline/variant mass fixture; no latest promotion."""
from pathlib import Path
import argparse,hashlib,json,os,resource,selectors,shlex,shutil,subprocess,sys,tempfile,time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.assembly.articulated import CanonicalRegistration
from ihm.native.mechanical_stream import SOURCE_FILES
from build_simbody_instance_mass_variant import validate
ROOT=Path(__file__).resolve().parents[1]
ENERGY=('muscle_metabolic_energy_j','signed_active_fiber_work_j','muscle_heat_energy_j','positive_active_fiber_work_j','external_work_j')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
class Session:
    def __init__(self,executable,inputs,out,env):
        out.mkdir();self.log=(out/'engine.log').open('w');self.other=(out/'other_stdout.log').open('w');self.commands=(out/'commands.txt').open('w');self.buffer=b''
        self.process=subprocess.Popen(['prlimit','--as=4294967296','--','nice','-n','10',str(executable),str(inputs),str(out),'supine','77.6122029'],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=self.log,env=env,cwd=out,bufsize=0)
        self.selector=selectors.DefaultSelector();self.selector.register(self.process.stdout,selectors.EVENT_READ)
        try:self.initial=self.read()
        except BaseException:self.close();raise
    def read(self):
        while True:
            while b'\n' not in self.buffer:
                if not self.selector.select(15):raise TimeoutError('bounded native response timeout')
                raw=os.read(self.process.stdout.fileno(),65536)
                if not raw:raise RuntimeError('native stream exited')
                self.buffer+=raw
            line,self.buffer=self.buffer.split(b'\n',1);line=line.decode()
            if not line.startswith('@IHM '):self.other.write(line+'\n');continue
            result=json.loads(line[5:])
            if 'error' in result:raise ValueError(result['error'])
            return result
    def request(self,line):self.commands.write(line+'\n');self.commands.flush();self.process.stdin.write((line+'\n').encode());self.process.stdin.flush();return self.read()
    def close(self):
        try:
            if self.process.poll() is None:self.process.stdin.write(b'close\n');self.process.stdin.flush();self.process.wait(timeout=3)
        except (OSError,subprocess.TimeoutExpired):self.process.kill();self.process.wait()
        finally:self.selector.close();self.process.stdin.close();self.process.stdout.close();self.log.close();self.other.close();self.commands.close()
def compare(a,b,path=''):
    if isinstance(a,dict):
        assert set(a)==set(b),path
        return max([compare(a[k],b[k],path+'/'+k) for k in a] or [0.])
    if isinstance(a,list):
        assert len(a)==len(b),path
        return max([compare(x,y,path) for x,y in zip(a,b)] or [0.])
    if isinstance(a,(float,int)) and not isinstance(a,bool):
        assert np.isfinite(a) and np.isfinite(b),path
        error=abs(a-b);assert error<=2e-10*max(1.,abs(a),abs(b)),(path,a,b)
        return error
    assert a==b,(path,a,b);return 0.
def main():
    p=argparse.ArgumentParser();p.add_argument('--variant',type=Path,default=Path('data/runtime/opensim/variants/instance_mass_v1'));p.add_argument('--run-native',action='store_true');a=p.parse_args()
    if not a.run_native:raise SystemExit('Requires coordinated compile + baseline/variant native slot')
    variant=(ROOT/a.variant).resolve();vm=json.loads((variant/'manifest.json').read_text());assert vm['complete'];validate(variant,vm);assert sha(variant/'libSimTKsimbody.so.3.9')==vm['library_sha256']
    out=Path(tempfile.mkdtemp(prefix='opensim-instance-mass-',dir=ROOT/'data/derived'));build=out/'build';build.mkdir();started=time.monotonic()
    names=('native_mechanical_stream.cpp','native_muscle_metabolism.h','native_static_pose.h','native_surface_foundation.h','native_bed_compression.h','native_opensim_mass_validation.h')
    frozen={str((ROOT/'scripts'/n).relative_to(ROOT)):sha(ROOT/'scripts'/n) for n in names}
    for n in names:(build/n).write_bytes((ROOT/'scripts'/n).read_bytes())
    source=build/'native_mechanical_stream.cpp';text=source.read_text();anchor='   if(command=="close")break;'
    assert text.count(anchor)==1
    text=text.replace('#include "native_surface_foundation.h"','#include "native_surface_foundation.h"\n#include "native_opensim_mass_validation.h"')
    text=text.replace(anchor,anchor+'\n   if(command=="mass_probe"){const auto result=ihm_mass_validation::probe(model,state);std::cout<<"@IHM "<<result<<std::endl;continue;}\n   if(command=="mass_transfer"){const auto result=ihm_mass_validation::transfer(model,state,in);std::cout<<"@IHM "<<result<<std::endl;continue;}')
    text=text.replace('state.invalidateAllCacheAtOrAbove(SimTK::Stage::Dynamics);emit("restored")','state.invalidateAllCacheAtOrAbove(SimTK::Stage::Instance);emit("restored")')
    text=text.replace('heat_energy=before.heat_energy;std::ostringstream o;', 'heat_energy=before.heat_energy;state.invalidateAllCacheAtOrAbove(SimTK::Stage::Instance);model.markControlsAsInvalid(state);std::ostringstream o;');source.write_text(text)
    runtime=ROOT/'data/runtime/opensim';flags=['-std=c++20','-O0','-DSWIG_PYTHON']
    for x in ('install/opensim/include','install/opensim/include/OpenSim','install/simbody/include/simbody'):flags+=['-isystem',str(runtime/x)]
    previous=shlex.split((runtime/'dynamics-adapter/build/CMakeFiles/native_opensim_dynamics.dir/link.txt').read_text());libraries=[v for v in previous if v.endswith('.so') or '.so.' in v or v.startswith('-Wl,') or v.startswith('-l')]
    for path in libraries:
        if Path(path).is_file():frozen[str(Path(path).relative_to(ROOT))]=sha(path)
    command=['prlimit','--as=4294967296','--','nice','-n','10','c++',*flags,str(source),'-o',str(build/'fixture'),*libraries,'-ldl'];dump(build/'command.json',command)
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    with (build/'compile.log').open('w') as log:r=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,env=env)
    if r.returncode:raise RuntimeError('fixture compile failed: '+str(build/'compile.log'))
    assert all(sha(ROOT/path)==digest for path,digest in frozen.items()),'frozen input changed during compile'
    inputs=out/'inputs';inputs.mkdir();registration_path=ROOT/'data/derived/mechanics/whole_body_arm26_v2/registration.json';registration=json.loads(registration_path.read_text());original=ROOT/'data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/example3DWalking'
    assert sha(ROOT/registration['catalog_path'])==registration['catalog_sha256']
    assert all(sha(ROOT/name)==digest for name,digest in registration['sources'].items())
    catalog=json.loads((ROOT/registration['catalog_path']).read_text())
    for n in SOURCE_FILES:(inputs/n).write_bytes(((ROOT/registration['model_path']) if n=='subject_walk_scaled.osim' else original/n).read_bytes())
    assert sha(inputs/'subject_walk_scaled.osim')==registration['model_sha256']
    input_hashes={n:sha(inputs/n) for n in SOURCE_FILES};dump(out/'input_identity.json',{'registration_sha256':sha(registration_path),'files':input_hashes,'canonical_sha256':sha(ROOT/'data/derived/canonical/mechanics.json')})
    libdirs=[runtime/'install/opensim/lib',runtime/'install/simbody/lib',*[runtime/'sysroot/usr/lib/aarch64-linux-gnu'/x for x in ('lapack','blas','')]]
    baseline={};variant_data={};loaded_hashes={};site=None
    for label,paths in [('baseline',libdirs),('variant',[variant,*libdirs])]:
        runenv=dict(env,LD_LIBRARY_PATH=':'.join(map(str,paths)))
        ld=subprocess.run(['ldd',str(build/'fixture')],capture_output=True,text=True,env=runenv);(out/(label+'_loaded_libraries.txt')).write_text(ld.stdout)
        selected=[Path(x) for x in ld.stdout.split() if x.startswith('/') and Path(x).is_file()];loaded_hashes[label]={str(x):sha(x) for x in selected}
        assert (str(variant/'libSimTKsimbody.so.3.9') in ld.stdout)==(label=='variant')
        session=Session(build/'fixture',inputs,out/label,runenv);rows=baseline if label=='baseline' else variant_data
        try:
            rows['initial']=session.initial;assert len(rows['initial']['muscles'])==92 and set(rows['initial']['muscles'])=={m['id'] for m in catalog}
            rows['initial_probe']=session.request('mass_probe')
            rows['warm']=session.request('advance .002 0 0');rows['warm_probe']=session.request('mass_probe');session.request('checkpoint warm')
            rows['zero']=session.request('mass_transfer torso 0 0 0 0 0 0 0');rows['after_zero_probe']=session.request('mass_probe')
            assert rows['warm_probe']==rows['after_zero_probe'],'explicit zero changed state/mass/M'
            rows['zero_step']=session.request('advance .002 0 0');rows['zero_step_probe']=session.request('mass_probe')
            if label=='baseline':continue
            session.request('restore warm')
            canonical=json.loads((ROOT/'data/derived/canonical/mechanics.json').read_text());reg=CanonicalRegistration(canonical,rows['initial']);stomach=next(e for e in canonical['entities'] if e['name']=='stomach')
            assert reg.rows[stomach['id']]['body']=='torso'
            station=(np.linalg.inv(np.array(rows['initial']['bodies']['torso']['transform_ground']))@np.linalg.inv(reg.global_map)@np.r_[stomach['centroid_m'],1.])[:3]
            b=rows['warm']['bodies']['torso'];rotation=np.array(b['transform_ground'])[:3,:3];material=np.array(b['origin_velocity_m_s'])+np.cross(b['angular_velocity_rad_s'],rotation@station);velocity=material+np.array([.05,0.,0.])
            site={'canonical_entity_id':stomach['id'],'body':'torso','station_m':station.tolist(),'incoming_velocity_source_m_s':velocity.tolist(),'delta_mass_kg':.5,'registration_basis':reg.rows[stomach['id']],'scope':'Synthetic localized payload at retained stomach centroid; no physiology inventory/action consumed'};dump(out/'site.json',site)
            line='mass_transfer torso .5 '+' '.join(map(str,[*station,*velocity]));rows['transfer']=session.request(line);rows['loaded_probe']=session.request('mass_probe');rows['loaded']=session.request('observe')
            assert rows['loaded_probe']['q']==rows['warm_probe']['q'] and rows['loaded_probe']['z']==rows['warm_probe']['z']
            assert rows['loaded']['metabolic_reference']==rows['warm']['metabolic_reference'] and all(rows['loaded'][k]==rows['warm'][k] for k in ENERGY)
            assert rows['transfer']['capture_dissipation_j']>=-1e-10
            changed=[k for k,v in rows['loaded_probe']['effective_body_mass_kg'].items() if v!=rows['warm_probe']['effective_body_mass_kg'][k]];assert changed==['torso'],changed
            assert np.linalg.norm(np.array(rows['loaded_probe']['M'])-np.array(rows['warm_probe']['M']))>0
            rows['continued']=session.request('advance .002 0 0');rows['continued_probe']=session.request('mass_probe')
            session.request('restore warm');restored=session.request('mass_probe');assert restored==rows['warm_probe'],'warm mass/state restore differs'
            replay_transfer=session.request(line);replay=session.request('advance .002 0 0');replay_probe=session.request('mass_probe')
            assert replay_transfer==rows['transfer'] and replay==rows['continued'] and replay_probe==rows['continued_probe'],'mass/muscle/ledger replay differs'
            try:session.request('mass_transfer torso .6 0 0 0 0 0 0')
            except ValueError as error:assert 'invalid bounded torso' in str(error)
            else:raise AssertionError('out-of-scope mass payload accepted')
            rows['after_rejection_probe']=session.request('mass_probe');assert rows['after_rejection_probe']==replay_probe
            assert rows['continued']['muscle_metabolic_energy_j']>rows['loaded']['muscle_metabolic_energy_j']
            for e in (rows['loaded'],rows['continued']):assert abs(e['muscle_metabolic_energy_j']-e['signed_active_fiber_work_j']-e['muscle_heat_energy_j'])<1e-12
        finally:session.close();dump(out/(label+'_frames.json'),rows)
    differences={key:compare(baseline[key],variant_data[key],key) for key in baseline}
    assert all(sha(path)==digest for hashes in loaded_hashes.values() for path,digest in hashes.items()),'loaded library changed'
    assert all(sha(inputs/n)==digest for n,digest in input_hashes.items())
    assert all(sha(ROOT/path)==digest for path,digest in frozen.items()),'frozen input changed during native fixture'
    assert sha(variant/'libSimTKsimbody.so.3.9')==vm['library_sha256']
    report={'passed':True,'muscles':92,'zero_mass_baseline_variant_max_abs_error':differences,'transfer':variant_data['transfer'],'site':site,'wall_s':time.monotonic()-started,'max_child_rss_kib':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,'loaded_library_sha256':loaded_hashes,'source_library_sha256':frozen,'fixture_source_sha256':sha(source),'variant_library_sha256':vm['library_sha256'],'scope':'Actual92-muscle assembly zero parity, local torso GI-site impulse, constraint/momentum/energy accounting and continuing State/muscle/ledger replay; no physiology inventory wiring or settled support claim'}
    dump(out/'report.json',report);print(json.dumps({'passed':True,'output':str(out),'transfer':report['transfer'],'zero_mass_max_abs_error':max(differences.values()),'wall_s':report['wall_s'],'max_child_rss_kib':report['max_child_rss_kib']},indent=2))
if __name__=='__main__':main()
