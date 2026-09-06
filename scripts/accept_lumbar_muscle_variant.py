"""One bounded native compile/load fixture; run only with coordinated heavy slot."""
from pathlib import Path
import csv,hashlib,json,os,resource,shlex,shutil,subprocess,tempfile,time
import numpy as np
from materialize_lumbar_muscle_variant import materialize,ROOT

def main():
    out=Path(tempfile.mkdtemp(prefix='lumbar-muscle-native-',dir=ROOT/'data/derived'))
    r=materialize(ROOT,out/'variant');runtime=ROOT/'data/runtime/opensim'
    source=ROOT/'scripts/native_lumbar_muscle_probe.cpp';retained=out/source.name;retained.write_bytes(source.read_bytes())
    frozen={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (source,ROOT/'scripts/accept_lumbar_muscle_variant.py',ROOT/'scripts/materialize_lumbar_muscle_variant.py')}
    previous=shlex.split((runtime/'dynamics-adapter/build/CMakeFiles/native_opensim_dynamics.dir/link.txt').read_text())
    libs=[v for v in previous if v.endswith('.so') or '.so.' in v or v.startswith('-Wl,') or v.startswith('-l')]
    for v in libs:
        if Path(v).is_file():frozen[str(Path(v).relative_to(ROOT))]=hashlib.sha256(Path(v).read_bytes()).hexdigest()
    flags=['-std=c++20','-O0','-DSWIG_PYTHON']
    for p in [runtime/'install/opensim/include',runtime/'install/opensim/include/OpenSim',runtime/'install/simbody/include/simbody']:flags+=['-isystem',str(p)]
    exe=out/'probe';command=[shutil.which('prlimit'),'--as=4294967296','--','nice','-n','10','c++',*flags,str(retained),'-o',str(exe),*libs]
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1');started=time.monotonic()
    env['LD_LIBRARY_PATH']=':'.join(str(p) for p in [runtime/'install/opensim/lib',runtime/'install/simbody/lib',runtime/'sysroot/usr/lib/aarch64-linux-gnu/lapack',runtime/'sysroot/usr/lib/aarch64-linux-gnu/blas',runtime/'sysroot/usr/lib/aarch64-linux-gnu'])
    (out/'command.json').write_text(json.dumps(command,indent=2)+'\n')
    print(str(out),flush=True)
    with (out/'compile.log').open('w') as log:subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=120)
    fitted=ROOT/'data/derived/locomotion/moco3d-torque-v1/subject_walk_scaled_FunctionBasedPathSet.xml'
    with (out/'probe.log').open('w') as log:subprocess.run([shutil.which('prlimit'),'--as=2147483648','--',str(exe),str(ROOT/r['model_path']),str(fitted),str(out/'observations.csv')],env=env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=60)
    rows=list(csv.DictReader((out/'observations.csv').open()));assert len(rows)==126
    audit=json.loads((ROOT/'data/research/lumbar_shoulder_coverage/audit.json').read_text());expected={m['candidate_name']:m for m in audit['lumbar_candidate']};coords=audit['registration']['coordinate_order']
    maxfd=maxbalance=0.
    for row in rows:
        for key in ('q','length','moment_arm','finite_difference','tendon_force','active_fiber_force','passive_fiber_force','cos_pennation','activation','fiber_velocity'):row[key]=float(row[key]);assert np.isfinite(row[key])
        maxfd=max(maxfd,abs(row['moment_arm']-row['finite_difference']))
        balance=abs(row['tendon_force']-(row['active_fiber_force']+row['passive_fiber_force'])*row['cos_pennation']);maxbalance=max(maxbalance,balance)
        assert balance<1e-5*max(1.,row['tendon_force']),row
        if row['pose']=='0':
            target=expected[row['muscle']];assert abs(row['length']-target['neutral_path_length_m'])<1e-10
            assert abs(row['moment_arm']-target['neutral_moment_arms_m'][coords.index(row['axis'])])<1e-9
    assert maxfd<1e-8
    assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in {**frozen,**r['sources']}.items())
    report=dict(schema='ihm.lumbar-native-acceptance.v1',passed=True,output_dir=str(out.relative_to(ROOT)),variant_registration=str((out/'variant/registration.json').relative_to(ROOT)),muscles=98,bodies=22,poses=7,observations=len(rows),max_virtual_work_moment_arm_error_m=maxfd,max_tendon_fiber_balance_error_N=maxbalance,activation_values=sorted(set(row['activation'] for row in rows)),native_time_advanced_s=0,default_promotion=False,source_hashes=frozen,executable_sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),wall_s=time.monotonic()-started,child_max_rss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
