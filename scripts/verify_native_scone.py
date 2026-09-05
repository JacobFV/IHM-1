"""Actual native SCONE acceptance; no fallback trajectories."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import re
import tempfile
import numpy as np
ROOT=Path(__file__).resolve().parents[1]


def read_sto(path):
    lines=path.read_text().splitlines(); start=lines.index('endheader')+1
    names=lines[start].split(); data=np.loadtxt(lines[start+1:])
    if data.ndim!=2 or data.shape[1]!=len(names) or not np.isfinite(data).all():raise ValueError('Invalid native trajectory')
    headers=dict(line.split('=',1) for line in lines[:start] if '=' in line)
    if int(headers['nRows'])!=len(data) or int(headers['nColumns'])!=len(names):raise ValueError('Trajectory header mismatch')
    if not np.all(np.diff(data[:,0])>0) or abs(data[0,0])>1e-12:raise ValueError('Invalid native time axis')
    return names,data


def inspect(out):
    files=list(out.glob('trajectory*.sto'))
    if len(files)!=1:raise ValueError('Expected exactly one native trajectory')
    names,data=read_sto(files[0]); execution=json.loads((out/'execution.json').read_text())
    if execution['exit_code']!=0 or abs(data[-1,0]-execution['native_simulation_time_s'])>.011:raise ValueError('Native time receipt mismatch')
    groups={key:[i for i,n in enumerate(names) if n.endswith(suffix)] for key,suffix in [('activation','/activation'),('excitation','.excitation'),('fiber_length','/fiber_length'),('mtu_force','.mtu_force')]}
    for key,idx in groups.items():
        if len(idx)!=14:raise ValueError('Missing fourteen-muscle telemetry: '+key)
        values=data[:,idx]
        if values.min() < -1e-7 or (key in ('activation','excitation') and values.max()>1+1e-7):raise ValueError('Muscle bounds: '+key)
        if key=='fiber_length' and values.min()<=0:raise ValueError('Nonpositive fibers')
    contacts=[i for i,n in enumerate(names) if 'contact' in n.lower() or '.grf' in n.lower()]
    if not contacts:raise ValueError('Missing contact telemetry')
    vertical=[i for i,n in enumerate(names) if n.endswith('.grf_norm_y')]
    if len(vertical)!=2 or data[:,vertical].min() < -1e-7 or data[:,vertical].max()<=0:raise ValueError('Invalid unilateral ground reaction')
    x=names.index('/jointset/ground_pelvis/pelvis_tx/value'); y=names.index('/jointset/ground_pelvis/pelvis_ty/value')
    log=(out/'runner.log').read_text()
    match=re.findall(r'step_count\s*=\s*([0-9]+)',log)
    if not match:raise ValueError('Missing native gait step counter')
    summary={'step_count':int(match[-1]),'duration_s':float(data[-1,0]),'rows':len(data),'columns':len(names),'pelvis_displacement_m':float(data[-1,x]-data[0,x]),'minimum_pelvis_height_m':float(data[:,y].min()),'contact_columns':[names[i] for i in contacts], 'muscle_ranges':{key:[float(data[:,idx].min()),float(data[:,idx].max())] for key,idx in groups.items()},'trajectory_sha256':hashlib.sha256(files[0].read_bytes()).hexdigest()}
    return names,data,summary


def main():
    p=argparse.ArgumentParser();p.add_argument('--output');args=p.parse_args()
    parent=ROOT/'data/derived/locomotion/scone_acceptance';parent.mkdir(parents=True,exist_ok=True)
    out=Path(args.output) if args.output else Path(tempfile.mkdtemp(prefix='run-',dir=parent))
    out.mkdir(parents=True,exist_ok=True)
    results={}; arrays={};checks={}
    invalid=subprocess.run([sys.executable,str(ROOT/'scripts/run_native_scone.py'),'--output',str(out/'invalid'),'--seconds','nan'],capture_output=True)
    checks['nonfinite_duration_rejected']=invalid.returncode!=0 and not (out/'invalid').exists()
    for label,payload in [('nonfinite','time x\n0 1\n1 nan'),('backward','time x\n0 1\n0 2')]:
        bad=out/(label+'.sto');bad.write_text('nRows=2\nnColumns=2\nendheader\n'+payload+'\n')
        try:read_sto(bad)
        except ValueError:checks[label+'_trajectory_rejected']=True
        else:checks[label+'_trajectory_rejected']=False
    for name,seconds,accuracy,extra in [('intact_30s',30,.002,[]),('refined_30s',30,.001,[]),('controller_off',30,.002,['--controller-off','--allow-early-termination']),('coarse_1s',1,.002,[]),('medium_1s',1,.001,[]),('fine_1s',1,.0005,[])]:
        command=[sys.executable,str(ROOT/'scripts/run_native_scone.py'),'--output',str(out/name),'--seconds',str(seconds),'--accuracy',str(accuracy),*extra]
        run=subprocess.run(command,cwd=ROOT,capture_output=True,text=True)
        (out/(name+'-invocation.log')).write_text(run.stdout+run.stderr)
        checks[name+'_process_success']=run.returncode==0
        try: names,data,summary=inspect(out/name);arrays[name]=(names,data);results[name]=summary
        except Exception as exc:results[name]={'error':str(exc)}
    if 'intact_30s' in arrays and 'controller_off' in arrays:
        n,a=arrays['intact_30s'];m,b=arrays['controller_off']
        # All physical state columns must start identically before control acts.
        states=[i for i,k in enumerate(n) if k.startswith('/')]
        other=[m.index(n[i]) for i in states]
        checks['same_initial_physical_state']=np.array_equal(a[0,states],b[0,other])
        pre=b[:,0]<1
        checks['same_preintervention_physical_state']=np.array_equal(a[:sum(pre),states],b[pre][:,other])
        checks['intact_completes_30s']=abs(a[-1,0]-30)<.011
        checks['controller_off_terminates_early']=b[-1,0]<5
        inputs=[i for i,k in enumerate(m) if k.endswith('.input')]
        checks['off_inputs_zero']=len(inputs)==14 and np.max(abs(b[b[:,0]>1.01][:,inputs]))<1e-12
    if all(k in arrays for k in ['coarse_1s','medium_1s','fine_1s']):
        n,a=arrays['coarse_1s'];_,b=arrays['medium_1s'];_,c=arrays['fine_1s']
        # Predeclared engineering acceptance, not a physiological tolerance:
        # first second max coordinate difference <2 cm translation/<0.1 rad angle.
        if a.shape==b.shape==c.shape:
            errors={}
            for label,cols,tolerance in [('translation',[i for i,k in enumerate(n) if k.endswith('/value') and ('pelvis_tx' in k or 'pelvis_ty' in k)],.02),('angle',[i for i,k in enumerate(n) if k.endswith('/value') and not ('pelvis_tx' in k or 'pelvis_ty' in k)],.1)]:
                errors[label]={'coarse_medium_max':float(np.max(abs(a[:,cols]-b[:,cols]))),'medium_fine_max':float(np.max(abs(b[:,cols]-c[:,cols]))),'tolerance':tolerance}
                checks[label+'_refinement_bound']=max(errors[label]['coarse_medium_max'],errors[label]['medium_fine_max'])<tolerance
            results['refinement']=errors
        else:checks['refinement_same_grid']=False
    if 'refined_30s' in arrays and 'intact_30s' in arrays:
        distance_a=results['intact_30s']['pelvis_displacement_m'];distance_b=results['refined_30s']['pelvis_displacement_m']
        relative=abs(distance_a-distance_b)/abs(distance_b)
        results['distance_refinement_relative']=relative
        checks['long_horizon_distance_refinement_5percent']=relative<.05
    negative=subprocess.run([sys.executable,str(ROOT/'scripts/run_native_scone.py'),'--output',str(out/'early_without_opt_in'),'--seconds','3','--controller-off'],capture_output=True,text=True)
    (out/'early_without_opt_in-invocation.log').write_text(negative.stdout+negative.stderr)
    checks['early_termination_requires_opt_in']=negative.returncode!=0 and (out/'early_without_opt_in/execution.json').exists() and not json.loads((out/'early_without_opt_in/execution.json').read_text())['completed_requested_duration']
    checks['all_runs_have_telemetry']=len(arrays)==6
    checks={key:bool(value) for key,value in checks.items()}
    report={'passed':all(checks.values()),'checks':checks,'results':results,'scope':'Planar shipped SCONE subject and optimized reflexes; numerical acceptance, not human validation.'}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(out);print(json.dumps(checks,indent=2))
    if not report['passed']:raise SystemExit(1)

if __name__=='__main__':main()
