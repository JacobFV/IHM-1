"""Read retained native acceptance and independently compare geodesic seeds."""
from pathlib import Path
import csv,json,hashlib
import numpy as np
from ellipsoid_geodesic_reference import solve_path
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/derived/stationary-wrap-candidate-4xo9l1c6'
report=json.loads((OUT/'report.json').read_text());assert report['candidate_gradient_gate_passed'] and not report['whole_source_range_certified']
rows={mode:list(csv.DictReader((OUT/(mode+'.csv')).open())) for mode in ('observed','candidate')}
trace={mode:[json.loads(line) for line in (OUT/(mode+'.jsonl')).read_text().splitlines()] for mode in ('observed','candidate')}
assert len(report['rejected_observations'])==16 and report['accepted_corrected_wrap_observations']==400
base={(r['case'],r['muscle']):r for r in rows['observed']};work=[];power=[];wrench=[]
for row in rows['candidate']:
    original=base[row['case'],row['muscle']]
    if not row['muscle'].startswith(('arm26_BRA_','arm26_BIClong_')):
        assert original['status']==row['status']
        if row['status']=='ok':assert abs(float(original['length'])-float(row['length']))<1e-12
    elif row['status']=='ok':
        values=[float(row[k]) for k in ['fd','moment_arm','lengthening_speed','unit_body_power','unit_resultant_force','unit_resultant_moment']];assert np.isfinite(values).all()
        fd,moment,speed,p,force,torque=values
        assert abs(fd-moment)<1e-6 and abs(speed+.07*fd)<1e-8 and abs(p+speed)<1e-10
        work.append(abs(speed+.07*fd));power.append(abs(p+speed));wrench.append((force,torque))
checks=[]
for side in ('r','l'):
    name=f'local_arm_add_{side}_0:arm26_BIClong_{side}'
    source=next(r for r in trace['observed'] if r['case']==name and r.get('wrapped'))
    actual=next(r for r in trace['candidate'] if r['case']==name and r.get('wrapped'))
    ref=solve_path(source['p1'],source['p2'],source['r1'],source['wrap_length'])
    error=max(np.linalg.norm(ref['r1']-actual['r1']),np.linalg.norm(ref['r2']-actual['r2']),abs(ref['arc_m']-actual['wrap_length']))
    assert error<1e-9
    repeats=[r for r in trace['candidate'] if r.get('wrapped') and r['case'].startswith('local') and r['case'].endswith('_0:arm26_BIClong_'+side)]
    span=max(np.linalg.norm(np.array(r['r1'])-actual['r1']) for r in repeats)
    assert span<1e-9
    checks.append({'side':side,'native_DOP853_reference_error_m':error,'repeated_seed_contact_span_m':span,'repeated_seed_cases':len(repeats)})
print(json.dumps({'passed':True,'native_corrected_observations':len(work),'maximum_speed_vs_length_derivative_error_m_s':max(work),'maximum_unit_tension_power_error_W':max(power),'maximum_unit_resultant_force_N':max(x[0] for x in wrench),'maximum_unit_resultant_moment_Nm':max(x[1] for x in wrench),'independent_geodesic_and_repeat_checks':checks},indent=2))
