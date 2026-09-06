"""Retained sided convergence classification; never relaxes the gradient gate."""
from pathlib import Path
import json,hashlib
ROOT=Path(__file__).resolve().parents[1]

def summarize():
    path=ROOT/'data/derived/native-wrap-work-ovst3fef/report.json';report=json.loads(path.read_text());rows=[]
    for coordinate in dict.fromkeys(r['coordinate'] for r in report['rows']):
        name=('arm26_BRA_' if coordinate.startswith('elbow') else 'arm26_BIClong_')+coordinate[-1]
        sizes=[]
        for h in (1e-3,1e-4,1e-5):
            pair=[next(m for m in r['muscles'] if m['name']==name) for r in report['rows'] if r['coordinate']==coordinate and abs(r['requested_step_rad'])==h]
            errors=[m['trapezoid_virtual_work_error_m_per_rad'] for m in pair]
            sizes.append(dict(step_rad=h,sided_work_errors_m_per_rad=errors,mean_work_error_m_per_rad=sum(errors)/2,sided_spread_m_per_rad=abs(errors[0]-errors[1]),maximum_effective_energy_length_work_error_nm=max(abs(m['effective_energy_length_work_error_nm']) for m in pair)))
        change=abs(sizes[-1]['mean_work_error_m_per_rad']-sizes[-2]['mean_work_error_m_per_rad'])
        rows.append(dict(coordinate=coordinate,muscle=name,steps=sizes,fine_mean_change_m_per_rad=change,
            persistent_above_original_length_gradient_gate=abs(sizes[-1]['mean_work_error_m_per_rad'])>1e-6,
            classification='Persistent fine-step sided force/path discrepancy; coarse sided transitions retained, not assigned a native branch ID'))
    return dict(native_receipt=str(path.relative_to(ROOT)),native_receipt_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),rows=rows,gradient_gate_passed=False,
        next_step='Do not use stored path-length primitive as complete native-force gradient. Either audit native wrap tangent construction independently or formulate a residual trust region using actual native generalized forces; no physical criteria change.')

if __name__=='__main__':print(json.dumps(summarize(),indent=2,allow_nan=False))
