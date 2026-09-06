"""Retain bounded free-thorax equation evidence; no time integration/native job."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from ihm.assembly.thoracic_mechanism import ThoracicMechanism
from ihm.assembly.thoracic_free_dynamics import ThoracicFreeDynamics

ROOT=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    if args.output.exists():raise ValueError('Require a fresh evidence directory')
    manifest=ROOT/'data/research/thoracic_mechanism/v2/manifest.json'
    paths=[manifest,ROOT/'ihm/assembly/thoracic_free_dynamics.py',ROOT/'ihm/assembly/thoracic_mechanism.py',Path(__file__).resolve()]
    captured={p:p.read_bytes() for p in paths}
    model=ThoracicMechanism(manifest);dynamics=ThoracicFreeDynamics(model)
    q=np.zeros(26);q[0]=.002;q[24]=.0003;q[25]=.001;u=np.linspace(-.03,.04,32)
    for k in model.locked:u[6+k]=0
    result=dynamics.evaluate(q,u);matrix=result['mass_matrix'];accel=result['velocity_derivative']
    eps=1e-5
    mdot=(model.kinetic(q+eps*u[6:])['mass_matrix']-model.kinetic(q-eps*u[6:])['mass_matrix'])/(2*eps)
    momentum=matrix@u;momentum_rate=mdot@u+matrix@accel
    cavity=model.cavity(q,pressure_pa=75);force=np.zeros(32);force[6:]=cavity['pressure_generalized_force']
    pressure=dynamics.solve(result,force)
    energy_error=float(u@matrix@accel+.5*u@mdot@u)
    pressure_error=float(u@matrix@pressure['velocity_derivative']+.5*u@mdot@u-u@force)
    linear=momentum_rate[:3]+np.cross(u[3:6],momentum[:3])
    angular=momentum_rate[3:6]+np.cross(u[3:6],momentum[3:6])+np.cross(u[:3],momentum[:3])
    if abs(energy_error)>2e-10 or abs(pressure_error)>2e-10 or np.max(np.abs(linear))>2e-9 or np.max(np.abs(angular))>2e-9:raise ValueError('Differential conservation check failed')
    if any(p.read_bytes()!=raw for p,raw in captured.items()):raise ValueError('Source changed during probe')
    report={'schema':'thoracic-free-dynamics-probe-v1','native_activation_allowed':False,
        'scope':'Differential source-only RHS evidence at one declared pose; not trajectory convergence or native acceptance',
        'receipts':[{'path':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw)} for p,raw in captured.items()],
        'q':q.tolist(),'velocity':u.tolist(),'velocity_convention':'u=(R^T tdot_world, omega_body, qdot_internal); Rdot=R[omega_body]x',
        'mass_kg':result['mass_kg'],'locked_internal_coordinates':sorted(model.locked),'active_mass_eigenvalues':result['active_eigenvalues'].tolist(),
        'inertial_bias':result['inertial_bias'].tolist(),'free_velocity_derivative':accel.tolist(),
        'mass_directional_difference_step_s':eps,'free_energy_rate_error_W':energy_error,
        'world_linear_momentum_rate_body_components_N':linear.tolist(),'world_angular_momentum_rate_body_components_Nm':angular.tolist(),
        'pressure_pa':75,'pressure_power_W':pressure['external_power_W'],'pressure_energy_rate_error_W':pressure_error,
        'active_equation_residual_max':float(np.max(np.abs(result['equation_residual_active']))),
        'ownership':'Prospective cervical-reduced torso only; no native chest compliance/ideal pressure drive transfer, no new mass/recoil/activation'}
    args.output.mkdir(parents=True);(args.output/'probe.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['mass_kg','free_energy_rate_error_W','pressure_energy_rate_error_W','active_equation_residual_max']}))

if __name__=='__main__':main()
