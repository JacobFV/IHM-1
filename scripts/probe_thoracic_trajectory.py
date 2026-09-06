"""Retain a 50 ms actual-source free trajectory refinement experiment."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from ihm.assembly.thoracic_mechanism import ThoracicMechanism
from ihm.assembly.thoracic_free_dynamics import ThoracicFreeDynamics
from ihm.assembly.thoracic_trajectory import FreeThoraxTrajectory,rotation_increment
ROOT=Path(__file__).resolve().parents[1]


def serial(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,dict):return {k:serial(v) for k,v in value.items()}
    return value


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    if args.output.exists():raise ValueError('Fresh evidence directory required')
    manifest=ROOT/'data/research/thoracic_mechanism/v2/manifest.json'
    paths=[manifest,*[ROOT/'ihm/assembly'/name for name in ['thoracic_mechanism.py','thoracic_free_dynamics.py','thoracic_trajectory.py']],Path(__file__).resolve()]
    captured={p:p.read_bytes() for p in paths}
    model=ThoracicMechanism(manifest);flow=FreeThoraxTrajectory(ThoracicFreeDynamics(model))
    q=np.zeros(26);q[0]=.002;q[24]=.0003;q[25]=.001;u=np.linspace(-.01,.02,32);u[:6]=[.2,.1,-.15,.6,-.4,.5]
    for k in model.locked:u[6+k]=0
    initial={'q':q,'velocity':u,'translation_m':np.array([.3,-.2,.4]),'rotation':rotation_increment([.1,.2,-.1])}
    start=flow.invariants(initial);coarse=flow.step(initial,.05);half=flow.step(initial,.025);fine=flow.step(half,.025)
    endpoints={'coarse':flow.invariants(coarse),'fine':flow.invariants(fine)};errors={}
    for key in ['kinetic_energy_J','world_linear_momentum_kg_m_per_s','world_angular_momentum_kg_m2_per_s']:
        a=float(np.linalg.norm(np.asarray(endpoints['coarse'][key])-start[key]));b=float(np.linalg.norm(np.asarray(endpoints['fine'][key])-start[key]))
        if b>=.4*a+2e-13 or b>=2e-5:raise ValueError('Refinement/conservation tolerance failed')
        errors[key]={'coarse_absolute_error':a,'fine_absolute_error':b,'fine_over_coarse':b/a if a else None}
    if any(p.read_bytes()!=raw for p,raw in captured.items()):raise ValueError('Source changed during probe')
    report={'schema':'thoracic-free-trajectory-probe-v1','native_activation_allowed':False,
        'scope':'One50msfreeinterval withone full vs twohalf Lie midpoint steps; local refinement evidence, no long-duration stability guarantee',
        'receipts':[{'path':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw)} for p,raw in captured.items()],
        'duration_s':.05,'step_sizes_s':[.05,.025],'initial_state':initial,'initial_invariants':start,
        'coarse_endpoint':coarse,'fine_midpoint':half,'fine_endpoint':fine,'endpoint_invariants':endpoints,'errors':errors,
        'locked_internal_coordinates':sorted(model.locked),'no_added_mass_recoil_gravity_or_force':True,
        'physical_locked_joint_reactions_available':False,'native_composition_precondition':'Compose9cervicalbodies7.418919568kg plus this20.235757397kg subsystem against native27.654676965kg torso, including parent cross-inertia'}
    args.output.mkdir(parents=True);(args.output/'probe.json').write_text(json.dumps(serial(report),indent=2)+'\n');print(json.dumps(errors))

if __name__=='__main__':main()
