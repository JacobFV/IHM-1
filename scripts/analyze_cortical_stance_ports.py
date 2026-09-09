"""Quantify exact static information/actuation losses in cortical port projections."""
from pathlib import Path
import argparse,json,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
import torch
from ihm.native.cortical_stance import load_cortical_stance


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--artifact',required=True);p.add_argument('--teacher',required=True);p.add_argument('--out',required=True);a=p.parse_args()
    torch.set_num_threads(1)
    policy,artifact=load_cortical_stance(a.artifact)
    with np.load(a.teacher,allow_pickle=False) as data:
        K=data['K'];Ad=data['Ad'];Bd=data['Bd'];model=str(data['model_sha256'].item())
    if model!=artifact['provenance']['model_sha256']:raise ValueError('Teacher native model differs')
    E=policy.encoder.detach().numpy().astype(float);scale=policy.state_scale.numpy().astype(float)
    D=policy.decoder.weight.detach().numpy().astype(float);D=D-D.mean(1,keepdims=True)
    Pe=np.linalg.pinv(E)@E;Pd=D@np.linalg.pinv(D)
    Ks=K*scale[None,:];null=Ks@(np.eye(len(scale))-Pe)
    fraction=lambda value:float(np.linalg.norm(value)/np.linalg.norm(Ks))
    def closed(gain):
        vals=np.linalg.eigvals(Ad-Bd@gain)
        return {'spectral_radius':float(abs(vals).max()),'modes_above_1_plus_1e8':int(np.sum(abs(vals)>1+1e-8))}
    variants={'exact_teacher':K,'sensor_projection_only':(Ks@Pe)/scale[None,:],
        'motor_projection_only':Pd@K,'both_port_projections':(Pd@Ks@Pe)/scale[None,:]}
    report={'schema':'ihm.cortical-stance-port-analysis.v1','artifact_sha256':artifact['artifact_sha256'],
        'state_width':len(scale),'sensory_ports':len(E),'sensory_projection_rank':int(np.linalg.matrix_rank(E)),
        'unobserved_instantaneous_state_dimensions':len(scale)-int(np.linalg.matrix_rank(E)),
        'motor_ports':D.shape[1],'motor_readout_rank':int(np.linalg.matrix_rank(D)),
        'teacher_gain_rank':int(np.linalg.matrix_rank(K)),
        'normalized_teacher_gain_fraction_in_sensor_nullspace':fraction(null),
        'normalized_teacher_gain_fraction_outside_motor_space':fraction(Ks-Pd@Ks),
        'normalized_teacher_gain_fraction_lost_by_both':fraction(Ks-Pd@Ks@Pe),
        'no_delay_native_closed_loop':{name:closed(gain) for name,gain in variants.items()},
        'scope':['Frobenius fractions in explicitlydeclared normalizednative state units',
            'Best orthogonal staticgain projections; not an E/I nonlinear simulation',
            'Randomprojection nullspace is instantaneous; persistent dynamical observability could recover some hidden states with an observer',
            'These numbers identify architecture limits and do not establish a standing controller']}
    path=Path(a.out);path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
