"""Retain bounded static/kinematic mechanism evidence; no native launch."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from ihm.assembly.thoracic_mechanism import ThoracicMechanism

ROOT=Path(__file__).resolve().parents[1]

def receipt(path):
    raw=path.read_bytes();return {'path':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw)}

def probe(output):
    output=Path(output).resolve()
    if output.exists() or not output.is_relative_to(ROOT):raise ValueError('Fresh workspace output required')
    path=ROOT/'data/research/thoracic_mechanism/v2/manifest.json';model=ThoracicMechanism(path)
    reference=model.kinetic(np.zeros(26));q=np.zeros(26);q[0]=.003;q[24]=.0005;q[25]=.002
    velocity=np.linspace(-.03,.04,32)
    for k in model.locked:velocity[6+k]=0
    moved=model.kinetic(q,velocity);cavity=model.cavity(q,pressure_pa=75.)
    fd=np.zeros(26);h=1e-7
    for k in range(26):
        if k in model.locked:continue
        step=np.zeros(26);step[k]=h
        fd[k]=(model.cavity(q+step)['geometric_volume_m3']-model.cavity(q-step)['geometric_volume_m3'])/(2*h)
    nodal_velocity=np.einsum('ijk,k->ij',cavity['material_jacobian'],velocity[6:])
    nodal_power=float(np.sum(cavity['pressure_nodal_forces']*nodal_velocity))
    expected_power=75*float(fd@velocity[6:]);orientation=[]
    for ident,m in model.material.items():
        x,_=model.material_state(ident,q);f=m['faces'];a=m['reference'][f];b=x[f]
        normal0=np.cross(a[:,1]-a[:,0],a[:,2]-a[:,0]);normal=np.cross(b[:,1]-b[:,0],b[:,2]-b[:,0])
        nondegenerate=np.linalg.norm(normal0,axis=1)>1e-20
        reversed_count=int(np.sum((np.einsum('ij,ij->i',normal0,normal)<=0)&nondegenerate))
        if reversed_count:orientation.append({'entity':ident,'reversed_triangles':reversed_count})
    native_root=ROOT/'data/raw/physiology/biogears/projects/biogears/libBiogears/src/engine'
    ownership=[]
    for file in [native_root/'Controller/BioGears.cpp',native_root/'Systems/Respiratory.cpp']:
        lines=file.read_text().splitlines();matches=[]
        for index,line in enumerate(lines,1):
            if any(word in line for word in ['LeftPleuralCavityToRespiratoryMuscle','RightPleuralCavityToRespiratoryMuscle','EnvironmentToRespiratoryMuscle']):matches.append({'line':index,'text':line.strip()})
        ownership.append(receipt(file)|{'path_matches':matches})
    prior=json.loads((ROOT/'data/research/cervical_inertia/v2/manifest.json').read_bytes())['residual_torso']
    result={'schema':'ihm.thoracic-mechanism-probe.v1','native_launched':False,'native_activation_allowed':False,
        'inputs':[receipt(path),receipt(Path(__file__).resolve()),receipt(ROOT/'ihm/assembly/thoracic_mechanism.py')],
        'reference':{'mass_kg':reference['mass_kg'],'center_m':reference['center_m'].tolist(),'inertia_kg_m2':reference['inertia_kg_m2'].tolist(),
            'mass_error_kg':reference['mass_kg']-prior['mass_kg'],'com_error_m':float(np.linalg.norm(reference['center_m']-prior['center_m'])),
            'inertia_error_kg_m2':float(np.linalg.norm(reference['inertia_kg_m2']-prior['inertia_kg_m2'])),
            'active_mass_min_eigenvalue':float(reference['active_eigenvalues'][0]),'eigenvalue_scope':'Declared mixed translational/rotational coordinates; positivity test, not unit-invariant conditioning'},
        'perturbed_q':q.tolist(),'test_velocity':velocity.tolist(),'active_generalized_indices':moved['active_indices'],
        'generalized_mass_matrix':moved['mass_matrix'].tolist(),'active_mass_eigenvalues':moved['active_eigenvalues'].tolist(),
        'kinetic_energy_J':moved['kinetic_energy_J'],'direct_material_energy_J':moved['direct_material_energy_J'],
        'kinetic_identity_error_J':moved['kinetic_energy_J']-moved['direct_material_energy_J'],
        'material_orientation_reversals':orientation,
        'cavity':{'reference_m3':model.reference_volume,'perturbed_m3':cavity['geometric_volume_m3'],
            'analytic_J':cavity['volume_jacobian_m3_per_coordinate'].tolist(),'finite_difference_J':fd.tolist(),
            'max_J_difference':float(np.max(np.abs(fd-cavity['volume_jacobian_m3_per_coordinate']))),
            'test_pressure_pa':75.,'pressure_sign':'Positive outward transmural mechanical pressure; no native driver sign conversion is assumed',
            'nodal_pressure_power_W':nodal_power,'finite_difference_pressure_power_W':expected_power,
            'pressure_work_error_W':nodal_power-expected_power,
            'pressure_resultant_force_N':cavity['pressure_resultant_force_N'].tolist(),'pressure_resultant_moment_Nm':cavity['pressure_resultant_moment_Nm'].tolist(),
            'gas_reference_volume_supplied':False},
        'native_ownership_receipts':ownership,
        'native_ownership_decision':'Keep bilateral chest/lung recoil and ideal native pressure drive. NoK/damping/activation is introduced; any later mechanical recoil/drive needs explicit replacement and clock/PV acceptance.',
        'acceptance_scope':'Source kinematics, material kinetic metric, reference mass ledger and virtual work only; no dynamic integration or physiological calibration'}
    output.mkdir(parents=True);(output/'probe.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return {'path':str(output),'mass_error_kg':result['reference']['mass_error_kg'],'kinetic_error_J':result['kinetic_identity_error_J'],
        'cavity_J_error':result['cavity']['max_J_difference'],'pressure_work_error_W':result['cavity']['pressure_work_error_W'],'orientation_reversals':orientation}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);print(json.dumps(probe(p.parse_args().output),indent=2))
