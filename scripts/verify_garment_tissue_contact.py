"""Causal acceptance for a generated-shorts panel / source-tissue experiment."""
from pathlib import Path
import sys,tempfile,json,argparse,hashlib
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))

def verify(seconds=.24,existing=None,refined=None):
    try:
        from ihm.assembly.garment_contact import run_garment_tissue_experiment
    except ImportError:raise AssertionError('Generated garment / source tissue physical experiment missing') from None
    root=Path(existing) if existing else Path(tempfile.mkdtemp(prefix='garment-tissue-',dir='data/derived'));results={}
    for name,friction,contact in [('coupled',True,True),('friction_off',False,True),('contact_off',False,False)]:
        r=json.loads((root/name/'report.json').read_text()) if existing else run_garment_tissue_experiment(root/name,seconds=seconds,friction=friction,contact=contact)
        for filename,digest in r['artifacts_sha256'].items():
            assert hashlib.sha256((root/name/filename).read_bytes()).hexdigest()==digest
        assert r['minimum_jacobian']>.8 and r['maximum_momentum_residual_ns']<1e-10
        assert abs(r['energy_audit_closure_j'])<1e-12
        assert abs(r['numerical_energy_defect_j'])<.02*r['initial_elastic_energy_j']
        raw=np.load(root/name/'trajectory.npz')
        force_residual=np.linalg.norm(raw['tissue_contact_force_n'].sum(axis=1)+raw['panel_contact_force_n'].sum(axis=1),axis=1)
        assert np.max(force_residual)<1e-10
        if contact:
            assert r['contact_resolutions']>0 and r['tissue_contact_impulse_norm_ns']>1e-6
            assert r['maximum_pair_impulse_residual_ns']<1e-10
            assert r['maximum_pair_work_balance_error_j']<1e-12
            assert abs(r['tissue_interface_work_j']+r['panel_interface_work_j']+r['contact_dissipation_j'])<1e-12
        else:assert r['contact_resolutions']==0
        results[name]=r
    coupled=np.array(results['coupled']['final_glans_displacement_m'])
    assert np.linalg.norm(coupled-results['contact_off']['final_glans_displacement_m'])>1e-6
    assert np.linalg.norm(coupled-results['friction_off']['final_glans_displacement_m'])>1e-8
    report={'status':'passed','cases':results,'output_dir':str(root)}
    if refined:
        fine=json.loads((Path(refined)/'report.json').read_text())
        a=np.load(root/'coupled/trajectory.npz');b=np.load(Path(refined)/'trajectory.npz')
        config=json.loads((root/'coupled/configuration.json').read_text());fine_config=json.loads((Path(refined)/'configuration.json').read_text())
        assert config['dt_s']==2*fine_config['dt_s']
        config.pop('dt_s');fine_config.pop('dt_s');assert config==fine_config
        for filename,digest in fine['artifacts_sha256'].items():
            assert hashlib.sha256((Path(refined)/filename).read_bytes()).hexdigest()==digest
        assert np.max(np.abs(a['time_s']-b['time_s']))<1e-10
        diff=np.linalg.norm(coupled-fine['final_glans_displacement_m'])
        assert diff<1e-5 and abs(fine['numerical_energy_defect_j'])<.6*abs(results['coupled']['numerical_energy_defect_j'])
        report['time_refinement']={'path':str(refined),'mean_glans_endpoint_difference_m':float(diff),
            'maximum_tissue_node_endpoint_difference_m':float(np.max(np.linalg.norm(a['tissue_positions_m'][-1]-b['tissue_positions_m'][-1],axis=1))),
            'maximum_panel_node_endpoint_difference_m':float(np.max(np.linalg.norm(a['panel_positions_m'][-1]-b['panel_positions_m'][-1],axis=1))),
            'fine_numerical_energy_defect_j':fine['numerical_energy_defect_j'],
            'fine_contact_dissipation_j':fine['contact_dissipation_j'],
            'interpretation':'Mean tissue response and energy defect improve; local cloth/contact trajectories are not established converged.'}
    # Rechecking retained runs creates a new receipt; original acceptance stays intact.
    destination=Path(tempfile.mkdtemp(prefix='garment-tissue-audit-',dir='data/derived')) if existing else root
    (destination/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(report,verification_path=str(destination/'verification.json')),indent=2))
    return report
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--seconds',type=float,default=.24);p.add_argument('--existing',type=Path);p.add_argument('--refined',type=Path);a=p.parse_args();verify(a.seconds,a.existing,a.refined)
