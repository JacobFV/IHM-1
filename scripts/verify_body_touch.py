"""Independent execution, recorded contact, causal delay and source receipt checks."""
from pathlib import Path
import sys,json,numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))

if __name__=='__main__':
    from ihm.app.experiments import read_experiment
    from ihm.human import ImplicitHuman
    data=read_experiment(ROOT,'forearm-touch');frames=data['frames'];clock=data['clock']
    loaded=data['states']['loaded']
    assert loaded['force_balance_residual_n']<1e-7 and loaded['minimum_jacobian']>.2
    assert loaded['indenter_reaction_n']>0 and loaded['maximum_penetration_m']<1e-10
    reference=np.array(data['states']['reference']['positions_m'])
    released=np.array(data['states']['released']['positions_m'])
    assert np.max(np.linalg.norm(released-reference,axis=1))<1e-7
    delay=data['receptors']['rapid']['delay_s']
    assert all(f['rapid_response']==0 and f['slow_response']==0 for f in frames if f['time_s']<=clock['onset_s']+delay)
    assert max(f['rapid_response'] for f in frames)>0
    assert min(f['rapid_response'] for f in frames if f['time_s']>clock['release_s']+delay)<0
    body=ImplicitHuman.open(ROOT)
    model=body.materialize('ibm-causal',sites_m=[data['anchor']['origin_m']])
    assert model.audit['package_sha256']==data['receptors']['rapid']['package_sha256']
    slow=body.materialize('ibm-causal',response_kind='slow',sites_m=[data['anchor']['origin_m']])
    assert slow.audit['implementation']=='mechanoreceptor_slow'
    transport=body.materialize('body-skin-transport')
    assert transport.to_dict()['native_blood_storage_connected'] is False
    print('Body-attached contact, release, causal signed adaptation/delay, source receipts and implicit materializers PASS')
