"""Actual elastic cloth-patch first crossing with transactional rejection."""
from pathlib import Path
import hashlib,json,sys,tempfile
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.assembly.clothing import Cloth
from ihm.assembly.swept_cloth import SweptEdgePair,step_first_impact,UnresolvedSweptContact

def patches():
    a=np.array([[-.01,0,.0000099],[.01,0,.0000099],[0,-.02,.0150099]])
    b=np.array([[0,-.01,0],[0,.01,0],[.02,0,-.015]])
    bodies={}
    for name,x in [('a',a),('b',b)]:
        body=Cloth(x,np.array([[0,1,2]]),areal_density_kg_m2=.18,edge_stiffness_n_m=10.)
        body.material_owner_ids=('synthetic-cloth-'+name,);bodies[name]=body
    bodies['a'].velocity_m_s[:]=[.02,0,-.1]
    return bodies

def run_case(dt,contact=True,friction=True):
    bodies=patches();trajectory=[];events=[];loss=0.;defect=0.;momentum=0.;angular=0.
    pair=SweptEdgePair('a',(0,1),'b',(0,1),.4 if friction else 0.,.3 if friction else 0.)
    initial=sum(b.elastic_forces(b.position_m)[1]+.5*np.sum(b.mass_kg[:,None]*b.velocity_m_s**2) for b in bodies.values())
    for i in range(round(.0001/dt)):
        result=step_first_impact(bodies,dt,contacts=[pair] if contact else [])
        loss+=result['contact_dissipation_j'];defect+=result['numerical_energy_defect_j']
        momentum=max(momentum,float(np.linalg.norm(result['momentum_residual_ns'])));angular=max(angular,float(np.linalg.norm(result['angular_impulse_residual_nms'])))
        if result['event']:events.append(result['event'])
        trajectory.append({'time_s':result['time_s'],'a':bodies['a'].position_m.tolist(),'b':bodies['b'].position_m.tolist(),
                           'velocity_a_m_s':bodies['a'].velocity_m_s.tolist(),'velocity_b_m_s':bodies['b'].velocity_m_s.tolist(),
                           'elastic_energy_j':result['elastic_energy_j'],'kinetic_energy_j':result['kinetic_energy_j']})
    return {'dt_s':dt,'events':events,'trajectory':trajectory,'contact_dissipation_j':loss,'numerical_energy_defect_j':defect,
            'maximum_momentum_residual_ns':momentum,'maximum_angular_impulse_residual_nms':angular,
            'energy_closure_j':result['total_energy_j']-initial+loss-defect,
            'final_elastic_energy_j':result['elastic_energy_j'],'masses_kg':{n:b.mass_kg.tolist() for n,b in bodies.items()}}

def verify():
    # Unsupported geometry after force prediction must leave every owner untouched.
    bodies=patches();before={n:(b.position_m.copy(),b.velocity_m_s.copy(),b.time_s) for n,b in bodies.items()}
    duplicate=[SweptEdgePair('a',(0,1),'b',(0,1))]*2
    try:step_first_impact(bodies,.0001,contacts=duplicate)
    except UnresolvedSweptContact:pass
    else:raise AssertionError('Unscheduled simultaneous impact accepted')
    for n,b in bodies.items():assert np.array_equal(b.position_m,before[n][0]) and np.array_equal(b.velocity_m_s,before[n][1]) and b.time_s==before[n][2]
    # A later owner rejecting final geometry must also leave earlier owners intact.
    original=bodies['b'].elastic_forces;calls=[0]
    def reject_final(x):
        calls[0]+=1
        if calls[0]==2:return np.full_like(x,np.nan),float('nan')
        return original(x)
    bodies['b'].elastic_forces=reject_final
    try:step_first_impact(bodies,.0001,contacts=[])
    except ValueError:pass
    else:raise AssertionError('Invalid final elastic state accepted')
    for n,b in bodies.items():assert np.array_equal(b.position_m,before[n][0]) and np.array_equal(b.velocity_m_s,before[n][1]) and b.time_s==before[n][2]
    results={name:run_case(dt,contact,friction) for name,dt,contact,friction in [
        ('coarse',.0001,True,True),('half',.00005,True,True),('quarter',.000025,True,True),
        ('contact_off',.000025,False,False),('friction_off',.000025,True,False)]}
    for name,r in results.items():
        assert r['maximum_momentum_residual_ns']<1e-12 and r['maximum_angular_impulse_residual_nms']<1e-12
        assert abs(r['energy_closure_j'])<1e-15
        if name!='contact_off':assert len(r['events'])==1 and r['contact_dissipation_j']>0 and r['final_elastic_energy_j']>0
        else:assert not r['events']
    positions=lambda name:np.array(results[name]['trajectory'][-1]['a'])
    assert np.max(np.abs(positions('quarter')-positions('contact_off')))>1e-8
    assert np.max(np.abs(positions('quarter')-positions('friction_off')))>1e-9
    for name in ('coarse','half'):
        assert abs(results[name]['events'][0]['body_time_s']-results['quarter']['events'][0]['body_time_s'])<1e-12
        assert np.max(np.abs(positions(name)-positions('quarter')))<1e-12
    root=Path(__file__).resolve().parents[1];destination=Path(tempfile.mkdtemp(prefix='swept-cloth-',dir=root/'data/derived'))
    sources=[Path(__file__).resolve(),root/'ihm/assembly/clothing.py',root/'ihm/assembly/swept_cloth.py',root/'ihm/assembly/swept_edge_contact.py']
    (destination/'inputs').mkdir()
    for p in sources:(destination/'inputs'/p.name).write_bytes(p.read_bytes())
    report={'passed':True,'results':results,'source_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
            'artifacts_sha256':{str(p.relative_to(destination)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((destination/'inputs').iterdir())},
            'scope':'Actual cloth masses and spring energies, one isolated crossing per run; refinement verifies event/first-impact parity, not sustained-contact relaxation or full garment containment.',
            'assumptions':{'areal_density_kg_m2':.18,'edge_stiffness_n_m':10.,'friction_static':.4,'friction_kinetic':.3,'evidence':'explicit engineering fixture priors'}}
    (destination/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'passed':True,'output_dir':str(destination),'event_time_s':results['quarter']['events'][0]['body_time_s'],'contact_loss_j':results['quarter']['contact_dissipation_j'],'final_elastic_energy_j':results['quarter']['final_elastic_energy_j'],'atomic_failure_test':True},indent=2))

if __name__=='__main__':verify()
