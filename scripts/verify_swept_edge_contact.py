"""Analytic swept crossing, impulse, friction and unsupported-case checks."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))

def verify():
    from ihm.assembly.swept_edge_contact import swept_edge_event,resolve_edge_impact
    a=np.array([[-.01,0,.001],[.01,0,.001]]);b=np.array([[0,-.01,0],[0,.01,0]])
    va=np.tile([0.,0,-.1],(2,1));vb=np.zeros((2,3));mass=np.ones(2)*.002
    event=swept_edge_event(a,va,b,vb,.02)
    assert event['status']=='impact' and abs(event['time_s']-.01)<1e-12
    assert np.allclose(event['weights_a'],[.5,.5]) and np.allclose(event['weights_b'],[.5,.5])
    result=resolve_edge_impact(a+va*.01,va,mass,b,vb,mass,event,friction_static=.4,friction_kinetic=.3)
    assert np.allclose(result['velocity_a_m_s'][:,2],-.05) and np.allclose(result['velocity_b_m_s'][:,2],-.05)
    assert np.isclose(result['normal_impulse_ns'],.0002) and np.isclose(result['dissipation_j'],1e-5)
    assert np.linalg.norm(result['momentum_residual_ns'])<1e-14
    assert np.linalg.norm(result['angular_impulse_residual_nms'])<1e-14
    assert abs(result['energy_residual_j'])<1e-14
    for reverse in (False,True):
        aa=a[::-1] if reverse else a;vva=va[::-1] if reverse else va
        e=swept_edge_event(aa,vva,b[::-1],vb,.02);assert e['status']=='impact' and abs(e['time_s']-.01)<1e-12
    e=swept_edge_event(a*np.array([1,1,-1]),-va,b,vb,.02)
    assert e['status']=='impact' and abs(e['time_s']-.01)<1e-12
    for speed,expected in ((.01,'sticking'),(.2,'sliding')):
        velocity=va+np.array([speed,0,0]);start=a-velocity*.0
        e=swept_edge_event(start,velocity,b,vb,.02);r=resolve_edge_impact(start+velocity*e['time_s'],velocity,mass,b,vb,mass,e,friction_static=.4,friction_kinetic=.3)
        assert r['friction_regime']==expected and r['dissipation_j']>=0 and abs(r['energy_residual_j'])<1e-14
        assert np.linalg.norm(r['angular_impulse_residual_nms'])<1e-14
    fixed=resolve_edge_impact(a+va*.01,va,mass,b,vb,mass,event,mobile_b=[False,False],friction_static=0,friction_kinetic=0)
    assert np.allclose(fixed['velocity_a_m_s'],0) and np.linalg.norm(fixed['support_impulse_ns'])>0
    assert np.linalg.norm(fixed['momentum_residual_ns'])<1e-14 and abs(fixed['energy_residual_j'])<1e-14
    moving=np.tile([0.,0,.02],(2,1));e=swept_edge_event(a,va,b,moving,.02);at=e['time_s']
    r=resolve_edge_impact(a+va*at,va,mass,b+moving*at,moving,mass,e,mobile_b=[False,False],friction_static=0,friction_kinetic=0)
    assert r['support_work_j']>0 and abs(r['energy_residual_j'])<1e-14
    endpoint=swept_edge_event(a,va,b+np.array([.01,0,0]),vb,.02)
    assert endpoint['status']=='impact' and np.allclose(endpoint['weights_a'],[0,1])
    changing=va+np.array([[-.01,-.02,0],[.01,.02,0]])
    e=swept_edge_event(a,changing,b,np.array([[.02,-.01,0],[-.02,.01,0]]),.02)
    assert e['status']=='impact' and abs(e['time_s']-.01)<1e-12
    translated=swept_edge_event(a+[1,2,3],va,b+[1,2,3],vb,.02)
    assert translated['status']=='impact' and abs(translated['time_s']-.01)<1e-12
    # Event time is independent of a shorter search interval or high speed.
    assert swept_edge_event(a,va,b,vb,.005)['status']=='no_impact'
    fast=swept_edge_event(a,va*100,b,vb,.02);assert abs(fast['time_s']-.0001)<1e-12
    for n in (1,2,4,8):
        found=[]
        for k in range(n):
            dt=.019/n;offset=k*dt;e=swept_edge_event(a+va*offset,va,b,vb,dt)
            if e['status']=='impact':found.append(offset+e['time_s'])
        assert len(found)==1 and abs(found[0]-.01)<1e-12
    assert swept_edge_event(a+np.array([1,0,0]),va,b,vb,.02)['status']=='no_impact'
    for aa,vv in ((a+va*.01,np.zeros_like(va)),(np.repeat(a[:1],2,axis=0),va)):
        assert swept_edge_event(aa,vv,b,vb,.02)['status']=='unresolved'
    before_a=a.copy();before_v=va.copy()
    for changes in ({'friction_static':-.1,'friction_kinetic':0},{'friction_static':.1,'friction_kinetic':.2},{'friction_static':0,'friction_kinetic':0,'mobile_a':[False,False],'mobile_b':[False,False]}):
        try:resolve_edge_impact(a+va*.01,va,mass,b,vb,mass,event,**changes)
        except ValueError:pass
        else:raise AssertionError('Invalid impact input accepted')
    assert np.array_equal(a,before_a) and np.array_equal(va,before_v)
    out=Path(__file__).resolve().parents[1]/'artifacts/verification/swept-edge-contact';out.mkdir(parents=True,exist_ok=True)
    report={'passed':True,'analytic_event_time_s':event['time_s'],'analytic_normal_impulse_ns':result['normal_impulse_ns'],
            'analytic_dissipation_j':result['dissipation_j'],'analytic_energy_residual_j':result['energy_residual_j'],
            'analytic_momentum_residual_ns':result['momentum_residual_ns'].tolist(),
            'analytic_angular_impulse_residual_nms':result['angular_impulse_residual_nms'].tolist(),
            'source_sha256':{str(p.relative_to(out.parents[2])):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__).resolve(),out.parents[2]/'ihm/assembly/swept_edge_contact.py')},
            'limitations':['Floating-point CCD, not certified exact predicates.','Initial/persistent, coplanar, parallel and grazing degeneracies remain explicit unresolved cases.','No self-contact scheduler, contact thickness or full garment containment.']}
    (out/'kernel-report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS swept/reversed/high-speed/refined events, analytic impulse/loss, friction, fixed support, angular ledger and explicit degeneracy')

if __name__=='__main__':verify()
