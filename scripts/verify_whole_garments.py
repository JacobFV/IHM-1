"""Whole-source garments: topology, force gradients, and actual supported dynamics."""
from pathlib import Path
import hashlib,json,sys,tempfile
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.assembly.whole_garment import materialize_garments,audit_pattern
from ihm.assembly.swept_cloth import step_first_impact

ROOT=Path(__file__).resolve().parents[1]
def main():
    tri=np.array([[0,1,2]])
    a=audit_pattern(np.array([[0.,0,0],[1,0,0],[0,1,0]]),tri)
    assert len(a['boundary_loops'])==1 and a['face_components']==1
    try:audit_pattern(np.array([[0.,0,0],[1,0,0],[0,1,0]]),np.array([[0,1,2],[0,1,2]]))
    except ValueError:pass
    else:raise AssertionError('Duplicate oriented faces accepted')
    garments,identity=materialize_garments(ROOT,areal_density_kg_m2=.18,edge_stiffness_n_m=12.)
    # Missing/mismatched source identities must not silently materialize bodies.
    with tempfile.TemporaryDirectory(prefix='whole-garment-negative-') as temporary:
        negative=Path(temporary)
        for path in ['data/derived/clothing/garments.json','app/src/clothing.js','data/derived/canonical/geometry/body-bp3d-FJ2810.json.gz']:
            target=negative/path;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes((ROOT/path).read_bytes())
        constructor=negative/'app/src/clothing.js';original=constructor.read_bytes();constructor.write_bytes(original+b'\n')
        try:materialize_garments(negative,areal_density_kg_m2=.18,edge_stiffness_n_m=12.)
        except ValueError:pass
        else:raise AssertionError('Changed garment constructor accepted')
        constructor.write_bytes(original);skin=negative/'data/derived/canonical/geometry/body-bp3d-FJ2810.json.gz';skin.write_bytes(skin.read_bytes()+b'\n')
        try:materialize_garments(negative,areal_density_kg_m2=.18,edge_stiffness_n_m=12.)
        except ValueError:pass
        else:raise AssertionError('Changed canonical skin accepted')
    reports={}
    for name,body in garments.items():
        assert len(body.pattern_audit['boundary_loops'])==({'shirt':4,'shorts':3}[name])
        x=body.position_m.copy();x[:,0]*=1.003
        force,energy=body.elastic_forces(x);direction=np.random.default_rng(123).normal(size=x.shape);direction/=np.linalg.norm(direction)
        eps=1e-7;gradient=(body.elastic_forces(x+eps*direction)[1]-body.elastic_forces(x-eps*direction)[1])/(2*eps)
        error=abs(gradient+np.sum(force*direction))
        assert error<1e-7 and np.linalg.norm(force.sum(0))<1e-9 and np.linalg.norm(np.cross(x,force).sum(0))<1e-9
        assert np.isclose(body.mass_kg.sum(),body.pattern_audit['area_m2']*.18,rtol=1e-14)
        cases={}
        for supported in [False,True]:
            bodies,_=materialize_garments(ROOT,areal_density_kg_m2=.18,edge_stiffness_n_m=12.);c=bodies[name]
            if supported:c.fixed_nodes=np.array(max(c.pattern_audit['boundary_loops'],key=lambda v:np.mean(c.position_m[v,1])))
            duration=.02;n=int(np.ceil(duration/min(.00005,c.max_explicit_dt_s)));dt=duration/n
            trajectory=[];reaction=np.zeros(3);defect=0.;max_momentum=0.
            for i in range(n):
                r=step_first_impact({name:c},dt,gravity_m_s2=(0,-9.81,0))
                reaction+=r['support_impulse_ns'];defect+=r['numerical_energy_defect_j'];max_momentum=max(max_momentum,float(np.linalg.norm(r['momentum_residual_ns'])))
                if i==n-1 or i%max(1,n//10)==0:trajectory.append({'time_s':c.time_s,'positions_m':c.position_m.tolist(),'velocities_m_s':c.velocity_m_s.tolist(),'elastic_energy_j':r['elastic_energy_j'],'kinetic_energy_j':r['kinetic_energy_j']})
            assert max_momentum<1e-10
            if supported:assert reaction[1]>0 and r['elastic_energy_j']>0 and np.array_equal(c.position_m[c.fixed_nodes],c.reference_position_m[c.fixed_nodes])
            else:
                assert np.linalg.norm(reaction)==0 and r['elastic_energy_j']<1e-20
                assert np.max(np.abs(c.velocity_m_s[:,1]+9.81*duration))<1e-10
            cases['supported' if supported else 'free']={'dt_s':dt,'trajectory':trajectory,'support_impulse_ns':reaction.tolist(),'numerical_energy_defect_j':defect,'max_momentum_residual_ns':max_momentum}
        reports[name]={'topology':body.pattern_audit,'mass_kg':float(body.mass_kg.sum()),'force_gradient_error_n':error,'cases':cases}
    out=Path(tempfile.mkdtemp(prefix='whole-garments-',dir=ROOT/'data/derived'));(out/'inputs').mkdir()
    paths=['ihm/assembly/whole_garment.py','ihm/assembly/clothing.py','ihm/assembly/swept_cloth.py','ihm/assembly/swept_edge_contact.py','scripts/verify_whole_garments.py','app/src/clothing.js','data/derived/clothing/garments.json','data/derived/canonical/geometry/body-bp3d-FJ2810.json.gz']
    sources={}
    for path in paths:
        p=ROOT/path;target=out/'inputs'/path;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(p.read_bytes());sources[path]=hashlib.sha256(p.read_bytes()).hexdigest()
    result={'passed':True,'identity':identity,'sources_sha256':sources,'garments':reports,'scope':'Whole elastic garment with engineering coefficients and prescribed upper-opening support; no body contact or self-contact, bending, calibrated cloth, or wearable containment claim.'}
    (out/'report.json').write_text(json.dumps(result)+'\n')
    print(json.dumps({'passed':True,'output_dir':str(out),'garments':{k:{'mass_kg':v['mass_kg'],'loops':len(v['topology']['boundary_loops']),'force_gradient_error_n':v['force_gradient_error_n'],'support_impulse_ns':v['cases']['supported']['support_impulse_ns']} for k,v in reports.items()}},indent=2))
if __name__=='__main__':main()
