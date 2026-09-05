"""Independent objective elastic attachment and actual atlas pairing checks."""
from pathlib import Path
import json,sys,tempfile,hashlib
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.assembly.garment_body_ports import BarycentricAttachment,attachment_forces,build_source_pairings,SourceSurface
from ihm.assembly.contact_dynamics import resolve_node_triangle_contact

def main():
    a=np.array([[.2,.3,.1]]);b=np.array([[0.,0,0],[1,0,0],[0,1,0]])
    port=BarycentricAttachment(0,(0,1,2),(.5,.2,.3),rest_length_m=.07,stiffness_n_m=12.)
    r=attachment_forces(a,b,[port]);assert np.allclose(r['garment_force_n'],[[0,0,-.36]])
    assert np.linalg.norm(r['paired_force_residual_n'])<1e-14 and np.linalg.norm(r['paired_torque_residual_nm'])<1e-14
    eps=1e-6
    for who,x in [('garment',a),('body',b)]:
        for i in range(len(x)):
            for j in range(3):
                plus=x.copy();minus=x.copy();plus[i,j]+=eps;minus[i,j]-=eps
                if who=='garment':up=attachment_forces(plus,b,[port])['energy_j'];um=attachment_forces(minus,b,[port])['energy_j']
                else:up=attachment_forces(a,plus,[port])['energy_j'];um=attachment_forces(a,minus,[port])['energy_j']
                assert abs((up-um)/(2*eps)+r[who+'_force_n'][i,j])<1e-8
    rot=np.array([[0,0,1],[0,1,0],[-1,0,0]]);shift=np.array([2,-3,.4]);rr=attachment_forces(a@rot.T+shift,b@rot.T+shift,[port])
    assert abs(rr['energy_j']-r['energy_j'])<1e-14 and np.allclose(rr['garment_force_n'],r['garment_force_n']@rot.T)
    for bad in [(-.1,.5,.6),(.1,.1,.1)]:
        try:BarycentricAttachment(0,(0,1,2),bad,rest_length_m=0,stiffness_n_m=12)
        except ValueError:pass
        else:raise AssertionError('Invalid barycentric weights accepted')
    surface=SourceSurface(b,[[0,1,2]])
    for point,expected,feature in [([.2,.3,.1],[.2,.3,0],'face'),([.4,-.2,0],[.4,0,0],'edge'),([-1,-1,0],[0,0,0],'vertex')]:
        q=surface.closest(point);assert np.allclose(q['body_point_m'],expected) and q['closest_feature']==feature
    # The prescribed-body limit of the existing finite-mass contact kernel is
    # exactly independent of its required but inactive target mass placeholders.
    contacts=[]
    for placeholder in [1e-9,1.,1e9]:
        contact=resolve_node_triangle_contact([[.2,.3,-.001]],[[.2,0,-.1]],[.001],b,np.zeros_like(b),np.full(3,placeholder),np.array([[0,1,2]]),
                    node_ids=np.array([0]),friction_static=.4,friction_kinetic=.3,search_distance_m=.01,mobile_b=np.zeros(3,bool))
        assert contact['contact_count']==1 and np.linalg.norm(contact['paired_impulse_residual_ns'])<1e-14
        assert np.array_equal(contact['position_b_m'],b) and np.array_equal(contact['velocity_b_m_s'],np.zeros_like(b))
        assert abs(contact['kinetic_transfer_a_j']+contact['dissipation_j'])<1e-14
        contacts.append(contact)
    for contact in contacts[1:]:
        assert np.array_equal(contact['velocity_a_m_s'],contacts[0]['velocity_a_m_s']) and np.array_equal(contact['impulse_b_ns'],contacts[0]['impulse_b_ns'])
    # Retained source pairing is geometric: edges/vertices are allowed as
    # attachment endpoints, without claiming their frictional contact solved.
    root=Path(__file__).resolve().parents[1];pairings=build_source_pairings(root)
    assert set(pairings['garments'])=={'shirt','shorts'}
    for name,g in pairings['garments'].items():
        assert g['support_pairings'] and g['contact_samples']
        for p in g['support_pairings']+g['contact_samples']:
            assert abs(sum(p['barycentric'])-1)<1e-12 and min(p['barycentric'])>=-1e-12
            assert p['reconstruction_error_m']<1e-12
    out=Path(tempfile.mkdtemp(prefix='garment-body-ports-',dir=root/'data/derived'));(out/'pairings.json').write_text(json.dumps(pairings,indent=2)+'\n')
    sources=['ihm/assembly/garment_body_ports.py','scripts/verify_garment_body_ports.py','ihm/assembly/whole_garment.py','ihm/assembly/clothing.py','ihm/assembly/contact_dynamics.py','app/src/clothing.js','data/derived/clothing/garments.json','data/derived/canonical/geometry/body-bp3d-FJ2810.json.gz']
    for source in sources:
        p=out/'inputs'/source;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes((root/source).read_bytes())
    report={'passed':True,'sources_sha256':{s:hashlib.sha256((root/s).read_bytes()).hexdigest() for s in sources},'pairings_sha256':hashlib.sha256((out/'pairings.json').read_bytes()).hexdigest(),'sample_energy_j':r['energy_j']}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'passed':True,'output_dir':str(out),'counts':{k:{'supports':len(v['support_pairings']),'contact_samples':len(v['contact_samples'])} for k,v in pairings['garments'].items()}},indent=2))
if __name__=='__main__':main()
