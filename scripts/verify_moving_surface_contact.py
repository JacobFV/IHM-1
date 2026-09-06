"""Small independent parity/conservation checks; no native process."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from ihm.assembly.garment_surface_contact import MovingSurfaceContact
from ihm.assembly.contact_dynamics import resolve_node_triangle_contact

def main():
    b=np.array([[0.,0,0],[1,0,0],[0,1,0]]);tri=np.array([[0,1,2]]);a=np.array([[.2,.2,-.001],[.3,.1,-.002],[1.001,0,-.0001]])
    v=np.array([[.2,0,-1.],[2,0,-.1],[0,0,-1.] ]);m=np.array([.002,.003,.004]);motion=np.array([.02,0,0]);dt=.01
    surface=MovingSurfaceContact(b,b+motion*dt,tri,dt,friction_static=.4,friction_kinetic=.3,search_distance_m=.004)
    r=surface.resolve(a,v,m,dt)
    reference=resolve_node_triangle_contact(a,v,m,b,np.tile(motion,(3,1)),np.ones(3),tri,node_ids=np.arange(3),friction_static=.4,friction_kinetic=.3,search_distance_m=.004,mobile_b=np.zeros(3,bool))
    assert np.allclose(r['positions_m'],reference['position_a_m']) and np.allclose(r['velocities_m_s'],reference['velocity_a_m_s'])
    assert np.allclose(r['body_reaction_impulses_ns'],reference['impulse_b_ns']) and r['unresolved_edge_contacts']==1
    delta=.5*np.sum(m[:,None]*(r['velocities_m_s']**2-v**2));loss=r['normal_impact_dissipation_j']+r['friction_dissipation_j']
    assert abs(delta-r['prescribed_surface_work_j']+loss)<1e-14 and np.linalg.norm(r['paired_impulse_residual_ns'])<1e-14
    assert np.linalg.norm(r['paired_contact_point_angular_impulse_residual_nms'])<1e-14
    # Motion-bound query must find a target moved beyond the reference search band.
    moving=MovingSurfaceContact(b,b+[0,0,.1],tri,dt,friction_static=0,friction_kinetic=0,search_distance_m=.004);moving.fraction=1
    q=moving.resolve([[.2,.2,.099]],[[0,0,0]],[.002],dt);assert q['contact_count']==1
    bounded=MovingSurfaceContact(b,b,tri,dt,friction_static=0,friction_kinetic=0,max_candidate_pairs=1)
    try:bounded.resolve(a,v,m,dt)
    except MemoryError:pass
    else:raise AssertionError('Contact query budget ignored')
    print('PASS moving surface parity, equal/opposite impulses, moving-boundary work, edge rejection, conservative motion bound and memory budget')
if __name__=='__main__':main()
