"""Continuum inertia and two-way finite-mass contact behavioral checks."""
from pathlib import Path
import sys,json,tempfile
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))

def verify():
    try:
        from ihm.assembly.contact_dynamics import DynamicTetrahedra, NodeTriangleContact, step_coupled, resolve_node_triangle_contact
    except ImportError:raise AssertionError('Dynamic coupled continuum contact API missing') from None
    from ihm.assembly.mechanics_backend import tetra_box
    # One finite-mass point strikes a finite-mass triangle: both move, the normal
    # impact loses the analytically expected relative kinetic energy.
    xa=np.array([[.2,.2,.0001]]);xb=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]])
    va=np.array([[.3,0,-1.]]);vb=np.zeros((3,3));ma=np.array([2.]);mb=np.ones(3)*3
    result=resolve_node_triangle_contact(xa+va*.001,va,ma,xb,vb,mb,np.array([[0,1,2]]),node_ids=np.array([0]),friction_static=.5,friction_kinetic=.3,search_distance_m=.01)
    pa=(result['velocity_a_m_s']-va)*ma[:,None];pb=(result['velocity_b_m_s']-vb)*mb[:,None]
    assert np.linalg.norm(pa.sum(axis=0)+pb.sum(axis=0))<1e-12
    assert np.linalg.norm(pb)>0 and result['contact_count']==1
    before=.5*np.sum(ma[:,None]*va**2)+.5*np.sum(mb[:,None]*vb**2)
    after=.5*np.sum(ma[:,None]*result['velocity_a_m_s']**2)+.5*np.sum(mb[:,None]*result['velocity_b_m_s']**2)
    assert abs(before-after-result['dissipation_j'])<1e-12
    assert 'kinetic_transfer_a_j' in result, 'Interface impulse work receipt missing'
    assert abs(result['kinetic_transfer_a_j']+result['kinetic_transfer_b_j']+result['dissipation_j'])<1e-12
    # Homogeneous strain uses the independently benchmarked native FEBio law.
    x,t=tetra_box((.02,.02,.02),(2,2,2))
    a=DynamicTetrahedra(x,t,mu_pa=500,lambda_pa=0,density_kg_m3=1000)
    y=x.copy();y[:,2]*=.9
    force,energy=a.elastic_forces(y)
    top=np.isclose(x[:,2],.02)
    assert abs(-force[top,2].sum()/(.02*.02)-500*(.9-1/.9))<1e-8
    native_reference=Path('data/derived/mechanics-reference/compression/benchmark.json')
    if native_reference.exists():
        native=json.loads(native_reference.read_text())
        assert native['status']=='passed' and native['returncode']==0
        assert np.max(np.abs(np.asarray(native['final_stress_Pa'])[:,-1]-(-force[top,2].sum()/(.02*.02))))<1e-7
    # Common rigid motion produces no strain energy and no artificial impulse.
    a.velocity_m_s[:]=(.1,.2,.3)
    r=step_coupled({'a':a},min(a.max_explicit_dt_s,1e-5),gravity_m_s2=(0,0,0))
    assert np.linalg.norm(r['momentum_residual_ns'])<1e-12
    assert abs(r['elastic_energy_j'])<1e-12
    assert abs(r['numerical_energy_defect_j'])<1e-12
    # A duplicate material identity must be rejected even if objects differ.
    a.material_owner_ids=('same-material',)
    duplicate=DynamicTetrahedra(x,t,mu_pa=500,lambda_pa=0,density_kg_m3=1000)
    duplicate.material_owner_ids=('same-material',)
    duplicate.time_s=a.time_s
    try:step_coupled({'a':a,'duplicate':duplicate},1e-5,gravity_m_s2=(0,0,0))
    except ValueError:pass
    else:raise AssertionError('Duplicate physical mass owners accepted')
    # Deformation must oscillate under continuum inertia with convergent energy
    # defect; density and spatial stiffness affect the computed response.
    x,t=tetra_box((.02,.02,.04),(2,2,4));fixed=np.flatnonzero(x[:,2]==0)
    def beam(dt,stiffness=1.,density=1000):
        centers=x[t].mean(axis=1);mu=np.where(centers[:,2]>.02,500*stiffness,500.)
        body=DynamicTetrahedra(x,t,mu_pa=mu,lambda_pa=0,density_kg_m3=density,fixed_nodes=fixed)
        body.position_m[:,0]+=.001*(x[:,2]/.04)**2
        initial=body.elastic_forces(body.position_m)[1];max_residual=0.;defect=0.
        for _ in range(round(.02/dt)):
            r=step_coupled({'beam':body},dt,gravity_m_s2=(0,0,0))
            max_residual=max(max_residual,np.linalg.norm(r['momentum_residual_ns']));defect+=r['numerical_energy_defect_j']
        assert body.minimum_jacobian()>.9 and max_residual<1e-10
        return body.position_m.copy(),{'energy_defect_j':defect,'initial_energy_j':initial,'final_tip_x_m':float(body.position_m[x[:,2]==.04,0].mean()),'momentum_residual_ns':max_residual}
    coarse,mc=beam(.0001);medium,mm=beam(.00005);fine,mf=beam(.000025)
    stiff,_=beam(.00005,stiffness=2);heavy,_=beam(.00005,density=2000)
    assert np.linalg.norm(fine-medium)<np.linalg.norm(medium-coarse)
    assert abs(mf['energy_defect_j'])<abs(mc['energy_defect_j'])
    assert np.linalg.norm(stiff-medium)>1e-6 and np.linalg.norm(heavy-medium)>1e-6
    # A failed trial must leave every body unchanged.
    saved=(a.position_m.copy(),a.velocity_m_s.copy(),a.time_s)
    bad=np.zeros_like(a.position_m);bad[-1]=float('nan')
    try:step_coupled({'a':a},1e-5,external_forces_n={'a':bad})
    except ValueError:pass
    else:raise AssertionError('Nonfinite trial accepted')
    assert np.array_equal(a.position_m,saved[0]) and np.array_equal(a.velocity_m_s,saved[1]) and a.time_s==saved[2]
    report={'status':'passed','beam_coarse':mc,'beam_medium':mm,'beam_fine':mf,
            'last_refinement_position_difference_m':float(np.linalg.norm(fine-medium)),
            'stiffness_position_difference_m':float(np.linalg.norm(stiff-medium)),
            'density_position_difference_m':float(np.linalg.norm(heavy-medium))}
    root=Path(tempfile.mkdtemp(prefix='contact-dynamics-',dir='data/derived'));(root/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(report,output_dir=str(root)),indent=2))
if __name__=='__main__':verify()
