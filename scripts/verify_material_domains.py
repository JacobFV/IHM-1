"""Positive source-shaped volumetric material partition and ownership checks."""
from pathlib import Path
import sys,json,tempfile,argparse,hashlib
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))

def verify():
    try:
        from ihm.assembly.material_domains import voxel_partition, MaterialOwnership, source_surface
    except ImportError:raise AssertionError('Heterogeneous material-domain API missing') from None
    from ihm.assembly.mechanics_backend import tetra_box
    # Two deliberately overlapping source boxes must create one owner per cell.
    def box(offset):
        x=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]],float)*.02+offset
        t=np.array([[0,2,1],[0,3,2],[4,5,6],[4,6,7],[0,1,5],[0,5,4],[1,2,6],[1,6,5],[2,3,7],[2,7,6],[3,0,4],[3,4,7]])
        return x,t
    a=box(np.zeros(3));b=box(np.array([.01,0,0]))
    p=voxel_partition([('a',*a),('b',*b)],spacing_m=.005)
    assert len(p['tetrahedra'])==6*96
    assert np.unique(p['cell_indices'],axis=0).shape[0]==len(p['cell_indices'])
    x=p['vertices_m'];t=p['tetrahedra']
    volume=np.linalg.det(np.swapaxes(x[t[:,1:]]-x[t[:,0,None]],1,2))/6
    assert volume.min()>0 and abs(volume.sum()-.03*.02*.02)<1e-15
    assert p['overlapping_source_cells']==32
    owner=MaterialOwnership({'a':.01,'b':.02})
    owner.claim('region',['a','b'])
    try:owner.claim('duplicated',['b'])
    except ValueError:pass
    else:raise AssertionError('Duplicate source material mass accepted')
    assert abs(owner.mass_kg('region')-.03)<1e-15
    print('Material partition, positive elements, overlap and ownership checks passed')

def verify_canonical(directory):
    from ihm.assembly.contact_dynamics import DynamicTetrahedra,step_coupled
    directory=Path(directory);manifest=json.loads((directory/'manifest.json').read_text())
    for name,digest in manifest['artifacts'].items():assert hashlib.sha256((directory/name).read_bytes()).hexdigest()==digest
    data=np.load(directory/'pelvic-domain.npz',allow_pickle=False)
    x=data['vertices_m'];t=data['tetrahedra'];materials=data['material_index'];fixed=data['fixed_nodes']
    assert len(np.unique(materials))==3 and manifest['minimum_reference_tetrahedron_volume_m3']>0
    glans=np.unique(t[materials==2]);owners=tuple(m['source_id'] for m in manifest['material_regions'])
    root=Path(tempfile.mkdtemp(prefix='canonical-material-dynamics-',dir='data/derived'))
    cases={}
    for name,dt,load in [('unloaded',.00005,0.),('loaded',.00005,.02),('refined',.000025,.02)]:
        b=DynamicTetrahedra(x,t,mu_pa=data['mu_pa'],lambda_pa=data['lambda_pa'],density_kg_m3=data['density_kg_m3'],fixed_nodes=fixed,material_owner_ids=owners)
        assert abs(b.mass_kg.sum()-manifest['mass_claim_kg'])<1e-12
        force=np.zeros_like(x);force[glans,1]=load/len(glans)
        snapshots=[x.copy()];times=[0.];audit=[];minj=1.;max_impulse_residual=0.;work=0.;defect=0.
        for i in range(round(.05/dt)):
            r=step_coupled({'pelvic-tissues':b},dt,external_forces_n={'pelvic-tissues':force},gravity_m_s2=(0,0,0))
            minj=min(minj,b.minimum_jacobian());max_impulse_residual=max(max_impulse_residual,np.linalg.norm(r['momentum_residual_ns']))
            work+=r['external_work_j'];defect+=r['numerical_energy_defect_j']
            if (i+1)%round(.005/dt)==0:snapshots.append(b.position_m.copy());times.append(b.time_s)
        assert minj>.9 and max_impulse_residual<1e-12
        if load:assert np.mean(b.position_m[glans,1]-x[glans,1])>.0001
        else:assert np.max(np.abs(b.position_m-x))<1e-12
        case={'minimum_jacobian':minj,'maximum_momentum_residual_ns':max_impulse_residual,
              'glans_mean_displacement_m':np.mean(b.position_m[glans]-x[glans],axis=0).tolist(),
              'maximum_displacement_m':float(np.linalg.norm(b.position_m-x,axis=1).max()),
              'external_work_j':work,'numerical_energy_defect_j':defect,'dt_s':dt,'total_glans_load_n':load}
        np.savez_compressed(root/(name+'.npz'),positions_m=np.array(snapshots),time_s=np.array(times),velocity_m_s=b.velocity_m_s,source_vertices_m=x,tetrahedra=t)
        cases[name]=case
    difference=np.linalg.norm(np.array(cases['loaded']['glans_mean_displacement_m'])-cases['refined']['glans_mean_displacement_m'])
    assert difference<1e-5
    report={'status':'passed','source_domain':str(directory.resolve()),'source_manifest_sha256':hashlib.sha256((directory/'manifest.json').read_bytes()).hexdigest(),
            'cases':cases,'glans_time_refinement_difference_m':float(difference),'scope':'source-shaped inertial deformation under an assumed load and proximal fixture; not calibrated tucking or garment contact'}
    (root/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(report,output_dir=str(root)),indent=2))
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--canonical');a=p.parse_args();verify()
    if a.canonical:verify_canonical(a.canonical)
