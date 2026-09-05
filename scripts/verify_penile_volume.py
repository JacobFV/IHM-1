"""Finite-element work conjugacy and actual source-volume time refinement."""
from pathlib import Path
import sys,json,argparse,tempfile,hashlib
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))

def verify_kernel():
    from ihm.assembly.penile_volume import PenileTetrahedra
    x=np.array([[0.,0,0],[.01,0,0],[0,.01,0],[0,0,.01],[.01,.01,.01]])
    t=np.array([[0,1,2,3],[1,2,3,4]])
    body=PenileTetrahedra(x,t,material_index=np.array([0,1]),
                          material_tissues=('corpus_cavernosum','corpus_spongiosum'),density_kg_m3=1000.)
    y=x@np.array([[1.04,.03,0],[0,.97,.01],[0,0,1.02]]).T
    force,energy=body.elastic_forces(y);numerical=np.zeros_like(y);h=1e-7
    for i,j in np.ndindex(y.shape):
        d=np.zeros_like(y);d[i,j]=h
        numerical[i,j]=-(body.elastic_forces(y+d)[1]-body.elastic_forces(y-d)[1])/(2*h)
    assert np.allclose(force,numerical,rtol=2e-6,atol=1e-9)
    assert np.linalg.norm(force.sum(axis=0))<1e-12
    assert np.linalg.norm(np.cross(y,force).sum(axis=0))<1e-12
    for labels,tissues in (([0,2],('corpus_cavernosum','corpus_spongiosum')),([0,1],('corpus_cavernosum','glans'))):
        try:PenileTetrahedra(x,t,material_index=np.array(labels),material_tissues=tissues,density_kg_m3=1000.)
        except ValueError:pass
        else:raise AssertionError('Unmapped or unsupported tissue accepted')
    print('PASS heterogeneous tetrahedral energy gradient, momentum/torque and explicit tissue scope')
    from ihm.human import ImplicitHuman
    human=ImplicitHuman.open();body=human.materialize('penile-volume')
    assert 'penile-volume' in human.describe()['materializations']
    assert len(body.region.tets)==5172 and body.evidence['canonical_handoff_applied'] is False
    assert body.material_owner_ids==('body-bp3d-FJ3132','body-bp3d-FJ3133')
    assert body.minimum_jacobian()>0.999999999
    print('PASS implicit body materializes the qualified source-shaped CC/CS volume')
    from ihm.assembly.penile_volume import source_cc_cs_domain
    with tempfile.TemporaryDirectory() as tmp:
        root=Path(tmp);parent=root/'domain';parent.mkdir()
        (parent/'manifest.json').write_text(json.dumps({'artifacts':{'../../outside': 'untrusted'}}))
        try:source_cc_cs_domain(parent,evidence_root=root)
        except ValueError as e:assert 'escapes' in str(e)
        else:raise AssertionError('Parent artifact path escaped evidence root')
    print('PASS material evidence path containment')

def verify_experiment(output,*,existing=False):
    from ihm.assembly.penile_volume import run_penile_volume,PenileTetrahedra,TISSUES
    out=Path(output)
    if not existing:
        if out.exists():raise ValueError('Choose a fresh experiment collection')
        out.mkdir(parents=True)
    results={};trajectories={}
    for name,law,dt in [('published','published',.00005),('half_step','published',.000025),
                        ('quarter_step','published',.0000125),('generic_prior','generic_prior',.000025)]:
        results[name]=(json.loads((out/name/'report.json').read_text()) if existing else run_penile_volume(out/name,law=law,dt_s=dt))
        print(name,{k:v for k,v in results[name].items() if k!='artifacts_sha256'},flush=True)
        with np.load(out/name/'trajectory.npz') as data:trajectories[name]={k:data[k] for k in data.files}
    for name,r in results.items():
        assert r['minimum_jacobian']>.8 and r['maximum_momentum_residual_ns']<1e-10
        assert abs(r['energy_audit_closure_j'])<1e-12
        assert abs(r['numerical_energy_defect_j'])<.03*r['external_work_j']
        for filename,digest in r['artifacts_sha256'].items():
            assert hashlib.sha256((out/name/filename).read_bytes()).hexdigest()==digest
        assert np.allclose(trajectories[name]['time_s'],trajectories['published']['time_s'],rtol=0,atol=1e-12)
        config=json.loads((out/name/'configuration.json').read_text())
        reference=json.loads((out/'published/configuration.json').read_text())
        expected_dt={'published':.00005,'half_step':.000025,'quarter_step':.0000125,'generic_prior':.000025}[name]
        assert config['dt_s']==expected_dt and config['law']==('generic_prior' if name=='generic_prior' else 'published')
        for key in ('law','dt_s','initial_tangent_dt_limit_s'):config.pop(key);reference.pop(key)
        assert config==reference,'Compared experiments changed more than law or timestep'
        with np.load(out/name/'domain.npz') as domain,np.load(out/'published/domain.npz') as reference_domain:
            assert domain.files==reference_domain.files and all(np.array_equal(domain[k],reference_domain[k]) for k in domain.files)
    a,b,c=(trajectories[k]['positions_m'] for k in ('published','half_step','quarter_step'))
    coarse=float(np.max(np.linalg.norm(a-b,axis=2)));fine=float(np.max(np.linalg.norm(b-c,axis=2)))
    assert fine<.7*coarse
    assert abs(results['quarter_step']['numerical_energy_defect_j'])<.7*abs(results['half_step']['numerical_energy_defect_j'])
    contrast=float(np.max(np.linalg.norm(b-trajectories['generic_prior']['positions_m'],axis=2)))
    assert contrast>10*fine and contrast>1e-6
    # Actual deformed source-volume virtual work, independent of the stepper.
    with np.load(out/'quarter_step/domain.npz') as data:
        body=PenileTetrahedra(data['vertices_m'],data['tetrahedra'],material_index=data['material_index'],
                              material_tissues=TISSUES,density_kg_m3=data['density_kg_m3'])
    state=c[-1];force,_=body.elastic_forces(state);rng=np.random.default_rng(1402);work_errors=[]
    for _ in range(3):
        direction=rng.normal(size=state.shape);direction/=np.max(np.abs(direction));h=1e-7
        derivative=(body.elastic_forces(state+h*direction)[1]-body.elastic_forces(state-h*direction)[1])/(2*h)
        conjugate=-float(np.sum(force*direction));work_errors.append(abs(derivative-conjugate))
        assert np.isclose(derivative,conjugate,rtol=1e-4,atol=1e-8)
    report={'passed':True,'output_dir':str(out),'maximum_full_trajectory_coarse_half_difference_m':coarse,
            'maximum_full_trajectory_half_quarter_difference_m':fine,
            'maximum_published_generic_trajectory_difference_m':contrast,
            'source_volume_virtual_work_gradient_errors_n':work_errors,
            'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'results':{k:{key:value for key,value in r.items() if key!='artifacts_sha256'} for k,r in results.items()},
            'interpretation':'Numerical time refinement and causal material-law contrast on the same source-derived CC/CS subdomain; no biological validation or spatial convergence claim.'}
    destination=Path(tempfile.mkdtemp(prefix='penile-volume-audit-',dir='artifacts/verification')) if existing else out
    (destination/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(report,verification_path=str(destination/'verification.json')),indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser();group=parser.add_mutually_exclusive_group();group.add_argument('--experiment',type=Path);group.add_argument('--existing',type=Path);args=parser.parse_args();verify_kernel()
    if args.experiment:verify_experiment(args.experiment)
    if args.existing:verify_experiment(args.existing,existing=True)
