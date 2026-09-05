"""Independent deformation, work conjugacy and evidence-scope checks."""
import json
from pathlib import Path
import numpy as np
from ihm.calibration.penile import penile_material

def main():
    checks=[]
    rotation=np.array([[0.,-1.,0.],[1.,0.,0.],[0.,0.,1.]])
    F=np.array([[1.13,.08,0.],[.02,.91,.03],[0.,0.,1.04]])
    for tissue in ('corpus_cavernosum','corpus_spongiosum','fascia','tunica_albuginea'):
        law=penile_material(tissue,**({'fiber_direction':[1,0,0]} if tissue=='tunica_albuginea' else {}))
        w,p=law.energy_piola(np.eye(3))
        assert abs(w)<1e-8 and np.max(abs(p))<1e-8
        w,p=law.energy_piola(F)
        wr,pr=law.energy_piola(rotation@F)
        assert np.allclose(w,wr,rtol=1e-9,atol=1e-8)
        assert np.allclose(pr,rotation@p,rtol=1e-9,atol=1e-7)
        numerical=np.zeros((3,3));step=1e-6
        for i,j in np.ndindex(3,3):
            d=np.zeros((3,3));d[i,j]=step
            numerical[i,j]=(law.energy_piola(F+d)[0]-law.energy_piola(F-d)[0])/(2*step)
        assert np.allclose(p,numerical,rtol=1e-5,atol=2e-3),(tissue,p,numerical)
        # Pure dilation has no distortional/fibre strain. This independently
        # tests the paper's different HGO and Ogden volumetric functions.
        a=1.02;J=a**3;_,piola=law.energy_piola(np.eye(3)*a)
        D=law.parameters['D_pa_inverse']
        pressure=(J-1/J)/D if tissue=='tunica_albuginea' else 2*(J-1)/D
        assert np.allclose(piola,np.eye(3)*pressure*J/a,rtol=1e-10,atol=1e-7)
        for invalid in (np.diag([-1.,1,1]),np.zeros((3,3)),np.full((3,3),np.nan)):
            try:law.energy_piola(invalid)
            except ValueError:pass
            else:raise AssertionError('invalid deformation accepted')
        checks.append(tissue+'_energy_stress_objectivity_and_volume')
    # Published Table3 independent initial shear/bulk targets, in Pa.
    for tissue,mu,bulk in [('corpus_cavernosum',22,2/.0018),('corpus_spongiosum',600,20000),('fascia',50,2500)]:
        law=penile_material(tissue);g=1e-5;f=np.eye(3);f[0,1]=g
        assert np.isclose(law.energy_piola(f)[1][0,1]/g,mu,rtol=1e-6)
        assert np.isclose(law.initial_moduli()['bulk_pa'],bulk)
    ta=penile_material('tunica_albuginea',fiber_direction=[1,0,0])
    along=np.diag([1.1,1/np.sqrt(1.1),1/np.sqrt(1.1)])
    across=np.diag([1/np.sqrt(1.1),1.1,1/np.sqrt(1.1)])
    assert ta.energy_piola(along)[0]>ta.energy_piola(across)[0]
    for tissue,kwargs in [('glans',{}),('skin',{}),('tunica_albuginea',{}),('tunica_albuginea',{'fiber_direction':[0,0,0]}),('fascia',{'fiber_direction':[1,0,0]})]:
        try:penile_material(tissue,**kwargs)
        except ValueError:pass
        else:raise AssertionError('unsupported material/fibre mapping accepted')
    checks+=['published_initial_moduli','explicit_anisotropic_fiber','unsupported_anatomy_rejected']
    from ihm.human import ImplicitHuman
    human=ImplicitHuman.open()
    assert 'penile-constitutive' in human.describe()['materializations']
    assert human.materialize('penile-constitutive',tissue='corpus_cavernosum').initial_moduli()['shear_pa']==22
    checks.append('implicit_body_evidence_bound_materialization')
    report={'passed':True,'checks':checks,'scope':'Published ex vivo fitted constitutive response; no new whole-organ or living-human validation'}
    out=Path(__file__).resolve().parents[1]/'artifacts/verification/penile-constitutive';out.mkdir(parents=True,exist_ok=True)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
