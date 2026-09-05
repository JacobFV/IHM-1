"""Numerical calibration, observable and source/holdout checks."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from ihm.calibration import bounded_fit, lateral_field, convert_field


def main():
    x=np.linspace(-2,2,30); y=3*x+7
    f=lambda p:p[0]*x+p[1]
    fit=bounded_fit(f,y,np.ones(30),[1,1],([-10,-10],[10,10]),parameter_scale=[10,10])
    np.testing.assert_allclose(fit['parameters'],[3,7],atol=1e-7)
    assert fit['rank']==2 and fit['identifiable']
    other=bounded_fit(lambda p:f(p)*1000,y*1000,np.ones(30)*1000,[1,1],([-10,-10],[10,10]),parameter_scale=[10,10])
    np.testing.assert_allclose(fit['parameters'],other['parameters'],atol=1e-6)
    np.testing.assert_allclose(fit['covariance'],other['covariance'],rtol=1e-5,atol=1e-8)
    deficient=bounded_fit(lambda p:np.ones(30)*(p[0]+p[1]),np.ones(30)*4,np.ones(30),[1,1],([-10,-10],[10,10]),parameter_scale=[10,10])
    assert deficient['rank']==1 and not deficient['identifiable'] and deficient['covariance'] is None
    np.testing.assert_allclose(convert_field([1,50],'mV/mm','V/m'),[1,50])
    np.testing.assert_allclose(lateral_field([.01,.02],[[0,0,0],[.001,0,0]],[(0,1)]),[-10])
    try: convert_field([1],'mV','V/m')
    except ValueError: pass
    else: raise AssertionError('membrane voltage accepted as field')
    from ihm.calibration.skin import fit_skin_observations, predict_skin
    import json
    rows=[json.loads(line) for line in (Path(__file__).resolve().parents[1]/'data/derived/integumentary/human-wound-observations.jsonl').read_text().splitlines()]
    result=fit_skin_observations(rows)
    assert not set(result['split']['train_subjects']) & set(result['split']['holdout_subjects'])
    assert len(result['observations'])==160
    assert len(result['split']['train_subjects'])+len(result['split']['holdout_subjects'])==40
    assert result['fit']['rank']==6
    changed=[dict(row) for row in rows]
    for row in changed:
        if row['participant_within_group']%5==0 and row['value'] is not None: row['value']+=10000
    alternative=fit_skin_observations(changed)
    np.testing.assert_allclose(result['fit']['parameters'],alternative['fit']['parameters'])
    prediction=predict_skin(result,'18-29','male','LA-1')
    assert np.isfinite(prediction['mean_V_per_m']) and prediction['mean_se_V_per_m']>=0
    from ihm.materialize.skin import SkinPatch, skin_model, electric_field
    from ihm.calibration import skin_surface_field
    patch=SkinPatch.line(); model=skin_model(patch)
    np.testing.assert_allclose(skin_surface_field(model,patch),electric_field(model,patch)['V_per_m'])
    print('calibration: identifiable recovery, unit rescaling, deficient rank and electric observable PASS')

if __name__=='__main__': main()
