"""Phenotypic human wound-field fit with subject-cluster uncertainty."""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.stats import t as student_t
from .fit import bounded_fit
from .observables import convert_field

FEATURES=['intercept','older_65_80','female','site_LA_1','site_LA_2','site_LL_1']
LIMITATIONS=[
    'Observation is lateral field above epidermis below stratum corneum; not transmembrane voltage, conductance or current.',
    'Additive age-group, sex and site associations are observational, not causal; omitted interactions and between-person heterogeneity remain.',
    'Per-person SEM describes usually three repeated scans, exact count unknown; it is retained but never treated as population SD or inverse-variance population precision.',
    'Subject keys combine source table age group, sex and participant row; repeated sites stay in one split. Cross-site row identity follows the supplied table extraction, not an independently released subject identifier.',
    'Young age category is 18-29 in table but 18-25 in methods/abstract; discrepancy retained and unresolved.',
    'Cluster sandwich uncertainty uses only training subjects and a small-sample correction; intervals are approximate with 32 training clusters.',
    'No heldout tuning, external cohort validation, microscopic ion channel identification or wound-healing outcome prediction.',
    'Field gradient observation operator binds apical extracellular potential to V/m, but human group fits do not calibrate the illustrative skin RC conductances.'
]


def _features(age,sex,site):
    if age not in ('18-29','65-80') or sex not in ('male','female') or site not in ('control','LA-1','LA-2','LL-1'):
        raise ValueError('prediction must use observed age groups, sex categories and sites')
    return np.array([1,age=='65-80',sex=='female',site=='LA-1',site=='LA-2',site=='LL-1'],float)


def _subject(row):
    return f"{row['age_group_table']}:{row['sex']}:{row['participant_within_group']}"


def predict_skin(result,age_group,sex,site):
    x=_features(age_group,sex,site); b=np.array(result['fit']['parameters']); cov=np.array(result['uncertainty']['cluster_covariance'])
    mean=float(x@b); se=float(np.sqrt(max(0,x@cov@x)))
    q=result['uncertainty']['t_critical_95']; spread=result['uncertainty']['training_residual_sd_V_per_m']
    half=q*np.sqrt(se**2+spread**2)
    return dict(mean_V_per_m=mean,mean_se_V_per_m=se,mean_interval_95_V_per_m=[mean-q*se,mean+q*se],
                individual_interval_95_V_per_m=[float(mean-half),float(mean+half)],
                interval_note='approximate pooled residual individual interval; heterogeneity may vary by group/site')


def fit_skin_observations(rows):
    if not rows: raise ValueError('no observations')
    observations=[]
    for row in rows:
        item=dict(row); item['subject_id']=_subject(row)
        item['split']='holdout' if int(row['participant_within_group'])%5==0 else 'train'
        item['included']=not row['missing'] and not row['excluded_by_authors']
        item['exclusion_reason']='missing' if row['missing'] else ('excluded_by_authors' if row['excluded_by_authors'] else None)
        if item['included']:
            if row.get('measurement')!='lateral field above epidermis, below stratum corneum':
                raise ValueError('observation does not match lateral field operator')
            item['field_V_per_m']=float(convert_field(row['value'],row['units']))
        else: item['field_V_per_m']=None
        observations.append(item)
    train=[r for r in observations if r['included'] and r['split']=='train']
    hold=[r for r in observations if r['included'] and r['split']=='holdout']
    if not train or not hold: raise ValueError('both training and holdout observations required')
    design=lambda data:np.array([_features(r['age_group_table'],r['sex'],r['site']) for r in data])
    x=design(train); y=np.array([r['field_V_per_m'] for r in train]); p=x.shape[1]
    # Equal row weights: scan SEM is not population heterogeneity. Unit scaling
    # below is solely numerical and confers no observational precision.
    numerical_scale=max(float(y.std()),1.)
    fit=bounded_fit(lambda b:x@b,y,np.full(len(y),numerical_scale),np.zeros(p),
        (np.full(p,-1000.),np.full(p,1000.)),parameter_scale=np.full(p,1000.))
    if not fit['identifiable'] or fit['active_bounds']: raise ValueError('phenotypic design deficient or bound constrained')
    fit['covariance']=None
    fit['covariance_assumption']='replaced by subject-cluster covariance; equal numerical weights do not encode known measurement precision'
    b=np.array(fit['parameters']); residual=y-x@b
    subjects=sorted({_subject(r) for r in train}); g=len(subjects); n=len(y)
    if g<=p or n<=p: raise ValueError('insufficient subject clusters for uncertainty')
    bread=np.linalg.inv(x.T@x); meat=np.zeros((p,p))
    for subject in subjects:
        mask=np.array([_subject(r)==subject for r in train]); score=x[mask].T@residual[mask]
        meat+=np.outer(score,score)
    cov=bread@meat@bread*(g/(g-1))*((n-1)/(n-p))
    result=dict(schema_version=1,id='human_wound_field_phenotype',kind='observational_additive_phenotype',
        quantity='lateral_field_above_epidermis',unit='V/m',feature_names=FEATURES,fit=fit,
        model_equation='intercept + older + female + site indicators; output V/m',
        parameter_bounds_note='broad numerical bounds +/-1000 V/m, not measured physiological priors',
        weighting='equal observations; supplied scan SEM retained as metadata only',
        split=dict(protocol='deterministic participant numbers divisible by five held out in every age/sex group',
            train_subjects=subjects,holdout_subjects=sorted({_subject(r) for r in hold}),
            train_observations=n,holdout_observations=len(hold)),
        uncertainty=dict(method='subject-cluster sandwich CR1, Student t G-1',cluster_count=g,
            cluster_covariance=cov.tolist(),parameter_se=np.sqrt(np.maximum(0,np.diag(cov))).tolist(),
            t_critical_95=float(student_t.ppf(.975,g-1)),training_residual_sd_V_per_m=float(np.sqrt(residual@residual/(n-p)))),
        observations=observations,limitations=LIMITATIONS)
    metrics={}
    for split,data in [('train',train),('holdout',hold)]:
        actual=np.array([r['field_V_per_m'] for r in data]); pred=design(data)@b
        errors=actual-pred
        metrics[split]=dict(n=len(data),rmse_V_per_m=float(np.sqrt(np.mean(errors**2))),mae_V_per_m=float(np.mean(abs(errors))),
            bias_prediction_minus_observed_V_per_m=float(np.mean(-errors)),
            training_mean_baseline_rmse_V_per_m=float(np.sqrt(np.mean((actual-y.mean())**2))))
    for item in observations:
        item['prediction']=predict_skin(result,item['age_group_table'],item['sex'],item['site'])
    result['metrics']=metrics
    result['groups']=[dict(age_group=age,sex=sex,site=site,**predict_skin(result,age,sex,site))
        for age in ('18-29','65-80') for sex in ('male','female') for site in ('control','LA-1','LA-2','LL-1')]
    return result


def build_skin_calibration(root):
    root=Path(root).resolve(); source=root/'data/derived/integumentary/human-wound-observations.jsonl'
    rows=[json.loads(line) for line in source.read_text().splitlines() if line.strip()]
    provenance=[]
    for name,expected in sorted({(r['source_file'],r['source_sha256']) for r in rows}):
        path=root/name; actual=hashlib.sha256(path.read_bytes()).hexdigest()
        if actual!=expected: raise ValueError('raw source hash mismatch: '+name)
        provenance.append(dict(path=name,sha256=actual,verified=True))
    result=fit_skin_observations(rows)
    result['provenance']=dict(observations_path=str(source.relative_to(root)),observations_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        raw_sources=provenance,references=sorted({r['reference'] for r in rows}),source_kind='human_measurement')
    result['observation_binding']=dict(model='ihm.materialize.skin',operator='ihm.calibration.skin_surface_field',
        expression='E_ab = (phi_apical[a] - phi_apical[b]) / distance_m',input_unit='V',output_unit='V/m',
        excluded_input='cell membrane_voltage',registration_status='no subject-matched electrode coordinates or patch geometry; quantity-level binding only')
    out=root/'data/derived/calibration'; out.mkdir(parents=True,exist_ok=True)
    (out/'skin-fit.json').write_text(json.dumps(result,separators=(',',':'),allow_nan=False))
    index=dict(schema_version=1,models=[dict(id=result['id'],path='data/derived/calibration/skin-fit.json',
        quantity=result['quantity'],unit=result['unit'],source_kind='human_measurement',metrics=result['metrics'],
        rank=result['fit']['rank'],parameters_count=len(FEATURES),groups=result['groups'],provenance=result['provenance'],limitations=LIMITATIONS)])
    (out/'index.json').write_text(json.dumps(index,separators=(',',':'),allow_nan=False))
    return result
