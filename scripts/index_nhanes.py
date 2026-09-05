"""index acquired measurements and estimate explicitly scoped population priors.

These estimates constrain population state distributions, never dynamic process
coefficients. Survey-weighted moments are not assay-error variances or SEs.
"""
import hashlib
import html
import json
from pathlib import Path
import re
import numpy as np
import pandas as pd
import pyreadstat
from ihm.forge.acquisition import sha256

RAW=Path('data/raw/population/nhanes-2017-2018')
OUT=Path('data/derived/population/nhanes-2017-2018')
CENSOR_FLAGS={'endocrine.insulin':'LBDINLC','immune.crp':'LBDHRPLC','blood.bilirubin':'LBDSTBLC'}
# file, variable, registered component, canonical unit, conversion scale, weight
BINDINGS=(
 ('GLU_J','LBDGLUSI','metabolic.glucose','mmol/L',1.,'WTSAF2YR'),
 ('INS_J','LBXIN','endocrine.insulin','mU/L',1.,'WTSAF2YR'),
 ('CBC_J','LBXHGB','blood.hemoglobin','g/L',10.,'WTMEC2YR'),
 ('CBC_J','LBXHCT','blood.hematocrit','1',.01,'WTMEC2YR'),
 ('CBC_J','LBXPLTSI','blood.platelets','1e9/L',1.,'WTMEC2YR'),
 ('CBC_J','LBDNENO','immune.neutrophils','cells/uL',1000.,'WTMEC2YR'),
 ('CBC_J','LBDMONO','immune.monocytes','cells/uL',1000.,'WTMEC2YR'),
 ('BIOPRO_J','LBDSALSI','blood.albumin','g/L',1.,'WTMEC2YR'),
 ('BIOPRO_J','LBDSBUSI','blood.urea','mmol/L',1.,'WTMEC2YR'),
 ('BIOPRO_J','LBDSCRSI','blood.creatinine','umol/L',1.,'WTMEC2YR'),
 ('BIOPRO_J','LBDSTBSI','blood.bilirubin','umol/L',1.,'WTMEC2YR'),
 ('BIOPRO_J','LBXSNASI','blood.sodium','mmol/L',1.,'WTMEC2YR'),
 ('BIOPRO_J','LBXSKSI','blood.potassium','mmol/L',1.,'WTMEC2YR'),
 ('BIOPRO_J','LBXSCLSI','blood.chloride','mmol/L',1.,'WTMEC2YR'),
 ('BIOPRO_J','LBXSC3SI','blood.total_co2','mmol/L',1.,'WTMEC2YR'),
 ('BIOPRO_J','LBDSCASI','blood.calcium','mmol/L',1.,'WTMEC2YR'),
 ('BIOPRO_J','LBDSPHSI','blood.phosphate','mmol/L',1.,'WTMEC2YR'),
 ('BIOPRO_J','LBDSTPSI','blood.total_protein','g/L',1.,'WTMEC2YR'),
 ('BIOPRO_J','LBDSUASI','blood.uric_acid','umol/L',1.,'WTMEC2YR'),
 ('BIOPRO_J','LBDSCHSI','blood.cholesterol','mmol/L',1.,'WTMEC2YR'),
 ('HSCRP_J','LBXHSCRP','immune.crp','mg/L',1.,'WTMEC2YR'),
 ('BPX_J','BPXPLS','vascular.pulse_rate','bpm',1.,'WTMEC2YR'),
 ('UCFLOW_J','URDFLOW1','urinary.collection_flow','mL/min',1.,'WTMEC2YR'),
)


def labels(path):
    if not path.exists(): return {}
    text=path.read_text()
    return {name:html.unescape(re.sub('<[^>]*>','',label)).strip()
            for name,label in re.findall(r'<h3[^>]+id="([^"]+)"[^>]*>(.*?)</h3>',text,re.S)}


def split(ids):
    return np.array([int(hashlib.sha256(('NHANES2017-2018:'+str(int(i))).encode()).hexdigest()[:8],16)%5!=0 for i in ids])


def moments(values,weights):
    mean=float(np.average(values,weights=weights)); var=float(np.average((values-mean)**2,weights=weights))
    order=np.argsort(values); sorted_x=np.asarray(values)[order]; sorted_w=np.asarray(weights)[order]
    q=np.interp([.025,.25,.5,.75,.975],(np.cumsum(sorted_w)-.5*sorted_w)/sorted_w.sum(),sorted_x)
    return {'mean':mean,'std':float(np.sqrt(var)),'variance':var,'quantiles':q.tolist(),
            'n':len(values),'weight_sum':float(np.sum(weights)),
            'kish_weight_effective_n':float(np.sum(weights)**2/np.sum(weights**2))}


def main():
    OUT.mkdir(parents=True,exist_ok=True); tables={}; inventory=[]; variable_count=0
    for path in sorted(RAW.glob('*.xpt')):
        frame,sas_metadata=pyreadstat.read_xport(str(path)); source_labels=sas_metadata.column_names_to_labels
        output=OUT/(path.stem+'.parquet'); frame.to_parquet(output,index=False)
        variables=[{'name':c,'label':source_labels.get(c),'nonmissing':int(frame[c].notna().sum()),
                    'dtype':str(frame[c].dtype)} for c in frame.columns]
        inventory.append({'id':path.stem,'raw_path':str(path),'sha256':sha256(path),
                          'rows':len(frame),'columns':len(frame.columns),'variables':variables,
                          'decoder':'pyreadstat '+pyreadstat.__version__,
                          'participant_unique':bool(frame.SEQN.is_unique) if 'SEQN' in frame else None,
                          'derived_path':str(output),'derived_sha256':sha256(output)})
        variable_count+=len(frame.columns)
        if path.stem in {b[0] for b in BINDINGS}|{'DEMO_J','BMX_J'}: tables[path.stem]=frame
    demo=tables['DEMO_J']; adult=demo.loc[demo.RIDAGEYR>=20].copy()
    adult['split']=np.where(split(adult.SEQN),'fit','holdout')
    # Keep survey design, age, sex and exact source participant IDs for every row.
    joined=adult[['SEQN','RIDAGEYR','RIAGENDR','SDMVPSU','SDMVSTRA','WTMEC2YR','split']].copy()
    fasting=tables['GLU_J'][['SEQN','WTSAF2YR']]
    joined=joined.merge(fasting,on='SEQN',how='left',validate='one_to_one')
    for filename,variable,component,unit,scale,weight in BINDINGS:
        if variable not in tables[filename]: raise ValueError(f'missing expected binding: {filename}/{variable}')
        value=tables[filename][['SEQN',variable]].rename(columns={variable:component})
        if component in CENSOR_FLAGS:
            flag=CENSOR_FLAGS[component]
            value[component+'__detection_flag']=tables[filename][flag]
        value[component]*=scale
        joined=joined.merge(value,on='SEQN',how='left',validate='one_to_one')
    joined.to_parquet(OUT/'adult-linked-measurements.parquet',index=False)
    estimates=[]
    for filename,variable,component,unit,scale,weight in BINDINGS:
        valid=joined.loc[joined[component].notna() & (joined[weight]>0)]
        for group,frame in [('adult_all',valid),('adult_male',valid.loc[valid.RIAGENDR==1]),('adult_female',valid.loc[valid.RIAGENDR==2])]:
            train=frame.loc[frame.split=='fit']; held=frame.loc[frame.split=='holdout']
            if len(train)<20 or len(held)<5: continue
            estimate=moments(train[component],train[weight]); holdout=moments(held[component],held[weight])
            estimates.append({'component':component,'unit':unit,'population':group,
                'source':'NHANES2017-2018','source_file':str(RAW/(filename+'.xpt')),
                'source_sha256':sha256(RAW/(filename+'.xpt')),'source_variable':variable,
                'source_label':labels(RAW/(filename+'.htm')).get(variable),'scale_to_canonical':scale,
                'reference':'https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/'+filename+'.htm',
                'weight':weight,'basis':'survey_weighted_empirical_population_state_prior',
                'calibrates_dynamics':False,
                'detection_limit_handling':{'policy':'retain CDC provided substituted values in moments; not exact uncensored measurements',
                    'flag_variable':CENSOR_FLAGS.get(component),
                    'fit_below_detection_n':int((train.get(component+'__detection_flag',pd.Series(dtype=float))==1).sum()),
                    'holdout_below_detection_n':int((held.get(component+'__detection_flag',pd.Series(dtype=float))==1).sum())},
                'estimate':estimate,'internal_holdout':holdout,
                'holdout_mean_difference':holdout['mean']-estimate['mean'],
                'uncertainty_semantics':'population dispersion, not measurement error or estimator standard error'})
    components=[b[2] for b in BINDINGS]
    complete=joined.loc[(joined.WTSAF2YR>0)].dropna(subset=components)
    train=complete.loc[complete.split=='fit']; held=complete.loc[complete.split=='holdout']
    x=train[components].to_numpy(float); w=train.WTSAF2YR.to_numpy(float)
    mean=np.average(x,axis=0,weights=w); delta=x-mean
    cov=(delta.T*w)@delta/w.sum()
    joint={'schema_version':1,'source':'NHANES2017-2018','population':'adults_age_20_plus_complete_cases_positive_fasting_weight',
        'components':components,'units':[b[3] for b in BINDINGS],'mean':mean.tolist(),'covariance':cov.tolist(),
        'decoder':'pyreadstat '+pyreadstat.__version__,
        'detection_limit_handling':{'policy':'CDC-provided substitution retained; censor flags preserved in linked measurements',
            'fit_below_detection_counts':{c:int((train[c+'__detection_flag']==1).sum()) for c in CENSOR_FLAGS}},
        'n_fit':len(train),'n_internal_holdout':len(held),'basis':'survey_weighted_empirical_population_state_prior',
        'parameter_calibration':False,'weight':'WTSAF2YR','complete_case_selection_bias_possible':True,
        'fit_subject_ids':[f'NHANES2017-2018:{int(i)}' for i in train.SEQN],
        'holdout_subject_ids':[f'NHANES2017-2018:{int(i)}' for i in held.SEQN],
        'bindings':[{'file':b[0],'variable':b[1],'component':b[2],'unit':b[3],'scale':b[4],
                     'sha256':sha256(RAW/(b[0]+'.xpt'))} for b in BINDINGS],
        'limitations':['cross-sectional state distribution does not identify dynamics',
          'below-detection substitutions retained: empirical moments do not fit a censored likelihood',
          'internal split of one survey; no external validation',
          'survey design standard errors not estimated; weights used for point moments',
          'adult cohort includes disease and treatment; not a healthy-only reference',
          'plasma/serum chemistry stays in blood components; no inferred interstitial equivalence',
          'pulse count and voided urine collection rate remain distinct from cardiac pacing and nephron flow',
          'Gaussian moment approximation may assign mass outside physiological bounds']}
    (OUT/'joint-population-prior.json').write_text(json.dumps(joint,indent=2,allow_nan=False)+'\n')
    with (OUT/'population-priors.jsonl').open('w') as stream:
        for row in estimates: stream.write(json.dumps(row,allow_nan=False)+'\n')
    report={'tables':inventory,'table_count':len(inventory),'variable_definitions':variable_count,
            'rows_across_tables_not_unique_people':sum(t['rows'] for t in inventory),
            'unique_participants':len(demo),'adults':len(adult),'state_prior_estimates':len(estimates),
            'decoder':'pyreadstat '+pyreadstat.__version__,
            'zero_weight_fasting_rows':int((tables['GLU_J'].WTSAF2YR==0).sum()),
            'zero_weight_MEC_rows':int((demo.WTMEC2YR==0).sum()),
            'joint_prior_fit_n':len(train),'joint_prior_holdout_n':len(held)}
    (OUT/'index.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='tables'}))


if __name__=='__main__': main()
