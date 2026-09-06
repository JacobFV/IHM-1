"""Conditional descriptive fit to Billich/Levitan1969 TableII isotonic means.

Infused sodium is the predictor, not instantaneous finite-lumen concentration.
This is an independent perfusion-response prior, not a native epithelial law.
"""
import json
from pathlib import Path
import numpy as np

DATA=Path(__file__).resolve().parents[2]/'data/research/gi_colon/table_ii_isotonic.json'
OUTPUTS=('water_ml_min','na_meq_hr','k_meq_hr','cl_meq_hr')

def fit():
    data=json.loads(DATA.read_text());rows=data['rows']
    x=np.array([r['infused_na_meq_l'] for r in rows]);y=np.array([[r[k] for k in OUTPUTS] for r in rows])
    design=np.c_[np.ones(len(x)),x]
    beta=np.linalg.lstsq(design,y,rcond=None)[0]
    held=[]
    for i in range(len(x)):
        mask=np.arange(len(x))!=i
        b=np.linalg.lstsq(design[mask],y[mask],rcond=None)[0]
        held.append((design[i]@b-y[i]).tolist())
    return dict(coefficients=beta.tolist(),outputs=OUTPUTS,leave_one_concentration_out_error=held,
                uncertainty='Equal-weight group-mean fit; SD descriptive, n repeated periods, no population posterior',
                temperature_k=data['temperature_k'],protocol=data['protocol'])

def response(infused_na_meq_l,*,flow_ml_hr,protocol,transfer_prior=False):
    if not np.isfinite([infused_na_meq_l,flow_ml_hr]).all() or infused_na_meq_l<0 or flow_ml_hr<=0:raise ValueError('Invalid perfusion condition')
    model=fit();matched=(0<=infused_na_meq_l<=150 and 820<=flow_ml_hr<=880 and protocol==model['protocol'])
    if not matched and not transfer_prior:raise ValueError('Outside acquired protocol; explicit transfer-prior opt-in required')
    rates=np.array([1,infused_na_meq_l])@np.array(model['coefficients'])
    return dict(rates=dict(zip(OUTPUTS,rates.tolist())),transfer_prior=not matched,
                temperature_k=None,temperature_status='unreported; thermal transfer uncalibrated',
                charge_equivalent_meq_hr=float(rates[1]+rates[2]-rates[3]),
                missing_charge_pathways='Bicarbonate/H and epithelial current not resolved; no chloride compensation')

def exchange(water_ml,ion_mol,duration_s,condition):
    """Paired finite lumen/serosal-bath exchange; no persistent native duplicate.

    Ion order Na,K,Cl. One common time fraction caps the entire signed empirical
    request when any donor exhausts. This preserves the requested rate ratios.
    """
    water=np.array(water_ml,float,copy=True);ions=np.array(ion_mol,float,copy=True)
    if water.shape!=(2,) or ions.shape!=(2,3) or not np.isfinite(water).all() or not np.isfinite(ions).all() or (water<0).any() or (ions<0).any() or not np.isfinite(duration_s) or not 0<=duration_s<=60:raise ValueError('Finite baths and bounded0–60s step required')
    r=condition['rates'];request=np.array([r['water_ml_min']/60,*[r[k]*1e-3/3600 for k in OUTPUTS[1:]]])*duration_s
    if not np.isfinite(request).all():raise ValueError('Nonfinite empirical rate request')
    pools=np.c_[water,ions];alpha=1.
    for j,q in enumerate(request):
        if q:alpha=min(alpha,pools[0 if q>0 else 1,j]/abs(q))
    moved=alpha*request
    result=pools.copy();result[0]-=moved;result[1]+=moved
    while (result<0).any():
        alpha=np.nextafter(alpha,0.);moved=alpha*request;result=pools.copy();result[0]-=moved;result[1]+=moved
    return dict(water_ml=result[:,0],ion_mol=result[:,1:],signed_lumen_to_serosa=moved,
                applied_time_fraction=alpha,uncompensated_charge_mol=float(moved[1]+moved[2]-moved[3]),
                native_stores_mutated=False)


def rectal_isotonic_observation():
    """Measured null uptake boundary, not a rectal permeability calibration."""
    return dict(water_absorption_ml_min=0.,na_absorption_meq_hr=0.,cl_absorption_meq_hr=0.,
                exposure_minutes_up_to=90,protocol='human rectal isotonic instillation/perfusion',
                source_doi='10.1136/gut.11.5.438',uncertainty=None,
                limitation='No measurable absorption reported; detection limit, exact composition and temperature not recovered. No zero-permeability assertion.')
