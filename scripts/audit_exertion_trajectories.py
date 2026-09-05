"""Independent reading of recorded physiology; no cached acceptance predicates."""
from pathlib import Path
import csv,json,sys,tempfile
import numpy as np
BASE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(BASE))
from scripts.audit_long_horizon_thermal import sha
from ihm.assembly.systemic_evidence import freeze_sources

def audit():
    source=BASE/'data/derived/systemic/exertion_v3'
    out=Path(tempfile.mkdtemp(prefix='exertion-v3-physiology-audit-',dir=BASE/'data/derived/audits'))
    records={};inputs={str(Path(__file__).relative_to(BASE)):sha(Path(__file__))}
    for protocol in ('rest','hydration','meal','exercise','meal_exercise'):
        path=source/protocol/'systemic.json';raw=source/protocol/'frames.jsonl'
        data=json.loads(path.read_text());frames=data['frames'];journal=[json.loads(line) for line in raw.read_text().splitlines()]
        assert frames==journal,'Assembled trajectory differs from append-only observation journal'
        times=np.array([f['time_s'] for f in frames]);assert np.array_equal(times,np.arange(0,3601,5))
        stats={};missing={};negative={}
        for key in frames[0]['values']:
            a=np.array([f['values'][key] if f['values'][key] is not None else np.nan for f in frames]);valid=np.isfinite(a)
            if not valid.all():missing[key]=times[~valid].tolist()
            if not valid.any():continue
            lo=int(np.nanargmin(a));hi=int(np.nanargmax(a))
            stats[key]=dict(initial=float(a[0]) if valid[0] else None,final=float(a[-1]) if valid[-1] else None,
                minimum=float(a[lo]),minimum_time_s=float(times[lo]),maximum=float(a[hi]),maximum_time_s=float(times[hi]))
            if key.endswith(('.mass_g','.volume_ml','.concentration_mg_per_dl')) or key.startswith(('stomach_','stored_')) or key in ('liver_glycogen_g','muscle_glycogen_g'):
                if (a[valid]<0).any():negative[key]=float(a[valid].min())
        v=lambda key:np.array([f['values'][key] for f in frames],float)
        glucose=v('Aorta.Glucose.concentration_mg_per_dl');ph=v('arterial_ph')
        def sampled_interval(mask):
            t=times[mask]
            return None if not len(t) else dict(first_sample_s=float(t[0]),last_sample_s=float(t[-1]),sample_count=len(t),sampled_support_seconds=int(mask[:-1].sum()*5))
        partition=v('metabolic_rate_w')-v('basal_metabolic_rate_w')-v('exercise_energy_demand_w')
        selected={str(int(t)):{k:frames[int(t/5)]['values'][k] for k in stats} for t in (0,5,600,1200,1800,1805,1835,2100,2400,2405,2700,3000,3600)}
        records[protocol]=dict(configuration=data['configuration'],actions=data['actions'],field_statistics=stats,missing=missing,
            negative_observed_local_stores=negative,sampled_glucose_below54=sampled_interval(glucose<54),
            sampled_ph_below735=sampled_interval(ph<7.35),sampled_ph_above745=sampled_interval(ph>7.45),
            maximum_requested_power_partition_residual_w=float(abs(partition).max()),
            integrated_requested_metabolic_energy_j=float(np.trapezoid(v('metabolic_rate_w'),times)),
            integrated_requested_exercise_energy_j=float(np.trapezoid(v('exercise_energy_demand_w'),times)),selected_times=selected,
            interpretation='Observed source-model trajectory. Clinical glucose/pH comparisons are screening flags, not a diagnosis; aortic simulated concentration is not a validated clinical plasma measurement.')
        inputs.update({str(p.relative_to(BASE)):sha(p) for p in (path,raw,source/protocol/'native/manifest.json')})
    records['meal_minus_hydration']={k:records['meal']['field_statistics'][k]['final']-records['hydration']['field_statistics'][k]['final'] for k in ('liver_glycogen_g','muscle_glycogen_g','metabolic_rate_w','Aorta.Glucose.concentration_mg_per_dl','insulin_synthesis_pmol_per_min')}
    for p in ('ihm/app/__init__.py','ihm/native/session.py','ihm/native/__init__.py','ihm/assembly/systemic.py','app/src/main.js',
              'data/runtime/physiology/variants/whole_body_integrity_skin_perfusion/Energy.cpp',
              'data/raw/physiology/biogears/projects/biogears/libBiogears/src/engine/Systems/Tissue.cpp',
              'data/raw/physiology/biogears/projects/biogears/libBiogears/src/engine/Systems/Hepatic.cpp'):
        inputs[p]=sha(BASE/p)
    freeze_sources(BASE,out,inputs)
    report=dict(status='audit_completed_with_physiological_readiness_blockers',records=records,source_hashes=inputs,
        findings=['Exercise and meal-exercise recorded glucose below54mg/dL.','Exercise, hydration and meal acid-base excursions need physiological/source review.','Meal thermogenesis is absent from requested metabolic power.','Rest fixture contains baseline gastric water/sodium/calcium.'],
        recommendation='Retain causal/model experiments; do not promote as healthy generic default. Separate numerical/causal checks from physiological plausibility and human calibration.')
    (out/'audit.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'output_dir':str(out),'status':report['status'],'minimum_glucose':{k:r['field_statistics']['Aorta.Glucose.concentration_mg_per_dl']['minimum'] for k,r in records.items() if 'field_statistics' in r}},indent=2))
    return out

if __name__=='__main__':audit()
