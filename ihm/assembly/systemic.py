"""Shared native meal, breathing, and exertion experiments; no shadow physiology."""
from dataclasses import dataclass, asdict
import json
import math
from pathlib import Path
from ihm.native import BASE, _sha
from ihm.native.session import NativeSession, SessionConfig, Meal, _ticks

PROTOCOLS = ('rest', 'hydration', 'meal', 'apnea', 'exercise', 'meal_exercise')


@dataclass(frozen=True)
class SystemicConfig:
    protocol: str = 'rest'
    seconds: float = 21600
    sample_interval_s: float = 30
    engine_variant: str = 'whole_body_integrity_renal'
    state_path: str | Path = BASE/'data/derived/canonical/native_baseline_v1/states/native_stabilized.xml'

    def __post_init__(self):
        if self.protocol not in PROTOCOLS:
            raise ValueError('Unknown systemic protocol')
        ticks = _ticks(self.seconds, 86400, 'seconds')
        cadence = _ticks(self.sample_interval_s, self.seconds, 'sample_interval_s')
        if ticks % cadence:
            raise ValueError('Duration must contain whole sample intervals')
        if self.protocol == 'apnea' and self.seconds < 180:
            raise ValueError('Apnea protocol requires at least 180 s')
        if self.protocol in ('exercise', 'meal_exercise') and self.seconds < 3000:
            raise ValueError('Exercise protocol requires at least 3000 s')
        SessionConfig(state_path=self.state_path, engine_variant=self.engine_variant, horizon_s=self.seconds)


def protocol_events(config):
    events = []
    if config.protocol in ('meal', 'meal_exercise'):
        events.append(dict(time_s=0, kind='meal', value=asdict(Meal(carbohydrate_g=60, protein_g=20,
            fat_g=20, sodium_g=1, calcium_mg=300, water_ml=500))))
    if config.protocol == 'hydration':
        events.append(dict(time_s=0, kind='meal', value=asdict(Meal(name='matched_hydration',
            sodium_g=1, calcium_mg=300, water_ml=500))))
    if config.protocol == 'apnea':
        events.extend([dict(time_s=30, kind='apnea', value=1), dict(time_s=90, kind='apnea', value=0)])
    if config.protocol in ('exercise', 'meal_exercise'):
        events.extend([dict(time_s=1800, kind='exercise', value=.15), dict(time_s=2400, kind='exercise', value=0)])
    return events


MECHANISMS = [
    dict(source='arterial_o2_mmhg + arterial_co2_mmhg', target='respiratory_request_cmh2o + respiratory_request_per_min',
         implementation='Nervous::ChemoreceptorFeedback', status='executing native reduced brainstem feedback'),
    dict(source='respiratory request', target='applied pressure, airflow, lung volume, blood gases',
         implementation='Respiratory::RespiratoryDriver + respiratory circuit + gas transport', status='executing native'),
    dict(source='stomach nutrients', target='SmallIntestineChyme',
         implementation='Gastrointestinal::DigestNutrient', status='executing lumped digestion; kinetics not independently calibrated'),
    dict(source='SmallIntestineChyme', target='SmallIntestineVasculature → portal path → LiverVasculature',
         implementation='Gastrointestinal::AbsorbNutrients + cardiovascular advection', status='executing native; lipid lymph bypass remains upstream limitation'),
    dict(source='circulating glucose', target='insulin/glucagon synthesis, liver/muscle glycogen',
         implementation='Endocrine + Hepatic + Tissue', status='executing shared native solute state'),
    dict(source='exercise demand', target='metabolism, gas demand, circulation, fatigue, thermal state',
         implementation='SEExercise + Energy + Tissue + Cardiovascular + Nervous', status='executing lumped demand; not computed skeletal muscle work'),
    dict(source='vascular fluid and solutes', target='glomerular filtration, urine storage/output',
         implementation='Renal + cardiovascular transport', status='executing native; renal transfer audit tracked separately'),
]


def _field(name):
    suffixes = [('concentration_mg_per_dl','mg/dL'), ('mass_g','g'), ('ml_per_min','mL/min'),
                ('l_per_min','L/min'), ('l_per_s','L/s'), ('pmol_per_min','pmol/min'),
                ('per_min','1/min'), ('mmhg','mmHg'), ('cmh2o','cmH2O'),
                ('_ml','mL'), ('_mg','mg'), ('_g','g'), ('_w','W'), ('_c','degC')]
    unit = next((unit for suffix,unit in suffixes if name.endswith(suffix)), '1')
    return dict(unit=unit, owner='BioGears', support=name.split('.')[0] if '.' in name else 'native systemic output',
                evidence='source-model state, not subject measurement', missing='null means absent/nonfinite native quantity')


def run_systemic(root, output_dir, config=None, **options):
    root = Path(root).resolve()
    config = config or SystemicConfig(**options)
    if not isinstance(config, SystemicConfig):
        raise ValueError('Expected SystemicConfig')
    out = Path(output_dir).resolve()
    if out.exists() and any(out.iterdir()):
        raise ValueError('Systemic experiment requires a fresh output directory')
    out.mkdir(parents=True, exist_ok=True)
    sources = ['ihm/assembly/systemic.py', 'ihm/native/session.py', 'scripts/native_body_ports.h',
               'scripts/native_biogears_stream.cpp', 'data/runtime/physiology/native_biogears_stream']
    from ihm.native import RUNTIME
    library=(RUNTIME/'biogears-build/outputs/Release/lib' if config.engine_variant=='upstream'
             else RUNTIME/'variants'/config.engine_variant)/'libbiogears.so.8.0.0'
    sources.extend([str(library.relative_to(root)), str(Path(config.state_path).resolve().relative_to(root))])
    source_hashes={p:_sha(root/p) for p in sources}
    (out/'source_receipts.json').write_text(json.dumps(source_hashes,indent=2)+'\n')
    events = protocol_events(config)
    total, stride = round(config.seconds*50), round(config.sample_interval_s*50)
    schedule = {round(e['time_s']*50): e for e in events}
    sample_ticks = set(range(0, total+1, stride))
    boundaries = sorted(sample_ticks | set(schedule))
    frames = []
    with NativeSession(SessionConfig(state_path=config.state_path, engine_variant=config.engine_variant,
                                    horizon_s=config.seconds), out/'native') as body:
        initial_state = body.save_state()
        current = 0
        for tick in boundaries:
            snapshot = body.step((tick-current)*.02) if tick > current else body.snapshot()
            current = tick
            if tick in sample_ticks:
                frame = dict(time_s=tick*.02, values=snapshot['values'])
                frames.append(frame)
                with (out/'frames.jsonl').open('a') as stream:
                    stream.write(json.dumps(frame, allow_nan=False)+'\n')
            if tick in schedule:
                event = schedule[tick]
                if event['kind']=='meal':
                    body.meal(Meal(**event['value']))
                else:
                    getattr(body, event['kind'])(event['value'])
        final_state = body.save_state()
    manifest = json.loads((out/'native/manifest.json').read_text())
    fields = {name:_field(name) for name in frames[0]['values']}
    stats = {}
    for name in fields:
        values = [frame['values'][name] for frame in frames]
        present = [v for v in values if v is not None]
        stats[name] = dict(initial=values[0], final=values[-1], minimum=min(present) if present else None,
                           maximum=max(present) if present else None, missing_count=len(values)-len(present))
    # Local stores are checked separately from concentrations and flows; no aggregate
    # sum double-counts native overlapping parent compartments.
    required_stores=['stomach_'+x for x in ('carbohydrate_g','protein_g','fat_g','sodium_g','calcium_mg','water_ml')]
    required_stores+=['liver_glycogen_g','muscle_glycogen_g','stored_protein_g','stored_fat_g']
    missing_required=[name for name in required_stores if name not in stats or stats[name]['missing_count']]
    negative = {name:stat['minimum'] for name,stat in stats.items()
                if ('.mass_g' in name or name in required_stores) and stat['minimum'] is not None and stat['minimum'] < -1e-10}
    result = dict(schema='ihm.systemic.v1', configuration=asdict(config), fields=fields, frames=frames,
                  actions=events, mechanism_edges=MECHANISMS, summary=stats,
                  checks=dict(sample_count=len(frames), expected_sample_count=total//stride+1,
                              negative_local_stores=negative, missing_required_stores=missing_required,
                              local_store_nonnegativity_passed=not negative and not missing_required),
                  runtime_sources=source_hashes, native_manifest=manifest,
                  initial_state_sha256=_sha(initial_state), final_state_sha256=_sha(final_state),
                  experimental_validation=False, limitations=[
                      'Generic source simulation, with no subject-specific calibration or parameter covariance.',
                      'Native gut is a lumped stomach/small intestine nutrient model; no complete chewing, peristalsis, colon or stool process.',
                      'Compartment storage changes mix local conversion with advection; they do not independently measure absorption flux.',
                      'Native metabolic conversions do not supply a closed elemental or exact first-law whole-body ledger.',
                      'Native brainstem requests execute; they are not the IBM cortical model or individual phrenic axons.',
                      'No online force/work return from canonical mechanics is included in this experiment.',
                      'Macro digestion rates and renal transport limitations remain separate source audits.'])
    (out/'systemic.json').write_text(json.dumps(result, indent=2, default=str, allow_nan=False)+'\n')
    return result


def verify_contrasts(results):
    """Require measured downstream changes; report magnitudes without clinical bounds."""
    checks = {}
    missing_controls=[]
    required={'meal':'hydration','apnea':'rest','exercise':'rest','meal_exercise':'meal'}
    for protocol,control in required.items():
        if protocol in results and control not in results:missing_controls.append(dict(protocol=protocol,control=control))
    for name, result in results.items():
        checks[name+'_local_stores'] = result['checks']['local_store_nonnegativity_passed']
    def changed(a,b,key,threshold):
        x,y=results[a],results[b]
        for identity in ('state_sha256','library_sha256','executable_sha256','patient_identity_input_sha256'):
            if not x['native_manifest'].get(identity) or x['native_manifest'][identity]!=y['native_manifest'].get(identity):
                raise ValueError('Causal contrasts require matching native '+identity)
        for identity in ('dependency_sha256','native_step_s'):
            if not x['native_manifest'].get(identity) or x['native_manifest'][identity]!=y['native_manifest'].get(identity):
                raise ValueError('Causal contrasts require matching native '+identity)
        if x['runtime_sources']!=y['runtime_sources']:
            raise ValueError('Causal contrasts require matching executing sources')
        if x['frames'][0]['values']!=y['frames'][0]['values']:
            raise ValueError('Causal contrasts require identical initial native observations')
        xa={event['time_s']:event for event in x['actions']}
        ya={event['time_s']:event for event in y['actions']}
        differences=[time for time in xa.keys()|ya.keys() if xa.get(time)!=ya.get(time)]
        if not differences:raise ValueError('Contrast actions are identical')
        intervention=min(differences)
        if [f['time_s'] for f in x['frames']] != [f['time_s'] for f in y['frames']]:
            raise ValueError('Causal contrasts require matching sample clocks')
        for u,v in zip(x['frames'],y['frames']):
            if u['time_s']<=intervention and u['values']!=v['values']:
                raise ValueError('Trajectories differ before the contrasted intervention')
        pairs=[(u['values'][key],v['values'][key]) for u,v in zip(x['frames'],y['frames']) if u['time_s']>intervention]
        if not pairs:raise ValueError('Contrast has no post-intervention samples')
        if any(u is None or v is None for u,v in pairs):
            raise ValueError('Causal contrast has missing native values')
        magnitude=max(abs(u-v) for u,v in pairs)
        checks[f'{a}_vs_{b}.{key}']=dict(maximum_absolute_difference=magnitude,passed=magnitude>threshold)
    if 'meal' in results and 'hydration' in results:
        for key, threshold in [('Aorta.Glucose.concentration_mg_per_dl', .01),
                              ('insulin_synthesis_pmol_per_min',.01), ('liver_glycogen_g',.01)]:
            changed('meal','hydration',key,threshold)
    if 'exercise' in results and 'rest' in results:
        changed('exercise','rest','metabolic_rate_w',1)
        changed('exercise','rest','oxygen_consumption_ml_per_min',1)
    if 'meal_exercise' in results and 'meal' in results:
        changed('meal_exercise','meal','metabolic_rate_w',1)
        changed('meal_exercise','meal','oxygen_consumption_ml_per_min',1)
    if 'apnea' in results and 'rest' in results:
        changed('apnea','rest','lung_volume_ml',50)
        changed('apnea','rest','arterial_co2_mmhg',.1)
    count=sum(isinstance(value,dict) for value in checks.values())
    return dict(checks=checks,missing_controls=missing_controls,causal_contrast_count=count,
                passed=bool(results) and count>0 and not missing_controls and all(v if isinstance(v,bool) else v['passed'] for v in checks.values()))
