"""Read-only native exercise demand audit; records known source defects."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import xml.etree.ElementTree as ET
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.native.session import NativeSession, SessionConfig
ROOT=Path(__file__).resolve().parents[1]


def saved_energy(path):
    tree=ET.parse(path)
    wanted={'ExerciseEnergyDemand','TotalMetabolicRate','BasalMetabolicRate','MaxWorkRate','TotalWorkRateLevel','AchievedExerciseLevel','LactateProductionRate','MuscleGlycogen','CoreTemperature','FatigueLevel'}
    values={'exercise_action_active':any(n.attrib.get('{http://www.w3.org/2001/XMLSchema-instance}type')=='ExerciseData' for n in tree.iter())}
    for node in tree.iter():
        name=node.tag.split('}')[-1]
        if name in wanted:values[name]=dict(node.attrib)
    return values


def main():
    p=argparse.ArgumentParser();p.add_argument('--variant',default='whole_body_integrity');args=p.parse_args()
    state=ROOT/'data/derived/canonical/native_baseline_v1/states/native_stabilized.xml'
    out=Path(tempfile.mkdtemp(prefix='exercise-energy-',dir=ROOT/'data/derived/audits'))
    sources=ROOT/'data/raw/physiology/biogears/projects/biogears/libBiogears/src'
    source_paths=[sources/'engine/Systems/Energy.cpp',sources/'engine/Systems/Tissue.cpp',sources/'cdm/patient/actions/SEExercise.cpp',sources/'cdm/scenario/SEPatientActionCollection.cpp']
    sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
    receipt={'state_path':str(state),'state_sha256':sha(state),'variant':args.variant,'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in source_paths}}
    records={}
    for mode in ('rest','continuous','stop_30s'):
        samples=[]
        with NativeSession(SessionConfig(state_path=state,engine_variant=args.variant,horizon_s=120),out/mode) as session:
            session.exercise(0 if mode=='rest' else .15)
            previous=0
            for target in (0,.02,1,5,10,30,30.02,31,60,120):
                observation=session.snapshot() if target==0 else session.step(round(target-previous,8))
                saved=session.save_state()
                samples.append({'elapsed_s':target,'values':observation['values'],'saved_energy':saved_energy(saved),'state_path':str(saved),'state_sha256':sha(saved)})
                if target==30 and mode=='stop_30s':
                    ack=session.exercise(0)
                    receipt['stop_ack']=ack
                    saved=session.save_state();receipt['immediate_post_stop_energy']=saved_energy(saved)
                previous=target
        records[mode]=samples
        receipt.setdefault('session_manifest_sha256',{})[mode]=sha(out/mode/'manifest.json')
    def val(mode,t,key):return next(x for x in records[mode] if x['elapsed_s']==t)['values'][key]
    def demand(mode,t):return float(next(x for x in records[mode] if x['elapsed_s']==t)['saved_energy']['ExerciseEnergyDemand']['value'])
    checks={'matched_initial':records['rest'][0]['values']==records['continuous'][0]['values']==records['stop_30s'][0]['values'],
        'same_pre_stop_state':all(a['values']==b['values'] for a,b in zip(records['continuous'][:6],records['stop_30s'][:6])),
        'stop_acknowledged_and_action_removed':receipt['stop_ack']['status']=='ok' and not receipt['immediate_post_stop_energy']['exercise_action_active'],
        'demand_units_w':all(s['saved_energy']['ExerciseEnergyDemand']['unit']=='W' for mode in records.values() for s in mode),
        'demand_accumulation_reproduced':demand('continuous',120)>2*demand('continuous',30),
        'stop_retains_demand_reproduced':demand('stop_30s',30)==demand('stop_30s',120),
        'source_unchanged':all(sha(ROOT/p)==v for p,v in receipt['source_sha256'].items())}
    report={'audit_reproduced':all(checks.values()),'checks':checks,'receipt':receipt,'samples':records,
            'status':'Known energy-demand consistency defects reproduced; not a model acceptance pass.'}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(out);print(json.dumps(checks,indent=2))
    for mode in records:
        print(mode,[(x['elapsed_s'],x['saved_energy'].get('ExerciseEnergyDemand'),x['values']['metabolic_rate_w'],x['values']['oxygen_consumption_ml_per_min'],x['values']['muscle_glycogen_g']) for x in records[mode]])
    if not report['audit_reproduced']:raise SystemExit(1)

if __name__=='__main__':main()
