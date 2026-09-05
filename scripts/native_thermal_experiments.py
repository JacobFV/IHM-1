#!/usr/bin/env python3
"""Matched-state hour thermal experiments; all source equations are retained."""
from pathlib import Path
import argparse,json,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.native import BASE,NativeConfig,run_native,summarize
SETTINGS={'control':(22,.5),'insulated':(22,1.5),'warm':(26,1.)}
def verify():
    base=BASE/'data/derived/physiology';report={}
    for case in SETTINGS:
        out=base/('native_hour_thermal_'+case)
        config=NativeConfig.from_dict(json.loads((out/'configuration.json').read_text()))
        summary=summarize(out,config)
        assert summary['rows']==3600 and summary['time_end_s']==3600
        report[case]={k:summary[k] for k in ('configuration','rows','csv_sha256','nutrition_context')}
        report[case]['quantities']={k:v for k,v in summary['summary'].items() if any(x in k for x in ('Temperature','BloodVolume','Metabolic','HeatLoss'))}
    assert report['insulated']['quantities']['CoreTemperature(degC)']['final']>report['control']['quantities']['CoreTemperature(degC)']['final']
    report['interpretation']='Insulation and room temperature reduce source-model cooling, but none of these declared protocols establishes hourly thermal equilibrium. Original equations retained; no clinical calibration. Default stomach water is a fluid input.'
    (base/'native_hour_thermal_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print('Three finite, complete 3600 s runs; insulation response verified; no equilibrium claim.')
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--case',choices=SETTINGS);parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    if args.verify:return verify()
    if not args.case:parser.error('--case or --verify required')
    ambient,clo=SETTINGS[args.case]
    config=NativeConfig(seconds=3600,sample_hz=1,state_path=str(BASE/'data/derived/physiology/native_hour_rest/states/native_stabilized.xml'),ambient_temperature_c=ambient,clothing_clo=clo)
    summary=run_native(config,BASE/'data/derived/physiology'/('native_hour_thermal_'+args.case))
    print(json.dumps({'configuration':summary['configuration'],'summary':summary['summary']},indent=2))
if __name__=='__main__':main()
