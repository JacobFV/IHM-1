#!/usr/bin/env python3
"""Verify preserved original failure and isolated source bounds correction artifacts."""
import json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.native import BASE,load_trajectory

def main():
    base=BASE/'data/derived/physiology'
    original=json.loads((base/'native_female_init_repro/execution.json').read_text())
    assert original['exit_code']==-6
    for name in ['native_saturation_bounds_StandardMale','native_saturation_bounds_StandardFemale','native_female_fluid_corrected']:
        summary=json.loads((base/name/'summary.json').read_text())
        assert summary['rows']==20 and summary['time_end_s']==2 and summary['all_finite']
        assert summary['engine_variant']=='saturation_bounds'
    male_equal=load_trajectory(base/'native_saturation_bounds_StandardMale')==load_trajectory(base/'native_saturation_control_StandardMale')
    assert male_equal
    report={'original_female_exit_code':-6,'corrected_female_rest_pass':True,'corrected_female_fluid_protocol_pass':True,'male_original_corrected_bitwise_csv_values_equal':male_equal,'external_clinical_validation':False}
    heat=base/'native_heatflux_corrected'
    if (heat/'summary.json').exists():
        a=load_trajectory(heat);b=load_trajectory(base/'native_heatflux_control')
        assert a['time_s']==b['time_s']
        assert all(a['values'][key]==b['values'][key] for key in a['values'] if key!='EvaporativeHeatLoss(W)')
        assert a['values']['EvaporativeHeatLoss(W)']!=b['values']['EvaporativeHeatLoss(W)']
        report['heatflux_patch_all_other_channels_bitwise_equal']=True
    (base/'native_saturation_variant_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
