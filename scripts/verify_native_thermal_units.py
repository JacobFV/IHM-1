#!/usr/bin/env python3
"""Native resistance-area regression with controlled synthetic state inputs."""
from pathlib import Path
import hashlib,json,math,re,sys,xml.etree.ElementTree as ET
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.native import BASE,NativeConfig,run_native,load_trajectory

def scalars(path):
    tree=ET.parse(path)
    for node in tree.iter():node.tag=node.tag.split('}')[-1]
    result={}
    for tag in ('SkinSurfaceArea','ClothingResistance','ConvectiveHeatTranferCoefficient','RadiativeHeatTranferCoefficient'):
        result[tag]=float(next(node for node in tree.iter(tag)).attrib['value'])
    for node in tree.iter('ThermalPath'):
        name=node.findtext('Name')
        if name in ('ClothingToEnvironment','ClothingToEnclosure'):
            result[name]=float(node.find('Resistance').attrib['value'])
    return result

def main():
    source=BASE/'data/derived/physiology/native_hour_rest/states/native_stabilized.xml'
    root=BASE/'data/derived/physiology/thermal_area_regression';root.mkdir(exist_ok=True)
    records={}
    for area in (1.,2.):
        text,count=re.subn(r'(<SkinSurfaceArea\b[^>]*value=")[^"]+',lambda match:match[1]+str(area),source.read_text())
        assert count==1
        state=root/f'controlled_area_{area}.xml';state.write_text(text)
        for variant in ('saturation_bounds_heatflux','saturation_bounds_heatflux_thermal_units'):
            out=root/f'{variant}_{area}'
            if not (out/'summary.json').exists():run_native(NativeConfig(seconds=.02,sample_hz=50,state_path=str(state),engine_variant=variant),out)
            values=scalars(out/'states/native_final.xml');records[f'{variant}_{area}']=values
            if variant.endswith('thermal_units'):
                for name,coefficient,factor in [('ClothingToEnvironment','ConvectiveHeatTranferCoefficient',.1),('ClothingToEnclosure','RadiativeHeatTranferCoefficient',5.)]:
                    expected=factor*max(values['ClothingResistance'],.1)/(area*values[coefficient])
                    assert math.isclose(values[name],expected,rel_tol=1e-12)
    for variant,expected_ratio in [('saturation_bounds_heatflux',2.),('saturation_bounds_heatflux_thermal_units',.5)]:
        a=records[f'{variant}_1.0'];b=records[f'{variant}_2.0']
        assert math.isclose(b['ClothingToEnvironment']/a['ClothingToEnvironment'],expected_ratio,rel_tol=1e-12)
    report={'source_state_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'controlled_change':'Only serialized SkinSurfaceArea changed to1 and2m²; synthetic numerical test, not physiological patient model.','records':records,'native_inverse_area_resistance_verified':True,'original_direct_area_defect_reproduced':True,'empirical_clothing_factors_retained':True}
    if '--hour' in sys.argv:
        report['hour_runs']={}
        for case,variant in [('control','saturation_bounds_heatflux'),('units','saturation_bounds_heatflux_thermal_units')]:
            out=BASE/'data/derived/physiology'/f'native_hour_thermal_{case}_heatflux'
            if not (out/'summary.json').exists():run_native(NativeConfig(seconds=3600,sample_hz=1,state_path=str(source),ambient_temperature_c=22,clothing_clo=.5,engine_variant=variant),out)
            summary=json.loads((out/'summary.json').read_text())
            assert summary['rows']==3600 and summary['time_end_s']==3600 and summary['engine_variant']==variant
            trajectory=load_trajectory(out)
            assert len(trajectory['time_s'])==3600 and trajectory['time_s'][-1]==3600
            assert hashlib.sha256((out/'native_multisystem.csv').read_bytes()).hexdigest()==summary['csv_sha256']
            report['hour_runs'][case]={'csv_sha256':summary['csv_sha256'],'core_temperature_degC':summary['summary']['CoreTemperature(degC)'],'skin_temperature_degC':summary['summary']['SkinTemperature(degC)'],'metabolic_power_W':summary['summary']['TotalMetabolicRate(W)']}
        report['hour_interpretation']='Dimensional correctness alone does not establish physiological validity; source empirical clothing factors retained, no fitted coefficients.'
    (root/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
