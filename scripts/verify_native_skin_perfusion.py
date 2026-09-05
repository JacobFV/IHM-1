"""Thermal regional perfusion must partition the single native skin inflow once."""
from pathlib import Path
import argparse,json,sys,tempfile,xml.etree.ElementTree as ET
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.audit_long_horizon_thermal import state_audit,tag,scalar,sha,BASE

def verify(path):
    root=ET.parse(path).getroot();values={tag(e):e.attrib for e in root.iter() if tag(e) in ('MeanSkinFlow','BloodDensity','BloodSpecificHeat')}
    assert values['MeanSkinFlow']['unit']=='mL/s' and values['BloodDensity']['unit']=='kg/m^3' and values['BloodSpecificHeat']['unit']=='J/K kg'
    flow=float(values['MeanSkinFlow']['value'])*1e-6;rho=float(values['BloodDensity']['value']);cp=float(values['BloodSpecificHeat']['value'])
    paths={next(c.text for c in e if tag(c)=='Name'):{tag(c):c.attrib for c in e if c.attrib} for e in root.iter() if tag(e)=='ThermalPath'}
    regional={name:1/(scalar(v['Resistance'],'resistance')*.42*rho*cp) for name,v in paths.items() if name.startswith('InternalCoreToInternal') and name.endswith('Skin')}
    ratio=sum(regional.values())/flow
    assert len(regional)==6 and abs(ratio-1)<1e-3,f'Regional thermal perfusion sums to {ratio:.9f} times native total skin inflow'
    thermal=state_audit(path)
    report={'status':'passed','source_path':str(Path(path)),'source_sha256':sha(path),'native_total_skin_flow_m3_s':flow,'regional_thermal_flow_m3_s':regional,'partition_ratio':ratio,'source_alpha_unchanged':.42,'thermal':thermal}
    out=Path(tempfile.mkdtemp(prefix='skin-perfusion-verification-',dir=BASE/'data/derived/audits'));(out/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'status':'passed','partition_ratio':ratio,'core_temperature_c':thermal['core_temperature_c'],'output_dir':str(out)},indent=2));return report
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('state',type=Path);verify(p.parse_args().state)
