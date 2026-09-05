"""Audit uninterrupted rest clocks, thermal energy and paired source identities."""
from pathlib import Path
import argparse,csv,json,sys,tempfile,xml.etree.ElementTree as ET
import numpy as np
BASE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(BASE))
from scripts.audit_long_horizon_thermal import state_audit,sha,tag,scalar
from ihm.assembly.systemic_evidence import resolve_sources
from ihm.assembly.native_environment_evidence import resolve_native_environment

def verify_case(case):
    case=Path(case).resolve()
    execution=json.loads((case/'execution.json').read_text());assert execution['returncode']==0
    config=json.loads((case/'configuration.json').read_text())
    sources=resolve_sources(BASE,case,execution['source_hashes']);environment=resolve_native_environment(BASE,case)
    initial=state_audit(case/'states/native_stabilized.xml');final=state_audit(case/'states/native_final.xml')
    rows=list(csv.DictReader((case/'native_multisystem.csv').open()));times=np.array([float(r['#Time(s)']) for r in rows])
    assert len(rows)==int(config['seconds']) and np.array_equal(times,np.arange(1,config['seconds']+1))
    assert np.isfinite(np.array([[float(v) for v in row.values()] for row in rows])).all()
    fixed_skin=sum(v for p,v in initial['metabolic_paths_w'].items() if p!='GroundToInternalCore')
    assert abs(fixed_skin-sum(v for p,v in final['metabolic_paths_w'].items() if p!='GroundToInternalCore'))<1e-9
    net=np.array([float(r['TotalMetabolicRate(W)'])+fixed_skin-sum(float(r[k]) for k in ('ConvectiveHeatLoss(W)','RadiativeHeatLoss(W)','EvaporativeHeatLoss(W)','RespirationHeatLoss(W)')) for r in rows])
    integrated=float(np.trapezoid(np.r_[initial['storage_rate_w'],net],np.r_[0,times]));change=final['stored_heat_j']-initial['stored_heat_j'];error=change-integrated
    assert abs(error)<max(5.,.002*abs(change)),'Sampled thermal flux and stored energy disagree'
    assert abs(final['equal_skin_temperature_equivalent_clo']-.5)<1e-10
    for label,path in [('Convective','ClothingToEnvironment'),('Radiative','ClothingToEnclosure')]:
        h=float(final['film_coefficients'][label+'HeatTranferCoefficient']['value']);area=float(final['environment']['SkinSurfaceArea']['value'])
        assert abs(final['resistance_k_w'][path]*h*area-1)<1e-5
    xml=ET.parse(case/'states/native_final.xml').getroot()
    vals={tag(e):e.attrib for e in xml.iter() if tag(e) in ('MeanSkinFlow','BloodDensity','BloodSpecificHeat')}
    flow=float(vals['MeanSkinFlow']['value'])*1e-6;rho=float(vals['BloodDensity']['value']);cp=float(vals['BloodSpecificHeat']['value'])
    paths={next(c.text for c in e if tag(c)=='Name'):{tag(c):c.attrib for c in e if c.attrib} for e in xml.iter() if tag(e)=='ThermalPath'}
    regional=[1/(scalar(v['Resistance'],'resistance')*.42*rho*cp) for name,v in paths.items() if name.startswith('InternalCoreToInternal') and name.endswith('Skin')]
    ratio=sum(regional)/flow;assert len(regional)==6
    expected=6 if config['variant']=='whole_body_integrity_thermal_boundary_v2' else 1
    assert abs(ratio-expected)<1e-3,'Regional thermal perfusion identity changed'
    core=np.array([float(r['CoreTemperature(degC)']) for r in rows])
    return dict(case=str(case.relative_to(BASE)),source_hashes=sources,environment_hashes=environment,configuration=config,
        initial=initial,final=final,regional_perfusion_sum_ratio=ratio,native_csv_sha256=sha(case/'native_multisystem.csv'),
        integrated_net_heat_j=integrated,stored_heat_change_j=change,sampled_integral_error_j=error,
        core_temperature_range_c=[float(core.min()),float(core.max())],
        late_linear_core_slope_c_per_hour=float(np.polyfit(times[-min(1800,len(times)):],core[-min(1800,len(core)):],1)[0]*3600),
        interpretation='Numerical and constitutive checks; no clinical calibration or all-system equilibrium implied.')

def main():
    p=argparse.ArgumentParser();p.add_argument('corrected',type=Path);p.add_argument('--control',type=Path);args=p.parse_args()
    cases={'corrected':verify_case(args.corrected)}
    if args.control:
        cases['control']=verify_case(args.control)
        a,b=[c['configuration'] for c in cases.values()]
        assert a['state_path'] and a['state_path']==b['state_path'],'Paired tests require the same initial state'
        key=str(Path(a['state_path']).resolve().relative_to(BASE))
        assert a['source_hashes'][key]==b['source_hashes'][key]
        assert a['seconds']==b['seconds']
        # Variant manifest objects must differ in exactly the regional-perfusion
        # Energy compilation unit; all other compiled inputs stay identical.
        manifests=[]
        for cfg in (a,b):
            receipt=json.loads((BASE/cases['corrected' if cfg is a else 'control']['case']/'frozen-sources.json').read_text())
            key=next(k for k in cfg['source_hashes'] if k.endswith('/manifest.json'))
            manifests.append(json.loads((BASE/receipt['copies'][key]['path']).read_text()))
        oa,ob=manifests[0]['object_sha256'],manifests[1]['object_sha256']
        sa,sb=set(oa),set(ob);assert len(sa-sb)==len(sb-sa)==1
        assert all(k.endswith('Energy.cpp.o') for k in sa^sb)
        assert all(oa[k]==ob[k] for k in sa&sb)
    out=Path(tempfile.mkdtemp(prefix='thermal-rest-verification-',dir=BASE/'data/derived/audits'))
    report=dict(status='passed',cases=cases,verifier_sha256=sha(__file__))
    (out/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status='passed',output_dir=str(out),core_temperature_c={k:v['final']['core_temperature_c'] for k,v in cases.items()},core_slope_c_hour={k:v['final']['core_temperature_rate_c_per_hr'] for k,v in cases.items()}),indent=2))

if __name__=='__main__':main()
