"""Read native thermal states and preserve explicit circuit energy/boundary evidence."""
from pathlib import Path
import argparse,csv,hashlib,json,tempfile,xml.etree.ElementTree as ET
import numpy as np
BASE=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def tag(e):return e.tag.split('}')[-1]
def scalar(e,quantity):
    value=float(e['value']);unit=e['unit']
    conversions={'power':{'W':1,'J/s':1,'kcal/hr':4184/3600,'kcal/day':4184/86400},
                 'capacity':{'J/K':1,'kcal/degC':4184,'kcal/K':4184},'resistance':{'K/W':1,'degC s/kcal':1/4184}}
    if quantity=='temperature':
        if unit not in ('K','degC'):raise ValueError('Unknown thermal temperature unit: '+unit)
        return value-273.15 if unit=='K' else value
    if unit not in conversions[quantity]:raise ValueError('Unknown '+quantity+' unit: '+unit)
    return value*conversions[quantity][unit]

def state_audit(path):
    root=ET.parse(path).getroot()
    def named(kind):
        return {next(c.text for c in e if tag(c)=='Name'):{tag(c):dict(c.attrib) if c.attrib else c.text for c in e} for e in root.iter() if tag(e)==kind}
    paths=named('ThermalPath');nodes=named('ThermalNode')
    def q(name,key='HeatTransferRate'):return scalar(paths[name][key],'power')
    capacity_paths=[p for p in paths if 'Capacitance' in paths[p]]
    metabolic_paths=[p for p in paths if p.startswith('GroundToInternal')]
    losses=['ClothingToEnvironment','ClothingToEnclosure','ExternalCoreToGround']+[p for p in paths if p.startswith('External') and p.endswith('SkinToGround')]
    storage=sum(q(p) for p in capacity_paths);inputs=sum(q(p) for p in metabolic_paths);outputs=sum(q(p) for p in losses)
    capacity={paths[p]['SourceNode']:scalar(paths[p]['Capacitance'],'capacity') for p in capacity_paths}
    temperatures={name:scalar(nodes[name]['Temperature'],'temperature') for name in capacity}
    if any(nodes[n]['Heat']['unit']!='J' for n in capacity):raise ValueError('Stored heat must be recorded in joules')
    clothing=[p for p in paths if p.startswith('External') and p.endswith('SkinToClothing')]
    resistance={p:scalar(paths[p]['Resistance'],'resistance') for p in clothing+['ClothingToEnvironment','ClothingToEnclosure']}
    environment=next(e for e in root if any(v=='BioGearsEnvironmentData' for v in e.attrib.values()))
    patient=next(e for e in root if tag(e)=='Patient')
    env={tag(e):dict(e.attrib) for e in environment.iter() if tag(e) in ('AmbientTemperature','ClothingResistance','MeanRadiantTemperature','AirVelocity','RelativeHumidity')}
    env.update({tag(e):dict(e.attrib) for e in patient if tag(e) in ('SkinSurfaceArea','Weight')})
    coefficients={tag(e):dict(e.attrib) for e in environment if tag(e) in ('ConvectiveHeatTranferCoefficient','RadiativeHeatTranferCoefficient')}
    area=float(env['SkinSurfaceArea']['value']);clo=float(env['ClothingResistance']['value'])
    if env['SkinSurfaceArea']['unit']!='m^2' or env['ClothingResistance']['unit']!='clo':raise ValueError('Unexpected area or clothing output units')
    req=1/sum(1/resistance[p] for p in clothing)
    report={'source_path':str(Path(path).resolve().relative_to(BASE)),'sha256':sha(path),
        'core_temperature_c':temperatures['InternalCore'],'node_temperature_c':temperatures,
        'capacities_j_k':capacity,'total_capacity_j_k':sum(capacity.values()),
        'stored_heat_j':sum(float(nodes[n]['Heat']['value']) for n in capacity),
        'storage_rate_w':storage,'metabolic_core_and_skin_w':inputs,'boundary_loss_w':outputs,
        'metabolic_paths_w':{p:q(p) for p in metabolic_paths},'film_coefficients':coefficients,
        'circuit_balance_residual_w':storage-inputs+outputs,'boundary_paths_w':{p:q(p) for p in losses},
        'core_storage_rate_w':q('InternalCoreToGround'),'core_temperature_rate_c_per_hr':q('InternalCoreToGround')/capacity['InternalCore']*3600,
        'resistance_k_w':resistance,'clothing_equivalent_resistance_k_w':req,
        'reported_clo':clo,'equal_skin_temperature_equivalent_clo':req*area/.155,'environment':env}
    assert abs(report['circuit_balance_residual_w'])<1e-3,'Thermal circuit imbalance requires investigation'
    return report

def main():
    p=argparse.ArgumentParser();p.add_argument('states',nargs='*',type=Path);args=p.parse_args()
    states=args.states or [BASE/'data/derived/canonical/native_baseline_v1/states/native_stabilized.xml',BASE/'data/derived/audits/fresh-energy-init-rqmr4of7/states/native_stabilized.xml']
    if not args.states:
        states+= [sorted(path.glob('*.xml'))[-1] for path in (BASE/'data/derived/systemic/six_hour_v3').glob('*/native/states')]
    report={'states':[state_audit(path) for path in states],'script_sha256':sha(__file__)}
    if not args.states:
        index=json.loads((BASE/'data/derived/thermal/index.json').read_text())
        report['jos3_held_comparators']=[{k:r[k] for k in ('id','configuration','final','audit')} for r in index['runs']]
        report['jos3_index_sha256']=sha(BASE/'data/derived/thermal/index.json')
        report['six_hour']={}
        initial=report['states'][0]
        for name in ('hydration','meal'):
            path=BASE/'data/derived/systemic/six_hour_v3'/name/'systemic.json';d=json.loads(path.read_text())
            frames=d['frames'];times=np.array([f['time_s'] for f in frames]);core=np.array([f['values']['core_temperature_c'] for f in frames])
            final=next(r for r in report['states'] if '/'+name+'/' in r['source_path'])
            report['six_hour'][name]={'source_sha256':sha(path),'initial_core_c':float(core[0]),'final_core_c':float(core[-1]),
                'first_sample_below_35c_s':float(times[np.flatnonzero(core<35)[0]]),
                'temperature_at_hour_c':[float(core[np.argmin(abs(times-3600*i))]) for i in range(7)],
                'stored_heat_change_j':final['stored_heat_j']-initial['stored_heat_j'],
                'mean_storage_rate_w':(final['stored_heat_j']-initial['stored_heat_j'])/21600,
                'interpretation':'Thermally invalid healthy-rest baseline despite local nutrient-store checks.'}
    out=Path(tempfile.mkdtemp(prefix='long-horizon-thermal-',dir=BASE/'data/derived/audits'));(out/'audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'output_dir':str(out),'states':[{k:r[k] for k in ('source_path','core_temperature_c','storage_rate_w','boundary_loss_w','circuit_balance_residual_w','equal_skin_temperature_equivalent_clo')} for r in report['states']]},indent=2))
if __name__=='__main__':main()
