"""Native boundary output versus independent dimensional budgets."""
from pathlib import Path
import argparse,csv,json,sys,tempfile
import numpy as np
BASE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(BASE))
from scripts.audit_long_horizon_thermal import sha
from ihm.assembly.systemic_evidence import resolve_sources
from ihm.assembly.native_environment_evidence import resolve_native_environment

def verify(directory,humidity=False):
    folder=Path(directory).resolve();execution=json.loads((folder/'execution.json').read_text());assert execution['returncode']==0
    sources=resolve_sources(BASE,folder,execution['source_hashes']);environment=resolve_native_environment(BASE,folder)
    rows=[{k:float(v) for k,v in r.items()} for r in csv.DictReader((folder/'evaporation.csv').open())]
    assert len(rows)==4*3*2*3*6 and np.isfinite([[v for v in r.values()] for r in rows]).all()
    fractions=[.36,.07,.092,.092,.193,.193];errors=[];derivatives=[]
    for r in rows:
        i=int(r['region']);area=r['area_m2'];sweat=r['sweat_kg_s']*r['latent_j_kg']
        # Area-only variant deliberately inherits the old whole-torso saturation
        # and unscaled ambient saturation. Humidity variant must use each actual
        # regional skin temperature and relative humidity exactly once.
        torso=next(a for a in rows if all(a[k]==r[k] for k in ('area_m2','rh','pattern','sweat_kg_s')) and a['region']==0)
        sp=r['skin_saturation_pa'] if humidity else torso['skin_saturation_pa']
        ambient=r['ambient_saturation_pa']*(r['rh'] if humidity else 1)
        capacity=max(0.,(sp-ambient)/1000*r['evap_coefficient_w_m2_kpa'])
        sweat_flux=min(sweat/area,capacity);wet=sweat_flux/capacity if capacity else 0
        expected=(sweat_flux+(1-wet)*.06*capacity)*area*fractions[i]
        errors.append(abs(expected-r['evaporation_w']))
    assert max(errors)<1e-8,f'Native regional evaporative law differs by {max(errors)} W'
    for area in sorted(set(r['area_m2'] for r in rows)):
        for rh in (0.,.5,1.):
            totals=[]
            for sweat in (0.,1e-6):
                subset=[r for r in rows if r['area_m2']==area and r['rh']==rh and r['pattern']==0 and r['sweat_kg_s']==sweat]
                totals.append(sum(r['evaporation_w'] for r in subset))
            latent=rows[0]['latent_j_kg']
            ratio=(totals[1]-totals[0])/(1e-6*latent*.94)
            assert abs(ratio-1)<1e-10,'Whole-body incremental latent heat was multiplied by area'
            derivatives.append(dict(area_m2=area,rh=rh,whole_body_increment_ratio=ratio))
    out=Path(tempfile.mkdtemp(prefix='native-evaporation-verification-',dir=BASE/'data/derived/audits'))
    report=dict(status='passed',humidity_and_local_temperature=humidity,maximum_heat_error_w=max(errors),latent_increment_checks=derivatives,
        sources=sources,environment_sources=environment,csv_sha256=sha(folder/'evaporation.csv'),verifier_sha256=sha(__file__),
        interpretation='Direct native constitutive probe, without advancing physiology. No fluid balance or exercise calibration claim.')
    (out/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(status='passed',maximum_error_w=max(errors),output_dir=str(out)),indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('directory',type=Path);p.add_argument('--humidity',action='store_true');args=p.parse_args();verify(args.directory,args.humidity)
