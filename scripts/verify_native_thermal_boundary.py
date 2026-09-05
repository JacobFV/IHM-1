"""Native area/insulation subdivision and dry boundary conservation checks."""
from pathlib import Path
import argparse,csv,subprocess,tempfile,json,hashlib,sys
import numpy as np
BASE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(BASE))

def verify_algebra():
    header=BASE/'scripts/native_thermal_boundary.h'
    assert header.exists(),'Native thermal boundary constitutive functions missing'
    out=Path(tempfile.mkdtemp(prefix='thermal-boundary-algebra-',dir=BASE/'data/derived/audits'))
    source=out/'test.cpp'
    source.write_text('''#include "'''+str(header)+'''"
#include <cassert>
#include <iostream>
int main() {
  const double A=1.9, I=.155*.5, floor=1e-8;
  double fractions[]={.36,.07,.092,.092,.193,.193};
  assert(ihm_thermal::segment_resistance(0.,A,.36,1e-100)>=1e-8/.36);
  for (double insulation: {0.,I,.155*2.}) {
    double conductance=0.;
    for (double f: fractions) conductance+=1/ihm_thermal::segment_resistance(insulation,A,f,floor);
    assert(std::abs(1/conductance-std::max(insulation/A,floor))<1e-12);
    const double film=ihm_thermal::film_resistance(A,4.,floor,1e12);
    assert(std::isfinite(10/(1/conductance+film)));
    assert(std::abs(film-1/(A*4))<1e-12);
  }
  // Subdividing an unequal region cannot create more clothing conductance.
  const double parent=1/ihm_thermal::segment_resistance(I,A,.36,floor);
  const double children=1/ihm_thermal::segment_resistance(I,A,.1,floor)+1/ihm_thermal::segment_resistance(I,A,.26,floor);
  assert(std::abs(parent-children)<1e-12);
  // Two distinct environmental temperatures: each boundary's opposite heat
  // transfer is counted once in the common algebraic clothing-node balance.
  const double rc=I/A, rv=ihm_thermal::film_resistance(A,4.,floor,1e12), rr=ihm_thermal::film_resistance(A,5.,floor,1e12);
  const double skin=34,air=22,radiant=20;
  const double cloth=(skin/rc+air/rv+radiant/rr)/(1/rc+1/rv+1/rr);
  const double input=(skin-cloth)/rc, conv=(cloth-air)/rv, rad=(cloth-radiant)/rr;
  assert(std::abs(input-conv-rad)<1e-10);
  bool rejected=false;
  try { ihm_thermal::segment_resistance(I,A,0,floor); } catch(const std::exception&) { rejected=true; }
  assert(rejected);
  std::cout << "PASS native unequal-area parallel/series insulation, finite nude boundary, and paired heat balance\\n";
}
''')
    subprocess.run(['c++','-std=c++17','-O2',str(source),'-o',str(out/'test')],check=True)
    result=subprocess.check_output([str(out/'test')],text=True)
    receipt={'status':'passed','output':result.strip(),'header_sha256':hashlib.sha256(header.read_bytes()).hexdigest(),'executable_sha256':hashlib.sha256((out/'test').read_bytes()).hexdigest()}
    (out/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(dict(receipt,output_dir=str(out)),indent=2))
    return out
def verify_runtime(directory,nude=None):
    from scripts.audit_long_horizon_thermal import state_audit
    from ihm.assembly.systemic_evidence import resolve_sources
    reports={}
    for name in ('source_boundary','corrected'):
        case=Path(directory)/name;execution=json.loads((case/'execution.json').read_text());assert execution['returncode']==0
        resolve_sources(BASE,case,execution['source_hashes'])
        initial=state_audit(case/'states/native_stabilized.xml');final=state_audit(case/'states/native_final.xml')
        rows=list(csv.DictReader((case/'native_multisystem.csv').open()));t=np.array([float(r['#Time(s)']) for r in rows]);assert t[-1]==3600
        assert np.isfinite(np.array([[float(v) for v in row.values()] for row in rows])).all()
        skin_source=sum(v for p,v in initial['metabolic_paths_w'].items() if p!='GroundToInternalCore')
        heat=np.array([float(r['TotalMetabolicRate(W)'])+skin_source-sum(float(r[k]) for k in ('ConvectiveHeatLoss(W)','RadiativeHeatLoss(W)','EvaporativeHeatLoss(W)','RespirationHeatLoss(W)')) for r in rows])
        integrated=float(np.trapezoid(np.r_[initial['storage_rate_w'],heat],np.r_[0,t]));change=final['stored_heat_j']-initial['stored_heat_j']
        error=change-integrated
        assert abs(error)<max(5.,.002*abs(change)), 'Sampled energy integral inconsistent with stored heat'
        if name=='corrected':
            assert abs(final['equal_skin_temperature_equivalent_clo']-.5)<1e-10
            for label,path in [('Convective','ClothingToEnvironment'),('Radiative','ClothingToEnclosure')]:
                coefficient=float(final['film_coefficients'][label+'HeatTranferCoefficient']['value']);area=float(final['environment']['SkinSurfaceArea']['value'])
                assert abs(final['resistance_k_w'][path]*coefficient*area-1)<1e-5
            assert final['core_temperature_c']>36 and final['core_temperature_c']>initial['core_temperature_c']-1
        reports[name]={'initial':initial,'final':final,'integrated_sampled_heat_j':integrated,'stored_heat_change_j':change,
                       'sampled_energy_integral_error_j':error,'native_csv_sha256':hashlib.sha256((case/'native_multisystem.csv').read_bytes()).hexdigest()}
    assert reports['corrected']['final']['core_temperature_c']>reports['source_boundary']['final']['core_temperature_c']+3
    if nude:
        nude=Path(nude);assert json.loads((nude/'execution.json').read_text())['returncode']==0
        final=state_audit(nude/'states/native_final.xml');assert final['reported_clo']==0
        assert abs(final['clothing_equivalent_resistance_k_w']-1e-8)<1e-15
        assert np.isfinite(final['boundary_loss_w']) and abs(final['circuit_balance_residual_w'])<1e-3
        reports['nude']=final
    out=Path(tempfile.mkdtemp(prefix='thermal-boundary-runtime-verification-',dir=BASE/'data/derived/audits'))
    report={'status':'passed','cases':reports,'interpretation':'Source-law correction and bounded one-hour execution verified; thermal equilibrium, long-horizon physiology and human agreement require separate checks.'}
    (out/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'status':'passed','output_dir':str(out),'final_core_c':{n:r['final']['core_temperature_c'] for n,r in reports.items() if 'final' in r}},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--runtime',type=Path);p.add_argument('--nude',type=Path);args=p.parse_args();verify_algebra()
    if args.runtime:verify_runtime(args.runtime,args.nude)
