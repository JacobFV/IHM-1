"""Actual native evaporation helper: area, latent energy, capacity and wetness."""
from pathlib import Path
import json,subprocess,tempfile,sys
BASE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(BASE))
from scripts.audit_long_horizon_thermal import sha

def verify():
    header=BASE/'scripts/native_sweat_evaporation.h'
    assert header.is_file(),'Native latent-energy / area correction missing'
    out=Path(tempfile.mkdtemp(prefix='sweat-evaporation-algebra-',dir=BASE/'data/derived/audits'))
    cpp=out/'test.cpp';cpp.write_text('''#include "'''+str(header)+'''"
#include <cassert>
#include <iostream>
int main(){
  double fractions[]={.36,.07,.092,.092,.193,.193};
  for(double area:{.7,1.,1.9012784297057,3.}){
    for(double latent:{0.,10.,200.,2000.}){
      double sweat=0,diffuse=0,unmet=0;
      for(double f:fractions){
        auto r=ihm_thermal::evaporation_budget(latent,area,f,100.,.06);
        assert(r.wetted_fraction>=0 && r.wetted_fraction<=1);
        assert(r.total_w>=0 && r.total_w<=100*area*f+1e-10);
        assert(std::abs(r.sweat_evaporated_w+r.sweat_unevaporated_latent_w-latent*f)<1e-10);
        sweat+=r.sweat_evaporated_w;diffuse+=r.diffusion_w;unmet+=r.sweat_unevaporated_latent_w;
      }
      assert(std::abs(sweat+unmet-latent)<1e-9);
      assert(std::abs(sweat-std::min(latent,100*area))<1e-9);
      if(latent==0)assert(std::abs(diffuse-.06*100*area)<1e-10);
    }
  }
  auto whole=ihm_thermal::evaporation_budget(20.,1.9,.36,100.,.06);
  auto a=ihm_thermal::evaporation_budget(20.,1.9,.1,100.,.06);
  auto b=ihm_thermal::evaporation_budget(20.,1.9,.26,100.,.06);
  assert(std::abs(whole.total_w-a.total_w-b.total_w)<1e-12);
  for(double potential:{0.,-10.}){
    auto r=ihm_thermal::evaporation_budget(20.,1.9,.36,potential,.06);
    assert(r.total_w==0 && r.sweat_unevaporated_latent_w==20*.36);
  }
  bool rejected=false;try{ihm_thermal::evaporation_budget(20,0,.36,100,.06);}catch(const std::exception&){rejected=true;}
  assert(rejected);
  std::cout<<"PASS total sweat latent energy, unequal areas, subdivision, zero sweat, finite wetness and evaporative capacity\\n";
}
''')
    subprocess.run(['c++','-std=c++17','-O2',str(cpp),'-o',str(out/'test')],check=True)
    result=subprocess.check_output([str(out/'test')],text=True)
    report={'status':'passed','header_sha256':sha(header),'executable_sha256':sha(out/'test'),'output':result.strip()}
    (out/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(report,output_dir=str(out)),indent=2));return out
if __name__=='__main__':verify()
