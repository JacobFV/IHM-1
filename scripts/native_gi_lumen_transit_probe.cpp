#include "native_gi_lumen_transit.h"
#include <biogears/engine/Controller/BioGears.h>
#include <iostream>
#include <cassert>
int main(){
 using namespace biogears;
 BioGears bg("probe.log");if(!bg.GetSubstances().LoadSubstanceDirectory())return 2;
 auto& mgr=bg.GetCompartments();auto& subs=bg.GetSubstances().GetSubstances();
 for(auto* s:subs)static_cast<SECompartmentManager&>(mgr).AddLiquidCompartmentSubstance(*s);
 std::array<SELiquidCompartment*,4> c={&mgr.CreateLiquidCompartment("FixtureExistingSI"),&mgr.CreateLiquidCompartment("FixtureNewColon"),&mgr.CreateLiquidCompartment("FixtureNewRectum"),&mgr.CreateLiquidCompartment("ExternalFecalLedger")};
 for(size_t i=0;i<4;i++){c[i]->GetVolume().SetValue(i==0?10.:0.,VolumeUnit::mL);for(size_t j=0;j<subs.size();j++){auto* q=c[i]->GetSubstanceQuantity(*subs[j]);q->GetMass().SetValue(i==0?(j+1)*1.e-8:0.,MassUnit::g);q->Balance(BalanceLiquidBy::Mass);}}
 auto r=ihm_gi::transit(c,{4,10,10},1.);
 assert(r.water_ml[0]==4&&r.water_ml[1]==0&&r.water_ml[2]==0);
 for(int step=0;step<4;step++) {
  ihm_gi::transit(c,{100,100,100},1.);
  double v=0;for(auto* x:c){assert(x->GetVolume(VolumeUnit::mL)>=0);v+=x->GetVolume(VolumeUnit::mL);}assert(std::abs(v-10)<1e-12);
  for(size_t j=0;j<subs.size();j++){double sum=0;for(auto* x:c){double m=x->GetSubstanceQuantity(*subs[j])->GetMass(MassUnit::g);assert(m>=0);sum+=m;}assert(std::abs(sum-(j+1)*1.e-8)<1e-18);}
 }
 assert(c[3]->GetVolume(VolumeUnit::mL)==10);
 // Dry unknown/dissolved residue is retained; no fabricated aqueous carrier.
 c[0]->GetSubstanceQuantity(*subs[0])->GetMass().SetValue(1.e-9,MassUnit::g);
 r=ihm_gi::transit(c,{1,1,1},1.);assert(r.mass_g[subs[0]][0]==0);
 bool rejected=false;auto alias=c;alias[1]=alias[0];try{ihm_gi::transit(alias,{1,1,1},1.);}catch(const std::runtime_error&){rejected=true;}assert(rejected);
 std::cout<<"PASS species="<<subs.size()<<" fecal_water_ml="<<c[3]->GetVolume(VolumeUnit::mL)<<" carrier_water_g=10 density_g_ml=1\n";
}
