#include "native_gi_mapped_transit.h"
#include <biogears/engine/Controller/BioGears.h>
#include <biogears/cdm/circuit/fluid/SEFluidCircuitCalculator.h>
#include <iostream>
int main(){
 using namespace biogears;
 BioGears bg("mapped.log");if(!bg.GetSubstances().LoadSubstanceDirectory())return 2;
 auto& mgr=bg.GetCompartments();auto& subs=bg.GetSubstances().GetSubstances();
 for(auto* s:subs)static_cast<SECompartmentManager&>(mgr).AddLiquidCompartmentSubstance(*s);
 auto& circuit=bg.GetCircuits().CreateFluidCircuit("IsolatedGIManualVolumeOwner");
 std::array<SELiquidCompartment*,4> c;std::array<SEFluidCircuitNode*,4> n;
 for(size_t i=0;i<4;i++){
  std::string name="GI_Mapped_"+std::to_string(i);n[i]=&circuit.CreateNode(name);c[i]=&mgr.CreateLiquidCompartment(name);
  n[i]->GetVolume().SetValue(i==0?10:0,VolumeUnit::mL);n[i]->GetNextVolume().SetValue(i==0?8:0,VolumeUnit::mL);
  n[i]->GetVolumeBaseline().SetValue(i==0?10:0,VolumeUnit::mL);c[i]->MapNode(*n[i]);
  n[i]->GetPressure().SetValue(0,PressureUnit::mmHg);n[i]->GetNextPressure().SetValue(0,PressureUnit::mmHg);
  for(size_t j=0;j<subs.size();j++)c[i]->GetSubstanceQuantity(*subs[j])->GetMass().SetValue(i==0?(j+1)*1e-8:0,MassUnit::g);
 }
 // Prospective SI volume8 differs from current10, as after upstream absorption.
 // Solute masses above are the corresponding post-absorption snapshot.
 SEFluidCircuitCalculator calculator(bg.GetLogger());ihm_gi_mapped::State state;
 auto r=state.transfer(c,n,{4,8,8},0);assert(r.water_ml[0]==4);
 assert(n[0]->GetVolume(VolumeUnit::mL)==10);assert(n[0]->GetNextVolume(VolumeUnit::mL)==4);
 assert(std::abs(c[0]->GetSubstanceQuantity(*subs[0])->GetMass(MassUnit::g)-.5e-8)<1e-20);
 bool rejected=false;try{state.save(c,n);}catch(...){rejected=true;}assert(rejected);
 calculator.PostProcess(circuit);state.settled(c,n);
 assert(n[0]->GetVolume(VolumeUnit::mL)==4);
 const auto saved=state.save(c,n);
 rejected=false;try{state.transfer(c,n,{0,0,0},0);}catch(...){rejected=true;}assert(rejected);
 state.transfer(c,n,{0,0,0},1);calculator.PostProcess(circuit);state.settled(c,n);
 assert(n[0]->GetVolume(VolumeUnit::mL)==4); // release/no-op has no stale flow
 state.load(saved,c,n);assert(state.consumed_epoch==0);assert(state.save(c,n)==saved);
 rejected=false;try{state.load(saved+"junk",c,n);}catch(...){rejected=true;}assert(rejected);assert(state.save(c,n)==saved);
 for(int epoch=1;epoch<5;epoch++){
  state.transfer(c,n,{100,100,100},epoch);calculator.PostProcess(circuit);state.settled(c,n);
  double v=0;for(auto* node:n)v+=node->GetVolume(VolumeUnit::mL);assert(std::abs(v-8)<1e-12);
  for(size_t j=0;j<subs.size();j++){double sum=0;for(auto* owner:c){double mass=owner->GetSubstanceQuantity(*subs[j])->GetMass(MassUnit::g);assert(mass>=0);sum+=mass;}assert(std::abs(sum-(j+1)*1e-8)<1e-18);}
 }
 assert(n[3]->GetVolume(VolumeUnit::mL)==8);
 std::cout<<"PASS mapped_species="<<subs.size()<<" fecal_ml=8 current_write_owner=native_PostProcess epoch="<<state.consumed_epoch<<" save_bytes="<<state.save(c,n).size()<<"\n";
}
