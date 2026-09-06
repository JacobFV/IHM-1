// Staged integration acceptance helper. Synthetic mechanical/action challenge, NOT a
// clinical burn simulation. Requires paired regenerated schema/core/CDM variant.
#pragma once
#include "native_tissue_burn_state.h"
#include <biogears/engine/Controller/BioGearsEngine.h>
#include <biogears/engine/Systems/Tissue.h>
#include <biogears/engine/Controller/BioGearsCompartments.h>
#include <biogears/engine/Controller/BioGearsCircuits.h>
#include <biogears/cdm/properties/SEProperties.h>
#include <biogears/cdm/patient/actions/SEBurnWound.h>
#include <biogears/cdm/patient/actions/SEEscharotomy.h>
#include <biogears/cdm/scenario/SEActionManager.h>
#include <biogears/cdm/scenario/SEPatientActionCollection.h>
#include <biogears/cdm/circuit/fluid/SEFluidCircuit.h>
#include "io/biogears/BioGearsPhysiology.h"
#include <iostream>
#include <memory>
namespace ihm_burn_acceptance {
using namespace biogears;
// Access control only: the pointer remains a Tissue member pointer and is
// invoked on the actual engine-owned Tissue, without casting its object layout.
struct BurnAccess:Tissue {using Tissue::CalculateCompartmentalBurn;};
inline void require(bool ok,const char* reason){if(!ok)throw std::runtime_error(reason);}
inline Tissue& tissue(BioGearsEngine& e){return dynamic_cast<Tissue&>(e.BioGears::GetTissue());}
inline ihm_burn::State history(BioGearsEngine& e){CDM::BioGearsTissueSystemData data;io::BiogearsPhysiology::Marshall(tissue(e),data);return ihm_burn::decode(data.IHMBurnHistory().get());}
inline void calculate(BioGearsEngine& e){auto method=&BurnAccess::CalculateCompartmentalBurn;(tissue(e).*method)();}
inline void challenge(BioGearsEngine& e,const std::vector<std::string>& regions,bool escharotomy){
 e.GetActions().GetPatientActions().RemoveBurnWound();
 SEBurnWound burn;burn.SetDegreeOfBurn(SEBurnDegree::Third);
 burn.RemoveCompartments();
 for(const auto& r:regions)burn.AddCompartment(r);
 burn.SetTotalBodySurfaceArea(.5); // Native setter caps TBSA against selected regions.
 require(e.ProcessAction(burn),"Burn fixture action rejected");
 auto& c=e.GetCircuits().GetActiveCardiovascularCircuit();
 c.GetNode("MuscleE3")->GetPressure().SetValue(20,PressureUnit::mmHg);
 for(const auto& r:regions){
  auto name=r=="Trunk"?std::string("Muscle"):r;
  c.GetNode(name+"1")->GetPressure().SetValue(20,PressureUnit::mmHg);
  auto* p=c.GetPath("Aorta1To"+name+"1");
  p->GetFlow().SetValue(1,VolumePerTimeUnit::mL_Per_s);
  p->GetNextResistance().SetValue(100,FlowResistanceUnit::mmHg_s_Per_mL);
 }
 if(escharotomy){require(regions.size()==1,"One escharotomy region required");SEEscharotomy action;action.SetLocation(regions.front());require(e.ProcessAction(action),"Escharotomy fixture action rejected");}
 calculate(e);
}
inline void require_equal(BioGearsEngine& direct,BioGearsEngine& resumed){
 require(ihm_burn::encode(history(direct))==ihm_burn::encode(history(resumed)),"All12 histories differ after continuation");
}
inline void require_accumulated(const ihm_burn::State& before,const ihm_burn::State& after){
 for(size_t i=0;i<5;++i)require(after.resistance_mmHg_s_per_mL[i]>before.resistance_mmHg_s_per_mL[i],"Resistance did not accumulate");
 require(after.syndrome_count>=before.syndrome_count+5,"Syndrome count did not evolve");
 require(after.baseline_ecf_mL>0,"Burn baseline did not initialize");
 if(before.baseline_ecf_mL>0)require(after.baseline_ecf_mL==before.baseline_ecf_mL,"Captured burn baseline changed");
}
inline void require_released(const ihm_burn::State& before,const ihm_burn::State& after,size_t region){
 require(region<5,"Invalid escharotomy region");
 require(after.escharotomy_flags==(before.escharotomy_flags|(1u<<region)),"Native escharotomy flags differ");
 require(after.resistance_mmHg_s_per_mL[region]==0,"Native escharotomy did not release history");
 for(size_t i=0;i<5;++i)if(i!=region)require(after.resistance_mmHg_s_per_mL[i]==before.resistance_mmHg_s_per_mL[i],"Unchallenged regional history changed");
 require(after.baseline_ecf_mL==before.baseline_ecf_mL,"Escharotomy changed burn baseline");
}
} // namespace ihm_burn_acceptance
