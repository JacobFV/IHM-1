// Public whole-engine acceptance helpers; synthetic inventory, not a dosing study.
#pragma once
#include <cmath>
#include <iomanip>
#include <map>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <biogears/engine/Controller/BioGearsEngine.h>
#include <biogears/engine/Controller/BioGearsSubstances.h>
#include <biogears/engine/Controller/BioGearsCompartments.h>
#include <biogears/cdm/system/physiology/SEGastrointestinalSystem.h>
#include <biogears/cdm/patient/actions/SESubstanceOralDose.h>
#include <biogears/cdm/compartment/fluid/SELiquidCompartment.h>
#include <biogears/cdm/compartment/substances/SELiquidSubstanceQuantity.h>
#include <biogears/cdm/properties/SEProperties.h>
namespace ihm_gi_acceptance {
using namespace biogears;
using Observation=std::map<std::string,double>;
inline void require(bool ok,const std::string& reason){if(!ok)throw std::runtime_error(reason);}
inline SESubstance& drug(BioGearsEngine& e){auto* s=e.BioGears::GetSubstances().GetSubstance("Fentanyl");require(s!=nullptr,"GI fixture requires native Fentanyl");return *s;}
inline SEDrugTransitState& transit(BioGearsEngine& e){auto* s=e.BioGears::GetGastrointestinal().GetDrugTransitState(&drug(e));require(s!=nullptr,"GI fixture requires initialized oral transit state");return *s;}
inline void request_tiny_oral_dose(BioGearsEngine& e){
 require(e.BioGears::GetGastrointestinal().GetDrugTransitStates().empty(),"GI fixture expects no preexisting oral inventory");
 SESubstanceOralDose dose(drug(e));dose.SetAdminRoute(SEOralAdministrationType::Gastrointestinal);
 dose.GetDose().SetValue(.001,MassUnit::ug);require(e.ProcessAction(dose),"Native oral action rejected");
 // Caller advances one real .02 s tick before seeding: Drugs PreProcess initializes CAT.
}
inline void seed_distinct_transit(BioGearsEngine& e){
 auto& s=transit(e);std::vector<double> solid(9),dissolved(9),enterocyte(8);
 // Distinct short decimal microgram inventories keep the challenge tiny. These
 // replace the initialized test inventory explicitly, not inferred dose history.
 for(size_t i=0;i<9;++i){solid[i]=(i+1)*.000001;dissolved[i]=(i+21)*.000001;}
 for(size_t i=0;i<8;++i)enterocyte[i]=(i+41)*.000001;
 require(s.SetLumenSolidMasses(solid,MassUnit::ug),"Solid setter failed");
 require(s.SetLumenDissolvedMasses(dissolved,MassUnit::ug),"Dissolved setter failed");
 require(s.SetEnterocyteMasses(enterocyte,MassUnit::ug),"Enterocyte setter failed");
 s.GetTotalMassMetabolized().SetValue(.000123,MassUnit::ug);
 s.GetTotalMassExcreted().SetValue(.000456,MassUnit::ug);
}
inline Observation capture(BioGearsEngine& e){
 Observation o;auto& s=transit(e);
 auto add=[&](const char* name,const std::vector<double>& v,size_t expected){require(v.size()==expected,"GI vector shape mismatch");for(size_t i=0;i<v.size();++i)o[std::string(name)+"/"+std::to_string(i)+"/ug"]=v[i];};
 add("solid",s.GetLumenSolidMasses(MassUnit::ug),9);
 add("dissolved",s.GetLumenDissolvedMasses(MassUnit::ug),9);
 add("enterocyte",s.GetEnterocyteMasses(MassUnit::ug),8);
 o["metabolized/ug"]=s.GetTotalMassMetabolized().GetValue(MassUnit::ug);
 o["excreted/ug"]=s.GetTotalMassExcreted().GetValue(MassUnit::ug);
 o["transit_state_count"]=e.BioGears::GetGastrointestinal().GetDrugTransitStates().size();
 auto* c=e.GetCompartments().GetLiquidCompartment("SmallIntestineVasculature");require(c!=nullptr,"Missing portal source compartment");
 auto* q=c->GetSubstanceQuantity(drug(e));require(q&&q->HasMass(),"Missing native portal drug mass");
 o["SmallIntestineVasculature/Fentanyl/ug"]=q->GetMass(MassUnit::ug);
 for(const auto& kv:o)require(std::isfinite(kv.second)&&kv.second>=0,"Invalid GI inventory: "+kv.first);
 return o;
}
inline bool equal(const Observation& a,const Observation& b,std::ostream& differences){
 bool ok=true;differences<<std::setprecision(17);
 for(const auto& kv:a){auto it=b.find(kv.first);if(it==b.end()){differences<<"GI missing "<<kv.first<<'\n';ok=false;}
 else if(it->second!=kv.second){differences<<"GI difference "<<kv.first<<" before="<<kv.second<<" after="<<it->second<<" delta="<<it->second-kv.second<<'\n';ok=false;}}
 for(const auto& kv:b)if(!a.count(kv.first)){differences<<"GI unexpected "<<kv.first<<'\n';ok=false;}
 return ok; // Exact comparison, no hidden tolerance or normalization.
}
inline void require_evolved(const Observation& before,const Observation& after){
 bool changed=false;for(const auto& kv:before)if(kv.first.find("solid/")==0||kv.first.find("dissolved/")==0||kv.first.find("enterocyte/")==0){auto it=after.find(kv.first);require(it!=after.end(),"Missing evolution key");changed|=it->second!=kv.second;}
 require(changed,"GI next tick did not evolve CAT inventory");
 require(after.at("metabolized/ug")>=before.at("metabolized/ug"),"GI metabolism cumulative decreased");
 require(after.at("excreted/ug")>=before.at("excreted/ug"),"GI excretion cumulative decreased");
}
} // namespace ihm_gi_acceptance
