#include "native_regional_species.h"
#include <biogears/engine/Controller/BioGears.h>
#include <iostream>
#include <iomanip>
using namespace biogears;
int main(){try{
 BioGears bg("regional_sweat_fixture.log");auto& mgr=bg.GetCompartments();
 auto& parent=mgr.CreateLiquidCompartment("ProbeSkin");std::array<SELiquidCompartment*,3> children{};
 for(size_t i=0;i<3;++i){children[i]=&mgr.CreateLiquidCompartment("ProbeSkin"+std::to_string(i));children[i]->GetVolume().SetValue(1,VolumeUnit::mL);parent.AddChild(*children[i]);}
 const std::array<std::string,3> names{"Sodium","Potassium","Chloride"};const std::array<int,3> valence{1,1,-1};
 for(size_t s=0;s<3;++s){auto* sub=bg.GetSubstances().GetSubstance(names[s]);static_cast<SECompartmentManager&>(mgr).AddLiquidCompartmentSubstance(*sub);for(size_t i=0;i<3;++i)children[i]->GetSubstanceQuantity(*sub)->GetMass().SetValue(i==0?1.:i==1?2.:7.,MassUnit::mg);}
 parent.StateChange();parent.Balance(BalanceLiquidBy::Mass);
 double residual=0.,charge_residual=0.;unsigned passed=0;
 for(size_t s=0;s<3;++s){auto* sub=bg.GetSubstances().GetSubstance(names[s]);auto& quantity=*parent.GetSubstanceQuantity(*sub);SEScalarMass waste;waste.SetValue(0,MassUnit::mg);
  ihm_regional::withdraw_skin_ion(parent,quantity,3.,waste,s,valence[s]);
  double mass=quantity.GetMass(MassUnit::mg),lost=waste.GetValue(MassUnit::mg);
  if(std::abs(mass-7.)>1e-12||lost!=3.||!quantity.GetMass().IsReadOnly())throw std::runtime_error("Sweat owning-leaf allocation failed");
  residual=std::max(residual,std::abs(mass+lost-10.));
  double charge=(mass+lost-10.)/sub->GetMolarMass(MassPerAmountUnit::mg_Per_mmol)*valence[s];charge_residual=std::max(charge_residual,std::abs(charge));++passed;
  ihm_regional::withdraw_skin_ion(parent,quantity,0.,waste,s,valence[s]);if(quantity.GetMass(MassUnit::mg)!=mass||waste.GetValue(MassUnit::mg)!=lost)throw std::runtime_error("Zero debit changed state");++passed;
  bool rejected=false;try{ihm_regional::withdraw_skin_ion(parent,quantity,8.,waste,s,valence[s]);}catch(const std::runtime_error&){rejected=true;}
  if(!rejected||quantity.GetMass(MassUnit::mg)!=mass||waste.GetValue(MassUnit::mg)!=lost)throw std::runtime_error("Overdraw not rejected before mutation");++passed;
  ihm_regional::withdraw_skin_ion(parent,quantity,mass,waste,s,valence[s]);if(quantity.GetMass(MassUnit::mg)!=0.||std::abs(waste.GetValue(MassUnit::mg)-10.)>1e-12)throw std::runtime_error("Exact donor depletion failed");++passed;
 }
 std::cout<<std::setprecision(17)<<"RESULT {\"passed\":true,\"checks\":"<<passed<<",\"maximum_mass_residual_mg\":"<<residual<<",\"maximum_paired_charge_residual_mmol\":"<<charge_residual<<",\"native_cdm_ownership\":true}\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
