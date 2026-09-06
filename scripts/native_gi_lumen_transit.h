#pragma once
#include <cassert>
#include <biogears/cdm/compartment/fluid/SELiquidCompartment.h>
#include <biogears/cdm/compartment/substances/SELiquidSubstanceQuantity.h>
#include <biogears/cdm/properties/SEProperties.h>
#include <array>
#include <map>
#include <set>
#include <stdexcept>
#include <cmath>
namespace ihm_gi {
struct Receipt { std::array<double,3> water_ml; std::map<biogears::SESubstance*,std::array<double,3>> mass_g; double fecal_carrier_water_g; };
// Isolated unmapped-native-compartment writer. Mapped NextVolume integration
// requires an explicit circuit adapter; silently changing current volume is forbidden.
inline Receipt transit(std::array<biogears::SELiquidCompartment*,4> owners,
                       std::array<double,3> request_ml,double carrier_density_g_ml) {
 using namespace biogears;
 if(!std::isfinite(carrier_density_g_ml)||carrier_density_g_ml<=0)throw std::runtime_error("density required");
 std::set<SELiquidCompartment*> unique;
 std::set<SESubstance*> species;
 std::array<double,4> volume;
 for(size_t i=0;i<4;i++) {
  auto* c=owners[i];if(!c||!unique.insert(c).second||c->HasChildren()||c->HasNodeMapping())throw std::runtime_error("distinct unmapped leaf owners required");
  volume[i]=c->GetVolume(VolumeUnit::mL);if(!std::isfinite(volume[i])||volume[i]<0)throw std::runtime_error("invalid volume");
  for(auto* q:c->GetSubstanceQuantities())species.insert(&q->GetSubstance());
 }
 Receipt result{};
 for(size_t i=0;i<3;i++) {if(!std::isfinite(request_ml[i])||request_ml[i]<0)throw std::runtime_error("invalid request");result.water_ml[i]=std::min(volume[i],request_ml[i]);}
 std::map<SESubstance*,std::array<double,4>> final_mass;
 for(auto* s:species) {
  std::array<double,4> mass;std::array<double,3> moved;
  for(size_t i=0;i<4;i++) {auto* q=owners[i]->GetSubstanceQuantity(*s);if(!q||!q->HasMass())throw std::runtime_error("missing species inventory");mass[i]=q->GetMass(MassUnit::g);if(!std::isfinite(mass[i])||mass[i]<0)throw std::runtime_error("invalid mass");}
  for(size_t i=0;i<3;i++)moved[i]=volume[i]>0?mass[i]*(result.water_ml[i]/volume[i]):0.;
  auto next=mass;for(size_t i=0;i<3;i++){next[i]-=moved[i];next[i+1]+=moved[i];}
  final_mass[s]=next;result.mass_g[s]=moved;
 }
 auto next_volume=volume;for(size_t i=0;i<3;i++){next_volume[i]-=result.water_ml[i];next_volume[i+1]+=result.water_ml[i];}
 // All validation precedes mutation. Complete known species union is copied,
 // never a hardcoded nutrient list. Fourth owner is external cumulative output.
 for(size_t i=0;i<4;i++)owners[i]->GetVolume().SetValue(next_volume[i],VolumeUnit::mL);
 for(auto& item:final_mass)for(size_t i=0;i<4;i++){auto* q=owners[i]->GetSubstanceQuantity(*item.first);q->GetMass().SetValue(item.second[i],MassUnit::g);q->Balance(BalanceLiquidBy::Mass);}
 result.fecal_carrier_water_g=result.water_ml[2]*carrier_density_g_ml;
 return result;
}
}
