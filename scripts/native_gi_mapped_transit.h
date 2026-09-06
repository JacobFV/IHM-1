#pragma once
#include <cassert>
#include <biogears/cdm/circuit/fluid/SEFluidCircuitNode.h>
#include <biogears/cdm/substance/SESubstance.h>
#include <sstream>
#include <iomanip>
#include <cstdint>
#include <biogears/cdm/compartment/fluid/SELiquidCompartment.h>
#include <biogears/cdm/compartment/substances/SELiquidSubstanceQuantity.h>
#include <biogears/cdm/properties/SEProperties.h>
#include <array>
#include <map>
#include <set>
#include <stdexcept>
#include <cmath>
namespace ihm_gi_mapped {
struct Receipt { std::array<double,3> water_ml; std::map<biogears::SESubstance*,std::array<double,3>> mass_g; double fecal_carrier_water_g; };
// Manual NextVolume owner; no transit flow paths may also move these volumes.
inline Receipt apply(std::array<biogears::SELiquidCompartment*,4> owners, std::array<biogears::SEFluidCircuitNode*,4> nodes,
                       std::array<double,3> request_ml,double carrier_density_g_ml) {
 using namespace biogears;
 if(!std::isfinite(carrier_density_g_ml)||carrier_density_g_ml<=0)throw std::runtime_error("density required");
 std::set<SELiquidCompartment*> unique;
 std::set<SESubstance*> species;
 std::array<double,4> volume;
 for(size_t i=0;i<4;i++) {
  auto* c=owners[i];if(!c||!unique.insert(c).second||c->HasChildren()||!c->HasNodeMapping())throw std::runtime_error("distinct mapped leaf owners required");
  auto& mapping=c->GetNodeMapping().GetNodes();if(mapping.size()!=1||mapping[0]!=nodes[i])throw std::runtime_error("binding mismatch");
  volume[i]=nodes[i]->GetNextVolume(VolumeUnit::mL);if(!std::isfinite(volume[i])||volume[i]<0)throw std::runtime_error("invalid volume");
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
 for(size_t i=0;i<4;i++)nodes[i]->GetNextVolume().SetValue(next_volume[i],VolumeUnit::mL);
 for(auto& item:final_mass)for(size_t i=0;i<4;i++){auto* q=owners[i]->GetSubstanceQuantity(*item.first);q->GetMass().SetValue(item.second[i],MassUnit::g);/* Balance only after native PostProcess commits prospective volume. */}
 result.fecal_carrier_water_g=result.water_ml[2]*carrier_density_g_ml;
 return result;
}
struct State {
 static void bindings(std::array<biogears::SELiquidCompartment*,4> owners,std::array<biogears::SEFluidCircuitNode*,4> nodes) {
  std::set<biogears::SELiquidCompartment*> cs;std::set<biogears::SEFluidCircuitNode*> ns;
  for(size_t i=0;i<4;i++)if(!owners[i]||!nodes[i]||!cs.insert(owners[i]).second||!ns.insert(nodes[i]).second||owners[i]->HasChildren()||!owners[i]->HasVolume()||owners[i]->GetNodeMapping().GetNodes()!=std::vector<biogears::SEFluidCircuitNode*>{nodes[i]})throw std::runtime_error("distinct exact mapped bindings required");
 }

 std::int64_t consumed_epoch=-1;bool pending=false;
 Receipt transfer(std::array<biogears::SELiquidCompartment*,4> owners,std::array<biogears::SEFluidCircuitNode*,4> nodes,std::array<double,3> requests,std::int64_t epoch) {
  bindings(owners,nodes);if(pending||epoch<0||epoch<=consumed_epoch)throw std::runtime_error("pending or repeated epoch");
  auto result=apply(owners,nodes,requests,1.);consumed_epoch=epoch;pending=true;return result;
 }
 void settled(std::array<biogears::SELiquidCompartment*,4> owners,std::array<biogears::SEFluidCircuitNode*,4> nodes) {
  using namespace biogears;bindings(owners,nodes);if(!pending)throw std::runtime_error("no pending transaction");
  for(auto* n:nodes)if(n->GetVolume(VolumeUnit::mL)!=n->GetNextVolume(VolumeUnit::mL))throw std::runtime_error("native PostProcess not committed");
  for(auto* c:owners)for(auto* q:c->GetSubstanceQuantities())q->Balance(BalanceLiquidBy::Mass);
  pending=false;
 }
 std::string save(std::array<biogears::SELiquidCompartment*,4> owners,std::array<biogears::SEFluidCircuitNode*,4> nodes)const {
  using namespace biogears;bindings(owners,nodes);if(pending)throw std::runtime_error("save only settled boundary");
  std::ostringstream out;out<<std::setprecision(17)<<"IHM_GI_TRANSIT_1 "<<consumed_epoch<<'\n';
  for(size_t i=0;i<4;i++){
   if(nodes[i]->GetVolume(VolumeUnit::mL)!=nodes[i]->GetNextVolume(VolumeUnit::mL))throw std::runtime_error("unsettled volume");
   out<<std::quoted(owners[i]->GetName())<<' '<<nodes[i]->GetVolume(VolumeUnit::mL)<<' '<<owners[i]->GetSubstanceQuantities().size()<<'\n';
   for(auto* q:owners[i]->GetSubstanceQuantities()){if(!q->HasMass()||!std::isfinite(q->GetMass(MassUnit::g))||q->GetMass(MassUnit::g)<0)throw std::runtime_error("invalid saved mass");out<<std::quoted(q->GetSubstance().GetName())<<' '<<q->GetMass(MassUnit::g)<<'\n';}
  }return out.str();
 }
 void load(const std::string& data,std::array<biogears::SELiquidCompartment*,4> owners,std::array<biogears::SEFluidCircuitNode*,4> nodes) {
  using namespace biogears;bindings(owners,nodes);if(pending)throw std::runtime_error("cannot load over pending work");
  std::istringstream in(data);std::string schema;std::int64_t epoch;in>>schema>>epoch;
  if(!in||schema!="IHM_GI_TRANSIT_1"||epoch< -1)throw std::runtime_error("invalid state header");
  std::array<double,4> volumes;std::array<std::vector<double>,4> masses;
  for(size_t i=0;i<4;i++){
   std::string name;size_t count;in>>std::quoted(name)>>volumes[i]>>count;
   if(!in||name!=owners[i]->GetName()||count!=owners[i]->GetSubstanceQuantities().size()||!std::isfinite(volumes[i])||volumes[i]<0)throw std::runtime_error("state owner/schema mismatch");
   for(auto* q:owners[i]->GetSubstanceQuantities()){double mass;in>>std::quoted(name)>>mass;if(!in||name!=q->GetSubstance().GetName()||!std::isfinite(mass)||mass<0)throw std::runtime_error("state species mismatch");masses[i].push_back(mass);}
  }
  in>>std::ws;if(!in.eof())throw std::runtime_error("trailing state data");
  for(size_t i=0;i<4;i++){nodes[i]->GetVolume().SetValue(volumes[i],VolumeUnit::mL);nodes[i]->GetNextVolume().SetValue(volumes[i],VolumeUnit::mL);size_t j=0;for(auto* q:owners[i]->GetSubstanceQuantities()){q->GetMass().SetValue(masses[i][j++],MassUnit::g);q->Balance(BalanceLiquidBy::Mass);}}
  consumed_epoch=epoch;
 }
};
}
