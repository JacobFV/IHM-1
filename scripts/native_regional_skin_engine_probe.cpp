// Standalone full-engine acceptance; never used by the production session API.
#include <cassert>
#include <biogears/cdm/properties/SEProperties.h>
#include <biogears/engine/BioGearsPhysiologyEngine.h>
#include <biogears/engine/Controller/BioGearsEngine.h>
#include <biogears/engine/Controller/BioGearsConfiguration.h>
#include <biogears/cdm/patient/SEPatient.h>
#include "native_regional_skin.h"
#include "native_regional_species.h"
#include "native_body_ports.h"
#include "native_tissue_ports.h"
#include <memory>
#include <iostream>
#include <iomanip>
#include <map>

using namespace biogears;
class RegionalEngineProbe : public BioGearsEngine {
public:
  using BioGearsEngine::BioGearsEngine;
  std::unique_ptr<ihm::NativeRegionalSkin> regional;
  // Generated verbatim from held BioGearsEngine.cpp with exactly two hook calls.
#include "native_regional_engine_step.inc"
};

int main(int argc,char** argv) {
 try {
  if(argc!=3)return 2;
  const std::string mode=argv[1];
  if(mode!="baseline"&&mode!="zero"&&mode!="load")return 3;
  auto engine=std::make_unique<RegionalEngineProbe>("regional_engine.log");
  if(!engine->LoadState(argv[2]))throw std::runtime_error("Native state load failed");
  engine->SetAutoTrackFlag(false);
  auto& compartments=engine->BioGears::GetCompartments();
  auto& owner=*compartments.GetLiquidCompartment("SkinTissueExtracellular");
  std::map<std::string,double> initial_mass;
  for(auto* q:owner.GetSubstanceQuantities())if(q->HasMass())initial_mass[q->GetSubstance().GetName()]=q->GetMass(MassUnit::ug);
  const double initial_volume=owner.GetVolume(VolumeUnit::mL);
  if(mode!="baseline")engine->regional=std::make_unique<ihm::NativeRegionalSkin>(engine->GetCircuits().GetActiveCardiovascularCircuit(),compartments,owner);
  const double origin=engine->GetSimulationTime(TimeUnit::s),dt=engine->GetTimeStep(TimeUnit::s);
  if(std::abs(dt-.02)>1e-12)throw std::runtime_error("Unexpected native timestep");
  auto output=[&](int tick,double prior_volume,const std::array<double,3>& prior_regions){
    auto values=body_ports(*engine);auto tissue=native_tissue_ports(*engine);values.insert(tissue.begin(),tissue.end());
    values["acceptance.skin.initial_volume_ml"]=initial_volume;
    values["acceptance.skin.aggregate_volume_ml"]=owner.GetVolume(VolumeUnit::mL);
    if(tick>0){
      const auto& ledger=*ihm_regional::ihm_regional_sweat_ledger();
      const std::array<std::string,3> ions{"Sodium","Potassium","Chloride"};
      for(size_t i=0;i<3;++i){const auto prefix="acceptance.sweat."+ions[i]+".";
        values[prefix+"calls"]=ledger.calls[i];values[prefix+"requested_mg"]=ledger.requested_mg[i];
        values[prefix+"donor_before_mg"]=ledger.owner_before_mg[i];values[prefix+"donor_after_mg"]=ledger.owner_after_mg[i];
        values[prefix+"native_waste_before_mg"]=ledger.waste_before_mg[i];values[prefix+"native_waste_after_mg"]=ledger.waste_after_mg[i];
        values[prefix+"paired_mass_residual_mg"]=ledger.residual_mg[i];values[prefix+"charge_to_waste_mmol"]=ledger.charge_to_waste_mmol[i];
      }
    }
    double ownership_volume=owner.GetVolume(VolumeUnit::mL),initial_mass_error=0.;
    if(engine->regional) {
      ownership_volume=0;
      if(owner.HasNodeMapping()||!owner.HasChildren()||owner.GetLeaves().size()!=3)throw std::runtime_error("Regional aggregate ownership broken");
      const auto& leaves=compartments.GetLiquidLeafCompartments();
      if(std::find(leaves.begin(),leaves.end(),&owner)!=leaves.end())throw std::runtime_error("Aggregate present in native owning leaves");
      for(size_t i=0;i<3;++i){auto* c=engine->regional->children[i];ownership_volume+=c->GetVolume(VolumeUnit::mL);if(std::find(leaves.begin(),leaves.end(),c)==leaves.end())throw std::runtime_error("Regional native leaf cache stale");}
      for(auto* graph:compartments.GetLiquidGraphs())if(graph->GetCompartment(owner.GetName())==&owner)throw std::runtime_error("Aggregate remains a native transporter vertex");
    }
    values["acceptance.skin.volume_ownership_residual_ml"]=ownership_volume-owner.GetVolume(VolumeUnit::mL);
    for(auto* q:owner.GetSubstanceQuantities())if(q->HasMass()) {
      const auto name=q->GetSubstance().GetName();const double parent_mass=q->GetMass(MassUnit::ug);double mass=parent_mass;
      if(!std::isfinite(parent_mass)||parent_mass<0)throw std::runtime_error("Invalid parent species mass: "+name);
      if(engine->regional){mass=0;for(auto* c:engine->regional->children){double m=c->GetSubstanceQuantity(q->GetSubstance())->GetMass(MassUnit::ug);if(!std::isfinite(m)||m<0)throw std::runtime_error("Invalid regional species mass: "+name);mass+=m;}}
      values["acceptance.skin.species."+name+".mass_ug"]=parent_mass;
      values["acceptance.skin.species."+name+".ownership_residual_ug"]=mass-parent_mass;
      if(tick==0)initial_mass_error=std::max(initial_mass_error,std::abs(parent_mass-initial_mass.at(name)));
    }
    if(tick==0){values["acceptance.skin.install_mass_residual_ug"]=initial_mass_error;values["acceptance.skin.install_volume_residual_ml"]=owner.GetVolume(VolumeUnit::mL)-initial_volume;}
    auto flow=[&](const std::string& name){auto* p=engine->GetCircuits().GetFluidPath(name);if(!p||!p->HasFlow())throw std::runtime_error("Missing native flow ledger path: "+name);return p->GetFlow(VolumePerTimeUnit::mL_Per_s);};
    if(tick>0) {
      double incidence=flow("SkinE2ToSkinE3")-flow("SkinE3ToSkinI")-flow("SkinE3ToSkinL1")-flow("SkinSweating");
      values["acceptance.skin.fluid_incidence_ml_per_s"]=incidence;
      values["acceptance.skin.fluid_step_residual_ml"]=owner.GetVolume(VolumeUnit::mL)-prior_volume-dt*incidence;
    }
    if(engine->regional)for(size_t i=0;i<3;++i) {
      auto& s=*engine->regional;auto* c=s.children[i];const auto prefix="acceptance.region."+s.regions[i]+".";
      values[prefix+"volume_ml"]=c->GetVolume(VolumeUnit::mL);values[prefix+"pressure_mmhg"]=s.nodes[i].at("SkinE3")->GetPressure(PressureUnit::mmHg);
      values[prefix+"external_pressure_pa"]=s.drives[i]->GetPressureSource(PressureUnit::Pa);
      values[prefix+"lymph_flow_ml_per_s"]=s.paths[i].at("SkinE3ToSkinL1")->HasFlow()?s.paths[i].at("SkinE3ToSkinL1")->GetFlow(VolumePerTimeUnit::mL_Per_s):0.;
      for(auto* q:c->GetSubstanceQuantities())if(q->HasMass())values[prefix+q->GetSubstance().GetName()+".mass_ug"]=q->GetMass(MassUnit::ug);
      if(tick>0){auto q=[&](const std::string& path){return s.paths[i].at(path)->GetFlow(VolumePerTimeUnit::mL_Per_s);};double incidence=q("SkinE2ToSkinE3")-q("SkinE3ToSkinI")-q("SkinE3ToSkinL1")-q("SkinSweating");values[prefix+"fluid_step_residual_ml"]=c->GetVolume(VolumeUnit::mL)-prior_regions[i]-dt*incidence;}
    }
    std::cout<<std::setprecision(17)<<"RESULT {\"mode\":\""<<mode<<"\",\"tick\":"<<tick<<",\"time_s\":"<<engine->GetSimulationTime(TimeUnit::s)<<",\"elapsed_s\":"<<engine->GetSimulationTime(TimeUnit::s)-origin<<",\"values\":{";
    bool first=true;for(const auto& [key,value]:values){if(!first)std::cout<<',';first=false;std::cout<<'"'<<key<<"\":";if(std::isfinite(value))std::cout<<value;else std::cout<<"null";}std::cout<<"}}"<<std::endl;
  };
  output(0,initial_volume,{});
  for(int tick=1;tick<=12;++tick){
    const double before=owner.GetVolume(VolumeUnit::mL);std::array<double,3> prior{};
    if(engine->regional){for(size_t i=0;i<3;++i)prior[i]=engine->regional->children[i]->GetVolume(VolumeUnit::mL);engine->regional->set_pressure_pa(0,mode=="load"&&tick>=6&&tick<=7?133.322387415:0.);}
    if(!engine->AdvanceModelTime(false))throw std::runtime_error("Native engine step failed");
    output(tick,before,prior);
  }
  return 0;
 }catch(const std::exception& e){std::cerr<<"REGIONAL_ENGINE_FAILURE: "<<e.what()<<'\n';return 1;}
}
