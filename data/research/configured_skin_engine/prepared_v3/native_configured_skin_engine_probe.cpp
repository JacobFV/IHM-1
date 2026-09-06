// Standalone full-engine acceptance; never used by the production session API.
#include <cassert>
#include <biogears/cdm/properties/SEProperties.h>
#include <biogears/engine/BioGearsPhysiologyEngine.h>
#include <biogears/engine/Controller/BioGearsEngine.h>
#include <biogears/engine/Controller/BioGearsConfiguration.h>
#include <biogears/cdm/patient/SEPatient.h>
#include "native_configured_regional_skin.h"
#include "native_regional_species.h"
#include "native_body_ports.h"
#include "native_tissue_ports.h"
#include <memory>
#include <iostream>
#include <iomanip>
#include <map>

using namespace biogears;
class ConfiguredRegionalEngineProbe : public BioGearsEngine {
public:
  using BioGearsEngine::BioGearsEngine;
  std::unique_ptr<ihm::NativeConfiguredRegionalSkin> regional;
  // Generated verbatim from held BioGearsEngine.cpp with exactly two hook calls.
#include "native_regional_engine_step.inc"
};

int main(int argc,char** argv) {
 try {
  if(argc!=3)return 2;
  const std::string mode=argv[1];
  if(mode!="baseline"&&mode!="zero"&&mode!="load")return 3;
  auto engine=std::make_unique<ConfiguredRegionalEngineProbe>("regional_engine.log");
  if(!engine->LoadState(argv[2]))throw std::runtime_error("Native state load failed");
  engine->SetAutoTrackFlag(false);
  auto& compartments=engine->BioGears::GetCompartments();
  auto& owner=*compartments.GetLiquidCompartment("SkinTissueExtracellular");
  std::map<std::string,double> initial_mass;
  for(auto* q:owner.GetSubstanceQuantities())if(q->HasMass())initial_mass[q->GetSubstance().GetName()]=q->GetMass(MassUnit::ug);
  const double initial_volume=owner.GetVolume(VolumeUnit::mL);
  if(mode!="baseline")engine->regional=std::make_unique<ihm::NativeConfiguredRegionalSkin>(engine->GetCircuits().GetActiveCardiovascularCircuit(),compartments,owner);
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
    // Observe the native owning-leaf inventory without introducing a second store.
    std::map<std::string,double> all_liquid_mass;
    double all_liquid_volume=0.;size_t unknown_mass=0,external_boundaries=0;
    const auto& owning_leaves=compartments.GetLiquidLeafCompartments();
    for(auto* leaf:owning_leaves) {
      if(leaf->HasChildren())throw std::runtime_error("Aggregate in native leaf inventory");
      const bool ambient=leaf->GetName()=="Ambient";
      if(ambient) {
        const auto& mapped=leaf->GetNodeMapping().GetNodes();
        if(mapped.size()!=1||mapped[0]->GetName()!="Ambient"||!leaf->HasVolume()||!std::isinf(leaf->GetVolume(VolumeUnit::mL))||leaf->GetVolume(VolumeUnit::mL)<0)
          throw std::runtime_error("Declared native Ambient boundary identity changed");
        ++external_boundaries;values["acceptance.boundary.Ambient.volume_is_positive_infinity"]=1;
        size_t unknown_boundary_mass=0;
        for(auto* q:leaf->GetSubstanceQuantities()) {
          const auto prefix="acceptance.boundary.Ambient.species."+q->GetSubstance().GetName()+".";
          if(!q->HasMass()){++unknown_boundary_mass;values[prefix+"mass_unknown"]=1;continue;}
          const double mass=q->GetMass(MassUnit::ug);
          if(std::isinf(mass)&&mass>0)values[prefix+"mass_is_positive_infinity"]=1;
          else if(std::isfinite(mass)&&mass>=0)values[prefix+"mass_ug"]=mass;
          else throw std::runtime_error("Invalid declared Ambient species mass");
        }
        values["acceptance.boundary.Ambient.unknown_mass_entries"]=unknown_boundary_mass;
        continue; // Named external reservoir is explicitly observed above, not an internal finite store.
      }
      if(leaf->HasVolume()){const double v=leaf->GetVolume(VolumeUnit::mL);if(!std::isfinite(v)||v<0)throw std::runtime_error("Invalid native owning liquid volume: "+leaf->GetName()+"="+std::to_string(v));all_liquid_volume+=v;}
      for(auto* q:leaf->GetSubstanceQuantities()) {
        if(!q->HasMass()){++unknown_mass;continue;}
        const double mass=q->GetMass(MassUnit::ug);
        if(!std::isfinite(mass)||mass<0)throw std::runtime_error("Invalid native owning species: "+leaf->GetName()+" / "+q->GetSubstance().GetName());
        all_liquid_mass[q->GetSubstance().GetName()]+=mass;
      }
    }
    values["acceptance.liquid_owners.external_boundary_count"]=external_boundaries;
    for(auto* link:compartments.GetActiveAerosolGraph().GetLinks()) {
      const bool source_ambient=link->GetSourceCompartment().GetName()=="Ambient";
      const bool target_ambient=link->GetTargetCompartment().GetName()=="Ambient";
      if(!source_ambient&&!target_ambient)continue;
      const auto key="acceptance.boundary.Ambient.active_aerosol_link."+link->GetName()+".";
      values[key+"source_is_Ambient"]=source_ambient?1:0;
      if(link->HasFlow()) {
        const double q=link->GetFlow(VolumePerTimeUnit::mL_Per_s);
        if(!std::isfinite(q))throw std::runtime_error("Invalid Ambient carrier flux");
        values[key+"carrier_flow_ml_per_s"]=q;
      }else values[key+"carrier_flow_unknown"]=1;
      values[key+"solute_flux_observed"]=0; // Do not fabricate mass flux from post-step concentrations.
    }
    values["acceptance.liquid_owners.count"]=owning_leaves.size();
    values["acceptance.liquid_owners.volume_ml"]=all_liquid_volume;
    values["acceptance.liquid_owners.unknown_mass_entries"]=unknown_mass;
    for(const auto& [name,mass]:all_liquid_mass)values["acceptance.liquid_owners.species."+name+".mass_ug"]=mass;
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
    for(const auto& name:ihm::NativeConfiguredRegionalSkin::path_names) {
      auto* path=engine->GetCircuits().GetFluidPath(name);
      if(!path)throw std::runtime_error("Missing original native law cache");
      const auto prefix="acceptance.original_path."+name+".";
      if(path->HasResistance())values[prefix+"resistance_mmhg_s_per_ml"]=path->GetResistance(FlowResistanceUnit::mmHg_s_Per_mL);
      if(path->HasCompliance())values[prefix+"compliance_ml_per_mmhg"]=path->GetCompliance(FlowComplianceUnit::mL_Per_mmHg);
      if(path->HasFlowSource())values[prefix+"flow_source_ml_per_s"]=path->GetFlowSource(VolumePerTimeUnit::mL_Per_s);
    }
    if(engine->regional)for(size_t i=0;i<3;++i) {
      auto& s=*engine->regional;auto* c=s.children[i];const auto prefix="acceptance.region."+s.regions[i]+".";
      values[prefix+"volume_ml"]=c->GetVolume(VolumeUnit::mL);values[prefix+"pressure_mmhg"]=s.nodes[i].at("SkinE3")->GetPressure(PressureUnit::mmHg);
      values[prefix+"external_pressure_pa"]=s.drives[i]->GetPressureSource(PressureUnit::Pa);
      values[prefix+"initial_fraction"]=s.fractions[i];
      for(const auto& [name,path]:s.paths[i]) {
        if(path->HasResistance())values[prefix+"path."+name+".resistance_mmhg_s_per_ml"]=path->GetResistance(FlowResistanceUnit::mmHg_s_Per_mL);
        if(path->HasCompliance())values[prefix+"path."+name+".compliance_ml_per_mmhg"]=path->GetCompliance(FlowComplianceUnit::mL_Per_mmHg);
        if(path->HasFlowSource())values[prefix+"path."+name+".flow_source_ml_per_s"]=path->GetFlowSource(VolumePerTimeUnit::mL_Per_s);
      }
      values[prefix+"lymph_flow_ml_per_s"]=s.paths[i].at("SkinE3ToSkinL1")->HasFlow()?s.paths[i].at("SkinE3ToSkinL1")->GetFlow(VolumePerTimeUnit::mL_Per_s):0.;
      for(auto* q:c->GetSubstanceQuantities())if(q->HasMass())values[prefix+q->GetSubstance().GetName()+".mass_ug"]=q->GetMass(MassUnit::ug);
      if(tick>0){auto q=[&](const std::string& path){return s.paths[i].at(path)->GetFlow(VolumePerTimeUnit::mL_Per_s);};double incidence=q("SkinE2ToSkinE3")-q("SkinE3ToSkinI")-q("SkinE3ToSkinL1")-q("SkinSweating");values[prefix+"fluid_step_residual_ml"]=c->GetVolume(VolumeUnit::mL)-prior_regions[i]-dt*incidence;}
    }
    std::cout<<std::setprecision(17)<<"RESULT {\"mode\":\""<<mode<<"\",\"configuration_sha256\":\""<<ihm::NativeConfiguredRegionalSkin::configuration_sha256<<"\",\"tick\":"<<tick<<",\"time_s\":"<<engine->GetSimulationTime(TimeUnit::s)<<",\"elapsed_s\":"<<engine->GetSimulationTime(TimeUnit::s)-origin<<",\"values\":{";
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
