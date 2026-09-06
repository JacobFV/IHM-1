// Regional native lifecycle, with the held engine method generated at build time.
#pragma once
#include <biogears/engine/Controller/BioGearsEngine.h>
#include <biogears/engine/Controller/BioGearsConfiguration.h>
#include <biogears/cdm/patient/SEPatient.h>
#include "native_body_ports.h"
#include "native_regional_skin.h"
#include "native_regional_species.h"
#include <memory>

class RegionalCoupledBioGearsEngine : public biogears::BioGearsEngine {
public:
  using BioGearsEngine::BioGearsEngine;
  std::unique_ptr<ihm::NativeRegionalSkin> regional;
  double external_pressure_pa=0.,generated_driver_pa=0.,applied_driver_pa=0.,external_work_j=0.;
  double prior_lung_m3=0.,initial_skin_volume_ml=0.,install_volume_residual_ml=0.;
  std::array<double,3> prior_region_ml{},fluid_step_residual_ml{},region_external_work_j{};
  std::map<std::string,double> install_mass_residual_ug;
  unsigned regional_steps=0;

  void install_regions() {
    if(regional)throw std::runtime_error("Regional topology already installed");
    auto& compartments=BioGears::GetCompartments();
    auto* owner=compartments.GetLiquidCompartment("SkinTissueExtracellular");
    if(!owner)throw std::runtime_error("Missing native Skin extracellular owner");
    initial_skin_volume_ml=owner->GetVolume(biogears::VolumeUnit::mL);
    std::map<std::string,double> before;
    for(auto* q:owner->GetSubstanceQuantities())if(q->HasMass())before[q->GetSubstance().GetName()]=q->GetMass(biogears::MassUnit::ug);
    regional=std::make_unique<ihm::NativeRegionalSkin>(GetCircuits().GetActiveCardiovascularCircuit(),compartments,*owner);
    install_volume_residual_ml=owner->GetVolume(biogears::VolumeUnit::mL)-initial_skin_volume_ml;
    for(auto* q:owner->GetSubstanceQuantities())if(q->HasMass())install_mass_residual_ug[q->GetSubstance().GetName()]=q->GetMass(biogears::MassUnit::ug)-before.at(q->GetSubstance().GetName());
  }
  void begin_coupled_step() {
    if(!regional)throw std::runtime_error("Regional topology is required before stepping");
    prior_lung_m3=GetRespiratorySystem()->GetTotalLungVolume(biogears::VolumeUnit::mL)*1e-6;
    for(size_t i=0;i<3;++i)prior_region_ml[i]=regional->children[i]->GetVolume(biogears::VolumeUnit::mL);
  }
  void coupled_after_preprocess() {
    auto* driver=GetCircuits().GetActiveRespiratoryCircuit().GetPath("EnvironmentToRespiratoryMuscle");
    generated_driver_pa=driver->GetNextPressureSource(biogears::PressureUnit::Pa);
    applied_driver_pa=generated_driver_pa+external_pressure_pa;
    if(external_pressure_pa!=0.)driver->GetNextPressureSource().SetValue(applied_driver_pa,biogears::PressureUnit::Pa);
    regional->after_preprocess();
  }
  void coupled_after_postprocess() {
    using namespace biogears;
    regional->after_postprocess();
    external_work_j-=external_pressure_pa*(GetRespiratorySystem()->GetTotalLungVolume(VolumeUnit::mL)*1e-6-prior_lung_m3);
    const double dt=GetTimeStep(TimeUnit::s);
    for(size_t i=0;i<3;++i) {
      auto flow=[&](const std::string& name){return regional->paths[i].at(name)->GetFlow(VolumePerTimeUnit::mL_Per_s);};
      const double delta=regional->children[i]->GetVolume(VolumeUnit::mL)-prior_region_ml[i];
      fluid_step_residual_ml[i]=delta-dt*(flow("SkinE2ToSkinE3")-flow("SkinE3ToSkinI")-flow("SkinE3ToSkinL1")-flow("SkinSweating"));
      region_external_work_j[i]-=regional->drives[i]->GetPressureSource(PressureUnit::Pa)*delta*1e-6;
    }
    ++regional_steps;
  }
  std::map<std::string,double> regional_ports() {
    using namespace biogears;
    if(!regional)throw std::runtime_error("Regional observation before topology installation");
    auto& compartments=BioGears::GetCompartments();auto& owner=regional->parent;
    if(owner.HasNodeMapping()||!owner.HasChildren()||owner.GetLeaves().size()!=3)throw std::runtime_error("Regional aggregate ownership broken");
    const auto& leaves=compartments.GetLiquidLeafCompartments();
    if(std::find(leaves.begin(),leaves.end(),&owner)!=leaves.end())throw std::runtime_error("Aggregate present in owning leaves");
    for(auto* graph:compartments.GetLiquidGraphs())if(graph->GetCompartment(owner.GetName())==&owner)throw std::runtime_error("Aggregate remains transporter vertex");
    std::map<std::string,double> values;
    const std::string prefix="tissue.regional_skin.";
    double sum=0,fluid_residual=0;
    for(size_t i=0;i<3;++i) {
      auto* child=regional->children[i];
      if(std::find(leaves.begin(),leaves.end(),child)==leaves.end())throw std::runtime_error("Regional owning leaf missing");
      const auto key=prefix+regional->regions[i]+".";
      const double volume=child->GetVolume(VolumeUnit::mL);
      if(!std::isfinite(volume)||volume<0)throw std::runtime_error("Invalid regional volume");
      sum+=volume;fluid_residual+=fluid_step_residual_ml[i];
      values[key+"fraction"]=regional->fractions[i];values[key+"volume_ml"]=volume;
      values[key+"pressure_mmhg"]=regional->nodes[i].at("SkinE3")->GetPressure(PressureUnit::mmHg);
      values[key+"requested_pressure_pa"]=regional->requested_pa[i];
      values[key+"applied_pressure_pa"]=regional->drives[i]->GetPressureSource(PressureUnit::Pa);
      values[key+"external_work_j"]=region_external_work_j[i];
      if(regional_steps)values[key+"fluid_step_residual_ml"]=fluid_step_residual_ml[i];
      for(const auto& name:regional->path_names) {
        auto* path=regional->paths[i].at(name);
        if(path->HasFlow())values[key+"path."+name+".flow_ml_per_s"]=path->GetFlow(VolumePerTimeUnit::mL_Per_s);
      }
      for(auto* q:child->GetSubstanceQuantities())if(q->HasMass()) {
        const double mass=q->GetMass(MassUnit::ug);
        if(!std::isfinite(mass)||mass<0)throw std::runtime_error("Invalid regional species mass");
        values[key+"species."+q->GetSubstance().GetName()+".mass_ug"]=mass;
      }
    }
    values[prefix+"aggregate.volume_ml"]=owner.GetVolume(VolumeUnit::mL);
    values[prefix+"aggregate.volume_ownership_residual_ml"]=sum-owner.GetVolume(VolumeUnit::mL);
    values[prefix+"aggregate.install_volume_residual_ml"]=install_volume_residual_ml;
    values[prefix+"aggregate.initial_volume_ml"]=initial_skin_volume_ml;
    values[prefix+"completed_native_steps"]=regional_steps;
    if(regional_steps)values[prefix+"aggregate.fluid_step_residual_ml"]=fluid_residual;
    for(auto* q:owner.GetSubstanceQuantities())if(q->HasMass()) {
      const auto name=q->GetSubstance().GetName();const double mass=q->GetMass(MassUnit::ug);double owning_sum=0.;
      if(!std::isfinite(mass)||mass<0)throw std::runtime_error("Invalid aggregate species mass");
      for(auto* child:regional->children)owning_sum+=child->GetSubstanceQuantity(q->GetSubstance())->GetMass(MassUnit::ug);
      const auto key=prefix+"aggregate.species."+name+".";
      values[key+"mass_ug"]=mass;values[key+"ownership_residual_ug"]=owning_sum-mass;
      values[key+"install_residual_ug"]=install_mass_residual_ug.at(name);
    }
    if(regional_steps) {
      const auto& ledger=*ihm_regional::ihm_regional_sweat_ledger();
      const std::array<std::string,3> ions{"Sodium","Potassium","Chloride"};
      for(size_t i=0;i<3;++i) {
        const auto key=prefix+"sweat."+ions[i]+".";
        values[key+"calls"]=ledger.calls[i];values[key+"requested_mg"]=ledger.requested_mg[i];
        values[key+"owner_before_mg"]=ledger.owner_before_mg[i];values[key+"owner_after_mg"]=ledger.owner_after_mg[i];
        values[key+"waste_before_mg"]=ledger.waste_before_mg[i];values[key+"waste_after_mg"]=ledger.waste_after_mg[i];
        values[key+"paired_mass_residual_mg"]=ledger.residual_mg[i];values[key+"charge_to_waste_mmol"]=ledger.charge_to_waste_mmol[i];
      }
    }
    return values;
  }
  // Exact held implementation with three narrow lifecycle insertions.
#include "native_regional_coupled_step.inc"
};
