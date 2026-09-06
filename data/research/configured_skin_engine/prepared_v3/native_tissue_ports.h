// Native-owned tissue exchange observations; no extra compartments or balances.
#pragma once
#include "native_body_ports.h"
#include <biogears/cdm/compartment/tissue/SETissueCompartment.h>
#include <biogears/cdm/circuit/fluid/SEFluidCircuitNode.h>
#include <biogears/engine/BioGearsPhysiologyEngine.h>
inline std::map<std::string,double> native_tissue_ports(BioGearsEngine& bg) {
  std::map<std::string,double> values;
  const double missing=std::numeric_limits<double>::quiet_NaN();
  auto& circuit=bg.GetCircuits().GetActiveCardiovascularCircuit();
  const auto liquid=[&](const std::string& key,const std::string& name){
    const auto* c=bg.GetCompartments().GetLiquidCompartment(name);
    values[key+".volume_ml"]=c&&c->HasVolume()?c->GetVolume(VolumeUnit::mL):missing;
    values[key+".pressure_mmhg"]=c&&c->HasPressure()?c->GetPressure(PressureUnit::mmHg):missing;
    for(const std::string sub : {"Albumin","Glucose","Oxygen","CarbonDioxide","Sodium","Potassium","Chloride"}) {
      const auto* substance=bg.GetSubstanceManager().GetSubstance(sub);
      const auto* q=c&&substance?c->GetSubstanceQuantity(*substance):nullptr;
      values[key+"."+sub+".mass_g"]=q&&q->HasMass()?q->GetMass(MassUnit::g):missing;
      values[key+"."+sub+".concentration_g_per_l"]=q&&q->HasConcentration()?q->GetConcentration(MassPerVolumeUnit::g_Per_L):missing;
      if(sub=="Sodium"||sub=="Potassium"||sub=="Chloride")
        values[key+"."+sub+".molarity_mmol_per_l"]=q&&q->HasMolarity()?q->GetMolarity(AmountPerVolumeUnit::mmol_Per_L):missing;
      if(sub=="Oxygen"||sub=="CarbonDioxide")
        values[key+"."+sub+".partial_pressure_mmhg"]=q&&q->HasPartialPressure()?q->GetPartialPressure(PressureUnit::mmHg):missing;
    }
  };
  for(const std::string organ : {"Fat","Bone","Brain","Gut","LeftKidney","RightKidney","Liver","LeftLung","RightLung","Muscle","Myocardium","Skin","Spleen"}) {
    const std::string prefix="tissue."+organ;
    liquid(prefix+".vascular",organ+"Vasculature");
    liquid(prefix+".extracellular",organ+"TissueExtracellular");
    liquid(prefix+".intracellular",organ+"TissueIntracellular");
    for(const std::string node : {organ+"E3",organ+"I",organ+"L1",organ+"L2"}) {
      const auto* n=circuit.GetNode(node);values["tissue.node."+node+".pressure_mmhg"]=n&&n->HasPressure()?n->GetPressure(PressureUnit::mmHg):missing;
    }
    for(const std::string path : {organ+"VTo"+organ+"E1",organ+"E1To"+organ+"E2",organ+"E2To"+organ+"E3",organ+"E3To"+organ+"I",organ+"E3To"+organ+"L1",organ+"L1To"+organ+"L2",organ+"ToLymphValve",organ+"E3ToGround"}) {
      const auto* p=circuit.GetPath(path);const std::string key="tissue.path."+path;
      values[key+".flow_ml_per_s"]=p&&p->HasFlow()?p->GetFlow(VolumePerTimeUnit::mL_Per_s):missing;
      values[key+".pressure_source_mmhg"]=p&&p->HasPressureSource()?p->GetPressureSource(PressureUnit::mmHg):missing;
      values[key+".compliance_ml_per_mmhg"]=p&&p->HasCompliance()?p->GetCompliance(FlowComplianceUnit::mL_Per_mmHg):missing;
    }
  }
  for(const std::string branch : {"LargeIntestine","SmallIntestine","Splanchnic"}) {
    const std::string name=branch+"VToGutE1";const auto* p=circuit.GetPath(name);
    values["tissue.path."+name+".pressure_source_mmhg"]=p&&p->HasPressureSource()?p->GetPressureSource(PressureUnit::mmHg):missing;
    values["tissue.path."+name+".flow_ml_per_s"]=p&&p->HasFlow()?p->GetFlow(VolumePerTimeUnit::mL_Per_s):missing;
  }
  liquid("tissue.Lymph","Lymph");liquid("tissue.VenaCava","VenaCava");
  for(const std::string name : {"LymphToVenaCava","SkinSweating"}) {
    const auto* p=circuit.GetPath(name);values["tissue.path."+name+".flow_ml_per_s"]=p&&p->HasFlow()?p->GetFlow(VolumePerTimeUnit::mL_Per_s):missing;
  }
  values["tissue.core_temperature_c"]=bg.GetEnergySystem()->GetCoreTemperature(TemperatureUnit::C);
  values["tissue.skin_temperature_c"]=bg.GetEnergySystem()->GetSkinTemperature(TemperatureUnit::C);
  return values;
}
