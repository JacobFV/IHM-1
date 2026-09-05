// Observational ports. All physiological state remains owned by BioGears.
#pragma once
#include <biogears/cdm/system/physiology/SEBloodChemistrySystem.h>
#include <biogears/cdm/system/physiology/SECardiovascularSystem.h>
#include <biogears/cdm/system/physiology/SERespiratorySystem.h>
#include <biogears/cdm/system/physiology/SENervousSystem.h>
#include <biogears/cdm/system/physiology/SEGastrointestinalSystem.h>
#include <biogears/cdm/system/physiology/SEEndocrineSystem.h>
#include <biogears/cdm/system/physiology/SEHepaticSystem.h>
#include <biogears/cdm/system/physiology/SETissueSystem.h>
#include <biogears/cdm/system/physiology/SEEnergySystem.h>
#include <biogears/cdm/system/physiology/SERenalSystem.h>
#include <biogears/cdm/patient/SENutrition.h>
#include <biogears/cdm/compartment/SECompartmentManager.h>
#include <biogears/cdm/compartment/fluid/SELiquidCompartment.h>
#include <biogears/cdm/compartment/substances/SELiquidSubstanceQuantity.h>
#include <biogears/cdm/substance/SESubstanceManager.h>
#include <biogears/cdm/circuit/fluid/SEFluidCircuit.h>
#include <biogears/cdm/circuit/fluid/SEFluidCircuitPath.h>
#include <biogears/engine/Controller/BioGearsEngine.h>
#include <biogears/engine/Controller/BioGearsCircuits.h>
#include <map>
#include <limits>
using namespace biogears;
inline std::map<std::string,double> body_ports(BioGearsEngine& bg) {
  std::map<std::string,double> v;
  const auto* cv=bg.GetCardiovascularSystem();
  const auto* r=bg.GetRespiratorySystem();
  const auto* b=bg.GetBloodChemistrySystem();
  const auto* n=bg.GetNervousSystem();
  const auto* gi=bg.GetGastrointestinalSystem();
  const auto* t=bg.GetTissueSystem();
  const auto* e=bg.GetEnergySystem();
  const auto* k=bg.GetRenalSystem();
  v["heart_rate_per_min"]=cv->GetHeartRate(FrequencyUnit::Per_min);
  v["arterial_pressure_mmhg"]=cv->GetArterialPressure(PressureUnit::mmHg);
  v["mean_arterial_pressure_mmhg"]=cv->GetMeanArterialPressure(PressureUnit::mmHg);
  v["cardiac_output_l_per_min"]=cv->GetCardiacOutput(VolumePerTimeUnit::L_Per_min);
  v["blood_volume_ml"]=cv->GetBloodVolume(VolumeUnit::mL);
  v["lung_volume_ml"]=r->GetTotalLungVolume(VolumeUnit::mL);
  v["respiration_rate_per_min"]=r->GetRespirationRate(FrequencyUnit::Per_min);
  v["respiratory_request_cmh2o"]=r->GetRespirationDriverPressure(PressureUnit::cmH2O);
  v["respiratory_request_per_min"]=r->GetRespirationDriverFrequency(FrequencyUnit::Per_min);
  v["respiratory_cycle_fraction"]=r->GetRespirationCyclePercentComplete();
  auto& circuit=bg.GetCircuits().GetActiveRespiratoryCircuit();
  v["respiratory_applied_cmh2o"]=circuit.GetPath("EnvironmentToRespiratoryMuscle")->GetPressureSource(PressureUnit::cmH2O);
  v["airway_flow_l_per_s"]=circuit.GetPath("MouthToTrachea")->GetFlow(VolumePerTimeUnit::L_Per_s);
  v["oxygen_saturation"]=b->GetOxygenSaturation();
  v["arterial_o2_mmhg"]=b->GetArterialOxygenPressure(PressureUnit::mmHg);
  v["arterial_co2_mmhg"]=b->GetArterialCarbonDioxidePressure(PressureUnit::mmHg);
  v["arterial_ph"]=b->GetArterialBloodPH();
  v["nervous_heart_rate_scale"]=n->GetHeartRateScale();
  v["nervous_heart_elastance_scale"]=n->GetHeartElastanceScale();
  v["nervous_compliance_scale"]=n->GetComplianceScale();
  v["nervous_muscle_resistance_scale"]=n->GetResistanceScaleMuscle();
  const auto* food=gi->GetStomachContents();
  const double missing=std::numeric_limits<double>::quiet_NaN();
  v["stomach_carbohydrate_g"]=food?food->GetCarbohydrate(MassUnit::g):missing;
  v["stomach_protein_g"]=food?food->GetProtein(MassUnit::g):missing;
  v["stomach_fat_g"]=food?food->GetFat(MassUnit::g):missing;
  v["stomach_sodium_g"]=food?food->GetSodium(MassUnit::g):missing;
  v["stomach_calcium_mg"]=food?food->GetCalcium(MassUnit::mg):missing;
  v["stomach_water_ml"]=food?food->GetWater(VolumeUnit::mL):missing;
  v["chyme_water_absorption_ml_per_min"]=gi->GetChymeAbsorptionRate(VolumePerTimeUnit::mL_Per_min);
  v["insulin_synthesis_pmol_per_min"]=bg.GetEndocrineSystem()->GetInsulinSynthesisRate(AmountPerTimeUnit::pmol_Per_min);
  v["glucagon_synthesis_pmol_per_min"]=bg.GetEndocrineSystem()->GetGlucagonSynthesisRate(AmountPerTimeUnit::pmol_Per_min);
  v["liver_glycogen_g"]=t->GetLiverGlycogen(MassUnit::g);
  v["muscle_glycogen_g"]=t->GetMuscleGlycogen(MassUnit::g);
  v["stored_protein_g"]=t->GetStoredProtein(MassUnit::g);
  v["stored_fat_g"]=t->GetStoredFat(MassUnit::g);
  v["oxygen_consumption_ml_per_min"]=t->GetOxygenConsumptionRate(VolumePerTimeUnit::mL_Per_min);
  v["co2_production_ml_per_min"]=t->GetCarbonDioxideProductionRate(VolumePerTimeUnit::mL_Per_min);
  v["respiratory_exchange_ratio"]=t->GetRespiratoryExchangeRatio();
  v["metabolic_rate_w"]=e->GetTotalMetabolicRate(PowerUnit::W);
  v["energy_deficit_w"]=e->GetEnergyDeficit(PowerUnit::W);
  v["core_temperature_c"]=e->GetCoreTemperature(TemperatureUnit::C);
  v["skin_temperature_c"]=e->GetSkinTemperature(TemperatureUnit::C);
  v["fatigue_fraction"]=e->GetFatigueLevel();
  v["urine_volume_ml"]=k->GetUrineVolume(VolumeUnit::mL);
  v["urine_production_ml_per_min"]=k->GetUrineProductionRate(VolumePerTimeUnit::mL_Per_min);
  v["glomerular_filtration_ml_per_min"]=k->GetGlomerularFiltrationRate(VolumePerTimeUnit::mL_Per_min);
  for(const std::string compartment : {"SmallIntestineChyme","SmallIntestineVasculature","LiverVasculature","Aorta","MuscleVasculature","Lymph","Bladder"}) {
    const auto* c=bg.GetCompartments().GetLiquidCompartment(compartment);
    v[compartment+".volume_ml"]=c?c->GetVolume(VolumeUnit::mL):missing;
    for(const std::string substance : {"Glucose","AminoAcids","Triacylglycerol","Sodium","Calcium","Lactate","Urea","Insulin","Glucagon"}) {
      const auto* s=bg.GetSubstanceManager().GetSubstance(substance);
      const auto* q=c&&s?c->GetSubstanceQuantity(*s):nullptr;
      v[compartment+"."+substance+".mass_g"]=q?q->GetMass(MassUnit::g):missing;
      v[compartment+"."+substance+".concentration_mg_per_dl"]=q?q->GetConcentration(MassPerVolumeUnit::mg_Per_dL):missing;
    }
  }
  return v;
}
