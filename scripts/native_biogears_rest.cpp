// Thin adapter for the upstream BioGears API; no physiological equations changed.
#include <cassert>
#include <biogears/cdm/engine/PhysiologyEngineTrack.h>
#include <biogears/cdm/properties/SEProperties.h>
#include <biogears/engine/BioGearsPhysiologyEngine.h>
#include <biogears/version.h>
#include <fstream>
#include <iostream>
#include <memory>
#include <cstdlib>
using namespace biogears;
int main(int argc, char** argv) {
  double seconds = argc > 1 ? std::stod(argv[1]) : 60.0;
  if (!(seconds > 0 && seconds <= 3600)) return 4;
  std::cout << "ENGINE_VERSION=" << full_version_string() << "\n";
  auto bg = CreateBioGearsEngine("native_engine.log");
  if (!bg->InitializeEngine("patients/StandardMale.xml")) return 2;
  std::cout << "STABILIZED_TIME_S=" << bg->GetSimulationTime(TimeUnit::s) << "\n";
  auto& dm = bg->GetEngineTrack()->GetDataRequestManager();
  dm.SetSamplesPerSecond(50);
  dm.CreatePhysiologyDataRequest().Set("ArterialPressure", PressureUnit::mmHg);
  dm.CreatePhysiologyDataRequest().Set("TotalLungVolume", VolumeUnit::mL);
  dm.CreatePhysiologyDataRequest().Set("HeartRate", FrequencyUnit::Per_min);
  dm.CreatePhysiologyDataRequest().Set("MeanArterialPressure", PressureUnit::mmHg);
  dm.CreatePhysiologyDataRequest().Set("SystolicArterialPressure", PressureUnit::mmHg);
  dm.CreatePhysiologyDataRequest().Set("DiastolicArterialPressure", PressureUnit::mmHg);
  dm.CreatePhysiologyDataRequest().Set("CardiacOutput", VolumePerTimeUnit::L_Per_min);
  dm.CreatePhysiologyDataRequest().Set("BloodVolume", VolumeUnit::mL);
  dm.CreatePhysiologyDataRequest().Set("RespirationRate", FrequencyUnit::Per_min);
  dm.CreatePhysiologyDataRequest().Set("TidalVolume", VolumeUnit::mL);
  dm.CreatePhysiologyDataRequest().Set("OxygenSaturation");
  dm.CreatePhysiologyDataRequest().Set("ArterialCarbonDioxidePressure", PressureUnit::mmHg);
  dm.CreatePhysiologyDataRequest().Set("ArterialOxygenPressure", PressureUnit::mmHg);
  dm.CreatePhysiologyDataRequest().Set("ArterialBloodPH");
  dm.CreatePhysiologyDataRequest().Set("GlomerularFiltrationRate", VolumePerTimeUnit::mL_Per_min);
  dm.CreatePhysiologyDataRequest().Set("RenalBloodFlow", VolumePerTimeUnit::mL_Per_min);
  dm.CreatePhysiologyDataRequest().Set("UrineProductionRate", VolumePerTimeUnit::mL_Per_min);
  dm.CreatePhysiologyDataRequest().Set("InsulinSynthesisRate", AmountPerTimeUnit::pmol_Per_min);
  dm.CreatePhysiologyDataRequest().Set("CoreTemperature", TemperatureUnit::C);
  dm.CreatePhysiologyDataRequest().Set("SkinTemperature", TemperatureUnit::C);
  dm.CreatePhysiologyDataRequest().Set("TotalMetabolicRate", PowerUnit::W);
  dm.CreatePhysiologyDataRequest().Set("SweatRate", MassPerTimeUnit::g_Per_s);
  dm.SetResultsFilename("native_multisystem.csv");
  bg->SaveStateToFile("native_stabilized.xml");
  if (!bg->AdvanceModelTime(seconds, TimeUnit::s, true)) return 3;
  bg->SaveStateToFile("native_final.xml");
  std::cout << "FINAL_TIME_S=" << bg->GetSimulationTime(TimeUnit::s) << "\n";
  return 0;
}
