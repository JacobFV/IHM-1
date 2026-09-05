// Thin adapter for the upstream BioGears API; no physiological equations changed.
#include <cassert>
#include <biogears/cdm/engine/PhysiologyEngineTrack.h>
#include <biogears/cdm/properties/SEProperties.h>
#include <biogears/engine/BioGearsPhysiologyEngine.h>
#include <biogears/version.h>
#include <fstream>
#include <cmath>
#include <vector>
#include <biogears/cdm/patient/actions/SEExercise.h>
#include <biogears/cdm/patient/actions/SEHemorrhage.h>
#include <biogears/cdm/patient/actions/SESubstanceCompoundInfusion.h>
#include <biogears/cdm/substance/SESubstanceManager.h>
#include <iostream>
#include <memory>
#include <cstdlib>
using namespace biogears;
int main(int argc, char** argv) {
  double seconds = argc > 1 ? std::stod(argv[1]) : 60.0;
  if (!(seconds > 0 && seconds <= 3600)) return 4;
  std::string patient = argc > 2 ? argv[2] : "StandardMale";
  std::string state = argc > 3 ? argv[3] : "-";
  double hz = argc > 4 ? std::stod(argv[4]) : 50;
  if (hz != 1 && hz != 2 && hz != 5 && hz != 10 && hz != 25 && hz != 50) return 4;
  struct Event { double time; std::string kind; double value; };
  std::vector<Event> events;
  if (argc > 5) {
    std::ifstream timeline(argv[5]); if (!timeline) return 4;
    Event e; double previous = -1;
    while (timeline >> e.time >> e.kind >> e.value) {
      if (!std::isfinite(e.time) || !std::isfinite(e.value) || e.time < previous || e.time < 0 || e.time >= seconds || e.value < 0 ||
          (e.kind != "exercise" && e.kind != "hemorrhage" && e.kind != "saline") || e.value > (e.kind == "exercise" ? .5 : 100)) return 4;
      events.push_back(e); previous = e.time;
    }
    if (!timeline.eof()) return 4;
  }
  std::cout << "ENGINE_VERSION=" << full_version_string() << "\n";
  auto bg = CreateBioGearsEngine("native_engine.log");
  if (state == "-" ? !bg->InitializeEngine("patients/" + patient + ".xml") : !bg->LoadState(state)) return 2;
  std::cout << "STABILIZED_TIME_S=" << bg->GetSimulationTime(TimeUnit::s) << "\n";
  auto& dm = bg->GetEngineTrack()->GetDataRequestManager();
  dm.SetSamplesPerSecond(hz);
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
  const double dt = bg->GetTimeStep(TimeUnit::s);
  if (std::abs(dt - .02) > 1e-12) return 4;
  const long total_steps = std::lround(seconds / dt);
  if (std::abs(total_steps * dt - seconds) > 1e-8) return 4;
  const long sample_stride = std::lround(1.0 / hz / dt);
  // Schedule both actions and tracking with integer native steps. The upstream
  // duration overload floors floating quotients, and its tracker clock can drift.
  bg->SetAutoTrackFlag(false);
  long elapsed_steps = 0;
  auto advance_to = [&](long target) {
    while (elapsed_steps < target) {
      if (!bg->AdvanceModelTime(false)) return false;
      ++elapsed_steps;
      if (elapsed_steps % sample_stride == 0)
        bg->GetEngineTrack()->TrackData(bg->GetSimulationTime(TimeUnit::s), true);
    }
    return true;
  };
  for (auto& e : events) {
    const long event_step = std::lround(e.time / dt);
    if (std::abs(event_step * dt - e.time) > 1e-8) return 4;
    if (!advance_to(event_step)) return 3;
    bool ok = false;
    if (e.kind == "exercise") {
      SEExercise::SEGeneric generic; generic.Intensity.SetValue(e.value);
      SEExercise action { generic }; ok = bg->ProcessAction(action);
    } else if (e.kind == "hemorrhage") {
      SEHemorrhage action; action.SetCompartment("RightLeg");
      action.GetInitialRate().SetValue(e.value, VolumePerTimeUnit::mL_Per_min);
      action.SetMCIS(); ok = bg->ProcessAction(action);
    } else {
      auto saline = bg->GetSubstanceManager().GetCompound("Saline"); if (!saline) return 5;
      SESubstanceCompoundInfusion action(*saline);
      action.GetBagVolume().SetValue(500, VolumeUnit::mL);
      action.GetRate().SetValue(e.value, VolumePerTimeUnit::mL_Per_min);
      ok = bg->ProcessAction(action);
    }
    if (!ok) return 5;
    std::cout << "ACTION_TIME_S=" << elapsed_steps * dt << " KIND=" << e.kind << " VALUE=" << e.value << "\n";
  }
  if (!advance_to(total_steps)) return 3;
  bg->SaveStateToFile("native_final.xml");
  std::cout << "FINAL_TIME_S=" << bg->GetSimulationTime(TimeUnit::s) << "\n";
  return 0;
}
