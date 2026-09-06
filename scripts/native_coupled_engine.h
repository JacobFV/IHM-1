// Optional load hook, preserving the pinned engine's step ordering and clocks.
#pragma once
#include <biogears/engine/Controller/BioGearsEngine.h>
#include <biogears/engine/Controller/BioGearsConfiguration.h>
#include <biogears/cdm/patient/SEPatient.h>
#include <biogears/cdm/circuit/fluid/SEFluidCircuitPath.h>
#include "native_body_ports.h"
#include "native_tissue_compression.h"
#include <memory>

class CoupledBioGearsEngine : public biogears::BioGearsEngine {
public:
  using BioGearsEngine::BioGearsEngine;
  bool coupling_enabled=false;
  double external_pressure_pa=0.,generated_driver_pa=0.,applied_driver_pa=0.,external_work_j=0.;
  std::unique_ptr<NativeTissueCompression> compression;

  bool AdvanceModelTime(bool appendDataTrack=false) override {
    using namespace biogears;
    if(!coupling_enabled)return BioGearsEngine::AdvanceModelTime(appendDataTrack);
    if(!IsReady() || m_Patient->IsEventActive(SEPatientEventType::IrreversibleState))return false;
    const double before_m3=GetRespiratorySystem()->GetTotalLungVolume(VolumeUnit::mL)*1e-6;
    PreProcess();
    auto* driver=GetCircuits().GetActiveRespiratoryCircuit().GetPath("EnvironmentToRespiratoryMuscle");
    generated_driver_pa=driver->GetNextPressureSource(PressureUnit::Pa);
    applied_driver_pa=generated_driver_pa+external_pressure_pa;
    // A zero boundary must not round-trip the native scalar through another
    // unit. Preserve its exact original representation for source parity.
    if(external_pressure_pa!=0.)driver->GetNextPressureSource().SetValue(applied_driver_pa,PressureUnit::Pa);
    if(compression)compression->after_preprocess();
    Process();PostProcess();
    if(compression)compression->after_postprocess();
    const double after_m3=GetRespiratorySystem()->GetTotalLungVolume(VolumeUnit::mL)*1e-6;
    external_work_j-=external_pressure_pa*(after_m3-before_m3);
    m_Patient->UpdateEvents(m_Configuration->GetTimeStep());
    m_CurrentTime->Increment(m_Configuration->GetTimeStep());
    m_SimulationTime->Increment(m_Configuration->GetTimeStep());
    const auto interval=1.0/GetEngineTrack()->GetDataRequestManager().GetSamplesPerSecond();
    m_timeSinceLastDataTrack+=m_Configuration->GetTimeStep(TimeUnit::s);
    if(m_timeSinceLastDataTrack>=interval){
      m_timeSinceLastDataTrack-=interval;
      if((m_isAutoTracking && EngineState::Active==m_State) || (m_State<EngineState::Active && m_areTrackingStabilization))
        GetEngineTrack()->TrackData(GetSimulationTime(TimeUnit::s),appendDataTrack);
    }
    return true;
  }
};
