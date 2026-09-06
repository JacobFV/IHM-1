bool AdvanceModelTime(bool appendDataTrack=false) override
{
  // TODO: I am starting to think this should be a protected function
  if (!IsReady()) {
    return false;
  }
  if (m_Patient->IsEventActive(SEPatientEventType::IrreversibleState)) {
    return false;
  }

  PreProcess();
  if(regional)regional->after_preprocess();
  Process();
  PostProcess();
  if(regional)regional->after_postprocess();

  m_Patient->UpdateEvents(m_Configuration->GetTimeStep());
  m_CurrentTime->Increment(m_Configuration->GetTimeStep());
  m_SimulationTime->Increment(m_Configuration->GetTimeStep());

  const auto sample_interval_s = 1.0 / GetEngineTrack()->GetDataRequestManager().GetSamplesPerSecond();
  m_timeSinceLastDataTrack += m_Configuration->GetTimeStep(TimeUnit::s);
  if (m_timeSinceLastDataTrack >= sample_interval_s) {
    m_timeSinceLastDataTrack -= sample_interval_s;
    if ((m_isAutoTracking && EngineState::Active == m_State) || (m_State < EngineState::Active && m_areTrackingStabilization)) {
      auto time = GetSimulationTime(TimeUnit::s);
      GetEngineTrack()->TrackData(time, appendDataTrack);
    }
  }

  return true;
}
