// Isolated whole-engine observation composition; never installed by a production adapter.
#pragma once
#include <cassert>
#include <biogears/engine/Controller/BioGearsEngine.h>
#include <biogears/engine/Controller/BioGearsConfiguration.h>
#include <biogears/cdm/patient/SEPatient.h>
#include "native_body_ports.h"
#include "native_signed_muscle_port.h"
#include "native_respiratory_work.h"
class RespiratoryEngineProbe : public biogears::BioGearsEngine {
public:
 using BioGearsEngine::BioGearsEngine;
 bool original=false,observing=false,terminal=false;
 std::uint64_t tick=0;
 double external_pa=0,generated_pa=0,applied_pa=0,lung_before_m3=0;
 double request_before_pa=0,request_after_pa=0,frequency_before_per_min=0,frequency_after_per_min=0,cycle_before=0,cycle_after=0;
 ihm_respiration::WorkObserver observer;
 void before_preprocess(){
  using namespace biogears;
  const auto* respiratory=GetRespiratorySystem();
  lung_before_m3=respiratory->GetTotalLungVolume(VolumeUnit::mL)*1e-6;
  request_before_pa=respiratory->GetRespirationDriverPressure(PressureUnit::Pa);
  frequency_before_per_min=respiratory->GetRespirationDriverFrequency(FrequencyUnit::Per_min);
  cycle_before=respiratory->GetRespirationCyclePercentComplete();
 }
 void after_preprocess(){
  using namespace biogears;
  auto& c=GetCircuits().GetActiveRespiratoryCircuit();auto* driver=c.GetPath("EnvironmentToRespiratoryMuscle");
  if(!driver)throw std::runtime_error("Missing native respiratory driver");
  generated_pa=driver->GetNextPressureSource(PressureUnit::Pa);applied_pa=generated_pa+external_pa;
  if(external_pa!=0)driver->GetNextPressureSource().SetValue(applied_pa,PressureUnit::Pa);
  request_after_pa=GetRespiratorySystem()->GetRespirationDriverPressure(PressureUnit::Pa);
  frequency_after_per_min=GetRespiratorySystem()->GetRespirationDriverFrequency(FrequencyUnit::Per_min);
  cycle_after=GetRespiratorySystem()->GetRespirationCyclePercentComplete();
  if(observing)observer.begin(c,tick,GetTimeStep(TimeUnit::s),generated_pa,external_pa,lung_before_m3);
 }
 void after_process(){
  if(observing)observer.solved(GetCircuits().GetActiveRespiratoryCircuit(),GetRespiratorySystem()->GetTotalLungVolume(biogears::VolumeUnit::mL)*1e-6);
 }
#include "native_respiratory_engine_step.inc"
 bool AdvanceModelTime(bool appendDataTrack=false)override{
  using namespace biogears;
  if(terminal)throw std::runtime_error("Respiratory probe is terminal");
  try{
   // All modes use the same zero signed delta; no respiratory work enters this record.
   ihm_signed::Record record;record.sequence=tick+1;record.start_s=GetSimulationTime(TimeUnit::s);record.end_s=record.start_s+GetTimeStep(TimeUnit::s);record.reference_id="respiratory-observer-zero-delta";
   ihm_signed::validate(record);ihm_signed::Scope scope(&record);
   if(!(original?BioGearsEngine::AdvanceModelTime(appendDataTrack):probe_advance(appendDataTrack)))throw std::runtime_error("Native respiratory probe advance failed");
   ihm_signed::finish(record);
   if(observing)observer.commit(tick+1);
   ++tick;return true;
  }catch(...){terminal=true;observer.fail();throw;}
 }
 std::map<std::string,double> work_values()const{
  if(!observing||tick==0)return {};
  const auto& r=observer.receipt();std::map<std::string,double> v{
   {"dt_s",r.dt_s},{"generated_pressure_pa",r.generated_pressure_pa},{"external_pressure_pa",r.external_pressure_pa},{"applied_pressure_pa",r.applied_pressure_pa},
   {"source_flow_m3_s",r.source_flow_m3_s},{"source_stroke_m3",r.source_stroke_m3},{"source_pressure_residual_pa",r.source_pressure_residual_pa},{"source_kcl_m3",r.source_kcl_m3},
   {"source_work_j",r.source_work_j},{"generated_work_j",r.generated_work_j},{"external_work_j",r.external_work_j},
   {"cumulative_source_work_j",r.cumulative_source_work_j},{"cumulative_generated_work_j",r.cumulative_generated_work_j},{"cumulative_external_work_j",r.cumulative_external_work_j},
   {"lung_before_m3",r.lung_before_m3},{"lung_after_m3",r.lung_after_m3},{"lung_delta_m3",r.lung_delta_m3},{"source_minus_lung_stroke_m3",r.source_minus_lung_stroke_m3},
   {"external_lung_proxy_work_j",r.external_lung_proxy_work_j},{"external_proxy_difference_j",r.external_proxy_difference_j},
   {"request_before_preprocess_pa",request_before_pa},{"request_after_preprocess_pa",request_after_pa},{"frequency_before_preprocess_per_min",frequency_before_per_min},{"frequency_after_preprocess_per_min",frequency_after_per_min},{"cycle_before_preprocess",cycle_before},{"cycle_after_preprocess",cycle_after}};
  for(size_t i=0;i<2;++i){const auto& x=r.chest[i];const std::string key=i==0?"left.":"right.";
   v[key+"pressure_before_pa"]=x.pressure_before_pa;v[key+"pressure_after_pa"]=x.pressure_after_pa;v[key+"compliance_before_m3_pa"]=x.compliance_before_m3_pa;v[key+"compliance_after_m3_pa"]=x.compliance_after_m3_pa;
   v[key+"flow_m3_s"]=x.flow_m3_s;v[key+"stroke_m3"]=x.stroke_m3;v[key+"constitutive_residual_m3"]=x.constitutive_residual_m3;
   v[key+"endpoint_work_j"]=x.endpoint_work_j;v[key+"storage_change_j"]=x.storage_change_j;v[key+"storage_residual_j"]=x.storage_residual_j;
  }return v;
 }
};
