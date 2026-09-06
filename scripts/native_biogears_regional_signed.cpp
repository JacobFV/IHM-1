// Persistent command/observation bridge; the native engine is the sole state owner.
#include <cassert>
#include <biogears/cdm/properties/SEProperties.h>
#include <biogears/engine/BioGearsPhysiologyEngine.h>
#include <biogears/cdm/patient/actions/SEApnea.h>
#include <biogears/cdm/patient/actions/SEConsumeNutrients.h>
#include <biogears/cdm/patient/actions/SEExercise.h>
#include <biogears/version.h>
#include <biogears/cdm/scenario/SEActionManager.h>
#include <biogears/cdm/scenario/SEPatientActionCollection.h>
#include "native_body_ports.h"
#include "native_regional_coupled_engine.h"
#include "native_tissue_ports.h"
#include "native_signed_muscle_port.h"
#include "native_intake_receipts.h"
#include <memory>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <regex>
int main(int argc, char** argv) {
 try {
  if(argc!=4) return 4;
  auto bg=std::make_unique<RegionalCoupledBioGearsEngine>("native_engine.log");
  if(std::string(argv[2])=="-" || !bg->LoadState(argv[2])) return 2;
  auto* engine=dynamic_cast<BioGearsEngine*>(bg.get());
  if(!engine || std::abs(bg->GetTimeStep(TimeUnit::s)-.02)>1e-12) return 6;
  bg->SetAutoTrackFlag(false);
  try {bg->install_regions();}
  catch(const std::exception& error){std::cerr<<"REGIONAL_NATIVE_FAILURE: "<<error.what()<<std::endl;return 15;}
  const long horizon=std::stol(argv[3]); long ticks=0, seq=0;
  const double origin=bg->GetSimulationTime(TimeUnit::s);
  ihm_signed::Record last_port;bool has_port=false;std::string fixed_reference;
  auto pending_meal=[&](){return engine->GetActions().GetPatientActions().HasConsumeNutrients();};
  ihm_intake::Ledger intake(ihm_intake::epoch());
  if(pending_meal())throw std::runtime_error("Loaded pending intake has no accepted command identity");
  auto held_intake=[&]()->std::optional<ihm_intake::Payload>{
    if(!pending_meal())return std::nullopt;
    if(engine->GetState()!=EngineState::Active)throw std::runtime_error("Intake consumer is not in the audited Active phase");
    auto* held=engine->GetActions().GetPatientActions().GetConsumeNutrients();
    if(!held || held->HasNutritionFile())throw std::runtime_error("Unbound native intake action payload");
    return ihm_intake::capture(held->GetNutrition(),MassUnit::kg,MassUnit::g,MassUnit::mg,VolumeUnit::mL);
  };
  auto before_native=[&](){
    intake.before(held_intake(),static_cast<std::uint64_t>(seq),static_cast<std::uint64_t>(ticks),bg->GetSimulationTime(TimeUnit::s));
    std::cerr<<"INTAKE_ATTEMPT ";intake.write_attempt(std::cerr);std::cerr<<std::endl;
  };
  std::cout<<"ENGINE_VERSION="<<full_version_string()<<'\n';
  auto emit=[&](const std::string& status){
    std::cout<<std::setprecision(17)<<"IHM\t{\"sequence\":"<<seq<<",\"status\":\""<<status<<"\",\"time_s\":"<<bg->GetSimulationTime(TimeUnit::s)<<",\"elapsed_s\":"<<ticks*.02<<",\"origin_s\":"<<origin<<",\"values\":{";
    bool first=true;
    auto observations=body_ports(*engine);
    observations["patient_weight_kg"]=engine->GetPatient().GetWeight(MassUnit::kg);
    observations.merge(native_tissue_ports(*engine));
    observations["coupling.external_pressure_pa"]=bg->external_pressure_pa;
    observations["coupling.generated_driver_pa"]=bg->generated_driver_pa;
    observations["coupling.applied_driver_pa"]=bg->applied_driver_pa;
    observations["coupling.external_work_j"]=bg->external_work_j;
    if(has_port) {
      observations["coupling.muscle_delta_m_w"]=last_port.delta_m_W;
      observations["coupling.muscle_delta_h_w"]=last_port.delta_h_W;
      observations["coupling.muscle_delta_w_w"]=last_port.delta_w_W;
      observations["coupling.muscle_allowed_decrement_w"]=last_port.allowed_decrement_W;
      observations["coupling.muscle_requested_kcal"]=last_port.muscle_requested_kcal;
      observations["coupling.muscle_unmet_kcal"]=last_port.unmet_muscle_kcal;
      observations["coupling.native_heat_w"]=last_port.native_heat_W;
      observations["coupling.effective_heat_w"]=last_port.effective_heat_W;
      observations["coupling.muscle_heat_count"]=last_port.heat_count;
      observations["coupling.muscle_tissue_count"]=last_port.tissue_count;
      const std::array<std::string,4> readers={"cardiovascular","endocrine","nervous_metabolic_fraction","nervous_tmr"};
      for(std::size_t i=0;i<readers.size();++i) {
        const auto key="coupling.effective_reader_"+readers[i];
        observations[key+"_count"]=last_port.effective_reader_count[i];
        // A skipped consumer has no observed demand; zero watts would imply a measurement.
        if(last_port.effective_reader_count[i])observations[key+"_w"]=last_port.effective_reader_W[i];
      }
    }
    observations.merge(bg->regional_ports());
    for(const auto& [key,value]:observations) {
      if(!first) std::cout<<','; first=false;
      std::cout<<'\"'<<key<<"\":";
      if(std::isfinite(value)) std::cout<<value; else std::cout<<"null";
    }
    std::cout<<"},\"pending_meal\":"<<(pending_meal()?"true":"false")<<",\"intake\":";
    intake.write_json(std::cout);std::cout<<"}"<<std::endl;
  };
  emit("ready"); std::string line;
  while(std::getline(std::cin,line)) {
    std::istringstream input(line); std::string op,extra; long requested;
    if(!(input>>requested>>op) || requested!=seq+1) return 7;
    ++seq; bool ok=true;
    if(op=="snapshot") { if(input>>extra) ok=false; }
    else if(op=="signed_step") {
      ihm_signed::Record record;std::string reference;
      if(!(input>>reference>>record.delta_m_W>>record.delta_h_W>>record.delta_w_W) || (input>>extra) || ticks>=horizon)ok=false;
      else {
        ok=std::regex_match(reference,std::regex("[a-f0-9]{64}")) && (fixed_reference.empty() || fixed_reference==reference);
        for(double value:{record.delta_m_W,record.delta_h_W,record.delta_w_W})ok=ok&&std::isfinite(value)&&std::abs(value)<=5000.;
        record.sequence=ticks+1;record.start_s=bg->GetSimulationTime(TimeUnit::s);record.end_s=record.start_s+.02;record.reference_id=reference;
        if(ok)try {ihm_signed::validate(record);}catch(const std::exception&){ok=false;}
        if(ok) {
          try {
            ihm_signed::Scope scope(&record);
            before_native();
            if(!bg->AdvanceModelTime())throw std::runtime_error("Native advancement failed before intake commit");
            ihm_signed::finish(record);
            intake.commit(pending_meal(),static_cast<std::uint64_t>(ticks+1),bg->GetSimulationTime(TimeUnit::s));
            ++ticks;fixed_reference=reference;last_port=record;has_port=true;
          } catch(const std::exception& error) {
            intake.terminal_failure();
            std::cerr<<std::setprecision(17)<<"SIGNED_NATIVE_FAILURE: "<<error.what()
              <<" delta_m_w="<<record.delta_m_W<<" delta_h_w="<<record.delta_h_W<<" delta_w_w="<<record.delta_w_W
              <<" allowed_decrement_w="<<record.allowed_decrement_W<<" heat_count="<<record.heat_count
              <<" tissue_count="<<record.tissue_count<<std::endl;
            return 14; // State may have changed: terminate, never claim rollback.
          }
        }
      }
    }
    else if(op=="step") {
      long count=0;
      if(!fixed_reference.empty() || !(input>>count) || count<1 || count>horizon-ticks || (input>>extra)) ok=false;
      else for(long i=0;i<count;++i) {
        before_native();
        if(!bg->AdvanceModelTime())throw std::runtime_error("Native advancement failed before intake commit");
        intake.commit(pending_meal(),static_cast<std::uint64_t>(ticks+1),bg->GetSimulationTime(TimeUnit::s));
        ++ticks;
      }
    } else if(op=="regional_skin_pressure") {
      std::string region;double value;
      if(!(input>>region>>value) || !std::isfinite(value) || value<0 || value>5000 || (input>>extra))ok=false;
      else {
        const auto& names=ihm::NativeRegionalSkin::regions;
        auto found=std::find(names.begin(),names.end(),region);
        if(found==names.end())ok=false;
        else bg->regional->set_pressure_pa(static_cast<size_t>(found-names.begin()),value);
      }
    } else if(op=="skin_compression") {
      ok=false; // Whole-Skin compression is incompatible with regional topology.
    } else if(op=="respiratory_load") {
      double value;
      if(!(input>>value) || !std::isfinite(value) || value < -5000. || value>5000. || (input>>extra))ok=false;
      else bg->external_pressure_pa=value;
    } else if(op=="apnea" || op=="exercise") {
      double value;
      if(!(input>>value) || !std::isfinite(value) || value<0 || value>(op=="exercise"?.5:1.) || (input>>extra)) ok=false;
      else if(op=="apnea") { SEApnea a; a.GetSeverity().SetValue(value); ok=bg->ProcessAction(a); }
      else {
        // Mutating GetGenericExercise() leaves the default action mode at NONE.
        SEExercise::SEGeneric generic; generic.Intensity.SetValue(value);
        SEExercise a { generic }; ok=bg->ProcessAction(a);
      }
    } else if(op=="meal") {
      std::string name; double c,p,f,s,ca,w;
      if(!(input>>name>>c>>p>>f>>s>>ca>>w) || (input>>extra)) ok=false;
      else {
        ok=!pending_meal() && std::regex_match(name,std::regex("[A-Za-z0-9_]{1,64}"));
        for(double x:{c,p,f,s,ca,w}) ok=ok && std::isfinite(x) && x>=0 && x<=10000;
        ok=ok && c+p+f+s+ca+w>0;
        if(ok) {
          SEConsumeNutrients a; auto& n=a.GetNutrition(); n.SetName(name);
          n.GetCarbohydrate().SetValue(c,MassUnit::g); n.GetProtein().SetValue(p,MassUnit::g);
          n.GetFat().SetValue(f,MassUnit::g); n.GetSodium().SetValue(s,MassUnit::g);
          n.GetCalcium().SetValue(ca,MassUnit::mg); n.GetWater().SetValue(w,VolumeUnit::mL);
          ok=bg->ProcessAction(a);
          if(ok) {
            auto held=held_intake();
            if(!held)throw std::runtime_error("Accepted native intake action was not retained");
            intake.accepted(static_cast<std::uint64_t>(seq),*held);
          } else if(pending_meal())throw std::runtime_error("Rejected intake left a native pending action");
        }
      }
    } else if(op=="save") {
      ok=false; // Regional topology and adapter ledgers are not serializable.
    } else if(op=="quit") {
      if(input>>extra) {emit("rejected");continue;}
      if(pending_meal()) {emit("rejected");continue;}
      emit("closed"); return 0;
    } else ok=false;
    emit(ok?"ok":"rejected");
  }
  if(pending_meal()) {
    std::cerr<<"UNCONSUMED_ACTION: stdin ended before the pending meal was consumed\n";
    return 9;
  }
  return 0;
 }catch(const std::exception& error){std::cerr<<"INTAKE_NATIVE_FAILURE: "<<error.what()<<std::endl;return 16;}
}

