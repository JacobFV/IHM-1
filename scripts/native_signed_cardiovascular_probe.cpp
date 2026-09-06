// Isolated whole-engine probe; no shared adapter or default library selection.
#include <cassert>
#include "native_coupled_engine.h"
#include "native_signed_muscle_port.h"
#include "native_signed_vascular_prior.h"
#include <biogears/engine/Systems/Energy.h>
#include <dlfcn.h>
#include <iostream>
#include <iomanip>
int main(int argc,char** argv){try{
 if(argc!=3)return 2;const std::string mode=argv[2];
 CoupledBioGearsEngine bg("engine.log");if(!bg.LoadState(argv[1]))throw std::runtime_error("load failed");
 bg.SetAutoTrackFlag(false);bg.coupling_enabled=true;
 using TraceGetter=ihm_signed_vascular::Trace*(*)();auto trace=reinterpret_cast<TraceGetter>(dlsym(RTLD_DEFAULT,"ihm_signed_vascular_trace"));
 for(int i=0;i<3;++i){
  ihm_signed::Record r;r.sequence=i+1;r.reference_id="explicit-isolated-vascular-probe";r.start_s=bg.GetSimulationTime(TimeUnit::s);r.end_s=r.start_s+.02;
  r.delta_m_W=i==0?(mode=="positive"?.01:mode=="negative"?-.01:0):0;r.delta_h_W=r.delta_m_W;
  if(mode=="disabled"){if(!bg.AdvanceModelTime())throw std::runtime_error("advance failed");}
  else{ihm_signed::validate(r);ihm_signed::Scope scope(&r);if(!bg.AdvanceModelTime())throw std::runtime_error("signed advance failed");}
  auto values=body_ports(bg);
  values["exercise_map_target_mmhg"]=bg.GetEnergy().GetExerciseMeanArterialPressureDelta(PressureUnit::mmHg);
  for(const auto& name:{"Aorta1ToMuscle1","Muscle1ToMuscle2","Muscle2ToVenaCava"}){
   auto* path=bg.GetCircuits().GetActiveCardiovascularCircuit().GetPath(name);if(!path)throw std::runtime_error("missing muscle path probe");
   values[std::string("has_resistance.")+name]=path->HasResistance()?1:0;
   values[std::string("has_cardiovascular_region.")+name]=path->HasCardiovascularRegion()?1:0;
   values[std::string("resistance.")+name]=path->HasResistance()?path->GetResistance(FlowResistanceUnit::mmHg_s_Per_mL):std::numeric_limits<double>::quiet_NaN();
   values[std::string("flow.")+name]=path->GetFlow(VolumePerTimeUnit::mL_Per_s);
  }
  std::cout<<std::setprecision(17)<<"IHM\t{\"step\":"<<i+1<<",\"delta_m_w\":"<<r.delta_m_W<<",\"reader_count\":"<<r.effective_reader_count[0]<<",\"reader_w\":"<<r.effective_reader_W[0]<<",\"heat_count\":"<<r.heat_count<<",\"tissue_count\":"<<r.tissue_count<<",\"unmet_kcal\":"<<r.unmet_muscle_kcal<<",\"values\":{";
  bool first=true;for(const auto& [key,value]:values){if(!first)std::cout<<',';first=false;std::cout<<'"'<<key<<"\":";if(std::isfinite(value))std::cout<<std::setprecision(17)<<value;else std::cout<<"null";}std::cout<<"},\"trace\":";
  if(trace){auto* t=trace();std::cout<<"{\"sequence\":"<<t->sequence<<",\"background_w\":"<<t->background_W<<",\"basal_w\":"<<t->basal_W<<",\"delta_w\":"<<t->delta_W<<",\"ratio\":"<<t->ratio<<",\"muscle_before\":"<<t->muscle_before<<",\"muscle_after\":"<<t->muscle_after<<",\"other_before\":"<<t->other_before<<",\"other_after\":"<<t->other_after<<'}';}else std::cout<<"null";
  std::cout<<'}'<<std::endl;
 }
 return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
