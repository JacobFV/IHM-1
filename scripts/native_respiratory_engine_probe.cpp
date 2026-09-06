#include "native_respiratory_engine_probe.h"
#include <iostream>
#include <iomanip>
#include <memory>
void values(const std::map<std::string,double>& v){bool first=true;std::cout<<'{';for(const auto& [k,x]:v){if(!first)std::cout<<',';first=false;std::cout<<'"'<<k<<"\":";if(std::isfinite(x))std::cout<<x;else std::cout<<"null";}std::cout<<'}';}
int main(int argc,char** argv){try{
 using namespace biogears;
 if(argc!=4)return 2;const std::string mode=argv[1];
 if(mode!="baseline"&&mode!="observer"&&mode!="pulse_control"&&mode!="pulse_observer")return 3;
 auto engine=std::make_unique<RespiratoryEngineProbe>("respiratory-engine.log");
 if(!engine->LoadState(argv[2]))throw std::runtime_error("Native state load failed");engine->SetAutoTrackFlag(false);
 engine->original=mode=="baseline";engine->observing=mode=="observer"||mode=="pulse_observer";
 if(std::abs(engine->GetTimeStep(TimeUnit::s)-.02)>1e-12)throw std::runtime_error("Unexpected native timestep");
 auto emit=[&]{auto v=body_ports(*engine);auto& c=engine->GetCircuits().GetActiveRespiratoryCircuit();auto* source=c.GetPath("EnvironmentToRespiratoryMuscle");
  v["probe.source_pressure_pa"]=source->GetPressureSource(PressureUnit::Pa);v["probe.source_flow_m3_s"]=source->HasFlow()?source->GetFlow(VolumePerTimeUnit::mL_Per_s)*1e-6:std::numeric_limits<double>::quiet_NaN();
  for(const std::string side:{"Left","Right"}){auto* p=c.GetPath(side+"PleuralCavityToRespiratoryMuscle");v["probe."+side+".pleural_pressure_pa"]=p->GetSourceNode().GetPressure(PressureUnit::Pa);v["probe."+side+".chest_compliance_m3_pa"]=p->GetCompliance(FlowComplianceUnit::m3_Per_Pa);}
  std::cout<<std::setprecision(17)<<"RESULT {\"mode\":\""<<mode<<"\",\"tick\":"<<engine->tick<<",\"time_s\":"<<engine->GetSimulationTime(TimeUnit::s)<<",\"values\":";values(v);std::cout<<",\"work\":";values(engine->work_values());std::cout<<"}\n";
 };
 emit();for(int step=1;step<=8;++step){engine->external_pa=(mode=="pulse_control"||mode=="pulse_observer")&&step>=3&&step<=4?10.:0.;engine->AdvanceModelTime(false);emit();}
 engine->SaveStateToFile(argv[3]);
 return 0;
}catch(const std::exception& e){std::cerr<<"RESPIRATORY_ENGINE_FAILURE "<<e.what()<<'\n';return 1;}}
