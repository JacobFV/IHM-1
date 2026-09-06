// Isolated native circuit solver fixture. No patient initialization or gas chemistry.
#include "native_respiratory_work.h"
#include <biogears/cdm/circuit/SECircuitManager.h>
#include <biogears/cdm/circuit/fluid/SEFluidCircuitCalculator.h>
#include <iostream>
#include <iomanip>
using namespace biogears;
void near(double a,double b,double tolerance=1e-10){if(!std::isfinite(a)||!std::isfinite(b)||std::abs(a-b)>tolerance)throw std::runtime_error("Respiratory work identity failed");}
template<class F>void rejects(F f){bool rejected=false;try{f();}catch(const std::exception&){rejected=true;}if(!rejected)throw std::runtime_error("Expected rejection");}
struct Fixture {
 Logger logger;SECircuitManager manager;SEFluidCircuit& c;SEFluidCircuitCalculator solver;
 Fixture():logger("respiratory-work.log"),manager(&logger),c(manager.CreateFluidCircuit("RespiratoryWorkFixture")),solver(&logger) {
  auto node=[&](std::string name){auto& n=c.CreateNode(name);n.GetPressure().SetValue(101325,PressureUnit::Pa);n.GetNextPressure().SetValue(101325,PressureUnit::Pa);return &n;};
  auto* environment=node("Environment");c.AddReferenceNode(*environment);auto* muscle=node("RespiratoryMuscle");
  c.CreatePath(*environment,*muscle,"EnvironmentToRespiratoryMuscle").GetPressureSourceBaseline().SetValue(0,PressureUnit::Pa);
  for(const std::string side:{"Left","Right"}) {
   auto* alveoli=node(side+"Alveoli");auto* connection=node(side+"PleuralConnection");auto* pleural=node(side+"PleuralCavity");
   alveoli->GetVolumeBaseline().SetValue(1000,VolumeUnit::mL);pleural->GetVolumeBaseline().SetValue(8.5,VolumeUnit::mL);
   c.CreatePath(*environment,*alveoli,side+"Airway").GetResistanceBaseline().SetValue(2,FlowResistanceUnit::cmH2O_s_Per_L);
   c.CreatePath(*alveoli,*connection,side+"LungCompliance").GetComplianceBaseline().SetValue(.1,FlowComplianceUnit::L_Per_cmH2O);
   c.CreatePath(*connection,*pleural,side+"PleuralConnectionTo"+side+"PleuralCavity");
   c.CreatePath(*pleural,*muscle,side+"PleuralCavityToRespiratoryMuscle").GetComplianceBaseline().SetValue(.1,FlowComplianceUnit::L_Per_cmH2O);
  }
  c.SetNextAndCurrentFromBaselines();c.StateChange();
 }
 double volume(bool next=false){double sum=0;for(const std::string name:{"LeftAlveoli","RightAlveoli"}){auto* n=c.GetNode(name);sum+=next?n->GetNextVolume(VolumeUnit::mL):n->GetVolume(VolumeUnit::mL);}return sum*1e-6;}
 void drive(double pa){c.GetPath("EnvironmentToRespiratoryMuscle")->GetNextPressureSource().SetValue(pa,PressureUnit::Pa);}
};
void emit(const ihm_respiration::WorkReceipt& r){
 std::cout<<std::setprecision(17)<<"CASE {\"start_tick\":"<<r.start_tick<<",\"end_tick\":"<<r.end_tick<<",\"dt_s\":"<<r.dt_s
 <<",\"generated_pa\":"<<r.generated_pressure_pa<<",\"external_pa\":"<<r.external_pressure_pa<<",\"applied_pa\":"<<r.applied_pressure_pa
 <<",\"source_flow_m3_s\":"<<r.source_flow_m3_s<<",\"source_stroke_m3\":"<<r.source_stroke_m3<<",\"lung_delta_m3\":"<<r.lung_delta_m3
 <<",\"source_work_j\":"<<r.source_work_j<<",\"generated_work_j\":"<<r.generated_work_j<<",\"external_work_j\":"<<r.external_work_j
 <<",\"cumulative_source_work_j\":"<<r.cumulative_source_work_j<<",\"chest\":[";
 for(size_t i=0;i<2;++i){const auto& x=r.chest[i];if(i)std::cout<<',';std::cout<<"{\"pressure_before_pa\":"<<x.pressure_before_pa<<",\"pressure_after_pa\":"<<x.pressure_after_pa<<",\"compliance_before_m3_pa\":"<<x.compliance_before_m3_pa<<",\"compliance_after_m3_pa\":"<<x.compliance_after_m3_pa<<",\"flow_m3_s\":"<<x.flow_m3_s<<",\"endpoint_work_j\":"<<x.endpoint_work_j<<",\"storage_change_j\":"<<x.storage_change_j<<",\"storage_residual_j\":"<<x.storage_residual_j<<'}';}
 std::cout<<"]}\n";
}
int main(){try{
 Fixture observed,control;ihm_respiration::WorkObserver observer;
 double max_kcl=0,max_constitutive=0,max_partition=0,max_parity=0,cumulative=0;unsigned cases=0;bool positive=false,negative=false;
 for(int i=0;i<12;++i){
  const double generated=i<6?-100.:(i<9?-50.:0.),external=(i==3?20.:0.);observed.drive(generated+external);control.drive(generated+external);
  observer.begin(observed.c,i,.02,generated,external,observed.volume());
  observed.solver.Process(observed.c,.02);control.solver.Process(control.c,.02);
  observer.solved(observed.c,observed.volume(true));
  rejects([&]{(void)observer.receipt();}); // Never publish a solved candidate before accepted commit.
  observed.solver.PostProcess(observed.c);control.solver.PostProcess(control.c);observer.commit(i+1);
  const auto& r=observer.receipt();emit(r);cumulative+=r.source_work_j;near(r.cumulative_source_work_j,cumulative,0);near(r.source_work_j,r.generated_work_j+r.external_work_j);
  max_kcl=std::max(max_kcl,std::abs(r.source_kcl_m3));max_partition=std::max(max_partition,std::abs(r.source_work_j-r.generated_work_j-r.external_work_j));
  max_parity=std::max(max_parity,std::abs(observed.volume()-control.volume()));near(observed.volume(),control.volume(),0);
  for(const auto& x:r.chest){max_constitutive=std::max(max_constitutive,std::abs(x.constitutive_residual_m3));near(x.endpoint_work_j-x.storage_change_j,x.storage_residual_j);near(x.storage_residual_j,.5*x.compliance_after_m3_pa*std::pow(x.pressure_after_pa-x.pressure_before_pa,2));}
  near(r.source_kcl_m3,0,1e-12);near(r.source_pressure_residual_pa,0,1e-8);near(r.source_stroke_m3,r.lung_delta_m3,1e-12);
  positive|=r.source_work_j>1e-10;negative|=r.source_work_j< -1e-10;++cases;
 }
 // Observations must retain a changed constitutive parameter, not call its energy residual heat.
 observed.drive(-75);observed.c.GetPath("LeftPleuralCavityToRespiratoryMuscle")->GetNextCompliance().SetValue(.12,FlowComplianceUnit::L_Per_cmH2O);
 observer.begin(observed.c,12,.02,-75,0,observed.volume());observed.solver.Process(observed.c,.02);observer.solved(observed.c,observed.volume(true));observed.solver.PostProcess(observed.c);observer.commit(13);
 emit(observer.receipt());
 if(observer.receipt().chest[0].compliance_after_m3_pa==observer.receipt().chest[0].compliance_before_m3_pa)throw std::runtime_error("Lost variable compliance");
 near(observer.receipt().chest[0].constitutive_residual_m3,0,1e-12);
 // Later failure must retain no accepted candidate and cannot be retried.
 ihm_respiration::WorkObserver failed;observed.drive(0);failed.begin(observed.c,13,.02,0,0,observed.volume());observed.solver.Process(observed.c,.02);failed.solved(observed.c,observed.volume(true));failed.fail();rejects([&]{failed.commit(14);});rejects([&]{(void)failed.receipt();});
 ihm_respiration::WorkObserver wrong;rejects([&]{wrong.begin(observed.c,13,.02,10,0,observed.volume());});
 if(!positive)throw std::runtime_error("No delivered work exercised");
 std::cout<<std::setprecision(17)<<"RESULT {\"passed\":true,\"patient_initialized\":false,\"cases\":"<<cases+1<<",\"max_kcl_m3\":"<<max_kcl<<",\"max_constitutive_m3\":"<<max_constitutive<<",\"max_partition_j\":"<<max_partition<<",\"max_observer_parity_m3\":"<<max_parity<<",\"positive_source_work\":"<<positive<<",\"negative_source_work\":"<<negative<<"}\n";
 return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
