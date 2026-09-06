// Read-only native respiratory source/chest work observer. No actuator or energy writer.
#pragma once
#include <cassert>
#include <biogears/cdm/circuit/fluid/SEFluidCircuit.h>
#include <biogears/cdm/circuit/fluid/SEFluidCircuitPath.h>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>

namespace ihm_respiration {
struct ChestWork {
 double pressure_before_pa=0,pressure_after_pa=0;
 double compliance_before_m3_pa=0,compliance_after_m3_pa=0;
 double flow_m3_s=0,stroke_m3=0,constitutive_residual_m3=0;
 double endpoint_work_j=0,storage_change_j=0,storage_residual_j=0;
};
struct WorkReceipt {
 std::uint64_t start_tick=0,end_tick=0;
 double dt_s=0,generated_pressure_pa=0,external_pressure_pa=0,applied_pressure_pa=0;
 double source_flow_m3_s=0,source_stroke_m3=0,source_pressure_residual_pa=0,source_kcl_m3=0;
 double source_work_j=0,generated_work_j=0,external_work_j=0;
 double cumulative_source_work_j=0,cumulative_generated_work_j=0,cumulative_external_work_j=0;
 double lung_before_m3=0,lung_after_m3=0,lung_delta_m3=0,source_minus_lung_stroke_m3=0;
 double external_lung_proxy_work_j=0,external_proxy_difference_j=0;
 std::array<ChestWork,2> chest{};
};
class WorkObserver {
 enum class Phase{idle,prepared,solved,failed};Phase phase_=Phase::idle;
 bool accepted_=false;WorkReceipt candidate_{},last_{};
 biogears::SEFluidCircuit* circuit_=nullptr;
 std::array<biogears::SEFluidCircuitPath*,3> paths_{};
 static double finite(double value){if(!std::isfinite(value))throw std::runtime_error("Nonfinite respiratory work observation");return value;}
 static void require(bool condition,const char* message){if(!condition)throw std::runtime_error(message);}
 void live()const{require(phase_!=Phase::failed,"Respiratory work observer failed terminally");}
public:
 void fail(){phase_=Phase::failed;}
 // Call after native PreProcess and after the sole pressure modifier, before Process.
 // generated_pa MUST be captured before that modifier; current lung volume is a separate observation.
 void begin(biogears::SEFluidCircuit& circuit,std::uint64_t tick,double dt,double generated_pa,double external_pa,double lung_m3){
  using namespace biogears;
  try{
   live();require(phase_==Phase::idle,"Respiratory work candidate already active");
   require(tick<std::numeric_limits<std::uint64_t>::max(),"Respiratory tick overflow");
   require(!accepted_||tick==last_.end_tick,"Noncontiguous respiratory work interval");
   require(finite(dt)>0&&finite(lung_m3)>0,"Invalid respiratory interval or lung volume");
   circuit_=&circuit;
   constexpr std::array<const char*,3> names{"EnvironmentToRespiratoryMuscle","LeftPleuralCavityToRespiratoryMuscle","RightPleuralCavityToRespiratoryMuscle"};
   for(size_t i=0;i<3;++i){paths_[i]=circuit.GetPath(names[i]);require(paths_[i]!=nullptr,"Required respiratory source/chest path missing");}
   auto& driver=*paths_[0];require(circuit.IsReferenceNode(driver.GetSourceNode())&&driver.GetTargetNode().GetName()=="RespiratoryMuscle","Unexpected respiratory source orientation");
   unsigned incident=0;for(auto* p:circuit.GetPaths())if(&p->GetSourceNode()==&driver.GetTargetNode()||&p->GetTargetNode()==&driver.GetTargetNode())++incident;
   require(incident==3,"Unsupported additional respiratory muscle boundary");
   require(driver.HasNextPressureSource()&&driver.NumberOfNextElements()==1,"Respiratory driver is not a sole pressure source");
   candidate_=WorkReceipt{};candidate_.start_tick=tick;candidate_.dt_s=dt;candidate_.lung_before_m3=lung_m3;
   candidate_.generated_pressure_pa=finite(generated_pa);candidate_.external_pressure_pa=finite(external_pa);
   candidate_.applied_pressure_pa=finite(driver.GetNextPressureSource(PressureUnit::Pa));
   require(std::abs(candidate_.applied_pressure_pa-generated_pa-external_pa)<=1e-9,"Unbound respiratory pressure modifier");
   for(size_t i=0;i<2;++i){auto& p=*paths_[i+1];auto& x=candidate_.chest[i];
    require(p.GetTargetNode().GetName()=="RespiratoryMuscle"&&p.GetSourceNode().GetName()==std::string(i==0?"Left":"Right")+"PleuralCavity","Unexpected chest orientation");
    require(p.HasCompliance()&&p.HasNextCompliance()&&p.NumberOfNextElements()==1&&!p.HasNextPolarizedState(),"Unsupported chest constitutive element");
    x.pressure_before_pa=finite(p.GetSourceNode().GetPressure(PressureUnit::Pa)-p.GetTargetNode().GetPressure(PressureUnit::Pa));
    x.compliance_before_m3_pa=finite(p.GetCompliance(FlowComplianceUnit::m3_Per_Pa));
    x.compliance_after_m3_pa=finite(p.GetNextCompliance(FlowComplianceUnit::m3_Per_Pa));
    require(x.compliance_before_m3_pa>0&&x.compliance_after_m3_pa>0,"Nonpositive chest compliance");
   }
   phase_=Phase::prepared;
  }catch(...){fail();throw;}
 }
 // Read solved NEXT values before native PostProcess copies/resets them.
 // Supplied lung volume must likewise be the solved native volume, not the old CURRENT nodes.
 void solved(biogears::SEFluidCircuit& circuit,double lung_after_m3){
  using namespace biogears;
  try{
   live();require(phase_==Phase::prepared&&circuit_==&circuit,"Respiratory solved phase or circuit mismatch");
   auto& r=candidate_;auto& driver=*paths_[0];
   require(driver.GetNextPressureSource(PressureUnit::Pa)==r.applied_pressure_pa,"Respiratory source changed during solve");
   r.source_flow_m3_s=finite(driver.GetNextFlow(VolumePerTimeUnit::mL_Per_s))*1e-6;
   r.source_stroke_m3=-r.source_flow_m3_s*r.dt_s;
   r.source_pressure_residual_pa=finite(driver.GetTargetNode().GetNextPressure(PressureUnit::Pa)-driver.GetSourceNode().GetNextPressure(PressureUnit::Pa))-r.applied_pressure_pa;
   r.source_work_j=r.applied_pressure_pa*r.source_flow_m3_s*r.dt_s;
   r.generated_work_j=r.generated_pressure_pa*r.source_flow_m3_s*r.dt_s;
   r.external_work_j=r.external_pressure_pa*r.source_flow_m3_s*r.dt_s;
   r.lung_after_m3=finite(lung_after_m3);require(lung_after_m3>0,"Invalid solved lung volume");
   r.lung_delta_m3=lung_after_m3-r.lung_before_m3;r.source_minus_lung_stroke_m3=r.source_stroke_m3-r.lung_delta_m3;
   r.external_lung_proxy_work_j=-r.external_pressure_pa*r.lung_delta_m3;r.external_proxy_difference_j=r.external_work_j-r.external_lung_proxy_work_j;
   r.source_kcl_m3=r.source_flow_m3_s*r.dt_s;
   for(size_t i=0;i<2;++i){auto& p=*paths_[i+1];auto& x=r.chest[i];
    require(p.GetNextCompliance(FlowComplianceUnit::m3_Per_Pa)==x.compliance_after_m3_pa,"Chest compliance changed during solve");
    x.pressure_after_pa=finite(p.GetSourceNode().GetNextPressure(PressureUnit::Pa)-p.GetTargetNode().GetNextPressure(PressureUnit::Pa));
    x.flow_m3_s=finite(p.GetNextFlow(VolumePerTimeUnit::mL_Per_s))*1e-6;x.stroke_m3=x.flow_m3_s*r.dt_s;
    x.constitutive_residual_m3=x.stroke_m3-(x.compliance_after_m3_pa*x.pressure_after_pa-x.compliance_before_m3_pa*x.pressure_before_pa);
    x.endpoint_work_j=x.pressure_after_pa*x.stroke_m3;
    x.storage_change_j=.5*(x.compliance_after_m3_pa*x.pressure_after_pa*x.pressure_after_pa-x.compliance_before_m3_pa*x.pressure_before_pa*x.pressure_before_pa);
    x.storage_residual_j=x.endpoint_work_j-x.storage_change_j;r.source_kcl_m3+=x.stroke_m3;
    finite(x.endpoint_work_j);finite(x.storage_change_j);finite(x.storage_residual_j);
   }
   finite(r.source_work_j);finite(r.generated_work_j);finite(r.external_work_j);phase_=Phase::solved;
  }catch(...){fail();throw;}
 }
 // Only after successful full native advancement AND signed muscle finish.
 void commit(std::uint64_t end_tick){try{live();require(phase_==Phase::solved&&end_tick==candidate_.start_tick+1,"Unaccepted respiratory work interval");candidate_.end_tick=end_tick;candidate_.cumulative_source_work_j=finite(last_.cumulative_source_work_j+candidate_.source_work_j);candidate_.cumulative_generated_work_j=finite(last_.cumulative_generated_work_j+candidate_.generated_work_j);candidate_.cumulative_external_work_j=finite(last_.cumulative_external_work_j+candidate_.external_work_j);last_=candidate_;accepted_=true;phase_=Phase::idle;}catch(...){fail();throw;}}
 const WorkReceipt& receipt()const{live();require(phase_==Phase::idle&&accepted_,"No committed respiratory work receipt");return last_;}
};
}
