// Tiny native circuit/compartment fixture: no patient stabilization or advancement.
#include "native_regional_skin.h"
#include <biogears/engine/Controller/BioGears.h>
#include <biogears/cdm/circuit/fluid/SEFluidCircuitCalculator.h>
#include <iomanip>
#include <iostream>
#include <memory>

using namespace biogears;
struct Fixture {
  BioGears bg;
  SEFluidCircuit& circuit;
  SELiquidCompartment& owner;
  SEFluidCircuitCalculator calculator;
  std::unique_ptr<ihm::NativeRegionalSkin> split;
  Fixture(const std::string& label,bool regional):bg(label+".log"),
    circuit(bg.GetCircuits().CreateFluidCircuit("RegionalFixture")),
    owner(bg.GetCompartments().CreateLiquidCompartment("SkinTissueExtracellular")),calculator(bg.GetLogger()) {
    if(!bg.GetSubstances().LoadSubstanceDirectory())throw std::runtime_error("Substance directory unavailable");
    for(const std::string sub:{"Albumin","Glucose","Sodium"})
      static_cast<SECompartmentManager&>(bg.GetCompartments()).AddLiquidCompartmentSubstance(*bg.GetSubstances().GetSubstance(sub));
    auto node=[&](std::string name,double pressure){auto& n=circuit.CreateNode(name);n.GetPressure().SetValue(pressure,PressureUnit::mmHg);n.GetNextPressure().SetValue(pressure,PressureUnit::mmHg);return &n;};
    auto* ground=node("Ground",0);circuit.AddReferenceNode(*ground);
    auto* vascular=node("SkinV",20);auto* lymph=node("Lymph",0);auto* intracellular=node("SkinI",10);
    auto pressure=[&](SEFluidCircuitNode& a,SEFluidCircuitNode& b,std::string name,double p){auto& path=circuit.CreatePath(a,b,name);path.GetPressureSourceBaseline().SetValue(p,PressureUnit::mmHg);return &path;};
    pressure(*ground,*vascular,"FixtureBloodPressure",20);pressure(*ground,*lymph,"FixtureLymphPressure",0);pressure(*ground,*intracellular,"FixtureIntracellularPressure",10);
    auto* e1=node("SkinE1",18);auto* e2=node("SkinE2",8);auto* e3=node("SkinE3",10);auto* l1=node("SkinL1",10);auto* l2=node("SkinL2",0);
    e3->GetVolumeBaseline().SetValue(100,VolumeUnit::mL);
    pressure(*vascular,*e1,"SkinVToSkinE1",-2);
    circuit.CreatePath(*e1,*e2,"SkinE1ToSkinE2").GetResistanceBaseline().SetValue(100,FlowResistanceUnit::mmHg_s_Per_mL);
    pressure(*e2,*e3,"SkinE2ToSkinE3",2);
    circuit.CreatePath(*e3,*ground,"SkinE3ToGround").GetComplianceBaseline().SetValue(2,FlowComplianceUnit::mL_Per_mmHg);
    circuit.CreatePath(*e3,*intracellular,"SkinE3ToSkinI").GetFlowSourceBaseline().SetValue(0,VolumePerTimeUnit::mL_Per_s);
    pressure(*e3,*l1,"SkinE3ToSkinL1",0);
    circuit.CreatePath(*l1,*l2,"SkinL1ToSkinL2").GetResistanceBaseline().SetValue(100,FlowResistanceUnit::mmHg_s_Per_mL);
    circuit.CreatePath(*l2,*lymph,"SkinToLymphValve").SetNextValve(SEOpenClosed::Closed);
    circuit.CreatePath(*e3,*ground,"SkinSweating").GetFlowSourceBaseline().SetValue(0,VolumePerTimeUnit::mL_Per_s);
    circuit.SetNextAndCurrentFromBaselines();circuit.StateChange();
    for(auto* n:{e1,e2,e3,l1,l2})owner.MapNode(*n);
    owner.StateChange();
    for(const std::string sub:{"Albumin","Glucose","Sodium"}) {
      auto* q=owner.GetSubstanceQuantity(*bg.GetSubstances().GetSubstance(sub));
      q->GetMass().SetValue(sub=="Albumin"?12345.67:765.43,MassUnit::ug);q->Balance(BalanceLiquidBy::Mass);
    }
    if(regional)split=std::make_unique<ihm::NativeRegionalSkin>(circuit,bg.GetCompartments(),owner);
  }
  void step(double pressure=0) {
    if(split){split->set_pressure_pa(0,pressure);split->after_preprocess();}
    calculator.Process(circuit,.02);calculator.PostProcess(circuit);
    if(split)split->after_postprocess();
    owner.Balance(BalanceLiquidBy::Mass);
  }
  double q() {return bg.GetCircuits().GetFluidPath("SkinE3ToSkinL1")->GetFlow(VolumePerTimeUnit::mL_Per_s);}
};

int main() {
 try {
  Fixture baseline("baseline",false),regional("regional",true);
  auto& split=*regional.split;
  auto near=[](double a,double b,double tolerance){if(!std::isfinite(a)||!std::isfinite(b)||std::abs(a-b)>tolerance)throw std::runtime_error("Native parity/ownership invariant failed: "+std::to_string(a)+" vs "+std::to_string(b));};
  if(regional.owner.HasNodeMapping()||!regional.owner.HasChildren()||regional.owner.GetLeaves().size()!=3)throw std::runtime_error("Parent retained node ownership");
  const auto& leaves=regional.bg.GetCompartments().GetLiquidLeafCompartments();
  if(std::find(leaves.begin(),leaves.end(),&regional.owner)!=leaves.end())throw std::runtime_error("Parent remained in global owning-leaf inventory");
  for(auto* child:split.children)if(std::find(leaves.begin(),leaves.end(),child)==leaves.end())throw std::runtime_error("Regional child missing from global owning-leaf inventory");
  if(regional.circuit.GetPath("SkinE3ToSkinL1")!=nullptr||regional.bg.GetCircuits().GetFluidPath("SkinE3ToSkinL1")==nullptr)throw std::runtime_error("Detached source-law lookup invariant");
  for(const std::string sub:{"Albumin","Glucose","Sodium"}) {
    auto* s=regional.bg.GetSubstances().GetSubstance(sub);double children=0;
    for(auto* child:split.children)children+=child->GetSubstanceQuantity(*s)->GetMass(MassUnit::ug);
    near(children,regional.owner.GetSubstanceQuantity(*s)->GetMass(MassUnit::ug),1e-10);
    near(children,baseline.owner.GetSubstanceQuantity(*baseline.bg.GetSubstances().GetSubstance(sub))->GetMass(MassUnit::ug),1e-10);
    if(!regional.owner.GetSubstanceQuantity(*s)->GetMass().IsReadOnly())throw std::runtime_error("Parent mass still writable");
  }
  double volume=0,next=0,base=0,compliance=0;
  for(size_t i=0;i<3;++i){auto* n=split.nodes[i].at("SkinE3");volume+=n->GetVolume(VolumeUnit::mL);next+=n->GetNextVolume(VolumeUnit::mL);base+=n->GetVolumeBaseline(VolumeUnit::mL);compliance+=split.paths[i].at("SkinE3ToGround")->GetComplianceBaseline(FlowComplianceUnit::mL_Per_mmHg);}
  near(volume,100,1e-12);near(next,100,1e-12);near(base,100,1e-12);near(compliance,2,1e-12);
  double max_volume=0,max_flow=0;
  for(int i=0;i<10;++i) {
    // Emulate a native law changing its NEXT pressure/C/R values. Detached
    // originals must remain writable and their commits must reach all regions.
    double cop=2.+i*.01;
    baseline.bg.GetCircuits().GetFluidPath("SkinE2ToSkinE3")->GetNextPressureSource().SetValue(cop,PressureUnit::mmHg);
    split.original_paths.at("SkinE2ToSkinE3")->GetNextPressureSource().SetValue(cop,PressureUnit::mmHg);
    baseline.step();regional.step();
    max_volume=std::max(max_volume,std::abs(baseline.owner.GetVolume(VolumeUnit::mL)-regional.owner.GetVolume(VolumeUnit::mL)));
    max_flow=std::max(max_flow,std::abs(baseline.q()-regional.q()));
    near(baseline.owner.GetVolume(VolumeUnit::mL),regional.owner.GetVolume(VolumeUnit::mL),1e-9);
    near(baseline.q(),regional.q(),1e-10);
  }
  regional.step(133.322387415);
  double p0=split.nodes[0].at("SkinE3")->GetPressure(PressureUnit::mmHg),p1=split.nodes[1].at("SkinE3")->GetPressure(PressureUnit::mmHg);
  if(std::abs(p0-p1)<.1)throw std::runtime_error("Local pressure failed to separate regional states");
  double q0=split.paths[0].at("SkinE3ToSkinL1")->GetFlow(VolumePerTimeUnit::mL_Per_s)/.2;
  double q1=split.paths[1].at("SkinE3ToSkinL1")->GetFlow(VolumePerTimeUnit::mL_Per_s)/.3;
  if(std::abs(q0-q1)<1e-5)throw std::runtime_error("Local pressure failed to separate regional drainage");
  regional.step(0);
  near(split.drives[0]->GetPressureSource(PressureUnit::Pa),0,0);
  std::cout<<"RESULT {\"passed\":true,\"regions\":3,\"zero_load_steps\":10,\"max_volume_error_ml\":"<<std::setprecision(17)<<max_volume<<",\"max_flow_error_ml_per_s\":"<<max_flow<<",\"local_pressure_difference_mmhg\":"<<p0-p1<<",\"local_specific_flow_difference_ml_per_s\":"<<q0-q1<<",\"patient_initialized\":false,\"protein_transport_executed\":false}\n";
  return 0;
 } catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
