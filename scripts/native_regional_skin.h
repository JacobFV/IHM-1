// Experimental native circuit split. No Python reservoirs or transport solver.
// Two engineering regions plus residual; original law objects remain detached
// from the active solve. Parent extracellular compartment becomes aggregate only.
#pragma once
#include <cassert>
#include <biogears/cdm/circuit/fluid/SEFluidCircuit.h>
#include <biogears/cdm/compartment/SECompartmentManager.h>
#include <biogears/cdm/compartment/fluid/SELiquidCompartment.h>
#include <biogears/cdm/compartment/substances/SELiquidSubstanceQuantity.h>
#include <biogears/cdm/properties/SEProperties.h>
#include <array>
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace ihm {
using namespace biogears;
class NativeRegionalSkin {
public:
  static constexpr std::array<double,3> fractions{.2,.3,.5};
  static inline const std::array<std::string,3> regions{"region_a","region_b","residual"};
  static inline const std::array<std::string,5> node_names{"SkinE1","SkinE2","SkinE3","SkinL1","SkinL2"};
  static inline const std::array<std::string,9> path_names{
    "SkinVToSkinE1","SkinE1ToSkinE2","SkinE2ToSkinE3","SkinE3ToGround",
    "SkinE3ToSkinI","SkinE3ToSkinL1","SkinL1ToSkinL2","SkinToLymphValve","SkinSweating"};
  std::map<std::string,SEFluidCircuitNode*> original_nodes;
  std::map<std::string,SEFluidCircuitPath*> original_paths;
  std::array<std::map<std::string,SEFluidCircuitNode*>,3> nodes;
  std::array<std::map<std::string,SEFluidCircuitPath*>,3> paths;
  std::array<SELiquidCompartment*,3> children{};
  std::array<SEFluidCircuitPath*,3> drives{};
  std::array<double,3> requested_pa{};
  SELiquidCompartment& parent;

  template<class Scalar,class Unit> static void assign(Scalar& target,double value,const Unit& unit) {
    bool ro=target.IsReadOnly();target.SetReadOnly(false);target.SetValue(value,unit);target.SetReadOnly(ro);
  }
  static std::array<double,3> partition(double total) {
    if(!std::isfinite(total))throw std::invalid_argument("Nonfinite regional allocation");
    double a=total*fractions[0],b=total*fractions[1];return {a,b,total-a-b};
  }
  static std::string name(size_t i,const std::string& native) {return "IHM_"+regions.at(i)+"_"+native;}

  NativeRegionalSkin(SEFluidCircuit& circuit,SECompartmentManager& compartments,SELiquidCompartment& owner)
    :parent(owner) {
    if(parent.HasChildren()||!parent.HasNodeMapping())throw std::invalid_argument("Expected native leaf with node mapping");
    auto mapped=parent.GetNodeMapping().GetNodes();
    if(mapped.size()!=node_names.size())throw std::invalid_argument("Unexpected Skin node ownership");
    for(const auto& n:node_names) {
      auto* node=circuit.GetNode(n);
      if(!node||std::find(mapped.begin(),mapped.end(),node)==mapped.end())throw std::invalid_argument("Skin node identity mismatch");
      original_nodes[n]=node;
    }
    auto* ground=circuit.GetNode("Ground");
    if(!ground)throw std::invalid_argument("Missing Ground");
    for(const auto& n:path_names) {
      auto* path=circuit.GetPath(n);
      if(!path)throw std::invalid_argument("Missing native Skin path: "+n);
      if(path->HasInertance()||path->HasNextInertance()||path->HasSwitch()||path->HasNextSwitch()||path->HasPolarizedState()||path->HasNextPolarizedState())
        throw std::invalid_argument("Unsupported native Skin path state");
      original_paths[n]=path;
    }
    for(auto* path:circuit.GetPaths()) {
      bool touches=original_nodes.count(path->GetSourceNode().GetName())||original_nodes.count(path->GetTargetNode().GetName());
      if(touches&&!original_paths.count(path->GetName()))throw std::invalid_argument("Unmapped regional boundary path: "+path->GetName());
    }
    for(size_t i=0;i<3;++i) {
      for(const auto& n:node_names)if(circuit.GetNode(name(i,n)))throw std::invalid_argument("Duplicate regional node");
      if(compartments.GetLiquidCompartment(name(i,"SkinTissueExtracellular")))throw std::invalid_argument("Duplicate regional compartment");
    }
    std::vector<std::pair<SESubstance*,double>> masses;
    for(auto* q:parent.GetSubstanceQuantities())if(q->HasMass()) {
      double mass=q->GetMass(MassUnit::ug);
      if(!std::isfinite(mass)||mass<0)throw std::invalid_argument("Invalid source species mass");
      masses.emplace_back(&q->GetSubstance(),mass);
    }
    // Snapshot each extensive scalar BEFORE changing the owner's node mapping.
    for(size_t i=0;i<3;++i) {
      auto& child=compartments.CreateLiquidCompartment(name(i,"SkinTissueExtracellular"));children[i]=&child;
      if(parent.HasWaterVolumeFraction())child.GetWaterVolumeFraction().SetValue(parent.GetWaterVolumeFraction().GetValue());
      if(parent.HasPH())child.GetPH().SetValue(parent.GetPH().GetValue());
      for(const auto& [n,source]:original_nodes) {
        auto& target=circuit.CreateNode(name(i,n));nodes[i][n]=&target;
#define IHM_NODE_INTENSIVE(Slot) if(source->Has##Slot())assign(target.Get##Slot(),source->Get##Slot(PressureUnit::mmHg),PressureUnit::mmHg)
        IHM_NODE_INTENSIVE(Pressure);IHM_NODE_INTENSIVE(NextPressure);
#undef IHM_NODE_INTENSIVE
#define IHM_NODE_EXTENSIVE(Slot) if(source->Has##Slot())assign(target.Get##Slot(),partition(source->Get##Slot(VolumeUnit::mL))[i],VolumeUnit::mL)
        IHM_NODE_EXTENSIVE(VolumeBaseline);IHM_NODE_EXTENSIVE(Volume);IHM_NODE_EXTENSIVE(NextVolume);
#undef IHM_NODE_EXTENSIVE
        child.MapNode(target);
      }
      auto& external=circuit.CreateNode(name(i,"ExternalPressure"));
      assign(external.GetPressure(),ground->GetPressure(PressureUnit::mmHg),PressureUnit::mmHg);
      assign(external.GetNextPressure(),ground->GetPressure(PressureUnit::mmHg),PressureUnit::mmHg);
      drives[i]=&circuit.CreatePath(*ground,external,name(i,"ExternalPressureDrive"));
      for(auto* scalar:{&drives[i]->GetPressureSourceBaseline(),&drives[i]->GetPressureSource(),&drives[i]->GetNextPressureSource()})assign(*scalar,0.,PressureUnit::Pa);
      for(const auto& [n,source]:original_paths) {
        auto* a=&source->GetSourceNode();auto* b=&source->GetTargetNode();
        if(nodes[i].count(a->GetName()))a=nodes[i].at(a->GetName());
        if(nodes[i].count(b->GetName()))b=nodes[i].at(b->GetName());
        if(n=="SkinE3ToGround")b=&external;
        auto& target=circuit.CreatePath(*a,*b,name(i,n));paths[i][n]=&target;
        copy_laws(*source,target,i);
#define IHM_FLOW(Slot) if(source->Has##Slot())assign(target.Get##Slot(),partition(source->Get##Slot(VolumePerTimeUnit::mL_Per_s))[i],VolumePerTimeUnit::mL_Per_s)
        IHM_FLOW(Flow);IHM_FLOW(NextFlow);
#undef IHM_FLOW
        if(source->HasValve())target.SetValve(source->GetValve());
        if(source->HasNextValve())target.SetNextValve(source->GetNextValve());
      }
      for(const auto& [sub,mass]:masses) {
        auto* q=child.GetSubstanceQuantity(*sub);
        if(!q)throw std::runtime_error("Active species missing in regional child");
        assign(q->GetMass(),partition(mass)[i],MassUnit::ug);
      }
    }
    // Native AddChild rejects node-mapped parents. Move mapping first; native
    // parent quantity getters subsequently derive read-only sums from children.
    for(auto* node:mapped)parent.GetNodeMapping().RemoveNode(*node);
    for(auto* child:children)parent.AddChild(*child);
    parent.StateChange();
    for(const auto& [n,path]:original_paths)circuit.RemovePath(*path);
    for(const auto& [n,node]:original_nodes)circuit.RemoveNode(*node);
    circuit.StateChange();
    for(auto* child:children)child->StateChange();
    parent.Balance(BalanceLiquidBy::Mass);
  }

  static void copy_laws(SEFluidCircuitPath& source,SEFluidCircuitPath& target,size_t i) {
#define IHM_LAW(Slot,Unit,Value) if(source.Has##Slot())assign(target.Get##Slot(),(Value),Unit)
#define IHM_RES(Slot) IHM_LAW(Slot,FlowResistanceUnit::mmHg_s_Per_mL,source.Get##Slot(FlowResistanceUnit::mmHg_s_Per_mL)/fractions[i])
    IHM_RES(ResistanceBaseline);IHM_RES(Resistance);IHM_RES(NextResistance);
#undef IHM_RES
#define IHM_COMP(Slot) IHM_LAW(Slot,FlowComplianceUnit::mL_Per_mmHg,partition(source.Get##Slot(FlowComplianceUnit::mL_Per_mmHg))[i])
    IHM_COMP(ComplianceBaseline);IHM_COMP(Compliance);IHM_COMP(NextCompliance);
#undef IHM_COMP
#define IHM_PRESS(Slot) IHM_LAW(Slot,PressureUnit::mmHg,source.Get##Slot(PressureUnit::mmHg))
    IHM_PRESS(PressureSourceBaseline);IHM_PRESS(PressureSource);IHM_PRESS(NextPressureSource);
#undef IHM_PRESS
#define IHM_SOURCE(Slot) IHM_LAW(Slot,VolumePerTimeUnit::mL_Per_s,partition(source.Get##Slot(VolumePerTimeUnit::mL_Per_s))[i])
    IHM_SOURCE(FlowSourceBaseline);IHM_SOURCE(FlowSource);IHM_SOURCE(NextFlowSource);
#undef IHM_SOURCE
#undef IHM_LAW
  }
  void set_pressure_pa(size_t region,double pressure) {
    if(region>=3||!std::isfinite(pressure)||pressure<0||pressure>5000)throw std::invalid_argument("Invalid regional pressure");
    requested_pa[region]=pressure;
  }
  void after_preprocess() {
    for(size_t i=0;i<3;++i) {
      for(const auto& [n,source]:original_paths)copy_laws(*source,*paths[i].at(n),i);
      assign(drives[i]->GetNextPressureSource(),requested_pa[i],PressureUnit::Pa);
    }
  }
  void after_postprocess() {
    // Detached originals are law/observation caches, never active fluid stores.
    for(const auto& [n,source]:original_paths) {
      double q=0,nq=0;bool has=true,next=true;
      for(size_t i=0;i<3;++i) {auto* p=paths[i].at(n);has&=p->HasFlow();next&=p->HasNextFlow();if(p->HasFlow())q+=p->GetFlow(VolumePerTimeUnit::mL_Per_s);if(p->HasNextFlow())nq+=p->GetNextFlow(VolumePerTimeUnit::mL_Per_s);}
      if(has)assign(source->GetFlow(),q,VolumePerTimeUnit::mL_Per_s);
      if(next)assign(source->GetNextFlow(),nq,VolumePerTimeUnit::mL_Per_s);
      // Mirror SECircuitCalculator's commit for detached native law pointers.
#define IHM_COMMIT(Slot,Unit) if(source->HasNext##Slot())Override<Unit>(source->GetNext##Slot(),source->Get##Slot())
      IHM_COMMIT(Resistance,FlowResistanceUnit);IHM_COMMIT(Compliance,FlowComplianceUnit);
      IHM_COMMIT(PressureSource,PressureUnit);IHM_COMMIT(FlowSource,VolumePerTimeUnit);
#undef IHM_COMMIT
    }
  }
};
}
