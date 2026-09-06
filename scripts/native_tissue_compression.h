// Optional mechanical pressure boundary. Native fluid/solute states remain sole owners.
// Install after LoadState; call after_preprocess before cardiovascular Process,
// and after_postprocess after the native circuit commits its step.
#pragma once
#include "native_tissue_ports.h"
#include <cmath>
#include <stdexcept>
class NativeTissueCompression {
  biogears::SEFluidCircuitPath *original_,*replacement_,*drive_;
  biogears::SEFluidCircuitNode *reference_;
  double requested_pa_=0;
public:
  explicit NativeTissueCompression(BioGearsEngine& engine,const std::string& organ="Skin") {
    if(organ!="Skin")throw std::invalid_argument("Compression topology currently audited only for native lumped Skin");
    auto& circuit=engine.GetCircuits().GetActiveCardiovascularCircuit();
    original_=circuit.GetPath(organ+"E3ToGround");
    if(!original_||!original_->HasCompliance()||!original_->HasNextCompliance()||!original_->HasComplianceBaseline()||original_->HasPressureSource()||original_->HasFlowSource())
      throw std::runtime_error("Missing or non-compliance native tissue boundary");
    auto* ground=circuit.GetNode("Ground");auto* tissue=circuit.GetNode(organ+"E3");
    if(!ground||!tissue||circuit.GetNode("IHM"+organ+"ExternalPressure"))throw std::runtime_error("Invalid or duplicate compression topology");
    reference_=&circuit.CreateNode("IHM"+organ+"ExternalPressure");
    reference_->GetPressure().SetValue(ground->GetPressure(PressureUnit::Pa),PressureUnit::Pa);
    reference_->GetNextPressure().SetValue(ground->GetPressure(PressureUnit::Pa),PressureUnit::Pa);
    drive_=&circuit.CreatePath(*ground,*reference_,"IHMGroundTo"+organ+"ExternalPressure");
    drive_->GetPressureSourceBaseline().SetValue(0,PressureUnit::Pa);drive_->GetPressureSource().SetValue(0,PressureUnit::Pa);drive_->GetNextPressureSource().SetValue(0,PressureUnit::Pa);
    replacement_=&circuit.CreatePath(*tissue,*reference_,"IHM"+organ+"E3ToExternalPressure");
    biogears::Override<FlowComplianceUnit>(original_->GetComplianceBaseline(),replacement_->GetComplianceBaseline());
    biogears::Override<FlowComplianceUnit>(original_->GetCompliance(),replacement_->GetCompliance());
    biogears::Override<FlowComplianceUnit>(original_->GetNextCompliance(),replacement_->GetNextCompliance());
    if(original_->HasFlow())biogears::Override<VolumePerTimeUnit>(original_->GetFlow(),replacement_->GetFlow());
    if(original_->HasNextFlow())biogears::Override<VolumePerTimeUnit>(original_->GetNextFlow(),replacement_->GetNextFlow());
    // Remove from active circuit, not manager: existing native pointers stay valid.
    circuit.RemovePath(*original_);circuit.StateChange();
  }
  void set_pressure_pa(double value) {
    if(!std::isfinite(value)||value<0)throw std::invalid_argument("External compression requires finite nonnegative Pa");
    requested_pa_=value;
  }
  void after_preprocess() {
    // Preserve all current native compliance changes and COP pressure sources.
    biogears::Override<FlowComplianceUnit>(original_->GetComplianceBaseline(),replacement_->GetComplianceBaseline());
    biogears::Override<FlowComplianceUnit>(original_->GetCompliance(),replacement_->GetCompliance());
    biogears::Override<FlowComplianceUnit>(original_->GetNextCompliance(),replacement_->GetNextCompliance());
    drive_->GetNextPressureSource().SetValue(requested_pa_,PressureUnit::Pa);
  }
  void after_postprocess() {
    // The detached native law object no longer receives circuit PostProcess.
    // Use the same per-scalar Override protocol as SECircuitCalculator,
    // preserving read-only flags and original stored units bit-for-bit.
    biogears::Override<FlowComplianceUnit>(replacement_->GetCompliance(),original_->GetCompliance());
    biogears::Override<FlowComplianceUnit>(replacement_->GetNextCompliance(),original_->GetNextCompliance());
    biogears::Override<VolumePerTimeUnit>(replacement_->GetFlow(),original_->GetFlow());
    biogears::Override<VolumePerTimeUnit>(replacement_->GetNextFlow(),original_->GetNextFlow());
  }
  std::map<std::string,double> ports() const {
    return {{"tissue.compression.Skin.requested_pa",requested_pa_},
      {"tissue.compression.Skin.applied_pa",drive_->GetPressureSource(PressureUnit::Pa)},
      {"tissue.compression.Skin.reference_pa",reference_->GetPressure(PressureUnit::Pa)},
      {"tissue.path.SkinE3ToGround.compliance_ml_per_mmhg",original_->GetCompliance(FlowComplianceUnit::mL_Per_mmHg)},
      {"tissue.path.SkinE3ToGround.flow_ml_per_s",replacement_->GetFlow(VolumePerTimeUnit::mL_Per_s)}};
  }
};
