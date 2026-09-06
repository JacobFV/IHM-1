// Source-native Umberger/Uchida muscle energetics. Adds probes, never body inertia.
#pragma once
#include <OpenSim/OpenSim.h>
#include <OpenSim/Simulation/Model/Umberger2010MuscleMetabolicsProbe.h>
#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

struct IHMMetabolicMusclePrior {
  double slow_twitch_fraction=0.5; // Held OpenSim generic parameter default.
  double provided_mass_kg=0.; // Zero selects native Fmax/tension*density*optimalLength.
};
struct IHMMetabolicMusclePower {
  std::string muscle_type;
  double metabolic_w,active_fiber_work_w,heat_w,analysis_mass_kg,slow_twitch_fraction;
};
struct IHMMetabolicPower {
  double total_muscle_metabolic_w=0.,active_fiber_work_w=0.,muscle_heat_w=0.,analysis_mass_kg=0.;
  std::map<std::string,IHMMetabolicMusclePower> muscles;
};

class IHMNativeMuscleMetabolism {
  OpenSim::Model* model_;
  OpenSim::Umberger2010MuscleMetabolicsProbe* probe_; // Owned by model.
  std::vector<const OpenSim::Muscle*> muscles_;
  static void valid(double x,const char* label) {
    if(!std::isfinite(x))throw std::runtime_error(std::string("Nonfinite native muscle metabolism: ")+label);
  }
public:
  // Install AFTER all source/added muscles, BEFORE final connections/initSystem.
  explicit IHMNativeMuscleMetabolism(OpenSim::Model& model,
      const std::map<std::string,IHMMetabolicMusclePrior>& priors={}) : model_(&model),probe_(nullptr) {
    std::set<std::string> names;
    for(const auto& muscle:model.getComponentList<OpenSim::Muscle>()) {
      const std::string type=muscle.getConcreteClassName();
      if(type!="Millard2012EquilibriumMuscle" && type!="Thelen2003Muscle")
        throw std::invalid_argument("Unaudited metabolic muscle implementation: "+type);
      if(!names.insert(muscle.getName()).second)throw std::invalid_argument("Duplicate metabolic muscle name");
      valid(muscle.getMaxIsometricForce(),"Fmax");valid(muscle.getOptimalFiberLength(),"optimal fiber length");
      if(muscle.getMaxIsometricForce()<=0||muscle.getOptimalFiberLength()<=0)throw std::invalid_argument("Nonpositive muscle morphology");
      muscles_.push_back(&muscle);
    }
    if(muscles_.empty())throw std::invalid_argument("No native muscles for metabolism");
    for(const auto& [name,p]:priors) {
      if(!names.count(name)||!std::isfinite(p.slow_twitch_fraction)||p.slow_twitch_fraction<0||p.slow_twitch_fraction>1||!std::isfinite(p.provided_mass_kg)||p.provided_mass_kg<0)
        throw std::invalid_argument("Unknown muscle or invalid metabolic prior");
    }
    auto owned=std::make_unique<OpenSim::Umberger2010MuscleMetabolicsProbe>(true,true,false,true);
    owned->setName("ihm_native_muscle_metabolics");owned->setOperation("value");owned->setGain(1.);
    owned->set_report_total_metabolics_only(false);
    owned->set_enforce_minimum_heat_rate_per_muscle(true);
    owned->set_aerobic_factor(1.5);owned->set_muscle_effort_scaling_factor(1.);
    owned->set_use_Bhargava_recruitment_model(true);
    owned->set_include_negative_mechanical_work(true);owned->set_forbid_negative_total_power(true);
    probe_=owned.get();model.addProbe(owned.release());
    for(const auto* muscle:muscles_) {
      const auto found=priors.find(muscle->getName());const auto p=found==priors.end()?IHMMetabolicMusclePrior{}:found->second;
      if(p.provided_mass_kg>0)probe_->addMuscle(muscle->getName(),p.slow_twitch_fraction,p.provided_mass_kg);
      else probe_->addMuscle(muscle->getName(),p.slow_twitch_fraction);
    }
  }
  IHMMetabolicPower sample(const SimTK::State& state) const {
    model_->realizeDynamics(state);
    const auto power=probe_->getProbeOutputs(state);
    if(power.size()!=static_cast<int>(muscles_.size()+2))throw std::runtime_error("Incomplete metabolic muscle coverage");
    valid(power[0],"total");valid(power[1],"basal");
    if(power[1]!=0.)throw std::runtime_error("Unexpected duplicate whole-body basal metabolic power");
    IHMMetabolicPower result;result.total_muscle_metabolic_w=power[0];double sum=0.;
    for(size_t i=0;i<muscles_.size();++i) {
      const auto& m=*muscles_[i];const double total=power[static_cast<int>(i)+2];
      const double force=m.getActiveFiberForce(state),velocity=m.getFiberVelocity(state),mass=probe_->getMuscleMass(m.getName());
      valid(total,"individual power");valid(force,"active force");valid(velocity,"fiber velocity");valid(mass,"analysis muscle mass");
      if(total<0||mass<=0)throw std::runtime_error("Out-of-domain native muscle metabolism");
      // This reproduces the source probe's active-force floor and signed work.
      const double work=-std::max(0.,force)*velocity,heat=total-work;
      result.muscles.emplace(m.getName(),IHMMetabolicMusclePower{m.getConcreteClassName(),total,work,heat,mass,probe_->getRatioSlowTwitchFibers(m.getName())});
      sum+=total;result.active_fiber_work_w+=work;result.muscle_heat_w+=heat;result.analysis_mass_kg+=mass;
    }
    if(std::abs(sum-power[0])>1e-10*std::max(1.,std::abs(power[0])))throw std::runtime_error("Metabolic total differs from native per-muscle sum");
    return result;
  }
};
