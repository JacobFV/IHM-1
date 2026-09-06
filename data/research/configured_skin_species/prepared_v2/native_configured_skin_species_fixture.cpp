// Bounded native species operators; no patient initialization or physiology advance.
#define main configured_circuit_reference_main
#include "native_configured_regional_skin_fixture.cpp"
#undef main
#include <biogears/engine/Systems/Diffusion.h>
#include <biogears/cdm/substance/SESubstanceTransport.h>
#include <biogears/cdm/compartment/fluid/SELiquidCompartmentLink.h>
#include <algorithm>
#include <limits>

// The held CDM library hides this template specialization. Instantiate the exact
// pinned native .inl implementation locally; do not replace its transport law.
template class biogears::SESubstanceTransporter<biogears::SELiquidTransportGraph,
  biogears::VolumePerTimeUnit,biogears::VolumeUnit,biogears::MassUnit,biogears::MassPerVolumeUnit>;

namespace biogears {
// Native Diffusion explicitly grants this test-access friend. No access macro or
// copied diffusion method; calls dispatch into the pinned graph_v2 shared library.
class BioGearsEngineTest {
public:
  static void lymph(BioGears& bg,SELiquidCompartment& source,SELiquidCompartment& target,
                    const SESubstance& sub,double dt) {
    auto& diffusion=bg.GetDiffusionCalculator();diffusion.m_dt_s=dt;
    diffusion.CalculatePassiveLymphDiffusion(source,target,sub);
  }
};
}

struct SpeciesSeed {const char* name;double concentration_ug_per_ml;};
#include "configured_species_seed.inc"

struct SpeciesFixture : Fixture {
  SELiquidCompartment& lymph;
  SELiquidCompartmentGraph& graph;
  SELiquidTransporter transporter;
  std::vector<SESubstance*> substances;
  std::vector<double> totals;
  SpeciesFixture(const std::string& name,bool regional):Fixture(name,regional),
    lymph(bg.GetCompartments().CreateLiquidCompartment("Lymph")),
    graph(bg.GetCompartments().CreateLiquidGraph("SelectedExtracellularLymph")),
    transporter(VolumePerTimeUnit::mL_Per_s,VolumeUnit::mL,MassUnit::ug,MassPerVolumeUnit::ug_Per_mL,bg.GetLogger()) {
    lymph.GetVolume().SetValue(10,VolumeUnit::mL);
    for(const auto& seed:species_seeds) {
      auto* sub=bg.GetSubstances().GetSubstance(seed.name);
      if(!sub)throw std::runtime_error("Missing native species: "+std::string(seed.name));
      static_cast<SECompartmentManager&>(bg.GetCompartments()).AddLiquidCompartmentSubstance(*sub);
      substances.push_back(sub);
      const double total=seed.concentration_ug_per_ml*100.;totals.push_back(total);
      auto set=[&](SELiquidCompartment& c,double mass) {
        auto* q=c.GetSubstanceQuantity(*sub);
        if(!q)throw std::runtime_error("Missing native liquid species quantity");
        q->GetMass().SetValue(mass,MassUnit::ug);q->Balance(BalanceLiquidBy::Mass);
      };
      if(split)for(size_t i=0;i<3;++i)set(*split->children[i],split->partition(total)[i]);
      else set(owner,total);
      set(lymph,0.);
    }
    owner.Balance(BalanceLiquidBy::Mass);
    graph.AddCompartment(lymph);
    auto edge=[&](SELiquidCompartment& source,SEFluidCircuitPath& path,const std::string& name) {
      graph.AddCompartment(source);
      auto& link=bg.GetCompartments().CreateLiquidLink(source,lymph,name);
      link.MapPath(path);graph.AddLink(link);
    };
    if(split)for(size_t i=0;i<3;++i)edge(*split->children[i],*split->paths[i].at("SkinE3ToSkinL1"),"SelectedDrain"+std::to_string(i));
    else edge(owner,*bg.GetCircuits().GetFluidPath("SkinE3ToSkinL1"),"SelectedDrain");
    graph.StateChange();
    if(split&&graph.GetCompartment(owner.GetName()))throw std::runtime_error("Nonowning parent is a transport vertex");
    if(graph.GetLinks().size()!=(split?3:1))throw std::runtime_error("Duplicate selected fluid edge");
  }
  double mass(SELiquidCompartment& c,SESubstance& s){return c.GetSubstanceQuantity(s)->GetMass(MassUnit::ug);}
  double verify() {
    double error=0;
    for(size_t j=0;j<substances.size();++j) {
      auto& sub=*substances[j];double owning=mass(owner,sub);
      if(split) {
        owning=0.;
        for(auto* child:split->children){const double value=mass(*child,sub);if(!std::isfinite(value)||value<0)throw std::runtime_error("Invalid owning leaf mass");owning+=value;}
        if(std::abs(owning-mass(owner,sub))>1e-7)throw std::runtime_error("Aggregate mass ownership mismatch");
        if(!owner.GetSubstanceQuantity(sub)->GetMass().IsReadOnly())throw std::runtime_error("Writable aggregate mass");
      }
      const double destination=mass(lymph,sub);
      if(!std::isfinite(owning)||owning<0||!std::isfinite(destination)||destination<0)throw std::runtime_error("Invalid donor/destination mass");
      const double residual=std::abs(owning+destination-totals[j]);error=std::max(error,residual);
      if(residual>1e-6)throw std::runtime_error("Species mass ledger failed: "+sub.GetName());
    }
    return error;
  }
  void transport(double dt=.02) {transporter.Transport(graph,dt);owner.Balance(BalanceLiquidBy::Mass);lymph.Balance(BalanceLiquidBy::Mass);}
  void albumin(double dt=.02) {
    BioGearsEngineTest::lymph(bg,owner,lymph,bg.GetSubstances().GetAlbumin(),dt);
    owner.Balance(BalanceLiquidBy::Mass);lymph.Balance(BalanceLiquidBy::Mass);
  }
};

int main() {
 try {
  double max_ledger=0,max_parity=0,loaded_difference=0,donor_after=0;
  SpeciesFixture baseline("species_baseline",false),regional("species_regional",true);
  for(int step=0;step<10;++step) {
    baseline.step();regional.step();baseline.transport();regional.transport();
    max_ledger=std::max({max_ledger,baseline.verify(),regional.verify()});
    for(size_t i=0;i<baseline.substances.size();++i) {
      double error=std::abs(baseline.mass(baseline.owner,*baseline.substances[i])-regional.mass(regional.owner,*regional.substances[i]));
      max_parity=std::max(max_parity,error);
      if(error>1e-6)throw std::runtime_error("Zero-load aggregate transported mass parity failed");
    }
  }
  auto& albumin=regional.bg.GetSubstances().GetAlbumin();
  double before0=regional.mass(*regional.split->children[0],albumin),before1=regional.mass(*regional.split->children[1],albumin);
  regional.step(133.322387415);regional.transport();max_ledger=std::max(max_ledger,regional.verify());
  loaded_difference=(before0-regional.mass(*regional.split->children[0],albumin))/regional.split->fractions[0]
    -(before1-regional.mass(*regional.split->children[1],albumin))/regional.split->fractions[1];
  if(std::abs(loaded_difference)<1e-8)throw std::runtime_error("No pressure-localized transported Albumin response");
  regional.step(0);regional.transport();max_ledger=std::max(max_ledger,regional.verify());
  if(regional.split->drives[0]->GetPressureSource(PressureUnit::Pa)!=0)throw std::runtime_error("Pressure release failed");
  // Separate phase: native physiological Albumin lymph operator on the aggregate.
  // It must not execute in addition to the generic transport phase on one state.
  SpeciesFixture passive_base("albumin_baseline",false),passive_regional("albumin_regional",true);
  for(int step=0;step<10;++step) {
    passive_base.step();passive_regional.step();passive_base.albumin();passive_regional.albumin();
    max_ledger=std::max({max_ledger,passive_base.verify(),passive_regional.verify()});
    const double error=std::abs(passive_base.mass(passive_base.owner,passive_base.bg.GetSubstances().GetAlbumin())-
                              passive_regional.mass(passive_regional.owner,passive_regional.bg.GetSubstances().GetAlbumin()));
    max_parity=std::max(max_parity,error);
    if(error>1e-6)throw std::runtime_error("Native passive Albumin zero-load parity failed");
  }
  passive_regional.step(133.322387415);passive_regional.albumin();max_ledger=std::max(max_ledger,passive_regional.verify());
  passive_regional.step(0);passive_regional.albumin();max_ledger=std::max(max_ledger,passive_regional.verify());
  // Actual native donor cap: a deliberately huge operator dt requests more Albumin
  // than available; no patient or circuit is advanced through this interval.
  passive_regional.albumin(1e12);max_ledger=std::max(max_ledger,passive_regional.verify());
  donor_after=passive_regional.mass(passive_regional.owner,passive_regional.bg.GetSubstances().GetAlbumin());
  if(donor_after>1e-7)throw std::runtime_error("Native donor cap did not exhaust available Albumin");
  passive_regional.albumin(1e12);max_ledger=std::max(max_ledger,passive_regional.verify());
  std::cout<<std::setprecision(17)<<"RESULT {\"passed\":true,\"species_count\":"<<species_seeds.size()
    <<",\"zero_load_steps_per_operator\":10,\"maximum_mass_ledger_residual_ug\":"<<max_ledger
    <<",\"maximum_zero_load_aggregate_mass_error_ug\":"<<max_parity
    <<",\"local_specific_albumin_transfer_difference_ug\":"<<loaded_difference
    <<",\"albumin_donor_after_cap_ug\":"<<donor_after
    <<",\"native_liquid_transporter_executed\":true,\"native_passive_albumin_executed\":true,\"pressure_release_passed\":true,\"patient_initialized\":false}\n";
  return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
