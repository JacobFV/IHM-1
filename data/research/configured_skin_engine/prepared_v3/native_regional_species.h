// Native owning-leaf sweat withdrawal. Audit records are not physiological stores.
#pragma once
#include <cassert>
#include <biogears/cdm/compartment/fluid/SELiquidCompartment.h>
#include <biogears/cdm/compartment/substances/SELiquidSubstanceQuantity.h>
#include <biogears/cdm/properties/SEProperties.h>
#include <biogears/cdm/substance/SESubstance.h>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>
namespace ihm_regional {
struct SweatLedger {
  std::array<unsigned,3> calls{};
  std::array<double,3> requested_mg{},owner_before_mg{},owner_after_mg{},waste_before_mg{},waste_after_mg{},residual_mg{},charge_to_waste_mmol{};
};
extern "C" SweatLedger* ihm_regional_sweat_ledger();
#ifdef IHM_REGIONAL_SPECIES_IMPLEMENTATION
extern "C" SweatLedger* ihm_regional_sweat_ledger(){static thread_local SweatLedger record;return &record;}
#endif
inline void withdraw_skin_ion(biogears::SELiquidCompartment& owner,biogears::SELiquidSubstanceQuantity& aggregate,
                             double requested_mg,biogears::SEScalarMass& waste,size_t index,int valence) {
  using namespace biogears;
  if(index>=3||owner.GetSubstanceQuantity(aggregate.GetSubstance())!=&aggregate||!std::isfinite(requested_mg)||requested_mg<0)
    throw std::invalid_argument("Invalid native sweat withdrawal receipt");
  const double before=aggregate.GetMass(MassUnit::mg),waste_before=waste.GetValue(MassUnit::mg);
  if(!std::isfinite(before)||before<0||requested_mg>before||!std::isfinite(waste_before)||waste.IsReadOnly())throw std::runtime_error("Native sweat withdrawal exceeds available owner or waste is unavailable");
  if(!owner.HasChildren()) {
    aggregate.GetMass().IncrementValue(-requested_mg,MassUnit::mg);
  } else if(requested_mg>0.) {
    std::vector<std::pair<SELiquidSubstanceQuantity*,double>> debits;double used=0.;
    const auto& leaves=owner.GetLeaves();
    if(leaves.empty())throw std::runtime_error("Native aggregate has no owning leaves");
    for(size_t i=0;i<leaves.size();++i) {
      auto* quantity=leaves[i]->GetSubstanceQuantity(aggregate.GetSubstance());
      if(!quantity||!quantity->HasMass())throw std::runtime_error("Missing native sweat leaf species");
      const double mass=quantity->GetMass(MassUnit::mg);
      const double debit=requested_mg==before?mass:i+1==leaves.size()?requested_mg-used:requested_mg*(mass/before);
      if(!std::isfinite(mass)||mass<0||debit<0||debit>mass||quantity->GetMass().IsReadOnly())throw std::runtime_error("Invalid native regional sweat allocation");
      debits.emplace_back(quantity,debit);used+=debit;
    }
    for(const auto& [quantity,debit]:debits)quantity->GetMass().IncrementValue(-debit,MassUnit::mg);
  }
  // Credit the native cumulative waste quantity with exactly the native request;
  // no Python sink, clipping, per-region extra waste, or read-only bypass.
  waste.IncrementValue(requested_mg,MassUnit::mg);
  auto& record=*ihm_regional_sweat_ledger();record.calls[index]++;
  record.requested_mg[index]=requested_mg;record.owner_before_mg[index]=before;
  record.owner_after_mg[index]=aggregate.GetMass(MassUnit::mg);record.waste_before_mg[index]=waste_before;
  record.waste_after_mg[index]=waste.GetValue(MassUnit::mg);
  record.residual_mg[index]=(record.owner_after_mg[index]-before)+(record.waste_after_mg[index]-waste_before);
  record.charge_to_waste_mmol[index]=valence*requested_mg/aggregate.GetSubstance().GetMolarMass(MassPerAmountUnit::mg_Per_mmol);
}
}
