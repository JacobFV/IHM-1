#pragma once
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ihm_thermal {
struct EvaporationBudget {
  double total_w, sweat_evaporated_w, diffusion_w;
  double sweat_unevaporated_latent_w, wetted_fraction;
};
// The source producer reports whole-body sweat kg/s. Its latent power is
// distributed with a uniform flux prior over regional skin area; it is not a
// power per unit area until divided by total skin area. Sweat remaining liquid
// has already left the source fluid circuit, and creates no extra heat sink.
inline EvaporationBudget evaporation_budget(double total_sweat_latent_w,
    double total_area_m2, double region_fraction, double capacity_w_m2,
    double diffusion_wetted_fraction) {
  if (!std::isfinite(total_sweat_latent_w) || total_sweat_latent_w<0 ||
      !std::isfinite(total_area_m2) || total_area_m2<=0 ||
      !std::isfinite(region_fraction) || region_fraction<=0 || region_fraction>1 ||
      !std::isfinite(capacity_w_m2) || !std::isfinite(diffusion_wetted_fraction) ||
      diffusion_wetted_fraction<0 || diffusion_wetted_fraction>1)
    throw std::invalid_argument("Invalid regional evaporation budget");
  // This boundary models outward evaporation only, not condensation.
  const double potential=std::max(0.,capacity_w_m2);
  const double sweat_flux=total_sweat_latent_w/total_area_m2;
  const double evaporated_flux=std::min(sweat_flux,potential);
  const double wet=potential>0 ? evaporated_flux/potential : 0.;
  const double area=total_area_m2*region_fraction;
  const double sweat=evaporated_flux*area;
  const double diffusion=(1-wet)*diffusion_wetted_fraction*potential*area;
  return {sweat+diffusion,sweat,diffusion,
    std::max(0.,total_sweat_latent_w*region_fraction-sweat),wet};
}
}
