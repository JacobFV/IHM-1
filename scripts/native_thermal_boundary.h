#pragma once
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ihm_thermal {
// rsi is the whole-body area-normalized clothing insulation, m^2 K/W.
// A uniform local insulation has R_i=rsi/(A*f_i); parallel branches recover rsi/A.
inline double segment_resistance(double rsi, double area, double fraction, double closed) {
  if (!std::isfinite(rsi) || rsi<0 || !std::isfinite(area) || area<=0 ||
      !std::isfinite(fraction) || fraction<=0 || fraction>1 || !std::isfinite(closed) || closed<=0)
    throw std::invalid_argument("Invalid area-normalized thermal insulation");
  // Native's 1e-100 K/W short defaults overflow/ill-condition the zero-clo
  // assembled circuit. 1e-8 K/W is a numerical short, not fitted insulation:
  // its relative dry-film error is <2e-7 under the retained rest conditions.
  return std::max(rsi/area,std::max(closed,1e-8))/fraction;
}
// Convection/radiation film coefficients already have units W/(m^2 K).
// Clothing resistance is a separate series layer; it must not multiply this film.
inline double film_resistance(double area, double coefficient, double closed, double open) {
  if (!std::isfinite(area) || area<=0 || !std::isfinite(coefficient) || coefficient<0 ||
      !std::isfinite(closed) || closed<=0 || !std::isfinite(open) || open<closed)
    throw std::invalid_argument("Invalid thermal film coefficient");
  return coefficient==0 ? open : std::max(1/(area*coefficient),closed);
}
}
