// Read-only source-internal observer. No Nervous object construction or mutation.
#pragma once
#include <biogears/engine/Systems/Nervous.h>
#include <map>
#include <string>
#include <cmath>
#include <stdexcept>
namespace ihm_afferents {
class Access : public biogears::Nervous {
public:
  // Standard protected member access produces a BASE member pointer. The actual
  // engine-owned Nervous object is never cast to this uninstantiated subclass.
  static double biogears::Nervous::* chemo(){return &Access::m_AfferentChemoreceptor_Hz;}
  static double biogears::Nervous::* carotid(){return &Access::m_AfferentBaroreceptorCarotid_Hz;}
  static double biogears::Nervous::* aortic(){return &Access::m_AfferentBaroreceptorAortic_Hz;}
  static double biogears::Nervous::* pulmonary(){return &Access::m_AfferentPulmonaryStretchReceptor_Hz;}
};
inline std::map<std::string,double> observe(const biogears::SENervousSystem& system){
  const auto* native=dynamic_cast<const biogears::Nervous*>(&system);
  if(!native)throw std::runtime_error("Native Nervous afferent owner unavailable");
  std::map<std::string,double> values{{"nervous.afferent.chemoreceptor_hz",native->*Access::chemo()},
    {"nervous.afferent.baroreceptor_carotid_hz",native->*Access::carotid()},
    {"nervous.afferent.baroreceptor_aortic_hz",native->*Access::aortic()},
    {"nervous.afferent.pulmonary_stretch_hz",native->*Access::pulmonary()}};
  for(const auto& entry:values)if(!std::isfinite(entry.second)||entry.second<0||entry.second>1000)
    throw std::runtime_error("Native afferent rate outside explicit Hz domain");
  return values;
}
}
