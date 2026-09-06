// Reuse the actual signed adapter and native processes; extend observation only.
#include <cassert>
#include <biogears/engine/BioGearsPhysiologyEngine.h>
#include "native_body_ports.h"
#include "native_nervous_afferents.h"
inline std::map<std::string,double> afferent_body_ports(biogears::BioGearsEngine& engine){
  auto values=body_ports(engine);auto afferents=ihm_afferents::observe(*engine.GetNervousSystem());
  values.insert(afferents.begin(),afferents.end());return values;
}
#define body_ports afferent_body_ports
#include "native_biogears_signed.cpp"
#undef body_ports
