// Public paired-DSO continuation challenge. No internal engine lifecycle override.
#include "native_gi_persistence_assertions.h"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>

namespace {
using biogears::BioGearsEngine;
void mappings() {
  std::ifstream maps("/proc/self/maps");
  std::string line;
  bool core=false, cdm=false;
  while (std::getline(maps,line)) {
    const bool a=line.find("/libbiogears.so.8.0.0")!=std::string::npos;
    const bool b=line.find("/libbiogears_cdm.so.8.0.0")!=std::string::npos;
    if (a || b) { std::cout<<"LOADED_LIBRARY "<<line<<'\n'; core|=a; cdm|=b; }
  }
  if (!core || !cdm) throw std::runtime_error("Missing paired DSO mappings");
}
std::unique_ptr<BioGearsEngine> load(const std::string& state,const std::string& log) {
  auto e=std::make_unique<BioGearsEngine>(log);
  if (!e->LoadState(state)) throw std::runtime_error("Native continuation state load rejected");
  e->SetAutoTrackFlag(false);
  if (e->GetTimeStep(biogears::TimeUnit::s)!=.02)
    throw std::runtime_error("Unexpected native timestep");
  return e;
}
void advance(BioGearsEngine& e) {
  if (!e.AdvanceModelTime(false)) throw std::runtime_error("Native next tick rejected");
}
void emit(const std::string& label,BioGearsEngine& e,bool populated) {
  std::cout<<std::setprecision(17)<<"SNAPSHOT {\"label\":\""<<label
    <<"\",\"time_s\":"<<e.GetSimulationTime(biogears::TimeUnit::s)
    <<",\"gi_count\":"<<e.BioGears::GetGastrointestinal().GetDrugTransitStates().size()
    <<",\"gi\":{";
  bool first=true;
  if (populated) for (const auto& [key,value]:ihm_gi_acceptance::capture(e)) {
    if (!first) std::cout<<',';
    first=false; std::cout<<'"'<<key<<"\":"<<value;
  }
  std::cout<<"}}\n";
}
}
int main(int argc,char** argv) {
  try {
    if (argc!=4) return 2;
    const std::string mode=argv[1];
    if (mode!="empty" && mode!="populated" && mode!="reject") return 2;
    mappings();
    auto direct=load(argv[2],"continuation-direct.log");
    if (mode=="reject") {
      std::cerr<<"UNEXPECTED_HISTORY_ACCEPTANCE\n";
      return 4;
    }
    const std::filesystem::path out=argv[3];
    const bool populated=mode=="populated";
    if (populated) {
      ihm_gi_acceptance::request_tiny_oral_dose(*direct);
      advance(*direct);
      ihm_gi_acceptance::seed_distinct_transit(*direct);
    }
    emit("direct_before",*direct,populated);
    direct->SaveStateToFile((out/"saved.xml").string());
    auto resumed=load((out/"saved.xml").string(),"continuation-resumed.log");
    emit("resumed_before",*resumed,populated);
    resumed->SaveStateToFile((out/"reloaded.xml").string());
    // Comparator failures are retained as evidence; only native failure stops
    // advancement. Neither owner is retried or restored after a touched failure.
    advance(*direct);
    emit("direct_after",*direct,populated);
    direct->SaveStateToFile((out/"direct-next.xml").string());
    advance(*resumed);
    emit("resumed_after",*resumed,populated);
    resumed->SaveStateToFile((out/"resumed-next.xml").string());
    return 0;
  } catch (const std::exception& e) {
    std::cerr<<"CONTINUATION_FAILURE "<<e.what()<<'\n';
    return 1;
  }
}
