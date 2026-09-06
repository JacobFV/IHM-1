// Public paired-DSO acceptance wrapper; original respiratory probe remains frozen.
#define main frozen_respiratory_probe_main
#include "native_respiratory_engine_probe.cpp"
#undef main
#include <fstream>
int main(int argc,char** argv){
 std::ifstream maps("/proc/self/maps");std::string line;bool core=false,cdm=false;
 while(std::getline(maps,line)){
  const bool is_core=line.find("/libbiogears.so.8.0.0")!=std::string::npos;
  const bool is_cdm=line.find("/libbiogears_cdm.so.8.0.0")!=std::string::npos;
  if(is_core||is_cdm){std::cout<<"LOADED_LIBRARY "<<line<<'\n';core|=is_core;cdm|=is_cdm;}
 }
 if(!core||!cdm){std::cerr<<"Missing actual paired DSO mapping\n";return 2;}
 return frozen_respiratory_probe_main(argc,argv);
}
