#include <cassert>
#include "Circuit.h"
#include <biogears/cdm/circuit/SECircuitManager.h>
#include <biogears/cdm/circuit/fluid/SEFluidCircuitNode.h>
#include <biogears/cdm/circuit/fluid/SEFluidCircuitPath.h>
#include <iostream>
using namespace biogears;
int main(){try{
 SECircuitManager original(nullptr),restored(nullptr);auto& a=original.CreateFluidNode("a");auto& b=original.CreateFluidNode("b");
 const SEResistancePathType regions[]={SEResistancePathType::Cerebral,SEResistancePathType::Extrasplanchnic,SEResistancePathType::Muscle,SEResistancePathType::Splanchnic,SEResistancePathType::Myocardium};
 for(int i=0;i<5;++i){auto& p=original.CreateFluidPath(a,b,"p"+std::to_string(i));p.SetCardiovascularRegion(regions[i]);p.GetResistanceBaseline().SetValue(1.25+i,FlowResistanceUnit::mmHg_s_Per_mL);}
 CDM::CircuitManagerData record;io::Circuit::Marshall(original,record);io::Circuit::UnMarshall(record,restored);
 for(int i=0;i<5;++i){auto* p=restored.GetFluidPath("p"+std::to_string(i));if(!p||!p->HasCardiovascularRegion()||p->GetCardiovascularRegion()!=regions[i]||p->GetResistanceBaseline(FlowResistanceUnit::mmHg_s_Per_mL)!=1.25+i)throw std::runtime_error("region roundtrip failed");}
 auto* p=original.GetFluidPath("p0");CDM::FluidCircuitPathData repeated;io::Circuit::Marshall(*p,repeated);if(!repeated.CardiovascularRegion().present())throw std::runtime_error("region absent on write");p->InvalidateCardiovascularRegion();io::Circuit::Marshall(*p,repeated);if(repeated.CardiovascularRegion().present())throw std::runtime_error("stale serialized region");
 for(int i=0;i<5;++i)original.GetFluidPath("p"+std::to_string(i))->InvalidateCardiovascularRegion();CDM::CircuitManagerData absent;io::Circuit::Marshall(original,absent);io::Circuit::UnMarshall(absent,restored);
 for(int i=0;i<5;++i)if(restored.GetFluidPath("p"+std::to_string(i))->HasCardiovascularRegion())throw std::runtime_error("stale native region");
 std::cout<<"{\"passed\":true,\"regions\":5,\"stale_field_clear\":true,\"resistance_baselines_preserved\":true}"<<std::endl;return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
