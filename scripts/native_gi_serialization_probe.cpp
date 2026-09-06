#include <cassert>
#include <cstdlib>
#include <iostream>
#include <biogears/engine/Controller/BioGears.h>
#include <biogears/engine/Systems/Gastrointestinal.h>
#include <biogears/cdm/properties/SEProperties.h>
#include "io/biogears/BioGearsPhysiology.h"
#include "io/cdm/Physiology.h"
static void* watched=nullptr;static int destroyed=0;
void operator delete(void* p)noexcept{if(p==watched){++destroyed;watched=nullptr;}std::free(p);}
void operator delete(void* p,std::size_t)noexcept{if(p==watched){++destroyed;watched=nullptr;}std::free(p);}
class FixtureGI:public biogears::Gastrointestinal {
public:explicit FixtureGI(biogears::BioGears& bg):Gastrointestinal(bg){}
protected:void SetUp()override{} // Only preparation hook omitted; actual codecs run.
};
int main(int argc,char**argv){
 using namespace biogears;std::string mode=argc>1?argv[1]:"roundtrip";
 BioGears bg("codec.log");if(!bg.GetSubstances().LoadSubstanceDirectory())return 2;
 auto* a=bg.GetSubstances().GetSubstances()[0];auto* b=bg.GetSubstances().GetSubstances()[1];
 FixtureGI source(bg);auto* state=source.NewDrugTransitState(a);
 std::vector<double> solid(9),dissolved(9),enterocyte(8);for(int i=0;i<9;i++){solid[i]=i+1;dissolved[i]=20+i;}for(int i=0;i<8;i++)enterocyte[i]=40+i;
 assert(state->SetLumenSolidMasses(solid,MassUnit::ug));assert(state->SetLumenDissolvedMasses(dissolved,MassUnit::ug));assert(state->SetEnterocyteMasses(enterocyte,MassUnit::ug));
 state->GetTotalMassMetabolized().SetValue(123,MassUnit::ug);state->GetTotalMassExcreted().SetValue(456,MassUnit::ug);
 if(mode=="nested") {CDM::DrugTransitStateData encoded;io::Physiology::Marshall(*state,encoded);SEDrugTransitState restored(*a);io::Physiology::UnMarshall(encoded,restored);bool good=restored.GetTotalMassExcreted().GetValue(MassUnit::ug)==456;std::cout<<"NESTED "<<good<<std::endl;return good?0:3;}
 if(mode=="cleanup") {watched=state;source.Invalidate();bool good=source.GetDrugTransitStates().empty()&&destroyed==1;watched=nullptr;std::cout<<"CLEANUP "<<good<<std::endl;return good?0:3;}
 CDM::BioGearsGastrointestinalSystemData encoded;io::BiogearsPhysiology::Marshall(source,encoded);
 if(encoded.DrugTransitStates().size()!=1){std::cout<<"OUTER_OMITTED "<<encoded.DrugTransitStates().size()<<std::endl;return 3;}
 io::BiogearsPhysiology::Marshall(source,encoded);assert(encoded.DrugTransitStates().size()==1);
 FixtureGI target(bg);auto* old=target.NewDrugTransitState(b);watched=old;destroyed=0;
 io::BiogearsPhysiology::UnMarshall(encoded,bg.GetSubstances(),target);
 assert(destroyed==1);watched=nullptr;assert(target.GetDrugTransitStates().size()==1);auto* restored=target.GetDrugTransitState(a);
 assert(restored->GetLumenSolidMasses(MassUnit::ug)==solid);assert(restored->GetLumenDissolvedMasses(MassUnit::ug)==dissolved);assert(restored->GetEnterocyteMasses(MassUnit::ug)==enterocyte);
 assert(restored->GetTotalMassMetabolized().GetValue(MassUnit::ug)==123);assert(restored->GetTotalMassExcreted().GetValue(MassUnit::ug)==456);
 watched=restored;destroyed=0;target.NewDrugTransitState(a);assert(destroyed==1);watched=nullptr;
 // Missing historical records are absent, never reconstructed as fake zero states.
 CDM::BioGearsGastrointestinalSystemData historical;io::BiogearsPhysiology::UnMarshall(historical,bg.GetSubstances(),target);assert(target.GetDrugTransitStates().empty());
 std::cout<<"PASS vectors=9,9,8 metabolized_ug=123 excreted_ug=456 cleanup=owned replacement=owned historical=absent\n";
}
