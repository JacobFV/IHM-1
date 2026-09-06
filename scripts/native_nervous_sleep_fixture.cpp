// Isolated actual Nervous owner IO/calculation acceptance; no whole-engine advance.
#include <cassert>
#include "native_nervous_sleep_state.h"
#include <biogears/engine/Controller/BioGears.h>
#include <biogears/engine/Controller/BioGearsSubstances.h>
#include <biogears/engine/Systems/Nervous.h>
#include <biogears/cdm/patient/SEPatient.h>
#include <biogears/cdm/properties/SEScalarTime.h>
#include <biogears/schema/biogears/BioGears.hxx>
#include "io/Serializer.h"
#include "io/biogears/BioGearsPhysiology.h"
#include <iostream>
#include <fstream>
#include <sstream>
using namespace biogears;
struct ClockOwner : BioGears {
 explicit ClockOwner(const std::string& runtime):BioGears("sleep-owner.log",runtime){}
 void clock(double seconds){m_SimulationTime->SetValue(seconds,TimeUnit::s);}
};
struct SleepOwner : Nervous {
 explicit SleepOwner(BioGears& bg):Nervous(bg){}
 ihm_sleep::State snapshot(){
  if(GetSleepState()!=SESleepState::Awake&&GetSleepState()!=SESleepState::Sleeping)throw std::runtime_error("Invalid fixture sleep mode");
  return {m_AttentionLapses,m_BiologicalDebt,m_ReactionTime_s,m_TiredTime_hr,GetWakeTime().GetValue(TimeUnit::min),GetSleepTime().GetValue(TimeUnit::min),GetSleepState()==SESleepState::Sleeping?1u:0u};
 }
 void calculate(){CalculateSleepEffects();}
 void extreme_debt(){m_BiologicalDebt=std::numeric_limits<double>::max();GetBiologicalDebt().SetValue(m_BiologicalDebt);}
};
CDM::BioGearsNervousSystemData& nervous(CDM::ObjectData& data){
 auto& state=dynamic_cast<CDM::BioGearsStateData&>(data);
 CDM::BioGearsNervousSystemData* found=nullptr;
 for(auto& system:state.System())if(auto* candidate=dynamic_cast<CDM::BioGearsNervousSystemData*>(&system)){
  if(found)throw std::runtime_error("Duplicate Nervous owner");found=candidate;
 }
 if(!found)throw std::runtime_error("Missing Nervous owner");return *found;
}
template<class F>void rejects(F action){bool rejected=false;try{action();}catch(const std::exception&){rejected=true;}if(!rejected)throw std::runtime_error("Expected rejection was absent");}
void require(bool condition,const char* reason){if(!condition)throw std::runtime_error(reason);}
int main(int argc,char** argv){try{
 if(argc!=6)throw std::runtime_error("Expected runtime legacy-state explicit-seed output-XML patient-sleep-minutes");
 ClockOwner bg(argv[1]);bg.GetPatient().GetSleepAmount().SetValue(std::stod(argv[5]),TimeUnit::min);
 require(bg.GetTimeStep().GetValue(TimeUnit::s)==.02,"Unexpected native timestep");
 auto legacy=biogears::Serializer::ReadFile(argv[2],bg.GetLogger());auto state=biogears::Serializer::ReadFile(argv[3],bg.GetLogger());
 require(bool(legacy)&&bool(state),"Native XML parse failed");auto& initial=nervous(*state);
 SleepOwner direct(bg),resumed(bg),uninitialized(bg);CDM::BioGearsNervousSystemData invalid;
 rejects([&]{io::BiogearsPhysiology::Marshall(uninitialized,invalid);});
 rejects([&]{io::BiogearsPhysiology::UnMarshall(nervous(*legacy),bg.GetSubstances(),uninitialized);});
 auto expected=ihm_sleep::decode(initial.IHMSleepState().get());
 io::BiogearsPhysiology::UnMarshall(initial,bg.GetSubstances(),direct);
 io::BiogearsPhysiology::UnMarshall(initial,bg.GetSubstances(),resumed);
 require(ihm_sleep::encode(direct.snapshot())==ihm_sleep::encode(expected),"Initial native owner differs from complete explicit seed");
 rejects([&]{io::BiogearsPhysiology::UnMarshall(nervous(*legacy),bg.GetSubstances(),direct);});
 require(ihm_sleep::encode(direct.snapshot())==ihm_sleep::encode(expected),"Missing-version preflight touched initialized sleep owner");
 SleepOwner extreme(bg);io::BiogearsPhysiology::UnMarshall(initial,bg.GetSubstances(),extreme);
 extreme.extreme_debt();ihm_sleep::validate(extreme.snapshot());
 bool overflow_rejected=false;try{extreme.calculate();}catch(const std::exception& e){overflow_rejected=std::string(e.what())=="Invalid/uninitialized native sleep state";}
 require(overflow_rejected,"Finite native ODE overflow was not rejected at output validation"); // Discard failed owner; no rollback.
 auto roundtrip=[&]{
  CDM::BioGearsNervousSystemData destination;
  io::BiogearsPhysiology::Marshall(resumed,destination);
  nervous(*state)=destination; // Fresh destination: native marshaller appends effector vectors.
  std::ostringstream xml;xml_schema::namespace_infomap map;map[""].name="uri:/mil/tatrc/physiology/datamodel";
  CDM::BioGearsState(xml,dynamic_cast<CDM::BioGearsStateData&>(*state),map);
  auto bytes=xml.str();auto restored=biogears::Serializer::ReadBuffer(reinterpret_cast<const XMLByte*>(bytes.data()),bytes.size(),bg.GetLogger());
  require(bool(restored),"Native saved XML reparse failed");
  io::BiogearsPhysiology::UnMarshall(nervous(*restored),bg.GetSubstances(),resumed);state=std::move(restored);
  require(ihm_sleep::encode(direct.snapshot())==ihm_sleep::encode(resumed.snapshot()),"Native sleep XML continuation bit mismatch");
  return bytes;
 };
 std::string saved=roundtrip(); // save immediately after load, before calculation
 for(int i=0;i<12;++i){bg.clock(i*.02);direct.calculate();resumed.calculate();saved=roundtrip();}
 auto final=direct.snapshot();ihm_sleep::validate(final);
 if(expected.wake_time_min/expected.sleep_time_min>3&&expected.sleep_state==0){
  require(final.tired_time_hr>expected.tired_time_hr,"Awake tired history failed to evolve");
  require(final.reaction_time_s>.3&&final.reaction_time_s<1,"Reaction time seconds contract failed");
 }else if(expected.sleep_state==0){require(final.reaction_time_s==expected.reaction_time_s,"Rested branch changed reaction-time units");}
 // Native loader must reject both future versions and conflicting public mirrors.
 auto valid_payload=nervous(*state).IHMSleepState().get();nervous(*state).IHMSleepState("IHM_SLEEP_V2:unsupported");
 rejects([&]{io::BiogearsPhysiology::UnMarshall(nervous(*state),bg.GetSubstances(),uninitialized);});
 nervous(*state).IHMSleepState(valid_payload);nervous(*state).AttentionLapses(9000.);
 rejects([&]{io::BiogearsPhysiology::UnMarshall(nervous(*state),bg.GetSubstances(),uninitialized);});
 std::ofstream output(argv[4]);output<<saved;require(bool(output),"Fixture XML write failed");
 std::cout<<"RESULT {\"passed\":true,\"native_sleep_calculations_per_continuation_owner\":12,\"native_xml_roundtrips\":13,\"whole_engine_advances\":0,\"finite_extreme_overflow_rejected\":true,\"overflow_calculation_attempts\":1,\"sleep_owner_exact\":\""<<ihm_sleep::encode(final)<<"\"}\n";
 return 0;
}catch(const std::exception& e){std::cerr<<"SLEEP_OWNER_FAILURE "<<e.what()<<'\n';return 1;}}
