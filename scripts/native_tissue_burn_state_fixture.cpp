// Staged only: compile with the matching regenerated schema/core variant.
// Actual native owner codecs; no patient initialization, physiology or time step.
#include "native_tissue_burn_state.h"
#include <cassert>
#include <iostream>
#include <biogears/engine/Controller/BioGears.h>
#include <biogears/engine/Systems/Tissue.h>
#include <biogears/cdm/properties/SEProperties.h>
#include "io/biogears/BioGearsPhysiology.h"
class FixtureTissue: public biogears::Tissue {
public:
 explicit FixtureTissue(biogears::BioGears& bg):Tissue(bg){}
 unsigned setup_calls=0;
 void set(const ihm_burn::State& s){
  ihm_burn::validate(s);
  m_trunkDeltaResistance_mmHg_s_Per_mL=s.resistance_mmHg_s_per_mL[0];
  m_leftArmDeltaResistance_mmHg_s_Per_mL=s.resistance_mmHg_s_per_mL[1];
  m_rightArmDeltaResistance_mmHg_s_Per_mL=s.resistance_mmHg_s_per_mL[2];
  m_leftLegDeltaResistance_mmHg_s_Per_mL=s.resistance_mmHg_s_per_mL[3];
  m_rightLegDeltaResistance_mmHg_s_Per_mL=s.resistance_mmHg_s_per_mL[4];
  m_compartmentSyndromeCount=s.syndrome_count;m_baselineECFluidVolume_mL=s.baseline_ecf_mL;
  m_trunkEscharotomy=s.escharotomy_flags&1;m_leftArmEscharotomy=s.escharotomy_flags&2;
  m_rightArmEscharotomy=s.escharotomy_flags&4;m_leftLegEscharotomy=s.escharotomy_flags&8;m_rightLegEscharotomy=s.escharotomy_flags&16;
 }
 ihm_burn::State get() const{return {{m_trunkDeltaResistance_mmHg_s_Per_mL,m_leftArmDeltaResistance_mmHg_s_Per_mL,m_rightArmDeltaResistance_mmHg_s_Per_mL,m_leftLegDeltaResistance_mmHg_s_Per_mL,m_rightLegDeltaResistance_mmHg_s_Per_mL},m_compartmentSyndromeCount,m_baselineECFluidVolume_mL,(m_trunkEscharotomy?1u:0u)|(m_leftArmEscharotomy?2u:0u)|(m_rightArmEscharotomy?4u:0u)|(m_leftLegEscharotomy?8u:0u)|(m_rightLegEscharotomy?16u:0u)};}
 void seed_base(){
  using namespace biogears;
  m_RestingPatientMass_kg=70;m_RestingFluidMass_kg=40;
  m_O2ConsumedRunningAverage_mL_Per_s.Sample(0);m_CO2ProducedRunningAverage_mL_Per_s.Sample(0);m_RespiratoryQuotientRunningAverage.Sample(0);m_FatigueRunningAverage.Sample(0);
  GetCarbonDioxideProductionRate().SetValue(0,VolumePerTimeUnit::mL_Per_s);GetOxygenConsumptionRate().SetValue(0,VolumePerTimeUnit::mL_Per_s);
  GetDehydrationFraction().SetValue(0);GetIntracellularFluidPH().SetValue(7);GetRespiratoryExchangeRatio().SetValue(0);
  GetExtracellularFluidVolume().SetValue(10000,VolumeUnit::mL);GetExtravascularFluidVolume().SetValue(40000,VolumeUnit::mL);GetIntracellularFluidVolume().SetValue(30000,VolumeUnit::mL);GetTotalBodyFluidVolume().SetValue(40000,VolumeUnit::mL);
  GetLiverInsulinSetPoint().SetValue(1,AmountPerVolumeUnit::mmol_Per_L);GetMuscleInsulinSetPoint().SetValue(1,AmountPerVolumeUnit::mmol_Per_L);GetFatInsulinSetPoint().SetValue(1,AmountPerVolumeUnit::mmol_Per_L);
  GetLiverGlucagonSetPoint().SetValue(1,MassPerVolumeUnit::g_Per_L);GetMuscleGlucagonSetPoint().SetValue(1,MassPerVolumeUnit::g_Per_L);GetFatGlucagonSetPoint().SetValue(1,MassPerVolumeUnit::g_Per_L);
  GetLiverGlycogen().SetValue(1,MassUnit::g);GetMuscleGlycogen().SetValue(1,MassUnit::g);GetStoredProtein().SetValue(1,MassUnit::g);GetStoredFat().SetValue(1,MassUnit::g);
 }
protected:
 // Deliberately overwrite to prove the codec restores AFTER preparation.
 // This is not validation of full native compartment/circuit pointer setup.
 void SetUp() override{++setup_calls;set(ihm_burn::fresh());}
};
int main(){
 using namespace biogears;
 BioGears bg("tissue-burn-codec.log");FixtureTissue source(bg),target(bg);
 const ihm_burn::State history{{.1,.2,.3,.4,.5},7,12345.678901234567,21};
 source.seed_base();source.set(history);CDM::BioGearsTissueSystemData saved;
 io::BiogearsPhysiology::Marshall(source,saved);assert(saved.IHMBurnHistory().get()==ihm_burn::encode(history));
 target.set(ihm_burn::fresh());io::BiogearsPhysiology::UnMarshall(saved,bg.GetSubstances(),target);
 assert(target.setup_calls==1);assert(ihm_burn::encode(target.get())==ihm_burn::encode(history));
 target.set(ihm_burn::fresh());io::BiogearsPhysiology::UnMarshall(saved,bg.GetSubstances(),target);assert(target.setup_calls==2);assert(ihm_burn::encode(target.get())==ihm_burn::encode(history));
 io::BiogearsPhysiology::UnMarshall(saved,bg.GetSubstances(),source);assert(ihm_burn::encode(source.get())==ihm_burn::encode(history));
 // Missing or malformed history must fail before invalidating a live owner.
 saved.IHMBurnHistory("IHM_BURN_V1:broken");bool rejected=false;
 try{io::BiogearsPhysiology::UnMarshall(saved,bg.GetSubstances(),target);}catch(const std::exception&){rejected=true;}
 assert(rejected&&target.setup_calls==2&&ihm_burn::encode(target.get())==ihm_burn::encode(history));
 CDM::BioGearsTissueSystemData missing;rejected=false;
 try{io::BiogearsPhysiology::UnMarshall(missing,bg.GetSubstances(),target);}catch(const std::exception&){rejected=true;}
 assert(rejected&&target.setup_calls==2&&ihm_burn::encode(target.get())==ihm_burn::encode(history));
 std::cout<<"PASS native Tissue owner codec; all12 histories; fresh/new/same owner; setup ordering; malformed/missing rejection. No whole-engine XML replay or dynamics tested.\n";
}
