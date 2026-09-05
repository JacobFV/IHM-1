// Direct boundary constitutive probe. It does not advance physiology or fit a
// patient: controlled area, temperature and sweat inputs exercise PreProcess.
#include <cassert>
#include <biogears/engine/BioGearsPhysiologyEngine.h>
#include <biogears/engine/Controller/BioGearsEngine.h>
#include <biogears/engine/Controller/BioGearsCircuits.h>
#include <biogears/engine/Systems/Environment.h>
#include <biogears/cdm/circuit/thermal/SEThermalCircuit.h>
#include <biogears/cdm/circuit/thermal/SEThermalCircuitPath.h>
#include <biogears/cdm/circuit/thermal/SEThermalCircuitNode.h>
#include <biogears/cdm/system/physiology/SEEnergySystem.h>
#include <biogears/cdm/system/environment/SEEnvironmentalConditions.h>
#include <biogears/cdm/patient/SEPatient.h>
#include <biogears/cdm/properties/SEProperties.h>
#include <biogears/cdm/utils/GeneralMath.h>
#include <fstream>
#include <iomanip>
#include <vector>
using namespace biogears;
int main(int argc,char** argv){
  if(argc!=3)return 4;
  auto engine=CreateBioGearsEngine("probe.log");if(!engine->LoadState(argv[1]))return 2;
  auto& body=dynamic_cast<BioGearsEngine&>(*engine);
  auto& env=dynamic_cast<Environment&>(body.BioGears::GetEnvironment());
  // An isolated constitutive probe supplies boundary states directly. No time
  // integration occurs and this controlled state is never saved as physiology.
  body.GetCircuits().SetReadOnlyThermal(false);
  auto& circuit=body.GetCircuits().GetTemperatureCircuit();
  const char* regions[]={"ExternalTorsoSkin","ExternalHeadSkin","ExternalLeftArmSkin","ExternalRightArmSkin","ExternalLeftLegSkin","ExternalRightLegSkin"};
  std::ofstream out(argv[2]);out<<std::setprecision(17);
  out<<"area_m2,rh,pattern,sweat_kg_s,latent_j_kg,region,skin_c,evaporation_w,evap_coefficient_w_m2_kpa,skin_saturation_pa,ambient_saturation_pa\n";
  for(double area:{.7,1.,1.9012784297057,3.})for(double rh:{0.,.5,1.})for(int pattern:{0,1})for(double sweat:{0.,1e-6,.001}){
    body.BioGears::GetPatient().GetSkinSurfaceArea().SetValue(area,AreaUnit::m2);
    env.GetConditions().GetClothingResistance().SetValue(0.,HeatResistanceAreaUnit::clo);
    env.GetConditions().GetRelativeHumidity().SetValue(rh);
    env.GetConditions().GetAmbientTemperature().SetValue(22.,TemperatureUnit::C);
    circuit.GetNode("Ambient")->GetTemperature().SetValue(22.,TemperatureUnit::C);
    for(int i=0;i<6;i++)circuit.GetNode(regions[i])->GetTemperature().SetValue(pattern?24.+2*i:34.,TemperatureUnit::C);
    body.BioGears::GetEnergy().GetSweatRate().SetValue(sweat,MassPerTimeUnit::kg_Per_s);
    env.PreProcess();
    double tk=295.15,latent=(-.1004*tk*tk+22.173*tk+46375)/.0180153;
    for(int i=0;i<6;i++){
      double skin=circuit.GetNode(regions[i])->GetTemperature(TemperatureUnit::C);
      double sp=Convert(GeneralMath::AntoineEquation(skin),PressureUnit::mmHg,PressureUnit::Pa);
      double ap=Convert(GeneralMath::AntoineEquation(22.),PressureUnit::mmHg,PressureUnit::Pa);
      auto path=circuit.GetPath(std::string(regions[i])+"ToGround");
      out<<area<<','<<rh<<','<<pattern<<','<<sweat<<','<<latent<<','<<i<<','<<skin<<','<<path->GetNextHeatSource(PowerUnit::W)<<','<<env.GetEvaporativeHeatTranferCoefficient(HeatConductancePerAreaUnit::W_Per_m2_K)<<','<<sp<<','<<ap<<'\n';
    }
  }
  return out.good()?0:5;
}
