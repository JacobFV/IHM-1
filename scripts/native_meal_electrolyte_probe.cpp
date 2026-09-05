// Controlled nutrition additions, actual native transport and acid-base outputs.
#include <cassert>
#include <biogears/cdm/properties/SEProperties.h>
#include <biogears/engine/BioGearsPhysiologyEngine.h>
#include <biogears/cdm/patient/actions/SEConsumeNutrients.h>
#include "native_body_ports.h"
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace biogears;
int main(int argc,char**argv){
  if(argc!=4)return 4;
  const std::string mode=argv[2];const int seconds=std::stoi(argv[3]);
  if(seconds<1||seconds>1200||(mode!="rest"&&mode!="water_only"&&mode!="na_only"&&mode!="nutrients_only"&&mode!="na_water"))return 4;
  auto bg=CreateBioGearsEngine("native_engine.log");if(!bg->LoadState(argv[1]))return 2;
  auto& engine=dynamic_cast<BioGearsEngine&>(*bg);bg->SetAutoTrackFlag(false);
  std::ofstream file("electrolytes.csv");file<<std::setprecision(17);bool first=true;
  auto sample=[&](double time){
    auto v=body_ports(engine);v["time_s"]=time;
    v["native_sid_mmol_l"]=bg->GetBloodChemistrySystem()->GetStrongIonDifference(AmountPerVolumeUnit::mmol_Per_L);
    for(const std::string name:{"Aorta","VenaCava","SmallIntestineChyme","SmallIntestineVasculature","Lymph","Bladder"}){
      const auto* c=bg->GetCompartments().GetLiquidCompartment(name);
      if(name=="Aorta"||name=="VenaCava")v[name+".ph"]=c->GetPH();
      for(const std::string substance:{"Sodium","Potassium","Chloride","Calcium","Lactate","Bicarbonate"}){
        const auto* s=bg->GetSubstanceManager().GetSubstance(substance);const auto* q=c&&s?c->GetSubstanceQuantity(*s):nullptr;
        v[name+"."+substance+".mmol_l"]=q?q->GetMolarity(AmountPerVolumeUnit::mmol_Per_L):std::numeric_limits<double>::quiet_NaN();
        v[name+"."+substance+".mass_g"]=q?q->GetMass(MassUnit::g):std::numeric_limits<double>::quiet_NaN();
      }
    }
    v["vc_sid_reconstructed_mmol_l"]=v["VenaCava.Sodium.mmol_l"]+v["VenaCava.Potassium.mmol_l"]-v["VenaCava.Chloride.mmol_l"]-v["VenaCava.Lactate.mmol_l"]-1.02;
    if(first){bool comma=false;for(const auto&[k,x]:v){if(comma)file<<',';file<<k;comma=true;}file<<'\n';first=false;}
    bool comma=false;for(const auto&[k,x]:v){if(comma)file<<',';file<<x;comma=true;}file<<'\n';file.flush();
  };
  sample(0);
  if(mode!="rest"){
    SEConsumeNutrients action;auto& n=action.GetNutrition();n.SetName(mode);
    n.GetCarbohydrate().SetValue(mode=="nutrients_only"?60:0,MassUnit::g);
    n.GetProtein().SetValue(mode=="nutrients_only"?20:0,MassUnit::g);
    n.GetFat().SetValue(mode=="nutrients_only"?20:0,MassUnit::g);
    n.GetSodium().SetValue(mode=="na_only"||mode=="na_water"?1:0,MassUnit::g);
    n.GetCalcium().SetValue(0,MassUnit::mg);
    n.GetWater().SetValue(mode=="water_only"||mode=="na_water"?500:0,VolumeUnit::mL);
    if(!bg->ProcessAction(action))return 3;
  }
  for(int tick=1;tick<=seconds*50;tick++){
    if(!bg->AdvanceModelTime())return 5;
    if(tick==1||tick%500==0||tick==seconds*50)sample(tick*.02);
  }
  bg->SaveStateToFile("final_state.xml");return file.good()?0:6;
}
