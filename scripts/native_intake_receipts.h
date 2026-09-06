// Observational external-intake ledger. No tissue/GI stores or native mass writers.
#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <optional>
#include <ostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unistd.h>
namespace ihm_intake {
inline const std::array<std::string,6> nutrient_names{"carbohydrate_g","protein_g","fat_g","sodium_g","calcium_mg","water_ml"};
inline bool token(const std::string& s) {
 return !s.empty()&&std::all_of(s.begin(),s.end(),[](unsigned char c){return (c>='a'&&c<='z')||(c>='A'&&c<='Z')||(c>='0'&&c<='9')||c=='_'||c=='-';});
}
inline std::string epoch() {
 std::random_device random;std::ostringstream out;out<<getpid()<<'-'<<std::hex;
 for(unsigned i=0;i<4;++i){out<<std::setw(8)<<std::setfill('0')<<random();}
 return out.str();
}
struct Payload {
 std::string name;
 double mass_kg=0;
 std::array<double,6> nutrients{};
 bool operator==(const Payload&) const=default;
 void validate() const {
  if(!token(name)||!std::isfinite(mass_kg)||mass_kg<=0)throw std::runtime_error("Invalid held native intake payload");
  for(double x:nutrients)if(!std::isfinite(x)||x<0)throw std::runtime_error("Invalid held native nutrient");
 }
};
// Units are passed explicitly from the native caller; this helper never estimates
// weight from configured amounts or from a stomach/vascular difference.
template<class Nutrition,class MassUnit,class VolumeUnit>
Payload capture(const Nutrition& n,const MassUnit& kg,const MassUnit& g,const MassUnit& mg,const VolumeUnit& ml) {
 if(!n.HasCarbohydrate()||!n.HasProtein()||!n.HasFat()||!n.HasSodium()||!n.HasCalcium()||!n.HasWater())
  throw std::runtime_error("Held native intake is missing explicit nutrient quantities");
 Payload p{n.GetName(),n.GetWeight(kg),{n.GetCarbohydrate(g),n.GetProtein(g),n.GetFat(g),n.GetSodium(g),n.GetCalcium(mg),n.GetWater(ml)}};
 p.validate();return p;
}
inline void quantities(std::ostream& out,const std::array<double,6>& amounts) {
 out<<'{';for(size_t i=0;i<6;++i){if(i)out<<',';out<<std::quoted(nutrient_names[i])<<':'<<amounts[i];}out<<'}';
}
inline void payload_json(std::ostream& out,const Payload& p) {
 out<<"{\"name\":"<<std::quoted(p.name)<<",\"mass_kg\":"<<p.mass_kg;
 for(size_t i=0;i<6;++i){out<<','<<std::quoted(nutrient_names[i])<<':'<<p.nutrients[i];}
 out<<'}';
}
struct Pending {std::uint64_t meal_sequence;Payload payload;};
struct Attempt {std::uint64_t advance_sequence,start_tick;double start_s;std::optional<Pending> pending;};
struct Consumed {std::uint64_t count,meal_sequence,advance_sequence,start_tick,end_tick;double start_s,end_s;Payload payload;};
class Ledger {
 std::string owner;
 std::uint64_t count=0,last_meal_sequence=0;
 double mass_kg=0;
 std::array<double,6> totals{};
 std::optional<Pending> pending;
 std::optional<Attempt> attempt;
 std::optional<Consumed> last;
 bool poisoned=false;
 [[noreturn]] void fail(const std::string& error){poisoned=true;throw std::runtime_error("Native intake receipt: "+error);}
 void healthy(){if(poisoned)throw std::runtime_error("Native intake receipt owner poisoned");}
public:
 explicit Ledger(std::string owner_epoch):owner(std::move(owner_epoch)){if(!token(owner))throw std::invalid_argument("Invalid intake owner epoch");}
 void accepted(std::uint64_t meal_sequence,const Payload& held) {
  healthy();held.validate();
  if(pending||attempt||meal_sequence==0||meal_sequence<=last_meal_sequence)fail("duplicate or overlapping accepted meal identity");
  pending=Pending{meal_sequence,held};last_meal_sequence=meal_sequence;
 }
 void before(const std::optional<Payload>& held,std::uint64_t advance_sequence,std::uint64_t start_tick,double start_s) {
  healthy();
  if(attempt||advance_sequence==0||!std::isfinite(start_s)||start_s<0)fail("invalid advance identity");
  if(bool(held)!=bool(pending))fail("held native action has no matching accepted identity");
  if(held){held->validate();if(*held!=pending->payload||advance_sequence<=pending->meal_sequence)fail("held native action changed before consumption");}
  attempt=Attempt{advance_sequence,start_tick,start_s,pending};
 }
 void commit(bool native_action_pending,std::uint64_t end_tick,double end_s) {
  healthy();
  if(!attempt||attempt->start_tick==std::numeric_limits<std::uint64_t>::max()||end_tick!=attempt->start_tick+1||
     !std::isfinite(end_s)||end_s<=attempt->start_s)fail("invalid successful native interval");
  if(native_action_pending)fail("unexpected native action after successful advance");
  if(attempt->pending) {
   if(count==std::numeric_limits<std::uint64_t>::max())fail("consumed count overflow");
   const auto& a=*attempt;const auto& p=a.pending->payload;
   const double next_mass=mass_kg+p.mass_kg;auto next_totals=totals;
   if(!std::isfinite(next_mass))fail("cumulative native intake mass overflow");
   for(size_t i=0;i<6;++i){next_totals[i]+=p.nutrients[i];if(!std::isfinite(next_totals[i]))fail("cumulative native nutrient overflow");}
   Consumed receipt{count+1,a.pending->meal_sequence,a.advance_sequence,a.start_tick,end_tick,a.start_s,end_s,p};
   last=std::move(receipt);++count;mass_kg=next_mass;totals=next_totals;pending.reset();
  }
  attempt.reset();
 }
 void terminal_failure(){poisoned=true;}
 void write_attempt(std::ostream& out) const {
  if(!attempt)throw std::runtime_error("No native intake attempt");
  out<<std::setprecision(17)<<"{\"owner_epoch\":"<<std::quoted(owner)<<",\"advance_sequence\":"<<attempt->advance_sequence
     <<",\"interval_start_tick\":"<<attempt->start_tick<<",\"native_start_s\":"<<attempt->start_s<<",\"pending\":";
  write_pending(out,attempt->pending);out<<'}';
 }
 static void write_pending(std::ostream& out,const std::optional<Pending>& item) {
  if(!item){out<<"null";return;}out<<"{\"meal_sequence\":"<<item->meal_sequence<<",\"native_payload\":";
  payload_json(out,item->payload);out<<'}';
 }
 void write_json(std::ostream& out) const {
  if(poisoned)throw std::runtime_error("No successful observation from poisoned intake owner");
  out<<std::setprecision(17)<<"{\"schema\":\"ihm.native-consumed-intake.v1\",\"owner_epoch\":"<<std::quoted(owner)
     <<",\"consumed_count\":"<<count<<",\"cumulative_mass_kg\":"<<mass_kg<<",\"cumulative_native\":";
  quantities(out,totals);out<<",\"pending\":";write_pending(out,pending);out<<",\"last_consumed\":";
  if(!last)out<<"null";
  else {const auto& r=*last;out<<"{\"consumed_count\":"<<r.count<<",\"meal_sequence\":"<<r.meal_sequence
      <<",\"advance_sequence\":"<<r.advance_sequence<<",\"interval_start_tick\":"<<r.start_tick<<",\"interval_end_tick\":"<<r.end_tick
      <<",\"native_start_s\":"<<r.start_s<<",\"native_end_s\":"<<r.end_s<<",\"native_payload\":";payload_json(out,r.payload);out<<'}';}
  out<<'}';
 }
};
}
