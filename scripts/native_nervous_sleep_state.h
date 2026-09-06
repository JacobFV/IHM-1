// Versioned exact sleep-owner payload for an isolated native schema correction.
#pragma once
#include <array>
#include <bit>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>
namespace ihm_sleep {
static_assert(sizeof(double)==8&&std::numeric_limits<double>::is_iec559);
struct State {
 double attention_lapses,biological_debt,reaction_time_s,tired_time_hr,wake_time_min,sleep_time_min;
 unsigned sleep_state; // Wire convention: 0 Awake, 1 Sleeping, independent of native enum encoding.
 std::array<double,6> values()const{return {attention_lapses,biological_debt,reaction_time_s,tired_time_hr,wake_time_min,sleep_time_min};}
};
inline void validate(const State& s){
 for(double x:s.values())if(!std::isfinite(x)||x<0)throw std::runtime_error("Invalid/uninitialized native sleep state");
 if(s.reaction_time_s<=0||s.sleep_time_min<=0||s.sleep_state>1)throw std::runtime_error("Invalid native sleep state domain");
}
inline State fresh(double sleep_amount_min){State s{3,0,.3,0,0,sleep_amount_min,0};validate(s);return s;}
inline std::string encode(const State& s){
 validate(s);std::ostringstream out;out.imbue(std::locale::classic());out<<"IHM_SLEEP_V1:"<<s.sleep_state<<std::hex<<std::setfill('0');
 for(double x:s.values())out<<':'<<std::setw(16)<<std::bit_cast<std::uint64_t>(x);
 return out.str();
}
inline State decode(const std::string& text){
 const std::string prefix="IHM_SLEEP_V1:";
 if(!text.starts_with(prefix)||text.size()!=prefix.size()+1+6*17)throw std::runtime_error("Missing/unknown native sleep state version or length");
 unsigned mode=text[prefix.size()]-'0';std::array<double,6> v{};size_t offset=prefix.size()+1;
 for(double& x:v){if(text[offset++]!=':')throw std::runtime_error("Invalid sleep state separator");std::uint64_t bits=0;
  const char* first=text.data()+offset;auto parsed=std::from_chars(first,first+16,bits,16);
  if(parsed.ec!=std::errc{}||parsed.ptr!=first+16)throw std::runtime_error("Invalid sleep state bits");
  x=std::bit_cast<double>(bits);offset+=16;
 }
 State s{v[0],v[1],v[2],v[3],v[4],v[5],mode};validate(s);
 if(encode(s)!=text)throw std::runtime_error("Noncanonical native sleep snapshot");
 return s;
}
// Held XSD serializer uses digits10=15, without the max-precision macro.
// The authoritative versioned payload is exact; public XML mirrors may have that
// specific formatting loss. Reject any other mismatch instead of hiding it.
inline bool public_compatible(double actual,double exact){
 if(actual==exact)return true;
 std::ostringstream text;text.imbue(std::locale::classic());text<<std::setprecision(std::numeric_limits<double>::digits10)<<exact;
 std::istringstream input(text.str());input.imbue(std::locale::classic());double rounded=0;input>>rounded;
 return input&&actual==rounded;
}
inline void require_public(double actual,double exact){if(!public_compatible(actual,exact))throw std::runtime_error("Native sleep XML mirror conflicts with exact owner snapshot");}
}
