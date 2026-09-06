// Exact versioned owner history for isolated Tissue serialization repair.
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
namespace ihm_burn {
static_assert(sizeof(double)==8 && std::numeric_limits<double>::is_iec559);
struct State {
 // Region order is trunk, left arm, right arm, left leg, right leg.
 std::array<double,5> resistance_mmHg_s_per_mL;
 double syndrome_count; // Existing native counter, not recomputed event count.
 double baseline_ecf_mL; // Zero only before first burn, not current ECF.
 unsigned escharotomy_flags; // Bits 0..4 in the same region order.
 std::array<double,7> values()const{return {resistance_mmHg_s_per_mL[0],resistance_mmHg_s_per_mL[1],resistance_mmHg_s_per_mL[2],resistance_mmHg_s_per_mL[3],resistance_mmHg_s_per_mL[4],syndrome_count,baseline_ecf_mL};}
};
inline void validate(const State& s){
 for(double x:s.values())if(!std::isfinite(x)||x<0)throw std::runtime_error("Invalid/uninitialized Tissue burn history");
 if(s.escharotomy_flags>31||std::floor(s.syndrome_count)!=s.syndrome_count||s.syndrome_count>9007199254740992.)throw std::runtime_error("Invalid Tissue burn flags/count");
 bool evolved=s.syndrome_count>0||s.escharotomy_flags!=0;
 for(double x:s.resistance_mmHg_s_per_mL)evolved=evolved||x>0;
 if(evolved&&s.baseline_ecf_mL<=0)throw std::runtime_error("Evolved burn history requires historical ECF baseline");
}
inline State fresh(){return State{{0,0,0,0,0},0,0,0};}
inline std::string encode(const State& s){
 validate(s);std::ostringstream out;out.imbue(std::locale::classic());
 out<<"IHM_BURN_V1:"<<std::hex<<std::setfill('0')<<std::setw(2)<<s.escharotomy_flags;
 for(double x:s.values())out<<':'<<std::setw(16)<<std::bit_cast<std::uint64_t>(x);
 return out.str();
}
inline State decode(const std::string& text){
 const std::string prefix="IHM_BURN_V1:";
 if(!text.starts_with(prefix)||text.size()!=prefix.size()+2+7*17)throw std::runtime_error("Missing/unknown Tissue burn history version or length");
 unsigned flags=0;const char* flag=text.data()+prefix.size();auto fp=std::from_chars(flag,flag+2,flags,16);
 if(fp.ec!=std::errc{}||fp.ptr!=flag+2)throw std::runtime_error("Invalid burn flags");
 std::array<double,7> values{};size_t offset=prefix.size()+2;
 for(double& x:values){if(text[offset++]!=':')throw std::runtime_error("Invalid burn separator");
  std::uint64_t bits=0;const char* first=text.data()+offset;auto p=std::from_chars(first,first+16,bits,16);
  if(p.ec!=std::errc{}||p.ptr!=first+16)throw std::runtime_error("Invalid burn history bits");
  x=std::bit_cast<double>(bits);offset+=16;
 }
 State s{{values[0],values[1],values[2],values[3],values[4]},values[5],values[6],flags};
 validate(s);if(encode(s)!=text)throw std::runtime_error("Noncanonical Tissue burn history");return s;
}
} // namespace ihm_burn
