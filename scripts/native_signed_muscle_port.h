#pragma once
// Diagnostic source variant ABI. The shared library owns the active pointer;
// adapters own Record lifetime. No native scalar or substance writer lives here.
#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
namespace ihm_signed {
struct Record {
  std::uint64_t sequence=0;
  double start_s=0,end_s=0,delta_m_W=0,delta_h_W=0,delta_w_W=0;
  std::string reference_id;
  unsigned heat_count=0,tissue_count=0;
  double native_heat_W=0,effective_heat_W=0,allowed_decrement_W=0;
  double muscle_base_kcal=0,obligatory_aa_kcal=0,mandatory_anaerobic_kcal=0;
  double muscle_requested_kcal=0,unmet_muscle_kcal=0;
  std::array<double,4> effective_reader_W{};
  std::array<unsigned,4> effective_reader_count{};
};
extern "C" Record* ihm_signed_active_record(Record*,bool);
#ifdef IHM_SIGNED_MUSCLE_IMPLEMENTATION
extern "C" Record* ihm_signed_active_record(Record* replacement,bool set) {
 static thread_local Record* active=nullptr;
 Record* old=active;if(set)active=replacement;return old;
}
#endif
inline Record* active(){return ihm_signed_active_record(nullptr,false);}
struct Scope {
 Record* prior;
 explicit Scope(Record* r):prior(ihm_signed_active_record(r,true)){}
 ~Scope(){ihm_signed_active_record(prior,true);}
 Scope(const Scope&)=delete;Scope& operator=(const Scope&)=delete;
};
inline void validate(const Record& r) {
 for(double x:{r.start_s,r.end_s,r.delta_m_W,r.delta_h_W,r.delta_w_W})
  if(!std::isfinite(x))throw std::runtime_error("signed muscle nonfinite record");
 if(r.sequence==0||r.reference_id.empty()||r.end_s<=r.start_s||r.heat_count||r.tissue_count)
  throw std::runtime_error("signed muscle invalid or consumed record");
 if(std::abs(r.delta_m_W-r.delta_h_W-r.delta_w_W)>1e-10*(1+std::abs(r.delta_m_W)+std::abs(r.delta_h_W)+std::abs(r.delta_w_W)))
  throw std::runtime_error("signed muscle chemical/heat/work mismatch");
}
inline double effective(double base,double watts_to_units,unsigned reader) {
 auto* r=active();if(!r)return base;
 const double result=base+r->delta_m_W*watts_to_units;
 if(!std::isfinite(result)||result<0)throw std::runtime_error("signed muscle effective demand outside domain");
 r->effective_reader_W.at(reader)=result/watts_to_units;r->effective_reader_count.at(reader)++;
 return r->delta_m_W==0?base:result;
}
inline double heat(double base) {
 auto* r=active();if(!r)return base;
 if(r->heat_count)throw std::runtime_error("signed muscle heat consumed twice");
 r->native_heat_W=base;r->effective_heat_W=base+r->delta_h_W;
 if(!std::isfinite(r->effective_heat_W)||r->effective_heat_W<0)
  throw std::runtime_error("signed muscle absolute thermal source outside domain");
 ++r->heat_count;return r->delta_h_W==0?base:r->effective_heat_W;
}
inline double budget(double discretionary_before_aa,double aa,double mandatory,double dt,double watts_to_kcal_per_s) {
 auto* r=active();if(!r)return 0;
 if(r->tissue_count)throw std::runtime_error("signed muscle tissue consumed twice");
 if(std::abs(dt-(r->end_s-r->start_s))>1e-10)throw std::runtime_error("signed muscle native interval mismatch");
 const double delta=r->delta_m_W*dt*watts_to_kcal_per_s;
 r->muscle_base_kcal=discretionary_before_aa+mandatory;
 r->obligatory_aa_kcal=aa;r->mandatory_anaerobic_kcal=mandatory;
 r->allowed_decrement_W=(discretionary_before_aa-aa)/(dt*watts_to_kcal_per_s);
 r->muscle_requested_kcal=r->muscle_base_kcal+delta;
 if(!std::isfinite(aa)||aa<0||!std::isfinite(discretionary_before_aa)||!std::isfinite(mandatory)||mandatory<0||
    (r->delta_m_W!=0 && discretionary_before_aa+delta-aa<0))
  throw std::runtime_error("signed muscle decrement exceeds native discretionary budget: requested_W="+std::to_string(r->delta_m_W)+" allowed_decrement_W="+std::to_string(r->allowed_decrement_W));
 ++r->tissue_count;return delta;
}
inline void finish(const Record& r){if(r.heat_count!=1||r.tissue_count!=1)throw std::runtime_error("signed muscle interval not consumed exactly once");}
}
