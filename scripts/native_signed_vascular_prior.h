#pragma once
// Experimental transfer of the native exercise resistance prior. No energy owner.
#include <cmath>
#include <stdexcept>
namespace ihm_signed_vascular {
struct Trace { double sequence=0,background_W=0,basal_W=0,delta_W=0,ratio=1,muscle_before=0,muscle_after=0,other_before=0,other_after=0; };
extern "C" Trace* ihm_signed_vascular_trace();
struct Factors { double muscle,other,ratio; };
inline Factors response(double background_W,double basal_W,double signed_delta_W) {
 if(!std::isfinite(background_W)||!std::isfinite(basal_W)||!std::isfinite(signed_delta_W)||background_W<0||basal_W<=0||background_W+signed_delta_W<0)
  throw std::runtime_error("signed vascular prior outside nonnegative effective-demand/positive basal domain");
 if(signed_delta_W==0)return {1,1,1};
 // q(F)=(1.5 F+5.5)/7 from Cardiovascular::MetabolicToneResponse.
 // Normalize against current native background, not zero absolute metabolism.
 const double ratio=1+1.5*signed_delta_W/(1.5*background_W+5.5*basal_W);
 if(!std::isfinite(ratio)||ratio<=0)throw std::runtime_error("invalid signed vascular prior ratio");
 return {1/(ratio*ratio),1/ratio,ratio};
}
}
