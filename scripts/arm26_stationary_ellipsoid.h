// Bounded native-seeded short geodesic; isolated candidate, no force fudge.
#pragma once
#include "arm26_wrap_candidate.h"
namespace arm26geodesic {
using SimTK::Vec3;
struct State {Vec3 x,v;};
inline State add(State x,State d,double h){return {x.x+h*d.x,x.v+h*d.v};}
inline State rhs(State s,Vec3 inv){Vec3 n(s.x[0]*inv[0],s.x[1]*inv[1],s.x[2]*inv[2]);double vv=0;for(int i=0;i<3;++i)vv+=s.v[i]*s.v[i]*inv[i];return {s.v,-vv/(~n*n)*n};}
struct Shot {State start,end;double residual[4];std::vector<State> points;};
inline Shot shoot(const double*z,Vec3 p,Vec3 q,Vec3 axes,Vec3 y0,Vec3 e1,Vec3 e2,int steps){
 double scale=std::min({axes[0],axes[1],axes[2]});Vec3 y=y0+z[0]*e1+z[1]*e2;y/=y.norm();Vec3 x,inv;
 for(int i=0;i<3;++i){x[i]=axes[i]*y[i];inv[i]=1/(axes[i]*axes[i]);}
 Vec3 n(x[0]*inv[0],x[1]*inv[1],x[2]*inv[2]);n/=n.norm();Vec3 incoming=x-p;incoming/=incoming.norm();Vec3 tangent=incoming-(~n*incoming)*n;tangent/=tangent.norm();
 State s{x,tangent};Shot result;result.start=s;result.points.push_back(s);double h=z[2]*scale/steps;
 for(int i=0;i<steps;++i){auto a=rhs(s,inv),b=rhs(add(s,a,h/2),inv),c=rhs(add(s,b,h/2),inv),d=rhs(add(s,c,h),inv);s.x+=h/6*(a.x+2*b.x+2*c.x+d.x);s.v+=h/6*(a.v+2*b.v+2*c.v+d.v);result.points.push_back(s);}
 result.end=s;result.residual[0]=~n*incoming;Vec3 mismatch=(s.x+z[3]*scale*s.v-q)/scale;for(int i=0;i<3;++i)result.residual[1+i]=mismatch[i];return result;
}
inline double norm(const Shot&s){double v=0;for(double e:s.residual)v=std::max(v,std::abs(e));return v;}
inline bool domain(const double*z){return std::abs(z[0])<.4&&std::abs(z[1])<.4&&z[2]>1e-4&&z[2]<1.81&&z[3]>1e-4&&z[3]<10.;}
inline void solve4(double a[4][5],double*x){for(int i=0;i<4;++i){int pivot=i;for(int k=i+1;k<4;++k)if(std::abs(a[k][i])>std::abs(a[pivot][i]))pivot=k;if(std::abs(a[pivot][i])<1e-12)throw std::runtime_error("singular geodesic shooting Jacobian");for(int j=i;j<5;++j)std::swap(a[i][j],a[pivot][j]);double d=a[i][i];for(int j=i;j<5;++j)a[i][j]/=d;for(int k=0;k<4;++k)if(k!=i){double v=a[k][i];for(int j=i;j<5;++j)a[k][j]-=v*a[i][j];}}for(int i=0;i<4;++i)x[i]=a[i][4];}
inline void correct(Vec3 p,Vec3 q,Vec3 axes,OpenSim::WrapResult&r){
 Vec3 pinned(.027559229367297551,.02204738349383804,.02204738349383804);if((axes-pinned).norm()>1e-12)throw std::runtime_error("unreviewed ellipsoid dimensions");
 double inside1=0,inside2=0;for(int i=0;i<3;++i){inside1+=p[i]*p[i]/(axes[i]*axes[i]);inside2+=q[i]*q[i]/(axes[i]*axes[i]);}if(std::min(inside1,inside2)<=1+1e-7)throw std::runtime_error("ellipsoid endpoints inside/near surface");
 const auto seed1=r.r1,seed2=r.r2;double scale=axes[1];if(r.wrap_path_length<=1e-6||r.wrap_path_length>=.04)throw std::runtime_error("ellipsoid outside short-branch domain");
 Vec3 y0;for(int i=0;i<3;++i)y0[i]=r.r1[i]/axes[i];y0/=y0.norm();int k=0;for(int i=1;i<3;++i)if(std::abs(y0[i])<std::abs(y0[k]))k=i;Vec3 axis(0);axis[k]=1;Vec3 e1=y0%axis;e1/=e1.norm();Vec3 e2=y0%e1;
 double z[4]={0,0,r.wrap_path_length/scale,(q-r.r2).norm()/scale};Shot current=shoot(z,p,q,axes,y0,e1,e2,128);bool converged=false;
 for(int iteration=0;iteration<18;++iteration){if(norm(current)<1e-11){converged=true;break;}double matrix[4][5];for(int column=0;column<4;++column){double plus[4],minus[4];for(int j=0;j<4;++j)plus[j]=minus[j]=z[j];plus[column]+=1e-5;minus[column]-=1e-5;auto a=shoot(plus,p,q,axes,y0,e1,e2,128),b=shoot(minus,p,q,axes,y0,e1,e2,128);for(int row=0;row<4;++row)matrix[row][column]=(a.residual[row]-b.residual[row])/2e-5;}
  for(int row=0;row<4;++row)matrix[row][4]=-current.residual[row];double step[4];solve4(matrix,step);bool accepted=false;
  for(int backtrack=0;backtrack<12;++backtrack){double trial[4];double alpha=std::ldexp(1.,-backtrack);for(int j=0;j<4;++j)trial[j]=z[j]+alpha*step[j];if(!domain(trial))continue;auto next=shoot(trial,p,q,axes,y0,e1,e2,128);if(norm(next)<norm(current)){for(int j=0;j<4;++j)z[j]=trial[j];current=next;accepted=true;break;}}
  if(!accepted)throw std::runtime_error("geodesic Newton branch did not decrease");
 }
 if(!converged&&norm(current)>=1e-11)throw std::runtime_error("geodesic shooting iteration cap");
 auto refined=shoot(z,p,q,axes,y0,e1,e2,256);if(norm(refined)>1e-8||(refined.end.x-current.end.x).norm()>1e-9)throw std::runtime_error("geodesic integration refinement failed");
 if((current.start.x-seed1).norm()>.002||(current.end.x-seed2).norm()>.002)throw std::runtime_error("geodesic left native contact branch neighborhood");
 for(const auto&s:current.points){double f=0;for(int i=0;i<3;++i)f+=s.x[i]*s.x[i]/(axes[i]*axes[i]);if(s.x[1]>=-1e-8||std::abs(f-1)>1e-8||std::abs(s.v.norm()-1)>1e-8)throw std::runtime_error("geodesic quadrant/surface/speed drift");}
 r.r1=current.start.x;r.r2=current.end.x;r.wrap_path_length=z[2]*scale;r.wrap_pts.setSize(0);for(int i=0;i<=128;i+=4)r.wrap_pts.append(current.points[i].x);
}
}
class StationaryArm26Ellipsoid:public OpenSim::WrapEllipsoid{
 OpenSim_DECLARE_CONCRETE_OBJECT(StationaryArm26Ellipsoid,OpenSim::WrapEllipsoid);
protected:int wrapLine(const SimTK::State&s,SimTK::Vec3&p,SimTK::Vec3&q,const OpenSim::PathWrap&w,OpenSim::WrapResult&r,bool&flag)const override{
 auto p0=p,q0=q;if(get_quadrant()!="-y")throw std::runtime_error("stationary ellipsoid requires retained minus-y quadrant");int result=WrapEllipsoid::wrapLine(s,p,q,w,r,flag);
 if(flag)arm26geodesic::correct(p0,q0,get_dimensions(),r);
 else {auto a=get_dimensions();SimTK::Vec3 x,d;for(int i=0;i<3;++i){x[i]=p0[i]/a[i];d[i]=(q0[i]-p0[i])/a[i];}double t=(~d*d)>0?std::clamp(-(~x*d)/(~d*d),0.,1.):0.;if((x+t*d).norm()<=1+1e-8)throw std::runtime_error("unsupported ellipsoid no-wrap branch");}
 arm26audit::record(*this,p0,q0,r,result,flag,flag);return result;
 }
};
