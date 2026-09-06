// Isolated experiment. No replacement of registered OpenSim base classes.
#pragma once
#include <OpenSim/OpenSim.h>
#include <OpenSim/Simulation/Wrap/WrapResult.h>
#include <fstream>
#include <iomanip>
#include <cmath>
namespace arm26audit {
inline std::ostream* trace=nullptr;inline std::string context;
inline void vec(std::ostream&o,const SimTK::Vec3&v){o<<'['<<v[0]<<','<<v[1]<<','<<v[2]<<']';}
inline void record(const OpenSim::WrapObject&w,const SimTK::Vec3&p,const SimTK::Vec3&q,const OpenSim::WrapResult&r,int status,bool flag,bool corrected){
 if(!trace)return;auto&o=*trace;o<<std::setprecision(17)<<"{\"case\":\""<<context<<"\",\"object\":\""<<w.getName()<<"\",\"type\":\""<<w.getConcreteClassName()<<"\",\"status\":"<<status<<",\"wrapped\":"<<(flag?"true":"false")<<",\"corrected\":"<<(corrected?"true":"false")<<",\"single_wrap\":"<<(r.singleWrap?"true":"false")<<",\"p1\":";vec(o,p);o<<",\"p2\":";vec(o,q);
 o<<",\"quadrant\":\""<<w.get_quadrant()<<"\",\"frame\":\""<<w.getFrame().getName()<<"\",\"transform\":[";
 const auto&x=w.getTransform();for(int i=0;i<3;++i){if(i)o<<',';o<<'[';for(int j=0;j<3;++j)o<<x.R().asMat33()(i,j)<<',';o<<x.p()[i]<<']';}o<<']';
 if(flag){o<<",\"r1\":";vec(o,r.r1);o<<",\"r2\":";vec(o,r.r2);o<<",\"wrap_length\":"<<r.wrap_path_length<<",\"wrap_points\":[";for(int i=0;i<r.wrap_pts.getSize();++i){if(i)o<<',';vec(o,r.wrap_pts[i]);}o<<']';
 if(std::isfinite(r.factor)&&r.factor>0){o<<",\"factor\":"<<r.factor<<",\"c1_unscaled\":";vec(o,r.c1/r.factor);o<<",\"sv_unscaled\":";vec(o,r.sv/r.factor);}}
 o<<"}\n";
}
inline void correctCylinder(const SimTK::Vec3&p,const SimTK::Vec3&q,double radius,double length,int status,bool flag,OpenSim::WrapResult&r){
 const double rho1=std::hypot(p[0],p[1]),rho2=std::hypot(q[0],q[1]);
 if(rho1<=radius+1e-9||rho2<=radius+1e-9)throw std::runtime_error("cylinder endpoint inside/near surface");
 if(!flag){SimTK::Vec3 d=q-p;double den=d[0]*d[0]+d[1]*d[1];double t=den?std::clamp(-(p[0]*d[0]+p[1]*d[1])/den,0.,1.):0.;if(std::hypot(p[0]+t*d[0],p[1]+t*d[1])<=radius+1e-9)throw std::runtime_error("unsupported no-wrap cap/branch case");return;}
 auto a=r.r1,b=r.r2;double ra=std::hypot(a[0],a[1]),rb=std::hypot(b[0],b[1]);
 if(std::abs(ra-radius)>1e-9||std::abs(rb-radius)>1e-9)throw std::runtime_error("native contacts not cylindrical");
 double d1=std::hypot(p[0]-a[0],p[1]-a[1]),d2=std::hypot(q[0]-b[0],q[1]-b[1]);
 if(std::abs((p[0]-a[0])*a[0]+(p[1]-a[1])*a[1])>1e-10||std::abs((q[0]-b[0])*b[0]+(q[1]-b[1])*b[1])>1e-10)throw std::runtime_error("native XY contacts not tangent");
 double arc2=r.wrap_path_length*r.wrap_path_length-(b[2]-a[2])*(b[2]-a[2]);if(arc2<=0)throw std::runtime_error("degenerate native arc");
 double theta=std::sqrt(arc2)/radius;
 double incoming=(a[0]-p[0])*(-a[1])+(a[1]-p[1])*a[0];double outgoing=(q[0]-b[0])*(-b[1])+(q[1]-b[1])*b[0];
 int sense=incoming>0?1:-1;if(incoming*outgoing<=0||theta<=1e-8||theta>=2*SimTK_PI-1e-8)throw std::runtime_error("unsupported winding branch");
 double c=std::cos(sense*theta),sn=std::sin(sense*theta);
 if(std::hypot(c*a[0]-sn*a[1]-b[0],sn*a[0]+c*a[1]-b[1])>1e-9)throw std::runtime_error("native arc winding does not join contacts");
 double horizontal=d1+radius*theta+d2,dz=q[2]-p[2];a[2]=p[2]+dz*d1/horizontal;b[2]=p[2]+dz*(d1+radius*theta)/horizontal;
 if(std::abs(a[2])>=length/2-1e-8||std::abs(b[2])>=length/2-1e-8)throw std::runtime_error("unsupported finite cylinder cap contact");
 r.r1=a;r.r2=b;r.wrap_path_length=std::hypot(radius*theta,b[2]-a[2]);r.wrap_pts.setSize(0);
 for(int i=0;i<=16;++i){double t=i/16.,cs=std::cos(sense*theta*t),ss=std::sin(sense*theta*t);r.wrap_pts.append(SimTK::Vec3(cs*a[0]-ss*a[1],ss*a[0]+cs*a[1],a[2]+t*(b[2]-a[2])));}
}
}
class ObservedArm26Cylinder:public OpenSim::WrapCylinder{
 OpenSim_DECLARE_CONCRETE_OBJECT(ObservedArm26Cylinder,OpenSim::WrapCylinder);
protected:int wrapLine(const SimTK::State&s,SimTK::Vec3&p,SimTK::Vec3&q,const OpenSim::PathWrap&w,OpenSim::WrapResult&r,bool&flag)const override{auto p0=p,q0=q;int result=WrapCylinder::wrapLine(s,p,q,w,r,flag);arm26audit::record(*this,p0,q0,r,result,flag,false);return result;}
};
class ExactArm26Cylinder:public OpenSim::WrapCylinder{
 OpenSim_DECLARE_CONCRETE_OBJECT(ExactArm26Cylinder,OpenSim::WrapCylinder);
protected:int wrapLine(const SimTK::State&s,SimTK::Vec3&p,SimTK::Vec3&q,const OpenSim::PathWrap&w,OpenSim::WrapResult&r,bool&flag)const override{auto p0=p,q0=q;if(get_quadrant()!="all")throw std::runtime_error("candidate requires retained all quadrant");int result=WrapCylinder::wrapLine(s,p,q,w,r,flag);arm26audit::correctCylinder(p0,q0,get_radius(),get_length(),result,flag,r);arm26audit::record(*this,p0,q0,r,result,flag,true);return result;}
};
class ObservedArm26Ellipsoid:public OpenSim::WrapEllipsoid{
 OpenSim_DECLARE_CONCRETE_OBJECT(ObservedArm26Ellipsoid,OpenSim::WrapEllipsoid);
protected:int wrapLine(const SimTK::State&s,SimTK::Vec3&p,SimTK::Vec3&q,const OpenSim::PathWrap&w,OpenSim::WrapResult&r,bool&flag)const override{auto p0=p,q0=q;int result=WrapEllipsoid::wrapLine(s,p,q,w,r,flag);arm26audit::record(*this,p0,q0,r,result,flag,false);return result;}
};
