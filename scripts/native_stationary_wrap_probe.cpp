#include "arm26_stationary_ellipsoid.h"
#include <OpenSim/Actuators/ModelOperators.h>
#include <sstream>
#include <map>
using namespace OpenSim;
int main(int argc,char**argv){try{
 if(argc!=6)throw std::runtime_error("model fitted cases metrics trace required");
 Object::registerType(StationaryArm26Ellipsoid());Object::registerType(ObservedArm26Cylinder());Object::registerType(ExactArm26Cylinder());Object::registerType(ObservedArm26Ellipsoid());Logger::setLevelString("warn");
 {
  auto fixture=[](double z){OpenSim::WrapResult r;double radius=.012,rho=.04;double x=radius*radius/rho,y=std::sqrt(radius*radius-x*x);r.r1=SimTK::Vec3(x,y,z);r.r2=SimTK::Vec3(-x,y,z);r.wrap_path_length=radius*std::atan2(2*x*y,y*y-x*x);return r;};
  int rejected=0;for(int k=0;k<4;++k){auto r=fixture(k==1?.04:0);SimTK::Vec3 p(.04,0,k==1?.04:0),q(-.04,0,k==1?.04:0);bool flag=true;if(k==0)p=SimTK::Vec3(0);if(k==2)flag=false;if(k==3)r.wrap_path_length+=2*SimTK_PI*.012;
   try{arm26audit::correctCylinder(p,q,.012,.04,1,flag,r);}catch(const std::exception&){++rejected;}}
  if(rejected!=4)throw std::runtime_error("cylinder rejection guard failed");std::cout<<"cylinder_boundary_rejections=4\n";
 }
 Model model(argv[1]);model.setUseVisualizer(false);model.set_assembly_accuracy(1e-10);Set<FunctionBasedPath> fitted{std::string(argv[2])};ModelFactory::replacePathsWithFunctionBasedPaths(model,fitted);model.finalizeConnections();auto initial=model.initSystem();
 std::ofstream metrics(argv[4]),trace(argv[5]);metrics<<std::setprecision(17)<<"case,coordinate,muscle,status,q,length,moment_arm,fd,unit_resultant_force,unit_resultant_moment,lengthening_speed,unit_body_power\n";
 std::ifstream cases(argv[3]);std::string line;
 while(std::getline(cases,line)){std::istringstream input(line);std::string label,coordinate,name;double value;input>>label>>coordinate;auto s=initial;while(input>>name>>value)model.updCoordinateSet().get(name).setValue(s,value,false);model.assemble(s);model.realizePosition(s);auto&c=model.updCoordinateSet().get(coordinate);
  for(const auto&m:model.getComponentList<Muscle>()){
   bool target=m.getName().find("arm26_BRA_")==0||m.getName().find("arm26_BIClong_")==0;
   try{
    arm26audit::context=label+":"+m.getName();arm26audit::trace=target?&trace:nullptr;double length=m.getLength(s);arm26audit::trace=nullptr;
    double moment=0,fd=0,resultant=0,torque=0,speed=0,power=0;
    if(target){moment=m.computeMomentArm(s,c);auto plus=s,minus=s;double q=c.getValue(s),h=1e-6;double d1=h,d2=-h;if(q-h<c.getRangeMin()){d1=h;d2=2*h;}else if(q+h>c.getRangeMax()){d1=-h;d2=-2*h;}c.setValue(plus,q+d1,false);c.setValue(minus,q+d2,false);model.realizePosition(plus);model.realizePosition(minus);if(d2==-d1)fd=-(m.getLength(plus)-m.getLength(minus))/(2*h);else fd=-(-3*length+4*m.getLength(plus)-m.getLength(minus))/(2*d1);
     SimTK::Vector_<SimTK::SpatialVec> bf(model.getMatterSubsystem().getNumBodies(),SimTK::SpatialVec(SimTK::Vec3(0),SimTK::Vec3(0)));SimTK::Vector mobility(s.getNU(),0.);m.getGeometryPath().addInEquivalentForces(s,1.,bf,mobility);SimTK::Vec3 totalF(0),totalM(0);
     for(int i=0;i<bf.size();++i){const auto&b=model.getMatterSubsystem().getMobilizedBody(SimTK::MobilizedBodyIndex(i));totalF+=bf[i][1];totalM+=bf[i][0]+b.getBodyOriginLocation(s)%bf[i][1];}resultant=totalF.norm();torque=totalM.norm();auto moving=s;moving.updU()=0;c.setSpeedValue(moving,.07);model.realizeVelocity(moving);speed=m.getLengtheningSpeed(moving);
     for(int i=0;i<bf.size();++i){const auto&b=model.getMatterSubsystem().getMobilizedBody(SimTK::MobilizedBodyIndex(i));const auto&v=b.getBodyVelocity(moving);power+=(~v[0]*bf[i][0])+(~v[1]*bf[i][1]);}power+=~mobility*moving.getU();
    }
    metrics<<label<<','<<coordinate<<','<<m.getName()<<",ok,"<<c.getValue(s)<<','<<length<<','<<moment<<','<<fd<<','<<resultant<<','<<torque<<','<<speed<<','<<power<<'\n';
   }catch(const std::exception&e){arm26audit::trace=nullptr;metrics<<label<<','<<coordinate<<','<<m.getName()<<",rejected,"<<c.getValue(s)<<",0,0,0,0,0,0,0\n";trace<<"{\"case\":\""<<label<<":"<<m.getName()<<"\",\"rejection\":\"";for(char v:std::string(e.what())){if(v=='\n'||v=='\r')trace<<' ';else if(v=='"'||v=='\\')trace<<' ';else trace<<v;}trace<<"\"}\n";}
  }
 }
 if(!metrics||!trace)throw std::runtime_error("output failure");return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}}
