#include "arm26_wrap_candidate.h"
#include <OpenSim/Actuators/ModelOperators.h>
#include <sstream>
#include <map>
using namespace OpenSim;
int main(int argc,char**argv){try{
 if(argc!=6)throw std::runtime_error("model fitted cases metrics trace required");
 Object::registerType(ObservedArm26Cylinder());Object::registerType(ExactArm26Cylinder());Object::registerType(ObservedArm26Ellipsoid());Logger::setLevelString("warn");
 Model model(argv[1]);model.setUseVisualizer(false);model.set_assembly_accuracy(1e-10);Set<FunctionBasedPath> fitted{std::string(argv[2])};ModelFactory::replacePathsWithFunctionBasedPaths(model,fitted);model.finalizeConnections();auto initial=model.initSystem();
 std::ofstream metrics(argv[4]),trace(argv[5]);metrics<<std::setprecision(17)<<"case,coordinate,muscle,status,q,length,moment_arm,fd,unit_resultant_force,unit_resultant_moment\n";
 std::ifstream cases(argv[3]);std::string line;
 while(std::getline(cases,line)){std::istringstream input(line);std::string label,coordinate,name;double value;input>>label>>coordinate;auto s=initial;while(input>>name>>value)model.updCoordinateSet().get(name).setValue(s,value,false);model.assemble(s);model.realizePosition(s);auto&c=model.updCoordinateSet().get(coordinate);
  for(const auto&m:model.getComponentList<Muscle>()){
   bool target=m.getName().find("arm26_BRA_")==0||m.getName().find("arm26_BIClong_")==0;
   try{
    arm26audit::context=label+":"+m.getName();arm26audit::trace=target?&trace:nullptr;double length=m.getLength(s);arm26audit::trace=nullptr;
    double moment=0,fd=0,resultant=0,torque=0;
    if(target){moment=m.computeMomentArm(s,c);auto plus=s,minus=s;double q=c.getValue(s),h=1e-6;double d1=h,d2=-h;if(q-h<c.getRangeMin()){d1=h;d2=2*h;}else if(q+h>c.getRangeMax()){d1=-h;d2=-2*h;}c.setValue(plus,q+d1,false);c.setValue(minus,q+d2,false);model.realizePosition(plus);model.realizePosition(minus);if(d2==-d1)fd=-(m.getLength(plus)-m.getLength(minus))/(2*h);else fd=-(-3*length+4*m.getLength(plus)-m.getLength(minus))/(2*d1);
     SimTK::Vector_<SimTK::SpatialVec> bf(model.getMatterSubsystem().getNumBodies(),SimTK::SpatialVec(SimTK::Vec3(0),SimTK::Vec3(0)));SimTK::Vector mobility(s.getNU(),0.);m.getGeometryPath().addInEquivalentForces(s,1.,bf,mobility);SimTK::Vec3 totalF(0),totalM(0);
     for(int i=0;i<bf.size();++i){const auto&b=model.getMatterSubsystem().getMobilizedBody(SimTK::MobilizedBodyIndex(i));totalF+=bf[i][1];totalM+=bf[i][0]+b.getBodyOriginLocation(s)%bf[i][1];}resultant=totalF.norm();torque=totalM.norm();
    }
    metrics<<label<<','<<coordinate<<','<<m.getName()<<",ok,"<<c.getValue(s)<<','<<length<<','<<moment<<','<<fd<<','<<resultant<<','<<torque<<'\n';
   }catch(const std::exception&e){arm26audit::trace=nullptr;metrics<<label<<','<<coordinate<<','<<m.getName()<<",rejected,"<<c.getValue(s)<<",0,0,0,0,0\n";trace<<"{\"case\":\""<<label<<":"<<m.getName()<<"\",\"rejection\":\"";for(char v:std::string(e.what())){if(v=='\n'||v=='\r')trace<<' ';else if(v=='"'||v=='\\')trace<<' ';else trace<<v;}trace<<"\"}\n";}
  }
 }
 if(!metrics||!trace)throw std::runtime_error("output failure");return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}}
