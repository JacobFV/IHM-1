// Standalone geometry/constitutive acceptance; never advances a patient.
#include <OpenSim/OpenSim.h>
#include <OpenSim/Actuators/ModelOperators.h>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace OpenSim;
int main(int argc,char**argv){try{
 if(argc!=4)throw std::runtime_error("model fitted_paths output required");
 Logger::setLevelString("warn");Model model(argv[1]);model.setUseVisualizer(false);
 Set<FunctionBasedPath> fitted{std::string(argv[2])};
 ModelFactory::replacePathsWithFunctionBasedPaths(model,fitted);
 model.finalizeConnections();auto initial=model.initSystem();
 int nm=0;for(const auto&m:model.getComponentList<Muscle>())++nm;
 if(nm!=98||model.getBodySet().getSize()!=22)throw std::runtime_error("coverage changed");
 std::ofstream out(argv[3]);out<<std::setprecision(17);
 out<<"pose,coordinate,q,muscle,axis,length,moment_arm,finite_difference,tendon_force,active_fiber_force,passive_fiber_force,cos_pennation,activation,fiber_velocity\n";
 const std::string coords[]={"lumbar_extension","lumbar_bending","lumbar_rotation"};
 for(int pose=0;pose<7;++pose){auto s=initial;int varied=(pose-1)/2;double q=pose==0?0:((pose%2)?-.1:.1);
  for(int k=0;k<3;++k)model.updCoordinateSet().get(coords[k]).setValue(s,pose&&k==varied?q:0,false);
  model.realizePosition(s);
  // Equilibrate at the unmodified native default activation, never select a
  // balancing activation. No time integration, external load or reserve solve.
  model.equilibrateMuscles(s);model.realizeDynamics(s);
  for(const auto&m:model.getComponentList<Muscle>())if(m.getName().find("gait2392_")==0){
   double length=m.getLength(s),force=m.getTendonForce(s);
   if(!std::isfinite(force)||force<0)throw std::runtime_error("invalid tensile force");
   for(int k=0;k<3;++k){auto&c=model.updCoordinateSet().get(coords[k]);double value=c.getValue(s),h=1e-6;
    auto plus=s,minus=s;c.setValue(plus,value+h,false);c.setValue(minus,value-h,false);model.realizePosition(plus);model.realizePosition(minus);
    double fd=-(m.getLength(plus)-m.getLength(minus))/(2*h),arm=m.computeMomentArm(s,c);
    if(std::abs(arm-fd)>1e-8)throw std::runtime_error("virtual work mismatch");
    out<<pose<<","<<(pose?coords[varied]:"neutral")<<","<<q<<","<<m.getName()<<","<<coords[k]<<","<<length<<","<<arm<<","<<fd<<","<<force<<","<<m.getActiveFiberForce(s)<<","<<m.getPassiveFiberForce(s)<<","<<std::cos(m.getPennationAngle(s))<<","<<m.getActivation(s)<<","<<m.getFiberVelocity(s)<<"\n";
   }
  }
 }
 if(!out)throw std::runtime_error("output failed");return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}}
