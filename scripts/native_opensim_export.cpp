// Thin official OpenSim API exporter. No mechanics equations or wrap algorithms replaced.
#include <OpenSim/OpenSim.h>
#include <OpenSim/Simulation/Wrap/PathWrapPoint.h>
#include <OpenSim/Common/About.h>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <filesystem>
using namespace OpenSim;
static void number(std::ostream& out,double x) { if(!std::isfinite(x)) throw std::runtime_error("Nonfinite engine mechanics"); out<<std::setprecision(17)<<x; }
static void str(std::ostream& out,const std::string& s) {out<<'"';for(char c:s){if(c=='"'||c=='\\')out<<'\\';if(c=='\n')out<<"\\n";else out<<c;}out<<'"';}
static void point(std::ostream& out,const SimTK::Vec3& v){out<<'[';for(int j=0;j<3;++j){if(j)out<<',';number(out,v[j]);}out<<']';}
int main(int argc,char** argv) {
  try {
    if(argc!=6)throw std::runtime_error("model output coordinate delta_rad activation required");
    const double delta=std::stod(argv[4]),activation=std::stod(argv[5]);
    if(!std::isfinite(delta)||!std::isfinite(activation)||std::abs(delta)>.15||activation<.01||activation>.5)throw std::runtime_error("Out of supported range");
    auto modelpath=std::filesystem::path(argv[1]);
    ModelVisualizer::addDirToGeometrySearchPaths((modelpath.parent_path()/"Geometry").string());
    ModelVisualizer::addDirToGeometrySearchPaths((modelpath.parent_path().parent_path().parent_path()/"Geometry").string());
    Model model(argv[1]);model.setUseVisualizer(false);model.set_assembly_accuracy(1e-12);
    auto& state=model.initSystem();
    const auto& target=model.getCoordinateSet().get(argv[3]);
    if(target.getMotionType()!=Coordinate::Rotational || target.isDependent(state) || target.getLocked(state))throw std::runtime_error("Expected unlocked independent rotational coordinate");
    const double requested=target.getDefaultValue()+delta;
    if(requested<target.getRangeMin() || requested>target.getRangeMax())throw std::runtime_error("Requested coordinate outside source range");
    target.setValue(state,requested,true);
    if(std::abs(target.getValue(state)-requested)>1e-8)throw std::runtime_error("Assembly failed requested coordinate tolerance");
    for(auto& muscle:model.updComponentList<Muscle>())muscle.setActivation(state,activation);
    model.equilibrateMuscles(state);model.realizeDynamics(state);
    std::ofstream out(argv[2]);if(!out)throw std::runtime_error("Cannot open output");
    out<<"{\"schema_version\":1,\"engine\":\"OpenSim\",\"version\":";str(out,GetVersion());
    out<<",\"activation\":";number(out,activation);out<<",\"perturbation_coordinate\":";str(out,target.getName());out<<",\"delta_rad\":";number(out,delta);
    out<<",\"coordinates\":{";bool first=true;
    std::vector<const Coordinate*> rotationCoords;
    for(int i=0;i<model.getCoordinateSet().getSize();++i){const auto& c=model.getCoordinateSet().get(i);if(!first)out<<',';first=false;str(out,c.getName());out<<':';number(out,c.getValue(state));if(c.getMotionType()==Coordinate::Rotational&&!c.isDependent(state)&&!c.getLocked(state))rotationCoords.push_back(&c);}
    out<<"},\"body_transforms_ground\":{";first=true;
    for(const auto& body:model.getComponentList<Body>()){
      if(!first)out<<',';first=false;str(out,body.getName());out<<": [";
      auto x=body.getTransformInGround(state);
      for(int i=0;i<4;++i){if(i)out<<',';out<<'[';for(int j=0;j<4;++j){if(j)out<<',';number(out,i==3?(j==3?1:0):(j==3?x.p()[i]:x.R().asMat33()(i,j)));}out<<']';}out<<']';
    }
    out<<"},\"muscles\":[";first=true;
    for(const auto& muscle:model.getComponentList<Muscle>()){
      if(!first)out<<',';first=false;out<<"{\"name\":";str(out,muscle.getName());
      const auto& path=muscle.getGeometryPath();
      out<<",\"length_m\":";number(out,path.getLength(state));
      out<<",\"tendon_force_N\":";number(out,muscle.getTendonForce(state));
      out<<",\"active_fiber_force_N\":";number(out,muscle.getActiveFiberForce(state));
      out<<",\"passive_fiber_force_N\":";number(out,muscle.getPassiveFiberForce(state));
      out<<",\"fiber_length_m\":";number(out,muscle.getFiberLength(state));
      out<<",\"fiber_velocity_m_s\":";number(out,muscle.getFiberVelocity(state));
      out<<",\"fiber_force_along_tendon_N\":";number(out,muscle.getActiveFiberForceAlongTendon(state)+muscle.getPassiveFiberForceAlongTendon(state));
      out<<",\"moment_arms_m\":{";bool fc=true;
      for(const auto* c:rotationCoords){if(!fc)out<<',';fc=false;str(out,c->getName());out<<':';number(out,path.computeMomentArm(state,*c));}
      out<<"},\"path_ground_m\":[";
      const auto& current=path.getCurrentPath(state);int wrapped=0;
      point(out,current[0]->getLocationInGround(state));
      // Identical traversal to upstream GeometryPath::generateDecorations.
      for(int i=1;i<current.getSize();++i){auto* wp=dynamic_cast<PathWrapPoint*>(current[i]);
        if(wp){const auto& points=wp->getWrapPath(state);const auto& x=wp->getParentFrame().getTransformInGround(state);if(points.getSize())++wrapped;for(int j=0;j<points.getSize();++j){out<<',';point(out,x*points[j]);}}
        else {out<<',';point(out,current[i]->getLocationInGround(state));}
      }
      out<<"],\"active_wrap_segments\":"<<wrapped<<",\"wrap_solved\":true";
      out<<",\"recomputed_lengths_m\":[";
      SimTK::State diagnostic=state;
      for(int k=0;k<3;++k){diagnostic.invalidateAllCacheAtOrAbove(SimTK::Stage::Position);model.realizePosition(diagnostic);if(k)out<<',';number(out,path.getLength(diagnostic));}
      out<<"],\"recomputed_moment_arm_m\":";number(out,path.computeMomentArm(diagnostic,target));
      out<<",\"path_nodes\":[";
      const auto& nodes=path.getCurrentPath(state);
      for(int k=0;k<nodes.getSize();++k){if(k)out<<',';out<<"{\"point_ground_m\":";point(out,nodes[k]->getLocationInGround(state));out<<",\"body\":";str(out,nodes[k]->getParentFrame().getName());out<<",\"wrap_object\":";str(out,nodes[k]->getWrapObject()?nodes[k]->getWrapObject()->getName():"");out<<'}';}
      out<<"],\"max_wrap_location_cache_error_m\":";
      double cache_error=0;
      for(int k=0;k<nodes.getSize();++k)if(nodes[k]->getWrapObject())cache_error=std::max(cache_error,(nodes[k]->getLocationInGround(state)-nodes[k]->getParentFrame().getTransformInGround(state)*nodes[k]->getLocation(state)).norm());
      number(out,cache_error);out<<'}';
    }
    out<<"],\"equilibrated\":true,\"force_scope\":\"source muscle equilibrium at explicit activation; external skeletal force balance not solved\"}\n";
    std::cout<<"EXPORTED="<<argv[2]<<"\n";return 0;
  }catch(const std::exception& e){std::cerr<<"EXPORT_ERROR="<<e.what()<<"\n";return 2;}
}
