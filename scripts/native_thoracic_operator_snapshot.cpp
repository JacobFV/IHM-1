// Read-only same-state Simbody tree operator extraction; no integration/forces.
#include <OpenSim/OpenSim.h>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <filesystem>

static void number(std::ostream& out,double x){if(!std::isfinite(x))throw std::runtime_error("Nonfinite operator snapshot");out<<std::setprecision(17)<<x;}
template<class V>static void vector(std::ostream& out,const V& v){out<<'[';for(int i=0;i<v.size();++i){if(i)out<<',';number(out,v[i]);}out<<']';}
template<class M>static void matrix(std::ostream& out,const M& m,int rows,int cols){out<<'[';for(int i=0;i<rows;++i){if(i)out<<',';out<<'[';for(int j=0;j<cols;++j){if(j)out<<',';number(out,m(i,j));}out<<']';}out<<']';}
int main(int argc,char** argv){try{
    if(argc!=3)throw std::runtime_error("model.osim fresh-output.json required");
    if(std::filesystem::exists(argv[2]))throw std::runtime_error("Refuse existing evidence output");
    OpenSim::Model model(argv[1]);model.setUseVisualizer(false);auto& state=model.initSystem();
    for(int i=0;i<state.getNU();++i)state.updU()[i]=.001*(i+1);
    const auto& system=model.getMultibodySystem();system.realize(state,SimTK::Stage::Velocity);
    const auto& matter=system.getMatterSubsystem();const auto& torso=model.getBodySet().get("torso");const auto index=torso.getMobilizedBodyIndex();
    SimTK::Matrix mass,jacobian;matter.calcM(state,mass);matter.calcFrameJacobian(state,index,SimTK::Vec3(0),jacobian);
    SimTK::Vector bias;matter.calcResidualForceIgnoringConstraints(state,SimTK::Vector(),SimTK::Vector_<SimTK::SpatialVec>(),SimTK::Vector(),bias);
    const auto frame_bias=matter.calcBiasForFrameJacobian(state,index,SimTK::Vec3(0));const auto transform=matter.getMobilizedBody(index).getBodyTransform(state);
    SimTK::Vector qdot;matter.multiplyByN(state,false,state.getU(),qdot);
    std::ofstream out(argv[2]);if(!out)throw std::runtime_error("Cannot open output");
    out<<"{\"schema\":\"thoracic-native-tree-operator-v1\",\"native_integrated\":false,\"native_constraints_projected\":false,\"zero_applied_force_inverse_dynamics\":true,\"q\":";vector(out,state.getQ());
    out<<",\"u\":";vector(out,state.getU());out<<",\"qdot\":";vector(out,qdot);
    out<<",\"mass_matrix\":";matrix(out,mass,mass.nrow(),mass.ncol());out<<",\"inertial_bias\":";vector(out,bias);
    out<<",\"torso_frame_jacobian_world_angular_first\":";matrix(out,jacobian,6,jacobian.ncol());
    out<<",\"torso_frame_bias_world_angular_first\":[";for(int part=0;part<2;++part)for(int k=0;k<3;++k){if(part||k)out<<',';number(out,frame_bias[part][k]);}out<<']';
    out<<",\"torso_rotation_body_to_world\":";matrix(out,transform.R().asMat33(),3,3);
    out<<",\"torso_translation_world_m\":";vector(out,transform.p());out<<",\"torso_mass_kg\":";number(out,torso.getMass());
    out<<",\"declared_native_constraint_components\":"<<model.getConstraintSet().getSize()<<",\"body_count\":"<<model.getBodySet().getSize()<<",\"source_muscle_count\":"<<model.getMuscles().getSize()<<"}\n";
    if(!out)throw std::runtime_error("Snapshot write failed");return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
