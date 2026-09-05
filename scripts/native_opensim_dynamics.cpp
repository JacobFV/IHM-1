// Native forward dynamics of the retained 3D walking subject. Coordinates are
// initialized once from measurements; subsequent states are freely integrated.
#include <OpenSim/OpenSim.h>
#include <OpenSim/Actuators/ModelOperators.h>
#include <OpenSim/Common/STOFileAdapter.h>
#include <fstream>
#include <iomanip>
#include <filesystem>
#include <cmath>
#include <algorithm>
using namespace OpenSim;
namespace fs=std::filesystem;
static void num(std::ostream& o,double x){if(!std::isfinite(x))throw std::runtime_error("Nonfinite native dynamics");o<<std::setprecision(17)<<x;}
static void str(std::ostream& o,const std::string& s){o<<'"';for(char c:s){if(c=='"'||c=='\\')o<<'\\';o<<c;}o<<'"';}
template<class V>static void vec(std::ostream& o,const V& v){o<<'[';for(int i=0;i<v.size();++i){if(i)o<<',';num(o,v[i]);}o<<']';}
// Bounded, instantaneous engineering velocity feedback. No adaptive-integrator
// callback mutates controller history. This is explicitly not a neural gait policy.
class BodyMuscleController final:public Controller{
    OpenSim_DECLARE_CONCRETE_OBJECT(BodyMuscleController,Controller);
public:
    OpenSim_DECLARE_PROPERTY(baseline,double,"Constant excitation, dimensionless");
    OpenSim_DECLARE_PROPERTY(velocity_gain,double,"Instantaneous positive velocity gain, dimensionless");
    BodyMuscleController(){constructProperty_baseline(.03);constructProperty_velocity_gain(0);}
    void computeControls(const SimTK::State& s,SimTK::Vector& controls)const override{
        const auto& socket=getSocket<Actuator>("actuators");
        for(int i=0;i<(int)socket.getNumConnectees();++i){
            const auto& m=dynamic_cast<const Muscle&>(socket.getConnectee(i));
            const double velocity=m.getLengtheningSpeed(s)/(m.getOptimalFiberLength()*m.getMaxContractionVelocity());
            const double u=std::clamp(get_baseline()+get_velocity_gain()*std::max(0.,velocity),.01,1.);
            m.addInControls(SimTK::Vector(1,u),controls);
        }
    }
};
int main(int argc,char** argv){try{
    if(argc!=12)throw std::runtime_error("source output seconds hz excitation gain dv stiffness friction accuracy max_step required");
    fs::path source=argv[1],out=argv[2];
    double duration=std::stod(argv[3]),hz=std::stod(argv[4]),excitation=std::stod(argv[5]),gain=std::stod(argv[6]),dv=std::stod(argv[7]),stiff=std::stod(argv[8]),friction=std::stod(argv[9]),accuracy=std::stod(argv[10]),maxstep=std::stod(argv[11]);
    // argc check is also enforced by the Python typed interface.
    Model model((source/"subject_walk_scaled.osim").string());model.setUseVisualizer(false);
    model.set_assembly_accuracy(1e-10);
    ForceSet passive((source/"subject_walk_scaled_ExpressionBasedCoordinateForceSet.xml").string());
    for(int i=0;i<passive.getSize();++i)model.addComponent(passive.get(i).clone());
    ContactGeometrySet geometry((source/"subject_walk_scaled_ContactGeometrySet.xml").string());
    for(int i=0;i<geometry.getSize();++i){auto* g=geometry.get(i).clone();if(g->getName()!="floor")g->upd_location()[1]+=.02;model.addContactGeometry(g);}
    ForceSet contacts((source/"subject_walk_scaled_ContactForceSet.xml").string());
    for(int i=0;i<contacts.getSize();++i){auto* f=dynamic_cast<SmoothSphereHalfSpaceForce*>(contacts.get(i).clone());
        f->set_stiffness(f->get_stiffness()*stiff);f->set_static_friction(f->get_static_friction()*friction);f->set_dynamic_friction(f->get_dynamic_friction()*friction);model.addComponent(f);}
    // Retain source Millard muscle and compliant tendon laws. Replace only paths
    // with the source fitted functions, avoiding known geometric-wrap instability.
    ModelFactory::replacePathsWithFunctionBasedPaths(model,(source/"subject_walk_scaled_FunctionBasedPathSet.xml").string());
    auto* control=new BodyMuscleController;control->setName("engineering_velocity_feedback");control->set_baseline(excitation);control->set_velocity_gain(gain);
    for(const auto& m:model.getComponentList<Muscle>())control->addActuator(m);
    model.addController(control);model.finalizeConnections();
    auto& initial=model.initSystem();
    TimeSeriesTable measured((source/"coordinates.sto").string());
    const auto& times=measured.getIndependentColumn();int row=(int)measured.getNearestRowIndexForTime(.48);
    const auto& labels=measured.getColumnLabels();
    for(const auto& c:model.getComponentList<Coordinate>()){
        if(c.isDependent(initial)||c.getLocked(initial))continue;
        auto label=c.getAbsolutePathString()+"/value";auto it=std::find(labels.begin(),labels.end(),label);
        if(it!=labels.end()){int column=(int)(it-labels.begin());c.setValue(initial,measured.getRowAtIndex(row)[column],false);c.setSpeedValue(initial,(measured.getRowAtIndex(row+1)[column]-measured.getRowAtIndex(row-1)[column])/(times[row+1]-times[row-1]));}
    }
    auto& pelvis=model.updCoordinateSet().get("pelvis_tx");pelvis.setSpeedValue(initial,pelvis.getSpeedValue(initial)+dv);
    model.assemble(initial);model.getMultibodySystem().projectU(initial,1e-10);
    for(const auto& m:model.getComponentList<Muscle>())m.setActivation(initial,excitation);
    model.equilibrateMuscles(initial);initial.setTime(0);model.realizeDynamics(initial);
    model.print((out/"assembled_model.osim").string());
    std::ofstream metadata(out/"model.json");
    metadata<<"{\"schema_version\":1,\"mass_kg\":";num(metadata,model.getTotalMass(initial));
    int muscles=0,bodies=0,ncontacts=0,prescribed=0;for(const auto& m:model.getComponentList<Muscle>())++muscles;for(const auto& b:model.getComponentList<Body>())++bodies;for(const auto& c:model.getComponentList<SmoothSphereHalfSpaceForce>())++ncontacts;
    metadata<<",\"coordinates\":[";bool first=true;
    for(const auto& c:model.getComponentList<Coordinate>()){if(!first)metadata<<',';first=false;prescribed+=c.isPrescribed(initial);metadata<<"{\"name\":";str(metadata,c.getName());metadata<<",\"unit\":";str(metadata,c.getMotionType()==Coordinate::Rotational?"rad":"m");metadata<<",\"dependent\":"<<(c.isDependent(initial)?"true":"false")<<",\"locked\":"<<(c.getLocked(initial)?"true":"false")<<'}';}
    metadata<<"],\"state_labels\":[";auto names=model.getStateVariableNames();for(int i=0;i<names.size();++i){if(i)metadata<<',';str(metadata,names[i]);}
    metadata<<"],\"muscle_count\":"<<muscles<<",\"body_count\":"<<bodies<<",\"contact_count\":"<<ncontacts<<",\"prescribed_coordinate_count\":"<<prescribed<<",\"reference_initial_time_s\":0.48,\"scope\":\"freely integrated source subject; fitted muscle paths; uncalibrated controller; no canonical registration or tissue coupling\"}\n";
    std::ofstream frames(out/"frames.jsonl");
    auto emit=[&](const SimTK::State& s){
        model.realizeAcceleration(s);frames<<"{\"time_s\":";num(frames,s.getTime());frames<<",\"q\":";vec(frames,s.getQ());frames<<",\"u\":";vec(frames,s.getU());frames<<",\"state\":";vec(frames,model.getStateVariableValues(s));
        frames<<",\"constraint_position_error\":";num(frames,s.getQErr().norm());frames<<",\"constraint_velocity_error\":";num(frames,s.getUErr().norm());
        frames<<",\"com_m\":";vec(frames,model.calcMassCenterPosition(s));frames<<",\"com_velocity_m_s\":";vec(frames,model.calcMassCenterVelocity(s));frames<<",\"kinetic_energy_J\":";num(frames,model.calcKineticEnergy(s));frames<<",\"potential_energy_J\":";num(frames,model.calcPotentialEnergy(s));
        SimTK::Vec3 total(0);frames<<",\"contacts\":[";bool first=true;
        for(const auto& f:model.getComponentList<SmoothSphereHalfSpaceForce>()){if(!first)frames<<',';first=false;auto wrench=f.getSphereForce(s);total+=wrench[1];frames<<"{\"name\":";str(frames,f.getName());frames<<",\"force_N\":";vec(frames,wrench[1]);frames<<",\"body_origin_moment_Nm\":";vec(frames,wrench[0]);frames<<",\"force_pair_residual_N\":";vec(frames,wrench[1]+f.getHalfSpaceForce(s)[1]);
            const auto& sphere=f.getConnectee<ContactSphere>("sphere");auto center=sphere.getFrame().findStationLocationInGround(s,sphere.get_location());
            frames<<",\"sphere_center_m\":";vec(frames,center);frames<<",\"penetration_m\":";num(frames,std::max(0.,sphere.getRadius()-center[1]));frames<<'}';}
        frames<<"],\"contact_force_N\":";vec(frames,total);
        frames<<",\"momentum_balance_residual_N\":";vec(frames,model.getTotalMass(s)*(model.calcMassCenterAcceleration(s)-model.getGravity())-total);
        frames<<",\"muscles\":[";first=true;
        for(const auto& m:model.getComponentList<Muscle>()){if(!first)frames<<',';first=false;frames<<"{\"name\":";str(frames,m.getName());frames<<",\"activation\":";num(frames,m.getActivation(s));frames<<",\"excitation\":";num(frames,m.getExcitation(s));frames<<",\"length_m\":";num(frames,m.getLength(s));frames<<",\"fiber_length_m\":";num(frames,m.getFiberLength(s));frames<<",\"fiber_velocity_m_s\":";num(frames,m.getFiberVelocity(s));frames<<",\"tendon_force_N\":";num(frames,m.getTendonForce(s));frames<<",\"mechanical_power_W\":";num(frames,-m.getTendonForce(s)*m.getLengtheningSpeed(s));frames<<'}';}
        frames<<"],\"body_transforms_ground\":{";first=true;
        for(const auto& b:model.getComponentList<Body>()){if(!first)frames<<',';first=false;str(frames,b.getName());auto x=b.getTransformInGround(s);frames<<":[";for(int i=0;i<4;++i){if(i)frames<<',';frames<<'[';for(int j=0;j<4;++j){if(j)frames<<',';num(frames,i==3?(j==3?1.:0.):(j==3?x.p()[i]:x.R().asMat33()(i,j)));}frames<<']';}frames<<']';}frames<<"}}\n";frames.flush();
    };
    emit(initial);Manager manager(model);manager.setIntegratorAccuracy(accuracy);manager.setIntegratorConstraintTolerance(1e-8);manager.setIntegratorMaximumStepSize(maxstep);manager.setIntegratorInternalStepLimit(1000000);manager.initialize(initial);
    bool completed=true;double final_time=0;
    for(int i=1;i<=std::lround(duration*hz);++i){
        const auto& s=manager.integrate(i/hz);emit(s);final_time=s.getTime();
        // This upright fixture only has foot contacts. Do not continue a fall
        // through the floor as though it had torso/hand/knee collision handling.
        if(model.calcMassCenterPosition(s)[1]<.55){completed=false;break;}
    }
    std::ofstream termination(out/"termination.json");
    termination<<"{\"completed\":"<<(completed?"true":"false")<<",\"reason\":";
    str(termination,completed?"requested_horizon_reached":"upright_plant_envelope_exceeded");
    termination<<",\"final_time_s\":";num(termination,final_time);termination<<",\"com_height_guard_m\":0.55,\"guard_evidence\":\"engineering restriction for feet-only collision coverage; not a calibrated fall threshold\"}\n";
    return 0;
}catch(const std::exception& e){std::cerr<<"DYNAMICS_ERROR="<<e.what()<<'\n';return 2;}}
