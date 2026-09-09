// Continuing native articulated plant. Source muscles, joints and inertia remain
// native; canonical geometry is registered by the Python materializer.
#include <OpenSim/OpenSim.h>
#include <OpenSim/Actuators/ModelOperators.h>
#include "native_muscle_metabolism.h"
#include "native_static_pose.h"
#include "native_surface_foundation.h"
#include "native_local_mass_port.h"
#include <filesystem>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>
#include <cmath>
namespace fs=std::filesystem;
using namespace OpenSim;
static void num(std::ostream& o,double v){if(!std::isfinite(v))throw std::runtime_error("nonfinite native state");o<<std::setprecision(17)<<v;}
static void str(std::ostream& o,const std::string& s){o<<'"';for(char c:s){if(c=='"'||c=='\\')o<<'\\';if(c=='\n')o<<"\\n";else o<<c;}o<<'"';}
static void vec(std::ostream& o,const SimTK::Vec3& v){o<<'[';for(int i=0;i<3;i++){if(i)o<<',';num(o,v[i]);}o<<']';}
static void matrix(std::ostream& o,const SimTK::Transform& x){o<<'[';for(int i=0;i<4;i++){if(i)o<<',';o<<'[';for(int j=0;j<4;j++){if(j)o<<',';num(o,i==3?(j==3?1.:0.):j==3?x.p()[i]:x.R().asMat33()(i,j));}o<<']';}o<<']';}
struct Load {std::string body;SimTK::Vec3 station,force;};
class PortForces final:public Force {
 OpenSim_DECLARE_CONCRETE_OBJECT(PortForces,Force);
public:
 std::vector<Load> loads;
 void computeForce(const SimTK::State& s,SimTK::Vector_<SimTK::SpatialVec>& f,SimTK::Vector&) const override {
  for(const auto& p:loads)applyForceToPoint(s,getModel().getBodySet().get(p.body),p.station,p.force,f);
 }
 double power(const SimTK::State& s) const {double p=0;for(const auto& l:loads)p+=dot(l.force,getModel().getBodySet().get(l.body).findStationVelocityInGround(s,l.station));return p;}
};
class ExcitationPorts final:public Controller {
 OpenSim_DECLARE_CONCRETE_OBJECT(ExcitationPorts,Controller);
public:
 std::map<std::string,double> values;
 void computeControls(const SimTK::State&,SimTK::Vector& controls) const override {
  const auto& socket=getSocket<Actuator>("actuators");
  for(int i=0;i<(int)socket.getNumConnectees();i++){const auto& a=socket.getConnectee(i);a.addInControls(SimTK::Vector(1,values.at(a.getName())),controls);}
 }
};
struct Saved {SimTK::State state;std::vector<Load> loads;std::map<std::string,double> excitations,torques;double work,positive_work,metabolic_energy,signed_work,heat_energy;ihm_mass_port::Ledger mass_port;};
int main(int argc,char** argv){try{
 if(argc!=5)throw std::runtime_error("source_dir output_dir environment(free|supine|upright) target_mass required");
 fs::path source=argv[1],out=argv[2];std::string environment=argv[3];
 if(environment!="free"&&environment!="supine"&&environment!="upright")throw std::runtime_error("unknown environment");
 Logger::setLevelString("warn");
 Model model((source/"subject_walk_scaled.osim").string());model.setUseVisualizer(false);model.set_assembly_accuracy(1e-10);
 double original_mass=0;for(const auto& b:model.getComponentList<Body>())original_mass+=b.getMass();
 double target_mass=std::stod(argv[4]);if(!std::isfinite(target_mass)||target_mass<=0)throw std::runtime_error("positive target mass required");
 double mass_scale=target_mass/original_mass;
 for(int i=0;i<model.updBodySet().getSize();i++){auto& b=model.updBodySet().get(i);b.setMass(b.getMass()*mass_scale);auto inertia=b.getInertia();inertia*=mass_scale;b.setInertia(inertia);}
 model.setGravity(environment=="free"?SimTK::Vec3(0):environment=="supine"?SimTK::Vec3(-9.81,0,0):SimTK::Vec3(0,-9.81,0));
 ForceSet passive((source/"subject_walk_scaled_ExpressionBasedCoordinateForceSet.xml").string());
 for(int i=0;i<passive.getSize();i++)model.addComponent(passive.get(i).clone());
 ModelFactory::replacePathsWithFunctionBasedPaths(model,(source/"subject_walk_scaled_FunctionBasedPathSet.xml").string());
 auto* excitation=new ExcitationPorts;excitation->setName("external_excitation_ports");
 for(const auto& m:model.getComponentList<Muscle>()){
  excitation->addActuator(m);auto* millard=dynamic_cast<const Millard2012EquilibriumMuscle*>(&m);
  auto* thelen=dynamic_cast<const Thelen2003Muscle*>(&m);
  excitation->values[m.getName()]=millard?millard->getDefaultActivation():thelen?thelen->getDefaultActivation():.01;
 }
 model.addController(excitation);
 // The source model declares CoordinateActuator torque ports for the lumbar and
 // for both arms (shoulder x3, elbow, forearm pronation), and until now NO
 // controller was connected to them, so their controls were identically zero for
 // the whole run.  The arm chain was therefore a passive rag doll: measured, it
 // reached 11.0 rad/s on pro_sup_r -- a coordinate the source passive force set
 // does not damp at all -- and dragged arm26_BIClong_l to 3.05 m/s of fiber
 // velocity against a 1.28 m/s maximum contraction velocity.  Outside the
 // force-velocity curve's domain the error controller collapses its step, which
 // is what stalls one 10 ms advance for minutes.  These ports are torque
 // actuators, NOT muscles: they are engineering drive on declared actuators, and
 // are reported separately from `muscles` so nothing can count them as muscle.
 auto* torque=new ExcitationPorts;torque->setName("declared_coordinate_actuator_ports");
 for(const auto& a:model.getComponentList<CoordinateActuator>()){torque->addActuator(a);torque->values[a.getName()]=0;}
 if(!torque->values.empty())model.addController(torque);
 auto* external=new PortForces;external->setName("external_material_point_forces");model.addForce(external);
 IHMNativeMuscleMetabolism metabolism(model);
 model.finalizeConnections();auto initial=model.initSystem();model.realizePosition(initial);
 ForceSet contact_templates((source/"subject_walk_scaled_ContactForceSet.xml").string());
 ihm_surface::Foundation* surface_foundation=nullptr;
 double support_plane=0;std::map<std::string,double> proxy_radius;
 if(environment=="upright"){
  ContactGeometrySet geo((source/"subject_walk_scaled_ContactGeometrySet.xml").string());
  for(int i=0;i<geo.getSize();i++)model.addContactGeometry(geo.get(i).clone());
  for(int i=0;i<contact_templates.getSize();i++)model.addComponent(contact_templates.get(i).clone());
  // The walking source only supplies foot contact. A falling trunk previously
  // passed through the floor until pelvis_ty hit its -1 m coordinate clamp.
  // Give the other segments unilateral contact, without anchoring their pose.
  // These inscribed COM spheres are engineering proxies inferred from inertia;
  // they are not an anatomical skin reconstruction or a balance controller.
  auto& floor=dynamic_cast<ContactHalfSpace&>(model.updContactGeometrySet().get("floor"));
  for(const auto& b:model.getComponentList<Body>()){
   // Talus has a tiny placeholder mass/inertia, not a geometric ellipsoid;
   // source foot contacts already represent the articulated foot assembly.
   if(b.getName().find("talus_")==0||b.getName().find("calcn_")==0||b.getName().find("toes_")==0)continue;
   const auto moments=b.getInertia().getMoments();double radius2=1e10;
   for(int k=0;k<3;k++)radius2=std::min(radius2,5*(moments[(k+1)%3]+moments[(k+2)%3]-moments[k])/(2*b.getMass()));
   if(!(radius2>0))throw std::runtime_error("source inertia cannot support whole-body contact radius");
   auto* sphere=new ContactSphere(std::sqrt(radius2),b.getMassCenter(),b);sphere->setName("fall_proxy_"+b.getName());model.addContactGeometry(sphere);
   auto* force=dynamic_cast<SmoothSphereHalfSpaceForce*>(contact_templates.get(0).clone());
   force->setName("fall_support_"+b.getName());force->connectSocket_sphere(*sphere);force->connectSocket_half_space(floor);model.addComponent(force);
  }
 }else if(environment=="supine" && fs::exists(source/"supine_surface_foundation.txt")){
  surface_foundation=new ihm_surface::Foundation;surface_foundation->setName("retained_skin_surface_foundation");
  surface_foundation->read((source/"supine_surface_foundation.txt").string(),model);
  if(fs::exists(source/"surface_sensor_indices.txt"))surface_foundation->select((source/"surface_sensor_indices.txt").string());
  if(fs::exists(source/"bed_compression.txt"))surface_foundation->readBed((source/"bed_compression.txt").string());
  support_plane=surface_foundation->plane;
  model.addForce(surface_foundation);
 }else if(environment=="supine"){
  // Engineering posterior contact proxies: uniform ellipsoid posterior semi-axis
  // derived from source mass/inertia, placed at retained segment COM. Contact
  // material law is explicitly transferred from the source foot model, not a
  // calibrated mattress or surface-anatomy reconstruction.
  support_plane=1e10;
  for(const auto& b:model.getComponentList<Body>()){
   auto moments=b.getInertia().getMoments();double radius2=5*(moments[1]+moments[2]-moments[0])/(2*b.getMass());
   if(!(radius2>0))throw std::runtime_error("source inertia cannot support ellipsoid posterior radius");
   double radius=std::sqrt(radius2);proxy_radius[b.getName()]=radius;
   support_plane=std::min(support_plane,b.findStationLocationInGround(initial,b.getMassCenter())[0]-radius);
  }
  auto* plane=new ContactHalfSpace(SimTK::Vec3(support_plane,0,0),SimTK::Vec3(0,0,SimTK::Pi),model.getGround());plane->setName("supine_plane");model.addContactGeometry(plane);
  for(const auto& b:model.getComponentList<Body>()){
   auto* sphere=new ContactSphere(proxy_radius.at(b.getName()),b.getMassCenter(),b);sphere->setName("posterior_"+b.getName());model.addContactGeometry(sphere);
   auto* force=dynamic_cast<SmoothSphereHalfSpaceForce*>(contact_templates.get(0).clone());
   force->setName("support_"+b.getName());force->connectSocket_sphere(*sphere);force->connectSocket_half_space(*plane);model.addComponent(force);
  }
 }
 model.finalizeConnections();SimTK::State state=model.initSystem();state.setTime(0);
 // Material registration belongs to the retained source pose. An optional
 // initialization must move that same body, not silently redefine its geometry.
 model.realizePosition(state);std::ostringstream registration_reference;registration_reference<<'{';bool first_reference=true;
 for(const auto& b:model.getComponentList<Body>()){
  if(!first_reference)registration_reference<<',';first_reference=false;str(registration_reference,b.getName());registration_reference<<":{\"transform_ground\":";
  matrix(registration_reference,b.getTransformInGround(state));registration_reference<<",\"mass_center_local_m\":";vec(registration_reference,b.getMassCenter());registration_reference<<'}';
 }registration_reference<<'}';
 const bool initial_pose_applied=fs::exists(source/"initial_pose.txt");
 if(initial_pose_applied){
  std::ifstream pose(source/"initial_pose.txt");std::string schema;int count;pose>>schema>>count;
  if(!pose||schema!="IHM_INITIAL_POSE_V1"||count<1||count>model.getCoordinateSet().getSize())throw std::runtime_error("invalid initial pose schema/count");
  std::map<std::string,double> coordinates;
  for(int i=0;i<count;i++){
   std::string name;double value;pose>>name>>value;
   if(!pose||!std::isfinite(value)||coordinates.count(name))throw std::runtime_error("invalid initial pose coordinate");
   const auto& c=model.getCoordinateSet().get(name);
   if(c.isDependent(state)||c.getLocked(state)||c.isPrescribed(state)||value<c.getRangeMin()||value>c.getRangeMax())throw std::runtime_error("initial pose requires independent free coordinates within source bounds");
   coordinates[name]=value;c.setValue(state,value,false);
  }
  std::string extra;if(pose>>extra)throw std::runtime_error("trailing initial pose data");
  model.assemble(state);state.updU()=0.;
  for(const auto& item:coordinates)if(std::abs(model.getCoordinateSet().get(item.first).getValue(state)-item.second)>1e-8)throw std::runtime_error("assembly changed requested initial pose");
 }
 for(const auto& m:model.getComponentList<Muscle>())m.setActivation(state,excitation->values.at(m.getName()));
 model.equilibrateMuscles(state);model.realizeDynamics(state);model.print((out/"assembled_model.osim").string());
 const char* mass_mode=std::getenv("IHM_INSTANCE_MASS_MODE");
 const bool mass_enabled=mass_mode&&std::string(mass_mode)=="1";
 const char* mass_reference=std::getenv("IHM_MASS_REFERENCE_ID");
 ihm_mass_port::Ledger mass_port(mass_enabled,mass_enabled&&mass_reference?mass_reference:"");
 const auto reference=metabolism.sample(state);
 double work=0,positive_work=0,metabolic_energy=0,signed_work=0,heat_energy=0;std::map<std::string,Saved> checkpoints;
 auto active_power=[&](){model.realizeDynamics(state);double power=0;for(const auto& m:model.getComponentList<Muscle>())power+=std::max(0.,-m.getActiveFiberForce(state)*m.getFiberVelocity(state));return power;};
 auto emit=[&](const std::string& kind){
  model.realizeAcceleration(state);auto metabolic=metabolism.sample(state);std::ostringstream o;o<<"{\"kind\":";str(o,kind);o<<",\"time_s\":";num(o,state.getTime());
  o<<",\"mass_transfer\":"<<mass_port.json();
  o<<",\"initial_pose_applied\":"<<(initial_pose_applied?"true":"false");
  o<<",\"registration_reference_bodies\":"<<registration_reference.str();
  o<<",\"mass_kg\":";num(o,model.getTotalMass(state));o<<",\"gravity_m_s2\":";vec(o,model.getGravity());
  o<<",\"kinetic_energy_j\":";num(o,model.calcKineticEnergy(state));o<<",\"potential_energy_j\":";num(o,model.calcPotentialEnergy(state));o<<",\"external_work_j\":";num(o,work);o<<",\"positive_active_fiber_work_j\":";num(o,positive_work);o<<",\"external_power_w\":";num(o,external->power(state));
  o<<",\"muscle_metabolic_energy_j\":";num(o,metabolic_energy);o<<",\"total_muscle_metabolic_w\":";num(o,metabolic.total_muscle_metabolic_w);o<<",\"signed_active_fiber_power_w\":";num(o,metabolic.active_fiber_work_w);o<<",\"muscle_heat_w\":";num(o,metabolic.muscle_heat_w);o<<",\"metabolic_analysis_mass_kg\":";num(o,metabolic.analysis_mass_kg);
  o<<",\"signed_active_fiber_work_j\":";num(o,signed_work);o<<",\"muscle_heat_energy_j\":";num(o,heat_energy);
  o<<",\"metabolic_reference\":{\"M0_w\":";num(o,reference.total_muscle_metabolic_w);o<<",\"W0_w\":";num(o,reference.active_fiber_work_w);o<<",\"H0_w\":";num(o,reference.total_muscle_metabolic_w-reference.active_fiber_work_w);o<<'}';
  o<<",\"constraint_position_error\":";num(o,state.getQErr().norm());o<<",\"constraint_velocity_error\":";num(o,state.getUErr().norm());
  o<<",\"original_source_mass_kg\":";num(o,original_mass);o<<",\"mass_scale\":";num(o,mass_scale);
  o<<",\"contact_model\":";str(o,surface_foundation?"retained_skin_foundation":environment=="upright"?"source_feet_and_inertia_inscribed_body_spheres":"source_sphere_proxies");
  o<<",\"support_plane_source_x_m\":";num(o,support_plane);bool first=true;
  std::ostringstream bodies;bodies<<'{';first=true;
  SimTK::Vector_<SimTK::SpatialVec> reactions;model.getMatterSubsystem().calcMobilizerReactionForces(state,reactions);
  for(const auto& b:model.getComponentList<Body>()){
   if(!first)bodies<<',';first=false;str(bodies,b.getName());bodies<<":{\"transform_ground\":";matrix(bodies,b.getTransformInGround(state));
   const auto effective=b.getMobilizedBody().getBodyMassProperties(state);
   bodies<<",\"mass_properties_basis\":\"effective Simbody State instance\",\"mass_kg\":";num(bodies,effective.getMass());bodies<<",\"mass_center_local_m\":";vec(bodies,effective.getMassCenter());bodies<<",\"inertia_moments_kg_m2\":";vec(bodies,effective.calcCentralInertia().getMoments());bodies<<",\"inertia_products_kg_m2\":";vec(bodies,effective.calcCentralInertia().getProducts());
   bodies<<",\"model_baseline_mass_properties\":{\"mass_kg\":";num(bodies,b.getMass());bodies<<",\"mass_center_local_m\":";vec(bodies,b.getMassCenter());bodies<<",\"inertia_moments_kg_m2\":";vec(bodies,b.getInertia().getMoments());bodies<<",\"inertia_products_kg_m2\":";vec(bodies,b.getInertia().getProducts());bodies<<'}';
   auto v=b.getMobilizedBody().getBodyVelocity(state);bodies<<",\"angular_velocity_rad_s\":";vec(bodies,v[0]);bodies<<",\"origin_velocity_m_s\":";vec(bodies,v[1]);
   auto reaction=reactions[b.getMobilizedBodyIndex()];bodies<<",\"joint_reaction_force_n\":";vec(bodies,reaction[1]);bodies<<",\"joint_reaction_moment_nm\":";vec(bodies,reaction[0]);bodies<<'}';
  }bodies<<'}';
  o<<",\"bodies\":"<<bodies.str();
  o<<",\"muscles\":{";first=true;
  for(const auto& m:model.getComponentList<Muscle>()){
   if(!first)o<<',';first=false;str(o,m.getName());o<<":{\"path_length_m\":";num(o,m.getLength(state));o<<",\"fiber_length_m\":";num(o,m.getFiberLength(state));o<<",\"fiber_velocity_m_s\":";num(o,m.getFiberVelocity(state));
   o<<",\"optimal_fiber_length_m\":";num(o,m.getOptimalFiberLength());o<<",\"tendon_force_n\":";num(o,m.getTendonForce(state));o<<",\"max_isometric_force_n\":";num(o,m.getMaxIsometricForce());
   o<<",\"activation\":";num(o,m.getActivation(state));o<<",\"excitation\":";num(o,m.getExcitation(state));o<<",\"muscle_type\":";str(o,m.getConcreteClassName());o<<",\"metabolic_power_w\":";num(o,metabolic.muscles.at(m.getName()).metabolic_w);o<<",\"metabolic_roundoff_tolerance_w\":";num(o,metabolic.muscles.at(m.getName()).roundoff_tolerance_w);o<<",\"metabolic_analysis_mass_kg\":";num(o,metabolic.muscles.at(m.getName()).analysis_mass_kg);o<<",\"sensor_basis\":\"native muscle states; retained source80 fitted paths or explicitly registered added geometry paths\"}";
  }o<<"},\"coordinate_actuators\":{";first=true;
  for(const auto& a:model.getComponentList<CoordinateActuator>()){
   if(!first)o<<',';first=false;str(o,a.getName());o<<":{\"coordinate\":";str(o,a.getCoordinate()->getName());
   o<<",\"command\":";num(o,torque->values.at(a.getName()));o<<",\"optimal_force_n_m\":";num(o,a.getOptimalForce());
   o<<",\"torque_n_m\":";num(o,torque->values.at(a.getName())*a.getOptimalForce());
   o<<",\"basis\":\"declared source CoordinateActuator driven by an explicit engineering torque port; not a muscle\"}";
  }o<<"},\"coordinates\":{";first=true;
  for(const auto& c:model.getComponentList<Coordinate>()){
   if(!first)o<<',';first=false;str(o,c.getName());o<<":{\"value\":";num(o,c.getValue(state));o<<",\"speed\":";num(o,c.getSpeedValue(state));o<<",\"unit\":";str(o,c.getMotionType()==Coordinate::Rotational?"rad":"m");o<<'}';
  }o<<"},\"contacts\":[";first=true;SimTK::Vec3 total(0),foot_r(0),foot_l(0);
  for(const auto& f:model.getComponentList<SmoothSphereHalfSpaceForce>()){
   if(!first)o<<',';first=false;auto wrench=f.getSphereForce(state);total+=wrench[1];
   const auto& sphere=f.getConnectee<ContactSphere>("sphere");auto center=sphere.getFrame().findStationLocationInGround(state,sphere.get_location());
   o<<"{\"name\":";str(o,f.getName());o<<",\"body_frame\":";str(o,sphere.getFrame().getName());o<<",\"force_n\":";vec(o,wrench[1]);o<<",\"moment_nm\":";vec(o,wrench[0]);
   o<<",\"paired_force_residual_n\":";vec(o,wrench[1]+f.getHalfSpaceForce(state)[1]);o<<",\"center_m\":";vec(o,center);o<<",\"radius_m\":";num(o,sphere.getRadius());o<<'}';
   auto n=sphere.getFrame().getName();if(n=="calcn_r"||n=="toes_r")foot_r+=wrench[1];if(n=="calcn_l"||n=="toes_l")foot_l+=wrench[1];
  }
  if(surface_foundation){
   const auto surface=surface_foundation->sample(state);total+=surface.force;
   for(const auto& item:surface.bodies){
    if(!first)o<<',';first=false;
    o<<"{\"name\":";str(o,"surface_"+item.first);o<<",\"body_frame\":";str(o,item.first);
    o<<",\"geometry_type\":\"retained_skin_foundation\",\"force_n\":";vec(o,item.second.force);
    o<<",\"moment_nm\":";vec(o,item.second.moment);o<<",\"penetration_m\":";num(o,item.second.maximum_penetration);
    o<<",\"contacting_points\":"<<item.second.contacting_points<<'}';
    if(item.first=="calcn_r"||item.first=="toes_r")foot_r+=item.second.force;
    if(item.first=="calcn_l"||item.first=="toes_l")foot_l+=item.second.force;
   }
  }
  SimTK::Vec3 ext(0);for(const auto& p:external->loads)ext+=p.force;
  o<<"]";
  if(surface_foundation){const auto surface=surface_foundation->sample(state);
   o<<",\"surface_foundation\":{\"elastic_energy_j\":";num(o,surface.energy);
   o<<",\"skin_elastic_energy_j\":";num(o,surface.skin_energy);o<<",\"bed_elastic_energy_j\":";num(o,surface.bed_energy);
   o<<",\"maximum_bed_indentation_m\":";num(o,surface.maximum_bed_indentation);o<<",\"maximum_total_approach_m\":";num(o,surface.maximum_total_approach);
   o<<",\"normal_rate_law\":";str(o,surface_foundation->bed_enabled?"conservative_measured_bed_skin_series":"transferred_skin_contact_dissipation");
   o<<",\"power_to_body_w\":";num(o,surface.power);o<<",\"dissipative_power_w\":";num(o,surface.dissipative_power);
   o<<",\"bed_force_n\":";vec(o,-surface.force);o<<",\"bed_moment_about_source_origin_nm\":";vec(o,surface.bed_moment);
   o<<",\"maximum_penetration_m\":";num(o,surface.maximum_penetration);
   o<<",\"sensor_points\":[";bool sensor_first=true;
   for(const auto& point:surface.observations){if(!sensor_first)o<<',';sensor_first=false;
    o<<"{\"quadrature_index\":"<<point.index<<",\"body_frame\":";str(o,point.body);
    o<<",\"point_source_m\":";vec(o,point.location);o<<",\"normal_source\":[-1,0,0],\"reaction_normal_source\":[1,0,0],\"force_n\":";vec(o,point.force);
    o<<",\"indentation_m\":";num(o,point.indentation);o<<",\"bed_indentation_m\":";num(o,point.bed_indentation);o<<",\"total_approach_m\":";num(o,point.total_approach);o<<",\"contact_area_m2\":";num(o,point.area);o<<'}';}
   o<<"]}";}
  o<<",\"contact_force_n\":";vec(o,total);o<<",\"momentum_balance_residual_n\":";vec(o,model.getTotalMass(state)*(model.calcMassCenterAcceleration(state)-model.getGravity())-total-ext);
  int axis=environment=="supine"?0:1;o<<",\"foot_contact_force_n\":{\"r\":";num(o,std::max(0.,foot_r[axis]));o<<",\"l\":";num(o,std::max(0.,foot_l[axis]));o<<"},\"environment\":";str(o,environment);o<<'}';
  std::cout<<"@IHM "<<o.str()<<std::endl;
 };
 emit("initialized");std::string line;
 while(std::getline(std::cin,line)){
  Saved before{state,external->loads,excitation->values,torque->values,work,positive_work,metabolic_energy,signed_work,heat_energy,mass_port};
  try{
   std::istringstream in(line);std::string command;in>>command;
   if(command=="close")break;
   if(command=="mass_transfer"){mass_port.apply(model,state,in);emit("mass_transferred");continue;}
   if(command=="moment_arms"){
    int n,m;in>>n;if(!in||n<1||n>1000)throw std::runtime_error("invalid moment arm muscle count");
    std::vector<std::string> muscles(n);for(auto& name:muscles){in>>name;if(!excitation->values.count(name))throw std::runtime_error("unknown moment arm muscle");}
    in>>m;if(!in||m<1||m>1000)throw std::runtime_error("invalid moment arm coordinate count");
    std::vector<std::string> coordinates(m);for(auto& name:coordinates){in>>name;model.getCoordinateSet().get(name);}
    std::string extra;if(!in||in>>extra)throw std::runtime_error("invalid moment arm query");
    model.realizePosition(state);std::ostringstream o;o<<"{\"kind\":\"moment_arms\",\"time_s\":";num(o,state.getTime());o<<",\"moment_arms_m\":{";
    bool first_muscle=true;for(const auto& name:muscles){if(!first_muscle)o<<',';first_muscle=false;str(o,name);o<<":{";bool first_coordinate=true;
     for(const auto& coordinate:coordinates){if(!first_coordinate)o<<',';first_coordinate=false;str(o,coordinate);o<<':';num(o,model.getMuscles().get(name).computeMomentArm(state,model.updCoordinateSet().get(coordinate)));}o<<'}';}
    o<<"}}";std::cout<<"@IHM "<<o.str()<<std::endl;continue;
   }
   if(command=="body_point"){
    std::string name,extra;SimTK::Vec3 station;in>>name;for(int k=0;k<3;++k)in>>station[k];
    if(!in||in>>extra||!station.isFinite())throw std::runtime_error("invalid body point query");
    model.realizeVelocity(state);const auto& body=model.getBodySet().get(name).getMobilizedBody();
    std::ostringstream o;o<<"{\"kind\":\"body_point\",\"body\":";str(o,name);o<<",\"station_m\":";vec(o,station);o<<",\"time_s\":";num(o,state.getTime());o<<",\"point_source_m\":";vec(o,body.findStationLocationInGround(state,station));o<<",\"velocity_source_m_s\":";vec(o,body.findStationVelocityInGround(state,station));o<<'}';std::cout<<"@IHM "<<o.str()<<std::endl;continue;
   }
   if(command=="observe"){emit("observed");continue;}
   if(command=="checkpoint"){std::string key;in>>key;if(key.empty()||checkpoints.size()>=64||checkpoints.count(key))throw std::runtime_error("invalid checkpoint id/capacity");checkpoints.emplace(key,before);std::cout<<"@IHM {\"kind\":\"checkpointed\"}"<<std::endl;continue;}
   if(command=="restore"){std::string key;in>>key;const auto& old=checkpoints.at(key);state=old.state;external->loads=old.loads;excitation->values=old.excitations;torque->values=old.torques;work=old.work;positive_work=old.positive_work;metabolic_energy=old.metabolic_energy;signed_work=old.signed_work;heat_energy=old.heat_energy;mass_port=old.mass_port;model.markControlsAsInvalid(state);state.invalidateAllCacheAtOrAbove(SimTK::Stage::Instance);emit("restored");continue;}
   if(command=="drop"){std::string key;in>>key;if(!checkpoints.erase(key))throw std::runtime_error("unknown checkpoint");std::cout<<"@IHM {\"kind\":\"dropped\"}"<<std::endl;continue;}
   if(command=="evaluate_static_pose"){const auto payload=ihm_static_pose::evaluate(model,state,in,environment,!external->loads.empty(),support_plane,surface_foundation);std::cout<<"@IHM "<<payload<<std::endl;continue;}
   if(command!="advance")throw std::runtime_error("unknown command");double dt;int count;in>>dt>>count;if(!in||!std::isfinite(dt)||dt<=0||dt>.02||count<0||count>10000)throw std::runtime_error("invalid native step/force count");
   model.realizePosition(state);external->loads.clear();
   for(int i=0;i<count;i++){std::string name;SimTK::Vec3 point,force;in>>name;for(int k=0;k<3;k++)in>>point[k];for(int k=0;k<3;k++)in>>force[k];if(!in||!point.isFinite()||!force.isFinite())throw std::runtime_error("invalid force port");const auto& body=model.getBodySet().get(name);external->loads.push_back({name,~body.getTransformInGround(state)*point,force});}
   in>>count;if(!in||count<0||count>1000)throw std::runtime_error("invalid excitation count");
   for(int i=0;i<count;i++){std::string name;double value;in>>name>>value;if(!in||!std::isfinite(value)||value<0||value>1||!excitation->values.count(name))throw std::runtime_error("invalid muscle excitation");excitation->values[name]=value;}
   if(in>>count){
    if(count<0||count>1000)throw std::runtime_error("invalid coordinate actuation count");
    for(int i=0;i<count;i++){std::string name;double value;in>>name>>value;if(!in||!std::isfinite(value)||value<-1||value>1||!torque->values.count(name))throw std::runtime_error("invalid coordinate actuation");torque->values[name]=value;}
   }
   in.clear();std::string extra;if(in>>extra)throw std::runtime_error("trailing command data");
   model.markControlsAsInvalid(state);state.invalidateAllCacheAtOrAbove(SimTK::Stage::Dynamics);model.realizeVelocity(state);double p0=external->power(state),active0=active_power();const auto metabolic0=metabolism.sample(state);
   Manager manager(model);manager.setIntegratorAccuracy(1e-7);manager.setIntegratorConstraintTolerance(1e-9);manager.setIntegratorMaximumStepSize(.0005);manager.setIntegratorInternalStepLimit(100000);manager.initialize(state);
   state=manager.integrate(state.getTime()+dt);model.realizeVelocity(state);work+=dt*.5*(p0+external->power(state));positive_work+=dt*.5*(active0+active_power());const auto metabolic1=metabolism.sample(state);
   // One endpoint quadrature owns chemical energy and signed active-fiber work.
   metabolic_energy+=dt*.5*(metabolic0.total_muscle_metabolic_w+metabolic1.total_muscle_metabolic_w);
   signed_work+=dt*.5*(metabolic0.active_fiber_work_w+metabolic1.active_fiber_work_w);
   heat_energy=metabolic_energy-signed_work;emit("advanced");
  }catch(const std::exception& error){state=before.state;external->loads=before.loads;excitation->values=before.excitations;torque->values=before.torques;work=before.work;positive_work=before.positive_work;metabolic_energy=before.metabolic_energy;signed_work=before.signed_work;heat_energy=before.heat_energy;mass_port=before.mass_port;state.invalidateAllCacheAtOrAbove(SimTK::Stage::Instance);model.markControlsAsInvalid(state);std::ostringstream o;o<<"{\"error\":";str(o,error.what());o<<'}';std::cout<<"@IHM "<<o.str()<<std::endl;}
 }
 return 0;
}catch(const std::exception& e){std::cerr<<"MECHANICAL_STREAM_ERROR="<<e.what()<<std::endl;return 2;}}
