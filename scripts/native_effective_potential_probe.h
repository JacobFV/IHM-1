#pragma once
// Numerical STATIC merit only. Never exported into physiological energy ledgers.
#include <OpenSim/OpenSim.h>
#include <lepton/Parser.h>
#include <lepton/ParsedExpression.h>
#include <functional>
#include "native_static_pose.h"
namespace ihm_effective {
using ihm_static_pose::number;using ihm_static_pose::quoted;using ihm_static_pose::array;
inline double integral(const std::function<double(double)>& f,double a,double b){
 if(a==b)return 0;if(a>b)return -integral(f,b,a);
 std::function<double(double,double,double,double,double,double,double,int)> recurse;
 recurse=[&](double l,double r,double fl,double fm,double fr,double whole,double tol,int depth){
  const double m=.5*(l+r),u=f(.5*(l+m)),v=f(.5*(m+r));
  const double left=(m-l)*(fl+4*u+fm)/6,right=(r-m)*(fm+4*v+fr)/6,error=left+right-whole;
  if(!std::isfinite(left+right))throw std::runtime_error("nonfinite static primitive");
  if(std::abs(error)<=15*tol)return left+right+error/15;
  if(depth==0)throw std::runtime_error("static primitive quadrature cap");
  return recurse(l,m,fl,u,fm,left,tol/2,depth-1)+recurse(m,r,fm,v,fr,right,tol/2,depth-1);
 };
 const double fa=f(a),fm=f(.5*(a+b)),fb=f(b);return recurse(a,b,fa,fm,fb,(b-a)*(fa+4*fm+fb)/6,1e-11,14);
}
inline std::string evaluate(OpenSim::Model& model,const SimTK::State& continuing,std::istream& input,const ihm_surface::Foundation* foundation){
 int count;if(!(input>>count)||count<1||count>31)throw std::runtime_error("invalid effective pose count");
 std::map<std::string,double> requested;
 for(int i=0;i<count;i++){
  std::string name;double q;if(!(input>>name>>q)||!std::isfinite(q)||requested.count(name))throw std::runtime_error("invalid effective pose input");
  const auto& c=model.getCoordinateSet().get(name);
  if(c.isDependent(continuing)||c.getLocked(continuing)||c.isPrescribed(continuing)||q<c.getRangeMin()||q>c.getRangeMax())throw std::runtime_error("effective pose outside free source coordinates");
  requested[name]=q;
 }
 std::string extra;if(input>>extra)throw std::runtime_error("trailing effective pose input");
 SimTK::State s=continuing;s.updU()=0.;for(const auto& item:requested)model.getCoordinateSet().get(item.first).setValue(s,item.second,false);
 model.assemble(s);s.updU()=0.;model.markControlsAsInvalid(s);model.equilibrateMuscles(s);model.realizeAcceleration(s);
 const auto& matter=model.getMatterSubsystem();const auto& system=model.getMultibodySystem();
 if(s.getNQ()!=s.getNU())throw std::runtime_error("non-square native coordinate velocity map");
 double nerror=0;for(int i=0;i<s.getNU();i++){SimTK::Vector u(s.getNU(),0.),qd;u[i]=1;matter.multiplyByN(s,false,u,qd);nerror=std::max(nerror,(qd-u).norm());}
 if(nerror>1e-12)throw std::runtime_error("nonidentity coordinate velocity map requires explicit pullback");
 SimTK::Vector tree,constrained,zero(s.getNU(),0.);matter.calcResidualForceIgnoringConstraints(s,system.getMobilityForces(s,SimTK::Stage::Dynamics),system.getRigidBodyForces(s,SimTK::Stage::Dynamics),zero,tree);
 matter.calcResidualForce(s,system.getMobilityForces(s,SimTK::Stage::Dynamics),system.getRigidBodyForces(s,SimTK::Stage::Dynamics),zero,s.getMultipliers(),constrained);
 std::vector<std::string> mobility(s.getNU());std::vector<const OpenSim::Coordinate*> independent;
 for(const auto& c:model.getComponentList<OpenSim::Coordinate>()){
  const auto& body=matter.getMobilizedBody(c.getBodyIndex());const int i=int(body.getFirstUIndex(s))+int(c.getMobilizerQIndex());
  if(i<0||i>=s.getNU()||!mobility[i].empty())throw std::runtime_error("ambiguous native coordinate map");mobility[i]=c.getName();
  if(!c.isDependent(s)&&!c.getLocked(s)&&!c.isPrescribed(s)){
   if(c.getValue(s)<c.getRangeMin()-1e-10||c.getValue(s)>c.getRangeMax()+1e-10)throw std::runtime_error("assembled effective q outside bounds");independent.push_back(&c);
  }
 }
 double gravity=0,joints=0,passive=0,active=0;
 SimTK::Vector_<SimTK::SpatialVec> gravity_bodies(matter.getNumBodies(),SimTK::SpatialVec(SimTK::Vec3(0),SimTK::Vec3(0))),contact_bodies=gravity_bodies;
 SimTK::Vector gravity_forces,contact_forces,joint_forces(s.getNU(),0.);
 for(const auto& b:model.getComponentList<OpenSim::Body>()){
  const auto com=b.findStationLocationInGround(s,b.getMassCenter());const auto force=b.getMass()*model.getGravity();gravity-=SimTK::dot(force,com);
  gravity_bodies[b.getMobilizedBodyIndex()]=SimTK::SpatialVec((com-b.getTransformInGround(s).p())%force,force);
 }
 matter.multiplyBySystemJacobianTranspose(s,gravity_bodies,gravity_forces);
 const auto surface=foundation?foundation->sample(s,false):ihm_surface::Sample{};
 for(const auto& item:surface.bodies){const auto& b=model.getBodySet().get(item.first);contact_bodies[b.getMobilizedBodyIndex()]=SimTK::SpatialVec(item.second.moment,item.second.force);}
 matter.multiplyBySystemJacobianTranspose(s,contact_bodies,contact_forces);
 std::ostringstream joint_json;joint_json<<'[';bool first=true;
 for(const auto& force:model.getComponentList<OpenSim::ExpressionBasedCoordinateForce>()){
  auto program=Lepton::Parser::parse(force.get_expression()).createProgram();const auto& c=model.getCoordinateSet().get(force.get_coordinate());
  const auto f=[&](double q){std::map<std::string,double> vars{{"q",q},{"qdot",0}};return program.evaluate(vars);};
  const double energy=-integral(f,0,c.getValue(s));joints+=energy;
  joint_forces[int(matter.getMobilizedBody(c.getBodyIndex()).getFirstUIndex(s))+int(c.getMobilizerQIndex())]+=f(c.getValue(s));
  if(!first)joint_json<<',';first=false;joint_json<<"{\"coordinate\":";quoted(joint_json,c.getName());joint_json<<",\"energy_j\":";number(joint_json,energy);joint_json<<",\"force\":";number(joint_json,f(c.getValue(s)));joint_json<<'}';
 }joint_json<<']';
 for(const auto& actuator:model.getComponentList<OpenSim::CoordinateActuator>())if(std::abs(actuator.getActuation(s))>1e-12)throw std::runtime_error("nonzero coordinate actuator outside effective merit");
 std::ostringstream muscles;muscles<<'[';first=true;int muscle_count=0;
 for(const auto& muscle:model.getComponentList<OpenSim::Muscle>()){
  const double activation=muscle.getActivation(s),lf=muscle.getFiberLength(s),lo=muscle.getOptimalFiberLength(),x=lf/lo;
  double minimum=0,primitive=0;
  if(const auto* thelen=dynamic_cast<const OpenSim::Thelen2003Muscle*>(&muscle)){
   minimum=thelen->getMinimumFiberLength();const double k=thelen->get_KshapeActive();primitive=std::sqrt(SimTK::Pi*k)/2*std::erf((x-1)/std::sqrt(k));
  }else if(const auto* millard=dynamic_cast<const OpenSim::Millard2012EquilibriumMuscle*>(&muscle)){
   minimum=millard->getMinimumFiberLength();primitive=integral([&](double y){return millard->getActiveForceLengthCurve().calcValue(y);},1.,x);
  }else throw std::runtime_error("unsupported muscle effective primitive");
  const double active_energy=muscle.getMaxIsometricForce()*activation*lo*primitive;
  const double passive_energy=muscle.getMusclePotentialEnergy(s);active+=active_energy;passive+=passive_energy;
  if(!first)muscles<<',';first=false;muscle_count++;muscles<<"{\"name\":";quoted(muscles,muscle.getName());muscles<<",\"type\":";quoted(muscles,muscle.getConcreteClassName());
  const auto field=[&](const char* key,double value){muscles<<',';quoted(muscles,key);muscles<<':';number(muscles,value);};
  field("length_m",muscle.getLength(s));field("fiber_length_m",lf);field("minimum_fiber_length_m",minimum);field("fiber_velocity_m_s",muscle.getFiberVelocity(s));field("activation",activation);
  field("cos_pennation",std::cos(muscle.getPennationAngle(s)));field("normalized_fiber_length",x);field("internal_stiffness_n_m",(muscle.getTendonStiffness(s)+muscle.getFiberStiffnessAlongTendon(s))/std::pow(std::cos(muscle.getPennationAngle(s)),2));field("tendon_force_n",muscle.getTendonForce(s));field("fiber_force_along_tendon_n",muscle.getFiberForceAlongTendon(s));field("fiso_n",muscle.getMaxIsometricForce());
  field("passive_energy_j",passive_energy);field("active_effective_energy_j",active_energy);muscles<<",\"moment_arms_m\":[";
  for(size_t i=0;i<independent.size();i++){if(i)muscles<<',';number(muscles,muscle.computeMomentArm(s,*independent[i]));}muscles<<"]}";
 }muscles<<']';if(muscle_count!=98)throw std::runtime_error("effective probe needs98 muscles");
 const double skin=surface.skin_energy,bed=surface.bed_energy;
 std::ostringstream out;out<<"{\"kind\":\"effective_pose_evaluated\",\"numerical_merit_only\":true,\"physical_time_advanced_s\":0,\"continuing_state_unchanged\":true,\"coordinate_velocity_map_error\":";number(out,nerror);
 out<<",\"gravity_generalized_forces\":";array(out,gravity_forces);out<<",\"contact_generalized_forces\":";array(out,contact_forces);out<<",\"joint_generalized_forces\":";array(out,joint_forces);
 out<<",\"q\":";array(out,s.getQ());out<<",\"udot\":";array(out,s.getUDot());out<<",\"tree_residual\":";array(out,tree);out<<",\"constrained_residual\":";array(out,constrained);
 out<<",\"mobility_names\":[";for(size_t i=0;i<mobility.size();i++){if(i)out<<',';quoted(out,mobility[i]);}out<<"],\"independent_names\":[";
 for(size_t i=0;i<independent.size();i++){if(i)out<<',';quoted(out,independent[i]->getName());}out<<"],\"independent_q\":[";for(size_t i=0;i<independent.size();i++){if(i)out<<',';number(out,independent[i]->getValue(s));}out<<']';
 const auto field=[&](const char* key,double value){out<<',';quoted(out,key);out<<':';number(out,value);};
 field("gravity_energy_j",gravity);field("joint_energy_j",joints);field("skin_energy_j",skin);field("bed_energy_j",bed);field("muscle_passive_energy_j",passive);field("muscle_active_effective_energy_j",active);field("effective_energy_j",gravity+joints+skin+bed+passive+active);
 field("constraint_position_error",s.getQErr().norm());field("constraint_velocity_error",s.getUErr().norm());field("constraint_acceleration_error",s.getUDotErr().norm());out<<",\"joints\":"<<joint_json.str()<<",\"muscles\":"<<muscles.str()<<'}';return out.str();
}
}
