#pragma once
// Isolated 92-muscle fixture commands. Never installed in default stream.
#include <dlfcn.h>
#include <sstream>
#include <iomanip>
namespace ihm_mass_validation {
using namespace SimTK;
inline void scalar(std::ostream& o,double x){if(!std::isfinite(x))throw std::runtime_error("mass probe nonfinite");o<<std::setprecision(17)<<x;}
inline void vector(std::ostream& o,const Vector& v){o<<'[';for(int i=0;i<v.size();++i){if(i)o<<',';scalar(o,v[i]);}o<<']';}
inline void matrix(std::ostream& o,const Matrix& m){o<<'[';for(int i=0;i<m.nrow();++i){if(i)o<<',';o<<'[';for(int j=0;j<m.ncol();++j){if(j)o<<',';scalar(o,m(i,j));}o<<']';}o<<']';}
inline std::string probe(OpenSim::Model& model,State& s){
 model.realizeVelocity(s);auto& matter=model.updMatterSubsystem();Matrix m,g;matter.calcM(s,m);matter.calcG(s,g);
 std::ostringstream o;o<<"{\"kind\":\"mass_probe\",\"time_s\":";scalar(o,s.getTime());
 o<<",\"M\":";matrix(o,m);o<<",\"G\":";matrix(o,g);o<<",\"q\":";vector(o,s.getQ());o<<",\"u\":";vector(o,s.getU());o<<",\"z\":";vector(o,s.getZ());
 o<<",\"effective_body_mass_kg\":{";bool first=true;
 for(const auto& b:model.getComponentList<OpenSim::Body>()){if(!first)o<<',';first=false;o<<'"'<<b.getName()<<"\":";scalar(o,b.getMobilizedBody().getBodyMass(s));}o<<"}}";return o.str();
}
inline std::string transfer(OpenSim::Model& model,State& s,std::istream& in){
 std::string name;double dm;Vec3 station,velocity;in>>name>>dm;
 for(int k=0;k<3;++k)in>>station[k];for(int k=0;k<3;++k)in>>velocity[k];std::string extra;
 if(!in||in>>extra||name!="torso"||!std::isfinite(dm)||std::abs(dm)>.5||!station.isFinite()||!velocity.isFinite())throw std::runtime_error("invalid bounded torso mass fixture input");
 // Explicit zero path changes no scalar, cache, mass, control or state variable.
 if(dm==0)return "{\"kind\":\"mass_transfer\",\"zero\":true}";
 using Setter=void(*)(SimbodyMatterSubsystem*,State*,int,const MassProperties*);
 auto setter=reinterpret_cast<Setter>(dlsym(RTLD_DEFAULT,"ihm_simbody_set_instance_mass"));
 if(!setter)throw std::runtime_error("mass variant setter unavailable");
 State before=s;
 try {
  model.realizeDynamics(s);auto& matter=model.updMatterSubsystem();const auto& body=model.getBodySet().get(name).getMobilizedBody();
  const auto old=body.getBodyMassProperties(s);const auto id=body.getMobilizedBodyIndex();
  if(old.getMass()+dm<=0)throw std::runtime_error("nonpositive body mass");
  if(dm<0&&(velocity-body.findStationVelocityInGround(s,station)).norm()>1e-10)throw std::runtime_error("outflow jet lacks propulsion contract");
  const Vector q=s.getQ(),u=s.getU(),z=s.getZ();const double time=s.getTime();
  const double k0=model.calcKineticEnergy(s),pe0=model.calcPotentialEnergy(s);const auto p0=matter.calcSystemMomentumAboutGroundOrigin(s);
  const Vec3 point=body.findStationLocationInGround(s,station),relative=body.getBodyRotation(s)*station;
  Vector old_momentum;matter.multiplyByM(s,u,old_momentum);
  Vector_<SpatialVec> body_impulse(matter.getNumBodies(),SpatialVec(Vec3(0),Vec3(0)));
  body_impulse[id]=SpatialVec(relative%(dm*velocity),dm*velocity);Vector input_impulse;
  matter.multiplyBySystemJacobianTranspose(s,body_impulse,input_impulse);
  const Inertia point_inertia(station,std::abs(dm));const auto inertia=dm>0?old.getInertia()+point_inertia:old.getInertia()-point_inertia;
  const MassProperties next(old.getMass()+dm,(old.getMass()*old.getMassCenter()+dm*station)/(old.getMass()+dm),inertia);
  setter(&matter,&s,int(id),&next);model.realizeDynamics(s);
  Vector unconstrained;matter.multiplyByMInv(s,old_momentum+input_impulse,unconstrained);
  Vector violation;matter.multiplyByG(s,unconstrained-u,violation);
  Vector multiplier,correction_impulse,du;matter.solveForConstraintImpulses(s,-violation,multiplier);
  matter.multiplyByGTranspose(s,multiplier,correction_impulse);matter.multiplyByMInv(s,correction_impulse,du);
  s.updU()=unconstrained+du;model.realizeDynamics(s);
  Vector momentum1;matter.multiplyByM(s,s.getU(),momentum1);
  const double generalized_residual=(momentum1-old_momentum-input_impulse-correction_impulse).norm();
  double constraint_work=0;for(int i=0;i<u.size();++i)constraint_work+=.5*(u[i]+s.getU()[i])*correction_impulse[i];
  const double k1=model.calcKineticEnergy(s),pe1=model.calcPotentialEnergy(s),kinetic=.5*dm*velocity.normSqr(),potential=-dm*dot(model.getGravity(),point);
  const auto residual=matter.calcSystemMomentumAboutGroundOrigin(s)-p0-SpatialVec(point%(dm*velocity),dm*velocity);
  const double constraint_error=s.getUErr().norm();
  if(generalized_residual>1e-9||constraint_error>1e-9||residual[0].norm()>1e-9||residual[1].norm()>1e-9||std::abs(pe1-pe0-potential)>1e-9)throw std::runtime_error("mass flux/constraint closure failure");
  if((s.getQ()-q).norm()!=0||(s.getZ()-z).norm()!=0||s.getTime()!=time)throw std::runtime_error("mass transfer reset continuous muscle/configuration state");
  std::ostringstream o;o<<"{\"kind\":\"mass_transfer\",\"zero\":false,\"delta_mass_kg\":";scalar(o,dm);
  o<<",\"kinetic_in_j\":";scalar(o,kinetic);o<<",\"potential_in_j\":";scalar(o,potential);
  o<<",\"constraint_impulse_work_j\":";scalar(o,constraint_work);o<<",\"capture_dissipation_j\":";scalar(o,k0+kinetic+constraint_work-k1);
  o<<",\"generalized_impulse_residual\":";scalar(o,generalized_residual);o<<",\"linear_momentum_residual_kg_m_s\":";scalar(o,residual[1].norm());o<<",\"angular_momentum_residual_kg_m2_s\":";scalar(o,residual[0].norm());
  o<<",\"constraint_velocity_error\":";scalar(o,constraint_error);o<<",\"constraint_multiplier_impulse\":";vector(o,multiplier);o<<'}';return o.str();
 }catch(...){s=before;s.invalidateAllCacheAtOrAbove(Stage::Instance);model.markControlsAsInvalid(s);throw;}
}
}
