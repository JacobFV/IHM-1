#pragma once
// Explicit local material boundary. Opt-in diagnostic State-mass variant only.
#include <OpenSim/OpenSim.h>
#include <dlfcn.h>
#include <map>
#include <sstream>
#include <iomanip>
#include <cctype>
#include <cstdint>
#include <cstdlib>
namespace ihm_mass_port {
using namespace SimTK;
inline void number(std::ostream& o,double x){if(!std::isfinite(x))throw std::runtime_error("nonfinite mass port ledger");o<<std::setprecision(17)<<x;}
inline void vector(std::ostream& o,const Vec3& v){o<<'[';for(int k=0;k<3;++k){if(k)o<<',';number(o,v[k]);}o<<']';}
inline bool token(const std::string& x){if(x.empty()||x.size()>128)return false;for(unsigned char c:x)if(!std::isalnum(c)&&c!='_'&&c!='-'&&c!='.'&&c!=':')return false;return true;}
struct Inventory {std::string body;Vec3 station;double mass_kg=0;};
struct Ledger {
 bool enabled=false;std::string reference_id;std::uint64_t last_sequence=0;
 std::map<std::string,Inventory> owners;std::string last_receipt="null";
 double kinetic_in_j=0,potential_in_j=0,capture_dissipation_j=0,constraint_work_j=0;
 explicit Ledger(bool on=false,const std::string& reference=""):enabled(on),reference_id(reference){
  if(enabled&&(!token(reference_id)||!dlsym(RTLD_DEFAULT,"ihm_simbody_set_instance_mass")))throw std::runtime_error("instance mass mode requires bound native variant/reference");
 }
 std::string json() const {
  std::ostringstream o;o<<"{\"enabled\":"<<(enabled?"true":"false")<<",\"reference_id\":\""<<reference_id<<"\",\"last_sequence\":"<<last_sequence<<",\"owners\":{";bool first=true;double total=0;
  for(const auto& [name,item]:owners){if(!first)o<<',';first=false;total+=item.mass_kg;o<<'"'<<name<<"\":{\"body\":\""<<item.body<<"\",\"mass_kg\":";number(o,item.mass_kg);o<<",\"station_m\":[";for(int k=0;k<3;++k){if(k)o<<',';number(o,item.station[k]);}o<<"]}";}
  o<<"},\"cumulative_kinetic_in_j\":";number(o,kinetic_in_j);o<<",\"cumulative_potential_in_j\":";number(o,potential_in_j);o<<",\"cumulative_capture_dissipation_j\":";number(o,capture_dissipation_j);o<<",\"cumulative_constraint_impulse_work_j\":";number(o,constraint_work_j);o<<",\"owned_payload_mass_kg\":";number(o,total);o<<",\"last_receipt\":"<<last_receipt<<",\"scope\":\"Explicit added local payload only; baseline GI/tissue mass remains excluded from this inventory; no physiology mass inference\"}";return o.str();
 }
 void apply(OpenSim::Model& model,State& s,std::istream& in) {
  if(!enabled)throw std::runtime_error("local mass port disabled");
  std::string reference,owner,body_name,sequence_text;double dm;Vec3 station,velocity;
  in>>reference>>sequence_text>>owner>>body_name>>dm;for(int k=0;k<3;++k)in>>station[k];for(int k=0;k<3;++k)in>>velocity[k];std::string extra;
  if(!in||in>>extra||reference!=reference_id||!token(owner)||body_name!="torso"||!std::isfinite(dm)||std::abs(dm)>.5||!station.isFinite()||!velocity.isFinite())throw std::runtime_error("invalid bound local mass command");
  if(sequence_text.empty()||sequence_text.size()>20||sequence_text.find_first_not_of("0123456789")!=std::string::npos)throw std::runtime_error("invalid mass sequence");
  const auto sequence=std::stoull(sequence_text);if(last_sequence==UINT64_MAX||sequence!=last_sequence+1)throw std::runtime_error("duplicate/out-of-order mass sequence");
  const auto found=owners.find(owner);const double inventory=found==owners.end()?0:found->second.mass_kg;
  if(found!=owners.end()&&(found->second.body!=body_name||found->second.station!=station))throw std::runtime_error("mass owner/site binding changed");
  if(found==owners.end()&&dm!=0&&owners.size()>=64)throw std::runtime_error("mass owner capacity exceeded");
  double total_inventory=0;for(const auto& item:owners)total_inventory+=item.second.mass_kg;
  if(total_inventory+dm>.5)throw std::runtime_error("total explicit payload exceeds validated 0.5 kg domain; splitting unsupported");
  if(inventory+dm<0)throw std::runtime_error("mass removal exceeds explicit owned inventory");
  if(dm==0){last_sequence=sequence;last_receipt="{\"zero\":true,\"sequence\":"+std::to_string(sequence)+"}";return;}
  const State before=s;const Ledger prior=*this;
  try {
   model.realizeDynamics(s);auto& matter=model.updMatterSubsystem();const auto& body=model.getBodySet().get(body_name).getMobilizedBody();const auto id=body.getMobilizedBodyIndex();
   const auto old=body.getBodyMassProperties(s);const Vector q=s.getQ(),u=s.getU(),z=s.getZ();const double time=s.getTime();
   if(dm<0&&(velocity-body.findStationVelocityInGround(s,station)).norm()>1e-10)throw std::runtime_error("outflow requires co-moving material; jet propulsion unsupported");
   const double k0=model.calcKineticEnergy(s),pe0=model.calcPotentialEnergy(s);const auto p0=matter.calcSystemMomentumAboutGroundOrigin(s);
   const Vec3 point=body.findStationLocationInGround(s,station),relative=body.getBodyRotation(s)*station;
   Vector old_momentum;matter.multiplyByM(s,u,old_momentum);Vector_<SpatialVec> spatial(matter.getNumBodies(),SpatialVec(Vec3(0),Vec3(0)));
   spatial[id]=SpatialVec(relative%(dm*velocity),dm*velocity);Vector input;matter.multiplyBySystemJacobianTranspose(s,spatial,input);
   const Inertia point_inertia(station,std::abs(dm));const auto inertia=dm>0?old.getInertia()+point_inertia:old.getInertia()-point_inertia;
   const MassProperties next(old.getMass()+dm,(old.getMass()*old.getMassCenter()+dm*station)/(old.getMass()+dm),inertia);
   using Setter=void(*)(SimbodyMatterSubsystem*,State*,int,const MassProperties*);
   auto setter=reinterpret_cast<Setter>(dlsym(RTLD_DEFAULT,"ihm_simbody_set_instance_mass"));setter(&matter,&s,int(id),&next);
#ifdef IHM_MASS_TEST_AFTER_SET
   IHM_MASS_TEST_AFTER_SET;
#endif
   model.realizeDynamics(s);Vector raw;matter.multiplyByMInv(s,old_momentum+input,raw);Vector violation;
   matter.multiplyByG(s,raw-u,violation);Vector lambda,correction,du;matter.solveForConstraintImpulses(s,-violation,lambda);matter.multiplyByGTranspose(s,lambda,correction);matter.multiplyByMInv(s,correction,du);
   s.updU()=raw+du;model.realizeDynamics(s);Vector momentum;matter.multiplyByM(s,s.getU(),momentum);
   const double generalized_error=(momentum-old_momentum-input-correction).norm();double constraint_work=0;
   for(int i=0;i<u.size();++i)constraint_work+=.5*(u[i]+s.getU()[i])*correction[i];
   const double k1=model.calcKineticEnergy(s),pe1=model.calcPotentialEnergy(s),kinetic=.5*dm*velocity.normSqr(),potential=-dm*dot(model.getGravity(),point),dissipation=k0+kinetic+constraint_work-k1;
   const auto residual=matter.calcSystemMomentumAboutGroundOrigin(s)-p0-SpatialVec(point%(dm*velocity),dm*velocity);
   for(double v:{generalized_error,constraint_work,k1,pe1,kinetic,potential,dissipation,residual.norm(),s.getUErr().norm()})if(!std::isfinite(v))throw std::runtime_error("nonfinite mass transfer result");
   if(generalized_error>1e-9||s.getUErr().norm()>1e-9||residual[0].norm()>1e-9||residual[1].norm()>1e-9||std::abs(pe1-pe0-potential)>1e-9||dissipation < -1e-9||(dm<0&&std::abs(dissipation)>1e-9))throw std::runtime_error("mass impulse/energy outside validated domain");
   if((s.getQ()-q).norm()!=0||(s.getZ()-z).norm()!=0||s.getTime()!=time)throw std::runtime_error("mass transfer changed preserved continuous state");
   owners[owner]={body_name,station,inventory+dm};last_sequence=sequence;
   kinetic_in_j+=kinetic;potential_in_j+=potential;capture_dissipation_j+=dissipation;constraint_work_j+=constraint_work;
   std::ostringstream o;o<<"{\"zero\":false,\"sequence\":"<<sequence<<",\"owner\":\""<<owner<<"\",\"body\":\""<<body_name<<"\",\"delta_mass_kg\":";number(o,dm);
   o<<",\"time_s\":";number(o,time);o<<",\"station_m\":";vector(o,station);o<<",\"point_source_m\":";vector(o,point);o<<",\"velocity_source_m_s\":";vector(o,velocity);o<<",\"body_mass_before_kg\":";number(o,old.getMass());o<<",\"body_mass_after_kg\":";number(o,next.getMass());
   o<<",\"kinetic_in_j\":";number(o,kinetic);o<<",\"potential_in_j\":";number(o,potential);o<<",\"capture_dissipation_j\":";number(o,dissipation);o<<",\"constraint_impulse_work_j\":";number(o,constraint_work);
   o<<",\"generalized_impulse_residual\":";number(o,generalized_error);o<<",\"linear_momentum_residual_kg_m_s\":";number(o,residual[1].norm());o<<",\"angular_momentum_residual_kg_m2_s\":";number(o,residual[0].norm());o<<",\"constraint_velocity_error\":";number(o,s.getUErr().norm());o<<'}';last_receipt=o.str();
  }catch(...){s=before;*this=prior;s.invalidateAllCacheAtOrAbove(Stage::Instance);model.markControlsAsInvalid(s);throw;}
 }
};
}
