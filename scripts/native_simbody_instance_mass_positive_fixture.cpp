// Positive acceptance for the separately linked State-instance mass variant.
#include "native_simbody_instance_mass.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
using namespace SimTK;
static void check(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
static bool near(double a,double b){return std::abs(a-b)<1e-11*std::max({1.,std::abs(a),std::abs(b)});}
struct Flux {double kinetic_in,potential_in,dissipation,momentum_error;};
static Flux transfer(MultibodySystem& sys,SimbodyMatterSubsystem& matter,State& s,
                     MobilizedBodyIndex id,double dm,const Vec3& station,const Vec3& velocity){
 check(matter.getNumConstraints()==0,"prototype needs unconstrained bodies");
 for(int i=1;i<matter.getNumBodies();++i)
  check(MobilizedBody::Free::isInstanceOf(matter.getMobilizedBody(MobilizedBodyIndex(i))),"prototype needs independent free bodies");
 State before=s;
 try {
  sys.realize(s,Stage::Velocity);const auto& body=matter.getMobilizedBody(id);
  const auto old=body.getBodyMassProperties(s);check(std::isfinite(dm)&&old.getMass()+dm>0&&station.isFinite()&&velocity.isFinite(),"invalid flux");
  if(dm<0)check((velocity-body.findStationVelocityInGround(s,station)).norm()<1e-12,"jet requires propulsion work");
  const auto p0=matter.calcSystemMomentumAboutGroundOrigin(s);const double k0=matter.calcKineticEnergy(s),v0=sys.calcPotentialEnergy(s);
  const Vec3 x=body.findStationLocationInGround(s,station),relative=body.getBodyRotation(s)*station;
  Vector momentum;matter.multiplyByM(s,s.getU(),momentum);
  Vector_<SpatialVec> impulse(matter.getNumBodies(),SpatialVec(Vec3(0),Vec3(0)));
  impulse[id]=SpatialVec(relative%(dm*velocity),dm*velocity);Vector incoming;
  matter.multiplyBySystemJacobianTranspose(s,impulse,incoming);
  const Inertia point(station,std::abs(dm));const auto inertia=dm>=0?old.getInertia()+point:old.getInertia()-point;
  const MassProperties next(old.getMass()+dm,(old.getMass()*old.getMassCenter()+dm*station)/(old.getMass()+dm),inertia);
  ihm_simbody_set_instance_mass(&matter,&s,int(id),&next);sys.realize(s,Stage::Position);
  Vector speed;matter.multiplyByMInv(s,momentum+incoming,speed);s.updU()=speed;sys.realize(s,Stage::Velocity);
  const auto p1=matter.calcSystemMomentumAboutGroundOrigin(s);const double k1=matter.calcKineticEnergy(s),v1=sys.calcPotentialEnergy(s);
  const double momentum_error=(p1-p0-SpatialVec(x%(dm*velocity),dm*velocity)).norm();check(momentum_error<1e-10,"momentum flux mismatch");
  const double kinetic=.5*dm*velocity.normSqr(),potential=dm*9.80665*x[1];check(near(v1-v0,potential),"gravity transport mismatch");
  return {kinetic,potential,k0+kinetic-k1,momentum_error};
 }catch(...){s=before;s.invalidateAllCacheAtOrAbove(Stage::Instance);throw;}
}
int main(){try{
 MultibodySystem sys;SimbodyMatterSubsystem matter(sys);GeneralForceSubsystem forces(sys);Force::Gravity gravity(forces,matter,Vec3(0,-9.80665,0));
 Body::Rigid carrier(MassProperties(10.,Vec3(0),Inertia(2.,3.,4.)));
 MobilizedBody::Free torso(matter.Ground(),Transform(),carrier,Transform());
 MobilizedBody::Free other(matter.Ground(),Transform(Vec3(2,0,0)),carrier,Transform());
 sys.realizeTopology();State s=sys.getDefaultState();s.setTime(3.25);sys.realize(s,Stage::Model);
 Vector u(s.getNU(),0.);u[0]=.2;u[1]=.3;u[2]=.1;u[3]=1.;u[4]=.5;s.updU()=u;sys.realize(s,Stage::Velocity);
 const State original=s;const Vector q0=s.getQ();const Vec3 station(.2,.3,-.1);
 const double k0=matter.calcKineticEnergy(s),v0=sys.calcPotentialEnergy(s);Matrix m0;matter.calcM(s,m0);
 const Vec3 material=torso.findStationVelocityInGround(s,station);const double expected_k=.25*material.normSqr();
 const auto old=torso.getBodyMassProperties(s);const MassProperties added(10.5,.5*station/10.5,old.getInertia()+Inertia(station,.5));
 ihm_simbody_set_instance_mass(&matter,&s,int(torso.getMobilizedBodyIndex()),&added);sys.realize(s,Stage::Velocity);Matrix m1;matter.calcM(s,m1);
 check(torso.getBodyMass(s)==10.5&&matter.calcSystemMass(s)==20.5&&other.getBodyMass(s)==10.,"mass getter/locality mismatch");
 check((m1-m0).norm()>0,"mass matrix unchanged");check(near(matter.calcKineticEnergy(s)-k0,expected_k),"kinetic delta mismatch");
 check(near(sys.calcPotentialEnergy(s)-v0,.5*9.80665*.3),"potential delta mismatch");
 const Vec3 expected_com=(10.5*torso.findMassCenterLocationInGround(s)+10.*other.findMassCenterLocationInGround(s))/20.5;
 check((matter.calcSystemMassCenterLocationInGround(s)-expected_com).norm()<1e-12,"COM owner mismatch");
 sys.realize(s,Stage::Acceleration);const auto acceleration=torso.getBodyAcceleration(s);
 check(std::isfinite(acceleration.norm()),"nonfinite changed-mass acceleration");
 check((matter.calcSystemMassCenterAccelerationInGround(s)-Vec3(0,-9.80665,0)).norm()<1e-10,"gravity/gyroscopic COM acceleration mismatch");
 check(s.getTime()==3.25&&(s.getQ()-q0).norm()==0&&(s.getU()-u).norm()==0,"setter changed physical state");
 // Concurrent original State remains unmodified by the instance setter.
 State independent=original;independent.invalidateAllCacheAtOrAbove(Stage::Instance);sys.realize(independent,Stage::Velocity);
 check(torso.getBodyMass(independent)==10.&&near(matter.calcKineticEnergy(independent),k0),"mass leaked between States");
 s=independent;const auto capture=transfer(sys,matter,s,torso.getMobilizedBodyIndex(),.5,station,Vec3(.1,-.2,.3));
 check(capture.dissipation>=0,"capture generated kinetic energy");const State captured=s;const Vector captured_u=s.getU();
 s=independent;s.invalidateAllCacheAtOrAbove(Stage::Instance);const auto replay=transfer(sys,matter,s,torso.getMobilizedBodyIndex(),.5,station,Vec3(.1,-.2,.3));
 check((s.getU()-captured_u).norm()==0&&replay.dissipation==capture.dissipation,"capture replay failed");
 const Vec3 departing=torso.findStationVelocityInGround(s,station);const auto release=transfer(sys,matter,s,torso.getMobilizedBodyIndex(),-.5,station,departing);
 check(near(release.dissipation,0.)&&(s.getU()-captured_u).norm()<1e-12&&torso.getBodyMass(s)==10.,"passive outflow mismatch");
 bool rejected=false;try{const MassProperties bad(-1.,Vec3(0),Inertia(1.));ihm_simbody_set_instance_mass(&matter,&s,int(torso.getMobilizedBodyIndex()),&bad);}catch(const std::exception&){rejected=true;}
 check(rejected&&torso.getBodyMass(s)==10.,"invalid mass accepted/mutated");
 std::cout<<std::setprecision(17)<<"{\"passed\":true,\"variable_mass_supported\":true,\"comoving_kinetic_in_j\":"<<expected_k<<",\"capture_dissipation_j\":"<<capture.dissipation<<",\"capture_momentum_error\":"<<capture.momentum_error<<",\"release_dissipation_j\":"<<release.dissipation<<",\"release_momentum_error\":"<<release.momentum_error<<",\"mass_matrix_change_norm\":"<<(m1-m0).norm()<<"}\n";
 return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
