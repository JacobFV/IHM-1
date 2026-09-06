// Falsification probe for private instance mass update; NOT a runtime mass port.
#include <Simbody.h>
#include "SimbodyMatterSubsystemRep.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
using namespace SimTK;
static void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
int main(){try{
 MultibodySystem system;SimbodyMatterSubsystem matter(system);
 Body::Rigid carrier(MassProperties(10.,Vec3(0),Inertia(2.,3.,4.)));
 MobilizedBody::Free torso(matter.Ground(),Transform(),carrier,Transform());
 MobilizedBody::Free other(matter.Ground(),Transform(Vec3(2,0,0)),carrier,Transform());
 system.realizeTopology();State s=system.getDefaultState();s.setTime(3.25);
 system.realize(s,Stage::Model);Vector u(s.getNU(),0.);u[0]=.2;u[1]=.3;u[2]=.1;u[3]=1.;u[4]=.5;s.updU()=u;
 system.realize(s,Stage::Velocity);State original=s;const Vector q0=s.getQ();
 Matrix mass0;matter.calcM(s,mass0);const double k0=matter.calcKineticEnergy(s);
 const auto old=torso.getBodyMassProperties(s);const Vec3 station(.2,.3,-.1);const double dm=.5;
 const Vec3 velocity=torso.findStationVelocityInGround(s,station);
 const double kinetic_flux=.5*dm*velocity.normSqr();
 const double potential_flux=dm*9.80665*torso.findStationLocationInGround(s,station)[1];
 const MassProperties next(old.getMass()+dm,(old.getMass()*old.getMassCenter()+dm*station)/(old.getMass()+dm),old.getInertia()+Inertia(station,dm));
 (void)next.calcCentralInertia();
 matter.updRep().updBodyMassProperties(s,torso.getMobilizedBodyIndex())=next;
 system.realize(s,Stage::Velocity);Matrix mass1;matter.calcM(s,mass1);
 const double k1=matter.calcKineticEnergy(s);
 const double matrix_change=(mass1-mass0).norm();
 require(matter.getRep().getBodyMassProperties(s,torso.getMobilizedBodyIndex()).getMass()==10.5&&other.getBodyMass(s)==10.,"internal instance mass not local");
 require(torso.getBodyMass(s)==10.,"public body getter changed: re-audit support");
 require(s.getTime()==3.25&&(s.getQ()-q0).norm()==0&&(s.getU()-u).norm()==0,"state time/configuration/speed changed");
 // A co-moving point payload MUST add kinetic energy without changing velocity.
 // Held implementation falsifies that contract: mass cache changes, M does not.
 require(kinetic_flux>0,"degenerate test velocity");
 require(matrix_change==0&&k1==k0,"native implementation changed: re-audit support");
 const State changed=s;
 s=original;s.invalidateAllCacheAtOrAbove(Stage::Instance);system.realize(s,Stage::Velocity);
 require(matter.getRep().getBodyMassProperties(s,torso.getMobilizedBodyIndex()).getMass()==10.&&(s.getU()-u).norm()==0,"State-copy restore failed");
 s=changed;s.invalidateAllCacheAtOrAbove(Stage::Instance);system.realize(s,Stage::Velocity);
 require(matter.getRep().getBodyMassProperties(s,torso.getMobilizedBodyIndex()).getMass()==10.5&&matter.calcKineticEnergy(s)==k1,"State-copy replay failed");
 // Inspect momentum operator's flux mismatch too; do not correct velocity using
 // a matrix demonstrated to retain the old inertia.
 original.invalidateAllCacheAtOrAbove(Stage::Instance);system.realize(original,Stage::Velocity);
 const auto p0=matter.calcSystemMomentumAboutGroundOrigin(original);
 const auto p1=matter.calcSystemMomentumAboutGroundOrigin(s);
 const Vec3 x=torso.findStationLocationInGround(s,station);
 const SpatialVec expected_flux(x%(dm*velocity),dm*velocity);
 const double momentum_flux_error=(p1-p0-expected_flux).norm();
 std::cout<<std::setprecision(17)<<"{\"passed\":true,\"variable_mass_supported\":false,\"public_body_mass_kg\":"<<torso.getBodyMass(s)<<",\"internal_instance_mass_kg\":"<<matter.getRep().getBodyMassProperties(s,torso.getMobilizedBodyIndex()).getMass()<<",\"reported_system_mass_kg\":"<<matter.calcSystemMass(s)<<",\"unchanged_other_mass_kg\":"<<other.getBodyMass(s)<<",\"mass_matrix_change_norm\":"<<matrix_change<<",\"actual_kinetic_change_j\":"<<k1-k0<<",\"required_comoving_kinetic_in_j\":"<<kinetic_flux<<",\"potential_in_j\":"<<potential_flux<<",\"momentum_flux_error\":"<<momentum_flux_error<<"}\n";
 return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
