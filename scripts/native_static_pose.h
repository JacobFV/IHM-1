#pragma once
// Isolated static residual evaluations. No assignment to the continuing State,
// no controller mutation, and no integration or accumulated-work updates.
#include <OpenSim/OpenSim.h>
#include "native_surface_foundation.h"
#include <cmath>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace ihm_static_pose {
inline void number(std::ostream& out, double value) {
    if (!std::isfinite(value)) throw std::runtime_error("nonfinite static pose response");
    out << std::setprecision(17) << value;
}
inline void quoted(std::ostream& out, const std::string& value) {
    out << '"';
    for (char c : value) {
        if (c == '"' || c == '\\') out << '\\';
        if (c == '\n') out << "\\n"; else out << c;
    }
    out << '"';
}
template<class Vector> inline void array(std::ostream& out, const Vector& values) {
    out << '[';
    for (int i=0; i<values.size(); ++i) { if (i) out << ','; number(out, values[i]); }
    out << ']';
}
inline std::string evaluate(OpenSim::Model& model, const SimTK::State& continuing,
                            std::istream& input, const std::string& environment,
                            bool has_external_loads, double support_plane, const ihm_surface::Foundation* foundation=nullptr) {
    if ((environment != "supine" && environment != "upright" && environment != "free") || has_external_loads)
        throw std::runtime_error("static pose requires a known environment with no external loads");
    int count;
    if (!(input >> count) || count < 1 || count > model.getCoordinateSet().getSize())
        throw std::runtime_error("invalid static pose coordinate count");
    std::map<std::string, double> requested;
    for (int i=0; i<count; ++i) {
        std::string name; double value;
        if (!(input >> name >> value) || !std::isfinite(value) || requested.count(name))
            throw std::runtime_error("invalid or duplicate static pose coordinate");
        const auto& coordinate = model.getCoordinateSet().get(name);
        if (coordinate.isDependent(continuing) || coordinate.getLocked(continuing) || coordinate.isPrescribed(continuing))
            throw std::runtime_error("static pose accepts independent free coordinates only");
        if (value < coordinate.getRangeMin() || value > coordinate.getRangeMax())
            throw std::runtime_error("static pose coordinate outside source bounds");
        requested.emplace(name,value);
    }
    // Optional suffix: activation_count muscle activation ... . The historical
    // coordinate-only request remains valid. Changes affect the copied State only.
    std::map<std::string,double> activations;
    input >> std::ws;
    if (input.peek() != std::char_traits<char>::eof()) {
        int muscle_count;
        if (!(input >> muscle_count) || muscle_count < 0 || muscle_count > model.getMuscles().getSize())
            throw std::runtime_error("invalid static activation count");
        for (int i=0;i<muscle_count;++i) {
            std::string name; double value;
            if (!(input >> name >> value) || !std::isfinite(value) || value<0 || value>1 || activations.count(name))
                throw std::runtime_error("invalid or duplicate static activation");
            model.getMuscles().get(name); // Check identity before evaluating anything.
            activations.emplace(name,value);
        }
        std::string extra;
        if (input >> extra) throw std::runtime_error("trailing static pose input");
    }
    SimTK::State candidate = continuing;
    candidate.updU() = 0.;
    for (const auto& value : requested)
        model.getCoordinateSet().get(value.first).setValue(candidate,value.second,false);
    model.assemble(candidate);
    candidate.updU() = 0.;
    model.markControlsAsInvalid(candidate);
    for(const auto& value:activations) model.getMuscles().get(value.first).setActivation(candidate,value.second);
    model.equilibrateMuscles(candidate);
    model.realizeAcceleration(candidate);
    // Assembly may alter independent coordinates. Do not let it evade bounds.
    for (const auto& coordinate : model.getComponentList<OpenSim::Coordinate>()) {
        if (!coordinate.isDependent(candidate) && !coordinate.getLocked(candidate) && !coordinate.isPrescribed(candidate)) {
            const double q = coordinate.getValue(candidate);
            if (q < coordinate.getRangeMin()-1e-10 || q > coordinate.getRangeMax()+1e-10)
                throw std::runtime_error("assembled static pose outside source bounds");
        }
    }
    const auto& system = model.getMultibodySystem();
    const auto& matter = model.getMatterSubsystem();
    SimTK::Vector tree_residual, constrained_residual;
    const SimTK::Vector zero(candidate.getNU(),0.);
    matter.calcResidualForceIgnoringConstraints(candidate,system.getMobilityForces(candidate,SimTK::Stage::Dynamics),
        system.getRigidBodyForces(candidate,SimTK::Stage::Dynamics),zero,tree_residual);
    matter.calcResidualForce(candidate,system.getMobilityForces(candidate,SimTK::Stage::Dynamics),
        system.getRigidBodyForces(candidate,SimTK::Stage::Dynamics),zero,candidate.getMultipliers(),constrained_residual);
    SimTK::Vector mass_udot;matter.multiplyByM(candidate,candidate.getUDot(),mass_udot);
    const double dynamic_identity=(constrained_residual+mass_udot).norm();
    std::ostringstream out;
    out << "{\"kind\":\"static_pose_evaluated\",\"physical_time_advanced_s\":0,\"continuing_state_unchanged\":true,\"time_s\":";
    number(out,candidate.getTime());
    out << ",\"activation_overrides\":{"; bool activation_first=true;
    for(const auto& value:activations){if(!activation_first)out<<',';activation_first=false;quoted(out,value.first);out<<":{\"requested\":";number(out,value.second);out<<",\"actual\":";number(out,model.getMuscles().get(value.first).getActivation(candidate));out<<'}';}
    out << "},\"q\":"; array(out,candidate.getQ());
    out << ",\"u\":"; array(out,candidate.getU());
    out << ",\"udot\":"; array(out,candidate.getUDot());
    std::vector<std::string> mobility_names(candidate.getNU());
    std::vector<int> mobility_rotational(candidate.getNU(),-1);
    for(const auto& coordinate:model.getComponentList<OpenSim::Coordinate>()){
        const auto& body=matter.getMobilizedBody(coordinate.getBodyIndex());
        if(body.getNumQ(candidate)!=body.getNumU(candidate))throw std::runtime_error("static force units need source coordinate/mobility mapping");
        const int index=int(body.getFirstUIndex(candidate))+int(coordinate.getMobilizerQIndex());
        if(index<0||index>=candidate.getNU()||!mobility_names[index].empty())throw std::runtime_error("ambiguous static mobility coordinate map");
        mobility_names[index]=coordinate.getName();mobility_rotational[index]=coordinate.getMotionType()==OpenSim::Coordinate::Rotational?1:0;
    }
    out<<",\"mobility_coordinate_names\":[";
    for(int i=0;i<candidate.getNU();i++){if(i)out<<',';if(mobility_names[i].empty())throw std::runtime_error("unmapped static mobility");quoted(out,mobility_names[i]);}
    out<<"],\"mobility_rotational\":";array(out,mobility_rotational);
    out << ",\"constraint_position_error\":"; number(out,candidate.getQErr().norm());
    out << ",\"constraint_velocity_error\":"; number(out,candidate.getUErr().norm());
    out << ",\"constraint_acceleration_error\":"; number(out,candidate.getUDotErr().norm());
    out<<",\"mass_times_udot\":";array(out,mass_udot);
    out<<",\"dynamic_residual_identity_norm\":";number(out,dynamic_identity);
    out << ",\"constraint_multipliers\":"; array(out,candidate.getMultipliers());
    out << ",\"tree_zero_acceleration_residual_mobility_force\":"; array(out,tree_residual);
    out << ",\"constrained_zero_acceleration_residual_mobility_force\":"; array(out,constrained_residual);
    out << ",\"residual_basis\":\"Simbody mobility order and mixed force/torque units; tree residual omits constraints; constrained residual uses native dynamic multipliers at this zero-speed pose\"";
    out << ",\"coordinates\":{";
    bool first=true;
    for (const auto& coordinate : model.getComponentList<OpenSim::Coordinate>()) {
        if (!first) out << ','; first=false;
        quoted(out,coordinate.getName()); out << ":{\"value\":"; number(out,coordinate.getValue(candidate));
        out << ",\"speed\":"; number(out,coordinate.getSpeedValue(candidate));
        out << ",\"acceleration\":"; number(out,coordinate.getAccelerationValue(candidate));
        out << ",\"rotational\":" << (coordinate.getMotionType()==OpenSim::Coordinate::Rotational?"true":"false");
        out << ",\"independent\":" << (!coordinate.isDependent(candidate)&&!coordinate.getLocked(candidate)&&!coordinate.isPrescribed(candidate)?"true":"false") << '}';
    }
    out << "},\"state_variables\":{";
    const auto names=model.getStateVariableNames();
    const auto values=model.getStateVariableValues(candidate);
    for (int i=0; i<names.getSize(); ++i) {
        if (i) out << ','; quoted(out,names[i]); out << ':'; number(out,values[i]);
    }
    out << "},\"mass_kg\":"; number(out,model.getTotalMass(candidate));
    out << ",\"gravity_m_s2\":"; array(out,model.getGravity());
    out << ",\"com_acceleration_m_s2\":"; array(out,model.calcMassCenterAcceleration(candidate));
    out << ",\"kinetic_energy_j\":"; number(out,model.calcKineticEnergy(candidate));
    out << ",\"contacts\":["; first=true;
    SimTK::Vec3 support(0); double penetration=0;
    for (const auto& force : model.getComponentList<OpenSim::SmoothSphereHalfSpaceForce>()) {
        const auto& sphere=force.getConnectee<OpenSim::ContactSphere>("sphere");
        const auto center=sphere.getFrame().findStationLocationInGround(candidate,sphere.get_location());
        const auto wrench=force.getSphereForce(candidate); support+=wrench[1];
        const double overlap=std::max(0.,support_plane+sphere.getRadius()-center[environment=="upright"?1:0]);
        penetration=std::max(penetration,overlap);
        if (!first) out << ','; first=false;
        out << "{\"name\":"; quoted(out,force.getName());
        out << ",\"force_n\":"; array(out,wrench[1]);
        out << ",\"body_frame\":"; quoted(out,sphere.getFrame().getName());
        out << ",\"center_m\":"; array(out,center);
        out << ",\"radius_m\":"; number(out,sphere.getRadius());
        out << ",\"signed_clearance_m\":"; number(out,center[environment=="upright"?1:0]-support_plane-sphere.getRadius());
        out << ",\"penetration_m\":"; number(out,overlap); out << '}';
    }
    if(foundation){const auto result=foundation->sample(candidate);support+=result.force;penetration=std::max(penetration,result.maximum_penetration);
        for(const auto& item:result.bodies){if(!first)out<<',';first=false;
            out<<"{\"name\":";quoted(out,"surface_"+item.first);out<<",\"force_n\":";array(out,item.second.force);
            out<<",\"penetration_m\":";number(out,item.second.maximum_penetration);out<<'}';}}
    out << "],\"contact_force_n\":"; array(out,support);
    out << ",\"maximum_penetration_m\":"; number(out,penetration);
    if(foundation){const auto result=foundation->sample(candidate);
        out<<",\"surface_foundation\":{\"elastic_energy_j\":";number(out,result.energy);
        out<<",\"skin_elastic_energy_j\":";number(out,result.skin_energy);out<<",\"bed_elastic_energy_j\":";number(out,result.bed_energy);
        out<<",\"maximum_bed_indentation_m\":";number(out,result.maximum_bed_indentation);out<<",\"maximum_total_approach_m\":";number(out,result.maximum_total_approach);
        out<<",\"bed_force_n\":";array(out,-result.force);out<<",\"bed_moment_about_source_origin_nm\":";array(out,result.bed_moment);out<<'}';}
    out << ",\"accepted_equilibrium\":false,\"scope\":\"Static candidate only; requires sustained forward verification and explicit new reference initialization\"}";
    return out.str();
}
} // namespace ihm_static_pose
