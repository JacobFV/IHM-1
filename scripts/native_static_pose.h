#pragma once
// Isolated static residual evaluations. No assignment to the continuing State,
// no controller mutation, and no integration or accumulated-work updates.
#include <OpenSim/OpenSim.h>
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
                            bool has_external_loads, double support_plane) {
    if (environment != "supine" || has_external_loads)
        throw std::runtime_error("static pose requires supine with no external loads");
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
    std::string extra;
    if (input >> extra) throw std::runtime_error("trailing static pose input");
    SimTK::State candidate = continuing;
    candidate.updU() = 0.;
    for (const auto& value : requested)
        model.getCoordinateSet().get(value.first).setValue(candidate,value.second,false);
    model.assemble(candidate);
    candidate.updU() = 0.;
    model.markControlsAsInvalid(candidate);
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
    std::ostringstream out;
    out << "{\"kind\":\"static_pose_evaluated\",\"physical_time_advanced_s\":0,\"continuing_state_unchanged\":true,\"time_s\":";
    number(out,candidate.getTime());
    out << ",\"q\":"; array(out,candidate.getQ());
    out << ",\"u\":"; array(out,candidate.getU());
    out << ",\"udot\":"; array(out,candidate.getUDot());
    out << ",\"constraint_position_error\":"; number(out,candidate.getQErr().norm());
    out << ",\"constraint_velocity_error\":"; number(out,candidate.getUErr().norm());
    out << ",\"constraint_acceleration_error\":"; number(out,candidate.getUDotErr().norm());
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
        const double overlap=std::max(0.,support_plane+sphere.getRadius()-center[0]);
        penetration=std::max(penetration,overlap);
        if (!first) out << ','; first=false;
        out << "{\"name\":"; quoted(out,force.getName());
        out << ",\"force_n\":"; array(out,wrench[1]);
        out << ",\"penetration_m\":"; number(out,overlap); out << '}';
    }
    out << "],\"contact_force_n\":"; array(out,support);
    out << ",\"maximum_penetration_m\":"; number(out,penetration);
    out << ",\"accepted_equilibrium\":false,\"scope\":\"Static candidate only; requires sustained forward verification and explicit new reference initialization\"}";
    return out.str();
}
} // namespace ihm_static_pose
