#pragma once
#include <Simbody.h>
// Separate diagnostic SimTKsimbody variant only. No layout change, no topology
// mutation. Rejects bodies with degenerate topology inertia in the first ABI.
extern "C" void ihm_simbody_set_instance_mass(
    SimTK::SimbodyMatterSubsystem*,SimTK::State*,int,const SimTK::MassProperties*);
