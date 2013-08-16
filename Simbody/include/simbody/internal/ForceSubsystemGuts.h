#ifndef SimTK_SIMBODY_FORCE_SUBSYSTEM_GUTS_H
#define SimTK_SIMBODY_FORCE_SUBSYSTEM_GUTS_H

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

/** @file
 * Define the extendable library-side implementation of the ForceSubsystem.
 */

#include "SimTKcommon.h"

#include "simbody/internal/common.h"
#include "simbody/internal/MultibodySystem.h"

namespace SimTK {

/// Public declaration of internals for ForceSubsystem extension
class ForceSubsystem::Guts : public Subsystem::Guts {
public:
    Guts(const String& name, const String& version) 
      : Subsystem::Guts(name,version)
    {
    }

    // Make sure the virtual destructor in Subsystem::Guts remains
    // virtual in this intermediate class.
    virtual ~Guts() { }

    // All the other Subsystem::Guts virtuals remain unresolved.

    // Return the MultibodySystem which owns this ForceSubsystem.
    const MultibodySystem& getMultibodySystem() const {
        return MultibodySystem::downcast(getSystem());
    }
    
    /// Get this subsystem's contribution to the potential energy.  The state must
    /// be at Dynamics stage or later.
    virtual Real calcPotentialEnergy(const State& state) const = 0;

    SimTK_DOWNCAST(ForceSubsystem::Guts, Subsystem::Guts);
};

// typedef ForceSubsystem::Guts ForceSubsystemRep;

} // namespace SimTK

#endif // SimTK_SIMBODY_FORCE_SUBSYSTEM_GUTS_H
