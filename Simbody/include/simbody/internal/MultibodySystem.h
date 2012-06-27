#ifndef SimTK_SIMBODY_MULTIBODY_SYSTEM_H_
#define SimTK_SIMBODY_MULTIBODY_SYSTEM_H_

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

#include "SimTKcommon.h"
#include "simbody/internal/common.h"

#include <vector>

namespace SimTK {

class SimbodyMatterSubsystem;
class ForceSubsystem;
class DecorationSubsystem;
class GeneralContactSubsystem;


/** The job of the MultibodySystem class is to coordinate the activities of 
various subsystems which can be part of a multibody system. We insist on 
having exactly one SimbodyMatterSubsystem, and we would like also to have:
    - one or more ForceSubsystems
    - a DecorationSubsystem for visualization
    - a GeneralContactSubsystem for contact geometry
There will also be a generic System-level "subsystem" for global variables.
**/
class SimTK_SIMBODY_EXPORT MultibodySystem : public System {
public:
    MultibodySystem();
    explicit MultibodySystem(SimbodyMatterSubsystem& m);

    // Steals ownership of the source; returns subsystem ID number.
    int addForceSubsystem(ForceSubsystem&);

    int setMatterSubsystem(SimbodyMatterSubsystem&);
    const SimbodyMatterSubsystem& getMatterSubsystem() const;
    SimbodyMatterSubsystem&       updMatterSubsystem();
    bool hasMatterSubsystem() const;

    int setDecorationSubsystem(DecorationSubsystem&);
    const DecorationSubsystem& getDecorationSubsystem() const;
    DecorationSubsystem&       updDecorationSubsystem();
    bool hasDecorationSubsystem() const;

    int setContactSubsystem(GeneralContactSubsystem&);
    const GeneralContactSubsystem& getContactSubsystem() const;
    GeneralContactSubsystem&       updContactSubsystem();
    bool hasContactSubsystem() const;


    /// Calculate the total potential energy of the system.  The state must
    /// be at Dynamics stage or later.
    const Real calcPotentialEnergy(const State&) const;
    /// Calculate the total kinetic energy of the system.  The state must
    /// be at Velocity stage or later.
    const Real calcKineticEnergy(const State&) const;
    /// Calculate the total energy of the system.  The state must
    /// be at Dynamics stage or later.
    Real calcEnergy(const State& s) const {
        return calcPotentialEnergy(s)+calcKineticEnergy(s);
    }

    // These methods are for use by our constituent subsystems to communicate 
    // with each other and with the MultibodySystem as a whole.

    // These cache entries belong to the global subsystem, which zeroes them at
    // the start of the corresponding stage. They are filled in by the force 
    // subsystems when they are realized to each stage. Forces are cumulative 
    // from stage to stage, so the Dynamics stage includes everything. That may
    // then be accessed by the matter subsystem in Acceleration stage to 
    // generate the accelerations.
    const Vector_<SpatialVec>& getRigidBodyForces(const State&, Stage) const;
    const Vector_<Vec3>&       getParticleForces (const State&, Stage) const;
    const Vector&              getMobilityForces (const State&, Stage) const;

    // These routines are for use by force subsystems during Dynamics stage.
    Vector_<SpatialVec>& updRigidBodyForces(const State&, Stage) const;
    Vector_<Vec3>&       updParticleForces (const State&, Stage) const;
    Vector&              updMobilityForces (const State&, Stage) const;

    // Private implementation.
    SimTK_PIMPL_DOWNCAST(MultibodySystem, System);
    class MultibodySystemRep& updRep();
    const MultibodySystemRep& getRep() const;
protected:
    explicit MultibodySystem(MultibodySystemRep*);
};


} // namespace SimTK

#endif // SimTK_SIMBODY_MULTIBODY_SYSTEM_H_
