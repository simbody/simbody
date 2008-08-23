#ifndef SimTK_SIMBODY_MULTIBODY_SYSTEM_H_
#define SimTK_SIMBODY_MULTIBODY_SYSTEM_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-8 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"

#include <vector>

namespace SimTK {

class SimbodyMatterSubsystem;
class ForceSubsystem;
class DecorationSubsystem;
class GeneralContactSubsystem;


/**
 * The job of the MultibodySystem class is to coordinate the activities of various
 * subsystems which can be part of a multibody system. We insist on having exactly one
 * MatterSubsystem, and we would like also to have:
 *    - one or more ForceSubsystems
 *    - an AnalyticGeometrySubsystem
 *    - a MassPropertiesSubsystem
 *    - a VisualizationSubsystem
 * There will also be a generic System-level "subsystem" for global variables.
 */
class SimTK_SIMBODY_EXPORT MultibodySystem : public System {
public:
    MultibodySystem();
    MultibodySystem(SimbodyMatterSubsystem& m);

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

    // These methods are for use by our constituent subsystems to communicate with
    // each other and with the MultibodySystem as a whole.

    // These cache entries belong to the global subsystem, which zeroes them at the
    // start of the corresponding stage. They are filled in by the force subsystems when
    // they are realized to each stage. Forces are cumulative from stage to stage,
    // so the Dynamics stage includes everything. That may then be accessed by the matter 
    // subsystem in Acceleration stage to generate the accelerations.
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

class SimTK_SIMBODY_EXPORT MultibodyDynamicsStudy : public Study {
public:
    MultibodyDynamicsStudy() { }
    MultibodyDynamicsStudy(const MultibodyDynamicsStudy&);
    MultibodyDynamicsStudy& operator=(const MultibodyDynamicsStudy&);
    ~MultibodyDynamicsStudy();

    MultibodyDynamicsStudy(const MultibodySystem&);

    const MultibodySystem& getMultibodySystem() const;

    void advanceTimeBy(const Real& h); //TODO

    SimTK_PIMPL_DOWNCAST(MultibodyDynamicsStudy, Study);
};


} // namespace SimTK

#endif // SimTK_SIMBODY_MULTIBODY_SYSTEM_H_
