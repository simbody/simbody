#ifndef SimTK_SIMBODY_MULTIBODY_SYSTEM_H_
#define SimTK_SIMBODY_MULTIBODY_SYSTEM_H_

/* Portions copyright (c) 2005-7 Stanford University and Michael Sherman.
 * Contributors:
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/System.h"
#include "simbody/internal/MatterSubsystem.h"

#include <vector>

namespace SimTK {

class AnalyticGeometry;
class DecorativeGeometry;

class MatterSubsystem;
class ForceSubsystem;
class DecorationSubsystem;


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

    // We inherit realize() from System, and add constraint projection here.
    // We are given a state whose continuous state variables y may violate
    // a set of constraints at position (q) and velocity (u) levels. In addition
    // we may be given a set of absolute error estimates y_err for y. This solver
    // performs two operations:
    //   (1) perform a least squares projection of y onto the constraint manifold,
    //       using the error test norm to define the least-squares direction
    //   (2) perform the same projection on y_err, returning a revised y_err which
    //       has a smaller norm
    // This routine returns true if any change was made to s or y_err, otherwise false.
    // 

    bool project(State& s, Vector& y_err, 
                 const Real& tol,               // must achieve this tolerance or better
                 const Real& dontProjectFac,    // skip projection if tol <= fac*tol
                 const Real& targetTol          // when projecting, try for this (<= tol)
                 ) const;


    // Steals ownership of the source; returns subsystem ID number.
    int setMatterSubsystem(SimbodyMatterSubsystem&);
    int addForceSubsystem(ForceSubsystem&);
    int setDecorationSubsystem(DecorationSubsystem&);
    const SimbodyMatterSubsystem& getMatterSubsystem() const;
    SimbodyMatterSubsystem&       updMatterSubsystem();
    const DecorationSubsystem& getDecorationSubsystem() const;
    DecorationSubsystem&       updDecorationSubsystem();

    // Responses available when the global subsystem is advanced to the indicated stage.
    const Real& getPotentialEnergy(const State&, Stage g=Stage::Dynamics) const;
    const Real& getKineticEnergy(const State&, Stage g=Stage::Dynamics) const;

    Real getEnergy(const State& s, Stage g=Stage::Dynamics) const {
        return getPotentialEnergy(s,g)+getKineticEnergy(s,g);
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
    Real&                updPotentialEnergy(const State&, Stage) const;
    Vector_<SpatialVec>& updRigidBodyForces(const State&, Stage) const;
    Vector_<Vec3>&       updParticleForces (const State&, Stage) const;
    Vector&              updMobilityForces (const State&, Stage) const;

    // This is for use by the matter subsystem while realizing Dynamics stage.
    Real& updKineticEnergy(const State&, Stage) const;

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
