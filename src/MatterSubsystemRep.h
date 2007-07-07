#ifndef SimTK_MATTER_SUBSYSTEM_REP_H_
#define SimTK_MATTER_SUBSYSTEM_REP_H_

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

/** @file
 * Define the private implementation MatterSubsystemRep of a MatterSubsystem, 
 * a still-abstract class derived from abstract base class SubsystemRep.
 */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/ForceSubsystem.h"
#include "simbody/internal/AnalyticGeometry.h"
#include "simbody/internal/DecorativeGeometry.h"

#include "SubsystemRep.h"

namespace SimTK {

class State;

class MatterSubsystemRep : public SubsystemRep {
public:
    MatterSubsystemRep(const String& name, const String& version)
      : SubsystemRep(name,version)
    {
        // Ground must be the 0th MobilizedBody; concrete MatterSubsystem
        // should be sure to call createGroundBody() after handle is set.
    }
    virtual ~MatterSubsystemRep() {
        invalidateSubsystemTopologyCache();
    }

    // Return the MultibodySystem which owns this MatterSubsystem.
    const MultibodySystem& getMultibodySystem() const {
        return MultibodySystem::downcast(getSystem());
    }

    const MatterSubsystem& getMyMatterSubsystemHandle() const {
        return MatterSubsystem::downcast(getMyHandle());
    }
    MatterSubsystem& updMyMatterSubsystemHandle() {
        return MatterSubsystem::updDowncast(updMyHandle());
    }


        // TOPOLOGY STAGE //

    virtual int getNBodies()      const = 0;    // includes ground, also # mobilizers (tree joints) +1
    virtual int getNParticles()   const {return 0;} // TODO
    virtual int getNMobilities()  const = 0;
    virtual int getNConstraints() const = 0;    // i.e., constraint elements (multiple equations)

    virtual MobilizedBodyId           getParent   (MobilizedBodyId) const = 0;
    virtual Array<MobilizedBodyId>    getChildren (MobilizedBodyId) const = 0;
    //virtual const Mobilizer& getMobilizer(MobilizedBodyId) const = 0;

    virtual const Transform& getDefaultMobilizerFrame(MobilizedBodyId) const = 0;
    virtual const Transform& getDefaultMobilizerFrameOnParent(MobilizedBodyId) const = 0;
    virtual const MassProperties& getDefaultBodyMassProperties(MobilizedBodyId) const = 0;

        // MODEL STAGE //

    // Access to Instance variables. In general Mobilizers and Constraints will define more
    // of these; here we just deal with variables that are always present.

    virtual const MassProperties& getBodyMassProperties    (const State&, MobilizedBodyId) const = 0;
    virtual const Transform&      getMobilizerFrame        (const State&, MobilizedBodyId) const = 0;
    virtual const Transform&      getMobilizerFrameOnParent(const State&, MobilizedBodyId) const = 0;

    virtual const Vector& getAllParticleMasses     (const State&) const {
        static const Vector v;
        return v;
    }

    // These update routines invalidate Stage::Instance.
    virtual MassProperties& updBodyMassProperties    (State&, MobilizedBodyId) const = 0;
    virtual Transform&      updMobilizerFrame        (State&, MobilizedBodyId) const = 0;
    virtual Transform&      updMobilizerFrameOnParent(State&, MobilizedBodyId) const = 0;
    virtual Vector&         updAllParticleMasses     (State&) const = 0;

    // Access to Position and Velocity variables. //

    // These report the start index and number of generalized coordinates q or generalized
    // speeds u associated with a body's mobilizer. These are indices into this subsystem's
    // Q or U allocation in the state, and can also be used other arrays which have the
    // same dimensions. For example, mobilizer force arrays have the same dimension as
    // the Us (that is, their length is the total mobility of the system).
    virtual void findMobilizerQs(const State& s, MobilizedBodyId body, int& qStart, int& nq) const = 0;
    virtual void findMobilizerUs(const State& s, MobilizedBodyId body, int& uStart, int& nu) const = 0;

    // A subset of the Q's, not including particle coordinates.
    virtual const Vector& getAllMobilizerCoords(const State&) const = 0;
    virtual const Vector& getAllMobilizerSpeeds(const State&) const = 0;

    // Invalidates Stage::Position.
    virtual Vector& updAllMobilizerCoords(State&) const = 0;
    // Invalidates Stage::Velocity.
    virtual Vector& updAllMobilizerSpeeds(State&) const = 0;

    virtual const Vector_<Vec3>&  getAllParticleLocations (const State&) const {
        static const Vector_<Vec3> v;
        return v;
    }
    virtual const Vector_<Vec3>&  getAllParticleVelocities(const State&) const {
        static const Vector_<Vec3> v;
        return v;
    }


    // Invalidate Stage::Position.
    virtual Vector_<Vec3>&  updAllParticleLocations (State&) const {
        static Vector_<Vec3> v;
        return v;
    }
    // Invalidate Stage::Velocity.
    virtual Vector_<Vec3>&  updAllParticleVelocities(State&) const {
        static Vector_<Vec3> v;
        return v;
    }

    // Access to Acceleration variables. //

    virtual const Vector&              getAllMobilizerAppliedForces(const State&) const = 0;
    virtual const Vector_<Vec3>&       getAllParticleAppliedForces (const State&) const = 0;
    virtual const Vector_<SpatialVec>& getAllBodyAppliedForces     (const State&) const = 0;

    // These update routines invalidate Stage::Acceleration.
    virtual Vector&              updAllMobilizerAppliedForces(State&) const = 0;
    virtual Vector_<Vec3>&       updAllParticleAppliedForces (State&) const = 0;
    virtual Vector_<SpatialVec>& updAllBodyAppliedForces     (State&) const = 0;


        // INSTANCE STAGE //
    virtual Real getTotalMass(const State&) const = 0;

        // POSITION, VELOCITY, ACCELERATION STAGES //
    virtual const Transform&  getBodyTransform(const State&, MobilizedBodyId) const = 0;
    virtual const SpatialVec& getBodyVelocity(const State&, MobilizedBodyId) const = 0;
    virtual const SpatialVec& getBodyAcceleration(const State&, MobilizedBodyId) const = 0; 

    virtual const Vector_<Vec3>&  getAllParticleAccelerations(const State&) const {
        static const Vector_<Vec3> v;
        return v;
    }

    virtual Real calcQConstraintNorm(const State&) const {
        return 0;
    }
    virtual Real calcUConstraintNorm(const State&) const {
        return 0;
    }
    virtual Real calcUDotConstraintNorm(const State&) const {
        return 0;
    }
    virtual bool projectQConstraints(State&, Vector& y_err, Real tol, Real targetTol) const {
        return false;
    }
    virtual bool projectUConstraints(State&, Vector& y_err, Real tol, Real targetTol) const {
        return false;
    }

    SimTK_DOWNCAST(MatterSubsystemRep, SubsystemRep);
private:

};

} // namespace SimTK

#endif // SimTK_MATTER_SUBSYSTEM_REP_H_
