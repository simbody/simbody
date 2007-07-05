#ifndef SimTK_SIMBODY_MULTIBODY_SYSTEM_REP_H_
#define SimTK_SIMBODY_MULTIBODY_SYSTEM_REP_H_

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
 * Define the private implementation of the MultibodySystem
 * class (a kind of System).
 */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/MolecularMechanicsSystem.h"
#include "simbody/internal/MatterSubsystem.h"
#include "simbody/internal/ForceSubsystem.h"
#include "simbody/internal/DuMMForceFieldSubsystem.h"
#include "simbody/internal/DecorationSubsystem.h"
#include "simbody/internal/AnalyticGeometry.h"
#include "simbody/internal/DecorativeGeometry.h"

#include "SystemRep.h"
#include "ForceSubsystemRep.h"
#include "MatterSubsystemRep.h"

#include <vector>

namespace SimTK {

class AnalyticGeometry;
class DecorativeGeometry;

/**
 * This is a global variable for the MultibodySystem, stored as a cache
 * entry in the system's State. Technically it is owned by the
 * MultibodySystemGlobalSubsystem, but it is manipulated directly by
 * the MultibodySystem.
 */
struct DynamicsCache {
    DynamicsCache()
        : potentialEnergy(NTraits<Real>::NaN), kineticEnergy(NTraits<Real>::NaN)
    { }
    Vector_<SpatialVec> rigidBodyForces;
    Vector_<Vec3>       particleForces;
    Vector              mobilityForces;
    Real                potentialEnergy;
    Real                kineticEnergy;
};

// Useless, but required by Value<T>.
inline std::ostream& operator<<(std::ostream& o, const DynamicsCache&) 
{assert(false);return o;}

/**
 * This is the subsystem used by a MultibodySystem to manage global state
 * calculations like forces and potential energy.
 */
class MultibodySystemGlobalSubsystemRep : public SubsystemRep {
    // Topological variables

    mutable int dynamicsCacheIndex;         // where in state to find our stuff

    const DynamicsCache& getDynamicsCache(const State& s) const {
        assert(subsystemTopologyHasBeenRealized());
        return Value<DynamicsCache>::downcast(
            getCacheEntry(s,dynamicsCacheIndex)).get();
    }
    DynamicsCache& updDynamicsCache(const State& s) const {
        assert(subsystemTopologyHasBeenRealized());
        return Value<DynamicsCache>::downcast(
            updCacheEntry(s,dynamicsCacheIndex)).upd();
    }
public:
    MultibodySystemGlobalSubsystemRep()
      : SubsystemRep("MultibodySystemGlobalSubsystem", "0.0.1"),
        dynamicsCacheIndex(-1)
    {
    }

    const MultibodySystem& getMultibodySystem() const {
        return MultibodySystem::downcast(getSystem());
    }

    const Vector_<SpatialVec>& getRigidBodyForces(const State& s) const {
        return getDynamicsCache(s).rigidBodyForces;
    }
    const Vector_<Vec3>& getParticleForces(const State& s) const {
        return getDynamicsCache(s).particleForces;
    }
    const Vector& getMobilityForces(const State& s) const {
        return getDynamicsCache(s).mobilityForces;
    }
    const Real& getPotentialEnergy(const State& s) const {
        return getDynamicsCache(s).potentialEnergy;
    }
    const Real& getKineticEnergy(const State& s) const {
        return getDynamicsCache(s).kineticEnergy;
    }
    
    Vector_<SpatialVec>& updRigidBodyForces(const State& s) const {
        return updDynamicsCache(s).rigidBodyForces;
    }
    Vector_<Vec3>& updParticleForces(const State& s) const {
        return updDynamicsCache(s).particleForces;
    }
    Vector& updMobilityForces(const State& s) const {
        return updDynamicsCache(s).mobilityForces;
    }
    Real& updPotentialEnergy(const State& s) const {
        return updDynamicsCache(s).potentialEnergy;
    }
    Real& updKineticEnergy(const State& s) const {
        return updDynamicsCache(s).kineticEnergy;
    }

    MultibodySystemGlobalSubsystemRep* cloneSubsystemRep() const {
        return new MultibodySystemGlobalSubsystemRep(*this);
    }

    // These override virtual methods from SubsystemRep.

    void realizeSubsystemTopologyImpl(State& s) const {
        const MultibodySystem& mbs = getMultibodySystem();
        dynamicsCacheIndex = s.allocateCacheEntry(getMySubsystemId(), 
            Stage::Dynamics, new Value<DynamicsCache>());
    }

    // At dynamics stage we make sure the force arrays are all the right length
    // and then initialize them and energy to zero. We expect the
    // force subsystems to then fill in forces and potential energy
    // as they are realized to this stage, and the matter subsystem to
    // fill in the kinetic energy.
    void realizeSubsystemDynamicsImpl(const State& s) const {
        DynamicsCache& dc = updDynamicsCache(s);
        dc.potentialEnergy = dc.kineticEnergy = 0;
        const MultibodySystem& mbs = getMultibodySystem();
        mbs.getMatterSubsystem().resetForces(dc.rigidBodyForces,
                                             dc.particleForces,
                                             dc.mobilityForces);
    }


    // no need for other realize() methods
    SimTK_DOWNCAST(MultibodySystemGlobalSubsystemRep, SubsystemRep);
};

class MultibodySystemGlobalSubsystem : public Subsystem {
public:
    MultibodySystemGlobalSubsystem() : Subsystem() {
        rep = new MultibodySystemGlobalSubsystemRep();
        rep->setMyHandle(*this);
    }

    SimTK_PIMPL_DOWNCAST(MultibodySystemGlobalSubsystem, Subsystem);
    const MultibodySystemGlobalSubsystemRep& getRep() const;
    MultibodySystemGlobalSubsystemRep&       updRep();
};


/**
 * The job of the MultibodySystem class is to coordinate the activities of a
 * MatterSubsystem and a set of ForceSubsystems.
 */
class MultibodySystemRep : public SystemRep {
public:
    MultibodySystemRep() 
      : SystemRep("MultibodySystem", "0.0.1")
    {
    }
    ~MultibodySystemRep() {
    }

    SubsystemId setGlobalSubsystem() {
        assert(!globalSub.isValid());
        MultibodySystemGlobalSubsystem glo;
        globalSub = adoptSubsystem(glo);
        return globalSub;
    }
    SubsystemId setMatterSubsystem(MatterSubsystem& m) {
        assert(!matterSub.isValid());
        matterSub = adoptSubsystem(m);
        return matterSub;
    }
    SubsystemId addForceSubsystem(ForceSubsystem& f) {
        forceSubs.push_back(adoptSubsystem(f));
        return forceSubs.back();
    }
    SubsystemId setDecorationSubsystem(DecorationSubsystem& d) {
        assert(!decorationSub.isValid());
        decorationSub = adoptSubsystem(d);
        return decorationSub;
    }

    const MatterSubsystem& getMatterSubsystem() const {
        assert(matterSub.isValid());
        return MatterSubsystem::downcast(getSubsystem(matterSub));
    }
    const ForceSubsystem& getForceSubsystem(SubsystemId id) const {
        return ForceSubsystem::downcast(getSubsystem(id));
    }
    const MultibodySystemGlobalSubsystem& getGlobalSubsystem() const {
        assert(globalSub.isValid());
        return MultibodySystemGlobalSubsystem::downcast(getSubsystem(globalSub));
    }
    bool hasDecorationSubsystem() const {return decorationSub.isValid();}
    const DecorationSubsystem& getDecorationSubsystem() const {
        assert(decorationSub.isValid());
        return DecorationSubsystem::downcast(getSubsystem(decorationSub));
    }

    MatterSubsystem& updMatterSubsystem() {
        assert(matterSub.isValid());
        return MatterSubsystem::updDowncast(updSubsystem(matterSub));
    }
    ForceSubsystem& updForceSubsystem(SubsystemId id) {
        return ForceSubsystem::updDowncast(updSubsystem(id));
    }
    MultibodySystemGlobalSubsystem& updGlobalSubsystem() {
        assert(globalSub.isValid());
        return MultibodySystemGlobalSubsystem::updDowncast(updSubsystem(globalSub));
    }
    DecorationSubsystem& updDecorationSubsystem() {
        assert(decorationSub.isValid());
        return DecorationSubsystem::updDowncast(updSubsystem(decorationSub));
    }
    // Global state variables dealing with interaction between forces & matter

    // Responses available when the global subsystem is advanced to Dynamics stage.
    const Vector_<SpatialVec>& getRigidBodyForces(const State& s) const {
        return getGlobalSubsystem().getRep().getRigidBodyForces(s);
    }
    const Vector_<Vec3>& getParticleForces(const State& s) const {
        return getGlobalSubsystem().getRep().getParticleForces(s);
    }
    const Vector& getMobilityForces(const State& s) const {
        return getGlobalSubsystem().getRep().getMobilityForces(s);
    }
    const Real& getPotentialEnergy(const State& s) const {
        return getGlobalSubsystem().getRep().getPotentialEnergy(s);
    }
    const Real& getKineticEnergy(const State& s) const {
        return getGlobalSubsystem().getRep().getKineticEnergy(s);
    }

    // Dynamics stage cache entries of the global subsystem.
    Vector_<SpatialVec>& updRigidBodyForces(const State& s) const {
        return getGlobalSubsystem().getRep().updRigidBodyForces(s);
    }
    Vector_<Vec3>& updParticleForces(const State& s) const {
        return getGlobalSubsystem().getRep().updParticleForces(s);
    }
    Vector& updMobilityForces(const State& s) const {
        return getGlobalSubsystem().getRep().updMobilityForces(s);
    }
    Real& updPotentialEnergy(const State& s) const {
        return getGlobalSubsystem().getRep().updPotentialEnergy(s);
    }
    Real& updKineticEnergy(const State& s) const {
        return getGlobalSubsystem().getRep().updKineticEnergy(s);
    }

    // Do all "q" projections before any "u" projections.
    bool project(State& s, Vector& y_err,
                 const Real& tol, const Real& dontProjectFac, 
                 const Real& targetTol) const 
    {
        bool anyChange = false;

        realize(s, Stage::Time);

        const MatterSubsystem& mech = getMatterSubsystem();

        mech.getRep().realizeSubsystemPosition(s);
        const Real qerr = mech.calcQConstraintNorm(s);
        if (dontProjectFac==0 || qerr > tol*dontProjectFac) {
            if (mech.projectQConstraints(s, y_err, tol, targetTol))
                anyChange = true;
        }

        realize(s, Stage::Position);  // realize the whole system now

        mech.getRep().realizeSubsystemVelocity(s);
        const Real uerr = mech.calcUConstraintNorm(s);
        if (dontProjectFac==0 || uerr > tol*dontProjectFac) {
            if (mech.projectUConstraints(s, y_err, tol, targetTol))
                anyChange = true;
        }

        realize(s, Stage::Velocity);  // realize the whole system now

        return anyChange;
    }

    // pure virtual
    MultibodySystemRep* cloneSystemRep() const {return new MultibodySystemRep(*this);}

    // Override the SystemRep default implementations for these virtual methods.
    void realizeTopologyImpl    (State& s)       const;
    void realizeModelImpl       (State& s)       const;
    void realizeInstanceImpl    (const State& s) const;
    void realizeTimeImpl        (const State& s) const;
    void realizePositionImpl    (const State& s) const;
    void realizeVelocityImpl    (const State& s) const;
    void realizeDynamicsImpl    (const State& s) const;
    void realizeAccelerationImpl(const State& s) const;
    void realizeReportImpl      (const State& s) const;

    SimTK_DOWNCAST(MultibodySystemRep, SystemRep);
private:
    SubsystemId  globalSub;             // index of global subsystem
    SubsystemId  matterSub;             // index of matter subsystems
    std::vector<SubsystemId> forceSubs; // indices of force subsystems
    SubsystemId  decorationSub;         // index of DecorationSubsystem if any, else -1
};


/**
 * This class is a kind of MultibodySystem which is required to have exactly
 * one DuMMForceFieldSubsystem.
 */
class MolecularMechanicsSystemRep : public MultibodySystemRep {
public:
    MolecularMechanicsSystemRep() : MultibodySystemRep()
    {
    }
    ~MolecularMechanicsSystemRep() {
    }

    int setMolecularMechanicsForceSubsystem(DuMMForceFieldSubsystem& mm) {
        assert(!molecularMechanicsSub.isValid());
        molecularMechanicsSub = addForceSubsystem(mm);
        return molecularMechanicsSub;
    }

    const DuMMForceFieldSubsystem& getMolecularMechanicsForceSubsystem() const {
        assert(molecularMechanicsSub.isValid());
        return DuMMForceFieldSubsystem::downcast(getSubsystem(molecularMechanicsSub));
    }
    DuMMForceFieldSubsystem& updMolecularMechanicsForceSubsystem() {
        assert(molecularMechanicsSub.isValid());
        return DuMMForceFieldSubsystem::updDowncast(updSubsystem(molecularMechanicsSub));
    }

    SimTK_DOWNCAST(MolecularMechanicsSystemRep, SystemRep);
private:
    SubsystemId molecularMechanicsSub;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_MULTIBODY_SYSTEM_REP_H_
