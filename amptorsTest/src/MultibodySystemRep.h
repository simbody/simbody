#ifndef SimTK_SIMBODY_MULTIBODY_SYSTEM_REP_H_
#define SimTK_SIMBODY_MULTIBODY_SYSTEM_REP_H_

/* Portions copyright (c) 2005-6 Stanford University and Michael Sherman.
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
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
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
    mutable bool built;

    const DynamicsCache& getDynamicsCache(const State& s) const {
        assert(built);
        return Value<DynamicsCache>::downcast(
            getCacheEntry(s,dynamicsCacheIndex)).get();
    }
    DynamicsCache& updDynamicsCache(const State& s) const {
        assert(built);
        return Value<DynamicsCache>::downcast(
            updCacheEntry(s,dynamicsCacheIndex)).upd();
    }
public:
    MultibodySystemGlobalSubsystemRep()
      : SubsystemRep("MultibodySystemGlobalSubsystem", "0.0.1"),
        dynamicsCacheIndex(-1), built(false)
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


    void realizeTopology(State& s) const {
        const MultibodySystem& mbs = getMultibodySystem();
        dynamicsCacheIndex = s.allocateCacheEntry(getMySubsystemIndex(), 
            Stage::Dynamics, new Value<DynamicsCache>());

        built = true;
    }

    // At dynamics stage we make sure the force arrays are all the right length
    // and then initialize them and energy to zero. We expect the
    // force subsystems to then fill in forces and potential energy
    // as they are realized to this stage, and the matter subsystem to
    // fill in the kinetic energy.
    void realizeDynamics(const State& s) const {
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
      : SystemRep("MultibodySystem", "0.0.1"), 
        built(false), globalSub(-1), matterSub(-1), decorationSub(-1)
    {
    }
    ~MultibodySystemRep() {
    }

    int setGlobalSubsystem() {
        assert(globalSub == -1);
        built = false;
        MultibodySystemGlobalSubsystem glo;
        globalSub = takeOverSubsystem(glo);
        return globalSub;
    }
    int setMatterSubsystem(MatterSubsystem& m) {
        assert(matterSub == -1);
        built = false;
        matterSub = takeOverSubsystem(m);
        return matterSub;
    }
    int addForceSubsystem(ForceSubsystem& f) {
        built = false;
        forceSubs.push_back(takeOverSubsystem(f));
        return forceSubs.back();
    }
    int setDecorationSubsystem(DecorationSubsystem& d) {
        assert(decorationSub == -1);
        built = false;
        decorationSub = takeOverSubsystem(d);
        return decorationSub;
    }

    const MatterSubsystem& getMatterSubsystem() const {
        assert(matterSub >= 0);
        return MatterSubsystem::downcast(getSubsystem(matterSub));
    }
    const ForceSubsystem& getForceSubsystem(int i) const {
        assert(i >= 0);
        return ForceSubsystem::downcast(getSubsystem(i));
    }
    const MultibodySystemGlobalSubsystem& getGlobalSubsystem() const {
        assert(globalSub >= 0);
        return MultibodySystemGlobalSubsystem::downcast(getSubsystem(globalSub));
    }
    bool hasDecorationSubsystem() const {return decorationSub >= 0;}
    const DecorationSubsystem& getDecorationSubsystem() const {
        assert(decorationSub >= 0);
        return DecorationSubsystem::downcast(getSubsystem(decorationSub));
    }

    MatterSubsystem& updMatterSubsystem() {
        assert(matterSub >= 0);
        return MatterSubsystem::updDowncast(updSubsystem(matterSub));
    }
    ForceSubsystem& updForceSubsystem(int i) {
        assert(i >= 0);
        return ForceSubsystem::updDowncast(updSubsystem(i));
    }
    MultibodySystemGlobalSubsystem& updGlobalSubsystem() {
        assert(globalSub >= 0);
        return MultibodySystemGlobalSubsystem::updDowncast(updSubsystem(globalSub));
    }
    DecorationSubsystem& updDecorationSubsystem() {
        assert(decorationSub >= 0);
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

        mech.realize(s, Stage::Position);
        const Real qerr = mech.calcQConstraintNorm(s);
        if (dontProjectFac==0 || qerr > tol*dontProjectFac) {
            if (mech.projectQConstraints(s, y_err, tol, targetTol))
                anyChange = true;
        }

        realize(s, Stage::Position);  // realize the whole system now

        mech.realize(s, Stage::Velocity);
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

    void realizeTopology(State& s) const {
        MultibodySystemRep& mutableThis = *const_cast<MultibodySystemRep*>(this);

        assert(globalSub != -1);
        assert(matterSub != -1);

        getGlobalSubsystem().realize(s, Stage::Topology);
        getMatterSubsystem().realize(s, Stage::Topology);
        for (int i=0; i < (int)forceSubs.size(); ++i)
            getForceSubsystem(forceSubs[i]).realize(s, Stage::Topology);

        if (hasDecorationSubsystem())
            getDecorationSubsystem().realize(s, Stage::Topology);

        mutableThis.built = true;
    }
    void realizeModel(State& s) const {
        getGlobalSubsystem().realize(s, Stage::Model);
        getMatterSubsystem().realize(s, Stage::Model);
        for (int i=0; i < (int)forceSubs.size(); ++i)
            getForceSubsystem(forceSubs[i]).realize(s, Stage::Model);
 
        if (hasDecorationSubsystem())
            getDecorationSubsystem().realize(s, Stage::Model);
   }
    void realizeInstance(const State& s) const {
        getGlobalSubsystem().realize(s, Stage::Instance);
        getMatterSubsystem().realize(s, Stage::Instance);
        for (int i=0; i < (int)forceSubs.size(); ++i)
            getForceSubsystem(forceSubs[i]).realize(s, Stage::Instance);
 
        if (hasDecorationSubsystem())
            getDecorationSubsystem().realize(s, Stage::Instance);
   }
    void realizeTime(const State& s) const {
        getGlobalSubsystem().realize(s, Stage::Time);
        getMatterSubsystem().realize(s, Stage::Time);
        for (int i=0; i < (int)forceSubs.size(); ++i)
            getForceSubsystem(forceSubs[i]).realize(s, Stage::Time);

        if (hasDecorationSubsystem())
            getDecorationSubsystem().realize(s, Stage::Time);
    }
    void realizePosition(const State& s) const {
        getGlobalSubsystem().realize(s, Stage::Position);
        getMatterSubsystem().realize(s, Stage::Position);
        for (int i=0; i < (int)forceSubs.size(); ++i)
            getForceSubsystem(forceSubs[i]).realize(s, Stage::Position);

        if (hasDecorationSubsystem())
            getDecorationSubsystem().realize(s, Stage::Position);
    }
    void realizeVelocity(const State& s) const {
        getGlobalSubsystem().realize(s, Stage::Velocity);
        getMatterSubsystem().realize(s, Stage::Velocity);
        for (int i=0; i < (int)forceSubs.size(); ++i)
            getForceSubsystem(forceSubs[i]).realize(s, Stage::Velocity);
 
        if (hasDecorationSubsystem())
            getDecorationSubsystem().realize(s, Stage::Velocity);
   }
    void realizeDynamics(const State& s) const {
        getGlobalSubsystem().realize(s, Stage::Dynamics);
        // note order: forces first (TODO: does that matter?)
        for (int i=0; i < (int)forceSubs.size(); ++i)
            getForceSubsystem(forceSubs[i]).realize(s, Stage::Dynamics);
        getMatterSubsystem().realize(s, Stage::Dynamics);

        if (hasDecorationSubsystem())
            getDecorationSubsystem().realize(s, Stage::Dynamics);
    }
    void realizeAcceleration(const State& s) const {
        getGlobalSubsystem().realize(s, Stage::Acceleration);
        // note order: forces first (TODO: does that matter?)
        for (int i=0; i < (int)forceSubs.size(); ++i)
            getForceSubsystem(forceSubs[i]).realize(s, Stage::Acceleration);
        getMatterSubsystem().realize(s, Stage::Acceleration);

        if (hasDecorationSubsystem())
            getDecorationSubsystem().realize(s, Stage::Acceleration);
    }

    SimTK_DOWNCAST(MultibodySystemRep, SystemRep);
private:
    bool built;
    int  globalSub;             // index of global subsystem
    int  matterSub;             // index of matter subsystems
    std::vector<int> forceSubs; // indices of force subsystems
    int  decorationSub;         // index of DecorationSubsystem if any, else -1
};


/**
 * This class is a kind of MultibodySystem which is required to have exactly
 * one DuMMForceFieldSubsystem.
 */
class MolecularMechanicsSystemRep : public MultibodySystemRep {
public:
    MolecularMechanicsSystemRep() 
      : MultibodySystemRep(), molecularMechanicsSub(-1)
    {
    }
    ~MolecularMechanicsSystemRep() {
    }

    int setMolecularMechanicsForceSubsystem(DuMMForceFieldSubsystem& mm) {
        assert(molecularMechanicsSub == -1);
        molecularMechanicsSub = addForceSubsystem(mm);
        return molecularMechanicsSub;
    }

    const DuMMForceFieldSubsystem& getMolecularMechanicsForceSubsystem() const {
        assert(molecularMechanicsSub >= 0);
        return DuMMForceFieldSubsystem::downcast(getSubsystem(molecularMechanicsSub));
    }
    DuMMForceFieldSubsystem& updMolecularMechanicsForceSubsystem() {
        assert(molecularMechanicsSub >= 0);
        return DuMMForceFieldSubsystem::updDowncast(updSubsystem(molecularMechanicsSub));
    }

    SimTK_DOWNCAST(MolecularMechanicsSystemRep, SystemRep);
private:
    int molecularMechanicsSub;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_MULTIBODY_SYSTEM_REP_H_
