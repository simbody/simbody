#ifndef SimTK_SIMBODY_MULTIBODY_SYSTEM_REP_H_
#define SimTK_SIMBODY_MULTIBODY_SYSTEM_REP_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-9 Stanford University and the Authors.         *
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

/** @file
 * Define the private implementation of the MultibodySystem
 * class (a kind of System).
 */

#include "SimTKcommon.h"
#include "SimTKcommon/internal/SystemGuts.h"

#include "simbody/internal/common.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/ForceSubsystem.h"
#include "simbody/internal/DecorationSubsystem.h"
#include "simbody/internal/GeneralContactSubsystem.h"

#include "simbody/internal/ForceSubsystemGuts.h"
#include "SimbodyMatterSubsystemRep.h"

#include <vector>

namespace SimTK {

class AnalyticGeometry;
class DecorativeGeometry;

/**
 * This is a set of global variables for the MultibodySystem, stored as a cache
 * entry in the system's State, one at each Stage. Technically it is owned by the
 * MultibodySystemGlobalSubsystem, but it is manipulated directly by
 * the MultibodySystem.
 *
 * The entry at a given Stage includes all the contributions from the previous
 * Stages. For example, if some forces are generated at Stage::Model, those
 * are used to initialize the ones at Stage::Instance. That way when we get
 * to the final set at Stage::Dynamics we can use it directly to produce
 * accelerations. This structure allows us to invalidate a higher Stage without
 * having to recalculate forces that were known at a lower Stage.
 */
struct ForceCacheEntry {
    ForceCacheEntry() 
    { }
    // default copy constructor, copy assignment, destructor

    Vector_<SpatialVec> rigidBodyForces;
    Vector_<Vec3>       particleForces;
    Vector              mobilityForces;

    void ensureAllocatedTo(int nRigidBodies, int nParticles, int nMobilities) {
        rigidBodyForces.resize(nRigidBodies);
        particleForces.resize(nParticles);
        mobilityForces.resize(nMobilities);
    }

    void setAllForcesToZero() {
        rigidBodyForces.setToZero();
        particleForces.setToZero();
        mobilityForces.setToZero();
    }

    // This is just an assignment but allows for some bugcatchers. All the
    // ForceCacheEntries at every stage are supposed to have the same dimensions.
    void initializeFromSimilarForceEntry(const ForceCacheEntry& src) {
        assert(src.rigidBodyForces.size() == rigidBodyForces.size());
        assert(src.particleForces.size()  == particleForces.size());
        assert(src.mobilityForces.size()  == mobilityForces.size());
        *this = src;
    }
};

// Useless, but required by Value<T>.
inline std::ostream& operator<<(std::ostream& o, const ForceCacheEntry&) 
{assert(false);return o;}

/**
 * This is the subsystem used by a MultibodySystem to manage global state
 * calculations like forces and potential energy.
 */
class MultibodySystemGlobalSubsystemRep : public Subsystem::Guts {
    // Topological variables

    static const int NumForceCacheEntries = (Stage::DynamicsIndex-Stage::ModelIndex+1);
    mutable CacheEntryIndex forceCacheIndices[NumForceCacheEntries]; // where in state to find our stuff

    const ForceCacheEntry& getForceCacheEntry(const State& s, Stage g) const {
        assert(subsystemTopologyHasBeenRealized());
        SimTK_STAGECHECK_RANGE(Stage::Model, g, Stage::Dynamics,
            "MultibodySystem::getForceCacheEntry()");

        return Value<ForceCacheEntry>::downcast(
            getCacheEntry(s,forceCacheIndices[g-Stage::Model])).get();
    }
    ForceCacheEntry& updForceCacheEntry(const State& s, Stage g) const {
        assert(subsystemTopologyHasBeenRealized());
        SimTK_STAGECHECK_RANGE(Stage::Model, g, Stage::Dynamics,
            "MultibodySystem::getForceCacheEntry()");

        return Value<ForceCacheEntry>::downcast(
            updCacheEntry(s,forceCacheIndices[g-Stage::Model])).upd();
    }
public:
    MultibodySystemGlobalSubsystemRep()
      : Subsystem::Guts("MultibodySystemGlobalSubsystem", "0.0.2")
    {
        invalidateSubsystemTopologyCache();
    }

    const MultibodySystem& getMultibodySystem() const {
        return MultibodySystem::downcast(getSystem());
    }

    const Vector_<SpatialVec>& getRigidBodyForces(const State& s, Stage g) const {
        return getForceCacheEntry(s,g).rigidBodyForces;
    }
    const Vector_<Vec3>& getParticleForces(const State& s, Stage g) const {
        return getForceCacheEntry(s,g).particleForces;
    }
    const Vector& getMobilityForces(const State& s, Stage g) const {
        return getForceCacheEntry(s,g).mobilityForces;
    }

    Vector_<SpatialVec>& updRigidBodyForces(const State& s, Stage g) const {
        return updForceCacheEntry(s,g).rigidBodyForces;
    }
    Vector_<Vec3>& updParticleForces(const State& s, Stage g) const {
        return updForceCacheEntry(s,g).particleForces;
    }
    Vector& updMobilityForces(const State& s, Stage g) const {
        return updForceCacheEntry(s,g).mobilityForces;
    }

    // These override virtual methods from Subsystem::Guts.

    // Use default copy constructor, but then clear out the cache indices
    // and invalidate topology.
    MultibodySystemGlobalSubsystemRep* cloneImpl() const {
        MultibodySystemGlobalSubsystemRep* p = 
            new MultibodySystemGlobalSubsystemRep(*this);
        for (int i=0; i<NumForceCacheEntries; ++i)
            p->forceCacheIndices[i].invalidate();
        p->invalidateSubsystemTopologyCache();
        return p;
    }

    // At Topology stage we just allocate some slots in the State to hold
    // the forces. We can't initialize the force arrays because we don't yet
    // know the problem size.
    int realizeSubsystemTopologyImpl(State& s) const {
        const MultibodySystem& mbs = getMultibodySystem();

        for (Stage g(Stage::Model); g<=Stage::Dynamics; ++g)
            forceCacheIndices[g-Stage::Model] = 
                allocateCacheEntry(s, g, new Value<ForceCacheEntry>());

        return 0;
    }

    // At Model stage we know the problem size, so we can allocate the
    // model stage forces (if necessary) and initialize them (to zero).
    int realizeSubsystemModelImpl(State& s) const {
        const MultibodySystem&        mbs    = getMultibodySystem();
        const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();

        ForceCacheEntry& modelForces = updForceCacheEntry(s, Stage::Model);
        modelForces.ensureAllocatedTo(matter.getNumBodies(),
                                      matter.getNumParticles(),
                                      matter.getNumMobilities());
        modelForces.setAllForcesToZero();

        return 0;
    }

    // We treat the other stages like Model except that we use the 
    // previous Stage's ForceCacheEntry to initialize this one.
    int realizeSubsystemInstanceImpl(const State& s) const {
        const MultibodySystem&        mbs    = getMultibodySystem();
        const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();

        const ForceCacheEntry& modelForces = getForceCacheEntry(s, Stage::Model);
        ForceCacheEntry& instanceForces = updForceCacheEntry(s, Stage::Instance);
        instanceForces.ensureAllocatedTo(matter.getNumBodies(),
                                         matter.getNumParticles(),
                                         matter.getNumMobilities());
        instanceForces.initializeFromSimilarForceEntry(modelForces);

        return 0;
    }

    int realizeSubsystemTimeImpl(const State& s) const {
        const MultibodySystem&        mbs    = getMultibodySystem();
        const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();

        const ForceCacheEntry& instanceForces = getForceCacheEntry(s, Stage::Instance);
        ForceCacheEntry& timeForces = updForceCacheEntry(s, Stage::Time);
        timeForces.ensureAllocatedTo(matter.getNumBodies(),
                                     matter.getNumParticles(),
                                     matter.getNumMobilities());
        timeForces.initializeFromSimilarForceEntry(instanceForces);

        return 0;
    }

    int realizeSubsystemPositionImpl(const State& s) const {
        const MultibodySystem&        mbs    = getMultibodySystem();
        const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();

        const ForceCacheEntry& timeForces = getForceCacheEntry(s, Stage::Time);
        ForceCacheEntry& positionForces = updForceCacheEntry(s, Stage::Position);
        positionForces.ensureAllocatedTo(matter.getNumBodies(),
                                         matter.getNumParticles(),
                                         matter.getNumMobilities());
        positionForces.initializeFromSimilarForceEntry(timeForces);

        return 0;
    }

    int realizeSubsystemVelocityImpl(const State& s) const {
        const MultibodySystem&        mbs    = getMultibodySystem();
        const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();

        const ForceCacheEntry& positionForces = getForceCacheEntry(s, Stage::Position);
        ForceCacheEntry& velocityForces = updForceCacheEntry(s, Stage::Velocity);
        velocityForces.ensureAllocatedTo(matter.getNumBodies(),
                                         matter.getNumParticles(),
                                         matter.getNumMobilities());
        velocityForces.initializeFromSimilarForceEntry(positionForces);
        
        return 0;
    }

    int realizeSubsystemDynamicsImpl(const State& s) const {
        const MultibodySystem&        mbs    = getMultibodySystem();
        const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();

        const ForceCacheEntry& velocityForces = getForceCacheEntry(s, Stage::Velocity);
        ForceCacheEntry& dynamicsForces = updForceCacheEntry(s, Stage::Dynamics);
        dynamicsForces.ensureAllocatedTo(matter.getNumBodies(),
                                         matter.getNumParticles(),
                                         matter.getNumMobilities());
        dynamicsForces.initializeFromSimilarForceEntry(velocityForces);

        return 0;
    }

    // no need for other realize() methods
    SimTK_DOWNCAST(MultibodySystemGlobalSubsystemRep, Subsystem::Guts);
};

class MultibodySystemGlobalSubsystem : public Subsystem {
public:
    MultibodySystemGlobalSubsystem() : Subsystem() {
        adoptSubsystemGuts(new MultibodySystemGlobalSubsystemRep());
    }

    SimTK_PIMPL_DOWNCAST(MultibodySystemGlobalSubsystem, Subsystem);
    const MultibodySystemGlobalSubsystemRep& getRep() const;
    MultibodySystemGlobalSubsystemRep&       updRep();
};


/**
 * The job of the MultibodySystem class is to coordinate the activities of a
 * MatterSubsystem and a set of ForceSubsystems.
 */
class MultibodySystemRep : public System::Guts {
public:
    MultibodySystemRep() 
        : System::Guts("MultibodySystem", "0.0.1")
    {
    }
    ~MultibodySystemRep() {
    }

    SubsystemIndex setGlobalSubsystem() {
        assert(!globalSub.isValid());
        MultibodySystemGlobalSubsystem glo;
        globalSub = adoptSubsystem(glo);
        return globalSub;
    }
    SubsystemIndex setMatterSubsystem(SimbodyMatterSubsystem& m) {
        assert(!matterSub.isValid());
        matterSub = adoptSubsystem(m);
        return matterSub;
    }
    SubsystemIndex addForceSubsystem(ForceSubsystem& f) {
        forceSubs.push_back(adoptSubsystem(f));
        return forceSubs.back();
    }
    SubsystemIndex setDecorationSubsystem(DecorationSubsystem& d) {
        assert(!decorationSub.isValid());
        decorationSub = adoptSubsystem(d);
        return decorationSub;
    }
    SubsystemIndex setContactSubsystem(GeneralContactSubsystem& c) {
        assert(!contactSub.isValid());
        contactSub = adoptSubsystem(c);
        return contactSub;
    }

    const SimbodyMatterSubsystem& getMatterSubsystem() const {
        assert(matterSub.isValid());
        return SimbodyMatterSubsystem::downcast(getSubsystem(matterSub));
    }
    const ForceSubsystem& getForceSubsystem(SubsystemIndex id) const {
        return ForceSubsystem::downcast(getSubsystem(id));
    }
    const MultibodySystemGlobalSubsystem& getGlobalSubsystem() const {
        assert(globalSub.isValid());
        return MultibodySystemGlobalSubsystem::downcast(getSubsystem(globalSub));
    }

    bool hasDecorationSubsystem() const {return decorationSub.isValid();}
    bool hasContactSubsystem() const {return contactSub.isValid();}
    bool hasMatterSubsystem() const {return matterSub.isValid();}
    bool hasGlobalSubsystem() const {return globalSub.isValid();}

    const DecorationSubsystem& getDecorationSubsystem() const {
        assert(decorationSub.isValid());
        return DecorationSubsystem::downcast(getSubsystem(decorationSub));
    }

    const GeneralContactSubsystem& getContactSubsystem() const {
        assert(contactSub.isValid());
        return GeneralContactSubsystem::downcast(getSubsystem(contactSub));
    }

    SimbodyMatterSubsystem& updMatterSubsystem() {
        assert(matterSub.isValid());
        return SimbodyMatterSubsystem::updDowncast(updSubsystem(matterSub));
    }
    ForceSubsystem& updForceSubsystem(SubsystemIndex id) {
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
    GeneralContactSubsystem& updContactSubsystem() {
        assert(contactSub.isValid());
        return GeneralContactSubsystem::updDowncast(updSubsystem(contactSub));
    }

    // Global state cache entries dealing with interaction between forces & 
    // matter.

    // Responses available when the global subsystem is advanced to the
    // indicated stage or higher.
    const Vector_<SpatialVec>& getRigidBodyForces(const State& s, Stage g) const {
        return getGlobalSubsystem().getRep().getRigidBodyForces(s,g);
    }
    const Vector_<Vec3>& getParticleForces(const State& s, Stage g) const {
        return getGlobalSubsystem().getRep().getParticleForces(s,g);
    }
    const Vector& getMobilityForces(const State& s, Stage g) const {
        return getGlobalSubsystem().getRep().getMobilityForces(s,g);
    }
    const Real calcPotentialEnergy(const State& s) const {
        Real pe = 0;
        for (int i = 0; i < (int) forceSubs.size(); ++i)
            pe += getForceSubsystem(forceSubs[i]).getRep().calcPotentialEnergy(s);
        return pe;
    }

    // These are the global subsystem cache entries at the indicated Stage.
    // This will reduce the stage of s to the previous stage.
    Vector_<SpatialVec>& updRigidBodyForces(const State& s, Stage g) const {
        return getGlobalSubsystem().getRep().updRigidBodyForces(s,g);
    }
    Vector_<Vec3>& updParticleForces(const State& s, Stage g) const {
        return getGlobalSubsystem().getRep().updParticleForces(s,g);
    }
    Vector& updMobilityForces(const State& s, Stage g) const {
        return getGlobalSubsystem().getRep().updMobilityForces(s,g);
    }

    // pure virtual
    MultibodySystemRep* cloneImpl() const {return new MultibodySystemRep(*this);}

    // Override the SystemRep default implementations for these virtual methods.
    int realizeTopologyImpl    (State&)       const;
    int realizeModelImpl       (State&)       const;
    int realizeInstanceImpl    (const State&) const;
    int realizeTimeImpl        (const State&) const;
    int realizePositionImpl    (const State&) const;
    int realizeVelocityImpl    (const State&) const;
    int realizeDynamicsImpl    (const State&) const;
    int realizeAccelerationImpl(const State&) const;
    int realizeReportImpl      (const State&) const;

    // Currently prescribe() and project() affect only the Matter subsystem.

    int prescribeImpl(State& s, Stage g) {
        const SimbodyMatterSubsystem& mech = getMatterSubsystem();
        mech.getRep().prescribe(s,g);
        return 0;
    }

    // Note that we do all "q" projections before any "u" projections.
    //
    // TODO: yWeights & ooTols are being ignored here but shouldn't be!
    int projectImpl(State& s, Real consAccuracy, const Vector& yWeights,
                    const Vector& ooTols, Vector& yErrest, System::ProjectOptions opts) const
    {
        bool anyChange = false;

        realize(s, Stage::Time);

        const SimbodyMatterSubsystem& mech = getMatterSubsystem();

        if (opts.hasAnyPositionOptions()) {
            mech.getRep().realizeSubsystemPosition(s);
            if (mech.projectQConstraints(s, consAccuracy, yWeights, ooTols, yErrest, opts))
                anyChange = true;
        }

        realize(s, Stage::Position);  // realize the whole system now

        if (opts.hasAnyVelocityOptions()) {
            mech.getRep().realizeSubsystemVelocity(s);
            if (mech.projectUConstraints(s, consAccuracy, yWeights, ooTols, yErrest, opts))
                anyChange = true;
        }

        realize(s, Stage::Velocity);  // realize the whole system now

        return 0;
    }

    /* TODO: not yet
    virtual Real calcTimescaleImpl(const State&) const;
    virtual int calcYUnitWeightsImpl(const State&, Vector& weights) const;
    virtual int calcYErrUnitTolerancesImpl(const State&, Vector& tolerances) const;
    virtual int handleEventsImpl
       (State&, EventCause, const Array<EventId>& eventIds,
        Real accuracy, const Vector& yWeights, const Vector& ooConstraintTols,
        Stage& lowestModified, bool& shouldTerminate) const;
    virtual int calcEventTriggerInfoImpl(const State&, Array<EventTriggerInfo>&) const;
    virtual int calcTimeOfNextScheduledEventImpl
        (const State&, Real& tNextEvent, Array<EventId>& eventIds) const;
    */

    SimTK_DOWNCAST(MultibodySystemRep, System::Guts);
private:
    SubsystemIndex         globalSub;       // index of global subsystem
    SubsystemIndex         matterSub;       // index of matter subsystems
    Array_<SubsystemIndex> forceSubs;       // indices of force subsystems
    SubsystemIndex         decorationSub;   // index of DecorationSubsystem if any, else -1
    SubsystemIndex         contactSub;      // index of contact subsystem if any, else -1
};


} // namespace SimTK

#endif // SimTK_SIMBODY_MULTIBODY_SYSTEM_REP_H_
