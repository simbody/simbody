#ifndef SimTK_SIMBODY_MULTIBODY_SYSTEM_REP_H_
#define SimTK_SIMBODY_MULTIBODY_SYSTEM_REP_H_

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

/* Here we define the private implementation of the MultibodySystem class 
(a kind of System). This is not part of the Simbody API. */

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


//==============================================================================
//                              FORCE CACHE ENTRY
//==============================================================================
/* This is a set of global variables for the MultibodySystem, stored as a cache
entry in the system's State, one at each Stage. Technically it is owned by the
MultibodySystemGlobalSubsystem, but it is manipulated directly by
the MultibodySystem.

The entry at a given Stage includes all the contributions from the previous
Stages. For example, if some forces are generated at Stage::Model, those
are used to initialize the ones at Stage::Instance. That way when we get
to the final set at Stage::Dynamics we can use it directly to produce
accelerations. This structure allows us to invalidate a higher Stage without
having to recalculate forces that were known at a lower Stage. */
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



//==============================================================================
//                   MULTIBODY SYSTEM GLOBAL SUBSYSTEM REP
//==============================================================================
/* This is the subsystem used by a MultibodySystem to manage global state
calculations like forces and potential energy. */
class MultibodySystemGlobalSubsystemRep : public Subsystem::Guts {
    // Topological variables

    static const int NumForceCacheEntries = (Stage::Dynamics-Stage::Model+1);
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

        return Value<ForceCacheEntry>::updDowncast(
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



//==============================================================================
//                          MULTIBODY SYSTEM REP
//==============================================================================
/* The job of the MultibodySystem class is to coordinate the activities of a
MatterSubsystem and a set of ForceSubsystems. */
class MultibodySystemRep : public System::Guts {
public:
    MultibodySystemRep();
    ~MultibodySystemRep() {}

    SubsystemIndex setGlobalSubsystem() {
        assert(!m_globalSub.isValid());
        MultibodySystemGlobalSubsystem glo;
        m_globalSub = adoptSubsystem(glo);
        return m_globalSub;
    }
    SubsystemIndex setMatterSubsystem(SimbodyMatterSubsystem& m) {
        assert(!m_matterSub.isValid());
        m_matterSub = adoptSubsystem(m);
        return m_matterSub;
    }
    SubsystemIndex addForceSubsystem(ForceSubsystem& f) {
        m_forceSubs.push_back(adoptSubsystem(f));
        return m_forceSubs.back();
    }
    SubsystemIndex setDecorationSubsystem(DecorationSubsystem& d) {
        assert(!m_decorationSub.isValid());
        m_decorationSub = adoptSubsystem(d);
        return m_decorationSub;
    }
    SubsystemIndex setContactSubsystem(GeneralContactSubsystem& c) {
        assert(!m_contactSub.isValid());
        m_contactSub = adoptSubsystem(c);
        return m_contactSub;
    }

    const SimbodyMatterSubsystem& getMatterSubsystem() const {
        assert(m_matterSub.isValid());
        return SimbodyMatterSubsystem::downcast(getSubsystem(m_matterSub));
    }
    const ForceSubsystem& getForceSubsystem(SubsystemIndex id) const {
        return ForceSubsystem::downcast(getSubsystem(id));
    }
    const MultibodySystemGlobalSubsystem& getGlobalSubsystem() const {
        assert(m_globalSub.isValid());
        return MultibodySystemGlobalSubsystem::downcast
                                                    (getSubsystem(m_globalSub));
    }

    bool hasDecorationSubsystem() const {return m_decorationSub.isValid();}
    bool hasContactSubsystem() const {return m_contactSub.isValid();}
    bool hasMatterSubsystem() const {return m_matterSub.isValid();}
    bool hasGlobalSubsystem() const {return m_globalSub.isValid();}

    const DecorationSubsystem& getDecorationSubsystem() const {
        assert(m_decorationSub.isValid());
        return DecorationSubsystem::downcast(getSubsystem(m_decorationSub));
    }

    const GeneralContactSubsystem& getContactSubsystem() const {
        assert(m_contactSub.isValid());
        return GeneralContactSubsystem::downcast(getSubsystem(m_contactSub));
    }

    SimbodyMatterSubsystem& updMatterSubsystem() {
        assert(m_matterSub.isValid());
        return SimbodyMatterSubsystem::updDowncast(updSubsystem(m_matterSub));
    }
    ForceSubsystem& updForceSubsystem(SubsystemIndex id) {
        return ForceSubsystem::updDowncast(updSubsystem(id));
    }
    MultibodySystemGlobalSubsystem& updGlobalSubsystem() {
        assert(m_globalSub.isValid());
        return MultibodySystemGlobalSubsystem::updDowncast
                                                    (updSubsystem(m_globalSub));
    }
    DecorationSubsystem& updDecorationSubsystem() {
        assert(m_decorationSub.isValid());
        return DecorationSubsystem::updDowncast(updSubsystem(m_decorationSub));
    }
    GeneralContactSubsystem& updContactSubsystem() {
        assert(m_contactSub.isValid());
        return GeneralContactSubsystem::updDowncast(updSubsystem(m_contactSub));
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
        for (int i = 0; i < (int)m_forceSubs.size(); ++i)
            pe += getForceSubsystem(m_forceSubs[i]).getRep()
                                                        .calcPotentialEnergy(s);
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

    EventId getImpactEventId() const {return m_impactEventId;}
    EventId getContactChangeEventId() const {return m_contactChangeEventId;}

    // pure virtual
    MultibodySystemRep* cloneImpl() const override
    {   return new MultibodySystemRep(*this); }

    // Override the SystemRep default implementations for these virtual methods.
    int realizeTopologyImpl    (State&)       const override;
    int realizeModelImpl       (State&)       const override;
    int realizeInstanceImpl    (const State&) const override;
    int realizeTimeImpl        (const State&) const override;
    int realizePositionImpl    (const State&) const override;
    int realizeVelocityImpl    (const State&) const override;
    int realizeDynamicsImpl    (const State&) const override;
    int realizeAccelerationImpl(const State&) const override;
    int realizeReportImpl      (const State&) const override;


    void multiplyByNImpl(const State& s, const Vector& u, 
                         Vector& dq) const override {
        const SimbodyMatterSubsystem& mech = getMatterSubsystem();
        mech.getRep().multiplyByN(s,false,u,dq);
    }
    void multiplyByNTransposeImpl(const State& s, const Vector& fq, 
                                  Vector& fu) const override {
        const SimbodyMatterSubsystem& mech = getMatterSubsystem();
        mech.getRep().multiplyByN(s,true,fq,fu);
    }
    void multiplyByNPInvImpl(const State& s, const Vector& dq, 
                             Vector& u) const override {
        const SimbodyMatterSubsystem& mech = getMatterSubsystem();
        mech.getRep().multiplyByNInv(s,false,dq,u);
    }
    void multiplyByNPInvTransposeImpl(const State& s, const Vector& fu, 
                                      Vector& fq) const override {
        const SimbodyMatterSubsystem& mech = getMatterSubsystem();
        mech.getRep().multiplyByNInv(s,true,fu,fq);
    }  

    // Currently prescribe() and project() affect only the Matter subsystem.
    bool prescribeQImpl(State& state) const override {
        const SimbodyMatterSubsystem& mech = getMatterSubsystem();
        return mech.getRep().prescribeQ(state);
    }
    bool prescribeUImpl(State& state) const override {
        const SimbodyMatterSubsystem& mech = getMatterSubsystem();
        return mech.getRep().prescribeU(state);
    }

    void projectQImpl(State& state, Vector& qErrEst, 
                      const ProjectOptions& options, 
                      ProjectResults& results) const override {
        const SimbodyMatterSubsystem& mech = getMatterSubsystem();
        mech.getRep().projectQ(state, qErrEst, options, results);
        realize(state, Stage::Position);  // realize the whole system now
    }
    void projectUImpl(State& state, Vector& uErrEst, 
                      const ProjectOptions& options, 
                      ProjectResults& results) const override {
        const SimbodyMatterSubsystem& mech = getMatterSubsystem();
        mech.getRep().projectU(state, uErrEst, options, results);
        realize(state, Stage::Velocity);  // realize the whole system now
    }

    void getFreeQIndexImpl
       (const State& s, Array_<SystemQIndex>& freeQs) const override {
        const SimbodyMatterSubsystem& matter = getMatterSubsystem();
        const SystemQIndex qStart = matter.getQStart(s);
        const Array_<QIndex>& matterFreeQs = matter.getFreeQIndex(s);
        freeQs.resize(matterFreeQs.size());
        for (unsigned i=0; i < matterFreeQs.size(); ++i)
            freeQs[i] = SystemQIndex(qStart + matterFreeQs[i]);
    }
    void getFreeUIndexImpl
       (const State& s, Array_<SystemUIndex>& freeUs) const override {
        const SimbodyMatterSubsystem& matter = getMatterSubsystem();
        const SystemUIndex uStart = matter.getUStart(s);
        const Array_<UIndex>& matterFreeUs = matter.getFreeUIndex(s);
        freeUs.resize(matterFreeUs.size());
        for (unsigned i=0; i < matterFreeUs.size(); ++i)
            freeUs[i] = SystemUIndex(uStart + matterFreeUs[i]);
    }

    SimTK_DOWNCAST(MultibodySystemRep, System::Guts);

private:
friend class MultibodySystem;

    SubsystemIndex         m_globalSub;       // index of global subsystem
    SubsystemIndex         m_matterSub;       // index of matter subsystems
    Array_<SubsystemIndex> m_forceSubs;       // indices of force subsystems
    SubsystemIndex         m_decorationSub;   // DecorationSubsystem indx if any
    SubsystemIndex         m_contactSub;      // Contact subsystem index if any

    // These are IDs of events used to activate and deactivate conditional 
    // constraints.
    EventId     m_impactEventId;
    EventId     m_contactChangeEventId;
};


} // namespace SimTK

#endif // SimTK_SIMBODY_MULTIBODY_SYSTEM_REP_H_
