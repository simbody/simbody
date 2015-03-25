/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-13 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors: Peter Eastman, Nabeel Allana, Chris Dembia                   *
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

/**@file
 *
 * Private implementation of GeneralForceSubsystem.
 */

#include "SimTKcommon.h"

#include "simbody/internal/common.h"
#include "simbody/internal/ForceSubsystem.h"
#include "simbody/internal/ForceSubsystemGuts.h"
#include "simbody/internal/GeneralForceSubsystem.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/MultibodySystem.h"

#include "ForceImpl.h"

#include <memory>

namespace SimTK {

// TODO create calcForceExecutor.
// TODO modify test case.
// TODO shouldBeParallelized.
// TODO choosing number of threads.
// TODO should we allow disabling parallelism?
        // The "+ 1" is for the non-parallel forces.
        /*
        const int numTasks = enabledParallelForces.size() + 1;
            calcForceTask.calcForceCachedAndNonCached(
                    enabledNonParallelForces, enabledParallelForces,
                    rigidBodyForces, particleForces, mobilityForces,
                    rigidBodyForceCache, particleForceCache,
                    mobilityForceCache);
            calcForceTask.calcForceNonCached(
                    enabledNonParallelForces, enabledParallelForces,
                    rigidBodyForces, particleForces, mobilityForces);
                    */

/** Calls calcForce() on all enabled forces. For an index of 0, this Task
 calculates the force for all non-parallelized forces. For all other indices,
 this Task executes a single parallelized force. */
class CalcForcesTask : public ParallelExecutor::Task {
public:

    enum Mode {
        All,
        CachedAndNonCached,
        NonCached
    };

    CalcForcesTask(ParallelExecutor& executor)
        : m_executor(executor), m_mode(All) {}

    void calcForceAll(
            const State& s,
            const Array_<Force*>& enabledNonParallelForces,
            const Array_<Force*>& enabledParallelForces,
            Vector_<SpatialVec>& rigidBodyForces,
            Vector_<Vec3>& particleForces,
            Vector& mobilityForces) {
        m_s = &s;
        m_enabledNonParallelForces = &enabledNonParallelForces;
        m_enabledParallelForces = &enabledParallelForces;
        m_rigidBodyForces = &rigidBodyForces;
        m_particleForces = &particleForces;
        m_mobilityForces = &mobilityForces;
        // The "+ 1" is for the non-parallelized forces.
        const int numTimes = enabledParallelForces.size() + 1;
        m_mode = All;
        m_executor.execute(*this, numTimes);
    }
    void calcForceCachedAndNonCached(
            const State& s,
            const Array_<Force*>& enabledNonParallelForces,
            const Array_<Force*>& enabledParallelForces,
            Vector_<SpatialVec>& rigidBodyForces,
            Vector_<Vec3>& particleForces,
            Vector& mobilityForces,
            Vector_<SpatialVec>& rigidBodyForceCache,
            Vector_<Vec3>& particleForceCache,
            Vector& mobilityForceCache) {
        m_s = &s;
        m_enabledNonParallelForces = &enabledNonParallelForces;
        m_enabledParallelForces = &enabledParallelForces;
        m_rigidBodyForces = &rigidBodyForces;
        m_particleForces = &particleForces;
        m_mobilityForces = &mobilityForces;
        m_rigidBodyForceCache = &rigidBodyForceCache;
        m_particleForceCache = &particleForceCache;
        m_mobilityForceCache = &mobilityForceCache;
        // The "+ 1" is for the non-parallelized forces.
        const int numTimes = enabledParallelForces.size() + 1;
        m_mode = CachedAndNonCached;
        m_executor.execute(*this, numTimes);
    }
    void calcForceNonCached(
            const State& s,
            const Array_<Force*>& enabledNonParallelForces,
            const Array_<Force*>& enabledParallelForces,
            Vector_<SpatialVec>& rigidBodyForces,
            Vector_<Vec3>& particleForces,
            Vector& mobilityForces) {
        m_s = &s;
        m_enabledNonParallelForces = &enabledNonParallelForces;
        m_enabledParallelForces = &enabledParallelForces;
        m_rigidBodyForces = &rigidBodyForces;
        m_particleForces = &particleForces;
        m_mobilityForces = &mobilityForces;
        // The "+ 1" is for the non-parallelized forces.
        const int numTimes = enabledParallelForces.size() + 1;
        m_mode = NonCached;
        m_executor.execute(*this, numTimes);
    }

    void execute(int index) override {
        switch (m_mode) {
        case All:
            if (index == 0) {
                // Process all non-parallel forces.
                for (Force* force : *m_enabledNonParallelForces) {
                    force->getImpl().calcForce(*m_s,
                            *m_rigidBodyForces, *m_particleForces,
                            *m_mobilityForces);
                }
            } else {
                // Process a single parallel force. Subtract 1 from index b/c
                // we use 0 for the non-parallel forces.
                const auto& impl =
                    m_enabledParallelForces->getElt(index-1)->getImpl();
                impl.calcForce(*m_s,
                        *m_rigidBodyForces, *m_particleForces,
                        *m_mobilityForces);
            }
            break;
        
        // TODO
        case CachedAndNonCached:
            if (index == 0) {
                // Process all non-parallel forces.
                for (Force* force : *m_enabledNonParallelForces) {
                    const auto& impl = force->getImpl();
                    if (impl.dependsOnlyOnPositions()) {
                        impl.calcForce(*m_s,
                                *m_rigidBodyForceCache, *m_particleForceCache,
                                *m_mobilityForceCache);
                    } else { // ordinary velocity dependent force
                        impl.calcForce(*m_s,
                                *m_rigidBodyForces, *m_particleForces,
                                *m_mobilityForces);
                    }
                }
            } else {
                // Process a single parallel force. Subtract 1 from index b/c
                // we use 0 for the non-parallel forces.
                const auto& impl =
                    m_enabledParallelForces->getElt(index-1)->getImpl();
                if (impl.dependsOnlyOnPositions()) {
                    impl.calcForce(*m_s,
                            *m_rigidBodyForceCache, *m_particleForceCache,
                            *m_mobilityForceCache);
                } else { // ordinary velocity dependent force
                    impl.calcForce(*m_s,
                            *m_rigidBodyForces, *m_particleForces,
                            *m_mobilityForces);
                }
            }
            break;

        // TODO
        case NonCached:
            if (index == 0) {
                // Process all non-parallel forces.
                for (Force* force : *m_enabledNonParallelForces) {
                    const auto& impl = force->getImpl();
                    if (!impl.dependsOnlyOnPositions()) {
                        impl.calcForce(*m_s,
                                *m_rigidBodyForces, *m_particleForces,
                                *m_mobilityForces);
                    }
                }
            } else {
                // Process a single parallel force. Subtract 1 from index b/c
                // we use 0 for the non-parallel forces.
                const auto& impl =
                    m_enabledParallelForces->getElt(index-1)->getImpl();
                    if (!impl.dependsOnlyOnPositions()) {
                        impl.calcForce(*m_s,
                                *m_rigidBodyForces, *m_particleForces,
                                *m_mobilityForces);
                    }
            }
            break;
        }
    }
    void initialize() override {
        // The ParallelExecutor is being executed (again). Must set these to 0
        // so that we can properly accumulate forces.
        // TODO Sherm uses "= Vec3(0)" etc. below; is that what I should do
        // here? Also, I'm not resizing; could the size change?
        m_rigidBodyForcesLocal.upd().resize(m_rigidBodyForces->size());
        m_rigidBodyForcesLocal.upd().setToZero();
        m_particleForcesLocal.upd().resize(m_particleForces->size());
        m_particleForcesLocal.upd().setToZero();
        m_mobilityForcesLocal.upd().resize(m_mobilityForces->size());
        m_mobilityForcesLocal.upd().setToZero();

        if (m_mode == CachedAndNonCached) {
            m_rigidBodyForceCacheLocal.upd().resize(m_rigidBodyForceCache->size());
            m_rigidBodyForceCacheLocal.upd().setToZero();
            m_particleForceCacheLocal.upd().resize(m_particleForceCache->size());
            m_particleForceCacheLocal.upd().setToZero();
            m_mobilityForceCacheLocal.upd().resize(m_mobilityForceCache->size());
            m_mobilityForceCacheLocal.upd().setToZero();
        }
    }
    void finish() override {
        // Add in this thread's contribution.
        *m_rigidBodyForces += m_rigidBodyForcesLocal.get();
        *m_particleForces += m_particleForcesLocal.get();
        *m_mobilityForces += m_mobilityForcesLocal.get();

        if (m_mode == CachedAndNonCached) {
            *m_rigidBodyForceCache += m_rigidBodyForceCacheLocal.get();
            *m_particleForceCache += m_particleForceCacheLocal.get();
            *m_mobilityForceCache += m_mobilityForceCacheLocal.get();
        }
    }

private:

    ParallelExecutor& m_executor;

    Mode m_mode;

    const State* m_s;

    // These are separated since we handle them differently.
    const Array_<Force*>* m_enabledNonParallelForces;
    const Array_<Force*>* m_enabledParallelForces;

    // We eventually add our thread-local result to these vectors.
    Vector_<SpatialVec>* m_rigidBodyForces;
    Vector_<Vec3>* m_particleForces;
    Vector* m_mobilityForces;

    Vector_<SpatialVec>* m_rigidBodyForceCache;
    Vector_<Vec3>* m_particleForceCache;
    Vector* m_mobilityForceCache;

    // These variables are local to a thread. They are set to their default
    // value when the threads are spawned.
    ThreadLocal<Vector_<SpatialVec>> m_rigidBodyForcesLocal;
    ThreadLocal<Vector_<Vec3>> m_particleForcesLocal;
    ThreadLocal<Vector> m_mobilityForcesLocal;

    ThreadLocal<Vector_<SpatialVec>> m_rigidBodyForceCacheLocal;
    ThreadLocal<Vector_<Vec3>> m_particleForceCacheLocal;
    ThreadLocal<Vector> m_mobilityForceCacheLocal;
};

// There is some tricky caching being done here for forces that have overridden
// dependsOnlyOnPositions() (and returned "true"). This is probably only worth
// doing for very expensive position-only forces like atomic force fields. We
// try not to incur any overhead if there are no such forces in the System.
// Note in particular that Force::Gravity does its own caching so doesn't make
// use of this service.
class GeneralForceSubsystemRep : public ForceSubsystem::Guts {
public:
    GeneralForceSubsystemRep()
     : ForceSubsystemRep("GeneralForceSubsystem", "0.0.1")
    {
        calcForcesExecutor = new ParallelExecutor;
        calcForcesTask = new CalcForcesTask(*calcForcesExecutor);
    }
    
    ~GeneralForceSubsystemRep() {
        // Delete in reverse order to be nice to heap system.
        for (int i = (int)forces.size()-1; i >= 0; --i)
            delete forces[i]; 
        delete calcForcesExecutor;
        delete calcForcesTask;
    }
    
    ForceIndex adoptForce(Force& force) {
        invalidateSubsystemTopologyCache();
        const ForceIndex index((int) forces.size());
        forces.push_back(new Force()); // grow
        Force& f = *forces.back(); // refer to the empty handle we just created
        force.disown(f); // transfer ownership to f
        return index;
    }
    
    int getNumForces() const {
        return forces.size();
    }
    
    const Force& getForce(ForceIndex index) const {
        assert(index >= 0 && index < forces.size());
        return *forces[index];
    }

    Force& updForce(ForceIndex index) {
        assert(index >= 0 && index < forces.size());
        return *forces[index];
    }
    
    bool isForceDisabled(const State& state, ForceIndex index) const {
        const Array_<bool>& forceEnabled = Value< Array_<bool> >::downcast
            (getDiscreteVariable(state, forceEnabledIndex));
        return !forceEnabled[index];
    }
    
    void setForceIsDisabled
       (State& s, ForceIndex index, bool shouldDisable) const 
    {
        Array_<bool>& forceEnabled = Value< Array_<bool> >::updDowncast
            (updDiscreteVariable(s, forceEnabledIndex));

        const bool shouldEnable = !shouldDisable; // sorry
        if (forceEnabled[index] != shouldEnable) {
            forceEnabled[index] = shouldEnable;

            // If we're caching position-dependent forces, make sure they are
            // marked invalid here.
            if (cachedForcesAreValidCacheIndex.isValid()) {
                Value<bool>::updDowncast
                   (updCacheEntry(s, cachedForcesAreValidCacheIndex)) = false;
            }
        }
    }

    void setNumberOfThreads(unsigned int numThreads) {
        delete calcForcesExecutor;
        delete calcForcesTask;
        calcForcesExecutor = new ParallelExecutor(numThreads);
        calcForcesTask = new CalcForcesTask(*calcForcesExecutor);
    }

    // These override default implementations of virtual methods in the 
    // Subsystem::Guts class.

    GeneralForceSubsystemRep* cloneImpl() const override
    {   return new GeneralForceSubsystemRep(*this); }

    int realizeSubsystemTopologyImpl(State& s) const  override {
        forceEnabledIndex.invalidate();
        cachedForcesAreValidCacheIndex.invalidate();
        rigidBodyForceCacheIndex.invalidate();
        mobilityForceCacheIndex.invalidate();
        particleForceCacheIndex.invalidate();

        // Some forces are disabled by default; initialize the enabled flags
        // accordingly. Also, see if we're going to need to do any caching
        // on behalf of any forces that don't depend on velocities. 
        Array_<bool> forceEnabled(getNumForces());
        bool someForceElementNeedsCaching = false;
        for (int i = 0; i < (int)forces.size(); ++i) {
            forceEnabled[i] = !(forces[i]->isDisabledByDefault());
            if (!someForceElementNeedsCaching)
                someForceElementNeedsCaching = 
                    forces[i]->getImpl().dependsOnlyOnPositions();
        }

        forceEnabledIndex = allocateDiscreteVariable(s, Stage::Instance, 
            new Value<Array_<bool> >(forceEnabled));

        // Note that we'll allocate these even if all the needs-caching 
        // elements are presently disabled. That way they'll be around when
        // the force gets enabled.
        if (someForceElementNeedsCaching) {
            cachedForcesAreValidCacheIndex = 
                allocateCacheEntry(s, Stage::Position, new Value<bool>());
            rigidBodyForceCacheIndex = allocateCacheEntry(s, Stage::Dynamics, 
                new Value<Vector_<SpatialVec> >());
            mobilityForceCacheIndex = allocateCacheEntry(s, Stage::Dynamics, 
                new Value<Vector>());
            particleForceCacheIndex = allocateCacheEntry(s, Stage::Dynamics, 
                new Value<Vector_<Vec3> >());
        }

        // We must realizeTopology() even if the force is disabled by default.
        for (int i = 0; i < (int) forces.size(); ++i)
            forces[i]->getImpl().realizeTopology(s);
        return 0;
    }

    // Forces must realizeModel() even if they are currently disabled.
    int realizeSubsystemModelImpl(State& s) const override {
        for (int i = 0; i < (int) forces.size(); ++i)
            forces[i]->getImpl().realizeModel(s);
        return 0;
    }

    // No need to realize Instance stage or later for force elements that are
    // currently disabled.
    int realizeSubsystemInstanceImpl(const State& s) const override {
        const Array_<bool>& enabled = Value<Array_<bool> >::downcast
            (getDiscreteVariable(s, forceEnabledIndex));
        for (int i = 0; i < (int) forces.size(); ++i)
            if (enabled[i]) forces[i]->getImpl().realizeInstance(s);
        return 0;
    }

    int realizeSubsystemTimeImpl(const State& s) const override {
        const Array_<bool>& enabled = Value<Array_<bool> >::downcast
            (getDiscreteVariable(s, forceEnabledIndex));
        for (int i = 0; i < (int) forces.size(); ++i)
            if (enabled[i]) forces[i]->getImpl().realizeTime(s);
        return 0;
    }

    int realizeSubsystemPositionImpl(const State& s) const override {
        const Array_<bool>& enabled = Value<Array_<bool> >::downcast
            (getDiscreteVariable(s, forceEnabledIndex));
        // If we're caching position-dependent forces, make sure they are
        // marked invalid here.
        if (cachedForcesAreValidCacheIndex.isValid()) {
            Value<bool>::updDowncast
               (updCacheEntry(s, cachedForcesAreValidCacheIndex)) = false;
        }
        for (int i = 0; i < (int) forces.size(); ++i)
            if (enabled[i]) forces[i]->getImpl().realizePosition(s);
        return 0;
    }

    int realizeSubsystemVelocityImpl(const State& s) const override {
        const Array_<bool>& enabled = Value<Array_<bool> >::downcast
            (getDiscreteVariable(s, forceEnabledIndex));
        for (int i = 0; i < (int) forces.size(); ++i)
            if (enabled[i]) forces[i]->getImpl().realizeVelocity(s);
        return 0;
    }

    int realizeSubsystemDynamicsImpl(const State& s) const override {
        const MultibodySystem&        mbs    = getMultibodySystem(); // my owner
        const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();

        // Get access to the discrete state variable that records whether each
        // force element is enabled.
        const Array_<bool>& forceEnabled = Value< Array_<bool> >::downcast
                                    (getDiscreteVariable(s, forceEnabledIndex));

        // Get access to System-global force cache arrays.
        Vector_<SpatialVec>&   rigidBodyForces = 
                                    mbs.updRigidBodyForces(s, Stage::Dynamics);
        Vector_<Vec3>&         particleForces  = 
                                    mbs.updParticleForces (s, Stage::Dynamics);
        Vector&                mobilityForces  = 
                                    mbs.updMobilityForces (s, Stage::Dynamics);

        // Track the enabled forces and whether they should be parallelized.
        // TODO is this expensive?
        Array_<Force*> enabledNonParallelForces;
        Array_<Force*> enabledParallelForces;
        // Avoid repeatedly allocating memory.
        enabledNonParallelForces.reserve(forces.size());
        enabledParallelForces.reserve(forces.size());
        for (int i = 0; i < (int) forces.size(); ++i) {
            if (forceEnabled[i]) {
                if (forces[i]->getImpl().shouldBeParallelized())
                    enabledParallelForces.push_back(forces[i]);
                else
                    enabledNonParallelForces.push_back(forces[i]);
            }
        }

        // Short circuit if we're not doing any caching here. Note that we're
        // checking whether the *index* is valid (i.e. does the cache entry
        // exist?), not the contents.
        if (!cachedForcesAreValidCacheIndex.isValid()) {

            // Call calcForce() on all Forces, in parallel.
            calcForcesTask->calcForceAll(s,
                    enabledNonParallelForces, enabledParallelForces,
                    rigidBodyForces, particleForces, mobilityForces);
            /*
            for (int i = 0; i < (int)forces.size(); ++i) {
                if (forceEnabled[i])
                    forces[i]->getImpl().calcForce
                        (s, rigidBodyForces, particleForces, mobilityForces);
            }
            */

            // Allow forces to do their own realization, but wait until all
            // forces have executed calcForce(). TODO: not sure if that is
            // necessary (sherm 20130716).
            for (int i = 0; i < (int)forces.size(); ++i)
                if (forceEnabled[i]) 
                    forces[i]->getImpl().realizeDynamics(s);
            return 0;
        }

        // OK, we're doing some caching. This is a little messier. 

        // Get access to subsystem force cache entries.
        bool& cachedForcesAreValid = Value<bool>::downcast
                          (updCacheEntry(s, cachedForcesAreValidCacheIndex));

        Vector_<SpatialVec>&
            rigidBodyForceCache = Value<Vector_<SpatialVec> >::downcast
                                 (updCacheEntry(s, rigidBodyForceCacheIndex));
        Vector_<Vec3>&
            particleForceCache  = Value<Vector_<Vec3> >::downcast
                                 (updCacheEntry(s, particleForceCacheIndex));
        Vector&
            mobilityForceCache  = Value<Vector>::downcast
                                 (updCacheEntry(s, mobilityForceCacheIndex));


        if (!cachedForcesAreValid) {
            // We need to calculate the velocity independent forces.
            rigidBodyForceCache.resize(matter.getNumBodies());
            rigidBodyForceCache = SpatialVec(Vec3(0), Vec3(0));
            particleForceCache.resize(matter.getNumParticles());
            particleForceCache = Vec3(0);
            mobilityForceCache.resize(matter.getNumMobilities());
            mobilityForceCache = 0;

            // Run through all the forces, accumulating directly into the
            // force arrays or indirectly into the cache as appropriate.
            calcForcesTask->calcForceCachedAndNonCached(s,
                    enabledNonParallelForces, enabledParallelForces,
                    rigidBodyForces, particleForces, mobilityForces,
                    rigidBodyForceCache, particleForceCache,
                    mobilityForceCache);

            /* TODO leaving this here in case we allow not parallelizing.
            for (int i = 0; i < (int) forces.size(); ++i) {
                if (!forceEnabled[i]) continue;
                const ForceImpl& impl = forces[i]->getImpl();
                if (impl.dependsOnlyOnPositions())
                    impl.calcForce(s, rigidBodyForceCache, particleForceCache, 
                                      mobilityForceCache);
                else // ordinary velocity dependent force
                    impl.calcForce(s, rigidBodyForces, particleForces, 
                                      mobilityForces);
            }
            */

            cachedForcesAreValid = true;
        } else {

            // Cache already valid; just need to do the non-cached ones (the
            // ones for which dependsOnlyOnPositions is false).
            calcForcesTask->calcForceNonCached(s,
                    enabledNonParallelForces, enabledParallelForces,
                    rigidBodyForces, particleForces, mobilityForces);

            /* TODO leaving this here in case we allow not parallelizing.
            for (int i = 0; i < (int) forces.size(); ++i) {
                if (!forceEnabled[i]) continue;
                const ForceImpl& impl = forces[i]->getImpl();
                if (!impl.dependsOnlyOnPositions())
                    impl.calcForce(s, rigidBodyForces, particleForces, 
                                      mobilityForces);
            }
            */
        }

        // Accumulate the values from the cache into the global arrays.
        rigidBodyForces += rigidBodyForceCache;
        particleForces += particleForceCache;
        mobilityForces += mobilityForceCache;
        
        // Allow forces to do their own Dynamics-stage realization. Note that
        // this *follows* all the calcForce() calls.
        for (int i = 0; i < (int) forces.size(); ++i)
            if (forceEnabled[i]) forces[i]->getImpl().realizeDynamics(s);
        return 0;
    }
    
    Real calcPotentialEnergy(const State& state) const override {
        const Array_<bool>& forceEnabled = Value<Array_<bool> >::downcast
           (getDiscreteVariable(state, forceEnabledIndex)).get();
        Real energy = 0;
        for (int i = 0; i < (int) forces.size(); ++i) {
            if (forceEnabled[i]) {
                const Force& f = *forces[i];
                energy += f.getImpl().calcPotentialEnergy(state);
            }
        }
        return energy;
    }

    int realizeSubsystemAccelerationImpl(const State& s) const override {
        const Array_<bool>& enabled = Value<Array_<bool> >::downcast
            (getDiscreteVariable(s, forceEnabledIndex));
        for (int i = 0; i < (int) forces.size(); ++i)
            if (enabled[i]) forces[i]->getImpl().realizeAcceleration(s);
        return 0;
    }

    int realizeSubsystemReportImpl(const State& s) const override {
        const Array_<bool>& enabled = Value<Array_<bool> >::downcast
            (getDiscreteVariable(s, forceEnabledIndex));
        for (int i = 0; i < (int) forces.size(); ++i)
            if (enabled[i]) forces[i]->getImpl().realizeReport(s);
        return 0;
    }

    int calcDecorativeGeometryAndAppendImpl
       (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const
       override
    {
        const Array_<bool>& enabled = Value<Array_<bool> >::downcast
            (getDiscreteVariable(s, forceEnabledIndex));
        for (int i = 0; i < (int) forces.size(); ++i)
            if (enabled[i]) 
                forces[i]->getImpl().calcDecorativeGeometryAndAppend(s,stage,geom);
        return 0;
    }

private:
    Array_<Force*>                  forces;

    // For parallel calculation of forces.
    // TODO figure out memory management.
    ParallelExecutor* calcForcesExecutor;
    CalcForcesTask* calcForcesTask;

        // TOPOLOGY "CACHE"
    // These indices must be filled in during realizeTopology and treated
    // as const thereafter.

    // This instance-stage variable holds a bool for each force element.
    mutable DiscreteVariableIndex   forceEnabledIndex;

    // This set of cache entries is allocated only if some force element
    // overrode dependsOnlyOnPositions().
    mutable CacheEntryIndex         cachedForcesAreValidCacheIndex;
    mutable CacheEntryIndex         rigidBodyForceCacheIndex;
    mutable CacheEntryIndex         mobilityForceCacheIndex;
    mutable CacheEntryIndex         particleForceCacheIndex;
};

    ///////////////////////////
    // GeneralForceSubsystem //
    ///////////////////////////


/*static*/ bool 
GeneralForceSubsystem::isInstanceOf(const ForceSubsystem& s) {
    return GeneralForceSubsystemRep::isA(s.getRep());
}
/*static*/ const GeneralForceSubsystem&
GeneralForceSubsystem::downcast(const ForceSubsystem& s) {
    assert(isInstanceOf(s));
    return static_cast<const GeneralForceSubsystem&>(s);
}
/*static*/ GeneralForceSubsystem&
GeneralForceSubsystem::updDowncast(ForceSubsystem& s) {
    assert(isInstanceOf(s));
    return static_cast<GeneralForceSubsystem&>(s);
}

const GeneralForceSubsystemRep& 
GeneralForceSubsystem::getRep() const {
    return SimTK_DYNAMIC_CAST_DEBUG<const GeneralForceSubsystemRep&>(ForceSubsystem::getRep());
}
GeneralForceSubsystemRep&       
GeneralForceSubsystem::updRep() {
    return SimTK_DYNAMIC_CAST_DEBUG<GeneralForceSubsystemRep&>(ForceSubsystem::updRep());
}

// Create Subsystem but don't associate it with any System. This isn't much use except
// for making std::vector's, which require a default constructor to be available.
GeneralForceSubsystem::GeneralForceSubsystem()
:   ForceSubsystem() 
{   adoptSubsystemGuts(new GeneralForceSubsystemRep()); }

GeneralForceSubsystem::GeneralForceSubsystem(MultibodySystem& mbs)
:   ForceSubsystem() 
{   adoptSubsystemGuts(new GeneralForceSubsystemRep());
    mbs.addForceSubsystem(*this); } // steal ownership

ForceIndex GeneralForceSubsystem::adoptForce(Force& force) 
{   return updRep().adoptForce(force); }

int GeneralForceSubsystem::getNumForces() const 
{   return getRep().getNumForces(); }

const Force& GeneralForceSubsystem::getForce(ForceIndex index) const 
{   return getRep().getForce(index); }

Force& GeneralForceSubsystem::updForce(ForceIndex index) 
{   return updRep().updForce(index); }

bool GeneralForceSubsystem::isForceDisabled
   (const State& state, ForceIndex index) const 
{   return getRep().isForceDisabled(state, index); }

void GeneralForceSubsystem::setForceIsDisabled
   (State& state, ForceIndex index, bool disabled) const 
{   getRep().setForceIsDisabled(state, index, disabled); }

void GeneralForceSubsystem::setNumberOfThreads(unsigned int numThreads)
{   updRep().setNumberOfThreads(numThreads); }

const MultibodySystem& GeneralForceSubsystem::getMultibodySystem() const
{   return MultibodySystem::downcast(getSystem()); }

} // namespace SimTK

