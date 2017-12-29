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
 * Contributors: Peter Eastman, Nabeel Allana, Chris Dembia, Thomas Lau       *
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
#include <iostream>
#include <exception>

#include "ForceImpl.h"

#include <memory>

//Threading constants used by CalcForcesTask
namespace {
using namespace SimTK;

const int NumNonParallelThreads = 1;
const int NonParallelForcesIndex = 0;

/* Base class for CalcForcesParallelTask and CalcForcesNonParallelTask - lays 
out common methods that will be implemented to suit the parallel/non-parallel
use cases*/
class CalcForcesTask : public ParallelExecutor::Task {
public:
    CalcForcesTask() = default;
    
    virtual CalcForcesTask* clone() const = 0;
    
    virtual void initializeAll(
            const Array_<Force*>& forces, const State& s,
            const Array_<ForceIndex>& enabledNonParallelForces,
            const Array_<ForceIndex>& enabledParallelForces,
            Vector_<SpatialVec>& rigidBodyForces,
            Vector_<Vec3>& particleForces,
            Vector& mobilityForces) = 0;
    virtual void initializeCachedAndNonCached(
            const Array_<Force*>& forces, const State& s,
            const Array_<ForceIndex>& enabledNonParallelForces,
            const Array_<ForceIndex>& enabledParallelForces,
            Vector_<SpatialVec>& rigidBodyForces,
            Vector_<Vec3>& particleForces,
            Vector& mobilityForces,
            Vector_<SpatialVec>& rigidBodyForceCache,
            Vector_<Vec3>& particleForceCache,
            Vector& mobilityForceCache) = 0;
    virtual void initializeNonCached(
            const Array_<Force*>& forces, const State& s,
            const Array_<ForceIndex>& enabledNonParallelForces,
            const Array_<ForceIndex>& enabledParallelForces,
            Vector_<SpatialVec>& rigidBodyForces,
            Vector_<Vec3>& particleForces,
            Vector& mobilityForces) = 0;
    
};
/*Calculates each enabled force's contribution in the MultibodySystem.
CalcForcesParallelTask allows force calculations to occur in parallel with
non-parallel forces being calculated on Thread NonParallelForcesIndex and all
other parallel forces being calculated on separate threads.*/

//Implementation of CalcForcesTask for parallel forces
class CalcForcesParallelTask : public CalcForcesTask {
public:
    //The different categories of force calculations that GeneralForceSubsystem
    //chooses between
    enum Mode {
        All,
        CachedAndNonCached,
        NonCached
    };
    
    CalcForcesParallelTask* clone() const override {
        return new CalcForcesParallelTask();
    }
    // Each different force calculation Mode requires different parameters
    //
    // Note: Execute MUST be called directly after CalcForceTask is initialized
    void initializeAll(
            const Array_<Force*>& forces, const State& s,
            const Array_<ForceIndex>& enabledNonParallelForces,
            const Array_<ForceIndex>& enabledParallelForces,
            Vector_<SpatialVec>& rigidBodyForces,
            Vector_<Vec3>& particleForces,
            Vector& mobilityForces) override
    {
        m_forces = &forces;
        m_state = &s;
        m_enabledNonParallelForces = &enabledNonParallelForces;
        m_enabledParallelForces = &enabledParallelForces;
        m_rigidBodyForces = &rigidBodyForces;
        m_particleForces = &particleForces;
        m_mobilityForces = &mobilityForces;
        
        m_mode = All;
    }
    void initializeCachedAndNonCached(
            const Array_<Force*>& forces, const State& s,
            const Array_<ForceIndex>& enabledNonParallelForces,
            const Array_<ForceIndex>& enabledParallelForces,
            Vector_<SpatialVec>& rigidBodyForces,
            Vector_<Vec3>& particleForces,
            Vector& mobilityForces,
            Vector_<SpatialVec>& rigidBodyForceCache,
            Vector_<Vec3>& particleForceCache,
            Vector& mobilityForceCache) override
    {
        m_forces = &forces;
        m_state = s;
        m_enabledNonParallelForces = &enabledNonParallelForces;
        m_enabledParallelForces = &enabledParallelForces;
        m_rigidBodyForces = &rigidBodyForces;
        m_particleForces = &particleForces;
        m_mobilityForces = &mobilityForces;
        m_rigidBodyForceCache = &rigidBodyForceCache;
        m_particleForceCache = &particleForceCache;
        m_mobilityForceCache = &mobilityForceCache;

        m_mode = CachedAndNonCached;
    }
    void initializeNonCached(
            const Array_<Force*>& forces, const State& s,
            const Array_<ForceIndex>& enabledNonParallelForces,
            const Array_<ForceIndex>& enabledParallelForces,
            Vector_<SpatialVec>& rigidBodyForces,
            Vector_<Vec3>& particleForces,
            Vector& mobilityForces) override
    {
        m_forces = &forces;
        m_state = s;
        m_enabledNonParallelForces = &enabledNonParallelForces;
        m_enabledParallelForces = &enabledParallelForces;
        m_rigidBodyForces = &rigidBodyForces;
        m_particleForces = &particleForces;
        m_mobilityForces = &mobilityForces;

        m_mode = NonCached;
    }
    
    // Set all of the local force contribution arrays to 0. This is so that
    // each local thread can sum up its force contribution to be later
    // added into the total force array.
    void initialize() override {
        m_rigidBodyForcesLocalStatic.resize(m_rigidBodyForces->size());
        m_rigidBodyForcesLocalStatic.setToZero();
        m_particleForcesLocalStatic.resize(m_particleForces->size());
        m_particleForcesLocalStatic.setToZero();
        m_mobilityForcesLocalStatic.resize(m_mobilityForces->size());
        m_mobilityForcesLocalStatic.setToZero();

        if (m_mode == CachedAndNonCached) {
            m_rigidBodyForceCacheLocalStatic.resize(m_rigidBodyForceCache->size());
            m_rigidBodyForceCacheLocalStatic.setToZero();
            m_particleForceCacheLocalStatic.resize(m_particleForceCache->size());
            m_particleForceCacheLocalStatic.setToZero();
            m_mobilityForceCacheLocalStatic.resize(m_mobilityForceCache->size());
            m_mobilityForceCacheLocalStatic.setToZero();
        }
    }
    
    // Calculate all enabled forces (taking into account mode and parallelism)
    void execute(int threadIndex) override {
      switch (m_mode) {
        case All:
            if (threadIndex == NonParallelForcesIndex) {
                // Process all non-parallel forces
                for (const auto& forceIndex : *m_enabledNonParallelForces) {
                    const auto force = m_forces.getRef()[forceIndex];
                    force->getImpl().calcForce(*m_state, m_rigidBodyForcesLocalStatic, m_particleForcesLocalStatic, m_mobilityForcesLocalStatic);
                }
            } else {
                // Process a single parallel force. Subtract 1 from index b/c
                // we use 0 for the non-parallel forces.
                const auto& forceIndex =
                        m_enabledParallelForces->getElt(threadIndex-1);
                const auto& impl = m_forces.getRef()[forceIndex]->getImpl();
                impl.calcForce(*m_state, m_rigidBodyForcesLocalStatic, m_particleForcesLocalStatic, m_mobilityForcesLocalStatic);

            }
            break;

        case CachedAndNonCached:
            if (threadIndex == NonParallelForcesIndex) {
                // Process all non-parallel forces.
                for (const auto& forceIndex : *m_enabledNonParallelForces) {
                    const auto& impl = m_forces.getRef()[forceIndex]->getImpl();
                    if (impl.dependsOnlyOnPositions()) {
                        impl.calcForce(*m_state, *m_rigidBodyForceCache, *m_particleForceCache, *m_mobilityForceCache);
                    } else { // ordinary velocity dependent force
                        impl.calcForce(*m_state, *m_rigidBodyForces, *m_particleForces, *m_mobilityForces);
                    }
                }
            } else {
                // Process a single parallel force. Subtract 1 from threadIndex b/c
                // we use 0 for the non-parallel forces.
                const auto& forceIndex =
                        m_enabledParallelForces->getElt(threadIndex-1);
                const auto& impl = m_forces.getRef()[forceIndex]->getImpl();
                if (impl.dependsOnlyOnPositions()) {
                    impl.calcForce(*m_state, m_rigidBodyForceCacheLocalStatic, m_particleForceCacheLocalStatic, m_mobilityForceCacheLocalStatic);
                } else { // ordinary velocity dependent force
                    impl.calcForce(*m_state, m_rigidBodyForcesLocalStatic, m_particleForcesLocalStatic, m_mobilityForcesLocalStatic);
                }
            }
            break;

        case NonCached:
            if (threadIndex == NonParallelForcesIndex) {
                // Process all non-parallel forces.
                for (const auto& forceIndex : *m_enabledNonParallelForces) {
                    const auto& impl = m_forces.getRef()[forceIndex]->getImpl();
                    if (!impl.dependsOnlyOnPositions()) {
                        impl.calcForce(*m_state,
                                *m_rigidBodyForces, *m_particleForces,
                                *m_mobilityForces);
                    }
                }
            } else {
                // Process a single parallel force. Subtract 1 from index b/c
                // we use 0 for the non-parallel forces.
                const auto& forceIndex =
                        m_enabledParallelForces->getElt(threadIndex-1);
                const auto& impl = m_forces.getRef()[forceIndex]->getImpl();
                if (!impl.dependsOnlyOnPositions()) {
                    impl.calcForce(*m_state,
                            m_rigidBodyForcesLocalStatic, m_particleForcesLocalStatic,
                            m_mobilityForcesLocalStatic);
                }
            }
            break;
        }
    }
    
    //Once a thread has finished its force calculations, we add in the thread's
    //contribution into the cached force arrays in the State
    void finish() override {
        *m_rigidBodyForces += m_rigidBodyForcesLocalStatic;
        *m_particleForces += m_particleForcesLocalStatic;
        *m_mobilityForces += m_mobilityForcesLocalStatic;

        if (m_mode == CachedAndNonCached) {
            *m_rigidBodyForceCache += m_rigidBodyForceCacheLocalStatic;
            *m_particleForceCache += m_particleForceCacheLocalStatic;
            *m_mobilityForceCache += m_mobilityForceCacheLocalStatic;
        }
    }
private:
    Mode m_mode;

    ReferencePtr<const Array_<Force*>> m_forces;
    ReferencePtr<const State> m_state;

    // Constant state-cache variables for the enabled parallel and non-parallel
    // forces.
    ReferencePtr<const Array_<ForceIndex>> m_enabledNonParallelForces;
    ReferencePtr<const Array_<ForceIndex>> m_enabledParallelForces;

    // ReferencePtrs that point to caches in the state; We eventually add our
    // thread-local result to these vectors.
    ReferencePtr<Vector_<SpatialVec>> m_rigidBodyForces;
    ReferencePtr<Vector_<Vec3>> m_particleForces;
    ReferencePtr<Vector> m_mobilityForces;

    ReferencePtr<Vector_<SpatialVec>> m_rigidBodyForceCache;
    ReferencePtr<Vector_<Vec3>> m_particleForceCache;
    ReferencePtr<Vector> m_mobilityForceCache;
    
    // These variables are local to a thread. They are set to their default
    // value when the threads are spawned. We use them to keep track of each
    // thread's contribution that we will later add in to the final state cache.
    static thread_local Vector_<SpatialVec> m_rigidBodyForcesLocalStatic;
    static thread_local Vector_<Vec3> m_particleForcesLocalStatic;
    static thread_local Vector m_mobilityForcesLocalStatic;

    static thread_local Vector_<SpatialVec> m_rigidBodyForceCacheLocalStatic;
    static thread_local Vector_<Vec3> m_particleForceCacheLocalStatic;
    static thread_local Vector m_mobilityForceCacheLocalStatic;
};

//local declarations of static member variables
/*static*/ thread_local Vector_<SpatialVec>
                    CalcForcesParallelTask::m_rigidBodyForcesLocalStatic;
/*static*/ thread_local Vector_<Vec3>
                    CalcForcesParallelTask::m_particleForcesLocalStatic;
/*static*/ thread_local Vector
                    CalcForcesParallelTask::m_mobilityForcesLocalStatic;

/*static*/ thread_local Vector_<SpatialVec>
                    CalcForcesParallelTask::m_rigidBodyForceCacheLocalStatic;
/*static*/ thread_local Vector_<Vec3>
                    CalcForcesParallelTask::m_particleForceCacheLocalStatic;
/*static*/ thread_local Vector
                    CalcForcesParallelTask::m_mobilityForceCacheLocalStatic;

/* Calculates each enabled force's contribution in the MultibodySystem. These
calculations occur on the main thread, without use of local thread variables.*/

//Implementation of CalcForcesTask for non-parallel forces
class CalcForcesNonParallelTask : public CalcForcesTask {
public:
    //The different categories of force calculations that GeneralForceSubsystem
    //chooses between
    enum Mode {
        All,
        CachedAndNonCached,
        NonCached
    };
    
    CalcForcesNonParallelTask* clone() const override {
        return new CalcForcesNonParallelTask();
    }
    
    // Each different force calculation Mode requires different parameters
    //
    // Note: Execute MUST be called directly after CalcForceTask is initialized
    void initializeAll(
            const Array_<Force*>& forces, const State& s,
            const Array_<ForceIndex>& enabledNonParallelForces,
            const Array_<ForceIndex>& enabledParallelForces,
            Vector_<SpatialVec>& rigidBodyForces,
            Vector_<Vec3>& particleForces,
            Vector& mobilityForces) override
    {
        m_forces = &forces;
        m_state = &s;
        m_enabledNonParallelForces = &enabledNonParallelForces;
        m_enabledParallelForces = &enabledParallelForces;
        m_rigidBodyForces = &rigidBodyForces;
        m_particleForces = &particleForces;
        m_mobilityForces = &mobilityForces;
        
        m_mode = All;
    }
    void initializeCachedAndNonCached(
            const Array_<Force*>& forces, const State& s,
            const Array_<ForceIndex>& enabledNonParallelForces,
            const Array_<ForceIndex>& enabledParallelForces,
            Vector_<SpatialVec>& rigidBodyForces,
            Vector_<Vec3>& particleForces,
            Vector& mobilityForces,
            Vector_<SpatialVec>& rigidBodyForceCache,
            Vector_<Vec3>& particleForceCache,
            Vector& mobilityForceCache) override
    {
        m_forces = &forces;
        m_state = s;
        m_enabledNonParallelForces = &enabledNonParallelForces;
        m_enabledParallelForces = &enabledParallelForces;
        m_rigidBodyForces = &rigidBodyForces;
        m_particleForces = &particleForces;
        m_mobilityForces = &mobilityForces;
        m_rigidBodyForceCache = &rigidBodyForceCache;
        m_particleForceCache = &particleForceCache;
        m_mobilityForceCache = &mobilityForceCache;

        m_mode = CachedAndNonCached;
    }
    void initializeNonCached(
            const Array_<Force*>& forces, const State& s,
            const Array_<ForceIndex>& enabledNonParallelForces,
            const Array_<ForceIndex>& enabledParallelForces,
            Vector_<SpatialVec>& rigidBodyForces,
            Vector_<Vec3>& particleForces,
            Vector& mobilityForces) override
    {
        m_forces = &forces;
        m_state = s;
        m_enabledNonParallelForces = &enabledNonParallelForces;
        m_enabledParallelForces = &enabledParallelForces;
        m_rigidBodyForces = &rigidBodyForces;
        m_particleForces = &particleForces;
        m_mobilityForces = &mobilityForces;

        m_mode = NonCached;
    }
    
    // Set all of the local force contribution arrays to 0. This is so that
    // each local thread can sum up it's force contribution to be later
    // added into the total force array.
    void initialize() override {
        m_rigidBodyForcesLocal.resize(m_rigidBodyForces->size());
        m_rigidBodyForcesLocal.setToZero();
        m_particleForcesLocal.resize(m_particleForces->size());
        m_particleForcesLocal.setToZero();
        m_mobilityForcesLocal.resize(m_mobilityForces->size());
        m_mobilityForcesLocal.setToZero();

        if (m_mode == CachedAndNonCached) {
            m_rigidBodyForceCacheLocal.resize(m_rigidBodyForceCache->size());
            m_rigidBodyForceCacheLocal.setToZero();
            m_particleForceCacheLocal.resize(m_particleForceCache->size());
            m_particleForceCacheLocal.setToZero();
            m_mobilityForceCacheLocal.resize(m_mobilityForceCache->size());
            m_mobilityForceCacheLocal.setToZero();
        }
    }
    
    // Calculate all enabled forces (taking into account mode and parallelism)
    void execute(int threadIndex) override {
      switch (m_mode) {
        case All:
            if (threadIndex == NonParallelForcesIndex) {
                // Process all non-parallel forces
                for (const auto& forceIndex : *m_enabledNonParallelForces) {
                    const auto force = m_forces.getRef()[forceIndex];
                    force->getImpl().calcForce(*m_state, m_rigidBodyForcesLocal,
                                  m_particleForcesLocal, m_mobilityForcesLocal);
                }
            }
            break;

        case CachedAndNonCached:
            if (threadIndex == NonParallelForcesIndex) {
                // Process all non-parallel forces.
                for (const auto& forceIndex : *m_enabledNonParallelForces) {
                    const auto& impl = m_forces.getRef()[forceIndex]->getImpl();
                    if (impl.dependsOnlyOnPositions()) {
                        impl.calcForce(*m_state, *m_rigidBodyForceCache,
                                  *m_particleForceCache, *m_mobilityForceCache);
                    } else { // ordinary velocity dependent force
                        impl.calcForce(*m_state, *m_rigidBodyForces,
                                          *m_particleForces, *m_mobilityForces);
                    }
                }
            }
            break;

        case NonCached:
            if (threadIndex == NonParallelForcesIndex) {
                // Process all non-parallel forces.
                for (const auto& forceIndex : *m_enabledNonParallelForces) {
                    const auto& impl = m_forces.getRef()[forceIndex]->getImpl();
                    if (!impl.dependsOnlyOnPositions()) {
                        impl.calcForce(*m_state,
                                *m_rigidBodyForces, *m_particleForces,
                                *m_mobilityForces);
                    }
                }
            }
            break;
        }
    }
    
    //Once a thread has finished it's force calculations, we add in the thread's
    //contribution into the cached force arrays in the State
    void finish() override {
        *m_rigidBodyForces += m_rigidBodyForcesLocal;
        *m_particleForces += m_particleForcesLocal;
        *m_mobilityForces += m_mobilityForcesLocal;

        if (m_mode == CachedAndNonCached) {
            *m_rigidBodyForceCache += m_rigidBodyForceCacheLocal;
            *m_particleForceCache += m_particleForceCacheLocal;
            *m_mobilityForceCache += m_mobilityForceCacheLocal;
        }
    }
private:
    Mode m_mode;

    ReferencePtr<const Array_<Force*>> m_forces;
    ReferencePtr<const State> m_state;

    // Const. state-cache variables for the enabled parallel and non-parallel
    // forces.
    ReferencePtr<const Array_<ForceIndex>> m_enabledNonParallelForces;
    ReferencePtr<const Array_<ForceIndex>> m_enabledParallelForces;

    // ReferencePtrs that point to caches in the state; We eventually add our
    // thread-local result to these vectors.
    ReferencePtr<Vector_<SpatialVec>> m_rigidBodyForces;
    ReferencePtr<Vector_<Vec3>> m_particleForces;
    ReferencePtr<Vector> m_mobilityForces;

    ReferencePtr<Vector_<SpatialVec>> m_rigidBodyForceCache;
    ReferencePtr<Vector_<Vec3>> m_particleForceCache;
    ReferencePtr<Vector> m_mobilityForceCache;
    
    // These variables belong to the main thread. We use them to keep track of
    // each thread's contribution that we will later add in to the final state
    // cache.
    Vector_<SpatialVec> m_rigidBodyForcesLocal;
    Vector_<Vec3> m_particleForcesLocal;
    Vector m_mobilityForcesLocal;

    Vector_<SpatialVec> m_rigidBodyForceCacheLocal;
    Vector_<Vec3> m_particleForceCacheLocal;
    Vector m_mobilityForceCacheLocal;
};
} //namespace

namespace SimTK{
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
        //The default number of threads is the physical number of processors
        //call setNumberOfThreads() if you want to override the thread count
        calcForcesExecutor = new ParallelExecutor();
    }

    ~GeneralForceSubsystemRep() {
        // Delete in reverse order to be nice to heap system.
        for (int i = (int) forces.size()-1; i >= 0; --i)
            delete forces[i];
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

    void setNumberOfThreads(unsigned numThreads) {
        SimTK_APIARGCHECK_ALWAYS(numThreads > 0, "GeneralForceSubsystemRep",
                    "setNumberOfThreads", "Number of threads must be positive");
        calcForcesExecutor = new ParallelExecutor(numThreads);
    }
    
    int getNumberOfThreads() const{
      return calcForcesExecutor->getMaxThreads();
    }

    // These override default implementations of virtual methods in the
    // Subsystem::Guts class.

    GeneralForceSubsystemRep* cloneImpl() const override
    {   return new GeneralForceSubsystemRep(*this); }

    int realizeSubsystemTopologyImpl(State& s) const  override {
        forceEnabledIndex.invalidate();
        enabledParallelForcesIndex.invalidate();
        enabledNonParallelForcesIndex.invalidate();
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
        
        
        // Track the enabled forces and whether they should be parallelized.
        Array_<ForceIndex> enabledNonParallelForces;
        Array_<ForceIndex> enabledParallelForces;

        // Avoid repeatedly allocating memory.
        enabledNonParallelForces.reserve(forces.size());
        enabledParallelForces.reserve(forces.size());

        for (int i = 0; i < (int) forces.size(); ++i) {
            if (forceEnabled[i])
            {
                if (forces[i]->getImpl().shouldBeParallelIfPossible())
                    enabledParallelForces.push_back(ForceIndex(i));
                else
                    enabledNonParallelForces.push_back(ForceIndex(i));
            }
        }

        enabledNonParallelForcesIndex = allocateCacheEntry(s, Stage::Instance,
                new Value<Array_<ForceIndex> >(enabledNonParallelForces));
        enabledParallelForcesIndex = allocateCacheEntry(s, Stage::Instance,
                new Value<Array_<ForceIndex> >(enabledParallelForces));

        //Determine whether the subsystem has parallel forces - if so, use the
        //parallel implementation of CalcForcesTask (even if those parallel
        //forces are not currently enabled - they could be enabled in the
        //future)
        //
        // *Consider the following example:*
        // Parallel Forces: A B C
        // NonParallel Forces: D E
        //
        // Enabled Forces at Time 0: D E
        // Enabled Forces at Time 1: A B C D E
        
        bool hasParallelForces = false;
        for(int x = 0; x < (int)forces.size(); ++x)
        {
            if (forces[x]->getImpl().shouldBeParallelIfPossible())
            {
                hasParallelForces = true;
                break;
            }
        }
        if (hasParallelForces)
            calcForcesTask = new CalcForcesParallelTask();
        else
            calcForcesTask = new CalcForcesNonParallelTask();
        
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
            
        //Update non/parallel enabled force arrays
        // Track the enabled forces and whether they should be parallelized.
        const Array_<bool>& forceEnabled = Value< Array_<bool> >::downcast
                                    (getDiscreteVariable(s, forceEnabledIndex));
        Array_<ForceIndex>& enabledNonParallelForces =
                Value< Array_<ForceIndex> >::
                    updDowncast(updCacheEntry(s,enabledNonParallelForcesIndex));
        Array_<ForceIndex>& enabledParallelForces =
                Value< Array_<ForceIndex> >::
                    updDowncast(updCacheEntry(s,enabledParallelForcesIndex));

        // Avoid repeatedly allocating memory.
        enabledParallelForces.resize(0);
        enabledNonParallelForces.resize(0);
        
        enabledNonParallelForces.reserve(forces.size());
        enabledParallelForces.reserve(forces.size());

        for (int i = 0; i < (int) forces.size(); ++i) {
            if (forceEnabled[i])
            {
                if (forces[i]->getImpl().shouldBeParallelIfPossible())
                    enabledParallelForces.push_back(ForceIndex(i));
                else
                    enabledNonParallelForces.push_back(ForceIndex(i));
            }
        }
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
        const MultibodySystem&        mbs    = getMultibodySystem();
        const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();

        // Get access to the discrete state variable that records whether each
        // force element is enabled.
        const Array_<bool>& forceEnabled = Value< Array_<bool> >::downcast
                                    (getDiscreteVariable(s, forceEnabledIndex));
        const Array_<ForceIndex>& enabledNonParallelForces =
                Value<Array_<ForceIndex>>::
                      downcast(getCacheEntry(s, enabledNonParallelForcesIndex));
        const Array_<ForceIndex>& enabledParallelForces =
                Value<Array_<ForceIndex>>::
                         downcast(getCacheEntry(s, enabledParallelForcesIndex));

        // Get access to System-global force cache arrays.
        Vector_<SpatialVec>&   rigidBodyForces =
                                    mbs.updRigidBodyForces(s, Stage::Dynamics);
        Vector_<Vec3>&         particleForces  =
                                    mbs.updParticleForces (s, Stage::Dynamics);
        Vector&                mobilityForces  =
                                    mbs.updMobilityForces (s, Stage::Dynamics);

        // Short circuit if we're not doing any caching here. Note that we're
        // checking whether the *index* is valid (i.e. does the cache entry
        // exist?), not the contents.
        if (!cachedForcesAreValidCacheIndex.isValid()) {
            // Call calcForce() on all Forces, in parallel.
            calcForcesTask->initializeAll(forces, s,
                    enabledNonParallelForces, enabledParallelForces,
                    rigidBodyForces, particleForces, mobilityForces);
            calcForcesExecutor->execute(calcForcesTask.updRef(),
                          enabledParallelForces.size() + NumNonParallelThreads);

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
        bool& cachedForcesAreValid = Value<bool>::updDowncast
                          (updCacheEntry(s, cachedForcesAreValidCacheIndex));

        Vector_<SpatialVec>&
            rigidBodyForceCache = Value<Vector_<SpatialVec> >::updDowncast
                                 (updCacheEntry(s, rigidBodyForceCacheIndex));
        Vector_<Vec3>&
            particleForceCache  = Value<Vector_<Vec3> >::updDowncast
                                 (updCacheEntry(s, particleForceCacheIndex));
        Vector&
            mobilityForceCache  = Value<Vector>::updDowncast
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
            calcForcesTask->initializeCachedAndNonCached(forces, s,
                                enabledNonParallelForces, enabledParallelForces,
                                rigidBodyForces, particleForces, mobilityForces,
                                rigidBodyForceCache, particleForceCache,
                                mobilityForceCache);
            calcForcesExecutor->execute(calcForcesTask.updRef(),
                          enabledParallelForces.size() + NumNonParallelThreads);
            cachedForcesAreValid = true;
        } else {
            // Cache already valid; just need to do the non-cached ones (the
            // ones for which dependsOnlyOnPositions is false).
            calcForcesTask->initializeNonCached(forces, s,
                               enabledNonParallelForces, enabledParallelForces,
                               rigidBodyForces, particleForces, mobilityForces);
            calcForcesExecutor->execute(calcForcesTask.updRef(),
                          enabledParallelForces.size() + NumNonParallelThreads);
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
    mutable ClonePtr<ParallelExecutor>               calcForcesExecutor;
    mutable ClonePtr<CalcForcesTask>                 calcForcesTask;
    
    // TOPOLOGY "CACHE"
    // These indices must be filled in during realizeTopology and treated
    // as const thereafter.

    // This instance-stage variable holds a bool for each force element.
    mutable DiscreteVariableIndex   forceEnabledIndex;
    
    //This set of cache entries stores an array of Force* elements for enabled
    //parallel and non-parallel forces
    mutable CacheEntryIndex   enabledParallelForcesIndex;
    mutable CacheEntryIndex   enabledNonParallelForcesIndex;

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

void GeneralForceSubsystem::setNumberOfThreads(unsigned numThreads)
{   updRep().setNumberOfThreads(numThreads); }

int GeneralForceSubsystem::getNumberOfThreads() const
{   return getRep().getNumberOfThreads(); }

const MultibodySystem& GeneralForceSubsystem::getMultibodySystem() const
{   return MultibodySystem::downcast(getSystem()); }

} // namespace SimTK
