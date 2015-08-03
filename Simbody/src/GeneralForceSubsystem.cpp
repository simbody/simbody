/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-14 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors: Peter Eastman, Nabeel Allana                                 *
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


namespace SimTK {

// There is some tricky caching being done here for forces that have overridden
// dependsOnlyOnPositions() (and returned "true"). This is probably only worth
// doing for very expensive position-only forces like atomic force fields. We
// try not to incur any overhead if there are no such forces in the System.
// Note in particular that Force::Gravity does its own caching so doesn't make
// use of this service.
class GeneralForceSubsystemRep : public ForceSubsystem::Guts {
public:
	GeneralForceSubsystemRep(int nThreads)
		: ForceSubsystemRep("GeneralForceSubsystem", "0.0.1")
    {
		assert(nThreads > 0);

		work_queue = new ParallelWorkQueue11(nThreads);

		rigidBodyForcesPool = new ObjectPool<Vector_<SpatialVec>>();
		particleForcesPool = new ObjectPool<Vector_<Vec3>>;
		mobilityForcesPool = new ObjectPool<Vector>();

		rigidBodyForcesCachePool = new ObjectPool<Vector_<SpatialVec>>();
		particleForcesCachePool = new ObjectPool<Vector_<Vec3>>;
		mobilityForcesCachePool = new ObjectPool<Vector>();

		rigidBodyForcesPool->initialize(nThreads);
		particleForcesPool->initialize(nThreads);
		mobilityForcesPool->initialize(nThreads);

		rigidBodyForcesCachePool->initialize(nThreads);
		particleForcesCachePool->initialize(nThreads);
		mobilityForcesCachePool->initialize(nThreads);

    }
    
    ~GeneralForceSubsystemRep() {
        // Delete in reverse order to be nice to heap system.
        for (int i = (int)forces.size()-1; i >= 0; --i)
            delete forces[i]; 

		delete work_queue;
		delete rigidBodyForcesPool;
		delete particleForcesPool;
		delete mobilityForcesPool;
		delete rigidBodyForcesCachePool;
		delete particleForcesCachePool;
		delete mobilityForcesCachePool;
    }

	void setNumberOfThreads(int nThreads) const
	{
		assert(nThreads > 0);

		work_queue->resize(nThreads);

		rigidBodyForcesPool->resize(nThreads);
		particleForcesPool->resize(nThreads);
		mobilityForcesPool->resize(nThreads);

		rigidBodyForcesCachePool->resize(nThreads);
		particleForcesCachePool->resize(nThreads);
		mobilityForcesCachePool->resize(nThreads);
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

		// zero out all entries of the force pools
		for (Vector_<SpatialVec>* vec : rigidBodyForcesPool->list())
		{
			vec->resize(matter.getNumBodies());
			*vec = SpatialVec(Vec3(0), Vec3(0));
		}
		for (Vector_<Vec3>* vec : particleForcesPool->list())
		{
			vec->resize(matter.getNumParticles());
			*vec = Vec3(0);
		}
		for (Vector* vec : mobilityForcesPool->list())
		{
			vec->resize(matter.getNumMobilities());
			*vec = 0;
		}

        // Short circuit if we're not doing any caching here. Note that we're
        // checking whether the *index* is valid (i.e. does the cache entry
        // exist?), not the contents.
        if (!cachedForcesAreValidCacheIndex.isValid()) {

			// push all forces that requested parallelism into a seperate entry
			// within the work queue, creating a lambda function to execute it
            for (int i = 0; i < (int)forces.size(); ++i) {
				if (forceEnabled[i])
				if (forces[i]->getImpl().shouldBeParallelized())
					work_queue->push(
					[this, &s, i]
				{
							Vector_<SpatialVec>* rigidBodyForce = this->rigidBodyForcesPool->get();
							Vector_<Vec3>* particleForce = this->particleForcesPool->get();
							Vector* mobilityForce = this->mobilityForcesPool->get();

							this->forces[i]->getImpl().calcForce
								(s, *rigidBodyForce, *particleForce, *mobilityForce);

							this->rigidBodyForcesPool->push(rigidBodyForce);
							this->particleForcesPool->push(particleForce);
							this->mobilityForcesPool->push(mobilityForce);
						}
					);
            }

			// push an entry to the queue containing all other force elements
			work_queue->push([this, &forceEnabled, &s] {
				
				Vector_<SpatialVec>* rigidBodyForce = this->rigidBodyForcesPool->get();
				Vector_<Vec3>* particleForce = this->particleForcesPool->get();
				Vector* mobilityForce = this->mobilityForcesPool->get();

				for (int i = 0; i < (int)this->forces.size(); ++i)
				if (forceEnabled[i])
				if (!forces[i]->getImpl().shouldBeParallelized())
					forces[i]->getImpl().calcForce
					(s, *rigidBodyForce, *particleForce, *mobilityForce);

				this->rigidBodyForcesPool->push(rigidBodyForce);
				this->particleForcesPool->push(particleForce);
				this->mobilityForcesPool->push(mobilityForce);

			});

			work_queue->execute();

			// sum up resulting forces from each pool
			for (Vector_<SpatialVec>* vec : rigidBodyForcesPool->list())
				rigidBodyForces += *vec;
			for (Vector_<Vec3>* vec : particleForcesPool->list())
				particleForces += *vec;
			for (Vector* vec : mobilityForcesPool->list())
				mobilityForces += *vec;

            // Allow forces to do their own realization, but wait until all
            // forces have executed calcForce(). TODO: not sure if that is
            // necessary (sherm 20130716).
            for (int i = 0; i < (int)forces.size(); ++i)
            if (forceEnabled[i]) 
			if (forces[i]->getImpl().shouldBeParallelized())
				work_queue->push([this, &s, i]
				{
					this->forces[i]->getImpl().realizeDynamics(s);
				}
				);
			
			work_queue->push([this, &forceEnabled, &s] {
						
				for (int i = 0; i < (int)this->forces.size(); ++i)
				if (forceEnabled[i])
				if (!forces[i]->getImpl().shouldBeParallelized())
					forces[i]->getImpl().realizeDynamics(s);
					
			});

			work_queue->execute();

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

			for (Vector_<SpatialVec>* vec : rigidBodyForcesCachePool->list())
			{
				vec->resize(matter.getNumBodies());
				*vec = SpatialVec(Vec3(0), Vec3(0));
			}
			for (Vector_<Vec3>* vec : particleForcesCachePool->list())
			{
				vec->resize(matter.getNumParticles());
				*vec = Vec3(0);
			}
			for (Vector* vec : mobilityForcesCachePool->list())
			{
				vec->resize(matter.getNumMobilities());
				*vec = 0;
			}

            rigidBodyForceCache.resize(matter.getNumBodies());
            rigidBodyForceCache = SpatialVec(Vec3(0), Vec3(0));
            particleForceCache.resize(matter.getNumParticles());
            particleForceCache = Vec3(0);
            mobilityForceCache.resize(matter.getNumMobilities());
            mobilityForceCache = 0;


			for (int i = 0; i < (int)forces.size(); ++i) {
				if (forceEnabled[i])
				if (forces[i]->getImpl().shouldBeParallelized())
					work_queue->push(
					[this, &s, i]
				{
					const ForceImpl& impl = forces[i]->getImpl();
					
					Vector_<SpatialVec>* rigidBodyForce;
					Vector_<Vec3>* particleForce;
					Vector* mobilityForce;

					if (impl.dependsOnlyOnPositions())
					{
						rigidBodyForce = this->rigidBodyForcesCachePool->get();
						particleForce = this->particleForcesCachePool->get();
						mobilityForce = this->mobilityForcesCachePool->get();
					}
					else
					{
						rigidBodyForce = this->rigidBodyForcesPool->get();
						particleForce = this->particleForcesPool->get();
						mobilityForce = this->mobilityForcesPool->get();
					}

					impl.calcForce
						(s, *rigidBodyForce, *particleForce, *mobilityForce);

					if (impl.dependsOnlyOnPositions())
					{
						this->rigidBodyForcesCachePool->push(rigidBodyForce);
						this->particleForcesCachePool->push(particleForce);
						this->mobilityForcesCachePool->push(mobilityForce);
					}
					else
					{
						this->rigidBodyForcesPool->push(rigidBodyForce);
						this->particleForcesPool->push(particleForce);
						this->mobilityForcesPool->push(mobilityForce);
					}
				}
				);
			}

			work_queue->push([this, &forceEnabled, &s] {

				Vector_<SpatialVec>* rigidBodyForce = this->rigidBodyForcesPool->get();
				Vector_<Vec3>* particleForce = this->particleForcesPool->get();
				Vector* mobilityForce = this->mobilityForcesPool->get();

				Vector_<SpatialVec>* rigidBodyCacheForce = this->rigidBodyForcesCachePool->get();
				Vector_<Vec3>* particleCacheForce = this->particleForcesCachePool->get();
				Vector* mobilityCacheForce = this->mobilityForcesCachePool->get();

				for (int i = 0; i < (int)this->forces.size(); ++i)
				if (forceEnabled[i])
				if (!forces[i]->getImpl().shouldBeParallelized())
				{
					const ForceImpl& impl = forces[i]->getImpl();
					if (impl.dependsOnlyOnPositions())
						impl.calcForce(s, *rigidBodyCacheForce, *particleCacheForce,
						*mobilityCacheForce);
					else // ordinary velocity dependent force
						impl.calcForce(s, *rigidBodyForce, *particleForce,
						*mobilityForce);
				}

				this->rigidBodyForcesPool->push(rigidBodyForce);
				this->particleForcesPool->push(particleForce);
				this->mobilityForcesPool->push(mobilityForce);

				this->rigidBodyForcesCachePool->push(rigidBodyCacheForce);
				this->particleForcesCachePool->push(particleCacheForce);
				this->mobilityForcesCachePool->push(mobilityCacheForce);

			});

			work_queue->execute();

			for (Vector_<SpatialVec>* vec : rigidBodyForcesCachePool->list())
				rigidBodyForceCache += *vec;
			for (Vector_<Vec3>* vec : particleForcesCachePool->list())
				particleForceCache += *vec;
			for (Vector* vec : mobilityForcesCachePool->list())
				mobilityForceCache += *vec;

			for (Vector_<SpatialVec>* vec : rigidBodyForcesPool->list())
				rigidBodyForces += *vec;
			for (Vector_<Vec3>* vec : particleForcesPool->list())
				particleForces += *vec;
			for (Vector* vec : mobilityForcesPool->list())
				mobilityForces += *vec;

			cachedForcesAreValid = true;
        } else {

			for (int i = 0; i < (int)forces.size(); ++i) {
				if (forceEnabled[i])
				if (!forces[i]->getImpl().dependsOnlyOnPositions())
				if (forces[i]->getImpl().shouldBeParallelized())
					work_queue->push(
					[this, &s, i]
				{
					Vector_<SpatialVec>* rigidBodyForce = this->rigidBodyForcesPool->get();
					Vector_<Vec3>* particleForce = this->particleForcesPool->get();
					Vector* mobilityForce = this->mobilityForcesPool->get();

					this->forces[i]->getImpl().calcForce
						(s, *rigidBodyForce, *particleForce, *mobilityForce);

					this->rigidBodyForcesPool->push(rigidBodyForce);
					this->particleForcesPool->push(particleForce);
					this->mobilityForcesPool->push(mobilityForce);
				}
				);
			}

			work_queue->push([this, &forceEnabled, &s] {

				Vector_<SpatialVec>* rigidBodyForce = this->rigidBodyForcesPool->get();
				Vector_<Vec3>* particleForce = this->particleForcesPool->get();
				Vector* mobilityForce = this->mobilityForcesPool->get();

				for (int i = 0; i < (int)this->forces.size(); ++i)
				if (forceEnabled[i])
				if (!forces[i]->getImpl().dependsOnlyOnPositions())
				if (!forces[i]->getImpl().shouldBeParallelized())
					forces[i]->getImpl().calcForce
					(s, *rigidBodyForce, *particleForce, *mobilityForce);

				this->rigidBodyForcesPool->push(rigidBodyForce);
				this->particleForcesPool->push(particleForce);
				this->mobilityForcesPool->push(mobilityForce);

			});

			work_queue->execute();

			for (Vector_<SpatialVec>* vec : rigidBodyForcesPool->list())
				rigidBodyForces += *vec;
			for (Vector_<Vec3>* vec : particleForcesPool->list())
				particleForces += *vec;
			for (Vector* vec : mobilityForcesPool->list())
				mobilityForces += *vec;

        }

        // Accumulate the values from the cache into the global arrays.
        rigidBodyForces += rigidBodyForceCache;
        particleForces += particleForceCache;
        mobilityForces += mobilityForceCache;
        
        // Allow forces to do their own Dynamics-stage realization. Note that
        // this *follows* all the calcForce() calls.

		for (int i = 0; i < (int)forces.size(); ++i)
		if (forceEnabled[i])
		if (forces[i]->getImpl().shouldBeParallelized())
			work_queue->push(
				[this, &s, i]
				{
					forces[i]->getImpl().realizeDynamics(s);
				}
			);

		work_queue->push([this, &forceEnabled, &s] {

			for (int i = 0; i < (int)this->forces.size(); ++i)
			if (forceEnabled[i])
			if (!forces[i]->getImpl().shouldBeParallelized())
				forces[i]->getImpl().realizeDynamics(s);

		});

		work_queue->execute();

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

    Array_<Force*>						forces;
	ParallelWorkQueue11*				work_queue;

	ObjectPool<Vector_<SpatialVec>>*	rigidBodyForcesPool;
	ObjectPool<Vector_<Vec3>>*			particleForcesPool;
	ObjectPool<Vector>*					mobilityForcesPool;

	ObjectPool<Vector_<SpatialVec>>*	rigidBodyForcesCachePool;
	ObjectPool<Vector_<Vec3>>*			particleForcesCachePool;
	ObjectPool<Vector>*					mobilityForcesCachePool;
    
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
{
	adoptSubsystemGuts(new GeneralForceSubsystemRep(ParallelExecutor::getNumProcessors()));
}

GeneralForceSubsystem::GeneralForceSubsystem(MultibodySystem& mbs)
:   ForceSubsystem() 
{
	adoptSubsystemGuts(new GeneralForceSubsystemRep(ParallelExecutor::getNumProcessors()));
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


void GeneralForceSubsystem::setNumberOfThreads(int nThreads) const
{ getRep().setNumberOfThreads(nThreads); }

const MultibodySystem& GeneralForceSubsystem::getMultibodySystem() const
{   return MultibodySystem::downcast(getSystem()); }

} // namespace SimTK

