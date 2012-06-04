/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-10 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors: Peter Eastman                                                *
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


class GeneralForceSubsystemRep : public ForceSubsystemRep {
public:
    GeneralForceSubsystemRep()
     : ForceSubsystemRep("GeneralForceSubsystem", "0.0.1")
    {
    }
    
    ~GeneralForceSubsystemRep() {
        // Delete in reverse order to be nice to heap system.
        for (int i = (int)forces.size()-1; i >= 0; --i)
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
       (State& state, ForceIndex index, bool disabled) const {
        Array_<bool>& forceEnabled = Value< Array_<bool> >::downcast
            (updDiscreteVariable(state, forceEnabledIndex));
        bool& forceValid = Value<bool>::downcast
            (updCacheEntry(state, forceValidCacheIndex));
        forceEnabled[index] = !disabled;
        forceValid = false;
    }

    
    // These override default implementations of virtual methods in the Subsystem::Guts
    // class.

    GeneralForceSubsystemRep* cloneImpl() const 
    {   return new GeneralForceSubsystemRep(*this); }

    int realizeSubsystemTopologyImpl(State& s) const {
        forceValidCacheIndex = allocateCacheEntry(s, Stage::Position, 
            new Value<bool>());
        rigidBodyForceCacheIndex = allocateCacheEntry(s, Stage::Dynamics, 
            new Value<Vector_<SpatialVec> >());
        mobilityForceCacheIndex = allocateCacheEntry(s, Stage::Dynamics, 
            new Value<Vector>());
        particleForceCacheIndex = allocateCacheEntry(s, Stage::Dynamics, 
            new Value<Vector_<Vec3> >());

        // Some forces are disabled by default; initialize the enabled flags
        // accordingly.
        Array_<bool> forceEnabled(getNumForces());
        for (int i = 0; i < (int) forces.size(); ++i)
            forceEnabled[i] = !(forces[i]->isDisabledByDefault());

        forceEnabledIndex = allocateDiscreteVariable(s, Stage::Instance, 
            new Value<Array_<bool> >(forceEnabled));

        // Note that we must realizeTopology() even if the force is disabled
        // by default.
        for (int i = 0; i < (int) forces.size(); ++i)
            forces[i]->getImpl().realizeTopology(s);

        return 0;
    }

    // Forces must realizeModel() even if they are currently disabled.
    int realizeSubsystemModelImpl(State& s) const {
        for (int i = 0; i < (int) forces.size(); ++i)
            forces[i]->getImpl().realizeModel(s);
        return 0;
    }

    // No need to realize Instance stage or later for force elements that are
    // currently disabled.
    int realizeSubsystemInstanceImpl(const State& s) const {
        const Array_<bool>& enabled = Value<Array_<bool> >::downcast
            (getDiscreteVariable(s, forceEnabledIndex));
        for (int i = 0; i < (int) forces.size(); ++i)
            if (enabled[i]) forces[i]->getImpl().realizeInstance(s);
        return 0;
    }

    int realizeSubsystemTimeImpl(const State& s) const {
        const Array_<bool>& enabled = Value<Array_<bool> >::downcast
            (getDiscreteVariable(s, forceEnabledIndex));
        for (int i = 0; i < (int) forces.size(); ++i)
            if (enabled[i]) forces[i]->getImpl().realizeTime(s);
        return 0;
    }

    int realizeSubsystemPositionImpl(const State& s) const {
        const Array_<bool>& enabled = Value<Array_<bool> >::downcast
            (getDiscreteVariable(s, forceEnabledIndex));
        Value<bool>::downcast(updCacheEntry(s, forceValidCacheIndex)) = false;
        for (int i = 0; i < (int) forces.size(); ++i)
            if (enabled[i]) forces[i]->getImpl().realizePosition(s);
        return 0;
    }

    int realizeSubsystemVelocityImpl(const State& s) const {
        const Array_<bool>& enabled = Value<Array_<bool> >::downcast
            (getDiscreteVariable(s, forceEnabledIndex));
        for (int i = 0; i < (int) forces.size(); ++i)
            if (enabled[i]) forces[i]->getImpl().realizeVelocity(s);
        return 0;
    }

    int realizeSubsystemDynamicsImpl(const State& s) const {

        const MultibodySystem&        mbs    = getMultibodySystem(); // my owner
        const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();

        // Get access to subsystem force cache entries.
        Vector_<SpatialVec>&    
            rigidBodyForceCache = Value<Vector_<SpatialVec> >::downcast
                                 (updCacheEntry(s, rigidBodyForceCacheIndex));
        Vector_<Vec3>&         
            particleForceCache  = Value<Vector_<Vec3> >::downcast
                                 (updCacheEntry(s, particleForceCacheIndex));
        Vector&                 
            mobilityForceCache  = Value<Vector>::downcast
                                 (updCacheEntry(s, mobilityForceCacheIndex));

        // Get access to subsystem discrete state and "is valid" cache.
        const Array_<bool>& 
            forceEnabled = Value<Array_<bool> >::downcast
                          (getDiscreteVariable(s, forceEnabledIndex));
        bool& forceValid = Value<bool>::downcast
                          (updCacheEntry(s, forceValidCacheIndex));

        if (!forceValid) {
            // We need to calculate the velocity independent forces.
            rigidBodyForceCache.resize(matter.getNumBodies());
            rigidBodyForceCache = SpatialVec(Vec3(0), Vec3(0));
            particleForceCache.resize(matter.getNumParticles());
            particleForceCache = Vec3(0);
            mobilityForceCache.resize(matter.getNumMobilities());
            mobilityForceCache = 0;
        }

        // Get access to System-global force cache arrays.
        Vector_<SpatialVec>&   rigidBodyForces = mbs.updRigidBodyForces(s, Stage::Dynamics);
        Vector_<Vec3>&         particleForces  = mbs.updParticleForces (s, Stage::Dynamics);
        Vector&                mobilityForces  = mbs.updMobilityForces (s, Stage::Dynamics);
        for (int i = 0; i < (int) forces.size(); ++i) {
            if (forceEnabled[i]) {
                const Force& f = *forces[i];
                if (!f.getImpl().dependsOnlyOnPositions())
                    f.getImpl().calcForce(s, rigidBodyForces, particleForces, mobilityForces);
                else if (!forceValid)
                    f.getImpl().calcForce(s, rigidBodyForceCache, particleForceCache, mobilityForceCache);
            }
        }

        // Accumulate the values from the cache into the global arrays.
        rigidBodyForces += rigidBodyForceCache;
        particleForces += particleForceCache;
        mobilityForces += mobilityForceCache;
        forceValid = true;
        
        // Allow forces to do their own realization.
        // TODO: sherm 100313: shouldn't this come before calcForce() calls?
        for (int i = 0; i < (int) forces.size(); ++i)
            if (forceEnabled[i]) forces[i]->getImpl().realizeDynamics(s);

        return 0;
    }
    
    Real calcPotentialEnergy(const State& state) const {
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

    int realizeSubsystemAccelerationImpl(const State& s) const {
        const Array_<bool>& enabled = Value<Array_<bool> >::downcast
            (getDiscreteVariable(s, forceEnabledIndex));
        for (int i = 0; i < (int) forces.size(); ++i)
            if (enabled[i]) forces[i]->getImpl().realizeAcceleration(s);
        return 0;
    }

    int realizeSubsystemReportImpl(const State& s) const {
        const Array_<bool>& enabled = Value<Array_<bool> >::downcast
            (getDiscreteVariable(s, forceEnabledIndex));
        for (int i = 0; i < (int) forces.size(); ++i)
            if (enabled[i]) forces[i]->getImpl().realizeReport(s);
        return 0;
    }

private:
    Array_<Force*>                  forces;
    
        // TOPOLOGY "CACHE"
    // This must be filled in during realizeTopology and treated
    // as const thereafter.
    mutable CacheEntryIndex         forceValidCacheIndex;
    mutable CacheEntryIndex         rigidBodyForceCacheIndex;
    mutable CacheEntryIndex         mobilityForceCacheIndex;
    mutable CacheEntryIndex         particleForceCacheIndex;
    mutable DiscreteVariableIndex   forceEnabledIndex;
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
    return reinterpret_cast<const GeneralForceSubsystem&>(s);
}
/*static*/ GeneralForceSubsystem&
GeneralForceSubsystem::updDowncast(ForceSubsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<GeneralForceSubsystem&>(s);
}

const GeneralForceSubsystemRep& 
GeneralForceSubsystem::getRep() const {
    return dynamic_cast<const GeneralForceSubsystemRep&>(ForceSubsystem::getRep());
}
GeneralForceSubsystemRep&       
GeneralForceSubsystem::updRep() {
    return dynamic_cast<GeneralForceSubsystemRep&>(ForceSubsystem::updRep());
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

const MultibodySystem& GeneralForceSubsystem::getMultibodySystem() const
{   return MultibodySystem::downcast(getSystem()); }

} // namespace SimTK

