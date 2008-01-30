/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-7 Stanford University and the Authors.         *
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

#include "simbody/internal/common.h"
#include "simbody/internal/UniformGravitySubsystem.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"

#include "ForceSubsystemRep.h"

#include <cmath>

namespace SimTK {

// 
// Define a uniform gravity field that affects all the matter in the system.
// Parameters exist for the gravity vector, zero height, and enable/disable.
// The MultibodySystem we're part of provides the memory into which we
// accumulate forces and potential energy.
//

class UniformGravitySubsystemRep : public ForceSubsystemRep {
public:
    UniformGravitySubsystemRep() 
      : ForceSubsystemRep("UniformGravitySubsystem", "0.0.1"), 
        instanceVarsIndex(-1), instanceCacheIndex(-1) { }

    explicit UniformGravitySubsystemRep(const Vec3& g, const Real& z=0)
      : ForceSubsystemRep("UniformGravitySubsystem", "0.0.1"), 
        defaultParameters(g,z),
        instanceVarsIndex(-1), instanceCacheIndex(-1) { } 

    // Access to state variables (parameters).
    const Vec3& getGravity(const State& s) const {return getParameters(s).gravity;}
    Vec3&       updGravity(State& s)       const {return updParameters(s).gravity;}

    const Real& getZeroHeight(const State& s) const {return getParameters(s).zeroHeight;}
    Real&       updZeroHeight(State& s)       const {return updParameters(s).zeroHeight;}

    bool  isEnabled(const State& s) const {return getParameters(s).enabled;}
    bool& updIsEnabled(State& s)    const {return updParameters(s).enabled;}

    // Responses (not available through the client-side handle class).
    const Real& getGravityMagnitude(const State& s) const {
        return getInstanceCache(s).gMagnitude;
    }
    const Real& getPEOffset(const State& s) const {
        return getInstanceCache(s).gz;
    }


    UniformGravitySubsystemRep* cloneImpl() const {return new UniformGravitySubsystemRep(*this);}

    int realizeSubsystemTopologyImpl(State& s) const;
    //   realizeSubsystemModelImpl() not needed
    int realizeSubsystemInstanceImpl(const State& s) const;
    //   realizeSubsystemTime, Position, VelocityImpl not needed
    int realizeSubsystemDynamicsImpl(const State& s) const;
    //   realizeSubsystemAccelerationImpl() not needed
    //   realizeSubsystemReportImpl() not needed


private:
    // State entries. TODO: these should be at a later stage
    // since they can't affect anything until potential energy.
    struct Parameters {
        Parameters() 
          : gravity(0), zeroHeight(0), enabled(true) { }
        Parameters(const Vec3& g, const Real& z) 
          : gravity(g), zeroHeight(z), enabled(true) { }

        Vec3 gravity;
        Real zeroHeight;
        bool enabled;
    };

    struct ParameterCache {
        Real gMagnitude;    // |g|, used to avoid work if g=0
        Real gz;            // precalculated PE offset: pe = m*(|g|h - |g|z)
    };

    // topological variables
    Parameters defaultParameters;

    // These must be filled in during realizeTopology and treated
    // as const thereafter.
    mutable int instanceVarsIndex;
    mutable int instanceCacheIndex;

    const Parameters& getParameters(const State& s) const {
        assert(subsystemTopologyHasBeenRealized());
        return Value<Parameters>::downcast(
            getDiscreteVariable(s,instanceVarsIndex)).get();
    }
    Parameters& updParameters(State& s) const {
        assert(subsystemTopologyHasBeenRealized());
        return Value<Parameters>::downcast(
            updDiscreteVariable(s,instanceVarsIndex)).upd();
    }

    const ParameterCache& getInstanceCache(const State& s) const {
        assert(subsystemTopologyHasBeenRealized());
        return Value<ParameterCache>::downcast(
            getCacheEntry(s,instanceCacheIndex)).get();
    }
    ParameterCache& updInstanceCache(const State& s) const {
        assert(subsystemTopologyHasBeenRealized());
        return Value<ParameterCache>::downcast(
            updCacheEntry(s,instanceCacheIndex)).upd();
    }

    friend std::ostream& operator<<(std::ostream& o, 
                         const UniformGravitySubsystemRep::Parameters&); 
    friend std::ostream& operator<<(std::ostream& o, 
                         const UniformGravitySubsystemRep::ParameterCache&);
};

// Useless, but required by Value<T>.
std::ostream& operator<<(std::ostream& o, 
                         const UniformGravitySubsystemRep::Parameters&) 
{assert(false);return o;}
std::ostream& operator<<(std::ostream& o, 
                         const UniformGravitySubsystemRep::ParameterCache&) 
{assert(false);return o;}


    ///////////////////////////////
    // UNIFORM GRAVITY SUBSYSTEM //
    ///////////////////////////////

/*static*/ bool 
UniformGravitySubsystem::isInstanceOf(const ForceSubsystem& s) {
    return UniformGravitySubsystemRep::isA(s.getRep());
}
/*static*/ const UniformGravitySubsystem&
UniformGravitySubsystem::downcast(const ForceSubsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const UniformGravitySubsystem&>(s);
}
/*static*/ UniformGravitySubsystem&
UniformGravitySubsystem::updDowncast(ForceSubsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<UniformGravitySubsystem&>(s);
}

const UniformGravitySubsystemRep& 
UniformGravitySubsystem::getRep() const {
    return dynamic_cast<const UniformGravitySubsystemRep&>(ForceSubsystem::getRep());
}
UniformGravitySubsystemRep&       
UniformGravitySubsystem::updRep() {
    return dynamic_cast<UniformGravitySubsystemRep&>(ForceSubsystem::updRep());
}

// Create Subsystem but don't associate it with any System. This isn't much use except
// for making std::vector's, which require a default constructor to be available.
UniformGravitySubsystem::UniformGravitySubsystem()
  : ForceSubsystem()
{
     adoptSubsystemGuts(new UniformGravitySubsystemRep());
}

UniformGravitySubsystem::UniformGravitySubsystem(MultibodySystem& mbs)
  : ForceSubsystem() 
{
    adoptSubsystemGuts(new UniformGravitySubsystemRep());
    mbs.addForceSubsystem(*this); // steal ownership
}

UniformGravitySubsystem::UniformGravitySubsystem(MultibodySystem& mbs, const Vec3& g, const Real& zeroHeight)
  : ForceSubsystem()
{
    adoptSubsystemGuts(new UniformGravitySubsystemRep(g,zeroHeight));
    mbs.addForceSubsystem(*this); // steal ownership
}

const Vec3& UniformGravitySubsystem::getGravity(const State& s) const {
    return getRep().getGravity(s);
}
Vec3& UniformGravitySubsystem::updGravity(State& s) const {
    return getRep().updGravity(s);
}
const Real& UniformGravitySubsystem::getZeroHeight(const State& s) const {
    return getRep().getZeroHeight(s);
}
Real& UniformGravitySubsystem::updZeroHeight(State& s) const {
    return getRep().updZeroHeight(s);
}
bool UniformGravitySubsystem::isEnabled(const State& s) const {
    return getRep().isEnabled(s);
}
bool& UniformGravitySubsystem::updIsEnabled(State& s) const {
    return getRep().updIsEnabled(s);
}

    ///////////////////////////////////
    // UNIFORM GRAVITY SUBSYSTEM REP //
    ///////////////////////////////////

int UniformGravitySubsystemRep::realizeSubsystemTopologyImpl(State& s) const {
    // Note that although these are *instance* variables, they are allocated as
    // part of the *topology*. That allows us to store the indices locally rather
    // than in the state.
    instanceVarsIndex = allocateDiscreteVariable(s, Stage::Instance, 
        new Value<Parameters>(defaultParameters));
    instanceCacheIndex = allocateCacheEntry(s, Stage::Instance,
        new Value<ParameterCache>());
    return 0;
}

// realizeModel() not needed since there are no modeling options here

int UniformGravitySubsystemRep::realizeSubsystemInstanceImpl(const State& s) const {
    // any values are acceptable
    ParameterCache& pc = updInstanceCache(s);
    pc.gMagnitude = getGravity(s).norm();
    pc.gz = pc.gMagnitude * getZeroHeight(s);
    return 0;
}

// realizeTime, Position, Velocity not needed

int UniformGravitySubsystemRep::realizeSubsystemDynamicsImpl(const State& s) const {
    if (!isEnabled(s) || getGravityMagnitude(s)==0)
        return 0; // nothing to do

    const Vec3& g   = getGravity(s);  // gravity is non zero
    const Real& gz  = getPEOffset(s); // amount to subtract from gh for pe

    const MultibodySystem&        mbs    = MultibodySystem::downcast(getSystem());
    const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();
    const int nBodies    = matter.getNBodies();
    const int nParticles = matter.getNParticles();

    Vector_<SpatialVec>& rigidBodyForces = mbs.updRigidBodyForces(s, Stage::Dynamics);
    Vector_<Vec3>&       particleForces  = mbs.updParticleForces (s, Stage::Dynamics);
    Real&                pe              = mbs.updPotentialEnergy(s, Stage::Dynamics);

    assert(rigidBodyForces.size() == nBodies);
    assert(particleForces.size() == nParticles);

    if (nParticles) {
        const Vector& m = matter.getAllParticleMasses(s);
        const Vector_<Vec3>& loc_G = matter.getAllParticleLocations(s);
        for (int i=0; i < nParticles; ++i) {
            pe -= m[i]*(~g*loc_G[i] + gz); // odd signs because height is in -g direction
            particleForces[i] += g * m[i];
        }
    }

    // no need to apply gravity to Ground!
    for (MobilizedBodyIndex i(1); i < nBodies; ++i) {
        const MassProperties& mprops = matter.getMobilizedBody(i).getBodyMassProperties(s);
        const Real&      m       = mprops.getMass();
        const Vec3&      com_B   = mprops.getMassCenter();
        const Transform& X_GB    = matter.getMobilizedBody(i).getBodyTransform(s);
        const Vec3       com_B_G = X_GB.R()*com_B;
        const Vec3       com_G   = X_GB.T() + com_B_G;
        const Vec3       frc_G   = m*g;

        pe -= m*(~g*com_G + gz); // odd signs because height is in -g direction
        rigidBodyForces[i] += SpatialVec(com_B_G % frc_G, frc_G); 
    }

    return 0;
}

} // namespace SimTK

