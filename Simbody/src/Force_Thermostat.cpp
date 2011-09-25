/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Michael Sherman                                                   *
 * Contributors: Chris Bruns                                                  *
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

#include "SimTKcommon.h"

#include "simbody/internal/common.h"
#include "simbody/internal/MobilizedBody.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/Force_Thermostat.h"

#include "ForceImpl.h"

namespace SimTK {

// Implementation class for Force::Thermostat.
class Force::ThermostatImpl : public ForceImpl {
public:
    ThermostatImpl(const SimbodyMatterSubsystem& matter, 
                   Real boltzmannsConstant, 
                   Real defBathTemp, 
                   Real defRelaxationTime,
                   int  defNumExcludedDofs)
    :   matter(matter), kB(boltzmannsConstant), 
        defaultNumChains(DefaultDefaultNumChains), 
        defaultBathTemp(defBathTemp),
        defaultRelaxationTime(defRelaxationTime), 
        defaultNumExcludedDofs(defNumExcludedDofs) {}

    ThermostatImpl* clone() const {return new ThermostatImpl(*this);}
    bool dependsOnlyOnPositions() const {return false;}

    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, 
                   Vector_<Vec3>& particleForces, Vector& mobilityForces) const;

    // Temperature does not contribute to potential energy.
    Real calcPotentialEnergy(const State& state) const {return 0;}

    void realizeTopology(State& state) const;
    void realizeModel(State& state) const;
    void realizeVelocity(const State& state) const;
    void realizeDynamics(const State& state) const;

    // Get/update the current number of chains.
    int getNumChains(const State& s) const {
        assert(dvNumChains.isValid());
        return Value<int>::downcast(getForceSubsystem().getDiscreteVariable(s, dvNumChains));
    }
    int& updNumChains(State& s) const {
        assert(dvNumChains.isValid());
        return Value<int>::updDowncast(getForceSubsystem().updDiscreteVariable(s, dvNumChains));
    }

    // Get/update the current number of excluded dofs.
    int getNumExcludedDofs(const State& s) const {
        assert(dvNumExcludedDofs.isValid());
        return Value<int>::downcast
            (getForceSubsystem().getDiscreteVariable(s, dvNumExcludedDofs));
    }
    int& updNumExcludedDofs(State& s) const {
        assert(dvNumExcludedDofs.isValid());
        return Value<int>::updDowncast
            (getForceSubsystem().updDiscreteVariable(s, dvNumExcludedDofs));
    }

    int getNumThermalDOFs(const State& s) const;

    // Get the auxiliary continuous state index of the 0'th thermostat state variable z.
    ZIndex getZ0Index(const State& s) const {
        assert(cacheZ0Index.isValid());
        return Value<ZIndex>::downcast(getForceSubsystem().getCacheEntry(s, cacheZ0Index));
    }
    ZIndex& updZ0Index(State& s) const {
        assert(cacheZ0Index.isValid());
        return Value<ZIndex>::updDowncast(getForceSubsystem().updCacheEntry(s, cacheZ0Index));
    }

    // Get the current bath temperature (i.e., the desired temperature).
    Real getBathTemp(const State& s) const {
        assert(dvBathTemp.isValid());
        return Value<Real>::downcast(getForceSubsystem().getDiscreteVariable(s, dvBathTemp));
    }
    Real& updBathTemp(State& s) const {
        assert(dvBathTemp.isValid());
        return Value<Real>::updDowncast(getForceSubsystem().updDiscreteVariable(s, dvBathTemp));
    }

    // Get the current bath temperature (i.e., the desired temperature).
    Real getRelaxationTime(const State& s) const {
        assert(dvRelaxationTime.isValid());
        return Value<Real>::downcast(getForceSubsystem().getDiscreteVariable(s, dvRelaxationTime));
    }
    Real& updRelaxationTime(State& s) const {
        assert(dvRelaxationTime.isValid());
        return Value<Real>::updDowncast(getForceSubsystem().updDiscreteVariable(s, dvRelaxationTime));
    }

    // Get the calculated momentum M*u (after Stage::Velocity).
    const Vector& getMomentum(const State& s) const {
        assert(cacheMomentumIndex.isValid());
        return Value<Vector>::downcast(getForceSubsystem().getCacheEntry(s, cacheMomentumIndex));
    }
    Vector& updMomentum(const State& s) const {
        assert(cacheMomentumIndex.isValid());
        return Value<Vector>::updDowncast(getForceSubsystem().updCacheEntry(s, cacheMomentumIndex));
    }

    // Get the calculated system kinetic energy ~u*M*u/2 (after Stage::Velocity).
    const Real& getKE(const State& s) const {
        assert(cacheKEIndex.isValid());
        return Value<Real>::downcast(getForceSubsystem().getCacheEntry(s, cacheKEIndex));
    }
    Real& updKE(const State& s) const {
        assert(cacheKEIndex.isValid());
        return Value<Real>::updDowncast(getForceSubsystem().updCacheEntry(s, cacheKEIndex));
    }

    Real getExternalWork(const State& s) const 
    {   return getForceSubsystem().getZ(s)[workZIndex]; }
    Real& updExternalWork(State& s) const 
    {   return getForceSubsystem().updZ(s)[workZIndex]; }

    Real& updWorkZDot(const State& s) const 
    {   return getForceSubsystem().updZDot(s)[workZIndex];  }

    Real calcExternalPower(const State& s) const;

    // Get the current value of one of the thermostat state variables.
    Real getZ(const State& s, int i) const {
        assert(0 <= i && i < 2*getNumChains(s));
        const ZIndex z0 = getZ0Index(s);
        return getForceSubsystem().getZ(s)[z0+i];
    }
    Real& updZ(State& s, int i) const {
        assert(0 <= i && i < 2*getNumChains(s));
        const ZIndex z0 = getZ0Index(s);
        return getForceSubsystem().updZ(s)[z0+i];
    }
    Real getZDot(const State& s, int i) const {
        assert(0 <= i && i < 2*getNumChains(s));
        const ZIndex z0 = getZ0Index(s);
        return getForceSubsystem().getZDot(s)[z0+i];
    }
    // State is const here because ZDot is a cache entry.
    Real& updZDot(const State& s, int i) const {
        assert(0 <= i && i < 2*getNumChains(s));
        const ZIndex z0 = getZ0Index(s);
        return getForceSubsystem().updZDot(s)[z0+i];
    }

    static const int DefaultDefaultNumChains = 3;
private:
    const SimbodyMatterSubsystem& matter;
    const Real  kB;     // Boltzmann's constant in compatible units

    // Topology-stage "state" variables.
    int         defaultNumChains;       // # chains in a new State
    int         defaultNumExcludedDofs; // # non-thermal rigid body dofs
    Real        defaultBathTemp;        // bath temperature
    Real        defaultRelaxationTime;  // relaxation time

    // These indices are Topology-stage "cache" variables.
    DiscreteVariableIndex dvNumChains;          // integer
    DiscreteVariableIndex dvNumExcludedDofs;    // integer
    DiscreteVariableIndex dvBathTemp;           // Real
    DiscreteVariableIndex dvRelaxationTime;     // Real
    CacheEntryIndex       cacheZ0Index;         // ZIndex
    CacheEntryIndex       cacheMomentumIndex;   // M*u
    CacheEntryIndex       cacheKEIndex;         // ~u*M*u/2
    ZIndex                workZIndex;           // power integral

friend class Force::Thermostat;
};

//-------------------------------- Thermostat ----------------------------------
//------------------------------------------------------------------------------

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Force::Thermostat, Force::ThermostatImpl, Force);

Force::Thermostat::Thermostat
   (GeneralForceSubsystem& forces, const SimbodyMatterSubsystem& matter,
    Real boltzmannsConstant, Real bathTemperature, Real relaxationTime,
    int numExcludedDofs) 
:   Force(new ThermostatImpl(matter, boltzmannsConstant, bathTemperature, 
                             relaxationTime, numExcludedDofs))
{
    SimTK_APIARGCHECK1_ALWAYS(boltzmannsConstant > 0, 
        "Force::Thermostat","ctor", "Illegal Boltzmann's constant %g.", boltzmannsConstant);
    SimTK_APIARGCHECK1_ALWAYS(bathTemperature > 0, 
        "Force::Thermostat","ctor", "Illegal bath temperature %g.", bathTemperature);
    SimTK_APIARGCHECK1_ALWAYS(relaxationTime > 0, 
        "Force::Thermostat","ctor", "Illegal relaxation time %g.", relaxationTime);
    SimTK_APIARGCHECK1_ALWAYS(0 <= numExcludedDofs && numExcludedDofs <= 6, 
        "Force::Thermostat","ctor", 
        "Illegal number of excluded rigid body dofs %d (must be 0-6).", 
        numExcludedDofs);

    updImpl().setForceSubsystem(forces, forces.adoptForce(*this));
}

Force::Thermostat& Force::Thermostat::
setDefaultNumChains(int numChains) {
    SimTK_APIARGCHECK1_ALWAYS(numChains > 0, 
        "Force::Thermostat","setDefaultNumChains", 
        "Illegal number of chains %d.", numChains);

    getImpl().invalidateTopologyCache();
    updImpl().defaultNumChains = numChains;
    return *this;
}

Force::Thermostat& Force::Thermostat::
setDefaultBathTemperature(Real bathTemperature) {
    SimTK_APIARGCHECK1_ALWAYS(bathTemperature > 0, 
        "Force::Thermostat","setDefaultBathTemperature", 
        "Illegal bath temperature %g.", bathTemperature);

    getImpl().invalidateTopologyCache();
    updImpl().defaultBathTemp = bathTemperature;
    return *this;
}

Force::Thermostat& Force::Thermostat::
setDefaultRelaxationTime(Real relaxationTime) {
    SimTK_APIARGCHECK1_ALWAYS(relaxationTime > 0, 
        "Force::Thermostat","setDefaultRelaxationTime", 
        "Illegal bath temperature %g.", relaxationTime);

    getImpl().invalidateTopologyCache();
    updImpl().defaultRelaxationTime = relaxationTime;
    return *this;
}

Force::Thermostat& Force::Thermostat::
setDefaultNumExcludedDofs(int numExcludedDofs) {
    SimTK_APIARGCHECK1_ALWAYS(0 <= numExcludedDofs && numExcludedDofs <= 6, 
        "Force::Thermostat","setDefaultNumExcludedDofs", 
        "Illegal number of excluded rigid body dofs %d (must be 0-6).", 
        numExcludedDofs);

    getImpl().invalidateTopologyCache();
    updImpl().defaultNumExcludedDofs = numExcludedDofs;
    return *this;
}

int Force::Thermostat::getDefaultNumChains() const {return getImpl().defaultNumChains;}
Real Force::Thermostat::getDefaultBathTemperature() const {return getImpl().defaultBathTemp;}
Real Force::Thermostat::getDefaultRelaxationTime() const {return getImpl().defaultRelaxationTime;}
int Force::Thermostat::getDefaultNumExcludedDofs() const {return getImpl().defaultNumExcludedDofs;}
Real Force::Thermostat::getBoltzmannsConstant() const {return getImpl().kB;}

const Force::Thermostat& Force::Thermostat::
setNumChains(State& s, int numChains) const {
    SimTK_APIARGCHECK1_ALWAYS(numChains > 0, 
        "Force::Thermostat","setNumChains", 
        "Illegal number of chains %d.", numChains);

    getImpl().updNumChains(s) = numChains;
    return *this;
}

const Force::Thermostat& Force::Thermostat::
setBathTemperature(State& s, Real bathTemperature) const {
    SimTK_APIARGCHECK1_ALWAYS(bathTemperature > 0, 
        "Force::Thermostat","setBathTemperature", 
        "Illegal bath temperature %g.", bathTemperature);

    getImpl().updBathTemp(s) = bathTemperature;
    return *this;
}

const Force::Thermostat& Force::Thermostat::
setRelaxationTime(State& s, Real relaxationTime) const {
    SimTK_APIARGCHECK1_ALWAYS(relaxationTime > 0, 
        "Force::Thermostat","setRelaxationTime", 
        "Illegal bath temperature %g.", relaxationTime);

    getImpl().updRelaxationTime(s) = relaxationTime;
    return *this;
}

const Force::Thermostat& Force::Thermostat::
setNumExcludedDofs(State& s, int numExcludedDofs) const {
    SimTK_APIARGCHECK1_ALWAYS(0 <= numExcludedDofs && numExcludedDofs <= 6, 
        "Force::Thermostat","setNumExcludedDofs", 
        "Illegal number of excluded rigid body dofs %d (must be 0-6).", 
        numExcludedDofs);

    getImpl().updNumExcludedDofs(s) = numExcludedDofs;
    return *this;
}

int Force::Thermostat::getNumChains(const State& s) const {return getImpl().getNumChains(s);}
Real Force::Thermostat::getBathTemperature(const State& s) const {return getImpl().getBathTemp(s);}
Real Force::Thermostat::getRelaxationTime(const State& s) const {return getImpl().getRelaxationTime(s);}
int Force::Thermostat::getNumExcludedDofs(const State& s) const {return getImpl().getNumExcludedDofs(s);}

void Force::Thermostat::initializeChainState(State& s) const {
    const ThermostatImpl& impl = getImpl();
    const int nChains = impl.getNumChains(s);
    for (int i=0; i < 2*nChains; ++i)
        impl.updZ(s, i) = 0;
}

void Force::Thermostat::setChainState(State& s, const Vector& z) const {
    const ThermostatImpl& impl = getImpl();
    const int nChains = impl.getNumChains(s);
    SimTK_APIARGCHECK2_ALWAYS(z.size() == 2*nChains,
        "Force::Thermostat", "setChainState", 
        "Number of values supplied (%d) didn't match twice the number of chains %d.", 
        z.size(), nChains);
    for (int i=0; i < 2*nChains; ++i)
        impl.updZ(s, i) = z[i];
}

Vector Force::Thermostat::getChainState(const State& s) const {
    const ThermostatImpl& impl = getImpl();
    const int nChains = impl.getNumChains(s);
    Vector out(2*nChains);
    for (int i=0; i < 2*nChains; ++i)
        out[i] = impl.getZ(s, i);
    return out;
}

// Cost is a divide and two flops, say about 20 flops.
Real Force::Thermostat::getCurrentTemperature(const State& s) const {
    const Real ke = getImpl().getKE(s); // Cached value for kinetic energy
    const Real kB = getImpl().kB;       // Boltzmann's constant
    const int  N  = getImpl().getNumThermalDOFs(s);
    return (2*ke) / (N*kB);
}

int Force::Thermostat::getNumThermalDofs(const State& s) const {
    return getImpl().getNumThermalDOFs(s);
}

// Bath energy is KEb + PEb where
//    KEb = 1/2 kT t^2 (N z0^2 + sum(zi^2))
//    PEb = kT (N s0 + sum(si))
// Cost is about 7 flops + 4 flops/chain, say about 25 flops.
Real Force::Thermostat::calcBathEnergy(const State& state) const {
    const ThermostatImpl& impl = getImpl();
    const int nChains = impl.getNumChains(state);
    const int  N = impl.getNumThermalDOFs(state);
    const Real kT = impl.kB * impl.getBathTemp(state);
    const Real t = impl.getRelaxationTime(state);

    Real zsqsum = N * square(impl.getZ(state,0));
    for (int i=1; i < nChains; ++i)
        zsqsum += square(impl.getZ(state,i));

    Real ssum = N * impl.getZ(state, nChains);
    for (int i=1; i < nChains; ++i)
        ssum += impl.getZ(state, nChains+i);

    const Real KEb = (kT/2) * t*t * zsqsum;
    const Real PEb = kT * ssum;

    return KEb + PEb;
}

Real Force::Thermostat::getExternalPower(const State& state) const {
    return getImpl().calcExternalPower(state);
}

Real Force::Thermostat::getExternalWork(const State& state) const {
    return getImpl().getExternalWork(state);
}

void Force::Thermostat::setExternalWork(State& state, Real w) const {
    getImpl().updExternalWork(state) = w;
}

// This is the number of thermal dofs. TODO: we're ignoring constraint 
// redundancy but we shouldn't be! That could result in negative dofs, so we'll 
// make sure that doesn't happen. But don't expect meaningful results
// in that case. Note that it is the acceleration-level constraints that
// matter; they remove dofs regardless of whether there is a corresponding
// velocity constraint.
int Force::ThermostatImpl::getNumThermalDOFs(const State& state) const {
    const int ndofs = state.getNU() - state.getNUDotErr()
                        - getNumExcludedDofs(state);
    const int N = std::max(1, ndofs);
    return N;
}

// This force produces only mobility forces, with 
//      f = -z0 * M * u
// Conveniently we already calculated the momentum M*u and cached it
// in realizeVelocity(). Additional cost here is 2*N flops.
void Force::ThermostatImpl::
calcForce(const State& state, Vector_<SpatialVec>&, Vector_<Vec3>&, 
          Vector& mobilityForces) const 
{
    const Vector& p = getMomentum(state);

    // Generate momentum-weighted forces and apply to mobilities.
    // This is 2*N flops.
    mobilityForces -= getZ(state, 0)*p;
}

// All the power generated by this force is external (to or from the
// thermal bath). So the power is
//      p = ~u * f = -z0 * (~u * M * u) = -z0 * (2*KE).
// We already calculated and cached KE in realizeVelocity() so power
// is practically free here (2 flops).
Real Force::ThermostatImpl::
calcExternalPower(const State& state) const {
    return -2 * getZ(state, 0) * getKE(state);
}

// Allocate and initialize state variables.
void Force::ThermostatImpl::realizeTopology(State& state) const {
    // Make these writable just here where we need to fill in the Topology "cache"
    // variables; after this they are const.
    Force::ThermostatImpl* mutableThis = const_cast<Force::ThermostatImpl*>(this);
    mutableThis->dvNumChains = 
        getForceSubsystem().allocateDiscreteVariable(state, Stage::Model, 
                                                     new Value<int>(defaultNumChains));
    mutableThis->dvBathTemp = 
        getForceSubsystem().allocateDiscreteVariable(state, Stage::Instance, 
                                                     new Value<Real>(defaultBathTemp));
    mutableThis->dvRelaxationTime = 
        getForceSubsystem().allocateDiscreteVariable(state, Stage::Instance, 
                                                     new Value<Real>(defaultRelaxationTime));
    mutableThis->dvNumExcludedDofs = 
        getForceSubsystem().allocateDiscreteVariable(state, Stage::Model, 
                                                     new Value<int>(defaultNumExcludedDofs));

    // This cache entry holds the auxiliary state index of our first
    // thermostat state variable. It is valid after realizeModel().
    mutableThis->cacheZ0Index = 
        getForceSubsystem().allocateCacheEntry(state, Stage::Model, 
                                               new Value<ZIndex>());

    // This cache entry holds the generalized momentum M*u. The vector
    // will be allocated to hold nu values.
    mutableThis->cacheMomentumIndex =
        getForceSubsystem().allocateCacheEntry(state, Stage::Velocity, 
                                               new Value<Vector>());

    // This cache entry holds the kinetic energy ~u*M*u/2.
    mutableThis->cacheKEIndex =
        getForceSubsystem().allocateCacheEntry(state, Stage::Velocity, 
                                               new Value<Real>(NaN));

    const Vector workZInit(1, Zero);
    mutableThis->workZIndex = 
        getForceSubsystem().allocateZ(state, workZInit);
}

// Allocate the chain state variables and bath energy variables.
// TODO: this should be done at Instance stage.
void Force::ThermostatImpl::realizeModel(State& state) const {
    const int nChains = getNumChains(state);
    const Vector zInit(2*nChains, Zero);
    updZ0Index(state) = getForceSubsystem().allocateZ(state, zInit);
}

// Calculate velocity-dependent terms, the internal coordinate
// momentum and the kinetic energy. This is the expensive part
// at about 125*N flops if all joints are 1 dof.
void Force::ThermostatImpl::realizeVelocity(const State& state) const {
    Vector& Mu = updMomentum(state);
    matter.multiplyByM(state, state.getU(), Mu); // <-- expensive ~123*N flops
    updKE(state) = (~state.getU() * Mu) / 2;     // 2*N flops
}

// Calculate time derivatives of the various state variables.
// This is just a fixed cost for the whole system, independent of
// size: 3 divides + a few flops per chain, maybe 100 flops total.
void Force::ThermostatImpl::realizeDynamics(const State& state) const {
    const Real t    = getRelaxationTime(state);
    const Real oot2 = 1 / square(t);
    const int  m    = getNumChains(state);

    // This is the desired kinetic energy per dof.
    const Real Eb = kB * getBathTemp(state) / 2;

    const int  N = getNumThermalDOFs(state);

    // This is the current average kinetic energy per dof.
    const Real E = getKE(state) / N;

    updZDot(state, 0) = (E/Eb - 1) * oot2;

    int Ndofs = N;  // only for z0
    for (int k=1; k < m; ++k) {
        const Real zk1 = getZ(state, k-1);
        const Real zk  = getZ(state, k);
        updZDot(state, k-1) -= zk1 * zk;
        updZDot(state, k) = Ndofs * square(zk1) - oot2;
        Ndofs = 1; // z1..m-1 control only 1 dof each
    }

    // Calculate sdot's for energy calculation.
    for (int k=0; k < m; ++k)
        updZDot(state, m+k) = getZ(state, k);

    updWorkZDot(state) = calcExternalPower(state);
}


} // namespace SimTK

