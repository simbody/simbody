/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2012-13 Stanford University and the Authors.        *
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

#include "SimTKcommon.h"

#include "simbody/internal/common.h"
#include "simbody/internal/MobilizedBody.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/CableSpring.h"

#include "ForceImpl.h"

using namespace SimTK;

//==============================================================================
//                          CABLE SPRING :: IMPL
//==============================================================================
// This is the hidden implementation class for a CableSpring force element.
class CableSpring::Impl : public ForceImpl {
friend class CableSpring;

    // Type of the discrete state variable that holds values for this
    // cable spring's changeable parameters in a State.
    struct InstanceVars {
        InstanceVars(Real defStiffness, Real defSlackLength, 
                     Real defDissipationCoef)
        :   k(defStiffness), L0(defSlackLength), c(defDissipationCoef) {}

        Real      k, L0, c;
    };

    // Type of the velocity-stage lazy cache entry that holds the spring
    // tension and instantaneous power if either has been requested.
    struct ForceCache {
        Real       f;               // scalar tension
        Real       powerLoss;
    };

    Impl(const CablePath&   path, 
         Real               stiffness,
         Real               slackLength,
         Real               dissipationCoef)
    :   path(path), defK(stiffness), defL0(slackLength), defC(dissipationCoef) 
    {
    }

    // Implementation of virtual methods from ForceImpl:
    Impl* clone() const override {return new Impl(*this);}
    bool dependsOnlyOnPositions() const override {return false;}

    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces,
                   Vector_<Vec3>& particleForces, Vector& mobilityForces) const
                   override 
    {   
        const ForceCache& forceCache = ensureForceCacheValid(state);
        path.applyBodyForces(state, forceCache.f, bodyForces); 
    }

    // We're not bothering to cache P.E. -- just recalculate it when asked.
    Real calcPotentialEnergy(const State& state) const override 
    {
        const InstanceVars& inst = getInstanceVars(state);
        const Real x = calcStretch(state, inst); // x >= 0
        const Real pe = inst.k * x*x / 2;
        return pe;
    }

    // Allocate the state variables and cache entry. 
    void realizeTopology(State& s) const override {
        // Allocate the discrete variable for instance parameters.
        const InstanceVars iv(defK,defL0,defC);
        instanceVarsIx = getForceSubsystem()
            .allocateDiscreteVariable(s, Stage::Instance, 
                                      new Value<InstanceVars>(iv));

        // Allocate a continuous variable to hold the integrated power loss.
        const Vector einit(1, Real(0));
        dissipatedEnergyIx = getForceSubsystem().allocateZ(s,einit);

        // Allocate a lazy cache entry for force and power loss.
        forceCacheIx = getForceSubsystem().allocateLazyCacheEntry(s,
            Stage::Velocity, new Value<ForceCache>());
    }

    // We must always evaluate the power if we're going to calculate its
    // integral.
    void realizeAcceleration(const State& s) const {
        const ForceCache& forceCache = ensureForceCacheValid(s);
        updDissipatedEnergyDeriv(s) = forceCache.powerLoss;
    }
        
    // Return the amount x by which the cable is stretched beyond its slack
    // length or zero if the cable is slack. Must be at stage Position.
    Real calcStretch(const State& state, const InstanceVars& inst) const {
        const Real x0 = path.getCableLength(state) - inst.L0; // signed
        return std::max(Real(0), x0); // return x >= 0
    }

    // Calculate tension f=f_stretch+f_rate, and power loss f_rate*xdot. See 
    // theory discussion for CableSpring class for an explanation. Must be at 
    // stage Velocity. 
    void calcTensionAndPowerLoss(const State& state, 
                                 Real& f, Real& powerLoss) const 
    {
        const InstanceVars& inst = getInstanceVars(state);
        const Real x = calcStretch(state, inst); // >= 0
        if (x == 0) {powerLoss = f = 0; return;} // cable spring is slack

        const Real f_stretch = inst.k * x;
        const Real xdot = path.getCableLengthDot(state);
        const Real diss = f_stretch*inst.c*xdot;
        const Real f_rate = std::max(-f_stretch, diss);
        f = f_stretch + f_rate;    // f=0 if dissipation too negative
        powerLoss = f_rate * xdot; // but can still dissipate power
    }

    // If state is at stage Velocity, we can calculate and store tension
    // in the cache if it hasn't already been calculated.
    const ForceCache& ensureForceCacheValid(const State& state) const {
        if (isForceCacheValid(state)) 
            return getForceCache(state);
        ForceCache& forceCache = updForceCache(state);
        calcTensionAndPowerLoss(state, forceCache.f, forceCache.powerLoss);
        markForceCacheValid(state);
        return forceCache;
    }

    const InstanceVars& getInstanceVars(const State& s) const
    {   return Value<InstanceVars>::downcast
           (getForceSubsystem().getDiscreteVariable(s,instanceVarsIx)); }
    InstanceVars& updInstanceVars(State& s) const
    {   return Value<InstanceVars>::updDowncast
           (getForceSubsystem().updDiscreteVariable(s,instanceVarsIx)); }

    const Real& getDissipatedEnergyVar(const State& s) const
    {   return getForceSubsystem().getZ(s)[dissipatedEnergyIx]; }
    Real& updDissipatedEnergyVar(State& s) const
    {   return getForceSubsystem().updZ(s)[dissipatedEnergyIx]; }
    Real& updDissipatedEnergyDeriv(const State& s) const
    {   return getForceSubsystem().updZDot(s)[dissipatedEnergyIx]; }

    const ForceCache& getForceCache(const State& s) const
    {   return Value<ForceCache>::downcast
            (getForceSubsystem().getCacheEntry(s,forceCacheIx)); }
    ForceCache& updForceCache(const State& s) const
    {   return Value<ForceCache>::updDowncast
            (getForceSubsystem().updCacheEntry(s,forceCacheIx)); }

    bool isForceCacheValid(const State& s) const
    {   return getForceSubsystem().isCacheValueRealized(s,forceCacheIx); }
    void markForceCacheValid(const State& s) const
    {   getForceSubsystem().markCacheValueRealized(s,forceCacheIx); }


//------------------------------------------------------------------------------
    // TOPOLOGY STATE
    const CablePath                 path; // shallow, reference counted copy
    Real                            defK, defL0, defC;

    // TOPOLOGY CACHE (Set only once in realizeTopology(); const thereafter.)
    mutable DiscreteVariableIndex   instanceVarsIx; // k, L0, c
    mutable ZIndex                  dissipatedEnergyIx;
    mutable CacheEntryIndex         forceCacheIx;
};


//==============================================================================
//                          CABLE SPRING DEFINITIONS
//==============================================================================

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(CableSpring, 
                                        CableSpring::Impl, Force);

CableSpring::CableSpring
   (GeneralForceSubsystem&  forces, 
    const CablePath&        path,
    Real                    defaultStiffness,
    Real                    defaultSlackLength,
    Real                    defaultDissipationCoef)
:   Force(new CableSpring::Impl(path, defaultStiffness, defaultSlackLength, 
                                defaultDissipationCoef))
{
    SimTK_ERRCHK_ALWAYS(defaultStiffness >= 0,
        "CableSpring constructor", "Spring constant must be nonnegative.");
    SimTK_ERRCHK_ALWAYS(defaultSlackLength >= 0,
        "CableSpring constructor", "Slack length must be nonnegative.");
    SimTK_ERRCHK_ALWAYS(defaultDissipationCoef >= 0,
        "CableSpring constructor", "Dissipation coefficient must be nonnegative.");
    updImpl().setForceSubsystem(forces, forces.adoptForce(*this));
}

CableSpring& CableSpring::
setDefaultStiffness(Real stiffness) {
    SimTK_ERRCHK_ALWAYS(stiffness >= 0,
        "CableSpring::setDefaultStiffness()", 
        "Spring constant must be nonnegative.");
    getImpl().invalidateTopologyCache();
    updImpl().defK = stiffness; 
    return *this;
}
CableSpring& CableSpring::
setDefaultSlackLength(Real slackLength) {
    SimTK_ERRCHK_ALWAYS(slackLength >= 0,
        "CableSpring::setDefaultSlackLength()", 
        "Slack length must be nonnegative.");
    getImpl().invalidateTopologyCache();
    updImpl().defL0 = slackLength; 
    return *this;
}
CableSpring& CableSpring::
setDefaultDissipationCoef(Real dissipationCoef) {
    SimTK_ERRCHK_ALWAYS(dissipationCoef >= 0,
        "CableSpring::setDefaultDissipationCoef()",
        "Dissipation coefficient must be nonnegative.");
    getImpl().invalidateTopologyCache();
    updImpl().defC = dissipationCoef; 
    return *this;
}

Real CableSpring::
getDefaultStiffness() const {return getImpl().defK;}
Real CableSpring::
getDefaultSlackLength() const {return getImpl().defL0;}
Real CableSpring::
getDefaultDissipationCoef() const {return getImpl().defC;}

const CableSpring& CableSpring::
setStiffness(State& state, Real stiffness) const {
    SimTK_ERRCHK_ALWAYS(stiffness >= 0,
        "CableSpring::setStiffness()",
        "Spring constant must be nonnegative.");
    getImpl().updInstanceVars(state).k = stiffness; 
    return *this;
}
const CableSpring& CableSpring::
setSlackLength(State& state, Real slackLength) const {
    SimTK_ERRCHK_ALWAYS(slackLength >= 0,
        "CableSpring::setSlackLength()",
        "Slack length must be nonnegative.");
    getImpl().updInstanceVars(state).L0 = slackLength; 
    return *this;
}
const CableSpring& CableSpring::
setDissipationCoef(State& state, Real dissipationCoef) const {
    SimTK_ERRCHK_ALWAYS(dissipationCoef >= 0,
        "CableSpring::setDissipationCoef()",
        "Dissipation coefficient must be nonnegative.");
    getImpl().updInstanceVars(state).c = dissipationCoef; 
    return *this;
}

Real CableSpring::
getStiffness(const State& state) const 
{   return getImpl().getInstanceVars(state).k; }
Real CableSpring::
getSlackLength(const State& state) const 
{   return getImpl().getInstanceVars(state).L0; }
Real CableSpring::
getDissipationCoef(const State& state) const 
{   return getImpl().getInstanceVars(state).c; }

Real CableSpring::
getLength(const State& state) const
{   return getImpl().path.getCableLength(state); }
Real CableSpring::
getLengthDot(const State& state) const
{   return getImpl().path.getCableLengthDot(state); }


Real CableSpring::
getTension(const State& s) const
{   getImpl().ensureForceCacheValid(s); 
    return getImpl().getForceCache(s).f; }

Real CableSpring::
getPotentialEnergy(const State& s) const
{   return getImpl().calcPotentialEnergy(s); }

Real CableSpring::
getPowerDissipation(const State& s) const
{   getImpl().ensureForceCacheValid(s); 
    return getImpl().getForceCache(s).powerLoss; }

Real CableSpring::
getDissipatedEnergy(const State& s) const
{   return getImpl().getDissipatedEnergyVar(s); }

void CableSpring::
setDissipatedEnergy(State& s, Real energy) const {
    SimTK_ERRCHK1_ALWAYS(energy >= 0,
        "CableSpring::setDissipatedEnergy()",
        "The initial value for the dissipated energy must be nonnegative"
        " but an attempt was made to set it to %g.", energy);
    getImpl().updDissipatedEnergyVar(s) = energy; 
}

const CablePath& CableSpring::
getCablePath() const
{   return getImpl().path; }


