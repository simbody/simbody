/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-13 Stanford University and the Authors.        *
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
#include "simbody/internal/Force_Gravity.h"

#include "ForceImpl.h"

namespace SimTK {

//==============================================================================
//                          FORCE :: GRAVITY IMPL
//==============================================================================
// This is the hidden implementation class for Force::Gravity.
class Force::GravityImpl : public ForceImpl {
friend class Force::Gravity;

    // These are settable parameters including gravity vector, zero height,
    // and which if any mobilized bodies are immune to gravity. Modifying this
    // variable invalidates Dynamics stage automatically, but every modification
    // must explicitly invalidate the force cache because that can be filled in
    // as early as Position stage.
    struct Parameters {
        Parameters(const UnitVec3& defDirection, 
                     Real defMagnitude, Real defZeroHeight, 
                     const Array_<bool,MobilizedBodyIndex>& defMobodIsImmune)
        :   d(defDirection), g(defMagnitude), z(defZeroHeight),
            mobodIsImmune(defMobodIsImmune) {}
        UnitVec3    d;
        Real        g, z;
        Array_<bool,MobilizedBodyIndex> mobodIsImmune; // [nb]
    };

    // The cache has a SpatialVec for each mobilized body, a Vec3 for each
    // particle [not used], and a scalar for potential energy. The SpatialVec 
    // corresponding to Ground is initialized to zero and stays that way.
    // This is a lazy-evaluated cache entry that can be calculated any time 
    // after Position stage. It is automatically invalidated if a Position
    // stage variable changes. But, it is also dependent on the parameter
    // values in the discrete Parameters variable and must be invalidated
    // explicitly if any of those change.
    struct ForceCache {
        ForceCache() {}
        void allocate(int nb, int np, bool initToZero)
        {   F_GB.resize(nb); f_GP.resize(np); 
            if (initToZero) setToZero(); else setToNaN(); }
        void setToZero() {F_GB.setToZero(); f_GP.setToZero(); pe=0;}
        void setToNaN()  
        {   F_GB.setToNaN(); F_GB[0]=SpatialVec(Vec3(0)); // Ground
            f_GP.setToNaN(); pe=NaN;}
        Vector_<SpatialVec> F_GB; // rigid body forces
        Vector_<Vec3>       f_GP; // particle forces
        Real                pe;   // total potential energy
    };

    // Constructor from a direction and magnitude.
    GravityImpl(const SimbodyMatterSubsystem&   matter,
                const UnitVec3&                 direction,
                Real                            magnitude,
                Real                            zeroHeight)
    :   matter(matter), defDirection(direction), defMagnitude(magnitude), 
        defZeroHeight(zeroHeight), 
        defMobodIsImmune(matter.getNumBodies(), false),
        numEvaluations(0)
    {   defMobodIsImmune.front() = true; } // Ground is always immune

    // Constructor from a gravity vector, which might have zero magnitude.
    // In that case we'll negate the default up direction in the System as the
    // default direction. That would only be noticeable if the magnitude were
    // later increased without giving a new direction.
    GravityImpl(const SimbodyMatterSubsystem&   matter,
                const Vec3&                     gravityVec,
                Real                            zeroHeight)
    :   matter(matter)
    {   // Invoke the general constructor using system up direction. 
        new (this) GravityImpl(matter, -matter.getSystem().getUpDirection(),
                               gravityVec.norm(), zeroHeight);
        if (defMagnitude > 0) 
            defDirection=UnitVec3(gravityVec/defMagnitude,true);
    }

    // Constructor from just gravity magnitude; use negative of System "up"
    // direction as the down direction here.
    GravityImpl(const SimbodyMatterSubsystem&   matter,
                Real                            magnitude,
                Real                            zeroHeight)
    :   matter(matter)
    {   // Invoke the general constructor using system up direction. 
        new (this) GravityImpl(matter, -matter.getSystem().getUpDirection(),
                               magnitude, zeroHeight);
    }

    void setMobodIsImmuneByDefault(MobilizedBodyIndex mbx, bool isImmune) {
        if (mbx == 0) return; // can't change Ground's innate immunity
        invalidateTopologyCache();
        if (defMobodIsImmune.size() < mbx+1)
            defMobodIsImmune.resize(mbx+1, false);
        defMobodIsImmune[mbx] = isImmune;
    }

    bool getMobodIsImmuneByDefault(MobilizedBodyIndex mbx) const {
        if (defMobodIsImmune.size() < mbx+1)
            return false;
        return defMobodIsImmune[mbx];
    }

    void setMobodIsImmune(State& state, MobilizedBodyIndex mbx,
                          bool isImmune) const {
        if (mbx == 0) return; // no messing with Ground
        invalidateForceCache(state);
        Parameters& p = updParameters(state);
        p.mobodIsImmune[mbx] = isImmune;
    }

    bool getMobodIsImmune(const State& state, MobilizedBodyIndex mbx) const {
        const Parameters& p = getParameters(state);
        return p.mobodIsImmune[mbx];
    }

    GravityImpl* clone() const OVERRIDE_11 {
        return new GravityImpl(*this);
    }

    // We are doing our own caching here, so don't override the 
    // dependsOnlyOnPositions() method which would cause the base class also
    // to cache the results.

    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces,
                   Vector_<Vec3>& particleForces, Vector& mobilityForces) const
                   OVERRIDE_11;
    Real calcPotentialEnergy(const State& state) const OVERRIDE_11;

    // Allocate the state variables and cache entries.
    void realizeTopology(State& s) const OVERRIDE_11;

    const Parameters& getParameters(const State& s) const
    {   return Value<Parameters>::downcast
           (getForceSubsystem().getDiscreteVariable(s,parametersIx)); }
    Parameters& updParameters(State& s) const
    {   return Value<Parameters>::updDowncast
           (getForceSubsystem().updDiscreteVariable(s,parametersIx)); }

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
    void invalidateForceCache(const State& s) const
    {   getForceSubsystem().markCacheValueNotRealized(s,forceCacheIx); }

    // This method calculates gravity forces if needed, and bumps the 
    // numEvaluations counter if it has to do any work.
    void ensureForceCacheValid(const State&) const;

    // TOPOLOGY STATE
    const SimbodyMatterSubsystem&   matter;
    UnitVec3                        defDirection;
    Real                            defMagnitude;
    Real                            defZeroHeight;
    Array_<bool,MobilizedBodyIndex> defMobodIsImmune;

    // TOPOLOGY CACHE
    DiscreteVariableIndex           parametersIx;
    CacheEntryIndex                 forceCacheIx;

    mutable long long               numEvaluations;
};


//==============================================================================
//                             FORCE :: GRAVITY
//==============================================================================

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Force::Gravity, 
                                        Force::GravityImpl, Force);

Force::Gravity::Gravity
   (GeneralForceSubsystem&          forces, 
    const SimbodyMatterSubsystem&   matter,
    const UnitVec3&                 defDirection,
    Real                            defMagnitude,
    Real                            defZeroHeight)
:   Force(new GravityImpl(matter,defDirection,defMagnitude,defZeroHeight))
{
    SimTK_ERRCHK1_ALWAYS(defMagnitude >= 0,
        "Force::Gravity::Gravity(downDirection,magnitude)",
        "The gravity magnitude g must be nonnegative but was specified as %g.",
        defMagnitude);
    SimTK_ERRCHK_ALWAYS(defDirection.isFinite(),
        "Force::Gravity::Gravity(downDirection,magnitude)",
        "A non-finite 'down' direction was received; did you specify a zero-"
        "length Vec3? The direction must be non-zero.");
    updImpl().setForceSubsystem(forces, forces.adoptForce(*this));
}

Force::Gravity::Gravity
   (GeneralForceSubsystem&          forces, 
    const SimbodyMatterSubsystem&   matter,
    const Vec3&                     defGravity)
:   Force(new GravityImpl(matter,defGravity,0))
{
    updImpl().setForceSubsystem(forces, forces.adoptForce(*this));
}

Force::Gravity::Gravity
   (GeneralForceSubsystem&          forces, 
    const SimbodyMatterSubsystem&   matter,
    Real                            defMagnitude)
:   Force(new GravityImpl(matter,defMagnitude,0))
{
    SimTK_ERRCHK1_ALWAYS(defMagnitude >= 0,
        "Force::Gravity::Gravity(magnitude)",
        "The gravity magnitude g must be nonnegative but was specified as %g.",
        defMagnitude);
    updImpl().setForceSubsystem(forces, forces.adoptForce(*this));
}

// Each of the setDefault methods must invalidate the topology cache.

Force::Gravity& Force::Gravity::
setDefaultBodyIsExcluded(MobilizedBodyIndex mobod, bool isExcluded) {
    // Invalidates topology cache.
    updImpl().setMobodIsImmuneByDefault(mobod, isExcluded);
    return *this;
}

Force::Gravity& Force::Gravity::
setDefaultGravityVector(const Vec3& gravity) {
    const Real g = gravity.norm();
    getImpl().invalidateTopologyCache();
    updImpl().defMagnitude = g;
    // Don't change the direction if the magnitude is zero.
    if (g > 0)
        updImpl().defDirection = UnitVec3(gravity/g, true);
    return *this;
}

Force::Gravity& Force::Gravity::
setDefaultDownDirection(const UnitVec3& down) {
    SimTK_ERRCHK_ALWAYS(down.isFinite(),
        "Force::Gravity::setDefaultDownDirection()",
        "A non-finite 'down' direction was received; did you specify a zero-"
        "length Vec3? The direction must be non-zero.");

    getImpl().invalidateTopologyCache();
    updImpl().defDirection = down;
    return *this;
}

Force::Gravity& Force::Gravity::
setDefaultMagnitude(Real g) {
    SimTK_ERRCHK1_ALWAYS(g >= 0,
        "Force::Gravity::setDefaultMagnitude()",
        "The gravity magnitude g must be nonnegative but was specified as %g.",
        g);

    getImpl().invalidateTopologyCache();
    updImpl().defMagnitude = g;
    return *this;
}

Force::Gravity& Force::Gravity::
setDefaultZeroHeight(Real zeroHeight) {
    getImpl().invalidateTopologyCache();
    updImpl().defZeroHeight = zeroHeight;
    return *this;
}

bool Force::Gravity::
getDefaultBodyIsExcluded(MobilizedBodyIndex mobod) const
{   return getImpl().getMobodIsImmuneByDefault(mobod); }
Vec3 Force::Gravity::
getDefaultGravityVector() const 
{   return getImpl().defDirection * getImpl().defMagnitude;}
const UnitVec3& Force::Gravity::
getDefaultDownDirection() const {return getImpl().defDirection;}
Real Force::Gravity::
getDefaultMagnitude() const {return getImpl().defMagnitude;}
Real Force::Gravity::
getDefaultZeroHeight() const {return getImpl().defZeroHeight;}


// These set routines must explicitly invalidate the force cache since only
// Dynamics stage is invalidated automatically by this state variable. Try
// not to do anything if the new value is the same as the old one.

const Force::Gravity& Force::Gravity::
setBodyIsExcluded(State& state, MobilizedBodyIndex mobod, 
                  bool isExcluded) const
{   
    const GravityImpl& impl = getImpl();
    SimTK_ERRCHK2_ALWAYS(mobod < impl.matter.getNumBodies(),
        "Force::Gravity::setBodyIsExcluded()",
        "Attemped to exclude mobilized body with index %d but only mobilized"
        " bodies with indices between 0 and %d exist in this System.", 
        (int)mobod, impl.matter.getNumBodies()-1);

    if (getBodyIsExcluded(state, mobod) != isExcluded) {
        // Invalidates force cache.
        impl.setMobodIsImmune(state, mobod, isExcluded);
        // The zero must be precalculated if the body is immune to gravity.
        SpatialVec& F = impl.updForceCache(state).F_GB[mobod];
        if (isExcluded) F.setToZero(); else F.setToNaN();
    }
    return *this;
}

const Force::Gravity& Force::Gravity::
setGravityVector(State& state, const Vec3& gravity) const {
    const Real     newg = gravity.norm();
    const UnitVec3 newd = newg > 0 ? UnitVec3(gravity/newg, true)
                                   : getDownDirection(state);
    if (getMagnitude(state) != newg || getDownDirection(state) != newd) {
        getImpl().invalidateForceCache(state);
        getImpl().updParameters(state).g = newg; 
        getImpl().updParameters(state).d = newd; 

        if (newg == 0) 
            getImpl().updForceCache(state).setToZero(); // must precalculate
    }
    return *this;
}

const Force::Gravity& Force::Gravity::
setDownDirection(State& state, const UnitVec3& down) const {
    SimTK_ERRCHK_ALWAYS(down.isFinite(),
        "Force::Gravity::setDownDirection()",
        "A non-finite 'down' direction was received; did you specify a zero-"
        "length Vec3? The direction must be non-zero.");

    if (getDownDirection(state) != down) {
        getImpl().invalidateForceCache(state);
        getImpl().updParameters(state).d = down; 
    }
    return *this;
}

const Force::Gravity& Force::Gravity::
setMagnitude(State& state, Real g) const {
    SimTK_ERRCHK1_ALWAYS(g >= 0,
        "Force::Gravity::setMagnitude()",
        "The gravity magnitude g must be nonnegative but was specified as %g.",
        g);

    if (getMagnitude(state) != g) {
        getImpl().invalidateForceCache(state);
        getImpl().updParameters(state).g = g; 

        if (g == 0) 
            getImpl().updForceCache(state).setToZero(); // must precalculate
    }
    return *this;
}

const Force::Gravity& Force::Gravity::
setZeroHeight(State& state, Real zeroHeight) const {
    if (getZeroHeight(state) != zeroHeight) {
        getImpl().invalidateForceCache(state);
        getImpl().updParameters(state).z = zeroHeight; 
    }
    return *this;
}

bool Force::Gravity::
getBodyIsExcluded(const State& state, MobilizedBodyIndex mobod) const
{   return getImpl().getMobodIsImmune(state, mobod); }


Vec3 Force::Gravity::
getGravityVector(const State& state) const
{   const GravityImpl::Parameters& p = getImpl().getParameters(state);
    return p.g*p.d; }

const UnitVec3& Force::Gravity::
getDownDirection(const State& state) const
{   return getImpl().getParameters(state).d; }


Real Force::Gravity::getMagnitude(const State& state) const
{   return getImpl().getParameters(state).g; }

Real Force::Gravity::getZeroHeight(const State& state) const
{   return getImpl().getParameters(state).z; }

Real Force::Gravity::
getPotentialEnergy(const State& s) const
{   getImpl().ensureForceCacheValid(s); 
    return getImpl().getForceCache(s).pe; }

const Vector_<SpatialVec>& Force::Gravity::
getBodyForces(const State& s) const
{   getImpl().ensureForceCacheValid(s);
    const GravityImpl::ForceCache& fc = getImpl().getForceCache(s);
    return fc.F_GB; }

const Vector_<Vec3>& Force::Gravity::
getParticleForces(const State& s) const
{   getImpl().ensureForceCacheValid(s);
    const GravityImpl::ForceCache& fc = getImpl().getForceCache(s);
    return fc.f_GP; }

long long Force::Gravity::
getNumEvaluations() const
{   return getImpl().numEvaluations; }

bool Force::Gravity::
isForceCacheValid(const State& state) const
{   return getImpl().isForceCacheValid(state); }

void Force::Gravity::
invalidateForceCache(const State& state) const 
{   getImpl().invalidateForceCache(state); }


//==============================================================================
//                            FORCE :: GRAVITY IMPL
//==============================================================================

//----------------------------- REALIZE TOPOLOGY -------------------------------
// Allocate the state variables and cache entries. The force cache is a lazy-
// evaluation entry - although it can be calculated any time after
// Stage::Position, it won't be unless someone asks for it. And if it is 
// evaluated early, it should not be re-evaluated when used as a force during
// Stage::Dynamics (via the calcForce() call).
//
// In addition, the force cache has a dependency on the parameter values that
// are stored in a discrete state variable. Changes to that variable 
// automatically invalidate Dynamics stage, but must be carefully managed also
// to invalidate the force cache here since it is only Position-dependent.
void Force::GravityImpl::
realizeTopology(State& s) const {
    GravityImpl* mThis = const_cast<GravityImpl*>(this);
    const int nb=matter.getNumBodies(), np=matter.getNumParticles();

    // In case more mobilized bodies were added after this Gravity element
    // was constructed, make room for the rest now. Earlier default immunity
    // settings are preserved.
    if (defMobodIsImmune.size() != nb)
        mThis->defMobodIsImmune.resize(nb, false);

    // Allocate a discrete state variable to hold parameters; see above comment.
    const Parameters p(defDirection,defMagnitude,defZeroHeight,
                       defMobodIsImmune); // initial value
    mThis->parametersIx = getForceSubsystem()
        .allocateDiscreteVariable(s, Stage::Dynamics, new Value<Parameters>(p));

    // Don't allocate force cache space yet since we would have to copy it.
    // Caution -- dependence on Parameters requires manual invalidation.
    mThis->forceCacheIx = getForceSubsystem().allocateLazyCacheEntry(s,
        Stage::Position, new Value<ForceCache>());

    // Now allocate the appropriate amount of space, and set to zero now 
    // any forces that we know will end up zero so we don't have to calculate
    // them at run time. Precalculated zeroes must be provided for any
    // immune elements, or for all elements if g==0, and this must be kept
    // up to date if there are runtime changes to the parameters.

    ForceCache& fc = updForceCache(s);
    if (defMagnitude == 0) {
        fc.allocate(nb, np, true); // initially zero since no gravity
    } else {
        fc.allocate(nb, np, false); // initially NaN except for Ground
        for (MobilizedBodyIndex mbx(1); mbx < nb; ++mbx)
            if (defMobodIsImmune[mbx])
                fc.F_GB[mbx] = SpatialVec(Vec3(0),Vec3(0));
        // This doesn't mean the ForceCache is valid yet.
    }      
}


//------------------------- ENSURE FORCE CACHE VALID ---------------------------
// This will also calculate potential energy since we can do it on the cheap 
// simultaneously with the force. Note that if the strength of gravity was set 
// to zero then we already zeroed out the forces and pe during realizeInstance()
// so all we have to do in that case is mark the cache valid now. Also, any
// immune bodies have their force set to zero in realizeInstance() so we don't
// have to do it again here.
void Force::GravityImpl::
ensureForceCacheValid(const State& state) const {
    if (isForceCacheValid(state)) return;

    SimTK_STAGECHECK_GE_ALWAYS(state.getSystemStage(), Stage::Position, 
        "Force::GravityImpl::ensureForceCacheValid()");

    const Parameters& p = getParameters(state);
    if (p.g == 0) { // no gravity
        markForceCacheValid(state); // zeroes must have been precalculated
        return;
    }

    // Gravity is non-zero and gravity forces are not up to date, so this counts
    // as an evaluation.
    ++numEvaluations;

    const Vec3 gravity      = p.g * p.d;
    const Real zeroPEOffset = p.g * p.z;
    ForceCache& fc = updForceCache(state);
    fc.pe = 0;

    const int nb = matter.getNumBodies();
    // Skip Ground since we know it is immune.
    for (MobilizedBodyIndex mbx(1); mbx < nb; ++mbx) {
        if (p.mobodIsImmune[mbx])
            continue; // don't apply gravity to this body; F already zero

        const MobilizedBody&     mobod  = matter.getMobilizedBody(mbx);
        const MassProperties&    mprops = mobod.getBodyMassProperties(state);
        const Transform&         X_GB   = mobod.getBodyTransform(state);

        Real        m       = mprops.getMass();
        const Vec3& p_CB    = mprops.getMassCenter(); // in B
        const Vec3  p_CB_G  = X_GB.R()*p_CB;          // exp. in G; 15 flops
        const Vec3  p_G_CB  = X_GB.p() + p_CB_G;      // meas. in G; 3 flops

        const Vec3  F_CB_G  = m*gravity; // force at mass center; 3 flops
        fc.F_GB[mbx] = SpatialVec(p_CB_G % F_CB_G, F_CB_G); // body frc; 9 flops 

        // odd signs here because height is in -gravity direction.
        fc.pe -= m*(~gravity*p_G_CB + zeroPEOffset); // 8 flops
    }

    const int np = matter.getNumParticles();
    if (np) {
        const Vector&        m    = matter.getAllParticleMasses(state);
        const Vector_<Vec3>& p_GP = matter.getAllParticleLocations(state);
        for (ParticleIndex px(0); px < np; ++px) {
            fc.f_GP[px] = m[px] * gravity;                     // 3 flops
            fc.pe -= m[px]*(~gravity*p_GP[px] + zeroPEOffset); // 8 flops
        }
    }

    markForceCacheValid(state);
}

//------------------------------- CALC FORCE -----------------------------------
void Force::GravityImpl::
calcForce(const State& state, Vector_<SpatialVec>& bodyForces, 
          Vector_<Vec3>& particleForces, Vector& mobilityForces) const 
{   ensureForceCacheValid(state);
    const ForceCache& fc = getForceCache(state);
    bodyForces     += fc.F_GB;
    particleForces += fc.f_GP; }


//-------------------------- CALC POTENTIAL ENERGY -----------------------------
// If the force was calculated, then the potential energy will already
// be valid. Otherwise we'll have to calculate it.
Real Force::GravityImpl::
calcPotentialEnergy(const State& state) const 
{   ensureForceCacheValid(state);
    const ForceCache& fc = getForceCache(state);
    return fc.pe; }


} // namespace SimTK

