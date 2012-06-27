/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-12 Stanford University and the Authors.        *
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

class Force::GravityImpl : public ForceImpl {
    struct InstanceVars {
        InstanceVars(const UnitVec3& defDirection, Real defMagnitude,
                     Real defZeroHeight)
        :   d(defDirection), g(defMagnitude), z(defZeroHeight) {}
        UnitVec3  d;
        Real      g, z;
    };

    // The cache has a SpatialVec for each mobilized body and a Vec3 for each
    // particle. The SpatialVec corresponding to Ground is initialized to zero
    // and stays that way.
    struct ForceCache {
        ForceCache() {}
        void allocate(int nb, int np)
        {   F_GB.resize(nb); f_GP.resize(np); setToNaN(); }
        void setToZero() {F_GB.setToZero(); f_GP.setToZero(); pe=0;}
        void setToNaN()  
        {   F_GB.setToNaN(); F_GB[0]=SpatialVec(Vec3(0)); // Ground
            f_GP.setToNaN(); pe=NaN;}
        Vector_<SpatialVec> F_GB; // rigid body forces
        Vector_<Vec3>       f_GP; // particle forces
        Real                pe;   // total potential energy
    };

public:
    GravityImpl(const SimbodyMatterSubsystem&   matter,
                const UnitVec3&                 direction,
                Real                            magnitude,
                Real                            zeroHeight)
    :   matter(matter), defDirection(direction), defMagnitude(magnitude), 
        defZeroHeight(zeroHeight) {}

    GravityImpl* clone() const {
        return new GravityImpl(*this);
    }
    bool dependsOnlyOnPositions() const {
        return true;
    }
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces,
                   Vector_<Vec3>& particleForces, Vector& mobilityForces) const;
    Real calcPotentialEnergy(const State& state) const;

    // Allocate the state variables and cache entries. The cached values are
    // lazy-evaluation entries - be sure to check whether they have already
    // been calculated; calculate them if not; and then mark them done.
    // They will be invalidated when the indicated stage has changed and
    // can be recalculated any time after that stage is realized.
    void realizeTopology(State& s) const {
        GravityImpl* mThis = const_cast<GravityImpl*>(this);

        const InstanceVars iv(defDirection,defMagnitude,defZeroHeight);
        mThis->instanceVarsIx = getForceSubsystem()
            .allocateDiscreteVariable(s, Stage::Instance, 
                                      new Value<InstanceVars>(iv));

        // Don't allocate force cache space yet since we have to copy
        // into the Value element.
        mThis->forceCacheIx = getForceSubsystem().allocateCacheEntry(s,
            Stage::Position, Stage::Infinity, new Value<ForceCache>());

        // Now allocate the appropriate amount of space.
        ForceCache& fc = updForceCache(s);
        fc.allocate(matter.getNumBodies(), matter.getNumParticles());
    }

    // If the magnitude of gravity was set to zero then we can calculate all
    // the forces on the affected bodies now -- they are zero!
    void realizeInstance(const State& s) const {
        const InstanceVars& iv = getInstanceVars(s);
        ForceCache&         fc = updForceCache(s);
        if (iv.g == 0) fc.setToZero();
        else           fc.setToNaN();
        // This doesn't mean the ForceCache is valid yet.
    }

private:
    const InstanceVars& getInstanceVars(const State& s) const
    {   return Value<InstanceVars>::downcast
           (getForceSubsystem().getDiscreteVariable(s,instanceVarsIx)); }
    InstanceVars& updInstanceVars(State& s) const
    {   return Value<InstanceVars>::updDowncast
           (getForceSubsystem().updDiscreteVariable(s,instanceVarsIx)); }

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
    void ensureForceCacheValid(const State&) const;

    // TOPOLOGY STATE
    const SimbodyMatterSubsystem&   matter;
    std::set<MobilizedBodyIndex>    affectedBodySet; // no dups or Ground
    UnitVec3                        defDirection;
    Real                            defMagnitude;
    Real                            defZeroHeight;

    // TOPOLOGY CACHE
    DiscreteVariableIndex           instanceVarsIx;
    CacheEntryIndex                 forceCacheIx;

friend class Force::Gravity;
friend std::ostream& operator<<(std::ostream&,const InstanceVars&);
friend std::ostream& operator<<(std::ostream&,const ForceCache&);
};

// These are required by Value<T>.
inline std::ostream& operator<<
   (std::ostream& o, const Force::GravityImpl::InstanceVars& iv)
{   assert(!"implemented"); return o; }
inline std::ostream& operator<<
   (std::ostream& o, const Force::GravityImpl::ForceCache& fc)
{   assert(!"implemented"); return o; }


//--------------------------------- Gravity ------------------------------------
//------------------------------------------------------------------------------

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
        "Force::Gravity::ctor(downDirection,magnitude)",
        "The gravity magnitude g must be nonnegative but was specified as %g.",
        defMagnitude);
    SimTK_ERRCHK_ALWAYS(defDirection.isFinite(),
        "Force::Gravity::ctor(downDirection,magnitude)",
        "A non-finite 'down' direction was received; did you specify a zero-"
        "length Vec3? The direction must be non-zero.");
    updImpl().setForceSubsystem(forces, forces.adoptForce(*this));
}

Force::Gravity::Gravity
   (GeneralForceSubsystem&          forces, 
    const SimbodyMatterSubsystem&   matter,
    const Vec3&                     defGravity)
:   Force(new GravityImpl(matter,UnitVec3(defGravity),defGravity.norm(),0))
{
    SimTK_ERRCHK_ALWAYS(defGravity.norm() > 0,
        "Force::Gravity::ctor(Vec3)",
        "This constructor requires a non-zero Vec3 as the gravity vector"
        " because it has to extract the gravity direction. If you want to"
        " create a Gravity force element for which the default gravity"
        " strength is zero, use the other constructor that allows strength"
        " and direction to be supplied separately.");
    updImpl().setForceSubsystem(forces, forces.adoptForce(*this));
}

Force::Gravity& Force::Gravity::
setDefaultGravityVector(const Vec3& gravity) {
    const Real g = gravity.norm();
    SimTK_ERRCHK_ALWAYS(g > 0,
        "Force::Gravity::setDefaultGravityVector()",
        "This method requires a non-zero Vec3 as the gravity vector"
        " because it has to determine the 'down' direction. If you want to"
        " set the default gravity strength to zero, use setDefaultMagnitude(0)"
        " instead of this method.");
    getImpl().invalidateTopologyCache();
    updImpl().defMagnitude = g;
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

Vec3 Force::Gravity::
getDefaultGravityVector() const 
{   return getImpl().defDirection * getImpl().defMagnitude;}
const UnitVec3& Force::Gravity::
getDefaultDownDirection() const {return getImpl().defDirection;}
Real Force::Gravity::
getDefaultMagnitude() const {return getImpl().defMagnitude;}
Real Force::Gravity::
getDefaultZeroHeight() const {return getImpl().defZeroHeight;}

const Force::Gravity& Force::Gravity::
setGravityVector(State& state, const Vec3& gravity) const {
    const Real g = gravity.norm();
    SimTK_ERRCHK_ALWAYS(g > 0,
        "Force::Gravity::setGravityVector()",
        "This method requires a non-zero Vec3 as the gravity vector"
        " because it has to separate the gravity direction and magnitude."
        " If you want to disable this Gravity force element in this State" 
        " use setMagnitude(0) instead which leaves the direction unchanged but"
        " sets the magnitude to zero.");
    getImpl().updInstanceVars(state).d = UnitVec3(gravity/g, true); 
    getImpl().updInstanceVars(state).g = g; 
    return *this;
}

Vec3 Force::Gravity::
getGravityVector(const State& state) const
{   const GravityImpl::InstanceVars& iv = getImpl().getInstanceVars(state);
    return iv.g*iv.d; }

const Force::Gravity& Force::Gravity::
setDownDirection(State& state, const UnitVec3& down) const {
    SimTK_ERRCHK_ALWAYS(down.isFinite(),
        "Force::Gravity::setDownDirection()",
        "A non-finite 'down' direction was received; did you specify a zero-"
        "length Vec3? The direction must be non-zero.");
    getImpl().updInstanceVars(state).d = down; 
    return *this;
}

const UnitVec3& Force::Gravity::getDownDirection(const State& state) const
{   return getImpl().getInstanceVars(state).d; }

const Force::Gravity& Force::Gravity::
setMagnitude(State& state, Real g) const {
    SimTK_ERRCHK1_ALWAYS(g >= 0,
        "Force::Gravity::setMagnitude()",
        "The gravity magnitude g must be nonnegative but was specified as %g.",
        g);
    getImpl().updInstanceVars(state).g = g; 
    return *this;
}

Real Force::Gravity::getMagnitude(const State& state) const
{   return getImpl().getInstanceVars(state).g; }

const Force::Gravity& Force::Gravity::
setZeroHeight(State& state, Real zeroHeight) const {
    getImpl().updInstanceVars(state).z = zeroHeight; 
    return *this;
}
Real Force::Gravity::getZeroHeight(const State& state) const
{   return getImpl().getInstanceVars(state).z; }

Real Force::Gravity::
getPotentialEnergy(const State& s) const
{   getImpl().ensureForceCacheValid(s); 
    return getImpl().getForceCache(s).pe; }

const SpatialVec& Force::Gravity::
getBodyForce(const State& s, MobilizedBodyIndex mbx) const
{   getImpl().ensureForceCacheValid(s);
    const GravityImpl::ForceCache& fc = getImpl().getForceCache(s);
    return fc.F_GB[mbx]; }

const Vec3& Force::Gravity::
getParticleForce(const State& s, ParticleIndex px) const
{   getImpl().ensureForceCacheValid(s);
    const GravityImpl::ForceCache& fc = getImpl().getForceCache(s);
    return fc.f_GP[px]; }


//------------------------------- GravityImpl ----------------------------------
//------------------------------------------------------------------------------

// This will also calculate potential energy since we can do it on the
// cheap simultaneously with the force. Note that if the strength of gravity
// was set to zero then we already zeroed out the forces and pe during
// realizeInstance() so all we have to do in that case is mark the cache valid
// now.
void Force::GravityImpl::
ensureForceCacheValid(const State& state) const {
    if (isForceCacheValid(state)) return;

    const InstanceVars& iv = getInstanceVars(state);
    if (iv.g == 0) {
        markForceCacheValid(state);
        return;
    }

    const Vec3 gravity      = iv.g * iv.d;
    const Real zeroPEOffset = iv.g * iv.z;
    ForceCache& fc = updForceCache(state);
    fc.pe = 0;

    const int nb = matter.getNumBodies();
    for (MobilizedBodyIndex mbx(1); mbx < nb; ++mbx) {
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

void Force::GravityImpl::
calcForce(const State& state, Vector_<SpatialVec>& bodyForces, 
          Vector_<Vec3>& particleForces, Vector& mobilityForces) const 
{   ensureForceCacheValid(state);
    const ForceCache& fc = getForceCache(state);
    bodyForces     += fc.F_GB;
    particleForces += fc.f_GP; }

// If the force was calculated, then the potential energy will already
// be valid. Otherwise we'll have to calculate it.
Real Force::GravityImpl::
calcPotentialEnergy(const State& state) const 
{   ensureForceCacheValid(state);
    const ForceCache& fc = getForceCache(state);
    return fc.pe; }


} // namespace SimTK

