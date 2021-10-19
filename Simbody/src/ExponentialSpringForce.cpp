/*-----------------------------------------------------------------------------
                               Simbody(tm)
-------------------------------------------------------------------------------
 Copyright (c) 2021 Authors.
 Authors: Frank C. Anderson
 Contributors:

 Licensed under the Apache License, Version 2.0 (the "License"); you may
 not use this file except in compliance with the License. You may obtain a
 copy of the License at http://www.apache.org/licenses/LICENSE-2.0.

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
 ----------------------------------------------------------------------------*/


#include "SimTKcommon.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/ForceSubsystem.h"
#include "simbody/internal/ForceSubsystemGuts.h"
#include "simbody/internal/ExponentialSpringForce.h"

namespace SimTK {

 //=============================================================================
 // Class Declaration for ExponentialSpringForceImpl
 //=============================================================================
class ExponentialSpringForceImpl : public ForceSubsystem::Guts {
public:
    // Constuctor for default parameters.
    //ExponentialSpringForceImpl(const Transform& floor,
    //    const MobilizedBody& body, const Vec3& station, Real mus, Real muk);
    // Constructor for customized parameters.
    ExponentialSpringForceImpl(const Transform& floor,
        const MobilizedBody& body, const Vec3& station, Real mus, Real muk,
        const ExponentialSpringParameters& params);
    Subsystem::Guts* cloneImpl() const override;
    int realizeSubsystemTopologyImpl(State& state) const override;
    int realizeSubsystemDynamicsImpl(const State& state) const override;
    int realizeSubsystemAccelerationImpl(const State& state) const override;
    Real calcPotentialEnergy(const State& state) const override;
    // Static coefficient of friction
    void setMuStatic(State& state, Real mus);
    const Real& getMuStatic(const State& state) const;
    // Kinetic coefficient of friction
    void setMuKinetic(State& state, Real muk) const;
    const Real& getMuKinetic(const State& state) const;
    // Spring Zero
    void resetSprZero(State& state) const;
    Vec3& updSprZero(State& state) const;
    const Vec3& getSprZero(const State& state) const;
    void updSprZeroInCache(const State& state, const Vec3& setpnt) const;
    Vec3 getSprZeroInCache(const State& state) const;
    // Sliding state
    void updSlidingDotInCache(const State& state, Real slidingDot) const;
    const Real getSlidingDotInCache(const State& state) const;
    void realizeSprZeroCache(const State& state) const;
    // Data
    ExponentialSpringData& updData(const State& state) const;
    const ExponentialSpringData& getData(const State& state) const;
    // Parameters
    void setParameters(const ExponentialSpringParameters& params);
    const ExponentialSpringParameters& getParameters() const;
    // STATIC METHODS
    static Real Sigma(Real t0, Real tau, Real t);
    static Real ClampAboveZero(Real value, Real max);
protected:
    ExponentialSpringParameters params;
private:
    Transform contactPlane;
    const MobilizedBody& body;
    Real defaultMus;
    Real defaultMuk;
    Vec3 station;
    Vec3 defaultSprZero;
    ExponentialSpringData defaultData;
    mutable DiscreteVariableIndex indexMus;
    mutable DiscreteVariableIndex indexMuk;
    mutable DiscreteVariableIndex indexSprZero;
    mutable CacheEntryIndex indexSprZeroInCache;
    mutable ZIndex indexZ;
    mutable CacheEntryIndex indexData;
};

} // namespace SimTK


using namespace SimTK;
using std::cout;
using std::endl;

//=============================================================================
// Class - ExponentialSpringData
//=============================================================================
//_____________________________________________________________________________
// Default Constructor
ExponentialSpringData::
ExponentialSpringData() :
    p_G(NaN), v_G(NaN), p(NaN), v(NaN),
    py(NaN), vy(NaN), pxz(NaN), vxz(NaN), pxz_G(NaN),
    fyElas(NaN), fyDamp(NaN), fy(NaN),
    mu(NaN), fxyLimit(NaN),
    fricElas(NaN), fricDamp(NaN), fric(NaN), fxy(NaN),
    f(NaN), f_G(NaN) {}
//_____________________________________________________________________________
// Copy Constructor
ExponentialSpringData::
ExponentialSpringData(const ExponentialSpringData& data) {
    operator=(data);
}
//_____________________________________________________________________________
// Assignment Operator
ExponentialSpringData&
ExponentialSpringData::
operator=(const ExponentialSpringData& data) {
    if(&data!=this) {
        p_G = data.p_G;
        v_G = data.v_G;
        p = data.p;
        v = data.v;
        py = data.py;
        vy = data.vy;
        pxz = data.pxz;
        vxz = data.vxz;
        pxz_G = data.pxz_G;
        fyElas = data.fyElas;
        fyDamp = data.fyDamp;
        fy = data.fy;
        mu = data.mu;
        fxyLimit = data.fxyLimit;
        fricElas = data.fricElas;
        fricDamp = data.fricDamp;
        fric = data.fric;
        fxy = data.fxy;
        f = data.f;
        f_G = data.f_G;
    }
    return *this;
}


//_____________________________________________________________________________
// Constructor
ExponentialSpringForceImpl::
ExponentialSpringForceImpl(const Transform& floor,
    const MobilizedBody& body, const Vec3& station, Real mus, Real muk,
    const ExponentialSpringParameters& params) :
    ForceSubsystem::Guts("ExponentialSpringForce", "0.0.1"),
    contactPlane(floor), body(body), station(station),
    defaultMus(mus), defaultMuk(muk), defaultSprZero(Vec3(0., 0., 0.))
{
    // Check for valid static coefficient
    if(defaultMus < 0.0) defaultMus = 0.0;

    // Check for valid kinetic coefficient
    if(defaultMuk < 0.0) defaultMuk = 0.0;
    if(defaultMuk > defaultMus) defaultMuk = defaultMus;

    // Assign the parameters
    this->params = params;
}
//_____________________________________________________________________________
// Clone
Subsystem::Guts*
ExponentialSpringForceImpl::
cloneImpl() const {
    return new ExponentialSpringForceImpl(contactPlane, body, station,
        defaultMus, defaultMuk, params);
}
//_____________________________________________________________________________
// Realize the system at the Topology Stage.
// This is where state variables and their associated indices are most often
// created, the Model stage being the other possibility for things like
// Q's (i.e., quaternions need 4 Q's).
int
ExponentialSpringForceImpl::
realizeSubsystemTopologyImpl(State& state) const {
    // Coefficients of friction: mus and muk
    // Both are treated as discrete states.  By doing so, mus and muk
    // can be changed during a simulation.
    indexMus = allocateDiscreteVariable(state,
        Stage::Dynamics, new Value<Real>(defaultMus));
    indexMuk = allocateDiscreteVariable(state,
        Stage::Dynamics, new Value<Real>(defaultMuk));

    // SprZero
    // The spring SprZero is a discrete variable that is auto updated during
    // the course of a simulation so that the frictional force generated
    // by a spring is consistent with the spring's theoretical limit
    // (mu*Fnormal).
    // Because the SprZero could potentially need updating after every
    // integration step, it is treated as an AutoUpdate Discrete Variable.

    // Index to the actual discrete variable
    // Changing the SprZero will invalidate the Dynamics Stage.
    // The SprZero held in Cache depends on realization through the Velocity
    // Stage.
    indexSprZero =
        allocateAutoUpdateDiscreteVariable(state, Stage::Dynamics,
            new Value<Vec3>(defaultSprZero),Stage::Velocity);
    // Index to the update value that is held in the Cache
    indexSprZeroInCache =
        getDiscreteVarUpdateIndex(state, indexSprZero);

    // Sliding
    // Sliding is a continuous state variable used to model the transition
    // of a spring between fixed and sliding.
    //        Sliding = 0        means fully fixed in place (lower bound)
    //        Sliding = 1        means fully sliding (upper bound)
    // Because Sliding is a continuous state variable (a Z variable in
    // Simbody terminology), changes in Sliding are made only by setting
    // the value of its time derivative, SlidingDot, which is updated in
    // realizeSubsystemAccelerationImpl().
    Real initialValue = 0.0;
    Vector zInit(1, initialValue);
    indexZ = allocateZ(state, zInit);

    // Data
    // Useful information that is computed during a simulation is organized
    // in class ExponentialSpringData and stored as a Cache Entry.
    indexData = allocateCacheEntry(state, Stage::Dynamics,
        new Value<ExponentialSpringData>(defaultData));

    return 0;
}
//_____________________________________________________________________________
// Realize compuations at the Dynamics Stage.  In other words, compute the
// forces modeled by this Subsystem.
//
// In the code below...
//
// "params" references the configurable parameters that govern the behavior
// of the exponential spring.  These can be changed by the user, but the
// System must be realized at the Topology Stage after any such change.
//
// "data" references the key data that are calculated and stored as a
// Cache Entry when the System is realized at the Dynamics Stage.
// These data can be retrieved during a simulation by a reporter or handler,
// for example.
// 
// Variables without a suffix are expressed in the floor frame.
// 
// Variables with the _G suffix are expressed in the ground frame.
//
// Most everything important happens in this one method.
int
ExponentialSpringForceImpl::
realizeSubsystemDynamicsImpl(const State& state) const {
    // Get current accumulated forces
    const MultibodySystem& system = MultibodySystem::downcast(getSystem());
    const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();
    Vector_<SpatialVec>& forces_G =
        system.updRigidBodyForces(state, Stage::Dynamics);

    // Retrieve a writable reference to the data cache entry.
    // Most computed quantities are stored in the data cache.
    ExponentialSpringData& data = updData(state);

    // Get position and velocity of the spring station in the ground frame
    data.p_G = body.findStationLocationInGround(state, station);
    data.v_G = body.findStationVelocityInGround(state, station);

    // Transform the position and velocity into the floor frame.
    data.p = contactPlane.shiftBaseStationToFrame(data.p_G);
    data.v = contactPlane.xformBaseVecToFrame(data.v_G);

    // Resolve into normal (y) and tangential parts (xz plane)
    // Normal (perpendicular to Floor)
    data.py = data.p[1];
    data.vy = data.v[1];
    // Tangent (in plane of Floor)
    data.pxz = data.p;    data.pxz[1] = 0.0;
    data.vxz = data.v;    data.vxz[1] = 0.0;
    // Not used to calculate force, but likely useful for visualization
    data.pxz_G = contactPlane.shiftFrameStationToBase(data.pxz);

    // Get all the parameters upfront
    Real d0, d1, d2;
    params.getShapeParameters(d0, d1, d2);
    Real kvNorm = params.getNormalViscosity();
    Real kpFric = params.getElasticity();
    Real kvFric = params.getViscosity();
    Real kTau = 1.0 / params.getSlidingTimeConstant();
    Real vSettle = params.getSettleVelocity();

    // Normal Force (perpendicular to floor) --------------------------------
    // Elastic Part
    data.fyElas = d1*std::exp(-d2*(data.py-d0));
    // Damping Part
    data.fyDamp = kvNorm * data.vy * data.fyElas;
    // Total
    data.fy = data.fyElas-data.fyDamp;
    // Don't allow the normal force to be negative or too large.
    data.fy = ClampAboveZero(data.fy, 100000.0);
    //if (data.fy < 0.0) data.fy = 0.0;

    // Friction (in the plane of floor) -------------------------------------
    // Get the sliding state.
    Real sliding = getZ(state)[indexZ];
    // Compute the maximum allowed frictional force based on the current
    // coefficient of friction.
    Real mus = getMuStatic(state);
    Real muk = getMuKinetic(state);
    data.mu = mus - sliding * (mus - muk);
    data.fxyLimit = data.mu * data.fy;
    // Access the SprZero from the State.
    Vec3 p0 = getSprZero(state);
    // The SprZero is always expressed in the Floor frame, so its y-component
    // should always be zero. The following statement shouldn't be necessary,
    // but rounding is possible
    p0[1] = 0.0;
    // Elastic part
    Vec3 r = data.pxz - p0;
    data.fricElas = -kpFric * r;
    Real fxyElas = data.fricElas.norm();
    // Viscous part (damping)
    data.fricDamp = -kvFric * data.vxz;
    // Total
    data.fric = data.fricElas + data.fricDamp;
    data.fxy = data.fric.norm();
    bool limitReached = false;
    if (data.fxy > data.fxyLimit) {
        data.fxy = data.fxyLimit;
        data.fric = data.fxy * data.fric.normalize();
        limitReached = true;
    }
    // If the spring is stretched beyond its limit, update the spring zero.
    // Note that no discontinuities in the friction force are introduced.
    // The spring zero is just made to be consistent with the limiting
    // frictional force.
    Vec3 p0New, fricElasNew;
    if (fxyElas > data.fxyLimit) {
        // Compute a new spring zero.
        fxyElas = data.fxyLimit;
        fricElasNew = fxyElas * data.fricElas.normalize();
        p0New = data.pxz + fricElasNew / kpFric;
        // Make sure that p0 is always in the plane of the floor.
        p0New[1] = 0.0;
        // Update the spring zero cache and mark the cache as realized.
        // Only place that the following two lines are called.
        updSprZeroInCache(state, p0New);
        markCacheValueRealized(state, indexSprZeroInCache);
    }

    // Update SlidingDot
    Real vMag = data.vxz.norm();
    Real slidingDot = 0.0;
    if (limitReached)  slidingDot = kTau*(1.0 - sliding);
    else if (vMag < vSettle)  slidingDot = -kTau*sliding;
    updSlidingDotInCache(state, slidingDot);

    // Total spring force expressed in the floor frame
    data.f = data.fric;        // The x and z components are friction.
    data.f[1] = data.fy;    // The y component is the normal force.

    // Transform the spring forces back to the Ground frame
    data.f_G = contactPlane.xformFrameVecToBase(data.f);

    // Apply the force
    body.applyForceToBodyPoint(state, station, data.f_G, forces_G);

    return 0;
}
//_____________________________________________________________________________
// Realize compuations at the Acceleration Stage,
// which means compute and update the time derivative of each continuous
// variable (Q, U, and Z).
// This spring only has one Z variable -- the Sliding state.
// Because all of the quantities needed for computing SlidingDot are
// available in realizeSubsystemDynamicsImpl(), SlidingDot is updated
// there.  This method does nothing.
int
ExponentialSpringForceImpl::
realizeSubsystemAccelerationImpl(const State& state) const {
    //Real sliding = getZ(state)[indexZ];
    //Real slidingDot = -sliding / 0.01;
    //updSlidingDotInCache(state, 0.0);
    return 0;
}
//_____________________________________________________________________________
// Calculate the potential energy stored in the spring.
// The System should be realized through Stage::Dynamics before a call to
// this method is made; an exception will be thrown otherwise.
Real
ExponentialSpringForceImpl::
calcPotentialEnergy(const State& state) const {
    const MultibodySystem& system = MultibodySystem::downcast(getSystem());
    const ExponentialSpringData& data = getData(state);

    // Strain energy in the normal direction (exponential spring)
    Real d0, d1, d2;
    params.getShapeParameters(d0, d1, d2);
    double energy = data.fyElas / d2;

    // Strain energy in the tangent plane (friction spring)
    Vec3 p0 = getSprZero(state);
    Vec3 r = data.pxz - p0;
    energy += 0.5 * params.getElasticity() * r.norm() * r.norm();

    return energy;
}

//_____________________________________________________________________________
// Update the coefficient of static friction for this spring.
// The specified value must obey the following constraints:
//
//        0.0 <= muk <= mus
//
// If mus is less than 0.0, mus is set to 0.0.
// If mus is less than mus, muk is set to mus.
void
ExponentialSpringForceImpl::
setMuStatic(State& state, Real mus) {
    // Keep mus greter than or equal to 0.0.
    if (mus < 0.0) mus = 0.0;
    Value<Real>::updDowncast(updDiscreteVariable(state, indexMus)) = mus;

    // Make sure muk is less than or equal to mus
    Real muk = getMuKinetic(state);
    if (muk > mus) {
        muk = mus;
        Value<Real>::updDowncast(updDiscreteVariable(state, indexMuk)) = muk;
    }
}
//_____________________________________________________________________________
// Get the coefficient of static friction for this spring.
const Real&
ExponentialSpringForceImpl::
getMuStatic(const State& state) const {
    return Value<Real>::downcast(getDiscreteVariable(state, indexMus));
}

//_____________________________________________________________________________
// Update the coefficient of kinetic friction for this spring.
// The specified value must obey the following constraints:
//
//       0.0 <= muk <= mus
//
// If muk is less than 0.0, muk is set to 0.0.
// If muk is greater than mus, mus is set to muk.
void
ExponentialSpringForceImpl::
setMuKinetic(State& state, Real muk) const {
    // Keep muk >= to zero.
    if (muk < 0.0) muk = 0.0;
    Value<Real>::updDowncast(updDiscreteVariable(state, indexMuk)) = muk;

    // Make sure mus is greater than or equal to muk
    Real mus = getMuStatic(state);
    if (muk > mus) {
        Value<Real>::updDowncast(updDiscreteVariable(state, indexMus)) = muk;
    }
}
//_____________________________________________________________________________
// Get the coefficient of kinetic friction for this spring.
const Real&
ExponentialSpringForceImpl::
getMuKinetic(const State& state) const {
    return Value<Real>::downcast(getDiscreteVariable(state, indexMuk));
}


//_____________________________________________________________________________
// Reset the zero of this spring.
//
// This method sets the spring zero to the point on the Floor that coincides
// with the Station that has been specified on the MobilizedBody.  In this
// process, the System is realized through the Position Stage.
void
ExponentialSpringForceImpl::
resetSprZero(State& state) const {
    // Realize through to the Position Stage
    const MultibodySystem& system = MultibodySystem::downcast(getSystem());
    system.realize(state, Stage::Position);

    // Get position of the spring station in the Ground frame
    Vec3 p_G = body.findStationLocationInGround(state, station);

    // Transform the position to the Floor frame.
    Vec3 p_F = contactPlane.shiftBaseStationToFrame(p_G);

    // Project into the plane of the Floor
    p_F[1] = 0.0;

    // Update the spring zero
    updSprZero(state) = p_F;
}
//_____________________________________________________________________________
// Update the SprZero of this exponential spring.
// Updating the SprZero will invalidate the state at the Dynamics Stage.
//
// The SprZero is the current zero of the spring in the floor plane.  It is
// the basis for computing the position-dependent portion of the frictional
// force exerted by the spring.
// The SprZero should be updated whenever the elastic part of the frictional
// force exeeds its theoretical limit (i.e., mu*Fnormal).  The new SprZero
// should be such that the elastic part generates the limit, not more or less.
// The updated value should not account for the size of the damping force.
Vec3&
ExponentialSpringForceImpl::
updSprZero(State& state) const {
    return Value<Vec3>::updDowncast(updDiscreteVariable(state, indexSprZero));
}
//_____________________________________________________________________________
// Get the SprZero of this exponential spring.
// The SprZero is the current zero of the spring in the ground plane.  It is
// the basis for computing the position-dependent portion of the frictional
// force exerted by the spring.
const Vec3&
ExponentialSpringForceImpl::
getSprZero(const State& state) const {
    return Value<Vec3>::downcast(getDiscreteVariable(state, indexSprZero));
}
//_____________________________________________________________________________
// Update the SprZero held in the Cache of this exponential spring.
// Updating the SprZero held in the Cache will not invalidate the State.
// Use this method.
void
ExponentialSpringForceImpl::
updSprZeroInCache(const State& state, const Vec3& setpoint) const {
    Value<Vec3>::updDowncast(updDiscreteVarUpdateValue(state, indexSprZero))
        = setpoint;
}
//_____________________________________________________________________________
// Get the SprZero of this exponential spring.
// The SprZero is the current zero of the spring in the ground plane.  It is
// the basis for computing the position-dependent portion of the frictional
// force exerted by the spring.
Vec3
ExponentialSpringForceImpl::
getSprZeroInCache(const State& state) const {
    return Value<Vec3>::downcast(
        getDiscreteVarUpdateValue(state,indexSprZero) );
}
//_____________________________________________________________________________
// Realize the SprZero Cache.
// There needs to be some assurance that the initial value of the SprZero
// is valid.  Therefore, the first time a getSprZeroInCache() is called
// a calle to realizeSprZeroCache() is also made.  Once the Cache for
// the SprZero has been realized, it may be repeatedly accessed with only
// the cost of the method call and the if statement.
//
// Not really using this method, but it might be useful to keep around
// in case I want a reminder of how cache access works.
void
ExponentialSpringForceImpl::
realizeSprZeroCache(const State& state) const {
    if (isCacheValueRealized(state, indexSprZeroInCache)) return;
    Vec3 sprZero = getSprZero(state);
    Real time = state.getTime();
    sprZero[0] = 0.01 * time;
    sprZero[1] = 0.0;
    sprZero[2] = 0.01 * time;
    updSprZeroInCache(state, sprZero);
    markCacheValueRealized(state, indexSprZeroInCache);
}
//_____________________________________________________________________________
// Update SlidingDot in the Cache of this exponential spring.
// Updating the value of SlidingDot in the Cache will not invalidate the State.
void
ExponentialSpringForceImpl::
updSlidingDotInCache(const State& state, Real slidingDot) const {
    updZDot(state)[indexZ] = slidingDot;
}
//_____________________________________________________________________________
// Get value of SlidingDot that is held in cashe.
// SlidingDot is the time derivative of the Sliding state variable.
const Real
ExponentialSpringForceImpl::
getSlidingDotInCache(const State& state) const {
    return getZDot(state)[indexZ];
}

//-----------------------------------------------------------------------------
// DATA
//-----------------------------------------------------------------------------
//_____________________________________________________________________________
// Retrieve a writable reference to the spring data contained in the Cache.
ExponentialSpringData&
ExponentialSpringForceImpl::
updData(const State& state) const {
    return Value<ExponentialSpringData>::updDowncast(updCacheEntry(state,
        indexData));
}
//_____________________________________________________________________________
// Retrieve a const reference to the spring data contained in the Cache.
const ExponentialSpringData&
ExponentialSpringForceImpl::
getData(const State& state) const {
    return Value<ExponentialSpringData>::downcast(getCacheEntry(state,
        indexData));
}


//-----------------------------------------------------------------------------
// PARAMETERS
//-----------------------------------------------------------------------------
//_____________________________________________________________________________
// Update the parameters for this exponential spring.
void
ExponentialSpringForceImpl::
setParameters(const ExponentialSpringParameters& params) {
    this->params = params;
    invalidateSubsystemTopologyCache();
}
//_____________________________________________________________________________
// Get the parameters for this exponential spring.
const ExponentialSpringParameters&
ExponentialSpringForceImpl::
getParameters() const {
    return params;
}

//_____________________________________________________________________________
// Sigma - a function that transitions smoothly from 0.0 to 1.0 or
// from 1.0 to 0.0.
//
//   f(t) = 1.0 / {1.0 + exp[(t - t0) / tau]}
//   t0 - time about which the transition is centered.F(t0) = 0.5.
//   tau - time constant modifying the rate of the transiton occurs.
//   tau < 0.0 generates a step up
//   tau > 0.0 generates a step down
//   A larger value of tau results in a more gradual transition.
//
// Step Up(negative tau)
//                    | f(t)
//                   1.0                         * ************
//                    |                  *
//                    |              *
//                   0.5 +
//                    |         *
//                    |    *
//  ***************---|-----------t0-------------------------  t
//                    |
//
// Step Down(positive tau)
//                    | f(t)
//  ***************  1.0
//                    |    *
//                    |         *
//                   0.5 +
//                    |               *
//                    |                   *
//  ------------------|-----------t0------------*************  t
//
Real
ExponentialSpringForceImpl::
Sigma(Real t0, Real tau, Real t) {
    Real x = (t - t0) / tau;
    Real s = 1.0 / (1.0 + std::exp(x));
    return s;
}
//_____________________________________________________________________________
// Clamp a value between zero and a maximum value.
Real
ExponentialSpringForceImpl::
ClampAboveZero(Real value, Real max) {
    if (value > max) return max;
    if (value < 0.0) return 0.0;
    return value;
}


//=============================================================================
// Class - ExponentialSpringForce
//=============================================================================
//_____________________________________________________________________________
// Constructor for default spring parameters.
//ExponentialSpringForce::
//ExponentialSpringForce(MultibodySystem& system,
//    const Transform& contactPlane,
//    const MobilizedBody& body,const Vec3& station,
//    Real mus, Real muk)
//{
//    adoptSubsystemGuts(
//        new ExponentialSpringForceImpl(contactPlane, body, station,
//            mus, muk));
//    system.addForceSubsystem(*this);
//}

//_____________________________________________________________________________
// Constructor.
ExponentialSpringForce::
ExponentialSpringForce(MultibodySystem& system,
    const Transform& contactPlane,
    const MobilizedBody& body, const Vec3& station,
    Real mus, Real muk, ExponentialSpringParameters params)
{
    adoptSubsystemGuts(
        new ExponentialSpringForceImpl(contactPlane, body, station,
            mus, muk, params));
    system.addForceSubsystem(*this);
}
//_____________________________________________________________________________
// Set new parameters for this exponential spring.
//
// Note that the underlying implementation (ExponentialSpringForceImpl) owns
// its own ExponentialSpringParameters instance.  When this method is called,
// the underlying implementation sets its parameters equal to the parameters
// sent in through the argument list by calling the assignment operator:
//
//        ExponentialSpringForceImple::params = params
// 
// @see ExponentialSpringParameters for the list of parameters.
void
ExponentialSpringForce::
setParameters(const ExponentialSpringParameters& params) {
    updImpl().setParameters(params);
}
//_____________________________________________________________________________
// Get the current parameters for this exponential spring.
// @see ExponentialSpringParameters for the list of parameters.
const ExponentialSpringParameters&
ExponentialSpringForce::
getParameters() const {
    return getImpl().getParameters();
}

//_____________________________________________________________________________
// Set the static coefficient of fricition
void
ExponentialSpringForce::
setMuStatic(State& state, Real mus) {
    updImpl().setMuStatic(state, mus);
}
//_____________________________________________________________________________
// Get the static coefficient of fricition
const Real&
ExponentialSpringForce::
getMuStatic(const State& state) const {
    return getImpl().getMuStatic(state);
}

//_____________________________________________________________________________
// Set the kinetic coefficient of fricition
void
ExponentialSpringForce::
setMuKinetic(State& state, Real muk) {
    updImpl().setMuKinetic(state, muk);
}
//_____________________________________________________________________________
// Get the kinetic coefficient of fricition
const Real&
ExponentialSpringForce::
getMuKinetic(const State& state) const {
    return getImpl().getMuKinetic(state);
}

//_____________________________________________________________________________
// Reset the spring zero.
// This method sets the spring zero to the point on the contact plane that
// coincides with the Station that has been specified on the MobilizedBody
// for which this exponential spring was constructed.
void
ExponentialSpringForce::
resetSpringZero(State& state) const {
    getImpl().resetSprZero(state);
}

//_____________________________________________________________________________
// Retrieve access to a const reference to the data cache.  The System
// should be realized through Stage::Dynamics before a call to this method;
// an exception will be thrown otherwise.
const ExponentialSpringData&
ExponentialSpringForce::
getData(const State& state) const {
    return getImpl().getData(state);
}

//_____________________________________________________________________________
// Get a reference to the underlying implementation that will allow changes
// to be made to underlying parameters and states.
ExponentialSpringForceImpl&
ExponentialSpringForce::
updImpl() {
    return dynamic_cast<ExponentialSpringForceImpl&>(updRep());
}
//_____________________________________________________________________________
// Get a reference to the underlying implementation that will allow
// access (but not change) to underlying parameters and states.
const ExponentialSpringForceImpl&
ExponentialSpringForce::
getImpl() const {
    return dynamic_cast<const ExponentialSpringForceImpl&>(getRep());
}


