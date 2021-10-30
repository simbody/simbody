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
// Struct ExponentialSpringData
//=============================================================================
/** ExponentialSpringData is an internal data structure used by the
implementation class ExponentialSpringForceImpl to store and retrieve
important quantities kept in the State's data cache.

Note - As originally coded, users could access class ExponentialSpringData
directly via a const reference. To make the exponential spring code less
brittle, however, class ExponentialSpringData was migrated to being an
internal data structure that is accessed directly only by class
ExponentialSpringForceImpl. The data stored in ExponentialSpringData can
now only be accessed by the end user via conventional assessor methods
in class ExponentialSpringForce.

Even though ExponentialSpringData is now used only internally, the
Doxygen-compliant comments originally written for the class might be helpful
to developers. I therefore kept the comments. They follow below:

ExponentialSpringData is a helper class that is used to store key
quantities associated with the ExponentialSpringForce Subsystem during a
simulation. An instance of this class serves as the data Cache Entry for the
ExponentialSpringForceImpl Subsystem. All of its member variables are
guaranteed to be calculated and set once the System has been realized
to Stage::Dynamics.

To understand what the quantities organized by this class represent, a basic
description of the contact problem that is solved, along with a description
of coordinate frame conventions, will be helpful.

Class ExponentialSpringForce computes and applies a contact force at a
specified point on a MobilizedBody (i.e., a Station) due to interaction of
that point with a specified contact plane. That plane is typically used to
model interactions with a floor, but need not be limited to this use case.
The contact plane can be rotated and displaced relative to the ground frame
and so can be used to model a wall or ramp, for example.

%Contact force computations are carried out in the frame of the contact plane.
The positive y-axis of the contact frame defines the normal of the
contact plane. The positive y-axis is the axis along which a repelling
normal force (modeled using an exponential) is applied. The x-axis and
z-axis of the contact frame are tangent to the contact plane. The friction
force will always lie in x-z plane.

Member variables with a "y" suffix (e.g., py, vy, or fy) indicate that
these quantities are directed normal to the contact plane.  Member varaibles
with an "xz" suffix (e.g., pxz, vxz, or fxz) indicate that these quantities
lie in the contact plane (or tangent to it) and are associated with the
friction force.

Member variables with a "_G" suffix are expressed in the Ground frame. Member
variables without a "_G" suffix are expressed the contact plane frame. */
struct ExponentialSpringData {
    /** Position of the body spring station in the ground frame. */
    Vec3 p_G;
    /** Velocity of the body spring station in the ground frame. */
    Vec3 v_G;
    /** Position of the body spring station in the frame of the contact
    plane. */
    Vec3 p;
    /** Velocity of the body spring station in the frame of the contact
    plane. */
    Vec3 v;
    /** Displacement of the body spring station normal to the floor expressed
    in the frame of the contact plane. */
    Real py;
    /** Velocity of the body spring station normal to the contact plane
    expressed in the frame of the contact plane. */
    Real vy;
    /** Position of the body spring station projected onto the contact plane
    expressed in the frame of the contact plane. */
    Vec3 pxz;
    /** Velocity of the body spring station in the contact plane expressed in
    the frame of the contact plane. */
    Vec3 vxz;
    /** Elastic force in the normal direction. */
    Real fyElas;
    /** Damping force in the normal direction. */
    Real fyDamp;
    /** Total normal force expressed in the frame of the contact plane. */
    Real fy;
    /** Instantaneous coefficient of friction. */
    Real mu;
    /** Limit of the frictional force. */
    Real fxzLimit;
    /** Elastic frictional force expressed in the frame of the contact plane.*/
    Vec3 fricElas;
    /** Damping frictional force expressed in the frame of the contact plane.*/
    Vec3 fricDamp;
    /** Total frictional force (elastic + damping) expressed in the frame of
    the contact plane. */
    Vec3 fric;
    /** Magnitude of the frictional force. */
    Real fxz;
    /** Resultant spring force (normal + friction) expressed in the floor
    frame. */
    Vec3 f;
    /** Resultant spring force (normal + friction) expressed in the ground
    frame. */
    Vec3 f_G;
};


 //============================================================================
 // Class ExponentialSpringForceImpl
 //============================================================================
class ExponentialSpringForceImpl : public ForceSubsystem::Guts {
public:
// Constructor
ExponentialSpringForceImpl(const Transform& floor,
const MobilizedBody& body, const Vec3& station, Real mus, Real muk,
const ExponentialSpringParameters& params) :
ForceSubsystem::Guts("ExponentialSpringForce", "0.0.1"),
contactPlane(floor), body(body), station(station),
defaultMus(mus), defaultMuk(muk), defaultSprZero(Vec3(0., 0., 0.)) {
    // Check for valid static coefficient
    if(defaultMus < 0.0) defaultMus = 0.0;
    // Check for valid kinetic coefficient
    if(defaultMuk < 0.0) defaultMuk = 0.0;
    if(defaultMuk > defaultMus) defaultMuk = defaultMus;
    // Assign the parameters
    this->params = params;
}

//-----------------------------------------------------------------------------
// Accessors
//-----------------------------------------------------------------------------
// SIMPLE
const Transform& getContactPlane() const { return contactPlane; }
const MobilizedBody& getBody() const { return body; }
const Vec3& getStation() const { return station; }

// TOPOLOGY PARAMETERS
const ExponentialSpringParameters& getParameters() const { return params; }
void setParameters(const ExponentialSpringParameters& params) {
    this->params = params;
    invalidateSubsystemTopologyCache(); }

// DATA CACHE
ExponentialSpringData& updData(const State& state) const {
    return Value<ExponentialSpringData>::updDowncast(updCacheEntry(state,
        indexData)); }
const ExponentialSpringData& getData(const State& state) const {
    return Value<ExponentialSpringData>::downcast(getCacheEntry(state,
        indexData)); }

// SLIDING STATE
Real getSliding(const State& state) const {
    return getZ(state)[indexZ]; }
Real getSlidingDotInCache(const State& state) const {
    return getZDot(state)[indexZ]; }
void updSlidingDotInCache(const State& state, Real slidingDot) const {
    // Does not invalidate the State.
    updZDot(state)[indexZ] = slidingDot; }

// SPRING ZERO
const Vec3& getSprZero(const State& state) const {
    return Value<Vec3>::downcast(getDiscreteVariable(state, indexSprZero)); }
Vec3& updSprZero(State& state) const {
    // Update occurs when the elastic force exceeds mu*N.
    return Value<Vec3>::updDowncast(updDiscreteVariable(state,indexSprZero)); }
Vec3 getSprZeroInCache(const State& state) const {
    return Value<Vec3>::downcast(
        getDiscreteVarUpdateValue(state, indexSprZero));
}void updSprZeroInCache(const State& state, const Vec3& setpoint) const {
    // Will not invalidate the State.
    Value<Vec3>::updDowncast(updDiscreteVarUpdateValue(state, indexSprZero))
        = setpoint; }

// STATIC COEFFICENT OF FRICTION
const Real& getMuStatic(const State& state) const {
    return Value<Real>::downcast(getDiscreteVariable(state, indexMus)); }
void setMuStatic(State& state, Real mus) {
    // Keep mus greter than or equal to 0.0.
    if(mus < 0.0) mus = 0.0;
    Value<Real>::updDowncast(updDiscreteVariable(state, indexMus)) = mus;
    // Make sure muk is less than or equal to mus
    Real muk = getMuKinetic(state);
    if(muk > mus) {
        muk = mus;
        Value<Real>::updDowncast(updDiscreteVariable(state, indexMuk)) = muk;
    }
}

// KINETIC COEFFICENT OF FRICTION
const Real& getMuKinetic(const State& state) const {
    return Value<Real>::downcast(getDiscreteVariable(state, indexMuk)); }
void setMuKinetic(State& state, Real muk) const {
    // Keep muk >= to zero.
    if(muk < 0.0) muk = 0.0;
    Value<Real>::updDowncast(updDiscreteVariable(state, indexMuk)) = muk;
    // Make sure mus is greater than or equal to muk
    Real mus = getMuStatic(state);
    if(muk > mus) {
        Value<Real>::updDowncast(updDiscreteVariable(state, indexMus)) = muk;
    }
}

//-----------------------------------------------------------------------------
// ForceSubsystem Methods (overrides of virtual methods)
//-----------------------------------------------------------------------------
//_____________________________________________________________________________
// Clone
Subsystem::Guts*
cloneImpl() const override {
    return new ExponentialSpringForceImpl(contactPlane, body, station,
        defaultMus, defaultMuk, params);
}
//_____________________________________________________________________________
// Topology - allocate state variables and the data cache.
int
realizeSubsystemTopologyImpl(State& state) const override {
    // Coefficients of friction: mus and muk
    indexMus = allocateDiscreteVariable(state,
        Stage::Dynamics, new Value<Real>(defaultMus));
    indexMuk = allocateDiscreteVariable(state,
        Stage::Dynamics, new Value<Real>(defaultMuk));
    // SprZero
    indexSprZero =
        allocateAutoUpdateDiscreteVariable(state, Stage::Dynamics,
            new Value<Vec3>(defaultSprZero), Stage::Velocity);
    indexSprZeroInCache =
        getDiscreteVarUpdateIndex(state, indexSprZero);
    // Sliding
    Real initialValue = 0.0;
    Vector zInit(1, initialValue);
    indexZ = allocateZ(state, zInit);
    // Data
    indexData = allocateCacheEntry(state, Stage::Dynamics,
        new Value<ExponentialSpringData>(defaultData));
    return 0;
}
//_____________________________________________________________________________
// Dynamics - compute the forces modeled by this Subsystem.
//
// "params" references the configurable parameters that govern the behavior
// of the exponential spring. These can be changed by the user, but the
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
// Most every calculation happens in this one method.
int
realizeSubsystemDynamicsImpl(const State& state) const override {
    // Get current accumulated forces
    const MultibodySystem& system = MultibodySystem::downcast(getSystem());
    const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();
    Vector_<SpatialVec>& forces_G =
        system.updRigidBodyForces(state, Stage::Dynamics);

    // Retrieve a writable reference to the data cache entry.
    // Most computed quantities are stored in the data cache.
    ExponentialSpringData& data = updData(state);

    // Get position and velocity of the spring station in Ground
    data.p_G = body.findStationLocationInGround(state, station);
    data.v_G = body.findStationVelocityInGround(state, station);

    // Transform the position and velocity into the contact frame.
    data.p = contactPlane.shiftBaseStationToFrame(data.p_G);
    data.v = contactPlane.xformBaseVecToFrame(data.v_G);

    // Resolve into normal (y) and tangential parts (xz plane)
    // Normal (perpendicular to Floor)
    data.py = data.p[1];
    data.vy = data.v[1];
    // Tangent (in plane of Floor)
    data.pxz = data.p;    data.pxz[1] = 0.0;
    data.vxz = data.v;    data.vxz[1] = 0.0;

    // Get all the parameters upfront
    Real d0, d1, d2;
    params.getShapeParameters(d0, d1, d2);
    Real kvNorm = params.getNormalViscosity();
    Real kpFric = params.getElasticity();
    Real kvFric = params.getViscosity();
    Real kTau = 1.0 / params.getSlidingTimeConstant();
    Real vSettle = params.getSettleVelocity();

    // Normal Force (perpendicular to contact plane) -------------------------
    // Elastic Part
    data.fyElas = d1 * std::exp(-d2 * (data.py - d0));
    // Damping Part
    data.fyDamp = -kvNorm * data.vy * data.fyElas;
    // Total
    data.fy = data.fyElas + data.fyDamp;
    // Don't allow the normal force to be negative or too large.
    // Note that conservation of energy will fail if bounds are enforced.
    data.fy = ClampAboveZero(data.fy, 1000000.0);

    // Friction (in the plane of contact plane) ------------------------------
    // Get the sliding state.
    Real sliding = getZ(state)[indexZ];
    // Compute the maximum allowed frictional force based on the current
    // coefficient of friction.
    Real mus = getMuStatic(state);
    Real muk = getMuKinetic(state);
    data.mu = mus - sliding * (mus - muk);
    data.fxzLimit = data.mu * data.fy;
    // Access the SprZero from the State.
    Vec3 p0 = getSprZero(state);
    // The SprZero is always expressed in the Floor frame, so its y-component
    // should always be zero. The following statement shouldn't be necessary,
    // but rounding is possible
    p0[1] = 0.0;
    // Elastic part
    Vec3 r = data.pxz - p0;
    data.fricElas = -kpFric * r;
    Real fxzElas = data.fricElas.norm();
    // Viscous part (damping)
    data.fricDamp = -kvFric * data.vxz;
    // Total
    data.fric = data.fricElas + data.fricDamp;
    data.fxz = data.fric.norm();
    bool limitReached = false;
    if(data.fxz > data.fxzLimit) {
        data.fxz = data.fxzLimit;
        data.fric = data.fxz * data.fric.normalize();
        limitReached = true;
    }
    // If the spring is stretched beyond its limit, update the spring zero.
    // Note that no discontinuities in the friction force are introduced.
    // The spring zero is just made to be consistent with the limiting
    // frictional force.
    Vec3 p0New, fricElasNew;
    if(fxzElas > data.fxzLimit) {
        // Compute a new spring zero.
        fxzElas = data.fxzLimit;
        fricElasNew = fxzElas * data.fricElas.normalize();
        p0New = data.pxz + fricElasNew / kpFric;
        // Make sure that p0 is always in the contact plane.
        p0New[1] = 0.0;
        // Update the spring zero cache and mark the cache as realized.
        // Only place that the following two lines are called.
        updSprZeroInCache(state, p0New);
        markCacheValueRealized(state, indexSprZeroInCache);
    }

    // Update SlidingDot
    Real vMag = data.vxz.norm();
    Real slidingDot = 0.0;
    if(limitReached)  slidingDot = kTau * (1.0 - sliding);
    else if(vMag < vSettle)  slidingDot = -kTau * sliding;
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
// Acceleration - compute and update the derivatives of continuous states.
// Because all of the quantities needed for computing SlidingDot are
// available in realizeSubsystemDynamicsImpl(), SlidingDot is updated
// there. This method at the moment does nothing.
int
realizeSubsystemAccelerationImpl(const State& state) const override {
    //Real sliding = getZ(state)[indexZ];
    //Real slidingDot = -sliding / 0.01;
    //updSlidingDotInCache(state, 0.0);
    return 0;
}
//_____________________________________________________________________________
// Potential Energy - calculate the potential energy stored in the spring.
// The System should be realized through Stage::Dynamics before a call to
// this method is made.
Real
calcPotentialEnergy(const State& state) const override {
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

//-----------------------------------------------------------------------------
// Utility and Static Methods
//-----------------------------------------------------------------------------
//_____________________________________________________________________________
// Project the body spring station onto the contact plane.
void
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
// Note - Not really using this method, but keeping it around in case I want
// a reminder of how cache access works.
// Realize the SprZero Cache.
// There needs to be some assurance that the initial value of the SprZero
// is valid. Therefore, the first time a getSprZeroInCache() is called
// a call to realizeSprZeroCache() is also made.  Once the Cache for
// the SprZero has been realized, it may be repeatedly accessed with only
// the cost of the method call and the if statement.
void
realizeSprZeroCache(const State& state) const {
    if(isCacheValueRealized(state, indexSprZeroInCache)) return;
    Vec3 sprZero = getSprZero(state);
    Real time = state.getTime();
    sprZero[0] = 0.01 * time;
    sprZero[1] = 0.0;
    sprZero[2] = 0.01 * time;
    updSprZeroInCache(state, sprZero);
    markCacheValueRealized(state, indexSprZeroInCache);
}
//_____________________________________________________________________________
// Clamp a value between zero and a maximum value.
static Real
ClampAboveZero(Real value, Real max) {
    if(value > max) return max;
    if(value < 0.0) return 0.0;
    return value;
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
static Real
Sigma(Real t0, Real tau, Real t) {
    Real x = (t - t0) / tau;
    Real s = 1.0 / (1.0 + std::exp(x));
    return s;
}

//-----------------------------------------------------------------------------
// Data Members
//-----------------------------------------------------------------------------
private:
    ExponentialSpringParameters params;
    Transform contactPlane;
    const MobilizedBody& body;
    Vec3 station;
    Real defaultMus;
    Real defaultMuk;
    Vec3 defaultSprZero;
    ExponentialSpringData defaultData;
    mutable DiscreteVariableIndex indexMus;
    mutable DiscreteVariableIndex indexMuk;
    mutable DiscreteVariableIndex indexSprZero;
    mutable CacheEntryIndex indexSprZeroInCache;
    mutable ZIndex indexZ;
    mutable CacheEntryIndex indexData;

};  // end of class ExponentialSpringForceImpl

} // namespace SimTK


using namespace SimTK;
using std::cout;
using std::endl;

//=============================================================================
// Class ExponentialSpringForce
//=============================================================================
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
// Get the Transform specifying the location and orientation of the Contact
// Plane.
const Transform&
ExponentialSpringForce::
getContactPlane() const {
    return getImpl().getContactPlane();
}
//_____________________________________________________________________________
// Get the point on the body that interacts with the contact plane and at
// which the contact force is applied.
const MobilizedBody&
ExponentialSpringForce::
getBody() const {
    return getImpl().getBody();
}
//_____________________________________________________________________________
// Get the point on the body that interacts with the contact plane and at
// which the contact force is applied.
const Vec3&
ExponentialSpringForce::
getStation() const {
    return getImpl().getStation();
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
Real
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
Real
ExponentialSpringForce::
getMuKinetic(const State& state) const {
    return getImpl().getMuKinetic(state);
}
//_____________________________________________________________________________
// Get the value of the Sliding state.
Real
ExponentialSpringForce::
getSliding(const State& state) const {
    return getImpl().getSliding(state);
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

//-----------------------------------------------------------------------------
// Spring Data Accessor Methods
//-----------------------------------------------------------------------------
//_____________________________________________________________________________
// Get the elastic part of the normal force.
Vec3
ExponentialSpringForce::
getNormalForceElasticPart(const State& state, bool inGround) const {
    Vec3 fyElas(0.);
    fyElas[1] = getImpl().getData(state).fyElas;
    if(inGround) fyElas = getContactPlane().xformFrameVecToBase(fyElas);
    return fyElas;
}
//_____________________________________________________________________________
// Get the damping part of the normal force.
Vec3
ExponentialSpringForce::
getNormalForceDampingPart(const State& state, bool inGround) const {
    Vec3 fyDamp(0.);
    fyDamp[1] = getImpl().getData(state).fyDamp;
    if(inGround) fyDamp = getContactPlane().xformFrameVecToBase(fyDamp);
    return fyDamp;
}
//_____________________________________________________________________________
// Get the magnitude of the normal force.
Vec3
ExponentialSpringForce::
getNormalForce(const State& state, bool inGround) const {
    Vec3 fy(0.);
    fy[1] = getImpl().getData(state).fy;
    if(inGround) fy = getContactPlane().xformFrameVecToBase(fy);
    return fy;
}
//_____________________________________________________________________________
// Get the instantaneous coefficient of friction.
Real
ExponentialSpringForce::
getMu(const State& state) const {
    return getImpl().getData(state).mu;
}
//_____________________________________________________________________________
// Get the friction limit.
Real
ExponentialSpringForce::
getFrictionForceLimit(const State& state) const {
    return getImpl().getData(state).fxzLimit;
}
//_____________________________________________________________________________
// Get the elastic part of the friction force.
Vec3
ExponentialSpringForce::
getFrictionForceElasticPart(const State& state, bool inGround) const {
    Vec3 fricElas = getImpl().getData(state).fricElas;;
    if(inGround) fricElas = getContactPlane().xformFrameVecToBase(fricElas);
    return fricElas;
}
//_____________________________________________________________________________
// Get the elastic part of the friction force.
Vec3
ExponentialSpringForce::
getFrictionForceDampingPart(const State& state, bool inGround) const {
    Vec3 fricDamp = getImpl().getData(state).fricDamp;;
    if(inGround) fricDamp = getContactPlane().xformFrameVecToBase(fricDamp);
    return fricDamp;
}
//_____________________________________________________________________________
// Get the total friction force.
Vec3
ExponentialSpringForce::
getFrictionForce(const State& state, bool inGround) const {
    Vec3 fric = getImpl().getData(state).fric;;
    if(inGround) fric = getContactPlane().xformFrameVecToBase(fric);
    return fric;
}

//_____________________________________________________________________________
// Get the spring force applied to the MobilizedBody.
Vec3
ExponentialSpringForce::
getForce(const State& state, bool inGround) const {
    Vec3 force;
    if(inGround) {
        force = getImpl().getData(state).f_G;
    } else {
        force = getImpl().getData(state).f;
    }
    return force;
}
//_____________________________________________________________________________
// Get the position of the spring station.
Vec3
ExponentialSpringForce::
getStationPosition(const State& state, bool inGround) const {
    Vec3 pos_B = getStation();
    Vec3 pos_G = getBody().findStationLocationInGround(state, pos_B);
    if(inGround) return pos_G;
    Vec3 pos = getContactPlane().shiftBaseStationToFrame(pos_G);
    return pos;
}
//_____________________________________________________________________________
// Get the velocity of the spring station.
Vec3
ExponentialSpringForce::
getStationVelocity(const State& state, bool inGround) const {
    Vec3 pos_B = getStation();
    Vec3 vel_G = getBody().findStationVelocityInGround(state, pos_B);
    if(inGround) return vel_G;
    Vec3 vel = getContactPlane().xformBaseVecToFrame(vel_G);
    return vel;
}
//_____________________________________________________________________________
// Get the position of the spring zero.
Vec3
ExponentialSpringForce::
getSpringZeroPosition(const State& state, bool inGround) const {
    Vec3 p0 = getImpl().getSprZero(state);
    if(inGround) {
        p0 = getContactPlane().shiftFrameStationToBase(p0);
    }
    return p0;
}



//-----------------------------------------------------------------------------
// Implmentation Accesssors
//-----------------------------------------------------------------------------
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


