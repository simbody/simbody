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
#include "simbody/internal/MobilizedBody_Ground.h"
#include "simbody/internal/ForceSubsystem.h"
#include "simbody/internal/ForceSubsystemGuts.h"
#include "simbody/internal/ExponentialSpringForce.h"


namespace SimTK {

// Begin anonymous namespace for ExponentialSpringData.
// Keeps ExponentialSpringData from being exposed in the SimTK namespace
// and limits the scope of ExponentialSpringData to this file.
namespace {

//=============================================================================
// Struct ExponentialSpringData
//=============================================================================
/* ExponentialSpringData is an internal data structure used by the
implementation class ExponentialSpringForceImpl to store and retrieve
important quantities kept in the State's data cache. End-user access to these
quantities is provided by ExponentialSpringForce's API.

Member variables with an "_G" suffix are expressed in the Ground frame.
Member variables with an "_P" suffix are expressed in the Contact Plane. */
struct ExponentialSpringData {
    // Position of the body station in the ground frame.
    Vec3 p_G;
    // Velocity of the body station in the ground frame.
    Vec3 v_G;
    // Position of the body station in the frame of the contact plane.
    Vec3 p_P;
    // Velocity of the body station in the frame of the contact plane.
    Vec3 v_P;
    // Displacement of the body station normal to the floor expressed in the
    // frame of the contact plane.
    Real pz;
    // Velocity of the body station normal to the contact plane expressed
    // in the frame of the contact plane.
    Real vz;
    // Position of the body station projected onto the contact plane
    // expressed in the frame of the contact plane.
    Vec3 pxy;
    // Velocity of the body station in the contact plane expressed in
    // the frame of the contact plane.
    Vec3 vxy;
    // Elastic force in the normal direction.
    Real fzElas;
    // Damping force in the normal direction.
    Real fzDamp;
    // Total normal force expressed in the frame of the contact plane.
    Real fz;
    // Instantaneous coefficient of friction.
    Real mu;
    // Limit of the friction force.
    Real fxyLimit;
    // Flag indicating if the friction limit was exceeded.
    bool limitReached;
    // Damping part of the friction force in Model 1.
    Vec3 fricDampMod1_P;
    // Total friction force in Model 1.
    Vec3 fricMod1_P;
    // Elastic part of the friction spring force in Model 2.
    Vec3 fricElasMod2_P;
    // Damping part of the friction spring force in Model 2.
    Vec3 fricDampMod2_P;
    // Total friction spring force in Model 2.
    Vec3 fricMod2_P;
    // Elastic friction force after blending.
    Vec3 fricElas_P;
    // Damping friction force after blending.
    Vec3 fricDamp_P;
    // Total friction force after blending.
    Vec3 fric_P;
    // Resultant force (normal + friction) expressed in the frame of the
    // contact frame. This is the force that will be applied to the body
    // after expressing it in the appropriate frame.
    Vec3 f_P;
    // Resultant force (normal + friction) expressed in the Ground frame.
    // This is the force applied to the body.
    Vec3 f_G;
};
} // End anonymous namespace for ExponentialSpringData


 //============================================================================
 // Class ExponentialSpringForceImpl
 //============================================================================
class ExponentialSpringForceImpl : public ForceSubsystem::Guts {
public:

// Flag for managing the Sliding state and SlidingDot.
// An auto update discrete state is used to store the sliding action
// for successive integration steps.
enum SlidingAction {
    Decay = 0,  // Decay all the way to fully fixed (Sliding = 0).
    Rise = 1,   // Rise all the way to fully slipping (Sliding = 1).
    Check = 2   // Check if a transition condition is met.
};

// Constructor
ExponentialSpringForceImpl(const Transform& XContactPlane,
const MobilizedBody& body, const Vec3& station,
Real mus, Real muk, const ExponentialSpringParameters& params) :
ForceSubsystem::Guts("ExponentialSpringForce", "0.0.1"),
X_GP(XContactPlane), body(body), station(station),
defaultMus(mus), defaultMuk(muk),
defaultSprZero(Vec3(0., 0., 0.)),
defaultSlidingAction(SlidingAction::Check), defaultSliding(1.0) {
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
const Transform& getContactPlaneTransform() const { return X_GP; }
const MobilizedBody& getBody() const { return body; }
const Vec3& getStation() const { return station; }

// TOPOLOGY PARAMETERS
const ExponentialSpringParameters& getParameters() const { return params; }
void setParameters(const ExponentialSpringParameters& params) {
    this->params = params;
    invalidateSubsystemTopologyCache();
}

// DATA CACHE
ExponentialSpringData& updData(const State& state) const {
    return Value<ExponentialSpringData>::updDowncast(updCacheEntry(state,
        indexData));
}
const ExponentialSpringData& getData(const State& state) const {
    return Value<ExponentialSpringData>::downcast(getCacheEntry(state,
        indexData));
}

// SLIDING STATE
void setSliding(State& state, Real sliding) {
    updZ(state)[indexZ] = clampInPlace(0.0, sliding, 1.0);
}
Real getSliding(const State& state) const {
    return getZ(state)[indexZ];
}
Real getSlidingDotInCache(const State& state) const {
    return getZDot(state)[indexZ];
}
void updSlidingDotInCache(const State& state, Real slidingDot) const {
    updZDot(state)[indexZ] = slidingDot; /* Doesn't invalidate the State. */
}

// SLIDING ACTION
SlidingAction getSlidingAction(const State& state) const {
    return Value<SlidingAction>::
        downcast(getDiscreteVariable(state, indexSlidingAction));
}
SlidingAction& updSlidingAction(State& state) const {
    return Value<SlidingAction>::
        updDowncast(updDiscreteVariable(state, indexSlidingAction));
}
SlidingAction getSlidingActionInCache(const State& state) const {
    return Value<SlidingAction>::
        downcast(getDiscreteVarUpdateValue(state, indexSlidingAction));
}
void updSlidingActionInCache(const State& state,
    SlidingAction action) const { // Will not invalidate the State.
    Value<SlidingAction>::updDowncast(
        updDiscreteVarUpdateValue(state, indexSlidingAction)) = action;
    markDiscreteVarUpdateValueRealized(state, indexSlidingAction);
}

// SPRING ZERO
const Vec3& getSprZero(const State& state) const {
    return Value<Vec3>::
        downcast(getDiscreteVariable(state, indexSprZero));
}
Vec3& updSprZero(State& state) const {
    return Value<Vec3>::
        updDowncast(updDiscreteVariable(state,indexSprZero));
}
Vec3 getSprZeroInCache(const State& state) const {
    return Value<Vec3>::downcast(
        getDiscreteVarUpdateValue(state, indexSprZero));
}
void updSprZeroInCache(const State& state, const Vec3& p0) const {
    // Will not invalidate the State.
    Value<Vec3>::updDowncast(
        updDiscreteVarUpdateValue(state, indexSprZero)) = p0;
    markDiscreteVarUpdateValueRealized(state, indexSprZero);
}

// STATIC COEFFICENT OF FRICTION
const Real& getMuStatic(const State& state) const {
    return Value<Real>::downcast(getDiscreteVariable(state, indexMus));
}
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
    return Value<Real>::downcast(getDiscreteVariable(state, indexMuk));
}
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
    return new ExponentialSpringForceImpl(X_GP, body, station,
        defaultMus, defaultMuk, params);
}
//_____________________________________________________________________________
// Topology - allocate state variables and the data cache.
int
realizeSubsystemTopologyImpl(State& state) const override {
    // Create a mutableThis
    ExponentialSpringForceImpl* mutableThis =
        const_cast<ExponentialSpringForceImpl*>(this);

    // Coefficients of friction: mus and muk
    mutableThis->indexMus = allocateDiscreteVariable(state,
        Stage::Dynamics, new Value<Real>(defaultMus));
    mutableThis->indexMuk = allocateDiscreteVariable(state,
        Stage::Dynamics, new Value<Real>(defaultMuk));

    // SprZero
    mutableThis->indexSprZero =
        allocateAutoUpdateDiscreteVariable(state, Stage::Dynamics,
            new Value<Vec3>(defaultSprZero), Stage::Dynamics);
    mutableThis->indexSprZeroInCache =
        getDiscreteVarUpdateIndex(state, indexSprZero);

    // SlidingAction
    mutableThis->indexSlidingAction =
        allocateAutoUpdateDiscreteVariable(state, Stage::Acceleration,
            new Value<SlidingAction>(defaultSlidingAction), Stage::Dynamics);
    mutableThis->indexSlidingActionInCache =
        getDiscreteVarUpdateIndex(state, indexSlidingAction);

    // Sliding
    const Vector zInit(1, defaultSliding);
    mutableThis->indexZ = allocateZ(state, zInit);

    // Data
    mutableThis->indexData = allocateCacheEntry(state, Stage::Dynamics,
        new Value<ExponentialSpringData>(defaultData));

    return 0;
}
//_____________________________________________________________________________
// Stage::Position - compute the positions needed by this Subsystem.
//
// "data" is a struct that stores the key quantities that are calculated and
// stored as cache entries. Values are updated in the data cache when the
// System is realized at the following stages: Position, Velocity, Dynamics,
// Acceleration. These data can be retrieved during a simulation by a reporter
// or handler, for example.
//
// Variables with a _P suffix are expressed in the frame of the contact plane.
//
// Variables with a _G suffix are expressed in the ground frame.
int
realizeSubsystemPositionImpl(const State& state) const override {
    // Retrieve a writable reference to the data cache entry.
    ExponentialSpringData& data = updData(state);
    // Get the position of the body station in Ground
    data.p_G = body.findStationLocationInGround(state, station);
    // Transform the position into the contact plane frame.
    data.p_P = ~X_GP * data.p_G;
    // Resolve into normal (z) and tangential parts (xy plane)
    // Normal (perpendicular to contact plane)
    data.pz = data.p_P[2];
    // Tangent (tangent to contact plane)
    data.pxy = data.p_P;    data.pxy[2] = 0.0;

    return 0;
}
//_____________________________________________________________________________
// Stage::Velocity - compute the velocities needed by this Subsystem.
//
// "data" is a struct that stores the key quantities that are calculated and
// stored as cache entries. Values are updated in the data cache when the
// System is realized at the following stages: Position, Velocity, Dynamics,
// Acceleration. These data can be retrieved during a simulation by a reporter
// or handler, for example.
//
// Variables with a _P suffix are expressed in the frame of the contact plane.
//
// Variables with a _G suffix are expressed in the ground frame.
int
realizeSubsystemVelocityImpl(const State& state) const override {
    // Retrieve a writable reference to the data cache entry.
    ExponentialSpringData& data = updData(state);
    // Get the velocity of the spring station in Ground
    data.v_G = body.findStationVelocityInGround(state, station);
    // Transform the velocity into the contact plane frame.
    data.v_P = ~X_GP.R() * data.v_G;
    // Resolve into normal (z) and tangential parts (xy plane)
    // Normal (perpendicular to contact plane)
    data.vz = data.v_P[2];
    // Tangent (tangent to contact plane)
    data.vxy = data.v_P;    data.vxy[2] = 0.0;

    return 0;
}
//_____________________________________________________________________________
// Stage::Dynamics - compute the forces modeled by this Subsystem.
//
// "params" references the configurable topology-stage parameters that govern
// the behavior of the exponential spring. These can be changed by the user,
// but the System must be realized at the Topology Stage after any such
// change.
//
// "data" is a struct that stores the key quantities that are calculated and
// stored as cache entries. Values are updated in the data cache when the
// System is realized at the following stages: Position, Velocity, Dynamics,
// Acceleration. These data can be retrieved during a simulation by a reporter
// or handler, for example.
//
// Variables with a _P suffix are expressed in the frame of the contact plane.
//
// Variables with a _G suffix are expressed in the ground frame.
int
realizeSubsystemDynamicsImpl(const State& state) const override {
    // Get current accumulated forces
    const MultibodySystem& system = MultibodySystem::downcast(getSystem());
    const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();
    Vector_<SpatialVec>& forces_G =
        system.updRigidBodyForces(state, Stage::Dynamics);

    // Perform the kinematic and force calculations.
    calcNormalForce(state);
    calcFrictionForceBlended(state);

    // Kinematic and force calculations were just stored in the data cache.
    ExponentialSpringData& data = updData(state);

    // Set total force expressed in the frame of the Contact Plane.
    data.f_P = data.fric_P;  // The x and y components are friction.
    data.f_P[2] = data.fz;   // The z component is the normal force.

    // Transform the total force to the Ground frame
    data.f_G = X_GP.R() * data.f_P;

    // Apply the force to the body and to Ground.
    // TODO(fcanderson) Add a test to see that Ground registers the
    // reaction of the force applied to the body.
    const MobilizedBody& ground = matter.getGround();
    ground.applyForceToBodyPoint(state, data.p_G, -data.f_G, forces_G);
    body.applyForceToBodyPoint(state, station, data.f_G, forces_G);

    return 0;
}
//_____________________________________________________________________________
// Calculate the normal force.
// The normal is defined by the z axis of the contact plane.
// Key quantities are saved in the data cache.
void
calcNormalForce(const State& state) const {
    // Retrieve a writable reference to the data cache entry.
    ExponentialSpringData& data = updData(state);
    // Get the relevant parameters upfront
    Real d0, d1, d2;
    params.getShapeParameters(d0, d1, d2);
    Real kvNorm = params.getNormalViscosity();
    // Normal Force (perpendicular to contact plane) -------------------------
    // Elastic Part
    data.fzElas = d1 * std::exp(-d2 * (data.pz - d0));
    // Damping Part
    data.fzDamp = -kvNorm * data.vz * data.fzElas;
    // Total
    data.fz = data.fzElas + data.fzDamp;
    // Don't allow the normal force to be negative or too large.
    // The upper limit can be justified as a crude model of material yielding.
    // Note that conservation of energy may fail if the material actually
    // yields.
    // The lower limit just means that the contact plane will not pull the
    // MoblizedBody down.
    // Make sure that any change in fz is accompanied by an adjustment
    // in fzElas and fzDamp so that 'fz = fzElas + fzDamp' remains true.
    if(data.fz < 0.0) {
        data.fz = 0.0;
        data.fzDamp = -data.fzElas;
    }
    if(data.fz > 100000.0) {
        data.fz = 100000.0;
        data.fzElas = data.fz - data.fzDamp;
    }
}
//_____________________________________________________________________________
// Calculate the friction force using the Blended Option.
// The friction plane is defined by the x and y axes of the contact plane.
// Key quantities are saved in the data cache.
void
calcFrictionForceBlended(const State& state) const {
    // Initializations
    ExponentialSpringData& data = updData(state);  // data cache
    Real kpFric = params.getFrictionElasticity();   // elasticity
    Real kvFric = params.getFrictionViscosity();    // viscosity
    Vec3 p0 = getSprZero(state);            // spring zero
    data.limitReached = false;              // true if force limit exceeded

    // Bound Sliding by 0.0 and 1.0.
    Real sliding = getZ(state)[indexZ];
    if(sliding < 0.0) sliding = 0.0;
    else if(sliding > 1.0) sliding = 1.0;

    // Compute max friction force based on the instantaneous mu.
    Real mus = getMuStatic(state);
    Real muk = getMuKinetic(state);
    data.mu = mus - sliding * (mus - muk);
    data.fxyLimit = data.mu * data.fz;

    // Friction limit is too small. Set all forces to 0.0.
    if(data.fxyLimit < SignificantReal) {
        Vec3 zero(0.0);
        data.fricMod1_P = data.fricDampMod1_P = zero;
        data.fricMod2_P = data.fricDampMod2_P = data.fricElasMod2_P = zero;
        data.fric_P = data.fricDamp_P = data.fricElas_P = zero;
        p0 = data.pxy;
        data.limitReached = true;

    // Friction limit is large enough for meaningful calculations.
    } else {
        // Model 1: Pure Damping (when Sliding = 1.0)
        // Friction is the result purely of damping (no elastic term).
        // If damping force is greater than data.fxyLimit, the damping force
        // is capped at data.fxyLimit.
        data.fricDampMod1_P = data.fricDampMod2_P = -kvFric * data.vxy;
        if(data.fricDampMod1_P.normSqr() > square(data.fxyLimit)) {
            data.fricDampMod1_P =
                data.fxyLimit * data.fricDampMod1_P.normalize();
            data.limitReached = true;
        }
        // Model 2: Damped Linear Spring (when Sliding = 0.0)
        // The elastic component prevents drift while maintaining reasonable
        // integrator step sizes, at least when compared to just increasing
        // the damping coefficient.
        data.fricElasMod2_P = -kpFric * (data.pxy - p0);
        data.fricMod2_P = data.fricElasMod2_P + data.fricDampMod2_P;
        Real fxyMod2 = data.fricMod2_P.norm();
        if(fxyMod2 > data.fxyLimit) {
            Real scale = data.fxyLimit / fxyMod2;
            data.fricElasMod2_P *= scale;
            data.fricDampMod2_P *= scale;
            data.fricMod2_P = data.fricElasMod2_P + data.fricDampMod2_P;
            data.limitReached = true;
        }
        // Blend Model 1 and Model 2 according to the Sliding state
        // As Sliding --> 1.0, Model 1 dominates
        // As Sliding --> 0.0, Model 2 dominates
        data.fricElas_P = data.fricElasMod2_P * (1.0 - sliding);
        data.fricDamp_P = data.fricDampMod2_P +
            (data.fricDampMod1_P - data.fricDampMod2_P) * sliding;
        data.fric_P = data.fricElas_P + data.fricDamp_P;
        p0 = data.pxy + data.fricElas_P / kpFric;  p0[2] = 0.0;
    }

    // Update the spring zero
    updSprZeroInCache(state, p0);
}
//_____________________________________________________________________________
// Stage::Acceleration - compute and update the derivatives of continuous,
// acceleration-dependent states.
//
// Two states are managed at this Stage: Sliding amd SlidingAction.
//
// Sliding is a continuous state (a "Z" in Simbody vocabuary), and its time
// derivative (SlidingDot) is set here. Sliding is bound between
// 0.0 (indicating that the body station is fixed in place or static) and
// 1.0 (indicating that the body station is moving or kinetic).
//
// SlidingAction is a discrete state used to trigger when Sliding should
// rise to 1.0 (kinetic - fully sliding) or decay to 0.0 (static - fully fixed
// in place). Some conditions for triggering a rise or decay depend on the
// acceleration of the body station. Thus, the need to manage these states
// at the Acceleration Stage.
int
realizeSubsystemAccelerationImpl(const State& state) const override {
    // Parameters
    Real kTau = 1.0 / params.getSlidingTimeConstant();
    Real vSettle = params.getSettleVelocity();
    Real aSettle = params.getSettleAcceleration();

    // Current Sliding State
    Real sliding = getZ(state)[indexZ];
    SlidingAction action = getSlidingAction(state);

    // Const reference to the data cache
    // Values were updated in realizeSubsystemDynamicsImpl() (see above).
    const ExponentialSpringData& data = getData(state);

    // Decision tree for managing SlidingDot. Two things happen:
    // 1) Assign target to 0.0 or 1.0.
    // 2) Assign next action to Check, Rise, or Decay.
    // Note that the reason for the 0.05 and 0.95 thresholds are because it
    // takes a LONG time for Sliding to decay all the way to 0.0 or rise all
    // the way to 1.0 (i.e., much longer than tau). Checking can resume when
    // Sliding gets reasonably close to its target.
    Real target = 1.0;
    if(action == SlidingAction::Check) {

        // Conditions for Rise (Sliding --> 1.0)
        // 1. limitReached = true, OR
        // 2. fz < SimTK::SignificantReal (not "touching" contact plane)
        if((data.limitReached) && (sliding < 0.95)) {
            action = SlidingAction::Rise;

        // Conditions for Decay (Sliding --> 0.0)
        // The requirement is basically static equilibrium.
        // 1. Friction limit not reached, AND
        // 2. |v| < vSettle  (Note: v must be small in ALL directions), AND
        // 3. |a| < aSettle  (Note: a must be small in ALL directions)
        } else {
            if(!data.limitReached && (data.v_G.normSqr() < vSettle * vSettle)
                && (sliding > 0.05)) {
                // Using another tier of the conditional to avoid computing
                // the acceleration if possible.
                // Acceleration takes 48 flops.
                // norm() takes 15 to 20 flops, but we'll use normSqr() which
                // only takes a 5-flop dot product.
                Vec3 a = body.findStationAccelerationInGround(state, station);
                if(a.normSqr() < aSettle * aSettle) {
                    target = 0.0;
                    action = SlidingAction::Decay;
                }
            }
        }

        // If the action is still StateAction::Check, finish transitioning
        // to the current target.
        if((action == SlidingAction::Check) && (sliding < 0.06))
            target = 0.0;

    // Rise
    } else if(action == SlidingAction::Rise) {
        if(sliding >= 0.95) action = SlidingAction::Check;

    // Decay
    } else if(action == SlidingAction::Decay) {
        target = 0.0;
        if(sliding <= 0.05) action = SlidingAction::Check;
    }

    // Update
    updSlidingActionInCache(state, action);
    Real slidingDot = kTau * (target - sliding);
    updSlidingDotInCache(state, slidingDot);

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
    double energy = data.fzElas / d2;
    // Strain energy in the tangent plane (friction spring)
    // Note that the updated spring zero (the one held in cache) needs to be
    // used, not the one in the state.
    // In the process of realizing to Stage::Dynamics, the spring zero is
    // changed when fxzElas > fxzLimit. This change is not reflected in the
    // state, just in the cache.
    Vec3 p0Cache = getSprZeroInCache(state);
    Vec3 r = data.pxy - p0Cache;
    energy += 0.5 * params.getFrictionElasticity() * r.normSqr();
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
    // Express the position in the contact plane.
    Vec3 p_P = ~X_GP * p_G;
    // Project onto the contact plane.
    p_P[2] = 0.0;
    // Update the spring zero
    updSprZero(state) = p_P;
}
//-----------------------------------------------------------------------------
// Data Members
//-----------------------------------------------------------------------------
private:
    ExponentialSpringParameters params;
    ExponentialSpringData defaultData;
    Transform X_GP;
    const MobilizedBody& body;
    Vec3 station;
    Real defaultMus;
    Real defaultMuk;
    Vec3 defaultSprZero;
    SlidingAction defaultSlidingAction;
    Real defaultSliding;
    DiscreteVariableIndex indexMus;
    DiscreteVariableIndex indexMuk;
    DiscreteVariableIndex indexSprZero;
    CacheEntryIndex indexSprZeroInCache;
    DiscreteVariableIndex indexSlidingAction;
    CacheEntryIndex indexSlidingActionInCache;
    ZIndex indexZ;
    CacheEntryIndex indexData;
}; // end of class ExponentialSpringForceImpl
}  // end of namespace SimTK


using namespace SimTK;
using std::cout;
using std::endl;

//=============================================================================
// Class ExponentialSpringForce
//=============================================================================
//_____________________________________________________________________________
ExponentialSpringForce::
ExponentialSpringForce(MultibodySystem& system,
    const Transform& XContactPlane,
    const MobilizedBody& body, const Vec3& station,
    Real mus, Real muk, ExponentialSpringParameters params)
{
    adoptSubsystemGuts(
        new ExponentialSpringForceImpl(XContactPlane, body, station,
            mus, muk, params));
    system.addForceSubsystem(*this);
}
//_____________________________________________________________________________
const Transform&
ExponentialSpringForce::
getContactPlaneTransform() const {
    return getImpl().getContactPlaneTransform();
}
//_____________________________________________________________________________
const MobilizedBody&
ExponentialSpringForce::
getBody() const {
    return getImpl().getBody();
}
//_____________________________________________________________________________
const Vec3&
ExponentialSpringForce::
getStation() const {
    return getImpl().getStation();
}

//_____________________________________________________________________________
// Note that the underlying implementation (ExponentialSpringForceImpl) owns
// its own ExponentialSpringParameters instance. When this method is called,
// the underlying implementation sets its parameters equal to the parameters
// sent in through the argument list by calling the assignment operator:
//
//        ExponentialSpringForceImple::params = params
void
ExponentialSpringForce::
setParameters(const ExponentialSpringParameters& params) {
    updImpl().setParameters(params);
}
//_____________________________________________________________________________
const ExponentialSpringParameters&
ExponentialSpringForce::
getParameters() const {
    return getImpl().getParameters();
}

//_____________________________________________________________________________
void
ExponentialSpringForce::
setMuStatic(State& state, Real mus) {
    updImpl().setMuStatic(state, mus);
}
//_____________________________________________________________________________
Real
ExponentialSpringForce::
getMuStatic(const State& state) const {
    return getImpl().getMuStatic(state);
}

//_____________________________________________________________________________
void
ExponentialSpringForce::
setMuKinetic(State& state, Real muk) {
    updImpl().setMuKinetic(state, muk);
}
//_____________________________________________________________________________
Real
ExponentialSpringForce::
getMuKinetic(const State& state) const {
    return getImpl().getMuKinetic(state);
}
//_____________________________________________________________________________
void
ExponentialSpringForce::
setSliding(State& state, Real sliding) {
    updImpl().setSliding(state, sliding);
}
//_____________________________________________________________________________
Real
ExponentialSpringForce::
getSliding(const State& state) const {
    return getImpl().getSliding(state);
}
//_____________________________________________________________________________
void
ExponentialSpringForce::
resetSpringZero(State& state) const {
    getImpl().resetSprZero(state);
}

//-----------------------------------------------------------------------------
// Data Cache Accessor Methods
//-----------------------------------------------------------------------------
//_____________________________________________________________________________
Vec3
ExponentialSpringForce::
getNormalForceElasticPart(const State& state, bool inGround) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(state), Stage::Dynamics,
        "ExponentialSpringForce::getNormalForceElasticPart");
    Vec3 fzElas(0.);
    fzElas[2] = getImpl().getData(state).fzElas;
    if(inGround) fzElas =
        getContactPlaneTransform().xformFrameVecToBase(fzElas);
    return fzElas;
}
//_____________________________________________________________________________
Vec3
ExponentialSpringForce::
getNormalForceDampingPart(const State& state, bool inGround) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(state), Stage::Dynamics,
        "ExponentialSpringForce::getNormalForceDampingPart");
    Vec3 fzDamp(0.);
    fzDamp[2] = getImpl().getData(state).fzDamp;
    if(inGround) fzDamp =
        getContactPlaneTransform().xformFrameVecToBase(fzDamp);
    return fzDamp;
}
//_____________________________________________________________________________
Vec3
ExponentialSpringForce::
getNormalForce(const State& state, bool inGround) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(state), Stage::Dynamics,
        "ExponentialSpringForce::getNormalForce");
    Vec3 fz(0.);
    fz[2] = getImpl().getData(state).fz;
    if(inGround) fz = getContactPlaneTransform().xformFrameVecToBase(fz);
    return fz;
}
//_____________________________________________________________________________
Real
ExponentialSpringForce::
getMu(const State& state) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(state), Stage::Dynamics,
        "ExponentialSpringForce::getMu");
    return getImpl().getData(state).mu;
}
//_____________________________________________________________________________
Real
ExponentialSpringForce::
getFrictionForceLimit(const State& state) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(state), Stage::Dynamics,
        "ExponentialSpringForce::getFrictionForceLimit");
    return getImpl().getData(state).fxyLimit;
}
//_____________________________________________________________________________
Vec3
ExponentialSpringForce::
getFrictionForceElasticPart(const State& state, bool inGround) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(state), Stage::Dynamics,
        "ExponentialSpringForce::getFrictionForceElasticPart");
    Vec3 fricElas = getImpl().getData(state).fricElas_P;
    if(inGround) fricElas =
        getContactPlaneTransform().xformFrameVecToBase(fricElas);
    return fricElas;
}
//_____________________________________________________________________________
Vec3
ExponentialSpringForce::
getFrictionForceDampingPart(const State& state, bool inGround) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(state), Stage::Dynamics,
        "ExponentialSpringForce::getFrictionForceDampingPart");
    Vec3 fricDamp = getImpl().getData(state).fricDamp_P;
    if(inGround) fricDamp =
        getContactPlaneTransform().xformFrameVecToBase(fricDamp);
    return fricDamp;
}
//_____________________________________________________________________________
Vec3
ExponentialSpringForce::
getFrictionForce(const State& state, bool inGround) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(state), Stage::Dynamics,
        "ExponentialSpringForce::getFrictionForce");
    Vec3 fric = getImpl().getData(state).fric_P;
    if(inGround) fric = getContactPlaneTransform().xformFrameVecToBase(fric);
    return fric;
}
//_____________________________________________________________________________
Vec3
ExponentialSpringForce::
getForce(const State& state, bool inGround) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(state), Stage::Dynamics,
        "ExponentialSpringForce::getForce");
    Vec3 force;
    if(inGround) {
        force = getImpl().getData(state).f_G;
    } else {
        force = getImpl().getData(state).f_P;
    }
    return force;
}
//_____________________________________________________________________________
Vec3
ExponentialSpringForce::
getStationPosition(const State& state, bool inGround) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(state), Stage::Position,
        "ExponentialSpringForce::getStationPosition");
    const ExponentialSpringData& data = getImpl().getData(state);
    if(inGround) return data.p_G;
    return data.p_P;
}
//_____________________________________________________________________________
Vec3
ExponentialSpringForce::
getStationVelocity(const State& state, bool inGround) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(state), Stage::Velocity,
        "ExponentialSpringForce::getStationVelocity");
    const ExponentialSpringData& data = getImpl().getData(state);
    if(inGround) return data.v_G;
    return data.v_P;
}
//_____________________________________________________________________________
// Only the updated value stored in the data cache for the friction spring
// zero (i.e., the elastic anchor point) is guaranteed to be consistent
// with the state. Therefore, in order for this method to give a correct
// value, the System must be realized through Stage::Dynamics prior to a
// call to this method.
Vec3
ExponentialSpringForce::
getFrictionSpringZeroPosition(const State& state, bool inGround) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(state), Stage::Dynamics,
        "ExponentialSpringForce::getFrictionSpringZeroPosition");
    Vec3 p0 = getImpl().getSprZeroInCache(state);
    if(inGround) {
        p0 = getContactPlaneTransform().shiftFrameStationToBase(p0);
    }
    return p0;
}

//-----------------------------------------------------------------------------
// Implementation Accesssors
//-----------------------------------------------------------------------------
//_____________________________________________________________________________
// Get a reference to the underlying implementation that will allow changes
// to be made to underlying member variables and states.
ExponentialSpringForceImpl&
ExponentialSpringForce::
updImpl() {
    return dynamic_cast<ExponentialSpringForceImpl&>(updRep());
}
//_____________________________________________________________________________
// Get a reference to the underlying implementation that will allow
// access, but no change, to underlying parameters and states.
const ExponentialSpringForceImpl&
ExponentialSpringForce::
getImpl() const {
    return dynamic_cast<const ExponentialSpringForceImpl&>(getRep());
}
