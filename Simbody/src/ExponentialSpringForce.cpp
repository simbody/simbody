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

Note that member variables with a "_G" suffix are expressed in the
Ground frame. Member variables without a "_G" suffix are expressed the
contact plane frame. */
struct ExponentialSpringData {
    // Position of the body station in the ground frame.
    Vec3 p_G;
    // Velocity of the body station in the ground frame.
    Vec3 v_G;
    // Position of the body station in the frame of the contact plane.
    Vec3 p;
    // Velocity of the body station in the frame of the contact plane.
    Vec3 v;
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
    Vec3 fricDampMod1;
    // Total friction force in Model 1.
    Vec3 fricMod1;
    // Elastic part of the friction spring force in Model 2.
    Vec3 fricElasMod2;
    // Damping part of the friction spring force in Model 2.
    Vec3 fricDampMod2;
    // Total friction spring force in Model 2.
    Vec3 fricMod2;
    // Elastic friction force after blending.
    Vec3 fricElas;
    // Damping friction force after blending.
    Vec3 fricDamp;
    // Total friction force after blending.
    Vec3 fric;
    // Resultant force (normal + friction) expressed in the frame of the
    // contact frame. This is the force that will be applied to the body
    // after expressing it in the appropriate frame.
    Vec3 f;
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
ExponentialSpringForceImpl(const Transform& floor,
const MobilizedBody& body, const Vec3& station, Real mus, Real muk,
const ExponentialSpringParameters& params) :
ForceSubsystem::Guts("ExponentialSpringForce", "0.0.1"),
contactPlane(floor), body(body), station(station),
defaultMus(mus), defaultMuk(muk), useBlended(true),
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
const Transform& getContactPlaneTransform() const { return contactPlane; }
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
void setSliding(State& state, Real sliding) {
    updZ(state)[indexZ] = clampInPlace(0.0, sliding, 1.0); }
Real getSliding(const State& state) const {
    return getZ(state)[indexZ]; }
Real getSlidingDotInCache(const State& state) const {
    return getZDot(state)[indexZ]; }
void updSlidingDotInCache(const State& state, Real slidingDot) const {
    updZDot(state)[indexZ] = slidingDot; /* Doesn't invalidate the State. */ }

// SLIDING ACTION
SlidingAction getSlidingAction(const State& state) const {
    return Value<SlidingAction>::
        downcast(getDiscreteVariable(state, indexSlidingAction)); }
SlidingAction& updSlidingAction(State& state) const {
    return Value<SlidingAction>::
        updDowncast(updDiscreteVariable(state, indexSlidingAction)); }
SlidingAction getSlidingActionInCache(const State& state) const {
    return Value<SlidingAction>::
        downcast(getDiscreteVarUpdateValue(state, indexSlidingAction)); }
void updSlidingActionInCache(const State& state,
    SlidingAction action) const {
    // Will not invalidate the State.
    Value<SlidingAction>::updDowncast(
        updDiscreteVarUpdateValue(state, indexSlidingAction)) = action; }

// SPRING ZERO
const Vec3& getSprZero(const State& state) const {
    return Value<Vec3>::
        downcast(getDiscreteVariable(state, indexSprZero)); }
Vec3& updSprZero(State& state) const {
    return Value<Vec3>::
        updDowncast(updDiscreteVariable(state,indexSprZero)); }
Vec3 getSprZeroInCache(const State& state) const {
    return Value<Vec3>::downcast(
        getDiscreteVarUpdateValue(state, indexSprZero)); }
void updSprZeroInCache(const State& state, const Vec3& p0) const {
    // Will not invalidate the State.
    Value<Vec3>::updDowncast(
        updDiscreteVarUpdateValue(state, indexSprZero)) = p0; }

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
// Dynamics - compute the forces modeled by this Subsystem.
//
// "params" references the configurable topology-stage parameters that govern
// the behavior of the exponential spring. These can be changed by the user,
// but the System must be realized at the Topology Stage after any such
// change.
//
// "data" references the key data that are calculated and stored as a
// Cache Entry when the System is realized at the Dynamics Stage.
// These data can be retrieved during a simulation by a reporter or handler,
// for example.
// 
// Variables without a suffix are expressed in the frame of the contact plane.
// 
// Variables with the _G suffix are expressed in the ground frame.
//
// Most every calculation happens in this one method, the calculations for
// setting SlidingDot being the notable exception. The conditions that must
// be met for transitioning to Sliding = 0 (fixed in place) include the
// acceleration of the body station. Therefore, SlidingDot are computed
// in realizeSubsystemAccelerationImpl().
int
realizeSubsystemDynamicsImpl(const State& state) const override {
    // Get current accumulated forces
    const MultibodySystem& system = MultibodySystem::downcast(getSystem());
    const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();
    Vector_<SpatialVec>& forces_G =
        system.updRigidBodyForces(state, Stage::Dynamics);

    // Perform the kinematic and force calculations.
    calcStationKinematics(state);
    calcNormalForce(state);
    calcFrictionForceBlended(state);
    //calcFrictionForceSpringOnly(state);
 
    // The kinematic and force calculations were just stored in the data cache.
    ExponentialSpringData& data = updData(state);

    // Set total spring force expressed in the frame of the Contact Plane.
    data.f = data.fric;     // The x and y components are friction.
    data.f[2] = data.fz;    // The z component is the normal force.

    // Transform the total force to the Ground frame
    data.f_G = contactPlane.xformFrameVecToBase(data.f);

    // Apply the force
    body.applyForceToBodyPoint(state, station, data.f_G, forces_G);

    return 0;
}
//_____________________________________________________________________________
// Calculate the spring station kinematics.
// The normal is defined by the z axis of the contact plane.
// The friction plane is defined by the x and y axes of the contact plane.
// Key quantities are saved in the data cache.
void
calcStationKinematics(const State& state) const {
    // Retrieve a writable reference to the data cache entry.
    ExponentialSpringData& data = updData(state);
    // Get position and velocity of the spring station in Ground
    data.p_G = body.findStationLocationInGround(state, station);
    data.v_G = body.findStationVelocityInGround(state, station);
    // Transform the position and velocity into the contact frame.
    data.p = contactPlane.shiftBaseStationToFrame(data.p_G);
    data.v = contactPlane.xformBaseVecToFrame(data.v_G);
    // Resolve into normal (y) and tangential parts (xz plane)
    // Normal (perpendicular to contact plane)
    data.pz = data.p[2];
    data.vz = data.v[2];
    // Tangent (tangent to contact plane)
    data.pxy = data.p;    data.pxy[2] = 0.0;
    data.vxy = data.v;    data.vxy[2] = 0.0;
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

    // A different take on Sliding... (experimental)
    // Sliding should really be a function of the velocity of the
    // MobilizedBody relative to the contact plane.
    // It doesn't really make sense for Sliding to have rise and decay
    // rates that are independent of the body dynamics.
    //Real vSettle = params.getSettleVelocity();
    //Real tau = params.getSlidingTimeConstant();
    //Real speed = data.vxy.norm();
    //sliding = sigma(0.01, -0.001, speed);
    //if(sliding < 0.0) sliding = 0.0;
    //else if(sliding > 1.0) sliding = 1.0;

    // Compute max friction force based on the instantaneous mu.
    Real mus = getMuStatic(state);
    Real muk = getMuKinetic(state);
    data.mu = mus - sliding * (mus - muk);
    data.fxyLimit = data.mu * data.fz;

    // Friction limit is too small. Set all forces to 0.0.
    if(data.fxyLimit < SignificantReal) {
        Vec3 zero(0.0);
        data.fricMod1 = data.fricDampMod1 = zero;
        data.fricMod2 = data.fricDampMod2 = data.fricElasMod2 = zero;
        data.fric = data.fricDamp = data.fricElas = zero;
        p0 = data.pxy;
        data.limitReached = true;

    // Friction limit is large enough for meaningful calculations.
    } else {
        // Model 1: Pure Damping (when Sliding = 1.0)
        // Friction is the result purely of damping (no elastic term).
        // If damping force is greater than data.fxyLimit, the damping force
        // is capped at data.fxyLimit.
        data.fricDampMod1 = data.fricDampMod2 = -kvFric * data.vxy;
        if(data.fricDampMod1.norm() > data.fxyLimit) {
            data.fricDampMod1 = data.fxyLimit * data.fricDampMod1.normalize();
            data.limitReached = true;
        }
        // Model 2: Damped Linear Spring (when Sliding = 0.0)
        // The elastic component prevents drift while maintaining reasonable
        // integrator step sizes, at least when compared to just increasing the
        // damping coefficient.
        data.fricElasMod2 = -kpFric * (data.pxy - p0);
        data.fricMod2 = data.fricElasMod2 + data.fricDampMod2;
        Real fxyMod2 = data.fricMod2.norm();
        if(fxyMod2 > data.fxyLimit) {
            Real scale = data.fxyLimit / fxyMod2;
            data.fricElasMod2 *= scale;
            data.fricDampMod2 *= scale;
            data.fricMod2 = data.fricElasMod2 + data.fricDampMod2;
            data.limitReached = true;
        }
        // Blend Model 1 and Model 2 according to the Sliding state
        // As Sliding --> 1.0, Model 1 dominates
        // As Sliding --> 0.0, Model 2 dominates
        data.fricElas = data.fricElasMod2 * (1.0 - sliding);
        data.fricDamp = data.fricDampMod2 +
            (data.fricDampMod1 - data.fricDampMod2) * sliding;
        data.fric = data.fricElas + data.fricDamp;
        p0 = data.pxy + data.fricElas / kpFric;  p0[2] = 0.0;
    }

    // Update the spring zero
    updSprZeroInCache(state, p0);
    markDiscreteVarUpdateValueRealized(state, indexSprZero);
}
//_____________________________________________________________________________
// Calculate the friction force using the Spring Only Option.
// Friction occurs in a plane that is parallel to the plane defined by the x
// and y axes of the contact plane.
// Key quantities are saved in the data cache.
void
calcFrictionForceSpringOnly(const State& state) const {
    // Retrieve a writable reference to the data cache entry.
    ExponentialSpringData& data = updData(state);
    // Get the friction parameters
    Real kpFric = params.getFrictionElasticity();
    Real kvFric = params.getFrictionViscosity();
    // Get the sliding state, which is bounded by 0.0 and 1.0.
    Real sliding = getZ(state)[indexZ];
    if(sliding < 0.0) sliding = 0.0;
    else if(sliding > 1.0) sliding = 1.0;
    // Compute the maximum allowed frictional force based on the current
    // coefficient of friction.
    Real mus = getMuStatic(state);
    Real muk = getMuKinetic(state);
    data.mu = mus - sliding * (mus - muk);
    data.fxyLimit = data.mu * data.fz;
    // Zero out stuff used for the Blended Option
    data.fricDampMod1 = data.fricDampMod2 = 0.0;
    data.fricElasMod2 = 0.0;
    data.fricMod1 = data.fricMod2 = 0.0;
    // Some initializations
    Vec3 p0 = getSprZero(state);
    data.limitReached = false;
    // fxyLimit is very small, so set the friction force to zero.
    if(data.fxyLimit < SignificantReal) {
        data.fricDamp = 0.0;
        data.fricElas = 0.0;
        data.fric = 0.0;
        p0 = data.pxy; p0[2] = 0.0;
        data.limitReached = true;
    // fxyLimit is large enough to do some calculations
    } else {
        // Elastic part
        data.fricElas = -kpFric * (data.pxy - p0); /*
        if(data.fricElas.norm() > data.fxyLimit) {
            data.fricElas = data.fxyLimit * data.fricElas.normalize();
            p0 = data.pxy + data.fricElas / kpFric;
            p0[2] = 0.0;
            data.limitReached = true;
        }*/
        // Damping part
        data.fricDamp = -kvFric * data.vxy; /*
        if(data.fricDamp.norm() > data.fxyLimit) {
            data.fricDamp = data.fxyLimit * data.fricDamp.normalize();
        }*/
        // Total
        data.fric = data.fricElas + data.fricDamp;
        Real fxy = data.fric.norm();
        if(fxy > data.fxyLimit) {
            Real scale = data.fxyLimit / fxy;
            data.fricElas *= scale;
            data.fricDamp *= scale;
            data.fric = data.fricElas + data.fricDamp;
            data.limitReached = true;
        }
    }
    // Update the spring zero based on the final elastic force.
    p0 = data.pxy + data.fricElas / kpFric;
    p0[2] = 0.0;
    updSprZeroInCache(state, p0);
    markCacheValueRealized(state, indexSprZeroInCache);
}

//_____________________________________________________________________________
// Acceleration - compute and update the derivatives of continuous states.
// The only such state in ExponentialSpringForce is the Sliding state.
int
realizeSubsystemAccelerationImpl(const State& state) const override {
    // Parameters
    Real kTau = 1.0 / params.getSlidingTimeConstant();
    Real vSettle = params.getSettleVelocity();
    Real aSettle = params.getSettleAcceleration();

    // Current Sliding State
    Real sliding = getZ(state)[indexZ];
    SlidingAction action = getSlidingAction(state);

    // Writable reference to the data cache
    // Values are updated during System::realize(Stage::Dynamics) (see above)
    const ExponentialSpringData& data = getData(state);

    // Decision tree for managing SlidingDot. Two things happen:
    // 1) Assign target to 0.0 or 1.0.
    // 2) Assign next action to Check, Rise, or Decay.
    // Note that the reason for the 0.05 and 0.95 thresholds are because it
    // taks a LONG time for Sliding to decay all the way to 0.0 or rise all
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
            if(!data.limitReached && (data.v_G.norm() < vSettle) &&
                (sliding > 0.05)) {
                // Using another tier of the conditional to avoid computing
                // the acceleration if possible.
                // Computing acceleration takes 48 flops.
                // Computing norm takes ?? flops.
                Vec3 a = body.findStationAccelerationInGround(state, station);
                if(a.norm() < aSettle) {
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
    markDiscreteVarUpdateValueRealized(state, indexSlidingAction);
    Real slidingDot = kTau * (target - sliding);
    updSlidingDotInCache(state, slidingDot);
    //std::cout << "t= " << state.getTime() <<
    //    "  s= " << sliding << ", " << slidingDot <<
    //    "  action= " << action << std::endl;

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
    energy += 0.5 * params.getFrictionElasticity() * r.norm() * r.norm();
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
    Vec3 p = contactPlane.shiftBaseStationToFrame(p_G);
    // Project onto the contact plane.
    p[2] = 0.0;
    // Update the spring zero
    updSprZero(state) = p;
}
//_____________________________________________________________________________
// Note - Not using this method, but keeping it around in case I want
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
    sprZero[1] = 0.01 * time; 
    sprZero[2] = 0.0;
    updSprZeroInCache(state, sprZero);
    markCacheValueRealized(state, indexSprZeroInCache);
}
//_____________________________________________________________________________
// sigma() - a function that transitions smoothly from 0.0 to 1.0 or
// from 1.0 to 0.0.
//
// Note that this function does essentially the same thing as SimTK::stepUp()
// and SimTK::stepDown(). I (FC Anderson) added this sigma() function for
// the sake of comparison to past implementations of an exponential spring
// contact model. Currently, sigma() is not being used in class
// ExponentialSpringForce; however, it might be of use in making future improvements to
// the ExponentialSpringForce class.
//
//   f(t) = 1.0 / {1.0 + exp[(t - t0) / tau]}
//   t0 - time about which the transition is centered.  F(t0) = 0.5.
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
sigma(Real t0, Real tau, Real t) {
    Real x = (t - t0) / tau;
    Real s = 1.0 / (1.0 + std::exp(x));
    return s;
}

//-----------------------------------------------------------------------------
// Data Members
//-----------------------------------------------------------------------------
private:
    bool useBlended;
    ExponentialSpringParameters params;
    ExponentialSpringData defaultData;
    Transform contactPlane;
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
getContactPlaneTransform() const {
    return getImpl().getContactPlaneTransform();
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
// Set the value of the Sliding state.
void
ExponentialSpringForce::
setSliding(State& state, Real sliding) {
    updImpl().setSliding(state, sliding);
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
    Vec3 fzElas(0.);
    fzElas[2] = getImpl().getData(state).fzElas;
    if(inGround) fzElas =
        getContactPlaneTransform().xformFrameVecToBase(fzElas);
    return fzElas;
}
//_____________________________________________________________________________
// Get the damping part of the normal force.
Vec3
ExponentialSpringForce::
getNormalForceDampingPart(const State& state, bool inGround) const {
    Vec3 fzDamp(0.);
    fzDamp[2] = getImpl().getData(state).fzDamp;
    if(inGround) fzDamp =
        getContactPlaneTransform().xformFrameVecToBase(fzDamp);
    return fzDamp;
}
//_____________________________________________________________________________
// Get the magnitude of the normal force.
Vec3
ExponentialSpringForce::
getNormalForce(const State& state, bool inGround) const {
    Vec3 fz(0.);
    fz[2] = getImpl().getData(state).fz;
    if(inGround) fz = getContactPlaneTransform().xformFrameVecToBase(fz);
    return fz;
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
    return getImpl().getData(state).fxyLimit;
}
//_____________________________________________________________________________
// Get the elastic part of the friction force.
Vec3
ExponentialSpringForce::
getFrictionForceElasticPart(const State& state, bool inGround) const {
    Vec3 fricElas = getImpl().getData(state).fricElas;;
    if(inGround) fricElas =
        getContactPlaneTransform().xformFrameVecToBase(fricElas);
    return fricElas;
}
//_____________________________________________________________________________
// Get the elastic part of the friction force.
Vec3
ExponentialSpringForce::
getFrictionForceDampingPart(const State& state, bool inGround) const {
    Vec3 fricDamp = getImpl().getData(state).fricDamp;;
    if(inGround) fricDamp =
        getContactPlaneTransform().xformFrameVecToBase(fricDamp);
    return fricDamp;
}
//_____________________________________________________________________________
// Get the total friction force.
Vec3
ExponentialSpringForce::
getFrictionForce(const State& state, bool inGround) const {
    Vec3 fric = getImpl().getData(state).fric;;
    if(inGround) fric = getContactPlaneTransform().xformFrameVecToBase(fric);
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
    Vec3 pos = getContactPlaneTransform().shiftBaseStationToFrame(pos_G);
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
    Vec3 vel = getContactPlaneTransform().xformBaseVecToFrame(vel_G);
    return vel;
}
//_____________________________________________________________________________
// Get the position of the spring zero.
// ? Should I be returning the spring zero that is stored in the Cache?
// YES! Only the updated value in the cache is guaranteed to be consistent
// with a state.
Vec3
ExponentialSpringForce::
getFrictionSpringZeroPosition(const State& state, bool inGround) const {
    Vec3 p0 = getImpl().getSprZeroInCache(state);
    if(inGround) {
        p0 = getContactPlaneTransform().shiftFrameStationToBase(p0);
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


