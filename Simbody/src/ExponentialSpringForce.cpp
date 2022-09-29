/*-----------------------------------------------------------------------------
                               Simbody(tm)
-------------------------------------------------------------------------------
 Copyright (c) 2021-22 Authors.
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

#include "simbody/internal/common.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/MobilizedBody_Ground.h"
#include "simbody/internal/GeneralForceSubsystem.h"
#include "simbody/internal/ExponentialSpringForce.h"

#include "ForceImpl.h"


namespace SimTK {

// Begin anonymous namespace for ExponentialSpringData.
// Keeps ExponentialSpringData from being exposed in the SimTK namespace
// and limits the scope of ExponentialSpringData to this file.
namespace {

//=============================================================================
// Struct ExponentialSpringData
//=============================================================================
/* ExponentialSpringData is internal data structure used by the implementation
class ExponentialSpringForceImpl to store and retrieve important quantities
kept in the State's data cache. End-user access to these quantities is
provided by ExponentialSpringForce's API.

The struct has three member structs (Pos, Vel, and Dyn) that each correspond
to a realization stage. Splitting up the data cache in this way allows data to
be accessed at the earliest appropriate stage.

Member variables with an "_G" suffix are expressed in the Ground frame.
Member variables with an "_P" suffix are expressed in the Contact Plane. */
struct ExponentialSpringData {
    struct Pos {
        // Position of the body station in the ground frame.
        Vec3 p_G{NaN};
        // Position of the body station in the frame of the contact plane.
        Vec3 p_P{NaN};
        // Displacement of the body station normal to the floor expressed in
        // the frame of the contact plane.
        Real pz{NaN};
        // Position of the body station projected onto the contact plane
        // expressed in the frame of the contact plane.
        Vec3 pxy{NaN};
    };
    struct Vel {
        // Velocity of the body station in the ground frame.
        Vec3 v_G{NaN};
        // Velocity of the body station in the frame of the contact plane.
        Vec3 v_P{NaN};
        // Velocity of the body station normal to the contact plane expressed
        // in the frame of the contact plane.
        Real vz{NaN};
        // Velocity of the body station in the contact plane expressed in
        // the frame of the contact plane.
        Vec3 vxy{NaN};
    };
    struct Dyn {
        // Elastic force in the normal direction.
        Real fzElas{NaN};
        // Damping force in the normal direction.
        Real fzDamp{NaN};
        // Total normal force expressed in the frame of the contact plane.
        Real fz{NaN};
        // Instantaneous coefficient of friction.
        Real mu{NaN};
        // Limit of the friction force.
        Real fxyLimit{NaN};
        // Flag indicating if the friction limit was exceeded.
        bool limitReached{false};
        // Damping part of the friction force in Model 1.
        Vec3 fricDampMod1_P{NaN};
        // Total friction force in Model 1.
        Vec3 fricMod1_P{NaN};
        // Elastic part of the friction spring force in Model 2.
        Vec3 fricElasMod2_P{NaN};
        // Damping part of the friction spring force in Model 2.
        Vec3 fricDampMod2_P{NaN};
        // Total friction spring force in Model 2.
        Vec3 fricMod2_P{NaN};
        // Elastic friction force after blending.
        Vec3 fricElas_P{NaN};
        // Damping friction force after blending.
        Vec3 fricDamp_P{NaN};
        // Total friction force after blending.
        Vec3 fric_P{NaN};
        // Resultant force (normal + friction) expressed in the frame of the
        // contact frame. This is the force that will be applied to the body
        // after expressing it in the appropriate frame.
        Vec3 f_P{NaN};
        // Resultant force (normal + friction) expressed in the Ground frame.
        // This is the force applied to the body.
        Vec3 f_G{NaN};
    };
};
} // End anonymous namespace for ExponentialSpringData


 //============================================================================
 // Class ExponentialSpringForceImpl
 //============================================================================
class ExponentialSpringForceImpl : public ForceImpl {
public:

    // Constructor
    ExponentialSpringForceImpl(const Transform& XContactPlane,
        const MobilizedBody& body, const Vec3& station,
        const ExponentialSpringParameters& params) :
        X_GP(XContactPlane), body(body), station(station),
        defaultAnchorPoint(Vec3(0., 0., 0.)), defaultSliding(1.0)
    {
        this->params = params;
    }

    //-------------------------------------------------------------------------
    // Accessors
    //-------------------------------------------------------------------------
    // CONSTRUCTOR CHOICES THAT ARE NOT CHANGEABLE
    const Transform& getContactPlaneTransform() const { return X_GP; }
    const MobilizedBody& getBody() const { return body; }
    const Vec3& getStation() const { return station; }

    // TOPOLOGY PARAMETERS
    const ExponentialSpringParameters& getParameters() const {return params;}
    void setParameters(const ExponentialSpringParameters& params) {
        this->params = params;
        invalidateTopologyCache();
    }

    // DATA CACHE
    // Position Stage
    ExponentialSpringData::Pos& updDataPos(const State& state) const {
        return Value<ExponentialSpringData::Pos>::updDowncast(
            updCacheEntry(state, indexDataPos));
    }
    const ExponentialSpringData::Pos& getDataPos(const State& state) const {
        return Value<ExponentialSpringData::Pos>::downcast(
            getCacheEntry(state, indexDataPos));
    }
    // Velocity Stage
    ExponentialSpringData::Vel& updDataVel(const State& state) const {
        return Value<ExponentialSpringData::Vel>::updDowncast(
            updCacheEntry(state, indexDataVel));
    }
    const ExponentialSpringData::Vel& getDataVel(const State& state) const {
        return Value<ExponentialSpringData::Vel>::downcast(
            getCacheEntry(state, indexDataVel));
    }
    // Dynamics Stage
    ExponentialSpringData::Dyn& updDataDyn(const State& state) const {
        return Value<ExponentialSpringData::Dyn>::updDowncast(
            updCacheEntry(state, indexDataDyn));
    }
    const ExponentialSpringData::Dyn& getDataDyn(const State& state) const {
        return Value<ExponentialSpringData::Dyn>::downcast(
            getCacheEntry(state, indexDataDyn));
    }

    // ELASTIC ANCHOR POINT
    const Vec3& getAnchorPoint(const State& state) const {
        return Value<Vec3>::downcast(
            getDiscreteVariable(state, indexAnchorPoint));
    }
    Vec3& updAnchorPoint(State& state) const {
        return Value<Vec3>::updDowncast(
            updDiscreteVariable(state, indexAnchorPoint));
    }
    Vec3 getAnchorPointInCache(const State& state) const {
        return Value<Vec3>::downcast(
            getDiscreteVarUpdateValue(state, indexAnchorPoint));
    }
    void updAnchorPointInCache(const State& state, const Vec3& p0) const {
        // Will not invalidate the State.
        Value<Vec3>::updDowncast(
            updDiscreteVarUpdateValue(state, indexAnchorPoint)) = p0;
        markDiscreteVarUpdateValueRealized(state, indexAnchorPoint);
    }

    // SLIDING
    const Real& getSliding(const State& state) const {
        return Value<Real>::downcast(
            getDiscreteVariable(state, indexSliding));
    }
    Real& updSliding(State& state) const {
        return Value<Real>::updDowncast(
            updDiscreteVariable(state, indexSliding));
    }
    Real getSlidingInCache(const State& state) const {
        return Value<Real>::downcast(
            getDiscreteVarUpdateValue(state, indexSliding));
    }
    void updSlidingInCache(const State& state, const Real& p0Speed) const {
        // Will not invalidate the State.
        Value<Real>::updDowncast(
            updDiscreteVarUpdateValue(state, indexSliding)) = p0Speed;
        markDiscreteVarUpdateValueRealized(state, indexSliding);
    }

    // STATIC COEFFICENT OF FRICTION
    const Real& getMuStatic(const State& state) const {
        return Value<Real>::downcast(
            getDiscreteVariable(state, indexMus));
    }
    void setMuStatic(State& state, Real mus) {
        // Keep mus greater than or equal to 0.0.
        if (mus < 0.0) mus = 0.0;
        Value<Real>::updDowncast(updDiscreteVariable(state, indexMus)) = mus;
        // Make sure muk is less than or equal to mus
        Real muk = getMuKinetic(state);
        if (muk > mus) {
            muk = mus;
            Value<Real>::updDowncast(updDiscreteVariable(state, indexMuk))
                = muk;
        }
    }

    // KINETIC COEFFICENT OF FRICTION
    const Real& getMuKinetic(const State& state) const {
        return Value<Real>::downcast(
            getDiscreteVariable(state, indexMuk));
    }
    void setMuKinetic(State& state, Real muk) const {
        // Keep muk >= to zero.
        if (muk < 0.0) muk = 0.0;
        Value<Real>::updDowncast(updDiscreteVariable(state, indexMuk)) = muk;
        // Make sure mus is greater than or equal to muk
        Real mus = getMuStatic(state);
        if (muk > mus) {
            Value<Real>::updDowncast(updDiscreteVariable(state, indexMus))
                = muk;
        }
    }

    //-------------------------------------------------------------------------
    // ForceImpl Methods (overrides of virtual methods)
    //-------------------------------------------------------------------------
    //_________________________________________________________________________
    // Clone
    ForceImpl*
    clone() const override {
        return new ExponentialSpringForceImpl(X_GP, body, station, params);
    }
    //_________________________________________________________________________
    // Topology - allocate state variables and the data cache.
    void
    realizeTopology(State& state) const override {
        // Create a mutableThis
        ExponentialSpringForceImpl* mutableThis =
            const_cast<ExponentialSpringForceImpl*>(this);

        // Get the GeneralForceSubsystem
        const GeneralForceSubsystem& fsub = getForceSubsystem();

        // Coefficients of friction: mus and muk
        mutableThis->indexMus = fsub.allocateDiscreteVariable(state,
            Stage::Dynamics, new Value<Real>(params.getInitialMuStatic()));
        mutableThis->indexMuk = fsub.allocateDiscreteVariable(state,
            Stage::Dynamics, new Value<Real>(params.getInitialMuKinetic()));

        // AnchorPoint (p0) (aka: friction spring zero)
        mutableThis->indexAnchorPoint =
            fsub.allocateAutoUpdateDiscreteVariable(state, Stage::Dynamics,
                new Value<Vec3>(defaultAnchorPoint), Stage::Dynamics);
        mutableThis->indexAnchorPointInCache =
            fsub.getDiscreteVarUpdateIndex(state, indexAnchorPoint);

        // Sliding (K --> for "Kinetic")
        mutableThis->indexSliding = 
            fsub.allocateAutoUpdateDiscreteVariable(state, Stage::Dynamics,
                new Value<Real>(defaultSliding), Stage::Dynamics);
        mutableThis->indexSlidingInCache =
            fsub.getDiscreteVarUpdateIndex(state, indexSliding);

        // Data
        mutableThis->indexDataPos = fsub.allocateCacheEntry(state,
            Stage::Position, new Value<ExponentialSpringData::Pos>());
        mutableThis->indexDataVel = fsub.allocateCacheEntry(state,
            Stage::Velocity, new Value<ExponentialSpringData::Vel>());
        mutableThis->indexDataDyn = fsub.allocateCacheEntry(state,
            Stage::Dynamics, new Value<ExponentialSpringData::Dyn>());
    }
    //_________________________________________________________________________
    // Stage::Position - compute the positions needed by this Force and store
    // them in the data cache.
    //
    // Variables with a _P suffix are expressed in the contact plane frame.
    //
    // Variables with a _G suffix are expressed in the ground frame.
    void
    realizePosition(const State& state) const override {
        // Retrieve a writable reference to the data cache entry.
        ExponentialSpringData::Pos& dataPos = updDataPos(state);
        // Get the position of the body station in Ground
        dataPos.p_G = body.findStationLocationInGround(state, station);
        // Transform the position into the contact plane frame.
        dataPos.p_P = ~X_GP * dataPos.p_G;
        // Resolve into normal (z) and tangential parts (xy plane)
        // Normal (perpendicular to contact plane)
        dataPos.pz = dataPos.p_P[2];
        // Tangent (tangent to contact plane)
        dataPos.pxy = dataPos.p_P;    dataPos.pxy[2] = 0.0;
    }
    //_________________________________________________________________________
    // Stage::Velocity - compute the velocities needed by this Force and store
    // them in the data cache.
    //
    // Variables with a _P suffix are expressed in the contact plane frame.
    //
    // Variables with a _G suffix are expressed in the ground frame.
    void
    realizeVelocity(const State& state) const override {
        // Retrieve a writable reference to the data cache entry.
        ExponentialSpringData::Vel& dataVel = updDataVel(state);
        // Get the velocity of the spring station in Ground
        dataVel.v_G = body.findStationVelocityInGround(state, station);
        // Transform the velocity into the contact plane frame.
        dataVel.v_P = ~X_GP.R() * dataVel.v_G;
        // Resolve into normal (z) and tangential parts (xy plane)
        // Normal (perpendicular to contact plane)
        dataVel.vz = dataVel.v_P[2];
        // Tangent (tangent to contact plane)
        dataVel.vxy = dataVel.v_P;    dataVel.vxy[2] = 0.0;
    }
    //_________________________________________________________________________
    // Calculate the normal force.
    // The normal is defined by the z axis of the contact plane.
    // Key quantities are saved in the data cache.
    void
    calcNormalForce(const State& state) const {
        // Retrieve references to data cache entries.
        const ExponentialSpringData::Pos& dataPos = getDataPos(state);
        const ExponentialSpringData::Vel& dataVel = getDataVel(state);
        ExponentialSpringData::Dyn& dataDyn = updDataDyn(state);
        // Get the relevant parameters upfront
        Real d0, d1, d2;
        params.getShapeParameters(d0, d1, d2);
        Real kvNorm = params.getNormalViscosity();
        // Elastic Part
        dataDyn.fzElas = d1 * std::exp(-d2 * (dataPos.pz - d0));
        // Damping Part
        dataDyn.fzDamp = -kvNorm * dataVel.vz * dataDyn.fzElas;
        // Total
        dataDyn.fz = dataDyn.fzElas + dataDyn.fzDamp;
        // Don't allow the normal force to be negative or too large.
        // The upper limit can be justified as a crude model of material
        // yielding. Note that conservation of energy may fail if the material
        // actually yields.
        // The lower limit just means that the contact plane will not pull the
        // MoblizedBody down.
        // Make sure that any change in fz is accompanied by an adjustment
        // in fzElas and fzDamp so that 'fz = fzElas + fzDamp' remains true.
        if (dataDyn.fz < 0.0) {
            dataDyn.fz = 0.0;
            dataDyn.fzDamp = -dataDyn.fzElas;
        }
        if (dataDyn.fz > params.getMaxNormalForce()) {
            dataDyn.fz = params.getMaxNormalForce();
            dataDyn.fzElas = dataDyn.fz - dataDyn.fzDamp;
        }
    }
    //_________________________________________________________________________
    // Calculate the friction force.
    // The friction plane is defined by the x and y axes of the contact plane.
    // Key quantities are saved in the data cache.
    void
    calcFrictionForce(const State& state) const {
        // Retrieve references to data cache entries.
        const ExponentialSpringData::Pos& dataPos = getDataPos(state);
        const ExponentialSpringData::Vel& dataVel = getDataVel(state);
        ExponentialSpringData::Dyn& dataDyn = updDataDyn(state);
        // Parameters
        Real kpFric = params.getFrictionElasticity();
        Real kvFric = params.getFrictionViscosity();
        // Other initializations
        Vec3 p0 = getAnchorPoint(state);
        Vec3 p0Last = p0;
        Real mus = getMuStatic(state);
        Real muk = getMuKinetic(state);
        dataDyn.limitReached = false;

        // Get the Sliding state (K). ("K" stands for "kinetic")
        // K is handled as an auto update discrete state. It is set at the end
        // of this method for the next time step based on the velocity of p0.
        // The value of K obtained on the line just below is the value set at
        // the end of the previous successful time step.
        Real K = getSliding(state);

        // Compute max friction force based on the instantaneous mu.
        dataDyn.mu = mus - K * (mus - muk);
        dataDyn.fxyLimit = dataDyn.mu * dataDyn.fz;

        // Friction limit is too small.
        // Set all forces to 0.0.
        // Set p0 = pxy.
        // Set limitReached to true.
        if (dataDyn.fxyLimit < SignificantReal) {
            Vec3 zero(0.0);
            dataDyn.fricMod1_P = dataDyn.fricDampMod1_P = zero;
            dataDyn.fricMod2_P = dataDyn.fricDampMod2_P =
                dataDyn.fricElasMod2_P = zero;
            dataDyn.fric_P = dataDyn.fricDamp_P = dataDyn.fricElas_P = zero;
            p0 = dataPos.pxy;
            dataDyn.limitReached = true;

        // Friction limit is large enough for meaningful calculations.
        } else {
            // Model 1: Pure Damping (when Sliding = 1.0)
            // Friction is the result purely of damping (no elastic term).
            // If damping force is greater than data.fxyLimit, the damping
            // force is capped at data.fxyLimit.
            Real fxyLimitSqr = square(dataDyn.fxyLimit);
            dataDyn.fricDampMod1_P =
                dataDyn.fricDampMod2_P = -kvFric*dataVel.vxy;
            if (dataDyn.fricDampMod1_P.normSqr() > fxyLimitSqr) {
                dataDyn.fricDampMod1_P =
                    dataDyn.fxyLimit * dataDyn.fricDampMod1_P.normalize();
                dataDyn.limitReached = true;
            }
            // Model 2: Damped Linear Spring (when Sliding = 0.0)
            // The elastic component prevents drift while maintaining
            // reasonable integrator step sizes, at least when compared to
            // just increasing the damping coefficient.
            dataDyn.fricElasMod2_P = -kpFric * (dataPos.pxy - p0);
            dataDyn.fricMod2_P = dataDyn.fricElasMod2_P +
                dataDyn.fricDampMod2_P;
            Real fxyMod2Sqr = dataDyn.fricMod2_P.normSqr();
            if (fxyMod2Sqr > fxyLimitSqr) {
                Real scale = dataDyn.fxyLimit / std::sqrt(fxyMod2Sqr);
                dataDyn.fricElasMod2_P *= scale;
                dataDyn.fricDampMod2_P *= scale;
                dataDyn.fricMod2_P =
                    dataDyn.fricElasMod2_P + dataDyn.fricDampMod2_P;
                dataDyn.limitReached = true;
            }
            // Blend Model 1 and Model 2 according to K
            // As K --> 1.0, Model 1 dominates
            // As K --> 0.0, Model 2 dominates
            dataDyn.fricElas_P = dataDyn.fricElasMod2_P * (1.0 - K);
            dataDyn.fricDamp_P = dataDyn.fricDampMod2_P +
                (dataDyn.fricDampMod1_P - dataDyn.fricDampMod2_P) * K;
            dataDyn.fric_P = dataDyn.fricElas_P + dataDyn.fricDamp_P;

            // Ensure p0 is consistent with the elastic component
            p0 = dataPos.pxy + dataDyn.fricElas_P / kpFric;  p0[2] = 0.0;
        }

        // Update the elastic anchor point in the data cache
        updAnchorPointInCache(state, p0);

        // Determine the next sliding state (K) of p0.
        K = 1.0;
        Real tDelta = state.getTime() -
            getDiscreteVarLastUpdateTime(state, indexAnchorPoint);
        if (tDelta > SignificantReal && dataDyn.fxyLimit > SignificantReal) {
            // Compute the average speed of the elastic anchor point.
            Vec3 p0Delta(0.0);
            p0Delta = p0 - p0Last;
            Real speed = p0Delta.norm() / tDelta;
            // Normalize speed by the settle velocity
            Real speedFrac = speed / params.getSettleVelocity();
            // Use a step up function to determine K
            K = stepUp(clamp(0.0, speedFrac, 1.0));
        }
        updSlidingInCache(state, K);
    }
    //_________________________________________________________________________
    // Force - calculate the forces and add them to the GeneralForceSubsystem.
    //
    // Note- In a GeneralizedForceSubsystem, calcForce() is called instead of
    // realizeDynamics().
    void
    calcForce(const State& state, Vector_<SpatialVec>& bodyForces,
        Vector_<Vec3>& particleForces,
        Vector& mobilityForces) const override {

        // Calculate the normal and friction forces.
        // Results are stored in the Dyn data cache.
        calcNormalForce(state);
        calcFrictionForce(state);

        // Totals and a transform need to be done for the force.
        ExponentialSpringData::Dyn& dataDyn = updDataDyn(state);
        dataDyn.f_P = dataDyn.fric_P;  // The x and y components are friction.
        dataDyn.f_P[2] = dataDyn.fz;   // The z component is the normal force.
        dataDyn.f_G = X_GP.R() * dataDyn.f_P;  // Transform to Ground

        // Add in the force to the body forces.
        body.applyForceToBodyPoint(state, station, dataDyn.f_G, bodyForces);

        // Add in the force to the Ground forces.
        // TODO(fcanderson) Add a test to see that Ground registers the
        // reaction of the force applied to the body.
        const MultibodySystem& system = MultibodySystem::downcast(
            getForceSubsystem().getSystem());
        const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();
        const MobilizedBody& ground = matter.getGround();
        // The point at which the force should be applied to Ground was
        // stored in the data cache during the call to realizePosition().
        const ExponentialSpringData::Pos& dataPos = getDataPos(state);
        ground.applyForceToBodyPoint(state, dataPos.p_G, -dataDyn.f_G,
            bodyForces);
    }
    //_________________________________________________________________________
    // Potential Energy - calculate the potential energy stored in the spring.
    // The System should be realized through Stage::Dynamics before a call to
    // this method is made.
    //
    // TODO(fcanderson) Correct potential energy calculation when the normal
    // force is capped at its maximum.
    Real
    calcPotentialEnergy(const State& state) const override {
        // Access the data cache
        const ExponentialSpringData::Pos& dataPos = getDataPos(state);
        const ExponentialSpringData::Dyn& dataDyn = getDataDyn(state);
        // Strain energy in the normal direction (exponential spring)
        Real d0, d1, d2;
        params.getShapeParameters(d0, d1, d2);
        double energy = dataDyn.fzElas / d2;
        // Strain energy in the tangent plane (friction spring)
        // Note that the updated anchor point (the one held in cache) needs
        // to be used, not the one in the state.
        // In the process of realizing to Stage::Dynamics, the anchor point is
        // changed when fxzElas > fxzLimit. This change is not reflected in
        // the state, just in the cache.
        Vec3 p0Cache = getAnchorPointInCache(state);
        Vec3 r = dataPos.pxy - p0Cache;
        energy += 0.5 * params.getFrictionElasticity() * r.normSqr();
        return energy;
    }

    //-------------------------------------------------------------------------
    // Utility and Static Methods
    //-------------------------------------------------------------------------
    //_________________________________________________________________________
    // Project the body spring station onto the contact plane.
    void
    resetAnchorPoint(State& state) const {
        // Realize through to the Position Stage
        const MultibodySystem& system = MultibodySystem::downcast(
            getForceSubsystem().getSystem());
        system.realize(state, Stage::Position);
        // Get position of the spring station in the Ground frame
        Vec3 p_G = body.findStationLocationInGround(state, station);
        // Express the position in the contact plane.
        Vec3 p_P = ~X_GP * p_G;
        // Project onto the contact plane.
        p_P[2] = 0.0;
        // Update the elastic anchor point
        updAnchorPoint(state) = p_P;
    }
private:
    //-------------------------------------------------------------------------
    // data members
    //-------------------------------------------------------------------------
    ExponentialSpringParameters params;
    Transform X_GP;
    const MobilizedBody& body;
    Vec3 station;
    Vec3 defaultAnchorPoint;
    Real defaultSliding;
    DiscreteVariableIndex indexMus;
    DiscreteVariableIndex indexMuk;
    DiscreteVariableIndex indexAnchorPoint;
    CacheEntryIndex indexAnchorPointInCache;
    DiscreteVariableIndex indexSliding;
    CacheEntryIndex indexSlidingInCache;
    CacheEntryIndex indexDataPos;
    CacheEntryIndex indexDataVel;
    CacheEntryIndex indexDataDyn;

    //-------------------------------------------------------------------------
    // pass-through methods to reduce calls to getForceSubsystem().
    //-------------------------------------------------------------------------
    // DiscreteVariable
    const AbstractValue& getDiscreteVariable(
        const State& s, DiscreteVariableIndex i) const {
        return getForceSubsystem().getDiscreteVariable(s, i);
    }
    AbstractValue& updDiscreteVariable(
        State& s, DiscreteVariableIndex i) const {
        return getForceSubsystem().updDiscreteVariable(s, i);
    }

    // DiscreteVarUpdateValue
    const AbstractValue& getDiscreteVarUpdateValue(
        const State& s, DiscreteVariableIndex i) const {
        return getForceSubsystem().getDiscreteVarUpdateValue(s, i);
    }
    Real getDiscreteVarLastUpdateTime(
        const State& s, DiscreteVariableIndex i) const {
        return getForceSubsystem().getDiscreteVarLastUpdateTime(s, i);
    }
    AbstractValue& updDiscreteVarUpdateValue(
        const State& s, DiscreteVariableIndex i) const {
        return getForceSubsystem().updDiscreteVarUpdateValue(s, i);
    }
    void markDiscreteVarUpdateValueRealized(
        const State& s, DiscreteVariableIndex i) const {
        getForceSubsystem().markDiscreteVarUpdateValueRealized(s, i);
    }

    // CacheEntry
    const AbstractValue& getCacheEntry(
        const State& s, CacheEntryIndex i) const {
        return getForceSubsystem().getCacheEntry(s, i);
    }
    AbstractValue& updCacheEntry(const State& s, CacheEntryIndex i) const {
        return getForceSubsystem().updCacheEntry(s, i);
    }

}; // end of class ExponentialSpringForceImpl
}  // end of namespace SimTK


using namespace SimTK;
using std::cout;
using std::endl;


SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(
    ExponentialSpringForce, ExponentialSpringForceImpl, Force);

//=============================================================================
// Class ExponentialSpringForce
//=============================================================================
//_____________________________________________________________________________
ExponentialSpringForce::
ExponentialSpringForce(GeneralForceSubsystem& forces,
    const Transform& XContactPlane,
    const MobilizedBody& body, const Vec3& station,
    ExponentialSpringParameters params) :
    Force(new ExponentialSpringForceImpl(XContactPlane, body, station, params))
{
    updImpl().setForceSubsystem(forces, forces.adoptForce(*this));
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
resetAnchorPoint(State& state) const {
    getImpl().resetAnchorPoint(state);
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
    fzElas[2] = getImpl().getDataDyn(state).fzElas;
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
    fzDamp[2] = getImpl().getDataDyn(state).fzDamp;
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
    fz[2] = getImpl().getDataDyn(state).fz;
    if(inGround) fz = getContactPlaneTransform().xformFrameVecToBase(fz);
    return fz;
}
//_____________________________________________________________________________
Real
ExponentialSpringForce::
getMu(const State& state) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(state), Stage::Dynamics,
        "ExponentialSpringForce::getMu");
    return getImpl().getDataDyn(state).mu;
}
//_____________________________________________________________________________
Real
ExponentialSpringForce::
getFrictionForceLimit(const State& state) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(state), Stage::Dynamics,
        "ExponentialSpringForce::getFrictionForceLimit");
    return getImpl().getDataDyn(state).fxyLimit;
}
//_____________________________________________________________________________
Vec3
ExponentialSpringForce::
getFrictionForceElasticPart(const State& state, bool inGround) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(state), Stage::Dynamics,
        "ExponentialSpringForce::getFrictionForceElasticPart");
    Vec3 fricElas = getImpl().getDataDyn(state).fricElas_P;
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
    Vec3 fricDamp = getImpl().getDataDyn(state).fricDamp_P;
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
    Vec3 fric = getImpl().getDataDyn(state).fric_P;
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
        force = getImpl().getDataDyn(state).f_G;
    } else {
        force = getImpl().getDataDyn(state).f_P;
    }
    return force;
}
//_____________________________________________________________________________
Vec3
ExponentialSpringForce::
getStationPosition(const State& state, bool inGround) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(state), Stage::Position,
        "ExponentialSpringForce::getStationPosition");
    const ExponentialSpringData::Pos& data = getImpl().getDataPos(state);
    if(inGround) return data.p_G;
    return data.p_P;
}
//_____________________________________________________________________________
Vec3
ExponentialSpringForce::
getStationVelocity(const State& state, bool inGround) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(state), Stage::Velocity,
        "ExponentialSpringForce::getStationVelocity");
    const ExponentialSpringData::Vel& dataVel = getImpl().getDataVel(state);
    if(inGround) return dataVel.v_G;
    return dataVel.v_P;
}

//-----------------------------------------------------------------------------
// Auto Update Discrete States
//-----------------------------------------------------------------------------
//_____________________________________________________________________________
// Only the updated value stored in the data cache for the elastic anchor
// point is guaranteed to be consistent with the state. Therefore, in order
// for this method to give a correct value, the System must be realized
// through Stage::Dynamics prior to a call to this method.
Vec3
ExponentialSpringForce::
getAnchorPointPosition(const State& state, bool inGround) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(state), Stage::Dynamics,
        "ExponentialSpringForce::getAnchorPointPosition");
    Vec3 p0 = getImpl().getAnchorPointInCache(state);
    if(inGround) {
        p0 = getContactPlaneTransform().shiftFrameStationToBase(p0);
    }
    return p0;
}
//_____________________________________________________________________________
// Only the updated value stored in the data cache for the sliding state is
// guaranteed to be consistent with the state. Therefore, in order for this
// method to give a correct value, the System must be realized through
// Stage::Dynamics prior to a call to this method.
Real
ExponentialSpringForce::
getSliding(const State& state) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(state), Stage::Dynamics,
        "ExponentialSpringForce::getSliding");
    return getImpl().getSlidingInCache(state);
}
