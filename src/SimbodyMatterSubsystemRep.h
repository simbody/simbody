#ifndef SimTK_SIMBODY_MATTER_SUBSYSTEM_REP_H_
#define SimTK_SIMBODY_MATTER_SUBSYSTEM_REP_H_

/* Portions copyright (c) 2005-7 Stanford University and Michael Sherman.
 * Contributors: Derived from NIH IVM code written by Charles Schwieters.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "SimTKcommon.h"
#include "SimTKcommon/internal/SubsystemGuts.h"

#include "simbody/internal/common.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/MobilizedBody.h"

#include "SimbodyTreeState.h"
#include "RigidBodyNode.h"

using namespace SimTK;

class RigidBodyNode;
class ConstraintNode;
class RBDistanceConstraint;
class RBStation;

namespace SimTK {
class SBModelVars;
class SBInstanceVars;
class SBTimeVars;
class SBPositionVars;
class SBVelocityVars;
class SBDynamicsVars;
class SBAccelerationVars;
class SBModelCache;
class SBInstanceCache;
class SBTimeCache;
class SBPositionCache;
class SBVelocityCache;
class SBDynamicsCache;
class SBAccelerationCache;
}

#include <cassert>
#include <vector>
#include <iostream>

typedef std::vector<const RigidBodyNode*>   RBNodePtrList;
typedef Vector_<SpatialVec>           SpatialVecList;

class IVM;
class LengthConstraints;

/*
 * The SimbodyMatterSubsystemRep class owns the tree of mobilizer-connected rigid bodies, called
 * RigidBodyNodes. The tree is stored by levels, with level 0 being ground, level 1
 * being bodies which are connected to ground (base bodies), level 2 connected to
 * level 1 and so on. Nodes at the same level are stored together in an array,
 * but the order does not reflect the logical tree structure; that is maintained
 * via parent & children pointers kept in the nodes.
 */
class SimbodyMatterSubsystemRep : public SimTK::Subsystem::Guts {
public:
    SimbodyMatterSubsystemRep() 
      : Subsystem::Guts("SimbodyMatterSubsystem", "0.5.5"),
        lConstraints(0)
    { 
        clearTopologyCache();
    }

    SimbodyMatterSubsystemRep(const SimbodyMatterSubsystemRep&);
    SimbodyMatterSubsystemRep& operator=(const SimbodyMatterSubsystemRep&);


        // IMPLEMENTATIONS OF SUBSYSTEM::GUTS VIRTUAL METHODS

    // Destructor is virtual
    ~SimbodyMatterSubsystemRep() {
        invalidateSubsystemTopologyCache();
        clearTopologyCache(); // should do cache before state
        clearTopologyState();
    }

    SimbodyMatterSubsystemRep* cloneImpl() const {
        return new SimbodyMatterSubsystemRep(*this);
    }

    int realizeSubsystemTopologyImpl    (State&) const;
    int realizeSubsystemModelImpl       (State&) const;
    int realizeSubsystemInstanceImpl    (const State&) const;
    int realizeSubsystemTimeImpl        (const State&) const;
    int realizeSubsystemPositionImpl    (const State&) const;
    int realizeSubsystemVelocityImpl    (const State&) const;
    int realizeSubsystemDynamicsImpl    (const State&) const;
    int realizeSubsystemAccelerationImpl(const State&) const;
    int realizeSubsystemReportImpl      (const State&) const;

    int calcDecorativeGeometryAndAppendImpl
       (const State& s, Stage stage, Array<DecorativeGeometry>& geom) const;

    // TODO: these are just unit weights and tolerances. They should be calculated
    // to be something more reasonable.

    int calcQUnitWeightsImpl(const State& s, Vector& weights) const {
        weights.resize(getNQ(s));
        weights = 1; // default says everyone's opinion is just as valid
        return 0;
    }
    int calcUUnitWeightsImpl(const State& s, Vector& weights) const {
        weights.resize(getNU(s));
        weights = 1;
        return 0;
    }
    int calcZUnitWeightsImpl(const State& s, Vector& weights) const {
        weights.resize(getNZ(s));
        weights = 1;
        return 0;
    }
    int calcQErrUnitTolerancesImpl(const State& s, Vector& tolerances) const {
        tolerances.resize(getNQErr(s));
        tolerances = 1;
        return 0;
    }
    int calcUErrUnitTolerancesImpl(const State& s, Vector& tolerances) const {
        tolerances.resize(getNUErr(s));
        tolerances = 1;
        return 0;
    }

        // END OF VIRTUALS.

    // Return the MultibodySystem which owns this MatterSubsystem.
    const MultibodySystem& getMultibodySystem() const {
        return MultibodySystem::downcast(getSystem());
    }


        // CONSTRUCTION STAGE //

    // The MatterSubsystemRep takes over ownership of the child
    // MobilizedBody handle (leaving child as a non-owner reference), and makes it
    // a child (outboard body) of the indicated parent. The new child body's id is
    // returned, and will be greater than the parent's id.

    MobilizedBodyId adoptMobilizedBody(MobilizedBodyId parentId, MobilizedBody& child);
    int getNumMobilizedBodies() const {return (int)mobilizedBodies.size();}

    const MobilizedBody& getMobilizedBody(MobilizedBodyId id) const {
        assert(id < (int)mobilizedBodies.size());
        assert(mobilizedBodies[id]);
        return *mobilizedBodies[id];
    }

    // Note that we do not invalidate the subsystem topology cache yet, even
    // though we're handing out a writable reference to a MobilizedBody here.
    // Otherwise every time someone references Ground, e.g., by calling
    // matterSubsys.Ground() the subsystem would have its topology marked
    // invalid. For this reason, and also because the main program normally
    // retains a writable reference to MobilizedBodies, it is essential that every 
    // non-const method of MobilizedBody and its many descendants mark the
    // subsystem topology invalid when called.
    MobilizedBody& updMobilizedBody(MobilizedBodyId id) {
        assert(id < (int)mobilizedBodies.size());
        assert(mobilizedBodies[id]);
        return *mobilizedBodies[id]; // topology not marked invalid yet
    }

    void createGroundBody();

    const MobilizedBody::Ground& getGround() const {
        return MobilizedBody::Ground::downcast(getMobilizedBody(MobilizedBodyId(0)));
    }
    MobilizedBody::Ground& updGround() {
        return MobilizedBody::Ground::updDowncast(updMobilizedBody(MobilizedBodyId(0)));
    }


    // Constraints are treated similarly to MobilizedBodies here.

    ConstraintId adoptConstraint(Constraint& child);
    int getNumConstraints() const {return (int)constraints.size();}

    const Constraint& getConstraint(ConstraintId id) const {
        assert(id < (int)constraints.size());
        assert(constraints[id]);
        return *constraints[id];
    }
    Constraint& updConstraint(ConstraintId id) {
        assert(id < (int)constraints.size());
        assert(constraints[id]);
        invalidateSubsystemTopologyCache();
        return *constraints[id];
    }

    // MatterSubsystemRep interface. These provide local implementations for
    // virtual methods of MatterSubsystemRep.

    // These counts can be obtained even during construction, where they
    // just return the current counts.
    // includes ground
    int getNBodies()      const {return mobilizedBodies.size();}
    int getNParticles()   const {return 0;} // TODO
    int getNMobilities()  const {return getTotalDOF();}
    int getNConstraints() const {return constraints.size();}
    MobilizedBodyId getParent(MobilizedBodyId) const;
    Array<MobilizedBodyId> getChildren(MobilizedBodyId) const;

    const MassProperties& getDefaultBodyMassProperties    (MobilizedBodyId b) const;
    const Transform&      getDefaultMobilizerFrame        (MobilizedBodyId b) const;
    const Transform&      getDefaultMobilizerFrameOnParent(MobilizedBodyId b) const;

    void findMobilizerQs(const State& s, MobilizedBodyId body, int& qStart, int& nq) const {
        const RigidBodyNode& n = getRigidBodyNode(body);
        qStart = n.getQIndex();
        nq     = n.getNQ(getModelVars(s));
    }

    void findMobilizerUs(const State& s, MobilizedBodyId body, int& uStart, int& nu) const {
        const RigidBodyNode& n = getRigidBodyNode(body);
        uStart = n.getUIndex();
        nu     = n.getDOF();
    }

    // Access to Instance variables. //

    const MassProperties& getBodyMassProperties(const State& s, MobilizedBodyId b) const {
        return getInstanceVars(s).bodyMassProperties[b];
    }
    const Transform& getMobilizerFrame(const State& s, MobilizedBodyId b) const {
        return getInstanceVars(s).outboardMobilizerFrames[b];
    }
    const Transform& getMobilizerFrameOnParent(const State& s, MobilizedBodyId b) const {
        return getInstanceVars(s).inboardMobilizerFrames[b];
    }

    MassProperties& updBodyMassProperties(State& s, MobilizedBodyId b) const {
        return updInstanceVars(s).bodyMassProperties[b];
    }
    Transform& updMobilizerFrame(State& s, MobilizedBodyId b) const {
        return updInstanceVars(s).outboardMobilizerFrames[b];
    }
    Transform& updMobilizerFrameOnParent(State& s, MobilizedBodyId b) const {
        return updInstanceVars(s).inboardMobilizerFrames[b];
    }

    const Vector& getAllParticleMasses(const State& s) const {
        return getInstanceVars(s).particleMasses;
    }
    Vector& updAllParticleMasses(State& s) const {
        return updInstanceVars(s).particleMasses;
    }

    Real getTotalMass(const State& s) const {
        return getInstanceCache(s).totalMass;
    }

    const Transform&  getBodyTransform(const State&, MobilizedBodyId) const;
    const SpatialVec& getBodyVelocity (const State&, MobilizedBodyId) const;

    // velocity dependent
    const SpatialVec& getCoriolisAcceleration     (const State&, MobilizedBodyId) const;
    const SpatialVec& getTotalCoriolisAcceleration(const State&, MobilizedBodyId) const;
    const SpatialVec& getGyroscopicForce          (const State&, MobilizedBodyId) const;
    const SpatialVec& getCentrifugalForces        (const State&, MobilizedBodyId) const;

    // PARTICLES TODO

    const Vector_<Vec3>&  getAllParticleLocations (const State&) const {
        static const Vector_<Vec3> v;
        return v;
    }
    const Vector_<Vec3>&  getAllParticleVelocities(const State&) const {
        static const Vector_<Vec3> v;
        return v;
    }
    // Invalidate Stage::Position.
    Vector_<Vec3>&  updAllParticleLocations (State&) const {
        static Vector_<Vec3> v;
        return v;
    }
    // Invalidate Stage::Velocity.
    Vector_<Vec3>&  updAllParticleVelocities(State&) const {
        static Vector_<Vec3> v;
        return v;
    }
    const Vector_<Vec3>&  getAllParticleAccelerations(const State&) const {
        static const Vector_<Vec3> v;
        return v;
    }

    // TODO: this is unweighted RMS norm
    Real calcQConstraintNorm(const State& s) const {
        const Vector& qerr = getQErr(s);
        return qerr.size() ? std::sqrt(qerr.normSqr()/qerr.size()) : 0.;
    }

    // TODO: this is unweighted, untimescaled RMS norm
    Real calcUConstraintNorm(const State& s) const {
        const Vector& uerr = getUErr(s);
        return uerr.size() ? std::sqrt(uerr.normSqr()/uerr.size()) : 0.;
    }

    // TODO: this is unweighted, untimescaled RMS norm
    Real calcUDotConstraintNorm(const State& s) const {
        const Vector& uderr = getUDotErr(s);
        return uderr.size() ? std::sqrt(uderr.normSqr()/uderr.size()) : 0.;
    }

    bool projectQConstraints(State& s, Vector& y_err, Real tol, Real targetTol) const {
        // TODO
        enforcePositionConstraints(s, tol, targetTol);
        return true;
    }
    bool projectUConstraints(State& s, Vector& y_err, Real tol, Real targetTol) const {
        // TODO
        enforceVelocityConstraints(s, tol, targetTol);
        return true;
    }


    Real calcKineticEnergy(const State&) const;


    /// Calculate the product J*X where J is the partial velocity Jacobian dV/du
    /// and X is a vector of SpatialVec's, one per body. See Eq. 76&77 in
    /// Schwieters' paper, and see 81a & b for a use of this routine to compute
    /// energy gradient in internal coordinates. In that case X=dE/dR, that is
    /// the gradient of the energy w.r.t. atomic positions, summed and shifted
    /// to body origins. There we are pretending dR/dq is the same as dV/du, which
    /// will be true if dq/dt = u, which works for all cases except quaternions.
    /// Schwieters handles that by using Euler angles for orientation coordinates
    /// when doing minimizations. But note that the routine works in terms of u,
    /// not q, so it produces a meaningful result in all cases, just not one that
    /// can be mapped directly back to quaternion coordinates. This is an O(n)
    /// operator which can be called after realizePosition().
    /// It has no effect on the cache.
    void calcInternalGradientFromSpatial(const State&, 
        const SpatialVecList& X, 
        Vector&               JX) const;

    // Given a set of body forces, return the equivalent set of mobilizer torques 
    // IGNORING CONSTRAINTS.
    // Must be in DynamicsStage so that articulated body inertias are available,
    // however, velocities are ignored. This operator has NO effect on the state
    // cache. It makes a single O(N) pass.
    void calcTreeEquivalentMobilityForces(const State&, 
        const Vector_<SpatialVec>& bodyForces,
        Vector&                    mobilityForces) const;

    void calcTreeAccelerations(const State& s,
        const Vector&              mobilityForces,
        const Vector_<SpatialVec>& bodyForces,
        Vector&                    netHingeForces,
        Vector_<SpatialVec>&       A_GB,
        Vector&                    udot) const; 

    void calcMInverseF(const State& s,
        const Vector&        f,
        Vector_<SpatialVec>& A_GB,
        Vector&              udot) const; 

    // Must be in Stage::Position to calculate qdot = Q*u.
    void calcQDot(const State& s,
        const Vector& u,
        Vector&       qdot) const;

    // Must be in Stage::Velocity to calculate qdotdot = Qdot*u + Q*udot.
    void calcQDotDot(const State& s,
        const Vector& udot,
        Vector&       qdotdot) const;

    void setDefaultModelValues       (const SBTopologyCache&, SBModelVars&)        const;
    void setDefaultInstanceValues    (const SBModelVars&,     SBInstanceVars&)     const;
    void setDefaultTimeValues        (const SBModelVars&,     SBTimeVars&)         const;
    void setDefaultPositionValues    (const SBModelVars&,     Vector& q)           const;
    void setDefaultVelocityValues    (const SBModelVars&,     Vector& u)           const;
    void setDefaultDynamicsValues    (const SBModelVars&,     SBDynamicsVars&)     const;
    void setDefaultAccelerationValues(const SBModelVars&,     SBAccelerationVars&) const;


    // A single Constraint may generate multiple of these constraint equations.
    int getNDistanceConstraints() const {return distanceConstraints.size();}

        // CALLABLE AFTER realizeTopology()

    int getTotalDOF()    const {assert(subsystemTopologyHasBeenRealized()); return DOFTotal;}
    int getTotalSqDOF()  const {assert(subsystemTopologyHasBeenRealized()); return SqDOFTotal;}
    int getTotalQAlloc() const {assert(subsystemTopologyHasBeenRealized()); return maxNQTotal;}

    int getNTopologicalPositionConstraintEquations()     const {
        assert(subsystemTopologyHasBeenRealized()); return nextQErrSlot;
    }
    int getNTopologicalVelocityConstraintEquations()     const {
        assert(subsystemTopologyHasBeenRealized()); return nextUErrSlot;
    }
    int getNTopologicalAccelerationConstraintEquations() const {
        assert(subsystemTopologyHasBeenRealized()); return nextMultSlot;
    }

    int getQIndex(MobilizedBodyId) const;
    int getQAlloc(MobilizedBodyId) const;
    int getUIndex(MobilizedBodyId) const;
    int getDOF   (MobilizedBodyId) const;

    // Modeling info.

    void setUseEulerAngles(State& s, bool useAngles) const;
    void setMobilizerIsPrescribed(State& s, MobilizedBodyId, bool prescribe) const;
    void setConstraintIsEnabled(State& s, ConstraintId constraint, bool enable) const;
    bool getUseEulerAngles(const State& s) const;
    bool isMobilizerPrescribed(const State& s, MobilizedBodyId) const;
    bool isConstraintEnabled(const State& s, ConstraintId constraint) const;

        // CALLABLE AFTER realizeModel()

    int  getNQuaternionsInUse(const State&) const;
    bool isUsingQuaternion(const State&, MobilizedBodyId) const;
    int  getQuaternionIndex(const State&, MobilizedBodyId) const; // -1 if none

    // Call after realizeDynamics()
    const SpatialMat& getArticulatedBodyInertia(const State& s, MobilizedBodyId) const;

    // Call after realizeAcceleration()
    const SpatialVec& getBodyAcceleration(const State& s, MobilizedBodyId) const;

    // This is a solver which generates internal velocities from spatial ones.
    void velFromCartesian(const Vector& pos, Vector& vel) {assert(false);/*TODO*/}

    void enforcePositionConstraints(State&, const Real& requiredTol, const Real& desiredTol) const;
    void enforceVelocityConstraints(State&, const Real& requiredTol, const Real& desiredTol) const;

    // Unconstrained (tree) dynamics 

    // articulated body inertias
    void calcArticulatedBodyInertias(const State&) const;

     // articulated body remainder forces
    void calcZ(const State&, 
        const Vector&              mobilityForces,
        const Vector_<SpatialVec>& bodyForces) const;
    void calcTreeAccel(const State&) const;                // accels with forces from last calcZ

    void fixVel0(State&, Vector& vel) const; // TODO -- yuck

    /// Part of constrained dynamics (TODO -- more to move here)
    void calcY(const State&) const;

    /// Pass in internal forces in T; they will be adjusted by this routine.
    void calcConstraintCorrectedInternalForces(const State&, Vector& T); 

    const RigidBodyNode& getRigidBodyNode(int nodeNum) const {
        const RigidBodyNodeIndex& ix = nodeNum2NodeMap[nodeNum];
        return *rbNodeLevels[ix.level][ix.offset];
    }
    //RigidBodyNode& updRigidBodyNode(int nodeNum) {
    //    const RigidBodyNodeIndex& ix = nodeNum2NodeMap[nodeNum];
    //    return *rbNodeLevels[ix.level][ix.offset];
    //}

    // Add a distance constraint equation and assign it a particular slot in
    // the qErr, uErr, and multiplier arrays.
    // Return the assigned distance constraint index for caller's use.
    int addOneDistanceConstraintEquation(
        const RBStation& s1, const RBStation& s2, const Real& d,
        int qerrIndex, int uerrIndex, int multIndex);

    //
    //    STATE ARCHEOLOGY
    //    These routines know where the bodies are buried (no pun intended).
    //

    // The TopologyCache in the State should be a copy of the one
    // we keep locally here. We always use our local copy rather than
    // this one except for checking that the State looks reasonable.
    const SBTopologyCache& getTopologyCache(const State& s) const {
        assert(subsystemTopologyHasBeenRealized() && topologyCacheIndex >= 0);
        return Value<SBTopologyCache>::downcast
            (s.getCacheEntry(getMySubsystemId(),topologyCacheIndex)).get();
    }
    SBTopologyCache& updTopologyCache(const State& s) const { //mutable
        assert(subsystemTopologyHasBeenRealized() && topologyCacheIndex >= 0);
        return Value<SBTopologyCache>::downcast
            (s.updCacheEntry(getMySubsystemId(),topologyCacheIndex)).upd();
    }

    const SBModelCache& getModelCache(const State& s) const {
        return Value<SBModelCache>::downcast
            (s.getCacheEntry(getMySubsystemId(),topologyCache.modelingCacheIndex)).get();
    }
    SBModelCache& updModelCache(const State& s) const { //mutable
        return Value<SBModelCache>::downcast
            (s.updCacheEntry(getMySubsystemId(),topologyCache.modelingCacheIndex)).upd();
    }

    const SBInstanceCache& getInstanceCache(const State& s) const {
        return Value<SBInstanceCache>::downcast
            (s.getCacheEntry(getMySubsystemId(),getModelCache(s).instanceCacheIndex)).get();
    }
    SBInstanceCache& updInstanceCache(const State& s) const { //mutable
        return Value<SBInstanceCache>::downcast
            (s.updCacheEntry(getMySubsystemId(),getModelCache(s).instanceCacheIndex)).upd();
    }

    const SBTimeCache& getTimeCache(const State& s) const {
        return Value<SBTimeCache>::downcast
            (s.getCacheEntry(getMySubsystemId(),getModelCache(s).timeCacheIndex)).get();
    }
    SBTimeCache& updTimeCache(const State& s) const { //mutable
        return Value<SBTimeCache>::downcast
            (s.updCacheEntry(getMySubsystemId(),getModelCache(s).timeCacheIndex)).upd();
    }

    const SBPositionCache& getPositionCache(const State& s) const {
        return Value<SBPositionCache>::downcast
            (s.getCacheEntry(getMySubsystemId(),getModelCache(s).qCacheIndex)).get();
    }
    SBPositionCache& updPositionCache(const State& s) const { //mutable
        return Value<SBPositionCache>::downcast
            (s.updCacheEntry(getMySubsystemId(),getModelCache(s).qCacheIndex)).upd();
    }

    const SBVelocityCache& getVelocityCache(const State& s) const {
        return Value<SBVelocityCache>::downcast
            (s.getCacheEntry(getMySubsystemId(),getModelCache(s).uCacheIndex)).get();
    }
    SBVelocityCache& updVelocityCache(const State& s) const { //mutable
        return Value<SBVelocityCache>::downcast
            (s.updCacheEntry(getMySubsystemId(),getModelCache(s).uCacheIndex)).upd();
    }

    const SBDynamicsCache& getDynamicsCache(const State& s) const {
        return Value<SBDynamicsCache>::downcast
            (s.getCacheEntry(getMySubsystemId(),getModelCache(s).dynamicsCacheIndex)).get();
    }
    SBDynamicsCache& updDynamicsCache(const State& s) const { //mutable
        return Value<SBDynamicsCache>::downcast
            (s.updCacheEntry(getMySubsystemId(),getModelCache(s).dynamicsCacheIndex)).upd();
    }

    const SBAccelerationCache& getAccelerationCache(const State& s) const {
        return Value<SBAccelerationCache>::downcast
            (s.getCacheEntry(getMySubsystemId(),getModelCache(s).accelerationCacheIndex)).get();
    }
    SBAccelerationCache& updAccelerationCache(const State& s) const { //mutable
        return Value<SBAccelerationCache>::downcast
            (s.updCacheEntry(getMySubsystemId(),getModelCache(s).accelerationCacheIndex)).upd();
    }


    const SBModelVars& getModelVars(const State& s) const {
        return Value<SBModelVars>::downcast
            (s.getDiscreteVariable(getMySubsystemId(),topologyCache.modelingVarsIndex)).get();
    }
    SBModelVars& updModelVars(State& s) const {
        return Value<SBModelVars>::downcast
            (s.updDiscreteVariable(getMySubsystemId(),topologyCache.modelingVarsIndex)).upd();
    }

    const SBInstanceVars& getInstanceVars(const State& s) const {
        return Value<SBInstanceVars>::downcast
            (s.getDiscreteVariable(getMySubsystemId(),getModelCache(s).instanceVarsIndex)).get();
    }
    SBInstanceVars& updInstanceVars(State& s) const {
        return Value<SBInstanceVars>::downcast
            (s.updDiscreteVariable(getMySubsystemId(),getModelCache(s).instanceVarsIndex)).upd();
    }

    const SBTimeVars& getTimeVars(const State& s) const {
        return Value<SBTimeVars>::downcast
            (s.getDiscreteVariable(getMySubsystemId(),getModelCache(s).timeVarsIndex)).get();
    }
    SBTimeVars& updTimeVars(State& s) const {
        return Value<SBTimeVars>::downcast
            (s.updDiscreteVariable(getMySubsystemId(),getModelCache(s).timeVarsIndex)).upd();
    }

    const SBPositionVars& getPositionVars(const State& s) const {
        return Value<SBPositionVars>::downcast
            (s.getDiscreteVariable(getMySubsystemId(),getModelCache(s).qVarsIndex)).get();
    }
    SBPositionVars& updPositionVars(State& s) const {
        return Value<SBPositionVars>::downcast
            (s.updDiscreteVariable(getMySubsystemId(),getModelCache(s).qVarsIndex)).upd();
    }

    const SBVelocityVars& getVelocityVars(const State& s) const {
        return Value<SBVelocityVars>::downcast
            (s.getDiscreteVariable(getMySubsystemId(),getModelCache(s).uVarsIndex)).get();
    }
    SBVelocityVars& updVelocityVars(State& s) const {
        return Value<SBVelocityVars>::downcast
            (s.updDiscreteVariable(getMySubsystemId(),getModelCache(s).uVarsIndex)).upd();
    }

    const SBDynamicsVars& getDynamicsVars(const State& s) const {
        return Value<SBDynamicsVars>::downcast
            (s.getDiscreteVariable(getMySubsystemId(),getModelCache(s).dynamicsVarsIndex)).get();
    }
    SBDynamicsVars& updDynamicsVars(State& s) const {
        return Value<SBDynamicsVars>::downcast
            (s.updDiscreteVariable(getMySubsystemId(),getModelCache(s).dynamicsVarsIndex)).upd();
    }

    const SBAccelerationVars& getAccelerationVars(const State& s) const {
        return Value<SBAccelerationVars>::downcast
            (s.getDiscreteVariable(getMySubsystemId(),getModelCache(s).accelerationVarsIndex)).get();
    }
    SBAccelerationVars& updAccelerationVars(State& s) const {
        return Value<SBAccelerationVars>::downcast
            (s.updDiscreteVariable(getMySubsystemId(),getModelCache(s).accelerationVarsIndex)).upd();
    }

    const SimbodyMatterSubsystem& getMySimbodyMatterSubsystemHandle() const {
        return SimbodyMatterSubsystem::downcast(getOwnerSubsystemHandle());
    }
    SimbodyMatterSubsystem& updMySimbodyMatterSubsystemHandle() {
        return SimbodyMatterSubsystem::updDowncast(updOwnerSubsystemHandle());
    }
private:
    ConstraintId addConstraintNode(ConstraintNode*&);

    void calcTreeForwardDynamicsOperator(const State&,
        const Vector&              mobilityForces,
        const Vector_<Vec3>&       particleForces,
        const Vector_<SpatialVec>& bodyForces,
        const Vector*              extraMobilityForces,
        const Vector_<SpatialVec>* extraBodyForces,
        SBAccelerationCache&       ac,
        Vector&                    udot,
        Vector&                    udotErr) const;

    void calcLoopForwardDynamicsOperator(const State&, 
        const Vector&              mobilityForces,
        const Vector_<Vec3>&       particleForces,
        const Vector_<SpatialVec>& bodyForces,
        SBAccelerationCache&       ac,
        Vector&                    udot,
        Vector&                    multipliers,
        Vector&                    udotErr) const;

    // Given a set of forces, calculate accelerations ignoring
    // constraints, and leave the results in the state cache. 
    // Must have already called realizeDynamics().
    // We also allow some extra forces to be supplied, with the intent
    // that these will be used to deal with internal forces generated
    // by constraints; set the pointers to zero if you don't have any
    // extras to pass in.
    void calcTreeForwardDynamics (const State& s,
        const Vector&              mobilityForces,
        const Vector_<Vec3>&       particleForces,
        const Vector_<SpatialVec>& bodyForces,
        const Vector*              extraMobilityForces,
        const Vector_<SpatialVec>* extraBodyForces) const;

    // Given a set of forces, calculate acclerations resulting from
    // those forces and enforcement of acceleration constraints, and update 
    // the state cache with the results.
    void calcLoopForwardDynamics(const State&,
        const Vector&              mobilityForces,
        const Vector_<Vec3>&       particleForces,
        const Vector_<SpatialVec>& bodyForces) const;

    friend std::ostream& operator<<(std::ostream&, const SimbodyMatterSubsystemRep&);
    friend class SimTK::SimbodyMatterSubsystem;

    struct RigidBodyNodeIndex {
        RigidBodyNodeIndex(int l, int o) : level(l), offset(o) { }
        int level, offset;
    };

    SimTK_DOWNCAST(SimbodyMatterSubsystemRep, Subsystem::Guts);

private:
        // TOPOLOGY "STATE VARIABLES"

    void clearTopologyState(); // note that this requires non-const access

    // The handles in this array are the owners of the MobilizedBodies after they
    // are adopted. The MobilizedBodyId (converted to int) is the index of a
    // MobilizedBody in this array.
    std::vector<MobilizedBody*> mobilizedBodies;

    // Constraints are treated similarly.
    std::vector<Constraint*>    constraints;

        // TOPOLOGY CACHE

    // Our realizeTopology method calls this after all bodies & constraints have been added,
    // to construct part of the topology cache below.
    void endConstruction();

    // The data members here are filled in when realizeTopology() is called.
    // The flag which remembers whether we have realized topology is in 
    // the SubsystemRep base class.
    // Note that a cache is treated as mutable, so the methods that manipulate
    // it are const.

    void clearTopologyCache() const {
        SimbodyMatterSubsystemRep& mthis = const_cast<SimbodyMatterSubsystemRep&>(*this);
        mthis.clearTopologyCache(); // call the non-const version
    }

    // This method should clear out all the data members below.
    void clearTopologyCache();

    // TODO: these state indices and counters should be deferred to realizeModel()
    // so we can have Model stage variables which change the number of state
    // slots needed.

    // Initialize to 0 at beginning of construction. These are for doling
    // out Q & U state variables to the nodes.
    int nextUSlot;
    int nextUSqSlot;
    int nextQSlot;

    // These are similarly for doling out slots for *topological* constraint
    // equations for position errors, velocity errors, and multipliers.
    // Note that quaternion normalization constraints are *not* topological
    // so aren't counted here (although we begin assigning them slots at
    // the end of the qErr's used up here).
    int nextQErrSlot;
    int nextUErrSlot;
    int nextMultSlot;

    int DOFTotal;   // summed over all nodes
    int SqDOFTotal; // sum of squares of ndofs per node
    int maxNQTotal; // sum of dofs with room for quaternions

    SBTopologyCache topologyCache;
    int topologyCacheIndex; // topologyCache is copied here in the State

    // This holds pointers to nodes and serves to map (level,offset) to nodeNum.
    Array<RBNodePtrList>      rbNodeLevels;
    // Map nodeNum to (level,offset).
    std::vector<RigidBodyNodeIndex> nodeNum2NodeMap;

    // These are the distant constraint equations generated by Constraints.
    Array<RBDistanceConstraint*> distanceConstraints;
    
    LengthConstraints* lConstraints;
};

std::ostream& operator<<(std::ostream&, const SimbodyMatterSubsystemRep&);

#endif // SimTK_SIMBODY_MATTER_SUBSYSTEM_REP_H_
