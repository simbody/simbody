#ifndef SimTK_SIMBODY_MATTER_SUBSYSTEM_REP_H_
#define SimTK_SIMBODY_MATTER_SUBSYSTEM_REP_H_

/* Portions copyright (c) 2005-6 Stanford University and Michael Sherman.
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
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "simbody/internal/common.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"

#include "SimbodyTreeState.h"
#include "MatterSubsystemRep.h"
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

typedef std::vector<RigidBodyNode*>   RBNodePtrList;
typedef Vector_<SpatialVec>           SpatialVecList;

class IVM;
class LengthConstraints;

/**
 * The SimbodyMatterSubsystemRep class owns the tree of mobilizer-connected rigid bodies, called
 * RigidBodyNodes. The tree is stored by levels, with level 0 being ground, level 1
 * being bodies which are connected to ground (base bodies), level 2 connected to
 * level 1 and so on. Nodes at the same level are stored together in an array,
 * but the order does not reflect the logical tree structure; that is maintained
 * via parent & children pointers kept in the nodes.
 * 
 * SimbodyMatterSubsystemRep is the owner of the RigidBodyNode objects (which are abstract), pointers to
 * which are stored in the tree.
 */
class SimbodyMatterSubsystemRep : public SimTK::MatterSubsystemRep {
public:
    SimbodyMatterSubsystemRep() 
      : MatterSubsystemRep("SimbodyMatterSubsystemRep", "0.5.3"), 
        nextUSlot(0), nextUSqSlot(0), nextQSlot(0), 
        nextQErrSlot(0), nextUErrSlot(0), nextMultSlot(0),
        DOFTotal(-1), SqDOFTotal(-1), maxNQTotal(-1), 
        built(false), topologyCacheIndex(-1), lConstraints(0) 
    { 
        addGroundNode(); 
    }

    SimbodyMatterSubsystemRep(const SimbodyMatterSubsystemRep&);
    SimbodyMatterSubsystemRep& operator=(const SimbodyMatterSubsystemRep&);
    ~SimbodyMatterSubsystemRep();

    // Create a new node, add it to the tree, and assign it
    // a node number, which is a regular labeling starting with node 0 which is ground.
    int addRigidBodyNode
        (RigidBodyNode&           parent,
         const MassProperties&    m,            // mass properties in body frame
         const Transform&         X_PMb,        // parent's frame for attaching this mobilizer
         const Transform&         X_BM,         // mobilizer frame M in body frame
         const Mobilizer&         mobilizer,
         int&                     nxtU,
         int&                     nxtUSq,
         int&                     nxtQ); 

    // Constrain stations on each of two distinct bodies to remain a
    // particular distance apart at all times. Distance must be
    // significantly greater than 0 so that this can be implemented as a
    // single constraint force acting along the instantaneous line between
    // the stations. Parent and child distinction here is meaningless.
    ConstraintId addConstantDistanceConstraint(const RigidBodyNode& parent, const Vec3& stationInP,
                                               const RigidBodyNode& child,  const Vec3& stationInC,
                                               const Real& distance);

    // Constrain stations on each of two distinct bodies to remain superimposed.
    // This restricts all translation but no rotation so adds three constraint
    // equations. Parent and child distinction here is meaningless.
    ConstraintId addCoincidentStationsConstraint(const RigidBodyNode& parent, const Vec3& stationInP,
                                                 const RigidBodyNode& child,  const Vec3& stationInC);

    // Constrain frames fixed to each of two distinct bodies to remain
    // superimposed. Parent and child here mean nothing! This adds six
    // constraint equations.
    ConstraintId addWeldConstraint(const RigidBodyNode& parent, const Transform& frameInP,
                                   const RigidBodyNode& child,  const Transform& frameInC);

    // Call this after all bodies & constraints have been added.
    void endConstruction(); // will set built==true

    // SubsystemRep interface
    SimbodyMatterSubsystemRep* cloneSubsystemRep() const {
        return new SimbodyMatterSubsystemRep(*this);
    }

    // MatterSubsystemRep interface. These provide local implementation sfor
    // virtual methods of MatterSubsystemRep.

    // These counts can be obtained even during construction, where they
    // just return the current counts.
    // includes ground
    int getNBodies()      const {return nodeNum2NodeMap.size();}
    int getNParticles()   const {return 0;} // TODO
    int getNMobilities()  const {return getTotalDOF();}
    int getNConstraints() const {return constraintNodes.size();}
    BodyId getParent(BodyId) const;
    Array<BodyId> getChildren(BodyId) const;

    const MassProperties& getDefaultBodyMassProperties    (BodyId b) const;
    const Transform&      getDefaultMobilizerFrame        (BodyId b) const;
    const Transform&      getDefaultMobilizerFrameOnParent(BodyId b) const;

    void findMobilizerQs(const State& s, BodyId body, int& qStart, int& nq) const {
        const RigidBodyNode& n = getRigidBodyNode(body);
        qStart = n.getQIndex();
        nq     = n.getNQ(getModelVars(s));
    }

    void findMobilizerUs(const State& s, BodyId body, int& uStart, int& nu) const {
        const RigidBodyNode& n = getRigidBodyNode(body);
        uStart = n.getUIndex();
        nu     = n.getDOF();
    }

    // Access to Instance variables. //

    const MassProperties& getBodyMassProperties(const State& s, BodyId b) const {
        return getInstanceVars(s).bodyMassProperties[b];
    }
    const Transform& getMobilizerFrame(const State& s, BodyId b) const {
        return getInstanceVars(s).outboardMobilizerFrames[b];
    }
    const Transform& getMobilizerFrameOnParent(const State& s, BodyId b) const {
        return getInstanceVars(s).inboardMobilizerFrames[b];
    }

    MassProperties& updBodyMassProperties(State& s, BodyId b) const {
        return updInstanceVars(s).bodyMassProperties[b];
    }
    Transform& updMobilizerFrame(State& s, BodyId b) const {
        return updInstanceVars(s).outboardMobilizerFrames[b];
    }
    Transform& updMobilizerFrameOnParent(State& s, BodyId b) const {
        return updInstanceVars(s).inboardMobilizerFrames[b];
    }

    const Vector& getAllParticleMasses(const State& s) const {
        return getInstanceVars(s).particleMasses;
    }
    Vector& updAllParticleMasses(State& s) const {
        return updInstanceVars(s).particleMasses;
    }

    const Vector& getAllMobilizerCoords(const State& s) const {
        return getQ(s); // TODO: return only the non-particle subset
    }
    Vector& updAllMobilizerCoords(State& s) const {
        return updQ(s); // TODO: return only the non-particle subset
    } 

    const Vector& getAllMobilizerSpeeds(const State& s) const {
        return getU(s); // TODO: return only the non-particle subset
    }

    Vector& updAllMobilizerSpeeds(State& s) const {
        return updU(s); // TODO: return only the non-particle subset
    } 

    // Access to Acceleration variables. //

    const Vector& getAllMobilizerAppliedForces(const State& s) const {
        return getAccelerationVars(s).appliedMobilityForces;
    }
    const Vector_<Vec3>& getAllParticleAppliedForces(const State& s) const {
        return getAccelerationVars(s).appliedParticleForces;
    }
    const Vector_<SpatialVec>& getAllBodyAppliedForces(const State& s) const {
        return getAccelerationVars(s).appliedRigidBodyForces;
    }

    // These update routines invalidate Stage::Acceleration.
    Vector& updAllMobilizerAppliedForces(State& s) const {
        return updAccelerationVars(s).appliedMobilityForces;
    }
    Vector_<Vec3>& updAllParticleAppliedForces(State& s) const{
        return updAccelerationVars(s).appliedParticleForces;
    }
    Vector_<SpatialVec>& updAllBodyAppliedForces(State& s) const {
        return updAccelerationVars(s).appliedRigidBodyForces;
    }

    Real getTotalMass(const State& s) const {
        return getInstanceCache(s).totalMass;
    }

    const Transform&  getBodyTransform(const State&, BodyId) const;
    const SpatialVec& getBodyVelocity (const State&, BodyId) const;

    // velocity dependent
    const SpatialVec& getCoriolisAcceleration     (const State&, BodyId) const;
    const SpatialVec& getTotalCoriolisAcceleration(const State&, BodyId) const;
    const SpatialVec& getGyroscopicForce          (const State&, BodyId) const;
    const SpatialVec& getCentrifugalForces        (const State&, BodyId) const;

    const Transform&  getMobilizerTransform(const State&, BodyId) const;
    const SpatialVec& getMobilizerVelocity (const State&, BodyId) const;

    void setMobilizerTransform  (State&, BodyId, const Transform& X_MbM) const;
    void setMobilizerRotation   (State&, BodyId, const Rotation&  R_MbM) const;
    void setMobilizerTranslation(State&, BodyId, const Vec3&      T_MbM,
                                 bool dontChangeOrientation) const;

    void setMobilizerVelocity       (State&, BodyId, const SpatialVec& V_MbM) const;
    void setMobilizerAngularVelocity(State&, BodyId, const Vec3&       w_MbM) const;
    void setMobilizerLinearVelocity (State&, BodyId, const Vec3&       v_MbM,
                                     bool dontChangeAngularVelocity) const;

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


    void realizeTopology    (State&) const;
    void realizeModel       (State&) const;
    void realizeInstance    (const State&) const;
    void realizeTime        (const State&) const;
    void realizePosition    (const State&) const;
    void realizeVelocity    (const State&) const;
    void realizeDynamics    (const State&) const;
    void realizeAcceleration(const State&) const;

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

    int getTotalDOF()    const {assert(built); return DOFTotal;}
    int getTotalSqDOF()  const {assert(built); return SqDOFTotal;}
    int getTotalQAlloc() const {assert(built); return maxNQTotal;}

    int getNTopologicalPositionConstraintEquations()     const {assert(built); return nextQErrSlot;}
    int getNTopologicalVelocityConstraintEquations()     const {assert(built); return nextUErrSlot;}
    int getNTopologicalAccelerationConstraintEquations() const {assert(built); return nextMultSlot;}

    int getQIndex(BodyId) const;
    int getQAlloc(BodyId) const;
    int getUIndex(BodyId) const;
    int getDOF   (BodyId) const;

    // Modeling info.

    void setUseEulerAngles(State& s, bool useAngles) const;
    void setMobilizerIsPrescribed(State& s, BodyId, bool prescribe) const;
    void setConstraintIsEnabled(State& s, int constraint, bool enable) const;
    bool getUseEulerAngles(const State& s) const;
    bool isMobilizerPrescribed(const State& s, BodyId) const;
    bool isConstraintEnabled(const State& s, int constraint) const;

        // CALLABLE AFTER realizeModel()

    int  getNQuaternionsInUse(const State&) const;
    bool isUsingQuaternion(const State&, BodyId) const;
    int  getQuaternionIndex(const State&, BodyId) const; // -1 if none

    // Call after realizeDynamics()
    const SpatialMat& getArticulatedBodyInertia(const State& s, BodyId) const;

    // Call after realizeAcceleration()
    const SpatialVec& getBodyAcceleration(const State& s, BodyId) const;

    /// This is a solver which generates internal velocities from spatial ones.
    void velFromCartesian(const Vector& pos, Vector& vel) {assert(false);/*TODO*/}

    void enforcePositionConstraints(State&, const Real& requiredTol, const Real& desiredTol) const;
    void enforceVelocityConstraints(State&, const Real& requiredTol, const Real& desiredTol) const;

    /// Unconstrained (tree) dynamics 
    void calcArticulatedBodyInertias(const State&) const;                // articulated body inertias
    void calcZ(const State&, const SpatialVecList& spatialForces) const; // articulated body remainder forces
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
    RigidBodyNode& updRigidBodyNode(int nodeNum) {
        const RigidBodyNodeIndex& ix = nodeNum2NodeMap[nodeNum];
        return *rbNodeLevels[ix.level][ix.offset];
    }

    bool isBuilt() const {return built;}

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
        assert(built && topologyCacheIndex >= 0);
        return Value<SBTopologyCache>::downcast
            (s.getCacheEntry(getMySubsystemIndex(),topologyCacheIndex)).get();
    }
    SBTopologyCache& updTopologyCache(const State& s) const { //mutable
        assert(built && topologyCacheIndex >= 0);
        return Value<SBTopologyCache>::downcast
            (s.updCacheEntry(getMySubsystemIndex(),topologyCacheIndex)).upd();
    }

    const SBModelCache& getModelCache(const State& s) const {
        return Value<SBModelCache>::downcast
            (s.getCacheEntry(getMySubsystemIndex(),topologyCache.modelingCacheIndex)).get();
    }
    SBModelCache& updModelCache(const State& s) const { //mutable
        return Value<SBModelCache>::downcast
            (s.updCacheEntry(getMySubsystemIndex(),topologyCache.modelingCacheIndex)).upd();
    }

    const SBInstanceCache& getInstanceCache(const State& s) const {
        return Value<SBInstanceCache>::downcast
            (s.getCacheEntry(getMySubsystemIndex(),getModelCache(s).instanceCacheIndex)).get();
    }
    SBInstanceCache& updInstanceCache(const State& s) const { //mutable
        return Value<SBInstanceCache>::downcast
            (s.updCacheEntry(getMySubsystemIndex(),getModelCache(s).instanceCacheIndex)).upd();
    }

    const SBTimeCache& getTimeCache(const State& s) const {
        return Value<SBTimeCache>::downcast
            (s.getCacheEntry(getMySubsystemIndex(),getModelCache(s).timeCacheIndex)).get();
    }
    SBTimeCache& updTimeCache(const State& s) const { //mutable
        return Value<SBTimeCache>::downcast
            (s.updCacheEntry(getMySubsystemIndex(),getModelCache(s).timeCacheIndex)).upd();
    }

    const SBPositionCache& getPositionCache(const State& s) const {
        return Value<SBPositionCache>::downcast
            (s.getCacheEntry(getMySubsystemIndex(),getModelCache(s).qCacheIndex)).get();
    }
    SBPositionCache& updPositionCache(const State& s) const { //mutable
        return Value<SBPositionCache>::downcast
            (s.updCacheEntry(getMySubsystemIndex(),getModelCache(s).qCacheIndex)).upd();
    }

    const SBVelocityCache& getVelocityCache(const State& s) const {
        return Value<SBVelocityCache>::downcast
            (s.getCacheEntry(getMySubsystemIndex(),getModelCache(s).uCacheIndex)).get();
    }
    SBVelocityCache& updVelocityCache(const State& s) const { //mutable
        return Value<SBVelocityCache>::downcast
            (s.updCacheEntry(getMySubsystemIndex(),getModelCache(s).uCacheIndex)).upd();
    }

    const SBDynamicsCache& getDynamicsCache(const State& s) const {
        return Value<SBDynamicsCache>::downcast
            (s.getCacheEntry(getMySubsystemIndex(),getModelCache(s).dynamicsCacheIndex)).get();
    }
    SBDynamicsCache& updDynamicsCache(const State& s) const { //mutable
        return Value<SBDynamicsCache>::downcast
            (s.updCacheEntry(getMySubsystemIndex(),getModelCache(s).dynamicsCacheIndex)).upd();
    }

    const SBAccelerationCache& getAccelerationCache(const State& s) const {
        return Value<SBAccelerationCache>::downcast
            (s.getCacheEntry(getMySubsystemIndex(),getModelCache(s).accelerationCacheIndex)).get();
    }
    SBAccelerationCache& updAccelerationCache(const State& s) const { //mutable
        return Value<SBAccelerationCache>::downcast
            (s.updCacheEntry(getMySubsystemIndex(),getModelCache(s).accelerationCacheIndex)).upd();
    }


    const SBModelVars& getModelVars(const State& s) const {
        return Value<SBModelVars>::downcast
            (s.getDiscreteVariable(getMySubsystemIndex(),topologyCache.modelingVarsIndex)).get();
    }
    SBModelVars& updModelVars(State& s) const {
        return Value<SBModelVars>::downcast
            (s.updDiscreteVariable(getMySubsystemIndex(),topologyCache.modelingVarsIndex)).upd();
    }

    const SBInstanceVars& getInstanceVars(const State& s) const {
        return Value<SBInstanceVars>::downcast
            (s.getDiscreteVariable(getMySubsystemIndex(),getModelCache(s).instanceVarsIndex)).get();
    }
    SBInstanceVars& updInstanceVars(State& s) const {
        return Value<SBInstanceVars>::downcast
            (s.updDiscreteVariable(getMySubsystemIndex(),getModelCache(s).instanceVarsIndex)).upd();
    }

    const SBTimeVars& getTimeVars(const State& s) const {
        return Value<SBTimeVars>::downcast
            (s.getDiscreteVariable(getMySubsystemIndex(),getModelCache(s).timeVarsIndex)).get();
    }
    SBTimeVars& updTimeVars(State& s) const {
        return Value<SBTimeVars>::downcast
            (s.updDiscreteVariable(getMySubsystemIndex(),getModelCache(s).timeVarsIndex)).upd();
    }

    const SBPositionVars& getPositionVars(const State& s) const {
        return Value<SBPositionVars>::downcast
            (s.getDiscreteVariable(getMySubsystemIndex(),getModelCache(s).qVarsIndex)).get();
    }
    SBPositionVars& updPositionVars(State& s) const {
        return Value<SBPositionVars>::downcast
            (s.updDiscreteVariable(getMySubsystemIndex(),getModelCache(s).qVarsIndex)).upd();
    }

    const SBVelocityVars& getVelocityVars(const State& s) const {
        return Value<SBVelocityVars>::downcast
            (s.getDiscreteVariable(getMySubsystemIndex(),getModelCache(s).uVarsIndex)).get();
    }
    SBVelocityVars& updVelocityVars(State& s) const {
        return Value<SBVelocityVars>::downcast
            (s.updDiscreteVariable(getMySubsystemIndex(),getModelCache(s).uVarsIndex)).upd();
    }

    const SBDynamicsVars& getDynamicsVars(const State& s) const {
        return Value<SBDynamicsVars>::downcast
            (s.getDiscreteVariable(getMySubsystemIndex(),getModelCache(s).dynamicsVarsIndex)).get();
    }
    SBDynamicsVars& updDynamicsVars(State& s) const {
        return Value<SBDynamicsVars>::downcast
            (s.updDiscreteVariable(getMySubsystemIndex(),getModelCache(s).dynamicsVarsIndex)).upd();
    }

    const SBAccelerationVars& getAccelerationVars(const State& s) const {
        return Value<SBAccelerationVars>::downcast
            (s.getDiscreteVariable(getMySubsystemIndex(),getModelCache(s).accelerationVarsIndex)).get();
    }
    SBAccelerationVars& updAccelerationVars(State& s) const {
        return Value<SBAccelerationVars>::downcast
            (s.updDiscreteVariable(getMySubsystemIndex(),getModelCache(s).accelerationVarsIndex)).upd();
    }

    
    // Access to our portion of State arrays.
    void setQ(State& s, const Vector& q) const {
        assert(q.size() == topologyCache.maxNQs);
        updQ(s) = q;
    }
    void setU(State& s, const Vector& u) const {
        assert(u.size() == topologyCache.nDOFs);
        updU(s) = u;
    }


private:
    void addGroundNode();
    ConstraintId addConstraintNode(ConstraintNode*&);

    // Given a forces in the state, calculate accelerations ignoring
    // constraints, and leave the results in the state. 
    // Must have already called realizeDynamics().
    // We also allow some extra forces to be supplied, with the intent
    // that these will be used to deal with internal forces generated
    // by constraints. 
    void calcTreeForwardDynamics (const State& s,
        const Vector*              extraMobilityForces,
        const Vector_<SpatialVec>* extraBodyForces) const;

    // Given a set of forces in the state, calculate acclerations resulting from
    // those forces and enforcement of acceleration constraints, and update the state.
    void calcLoopForwardDynamics(const State&) const;

    friend std::ostream& operator<<(std::ostream&, const SimbodyMatterSubsystemRep&);
    friend class SimTK::SimbodyMatterSubsystem;

    struct RigidBodyNodeIndex {
        RigidBodyNodeIndex(int l, int o) : level(l), offset(o) { }
        int level, offset;
    };

    SimTK_DOWNCAST(SimbodyMatterSubsystemRep, SubsystemRep);

private:
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

    // set by endConstruction
    int DOFTotal;   // summed over all nodes
    int SqDOFTotal; // sum of squares of ndofs per node
    int maxNQTotal; // sum of dofs with room for quaternions
    bool built;

    // set by realizeTopology
    SBTopologyCache topologyCache;
    int topologyCacheIndex; // topologyCache is copied here in the State

    // This holds pointers to nodes and serves to map (level,offset) to nodeNum.
    Array<RBNodePtrList>      rbNodeLevels;
    // Map nodeNum to (level,offset).
    std::vector<RigidBodyNodeIndex> nodeNum2NodeMap;

    // This holds pointers to the abstract constraint nodes which correspond
    // the the user's idea of constraints in a manner analogous to the
    // linked bodies represented by RigidBodyNodes. Each of these may generate
    // several constraint equations.
    Array<ConstraintNode*>       constraintNodes;
    Array<RBDistanceConstraint*> distanceConstraints;
    
    LengthConstraints* lConstraints;

};

std::ostream& operator<<(std::ostream&, const SimbodyMatterSubsystemRep&);

#endif // SimTK_SIMBODY_MATTER_SUBSYSTEM_REP_H_
