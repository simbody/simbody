#ifndef SimTK_SIMBODY_MATTER_SUBSYSTEM_REP_H_
#define SimTK_SIMBODY_MATTER_SUBSYSTEM_REP_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-7 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
 * Contributors: Derived from NIH IVM code written by Charles Schwieters      *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "SimTKcommon.h"
#include "SimTKcommon/internal/SubsystemGuts.h"

#include "simbody/internal/common.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/SimbodyMatterSubtree.h"
#include "simbody/internal/MobilizedBody.h"

#include "SimbodyTreeState.h"
#include "RigidBodyNode.h"

#include <set>
#include <map>
#include <vector>
#include <utility> // std::pair
using std::pair;


#include <cassert>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

class RigidBodyNode;
class RBDistanceConstraint;
class RBStation;
class RBDirection;

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

using namespace SimTK;

typedef std::vector<const RigidBodyNode*>   RBNodePtrList;
typedef Vector_<SpatialVec>           SpatialVecList;

class IVM;
class LengthConstraints;

/*
 * A CoupledConstraintSet is a set of Simbody Constraints which must be
 * handled simultaneously for purposes of a particular computation. We maintain a
 * Subtree here which encompasses all the included Constraints. This
 * Subtree has as its terminal bodies the set of all constrained bodies
 * from any of the included Constraints. The ancestor body will be
 * the outmost common ancestor of all those constrained bodies, or 
 * equivalently the outmost common ancestor of all the Constraints'
 * ancestor bodies.
 *
 * Operators and responses here set up an environment in which the
 * individual constraint error functions are evaluated. The environments
 * are: 
 *   - same as global System state (i.e., this is just a response)
 *   - all mobility variables from State, except for one which is perturbed (q,u,udot)
 *   - all mobility variables are zero (u,udot)
 *   - all mobility variables are zero except for one (u,udot)
 */

class CoupledConstraintSet {
public:
    CoupledConstraintSet() : topologyRealized(false) { }
    void addConstraint(ConstraintIndex cid) {
        topologyRealized = false;
        constraints.insert(cid);
    }
    
    void mergeInConstraintSet(const CoupledConstraintSet& src) {
        topologyRealized = false;
        for (int i=0; i < (int)src.getCoupledConstraints().size(); ++i)
            constraints.insert(src.getCoupledConstraints()[i]);
    }

    void realizeTopology(const SimbodyMatterSubsystem& matter) {
        coupledConstraints.resize(constraints.size());
        std::set<ConstraintIndex>::const_iterator i = constraints.begin();
        for (int nxt=0; i != constraints.end(); ++i, ++nxt)
            coupledConstraints[nxt] = *i;
        topologyRealized = true;
    }

    const std::vector<ConstraintIndex>& getCoupledConstraints() const {
        assert(topologyRealized);
        return coupledConstraints;
    }

private:
    // TOPOLOGY STATE
    std::set<ConstraintIndex> constraints;

    // TOPOLOGY CACHE
    bool topologyRealized;

    // Sorted in nondecreasing order of ancestor MobilizedBodyIndex.
    std::vector<ConstraintIndex>    coupledConstraints;
    SimbodyMatterSubtree coupledSubtree; // with the new ancestor
};

    //////////////////////////////////
    // SIMBODY MATTER SUBSYSTEM REP //
    //////////////////////////////////

/*
 * The SimbodyMatterSubsystemRep class owns the tree of MobilizedBodies and their
 * associated computational form inherited from IVM, called RigidBodyNodes. 
 * Here we store references to the RigidBodyNodes in a tree structure organized
 * in "levels" with level 0 being Ground, level 1
 * being bodies which are connected to Ground (base bodies), level 2 connected to
 * level 1 and so on. Nodes at the same level are stored together in an array,
 * but the order does not reflect the logical tree structure; that is maintained
 * via parent & children pointers kept in the nodes.
 *
 * Access to mobilized body information is requested via MobilizedBodyIndex's here,
 * which are small integer indices, rather than via MobilizedBody objects. These
 * integers are used both to index MobilizedBodies and their corresponding
 * RigidBodyNodes. You can get the abstract MobilizedBody corresponding to a
 * MobilizedBodyIndex if you want one, but you should prefer MobilizedBodyIndex's when
 * working at the "Rep" level here.
 */
class SimbodyMatterSubsystemRep : public SimTK::Subsystem::Guts {
public:
    SimbodyMatterSubsystemRep() 
      : Subsystem::Guts("SimbodyMatterSubsystem", "0.7.1"),
        lConstraints(0)
    { 
        clearTopologyCache();
    }

    SimbodyMatterSubsystemRep(const SimbodyMatterSubsystemRep&);
    SimbodyMatterSubsystemRep& operator=(const SimbodyMatterSubsystemRep&);


        // IMPLEMENTATIONS OF SUBSYSTEM::GUTS VIRTUAL METHODS

    // Subsystem::Guts destructor is virtual
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
       (const State& s, Stage stage, std::vector<DecorativeGeometry>& geom) const;

    // TODO: these are just unit weights and tolerances. They should be calculated
    // to be something more reasonable.

    int calcQUnitWeightsImpl(const State& s, Vector& weights) const;
    int calcUUnitWeightsImpl(const State& s, Vector& weights) const;
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

        // END OF SUBSYSTEM::GUTS VIRTUALS.

    // Return the MultibodySystem which owns this MatterSubsystem.
    const MultibodySystem& getMultibodySystem() const {
        return MultibodySystem::downcast(getSystem());
    }


        // CONSTRUCTION STAGE //

    // The SimbodyMatterSubsystemRep takes over ownership of the child
    // MobilizedBody handle (leaving child as a non-owner reference), and makes it
    // a child (outboard body) of the indicated parent. The new child body's id is
    // returned, and will be greater than the parent's id.

    MobilizedBodyIndex adoptMobilizedBody(MobilizedBodyIndex parentIndex, MobilizedBody& child);
    int getNumMobilizedBodies() const {return (int)mobilizedBodies.size();}

    const MobilizedBody& getMobilizedBody(MobilizedBodyIndex id) const {
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
    MobilizedBody& updMobilizedBody(MobilizedBodyIndex id) {
        assert(id < (int)mobilizedBodies.size());
        assert(mobilizedBodies[id]);
        return *mobilizedBodies[id]; // topology not marked invalid yet
    }

    void createGroundBody();

    const MobilizedBody::Ground& getGround() const {
        return MobilizedBody::Ground::downcast(getMobilizedBody(MobilizedBodyIndex(0)));
    }
    MobilizedBody::Ground& updGround() {
        return MobilizedBody::Ground::updDowncast(updMobilizedBody(MobilizedBodyIndex(0)));
    }


    // Constraints are treated similarly to MobilizedBodies here.

    ConstraintIndex adoptConstraint(Constraint& child);
    int getNumConstraints() const {return (int)constraints.size();}

    const Constraint& getConstraint(ConstraintIndex id) const {
        assert(id < (int)constraints.size());
        assert(constraints[id]);
        return *constraints[id];
    }
    Constraint& updConstraint(ConstraintIndex id) {
        assert(id < (int)constraints.size());
        assert(constraints[id]);
        return *constraints[id];
    }

    // Topology stage cache entry
    AncestorConstrainedBodyPoolIndex allocateNextAncestorConstrainedBodyPoolSlot() const {
        const AncestorConstrainedBodyPoolIndex nxt = nextAncestorConstrainedBodyPoolSlot;
        // Make this mutable briefly.
        ++ const_cast<SimbodyMatterSubsystemRep*>(this)->nextAncestorConstrainedBodyPoolSlot;
        return nxt;
    }


    // These counts can be obtained even during construction, where they
    // just return the current counts.
    // NBodies includes ground.
    int getNBodies()      const {return mobilizedBodies.size();}
    int getNParticles()   const {return 0;} // TODO
    int getNMobilities()  const {return getTotalDOF();}
    int getNConstraints() const {return constraints.size();}
    MobilizedBodyIndex getParent(MobilizedBodyIndex) const;
    std::vector<MobilizedBodyIndex> getChildren(MobilizedBodyIndex) const;

    const MassProperties& getDefaultBodyMassProperties    (MobilizedBodyIndex b) const;
    const Transform&      getDefaultMobilizerFrame        (MobilizedBodyIndex b) const;
    const Transform&      getDefaultMobilizerFrameOnParent(MobilizedBodyIndex b) const;

    void findMobilizerQs(const State& s, MobilizedBodyIndex body, QIndex& qStart, int& nq) const {
        const RigidBodyNode& n = getRigidBodyNode(body);
        qStart = n.getQIndex();
        nq     = n.getNQInUse(getModelVars(s));
    }

    void findMobilizerUs(const State& s, MobilizedBodyIndex body, UIndex& uStart, int& nu) const {
        const RigidBodyNode& n = getRigidBodyNode(body);
        uStart = n.getUIndex();
        nu     = n.getDOF();
    }

    // Access to Instance variables. //

    const MassProperties& getBodyMassProperties(const State& s, MobilizedBodyIndex b) const {
        return getInstanceVars(s).bodyMassProperties[b];
    }
    const Transform& getMobilizerFrame(const State& s, MobilizedBodyIndex b) const {
        return getInstanceVars(s).outboardMobilizerFrames[b];
    }
    const Transform& getMobilizerFrameOnParent(const State& s, MobilizedBodyIndex b) const {
        return getInstanceVars(s).inboardMobilizerFrames[b];
    }

    MassProperties& updBodyMassProperties(State& s, MobilizedBodyIndex b) const {
        return updInstanceVars(s).bodyMassProperties[b];
    }
    Transform& updMobilizerFrame(State& s, MobilizedBodyIndex b) const {
        return updInstanceVars(s).outboardMobilizerFrames[b];
    }
    Transform& updMobilizerFrameOnParent(State& s, MobilizedBodyIndex b) const {
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

    // Extract position, velocity, and acceleration information for MobilizedBodies out of the 
    // State cache. Normally the State already has to have been realized to the Position, Velocity,
    // or Acceleration stage resp., but if you are calling one of these *during* realization of
    // that stage (say, from a Constraint error-calculation routine) then set the indicated flag
    // and the routine will allow you to retrieve the cache value one stage earlier. In that case you are on
    // the honor system not to reference the value before it is available.
    const Transform&  getBodyTransform   (const State&, MobilizedBodyIndex, bool realizingPosition=false) const;
    const SpatialVec& getBodyVelocity    (const State&, MobilizedBodyIndex, bool realizingVelocity=false) const;
    const SpatialVec& getBodyAcceleration(const State&, MobilizedBodyIndex, bool realizingAcceleration=false) const;

    // velocity dependent
    const SpatialVec& getCoriolisAcceleration     (const State&, MobilizedBodyIndex) const;
    const SpatialVec& getTotalCoriolisAcceleration(const State&, MobilizedBodyIndex) const;
    const SpatialVec& getGyroscopicForce          (const State&, MobilizedBodyIndex) const;
    const SpatialVec& getCentrifugalForces        (const State&, MobilizedBodyIndex) const;

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

    bool projectQConstraints(State& s, Real consAccuracy, const Vector& yWeights,
							 const Vector& ooTols, Vector& yErrest, System::ProjectOptions opts) const
	{
        // TODO
        enforcePositionConstraints(s, consAccuracy, yWeights, ooTols, yErrest, opts);
        return true;
    }
    bool projectUConstraints(State& s, Real consAccuracy, const Vector& yWeights,
							 const Vector& ooTols, Vector& yErrest, System::ProjectOptions opts) const
	{ 
		// TODO
        enforceVelocityConstraints(s, consAccuracy, yWeights, ooTols, yErrest, opts);
        return true;
    }


        // OPERATORS //

    Real calcKineticEnergy(const State&) const;

    // Calculate the product J*v where J is the kinematic Jacobian dV/du=~Phi*~H and
    // v is a vector in mobility space (internal coordinates). If v==u, that is,
    // it is a vector of generalized speeds, then the result will be V_GB, the
    // spatial velocity of each body.
    // This is an O(N) operator which can be called once the State is realized
    // to Stage::Position or higher. Because this is an operator, there is
    // no effect on the State cache.
    void calcSpatialKinematicsFromInternal(const State&,
        const Vector&        v,
        Vector_<SpatialVec>& Jv) const;

    // Calculate the product ~J*X where J is the partial velocity Jacobian dV/du
    // (~J=H*Phi)and X is a vector of SpatialVec's, one per body. See Eq. 76&77 in
    // Schwieters' paper, and see 81a & b for a use of this routine to compute
    // energy gradient in internal coordinates. In that case X=dE/dR, that is
    // the gradient of the energy w.r.t. atomic positions, summed and shifted
    // to body origins. There we are pretending dR/dq is the same as dV/du, which
    // will be true if dq/dt = u, which works for all cases except quaternions.
    // Schwieters handles that by using Euler angles for orientation coordinates
    // when doing minimizations. But note that the routine works in terms of u,
    // not q, so it produces a meaningful result in all cases, just not one that
    // can be mapped directly back to quaternion coordinates. This is an O(n)
    // operator which can be called after realizePosition().
    // Because this is an operator, there is no effect on the State cache.
    void calcInternalGradientFromSpatial(const State&, 
        const Vector_<SpatialVec>& X, 
        Vector&                    JX) const;

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


    // Must be in Stage::Position to calculate out_q = Q(q)*in_u (e.g., qdot=Q*u)
    // or out_u = in_q * Q(q). Note that one of "in" and "out" is always "q-like" while
    // the other is "u-like", but which is which changes if the matrix is on the right.
    // This is an O(N) operator since Q is block diagonal.
    void multiplyByQMatrix(const State& s, bool matrixOnRight, const Vector& in, Vector& out) const;

    // Must be in Stage::Position to calculate out_u = QInv(q)*in_q (e.g., u=QInv*qdot)
    // or out_q = in_u * QInv(q). Note that one of "in" and "out" is always "q-like" while
    // the other is "u-like", but which is which changes if the matrix is on the right.
    // This is an O(N) operator since QInv is block diagonal.
    void multiplyByQMatrixInverse(const State& s, bool matrixOnRight, const Vector& in, Vector& out) const;

    // Must be in Stage::Position to calculate qdot = Q*u.
    void calcQDot(const State& s,
        const Vector& u,
        Vector&       qdot) const;

    // Must be in Stage::Velocity to calculate qdotdot = Qdot*u + Q*udot.
    void calcQDotDot(const State& s,
        const Vector& udot,
        Vector&       qdotdot) const;

        // MOBILIZER OPERATORS //

    // These operators deal with an isolated mobilizer and are thus independent of 
    // any other generalized coordinates or speeds.

    // State must be realized to Stage::Position, so that we can extract Q(q) from it to calculate
    // qdot=Q(q)*u for this mobilizer.
    void calcMobilizerQDotFromU(const State&, MobilizedBodyIndex, int nu, const Real* u, 
                                int nq, Real* qdot) const;

    // State must be realized to Stage::Velocity, so that we can extract Q(q), QDot(q,u), and u from it to calculate
    // qdotdot=Q(q)*udot + QDot(q,u)*u for this mobilizer.
    void calcMobilizerQDotDotFromUDot(const State&, MobilizedBodyIndex, int nu, const Real* udot, 
                                      int nq, Real* qdotdot) const;

    // State must be realized through Stage::Instance. Neither the State nor its
    // cache are modified by this method, since it is an operator.
    // The number of q's is passed in as a sanity check, to make sure the caller
    // and the called mobilizer agree on the generalized coordinates.
    // Returns X_FM(q).
    Transform calcMobilizerTransformFromQ(const State&, MobilizedBodyIndex, int nq, const Real* q) const;

    // State must be realized through Stage::Position. Neither the State nor its
    // cache are modified by this method, since it is an operator.
    // The number of u's is passed in as a sanity check, to make sure the caller
    // and the called mobilizer agree on the generalized speeds.
    // Returns V_FM(q,u)=H_FM(q)*u, where the q dependency is extracted from the State via
    // the hinge transition matrix H_FM(q).
    SpatialVec calcMobilizerVelocityFromU(const State&, MobilizedBodyIndex, int nu, const Real* u) const;

    // State must be realized through Stage::Velocity. Neither the State nor its
    // cache are modified by this method, since it is an operator.
    // The number of u's (and udot's) is passed in as a sanity check, to make sure the caller
    // and the called mobilizer agree on the generalized accelerations.
    // Returns A_FM(q,u,udot)=H_FM(q)*udot + HDot_FM(q,u)*u where the q and u dependencies
    // are extracted from the State via H_FM(q), and HDot_FM(q,u).
    SpatialVec calcMobilizerAccelerationFromUDot(const State&, MobilizedBodyIndex, int nu, const Real* udot) const;

    // These perform the same computations as above but then transform the results so that they
    // relate the child body's frame B to its parent body's frame P, rather than the M and F frames
    // which are attached to B and P respectively but differ by a constant transform.
    Transform  calcParentToChildTransformFromQ(const State& s, MobilizedBodyIndex mb, int nq, const Real* q) const;
    SpatialVec calcParentToChildVelocityFromU (const State& s, MobilizedBodyIndex mb, int nu, const Real* u) const;
    SpatialVec calcParentToChildAccelerationFromUDot(const State& s, MobilizedBodyIndex mb, int nu, const Real* udot) const;


    void setDefaultModelValues       (const SBTopologyCache&, SBModelVars&)        const;
    void setDefaultInstanceValues    (const SBModelVars&,     SBInstanceVars&)     const;
    void setDefaultTimeValues        (const SBModelVars&,     SBTimeVars&)         const;
    void setDefaultPositionValues    (const SBModelVars&,     Vector& q)           const;
    void setDefaultVelocityValues    (const SBModelVars&,     Vector& u)           const;
    void setDefaultDynamicsValues    (const SBModelVars&,     SBDynamicsVars&)     const;
    void setDefaultAccelerationValues(const SBModelVars&,     SBAccelerationVars&) const;


    // A single Constraint may generate multiple of these constraint equations.
    // TODO: this is part of the obsolete IVM constraint system
    int getNDistanceConstraints()     const {return distanceConstraints.size();}

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

    int getQIndex(MobilizedBodyIndex) const;
    int getQAlloc(MobilizedBodyIndex) const;
    int getUIndex(MobilizedBodyIndex) const;
    int getDOF   (MobilizedBodyIndex) const;

    // Modeling info.

    void setUseEulerAngles(State& s, bool useAngles) const;
    void setMobilizerIsPrescribed(State& s, MobilizedBodyIndex, bool prescribe) const;
    void setConstraintIsDisabled(State& s, ConstraintIndex constraint, bool disabled) const;
    bool getUseEulerAngles(const State& s) const;
    bool isMobilizerPrescribed(const State& s, MobilizedBodyIndex) const;
    bool isConstraintDisabled(const State& s, ConstraintIndex constraint) const;
    void convertToEulerAngles(const State& inputState, State& outputState) const;
    void convertToQuaternions(const State& inputState, State& outputState) const;

        // CALLABLE AFTER realizeModel()

    int  getNQuaternionsInUse(const State&) const;						 // mquat
    bool isUsingQuaternion(const State&, MobilizedBodyIndex) const;
    QuaternionPoolIndex getQuaternionPoolIndex(const State&, MobilizedBodyIndex) const; // Invalid if none
    AnglePoolIndex      getAnglePoolIndex     (const State&, MobilizedBodyIndex) const; // Invalid if none

	// Note that although holonomic constraints are position-level constraints, they
	// do *not* include quaternion constraints (although the state's QErr vector does
	// include both). The total number of position-level constraints is thus
	// getNHolonomicConstraintEquationsInUse()+getNQuaternionsInUse()==mp+mquat.

	int getNHolonomicConstraintEquationsInUse       (const State&) const; // mh
	int getNNonholonomicConstraintEquationsInUse    (const State&) const; // mn
	int getNAccelerationOnlyConstraintEquationsInUse(const State&) const; // ma

	void calcHolonomicConstraintMatrixPQInverse(const State&, Matrix&) const; // mh X nq
	void calcHolonomicVelocityConstraintMatrixP(const State&, Matrix&) const; // mh X nu
	void calcHolonomicVelocityConstraintMatrixPt(const State&, Matrix&) const; // nu X mh
	void calcNonholonomicConstraintMatrixV     (const State&, Matrix&) const; // mn X nu
	void calcNonholonomicConstraintMatrixVt    (const State&, Matrix&) const; // nu X mn
	void calcAccelerationOnlyConstraintMatrixA (const State&, Matrix&) const; // ma X nu
	void calcAccelerationOnlyConstraintMatrixAt(const State&, Matrix&) const; // nu X ma

    // Treating all constraints together, given a comprehensive set of multipliers lambda,
    // generate the complete set of body and mobility forces applied by all the 
    // constraints.
	void calcConstraintForcesFromMultipliers
      (const State& s, const Vector& lambda,
	   Vector_<SpatialVec>& bodyForcesInG,
	   Vector&              mobilityForces) const;

    // Call after realizeDynamics()
    const SpatialMat& getArticulatedBodyInertia(const State& s, MobilizedBodyIndex) const;


    // This is a solver which generates internal velocities from spatial ones.
    void velFromCartesian(const Vector& pos, Vector& vel) {assert(false);/*TODO*/}

    void enforcePositionConstraints(State& s, Real consAccuracy, const Vector& yWeights,
									const Vector& ooTols, Vector& yErrest, System::ProjectOptions) const;
    void enforceVelocityConstraints(State& s, Real consAccuracy, const Vector& yWeights,
									const Vector& ooTols, Vector& yErrest, System::ProjectOptions) const;

    // Unconstrained (tree) dynamics 

    // articulated body inertias
    void calcArticulatedBodyInertias(const State&) const;

     // articulated body remainder forces
    void calcZ(const State&, 
        const Vector&              mobilityForces,
        const Vector_<SpatialVec>& bodyForces) const;
    void calcTreeAccel(const State&) const;                // accels with forces from last calcZ

    void fixVel0(State&, Vector& vel) const; // TODO -- yuck

    // Part of constrained dynamics (TODO -- more to move here)
    void calcY(const State&) const;

    // Pass in internal forces in T; they will be adjusted by this routine.
    void calcConstraintCorrectedInternalForces(const State&, Vector& T); 

    const RigidBodyNode& getRigidBodyNode(int nodeNum) const {
        const RigidBodyNodeIndex& ix = nodeNum2NodeMap[nodeNum];
        return *rbNodeLevels[ix.level][ix.offset];
    }

    // Add a distance constraint equation and assign it a particular slot in
    // the qErr, uErr, and multiplier arrays.
    // Return the assigned distance constraint index for caller's use.
    //TODO: OBSOLETE
    int addOneDistanceConstraintEquation(
        const RBStation& s1, const RBStation& s2, const Real& d);

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
            (s.getCacheEntry(getMySubsystemIndex(),topologyCacheIndex)).get();
    }
    SBTopologyCache& updTopologyCache(const State& s) const { //mutable
        assert(subsystemTopologyHasBeenRealized() && topologyCacheIndex >= 0);
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

    const SBPositionCache& getPositionCache(const State& s, bool realizingPosition=false) const {
        const AbstractValue& cacheEntry = 
            realizingPosition ? (const AbstractValue&)s.updCacheEntry(getMySubsystemIndex(),getModelCache(s).qCacheIndex)
                              : s.getCacheEntry(getMySubsystemIndex(),getModelCache(s).qCacheIndex);
        return Value<SBPositionCache>::downcast(cacheEntry).get();
    }
    SBPositionCache& updPositionCache(const State& s) const { //mutable
        return Value<SBPositionCache>::downcast
            (s.updCacheEntry(getMySubsystemIndex(),getModelCache(s).qCacheIndex)).upd();
    }

    const SBVelocityCache& getVelocityCache(const State& s, bool realizingVelocity=false) const {
        const AbstractValue& cacheEntry = 
            realizingVelocity ? (const AbstractValue&)s.updCacheEntry(getMySubsystemIndex(),getModelCache(s).uCacheIndex)
                              : s.getCacheEntry(getMySubsystemIndex(),getModelCache(s).uCacheIndex);
        return Value<SBVelocityCache>::downcast(cacheEntry).get();
    }
    SBVelocityCache& updVelocityCache(const State& s) const { //mutable
        return Value<SBVelocityCache>::downcast
            (s.updCacheEntry(getMySubsystemIndex(),getModelCache(s).uCacheIndex)).upd();
    }

    const SBDynamicsCache& getDynamicsCache(const State& s, bool realizingDynamics=false) const {
        const AbstractValue& cacheEntry = 
            realizingDynamics ? (const AbstractValue&)s.updCacheEntry(getMySubsystemIndex(),getModelCache(s).dynamicsCacheIndex)
                              : s.getCacheEntry(getMySubsystemIndex(),getModelCache(s).dynamicsCacheIndex);
        return Value<SBDynamicsCache>::downcast(cacheEntry).get();
    }
    SBDynamicsCache& updDynamicsCache(const State& s) const { //mutable
        return Value<SBDynamicsCache>::downcast
            (s.updCacheEntry(getMySubsystemIndex(),getModelCache(s).dynamicsCacheIndex)).upd();
    }

    const SBAccelerationCache& getAccelerationCache(const State& s, bool realizingAcceleration=false) const {
        const AbstractValue& cacheEntry = 
            realizingAcceleration ? (const AbstractValue&)s.updCacheEntry(getMySubsystemIndex(),getModelCache(s).accelerationCacheIndex)
                                  : s.getCacheEntry(getMySubsystemIndex(),getModelCache(s).accelerationCacheIndex);
        return Value<SBAccelerationCache>::downcast(cacheEntry).get();
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

    const SimbodyMatterSubsystem& getMySimbodyMatterSubsystemHandle() const {
        return SimbodyMatterSubsystem::downcast(getOwnerSubsystemHandle());
    }
    SimbodyMatterSubsystem& updMySimbodyMatterSubsystemHandle() {
        return SimbodyMatterSubsystem::updDowncast(updOwnerSubsystemHandle());
    }
    bool getShowDefaultGeometry() const;
    void setShowDefaultGeometry(bool show);

private:
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
    
    // Recursively calculate the weights for Q's or U's.
    void calcQUnitWeightsRecursively(const State& s, State& tempState, Vector& weights, Vec6& bounds, const RigidBodyNode& body) const;
    void calcUUnitWeightsRecursively(const State& s, State& tempState, Vector& weights, Vec6& bounds, const RigidBodyNode& body) const;

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
    // are adopted. The MobilizedBodyIndex (converted to int) is the index of a
    // MobilizedBody in this array.
    std::vector<MobilizedBody*> mobilizedBodies;

    // Constraints are treated similarly.
    std::vector<Constraint*>    constraints;

    // Our realizeTopology method calls this after all bodies & constraints have been added,
    // to construct part of the topology cache below.
    void endConstruction(State&);

        // TOPOLOGY CACHE

    // The data members here are filled in when realizeTopology() is called.
    // The flag which remembers whether we have realized topology is in 
    // the Subsystem::Guts base class.
    // Note that a cache is treated as mutable, so the methods that manipulate
    // it are const.

    void clearTopologyCache() const {
        SimbodyMatterSubsystemRep& mthis = const_cast<SimbodyMatterSubsystemRep&>(*this);
        mthis.clearTopologyCache(); // call the non-const version
    }

    // This method should clear out all the data members below.
    void clearTopologyCache();

		// Mobilized bodies and their rigid body nodes

    // This holds pointers to nodes and serves to map (level,offset) to nodeNum.
    std::vector<RBNodePtrList>      rbNodeLevels;
    // Map nodeNum to (level,offset).
    std::vector<RigidBodyNodeIndex> nodeNum2NodeMap;

		// Constraints

    // Here we sort the above constraints by branch (ancestor's base body), then by
    // level within that branch. That is, each constraint is addressed
    // by three indices [branch][levelOfAncestor][offset]
    // where offset is an arbitrary unique integer assigned to all
    // the constraints on the same branch with the same level of ancestor.
    std::vector< std::vector< std::vector<ConstraintIndex> > > branches;

    // Partition the constraints into groups which are coupled by constraints at
    // the indicated level. Only Constraints which generate holonomic constraint
    // equations (mp>0) are coupled at the position level. Constraints with
    // *either* holonomic or nonholonomic (mv>0) constraint equations can be
    // coupled at the velocity level because the derivatives of the holonomic
    // constraints are velocity constraints. And any Constraints can be coupled
    // at the accleration level.
    //
    // Each ConstraintSet contains a Subtree comprising all the relevant
    // mobilized bodies for that set. Note that this is *kinematic* coupling;
    // for dynamics there is additional coupling brought in by the mass
    // matrix in the computation of (G M^-1 G^T).

    // sorted in nondecreasing order of ancestor MobilizedBodyIndex
    std::vector<CoupledConstraintSet> positionCoupledConstraints;     // for P
    std::vector<CoupledConstraintSet> velocityCoupledConstraints;     // for PV
    std::vector<CoupledConstraintSet> accelerationCoupledConstraints; // for G=PVA

    // This further partitions the accelerationCoupledConstraints into
    // larger groups comprised of accelerationCoupledConstraints whose ancestor
    // bodies' inboard paths share a body other than ground. That is, we couple
    // the constraints unless they involve completely disjoint grounded subtrees.
    // These groups correspond to disjoint blocks in M (as
    // well as G, of course), so we can decouple them for the (G M^-1 G^T) 
    // calculation.
    // TODO: acceleration constraints should instead be dealt with recursively,
    // based on *kinematic* coupling; where the kinematically coupled groups are
    // used to modify the articulated body inertias.
    std::vector<CoupledConstraintSet> dynamicallyCoupledConstraints;


    // TODO: these state indices and counters should be deferred to realizeModel()
    // so we can have Model stage variables which change the number of state
    // slots needed.

    // Initialize to 0 at beginning of construction. These are for doling
    // out Q & U state variables to the nodes.
    UIndex        nextUSlot;
    USquaredIndex nextUSqSlot;
    QIndex        nextQSlot;

    // These are similarly for doling out slots for *topological* constraint
    // equations for position errors, velocity errors, and multipliers.
    // Note that quaternion normalization constraints are *not* topological
    // so aren't counted here (although we begin assigning them slots at
    // the end of the qErr's used up here).
    int nextQErrSlot;
    int nextUErrSlot;
    int nextMultSlot;

    // Dole out slots for all Constrained Bodies which belong to Constraints
    // whose Ancestor is not Ground (except for the Ancestor bodies themselves).
    AncestorConstrainedBodyPoolIndex nextAncestorConstrainedBodyPoolSlot;

    int DOFTotal;   // summed over all nodes
    int SqDOFTotal; // sum of squares of ndofs per node
    int maxNQTotal; // sum of dofs with room for quaternions

    SBTopologyCache topologyCache;
    int topologyCacheIndex; // topologyCache is copied here in the State



    // These are the distant constraint equations and point-in-plane
    // constraint equations generated by Constraints.
	// TODO: these are OBSOLETE; part of the old IVM constraint system
    std::vector<RBDistanceConstraint*>     distanceConstraints;
    LengthConstraints* lConstraints;
    
    // Specifies whether default decorative geometry should be shown.
    bool showDefaultGeometry;
};

std::ostream& operator<<(std::ostream&, const SimbodyMatterSubsystemRep&);


#endif // SimTK_SIMBODY_MATTER_SUBSYSTEM_REP_H_
