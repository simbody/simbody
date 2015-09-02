#ifndef SimTK_SIMBODY_MATTER_SUBSYSTEM_REP_H_
#define SimTK_SIMBODY_MATTER_SUBSYSTEM_REP_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-15 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors: Derived from IVM code written by Charles Schwieters          *
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
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/SimbodyMatterSubtree.h"
#include "simbody/internal/MobilizedBody.h"
#include "simbody/internal/MobilizedBody_Ground.h"

#include "SimbodyTreeState.h"
#include "RigidBodyNode.h"

#include <set>
#include <map>
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

using namespace SimTK;

typedef Array_<const RigidBodyNode*>   RBNodePtrList;
typedef Vector_<SpatialVec>            SpatialVecList;

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
        coupledConstraints.resize((unsigned)constraints.size());
        std::set<ConstraintIndex>::const_iterator i = constraints.begin();
        for (int nxt=0; i != constraints.end(); ++i, ++nxt)
            coupledConstraints[nxt] = *i;
        topologyRealized = true;
    }

    const Array_<ConstraintIndex>& getCoupledConstraints() const {
        assert(topologyRealized);
        return coupledConstraints;
    }

private:
    // TOPOLOGY STATE
    std::set<ConstraintIndex> constraints;

    // TOPOLOGY CACHE
    bool topologyRealized;

    // Sorted in nondecreasing order of ancestor MobilizedBodyIndex.
    Array_<ConstraintIndex>    coupledConstraints;
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
      : Subsystem::Guts("SimbodyMatterSubsystem", "0.7.1")
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
       (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const;

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
    // MobilizedBody handle (leaving child as a non-owner reference), and 
    // makes it a child (outboard body) of the indicated parent. The new 
    // child body's id is returned, and will be greater than the parent's id.

    MobilizedBodyIndex adoptMobilizedBody(MobilizedBodyIndex parentIndex, MobilizedBody& child);
    int getNumMobilizedBodies() const {return (int)mobilizedBodies.size();}

    const MobilizedBody& getMobilizedBody(MobilizedBodyIndex ix) const {
        SimTK_INDEXCHECK(ix, (int)mobilizedBodies.size(),
                         "SimbodyMatterSubsystem::getMobilizedBody()");
        assert(mobilizedBodies[ix]);
        return *mobilizedBodies[ix];
    }

    // Note that we do not invalidate the subsystem topology cache yet, even
    // though we're handing out a writable reference to a MobilizedBody here.
    // Otherwise every time someone references Ground, e.g., by calling
    // matterSubsys.Ground() the subsystem would have its topology marked
    // invalid. For this reason, and also because the main program normally
    // retains a writable reference to MobilizedBodies, it is essential that 
    // every non-const method of MobilizedBody and its many descendants mark 
    // the subsystem topology invalid when called.
    MobilizedBody& updMobilizedBody(MobilizedBodyIndex ix) {
        SimTK_INDEXCHECK(ix, (int)mobilizedBodies.size(),
                         "SimbodyMatterSubsystem::updMobilizedBody()");
        assert(mobilizedBodies[ix]);
        return *mobilizedBodies[ix]; // topology not marked invalid yet
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

    const Constraint& getConstraint(ConstraintIndex ix) const {
        SimTK_INDEXCHECK(ix, (int)constraints.size(),
                         "SimbodyMatterSubsystem::getConstraint()");
        assert(constraints[ix]);
        return *constraints[ix];
    }
    Constraint& updConstraint(ConstraintIndex ix) {
        SimTK_INDEXCHECK(ix, (int)constraints.size(),
                         "SimbodyMatterSubsystem::updConstraint()");
        assert(constraints[ix]);
        return *constraints[ix];
    }

    UnilateralContactIndex adoptUnilateralContact(UnilateralContact*);
    int getNumUnilateralContacts() const 
    {   return (int)uniContacts.size(); }
    const UnilateralContact& 
    getUnilateralContact(UnilateralContactIndex ix) const {
        SimTK_INDEXCHECK(ix, (int)uniContacts.size(),
                         "SimbodyMatterSubsystem::getUnilateralContact()");
        assert(uniContacts[ix]);
        return *uniContacts[ix];
    }
    UnilateralContact& 
    updUnilateralContact(UnilateralContactIndex ix) {
        SimTK_INDEXCHECK(ix, (int)uniContacts.size(),
                         "SimbodyMatterSubsystem::updUnilateralContact()");
        assert(uniContacts[ix]);
        return *uniContacts[ix];
    }

    StateLimitedFrictionIndex adoptStateLimitedFriction(StateLimitedFriction*);
    int getNumStateLimitedFrictions() const 
    {   return (int)stateLtdFriction.size(); }

    const StateLimitedFriction& 
    getStateLimitedFriction(StateLimitedFrictionIndex ix) const {
        SimTK_INDEXCHECK(ix, (int)stateLtdFriction.size(),
                         "SimbodyMatterSubsystem::getStateLimitedFriction()");
        assert(stateLtdFriction[ix]);
        return *stateLtdFriction[ix];
    }
    StateLimitedFriction& 
    updStateLimitedFriction(StateLimitedFrictionIndex ix) {
        SimTK_INDEXCHECK(ix, (int)stateLtdFriction.size(),
                         "SimbodyMatterSubsystem::updStateLimitedFriction()");
        assert(stateLtdFriction[ix]);
        return *stateLtdFriction[ix];
    }

    // Topology stage cache entry
    AncestorConstrainedBodyPoolIndex 
    allocateNextAncestorConstrainedBodyPoolSlot() const {
        const AncestorConstrainedBodyPoolIndex nxt = nextAncestorConstrainedBodyPoolSlot;
        // Make this mutable briefly.
        ++ const_cast<SimbodyMatterSubsystemRep*>(this)->nextAncestorConstrainedBodyPoolSlot;
        return nxt;
    }


    // These counts can be obtained even during construction, where they
    // just return the current counts.
    // NumBodies includes ground.
    int getNumBodies()      const {return mobilizedBodies.size();}
    int getNumParticles()   const {return 0;} // TODO
    int getNumMobilities()  const {return getTotalDOF();}
    int getNumConstraints() const {return constraints.size();}
    MobilizedBodyIndex getParent(MobilizedBodyIndex) const;
    Array_<MobilizedBodyIndex> getChildren(MobilizedBodyIndex) const;

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

    // Extract position, velocity, and acceleration information for 
    // MobilizedBodies out of the State cache, accessing only the Tree cache
    // entries at each level, which should already have been marked valid.
    const Transform&  getBodyTransform   (const State&, MobilizedBodyIndex) const;
    const SpatialVec& getBodyVelocity    (const State&, MobilizedBodyIndex) const;
    const SpatialVec& getBodyAcceleration(const State&, MobilizedBodyIndex) const;

    // Call at Position stage or later. If necessary, composite body inertias 
    // will be realized first.
    const Array_<SpatialInertia,MobilizedBodyIndex>& 
    getCompositeBodyInertias(const State& s) const {
        realizeCompositeBodyInertias(s);
        return getCompositeBodyInertiaCache(s).compositeBodyInertia;
    }

    // Call at Position stage or later. If necessary, articulated body 
    // inertias will be realized first.
    const Array_<ArticulatedInertia,MobilizedBodyIndex>& 
    getArticulatedBodyInertias(const State& s) const {
        realizeArticulatedBodyInertias(s);
        return getArticulatedBodyInertiaCache(s).articulatedBodyInertia;
    }

    const Array_<ArticulatedInertia,MobilizedBodyIndex>& 
    getArticulatedBodyInertiasPlus(const State& s) const {
        realizeArticulatedBodyInertias(s);
        return getArticulatedBodyInertiaCache(s).pPlus;
    }

    // Call at Acceleration stage only.
    const Array_<SpatialVec,MobilizedBodyIndex>& 
    getArticulatedBodyForces(const State& s) const {
        return getTreeAccelerationCache(s).z;
    }
    const Array_<SpatialVec,MobilizedBodyIndex>& 
    getArticulatedBodyForcesPlus(const State& s) const {
        return getTreeAccelerationCache(s).zPlus;
    }

    // velocity dependent
    const SpatialVec& getMobilizerCoriolisAcceleration(const State&, MobilizedBodyIndex) const;
    const SpatialVec& getTotalCoriolisAcceleration    (const State&, MobilizedBodyIndex) const;
    const SpatialVec& getGyroscopicForce              (const State&, MobilizedBodyIndex) const;
    const SpatialVec& getMobilizerCentrifugalForces   (const State&, MobilizedBodyIndex) const;
    const SpatialVec& getTotalCentrifugalForces       (const State&, MobilizedBodyIndex) const;

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

    // This is a solver that sets continuous state variables q or u
    // (depending on stage) to their prescribed values that will already
    // have been computed. Note that prescribed udot=udot(t,q,u) is not
    // dealt with here because it does not involve a state change.
    // Returns true if it makes any changes.
    bool prescribeQ(State& s) const;
    bool prescribeU(State& s) const;

    // Project quaternions onto their constraint manifold by normalizing
    // them. Also removes any error component along the length of the
    // quaternions if given a qErrest. Returns true if any change was made.
    // This is used by projectQ().
    bool normalizeQuaternions(State& s, Vector& qErrest) const;

    int projectQ(State& s, Vector& qErrest, 
                 const ProjectOptions& opts,
                 ProjectResults& results) const;
    int projectU(State& s, Vector& uErrest, 
                 const ProjectOptions& opts,
                 ProjectResults& results) const;

        // REALIZATIONS //


    // Call at Instance Stage or later.
    void realizePositionKinematics(const State&) const;

    // Call at Instance + PositionKinematics Stage or later.
    void realizeVelocityKinematics(const State&) const;

    // Call at Position Stage or later.
    void realizeCompositeBodyInertias(const State&) const;

    // Call at Position Stage or later.
    void realizeArticulatedBodyInertias(const State&) const;

    bool isPositionKinematicsRealized(const State&) const;
    bool isVelocityKinematicsRealized(const State&) const;

    // These are just used in timing tests.
    void invalidatePositionKinematics(const State&) const;
    void invalidateVelocityKinematics(const State&) const;
    void invalidateCompositeBodyInertias(const State&) const;
    void invalidateArticulatedBodyInertias(const State&) const;

        // OPERATORS //

    Real calcKineticEnergy(const State&) const;

    void calcCompositeBodyInertias(const State&,
        Array_<SpatialInertia,MobilizedBodyIndex>& R) const;

    // Calculate the product J*v where J is the kinematic Jacobian 
    // dV/du=~Phi*~H (Schwieters' and Jain's terminology; our H is transposed
    // from theirs), and  v is a vector in mobility space (internal 
    // coordinates). If v==u, that is, it is a vector of generalized speeds, 
    // then the result will be V_GB, the spatial velocity of each body.
    // This is an O(N) operator which can be called once the State is realized
    // to Stage::Position or higher. Because this is an operator, there is
    // no effect on the State cache.
    void multiplyBySystemJacobian(const State&,
        const Vector&        v,
        Vector_<SpatialVec>& Jv) const;

    // Calculate the product ~J*X where J is the partial velocity Jacobian 
    // dV/du (~J=H*Phi)and X is a vector of force-space SpatialVec's, one per 
    // body. See Eq. 76&77 in Schwieters' paper, and see 81a & b for a use of 
    // this routine to compute energy gradient in internal coordinates. In that
    // case X=dE/dR, that is the gradient of the energy w.r.t. atomic positions, 
    // summed and shifted to body origins. There we are pretending dR/dq is the
    // same as dV/du, which will be true if dq/dt = u. In general, we have
    // dR/dq = (dV/du)*N^-1, where dq/dt = N*u (i.e., N= d qdot/du). But note 
    // that this method works in terms of u, not q, so it produces a meaningful
    // result in all cases, just not one that can be mapped directly back to 
    // generalized coordinates q. This is an O(n) operator which can be called 
    // after realizePosition(). Because this is an operator, there is no effect
    // on the State cache.
    void multiplyBySystemJacobianTranspose(const State&, 
        const Vector_<SpatialVec>& X, 
        Vector&                    JtX) const;

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
        const Array_<Real>&        presUDots, // packed
        Vector&                    netHingeForces,
        Array_<SpatialVec,MobilizedBodyIndex>& abForcesZ, 
        Array_<SpatialVec,MobilizedBodyIndex>& abForcesZPlus, 
        Vector_<SpatialVec>&       A_GB,
        Vector&                    udot, // in/out (in for prescribed udots)
        Vector&                    qdotdot,
        Vector&                    tau) const; 

    // Multiply by the mass matrix in O(n) time.
    void multiplyByM(const State& s,
        const Vector&             a,
        Vector&                   Ma) const;

    // Multiply by the mass matrix inverse in O(n) time. Works only with the
    // non-prescribed submatrix Mrr of M; entries f_p in f are not accessed,
    // and entries MInvf_p in MInvf are not written.
    void multiplyByMInv(const State&    s,
        const Vector&                   f,
        Vector&                         MInvf) const; 

    // Calculate the mass matrix in O(n^2) time. State must have already
    // been realized to Position stage. M must be resizeable or already the
    // right size (nXn). The result is symmetric but the entire matrix is
    // filled in.
    void calcM(const State& s, Matrix& M) const;

    // Calculate the mass matrix inverse in O(n^2) time. State must have already
    // been realized to Position stage. MInv must be resizeable or already the
    // right size (nXn). The result is symmetric but the entire matrix is
    // filled in. Only the non-prescribed block Mrr is inverted; other elements
    // are not written.
    void calcMInv(const State& s, Matrix& MInv) const;

    void calcTreeResidualForces(const State&,
        const Vector&               appliedMobilityForces,
        const Vector_<SpatialVec>&  appliedBodyForces,
        const Vector&               knownUdot,
        Vector_<SpatialVec>&        A_GB,
        Vector&                     residualMobilityForces) const;



    // Must be in Stage::Position to calculate out_q = N(q)*in_u (e.g., qdot=N*u)
    // or out_u = in_q * N(q). Note that one of "in" and "out" is always "q-like" while
    // the other is "u-like", but which is which changes if the matrix is on the right.
    // This is an O(n) operator since N is block diagonal.
    void multiplyByN(const State& s, bool matrixOnRight, const Vector& in, Vector& out) const;

    // Must be in Stage::Position to calculate out_u = NInv(q)*in_q (e.g., u=NInv*qdot)
    // or out_q = in_u * NInv(q). Note that one of "in" and "out" is always "q-like" while
    // the other is "u-like", but which is which changes if the matrix is on the right.
    // This is an O(n) operator since NInv is block diagonal.
    void multiplyByNInv(const State& s, bool matrixOnRight, const Vector& in, Vector& out) const;

    // Must be in Stage::Velocity to calculate out_q = NDot(q,u)*in_u
    // or out_u = in_q * NDot(q,u). Note that one of "in" and "out" is always 
    // "q-like" while the other is "u-like", but which is which changes if the 
    // matrix is on the right. This is an O(n) operator since NDot is block diagonal.
    void multiplyByNDot(const State& s, bool matrixOnRight, const Vector& in, Vector& out) const;

    // Must be in Stage::Position to calculate qdot = N*u.
    void calcQDot(const State& s,
        const Vector& u,
        Vector&       qdot) const;

    // Must be in Stage::Velocity to calculate qdotdot = Ndot*u + N*udot.
    void calcQDotDot(const State& s,
        const Vector& udot,
        Vector&       qdotdot) const;

        // MOBILIZER OPERATORS //

    // These operators deal with an isolated mobilizer and are thus independent of 
    // any other generalized coordinates or speeds.

    // State must be realized to Stage::Position, so that we can extract N(q) from it to calculate
    // qdot=N(q)*u for this mobilizer.
    void calcMobilizerQDotFromU(const State&, MobilizedBodyIndex, int nu, const Real* u, 
                                int nq, Real* qdot) const;

    // State must be realized to Stage::Velocity, so that we can extract N(q), NDot(q,u), and u from it to calculate
    // qdotdot=N(q)*udot + NDot(q,u)*u for this mobilizer.
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

        // CALLABLE AFTER realizeTopology()

    int getTotalDOF()    const {assert(subsystemTopologyHasBeenRealized()); return DOFTotal;}
    int getTotalSqDOF()  const {assert(subsystemTopologyHasBeenRealized()); return SqDOFTotal;}
    int getTotalQAlloc() const {assert(subsystemTopologyHasBeenRealized()); return maxNQTotal;}

    int getNumTopologicalPositionConstraintEquations()     const {
        assert(subsystemTopologyHasBeenRealized()); return nextQErrSlot;
    }
    int getNumTopologicalVelocityConstraintEquations()     const {
        assert(subsystemTopologyHasBeenRealized()); return nextUErrSlot;
    }
    int getNumTopologicalAccelerationConstraintEquations() const {
        assert(subsystemTopologyHasBeenRealized()); return nextMultSlot;
    }

    int getQIndex(MobilizedBodyIndex) const;
    int getQAlloc(MobilizedBodyIndex) const;
    int getUIndex(MobilizedBodyIndex) const;
    int getDOF   (MobilizedBodyIndex) const;

    // Modeling info.

    void setUseEulerAngles(State& s, bool useAngles) const;
    void setConstraintIsDisabled(State& s, ConstraintIndex constraint, bool disabled) const;
    bool getUseEulerAngles(const State& s) const;
    bool isConstraintDisabled(const State& s, ConstraintIndex constraint) const;
    void convertToEulerAngles(const State& inputState, State& outputState) const;
    void convertToQuaternions(const State& inputState, State& outputState) const;

        // CALLABLE AFTER realizeModel()

    int  getNumQuaternionsInUse(const State&) const;                // mquat
    bool isUsingQuaternion(const State&, MobilizedBodyIndex) const;
    QuaternionPoolIndex getQuaternionPoolIndex(const State&, MobilizedBodyIndex) const; // Invalid if none

    // Note that although holonomic constraints are position-level constraints, they
    // do *not* include quaternion constraints (although the state's QErr vector does
    // include both). The total number of position-level constraints is thus
    // getNumHolonomicConstraintEquationsInUse()+getNumQuaternionsInUse()==mp+mquat.

    int getNumHolonomicConstraintEquationsInUse       (const State&) const; // mh
    int getNumNonholonomicConstraintEquationsInUse    (const State&) const; // mn
    int getNumAccelerationOnlyConstraintEquationsInUse(const State&) const; // ma

    void calcHolonomicConstraintMatrixPNInv    (const State&, Matrix&) const; // mh X nq
    void calcHolonomicVelocityConstraintMatrixP(const State&, Matrix&) const; // mh X nu
    void calcHolonomicVelocityConstraintMatrixPt(const State&, Matrix&) const; // nu X mh
    void calcNonholonomicConstraintMatrixV     (const State&, Matrix&) const; // mn X nu
    void calcNonholonomicConstraintMatrixVt    (const State&, Matrix&) const; // nu X mn
    void calcAccelerationOnlyConstraintMatrixA (const State&, Matrix&) const; // ma X nu
    void calcAccelerationOnlyConstraintMatrixAt(const State&, Matrix&) const; // nu X ma

    void calcMobilizerReactionForces(const State& s, Vector_<SpatialVec>& forces) const;
    // This alternative is for debugging and testing; it is slow but should
    // produce the same answers as calcMobilizerReactionForces().
    void calcMobilizerReactionForcesUsingFreebodyMethod(const State& s, Vector_<SpatialVec>& forces) const;



    // Constraint multipliers are known to the State object.
    const Vector& getConstraintMultipliers(const State& s) const
    {   return this->getMultipliers(s); }

    // But taus (prescribed motion "multipliers") are internal.
    const Vector& getMotionMultipliers(const State& s) const {
        const SBTreeAccelerationCache& tac = getTreeAccelerationCache(s);
        return tac.presMotionForces;
    }

    Vector calcMotionErrors(const State& state, const Stage& stage) const;

    void findMotionForces(const State&         s,
                          Vector&              mobilityForces) const;
    void findConstraintForces(const State&         s, 
                              Vector_<SpatialVec>& bodyForcesInG,
                              Vector&              mobilityForces) const;

    Real calcMotionPower(const State& s) const;
    Real calcConstraintPower(const State& s) const;

    // Treating all constraints together, given a comprehensive set of 
    // multipliers lambda, generate the complete set of body and mobility 
    // forces applied by all the constraints. Return the
    // individual Constraints' force contributions in the final two arguments.
    void calcConstraintForcesFromMultipliers
      (const State&         state, 
       const Vector&        lambda,
       Vector_<SpatialVec>& bodyForcesInG,
       Vector&              mobilityForces,
       Array_<SpatialVec>&  constrainedBodyForcesInG,
       Array_<Real>&        contraintMobilityForces) const;

    // Call this signature if you don't care about the individual constraint
    // contributions.
    void calcConstraintForcesFromMultipliers
      (const State&         state, 
       const Vector&        lambda,
       Vector_<SpatialVec>& bodyForcesInG,
       Vector&              mobilityForces) const
    {
        const SBInstanceCache& ic = getInstanceCache(state);
        const int ncb = ic.totalNConstrainedBodiesInUse;
        const int ncu = ic.totalNConstrainedUInUse;
        Array_<SpatialVec> constrainedBodyForcesInG(ncb);
        Array_<Real>       constraintMobilityForces(ncu);
        calcConstraintForcesFromMultipliers(state,lambda,bodyForcesInG,
            mobilityForces,constrainedBodyForcesInG,constraintMobilityForces);
    }

    // Form the product 
    //    fu = [ ~P ~V ~A ] * lambda
    // with all or a subset of P,V,A included. The multiplier-like vector 
    // must have length m=mp+mv+ma always. This is an O(n+m) method.
    void multiplyByPVATranspose(const State&     state,
                                bool             includeP,
                                bool             includeV,
                                bool             includeA,
                                const Vector&    lambda,
                                Vector&          fu) const;

    // Explicitly form the u-space constraint Jacobian transpose 
    // ~G=[~P ~V ~A] or selected submatrices of it. Performance is best if the 
    // output matrix has columns stored contiguously in memory, but this method
    // will work anyway, in that case using a contiguous temporary for column 
    // calculations and then copying out into the result. The matrix will be
    // resized as necessary to nu X (mp+mv+ma) where the constraint dimensions
    // can be zero if that submatrix is not selected.
    void calcPVATranspose(  const State&     state,
                            bool             includeP,
                            bool             includeV,
                            bool             includeA,
                            Matrix&          PVAt) const;

    // Form the product 
    //    fq = [ ~Pq ] * lambdap
    // The multiplier-like vector must have length m=mp always.
    // This is equivalent to ~(P*N^-1)*lambdap = ~N^-1 * (~P * lambdap).
    // This is an O(n+mp) method.
    void multiplyByPqTranspose(const State&     state,
                               const Vector&    lambdap,
                               Vector&          fq) const;

    // Explicitly form the q-space holonomic constraint Jacobian transpose
    // Pqt (= ~(N^-1) * ~P). Performance is best if the output matrix has 
    // columns stored contiguously in memory, but this method will work anyway,
    // in that case using a contiguous temporary for column calculations and 
    // then copying out into the result. The matrix will be resized as 
    // necessary to nq X mp.
    void calcPqTranspose(   const State&     state,
                            Matrix&          Pqt) const;

    // Calculate the bias vector from the constraint error
    // equations used in multiplyByPVA. Here bias is what you would get
    // when ulike==0. The output Vector must use contiguous storage. It will 
    // be resized if necessary to length m=mp+mv+ma.
    void calcBiasForMultiplyByPVA(const State& state,
                                  bool         includeP,
                                  bool         includeV,
                                  bool         includeA,
                                  Vector&      bias) const;

    // Calculate the bias vector from the acceleration constraint error
    // equations used. Here bias is what you would get from paerr, vaerr,
    // and aerr when udot==0. This is different than calcBiasForMultiplyByPVA()
    // because that method uses pverr for holonomic constraints rather than
    // paerr. The output Vector must use contiguous storage. It will 
    // be resized if necessary to length m=mp+mv+ma.
    void calcBiasForAccelerationConstraints(const State& state,
                                            bool         includeP,
                                            bool         includeV,
                                            bool         includeA,
                                            Vector&      bias) const;

    void calcBiasForMultiplyByPq(const State& state,
                                 Vector&      bias) const
    {   calcBiasForMultiplyByPVA(state,true,false,false,bias); }

    // Given a bias calculated by the above method using the same settings
    // for the "include" flags, form the product 
    //           [ P ]
    //    PVAu = [ V ] * ulike
    //           [ A ]
    // with all or a subset of P,V,A included. The u-like vector must have
    // length nu always. This is an O(n+m) method.
    void multiplyByPVA(const State&     state,
                       bool             includeP,
                       bool             includeV,
                       bool             includeA,
                       const Vector&    bias,
                       const Vector&    ulike,
                       Vector&          PVAu) const;

    // Explicitly form the u-space constraint Jacobian G=[P;V;A] or 
    // selected submatrices of it. Performance is best if the output matrix 
    // has columns stored contiguously in memory, but this method will work 
    // anyway, in that case using a contiguous temporary for column 
    // calculations and then copying out into the result. The matrix will be
    // resized as necessary to (mp+mv+ma) X nu where the constraint dimensions
    // can be zero if that submatrix is not selected.
    void calcPVA(   const State&     state,
                    bool             includeP,
                    bool             includeV,
                    bool             includeA,
                    Matrix&          PVA) const;

    // Given a bias calculated by the above method using just includeP=true
    // (or the leading bias_p segment of a complete bias vector), form the
    // product PqXqlike = Pq*qlike (= P*N^-1*qlike). The q-like vector must 
    // have length nq always and all Vectors must use contiguous storage.
    // This is an O(nq+mp) method.
    void multiplyByPq(  const State&   state,
                        const Vector&  bias_p,
                        const Vector&  qlike,
                        Vector&        PqXqlike) const;

    // Explicitly form the q-space holonomic constraint Jacobian Pq (= P*N^-1).
    // Performance is best if the output matrix has columns stored contiguously
    // in memory, but this method will work anyway, in that case using a 
    // contiguous temporary for column calculations and then copying out into 
    // the result. The matrix will be resized as necessary to mp X nq.
    void calcPq(    const State&     state,
                    Matrix&          Pq) const;

    // Calculate the mXm "projected mass matrix" G * M^-1 * G^T. By using
    // a combination of O(n) operators we can calculate this in O(m*n) time.
    // The method requires only O(n) memory also, except for the mXm result.
    // State must be realized through Velocity stage unless all constraints
    // are holonomic in which case Position will do. This matrix is used
    // when solving for Lagrange multipliers: (G M^-1 G^T) lambda = aerr
    // gives values for lambda that elimnate aerr.
    // Performance is best if the output matrix has columns stored 
    // contiguously in memory, but this method will work anyway, in that case
    // using a contiguous temporary for column calculations and then copying
    // out into the result.
    void calcGMInvGt(const State&   state,
                     Matrix&        GMInvGt) const;

    // Use factored GMInvGt to solve GMinvGt*impulse=deltaV. The main benefit
    // of this method is that it promises to use the same method Simbody does
    // to deal with constraint redundancies.
    void solveForConstraintImpulses(const State&     state,
                                    const Vector&    deltaV,
                                    Vector&          impulse) const;

    // Given an array of nu udots, return nb body accelerations in G (including
    // Ground as the 0th body with A_GB[0]=0). The returned accelerations are
    // A = J*udot + Jdot*u, with the Jdot*u (coriolis acceleration) term
    // extracted from the supplied state, which must have been realized to
    // Velocity stage.
    // The input and output Vectors must use contiguous storage.
    void calcBodyAccelerationFromUDot(const State&          state,
                                      const Vector_<Real>&  knownUDot,
                                      Vector_<SpatialVec>&  A_GB) const;

    // Given an array of nu udots and already-calculated corresponding 
    // qdotdot=N*udot + NDot*u, and the set of corresponding
    // Ground-relative body accelerations, compute the constraint
    // acceleration errors that result due to the constraints currently active
    // in the given state. All acceleration-level constraints are included:
    // holonomic second derivatives, nonholonomic first derivatives, and 
    // acceleratin-only constraints. This is a pure operator and does not 
    // affect the state or state cache. Vectors must use contiguous data.
    // State must have been realized through Velocity stage.
    void calcConstraintAccelerationErrors
       (const State&                state,
        const Vector_<SpatialVec>&  A_GB,
        const Vector_<Real>&        udot,
        const Vector_<Real>&        qdotdot,
        Vector&                     pvaerr) const;

    // This is a solver which generates internal velocities from spatial ones.
    void velFromCartesian(const Vector& pos, Vector& vel) {assert(false);/*TODO*/}

    void enforcePositionConstraints(State& s, Real consAccuracy, const Vector& yWeights,
                                    const Vector& ooTols, Vector& yErrest, ProjectOptions) const;
    void enforceVelocityConstraints(State& s, Real consAccuracy, const Vector& yWeights,
                                    const Vector& ooTols, Vector& yErrest, ProjectOptions) const;

    // Unconstrained (tree) dynamics methods for use during realization.

    // Part of OLD constrained dynamics; TODO: may be useful in op space inertia calcs.
    void realizeY(const State&) const;

    const RigidBodyNode& getRigidBodyNode(MobilizedBodyIndex nodeNum) const {
        const RigidBodyNodeIndex& ix = nodeNum2NodeMap[nodeNum];
        return *rbNodeLevels[ix.level][ix.offset];
    }

    //
    //    STATE ARCHEOLOGY
    //    These routines know where the bodies are buried (no pun intended).
    //

    // This first one just provides access to the local topologyCache variable
    // that is stored directly in the SimbodyMatterSubsystemRep. This is the one
    // that counts - there is a copy in the State but it is just for sanity
    // checking.
    const SBTopologyCache& getMatterTopologyCache() const 
    {   return topologyCache; }

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
            (s.getCacheEntry(getMySubsystemIndex(),topologyCache.instanceCacheIndex)).get();
    }
    SBInstanceCache& updInstanceCache(const State& s) const { //mutable
        return Value<SBInstanceCache>::downcast
            (s.updCacheEntry(getMySubsystemIndex(),topologyCache.instanceCacheIndex)).upd();
    }

    const SBTimeCache& getTimeCache(const State& s) const {
        return Value<SBTimeCache>::downcast
            (s.getCacheEntry(getMySubsystemIndex(),topologyCache.timeCacheIndex)).get();
    }
    SBTimeCache& updTimeCache(const State& s) const { //mutable
        return Value<SBTimeCache>::downcast
            (s.updCacheEntry(getMySubsystemIndex(),topologyCache.timeCacheIndex)).upd();
    }

    const SBTreePositionCache& getTreePositionCache(const State& state) const {
        return Value<SBTreePositionCache>::downcast
           (state.getCacheEntry(getMySubsystemIndex(),
                                topologyCache.treePositionCacheIndex)).get();
    }

    SBTreePositionCache& updTreePositionCache(const State& s) const { //mutable
        return Value<SBTreePositionCache>::downcast
            (s.updCacheEntry(getMySubsystemIndex(),topologyCache.treePositionCacheIndex)).upd();
    }

    const SBConstrainedPositionCache& getConstrainedPositionCache(const State& s) const {
        return Value<SBConstrainedPositionCache>::downcast
            (s.getCacheEntry(getMySubsystemIndex(),topologyCache.constrainedPositionCacheIndex)).get();
    }
    SBConstrainedPositionCache& updConstrainedPositionCache(const State& s) const { //mutable
        return Value<SBConstrainedPositionCache>::downcast
            (s.updCacheEntry(getMySubsystemIndex(),topologyCache.constrainedPositionCacheIndex)).upd();
    }

    const SBCompositeBodyInertiaCache& getCompositeBodyInertiaCache(const State& s) const {
        return Value<SBCompositeBodyInertiaCache>::downcast
            (getCacheEntry(s,topologyCache.compositeBodyInertiaCacheIndex));
    }
    SBCompositeBodyInertiaCache& updCompositeBodyInertiaCache(const State& s) const { //mutable
        return Value<SBCompositeBodyInertiaCache>::updDowncast
            (updCacheEntry(s,topologyCache.compositeBodyInertiaCacheIndex));
    }

    const SBArticulatedBodyInertiaCache& getArticulatedBodyInertiaCache(const State& s) const {
        return Value<SBArticulatedBodyInertiaCache>::downcast
            (getCacheEntry(s,topologyCache.articulatedBodyInertiaCacheIndex));
    }
    SBArticulatedBodyInertiaCache& updArticulatedBodyInertiaCache(const State& s) const { //mutable
        return Value<SBArticulatedBodyInertiaCache>::updDowncast
            (updCacheEntry(s,topologyCache.articulatedBodyInertiaCacheIndex));
    }

    const SBTreeVelocityCache& getTreeVelocityCache(const State& state) const {
        return Value<SBTreeVelocityCache>::downcast
           (state.getCacheEntry(getMySubsystemIndex(),
                                topologyCache.treeVelocityCacheIndex)).get();
    }

    SBTreeVelocityCache& updTreeVelocityCache(const State& s) const { //mutable
        return Value<SBTreeVelocityCache>::downcast
            (s.updCacheEntry(getMySubsystemIndex(),topologyCache.treeVelocityCacheIndex)).upd();
    }

    const SBConstrainedVelocityCache& getConstrainedVelocityCache(const State& s) const {
        return Value<SBConstrainedVelocityCache>::downcast
            (s.getCacheEntry(getMySubsystemIndex(),topologyCache.constrainedVelocityCacheIndex)).get();
    }
    SBConstrainedVelocityCache& updConstrainedVelocityCache(const State& s) const { //mutable
        return Value<SBConstrainedVelocityCache>::downcast
            (s.updCacheEntry(getMySubsystemIndex(),topologyCache.constrainedVelocityCacheIndex)).upd();
    }

    const SBDynamicsCache& getDynamicsCache(const State& s, bool realizingDynamics=false) const {
        const AbstractValue& cacheEntry = 
            realizingDynamics ? (const AbstractValue&)s.updCacheEntry(getMySubsystemIndex(),topologyCache.dynamicsCacheIndex)
                              : s.getCacheEntry(getMySubsystemIndex(),topologyCache.dynamicsCacheIndex);
        return Value<SBDynamicsCache>::downcast(cacheEntry).get();
    }
    SBDynamicsCache& updDynamicsCache(const State& s) const { //mutable
        return Value<SBDynamicsCache>::downcast
            (s.updCacheEntry(getMySubsystemIndex(),topologyCache.dynamicsCacheIndex)).upd();
    }

    const SBTreeAccelerationCache& getTreeAccelerationCache(const State& s) const {
        return Value<SBTreeAccelerationCache>::downcast
            (s.getCacheEntry(getMySubsystemIndex(),topologyCache.treeAccelerationCacheIndex)).get();
    }
    SBTreeAccelerationCache& updTreeAccelerationCache(const State& s) const { //mutable
        return Value<SBTreeAccelerationCache>::downcast
            (s.updCacheEntry(getMySubsystemIndex(),topologyCache.treeAccelerationCacheIndex)).upd();
    }

    const SBConstrainedAccelerationCache& getConstrainedAccelerationCache(const State& s) const {
        return Value<SBConstrainedAccelerationCache>::downcast
            (s.getCacheEntry(getMySubsystemIndex(),topologyCache.constrainedAccelerationCacheIndex)).get();
    }
    SBConstrainedAccelerationCache& updConstrainedAccelerationCache(const State& s) const { //mutable
        return Value<SBConstrainedAccelerationCache>::downcast
            (s.updCacheEntry(getMySubsystemIndex(),topologyCache.constrainedAccelerationCacheIndex)).upd();
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
            (s.getDiscreteVariable(getMySubsystemIndex(),topologyCache.topoInstanceVarsIndex)).get();
    }
    SBInstanceVars& updInstanceVars(State& s) const {
        return Value<SBInstanceVars>::downcast
            (s.updDiscreteVariable(getMySubsystemIndex(),topologyCache.topoInstanceVarsIndex)).upd();
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

    const SimbodyMatterSubsystem& getMySimbodyMatterSubsystemHandle() const 
    {   return SimbodyMatterSubsystem::downcast(getOwnerSubsystemHandle()); }
    SimbodyMatterSubsystem& updMySimbodyMatterSubsystemHandle() 
    {   return SimbodyMatterSubsystem::updDowncast(updOwnerSubsystemHandle()); }
    bool getShowDefaultGeometry() const;
    void setShowDefaultGeometry(bool show);

    void calcTreeForwardDynamicsOperator(const State&,
        const Vector&                   mobilityForces,
        const Vector_<Vec3>&            particleForces,
        const Vector_<SpatialVec>&      bodyForces,
        const Vector*                   extraMobilityForces,
        const Vector_<SpatialVec>*      extraBodyForces,
        SBTreeAccelerationCache&        tac,    // kinematics & prescribed forces into here
        Vector&                         udot,   // in/out (in for prescribed udot)
        Vector&                         qdotdot,
        Vector&                         udotErr) const;

    void calcLoopForwardDynamicsOperator(const State&, 
        const Vector&                   mobilityForces,
        const Vector_<Vec3>&            particleForces,
        const Vector_<SpatialVec>&      bodyForces,
        SBTreeAccelerationCache&        tac,    // kinematics & prescribed forces into here
        SBConstrainedAccelerationCache& cac,    // constraint forces go here
        Vector&                         udot,   // in/out (in for prescribed udot)
        Vector&                         qdotdot,
        Vector&                         multipliers,
        Vector&                         udotErr) const;

    // Given a set of forces, calculate accelerations ignoring
    // constraints, and leave the results in the state cache. 
    // Must have already called realizeDynamics().
    // We also allow some extra forces to be supplied, with the intent
    // that these will be used to deal with internal forces generated
    // by constraints; set the pointers to zero if you don't have any
    // extras to pass in.
    void realizeTreeForwardDynamics (const State& s,
        const Vector&              mobilityForces,
        const Vector_<Vec3>&       particleForces,
        const Vector_<SpatialVec>& bodyForces,
        const Vector*              extraMobilityForces,
        const Vector_<SpatialVec>* extraBodyForces) const;

    // Given a set of forces, calculate acclerations resulting from
    // those forces and enforcement of acceleration constraints, and update 
    // the state cache with the results.
    void realizeLoopForwardDynamics(const State&,
        const Vector&              mobilityForces,
        const Vector_<Vec3>&       particleForces,
        const Vector_<SpatialVec>& bodyForces) const;

    // calc ~(Tp Pq Wq^-1)_r (nfq X mp)
    void calcWeightedPqrTranspose(   
        const State&     state,
        const Vector&    Tp,    // 1/perr tols
        const Vector&    Wqinv, // 1/q weights
        Matrix&          Pqrt) const;

    // calc ~(Tp P Wu^-1)
    //       (Tv V Wu^-1)_r (nfu X (mp+mv))
    void calcWeightedPVrTranspose(
        const State&     s,
        const Vector&    Tpv,   // 1/verr tols
        const Vector&    Wuinv, // 1/u weights
        Matrix&          PVrt) const;

    const Array_<QIndex>& getFreeQIndex(const State& state) const;
    const Array_<QIndex>& getPresQIndex(const State& state) const;
    const Array_<QIndex>& getZeroQIndex(const State& state) const;

    const Array_<UIndex>& getFreeUIndex(const State& state) const;
    const Array_<UIndex>& getPresUIndex(const State& state) const;
    const Array_<UIndex>& getZeroUIndex(const State& state) const;

    const Array_<UIndex>& getFreeUDotIndex(const State& state) const;
    const Array_<UIndex>& getKnownUDotIndex(const State& state) const;

    // Output must already be sized for number of free q's nfq.
    // Input must be size nq.
    void packFreeQ(const State& s, const Vector& allQ,
                   Vector& packedFreeQ) const;

    // For efficiency, you must provide an output array of the right size nq.
    // This method *will not* touch the prescribed slots in the output so
    // if you want them zero make sure you do it yourself.
    void unpackFreeQ(const State& s, const Vector& packedFreeQ,
                     Vector& unpackedFreeQ) const;

    // Given a q-like array with nq entries, write zeroes onto the entries
    // corresponding to known (prescribed) q's. The result looks like
    // a properly-zeroed unpackedFreeQ.
    void zeroKnownQ(const State& s, Vector& qlike) const;

    // Output must already be sized for number of free u's nfu.
    // Input must be size nu.
    void packFreeU(const State& s, const Vector& allU,
                   Vector& packedFreeU) const;

    // For efficiency, you must provide an output array of the right size nu.
    // This method *will not* touch the prescribed slots in the output so
    // if you want them zero make sure you do it yourself.
    void unpackFreeU(const State& s, const Vector& packedFreeU,
                     Vector& unpackedFreeU) const;

    // Given a u-like array with nu entries, write zeroes onto the entries
    // corresponding to known (prescribed) u's. The result looks like
    // a properly-zeroed unpackedFreeU.
    void zeroKnownU(const State& s, Vector& ulike) const;

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
    Array_<MobilizedBody*,MobilizedBodyIndex>               mobilizedBodies;
    // Constraints are treated similarly.
    Array_<Constraint*,ConstraintIndex>                     constraints;

    Array_<UnilateralContact*,UnilateralContactIndex>       uniContacts;
    Array_<StateLimitedFriction*,StateLimitedFrictionIndex> stateLtdFriction;

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
    Array_<RBNodePtrList>      rbNodeLevels;
    // Map nodeNum (a.k.a. MobilizedBodyIndex) to (level,offset).
    Array_<RigidBodyNodeIndex,MobilizedBodyIndex> nodeNum2NodeMap;

        // Constraints

    // Here we sort the above constraints by branch (ancestor's base body), then by
    // level within that branch. That is, each constraint is addressed
    // by three indices [branch][levelOfAncestor][offset]
    // where offset is an arbitrary unique integer assigned to all
    // the constraints on the same branch with the same level of ancestor.
    Array_< Array_< Array_<ConstraintIndex> > > branches;

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
    Array_<CoupledConstraintSet> positionCoupledConstraints;     // for P
    Array_<CoupledConstraintSet> velocityCoupledConstraints;     // for PV
    Array_<CoupledConstraintSet> accelerationCoupledConstraints; // for G=PVA

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
    Array_<CoupledConstraintSet> dynamicallyCoupledConstraints;


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
    CacheEntryIndex topologyCacheIndex; // topologyCache is copied here in the State
    
    // Specifies whether default decorative geometry should be shown.
    bool showDefaultGeometry;
};

std::ostream& operator<<(std::ostream&, const SimbodyMatterSubsystemRep&);


#endif // SimTK_SIMBODY_MATTER_SUBSYSTEM_REP_H_
