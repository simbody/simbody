#ifndef SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_H_
#define SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-9 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *    Charles Schwieters (NIH): wrote the public domain IVM code from which   *
 *                              this was derived.                             *
 *    Peter Eastman: wrote the Euler Angle<->Quaternion conversion            *
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


/**@file
 * This file contains just the templatized class RigidBodyNodeSpec<dof> which
 * is used as the base class for most of the RigidBodyNodes (that is, the
 * implementations of Mobilizers). The only exceptions are nodes whose 
 * mobilizers provide no degrees of freedom -- Ground and Weld.
 *
 * This file contains all the multibody mechanics method declarations that 
 * involve a single body and its mobilizer (inboard joint), that is, one node 
 * in the multibody tree. These methods constitute the inner loops of the 
 * multibody calculations, and much suffering is undergone here to make them 
 * run fast. In particular most calculations are templatized by the number of 
 * mobilities, so that compile-time sizes are known for everything.
 *
 * Most methods here expect to be called in a particular order during traversal 
 * of the tree -- either base to tip or tip to base.
 */

#include "SimbodyMatterSubsystemRep.h"
#include "RigidBodyNode.h"

/**
 * This still-abstract class is a skeleton implementation of a built-in 
 * mobilizer, with the number of mobilities (dofs, u's) specified as a template 
 * parameter. That way all the code that is common except for the dimensionality
 * of the mobilizer can be written once, and the compiler generates specific 
 * implementations for each of the six possible dimensionalities (1-6 
 * mobilities). Each implementation works only on fixed size Vec<> and Mat<> 
 * types, so can use very high speed inline operators.
 */
template<int dof, bool noR_FM=false, bool noX_MB=false, bool noR_PF=false>
class RigidBodyNodeSpec : public RigidBodyNode {
public:

RigidBodyNodeSpec(const MassProperties& mProps_B,
                  const Transform&      X_PF,
                  const Transform&      X_BM,
                  UIndex&               nextUSlot,
                  USquaredIndex&        nextUSqSlot,
                  QIndex&               nextQSlot,
                  QDotHandling          qdotHandling,
                  QuaternionUse         quaternionUse,
                  bool                  isReversed)
:   RigidBodyNode(mProps_B, X_PF, X_BM, 
                  qdotHandling, quaternionUse, isReversed)
{
    // don't call any virtual methods in here!
    uIndex   = nextUSlot;
    uSqIndex = nextUSqSlot;
    qIndex   = nextQSlot;
}

void updateSlots(UIndex& nextUSlot, USquaredIndex& nextUSqSlot, QIndex& nextQSlot) {
    // OK to call virtual method here.
    nextUSlot   += getDOF();
    nextUSqSlot += getDOF()*getDOF();
    nextQSlot   += getMaxNQ();
}

virtual ~RigidBodyNodeSpec() {}

// This is the type of the joint transition matrix H, but our definition
// of H is transposed from Jain's or Schwieters'. That is, what we're 
// calling H here is what Jain calls H* and Schwieters calls H^T. So
// our H matrix is 6 x nu, or more usefully it is 2 rows of Vec3's.
// The first row define how u's contribute to angular velocities; the
// second defines how u's contribute to linear velocities.
typedef Mat<2,dof,Vec3> HType;

// Provide default implementations for setQToFitTransformImpl() and setQToFitVelocityImpl() 
// which are implemented using the rotational and translational quantity routines. These assume
// that the rotational and translational coordinates are independent, with rotation handled
// first and then left alone. If a mobilizer type needs to deal with rotation and
// translation simultaneously, it should provide a specific implementation for these two routines.
// *Each* mobilizer must implement setQToFit{Rotation,Translation} and 
// setUToFit{AngularVelocity,LinearVelocity}; there are no defaults.

virtual void setQToFitTransformImpl(const SBStateDigest& sbs, const Transform& X_FM, Vector& q) const {
    setQToFitRotationImpl   (sbs,X_FM.R(),q);
    setQToFitTranslationImpl(sbs,X_FM.p(),q);
}

virtual void setUToFitVelocityImpl(const SBStateDigest& sbs, const Vector& q, const SpatialVec& V_FM, Vector& u) const {
    setUToFitAngularVelocityImpl(sbs,q,V_FM[0],u);
    setUToFitLinearVelocityImpl (sbs,q,V_FM[1],u);
}

// The following routines calculate joint-specific position kinematic
// quantities. They may assume that *all* position kinematics (not just
// joint-specific) has been done for the parent, and that the position
// state variables q are available. Each routine may also assume that the
// previous routines have been called, in the order below.
// The routines are structured as operators -- they use the State but
// do not change anything in it, including the cache. Instead they are
// passed arguments to write their results into. In practice, these
// arguments will typically be in the State cache (see below).


// This mandatory routine calculates the joint transition matrix H_FM, giving
// the change of velocity induced by the generalized speeds u for this 
// mobilizer, expressed in the mobilizer's inboard "fixed" frame F (attached to
// this body's parent). 
// This method constitutes the *definition* of the generalized speeds for
// a particular joint.
// This routine can depend on X_FM having already
// been calculated and available in the PositionCache but must NOT depend
// on any quantities involving Ground or other bodies.
// Note: this calculates the H matrix *as defined*; we might reverse
// the frames in use. So we call this H_F0M0 while the possibly reversed
// version is H_FM. Be sure to use X_F0M0 in your calculation for H_F0M0
// if you need access to the cross-mobilizer transform.
// CAUTION: our definition of H is transposed from Jain and Schwieters.
virtual void calcAcrossJointVelocityJacobian(
    const SBStateDigest& sbs,
    HType&               H_F0M0) const=0;

// This mandatory routine calculates the time derivative taken in F of
// the matrix H_FM (call it HDot_FM). This is zero if the generalized
// speeds are all defined in terms of the F frame, which is often the case.
// This routine can depend on X_FM and H_FM being available already in
// the state PositionCache, and V_FM being already in the state VelocityCache.
// However, it must NOT depend on any quantities involving Ground or other bodies.
// Note: this calculates the HDot matrix *as defined*; we might reverse
// the frames in use. So we call this HDot_F0M0 while the possibly reversed
// version is HDot_FM. Be sure to use X_F0M0 and V_F0M0 in your calculation
// of HDot_F0M0 if you need access to the cross-mobilizer transform and/or velocity.
// CAUTION: our definition of H is transposed from Jain and Schwieters.
virtual void calcAcrossJointVelocityJacobianDot(
    const SBStateDigest& sbs,
    HType&               HDot_F0M0) const = 0;

// We allow a mobilizer to be defined in reverse when that is more
// convenient. That is, the H matrix can be defined by giving H_F0M0=H_MF and 
// HDot_F0M0=HDot_MF instead of H_FM and HDot_FM. In that case the 
// following methods are called instead of the above; the default 
// implementation postprocesses the output from the above methods, but a 
// mobilizer can override it if it knows how to get the job done faster.
virtual void calcReverseMobilizerH_FM(
    const SBStateDigest& sbs,
    HType&               H_FM) const;

virtual void calcReverseMobilizerHDot_FM(
    const SBStateDigest& sbs,
    HType&               HDot_FM) const;

// This routine is NOT joint specific, but cannot be called until the across-joint
// transform X_FM has been calculated and is available in the State cache.
void calcBodyTransforms(
    const SBTreePositionCache&  pc, 
    Transform&                  X_PB, 
    Transform&                  X_GB) const 
{
    const Transform& X_MB = getX_MB();   // fixed
    const Transform& X_PF = getX_PF();   // fixed
    const Transform& X_FM = getX_FM(pc); // just calculated
    const Transform& X_GP = getX_GP(pc); // already calculated

    const Transform X_FB = (noX_MB ? X_FM : X_FM*X_MB); // 63 flops
    X_PB = (noR_PF ? Transform(X_FB.R(), X_PF.p()+X_FB.p()) : X_PF*X_FB); // 63 flops
    X_GB = X_GP * X_PB;        //  63 flops
}

// Same for all mobilizers. The return matrix here is precisely the 
// one used by Jain and Schwieters, but transposed.
void calcParentToChildVelocityJacobianInGround(
    const SBModelVars&          mv,
    const SBTreePositionCache&  pc, 
    HType&                      H_PB_G) const;

// Same for all mobilizers. This is the time derivative of 
// the matrix H_PB_G above, with the derivative taken in the 
// Ground frame.
void calcParentToChildVelocityJacobianInGroundDot(
    const SBModelVars&          mv,
    const SBTreePositionCache&  pc, 
    const SBTreeVelocityCache&  vc, 
    HType&                      HDot_PB_G) const;

void realizeModel(SBStateDigest& sbs) const 
{
}

void realizeInstance(const SBStateDigest& sbs) const
{
}

// Set a new configuration and calculate the consequent kinematics.
// Must call base-to-tip.
void realizePosition(const SBStateDigest& sbs) const 
{
    const SBModelVars&      mv   = sbs.getModelVars();
    const SBModelCache&     mc   = sbs.getModelCache();
    const SBInstanceCache&  ic   = sbs.getInstanceCache();
    const Vector&           allQ = sbs.getQ();
    SBTreePositionCache&    pc   = sbs.updTreePositionCache();
    Vector&                 allQErr = sbs.updQErr();

    // Mobilizer specific.
    
    const SBModelCache::PerMobilizedBodyModelInfo& mbInfo = getModelInfo(mc);

    // First perform precalculations on these new q's, such as stashing away
    // sines and cosines of angles. We'll put all the results into the State's
    // position cache for later use. If there are local constraints on the q's
    // (currently that means quaternion normalization constraints), then we
    // calculate those here and put the result in the appropriate slot of qerr.

    const int nq=mbInfo.nQInUse, nqpool=mbInfo.nQPoolInUse,
              nqerr=(mbInfo.hasQuaternionInUse ? 1 : 0);
    const Real* q0    = nq     ? &allQ[mbInfo.firstQIndex]                : 0;
    Real*       qpool0= nqpool ? &pc.mobilizerQCache[mbInfo.startInQPool] : 0;
    Real*       qerr0 = nqerr  ? &allQErr[ic.firstQuaternionQErrSlot
                                          + mbInfo.quaternionPoolIndex]     
                               : 0;
    performQPrecalculations(sbs, q0, nq, qpool0, nqpool, qerr0, nqerr);

    // Now that we've done the necessary precalculations, calculate the cross-
    // mobilizer transform X_FM without recalculating anything. For reversed 
    // mobilizers we have to do our own reversal here because mobilizers don't
    // handle that themselves.

    if (isReversed()) {
        Transform X_MF;
        calcX_FM(sbs, q0, nq, qpool0, nqpool, X_MF);
        updX_FM(pc) = ~X_MF;
    } else 
        calcX_FM(sbs, q0, nq, qpool0, nqpool, updX_FM(pc));

    // With X_FM in the cache, and X_GP for the parent already calculated (we're doing
    // an outward pass), we can calculate X_PB and X_GB now.
    calcBodyTransforms(pc, updX_PB(pc), updX_GB(pc));

    // Here we do allow the mobilizer to calculate the reversed H matrix, but the
    // default implementation of the reversed method just calls the forward method
    // and then reverses it. For some mobilizers that is unreasonably expensive.
    // REMINDER: our H matrix definition is transposed from Jain and Schwieters.
    if (isReversed()) calcReverseMobilizerH_FM       (sbs, updH_FM(pc));
    else              calcAcrossJointVelocityJacobian(sbs, updH_FM(pc));

    // Here we're using the cross-mobilizer hinge matrix H_FM that we just 
    // calculated to compute H(==H_PB_G), the equivalent hinge matrix between 
    // the parent's body frame P and child's body frame B, but expressed in
    // Ground. (F is fixed on P and M is fixed on B.)
    calcParentToChildVelocityJacobianInGround(mv,pc, updH(pc));

    // Mobilizer independent.
    calcJointIndependentKinematicsPos(pc);
}

// Set new velocities for the current configuration, and calculate
// all the velocity-dependent terms. Must call base-to-tip.
// This routine may assume that *all* position 
// kinematics (not just joint-specific) has been done for this node,
// that all velocity kinematics has been done for the parent, and
// that the velocity state variables (u) are available. The
// quantities that must be computed are:
//   V_FM   relative velocity of B's M frame in P's F frame, 
//             expressed in F (note: this is also V_PM_F since
//             F is fixed on P). (This is H_FM*u.)
//   V_PB_G  relative velocity of B in P (==H*u), expr. in G
//   HDot_FM time derivative of hinge matrix H_FM
//   HDot    time derivative of H (==H_PB_G) hinge matrix relating body frames
//   VD_PB_G acceleration remainder term HDot*u, expr. in G
// The code is the same for all joints, although parametrized by ndof.
void realizeVelocity(const SBStateDigest& sbs) const 
{
    const SBModelVars&          mv = sbs.getModelVars();
    const SBTreePositionCache&  pc = sbs.getTreePositionCache();
    SBTreeVelocityCache&        vc = sbs.updTreeVelocityCache();
    const Vector&               allU = sbs.getU();
    const Vec<dof>&             u = fromU(allU);

    // Mobilizer specific.
    calcQDot(sbs, &allU[uIndex], &sbs.updQDot()[qIndex]);

    updV_FM(vc)    = getH_FM(pc) * u;   // 6*dof flops
    updV_PB_G(vc)  = getH(pc)    * u;   // 6*dof flops

    // REMINDER: our H matrix definition is transposed from Jain and Schwieters.
    if (isReversed()) calcReverseMobilizerHDot_FM       (sbs, updHDot_FM(vc));
    else              calcAcrossJointVelocityJacobianDot(sbs, updHDot_FM(vc));

    // Here we're using the cross-mobilizer hinge matrix derivative HDot_FM 
    // that we just calculated to compute HDot(==HDot_PB_G), the derivative
    // of H (==H_PB_G) between the parent's body frame P and child's body 
    // frame B. The derivative is taken in the Ground frame and the result
    // is expressed in Ground. (F is fixed on P and M is fixed on B.)
    calcParentToChildVelocityJacobianInGroundDot(mv,pc,vc, updHDot(vc));
    updVD_PB_G(vc) = getHDot(vc) * u;   // 6*dof flops

    // Mobilizer independent.
    calcJointIndependentKinematicsVel(pc,vc);
}

void realizeDynamics(const SBArticulatedBodyInertiaCache&   abc,
                     const SBStateDigest&                   sbs) const 
{
    const SBTreePositionCache&  pc = sbs.getTreePositionCache();
    const SBTreeVelocityCache&  vc = sbs.getTreeVelocityCache();
    SBDynamicsCache&            dc = sbs.updDynamicsCache();

    // Mobilizer-specific.
    // None.

    // Mobilizer independent.
    calcJointIndependentDynamicsVel(pc,abc,vc,dc);
}

void realizeAcceleration(const SBStateDigest& sbs) const
{
}

void realizeReport(const SBStateDigest& sbs) const
{
}

// Use base class implementation of calcCompositeBodyInertiasInward() since
// that is independent of mobilities.

// This is a dynamics-stage calculation and must be called tip-to-base (inward).
void realizeArticulatedBodyInertiasInward(
    const SBInstanceCache&          ic,
    const SBTreePositionCache&      pc,
    SBArticulatedBodyInertiaCache&  abc) const;

// calcJointIndependentDynamicsVel() must be called after ArticulatedBodyInertias.

// This dynamics-stage calculation is needed for handling constraints. It
// must be called base-to-tip (outward);
void realizeYOutward(
    const SBInstanceCache&                  ic,
    const SBTreePositionCache&              pc,
    const SBArticulatedBodyInertiaCache&    abc,
    SBDynamicsCache&                        dc) const;

void realizeZ(
    const SBTreePositionCache&              pc,
    const SBArticulatedBodyInertiaCache&    abc,
    const SBTreeVelocityCache&              vc,
    const SBDynamicsCache&                  dc,
    SBTreeAccelerationCache&                ac,
    const Real*                             mobilityForces,
    const SpatialVec*                       bodyForces) const;

void realizeAccel(
    const SBTreePositionCache&              pc,
    const SBArticulatedBodyInertiaCache&    abc,
    const SBTreeVelocityCache&              vc,
    const SBDynamicsCache&                  dc,
    SBTreeAccelerationCache&                ac,
    Real*                                   udot) const;

// These routines give each node a chance to set appropriate defaults in a piece
// of the state corresponding to a particular stage. Default implementations here
// assume non-ball joint; override if necessary.
virtual void setMobilizerDefaultModelValues   (const SBTopologyCache&, SBModelVars&)  const {}
virtual void setMobilizerDefaultInstanceValues(const SBModelVars&, SBInstanceVars&) const {}
virtual void setMobilizerDefaultTimeValues    (const SBModelVars&, SBTimeVars&)    const {}

virtual void setMobilizerDefaultPositionValues(const SBModelVars& s, Vector& q) const 
{
    toQ(q) = 0.;
}
virtual void setMobilizerDefaultVelocityValues(const SBModelVars&, Vector& u) const 
{
    toU(u) = 0.;
}
virtual void setMobilizerDefaultDynamicsValues(const SBModelVars&, SBDynamicsVars&) const {}
virtual void setMobilizerDefaultAccelerationValues(const SBModelVars&, 
                                                   SBDynamicsVars& v) const {}

int          getDOF()            const {return dof;}
virtual int  getMaxNQ()          const {
    assert(quaternionUse == QuaternionIsNeverUsed);
    return dof; // maxNQ can be larger than dof if there's a quaternion
}

// Default method must be overridden if there is a chance that NQ
// might not be the same as NU.
virtual int  getNQInUse(const SBModelVars&) const {
    assert(quaternionUse == QuaternionIsNeverUsed); 
    return dof; // DOF <= NQ <= maxNQ
}

// Currently NU is always just the Mobilizer's DOFs (the template argument)
// Later we may want to offer modeling options to lock joints, or perhaps
// break them.
virtual int  getNUInUse(const SBModelVars&) const {
    return dof;
}

// Default method must be overridden if this mobilizer might use a
// quaternion under some conditions.
virtual bool isUsingQuaternion(const SBStateDigest&, MobilizerQIndex& startOfQuaternion) const {
    assert(quaternionUse == QuaternionIsNeverUsed);
    startOfQuaternion.invalidate();
    return false;
}

// You must override the next few methods if qdot might not be the same
// as u under *any* circumstance.

// This method should calculate qdot=N*u, where N=N(q) is the kinematic
// coupling matrix. State digest should be at Stage::Position.
virtual void calcQDot(const SBStateDigest&, 
                                     const Real* u, Real* qdot) const 
{
    assert(qdotHandling == QDotIsAlwaysTheSameAsU);
    Vec<dof>::updAs(qdot) = Vec<dof>::getAs(u); // default says qdot=u
}

// This method should calculate u=N^-1*qdot, where N=N(q).
// State digest should be at Stage::Position.
virtual void calcLocalUFromLocalQDot(const SBStateDigest&, 
                                     const Real* qdot, Real* u) const 
{
    assert(qdotHandling == QDotIsAlwaysTheSameAsU);
    Vec<dof>::updAs(u) = Vec<dof>::getAs(qdot); // default says u=qdot
}

// This method should calculate qdotdot=N*udot + NDot*u, where N=N(q),
// NDot=NDot(q,u). State digest should be at Stage::Velocity.
virtual void calcQDotDot(const SBStateDigest&, 
                                           const Real* udot, Real* qdotdot) const
{
    assert(qdotHandling == QDotIsAlwaysTheSameAsU);
    Vec<dof>::updAs(qdotdot) = Vec<dof>::getAs(udot); // default: qdotdot=udot
}

// This method should calculate udot=N^-1*(qdotdot-NDot*u), where N=N(q),
// NDot=NDot(q,u). State digest should be at Stage::Velocity.
virtual void calcLocalUDotFromLocalQDotDot(const SBStateDigest&, 
                                           const Real* qdotdot, Real* udot) const
{
    assert(qdotHandling == QDotIsAlwaysTheSameAsU);
    Vec<dof>::updAs(udot) = Vec<dof>::getAs(qdotdot); // default: udot=qdotdot
}




// State digest should be at Stage::Position for calculating N (the matrix that
// maps generalized speeds u to coordinate derivatives qdot, qdot=Nu). The 
// default implementation assumes nq==nu and the nqXnu block of N corresponding 
// to this mobilizer is identity. Then either operation (regardless of side) 
// just copies nu numbers from in to out.
//
// THIS MUST BE OVERRIDDEN by any mobilizer for which nq != nu, or qdot != u.
virtual void multiplyByN(const SBStateDigest&, bool matrixOnRight,  
                         const Real* in, Real* out) const
{
    assert(qdotHandling == QDotIsAlwaysTheSameAsU);
    Vec<dof>::updAs(out) = Vec<dof>::getAs(in);
}

virtual void multiplyByNInv(const SBStateDigest&, bool matrixOnRight, 
                            const Real* in, Real* out) const
{
    assert(qdotHandling == QDotIsAlwaysTheSameAsU);
    Vec<dof>::updAs(out) = Vec<dof>::getAs(in);
}

// State digest should be at Stage::Velocity for calculating NDot (the matrix 
// that is used in mapping generalized accelerations udot to coordinate 2nd 
// derivatives qdotdot, qdotdot=N udot + NDot u. The default implementation 
// assumes nq==nu and the nqXnu block of N corresponding to this mobilizer is 
// identity, so NDot is an nuXnu block of zeroes. Then either operation 
// (regardless of side) just copies nu zeros to out.
//
// THIS MUST BE OVERRIDDEN by any mobilizer for which nq != nu, or qdot != u.
virtual void multiplyByNDot(const SBStateDigest&, bool matrixOnRight, 
                            const Real* in, Real* out) const
{
    assert(qdotHandling == QDotIsAlwaysTheSameAsU);
    Vec<dof>::updAs(out) = 0;
}


// No default implementations here for:
// calcMobilizerTransformFromQ
// calcMobilizerVelocityFromU
// calcMobilizerAccelerationFromUDot
// calcParentToChildTransformFromQ
// calcParentToChildVelocityFromU
// calcParentToChildAccelerationFromUDot


virtual void setVelFromSVel(
    const SBTreePositionCache&  pc, 
    const SBTreeVelocityCache&  vc,
    const SpatialVec&           sVel, 
    Vector&                     u) const;

// Return true if any change is made to the q vector.
virtual bool enforceQuaternionConstraints(
    const SBStateDigest&    sbs,
    Vector&                 q,
    Vector&                 qErrest) const 
{
    assert(quaternionUse == QuaternionIsNeverUsed);
    return false;
}

void convertToEulerAngles(const Vector& inputQ, Vector& outputQ) const {
    // The default implementation just copies Q.  Subclasses may override this.
    assert(quaternionUse == QuaternionIsNeverUsed);
    toQ(outputQ) = fromQ(inputQ);
}

void convertToQuaternions(const Vector& inputQ, Vector& outputQ) const {
    // The default implementation just copies Q.  Subclasses may override this.
    assert(quaternionUse == QuaternionIsNeverUsed);
    toQ(outputQ) = fromQ(inputQ);
}

// Get a column of H_PB_G, which is what Jain calls H* and Schwieters calls H^T.
const SpatialVec& getHCol(const SBTreePositionCache& pc, int j) const {
    return getH(pc)(j);
}

// Get a column of H_FM the local cross-mobilizer hinge matrix expressed in the
// parent (inboard) mobilizer frame F.
const SpatialVec& getH_FMCol(const SBTreePositionCache& pc, int j) const {
    return getH_FM(pc)(j);
}

// Access to body-oriented state and cache entries is the same for all nodes,
// and joint oriented access is almost the same but parametrized by dof. There is a special
// case for quaternions because they use an extra state variable, and although we don't
// have to we make special scalar routines available for 1-dof joints. Note that all State access
// routines are inline, not virtual, so the cost is just an indirection and an index.
//
// TODO: these inner-loop methods probably shouldn't be indexing a Vector, which requires
// several indirections. Instead, the top-level caller should find the Real* data contained in the
// Vector and then pass that to the RigidBodyNode routines which will call these ones.

// General joint-dependent select-my-goodies-from-the-pool routines.
const Vec<dof>&     fromQ  (const Real* q)   const {return Vec<dof>::getAs(&q[qIndex]);}
Vec<dof>&           toQ    (      Real* q)   const {return Vec<dof>::updAs(&q[qIndex]);}
const Vec<dof>&     fromU  (const Real* u)   const {return Vec<dof>::getAs(&u[uIndex]);}
Vec<dof>&           toU    (      Real* u)   const {return Vec<dof>::updAs(&u[uIndex]);}
const Mat<dof,dof>& fromUSq(const Real* uSq) const {return Mat<dof,dof>::getAs(&uSq[uSqIndex]);}
Mat<dof,dof>&       toUSq  (      Real* uSq) const {return Mat<dof,dof>::updAs(&uSq[uSqIndex]);}

// Same, but specialized for the common case where dof=1 so everything is scalar.
const Real& from1Q  (const Real* q)   const {return q[qIndex];}
Real&       to1Q    (      Real* q)   const {return q[qIndex];}
const Real& from1U  (const Real* u)   const {return u[uIndex];}
Real&       to1U    (      Real* u)   const {return u[uIndex];}
const Real& from1USq(const Real* uSq) const {return uSq[uSqIndex];}
Real&       to1USq  (      Real* uSq) const {return uSq[uSqIndex];}

// Same, specialized for quaternions. We're assuming that the quaternion is stored
// in the *first* four coordinates, if there are more than four altogether.
const Vec4& fromQuat(const Real* q) const {return Vec4::getAs(&q[qIndex]);}
Vec4&       toQuat  (      Real* q) const {return Vec4::updAs(&q[qIndex]);}

// Extract a Vec3 from a Q-like or U-like object, beginning at an offset from the qIndex or uIndex.
const Vec3& fromQVec3(const Real* q, int offs) const {return Vec3::getAs(&q[qIndex+offs]);}
Vec3&       toQVec3  (      Real* q, int offs) const {return Vec3::updAs(&q[qIndex+offs]);}
const Vec3& fromUVec3(const Real* u, int offs) const {return Vec3::getAs(&u[uIndex+offs]);}
Vec3&       toUVec3  (      Real* u, int offs) const {return Vec3::updAs(&u[uIndex+offs]);}

// These provide an identical interface for when q,u are given as Vectors rather than Reals.

const Vec<dof>&     fromQ  (const Vector& q)   const {return fromQ(&q[0]);} // convert to array of Real
Vec<dof>&           toQ    (      Vector& q)   const {return toQ  (&q[0]);}
const Vec<dof>&     fromU  (const Vector& u)   const {return fromU(&u[0]);}
Vec<dof>&           toU    (      Vector& u)   const {return toU  (&u[0]);}
const Mat<dof,dof>& fromUSq(const Vector& uSq) const {return fromUSq(&uSq[0]);}
Mat<dof,dof>&       toUSq  (      Vector& uSq) const {return toUSq  (&uSq[0]);}

const Real& from1Q  (const Vector& q)   const {return from1Q  (&q[0]);}
Real&       to1Q    (      Vector& q)   const {return to1Q    (&q[0]);}
const Real& from1U  (const Vector& u)   const {return from1U  (&u[0]);}
Real&       to1U    (      Vector& u)   const {return to1U    (&u[0]);}
const Real& from1USq(const Vector& uSq) const {return from1USq(&uSq[0]);}
Real&       to1USq  (      Vector& uSq) const {return to1USq  (&uSq[0]);}

// Same, specialized for quaternions. We're assuming that the quaternion comes first in the coordinates.
const Vec4& fromQuat(const Vector& q) const {return fromQuat(&q[0]);}
Vec4&       toQuat  (      Vector& q) const {return toQuat  (&q[0]);}

// Extract a Vec3 from a Q-like or U-like object, beginning at an offset from the qIndex or uIndex.
const Vec3& fromQVec3(const Vector& q, int offs) const {return fromQVec3(&q[0], offs);}
Vec3&       toQVec3  (      Vector& q, int offs) const {return toQVec3  (&q[0], offs);}
const Vec3& fromUVec3(const Vector& u, int offs) const {return fromUVec3(&u[0], offs);}
Vec3&       toUVec3  (      Vector& u, int offs) const {return toUVec3  (&u[0], offs);}

// Applications of the above extraction routines to particular interesting items in the State. Note
// that you can't use these for quaternions since they extract "dof" items.

// Cache entries (cache is mutable in a const State)

    // Position


// CAUTION: our H definition is transposed from Jain and Schwieters.
const HType& getH_FM(const SBTreePositionCache& pc) const
  { return HType::getAs(&pc.storageForH_FM[2*uIndex]); }
HType&       updH_FM(SBTreePositionCache& pc) const
  { return HType::updAs(&pc.storageForH_FM[2*uIndex]); }

// "H" here should really be H_PB_G, that is, cross joint transition
// matrix relating parent and body frames, but expressed in Ground.
// CAUTION: our H definition is transposed from Jain and Schwieters.
const HType& getH(const SBTreePositionCache& pc) const
  { return HType::getAs(&pc.storageForH[2*uIndex]); }
HType&       updH(SBTreePositionCache& pc) const
  { return HType::updAs(&pc.storageForH[2*uIndex]); }

    // Velocity

// CAUTION: our H definition is transposed from Jain and Schwieters.
const HType& getHDot_FM(const SBTreeVelocityCache& vc) const
  { return HType::getAs(&vc.storageForHDot_FM[2*uIndex]); }
HType&       updHDot_FM(SBTreeVelocityCache& vc) const
  { return HType::updAs(&vc.storageForHDot_FM[2*uIndex]); }

// CAUTION: our H definition is transposed from Jain and Schwieters.
const HType& getHDot(const SBTreeVelocityCache& vc) const
  { return HType::getAs(&vc.storageForHDot[2*uIndex]); }
HType&       updHDot(SBTreeVelocityCache& vc) const
  { return HType::updAs(&vc.storageForHDot[2*uIndex]); }

    // Dynamics

// These are calculated with articulated body inertias.
const Mat<dof,dof>& getD(const SBArticulatedBodyInertiaCache& abc) const 
{return fromUSq(abc.storageForD);}
Mat<dof,dof>&       updD(SBArticulatedBodyInertiaCache&       abc) const 
{return toUSq  (abc.storageForD);}

const Mat<dof,dof>& getDI(const SBArticulatedBodyInertiaCache& abc) const 
{return fromUSq(abc.storageForDI);}
Mat<dof,dof>&       updDI(SBArticulatedBodyInertiaCache&       abc) const 
{return toUSq  (abc.storageForDI);}

const Mat<2,dof,Vec3>& getG(const SBArticulatedBodyInertiaCache& abc) const
  { return Mat<2,dof,Vec3>::getAs(&abc.storageForG[2*uIndex]); }
Mat<2,dof,Vec3>&       updG(SBArticulatedBodyInertiaCache&       abc) const
  { return Mat<2,dof,Vec3>::updAs(&abc.storageForG[2*uIndex]); }

const Vec<dof>& getTau(const SBInstanceCache& ic, const Real* tau) const {
    const PresForcePoolIndex tauIx = ic.getMobodInstanceInfo(nodeNum).firstPresForce;
    assert(tauIx.isValid());
    return Vec<dof>::getAs(&tau[tauIx]);
}
Vec<dof>& updTau(const SBInstanceCache& ic, Real* tau) const {
    const PresForcePoolIndex tauIx = ic.getMobodInstanceInfo(nodeNum).firstPresForce;
    assert(tauIx.isValid());
    return Vec<dof>::updAs(&tau[tauIx]);
}

    // Acceleration

const Vec<dof>&   getEpsilon (const SBTreeAccelerationCache& ac) const {return fromU (ac.epsilon);}
Vec<dof>&         updEpsilon (SBTreeAccelerationCache&       ac) const {return toU   (ac.epsilon);}
const Real&       get1Epsilon(const SBTreeAccelerationCache& ac) const {return from1U(ac.epsilon);}
Real&             upd1Epsilon(SBTreeAccelerationCache&       ac) const {return to1U  (ac.epsilon);}


void multiplyBySystemJacobian(
    const SBTreePositionCache&  pc,
    const Real*                 v,
    SpatialVec*                 Jv) const;

void multiplyBySystemJacobianTranspose(
    const SBTreePositionCache&  pc, 
    SpatialVec*                 zTmp,
    const SpatialVec*           X, 
    Real*                       JtX) const;

void calcEquivalentJointForces(
    const SBTreePositionCache&  pc,
    const SBDynamicsCache&      dc,
    const SpatialVec*           bodyForces,
    SpatialVec*                 allZ,
    Real*                       jointForces) const;

void calcUDotPass1Inward(
    const SBInstanceCache&      ic,
    const SBTreePositionCache&  pc,
    const SBArticulatedBodyInertiaCache&,
    const SBDynamicsCache&      dc,
    const Real*                 jointForces,
    const SpatialVec*           bodyForces,
    const Real*                 allUDot,
    SpatialVec*                 allZ,
    SpatialVec*                 allGepsilon,
    Real*                       allEpsilon) const;

void calcUDotPass2Outward(
    const SBInstanceCache&      ic,
    const SBTreePositionCache&  pc,
    const SBArticulatedBodyInertiaCache&,
    const SBTreeVelocityCache&  vc,
    const SBDynamicsCache&      dc,
    const Real*                 epsilonTmp,
    SpatialVec*                 allA_GB,
    Real*                       allUDot,
    Real*                       allTau) const;

void calcMInverseFPass1Inward(
    const SBInstanceCache&      ic,
    const SBTreePositionCache&  pc,
    const SBArticulatedBodyInertiaCache&,
    const SBDynamicsCache&      dc,
    const Real*                 f,
    SpatialVec*                 allZ,
    SpatialVec*                 allGepsilon,
    Real*                       allEpsilon) const;

void calcMInverseFPass2Outward(
    const SBInstanceCache&      ic,
    const SBTreePositionCache&  pc,
    const SBArticulatedBodyInertiaCache&,
    const SBDynamicsCache&      dc,
    const Real*                 epsilonTmp,
    SpatialVec*                 allA_GB,
    Real*                       allUDot) const;

void calcInverseDynamicsPass1Outward(
    const SBTreePositionCache&  pc,
    const SBTreeVelocityCache&  vc,
    const Real*                 allUDot,
    SpatialVec*                 allA_GB) const;

void calcInverseDynamicsPass2Inward(
    const SBTreePositionCache&  pc,
    const SBTreeVelocityCache&  vc,
    const SpatialVec*           allA_GB,
    const Real*                 jointForces,
    const SpatialVec*           bodyForces,
    SpatialVec*                 allFTmp,
    Real*                       allTau) const; 

void calcMVPass1Outward(
    const SBTreePositionCache&  pc,
    const Real*                 allUDot,
    SpatialVec*                 allA_GB) const;
void calcMVPass2Inward(
    const SBTreePositionCache&  pc,
    const SpatialVec*           allA_GB,
    SpatialVec*                 allFTmp,
    Real*                       allTau) const;

};

#endif // SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_H_
