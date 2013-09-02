#ifndef SimTK_SIMBODY_RIGID_BODY_NODE_H_
#define SimTK_SIMBODY_RIGID_BODY_NODE_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-12 Stanford University and the Authors.        *
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

#include "simbody/internal/common.h"
#include "SimbodyTreeState.h"

#include <cassert>

using SimTK::Vector;
using SimTK::Vector_;
using SimTK::Vec3;
using SimTK::SpatialVec;
using SimTK::Transform;
using SimTK::Rotation;
using SimTK::Inertia;
using SimTK::MassProperties;
using SimTK::Array_;

/* This is an abstract class representing the *computational* form of a body 
and its (generic) mobilizer, that is, the joint connecting it to its parent. 
Concrete classes are derived from this one to represent each specific type of 
mobilizer, with the emphasis being on fast calculation since node traversals 
are inner loops of all O(N) multibody algorithms. In particular, while the 
generic RigidBodyNode has a variable number of mobilities, each concrete node
has a compile-time known size so can perform inline floating point operations 
using stack-allocated variables.

Every body has a body frame B, and a unique inboard mobilizer frame M fixed to
that body. For convenience, we refer to the body frame of a body's unique 
parent body P as the P frame. There is a frame F fixed to P which is where 
B's inboard mobilizer attaches. The transform X_FM(q) tracks the across-
mobilizer change in configuration induced by the generalized coordinates q. 
When all the mobilizer coordinates are 0 (=1000 for quaternions), M and F take
on their "reference configuration" relationship. Usually M==F in the reference
configuration, but sometimes only the axes are aligned (X_FM(0).R()=Identity) 
while the origins are separated (X_FM(0).p() != 0); the ellipsoid joint is an 
example.

The mobilizer frame M is fixed with respect to B, and F is fixed with
respect to P. In some cases M and B or F and P will be the same, but not always.
The constant (TODO: Instance stage) transforms X_BM and X_PF provide the 
configuration of the mobilizer frames with respect to the body and parent 
frames. With these definitions we can easily calculate X_PB as 
X_PB = X_PF*X_FM*X_MB, where X_FM is the q-dependent cross-mobilizer transform
calculated at Position stage.

RigidBodyNodes know how to extract and deposit their own information from and
to the Simbody state variables and cache entries, but they don't know anything
about the State class, stages, etc. Instead they depend on being given 
appropriate access by the caller, whose job it is to mine the State.

The RigidBodyNode abstract base class defines the interface a concrete 
mobilized body node must implement. The nodes must implement two kinds of 
interface methods:
  - mobilizer-specific local computations
  - single-node contributions to multibody tree-sweeping algorithms
The latter methods are called in a well-defined sequence, sweeping either 
inwards or outwards so that a node can depend on either its child or parent 
computations (respectively) having already been completed.

Here is the spec (TODO):

Mobilizer-specific

  Mobility-dependent local kinematics calculations
  ------------------------------------------------
  calcQDot
  calcQDotDot
  calcMobilizerTransformFromQ
  calcMobilizerVelocityFromU
  calcMobilizerAccelerationFromUDot


  Intialization of state variables
  -------------------------------------
  setMobilizerDefaultModelValues
  setMobilizerDefaultInstanceValues
  setMobilizerDefaultTimeValues
  setMobilizerDefaultPositionValues
  setMobilizerDefaultVelocityValues
  setMobilizerDefaultDynamicsValues
  setMobilizerDefaultAccelerationValues

  Local initial condition setting (solvers)
  -----------------------------------------
  setQToFitTransform
  setQToFitRotation
  setQToFitTranslation
  setUToFitVelocity
  setUToFitAngularVelocity
  setUToFitLinearVelocity
  enforceQuaternionConstraints

  Bookkeeping
  -----------
  getName
  getNumU
  getNumQ
  isUsingQuaternion
*/
class RigidBodyNode {
public:

// Every concrete RigidBodyNode sets this property on construction so we know
// whether it can use simple default implementations for kinematic equations
// (those involving the N matrix), which rely on qdot=u.
enum QDotHandling {
    QDotIsAlwaysTheSameAsU, // can use default implementations
    QDotMayDifferFromU      // must supply custom implementation
};

// This is a fixed property of a given concrete RigidBodyNode.
bool isQDotAlwaysTheSameAsU() const
{   return qdotHandling==QDotIsAlwaysTheSameAsU; }

// This property tells us whether a RigidBodyNode may use a quaternion for
// some set of modeling options. If so, we know that nq > nu, and it must
// also be the case that QDotMayDifferFromU. (But note that there may be
// mobilizers for which qdot != u but there is no quaternion.)
enum QuaternionUse {
    QuaternionIsNeverUsed,
    QuaternionMayBeUsed
};

// This is a fixed property of a given concrete RigidBodyNode.
bool isQuaternionEverUsed() const
{   return quaternionUse==QuaternionMayBeUsed; }

virtual ~RigidBodyNode() {}

    // MOBILIZER-SPECIFIC VIRTUAL METHODS //

virtual const char* type()     const {return "unknown";}
virtual int  getDOF()   const=0; //number of independent dofs
virtual int  getMaxNQ() const=0; //dofs plus extra quaternion coordinate if any

virtual int getNQInUse(const SBModelVars&) const=0; //actual number of q's
virtual int getNUInUse(const SBModelVars&) const=0; // actual number of u's

// This depends on the mobilizer type and modeling options. If it returns
// true it also returns the first generalized coordinate of the quaternion
// (numbering locally), which is assumed to occupy that generalized coordinate
// and the next three contiguous ones. Quaternion-using mobilizers should be 
// assigned a slot in the quaternion pool.
virtual bool isUsingQuaternion(const SBStateDigest&, 
                               MobilizerQIndex& startOfQuaternion) const=0;

// This depends on the mobilizer type and modeling options. This is the amount
// of position-cache storage this mobilizer wants us to set aside for
// precalculations involving its q's, in units of number of Reals. The meaning
// of the entries in this pool is known only to the node.
virtual int calcQPoolSize(const SBModelVars&) const = 0;

// This operator is the first step in realizePosition() for this RBNode. Taking
// modeling information from the supplied state digest, perform calculations
// that can be done knowing only the current modeling parameters and the
// values of the q's. The results go into the posCache and qErr arrays which
// have already been set to the first entry belonging exclusively to this
// RBNode. DO NOT write directly into the provided State. (Typically 
// the output arrays will be referring to cache entries in the same State, but
// they don't have to.) Either of those output arrays may be null, and the 
// expected lengths are passed in for consistency checking -- at least in Debug
// mode you should make sure they are what you're expecting.
//
// The only mandatory calculation here is that this routine *must* calculate 
// constraint errors if there are any mobilizer-local constraints among
// the q's -- typically that is the quaternion normalization error where
// qerr=|q|-1. For many mobilizers, it is a good idea to precalculate
// expensive operations like sin(q), cos(q), or 1/norm(q) so that these don't
// need to be recalculated when forming X_FM, N, NInv, and NDot later.
// Note that if you want to save such precalculations, you must have asked
// for space to be set aside during realizeModel().
//
// The default implementation here checks the parameters and then does
// nothing.
virtual void performQPrecalculations(const SBStateDigest& sbs,
                                     const Real* q,  int nq,
                                     Real* posCache, int nPosCache,
                                     Real* qErr,     int nQErr) const=0;

// This mandatory routine calculates the across-joint transform X_FM generated
// by the supplied q values. This may depend on sines & cosines or normalized
// quaternions already having been calculated and placed into the supplied
// posCache. If you just have q's, call the calcAcrossJointTransform() method
// instead because it knows how to fill in the posCache if need be.
//
// This method constitutes the *definition* of the generalized coordinates for
// a particular mobilizer.
//
// Note: this calculates the transform between the *as defined* frames; we 
// might reverse the frames in use. So we call this X_F0M0 while the possibly 
// reversed version is X_FM.
virtual void calcX_FM(const SBStateDigest& sbs,
                      const Real* q,        int nq,
                      const Real* posCache, int nPos,
                      Transform&  X_F0M0) const = 0;

// This is a pure operator form of calcX_FM(). The StateDigest must have been
// realized to model stage. The result depends only on the passed-in q,
// not anything in the State beyond model stage.
void calcAcrossJointTransform(
    const SBStateDigest& sbs,
    const Real* q, int nq,
    Transform&  X_F0M0) const
{
    const int poolz = getQPoolSize(sbs.getModelCache());
    Real* qPool = poolz ? new Real[poolz] : 0;
    Real qErr;
    const int nQErr = isQuaternionInUse(sbs.getModelCache()) ? 1 : 0;
    performQPrecalculations(sbs, q, nq, qPool, poolz, &qErr, nQErr);
    calcX_FM(sbs, q, nq, qPool, poolz, X_F0M0);
    delete[] qPool;
}

// An alternate signature for when you have all the Q's together; this
// method will find just ours and call the above method.
void calcAcrossJointTransform(
    const SBStateDigest& sbs,
    const Vector&        q,
    Transform&           X_F0M0) const
{
    const SBModelCache mc = sbs.getModelCache();
    const int   nq = getNumQInUse(mc);
    const Real* qp = nq ? &q[getFirstQIndex(mc)] : (Real*)0;
    calcAcrossJointTransform(sbs, qp, nq, X_F0M0);
}


// This operator pulls N(q) from the StateDigest if necessary and calculates 
// qdot=N(q)*u from the supplied argument. For many mobilizers it 
// can simply copy u to qdot without referencing the state at all.
virtual void calcQDot
   (const SBStateDigest&,    const Real* u,    Real* qdot)    const 
{   SimTK_THROW2(Exception::UnimplementedVirtualMethod, "RigidBodeNode", "calcQDot");}
// This operator pulls N(q) and NDot(q,u) from the StateDigest if necessary
// and calculates qdotdot=N*udot + NDot*u from the supplied argument.
// For many mobilizers it can simply copy udot to qdotdot without referencing
// the state at all.
virtual void calcQDotDot
   (const SBStateDigest&,    const Real* udot, Real* qdotdot) const 
{   SimTK_THROW2(Exception::UnimplementedVirtualMethod, "RigidBodeNode", "calcQDotDot");}
virtual Transform  calcMobilizerTransformFromQ          
   (const SBStateDigest&,    const Real* q)    const 
{   SimTK_THROW2(Exception::UnimplementedVirtualMethod, "RigidBodeNode", "calcMobilizerTransformFromQ");}
virtual SpatialVec calcMobilizerVelocityFromU           
   (const SBStateDigest&,    const Real* u)    const 
{   SimTK_THROW2(Exception::UnimplementedVirtualMethod, "RigidBodeNode", "calcMobilizerVelocityFromU");}
virtual SpatialVec calcMobilizerAccelerationFromUDot    
   (const SBStateDigest&,    const Real* udot) const 
{   SimTK_THROW2(Exception::UnimplementedVirtualMethod, "RigidBodeNode", "calcMobilizerAccelerationFromUDot");}
virtual Transform  calcParentToChildTransformFromQ      
   (const SBStateDigest&,    const Real* q)    const 
{   SimTK_THROW2(Exception::UnimplementedVirtualMethod, "RigidBodeNode", "calcParentToChildTransformFromQ");}
virtual SpatialVec calcParentToChildVelocityFromU       
   (const SBStateDigest&,    const Real* u)    const 
{   SimTK_THROW2(Exception::UnimplementedVirtualMethod, "RigidBodeNode", "calcParentToChildVelocityFromU");}
virtual SpatialVec calcParentToChildAccelerationFromUDot
   (const SBStateDigest&,    const Real* udot) const 
{   SimTK_THROW2(Exception::UnimplementedVirtualMethod, "RigidBodeNode", "calcParentToChildAccelerationFromUDot");}

// Operators involving kinematics matrix N and related matrices NInv and
// NDot. These methods perform operations involving just the block on the 
// diagonal of the matrix that corresponds to this mobilizer. These blocks
// can be rectangular, in which case the u-dimension is always the number
// of mobilizers dofs (generalized speeds) but the q-dimension may depend on 
// modeling options (specifically, whether the mobilizer orientation is
// modeled with 4 quaternions or 3 Euler angles). These normally treat the
// input as a column vector and multiply with the matrix on the left.
// Optionally they will treat the input as a row vector and multiply with 
// the matrix on the right. (The latter is equivalent to multiplication on 
// the left by the transpose of the matrix.)
//
//   multiplyByN
//     Left   out_q = N * in_u
//     Right  out_u = in_q * N
//
//   multiplyByNInv
//     Left  out_u = inv(N) * in_q
//     Right out_q = in_u * inv(N)
//
//   multiplyByNDot
//     Left   out_q = NDot * in_u
//     Right  out_u = in_q * NDot
//      (not invertible but inverse isn't needed; see below)
//
// where 'in_q' and 'in_u' etc. indicate q-like and u-like vectors or row 
// vectors, not the actual state variables of the same names. Note that the
// interpretation of the in and out arrays is different depending on whether 
// we're multiplying on the left or right; which is q-like and which is
// u-like reverses.
//
// Note that N=N(q), NInv=Ninv(q), and NDot=NDot(q,u) where now I am talking 
// about the real set of generalized coordinates q and generalized speeds u
// but *just for this mobilizer*, whose values are taken from the supplied 
// StateDigest.
//
// In typical usage these are used to calculate qdot=N*u, 
// qdotdot=N*udot + NDot*u. These can be inverted to calculate u=inv(N)*qdot, 
// and udot=inv(N)*(qdotdot - NDot*u); note that inv(NDot) is not needed (a 
// good thing since it is likely to be singular!).
//
// WARNING: these routines do not normalize quaternions before using them, 
// so the resulting matrices are scaled by the quaternion norm (if it 
// isn't 1). That's a feature not a bug! This is important for correct 
// calculation of qdots because otherwise they would not be the correct 
// derivatives for unnormalized q's. 

virtual void multiplyByN   (const SBStateDigest&, bool matrixOnRight, 
                            const Real* in, Real* out) const = 0;
virtual void multiplyByNInv(const SBStateDigest&, bool matrixOnRight,
                            const Real* in, Real* out) const = 0;
virtual void multiplyByNDot(const SBStateDigest&, bool matrixOnRight,
                            const Real* in, Real* out) const = 0;

// This will do nothing unless the mobilizer is using a quaternion. Otherwise it
// will normalize its quaternion in q, and if qErrest has non-zero length then
// it will remove the component of the error estimate which was along the 
// direction of the quaternion, since that error will now be zero. That is, 
// we'll set
//     q = q_fixed = q/|q|
// and qErrest -= dot(qErrest,q_fixed)*q_fixed
virtual bool enforceQuaternionConstraints(
    const SBStateDigest& sbs,
    Vector&            q,
    Vector&            qErrest) const=0;

// Convert from quaternion to Euler angle representations.
virtual void convertToEulerAngles(const Vector& inputQ, Vector& outputQ) const 
{   SimTK_THROW2(Exception::UnimplementedVirtualMethod, "RigidBodeNode", "convertToEulerAngles");};
// Convert from Euler angle to quaternion representations.
virtual void convertToQuaternions(const Vector& inputQ, Vector& outputQ) const 
{   SimTK_THROW2(Exception::UnimplementedVirtualMethod, "RigidBodeNode", "convertToQuaternions");};

// Called after Model variables are allocated by realizeTopology()
virtual void setMobilizerDefaultModelValues
   (const SBTopologyCache&, SBModelVars&)        const {}

// All the rest are called right after realizeModel() since that's when all the
// remaining state variables are allocated.
virtual void setMobilizerDefaultInstanceValues    
   (const SBModelVars&,     SBInstanceVars&)     const {}
virtual void setMobilizerDefaultTimeValues        
   (const SBModelVars&,     SBTimeVars&)         const {}
virtual void setMobilizerDefaultPositionValues    
   (const SBModelVars&,     Vector& q)           const {}
virtual void setMobilizerDefaultVelocityValues    
   (const SBModelVars&,     Vector& u)           const {}
virtual void setMobilizerDefaultDynamicsValues    
   (const SBModelVars&,     SBDynamicsVars&)     const {}
virtual void setMobilizerDefaultAccelerationValues
   (const SBModelVars&,     SBAccelerationVars&) const {}

// These attempt to set the mobilizer's internal configuration or velocity
// to a specified value. This is intended to be a fast, local calculation that 
// produces an answer to machine precision *if* the mobilizer can represent the
// given quantity exactly. The answer is returned in the appropriate slots of a
// caller-provided "q-like" or "u-like" Vector; the implementation must not 
// look at or change any other slots.
//
// It is OK for the implementation to use the current values of the coordinates
// or speeds, and it is required to preserve any of these that are not needed 
// to satisfy the request. If the mobilizer can't satisfy the request to 
// machine precision it should just do nothing or the best it can, with the 
// only general rule being that it shouldn't make things worse. In particular,
// it does not need to work hard on an approximate solution.
//
// Note: these are non-virtual wrappers which arrange to reverse the request
// for reversed mobilizers, so that the mobilizers themselves do not need to
// know they have been reversed. The corresponding pure virtuals are protected.

void setQToFitTransform
   (const SBStateDigest& sbs, const Transform& X_FM, Vector& q) const
{   setQToFitTransformImpl(sbs, isReversed() ? Transform(~X_FM) : X_FM, q); }

void setQToFitRotation
   (const SBStateDigest& sbs, const Rotation& R_FM, Vector& q) const
{   setQToFitRotationImpl(sbs, isReversed() ? Rotation(~R_FM) : R_FM, q); }

// Reversing a translation requires that we obtain the current orientation,
// which we'll do assuming there is something meaningful in the rotational
// part of the q's.
void setQToFitTranslation
   (const SBStateDigest& sbs, const Vec3& p_FM, Vector& q) const
{   
    if (!isReversed()) {setQToFitTranslationImpl(sbs, p_FM, q); return;}

    Transform X_MF; // note reversal of frames
    calcAcrossJointTransform(sbs, q, X_MF);
    setQToFitTranslationImpl(sbs, X_MF.R()*(-p_FM), q);
}

void setUToFitVelocity
   (const SBStateDigest& sbs, const Vector& q, const SpatialVec& V_FM, Vector& u) const
{   if (!isReversed()) {setUToFitVelocityImpl(sbs, q, V_FM, u); return;}
    Transform X_MF; // note reversal of frames
    calcAcrossJointTransform(sbs, q, X_MF);
    setUToFitVelocityImpl(sbs, q, reverseSpatialVelocity(~X_MF,V_FM), u);
}

void setUToFitAngularVelocity
   (const SBStateDigest& sbs, const Vector& q, const Vec3& w_FM, Vector& u)       const
{   if (!isReversed()) {setUToFitAngularVelocityImpl(sbs, q, w_FM, u); return;}
    Transform X_MF; // note reversal of frames
    calcAcrossJointTransform(sbs, q, X_MF);
    setUToFitAngularVelocityImpl(sbs, q, reverseAngularVelocity(~X_MF.R(),w_FM), u);
}

void setUToFitLinearVelocity
   (const SBStateDigest& sbs, const Vector& q, const Vec3& v_FM, Vector& u) const
{   if (!isReversed()) {setUToFitLinearVelocityImpl(sbs, q, v_FM, u); return;}
    Transform X_MF; // note reversal of frames
    calcAcrossJointTransform(sbs, q, X_MF);
    //TODO: we have to assume angular velocity is zero here. Should be some
    // way to get the joint to compute w_FM.
    const SpatialVec V_FM( Vec3(0), v_FM );
    setUToFitLinearVelocityImpl(sbs, q, reverseSpatialVelocity(~X_MF,V_FM)[1], u);
}


    // VIRTUAL METHODS FOR SINGLE-NODE OPERATOR CONTRIBUTIONS //

virtual void realizeModel(
    SBStateDigest& sbs) const=0;

virtual void realizeInstance(
    const SBStateDigest& sbs) const=0;

// The multibody tree cannot have any time dependence; hence no
// realizeTime().

// Introduce new values for generalized coordinates and calculate
// all the position-dependent kinematic terms, including position
// constraint errors. Must be called base to tip.
virtual void realizePosition(
    const SBStateDigest&      sbs) const=0;

// Introduce new values for generalized speeds and calculate
// all the velocity-dependent kinematic terms. Assumes realizePosition()
// has already been called on all nodes. Must be called base to tip.
virtual void realizeVelocity(
    const SBStateDigest&         sbs) const=0;

// Calculate base-to-tip velocity-dependent terms which will be used
// in Dynamics stage operators. Assumes realizeVelocity()
// has already been called on all nodes, as well as any Dynamics 
// stage calculations which must go tip-to-base (e.g. articulated
// body inertias).
virtual void realizeDynamics(
    const SBArticulatedBodyInertiaCache&    abc,
    const SBStateDigest&                    sbs) const=0;

// There is no realizeAcceleration().

virtual void realizeReport(
    const SBStateDigest&         sbs) const=0;

virtual void realizeArticulatedBodyInertiasInward(
    const SBInstanceCache&          ic,
    const SBTreePositionCache&      pc,
    SBArticulatedBodyInertiaCache&  abc) const=0;

virtual void realizeYOutward(
    const SBInstanceCache&                ic,
    const SBTreePositionCache&            pc,
    const SBArticulatedBodyInertiaCache&  abc,
    SBDynamicsCache&                      dc) const=0;

// This has a default implementation that is good for everything
// but Ground.
virtual void calcCompositeBodyInertiasInward(
    const SBTreePositionCache&                  pc,
    Array_<SpatialInertia,MobilizedBodyIndex>&  R) const;

virtual void multiplyBySystemJacobian(
    const SBTreePositionCache&  pc,
    const Real*                 v,
    SpatialVec*                 Jv) const
  { SimTK_THROW2(Exception::UnimplementedVirtualMethod, "RigidBodeNode", 
    "multiplyBySystemJacobian"); }

virtual void multiplyBySystemJacobianTranspose(
    const SBTreePositionCache&  pc, 
    SpatialVec*                 zTmp,
    const SpatialVec*           X, 
    Real*                       JtX) const
  { SimTK_THROW2(Exception::UnimplementedVirtualMethod, "RigidBodeNode", 
    "multiplyBySystemJacobianTranspose"); }

virtual void calcEquivalentJointForces(
    const SBTreePositionCache&  pc,
    const SBDynamicsCache&      dc,
    const SpatialVec*           bodyForces,
    SpatialVec*                 allZ,
    Real*                       jointForces) const
  { SimTK_THROW2(Exception::UnimplementedVirtualMethod, "RigidBodeNode", "calcEquivalentJointForces"); }

virtual void calcUDotPass1Inward(
    const SBInstanceCache&                  ic,
    const SBTreePositionCache&              pc,
    const SBArticulatedBodyInertiaCache&    abc,
    const SBDynamicsCache&                  dc,
    const Real*                             jointForces,
    const SpatialVec*                       bodyForces,
    const Real*                             allUDot,
    SpatialVec*                             allZ,
    SpatialVec*                             allGepsilon,
    Real*                                   allEpsilon) const=0;
virtual void calcUDotPass2Outward(
    const SBInstanceCache&                  ic,
    const SBTreePositionCache&              pc,
    const SBArticulatedBodyInertiaCache&    abc,
    const SBTreeVelocityCache&              vc,
    const SBDynamicsCache&                  dc,
    const Real*                             epsilonTmp,
    SpatialVec*                             allA_GB,
    Real*                                   allUDot,
    Real*                                   allTau) const=0;

virtual void multiplyByMInvPass1Inward(
    const SBInstanceCache&                  ic,
    const SBTreePositionCache&              pc,
    const SBArticulatedBodyInertiaCache&    abc,
    const Real*                             f,
    SpatialVec*                             allZ,
    SpatialVec*                             allGepsilon,
    Real*                                   allEpsilon) const=0;
virtual void multiplyByMInvPass2Outward(
    const SBInstanceCache&                  ic,
    const SBTreePositionCache&              pc,
    const SBArticulatedBodyInertiaCache&    abc,
    const Real*                             epsilonTmp,
    SpatialVec*                             allA_GB,
    Real*                                   allUDot) const=0;

// Also serves as pass 1 for inverse dynamics.
virtual void calcBodyAccelerationsFromUdotOutward(
    const SBTreePositionCache&  pc,
    const SBTreeVelocityCache&  vc,
    const Real*                 allUDot,
    SpatialVec*                 allA_GB) const
  { SimTK_THROW2(Exception::UnimplementedVirtualMethod, "RigidBodeNode", "calcBodyAccelerationsFromUdotOutward"); }
virtual void calcInverseDynamicsPass2Inward(
    const SBTreePositionCache&  pc,
    const SBTreeVelocityCache&  vc,
    const SpatialVec*           allA_GB,
    const Real*                 jointForces,
    const SpatialVec*           bodyForces,
    SpatialVec*                 allFTmp,
    Real*                       allTau) const
  { SimTK_THROW2(Exception::UnimplementedVirtualMethod, "RigidBodeNode", "calcInverseDynamicsPass2Inward"); }

virtual void multiplyByMPass1Outward(
    const SBTreePositionCache&  pc,
    const Real*                 allUDot,
    SpatialVec*                 allA_GB) const
  { SimTK_THROW2(Exception::UnimplementedVirtualMethod, "RigidBodeNode", "multiplyByMPass1Outward"); }
virtual void multiplyByMPass2Inward(
    const SBTreePositionCache&  pc,
    const SpatialVec*           allA_GB,
    SpatialVec*                 allFTmp,
    Real*                       allTau) const
  { SimTK_THROW2(Exception::UnimplementedVirtualMethod, "RigidBodeNode", "multiplyByMPass2Inward"); }


virtual void setVelFromSVel(const SBStateDigest&,
                            const SpatialVec&, Vector& u) const {SimTK_THROW2(Exception::UnimplementedVirtualMethod, "RigidBodeNode", "setVelFromSVel");}

// Note that this requires columns of H to be packed like SpatialVec.
virtual const SpatialVec& getHCol(const SBTreePositionCache&, int j) const 
{SimTK_THROW2(Exception::UnimplementedVirtualMethod, "RigidBodeNode", "getHCol");}

virtual const SpatialVec& getH_FMCol(const SBTreePositionCache&, int j) const 
{SimTK_THROW2(Exception::UnimplementedVirtualMethod, "RigidBodeNode", "getH_FMCol");}

//TODO (does this even belong here?)
virtual void velFromCartesian() {}


    // BASE CLASS METHODS //

RigidBodyNode& operator=(const RigidBodyNode&);

static RigidBodyNode* createGroundNode();

// Register the passed-in node as a child of this one.
void addChild(RigidBodyNode* child);
void setParent(RigidBodyNode* p) {parent=p;}
void setNodeNum(MobilizedBodyIndex n) {nodeNum=n;}
void setLevel(int i)   {level=i;}

    // TOPOLOGICAL INFO: no State needed

RigidBodyNode* getParent() const {return parent;}
int            getNumChildren()  const {return (int)children.size();}
RigidBodyNode* getChild(int i) const {return (i<(int)children.size()?children[i]:0);}

bool           isReversed() const {return reversed;}

// Return this node's level, that is, how many ancestors separate it from
// the Ground node at level 0. Level 1 nodes (directly connected to the
// Ground node) are called 'base' nodes.
int  getLevel() const  {return level;}

// This is the unique "body index" for this (body,mobilizer) node. It is used to index
// arrays of body quantities.
MobilizedBodyIndex getNodeNum() const {return nodeNum;}

bool isGroundNode() const { return level==0; }
bool isBaseNode()   const { return level==1; }

// TODO: these should come from the model cache.
UIndex getUIndex() const {return uIndex;}
QIndex getQIndex() const {return qIndex;}


// Access routines for plucking the right per-body data from the pool in the State.
const Transform&  fromB(const Array_<Transform,MobilizedBodyIndex>& x) const {return x[nodeNum];}
const PhiMatrix&  fromB(const Array_<PhiMatrix,MobilizedBodyIndex>& p) const {return p[nodeNum];}
const MassProperties& fromB(const Array_<MassProperties,MobilizedBodyIndex>& m) const {return m[nodeNum];}
const Inertia&    fromB(const Array_<Inertia,MobilizedBodyIndex>&   i) const {return i[nodeNum];}
const UnitInertia& fromB(const Array_<UnitInertia,MobilizedBodyIndex>&   i) const {return i[nodeNum];}
int               fromB(const Array_<int,MobilizedBodyIndex>&       i) const {return i[nodeNum];}
const SpatialVec& fromB(const Array_<SpatialVec,MobilizedBodyIndex>&    v) const {return v[nodeNum];}
const SpatialMat& fromB(const Array_<SpatialMat,MobilizedBodyIndex>&    m) const {return m[nodeNum];}
const SpatialInertia& fromB(const Array_<SpatialInertia,MobilizedBodyIndex>& m) const {return m[nodeNum];}
const ArticulatedInertia& fromB(const Array_<ArticulatedInertia,MobilizedBodyIndex>& m) const {return m[nodeNum];}
const Vec3&       fromB(const Array_<Vec3,MobilizedBodyIndex>&          v) const {return v[nodeNum];}


Transform&  toB(Array_<Transform,MobilizedBodyIndex>& x) const {return x[nodeNum];}
PhiMatrix&  toB(Array_<PhiMatrix,MobilizedBodyIndex>& p) const {return p[nodeNum];}
MassProperties& toB(Array_<MassProperties,MobilizedBodyIndex>& m) const {return m[nodeNum];}
Inertia&    toB(Array_<Inertia,MobilizedBodyIndex>&   i) const {return i[nodeNum];}
UnitInertia& toB(Array_<UnitInertia,MobilizedBodyIndex>&   i) const {return i[nodeNum];}
int&        toB(Array_<int,MobilizedBodyIndex>&       i) const {return i[nodeNum];}
SpatialVec& toB(Array_<SpatialVec,MobilizedBodyIndex>&    v) const {return v[nodeNum];}
SpatialMat& toB(Array_<SpatialMat,MobilizedBodyIndex>&    m) const {return m[nodeNum];}
SpatialInertia& toB(Array_<SpatialInertia,MobilizedBodyIndex>& m) const {return m[nodeNum];}
ArticulatedInertia& toB(Array_<ArticulatedInertia,MobilizedBodyIndex>& m) const {return m[nodeNum];}
Vec3&       toB(Array_<Vec3,MobilizedBodyIndex>&          v) const {return v[nodeNum];}

// Elementwise access to Vectors is relatively expensive; use sparingly.
const SpatialVec& fromB(const Vector_<SpatialVec>&    v) const {return v[nodeNum];}
const SpatialMat& fromB(const Vector_<SpatialMat>&    m) const {return m[nodeNum];}
const Vec3&       fromB(const Vector_<Vec3>&          v) const {return v[nodeNum];}
SpatialMat& toB(Vector_<SpatialMat>&    m) const {return m[nodeNum];}
Vec3&       toB(Vector_<Vec3>&          v) const {return v[nodeNum];}
SpatialVec& toB(Vector_<SpatialVec>&    v) const {return v[nodeNum];}

    // MODELING INFO
bool getUseEulerAngles(const SBModelVars& mv) const {return mv.useEulerAngles;}

// Find cache resources allocated to this RigidBodyNode.
const SBModelPerMobodInfo& 
getModelInfo(const SBModelCache& mc) const
{   return mc.getMobodModelInfo(nodeNum); }

QIndex getFirstQIndex(const SBModelCache& mc) const
{   return getModelInfo(mc).firstQIndex; }
int getNumQInUse(const SBModelCache& mc) const
{   return getModelInfo(mc).nQInUse; }

UIndex getFirstUIndex(const SBModelCache& mc) const
{   return getModelInfo(mc).firstUIndex; }
int getNumUInUse(const SBModelCache& mc) const
{   return getModelInfo(mc).nUInUse; }

// Return value will be invalid if there are no quaternions currently in use here.
QuaternionPoolIndex getQuaternionPoolIndex(const SBModelCache& mc) const 
{   return getModelInfo(mc).quaternionPoolIndex; }
bool isQuaternionInUse(const SBModelCache& mc) const
{   return getModelInfo(mc).hasQuaternionInUse; }

// Return value will be invalid if this node is not currently using any space 
// in the q-pool position cache entry.
MobodQPoolIndex getQPoolIndex(const SBModelCache& mc) const 
{   return getModelInfo(mc).startInQPool; }
int getQPoolSize(const SBModelCache& mc) const
{   return getModelInfo(mc).nQPoolInUse; }

// Return a pointer to the first slot in the q cache pool that is
// assigned to this node, or null if there are none.
const Real* getQPool(const SBModelCache& mc,
                       const SBTreePositionCache& pc) const
{   const MobodQPoolIndex ix = getQPoolIndex(mc);
    return ix.isValid() ? &pc.mobilizerQCache[ix] : (Real*)0; }

    // INSTANCE INFO

// TODO: These ignore State currently since they aren't parametrizable.
const MassProperties& getMassProperties_OB_B() const {return massProps_B;}
const Real&           getMass          () const {return massProps_B.getMass();}
const Vec3&           getCOM_B         () const {return massProps_B.getMassCenter();}
const UnitInertia&    getUnitInertia_OB_B() const {return massProps_B.getUnitInertia();}
const Transform&      getX_BM          () const {return X_BM;}
const Transform&      getX_PF          () const {return X_PF;}

// These are calculated on construction.
// TODO: principal axes
const Inertia&        getInertia_CB_B  () const {return inertia_CB_B;}
const Transform&      getX_MB          () const {return X_MB;}

// This says whether this mobilizer has prescribed *acceleration*, and if so
// whether the acceleration is known to be zero.
bool isUDotKnown(const SBInstanceCache& ic) const 
{   return ic.getMobodInstanceInfo(nodeNum).udotMethod!=Motion::Free; }
bool isUDotKnownToBeZero(const SBInstanceCache& ic) const 
{   return ic.getMobodInstanceInfo(nodeNum).udotMethod==Motion::Zero; }

    // POSITION INFO

// Extract from the cache  X_FM, the cross-mobilizer transformation matrix giving the configuration
// of this body's mobilizer frame M, measured from and expressed in the corresponding outboard
// mobilizer frame F attached to the parent. This transformation is defined to be zero (that is, F=M)
// in the reference configuration where the joint coordinates are all 0 (or 1,0,0,0 for quaternions).
// This is NOT a spatial (ground frame) transformation.
const Transform& getX_FM(const SBTreePositionCache& pc) const {return fromB(pc.bodyJointInParentJointFrame);}
Transform&       updX_FM(SBTreePositionCache&       pc) const {return toB  (pc.bodyJointInParentJointFrame);}

// Return X_F0M0, the cross-joint mobilizer transform *as it was defined* in the mobilizer 
// specification. If the mobilizer has been reversed, then X_F0M0=~X_FM, otherwise it is 
// just X_FM.
Transform findX_F0M0(const SBTreePositionCache& pc) const 
{   return isReversed() ? Transform(~getX_FM(pc))   // 18 flops
                        : getX_FM(pc); }

// Extract from the cache  X_PB, the cross-joint transformation matrix giving the configuration
// of this body's frame B measured from and expressed in its *parent* frame P. Thus this is NOT
// a spatial (ground frame) transformation.
const Transform& getX_PB(const SBTreePositionCache& pc) const {return fromB(pc.bodyConfigInParent);}
Transform&       updX_PB(SBTreePositionCache&       pc) const {return toB  (pc.bodyConfigInParent);}

// Extract from the cache X_GB, the transformation matrix giving the spatial configuration of this
// body's frame B measured from and expressed in ground. This consists of a rotation matrix
// R_GB, and a ground-frame vector r_OG_OB from ground's origin to the origin point of frame B.
const Transform& getX_GB(const SBTreePositionCache& pc) const {
    return fromB(pc.bodyConfigInGround);
}
Transform& updX_GB(SBTreePositionCache& pc) const {
    return toB(pc.bodyConfigInGround);
}

// Extract from the cache the body-to-parent shift matrix "phi". 
const PhiMatrix& getPhi(const SBTreePositionCache& pc) const {return fromB(pc.bodyToParentShift);}
PhiMatrix&       updPhi(SBTreePositionCache&       pc) const {return toB  (pc.bodyToParentShift);}

// Extract this body's spatial inertia matrix from the cache. This contains 
// the mass properties measured from (and about) the body frame origin, but 
// expressed in the Ground frame.
const SpatialInertia& getMk_G(const SBTreePositionCache& pc) const {return fromB(pc.bodySpatialInertiaInGround);}
SpatialInertia&       updMk_G(SBTreePositionCache&       pc) const {return toB  (pc.bodySpatialInertiaInGround);}

// Extract from the cache the location of the body's center of mass, measured from the ground
// origin and expressed in Ground.
const Vec3& getCOM_G(const SBTreePositionCache& pc) const {return fromB(pc.bodyCOMInGround);}
Vec3&       updCOM_G(SBTreePositionCache&       pc) const {return toB  (pc.bodyCOMInGround);}

// Extract from the cache the vector from body B's origin to its center of 
// mass, reexpressed in Ground.
const Vec3& getCB_G(const SBTreePositionCache& pc) const 
{   return getMk_G(pc).getMassCenter(); }

// Extract from the cache the body's inertia about the body origin OB, but 
// reexpressed in Ground.
const UnitInertia& getUnitInertia_OB_G(const SBTreePositionCache& pc) const 
{   return getMk_G(pc).getUnitInertia(); }

// Extract from the cache the spatial (ground-relative) location and orientation of this body's
// *parent's* body frame P.
const Transform& getX_GP(const SBTreePositionCache& pc) const {assert(parent); return parent->getX_GB(pc);}

        // VELOCITY INFO

const SpatialVec& getV_FM(const SBTreeVelocityCache& vc) const {return fromB(vc.mobilizerRelativeVelocity);}
SpatialVec&       updV_FM(SBTreeVelocityCache&       vc) const {return toB  (vc.mobilizerRelativeVelocity);}

// Given V_AB, the spatial velocity of frame B in A expressed in A, return V_BA, the spatial 
// velocity of frame A in B, expressed in B. The reversal also requires knowing X_AB, the spatial 
// position and orientation of B in A. Theory:
//      V_BA = -R_BA * [        w_AB        ]
//                     [ v_AB + p_AB x w_AB ]
// This costs 45 flops. See reverseAngularVelocity() if that's all you need; it's a lot cheaper.
static SpatialVec reverseSpatialVelocity(const Transform& X_AB, const SpatialVec& V_AB) {
    const Rotation& R_AB = X_AB.R();
    const Vec3&     p_AB = X_AB.p();
    const Vec3&     w_AB = V_AB[0];
    const Vec3&     v_AB = V_AB[1];

    return ~R_AB * SpatialVec( -w_AB,   (w_AB % p_AB) - v_AB );
}

// Given w_AB, the angular velocity of frame B in A expressed in A, return w_BA, the angular 
// velocity of frame A in B, expressed in B. The reversal also requires knowing R_AB, the 
// orientation of B in A. Theory:
//      w_BA = -R_BA*w_AB
// This is just a subset of what reverseSpatialVelocity does, but it is cheap and
// often all you need is the angular velocity.
// This costs 18 flops.
static Vec3 reverseAngularVelocity(const Rotation& R_AB, const Vec3& w_AB) {
    return ~R_AB * (-w_AB);
}

// Return V_F0M0, the cross-joint mobilizer velocity *as it was defined* in the mobilizer 
// specification. If the mobilizer has not been reversed, this is just V_FM.
// 45 flops if reversed.
SpatialVec findV_F0M0(const SBTreePositionCache& pc, const SBTreeVelocityCache& vc) const {
    return isReversed() ? reverseSpatialVelocity(getX_FM(pc), getV_FM(vc))
                        : getV_FM(vc);
}

// Return w_F0M0, the cross-joint mobilizer angular velocity *as it was defined* 
// in the mobilizer specification. If the mobilizer has not been reversed, this is just w_FM.
// This is a useful and much cheaper subset of findV_F0M0.
// 18 flops if reversed.
Vec3 find_w_F0M0(const SBTreePositionCache& pc, const SBTreeVelocityCache& vc) const {
    return isReversed() ? reverseAngularVelocity(getX_FM(pc).R(), getV_FM(vc)[0])
                        : getV_FM(vc)[0];
}

// Extract from the cache V_GB, the spatial velocity of this body's frame B measured in and
// expressed in ground. This contains the angular velocity of B in G, and the linear velocity
// of B's origin point OB in G, with both vectors expressed in G.
const SpatialVec& getV_GB   (const SBTreeVelocityCache& vc) const {return fromB(vc.bodyVelocityInGround);}
SpatialVec&       updV_GB   (SBTreeVelocityCache&       vc) const {return toB  (vc.bodyVelocityInGround);}

// Extract from the cache V_PB_G, the *spatial* velocity of this body's frame B, that is the
// cross-joint velocity measured with respect to the parent frame, but then expressed in the
// *ground* frame. This contains the angular velocity of B in P, and the linear velocity
// of B's origin point OB in P, with both vectors expressed in *G*.
const SpatialVec& getV_PB_G (const SBTreeVelocityCache& vc) const {return fromB(vc.bodyVelocityInParent);}
SpatialVec&       updV_PB_G (SBTreeVelocityCache&       vc) const {return toB  (vc.bodyVelocityInParent);}

const SpatialVec& getV_GP(const SBTreeVelocityCache& vc) const {assert(parent); return parent->getV_GB(vc);}

const SpatialVec& getSpatialVel   (const SBTreeVelocityCache& vc) const {return getV_GB(vc);}
const Vec3&       getSpatialAngVel(const SBTreeVelocityCache& vc) const {return getV_GB(vc)[0];}
const Vec3&       getSpatialLinVel(const SBTreeVelocityCache& vc) const {return getV_GB(vc)[1];}

// Extract from the cache VD_PB_G, the *spatial* velocity derivative remainder term
// generated by H_PB_G_Dot*u, where H_PB_G_Dot = d/dt H_PB_G with the derivative taken
// in the Ground frame. This is used in calculation of coriolis acceleration.
// CAUTION: our definition for the H matrix is transposed from those used by Jain and Schwieters.
const SpatialVec& getVD_PB_G (const SBTreeVelocityCache& vc) const 
    {return fromB(vc.bodyVelocityInParentDerivRemainder);}
SpatialVec&       updVD_PB_G (SBTreeVelocityCache&       vc) const 
    {return toB  (vc.bodyVelocityInParentDerivRemainder);}

const SpatialVec& getGyroscopicForce(const SBTreeVelocityCache& vc) const {return fromB(vc.gyroscopicForces);}
SpatialVec&       updGyroscopicForce(SBTreeVelocityCache&       vc) const {return toB  (vc.gyroscopicForces);}

const SpatialVec& getMobilizerCoriolisAcceleration(const SBTreeVelocityCache& vc) const {return fromB(vc.mobilizerCoriolisAcceleration);}
SpatialVec&       updMobilizerCoriolisAcceleration(SBTreeVelocityCache&       vc) const {return toB  (vc.mobilizerCoriolisAcceleration);}

const SpatialVec& getTotalCoriolisAcceleration(const SBTreeVelocityCache& vc) const {return fromB(vc.totalCoriolisAcceleration);}
SpatialVec&       updTotalCoriolisAcceleration(SBTreeVelocityCache&       vc) const {return toB  (vc.totalCoriolisAcceleration);}

    // DYNAMICS INFO

// Composite body inertias.
const SpatialInertia& getR(const SBCompositeBodyInertiaCache& cbc) const {return fromB(cbc.compositeBodyInertia);}
SpatialInertia&       updR(SBCompositeBodyInertiaCache&       cbc) const {return toB  (cbc.compositeBodyInertia);}

// Articulated body inertias and related calculations
const ArticulatedInertia& getP(const SBArticulatedBodyInertiaCache& abc) const {return fromB(abc.articulatedBodyInertia);}
ArticulatedInertia&       updP(SBArticulatedBodyInertiaCache&       abc) const {return toB  (abc.articulatedBodyInertia);}

const ArticulatedInertia& getPPlus(const SBArticulatedBodyInertiaCache& abc) const {return fromB(abc.pPlus);}
ArticulatedInertia&       updPPlus(SBArticulatedBodyInertiaCache&       abc) const {return toB  (abc.pPlus);}

// Others

const SpatialVec& getMobilizerCentrifugalForces(const SBDynamicsCache& dc) const {return fromB(dc.mobilizerCentrifugalForces);}
SpatialVec&       updMobilizerCentrifugalForces(SBDynamicsCache&       dc) const {return toB  (dc.mobilizerCentrifugalForces);}

const SpatialVec& getTotalCentrifugalForces(const SBDynamicsCache& dc) const {return fromB(dc.totalCentrifugalForces);}
SpatialVec&       updTotalCentrifugalForces(SBDynamicsCache&       dc) const {return toB  (dc.totalCentrifugalForces);}

const SpatialVec& getZ(const SBTreeAccelerationCache& tac) const {return fromB(tac.z);}
SpatialVec&       updZ(SBTreeAccelerationCache&       tac) const {return toB  (tac.z);}

const SpatialVec& getZPlus(const SBTreeAccelerationCache& tac) const {return fromB(tac.zPlus);}
SpatialVec&       updZPlus(SBTreeAccelerationCache&       tac) const {return toB  (tac.zPlus);}


const SpatialMat& getY(const SBDynamicsCache& dc) const {return fromB(dc.Y);}
SpatialMat&       updY(SBDynamicsCache&       dc) const {return toB  (dc.Y);}

    // ACCELERATION INFO

// Extract from the cache A_GB, the spatial acceleration of this body's frame B measured in and
// expressed in ground. This contains the inertial angular acceleration of B in G, and the
// linear acceleration of B's origin point OB in G, with both vectors expressed in G.
const SpatialVec& getA_GB (const SBTreeAccelerationCache& ac) const {return fromB(ac.bodyAccelerationInGround);}
SpatialVec&       updA_GB (SBTreeAccelerationCache&       ac) const {return toB  (ac.bodyAccelerationInGround);}

const SpatialVec& getSpatialAcc   (const SBTreeAccelerationCache& ac) const {return getA_GB(ac);}
const Vec3&       getSpatialAngAcc(const SBTreeAccelerationCache& ac) const {return getA_GB(ac)[0];}
const Vec3&       getSpatialLinAcc(const SBTreeAccelerationCache& ac) const {return getA_GB(ac)[1];}



// These are called just after new state variables are allocated,
// in case there are any node-specific default values.
// We can handle the "body" variables (like mass) here, but we have to forward the
// request to the mobilizers to handle their own variables. At the Position
// stage, for example, mobilizers which use quaternions will set the default ball 
// joint q's to 1,0,0,0.
// Most of these will use the default implementations here, i.e. do nothing.

// Called after Model variables are allocated by realizeTopology()
void setNodeDefaultModelValues(const SBTopologyCache& tc, SBModelVars& mv) const {
    // no body model variables at the moment; TODO: should forward to Body type
    setMobilizerDefaultModelValues(tc, mv);
}

// All the rest are called right after realizeModel() since that's when all the
// remaining state variables are allocated.
void setNodeDefaultInstanceValues(const SBModelVars& mv, SBInstanceVars& iv) const {
    // mass properties, inb and outb frame are handled here
    toB(iv.bodyMassProperties)      = getMassProperties_OB_B();
    toB(iv.outboardMobilizerFrames) = getX_BM();
    toB(iv.inboardMobilizerFrames)  = getX_PF();
    setMobilizerDefaultInstanceValues(mv,iv);
}
void setNodeDefaultTimeValues(const SBModelVars& mv, SBTimeVars& tv)  const {
    // no body model variables at the moment; TODO: should forward to Body type
    setMobilizerDefaultTimeValues(mv, tv);
}
void setNodeDefaultPositionValues(const SBModelVars& mv, Vector& q) const {
    // no body model variables at the moment; TODO: should forward to Body type
    setMobilizerDefaultPositionValues(mv, q);
}
void setNodeDefaultVelocityValues(const SBModelVars& mv, Vector& u) const {
     // no body model variables at the moment; TODO: should forward to Body type
    setMobilizerDefaultVelocityValues(mv, u);
}
void setNodeDefaultDynamicsValues(const SBModelVars& mv, SBDynamicsVars& dv) const {
    // no body model variables at the moment; TODO: should forward to Body type
    setMobilizerDefaultDynamicsValues(mv, dv);
}
void setNodeDefaultAccelerationValues(const SBModelVars& mv, SBAccelerationVars& av) const {
     // no body model variables at the moment; TODO: should forward to Body type
    setMobilizerDefaultAccelerationValues(mv, av);
}


// Calculate kinetic energy (from spatial quantities only).
Real calcKineticEnergy(
    const SBTreePositionCache& pc,
    const SBTreeVelocityCache& vc) const;   

// Calculate all spatial configuration quantities, assuming availability of
// joint-specific relative quantities.
void calcJointIndependentKinematicsPos(
    SBTreePositionCache& pc) const;

// Calcluate all spatial velocity quantities, assuming availability of
// joint-specific relative quantities and all position kinematics.
void calcJointIndependentKinematicsVel(
    const SBTreePositionCache& pc,
    SBTreeVelocityCache&       vc) const;

// Calculate velocity-dependent quantities which will be needed for
// computing accelerations.
void calcJointIndependentDynamicsVel(
    const SBTreePositionCache&              pc,
    const SBArticulatedBodyInertiaCache&    abc,
    const SBTreeVelocityCache&              vc,
    SBDynamicsCache&                        dc) const;


protected:
// This is the constructor for the abstract base type for use by the derived
// concrete types in their constructors.
RigidBodyNode(const MassProperties& mProps_B,
              const Transform&      xform_PF,
              const Transform&      xform_BM,
              QDotHandling          qdotType,
              QuaternionUse         quatUse,
              bool                  reverse=false)
  : parent(0), children(), level(-1),
    massProps_B(mProps_B), 
    inertia_CB_B(mProps_B.isFinite()
                 ? mProps_B.calcCentralInertia()
                 : (mProps_B.isInf() ? Inertia(Infinity) : Inertia())),
    X_BM(xform_BM), X_PF(xform_PF), X_MB(~xform_BM),
    qdotHandling(qdotType), quaternionUse(quatUse), reversed(reverse)
{
    // If a quaternion might be used, it can't possibly be true that qdot is
    // always the same as u.
    assert(!(quatUse==QuaternionMayBeUsed && qdotType==QDotIsAlwaysTheSameAsU));
}

// These have wrappers above which call them and should not be
// called directly from outside this class and its descendents. Note that
// the requests are couched in terms of the way the mobilizer was
// *defined*, not as it is after reversal (hence F0 instead of F and
// M0 instead of M). The q's and u's have identical meanings for both
// forward and reverse mobilizers.
virtual void setQToFitTransformImpl
   (const SBStateDigest&, const Transform& X_F0M0, Vector& q) const = 0;
virtual void setQToFitRotationImpl
   (const SBStateDigest&, const Rotation& R_F0M0, Vector& q) const = 0;
virtual void setQToFitTranslationImpl
   (const SBStateDigest&, const Vec3& p_F0M0, Vector& q) const = 0;

virtual void setUToFitVelocityImpl
   (const SBStateDigest&, const Vector& q, const SpatialVec& V_F0M0, Vector& u) const = 0;
virtual void setUToFitAngularVelocityImpl
   (const SBStateDigest&, const Vector& q, const Vec3& w_F0M0, Vector& u)       const = 0;
virtual void setUToFitLinearVelocityImpl
   (const SBStateDigest&, const Vector& q, const Vec3& v_F0M0, Vector& u) const = 0;

typedef Array_<RigidBodyNode*> RigidBodyNodeList;

QIndex              qIndex;   // index into internal coord pos array
UIndex              uIndex;   // index into internal coord vel,acc arrays
USquaredIndex       uSqIndex; // index into array of DOF^2 objects
QuaternionPoolIndex quaternionIndex; // if this mobilizer has a quaternion, this is our slot

RigidBodyNode*     parent; 
RigidBodyNodeList  children;
int                level;        //how far from base 
MobilizedBodyIndex nodeNum;      //unique ID number in SimbodyMatterSubsystemRep

// These are the default body properties, all supplied or calculated on
// construction. TODO: they should be 
// (optionally?) overrideable by Instance-level state variable entries.

// This is the mass, center of mass, and inertia as supplied at construction.
// Here the inertia is taken about the B origin OB.
const MassProperties massProps_B;

// This is the supplied inertia, shifted to the center of mass. It is still
// a constant expressed in B, but is taken about the COM.
const Inertia   inertia_CB_B;

// Orientation and location of inboard mobilizer frame M, measured
// and expressed in body frame B.
const Transform X_BM; 
const Transform X_MB; // inverse of X_BM, calculated on construction

// This is set when we attach this node to its parent in the tree. This is the
// configuration of the parent's outboard mobilizer attachment frame corresponding
// to body B (F) measured from and expressed in the parent frame P. It is 
// a constant in frame P. TODO: make it parameterizable.
const Transform X_PF;

// Concrete RigidBodyNodes should set this flag on construction to indicate whether they can guarantee
// that their mobilizer's qdots are just the generalized speeds u, for all possible
// modeling options. That is the same as saying nq=nu and qdot=u, or that the block of the kinematic
// matrix N corresponding to this node is an nuXnu identity matrix; NInv is the same; and
// NDot is nuXnu zero. These conditions allows the use of default implementations of many
// methods which would otherwise have to be overridden. The default implementations assert
// that qdotHandling==QDotIsAlwaysTheSameAsU.
const QDotHandling qdotHandling;

// Concrete RigidBodyNodes should set this flag on construction to indicate whether they
// can guarantee never to use a quaternion in their generailized coordinates for any set
// of modeling options.
const QuaternionUse quaternionUse;

// This is true if the mobilizer is specified in the reverse direction,
// that is from child to parent.
const bool reversed;


};

#endif // SimTK_SIMBODY_RIGID_BODY_NODE_H_
