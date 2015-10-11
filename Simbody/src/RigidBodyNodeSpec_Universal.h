#ifndef SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_UNIVERSAL_H_
#define SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_UNIVERSAL_H_

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

/**@file
 * Define the RigidBodyNode that implements a Universal mobilizer, also known
 * as a UJoint or Hooke's joint.
 */

#include "SimbodyMatterSubsystemRep.h"
#include "RigidBodyNode.h"
#include "RigidBodyNodeSpec.h"

    // UNIVERSAL (U-JOINT, HOOKE'S JOINT) //

// This is a Universal Joint (U-Joint), also known as a Hooke's joint. This is
// identical to the joint that connects pieces of a driveshaft together under
// a car. Physically, you can think of this as a parent body P, hinged to a
// massless cross-shaped coupler, which is then hinged to the child body B. The
// massless coupler doesn't actually appear in the model. Instead, we use a
// body-fixed 1-2 Euler rotation sequence for orientation, which has the same
// effect: starting with frames B and P aligned (when q0=q1=0), rotate frame
// B about the Px(=Bx) axis by q0; then, rotate frame B further about the new
// By(!=Py) by q1. For generalized speeds u we use the Euler angle derivatives
// qdot, which are *not* the same as angular velocity components because u0 is
// a rotation rate around Px(!=Bx any more) while u1 is a rotation rate
// about By.
//
// To summarize,
//    q's: a two-angle body-fixed rotation sequence about x, then new y
//    u's: time derivatives of the q's
//
// Note that the U-Joint degrees of freedom relating the parent's F frame to
// the child's M frame are about x and y, with the "long" axis of the
// driveshaft along z.
template<bool noX_MB, bool noR_PF>
class RBNodeUJoint : public RigidBodyNodeSpec<2, false, noX_MB, noR_PF> {
public:
typedef typename RigidBodyNodeSpec<2, false, noX_MB, noR_PF>::HType HType;
virtual const char* type() { return "ujoint"; }

RBNodeUJoint(const MassProperties& mProps_B,
                const Transform&      X_PF,
                const Transform&      X_BM,
                bool                  isReversed,
                UIndex&               nextUSlot,
                USquaredIndex&        nextUSqSlot,
                QIndex&               nextQSlot)
:   RigidBodyNodeSpec<2, false, noX_MB, noR_PF>(mProps_B,X_PF,X_BM,nextUSlot,nextUSqSlot,nextQSlot,
                         RigidBodyNode::QDotIsAlwaysTheSameAsU, RigidBodyNode::QuaternionIsNeverUsed,
                         isReversed)
{
    this->updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
}

void setQToFitRotationImpl(const SBStateDigest& sbs, const Rotation& R_FM,
                           Vector& q) const {
    // The only rotations this joint can handle are about Mx and My.
    // TODO: isn't there a better way to come up with "the rotation around x&y
    // that best approximates a rotation R"? Here we're just hoping that the
    // supplied rotation matrix can be decomposed into (x,y) rotations.
    const Vec2 angles12 = R_FM.convertRotationToBodyFixedXY();
    this->toQ(q) = angles12;
}

void setQToFitTranslationImpl(const SBStateDigest& sbs, const Vec3& p_FM,
                              Vector& q) const {
    // M and F frame origins are always coincident for this mobilizer so there
    // is no way to create a translation by rotating. So the only translation
    // we can represent is 0.
}

// We can only express angular velocity that can be produced with our
// generalized speeds which are Fx and My rotations rates. So we'll take the
// supplied angular velocity expressed in F, project it on Fx and use that as
// the first generalized speed. Then take whatever angular velocity is
// unaccounted for, express it in M, and project onto My and use that as the
// second generalized speed.
void setUToFitAngularVelocityImpl(const SBStateDigest& sbs, const Vector& q,
                                  const Vec3& w_FM, Vector& u) const {
    const Rotation R_FM =
        Rotation(BodyRotationSequence, this->fromQ(q)[0], XAxis, this->fromQ(q)[1], YAxis);  // body fixed 1-2 sequence
    const Vec3 wyz_FM_M = ~R_FM*Vec3(0,w_FM[1],w_FM[2]);
    this->toU(u) = Vec2(w_FM[0], wyz_FM_M[1]);
}

void setUToFitLinearVelocityImpl
    (const SBStateDigest& sbs, const Vector&, const Vec3& v_FM, Vector& u) const
{
    // M and F frame origins are always coincident for this mobilizer so there
    // is no way to create a linear velocity by rotating. So the only linear
    // velocity we can represent is 0.
}

enum {PoolSize=4}; // number of Reals
enum {CosQ=0, SinQ=2};
// We want space for cos(q01) and sin(q01).
int calcQPoolSize(const SBModelVars&) const
{   return PoolSize; }

// Precalculation of sin/cos costs around 100 flops.
void performQPrecalculations(const SBStateDigest& sbs,
                                const Real* q, int nq,
                                Real* qCache,  int nQCache,
                                Real* qErr,    int nQErr) const
{
    assert(q && nq==2 && qCache && nQCache==PoolSize && nQErr==0);
    Vec2::updAs(&qCache[CosQ]) = Vec2(std::cos(q[0]),std::cos(q[1]));
    Vec2::updAs(&qCache[SinQ]) = Vec2(std::sin(q[0]),std::sin(q[1]));
}

// TODO: should use precalculated sin/cos but doesn't.
void calcX_FM(const SBStateDigest& sbs,
                const Real* q,      int nq,
                const Real* qCache, int nQCache,
                Transform&  X_FM) const
{
    assert(q && nq==2 && qCache && nQCache==PoolSize);
    // TODO: use qCache
    // We're only updating the orientation here because a U-joint
    // can't translate. This is a body fixed X-Y sequence.
    X_FM.updR() = Rotation(BodyRotationSequence, q[0], XAxis, q[1], YAxis);
    X_FM.updP() = 0;
}

// The generalized speeds for this 2-dof rotational joint are the time
// derivatives of the body-fixed 1-2 rotation sequence defining the
// orientation. That is, the first speed is just a rotation rate about Fx. The
// second is a rotation rate about the current My, so we have to transform it
// into F to make H_FM uniformly expressed in F.
void calcAcrossJointVelocityJacobian(
    const SBStateDigest& sbs,
    HType&               H_FM) const
{
    // "upd" because we're realizing positions now
    const SBTreePositionCache& pc = sbs.updTreePositionCache();
    const Transform  X_F0M0 = this->findX_F0M0(pc);

    // Dropping the 0's here.
    const Rotation& R_FM = X_F0M0.R();

    H_FM(0) = SpatialVec(  Vec3(1,0,0) , Vec3(0) );
    H_FM(1) = SpatialVec(    R_FM.y()  , Vec3(0) );
}

// Since the second row of the Jacobian H_FM above is not constant in F,
// its time derivative is non zero. Here we use the fact that for
// a vector r_B_A fixed in a moving frame B but expressed in another frame A,
// its time derivative in A is the angular velocity of B in A crossed with
// the vector, i.e., d_A/dt r_B_A = w_AB % r_B_A.
void calcAcrossJointVelocityJacobianDot(
    const SBStateDigest& sbs,
    HType&               HDot_FM) const
{
    const SBTreePositionCache& pc = sbs.getTreePositionCache();
    // "upd" because we're realizing velocities now
    const SBTreeVelocityCache& vc = sbs.updTreeVelocityCache();

    const Transform  X_F0M0 = this->findX_F0M0(pc);

    // Dropping the 0's here.
    const Rotation& R_FM = X_F0M0.R();
    const Vec3      w_FM = this->find_w_F0M0(pc,vc); // angular velocity of M in F

    HDot_FM(0) = SpatialVec(     Vec3(0)     , Vec3(0) );
    HDot_FM(1) = SpatialVec( w_FM % R_FM.y() , Vec3(0) );
}

};




#endif // SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_UNIVERSAL_H_

