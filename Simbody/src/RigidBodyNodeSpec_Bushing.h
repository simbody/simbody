#ifndef SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_BUSHING_H_
#define SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_BUSHING_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2012 Stanford University and the Authors.           *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
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
 * Define the RigidBodyNode that implements a Bushing mobilizer.
 */

#include "SimbodyMatterSubsystemRep.h"
#include "RigidBodyNode.h"
#include "RigidBodyNodeSpec.h"


    // BUSHING //

// Bushing joint. This provides a full six degrees of freedom but with
// orientation represented using Euler angles and thus not suited to
// arbitrarily large rotations.
// The generalized coordinates q are (in this order):
//   * 3 x-y-z (1-2-3) body fixed Euler angles (that is, fixed in M)
//   * 3 measure numbers of position vector p_FM
// and generalized speeds u are:
//   * u = qdot, that is, the x-y-z body fixed Euler angle derivatives
//     followed by the velocity vector v_FM
//
// NOTE: This joint has a singularity when the middle angle (q[1]) is near
// +/-90 degrees. In most cases you should use a Free joint instead, which by
// default uses a quaternion as its generalized coordinates to avoid this
// singularity. A modeling option allows the Free to be switched to use Euler
// angles when necessary.

template<bool noX_MB, bool noR_PF>
class RBNodeBushing : public RigidBodyNodeSpec<6, false, noX_MB, noR_PF> {
public:

typedef typename RigidBodyNodeSpec<6, false, noX_MB, noR_PF>::HType HType;
virtual const char* type() { return "bushing"; }

RBNodeBushing(const MassProperties& mProps_B,
              const Transform&      X_PF,
              const Transform&      X_BM,
              bool                  isReversed,
              UIndex&               nextUSlot,
              USquaredIndex&        nextUSqSlot,
              QIndex&               nextQSlot)
:   RigidBodyNodeSpec<6, false, noX_MB, noR_PF>
       (mProps_B,X_PF,X_BM,
        nextUSlot,nextUSqSlot,nextQSlot,
        RigidBodyNode::QDotIsAlwaysTheSameAsU,
        RigidBodyNode::QuaternionIsNeverUsed,
        isReversed)
{
    this->updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
}

    // Implementations of virtual methods.

// This has a default implementation but it rotates first then translates,
// which works fine for the normal bushing but produces wrong behavior when
// the mobilizer is reversed.
void setQToFitTransformImpl(const SBStateDigest& sbs, const Transform& X_FM,
                            Vector& q) const override
{
    setQToFitTranslationImpl(sbs, X_FM.p(), q); // see below
    setQToFitRotationImpl(sbs, X_FM.R(), q);
}

void setQToFitRotationImpl(const SBStateDigest& sbs, const Rotation& R_FM,
                           Vector& q) const override {
    this->toQVec3(q,0) = R_FM.convertRotationToBodyFixedXYZ();
}

void setQToFitTranslationImpl(const SBStateDigest& sbs, const Vec3& p_FM,
                              Vector& q) const override {
    this->toQVec3(q,3) = p_FM; // skip the 3 Euler angles
}

// Given the angular velocity of M in F, expressed in F, compute the Euler
// angle derivatives qdot that would produce that angular velocity, and
// return u=qdot.
void setUToFitAngularVelocityImpl
   (const SBStateDigest& sbs, const Vector& q, const Vec3& w_FM,
    Vector& u) const override
{
    const Vec2 cosxy(std::cos(q[0]), std::cos(q[1]));
    const Vec2 sinxy(std::sin(q[0]), std::sin(q[1]));
    const Real oocosy = 1 / cosxy[1];
    const Vec3 qdot =
        Rotation::convertAngVelInParentToBodyXYZDot(cosxy,sinxy,oocosy,w_FM);
    this->toUVec3(u,0) = qdot;
}

void setUToFitLinearVelocityImpl
   (const SBStateDigest& sbs, const Vector&, const Vec3& v_FM,
    Vector& u) const override
{
    this->toUVec3(u,3) = v_FM;
}

// We want to cache cos and sin for each angle, and also 1/cos of the middle
// angle will be handy to have around.
enum {PoolSize=7};
// cos x,y,z sin x,y,z 1/cos(y)
enum {CosQ=0, SinQ=3, OOCosQy=6};
int calcQPoolSize(const SBModelVars& mv) const override
{   return PoolSize; }

// This is expensive since we have three sin/cos computations and a divide
// to do, approx. 150 flops. But we hope to re-use these calculations several
// times before we're done with a complete realization.
void performQPrecalculations(const SBStateDigest& sbs,
                             const Real* q,      int nq,
                             Real*       qCache, int nQCache,
                             Real*       qErr,   int nQErr) const override
{
    assert(q && nq==6 && qCache && nQCache==PoolSize && nQErr==0);

    const Real cy = std::cos(q[1]);
    Vec3::updAs(&qCache[CosQ]) =
        Vec3(std::cos(q[0]), cy, std::cos(q[2]));
    Vec3::updAs(&qCache[SinQ]) =
        Vec3(std::sin(q[0]), std::sin(q[1]), std::sin(q[2]));
    qCache[OOCosQy] = 1/cy; // trouble at 90 or 270 (-90) degrees
}

// Because of the precalculations we can calculate the cross-mobilizer
// transform in only 18 flops.
void calcX_FM(const SBStateDigest& sbs,
              const Real* q,      int nq,
              const Real* qCache, int nQCache,
              Transform&  X_F0M0) const override
{
    assert(q && nq==6 && qCache && nQCache==PoolSize);

    X_F0M0.updR().setRotationToBodyFixedXYZ // 18 flops
        (Vec3::getAs(&qCache[CosQ]), Vec3::getAs(&qCache[SinQ]));
    X_F0M0.updP() = Vec3::getAs(&q[3]); // a012 x y z
}

// Generalized speeds are the Euler angle derivatives and the velocity
// vector v_FM. The H_FM matrix maps those into the angular velocity of
// M in F, expressed in F, and the linear velocity
//  [w_FM]   [Hw_FM]
//  [    ] = [     ] * qdot      (where v_FM=qdot(3,3))
//  [v_FM]   [ 0 I ]
// Using the precalculations this requires only 3 flops, although there is
// a fair bit of poking around in memory required.
void calcAcrossJointVelocityJacobian(
    const SBStateDigest& sbs,
    HType&               H_FM) const override
{
    const SBModelCache&        mc = sbs.getModelCache();
    // Use "upd" here because we're realizing positions now.
    const SBTreePositionCache& pc = sbs.updTreePositionCache();
    const Real*                pool = this->getQPool(mc, pc);

    const Real c0 = pool[CosQ], c1 = pool[CosQ+1];
    const Real s0 = pool[SinQ], s1 = pool[SinQ+1];

    // Fill in columns of H_FM. See Rotation::calcNInvForBodyXYZInParentFrame().
    H_FM(0) = SpatialVec(Vec3(1,     0,    0),      Vec3(0)  ); // rotations
    H_FM(1) = SpatialVec(Vec3(0,    c0,    s0),     Vec3(0)  );
    H_FM(2) = SpatialVec(Vec3(s1, -s0*c1, c0*c1),   Vec3(0)  );
    H_FM(3) = SpatialVec(         Vec3(0),        Vec3(1,0,0)); // translations
    H_FM(4) = SpatialVec(         Vec3(0),        Vec3(0,1,0));
    H_FM(5) = SpatialVec(         Vec3(0),        Vec3(0,0,1));
}

// Differentiate H_FM to get HDot_FM. Note that this depends on qdot:
//    d/dt cos(q0) = -sin(q0)*qdot0, etc.
void calcAcrossJointVelocityJacobianDot(
    const SBStateDigest& sbs,
    HType&               HDot_FM) const override
{
    const SBModelCache&        mc = sbs.getModelCache();
    const SBTreePositionCache& pc = sbs.getTreePositionCache();

    const Real* pool = this->getQPool(mc, pc);
    const Real c0 = pool[CosQ], c1 = pool[CosQ+1];
    const Real s0 = pool[SinQ], s1 = pool[SinQ+1];

    // Use "upd" here because we're realizing velocities now.
    const Vec6& qdot = this->fromQ(sbs.updQDot());
    const Real qd0 = qdot[0], qd1 = qdot[1];

    const Real dc0 = -s0*qd0, dc1 = -s1*qd1; // derivatives of c0,c1,s0,s1
    const Real ds0 =  c0*qd0, ds1 =  c1*qd1;

    // Compare with H_FM above.
    HDot_FM(0) = SpatialVec(Vec3(0,         0,              0),       Vec3(0));
    HDot_FM(1) = SpatialVec(Vec3(0,        dc0,            ds0),      Vec3(0));
    HDot_FM(2) = SpatialVec(Vec3(ds1, -ds0*c1-s0*dc1, dc0*c1+c0*dc1), Vec3(0));
    HDot_FM(3) = SpatialVec(             Vec3(0),                     Vec3(0));
    HDot_FM(4) = SpatialVec(             Vec3(0),                     Vec3(0));
    HDot_FM(5) = SpatialVec(             Vec3(0),                     Vec3(0));
}

// Can use default for calcQDot, multiplyByN, etc., since qdot==u for Bushing
// mobilizer.

};



#endif // SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_BUSHING_H_

