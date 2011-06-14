#ifndef SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_GIMBAL_H_
#define SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_GIMBAL_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-11 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *    Charles Schwieters (NIH): wrote the public domain IVM code from which   *
 *                              this was derived.                             *
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
 * Define the RigidBodyNode that implements a Gimbal mobilizer.
 */

#include "SimbodyMatterSubsystemRep.h"
#include "RigidBodyNode.h"
#include "RigidBodyNodeSpec.h"


    // GIMBAL //

// Gimbal joint. This provides three degrees of rotational freedom,  i.e.,
// unrestricted orientation of the body's M frame in the parent's F frame.
// The generalized coordinates are:
//   * 3 X-Y-Z (1-2-3) body fixed Euler angles (that is, fixed in M)
// and generalized speeds are:
//   * angular velocity w_FM as a vector expressed in the F (parent) frame.
// Thus rotational qdots have to be derived from the generalized speeds to
// be turned into 3 Euler angle derivatives.
//
// NOTE: This joint has a singularity when the middle angle is near 90 degrees.
// In most cases you should use a Ball joint instead, which by default uses
// a quaternion as its generalized coordinates to avoid this singularity. A
// modeling option allows the Ball to be switched to use Euler angles when 
// convenient.

template<bool noX_MB, bool noR_PF>
class RBNodeGimbal : public RigidBodyNodeSpec<3, false, noX_MB, noR_PF> {
public:

typedef typename RigidBodyNodeSpec<3, false, noX_MB, noR_PF>::HType HType;
virtual const char* type() { return "gimbal"; }

RBNodeGimbal( const MassProperties& mProps_B,
              const Transform&      X_PF,
              const Transform&      X_BM,
              bool                  isReversed,
              UIndex&               nextUSlot,
              USquaredIndex&        nextUSqSlot,
              QIndex&               nextQSlot)
  : RigidBodyNodeSpec<3, false, noX_MB, noR_PF>(mProps_B,X_PF,X_BM,nextUSlot,nextUSqSlot,nextQSlot,
                         RigidBodyNode::QDotMayDifferFromU, RigidBodyNode::QuaternionIsNeverUsed, isReversed)
{
    this->updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
}

void setQToFitRotationImpl(const SBStateDigest& sbs, const Rotation& R_FM,
                           Vector& q) const {
    this->toQ(q) = R_FM.convertRotationToBodyFixedXYZ();
}

void setQToFitTranslationImpl(const SBStateDigest& sbs, const Vec3& p_FM, 
                              Vector& q) const {
    // M and F frame origins are always coincident for this mobilizer so 
    // there is no way to create a translation by rotating. So the only 
    // translation we can represent is 0.
}

void setUToFitAngularVelocityImpl
   (const SBStateDigest& sbs, const Vector&, const Vec3& w_FM,
    Vector& u) const 
{
    this->toU(u) = w_FM; // relative ang. vel. always used as generalized speeds
}

void setUToFitLinearVelocityImpl
   (const SBStateDigest& sbs, const Vector&, const Vec3& v_FM, 
    Vector& u) const
{
    // M and F frame origins are always coincident for this mobilizer so 
    // there is no way to create a linear velocity by rotating. So the 
    // only linear velocity we can represent is 0.
}

// We want to cache cos and sin for each angle, and also 1/cos of the middle 
// angle will be handy to have around.
enum {PoolSize=7};
// cos x,y,z sin x,y,z 1/cos(y)
enum {CosQ=0, SinQ=3, OOCosQy=6};
int calcQPoolSize(const SBModelVars& mv) const
{   return PoolSize; }

// This is expensive since we have three sin/cos computations and a divide
// to do, approx. 150 flops. But we hope to re-use these calculations several
// times before we're done with a complete realization.
void performQPrecalculations(const SBStateDigest& sbs,
                             const Real* q,      int nq,
                             Real*       qCache, int nQCache,
                             Real*       qErr,   int nQErr) const
{
    assert(q && nq==3 && qCache && nQCache==PoolSize && nQErr==0);

    const Real cy = std::cos(q[1]);
    Vec3::updAs(&qCache[CosQ]) =
        Vec3(std::cos(q[0]), cy, std::cos(q[2]));
    Vec3::updAs(&qCache[SinQ]) =
        Vec3(std::sin(q[0]), std::sin(q[1]), std::sin(q[2]));
    qCache[OOCosQy] = 1/cy; // trouble at 90 or 270 degrees
}

// Because of the precalculations we can calculate the cross-mobilizer
// transform in only 18 flops.
void calcX_FM(const SBStateDigest& sbs,
              const Real* q,      int nq,
              const Real* qCache, int nQCache,
              Transform&  X_F0M0) const
{
    assert(q && nq==3 && qCache && nQCache==PoolSize);

    X_F0M0.updR().setRotationToBodyFixedXYZ // 18 flops
        (Vec3::getAs(&qCache[CosQ]), Vec3::getAs(&qCache[SinQ]));
    X_F0M0.updP() = 0.; // This joint can't translate.
}

// Generalized speeds are the angular velocity expressed in F (the parent), 
// so they cause rotations around F x,y,z axes respectively.
void calcAcrossJointVelocityJacobian(
    const SBStateDigest& sbs,
    HType&               H_FM) const
{
    H_FM(0) = SpatialVec( Vec3(1,0,0), Vec3(0) );
    H_FM(1) = SpatialVec( Vec3(0,1,0), Vec3(0) );
    H_FM(2) = SpatialVec( Vec3(0,0,1), Vec3(0) );
}

// The derivative of constant matrix H_FM is not too interesting.
void calcAcrossJointVelocityJacobianDot(
    const SBStateDigest& sbs,
    HType&               HDot_FM) const
{
    HDot_FM(0) = SpatialVec( Vec3(0), Vec3(0) );
    HDot_FM(1) = SpatialVec( Vec3(0), Vec3(0) );
    HDot_FM(2) = SpatialVec( Vec3(0), Vec3(0) );
}


// Here we again use the pooled trig evaluations to calculate qdot=N(q)*u
// on the cheap (10 flops).
void calcQDot(const SBStateDigest& sbs, const Real* u, 
                             Real* qdot) const {
    assert(sbs.getStage() >= Stage::Position);
    assert(u && qdot);

    const SBModelCache&        mc   = sbs.getModelCache();
    const SBTreePositionCache& pc   = sbs.getTreePositionCache();
    const Real*                pool = this->getQPool(mc, pc);

    const Vec3& w_FM = Vec3::getAs(u);
    Vec3::updAs(qdot) = Rotation::convertAngVelInParentToBodyXYZDot(
         Vec2::getAs(&pool[CosQ]), Vec2::getAs(&pool[SinQ]), // only x,y
         pool[OOCosQy], w_FM);
}


// Compute out_q = N * in_u
//   or    out_u = in_q * N
// About 10 flops.
void multiplyByN(const SBStateDigest& sbs, bool matrixOnRight, 
                 const Real* in, Real* out) const
{
    assert(sbs.getStage() > Stage::Position);
    assert(in && out);

    const SBModelCache&        mc   = sbs.getModelCache();
    const SBTreePositionCache& pc   = sbs.getTreePositionCache();
    const Real*                pool = this->getQPool(mc, pc);

    // Aliases.
    const Vec2& cosxy  = Vec2::getAs(&pool[CosQ]); // only need x,y angles
    const Vec2& sinxy  = Vec2::getAs(&pool[SinQ]);
    const Real& oocosy = pool[OOCosQy];
    const Vec3& in_v   = Vec3::getAs(in);
    Vec3&       out_v  = Vec3::updAs(out);

    if (matrixOnRight) { // out_u = in_q * N = ~N * in_q (9 flops)
        // Transpose on left is same as untransposed on right.
        out_v = Rotation::multiplyByBodyXYZ_NT_P(cosxy, sinxy, oocosy, in_v);
    } else {             // out_q = N*in_u (10 flops)
        out_v = Rotation::multiplyByBodyXYZ_N_P(cosxy, sinxy, oocosy, in_v);
    }
}

// Compute out_u = inv(N) * in_q
//   or    out_q = in_u * inv(N)
// About 10 flops.
void multiplyByNInv(const SBStateDigest& sbs, bool matrixOnRight, 
                    const Real* in, Real* out) const
{
    assert(sbs.getStage() > Stage::Position);
    assert(in && out);

    const SBModelCache&        mc   = sbs.getModelCache();
    const SBTreePositionCache& pc   = sbs.getTreePositionCache();
    const Real*                pool = this->getQPool(mc, pc);

    // Aliases.
    const Vec2& cosxy  = Vec2::getAs(&pool[CosQ]); // only need x,y angles
    const Vec2& sinxy  = Vec2::getAs(&pool[SinQ]);
    const Vec3& in_v   = Vec3::getAs(in);
    Vec3&       out_v  = Vec3::updAs(out);

    if (matrixOnRight) { // out_q = in_u * inv(N) = ~inv(N) * in_u  (10 flops)
        // transpose on left is same as untransposed on right
        out_v = Rotation::multiplyByBodyXYZ_NInvT_P(cosxy, sinxy, in_v);
    } else {             // out_u = inv(N) * in_q                   (9 flops)
        out_v = Rotation::multiplyByBodyXYZ_NInv_P(cosxy, sinxy, in_v);
    }
}


// Compute out_q = NDot * in_u
//   or    out_u = in_q * NDot
// Either way, 36 flops. Note that it is much cheaper to calculate
// qdotdot directly than to do it using N and NDot; see below.
void multiplyByNDot(const SBStateDigest& sbs, bool matrixOnRight, 
                    const Real* in, Real* out) const
{
    assert(sbs.getStage() > Stage::Velocity);
    assert(in && out);

    const SBModelCache&        mc   = sbs.getModelCache();
    const SBTreePositionCache& pc   = sbs.getTreePositionCache();
    const Real*                pool = this->getQPool(mc, pc);

    // Aliases.
    const Vec2& cosxy  = Vec2::getAs(&pool[CosQ]); // only need x,y angles
    const Vec2& sinxy  = Vec2::getAs(&pool[SinQ]);
    const Real& oocosy = pool[OOCosQy];
    const Vec3& qdot   = this->fromQ(sbs.getQDot());

    // We don't have a nice multiply-by routine here so just get the NDot
    // matrix and use it (21 flops).
    const Mat33 NDot_F = 
        Rotation::calcNDotForBodyXYZInParentFrame(cosxy, sinxy, oocosy, qdot);

    if (matrixOnRight) {    // out_u = in_q * NDot (15 flops)
        Row3::updAs(out) = Row3::getAs(in) * NDot_F;
    } else {                // out_q = NDot * in_u (15 flops)
        Vec3::updAs(out) = NDot_F * Vec3::getAs(in);
    }
}

// Calculate qdotdot from udot as qdotdot = N(q)*udot + NDot(q,qdot)*u, but
// faster. Here we assume that qdot=N*u has already been calculated and this 
// allows us to calculate qdotdot in only 22 flops, which is faster than
// calculating NDot*v alone using the multiplyByNDot operator.
void calcQDotDot(const SBStateDigest& sbs, 
                                   const Real* udot, Real* qdotdot) const {
    assert(sbs.getStage() > Stage::Velocity);
    assert(udot && qdotdot);

    const SBModelCache&        mc   = sbs.getModelCache();
    const SBTreePositionCache& pc   = sbs.getTreePositionCache();
    const Real*                pool = this->getQPool(mc, pc);

    // Aliases.
    const Vec2& cosxy  = Vec2::getAs(&pool[CosQ]); // only need x,y angles
    const Vec2& sinxy  = Vec2::getAs(&pool[SinQ]);
    const Real& oocosy = pool[OOCosQy];
    const Vec3& qdot = this->fromQ(sbs.getQDot());

    const Vec3& b_FM      = Vec3::getAs(udot); // = w_FM_dot (angular accel.)
    Vec3&       qdotdot_v = Vec3::updAs(qdotdot);

    // 22 flops.
    qdotdot_v = Rotation::convertAngAccInParentToBodyXYZDotDot
                                        (cosxy, sinxy, oocosy, qdot, b_FM);
}

};



#endif // SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_GIMBAL_H_

