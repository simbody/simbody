#ifndef SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_BALL_H_
#define SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_BALL_H_

/* -------------------------------------------------------------------------- *
 *                              SimTK Simbody(tm)                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
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
 * Define the RigidBodyNode that implements a Ball mobilizer, also known as an
 * Orientation or Spherical joint.
 */

#include "SimbodyMatterSubsystemRep.h"
#include "RigidBodyNode.h"
#include "RigidBodyNodeSpec.h"

// ORIENTATION (BALL) //

// Ball joint. This provides three degrees of rotational freedom,  i.e.,
// unrestricted orientation of the body's M frame in the parent's F frame.
// The generalized coordinates are:
//   * 4 quaternions or 3 X-Y-Z (1-2-3) body fixed Euler angles (that is, 
//     fixed in M)
// and generalized speeds are:
//   * angular velocity w_FM as a vector expressed in the F (parent) frame.
// Thus rotational qdots have to be derived from the generalized speeds to
// be turned into either 4 quaternion derivatives or 3 Euler angle derivatives.
//
// NOTE: An XYZ Euler angle sequence has a singularity when the middle angle
// is at 90 or 270 degrees; quaternions are never singular.
class RBNodeBall : public RigidBodyNodeSpec<3> {
public:

virtual const char* type() { return "ball"; }

RBNodeBall(const MassProperties& mProps_B,
           const Transform&      X_PF,
           const Transform&      X_BM,
           bool                  isReversed,
           UIndex&               nextUSlot,
           USquaredIndex&        nextUSqSlot,
           QIndex&               nextQSlot)
  : RigidBodyNodeSpec<3>(mProps_B,X_PF,X_BM,nextUSlot,nextUSqSlot,nextQSlot,
                         QDotMayDifferFromU, QuaternionMayBeUsed, isReversed)
{
    updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
}

void setQToFitRotationImpl(const SBStateDigest& sbs, const Rotation& R_FM,
                          Vector& q) const 
{
    if (getUseEulerAngles(sbs.getModelVars()))
        toQ(q)    = R_FM.convertRotationToBodyFixedXYZ();
    else
        toQuat(q) = R_FM.convertRotationToQuaternion().asVec4();
}

void setQToFitTranslationImpl(const SBStateDigest& sbs, const Vec3& p_FM, Vector& q) const {
    // M and F frame origins are always coincident for this mobilizer so there is no
    // way to create a translation by rotating. So the only translation we can represent is 0.
}

void setUToFitAngularVelocityImpl(const SBStateDigest& sbs, const Vector&, const Vec3& w_FM,
                                 Vector& u) const
{
        toU(u) = w_FM; // relative angular velocity always used as generalized speeds
}

void setUToFitLinearVelocityImpl
   (const SBStateDigest& sbs, const Vector&, const Vec3& v_FM, Vector& u) const
{
    // M and F frame origins are always coincident for this mobilizer so there is no
    // way to create a linear velocity by rotating. So the only linear velocity
    // we can represent is 0.
}

// When we're using Euler angles, we're going to want to cache cos and sin for
// each angle, and also 1/cos of the middle angle will be handy to have around.
// When we're using quaternions, we'll cache 1/|q| for convenient normalizing
// of quaternions.
enum {AnglePoolSize=7, QuatPoolSize=1};
// These are indices into the pool: cos x,y,z sin x,y,z 1/cos(y)
enum {AngleCosQ=0, AngleSinQ=3, AngleOOCosQy=6};
enum {QuatOONorm=0};
int calcQPoolSize(const SBModelVars& mv) const
{   return getUseEulerAngles(mv) ? AnglePoolSize : QuatPoolSize; }

// This is expensive in Euler angle mode due to the three sin/cos computations
// required (about 150 flops). In quaternion mode it is much cheaper (about
// 25 flops).
void performQPrecalculations(const SBStateDigest& sbs,
                             const Real* q,      int nq,
                             Real*       qCache, int nQCache,
                             Real*       qErr,   int nQErr) const
{
    if (getUseEulerAngles(sbs.getModelVars())) {
        assert(q && nq==3 && qCache && nQCache==AnglePoolSize && nQErr==0);

        const Real cy = std::cos(q[1]);
        Vec3::updAs(&qCache[AngleCosQ]) =
            Vec3(std::cos(q[0]), cy, std::cos(q[2]));
        Vec3::updAs(&qCache[AngleSinQ]) =
            Vec3(std::sin(q[0]), std::sin(q[1]), std::sin(q[2]));
        qCache[AngleOOCosQy] = 1/cy; // trouble at 90 or 270 degrees
    } else {
        assert(q && nq==4 && qCache && nQCache==QuatPoolSize && 
               qErr && nQErr==1);

        const Real quatLen = Vec4::getAs(q).norm();
        qErr[0] = quatLen - Real(1);    // normalization error
        qCache[QuatOONorm] = 1/quatLen; // save for later
    }
}

// Because of the precalculations above we can calculate the cross-mobilizer
// transform in 18 flops (Euler angles) or 33 flops (quaternions).
void calcX_FM(const SBStateDigest& sbs,
              const Real* q,      int nq,
              const Real* qCache, int nQCache,
              Transform&  X_F0M0) const
{
    // Rotation.
    if (getUseEulerAngles(sbs.getModelVars())) {
        assert(q && nq==3 && qCache && nQCache==AnglePoolSize);
        X_F0M0.updR().setRotationToBodyFixedXYZ // 18 flops
           (Vec3::getAs(&qCache[AngleCosQ]), Vec3::getAs(&qCache[AngleSinQ]));
    } else {
        assert(q && nq==4 && qCache && nQCache==QuatPoolSize);
        // Must use a normalized quaternion to generate the rotation matrix.
        // Here we normalize with just 4 flops using precalculated 1/norm(q).
        const Quaternion quat(Vec4::getAs(q)*qCache[QuatOONorm], true); 
        X_F0M0.updR().setRotationFromQuaternion(quat); // 29 flops
    }

    // Translation.
    X_F0M0.updP() = 0; // This joint can't translate.
}

// Generalized speeds are the angular velocity expressed in F, so they
// cause rotations around F x,y,z axes respectively.
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


// Calculate qdot=N(q)*u. Precalculations make this fast in Euler angle
// mode (10 flops); quaternion mode is 27 flops.
// TODO: we're expecting that there are 4 qdots even in Euler angle mode
// so we have to zero out the last one in that case.
void calcQDot(const SBStateDigest& sbs, const Real* u, 
                             Real* qdot) const {
    assert(sbs.getStage() >= Stage::Position);
    assert(u && qdot);

    const SBModelVars& mv = sbs.getModelVars();

    // Generalized speeds are the same in either mode -- angular velocity
    // of F in M (or B in P), expressed in F (fixed on parent) frame.
    const Vec3& w_FM = Vec3::getAs(u);

    if (getUseEulerAngles(mv)) {
        // Euler angle mode (10 flops)
        qdot[3] = 0; // TODO: kludge, clear unused element

        const SBModelCache&        mc   = sbs.getModelCache();
        const SBTreePositionCache& pc   = sbs.getTreePositionCache();
        const Real*                pool = getQPool(mc, pc);

        Vec3::updAs(qdot) = Rotation::convertAngVelInParentToBodyXYZDot(
             Vec2::getAs(&pool[AngleCosQ]), Vec2::getAs(&pool[AngleSinQ]), //x,y
             pool[AngleOOCosQy], w_FM);
    } else {
        // Quaternion mode (27 flops)
        const Vec4& quat  = fromQuat(sbs.getQ()); // unnormalized
        Vec4::updAs(qdot) = Rotation::convertAngVelToQuaternionDot(quat, w_FM);
    }
}

// CAUTION: we do not zero the unused 4th element of q for Euler angles; it
// is up to the caller to do that if it is necessary.
// Compute out_q = N * in_u
//   or    out_u = in_q * N
// The u quantity is a vector v_F in the parent frame (e.g. angular 
// velocity of child in parent), while the q quantity is the derivative 
// of a quaternion representing R_FM, i.e. orientation of child in parent.
// Cost: Euler angle mode - 10 flops, Quaternion mode - 27 flops.
void multiplyByN(const SBStateDigest& sbs, bool matrixOnRight, 
                 const Real* in, Real* out) const
{
    assert(sbs.getStage() >= Stage::Position);
    assert(in && out);

    const SBModelVars& mv = sbs.getModelVars();

    if (getUseEulerAngles(mv)) {
        // Euler angle mode (10 flops).

        const SBModelCache&        mc   = sbs.getModelCache();
        const SBTreePositionCache& pc   = sbs.getTreePositionCache();
        const Real*                pool = getQPool(mc, pc);

        // Aliases.
        const Vec2& cosxy  = Vec2::getAs(&pool[AngleCosQ]); // only need x,y
        const Vec2& sinxy  = Vec2::getAs(&pool[AngleSinQ]);
        const Real& oocosy = pool[AngleOOCosQy];
        const Vec3& in_v   = Vec3::getAs(in);
        Vec3&       out_v  = Vec3::updAs(out);

        if (matrixOnRight) { // out_u = in_q * N = ~N * in_q (9 flops)
            // Transpose on left is same as untransposed on right.
            out_v = Rotation::multiplyByBodyXYZ_NT_P(cosxy, sinxy, oocosy, in_v);
        } else {             // out_q = N*in_u (10 flops)
            out_v = Rotation::multiplyByBodyXYZ_N_P(cosxy, sinxy, oocosy, in_v);
        }
    } else {
        // Quaternion mode (27 flops).
        const Vec4& quat  = fromQuat(sbs.getQ()); // unnormalized
        const Mat43 N_F = // This method returns N for the parent frame
            Rotation::calcUnnormalizedNForQuaternion(quat); // 7 flops
        if (matrixOnRight) 
                Row3::updAs(out) = Row4::getAs(in) * N_F; // 20 flops
        else    Vec4::updAs(out) = N_F * Vec3::getAs(in);
    }
}

// Compute out_u = inv(N) * in_q
//   or    out_q = in_u * inv(N)
// Cost: Euler angle mode ~ 10 flops, Quaternion mode ~ 27 flops.
void multiplyByNInv(const SBStateDigest& sbs, bool matrixOnRight,
                    const Real* in, Real* out) const
{
    assert(sbs.getStage() >= Stage::Position);
    assert(in && out);

    const SBModelVars& mv = sbs.getModelVars();

    if (getUseEulerAngles(mv)) {
        // Euler angle mode (10 flops).
        const SBModelCache&        mc   = sbs.getModelCache();
        const SBTreePositionCache& pc   = sbs.getTreePositionCache();
        const Real*                pool = getQPool(mc, pc);

        // Aliases.
        const Vec2& cosxy  = Vec2::getAs(&pool[AngleCosQ]); // only need x,y
        const Vec2& sinxy  = Vec2::getAs(&pool[AngleSinQ]);
        const Vec3& in_v   = Vec3::getAs(in);
        Vec3&       out_v  = Vec3::updAs(out);

        if (matrixOnRight) { // out_q = in_u * inv(N) = ~inv(N) * in_u  (10 flops)
            // transpose on left is same as untransposed on right
            out_v = Rotation::multiplyByBodyXYZ_NInvT_P(cosxy, sinxy, in_v);
        } else {             // out_u = inv(N) * in_q                   (9 flops)
            out_v = Rotation::multiplyByBodyXYZ_NInv_P(cosxy, sinxy, in_v);
        }
    } else {
        // Quaternion mode (27 flops).
        const Vec4& quat  = fromQuat(sbs.getQ()); // unnormalized
        const Mat34 NInv_F = // Returns NInv for parent frame.
            Rotation::calcUnnormalizedNInvForQuaternion(quat);  // 7 flops
        if (matrixOnRight) 
                Row4::updAs(out) = Row3::getAs(in) * NInv_F; // 20 flops
        else    Vec3::updAs(out) = NInv_F * Vec4::getAs(in); // 21 flops
    }
}

// Compute out_q = NDot * in_u
//   or    out_u = in_q * NDot
// Note that it is much cheaper to calculate
// qdotdot directly than to do it using N and NDot; see below.
// Cost: Euler angle mode - 36 flops, Quaternion mode - 27 flops.
void multiplyByNDot(const SBStateDigest& sbs, bool matrixOnRight,
                    const Real* in, Real* out) const
{
    assert(sbs.getStage() > Stage::Velocity);
    assert(in && out);

    const SBModelVars& mv = sbs.getModelVars();

    if (getUseEulerAngles(mv)) {
        // Euler angle mode (36 flops).
        const SBModelCache&        mc   = sbs.getModelCache();
        const SBTreePositionCache& pc   = sbs.getTreePositionCache();
        const Real*                pool = getQPool(mc, pc);

        // Aliases.
        const Vec2& cosxy  = Vec2::getAs(&pool[AngleCosQ]); // only need x,y
        const Vec2& sinxy  = Vec2::getAs(&pool[AngleSinQ]);
        const Real& oocosy = pool[AngleOOCosQy];
        const Vec3& qdot   = fromQ(sbs.getQDot());

        // We don't have a nice multiply-by routine here so just get the NDot
        // matrix and use it (21 flops).
        const Mat33 NDot_F = 
            Rotation::calcNDotForBodyXYZInParentFrame(cosxy, sinxy, oocosy, qdot);

        if (matrixOnRight) {    // out_u = in_q * NDot (15 flops)
            Row3::updAs(out) = Row3::getAs(in) * NDot_F;
        } else {                // out_q = NDot * in_u (15 flops)
            Vec3::updAs(out) = NDot_F * Vec3::getAs(in);
        }
    } else {
        // Quaternion mode (27 flops).
        const Vec4& qdot  = fromQuat(sbs.getQDot()); // unnormalized
        const Mat43 NDot_F = // This method returns NDot for the parent frame
            Rotation::calcUnnormalizedNDotForQuaternion(qdot); // 7 flops
        if (matrixOnRight) 
                Row3::updAs(out) = Row4::getAs(in) * NDot_F;  // 21 flops
        else    Vec4::updAs(out) = NDot_F * Vec3::getAs(in);  // 20 flops
    }
}

// Calculate qdotdot from udot as qdotdot = N(q)*udot + NDot(q,qdot)*u, but
// faster. Here we assume that qdot=N*u has already been calculated, which
// is very helpful in Euler angle mode. In either mode it is a very bad
// idea to calculate this using explicit N and NDot matrices; we can save
// a lot of time by using more specialized methods.
// Cost: Euler angle mode - 22 flops, Quaternion mode - 41 flops.
void calcQDotDot(const SBStateDigest& sbs, const Real* udot, 
                                   Real* qdotdot) const {
    assert(sbs.getStage() > Stage::Velocity);
    assert(udot && qdotdot);

    const SBModelVars& mv = sbs.getModelVars();

    // Angular acceleration is the same in either mode.
    const Vec3& b_FM = Vec3::getAs(udot);   // a.k.a. w_FM_dot

    if (getUseEulerAngles(mv)) {
        // Euler angle mode (22 flops).
        Vec4::updAs(qdotdot) = 0; // TODO: kludge, clear unused element

        const SBModelCache&        mc   = sbs.getModelCache();
        const SBTreePositionCache& pc   = sbs.getTreePositionCache();
        const Real*                pool = getQPool(mc, pc);

        // Aliases.
        const Vec2& cosxy     = Vec2::getAs(&pool[AngleCosQ]); // only need x,y
        const Vec2& sinxy     = Vec2::getAs(&pool[AngleSinQ]);
        const Real& oocosy    = pool[AngleOOCosQy];
        const Vec3& qdot      = fromQ(sbs.getQDot());
        Vec3&       qdotdot_v = Vec3::updAs(qdotdot);

        // 22 flops.
        qdotdot_v = Rotation::convertAngAccInParentToBodyXYZDotDot
                                            (cosxy, sinxy, oocosy, qdot, b_FM);
    } else {
        // Quaternion mode (41 flops).
        const Vec4& quat  = fromQuat(sbs.getQ()); // unnormalized
        const Vec3& w_FM = fromU(sbs.getU());
        Vec4::updAs(qdotdot) = 
            Rotation::convertAngVelDotToQuaternionDotDot(quat,w_FM,b_FM);
    }
}

int getMaxNQ()              const {return 4;}
int getNQInUse(const SBModelVars& mv) const {
    return getUseEulerAngles(mv) ? 3 : 4;
} 
bool isUsingQuaternion(const SBStateDigest& sbs, MobilizerQIndex& startOfQuaternion) const {
    if (getUseEulerAngles(sbs.getModelVars())) {
        startOfQuaternion.invalidate(); 
        return false;
    }
    startOfQuaternion = MobilizerQIndex(0); // quaternion comes first
    return true;
}

void setMobilizerDefaultPositionValues(
    const SBModelVars& mv,
    Vector&            q) const 
{
    if (getUseEulerAngles(mv)) {
        //TODO: kludge
        toQuat(q) = Vec4(0); // clear unused element
        toQ(q) = 0.;
    }
    else toQuat(q) = Vec4(1.,0.,0.,0.);
}

bool enforceQuaternionConstraints(
    const SBStateDigest& sbs,
    Vector&             q,
    Vector&             qErrest) const 
{
    if (getUseEulerAngles(sbs.getModelVars())) 
        return false;   // no change

    Vec4& quat = toQuat(q);
    quat = quat / quat.norm();

    if (qErrest.size()) {
        Vec4& qerr = toQuat(qErrest);
        qerr -= dot(qerr,quat) * quat;
    }

    return true;
}

void convertToEulerAngles(const Vector& inputQ, Vector& outputQ) const {
    toQuat(outputQ) = Vec4(0); // clear unused element
    toQ(outputQ) = Rotation(Quaternion(fromQuat(inputQ))).convertRotationToBodyFixedXYZ();
}

void convertToQuaternions(const Vector& inputQ, Vector& outputQ) const {
    Rotation rot;
    rot.setRotationToBodyFixedXYZ(fromQ(inputQ));
    toQuat(outputQ) = rot.convertRotationToQuaternion().asVec4();
}


};


#endif // SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_BALL_H_

