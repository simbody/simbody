#ifndef SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_BALL_H_
#define SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_BALL_H_

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
//   * 4 quaternions or 3 1-2-3 body fixed Euler angles (that is, fixed in M)
// and generalized speeds are:
//   * angular velocity w_FM as a vector expressed in the F frame.
// Thus rotational qdots have to be derived from the generalized speeds to
// be turned into either 4 quaternion derivatives or 3 Euler angle derivatives.
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

// This is required for all mobilizers.
bool isUsingAngles(const SBStateDigest& sbs, 
                   MobilizerQIndex&     startOfAngles, 
                   int&                 nAngles) const 
{
    // Ball joint has three angular coordinates when Euler angles are being used, 
    // none when quaternions are being used.
    if (!getUseEulerAngles(sbs.getModelVars())) 
    {   startOfAngles.invalidate(); nAngles=0; return false; }

    startOfAngles = MobilizerQIndex(0);
    nAngles = 3;
    return true;
}

// Precalculate sines and cosines.
void calcJointSinCosQNorm(
    const SBModelVars&      mv,
    const SBModelCache&     mc,
    const SBInstanceCache&  ic,
    const Vector&           q, 
    Vector&                 sine, 
    Vector&                 cosine, 
    Vector&                 qErr,
    Vector&                 qnorm) const
{
    const SBModelCache::PerMobilizedBodyModelInfo& 
        bInfo = mc.getMobilizedBodyModelInfo(nodeNum);

    if (getUseEulerAngles(mv)) {
        const Vec3& a = fromQ(q); // angular coordinates
        toQ(sine)   = Vec3(std::sin(a[0]), std::sin(a[1]), std::sin(a[2]));
        toQ(cosine) = Vec3(std::cos(a[0]), std::cos(a[1]), std::cos(a[2]));
        // no quaternions
    } else {
        // no angles
        const Vec4& quat = fromQuat(q); // unnormalized quaternion from state
        const Real  quatLen = quat.norm();
        assert(bInfo.hasQuaternionInUse && bInfo.quaternionPoolIndex.isValid());
        qErr[ic.firstQuaternionQErrSlot+bInfo.quaternionPoolIndex] = 
            quatLen - Real(1);
        toQuat(qnorm) = quat / quatLen;
    }
}

// Calculate X_FM.
void calcAcrossJointTransform(
    const SBStateDigest& sbs,
    const Vector&        q,
    Transform&           X_FM) const
{
    const SBModelVars& mv = sbs.getModelVars();
    X_FM.updP() = 0.; // This joint can't translate.
    if (getUseEulerAngles(mv))
        X_FM.updR().setRotationToBodyFixedXYZ( fromQ(q) );
    else {
        // TODO: should use qnorm pool; this performs costly normalize
        X_FM.updR().setRotationFromQuaternion( Quaternion(fromQuat(q)) );
    }
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

void calcAcrossJointVelocityJacobianDot(
    const SBStateDigest& sbs,
    HType&               HDot_FM) const
{
    HDot_FM(0) = SpatialVec( Vec3(0), Vec3(0) );
    HDot_FM(1) = SpatialVec( Vec3(0), Vec3(0) );
    HDot_FM(2) = SpatialVec( Vec3(0), Vec3(0) );
}

void calcQDot(
    const SBStateDigest&   sbs,
    const Vector&          u, 
    Vector&                qdot) const 
{
    const SBModelVars& mv = sbs.getModelVars();
    const SBTreePositionCache& pc = sbs.getTreePositionCache();
    const Vec3& w_FM = fromU(u); // angular velocity of M in F 
    if (getUseEulerAngles(mv)) {
        toQuat(qdot) = Vec4(0); // TODO: kludge, clear unused element
        const Rotation& R_FM = getX_FM(pc).R();
        toQ(qdot) = Rotation::convertAngVelInBodyFrameToBodyXYZDot(fromQ(sbs.getQ()),
                                    ~R_FM*w_FM); // need w in *body*, not parent
    } else
        toQuat(qdot) = Rotation::convertAngVelToQuaternionDot(fromQuat(sbs.getQ()),w_FM);
}

void calcLocalQDotFromLocalU(const SBStateDigest& sbs, const Real* u, Real* qdot) const {
    assert(sbs.getStage() >= Stage::Position);
    assert(u && qdot);

    const SBModelVars&          mv   = sbs.getModelVars();
    const SBTreePositionCache&  pc   = sbs.getTreePositionCache();
    const Vector&               allQ = sbs.getQ();

    const Vec3&            w_FM = Vec3::getAs(u);

    if (getUseEulerAngles(mv)) {
        Vec4::updAs(qdot) = 0; // TODO: kludge, clear unused element
        const Rotation& R_FM = getX_FM(pc).R();
        Vec3::updAs(qdot) = Rotation::convertAngVelInBodyFrameToBodyXYZDot(fromQ(allQ),
                                         ~R_FM*w_FM); // need w in *body*, not parent
    } else
        Vec4::updAs(qdot) = Rotation::convertAngVelToQuaternionDot(fromQuat(allQ),w_FM);
}

// CAUTION: we do not zero the unused 4th element of q for Euler angles; it
// is up to the caller to do that if it is necessary.
// Compute out_q = N * in_u
//   or    out_u = in_q * N
// The u quantity is a vector v_F in the parent frame (e.g. angular 
// velocity of child in parent), while the q quantity is the derivative 
// of a quaternion representing R_FM, i.e. orientation of child in parent.
void multiplyByN(const SBStateDigest& sbs, bool matrixOnRight, 
                 const Real* in, Real* out) const
{
    assert(sbs.getStage() >= Stage::Position);
    assert(in && out);

    const SBModelVars&          mv   = sbs.getModelVars();
    const SBTreePositionCache&  pc   = sbs.getTreePositionCache();
    const Vector&               allQ = sbs.getQ();

    if (getUseEulerAngles(mv)) {
        const Rotation& R_FM = getX_FM(pc).R();
        const Mat33 N_FM = // N here mixes parent and body frames
            Rotation::calcNForBodyXYZInBodyFrame(fromQ(allQ)) ;
        if (matrixOnRight) 
                Row3::updAs(out) = (Row3::getAs(in)*N_FM) * ~R_FM;
        else    Vec3::updAs(out) = N_FM * (~R_FM*Vec3::getAs(in));
    } else {
        // Quaternion
        const Mat43 N_FF = // This method returns N in the parent frame
            Rotation::calcUnnormalizedNForQuaternion(fromQuat(allQ));
        if (matrixOnRight) 
                Row3::updAs(out) = Row4::getAs(in) * N_FF;
        else    Vec4::updAs(out) = N_FF * Vec3::getAs(in);
    }
}

// Compute out_u = inv(N) * in_q
//   or    out_q = in_u * inv(N)
void multiplyByNInv(const SBStateDigest& sbs, bool matrixOnRight,
                    const Real* in, Real* out) const
{
    assert(sbs.getStage() >= Stage::Position);
    assert(in && out);

    const SBModelVars&          mv   = sbs.getModelVars();
    const SBTreePositionCache&  pc   = sbs.getTreePositionCache();
    const Vector&               allQ = sbs.getQ();

    if (getUseEulerAngles(mv)) {
        const Rotation& R_FM = getX_FM(pc).R();
        const Mat33 NInv_MF = // NInv here mixes parent and body frames
            Rotation::calcNInvForBodyXYZInBodyFrame(fromQ(allQ));
        if (matrixOnRight) 
                Row3::updAs(out) = (Row3::getAs(in)*R_FM) * NInv_MF;
        else    Vec3::updAs(out) = R_FM * (NInv_MF*Vec3::getAs(in));
    } else {
        // Quaternion
        const Mat34 NInv_FF = // This method returns NInv in  parent frame
            Rotation::calcUnnormalizedNInvForQuaternion(fromQuat(allQ));
        if (matrixOnRight) 
                Row4::updAs(out) = Row3::getAs(in) * NInv_FF;
        else    Vec3::updAs(out) = NInv_FF * Vec4::getAs(in);
    }
}

// Compute out_q = NDot * in_u
//   or    out_u = in_q * NDot
void multiplyByNDot(const SBStateDigest& sbs, bool matrixOnRight,
                    const Real* in, Real* out) const
{
    assert(sbs.getStage() >= Stage::Velocity);
    assert(in && out);

    const SBModelVars&          mv      = sbs.getModelVars();
    const SBTreePositionCache&  pc      = sbs.getTreePositionCache();
    const Vector&               allQ    = sbs.getQ();
    const Vector&               allQDot = sbs.getQDot();

    if (getUseEulerAngles(mv)) {
        const Rotation& R_FM = getX_FM(pc).R();
        const Mat33 NDot_FM = // NDot here mixes parent and body frames
            Rotation::calcNDotForBodyXYZInBodyFrame(fromQ(allQ), 
                                                    fromQ(allQDot)) ;
        if (matrixOnRight) 
                Row3::updAs(out) = (Row3::getAs(in)*NDot_FM) * ~R_FM;
        else    Vec3::updAs(out) = NDot_FM * (~R_FM*Vec3::getAs(in));
    } else {
        // Quaternion
        const Mat43 NDot_FF = // This method returns NDot in the parent frame
            Rotation::calcUnnormalizedNDotForQuaternion(fromQuat(allQDot));
        if (matrixOnRight) 
                Row3::updAs(out) = Row4::getAs(in) * NDot_FF;
        else    Vec4::updAs(out) = NDot_FF * Vec3::getAs(in);
    }
}

void calcQDotDot(
    const SBStateDigest&   sbs,
    const Vector&          udot, 
    Vector&                qdotdot) const 
{
    const SBModelVars&          mv = sbs.getModelVars();
    const SBTreePositionCache&  pc = sbs.getTreePositionCache();
    // angular velocity of M in F, expr in F
    const Vec3& w_FM     = fromU(sbs.getU()); 
    const Vec3& w_FM_dot = fromU(udot);

    if (getUseEulerAngles(mv)) {
        toQuat(qdotdot) = Vec4(0); // TODO: kludge, clear unused element
        const Rotation& R_FM = getX_FM(pc).R();
        toQ(qdotdot)    = Rotation::convertAngVelDotInBodyFrameToBodyXYZDotDot
                              (fromQ(sbs.getQ()), ~R_FM*w_FM, ~R_FM*w_FM_dot);
    } else
        toQuat(qdotdot) = Rotation::convertAngVelDotToQuaternionDotDot
                              (fromQuat(sbs.getQ()),w_FM,w_FM_dot);
}

void calcLocalQDotDotFromLocalUDot(const SBStateDigest& sbs, const Real* udot, Real* qdotdot) const {
    assert(sbs.getStage() >= Stage::Velocity);
    assert(udot && qdotdot);

    const SBModelVars&          mv   = sbs.getModelVars();
    const SBTreePositionCache&  pc   = sbs.getTreePositionCache();
    const Vector&               allQ = sbs.getQ();
    const Vector&               allU = sbs.getU();

    const Vec3& w_FM     = fromU(allU);
    const Vec3& w_FM_dot = Vec3::getAs(udot);

    if (getUseEulerAngles(mv)) {
        Vec4::updAs(qdotdot) = 0; // TODO: kludge, clear unused element
        const Rotation& R_FM = getX_FM(pc).R();
        Vec3::updAs(qdotdot) = Rotation::convertAngVelDotInBodyFrameToBodyXYZDotDot
                                    (fromQ(allQ), ~R_FM*w_FM, ~R_FM*w_FM_dot);
    } else
        Vec4::updAs(qdotdot) = Rotation::convertAngVelDotToQuaternionDotDot
                                    (fromQuat(allQ),w_FM,w_FM_dot);
}

void copyQ(
    const SBModelVars& mv, 
    const Vector&      qIn, 
    Vector&            q) const 
{
    if (getUseEulerAngles(mv))
        toQ(q) = fromQ(qIn);
    else
        toQuat(q) = fromQuat(qIn);
}

int getMaxNQ()              const {return 4;}
int getNQInUse(const SBModelVars& mv) const {
    return getUseEulerAngles(mv) ? 3 : 4;
} 
bool isUsingQuaternion(const SBStateDigest& sbs, MobilizerQIndex& startOfQuaternion) const {
    if (getUseEulerAngles(sbs.getModelVars())) {startOfQuaternion.invalidate(); return false;}
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

