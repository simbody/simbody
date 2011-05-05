#ifndef SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_FREE_H_
#define SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_FREE_H_

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
 * Define the RigidBodyNode that implements a Free mobilizer.
 */

#include "SimbodyMatterSubsystemRep.h"
#include "RigidBodyNode.h"
#include "RigidBodyNodeSpec.h"


    // FREE //

// Free joint. This provides six degrees of freedom, three rotational and
// three translational. The rotation is like the ball joint above; the
// translation is like the Cartesian joint above.
// TODO: to get this to work I had to make the translations be in the outboard
// frame (M, not F). So currently the generalized coordinates are:
//   * 4 quaternions or 3 1-2-3 body fixed Euler angles (that is, fixed in M)
//   * translation from OF to OM as a 3-vector in the outboard body 
//        mobilizer (M) frame
// and generalized speeds are:
//   * angular velocity w_FM as a vector expressed in the F frame
//   * linear velocity of the M origin in F (v_FM), expressed in M
// Thus translational qdots are just generalized speeds, but rotational
// qdots have to be derived from the generalized speeds to be turned into
// either 4 quaternion derivatives or 3 Euler angle derivatives.
//   
class RBNodeFree : public RigidBodyNodeSpec<6> {
public:

virtual const char* type() { return "free"; }

RBNodeFree(const MassProperties& mProps_B,
           const Transform&      X_PF,
           const Transform&      X_BM,
           bool                  isReversed,
           UIndex&               nextUSlot,
           USquaredIndex&        nextUSqSlot,
           QIndex&               nextQSlot)
  : RigidBodyNodeSpec<6>(mProps_B,X_PF,X_BM,nextUSlot,nextUSqSlot,nextQSlot,
                         QDotMayDifferFromU, QuaternionMayBeUsed, isReversed)
{
    updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
}

void setQToFitRotationImpl(const SBStateDigest& sbs, const Rotation& R_FM,
                          Vector& q) const 
{
    if (getUseEulerAngles(sbs.getModelVars()))
        toQVec3(q,0) = R_FM.convertRotationToBodyFixedXYZ();
    else
        toQuat(q) = R_FM.convertRotationToQuaternion().asVec4();
}

// The user gives us the translation vector from OF to OM as a vector expressed
// in F, which is what we use as translational generalized coordinates. Also, 
// with a free joint we never have to change orientation coordinates in order 
// to achieve a translation.
void setQToFitTranslationImpl(const SBStateDigest& sbs, const Vec3& p_FM, Vector& q) const {
    if (getUseEulerAngles(sbs.getModelVars()))
        toQVec3(q,3) = p_FM; // skip the 3 Euler angles
    else
        toQVec3(q,4) = p_FM; // skip the 4 quaternions
}

// Our 3 rotational generalized speeds are just the angular velocity vector of 
// M in F, expressed in F, which is exactly what the user provides here.
void setUToFitAngularVelocityImpl(const SBStateDigest& sbs, 
                                  const Vector& q, const Vec3& w_FM,
                                  Vector& u) const
{
    toUVec3(u,0) = w_FM; // relative ang. vel. always used as generalized speeds
}

// Our 3 translational generalized speeds are the linear velocity of M's origin 
// in F, expressed in F, which is just what the user gives us.
void setUToFitLinearVelocityImpl
   (const SBStateDigest& sbs, const Vector& q, const Vec3& v_FM, Vector& u) const
{
    toUVec3(u,3) = v_FM;
}

// When we're using Euler angles, we're going to want to cache cos and sin for
// each angle, and also 1/cos of the middle angle will be handy to have around.
// When we're using quaternions, we'll cache 1/|q| for convenient normalizing
// of quaternions.
enum {AnglePoolSize=7, QuatPoolSize=1};
// cos x,y,z sin x,y,z 1/cos(y)
enum {AngleCosQ=0, AngleSinQ=3, AngleOOCosQy=6};
enum {QuatOONorm=0};
int calcQPoolSize(const SBModelVars& mv) const
{   return getUseEulerAngles(mv) ? AnglePoolSize : QuatPoolSize; }

void performQPrecalculations(const SBStateDigest& sbs,
                             const Real* q,      int nq,
                             Real*       qCache, int nQCache,
                             Real*       qErr,   int nQErr) const
{
    if (getUseEulerAngles(sbs.getModelVars())) {
        assert(q && nq==6 && qCache && nQCache==AnglePoolSize && nQErr==0);
        const Real cy = std::cos(q[1]);
        Vec3::updAs(&qCache[AngleCosQ]) =
            Vec3(std::cos(q[0]), cy, std::cos(q[2]));
        Vec3::updAs(&qCache[AngleSinQ]) =
            Vec3(std::sin(q[0]), std::sin(q[1]), std::sin(q[2]));
        qCache[AngleOOCosQy] = 1/cy; // trouble at 90 degrees
    } else {
        assert(q && nq==7 && qCache && nQCache==QuatPoolSize && 
               qErr && nQErr==1);
        const Real quatLen = Vec4::getAs(q).norm();
        qErr[0] = quatLen - Real(1);    // normalization error
        qCache[QuatOONorm] = 1/quatLen; // save for later
    }
}

void calcX_FM(const SBStateDigest& sbs,
              const Real* q,      int nq,
              const Real* qCache, int nQCache,
              Transform&  X_F0M0) const
{
    if (getUseEulerAngles(sbs.getModelVars())) {
        assert(q && nq==6 && qCache && nQCache==AnglePoolSize);
        X_F0M0.updR().setRotationToBodyFixedXYZ // 18 flops
           (Vec3::getAs(&qCache[AngleCosQ]), Vec3::getAs(&qCache[AngleSinQ]));
        X_F0M0.updP() = Vec3::getAs(&q[3]); // a012 x y z
    } else {
        assert(q && nq==7 && qCache && nQCache==QuatPoolSize);
        // Must use a normalized quaternion to generate the rotation matrix.
        // Here we normalize with just 4 flops using precalculated 1/norm(q).
        const Quaternion quat(Vec4::getAs(q)*qCache[QuatOONorm], true); 
        X_F0M0.updR().setRotationFromQuaternion(quat); // 29 flops
        X_F0M0.updP() = Vec3::getAs(&q[4]); // q0123 x y z
    }
}


// The generalized speeds for this 6-dof ("free") joint are 
//   (1) the angular velocity of M in the F frame, expressed in F, and
//   (2) the (linear) velocity of M's origin in F, expressed in F.
void calcAcrossJointVelocityJacobian(
    const SBStateDigest& sbs,
    HType&               H_FM) const
{
    H_FM(0) = SpatialVec( Vec3(1,0,0),   Vec3(0)   );  // rotations
    H_FM(1) = SpatialVec( Vec3(0,1,0),   Vec3(0)   );
    H_FM(2) = SpatialVec( Vec3(0,0,1),   Vec3(0)   );

    H_FM(3) = SpatialVec(   Vec3(0),   Vec3(1,0,0) );  // translations
    H_FM(4) = SpatialVec(   Vec3(0),   Vec3(0,1,0) );
    H_FM(5) = SpatialVec(   Vec3(0),   Vec3(0,0,1) );
}

// Since the Jacobian above is constant in F, its derivative in F is 0.
void calcAcrossJointVelocityJacobianDot(
    const SBStateDigest& sbs,
    HType&               HDot_FM) const
{
    HDot_FM(0) = SpatialVec( Vec3(0), Vec3(0) );
    HDot_FM(1) = SpatialVec( Vec3(0), Vec3(0) );
    HDot_FM(2) = SpatialVec( Vec3(0), Vec3(0) );

    HDot_FM(3) = SpatialVec( Vec3(0), Vec3(0) );
    HDot_FM(4) = SpatialVec( Vec3(0), Vec3(0) );
    HDot_FM(5) = SpatialVec( Vec3(0), Vec3(0) );
}

// CAUTION: we do not zero the unused 4th element of q for Euler angles; it
// is up to the caller to do that if it is necessary.
// Compute out_q = N * in_u
//   or    out_u = in_q * N
// The u quantity is a vector v_F in the parent frame (e.g. angular 
// velocity of child in parent), while the q quantity is the derivative 
// of a quaternion or Euler sequence representing R_FM, i.e. orientation 
// of child in parent.
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
            Rotation::calcNForBodyXYZInBodyFrame(fromQVec3(allQ,0)) ;
        if (matrixOnRight) 
                Row3::updAs(out) = (Row3::getAs(in)*N_FM) * ~R_FM;
        else    Vec3::updAs(out) = N_FM * (~R_FM*Vec3::getAs(in));
        // Translational part of N block is identity.
        Vec3::updAs(out+3) = Vec3::getAs(in+3);
    } else {
        // Quaternion
        const Mat43 N_FF = // This method returns N in the parent frame
            Rotation::calcUnnormalizedNForQuaternion(fromQuat(allQ));
        // Translational part of N is identity so just copy in to out.
        if (matrixOnRight) {
            Row3::updAs(out) = Row4::getAs(in) * N_FF;
            Row3::updAs(out+3) = Row3::getAs(in+4);     // translation
        } else {
            Vec4::updAs(out) = N_FF * Vec3::getAs(in);
            Vec3::updAs(out+4) = Vec3::getAs(in+3);     // translation
        }
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
            Rotation::calcNInvForBodyXYZInBodyFrame(fromQVec3(allQ,0));
        if (matrixOnRight) 
                Row3::updAs(out) = (Row3::getAs(in)*R_FM) * NInv_MF;
        else    Vec3::updAs(out) = R_FM * (NInv_MF*Vec3::getAs(in));
        // Translational part of NInv block is identity.
        Vec3::updAs(out+3) = Vec3::getAs(in+3);
    } else {
        // Quaternion
        const Mat34 NInv_FF = // This method returns NInv in  parent frame
            Rotation::calcUnnormalizedNInvForQuaternion(fromQuat(allQ));
        // Translational part of NInv is identity so just copy in to out.
        if (matrixOnRight) {
            Row4::updAs(out) = Row3::getAs(in) * NInv_FF;
            Row3::updAs(out+4) = Row3::getAs(in+3); // translation
        } else { 
            Vec3::updAs(out) = NInv_FF * Vec4::getAs(in);
            Vec3::updAs(out+3) = Vec3::getAs(in+4); // translation
        }
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
            Rotation::calcNDotForBodyXYZInBodyFrame(fromQVec3(allQ,0), 
                                                    fromQVec3(allQDot,0));
        if (matrixOnRight) 
                Row3::updAs(out) = (Row3::getAs(in)*NDot_FM) * ~R_FM;
        else    Vec3::updAs(out) = NDot_FM * (~R_FM*Vec3::getAs(in));
        // Translational part of NDot block is zero.
        Vec3::updAs(out+3) = 0;
    } else {
        // Quaternion
        const Mat43 NDot_FF = // This method returns NDot in the parent frame
            Rotation::calcUnnormalizedNDotForQuaternion(fromQuat(allQDot));
        // Translational part of NDot is zero so set out to zero.
        if (matrixOnRight) {
            Row3::updAs(out) = Row4::getAs(in) * NDot_FF;
            Row3::updAs(out+3) = 0;     // translation
        } else {
            Vec4::updAs(out) = NDot_FF * Vec3::getAs(in);
            Vec3::updAs(out+4) = 0;     // translation
        }
    }
}

void calcQDot(
    const SBStateDigest&   sbs,
    const Vector&          u,
    Vector&                qdot) const
{
    const SBModelVars&          mv = sbs.getModelVars();
    const SBTreePositionCache&  pc = sbs.getTreePositionCache();
    const Vec3& w_FM = fromUVec3(u,0); // Angular velocity in F
    const Vec3& v_FM = fromUVec3(u,3); // Linear velocity in F
    if (getUseEulerAngles(mv)) {
        const Rotation& R_FM = getX_FM(pc).R();
        const Vec3& theta = fromQVec3(sbs.getQ(),0); // Euler angles
        toQVec3(qdot,0) = Rotation::convertAngVelInBodyFrameToBodyXYZDot(theta,
                                        ~R_FM*w_FM); // need w in *body*, not parent
        toQVec3(qdot,4) = Vec3(0); // TODO: kludge, clear unused element
        toQVec3(qdot,3) = v_FM;
    } else {
        const Vec4& quat = fromQuat(sbs.getQ());
        toQuat (qdot)   = Rotation::convertAngVelToQuaternionDot(quat,w_FM);
        toQVec3(qdot,4) = v_FM;
    }
}

void calcQDotDot(
    const SBStateDigest&   sbs,
    const Vector&          udot, 
    Vector&                qdotdot) const 
{
    const SBModelVars&          mv = sbs.getModelVars();
    const SBTreePositionCache&  pc = sbs.getTreePositionCache();
    const Vec3& w_FM     = fromUVec3(sbs.getU(),0); // angular velocity of M in F
    const Vec3& v_FM     = fromUVec3(sbs.getU(),3); // linear velocity of M in F, expressed in F
    const Vec3& w_FM_dot = fromUVec3(udot,0);
    const Vec3& v_FM_dot = fromUVec3(udot,3);
    if (getUseEulerAngles(mv)) {
        const Rotation& R_FM = getX_FM(pc).R();
        const Vec3& theta  = fromQVec3(sbs.getQ(),0); // Euler angles
        toQVec3(qdotdot,0) = Rotation::convertAngVelDotInBodyFrameToBodyXYZDotDot
                                         (theta, ~R_FM*w_FM, ~R_FM*w_FM_dot);
        toQVec3(qdotdot,4) = Vec3(0); // TODO: kludge, clear unused element
        toQVec3(qdotdot,3) = v_FM_dot;
    } else {
        const Vec4& quat  = fromQuat(sbs.getQ());
        toQuat(qdotdot)   = Rotation::convertAngVelDotToQuaternionDotDot
                                         (quat,w_FM,w_FM_dot);
        toQVec3(qdotdot,4) = v_FM_dot;
    }
}

void copyQ(const SBModelVars& mv, const Vector& qIn, Vector& q) const {
    if (getUseEulerAngles(mv))
        toQ(q) = fromQ(qIn);
    else {
        toQuat(q)    = fromQuat(qIn);
        toQVec3(q,4) = fromQVec3(qIn,4);
    }
}

int  getMaxNQ()                   const {return 7;}
int  getNQInUse(const SBModelVars& mv) const {return getUseEulerAngles(mv) ? 6 : 7;} 
bool isUsingQuaternion(const SBStateDigest& sbs, MobilizerQIndex& startOfQuaternion) const {
    if (getUseEulerAngles(sbs.getModelVars())) {startOfQuaternion.invalidate(); return false;}
    startOfQuaternion = MobilizerQIndex(0); // quaternion comes first
    return true;
}

void setMobilizerDefaultPositionValues(const SBModelVars& mv, Vector& q) const 
{
    if (getUseEulerAngles(mv)) {
        toQVec3(q,4) = Vec3(0); // TODO: kludge, clear unused element
        toQ(q) = 0.;
    } else {
        toQuat(q) = Vec4(1.,0.,0.,0.);
        toQVec3(q,4) = 0.;
    }
}

bool enforceQuaternionConstraints(
    const SBStateDigest& sbs, 
    Vector&             q,
    Vector&             qErrest) const
{
    if (getUseEulerAngles(sbs.getModelVars())) 
        return false; // no change

    Vec4& quat = toQuat(q);
    quat = quat / quat.norm();

    if (qErrest.size()) {
        Vec4& qerr = toQuat(qErrest);
        qerr -= dot(qerr,quat) * quat;
    }

    return true;
}

void convertToEulerAngles(const Vector& inputQ, Vector& outputQ) const {
    toQVec3(outputQ, 4) = Vec3(0); // clear unused element
    toQVec3(outputQ, 3) = fromQVec3(inputQ, 4);
    toQVec3(outputQ, 0) = Rotation(Quaternion(fromQuat(inputQ))).convertRotationToBodyFixedXYZ();
}

void convertToQuaternions(const Vector& inputQ, Vector& outputQ) const {
    toQVec3(outputQ, 4) = fromQVec3(inputQ, 3);
    Rotation rot;
    rot.setRotationToBodyFixedXYZ(fromQVec3(inputQ, 0));
    toQuat(outputQ) = rot.convertRotationToQuaternion().asVec4();
}


};


 
#endif // SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_FREE_H_

