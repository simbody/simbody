#ifndef SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_LINEORIENTATION_H_
#define SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_LINEORIENTATION_H_

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
 * Contributors: Paul Mitiguy                                                 *
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
 * Define the RigidBodyNode that implements a LineOrientation mobilizer.
 */

#include "SimbodyMatterSubsystemRep.h"
#include "RigidBodyNode.h"
#include "RigidBodyNodeSpec.h"


    // LINE ORIENTATION //

// LineOrientation joint. Like a Ball joint, this provides full rotational
// freedom, but for a degenerate body which is thin (inertialess) along its
// own z axis. These arise in molecular modeling for linear molecules formed
// by pairs of atoms, or by multiple atoms in a linear arrangement like
// carbon dioxide (CO2) whose structure is O=C=O in a straight line. We are
// assuming that there is no meaning to a rotation about the linear axis,
// so free orientation requires just *two* degrees of freedom, not *three*
// as is required for general rigid bodies. And in fact we can get away with
// just two generalized speeds. But so far, no one has been able to come up
// with a way to manage with only two generalized coordinates, so this joint
// has the same q's as a regular Orientation (Ball) joint: either a quaternion
// for unconditional stability, or a three-angle (body fixed 1-2-3)
// Euler sequence which will be dynamically singular when the middle (y) axis
// is 90 degrees. Use the Euler sequence only for small motions or for 
// kinematics problems (and note that only the first two are meaningful).
//
// To summarize, the generalized coordinates are:
//   * 4 quaternions or 3 1-2-3 body fixed Euler angles (that is, fixed in M)
// and generalized speeds are:
//   * the x,y components of the angular velocity w_FM_M, that is, the angular
//     velocity of M in F expressed in M (where we want wz=0).
//     NOTE: THAT IS A DIFFERENT FRAME THAN IS USED FOR BALL AND GIMBAL
// Thus the qdots have to be derived from the generalized speeds to
// be turned into either 4 quaternion derivatives or 3 Euler angle derivatives.
class RBNodeLineOrientation : public RigidBodyNodeSpec<2> {
public:

virtual const char* type() { return "lineOrientation"; }

RBNodeLineOrientation(const MassProperties& mProps_B,
                      const Transform&      X_PF,
                      const Transform&      X_BM,
                      bool                  isReversed,
                      UIndex&               nextUSlot,
                      USquaredIndex&        nextUSqSlot,
                      QIndex&               nextQSlot)
  : RigidBodyNodeSpec<2>(mProps_B,X_PF,X_BM,nextUSlot,nextUSqSlot,nextQSlot,
                         QDotMayDifferFromU, QuaternionMayBeUsed, isReversed)
{
    updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
}

void setQToFitRotationImpl(const SBStateDigest& sbs, const Rotation& R_FM,
                          Vector& q) const 
{
    if (getUseEulerAngles(sbs.getModelVars()))
        toQVec3(q,0)    = R_FM.convertRotationToBodyFixedXYZ();
    else
        toQuat(q) = R_FM.convertRotationToQuaternion().asVec4();
}

void setQToFitTranslationImpl(const SBStateDigest& sbs, const Vec3& p_FM, 
                              Vector& q) const {
    // M and F frame origins are always coincident for this mobilizer so there is no
    // way to create a translation by rotating. So the only translation we can represent is 0.
}

void setUToFitAngularVelocityImpl(const SBStateDigest& sbs, 
                                  const Vector& q, const Vec3& w_FM,
                                  Vector& u) const
{
    Rotation R_FM;
    if (getUseEulerAngles(sbs.getModelVars()))
        R_FM.setRotationToBodyFixedXYZ( fromQVec3(q,0) );
    else {
        // TODO: should use qnorm pool
        R_FM.setRotationFromQuaternion( Quaternion(fromQuat(q)) ); // normalize
    }
    const Vec3 w_FM_M = ~R_FM*w_FM;
    toU(u) = Vec2(w_FM_M[0], w_FM_M[1]); // (x,y) of relative angular velocity always used as generalized speeds
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
        assert(q && nq==3 && qCache && nQCache==AnglePoolSize && nQErr==0);
        const Real cy = std::cos(q[1]);
        Vec3::updAs(&qCache[AngleCosQ]) =
            Vec3(std::cos(q[0]), cy, std::cos(q[2]));
        Vec3::updAs(&qCache[AngleSinQ]) =
            Vec3(std::sin(q[0]), std::sin(q[1]), std::sin(q[2]));
        qCache[AngleOOCosQy] = 1/cy; // trouble at 90 degrees
    } else {
        assert(q && nq==4 && qCache && nQCache==QuatPoolSize && 
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
    X_F0M0.updP() = 0.; // This joint can't translate.
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
}

// The generalized speeds for this 2-dof rotational joint are the x and y
// components of the angular velocity of M in the F frame, expressed in the *M*
// frame.
void calcAcrossJointVelocityJacobian(
    const SBStateDigest& sbs,
    HType&               H_FM) const
{
    const SBTreePositionCache& pc = sbs.updTreePositionCache(); // "upd" because we're realizing positions now
    const Transform X_F0M0 = findX_F0M0(pc);

    // Dropping the 0's here.
    const Rotation& R_FM = X_F0M0.R();
    const Vec3&     Mx_F = R_FM.x(); // M's x axis, expressed in F
    const Vec3&     My_F = R_FM.y(); // M's y axis, expressed in F

    H_FM(0) = SpatialVec( Mx_F, Vec3(0) );
    H_FM(1) = SpatialVec( My_F, Vec3(0) );
}

// Since the Jacobian above is not constant in F,
// its time derivative is non zero. Here we use the fact that for
// a vector r_B_A fixed in a moving frame B but expressed in another frame A,
// its time derivative in A is the angular velocity of B in A crossed with
// the vector, i.e., d_A/dt r_B_A = w_AB % r_B_A.
void calcAcrossJointVelocityJacobianDot(
    const SBStateDigest& sbs,
    HType&               HDot_FM) const
{
    const SBTreePositionCache& pc = sbs.getTreePositionCache();
    const SBTreeVelocityCache& vc = sbs.updTreeVelocityCache(); // "upd" because we're realizing velocities now
    const Transform  X_F0M0 = findX_F0M0(pc);

    // Dropping the 0's here.
    const Rotation& R_FM = X_F0M0.R();
    const Vec3&     Mx_F = R_FM.x(); // M's x axis, expressed in F
    const Vec3&     My_F = R_FM.y(); // M's y axis, expressed in F

    const Vec3      w_FM = find_w_F0M0(pc,vc); // angular velocity of M in F

    HDot_FM(0) = SpatialVec( w_FM % Mx_F, Vec3(0) );
    HDot_FM(1) = SpatialVec( w_FM % My_F, Vec3(0) );
}

// CAUTION: we do not zero the unused 4th element of q for Euler angles; it
// is up to the caller to do that if it is necessary.
// Compute out_q = N * in_u
//   or    out_u = in_q * N
// The u quantity is a vector v_M in the child body frame (e.g. angular 
// velocity of child in parent, but expressed in child), while the q 
// quantity is the derivative of a quaternion or Euler sequence representing 
// R_FM, i.e. orientation of child in parent.
void multiplyByN(const SBStateDigest& sbs, bool matrixOnRight,
                 const Real* in, Real* out) const
{
    assert(sbs.getStage() >= Stage::Position);
    assert(in && out);

    const SBModelVars&          mv   = sbs.getModelVars();
    const SBTreePositionCache&  pc   = sbs.getTreePositionCache();
    const Vector&               allQ = sbs.getQ();

    if (getUseEulerAngles(mv)) {
        const Mat32 N_FM = // N here mixes parent and body frames (convenient)
            Rotation::calcNForBodyXYZInBodyFrame(fromQVec3(allQ,0))
                        .getSubMat<3,2>(0,0); // drop 3rd column
        if (matrixOnRight) Row2::updAs(out) = Row3::getAs(in) * N_FM;
        else               Vec3::updAs(out) = N_FM * Vec2::getAs(in);
    } else {
        // Quaternion: N block is only available expecting angular velocity 
        // in the parent frame F, but we have it in M for this joint so we
        // have to calculate N_FM = N_FF*R_FM.
        const Rotation& R_FM = getX_FM(pc).R();
        const Mat42 N_FM = (Rotation::calcUnnormalizedNForQuaternion(fromQuat(allQ))*R_FM)
                            .getSubMat<4,2>(0,0); // drop 3rd column
        if (matrixOnRight) Row2::updAs(out) = Row4::getAs(in) * N_FM;
        else               Vec4::updAs(out) = N_FM * Vec2::getAs(in);
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
        const Mat23 NInv_MF = Rotation::calcNInvForBodyXYZInBodyFrame(fromQVec3(allQ,0))
                                .getSubMat<2,3>(0,0);   // drop 3rd row
        if (matrixOnRight) Row3::updAs(out) = Row2::getAs(in) * NInv_MF;
        else               Vec2::updAs(out) = NInv_MF * Vec3::getAs(in);
    } else {
        // Quaternion: NInv block is only available expecting angular 
        // velocity in the parent frame F, but we have it in M for 
        // this joint so we have to calculate NInv_MF = R_MF*NInv_FF.
        const Rotation& R_FM = getX_FM(pc).R();
        const Mat24 NInv_MF = (~R_FM*Rotation::calcUnnormalizedNInvForQuaternion(fromQuat(allQ)))
                                .getSubMat<2,4>(0,0);   // drop 3rd row
        if (matrixOnRight) Row4::updAs(out) = Row2::getAs(in) * NInv_MF;
        else               Vec2::updAs(out) = NInv_MF * Vec4::getAs(in);
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
        const Mat32 NDot_FM = // NDot here mixes parent and body frames (convenient)
            Rotation::calcNDotForBodyXYZInBodyFrame(fromQVec3(allQ,0), 
                                                    fromQVec3(allQDot,0))
                        .getSubMat<3,2>(0,0); // drop 3rd column
        if (matrixOnRight) Row2::updAs(out) = Row3::getAs(in) * NDot_FM;
        else               Vec3::updAs(out) = NDot_FM * Vec2::getAs(in);
    } else {
        // Quaternion: NDot block is only available expecting angular velocity 
        // in the parent frame F, but we have it in M for this joint so we
        // have to calculate NDot_FM = NDot_FF*R_FM.
        const Rotation& R_FM = getX_FM(pc).R();
        const Mat42 NDot_FM = 
           (Rotation::calcUnnormalizedNDotForQuaternion(fromQuat(allQDot))*R_FM)
                        .getSubMat<4,2>(0,0); // drop 3rd column
        if (matrixOnRight) Row2::updAs(out) = Row4::getAs(in) * NDot_FM;
        else               Vec4::updAs(out) = NDot_FM * Vec2::getAs(in);
    }
}

void calcQDot(
    const SBStateDigest&   sbs,
    const Real*            u, 
    Real*                  qdot) const 
{
    const SBModelVars& mv = sbs.getModelVars();
    const SBTreePositionCache& pc = sbs.getTreePositionCache();
    const Vec3 w_FM_M = Vec2::getAs(u).append1(0); // angular velocity of M in F, exp in M (with wz=0) 
    if (getUseEulerAngles(mv)) {
        Vec3::updAs(qdot) = Rotation::convertAngVelInBodyFrameToBodyXYZDot(fromQVec3(sbs.getQ(),0),
                                    w_FM_M); // need w in *body*, not parent
        qdot[3] = 0;
    } else {
        const Rotation& R_FM = getX_FM(pc).R();
        Vec4::updAs(qdot) = Rotation::convertAngVelToQuaternionDot(fromQuat(sbs.getQ()),
                                    R_FM*w_FM_M); // need w in *parent* frame
    }
}

void calcQDotDot(
    const SBStateDigest&   sbs,
    const Real*            udot, 
    Real*                  qdotdot) const 
{
    const SBModelVars&          mv = sbs.getModelVars();
    const SBTreePositionCache&  pc = sbs.getTreePositionCache();
    const Vec3 w_FM_M     = fromU(sbs.getU()).append1(0); // angular velocity of M in F, exp in M (with wz=0)
    const Vec3 w_FM_M_dot = Vec2::getAs(udot).append1(0);

    if (getUseEulerAngles(mv)) {
        Vec3::updAs(qdotdot) = Rotation::convertAngVelDotInBodyFrameToBodyXYZDotDot
                                   (fromQVec3(sbs.getQ(),0), w_FM_M, w_FM_M_dot); // body frame
        qdotdot[3] = 0;
    } else {
        const Rotation& R_FM = getX_FM(pc).R();
        Vec4::updAs(qdotdot) = Rotation::convertAngVelDotToQuaternionDotDot
                              (fromQuat(sbs.getQ()),R_FM*w_FM_M,R_FM*w_FM_M_dot); // parent frame
    }
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
    toQVec3(outputQ, 4) = Vec3(0); // clear unused element
    toQVec3(outputQ, 2) = fromQVec3(inputQ, 3);
    toQVec3(outputQ, 0) = Rotation(Quaternion(fromQuat(inputQ))).convertRotationToBodyFixedXYZ();
}

void convertToQuaternions(const Vector& inputQ, Vector& outputQ) const {
    toQVec3(outputQ, 4) = fromQVec3(inputQ, 3);
    Rotation rot;
    rot.setRotationToBodyFixedXYZ(fromQVec3(inputQ, 0));
    toQuat(outputQ) = rot.convertRotationToQuaternion().asVec4();
}


};



#endif // SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_LINEORIENTATION_H_

