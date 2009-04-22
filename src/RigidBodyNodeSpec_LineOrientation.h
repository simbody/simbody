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
// is 90 degrees. Use the Euler sequence only for small motions or for kinematics
// problems (and note that only the first two are meaningful).
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

    void setQToFitTranslationImpl(const SBStateDigest& sbs, const Vec3& p_FM, Vector& q) const {
        // M and F frame origins are always coincident for this mobilizer so there is no
        // way to create a translation by rotating. So the only translation we can represent is 0.
    }

    void setUToFitAngularVelocityImpl(const SBStateDigest& sbs, const Vector& q, const Vec3& w_FM,
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

    // This is required for all mobilizers.
    bool isUsingAngles(const SBStateDigest& sbs, MobilizerQIndex& startOfAngles, int& nAngles) const {
        // LineOrientation joint has three angular coordinates when Euler angles are being used, 
        // none when quaternions are being used.
        if (!getUseEulerAngles(sbs.getModelVars())) {startOfAngles.invalidate(); nAngles=0; return false;} 
        startOfAngles = MobilizerQIndex(0);
        nAngles = 3;
        return true;
    }

    // Precalculate sines and cosines.
    void calcJointSinCosQNorm(
        const SBModelVars&  mv,
        const SBModelCache& mc,
        const SBInstanceCache& ic,
        const Vector&       q, 
        Vector&             sine, 
        Vector&             cosine, 
        Vector&             qErr,
        Vector&             qnorm) const
    {
        const SBModelCache::PerMobilizedBodyModelInfo& bInfo = mc.getMobilizedBodyModelInfo(nodeNum);

        if (getUseEulerAngles(mv)) {
            const Vec3& a = fromQVec3(q,0); // angular coordinates
            toQVec3(sine,0)   = Vec3(std::sin(a[0]), std::sin(a[1]), std::sin(a[2]));
            toQVec3(cosine,0) = Vec3(std::cos(a[0]), std::cos(a[1]), std::cos(a[2]));
            // no quaternions
        } else {
            // no angles
            const Vec4& quat = fromQuat(q); // unnormalized quaternion from state
            const Real  quatLen = quat.norm();
            assert(bInfo.hasQuaternionInUse && bInfo.quaternionPoolIndex.isValid());
            qErr[ic.firstQuaternionQErrSlot+bInfo.quaternionPoolIndex] = quatLen - Real(1);
            toQuat(qnorm) = quat / quatLen;
        }
    }

    // Calculate X_F0M0.
    void calcAcrossJointTransform(
        const SBStateDigest& sbs,
        const Vector&        q,
        Transform&           X_F0M0) const
    {
        const SBModelVars& mv = sbs.getModelVars();
        X_F0M0.updP() = 0.; // This joint can't translate.
        if (getUseEulerAngles(mv))
            X_F0M0.updR().setRotationToBodyFixedXYZ( fromQVec3(q,0) );
        else {
            // TODO: should use qnorm pool
            X_F0M0.updR().setRotationFromQuaternion( Quaternion(fromQuat(q)) ); // normalize
        }
    }

    // The generalized speeds for this 2-dof rotational joint are the x and y
    // components of the angular velocity of M in the F frame, expressed in the *M*
    // frame.
    void calcAcrossJointVelocityJacobian(
        const SBStateDigest& sbs,
        HType&               H_FM) const
    {
        const SBPositionCache& pc = sbs.updPositionCache(); // "upd" because we're realizing positions now
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
        const SBPositionCache& pc = sbs.getPositionCache();
        const SBVelocityCache& vc = sbs.updVelocityCache(); // "upd" because we're realizing velocities now
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
    void multiplyByN(const SBStateDigest& sbs, bool useEulerAnglesIfPossible, const Real* q,
                          bool matrixOnRight, const Real* in, Real* out) const
    {
        assert(sbs.getStage() >= Stage::Model);
        assert(q && in && out);

        if (useEulerAnglesIfPossible) {
            const Mat32    N = Rotation::calcQBlockForBodyXYZInBodyFrame(Vec3::getAs(q))
                                    .getSubMat<3,2>(0,0); // drop 3rd column
            if (matrixOnRight) Row2::updAs(out) = Row3::getAs(in) * N;
            else               Vec3::updAs(out) = N * Vec2::getAs(in);
        } else {
            // Quaternion: N block is only available expecting angular velocity in the
            // parent frame F, but we have it in M for this joint.
            const Rotation R_FM(Quaternion(Vec4::getAs(q)));
            const Mat42 N = (Rotation::calcUnnormalizedQBlockForQuaternion(Vec4::getAs(q))*R_FM)
                                .getSubMat<4,2>(0,0); // drop 3rd column
            if (matrixOnRight) Row2::updAs(out) = Row4::getAs(in) * N;
            else               Vec4::updAs(out) = N * Vec2::getAs(in);
        }
    }

    // Compute out_u = inv(N) * in_q
    //   or    out_q = in_u * inv(N)
    void multiplyByNInv(const SBStateDigest& sbs, bool useEulerAnglesIfPossible, const Real* q,
                             bool matrixOnRight, const Real* in, Real* out) const
    {
        assert(sbs.getStage() >= Stage::Position);
        assert(in && out);

        if (useEulerAnglesIfPossible) {
            const Mat23    NInv = Rotation::calcQInvBlockForBodyXYZInBodyFrame(Vec3::getAs(q))
                                        .getSubMat<2,3>(0,0); // drop 3rd row
            if (matrixOnRight) Row3::updAs(out) = Row2::getAs(in) * NInv;
            else               Vec2::updAs(out) = NInv * Vec3::getAs(in);
        } else {
            // Quaternion: QInv block is only available expecting angular velocity in the
            // parent frame F, but we have it in M for this joint.
            const Rotation R_FM(Quaternion(Vec4::getAs(q)));
            const Mat24 NInv = (~R_FM*Rotation::calcUnnormalizedQInvBlockForQuaternion(Vec4::getAs(q)))
                                    .getSubMat<2,4>(0,0);   // drop 3rd row
            if (matrixOnRight) Row4::updAs(out) = Row2::getAs(in) * NInv;
            else               Vec2::updAs(out) = NInv * Vec4::getAs(in);
        }
    }

    void calcQDot(
        const SBStateDigest&   sbs,
        const Vector&          u, 
        Vector&                qdot) const 
    {
        const SBModelVars& mv = sbs.getModelVars();
        const SBPositionCache& pc = sbs.getPositionCache();
        const Vec3 w_FM_M = fromU(u).append1(0); // angular velocity of M in F, exp in M (with wz=0) 
        if (getUseEulerAngles(mv)) {
            toQuat(qdot)    = Vec4(0); // TODO: kludge, clear unused element
            toQVec3(qdot,0) = Rotation::convertAngVelToBodyFixed123Dot(fromQVec3(sbs.getQ(),0),
                                        w_FM_M); // need w in *body*, not parent
        } else {
            const Rotation& R_FM = getX_FM(pc).R();
            toQuat(qdot) = Rotation::convertAngVelToQuaternionDot(fromQuat(sbs.getQ()),
                                        R_FM*w_FM_M); // need w in *parent* frame
        }
    }
 
    void calcQDotDot(
        const SBStateDigest&   sbs,
        const Vector&          udot, 
        Vector&                qdotdot) const 
    {
        const SBModelVars& mv = sbs.getModelVars();
        const SBPositionCache& pc = sbs.getPositionCache();
        const Vec3 w_FM_M     = fromU(sbs.getU()).append1(0); // angular velocity of M in F, exp in M (with wz=0)
        const Vec3 w_FM_M_dot = fromU(udot).append1(0);

        if (getUseEulerAngles(mv)) {
            toQuat(qdotdot)    = Vec4(0); // TODO: kludge, clear unused element
            toQVec3(qdotdot,0) = Rotation::convertAngVelDotToBodyFixed123DotDot
                                       (fromQVec3(sbs.getQ(),0), w_FM_M, w_FM_M_dot); // body frame
        } else {
            const Rotation& R_FM = getX_FM(pc).R();
            toQuat(qdotdot) = Rotation::convertAngVelDotToQuaternionDotDot
                                  (fromQuat(sbs.getQ()),R_FM*w_FM_M,R_FM*w_FM_M_dot); // parent frame
        }
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

