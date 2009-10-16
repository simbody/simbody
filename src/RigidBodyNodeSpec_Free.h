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
//   * translation from OF to OM as a 3-vector in the outboard body mobilizer (M) frame
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

    // The user gives us the translation vector from OF to OM as a vector expressed in F, which
    // is what we use as translational generalized coordinates. Also, with a free joint 
    // we never have to change orientation coordinates in order to achieve a translation.
    void setQToFitTranslationImpl(const SBStateDigest& sbs, const Vec3& p_FM, Vector& q) const {
        if (getUseEulerAngles(sbs.getModelVars()))
            toQVec3(q,3) = p_FM; // skip the 3 Euler angles
        else
            toQVec3(q,4) = p_FM; // skip the 4 quaternions
    }

    // Our 3 rotational generalized speeds are just the angular velocity vector of M in F,
    // expressed in F, which is exactly what the user provides here.
    void setUToFitAngularVelocityImpl(const SBStateDigest& sbs, const Vector& q, const Vec3& w_FM,
                                     Vector& u) const
    {
        toUVec3(u,0) = w_FM; // relative angular velocity always used as generalized speeds
    }

    // Our 3 translational generalized speeds are the linear velocity of M's origin in F,
    // expressed in F, which is just what the user gives us.
    void setUToFitLinearVelocityImpl
       (const SBStateDigest& sbs, const Vector& q, const Vec3& v_FM, Vector& u) const
    {
        toUVec3(u,3) = v_FM;
    }

    // This is required for all mobilizers.
    bool isUsingAngles(const SBStateDigest& sbs, MobilizerQIndex& startOfAngles, int& nAngles) const {
        // Free joint has three angular coordinates when Euler angles are being used, 
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
            const Vec3& a = fromQ(q).getSubVec<3>(0); // angular coordinates
            toQ(sine).updSubVec<3>(0)   = Vec3(std::sin(a[0]), std::sin(a[1]), std::sin(a[2]));
            toQ(cosine).updSubVec<3>(0) = Vec3(std::cos(a[0]), std::cos(a[1]), std::cos(a[2]));
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

    // Calculate X_FM.
    void calcAcrossJointTransform(
        const SBStateDigest& sbs,
        const Vector&        q,
        Transform&           X_FM) const 
    {
        const SBModelVars& mv = sbs.getModelVars();
        if (getUseEulerAngles(mv)) {
            X_FM.updR().setRotationToBodyFixedXYZ( fromQVec3(q,0) );
            X_FM.updP() = fromQVec3(q,3); // translation is in F already
        } else {
            X_FM.updR().setRotationFromQuaternion( Quaternion(fromQuat(q)) ); // normalize
            X_FM.updP() = fromQVec3(q,4); // translation is in F already
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
    void multiplyByN(const SBStateDigest& sbs, bool useEulerAnglesIfPossible, const Real* q,
                          bool matrixOnRight, const Real* in, Real* out) const
    {
        assert(sbs.getStage() >= Stage::Model);
        assert(q && in && out);

        if (useEulerAnglesIfPossible) {
            // TODO: it's annoying that this Q block is only available in the Body (M) frame,
            // because this mobilizer uses angular velocity in the Parent (F) frame
            // as generalized speeds. So we have to do an expensive conversion here. It
            // would be just as easy to compute this matrix in the Parent frame in 
            // the first place.
            const Rotation R_FM(BodyRotationSequence, q[0], XAxis, q[1], YAxis, q[2], ZAxis);
            const Mat33    N = Rotation::calcQBlockForBodyXYZInBodyFrame(Vec3::getAs(q)) * ~R_FM;
            if (matrixOnRight) Row3::updAs(out) = Row3::getAs(in) * N;
            else               Vec3::updAs(out) = N * Vec3::getAs(in);
            // translational part of Q block is identity
            Vec3::updAs(out+3) = Vec3::getAs(in+3);
        } else {
            // Quaternion
            const Mat43 N = Rotation::calcUnnormalizedQBlockForQuaternion(Vec4::getAs(q));
            if (matrixOnRight) {
                Row3::updAs(out)   = Row4::getAs(in) * N;
                Row3::updAs(out+3) = Row3::getAs(in+4); // translational part of N block is identity
            } else { // matrix on left
                Vec4::updAs(out) = N * Vec3::getAs(in);
                Vec3::updAs(out+4) = Vec3::getAs(in+3); // translational part of N block is identity
            }
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
            // TODO: see above regarding the need for this R_FM kludge
            const Rotation R_FM(BodyRotationSequence, q[0], XAxis, q[1], YAxis, q[2], ZAxis);
            const Mat33    NInv(R_FM*Rotation::calcQInvBlockForBodyXYZInBodyFrame(Vec3::getAs(q)));
            if (matrixOnRight) Row3::updAs(out) = Row3::getAs(in) * NInv;
            else               Vec3::updAs(out) = NInv * Vec3::getAs(in);
            // translational part of NInv block is identity
            Vec3::updAs(out+3) = Vec3::getAs(in+3);
        } else {           
            // Quaternion
            const Mat34 NInv = Rotation::calcUnnormalizedQInvBlockForQuaternion(Vec4::getAs(q));
            if (matrixOnRight) {
                Row4::updAs(out) = Row3::getAs(in) * NInv;
                Row3::updAs(out+4) = Row3::getAs(in+3); // translational part of NInv block is identity
            } else { // matrix on left
                Vec3::updAs(out) = NInv * Vec4::getAs(in);
                Vec3::updAs(out+3) = Vec3::getAs(in+4); // translational part of NInv block is identity
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
            toQVec3(qdot,0) = Rotation::convertAngVelToBodyFixed123Dot(theta,
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
            toQVec3(qdotdot,0) = Rotation::convertAngVelDotToBodyFixed123DotDot
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

