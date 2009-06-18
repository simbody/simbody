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
 * Define the RigidBodyNode that implements a Gimbal mobilizer.
 */

#include "SimbodyMatterSubsystemRep.h"
#include "RigidBodyNode.h"
#include "RigidBodyNodeSpec.h"


    // GIMBAL //

// Gimbal joint. This provides three degrees of rotational freedom,  i.e.,
// unrestricted orientation of the body's M frame in the parent's F frame.
// The generalized coordinates are:
//   * 3 1-2-3 body fixed Euler angles (that is, fixed in M)
// and generalized speeds are:
//   * angular velocity w_FM as a vector expressed in the F frame.
// Thus rotational qdots have to be derived from the generalized speeds to
// be turned into 3 Euler angle derivatives.
//
// NOTE: This joint has a singularity when the middle angle is near 90 degrees.
// In most cases you should use a Ball joint instead, which by default uses
// a quaternion as its generalized coordinates to avoid this singularity.

class RBNodeGimbal : public RigidBodyNodeSpec<3> {
public:
    virtual const char* type() { return "gimbal"; }

    RBNodeGimbal( const MassProperties& mProps_B,
                  const Transform&      X_PF,
                  const Transform&      X_BM,
                  bool                  isReversed,
                  UIndex&               nextUSlot,
                  USquaredIndex&        nextUSqSlot,
                  QIndex&               nextQSlot)
      : RigidBodyNodeSpec<3>(mProps_B,X_PF,X_BM,nextUSlot,nextUSqSlot,nextQSlot,
                             QDotMayDifferFromU, QuaternionIsNeverUsed, isReversed)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }

    void setQToFitRotationImpl(const SBStateDigest& sbs, const Rotation& R_FM,
                           Vector& q) const 
    {
        toQ(q) = R_FM.convertRotationToBodyFixedXYZ();
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
    bool isUsingAngles(const SBStateDigest& sbs, MobilizerQIndex& startOfAngles, int& nAngles) const {
        // Gimbal joint has three angular coordinates.
        startOfAngles = MobilizerQIndex(0); nAngles=3; 
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
        const Vec3& a = fromQ(q); // angular coordinates
        toQ(sine)   = Vec3(std::sin(a[0]), std::sin(a[1]), std::sin(a[2]));
        toQ(cosine) = Vec3(std::cos(a[0]), std::cos(a[1]), std::cos(a[2]));
        // no quaternions
    }

    // Calculate X_FM.
    void calcAcrossJointTransform(
        const SBStateDigest& sbs,
        const Vector&        q,
        Transform&           X_FM) const
    {
        X_FM.updP() = 0.; // This joint can't translate.
        X_FM.updR().setRotationToBodyFixedXYZ( fromQ(q) );
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
        const SBPositionCache& pc = sbs.getPositionCache();
        const Vec3& w_FM = fromU(u); // angular velocity of M in F 
        const Rotation& R_FM = getX_FM(pc).R();
        toQ(qdot) = Rotation::convertAngVelToBodyFixed123Dot(fromQ(sbs.getQ()),
                                    ~R_FM*w_FM); // need w in *body*, not parent
    }

    void calcLocalQDotFromLocalU(const SBStateDigest& sbs, const Real* u, Real* qdot) const {
        assert(sbs.getStage() >= Stage::Position);
        assert(u && qdot);

        const SBModelVars&     mv   = sbs.getModelVars();
        const SBPositionCache& pc   = sbs.getPositionCache();
        const Vector&          allQ = sbs.getQ();

        const Vec3&            w_FM = Vec3::getAs(u);

        const Rotation& R_FM = getX_FM(pc).R();
        Vec3::updAs(qdot) = Rotation::convertAngVelToBodyFixed123Dot(fromQ(allQ),
                                         ~R_FM*w_FM); // need w in *body*, not parent
    }

    // Compute out_q = N * in_u
    //   or    out_u = in_q * N
    void multiplyByN(const SBStateDigest& sbs, bool useEulerAnglesIfPossible, const Real* q,
                          bool matrixOnRight, const Real* in, Real* out) const
    {
        assert(sbs.getStage() >= Stage::Model);
        assert(q && in && out);

        // TODO: it's annoying that this N block is only available in the Body (M) frame,
        // because this mobilizer uses angular velocity in the Parent (F) frame
        // as generalized speeds. So we have to do an expensive conversion here. It
        // would be just as easy to compute this matrix in the Parent frame in 
        // the first place.
        const Rotation R_FM(BodyRotationSequence, q[0], XAxis, q[1], YAxis, q[2], ZAxis);
        const Mat33    N = Rotation::calcQBlockForBodyXYZInBodyFrame(Vec3::getAs(q)) * ~R_FM;
        if (matrixOnRight) Row3::updAs(out) = Row3::getAs(in) * N;
        else               Vec3::updAs(out) = N * Vec3::getAs(in);
    }

    // Compute out_u = inv(N) * in_q
    //   or    out_q = in_u * inv(N)
    void multiplyByNInv(const SBStateDigest& sbs, bool useEulerAnglesIfPossible, const Real* q,
                             bool matrixOnRight, const Real* in, Real* out) const
    {
        assert(sbs.getStage() >= Stage::Position);
        assert(in && out);

        // TODO: see above regarding the need for this R_FM kludge
        const Rotation R_FM(BodyRotationSequence, q[0], XAxis, q[1], YAxis, q[2], ZAxis);
        const Mat33    NInv(R_FM*Rotation::calcQInvBlockForBodyXYZInBodyFrame(Vec3::getAs(q)));
        if (matrixOnRight) Row3::updAs(out) = Row3::getAs(in) * NInv;
        else               Vec3::updAs(out) = NInv * Vec3::getAs(in);
    }
 
    void calcQDotDot(
        const SBStateDigest&   sbs,
        const Vector&          udot, 
        Vector&                qdotdot) const 
    {
        const SBPositionCache& pc = sbs.getPositionCache();
        const Vec3& w_FM     = fromU(sbs.getU()); // angular velocity of J in Jb, expr in Jb
        const Vec3& w_FM_dot = fromU(udot);

        const Rotation& R_FM = getX_FM(pc).R();
        toQ(qdotdot)    = Rotation::convertAngVelDotToBodyFixed123DotDot
                              (fromQ(sbs.getQ()), ~R_FM*w_FM, ~R_FM*w_FM_dot);
    }

    void calcLocalQDotDotFromLocalUDot(const SBStateDigest& sbs, const Real* udot, Real* qdotdot) const {
        assert(sbs.getStage() >= Stage::Velocity);
        assert(udot && qdotdot);

        const SBModelVars&     mv   = sbs.getModelVars();
        const SBPositionCache& pc   = sbs.getPositionCache();
        const Vector&          allQ = sbs.getQ();
        const Vector&          allU = sbs.getU();

        const Vec3& w_FM     = fromU(allU);
        const Vec3& w_FM_dot = Vec3::getAs(udot);

        const Rotation& R_FM = getX_FM(pc).R();
        Vec3::updAs(qdotdot) = Rotation::convertAngVelDotToBodyFixed123DotDot
                                    (fromQ(allQ), ~R_FM*w_FM, ~R_FM*w_FM_dot);
    }

    void getInternalForce(const SBAccelerationCache&, Vector&) const {
        assert(false); // TODO: decompose cross-joint torque into 123 gimbal torques
        /* OLD BALL CODE:
        Vector& f = s.cache->netHingeForces;
        //dependency: calcR_PB must be called first
        assert( useEuler );

        const Vec<3,Vec2>& scq = getSinCosQ(s);
        const Real sPhi   = scq[0][0], cPhi   = scq[0][1];
        const Real sTheta = scq[1][0], cTheta = scq[1][1];
        const Real sPsi   = scq[2][0], cPsi   = scq[2][1];

        Vec3 torque = forceInternal;
        const Mat33 M( 0.          , 0.          , 1.    ,
                      -sPhi        , cPhi        , 0.    ,
                       cPhi*cTheta , sPhi*cTheta ,-sTheta );
        Vec3 eTorque = RigidBodyNode::DEG2RAD * M * torque;

        Vec3::updAs(&v[uIndex]) = eTorque;
        */
    }
};



#endif // SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_GIMBAL_H_

