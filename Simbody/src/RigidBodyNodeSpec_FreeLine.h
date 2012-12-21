#ifndef SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_FREELINE_H_
#define SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_FREELINE_H_

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
 * Contributors: Paul Mitiguy                                                 *
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
 * Define the RigidBodyNode that implements a Translation mobilizer, also
 * known as a Cartesian joint.
 */

#include "SimbodyMatterSubsystemRep.h"
#include "RigidBodyNode.h"
#include "RigidBodyNodeSpec.h"


    // FREE LINE //

// FreeLine joint. Like a Free joint, this provides full rotational and
// translational freedom, but for a degenerate body which is thin (inertialess)
// along its own z axis. These arise in molecular modeling for linear molecules 
// formed by pairs of atoms, or by multiple atoms in a linear arrangement like
// carbon dioxide (CO2) whose structure is O=C=O in a straight line. We are
// assuming that there is no meaning to a rotation about the linear axis,
// so free orientation requires just *two* degrees of freedom, not *three*
// as is required for general rigid bodies. And in fact we can get away with
// just two rotational generalized speeds so this joint provides only 5 
// mobilities.
//
// But so far, no one has been able to come up with a way to manage with only 
// two rotational generalized *coordinates*, so this joint has the same q's as 
// a regular Free joint: either a quaternion for unconditional stability, or a 
// three-angle (body fixed 1-2-3) Euler sequence which will be dynamically 
// singular when the middle (y) axis is 90 degrees. Use the Euler sequence only 
// for small motions or for kinematics problems (and note that only the first 
// two are meaningful). Translations here are treated exactly as for a Free 
// joint (or for a Cartesian joint for that matter).
//
// To summarize, the generalized coordinates are:
//   * 4 quaternions or 3 1-2-3 body fixed Euler angles (that is, fixed in M)
//   * 3 components of the translation vector p_FM (that is, vector from origin
//     of F to origin of M, expressed in F)
// and generalized speeds are:
//   * the x,y components of the angular velocity w_FM_M, that is, the angular
//     velocity of M in F expressed in *M* (where we want wz=0).
//   * 3 components of the linear velocity of origin of M in F, expressed in F.
//     NOTE: THAT IS NOT THE SAME FRAME AS FOR A FREE JOINT
// Thus the qdots have to be derived from the generalized speeds to
// be turned into either 4 quaternion derivatives or 3 Euler angle derivatives.
template<bool noX_MB, bool noR_PF>
class RBNodeFreeLine : public RigidBodyNodeSpec<5, false, noX_MB, noR_PF> {
public:

typedef typename RigidBodyNodeSpec<5, false, noX_MB, noR_PF>::HType HType;
virtual const char* type() { return "full"; }

RBNodeFreeLine(const MassProperties& mProps_B,
               const Transform&      X_PF,
               const Transform&      X_BM,
               bool                  isReversed,
               UIndex&               nextUSlot,
               USquaredIndex&        nextUSqSlot,
               QIndex&               nextQSlot)
  : RigidBodyNodeSpec<5, false, noX_MB, noR_PF>(mProps_B,X_PF,X_BM,nextUSlot,nextUSqSlot,nextQSlot,
                         RigidBodyNode::QDotMayDifferFromU, RigidBodyNode::QuaternionMayBeUsed, isReversed)
{
    this->updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
}

    // Implementations of virtual methods.

// This has a default implementation but it rotates first then translates,
// which works fine for the normal FreeLine joint but produces wrong behavior
// when the mobilizer is reversed.
void setQToFitTransformImpl(const SBStateDigest& sbs, const Transform& X_FM, 
                            Vector& q) const OVERRIDE_11 
{
    setQToFitTranslationImpl(sbs, X_FM.p(), q); // see below
    setQToFitRotationImpl(sbs, X_FM.R(), q);
}

void setQToFitRotationImpl(const SBStateDigest& sbs, const Rotation& R_FM,
                          Vector& q) const 
{
    if (this->getUseEulerAngles(sbs.getModelVars()))
        this->toQVec3(q,0) = R_FM.convertRotationToBodyFixedXYZ();
    else
        this->toQuat(q) = R_FM.convertRotationToQuaternion().asVec4();
}

// The user gives us the translation vector from OF to OM as a vector expressed 
// in F. With a free joint we never have to *change* orientation coordinates in 
// order to achieve a translation.
void setQToFitTranslationImpl(const SBStateDigest& sbs, const Vec3& p_FM, Vector& q) const {
    if (this->getUseEulerAngles(sbs.getModelVars()))
        this->toQVec3(q,3) = p_FM; // skip the 3 Euler angles
    else
        this->toQVec3(q,4) = p_FM; // skip the 4 quaternions
}

// Our 2 rotational generalized speeds are just the (x,y) components of the
// angular velocity vector of M in F, expressed in *M*.
// Note: a quaternion from a state is not necessarily normalized so can't be 
// used directly as though it were a set of Euler parameters; it must be 
// normalized first.
void setUToFitAngularVelocityImpl(const SBStateDigest& sbs, 
                                  const Vector& q, const Vec3& w_FM,
                                  Vector& u) const
{
    Rotation R_FM;
    if (this->getUseEulerAngles(sbs.getModelVars()))
        R_FM.setRotationToBodyFixedXYZ( this->fromQVec3(q,0) );
    else {
        // can't use qnorm pool here since state hasn't been 
        // realized to position stage yet; q's can be anything
        R_FM.setRotationFromQuaternion( Quaternion(this->fromQuat(q)) ); // normalize
    }
    const Vec3 w_FM_M = ~R_FM*w_FM;
    // (x,y) of relative angular velocity always used as generalized speeds
    Vec2::updAs(&this->toU(u)[0]) = Vec2(w_FM_M[0], w_FM_M[1]); 
}

// Our 3 translational generalized speeds are the linear velocity of M's origin
// in F, expressed in F. The user gives us that same vector.
void setUToFitLinearVelocityImpl
   (const SBStateDigest& sbs, const Vector& q, const Vec3& v_FM, Vector& u) const
{
    this->toUVec3(u,2) = v_FM;
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
{   return this->getUseEulerAngles(mv) ? AnglePoolSize : QuatPoolSize; }

void performQPrecalculations(const SBStateDigest& sbs,
                             const Real* q,      int nq,
                             Real*       qCache, int nQCache,
                             Real*       qErr,   int nQErr) const
{
    if (this->getUseEulerAngles(sbs.getModelVars())) {
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
    if (this->getUseEulerAngles(sbs.getModelVars())) {
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


// The generalized speeds for this 5-dof ("free line") joint are 
//   (1) the (x,y) components of angular velocity of M in the F frame, 
//         expressed in *M*, and
//   (2) the (linear) velocity of M's origin in F, expressed in F.
void calcAcrossJointVelocityJacobian(
    const SBStateDigest& sbs,
    HType&               H_FM) const
{
    const SBTreePositionCache& pc = sbs.updTreePositionCache(); // "upd" because we're realizing positions now
    const Transform  X_F0M0 = this->findX_F0M0(pc);

    // Dropping the 0's here.
    const Rotation& R_FM = X_F0M0.R();
    const Vec3&     Mx_F = R_FM.x(); // M's x axis, expressed in F
    const Vec3&     My_F = R_FM.y(); // M's y axis, expressed in F

    H_FM(0) = SpatialVec( Mx_F, Vec3(0) );  // x,y ang. vel. in M, exp. in F
    H_FM(1) = SpatialVec( My_F, Vec3(0) );

    H_FM(2) = SpatialVec( Vec3(0), Vec3(1,0,0) );   // translations in F
    H_FM(3) = SpatialVec( Vec3(0), Vec3(0,1,0) );
    H_FM(4) = SpatialVec( Vec3(0), Vec3(0,0,1) );
}

// Since the first two rows of the Jacobian above are not constant in F,
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
    const Transform  X_F0M0 = this->findX_F0M0(pc);

    // Dropping the 0's here.
    const Rotation& R_FM = X_F0M0.R();
    const Vec3&     Mx_F = R_FM.x(); // M's x axis, expressed in F
    const Vec3&     My_F = R_FM.y(); // M's y axis, expressed in F

    const Vec3      w_FM = this->find_w_F0M0(pc,vc); // angular velocity of M in F

    HDot_FM(0) = SpatialVec( w_FM % Mx_F, Vec3(0) );
    HDot_FM(1) = SpatialVec( w_FM % My_F, Vec3(0) );

    // For translation in F.
    HDot_FM(2) = SpatialVec( Vec3(0), Vec3(0) );
    HDot_FM(3) = SpatialVec( Vec3(0), Vec3(0) );
    HDot_FM(4) = SpatialVec( Vec3(0), Vec3(0) );
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

    if (this->getUseEulerAngles(mv)) {
        const Vec3& q = this->fromQVec3(allQ,0);
        const Mat32 N_FM = // N here mixes parent and body frames (convenient)
            Rotation::calcNForBodyXYZInBodyFrame(q)
                        .getSubMat<3,2>(0,0); // drop 3rd column
        // Translational part of N is identity so just copy in to out.
        if (matrixOnRight) {
            Row2::updAs(out)   = Row3::getAs(in) * N_FM;
            Row3::updAs(out+2) = Row3::getAs(in+3);// translation
        } else {
            Vec3::updAs(out)   = N_FM * Vec2::getAs(in);        
            Vec3::updAs(out+3) = Vec3::getAs(in+2);// translation
        }
    } else {
        // Quaternion: N block is only available expecting angular velocity 
        // in the parent frame F, but we have it in M for this joint so we
        // have to calculate N_FM = N_FF*R_FM.
        const Rotation& R_FM = this->getX_FM(pc).R();
        const Vec4& q = this->fromQuat(allQ);
        const Mat42 N_FM = (Rotation::calcUnnormalizedNForQuaternion(q)*R_FM)
                            .getSubMat<4,2>(0,0); // drop 3rd column
        // Translational part of N is identity so just copy in to out.
        if (matrixOnRight) {
            Row2::updAs(out)   = Row4::getAs(in) * N_FM;
            Row3::updAs(out+2) = Row3::getAs(in+4); // translation
        } else { // matrix on left
            Vec4::updAs(out)   = N_FM * Vec2::getAs(in);
            Vec3::updAs(out+4) = Vec3::getAs(in+2); // translation
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

    if (this->getUseEulerAngles(mv)) {
        const Vec3& q = this->fromQVec3(allQ,0);
        const Mat23 NInv_MF = Rotation::calcNInvForBodyXYZInBodyFrame(q)
                                .getSubMat<2,3>(0,0);   // drop 3rd row
        // Translational part of NInv block is identity.
        if (matrixOnRight) {
            Row3::updAs(out)   = Row2::getAs(in) * NInv_MF;
            Row3::updAs(out+3) = Row3::getAs(in+2); // translation
        } else {
            Vec2::updAs(out)   = NInv_MF * Vec3::getAs(in);
            Vec3::updAs(out+2) = Vec3::getAs(in+3); // translation
        }
    } else {           
        // Quaternion: NInv block is only available expecting angular 
        // velocity in the parent frame F, but we have it in M for 
        // this joint so we have to calculate NInv_MF = R_MF*NInv_FF.
        const Rotation& R_FM = this->getX_FM(pc).R();
        const Vec4& q = this->fromQuat(allQ);
        const Mat24 NInv_MF = (~R_FM*Rotation::calcUnnormalizedNInvForQuaternion(q))
                                .getSubMat<2,4>(0,0);   // drop 3rd row
        // Translational part of NInv block is identity.
        if (matrixOnRight) {
            Row4::updAs(out)   = Row2::getAs(in) * NInv_MF;
            Row3::updAs(out+4) = Row3::getAs(in+2); // translation
        } else { // matrix on left
            Vec2::updAs(out)   = NInv_MF * Vec4::getAs(in);
            Vec3::updAs(out+2) = Vec3::getAs(in+4); // translation
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

    if (this->getUseEulerAngles(mv)) {
        const Vec3& q = this->fromQVec3(allQ,0);
        const Vec3& qdot = this->fromQVec3(allQDot,0);
        const Mat32 NDot_FM = // NDot here mixes parent and body frames (convenient)
            Rotation::calcNDotForBodyXYZInBodyFrame(q, qdot)
                        .getSubMat<3,2>(0,0); // drop 3rd column
        // Translational part of NDot is zero so set out to zero.
        if (matrixOnRight) {
            Row2::updAs(out)   = Row3::getAs(in) * NDot_FM;
            Row3::updAs(out+2) = 0; // translation
        } else {
            Vec3::updAs(out)   = NDot_FM * Vec2::getAs(in);        
            Vec3::updAs(out+3) = 0; // translation
        }
    } else {
        // Quaternion: NDot block is only available expecting angular velocity 
        // in the parent frame F, but we have it in M for this joint so we
        // have to calculate NDot_FM = NDot_FF*R_FM.
        const Rotation& R_FM = this->getX_FM(pc).R();
        const Vec4& qdot = this->fromQuat(allQDot);
        const Mat42 NDot_FM = 
           (Rotation::calcUnnormalizedNDotForQuaternion(qdot)*R_FM)
                        .getSubMat<4,2>(0,0); // drop 3rd column
        // Translational part of N is identity so just copy in to out.
        if (matrixOnRight) {
            Row2::updAs(out)   = Row4::getAs(in) * NDot_FM;
            Row3::updAs(out+2) = 0; // translation
        } else { // matrix on left
            Vec4::updAs(out)   = NDot_FM * Vec2::getAs(in);
            Vec3::updAs(out+4) = 0; // translation
        }
    }
}


void calcQDot(
    const SBStateDigest&   sbs,
    const Real*            u,
    Real*                  qdot) const
{
    const SBModelVars&          mv = sbs.getModelVars();
    const SBTreePositionCache&  pc = sbs.getTreePositionCache();
    const Vec3  w_FM_M = Vec3(u[0], u[1], 0); // Angular velocity in M
    const Vec3& v_FM   = Vec3::getAs(&u[2]);  // Linear velocity in F

    if (this->getUseEulerAngles(mv)) {
        const Vec3& theta = this->fromQVec3(sbs.getQ(),0); // Euler angles
        Vec3::updAs(qdot) = Rotation::convertAngVelInBodyFrameToBodyXYZDot(theta,
                                        w_FM_M); // need w in *body*, not parent
        Vec3::updAs(&qdot[3]) = v_FM;
        qdot[6] = 0;
    } else {
        const Rotation& R_FM = this->getX_FM(pc).R();
        const Vec4& quat = this->fromQuat(sbs.getQ());
        Vec4::updAs(qdot)   = Rotation::convertAngVelToQuaternionDot(quat,
                                        R_FM*w_FM_M); // need w in *parent* frame here
        Vec3::updAs(&qdot[4]) = v_FM;
    }
}

void calcQDotDot(
    const SBStateDigest&   sbs,
    const Real*            udot, 
    Real*                  qdotdot) const 
{
    const SBModelVars&          mv = sbs.getModelVars();
    const SBTreePositionCache&  pc = sbs.getTreePositionCache();
    const Vec3  w_FM_M     = Vec3(this->fromU(sbs.getU())[0], this->fromU(sbs.getU())[1], 0); // Angular velocity of M in F, exp. in M
    const Vec3& v_FM       = this->fromUVec3(sbs.getU(),2); // linear velocity of M in F, expressed in M
    const Vec3  w_FM_M_dot = Vec3(udot[0], udot[1], 0);
    const Vec3& v_FM_dot   = Vec3::getAs(&udot[2]);

    if (this->getUseEulerAngles(mv)) {
        const Vec3& theta  = this->fromQVec3(sbs.getQ(),0); // Euler angles
        Vec3::updAs(qdotdot) = Rotation::convertAngVelDotInBodyFrameToBodyXYZDotDot
                                         (theta, w_FM_M, w_FM_M_dot); // needed in body frame here
        Vec3::updAs(&qdotdot[3]) = v_FM_dot;
        qdotdot[6] = 0;
    } else {
        const Rotation& R_FM = this->getX_FM(pc).R();
        const Vec4& quat  = this->fromQuat(sbs.getQ());
        Vec4::updAs(qdotdot)   = Rotation::convertAngVelDotToQuaternionDotDot
                                         (quat,R_FM*w_FM_M,R_FM*w_FM_M_dot); // needed in parent frame
        Vec3::updAs(&qdotdot[4]) = v_FM_dot;
    }
}

int  getMaxNQ()                   const {return 7;}
int  getNQInUse(const SBModelVars& mv) const {return this->getUseEulerAngles(mv) ? 6 : 7;} 
bool isUsingQuaternion(const SBStateDigest& sbs, MobilizerQIndex& startOfQuaternion) const {
    if (this->getUseEulerAngles(sbs.getModelVars())) {startOfQuaternion.invalidate(); return false;}
    startOfQuaternion = MobilizerQIndex(0); // quaternion comes first
    return true;
}

void setMobilizerDefaultPositionValues(const SBModelVars& mv, Vector& q) const 
{
    if (this->getUseEulerAngles(mv)) {
        this->toQVec3(q,4) = Vec3(0); // TODO: kludge, clear unused element
        this->toQ(q) = 0.;
    } else {
        this->toQuat(q) = Vec4(1.,0.,0.,0.);
        this->toQVec3(q,4) = 0.;
    }
}

bool enforceQuaternionConstraints(
    const SBStateDigest& sbs, 
    Vector&             q,
    Vector&             qErrest) const 
{
    if (this->getUseEulerAngles(sbs.getModelVars())) 
        return false; // no change

    Vec4& quat = this->toQuat(q);
    quat = quat / quat.norm();

    if (qErrest.size()) {
        Vec4& qerr = this->toQuat(qErrest);
        qerr -= dot(qerr,quat) * quat;
    }

    return true;
}

void convertToEulerAngles(const Vector& inputQ, Vector& outputQ) const {
    this->toQVec3(outputQ, 4) = Vec3(0); // clear unused element
    this->toQVec3(outputQ, 3) = this->fromQVec3(inputQ, 4);
    this->toQVec3(outputQ, 0) = Rotation(Quaternion(this->fromQuat(inputQ))).convertRotationToBodyFixedXYZ();
}

void convertToQuaternions(const Vector& inputQ, Vector& outputQ) const {
    this->toQVec3(outputQ, 4) = this->fromQVec3(inputQ, 3);
    Rotation rot;
    rot.setRotationToBodyFixedXYZ(this->fromQVec3(inputQ, 0));
    this->toQuat(outputQ) = rot.convertRotationToQuaternion().asVec4();
}


};




#endif // SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_FREELINE_H_

