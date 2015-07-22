#ifndef SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_POLARCOORDS_H_
#define SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_POLARCOORDS_H_

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
 * Contributors: Derived from IVM code written by Charles Schwieters          *
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
 * Define the RigidBodyNode that implements a PolarCoords mobilizer, also known
 * as a BendStretch joint.
 */

#include "SimbodyMatterSubsystemRep.h"
#include "RigidBodyNode.h"
#include "RigidBodyNodeSpec.h"


    // BEND-STRETCH //

// This is a "bend-stretch" joint, meaning one degree of rotational freedom
// about a particular axis, and one degree of translational freedom
// along a perpendicular axis. The z axis of the parent's F frame is
// used for rotation (and that is always aligned with the M frame z axis).
// The x axis of the *M* frame is used for translation; that is, first
// we rotate around z, which moves M's x with respect to F's x. Then
// we slide along the rotated x axis. The two generalized coordinates are the
// rotation and the translation, in that order.
template<bool noX_MB, bool noR_PF>
class RBNodeBendStretch : public RigidBodyNodeSpec<2, false, noX_MB, noR_PF> {
public:
typedef typename RigidBodyNodeSpec<2, false, noX_MB, noR_PF>::HType HType;
virtual const char* type() { return "bendstretch"; }

RBNodeBendStretch(const MassProperties&   mProps_B,
                    const Transform&      X_PF,
                    const Transform&      X_BM,
                    bool                  isReversed,
                    UIndex&               nextUSlot,
                    USquaredIndex&        nextUSqSlot,
                    QIndex&               nextQSlot)
:   RigidBodyNodeSpec<2, false, noX_MB, noR_PF>(mProps_B,X_PF,X_BM,nextUSlot,nextUSqSlot,nextQSlot,
                         RigidBodyNode::QDotIsAlwaysTheSameAsU, RigidBodyNode::QuaternionIsNeverUsed,
                         isReversed)
{
    this->updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
}


void setQToFitRotationImpl(const SBStateDigest& sbs, const Rotation& R_FM,
                           Vector& q) const {
    // The only rotation our bend-stretch joint can handle is about z.
    // TODO: this code is bad -- see comments for Torsion joint above.
    const Vec3 angles123 = R_FM.convertRotationToBodyFixedXYZ();
    this->toQ(q)[0] = angles123[2];
}

void setQToFitTranslationImpl(const SBStateDigest& sbs, const Vec3& p_FM,
                              Vector& q) const {
    // We can represent any translation that puts the M origin in the x-y plane
    // of F, by a suitable rotation around z followed by translation along x.
    const Vec2 r = p_FM.getSubVec<2>(0); // (rx, ry)
    const Real d = r.norm();

    // If there is no translation worth mentioning, we'll leave the rotational
    // coordinate alone, otherwise rotate so M's x axis is aligned with r.
    if (d >= 4*Eps) {
        const Real angle = std::atan2(r[1],r[0]);
        this->toQ(q)[0] = angle;
        this->toQ(q)[1] = d;
    } else
        this->toQ(q)[1] = 0;
}

void setUToFitAngularVelocityImpl(const SBStateDigest& sbs, const Vector& q,
                                  const Vec3& w_FM, Vector& u) const {
    // We can only represent an angular velocity along z with this joint.
    this->toU(u)[0] = w_FM[2];
}

// If the translational coordinate is zero, we can only represent a linear
// velocity of Mo in F which is along M's current x axis direction. Otherwise,
// we can represent any velocity in the x-y plane by introducing angular
// velocity about z. We can never represent a linear velocity along z.
void setUToFitLinearVelocityImpl
   (const SBStateDigest& sbs, const Vector& q, const Vec3& v_FM, Vector& u) const
{
    // Decompose the requested v into "along Mx" and "along My" components.
    const Rotation R_FM = Rotation( this->fromQ(q)[0], ZAxis ); // =[ Mx My Mz ] in F
    const Vec3 v_FM_M = ~R_FM*v_FM; // re-express in M frame

    this->toU(u)[1] = v_FM_M[0]; // velocity along Mx we can represent directly

    const Real x = this->fromQ(q)[1]; // translation along Mx (signed)
    if (std::abs(x) < SignificantReal) {
        // No translation worth mentioning; we can only do x velocity, which
        // we just set above.
        return;
    }

    // significant translation
    this->toU(u)[0] = v_FM_M[1] / x; // set angular velocity about z to produce vy
}

enum {PoolSize=2}; // number of Reals
enum {CosQ=0, SinQ=1};
// We want space for cos(q0) and sin(q0).
int calcQPoolSize(const SBModelVars&) const
{   return PoolSize; }

// Precalculation of sin/cos costs around 50 flops.
void performQPrecalculations(const SBStateDigest& sbs,
                             const Real* q, int nq,
                             Real* qCache,  int nQCache,
                             Real* qErr,    int nQErr) const
{
    assert(q && nq==2 && qCache && nQCache==PoolSize && nQErr==0);
    qCache[CosQ] = std::cos(q[0]);
    qCache[SinQ] = std::sin(q[0]);
}

// This is about 15 flops.
void calcX_FM(const SBStateDigest& sbs,
                const Real* q,      int nq,
                const Real* qCache, int nQCache,
                Transform&  X_FM) const
{
    assert(q && nq==2 && qCache && nQCache==PoolSize);
    X_FM.updR().setRotationFromAngleAboutZ(qCache[CosQ], qCache[SinQ]);
    // Translational coordinate is in M frame so must be reexpressed in F.
    X_FM.updP() = X_FM.R()*Vec3(q[1],0,0); // 15 flops; could be cheaper
}

// The generalized speeds for this bend-stretch joint are (1) the angular
// velocity of M in the F frame, about F's z axis, expressed in F, and
// (2) the (linear) velocity of M's origin in F, along *M's* current x axis
// (that is, after rotation about z). (The z axis is also constant in M for
// this joint.)
// 9 flops (could be cheaper but who cares?)
void calcAcrossJointVelocityJacobian(
    const SBStateDigest& sbs,
    HType&               H_FM) const
{
    // use "upd" because we're realizing positions now
    const SBTreePositionCache& pc = sbs.updTreePositionCache();
    const Transform X_F0M0 = this->findX_F0M0(pc);
    const Rotation& R_F0M0 = X_F0M0.R();

    // Dropping the 0's here.
    const Vec3&     p_FM = X_F0M0.p();
    const Vec3&     Mx_F = R_F0M0.x(); // M's x axis, expressed in F

    H_FM(0) = SpatialVec( Vec3(0,0,1), Vec3(0,0,1) % p_FM   );
    H_FM(1) = SpatialVec( Vec3(0),            Mx_F          );
}

// Since the the Jacobian above is not constant in F,
// its time derivative is non zero. Here we use the fact that for
// a vector r_B_A fixed in a moving frame B but expressed in another frame A,
// its time derivative in A is the angular velocity of B in A crossed with
// the vector, i.e., d_A/dt r_B_A = w_AB % r_B_A.
// 18 flops (could be cheaper).
void calcAcrossJointVelocityJacobianDot(
    const SBStateDigest& sbs,
    HType&               HDot_FM) const
{
    const SBTreePositionCache& pc = sbs.getTreePositionCache();
    const SBTreeVelocityCache& vc = sbs.updTreeVelocityCache(); // use "upd" because we're realizing velocities now

    const Transform  X_F0M0 = this->findX_F0M0(pc);
    const Rotation&  R_F0M0 = X_F0M0.R();
    const SpatialVec V_F0M0 = this->findV_F0M0(pc,vc);

    // Dropping the 0's here.
    const Vec3&     Mx_F = R_F0M0.x(); // M's x axis, expressed in F
    const Vec3&     w_FM = V_F0M0[0]; // angular velocity of M in F
    const Vec3&     v_FM = V_F0M0[1]; // linear velocity of OM in F

    HDot_FM(0) = SpatialVec( Vec3(0), Vec3(0,0,1) % v_FM );
    HDot_FM(1) = SpatialVec( Vec3(0),    w_FM % Mx_F );
}

};


#endif // SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_POLARCOORDS_H_

