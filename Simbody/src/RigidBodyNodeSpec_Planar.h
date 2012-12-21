#ifndef SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_PLANAR_H_
#define SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_PLANAR_H_

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
 * Define the RigidBodyNode that implements a Planar mobilizer.
 */

#include "SimbodyMatterSubsystemRep.h"
#include "RigidBodyNode.h"
#include "RigidBodyNodeSpec.h"


    // PLANAR //

// This provides free motion (translation and rotation) in a plane. We use
// the 2d coordinate system formed by the x,y axes of F as the translations,
// and the common z axis of F and M as the rotational axis. The generalized
// coordinates are theta,x,y interpreted as rotation around z and translation
// along the (space fixed) Fx and Fy axes.
template<bool noX_MB, bool noR_PF>
class RBNodePlanar : public RigidBodyNodeSpec<3, false, noX_MB, noR_PF> {
public:
typedef typename RigidBodyNodeSpec<3, false, noX_MB, noR_PF>::HType HType;
virtual const char* type() { return "planar"; }

RBNodePlanar(const MassProperties&    mProps_B,
                const Transform&      X_PF,
                const Transform&      X_BM,
                bool                  isReversed,
                UIndex&               nextUSlot,
                USquaredIndex&        nextUSqSlot,
                QIndex&               nextQSlot)
:   RigidBodyNodeSpec<3, false, noX_MB, noR_PF>(mProps_B,X_PF,X_BM,nextUSlot,nextUSqSlot,nextQSlot,
                         RigidBodyNode::QDotIsAlwaysTheSameAsU, RigidBodyNode::QuaternionIsNeverUsed, 
                         isReversed)
{
    this->updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
}

    // Implementations of virtual methods.

// This has a default implementation but it rotates first then translates,
// which works fine for the normal Planar joint but produces wrong behavior
// when the mobilizer is reversed.
void setQToFitTransformImpl(const SBStateDigest& sbs, const Transform& X_FM, 
                            Vector& q) const OVERRIDE_11 
{
    setQToFitTranslationImpl(sbs, X_FM.p(), q); // see below
    setQToFitRotationImpl(sbs, X_FM.R(), q);
}

void setQToFitRotationImpl(const SBStateDigest& sbs, const Rotation& R_FM, Vector& q) const {
    // The only rotation our planar joint can handle is about z.
    // TODO: should use 321 to deal with singular configuration (angle2==pi/2) better;
    // in that case 1 and 3 are aligned and the conversion routine allocates all the
    // rotation to whichever comes first.
    // TODO: isn't there a better way to come up with "the rotation around z that
    // best approximates a rotation R"?
    const Vec3 angles123 = R_FM.convertRotationToBodyFixedXYZ();
    this->toQ(q)[0] = angles123[2];
}
void setQToFitTranslationImpl(const SBStateDigest& sbs, const Vec3&  p_FM, Vector& q) const {
    // Ignore translation in the z direction.
    this->toQ(q)[1] = p_FM[0]; // x
    this->toQ(q)[2] = p_FM[1]; // y
}

void setUToFitAngularVelocityImpl(const SBStateDigest& sbs, const Vector&, const Vec3& w_FM, Vector& u) const {
    // We can represent the z angular velocity exactly, but nothing else.
    this->toU(u)[0] = w_FM[2];
}
void setUToFitLinearVelocityImpl
    (const SBStateDigest& sbs, const Vector&, const Vec3& v_FM, Vector& u) const
{
    // Ignore translational velocity in the z direction.
    this->toU(u)[1] = v_FM[0]; // x
    this->toU(u)[2] = v_FM[1]; // y
}

enum {PoolSize=2};
enum {CosQ=0, SinQ=1};
// We want space for cos(q) and sin(q).
int calcQPoolSize(const SBModelVars&) const
{   return PoolSize; }

// Precalculation of sin/cos costs around 50 flops.
void performQPrecalculations(const SBStateDigest& sbs,
                             const Real* q, int nq,
                             Real* qCache,  int nQCache,
                             Real* qErr,    int nQErr) const
{
    assert(q && nq==3 && qCache && nQCache==PoolSize && nQErr==0);
    qCache[CosQ] = std::cos(q[0]);
    qCache[SinQ] = std::sin(q[0]);
}

// This is nearly free since we already calculated sin/cos.
void calcX_FM(const SBStateDigest& sbs,
                const Real* q,      int nq,
                const Real* qCache, int nQCache,
                Transform&  X_FM) const
{
    assert(q && nq==3 && qCache && nQCache==PoolSize);
    // Rotational q is about common z axis, translational q's along 
    // Fx and Fy.
    X_FM.updR().setRotationFromAngleAboutZ(qCache[CosQ], qCache[SinQ]);
    X_FM.updP() = Vec3(q[1], q[2], 0);
}

// The rotational generalized speed is about the common z axis; translations
// are along Fx and Fy so all axes are constant in F.
void calcAcrossJointVelocityJacobian(
    const SBStateDigest& sbs,
    HType&               H_FM) const
{
    H_FM(0) = SpatialVec( Vec3(0,0,1),   Vec3(0) );
    H_FM(1) = SpatialVec(   Vec3(0),   Vec3(1,0,0) );
    H_FM(2) = SpatialVec(   Vec3(0),   Vec3(0,1,0) );
}

// Since the Jacobian above is constant in F, its time derivative is zero.
void calcAcrossJointVelocityJacobianDot(
    const SBStateDigest& sbs,
    HType&               HDot_FM) const
{
    HDot_FM(0) = SpatialVec( Vec3(0), Vec3(0) );
    HDot_FM(1) = SpatialVec( Vec3(0), Vec3(0) );
    HDot_FM(2) = SpatialVec( Vec3(0), Vec3(0) );
}

};


#endif // SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_PLANAR_H_

