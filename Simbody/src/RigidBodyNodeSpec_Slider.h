#ifndef SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_SLIDER_H_
#define SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_SLIDER_H_

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
 * Define the RigidBodyNode that implements a Slider mobilizer, also known as a
 * Sliding or Prismatic joint.
 */

#include "SimbodyMatterSubsystemRep.h"
#include "RigidBodyNode.h"
#include "RigidBodyNodeSpec.h"

    // SLIDING (PRISMATIC) //

// Sliding joint (1 dof translation). The translation is along the x
// axis of the parent body's F frame, with M=F when the coordinate
// is zero and the orientation of M in F frozen at 0 forever.
template<bool noX_MB, bool noR_PF>
class RBNodeSlider : public RigidBodyNodeSpec<1, true, noX_MB, noR_PF> {
public:
typedef typename RigidBodyNodeSpec<1, false, noX_MB, noR_PF>::HType HType;
virtual const char* type() { return "slider"; }

RBNodeSlider(const MassProperties&    mProps_B,
                const Transform&      X_PF,
                const Transform&      X_BM,
                bool                  isReversed,
                UIndex&               nextUSlot,
                USquaredIndex&        nextUSqSlot,
                QIndex&               nextQSlot)
:   RigidBodyNodeSpec<1, true, noX_MB, noR_PF>(mProps_B,X_PF,X_BM,nextUSlot,nextUSqSlot,nextQSlot,
                         RigidBodyNode::QDotIsAlwaysTheSameAsU, RigidBodyNode::QuaternionIsNeverUsed, 
                         isReversed)
{
    this->updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
}

    // Implementations of virtual methods.

void setQToFitRotationImpl(const SBStateDigest& sbs, const Rotation& R_FM, 
                           Vector& q) const {
    // The only rotation a slider can represent is identity.
}

void setQToFitTranslationImpl(const SBStateDigest& sbs, const Vec3& p_FM, 
                              Vector& q) const {
    // We can only represent the x coordinate with this joint.
    this->to1Q(q) = p_FM[0];
}

void setUToFitAngularVelocityImpl(const SBStateDigest& sbs, const Vector&, 
                                  const Vec3& w_FM, Vector& u) const {
    // The only angular velocity a slider can represent is zero.
}

void setUToFitLinearVelocityImpl(const SBStateDigest& sbs, const Vector&, 
                                 const Vec3& v_FM, Vector& u) const
{
    // We can only represent a velocity along x with this joint.
    this->to1U(u) = v_FM[0];
}

// A sliding joint doesn't need to cache any q calculations.
int calcQPoolSize(const SBModelVars&) const
{   return 0; }

// Nothing to precalculate.
void performQPrecalculations(const SBStateDigest& sbs,
                                const Real* q, int nq,
                                Real* qCache,  int nQCache,
                                Real* qErr,    int nQErr) const
{
    assert(q && nq==1 && nQCache==0 && nQErr==0);
}

// This is free.
void calcX_FM(const SBStateDigest& sbs,
                const Real* q,      int nq,
                const Real* qCache, int nQCache,
                Transform&  X_FM) const
{
    assert(q && nq==1 && nQCache==0);
    // Translation vector q is expressed in F (and M since they have same 
    // orientation). A sliding joint can't change orientation, and only 
    // translates along x. 
    X_FM = Transform(Rotation(), Vec3(q[0],0,0));
}

// The generalized speed is the velocity of M's origin in the F frame,
// along F's x axis, expressed in F.
void calcAcrossJointVelocityJacobian(
    const SBStateDigest& sbs,
    HType&               H_FM) const
{
    H_FM(0) = SpatialVec( Vec3(0), Vec3(1,0,0) );
}

// Since the Jacobian above is constant in F, its time derivative is zero.
void calcAcrossJointVelocityJacobianDot(
    const SBStateDigest& sbs,
    HType&               HDot_FM) const
{
    HDot_FM(0) = SpatialVec( Vec3(0), Vec3(0) );
}

// Override the computation of reverse-H for this simple mobilizer.
void calcReverseMobilizerH_FM(
    const SBStateDigest& sbs,
    HType&               H_FM) const
{
    H_FM(0) = SpatialVec( Vec3(0), Vec3(-1,0,0) );
}

// Override the computation of reverse-HDot for this simple mobilizer.
void calcReverseMobilizerHDot_FM(
    const SBStateDigest& sbs,
    HType&               HDot_FM) const
{
    HDot_FM(0) = SpatialVec( Vec3(0), Vec3(0) );
}

};


#endif // SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_SLIDER_H_

