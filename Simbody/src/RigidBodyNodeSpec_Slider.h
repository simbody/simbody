#ifndef SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_SLIDER_H_
#define SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_SLIDER_H_

/* -------------------------------------------------------------------------- *
 *                              SimTK Simbody(tm)                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-11 Stanford University and the Authors.        *
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
class RBNodeSlider : public RigidBodyNodeSpec<1, true> {
public:
virtual const char* type() { return "slider"; }

RBNodeSlider(const MassProperties&    mProps_B,
                const Transform&      X_PF,
                const Transform&      X_BM,
                bool                  isReversed,
                UIndex&               nextUSlot,
                USquaredIndex&        nextUSqSlot,
                QIndex&               nextQSlot)
:   RigidBodyNodeSpec<1, true>(mProps_B,X_PF,X_BM,nextUSlot,nextUSqSlot,nextQSlot,
                         QDotIsAlwaysTheSameAsU, QuaternionIsNeverUsed, 
                         isReversed)
{
    updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
}

    // Implementations of virtual methods.

void setQToFitRotationImpl(const SBStateDigest& sbs, const Rotation& R_FM, 
                           Vector& q) const {
    // The only rotation a slider can represent is identity.
}

void setQToFitTranslationImpl(const SBStateDigest& sbs, const Vec3& p_FM, 
                              Vector& q) const {
    // We can only represent the x coordinate with this joint.
    to1Q(q) = p_FM[0];
}

void setUToFitAngularVelocityImpl(const SBStateDigest& sbs, const Vector&, 
                                  const Vec3& w_FM, Vector& u) const {
    // The only angular velocity a slider can represent is zero.
}

void setUToFitLinearVelocityImpl(const SBStateDigest& sbs, const Vector&, 
                                 const Vec3& v_FM, Vector& u) const
{
    // We can only represent a velocity along x with this joint.
    to1U(u) = v_FM[0];
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

