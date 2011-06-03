#ifndef SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_TRANSLATION_H_
#define SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_TRANSLATION_H_

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
 * Define the RigidBodyNode that implements a Translation mobilizer, also
 * known as a Cartesian joint.
 */

#include "SimbodyMatterSubsystemRep.h"
#include "RigidBodyNode.h"
#include "RigidBodyNodeSpec.h"

    // TRANSLATION (CARTESIAN) //

// Translate (Cartesian) joint. This provides three degrees of translational 
// freedom which is suitable (e.g.) for connecting a free atom to ground. The 
// Cartesian directions are the axes of the parent body's F frame, with M=F 
// when all 3 coords are 0, and the orientation of M in F is 0 (identity) 
// forever.
template<bool noX_MB, bool noR_PF>
class RBNodeTranslate : public RigidBodyNodeSpec<3, true, noX_MB, noR_PF> {
public:
typedef typename RigidBodyNodeSpec<3, true, noX_MB, noR_PF>::HType HType;
virtual const char* type() { return "translate"; }

RBNodeTranslate(const MassProperties& mProps_B,
                const Transform&      X_PF,
                const Transform&      X_BM,
                bool                  isReversed,
                UIndex&               nextUSlot,
                USquaredIndex&        nextUSqSlot,
                QIndex&               nextQSlot)
:   RigidBodyNodeSpec<3, true, noX_MB, noR_PF>(mProps_B,X_PF,X_BM,nextUSlot,nextUSqSlot,nextQSlot,
                         RigidBodyNode::QDotIsAlwaysTheSameAsU, RigidBodyNode::QuaternionIsNeverUsed, 
                         isReversed)
{
    this->updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
}

    // Implementations of virtual methods.

void setQToFitRotationImpl(const SBStateDigest& sbs, const Rotation& R_FM, 
                           Vector& q) const {
    // the only rotation this mobilizer can represent is identity
}
void setQToFitTranslationImpl(const SBStateDigest& sbs, const Vec3&  p_FM, 
                              Vector& q) const {
    // here's what this joint is really good at!
    this->toQ(q) = p_FM;
}

void setUToFitAngularVelocityImpl(const SBStateDigest& sbs, const Vector&, 
                                  const Vec3& w_FM, Vector& u) const {
    // The only angular velocity this can represent is zero.
}
void setUToFitLinearVelocityImpl
   (const SBStateDigest& sbs, const Vector&, const Vec3& v_FM, Vector& u) const
{
    // linear velocity is in a Cartesian joint's sweet spot
    this->toU(u) = v_FM;
}

// A translation joint doesn't need to cache any q calculations.
int calcQPoolSize(const SBModelVars&) const
{   return 0; }

// Nothing to precalculate.
void performQPrecalculations(const SBStateDigest& sbs,
                             const Real* q, int nq,
                             Real* qCache,  int nQCache,
                             Real* qErr,    int nQErr) const
{
    assert(q && nq==3 && nQCache==0 && nQErr==0);
}

// This is free.
void calcX_FM(const SBStateDigest& sbs,
              const Real* q,      int nq,
              const Real* qCache, int nQCache,
              Transform&  X_FM) const
{
    assert(q && nq==3 && nQCache==0);
    // Translation vector q is expressed in F (and M since they have same 
    // orientation). A translation (Cartesian) joint can't change orientation. 
    X_FM = Transform(Rotation(), Vec3::getAs(&q[0]));
}

// Generalized speeds together are the velocity of M's origin in the F frame,
// expressed in F. So individually they produce velocity along F's x,y,z
// axes respectively.
void calcAcrossJointVelocityJacobian(
    const SBStateDigest& sbs,
    HType&               H_FM) const
{
    H_FM(0) = SpatialVec( Vec3(0), Vec3(1,0,0) );
    H_FM(1) = SpatialVec( Vec3(0), Vec3(0,1,0) );
    H_FM(2) = SpatialVec( Vec3(0), Vec3(0,0,1) );
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

// Override the computation of reverse-H for this simple mobilizer.
void calcReverseMobilizerH_FM(
    const SBStateDigest& sbs,
    HType&               H_FM) const
{
    H_FM(0) = SpatialVec( Vec3(0), Vec3(-1, 0, 0) );
    H_FM(1) = SpatialVec( Vec3(0), Vec3( 0,-1, 0) );
    H_FM(2) = SpatialVec( Vec3(0), Vec3( 0, 0,-1) );
}

// Override the computation of reverse-HDot for this simple mobilizer.
void calcReverseMobilizerHDot_FM(
    const SBStateDigest& sbs,
    HType&               HDot_FM) const
{
    HDot_FM(0) = SpatialVec( Vec3(0), Vec3(0) );
    HDot_FM(1) = SpatialVec( Vec3(0), Vec3(0) );
    HDot_FM(2) = SpatialVec( Vec3(0), Vec3(0) );
}

};


#endif // SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_TRANSLATION_H_

