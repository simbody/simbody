#ifndef SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_PIN_H_
#define SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_PIN_H_

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
 * Define the RigidBodyNode that implements a Pin mobilizer, also known as a
 * Torsion or Revolute joint.
 */

#include "SimbodyMatterSubsystemRep.h"
#include "RigidBodyNode.h"
#include "RigidBodyNodeSpec.h"


    // PIN (TORSION) //

// This is a "pin" or "torsion" or "revolute" joint, meaning one degree of 
// rotational freedom about a particular axis, the z axis of the parent's F 
// frame, which is aligned forever with the z axis of the body's M frame. In 
// addition, the origin points Mo of M and Fo of F are identical forever.
template<bool noX_MB, bool noR_PF>
class RBNodeTorsion : public RigidBodyNodeSpec<1, false, noX_MB, noR_PF> {
public:
virtual const char* type() { return "torsion"; }
typedef typename RigidBodyNodeSpec<1, false, noX_MB, noR_PF>::HType HType;

RBNodeTorsion(const MassProperties&   mProps_B,
                const Transform&      X_PF,
                const Transform&      X_BM,
                bool                  isReversed,
                UIndex&               nextUSlot,
                USquaredIndex&        nextUSqSlot,
                QIndex&               nextQSlot)
:   RigidBodyNodeSpec<1, false, noX_MB, noR_PF>(mProps_B,X_PF,X_BM,nextUSlot,nextUSqSlot,nextQSlot,
                         RigidBodyNode::QDotIsAlwaysTheSameAsU, RigidBodyNode::QuaternionIsNeverUsed, 
                         isReversed)
{
    this->updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
}

void setQToFitRotationImpl(const SBStateDigest& sbs, const Rotation& R_FM, 
                           Vector& q) const {
    // The only rotation our pin joint can handle is about z.
    // TODO: should use 321 to deal with singular configuration (angle2==pi/2) 
    // better; in that case 1 and 3 are aligned and the conversion routine 
    // allocates all the rotation to whichever comes first.
    // TODO: isn't there a better way to come up with "the rotation around z 
    // that best approximates a rotation R"?
    const Vec3 angles123 = R_FM.convertRotationToBodyFixedXYZ();
    this->to1Q(q) = angles123[2];
}

void setQToFitTranslationImpl(const SBStateDigest& sbs, const Vec3& p_FM, 
                              Vector& q) const {
    // M and F frame origins are always coincident for this mobilizer so there 
    // is no way to create a translation by rotating. So the only translation 
    // we can represent is 0.
}

void setUToFitAngularVelocityImpl(const SBStateDigest& sbs, const Vector&, 
                                  const Vec3& w_FM, Vector& u) const {
    // We can only represent an angular velocity along z with this joint.
    this->to1U(u) = w_FM[2]; // project angular velocity onto z axis
}

void setUToFitLinearVelocityImpl
    (const SBStateDigest& sbs, const Vector&, const Vec3& v_FM, Vector& u) const
{
    // M and F frame origins are always coincident for this mobilizer so there 
    // is no way to create a linear velocity by rotating. So the only linear 
    // velocity we can represent is 0.
}

enum {PoolSize=2}; // number of Reals
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
    assert(q && nq==1 && qCache && nQCache==PoolSize && nQErr==0);
    qCache[CosQ] = std::cos(q[0]);
    qCache[SinQ] = std::sin(q[0]);
}

// This is nearly free since we already calculated sin/cos.
void calcX_FM(const SBStateDigest& sbs,
                const Real* q,      int nq,
                const Real* qCache, int nQCache,
                Transform&  X_FM) const
{
    assert(q && nq==1 && qCache && nQCache==PoolSize);
    X_FM.updR().setRotationFromAngleAboutZ(qCache[CosQ], qCache[SinQ]);
    X_FM.updP() = 0; // can't translate
}

// The generalized speed is the angular velocity of M in the F frame,
// about F's z axis, expressed in F. (This axis is also constant in M.)
void calcAcrossJointVelocityJacobian(
    const SBStateDigest& sbs,
    HType&               H_FM) const
{
    H_FM(0) = SpatialVec( Vec3(0,0,1), Vec3(0) );
}


// Since H_FM above is constant in F, its time derivative in F is zero.
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
    H_FM(0) = SpatialVec( Vec3(0,0,-1), Vec3(0) );
}

// Override the computation of reverse-HDot for this simple mobilizer.
void calcReverseMobilizerHDot_FM(
    const SBStateDigest& sbs,
    HType&               HDot_FM) const
{
    // doesn't get better than this!
    HDot_FM(0) = SpatialVec( Vec3(0), Vec3(0) ); 
}

};


#endif // SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_PIN_H_

