#ifndef SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_CYLINDER_H_
#define SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_CYLINDER_H_

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
 * Define the RigidBodyNode that implements a Cylinder mobilizer.
 */

#include "SimbodyMatterSubsystemRep.h"
#include "RigidBodyNode.h"
#include "RigidBodyNodeSpec.h"


    // CYLINDER //

// This is a "cylinder" joint, meaning one degree of rotational freedom
// about a particular axis, and one degree of translational freedom
// along the same axis. For molecules you can think of this as a combination
// of torsion and bond stretch. The axis used is the z axis of the parent's
// F frame, which is aligned forever with the z axis of the body's M frame.
// In addition, the origin points of M and F are separated only along the
// z axis; i.e., they have the same x & y coords in the F frame. The two
// generalized coordinates are the rotation and the translation, in that order.
class RBNodeCylinder : public RigidBodyNodeSpec<2> {
public:
    virtual const char* type() { return "cylinder"; }

    RBNodeCylinder(const MassProperties& mProps_B,
                   const Transform&      X_PF,
                   const Transform&      X_BM,
                   bool                  isReversed,
                   UIndex&               nextUSlot,
                   USquaredIndex&        nextUSqSlot,
                   QIndex&               nextQSlot)
      : RigidBodyNodeSpec<2>(mProps_B,X_PF,X_BM,nextUSlot,nextUSqSlot,nextQSlot,
                             QDotIsAlwaysTheSameAsU, QuaternionIsNeverUsed, isReversed)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }


    void setQToFitRotationImpl(const SBStateDigest& sbs, const Rotation& R_FM, Vector& q) const {
        // The only rotation our cylinder joint can handle is about z.
        // TODO: this code is bad -- see comments for Torsion joint above.
        const Vec3 angles123 = R_FM.convertRotationToBodyFixedXYZ();
        toQ(q)[0] = angles123[2];
    }

    void setQToFitTranslationImpl(const SBStateDigest& sbs, const Vec3& p_FM, Vector& q) const {
        // Because the M and F origins must lie along their shared z axis, there is no way to
        // create a translation by rotating around z. So the only translation we can represent
        // is that component which is along z.
        toQ(q)[1] = p_FM[2];
    }

    void setUToFitAngularVelocityImpl(const SBStateDigest& sbs, const Vector&, const Vec3& w_FM, Vector& u) const {
        // We can only represent an angular velocity along z with this joint.
        toU(u)[0] = w_FM[2];
    }

    void setUToFitLinearVelocityImpl
       (const SBStateDigest& sbs, const Vector&, const Vec3& v_FM, Vector& u) const
    {
        // Because the M and F origins must lie along their shared z axis, there is no way to
        // create a linear velocity by rotating around z. So the only linear velocity we can represent
        // is that component which is along z.
        toU(u)[1] = v_FM[2];
    }

    // This is required for all mobilizers.
    bool isUsingAngles(const SBStateDigest& sbs, MobilizerQIndex& startOfAngles, int& nAngles) const {
        // Cylinder joint has one angular coordinate, which comes first.
        startOfAngles = MobilizerQIndex(0); nAngles=1; 
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
        const Real& angle = fromQ(q)[0];
        toQ(sine)[0]    = std::sin(angle);
        toQ(cosine)[0]  = std::cos(angle);
        // no quaternions
    }

    // Calculate X_FM.
    void calcAcrossJointTransform(
        const SBStateDigest& sbs,
        const Vector&        q,
        Transform&           X_FM) const
    {
        const Vec2& coords  = fromQ(q);

        X_FM.updR().setRotationFromAngleAboutZ(coords[0]);
        X_FM.updP() = Vec3(0,0,coords[1]);
    }


    // The generalized speeds are (1) the angular velocity of M in the F frame,
    // about F's z axis, expressed in F, and (2) the velocity of M's origin
    // in F, along F's z axis. (The z axis is also constant in M for this joint.)
    void calcAcrossJointVelocityJacobian(
        const SBStateDigest& sbs,
        HType&               H_FM) const
    {
        H_FM(0) = SpatialVec( Vec3(0,0,1), Vec3(0)     );
        H_FM(1) = SpatialVec( Vec3(0),     Vec3(0,0,1) );
    }

    // Since the Jacobian above is constant in F, its time derivative in F is zero.
    void calcAcrossJointVelocityJacobianDot(
        const SBStateDigest& sbs,
        HType&               HDot_FM) const
    {
        HDot_FM(0) = SpatialVec( Vec3(0), Vec3(0) );
        HDot_FM(1) = SpatialVec( Vec3(0), Vec3(0) );
    }

    // Override the computation of reverse-H for this simple mobilizer.
    void calcReverseMobilizerH_FM(
        const SBStateDigest& sbs,
        HType&               H_FM) const
    {
        H_FM(0) = SpatialVec( Vec3(0,0,-1), Vec3(0)     );
        H_FM(1) = SpatialVec( Vec3(0),      Vec3(0,0,-1) );
    }

    // Override the computation of reverse-HDot for this simple mobilizer.
    void calcReverseMobilizerHDot_FM(
        const SBStateDigest& sbs,
        HType&               HDot_FM) const
    {
        HDot_FM(0) = SpatialVec( Vec3(0), Vec3(0) );
        HDot_FM(1) = SpatialVec( Vec3(0), Vec3(0) );
    }
};



#endif // SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_CYLINDER_H_

