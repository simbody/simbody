#ifndef SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_SCREW_H_
#define SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_SCREW_H_

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
 * Define the RigidBodyNode that implements a Screw mobilizer.
 */

#include "SimbodyMatterSubsystemRep.h"
#include "RigidBodyNode.h"
#include "RigidBodyNodeSpec.h"


    // SCREW //

// This is a one-dof "screw" joint, meaning one degree of rotational freedom
// about a particular axis, coupled to translation along that same axis.
// Here we use the common z axis of the F and M frames, which remains
// aligned forever. 
// For the generalized coordinate q, we use the rotation angle. For the
// generalized speed u we use the rotation rate, which is also the
// angular velocity of M in F (about the z axis). We compute the
// translational position as pitch*q, and the translation rate as pitch*u.
class RBNodeScrew : public RigidBodyNodeSpec<1> {
    Real pitch;
public:
    virtual const char* type() { return "screw"; }

    RBNodeScrew(const MassProperties& mProps_B,
                const Transform&      X_PF,
                const Transform&      X_BM,
                Real                  p,  // the pitch
                bool                  isReversed,
                UIndex&               nextUSlot,
                USquaredIndex&        nextUSqSlot,
                QIndex&               nextQSlot)
      : RigidBodyNodeSpec<1>(mProps_B,X_PF,X_BM,nextUSlot,nextUSqSlot,nextQSlot,
                             QDotIsAlwaysTheSameAsU, QuaternionIsNeverUsed, isReversed),
        pitch(p)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }

    void setQToFitRotationImpl(const SBStateDigest& sbs, const Rotation& R_FM, Vector& q) const {
        // The only rotation our screw joint can handle is about z.
        // TODO: should use 321 to deal with singular configuration (angle2==pi/2) better;
        // in that case 1 and 3 are aligned and the conversion routine allocates all the
        // rotation to whichever comes first.
        // TODO: isn't there a better way to come up with "the rotation around z that
        // best approximates a rotation R"?
        const Vec3 angles123 = R_FM.convertRotationToBodyFixedXYZ();
        to1Q(q) = angles123[2];
    }

    void setQToFitTranslationImpl(const SBStateDigest& sbs, const Vec3& p_FM, Vector& q) const {
        to1Q(q) = p_FM[2]/pitch;
    }

    void setUToFitAngularVelocityImpl(const SBStateDigest& sbs, const Vector&, const Vec3& w_FM, Vector& u) const {
        // We can only represent an angular velocity along z with this joint.
        to1U(u) = w_FM[2]; // project angular velocity onto z axis
    }

    void setUToFitLinearVelocityImpl
       (const SBStateDigest& sbs, const Vector&, const Vec3& v_FM, Vector& u) const
    {
        to1U(u) = v_FM[2]/pitch;
    }

    // This is required for all mobilizers.
    bool isUsingAngles(const SBStateDigest& sbs, MobilizerQIndex& startOfAngles, int& nAngles) const {
        // We're currently using an angle as the generalized coordinate for the screw joint
        // but could just as easily have used translation or some non-physical coordinate. It
        // might make sense to offer a Model stage option to set the coordinate meaning.
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
        const Real& angle = from1Q(q); // angular coordinate
        to1Q(sine)    = std::sin(angle);
        to1Q(cosine)  = std::cos(angle);
        // no quaternions
    }

    // Calculate X_FM.
    void calcAcrossJointTransform(
        const SBStateDigest& sbs,
        const Vector&        q,
        Transform&           X_FM) const
    {
        const Real& theta  = from1Q(q);    // angular coordinate

        X_FM.updR().setRotationFromAngleAboutZ(theta);
        X_FM.updP() = Vec3(0,0,theta*pitch);
    }


    // The generalized speed is the angular velocity of M in the F frame,
    // about F's z axis, expressed in F. (This axis is also constant in M.)
    void calcAcrossJointVelocityJacobian(
        const SBStateDigest& sbs,
        HType&               H_FM) const
    {
        H_FM(0) = SpatialVec( Vec3(0,0,1), Vec3(0,0,pitch) );
    }

    // Since the Jacobian above is constant in F, its time derivative in F is zero.
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
        H_FM(0) = SpatialVec( Vec3(0,0,-1), Vec3(0,0,-pitch) );
    }

    // Override the computation of reverse-HDot for this simple mobilizer.
    void calcReverseMobilizerHDot_FM(
        const SBStateDigest& sbs,
        HType&               HDot_FM) const
    {
        HDot_FM(0) = SpatialVec( Vec3(0), Vec3(0) );
    }
};



#endif // SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_SCREW_H_

