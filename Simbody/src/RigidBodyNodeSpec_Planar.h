#ifndef SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_PLANAR_H_
#define SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_PLANAR_H_

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
class RBNodePlanar : public RigidBodyNodeSpec<3> {
public:
    virtual const char* type() { return "planar"; }

    RBNodePlanar(const MassProperties& mProps_B,
                 const Transform&      X_PF,
                 const Transform&      X_BM,
                 bool                  isReversed,
                 UIndex&               nextUSlot,
                 USquaredIndex&        nextUSqSlot,
                 QIndex&               nextQSlot)
      : RigidBodyNodeSpec<3>(mProps_B,X_PF,X_BM,nextUSlot,nextUSqSlot,nextQSlot,
                             QDotIsAlwaysTheSameAsU, QuaternionIsNeverUsed, isReversed)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }

        // Implementations of virtual methods.

    void setQToFitRotationImpl(const SBStateDigest& sbs, const Rotation& R_FM, Vector& q) const {
        // The only rotation our planar joint can handle is about z.
        // TODO: should use 321 to deal with singular configuration (angle2==pi/2) better;
        // in that case 1 and 3 are aligned and the conversion routine allocates all the
        // rotation to whichever comes first.
        // TODO: isn't there a better way to come up with "the rotation around z that
        // best approximates a rotation R"?
        const Vec3 angles123 = R_FM.convertRotationToBodyFixedXYZ();
        toQ(q)[0] = angles123[2];
    }
    void setQToFitTranslationImpl(const SBStateDigest& sbs, const Vec3&  p_FM, Vector& q) const {
        // Ignore translation in the z direction.
        toQ(q)[1] = p_FM[0]; // x
        toQ(q)[2] = p_FM[1]; // y
    }

    void setUToFitAngularVelocityImpl(const SBStateDigest& sbs, const Vector&, const Vec3& w_FM, Vector& u) const {
        // We can represent the z angular velocity exactly, but nothing else.
        toU(u)[0] = w_FM[2];
    }
    void setUToFitLinearVelocityImpl
       (const SBStateDigest& sbs, const Vector&, const Vec3& v_FM, Vector& u) const
    {
        // Ignore translational velocity in the z direction.
        toU(u)[1] = v_FM[0]; // x
        toU(u)[2] = v_FM[1]; // y
    }

    // This is required for all mobilizers.
    bool isUsingAngles(const SBStateDigest& sbs, MobilizerQIndex& startOfAngles, int& nAngles) const {
        // Planar joint has one angular coordinate, which comes first.
        startOfAngles = MobilizerQIndex(0); nAngles=1; 
        return true;
    }

    // This is required but does nothing here since there are no rotations for this joint.
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
        const Real& angle = fromQ(q)[0]; // angular coordinate
        to1Q(sine)    = std::sin(angle);
        to1Q(cosine)  = std::cos(angle);
        // no quaternions
    }

    // Calculate X_FM.
    void calcAcrossJointTransform(
        const SBStateDigest& mv,
        const Vector&        q,
        Transform&           X_FM) const
    {
        // Rotational q is about common z axis, translational q's along Fx and Fy.
        X_FM = Transform(Rotation( fromQ(q)[0], ZAxis ), 
                          Vec3(fromQ(q)[1], fromQ(q)[2], 0));
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

