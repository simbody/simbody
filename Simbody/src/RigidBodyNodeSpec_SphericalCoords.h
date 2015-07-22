#ifndef SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_SPHERICALCOORDS_H_
#define SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_SPHERICALCOORDS_H_

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
 * Define the RigidBodyNode that implements a SphericalCoods mobilizer.
 */

#include "SimbodyMatterSubsystemRep.h"
#include "RigidBodyNode.h"
#include "RigidBodyNodeSpec.h"

/* This is a mobilizer whose generalized coordinates are directly interpretable
as a spherical coordinate system in this order: azimuth (longitude), zenith
(latitude-90), radius. In the reference configuration (all coordinates zero)
the parent's F frame and child's M frame are aligned, with the z axes pointing
to the north pole. The equator is in the common x-y plane. Then the azimuth is
a rotation of M about the common z axis. Next, the zenith angle is a rotation
about My, the *new* child body y axis. This angle measures the deviation of Mz
from Fz, that is, 0 when Mz is pointing to F's north pole; Pi (180 degrees)
when Mz points to the south pole. The final coordinate is a translation along
a selected axis (Mx or Mz; default is Mz). Whenever the translation is zero,
the F and M origins are coincident.

This mobilizer can be used to connect an inertialess particle to its parent in
an unrestricted way that corresponds directly to (torsion, bend, stretch)
coordinates, with two caveats: (1) the joint is singular whenever the stretch
is zero, and (2) the joint is singular when the bend angle is 0 or Pi, because
then the torsion rotates the particle about its center. The joint is never
singular if the child body has a full inertia matrix.

Because there are an assortment of useful conventions which are minor
variations on the spherical coordinate theme, we allow some flexibility in
defining the coordinates here. Although the rotation axes are fixed in a 3-2
body-fixed sequence, the third (sliding) axis can be selected as Mx or Mz.
Any of the three coordinates can be defined to use a reversed sign convention,
and you can "bake in" fixed offsets for the angles. So we define the mobilizer
in terms of aximuth angle a, zenith angle z, and stretch length d and define
these in terms of the generalized coordinates like this:
<pre>
     a = q0*s0 + a0  (about Mz==Fz)
     z = q1*s1 + z0  (about My)
     d = q2*s2       (along v from {Mx,Mz})
</pre>
In all cases we have qi_dot = ui.

With this you can define a "geographical" coordinate system where Mx is the
Greenwich line, a is latitude and z longitude (with north positive):
<pre>
     v  = Mx
     s0 =  1, a0 = 0
     s1 = -1, z0 = 0
     s2 =  1
</pre>
If you want the translation direction to be in Mz (the default) but would like
q1=0 to mean the equatorial position rather than the (possibly singular) north
pole which should be -90, define
<pre>
     v  = Mz
     s0 = 1, a0 = 0
     s1 = 1, z0 = Pi/2
     s2 = 1
</pre>
One common convention for atomic (torsion,bend,stretch) uses the default
spherical coordinate system but the final stretch is along the -z direction.
For that, take all defaults but set s2=-1. */

template<bool noX_MB, bool noR_PF>
class RBNodeSphericalCoords : public RigidBodyNodeSpec<3, false, noX_MB, noR_PF> {
public:
typedef typename RigidBodyNodeSpec<3, false, noX_MB, noR_PF>::HType HType;
virtual const char* type() { return "spherical coords"; }

RBNodeSphericalCoords(const MassProperties&   mProps_B,
                        const Transform&      X_PF,
                        const Transform&      X_BM,
                        Real                  azimuthOffset, // radians
                        bool                  azimuthNegated,
                        Real                  zenithOffset,
                        bool                  zenithNegated,
                        CoordinateAxis        translationAxis, // only x or z
                        bool                  translationNegated,
                        bool                  isReversed,
                        UIndex&               nextUSlot,
                        USquaredIndex&        nextUSqSlot,
                        QIndex&               nextQSlot)
:   RigidBodyNodeSpec<3, false, noX_MB, noR_PF>(mProps_B,X_PF,X_BM,nextUSlot,nextUSqSlot,nextQSlot,
                         RigidBodyNode::QDotIsAlwaysTheSameAsU, RigidBodyNode::QuaternionIsNeverUsed,
                         isReversed),
    az0(azimuthOffset), ze0(zenithOffset), axisT(translationAxis),
    signAz(Real(azimuthNegated ? -1 : 1)), signZe(Real(zenithNegated ? -1 : 1)),
    signT(Real(translationNegated ? -1 : 1))
{
    SimTK_ASSERT_ALWAYS( translationAxis != YAxis,
        "RBNodeSphericalCoords: translation axis must be x or z; not y.");
    this->updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
}

void setQToFitRotationImpl(const SBStateDigest& sbs, const Rotation& R_FM,
                           Vector& q) const {
    // The only rotations this joint can handle are about Mz and My.
    // TODO: isn't there a better way to come up with "the rotation around z&y
    // that best approximates a rotation R"? Here we're just hoping that the
    // supplied rotation matrix can be decomposed into (z,y) rotations.
    const Vec2 angles = R_FM.convertTwoAxesRotationToTwoAngles
                                        (BodyRotationSequence, ZAxis, YAxis);
    this->toQ(q)[0] = signAz * (angles[0] - az0);
    this->toQ(q)[1] = signZe * (angles[1] - ze0);
}

void setQToFitTranslationImpl(const SBStateDigest& sbs, const Vec3& p_FM,
                              Vector& q) const {
    // We're not going to attempt any rotations to fit the translation. We'll
    // find out which way the translation axis is currently pointing and move
    // along it.
    const Rotation R_FM = calcR_FM(q);
    // note: rows of R_FM are columns of R_MF.
    this->toQ(q)[2] = signT * dot(p_FM, R_FM[axisT]);
}

// We can only express angular velocity that can be produced with our
// generalized speeds which are Fz(=Mz) and My rotation rates. So we'll take
// the supplied angular velocity expressed in F, project it on Fz and use that
// as the first generalized speed. Then take whatever angular velocity is
// unaccounted for, express it in M, and project onto My and use that as the
// second generalized speed.
void setUToFitAngularVelocityImpl(const SBStateDigest& sbs, const Vector& q,
                                  const Vec3& w_FM, Vector& u) const {
    const Rotation R_FM = calcR_FM(q);
    this->toU(u)[0] = w_FM[2];
    // Calculate the remaining angular velocity w-wz and re-express in M.
    const Vec3 wxy_FM_M = ~R_FM*Vec3(w_FM[0],w_FM[1],0);
    this->toU(u)[1] = wxy_FM_M[1]; // can only deal with the y angular velocity now
}

// Although we could try to use angular velocity to affect linear velocity,
// we're going to restrain ourselves here and only change the translational
// velocity along the final axis.
void setUToFitLinearVelocityImpl
   (const SBStateDigest& sbs, const Vector& q, const Vec3& v_FM, Vector& u) const
{
    const Rotation R_FM = calcR_FM(q);
    this->toU(u)[2] = dot(v_FM, R_FM[axisT]); // i.e., v_FM*Mx_F or v_FM*Mz_F.
}

enum {PoolSize=4}; // number of Reals
enum {CosQ=0, SinQ=2};
// We want space for cos(q01) and sin(q01).
int calcQPoolSize(const SBModelVars&) const
{   return PoolSize; }

// Precalculate sines and cosines. Note: these are sines and cosines of q,
// not az and el.
// az = s*q + az0, so
// sin(az) = sin(s*q + az0)
//    = sin(s*q)cos(az0) + cos(s*q)sin(az0)
//    = s*sin(q)cos(az0) + cos(q)sin(az0)
// cos(az) = cos(s*q + az0)
//    = cos(s*q)cos(az0) - sin(s*q)sin(az0)
//    = cos(q)cos(az0) - s*sin(q)sin(az0)
// and similarly for the zenith angle. We can precalculate sine and cosine of
// the offsets so we need only 4 flops to get sin or cos of the actual angles
// given sin and cos of the generalized coordinates.
// Cost here is about 100 flops.
void performQPrecalculations(const SBStateDigest& sbs,
                             const Real* q, int nq,
                             Real* qCache,  int nQCache,
                             Real* qErr,    int nQErr) const
{
    assert(q && nq==3 && qCache && nQCache==PoolSize && nQErr==0);
    Vec2::updAs(&qCache[CosQ]) = Vec2(std::cos(q[0]),std::cos(q[1]));
    Vec2::updAs(&qCache[SinQ]) = Vec2(std::sin(q[0]),std::sin(q[1]));
}

// TODO: should be using sin/cos cache but isn't.
void calcX_FM(const SBStateDigest& sbs,
              const Real* q,      int nq,
              const Real* qCache, int nQCache,
              Transform&  X_FM) const
{
    assert(q && nq==3 && qCache && nQCache==PoolSize);
    // Calculate azimuth and zenith angles from the first two qs.
    const Vec2 azZe = calcAzZe(q[0], q[1]);

    // Calculate R_FM from q's. This is always a body fixed Z-Y (3-2)
    // sequence although the angles may be negated and shifted from the
    // generalized coordinates.
    // TODO: use qCache for speed.
    X_FM.updR() = Rotation( BodyRotationSequence,
                            azZe[0], ZAxis, azZe[1], YAxis );

    const Real t = signT * q[2];
    X_FM.updP() = t * X_FM.R()(axisT); // i.e., t*Mx or t*Mz, expressed in F
}

// The generalized speeds for this joint are the time derivatives of the
// body-fixed 3-2 rotation sequence defining the orientation, and the
// translation velocity. That is, the first speed is just a rotation rate about
// Fx. The second is a rotation rate about the current My, so we have to
// transform it into F to make H_FM uniformly expressed in F. The third axis is
// either Mx or My and is likewise transformed.
// 18 flops normally; 36 if reversed.
void calcAcrossJointVelocityJacobian(
    const SBStateDigest& sbs,
    HType&               H_FM) const
{
    // "upd" because we're realizing positions now
    const SBTreePositionCache& pc = sbs.updTreePositionCache();
    const Transform  X_F0M0 = this->findX_F0M0(pc);   // 18 flops if reversed

    // Dropping the 0's here.
    const Rotation& R_FM = X_F0M0.R();
    const Vec3&     p_FM = X_F0M0.p();

    const Vec3 sFz = Vec3(0,0,signAz);  // signed Fz for 1st rotation
    const Vec3 sMy = signZe*R_FM.y();   // signed My for 2nd rotation (3 flops)
    const Vec3 sMt = signT*R_FM(axisT); // signed Mz or Mx for trans. (3 flops)

    // We need sFz % p_FM but sFz is the z or -z axis so this is a very
    // cheap calculation.
    const Vec3 sFzXp_FM = Vec3( -sFz[2]*p_FM[1], sFz[2]*p_FM[0], 0); // 3 flops

    H_FM(0) = SpatialVec(  sFz    ,  sFzXp_FM );
    H_FM(1) = SpatialVec(  sMy    , sMy % p_FM );                    // 9 flops
    H_FM(2) = SpatialVec( Vec3(0) ,    sMt );
}

// Since the last two rows of H_FM above are not constant in F, their time
// derivatives are non zero. Here we use the fact that for a vector r_B_A
// fixed in a moving frame B but expressed in another frame A, its time
// derivative in A is the angular velocity of B in A crossed with the vector,
// i.e., d_A/dt r_B_A = w_AB % r_B_A. 48 flops normally; 111 if reversed.
void calcAcrossJointVelocityJacobianDot(
    const SBStateDigest& sbs,
    HType&               HDot_FM) const
{
    const SBTreePositionCache& pc = sbs.getTreePositionCache();
    // "upd" because we're realizing velocities now
    const SBTreeVelocityCache& vc = sbs.updTreeVelocityCache();

    const Transform  X_F0M0 = this->findX_F0M0(pc);       // 18 flops if reversed
    const SpatialVec V_F0M0 = this->findV_F0M0(pc,vc);    // 45 flops if reversed

    // Dropping the 0's here.
    const Rotation& R_FM = X_F0M0.R();
    const Vec3&     p_FM = X_F0M0.p();

    const Vec3  sFz = Vec3(0,0,signAz);  // signed Fz for 1st rotation
    const Vec3  sMy = signZe*R_FM.y();   // signed My for 2nd rotation (3 flops)
    const Vec3  sMt = signT*R_FM(axisT); // signed Mz or Mx for trans. (3 flops)

    const Vec3  w_FM = V_F0M0[0]; // angular velocity of M in F
    const Vec3& v_FM = V_F0M0[1]; // linear velocity of OM in F

    const Vec3  dsMy = w_FM % sMy;  // time derivative of sMy         (9 flops)
    const Vec3  dsMt = w_FM % sMt;  // time derivative of sMz or sMx  (9 flops)

    // We need d/dt (sFz % p_FM) = (sFz % v_FM). But since sFz is the
    // z or -z axis  this is a very cheap calculation.
    const Vec3  sFzXv_FM = Vec3( -sFz[2]*v_FM[1], sFz[2]*v_FM[0], 0); // 3 flops

    HDot_FM(0) = SpatialVec( Vec3(0) ,         sFzXv_FM );
    HDot_FM(1) = SpatialVec(  dsMy   , dsMy % p_FM + sMy % v_FM );  // 21 flops
    HDot_FM(2) = SpatialVec( Vec3(0) ,           dsMt );
}

private:
const Real              az0, ze0;               // angle offsets
const CoordinateAxis    axisT;                  // translation axis (X or Z)
const Real              signAz, signZe, signT;  // 1 or -1

// Calculate azimuth and zenith angles from the first two qs.
Vec2 calcAzZe(Real q0, Real q1) const {
    const Real az = signAz * q0 + az0;
    const Real ze = signZe * q1 + ze0;
    return Vec2(az,ze);
}

// Calculate R_FM from q's. This is always a body fixed 3-2 sequence
// although the angles may be negated and shifted from the generalized
// coordinates.
Rotation calcR_FM(const Vector& q) const {
    const Vec2 azZe = calcAzZe(this->fromQ(q)[0], this->fromQ(q)[1]);
    return Rotation( BodyRotationSequence, azZe[0], ZAxis, azZe[1], YAxis );
}

};


#endif // SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_SPHERICALCOORDS_H_

