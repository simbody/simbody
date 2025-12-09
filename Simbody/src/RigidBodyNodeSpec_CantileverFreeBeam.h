#ifndef SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_CANTILEVER_FREE_BEAM_H_
#define SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_CANTILEVER_FREE_BEAM_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-25 Stanford University and the Authors.      *
 * Authors: Nicholas Bianco                                                   *
 * Contributors:                                                              *
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
 * Define the RigidBodyNode that implements a CantileverFreeBeam mobilizer.
 */

#include "SimbodyMatterSubsystemRep.h"
#include "RigidBodyNode.h"
#include "RigidBodyNodeSpec.h"

    // CANTILEVER FREE BEAM //

// CANTILEVER FREE BEAM mobilizer. This provides three degrees of rotational
// freedom of the body's M frame in the parent's F frame, along with coordinated
// translation that defines the position of the M frame origin based on the
// deflection of a cantilever beam subject to a point load at the end of the
// beam.
//
// Unlike most joints, the reference configuration (i.e., X_FM when q=0) is
// not the identity transform. Instead, although the frames are aligned the
// M frame origin is offset from F along their shared +z axis, so that it lies
// at the point (0,0,L) where L is the length of the undeflected beam.
//
// The generalized coordinates q are:
//   * 3 X-Y-Z (1-2-3) body-fixed Euler angles (that is, fixed in M)
//     Angles 1 and 2 induce deflections in the -Fy and Fx directions,
//     respectively, according to the equations for a endpoint-loaded
//     cantilever beam derived from Euler-Bernoulli beam theory. Angle 3 is the
//     spin of the beam about the Mz axis, which is always tangent to the beam
//     at the beam's endpoint.
//
// and generalized speeds u are:
//   * u = qdot, that is, the X-Y-Z body fixed Euler angle derivatives
//
// NOTE: This mobilizer has a singularity when the middle angle (q[1]) is near
// +/-90 degrees.

template<bool noX_MB, bool noR_PF>
class RBNodeCantileverFreeBeam :
        public RigidBodyNodeSpec<3, false, noX_MB, noR_PF> {
    Real length;  // length of the beam
    Real deflectionCoefficient;
    Real displacementCoefficient;
    Real displacementSpeedCoefficient;
public:

typedef typename RigidBodyNodeSpec<3, false, noX_MB, noR_PF>::HType HType;
virtual const char* type() { return "cantilever free beam"; }

RBNodeCantileverFreeBeam(const MassProperties& mProps_B,
                         const Transform&      X_PF,
                         const Transform&      X_BM,
                         const Real&           length,
                         bool                  isReversed,
                         UIndex&               nextUSlot,
                         USquaredIndex&        nextUSqSlot,
                         QIndex&               nextQSlot)
  : RigidBodyNodeSpec<3, false, noX_MB, noR_PF>(
        mProps_B, X_PF, X_BM,
        nextUSlot, nextUSqSlot, nextQSlot,
        RigidBodyNode::QDotIsAlwaysTheSameAsU,
        RigidBodyNode::QuaternionIsNeverUsed,
        isReversed),
    length(length)
{
    // Multiplying this term by the beam deflection angle gives the beam
    // deflection, which is the absolute value of the beam's end point position
    // in the Fx and Fy directions. This coefficient may also be used to
    // calculate the beam deflection speed.
    deflectionCoefficient = (2.0 / 3.0) * length;

    // Multiplying this term by the beam deflection angle squared gives the beam
    // displacement, which can be subtracted from the beam length to give the
    // Fz-position of the beam's endpoint.
    displacementCoefficient = (4.0 / 15.0) * length;

    this->updateSlots(nextUSlot, nextUSqSlot, nextQSlot);
}

    // Implementations of virtual methods.

void setQToFitRotationImpl(const SBStateDigest& sbs, const Rotation& R_FM,
    Vector& q) const override
{
    this->toQ(q) = R_FM.convertRotationToBodyFixedXYZ();
}

// We can't represent arbitrary translations with a joint that has only
// rotational coordinates, but we can find a set of rotations to obtain a
// translation in the direction of the requested translation. Based on the
// implementation from RigidBodyNodeSpec_Ellipsoid.h.
void setQToFitTranslationImpl(const SBStateDigest& sbs, const Vec3& p_FM,
                              Vector& q) const override {
   if (p_FM.norm() < Eps) return;

   const UnitVec3 e(p_FM); // direction from F origin towards desired M origin
   const Real latitude  = std::atan2(-e[1],e[2]); // project onto F's yz plane
   const Real longitude = std::atan2( e[0],e[2]); // project onto F's xz plane

   // Calculate the current value of the spin coordinate (3rd Euler angle).
   Real spin = this->fromQ(q)[2];

   // Calculate the desired rotation, which is a space-fixed 1-2 sequence for
   // latitude and longitude, followed by a body fixed rotation for spin.
   const Rotation R_FM =
       Rotation(SpaceRotationSequence, latitude, XAxis, longitude, YAxis) *
       Rotation(spin, ZAxis);

   this->toQ(q) = R_FM.convertRotationToBodyFixedXYZ();
}

// Given the angular velocity of M in F, expressed in F, compute the Euler
// angle derivatives qdot that would produce that angular velocity, and
// return u=qdot.
void setUToFitAngularVelocityImpl(const SBStateDigest& sbs, const Vector& q,
    const Vec3& w_FM, Vector& u) const override
{
    const Vec2 cosxy(std::cos(q[0]), std::cos(q[1]));
    const Vec2 sinxy(std::sin(q[0]), std::sin(q[1]));
    const Real oocosy = 1 / cosxy[1];
    const Vec3 qdot =
        Rotation::convertAngVelInParentToBodyXYZDot(cosxy,sinxy,oocosy,w_FM);
    this->toU(u) = qdot;
}

// The linear velocity of the beam's endpoint does not uniquely determine qdot,
// but we know that is a linear function of qdot with known q's. Therefore, we
// can estimate qdot via least squares.
void setUToFitLinearVelocityImpl(const SBStateDigest& sbs, const Vector& q,
    const Vec3& v_FM, Vector& u) const override
{
    Real q0 = this->fromQ(q)[0];
    Real q1 = this->fromQ(q)[1];

    Matrix m(3, 2, 0.0);

    // The y-component of qdot induces a positive Fx speed.
    m(0, 1) = deflectionCoefficient;

    // The x-component of qdot induces a negative Fy speed.
    m(1, 0) = -deflectionCoefficient;

    // The x- and y-components of qdot induce a negative Fz speed.
    m(2, 0) = -displacementSpeedCoefficient * q0;
    m(2, 1) = -displacementSpeedCoefficient * q1;

    // Solve for qdot0 and qdot1 via least squares.
    FactorQTZ qtz(m);
    Vector b(3, 0.0);
    b[0] = v_FM[0];
    b[1] = v_FM[1];
    b[2] = v_FM[2];
    Vector x(2, 0.0);
    qtz.solve(b, x);

    // Return the the first two elements of qdot we solved for above and leave
    // the third element unchanged.
    Vec3 qdot(x[0], x[1], this->fromU(u)[2]);
    this->toU(u) = qdot;
}

// We want to cache cos and sin for each angle, and also 1/cos of the middle
// angle will be handy to have around.
enum {PoolSize=7};
// cos x,y,z sin x,y,z 1/cos(y)
enum {CosQ=0, SinQ=3, OOCosQy=6};
int calcQPoolSize(const SBModelVars& mv) const override
{   return PoolSize; }

// This is expensive since we have three sin/cos computations and a divide
// to do, approx. 150 flops. But we hope to re-use these calculations several
// times before we're done with a complete realization.
void performQPrecalculations(const SBStateDigest& sbs,
                             const Real* q,      int nq,
                             Real*       qCache, int nQCache,
                             Real*       qErr,   int nQErr) const override
{
    assert(q && nq==3 && qCache && nQCache==PoolSize && nQErr==0);
    const Real cy = std::cos(q[1]);
    Vec3::updAs(&qCache[CosQ]) =
        Vec3(std::cos(q[0]), cy, std::cos(q[2]));
    Vec3::updAs(&qCache[SinQ]) =
        Vec3(std::sin(q[0]), std::sin(q[1]), std::sin(q[2]));
    qCache[OOCosQy] = 1/cy; // trouble at 90 or 270 (-90) degrees
}

void calcX_FM(const SBStateDigest& sbs,
              const Real* q,      int nq,
              const Real* qCache, int nQCache,
              Transform&  X_F0M0) const override
{
    assert(q && nq==3 && qCache && nQCache==PoolSize);

    X_F0M0.updR().setRotationToBodyFixedXYZ // 18 flops
        (Vec3::getAs(&qCache[CosQ]), Vec3::getAs(&qCache[SinQ]));

    const Real& q0 = Vec3::getAs(q)[0];
    const Real& q1 = Vec3::getAs(q)[1];
    X_F0M0.updP() = Vec3(
        q1 * deflectionCoefficient,
        -q0 * deflectionCoefficient,
        length - displacementCoefficient * (q0*q0 + q1*q1)
    );
}

// Generalized speeds are the Euler angle derivatives.
// Using the precalculations this requires only 5 flops, although there is
// a fair bit of poking around in memory required.
void calcAcrossJointVelocityJacobian(const SBStateDigest& sbs,
                                     HType& H_FM) const override
{
    const SBModelCache&        mc = sbs.getModelCache();
    // Use "upd" here because we're realizing positions now.
    const SBTreePositionCache& pc = sbs.updTreePositionCache();
    const Real*                pool = this->getQPool(mc, pc);

    const Real c0 = pool[CosQ], c1 = pool[CosQ+1];
    const Real s0 = pool[SinQ], s1 = pool[SinQ+1];
    const Vec3& q = this->fromQ(sbs.getQ());

    // Fill in columns of H_FM. See Rotation::calcNInvForBodyXYZInParentFrame().
    // Include the contributions of the Euler angle derivatives to the linear
    // velocity of the beam's endpoint.
    H_FM(0) = SpatialVec(Vec3(1,     0,    0),    Vec3(0,     -deflectionCoefficient, -2.0*displacementCoefficient*q[0]) );
    H_FM(1) = SpatialVec(Vec3(0,    c0,    s0),   Vec3(deflectionCoefficient, 0,      -2.0*displacementCoefficient*q[1]) );
    H_FM(2) = SpatialVec(Vec3(s1, -s0*c1, c0*c1), Vec3(0) );
}

// Differentiate H_FM to get HDot_FM. Note that this depends on qdot:
//    d/dt cos(q0) = -sin(q0)*qdot0, etc.
void calcAcrossJointVelocityJacobianDot(
    const SBStateDigest& sbs,
    HType&               HDot_FM) const override
{
    const SBModelCache&        mc = sbs.getModelCache();
    const SBTreePositionCache& pc = sbs.getTreePositionCache();

    const Real* pool = this->getQPool(mc, pc);
    const Real c0 = pool[CosQ], c1 = pool[CosQ+1];
    const Real s0 = pool[SinQ], s1 = pool[SinQ+1];

    // Use "upd" here because we're realizing velocities now.
    const Vec3& qdot = this->fromQ(sbs.updQDot());
    const Real qd0 = qdot[0], qd1 = qdot[1];

    const Real dc0 = -s0*qd0, dc1 = -s1*qd1; // derivatives of c0,c1,s0,s1
    const Real ds0 =  c0*qd0, ds1 =  c1*qd1;

    // Compare with H_FM above.
    HDot_FM(0) = SpatialVec(Vec3(0,         0,              0),       Vec3(0,  0,  -2.0*displacementCoefficient*qd0));
    HDot_FM(1) = SpatialVec(Vec3(0,        dc0,            ds0),      Vec3(0,  0,  -2.0*displacementCoefficient*qd1));
    HDot_FM(2) = SpatialVec(Vec3(ds1, -ds0*c1-s0*dc1, dc0*c1+c0*dc1), Vec3(0));
}

// Can use default for calcQDot, multiplyByN, etc., since qdot==u for
// CantileverFreeBeam mobilizer.

};

#endif // SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_CANTILEVER_FREE_BEAM_H_
