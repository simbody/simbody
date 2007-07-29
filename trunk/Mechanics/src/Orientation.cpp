/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmatrix(tm)                       *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-7 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
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
 * Implementations of non-inline methods of classes dealing
 * with orientation.
 */

#include "SimTKcommon/basics.h"
#include "SimTKcommon/internal/Scalar.h"
#include "SimTKcommon/internal/SmallMatrix.h"
#include "SimTKcommon/internal/Orientation.h"

#include <iostream>

namespace SimTK {

// Calculate a rotation matrix R_GF which defines the F
// coordinate frame by rotating the G frame z axis into alignment 
// with the passed-in zF vector (that is, zF is the z axis
// of R_GF expressed in G. The result is not unique.
// TODO: (sherm) I think this can be done with one sqrt and no trig functions.
// Just create a vector pointing in a substantially different
// direction than zF (e.g., move its largest coord), then cross it with zF to get the first 
// perpendicular. Normalize that, cross again and you have a frame.
// Must be careful about signs to get a right-handed set.
// TODO: (sherm) Uh, shouldn't the 3rd column of the result just be zF?
Rotation::Rotation(const UnitVec<1>& zF) {
    // Use the individual measure numbers (i.e., projections of zF
    // onto the G coordinate axes) to calculate spherical coordinates,
    // considering an equatorial x-y plane with z North.
    const Real ct=zF[2];    // cos(theta)
    const Real theta = std::acos( ct );             // zenith (90-elevation)
    const Real psi   = std::atan2( zF[0] , zF[1] ); // 90-azimuth
    const Real st=std::sin(theta), cp=std::cos(psi), sp=std::sin(psi);

    // This is a space fixed 1-2-3 sequence with angles
    // a1=-theta, a2=0, a3=-psi. That is, to get from G to F
    // first rotate by -theta around the G frame x axis, 
    // then rotate by -psi around the G frame z axis.

    BaseMat::operator=(BaseMat( cp, ct*sp, sp*st,
                               -sp, ct*cp, cp*st,
                                0.,   -st,    ct));
}

Rotation::Rotation(const Quaternion& q) {
    setToQuaternion(q);
}

// Orthogonalize the supplied matrix to make it a legitimate rotation.
// We do this by conversion to quaternions and back, although that may
// not produce the closest rotation.
Rotation::Rotation(const Mat33& m) {
    const Rotation bad(m, true);
    const Quaternion good(bad);
    this->setToQuaternion(good);
}

/*static*/Rotation Rotation::aboutAxis
   (const Real& angInRad, const UnitVec3& axis)
{
    Quaternion q; q.setToAngleAxis(angInRad, axis);
    return Rotation(q);
}

// Set this Rotation to represent the same rotation as
// the passed-in quaternion. The 0th element is the quaternion
// scalar.
void Rotation::setToQuaternion(const Quaternion& q) {
    const Real q00=q[0]*q[0], q11=q[1]*q[1], q22=q[2]*q[2], q33=q[3]*q[3];
    const Real q01=q[0]*q[1], q02=q[0]*q[2], q03=q[0]*q[3];
    const Real q12=q[1]*q[2], q13=q[1]*q[3], q23=q[2]*q[3];

    BaseMat::operator=( 
        BaseMat(q00+q11-q22-q33,   2.*(q12-q03)  ,   2.*(q13+q02),
                 2.*(q12+q03)  ,  q00-q11+q22-q33,   2.*(q23-q01),
                 2.*(q13-q02)  ,   2.*(q23+q01)  , q00-q11-q22+q33));
}

// Set this Rotation to represent a rotation of +q0 about
// the body frame's Z axis, followed by a rotation of +q1 about
// the body frame's NEW Y axis, followed by a rotation of +q3
// about the body frame's NEW X axis.
// See Kane, Spacecraft Dynamics, pg. 423, body-three: 3-2-1.
void Rotation::setToBodyFixed321(const Vec3& q) {
    const Real sq0 = std::sin(q[0]), cq0 = std::cos(q[0]);
    const Real sq1 = std::sin(q[1]), cq1 = std::cos(q[1]);
    const Real sq2 = std::sin(q[2]), cq2 = std::cos(q[2]);
    BaseMat::operator=( 
        BaseMat( cq0*cq1 , cq0*sq1*sq2-sq0*cq2 , cq0*sq1*cq2+sq0*sq2,
                 sq0*cq1 , sq0*sq1*sq2+cq0*cq2 , sq0*sq1*cq2-cq0*sq2,
                  -sq1   ,       cq1*sq2       ,      cq1*cq2        ));
}

// Set this Rotation to represent a rotation of +q0 about
// the body frame's X axis, followed by a rotation of +q1 about
// the body frame's NEW Y axis, followed by a rotation of +q3
// about the body frame's NEW Z axis.
// See Kane, Spacecraft Dynamics, pg. 423, body-three: 1-2-3.
void Rotation::setToBodyFixed123(const Vec3& q) {
    const Real sq0 = std::sin(q[0]), cq0 = std::cos(q[0]);
    const Real sq1 = std::sin(q[1]), cq1 = std::cos(q[1]);
    const Real sq2 = std::sin(q[2]), cq2 = std::cos(q[2]);
    BaseMat::operator=( 
        BaseMat(      cq1*cq2        ,       -cq1*sq2       ,  sq1    ,
                 sq0*sq1*cq2+cq0*sq2 , -sq0*sq1*sq2+cq0*cq2 , -sq0*cq1,
                -cq0*sq1*cq2+sq0*sq2 ,  cq0*sq1*sq2+sq0*cq2 ,  cq0*cq1 ));
}


// Convert this Rotation matrix to the equivalent quaternion. This
// is tricky to do without numerical errors. We use a modification
// of Richard Spurrier's method. See Spurrier, R.A., "Comment
// on 'Singularity-Free Extraction of a Quaternion from a
// Direction-Cosine Matrix'", J. Spacecraft and Rockets, 15(4):255, 1977.
// Our modification avoids all but one square root and divide.
// In each of the four cases we compute 4q[m]*q where m is the "max"
// element, with m=0 if the trace is larger than any diagonal or
// m=i if the i,i element is the largest diagonal and larger than
// the trace. Then when we normalize at the end the scalar 4q[m] 
// evaporates leaving us with a perfectly normalized quaternion.
// 
// The returned quaternion can be interpreted as a rotation angle
// a about a unit vector v=[vx vy vz] like this:
//    q = [ cos(a/2) sin(a/2)*v ]
// We canonicalize the returned quaternion by insisting that
// cos(a/2) >= 0, meaning that -180 < a <= 180.
//
Quaternion Rotation::convertToQuaternion() const {
    const Mat33& R = *this; // upcast
    const Real tr = R.trace();

    Vec4 q; // order is (cos sinx siny sinz)

    if (tr >= R(0,0) && tr >= R(1,1) && tr >= R(2,2)) {
        // trace is larger than any diagonal
        q[0] = 1 + tr;          // scalar element
        q[1] = R(2,1)-R(1,2);
        q[2] = R(0,2)-R(2,0);
        q[3] = R(1,0)-R(0,1);
    } else if (R(0,0) >= R(1,1) && R(0,0) >= R(2,2)) {
        // 0,0 element is largest
        q[0] = R(2,1)-R(1,2);   // scalar element
        q[1] = 1 - (tr - 2*R(0,0));
        q[2] = R(0,1)+R(1,0);
        q[3] = R(0,2)+R(2,0);
    } else if (R(1,1) >= R(2,2)) {
        // 1,1 element is largest
        q[0] = R(0,2)-R(2,0);   // scalar element
        q[1] = R(0,1)+R(1,0);
        q[2] = 1 - (tr - 2*R(1,1));
        q[3] = R(1,2)+R(2,1);
    } else {
        // 2,2 element is largest
        q[0] = R(1,0)-R(0,1);   // scalar element
        q[1] = R(0,2)+R(2,0);
        q[2] = R(1,2)+R(2,1);
        q[3] = 1 - (tr - 2*R(2,2));
    }
    Real scale = q.norm(); 
    if (q[0] < 0) scale = -scale; // canonicalize
    return Quaternion(q/scale, true); // prevent re-normalization
}

// Here is a rotation matrix made from a body012 (that is, 123) Euler 
// angle sequence (q0,q1,q2):
//
//         cq1*cq2        ,       -cq1*sq2       ,  sq1    ,
//    sq0*sq1*cq2+cq0*sq2 , -sq0*sq1*sq2+cq0*cq2 , -sq0*cq1,
//   -cq0*sq1*cq2+sq0*sq2 ,  cq0*sq1*sq2+sq0*cq2 ,  cq0*cq1
//
// We'll return Euler angles in the range -pi <= q0,q2 <= pi,
// and -pi/2 <= q1 <= pi/2.
//
// This computation is singular when the middle angle, q1, is +/-90 degrees.
// In that configuration the 1st & 3rd axes are aligned so we
// can arbitarily choose how to split that rotation between them.
// We'll choose q2=0 so that we'll have (q0, +/- pi/2,0).
//
// TODO: this routine will produce a pointing error of around
// sqrt(eps) radians for q1 within about 1e-8 of pi/2, that is,
// when sin(q1)==1. I think a more complicated algorithm, 
// perhaps iterative, could reduce the pointing error to machine
// precision. I tried making use of the fact that cos(q1) is well
// behaved there, using the 00,01 and 12,22 pairs to get q0 and q2
// (look at the matrix above). But since cos(q1) is close to zero
// when q1 is near pi/2, those terms can be indistinguishable
// from noise when sin(q1) is exactly 1, so the computation was
// too sensitive to junk there. Those terms are useful as soon
// as sin(q1) is anything but 1, because cos(q1) is already sqrt(eps)
// by then (see below). I also tried an iterative algorithm with
// a numerical Jacobian but it appears that an analytic (or complex
// step?) one is required due to the extreme sensitivity near pi/2.
// For now we'll just live with the occasional 1e-8 radian pointing
// error in conversion to singular Euler sequence, which really 
// isn't that bad! 
// (sherm 070415)
Vec3 Rotation::convertToBodyFixed123() const {
    const Rotation& R = *this;
    Real q0, q1, q2;

    const Real sq1 = R[0][2];

    // Any angle q1=pi/2 +/- sqrt(eps) will give sin(q1)==1, because
    // sin() is flat there. (Careful: sq1 might be a little greater
    // then 1 due to noise.) Numerically, we can get very close to 1
    // before we have to treat this as singular, because our calculation
    // of cos(q1) below will be well behaved.
    if (1-std::abs(sq1) <= 4*Eps)
    {
        // sq1==1, so we'll assume cq1==0 and we're going to set the last
        // angle to zero making sq2=0 and cq2=1, vastly simplifying the
        // middle column above.
        q0 = std::atan2(R[2][1],R[1][1]);
        q1 = sq1 > 0 ? Pi/2 : -Pi/2;
        q2 = 0;
    } else {
        // cq1 isn't zero; in fact it isn't much smaller than sqrt(eps)
        const Vec2 v0(R[2][2], -R[1][2]); // used for getting q0
        q0 = std::atan2(v0[1],v0[0]);

        const Vec2 v2(R[0][0], -R[0][1]); // used for getting q2
        q2 = std::atan2(v2[1],v2[0]);

        // Rather than using asin(sq1), we'll try to get a decent estimate of cq1
        // and then use the better behaved atan2(sq1,cq1). There are four terms
        // from which we can get cq1 in the rotation matrix (see above), two 
        // dependent on q0 and two dependent on q2. In either case squaring and
        // adding the two terms leaves us with just cq1^2, and since for
        // q1 in range -pi/2:pi/2, cos(q1) >= 0 we can just take the square root.
        // (Otherwise we'd lose a sign here.) To reduce noise sensitivity we'll
        // do it both ways and average the results.
        const Real cq1a = v0.norm(), cq1b = v2.norm();
        q1 = std::atan2(sq1, 0.5*(cq1a+cq1b));
    }

    return Vec3(q0,q1,q2);
}

// Here is a rotation matrix made from a body01 (that is, 12 or xy) Euler 
// angle sequence (q0,q1):
//
//        cy   ,    0   ,   sy   ,
//       sy*sx ,   cx   , -cy*sx ,
//      -sy*cx ,   sx   ,  cy*cx 
//
// The results can also be interpreted as a space10 (21, yx) sequence,
// meaning rotate by q1 about y, the q0 about the old x.
//
// This routine assumes that the current Rotation has the above structure,
// that is, that it is a rotation that could have resulted from the desired
// two-angle sequence. If not, the results will not be very meaningful.
//
// We'll return Euler angles in the range -pi <= q0,q1 <= pi. There are
// no singular configurations.
//
Vec2 Rotation::convertToBodyFixed12() const {
    const Rotation& R = *this;

    const Real q0 = std::atan2(R[2][1], R[1][1]);
    const Real q1 = std::atan2(R[0][2], R[0][0]);

    return Vec2(q0,q1);
}

// Here is a rotation matrix made from a space01 (that is, 12 or xy) Euler 
// angle sequence (q0,q1):
//
//        cy   ,  sy*sx ,  sy*cx ,
//         0   ,   cx   ,  -sx   ,
//       -sy   ,  cy*sx ,  cy*cx 
//
// The results can also be interpreted as a body10 (21, yx) sequence,
// meaning rotate by q1 about y, the q0 about the new x.
//
// This routine assumes that the current Rotation has the above structure,
// that is, that it is a rotation that could have resulted from the desired
// two-angle sequence. If not, the results will not be very meaningful.
//
// We'll return Euler angles in the range -pi <= q0,q1 <= pi. There are
// no singular configurations.
//
Vec2 Rotation::convertToSpaceFixed12() const {
    const Rotation& R = *this;

    const Real q0 = std::atan2(-R[1][2], R[1][1]);
    const Real q1 = std::atan2(-R[2][0], R[0][0]);

    return Vec2(q0,q1);
}

// The lazy but correct way. I have not seen a roundoff-safe method
// for direct conversion.
Vec4 Rotation::convertToAngleAxis() const {
    return convertToQuaternion().convertToAngleAxis();
}

// Call the current Rotation R_GB. We are passed R_GX and want to
// know if these are the same rotation, to within a specified "pointing error"
// in radians. We'll interpret that to mean that R_BX, when converted to
// angle-axis format, should have a rotation angle of no more than the given
// acceptable pointing error.
// We require 0 < pointingError < pi. 0 would require perfect accuracy; 180 would
// answer true for all possible rotation matrices.
bool Rotation::isSameRotationToWithinAngle(const Rotation& R_GX, Real okPointingError) const {
    assert(0 < okPointingError && okPointingError < Pi);
    const Rotation R_BX = ~(*this)*R_GX;
    const Real pointingError = std::abs(R_BX.convertToAngleAxis()[0]);
    return pointingError <= okPointingError;
}

// Call the current Rotation R_GB. We are passed R_GX and want to
// know if these are the same rotation, to within machine precision.
// We'll interpret that to mean that R_BX, when converted to
// angle-axis format, should have a rotation angle of no more than
// machine epsilon^(7/8), i.e. around 1e-14 radians in double precision.
bool Rotation::isSameRotationToMachinePrecision(const Rotation& R_GX) const {
    return isSameRotationToWithinAngle(R_GX, SignificantReal);
}


Quaternion::Quaternion(const Rotation& r) 
  : BaseVec(r.convertToQuaternion())
{ 
}

// DO NOT do the obvious acos(q[0]) to get the rotation angle!!!
// You will get numerical garbage anywhere near zero, and I don't
// mean all that near!
Vec4 Quaternion::convertToAngleAxis() const {
    const Real& ca2  = (*this)[0];      // cos(a/2)
    const Vec3& sa2v = getSubVec<3>(1); // sin(a/2)*v
    Real        sa2  = sa2v.norm();     // always >= 0

    // TODO: what is the right value to use here?? Norms can be
    // much less than eps and still OK -- this is 1e-32 in double.
    if (sa2 < square(Eps))
        return Vec4(0,1,0,0); // no rotation, x axis

    Vec4 av;

    // atan2 is numerically perfect. Since sa2>=0, it will return
    // an angle between 0 and pi, but because of the factor of 2
    // here we'll get angles between 0 and 2pi which we want to
    // pull into the -pi < a <= pi range (that is, instead of rotating say,
    // 359 degrees clockwise, rotate -1 degree counterclockwise.
    av[0] = 2*std::atan2(sa2,ca2);
    if (av[0] > Pi) av[0] -= 2*Pi;
    av.updSubVec<3>(1) = sa2v/sa2;
    return av;
}

// av = [ a vx vy vz ]
// If |a| < machine precision we'll treat this as zero rotation
// which will produce quaternion q=[1 0 0 0].
// Otherwise we'll insist that v have length at least machine precision,
// return NaN if not, and otherwise normalize it and use the result as the
// rotation axis.
void Quaternion::setToAngleAxis(const Vec4& av) {
    const Real eps = std::numeric_limits<Real>::epsilon();
    const Real& a = av[0];  // the angle
    if (std::fabs(a) < eps) {
        BaseVec::operator=(Vec4(1,0,0,0));
        return;
    }
    const Vec3& vIn = av.getSubVec<3>(1);
    const Real vnorm = vIn.norm();
    if (vnorm < eps) setToNaN();
    else setToAngleAxis(a, UnitVec3(vIn/vnorm, true));
}

// Here there can be no problems. The angle can be anything, but
// the quaternion will effectively reduce it to the -pi < a <= pi
// range, meaning the scalar part of the quaternion (cos(a/2))
// will be nonnegative.
void Quaternion::setToAngleAxis(const Real& a, const UnitVec3& v) {
    Real ca2 = std::cos(0.5*a), sa2 = std::sin(0.5*a);

    // If ca2 < 0 we have 90 < |a/2| < 180. We can move that
    // to 0 < |a/2| < 90 by adding or subtracting 180 to a/2,
    // which is (a +/- 360), changing nothing. That will have
    // the effect of negating both ca2 & sa2.
    if (ca2 < 0) 
        ca2 = -ca2, sa2 = -sa2;
    // OK, we now have -90 <= a/2 <= 90, so -180 <= a <= 180.
    (*this)[0] = ca2;
    (*this).updSubVec<3>(1) = sa2*v;
}


std::ostream& operator<<(std::ostream& o, const Rotation& m) {
    return o << m.asMat33();
}

std::ostream& operator<<(std::ostream& o, const Transform& x) {
    return o << x.asMat34() << Row4(0,0,0,1) << std::endl;
}

} // namespace SimTK

