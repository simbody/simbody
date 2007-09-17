//-----------------------------------------------------------------------------
// File:     Quaternion.h
// Class:    Quaternion
// Parent:   Vec4
// Purpose:  Quaternion (Euler parameters) for representing orientation
//-----------------------------------------------------------------------------

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
#include "SimTKcommon/internal/Orientation.h"

//-------------------------------------------------------------------
namespace SimTK {


Quaternion::Quaternion( const Rotation& r ) : BaseVec(r.convertToQuaternion()) { }


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
void  Quaternion::setToAngleAxis( const Vec4& av ) {
    const Real eps = std::numeric_limits<Real>::epsilon();
    const Real& a = av[0];  // the angle
    if (std::fabs(a) < eps) {
        BaseVec::operator=( Vec4(1,0,0,0) );
        return;
    }
    const Vec3& vIn = av.getSubVec<3>(1);
    const Real vnorm = vIn.norm();
    if (vnorm < eps) setToNaN();
    else setToAngleAxis( a, UnitVec3(vIn/vnorm, true) );
}


// Here there can be no problems. The angle can be anything, but
// the quaternion will effectively reduce it to the -pi < a <= pi
// range, meaning the scalar part of the quaternion (cos(a/2))
// will be nonnegative.
void  Quaternion::setToAngleAxis( const Real& a, const UnitVec3& v ) {
    Real ca2 = std::cos(0.5*a), sa2 = std::sin(0.5*a);

    // If ca2 < 0 we have 90 < |a/2| < 180. We can move that
    // to 0 < |a/2| < 90 by adding or subtracting 180 to a/2,
    // which is (a +/- 360), changing nothing. That will have
    // the effect of negating both ca2 & sa2.
    if( ca2 < 0 )
        ca2 = -ca2, sa2 = -sa2;
    // OK, we now have -90 <= a/2 <= 90, so -180 <= a <= 180.
    (*this)[0] = ca2;
    (*this).updSubVec<3>(1) = sa2*v;
}



//------------------------------------------------------------------------------
}  // End of namespace SimTK


