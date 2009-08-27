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
 * Portions copyright (c) 2005-9 Stanford University and the Authors.         *
 * Authors: Michael Sherman and Paul Mitiguy                                  *
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
 * Implementations of non-inline methods of classes dealing with quaternions.
 */

#include "SimTKcommon/basics.h"
#include "SimTKcommon/internal/Quaternion.h"
#include "SimTKcommon/internal/Rotation.h"

//------------------------------------------------------------------------------
namespace SimTK {


//------------------------------------------------------------------------------
// Constructs a canonical quaternion from a rotation matrix (cost is ~60 flops).
//------------------------------------------------------------------------------
template <class P>
Quaternion_<P>::Quaternion_(const Rotation_<P>& r) 
:   Vec<4,P>(r.convertRotationToQuaternion()) {}


//------------------------------------------------------------------------------
// Returns [ a vx vy vz ] with (a,v) in canonical form, i.e., -180 < a <= 180 and |v|=1.
// The cost of this operation is roughly one atan2, one sqrt, and one divide (about 100 flops).
//------------------------------------------------------------------------------
template <class P> Vec<4,P>
Quaternion_<P>::convertQuaternionToAngleAxis() const {
    const RealP& ca2  = (*this)[0];       // cos(a/2)
    const Vec3P& sa2v = this->getSubVec<3>(1);  // sin(a/2) * v
    RealP        sa2  = sa2v.norm();      // sa2 is always >= 0

    const RealP Eps = NTraits<P>::getEps();
    const RealP Pi  = NTraits<P>::getPi();

    // TODO: what is the right value to use here?? Norms can be
    // much less than eps and still OK -- this is 1e-32 in double.
    if( sa2 < square(Eps) )  return Vec4P(0,1,0,0); // no rotation, x axis

    // Use atan2.  Do NOT just use acos(q[0]) to calculate the rotation angle!!!
    // Otherwise results are numerical garbage anywhere near zero (or less near).
    RealP angle = 2 * std::atan2(sa2,ca2);

    // Since sa2>=0, atan2 returns a value between 0 and pi, which is then
    // multiplied by 2 which means the angle is between 0 and 2pi.
    // We want an angle in the range:  -pi < angle <= pi range.
    // E.g., instead of rotating 359 degrees clockwise, rotate -1 degree counterclockwise.
    if( angle > Pi ) angle -= 2*Pi;

    // Normalize the axis part of the return value
    const Vec3P axis = sa2v / sa2;

	// Return (angle/axis)
    return Vec4P( angle, axis[0], axis[1], axis[2] );
}


//-------------------------------------------------------------------
template <class P> void
Quaternion_<P>::setQuaternionFromAngleAxis( const Vec4P& av ) {
    // av = [ a vx vy vz ]
    // If |a| < machine precision,  we treat as a zero rotation which produces quaternion q=[1 0 0 0].
    const RealP eps = std::numeric_limits<RealP>::epsilon();
    const RealP& a = av[0];  // the angle
    if( std::fabs(a) < eps ) { Vec4P::operator=( Vec4P(1,0,0,0) );  return; }

    // The vector v must have length at least machine precision (or return NaN).
    const Vec3P& vIn = av.getSubVec<3>(1);
    const RealP vnorm = vIn.norm();
    if( vnorm < eps ) setQuaternionToNaN();

    // Otherwise, the vector v is normalized and used as the rotation axis.
    // Note: The cost of this method is 120 flops, including normalization (about 40 flops)
    //       and the sine and cosine calculations (80 flops) used in the next method.
    else setToAngleAxis( a, UnitVec<P,1>(vIn/vnorm, true) );
}


//-------------------------------------------------------------------
template <class P> void
Quaternion_<P>::setQuaternionFromAngleAxis( const RealP& a, const UnitVec<P,1>& v ) {
    /// The cost of this method is approximately 80 flops (one sin and one cos).
    RealP ca2 = std::cos(a/2), sa2 = std::sin(a/2);

    // Multiplying an entire quaternion by -1 produces the same Rotation matrix
    // (each element of the Rotation element involves the product of two quaternion elements).
    // The canonical form is to make the first element of the quaternion positive.
    if( ca2 < 0 ) { ca2 = -ca2; sa2 = -sa2; }
    (*this)[0] = ca2;
    (*this).updSubVec<3>(1) = sa2*v;
}

// Instantiate now to catch bugs.
template class Quaternion_<float>;
template class Quaternion_<double>;

//------------------------------------------------------------------------------
}  // End of namespace SimTK


