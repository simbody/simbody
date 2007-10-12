//-----------------------------------------------------------------------------
// File:     Quaternion.h
// Class:    Quaternion 
// Parent:   Vec4
// Purpose:  Quaternion (Euler parameters) for representing orientation
//-----------------------------------------------------------------------------
#ifndef SimTK_QUATERNION_H 
#define SimTK_QUATERNION_H 

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmatrix(tm)                       *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-7 Stanford University and the Authors.         *
 * Authors: Michael Sherman and Paul Mitiguy                                  *
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

//-----------------------------------------------------------------------------
#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/Constants.h"
#include "SimTKcommon/internal/Scalar.h"
#include "SimTKcommon/internal/SmallMatrix.h"
#include "SimTKcommon/internal/UnitVec.h"
//-----------------------------------------------------------------------------
#include <iosfwd>  // Forward declaration of iostream
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
namespace SimTK {

//-----------------------------------------------------------------------------
// Forward declarations
class Rotation;


//-----------------------------------------------------------------------------
/**
 * A Quaternion is a Vec4 with the following behavior:
 *   - its length is always 1 (or else it is all NaN)
 *   - it is equivalent to an angle/axis rotation for
 *     angle a, axis unit vector v, as:  q = [ cos(a/2) sin(a/2)*v ]
 * A quaternion is in "canonical form" when its first element is nonnegative. 
 * This corresponds to rotation angles in the range -180 < a <= 180 degrees. 
 * Quaternions are not required to be in canonical form (e.g., during numerical integration).
 * When appropriate, they are put in canonical form.
 *
 * Conversion from quaternion to (angle,axis) form is handled here also. 
 * (angle,axis) is in canonical form when -180 < angle <= 180 and |axis|=1.
 * However, (angle,axis) is meaningful for any angle and for any axis where |axis| > 0.
 */
//-----------------------------------------------------------------------------
class Quaternion : public Vec4 {
public:
    typedef Vec4 BaseVec;

    /// Default constructor produces the ZeroRotation quaternion [1 0 0 0].
    Quaternion() : Vec4(1,0,0,0) { }

    /// Construct a quaternion from a Vec4 v which is interpreted as a
    /// quaternion that needs normalization [NOT (angle,axis)].
    /// If the passed-in vector v is *exactly* zero, the quaternion is set to [1 0 0 0].
    /// If the length of v is 0 < len < eps (eps being machine tolerance),
    /// the quaternion is set to NaN (treated as an error). 
    /// Otherwise, the quaternion is set by normalizing v (40 flops).
    /// The constructed quaternion is NOT put in canonical form.
    explicit Quaternion( const Vec4& v ) {
        const Real eps = std::numeric_limits<Real>::epsilon();
        const Real len = v.norm();
        if      (len == 0)  setQuaternionToZeroRotation();
        else if (len < eps) setQuaternionToNaN();
        else BaseVec::operator=( v/len );
    }

    /// Constructs a canonical quaternion from a rotation matrix (cost is about 60 flops).
    SimTK_SimTKCOMMON_EXPORT explicit Quaternion( const Rotation& );

    /// Zero-cost copy constructor just copies the source without conversion to canonical form or normalization.
    Quaternion( const Quaternion& q ) : BaseVec(q) { }

    /// Zero-cost copy assignment just copies the source without conversion to canonical form or normalization.
    Quaternion& operator=( const Quaternion& q ) { BaseVec::operator=( q.asVec4() );  return *this; }
    
    /// The ZeroRotation quaternion is [1 0 0 0].
    void setQuaternionToZeroRotation()  { BaseVec::operator=( Vec4(1,0,0,0) ); }
    void setToZero()                    { setQuaternionToZeroRotation(); }

    /// This is the only exception to the "must be normalized" rule for quaternions. 
    /// Note: Unlike naked Vec4's, Quaternions do not start out NaN even in Debug mode.
    /// The default constructor sets the Quaternion to "zero rotation" instead.
    void setQuaternionToNaN() { BaseVec::setToNaN(); }
    void setToNaN()           { setQuaternionToNaN(); }

    /// Resulting 4-vector is [ a vx vy vz ] with (a,v) in canonical form.
    /// That is, -180 < a <= 180 and |v|=1. The cost of this operation is
    /// roughly one atan2, one sqrt, and one divide, say about 100 flops.
    SimTK_SimTKCOMMON_EXPORT Vec4 convertQuaternionToAngleAxis() const;
    Vec4 convertToAngleAxis() const  { return convertQuaternionToAngleAxis(); }

    /// The quaternion that is set by this method has a non-negative first element (canonical form).
    /// If the "axis" portion of av is a zero vector, the quaternion is set to all-NaN.
    SimTK_SimTKCOMMON_EXPORT void setQuaternionFromAngleAxis( const Vec4& av );
    void setToAngleAxis( const Vec4& av ) { setQuaternionFromAngleAxis(av); }

    /// The quaternion that is set by this method has a non-negative first element (canonical form).
    SimTK_SimTKCOMMON_EXPORT void setQuaternionFromAngleAxis( const Real& a, const UnitVec3& v );
    void setToAngleAxis( const Real& a, const UnitVec3& v )  { setQuaternionFromAngleAxis(a,v); }

    /// Upcast this Quaternion to its parent class, a Vec4. This is inline
    /// and should generate no code. You can do the same thing with static_cast
    /// if you prefer. Zero cost.
    const BaseVec& asVec4() const  { return *static_cast<const BaseVec*>(this); }

    /// Use this method only if you are *sure* v is normalized to 1.0.
    /// This zero cost method is faster than the Quaternion(Vec4) constructor 
    /// which normalizes the Vec4. The second argument forces the compiler to call 
    //  the fast constructor; it is otherwise ignored. 
    /// By convention, set the second argument to "true". 
    Quaternion( const Vec4& v, bool ) : Vec4(v) {}
};


//------------------------------------------------------------------------------
}  // End of namespace SimTK

//--------------------------------------------------------------------------
#endif // SimTK_QUATERNION_H_
//--------------------------------------------------------------------------

