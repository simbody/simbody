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
    /// Default constructor produces the ZeroRotation quaternion [1 0 0 0] (not NaN - even in debug mode).
    Quaternion() : Vec4(1,0,0,0) { }

    /// Zero-cost copy constructor just copies the source without conversion to canonical form or normalization.
    Quaternion( const Quaternion& q ) : Vec4(q) {}

    /// Zero-cost copy assignment just copies the source without conversion to canonical form or normalization.
    Quaternion& operator=( const Quaternion& q ) { Vec4::operator=( q.asVec4() );  return *this; }
    
    /// Construct a quaternion and normalize it 
    Quaternion( Real e0, Real e1, Real e2, Real e3 ) : Vec4(e0,e1,e2,e3)  { normalizeThis(); }
    explicit Quaternion( const Vec4& q ) : Vec4(q)                        { normalizeThis(); }

    /// Constructs a canonical quaternion from a rotation matrix (cost is about 60 flops).
    SimTK_SimTKCOMMON_EXPORT explicit Quaternion( const Rotation& );

    /// The ZeroRotation quaternion is [1 0 0 0].
    /// Note: Default constructor is ZeroRotation (unlike Vec4 which start as NaN in Debug mode).
    void setQuaternionToZeroRotation()  { Vec4::operator=( Vec4(1,0,0,0) ); }
    void setQuaternionToNaN()           { Vec4::setToNaN(); }

    /// The quaternion that is set by this method has a non-negative first element (canonical form).
    /// If the "axis" portion of av is a zero vector, the quaternion is set to all-NaN.
    SimTK_SimTKCOMMON_EXPORT void  setQuaternionFromAngleAxis( const Vec4& av );
    SimTK_SimTKCOMMON_EXPORT void  setQuaternionFromAngleAxis( const Real& a, const UnitVec3& v );

    /// Returns [ a vx vy vz ] with (a,v) in canonical form, i.e., -180 < a <= 180 and |v|=1. 
    SimTK_SimTKCOMMON_EXPORT Vec4  convertQuaternionToAngleAxis() const;

    /// Zero-cost cast of a Quaternion to a Vec4.
    const Vec4&  asVec4() const  { return *static_cast<const Vec4*>(this); }

    /// Normalize an already constructed quaternion.
    /// If the quaternion is *exactly* zero, set it to [1 0 0 0].
    /// If its magnitude is:  0 < magnitude < epsilon  (epsilon is machine tolerance), set it to NaN (treated as an error). 
    /// Otherwise, normalize the quaternion which costs about 40 flops.
    /// The quaternion is NOT put in canonical form.
    Quaternion&  normalizeThis() { 
        const Real epsilon = std::numeric_limits<Real>::epsilon();
        const Real magnitude = Vec4::norm();
        if(      magnitude == 0      )  setQuaternionToZeroRotation();
        else if( magnitude < epsilon )  setQuaternionToNaN();
		else (*this) *= (1.0/magnitude);
        return *this;
    }

    /// Use this constructor only if you are *sure* v is normalized to 1.0.
    /// This zero cost method is faster than the Quaternion(Vec4) constructor which normalizes the Vec4. 
    /// The second argument forces the compiler to call the fast constructor; it is otherwise ignored. 
    /// By convention, set the second argument to "true". 
    Quaternion( const Vec4& v, bool ) : Vec4(v) {}

//----------------------------------------------------------------------------------------------------
// The following code is obsolete - it is here temporarily for backward compatibility (Mitiguy 10/13/2007)
//----------------------------------------------------------------------------------------------------
public:
    Vec4  convertToAngleAxis() const                          { return convertQuaternionToAngleAxis(); }
    void  setToAngleAxis( const Vec4& av )                    { setQuaternionFromAngleAxis(av); }
    void  setToAngleAxis( const Real& a, const UnitVec3& v )  { setQuaternionFromAngleAxis(a,v); }
    void setToNaN()                                           { setQuaternionToNaN(); }
    void setToZero()                                          { setQuaternionToZeroRotation(); }

};


//------------------------------------------------------------------------------
}  // End of namespace SimTK

//--------------------------------------------------------------------------
#endif // SimTK_QUATERNION_H_
//--------------------------------------------------------------------------

