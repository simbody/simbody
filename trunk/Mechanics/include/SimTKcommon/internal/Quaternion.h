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
 * Authors: Michael Sherman                                                   *
 * Contributors: Paul Mitiguy                                                 *
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
#include "SimTKcommon/internal/Orientation.h"

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
 *     angle a, axis unit vector v, like this:
 *        q=[ cos(a/2) sin(a/2)*v ]
 * We consider a quaternion to be in "canonical form" when its
 * first element is nonnegative. That corresponds to rotation
 * angles in the range -180 < a <= 180 degrees. We don't require
 * quaternions to be in canonical form; continuity during integration
 * requires them to range more widely. However, when we're creating
 * them from scratch and have a choice, we'll return them in 
 * canonical form.
 *
 * The (angle,axis) form is handled here also. When these are in
 * their canonical form, they have -180 < angle <= 180 and |axis|=1.
 * However, (angle,axis) is meaningful for any value of angle and
 * for any axis where |axis| > 0.
 */
//-----------------------------------------------------------------------------
class Quaternion : public Vec4 {
public:
    typedef Vec4 BaseVec;

    /// Default constructor produces a 0-rotation (identity) quaternion
    Quaternion() : Vec4(1,0,0,0) { }

    /// Initialize this quaternion from an unnormalized quaternion
    /// stored in a Vec4. The argument is *not* interpreted as an
    /// (angle,axis) -- it is just a slightly mangled quaternion.
    /// If the passed-in vector is *exactly* zero, we will assume
    /// it indicates a "zero rotation" and set the Quaternion to
    /// [1 0 0 0]. If the length is 0 < len < eps (eps being machine
    /// tolerance) we consider that an error condition and set the
    /// Quaternion to NaN. Otherwise we normalize it and return.
    /// The constructed quaternion is NOT put in canonical form -- it is
    /// as close to the original as possible.
    /// Because of the normalization, this costs about 40 flops.
    explicit Quaternion( const Vec4& v ) {
        const Real eps = std::numeric_limits<Real>::epsilon();
        const Real len = v.norm();
        if      (len == 0)  setToZero();
        else if (len < eps) setToNaN();
        else BaseVec::operator=( v/len );
    }

    /// Initialize this quaternion from a rotation matrix. The result
    /// will be in canonical form. The cost is about 60 flops.
    SimTK_SimTKCOMMON_EXPORT explicit Quaternion( const Rotation& );

    /// Copy constructor copies the source as-is; it does not 
    /// convert to canonical form, or normalize, or anything else.
    /// Zero cost.
    Quaternion( const Quaternion& q ) : BaseVec(q) { }

    /// Copy assignment copies the source as-is; it does not 
    /// convert to canonical form, or normalize, or anything else.
    /// Zero cost.
    Quaternion& operator=( const Quaternion& q ) { BaseVec::operator=( q.asVec4() );  return *this; }
    
    /// By zero here we mean "zero rotation", i.e., an identity rotation
    /// represented as [1 0 0 0]. This is in canonical form; [-1 0 0 0] would
    /// mean the same thing.
    void setToZero() { BaseVec::operator=( Vec4(1,0,0,0) ); }

    /// This is the only exception to the "must be normalized" rule for
    /// quaternions -- all elements are set to NaN. Note that unlike 
    /// naked Vec4's, Quaternions do not start out NaN even in Debug mode.
    /// The default constructor sets the Quaternion to "zero rotation" instead.
    void setToNaN() { BaseVec::setToNaN(); }

    /// Resulting 4-vector is [ a vx vy vz ] with (a,v) in canonical form.
    /// That is, -180 < a <= 180 and |v|=1. The cost of this operation is
    /// roughly one atan2, one sqrt, and one divide, say about 100 flops.
    SimTK_SimTKCOMMON_EXPORT Vec4 convertToAngleAxis() const;

    /// Assign the current quaternion to the rotation represented by the 
    /// passed-in (angle,axis) form. The resulting quaternion will be in
    /// canonical form regardless of the condition of the (angle,axis) input.
    /// The "axis" will be normalized here unless it has zero length on
    /// entry, in which case the quaternion will be all-NaN.
    /// Cost is a normalization, a sin and a cos, or about 120 flops.
    SimTK_SimTKCOMMON_EXPORT void setToAngleAxis( const Vec4& av );

    /// Assign the current quaternion to the rotation represented by the 
    /// passed-in (angle,unitVector) form. The resulting quaternion will be in
    /// canonical form regardless of the condition of the (angle,axis) input.
    /// This can't fail for any angle since we know we have a good axis.
    /// Cost is one sin, one cos or roughly 80 flops.
    SimTK_SimTKCOMMON_EXPORT void setToAngleAxis( const Real& a, const UnitVec3& v );

    /// Upcast this Quaternion to its parent class, a Vec4. This is inline
    /// and should generate no code. You can do the same thing with static_cast
    /// if you prefer. Zero cost.
    const BaseVec& asVec4() const  { return *static_cast<const BaseVec*>(this); }

    /// Don't use this unless you are *sure* this is already normalized! This
    /// is much faster than the normal Quaternion(Vec4) constructor which
    /// expects the Vec4 to need cleaning up. The second argument here is just
    /// to allow you to force a call to the fast constructor; it is otherwise
    /// ignored. By convention however, you should call this with the second
    /// argument set to "true". Zero cost.
    Quaternion( const Vec4& v, bool ) : Vec4(v) { }
};


//------------------------------------------------------------------------------
}  // End of namespace SimTK

//--------------------------------------------------------------------------
#endif // SimTK_QUATERNION_H_
//--------------------------------------------------------------------------

