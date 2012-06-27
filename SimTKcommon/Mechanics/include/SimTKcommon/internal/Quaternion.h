#ifndef SimTK_SimTKCOMMON_QUATERNION_H 
#define SimTK_SimTKCOMMON_QUATERNION_H 

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman, Paul Mitiguy                                     *
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

//-----------------------------------------------------------------------------
#include "SimTKcommon/SmallMatrix.h"
#include "SimTKcommon/internal/UnitVec.h"
//-----------------------------------------------------------------------------
#include <iosfwd>  // Forward declaration of iostream
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
namespace SimTK {

//-----------------------------------------------------------------------------
// Forward declarations
template <class P> class Rotation_;
template <class P> class Quaternion_;

typedef Quaternion_<Real>    Quaternion;
typedef Quaternion_<float>  fQuaternion;
typedef Quaternion_<double> dQuaternion;

//-----------------------------------------------------------------------------
/**
 * A Quaternion is a Vec4 with the following behavior:
 *   - its length is always 1 (or else it is all NaN)
 *   - it is equivalent to an angle/axis rotation for
 *     angle a, axis unit vector v, as:  q = [ cos(a/2) sin(a/2)*v ]
 * A quaternion is in "canonical form" when its first element is nonnegative. 
 * This corresponds to rotation angles in the range -180 < a <= 180 degrees. 
 * Quaternions are not required to be in canonical form (e.g., during numerical
 * integration). When appropriate, they are put in canonical form.
 *
 * Conversion from quaternion to (angle,axis) form is handled here also. 
 * (angle,axis) is in canonical form when -180 < angle <= 180 and |axis|=1.
 * However, (angle,axis) is meaningful for any angle and for any axis where 
 * |axis| > 0.
 */
//-----------------------------------------------------------------------------
template <class P>
class Quaternion_ : public Vec<4,P> {
    typedef P           RealP;
    typedef Vec<3,P>    Vec3P;
    typedef Vec<4,P>    Vec4P;
public:
    /// Default constructor produces the ZeroRotation quaternion [1 0 0 0] 
    /// (not NaN - even in debug mode).
    Quaternion_() :  Vec4P(1,0,0,0) { }

    /// Zero-cost copy constructor just copies the source without conversion to 
    /// canonical form or normalization.
    Quaternion_(const Quaternion_& q) :  Vec4P(q) {}

    /// Zero-cost copy assignment just copies the source without conversion to 
    /// canonical form or normalization.
    Quaternion_& operator=( const Quaternion_& q ) 
    {   Vec4P::operator=(q.asVec4()); return *this; }
    
    /// Construct a quaternion from four scalars and normalize the result,
    /// which costs about 40 flops.
    Quaternion_( RealP e0, RealP e1, RealP e2, RealP e3 ) : Vec4P(e0,e1,e2,e3) 
    {   normalizeThis(); }
    /// Construct a quaternion from a 4-vector and normalize the result,
    /// which costs about 40 flops.
    explicit Quaternion_( const Vec4P& q ) : Vec4P(q) 
    {   normalizeThis(); }

    /// Constructs a canonical quaternion from a rotation matrix (cost is 
    /// about 60 flops).
    SimTK_SimTKCOMMON_EXPORT explicit Quaternion_(const Rotation_<P>&);

    /// The ZeroRotation quaternion is [1 0 0 0]. Note: Default constructor 
    /// is ZeroRotation (unlike Vec4P which start as NaN in Debug mode).
    void setQuaternionToZeroRotation()  {Vec4P::operator=( Vec4P(1,0,0,0) );}
    /// Set quaternion to all-NaN. Note that this is not the same as produced
    /// by default construction, even in Debug mode -- default construction
    /// always produces an identity rotation of [1 0 0 0].
    void setQuaternionToNaN()           {Vec4P::setToNaN();}

    /// Set this quaternion from an angle-axis rotation packed into a 4-vector
    /// as [a vx vy vz]. The result will be put in canonical form, i.e., 
    /// it will have a non-negative first element. If the "axis" portion of av 
    /// is a zero vector on input, the quaternion is set to all-NaN.
    SimTK_SimTKCOMMON_EXPORT void  setQuaternionFromAngleAxis(const Vec4P& av);
    /// Set this quaternion from an angle-axis rotation provided as an angle a
    /// and a separate unit vector [vx vy vz]. The result will be put in 
    /// canonical form, i.e., it will have a non-negative first element.
    SimTK_SimTKCOMMON_EXPORT void  setQuaternionFromAngleAxis(const RealP& a, const UnitVec<P,1>& v);

    /// Returns [ a vx vy vz ] with (a,v) in canonical form, i.e., 
    /// -180 < a <= 180 and |v|=1. 
    SimTK_SimTKCOMMON_EXPORT Vec4P convertQuaternionToAngleAxis() const;

    /// Zero-cost cast of a Quaternion_ to its underlying Vec4; this is \e not
    /// converted to axis-angle form.
    const Vec4P& asVec4() const  {return *static_cast<const Vec4P*>(this);}

    /// Normalize an already constructed quaternion in place; but do you really
    /// need to do this? Quaternions should be kept normalized at all times. 
    /// One of the advantages of using them is that you don't have to check if 
    /// they are normalized or renormalize them. However, under some situations
    /// they do need renormalization, but it is costly if you don't actually 
    /// need it. If the quaternion is \e exactly zero, set it to [1 0 0 0]. If 
    /// its magnitude is 0 < magnitude < epsilon  (epsilon is machine 
    /// tolerance), set it to NaN (treated as an error). Otherwise, normalize 
    /// the quaternion which costs about 40 flops. The quaternion is NOT put 
    /// in canonical form.
    Quaternion_& normalizeThis() { 
        const RealP epsilon = std::numeric_limits<RealP>::epsilon();
        const RealP magnitude = Vec4P::norm();
        if(      magnitude == 0      )  setQuaternionToZeroRotation();
        else if( magnitude < epsilon )  setQuaternionToNaN();
        else (*this) *= (1/magnitude);
        return *this;
    }

    /// Return a normalized copy of this quaternion; but do you really need to
    /// do this? Quaternions should be kept normalized at all times. One of
    /// the advantages of using them is that you don't have to check if they
    /// are normalized or renormalize them. However, under some situations they
    /// do need renormalization, but it is costly if you don't actually need it.
    /// @see normalizeThis() for details.
    Quaternion_ normalize() const {return Quaternion_(*this).normalizeThis();}

    /// Use this constructor only if you are *sure* v is normalized to 1.0.
    /// This zero cost method is faster than the Quaternion(Vec4) constructor 
    /// which normalizes the Vec4. The second argument forces the compiler to 
    /// call the fast constructor; it is otherwise ignored. By convention, set 
    /// the second argument to "true". 
    Quaternion_(const Vec4P& v, bool) : Vec4P(v) {}

//----------------------------------------------------------------------------------------------------
// The following code is obsolete - it is here temporarily for backward compatibility (Mitiguy 10/13/2007)
//----------------------------------------------------------------------------------------------------
private:
    Vec4P convertToAngleAxis() const                          { return convertQuaternionToAngleAxis();}
    void setToAngleAxis(const Vec4P& av)                      { setQuaternionFromAngleAxis(av);}
    void setToAngleAxis(const RealP& a, const UnitVec<P,1>& v){ setQuaternionFromAngleAxis(a,v);}
    void setToNaN()                                           { setQuaternionToNaN(); }
    void setToZero()                                          { setQuaternionToZeroRotation();}

};


//------------------------------------------------------------------------------
}  // End of namespace SimTK

//--------------------------------------------------------------------------
#endif // SimTK_SimTKCOMMON_QUATERNION_H_
//--------------------------------------------------------------------------

