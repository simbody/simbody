//-----------------------------------------------------------------------------
// File:     UnitVec.h
// Classes:  UnitVec and UnitRow
// Parents:  Vec and Row
// Purpose:  Unit vector class (pure direction - magnitude is always 1.0)
//-----------------------------------------------------------------------------
#ifndef SimTK_UNITVEC_H 
#define SimTK_UNITVEC_H 

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmatrix(tm)                       *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-9 Stanford University and the Authors.         *
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
#include "SimTKcommon/SmallMatrix.h"
//-----------------------------------------------------------------------------
#include <iosfwd>  // Forward declaration of iostream
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
namespace SimTK {

//-----------------------------------------------------------------------------
// Forward declarations. These are templatized by precision P and and stride S
// but always have length 3. TODO: this should be generalized to other lengths.
template <class P, int S> class UnitVec;
template <class P, int S> class UnitRow;

// UnitVec3 is more intelligible name for UnitVec<Real,1>.
typedef UnitVec<Real,1> UnitVec3;


//-----------------------------------------------------------------------------
/**
 * This class is a Vec3 plus an ironclad guarantee either that:
 *      - the length is one (to within a very small tolerance), or
 *      - all components are NaN.
 */
//-----------------------------------------------------------------------------
template <class P, int S>
class UnitVec : public Vec<3,P,S> {
    typedef P   RealP;
public:
    typedef Vec<3,P,S>      BaseVec;
    typedef UnitRow<P,S>    TransposeType;

    /// Default constructor initializes to all-NaN even in Release mode so that
    /// we maintain the above-promised contract.
    UnitVec() : BaseVec(NTraits<P>::getNaN()) {}

    /// Copy constructor does not require normalization since we know the
    /// source is a unit vector.
    UnitVec(const UnitVec& u) 
    :   BaseVec( static_cast<const BaseVec&>(u) ) {}

    /// Automatic conversion from UnitVec with different stride; no computation
    /// required.
    template <int S2> UnitVec(const UnitVec<P,S2>& u) 
    :   BaseVec( static_cast<const typename UnitVec<P,S2>::BaseVec&>(u) ) {}

    /// Explicit conversion from Vec to UnitVec, requiring expensive normalization.
    explicit UnitVec(const BaseVec& v) : BaseVec(v/v.norm()) {}
    /// Explicit conversion from Vec of any stride to this UnitVec, requiring 
    /// expensive normalization.
    template <int S2> 
    explicit UnitVec(const Vec<3,P,S2>& v) : BaseVec(v/v.norm())  {}

    /// Create a unit vector in the direction of the vector (x,y,z) whose measure
    /// numbers are supplied -- this requires an expensive normalization since
    /// we don't know that the supplied vector is normalized.
    UnitVec(const RealP& x, const RealP& y, const RealP& z) : BaseVec(x,y,z)  
    {   static_cast<BaseVec&>(*this) /= BaseVec::norm(); }

    /// Construct a unit axis vector 100 010 001 given 0,1, or 2.
    explicit UnitVec(int axis) : BaseVec(0) 
    {   assert(0 <= axis && axis <= 2);
        BaseVec::operator[](axis) = 1; }

    /// Copy assignment does not require normalization.
    UnitVec& operator=(const UnitVec& u) 
    {   BaseVec::operator=(static_cast<const BaseVec&>(u)); 
        return *this; }

    /// Copy assignment from a UnitVec whose stride differs from this one; no
    /// normalization required.
    template <int S2> UnitVec& operator=(const UnitVec<P,S2>& u) 
    {   BaseVec::operator=(static_cast<const typename UnitVec<P,S2>::BaseVec&>(u));
        return *this; }

    /// Return a reference to the underlying Vec3 (no copying here).
    const BaseVec& asVec3() const {return static_cast<const BaseVec&>(*this);}

    // Override Vec3 methods which preserve length. These return a 
    // packed UnitVec regardless of our stride.

    /// Returns a new unit vector pointing in the opposite direction from this one;
    /// does \e not modify this UnitVec object. Cost is 3 flops.
    UnitVec<P,1> negate()    const {return UnitVec<P,1>(-asVec3(),true);}
    /// Returns a new unit vector pointing in the opposite direction from this one.
    /// Cost is 3 flops.
    UnitVec<P,1> operator-() const {return negate();}

    /// Return a const reference to this unit vector re-expressed as a unit row; no
    /// computational cost.
    const TransposeType& operator~() const {return *reinterpret_cast<const TransposeType*>(this);}
    /// Return a writable reference to this unit vector re-expressed as a unit row; no
    /// computational cost.
    TransposeType& operator~() {return *reinterpret_cast<TransposeType*>(this);}

    // We have to define these here so that the non-const ones won't be
    // inherited. We don't trust anyone to write on one element of a UnitVec!

    /// Return one element of this unit vector as a const reference; there is no 
    /// corresponding writable index function since changing a single element of
    /// a unit vector would violate the contract that it has unit length at all times.
    const RealP&  operator[](int i) const  { return BaseVec::operator[](i); }
    /// Return one element of this unit vector as a const reference; there is no 
    /// corresponding writable index function since changing a single element of
    /// a unit vector would violate the contract that it has unit length at all times.
    const RealP&  operator()(int i) const  { return BaseVec::operator()(i); }

    /// Return a new unit vector whose measure numbers are the absolute values
    /// of the ones here. This will still have unit length but will be
    /// a reflection of this unit vector into the first octant (+x,+y,+z).
    /// Note that we are returning the packed form of UnitVec regardless
    /// of our stride here.
    UnitVec<P,1> abs() const {return UnitVec<P,1>( asVec3().abs(), true );}

    /// Return a new unit vector perpendicular to this one but otherwise
    /// arbitrary. Some care is taken to ensure good numerical conditioning
    /// for the result regardless of what goes in. Cost is about 50 flops.
    inline UnitVec<P,1> perp() const;

    /// This constructor is only for our friends whom we trust to
    /// give us an already-normalized vector which we simply accept as
    /// normalized without checking.
    UnitVec(const BaseVec& v, bool) : BaseVec(v) {}
    /// This constructor is only for our friends whom we trust to
    /// give us an already-normalized vector which we simply accept as
    /// normalized without checking (this version accepts an input
    /// vector of any stride).
    template <int S2> UnitVec(const Vec<3,RealP,S2>& v, bool) : BaseVec(v) { }
};


template <class P, int S> inline UnitVec<P,1> 
UnitVec<P,S>::perp() const {
    // Choose the coordinate axis which makes the largest angle
    // with this vector, that is, has the "least u" along it.
    const UnitVec<P,1> u(abs());    // reflect to first octant
    const int minAxis = u[0] <= u[1] ? (u[0] <= u[2] ? 0 : 2)
                                     : (u[1] <= u[2] ? 1 : 2);
    // Cross returns a Vec3 result which is then normalized.
    return UnitVec<P,1>( *this % UnitVec<P,1>(minAxis) );
}


//-----------------------------------------------------------------------------
/**
 * This type is used for the transpose of UnitVec, and as the returned row
 * type of a Rotation. Don't construct these directly.
 */
//-----------------------------------------------------------------------------
template <class P, int S>
class UnitRow : public Row<3,P,S> {
    typedef P   RealP;
public:
    typedef Row<3,P,S>      BaseRow;
    typedef UnitVec<P,S>    TransposeType;

    UnitRow() : BaseRow(NTraits<P>::getNaN()) { }

    /// Copy constructor does not require normalization.
    UnitRow(const UnitRow& u) 
    :   BaseRow(static_cast<const BaseRow&>(u)) {}

    /// Implicit conversion from UnitRow with different stride; no
    /// normalization required.
    template <int S2> UnitRow(const UnitRow<P,S2>& u)
    :   BaseRow(static_cast<const typename UnitRow<P,S2>::BaseRow&>(u)) { }

    /// Copy assignment does not require normalization.
    UnitRow& operator=(const UnitRow& u) 
    {   BaseRow::operator=(static_cast<const BaseRow&>(u)); 
        return *this; }

    /// Copy assignment from UnitRow with different stride; no computation needed.
    template <int S2> UnitRow& operator=(const UnitRow<P,S2>& u) 
    {   BaseRow::operator=(static_cast<const typename UnitRow<P,S2>::BaseRow&>(u));
        return *this; }

    /// Explicit conversion from Row to UnitRow, requiring expensive normalization.
    explicit UnitRow(const BaseRow& v) : BaseRow(v/v.norm()) {}
    /// Explicit conversion from Row of any stride to UnitRow, requiring expensive 
    /// normalization.
    template <int S2> 
    explicit UnitRow(const Row<3,P,S2>& v) : BaseRow(v/v.norm()) {}

    /// Create a unit row from explicitly specified measure numbers (x,y,z); 
    /// requires expensive normalization.
    UnitRow(const RealP& x, const RealP& y, const RealP& z)
    :   BaseRow(x,y,z)
    {   static_cast<BaseRow&>(*this) /= BaseRow::norm(); }

    /// Create a unit axis vector 100 010 001 given 0, 1, or 2.
    explicit UnitRow(int axis) : BaseRow(0) 
    {   assert(0 <= axis && axis <= 2);
        BaseRow::operator[](axis) = 1; }

    /// Return a const reference to the Row3 underlying this UnitRow.
    const BaseRow& asRow3() const  {return static_cast<const BaseRow&>(*this);}

    // Override Row3 methods which preserve length. These return the 
    // packed UnitRow regardless of our stride.

    /// Returns a new unit vector pointing in the opposite direction from this one;
    /// does \e not modify this UnitVec object. Cost is 3 flops.
    UnitRow<P,1> negate()    const  { return UnitRow<P,1>(-asRow3(),true); }
    /// Returns a new unit vector pointing in the opposite direction from this one.
    /// Cost is 3 flops.
    UnitRow<P,1> operator-() const  { return negate();}

    /// Return a const reference to this UnitRow reinterpreted as a UnitVec; no
    /// computation requires since this is just a type cast.
    const TransposeType&  operator~() const {return *reinterpret_cast<const TransposeType*>(this);}
    /// Return a writable reference to this UnitRow reinterpreted as a UnitVec; no
    /// computation requires since this is just a type cast.
    TransposeType& operator~() {return *reinterpret_cast<TransposeType*>(this);}

    // We have to define these here so that the non-const ones won't be
    // inherited. We don't trust anyone to write on one element of a UnitRow!

    /// Return one element of this unit row as a const reference; there is no 
    /// corresponding writable index function since changing a single element of
    /// a unit vector would violate the contract that it has unit length at all times.
    const RealP&  operator[](int i) const  { return BaseRow::operator[](i); }
    /// Return one element of this unit row as a const reference; there is no 
    /// corresponding writable index function since changing a single element of
    /// a unit vector would violate the contract that it has unit length at all times.
    const RealP&  operator()(int i) const  { return BaseRow::operator()(i); }

    /// Return a new UnitRow whose measure numbers are the absolute values
    /// of the ones here. This will still have unit length but will be
    /// a reflection of this unit vector into the first octant (+x,+y,+z).
    /// Note that we are returning the packed form of UnitRow regardless
    /// of our stride here.
    UnitRow<P,1> abs() const {return UnitRow<P,1>(asRow3().abs(),true);}

    /// Return a new UnitRow perpendicular to this one but otherwise
    /// arbitrary. Some care is taken to ensure good numerical conditioning
    /// for the result regardless of what goes in. Cost is about 50 flops.
    inline UnitRow<P,1> perp() const;

    // This constructor is only for our friends whom we trust to
    // give us an already-normalized vector.
    UnitRow( const BaseRow& v, bool ) : BaseRow(v) { }
    template <int S2> UnitRow( const Row<3,P,S2>& v, bool ) : BaseRow(v) { }
};

template <class P, int S>
inline UnitRow<P,1> UnitRow<P,S>::perp() const {
    // Choose the coordinate axis which makes the largest angle
    // with this vector, that is, has the "least u" along it.
    const UnitRow<P,1> u(abs());    // reflect to first octant
    const int minAxis = u[0] <= u[1] ? (u[0] <= u[2] ? 0 : 2)
                                     : (u[1] <= u[2] ? 1 : 2);
    // Cross returns a Row3 result which is then normalized.
    return UnitRow<P,1>(*this % UnitRow<P,1>(minAxis));
}



//------------------------------------------------------------------------------
}  // End of namespace SimTK

//--------------------------------------------------------------------------
#endif // SimTK_UNITVEC_H_
//--------------------------------------------------------------------------


