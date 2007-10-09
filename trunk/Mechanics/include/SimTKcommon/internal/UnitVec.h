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

//-----------------------------------------------------------------------------
#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/Constants.h"
#include "SimTKcommon/internal/Scalar.h"
#include "SimTKcommon/internal/SmallMatrix.h"
//-----------------------------------------------------------------------------
#include <iosfwd>  // Forward declaration of iostream
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
namespace SimTK {

//-----------------------------------------------------------------------------
// Forward declarations
template <int S> class UnitRow;


//-----------------------------------------------------------------------------
/**
 * This class is a Vec3 plus an ironclad guarantee either that:
 *      - the length is one (to within a very small tolerance), or
 *      - all components are NaN.
 */
//-----------------------------------------------------------------------------
template <int S>
class UnitVec : public Vec<3,Real,S> {
public:
    typedef Vec<3,Real,S> BaseVec;
    typedef UnitRow<S>    TransposeType;

    UnitVec() : BaseVec(NaN) { }

    // Copy constructor.
    UnitVec(const UnitVec& u) 
      : BaseVec(static_cast<const BaseVec&>(u)) { }

    // Automatic conversion from UnitVec with different stride.
    template <int S2> UnitVec(const UnitVec<S2>& u)
        : BaseVec(static_cast<const typename UnitVec<S2>::BaseVec&>(u)) { }

    // Explicit conversion from Vec to UnitVec, requiring expensive normalization.
    explicit UnitVec(const BaseVec& v) : BaseVec(v/v.norm()) { }
    template <int S2> explicit UnitVec(const Vec<3,Real,S2>& v)
        : BaseVec(v/v.norm()) { }

    UnitVec(const Real& x, const Real& y, const Real& z) : BaseVec(x,y,z) {
        static_cast<BaseVec&>(*this) /= BaseVec::norm();
    }

    // Create a unit axis vector 100 010 001
    explicit UnitVec(int axis) : BaseVec(0) {
        assert(0 <= axis && axis <= 2);
        BaseVec::operator[](axis) = 1.;
    }

    UnitVec& operator=(const UnitVec& u) {
        BaseVec::operator=(static_cast<const BaseVec&>(u)); 
        return *this;
    }
    template <int S2> UnitVec& operator=(const UnitVec<S2>& u) {
        BaseVec::operator=(static_cast<const typename UnitVec<S2>::BaseVec&>(u));
        return *this;
    }

    const BaseVec& asVec3() const {return static_cast<const BaseVec&>(*this);}

    // Override Vec3 methods which preserve length. These return the 
    // packed UnitVec regardless of our stride.
    UnitVec<1> negate()    const {return UnitVec<1>(-asVec3(),true);}
    UnitVec<1> operator-() const {return negate();}

    const TransposeType&  operator~() const  { return *reinterpret_cast<const TransposeType*>(this); }
    TransposeType&        operator~()        { return *reinterpret_cast<TransposeType*>(this); }

    // We have to define these here so that the non-const ones won't be
    // inherited. We don't trust anyone to write on one element of a UnitVec!
    const Real&  operator[](int i) const  { return BaseVec::operator[](i); }
    const Real&  operator()(int i) const  { return BaseVec::operator()(i); }

    // Return a vector whose measure numbers are the absolute values
    // of the ones here. This will still have unit length but will be
    // a reflection of this unit vector into the first octant (+x,+y,+z).
    // Note that we are returning the packed form of UnitVec regardless
    // of our stride here.
    UnitVec<1>  abs() const  { return UnitVec<1>( asVec3().abs(), true ); }

    // Return a unit vector perpendicular to this one (arbitrary).
    inline UnitVec<1>  perp() const;

    // This constructor is only for our friends whom we trust to
    // give us an already-normalized vector.
    UnitVec( const BaseVec& v, bool ) : BaseVec(v) { }
    template <int S2>  UnitVec( const Vec<3,Real,S2>& v, bool ) : BaseVec(v) { }
};


template <int S>
inline UnitVec<1> UnitVec<S>::perp() const {
    // Choose the coordinate axis which makes the largest angle
    // with this vector, that is, has the "least u" along it.
    const UnitVec<1> u(abs());    // reflect to first octant
    const int minAxis = u[0] <= u[1] ? (u[0] <= u[2] ? 0 : 2)
                                     : (u[1] <= u[2] ? 1 : 2);
    // Cross returns a Vec3 result which is then normalized.
    return UnitVec<1>( *this % UnitVec<1>(minAxis) );
}


//-----------------------------------------------------------------------------
/**
 * This type is used for the transpose of UnitVec, and as the returned row
 * type of a Rotation. Don't construct these directly.
 */
//-----------------------------------------------------------------------------
template <int S>
class UnitRow : public Row<3,Real,S> {
public:
    typedef Row<3,Real,S> BaseRow;
    typedef UnitVec<S>    TransposeType;

    UnitRow() : BaseRow(NaN) { }

    // Copy constructor.
    UnitRow(const UnitRow& u) 
      : BaseRow(static_cast<const BaseRow&>(u)) { }

    // Automatic conversion from UnitRow with different stride.
    template <int S2> UnitRow(const UnitRow<S2>& u)
        : BaseRow(static_cast<const typename UnitRow<S2>::BaseRow&>(u)) { }

    // Copy assignment.
    UnitRow& operator=(const UnitRow& u) {
        BaseRow::operator=(static_cast<const BaseRow&>(u)); 
        return *this;
    }
    // Assignment from UnitRow with different stride.
    template <int S2> UnitRow& operator=(const UnitRow<S2>& u) {
        BaseRow::operator=(static_cast<const typename UnitRow<S2>::BaseRow&>(u));
        return *this;
    }

    // Explicit conversion from Row to UnitRow, requiring expensive normalization.
    explicit UnitRow(const BaseRow& v) : BaseRow(v/v.norm()) { }
    template <int S2> explicit UnitRow(const Row<3,Real,S2>& v)
        : BaseRow(v/v.norm()) { }

    UnitRow(const Real& x, const Real& y, const Real& z) : BaseRow(x,y,z) {
        static_cast<BaseRow&>(*this) /= BaseRow::norm();
    }

    // Create a unit axis vector 100 010 001
    explicit UnitRow(int axis) : BaseRow(0) {
        assert(0 <= axis && axis <= 2);
        BaseRow::operator[](axis) = 1.;
    }

    const BaseRow&  asRow3() const  { return static_cast<const BaseRow&>(*this); }

    // Override Row3 methods which preserve length. These return the 
    // packed UnitRow regardless of our stride.
    UnitRow<1>  negate()    const  { return UnitRow<1>(-asRow3(),true); }
    UnitRow<1>  operator-() const  { return negate();}

    const TransposeType&  operator~() const { return *reinterpret_cast<const TransposeType*>(this); }
    TransposeType&        operator~()       { return *reinterpret_cast<TransposeType*>(this); }

    // We have to define these here so that the non-const ones won't be
    // inherited. We don't trust anyone to write on one element of a UnitRow!
    const Real&  operator[](int i) const  { return BaseRow::operator[](i); }
    const Real&  operator()(int i) const  { return BaseRow::operator()(i); }

    // Return a vector whose measure numbers are the absolute values
    // of the ones here. This will still have unit length but will be
    // a reflection of this unit vector into the first octant (+x,+y,+z).
    // Note that we are returning the packed form of UnitVec regardless
    // of our stride here.
    UnitRow<1>  abs() const  { return UnitRow<1>(asRow3().abs(),true); }

    // Return a unit row vector perpendicular to this one (arbitrary).
    inline UnitRow<1>  perp() const;

    // This constructor is only for our friends whom we trust to
    // give us an already-normalized vector.
    UnitRow( const BaseRow& v, bool ) : BaseRow(v) { }
    template <int S2> UnitRow( const Row<3,Real,S2>& v, bool ) : BaseRow(v) { }
};

template <int S>
inline UnitRow<1> UnitRow<S>::perp() const {
    // Choose the coordinate axis which makes the largest angle
    // with this vector, that is, has the "least u" along it.
    const UnitRow<1> u(abs());    // reflect to first octant
    const int minAxis = u[0] <= u[1] ? (u[0] <= u[2] ? 0 : 2)
                                     : (u[1] <= u[2] ? 1 : 2);
    // Cross returns a Row3 result which is then normalized.
    return UnitRow<1>(*this % UnitRow<1>(minAxis));
}


// UnitVec3 is more intelligible name for UnitVec<1> now that UnitVec class is defined
typedef UnitVec<1> UnitVec3;


//------------------------------------------------------------------------------
}  // End of namespace SimTK

//--------------------------------------------------------------------------
#endif // SimTK_UNITVEC_H_
//--------------------------------------------------------------------------


