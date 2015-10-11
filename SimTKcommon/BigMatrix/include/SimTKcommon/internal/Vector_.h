#ifndef SimTK_SIMMATRIX_VECTOR_H_
#define SimTK_SIMMATRIX_VECTOR_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-13 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
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

/** @file
Define the SimTK::Vector_ class that is part of Simbody's BigMatrix toolset. **/

namespace SimTK {

//==============================================================================
//                                 VECTOR
//==============================================================================
/** @brief This is the vector class intended to appear in user code for large,
variable size column vectors.

@ingroup MatVecUtilities

More commonly, the typedef @ref SimTK::Vector "Vector" is used instead; that is
just an abbreviation for \c Vector_<Real>.

A %Vector_ can be a fixed-size view of someone else's data, or can be a
resizable data owner itself, although of course it will always have just one
column.

@see Vec for handling of small, fixed-size vectors with no runtime overhead.
@see Matrix_ for variable %size, two-dimensional matrix.
**/
template <class ELT> class Vector_ : public VectorBase<ELT> {
    typedef typename CNT<ELT>::Scalar       S;
    typedef typename CNT<ELT>::Number       Number;
    typedef typename CNT<ELT>::StdNumber    StdNumber;
    typedef typename CNT<ELT>::TNeg         ENeg;
    typedef VectorBase<ELT>                 Base;
    typedef VectorBase<ENeg>                BaseNeg;
public:
    // Uses default destructor.

    /** @name                  Owner constructors
    These constructors create a Vector_ object that allocates new data on the
    heap and serves as the owner of that data. Any existing data that is
    supplied in these constructors is \e copied, not \e shared. **/
    /**@{**/

    /** Default constructor creates a 0x1, reallocatable vector. **/
    Vector_() : Base() {}

    /** Copy constructor is deep, that is the source %Vector_ is copied, not
    referenced. **/
    Vector_(const Vector_& src) : Base(src) {}

    /** This copy constructor serves as an implicit conversion from objects of
    the base class type (for example, a VectorView_<ELT>), to objects of this
    Vector_<ELT> type. Note that the source object is copied, not
    referenced. **/
    Vector_(const Base& src) : Base(src) {}    // e.g., VectorView
    /** This copy constructor serves as an implicit conversion from objects of
    the base class but with negated elements, to objects of this Vector_<ELT>
    type.  Note that the source object is copied, not referenced. **/
    Vector_(const BaseNeg& src) : Base(src) {}

    /** Construct a %Vector_ with a given preallocated size, but with
    uninitialized values. In Debug builds the elements will be initialized to
    NaN as a debugging aid, but in Release (optimized) builds they will be
    uninitialized garbage so that we don't waste time setting them. **/
    explicit Vector_(int m) : Base(m) { }
    /** Construct an owner %Vector_ of a given size \a m and initialize it from
    a C++ array of \a m ELT values pointed to by \a cppInitialValues. Note that
    there is no way to check that the correct number of elements has been
    provided; make sure you have supplied enough of them. **/
    Vector_(int m, const ELT* cppInitialValues) : Base(m, cppInitialValues) {}
    /** Construct an owner %Vector_ of a given size \a m and initialize all the
    elements to the given ELT value \a initialValue. **/
    Vector_(int m, const ELT& initialValue) : Base(m, initialValue) {}

    /** This constructor creates a new owner %Vector_ that is the same length
    and has a copy of the contents of a given fixed-size Vec. **/
    template <int M>
    explicit Vector_(const Vec<M,ELT>& v) : Base(M) {
        for (int i = 0; i < M; ++i)
            this->updElt(i, 0) = v(i);
    }
    /**@}**/

    /** @name                  View constructors
    These constructors create a %Vector_ object that is a view of existing
    data that is owned by some other object. Any existing data that is
    supplied in these constructors is \e shared, not \e copied, so writing
    to the created %Vector_ view modifies the original data. **/
    /**@{**/

    /** Construct a %Vector_ which shares read-only, borrowed space with assumed
    element-to-element stride equal to the C++ element spacing. Last parameter
    is a required dummy to avoid overload conflicts when ELT=S; pass it as
    "true". **/
    Vector_(int m, const S* cppData, bool)
        : Base(m, Base::CppNScalarsPerElement, cppData) {}
    /** Construct a %Vector_ which shares writable, borrowed space with assumed
    element-to-element stride equal to the C++ element spacing. Last parameter
    is a required dummy to avoid overload conflicts when ELT=S; pass it as
    "true". **/
    Vector_(int m, S* cppData, bool)
        : Base(m, Base::CppNScalarsPerElement, cppData) {}

    /** Borrowed-space read-only construction with explicit stride supplied as
    "number of scalars between elements". Last parameter is a dummy to avoid
    overload conflicts; pass it as "true". **/
    Vector_(int m, int stride, const S* data, bool)
        : Base(m, stride, data) {}
    /** Borrowed-space writable construction with explicit stride supplied as
    "number of scalars between elements". Last parameter is a dummy to avoid
    overload conflicts; pass it as "true". **/
    Vector_(int m, int stride, S* data, bool)
        : Base(m, stride, data) {}
    /**@}**/

    /** @name              Assignment operators **/
    /**@{**/
    /** Copy assignment is deep and can be reallocating if this %Vector_ is an
    owner. Otherwise the source is copied into the destination view. **/
    Vector_& operator=(const Vector_& src)
    {   Base::operator=(src); return*this; }

    /** Like copy assignment but the source can be any object of the base
    type, even with a different element type \a EE as long as that element
    type is assignment-compatible to an element of our type \a ELT. **/
    template <class EE> Vector_& operator=(const VectorBase<EE>& src)
    {   Base::operator=(src); return*this; }

    /** Set all elements of this %Vector_ to the same value \a v. **/
    Vector_& operator=(const ELT& v) { Base::operator=(v); return *this; }

    /** Add in a conforming vector of the base type even if it has a different
    element type \a EE, provided that element type is add-compatible with
    our element type \a ELT, meaning that ELT+=EE is allowed. **/
    template <class EE> Vector_& operator+=(const VectorBase<EE>& m)
    {   Base::operator+=(m); return*this; }
    /** In-place subtract of a conforming vector of the base type even if it has
    a different element type \a EE, provided that element type is
    subtract-compatible with our element type \a ELT, meaning that ELT-=EE is
    allowed. **/
    template <class EE> Vector_& operator-=(const VectorBase<EE>& m)
    {   Base::operator-=(m); return*this; }

    /** In-place multiply of each element of this %Vector_ by a scalar \a t.
    Returns a reference to the now-modified %Vector_. **/
    Vector_& operator*=(const StdNumber& t) { Base::operator*=(t); return *this; }
    /** In-place divide of each element of this %Vector_ by a scalar \a t.
    Returns a reference to the now-modified %Vector_. **/
    Vector_& operator/=(const StdNumber& t) { Base::operator/=(t); return *this; }

    /** In-place add to each element of this %Vector_ the same given value \a b.
    Returns a reference to the now-modified %Vector_. This is the same
    operation as elementwiseAddScalarInPlace(). **/
    Vector_& operator+=(const ELT& b)
    {   this->elementwiseAddScalarInPlace(b); return *this; }
    /** In-place subtract from each element of this %Vector_ the same given
    value \a b. Returns a reference to the now-modified %Vector_. This is the
    same operation as elementwiseSubtractScalarInPlace(). **/
    Vector_& operator-=(const ELT& b)
    {   this->elementwiseSubtractScalarInPlace(b); return *this; }
    /**@}**/

    /** @name             Operator equivalents for scripting
    Functions to be used for Scripting in MATLAB and languages that do not
    support operator overloading. **/
    /**@{**/
    /** toString() returns a string representation of the %Vector_. Please refer
    to operator<< for details. **/
    std::string toString() const {
        std::stringstream stream;
        stream << (*this) ;
        return stream.str();
    }
    /** Variant of operator[] that's scripting friendly to get ith element. **/
    const ELT& get(int i) const { return (*this)[i]; }
    /** Variant of operator[] that's scripting friendly to set ith element. **/
    void set(int i, const ELT& value)  { (*this)[i]=value; }
    /**@}**/

private:
    // NO DATA MEMBERS ALLOWED
};

} //namespace SimTK

#endif // SimTK_SIMMATRIX_VECTOR_H_
