#ifndef SimTK_SIMMATRIX_ROWVECTOR_H_
#define SimTK_SIMMATRIX_ROWVECTOR_H_

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
Define the SimTK::RowVector_ class that is part of Simbody's BigMatrix
toolset. **/

namespace SimTK {

//==============================================================================
//                                ROW VECTOR
//==============================================================================
/** @brief Represents a variable size row vector; much less common than the
column vector type Vector_.

@ingroup MatVecUtilities

Row vectors are much less commonly used than column vectors; they mostly arise
implicitly as the type of a transposed column vector (represented by Simbody's
Vector_ class). However, you want to use rows this is the class intended to
appear in user code. It can be a fixed-size view of someone else's data, or can
be a resizable data owner itself, although of course it will always have just
one row.

@see Row for handling of small, fixed-size row vectors with no runtime overhead
@see Matrix_ for variable %size, two-dimensional matrix.
@see RowVectorView_, RowVectorBase
**/
template <class ELT> class RowVector_ : public RowVectorBase<ELT> {
    typedef typename CNT<ELT>::Scalar       S;
    typedef typename CNT<ELT>::Number       Number;
    typedef typename CNT<ELT>::StdNumber    StdNumber;
    typedef typename CNT<ELT>::TNeg         ENeg;

    typedef RowVectorBase<ELT>              Base;
    typedef RowVectorBase<ENeg>             BaseNeg;
public:
    RowVector_() : Base() {}   // 1x0 reallocatable
    // Uses default destructor.

    // Copy constructor is deep.
    RowVector_(const RowVector_& src) : Base(src) {}

    // Implicit conversions.
    RowVector_(const Base& src) : Base(src) {}    // e.g., RowVectorView
    RowVector_(const BaseNeg& src) : Base(src) {}

    // Copy assignment is deep and can be reallocating if this RowVector
    // has no View.
    RowVector_& operator=(const RowVector_& src) {
        Base::operator=(src); return*this;
    }


    explicit RowVector_(int n) : Base(n) { }
    RowVector_(int n, const ELT* cppInitialValues) : Base(n, cppInitialValues) {}
    RowVector_(int n, const ELT& initialValue)     : Base(n, initialValue) {}

    /// Construct a Vector which uses borrowed space with assumed
    /// element-to-element stride equal to the C++ element spacing.
    /// Last parameter is a dummy to avoid overload conflicts when ELT=S;
    /// pass it as "true".
    RowVector_(int n, const S* cppData, bool): Base(n, Base::CppNScalarsPerElement, cppData) {}
    RowVector_(int n,       S* cppData, bool): Base(n, Base::CppNScalarsPerElement, cppData) {}

    /// Borrowed-space construction with explicit stride supplied as
    /// "number of scalars between elements". Last parameter is a
    /// dummy to avoid overload conflicts; pass it as "true".
    RowVector_(int n, int stride, const S* data, bool) : Base(n, stride, data) {}
    RowVector_(int n, int stride,       S* data, bool) : Base(n, stride, data) {}

    /// Convert a Row to a RowVector_.
    template <int M>
    explicit RowVector_(const Row<M,ELT>& v) : Base(M) {
        for (int i = 0; i < M; ++i)
            this->updElt(0, i) = v(i);
    }

    RowVector_& operator=(const ELT& v) { Base::operator=(v); return *this; }

    template <class EE> RowVector_& operator=(const RowVectorBase<EE>& b)
      { Base::operator=(b); return*this; }
    template <class EE> RowVector_& operator+=(const RowVectorBase<EE>& b)
      { Base::operator+=(b); return*this; }
    template <class EE> RowVector_& operator-=(const RowVectorBase<EE>& b)
      { Base::operator-=(b); return*this; }

    RowVector_& operator*=(const StdNumber& t) { Base::operator*=(t); return *this; }
    RowVector_& operator/=(const StdNumber& t) { Base::operator/=(t); return *this; }
    RowVector_& operator+=(const ELT& b) { this->elementwiseAddScalarInPlace(b); return *this; }
    RowVector_& operator-=(const ELT& b) { this->elementwiseSubtractScalarInPlace(b); return *this; }

private:
    // NO DATA MEMBERS ALLOWED
};

} //namespace SimTK

#endif // SimTK_SIMMATRIX_ROWVECTOR_H_
