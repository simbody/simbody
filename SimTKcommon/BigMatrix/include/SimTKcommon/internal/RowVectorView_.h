#ifndef SimTK_SIMMATRIX_ROWVECTORVIEW_H_
#define SimTK_SIMMATRIX_ROWVECTORVIEW_H_

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
Define the SimTK::RowVectorView_ class that is part of Simbody's BigMatrix
toolset. **/

namespace SimTK {

//==============================================================================
//                              ROW VECTOR VIEW
//==============================================================================
/** @brief (Advanced) This class is identical to RowVector_ except that it has
shallow (reference) copy and assignment semantics.

Despite the name, this may be an owner if a RowVector_ is recast to a
%RowVectorView_. However, there are no owner constructors for %RowVectorView_.
@see RowVector_, RowVectorBase, MatrixView_ **/
template <class ELT> class RowVectorView_ : public RowVectorBase<ELT> {
    typedef RowVectorBase<ELT>                              Base;
    typedef typename CNT<ELT>::Scalar                       S;
    typedef typename CNT<ELT>::Number                       Number;
    typedef typename CNT<ELT>::StdNumber                    StdNumber;
    typedef RowVectorView_<ELT>                             T;
    typedef RowVectorView_< typename CNT<ELT>::TNeg >       TNeg;
    typedef VectorView_< typename CNT<ELT>::THerm >         THerm;
public:
    // Default construction is suppressed.
    // Uses default destructor.

    // Create a RowVectorView_ handle using a given helper rep.
    explicit RowVectorView_(MatrixHelperRep<S>* hrep) : Base(hrep) {}

    // Copy constructor is shallow. CAUTION: despite const argument, this preserves writability
    // if it was present in the source. This is necessary to allow temporary views to be
    // created and used as lvalues.
    RowVectorView_(const RowVectorView_& r)
      : Base(const_cast<MatrixHelper<S>&>(r.getHelper()), typename MatrixHelper<S>::ShallowCopy()) { }

    // Copy assignment is deep but not reallocating.
    RowVectorView_& operator=(const RowVectorView_& r) {
        Base::operator=(r); return *this;
    }

    // Ask for shallow copy
    explicit RowVectorView_(const MatrixHelper<S>& h) : Base(h, typename MatrixHelper<S>::ShallowCopy()) { }
    explicit RowVectorView_(MatrixHelper<S>&       h) : Base(h, typename MatrixHelper<S>::ShallowCopy()) { }

    RowVectorView_& operator=(const Base& b) { Base::operator=(b); return *this; }

    RowVectorView_& operator=(const ELT& v) { Base::operator=(v); return *this; }

    template <class EE> RowVectorView_& operator=(const RowVectorBase<EE>& m)
      { Base::operator=(m); return*this; }
    template <class EE> RowVectorView_& operator+=(const RowVectorBase<EE>& m)
      { Base::operator+=(m); return*this; }
    template <class EE> RowVectorView_& operator-=(const RowVectorBase<EE>& m)
      { Base::operator-=(m); return*this; }

    RowVectorView_& operator*=(const StdNumber& t) { Base::operator*=(t); return *this; }
    RowVectorView_& operator/=(const StdNumber& t) { Base::operator/=(t); return *this; }
    RowVectorView_& operator+=(const ELT& b) { this->elementwiseAddScalarInPlace(b); return *this; }
    RowVectorView_& operator-=(const ELT& b) { this->elementwiseSubtractScalarInPlace(b); return *this; }

private:
    // NO DATA MEMBERS ALLOWED
    RowVectorView_(); // default construction suppressed (what is it a view of?)
};

} //namespace SimTK

#endif // SimTK_SIMMATRIX_ROWVECTORVIEW_H_
