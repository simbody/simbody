#ifndef SimTK_SIMMATRIX_MATRIXVIEW_H_
#define SimTK_SIMMATRIX_MATRIXVIEW_H_

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
Define the SimTK::MatrixView_ class that is part of Simbody's BigMatrix
toolset. **/

namespace SimTK {

//==============================================================================
//                               MATRIX VIEW
//==============================================================================
/** @brief (Advanced) This class is identical to Matrix_ except that it has
shallow (reference) copy and assignment semantics.

Despite the name, this may be an owner if a Matrix_ is recast to a MatrixView_.
However, there are no owner constructors for MatrixView_.

@see Matrix_, MatrixBase, VectorView_ **/
template <class ELT> class MatrixView_ : public MatrixBase<ELT> {
    typedef MatrixBase<ELT>                 Base;
    typedef typename CNT<ELT>::Scalar       S;
    typedef typename CNT<ELT>::StdNumber    StdNumber;
public:
    // Default construction is suppressed.
    // Uses default destructor.

    // Create a MatrixView_ handle using a given helper rep.
    explicit MatrixView_(MatrixHelperRep<S>* hrep) : Base(hrep) {}

    // Copy constructor is shallow. CAUTION: despite const argument, this preserves writability
    // if it was present in the source. This is necessary to allow temporary views to be
    // created and used as lvalues.
    MatrixView_(const MatrixView_& m)
      : Base(MatrixCommitment(),
             const_cast<MatrixHelper<S>&>(m.getHelper()),
             typename MatrixHelper<S>::ShallowCopy()) {}

    // Copy assignment is deep but not reallocating.
    MatrixView_& operator=(const MatrixView_& m) {
        Base::operator=(m); return *this;
    }

    // Ask for shallow copy
    MatrixView_(const MatrixHelper<S>& h) : Base(MatrixCommitment(), h, typename MatrixHelper<S>::ShallowCopy()) { }
    MatrixView_(MatrixHelper<S>&       h) : Base(MatrixCommitment(), h, typename MatrixHelper<S>::ShallowCopy()) { }

    MatrixView_& operator=(const Matrix_<ELT>& v)     { Base::operator=(v); return *this; }
    MatrixView_& operator=(const ELT& e)              { Base::operator=(e); return *this; }

    template <class EE> MatrixView_& operator=(const MatrixBase<EE>& m)
      { Base::operator=(m); return *this; }
    template <class EE> MatrixView_& operator+=(const MatrixBase<EE>& m)
      { Base::operator+=(m); return *this; }
    template <class EE> MatrixView_& operator-=(const MatrixBase<EE>& m)
      { Base::operator-=(m); return *this; }

    MatrixView_& operator*=(const StdNumber& t) { Base::operator*=(t); return *this; }
    MatrixView_& operator/=(const StdNumber& t) { Base::operator/=(t); return *this; }
    MatrixView_& operator+=(const ELT& r)       { this->updDiag() += r; return *this; }
    MatrixView_& operator-=(const ELT& r)       { this->updDiag() -= r; return *this; }

    operator const Matrix_<ELT>&() const { return *reinterpret_cast<const Matrix_<ELT>*>(this); }
    operator Matrix_<ELT>&()             { return *reinterpret_cast<Matrix_<ELT>*>(this); }

private:
    // NO DATA MEMBERS ALLOWED
    MatrixView_(); // default constructor suppressed (what's it a view of?)
};



} //namespace SimTK

#endif // SimTK_SIMMATRIX_MATRIXVIEW_H_
