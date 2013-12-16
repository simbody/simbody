#ifndef SimTK_SIMMATRIX_VECTORVIEW_H_
#define SimTK_SIMMATRIX_VECTORVIEW_H_

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
Define the SimTK::VectorView_ class that is part of Simbody's BigMatrix 
toolset. **/

namespace SimTK {

//==============================================================================
//                                VECTOR VIEW
//==============================================================================
/** @brief (Advanced) This class is identical to Vector_ except that it has 
shallow (reference) copy and assignment semantics. 

Despite the name, this may be an owner if a Vector_ is recast to a %VectorView_.
However, there are no owner constructors for %VectorView_. 
@see Vector_, VectorBase, MatrixView_ **/
template <class ELT> class VectorView_ : public VectorBase<ELT> {
    typedef VectorBase<ELT>                             Base;
    typedef typename CNT<ELT>::Scalar                   S;
    typedef typename CNT<ELT>::Number                   Number;
    typedef typename CNT<ELT>::StdNumber                StdNumber;
    typedef VectorView_<ELT>                            T;
    typedef VectorView_< typename CNT<ELT>::TNeg >      TNeg;
    typedef RowVectorView_< typename CNT<ELT>::THerm >  THerm;
public:
    // Default construction is suppressed.
    // Uses default destructor.

    // Create a VectorView_ handle using a given helper rep. 
    explicit VectorView_(MatrixHelperRep<S>* hrep) : Base(hrep) {}

    // Copy constructor is shallow. CAUTION: despite const argument, this preserves writability
    // if it was present in the source. This is necessary to allow temporary views to be
    // created and used as lvalues.
    VectorView_(const VectorView_& v) 
      : Base(const_cast<MatrixHelper<S>&>(v.getHelper()), typename MatrixHelper<S>::ShallowCopy()) { }

    // Copy assignment is deep but not reallocating.
    VectorView_& operator=(const VectorView_& v) {
        Base::operator=(v); return *this;
    }

    // Ask for shallow copy    
    explicit VectorView_(const MatrixHelper<S>& h) : Base(h, typename MatrixHelper<S>::ShallowCopy()) { }
    explicit VectorView_(MatrixHelper<S>&       h) : Base(h, typename MatrixHelper<S>::ShallowCopy()) { }
    
    VectorView_& operator=(const Base& b) { Base::operator=(b); return *this; }

    VectorView_& operator=(const ELT& v) { Base::operator=(v); return *this; } 

    template <class EE> VectorView_& operator=(const VectorBase<EE>& m)
      { Base::operator=(m); return*this; }
    template <class EE> VectorView_& operator+=(const VectorBase<EE>& m)
      { Base::operator+=(m); return*this; }
    template <class EE> VectorView_& operator-=(const VectorBase<EE>& m)
      { Base::operator-=(m); return*this; }

    VectorView_& operator*=(const StdNumber& t) { Base::operator*=(t); return *this; }
    VectorView_& operator/=(const StdNumber& t) { Base::operator/=(t); return *this; }
    VectorView_& operator+=(const ELT& b) { this->elementwiseAddScalarInPlace(b); return *this; }
    VectorView_& operator-=(const ELT& b) { this->elementwiseSubtractScalarInPlace(b); return *this; }

private:
    // NO DATA MEMBERS ALLOWED
    VectorView_(); // default construction suppressed (what's it a View of?)
};

} //namespace SimTK

#endif // SimTK_SIMMATRIX_VECTORVIEW_H_
