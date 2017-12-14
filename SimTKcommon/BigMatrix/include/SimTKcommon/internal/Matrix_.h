#ifndef SimTK_SIMMATRIX_MATRIX_H_
#define SimTK_SIMMATRIX_MATRIX_H_

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
Define the SimTK::Matrix_ class that is part of Simbody's BigMatrix toolset. **/

namespace SimTK {

//==============================================================================
//                                MATRIX
//==============================================================================
/** @brief This is the matrix class intended to appear in user code for large,
variable size matrices.

@ingroup MatVecUtilities

More commonly, the typedef @ref SimTK::Matrix "Matrix" is used instead; that is
just an abbreviation for \c Matrix_<Real>.

A %Matrix_ can be a fixed-size view of someone else's data, or can be a 
resizable data owner itself.

@see Mat for handling of small, fixed-size matrices with no runtime overhead.
@see Vector_ for variable size, one-dimensional column vector.
@see MatrixView_, MatrixBase
**/
//------------------------------------------------------------------------------
template <class ELT> class Matrix_ : public MatrixBase<ELT> {
    typedef typename CNT<ELT>::Scalar       S;
    typedef typename CNT<ELT>::Number       Number;
    typedef typename CNT<ELT>::StdNumber    StdNumber;

    typedef typename CNT<ELT>::TNeg         ENeg;
    typedef typename CNT<ELT>::THerm        EHerm;

    typedef MatrixBase<ELT>     Base;
    typedef MatrixBase<ENeg>    BaseNeg;
    typedef MatrixBase<EHerm>   BaseHerm;

    typedef MatrixView_<ELT>    TView;
    typedef Matrix_<ENeg>       TNeg;

public:
    Matrix_() : Base() { }
    explicit Matrix_(const MatrixCommitment& mc) : Base(mc) {}

    // Copy constructor is deep.
    Matrix_(const Matrix_& src) : Base(src) { }

    // Assignment is a deep copy and will also allow reallocation if this Matrix
    // doesn't have a view.
    Matrix_& operator=(const Matrix_& src) { 
        Base::operator=(src); return *this;
    }

    // Force a deep copy of the view or whatever this is.
    // Note that this is an implicit conversion.
    Matrix_(const Base& v) : Base(v) {}   // e.g., MatrixView

    // Allow implicit conversion from a source matrix that
    // has a negated version of ELT.
    Matrix_(const BaseNeg& v) : Base(v) {}

    // TODO: implicit conversion from conjugate. This is trickier
    // since real elements are their own conjugate so you'll get
    // duplicate methods defined from Matrix_(BaseHerm) and Matrix_(Base).

    Matrix_(int m, int n) : Base(MatrixCommitment(), m, n) {}

    Matrix_(int m, int n, const ELT* cppInitialValuesByRow) 
    :   Base(MatrixCommitment(), m, n, cppInitialValuesByRow) {}
    Matrix_(int m, int n, const ELT& initialValue) 
    :   Base(MatrixCommitment(), m, n, initialValue) {}
    
    Matrix_(int m, int n, int leadingDim, const S* data) // read only
    :   Base(MatrixCommitment(), MatrixCharacter::LapackFull(m,n), 
             leadingDim, data) {}
    Matrix_(int m, int n, int leadingDim, S* data) // writable
    :   Base(MatrixCommitment(), MatrixCharacter::LapackFull(m,n), 
             leadingDim, data) {}
    
    /// Convert a Mat to a Matrix_.
    template <int M, int N, int CS, int RS>
    explicit Matrix_(const Mat<M,N,ELT,CS,RS>& mat)
    :   Base(MatrixCommitment(), M, N)
    {   for (int i = 0; i < M; ++i)
            for (int j = 0; j < N; ++j)
                this->updElt(i, j) = mat(i, j); }

    Matrix_& operator=(const ELT& v) { Base::operator=(v); return *this; }

    template <class EE> Matrix_& operator=(const MatrixBase<EE>& m)
      { Base::operator=(m); return*this; }
    template <class EE> Matrix_& operator+=(const MatrixBase<EE>& m)
      { Base::operator+=(m); return*this; }
    template <class EE> Matrix_& operator-=(const MatrixBase<EE>& m)
      { Base::operator-=(m); return*this; }

    Matrix_& operator*=(const StdNumber& t) { Base::operator*=(t); return *this; }
    Matrix_& operator/=(const StdNumber& t) { Base::operator/=(t); return *this; }
    Matrix_& operator+=(const ELT& r)       { this->updDiag() += r; return *this; }
    Matrix_& operator-=(const ELT& r)       { this->updDiag() -= r; return *this; }  

    const TNeg& negate()    const {return *reinterpret_cast<const TNeg*>(this); }
    TNeg&       updNegate()       {return *reinterpret_cast<TNeg*>(this); }

    const TNeg& operator-() const {return negate();}
    TNeg&       operator-()       {return updNegate();}
   
    // Functions to be used for Scripting in MATLAB and languages that do not support operator overloading
    /** toString() returns a string representation of the Matrix_. Please refer to operator<< for details. **/
    std::string toString() const {
        std::stringstream stream;
        stream <<  (*this) ;
        return stream.str(); 
    }
    /** Variant of indexing operator that's scripting friendly to get entry (i, j) **/
    const ELT& get(int i,int j) const { return this->getElt(i,j); }
    /** Variant of indexing operator that's scripting friendly to set entry (i, j) **/
    void       set(int i,int j, const ELT& value)       { this->updElt(i,j)=value; }

private:
    // NO DATA MEMBERS ALLOWED
};

} //namespace SimTK

#endif // SimTK_SIMMATRIX_MATRIX_H_
