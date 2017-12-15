#ifndef SimTK_SIMMATRIX_MATRIXBASE_H_
#define SimTK_SIMMATRIX_MATRIXBASE_H_

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
Define the SimTK::MatrixBase class that is part of Simbody's BigMatrix 
toolset. **/

namespace SimTK {

//==============================================================================
//                               MATRIX BASE
//==============================================================================
/** @brief This is the common base class for Simbody's Vector_ and Matrix_
classes for handling large, variable-sized vectors and matrices.

%MatrixBase does not normally appear in user programs. Instead classes
Vector_ and Matrix_ are used, or more commonly the typedefs Vector and Matrix
which are abbreviations for Vector_<Real> and Matrix_<Real>.

%Matrix base is a variable-size 2d matrix of Composite Numerical Type (CNT) 
elements. This is a container of such elements, it is NOT a Composite Numerical
Type itself. 

<h2>Implementation</h2>
MatrixBase<ELT> uses MatrixHelper<S> for implementation, where S 
is ELT::Scalar, that is, the underlying float, double,
complex<float>, negator<conjugate<double>>, 
etc. from which ELT is constructed. This is a finite set of which all
members are explicitly instantiated in the implementation code, so 
clients don't have to know how anything is implemented.

MatrixBase is the only class in the Matrix/Vector family which has any
data members (it has exactly one MatrixHelper, which itself consists only
of a single pointer to an opaque class). Thus all other objects
in this family (that is, derived from MatrixBase) are exactly the same
size in memory and may be "reinterpreted" as appropriate. For example,
a Vector may be reinterpreted as a Matrix or vice versa, provided runtime
requirements are met (e.g., exactly 1 column).

Unlike the small matrix classes, very little is encoded in the type.
Only the element type, and matrix vs. vector vs. row are in the type;
everything else like shape, storage layout, and writability are handled
at run time. **/
//  ----------------------------------------------------------------------------
template <class ELT> class MatrixBase {  
public:
    // These typedefs are handy, but despite appearances you cannot 
    // treat a MatrixBase as a composite numerical type. That is,
    // CNT<MatrixBase> will not compile, or if it does it won't be
    // meaningful.

    typedef ELT                                 E;
    typedef typename CNT<E>::TNeg               ENeg;
    typedef typename CNT<E>::TWithoutNegator    EWithoutNegator;
    typedef typename CNT<E>::TReal              EReal;
    typedef typename CNT<E>::TImag              EImag;
    typedef typename CNT<E>::TComplex           EComplex;
    typedef typename CNT<E>::THerm              EHerm;       
    typedef typename CNT<E>::TPosTrans          EPosTrans;

    typedef typename CNT<E>::TAbs               EAbs;
    typedef typename CNT<E>::TStandard          EStandard;
    typedef typename CNT<E>::TInvert            EInvert;
    typedef typename CNT<E>::TNormalize         ENormalize;
    typedef typename CNT<E>::TSqHermT           ESqHermT;
    typedef typename CNT<E>::TSqTHerm           ESqTHerm;

    typedef typename CNT<E>::Scalar             EScalar;
    typedef typename CNT<E>::Number             ENumber;
    typedef typename CNT<E>::StdNumber          EStdNumber;
    typedef typename CNT<E>::Precision          EPrecision;
    typedef typename CNT<E>::ScalarNormSq       EScalarNormSq;

    typedef EScalar    Scalar;        // the underlying Scalar type
    typedef ENumber    Number;        // negator removed from Scalar
    typedef EStdNumber StdNumber;     // conjugate goes to complex
    typedef EPrecision Precision;     // complex removed from StdNumber
    typedef EScalarNormSq  ScalarNormSq;      // type of scalar^2

    typedef MatrixBase<ENeg>             TNeg;
    typedef MatrixBase<EWithoutNegator>  TWithoutNegator;
    typedef MatrixBase<EReal>            TReal;
    typedef MatrixBase<EImag>            TImag;
    typedef MatrixBase<EComplex>         TComplex;
    typedef MatrixBase<EHerm>            THerm; 
    typedef MatrixBase<E>                TPosTrans;

    typedef MatrixBase<EAbs>             TAbs;
    typedef MatrixBase<EStandard>        TStandard;
    typedef MatrixBase<EInvert>          TInvert;
    typedef MatrixBase<ENormalize>       TNormalize;
    typedef MatrixBase<ESqHermT>         TSqHermT;  // ~Mat*Mat
    typedef MatrixBase<ESqTHerm>         TSqTHerm;  // Mat*~Mat

    const MatrixCommitment& getCharacterCommitment() const {return helper.getCharacterCommitment();}
    const MatrixCharacter& getMatrixCharacter()     const {return helper.getMatrixCharacter();}

    /// Change the handle commitment for this matrix handle; only allowed if the 
    /// handle is currently clear.
    void commitTo(const MatrixCommitment& mc)
    {   helper.commitTo(mc); }

    // This gives the resulting matrix type when (m(i,j) op P) is applied to each element.
    // It will have element types which are the regular composite result of E op P.
    template <class P> struct EltResult { 
        typedef MatrixBase<typename CNT<E>::template Result<P>::Mul> Mul;
        typedef MatrixBase<typename CNT<E>::template Result<P>::Dvd> Dvd;
        typedef MatrixBase<typename CNT<E>::template Result<P>::Add> Add;
        typedef MatrixBase<typename CNT<E>::template Result<P>::Sub> Sub;
    };

    /// Return the number of rows m in the logical shape of this matrix.
    int  nrow() const {return helper.nrow();}
    /// Return the number of columns n in the logical shape of this matrix.
    int  ncol() const {return helper.ncol();}

    /// Return the number of elements in the \e logical shape of this matrix.
    /// This has nothing to do with how many elements are actually stored;
    /// it is simply the product of the logical number of rows and columns,
    /// that is, nrow()*ncol(). Note that although each dimension is limited
    /// to a 32 bit size, the product of those dimensions may be > 32 bits 
    /// on a 64 bit machine so the return type may be larger than that of
    /// nrow() and ncol().
    ptrdiff_t nelt() const {return helper.nelt();}

    /// Return true if either dimension of this Matrix is resizable.
    bool isResizeable() const {return getCharacterCommitment().isResizeable();}

    enum { 
        NScalarsPerElement    = CNT<E>::NActualScalars,
        CppNScalarsPerElement = sizeof(E) / sizeof(Scalar)
    };
  
    /// The default constructor builds a 0x0 matrix managed by a helper that
    /// understands how many scalars there are in one of our elements but is
    /// otherwise uncommitted.
    MatrixBase() : helper(NScalarsPerElement,CppNScalarsPerElement) {}

    /// This constructor allocates the default matrix a completely uncommitted
    /// matrix commitment, given particular initial dimensions.
    MatrixBase(int m, int n) 
    :   helper(NScalarsPerElement,CppNScalarsPerElement,MatrixCommitment(),m,n) {}

    /// This constructor takes a handle commitment and allocates the default
    /// matrix for that kind of commitment. If a dimension is set to a 
    /// particular (unchangeable) value in the commitment then the initial
    /// allocation will use that value. Unlocked dimensions are given the
    /// smallest value consistent with other committed attributes, typically 0.
    explicit MatrixBase(const MatrixCommitment& commitment) 
    :   helper(NScalarsPerElement,CppNScalarsPerElement,commitment) {}


    /// This constructor takes a handle commitment and allocates the default
    /// matrix for that kind of commitment given particular initial minimum
    /// dimensions, which cannot be larger than those permitted by the 
    /// commitment.
    MatrixBase(const MatrixCommitment& commitment, int m, int n) 
    :   helper(NScalarsPerElement,CppNScalarsPerElement,commitment,m,n) {}

    /// Copy constructor is a deep copy (not appropriate for views!).    
    MatrixBase(const MatrixBase& b)
      : helper(b.helper.getCharacterCommitment(), 
               b.helper, typename MatrixHelper<Scalar>::DeepCopy()) { }

    /// Implicit conversion from matrix with negated elements (otherwise this
    /// is just like the copy constructor.
    MatrixBase(const TNeg& b)
      : helper(b.helper.getCharacterCommitment(),
               b.helper, typename MatrixHelper<Scalar>::DeepCopy()) { }
    
    /// Copy assignment is a deep copy but behavior depends on type of lhs: if 
    /// view, rhs must match. If owner, we reallocate and copy rhs.
    MatrixBase& copyAssign(const MatrixBase& b) {
        helper.copyAssign(b.helper);
        return *this;
    }
    MatrixBase& operator=(const MatrixBase& b) { return copyAssign(b); }


    /// View assignment is a shallow copy, meaning that we disconnect the MatrixBase 
    /// from whatever it used to refer to (destructing as necessary), then make it a new view
    /// for the data descriptor referenced by the source.
    /// CAUTION: we always take the source as const, but that is ignored in 
    /// determining whether the resulting view is writable. Instead, that is
    /// inherited from the writability status of the source. We have to do this
    /// in order to allow temporary view objects to be writable -- the compiler
    /// creates temporaries like m(i,j,m,n) as const.
    MatrixBase& viewAssign(const MatrixBase& src) {
        helper.writableViewAssign(const_cast<MatrixHelper<Scalar>&>(src.helper));
        return *this;
    }

    // default destructor

    /// Initializing constructor with all of the initially-allocated elements
    /// initialized to the same value. The given dimensions are treated as
    /// minimum dimensions in case the commitment requires more. So it is 
    /// always permissible to set them both to 0 in which case you'll get
    /// the smallest matrix that satisfies the commitment, with each of its
    /// elements (if any) set to the given initial value.
    MatrixBase(const MatrixCommitment& commitment, int m, int n, const ELT& initialValue) 
    :   helper(NScalarsPerElement, CppNScalarsPerElement, commitment, m, n)
    {   helper.fillWith(reinterpret_cast<const Scalar*>(&initialValue)); }  

    /// Initializing constructor with the initially-allocated elements
    /// initialized from a C++ array of elements, which is provided in
    /// <i>row major</i> order. The given dimensions are treated as
    /// minimum dimensions in case the commitment requires more. The
    /// array is presumed to be long enough to supply a value for each
    /// element. Note that C++ packing for elements may be different than
    /// Simmatrix packing of the same elements (Simmatrix packs them
    /// more tightly in some cases). So you should not use this constructor
    /// to copy elements from one Simmatrix matrix to another; this is
    /// exclusively for initializing a Simmatrix from a C++ array.
    MatrixBase(const MatrixCommitment& commitment, int m, int n, 
               const ELT* cppInitialValuesByRow) 
    :   helper(NScalarsPerElement, CppNScalarsPerElement, commitment, m, n)
    {   helper.copyInByRowsFromCpp(reinterpret_cast<const Scalar*>(cppInitialValuesByRow)); }
     
    /// @name           Matrix view of pre-exising data
    ///
    /// Non-resizeable view of someone else's already-allocated 
    /// memory of a size and storage type indicated by the supplied
    /// MatrixCharacter. The \a spacing argument has different interpretations
    /// depending on the storage format. Typically it is the leading
    /// dimension for Lapack-style full storage or stride for a vector.
    /// Spacing is in units like "number of scalars between elements" or
    /// "number of scalars between columns" so it can be used to deal
    /// with C++ packing vs. Simmatrix packing if necessary.
    /// @{

    /// Construct a read-only view of pre-existing data.
    MatrixBase(const MatrixCommitment& commitment, 
               const MatrixCharacter&  character, 
               int spacing, const Scalar* data) // read only data
    :   helper(NScalarsPerElement, CppNScalarsPerElement, 
               commitment, character, spacing, data) {}  

    /// Construct a writable view of pre-existing data.
    MatrixBase(const MatrixCommitment& commitment, 
               const MatrixCharacter&  character, 
               int spacing, Scalar* data) // writable data
    :   helper(NScalarsPerElement, CppNScalarsPerElement, 
               commitment, character, spacing, data) {}  
    /// @}
        
    // Create a new MatrixBase from an existing helper. Both shallow and deep copies are possible.
    MatrixBase(const MatrixCommitment& commitment, 
               MatrixHelper<Scalar>&   source, 
               const typename MatrixHelper<Scalar>::ShallowCopy& shallow) 
    :   helper(commitment, source, shallow) {}
    MatrixBase(const MatrixCommitment&      commitment, 
               const MatrixHelper<Scalar>&  source, 
               const typename MatrixHelper<Scalar>::ShallowCopy& shallow) 
    :   helper(commitment, source, shallow) {}
    MatrixBase(const MatrixCommitment&      commitment, 
               const MatrixHelper<Scalar>&  source, 
               const typename MatrixHelper<Scalar>::DeepCopy& deep)    
    :   helper(commitment, source, deep) {}

    /// This restores the MatrixBase to the state it would be in had it
    /// been constructed specifying only its handle commitment. The size will
    /// have been reduced to the smallest size consistent with the commitment.
    void clear() {helper.clear();}

    MatrixBase& operator*=(const StdNumber& t)  { helper.scaleBy(t);              return *this; }
    MatrixBase& operator/=(const StdNumber& t)  { helper.scaleBy(StdNumber(1)/t); return *this; }
    MatrixBase& operator+=(const MatrixBase& r) { helper.addIn(r.helper);         return *this; }
    MatrixBase& operator-=(const MatrixBase& r) { helper.subIn(r.helper);         return *this; }  

    template <class EE> MatrixBase(const MatrixBase<EE>& b)
      : helper(MatrixCommitment(),b.helper, typename MatrixHelper<Scalar>::DeepCopy()) { }

    template <class EE> MatrixBase& operator=(const MatrixBase<EE>& b) 
      { helper = b.helper; return *this; }
    template <class EE> MatrixBase& operator+=(const MatrixBase<EE>& b) 
      { helper.addIn(b.helper); return *this; }
    template <class EE> MatrixBase& operator-=(const MatrixBase<EE>& b) 
      { helper.subIn(b.helper); return *this; }

    /// Matrix assignment to an element sets only the *diagonal* elements to
    /// the indicated value; everything else is set to zero. This is particularly
    /// useful for setting a Matrix to zero or to the identity; for other values
    /// it creates a Matrix which acts like the scalar. That is, if the scalar
    /// is s and we do M=s, then multiplying another Matrix B by the resulting 
    /// diagonal matrix M gives the same result as multiplying B by s. That is
    /// (M=s)*B == s*B.
    ///
    /// NOTE: this must be overridden for Vector and RowVector since then scalar
    /// assignment is defined to copy the scalar to every element.
    MatrixBase& operator=(const ELT& t) { 
        setToZero(); updDiag().setTo(t); 
        return *this;
    }

    /// Set M's diagonal elements to a "scalar" value S, and all off-diagonal
    /// elements to zero. S can be any type which is assignable to an element
    /// of type E. This is the same as the Matrix assignment operator M=S for
    /// a scalar type S. It is overriden for Vector and Row types to behave
    /// as elementwiseScalarAssign.
    template <class S> inline MatrixBase&
    scalarAssign(const S& s) {
        setToZero(); updDiag().setTo(s);
        return *this;
    }

    /// Add a scalar to M's diagonal. This is the same as the Matrix +=
    /// operator. This is overridden for Vector and Row types to behave
    /// as elementwiseAddScalarInPlace.
    template <class S> inline MatrixBase&
    scalarAddInPlace(const S& s) {
        updDiag().elementwiseAddScalarInPlace(s);
    }


    /// Subtract a scalar from M's diagonal. This is the same as the Matrix -=
    /// operator. This is overridden for Vector and Row types to behave
    /// as elementwiseSubtractScalarInPlace.
    template <class S> inline MatrixBase&
    scalarSubtractInPlace(const S& s) {
        updDiag().elementwiseSubtractScalarInPlace(s);
    }

    /// Set M(i,i) = S - M(i,i), M(i,j) = -M(i,j) for i!=j. This is overridden
    /// for Vector and Row types to behave as elementwiseSubtractFromScalarInPlace.
    template <class S> inline MatrixBase&
    scalarSubtractFromLeftInPlace(const S& s) {
        negateInPlace();
        updDiag().elementwiseAddScalarInPlace(s); // yes, add
    }

    /// Set M(i,j) = M(i,j)*S for some "scalar" S. Actually S can be any
    /// type for which E = E*S makes sense. That is, S must be conformant
    /// with E and it must be possible to store the result back in an E.
    /// This is the *= operator for M *= S and behaves the same way for
    /// Matrix, Vector, and RowVector: every element gets multiplied in
    /// place on the right by S.
    template <class S> inline MatrixBase&
    scalarMultiplyInPlace(const S&);

    /// Set M(i,j) = S * M(i,j) for some "scalar" S. This is the same
    /// as the above routine if S really is a scalar, but for S a more
    /// complicated CNT it will be different.
    template <class S> inline MatrixBase&
    scalarMultiplyFromLeftInPlace(const S&);

    /// Set M(i,j) = M(i,j)/S for some "scalar" S. Actually S can be any
    /// type for which E = E/S makes sense. That is, S^-1 must be conformant
    /// with E and it must be possible to store the result back in an E.
    /// This is the /= operator for M /= S and behaves the same way for
    /// Matrix, Vector, and RowVector: every element gets divided in
    /// place on the right by S.
    template <class S> inline MatrixBase&
    scalarDivideInPlace(const S&);

    /// Set M(i,j) = S/M(i,j) for some "scalar" S. Actually S can be any
    /// type for which E = S/E makes sense. That is, S must be conformant
    /// with E^-1 and it must be possible to store the result back in an E.
    template <class S> inline MatrixBase&
    scalarDivideFromLeftInPlace(const S&);


    /// M = diag(r) * M; r must have nrow() elements.
    /// That is, M[i] *= r[i].
    template <class EE> inline MatrixBase& 
    rowScaleInPlace(const VectorBase<EE>&);

    /// Return type is a new matrix which will have the same dimensions as 'this' but
    /// will have element types appropriate for the elementwise multiply being performed.
    template <class EE> inline void 
    rowScale(const VectorBase<EE>& r, typename EltResult<EE>::Mul& out) const;

    template <class EE> inline typename EltResult<EE>::Mul 
    rowScale(const VectorBase<EE>& r) const {
        typename EltResult<EE>::Mul out(nrow(), ncol()); rowScale(r,out); return out;
    }

    /// M = M * diag(c); c must have ncol() elements.
    /// That is, M(j) *= c[j].
    template <class EE> inline MatrixBase& 
    colScaleInPlace(const VectorBase<EE>&);

    template <class EE> inline void 
    colScale(const VectorBase<EE>& c, typename EltResult<EE>::Mul& out) const;

    template <class EE> inline typename EltResult<EE>::Mul
    colScale(const VectorBase<EE>& c) const {
        typename EltResult<EE>::Mul out(nrow(), ncol()); colScale(c,out); return out;
    }


    /// M = diag(r) * M * diag(c); r must have nrow() elements;  must have ncol() elements.
    /// That is, M(i,j) *= r[i]*c[j].
    /// Having a combined row & column scaling operator means we can go through the matrix
    /// memory once instead of twice.
    template <class ER, class EC> inline MatrixBase& 
    rowAndColScaleInPlace(const VectorBase<ER>& r, const VectorBase<EC>& c);

    template <class ER, class EC> inline void 
    rowAndColScale(const VectorBase<ER>& r, const VectorBase<EC>& c, 
                   typename EltResult<typename VectorBase<ER>::template EltResult<EC>::Mul>::Mul& out) const;

    template <class ER, class EC> inline typename EltResult<typename VectorBase<ER>::template EltResult<EC>::Mul>::Mul
    rowAndColScale(const VectorBase<ER>& r, const VectorBase<EC>& c) const {
        typename EltResult<typename VectorBase<ER>::template EltResult<EC>::Mul>::Mul 
            out(nrow(), ncol()); 
        rowAndColScale(r,c,out); return out;
    }

    /// Set M(i,j)=s for every element of M and some value s. This requires only
    /// that s be assignment compatible with M's elements; s doesn't
    /// actually have to be a scalar. Note that for Matrix types this behavior
    /// is different than scalar assignment, which puts the scalar only on M's
    /// diagonal and sets the rest of M to zero. For Vector and RowVector types,
    /// this operator is identical to the normal assignment operator and
    /// scalarAssignInPlace() method which also assign the scalar to every element.
    template <class S> inline MatrixBase&
    elementwiseAssign(const S& s);

    /// Overloaded to allow an integer argument, which is converted to Real.
    MatrixBase& elementwiseAssign(int s)
    {   return elementwiseAssign<Real>(Real(s)); }

    /// Set M(i,j) = M(i,j)^-1.
    MatrixBase& elementwiseInvertInPlace();

    void elementwiseInvert(MatrixBase<typename CNT<E>::TInvert>& out) const;

    MatrixBase<typename CNT<E>::TInvert> elementwiseInvert() const {
        MatrixBase<typename CNT<E>::TInvert> out(nrow(), ncol());
        elementwiseInvert(out);
        return out;
    }

    /// Set M(i,j)+=s for every element of M and some value s. This requires that s be
    /// conformant with M's elements (of type E) and that the result can
    /// be stored in an E. For Matrix types this behavior is different than
    /// the normal += or scalarAddInPlace() operators, which add the scalar
    /// only to the Matrix diagonal. For Vector and RowVector, this operator
    /// is identical to += and scalarAddInPlace() which also add the scalar
    /// to every element.
    template <class S> inline MatrixBase&
    elementwiseAddScalarInPlace(const S& s);

    template <class S> inline void
    elementwiseAddScalar(const S& s, typename EltResult<S>::Add&) const;

    template <class S> inline typename EltResult<S>::Add
    elementwiseAddScalar(const S& s) const {
        typename EltResult<S>::Add out(nrow(), ncol());
        elementwiseAddScalar(s,out);
        return out;
    }

    /// Set M(i,j)-=s for every element of M and some value s. This requires that s be
    /// conformant with M's elements (of type E) and that the result can
    /// be stored in an E. For Matrix types this behavior is different than
    /// the normal -= or scalarSubtractInPlace() operators, which subtract the scalar
    /// only from the Matrix diagonal. For Vector and RowVector, this operator
    /// is identical to -= and scalarSubtractInPlace() which also subtract the scalar
    /// from every element.
    template <class S> inline MatrixBase&
    elementwiseSubtractScalarInPlace(const S& s);

    template <class S> inline void
    elementwiseSubtractScalar(const S& s, typename EltResult<S>::Sub&) const;

    template <class S> inline typename EltResult<S>::Sub
    elementwiseSubtractScalar(const S& s) const {
        typename EltResult<S>::Sub out(nrow(), ncol());
        elementwiseSubtractScalar(s,out);
        return out;
    }

    /// Set M(i,j) = s - M(i,j) for every element of M and some value s. This requires that s be
    /// conformant with M's elements (of type E) and that the result can
    /// be stored in an E. For Matrix types this behavior is different than
    /// the scalarSubtractFromLeftInPlace() operator, which subtracts only the diagonal
    /// elements of M from s, while simply negating the off diagonal elements.
    /// For Vector and RowVector, this operator
    /// is identical to scalarSubtractFromLeftInPlace() which also subtracts every
    /// element of M from the scalar.
    template <class S> inline MatrixBase&
    elementwiseSubtractFromScalarInPlace(const S& s);

    template <class S> inline void
    elementwiseSubtractFromScalar(
        const S&, 
        typename MatrixBase<S>::template EltResult<E>::Sub&) const;

    template <class S> inline typename MatrixBase<S>::template EltResult<E>::Sub
    elementwiseSubtractFromScalar(const S& s) const {
        typename MatrixBase<S>::template EltResult<E>::Sub out(nrow(), ncol());
        elementwiseSubtractFromScalar<S>(s,out);
        return out;
    }

    /// M(i,j) *= R(i,j); R must have same dimensions as this.
    template <class EE> inline MatrixBase& 
    elementwiseMultiplyInPlace(const MatrixBase<EE>&);

    template <class EE> inline void 
    elementwiseMultiply(const MatrixBase<EE>&, typename EltResult<EE>::Mul&) const;

    template <class EE> inline typename EltResult<EE>::Mul 
    elementwiseMultiply(const MatrixBase<EE>& m) const {
        typename EltResult<EE>::Mul out(nrow(), ncol()); 
        elementwiseMultiply<EE>(m,out); 
        return out;
    }

    /// M(i,j) = R(i,j) * M(i,j); R must have same dimensions as this.
    template <class EE> inline MatrixBase& 
    elementwiseMultiplyFromLeftInPlace(const MatrixBase<EE>&);

    template <class EE> inline void 
    elementwiseMultiplyFromLeft(
        const MatrixBase<EE>&, 
        typename MatrixBase<EE>::template EltResult<E>::Mul&) const;

    template <class EE> inline typename MatrixBase<EE>::template EltResult<E>::Mul 
    elementwiseMultiplyFromLeft(const MatrixBase<EE>& m) const {
        typename EltResult<EE>::Mul out(nrow(), ncol()); 
        elementwiseMultiplyFromLeft<EE>(m,out); 
        return out;
    }

    /// M(i,j) /= R(i,j); R must have same dimensions as this.
    template <class EE> inline MatrixBase& 
    elementwiseDivideInPlace(const MatrixBase<EE>&);

    template <class EE> inline void 
    elementwiseDivide(const MatrixBase<EE>&, typename EltResult<EE>::Dvd&) const;

    template <class EE> inline typename EltResult<EE>::Dvd 
    elementwiseDivide(const MatrixBase<EE>& m) const {
        typename EltResult<EE>::Dvd out(nrow(), ncol()); 
        elementwiseDivide<EE>(m,out); 
        return out;
    }

    /// M(i,j) = R(i,j) / M(i,j); R must have same dimensions as this.
    template <class EE> inline MatrixBase& 
    elementwiseDivideFromLeftInPlace(const MatrixBase<EE>&);

    template <class EE> inline void 
    elementwiseDivideFromLeft(
        const MatrixBase<EE>&,
        typename MatrixBase<EE>::template EltResult<E>::Dvd&) const;

    template <class EE> inline typename MatrixBase<EE>::template EltResult<EE>::Dvd 
    elementwiseDivideFromLeft(const MatrixBase<EE>& m) const {
        typename MatrixBase<EE>::template EltResult<E>::Dvd out(nrow(), ncol()); 
        elementwiseDivideFromLeft<EE>(m,out); 
        return out;
    }

    /// Fill every element in current allocation with given element (or NaN or 0).
    MatrixBase& setTo(const ELT& t) {helper.fillWith(reinterpret_cast<const Scalar*>(&t)); return *this;}
    MatrixBase& setToNaN() {helper.fillWithScalar(CNT<StdNumber>::getNaN()); return *this;}
    MatrixBase& setToZero() {helper.fillWithScalar(StdNumber(0)); return *this;}
 
    // View creating operators.   
    inline RowVectorView_<ELT> row(int i) const;   // select a row
    inline RowVectorView_<ELT> updRow(int i);
    inline VectorView_<ELT>    col(int j) const;   // select a column
    inline VectorView_<ELT>    updCol(int j);

    RowVectorView_<ELT> operator[](int i) const {return row(i);}
    RowVectorView_<ELT> operator[](int i)       {return updRow(i);}
    VectorView_<ELT>    operator()(int j) const {return col(j);}
    VectorView_<ELT>    operator()(int j)       {return updCol(j);}
     
    // Select a block.
    inline MatrixView_<ELT> block(int i, int j, int m, int n) const;
    inline MatrixView_<ELT> updBlock(int i, int j, int m, int n);

    MatrixView_<ELT> operator()(int i, int j, int m, int n) const
      { return block(i,j,m,n); }
    MatrixView_<ELT> operator()(int i, int j, int m, int n)
      { return updBlock(i,j,m,n); }
 
    // Hermitian transpose.
    inline MatrixView_<EHerm> transpose() const;
    inline MatrixView_<EHerm> updTranspose();

    MatrixView_<EHerm> operator~() const {return transpose();}
    MatrixView_<EHerm> operator~()       {return updTranspose();}

    /// Select main diagonal (of largest leading square if rectangular) and
    /// return it as a read-only view of the diagonal elements of this Matrix.
    inline VectorView_<ELT> diag() const;
    /// Select main diagonal (of largest leading square if rectangular) and
    /// return it as a writable view of the diagonal elements of this Matrix.
    inline VectorView_<ELT> updDiag();
    /// This non-const version of diag() is an alternate name for updDiag()
    /// available for historical reasons.
    VectorView_<ELT> diag() {return updDiag();}

    // Create a view of the real or imaginary elements. TODO
    //inline MatrixView_<EReal> real() const;
    //inline MatrixView_<EReal> updReal();
    //inline MatrixView_<EImag> imag() const;
    //inline MatrixView_<EImag> updImag();

    // Overload "real" and "imag" for both read and write as a nod to convention. TODO
    //MatrixView_<EReal> real() {return updReal();}
    //MatrixView_<EReal> imag() {return updImag();}

    // TODO: this routine seems ill-advised but I need it for the IVM port at the moment
    TInvert invert() const {  // return a newly-allocated inverse; dump negator 
        TInvert m(*this);
        m.helper.invertInPlace();
        return m;   // TODO - bad: makes an extra copy
    }

    void invertInPlace() {helper.invertInPlace();}

    /// Matlab-compatible debug output.
    void dump(const char* msg=0) const {
        helper.dump(msg);
    }

    /// Element selection for stored elements. These are the fastest element access
    /// methods but may not be able to access all elements of the logical matrix when
    /// some of its elements are not stored in memory. For example, a Hermitian matrix
    /// stores only half its elements and other ones have to be calculated by conjugation
    /// if they are to be returned as type ELT. (You can get them for free by recasting
    /// the matrix so that the elements are reinterpreted as conjugates.) If you want
    /// to guarantee that you can access the value of every element of a matrix, stored or not,
    /// use getAnyElt() instead.
    const ELT& getElt(int i, int j) const { return *reinterpret_cast<const ELT*>(helper.getElt(i,j)); }
    ELT&       updElt(int i, int j)       { return *reinterpret_cast<      ELT*>(helper.updElt(i,j)); }

    const ELT& operator()(int i, int j) const {return getElt(i,j);}
    ELT&       operator()(int i, int j)       {return updElt(i,j);}

    /// This returns a copy of the element value for any position in the logical matrix,
    /// regardless of whether it is stored in memory. If necessary the element's value
    /// is calculated. This is much slower than getElt() but less restrictive.
    /// @see getElt()
    void getAnyElt(int i, int j, ELT& value) const
    {   helper.getAnyElt(i,j,reinterpret_cast<Scalar*>(&value)); }
    ELT getAnyElt(int i, int j) const {ELT e; getAnyElt(i,j,e); return e;}

    /// Scalar norm square is sum( squares of all scalars ). Note that this
    /// is not very useful unless the elements are themselves scalars.
    // TODO: very slow! Should be optimized at least for the case
    //       where ELT is a Scalar.
    ScalarNormSq scalarNormSqr() const {
        const int nr=nrow(), nc=ncol();
        ScalarNormSq sum(0);
        for(int j=0;j<nc;++j) 
            for (int i=0; i<nr; ++i)
                sum += CNT<E>::scalarNormSqr((*this)(i,j));
        return sum;
    }

    /// abs() is elementwise absolute value; that is, the return value has the same
    /// dimension as this Matrix but with each element replaced by whatever it thinks
    /// its absolute value is.
    // TODO: very slow! Should be optimized at least for the case
    //       where ELT is a Scalar.
    void abs(TAbs& mabs) const {
        const int nr=nrow(), nc=ncol();
        mabs.resize(nr,nc);
        for(int j=0;j<nc;++j) 
            for (int i=0; i<nr; ++i)
                mabs(i,j) = CNT<E>::abs((*this)(i,j));
    }

    /// abs() with the result as a function return. More convenient than the 
    /// other abs() member function, but may involve an additional copy of the 
    /// matrix.
    TAbs abs() const { TAbs mabs; abs(mabs); return mabs; }

    /// Return a Matrix of the same shape and contents as this one but
    /// with the element type converted to one based on the standard
    /// C++ scalar types: float, double, complex<float>,
    /// complex<double>. That is, negator<>
    /// and conjugate<> are eliminated from the element type by 
    /// performing any needed negations computationally.
    /// Note that this is actually producing a new matrix with new data;
    /// you can also do this for free by reinterpreting the current
    /// matrix as a different type, if you don't mind looking at
    /// shared data.
    TStandard standardize() const {
        const int nr=nrow(), nc=ncol();
        TStandard mstd(nr, nc);
        for(int j=0;j<nc;++j) 
            for (int i=0; i<nr; ++i)
                mstd(i,j) = CNT<E>::standardize((*this)(i,j));
        return mstd;
    }

    /// This is the scalar Frobenius norm, and its square. Note: if this is a 
    /// Matrix then the Frobenius norm is NOT the same as the 2-norm, although
    /// they are equivalent for Vectors.
    ScalarNormSq normSqr() const { return scalarNormSqr(); }
    // TODO -- not good; unnecessary overflow
    typename CNT<ScalarNormSq>::TSqrt 
        norm() const { return CNT<ScalarNormSq>::sqrt(scalarNormSqr()); }

    /// We only allow RMS norm if the elements are scalars. If there are no 
    /// elements in this Matrix, we'll define its RMS norm to be 0, although 
    /// NaN might be a better choice.
    typename CNT<ScalarNormSq>::TSqrt 
    normRMS() const {
        if (!CNT<ELT>::IsScalar)
            SimTK_THROW1(Exception::Cant, "normRMS() only defined for scalar elements");
        if (nelt() == 0)
            return typename CNT<ScalarNormSq>::TSqrt(0);
        return CNT<ScalarNormSq>::sqrt(scalarNormSqr()/nelt());
    }

    /// Form the column sums of this matrix, returned as a RowVector.
    RowVector_<ELT> colSum() const {
        const int cols = ncol();
        RowVector_<ELT> row(cols);
        for (int j = 0; j < cols; ++j)
            helper.colSum(j, reinterpret_cast<Scalar*>(&row[j]));
        return row;
    }
    /// Alternate name for colSum(); behaves like the Matlab function sum().
    RowVector_<ELT> sum() const {return colSum();}

    /// Form the row sums of this matrix, returned as a Vector.
    Vector_<ELT> rowSum() const {
        const int rows = nrow();
        Vector_<ELT> col(rows);
        for (int i = 0; i < rows; ++i)
            helper.rowSum(i, reinterpret_cast<Scalar*>(&col[i]));
        return col;
    }

    //TODO: make unary +/- return a self-view so they won't reallocate?
    const MatrixBase& operator+() const {return *this; }
    const TNeg&       negate()    const {return *reinterpret_cast<const TNeg*>(this); }
    TNeg&             updNegate()       {return *reinterpret_cast<TNeg*>(this); }

    const TNeg&       operator-() const {return negate();}
    TNeg&             operator-()       {return updNegate();}

    MatrixBase& negateInPlace() {(*this) *= EPrecision(-1); return *this;}
 
    /// Change the size of this matrix. This is only allowed for owner matrices. The
    /// current storage format is retained, but all the data is lost. If you want
    /// to keep the old data, use resizeKeep().
    /// @see resizeKeep()
    MatrixBase& resize(int m, int n)     { helper.resize(m,n); return *this; }
    /// Change the size of this matrix, retaining as much of the old data as will
    /// fit. This is only allowed for owner matrices. The
    /// current storage format is retained, and the existing data is copied
    /// into the new memory to the extent that it will fit.
    /// @see resize()
    MatrixBase& resizeKeep(int m, int n) { helper.resizeKeep(m,n); return *this; }

    // This prevents shape changes in a Matrix that would otherwise allow it. No harm if is
    // are called on a Matrix that is locked already; it always succeeds.
    void lockShape() {helper.lockShape();}

    // This allows shape changes again for a Matrix which was constructed to allow them
    // but had them locked with the above routine. No harm if this is called on a Matrix
    // that is already unlocked, but it is not allowed to call this on a Matrix which
    // *never* allowed resizing. An exception will be thrown in that case.
    void unlockShape() {helper.unlockShape();}
    
    // An assortment of handy conversions
    const MatrixView_<ELT>& getAsMatrixView() const { return *reinterpret_cast<const MatrixView_<ELT>*>(this); }
    MatrixView_<ELT>&       updAsMatrixView()       { return *reinterpret_cast<      MatrixView_<ELT>*>(this); } 
    const Matrix_<ELT>&     getAsMatrix()     const { return *reinterpret_cast<const Matrix_<ELT>*>(this); }
    Matrix_<ELT>&           updAsMatrix()           { return *reinterpret_cast<      Matrix_<ELT>*>(this); }
         
    const VectorView_<ELT>& getAsVectorView() const 
      { assert(ncol()==1); return *reinterpret_cast<const VectorView_<ELT>*>(this); }
    VectorView_<ELT>&       updAsVectorView()       
      { assert(ncol()==1); return *reinterpret_cast<      VectorView_<ELT>*>(this); } 
    const Vector_<ELT>&     getAsVector()     const 
      { assert(ncol()==1); return *reinterpret_cast<const Vector_<ELT>*>(this); }
    Vector_<ELT>&           updAsVector()           
      { assert(ncol()==1); return *reinterpret_cast<      Vector_<ELT>*>(this); }
    const VectorBase<ELT>& getAsVectorBase() const 
      { assert(ncol()==1); return *reinterpret_cast<const VectorBase<ELT>*>(this); }
    VectorBase<ELT>&       updAsVectorBase()       
      { assert(ncol()==1); return *reinterpret_cast<      VectorBase<ELT>*>(this); } 
                
    const RowVectorView_<ELT>& getAsRowVectorView() const 
      { assert(nrow()==1); return *reinterpret_cast<const RowVectorView_<ELT>*>(this); }
    RowVectorView_<ELT>&       updAsRowVectorView()       
      { assert(nrow()==1); return *reinterpret_cast<      RowVectorView_<ELT>*>(this); } 
    const RowVector_<ELT>&     getAsRowVector()     const 
      { assert(nrow()==1); return *reinterpret_cast<const RowVector_<ELT>*>(this); }
    RowVector_<ELT>&           updAsRowVector()           
      { assert(nrow()==1); return *reinterpret_cast<      RowVector_<ELT>*>(this); }        
    const RowVectorBase<ELT>& getAsRowVectorBase() const 
      { assert(nrow()==1); return *reinterpret_cast<const RowVectorBase<ELT>*>(this); }
    RowVectorBase<ELT>&       updAsRowVectorBase()       
      { assert(nrow()==1); return *reinterpret_cast<      RowVectorBase<ELT>*>(this); } 

    // Access to raw data. We have to return the raw data
    // pointer as pointer-to-scalar because we may pack the elements tighter
    // than a C++ array would.

    /// This is the number of consecutive scalars used to represent one
    /// element of type ELT. This may be fewer than C++ would use for the
    /// element, since it may introduce some padding.
    int getNScalarsPerElement()  const {return NScalarsPerElement;}

    /// This is like sizeof(ELT), but returning the number of bytes \e we use
    /// to store the element which may be fewer than what C++ would use. We store
    /// these packed elements adjacent to one another in memory.
    int getPackedSizeofElement() const {return NScalarsPerElement*sizeof(Scalar);}

    bool hasContiguousData() const {return helper.hasContiguousData();}
    ptrdiff_t getContiguousScalarDataLength() const {
        return helper.getContiguousDataLength();
    }
    const Scalar* getContiguousScalarData() const {
        return helper.getContiguousData();
    }
    Scalar* updContiguousScalarData() {
        return helper.updContiguousData();
    }
    void replaceContiguousScalarData(Scalar* newData, ptrdiff_t length, bool takeOwnership) {
        helper.replaceContiguousData(newData,length,takeOwnership);
    }
    void replaceContiguousScalarData(const Scalar* newData, ptrdiff_t length) {
        helper.replaceContiguousData(newData,length);
    }
    void swapOwnedContiguousScalarData(Scalar* newData, ptrdiff_t length, Scalar*& oldData) {
        helper.swapOwnedContiguousData(newData,length,oldData);
    }

    /// Helper rep-stealing constructor. We take over ownership of this rep here. Note
    /// that this \e defines the handle commitment for this handle. This is intended
    /// for internal use only -- don't call this constructor unless you really 
    /// know what you're doing.
    explicit MatrixBase(MatrixHelperRep<Scalar>* hrep) : helper(hrep) {}

protected:
    const MatrixHelper<Scalar>& getHelper() const {return helper;}
    MatrixHelper<Scalar>&       updHelper()       {return helper;}

private:
    MatrixHelper<Scalar> helper; // this is just one pointer

    template <class EE> friend class MatrixBase;

    // ============================= Unimplemented =============================
    // This routine is useful for implementing friendlier Matrix expressions and operators.
    // It maps closely to the Level-3 BLAS family of pxxmm() routines like sgemm(). The
    // operation performed assumes that "this" is the result, and that "this" has 
    // already been sized correctly to receive the result. We'll compute
    //     this = beta*this + alpha*A*B
    // If beta is 0 then "this" can be uninitialized. If alpha is 0 we promise not
    // to look at A or B. The routine can work efficiently even if A and/or B are transposed
    // by their views, so an expression like
    //        C += s * ~A * ~B
    // can be performed with the single equivalent call
    //        C.matmul(1., s, Tr(A), Tr(B))
    // where Tr(A) indicates a transposed view of the original A.
    // The ultimate efficiency of this call depends on the data layout and views which
    // are used for the three matrices.
    // NOTE: neither A nor B can be the same matrix as 'this', nor views of the same data
    // which would expose elements of 'this' that will be modified by this operation.
    template <class ELT_A, class ELT_B>
    MatrixBase& matmul(const StdNumber& beta,   // applied to 'this'
                       const StdNumber& alpha, const MatrixBase<ELT_A>& A, const MatrixBase<ELT_B>& B)
    {
        helper.matmul(beta,alpha,A.helper,B.helper);
        return *this;
    }

};

} //namespace SimTK

#endif // SimTK_SIMMATRIX_MATRIXBASE_H_
