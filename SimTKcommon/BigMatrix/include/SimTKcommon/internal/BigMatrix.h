#ifndef SimTK_SIMMATRIX_BIGMATRIX_H_
#define SimTK_SIMMATRIX_BIGMATRIX_H_

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

/** @file
 * This file defines the client side of the SimTK::Matrix classes, which
 * hold medium to large, variable-sized matrices whose elements are packed 
 * SimTK "Composite Numerical Types" (CNTs). Unlike CNTs, the implemention here 
 * is opaque, and almost all properties are captured in the implementation at 
 * run time rather than in the type at compile time. 
 *
 * Every Matrix consists logically of three pieces: 
 *  - the matrix handle 
 *  - the matrix helper
 *  - and the matrix data. 
 *
 * They are organized like this:
 * <pre>
 *      ------------            ------------
 *     |  Handle<E> | -------> |            |
 *      ------------  <------- | Helper<S>  |
 *                             |            |
 *                             |            |          --------~ ~--
 *                             |            | ------> | Data<S> ... |
 *                              ------------           --------~ ~--
 * </pre>
 * The handle is the object actually appearing in SimTK API user programs.
 * It always consists of just a single pointer, pointing to a library-side
 * "helper" object whose implementation is opaque. The handle is templatized
 * by the user's element type, which may be any packed composite numerical
 * type, including scalar types like \c float and \c complex<double>, but also
 * composite types such as \c Vec3 or \c Mat<2,2,Mat<3,3>>. A Matrix handle
 * owns the helper to which it points and must destruct the helper when
 * the handle's destructor is called.
 *
 * The helper, on the other hand, is parameterized only by the underlying scalar 
 * type. There are exactly 12 SimTK scalar types, so all can be instantiated on 
 * the library side leaving the implementation opaque and thus flexible from
 * release to release without compromising binary compatibility. (The scalar 
 * types are: the four C++ standard types float and double, 
 * complex<float>, and complex<double>; the SimTK numbers conjugate<float> and 
 * conjugate<double>; and negator<> types templatized by any of the six
 * numeric types.) The helper contains several kinds of information:
 *  - the underlying scalar type S (as its template parameter)
 *  - the number of scalars in the handle's logical element type E
 *  - whether this is an owner matrix, or just a view
 *  - the handle "commitment"; defining the range of matrix characteristics
 *      to which that handle may refer
 *  - the actual characteristics of the matrix currently represented by
 *      the helper
 *  - a virtual function table full of methods which are aware of the
 *      logical structure of the Matrix and the physical structure of 
 *      the data to support operations such as element indexing
 *  - a pointer to the underlying data, which may be shared with other 
 *      helpers
 *
 * The data itself consists only of scalars
 * S of the same type as the helper's template argument, but different 
 * helpers can look at the same data differently. For examples, when the
 * elements are composite consisting of k scalars, the helper will provide a 
 * view of the data in which its scalars are interpreted in groups of k.
 * Many other reinterpretations of the data are possible and useful, such
 * as a real-valued helper viewing only the real or imaginary part of
 * complex data, or a helper which views the data order as though it were
 * transposed.
 *
 * At most \e one matrix helper owns the matrix data and is responsible
 * for deleting that data when no longer needed. That is called an "owner"
 * helper and its associated handle is an owner handle. Normally the owner 
 * is the handle (and helper) that allocated the data, and
 * in most cases an owner can resize the data at will. Many other handles
 * may reference the same data; those non-owner handles are called "views".
 * Every view may present a different picture of the underlying data. The
 * default view is "whole" meaning that all the elements of the data are 
 * visible, and appear in their normal order. A "transpose" view also shows 
 * all the elements but the matrix dimensions and indices are reversed. 
 * Other common views are "block" to select a sub-block of a matrix, and 
 * "diagonal" which shows only the diagonal of a matrix (as a vector). 
 *
 * NOTE: Destruction of an owner destructs the data it owns
 * *regardless* of the presence of other views into that data! I.e., these
 * are not reference counted. TODO: should we change that?
 * 
 * In some cases there may be no owner helper for a particular piece of 
 * matrix data. That occurs when pre-existing memory, such as a Fortran
 * array, is used to construct a Matrix. In that case all the helpers are
 * views and the data will persist after the destruction of the last
 * referencing helper.
 *                 
 * A Matrix that is the owner of its data will be resized whenever
 * necessary, unless you take active steps to prevent that. For example, if
 * you declare a Vector, the number of rows can resize but the number of
 * columns will be locked at 1. A RowVector does the reverse. You can also
 * explicitly lock the number of rows and/or columns of a matrix to prevent
 * unwanted resizes.
 *
 * Here are the classes and short descriptions:
 * <pre>
 *   MatrixHelper<S>  interface to the opaque implementation, templatized
 *                      by scalar type only
 *   MatrixBase<CNT>  fully templatized client, contains a MatrixHelper
 * </pre>
 *
 * The rest are dataless classes all of which can be interconverted just
 * by recasting. Every one of these classes has a default conversion to
 * type Matrix_<same element type>, so users can write functions that expect
 * a Matrix argument and pass it a Vector, RowVectorView, or whatever.
 *
 * <pre>
 *   VectorBase<CNT>    these are derived from MatrixBase and add no new data,
 *   RowVectorBase<CNT> but change some of the operators and other methods to
 *                        be appropriate for 1d data.
 *
 *   Matrix_<CNT>      2d owner class     (a MatrixBase<CNT>)
 *   Vector_<CNT>      column owner class (a VectorBase<CNT>)
 *   RowVector_<CNT>   row owner class    (a RowVectorBase<CNT>)
 * </pre>
 *
 * Views are exactly the same as the corresponding owner class, but with
 * shallow construction and assignment semantics.
 *
 * <pre>
 *   MatrixView_<CNT>, VectorView_<CNT>, RowVectorView_<CNT>
 * </pre>
 *
 * Dead matrices are owners that are about to be destructed. Anything
 * they own may be taken from them, including the helper and/or
 * the data. This is a very effective performance trick for sequences
 * of operations since it eliminates most of the need for allocating and
 * deallocating temporaries.
 *
 * <pre>
 *   DeadMatrix_<CNT>, DeadVector_<CNT>, DeadRowVector_<CNT>
 * </pre>
 *
 */
// TODO: matrix expression templates for delaying operator execution.

#include "SimTKcommon/Scalar.h"
#include "SimTKcommon/SmallMatrix.h"

#include "SimTKcommon/internal/MatrixHelper.h"
#include "SimTKcommon/internal/MatrixCharacteristics.h"

#include <iostream>
#include <cassert>
#include <complex>
#include <cstddef>
#include <limits>

namespace SimTK {

template <class ELT>    class MatrixBase;
template <class ELT>    class VectorBase;
template <class ELT>    class RowVectorBase;

template <class T=Real> class MatrixView_;
template <class T=Real> class DeadMatrixView_;
template <class T=Real> class Matrix_;
template <class T=Real> class DeadMatrix_;

template <class T=Real> class VectorView_;
template <class T=Real> class DeadVectorView_;
template <class T=Real> class Vector_;
template <class T=Real> class DeadVector_;

template <class T=Real> class RowVectorView_;
template <class T=Real> class DeadRowVectorView_;
template <class T=Real> class RowVector_;
template <class T=Real> class DeadRowVector_;

template <class ELT, class VECTOR_CLASS> class VectorIterator;

//  -------------------------------- MatrixBase --------------------------------
/// Variable-size 2d matrix of Composite Numerical Type (ELT) elements. This is
/// a container of such elements, it is NOT a Composite Numerical Type itself.
/// MatrixBase<ELT> uses MatrixHelper<S> for implementation, where S is 
/// ELT::Scalar, that is, the underlying float, double, long double,
/// complex<float>, negator<conjugate<double>>, 
/// etc. from which ELT is constructed. This is a finite set of which all
/// members are explicitly instantiated in the implementation code, so 
/// clients don't have to know how anything is implemented.
/// 
/// MatrixBase is the only class in the Matrix/Vector family which has any
/// data members (it has exactly one MatrixHelper, which itself consists only
/// of a single pointer to an opaque class). Thus all other objects
/// in this family (that is, derived from MatrixBase) are exactly the same
/// size in memory and may be "reinterpreted" as appropriate. For example,
/// a Vector may be reinterpreted as a Matrix or vice versa, provided runtime
/// requirements are met (e.g., exactly 1 column).
///
/// Unlike the small matrix classes, very little is encoded in the type.
/// Only the element type, and matrix vs. vector vs. row are in the type;
/// everything else like shape, storage layout, and writability are handled
/// at run time. 
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

    typedef MatrixBase<E>                T;
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
 
    // View creating operators. TODO: these should be DeadMatrixViews  
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

    /// abs() with the result as a function return. More convenient than the other
    /// abs() member function, but may involve an additional copy of the matrix.
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

    /// This is the scalar Frobenius norm, and its square. Note: if this is a Matrix then the Frobenius
    /// norm is NOT the same as the 2-norm, although they are equivalent for Vectors.
    ScalarNormSq normSqr() const { return scalarNormSqr(); }
    // TODO -- not good; unnecessary overflow
    typename CNT<ScalarNormSq>::TSqrt 
        norm() const { return CNT<ScalarNormSq>::sqrt(scalarNormSqr()); }

    /// We only allow RMS norm if the elements are scalars. If there are no elements in this Matrix,
    /// we'll define its RMS norm to be 0, although NaN might be a better choice.
    typename CNT<ScalarNormSq>::TSqrt 
    normRMS() const {
        if (!CNT<ELT>::IsScalar)
            SimTK_THROW1(Exception::Cant, "normRMS() only defined for scalar elements");
        if (nelt() == 0)
            return typename CNT<ScalarNormSq>::TSqrt(0);
        return CNT<ScalarNormSq>::sqrt(scalarNormSqr()/nelt());
    }

    /// Form the column sums of this matrix, returned as a RowVector.
    RowVectorBase<ELT> sum() const {
        const int cols = ncol();
        RowVectorBase<ELT> row(cols);
        for (int i = 0; i < cols; ++i)
            helper.colSum(i, reinterpret_cast<Scalar*>(&row[i]));
        return row;
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
};



//  -------------------------------- VectorBase --------------------------------
/// This is a dataless rehash of the MatrixBase class to specialize it for Vectors.
/// This mostly entails overriding a few of the methods. Note that all the MatrixBase
/// operations remain available if you static_cast<> this up to a MatrixBase.
//  ----------------------------------------------------------------------------
template <class ELT> class VectorBase : public MatrixBase<ELT> {
    typedef MatrixBase<ELT>                             Base;
    typedef typename CNT<ELT>::Scalar                   Scalar;
    typedef typename CNT<ELT>::Number                   Number;
    typedef typename CNT<ELT>::StdNumber                StdNumber;
    typedef VectorBase<ELT>                             T;
    typedef VectorBase<typename CNT<ELT>::TAbs>         TAbs;
    typedef VectorBase<typename CNT<ELT>::TNeg>         TNeg;
    typedef RowVectorView_<typename CNT<ELT>::THerm>    THerm;
public:  
    //  ------------------------------------------------------------------------
    /// @name       VectorBase "owner" construction
    ///
    /// These constructors create new VectorBase objects which own their
    /// own data and are (at least by default) resizable. The resulting matrices
    /// are m X 1 with the number of columns locked at 1. If there is any data
    /// allocated but not explicitly initialized, that data will be uninitialized
    /// garbage in Release builds but will be initialized to NaN (at a performance
    /// cost) in Debug builds.
    /// @{

    /// Default constructor makes a 0x1 matrix locked at 1 column; you can
    /// provide an initial allocation if you want.
    explicit VectorBase(int m=0) : Base(MatrixCommitment::Vector(), m, 1) {}

    /// Copy constructor is a deep copy (not appropriate for views!). That
    /// means it creates a new, densely packed vector whose elements are
    /// initialized from the source object.
    VectorBase(const VectorBase& source) : Base(source) {}

    /// Implicit conversion from compatible vector with negated elements.
    VectorBase(const TNeg& source) : Base(source) {}

    /// Construct an owner vector of length m, with each element initialized to
    /// the given value.
    VectorBase(int m, const ELT& initialValue)
    :   Base(MatrixCommitment::Vector(),m,1,initialValue) {}  

    /// Construct an owner vector of length m, with the elements initialized sequentially
    /// from a C++ array of elements which is assumed to be of length m. Note that we
    /// are expecting C++ packing; don't use this to initialize one Simmatrix vector
    /// from another because Simmatrix may pack its elements more densely than C++.
    VectorBase(int m, const ELT* cppInitialValues)
    :   Base(MatrixCommitment::Vector(),m,1,cppInitialValues) {}
    /// @}

    //  ------------------------------------------------------------------------
    /// @name       VectorBase construction from pre-existing data
    ///
    /// Construct a non-resizeable, VectorBase view of externally supplied data. Note that
    /// stride should be interpreted as "the number of scalars between elements" and
    /// for composite elements may have a different value if the source is a C++ array 
    /// of elements vs. a Simmatrix packed data array. We provide constructors for
    /// both read-only and writable external data.
    /// @{

    /// Construct a read-only view of existing data.
    VectorBase(int m, int stride, const Scalar* s)
    :   Base(MatrixCommitment::Vector(m), MatrixCharacter::Vector(m),stride,s) { }
    /// Construct a writable view into existing data.
    VectorBase(int m, int stride, Scalar* s)
    :   Base(MatrixCommitment::Vector(m), MatrixCharacter::Vector(m),stride,s) { }
    /// @}
        
    //  ------------------------------------------------------------------------
    /// @name       VectorBase construction from an existing Helper.
    ///
    /// Create a new VectorBase from an existing helper. Both shallow (view) and deep 
    /// copies are possible. For shallow copies, there is a constructor providing a read-only
    /// view of the original data and one providing a writable view into the original data.
    /// @{

    /// Construct a writable view into the source data.
    VectorBase(MatrixHelper<Scalar>& h, const typename MatrixHelper<Scalar>::ShallowCopy& s) 
    :   Base(MatrixCommitment::Vector(), h,s) { }
    /// Construct a read-only view of the source data.
    VectorBase(const MatrixHelper<Scalar>& h, const typename MatrixHelper<Scalar>::ShallowCopy& s) 
    :   Base(MatrixCommitment::Vector(), h,s) { }
    /// Construct a new owner vector initialized with the data from the source.
    VectorBase(const MatrixHelper<Scalar>& h, const typename MatrixHelper<Scalar>::DeepCopy& d)    
    :   Base(MatrixCommitment::Vector(), h,d) { }
    /// @}

    // This gives the resulting vector type when (v[i] op P) is applied to each element.
    // It will have element types which are the regular composite result of ELT op P.
    template <class P> struct EltResult { 
        typedef VectorBase<typename CNT<ELT>::template Result<P>::Mul> Mul;
        typedef VectorBase<typename CNT<ELT>::template Result<P>::Dvd> Dvd;
        typedef VectorBase<typename CNT<ELT>::template Result<P>::Add> Add;
        typedef VectorBase<typename CNT<ELT>::template Result<P>::Sub> Sub;
    };

    /// Copy assignment is deep copy but behavior depends on type of lhs: if view, rhs
    /// must match. If owner, we reallocate and copy rhs.
    VectorBase& operator=(const VectorBase& b) {
        Base::operator=(b); return *this;
    }

    // default destructor


    VectorBase& operator*=(const StdNumber& t)  { Base::operator*=(t); return *this; }
    VectorBase& operator/=(const StdNumber& t)  { Base::operator/=(t); return *this; }
    VectorBase& operator+=(const VectorBase& r) { Base::operator+=(r); return *this; }
    VectorBase& operator-=(const VectorBase& r) { Base::operator-=(r); return *this; }  


    template <class EE> VectorBase& operator=(const VectorBase<EE>& b) 
      { Base::operator=(b);  return *this; } 
    template <class EE> VectorBase& operator+=(const VectorBase<EE>& b) 
      { Base::operator+=(b); return *this; } 
    template <class EE> VectorBase& operator-=(const VectorBase<EE>& b) 
      { Base::operator-=(b); return *this; } 


    /// Fill current allocation with copies of element. Note that this is not the 
    /// same behavior as assignment for Matrices, where only the diagonal is set (and
    /// everything else is set to zero.)
    VectorBase& operator=(const ELT& t) { Base::setTo(t); return *this; }  

	/// There's only one column here so it's a bit wierd to use rowScale rather than
	/// elementwiseMultiply, but there's nothing really wrong with it. Using colScale
	/// would be really wacky since it is the same as a scalar multiply. We won't support
	/// colScale here except through inheritance where it will not be much use.
	template <class EE> VectorBase& rowScaleInPlace(const VectorBase<EE>& v)
	  { Base::template rowScaleInPlace<EE>(v); return *this; }
	template <class EE> inline void rowScale(const VectorBase<EE>& v, typename EltResult<EE>::Mul& out) const
	  { Base::rowScale(v,out); }
	template <class EE> inline typename EltResult<EE>::Mul rowScale(const VectorBase<EE>& v) const
	  { typename EltResult<EE>::Mul out(nrow()); Base::rowScale(v,out); return out; }

    /// Set this[i] = this[i]^-1.
    VectorBase& elementwiseInvertInPlace() {
        Base::elementwiseInvertInPlace();
        return *this;
    }

    /// Set supplied out[i] = this[i]^-1
    void elementwiseInvert(VectorBase<typename CNT<ELT>::TInvert>& out) const {
        Base::elementwiseInvert(out);
    }

    /// Return out[i]=this[i]^-1 as function return.
    VectorBase<typename CNT<ELT>::TInvert> elementwiseInvert() const {
        VectorBase<typename CNT<ELT>::TInvert> out(nrow());
        Base::elementwiseInvert(out);
        return out;
    }

        // elementwise multiply
	template <class EE> VectorBase& elementwiseMultiplyInPlace(const VectorBase<EE>& r)
	  { Base::template elementwiseMultiplyInPlace<EE>(r); return *this; }
	template <class EE> inline void elementwiseMultiply(const VectorBase<EE>& v, typename EltResult<EE>::Mul& out) const
	  { Base::template elementwiseMultiply<EE>(v,out); }
	template <class EE> inline typename EltResult<EE>::Mul elementwiseMultiply(const VectorBase<EE>& v) const
	  { typename EltResult<EE>::Mul out(nrow()); Base::template elementwiseMultiply<EE>(v,out); return out; }

        // elementwise multiply from left
	template <class EE> VectorBase& elementwiseMultiplyFromLeftInPlace(const VectorBase<EE>& r)
	  { Base::template elementwiseMultiplyFromLeftInPlace<EE>(r); return *this; }
	template <class EE> inline void 
    elementwiseMultiplyFromLeft(
        const VectorBase<EE>& v, 
        typename VectorBase<EE>::template EltResult<ELT>::Mul& out) const
	{ 
        Base::template elementwiseMultiplyFromLeft<EE>(v,out);
    }
	template <class EE> inline typename VectorBase<EE>::template EltResult<ELT>::Mul 
    elementwiseMultiplyFromLeft(const VectorBase<EE>& v) const
	{ 
        typename VectorBase<EE>::template EltResult<ELT>::Mul out(nrow()); 
        Base::template elementwiseMultiplyFromLeft<EE>(v,out); 
        return out;
    }

        // elementwise divide
	template <class EE> VectorBase& elementwiseDivideInPlace(const VectorBase<EE>& r)
	  { Base::template elementwiseDivideInPlace<EE>(r); return *this; }
	template <class EE> inline void elementwiseDivide(const VectorBase<EE>& v, typename EltResult<EE>::Dvd& out) const
	  { Base::template elementwiseDivide<EE>(v,out); }
	template <class EE> inline typename EltResult<EE>::Dvd elementwiseDivide(const VectorBase<EE>& v) const
	  { typename EltResult<EE>::Dvd out(nrow()); Base::template elementwiseDivide<EE>(v,out); return out; }

        // elementwise divide from left
	template <class EE> VectorBase& elementwiseDivideFromLeftInPlace(const VectorBase<EE>& r)
	  { Base::template elementwiseDivideFromLeftInPlace<EE>(r); return *this; }
	template <class EE> inline void 
    elementwiseDivideFromLeft(
        const VectorBase<EE>& v, 
        typename VectorBase<EE>::template EltResult<ELT>::Dvd& out) const
	{ 
        Base::template elementwiseDivideFromLeft<EE>(v,out);
    }
	template <class EE> inline typename VectorBase<EE>::template EltResult<ELT>::Dvd 
    elementwiseDivideFromLeft(const VectorBase<EE>& v) const
	{ 
        typename VectorBase<EE>::template EltResult<ELT>::Dvd out(nrow()); 
        Base::template elementwiseDivideFromLeft<EE>(v,out); 
        return out;
    }

    // Implicit conversions are allowed to Vector or Matrix, but not to RowVector.   
    operator const Vector_<ELT>&()     const { return *reinterpret_cast<const Vector_<ELT>*>(this); }
    operator       Vector_<ELT>&()           { return *reinterpret_cast<      Vector_<ELT>*>(this); }
    operator const VectorView_<ELT>&() const { return *reinterpret_cast<const VectorView_<ELT>*>(this); }
    operator       VectorView_<ELT>&()       { return *reinterpret_cast<      VectorView_<ELT>*>(this); }
    
    operator const Matrix_<ELT>&()     const { return *reinterpret_cast<const Matrix_<ELT>*>(this); }
    operator       Matrix_<ELT>&()           { return *reinterpret_cast<      Matrix_<ELT>*>(this); } 
    operator const MatrixView_<ELT>&() const { return *reinterpret_cast<const MatrixView_<ELT>*>(this); }
    operator       MatrixView_<ELT>&()       { return *reinterpret_cast<      MatrixView_<ELT>*>(this); } 


    // size() for Vectors is Base::nelt() but returns int instead of ptrdiff_t.
	int size() const { 
		assert(Base::nelt() <= (ptrdiff_t)std::numeric_limits<int>::max()); 
		assert(Base::ncol()==1);
		return (int)Base::nelt();
	}
	int       nrow() const {assert(Base::ncol()==1); return Base::nrow();}
	int       ncol() const {assert(Base::ncol()==1); return Base::ncol();}
    ptrdiff_t nelt() const {assert(Base::ncol()==1); return Base::nelt();}

    // Override MatrixBase operators to return the right shape
    TAbs abs() const {TAbs result; Base::abs(result); return result;}
    
    // Override MatrixBase indexing operators          
    const ELT& operator[](int i) const {return *reinterpret_cast<const ELT*>(Base::getHelper().getElt(i));}
    ELT&       operator[](int i)       {return *reinterpret_cast<ELT*>      (Base::updHelper().updElt(i));}
    const ELT& operator()(int i) const {return *reinterpret_cast<const ELT*>(Base::getHelper().getElt(i));}
    ELT&       operator()(int i)       {return *reinterpret_cast<ELT*>      (Base::updHelper().updElt(i));}
         
    // Block (contiguous subvector) view creation      
    VectorView_<ELT> operator()(int i, int m) const {return Base::operator()(i,0,m,1).getAsVectorView();}
    VectorView_<ELT> operator()(int i, int m)       {return Base::operator()(i,0,m,1).updAsVectorView();}

    // Indexed view creation (arbitrary subvector). Indices must be monotonically increasing.
    VectorView_<ELT> index(const Array_<int>& indices) const {
        MatrixHelper<Scalar> h(Base::getHelper().getCharacterCommitment(), Base::getHelper(), indices);
        return VectorView_<ELT>(h);
    }
    VectorView_<ELT> updIndex(const Array_<int>& indices) {
        MatrixHelper<Scalar> h(Base::getHelper().getCharacterCommitment(), Base::updHelper(), indices);
        return VectorView_<ELT>(h);
    }

    VectorView_<ELT> operator()(const Array_<int>& indices) const {return index(indices);}
    VectorView_<ELT> operator()(const Array_<int>& indices)       {return updIndex(indices);}
 
    // Hermitian transpose.
    THerm transpose() const {return Base::transpose().getAsRowVectorView();}
    THerm updTranspose()    {return Base::updTranspose().updAsRowVectorView();}

    THerm operator~() const {return transpose();}
    THerm operator~()       {return updTranspose();}

    const VectorBase& operator+() const {return *this; }

    // Negation

    const TNeg& negate()    const {return *reinterpret_cast<const TNeg*>(this); }
    TNeg&       updNegate()       {return *reinterpret_cast<TNeg*>(this); }

    const TNeg& operator-() const {return negate();}
    TNeg&       operator-()       {return updNegate();}

    VectorBase& resize(int m)     {Base::resize(m,1); return *this;}
    VectorBase& resizeKeep(int m) {Base::resizeKeep(m,1); return *this;}

	//TODO: this is not re-locking the number of columns at 1.
	void clear() {Base::clear(); Base::resize(0,1);}

    ELT sum() const {ELT s; Base::getHelper().sum(reinterpret_cast<Scalar*>(&s)); return s; } // add all the elements        
    VectorIterator<ELT, VectorBase<ELT> > begin() {
        return VectorIterator<ELT, VectorBase<ELT> >(*this, 0);
    }
    VectorIterator<ELT, VectorBase<ELT> > end() {
        return VectorIterator<ELT, VectorBase<ELT> >(*this, size());
    }

protected:
    // Create a VectorBase handle using a given helper rep. 
    explicit VectorBase(MatrixHelperRep<Scalar>* hrep) : Base(hrep) {}

private:
    // NO DATA MEMBERS ALLOWED
};



//  ------------------------------- RowVectorBase ------------------------------
/// This is a dataless rehash of the MatrixBase class to specialize it for RowVectors.
/// This mostly entails overriding a few of the methods. Note that all the MatrixBase
/// operations remain available if you static_cast<> this up to a MatrixBase.
//  ----------------------------------------------------------------------------
template <class ELT> class RowVectorBase : public MatrixBase<ELT> {
    typedef MatrixBase<ELT>                             Base;
    typedef typename CNT<ELT>::Scalar                   Scalar;
    typedef typename CNT<ELT>::Number                   Number;
    typedef typename CNT<ELT>::StdNumber                StdNumber;
    typedef RowVectorBase<ELT>                          T;
    typedef RowVectorBase<typename CNT<ELT>::TAbs>      TAbs;
    typedef RowVectorBase<typename CNT<ELT>::TNeg>      TNeg;
    typedef VectorView_<typename CNT<ELT>::THerm>       THerm;
public: 
    //  ------------------------------------------------------------------------
    /// @name       RowVectorBase "owner" construction
    ///
    /// These constructors create new RowVectorBase objects which own their
    /// own data and are (at least by default) resizable. The resulting matrices
    /// are 1 x n with the number of rows locked at 1. If there is any data
    /// allocated but not explicitly initialized, that data will be uninitialized
    /// garbage in Release builds but will be initialized to NaN (at a performance
    /// cost) in Debug builds.
    /// @{

    /// Default constructor makes a 1x0 matrix locked at 1 row; you can
    /// provide an initial allocation if you want.
    explicit RowVectorBase(int n=0) : Base(MatrixCommitment::RowVector(), 1, n) {}
    
    /// Copy constructor is a deep copy (not appropriate for views!). That
    /// means it creates a new, densely packed vector whose elements are
    /// initialized from the source object.    
    RowVectorBase(const RowVectorBase& source) : Base(source) {}

    /// Implicit conversion from compatible row vector with negated elements.
    RowVectorBase(const TNeg& source) : Base(source) {}

    /// Construct an owner row vector of length n, with each element initialized to
    /// the given value.
    RowVectorBase(int n, const ELT& initialValue)
    :   Base(MatrixCommitment::RowVector(),1,n,initialValue) {}  

    /// Construct an owner vector of length n, with the elements initialized sequentially
    /// from a C++ array of elements which is assumed to be of length n. Note that we
    /// are expecting C++ packing; don't use this to initialize one Simmatrix vector
    /// from another because Simmatrix may pack its elements more densely than C++.
    RowVectorBase(int n, const ELT* cppInitialValues)
    :   Base(MatrixCommitment::RowVector(),1,n,cppInitialValues) {}
    /// @}

    //  ------------------------------------------------------------------------
    /// @name       RowVectorBase construction from pre-existing data
    ///
    /// Construct a non-resizeable, RowVectorBase view of externally supplied data. Note that
    /// stride should be interpreted as "the number of scalars between elements" and
    /// for composite elements may have a different value if the source is a C++ array 
    /// of elements vs. a Simmatrix packed data array. We provide constructors for
    /// both read-only and writable external data.
    /// @{

    /// Construct a read-only view of existing data.
    RowVectorBase(int n, int stride, const Scalar* s)
    :   Base(MatrixCommitment::RowVector(n), MatrixCharacter::RowVector(n),stride,s) { }
    /// Construct a writable view into existing data.
    RowVectorBase(int n, int stride, Scalar* s)
    :   Base(MatrixCommitment::RowVector(n), MatrixCharacter::RowVector(n),stride,s) { }
    /// @}

    //  ------------------------------------------------------------------------
    /// @name       RowVectorBase construction from an existing Helper.
    ///
    /// Create a new RowVectorBase from an existing helper. Both shallow (view) and deep 
    /// copies are possible. For shallow copies, there is a constructor providing a read-only
    /// view of the original data and one providing a writable view into the original data.
    /// @{

    /// Construct a writable view into the source data.
    RowVectorBase(MatrixHelper<Scalar>& h, const typename MatrixHelper<Scalar>::ShallowCopy& s) 
    :   Base(MatrixCommitment::RowVector(), h,s) { }
    /// Construct a read-only view of the source data.
    RowVectorBase(const MatrixHelper<Scalar>& h, const typename MatrixHelper<Scalar>::ShallowCopy& s) 
    :   Base(MatrixCommitment::RowVector(), h,s) { }
    /// Construct a new owner vector initialized with the data from the source.
    RowVectorBase(const MatrixHelper<Scalar>& h, const typename MatrixHelper<Scalar>::DeepCopy& d)    
    :   Base(MatrixCommitment::RowVector(), h,d) { }
    /// @}

	// This gives the resulting rowvector type when (r(i) op P) is applied to each element.
    // It will have element types which are the regular composite result of ELT op P.
    template <class P> struct EltResult { 
        typedef RowVectorBase<typename CNT<ELT>::template Result<P>::Mul> Mul;
        typedef RowVectorBase<typename CNT<ELT>::template Result<P>::Dvd> Dvd;
        typedef RowVectorBase<typename CNT<ELT>::template Result<P>::Add> Add;
        typedef RowVectorBase<typename CNT<ELT>::template Result<P>::Sub> Sub;
    };

    /// Copy assignment is deep copy but behavior depends on type of lhs: if view, rhs
    /// must match. If owner, we reallocate and copy rhs.
    RowVectorBase& operator=(const RowVectorBase& b) {
        Base::operator=(b); return *this;
    }

    // default destructor

    RowVectorBase& operator*=(const StdNumber& t)     {Base::operator*=(t); return *this;}
    RowVectorBase& operator/=(const StdNumber& t)     {Base::operator/=(t); return *this;}
    RowVectorBase& operator+=(const RowVectorBase& r) {Base::operator+=(r); return *this;}
    RowVectorBase& operator-=(const RowVectorBase& r) {Base::operator-=(r); return *this;}  

    template <class EE> RowVectorBase& operator=(const RowVectorBase<EE>& b) 
      { Base::operator=(b);  return *this; } 
    template <class EE> RowVectorBase& operator+=(const RowVectorBase<EE>& b) 
      { Base::operator+=(b); return *this; } 
    template <class EE> RowVectorBase& operator-=(const RowVectorBase<EE>& b) 
      { Base::operator-=(b); return *this; } 

    // default destructor
 
    /// Fill current allocation with copies of element. Note that this is not the 
    /// same behavior as assignment for Matrices, where only the diagonal is set (and
    /// everything else is set to zero.)
    RowVectorBase& operator=(const ELT& t) { Base::setTo(t); return *this; } 

	/// There's only one row here so it's a bit wierd to use colScale rather than
	/// elementwiseMultiply, but there's nothing really wrong with it. Using rowScale
	/// would be really wacky since it is the same as a scalar multiply. We won't support
	/// rowScale here except through inheritance where it will not be much use.
	template <class EE> RowVectorBase& colScaleInPlace(const VectorBase<EE>& v)
	  { Base::template colScaleInPlace<EE>(v); return *this; }
	template <class EE> inline void colScale(const VectorBase<EE>& v, typename EltResult<EE>::Mul& out) const
	  { return Base::template colScale<EE>(v,out); }
	template <class EE> inline typename EltResult<EE>::Mul colScale(const VectorBase<EE>& v) const
	  { typename EltResult<EE>::Mul out(ncol()); Base::template colScale<EE>(v,out); return out; }


        // elementwise multiply
	template <class EE> RowVectorBase& elementwiseMultiplyInPlace(const RowVectorBase<EE>& r)
	  { Base::template elementwiseMultiplyInPlace<EE>(r); return *this; }
	template <class EE> inline void elementwiseMultiply(const RowVectorBase<EE>& v, typename EltResult<EE>::Mul& out) const
	  { Base::template elementwiseMultiply<EE>(v,out); }
	template <class EE> inline typename EltResult<EE>::Mul elementwiseMultiply(const RowVectorBase<EE>& v) const
	  { typename EltResult<EE>::Mul out(nrow()); Base::template elementwiseMultiply<EE>(v,out); return out; }

        // elementwise multiply from left
	template <class EE> RowVectorBase& elementwiseMultiplyFromLeftInPlace(const RowVectorBase<EE>& r)
	  { Base::template elementwiseMultiplyFromLeftInPlace<EE>(r); return *this; }
	template <class EE> inline void 
    elementwiseMultiplyFromLeft(
        const RowVectorBase<EE>& v, 
        typename RowVectorBase<EE>::template EltResult<ELT>::Mul& out) const
	{ 
        Base::template elementwiseMultiplyFromLeft<EE>(v,out);
    }
	template <class EE> inline typename RowVectorBase<EE>::template EltResult<ELT>::Mul 
    elementwiseMultiplyFromLeft(const RowVectorBase<EE>& v) const
	{ 
        typename RowVectorBase<EE>::template EltResult<ELT>::Mul out(nrow()); 
        Base::template elementwiseMultiplyFromLeft<EE>(v,out); 
        return out;
    }

        // elementwise divide
	template <class EE> RowVectorBase& elementwiseDivideInPlace(const RowVectorBase<EE>& r)
	  { Base::template elementwiseDivideInPlace<EE>(r); return *this; }
	template <class EE> inline void elementwiseDivide(const RowVectorBase<EE>& v, typename EltResult<EE>::Dvd& out) const
	  { Base::template elementwiseDivide<EE>(v,out); }
	template <class EE> inline typename EltResult<EE>::Dvd elementwiseDivide(const RowVectorBase<EE>& v) const
	  { typename EltResult<EE>::Dvd out(nrow()); Base::template elementwiseDivide<EE>(v,out); return out; }

        // elementwise divide from left
	template <class EE> RowVectorBase& elementwiseDivideFromLeftInPlace(const RowVectorBase<EE>& r)
	  { Base::template elementwiseDivideFromLeftInPlace<EE>(r); return *this; }
	template <class EE> inline void 
    elementwiseDivideFromLeft(
        const RowVectorBase<EE>& v, 
        typename RowVectorBase<EE>::template EltResult<ELT>::Dvd& out) const
	{ 
        Base::template elementwiseDivideFromLeft<EE>(v,out);
    }
	template <class EE> inline typename RowVectorBase<EE>::template EltResult<ELT>::Dvd 
    elementwiseDivideFromLeft(const RowVectorBase<EE>& v) const
	{ 
        typename RowVectorBase<EE>::template EltResult<ELT>::Dvd out(nrow()); 
        Base::template elementwiseDivideFromLeft<EE>(v,out); 
        return out;
    }

    // Implicit conversions are allowed to RowVector or Matrix, but not to Vector.   
    operator const RowVector_<ELT>&()     const {return *reinterpret_cast<const RowVector_<ELT>*>(this);}
    operator       RowVector_<ELT>&()           {return *reinterpret_cast<      RowVector_<ELT>*>(this);}
    operator const RowVectorView_<ELT>&() const {return *reinterpret_cast<const RowVectorView_<ELT>*>(this);}
    operator       RowVectorView_<ELT>&()       {return *reinterpret_cast<      RowVectorView_<ELT>*>(this);}
    
    operator const Matrix_<ELT>&()     const {return *reinterpret_cast<const Matrix_<ELT>*>(this);}
    operator       Matrix_<ELT>&()           {return *reinterpret_cast<      Matrix_<ELT>*>(this);} 
    operator const MatrixView_<ELT>&() const {return *reinterpret_cast<const MatrixView_<ELT>*>(this);}
    operator       MatrixView_<ELT>&()       {return *reinterpret_cast<      MatrixView_<ELT>*>(this);} 
    

    // size() for RowVectors is Base::nelt() but returns int instead of ptrdiff_t.
	int size() const { 
		assert(Base::nelt() <= (ptrdiff_t)std::numeric_limits<int>::max()); 
		assert(Base::nrow()==1);
		return (int)Base::nelt();
	}
	int       nrow() const {assert(Base::nrow()==1); return Base::nrow();}
	int       ncol() const {assert(Base::nrow()==1); return Base::ncol();}
	ptrdiff_t nelt() const {assert(Base::nrow()==1); return Base::nelt();}

    // Override MatrixBase operators to return the right shape
    TAbs abs() const {
        TAbs result; Base::abs(result); return result;
    }

    // Override MatrixBase indexing operators          
    const ELT& operator[](int j) const {return *reinterpret_cast<const ELT*>(Base::getHelper().getElt(j));}
    ELT&       operator[](int j)       {return *reinterpret_cast<ELT*>      (Base::updHelper().updElt(j));}
    const ELT& operator()(int j) const {return *reinterpret_cast<const ELT*>(Base::getHelper().getElt(j));}
    ELT&       operator()(int j)       {return *reinterpret_cast<ELT*>      (Base::updHelper().updElt(j));}
         
    // Block (contiguous subvector) creation      
    RowVectorView_<ELT> operator()(int j, int n) const {return Base::operator()(0,j,1,n).getAsRowVectorView();}
    RowVectorView_<ELT> operator()(int j, int n)       {return Base::operator()(0,j,1,n).updAsRowVectorView();}

    // Indexed view creation (arbitrary subvector). Indices must be monotonically increasing.
    RowVectorView_<ELT> index(const Array_<int>& indices) const {
        MatrixHelper<Scalar> h(Base::getHelper().getCharacterCommitment(), Base::getHelper(), indices);
        return RowVectorView_<ELT>(h);
    }
    RowVectorView_<ELT> updIndex(const Array_<int>& indices) {
        MatrixHelper<Scalar> h(Base::getHelper().getCharacterCommitment(), Base::updHelper(), indices);
        return RowVectorView_<ELT>(h);
    }

    RowVectorView_<ELT> operator()(const Array_<int>& indices) const {return index(indices);}
    RowVectorView_<ELT> operator()(const Array_<int>& indices)       {return updIndex(indices);}
 
    // Hermitian transpose.
    THerm transpose() const {return Base::transpose().getAsVectorView();}
    THerm updTranspose()    {return Base::updTranspose().updAsVectorView();}

    THerm operator~() const {return transpose();}
    THerm operator~()       {return updTranspose();}

    const RowVectorBase& operator+() const {return *this; }

    // Negation

    const TNeg& negate()    const {return *reinterpret_cast<const TNeg*>(this); }
    TNeg&       updNegate()       {return *reinterpret_cast<TNeg*>(this); }

    const TNeg& operator-() const {return negate();}
    TNeg&       operator-()       {return updNegate();}

    RowVectorBase& resize(int n)     {Base::resize(1,n); return *this;}
    RowVectorBase& resizeKeep(int n) {Base::resizeKeep(1,n); return *this;}

	//TODO: this is not re-locking the number of rows at 1.
	void clear() {Base::clear(); Base::resize(1,0);}

    ELT sum() const {ELT s; Base::getHelper().sum(reinterpret_cast<Scalar*>(&s)); return s; } // add all the elements        
    VectorIterator<ELT, RowVectorBase<ELT> > begin() {
        return VectorIterator<ELT, RowVectorBase<ELT> >(*this, 0);
    }
    VectorIterator<ELT, RowVectorBase<ELT> > end() {
        return VectorIterator<ELT, RowVectorBase<ELT> >(*this, size());
    }

protected:
    // Create a RowVectorBase handle using a given helper rep. 
    explicit RowVectorBase(MatrixHelperRep<Scalar>* hrep) : Base(hrep) {}

private:
    // NO DATA MEMBERS ALLOWED
};



//  ------------------------------- MatrixView_ --------------------------------
/// This class is identical to a Matrix_; it is used only to manage the C++ rules
/// for when copy constructors are called by introducing a separate type to
/// prevent certain allowed optimizations from occuring when we don't want them.
/// Despite the name, this may be an owner if a Matrix_ is recast to a MatrixView_.
/// However, there are no owner constructors for MatrixView_. 
//  ----------------------------------------------------------------------------
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

    // Copy construction and copy assignment from a DeadMatrixView steals the helper.
    MatrixView_(DeadMatrixView_<ELT>&);
    MatrixView_& operator=(DeadMatrixView_<ELT>&);

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



//  ----------------------------- DeadMatrixView_ ------------------------------
/// This is a MatrixView_ with the additional property that we are about to delete it.
/// If this is the source for an assignment or copy construction, the destination is
/// free to steal the helper and/or the underlying data.
//  ----------------------------------------------------------------------------
template <class ELT> class DeadMatrixView_ : public MatrixView_<ELT> {
    typedef MatrixView_<ELT>                Base;
    typedef typename CNT<ELT>::Scalar       S;
    typedef typename CNT<ELT>::StdNumber    StdNumber;
public:
    // Default construction is suppressed.
    // Uses default destructor.
    
    // All functionality is passed through to MatrixView_.
    explicit DeadMatrixView_(MatrixHelperRep<S>* hrep) : Base(hrep) {}
    DeadMatrixView_(const Base& m) : Base(m) {}
    DeadMatrixView_& operator=(const Base& m) {
        Base::operator=(m); return *this;
    }

    // Ask for shallow copy    
    DeadMatrixView_(const MatrixHelper<S>& h) : Base(h) {}
    DeadMatrixView_(MatrixHelper<S>&       h) : Base(h) {}

    DeadMatrixView_& operator=(const Matrix_<ELT>& v)     { Base::operator=(v); return *this; }
    DeadMatrixView_& operator=(const ELT& e)              { Base::operator=(e); return *this; }

    template <class EE> DeadMatrixView_& operator=(const MatrixBase<EE>& m)
      { Base::operator=(m); return *this; }
    template <class EE> DeadMatrixView_& operator+=(const MatrixBase<EE>& m)
      { Base::operator+=(m); return *this; }
    template <class EE> DeadMatrixView_& operator-=(const MatrixBase<EE>& m)
      { Base::operator-=(m); return *this; }

    DeadMatrixView_& operator*=(const StdNumber& t) { Base::operator*=(t); return *this; }
    DeadMatrixView_& operator/=(const StdNumber& t) { Base::operator/=(t); return *this; }
    DeadMatrixView_& operator+=(const ELT& r)       { this->updDiag() += r; return *this; }
    DeadMatrixView_& operator-=(const ELT& r)       { this->updDiag() -= r; return *this; }  

private:
    // NO DATA MEMBERS ALLOWED
    DeadMatrixView_(); // default constructor suppressed (what's it a view of?)
};

template <class ELT> inline 
MatrixView_<ELT>::MatrixView_(DeadMatrixView_<ELT>& dead) 
:   Base(dead.updHelper().stealRep()) {}

template <class ELT> inline MatrixView_<ELT>& 
MatrixView_<ELT>::operator=(DeadMatrixView_<ELT>& dead) {
    if (Base::getHelper().getCharacterCommitment().isSatisfiedBy(dead.getMatrixCharacter()))
        Base::updHelper().replaceRep(dead.updHelper().stealRep());
    else
        Base::operator=(dead);
    return *this;
}


//  ---------------------------------- Matrix_ ---------------------------------
/// This is the Matrix class intended to appear in user code. It can be a 
/// fixed-size view of someone else's data, or can be a resizable data owner itself.
//  ----------------------------------------------------------------------------
template <class ELT> class Matrix_ : public MatrixBase<ELT> {
    typedef typename CNT<ELT>::Scalar       S;
    typedef typename CNT<ELT>::Number       Number;
    typedef typename CNT<ELT>::StdNumber    StdNumber;

    typedef typename CNT<ELT>::TNeg         ENeg;
    typedef typename CNT<ELT>::THerm        EHerm;

    typedef MatrixBase<ELT>     Base;
    typedef MatrixBase<ENeg>    BaseNeg;
    typedef MatrixBase<EHerm>   BaseHerm;

    typedef Matrix_<ELT>        T;
    typedef MatrixView_<ELT>    TView;
    typedef Matrix_<ENeg>       TNeg;

public:
    Matrix_() : Base() { }
    Matrix_(const MatrixCommitment& mc) : Base(mc) {}

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
   
private:
    // NO DATA MEMBERS ALLOWED
};



//  -------------------------------- VectorView_ -------------------------------
/// This class is identical to a Vector_; it is used only to manage the C++ rules
/// for when copy constructors are called by introducing a separate type to
/// prevent certain allowed optimizations from occuring when we don't want them.
/// Despite the name, this may be an owner if a Vector_ is recast to a VectorView_.
/// However, there are no owner constructors for VectorView_. 
//  ----------------------------------------------------------------------------
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
    VectorView_& operator+=(const ELT& b) { elementwiseAddScalarInPlace(b); return *this; } 
    VectorView_& operator-=(const ELT& b) { elementwiseSubtractScalarInPlace(b); return *this; } 

private:
    // NO DATA MEMBERS ALLOWED
    VectorView_(); // default construction suppressed (what's it a View of?)
};



//  ---------------------------------- Vector_ ---------------------------------
/// This is the Vector class intended to appear in user code. It can be a 
/// fixed-size view of someone else's data, or can be a resizable data owner 
/// itself, although of course it will always have just one column.
//  ----------------------------------------------------------------------------
template <class ELT> class Vector_ : public VectorBase<ELT> {
    typedef typename CNT<ELT>::Scalar       S;
    typedef typename CNT<ELT>::Number       Number;
    typedef typename CNT<ELT>::StdNumber    StdNumber;
    typedef typename CNT<ELT>::TNeg         ENeg;
    typedef VectorBase<ELT>                 Base;
    typedef VectorBase<ENeg>                BaseNeg;
public:
    Vector_() : Base() {}  // 0x1 reallocatable
    // Uses default destructor.

    // Copy constructor is deep.
    Vector_(const Vector_& src) : Base(src) {}

    // Implicit conversions.
    Vector_(const Base& src) : Base(src) {}    // e.g., VectorView
    Vector_(const BaseNeg& src) : Base(src) {}

    // Copy assignment is deep and can be reallocating if this Vector
    // has no View.
    Vector_& operator=(const Vector_& src) {
        Base::operator=(src); return*this;
    }


    explicit Vector_(int m) : Base(m) { }
    Vector_(int m, const ELT* cppInitialValues) : Base(m, cppInitialValues) {}
    Vector_(int m, const ELT& initialValue)     : Base(m, initialValue) {}

    /// Construct a Vector which uses borrowed space with assumed
    /// element-to-element stride equal to the C++ element spacing.
    /// Last parameter is a dummy to avoid overload conflicts when ELT=S;
    /// pass it as "true".
    Vector_(int m, const S* cppData, bool): Base(m, Base::CppNScalarsPerElement, cppData) {}
    Vector_(int m,       S* cppData, bool): Base(m, Base::CppNScalarsPerElement, cppData) {}

    /// Borrowed-space construction with explicit stride supplied as
    /// "number of scalars between elements". Last parameter is a 
    /// dummy to avoid overload conflicts; pass it as "true".
    Vector_(int m, int stride, const S* data, bool) : Base(m, stride, data) {}
    Vector_(int m, int stride,       S* data, bool) : Base(m, stride, data) {}

    /// Convert a Vec to a Vector_.
    template <int M>
    explicit Vector_(const Vec<M,ELT>& v) : Base(M) {
        for (int i = 0; i < M; ++i)
            this->updElt(i, 0) = v(i);
    }

    Vector_& operator=(const ELT& v) { Base::operator=(v); return *this; } 

    template <class EE> Vector_& operator=(const VectorBase<EE>& m)
      { Base::operator=(m); return*this; }
    template <class EE> Vector_& operator+=(const VectorBase<EE>& m)
      { Base::operator+=(m); return*this; }
    template <class EE> Vector_& operator-=(const VectorBase<EE>& m)
      { Base::operator-=(m); return*this; }

    Vector_& operator*=(const StdNumber& t) { Base::operator*=(t); return *this; }
    Vector_& operator/=(const StdNumber& t) { Base::operator/=(t); return *this; }
    Vector_& operator+=(const ELT& b) { elementwiseAddScalarInPlace(b); return *this; } 
    Vector_& operator-=(const ELT& b) { elementwiseSubtractScalarInPlace(b); return *this; } 
 
private:
    // NO DATA MEMBERS ALLOWED
};



//  ------------------------------ RowVectorView_ ------------------------------
/// This class is identical to a RowVector_; it is used only to manage the C++ 
/// rules for when copy constructors are called by introducing a separate type to
/// prevent certain allowed optimizations from occuring when we don't want them.
/// Despite the name, this may be an owner if a RowVector_ is recast to a 
/// RowVectorView_. However, there are no owner constructors for RowVectorView_. 
//  ----------------------------------------------------------------------------
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
    RowVectorView_& operator+=(const ELT& b) { elementwiseAddScalarInPlace(b); return *this; } 
    RowVectorView_& operator-=(const ELT& b) { elementwiseSubtractScalarInPlace(b); return *this; } 

private:
    // NO DATA MEMBERS ALLOWED
    RowVectorView_(); // default construction suppressed (what is it a view of?)
};



//  -------------------------------- RowVector_ --------------------------------
/// RowVectors are much less common than Vectors. However, if a Simmatrix user 
/// wants one, this is the class intended to appear in user code. It can be a 
/// fixed-size view of someone else's data, or can be a resizable data owner 
/// itself, although of course it will always have just one row.
//  ----------------------------------------------------------------------------
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
    RowVector_& operator+=(const ELT& b) { elementwiseAddScalarInPlace(b); return *this; } 
    RowVector_& operator-=(const ELT& b) { elementwiseSubtractScalarInPlace(b); return *this; } 

private:
    // NO DATA MEMBERS ALLOWED
};



//  ------------------------ MatrixBase definitions ----------------------------

template <class ELT> inline MatrixView_<ELT> 
MatrixBase<ELT>::block(int i, int j, int m, int n) const { 
    SimTK_INDEXCHECK(i,nrow()+1,"MatrixBase::block()");
    SimTK_INDEXCHECK(j,ncol()+1,"MatrixBase::block()");
    SimTK_SIZECHECK(i+m,nrow(),"MatrixBase::block()");
    SimTK_SIZECHECK(j+n,ncol(),"MatrixBase::block()");

    MatrixHelper<Scalar> h(MatrixCommitment(),helper,i,j,m,n);    
    return MatrixView_<ELT>(h.stealRep()); 
}
    
template <class ELT> inline MatrixView_<ELT>
MatrixBase<ELT>::updBlock(int i, int j, int m, int n) { 
    SimTK_INDEXCHECK(i,nrow()+1,"MatrixBase::updBlock()");
    SimTK_INDEXCHECK(j,ncol()+1,"MatrixBase::updBlock()");
    SimTK_SIZECHECK(i+m,nrow(),"MatrixBase::updBlock()");
    SimTK_SIZECHECK(j+n,ncol(),"MatrixBase::updBlock()");

    MatrixHelper<Scalar> h(MatrixCommitment(),helper,i,j,m,n);        
    return MatrixView_<ELT>(h.stealRep()); 
}

template <class E> inline MatrixView_<typename CNT<E>::THerm>
MatrixBase<E>::transpose() const { 
    MatrixHelper<typename CNT<Scalar>::THerm> 
        h(MatrixCommitment(),
          helper, typename MatrixHelper<typename CNT<Scalar>::THerm>::TransposeView());
    return MatrixView_<typename CNT<E>::THerm>(h.stealRep()); 
}
    
template <class E> inline MatrixView_<typename CNT<E>::THerm>
MatrixBase<E>::updTranspose() {     
    MatrixHelper<typename CNT<Scalar>::THerm> 
        h(MatrixCommitment(),
          helper, typename MatrixHelper<typename CNT<Scalar>::THerm>::TransposeView());
    return MatrixView_<typename CNT<E>::THerm>(h.stealRep()); 
}

template <class E> inline VectorView_<E>
MatrixBase<E>::diag() const { 
    MatrixHelper<Scalar> h(MatrixCommitment::Vector(),
                           helper, typename MatrixHelper<Scalar>::DiagonalView());
    return VectorView_<E>(h.stealRep()); 
}
    
template <class E> inline VectorView_<E>
MatrixBase<E>::updDiag() {     
    MatrixHelper<Scalar> h(MatrixCommitment::Vector(),
                           helper, typename MatrixHelper<Scalar>::DiagonalView());
    return VectorView_<E>(h.stealRep());
}

template <class ELT> inline VectorView_<ELT> 
MatrixBase<ELT>::col(int j) const { 
    SimTK_INDEXCHECK(j,ncol(),"MatrixBase::col()");

    MatrixHelper<Scalar> h(MatrixCommitment::Vector(),
                           helper,0,j,nrow(),1);    
    return VectorView_<ELT>(h.stealRep()); 
}
    
template <class ELT> inline VectorView_<ELT>
MatrixBase<ELT>::updCol(int j) {
    SimTK_INDEXCHECK(j,ncol(),"MatrixBase::updCol()");

    MatrixHelper<Scalar> h(MatrixCommitment::Vector(),
                           helper,0,j,nrow(),1);        
    return VectorView_<ELT>(h.stealRep()); 
}

template <class ELT> inline RowVectorView_<ELT> 
MatrixBase<ELT>::row(int i) const { 
    SimTK_INDEXCHECK(i,nrow(),"MatrixBase::row()");

    MatrixHelper<Scalar> h(MatrixCommitment::RowVector(),
                           helper,i,0,1,ncol());    
    return RowVectorView_<ELT>(h.stealRep()); 
}
    
template <class ELT> inline RowVectorView_<ELT>
MatrixBase<ELT>::updRow(int i) { 
    SimTK_INDEXCHECK(i,nrow(),"MatrixBase::updRow()");

    MatrixHelper<Scalar> h(MatrixCommitment::RowVector(),
                           helper,i,0,1,ncol());        
    return RowVectorView_<ELT>(h.stealRep()); 
}

// M = diag(v) * M; v must have nrow() elements.
// That is, M[i] *= v[i].
template <class ELT> template <class EE> inline MatrixBase<ELT>& 
MatrixBase<ELT>::rowScaleInPlace(const VectorBase<EE>& v) {
	assert(v.nrow() == nrow());
	for (int i=0; i < nrow(); ++i)
		(*this)[i] *= v[i];
	return *this;
}

template <class ELT> template <class EE> inline void
MatrixBase<ELT>::rowScale(const VectorBase<EE>& v, typename MatrixBase<ELT>::template EltResult<EE>::Mul& out) const {
	assert(v.nrow() == nrow());
	out.resize(nrow(), ncol());
    for (int j=0; j<ncol(); ++j)
	    for (int i=0; i<nrow(); ++i)
		   out(i,j) = (*this)(i,j) * v[i];
}

// M = M * diag(v); v must have ncol() elements
// That is, M(i) *= v[i]
template <class ELT> template <class EE>  inline MatrixBase<ELT>& 
MatrixBase<ELT>::colScaleInPlace(const VectorBase<EE>& v) {
    assert(v.nrow() == ncol());
	for (int j=0; j < ncol(); ++j)
		(*this)(j) *= v[j];
	return *this;
}

template <class ELT> template <class EE> inline void
MatrixBase<ELT>::colScale(const VectorBase<EE>& v, typename MatrixBase<ELT>::template EltResult<EE>::Mul& out) const {
	assert(v.nrow() == ncol());
	out.resize(nrow(), ncol());
    for (int j=0; j<ncol(); ++j)
	    for (int i=0; i<nrow(); ++i)
		   out(i,j) = (*this)(i,j) * v[j];
}


// M(i,j) *= r[i]*c[j]; r must have nrow() elements; c must have ncol() elements
template <class ELT> template <class ER, class EC> inline MatrixBase<ELT>& 
MatrixBase<ELT>::rowAndColScaleInPlace(const VectorBase<ER>& r, const VectorBase<EC>& c) {
	assert(r.nrow()==nrow() && c.nrow()==ncol());
    for (int j=0; j<ncol(); ++j)
	    for (int i=0; i<nrow(); ++i)
			(*this)(i,j) *= (r[i]*c[j]);
    return *this;
}

template <class ELT> template <class ER, class EC> inline void
MatrixBase<ELT>::rowAndColScale(
    const VectorBase<ER>& r, 
    const VectorBase<EC>& c,
    typename EltResult<typename VectorBase<ER>::template EltResult<EC>::Mul>::Mul& 
                          out) const
{
	assert(r.nrow()==nrow() && c.nrow()==ncol());
    out.resize(nrow(), ncol());
    for (int j=0; j<ncol(); ++j)
	    for (int i=0; i<nrow(); ++i)
			out(i,j) = (*this)(i,j) * (r[i]*c[j]);
}


// Set M(i,j) = M(i,j)^-1.
template <class ELT> inline MatrixBase<ELT>& 
MatrixBase<ELT>::elementwiseInvertInPlace() {
    const int nr=nrow(), nc=ncol();
    for (int j=0; j<nc; ++j)
        for (int i=0; i<nr; ++i) {
            ELT& e = updElt(i,j);
            e = CNT<ELT>::invert(e);
        }
    return *this;
}

template <class ELT> inline void
MatrixBase<ELT>::elementwiseInvert(MatrixBase<typename CNT<E>::TInvert>& out) const {
    const int nr=nrow(), nc=ncol();
    out.resize(nr,nc);
    for (int j=0; j<nc; ++j)
        for (int i=0; i<nr; ++i)
            out(i,j) = CNT<ELT>::invert((*this)(i,j));
}

// M(i,j) += s
template <class ELT> template <class S> inline MatrixBase<ELT>& 
MatrixBase<ELT>::elementwiseAddScalarInPlace(const S& s) {
    for (int j=0; j<ncol(); ++j)
	    for (int i=0; i<nrow(); ++i)
			(*this)(i,j) += s;
    return *this;
}

template <class ELT> template <class S> inline void 
MatrixBase<ELT>::elementwiseAddScalar(
    const S& s,
    typename MatrixBase<ELT>::template EltResult<S>::Add& out) const
{
    const int nr=nrow(), nc=ncol();
    out.resize(nr,nc);
    for (int j=0; j<nc; ++j)
	    for (int i=0; i<nr; ++i)
			out(i,j) = (*this)(i,j) + s;
}

// M(i,j) -= s
template <class ELT> template <class S> inline MatrixBase<ELT>& 
MatrixBase<ELT>::elementwiseSubtractScalarInPlace(const S& s) {
    for (int j=0; j<ncol(); ++j)
	    for (int i=0; i<nrow(); ++i)
			(*this)(i,j) -= s;
    return *this;
}

template <class ELT> template <class S> inline void 
MatrixBase<ELT>::elementwiseSubtractScalar(
    const S& s,
    typename MatrixBase<ELT>::template EltResult<S>::Sub& out) const
{
    const int nr=nrow(), nc=ncol();
    out.resize(nr,nc);
    for (int j=0; j<nc; ++j)
	    for (int i=0; i<nr; ++i)
			out(i,j) = (*this)(i,j) - s;
}

// M(i,j) = s - M(i,j)
template <class ELT> template <class S> inline MatrixBase<ELT>& 
MatrixBase<ELT>::elementwiseSubtractFromScalarInPlace(const S& s) {
    const int nr=nrow(), nc=ncol();
    for (int j=0; j<nc; ++j)
        for (int i=0; i<nr; ++i) {
            ELT& e = updElt(i,j);
			e = s - e;
        }
    return *this;
}

template <class ELT> template <class S> inline void 
MatrixBase<ELT>::elementwiseSubtractFromScalar(
    const S& s,
    typename MatrixBase<S>::template EltResult<ELT>::Sub& out) const
{
    const int nr=nrow(), nc=ncol();
    out.resize(nr,nc);
    for (int j=0; j<nc; ++j)
	    for (int i=0; i<nr; ++i)
			out(i,j) = s - (*this)(i,j);
}

// M(i,j) *= R(i,j); R must have same dimensions as this
template <class ELT> template <class EE> inline MatrixBase<ELT>& 
MatrixBase<ELT>::elementwiseMultiplyInPlace(const MatrixBase<EE>& r) {
    const int nr=nrow(), nc=ncol();
	assert(r.nrow()==nr && r.ncol()==nc);
    for (int j=0; j<nc; ++j)
	    for (int i=0; i<nr; ++i)
			(*this)(i,j) *= r(i,j);
    return *this;
}

template <class ELT> template <class EE> inline void 
MatrixBase<ELT>::elementwiseMultiply(
    const MatrixBase<EE>& r, 
    typename MatrixBase<ELT>::template EltResult<EE>::Mul& out) const
{
    const int nr=nrow(), nc=ncol();
	assert(r.nrow()==nr && r.ncol()==nc);
	out.resize(nr,nc);
    for (int j=0; j<nc; ++j)
	    for (int i=0; i<nr; ++i)
			out(i,j) = (*this)(i,j) * r(i,j);
}

// M(i,j) = R(i,j) * M(i,j); R must have same dimensions as this
template <class ELT> template <class EE> inline MatrixBase<ELT>& 
MatrixBase<ELT>::elementwiseMultiplyFromLeftInPlace(const MatrixBase<EE>& r) {
    const int nr=nrow(), nc=ncol();
	assert(r.nrow()==nr && r.ncol()==nc);
    for (int j=0; j<nc; ++j)
        for (int i=0; i<nr; ++i) {
            ELT& e = updElt(i,j);
			e = r(i,j) * e;
        }
    return *this;
}

template <class ELT> template <class EE> inline void 
MatrixBase<ELT>::elementwiseMultiplyFromLeft(
    const MatrixBase<EE>& r, 
    typename MatrixBase<EE>::template EltResult<ELT>::Mul& out) const
{
    const int nr=nrow(), nc=ncol();
	assert(r.nrow()==nr && r.ncol()==nc);
	out.resize(nr,nc);
    for (int j=0; j<nc; ++j)
	    for (int i=0; i<nr; ++i)
			out(i,j) =  r(i,j) * (*this)(i,j);
}

// M(i,j) /= R(i,j); R must have same dimensions as this
template <class ELT> template <class EE> inline MatrixBase<ELT>& 
MatrixBase<ELT>::elementwiseDivideInPlace(const MatrixBase<EE>& r) {
    const int nr=nrow(), nc=ncol();
	assert(r.nrow()==nr && r.ncol()==nc);
    for (int j=0; j<nc; ++j)
	    for (int i=0; i<nr; ++i)
			(*this)(i,j) /= r(i,j);
    return *this;
}

template <class ELT> template <class EE> inline void 
MatrixBase<ELT>::elementwiseDivide(
    const MatrixBase<EE>& r,
    typename MatrixBase<ELT>::template EltResult<EE>::Dvd& out) const
{
    const int nr=nrow(), nc=ncol();
	assert(r.nrow()==nr && r.ncol()==nc);
	out.resize(nr,nc);
    for (int j=0; j<nc; ++j)
	    for (int i=0; i<nr; ++i)
			out(i,j) = (*this)(i,j) / r(i,j);
}
// M(i,j) = R(i,j) / M(i,j); R must have same dimensions as this
template <class ELT> template <class EE> inline MatrixBase<ELT>& 
MatrixBase<ELT>::elementwiseDivideFromLeftInPlace(const MatrixBase<EE>& r) {
    const int nr=nrow(), nc=ncol();
	assert(r.nrow()==nr && r.ncol()==nc);
    for (int j=0; j<nc; ++j)
        for (int i=0; i<nr; ++i) {
            ELT& e = updElt(i,j);
			e = r(i,j) / e;
        }
    return *this;
}

template <class ELT> template <class EE> inline void 
MatrixBase<ELT>::elementwiseDivideFromLeft(
    const MatrixBase<EE>& r,
    typename MatrixBase<EE>::template EltResult<ELT>::Dvd& out) const
{
    const int nr=nrow(), nc=ncol();
	assert(r.nrow()==nr && r.ncol()==nc);
	out.resize(nr,nc);
    for (int j=0; j<nc; ++j)
	    for (int i=0; i<nr; ++i)
			out(i,j) = r(i,j) / (*this)(i,j);
}

/*
template <class ELT> inline MatrixView_< typename CNT<ELT>::TReal > 
MatrixBase<ELT>::real() const { 
    if (!CNT<ELT>::IsComplex) { // known at compile time
        return MatrixView_< typename CNT<ELT>::TReal >( // this is just ELT
            MatrixHelper(helper,0,0,nrow(),ncol()));    // a view of the whole matrix
    }
    // Elements are complex -- helper uses underlying precision (real) type.
    MatrixHelper<Precision> h(helper,typename MatrixHelper<Precision>::RealView);    
    return MatrixView_< typename CNT<ELT>::TReal >(h); 
}
*/


//  ----------------------------------------------------------------------------
/// @name Global operators involving Matrix objects
/// These operators take MatrixBase arguments and produce Matrix_ results.
/// @{

// + and - allow mixed element types, but will fail to compile if the elements aren't
// compatible. At run time these will fail if the dimensions are incompatible.
template <class E1, class E2>
Matrix_<typename CNT<E1>::template Result<E2>::Add>
operator+(const MatrixBase<E1>& l, const MatrixBase<E2>& r) {
    return Matrix_<typename CNT<E1>::template Result<E2>::Add>(l) += r;
}

template <class E>
Matrix_<E> operator+(const MatrixBase<E>& l, const typename CNT<E>::T& r) {
    return Matrix_<E>(l) += r;
}

template <class E>
Matrix_<E> operator+(const typename CNT<E>::T& l, const MatrixBase<E>& r) {
    return Matrix_<E>(r) += l;
}

template <class E1, class E2>
Matrix_<typename CNT<E1>::template Result<E2>::Sub>
operator-(const MatrixBase<E1>& l, const MatrixBase<E2>& r) {
    return Matrix_<typename CNT<E1>::template Result<E2>::Sub>(l) -= r;
}

template <class E>
Matrix_<E> operator-(const MatrixBase<E>& l, const typename CNT<E>::T& r) {
    return Matrix_<E>(l) -= r;
}

template <class E>
Matrix_<E> operator-(const typename CNT<E>::T& l, const MatrixBase<E>& r) {
    Matrix_<E> temp(r.nrow(), r.ncol());
    temp = l;
    return (temp -= r);
}

// Scalar multiply and divide. You might wish the scalar could be
// a templatized type "E2", but that would create horrible ambiguities since
// E2 would match not only scalar types but everything else including
// matrices.
template <class E> Matrix_<E>
operator*(const MatrixBase<E>& l, const typename CNT<E>::StdNumber& r) 
  { return Matrix_<E>(l)*=r; }

template <class E> Matrix_<E>
operator*(const typename CNT<E>::StdNumber& l, const MatrixBase<E>& r) 
  { return Matrix_<E>(r)*=l; }

template <class E> Matrix_<E>
operator/(const MatrixBase<E>& l, const typename CNT<E>::StdNumber& r) 
  { return Matrix_<E>(l)/=r; }

// Handle ints explicitly.
template <class E> Matrix_<E>
operator*(const MatrixBase<E>& l, int r) 
  { return Matrix_<E>(l)*= typename CNT<E>::StdNumber(r); }

template <class E> Matrix_<E>
operator*(int l, const MatrixBase<E>& r) 
  { return Matrix_<E>(r)*= typename CNT<E>::StdNumber(l); }

template <class E> Matrix_<E>
operator/(const MatrixBase<E>& l, int r) 
  { return Matrix_<E>(l)/= typename CNT<E>::StdNumber(r); }

/// @}

/// @name Global operators involving Vector objects
/// These operators take VectorBase arguments and produce Vector_ results.
/// @{

template <class E1, class E2>
Vector_<typename CNT<E1>::template Result<E2>::Add>
operator+(const VectorBase<E1>& l, const VectorBase<E2>& r) {
    return Vector_<typename CNT<E1>::template Result<E2>::Add>(l) += r;
}
template <class E>
Vector_<E> operator+(const VectorBase<E>& l, const typename CNT<E>::T& r) {
    return Vector_<E>(l) += r;
}
template <class E>
Vector_<E> operator+(const typename CNT<E>::T& l, const VectorBase<E>& r) {
    return Vector_<E>(r) += l;
}
template <class E1, class E2>
Vector_<typename CNT<E1>::template Result<E2>::Sub>
operator-(const VectorBase<E1>& l, const VectorBase<E2>& r) {
    return Vector_<typename CNT<E1>::template Result<E2>::Sub>(l) -= r;
}
template <class E>
Vector_<E> operator-(const VectorBase<E>& l, const typename CNT<E>::T& r) {
    return Vector_<E>(l) -= r;
}
template <class E>
Vector_<E> operator-(const typename CNT<E>::T& l, const VectorBase<E>& r) {
    Vector_<E> temp(r.size());
    temp = l;
    return (temp -= r);
}

// Scalar multiply and divide.

template <class E> Vector_<E>
operator*(const VectorBase<E>& l, const typename CNT<E>::StdNumber& r) 
  { return Vector_<E>(l)*=r; }

template <class E> Vector_<E>
operator*(const typename CNT<E>::StdNumber& l, const VectorBase<E>& r) 
  { return Vector_<E>(r)*=l; }

template <class E> Vector_<E>
operator/(const VectorBase<E>& l, const typename CNT<E>::StdNumber& r) 
  { return Vector_<E>(l)/=r; }

// Handle ints explicitly
template <class E> Vector_<E>
operator*(const VectorBase<E>& l, int r) 
  { return Vector_<E>(l)*= typename CNT<E>::StdNumber(r); }

template <class E> Vector_<E>
operator*(int l, const VectorBase<E>& r) 
  { return Vector_<E>(r)*= typename CNT<E>::StdNumber(l); }

template <class E> Vector_<E>
operator/(const VectorBase<E>& l, int r) 
  { return Vector_<E>(l)/= typename CNT<E>::StdNumber(r); }

// These are fancier "scalars"; whether they are allowed depends on
// whether the element type and the CNT are compatible.

// Vector * Vec
template <class E1, int M, class E2, int S> 
Vector_<typename CNT<E1>::template Result< Vec<M,E2,S> >::Mul>
operator*(const VectorBase<E1>& v, const Vec<M,E2,S>& s) {
    Vector_<typename CNT<E1>::template Result< Vec<M,E2,S> >::Mul> res(v.nrow());
    for (int i=0; i < v.nrow(); ++i)
        res[i] = v[i]*s; 
    return res;
}

// Vec * Vector
template <class E1, int M, class E2, int S> 
Vector_<typename CNT< Vec<M,E2,S> >::template Result<E1>::Mul>
operator*(const Vec<M,E2,S>& s, const VectorBase<E1>& v) {
    Vector_<typename CNT< Vec<M,E2,S> >::template Result<E1>::Mul> res(v.nrow());
    for (int i=0; i < v.nrow(); ++i)
        res[i] = s*v[i]; 
    return res;
}

// Vector * Row
template <class E1, int N, class E2, int S> 
Vector_<typename CNT<E1>::template Result< Row<N,E2,S> >::Mul>
operator*(const VectorBase<E1>& v, const Row<N,E2,S>& s) {
    Vector_<typename CNT<E1>::template Result< Row<N,E2,S> >::Mul> res(v.nrow());
    for (int i=0; i < v.nrow(); ++i)
        res[i] = v[i]*s; 
    return res;
}

// Row * Vector
template <class E1, int N, class E2, int S> 
Vector_<typename CNT< Row<N,E2,S> >::template Result<E1>::Mul>
operator*(const Row<N,E2,S>& s, const VectorBase<E1>& v) {
    Vector_<typename CNT< Row<N,E2,S> >::template Result<E1>::Mul> res(v.nrow());
    for (int i=0; i < v.nrow(); ++i)
        res[i] = s*v[i]; 
    return res;
}

// Vector * Mat
template <class E1, int M, int N, class E2, int S1, int S2> 
Vector_<typename CNT<E1>::template Result< Mat<M,N,E2,S1,S2> >::Mul>
operator*(const VectorBase<E1>& v, const Mat<M,N,E2,S1,S2>& s) {
    Vector_<typename CNT<E1>::template Result< Mat<M,N,E2,S1,S2> >::Mul> res(v.nrow());
    for (int i=0; i < v.nrow(); ++i)
        res[i] = v[i]*s; 
    return res;
}

// Mat * Vector
template <class E1, int M, int N, class E2, int S1, int S2> 
Vector_<typename CNT< Mat<M,N,E2,S1,S2> >::template Result<E1>::Mul>
operator*(const Mat<M,N,E2,S1,S2>& s, const VectorBase<E1>& v) {
    Vector_<typename CNT< Mat<M,N,E2,S1,S2> >::template Result<E1>::Mul> res(v.nrow());
    for (int i=0; i < v.nrow(); ++i)
        res[i] = s*v[i]; 
    return res;
}

// Vector * SymMat
template <class E1, int M, class E2, int S> 
Vector_<typename CNT<E1>::template Result< SymMat<M,E2,S> >::Mul>
operator*(const VectorBase<E1>& v, const SymMat<M,E2,S>& s) {
    Vector_<typename CNT<E1>::template Result< SymMat<M,E2,S> >::Mul> res(v.nrow());
    for (int i=0; i < v.nrow(); ++i)
        res[i] = v[i]*s; 
    return res;
}

// SymMat * Vector
template <class E1, int M, class E2, int S> 
Vector_<typename CNT< SymMat<M,E2,S> >::template Result<E1>::Mul>
operator*(const SymMat<M,E2,S>& s, const VectorBase<E1>& v) {
    Vector_<typename CNT< SymMat<M,E2,S> >::template Result<E1>::Mul> res(v.nrow());
    for (int i=0; i < v.nrow(); ++i)
        res[i] = s*v[i]; 
    return res;
}

/// @}

/// @name Global operators involving RowVector objects
/// These operators take RowVectorBase arguments and produce RowVector_ results.
/// @{

template <class E1, class E2>
RowVector_<typename CNT<E1>::template Result<E2>::Add>
operator+(const RowVectorBase<E1>& l, const RowVectorBase<E2>& r) {
    return RowVector_<typename CNT<E1>::template Result<E2>::Add>(l) += r;
}
template <class E>
RowVector_<E> operator+(const RowVectorBase<E>& l, const typename CNT<E>::T& r) {
    return RowVector_<E>(l) += r;
}
template <class E>
RowVector_<E> operator+(const typename CNT<E>::T& l, const RowVectorBase<E>& r) {
    return RowVector_<E>(r) += l;
}
template <class E1, class E2>
RowVector_<typename CNT<E1>::template Result<E2>::Sub>
operator-(const RowVectorBase<E1>& l, const RowVectorBase<E2>& r) {
    return RowVector_<typename CNT<E1>::template Result<E2>::Sub>(l) -= r;
}
template <class E>
RowVector_<E> operator-(const RowVectorBase<E>& l, const typename CNT<E>::T& r) {
    return RowVector_<E>(l) -= r;
}
template <class E>
RowVector_<E> operator-(const typename CNT<E>::T& l, const RowVectorBase<E>& r) {
    RowVector_<E> temp(r.size());
    temp = l;
    return (temp -= r);
}

// Scalar multiply and divide 

template <class E> RowVector_<E>
operator*(const RowVectorBase<E>& l, const typename CNT<E>::StdNumber& r) 
  { return RowVector_<E>(l)*=r; }

template <class E> RowVector_<E>
operator*(const typename CNT<E>::StdNumber& l, const RowVectorBase<E>& r) 
  { return RowVector_<E>(r)*=l; }

template <class E> RowVector_<E>
operator/(const RowVectorBase<E>& l, const typename CNT<E>::StdNumber& r) 
  { return RowVector_<E>(l)/=r; }

// Handle ints explicitly.
template <class E> RowVector_<E>
operator*(const RowVectorBase<E>& l, int r) 
  { return RowVector_<E>(l)*= typename CNT<E>::StdNumber(r); }

template <class E> RowVector_<E>
operator*(int l, const RowVectorBase<E>& r) 
  { return RowVector_<E>(r)*= typename CNT<E>::StdNumber(l); }

template <class E> RowVector_<E>
operator/(const RowVectorBase<E>& l, int r) 
  { return RowVector_<E>(l)/= typename CNT<E>::StdNumber(r); }


// These are fancier "scalars"; whether they are allowed depends on
// whether the element type and the CNT are compatible.

// RowVector * Vec
template <class E1, int M, class E2, int S> 
RowVector_<typename CNT<E1>::template Result< Vec<M,E2,S> >::Mul>
operator*(const RowVectorBase<E1>& v, const Vec<M,E2,S>& s) {
    RowVector_<typename CNT<E1>::template Result< Vec<M,E2,S> >::Mul> res(v.ncol());
    for (int i=0; i < v.ncol(); ++i)
        res[i] = v[i]*s; 
    return res;
}

// Vec * RowVector
template <class E1, int M, class E2, int S> 
RowVector_<typename CNT< Vec<M,E2,S> >::template Result<E1>::Mul>
operator*(const Vec<M,E2,S>& s, const RowVectorBase<E1>& v) {
    RowVector_<typename CNT< Vec<M,E2,S> >::template Result<E1>::Mul> res(v.ncol());
    for (int i=0; i < v.ncol(); ++i)
        res[i] = s*v[i]; 
    return res;
}

// RowVector * Row
template <class E1, int N, class E2, int S> 
RowVector_<typename CNT<E1>::template Result< Row<N,E2,S> >::Mul>
operator*(const RowVectorBase<E1>& v, const Row<N,E2,S>& s) {
    RowVector_<typename CNT<E1>::template Result< Row<N,E2,S> >::Mul> res(v.ncol());
    for (int i=0; i < v.ncol(); ++i)
        res[i] = v[i]*s; 
    return res;
}

// Row * RowVector
template <class E1, int N, class E2, int S> 
RowVector_<typename CNT< Row<N,E2,S> >::template Result<E1>::Mul>
operator*(const Row<N,E2,S>& s, const RowVectorBase<E1>& v) {
    RowVector_<typename CNT< Row<N,E2,S> >::template Result<E1>::Mul> res(v.ncol());
    for (int i=0; i < v.ncol(); ++i)
        res[i] = s*v[i]; 
    return res;
}

// RowVector * Mat
template <class E1, int M, int N, class E2, int S1, int S2> 
RowVector_<typename CNT<E1>::template Result< Mat<M,N,E2,S1,S2> >::Mul>
operator*(const RowVectorBase<E1>& v, const Mat<M,N,E2,S1,S2>& s) {
    RowVector_<typename CNT<E1>::template Result< Mat<M,N,E2,S1,S2> >::Mul> res(v.ncol());
    for (int i=0; i < v.ncol(); ++i)
        res[i] = v[i]*s; 
    return res;
}

// Mat * RowVector
template <class E1, int M, int N, class E2, int S1, int S2> 
RowVector_<typename CNT< Mat<M,N,E2,S1,S2> >::template Result<E1>::Mul>
operator*(const Mat<M,N,E2,S1,S2>& s, const RowVectorBase<E1>& v) {
    RowVector_<typename CNT< Mat<M,N,E2,S1,S2> >::template Result<E1>::Mul> res(v.ncol());
    for (int i=0; i < v.ncol(); ++i)
        res[i] = s*v[i]; 
    return res;
}

// RowVector * SymMat
template <class E1, int M, class E2, int S> 
RowVector_<typename CNT<E1>::template Result< SymMat<M,E2,S> >::Mul>
operator*(const RowVectorBase<E1>& v, const SymMat<M,E2,S>& s) {
    RowVector_<typename CNT<E1>::template Result< SymMat<M,E2,S> >::Mul> res(v.ncol());
    for (int i=0; i < v.ncol(); ++i)
        res[i] = v[i]*s; 
    return res;
}

// SymMat * RowVector
template <class E1, int M, class E2, int S> 
RowVector_<typename CNT< SymMat<M,E2,S> >::template Result<E1>::Mul>
operator*(const SymMat<M,E2,S>& s, const RowVectorBase<E1>& v) {
    RowVector_<typename CNT< SymMat<M,E2,S> >::template Result<E1>::Mul> res(v.ncol());
    for (int i=0; i < v.ncol(); ++i)
        res[i] = s*v[i]; 
    return res;
}

/// @}


/// @name Global operators involving mixed matrix, vector, and row vector objects
/// These operators take MatrixBase, VectorBase, and RowVectorBase arguments
/// and produce Matrix_, Vector_, and RowVector_ results.
/// @{

    // TODO: these should use LAPACK!

// Dot product
template <class E1, class E2> 
typename CNT<E1>::template Result<E2>::Mul
operator*(const RowVectorBase<E1>& r, const VectorBase<E2>& v) {
    assert(r.ncol() == v.nrow());
    typename CNT<E1>::template Result<E2>::Mul sum(0);
    for (int j=0; j < r.ncol(); ++j)
        sum += r(j) * v[j];
    return sum;
}

template <class E1, class E2> 
Vector_<typename CNT<E1>::template Result<E2>::Mul>
operator*(const MatrixBase<E1>& m, const VectorBase<E2>& v) {
    assert(m.ncol() == v.nrow());
    Vector_<typename CNT<E1>::template Result<E2>::Mul> res(m.nrow());
    for (int i=0; i< m.nrow(); ++i)
        res[i] = m[i]*v;
    return res;
}

template <class E1, class E2> 
Matrix_<typename CNT<E1>::template Result<E2>::Mul>
operator*(const MatrixBase<E1>& m1, const MatrixBase<E2>& m2) {
    assert(m1.ncol() == m2.nrow());
    Matrix_<typename CNT<E1>::template Result<E2>::Mul> 
        res(m1.nrow(),m2.ncol());

    for (int j=0; j < res.ncol(); ++j)
        for (int i=0; i < res.nrow(); ++i)
            res(i,j) = m1[i] * m2(j);

    return res;
}

/// @}

// This "private" static method is used to implement VectorView's 
// fillVectorViewFromStream() and Vector's readVectorFromStream() 
// namespace-scope static methods, which are in turn used to implement 
// VectorView's and 
// Vector's stream extraction operators ">>". This method has to be in the 
// header file so that we don't need to pass streams through the API, but it 
// is not intended for use by users and has no Doxygen presence, unlike 
// fillArrayFromStream() and readArrayFromStream() and (more commonly)
// the extraction operators.
template <class T> static inline 
std::istream& readVectorFromStreamHelper
   (std::istream& in, bool isFixedSize, Vector_<T>& out)
{
    // If already failed, bad, or eof, set failed bit and return without 
    // touching the Vector.
    if (!in.good()) {in.setstate(std::ios::failbit); return in;}

    // If the passed-in Vector isn't resizeable, then we have to treat it as
    // a fixed size VectorView regardless of the setting of the isFixedSize
    // argument.
    if (!out.isResizeable())
        isFixedSize = true; // might be overriding the argument here

    // numRequired will be ignored unless isFixedSize==true.
    const int numRequired = isFixedSize ? out.size() : 0;

    if (!isFixedSize)
        out.clear(); // We're going to replace the entire contents of the Array.

    // Skip initial whitespace. If that results in eof this may be a successful
    // read of a 0-length, unbracketed Vector. That is OK for either a
    // variable-length Vector or a fixed-length VectorView of length zero.
    std::ws(in); if (in.fail()) return in;
    if (in.eof()) {
        if (isFixedSize && numRequired != 0)
            in.setstate(std::ios_base::failbit); // zero elements not OK
        return in;
    }
    
    // Here the stream is good and the next character is non-white.
    assert(in.good());

    // Use this for raw i/o (peeks and gets).
    typename       std::iostream::int_type ch;
    const typename std::iostream::int_type EOFch = 
        std::iostream::traits_type::eof();

    // First we'll look for the optional "~". If found, the brackets become
    // required.
    bool tildeFound = false;
    ch = in.peek(); if (in.fail()) return in;
    assert(ch != EOFch); // we already checked above
    if ((char)ch == '~') {
        tildeFound = true;
        in.get(); // absorb the tilde
        // Eat whitespace after the tilde to see what's next.
        if (in.good()) std::ws(in);
        // If we hit eof after the tilde we don't like the formatting.
        if (!in.good()) {in.setstate(std::ios_base::failbit); return in;}
    }

    // Here the stream is good, the next character is non-white, and we
    // might have seen a tilde.
    assert(in.good());

    // Now see if the sequence is bare or surrounded by (), or [].
    bool lookForCloser = true;
    char openBracket, closeBracket;
    ch = in.peek(); if (in.fail()) return in;
    assert(ch != EOFch); // we already checked above

    openBracket = (char)ch;
    if      (openBracket=='(') {in.get(); closeBracket = ')';}
    else if (openBracket=='[') {in.get(); closeBracket = ']';}
    else lookForCloser = false;

    // If we found a tilde, the opening bracket was mandatory. If we didn't
    // find one then we reject the formatting.
    if (tildeFound && !lookForCloser)
    {   in.setstate(std::ios_base::failbit); return in;}

    // If lookForCloser is true, then closeBracket contains the terminating
    // delimiter, otherwise we're not going to quit until eof.

    // Eat whitespace after the opening bracket to see what's next.
    if (in.good()) std::ws(in);

    // If we're at eof now it must be because the open bracket was the
    // last non-white character in the stream, which is an error.
    if (!in.good()) {
        if (in.eof()) {
            assert(lookForCloser); // or we haven't read anything that could eof
            in.setstate(std::ios::failbit);
        }
        return in;
    }

    // istream is good and next character is non-white; ready to read first
    // value or terminator.

    // We need to figure out whether the elements are space- or comma-
    // separated and then insist on consistency.
    bool commaOK = true, commaRequired = false;
    bool terminatorSeen = false;
    int nextIndex = 0;
    while (true) {
        char c;

        // Here at the top of this loop, we have already successfully read 
        // n=nextIndex values of type T. For fixed-size reads, it might be
        // the case that n==numRequired already, but we still may need to
        // look for a closing bracket before we can declare victory.
        // The stream is good() (not at eof) but it might be the case that 
        // there is nothing but white space left; we don't know yet because
        // if we have satisfied the fixed-size count and are not expecting
        // a terminator then we should quit without absorbing the trailing
        // white space.
        assert(in.good());

        // Look for closing bracket before trying to read value.
        if (lookForCloser) {
            // Eat white space to find the closing bracket.
            std::ws(in); if (!in.good()) break; // eof?
            ch = in.peek(); assert(ch != EOFch);
            if (!in.good()) break;
            c = (char)ch;
            if (c == closeBracket) {   
                in.get(); // absorb the closing bracket
                terminatorSeen = true; 
                break; 
            }
            // next char not a closing bracket; fall through
        }

        // We didn't look or didn't find a closing bracket. The istream is good 
        // but we might be looking at white space.

        // If we already got all the elements we want, break for final checks.
        if (isFixedSize && (nextIndex == numRequired))
            break; // that's a full count.

        // Look for comma before value, except the first time.
        if (commaOK && nextIndex != 0) {
            // Eat white space to find the comma.
            std::ws(in); if (!in.good()) break; // eof?
            ch = in.peek(); assert(ch != EOFch);
            if (!in.good()) break;
            c = (char)ch;
            if (c == ',') {
                in.get(); // absorb comma
                commaRequired = true; // all commas from now on
            } else { // next char not a comma
                if (commaRequired) // bad, e.g.: v1, v2, v3 v4 
                {   in.setstate(std::ios::failbit); break; }
                else commaOK = false; // saw: v1 v2 (no commas now)
            }
            if (!in.good()) break; // might be eof
        }

        // No closing bracket yet; don't have enough elements; skipped comma 
        // if any; istream is good; might be looking at white space.
        assert(in.good());

        // Now read in an element of type T.
        // The extractor T::operator>>() will ignore leading white space.
        if (!isFixedSize)
            out.resizeKeep(out.size()+1); // grow by one (default consructed)
        in >> out[nextIndex]; if (in.fail()) break;
        ++nextIndex;

        if (!in.good()) break; // might be eof
    }

    // We will get here under a number of circumstances:
    //  - the fail bit is set in the istream, or
    //  - we reached eof
    //  - we saw a closing brace
    //  - we got all the elements we wanted (for a fixed-size read)
    // Note that it is possible that we consumed everything except some
    // trailing white space (meaning we're not technically at eof), but
    // for consistency with built-in operator>>()'s we won't try to absorb
    // that trailing white space.

    if (!in.fail()) {
        if (lookForCloser && !terminatorSeen)
            in.setstate(std::ios::failbit); // missing terminator

        if (isFixedSize && nextIndex != numRequired)
            in.setstate(std::ios::failbit); // wrong number of values
    }

    return in;
}



//------------------------------------------------------------------------------
//                          RELATED GLOBAL OPERATORS
//------------------------------------------------------------------------------
// These are logically part of the Matrix_<T> class but are not actually 
// class members; that is, they are in the SimTK namespace.

/**@name             Matrix_<T> serialization and I/O
These methods are at namespace scope but are logically part of the Vector
classes. These deal with reading and writing Vectors from and to streams,
which places an additional requirement on the element type T: the element 
must support the same operation you are trying to do on the Vector as a 
whole. **/
/*@{*/

/** Output a human readable representation of a Vector to an std::ostream
(like std::cout). The format is ~[ \e elements ] where \e elements is a 
space-separated list of the Vector's contents output by invoking the "<<" 
operator on the elements. This function will not compile if the element type 
does not support the "<<" operator. No newline is issued before
or after the output. @relates Vector_ **/
template <class T> inline std::ostream&
operator<<(std::ostream& o, const VectorBase<T>& v)
{ o << "~["; for(int i=0;i<v.size();++i) o<<(i>0?" ":"")<<v[i]; 
    return o << "]"; }

/** Output a human readable representation of a RowVector to an std::ostream
(like std::cout). The format is [ \e elements ] where \e elements is a 
space-separated list of the RowVector's contents output by invoking the "<<" 
operator on the elements. This function will not compile if the element type 
does not support the "<<" operator. No newline is issued before
or after the output. @relates RowVector_ **/
template <class T> inline std::ostream&
operator<<(std::ostream& o, const RowVectorBase<T>& v)
{ o << "["; for(int i=0;i<v.size();++i) o<<(i>0?" ":"")<<v[i]; 
    return o << "]"; }

/** Output a human readable representation of a Matrix to an std::ostream
(like std::cout). The format is one row per line, with each row output as
[ \e elements ] where \e elements is a 
space-separated list of the row's contents output by invoking the "<<" 
operator on the elements. This function will not compile if the element type 
does not support the "<<" operator. A newline is issued before each row and
at the end. @relates Matrix_ **/
template <class T> inline std::ostream&
operator<<(std::ostream& o, const MatrixBase<T>& m) {
    for (int i=0;i<m.nrow();++i)
        o << std::endl << m[i];
    if (m.nrow()) o << std::endl;
    return o; 
}


/** Read in a Vector_<T> from a stream, as a sequence of space-separated or
comma-separated values optionally surrounded by parentheses (), or square 
brackets [], or the "transposed" ~() or ~[]. In the case that the transpose
operator is present, the parentheses or brackets are required, otherwise they
are optional. We will continue to read elements of 
type T from the stream until we find a reason to stop, using type T's stream 
extraction operator>>() to read in each element and resizing the Vector as
necessary. If the data is bracketed, we'll read until we hit the closing 
bracket. If it is not bracketed, we'll read until we hit eof() or get an error
such as the element extractor setting the stream's fail bit due to bad 
formatting. On successful return, the stream will be positioned right after 
the final read-in element or terminating bracket, and the stream's status will 
be good() or eof(). We will not consume trailing whitespace after bracketed 
elements; that means the stream might actually be empty even if we don't 
return eof(). If you want to know whether there is anything else in the 
stream, follow this call with the STL whitespace skipper std::ws() like this:
@code
    if (readVectorFromStream(in,vec) && !in.eof()) 
        std::ws(in); // might take us to eof
    if (in.fail()) {...} // probably a formatting error
    else {
        // Here if the stream is good() then there is more to read; if the
        // stream got used up the status is guaranteed to be eof().
    }
@endcode
A compilation error will occur if you try to use this method on an Vector_<T>
for a type T for which there is no stream extraction operator>>(). 
@note If you want to fill a resizeable Vector_<T> with a fixed amount of data 
from the stream, resize() the Vector to the appropriate length and then use 
fillVectorFromStream() instead. @see fillVectorFromStream()
@relates Vector_ **/
template <class T> static inline 
std::istream& readVectorFromStream(std::istream& in, Vector_<T>& out)
{   return readVectorFromStreamHelper<T>(in, false /*variable sizez*/, out); }



/** Read in a fixed number of elements from a stream into a Vector. We expect 
to read in exactly size() elements of type T, using type T's stream extraction 
operator>>(). This will stop reading when we've read size() elements, or set 
the fail bit in the stream if we run out of elements or if any element's 
extract operator sets the fail bit. On successful return, all size() elements 
will have been set, the stream will be positioned right after the final 
read-in element or terminating bracket, and the stream's status will be good()
or eof(). We will not consume trailing whitespace after reading all the 
elements; that means the stream might actually be empty even if we don't 
return eof(). If you want to know whether there is anything else in the 
stream, follow this call with std::ws() like this:
@code
    if (fillVectorFromStream(in,vec))
        if (!in.eof()) std::ws(in); // might take us to eof
    if (in.fail()) {...} // deal with I/O or formatting error
    // Here if the stream is good() then there is more to read; if the
    // stream got used up the status is guaranteed to be eof().
@endcode
A compilation error will occur if you try to use this method on a Vector_<T>
for a type T for which there is no stream extraction operator>>().
@note If you want to read in a variable number of elements and have the 
Vector_<T> resized as needed, use readVectorFromStream() instead.
@see readVectorFromStream()
@relates Vector_ **/
template <class T> static inline 
std::istream& fillVectorFromStream(std::istream& in, Vector_<T>& out)
{   return readVectorFromStreamHelper<T>(in, true /*fixed size*/, out); }

/** Read in a fixed number of elements from a stream into an VectorView. See
fillVectorFromStream() for more information; this works the same way.
@see fillVectorFromStream()
@relates VectorView_ **/
template <class T> static inline 
std::istream& fillVectorViewFromStream(std::istream& in, VectorView_<T>& out)
{   return readVectorFromStreamHelper<T>(in, true /*fixed size*/, out); }


/** Read Vector_<T> from a stream as a sequence of space- or comma-separated
values of type T, optionally delimited by parentheses, or brackets, and
preceded by "~". The Vector_<T> may be an owner (variable size) or a view 
(fixed size n). In the case of an owner, we'll read all the elements in 
brackets or until eof if there are no brackets. In the case of a view, there 
must be exactly n elements in brackets, or if there are no brackets we'll 
consume exactly n elements and then stop. Each element is read in with its 
own operator ">>" so this won't work if no such operator is defined for 
type T. @relates Vector_ **/
template <class T> inline
std::istream& operator>>(std::istream& in, Vector_<T>& out) 
{   return readVectorFromStream<T>(in, out); }

/** Read a (fixed size n) VectorView_<T> from a stream as a sequence of space- 
or comma-separated values of type T, optionally delimited by parentheses or 
square brackets, and preceded by "~". If there are no delimiters then we will 
read size() values and then stop. Otherwise, there must be exactly size() 
values within the brackets. Each element is read in with its own 
operator ">>" so  this won't work if no such operator is defined for type T.
@relates VectorView_ **/
template <class T> inline
std::istream& operator>>(std::istream& in, VectorView_<T>& out) 
{   return fillVectorViewFromStream<T>(in, out); }

/*@}                     End of Matrix serialization. **/

// Friendly abbreviations for vectors and matrices with scalar elements.

typedef Vector_<Real>           Vector;
typedef Vector_<float>          fVector;
typedef Vector_<double>         dVector;
typedef Vector_<Complex>        ComplexVector;
typedef Vector_<fComplex>       fComplexVector;
typedef Vector_<dComplex>       dComplexVector;

typedef VectorView_<Real>       VectorView;
typedef VectorView_<float>      fVectorView;
typedef VectorView_<double>     dVectorView;
typedef VectorView_<Complex>    ComplexVectorView;
typedef VectorView_<fComplex>   fComplexVectorView;
typedef VectorView_<dComplex>   dComplexVectorView;

typedef RowVector_<Real>        RowVector;
typedef RowVector_<float>       fRowVector;
typedef RowVector_<double>      dRowVector;
typedef RowVector_<Complex>     ComplexRowVector;
typedef RowVector_<fComplex>    fComplexRowVector;
typedef RowVector_<dComplex>    dComplexRowVector;

typedef RowVectorView_<Real>    RowVectorView;
typedef RowVectorView_<float>   fRowVectorView;
typedef RowVectorView_<double>  dRowVectorView;
typedef RowVectorView_<Complex> ComplexRowVectorView;
typedef RowVectorView_<fComplex> fComplexRowVectorView;
typedef RowVectorView_<dComplex> dComplexRowVectorView;

typedef Matrix_<Real>           Matrix;
typedef Matrix_<float>          fMatrix;
typedef Matrix_<double>         dMatrix;
typedef Matrix_<Complex>        ComplexMatrix;
typedef Matrix_<fComplex>       fComplexMatrix;
typedef Matrix_<dComplex>       dComplexMatrix;

typedef MatrixView_<Real>       MatrixView;
typedef MatrixView_<float>      fMatrixView;
typedef MatrixView_<double>     dMatrixView;
typedef MatrixView_<Complex>    ComplexMatrixView;
typedef MatrixView_<fComplex>   fComplexMatrixView;
typedef MatrixView_<dComplex>   dComplexMatrixView;


/**
 * This is an iterator for iterating over the elements of a Vector.
 */

template <class ELT, class VECTOR_CLASS>
class VectorIterator {
public:
    typedef ELT value_type;
    typedef int difference_type;
    typedef ELT& reference;
    typedef ELT* pointer;
    typedef std::random_access_iterator_tag iterator_category;
    VectorIterator(VECTOR_CLASS& vector, int index) : vector(vector), index(index) {
    }
    VectorIterator(const VectorIterator& iter) : vector(iter.vector), index(iter.index) {
    }
    VectorIterator& operator=(const VectorIterator& iter) {
        vector = iter.vector;
        index = iter.index;
        return *this;
    }
    ELT& operator*() {
        assert (index >= 0 && index < vector.size());
        return vector[index];
    }
    ELT& operator[](int i) {
        assert (i >= 0 && i < vector.size());
        return vector[i];
    }
    VectorIterator operator++() {
        assert (index < vector.size());
        ++index;
        return *this;
    }
    VectorIterator operator++(int) {
        assert (index < vector.size());
        VectorIterator current = *this;
        ++index;
        return current;
    }
    VectorIterator operator--() {
        assert (index > 0);
        --index;
        return *this;
    }
    VectorIterator operator--(int) {
        assert (index > 0);
        VectorIterator current = *this;
        --index;
        return current;
    }
    bool operator<(VectorIterator iter) const {
        return (index < iter.index);
    }
    bool operator>(VectorIterator iter) const {
        return (index > iter.index);
    }
    bool operator<=(VectorIterator iter) const {
        return (index <= iter.index);
    }
    bool operator>=(VectorIterator iter) const {
        return (index >= iter.index);
    }
    int operator-(VectorIterator iter) const {
        return (index - iter.index);
    }
    VectorIterator operator-(int n) const {
        return VectorIterator(vector, index-n);
    }
    VectorIterator operator+(int n) const {
        return VectorIterator(vector, index+n);
    }
    bool operator==(VectorIterator iter) const {
        return (index == iter.index);
    }
    bool operator!=(VectorIterator iter) const {
        return (index != iter.index);
    }
private:
    VECTOR_CLASS& vector;
    int index;
};

} //namespace SimTK

#endif //SimTK_SIMMATRIX_BIGMATRIX_H_
