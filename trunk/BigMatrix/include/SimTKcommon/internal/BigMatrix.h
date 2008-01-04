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

/** @file
 * This file defines the client side of the SimTK::Matrix classes, which
 * hold large, variable-sized matrices whose elements are packed SimTK "Composite
 * Numerical Types" (CNTs). Unlike CNTs, the implemention here is opaque,
 * and almost all properties are captured in the implementation at run time
 * rather than in the type at compile time. 
 *
 * Every Matrix consists logically of two pieces: data descriptor and view (which
 * we call an "element filter"). Many matrices can refer to the same data descriptor,
 * but each presents its own view of that data. The default view is transparent,
 * that is, all the data managed by the descriptor is visible. A single Matrix
 * owns a data descriptor, although there can be many other views. 
 *
 * NOTE: Destruction of a data descriptor's owner destructs the
 * data descriptor *regardless* of the presences of other views! I.e., these
 * are not reference counted. TODO: should we change that?
 * 
 * Data descriptors refer to the actual in-memory data containing the 
 * floating point values. Typically the data descriptor owns the data, so
 * they disappear together. However, when working with pre-existing
 * arrays (for example when dealing with Fortran or C), the data descriptor
 * will describe the data layout without owning it so the data persists
 * after the death of the descriptor. In that case it is fine to have
 * multiple data descriptors for the same data.
 *
 *        Matrix handle 
 *            references & may own: MatrixRep (private implementation)
 *        MatrixRep
 *            owns:       view
 *            references & may own:   data descriptor
 *                                    references & may own: data
 *                 
 * A Matrix which is the owner of its data descriptor will be resized whenever
 * necessary, unless you take active steps to prevent that. For example, if
 * you declare a Vector, the number of rows can resize but the number of
 * columns will be locked at 1. A RowVector does the reverse. You can also
 * explicitly lock the number of rows and/or columns of a matrix to prevent
 * unwanted resizes.
 *
 * Here are the classes and short descriptions:
 *   MatrixHelper<S>  interface to the opaque implementation, templatized
 *                      by scalar type only
 *   MatrixBase<CNT>  fully templatized client, contains a MatrixHelper
 *
 *   The rest are dataless classes all of which can be interconverted just
 *   by recasting. Every one of these classes has a default conversion to
 *   type Matrix_<same element type>, so users can write functions that expect
 *   a Matrix argument and pass it a Vector, RowVectorView, or whatever.
 *
 *   VectorBase<CNT>    these are derived from MatrixBase and add no new data,
 *   RowVectorBase<CNT> but change some of the operators and other methods to
 *                        be appropriate for 1d data.
 *
 *   Matrix<CNT>      2d owner class     (a MatrixBase<CNT>)
 *   Vector<CNT>      column owner class (a VectorBase<CNT>)
 *   RowVector<CNT>   row owner class    (a RowVectorBase<CNT>)
 *
 *   Views are exactly the same as the corresponding owner class, but with
 *   shallow construction and assignment semantics.
 *
 *   MatrixView<CNT>, VectorView<CNT>, RowVectorView<CNT>
 *
 *   Dead matrices are owners which are about to be destructed. Anything
 *   they own may be taken from them, including the data descriptor and/or
 *   the element filter. This is a very effective performance trick for sequences
 *   of operations since it eliminates most of the need for allocating and
 *   deallocating temporaries.
 *
 *   DeadMatrix<CNT>, DeadVector<CNT>, DeadRowVector<CNT>
 *
 * TODO: matrix expression templates for delaying operator execution.
 */

#include "SimTKcommon/Scalar.h"
#include "SimTKcommon/SmallMatrix.h"

#include "SimTKcommon/internal/MatrixHelper.h"

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
template <class T=Real> class TmpMatrixView_;
template <class T=Real> class Matrix_;
template <class T=Real> class DeadMatrix_;

template <class T=Real> class VectorView_;
template <class T=Real> class TmpVectorView_;
template <class T=Real> class Vector_;
template <class T=Real> class DeadVector_;

template <class T=Real> class RowVectorView_;
template <class T=Real> class TmpRowVectorView_;
template <class T=Real> class RowVector_;
template <class T=Real> class DeadRowVector_;
 
/**
 * Variable-size 2d matrix of Composite Numerical Type (ELT) elements. This is
 * a container of such elements, it is NOT a Composite Numerical Type itself.
 * MatrixBase<ELT> uses MatrixHelper<S> for implementation, where S is 
 * ELT::Scalar, that is, the underlying float, double, long double,
 * complex<float>, negator<conjugate<long double>>, 
 * etc. from which ELT is constructed. This is a finite set of which all
 * members are explicitly instantiated in the implementation code, so 
 * clients don't have to know how anything is implemented.
 * 
 * MatrixBase is the only class in the Matrix/Vector family which has any
 * data members (it has exactly one MatrixHelper). Thus all other objects
 * in this family (that is, derived from MatrixBase) are exactly the same
 * size in memory and may be "reinterpreted" as appropriate. For example,
 * a Vector may be reinterpreted as a Matrix or vice versa, provided runtime
 * requirements are met (e.g., exactly 1 column).
 *
 * Unlike the small matrix classes, very little is encoded in the type.
 * Only the element type, and matrix vs. vector vs. row are in the type;
 * everything else like shape, storage layout, and writability are handled
 * at run time. 
 */
template <class ELT> class MatrixBase {  
public:
    // These typedefs are handy, but despite appearances you cannot 
    // treat a MatrixBase as a composite numerical type. That is,
    // CNT<MatrixBase> will not compile, or if it does it won't be
    // meaningful.

    typedef ELT                        E;
    typedef typename CNT<E>::TNeg      ENeg;
    typedef typename CNT<E>::TAbs      EAbs;
    typedef typename CNT<E>::TReal     EReal;
    typedef typename CNT<E>::TImag     EImag;
    typedef typename CNT<E>::TComplex  EComplex;
    typedef typename CNT<E>::THerm     EHerm;       
    typedef typename CNT<E>::TPosTrans EPosTrans;
    typedef typename CNT<E>::TSqHermT  ESqHermT;
    typedef typename CNT<E>::TSqTHerm  ESqTHerm;

    typedef typename CNT<E>::Scalar    Scalar;        // the underlying Scalar type
    typedef typename CNT<E>::Number    Number;        // negator removed from Scalar
    typedef typename CNT<E>::StdNumber StdNumber;     // conjugate goes to complex
    typedef typename CNT<E>::Precision Precision;     // complex removed from StdNumber

    typedef typename CNT<ELT>::ScalarSq  ScalarSq;

    typedef MatrixBase<E>                T;
    typedef MatrixBase<ENeg>             TNeg;
    typedef MatrixBase<EAbs>             TAbs;
    typedef MatrixBase<EReal>            TReal;
    typedef MatrixBase<EImag>            TImag;
    typedef MatrixBase<EComplex>         TComplex;
    typedef MatrixBase<EHerm>            THerm; 
    typedef MatrixBase<E>                TPosTrans;
    typedef MatrixBase<ESqHermT>         TSqHermT;  // ~Mat*Mat
    typedef MatrixBase<ESqTHerm>         TSqTHerm;  // Mat*~Mat

    void setMatrixStructure(MatrixStructures::Structure structure) {
        helper.setMatrixStructure( structure);  // default Uncommitted
    }
    MatrixStructures::Structure getMatrixStructure() const {
        return helper.getMatrixStructure();
    }
    void setMatrixShape(MatrixShapes::Shape shape) {
        helper.setMatrixShape( shape);          // default Uncommitted
    }
    MatrixShapes::Shape getMatrixShape() const  {
        return helper.getMatrixShape();
    }
    void setMatrixSparsity(MatrixSparseFormats::Sparsity sparsity) {
        helper.setMatrixSparsity( sparsity);    // default Uncommitted
    }
    MatrixSparseFormats::Sparsity getMatrixSparsity() const  {
        return helper.getMatrixSparsity();
    }
    void setMatrixStorage(MatrixStorageFormats::Storage storage) {
        helper.setMatrixStorage( storage);      // default Uncommitted
    }
    MatrixStorageFormats::Storage getMatrixStorage() const  {
        return helper.getMatrixStorage();
    }

    void setMatrixCondition(MatrixConditions::Condition condition) {
        helper.setMatrixCondition(condition);  //  default Uncommitted
    }
    MatrixConditions::Condition getMatrixCondition() const  {
        return helper.getMatrixCondition();
    }


    // This gives the resulting matrix type when (m(i,j) op P) is applied to each element.
    // It will have element types which are the regular composite result of E op P.
    template <class P> struct EltResult { 
        typedef MatrixBase<typename CNT<E>::template Result<P>::Mul> Mul;
        typedef MatrixBase<typename CNT<E>::template Result<P>::Dvd> Dvd;
        typedef MatrixBase<typename CNT<E>::template Result<P>::Add> Add;
        typedef MatrixBase<typename CNT<E>::template Result<P>::Sub> Sub;
    };
       
    // Product of dimensions may be > 32 bits.
    long size() const { return helper.size(); }
    int  nrow() const { return helper.nrow(); }
    int  ncol() const { return helper.ncol(); }

    // Scalar norm square is sum( squares of all scalars )
    // TODO: very slow! Should be optimized at least for the case
    //       where ELT is a Scalar.
    ScalarSq scalarNormSqr() const {
        const int nr=nrow(), nc=ncol();
        ScalarSq sum(0);
        for(int j=0;j<nc;++j) 
            for (int i=0; i<nr; ++i)
                sum += CNT<E>::scalarNormSqr((*this)(i,j));
        return sum;
    }

    // abs() is elementwise absolute value; that is, the return value has the same
    // dimension as this Matrix but with each element replaced by whatever it thinks
    // its absolute value is.
    // TODO: very slow! Should be optimized at least for the case
    //       where ELT is a Scalar.
    void abs(TAbs& mabs) const {
        const int nr=nrow(), nc=ncol();
        mabs.resize(nr,nc);
        for(int j=0;j<nc;++j) 
            for (int i=0; i<nr; ++i)
                mabs(i,j) = CNT<E>::abs((*this)(i,j));
    }

    TAbs abs() const { TAbs mabs; abs(mabs); return mabs; }

    enum { 
        NScalarsPerElement    = CNT<E>::NActualScalars,
        CppNScalarsPerElement = sizeof(E) / sizeof(Scalar)
    };
  
    MatrixBase() : helper(NScalarsPerElement,CppNScalarsPerElement)  { }

    // This restores the MatrixBase to its just-constructed state.
    void clear() {
        helper.clear();
    }

    
    // Copy constructor is a deep copy (not appropriate for views!).    
    MatrixBase(const MatrixBase& b)
      : helper(b.helper, typename MatrixHelper<Scalar>::DeepCopy()) { }
    
    // Copy assignment is a deep copy but behavior depends on type of lhs: if view, rhs
    // must match. If owner, we reallocate and copy rhs.
    MatrixBase& copyAssign(const MatrixBase& b) {
        helper.copyAssign(b.helper);
        return *this;
    }
    MatrixBase& operator=(const MatrixBase& b) { return copyAssign(b); }

    // View assignment is a shallow copy, meaning that we disconnect the MatrixBase 
    // from whatever it used to refer to (destructing as necessary), then make it a new view
    // for the data descriptor referenced by the source.
    // CAUTION: we always take the source as const, but that is ignored in 
    // determining whether the resulting view is writable. Instead, that is
    // inherited from the writability status of the source. We have to do this
    // in order to allow temporary view objects to be writable -- the compiler
    // creates temporaries like m(i,j,m,n) as const.

    MatrixBase& viewAssign(const MatrixBase& src) {
        helper.writableViewAssign(const_cast<MatrixHelper<Scalar>&>(src.helper));
        return *this;
    }

    // default destructor


    MatrixBase(int m, int n, bool lockNrow=false, bool lockNcol=false) 
      : helper(NScalarsPerElement, CppNScalarsPerElement, m, n, lockNrow, lockNcol)  { }
        
    // Non-resizeable owner sharing pre-allocated raw data
    MatrixBase(int m, int n, int leadingDim, const Scalar* s)
      : helper(NScalarsPerElement, CppNScalarsPerElement, m, n, leadingDim, s) { }    // read only
    MatrixBase(int m, int n, int leadingDim, Scalar* s)
      : helper(NScalarsPerElement, CppNScalarsPerElement, m, n, leadingDim, s) { }    // writable
            
    MatrixBase(int m, int n, const ELT& t) 
      : helper(NScalarsPerElement, CppNScalarsPerElement, m, n)
      { helper.fillWith(reinterpret_cast<const Scalar*>(&t)); }    
    MatrixBase(int m, int n, bool lockNrow, bool lockNcol, const ELT& t) 
      : helper(NScalarsPerElement, CppNScalarsPerElement, m, n, lockNrow, lockNcol)
      { helper.fillWith(reinterpret_cast<const Scalar*>(&t)); }
        
    MatrixBase(int m, int n, const ELT* p) 
      : helper(NScalarsPerElement, CppNScalarsPerElement, m, n)
      { helper.copyInByRows(reinterpret_cast<const Scalar*>(p)); }
    MatrixBase(int m, int n,  bool lockNrow, bool lockNcol, const ELT* p) 
      : helper(NScalarsPerElement, CppNScalarsPerElement, m, n, lockNrow, lockNcol)
      { helper.copyInByRows(reinterpret_cast<const Scalar*>(p)); }
        
    // Create a new MatrixBase from an existing helper. Both shallow and deep copies are possible.
    MatrixBase(MatrixHelper<Scalar>&       h, const typename MatrixHelper<Scalar>::ShallowCopy& s) : helper(h,s) { }
    MatrixBase(const MatrixHelper<Scalar>& h, const typename MatrixHelper<Scalar>::ShallowCopy& s) : helper(h,s) { }
    MatrixBase(const MatrixHelper<Scalar>& h, const typename MatrixHelper<Scalar>::DeepCopy& d)    : helper(h,d) { }


    MatrixBase& operator*=(const StdNumber& t)  { helper.scaleBy(t);              return *this; }
    MatrixBase& operator/=(const StdNumber& t)  { helper.scaleBy(StdNumber(1)/t); return *this; }
    MatrixBase& operator+=(const MatrixBase& r) { helper.addIn(r.helper);         return *this; }
    MatrixBase& operator-=(const MatrixBase& r) { helper.subIn(r.helper);         return *this; }  

    template <class EE> MatrixBase(const MatrixBase<EE>& b)
      : helper(b.helper, typename MatrixHelper<Scalar>::DeepCopy()) { }

    template <class EE> MatrixBase& operator=(const MatrixBase<EE>& b) 
      { helper = b.helper; return *this; }
    template <class EE> MatrixBase& operator+=(const MatrixBase<EE>& b) 
      { helper.addIn(b.helper); return *this; }
    template <class EE> MatrixBase& operator-=(const MatrixBase<EE>& b) 
      { helper.subIn(b.helper); return *this; }

    // Matrix assignment to an element sets only the *diagonal* elements to
    // the indicated value; everything else is set to zero. This is particularly
    // useful for setting a Matrix to zero or to the identity; for other values
    // it create as Matrix which acts like the scalar. That is, if the scalar
    // is s and we do M=s, then multiplying another Matrix B by the resulting 
    // diagonal matrix M gives the same result as multiplying B by s. That is
    // (M=s)*B == s*B.
    //
    // NOTE: this must be overridden for Vector and RowVector since then scalar
    // assignment is defined to copy the scalar to every element.
    MatrixBase& operator=(const ELT& t) { 
        setToZero(); updDiag().setTo(t); 
        return *this;
    }

	// M = diag(v) * M; v must have nrow() elements.
	// That is, M[i] *= v[i].
	template <class EE> inline MatrixBase& rowScaleInPlace(const VectorBase<EE>&);

	// Return type is a new matrix which will have the same dimensions as 'this' but
	// will have element types appropriate for the elementwise multiply being performed.
	template <class EE> inline void rowScale(const VectorBase<EE>&    v, typename EltResult<EE>::Mul& out) const;
	template <class EE> inline typename EltResult<EE>::Mul rowScale(const VectorBase<EE>& v) const
      { EltResult<EE>::Mul out(nrow(), ncol()); rowScale(v,out); return out; }

	// M = diag(v)^-1 * M; v must have nrow() elements.
	// That is, M[i] /= v[i].
	template <class EE> inline MatrixBase& rowUnscaleInPlace(const VectorBase<EE>&);
	template <class EE> inline void rowUnscale(const VectorBase<EE>& v, typename EltResult<EE>::Dvd& out) const;
	template <class EE> inline typename EltResult<EE>::Dvd rowUnscale(const VectorBase<EE>& v) const
      { EltResult<EE>::Dvd out(nrow(), ncol()); rowUnscale(v,out); return out; }

	// M = M * diag(v); v must have ncol() elements
	// That is, M(i) *= v[i]
	template <class EE> inline MatrixBase& colScaleInPlace(const VectorBase<EE>&);
	template <class EE> inline void colScale(const VectorBase<EE>& v, typename EltResult<EE>::Mul& out) const;
	template <class EE> inline typename EltResult<EE>::Mul colScale(const VectorBase<EE>& v) const
      { EltResult<EE>::Mul out(nrow(), ncol()); colScale(v,out); return out; }

	// M = M * diag(v)^-1; v must have ncol() elements
	// That is, M(i) /= v[i]
	template <class EE> inline MatrixBase& colUnscaleInPlace(const VectorBase<EE>&);
	template <class EE> inline void colUnscale(const VectorBase<EE>& v, typename EltResult<EE>::Dvd& out) const;
	template <class EE> inline typename EltResult<EE>::Dvd colUnscale(const VectorBase<EE>& v) const
      { EltResult<EE>::Dvd out(nrow(), ncol()); colUnscale(v,out); return out; }

	// M(i,j) *= R(i,j); R must have same dimensions as this
	template <class EE> inline MatrixBase& elementwiseMultiplyInPlace(const MatrixBase<EE>&);
	template <class EE> inline void elementwiseMultiply(const MatrixBase<EE>&, typename EltResult<EE>::Mul&) const;
	template <class EE> inline typename EltResult<EE>::Mul elementwiseMultiply(const MatrixBase<EE>&) const
      { EltResult<EE>::Mul out(nrow(), ncol()); elementwiseMultiply(v,out); return out; }

	// M(i,j) /= R(i,j); R must have same dimensions as this
	template <class EE> inline MatrixBase& elementwiseDivideInPlace(const MatrixBase<EE>& r);
	template <class EE> inline void elementwiseDivide(const MatrixBase<EE>&, typename EltResult<EE>::Dvd&) const;
	template <class EE> inline typename EltResult<EE>::Mul elementwiseDivide(const MatrixBase<EE>&) const
      { EltResult<EE>::Dvd out(nrow(), ncol()); elementwiseDivide(v,out); return out; }

    // fill every element in current allocation with given element (or NaN or 0)
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

    // Select matrix diagonal (of largest leading square if rectangular).
    inline VectorView_<ELT> diag() const;
    inline VectorView_<ELT> updDiag();

    // Create a view of the real or imaginary elements. TODO
    //inline MatrixView_<EReal> real() const;
    //inline MatrixView_<EReal> updReal();
    //inline MatrixView_<EImag> imag() const;
    //inline MatrixView_<EImag> updImag();

    // Overload "real" and "imag" for both read and write as a nod to convention. TODO
    //MatrixView_<EReal> real() {return updReal();}
    //MatrixView_<EReal> imag() {return updImag();}

    // TODO: this routine seems ill-advised but I need it for the IVM port at the moment
    Matrix_<StdNumber> invert() const {  // return a newly-allocated inverse; dump negator 
        Matrix_<StdNumber> m(*this);
        m.helper.invertInPlace();
        return m;   // TODO - bad: makes an extra copy
    }

    // Matlab-compatible debug output.
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

    // Element selection      
    const ELT& operator()(int i, int j) const { return *reinterpret_cast<const ELT*>(helper.getElt(i,j)); }
    ELT&       operator()(int i, int j)       { return *reinterpret_cast<      ELT*>(helper.updElt(i,j)); }

    // This is the scalar Frobenius norm, and its square. Note: if this is a Matrix then the Frobenius
    // norm is NOT the same as the 2-norm, although they are equivalent for Vectors.
    ScalarSq normSqr() const { return scalarNormSqr(); }
    ScalarSq norm()    const { return std::sqrt(scalarNormSqr()); } // TODO -- not good; unnecessary overflow

    // We only allow RMS norm if the elements are scalars. If there are no elements in this Matrix,
    // we'll define its RMS norm to be 0, although NaN might be a better choice.
    ScalarSq normRMS() const {
        if (!CNT<ELT>::IsScalar)
            SimTK_THROW1(Exception::Cant, "normRMS() only defined for scalar elements");
        const long nelt = (long)nrow()*(long)ncol();
        if (nelt == 0)
            return ScalarSq(0);
        return std::sqrt(scalarNormSqr()/nelt);
    }


    ELT sum() const { ELT s; helper.sum(reinterpret_cast<Scalar*>(&s)); return s; } // add all the elements        

    //TODO: make unary +/- return a self-view so they won't reallocate?
    const MatrixBase& operator+() const {return *this; }
    const TNeg&       negate()    const {return *reinterpret_cast<const TNeg*>(this); }
    TNeg&             updNegate()       {return *reinterpret_cast<TNeg*>(this); }

    const TNeg&       operator-() const {return negate();}
    TNeg&             operator-()       {return updNegate();}
 
    MatrixBase& resize(int m, int n)     { helper.resize(m,n); return *this; }
    MatrixBase& resizeKeep(int m, int n) { helper.resizeKeep(m,n); return *this; }

    // These prevent shape changes in a Matrix that would otherwise allow it. No harm if they
    // are called on a Matrix that is locked already; they always succeed.
    void lockNRows() {helper.lockNRows();}
    void lockNCols() {helper.lockNCols();}
    void lockShape() {helper.lockShape();}

    // These allow shape changes again for a Matrix which was constructed to allow them
    // but had them locked with the above routines. No harm if these are called on a Matrix
    // that is already unlocked, but it is not allowed to call them on a Matrix which
    // *never* allowed resizing. An exception will be thrown in that case.
    void unlockNRows() {helper.unlockNRows();}
    void unlockNCols() {helper.unlockNCols();}
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

    // This is the number of consecutive scalars used to represent one
    // element of type ELT. This may be fewer than C++ would use for the
    // element, since it may introduce some padding.
    int getNScalarsPerElement()  const {return NScalarsPerElement;}

    // This is like sizeof(ELT), but returning the number of bytes we use
    // to store the element which may be fewer than what C++ would use.
    int getPackedSizeofElement() const {return NScalarsPerElement*sizeof(Scalar);}

    bool hasContiguousData() const {return helper.hasContiguousData();}
    long getContiguousScalarDataLength() const {
        return helper.getContiguousDataLength();
    }
    const Scalar* getContiguousScalarData() const {
        return helper.getContiguousData();
    }
    Scalar* updContiguousScalarData() {
        return helper.updContiguousData();
    }
    void replaceContiguousScalarData(Scalar* newData, long length, bool takeOwnership) {
        helper.replaceContiguousData(newData,length,takeOwnership);
    }
    void replaceContiguousScalarData(const Scalar* newData, long length) {
        helper.replaceContiguousData(newData,length);
    }
    void swapOwnedContiguousScalarData(Scalar* newData, int length, Scalar*& oldData) {
        helper.swapOwnedContiguousData(newData,length,oldData);
    }
protected:
    MatrixHelper<Scalar> helper; // this is just one pointer

    template <class EE> friend class MatrixBase;
};

/**
 * This is a dataless rehash of the MatrixBase class to specialize it for Vectors.
 * This mostly entails overriding a few of the methods. Note that all the MatrixBase
 * operations remain available if you static_cast<> this up to a MatrixBase.
 */
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
    VectorBase() : Base() { }

    // Copy constructor is a deep copy (not appropriate for views!).
    VectorBase(const VectorBase& b) : Base(b) { }
    
    // Copy assignment is deep copy but behavior depends on type of lhs: if view, rhs
    // must match. If owner, we reallocate and copy rhs.
    VectorBase& operator=(const VectorBase& b) {
        Base::operator=(b); return *this;
    }

    // default destructor


    explicit VectorBase(int m, bool lockNrow=false)
      : Base(m,1,lockNrow,true) { }
        
    // Non-resizeable owner sharing pre-allocated raw data
    VectorBase(int m, int leadingDim, const Scalar* s)
      : Base(m,1,leadingDim,s) { }
    VectorBase(int m, int leadingDim, Scalar* s)
      : Base(m,1,leadingDim,s) { }
            
    VectorBase(int m, const ELT& t) : Base(m,1,t) { }  
    VectorBase(int m, bool lockNrow, const ELT& t) 
      : Base(m,1,lockNrow,true,t) { }
        
    VectorBase(int m, const ELT* p) : Base(m,1,p) { }
    VectorBase(int m, bool lockNrow, const ELT* p) 
      : Base(m, 1, lockNrow, true, p) { }
        
    // Create a new VectorBase from an existing helper. Both shallow and deep copies are possible.
    VectorBase(MatrixHelper<Scalar>&       h, const typename MatrixHelper<Scalar>::ShallowCopy& s) : Base(h,s) { }
    VectorBase(const MatrixHelper<Scalar>& h, const typename MatrixHelper<Scalar>::ShallowCopy& s) : Base(h,s) { }
    VectorBase(const MatrixHelper<Scalar>& h, const typename MatrixHelper<Scalar>::DeepCopy& d)    : Base(h,d) { }

    // This gives the resulting vector type when (v[i] op P) is applied to each element.
    // It will have element types which are the regular composite result of ELT op P.
    template <class P> struct EltResult { 
        typedef VectorBase<typename CNT<ELT>::template Result<P>::Mul> Mul;
        typedef VectorBase<typename CNT<ELT>::template Result<P>::Dvd> Dvd;
        typedef VectorBase<typename CNT<ELT>::template Result<P>::Add> Add;
        typedef VectorBase<typename CNT<ELT>::template Result<P>::Sub> Sub;
    };

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


    // Fill current allocation with copies of element. Note that this is not the 
    // same behavior as assignment for Matrices, where only the diagonal is set (and
    // everything else is set to zero.)
    VectorBase& operator=(const ELT& t) { Base::setTo(t); return *this; }  

	// There's only one column here so it's a bit wierd to use rowScale rather than
	// elementwiseMultiply, but there's nothing really wrong with it. Using colScale
	// would be really wacky since it is the same as a scalar multiply. We won't support
	// colScale here except through inheritance where it will not be much use.
	template <class EE> VectorBase& rowScaleInPlace(const VectorBase<EE>& v)
	  { Base::rowScaleInPlace<EE>(v); return *this; }
	template <class EE> inline void rowScale(const VectorBase<EE>& v, typename EltResult<EE>::Mul& out) const
	  { return Base::rowScale(v,out); }
	template <class EE> inline typename EltResult<EE>::Mul rowScale(const VectorBase<EE>& v) const
	  { typename EltResult<EE>::Mul out(nrow()); Base::rowScale(v,out); return out; }

	template <class EE> VectorBase& rowUnscaleInPlace(const VectorBase<EE>& v)
	  { Base::rowUnscaleInPlace<EE>(v); return *this; }
	template <class EE> inline void rowUnscale(const VectorBase<EE>& v, typename EltResult<EE>::Dvd& out) const
	  { return Base::rowUnscale(v,out); }
	template <class EE> inline typename EltResult<EE>::Dvd rowUnscale(const VectorBase<EE>& v) const
	  { typename EltResult<EE>::Dvd out(nrow()); Base::rowUnscale(v,out); return out; }

	template <class EE> VectorBase& elementwiseMultiplyInPlace(const VectorBase<EE>& r)
	  { Base::elementwiseMultiplyInPlace<EE>(r); return *this; }
	template <class EE> inline void elementwiseMultiply(const VectorBase<EE>& v, typename EltResult<EE>::Mul& out) const
	  { return Base::elementwiseMultiply(v,out); }
	template <class EE> inline typename EltResult<EE>::Mul elementwiseMultiply(const VectorBase<EE>& v) const
	  { typename EltResult<EE>::Mul out(nrow()); Base::elementwiseMultiply(v,out); return out; }

	template <class EE> VectorBase& elementwiseDivideInPlace(const VectorBase<EE>& r)
	  { Base::elementwiseDivideInPlace<EE>(r); return *this; }
	template <class EE> inline void elementwiseDivide(const VectorBase<EE>& v, typename EltResult<EE>::Dvd& out) const
	  { return Base::elementwiseDivide(v,out); }
	template <class EE> inline typename EltResult<EE>::Mul elementwiseDivide(const VectorBase<EE>& v) const
	  { typename EltResult<EE>::Dvd out(nrow()); Base::elementwiseDivide(v,out); return out; }

    // Implicit conversions are allowed to Vector or Matrix, but not to RowVector.   
    operator const Vector_<ELT>&()     const { return *reinterpret_cast<const Vector_<ELT>*>(this); }
    operator       Vector_<ELT>&()           { return *reinterpret_cast<      Vector_<ELT>*>(this); }
    operator const VectorView_<ELT>&() const { return *reinterpret_cast<const VectorView_<ELT>*>(this); }
    operator       VectorView_<ELT>&()       { return *reinterpret_cast<      VectorView_<ELT>*>(this); }
    
    operator const Matrix_<ELT>&()     const { return *reinterpret_cast<const Matrix_<ELT>*>(this); }
    operator       Matrix_<ELT>&()           { return *reinterpret_cast<      Matrix_<ELT>*>(this); } 
    operator const MatrixView_<ELT>&() const { return *reinterpret_cast<const MatrixView_<ELT>*>(this); }
    operator       MatrixView_<ELT>&()       { return *reinterpret_cast<      MatrixView_<ELT>*>(this); } 


    // Override MatrixBase size() for Vectors to return int instead of long.
    int size() const {
        assert(Base::size() <= std::numeric_limits<int>::max()); 
        return (int)Base::size();
    }

    // Override MatrixBase operators to return the right shape
    TAbs abs() const {TAbs result; Base::abs(result); return result;}
    
    // Override MatrixBase indexing operators          
    const ELT& operator[](int i) const {return Base::operator()(i,0);}
    ELT&       operator[](int i)       {return Base::operator()(i,0);}
    const ELT& operator()(int i) const {return Base::operator()(i,0);}
    ELT&       operator()(int i)       {return Base::operator()(i,0);}
         
    // View creation      
    VectorView_<ELT> operator()(int i, int m) const {return Base::operator()(i,0,m,1).getAsVectorView();}
    VectorView_<ELT> operator()(int i, int m)       {return Base::operator()(i,0,m,1).updAsVectorView();}
 
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
private:
    // NO DATA MEMBERS ALLOWED
};


/**
 * This is a dataless rehash of the MatrixBase class to specialize it for RowVectors.
 * This mostly entails overriding a few of the methods. Note that all the MatrixBase
 * operations remain available if you static_cast<> this up to a MatrixBase.
 */
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
    RowVectorBase() : Base() { }
    
    // Copy constructor is a deep copy (not appropriate for views!).    
    RowVectorBase(const RowVectorBase& b) : Base(b) { }

    // Copy assignment is deep copy but behavior depends on type of lhs: if view, rhs
    // must match. If owner, we reallocate and copy rhs.
    RowVectorBase& operator=(const RowVectorBase& b) {
        Base::operator=(b); return *this;
    }

    // default destructor

    explicit RowVectorBase(int n, bool lockNcol=false)
      : Base(1,n,true,lockNcol) { }
        
    // Non-resizeable owner sharing pre-allocated raw data
    RowVectorBase(int n, int leadingDim, const Scalar* s)
      : Base(1,n,leadingDim,s) { }
    RowVectorBase(int n, int leadingDim, Scalar* s)
      : Base(1,n,leadingDim,s) { }
            
    RowVectorBase(int n, const ELT& t) : Base(1,n,t) { }  
    RowVectorBase(int n, bool lockNcol, const ELT& t) 
      : Base(1,n,true,lockNcol,t) { }
        
    RowVectorBase(int n, const ELT* p) : Base(1,n,p) { }
    RowVectorBase(int n, bool lockNcol, const ELT* p) 
      : Base(1, n, true, lockNcol, p) { }
        
    // Create a new RowVectorBase from an existing helper. Both shallow and deep copies are possible.
    RowVectorBase(MatrixHelper<Scalar>&       h, const typename MatrixHelper<Scalar>::ShallowCopy& s) : Base(h,s) { }
    RowVectorBase(const MatrixHelper<Scalar>& h, const typename MatrixHelper<Scalar>::ShallowCopy& s) : Base(h,s) { }
    RowVectorBase(const MatrixHelper<Scalar>& h, const typename MatrixHelper<Scalar>::DeepCopy& d)    : Base(h,d) { }

	// This gives the resulting rowvector type when (r(i) op P) is applied to each element.
    // It will have element types which are the regular composite result of ELT op P.
    template <class P> struct EltResult { 
        typedef RowVectorBase<typename CNT<ELT>::template Result<P>::Mul> Mul;
        typedef RowVectorBase<typename CNT<ELT>::template Result<P>::Dvd> Dvd;
        typedef RowVectorBase<typename CNT<ELT>::template Result<P>::Add> Add;
        typedef RowVectorBase<typename CNT<ELT>::template Result<P>::Sub> Sub;
    };

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
 
    // Fill current allocation with copies of element. Note that this is not the 
    // same behavior as assignment for Matrices, where only the diagonal is set (and
    // everything else is set to zero.)
    RowVectorBase& operator=(const ELT& t) { Base::setTo(t); return *this; } 

	// There's only one row here so it's a bit wierd to use colScale rather than
	// elementwiseMultiply, but there's nothing really wrong with it. Using rowScale
	// would be really wacky since it is the same as a scalar multiply. We won't support
	// rowScale here except through inheritance where it will not be much use.

	template <class EE> RowVectorBase& colScaleInPlace(const VectorBase<EE>& v)
	  { Base::colScaleInPlace<EE>(v); return *this; }
	template <class EE> inline void colScale(const VectorBase<EE>& v, typename EltResult<EE>::Mul& out) const
	  { return Base::colScale(v,out); }
	template <class EE> inline typename EltResult<EE>::Mul colScale(const VectorBase<EE>& v) const
	  { typename EltResult<EE>::Mul out(ncol()); Base::colScale(v,out); return out; }

	template <class EE> RowVectorBase& colUnscaleInPlace(const VectorBase<EE>& v)
	  { Base::colUnscaleInPlace<EE>(v); return *this; }
	template <class EE> inline void colUnscale(const VectorBase<EE>& v, typename EltResult<EE>::Dvd& out) const
	  { return Base::colUnscale(v,out); }
	template <class EE> inline typename EltResult<EE>::Dvd colUnscale(const VectorBase<EE>& v) const
	  { typename EltResult<EE>::Dvd out(ncol()); Base::colUnscale(v,out); return out; }

	template <class EE> RowVectorBase& elementwiseMultiplyInPlace(const VectorBase<EE>& r)
	  { Base::elementwiseMultiplyInPlace<EE>(r); return *this; }
	template <class EE> inline void elementwiseMultiply(const VectorBase<EE>& v, typename EltResult<EE>::Mul& out) const
	  { return Base::elementwiseMultiply(v,out); }
	template <class EE> inline typename EltResult<EE>::Mul elementwiseMultiply(const VectorBase<EE>& v) const
	  { typename EltResult<EE>::Mul out(ncol()); Base::elementwiseMultiply(v,out); return out; }

	template <class EE> RowVectorBase& elementwiseDivideInPlace(const VectorBase<EE>& r)
	  { Base::elementwiseDivideInPlace<EE>(r); return *this; }
	template <class EE> inline void elementwiseDivide(const VectorBase<EE>& v, typename EltResult<EE>::Dvd& out) const
	  { return Base::elementwiseDivide(v,out); }
	template <class EE> inline typename EltResult<EE>::Mul elementwiseDivide(const VectorBase<EE>& v) const
	  { typename EltResult<EE>::Dvd out(ncol()); Base::elementwiseDivide(v,out); return out; }

    // Implicit conversions are allowed to RowVector or Matrix, but not to Vector.   
    operator const RowVector_<ELT>&()     const {return *reinterpret_cast<const RowVector_<ELT>*>(this);}
    operator       RowVector_<ELT>&()           {return *reinterpret_cast<      RowVector_<ELT>*>(this);}
    operator const RowVectorView_<ELT>&() const {return *reinterpret_cast<const RowVectorView_<ELT>*>(this);}
    operator       RowVectorView_<ELT>&()       {return *reinterpret_cast<      RowVectorView_<ELT>*>(this);}
    
    operator const Matrix_<ELT>&()     const {return *reinterpret_cast<const Matrix_<ELT>*>(this);}
    operator       Matrix_<ELT>&()           {return *reinterpret_cast<      Matrix_<ELT>*>(this);} 
    operator const MatrixView_<ELT>&() const {return *reinterpret_cast<const MatrixView_<ELT>*>(this);}
    operator       MatrixView_<ELT>&()       {return *reinterpret_cast<      MatrixView_<ELT>*>(this);} 
    

    // Override MatrixBase size() for RowVectors to return int instead of long.
    int size() const { 
        assert(Base::size() <= std::numeric_limits<int>::max()); 
        return (int)Base::size();
    }

    // Override MatrixBase operators to return the right shape
    TAbs abs() const {
        TAbs result; Base::abs(result); return result;
    }

    // Override MatrixBase indexing operators          
    const ELT& operator[](int j) const {return Base::operator()(0,j);}
    ELT&       operator[](int j)       {return Base::operator()(0,j);}
    const ELT& operator()(int j) const {return Base::operator()(0,j);}
    ELT&       operator()(int j)       {return Base::operator()(0,j);}
         
    // View creation      
    RowVectorView_<ELT> operator()(int j, int n) const {return Base::operator()(0,j,1,n).getAsRowVectorView();}
    RowVectorView_<ELT> operator()(int j, int n)       {return Base::operator()(0,j,1,n).updAsRowVectorView();}
 
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
private:
    // NO DATA MEMBERS ALLOWED
};

/**
 * XXX not ready for prime time
 */	
template <class ELT> class TmpMatrixView_ : public MatrixBase<ELT> {
    typedef MatrixBase<ELT> Base;
    typedef typename MatrixBase<ELT>::Scalar Scalar;
public:
    TmpMatrixView_() : Base() { }
    TmpMatrixView_(int m, int n) : Base(m,n) { }
    explicit TmpMatrixView_(const MatrixHelper<Scalar>& h) : Base(h) { }
    
    operator const Matrix_<ELT>&() const { return *reinterpret_cast<const Matrix_<ELT>*>(this); }
    operator Matrix_<ELT>&()             { return *reinterpret_cast<Matrix_<ELT>*>(this); }
    
    TmpMatrixView_* clone() const { return new TmpMatrixView_(*this); } 
private:
    // NO DATA MEMBERS ALLOWED
};

/**
 * This class is identical to a Matrix_; it is used only to manage the C++ rules
 * for when copy constructors are called by introducing a separate type to
 * prevent certain allowed optimizations from occuring when we don't want them.
 * Despite the name, this may be an owner if a Matrix_ is recast to a MatrixView_.
 * However, there are no owner constructors for MatrixView_. 
 */            
template <class ELT> class MatrixView_ : public MatrixBase<ELT> {
    typedef MatrixBase<ELT>                 Base;
    typedef typename CNT<ELT>::Scalar       S;
    typedef typename CNT<ELT>::StdNumber    StdNumber;
public:
    // Default construction is suppressed.
    // Uses default destructor.

    // Copy constructor is shallow. CAUTION: despite const argument, this preserves writability
    // if it was present in the source. This is necessary to allow temporary views to be
    // created and used as lvalues.
    MatrixView_(const MatrixView_& m) 
      : Base(const_cast<MatrixHelper<S>&>(m.helper), typename MatrixHelper<S>::ShallowCopy()) { }

    // Copy assignment is deep but not reallocating.
    MatrixView_& operator=(const MatrixView_& m) {
        Base::operator=(m); return *this;
    }

    // Ask for shallow copy    
    MatrixView_(const MatrixHelper<S>& h) : Base(h, typename MatrixHelper<S>::ShallowCopy()) { }
    MatrixView_(MatrixHelper<S>&       h) : Base(h, typename MatrixHelper<S>::ShallowCopy()) { }

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

    operator const Matrix_<ELT>&() const { return *reinterpret_cast<const Matrix_<ELT>*>(this); }
    operator Matrix_<ELT>&()             { return *reinterpret_cast<Matrix_<ELT>*>(this); }      

private:
    // NO DATA MEMBERS ALLOWED
    MatrixView_(); // default constructor suppressed (what's it a view of?)
};

template <class ELT> inline MatrixView_<ELT> 
MatrixBase<ELT>::block(int i, int j, int m, int n) const { 
    SimTK_INDEXCHECK(0,i,nrow()+1,"MatrixBase::block()");
    SimTK_INDEXCHECK(0,j,ncol()+1,"MatrixBase::block()");
    SimTK_SIZECHECK(i+m,nrow(),"MatrixBase::block()");
    SimTK_SIZECHECK(j+n,ncol(),"MatrixBase::block()");

    MatrixHelper<Scalar> h(helper,i,j,m,n);    
    return MatrixView_<ELT>(h); 
}
    
template <class ELT> inline MatrixView_<ELT>
MatrixBase<ELT>::updBlock(int i, int j, int m, int n) { 
    SimTK_INDEXCHECK(0,i,nrow()+1,"MatrixBase::updBlock()");
    SimTK_INDEXCHECK(0,j,ncol()+1,"MatrixBase::updBlock()");
    SimTK_SIZECHECK(i+m,nrow(),"MatrixBase::updBlock()");
    SimTK_SIZECHECK(j+n,ncol(),"MatrixBase::updBlock()");

    MatrixHelper<Scalar> h(helper,i,j,m,n);        
    return MatrixView_<ELT>(h); 
}

template <class E> inline MatrixView_<typename CNT<E>::THerm>
MatrixBase<E>::transpose() const { 
    MatrixHelper<typename CNT<Scalar>::THerm> 
        h(helper, typename MatrixHelper<typename CNT<Scalar>::THerm>::TransposeView());
    return MatrixView_<typename CNT<E>::THerm>(h); 
}
    
template <class E> inline MatrixView_<typename CNT<E>::THerm>
MatrixBase<E>::updTranspose() {     
    MatrixHelper<typename CNT<Scalar>::THerm> 
        h(helper, typename MatrixHelper<typename CNT<Scalar>::THerm>::TransposeView());
    return MatrixView_<typename CNT<E>::THerm>(h); 
}

template <class E> inline VectorView_<E>
MatrixBase<E>::diag() const { 
    MatrixHelper<Scalar> h(helper, typename MatrixHelper<Scalar>::DiagonalView());
    return VectorView_<E>(h); 
}
    
template <class E> inline VectorView_<E>
MatrixBase<E>::updDiag() {     
    MatrixHelper<Scalar> h(helper, typename MatrixHelper<Scalar>::DiagonalView());
    return VectorView_<E>(h);
}

template <class ELT> inline VectorView_<ELT> 
MatrixBase<ELT>::col(int j) const { 
    SimTK_INDEXCHECK(0,j,ncol(),"MatrixBase::col()");

    MatrixHelper<Scalar> h(helper,0,j,nrow(),1);    
    return VectorView_<ELT>(h); 
}
    
template <class ELT> inline VectorView_<ELT>
MatrixBase<ELT>::updCol(int j) {
    SimTK_INDEXCHECK(0,j,ncol(),"MatrixBase::updCol()");

    MatrixHelper<Scalar> h(helper,0,j,nrow(),1);        
    return VectorView_<ELT>(h); 
}

template <class ELT> inline RowVectorView_<ELT> 
MatrixBase<ELT>::row(int i) const { 
    SimTK_INDEXCHECK(0,i,nrow(),"MatrixBase::row()");

    MatrixHelper<Scalar> h(helper,i,0,1,ncol());    
    return RowVectorView_<ELT>(h); 
}
    
template <class ELT> inline RowVectorView_<ELT>
MatrixBase<ELT>::updRow(int i) { 
    SimTK_INDEXCHECK(0,i,nrow(),"MatrixBase::updRow()");

    MatrixHelper<Scalar> h(helper,i,0,1,ncol());        
    return RowVectorView_<ELT>(h); 
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
MatrixBase<ELT>::rowScale(const VectorBase<EE>& v, typename MatrixBase<ELT>::EltResult<EE>::Mul& out) const {
	assert(v.nrow() == nrow());
	out.resize(nrow(), ncol());
    for (int j=0; j<ncol(); ++j)
	    for (int i=0; i<nrow(); ++i)
		   out(i,j) = (*this)(i,j) * v[i];
}

// M = diag(v)^-1 * M; v must have nrow() elements.
// That is, M[i] /= v[i].
template <class ELT> template <class EE>  inline MatrixBase<ELT>& 
MatrixBase<ELT>::rowUnscaleInPlace(const VectorBase<EE>& v) {
	assert(v.nrow() == nrow());
	for (int i=0; i < nrow(); ++i)
		(*this)[i] /= v[i];
	return *this;
}

template <class ELT> template <class EE> inline void
MatrixBase<ELT>::rowUnscale(const VectorBase<EE>& v, typename MatrixBase<ELT>::EltResult<EE>::Dvd& out) const {
	assert(v.nrow() == nrow());
	out.resize(nrow(), ncol());
    for (int j=0; j<ncol(); ++j)
	    for (int i=0; i<nrow(); ++i)
		   out(i,j) = (*this)(i,j) / v[i];
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
MatrixBase<ELT>::colScale(const VectorBase<EE>& v, typename MatrixBase<ELT>::EltResult<EE>::Mul& out) const {
	assert(v.nrow() == ncol());
	out.resize(nrow(), ncol());
    for (int j=0; j<ncol(); ++j)
	    for (int i=0; i<nrow(); ++i)
		   out(i,j) = (*this)(i,j) * v[j];
}

// M = M * diag(v)^-1; v must have ncol() elements
// That is, M(i) /= v[i]
template <class ELT> template <class EE> inline MatrixBase<ELT>& 
MatrixBase<ELT>::colUnscaleInPlace(const VectorBase<EE>& v) {
    assert(v.nrow() == ncol());
	for (int j=0; j < ncol(); ++j)
		(*this)(j) /= v[j];
	return *this;
}

template <class ELT> template <class EE> inline void
MatrixBase<ELT>::colUnscale(const VectorBase<EE>& v, typename MatrixBase<ELT>::EltResult<EE>::Dvd& out) const {
	assert(v.nrow() == ncol());
	out.resize(nrow(), ncol());
    for (int j=0; j<ncol(); ++j)
	    for (int i=0; i<nrow(); ++i)
		   out(i,j) = (*this)(i,j) / v[j];
}

// M(i,j) *= R(i,j); R must have same dimensions as this
template <class ELT> template <class EE> inline MatrixBase<ELT>& 
MatrixBase<ELT>::elementwiseMultiplyInPlace(const MatrixBase<EE>& r) {
	assert(r.nrow()==nrow() && r.ncol()==ncol());
    for (int j=0; j<ncol(); ++j)
	    for (int i=0; i<nrow(); ++i)
			(*this)(i,j) *= r(i,j);
}

template <class ELT> template <class EE> inline void 
MatrixBase<ELT>::elementwiseMultiply(const MatrixBase<EE>& r, typename EltResult<EE>::Mul& out) const {
	assert(r.nrow()==nrow() && r.ncol()==ncol());
	out.resize(nrow(),ncol());
    for (int j=0; j<ncol(); ++j)
	    for (int i=0; i<nrow(); ++i)
			out(i,j) = (*this)(i,j) * r(i,j);
}

// M(i,j) /= R(i,j); R must have same dimensions as this
template <class ELT> template <class EE> inline MatrixBase<ELT>& 
MatrixBase<ELT>::elementwiseDivideInPlace(const MatrixBase<EE>& r) {
	assert(r.nrow()==nrow() && r.ncol()==ncol());
    for (int j=0; j<ncol(); ++j)
	    for (int i=0; i<nrow(); ++i)
			(*this)(i,j) /= r(i,j);
}

template <class ELT> template <class EE> inline void 
MatrixBase<ELT>::elementwiseDivide(const MatrixBase<EE>& r, typename EltResult<EE>::Dvd& out) const {
	assert(r.nrow()==nrow() && r.ncol()==ncol());
	out.resize(nrow(),ncol());
    for (int j=0; j<ncol(); ++j)
	    for (int i=0; i<nrow(); ++i)
			out(i,j) = (*this)(i,j) / r(i,j);
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


/**
 * This is the class intended to appear in user code. It can be a fixed-size view
 * of someone else's data, or can be a resizable data owner itself.
 */
template <class ELT> class Matrix_ : public MatrixBase<ELT> {
    typedef MatrixBase<ELT> Base;
    typedef typename CNT<ELT>::Scalar            S;
    typedef typename CNT<ELT>::Number            Number;
    typedef typename CNT<ELT>::StdNumber         StdNumber;

    typedef Matrix_<ELT>                         T;
    typedef MatrixView_<ELT>                     TView;
    typedef Matrix_< typename CNT<ELT>::TNeg >   TNeg;

public:
    Matrix_() : Base() { }

    // Copy constructor is deep.
    Matrix_(const Matrix_& src) : Base(src) { }

    // Assignment is a deep copy and will also allow reallocation if this Matrix
    // doesn't have a view.
    Matrix_& operator=(const Matrix_& src) { 
        Base::operator=(src); return *this;
    }

    // Force a deep copy of the view or whatever this is.    
    explicit Matrix_(const Base& v) : Base(v) { }   // e.g., MatrixView

    Matrix_(int m, int n) : Base(m,n) { }
    Matrix_(int m, int n, const ELT* initValsByRow) : Base(m,n,initValsByRow) { }
    Matrix_(int m, int n, const ELT& ival) : Base(m,n,ival) { }
    
    Matrix_(int m, int n, int leadingDim, const S* s): Base(m,n,leadingDim,s) { }
    Matrix_(int m, int n, int leadingDim,       S* s): Base(m,n,leadingDim,s) { }

    Matrix_& operator=(const ELT& v) { Base::operator=(v); return *this; }

    template <class EE> Matrix_& operator=(const MatrixBase<EE>& m)
      { Base::operator=(m); return*this; }
    template <class EE> Matrix_& operator+=(const MatrixBase<EE>& m)
      { Base::operator+=(m); return*this; }
    template <class EE> Matrix_& operator-=(const MatrixBase<EE>& m)
      { Base::operator-=(m); return*this; }

    Matrix_& operator*=(const StdNumber& t) { Base::operator*=(t); return *this; }
    Matrix_& operator/=(const StdNumber& t) { Base::operator/=(t); return *this; }

    const TNeg& negate()    const {return *reinterpret_cast<const TNeg*>(this); }
    TNeg&       updNegate()       {return *reinterpret_cast<TNeg*>(this); }

    const TNeg& operator-() const {return negate();}
    TNeg&       operator-()       {return updNegate();}
   
private:
    // NO DATA MEMBERS ALLOWED
};

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

    // Copy constructor is shallow. CAUTION: despite const argument, this preserves writability
    // if it was present in the source. This is necessary to allow temporary views to be
    // created and used as lvalues.
    VectorView_(const VectorView_& v) 
      : Base(const_cast<MatrixHelper<S>&>(v.helper), typename MatrixHelper<S>::ShallowCopy()) { }

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

private:
    // NO DATA MEMBERS ALLOWED
    VectorView_(); // default construction suppressed (what's it a View of?)
};

/**
 * XXX not ready for prime time
 */    
template <class ELT> class TmpVectorViewT : public VectorBase<ELT> {
    typedef VectorBase<ELT>             Base;
    typedef typename CNT<ELT>::Scalar   S;
public:
    TmpVectorViewT() : Base() { }
    explicit TmpVectorViewT(int m) : Base(m,1,false) { }
    explicit TmpVectorViewT(const MatrixHelper<S>& h) : Base(h) { }
    
    operator const Vector_<ELT>&() const { return *reinterpret_cast<const Vector_<ELT>*>(this); }
    operator Vector_<ELT>&()             { return *reinterpret_cast<Vector_<ELT>*>(this); }
    
    TmpVectorViewT* clone() const { return new TmpVectorViewT(*this); } 
private:
    // NO DATA MEMBERS ALLOWED
};

template <class ELT> class Vector_ : public VectorBase<ELT> {
    typedef VectorBase<ELT>                 Base;
    typedef typename CNT<ELT>::Scalar       S;
    typedef typename CNT<ELT>::Number       Number;
    typedef typename CNT<ELT>::StdNumber    StdNumber;
public:
    Vector_() : Base() { }  // 0x1 TODO, reallocatable
    // Uses default destructor.

    // Copy constructor is deep.
    Vector_(const Vector_& src) : Base(src) { }

    // Copy assignment is deep and can be reallocating if this Vector
    // has no View.
    Vector_& operator=(const Vector_& src) {
        Base::operator=(src); return*this;
    }

    explicit Vector_(const Base& src) : Base(src) { }    // e.g., VectorView

    explicit Vector_(int m) : Base(m,false) { }
    Vector_(int m, const ELT* initVals) : Base(m,false,initVals) { }
    Vector_(int m, const ELT& ival) : Base(m,false,ival) { }

    // Construct a Vector which uses borrowed space.
    // Last parameter is a dummy to avoid overload conflicts when ELT=S.    
    Vector_(int m, const S* s, bool): Base(m,m,s) { }
    Vector_(int m,       S* s, bool): Base(m,m,s) { }
    
    Vector_& operator=(const ELT& v) { Base::operator=(v); return *this; } 

    template <class EE> Vector_& operator=(const VectorBase<EE>& m)
      { Base::operator=(m); return*this; }
    template <class EE> Vector_& operator+=(const VectorBase<EE>& m)
      { Base::operator+=(m); return*this; }
    template <class EE> Vector_& operator-=(const VectorBase<EE>& m)
      { Base::operator-=(m); return*this; }

    Vector_& operator*=(const StdNumber& t) { Base::operator*=(t); return *this; }
    Vector_& operator/=(const StdNumber& t) { Base::operator/=(t); return *this; }
 
private:
    // NO DATA MEMBERS ALLOWED
};


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

    // Copy constructor is shallow. CAUTION: despite const argument, this preserves writability
    // if it was present in the source. This is necessary to allow temporary views to be
    // created and used as lvalues.
    RowVectorView_(const RowVectorView_& r) 
      : Base(const_cast<MatrixHelper<S>&>(r.helper), typename MatrixHelper<S>::ShallowCopy()) { }

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

private:
    // NO DATA MEMBERS ALLOWED
    RowVectorView_(); // default construction suppressed (what is it a view of?)
};

/**
 * XXX not ready for prime time
 */    
template <class ELT> class TmpRowVectorView_ : public RowVectorBase<ELT> {
    typedef RowVectorBase<ELT>          Base;
    typedef typename CNT<ELT>::Scalar   S;
public:
    TmpRowVectorView_() : Base() { }
    TmpRowVectorView_(int n) : Base(n,false) { }
    TmpRowVectorView_(const MatrixHelper<S>& h) : Base(h) { }
    
    operator const RowVector_<ELT>&() const { return *reinterpret_cast<const RowVector_<ELT>*>(this); }
    operator RowVector_<ELT>&()             { return *reinterpret_cast<RowVector_<ELT>*>(this); }
    
    TmpRowVectorView_* clone() const { return new TmpRowVectorView_(*this); } 
private:
    // NO DATA MEMBERS ALLOWED
};


template <class ELT> class RowVector_ : public RowVectorBase<ELT> {
    typedef RowVectorBase<ELT>              Base;
    typedef typename CNT<ELT>::Scalar       S;
    typedef typename CNT<ELT>::Number       Number;
    typedef typename CNT<ELT>::StdNumber    StdNumber;
public:
    RowVector_() : Base() { }   // 1x0 TODO, reallocatable
    // Uses default destructor.

    // Copy constructor is deep.
    RowVector_(const RowVector_& src) : Base(src) { }

    // Copy assignment is deep and can be reallocating if this RowVector
    // has no View.
    RowVector_& operator=(const RowVector_& src) {
        Base::operator=(src); return*this;
    }

    explicit RowVector_(const Base& src) : Base(src) { }    // e.g., RowVectorView

    explicit RowVector_(int n) : Base(n,false) { }
    RowVector_(int n, const ELT* initVals) : Base(n,false,initVals) { }
    RowVector_(int n, const ELT& ival)     : Base(n,false,ival) { }

    // Construct a RowVector which uses borrowed space.
    // Last parameter is a dummy to avoid overload conflicts when ELT=S.    
    RowVector_(int n, const S* s, bool): Base(n,1,s) { }
    RowVector_(int n,       S* s, bool): Base(n,1,s) { }
    
    RowVector_& operator=(const ELT& v) { Base::operator=(v); return *this; } 

    template <class EE> RowVector_& operator=(const RowVectorBase<EE>& b)
      { Base::operator=(b); return*this; }
    template <class EE> RowVector_& operator+=(const RowVectorBase<EE>& b)
      { Base::operator+=(b); return*this; }
    template <class EE> RowVector_& operator-=(const RowVectorBase<EE>& b)
      { Base::operator-=(b); return*this; }

    RowVector_& operator*=(const StdNumber& t) { Base::operator*=(t); return *this; }
    RowVector_& operator/=(const StdNumber& t) { Base::operator/=(t); return *this; }
 
private:
    // NO DATA MEMBERS ALLOWED
};

    // GLOBAL OPERATORS: Matrix_

// + and - allow mixed element types, but will fail to compile if the elements aren't
// compatible. At run time these will fail if the dimensions are incompatible.
template <class E1, class E2>
Matrix_<typename CNT<E1>::template Result<E2>::Add>
operator+(const Matrix_<E1>& l, const Matrix_<E2>& r) {
    return Matrix_<typename CNT<E1>::template Result<E2>::Add>(l) += r;
}
template <class E1, class E2>
Matrix_<typename CNT<E1>::template Result<E2>::Sub>
operator-(const Matrix_<E1>& l, const Matrix_<E2>& r) {
    return Matrix_<typename CNT<E1>::template Result<E2>::Sub>(l) -= r;
}

// Scalar multiply and divide. You might wish the scalar could be
// a templatized type "E2", but that would create horrible ambiguities since
// E2 would match not only scalar types but everything else including
// matrices.
template <class E> Matrix_<E>
operator*(const Matrix_<E>& l, const typename CNT<E>::StdNumber& r) 
  { return Matrix_<E>(l)*=r; }

template <class E> Matrix_<E>
operator*(const typename CNT<E>::StdNumber& l, const Matrix_<E>& r) 
  { return Matrix_<E>(r)*=l; }

template <class E> Matrix_<E>
operator/(const Matrix_<E>& l, const typename CNT<E>::StdNumber& r) 
  { return Matrix_<E>(l)/=r; }

    // GLOBAL OPERATORS: Vector_

template <class E1, class E2>
Vector_<typename CNT<E1>::template Result<E2>::Add>
operator+(const Vector_<E1>& l, const Vector_<E2>& r) {
    return Vector_<typename CNT<E1>::template Result<E2>::Add>(l) += r;
}
template <class E1, class E2>
Vector_<typename CNT<E1>::template Result<E2>::Sub>
operator-(const Vector_<E1>& l, const Vector_<E2>& r) {
    return Vector_<typename CNT<E1>::template Result<E2>::Sub>(l) -= r;
}

template <class E> Vector_<E>
operator*(const Vector_<E>& l, const typename CNT<E>::StdNumber& r) 
  { return Vector_<E>(l)*=r; }

template <class E> Vector_<E>
operator*(const typename CNT<E>::StdNumber& l, const Vector_<E>& r) 
  { return Vector_<E>(r)*=l; }

template <class E> Vector_<E>
operator/(const Vector_<E>& l, const typename CNT<E>::StdNumber& r) 
  { return Vector_<E>(l)/=r; }

    // GLOBAL OPERATORS: RowVector_

template <class E1, class E2>
RowVector_<typename CNT<E1>::template Result<E2>::Add>
operator+(const RowVector_<E1>& l, const RowVector_<E2>& r) {
    return RowVector_<typename CNT<E1>::template Result<E2>::Add>(l) += r;
}
template <class E1, class E2>
RowVector_<typename CNT<E1>::template Result<E2>::Sub>
operator-(const RowVector_<E1>& l, const RowVector_<E2>& r) {
    return RowVector_<typename CNT<E1>::template Result<E2>::Sub>(l) -= r;
}

template <class E> RowVector_<E>
operator*(const RowVector_<E>& l, const typename CNT<E>::StdNumber& r) 
  { return RowVector_<E>(l)*=r; }

template <class E> RowVector_<E>
operator*(const typename CNT<E>::StdNumber& l, const RowVector_<E>& r) 
  { return RowVector_<E>(r)*=l; }

template <class E> RowVector_<E>
operator/(const RowVector_<E>& l, const typename CNT<E>::StdNumber& r) 
  { return RowVector_<E>(l)/=r; }


    // GLOBAL OPERATORS: mixed

    // TODO: these should use LAPACK!

// Dot product
template <class E1, class E2> 
typename CNT<E1>::template Result<E2>::Mul
operator*(const RowVector_<E1>& r, const Vector_<E2>& v) {
    assert(r.ncol() == v.nrow());
    typename CNT<E1>::template Result<E2>::Mul sum(0);
    for (int j=0; j < r.ncol(); ++j)
        sum += r(j) * v[j];
    return sum;
}
template <class E1, class E2> 
typename CNT<E1>::template Result<E2>::Mul
operator*(const RowVectorView_<E1>& r, const Vector_<E2>& v) {
    return r.getAsRowVector()*v;
}
template <class E1, class E2> 
typename CNT<E1>::template Result<E2>::Mul
operator*(const RowVector_<E1>& r, const VectorView_<E2>& v) {
    return r*v.getAsVector();
}
template <class E1, class E2> 
typename CNT<E1>::template Result<E2>::Mul
operator*(const RowVectorView_<E1>& r, const VectorView_<E2>& v) {
    return r.getAsRowVector()*v.getAsVector();
}

template <class E1, class E2> 
Vector_<typename CNT<E1>::template Result<E2>::Mul>
operator*(const Matrix_<E1>& m, const Vector_<E2>& v) {
    assert(m.ncol() == v.nrow());
    Vector_<typename CNT<E1>::template Result<E2>::Mul> res(m.nrow());
    for (int i=0; i< m.nrow(); ++i)
        res[i] = m[i]*v;
    return res;
}
template <class E1, class E2> 
Vector_<typename CNT<E1>::template Result<E2>::Mul>
operator*(const MatrixView_<E1>& m, const Vector_<E2>& v) {
    return m.getAsMatrix()*v;
}
template <class E1, class E2> 
Vector_<typename CNT<E1>::template Result<E2>::Mul>
operator*(const Matrix_<E1>& m, const VectorView_<E2>& v) {
    return m*v.getAsVector();
}
template <class E1, class E2> 
Vector_<typename CNT<E1>::template Result<E2>::Mul>
operator*(const MatrixView_<E1>& m, const VectorView_<E2>& v) {
    return m.getAsMatrix()*v.getAsVector();
}

template <class E1, class E2> 
Matrix_<typename CNT<E1>::template Result<E2>::Mul>
operator*(const Matrix_<E1>& m1, const Matrix_<E2>& m2) {
    assert(m1.ncol() == m2.nrow());
    Matrix_<typename CNT<E1>::template Result<E2>::Mul> 
        res(m1.nrow(),m2.ncol());

    for (int j=0; j < res.ncol(); ++j)
        for (int i=0; i < res.nrow(); ++i)
            res(i,j) = m1[i] * m2(j);

    return res;
}
template <class E1, class E2> 
Matrix_<typename CNT<E1>::template Result<E2>::Mul>
operator*(const MatrixView_<E1>& m1, const Matrix_<E2>& m2) {
    return m1.getAsMatrix()*m2;
}
template <class E1, class E2> 
Matrix_<typename CNT<E1>::template Result<E2>::Mul>
operator*(const Matrix_<E1>& m1, const MatrixView_<E2>& m2) {
    return m1*m2.getAsMatrix();
}
template <class E1, class E2> 
Matrix_<typename CNT<E1>::template Result<E2>::Mul>
operator*(const MatrixView_<E1>& m1, const MatrixView_<E2>& m2) {
    return m1.getAsMatrix()*m2.getAsMatrix();
}
    // GLOBAL OPERATORS: I/O

template <class T> inline std::ostream&
operator<<(std::ostream& o, const VectorBase<T>& v)
{ o << "~["; for(int i=0;i<v.size();++i) o<<(i>0?" ":"")<<v[i]; 
    return o << "]"; }

template <class T> inline std::ostream&
operator<<(std::ostream& o, const RowVectorBase<T>& v)
{ o << "["; for(int i=0;i<v.size();++i) o<<(i>0?" ":"")<<v[i]; 
    return o << "]"; }

template <class T> inline std::ostream&
operator<<(std::ostream& o, const MatrixBase<T>& m) {
    for (int i=0;i<m.nrow();++i)
        o << std::endl << m[i];
    if (m.nrow()) o << std::endl;
    return o; 
}


// Friendly abbreviations
typedef Vector_<Real>                           Vector;     // default precision
typedef Vector_<float>                          FVector;
typedef Vector_<double>                         DVector;
typedef Vector_<long double>                    LVector;

typedef Vector_<Complex>                        CVector;    // default precision
typedef Vector_<std::complex<float> >           FCVector;
typedef Vector_<std::complex<double> >          DCVector;
typedef Vector_<std::complex<long double> >     LCVector;

typedef VectorView_<Real>                       VectorView;
typedef VectorView_<float>                      FVectorView;
typedef VectorView_<double>                     DVectorView;
typedef VectorView_<long double>                LVectorView;

typedef VectorView_<Complex>                    CVectorView;
typedef VectorView_<std::complex<float> >       FCVectorView;
typedef VectorView_<std::complex<double> >      DCVectorView;
typedef VectorView_<std::complex<long double> > LCVectorView;

typedef RowVector_<Real>                           RowVector;     // default precision
typedef RowVector_<float>                          FRowVector;
typedef RowVector_<double>                         DRowVector;
typedef RowVector_<long double>                    LRowVector;

typedef RowVector_<Complex>                        CRowVector;    // default precision
typedef RowVector_<std::complex<float> >           FCRowVector;
typedef RowVector_<std::complex<double> >          DCRowVector;
typedef RowVector_<std::complex<long double> >     LCRowVector;

typedef RowVectorView_<Real>                       RowVectorView;
typedef RowVectorView_<float>                      FRowVectorView;
typedef RowVectorView_<double>                     DRowVectorView;
typedef RowVectorView_<long double>                LRowVectorView;

typedef RowVectorView_<Complex>                    CRowVectorView;
typedef RowVectorView_<std::complex<float> >       FCRowVectorView;
typedef RowVectorView_<std::complex<double> >      DCRowVectorView;
typedef RowVectorView_<std::complex<long double> > LCRowVectorView;

typedef Matrix_<Real>                           Matrix;
typedef Matrix_<float>                          FMatrix;
typedef Matrix_<double>                         DMatrix;
typedef Matrix_<long double>                    LMatrix;

typedef Matrix_<Complex>                        CMatrix;
typedef Matrix_<std::complex<float> >           FCMatrix;
typedef Matrix_<std::complex<double> >          DCMatrix;
typedef Matrix_<std::complex<long double> >     LCMatrix;

typedef MatrixView_<Real>                       MatrixView;
typedef MatrixView_<float>                      FMatrixView;
typedef MatrixView_<double>                     DMatrixView;
typedef MatrixView_<long double>                LMatrixView;

typedef MatrixView_<Complex>                    CMatrixView;
typedef MatrixView_<std::complex<float> >       FCMatrixView;
typedef MatrixView_<std::complex<double> >      DCMatrixView;
typedef MatrixView_<std::complex<long double> > LCMatrixView;


    
// Not all combinations of structure/sparsity/storage/condition
// are allowed. 

class MatrixShape {
public:

    MatrixShape() : shape(MatrixShapes::Uncommitted) { }
    // implicit conversion
    MatrixShape(MatrixShapes::Shape s) : shape(s) { }
    MatrixShapes::Shape shape;
};

class MatrixSize {
public:
    enum Freedom {
        Uncommitted = 0x00,    // nrow, ncol variable
        FixedNRows  = 0x01,    // can't vary nrows
        FixedNCols  = 0x02,    // can't vary ncols
        Fixed       = FixedNRows | FixedNCols
    };

    MatrixSize() 
      : freedom(Uncommitted), nrow(0), ncol(0) { }

    MatrixSize(Freedom f, long nr, long nc) 
      : freedom(f), nrow(nr), ncol(nc) { }

    Freedom freedom;
    long nrow;
    long ncol;
};

class MatrixStructure {
public:

    MatrixStructure() : structure(MatrixStructures::Uncommitted) { }
    
    // implicit conversion
    MatrixStructure(MatrixStructures::Structure ms) : structure(ms) { }

    MatrixStructures::Structure structure;
};

class MatrixSparsity {
public:

    // If Banded, how stuck are we on a particular bandwidth?
    enum BandwidthFreedom {
        Free        = 0x00,    // upper & lower are free
        FixedUpper  = 0x01,
        FixedLower  = 0x02,
        Fixed       = FixedUpper | FixedLower
    };

    MatrixSparsity() 
      : sparsity(MatrixSparseFormats::Uncommitted), lowerBandwidth(-1), upperBandwidth(-1)
    {
    }

    MatrixSparsity(int lower, int upper) 
      : sparsity(MatrixSparseFormats::Banded), lowerBandwidth(lower), upperBandwidth(upper)
    {
        assert(lower >= 0 && upper >= 0);
    }

    MatrixSparseFormats::Sparsity sparsity;
    int lowerBandwidth;
    int upperBandwidth;
};

class MatrixStorage {
public:

    enum Position {
        Lower,              // lower is default
        Upper
    };

    // OR-able
    enum Assumptions {
        None         = 0x00,
        UnitDiagonal = 0x01
    };

    MatrixStorage() 
      : storage(MatrixStorageFormats::Uncommitted), position(Lower), assumptions(None)
    {
    }
    
    // also serves as implicit conversion from Storage type
    MatrixStorage(MatrixStorageFormats::Storage s, Position p=Lower, Assumptions a=None) 
        : storage(s), position(p), assumptions(a)
    {
    }

    MatrixStorageFormats::Storage     storage;
    Position    position;
    Assumptions assumptions;

    // All the 2d formats allow a leading dimension larger
    // than the number of rows, producing a gap between
    // each column.
    int leadingDimension;

    // 1d formats allow spacing between elements. Stride==1
    // means packed.
    int stride;
};

class MatrixCondition {
public:

    MatrixCondition() : condition(MatrixConditions::Uncommitted) { }
    // implicit conversion
    MatrixCondition(MatrixConditions::Condition c) : condition(c) { }

    MatrixConditions::Condition condition;
};
} //namespace SimTK

#endif //SimTK_SIMMATRIX_BIGMATRIX_H_
