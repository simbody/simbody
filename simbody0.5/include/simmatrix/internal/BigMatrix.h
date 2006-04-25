#ifndef _SimTK_SIMMATRIX_BIGMATRIX_H_
#define _SimTK_SIMMATRIX_BIGMATRIX_H_

/* Copyright (c) 2005-6 Stanford University and Michael Sherman.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "SimTKcommon.h"
#include "simmatrix/internal/common.h"
#include "simmatrix/internal/Scalar.h"
#include "simmatrix/internal/SmallMatrix.h"

#include <iostream>
#include <cassert>
#include <complex>
#include <cstddef>

namespace simtkimpl {
class MatrixViewImpl;
template <class S> class MatrixDataImpl;  
}

namespace SimTK {

template <class ELT>    class MatrixBase;

template <class T=Real> class MatrixView_;
template <class T=Real> class TmpMatrixView_;
template <class T=Real> class Matrix_;

template <class T=Real> class VectorView_;
template <class T=Real> class TmpVectorView_;
template <class T=Real> class Vector_;

template <class T=Real> class RowVectorView_;
template <class T=Real> class TmpRowVectorView_;
template <class T=Real> class RowVector_;

template <class S> class TmpMatrixHelper;

/**
 * Scalar-type templatized helper class for the more general, numerical-type
 * templatized class MatrixBase<ELT>. This will be instantiated once each for 
 * float, double, long double and the associated complex and conjugate types,
 * and their Negators (a total of 18 types). Element size is
 * dealt with at run time; otherwise the helper knows nothing about the 
 * structure of the elements.
 * 
 * The constructors for numerical types should not initialize the numerical
 * values. We take advantage of that here -- this class assumes it can simply
 * allocate the appropriate amount of data as an array of the underlying
 * scalar type, with no implicit initialization. We'll fill uninitialized data
 * with NaNs when debugging or if specifically requested; otherwise it will
 * contain garbage initially.
 * 
 * Note that the implementation here is critically opaque: there is a 32 bit
 * element size, a 32 bit flags word, and two pointers (view & data) which
 * will be 32 or 64 bits depending on the machine we're on. That's it, and
 * that is exactly how it will be until hell freezes over. The meaning of
 * the flags and what the pointers point to will change frequently until
 * perfection is achieved! However, binary compatibility will be maintained 
 * regardless.
 */
template <class S> 
class SimTK_SIMBODY_API MatrixHelper {
    typedef typename CNT<S>::Number     Number;     // strips off negator from S
    typedef typename CNT<S>::StdNumber  StdNumber;  // turns conjugate to complex
    typedef typename CNT<S>::Precision  Precision;  // strips off complex from StdNumber
public:     
    // no default constructor, copy constructor suppressed

    // Local types for directing us to the right constructor at compile time.
    class ShallowCopy { };
    class DeepCopy    { };
    class TransposeView   { };
    
    // Matrix owner constructors
    
    // 0x0, fully resizable
    explicit MatrixHelper(unsigned int esz) : eltSize(esz),flags(0),view(0),data(0) {
        SimTK_ASSERT1(esz > 0, "MatrixHelper constructor was given a bad element size %u", esz); 
    }
    
    // (m*esz) x n, resizable with optional restrictions 
    MatrixHelper(unsigned int esz, int m, int n, bool lockNrow=false, bool lockNcol=false);
    
    // This is a non-resizeable owner (of the data descriptor) but shares pre-existing raw data.
    MatrixHelper(unsigned int esz, int m, int n, int leadingDim, const S*); // read only
    MatrixHelper(unsigned int esz, int m, int n, int leadingDim, S*);       // writable
   
    // Matrix view constructors, read only and writable. Use these for submatrices,
    // rows, and columns. Indices are by *element* and these may consist of multiple
    // scalars of type template parameter S.
    MatrixHelper(const MatrixHelper&, int i, int j, int nrow, int ncol);
    MatrixHelper(MatrixHelper&,       int i, int j, int nrow, int ncol); 

    MatrixHelper(const MatrixHelper<typename CNT<S>::THerm>&, 
                 const TransposeView&);    // a read only transposed view
    MatrixHelper(MatrixHelper<typename CNT<S>::THerm>&,       
                 const TransposeView&);    // a writable transposed view

    // Copy an existing MatrixHelper using deep or shallow copy as requested. 
    // For shallow copy, const form loses writability, non-const retains same
    // writability status as source.    
    MatrixHelper(const MatrixHelper&, const ShallowCopy&);
    MatrixHelper(MatrixHelper&,       const ShallowCopy&);
    MatrixHelper(const MatrixHelper&, const DeepCopy&);

    MatrixHelper(const MatrixHelper<typename CNT<S>::TNeg>&, const DeepCopy&);

    // Transfer view and data from source; always a deep copy and practically free.
    // The source helper is left empty with neither view nor data.    
    inline explicit MatrixHelper(TmpMatrixHelper<S>&);

    
    // Behavior of assignment depends on whether this is an owner or view. If
    // it's an owner it is resized and ends up a new, dense copy of rhs. If
    // it's a view, sizes must match and we copy elements from rhs to corresponding
    // elements of lhs.        
    MatrixHelper& operator=(const MatrixHelper&); 
    MatrixHelper& operator=(const MatrixHelper<typename CNT<S>::TNeg>&);


    // Using *element* indices, obtain a pointer to the beginning of a particular
    // element. This is always a slow operation compared to raw array access;
    // use sparingly.    
    const S* getElt(int i, int j) const;
    S*       updElt(int i, int j);

    
    // Add up all the *elements*, returning the answer in the element given
    // by pointer to its first scalar.
    void sum(S* eltp) const;
    void colSum(int j, S* eltp) const;
    void rowSum(int i, S* eltp) const;
        
    // addition and subtraction (+= and -=)
    void addIn(const MatrixHelper&);   
    void addIn(const MatrixHelper<typename CNT<S>::TNeg>&);   
    void subIn(const MatrixHelper&); 
    void subIn(const MatrixHelper<typename CNT<S>::TNeg>&); 
    
    void fillWith(const S* eltp);
    void copyInByRows(const S* elts);
    
        // Scalar operations //

    // Fill every element with repeated copies of a single scalar value.
    void fillWithScalar(const StdNumber&);
            
    // Scalar multiply (and divide). This is useful no matter what the
    // element structure and will produce the correct result.
    void scaleBy(const StdNumber&);
    
    // Sums of scalars, regardless of element structure. Not much use except
    // when the elements are themselves scalars.
    S scalarSum() const;
    S scalarColSum(int j) const;
    S scalarRowSum(int i) const; 

    // This is only allowed for a matrix of real or complex or neg of those,
    // which is square, well-conditioned, and for which we have no view,
    // and element size 1.
    void invertInPlace();

    // See comment in MatrixBase::matmul for an explanation.
    template <class SA, class SB>
    void matmul(const StdNumber& beta,   // applied to 'this'
                const StdNumber& alpha, const MatrixHelper<SA>& A, const MatrixHelper<SB>& B);
    
    // Here we use element indices but return the first scalar we find there.
    // This is useless except when the elements are scalars.
    const S& getScalar(int i, int j) const;
    S&       updScalar(int i, int j);      
    
        // Bookkeeping //
        
    int nrow() const;
    int ncol() const;
    ptrdiff_t size() const; // nrow*ncol        
    void resize(int m, int n);
    ~MatrixHelper();


protected:
    const unsigned int            eltSize;  // # scalars per element
    unsigned int                  flags;    // meaning is opaque
    simtkimpl::MatrixViewImpl*    view;     // details are opaque
    simtkimpl::MatrixDataImpl<S>* data;

    MatrixHelper(unsigned int esz, unsigned int f, 
                 simtkimpl::MatrixViewImpl* v, simtkimpl::MatrixDataImpl<S>* d)
      : eltSize(esz), flags(f), view(v), data(d) {
    }

    template <class T> friend class MatrixHelper;

private:
    unsigned int getElementSize() const { return eltSize; }
    bool isDataOwned() const;
        
    // Suppress copy constructor.
    MatrixHelper(const MatrixHelper&);
};

/**
 * This is just a MatrixHelper which is about to be thrown away. We can
 * thus safely steal the pointers out of it rather than make copies. Saves
 * heap allocation.
 */
template <class S> class TmpMatrixHelper : public MatrixHelper<S> {
public:
    // view construction only
    
    // read only temporary view
    TmpMatrixHelper(const MatrixHelper<S>& h, int i, int j, int nrow, int ncol)
        : MatrixHelper<S>(h,i,j,nrow,ncol) { }

    // writable temporary view        
    TmpMatrixHelper(MatrixHelper<S>& h,int i, int j, int nrow, int ncol)   
        : MatrixHelper<S>(h,i,j,nrow,ncol) { }

    // copy constructor steals space; costs nothing            
    TmpMatrixHelper(TmpMatrixHelper& th) : MatrixHelper<S>(th.eltSize) 
    	{ th.steal(this->view,this->data); }    
  
    // default destructor
private:
    // NO DATA MEMBERS ALLOWED
    friend class MatrixHelper<S>;   
    void steal(simtkimpl::MatrixViewImpl*& v, simtkimpl::MatrixDataImpl<S>*& d) 
    	{ v=this->view; d=this->data; this->view=0; this->data=0; }
    TmpMatrixHelper& operator=(const TmpMatrixHelper&); // suppress
};

template <class S> inline
MatrixHelper<S>::MatrixHelper(TmpMatrixHelper<S>& th) 
    : eltSize(th.eltSize), flags(th.flags) { th.steal(view,data); }    
 
/**
 * Variable-size 2d matrix of Composite Numerical Type (ELT) elements. This is
 * a container of such elements, it is NOT a Composite Numerical Type itself.
 * MatrixBase<ELT> uses MatrixHelper<S> for implementation, where S is 
 * ELT::Scalar, that is, the underlying float, double, long double,
 * complex<float>, negator<conjugate<long double>>, 
 * etc. from which ELT is constructed (this is a finite set).
 * 
 * MatrixBase is the only class in the Matrix/Vector family which has any
 * data members (it has exactly one MatrixHelper). Thus all other objects
 * in this family are exactly the same size and may be "reinterpreted" as
 * appropriate. For example, a Vector may be reinterpreted as a Matrix or
 * vice versa, provided runtime requirements are met.
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

    typedef MatrixBase<E>                TRow;  // different shape, of course!
    typedef MatrixBase<E>                TCol;
    typedef MatrixBase<E>                TDiag;

    // This gives the resulting matrix type when (m(i,j) op P) is applied to each element.
    // It will have element types which are the regular composite result of E op P.
    template <class P> struct EltResult { 
        typedef MatrixBase<typename CNT<E>::template Result<P>::Mul> Mul;
        typedef MatrixBase<typename CNT<E>::template Result<P>::Dvd> Dvd;
        typedef MatrixBase<typename CNT<E>::template Result<P>::Add> Add;
        typedef MatrixBase<typename CNT<E>::template Result<P>::Sib> Sub;
    };
       
    ptrdiff_t size() const { return helper.size(); }
    int       nrow() const { return helper.nrow(); }
    int       ncol() const { return helper.ncol(); }

    // Scalar norm square is sum( squares of all scalars )
    // TODO: very slow! Should be optimized at least for the case
    //       where ELT is a Scalar.
    ScalarSq scalarNormSqr() const { 
        ScalarSq sum(0);
        for(int j=0;j<ncol();++j) 
            for (int i=0; i<nrow(); ++i)
                sum += CNT<E>::scalarNormSqr((*this)(i,j));
        return sum;
    }

    // abs() is elementwise absolute value; that is, the return value has the same
    // dimension as this Matrix but with each element replaced by whatever it thinks
    // its absolute value is.
    // TODO: very slow! Should be optimized at least for the case
    //       where ELT is a Scalar.
    TAbs abs() const { 
        TAbs mabs;
        for(int j=0;j<ncol();++j) 
            for (int i=0; i<nrow(); ++i)
                mabs(i,j) = CNT<E>::abs((*this)(i,j));
        return mabs;
    }

    enum { 
        NScalarsPerElement = CNT<E>::NActualScalars
    };
  
    MatrixBase() : helper(NScalarsPerElement)  { }
    
    // Copy constructor is a deep copy (not appropriate for views!).    
    MatrixBase(const MatrixBase& b) : helper(b.helper, typename MatrixHelper<Scalar>::DeepCopy()) { }
    
    // Copy assignment is a deep copy but behavior depends on type of lhs: if view, rhs
    // must match. If owner, we reallocate and copy rhs.
    MatrixBase& operator=(const MatrixBase& b) { helper=b.helper; return *this; }

    // default destructor


    MatrixBase(int m, int n, bool lockNrow=false, bool lockNcol=false) 
      : helper(NScalarsPerElement, m, n, lockNrow, lockNcol)  { }
        
    // Non-resizeable owner sharing pre-allocated raw data
    MatrixBase(int m, int n, int leadingDim, const Scalar* s)
      : helper(NScalarsPerElement, m, n, leadingDim, s) { }    // read only
    MatrixBase(int m, int n, int leadingDim, Scalar* s)
      : helper(NScalarsPerElement, m, n, leadingDim, s) { }    // writable
            
    MatrixBase(int m, int n, const ELT& t) : helper(NScalarsPerElement, m, n)
      { helper.fillWith(reinterpret_cast<const Scalar*>(&t)); }    
    MatrixBase(int m, int n, bool lockNrow, bool lockNcol, const ELT& t) 
      : helper(NScalarsPerElement, m, n, lockNrow, lockNcol)
      { helper.fillWith(reinterpret_cast<const Scalar*>(&t)); }
        
    MatrixBase(int m, int n, const ELT* p) : helper(NScalarsPerElement, m, n)
      { helper.copyInByRows(reinterpret_cast<const Scalar*>(p)); }
    MatrixBase(int m, int n,  bool lockNrow, bool lockNcol, const ELT* p) 
      : helper(NScalarsPerElement, m, n, lockNrow, lockNcol)
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

    // fill current allocation with copies of element
    MatrixBase& operator=(const ELT& t) 
        { helper.fillWith(reinterpret_cast<const Scalar*>(&t)); return *this; }    

    void setToNaN() {helper.fillWithScalar(CNT<StdNumber>::getNaN());}
    void setToZero() {helper.fillWithScalar(StdNumber(0));}
 
    // View creating operators. TODO: these should be TmpMatrixViews  
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

    // TODO: this routine seems ill-advised but I need it for the IVM port at the moment
    Matrix_<StdNumber> invert() const {  // return a newly-allocated inverse; dump negator 
        Matrix_<StdNumber> m(*this);
        m.helper.invertInPlace();
        return m;   // TODO - bad: makes an extra copy
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
    // NOTE: neither A nor B can be the same matrix as 'this', or views of the same data
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
    ScalarSq norm()    const { return std::sqrt(scalarNormSqr()); } // TODO -- not good

    ELT sum() const { ELT s; helper.sum(reinterpret_cast<Scalar*>(&s)); return s; } // add all the elements        

    //TODO: make unary +/- return a self-view so they won't reallocate?
    const MatrixBase& operator+() const {return *this; }
    const TNeg&       negate()    const {return *reinterpret_cast<const TNeg*>(this); }
    TNeg&             updNegate()       {return *reinterpret_cast<TNeg*>(this); }

    const TNeg&       operator-() const {return negate();}
    TNeg&             operator-()       {return updNegate();}
 
    void resize(int m, int n) { helper.resize(m,n); }
    
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
                
    const RowVectorView_<ELT>& getAsRowVectorView() const 
      { assert(nrow()==1); return *reinterpret_cast<const RowVectorView_<ELT>*>(this); }
    RowVectorView_<ELT>&       updAsRowVectorView()       
      { assert(nrow()==1); return *reinterpret_cast<      RowVectorView_<ELT>*>(this); } 
    const RowVector_<ELT>&     getAsRowVector()     const 
      { assert(nrow()==1); return *reinterpret_cast<const RowVector_<ELT>*>(this); }
    RowVector_<ELT>&           updAsRowVector()           
      { assert(nrow()==1); return *reinterpret_cast<      RowVector_<ELT>*>(this); }        
       
protected:
    MatrixHelper<Scalar>    helper;

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
    typedef VectorBase< typename CNT<ELT>::TNeg >       TNeg;
    typedef RowVectorView_< typename CNT<ELT>::THerm >  THerm;
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


    // fill current allocation with copies of element
    VectorBase& operator=(const ELT& t) { Base::operator=(t); return *this; }  

    // Implicit conversions are allowed to Vector or Matrix, but not to RowVector.   
    operator const Vector_<ELT>&()     const { return *reinterpret_cast<const Vector_<ELT>*>(this); }
    operator       Vector_<ELT>&()           { return *reinterpret_cast<      Vector_<ELT>*>(this); }
    operator const VectorView_<ELT>&() const { return *reinterpret_cast<const VectorView_<ELT>*>(this); }
    operator       VectorView_<ELT>&()       { return *reinterpret_cast<      VectorView_<ELT>*>(this); }
    
    operator const Matrix_<ELT>&()     const { return *reinterpret_cast<const Matrix_<ELT>*>(this); }
    operator       Matrix_<ELT>&()           { return *reinterpret_cast<      Matrix_<ELT>*>(this); } 
    operator const MatrixView_<ELT>&() const { return *reinterpret_cast<const MatrixView_<ELT>*>(this); }
    operator       MatrixView_<ELT>&()       { return *reinterpret_cast<      MatrixView_<ELT>*>(this); } 
    
    // Override MatrixBase indexing operators          
    const ELT& operator[](int i) const { return Base::operator()(i,0); }
    ELT&       operator[](int i)       { return Base::operator()(i,0); }
    const ELT& operator()(int i) const { return Base::operator()(i,0); }
    ELT&       operator()(int i)       { return Base::operator()(i,0); }
         
    // View creation      
    VectorView_<ELT> operator()(int i, int m) const { return Base::operator()(i,0,m,1).getAsVectorView(); }
    VectorView_<ELT> operator()(int i, int m)       { return Base::operator()(i,0,m,1).updAsVectorView(); }
 
    // Hermitian transpose.
    THerm transpose() const {return Base::transpose().getAsRowVectorView();}
    THerm updTranspose()    {return Base::transpose().updAsRowVectorView();}

    THerm operator~() const {return transpose();}
    THerm operator~()       {return updTranspose();}

    const VectorBase& operator+() const {return *this; }

    // Negation

    const TNeg& negate()    const {return *reinterpret_cast<const TNeg*>(this); }
    TNeg&       updNegate()       {return *reinterpret_cast<TNeg*>(this); }

    const TNeg& operator-() const {return negate();}
    TNeg&       operator-()       {return updNegate();}

    void resize(int m) { Base::resize(m,1); }
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
    typedef RowVectorBase< typename CNT<ELT>::TNeg >    TNeg;
    typedef VectorView_< typename CNT<ELT>::THerm >     THerm;
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

    RowVectorBase& operator*=(const StdNumber& t)     { Base::operator*=(t); return *this; }
    RowVectorBase& operator/=(const StdNumber& t)     { Base::operator/=(t); return *this; }
    RowVectorBase& operator+=(const RowVectorBase& r) { Base::operator+=(r); return *this; }
    RowVectorBase& operator-=(const RowVectorBase& r) { Base::operator-=(r); return *this; }  

    template <class EE> RowVectorBase& operator=(const RowVectorBase<EE>& b) 
      { Base::operator=(b);  return *this; } 
    template <class EE> RowVectorBase& operator+=(const RowVectorBase<EE>& b) 
      { Base::operator+=(b); return *this; } 
    template <class EE> RowVectorBase& operator-=(const RowVectorBase<EE>& b) 
      { Base::operator-=(b); return *this; } 

    // default destructor

    // fill current allocation with copies of element
    RowVectorBase& operator=(const ELT& t) { Base::operator=(t); return *this; }  

    // Implicit conversions are allowed to RowVector or Matrix, but not to Vector.   
    operator const RowVector_<ELT>&()     const { return *reinterpret_cast<const RowVector_<ELT>*>(this); }
    operator       RowVector_<ELT>&()           { return *reinterpret_cast<      RowVector_<ELT>*>(this); }
    operator const RowVectorView_<ELT>&() const { return *reinterpret_cast<const RowVectorView_<ELT>*>(this); }
    operator       RowVectorView_<ELT>&()       { return *reinterpret_cast<      RowVectorView_<ELT>*>(this); }
    
    operator const Matrix_<ELT>&()     const { return *reinterpret_cast<const Matrix_<ELT>*>(this); }
    operator       Matrix_<ELT>&()           { return *reinterpret_cast<      Matrix_<ELT>*>(this); } 
    operator const MatrixView_<ELT>&() const { return *reinterpret_cast<const MatrixView_<ELT>*>(this); }
    operator       MatrixView_<ELT>&()       { return *reinterpret_cast<      MatrixView_<ELT>*>(this); } 
    
    // Override MatrixBase indexing operators          
    const ELT& operator[](int j) const { return Base::operator()(0,j); }
    ELT&       operator[](int j)       { return Base::operator()(0,j); }
    const ELT& operator()(int j) const { return Base::operator()(0,j); }
    ELT&       operator()(int j)       { return Base::operator()(0,j); }
         
    // View creation      
    RowVectorView_<ELT> operator()(int j, int n) const { return Base::operator()(0,j,1,n).getAsRowVectorView(); }
    RowVectorView_<ELT> operator()(int j, int n)       { return Base::operator()(0,j,1,n).updAsRowVectorView(); }
 
    // Hermitian transpose.
    THerm transpose() const {return Base::transpose().getAsVectorView();}
    THerm updTranspose()    {return Base::transpose().updAsVectorView();}

    THerm operator~() const {return transpose();}
    THerm operator~()       {return updTranspose();}

    const RowVectorBase& operator+() const {return *this; }

    // Negation

    const TNeg& negate()    const {return *reinterpret_cast<const TNeg*>(this); }
    TNeg&       updNegate()       {return *reinterpret_cast<TNeg*>(this); }

    const TNeg& operator-() const {return negate();}
    TNeg&       operator-()       {return updNegate();}

    void resize(int n) { Base::resize(1,n); }
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
    explicit TmpMatrixView_(const TmpMatrixHelper<Scalar>& h) : Base(h) { }
    
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
    explicit TmpVectorViewT(const TmpMatrixHelper<S>& h) : Base(h) { }
    
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
    TmpRowVectorView_(const TmpMatrixHelper<S>& h) : Base(h) { }
    
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

    
} //namespace SimTK

#endif //_SimTK_SIMMATRIX_BIGMATRIX_H_
