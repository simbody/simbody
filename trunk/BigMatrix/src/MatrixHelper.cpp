/* Portions copyright (c) 2005-6 Stanford University and Michael Sherman.
 * Contributors:
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

#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/internal/ExceptionMacros.h"
#include "SimTKcommon/internal/Scalar.h"
#include "SimTKcommon/internal/SmallMatrix.h"
#include "SimTKcommon/internal/MatrixHelper.h"

#include "ElementFilter.h"
#include "DataDescriptor.h"

#include <iostream>
#include <cstdio>

namespace SimTK {

/**
 * This is the private implementation of the MatrixHelper<S> class. Note that this
 * must be explicitly instantiated in library-side source code for all possible
 * values of the template parameter S, which is just the set of all SimTK scalar
 * types.
 */
template <class S> 
class MatrixHelperRep {
    typedef typename CNT<S>::Number     Number;     // strips off negator from S
    typedef typename CNT<S>::StdNumber  StdNumber;  // turns conjugate to complex
    typedef typename CNT<S>::Precision  Precision;  // strips off complex from StdNumber
public:     
    // no default constructor, copy constructor suppressed

    // Local types for directing us to the right constructor at compile time.
    class ShallowCopy   { };
    class DeepCopy      { };
    class TransposeView { };
    
    // Matrix owner constructors
    
    // 0x0, fully resizable
    explicit MatrixHelperRep(int esz, int cppEsz) 
      : eltSize(esz),cppEltSize(cppEsz),view(0),data(0),myHandle(0) 
    {
        SimTK_ASSERT2(esz > 0 && cppEsz >= esz, 
            "MatrixHelper constructor was given a bad element size esz=%d, cppEsz=%d", esz, cppEsz); 
    }

    // Restore helper to its just-constructed state. We forget everything except
    // the element size, which can never change.
    void clear();

    // (m*esz) x n, resizable with optional restrictions 
    MatrixHelperRep(int esz, int cppEsz, int m, int n, bool lockNrow, bool lockNcol);
    
    // This is a non-resizeable owner (of the data descriptor) but shares pre-existing raw data.
    MatrixHelperRep(int esz, int cppEsz, int m, int n, int leadingDim, const S*); // read only
    MatrixHelperRep(int esz, int cppEsz, int m, int n, int leadingDim, S*);       // writable

    // Behavior of copy assignment depends on whether this is an owner or view. If
    // it's an owner it is resized and ends up a new, dense copy of rhs. If
    // it's a view, sizes must match and we copy elements from rhs to corresponding
    // elements of lhs.  
    MatrixHelperRep& copyAssign(const MatrixHelper<S>&);
    MatrixHelperRep& negatedCopyAssign(const MatrixHelper<typename CNT<S>::TNeg>&);

    // View assign always disconnects this helper from whatever view & data
    // it used to have (destructing as appropriate) and then makes it a new view
    // of the source. Writability is lost if the source is const, otherwise 
    // writability is inherited from the source.
    MatrixHelperRep& readOnlyViewAssign(const MatrixHelper<S>&);
    MatrixHelperRep& writableViewAssign(MatrixHelper<S>&);

    // These are like constructors for heap allocated MatrixHelperReps
    MatrixHelperRep* createReadOnlyShallowCopy() const;
    MatrixHelperRep* createWritableShallowCopy();

    MatrixHelperRep*                        createDeepCopy() const;
    MatrixHelperRep<typename CNT<S>::TNeg>* createNegatedDeepCopy() const;

    MatrixHelperRep* createReadOnlyBlockView(int i, int j, int m, int n) const;
    MatrixHelperRep* createWritableBlockView(int i, int j, int m, int n);

    MatrixHelperRep<typename CNT<S>::THerm>* createReadOnlyTransposedView() const;
    MatrixHelperRep<typename CNT<S>::THerm>* createWritableTransposedView();

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
    void addIn(const MatrixHelper<S>&);   
    void addIn(const MatrixHelper<typename CNT<S>::TNeg>&);   
    void subIn(const MatrixHelper<S>&); 
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

    void dump(const char* msg) const;

    // See comment in MatrixBase::matmul for an explanation.
    template <class SA, class SB>
    void matmul(const StdNumber& beta,   // applied to 'this'
                const StdNumber& alpha, const MatrixHelper<SA>& A, const MatrixHelper<SB>& B);
    
    // Here we use element indices but return the first scalar we find there.
    // This is useless except when the elements are scalars.
    const S& getScalar(int i, int j) const;
    S&       updScalar(int i, int j);      
    
        // Bookkeeping //
        
    int  nrow() const;
    int  ncol() const;
    long size() const; // nrow*ncol  

    void resize(int m, int n);
    void resizeKeep(int m, int n);

    void lockNRows();
    void lockNCols();
    void lockShape();
    void unlockNRows();
    void unlockNCols();
    void unlockShape();

    ~MatrixHelperRep();


    void setMyHandle(MatrixHelper<S>& h) {myHandle = &h;}
    const MatrixHelper<S>& getMyHandle() const {assert(myHandle); return *myHandle;}
    void clearMyHandle() {myHandle=0;}

private:
    // Matrix handle commitments.

    // We are always committed to a particular size element (size is minimum
    // # scalars required). We also make note here of the size that C++
    // allocates for one of these for use in communicating with ordinary
    // arrays when that is necessary.
    const int eltSize;      // # scalars per packed element
    const int cppEltSize;   // # scalars when C++ pads out the element

    // These are by default "Uncommitted", meaning we're happy to take on
    // matrices in any format. Otherwise the settings here limit what we'll
    // find acceptable in view & data.
    // Note: don't look at these to find out anything about the current
    // matrix. These are used only when the matrix is being created or
    // changed structurally. If there is a view, then that defines the
    // current matrix properties. Otherwise, the data descriptor does so.
    MatrixShape     shapeCommitment;
    MatrixSize      sizeCommitment;
    MatrixStructure structureCommitment;
    MatrixSparsity  sparsityCommitment;
    MatrixStorage   storageCommitment;
    MatrixCondition conditionCommitment;


    ElementFilter*     view;    // details are opaque
    DataDescriptor<S>* data;    // TODO: should be <StdNumber>?

    MatrixHelper<S>* myHandle; // point back to the owner handle
    friend class MatrixHelper<S>;

    // suppress assignment & copy
    MatrixHelperRep& operator=(const MatrixHelperRep&);
    MatrixHelperRep(const MatrixHelperRep&);

    MatrixHelperRep(int esz, ElementFilter* v, DataDescriptor<S>* d)
      : eltSize(esz), view(v), data(d) {
    }

    template <class T> friend class MatrixHelperRep;

    int getElementSize() const { return eltSize; }
    bool isDataOwned() const;
};

    ///////////////////
    // MATRIX HELPER //
    ///////////////////

// These are all pass-throughs to the MatrixHelperRep class.

template <class S> void
MatrixHelper<S>::clear() {
    rep->clear();
}

// Destructor. If we own the MatrixHelperRep it will point back here.
template <class S> 
MatrixHelper<S>::~MatrixHelper() {
    if (&rep->getMyHandle() == this)
        delete rep;
    rep=0;
}

template <class S> 
MatrixHelper<S>::MatrixHelper(int esz, int cppEsz) : rep(0) { 
    rep = new MatrixHelperRep<S>(esz,cppEsz);
    rep->setMyHandle(*this);
}
template <class S> 
MatrixHelper<S>::MatrixHelper(int esz, int cppEsz, int m, int n, bool lockNrow, bool lockNcol) : rep(0) { 
    rep = new MatrixHelperRep<S>(esz,cppEsz,m,n,lockNrow,lockNcol);
    rep->setMyHandle(*this);
}
// These are non-resizeable owners (of the data descriptor) but share pre-existing raw data.
// We have read-only and writable varieties.
template <class S> 
MatrixHelper<S>::MatrixHelper(int esz, int cppEsz, int m, int n, int leadingDim, const S* s) : rep(0) {
    rep = new MatrixHelperRep<S>(esz,cppEsz,m,n,leadingDim,s);
    rep->setMyHandle(*this);
}
template <class S> 
MatrixHelper<S>::MatrixHelper(int esz, int cppEsz, int m, int n, int leadingDim, S* s) : rep(0) {
    rep = new MatrixHelperRep<S>(esz,cppEsz,m,n,leadingDim,s);
    rep->setMyHandle(*this);
}

// copy constructor is suppressed; these are closest things but require
// specification of deep or shallow copy

// Create a read-only view of existing data.
template <class S>
MatrixHelper<S>::MatrixHelper(const MatrixHelper& h, const ShallowCopy&) : rep(0) {
    rep = h.rep->createReadOnlyShallowCopy();
    rep->setMyHandle(*this);
}
// Create a (possibly) writable view of existing data. 
template <class S>
MatrixHelper<S>::MatrixHelper(MatrixHelper& h, const ShallowCopy&) : rep(0) {
    rep = h.rep->createWritableShallowCopy();
    rep->setMyHandle(*this);
}
// Deep copy. "This" will have no view, and the result is always writable, packed storage.
template <class S>
MatrixHelper<S>::MatrixHelper(const MatrixHelper& h, const DeepCopy&) : rep(0) {
    rep = h.rep->createDeepCopy();
    rep->setMyHandle(*this);
}
// Negated deep copy. We'll get a brand new, viewless, writable, packed copy with 
// the same values as the original (duh, that's what "copy" means) -- BUT, the
// physical floating point representations will all have been negated since we're
// copying a Matrix whose elements are of type negator<S> while ours are type S.
template <class S>
MatrixHelper<S>::MatrixHelper(const MatrixHelper<typename CNT<S>::TNeg>& h, const DeepCopy&) : rep(0) {
    rep = h.rep->createNegatedDeepCopy();
    rep->setMyHandle(*this);
}
// construct read-only view
template <class S> 
MatrixHelper<S>::MatrixHelper(const MatrixHelper& h, int i, int j, int m, int n) : rep(0) {
    rep = h.rep->createReadOnlyBlockView(i,j,m,n);
    rep->setMyHandle(*this);
}
// construct writable view
template <class S> 
MatrixHelper<S>::MatrixHelper(MatrixHelper& h, int i, int j, int m, int n) : rep(0) {
    rep = h.rep->createWritableBlockView(i,j,m,n);
    rep->setMyHandle(*this);
}
// Construct read only transposed view of passed-in helper.
template <class S>
MatrixHelper<S>::MatrixHelper(const MatrixHelper<typename CNT<S>::THerm>& h,
                              const TransposeView&) : rep(0) {
    rep = h.rep->createReadOnlyTransposedView();
    rep->setMyHandle(*this);
}

// Construct (possibly) writable transposed view of passed-in helper.
template <class S>
MatrixHelper<S>::MatrixHelper(MatrixHelper<typename CNT<S>::THerm>& h,
                              const TransposeView&) : rep(0) {
    rep = h.rep->createWritableTransposedView();
    rep->setMyHandle(*this);
}

template <class S> MatrixHelper<S>&
MatrixHelper<S>::copyAssign(const MatrixHelper& h) {
    rep->copyAssign(h);
    return *this;
}
// Like copy assign but the source has elements with the opposite interpretation
// of sign. Note that the result still has the original value, but the in-memory
// representation will be different.
template <class S> MatrixHelper<S>&
MatrixHelper<S>::negatedCopyAssign(const MatrixHelper<typename CNT<S>::TNeg>& h) {
    rep->negatedCopyAssign(h);
    return *this;
}
// Create a read-only view of existing data. Element size of source and destination
// must match.
template <class S> MatrixHelper<S>&
MatrixHelper<S>::readOnlyViewAssign(const MatrixHelper& h) {
    rep->readOnlyViewAssign(h);
    return *this;
}
// Create a (possibly) writable view of existing data. 
template <class S> MatrixHelper<S>&
MatrixHelper<S>::writableViewAssign(MatrixHelper& h) {
    rep->writableViewAssign(h);
    return *this;
}
template <class S> void
MatrixHelper<S>::scaleBy(const typename CNT<S>::StdNumber& s) {
    rep->scaleBy(s);
}
template <class S> void
MatrixHelper<S>::addIn(const MatrixHelper& h) {
    rep->addIn(h);
}
template <class S> void
MatrixHelper<S>::addIn(const MatrixHelper<typename CNT<S>::TNeg>& h) {
    subIn(reinterpret_cast<const MatrixHelper<S>&>(h));
}   
template <class S> void
MatrixHelper<S>::subIn(const MatrixHelper& h) {
    rep->subIn(h);
}
template <class S> void
MatrixHelper<S>::subIn(const MatrixHelper<typename CNT<S>::TNeg>& h) {
    addIn(reinterpret_cast<const MatrixHelper<S>&>(h));
}
template <class S> void
MatrixHelper<S>::fillWith(const S* eltp) {
    rep->fillWith(eltp);
}
template <class S> void 
MatrixHelper<S>::copyInByRows(const S* elts) {
    rep->copyInByRows(elts);
}

template <class S> void
MatrixHelper<S>::invertInPlace() {
    rep->invertInPlace();
}

template <class S> void
MatrixHelper<S>::dump(const char* msg) const {
    rep->dump(msg);
}

template <class S> const S* 
MatrixHelper<S>::getElt(int i, int j) const {
    return rep->getElt(i,j);
}
template <class S> S* 
MatrixHelper<S>::updElt(int i, int j) {
    return rep->updElt(i,j);
}
template <class S> const S& 
MatrixHelper<S>::getScalar(int i, int j) const {
    return rep->getScalar(i,j);
} 

template <class S> S& 
MatrixHelper<S>::updScalar(const int i, const int j) {
    return rep->updScalar(i,j);
}
template <class S> int 
MatrixHelper<S>::nrow() const {return rep->nrow();}
template <class S> int 
MatrixHelper<S>::ncol() const {return rep->ncol();} 
template <class S> long 
MatrixHelper<S>::size() const {return rep->size();}

template <class S> void 
MatrixHelper<S>::resize(int m, int n) {
    rep->resize(m,n);         
} 

template <class S> void 
MatrixHelper<S>::resizeKeep(int m, int n) {
    rep->resizeKeep(m,n);     
} 
 
// "Lock" operations are only allowed on owners, since the caller is
// probably confused if they are calling this on a view. However,
// it is fine to lock a dimension of an unresizable owner, just
// don't try to unlock it!
template <class S> void MatrixHelper<S>::lockNRows() {rep->lockNRows();}
template <class S> void MatrixHelper<S>::lockNCols() {rep->lockNCols();}
template <class S> void MatrixHelper<S>::lockShape() {rep->lockShape();}

//
// These are very fussy. You can only call them on data descriptor owners, and
// then only on those which were originally unlocked in the indicated dimension.
// For example, you can lock and unlock the number of rows in a Vector, but
// you can never "unlock" the number of columns which was frozen at 1 on construction.
//
template <class S> void MatrixHelper<S>::unlockNRows() {rep->unlockNRows();}
template <class S> void MatrixHelper<S>::unlockNCols() {rep->unlockNCols();}
template <class S> void MatrixHelper<S>::unlockShape() {rep->unlockShape();}

template <class S> void
MatrixHelper<S>::sum(S* const answer) const {
    rep->sum(answer);
}     
template <class S> void
MatrixHelper<S>::colSum(int j, S* const answer) const {
    rep->colSum(j,answer);
}
template <class S> void
MatrixHelper<S>::rowSum(int i, S* const answer) const {
    rep->rowSum(i,answer);
}
template <class S> void
MatrixHelper<S>::fillWithScalar(const typename CNT<S>::StdNumber& s) {
    rep->fillWithScalar(s);
}
template <class S> S
MatrixHelper<S>::scalarColSum(int j) const {
    return rep->scalarColSum(j);
}
template <class S> S
MatrixHelper<S>::scalarRowSum(int i) const {
    return rep->scalarRowSum(i);
}
template <class S> S
MatrixHelper<S>::scalarSum() const {
    return rep->scalarSum();
}     

template <class S> bool 
MatrixHelper<S>::hasContiguousData() const {
    return rep->view==0 && rep->data && rep->data->hasContiguousData();
}
template <class S> long 
MatrixHelper<S>::getContiguousDataLength() const {
    assert(hasContiguousData());
    return rep->data->getContiguousDataLength();
}
template <class S> const S* 
MatrixHelper<S>::getContiguousData() const {
    assert(hasContiguousData());
    return rep->data->getContiguousData();
}
template <class S> S* 
MatrixHelper<S>::updContiguousData() {
    assert(hasContiguousData());
    return rep->data->updContiguousData();
}
template <class S> void 
MatrixHelper<S>::replaceContiguousData(S* newData, long length, bool takeOwnership) {
    assert(hasContiguousData());
    return rep->data->replaceContiguousData(newData,length,takeOwnership);
}
template <class S> void 
MatrixHelper<S>::replaceContiguousData(const S* newData, long length) {
    assert(hasContiguousData());
    return rep->data->replaceContiguousData(newData,length);
}
template <class S> void 
MatrixHelper<S>::swapOwnedContiguousData(S* newData, int length, S*& oldData) {
    assert(hasContiguousData());
    return rep->data->swapOwnedContiguousData(newData,length,oldData);
}


    ///////////////////////
    // MATRIX HELPER REP //
    ///////////////////////

// These are the non-inline implementations for all the MatrixHelper pass-throughs,
// plus other routines which have no client-side equivalents.

// Element size and handle remain unchanged by clear().
template <class S> void
MatrixHelperRep<S>::clear() {
    if (view==0) delete data;   // i.e., owner
    data=0;
    delete view; view=0;
}

template <class S> 
MatrixHelperRep<S>::MatrixHelperRep(int esz, int cppEsz, int m, int n, bool lockNrow, bool lockNcol)
  : eltSize(esz), cppEltSize(cppEsz), view(0), data(0), myHandle(0)
{ 
    SimTK_ASSERT2(esz > 0 && cppEsz >= esz, 
        "MatrixHelper constructor was given a bad element size esz=%d, cppEsz=%d", esz, cppEsz); 

    data = new DataDescriptor<S>(esz,m,n,lockNrow,lockNcol);
}

// These are non-resizeable owners (of the data descriptor) but share pre-existing raw data.
// We have read-only and writable varieties.
template <class S> 
MatrixHelperRep<S>::MatrixHelperRep(int esz, int cppEsz, int m, int n, int leadingDim, const S* s)
  : eltSize(esz), cppEltSize(cppEsz), view(0), data(0), myHandle(0)
{ 
    SimTK_ASSERT2(esz > 0 && cppEsz >= esz, 
        "MatrixHelper constructor was given a bad element size esz=%d, cppEsz=%d", esz, cppEsz); 

    data = new DataDescriptor<S>(esz,m,n,leadingDim,s);
}
template <class S> 
MatrixHelperRep<S>::MatrixHelperRep(int esz, int cppEsz, int m, int n, int leadingDim, S* s)
  : eltSize(esz), cppEltSize(cppEsz), view(0), data(0), myHandle(0)
{ 
    SimTK_ASSERT2(esz > 0 && cppEsz >= esz, 
        "MatrixHelper constructor was given a bad element size esz=%d, cppEsz=%d", esz, cppEsz);

    data = new DataDescriptor<S>(esz,m,n,leadingDim,s);
}

// copy constructor is suppressed; these are closest things but require
// specification of deep or shallow copy

// Create a read-only view of existing data.
template <class S> MatrixHelperRep<S>*
MatrixHelperRep<S>::createReadOnlyShallowCopy() const
{
    assert(data);

    MatrixHelperRep* copy = new MatrixHelperRep(eltSize,cppEltSize);

    copy->view = view 
        ? new ElementFilter(*view) // loses writability
        : new ElementFilter(false,data->nrowElt(eltSize),data->ncolElt(eltSize),
                             ElementFilter::Indexer(0,0));
    copy->data = data;
    return copy;
}

// Create a (possibly) writable view of existing data. 
template <class S> MatrixHelperRep<S>*
MatrixHelperRep<S>::createWritableShallowCopy() 
{
    assert(data);

    MatrixHelperRep* copy = new MatrixHelperRep(eltSize,cppEltSize);

    copy->view = view 
        ? new ElementFilter(*view, true)    // keep writability if possible
        : new ElementFilter(true,data->nrowElt(eltSize),data->ncolElt(eltSize),
                             ElementFilter::Indexer(0,0));
    copy->data = data;
    return copy;
}

// Deep copy. Copy will have no view, and the result is always writable, packed storage.
template <class S> MatrixHelperRep<S>*
MatrixHelperRep<S>::createDeepCopy() const
{
    MatrixHelperRep* copy = new MatrixHelperRep(eltSize,cppEltSize);

    if (data)
        copy->data = view 
            ? new DataDescriptor<S>(eltSize,*view,*data)
            : new DataDescriptor<S>(*data);

    return copy;
}

template <class S> MatrixHelperRep<typename CNT<S>::TNeg>*
MatrixHelperRep<S>::createNegatedDeepCopy() const
{
    MatrixHelperRep<typename CNT<S>::TNeg>* copy = 
        new MatrixHelperRep<typename CNT<S>::TNeg>(eltSize,cppEltSize);

    if (data) {
        const DataDescriptor<typename CNT<S>::TNeg>* negSrc = 
            reinterpret_cast<const DataDescriptor<typename CNT<S>::TNeg>*>(data);
        copy->data = view 
            ? new DataDescriptor<typename CNT<S>::TNeg>(eltSize,*view,*negSrc)
            : new DataDescriptor<typename CNT<S>::TNeg>(*negSrc);
        copy->scaleBy(typename CNT<S>::StdNumber(-1));
    }

    return copy;
}

template <class S> MatrixHelperRep<S>*
MatrixHelperRep<S>::createReadOnlyBlockView(int i, int j, int m, int n) const
{
    MatrixHelperRep* copy = new MatrixHelperRep(eltSize,cppEltSize);

    const ElementFilter::Indexer ix(i,j);

    copy->view = view 
        ? new ElementFilter(*view,false,m,n,ix)
        : new ElementFilter(false,m,n,ix);

    copy->data = data;
    return copy;
}


template <class S> MatrixHelperRep<S>*
MatrixHelperRep<S>::createWritableBlockView(int i, int j, int m, int n)
{
    MatrixHelperRep* copy = new MatrixHelperRep(eltSize,cppEltSize);

    const ElementFilter::Indexer ix(i,j);

    copy->view = view 
        ? new ElementFilter(*view,true,m,n,ix)
        : new ElementFilter(true,m,n,ix);

    copy->data = data;
    return copy;
}

template <class S> MatrixHelperRep<typename CNT<S>::THerm>*
MatrixHelperRep<S>::createReadOnlyTransposedView() const
{
    MatrixHelperRep<typename CNT<S>::THerm>* copy = 
        new MatrixHelperRep<typename CNT<S>::THerm>(eltSize,cppEltSize);

    const ElementFilter::Indexer ix(ElementFilter::Indexer(0,0).transpose());

    copy->view = view
        ? new ElementFilter(*view,false,ncol(),nrow(),ix)
        : new ElementFilter(false,ncol(),nrow(),ix);

    copy->data = reinterpret_cast<DataDescriptor<typename CNT<S>::THerm>*>(data);
    return copy;
}

template <class S> MatrixHelperRep<typename CNT<S>::THerm>*
MatrixHelperRep<S>::createWritableTransposedView()
{
    MatrixHelperRep<typename CNT<S>::THerm>* copy = 
        new MatrixHelperRep<typename CNT<S>::THerm>(eltSize,cppEltSize);

    const ElementFilter::Indexer ix(ElementFilter::Indexer(0,0).transpose());

    copy->view = view
        ? new ElementFilter(*view,true,ncol(),nrow(),ix)
        : new ElementFilter(true,ncol(),nrow(),ix);

    copy->data = reinterpret_cast<DataDescriptor<typename CNT<S>::THerm>*>(data);
    return copy;
}

// assignment depends on whether lhs (this) is a view
template <class S> MatrixHelperRep<S>&
MatrixHelperRep<S>::copyAssign(const MatrixHelper<S>& h) {
    const MatrixHelperRep& hrep = *h.rep;
    if (&hrep == this) return *this;

    assert(getElementSize()==hrep.getElementSize());

    // OK if we can reassign, or if the sizes match, or if
    // this is an empty Vector (0x1) or RowVector (1x0) and
    // the source is just an empty handle (0x0).
    assert(view==0 
           || (nrow()==hrep.nrow() && ncol()==hrep.ncol())
           || (nrow()*ncol()==0 && hrep.nrow()==0 && hrep.ncol()==0));

    if (view == 0) { // 'this' is the owner
        resize(hrep.nrow(), hrep.ncol());
        if (hrep.data) {
            if (hrep.view) data->copyFromCompatibleSourceView(getElementSize(),*hrep.view,*hrep.data);
            else data->copyFromCompatibleSourceData(*hrep.data);
        }
    } else {
        // "this" has a view. TODO: this is really bad!
        for (int j=0; j<ncol(); ++j)
            for (int i=0; i<nrow(); ++i) 
                DataDescriptor<S>::copyElement(getElementSize(),
                                               updElt(i,j), hrep.getElt(i,j));
    }

    return *this;
}

template <class S> MatrixHelperRep<S>&
MatrixHelperRep<S>::negatedCopyAssign(const MatrixHelper<typename CNT<S>::TNeg>& src) {
    // TODO: avoid making two passes through the data here
    copyAssign(reinterpret_cast<const MatrixHelper<S>&>(src));
    scaleBy(typename CNT<S>::StdNumber(-1));
    return *this;
}

// Create a read-only view of existing data. Element size of source and destination
// must match.
template <class S> MatrixHelperRep<S>&
MatrixHelperRep<S>::readOnlyViewAssign(const MatrixHelper<S>& h)
{
    const MatrixHelperRep& hrep = *h.rep;

    assert((eltSize == hrep.eltSize) && hrep.data);
    clear();

    view = hrep.view 
        ? new ElementFilter(*hrep.view) // loses writability
        : new ElementFilter(false,hrep.data->nrowElt(eltSize),hrep.data->ncolElt(eltSize),
                             ElementFilter::Indexer(0,0));
    data = hrep.data;
    return *this;
}

// Create a (possibly) writable view of existing data. 
template <class S> MatrixHelperRep<S>&
MatrixHelperRep<S>::writableViewAssign(MatrixHelper<S>& h) 
{
    MatrixHelperRep& hrep = *h.rep;

    assert((eltSize == hrep.eltSize) && hrep.data);
    clear();

    view = hrep.view 
        ? new ElementFilter(*hrep.view, true)    // keep writability if possible
        : new ElementFilter(true,hrep.data->nrowElt(eltSize),hrep.data->ncolElt(eltSize),
                             ElementFilter::Indexer(0,0));
    data = hrep.data;
    return *this;
}

template <class S> void
MatrixHelperRep<S>::scaleBy(const typename CNT<S>::StdNumber& s) {
    if (!view) { data->scaleBy(s); return; }
    // XXX -- really, really bad! Optimize for contiguous data!
    const int sz = getElementSize();
    for (int j=0; j<view->ncol(); ++j)
        for (int i=0; i<view->nrow(); ++i) 
            DataDescriptor<S>::scaleElement(sz,updElt(i,j),s);
}  
     
template <class S> void
MatrixHelperRep<S>::addIn(const MatrixHelper<S>& h) {
    const MatrixHelperRep& hrep = *h.rep;

    assert(nrow()==hrep.nrow() && ncol()==hrep.ncol());
    assert(getElementSize()==hrep.getElementSize());
    // XXX -- really, really bad! Optimize for contiguous data, missing views, etc.!
    const int sz = getElementSize();
    for (int j=0; j<ncol(); ++j)
        for (int i=0; i<nrow(); ++i)
            DataDescriptor<S>::addToElement(sz,updElt(i,j),hrep.getElt(i,j));
} 
     
template <class S> void
MatrixHelperRep<S>::subIn(const MatrixHelper<S>& h) {
    const MatrixHelperRep& hrep = *h.rep;

    assert(nrow()==hrep.nrow() && ncol()==hrep.ncol());
    assert(getElementSize()==hrep.getElementSize());
    // XXX -- really, really bad! Optimize for contiguous data, missing views, etc.!
    const int sz = getElementSize();
    for (int j=0; j<ncol(); ++j)
        for (int i=0; i<nrow(); ++i)
            DataDescriptor<S>::subFromElement(sz,updElt(i,j),hrep.getElt(i,j));
}  

template <class S> void
MatrixHelperRep<S>::fillWith(const S* eltp) {
    // XXX -- really, really bad! Optimize for contiguous data, missing views, etc.!
    const int sz = getElementSize();
    for (int j=0; j<ncol(); ++j)
        for (int i=0; i<nrow(); ++i)
            DataDescriptor<S>::copyElement(sz,updElt(i,j),eltp);
} 

template <class S> void 
MatrixHelperRep<S>::copyInByRows(const S* elts) {
    // XXX -- really, really bad! Optimize for contiguous data, missing views, etc.!
    const int sz = getElementSize();
    for (int j=0; j<ncol(); ++j)
        for (int i=0; i<nrow(); ++i)
            DataDescriptor<S>::copyElement(sz,updElt(i,j),elts+i*ncol()+j);
}

template <class S> 
MatrixHelperRep<S>::~MatrixHelperRep()
{
    if (view==0) delete data;
    delete view;
}

template <class S> void
MatrixHelperRep<S>::invertInPlace() {
    assert(view==0 && eltSize == 1);
    if (data)
        data->invertInPlace();
}

template <class S> static void
dumpElt(const S* p, int sz) {
    if (sz > 1) std::cout << "{";
    for (int k=0; k<sz; ++k) {
        if (k>0) std::cout << " ";
        std::cout << *p++;
    }
    if (sz > 1) std::cout << "}";
}

template <class S> void
MatrixHelperRep<S>::dump(const char* msg) const {
    if (msg) 
        std::cout << std::string(msg) << std::endl;
    std::cout << "Matrix " << nrow() << " X " << ncol() << " "
              << getElementSize() << "-scalar entries:" << std::endl;
    if (nrow()*ncol() == 0) {
        std::cout << "<EMPTY>" << std::endl;
        return;
    }

    const std::streamsize oldsz = std::cout.precision(20); 
    for (int i=0; i<nrow(); ++i) {
        for (int j=0; j<ncol(); ++j) {
            if (j>0) std::cout << "\t";
            dumpElt(getElt(i,j), getElementSize());
        }
        std::cout << std::endl;
    }
    std::cout.precision(oldsz);

}

template <class S> const S* 
MatrixHelperRep<S>::getElt(int i, int j) const {
    assert(0<=i && i<nrow() && 0<=j && j<ncol());
    if (view) return data->getElt(getElementSize(),view->r(i,j), view->c(i,j));
    else return data->getElt(getElementSize(),i,j);
}       
    
template <class S> S* 
MatrixHelperRep<S>::updElt(int i, int j) {
    assert(0<=i && i<nrow() && 0<=j && j<ncol());
    if (view) {
        assert(view->isViewWritable());
        return data->updElt(getElementSize(), view->r(i,j), view->c(i,j));
    }
    return data->updElt(getElementSize(),i,j);
}

template <class S> const S& 
MatrixHelperRep<S>::getScalar(const int i, const int j) const {
    assert(0<=i && i<nrow() && 0<=j && j<ncol());
    if (view) return data->getScalar(view->r(i,j), view->c(i,j));
    else return data->getScalar(i,j);
} 

template <class S> S& 
MatrixHelperRep<S>::updScalar(const int i, const int j) {
    assert(0<=i && i<nrow() && 0<=j && j<ncol());
    if (view) {
        assert(view->isViewWritable());
        return data->updScalar(view->r(i,j), view->c(i,j));
    }
    return data->updScalar(i,j);
}

template <class S> int 
MatrixHelperRep<S>::nrow() const {return view?view->nrow():(data?data->nrowElt(eltSize):0);}

template <class S> int 
MatrixHelperRep<S>::ncol() const {return view?view->ncol():(data?data->ncolElt(eltSize):0);} 

template <class S> long 
MatrixHelperRep<S>::size() const {return view?view->size():(data?data->sizeElt(eltSize):0);}

template <class S> bool
MatrixHelperRep<S>::isDataOwned() const {return view==0;} 

template <class S> void 
MatrixHelperRep<S>::resize(int m, int n) {
    if (nrow()==m && ncol()==n) return;
    assert(isDataOwned());
    if (data) data->resize(getElementSize(),m,n);
    else data=new DataDescriptor<S>(getElementSize(),m,n,false,false);          
} 

template <class S> void 
MatrixHelperRep<S>::resizeKeep(int m, int n) {
    if (nrow()==m && ncol()==n) return;
    assert(isDataOwned());
    if (data) data->resizeKeep(getElementSize(),m,n);
    else data=new DataDescriptor<S>(getElementSize(),m,n,false,false);          
} 
 
// "Lock" operations are only allowed on owners, since the caller is
// probably confused if they are calling this on a view. However,
// it is fine to lock a dimension of an unresizable owner, just
// don't try to unlock it!
template <class S> void MatrixHelperRep<S>::lockNRows() {
    assert(data && isDataOwned());
    data->lockNRows();
}
template <class S> void MatrixHelperRep<S>::lockNCols() {
    assert(data && isDataOwned());
    data->lockNCols();
}
template <class S> void MatrixHelperRep<S>::lockShape() {
    assert(data && isDataOwned());
    data->lockNRows(); 
    data->lockNCols();
}

//
// These are very fussy. You can only call them on data descriptor owners, and
// then only on those which were originally unlocked in the indicated dimension.
// For example, you can lock and unlock the number of rows in a Vector, but
// you can never "unlock" the number of columns which was frozen at 1 on construction.
//
template <class S> void MatrixHelperRep<S>::unlockNRows() {
    assert(data && isDataOwned());
    data->unlockNRows();
}
template <class S> void MatrixHelperRep<S>::unlockNCols() {
    assert(data && isDataOwned());
    data->unlockNCols();
}
template <class S> void MatrixHelperRep<S>::unlockShape() {
    assert(data && isDataOwned());
    data->unlockNRows(); 
    data->unlockNCols();
}

template <class S> void
MatrixHelperRep<S>::sum(S* const answer) const {
    const int sz = getElementSize();
    if (sz==1) *answer = scalarSum();
    else {
        S* csum = new S[sz];
        DataDescriptor<S>::fillElement(sz, answer, typename CNT<S>::StdNumber(0));
        for (int j=0; j<ncol(); ++j) {
            colSum(j, csum);
            DataDescriptor<S>::addToElement(sz, answer, csum); // answer+=csum
        }
        delete csum;
    }
}        

template <class S> void
MatrixHelperRep<S>::colSum(const int j, S* const answer) const {
    const int sz = getElementSize();
    if (sz==1) *answer = scalarColSum(j);
    else {
        DataDescriptor<S>::fillElement(sz, answer, typename CNT<S>::StdNumber(0));
        for (int i=0; i<nrow(); ++i)
            DataDescriptor<S>::addToElement(sz, answer, getElt(i,j));
    }
} 

template <class S> void
MatrixHelperRep<S>::rowSum(const int i, S* const answer) const {
    const int sz = getElementSize();
    if (sz==1) *answer = scalarRowSum(i);
    else {
        DataDescriptor<S>::fillElement(sz, answer, typename CNT<S>::StdNumber(0));
        for (int j=0; j<ncol(); ++j)
            DataDescriptor<S>::addToElement(sz, answer, getElt(i,j));
    }
} 

template <class S> void
MatrixHelperRep<S>::fillWithScalar(const typename CNT<S>::StdNumber& s) {
    // XXX -- really, really bad! Optimize for contiguous data, missing views, etc.!
    const int sz = getElementSize();
    for (int j=0; j<ncol(); ++j)
        for (int i=0; i<nrow(); ++i)
            DataDescriptor<S>::fillElement(sz,updElt(i,j),s);
} 

template <class S> S
MatrixHelperRep<S>::scalarColSum(const int j) const {
    assert(getElementSize()==1);
    S sum = S(0);
    for (int i=0; i<nrow(); ++i) sum += getScalar(i,j);
    return sum;
} 

template <class S> S
MatrixHelperRep<S>::scalarRowSum(const int i) const {
    assert(getElementSize()==1);
    S sum = S(0);
    for (int j=0; j<ncol(); ++j) sum += getScalar(i,j);
    return sum;
} 

template <class S> S
MatrixHelperRep<S>::scalarSum() const {
    assert(getElementSize()==1);
    S sum = S(0);
    for (int j=0; j<ncol(); ++j) sum += scalarColSum(j);
    return sum;
}     

// Instantiations for each of the 18 Scalar types. (These will instantiate
// the MatrixHelperRep<S> classes also.

template class MatrixHelper< float >;
template class MatrixHelper< double >;
template class MatrixHelper< long double >;
template class MatrixHelper< std::complex<float> >;
template class MatrixHelper< std::complex<double> >;
template class MatrixHelper< std::complex<long double> >;
template class MatrixHelper< conjugate<float> >;
template class MatrixHelper< conjugate<double> >;
template class MatrixHelper< conjugate<long double> >;

template class MatrixHelper<negator< float > >;
template class MatrixHelper<negator< double > >;
template class MatrixHelper<negator< long double > >;
template class MatrixHelper<negator< std::complex<float> > >;
template class MatrixHelper<negator< std::complex<double> > >;
template class MatrixHelper<negator< std::complex<long double> > >;
template class MatrixHelper<negator< conjugate<float> > >;
template class MatrixHelper<negator< conjugate<double> > >;
template class MatrixHelper<negator< conjugate<long double> > >;
   
} // namespace SimTK   
