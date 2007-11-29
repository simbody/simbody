#ifndef SimTK_SIMMATRIX_DATA_DESCRIPTOR_H_
#define SimTK_SIMMATRIX_DATA_DESCRIPTOR_H_

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

/**@file
 * Here we declare the DataDescriptor base class. TODO: this will be abstract
 * to support different data formats. This class is part of the private
 * implementation of the MatrixHandle class and is not visible to client
 * side code.
 */

#include "SimTKcommon/Scalar.h"
#include "SimTKcommon/SmallMatrix.h"

#include "SimTKcommon/internal/BigMatrix.h"

#include "ElementFilter.h"


namespace SimTK {


/**
 * This class describes the actual memory used to hold data values for Matrix and
 * Vector objects. The data member here points to the physical memory, which may
 * be heap that is owned by this object or may be memory owned by someone else,
 * such as a Fortran or C array that we're trying to access in a civilized way.
 * 
 * If we're the owner, then the data may be resizable. Normally it is fully
 * resizable, but when the higher level client is a Vector or RowVector it
 * will be resizable only in number of rows or number of columns, respectively.
 *
 * Note that although a Matrix can contain any Composite Numerical Type (for
 * example Mat<5,5>) this class is templatized only by *scalar* type, that is, 
 * float, double, long double, the corresponding complex and conjugate types,
 * and negators of these types. It is thus instantiated explicitly exactly
 * 18 times in the source file corresponding to this header, and that's
 * it! We are ignorant here of the size and structure of the individual elements;
 * that is considered a matter for the client-side Matrix class.
 * 
 * Here we deal with memory consisting of contiguous scalar values, currently
 * limited to column-ordered storage, possibly with space between the columns
 * as is typical with LAPACK-style storage (that is, the "leading dimension"
 * may be larger than the number of rows). We are happy to perform "element"
 * services here if someone passes in the element size; in that case we 
 * interpret our m scalar rows to consist of m/elementSize elements instead.
 * Very, very bad things will happen to anyone who tries to use an element
 * size which does not divide evenly into m. In kind moods, we may verify
 * that here via assert() but don't count on it. Note that the number of
 * columns is not affected by the alleged element size.
 */
template <class S> class DataDescriptor {
    typedef typename CNT<S>::TReal      TReal;
    typedef typename CNT<S>::TNeg       TNeg;
    typedef typename CNT<S>::THerm      THerm;
    typedef typename CNT<S>::Number     Number;     // i.e., without negator<>
    typedef typename CNT<S>::StdNumber  StdNumber;  // i.e., without conjugate<>
public:
    /// Construct an empty (0x0), writable, fully resizable owner.
    DataDescriptor() 
        : rowsLockedOnConstruction(false), colsLockedOnConstruction(false), m(0), n(0), 
          owner(true), writable(true), leadingDim(0), data(0) 
    { 
        rowsCurrentlyLocked = rowsLockedOnConstruction;
        colsCurrentlyLocked = colsLockedOnConstruction;
    }
 
    /// Construct and pre-allocate a resizable (mXsz)Xn owner with optional locking of
    /// the number of rows or columns (or both) at their initial values. No
    /// heap allocation occurs if one of the dimensions is 0, as is common
    /// for initial creation of a row or column.       
    DataDescriptor(int sz, int mm, int nn, bool mlock, bool nlock)
        : rowsLockedOnConstruction(mlock), colsLockedOnConstruction(nlock), m(mm*sz), n(nn), 
          owner(true), writable(true), leadingDim(mm*sz),  
          data(mm*nn*sz ? new S[mm*nn*sz] : 0)
    {
        assert(sz > 0 && mm >= 0 && nn >= 0);
        rowsCurrentlyLocked = rowsLockedOnConstruction;
        colsCurrentlyLocked = colsLockedOnConstruction;
    #ifndef NDEBUG
        fillWithScalar(CNT<StdNumber>::getNaN());
    #endif

    }

    /// Create a read-only descriptor for someone else's data. We're assuming it is 
    /// organized as nn columns of mm elements of size sz scalars, with leading
    /// dimension ldim *elements*. This is NOT resizable.
    DataDescriptor(int sz, int mm, int nn, int ldim, const S* s)
        : rowsLockedOnConstruction(true), colsLockedOnConstruction(true), m(mm*sz), n(nn), 
          owner(false), writable(false), leadingDim(ldim*sz),  
          data(const_cast<S*>(s))
    { 
        assert(sz > 0 && mm >= 0 && nn >= 0 && ldim >= 0);
        rowsCurrentlyLocked = rowsLockedOnConstruction;
        colsCurrentlyLocked = colsLockedOnConstruction;
    }
          
    /// Shared data constructor allowing writing. 
    DataDescriptor(int sz, int mm, int nn, int ldim, S* s)
        : rowsLockedOnConstruction(true), colsLockedOnConstruction(true), m(mm*sz), n(nn), 
          owner(false), writable(true), leadingDim(ldim*sz),  
          data(s) 
    { 
        assert(sz > 0 && mm >= 0 && nn >= 0 && ldim >= 0);
        rowsCurrentlyLocked = rowsLockedOnConstruction;
        colsCurrentlyLocked = colsLockedOnConstruction;
    }
         
    /// Caution: this is a *viewless* deep copy of the whole of the source. We'll 
    /// inherit rowsLocked, etc. but of course we are now the owner regardless of
    /// whether the source was, and leadingDim becomes just the number of rows.
    /// We also have write access since this is new data.
    DataDescriptor(const DataDescriptor& s)
        : rowsLockedOnConstruction(s.rowsCurrentlyLocked), 
          colsLockedOnConstruction(s.colsCurrentlyLocked), m(s.m), n(s.n), 
          owner(true), writable(true), leadingDim(s.m),
          data(s.m*s.n ? new S[s.m*s.n] : 0) 
    {
        rowsCurrentlyLocked = rowsLockedOnConstruction;
        colsCurrentlyLocked = colsLockedOnConstruction;

        copyFromCompatibleSourceData(s);
    }

    void copyFromCompatibleSourceData(const DataDescriptor& src) {
        assert(writable);
        assert(m == src.m && n == src.n);
        copyInDataByColumn(src.data, src.leadingDim);
    }

    void copyFromCompatibleSourceView(int esz, const ElementFilter& view, const DataDescriptor& src) {
        assert(esz > 0); assert(writable); 
        assert(nrowElt(esz) == view.nrow() && ncolElt(esz) == view.ncol());

        // bad! XXX
        for (int j=0; j < view.ncol(); ++j)
            for (int i=0; i < view.nrow(); ++i)
                copyInOneElt(esz,i,j,src.getElt(esz,view.r(i,j), view.c(i,j))); 
    }


    /// Construction of new "viewless" data from a view of existing data. We 
    /// will NOT inherit rowsLocked & colsLocked since the view might be
    /// transposing -- reset the locks if you care. Of course we are the
    /// owner of the new (writable) data, which will be stored compactly regardless
    /// of how the original was stored.
    DataDescriptor(int esz, const ElementFilter& v, const DataDescriptor& s)
        : rowsLockedOnConstruction(false), colsLockedOnConstruction(false), 
          m(v.nrow()*esz), n(v.ncol()), 
          owner(true), writable(true), leadingDim(v.nrow()*esz),  
          data(v.size() ? new S[v.size()*esz] : 0) 
    {
        assert(esz > 0);
        rowsCurrentlyLocked = rowsLockedOnConstruction;
        colsCurrentlyLocked = colsLockedOnConstruction;

        copyFromCompatibleSourceView(esz,v,s);
    }      
                   
    ~DataDescriptor() { if (owner) delete[] data; }
    
    const S* getElt(int eltsz, int i, int j) const 
        { assert(i*eltsz<m&&j<n); return &data[eindx(eltsz,i,j)]; }
    S*       updElt(int eltsz, int i, int j)       
        { assert(i*eltsz<m&&j<n); assert(writable); return &data[eindx(eltsz,i,j)]; }
    
    const S& getScalar(int i, int j) const 
        { assert(i<m&&j<n); return data[sindx(i,j)]; }
    S&       updScalar(int i, int j)       
        { assert(i<m&&j<n); assert(writable); return data[sindx(i,j)]; }
    
    int nrowElt(int esz) const { return m/esz; } 
    int ncolElt(int esz) const { return n; }
    int sizeElt(int esz) const { return m*n/esz; }
        
    void resize(int sz, int mm, int nn) {
        if (mm==nrowElt(sz) && nn==ncolElt(sz)) return;
        assert(owner 
                && (!rowsCurrentlyLocked||mm==nrowElt(sz)) 
                && (!colsCurrentlyLocked||nn==ncolElt(sz)));
        delete[] data; m=mm*sz; n=nn; leadingDim=m;
        data = new S[m*n];

        #ifndef NDEBUG
            fillWithScalar(CNT<StdNumber>::getNaN());
        #endif
    }

    void resizeKeep(int sz, int mm, int nn) {
        if (mm==nrowElt(sz) && nn==ncolElt(sz)) return;
        assert(owner 
                && (!rowsCurrentlyLocked||mm==nrowElt(sz)) 
                && (!colsCurrentlyLocked||nn==ncolElt(sz)));
        S* oldData = data;
        const int oldM = m, oldN = n, oldLeadingDim = leadingDim;
        m=mm*sz; n=nn; leadingDim=m;
        data = new S[m*n];

        #ifndef NDEBUG
            fillWithScalar(CNT<StdNumber>::getNaN());
        #endif

        const int colsToCopy = std::min(n,oldN);
        const int rowsToCopy = std::min(m,oldM);
        for (int j=0; j<colsToCopy; ++j) {
            S* const colBegin = &data[(long)j*leadingDim];
            const S* const oldColBegin = &oldData[(long)j*oldLeadingDim];
            for (int i=0; i<rowsToCopy; ++i) 
                colBegin[i] = oldColBegin[i];
        }
        delete[] oldData;
    }

    void lockNRows() {rowsCurrentlyLocked=true;}
    void lockNCols() {colsCurrentlyLocked=true;}

    // These die if resizing wasn't allowed at construction.
    void unlockNRows() {
        assert(!rowsLockedOnConstruction);
        rowsCurrentlyLocked=rowsLockedOnConstruction;
    }
    void unlockNCols() {
        assert(!colsLockedOnConstruction);
        colsCurrentlyLocked=colsLockedOnConstruction;
    }
    
    /// Multiply every data element by a scalar. Conveniently this can be done
    /// ignoring the element structure. Note the Scalar is a *number*, that is,
    /// the negator<>, if any, has been removed and the sign adjusted accordingly.
    /// TODO: should be a BLAS call.
    void scaleBy(const typename CNT<S>::StdNumber& s) {
        for (int j=0; j<n; ++j) {
            S* const colBeg = &data[sindx(0,j)];
            for (int i=0; i<m; ++i)
                colBeg[i] *= s;
        }
    }
 
    /// Multiply one of our data elements by a scalar. These are element indices.
    void scaleEltBy(int sz, int i, int j, const typename CNT<S>::StdNumber& s) {
        scaleElement(sz,updElt(sz,i,j),s);
    }            
    /// Add an element (scalarwise) to one of our data elements. These are element indices.
    void addToOneElt(int sz, int i, int j, const S* eltp) {
        addToElement(sz,updElt(sz,i,j),eltp);
    }        
    /// Subtract an element (scalarwise) from one of our data elements. These are element indices.
    void subFromOneElt(int sz, int i, int j, const S* eltp) {
        subFromElement(sz,updElt(sz,i,j),eltp);
    }  
    /// Copy in one element onto one of ours. These are element indices.
    void copyInOneElt(int sz, int i, int j, const S* eltp) {
        copyElement(sz,updElt(sz,i,j),eltp);
    }
                        
    /// Fill every currently allocated element with repeated copies of a *scalar* value.
    /// This probably makes sense only when the value is 0 or NaN since it ignores
    /// the internal structure of the elements.
    void fillWithScalar(const typename CNT<S>::StdNumber& val);
    
    /// Fill every currently allocated element with a copy of the passed-in element;
    void fillWithElement(int sz, const S* eltp);

    /// Given a pointer to densely packed elements in row order, copy them into
    /// the currently allocated space. The size is determined by the current
    /// allocation; all hell will break loose if you don't provide (m/sz)Xn elements.
    void copyInElementsByRow(const int sz, const S* eltp) {
        assert((m%sz)==0);
        for (int j=0; j<n; ++j)
            copyInOneEltColumn(sz, j, eltp + j*sz, m);
    }

    /// Given a pointer to densely packed *elements* in column order (possibly with
    /// a leading dimension greater than m/sz elements), copy them into
    /// the currently allocated space. The size is determined by the current
    /// allocation; all hell will break loose if you don't provide (m/sz)*n elements.
    void copyInElementsByColumn(const int sz, const S* eltp, int eldim=0) {
        assert((m%sz)==0);
        const int sldim = eldim ? eldim*sz : m;  // get leading dimension in scalars
        for (int j=0; j<n; ++j)
            copyInOneEltColumn(sz, j, eltp + (long)j*sldim, 1);
    }
    
    /// Here we are given a pointer to a column's worth of *elements*, possibly
    /// with regular gaps between the elements. The idea is to copy them
    /// into a particular column of this data. Estride==1 means the source elements
    /// are packed together.    
    void copyInOneEltColumn(int sz, int j, const S* eltp, int estride=1);
    
    /// Here we have a pointer to a columns worth of *scalars*, possibly with
    /// gaps of (stride-1) scalars separating them. Copy them into our column j.
    void copyInOneScalarColumn(int j, const S* sp, int stride=1); 
    
    /// Here we have a pointer to a whole matrix of scalars, of the same
    /// dimension (allegedly) as this one. There may be an arbitrary
    /// leading dimension, however (given in scalars). We use the leading dimension
    /// to find the beginning of each source column and then copy it onto
    /// the corresponding column of our data.
    void copyInDataByColumn(const S* data, int ldim=0); // raw data copy
    
    // Single element manipulation: VERY SLOW, use sparingly.    
    static void copyElement(const int sz, S* dest, const S* src)
        { for (int k=0; k<sz; ++k) dest[k] = src[k]; }
    static void fillElement(const int sz, S* dest, const typename CNT<S>::StdNumber& src)
        { for (int k=0; k<sz; ++k) dest[k] = src; }   
    static void addToElement(const int sz, S* dest, const S* src)
        { for (int k=0; k<sz; ++k) dest[k] += src[k]; }        
    static void subFromElement(const int sz, S* dest, const S* src)
        { for (int k=0; k<sz; ++k) dest[k] -= src[k]; }
    static void scaleElement(const int sz, S* dest, const typename CNT<S>::StdNumber& s)
        { for (int k=0; k<sz; ++k) dest[k] *= s; }

    // Assuming this is a square, well-conditioned matrix and we have write access,
    // this will replace the data with the data corresponding to the inverse of
    // the matrix. This ignores element size and considers the data as scalar.
    void invertInPlace();

    bool hasContiguousData() const {return leadingDim==m;}
    long getContiguousDataLength() const {
        assert(hasContiguousData());
        return nScalars();
    }
    const S* getContiguousData() const {
        assert(hasContiguousData());
        return data;
    }

    S* updContiguousData() {
        assert(hasContiguousData());
        return data;
    }

    // Replace the current data pointer, which must point to contiguous,
    // elements, with the passed-in pointer. 
    void replaceContiguousData(S* newData, long length, 
                               bool takeOwnership) 
    {
        assert(hasContiguousData() && length == nScalars());
        if (owner) delete[] data;
        data = newData;
        owner = takeOwnership;
        writable = true;
    }

    // Replace the current data pointer, which must point to contiguous
    // elements, with the passed-in pointer to read-only data. We *do not*
    // take ownership of the new data, even if we previously owned our data.
    void replaceContiguousData(const S* newData, long length) 
    {
        replaceContiguousData(const_cast<S*>(newData), length, false);
        writable = false;
    }


    // Given a pointer to a heap allocation of 'length' contiguous scalars,
    // swap it with our own heap allocation of the same length, taking ownership
    // of the new data and relinquishing ownership of the old data.
    void swapOwnedContiguousData(S* newData, long length, S*& oldData) 
    {
        assert(owner && hasContiguousData() && length == nScalars());
        oldData = data;
        data = newData;
    }
                                    
private:    
    const bool   rowsLockedOnConstruction;     // can't ever change m
    const bool   colsLockedOnConstruction;     // can't ever change n 

    bool rowsCurrentlyLocked;   // lock always allowed, unlock only if unlocked
    bool colsCurrentlyLocked;   //   when first constructed.
          
    int  m, n;     // logical shape of data m=nrows, n=ncols, in *scalars*
    
    // Physical data. Currently this just a LAPACK-style dense matrix, in column
    // order, with a leading dimension >= m.
    bool         owner;
    bool         writable;      // if not owner, we may only have read-only access
    int          leadingDim;    // in *scalars*, not *elements*
    S*           data;

    // Although we restrict any single dimension to what will fit in an int
    // (about 2 billion scalars) we allow the total number of scalars m*n
    // to be larger on a 64 bit machine and return the value as a long
    // which should be big enough on 64 bit machines to address anything that will
    // fit in memory. Note that we're only calculating the number of meaningful
    // scalars; we're not including any dead space introduced by leadingDim.
    long nScalars() const { return (long)m*n; }
    
    // These are linear indices into the data and can thus be very large, so use
    // long here too.
    long sindx(int i, int j) const { 
        return (long)j*leadingDim + i; 
    }
    long eindx(int sz, int i, int j) const { 
        return (long)j*leadingDim + i*sz; 
    }
     
private:
    // XXX suppress assignment
    DataDescriptor& operator=(const DataDescriptor&);
};  

#ifdef NOTYET

    /////////////////////////////////////
    // MEMORY LAYOUT CLASS DECLARATION //
    /////////////////////////////////////

// This is the abstract interface that describes the services which
// must be performed by a concrete data descriptor. It is templatized
// by a *standard* number, i.e., float, double, long double, 
// complex<float>, complex<double>, complex<long double> and nothing
// else. Each concrete memory layout class is expected to provide
// routines which can work with all the others, but only when
// templatized by the same type.
template <class S>
class MemoryLayout {
    typedef NTraits<S>::ScalarSq EAbs;
public:
    MemoryLayout() 
      : m(0), n(0), owner(true), writable(true), 
        transposed(false), alloc(0), data(0) { }
    virtual ~MemoryLayout() { 
        char* dataAsChar = reinterpret_cast<char*>(data);
        data = 0;
        delete[] dataAsChar;
    }

    // Computational routines, where "this" is the result.

    // this = a*this. Guarantee: won't look at 'this' if a==0.
    void selfScale(const S& a) {
        if (a == S(0)) concreteSetSelfToZero();
        else           concreteSelfScale(a);
    }

    // this = a*this + b*Y.
    // Restriction: Y is same shape as 'this' (enforced here)
    // Restriction: 'this' can contain the result without loss (concrete enforces)
    // Guarantee: won't look at this if a==0 (handled here)
    // Guarantee: won't look at Y if b==0 (handled here)
    void axpby(const S& a, const S& b, const MemoryLayout& Y) {
        assert(isSameShape(Y));
        if      (b==S(0)) selfScale(a);
        else if (a==S(0)) concreteCopyScaled(b,Y);
        else              concreteAxpby(a,b,Y);
    }

    // this = a*this + b*X*Y
    // Restriction: this,X,Y have conforming shapes (enforced here)
    // Restriction: 'this' can contain the result without loss (concrete enforces)
    // Guarantee: won't look at this if a==0 (concrete handles)
    // Guarantee: won't look at Y if b==0 (handled here)
    void matmul(const S& a, const S& b, 
                const MemoryLayout& X, const MemoryLayout& Y) {
        assert(X.nrow==nrow && Y.ncol==ncol && X.ncol==Y.nrow);
        if (b==S(0)) selfScale(a);
        else concreteMatmul(a,b,X,Y);
    }

    // The matrix 1-norm is the largest column sum of abs values. (This is just
    // the largest absolute value in a RowVector.)
    EAbs matrixNorm1() const {
        if (isEmpty()) return CNT<EAbs>::getNaN();  // not defined
        EAbs maxColSumAbs;
        if (nrow==1) {
            (void)concreteFindRowMaxAbs(0, &maxColSumAbs);
            return maxColSumAbs;
        }

        maxColSumAbs = concreteColSumAbs(0);
        for (c=1; c < ncol; ++c) {
            const EAbs colSumAbs = concreteColSumAbs(c);
            if (colSumAbs > maxColSumAbs)
                maxColSumAbs = colSumAbs;
        }
        return maxColSumAbs;
    }

    // The matrix infinity-norm is the largest row sum of abs values. (This is
    // just the largest absolute value in a Vector.)
    EAbs matrixNormInf() const {
        if (isEmpty()) return CNT<EAbs>::getNaN();  // not defined
        EAbs maxRowSumAbs;
        if (ncol==1) {
            (void)concreteFindColMaxAbs(0, &maxRowSumAbs);
            return maxRowSumAbs;
        }

        maxRowSumAbs = concreteRowSumAbs(0);
        for (r=1; r < nrow; ++r) {
            const EAbs rowSumAbs = concreteRowSumAbs(r);
            if (rowSumAbs > maxRowSumAbs)
                maxRowSumAbs = rowSumAbs;
        }
        return maxRowSumAbs;
    }

    // Vector norms are defined only if the matrix has at least one dimension
    // which is 1. Unlike matrix norms, they are the same for 1xn or nx1.

    // The vector 1-norm is the sum of the absolute values of the elements.
    EAbs vectorNorm1() const {
        if (ncol==1 && nrow>0) return concreteColSumAbs(0);
        if (nrow==1 && ncol>0) return concreteRowSumAbs(0);
        return CNT<S>::getNaN();
    }

    // The vector infinity-norm is the element of maximum absolute value.
    EAbs vectorNormInf() const {
        EAbs mx;
        if      (ncol==1 && nrow>0) (void)concreteFindColMaxAbs(0,&mx);
        else if (nrow==1 && ncol>0) (void)concreteFindRowMaxAbs(0,&mx);
        else mx = CNT<S>::getNaN();
        return mx;
    }

    void findMaxAbs(long& row, long& col, S* v=0) const {
        if (isEmpty()) {
            row = col = -1;
            if (v) *v = CNT<S>::getNaN();
            return;
        }
        concreteFindMaxAbs(row,col,v);
    }

    long findColMaxAbs(long col, S* v=0) const {
        if (isEmpty()) {
            if (v) *v = CNT<S>::getNaN();
            return -1;
        }
        return concreteFindColMaxAbs(col,v);
    }

    long findRowMaxAbs(long row, S* v=0) const {
        if (isEmpty()) {
            if (v) *v = CNT<S>::getNaN();
            return -1;
        }
        return concreteFindRowMaxAbs(col,v);
    }

        // VIRTUAL METHODS FOR CONCRETE CLASSES TO IMPLEMENT

    // No need to optimize for a==0.
    virtual void concreteSelfScale(const S& a) {
        // default implementation
        Lapack::scal<S>(nalloc,a,data,1); 
    }

    virtual void concreteSetSelfToZero() {
        // default implementation
        for (long i=0; i<nalloc; ++i)
            data[i] = S(0);
    }

    // Copy Y into this. Will not be called if b==0 or Y has wrong shape.
    // Concrete implementation must not look at current contents of 'this'.
    // 'this' must be able to hold the result losslessly.
    virtual void concreteCopyScaled(const S& b, const MemoryLayout& Y) = 0;

    // this = Y, selecting only from elements represented in Y
    virtual void copyRepresentedElements(const MemoryLayout& Y) = 0;

    // this = a*this + b*Y.
    // Will not be called if a==0 or b==0 or Y has wrong shape.
    virtual void concreteAxpby(const S& a, const S& b, const MemoryLayout& Y) = 0;

    // this = a*this + b*X*Y
    // will not be called if b==0.
    // Concrete implementation must guarantee not to look at 'this' if a==0.
    virtual void concreteMatmul(const S& a, const S& b, 
                                const MemoryLayout& X, const MemoryLayout& Y) = 0;

    // A concrete data descriptor should override this if it can do better,
    // otherwise we'll cobble it together from the column operations.
    // This routine will NOT be called if the matrix is empty.
    virtual void concreteFindMaxAbs(long& row, long& col, S* v=0) const {
        EAbs maxAbs;
        col=0; row=concreteFindColMaxAbs(col, &maxAbs);
        for (long c=1; c < ncol; ++c) {
            EAbs colMaxAbs;
            const long r = concreteFindColMaxAbs(c, &colMaxAbs);
            if (colMaxAbs > maxAbs) {
                maxAbs = colMaxAbs;
                col=c; row=r;
            }
        }
        if (v) *v = maxAbs;
    }

    // SEARCH (returns index, optionally value also)
    // The concrete routines will NOT be called on an empty matrix.
    virtual long concreteFindColMaxAbs(long col, S* v=0) const = 0;
    virtual long concreteFindRowMaxAbs(long row, S* v=0) const = 0;

    virtual long concreteFindColMinAbs(long col, S* v=0) const = 0;
    virtual long concreteFindRowMinAbs(long row, S* v=0) const = 0;

    // NORMS

    virtual EAbs concreteColSumAbs(long col) const = 0;
    virtual EAbs concreteRowSumAbs(long row) const = 0;

    // Frobenius norm: sqrt(sum( e_ij^2 )) for every element in logical
    // matrix represented by 'this'. This is also the 2-norm for vectors,
    // but not for 2d matrices.
    virtual EAbs normFro() const = 0;

protected:
    // Always allocate chars to ensure no initialization.
    void allocateUninitializedScalars(long ns) {
        assert(data==0);
        nalloc = ns;
        char* dataAsChar = new char[nalloc * sizeof(S)];
        data = reinterpret_cast<S*>(dataAsChar);
        #ifndef NDEBUG
            for (int i=0; i<ns; ++i)
                data[i] = NTraits<S>::getNaN();
        #endif
        owner = writable = true;
    }

private:
    bool isEmpty() const {return nrow==0 || ncol==0;}
    bool isSameShape(const MemoryLayout& m) const {
        return m.nrow==nrow && m.ncol==ncol;
    }

    long nrow, ncol; // logical dimensions of matrix
    bool owner;
    bool writable;

    bool transposed; // data is stored as Hermitian transpose of logical matrix
    long nalloc;     // # scalars S allocated in data
    S*   data;
};

    ///////////////////////////////////////////////
    // CONCRETE MEMORY LAYOUT CLASS DECLARATIONS //
    ///////////////////////////////////////////////

template <class S>
class LapackTriangularInConventional : public MemoryLayout<S> {
    typedef MemoryLayout<S> Base;
public:
    LapackTriangularInConventional(); // TODO

    bool isEltStored(long i, long j) {
        return (i<j && upper) || (i>j && !upper) || (i==j && !assumeUnitDiagonal);
    }

    const S& getStoredElt(long i, long j) {
        assert(isEltStored(i,j));
        return data[sindx(i,j)];
    }
    const S& updStoredElt(long i, long j) {
        assert(isEltStored(i,j));
        return data[sindx(i,j)];
    }

    const S& getAnyElt(long i, long j) {
        static const S zero = S(0);
        static const S one  = S(1);
        if (upper) {
            if (i < j) return data[sindx(i,j)];
        } else { // lower
            if (i > j) return data[sindx(i,j)];
        }
        if (i==j) return assumeUnitDiagonal ? one : data[sindx(i,j)];
        return zero;
    }
    
private:
    long sindx(long i, long j) const {return j*leadingDimension+i;}
    bool upper;
    bool assumeUnitDiagonal;
    long leadingDimension;
};

/**
 * This is the description of Lapack conventional matrix storage.
 * This is column-oriented storage with a leading dimension, suited
 * for passing directly to Lapack & Blas routines.
 */
template <class S> 
class LapackConventional : public MemoryLayout<S> {
    typedef MemoryLayout<S> Base;
public:
    LapackConventional(long nr, long nc)
      : Base(nr,nc) // owner, writable, resizeable
    {
        data = allocateUninitializedScalars(nr*nc);
        leadingDimension = nr;
    }

    LapackConventional(long nr, long nc, long lda, const S* d)
      : Base(nr,nc,d) // not owner, not writable, not resizeable
    {
        assert(lda >= nr);
        leadingDimension = lda;
    }

    LapackConventional(long nr, long nc, long lda, S* d)
      : Base(nr,nc,d) // not owner, writable, not resizeable
    {
        assert(lda >= nr);
        leadingDimension = lda;
    }

    // Map i,j scalar indices to data index
    long sindx(long i, long j) const {
        return j*leadingDimension + i;
    }

private:
    long leadingDimension;
};

/**
 * This is the description of a Lapack 1-d vector. It provides for
 * a starting address and a stride (1 means packed). You can feed the
 * raw data directly to BLAS routines.
 *
 * This is agnostic about whether the matrix itself is a column or
 * row vector.
 */
template <class S>
class LapackVector : public MemoryLayout<S> {
    typedef MemoryLayout<S> Base;
public:
    explicit LapackVector(long nr, long nc)
      : Base(nr,nc) // owner, writable, resizeable
    {
        data = allocateUninitializedScalars(nr*nc);
        leadingDimension = nr;
    }

    LapackVector(long nr, long nc, long lda, const S* d)
      : Base(nr,nc,d) // not owner, not writable, not resizeable
    {
        assert(lda >= nr);
        leadingDimension = lda;
    }

    LapackConventional(long nr, long nc, long lda, S* d)
      : Base(nr,nc,d) // not owner, writable, not resizeable
    {
        assert(lda >= nr);
        leadingDimension = lda;
    }

    // Map i,j scalar indices to data index
    long sindx(long i, long j) const {
        return data + j*leadingDimension + i;
    }
private:
    long stride;
};

#endif
   
} // namespace SimTK   


#endif // SimTK_SIMMATRIX_DATA_DESCRIPTOR_H_
