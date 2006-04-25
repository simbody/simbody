#ifndef SIMTKIMPL_SIMMATRIX_BIGMATRIXIMPL_H_
#define SIMTKIMPL_SIMMATRIX_BIGMATRIXIMPL_H_

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

/**@file
 * These are declarations for some internal classes needed for 
 * implementation of the BigMatrix objects.
 */

#include "SimTKcommon.h"
#include "simmatrix/internal/common.h"
#include "simmatrix/internal/Scalar.h"
#include "simmatrix/internal/SmallMatrix.h"
#include "simmatrix/internal/BigMatrix.h"


namespace simtkimpl {

using namespace SimTK;

/**
 * Describe a subset of the elements of a MatrixDataImpl object referenced by
 * the same MatrixHelper that references this MatrixViewImpl object. Note that
 * this is done entirely in terms of "elements" composed of one or more
 * scalars. The MatrixDataImpl object contains only scalars. It is guaranteed
 * that the storage for an element is composed of consecutive scalars, but 
 * otherwise we don't know anything about them. Note that different views of
 * the same data can claim elements of different sizes as long as the underlying
 * scalar types match.
 * 
 * A MatrixBase object does not have to contain a MatrixViewImpl object if it
 * permits unfettered access to all the data elements. In that case all MatrixBase 
 * operations are forward directly to the MatrixDataImpl along with the
 * element size we think we're looking at; otherwise the MatrixViewImpl
 * object must serve as an intermediary. 
 * 
 * This object expects to be intimately coupled with a MatrixDataImpl object
 * via the MatrixHelper. Don't move it around separately!
 */
class MatrixViewImpl {
public:
    /**
     * This class handles the mundane details of mapping from logical, element-oriented
     * view indices to physical but still element-oriented indices of the 
     * MatrixDataImpl object, which handles the final mapping of indices to 
     * memory address.
     */        
    class Indexer {
    public:
        // All arguments refer to *elements*, not *scalars*.
        Indexer(int r, int c, 
                ptrdiff_t drdx_=1, ptrdiff_t drdy_=0, ptrdiff_t dcdx_=0, ptrdiff_t dcdy_=1)
            : r0(r), c0(c), 
              drdx(drdx_), drdy(drdy_), dcdx(dcdx_), dcdy(dcdy_) { }
 
        // Compose a physical indexer (old) with a relative one (that is, an
        // indexer relative to this view) to produce a new physical one 
        // suitable for use on the original data.          
        Indexer(const Indexer& old, const Indexer& ix)
            : r0(old.row(ix.r0,ix.c0)), c0(old.col(ix.r0,ix.c0)),
              drdx(old.drdx*ix.drdx + old.drdy*ix.dcdx), 
              drdy(old.drdy*ix.dcdy + old.drdx*ix.drdy),
              dcdx(old.dcdx*ix.drdx + old.dcdy*ix.dcdx), 
              dcdy(old.dcdy*ix.dcdy + old.dcdx*ix.drdy) { }

        // Given element index (x,y) relative to this view, return element
        // indices (row,col) from which the data object can retrieve the
        // desired element.        
        int row(int x, int y) const { return r0 + drdx*x + drdy*y; }
        int col(int x, int y) const { return c0 + dcdx*x + dcdy*y; }

        // Make an indexer just like this one but with the roles of x and y reversed.
        // Still has same (0,0) element.
        Indexer transpose() {
            return Indexer(r0,c0,drdy,drdx,dcdy,dcdx);
        }
        
    private:             
        int r0,c0;        // indices of (0,0) element
        int drdx, drdy;   // row selection
        int dcdx, dcdy;   // column selection
        
        // no default construction
        Indexer();
    }; 

    // Like copy constructor but can *reduce* writability. Can't create writability
    // where none was permitted before, however.    
    MatrixViewImpl(MatrixViewImpl& v, bool wrt)
        : writable(v.writable && wrt), nr(v.nr), nc(v.nc), indexer(v.indexer) { }
        
    // Copy constructor takes a const MatrixViewImpl and loses writability.
    MatrixViewImpl(const MatrixViewImpl& v)
        : writable(false), nr(v.nr), nc(v.nc), indexer(v.indexer) { }


    MatrixViewImpl(bool wrt, int m, int n, const Indexer& ix)
        : writable(wrt), nr(m), nc(n), indexer(ix) { }
        
    // Offset constructor -- combine an old one and new instructions to get
    // another MatrixViewImpl suitable for use with the original data.
    MatrixViewImpl(const MatrixViewImpl& old, bool wrt, int m, int n,
                   const Indexer& ix)
        : writable(wrt), nr(m), nc(n), indexer(old.indexer, ix) { } 
        
    int nrow() const { return nr; }
    int ncol() const { return nc; } 
    size_t size() const { return nr*nc; }    
        
    bool isViewWritable() const { return writable; }
    
    int r(int i, int j) const { assert(i<nr&&j<nc); return indexer.row(i,j); }
    int c(int i, int j) const { assert(i<nr&&j<nc); return indexer.col(i,j); }
    void  rc(int i, int j, int& r, int& c) const 
        { assert(i<nr&&j<nc); r=indexer.row(i,j); c=indexer.col(i,j); }   
                       
private:
    bool    writable;           // does this view allow writing?
    int     nr,nc;              // logical shape
    Indexer indexer;            // view->data mapping
};

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
 * that is considered a matter for the View.
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
template <class S> class MatrixDataImpl {
    typedef typename CNT<S>::TReal      TReal;
    typedef typename CNT<S>::TNeg       TNeg;
    typedef typename CNT<S>::THerm      THerm;
    typedef typename CNT<S>::Number     Number;     // i.e., without negator<>
    typedef typename CNT<S>::StdNumber  StdNumber;  // i.e., without conjugate<>
public:
    /// Construct an empty (0x0), fully resizable owner.
    MatrixDataImpl() 
        : rowsLocked(false), colsLocked(false), m(0), n(0), 
          owner(true), writable(true), leadingDim(0), data(0) 
    { 
    }
 
    /// Construct and pre-allocate a resizable (mXsz)Xn owner with optional locking of
    /// the number of rows or columns (or both) at their initial values. No
    /// heap allocation occurs if one of the dimensions is 0, as is common
    /// for initial creation of a row or column.       
    MatrixDataImpl(int sz, int mm, int nn, bool mlock, bool nlock)
        : rowsLocked(mlock), colsLocked(nlock), m(mm*sz), n(nn), 
          owner(true), writable(true), leadingDim(mm*sz),  
          data(mm*nn*sz ? new S[mm*nn*sz] : 0)
    {
        assert(sz > 0 && mm >= 0 && nn >= 0);
    #ifndef NDEBUG
        fillWithScalar(CNT<StdNumber>::getNaN());
    #endif
    }

    /// Create a read-only descriptor for someone else's data. We're assuming it is 
    /// organized as nn columns of mm elements of size sz scalars, with leading
    /// dimension ldim *elements*. This is NOT resizable.
    MatrixDataImpl(int sz, int mm, int nn, int ldim, const S* s)
        : rowsLocked(true), colsLocked(true), m(mm*sz), n(nn), 
          owner(false), writable(false), leadingDim(ldim*sz),  
          data(const_cast<S*>(s))
    { 
        assert(sz > 0 && mm >= 0 && nn >= 0 && ldim >= 0);
    }
          
    /// Shared data constructor allowing writing. 
    MatrixDataImpl(int sz, int mm, int nn, int ldim, S* s)
        : rowsLocked(true), colsLocked(true), m(mm*sz), n(nn), 
          owner(false), writable(true), leadingDim(ldim*sz),  
          data(s) 
    { 
        assert(sz > 0 && mm >= 0 && nn >= 0 && ldim >= 0);
    }
         
    /// Caution: this is a *viewless* deep copy of the whole of the source. We'll 
    /// inherit rowsLocked, etc. but of course we are now the owner regardless of
    /// whether the source was, and leadingDim becomes just the number of rows.
    /// We also have write access since this is new data.
    MatrixDataImpl(const MatrixDataImpl& s)
        : rowsLocked(s.rowsLocked), colsLocked(s.colsLocked), m(s.m), n(s.n), 
          owner(true), writable(true), leadingDim(s.m),
          data(s.m*s.n ? new S[s.m*s.n] : 0) 
    {
        copyInDataByColumn(s.data,s.leadingDim); 
    }    

    /// Construction of new "viewless" data from a view of existing data. We 
    /// will NOT inherit rowsLocked & colsLocked since the view might be
    /// transposing -- reset the locks if you care. Of course we are the
    /// owner of the new (writable) data, which will be stored compactly regardless
    /// of how the original was stored.
    MatrixDataImpl(int esz, const MatrixViewImpl& v, const MatrixDataImpl& s)
        : rowsLocked(false), colsLocked(false), m(v.nrow()*esz), n(v.ncol()), 
          owner(true), writable(true), leadingDim(v.nrow()*esz),  
          data(v.size() ? new S[v.size()*esz] : 0) 
    {
        assert(esz > 0);

        // bad! XXX
        for (int j=0; j < v.ncol(); ++j)
            for (int i=0; i < v.nrow(); ++i)
                copyInOneElt(esz,i,j,s.getElt(esz,v.r(i,j), v.c(i,j))); 
    }      
                   
    ~MatrixDataImpl() { if (owner) delete data; }
    
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
        assert(owner && (!rowsLocked||mm==nrowElt(sz)) && (!colsLocked||nn==ncolElt(sz)));
        delete data; m=mm*sz; n=nn; leadingDim=m;
        data = new S[m*n];
    }
    
    /// Multiply every data element by a scalar. Conveniently this can be done
    /// ignoring the element structure. Note the Scalar is a *number*, that is,
    /// the negator<>, if any, has been removed and the sign adjusted accordingly.
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
            copyInOneEltColumn(sz, j, eltp + j*sldim, 1);
    }
    
    /// Here we are given a pointer to a column's worth of *elements*, possibly
    /// with regular gaps between the elements. The idea is to copy them
    /// into a particular column of this data. Stride==1 means the source elements
    /// are packed together.    
    void copyInOneEltColumn(int sz, int j, const S* eltp, int estride=1);
    
    /// Here we have a pointer to a columns worth of *scalars*, possibly with
    /// gaps of (stride-1) scalars separating them. Copy them into our column j.
    void copyInOneScalarColumn(int j, const S* sp, int stride=1); 
    
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
                                            
private:    
    const bool   rowsLocked;     // can't change m
    const bool   colsLocked;     // can't change n 
          
    int          m, n;     // logical shape of data m=nrows, n=ncols, in *elements*
    
    // Physical data. Currently this just a LAPACK-style dense matrix, in column
    // order, with a leading dimension >= m.
    bool         owner;
    bool         writable;    // if not owner, we may only have read-only access
    int          leadingDim;
    S*           data;
    
    // These are linear indices into the data and can thus be very large, so use
    // ptrdiff_t which is big enough on 64 bit machines to address everything.
    ptrdiff_t sindx(int i, int j)         const { return j*leadingDim + i; }
    ptrdiff_t eindx(int sz, int i, int j) const { return j*leadingDim + i*sz; }
     
private:
    // XXX suppress assignment
    MatrixDataImpl& operator=(const MatrixDataImpl&);
};       
   
} // namespace simtkimpl   


#endif //SIMTKIMPL_SIMMATRIX_BIGMATRIXIMPL_H_
