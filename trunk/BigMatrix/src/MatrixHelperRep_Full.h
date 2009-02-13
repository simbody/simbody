#ifndef SimTK_SimTKCOMMON_MATRIX_HELPER_REP_FULL_H_
#define SimTK_SimTKCOMMON_MATRIX_HELPER_REP_FULL_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmatrix(tm)                       *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-9 Stanford University and the Authors.         *
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

#include "SimTKcommon/Scalar.h"
#include "SimTKcommon/SmallMatrix.h"

#include "MatrixHelperRep.h"

#include <cstddef>

namespace SimTK {


//------------------------------- FullHelper -----------------------------------
//
// This abstract class represents a matrix for which every one
// of the mXn elements is stored explicitly in data. 
// Derived classes provide fast scalar elements, column- or row-order storage,
// or fancier indexing schemes.
//------------------------------------------------------------------------------
template <class S>
class FullHelper : public MatrixHelperRep<S> {
    typedef FullHelper<S>           This;
    typedef MatrixHelperRep<S>      Base;
    typedef FullHelper<SNeg>        ThisNeg;
    typedef FullHelper<SHerm>       ThisHerm;
public:
    // Allocate new memory for a full, contiguous storage, owner matrix. The
    // leading dimension will be either nr or nc depending on whether
    // this is column- or row-order storage. Note that nr and nc are
    // in elements while leadingDim is in scalars.
    This(int esz, int cppesz, int nr, int nc, int ldim)
    :   Base(esz,cppesz), m_leadingDim(ldim)
    {
        assert(m_leadingDim==nr*m_eltSize || m_leadingDim==nc*m_eltSize);
        m_owner     = true;
        m_writable  = true;
        allocateData(nr,nc);
        m_actual.setStructure(MatrixStructure::Full);
        m_actual.setActualSize(nr,nc);
    }

    // Use someone else's memory, which we assume to be the right size. 
    This(int esz, int cppesz, int nr, int nc, int ldim, const S* shared, bool canWrite)
    :   Base(esz,cppesz), m_leadingDim(ldim)
    {
        m_owner     = false;
        m_writable  = canWrite;
        setData(const_cast<S*>(shared));
        m_actual.setStructure(MatrixStructure::Full);
        m_actual.setActualSize(nr,nc);
    }

    int getLeadingDim() const {return m_leadingDim;}

    void getAnyElt_(int i, int j, S* value) const 
    {   copyElt(value, getElt_(i,j)); }

    // The meaning of "Full" is that every element is stored in memory somewhere.
    bool eltIsStored_(int, int) const {return true;}

    // Just changing the return type here.
    virtual FullHelper* cloneHelper_() const = 0;

    // A deep copy of a full matrix always results in contiguous storage, in
    // row or column order depending on what is cheaper.
    virtual FullHelper* createDeepCopy_() const = 0;

    // This serves for all block views of Full matrices because selecting
    // a block doesn't change the type of handle we need.
    FullHelper* createBlockView_(const EltBlock& block) {
        FullHelper* p = cloneHelper_();
        p->m_data = updElt_(block.row0(), block.col0());
        return p;
    }

protected:
    int m_leadingDim; // in scalars

    // Note that these indexers can take signed indices and produce signed
    // data offsets. Here we don't know whether the data is stored in row
    // or column order, so the indices are in terms of "fast" and "slow"
    // where is "fast" is the one that moves consecutively through memory.

    // First is for use with composite elements, second is for scalar elements.
    ptrdiff_t eltIx   (int fast, int slow) const {return (ptrdiff_t)slow*m_leadingDim + fast*m_eltSize;}
    ptrdiff_t scalarIx(int fast, int slow) const {return (ptrdiff_t)slow*m_leadingDim + fast;}

    bool isContiguousElt   (int nFast) const {return m_leadingDim == nFast*m_eltSize;}
    bool isContiguousScalar(int nFast) const {return m_leadingDim == nFast;}

    S scalarColSum(int j) const {
        S csum = S(0);
        for (int i=0; i<nrow(); ++i) csum += *getElt_(i,j);
        return csum;
    }

    S scalarRowSum(int i) const {
        S rsum = S(0);
        for (int j=0; j<ncol(); ++j) rsum += *getElt_(i,j);
        return rsum;
    }
};

//----------------------------- RegularFullHelper ------------------------------
//
// This is a Full matrix whose elements are regularly spaced in memory, meaning
// that a change in row index i creates a predictable memory offset independent
// of j, and similarly a change in column index j creates a predictable offset
// independent of i.
//------------------------------------------------------------------------------
template <class S>
class RegularFullHelper : public FullHelper<S> {
    typedef RegularFullHelper<S>    This;
    typedef FullHelper<S>           Base;
public:
    This(int esz, int cppesz, int nr, int nc, int ldim)
    :   Base(esz,cppesz,nr,nc,ldim) {}
    This(int esz, int cppesz, int nr, int nc, int ldim, const S* shared, bool canWrite)
    :   Base(esz,cppesz,nr,nc,ldim,shared,canWrite) {}

    bool hasRegularData_() const {return true;}

    // These implementations work for any RegularFullHelper.
    This*            createRegularView_(const EltBlock&, const EltIndexer&);
    VectorHelper<S>* createDiagonalView_();
    VectorHelper<S>* createColumnView_(int j, int i, int m);
    VectorHelper<S>* createRowView_(int i, int j, int n);

    // This is a new pure virtual that all RegularFullHelpers must provide.
    // Some will return their actual indexers, others will calculate one,
    // which is why we don't return a reference here.
    virtual EltIndexer getEltIndexer() const = 0;

protected:
    // This is for use by scalar col- and row- ordered matrices only. The same
    // code works because transpose(inv(m))==inv(transpose(m)).
    void lapackInvertInPlace() {
        assert(m_eltSize==1 && nrow()==ncol()); // should have been checked already
        const int m = nrow();
        StdNumber* rawMem = reinterpret_cast<StdNumber*>(m_data);
        std::vector<int> ipiv(m);
        int info;
        Lapack::getrf<StdNumber>(m,m,rawMem,m_leadingDim,&ipiv[0],info);
        assert(info==0);

        // Calculate optimal size for work
        StdNumber workSz;
        Lapack::getri<StdNumber>(m,rawMem,m_leadingDim,&ipiv[0],&workSz,-1,info);
        const int wsz = (int)CNT<StdNumber>::real(workSz);

        std::vector<StdNumber> work(wsz);
        Lapack::getri<StdNumber>(m,rawMem,m_leadingDim,&ipiv[0],&work[0],wsz,info);
        assert(info==0);
    }
};

// Full, column order, composite. Row i is the fast dimension here.
template <class S>
class FullColOrderEltHelper : public RegularFullHelper<S> {
    typedef FullColOrderEltHelper<S>    This;
    typedef RegularFullHelper<S>        Base;
public:
    // Leading dimension is # rows for contiguous column-ordered memory.
    This(int esz, int cppesz, int nr, int nc)
    :   Base(esz, cppesz, nr, nc, nr*esz) {   
        m_actual.setStorage(MatrixStorage(MatrixStorage::Full, MatrixStorage::NoPlacement,
                                          MatrixStorage::ColumnOrder, MatrixStorage::NoDiag)); 
    }

    This(int esz, int cppesz, int nr, int nc, int ldim, S* shared, bool canWrite)
    :   Base(esz, cppesz, nr, nc, ldim, shared, canWrite) {
        assert(ldim>=nr*esz);
        m_actual.setStorage(MatrixStorage(MatrixStorage::Full, MatrixStorage::NoPlacement,
                                          MatrixStorage::ColumnOrder, MatrixStorage::NoDiag));
    }

    bool preferRowOrder_() const {return false;}

    // These virtual methods should be overridden in the derived scalar class.
    virtual const S*    getElt_(int i, int j) const {return m_data + eltIx(i,j);}
    virtual S*          updElt_(int i, int j)       {return m_data + eltIx(i,j);}
    virtual bool        hasContiguousData_()  const {return isContiguousElt(nrow());}
    virtual This*       cloneHelper_()        const {return new This(*this);}

    // This implementation will return a FullRowOrderEltHelper. Override for scalars.
    virtual RegularFullHelper<S>* createTransposeView_();

    // These implementations are fine for both composite and scalar class.

    // Regular spacing as (dfast/di, dfast/dj) and (dslow/di, dslow/dj), with i,j in
    // elements.
    EltIndexer getEltIndexer() const {return EltIndexer(1,0,0,1);}

    This* createDeepCopy_() const {
        This* p = cloneHelper_();
        p->m_writable = true;
        p->m_owner = true;
        p->allocateData(nelt());
        if (hasContiguousData_())
            std::memcpy(p->m_data, m_data, nelt()*m_eltSize*sizeof(S));
        else for (int j=0; j < ncol(); ++j)
            std::memcpy(p->updElt_(0,j), getElt_(0,j), nrow()*m_eltSize*sizeof(S));
        return p;
    }


    // OK for any size elements.
    void resize_(int m, int n) {
        clearData();
        allocateData(m,n);
        m_leadingDim = m * m_eltSize; // number of scalars in a column
    }

    // OK for any size elements.
    void resizeKeep_(int m, int n) {
        const int newLeadingDim = m * m_eltSize; // number of scalars in a column
        S* const newData = allocateMemory(m,n);
        const int colsToCopy = std::min(n, ncol());
        const int rowsToCopy = std::min(m, nrow()); // in elements
        for (int j=0; j < colsToCopy; ++j) {
            S*       const dest = newData + (ptrdiff_t)j*newLeadingDim;
            const S* const src  = m_data  + (ptrdiff_t)j*m_leadingDim;
            std::memcpy(dest, src, rowsToCopy*m_eltSize*sizeof(S));
        }
        clearData();
        setData(newData);
        m_leadingDim = newLeadingDim;
    }
};

// Full, column order, scalar, final.
template <class S>
class FullColOrderScalarHelper : public FullColOrderEltHelper<S> {
    typedef FullColOrderScalarHelper<S> This;
    typedef FullColOrderEltHelper<S>    Base;
public:
    // Leading dimension is # rows for contiguous column-ordered memory.
    This(int nr, int nc)
    :   Base(1, 1, nr, nc) {}
    This(int nr, int nc, int ldim, S* shared, bool canWrite)
    :   Base(1, 1, nr, nc, ldim, shared, canWrite) {}

    // For speed, these override the Base implementations for composite elements.
    const S*    getElt_(int i, int j) const {return m_data + scalarIx(i,j);}
    S*          updElt_(int i, int j)       {return m_data + scalarIx(i,j);}
    bool        hasContiguousData_()  const {return isContiguousScalar(nrow());}
    This*       cloneHelper_()        const {return new This(*this);}

    // This implementation will return a FullRowOrderScalarHelper.
    RegularFullHelper<S>* createTransposeView_();

    void colSum_(int j, S* csum) const {*csum = scalarColSum(j);}
    void rowSum_(int i, S* rsum) const {*rsum = scalarRowSum(i);}
    // Sum element column by column to avoid cache faults.
    void sum_(S* esum) const {
        *esum = S(0);
        for (int j=0; j<ncol(); ++j) *esum += scalarColSum(j);
    }

    void invertInPlace_() {lapackInvertInPlace();}
};

// Full, row order, composite. Column j is the fast dimension here.
template <class S>
class FullRowOrderEltHelper : public RegularFullHelper<S> {
    typedef FullRowOrderEltHelper<S>    This;
    typedef RegularFullHelper<S>        Base;
public:
    // Leading dimension is # cols for contiguous row-ordered memory.
    This(int esz, int cppesz, int nr, int nc)
    :   Base(esz, cppesz, nr, nc, nc*esz) {   
        m_actual.setStorage(MatrixStorage(MatrixStorage::Full, MatrixStorage::NoPlacement,
                                          MatrixStorage::RowOrder, MatrixStorage::NoDiag)); 
    }

    This(int esz, int cppesz, int nr, int nc, int ldim, S* shared, bool canWrite)
    :   Base(esz, cppesz, nr, nc, ldim, shared, canWrite) {
        assert(ldim>=nc*esz);
        m_actual.setStorage(MatrixStorage(MatrixStorage::Full, MatrixStorage::NoPlacement,
                                          MatrixStorage::RowOrder, MatrixStorage::NoDiag));
    }

    bool preferRowOrder_() const {return true;}

    // These virtual methods should be overridden in the derived scalar class.
    virtual const S*    getElt_(int i, int j) const {return m_data + eltIx(j,i);}
    virtual S*          updElt_(int i, int j)       {return m_data + eltIx(j,i);}
    virtual bool        hasContiguousData_()  const {return isContiguousElt(ncol());}
    virtual This*       cloneHelper_()        const {return new This(*this);}

    // This implementation will return a FullColOrderEltHelper. Override for scalars.
    virtual RegularFullHelper<S>* createTransposeView_();

    // These implementations are fine for both composite and scalar class.
    EltIndexer getEltIndexer() const {return EltIndexer(0,1,1,0);}

    This* createDeepCopy_() const {
        This* p = cloneHelper_();
        p->m_writable = true;
        p->m_owner = true;
        p->allocateData(nelt());
        if (hasContiguousData_())
            std::memcpy(p->m_data, m_data, nelt()*m_eltSize*sizeof(S));
        else for (int i=0; i < nrow(); ++i)
            std::memcpy(p->updElt_(i,0), getElt_(i,0), ncol()*m_eltSize*sizeof(S));
        return p;
    }

    void resize_(int m, int n) {
        clearData();
        allocateData(m,n);
        m_leadingDim = n * m_eltSize;   // number of scalars in a row
    }

    // OK for any size elements.
    void resizeKeep_(int m, int n) {
        const int newLeadingDim = n * m_eltSize; // number of scalars in a row
        S* const newData = allocateMemory(m,n);
        const int colsToCopy = std::min(n, ncol()); // in elements
        const int rowsToCopy = std::min(m, nrow());
        for (int i=0; i < rowsToCopy; ++i) {
            S*       const dest = newData + (ptrdiff_t)i*newLeadingDim;
            const S* const src  = m_data  + (ptrdiff_t)i*m_leadingDim;
            std::memcpy(dest, src, colsToCopy*m_eltSize*sizeof(S));
        }
        clearData();
        setData(newData);
        m_leadingDim = newLeadingDim;
    }
};

// Full, row order, scalar.
template <class S>
class FullRowOrderScalarHelper : public FullRowOrderEltHelper<S> {
    typedef FullRowOrderScalarHelper<S> This;
    typedef FullRowOrderEltHelper<S>    Base;
public:
    // Leading dimension is # cols for contiguous row-ordered memory.
    This(int nr, int nc)
    :   Base(1, 1, nr, nc) {}
    This(int nr, int nc, int ldim, S* shared, bool canWrite)
    :   Base(1, 1, nr, nc, ldim, shared, canWrite) {}

    // For speed, these override the Base implementations for composite elements.
    const S*    getElt_(int i, int j) const {return m_data + scalarIx(j,i);}
    S*          updElt_(int i, int j)       {return m_data + scalarIx(j,i);}
    bool        hasContiguousData_()  const {return isContiguousScalar(ncol());}
    This*       cloneHelper_()         const {return new This(*this);}

    // This implementation will return a FullColOrderScalarHelper.
    RegularFullHelper<S>* createTransposeView_();

    void colSum_(int j, S* csum) const {*csum = scalarColSum(j);}
    void rowSum_(int i, S* rsum) const {*rsum = scalarRowSum(i);}
    // Sum element row by row to avoid cache faults.
    void sum_(S* esum) const {
        *esum = S(0);
        for (int i=0; i<nrow(); ++i) *esum += scalarRowSum(i);
    }

    void invertInPlace_() {lapackInvertInPlace();}
};

template <class S> inline RegularFullHelper<S>*
FullColOrderEltHelper<S>::createTransposeView_() {
    FullRowOrderEltHelper<S>* p = 
        new FullRowOrderEltHelper<S>(m_eltSize,m_cppEltSize, ncol(), nrow(), m_leadingDim, m_data, m_writable);
    return p;
}

template <class S> inline RegularFullHelper<S>*
FullColOrderScalarHelper<S>::createTransposeView_() {
    FullRowOrderScalarHelper<S>* p = 
        new FullRowOrderScalarHelper<S>(ncol(), nrow(), m_leadingDim, m_data, m_writable);
    return p;
}

template <class S> inline RegularFullHelper<S>*
FullRowOrderEltHelper<S>::createTransposeView_() {
    FullColOrderEltHelper<S>* p = 
        new FullColOrderEltHelper<S>(m_eltSize,m_cppEltSize, ncol(), nrow(), m_leadingDim, m_data, m_writable);
    return p;
}

template <class S> inline RegularFullHelper<S>*
FullRowOrderScalarHelper<S>::createTransposeView_() {
    FullColOrderScalarHelper<S>* p = 
        new FullColOrderScalarHelper<S>(ncol(), nrow(), m_leadingDim, m_data, m_writable);
    return p;
}

// Full, indexed, composite.
template <class S>
class FullIndexedEltHelper : public RegularFullHelper<S> {
    typedef FullIndexedEltHelper<S>     This;
    typedef RegularFullHelper<S>        Base;
public:
    // Construct a new view of any regularly-spaced full matrix, allowing
    // composite elements. Source is taken non-const but writability is only
    // granted if the source had it.
    // Note that the (r0,c0) element may be one element off the end of the
    // matrix as long as the corresponding dimension is 0.
    This(RegularFullHelper<S>& src, const EltBlock& block, const EltIndexer& ix)
    :   Base(src.getElementSize(), src.getCppElementSize(), 
             block.nrow(), block.ncol(), src.getLeadingDim(), 
             src.getElt(block.row0(), block.col0()), src.isWritable()), 
        m_indexer(src.getEltIndexer().postIndexBy(ix))
    {
        SimTK_SIZECHECK(block.row0(), src.nrow()-block.nrow(), "FullIndexedEltHelper ctor");
        SimTK_SIZECHECK(block.col0(), src.ncol()-block.ncol(), "FullIndexedEltHelper ctor"); 
    }

        // Implementations good for both composite and scalar elements. //

    bool        hasContiguousData_()  const {return false;}
    EltIndexer  getEltIndexer()       const {return m_indexer;}

    // Prefer to go through the data a row at a time if adjacent column
    // elements in the same row are closer together than adjacent row
    // elements in the same column.
    bool preferRowOrder_() const {return ixEltIx(0,1) < ixEltIx(1,0);}

        // Implementations that should be overridden for scalar elements. //

    virtual const S*    getElt_(int i, int j) const {return m_data + ixEltIx(i,j);}
    virtual S*          updElt_(int i, int j)       {return m_data + ixEltIx(i,j);}
    virtual This*       cloneHelper_()        const {return new This(*this);}

    // This implementation works for this class and the scalar derived class.
    This* createTransposeView_() {
        This* p = cloneHelper_();
        p->m_indexer = m_indexer.transpose();
        p->m_actual.setActualSize(ncol(), nrow());
        return p;
    }

    // A deep copy of an indexed matrix eliminates the index, producing a
    // full, column-ordered matrix with contiguous storage. 
    // TODO: should create col- or row-order depending on which way is
    // fastest to travel through the indexed matrix.
    virtual RegularFullHelper<S>* createDeepCopy_() const {
        if (preferRowOrder_()) {
            FullRowOrderEltHelper<S>* p = 
                new FullRowOrderEltHelper<S>(m_eltSize,m_cppEltSize,nrow(),ncol());
            for (int i=0; i < nrow(); ++i) {
                S* dest = p->updElt_(i,0);   // start of a dense row
                for (int j=0; j < ncol(); ++j, dest += m_eltSize)
                    copyElt(dest, getElt_(i,j));
            }
            return p;
        } else {
            FullColOrderEltHelper<S>* p = 
                new FullColOrderEltHelper<S>(m_eltSize,m_cppEltSize,nrow(),ncol());
            for (int j=0; j < ncol(); ++j) {
                S* dest = p->updElt_(0,j);   // start of a dense column
                for (int i=0; i < nrow(); ++i, dest += m_eltSize)
                    copyElt(dest, getElt_(i,j));
            }
            return p;
        }
    }

protected:
    EltIndexer m_indexer;

    int row(int i, int j) const {return m_indexer.row(i,j);}
    int col(int i, int j) const {return m_indexer.col(i,j);}

    ptrdiff_t ixEltIx(int i, int j) const {return eltIx(row(i,j),col(i,j));}
};

// Full, indexed, scalar, final.
template <class S>
class FullIndexedScalarHelper : public FullIndexedEltHelper<S> {
    typedef FullIndexedScalarHelper<S>  This;
    typedef FullIndexedEltHelper<S>     Base;
public:
    // Construct a new view of any regularly-spaced full matrix, as long as that
    // matrix has scalar elements.
    This(RegularFullHelper<S>& src, const EltBlock& block, const EltIndexer& ix)
    :   Base(src, block, ix)
    {
        SimTK_ASSERT_ALWAYS(src.getElementSize()==1, 
            "FullIndexedScalarHelper ctor: source must have scalar elements");
    }

    // For speed, these override the composite-element implementatiosn.
    const S*    getElt_(int i, int j) const {return m_data + ixScalarIx(i,j);}
    S*          updElt_(int i, int j)       {return m_data + ixScalarIx(i,j);}
    This*       cloneHelper_()         const {return new This(*this);}

    RegularFullHelper<S>* createDeepCopy_() const {
        if (preferRowOrder_()) {
            FullRowOrderScalarHelper<S>* p = 
                new FullRowOrderScalarHelper<S>(nrow(),ncol());
            for (int i=0; i < nrow(); ++i) {
                S* dest = p->updElt_(i,0);   // start of a dense row
                for (int j=0; j < ncol(); ++j, ++dest)
                    *dest = *getElt_(i,j);
            }
            return p;
        } else {
            FullColOrderScalarHelper<S>* p = 
                new FullColOrderScalarHelper<S>(nrow(),ncol());
            for (int j=0; j < ncol(); ++j) {
                S* dest = p->updElt_(0,j);   // start of a dense column
                for (int i=0; i < nrow(); ++i, ++dest)
                    *dest = *getElt_(i,j);
            }
            return p;
        }
    }

private:
    ptrdiff_t ixScalarIx(int i, int j) const {return scalarIx(row(i,j),col(i,j));}
};

template <class S> inline RegularFullHelper<S>*
RegularFullHelper<S>::createRegularView_(const EltBlock& block, const EltIndexer& ix) {
    RegularFullHelper<S>* p;
    if (m_eltSize==1) p = new FullIndexedScalarHelper<S>(*this, block, ix);
    else              p = new FullIndexedEltHelper<S>(*this, block, ix);
    return p;
}




} // namespace SimTK   


#endif // SimTK_SimTKCOMMON_MATRIX_HELPER_REP_FULL_H_
