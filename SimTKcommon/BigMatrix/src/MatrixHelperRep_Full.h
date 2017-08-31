#ifndef SimTK_SimTKCOMMON_MATRIX_HELPER_REP_FULL_H_
#define SimTK_SimTKCOMMON_MATRIX_HELPER_REP_FULL_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
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

#include "SimTKcommon/Scalar.h"
#include "SimTKcommon/SmallMatrix.h"

#include "MatrixHelperRep.h"

#include <cstddef>
#include <cstring>

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
    typedef typename CNT<S>::TNeg   SNeg;
    typedef typename CNT<S>::THerm  SHerm;
    typedef FullHelper<SNeg>        ThisNeg;
    typedef FullHelper<SHerm>       ThisHerm;
public:
    // Allocate new memory for a full, contiguous storage, owner matrix. The
    // leading dimension will be either nr or nc depending on whether
    // this is column- or row-order storage. Note that nr and nc are
    // in elements while leadingDim is in scalars.
    FullHelper(int esz, int cppesz, int nr, int nc, int ldim)
    :   Base(esz,cppesz), m_leadingDim(ldim)
    {
        assert(   m_leadingDim==nr*this->m_eltSize 
               || m_leadingDim==nc*this->m_eltSize);
        // The "this->" (or Base::) is required here (by gcc and the standard, 
        // though not VC++ 9) to delay lookup of these non-tempatized members
        // until instantiation. (Google "two-stage name lookup".)
        this->m_owner     = true;
        this->m_writable  = true;
        this->allocateData(nr,nc);
        this->m_actual.setStructure(MatrixStructure::Full);
        this->m_actual.setActualSize(nr,nc);
    }

    // Use someone else's memory, which we assume to be the right size. 
    FullHelper(int esz, int cppesz, int nr, int nc, int ldim, 
               const S* shared, bool canWrite)
    :   Base(esz,cppesz), m_leadingDim(ldim)
    {
        this->m_owner     = false;
        this->m_writable  = canWrite;
        this->setData(const_cast<S*>(shared));
        this->m_actual.setStructure(MatrixStructure::Full);
        this->m_actual.setActualSize(nr,nc);
    }

    int getLeadingDim() const {return m_leadingDim;}

    void getAnyElt_(int i, int j, S* value) const 
    {   this->copyElt(value, this->getElt_(i,j)); }

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
        p->m_data = this->updElt_(block.row0(), block.col0());
        return p;
    }

protected:
    int m_leadingDim; // in scalars

    // Note that these indexers can take signed indices and produce signed
    // data offsets. Here we don't know whether the data is stored in row
    // or column order, so the indices are in terms of "fast" and "slow"
    // where "fast" is the one that moves consecutively through memory.

    // First is for use with composite elements, second is for scalar elements.
    ptrdiff_t eltIx   (int fast, int slow) const 
    {   return (ptrdiff_t)slow*m_leadingDim + fast*this->m_eltSize; }
    ptrdiff_t scalarIx(int fast, int slow) const 
    {   return (ptrdiff_t)slow*m_leadingDim + fast; }

    bool isContiguousElt   (int nFast) const 
    {   return m_leadingDim == nFast*this->m_eltSize; }
    bool isContiguousScalar(int nFast) const 
    {   return m_leadingDim == nFast; }

    S scalarColSum(int j) const {
        S csum = S(0);
        for (int i=0; i<this->nrow(); ++i) csum += *this->getElt_(i,j);
        return csum;
    }

    S scalarRowSum(int i) const {
        S rsum = S(0);
        for (int j=0; j<this->ncol(); ++j) rsum += *this->getElt_(i,j);
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
    RegularFullHelper(int esz, int cppesz, int nr, int nc, int ldim)
    :   Base(esz,cppesz,nr,nc,ldim) {}
    RegularFullHelper(int esz, int cppesz, int nr, int nc, int ldim, 
                      const S* shared, bool canWrite)
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
    typedef typename CNT<S>::StdNumber StdNumber;
    // This is for use by scalar col- and row- ordered matrices only. The same
    // code works because transpose(inv(m))==inv(transpose(m)).
    void lapackInvertInPlace() {
        // should have been checked already
        assert(this->m_eltSize==1 && this->nrow()==this->ncol()); 
        const int m = this->nrow();
        StdNumber* rawMem = reinterpret_cast<StdNumber*>(this->m_data);
        Array_<int> ipiv(m);
        int info;
        Lapack::getrf<StdNumber>(m,m,rawMem,this->m_leadingDim,&ipiv[0],info);
        assert(info==0);

        // Calculate optimal size for work
        StdNumber workSz;
        Lapack::getri<StdNumber>(m,rawMem,this->m_leadingDim,&ipiv[0],
                                 &workSz,-1,info);
        const int wsz = (int)CNT<StdNumber>::real(workSz);

        Array_<StdNumber> work(wsz);
        Lapack::getri<StdNumber>(m,rawMem,this->m_leadingDim,&ipiv[0],
                                 &work[0],wsz,info);
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
    FullColOrderEltHelper(int esz, int cppesz, int nr, int nc)
    :   Base(esz, cppesz, nr, nc, nr*esz) {   
        this->m_actual.setStorage
           (MatrixStorage(MatrixStorage::Full, MatrixStorage::NoPlacement,
                          MatrixStorage::ColumnOrder, MatrixStorage::NoDiag)); 
    }

    FullColOrderEltHelper(int esz, int cppesz, int nr, int nc, int ldim, 
                          S* shared, bool canWrite)
    :   Base(esz, cppesz, nr, nc, ldim, shared, canWrite) {
        assert(ldim>=nr*esz);
        this->m_actual.setStorage
           (MatrixStorage(MatrixStorage::Full, MatrixStorage::NoPlacement,
                          MatrixStorage::ColumnOrder, MatrixStorage::NoDiag));
    }

    // This should be processed a column at a time if possible.
    bool preferRowOrder_() const {return false;}

    // These virtual methods should be overridden in the derived scalar class.
    virtual const S*    getElt_(int i, int j) const 
    {   return this->m_data + this->eltIx(i,j); }
    virtual S*          updElt_(int i, int j)       
    {   return this->m_data + this->eltIx(i,j); }
    virtual bool        hasContiguousData_()  const 
    {   return this->isContiguousElt(this->nrow()); }
    virtual This*       cloneHelper_()        const 
    {   return new This(*this); }

    // This implementation will return a FullRowOrderEltHelper. Override for 
    // scalars.
    virtual RegularFullHelper<S>* createTransposeView_();

    // These implementations are fine for both composite and scalar class.

    // Regular spacing as (dfast/di, dfast/dj) and (dslow/di, dslow/dj), 
    // with i,j in elements.
    EltIndexer getEltIndexer() const {return EltIndexer(1,0,0,1);}

    This* createDeepCopy_() const {
        This* p = cloneHelper_();
        p->m_leadingDim = this->nrow()*this->m_eltSize; // packed now
        p->m_writable = true;
        p->m_owner = true;
        p->allocateData(this->nelt());
        if (hasContiguousData_())
            std::copy(this->m_data, this->m_data + 
                      Base::nelt()*this->m_eltSize, p->m_data);
        else for (int j=0; j < this->ncol(); ++j)
            std::copy(getElt_(0,j), getElt_(0,j) + 
                      this->nrow()*this->m_eltSize, p->updElt_(0,j));
        return p;
    }


    // OK for any size elements.
    void resize_(int m, int n) {
        this->clearData();
        this->allocateData(m,n);
        this->m_leadingDim = m * this->m_eltSize; // # scalars in a column
    }

    // OK for any size elements.
    void resizeKeep_(int m, int n) {
        const int newLeadingDim = m * this->m_eltSize; // # scalars in a column
        S* const newData = this->allocateMemory(m,n);
        const int colsToCopy = std::min(n, this->ncol());
        const int rowsToCopy = std::min(m, this->nrow()); // in elements
        for (int j=0; j < colsToCopy; ++j) {
            S*       const dest = newData + (ptrdiff_t)j*newLeadingDim;
            const S* const src  = this->m_data + (ptrdiff_t)j*this->m_leadingDim;
            std::copy(src, src+rowsToCopy*this->m_eltSize, dest);
        }
        this->clearData();
        this->setData(newData);
        this->m_leadingDim = newLeadingDim;
    }
};

// Full, column order, scalar, final.
template <class S>
class FullColOrderScalarHelper : public FullColOrderEltHelper<S> {
    typedef FullColOrderScalarHelper<S> This;
    typedef FullColOrderEltHelper<S>    Base;
public:
    // Leading dimension is # rows for contiguous column-ordered memory.
    FullColOrderScalarHelper(int nr, int nc)
    :   Base(1, 1, nr, nc) {}
    FullColOrderScalarHelper(int nr, int nc, int ldim, S* shared, bool canWrite)
    :   Base(1, 1, nr, nc, ldim, shared, canWrite) {}

    // For speed, these override the Base implementations for composite elements.
    const S*    getElt_(int i, int j) const 
    {   return this->m_data + this->scalarIx(i,j); }
    S*          updElt_(int i, int j)       
    {   return this->m_data + this->scalarIx(i,j); }
    bool        hasContiguousData_()  const 
    {   return this->isContiguousScalar(this->nrow()); }
    This*       cloneHelper_()        const 
    {   return new This(*this); }

    // This implementation will return a FullRowOrderScalarHelper.
    RegularFullHelper<S>* createTransposeView_();

    void colSum_(int j, S* csum) const {*csum = this->scalarColSum(j);}
    void rowSum_(int i, S* rsum) const {*rsum = this->scalarRowSum(i);}
    // Sum element column by column to avoid cache faults.
    void sum_(S* esum) const {
        *esum = S(0);
        for (int j=0; j<this->ncol(); ++j) *esum += this->scalarColSum(j);
    }

    void invertInPlace_() {this->lapackInvertInPlace();}
};

// Full, row order, composite. Column j is the fast dimension here.
template <class S>
class FullRowOrderEltHelper : public RegularFullHelper<S> {
    typedef FullRowOrderEltHelper<S>    This;
    typedef RegularFullHelper<S>        Base;
public:
    // Leading dimension is # cols for contiguous row-ordered memory.
    FullRowOrderEltHelper(int esz, int cppesz, int nr, int nc)
    :   Base(esz, cppesz, nr, nc, nc*esz) {   
        this->m_actual.setStorage
           (MatrixStorage(MatrixStorage::Full, MatrixStorage::NoPlacement,
                          MatrixStorage::RowOrder, MatrixStorage::NoDiag)); 
    }

    FullRowOrderEltHelper(int esz, int cppesz, int nr, int nc, int ldim, 
                          S* shared, bool canWrite)
    :   Base(esz, cppesz, nr, nc, ldim, shared, canWrite) {
        assert(ldim>=nc*esz);
        this->m_actual.setStorage
           (MatrixStorage(MatrixStorage::Full, MatrixStorage::NoPlacement,
                          MatrixStorage::RowOrder, MatrixStorage::NoDiag));
    }

    // This should be processed a row at a time if possible.
    bool preferRowOrder_() const {return true;}

    // These virtual methods should be overridden in the derived scalar class.
    virtual const S*    getElt_(int i, int j) const 
    {   return this->m_data + this->eltIx(j,i); }
    virtual S*          updElt_(int i, int j)       
    {   return this->m_data + this->eltIx(j,i); }
    virtual bool        hasContiguousData_()  const 
    {   return this->isContiguousElt(this->ncol()); }
    virtual This*       cloneHelper_()        const 
    {   return new This(*this); }

    // This implementation will return a FullColOrderEltHelper. Override for 
    // scalars.
    virtual RegularFullHelper<S>* createTransposeView_();

    // These implementations are fine for both composite and scalar class.
    EltIndexer getEltIndexer() const {return EltIndexer(0,1,1,0);}

    This* createDeepCopy_() const {
        This* p = cloneHelper_();
        p->m_leadingDim = this->ncol()*this->m_eltSize; // packed now
        p->m_writable = true;
        p->m_owner = true;
        p->allocateData(this->nelt());
        if (hasContiguousData_())
            std::copy(this->m_data, this->m_data + 
                      this->nelt()*this->m_eltSize, p->m_data);
        else for (int i=0; i < this->nrow(); ++i)
            std::copy(this->getElt_(i,0), this->getElt_(i,0) + 
                      this->ncol()*this->m_eltSize, p->updElt_(i,0));
        return p;
    }

    void resize_(int m, int n) {
        this->clearData();
        this->allocateData(m,n);
        this->m_leadingDim = n * this->m_eltSize;   // # scalars in a row
    }

    // OK for any size elements.
    void resizeKeep_(int m, int n) {
        const int newLeadingDim = n * this->m_eltSize; // # scalars in a row
        S* const newData = this->allocateMemory(m,n);
        const int colsToCopy = std::min(n, this->ncol()); // in elements
        const int rowsToCopy = std::min(m, this->nrow());
        for (int i=0; i < rowsToCopy; ++i) {
            S*       const dest = newData + (ptrdiff_t)i*newLeadingDim;
            const S* const src  = this->m_data + (ptrdiff_t)i*this->m_leadingDim;
            std::copy(src, src+colsToCopy*this->m_eltSize, dest);
        }
        this->clearData();
        this->setData(newData);
        this->m_leadingDim = newLeadingDim;
    }
};

// Full, row order, scalar, final.
template <class S>
class FullRowOrderScalarHelper : public FullRowOrderEltHelper<S> {
    typedef FullRowOrderScalarHelper<S> This;
    typedef FullRowOrderEltHelper<S>    Base;
public:
    // Leading dimension is # cols for contiguous row-ordered memory.
    FullRowOrderScalarHelper(int nr, int nc)
    :   Base(1, 1, nr, nc) {}
    FullRowOrderScalarHelper(int nr, int nc, int ldim, S* shared, bool canWrite)
    :   Base(1, 1, nr, nc, ldim, shared, canWrite) {}

    // For speed, these override the Base implementations for composite elements.
    const S*    getElt_(int i, int j) const 
    {   return this->m_data + this->scalarIx(j,i); }
    S*          updElt_(int i, int j)       
    {   return this->m_data + this->scalarIx(j,i); }
    bool        hasContiguousData_()  const 
    {   return this->isContiguousScalar(this->ncol()); }
    This*       cloneHelper_()        const 
    {   return new This(*this); }

    // This implementation will return a FullColOrderScalarHelper.
    RegularFullHelper<S>* createTransposeView_();

    void colSum_(int j, S* csum) const {*csum = this->scalarColSum(j);}
    void rowSum_(int i, S* rsum) const {*rsum = this->scalarRowSum(i);}
    // Sum element row by row to avoid cache faults.
    void sum_(S* esum) const {
        *esum = S(0);
        for (int i=0; i<this->nrow(); ++i) *esum += this->scalarRowSum(i);
    }

    void invertInPlace_() {this->lapackInvertInPlace();}
};


// These definitions for inline methods had to wait until other classes
// were defined.

template <class S> inline RegularFullHelper<S>*
FullColOrderEltHelper<S>::createTransposeView_() {
    FullRowOrderEltHelper<S>* p = 
        new FullRowOrderEltHelper<S>(this->m_eltSize, this->m_cppEltSize, 
                                     this->ncol(), this->nrow(), 
                                     this->m_leadingDim, this->m_data, this->m_writable);
    return p;
}

template <class S> inline RegularFullHelper<S>*
FullColOrderScalarHelper<S>::createTransposeView_() {
    FullRowOrderScalarHelper<S>* p = 
        new FullRowOrderScalarHelper<S>(this->ncol(), this->nrow(), 
                                        this->m_leadingDim, this->m_data, this->m_writable);
    return p;
}

template <class S> inline RegularFullHelper<S>*
FullRowOrderEltHelper<S>::createTransposeView_() {
    FullColOrderEltHelper<S>* p = 
        new FullColOrderEltHelper<S>(this->m_eltSize, this->m_cppEltSize, 
                                     this->ncol(), this->nrow(), 
                                     this->m_leadingDim, this->m_data, this->m_writable);
    return p;
}

template <class S> inline RegularFullHelper<S>*
FullRowOrderScalarHelper<S>::createTransposeView_() {
    FullColOrderScalarHelper<S>* p = 
        new FullColOrderScalarHelper<S>(this->ncol(), this->nrow(), 
                                        this->m_leadingDim, this->m_data, this->m_writable);
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
    FullIndexedEltHelper(RegularFullHelper<S>& src, const EltBlock& block, const EltIndexer& ix)
    :   Base(src.getEltSize(), src.getCppEltSize(), 
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

    virtual const S*    getElt_(int i, int j) const {return this->m_data + ixEltIx(i,j);}
    virtual S*          updElt_(int i, int j)       {return this->m_data + ixEltIx(i,j);}
    virtual This*       cloneHelper_()        const {return new This(*this);}

    // This implementation works for this class and the scalar derived class.
    This* createTransposeView_() {
        This* p = cloneHelper_();
        p->m_indexer = m_indexer.transpose();
        p->m_actual.setActualSize(this->ncol(), this->nrow());
        return p;
    }

    // A deep copy of an indexed matrix eliminates the index, producing a
    // full, column- or row-ordered matrix with contiguous storage. 
    // We choose the destination storage order to match whichever way is
    // fastest to travel through the indexed matrix, since then we can
    // traverse both matrices in their fastest order.
    virtual RegularFullHelper<S>* createDeepCopy_() const {
        if (preferRowOrder_()) {
            FullRowOrderEltHelper<S>* p = 
                new FullRowOrderEltHelper<S>(this->m_eltSize,this->m_cppEltSize,
                                             this->nrow(),this->ncol());
            for (int i=0; i < this->nrow(); ++i) {
                S* dest = p->updElt_(i,0);   // start of a dense row
                for (int j=0; j < this->ncol(); ++j, dest += this->m_eltSize)
                    this->copyElt(dest, this->getElt_(i,j));
            }
            return p;
        } else {
            FullColOrderEltHelper<S>* p = 
                new FullColOrderEltHelper<S>(this->m_eltSize,this->m_cppEltSize,
                                             this->nrow(), this->ncol());
            for (int j=0; j < this->ncol(); ++j) {
                S* dest = p->updElt_(0,j);   // start of a dense column
                for (int i=0; i < this->nrow(); ++i, dest += this->m_eltSize)
                    this->copyElt(dest, this->getElt_(i,j));
            }
            return p;
        }
    }

protected:
    EltIndexer m_indexer;

    int row(int i, int j) const {return m_indexer.row(i,j);}
    int col(int i, int j) const {return m_indexer.col(i,j);}

    ptrdiff_t ixEltIx(int i, int j) const {return this->eltIx(row(i,j),col(i,j));}
};

// Full, indexed, scalar, final.
template <class S>
class FullIndexedScalarHelper : public FullIndexedEltHelper<S> {
    typedef FullIndexedScalarHelper<S>  This;
    typedef FullIndexedEltHelper<S>     Base;
public:
    // Construct a new view of any regularly-spaced full matrix, as long as that
    // matrix has scalar elements.
    FullIndexedScalarHelper(RegularFullHelper<S>& src, const EltBlock& block, 
                            const EltIndexer& ix)
    :   Base(src, block, ix)
    {
        SimTK_ASSERT_ALWAYS(src.getEltSize()==1, 
            "FullIndexedScalarHelper ctor: source must have scalar elements");
    }

    // For speed, these override the composite-element implementatiosn.
    const S* getElt_(int i, int j) const 
    {   return this->m_data + ixScalarIx(i,j); }
    S*       updElt_(int i, int j)       
    {   return this->m_data + ixScalarIx(i,j); }

    This* cloneHelper_() const {return new This(*this);}

    RegularFullHelper<S>* createDeepCopy_() const {
        if (this->preferRowOrder_()) {
            FullRowOrderScalarHelper<S>* p = 
                new FullRowOrderScalarHelper<S>(this->nrow(),this->ncol());
            for (int i=0; i < this->nrow(); ++i) {
                S* dest = p->updElt_(i,0);   // start of a dense row
                for (int j=0; j < this->ncol(); ++j, ++dest)
                    *dest = *this->getElt_(i,j);
            }
            return p;
        } else {
            FullColOrderScalarHelper<S>* p = 
                new FullColOrderScalarHelper<S>(this->nrow(),this->ncol());
            for (int j=0; j < this->ncol(); ++j) {
                S* dest = p->updElt_(0,j);   // start of a dense column
                for (int i=0; i < this->nrow(); ++i, ++dest)
                    *dest = *this->getElt_(i,j);
            }
            return p;
        }
    }

private:
    ptrdiff_t ixScalarIx(int i, int j) const 
    {   return this->scalarIx(this->row(i,j),this->col(i,j)); }
};

template <class S> inline RegularFullHelper<S>*
RegularFullHelper<S>::createRegularView_(const EltBlock& block, 
                                         const EltIndexer& ix) {
    RegularFullHelper<S>* p;
    if (this->m_eltSize==1) p = new FullIndexedScalarHelper<S>(*this, block, ix);
    else                    p = new FullIndexedEltHelper<S>(*this, block, ix);
    return p;
}




} // namespace SimTK   


#endif // SimTK_SimTKCOMMON_MATRIX_HELPER_REP_FULL_H_
