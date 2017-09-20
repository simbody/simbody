#ifndef SimTK_SimTKCOMMON_MATRIX_HELPER_REP_VECTOR_H_
#define SimTK_SimTKCOMMON_MATRIX_HELPER_REP_VECTOR_H_

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

namespace SimTK {


/*-------------------------------- VectorHelper --------------------------------

This abstract class represents a 1d matrix, a.k.a. a Vector or a RowVector
(covector). However, its use is not limited to the SimTK::Vector and and
SimTK::RowVector classes; any skinny Matrix or skinny slice of a fatter
Matrix might use a VectorHelper since it provides faster access to memory
than a generic 2d helper would.
    
Most operations don't care whether this is a column or a row,
however we have to know so that we can support matrix operations on the
vector when necessary. For example, getElt(i,j) still has to work (as long
as the appropriate one of i or j is 0), even though getElt(i) is more
efficient and direction-agnostic.

The most common layout is that all elements are stored, either consecutively
in memory or with a regular stride. Vectors constructed from a larger pool 
of stored data via indexing are also important and need to be implemented 
efficiently.

TODO:  However, there are several other important
layouts that arise most commonly from row and column selections performed
on non-full matrices, like triangular or symmetric matrices. Supporting such
selections allows simple (if inefficient) implementations of operations on
mixed types of matrices by breaking them into row and column operations.
Vectors selected in this way can be repetitions of the same element (often
zero), negations or conjugations of stored elements, and may sometimes have
a single "distinguished" element whose value is known (this occurs for example
when crossing a non-stored unit diagonal).

TODO: Finally, we allow a Vector to be formed of a composition of smaller 
Vectors. Then the whole vector can be accessed by element or more efficiently
by segments.
------------------------------------------------------------------------------*/
template <class S>
class VectorHelper : public MatrixHelperRep<S> {
    typedef VectorHelper<S>         This;
    typedef MatrixHelperRep<S>      Base;
    typedef typename CNT<S>::TNeg   SNeg;
    typedef typename CNT<S>::THerm  SHerm;
    typedef VectorHelper<SNeg>      ThisNeg;
    typedef VectorHelper<SHerm>     ThisHerm;
public:
    VectorHelper(int esz, int cppesz, bool isRow) 
    :   Base(esz,cppesz),m_row(isRow) 
    {
    }


    // A deep copy of a Vector will always return another Vector, so we'll 
    // change the return type here.
    virtual This* createDeepCopy_() const = 0;

    // Just changing the return type here.
    virtual This* cloneHelper_() const = 0;

    bool preferRowOrder_() const {return m_row;}


    // Forward the two-index operations to the appropriate one-index operation.

    // One of the indices must be zero.
    bool eltIsStored_(int i, int j)           const {return eltIsStored_(i+j);}
    const S* getElt_ (int i, int j)           const {return getElt_(i+j);}
    S*       updElt_ (int i, int j)                 {return updElt_(i+j);}
    void getAnyElt_  (int i, int j, S* value) const {getAnyElt_(i+j, value);}


    // These are mandatory for vectors.
    virtual bool eltIsStored_(int i)           const = 0;
    virtual const S* getElt_ (int i)           const = 0;
    virtual S*       updElt_ (int i)                 = 0;
    virtual void getAnyElt_  (int i, S* value) const = 0;

    // (Positional) transpose view is identical to this one except that we'll 
    // call it a row rather than a column or vice versa.
    This* createTransposeView_() {
        This* p = cloneHelper_();
        p->m_data = this->m_data;
        p->m_row = !m_row;
        p->m_actual.updStorage().setOrder
           (p->m_row ? MatrixStorage::RowOrder : MatrixStorage::ColumnOrder);
        return p;
    }

    // Not sure if this should ever be supported.
    MatrixHelperRep<S>*  createRegularView_(const EltBlock&, const EltIndexer&) {
        SimTK_ERRCHK_ALWAYS(!"not implemented", "VectorHelper::createRegularView_()", "not implemented");
        return 0;
    }


protected:
    bool    m_row; // true if this should be a row when treated as a matrix
};

//----------------------------- FullVectorHelper -------------------------------
/// All elements of the Vector are stored. The simplest form of this has the
/// data pointing to the Vector's 0th element with all the rest following
/// consecutively in memory. Derived classes add stride or indexing and 
/// optimize for scalar elements.
//------------------------------------------------------------------------------
template <class S>
class FullVectorHelper : public VectorHelper<S> {
    typedef FullVectorHelper<S>  This;
    typedef VectorHelper<S>      Base;
public:
    FullVectorHelper(int esz, int cppesz, bool isRow) 
    :   Base(esz,cppesz,isRow) {}

    // This will always produce a 1-element "contiguous" column vector.
    VectorHelper<S>* createDiagonalView_();

    // Source matches size and shape of this row or column.
    void copyInFromCompatibleSource_(const MatrixHelperRep<S>& source) {
        if (this->getEltSize() == 1) {
            // The elements are scalars, so we can copy them directly.
            if (this->nrow() == 1) // a row vector
                for (int j=0; j<this->ncol(); ++j)
                    *this->updElt_(j) = *source.getElt(0,j);
            else // a column vector
                for (int i=0; i<this->nrow(); ++i)
                    *this->updElt_(i) = *source.getElt(i,0);
        }
        else {  // Here the elements are not scalars.
            if (this->nrow() == 1) // a row vector
                for (int j=0; j<this->ncol(); ++j)
                    this->copyElt(this->updElt_(j), source.getElt(0,j));
            else // a column vector
                for (int i=0; i<this->nrow(); ++i)
                    this->copyElt(this->updElt_(i), source.getElt(i,0));
        }
    }
};

//-------------------------- ContiguousVectorHelper ----------------------------
/// All elements of the Vector are stored, and they are contiguous in memory.
/// This class handles general elements; a derived class handles the special
/// case of scalar elements by overriding some of the methods for speed.
//------------------------------------------------------------------------------
template <class S>
class ContiguousVectorHelper : public FullVectorHelper<S> {
    typedef ContiguousVectorHelper<S>   This;
    typedef FullVectorHelper<S>         Base;
public:
    // Allocate our own space.
    ContiguousVectorHelper(int esz, int cppesz, int n, bool isRow) 
    :   Base(esz,cppesz,isRow)
    {
        this->m_owner     = true;
        this->m_writable  = true;
        this->allocateData(n);
        this->m_actual.setStructure(MatrixStructure::Matrix1d);
        this->m_actual.setStorage(
            MatrixStorage(MatrixStorage::Vector, MatrixStorage::NoPlacement, 
                          isRow ? MatrixStorage::RowOrder : MatrixStorage::ColumnOrder, 
                          MatrixStorage::NoDiag));
        this->m_actual.setActualSize(isRow?1:n, isRow?n:1); // apparent size; sets Outline
    }

    // Use someone else's memory, which we assume to be the right size.
    // We take care of stride elsewhere.
    ContiguousVectorHelper(int esz, int cppesz, int n, bool isRow, 
                           const S* shared, bool canWrite) 
    :   Base(esz,cppesz,isRow)
    {        
        this->m_owner     = false;
        this->m_writable  = canWrite;
        this->setData(const_cast<S*>(shared));
        this->m_actual.setStructure(MatrixStructure::Matrix1d);
        this->m_actual.setStorage(
            MatrixStorage(MatrixStorage::Vector, MatrixStorage::NoPlacement, 
                          isRow ? MatrixStorage::RowOrder 
                                : MatrixStorage::ColumnOrder, 
                          MatrixStorage::NoDiag));
        // apparent size; sets Outline
        this->m_actual.setActualSize(isRow?1:n, isRow?n:1); 
    }

    virtual This* cloneHelper_() const {return new This(*this);}

    // A block view of a full, contiguous row/column is either 
    // (1) a smaller full, contiguous row/column, or (2) a zero-width
    // slice. In the latter case it may no longer be a row or column so we
    // have to switch helper types to a Full mX0 or 0Xn matrix.
    MatrixHelperRep<S>* createBlockView_(const EltBlock& block) {
        const int m=block.nrow(), n=block.ncol();
        if (m && n) { // normal case; same orientation but smaller
            This* p = cloneHelper_();
            p->m_data = updElt_(block.row0() + block.col0());
            return p;
        }
        // Here we know at least one of m or n is zero.
        if (m==1 || n==1) { // i.e., (1,0) or (0,1)
            This* p = cloneHelper_(); // still 1d "contiguous"
            p->m_data = 0;
            p->m_row = (m==1); // regardless of what it was
            return p;
        }

        // Here one of m,n is zero and the other is > 1; not 1d any more.
        RegularFullHelper<S>* p = 0;
        if (this->getEltSize()==1)
             p = new FullColOrderScalarHelper<S>(m,n,m,(S*)0,this->m_writable);
        else p = new FullColOrderEltHelper<S>(this->getEltSize(), 
                                              this->getCppEltSize(),
                                              m,n,m,(S*)0,this->m_writable);
        // Called shared-data constructor so p is non-owner.
        return p;
    }

    // This creates an mx1 column vector.
    This* createColumnView_(int j, int i, int m) {
        This* p = cloneHelper_();
        p->m_data = m > 0 ? updElt_(i+j) : 0;
        p->m_row = false; // even if this was a row
        p->m_actual.updStorage().setOrder(MatrixStorage::ColumnOrder);
        return p;
    }

    // This creates a 1xn row vector.
    This* createRowView_(int i, int j, int n) {
        This* p = cloneHelper_();
        p->m_data = n > 0 ? updElt_(i+j) : 0;
        p->m_row = true; // even if this was a column
        p->m_actual.updStorage().setOrder(MatrixStorage::RowOrder);
        return p;
    }

    // Override for strided or indexed data.
    virtual bool hasContiguousData_() const {return true;}
    // Override for indexed data.
    virtual bool hasRegularData_() const {return true;}

    const S* getElt_ (int i)           const {return this->m_data + i*this->m_eltSize;}
    S*       updElt_ (int i)                 {return this->m_data + i*this->m_eltSize;}
    bool eltIsStored_(int i)           const {return true;}

    // Every element is stored so this just forwards to getElt(i).
    void getAnyElt_(int i, S* value) const 
    {   this->copyElt(value, getElt_(i)); }

    // OK for any size elements that are packed contiguously.
    This* createDeepCopy_() const {
        This* p = cloneHelper_();
        p->m_writable = true;
        p->m_owner = true;
        p->allocateData(this->nelt());
        std::copy(this->m_data, this->m_data + 
                  this->length()*this->m_eltSize, p->m_data);
        return p;
    }

    // One of the lengths must be 1.
    void resize_     (int m, int n)                 {resize_(m*n);}
    void resizeKeep_ (int m, int n)                 {resizeKeep_(m*n);}

    // OK for any size elements.
    void resize_(int n) {
        this->clearData();
        this->allocateData(n);
    }

    // OK for any size elements that are packed contiguously.
    void resizeKeep_(int n) {
        S* const newData = this->allocateMemory(n);
        const int nToCopy = std::min(n, this->length());
        std::copy(this->m_data, this->m_data + 
                  nToCopy*this->m_eltSize, newData);
        this->clearData();
        this->setData(newData);
    }

    void copyInFromCompatibleSource_(const MatrixHelperRep<S>& source) {
        if (source.hasContiguousData() && this->nScalars())
            std::copy(source.getElt(0,0), source.getElt(0,0) + 
                      this->nScalars(), this->m_data);
        else
            FullVectorHelper<S>::copyInFromCompatibleSource_(source);
    }
};


//----------------------- ContiguousVectorScalarHelper -------------------------
/// All elements of the Vector are stored, and they are contiguous in memory,
/// and the elements are scalars. This inherits most functionality from the 
/// contiguous general-element case.
//------------------------------------------------------------------------------
template <class S>
class ContiguousVectorScalarHelper : public ContiguousVectorHelper<S> {
    typedef ContiguousVectorScalarHelper<S>     This;
    typedef ContiguousVectorHelper<S>           Base;
public:
    // Allocate our own space.
    ContiguousVectorScalarHelper(int n, bool isRow) : Base(1,1,n,isRow) {}
    // Use someone else's memory.
    ContiguousVectorScalarHelper(int n, bool isRow, const S* shared, bool canWrite) 
    :   Base(1,1,n,isRow,shared,canWrite) {}

    This* cloneHelper_() const {return new This(*this);}

    const S* getElt_ (int i) const {return this->m_data + i;}
    S*       updElt_ (int i)       {return this->m_data + i;}

    // Every element is stored so this just forwards to getElt(i).
    void getAnyElt_(int i, S* value) const {*value = *getElt_(i);}
};


template <class S> inline VectorHelper<S>* 
FullVectorHelper<S>::createDiagonalView_() {
    VectorHelper<S>* p = 0;
    const int nDiags = std::min(this->length(), 1); // 0 or 1
    S* data = nDiags ? this->updElt_(0) : 0;

    p = (this->m_eltSize==1) 
        ? new ContiguousVectorScalarHelper<S>(nDiags, false, data, false)
        : new ContiguousVectorHelper<S>(this->m_eltSize, this->m_cppEltSize, nDiags,
                                        false, data, false);
    return p;
}


//---------------------------- StridedVectorHelper -----------------------------
/// This is a vector with regularly-spaced but non-contiguous elements. This
/// is only for views so does not include an implementation for resizing.
//------------------------------------------------------------------------------
template <class S>
class StridedVectorHelper : public FullVectorHelper<S> {
    typedef StridedVectorHelper<S>  This;
    typedef FullVectorHelper<S>     Base;
public:
    /// Use someone else's memory, which we assume to be the right size. Note
    /// that stride is given in elements, with stride==1 meaning the elements
    /// are packed contiguously, in which case this is the wrong helper class to use.
    StridedVectorHelper(int esz, int cppesz, int n, bool isRow, 
         int strideInElements, const S* shared, bool canWrite) 
    :   Base(esz,cppesz,isRow), m_spacing((ptrdiff_t)strideInElements * esz)
    {        
        SimTK_ASSERT1(strideInElements >= 2, 
            "StridedVectorHelper::ctor(): illegal stride %d", strideInElements);
        this->m_owner     = false;
        this->m_writable  = canWrite;
        this->setData(const_cast<S*>(shared));
        this->m_actual.setStructure(MatrixStructure::Matrix1d);
        this->m_actual.setStorage(
            MatrixStorage(MatrixStorage::Vector, MatrixStorage::NoPlacement, 
                          isRow ? MatrixStorage::RowOrder : MatrixStorage::ColumnOrder, 
                          MatrixStorage::NoDiag));
        this->m_actual.setActualSize(isRow?1:n, isRow?n:1); // apparent size; sets Outline
    }

    virtual This* cloneHelper_() const {return new This(*this);}

    // A block view of a strided vector is usually a smaller vector with 
    // identical stride. No stride needed if there are fewer than two 
    // elements, though. Also, this could be an mX0 or 0Xn slice which is
    // no longer a Vector unless m==1 or n==1.
     MatrixHelperRep<S>* createBlockView_(const EltBlock& block) {
        const int m=block.nrow(), n=block.ncol();
        if ((m==0 && n!=1) || (n==0 && m!=1)) {
            // One or both dimensions is 0 and the other is not 1, so this
            // is no longer a 1d object.
            RegularFullHelper<S>* p = 0;
            if (this->getEltSize()==1)
                 p = new FullColOrderScalarHelper<S>(m,n,m,(S*)0,this->m_writable);
            else p = new FullColOrderEltHelper<S>(this->getEltSize(), 
                                                  this->getCppEltSize(),
                                                  m,n,m,(S*)0,this->m_writable);
            // Called shared-data constructor so p is non-owner.
            return p;
        }

        // At least one of m,n is a 1. Could still be 1x0 or 0x1.

        VectorHelper<S>* p = 0;
        const int start = block.row0() + block.col0(); // one of those is zero
        const int length = m*n;  // one of those is one
        S* data = length ? updElt_(start) : 0;

        if (length <= 1) {
            p = (this->m_eltSize==1) 
                ? new ContiguousVectorScalarHelper<S>(length, false, data, false)
                : new ContiguousVectorHelper<S>(this->m_eltSize, this->m_cppEltSize, length,
                                                false, data, false);
            // called a shared-data constructor, so p is non-owner
            return p;
        }

        p = cloneHelper_();
        p->setData(data);
        return p;
    }

    // Row and column view are like block view but without the possibility
    // of a non-Vector result.
    VectorHelper<S>* createColumnView_(int j, int i, int m) {
        VectorHelper<S>* p = 0;
        const int start = i+j; // one of those is zero
        S* data = m ? updElt_(start) : 0;

        if (m <= 1) {
            p = (this->m_eltSize==1) 
                ? new ContiguousVectorScalarHelper<S>(m, false, data, false)
                : new ContiguousVectorHelper<S>(this->m_eltSize, this->m_cppEltSize, m, 
                                                false, data, false);
            return p;
        }

        // length is > 1 so this must already be a column
        assert(j==0);
        p = cloneHelper_();
        p->setData(data);
        return p;
    }

    VectorHelper<S>* createRowView_(int i, int j, int n) {
        VectorHelper<S>* p = 0;
        const int start = i+j; // one of those is zero
        S* data = n ? updElt_(start) : 0;

        if (n <= 1) {
            p = (this->m_eltSize==1) 
                ? new ContiguousVectorScalarHelper<S>(n, true, data, false)
                : new ContiguousVectorHelper<S>(this->m_eltSize, this->m_cppEltSize, n,
                                                true, data, false);
            return p;
        }

        // length is > 1 so this must already be a row
        assert(i==0);
        p = cloneHelper_();
        p->setData(data);
        return p;
    }


    bool hasContiguousData_() const {return false;}
    bool hasRegularData_()    const {return true;}

    bool eltIsStored_(int i)           const {return true;}
    const S* getElt_ (int i)           const {return this->m_data + i*m_spacing;}
    S*       updElt_ (int i)                 {return this->m_data + i*m_spacing;}

    // Every element is stored so this just forwards to getElt(i).
    void getAnyElt_(int i, S* value) const 
    {   this->copyElt(value, getElt_(i)); }

    /// A deep copy of a strided vector produces a contiguous (stride==1)
    /// vector containing the same number of elements.
    FullVectorHelper<S>* createDeepCopy_() const {
        ContiguousVectorHelper<S>* p = 
            new ContiguousVectorHelper<S>(this->m_eltSize, this->m_cppEltSize, 
                                          this->length(), this->m_row);
        for (int i=0; i < this->length(); ++i)
            this->copyElt(p->updData() + i*this->m_eltSize, this->getData() + i*m_spacing);
        return p;
    }

protected:
    ptrdiff_t m_spacing; // m_eltsize*stride (spacing in scalars)
};


//------------------------- StridedVectorScalarHelper --------------------------
/// This is a vector with regularly-spaced but non-contiguous elements, where
/// the elements are scalars. Most functionality is inherited from the general-
/// element case but we specialize a few methods here for speed. This
/// is only for views so does not include an implementation for resizing.
//------------------------------------------------------------------------------
template <class S>
class StridedVectorScalarHelper : public StridedVectorHelper<S> {
    typedef StridedVectorScalarHelper<S>    This;
    typedef StridedVectorHelper<S>          Base;
public:
    /// Use someone else's memory, which we assume to be the right size. Note
    /// that stride is given in elements, with stride==1 meaning the elements
    /// are packed contiguously, in which case this is the wrong helper class to use.
    StridedVectorScalarHelper(int n, bool isRow, int stride, const S* shared, bool canWrite) 
    :   Base(1,1,n,isRow,stride,shared,canWrite) {}

    This* cloneHelper_() const {return new This(*this);}

    // Every element is stored so this just forwards to getElt(i).
    void getAnyElt_(int i, S* value) const {*value = *this->getElt_(i);} 

    /// A deep copy of a strided vector produces a contiguous (stride==1)
    /// vector containing the same number of elements.
    FullVectorHelper<S>* createDeepCopy_() const {
        ContiguousVectorScalarHelper<S>* p = 
            new ContiguousVectorScalarHelper<S>(this->length(), this->m_row);
        const int nToCopy = this->length();
        for (int i=0; i < nToCopy; ++i)
            p->updData()[i] = this->getData()[i*this->m_spacing];
        return p;
    }
};

//--------------------------- IndexedVectorHelper ------------------------------
/// All elements of the Vector are stored, but they are selected irregularly
/// from the available data. This is only for views so does not include an 
/// implementation for resizing, and there is no constructor for data allocation.
/// This class handles general elements; a derived class handles the special
/// case of scalar elements by overriding some of the methods for speed.
/// TODO: we only allow 32 bit indices, although 64 bit indices
/// are possible (just need another helper class for "long long" indices).
//------------------------------------------------------------------------------
template <class S>
class IndexedVectorHelper : public FullVectorHelper<S> {
    typedef IndexedVectorHelper<S>  This;
    typedef FullVectorHelper<S>     Base;
public:
    // Use someone else's memory, which we assume to be the right size. Here the
    // indices must be in terms of elements. We'll store them internally in terms
    // of scalars instead. We insist here that the indices are all nonnegative and
    // in monotonically increasing order.
    IndexedVectorHelper(int esz, int cppesz, int n, bool isRow, 
         const int* eltIndices, const S* shared, bool canWrite) 
    :   Base(esz,cppesz,isRow), m_scalarIndices(0)
    {        
        this->m_owner     = false;
        this->m_writable  = canWrite;
        this->setData(const_cast<S*>(shared));
        this->m_actual.setStructure(MatrixStructure::Matrix1d);
        this->m_actual.setStorage(
            MatrixStorage(MatrixStorage::Vector, MatrixStorage::NoPlacement, 
                          isRow ? MatrixStorage::RowOrder : MatrixStorage::ColumnOrder, 
                          MatrixStorage::NoDiag));
        this->m_actual.setActualSize(isRow?1:n, isRow?n:1); // apparent size; sets Outline

        if (n) {
            m_scalarIndices = new int[n];
            for (int i=0; i<n; ++i) {
                SimTK_ERRCHK(i==0 || eltIndices[i]>eltIndices[i-1], "IndexedVectorHelper::ctor()",
                    "Indices must be in monotonically ascending order.");
                m_scalarIndices[i] = this->m_eltSize*eltIndices[i];
            }
        }
    }

    // Copy constructor must copy indices.
    IndexedVectorHelper(const This& src) : Base(src), m_scalarIndices(0) {
        if (src.length()) {
            m_scalarIndices = new int[src.length()];
            std::copy(src.m_scalarIndices, src.m_scalarIndices +
                      src.length(), m_scalarIndices);
        }
    }


    // (Virtual) destructor must delete indices.
    ~IndexedVectorHelper() {delete[] m_scalarIndices;}

    // Note that cloning the helper also copies the indices.
    virtual This* cloneHelper_() const {return new This(*this);}

    // A block view of an indexed vector is a shorter indexed vector.
    // Also, this could be an mX0 or 0Xn slice which is
    // no longer a Vector unless m==1 or n==1.
     MatrixHelperRep<S>* createBlockView_(const EltBlock& block) {
        const int m=block.nrow(), n=block.ncol();
        if ((m==0 && n!=1) || (n==0 && m!=1)) {
            // One or both dimensions is 0 and the other is not 1, so this
            // is no longer a 1d object.
            RegularFullHelper<S>* p = 0;
            if (this->getEltSize()==1)
                 p = new FullColOrderScalarHelper<S>(m,n,m,(S*)0,this->m_writable);
            else p = new FullColOrderEltHelper<S>(this->getEltSize(), 
                                                  this->getCppEltSize(),
                                                  m,n,m,(S*)0,this->m_writable);
            // Called a shared-data constructor so p is non-owner.
            return p;
        }

        // At least one of the dimensions is a 1. Could still be 1x0 or 0x1.

        const int start = block.row0() + block.col0(); // one of these is zero
        const int length = block.nrow()*block.ncol(); // one of these is 1

        if (length <= 1) {
            // No need for indices; this is contiguous now.
            S* data = length ? updElt_(start) : 0;
            ContiguousVectorHelper<S>* p = (this->m_eltSize==1) 
                ? new ContiguousVectorScalarHelper<S>
                        (length, false, data, false)
                : new ContiguousVectorHelper<S>
                        (this->m_eltSize, this->m_cppEltSize, length,
                         false, data, false);
            // Called a shared-data constructor so p is non-owner.
            return p;
        }

        // Still an indexed vector.
 
        This* p = new This(*this, true); // don't copy the indices
        p->m_data = this->m_data;
        p->m_scalarIndices = new int[length];
        std::copy(m_scalarIndices+start, m_scalarIndices+start + 
                  length, p->m_scalarIndices);
        return p;
    }

    VectorHelper<S>* createColumnView_(int j, int i, int m) {
        if (m <= 1) {
            S* data = m ? updElt_(i+j) : 0; // one of i or j is 0
            VectorHelper<S>* p = (this->m_eltSize==1) 
                ? (VectorHelper<S>*)new ContiguousVectorScalarHelper<S>(m, false, data, false)
                : (VectorHelper<S>*)new ContiguousVectorHelper<S>(this->m_eltSize, this->m_cppEltSize, 
                                                                  m, false, data, false);
            return p;
        }

        // This must already be a column or we couldn't make a >1 element column here.
        assert(j==0);
        This* p = new This(*this, true); // don't copy the indices
        p->setData(this->m_data); // leaving the indices the same, so data starts at 0
        p->m_scalarIndices = new int[m];
        std::copy(m_scalarIndices+i, m_scalarIndices+i+m, p->m_scalarIndices);
        return p;
    }


    VectorHelper<S>* createRowView_(int i, int j, int n) {
        if (n<= 1) {
            S* data = n ? updElt_(i+j) : 0; // one of i or j is 0
            VectorHelper<S>* p = (this->m_eltSize==1) 
                ? (VectorHelper<S>*)new ContiguousVectorScalarHelper<S>(n, true, data, false)
                : (VectorHelper<S>*)new ContiguousVectorHelper<S>(this->m_eltSize, this->m_cppEltSize, 
                                                                  n, true, data, false);
            return p;
        }

        // This must already be a row or we couldn't make a >1 element row here.
        assert(i==0);
        This* p = new This(*this, true); // don't copy the indices
        p->setData(this->m_data); // leaving the indices the same, so data starts at 0
        p->m_scalarIndices = new int[n];
        std::copy(m_scalarIndices+j, m_scalarIndices+j+n, p->m_scalarIndices);
        return p;
    }

    bool hasContiguousData_() const {return false;}
    bool hasRegularData_()    const {return false;}

    const S* getElt_ (int i)           const {return this->m_data + m_scalarIndices[i];}
    S*       updElt_ (int i)                 {return this->m_data + m_scalarIndices[i];}
    bool eltIsStored_(int i)           const {return true;}

    // Every element is stored so this just forwards to getElt(i).
    void getAnyElt_(int i, S* value) const 
    {   this->copyElt(value, getElt_(i)); }

    // A deep copy of an indexed vector produces a contiguous, non-indexed vector.
    ContiguousVectorHelper<S>* createDeepCopy_() const {
        ContiguousVectorHelper<S>* p = 
            new ContiguousVectorHelper<S>(this->m_eltSize, this->m_cppEltSize, 
                                          this->length(), this->m_row);
        for (int i=0; i<this->length(); ++i)
            this->copyElt(p->updElt_(i), getElt_(i));
        return p;
    }

protected:
    int*  m_scalarIndices;

private:
    // This is like a copy constructor, but it does not copy the indices. The second
    // parameter is a dummy.
    IndexedVectorHelper(const This& src, bool) : Base(src), m_scalarIndices(0) {}
};



} // namespace SimTK   


#endif // SimTK_SimTKCOMMON_MATRIX_HELPER_REP_VECTOR_H_
