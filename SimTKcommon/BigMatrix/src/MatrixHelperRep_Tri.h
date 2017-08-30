#ifndef SimTK_SimTKCOMMON_MATRIX_HELPER_REP_TRI_H_
#define SimTK_SimTKCOMMON_MATRIX_HELPER_REP_TRI_H_

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


//--------------------------------- TriHelper ----------------------------------

/// This abstract class represents a square matrix for which only half the elements 
/// need to be stored in memory, and that half is either the lower or upper triangle.
/// The remaining half must have elements whose values can be determined from
/// the stored ones. There are several cases:
///
///     - the missing elements are all zero (a triangular matrix)
///     - the missing (i,j)th element is the same as the stored (j,i)th element
///       (a positionally symmetric matrix)
///     - the missing (i,j)th element is the generalized hermitian transpose
///       of the stored (j,i)th element (a generalized hermitian matrix)
///     - the missing elements are the negatives of the above (skew-symmetric
///       or skew-hermitian matrix)
///
/// In addition, we allow for the possibility that the diagonal
/// elements are all the same, known value and don't need to be stored; we'll
/// call that a "known diagonal" matrix. The known diagonal capability permits 
/// two full triangle-shaped matrices to be packed in the space of a single 
/// full matrix provided that at least one of them is a known diagonal matrix.
/// Typically the known diagonal value is one (unit diagonal), but it can have
/// other values (e.g. a skew symmetric matrix has a zero diagonal).
///
/// Read-only references to any of the stored elements, and to the known diagonal
/// elements if any, can be obtained through the usual getElt_(i,j) virtual.
/// Writable references to the stored elements (but not a diagonal element of
/// a known diagonal matrix) can be obtained with the updElt_(i,j) virtual.
///
/// In any matrix where the (i,j)th element must be derived
/// from the (j,i)th via some non-identity transformation, the diagonal
/// elements are required to be invariant under that transformation. For a
/// skew-symmetric matrix the diagonals must be zero to make m(i,j)==-m(i,j).
/// For a hermitian matrix, the diagonals must be real so that m(i,j)==conj(m(i,j)).
/// For a skew-hermitian matrix the diagonals must be imaginary so that
/// m(i,j)==-conj(m(i,j)). In general the diagonals must be recursively
/// self-hermitian-transposes for composite elements in hermitian matrices.
///
/// Although the interesting part of the matrix is square, we allow that to be
/// the upper left square of a rectangular matrix, where all the additional
/// elements are zero (that's called "trapezoidal" although it might really
/// be a trapezoid because the triangle part may be flipped. In any case if 
/// the matrix has dimension m X n, the leading square has dimension min(m,n).

//------------------------------------------------------------------------------
template <class S>
class TriHelper : public MatrixHelperRep<S> {
    typedef TriHelper<S>            This;
    typedef MatrixHelperRep<S>      Base;
    typedef typename CNT<S>::TNeg   SNeg;
    typedef typename CNT<S>::THerm  SHerm;
    typedef TriHelper<SNeg>         ThisNeg;
    typedef TriHelper<SHerm>        ThisHerm;
public:
    TriHelper(int esz, int cppesz) : Base(esz,cppesz) {}

    // Override in derived classes with assumed-unit diagonals.
    virtual bool hasUnitDiagonal() const {return false;}

    // Just changing the return type here.
    virtual This* cloneHelper_() const = 0;

    // A deep copy of a Tri matrix will always return another
    // Tri matrix.
    virtual This* createDeepCopy_() const = 0;


    // TODO: this should allow at least views which are Triangular, meaning that
    // the upper left corner should be on the diagonal.
    MatrixHelperRep<S>* createBlockView_(const EltBlock& block) {
        SimTK_ERRCHK_ALWAYS(!"implemented", "TriHelper::createBlockView_()", "not implemented");
        return 0;
    }
    // TODO: the transpose has conjugate elements, which we would normally handle
    // by typecasting to CNT<S>::THerm. However, here we're being asked to make a
    // view without typecasting, so we have to do it with the elements we've got.
    // No problem: turn a "lower" into an "upper"; the missing elements (the 
    // conjugated ones) will have moved to their proper locations.
    MatrixHelperRep<S>* createTransposeView_() {
        SimTK_ERRCHK_ALWAYS(!"not implemented", "TriHelper::createTransposeView_()", "not implemented");
        return 0;
    }

    // Not sure if this should ever be supported.
    MatrixHelperRep<S>*  createRegularView_(const EltBlock&, const EltIndexer&) {
        SimTK_ERRCHK_ALWAYS(!"not implemented", "TriHelper::createRegularView_()", "not implemented");
        return 0;
    }

    VectorHelper<S>* createColumnView_(int,int,int){
        SimTK_ERRCHK_ALWAYS(!"not implemented", "TriHelper::createColumnView_()", "not implemented");
        return 0;
    }
    VectorHelper<S>* createRowView_(int,int,int){
        SimTK_ERRCHK_ALWAYS(!"not implemented", "TriHelper::createRowView_()", "not implemented");
        return 0;
    }
protected:
};


//----------------------------- TriInFullHelper --------------------------------
/// We are assuming non-packed storage here: the triangle matrix is stored
/// in a full matrix so the stored elements are indexed exactly as they would
/// be if the matrix were full. Two triangle matrices can be stored
/// in the available memory, one in the lower triangle and one in the upper,
/// provided that at least one has known diagonal elements that don't need to be
/// stored (that's still not considered "packed").
//------------------------------------------------------------------------------
template <class S>
class TriInFullHelper : public TriHelper<S> {
    typedef typename CNT<S>::StdNumber  StdNumber;
    typedef TriInFullHelper<S>          This;
    typedef TriHelper<S>                Base;
public:
    // If we're allocating the space, we'll just allocate a square even if the
    // matrix is trapezoidal.
    TriInFullHelper(int esz, int cppesz, int m, int n, 
         bool triangular, bool hermitian, bool skew, bool rowOrder) 
    :   Base(esz,cppesz), m_minmn(std::min(m,n)), m_leadingDim(m_minmn*esz), 
        m_triangular(triangular), m_hermitian(hermitian), 
        m_skew(skew), m_rowOrder(rowOrder), m_unstored(new S[NumUnstored*esz])
    {
        this->m_owner     = true;
        this->m_writable  = true;
        this->allocateData(m_minmn, m_minmn);
        this->m_actual.setActualSize(m, n); // the apparent size
    }

    // Use someone else's memory, which we assume to be the right size. 
    TriInFullHelper(int esz, int cppesz, int m, int n,
         bool triangular, bool hermitian, bool skew, bool rowOrder,
         int ldim, const S* shared, bool canWrite) 
    :   Base(esz,cppesz), m_minmn(std::min(m,n)), m_leadingDim(ldim), 
        m_triangular(triangular), m_hermitian(hermitian), 
        m_skew(skew), m_rowOrder(rowOrder), m_unstored(new S[NumUnstored*esz])
    {        
        assert(m_leadingDim >= m_minmn*this->m_eltSize);
        this->m_owner     = false;
        this->m_writable  = canWrite;
        this->setData(const_cast<S*>(shared));
        this->m_actual.setActualSize(m, n);
    }

    ~TriInFullHelper() {delete[] m_unstored;}

    bool preferRowOrder_() const {return m_rowOrder;}

    // In full storage, Tri matrix elements are regularly spaced.
    bool hasRegularData_() const {return true;}

    // Just changing the return type here: a cloned Tri-in-full-storage
    // handle will always be another Tri-in-full-storage handle.
    virtual This* cloneHelper_() const = 0;

    // In full storage, a Tri matrix is not stored contiguously.
    bool hasContiguousData_() const {return false;}


    // Compute the missing elements.
    void getAnyElt_(int i, int j, S* value) const { 
        if (i==j || this->eltIsStored_(i,j)) {
            this->copyElt(value, this->getElt_(i,j));
            return; 
        }
        // Missing elements are all zero for triangular matrices, or if
        // either dimension is outside the square part of any triangular-storage
        // matrix (i.e., the matrix is trapezoidal).
        if (m_triangular || i >= m_minmn || j >= m_minmn) {
            this->fillElt(value, StdNumber(0));
            return;
        }

        // Symmetric or hermitian. The returned (i,j)th element is some
        // function of the stored (j,i)th element.
        const S* e = this->getElt_(j,i);

        if (m_hermitian) {
            if (m_skew) this->copyAndNegConjugateElt(value, e);
            else this->copyAndConjugateElt(value, e);
        } else {
            // symmetric but not hermitian
            if (m_skew) this->copyAndNegateElt(value, e);
            else this->copyElt(value, e); // no change to stored element
        }
    }

    VectorHelper<S>* createDiagonalView_();

    virtual bool isUpper() const {return false;}
    virtual bool hasKnownDiagonal() const {return false;}

    This* createDeepCopy_() const {
        This* p = cloneHelper_();
        p->m_writable = true;
        p->m_owner = true;
        p->allocateData(m_minmn, m_minmn);
        p->m_minmn = m_minmn;
        p->m_leadingDim = m_minmn; // might be different than source

        const bool shortFirst = (m_rowOrder&&!isUpper()) || (!m_rowOrder&&isUpper());
        const int  known = hasKnownDiagonal() ? 1 : 0;
        const int  eltSize = this->m_eltSize;

        // Copy 1/2 of a square of this size, possibly excluding diagonals.
        const int nToCopy = m_minmn;
        S*       const dest = p->m_data    + known*this->m_eltSize; // skip diag if known
        const S* const src  = this->m_data + known*this->m_eltSize;
        if (shortFirst) {
            // column or row begins at (k,j) and has length (j+1)-k
            int lengthElt = eltSize; // 1st row/col has 1 element to copy
            // skip 0-length row or col if diag is known
            for (ptrdiff_t j=known; j < nToCopy; ++j, lengthElt += eltSize) {
                std::copy(src+j*this->m_leadingDim, src+j*this->m_leadingDim +
                          lengthElt, dest+j*p->m_leadingDim);
            }
        } else { // longFirst
            // column or row begins at (j+k,j), length m-j-k
            int startInScalars = 0; // data was already shifted by 1 if needed
            int lengthElt = (nToCopy-known)*eltSize;
            for (ptrdiff_t j=0; j < nToCopy-known; ++j) {
                std::copy(src+j*m_leadingDim + startInScalars,
                          src+j*m_leadingDim + startInScalars + lengthElt,
                          dest+j*p->m_leadingDim + startInScalars);
                startInScalars += this->m_eltSize;
                lengthElt -= eltSize;
            }
        }
        return p;
    }

    // OK for any size elements. Note that we only allocate a square to hold
    // the data even if the matrix is trapezoidal. The rest is known to be zero.
    // There is no need to reallocate if min(m,n) doesn't change -- that just
    // means more zeroes so the size will be changed.
    void resize_(int m, int n) {
        const int newMinmn = std::min(m,n);
        if (newMinmn == m_minmn) return;
        this->clearData();
        m_minmn = newMinmn;
        this->allocateData(m_minmn, m_minmn);
        m_leadingDim = m_minmn*this->m_eltSize; // number of scalars in a column
    }

    // OK for any size elements. We'll move along the fast direction; columns
    // for column-ordered and rows for row-ordered storage. There are two cases:
    //   (upper, col order) and (lower, row order) 
    //      -- data starts at 0 for each column [row] and get longer (shortFirst)
    //   (upper, row order) and (lower, col order) 
    //      -- data starts at the diagonal and gets shorter (longFirst)
    void resizeKeep_(int m, int n) {
        const int newMinmn = std::min(m,n);
        if (newMinmn == m_minmn) return;

        const bool shortFirst = (m_rowOrder&&!isUpper()) || (!m_rowOrder&&isUpper());
        const int  newLeadingDim = newMinmn;
        const int  known = hasKnownDiagonal() ? 1 : 0;
        const int  eltSize = this->m_eltSize;

        // Copy 1/2 of a square of this size, possibly excluding diagonals.
        S* const newData = this->allocateMemory(newMinmn, newMinmn);
        const int nToCopy = std::min(m_minmn, newMinmn);
        S*       const dest = newData + known*this->m_eltSize; // skip diag if known
        const S* const src  = this->m_data  + known*this->m_eltSize;
        if (shortFirst) {
            // column or row begins at (k,j) and has length (j+1)-k
            int lengthElt = eltSize; // 1st row/col has 1 element to copy
            // skip 0-length row or col if diag is known
            for (ptrdiff_t j=known; j < nToCopy; ++j, lengthElt += eltSize) {
                std::copy(src+j*this->m_leadingDim, src+j*this->m_leadingDim +
                          lengthElt, dest+j*newLeadingDim);
            }
        } else { // longFirst
            // column or row begins at (j+k,j), length m-j-k
            int startInScalars = 0; // data was already shifted by 1 if needed
            int lengthElt = (nToCopy-known)*eltSize;
            for (ptrdiff_t j=0; j < nToCopy-known; ++j) {
                std::copy(src+j*this->m_leadingDim + startInScalars,
                          src+j*this->m_leadingDim + startInScalars + lengthElt,
                          dest+j*newLeadingDim + startInScalars);
                startInScalars += this->m_eltSize;
                lengthElt -= eltSize;
            }
        }
        this->clearData();
        this->setData(newData);
        m_minmn      = newMinmn;
        m_leadingDim = newLeadingDim;
    }

protected:
    int     m_minmn;        // dimension of the square part (min(m,n))
    int     m_leadingDim;   // in scalars; can be row- or column-ordered.

    // These should be implicit in the subclasses but for now they're here 
    bool    m_triangular;   // true if triangular, false if symmetric or hermitian
    bool    m_hermitian;    // m(i,j) = hermitianTranspose(m(j,i))
    bool    m_skew;         // m(i,j) = -op(m(j,i))
    bool    m_rowOrder;

    // Allocate heap space here for unstored elements. Don't forget to multiply
    // by element size when accessing.
    static const int NumUnstored   = 3; // length is 3*element size in scalars
    static const int UnstoredZero  = 0; // first element holds a zero
    static const int UnstoredDummy = 1; // second element is "write only" garbage
    static const int UnstoredDiag  = 2; // this is the known diagonal if any
    S*      m_unstored;

private:

};

// Concrete class: Tri matrix in the upper triangle of full storage,
// diagonals are stored, composite elements.
// An upper triangle has i <= j.
template <class S>
class TriInFullUpperHelper : public TriInFullHelper<S> {
    typedef TriInFullUpperHelper<S> This;
    typedef TriInFullHelper<S>      Base;
public:
    TriInFullUpperHelper(int esz, int cppesz, int m, int n, 
         bool triangular, bool hermitian, bool skew, bool rowOrder) 
    :   Base(esz,cppesz,m,n,triangular,hermitian,skew,rowOrder) 
    {   this->m_actual.setStructure(calcUpperTriStructure()); }

    // Use someone else's memory, which we assume to be the right size. 
    TriInFullUpperHelper(int esz, int cppesz, int m, int n,
         bool triangular, bool hermitian, bool skew, bool rowOrder,
         int ldim, const S* shared, bool canWrite) 
    :   Base(esz,cppesz,m,n,triangular,hermitian,skew,rowOrder,ldim,shared,canWrite)
    {   this->m_actual.setStructure(calcUpperTriStructure()); }

    bool isUpper() const {return true;}

    virtual This* cloneHelper_() const {return new This(*this);}

    // override for unit diagonal
    virtual bool eltIsStored_(int i, int j) const {return i <= j && j < this->m_minmn;}

    // override for unit diagonals and scalar elements
    virtual const S* getElt_(int i, int j) const {
        SimTK_ERRCHK2(i <= j && j < this->m_minmn, "TriInFullUpperHelper::getElt_()",
            "Element index was (i,j)=(%d,%d) which refers to an element which is\n"
            " not available. Only the upper triangle and diagonal of this matrix are stored.",
            i, j);
        if (this->m_rowOrder) std::swap(i,j);
        return this->m_data + j*this->m_leadingDim + i*this->m_eltSize;
    }

    // override for unit diagonals and scalar elements
    virtual S* updElt_(int i, int j) {
        SimTK_ERRCHK2(i <= j && j < this->m_minmn, "TriInFullUpperHelper::updElt_()",
            "Element index was (i,j)=(%d,%d) which refers to an element which is\n"
            " not stored. Only the upper triangle and diagonal of this matrix are stored.",
            i, j);
        if (this->m_rowOrder) std::swap(i,j);
        return this->m_data + j*this->m_leadingDim + i*this->m_eltSize;
    }

private:
    MatrixStructure calcUpperTriStructure() const {
        MatrixStructure ms;
        ms.setPosition(MatrixStructure::Upper);
        ms.setDiagValue(this->hasKnownDiagonal() ? MatrixStructure::UnitDiag : MatrixStructure::StoredDiag);
        if (this->m_triangular) ms.setStructure(MatrixStructure::Triangular);
        else {
            if (this->m_hermitian) 
                ms.setStructure(this->m_skew ? MatrixStructure::SkewHermitian 
                                             : MatrixStructure::Hermitian);
            else // symmetric
                ms.setStructure(this->m_skew ? MatrixStructure::SkewSymmetric 
                                             : MatrixStructure::Symmetric);
        }
        return ms;
    }
};


// Concrete class: Tri matrix in the upper triangle of full storage,
// diagonals are known, composite elements.
template <class S>
class TriInFullUpperKnownDiagHelper : public TriInFullUpperHelper<S> {
    typedef TriInFullUpperKnownDiagHelper<S>    This;
    typedef TriInFullUpperHelper<S>             Base;
public:
    TriInFullUpperKnownDiagHelper(int esz, int cppesz, int m, int n, 
         bool triangular, bool hermitian, bool skew, bool rowOrder, const S* knownDiag) 
    :   Base(esz,cppesz,m,n,triangular,hermitian,skew,rowOrder) 
    {
        copyElt(&this->m_unknown[this->UnstoredDiag*this->m_eltSize], knownDiag);
    }

    TriInFullUpperKnownDiagHelper(int esz, int cppesz, int m, int n,
         bool triangular, bool hermitian, bool skew, bool rowOrder, const S* knownDiag,
         int ldim, const S* shared, bool canWrite) 
    :   Base(esz,cppesz,m,n,triangular,hermitian,skew,rowOrder,ldim,shared,canWrite) 
    {
        copyElt(&this->m_unknown[this->UnstoredDiag*this->m_eltSize], knownDiag);
    }

    bool hasKnownDiagonal() const {return true;}

    This* cloneHelper_() const {return new This(*this);}

    // We're overriding since only i<j is stored.
    bool eltIsStored_(int i, int j) const {return i<j && j < this->m_minmn;}

    // These should be overwritten for scalars, although they will work as is.
    virtual const S* getElt_(int i, int j) const {
        SimTK_ERRCHK2(i<=j && j < this->m_minmn, "TriInFullUpperKnownDiagHelper::getElt_()",
            "Element index (i,j)=(%d,%d) which refers to an element which is\n"
            " not available. Only the upper triangle of this matrix is stored.",
            i, j);
        if (i==j) return &this->m_unknown[this->UnstoredDiag*this->m_eltSize];
        if (this->m_rowOrder) std::swap(i,j);
        return this->m_data + j*this->m_leadingDim + i*this->m_eltSize;
    }

    // override for unit diagonals and scalar elements
    virtual S* updElt_(int i, int j) {
        SimTK_ERRCHK2(i<j && j < this->m_minmn, "TriInFullUpperKnownDiagHelper::updElt_()",
            "Element index (i,j)=(%d,%d) which refers to the lower triangle, but\n"
            " only the upper triangle (i<j) of this matrix are stored.", i, j);
        if (this->m_rowOrder) std::swap(i,j);
        return this->m_data + j*this->m_leadingDim + i*this->m_eltSize;
    }
};





} // namespace SimTK   


#endif // SimTK_SimTKCOMMON_MATRIX_HELPER_REP_TRI_H_
