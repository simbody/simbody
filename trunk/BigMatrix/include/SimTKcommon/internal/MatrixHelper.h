#ifndef SimTK_SIMMATRIX_MATRIX_HELPER_H_
#define SimTK_SIMMATRIX_MATRIX_HELPER_H_

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

/** @file
 *
 * Here we define class MatrixHelper<S>, the scalar-type templatized helper
 * class for the more general, composite numerical type-templatized 
 * class MatrixBase<ELT>. The helper class is not intended to appear
 * directly in user programs; it is client-side code used by the client-side
 * Matrix<CNT>, Vector<CNT>, etc. templates to reduce the infinite set
 * of possible CNT templates to a finite set of scalar templates which
 * can then have hidden implementations in the Simmatrix library.
 * The hidden implementation class will be instantiated once each for 
 * float, double, long double and the associated complex and conjugate types,
 * and their negators (a total of 18 types). Element size is
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
 * Note that this is just a templatized handle class. The implementation
 * is private, in the undefined class MatrixHelperRep<S>.
 */

#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/internal/Scalar.h"

#include <iostream>
#include <cassert>
#include <complex>
#include <cstddef>

namespace SimTK {

// Not all combinations of structure/sparsity/storage/condition
// are allowed. 

class MatrixShape {
public:
    enum Shape {
        Uncommitted,
        General,        // 2d rectangular matrix
        Square,         // restriction to nrow==ncol
        Vector,         // 1d column vector
        RowVector       // 1d row vector
    };

    MatrixShape() : shape(Uncommitted) { }
    // implicit conversion
    MatrixShape(Shape s) : shape(s) { }
    Shape shape;
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
    enum Structure {
        Uncommitted,
        Full,       // all elements are significant
        Diagonal,
        Symmetric,  // also means Hermitian
        Triangular,
        QuasiTriangular,
        Hessenberg,
        Permutation
    };

    MatrixStructure() : structure(Uncommitted) { }
    
    // implicit conversion
    MatrixStructure(Structure ms) : structure(ms) { }

    Structure structure;
};

class MatrixSparsity {
public:
    enum Sparsity {
        Uncommitted,
        Full,
        Banded // needs upper & lower bandwidth
    };

    // If Banded, how stuck are we on a particular bandwidth?
    enum BandwidthFreedom {
        Free        = 0x00,    // upper & lower are free
        FixedUpper  = 0x01,
        FixedLower  = 0x02,
        Fixed       = FixedUpper | FixedLower
    };

    MatrixSparsity() 
      : sparsity(Uncommitted), lowerBandwidth(-1), upperBandwidth(-1)
    {
    }

    MatrixSparsity(int lower, int upper) 
      : sparsity(Banded), lowerBandwidth(lower), upperBandwidth(upper)
    {
        assert(lower >= 0 && upper >= 0);
    }

    Sparsity sparsity;
    int lowerBandwidth;
    int upperBandwidth;
};

class MatrixStorage {
public:
    enum Storage {
        Uncommitted,
        Conventional,
        Packed,             // triangular or symmetric
        HouseholderProduct, // orthogonal only
        PivotArray          // pivot only
    };

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
      : storage(Uncommitted), position(Lower), assumptions(None)
    {
    }
    
    // also serves as implicit conversion from Storage type
    MatrixStorage(Storage s, Position p=Lower, Assumptions a=None) 
        : storage(s), position(p), assumptions(a)
    {
    }

    Storage     storage;
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
    enum Condition {
        Uncommitted,
        Unknown,
        Orthogonal,
        PositiveDefinite,
        WellConditioned,    // implies full rank
        FullRank,           // but might have bad conditioning
        Singular            // implies possible bad conditioning
    };

    MatrixCondition() : condition(Uncommitted) { }
    // implicit conversion
    MatrixCondition(Condition c) : condition(c) { }

    Condition condition;
};

template <class S> class MatrixHelperRep;

template <class S> 
class SimTK_SimTKCOMMON_EXPORT MatrixHelper {
    typedef typename CNT<S>::Number     Number;     // strips off negator from S
    typedef typename CNT<S>::StdNumber  StdNumber;  // turns conjugate to complex
    typedef typename CNT<S>::Precision  Precision;  // strips off complex from StdNumber
public:     
    // no default constructor, copy constructor suppressed

    // Destructor eliminates MatrixHelperRep object if "this" owns it.
    ~MatrixHelper();

    // Local types for directing us to the right constructor at compile time.
    class ShallowCopy   { };
    class DeepCopy      { };
    class TransposeView { };
    
    // Matrix owner constructors
    
    // 0x0, fully resizable
    explicit MatrixHelper(int esz, int cppEsz);

    // Restore helper to its just-constructed state. We forget everything except
    // the element size, which can never change.
    void clear();
    
    // (m*esz) x n, resizable with optional restrictions 
    MatrixHelper(int esz, int cppEsz, int m, int n, bool lockNrow=false, bool lockNcol=false);
    
    // This is a non-resizeable owner (of the data descriptor) but shares pre-existing raw data.
    MatrixHelper(int esz, int cppEsz, int m, int n, int leadingDim, const S*); // read only
    MatrixHelper(int esz, int cppEsz, int m, int n, int leadingDim, S*);       // writable
   
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

    // Construct from negated form, performing the required multiplications.
    MatrixHelper(const MatrixHelper<typename CNT<S>::TNeg>&, const DeepCopy&);
   
    // Behavior of copy assignment depends on whether this is an owner or view. If
    // it's an owner it is resized and ends up a new, dense copy of rhs. If
    // it's a view, sizes must match and we copy elements from rhs to corresponding
    // elements of lhs.  
    MatrixHelper& copyAssign(const MatrixHelper&);
    MatrixHelper& negatedCopyAssign(const MatrixHelper<typename CNT<S>::TNeg>&);

    // View assign always disconnects this helper from whatever view & data
    // it used to have (destructing as appropriate) and then makes it a new view
    // of the source. Writability is lost if the source is const, otherwise 
    // writability is inherited from the source.
    MatrixHelper& readOnlyViewAssign(const MatrixHelper&);
    MatrixHelper& writableViewAssign(MatrixHelper&);

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

    void dump(const char* msg=0) const; // For debugging -- comes out in a way you can feed to Matlab

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

    // Access to raw data. For now this is only allowed if there is no view
    // and the raw data is contiguous.
    bool hasContiguousData() const;
    long getContiguousDataLength() const;
    const S* getContiguousData() const;
    S* updContiguousData();
    void replaceContiguousData(S* newData, long length, bool takeOwnership);
    void replaceContiguousData(const S* newData, long length);
    void swapOwnedContiguousData(S* newData, int length, S*& oldData);

    // Pointer to the private implementation object.
    class MatrixHelperRep<S>* rep;

private:
    // Suppress copy constructor.
    MatrixHelper(const MatrixHelper&);
};
     
} //namespace SimTK

#endif // SimTK_SIMMATRIX_MATRIX_HELPER_H_
