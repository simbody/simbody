#ifndef SimTK_SIMMATRIX_MATRIX_HELPER_H_
#define SimTK_SIMMATRIX_MATRIX_HELPER_H_

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

#include "SimTKcommon/Scalar.h"

#include <iostream>
#include <cassert>
#include <complex>
#include <cstddef>

namespace SimTK {


template <class S> class MatrixHelperRep;

/**
 * Class of enums used to communicate  various attributes of matrices that
 * affect which alogrithms are used for factoring, solving etc.
 */
class MatrixStructures {
    public:
    enum Structure {
        Uncommitted        = 0,
        Full               = 1,  // all elements are significant
        Diagonal           = 2,
        Symmetric          = 3,  // also means Hermitian
        Triangular         = 4,
        QuasiTriangular    = 5,
        Hessenberg         = 6,
        Permutation        = 7,
        TriDiagonal        = 8
    };
};

class MatrixShapes {
    public:
    enum Shape {
        Uncommitted        = 0,
        General            = 1,  // 2d rectangular matrix
        Square             = 2,  // restriction to nrow==ncol
        Vector             = 3,  // 1d column vector
        RowVector          = 4   // 1d row vector
    };
};

class MatrixSparseFormats {
    public:
    enum Sparsity {
        Uncommitted        = 0,
        Full               = 1,
        Banded             = 2 // needs upper & lower bandwidth
    };
};
class MatrixStorageFormats {
    public:
    enum Storage {
        Uncommitted        = 0,
        Conventional       = 1,
        Packed             = 2, // triangular or symmetric
        HouseholderProduct = 3, // orthogonal only
        PivotArray         = 4  // pivot only
    };
};
class MatrixConditions {
    public:
    enum Condition {
        Uncommitted        = 0,
        Unknown            = 1,
        Orthogonal         = 2,
        PositiveDefinite   = 3,
        WellConditioned    = 4, // implies full rank
        FullRank           = 5, // but might have bad conditioning
        Singular           = 6  // implies possible bad conditioning
    };
};

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
    class DiagonalView  { };
    
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

    // These are the constructors for making a diagonal view of a Matrix.
    MatrixHelper(const MatrixHelper&, const DiagonalView&); // a read only diagonal view
    MatrixHelper(MatrixHelper&, const DiagonalView&);       // a writable diagonal view

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

    void setMatrixStructure(MatrixStructures::Structure );
    MatrixStructures::Structure getMatrixStructure() const;
    void setMatrixShape(MatrixShapes::Shape );
    MatrixShapes::Shape getMatrixShape() const;
    void setMatrixSparsity(MatrixSparseFormats::Sparsity );
    MatrixSparseFormats::Sparsity getMatrixSparsity() const;
    void setMatrixStorage(MatrixStorageFormats::Storage );
    MatrixStorageFormats::Storage getMatrixStorage() const;
    void setMatrixCondition(MatrixConditions::Condition );
    MatrixConditions::Condition getMatrixCondition() const;

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
