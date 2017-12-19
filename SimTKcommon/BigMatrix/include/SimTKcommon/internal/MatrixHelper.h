#ifndef SimTK_SIMMATRIX_MATRIX_HELPER_H_
#define SimTK_SIMMATRIX_MATRIX_HELPER_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-12 Stanford University and the Authors.        *
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

/** @file
 *
 * Here we declare the detailed interface to the Simmatrix classes. This is
 * still client-side code although nothing here should need to appear as part 
 * of the user documentation. Code here primarily deals with the interface
 * between the arbitrarily-templatized user visible Matrix, Vector, and
 * RowVector classes and the finitely templatized implementation classes,
 * which are available only for SimTK scalar types.
 */

#include "SimTKcommon/Scalar.h"

#include <iostream>
#include <cassert>
#include <complex>
#include <cstddef>
#include <utility> // for std::pair

namespace SimTK {


template <class S> class MatrixHelperRep;
class MatrixCharacter;
class MatrixCommitment;

//  ------------------------------ MatrixHelper --------------------------------

/// Here we define class MatrixHelper<S>, the scalar-type templatized helper
/// class for the more general, composite numerical type-templatized 
/// class MatrixBase<ELT>. The helper class is not intended to appear
/// directly in user programs; it is client-side code used by the client-side
/// Matrix<CNT>, Vector<CNT>, etc. templates to reduce the infinite set
/// of possible CNT templates to a finite set of scalar templates which
/// can then have hidden implementations in the Simmatrix library.
/// The hidden implementation class will be instantiated once each for 
/// float, double and the associated complex and conjugate types,
/// and their negators (a total of 12 types). Element size is
/// dealt with at run time; otherwise the helper knows nothing about the 
/// structure of the elements.
/// 
/// The constructors for numerical types should not initialize the numerical
/// values. We take advantage of that here -- this class assumes it can simply
/// allocate the appropriate amount of data as an array of the underlying
/// scalar type, with no implicit initialization. We'll fill uninitialized data
/// with NaNs when debugging or if specifically requested; otherwise it will
/// contain garbage initially.
/// 
/// Note that this is just a templatized handle class. The implementation
/// is private, in the undefined class MatrixHelperRep<S>.

//  ----------------------------------------------------------------------------
template <class S> 
class SimTK_SimTKCOMMON_EXPORT MatrixHelper {
    typedef MatrixHelper<S>                         This;
    typedef MatrixHelper<typename CNT<S>::TNeg>     ThisNeg;
    typedef MatrixHelper<typename CNT<S>::THerm>    ThisHerm;
public: 
    typedef typename CNT<S>::Number     Number;     // strips off negator from S
    typedef typename CNT<S>::StdNumber  StdNumber;  // turns conjugate to complex
    typedef typename CNT<S>::Precision  Precision;  // strips off complex from StdNumber

    // no default constructor
    // copy constructor suppressed

    // Destructor eliminates MatrixHelperRep object if "this" owns it.
    ~MatrixHelper() {deleteRepIfOwner();}

    // Local types for directing us to the right constructor at compile time.
    class ShallowCopy   { };
    class DeepCopy      { };
    class TransposeView { };
    class DiagonalView  { };
    
        //////////////////////////
        // "Owner" constructors //
        //////////////////////////
    
    // 0x0, fully resizable, fully uncommitted.
    MatrixHelper(int esz, int cppEsz);

    // Default allocation for the given commitment.
    MatrixHelper(int esz, int cppEsz, const MatrixCommitment&);

    // Allocate a matrix that satisfies a given commitment and has a
    // particular initial size.
    MatrixHelper(int esz, int cppEsz, const MatrixCommitment&, int m, int n);

    // Copy constructor that produces a new owner whose logical shape and contents are
    // the same as the source, but with a possibly better storage layout. Data will
    // be contiguous in the copy regardless of how spread out it was in the source.
    // The second argument is just to disambiguate this constructor from similar ones.
    MatrixHelper(const MatrixCommitment&, const MatrixHelper& source, const DeepCopy&);

    // This has the same semantics as the previous DeepCopy constructor, except that
    // the source has negated elements with respect to S. The resulting logical shape
    // and logical values are identical to the source, meaning that the negation must
    // actually be performed here, using one flop for every meaningful source scalar.
    MatrixHelper(const MatrixCommitment&, const MatrixHelper<typename CNT<S>::TNeg>& source, const DeepCopy&);
   
    
        //////////////////////////////////
        // External "View" constructors //
        //////////////////////////////////

    // These constructors produce a matrix which provides a view of externally-allocated
    // data, which is known to have the given MatrixCharacter. There is also provision
    // for a "spacing" parameter which defines gaps in the supplied data, although
    // the interpretation of that parameter depends on the MatrixCharacter (typically
    // it is the leading dimension for a matrix or the stride for a vector). Note
    // that spacing is interpreted as "number of scalars between adjacent elements"
    // which has the usual Lapack interpretation if the elements are scalars but
    // can be used for C++ vs. Simmatrix packing differences for composite elements.
    // The resulting handle is *not* the owner of the data, so destruction of the 
    // handle does not delete the data.

    // Create a read-only view into existing data.
    MatrixHelper(int esz, int cppEsz, const MatrixCommitment&,
                 const MatrixCharacter&, int spacing, const S* data);
    // Create a writable view into existing data.
    MatrixHelper(int esz, int cppEsz, const MatrixCommitment&,
                 const MatrixCharacter&, int spacing, S* data);


        ////////////////////////////////
        // Matrix "View" constructors //
        ////////////////////////////////

    // Matrix view constructors, read only and writable. Use these for submatrices,
    // rows, and columns. Indices are by *element* and these may consist of multiple
    // scalars of type template parameter S.

    // "Block" views
    MatrixHelper(const MatrixCommitment&, const MatrixHelper&, int i, int j, int nrow, int ncol);
    MatrixHelper(const MatrixCommitment&, MatrixHelper&,       int i, int j, int nrow, int ncol); 

    // "Transpose" views (note that this is Hermitian transpose; i.e., element
    // type is conjugated).
    MatrixHelper(const MatrixCommitment&, const MatrixHelper<typename CNT<S>::THerm>&, 
                 const TransposeView&);    // a read only transposed view
    MatrixHelper(const MatrixCommitment&, MatrixHelper<typename CNT<S>::THerm>&,       
                 const TransposeView&);    // a writable transposed view

    // "Diagonal" views.
    MatrixHelper(const MatrixCommitment&, const MatrixHelper&, const DiagonalView&); // a read only diagonal view
    MatrixHelper(const MatrixCommitment&, MatrixHelper&,       const DiagonalView&); // a writable diagonal view

    // "Indexed" view of a 1d matrix.
    MatrixHelper(const MatrixCommitment&, const MatrixHelper&, 
                 int n, const int* indices);
    MatrixHelper(const MatrixCommitment&, MatrixHelper&, 
                 int n, const int* indices);

    // These invoke the previous constructors but with friendlier index source.
    MatrixHelper(const MatrixCommitment& mc, const MatrixHelper& h, 
                 const Array_<int>& indices)
    {   new (this) MatrixHelper(mc, h, (int)indices.size(), indices.cbegin()); }
    MatrixHelper(const MatrixCommitment& mc, MatrixHelper& h, 
                 const Array_<int>& indices)
    {   new (this) MatrixHelper(mc, h, (int)indices.size(), indices.cbegin()); }

    // "Copy" an existing MatrixHelper by making a new view into the same data. 
    // The const form loses writability, non-const retains same writability status 
    // as source. If the source is already a view then the destination will have
    // an identical element filter so that the logical shape and values are the
    // same for both source and copy. (The second argument is used just to 
    // disambiguate this constructor from similar ones.)
    MatrixHelper(const MatrixCommitment&, const MatrixHelper& source, const ShallowCopy&);
    MatrixHelper(const MatrixCommitment&, MatrixHelper&       source, const ShallowCopy&);

        /////////////////////
        // Copy assignment //
        /////////////////////

    // Behavior of copy assignment depends on whether "this" is an owner or view. If
    // it's an owner it is resized and ends up a new, dense copy of "source" just as
    // for the DeepCopy constructor above. If "this" is a writable view, sizes must match 
    // and we copy elements from "source" to corresponding elements of "this". If
    // "this" is not writable then the operation will fail.
    MatrixHelper& copyAssign(const MatrixHelper& source);

    // Same as copyAssign() except the source has element types which are logically
    // negated from the element types of "this". Since the result must have the
    // same value as the source, this requires that all the copied elements are
    // actually negated at a cost of one flop per scalar.
    MatrixHelper& negatedCopyAssign(const MatrixHelper<typename CNT<S>::TNeg>&);

        /////////////////////
        // View assignment //
        /////////////////////

    // View assign always disconnects this helper from whatever view & data
    // it used to have (destructing as appropriate) and then makes it a new view
    // of the source. Writability is lost if the source is const, otherwise 
    // writability is inherited from the source.
    MatrixHelper& readOnlyViewAssign(const MatrixHelper& source);
    MatrixHelper& writableViewAssign(MatrixHelper&       source);

    // Restore helper to its just-constructed state. We forget everything except
    // the element size and handle commitment, which can never change. The
    // allocated helper will be the same as if we had just default-constructed
    // using the current commitment.
    void clear();

    // Return true if there is currently no data associated with this handle.
    bool isClear() const;

    // Using *element* indices, obtain a pointer to the beginning of a particular
    // element. This is always a slow operation compared to raw array access;
    // use sparingly. These will fail if the indices are outside the stored
    // portion of the matrix. Use getAnyElt() if you want something that always
    // works.
    const S* getElt(int i, int j) const;
    S*       updElt(int i, int j);

    // Faster for 1-d matrices (vectors) if you know you have one.
    const S* getElt(int i) const;
    S*       updElt(int i);

    // This returns the correct value for any element within the logical dimensions
    // of the matrix. In some cases it has to compute the value; in all cases
    // it has to copy it.
    void getAnyElt(int i, int j, S* value) const;

    // Faster for 1-d matrices (vectors) if you know you have one.
    void getAnyElt(int i, S* value) const;

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

    // Fill all our stored data with copies of the same supplied element.
    void fillWith(const S* eltp);

    // We're copying data from a C++ row-oriented matrix into our general
    // Matrix. In addition to the row ordering, C++ may use different spacing
    // for elements than Simmatrix does.
    void copyInByRowsFromCpp(const S* elts);
    
        // Scalar operations //

    // Fill every element with repeated copies of a single scalar value.
    void fillWithScalar(const StdNumber&);
            
    // Scalar multiply (and divide). This is useful no matter what the
    // element structure and will produce the correct result.
    void scaleBy(const StdNumber&);

    // This is only allowed for a matrix of real or complex or neg of those,
    // which is square, well-conditioned, and for which we have no view,
    // and element size 1.
    void invertInPlace();

    void dump(const char* msg=0) const; // For debugging -- comes out in a way you can feed to Matlab

        // Bookkeeping //
    
    // This is the number of logical *elements* in each column of this matrix; i.e., m.
    int    nrow() const;
    // This is the number of *elements* in each row of this matrix; i.e., n.
    int    ncol() const;
    // This is the total number of *elements* in the logical shape of this matrix, i.e. m*n.
    ptrdiff_t nelt() const; // nrow*ncol  
    // This is the number of elements if this is a 1d matrix (vector or rowvector). That is,
    // it is the same as one of nrow() or ncol(); the other must be 1. It is also the
    // same as nelt() but limited to fit in 32 bits.
    int length() const;

    // Change the logical size of the underlying memory for this Matrix to m X n, forgetting
    // everything that used to be there. This will fail if it would have to resize any
    // non-resizable dimension. However, it will succeed even on non-resizable matrices and
    // views provided the existing dimensions are already correct. If no size change is made,
    // no action is taken and the original data is still accessible.
    void resize(int m, int n);

    // Same as resize() except as much of the original data as will fit remains in place at
    // the same (i,j) location it had before. This may require copying the elements after
    // allocating new space. As for resize(), this is allowed for any Matrix whose dimensions
    // are already right, even if that Matrix is not resizable.
    void resizeKeep(int m, int n);

    void lockShape();
    void unlockShape();

    const MatrixCommitment& getCharacterCommitment() const;
    const MatrixCharacter&  getMatrixCharacter() const;
    void commitTo(const MatrixCommitment&);

    // Access to raw data. For now this is only allowed if there is no view
    // and the raw data is contiguous.
    bool      hasContiguousData() const;
    ptrdiff_t getContiguousDataLength() const;
    const S*  getContiguousData() const;
    S*        updContiguousData();

    void replaceContiguousData(S* newData, ptrdiff_t length, bool takeOwnership);
    void replaceContiguousData(const S* newData, ptrdiff_t length);
    void swapOwnedContiguousData(S* newData, ptrdiff_t length, S*& oldData);

    const MatrixHelperRep<S>& getRep() const {assert(rep); return *rep;}
    MatrixHelperRep<S>&       updRep()       {assert(rep); return *rep;}
    void setRep(MatrixHelperRep<S>* hrep)    {assert(!rep); rep = hrep;}
    MatrixHelperRep<S>* stealRep() 
    {   assert(rep); MatrixHelperRep<S>* stolen=rep; rep=0; return stolen; }

    void deleteRepIfOwner();
    void replaceRep(MatrixHelperRep<S>*);

    // Rep-stealing constructor. We're taking ownership of the supplied rep.
    // Internal use only!
    explicit MatrixHelper(MatrixHelperRep<S>*);

private:
    // Pointer to the private implementation object. This is the only
    // allowable data member in this class.
    class MatrixHelperRep<S>* rep;

    // Suppress copy constructor.
    MatrixHelper(const MatrixHelper&);

    // ============================= Unimplemented =============================
    // See comment in MatrixBase::matmul for an explanation.
    template <class SA, class SB>
    void matmul(const StdNumber& beta,   // applied to 'this'
                const StdNumber& alpha, const MatrixHelper<SA>& A, const MatrixHelper<SB>& B);    
    
friend class MatrixHelper<typename CNT<S>::TNeg>;
friend class MatrixHelper<typename CNT<S>::THerm>;
};
     
} //namespace SimTK

#endif // SimTK_SIMMATRIX_MATRIX_HELPER_H_
