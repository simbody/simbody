#ifndef SimTK_SimTKCOMMON_MATRIX_HELPER_REP_H_
#define SimTK_SimTKCOMMON_MATRIX_HELPER_REP_H_

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

#include "SimTKcommon/Scalar.h"
#include "SimTKcommon/SmallMatrix.h"

#include "ElementFilter.h"

#include <cstddef>

namespace SimTK {
template <class S> class MatrixHelper;      // the helper handle class
template <class S> class MatrixHelperRep;   // implementation base class

template <class S> class FullHelper;        // all elements stored
template <class S> class TriHelper;         // half the elements are stored
template <class S> class VectorHelper;      // 1-d array of elements

/* --------------------------- MatrixHelperRep ---------------------------------
 *
 * This is the private implementation of the MatrixHelper<S> class. Note that
 * this must be explicitly instantiated in library-side source code for all 
 * possible values of the template parameter S, which is just the set of all 
 * SimTK scalar types.
 *
 * Every MatrixHelperRep has the following members:
 *  - handle properties
 *  - matrix character commitment
 *  - actual matrix character
 *  - pointer to memory
 *  - pointer to memory's reference counter
 *
 * Concrete derivatives of MatrixHelperRep contain additional information 
 * used to map logical element indices to physical memory.
 *
 * There are some extremely performance-sensitive aspects to this object, 
 * especially with regard to selecting individual elements of a matrix. In order
 * to get here we have already had to follow the pointer in the Matrix handle.
 * We would like to be able to return an element in the fewest additional 
 * instructions possible. This is made difficult by the variety of run time 
 * features we want to support here:
 *  - writable vs. read-only access to data
 *  - row- and column-ordered data
 *  - scalar elements vs. composite ones
 *  - 1d Vector and RowVector objects vs. 2d Matrix objects
 *  - various ways of selecting meaningful elements from memory based
 *    on regular and irregular spacing
 *
 * If we had to check flags for each of these possibilities before finally
 * indexing the element of interest, access would be substantially slowed. We
 * are particularly concerned about performance for the most common case
 * of regularly-spaced scalar elements. Instead of checking flags at run time, 
 * we want to make use of the virtual function table so that the all of the 
 * above tests can be avoided with a single indirect call through the virtual 
 * function branch table. At run time the compiler will then translate a call 
 * like "rep->getElt(i,j)" to "call rep->vft[getElt](i,j)" where the particular
 * function called is the right one given this rep's particular combination of 
 * all the possibilities mentioned above, which is encapsulated in the concrete
 * type of the MatrixHelperRep object. Attempts to call update methods of 
 * non-writable reps will call directly to methods which throw an appropriate 
 * error condition.
 *
 * Sorting out all these issues at compile time requires a large number of
 * small classes. All classes are templatized by Scalar type <S> and explicitly 
 * instantiated for all the SimTK scalar types. 
 *
 * At handle construction, a suitable default helper is assigned with a minimal 
 * handle commitment derived from the handle type -- element size, column or 
 * row outline, and whether the handle must be a view. Everything else will be 
 * uncommitted. If a later assignment requires a different concrete helper, the
 * original one will be replaced but the commitments will remain unchanged. 
 * Note that MatrixHelpers cannot be shared; each is paired with just one 
 * handle.
 * ---------------------------------------------------------------------------*/
template <class S> 
class MatrixHelperRep {
    typedef MatrixHelperRep<S>                      This;
    typedef MatrixHelperRep<typename CNT<S>::TNeg>  ThisNeg;
    typedef MatrixHelperRep<typename CNT<S>::THerm> ThisHerm;
public:
    typedef typename CNT<S>::Number     Number;     // strips off negator from S
    typedef typename CNT<S>::StdNumber  StdNumber;  // turns conjugate to complex
    typedef typename CNT<S>::Precision  Precision;  // strips off complex from StdNumber
    typedef typename CNT<S>::TNeg       SNeg;       // the negated version of S
    typedef typename CNT<S>::THerm      SHerm;      // the conjugate version of S

    // Constructors are protected; for use by derived classes only.

    // This is the MatrixHelper factory method for owner matrices. Given a 
    // particular matrix character, this method will deliver the best 
    // MatrixHelperRep subclass with those characteristics.
    static This* createOwnerMatrixHelperRep(int esz, int cppEsz, 
                                            const MatrixCharacter&);

    // This is the MatrixHelper factory method for matrices providing a view
    // into externally-allocated data. The supplied MatrixCharacter describes
    // the actual properties of the external data as we want it to appear
    // through this view. Const and non-const methods are provided. There is
    // provision for a "spacing" argument for the external data; this may
    // mean different things. Typically it is the leading dimensions for
    // full-storage matrices, and the stride for vectors. We'll return the
    // best MatrixHelperRep we can come up with for these characteristics; it
    // will not be an owner rep.
    static This* createExternalMatrixHelperRep
       (int esz, int cppEsz, const MatrixCharacter& actual, 
        int spacing, S* data, bool canWrite=true);

    // This is the factory for const access to externally-allocated storage.
    static This* createExternalMatrixHelperRep
       (int esz, int cppEsz, const MatrixCharacter& actual, int spacing, 
        const S* data)
    {   return createExternalMatrixHelperRep(esz, cppEsz, actual, spacing,
                                             const_cast<S*>(data), false); }

    // Here we have determined that this source is acceptable for the data
    // already allocated in this handle. In particular both dimensions match.
    void copyInFromCompatibleSource(const MatrixHelperRep<S>& source) {
        if (!m_writable)
            SimTK_THROW1(Exception::OperationNotAllowedOnNonconstReadOnlyView, 
                         "assignment");
        copyInFromCompatibleSource_(source);
    }

    // Same as above except the source is negated as it is copied in.
    void negatedCopyInFromCompatibleSource(const ThisNeg& source) {
        copyInFromCompatibleSource
           (reinterpret_cast<const MatrixHelperRep<S>&>(source));
        scaleBy(StdNumber(-1));
    }

    // Fill every element with repeated copies of a single scalar value.
    // The default implementation is very slow.
    void fillWithScalar(const StdNumber& scalar) {
        if (!m_writable)
            SimTK_THROW1(Exception::OperationNotAllowedOnNonconstReadOnlyView, 
                         "fillWithScalar()");
        fillWithScalar_(scalar);
    }

    // These are used in copy constructors. Hence the resulting value must
    // be the same as the source value, even if we're removing a negator<> from
    // the elements. So the createNegatedDeepCopy actually has to negate every 
    // element with (yuck) floating point operations. Note that the returned 
    // Matrix is always an owner, writable, and has whatever is the most 
    // efficient storage type for the source data.
    This* createDeepCopy() const {return createDeepCopy_();}
    ThisNeg* createNegatedDeepCopy() const {
        ThisNeg* p = reinterpret_cast<const ThisNeg*>(this)->createDeepCopy_();
        p->scaleBy(StdNumber(-1));
        return p;
    }

    This* createWholeView(bool wantToWrite) const {
        This* p       = cloneHelper_();
        p->m_owner    = false;
        p->m_writable = m_writable && wantToWrite; 
        p->m_data     = m_data;
        p->m_actual   = m_actual;
        return p;
    }


    ThisHerm* createHermitianView(bool wantToWrite) const {
        const ThisHerm* thisHerm = reinterpret_cast<const ThisHerm*>(this);
        ThisHerm* p = const_cast<ThisHerm*>(thisHerm)->createTransposeView_();
        p->m_writable = m_writable && wantToWrite;
        p->m_actual.setActualSize(ncol(),nrow());
        return p;
    }


    This* createBlockView(const EltBlock& block, bool wantToWrite) const {
        SimTK_SIZECHECK(block.row0(), nrow()-block.nrow(), 
                        "MatrixHelperRep::createBlockView()");
        SimTK_SIZECHECK(block.col0(), ncol()-block.ncol(), 
                        "MatrixHelperRep::createBlockView()");
        This* p = const_cast<This*>(this)->createBlockView_(block);
        p->m_writable = m_writable && wantToWrite;
        p->m_actual.setActualSize(block.nrow(), block.ncol());
        return p;
    }
    This* createRegularView(const EltBlock& block, const EltIndexer& ix, 
                            bool wantToWrite) const {
        SimTK_SIZECHECK(block.row0(), nrow()-ix.row(block.nrow(), block.ncol()), 
                        "MatrixHelperRep::createRegularView()");
        SimTK_SIZECHECK(block.col0(), ncol()-ix.col(block.nrow(), block.ncol()), 
                        "MatrixHelperRep::createRegularView()");
        This* p = const_cast<This*>(this)->createRegularView_(block, ix);
        p->m_writable = m_writable && wantToWrite;
        p->m_actual.setActualSize(block.nrow(), block.ncol());
        return p;
    }

        // 1-d views //

    VectorHelper<S>* createDiagonalView(bool wantToWrite) const {
        VectorHelper<S>* p = const_cast<This*>(this)->createDiagonalView_();
        p->m_writable = m_writable && wantToWrite;
        p->m_actual.setActualSize(std::min(nrow(),ncol()), 1); // a column
        return p;
    }
    // Select column j or part of column j.
    VectorHelper<S>* createColumnView(int j, int i, int m, bool wantToWrite) const {
        SimTK_INDEXCHECK(j, ncol(), "MatrixHelperRep::createColumnView()");
        SimTK_SIZECHECK(i, nrow()-m, "MatrixHelperRep::createColumnView()");
        VectorHelper<S>* p = const_cast<This*>(this)->createColumnView_(j,i,m);
        p->m_writable = m_writable && wantToWrite;
        p->m_actual.setActualSize(m, 1); // a column
        return p;
    }
    // Select row i or part of row i.
    VectorHelper<S>* createRowView(int i, int j, int n, bool wantToWrite) const {
        SimTK_INDEXCHECK(i, nrow(), "MatrixHelperRep::createRowView()");
        SimTK_SIZECHECK(j, ncol()-n, "MatrixHelperRep::createRowView()");
        VectorHelper<S>* p = const_cast<This*>(this)->createRowView_(i,j,n);
        p->m_writable = m_writable && wantToWrite;
        p->m_actual.setActualSize(1, n); // a row
        return p;
    }

    // Is it more efficient to access this Matrix in row order rather than
    // the default column order? If so, you should try to operate on
    // all columns of each row before moving on to the next row. Otherwise,
    // try to operate on all rows of each column before moving on
    // to the next column.
    bool preferRowOrder() const {return preferRowOrder_();}

    // Is the memory that we ultimately reference organized contiguously?
    bool hasContiguousData() const {return hasContiguousData_();}

    // Using *element* indices, obtain a pointer to the beginning of a 
    // particular element. This is always a slow operation compared to raw 
    // array access; use sparingly.
    //
    // These inline base class interface routines provide Debug-mode range 
    // checking but evaporate completely in Release mode so that the underlying
    // virtual method is called directly.
    const S* getElt(int i, int j) const {
        SimTK_INDEXCHECK(i, nrow(), "MatrixHelperRep::getElt(i,j)");
        SimTK_INDEXCHECK(j, ncol(), "MatrixHelperRep::getElt(i,j)");
        return getElt_(i,j);
    }
    S* updElt(int i, int j) {
        SimTK_INDEXCHECK(i, nrow(), "MatrixHelperRep::updElt(i,j)");
        SimTK_INDEXCHECK(j, ncol(), "MatrixHelperRep::updElt(i,j)");
        SimTK_ERRCHK(m_writable, "MatrixHelperRep::updElt()", 
                     "Matrix not writable.");
        return updElt_(i,j);
    }

    void getAnyElt(int i, int j, S* value) const {
        SimTK_INDEXCHECK(i, nrow(), "MatrixHelperRep::getAnyElt(i,j)");
        SimTK_INDEXCHECK(j, ncol(), "MatrixHelperRep::getAnyElt(i,j)");
        SimTK_ERRCHK(value, "MatrixHelperRep::getAnyElt()", 
            "The value return pointer must not be null.");
        return getAnyElt_(i,j,value);
    }

    const S* getElt(int i) const {
        SimTK_INDEXCHECK(i, length(), "MatrixHelperRep::getElt(i)");
        return getElt_(i);
    }
    S* updElt(int i) {
        SimTK_INDEXCHECK(i, length(), "MatrixHelperRep::updElt(i)");
        SimTK_ERRCHK(m_writable, "MatrixHelperRep::updElt(i)", 
                     "Matrix not writable.");
        return updElt_(i);
    }

    void getAnyElt(int i, S* value) const {
        SimTK_INDEXCHECK(i, length(), "MatrixHelperRep::getAnyElt(i)");
        SimTK_ERRCHK(value, "MatrixHelperRep::getAnyElt()", 
            "The value return pointer must not be null.");
        return getAnyElt_(i,value);
    }

    // Many matrix storage formats assume values for some of their elements
    // and thus those elements are not stored. Generic algorithms must avoid
    // writing on such elements, so we expect to be able to query the concrete
    // class to find out if particular elements actually have a memory
    // representation.
    //
    // Note: for symmetric matrices only one of (i,j) and (j,i) should return
    // true here, and for unit diagonals (i,i) will return false.
    bool eltIsStored(int i, int j) const {
        SimTK_INDEXCHECK(i, nrow(), "MatrixHelperRep::eltIsStored(i,j)");
        SimTK_INDEXCHECK(j, ncol(), "MatrixHelperRep::eltIsStored(i,j)");
        return eltIsStored_(i,j);
    }
    bool eltIsStored(int i) const {
        SimTK_INDEXCHECK(i, length(), "MatrixHelperRep::eltIsStored(i)");
        return eltIsStored_(i);
    }

    // Add up all the *elements*, returning the answer in the element given
    // by pointer to its first scalar.
    void sum(S* eltp) const {sum_(eltp);}
    void colSum(int j, S* eltp) const {colSum_(j,eltp);}
    void rowSum(int i, S* eltp) const {rowSum_(i,eltp);}



    // Resizing is permitted only when one of these two conditions is met:
    //     1. The matrix already has the indicated dimensions,
    //  OR 2.  (a) this handle is the owner of the data,
    //     AND (b) handle commitment permits resizing of at least one dimension,
    //     AND (c) data has not been locked.
    void resize(int m, int n, bool keep) {
        SimTK_SIZECHECK_NONNEG(m, "MatrixHelperRep::resize()");
        SimTK_SIZECHECK_NONNEG(n, "MatrixHelperRep::resize()");
        if (m==nrow() && n==ncol()) return;
        if (!m_owner)
            SimTK_THROW1(Exception::OperationNotAllowedOnView, "resize()");
        // owner
        if (m_handleIsLocked)
            SimTK_THROW1(Exception::Cant, 
                "resize(), because owner handle is locked.");
        if (!m_commitment.isSizeOK(m,n))
            SimTK_THROW1(Exception::Cant, 
                "resize(), because handle commitment doesn't allow this size.");
        if (keep) resizeKeep_(m,n);
        else      resize_(m,n);
        m_actual.setActualSize(m,n);
    }


    // addition and subtraction (+= and -=)
    void addIn(const MatrixHelper<S>&);   
    void addIn(const MatrixHelper<typename CNT<S>::TNeg>&);   
    void subIn(const MatrixHelper<S>&); 
    void subIn(const MatrixHelper<typename CNT<S>::TNeg>&); 
    
    // Fill all our stored data with copies of the same supplied element.
    void fillWith(const S* eltp);

    // We're copying data from a C++ row-oriented matrix into our general
    // Matrix. In addition to the row ordering, C++ may use different spacing
    // for elements than Simmatrix does. Lucky we know that value!
    void copyInByRowsFromCpp(const S* elts);
            
    // Scalar multiply (and divide). This is useful no matter what the
    // element structure and will produce the correct result.
    void scaleBy(const StdNumber&);

    // If this is a full or symmetric square matrix with scalar elements,
    // it can be inverted in place (although that's not the best way to
    // solve linear equations!).
    void invertInPlace() {
        if (!m_writable)
            SimTK_THROW1(Exception::OperationNotAllowedOnNonconstReadOnlyView, 
                         "invertInPlace()");
        if (nrow() != ncol())
            SimTK_THROW1(Exception::Cant, 
                         "invertInPlace(): matrix must be square");
        invertInPlace_();
    }

    void dump(const char* msg) const;

    // See comment in MatrixBase::matmul for an explanation.
    template <class SA, class SB>
    void matmul(const StdNumber& beta,   // applied to 'this'
                const StdNumber& alpha, const MatrixHelper<SA>& A, const MatrixHelper<SB>& B);
    
        // Bookkeeping //


    // Dimensions come from the "actual" Matrix character. These are the *logical*
    // dimensions, in elements (not scalars), and have nothing to do with how much
    // memory is being used to store these elements.
    int       nrow()     const {return m_actual.nrow();}
    int       ncol()     const {return m_actual.ncol();}

    // This is only for 1D vectors; their lengths must fit in an "int". For
    // a 2D matrix use nelt() to get the total number of elements.
    int length() const {assert(nrow()==1||ncol()==1); return nrow()*ncol();}

    ptrdiff_t nelt()     const {return ptrdiff_t(nrow()) * ptrdiff_t(ncol());}
    ptrdiff_t nScalars() const {return nelt()*m_eltSize;}

    // Does this handle own the data descriptor it points to?
    int     getEltSize()        const {return m_eltSize;}
    int     getCppEltSize()     const {return m_cppEltSize;}
    bool    isWritable()        const {return m_writable;}
    bool    isOwner()           const {return m_owner;}
    bool    isClear()           const {return m_data==0;}

    void lockHandle()           {m_handleIsLocked=true;}
    void unlockHandle()         {m_handleIsLocked=false;}
    bool isHandleLocked() const {return m_handleIsLocked;}


    const MatrixCommitment& getCharacterCommitment() const {return m_commitment;}
    const MatrixCharacter&  getMatrixCharacter()  const {return m_actual;}

    // Get a pointer to the raw data viewed various ways.
    const S*     getData()     const {assert(m_data); return m_data;}
    S*           updData()           {assert(m_data); return m_data;}

    const SNeg*  getDataNeg()  const 
    {   return reinterpret_cast<const SNeg*>(getData()); }
    SNeg*        updDataNeg()        
    {   return reinterpret_cast<SNeg*>(updData()); }
    const SHerm* getDataHerm() const 
    {   return reinterpret_cast<const SHerm*>(getData()); }
    SHerm*       updDataHerm()       
    {   return reinterpret_cast<SHerm*>(updData()); }

    // Delete the data if necessary, and leave the data pointer null. If this
    // is not an owner handle, we just clear the pointer. If it is an owner
    // handle, we check that it isn't locked and if it isn't we'll return the
    // allocated data to the heap. It is an error to attempt to clear the
    // data while the handle is locked.
    void clearData() {
        if (m_handleIsLocked) {
            SimTK_THROW1(Exception::Cant, 
                "MatrixHelperRep::clearData(): handle is locked.");
            return;
        }
        if (m_owner)
            deleteAllocatedMemory(m_data);
        m_data = 0;
    }

    // Given a number of elements (not scalars), allocate
    // just enough memory to hold that many elements in packed storage.
    void allocateData(ptrdiff_t nelt) {
        assert(nelt >= 0);
        assert(m_owner && m_data==0);
        if (m_handleIsLocked) {
            SimTK_THROW1(Exception::Cant, 
                "MatrixHelperRep::allocateData(): handle is locked.");
            return;
        }
        m_data = allocateMemory(nelt);
    }

    // Given dimensions in number of elements (not scalars), allocate
    // just enough memory to hold m*n elements in packed storage.
    void allocateData(int m, int n) 
    {   assert(m>=0 && n>=0);
        allocateData(ptrdiff_t(m) * ptrdiff_t(n)); }

    // Allocate new heap space to hold nelt densely-packed elements each 
    // composed of m_eltSize scalars of type S. We actually allocate an array of
    // the underlying Precision type (float or double) to avoid any default 
    // construction of more complicated elements like complex. If we're in Debug
    // mode, we'll initialize the resulting data to NaN, otherwise we won't 
    // touch it. If nelt is zero we return a null pointer.
    S* allocateMemory(ptrdiff_t nElt) const {
        assert(nElt >= 0);
        if (nElt==0) 
            return 0;
        assert(sizeof(S) % sizeof(Precision) == 0);
        const ptrdiff_t nPrecPerElt = (sizeof(S)/sizeof(Precision))*m_eltSize;
        const ptrdiff_t nPrec       = nElt * nPrecPerElt;
        Precision* p = new Precision[nPrec];
        #ifndef NDEBUG
            const Precision nan = CNT<Precision>::getNaN();
            for (ptrdiff_t i=0; i < nPrec; ++i)
                p[i] = nan;
        #endif
        return reinterpret_cast<S*>(p);
    }
    // Allocate enough memory to hold m*n elements. m and n are ints, but their
    // product may not fit in an int, so we use ptrdiff_t which will be a 64-bit
    // signed integer on a 64 bit machine.
    S* allocateMemory(int m, int n) const 
    {   assert(m>=0 && n>=0);
        return allocateMemory(ptrdiff_t(m) * ptrdiff_t(n)); }

    // Use this method to delete help space that you allocated using 
    // allocateNewData() above. We recast it back to the form in which it was 
    // allocated before deleting which will keep the heap system happy and also 
    // prevent the calling of any element destructor that might be present for 
    // the fancier scalar types.
    static void deleteAllocatedMemory(S* mem) {
        Precision* p = reinterpret_cast<Precision*>(mem);
        delete[] p;
    }

    // Use setData only when there isn't already data in this handle. If this
    // is an owner handle we're taking over responsibility for the heap space.
    void setData(S* datap) {assert(!m_data); m_data = datap;}
    void setOwner(bool isOwner) {m_owner=isOwner;}

    // Single element manipulation: VERY SLOW, use sparingly.    
    void copyElt(S* dest, const S* src) const
    {   for (int k=0; k<m_eltSize; ++k) dest[k] = src[k]; }
    void copyAndScaleElt(S* dest, const StdNumber& s, const S* src) const
    {   for (int k=0; k<m_eltSize; ++k) dest[k] = s*src[k]; }
    void copyAndNegateElt(S* dest, const S* src) const
    {   for (int k=0; k<m_eltSize; ++k) dest[k] = CNT<S>::negate(src[k]); }
    void copyAndConjugateElt(S* dest, const S* src) const
    {   for (int k=0; k<m_eltSize; ++k) dest[k] = CNT<S>::transpose(src[k]); }
    void copyAndNegConjugateElt(S* dest, const S* src) const
    {   for (int k=0; k<m_eltSize; ++k) dest[k] = -CNT<S>::transpose(src[k]); }

    void fillElt(S* dest, const StdNumber& src) const
    {   for (int k=0; k<m_eltSize; ++k) dest[k] = src; }   
    void addToElt(S* dest, const S* src) const
    {   for (int k=0; k<m_eltSize; ++k) dest[k] += src[k]; }        
    void subFromElt(S* dest, const S* src) const
    {   for (int k=0; k<m_eltSize; ++k) dest[k] -= src[k]; }
    void scaleElt(S* dest, const StdNumber& s) const
    {   for (int k=0; k<m_eltSize; ++k) dest[k] *= s; }

protected:
    //--------------------------------------------------------------------------
    //                        ABSTRACT INTERFACE
    // This is the complete set of virtual methods which can be overridden
    // by classes derived from MatrixHelperRep. All the virtual methods have
    // names ending in an underscore "_". The base class provides an inline
    // interface method of the same name but without the underscore. That 
    // method performs operations common to all the implementations, such
    // as checking arguments and verifying that the handle permits the
    // operation. It then transfers control to the virtual method, the 
    // various implementations of which do not have to repeat the common
    // work.
    //--------------------------------------------------------------------------

        // DESTRUCTOR
        // Any concrete class that uses additional heap space must provide
        // a destructor.

    virtual ~MatrixHelperRep();

        // PURE VIRTUALS
        // All concrete classes must provide implementations.

    virtual MatrixHelperRep* createDeepCopy_() const = 0;

    // Make a clone of the current helper, except that the data and
    // myHandle pointer are left null.
    virtual This* cloneHelper_() const = 0;


    virtual This*   createBlockView_(const EltBlock&) = 0;
    virtual This*   createTransposeView_() = 0;
    virtual This*   createRegularView_(const EltBlock&, const EltIndexer&) = 0;

    virtual VectorHelper<S>* createDiagonalView_() = 0;
    virtual VectorHelper<S>* createColumnView_(int j, int i, int m) = 0;
    virtual VectorHelper<S>* createRowView_(int i, int j, int n) = 0;

    virtual bool preferRowOrder_() const = 0;
    virtual bool hasContiguousData_() const = 0;
    virtual bool hasRegularData_() const = 0;
    virtual bool eltIsStored_(int i, int j) const = 0;
    virtual const S* getElt_(int i, int j) const = 0;
    virtual S*       updElt_(int i, int j)       = 0;
    virtual void getAnyElt_(int i, int j, S* value) const = 0;

        // OPTIONAL FUNCTIONALITY
        // Optional methods. No default implementation. A derived class can 
        // supply these, but if it doesn't then this functionality will not 
        // be available for matrices which use that class as a helper.

    virtual void resize_(int m, int n) {
        SimTK_THROW1(Exception::Cant, 
            "resize_() not implemented for this kind of matrix");
    }
    virtual void resizeKeep_(int m, int n) {
        SimTK_THROW1(Exception::Cant, 
            "resizeKeep_() not implemented for this kind of matrix");
    }

    // If this gets called we know the matrix is writable and square.
    virtual void invertInPlace_()  {
        SimTK_THROW1(Exception::Cant, 
            "invertInPlace_() not implemented for this kind of matrix");
    }

        // One-index versions of above two-index methods for use in Vector
        // helpers.
    virtual bool eltIsStored_(int i) const {
        SimTK_THROW1(Exception::Cant, 
            "One-index eltIsStored_() not available for 2D matrices");
        return false;
    }
    virtual const S* getElt_(int i) const {
        SimTK_THROW1(Exception::Cant, 
            "One-index getElt_() not available for 2D matrices");
        return 0;
    }

    virtual S* updElt_(int i) {
        SimTK_THROW1(Exception::Cant, 
            "One-index updElt_() not available for 2D matrices");
        return 0;
    }

    virtual void getAnyElt_(int i, S* value) const {
        SimTK_THROW1(Exception::Cant, 
            "One-index getAnyElt_() not available for 2D matrices");
    }

    virtual void resize_(int n) {
        SimTK_THROW1(Exception::Cant, 
            "One-index resize_() not available for 2D matrices");
    }

    virtual void resizeKeep_(int n) {
        SimTK_THROW1(Exception::Cant, 
            "One-index resizeKeep_() not available for 2D matrices");
    }



        // VIRTUALS WITH DEFAULT IMPLEMENTATIONS
        // This functionality is required of all MatrixHelpers, but there is
        // a base class default implementation here, slow but functional.
        // In many cases a concrete class can do much better because of its 
        // intimate knowledge of the data layout; override if you can.


    // Overridable method to implement copyInFromCompatibleSource().
    // The default implementation works but is very slow.
    virtual void copyInFromCompatibleSource_(const MatrixHelperRep<S>& source) {
        if (preferRowOrder_()) 
            for (int i=0; i<nrow(); ++i)
                for (int j=0; j<ncol(); ++j) {
                    if (eltIsStored_(i,j))
                        copyElt(updElt_(i,j), source.getElt_(i,j));
                }
        else // column order (rows vary fastest)
            for (int j=0; j<ncol(); ++j)
                for (int i=0; i<nrow(); ++i) {
                    if (eltIsStored_(i,j))
                        copyElt(updElt_(i,j), source.getElt_(i,j));
                }
    }

    // Overridable method to implement fillWithScalar().
    // The default implementation works but is very slow.
    virtual void fillWithScalar_(const StdNumber& scalar) {
        if (preferRowOrder_())
            for (int i=0; i<nrow(); ++i)
                for (int j=0; j<ncol(); ++j) {
                    if (eltIsStored_(i,j))
                        fillElt(updElt_(i,j), scalar);
                }
        else // column order (rows vary fastest)
            for (int j=0; j<ncol(); ++j)
                for (int i=0; i<nrow(); ++i) {
                    if (eltIsStored_(i,j))
                        fillElt(updElt_(i,j), scalar);
                }
    }

    virtual void colSum_(int j, S* csum) const {
        fillElt(csum, 0);
        for (int i=0; i < nrow(); ++i)
            addToElt(csum, getElt_(i,j));
    }

    virtual void rowSum_(int i, S* rsum) const {
        fillElt(rsum, 0);
        for (int j=0; j < ncol(); ++j)
            addToElt(rsum, getElt_(i,j));
    }
    virtual void sum_(S* esum) const {
        fillElt(esum, 0);
        S* tsum = new S[m_eltSize]; // temporary variable for row or col sums
        if (preferRowOrder_()) // i.e., row sums are cheaper
            for (int i=0; i < nrow(); ++i) 
            {   rowSum_(i, tsum); addToElt(esum, tsum); }
        else // col sums are cheaper
            for (int j=0; j < ncol(); ++j) 
            {   colSum_(j, tsum); addToElt(esum, tsum); }
        delete[] tsum;
    }



protected:
    MatrixHelperRep(int esz, int cppesz) 
    :   m_data(0), m_actual(), m_writable(false), 
        m_eltSize(esz), m_cppEltSize(cppesz), 
        m_canBeOwner(true), m_owner(false), 
        m_handleIsLocked(false), m_commitment(), m_handle(0) {}

    MatrixHelperRep(int esz, int cppesz, const MatrixCommitment& commitment) 
    :   m_data(0), m_actual(), m_writable(false),
        m_eltSize(esz), m_cppEltSize(cppesz),
        m_canBeOwner(true), m_owner(false), 
        m_handleIsLocked(false), m_commitment(commitment), m_handle(0) {}

    // Copy constructor copies just the base class members, and *not* the data.
    // We assume we don't have write access and that we aren't going to be
    // the data owner. The handle character commitment and actual character are
    // copied from the source, but may need to be changed by the caller.
    MatrixHelperRep(const MatrixHelperRep& src)
    :   m_data(0), m_actual(src.m_actual), m_writable(false),
        m_eltSize(src.m_eltSize), m_cppEltSize(src.m_cppEltSize),  
        m_canBeOwner(true), m_owner(false), 
        m_handleIsLocked(false), m_commitment(src.m_commitment), m_handle(0) {}

        // Properties of the actual matrix //

    // Raw memory holding the elements of this matrix, as an array of Scalars.
    S*                  m_data;
    // The actual characteristics of the matrix as seen through this handle.
    MatrixCharacter     m_actual;
    // Whether we have write access to the data through the pointer above. 
    // If not, the data is treated as though it were "const".
    bool                m_writable;

        // Properties of the handle //

    const int           m_eltSize;
    const int           m_cppEltSize;

    bool                m_canBeOwner;
    bool                m_owner;
    bool                m_handleIsLocked; // temporarily prevent resize of owner

    /// All commitments are by default "Uncommitted", meaning we're happy to
    /// take on matrices in any format, provided that the element types match. 
    /// Otherwise the settings here limit what we'll find acceptable as actual data.
    /// Note: don't look at these to find out anything about the current
    /// matrix. These are used only when the matrix is being created or
    /// changed structurally.
    MatrixCommitment    m_commitment;


private:
    // Point back to the owner handle of this rep.
    MatrixHelper<S>*    m_handle;

    // suppress assignment
    MatrixHelperRep& operator=(const MatrixHelperRep&);

    // For use by our friend, MatrixHelper<S>.
    void                   setMyHandle(MatrixHelper<S>& h) {m_handle = &h;}
    const MatrixHelper<S>& getMyHandle() const {assert(m_handle); return *m_handle;}
    void                   clearMyHandle() {m_handle=0;}


friend class MatrixHelperRep<typename CNT<S>::TNeg>;
friend class MatrixHelperRep<typename CNT<S>::THerm>;
friend class MatrixHelper<S>;
};


} // namespace SimTK   


#endif // SimTK_SimTKCOMMON_MATRIX_HELPER_REP_H_
