/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-15 Stanford University and the Authors.        *
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

// Avoid annoying deprecated warnings regarding std::copy during instantiations.
#ifdef _MSC_VER
#pragma warning(disable:4996)
#endif

#include "SimTKcommon/Scalar.h"
#include "SimTKcommon/SmallMatrix.h"
#include "SimTKcommon/TemplatizedLapack.h"

#include "SimTKcommon/internal/MatrixHelper.h"
#include "SimTKcommon/internal/MatrixCharacteristics.h"

#include "MatrixHelperRep.h"
#include "MatrixHelperRep_Full.h"
#include "MatrixHelperRep_Tri.h"
#include "MatrixHelperRep_Vector.h"

#include <iostream>
#include <cstdio>

namespace SimTK {


//----------------------------- MATRIX HELPER ----------------------------------
//
// Implementations of MatrixHelper methods. Most are just pass throughs to
// MatrixHelperRep, but some may involve replacing the current MatrixHelperRep
// with a different one.
//------------------------------------------------------------------------------

template <class S> void 
MatrixHelper<S>::deleteRepIfOwner() {
    if (rep && &rep->getMyHandle() == this)
        delete rep;
    rep=0;
}

// Rep replacement. Delete the existing rep if we're the owner of it.
// We are going to be the owner of the new rep regardless.
template <class S> void 
MatrixHelper<S>::replaceRep(MatrixHelperRep<S>* hrep) {
    deleteRepIfOwner();
    rep=hrep;
    rep->setMyHandle(*this);
}

// Space-stealing constructor. We're going to take over ownership
// of this rep.
template <class S> 
MatrixHelper<S>::MatrixHelper(MatrixHelperRep<S>* hrep) : rep(hrep)
{   if (rep) rep->setMyHandle(*this); }

template <class S> const MatrixCommitment& 
MatrixHelper<S>::getCharacterCommitment() const
{   return getRep().getCharacterCommitment(); }

template <class S> const MatrixCharacter& 
MatrixHelper<S>::getMatrixCharacter() const
{   return getRep().getMatrixCharacter(); }



// This is the constructor for a Matrix in which only the handle commitment has been
// supplied. The allocated matrix will have the smallest size that satisfies the
// commitment, typically 0x0 or 0x1.
template <class S> 
MatrixHelper<S>::MatrixHelper(int esz, int cppEsz, const MatrixCommitment& mc) : rep(0) {
    // Determine the best actual matrix to allocate to satisfy this commitment.
    const MatrixCharacter actual = mc.calcDefaultCharacter(0,0); 

    rep = MatrixHelperRep<S>::createOwnerMatrixHelperRep(esz, cppEsz, actual);
    assert(rep);
    rep->setMyHandle(*this);

    rep->m_commitment = mc;
}

// This is effectively the default constructor. It creates a writable, fully resizable 0x0
// matrix, with the handle committed only to the given element size.
// Just calls the above constructor with a default commitment.
template <class S> 
MatrixHelper<S>::MatrixHelper(int esz, int cppEsz) : rep(0) {
    new (this) MatrixHelper(esz, cppEsz, MatrixCommitment());
}

// This is the constructor for a Matrix in which the handle commitment has been
// supplied, along with an initial allocation size. Provided the size satisfies the
// commitment, the resulting matrix will have that size.
template <class S> 
MatrixHelper<S>::MatrixHelper
   (int esz, int cppEsz, const MatrixCommitment& commitment, 
    int m, int n) : rep(0) 
{
    SimTK_ERRCHK2(commitment.isSizeOK(m,n),  "MatrixHelper::ctor()", 
        "The initial size allocation %s x %s didn't satisfy the supplied commitment.",
        m, n);

    // Determine the best actual matrix to allocate to satisfy this commitment.
    const MatrixCharacter actual = commitment.calcDefaultCharacter(m,n);

    rep = MatrixHelperRep<S>::createOwnerMatrixHelperRep(esz, cppEsz, actual);
    assert(rep);
    rep->setMyHandle(*this);

    rep->m_commitment = commitment;
}

// Create a read-only view into existing data.
template <class S> 
MatrixHelper<S>::MatrixHelper
   (int esz, int cppEsz, const MatrixCommitment& commitment,
    const MatrixCharacter& actual, int spacing, const S* data) : rep(0) 
{
    SimTK_ERRCHK(commitment.isSatisfiedBy(actual), "MatrixHelper::ctor(external data)",
    "The supplied actual matrix character for the external data did not "
    "satisfy the specified handle character commitment.");

    rep = MatrixHelperRep<S>::createExternalMatrixHelperRep
                (esz, cppEsz, actual, spacing, const_cast<S*>(data), false); 
    assert(rep);
    rep->setMyHandle(*this);

    rep->m_commitment = commitment;
}

// Create a writable view into existing data.
template <class S> 
MatrixHelper<S>::MatrixHelper
   (int esz, int cppEsz, const MatrixCommitment& commitment, 
    const MatrixCharacter& actual, int spacing, S* data) : rep(0) 
{
    SimTK_ERRCHK(commitment.isSatisfiedBy(actual), "MatrixHelper::ctor(external data)",
        "The supplied actual matrix character for the external data did not "
        "satisfy the specified handle character commitment.");

    rep = MatrixHelperRep<S>::createExternalMatrixHelperRep
                (esz, cppEsz, actual, spacing, data, true); 
    assert(rep);
    rep->setMyHandle(*this);

    rep->m_commitment = commitment;
}

// clear() restores this matrix to the state it would be in if it were constructed
// using its current character commitment. We'll replace the current HelperRep with
// a fresh one. Note that this is more expensive than resize(0,0) which doesn't
// require replacement of the HelperRep.
template <class S> void
MatrixHelper<S>::clear() {
    // Determine the best actual matrix to allocate to satisfy this commitment.
    const MatrixCharacter actual = getCharacterCommitment().calcDefaultCharacter(0,0); 
    const int             esz    = getRep().getEltSize();
    const int             cppEsz = getRep().getCppEltSize();

    MatrixHelperRep<S>* newRep = 
        MatrixHelperRep<S>::createOwnerMatrixHelperRep(esz, cppEsz, actual);
    assert(newRep);

    newRep->m_commitment = getRep().m_commitment;
    if (!getRep().m_canBeOwner) {
        newRep->m_canBeOwner = false;
        newRep->m_owner      = false;
    }

    replaceRep(newRep);
}

template <class S> bool
MatrixHelper<S>::isClear() const {return getRep().isClear();}


template <class S> void
MatrixHelper<S>::commitTo(const MatrixCommitment& mc) {
        SimTK_ERRCHK_ALWAYS(isClear(), "MatrixHelper::commitTo()",
            "You can only replace the character commitment in an empty Matrix handle."
            "  Call clear() first to empty the handle.");
        updRep().m_commitment = mc;
        clear(); // reallocate to satisfy the new commitment
}

// Copy constructor is suppressed; these are closest things but require
// specification of deep or shallow copy.

// Create a read-only view of existing data.
template <class S>
MatrixHelper<S>::MatrixHelper(const MatrixCommitment& mc, const MatrixHelper& h, const ShallowCopy&) : rep(0) {
    SimTK_ERRCHK(mc.isSatisfiedBy(h.getMatrixCharacter()), "MatrixHelper::ctor(const,shallow)",
        "The actual matrix character of the source did not "
        "satisfy the specified handle character commitment.");
    rep = h.rep->createWholeView(false);
    rep->setMyHandle(*this);
    rep->m_commitment = mc;
}

// Create a (possibly) writable view of existing data. 
template <class S>
MatrixHelper<S>::MatrixHelper(const MatrixCommitment& mc, MatrixHelper& h, const ShallowCopy&) : rep(0) {
    SimTK_ERRCHK(mc.isSatisfiedBy(h.getMatrixCharacter()), "MatrixHelper::ctor(writable,shallow)",
        "The actual matrix character of the source did not "
        "satisfy the specified handle character commitment.");
    rep = h.rep->createWholeView(true);
    rep->setMyHandle(*this);
    rep->m_commitment = mc;
}

// Deep copy. "This" be always writable and in densely packed storage.
template <class S>
MatrixHelper<S>::MatrixHelper(const MatrixCommitment& mc, const MatrixHelper& h, const DeepCopy&) : rep(0) {
    SimTK_ERRCHK(mc.isSatisfiedBy(h.getMatrixCharacter()), "MatrixHelper::ctor(deep)",
        "The actual matrix character of the source did not "
        "satisfy the specified handle character commitment.");
    rep = h.rep->createDeepCopy();
    rep->setMyHandle(*this);
    rep->m_commitment = mc;
}

// Negated deep copy. We'll get a brand new, filterless, writable, packed copy with 
// the same values as the original (duh, that's what "copy" means) -- BUT, the
// physical floating point representations will all have been negated since we're
// copying a Matrix whose elements are of type negator<S> while ours are type S (or
// vice versa).
template <class S>
MatrixHelper<S>::MatrixHelper(const MatrixCommitment& mc, const MatrixHelper<typename CNT<S>::TNeg>& h, const DeepCopy&) : rep(0) {
    SimTK_ERRCHK(mc.isSatisfiedBy(h.getMatrixCharacter()), "MatrixHelper::ctor(negated,deep)",
        "The actual matrix character of the source did not "
        "satisfy the specified handle character commitment.");
    rep = h.rep->createNegatedDeepCopy();
    rep->setMyHandle(*this);
    rep->m_commitment = mc;
}

// construct read-only block view
template <class S> 
MatrixHelper<S>::MatrixHelper(const MatrixCommitment& mc, const MatrixHelper& h, int i, int j, int m, int n) : rep(0) {
    if (n==1)
        rep = h.rep->createColumnView(j,i,m,false); // column j, from i to i+m-1
    else if (m==1)
        rep = h.rep->createRowView(i,j,n,false); // row i, from j to j+n-1
    else 
        rep = h.rep->createBlockView(EltBlock(i,j,m,n), false);

    SimTK_ERRCHK(mc.isSatisfiedBy(rep->getMatrixCharacter()), "MatrixHelper::ctor(const,block)",
        "The actual matrix character of the source block did not "
        "satisfy the specified handle character commitment.");

    rep->setMyHandle(*this);
    rep->m_commitment = mc;
}

// construct (possibly) writable block view
template <class S> 
MatrixHelper<S>::MatrixHelper(const MatrixCommitment& mc, MatrixHelper& h, int i, int j, int m, int n) : rep(0) {
    if (n==1)
        rep = h.rep->createColumnView(j,i,m,true); // column j, from i to i+m-1
    else if (m==1)
        rep = h.rep->createRowView(i,j,n,true); // row i, from j to j+n-1
    else 
        rep = h.rep->createBlockView(EltBlock(i,j,m,n),true);

    SimTK_ERRCHK(mc.isSatisfiedBy(rep->getMatrixCharacter()), "MatrixHelper::ctor(writable,block)",
        "The actual matrix character of the source block did not "
        "satisfy the specified handle character commitment.");

    rep->setMyHandle(*this);
    rep->m_commitment = mc;
}

// Construct read only transposed view of passed-in helper.
template <class S>
MatrixHelper<S>::MatrixHelper(const MatrixCommitment& mc, const MatrixHelper<typename CNT<S>::THerm>& h,
                              const TransposeView&) : rep(0) {
    rep = h.rep->createHermitianView(false);

    SimTK_ERRCHK(mc.isSatisfiedBy(rep->getMatrixCharacter()), "MatrixHelper::ctor(const,transpose)",
        "The actual matrix character of the transposed source did not "
        "satisfy the specified handle character commitment.");

    rep->setMyHandle(*this);
    rep->m_commitment = mc;
}

// Construct (possibly) writable transposed view of passed-in helper.
template <class S>
MatrixHelper<S>::MatrixHelper(const MatrixCommitment& mc, MatrixHelper<typename CNT<S>::THerm>& h,
                              const TransposeView&) : rep(0) {
    rep = h.rep->createHermitianView(true);

    SimTK_ERRCHK(mc.isSatisfiedBy(rep->getMatrixCharacter()), "MatrixHelper::ctor(writable,transpose)",
        "The actual matrix character of the transposed source did not "
        "satisfy the specified handle character commitment.");

    rep->setMyHandle(*this);
    rep->m_commitment = mc;
}

// Construct read only diagonal view of passed-in helper.
template <class S>
MatrixHelper<S>::MatrixHelper(const MatrixCommitment& mc, const MatrixHelper& h,
                              const DiagonalView&) : rep(0) {
    rep = h.rep->createDiagonalView(false);

    SimTK_ERRCHK(mc.isSatisfiedBy(rep->getMatrixCharacter()), "MatrixHelper::ctor(const,diagonal)",
        "The actual matrix character of the diagonal of the source did not "
        "satisfy the specified handle character commitment.");

    rep->setMyHandle(*this);
    rep->m_commitment = mc;
}

// Construct (possibly) writable diagonal view of passed-in helper.
template <class S>
MatrixHelper<S>::MatrixHelper(const MatrixCommitment& mc, MatrixHelper& h,
                              const DiagonalView&) : rep(0) {
    rep = h.rep->createDiagonalView(true);

    SimTK_ERRCHK(mc.isSatisfiedBy(rep->getMatrixCharacter()), 
        "MatrixHelper::ctor(writable,diagonal)",
        "The actual matrix character of the diagonal of the source did not "
        "satisfy the specified handle character commitment.");

    rep->setMyHandle(*this);
    rep->m_commitment = mc;
}

// Construct a read-only indexed view of a Vector.
template <class S>
MatrixHelper<S>::MatrixHelper(const MatrixCommitment& mc, const MatrixHelper& h,
                              int n, const int* ix) : rep(0) {
    const VectorHelper<S>& vh = 
        SimTK_DYNAMIC_CAST_DEBUG<const VectorHelper<S>&>(*h.rep);
    rep = new IndexedVectorHelper<S>(vh.getEltSize(), vh.getCppEltSize(), n,
                                     vh.preferRowOrder_(), ix, n ? vh.getElt_(0) : 0, false);

    SimTK_ERRCHK(mc.isSatisfiedBy(rep->getMatrixCharacter()), 
        "MatrixHelper::ctor(const,indexed)",
        "The actual matrix character of the indexed source did not "
        "satisfy the specified handle character commitment.");

    rep->setMyHandle(*this);
    rep->m_commitment = mc;
}

// Construct a writable indexed view of a Vector.
template <class S>
MatrixHelper<S>::MatrixHelper(const MatrixCommitment& mc, MatrixHelper& h, 
                              int n, const int* ix) : rep(0) {
    VectorHelper<S>& vh = SimTK_DYNAMIC_CAST_DEBUG<VectorHelper<S>&>(*h.rep);
    rep = new IndexedVectorHelper<S>(vh.getEltSize(), vh.getCppEltSize(), n,
                                     vh.preferRowOrder_(), ix, n ? vh.updElt_(0) : 0, true);

    SimTK_ERRCHK(mc.isSatisfiedBy(rep->getMatrixCharacter()), 
        "MatrixHelper::ctor(writable,indexed)",
        "The actual matrix character of the indexed source did not "
        "satisfy the specified handle character commitment.");

    rep->setMyHandle(*this);
    rep->m_commitment = mc;
}

// Perform deep assignment. "This" gets reallocated if necessary if it's an
// owner, but if this is a view the source and destination sizes must match.
template <class S> MatrixHelper<S>&
MatrixHelper<S>::copyAssign(const MatrixHelper& h) {
    if (rep->isOwner())
        rep->resize(h.nrow(), h.ncol(), false);
    else {
        // OK if the sizes match.
        assert(nrow()==h.nrow() && ncol()==h.ncol());
    }

    rep->copyInFromCompatibleSource(*h.rep);
    return *this;
}
// Like copy assign but the source has elements with the opposite interpretation
// of sign. Note that the result still has the original value, but the in-memory
// representation will be different, meaning flops will be burned here.
template <class S> MatrixHelper<S>&
MatrixHelper<S>::negatedCopyAssign(const MatrixHelper<typename CNT<S>::TNeg>& h) {
    if (rep->isOwner())
        rep->resize(h.nrow(), h.ncol(), false);
    else {
        // OK if the sizes match.
        assert(nrow()==h.nrow() && ncol()==h.ncol());
    }

    rep->negatedCopyInFromCompatibleSource(*h.rep);
    return *this;
}

// Create a read-only view of existing data. Element size of source and destination
// must match. The result will have no element filter unless the source already does.
template <class S> MatrixHelper<S>&
MatrixHelper<S>::readOnlyViewAssign(const MatrixHelper& h) {
    SimTK_ERRCHK(getCharacterCommitment().isSatisfiedBy(h.getMatrixCharacter()),
        "MatrixHelper<S>::readOnlyViewAssign()",
        "The source matrix does not satisfy this handle's character commitment.");

    MatrixHelperRep<S>* newRep = h.getRep().createWholeView(false);
    newRep->m_commitment = getRep().m_commitment;   // preserve commitment
    newRep->m_canBeOwner = getRep().m_canBeOwner;
    replaceRep(newRep);
    return *this;
}
// Create a (possibly) writable view of existing data. 
template <class S> MatrixHelper<S>&
MatrixHelper<S>::writableViewAssign(MatrixHelper& h) {
    SimTK_ERRCHK(getCharacterCommitment().isSatisfiedBy(h.getMatrixCharacter()),
        "MatrixHelper<S>::writableViewAssign()",
        "The source matrix does not satisfy this handle's character commitment.");

    MatrixHelperRep<S>* newRep = h.getRep().createWholeView(true);
    newRep->m_commitment = getRep().m_commitment;   // preserve commitment
    newRep->m_canBeOwner = getRep().m_canBeOwner;
    replaceRep(newRep);
    return *this;
}

template <class S> void
MatrixHelper<S>::scaleBy(const typename CNT<S>::StdNumber& s) {
    rep->scaleBy(s);
}
template <class S> void
MatrixHelper<S>::addIn(const MatrixHelper& h) {
    rep->addIn(h);
}
template <class S> void
MatrixHelper<S>::addIn(const MatrixHelper<typename CNT<S>::TNeg>& h) {
    subIn(reinterpret_cast<const MatrixHelper<S>&>(h));
}   
template <class S> void
MatrixHelper<S>::subIn(const MatrixHelper& h) {
    rep->subIn(h);
}
template <class S> void
MatrixHelper<S>::subIn(const MatrixHelper<typename CNT<S>::TNeg>& h) {
    addIn(reinterpret_cast<const MatrixHelper<S>&>(h));
}
template <class S> void
MatrixHelper<S>::fillWith(const S* eltp) {
    rep->fillWith(eltp);
}
template <class S> void 
MatrixHelper<S>::copyInByRowsFromCpp(const S* elts) {
    rep->copyInByRowsFromCpp(elts);
}

template <class S> void
MatrixHelper<S>::invertInPlace() {
    rep->invertInPlace();
}

template <class S> void
MatrixHelper<S>::dump(const char* msg) const {
    rep->dump(msg);
}

template <class S> const S* 
MatrixHelper<S>::getElt(int i, int j) const {return rep->getElt(i,j);}
template <class S> S* 
MatrixHelper<S>::updElt(int i, int j)       {return rep->updElt(i,j);}
template <class S> void 
MatrixHelper<S>::getAnyElt(int i, int j, S* value) const
{   return rep->getAnyElt(i,j,value); }

template <class S> const S* 
MatrixHelper<S>::getElt(int i) const {return rep->getElt(i);}
template <class S> S* 
MatrixHelper<S>::updElt(int i)       {return rep->updElt(i);}
template <class S> void 
MatrixHelper<S>::getAnyElt(int i, S* value) const
{   return rep->getAnyElt(i,value); }

template <class S> int 
MatrixHelper<S>::nrow() const {
    return rep->nrow();
}
template <class S> int 
MatrixHelper<S>::ncol() const {
    return rep->ncol();
} 
template <class S> ptrdiff_t 
MatrixHelper<S>::nelt() const {   
    return rep->nelt(); 
}
template <class S> int 
MatrixHelper<S>::length() const {
    return rep->length();
}

template <class S> void 
MatrixHelper<S>::resize    (int m, int n) {rep->resize(m,n,false);} 
template <class S> void 
MatrixHelper<S>::resizeKeep(int m, int n) {rep->resize(m,n,true);} 
 
template <class S> void MatrixHelper<S>::lockShape() {rep->lockHandle();}
template <class S> void MatrixHelper<S>::unlockShape() {rep->unlockHandle();}

template <class S> void
MatrixHelper<S>::sum(S* const answer) const {
    rep->sum(answer);
}     
template <class S> void
MatrixHelper<S>::colSum(int j, S* const answer) const {
    rep->colSum(j,answer);
}
template <class S> void
MatrixHelper<S>::rowSum(int i, S* const answer) const {
    rep->rowSum(i,answer);
}
template <class S> void
MatrixHelper<S>::fillWithScalar(const typename CNT<S>::StdNumber& s) {
    rep->fillWithScalar(s);
}

template <class S> bool 
MatrixHelper<S>::hasContiguousData() const {
    return rep->hasContiguousData();
}
template <class S> ptrdiff_t 
MatrixHelper<S>::getContiguousDataLength() const {
    assert(hasContiguousData());
    return rep->nScalars();
}
template <class S> const S* 
MatrixHelper<S>::getContiguousData() const {
    assert(hasContiguousData());
    return rep->m_data;
}
template <class S> S* 
MatrixHelper<S>::updContiguousData() {
    assert(hasContiguousData());
    return rep->m_data;
}
template <class S> void 
MatrixHelper<S>::replaceContiguousData(S* newData, ptrdiff_t length, bool takeOwnership) {
    assert(length == getContiguousDataLength());
    if (rep->m_owner) {
        MatrixHelperRep<S>::deleteAllocatedMemory(rep->m_data); 
        rep->m_data=0;
    }
    rep->m_data = newData;
    rep->m_owner = takeOwnership;
}
template <class S> void 
MatrixHelper<S>::replaceContiguousData(const S* newData, ptrdiff_t length) {
    replaceContiguousData(const_cast<S*>(newData), length, false);
    rep->m_writable = false;
}
template <class S> void 
MatrixHelper<S>::swapOwnedContiguousData(S* newData, ptrdiff_t length, S*& oldData) {
    assert(length == getContiguousDataLength());
    assert(rep->m_owner);
    oldData = rep->m_data;
    rep->m_data = newData;
}




//----------------------------- MATRIX HELPER REP ------------------------------
//
//------------------------------------------------------------------------------

template <class S> MatrixHelperRep<S>*
MatrixHelperRep<S>::createOwnerMatrixHelperRep(int esz, int cppEsz, const MatrixCharacter& want) {
    SimTK_ASSERT2((esz==1&&cppEsz==1) || (esz>1&&cppEsz>=esz),
        "MatrixHelperRep::createOwnerMatrixHelperRep(): bad element size esz=%d, cppEsz=%d", esz, cppEsz);

    MatrixHelperRep<S>* rep = 0;

    const int nr = want.nrow();
    const int nc = want.ncol();
    const MatrixStructure& structure = want.getStructure();
    const MatrixStorage&   storage   = want.getStorage();
    const bool rowOrder = (storage.getOrder() == MatrixStorage::RowOrder);

    switch (structure.getStructure()) {

      case MatrixStructure::Full:
        if (rowOrder)
            rep = (esz == 1 ? (RegularFullHelper<S>*)new FullRowOrderScalarHelper<S>(nr,nc) 
                            : (RegularFullHelper<S>*)new FullRowOrderEltHelper<S>(esz, cppEsz,nr,nc));
        else
            rep = (esz == 1 ? (RegularFullHelper<S>*)new FullColOrderScalarHelper<S>(nr,nc) 
                            : (RegularFullHelper<S>*)new FullColOrderEltHelper<S>(esz, cppEsz,nr,nc));
        break;

      case MatrixStructure::Triangular:
      case MatrixStructure::Symmetric:
      case MatrixStructure::Hermitian:
      case MatrixStructure::SkewSymmetric:
      case MatrixStructure::SkewHermitian: {
        const bool triangular = structure.getStructure() == MatrixStructure::Triangular;
        const bool hermitian  =    structure.getStructure() == MatrixStructure::Hermitian
                                || structure.getStructure() == MatrixStructure::SkewHermitian;
        const bool skew       =    structure.getStructure() == MatrixStructure::SkewSymmetric
                                || structure.getStructure() == MatrixStructure::SkewHermitian;

        if (/*storage.getPlacement() == MatrixStorage::Upper*/true)
            rep = (esz == 1 ? (TriInFullHelper<S>*)new TriInFullUpperHelper<S>(1,1,nr,nc,
                                                            triangular,hermitian,skew,rowOrder) 
                            : (TriInFullHelper<S>*)new TriInFullUpperHelper<S>(esz,cppEsz,nr,nc,
                                                            triangular,hermitian,skew,rowOrder));
        //else
        //  rep = (esz == 1 ? (TriInFullHelper<S>*)new TriInFullLowerHelper<S>(1,1,nr,nc,
         //                                               false,true,false,rowOrder) 
         //                 : (TriInFullHelper<S>*)new TriInFullLowerHelper<S>(esz,cppEsz,nr,nc,
         //                                               false,true,false,rowOrder));
        break;
      }

      case MatrixStructure::Matrix1d: {
        assert(nr==1 || nc==1);
        const int length = nr*nc;
        rep = (esz==1 ? (FullVectorHelper<S>*)new ContiguousVectorScalarHelper<S>(length,rowOrder)
                      : (FullVectorHelper<S>*)new ContiguousVectorHelper<S>(esz,cppEsz,length,rowOrder));
        break;
      }
                                        
      default:
          SimTK_ERRCHK1(!"not implemented", "MatrixHelperRep::createOwnerMatrixHelperRep()",
              "Matrix structure commitment %s not implemented yet.", 
              MatrixStructure::name(structure.getStructure()));
    }

    return rep;
}


template <class S> MatrixHelperRep<S>*
MatrixHelperRep<S>::createExternalMatrixHelperRep
   (int esz, int cppEsz, const MatrixCharacter& actual,
    int spacing, S* data, bool canWrite)
{
    SimTK_ASSERT2((esz==1&&cppEsz==1) || (esz>1&&cppEsz>=esz),
        "MatrixHelperRep::createExternalMatrixHelperRep(): bad element size esz=%d, cppEsz=%d", esz, cppEsz);

    MatrixHelperRep<S>* rep = 0;

    const int nr = actual.nrow();
    const int nc = actual.ncol();
    const MatrixStructure& structure = actual.getStructure();
    const MatrixStorage&   storage   = actual.getStorage();
    const bool rowOrder = (storage.getOrder() == MatrixStorage::RowOrder);

    switch (structure.getStructure()) {

      case MatrixStructure::Full:
        if (rowOrder)
            rep = (esz == 1 ? (RegularFullHelper<S>*)new FullRowOrderScalarHelper<S>(nr,nc,
                                                            spacing,data,canWrite) 
                            : (RegularFullHelper<S>*)new FullRowOrderEltHelper<S>(esz, cppEsz,nr,nc,
                                                            spacing,data,canWrite));
        else
            rep = (esz == 1 ? (RegularFullHelper<S>*)new FullColOrderScalarHelper<S>(nr,nc,
                                                            spacing,data,canWrite) 
                            : (RegularFullHelper<S>*)new FullColOrderEltHelper<S>(esz, cppEsz,nr,nc,
                                                            spacing,data,canWrite));
        break;

      case MatrixStructure::Triangular:
      case MatrixStructure::Symmetric:
      case MatrixStructure::Hermitian:
      case MatrixStructure::SkewSymmetric:
      case MatrixStructure::SkewHermitian: {
        const bool triangular = structure.getStructure() == MatrixStructure::Triangular;
        const bool hermitian  =    structure.getStructure() == MatrixStructure::Hermitian
                                || structure.getStructure() == MatrixStructure::SkewHermitian;
        const bool skew       =    structure.getStructure() == MatrixStructure::SkewSymmetric
                                || structure.getStructure() == MatrixStructure::SkewHermitian;

        if (/*storage.getPlacement() == MatrixStorage::Upper*/true)
            rep = (esz == 1 ? (TriInFullHelper<S>*)new TriInFullUpperHelper<S>(1,1,nr,nc,
                                                            triangular,hermitian,skew,rowOrder,
                                                            spacing,data,canWrite) 
                            : (TriInFullHelper<S>*)new TriInFullUpperHelper<S>(esz,cppEsz,nr,nc,
                                                            triangular,hermitian,skew,rowOrder,
                                                            spacing,data,canWrite));
        //else
        //  rep = (esz == 1 ? (TriInFullHelper<S>*)new TriInFullLowerHelper<S>(1,1,nr,nc,
         //                                               false,true,false,rowOrder) 
         //                 : (TriInFullHelper<S>*)new TriInFullLowerHelper<S>(esz,cppEsz,nr,nc,
         //                                               false,true,false,rowOrder));
        break;
      }

      case MatrixStructure::Matrix1d: {
        assert(nr==1 || nc==1);
        const int length = nr*nc;
        const int strideInElements = spacing/esz;
        if (strideInElements > 1)
            rep = (esz==1 ? (FullVectorHelper<S>*)new StridedVectorScalarHelper<S>(length,rowOrder,
                                                            strideInElements,data,canWrite)
                          : (FullVectorHelper<S>*)new StridedVectorHelper<S>(esz,cppEsz,length,rowOrder,
                                                            strideInElements,data,canWrite));
        else 
            rep = (esz==1 ? (FullVectorHelper<S>*)new ContiguousVectorScalarHelper<S>(length,rowOrder,
                                                            data,canWrite)
                          : (FullVectorHelper<S>*)new ContiguousVectorHelper<S>(esz,cppEsz,length,rowOrder,
                                                            data,canWrite));
        break;
      }
                                        
      default:
          SimTK_ERRCHK1(!"not implemented", "MatrixHelperRep::createOwnerMatrixHelperRep()",
              "Matrix structure commitment %s not implemented yet.", 
              MatrixStructure::name(structure.getStructure()));
    }

    return rep;
}

template <class S> void
MatrixHelperRep<S>::scaleBy(const typename CNT<S>::StdNumber& s) {
    // XXX -- really, really bad! Optimize for contiguous data!
    for (int j=0; j<ncol(); ++j)
        for (int i=0; i<nrow(); ++i) 
            scaleElt(updElt(i,j),s);
}  
     
template <class S> void
MatrixHelperRep<S>::addIn(const MatrixHelper<S>& h) {
    const MatrixHelperRep& hrep = h.getRep();

    assert(nrow()==hrep.nrow() && ncol()==hrep.ncol());
    assert(getEltSize()==hrep.getEltSize());
    // XXX -- really, really bad! Optimize for contiguous data, missing views, etc.!
    for (int j=0; j<ncol(); ++j)
        for (int i=0; i<nrow(); ++i)
            addToElt(updElt(i,j),hrep.getElt(i,j));
} 
template <class S> void
MatrixHelperRep<S>::addIn(const MatrixHelper<typename CNT<S>::TNeg>& nh) {
    subIn(reinterpret_cast<const MatrixHelper<S>&>(nh));
}
     
template <class S> void
MatrixHelperRep<S>::subIn(const MatrixHelper<S>& h) {
    const MatrixHelperRep& hrep = h.getRep();

    assert(nrow()==hrep.nrow() && ncol()==hrep.ncol());
    assert(getEltSize()==hrep.getEltSize());
    // XXX -- really, really bad! Optimize for contiguous data, missing views, etc.!
    for (int j=0; j<ncol(); ++j)
        for (int i=0; i<nrow(); ++i)
            subFromElt(updElt(i,j),hrep.getElt(i,j));
}  
template <class S> void
MatrixHelperRep<S>::subIn(const MatrixHelper<typename CNT<S>::TNeg>& nh) {
    addIn(reinterpret_cast<const MatrixHelper<S>&>(nh));
}

template <class S> void
MatrixHelperRep<S>::fillWith(const S* eltp) {
    if (hasContiguousData()) {
        const ptrdiff_t len = nelt();
        if (getEltSize() == 1)
            for (ptrdiff_t i = 0; i < len; i++)
                m_data[i] = *eltp;
        else
            for (ptrdiff_t i = 0; i < len; i++)
                for (int j = 0; j < getEltSize(); j++)
                    m_data[i*getEltSize()+j] = eltp[j];
    }
    else {
        for (int j=0; j<ncol(); ++j)
            for (int i=0; i<nrow(); ++i)
                copyElt(updElt(i,j),eltp);
    }
} 

// We're copying data from a C++ row-oriented matrix into our general
// Matrix. In addition to the row ordering, C++ may use different spacing
// for elements than Simmatrix does. Lucky we know that value!
template <class S> void 
MatrixHelperRep<S>::copyInByRowsFromCpp(const S* elts) {
    const int cppRowSz = m_cppEltSize * ncol();
    // XXX -- really, really bad! Optimize for contiguous data, missing views, etc.!
    for (int j=0; j<ncol(); ++j)
        for (int i=0; i<nrow(); ++i)
            copyElt(updElt(i,j), elts + i*cppRowSz + j*m_cppEltSize);
}

template <class S> 
MatrixHelperRep<S>::~MatrixHelperRep()
{
    if (isOwner()) 
        deleteAllocatedMemory(m_data);
}

template <class S> static void
dumpElt(const S* p, int sz) {
    if (sz > 1) std::cout << "{";
    for (int k=0; k<sz; ++k) {
        if (k>0) std::cout << " ";
        std::cout << *p++;
    }
    if (sz > 1) std::cout << "}";
}

template <class S> void
MatrixHelperRep<S>::dump(const char* msg) const {
    if (msg) 
        std::cout << std::string(msg) << std::endl;
    std::cout << "Matrix " << nrow() << " X " << ncol() << " "
              << getEltSize() << "-scalar entries:" << std::endl;
    if (nrow()*ncol() == 0) {
        std::cout << "<EMPTY>" << std::endl;
        return;
    }

    S* elt = new S[getEltSize()];
    const std::streamsize oldsz = std::cout.precision(20); 
    for (int i=0; i<nrow(); ++i) {
        for (int j=0; j<ncol(); ++j) {
            if (j>0) std::cout << "\t";
            getAnyElt(i,j,elt);
            dumpElt(elt, getEltSize());
        }
        std::cout << std::endl;
    }
    std::cout.precision(oldsz);
    delete[] elt;
}

//----------------------------- RegularFullHelper ------------------------------
//------------------------------------------------------------------------------

template <class S> VectorHelper<S>* 
RegularFullHelper<S>::createDiagonalView_() {
    VectorHelper<S>* p = 0;
    const int length = std::min(this->nrow(), this->ncol());
    S*        data   = length ? this->updElt_(0,0) : 0;

    const int strideInScalars = length > 1 ? int(this->getElt_(1,1) - this->getElt_(0,0)) 
                                           : this->m_eltSize;
    const int strideInElements = strideInScalars / this->m_eltSize;

    // No need for a stride if there's 0 or 1 element, or if the stride is 1. TODO: scalar helper
    if (strideInElements == 1) {
        p = (this->m_eltSize==1) 
            ? (VectorHelper<S>*)new ContiguousVectorScalarHelper<S>(length, false, data, false)
            : (VectorHelper<S>*)new ContiguousVectorHelper<S>(this->m_eltSize, this->m_cppEltSize, 
                                                              length, false, data, false);
        return p;
    }

    p = (this->m_eltSize==1)
        ? (VectorHelper<S>*)new StridedVectorScalarHelper<S>(length, false, strideInElements, data, false)
        : (VectorHelper<S>*)new StridedVectorHelper<S>(this->m_eltSize, this->m_cppEltSize, 
                                                       length, false, strideInElements, data, false);
    return p;
}

template <class S> VectorHelper<S>* 
RegularFullHelper<S>::createColumnView_(int j, int i, int m) {
    VectorHelper<S>* p = 0;
    S* data = m ? this->updElt_(i,j) : 0;

    const int strideInScalars = m > 1 ? int(this->getElt_(i+1,j) - this->getElt_(i,j)) : this->m_eltSize;
    const int strideInElements = strideInScalars / this->m_eltSize;

    // No need for a stride if there's 0 or 1 element, or if the stride is 1. TODO: scalar helper
    if (strideInElements == 1) {
        p = (this->m_eltSize==1) 
            ? (VectorHelper<S>*)new ContiguousVectorScalarHelper<S>(m, false, data, false)
            : (VectorHelper<S>*)new ContiguousVectorHelper<S>(this->m_eltSize, this->m_cppEltSize, 
                                                              m, false, data, false);
        return p;
    }

    p = (this->m_eltSize==1)
        ? (VectorHelper<S>*)new StridedVectorScalarHelper<S>(m, false, strideInElements, data, false)
        : (VectorHelper<S>*)new StridedVectorHelper<S>(this->m_eltSize, this->m_cppEltSize, 
                                                       m, false, strideInElements, data, false);
    return p;
}

template <class S> VectorHelper<S>* 
RegularFullHelper<S>::createRowView_(int i, int j, int n) {
    VectorHelper<S>* p = 0;
    S* data = n ? this->updElt_(i,j) : 0;

    const int strideInScalars = n > 1 ? int(this->getElt_(i,j+1) - this->getElt_(i,j)) 
                                      : this->m_eltSize;
    const int strideInElements = strideInScalars / this->m_eltSize;

    // No need for a stride if there's 0 or 1 element, or if the stride is 1. TODO: scalar helper
    if (strideInElements == 1) {
        p = (this->m_eltSize==1) 
            ? (VectorHelper<S>*)new ContiguousVectorScalarHelper<S>(n, true, data, false)
            : (VectorHelper<S>*)new ContiguousVectorHelper<S>(this->m_eltSize, this->m_cppEltSize,
                                                              n, true, data, false);
        return p;
    }

    p = (this->m_eltSize==1)
        ? (VectorHelper<S>*)new StridedVectorScalarHelper<S>(n, true, strideInElements, data, false)
        : (VectorHelper<S>*)new StridedVectorHelper<S>(this->m_eltSize, this->m_cppEltSize,
                                                       n, true, strideInElements, data, false);
    return p;
}


//------------------------------ TriInFullHelper -------------------------------
//------------------------------------------------------------------------------

template <class S> VectorHelper<S>* 
TriInFullHelper<S>::createDiagonalView_() {
    SimTK_ERRCHK_ALWAYS(!hasKnownDiagonal(), "TriInFullHelper::createDiagonalView_()", 
        "Diagonal view of a known-diagonal matrix is not yet implemented. Sorry.");

    VectorHelper<S>* p = 0;
    const int length = std::min(this->nrow(), this->ncol());
    S*        data   = length ? this->updElt_(0,0) : 0;

    const int strideInScalars = length > 1 ? int(this->getElt_(1,1) - this->getElt_(0,0)) 
                                           : this->m_eltSize;
    const int strideInElements = strideInScalars / this->m_eltSize;

    // No need for a stride if there's 0 or 1 element, or if the stride is 1. TODO: scalar helper
    if (strideInElements == 1) {
        p = (this->m_eltSize==1) 
            ? (VectorHelper<S>*)new ContiguousVectorScalarHelper<S>(length, false, data, false)
            : (VectorHelper<S>*)new ContiguousVectorHelper<S>(this->m_eltSize, this->m_cppEltSize, 
                                                              length, false, data, false);
        return p;
    }

    p = (this->m_eltSize==1)
        ? (VectorHelper<S>*)new StridedVectorScalarHelper<S>(length, false, strideInElements, data, false)
        : (VectorHelper<S>*)new StridedVectorHelper<S>(this->m_eltSize, this->m_cppEltSize, 
                                                       length, false, strideInElements, data, false);
    return p;
}

//------------------------------ INSTANTIATIONS --------------------------------
// Instantiations for each of the 12 Scalar types.
// It isn't actually necessary to instantiate these classes since they
// are only used internally. However, this is a good way to make sure that
// everything supplied at least compiles. Otherwise you won't find errors
// until all the code is actually used.

#define INSTANTIATE(Helper)         \
template class Helper< float >;     \
template class Helper< double >;    \
template class Helper< std::complex<float> >;   \
template class Helper< std::complex<double> >;  \
template class Helper< conjugate<float> >;      \
template class Helper< conjugate<double> >;     \
template class Helper< negator< float > >;      \
template class Helper< negator< double > >;     \
template class Helper< negator< std::complex<float> > >;    \
template class Helper< negator< std::complex<double> > >;   \
template class Helper< negator< conjugate<float> > >;       \
template class Helper< negator< conjugate<double> > >

INSTANTIATE(MatrixHelper);
INSTANTIATE(MatrixHelperRep);

INSTANTIATE(FullHelper);
INSTANTIATE(RegularFullHelper);

INSTANTIATE(TriHelper);

INSTANTIATE(VectorHelper);
INSTANTIATE(FullVectorHelper);
INSTANTIATE(ContiguousVectorHelper);
INSTANTIATE(ContiguousVectorScalarHelper);
INSTANTIATE(StridedVectorHelper);
INSTANTIATE(StridedVectorScalarHelper);
INSTANTIATE(IndexedVectorHelper);

} // namespace SimTK   
