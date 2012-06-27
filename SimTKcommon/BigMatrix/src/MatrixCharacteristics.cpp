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

#include "SimTKcommon/internal/MatrixCharacteristics.h"

#include <iostream>

namespace SimTK {


//  ------------------------------ MatrixStructure -----------------------------
const char* 
MatrixStructure::name(Structure structure) {
    switch (structure) {
      case NoStructure:        return "No Structure";
      case Matrix1d:           return "Matrix1d";
      case Zero:               return "Zero";
      case Identity:           return "Identity";
      case Permutation:        return "Permutation";
      case RepeatedDiagonal:   return "RepeatedDiagonal";
      case Diagonal:           return "Diagonal";
      case BiDiagonal:         return "BiDiagonal";
      case TriDiagonal:        return "TriDiagonal";
      case BandedSymmetric:    return "BandedSymmetric";
      case BandedHermitian:    return "BandedHermitian";
      case Banded:             return "Banded";
      case Triangular:         return "Triangular";
      case QuasiTriangular:    return "QuasiTriangular";
      case Hessenberg:         return "Hessenberg";
      case Symmetric:          return "Symmetric";
      case Hermitian:          return "Hermitian";
      case SkewSymmetric:      return "SkewSymmetric";
      case SkewHermitian:      return "SkewHermitian";
      case Full:               return "Full";
      default: return "INVALID MatrixStructure::Structure";
    };
}

const char*
MatrixStructure::name(Position position) {
    switch (position) {
      case NoPosition:  return "No Position";
      case Upper:       return "Upper";
      case Lower:       return "Lower";
      default: return "INVALID MatrixStructure::Position";
    };
}

const char*
MatrixStructure::name(DiagValue diagValue) {
    switch (diagValue) {
      case NoDiagValue: return "No DiagValue";
      case StoredDiag:  return "StoredDiag";
      case ZeroDiag:    return "ZeroDiag";
      case UnitDiag:    return "UnitDiag";
      default: return "INVALID MatrixStructure::DiagValue";
    };
}

MatrixStructure::Mask
MatrixStructure::mask() const {
    return Mask(calcStructureMask(getStructure()),
                calcPositionMask(getStructure()),
                calcDiagValueMask(getStructure()));
}

MatrixStructure::StructureMask
MatrixStructure::calcStructureMask(Structure structure) {
    switch (structure) {
      case NoStructure:        return UncommittedStructure;
      case Matrix1d:           return Matrix1d;
      case Zero:               return Zero;
      case Identity:           return Identity;
      case Permutation:        return Permutation | calcStructureMask(Identity);
      case RepeatedDiagonal:   return RepeatedDiagonal | calcStructureMask(Identity)
                                                       | calcStructureMask(Zero);
      case Diagonal:           return Diagonal    | calcStructureMask(RepeatedDiagonal);
      case BiDiagonal:         return BiDiagonal  | calcStructureMask(Diagonal);
      case TriDiagonal:        return TriDiagonal | calcStructureMask(BiDiagonal);
      case BandedSymmetric:    return BandedSymmetric | calcStructureMask(Diagonal);
      case BandedHermitian:    return BandedHermitian | calcStructureMask(Diagonal);
      case Banded:             return Banded | calcStructureMask(BandedHermitian)
                                             | calcStructureMask(BandedSymmetric)
                                             | calcStructureMask(TriDiagonal);
      case Triangular:         return Triangular      | calcStructureMask(BiDiagonal);
      case QuasiTriangular:    return QuasiTriangular | calcStructureMask(Triangular);
      case Hessenberg:         return Hessenberg | calcStructureMask(QuasiTriangular)
                                                 | calcStructureMask(TriDiagonal);
      case Symmetric:          return Symmetric | calcStructureMask(BandedSymmetric);
      case Hermitian:          return Hermitian | calcStructureMask(BandedHermitian);
      case SkewSymmetric:      return SkewSymmetric;
      case SkewHermitian:      return SkewHermitian;
      case Full:               return Full | calcStructureMask(Hermitian) | calcStructureMask(SkewHermitian)
                                           | calcStructureMask(Symmetric) | calcStructureMask(SkewSymmetric)
                                           | calcStructureMask(Hessenberg) 
                                           | calcStructureMask(Banded)
                                           | calcStructureMask(Permutation)
                                           | calcStructureMask(Matrix1d);

      default: SimTK_ASSERT1_ALWAYS(!"invalid MatrixStructure::Structure",
                   "MatrixStructure::mask(): value 0x%x is invalid.", (unsigned)structure);
               return NoStructure;
    };
}

MatrixStructure::PositionMask
MatrixStructure::calcPositionMask(Structure structure) {
    switch (structure) {
      case BiDiagonal:
      case Triangular:
      case QuasiTriangular:
      case Hessenberg:
        return AnyPosition; // i.e., upper or lower
      default:
        return NoPosition;
    };
}

MatrixStructure::DiagValueMask
MatrixStructure::calcDiagValueMask(Structure structure) {
    switch (structure) {
      case Triangular:
      case Symmetric:
      case Hermitian:
        return AnyDiagValue;

      default:
        return StoredDiag;
    };
}


// ------------------------------ MatrixStorage --------------------------------
const char*
MatrixStorage::name(Packing p) {
    switch (p) {
      case NoPacking: return "No Packing";
      case Full:      return "Full";
      case TriInFull: return "TriInFull";
      case TriPacked: return "TriPacked";
      case Banded:    return "Banded";
      case Vector:    return "Vector";
      case Scalar:    return "Scalar";
      case Permutation: return "Permutation";
      default: return "INVALID Packing";
    };
}
const char*
MatrixStorage::name(Placement p) {
    switch (p) {
      case NoPlacement: return "No Placement";
      case Lower:       return "Lower";
      case Upper:       return "Upper";
      default: return "INVALID Placement";
    };
}
const char*
MatrixStorage::name(Order o) {
    switch (o) {
      case NoOrder:     return "No Order";
      case ColumnOrder: return "ColumnOrder";
      case RowOrder:    return "RowOrder";
      default: return "INVALID Order";
    };
}
const char*
MatrixStorage::name(Diagonal d) {
    switch (d) {
      case NoDiag:      return "No Diag";
      case StoredDiag:  return "StoredDiag";
      case AssumedDiag: return "AssumedDiag";
      default: return "INVALID Diagonal";
    };
}

MatrixStorage 
MatrixStorage::calcDefaultStorage(const MatrixStructure& structure,
                                  const MatrixOutline&   outline) 
{
    // Just to be nice we'll default to Lower storage for lower-triangular
    // matrices and Upper storage for upper-triangular ones. (There is no
    // real necessity to do that, it's just more pleasant.)
    Placement placement = Lower;
    if (structure.getPosition()==MatrixStructure::Upper) placement = Upper;

    // Default order is Lapack-standard column order except for rows.
    // TODO: probably should use row order for "Wide" outlines too, but need
    // to check if that will limit us too much with Lapack.
    MatrixStorage::Order order = ColumnOrder;
    if (outline.getOutline()==MatrixOutline::Row) 
        order = RowOrder;

    Diagonal diagonal = StoredDiag;
    if (structure.getDiagValue() && structure.getDiagValue()!=MatrixStructure::StoredDiag)
        diagonal = AssumedDiag;

    switch (structure.getStructure()) {
      case MatrixStructure::NoStructure:        return MatrixStorage();
      case MatrixStructure::Matrix1d:           return MatrixStorage(Vector, NoPlacement, order, NoDiag);
      case MatrixStructure::Zero:               return MatrixStorage(Scalar); // don't specify anything else
      case MatrixStructure::Identity:           return MatrixStorage(Scalar);
      case MatrixStructure::Permutation:        return MatrixStorage(Permutation);
      case MatrixStructure::RepeatedDiagonal:   return MatrixStorage(Scalar);
      case MatrixStructure::Diagonal:           return MatrixStorage(Vector);
      case MatrixStructure::BiDiagonal:         return MatrixStorage(Banded);
      case MatrixStructure::TriDiagonal:        return MatrixStorage(Banded);
      case MatrixStructure::BandedSymmetric:    return MatrixStorage(Banded, placement);
      case MatrixStructure::BandedHermitian:    return MatrixStorage(Banded, placement);
      case MatrixStructure::Banded:             return MatrixStorage(Banded);
      case MatrixStructure::Triangular:         return MatrixStorage(TriInFull, placement, order, diagonal);
      case MatrixStructure::QuasiTriangular:    return MatrixStorage(Full, NoPlacement, order, NoDiag);
      case MatrixStructure::Hessenberg:         return MatrixStorage(Full, NoPlacement, order, NoDiag);
      case MatrixStructure::Symmetric:          return MatrixStorage(TriInFull, placement, order, diagonal);
      case MatrixStructure::Hermitian:          return MatrixStorage(TriInFull, placement, order, diagonal);
      case MatrixStructure::SkewSymmetric:      return MatrixStorage(TriInFull, placement, order, diagonal);
      case MatrixStructure::SkewHermitian:      return MatrixStorage(TriInFull, placement, order, diagonal);
      case MatrixStructure::Full:               return MatrixStorage(Full, NoPlacement, order, NoDiag);
      default: SimTK_ASSERT1_ALWAYS(!"invalid MatrixStructure::Structure",
                   "MatrixStorage::getDefaultStorage(): value 0x%x is invalid.", structure.getStructure());
               return MatrixStorage();
    };
}


// ------------------------------- MatrixOutline -------------------------------
const char* 
MatrixOutline::name(Outline a) {
    switch(a) {
      case NoOutline:   return "No Outline";
      case Scalar:      return "Scalar";
      case Column:      return "Column";
      case Row:         return "Row";
      case Square:      return "Square";
      case Wide:        return "Wide";
      case Tall:        return "Tall";
      case Rectangular: return "Rectangular";
      default:          return "INVALID MatrixOutline::Outline";
    };
}


MatrixOutline::OutlineMask 
MatrixOutline::calcMask(Outline outline) {
    switch(outline) {
      case NoOutline:   return UncommittedOutline;
      case Scalar:      return Scalar;
      case Column:      return Column       | calcMask(Scalar);
      case Row:         return Row          | calcMask(Scalar);
      case Square:      return Square       | calcMask(Scalar);
      case Wide:        return Wide         | calcMask(Square)  | calcMask(Row);
      case Tall:        return Tall         | calcMask(Square)  | calcMask(Column);
      case Rectangular: return Rectangular  | calcMask(Wide)    | calcMask(Tall);
      default: SimTK_ASSERT1_ALWAYS(!"invalid MatrixOutline::Outline",
                   "MatrixOutline::calcMask(): value 0x%x is invalid.", outline);
          return NoOutline;
    };
}

bool 
MatrixOutline::isSizeOK(int m, int n) const {
    SimTK_ASSERT2(m>=0 && n>=0,
        "MatrixOutline::isSizeOK(): sizes must be non-negative; got %dx%d", m, n);
    switch(getOutline()) {
      case Scalar:      return m==1 && n==1;
      case Column:      return n==1;
      case Row:         return m==1;
      case Square:      return m==n;
      case Wide:        return m<=n;
      case Tall:        return m>=n;
      case Rectangular: return true;
      default: SimTK_ASSERT1_ALWAYS(!"invalid MatrixOutline::Outline",
                   "MatrixOutline::isSizeOK(): value 0x%x is invalid.", getOutline());
               return false;
    };
}

void 
MatrixOutline::getMinimumSize(int& m, int& n) const {
         if (outline==Scalar) m=n=1;
    else if (outline==Column) m=0, n=1;
    else if (outline==Row)    m=1, n=0;
    else                      m=n=0;
}

MatrixOutline
MatrixOutline::calcFromSize(int m, int n) {
    Outline outline;

         if (m==1) outline = n==1 ? Scalar : Row;
    else if (n==1) outline = Column; // m!=1
    else if (m==n) outline = Square;
    else if (m<n)  outline = Wide;
    else           outline = Tall;

    return MatrixOutline(outline);
}


// ----------------------------- MatrixCondition -------------------------------
const char* 
MatrixCondition::name(Condition c) {
    switch(c) {
      case UnknownCondition:    return "UnknownCondition";
      case Orthogonal:          return "Orthogonal";
      case PositiveDefinite:    return "PositiveDefinite";
      case WellConditioned:     return "WellConditioned";
      case FullRank:            return "FullRank";
      case Singular:            return "Singular";
      default: return "INVALID MatrixCondition::Condition";
    };
}

MatrixCondition::ConditionMask 
MatrixCondition::calcMask(Condition c) {
    switch(c) {
      case UnknownCondition:    return UncommittedCondition;
      case Orthogonal:          return Orthogonal;
      case PositiveDefinite:    return PositiveDefinite;
      case WellConditioned:     return WellConditioned
                                       | calcMask(PositiveDefinite)
                                       | calcMask(Orthogonal);
      case FullRank:            return FullRank
                                       | calcMask(WellConditioned);
      case Singular:            return UncommittedCondition; // can't get worse than this
      default: SimTK_ASSERT1_ALWAYS(!"invalid MatrixCondition::Condition",
                   "MatrixCondition::calcMask(Condition): value 0x%x is invalid.", c);
          return UnknownCondition;
    };
}

const char* 
MatrixCondition::name(Diagonal d) {
    switch(d) {
      case UnknownDiagonal:     return "UnknownDiagonal";
      case ZeroDiagonal:        return "ZeroDiagonal";
      case OneDiagonal:         return "OneDiagonal";
      case RealDiagonal:        return "RealDiagonal";
      case ImaginaryDiagonal:   return "ImaginaryDiagonal";
      default: return "INVALID MatrixCondition::Diagonal";
    };
}

MatrixCondition::DiagonalMask 
MatrixCondition::calcMask(Diagonal d) {
    switch(d) {
      case UnknownDiagonal:     return UncommittedDiagonal;
      case ZeroDiagonal:        return ZeroDiagonal;
      case OneDiagonal:         return OneDiagonal;
      case RealDiagonal:        return RealDiagonal
                                       | calcMask(ZeroDiagonal)
                                       | calcMask(OneDiagonal);
      case ImaginaryDiagonal:   return ImaginaryDiagonal
                                       | calcMask(ZeroDiagonal);
      default: SimTK_ASSERT1_ALWAYS(!"invalid MatrixCondition::Diagonal",
                   "MatrixCondition::calcMask(Diagonal): value 0x%x is invalid.", d);
          return UnknownDiagonal;
    };
}


// ----------------------------- MatrixCharacter -------------------------------
std::ostream& operator<<(std::ostream& o, const MatrixCharacter& actual) {
    o << "MatrixCharacter:";
    o << "\n  size=" << actual.nrow() << "x" << actual.ncol();
    o << "\n  structure=" << actual.getStructure().name();
    o << "\n  outline=" << actual.getOutline().name();
    o << "\n  storage=" << actual.getStorage().name();
    o << "\n  condition=" << actual.getCondition().name();
    return o << std::endl;
}


// ----------------------------- MatrixCommitment ------------------------------
MatrixCharacter
MatrixCommitment::calcDefaultCharacter(int m, int n) const {
    MatrixCharacter actual; // nothing specified yet

    // Work on the size first.
    int nr, nc;
    outline.getMinimumSize(nr, nc);

    nr = std::max(nr, getDefaultNumRows());
    nc = std::max(nc, getDefaultNumCols());

    nr = std::max(nr, m);
    nc = std::max(nc, n);

    actual.setActualSize(nr,nc);  // calculate outline also

    SimTK_ERRCHK3(isOutlineOK(actual.getOutline()) && isSizeOK(actual.getSize()),
        "MatrixCommitment::calcDefaultCharacter()",
        "Outline commitment %s and required initial size %d x %d were inconsistent.",
        outline.name().c_str(), nr, nc);

    // Size and outline are now set in actual.

    if (isStructureCommitted()) 
        actual.setStructure(structure);
    else {
        switch (getOutlineCommitment().getOutline()) {
          case MatrixOutline::Scalar: 
            actual.setStructure(MatrixStructure::RepeatedDiagonal);
            break;
          case MatrixOutline::Column:
          case MatrixOutline::Row:
            actual.setStructure(MatrixStructure::Matrix1d);
            break;
          default:
            actual.setStructure(MatrixStructure::Full);
        };
    }
    actual.updStructure().setMissingAttributes();
    const MatrixStructure& astruct = actual.getStructure();
    const MatrixOutline&   aoutline = actual.getOutline();

    // Size, outline, and structure are now set in actual.

    if (isStorageCommitted()) 
        actual.setStorage(storage);
    else {
        // Default storage format depends on structure type and outline.
        actual.setStorage(MatrixStorage::calcDefaultStorage(astruct, aoutline));
    }
    actual.updStorage().setMissingAttributes();

    // We don't expect condition to be set as a commitment usually.
    if (isConditionCommitted()) actual.setCondition(condition);
    else actual.setCondition(MatrixCondition());

    return actual;
}


std::ostream& operator<<(std::ostream& o, const MatrixCommitment& commit) {
    const MatrixCharacter::Mask::SizeMask Uncommitted = MatrixCharacter::Mask::SizeUncommitted;
    o << "MatrixCommitment:";
    o << "\n  size=";
    if (commit.getNumRowsMask()==Uncommitted) o << "Any";
    else o << commit.getNumRowsMask();
    o << " x ";
    if (commit.getNumColsMask()==Uncommitted) o << "Any";
    else o << commit.getNumColsMask();
    o << "\n  structure=" << commit.getStructureCommitment().name();
    o << "\n  outline="   << commit.getOutlineCommitment().name();
    o << "\n  storage="   << commit.getStorageCommitment().name();
    o << "\n  condition=" << commit.getConditionCommitment().name();
    return o << std::endl;
}


} // namespace SimTK   
