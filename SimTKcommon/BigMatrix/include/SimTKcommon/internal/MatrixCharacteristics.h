#ifndef SimTK_SIMMATRIX_MATRIX_CHARACTERISTICS_H_
#define SimTK_SIMMATRIX_MATRIX_CHARACTERISTICS_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-13 Stanford University and the Authors.        *
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
 * Here we declare the classes needed for managing the properties of matrices,
 * which we call Matrix Characteristics.
 */

#include "SimTKcommon/Scalar.h"

#include <iostream>
#include <cassert>
#include <complex>
#include <cstddef>
#include <utility> // for std::pair

namespace SimTK {


class MatrixStructure;
class MatrixStorage;
class MatrixOutline;
class MatrixCondition;
class MatrixCharacter;
class MatrixCommitment;


//  ------------------------------ MatrixStructure -----------------------------
/// Matrix "structure" refers to an inherent mathematical (or at least 
/// algorithmic) characteristic of the matrix rather than a storage strategy. 
/// Symmetry is the clearest example of this; it is far more significant
/// mathematically than just a way to save storage and reduce operation count.
//  ----------------------------------------------------------------------------
class SimTK_SimTKCOMMON_EXPORT MatrixStructure {
public:
    enum Structure {
        NoStructure      = 0x00000000,  ///< unspecified structure
        Matrix1d         = 0x00000001,  ///< a 1-d matrix (could be row or column)
        Zero             = 0x00000002,  ///< a matrix of all zeroes
        Identity         = 0x00000004,  ///< diagonal matrix with repeated 1's
        Permutation      = 0x00000008,  ///< permutation of an identity matrix
        RepeatedDiagonal = 0x00000010,  ///< diagonal matrix with repeated element
        Diagonal         = 0x00000020,  ///< diagonal matrix with arbitrary elements
        BiDiagonal       = 0x00000040,  ///< diagonal plus one upper or lower band
        TriDiagonal      = 0x00000080,  ///< diagonal plus one upper and one lower band
        BandedSymmetric  = 0x00000100,  ///< diagonal plus adjacent symmetric bands
        BandedHermitian  = 0x00000200,  ///< diagonal plus adjacent conjugate bands
        Banded           = 0x00000400,  ///< diagonal plus upper and lower bands
        Triangular       = 0x00000800,  ///< diagonal plus all upper or all lower bands
        QuasiTriangular  = 0x00001000,  ///< triangular but with 2x2 blocks on the diagonal
        Hessenberg       = 0x00002000,  ///< triangular plus one band above/below diagonal
        Symmetric        = 0x00004000,  ///< symmetric: elt(i,j)==elt(j,i)
        Hermitian        = 0x00008000,  ///< hermitian: elt(i,j)==conjugate(elt(j,i))
        SkewSymmetric    = 0x00010000,  ///< skew symmetric: elt(i,j) == -elt(j,i)
        SkewHermitian    = 0x00020000,  ///< skew hermitian: elt(i,j) == -conjugate(elt(j,i))
        Full             = 0x00040000   ///< full mxn matrix, all elements distinct
    };
    static const char* name(Structure);

    typedef unsigned int StructureMask; // 32 bits
    static const StructureMask AnyStructure         = 0x0007ffffU; // see above
    static const StructureMask UncommittedStructure = 0xffffffffU;
    static StructureMask calcStructureMask(Structure);

    /// For triangular matrices, we have to know which triangle 
    /// we're talking about. Don't confuse this with MatrixStorage::Placement
    /// which has to do with where we put it in memory, *not* what
    /// matrix is being represented.
    enum Position {
        NoPosition = 0x0000,
        Lower      = 0x0001,  // matrix is lower triangular (default)
        Upper      = 0x0002   // matrix is upper triangular
    };
    static const char* name(Position);

    typedef unsigned short PositionMask; // 16 bits
    static const PositionMask AnyPosition         = 0x0003U;  // see above
    static const PositionMask UncommittedPosition = 0xffffU;
    static PositionMask calcPositionMask(Structure);

    /// For triangular, symmetric, and hermitian matrices the diagonal elements
    /// may have a single, assumed value rather than being stored in memory.
    /// This specifies the value. Don't confuse this with the similar
    /// MatrixStorage type which simply says whether there is a known value
    /// rather than stored values, *not* what that value is.
    enum DiagValue {
        NoDiagValue = 0x0000,
        StoredDiag  = 0x0001, // could be anything (default)
        ZeroDiag    = 0x0002, // zero (e.g. for skew matrices)
        UnitDiag    = 0x0004  // unit (one) diagonal is used frequently by Lapack
    };
    static const char* name(DiagValue);

    typedef unsigned short DiagValueMask; // 16 bits
    static const DiagValueMask AnyDiagValue         = 0x0003U;
    static const DiagValueMask UncommittedDiagValue = 0xffffU;
    static DiagValueMask calcDiagValueMask(Structure);

    MatrixStructure& setMissingAttributes() {
        if (structure == NoStructure)
            structure = Full;
        if (position == NoPosition)
            position = Lower;
        if (diagValue == NoDiagValue)
            diagValue = StoredDiag;
        return *this;
    }

    std::string name() const {
        return std::string(name(getStructure())) 
            + "|" + std::string(name(getPosition()))
            + "|" + std::string(name(getDiagValue()));
    }

    struct Mask {
        Mask() {setToUncommitted();}
        Mask(StructureMask sm, PositionMask pm, DiagValueMask dm)
        :   structure(sm), position(pm), diagValue(dm) {}
        Mask& setToUncommitted() 
        {   structure=UncommittedStructure; position=UncommittedPosition; 
            diagValue=UncommittedDiagValue; return *this; }
        bool isUncommitted() const
        {   return structure==UncommittedStructure && position==UncommittedPosition 
                && diagValue==UncommittedDiagValue; }
        bool isSatisfiedBy(Structure str, Position pos, DiagValue diag) const
        {   return ((StructureMask)str&structure)==(StructureMask)str 
                && ((PositionMask)pos&position)==(PositionMask)pos
                && ((DiagValueMask)diag&diagValue)==(DiagValueMask)diag; }
        bool isSatisfiedBy(const MatrixStructure& actual) const
        {  return isSatisfiedBy(actual.getStructure(), actual.getPosition(), 
                                actual.getDiagValue()); }

        StructureMask  structure;
        PositionMask   position;
        DiagValueMask  diagValue;
    };

    MatrixStructure() {setToNone();}

    /// This constructor is also an implicit conversion from the Structure enum
    /// to a MatrixStructure object which does not specify Position or DiagValue.
    MatrixStructure(Structure s, Position p=NoPosition, DiagValue d=NoDiagValue)
        :   structure(s), position(p), diagValue(d) {} 

    /// Given a Structure commitment, which more-restrictive Structures will
    /// still satisfy this commitment? Returned value is a mask with a bit
    /// set for every Structure that is satisfactory. For example, if the
    /// commitment is "Banded", "Diagonal" is also acceptable.
    Mask mask() const;

    Structure   getStructure() const {return structure;}
    Position    getPosition()  const {return position;}
    DiagValue   getDiagValue() const {return diagValue;}

    MatrixStructure& setStructure(Structure s) {structure=s; return *this;}
    MatrixStructure& setPosition (Position  p) {position=p; return *this;}
    MatrixStructure& setDiagValue(DiagValue d) {diagValue=d; return *this;}

    MatrixStructure& set(Structure s, Position p, DiagValue d)
    {   structure=s; position=p; diagValue=d; return *this; }

    MatrixStructure& setToNone() 
    {   structure=NoStructure; position=NoPosition; 
        diagValue=NoDiagValue; return *this; }

private:
    Structure  structure:32;
    Position   position:16;
    DiagValue  diagValue:16;
};


//  ------------------------------ MatrixStorage -------------------------------
/// Matrix "storage" refers to the physical layout of data in the computer’s 
/// memory. Whenever possible we attempt to store data in a format that enables 
/// use of special high performance methods, such as those available in the 
/// SimTK LAPACK/BLAS implementation.
//  ----------------------------------------------------------------------------
class SimTK_SimTKCOMMON_EXPORT MatrixStorage {
public:
    enum Packing {
        NoPacking    = 0x0000,
        Full         = 0x0001,  // full storage layout
        TriInFull    = 0x0002,  // a triangular piece of a full storage layout
        TriPacked    = 0x0004,  // triangle packed into minimal storage, at performance cost
        Banded       = 0x0008,  // a packed, banded storage format
        Vector       = 0x0010,  // a possibly-strided or scattered vector
        Scalar       = 0x0020,  // a single scalar is stored
        Permutation  = 0x0040   // a permuted identity matrix
    };
    static const char* name(Packing);
    typedef unsigned short PackingMask;
    static const PackingMask AllPacking = 0x007fU; // see above
    static const PackingMask UncommittedPacking = 0xffffU;

    enum Placement {
        NoPlacement  = 0x0000,
        Lower        = 0x0001,  // stored in lower triangle of full matrix
        Upper        = 0x0002,  // stored in upper triangle of full matrix
    };
    static const char* name(Placement);
    typedef unsigned short PlacementMask;
    static const PlacementMask AllPlacement = 0x0003U; // see above
    static const PlacementMask UncommittedPlacement = 0xffffU;

    enum Order {
        NoOrder      = 0x0000,
        ColumnOrder  = 0x0001,  // matrix is stored by columns
        RowOrder     = 0x0002,  // matrix is stored by rows
    };
    static const char* name(Order);
    typedef unsigned short OrderMask;
    static const OrderMask AllOrder = 0x03U; // see above
    static const OrderMask UncommittedOrder = 0xffU;

    enum Diagonal {
        NoDiag       = 0x0000,
        StoredDiag   = 0x0001,  // matrix diagonal is stored
        AssumedDiag  = 0x0002   // matrix diagonal is not stored but has known value
    };
    static const char* name(Diagonal);
    typedef unsigned short DiagonalMask;
    static const DiagonalMask AllDiagonal = 0x0003U; // see above
    static const DiagonalMask UncommittedDiagonal = 0xffffU;

    /// Use this class to represent sets of acceptable values for each of 
    /// the storage attributes (packing, position, order, diagonal).
    struct Mask {
        Mask()
        :   packing(UncommittedPacking), placement(UncommittedPlacement), 
            order(UncommittedOrder), diagonal(UncommittedDiagonal) {}
        Mask(PackingMask pkm, PlacementMask plm, OrderMask om, DiagonalMask dm)
        :   packing(pkm), placement(plm), order(om), diagonal(dm) {}
        Mask& setToUncommitted()
        {   packing=UncommittedPacking; placement=UncommittedPlacement; 
            order=UncommittedOrder;     diagonal=UncommittedDiagonal; return *this; }
        bool isUncommitted() const 
        {   return packing==UncommittedPacking && placement==UncommittedPlacement 
                && order==UncommittedOrder     && diagonal==UncommittedDiagonal; }        
        bool isSatisfiedBy(Packing pack, Placement place, Order ord, Diagonal diag) const
        {   return ((PackingMask)pack    & packing)   == (PackingMask)  pack 
                && ((PlacementMask)place & placement) == (PlacementMask)place
                && ((OrderMask)ord       & order)     == (OrderMask)    ord     
                && ((DiagonalMask)diag   & diagonal)  == (DiagonalMask) diag; }
        bool isSatisfiedBy(const MatrixStorage& actual) const
        {   return isSatisfiedBy(actual.getPacking(), actual.getPlacement(), 
                                 actual.getOrder(),   actual.getDiagonal());}      

        PackingMask   packing;
        PlacementMask placement;
        OrderMask     order;
        DiagonalMask  diagonal;
    };

    static MatrixStorage calcDefaultStorage(const MatrixStructure&,
                                            const MatrixOutline&);

    std::string name() const {
        return std::string(name(getPacking()))
            + "|" + std::string(name(getPlacement()))
            + "|" + std::string(name(getOrder()))
            + "|" + std::string(name(getDiagonal()));
    }

    /// Calculate the commitment mask associated with specifying "this" set
    /// of storage attributes as a commitment. Here the mask will either be
    /// fully uncommitted or set to a specific value for each attribute; they
    /// are all mutually exclusive.
    Mask mask() const {
        Mask ms; // initially uncommitted
        if (packing)   ms.packing   = (PackingMask)packing;
        if (placement) ms.placement = (PlacementMask)placement;
        if (order)     ms.order     = (OrderMask)order;
        if (diagonal)  ms.diagonal  = (DiagonalMask)diagonal;
        return ms;
    }

    /// Default constructor leaves all fields unspecified.
    MatrixStorage() 
    :   packing(NoPacking), placement(NoPlacement), order(NoOrder), diagonal(NoDiag) {}

    /// This constructor is also an implicit conversion from the Packing enum to a
    /// MatrixStorage object which does not contain any specification for placement,
    /// order, or storage of diagonal elements.
    MatrixStorage(Packing pk, Placement pl=NoPlacement, Order o=NoOrder, Diagonal d=NoDiag)
    :   packing(pk), placement(pl), order(o), diagonal(d) {}

    /// This constructor is for the common case of just packing and order, with
    /// no particular placement and a stored diagonal.
    MatrixStorage(Packing pk, Order o)
    :   packing(pk), placement(NoPlacement), order(o), diagonal(StoredDiag) {}

    /// Assuming this is an actual matrix description, set any unspecified attributes
    /// to appropriate defaults to match the specified packing.
    MatrixStorage& setMissingAttributes() {
        if (packing==NoPacking) 
            packing = Full;
        if (placement==NoPlacement)
            placement = Lower;
        if (order==NoOrder)
            order = ColumnOrder;
        if (diagonal==NoDiag)
            diagonal = StoredDiag;
        return *this;
    }

    /// Restore this object to its default-constructed state of "none".
    MatrixStorage& setToNone()
    {   packing=NoPacking; placement=NoPlacement; 
        order=NoOrder;     diagonal=NoDiag; return *this; }

    MatrixStorage& setPacking(Packing p)     {packing   = p; return *this;}
    MatrixStorage& setPlacement(Placement p) {placement = p; return *this;}
    MatrixStorage& setOrder(Order o)         {order     = o; return *this;}
    MatrixStorage& setDiagonal(Diagonal d)   {diagonal  = d; return *this;}

    Packing   getPacking()   const {return packing;}
    Placement getPlacement() const {return placement;}
    Order     getOrder()     const {return order;}
    Diagonal  getDiagonal()  const {return diagonal;}

private:
    Packing   packing:16;
    Placement placement:16;
    Order     order:16;
    Diagonal  diagonal:16;
};


//  ------------------------------- MatrixOutline ------------------------------
/// Matrix "outline" refers to the characteristic relationship between the number
/// of rows and columns of a matrix, without necessarily specifying the absolute
/// dimensions.
///
/// There are seven possible outline attributes:
/// - Rectangular (any m x n)
/// - Tall (m >= n)
/// - Wide (m <= n)
/// - Square (m == n)
/// - Column (m x 1)
/// - Row (1 x n)
/// - Scalar (1 x 1)
/// 
/// A matrix handle with no outline commitment can hold a general (rectangular) 
/// matrix, which of course includes all the other outlines as well. Vector 
/// handles are always committed to column outline, RowVector handles to row 
/// outline. Scalar matrices are also rows, columns, and square and some
/// rows and columns are square (if they are 1x1).
///
/// Note that certain outlines imply fixed sizes of one or both dimensions.
//  ----------------------------------------------------------------------------
class SimTK_SimTKCOMMON_EXPORT MatrixOutline {
public:
    enum Outline {
        NoOutline   = 0x0000,
        Scalar      = 0x0001,    // 1x1
        Column      = 0x0002,    // mx1, m != 1
        Row         = 0x0004,    // 1xn, n != 1
        Square      = 0x0008,    // mxn, m == n
        Wide        = 0x0010,    // mxn, m < n
        Tall        = 0x0020,    // mxn, m > n
        Rectangular = 0x0040     // mxn
    };
    static const char* name(Outline);

    typedef unsigned short OutlineMask;
    static const OutlineMask AnyOutline  = 0x007fU; // see above
    static const OutlineMask UncommittedOutline = 0xffffU;

    struct Mask {
        Mask() : outline(UncommittedOutline) {}
        explicit Mask(OutlineMask mask) : outline(mask) {}
        Mask& setToUncommitted() {outline=UncommittedOutline; return *this;}
        bool isUncommitted() const {return outline==UncommittedOutline;}
        bool isSatisfiedBy(const MatrixOutline& actual) const 
        {  return ((OutlineMask)actual.outline & outline) == (OutlineMask)actual.outline; }

        OutlineMask outline;
    };

    std::string name() const {return std::string(name(getOutline()));}

    /// Default constructor produces an object containing no outline specification.
    /// If used as a commitment it won't accept \e any outline!
    MatrixOutline() : outline(NoOutline) {}

    /// This is an implicit conversion from the Outline enum to a MatrixOutline object.
    MatrixOutline(Outline outline) : outline(outline) {}

    /// Set the outline back to its default-constructed value of "none".
    MatrixOutline& setToNone() {outline=NoOutline; return *this;}

    /// Compute a mask of acceptable Outline values given a particular
    /// value specified as a commitment. For example, if the commitment 
    /// is "Wide" then Square, Row, and Scalar outlines are also acceptable.
    static OutlineMask calcMask(Outline);

    /// When "this" outline is used as a commitment, it represents a mask of
    /// acceptable outlines. Calculate and return that mask.
    Mask mask() const {return Mask(calcMask(getOutline()));}

    /// Determine if the proposed shape satisfies this outline.
    bool isSizeOK(int m, int n) const;

    /// Return the minimum shape that will satisfy this outline.
    void getMinimumSize(int& m, int& n) const;

    /// Determine the outline from given actual dimensions.
    static MatrixOutline calcFromSize(int m, int n);

    /// Return the outline value stored in this MatrixOutline object.
    Outline getOutline() const {return outline;}

private:
    Outline outline:16;
};



//  ---------------------------- MatrixCondition -------------------------------
/// Matrix "condition" is a statement about the numerical characteristics of a 
/// Matrix. It can be set as a result of an operation, or by a knowledgeable user. 
/// Simmatrix is entitled to rely on the correctness of these assertions, 
/// although it will try to check them if requested or when it is possible to do 
/// so without sacrificing performance.
/// In most cases the condition will not be set at all, which we interpret to
/// mean "unknown condition" in an actual character, and "uncommitted" in
/// a character commitment.
//  ----------------------------------------------------------------------------
class SimTK_SimTKCOMMON_EXPORT MatrixCondition {
public:
    enum Condition {
        UnknownCondition = 0x0000,
        Orthogonal       = 0x0001, // implies well conditioned
        PositiveDefinite = 0x0002, // implies well conditioned
        WellConditioned  = 0x0004, // implies full rank
        FullRank         = 0x0008, // but might have bad conditioning
        Singular         = 0x0010  // implies possible bad conditioning 
    };
    static const char* name(Condition);

    typedef unsigned short ConditionMask;   // 16 bits in mask
    static const ConditionMask AnyCondition          = 0x001fU;  // see above
    static const ConditionMask UncommittedCondition  = 0xffffU;

    enum Diagonal {
        UnknownDiagonal   = 0x0000,   ///< we don't know the condition of the diagonal elements
        ZeroDiagonal      = 0x0001,   ///< all diagonal elements are zero (also pure real and pure imaginary)
        OneDiagonal       = 0x0002,   ///< all diagonal elements are one (and real)
        RealDiagonal      = 0x0004,   ///< all diagonal elements are pure real
        ImaginaryDiagonal = 0x0008    ///< all diagonal elements are pure imaginary
    };
    static const char* name(Diagonal);

    typedef unsigned short DiagonalMask;   // 16 bits in mask
    static const DiagonalMask AnyDiagonal          = 0x000fU;  // see above
    static const DiagonalMask UncommittedDiagonal  = 0xffffU;

    /// Use this class to represent a set of acceptable Condition values.
    struct Mask {
        Mask() : condition(UncommittedCondition), diagonal(UncommittedDiagonal) {}
        Mask(ConditionMask cmask, DiagonalMask dmask) : condition(cmask), diagonal(dmask) {}
        Mask& setToUncommitted()    
        {   condition=UncommittedCondition; diagonal=UncommittedDiagonal; return *this;}
        bool isUncommitted() const 
        {   return condition==UncommittedCondition && diagonal==UncommittedDiagonal;}
        bool isSatisfiedBy(const MatrixCondition& actual) const 
        {   return ((ConditionMask)actual.condition & condition) == (ConditionMask)actual.condition
                && ((DiagonalMask) actual.diagonal  & diagonal)  == (DiagonalMask)actual.diagonal; }

        ConditionMask   condition;
        DiagonalMask    diagonal;
    };

    std::string name() const 
    {   return std::string(name(getCondition())) + "|" + std::string(name(getDiagonal()));}

    /// The default constructor sets the condition to Unknown, which is typically
    /// where it remains.
    MatrixCondition() : condition(UnknownCondition), diagonal(UnknownDiagonal) {}

    /// This is an implicit conversion from the Condition enum to a
    /// MatrixCondition object.
    MatrixCondition(Condition cond, Diagonal diag=UnknownDiagonal) 
    :   condition(cond), diagonal(diag) {}

    /// Restore to default-constructed state of "none".
    MatrixCondition& setToNone() {condition=UnknownCondition; diagonal=UnknownDiagonal; return *this;}

    /// Given a particular Condition provided as a commitment, calculate
    /// the mask of all Condition values that would satisfy that commitment.
    /// For example, if the commitment is "WellConditioned" then "Orthogonal"
    /// and "PositiveDefinite" also qualify since they are automatically
    /// well conditioned. If the Condition is "Unknown" we'll return an Uncommitted mask.
    static ConditionMask calcMask(Condition);

    /// Given a particular Diagonal condition provided as a commitment, calculate
    /// the mask of all Diagonal conditions that would satisfy that commitment.
    /// For example, if the commitment is "RealDiagonal" then "ZeroDiagonal"
    /// and "OneDiagonal" would also be acceptable. If the Diagonal condition
    /// is specified as "UnknownDiagonal" then we'll return an Uncommitted mask.
    static DiagonalMask calcMask(Diagonal);

    /// Return the commitment mask corresponding to use of "this" condition
    /// as a commitment.
    Mask mask() const 
    {   return Mask(calcMask(getCondition()), calcMask(getDiagonal())); }

    Condition getCondition() const {return condition;}
    Diagonal  getDiagonal()  const {return diagonal;}

    MatrixCondition& setCondition(Condition c) {condition=c; return *this;}
    MatrixCondition& setDiagonal (Diagonal d)  {diagonal=d; return *this;}

private:
    Condition condition:16;
    Diagonal  diagonal:16;
};



//  ------------------------------ MatrixCharacter -----------------------------
/** A MatrixCharacter is a set containing a value for each of the matrix 
characteristics except element type, which is part of the templatized
declaration of a Matrix_, Vector_, or RowVector_ handle. MatrixCharacters are
used both as the handle "commitment", setting restrictions on what kinds
of matrices a handle can reference, and as the "facts on the ground" current
character of the matrix being referenced. The current character must always
satisfy the character commitment.

%Matrix characteristics are specifications of particular aspects of matrices:
 - Element type
 - Size
 - Structure
 - Storage format
 - Outline
 - Conditioning

Collectively, the set of values for the above properties is called a <i>matrix
character</i>. A matrix character can be used to describe an existing matrix,
or a character mask can be used to describe the range of characteristics that
a matrix handle may support. The  character mask describing the acceptable
matrices for a matrix handle is called the handle's <i>character commitment</i>
or just the <i>handle commitment</i>. The character describing an existing 
matrix is called the <i>actual character</i> of that matrix. Thus there are 
always two sets of characteristics associated with a matrix: the handle's 
commitment, and the actual character of the matrix to which the handle currently
refers. The actual character must always \e satisfy the character commitment.

When a handle presents a view into another handle's data, it is the 
characteristics of the matrix as seen through the view that must satisfy the 
handle's character commitment. So for example, a view showing one column of a 
full matrix satisfies a "column" outline commitment.

Element type for a matrix handle is always determined at compile time via the
template argument used in the declaration. For example, a matrix handle declared
\c Matrix_<Vec3> can only hold matrices whose elements are \c Vec3s. Also,
recall that \c Matrix is an abbreviation for \c Matrix_<Real> so that
declaration commits the matrix handle to Real-element matrices. Element type
is the only matrix characteristic for which no matrix handle can remain
uncommitted. However, different handles can provide views of the same data
through which that data is seen to contain different element types.

Each matrix characteristic other than sizes is represented by a class
defining one or more enumerated types, where individual characteristics are
assigned a single bit. Then an appropriate mask type (an unsigned integral
type) is defined which can represent a set of allowable characteristics.
The actual character of a matrix is represented via enumeration values; the
character commitment is represented by the compatible masks. The operation
of determining whether a particular actual character satisfies a handle
commitment can then be performed very quickly via bitwise logical operations.
**/
class SimTK_SimTKCOMMON_EXPORT MatrixCharacter {
public:
    /// Default constructor sets lengths to zero and the other characteristics
    /// to "none specified".
    MatrixCharacter() : nr(0), nc(0), lband(0), uband(0) {}

    // Some handy predefined MatrixCharacters.
    class LapackFull;
    class Vector;
    class RowVector;

    /// Restore this MatrixCharacter to its default-constructed state of "none".
    MatrixCharacter& setToNone() {
        nr=nc=lband=uband=0;
        structure.setToNone(); outline.setToNone();
        storage.setToNone();   condition.setToNone();
        return *this;
    }

    /// These are dimensions of the logical matrix and have nothing to do with 
    /// how much storage may be used to hold the elements.
    int                nrow()       const {return nr;}
    int                ncol()       const {return nc;}
    std::pair<int,int> getSize()    const {return std::pair<int,int>(nrow(),ncol());}
    ptrdiff_t          nelt()       const {return (ptrdiff_t)nrow() * (ptrdiff_t)ncol();}

    int                getLowerBandwidth() const {return lband;}
    int                getUpperBandwidth() const {return uband;}
    std::pair<int,int> getBandwidth()      const 
    {   return std::pair<int,int>(getLowerBandwidth(), getUpperBandwidth()); }

    const MatrixStructure&  getStructure() const {return structure;}
    const MatrixStorage&    getStorage()   const {return storage;}
    const MatrixOutline&    getOutline()   const {return outline;}
    const MatrixCondition&  getCondition() const {return condition;}

    MatrixStructure&  updStructure() {return structure;}
    MatrixStorage&    updStorage()   {return storage;}
    MatrixOutline&    updOutline()   {return outline;}
    MatrixCondition&  updCondition() {return condition;}

    MatrixCharacter& setStructure(const MatrixStructure& sa)  {structure = sa; return *this;}
    MatrixCharacter& setStorage  (const MatrixStorage&   sa)  {storage   = sa; return *this;}
    MatrixCharacter& setOutline  (const MatrixOutline&   oa)  {outline   = oa; return *this;}
    MatrixCharacter& setCondition(const MatrixCondition& ca)  {condition = ca; return *this;}


    /// Set the actual size and update the outline to match.
    MatrixCharacter& setActualSize(int m, int n)
    {   setSize(m,n); outline = MatrixOutline::calcFromSize(m,n); return *this; }
    MatrixCharacter& setActualNumRows(int m)
    {   setNumRows(m); outline = MatrixOutline::calcFromSize(m,ncol()); return *this; }
    MatrixCharacter& setActualNumCols(int n)
    {   setNumCols(n); outline = MatrixOutline::calcFromSize(nrow(),n); return *this; }

    MatrixCharacter& setBandwidth(int lb, int ub) {
        assert(lb>=0 && lb>=0);
        lband = lb; uband = ub;
        return *this;
    }
    MatrixCharacter& setLowerBandwidth(int lb) {
        assert(lb>=0);
        lband = lb;
        return *this;
    }
    MatrixCharacter& setUpperBandwidth(int ub) {
        assert(ub>=0);
        uband = ub;
        return *this;
    }

    class Mask; // defined below

protected:
    MatrixCharacter(int m, int n,
                    int lb, int ub,
                    MatrixStructure structure,
                    MatrixStorage   storage,
                    MatrixCondition condition)
    :   nr(m), nc(n), lband(lb), uband(ub),
        structure(structure), storage(storage), 
        outline(MatrixOutline::calcFromSize(m,n)), 
        condition(condition) {}


    int              nr,            ///< actual number of rows
                     nc;            ///< actual number of columns
    int              lband,         ///< actual lower bandwidth, if banded
                     uband;         ///< actual upper bandwidth, if banded
    MatrixStructure  structure;
    MatrixStorage    storage;
    MatrixOutline    outline;
    MatrixCondition  condition;

private:
    // These are private because they don't set the outline as well.
    MatrixCharacter& setSize(int m, int n) 
    {   assert(m>=0 && n>=0); nr = m; nc = n; return *this; }
    MatrixCharacter& setNumRows(int m) 
    {   assert(m>=0); nr = m; return *this; }
    MatrixCharacter& setNumCols(int n) 
    {   assert(n>=0); nc = n; return *this; }
};

/// Output a textual description of a MatrixCharacter; handy for debugging.
/// @relates SimTK::MatrixCharacter
SimTK_SimTKCOMMON_EXPORT std::ostream& 
operator<<(std::ostream& o, const MatrixCharacter&);

/// Predefined MatrixCharacter for an ordinary Lapack-style full matrix
/// of a particular dimension m x n (nrows X ncols). Note that the storage
/// format allows for a "leading dimension" larger than m, but that the 
/// leading dimension is not considered part of the matrix character. It is 
/// dealt with separately.
class MatrixCharacter::LapackFull : public MatrixCharacter {
public:
    LapackFull(int m, int n)
    :   MatrixCharacter(m,n,0,0,
            MatrixStructure(MatrixStructure::Full),
            MatrixStorage(MatrixStorage::Full,MatrixStorage::ColumnOrder),
            MatrixCondition()) {}                   
};

/// Predefined MatrixCharacter for an ordinary column vector of a particular
/// size. Note that standard vector storage allows for a "stride" greater than
/// 1, but that is not considered part of the matrix character. Stride is
/// dealt with separately.
class MatrixCharacter::Vector : public MatrixCharacter {
public:
    Vector(int m)
    :   MatrixCharacter(m,1,0,0,
            MatrixStructure(MatrixStructure::Matrix1d),
            MatrixStorage(MatrixStorage::Vector,MatrixStorage::ColumnOrder),
            MatrixCondition()) {}                   
};

/// Predefined MatrixCharacter for an ordinary row vector of a particular
/// size. Note that standard vector storage allows for a "stride" greater than
/// 1, but that is not considered part of the matrix character. Stride is
/// dealt with separately.
class MatrixCharacter::RowVector : public MatrixCharacter {
public:
    RowVector(int n)
    :   MatrixCharacter(1,n,0,0,
            MatrixStructure(MatrixStructure::Matrix1d),
            MatrixStorage(MatrixStorage::Vector,MatrixStorage::RowOrder),
            MatrixCondition()) {}                   
};

//  -------------------------- MatrixCharacter::Mask ---------------------------
/// This class collects masks of each characteristic type for representing sets
/// of accceptable characteristics.
//  ----------------------------------------------------------------------------
class MatrixCharacter::Mask {
public:
    Mask() {setToUncommitted();}

    typedef unsigned int SizeMask;
    static const SizeMask SizeUncommitted = 0xffffffffU;

    bool isResizeable()      const {return nr==SizeUncommitted || nc==SizeUncommitted;}
    bool isFullyResizeable() const {return nr==SizeUncommitted && nc==SizeUncommitted;}
    bool isNumRowsLocked()   const {return nr!=SizeUncommitted;}
    bool isNumColsLocked()   const {return nc!=SizeUncommitted;}

    unsigned int getNumRowsMask() const {return nr;}
    unsigned int getNumColsMask() const {return nc;}
    unsigned int getLowerBandwidthMask() const {return lband;}
    unsigned int getUpperBandwidthMask() const {return uband;}

    int getDefaultNumRows() const {return isNumRowsLocked() ? nr : 0;}
    int getDefaultNumCols() const {return isNumColsLocked() ? nc : 0;}

    bool isLowerBandwidthLocked()  const {return lband!=SizeUncommitted;}
    bool isUpperBandwidthLocked()  const {return uband!=SizeUncommitted;}
    int getDefaultLowerBandwidth() const {return isLowerBandwidthLocked() ? lband : 0;}
    int getDefaultUpperBandwidth() const {return isUpperBandwidthLocked() ? uband : 0;}

    /// Set all bits to one ("Uncommitted").
    Mask& setToUncommitted() {
        nr=nc=lband=uband=SizeUncommitted;
        structure.setToUncommitted(); storage.setToUncommitted();
        outline.setToUncommitted();   condition.setToUncommitted();
        return *this;
    }

    /// Return if all fields are set to "Uncommitted" (all bits are one).
    bool isUncommitted() const {
        return nr==SizeUncommitted       && nc==SizeUncommitted 
            && lband==SizeUncommitted    && uband==SizeUncommitted
            && structure.isUncommitted() && storage.isUncommitted()
            && outline.isUncommitted()   && condition.isUncommitted();
    }

    /// Check whether an actual matrix character satisfies this matrix commitment.
    bool isSatisfiedBy(const MatrixCharacter& actual) const {
        return isSizeOK(actual.nr, actual.nc) 
            && isBandwidthOK(actual.lband, actual.uband)
            && structure.isSatisfiedBy(actual.getStructure())
            && storage.isSatisfiedBy(actual.getStorage())
            && outline.isSatisfiedBy(actual.getOutline())
            && condition.isSatisfiedBy(actual.getCondition());
    }

    /// Check whether an actual size satisfies the size commitment.
    bool isSizeOK(int m, int n) const 
    {   return ((SizeMask)m & nr)      == (SizeMask)m
            && ((SizeMask)n & nc)      == (SizeMask)n; }

    /// Check whether an actual bandwidth satisfies the bandwidth commitment.
    /// (If the matrix isn't banded any bandwidth will be OK.)
    bool isBandwidthOK(int lower, int upper) const 
    {   return ((SizeMask)lower & lband) == (SizeMask)lower
            && ((SizeMask)upper & uband) == (SizeMask)upper; }

    SizeMask                nr,         ///< number of rows
                            nc;         ///< number of columns
    SizeMask                lband,      ///< lower bandwidth, if banded
                            uband;      ///< upper bandwidth, if banded
    MatrixStructure::Mask   structure;
    MatrixStorage::Mask     storage;
    MatrixOutline::Mask     outline;
    MatrixCondition::Mask   condition;

friend class MatrixCommitment;
};

//  ----------------------------- MatrixCommitment -----------------------------

/// A MatrixCommitment provides a set of acceptable matrix characteristics.
/// Since we define the characteristics each with its own bit, a commitment
/// is implemented as a set of masks with bits set corresponding to the
/// acceptable characteristics.

//  ----------------------------------------------------------------------------
class SimTK_SimTKCOMMON_EXPORT MatrixCommitment {
public:
    MatrixCommitment() {} // set commitments to "none" and masks to "uncommitted"

    /// This is an implicit conversion from a MatrixStructure specification to 
    /// a MatrixCommitment with storage, outline, and condition uncommitted.
    MatrixCommitment(const MatrixStructure& str)
    { new (this) MatrixCommitment(str, MatrixStorage(), MatrixOutline(), MatrixCondition());}

    class Vector;
    class RowVector;
    class Triangular;
    class Symmetric;
    class Hermitian;
    class SkewSymmetric;
    class SkewHermitian;

    MatrixCommitment& commitSize(int m, int n) 
    {   commitNumRows(m); commitNumCols(n); return *this; }
    MatrixCommitment& commitNumRows(int m) 
    {   SimTK_SIZECHECK_NONNEG(m, "MatrixCommitment::commitNumRows()");
        masks.nr = m; return *this; }
    MatrixCommitment& commitNumCols(int n)  
    {   SimTK_SIZECHECK_NONNEG(n, "MatrixCommitment::commitNumCols()");
        masks.nc = n; return *this; }

    MatrixCommitment& commitBandwidth(int lb, int ub) 
    {  commitLowerBandwidth(lb); commitUpperBandwidth(ub); return *this;}
    MatrixCommitment& commitLowerBandwidth(int lb)
    {   SimTK_SIZECHECK_NONNEG(lb, "MatrixCommitment::commitLowerBandwidth()");
        masks.lband = lb; return *this; }
    MatrixCommitment& commitUpperBandwidth(int ub)
    {   SimTK_SIZECHECK_NONNEG(ub, "MatrixCommitment::commitUpperBandwidth()");
        masks.uband = ub; return *this; }

    MatrixCommitment& commitStructure(const MatrixStructure& s) 
    {   structure=s; masks.structure=s.mask(); return *this; }
    MatrixCommitment& commitStorage  (const MatrixStorage&   s) 
    {   storage=s;   masks.storage  =s.mask(); return *this; }
    MatrixCommitment& commitOutline  (const MatrixOutline&   o) 
    {   outline=o;   masks.outline  =o.mask(); return *this; }
    MatrixCommitment& commitCondition(const MatrixCondition& c) 
    {   condition=c; masks.condition=c.mask(); return *this; }

    /// For any handle commitment, we can calculate a "best character" for
    /// an allocation that satisfies the commitment, optionally with an 
    /// initial allocation size. Typically it is the least-restrictive 
    /// actual character that satisfies the commitment. For
    /// example, if the commitment is Triangular we'll allocate a full triangular
    /// matrix, not a banded one, or a symmetric one, or for that matter an
    /// identity matrix, all of which would satisfy the commitment. 
    /// The supplied sizes are used as minima -- if the commitment requires
    /// a larger minimum size you'll get that. For example, if you specify
    /// 0x0 but you're committed to a Column outline, you'll get 0x1.
    MatrixCharacter calcDefaultCharacter(int minNumRows, int minNumCols) const;

    /// These report the commitment as it was specified.
    const MatrixStructure&  getStructureCommitment() const {return structure;}
    const MatrixStorage&    getStorageCommitment()   const {return storage;}
    const MatrixOutline&    getOutlineCommitment()   const {return outline;}
    const MatrixCondition&  getConditionCommitment() const {return condition;}

    /// These report the masks of acceptable values generated from the commitment.
    const MatrixStructure::Mask&   getStructureMask() const {return masks.structure;}
    const MatrixStorage::Mask&     getStorageMask()   const {return masks.storage;}
    const MatrixOutline::Mask&     getOutlineMask()   const {return masks.outline;}
    const MatrixCondition::Mask&   getConditionMask() const {return masks.condition;}

    MatrixCharacter::Mask::SizeMask getNumRowsMask() const {return masks.nr;}
    MatrixCharacter::Mask::SizeMask getNumColsMask() const {return masks.nc;}
    MatrixCharacter::Mask::SizeMask getLowerBandwidthMask() const {return masks.lband;}
    MatrixCharacter::Mask::SizeMask getUpperBandwidthMask() const {return masks.uband;}

    int getDefaultNumRows() const {return masks.getDefaultNumRows();}
    int getDefaultNumCols() const {return masks.getDefaultNumRows();}

    bool isSizeOK(int m, int n) const {return masks.isSizeOK(m,n);} 
    bool isSizeOK(const std::pair<int,int>& mn) const
    {   return isSizeOK(mn.first, mn.second); }

    bool isBandwidthOK(int lower, int upper) const {return masks.isBandwidthOK(lower,upper);} 

    bool isSatisfiedBy(const MatrixCharacter& actual) const 
    {   return masks.isSatisfiedBy(actual); }
    bool isStructureOK(const MatrixStructure& s) const
    {   return getStructureMask().isSatisfiedBy(s); }
    bool isStorageOK(const MatrixStorage& s) const
    {   return getStorageMask().isSatisfiedBy(s); }
    bool isOutlineOK(const MatrixOutline& o) const
    {   return getOutlineMask().isSatisfiedBy(o); }
    bool isConditionOK(const MatrixCondition& c) const
    {   return getConditionMask().isSatisfiedBy(c); }

    bool isResizeable()      const {return masks.isResizeable();}
    bool isFullyResizeable() const {return masks.isFullyResizeable();;}
    bool isNumRowsLocked()  const {return masks.isNumRowsLocked();}
    bool isNumColsLocked()  const {return masks.isNumColsLocked();}

    bool isStructureCommitted() const 
    {   return !getStructureMask().isUncommitted(); }
    bool isStorageCommitted()   const 
    {   return !getStorageMask().isUncommitted();}
    bool isOutlineCommitted()   const 
    {   return !getOutlineMask().isUncommitted(); }
    bool isConditionCommitted() const 
    {   return !getConditionMask().isUncommitted();}

    /// Set commitment s to "none" and masks to "uncommitted" for all characteristics.
    void clear() {
        structure.setToNone(); 
        storage.setToNone(); 
        outline.setToNone(); 
        condition.setToNone();
        masks.setToUncommitted();
    }

protected:
    MatrixCommitment(const MatrixStructure& structure,
                     const MatrixStorage&   storage,                    
                     const MatrixOutline&   outline,
                     const MatrixCondition& condition)
    :   structure(structure), storage(storage), 
        outline(outline), condition(condition),
        masks() // set to all 1's 
    {
        if (outline.getOutline()==MatrixOutline::Scalar) commitSize(1,1);
        else if (outline.getOutline()==MatrixOutline::Column) commitNumCols(1);
        else if (outline.getOutline()==MatrixOutline::Row) commitNumRows(1);

        masks.structure = structure.mask();
        masks.storage   = storage.mask();
        masks.outline   = outline.mask();
        masks.condition = condition.mask(); 
    }
    
    /// These are the commitments as specified. They are used to fill in
    /// the corresponding bitmasks of acceptable characteristics below.
    MatrixStructure         structure;
    MatrixStorage           storage;
    MatrixOutline           outline;
    MatrixCondition         condition;

    /// These are the bitmasks of acceptable characteristics which
    /// would satisfy the above-specified commitments.
    MatrixCharacter::Mask   masks;
};


/// This is the default commitment for a column vector.
class MatrixCommitment::Vector : public MatrixCommitment {
public:
    /// Commit to a resizeable column vector.
    Vector()
    :   MatrixCommitment
        (   MatrixStructure(MatrixStructure::Matrix1d), 
            MatrixStorage(),
            MatrixOutline(MatrixOutline::Column), 
            MatrixCondition())
    {
    }
    /// Commit to a column vector of a particular length.
    explicit Vector(int m)
    :   MatrixCommitment
        (   MatrixStructure(MatrixStructure::Matrix1d), 
            MatrixStorage(),
            MatrixOutline(MatrixOutline::Column), 
            MatrixCondition())
    {
        commitNumRows(m);
    }
};

/// This is the default commitment for a row vector.
class MatrixCommitment::RowVector : public MatrixCommitment {
public:
    /// Commit to a resizeable row vector.
    RowVector()
    :   MatrixCommitment
        (   MatrixStructure(MatrixStructure::Matrix1d), 
            MatrixStorage(),
            MatrixOutline(MatrixOutline::Row), 
            MatrixCondition())
    {
    }
    /// Commit to a row vector of a particular length.
    explicit RowVector(int n)
    :   MatrixCommitment
        (   MatrixStructure(MatrixStructure::Matrix1d), 
            MatrixStorage(),
            MatrixOutline(MatrixOutline::Row), 
            MatrixCondition())
    {
        commitNumCols(n);
    }
};

/// This is the default commitment for a triangular matrix.
class MatrixCommitment::Triangular : public MatrixCommitment {
public:
    Triangular()
    :   MatrixCommitment(MatrixStructure::Triangular, MatrixStorage(),
                         MatrixOutline(), MatrixCondition())
    {
    }
};

/// This is the default commitment for a symmetric (*not* Hermitian) matrix.
/// There are no restrictions on the underlying elements or diagonals.
class MatrixCommitment::Symmetric : public MatrixCommitment {
public:
    Symmetric()
    :   MatrixCommitment(MatrixStructure::Symmetric, MatrixStorage(),
                         MatrixOutline(), MatrixCondition())
    {
    }
};

/// This is the default commitment for a Hermitian (*not* symmetric) matrix.
/// Diagonal elements must be real since they have to serve as their own
/// conjugates.
class MatrixCommitment::Hermitian : public MatrixCommitment {
public:
    Hermitian()
    :   MatrixCommitment
        (   MatrixStructure::Hermitian, 
            MatrixStorage(),
            MatrixOutline(), 
            MatrixCondition().setDiagonal(MatrixCondition::RealDiagonal))
    {
    }
};

/// This is the default commitment for skew symmetric (*not* skew Hermitian) 
/// matrix. Diagonal elements must be all-zero since they have to be their
/// own negation. Otherwise any elements are acceptable.
class MatrixCommitment::SkewSymmetric : public MatrixCommitment {
public:
    SkewSymmetric()
    :   MatrixCommitment
        (   MatrixStructure::SkewSymmetric, 
            MatrixStorage(),
            MatrixOutline(), 
            MatrixCondition().setDiagonal(MatrixCondition::ZeroDiagonal))
    {
    }
};

/// This is the default commitment for a skew Hermitian (*not* skew symmetric) 
/// matrix. Diagonal elements must be pure imaginary since when conjugated they
/// must be their own negation for a skew matrix.
class MatrixCommitment::SkewHermitian : public MatrixCommitment {
public:
    SkewHermitian()
    :   MatrixCommitment
        (   MatrixStructure::SkewHermitian, 
            MatrixStorage(),
            MatrixOutline(), 
            MatrixCondition().setDiagonal(MatrixCondition::ImaginaryDiagonal))
    {
    }
};

/// Output a textual description of a MatrixCommitment; handy for debugging.
/// @relates SimTK::MatrixCommitment
SimTK_SimTKCOMMON_EXPORT std::ostream& 
operator<<(std::ostream& o, const MatrixCommitment&);
     
} //namespace SimTK

#endif // SimTK_SIMMATRIX_MATRIX_CHARACTERISTICS_H_
