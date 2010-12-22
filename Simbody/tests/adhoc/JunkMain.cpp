/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-7 Stanford University and the Authors.         *
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

/**@file
 * This is an outer block for simulating ??? in various ways with Simbody.
 * This is about testing Simbody, *not* studying ???!
 */

#include "SimTKsimbody.h"

#include <string>
#include <iostream>
#include <exception>
#include <cmath>
#include <vector>
#include <map>
#include <set>
using std::cout;
using std::endl;

int main() { return 0; }
#ifdef NOTDEF

using namespace SimTK;

static const Real Deg2Rad = (Real)SimTK_DEGREE_TO_RADIAN;

typedef int BioType;
static const BioType InvalidBioType = -1;
static const BioType EthaneC = 1;
static const BioType EthaneH = 2;

static const BioType MethaneC = 3;
static const BioType MethaneH = 4;

typedef int Element;
static const Element InvalidElement = -1;
static const Element Carbon = 6;
static const Element Hydrogen = 1;

// These are the Ids in the shared atom/bond pool.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(AtomIndex)
SimTK_DEFINE_UNIQUE_INDEX_TYPE(BondIndex)

// These are the Ids local to a Compound.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(CompoundAtomIndex)
SimTK_DEFINE_UNIQUE_INDEX_TYPE(CompoundBondIndex)
SimTK_DEFINE_UNIQUE_INDEX_TYPE(SubcompoundIndex)

typedef String AtomName;
typedef String BondName;
typedef String CompoundName;

class Compound;
class Atom {
public:
    Atom() : element(InvalidElement), biotype(InvalidBioType), location(NaN) { }

    Atom(String nm, Element elt, BioType bt, const Vec3& loc)
        : name(nm), element(elt), biotype(bt), location(loc)
    { 
    }

    Atom& setElement(Element);
    Atom& setBioType(BioType);
    Atom& setName(AtomName);
    Atom& setPDBName(AtomName);

    bool isInCompound() const {return myCompound != 0 && myId.isValid();}

    // If an Atom is owned by a Compound these won't blow up.
    AtomIndex getAtomIndex() const {assert(isInCompound()); return myId;}
    const Compound& getCompound() const {assert(isInCompound()); return *myCompound;}
    Compound&       updCompound() const {assert(isInCompound()); return *myCompound;}
private:
    friend class Compound;
    Compound* myCompound;
    AtomIndex    myId;

    String    name;
    Element   element;
    BioType   biotype;
    Vec3      location;

    // Only valid if containing Compound has had its topology realized.
    std::vector<AtomIndex> bond12; // directly bonded atoms

    void setCompound(Compound& c, AtomIndex id) {
        assert(!myCompound && !myId.isValid());
        myCompound = &c;
        myId = id;
    }
};

class Bond {
public:
    // Bonds are not directional; we'll sort these by AtomIndex.
    Bond(AtomIndex a, AtomIndex b) 
      : lowId(std::min(a,b)), highId(std::max(a,b))
    {
        assert(a != b);
    }

    AtomIndex getAtom1() const {return lowId;}
    AtomIndex getAtom2() const {return highId;}

    bool isInCompound() const {return myCompound != 0 && myId.isValid();}

    // If an Bond is owned by a Compound these won't blow up.
    BondIndex getBondIndex() const {assert(isInCompound()); return myId;}
    const Compound& getCompound() const {assert(isInCompound()); return *myCompound;}
    Compound&       updCompound() const {assert(isInCompound()); return *myCompound;}
private:
    friend class Compound;
    Compound* myCompound;
    BondIndex    myId;

    AtomIndex lowId;
    AtomIndex highId;

    void setCompound(Compound& c, BondIndex id) {
        assert(!myCompound && !myId.isValid());
        myCompound = &c;
        myId = id;
    }
};

// This class contains the atoms and bonds for the physical objects
// we're building, as well as room for annotations
// that can be used when we produce a computational model of them. Annotations
// include Tinker biotype classification and PDB symbol for atoms, default
// spatial locations, and how to chop the objects up into rigid bodies, mobilizers
// and constraints for modeling with Simbody.
//
// A bond connects two distinct atoms from the pool, and only a single bond
// may exist between a given pair of atoms. The pool can be partioned
// into disjoint sets of atoms and bonds, where the
// atoms are all connected directly or indirectly by the bonds. Depending on 
// context these are called "chains" or "molecules". Sets with just one atom
// and no bonds are "free atoms".
// 
// Only the top level Compound or System has ownership of the AtomPool. All subcompounds
// reference the same pool.
class AtomPool {
public:
    AtomPool() : postConstructionAnalysisDone(false) { }

    Array<Atom> atoms;
    Array<Bond> bonds;

    // post-construction analysis data structures go here
    bool postConstructionAnalysisDone;
    Array<Vec3> configuration; // one per atom
};

// These are internal to Compound, indexed by CompoundAtomIndex, and refer to AtomIndexs in the pool.
class AtomInfo {
public:
    AtomInfo(const String& nm, AtomIndex id, const Vec3& loc)
      : name(nm), poolAtomIndex(id), location(loc)
    {
    }

    AtomIndex  poolAtomIndex;
    String  name;     // default name -- might be synonyms too
    Vec3    location; // in compound frame C
};

// These are indexed by CompoundBondIndex and refer to BondIndexs in the pool.
class BondInfo {
public:
    BondInfo(const String& nm, BondIndex id)
        : name(nm), poolBondIndex(id)
    { 
    }
        
    BondIndex  poolBondIndex;
    String  name;     // default name -- might be synonyms too
};

// These are indexed by SubcompoundIndex and actually *contain* the
// subcompounds owned by this parent compound.
class SubcompoundInfo {
public:
    SubcompoundInfo(CompoundRep* rep, const String& nm, const Transform& xform_PC)
        : subcompound(rep), name(nm), X_PC(xform_PC)
    { 
        rep->myHandle = &subcompound;
    }

    Compound  subcompound;
    String    name; // default name -- might be synonyms too
    Transform X_PC; // in parent compound frame P
};

// An AtomPool represents the unstructured, physical objects in the molecular
// simulation: atoms and bonds. We'll call these "pool atoms" and "pool bonds"
// to distinguish them from similar Compound objects. The top level Compound in
// the Compound tree owns the AtomPool; all sub-Compounds in a tree inherit
// the AtomPool from their parent compound.
//
//    pool atom
//    - element (or "deleted")
//    - valence
//    - formal charge
//    - list of pool bonds referencing this pool atom
//    - pathname of the Atom compound which owns this pool atom
//    - misc atom properties (e.g., pdb name)
//    - default position in the Ground frame
//
//    pool bond
//    - "deleted" flag
//    - pool atoms 1 & 2
//    - pathname of the Bond compound which owns this pool bond
//    - misc bond properties (e.g., rotatable?)
// 
// The ordering of pool atoms and bonds is insignificant, however their
// indices are stable even if atoms and bonds are deleted.
//
// -------------------
//
// A Compound is a structured interpretation of a set of atoms and bonds. For example, 
// a Compound can be used to express the idea that
// four atoms consisting of a carbon bonded to 3 hydrogens constitutes
// a "methyl" with atom names C, H1, H2, and H3, and to give default geometry
// for this structure. Compounds can recursively contain
// sub-Compounds.
//
// Positioning
// -----------
// Compounds are positioned relative to one another. In the absence of bonds, compounds
// are positioned relative to their parent compound.
// 
//
// Every Compound provides a unique "compound reference frame" C used for positioning
// geometric properties of the compound (including subcompounds), and a property
// list for non-geometric information.
//
// Compound "terminals" are the basic compounds from which all other compounds
// are built up. The four terminals are: Atoms, BondCenters, Bonds, and Properties.
//
// ATOM
//
// An Atom is a Compound consisting of
//   (1) an element type, valence (# bonds) nb, and formal charge
//   (2) an ordered list of nb BondCenters with orientations relative to 
//         the atom (compound) reference frame C, but with origins identical
//         to the atom frame's origin
//   (3) names for each of the BondCenters
//
// BOND
//
// A Bond consists of 
//   (1) an ordered pair of references to BondCenters, one "fixed" (inboard)
//         and one "moving" (outboard). If both are inboard, this bond 
//         closes a loop
//   (2) local names for those bond centers
//
//   A Compound can have at most one inboard bond; a rotation about that bond
//   repositions the Compound but otherwise it is fixed with respect to its 
//   parent Compound. If a Compound has *no* inboard bond, it is the base
//   body of some molecule and can be repositioned freely.
//
// BOND CENTER
//
// A BondCenter is a Compound which can exist only as a subcompound of an Atom (although
// higher-level compounds can contain named *references* to BondCenters).
// It contains the following elements:
//   (1) A reference to its parent Atom.
//   (2) Like any Compound, a BondCenter has a reference frame C. But here C's
//         origin must be coincident with the parent Atom's origin. C's x axis
//         gives the bond direction with +x pointing away from the parent Atom
//         towards its bonded partner. See below for more detail.
//   (3) Properties, such as whether it can be considered rotatable for
//         positioning purposes, whether it can take only a finite number
//         of orientations, or whether it is rigid.
//   (4) "Bonded" status: has this BondCenter been connected to another one?
//       - if so, by which Bond
//       - is this the "fixed" (inboard) side of the bond, or the "moving" 
//         (outboard) side?
//
//   When two BondCenters are connected by a bond, the x axes are made antiparallel
//   while the y axes are aligned, and the origins are separated along x by
//   a given bond length. In the case of a rotatable bond, a default angle
//   can be given which will result in the "moving" bond center's y axis being
//   rotated counterclockwise about the "fixed" bond center's x axis, so that
//   the y axes are separated by the given angle.
//
// PROPERTY
//
//   A Property is a Compound consisting only of a value. Its parent compound gives
//   it a name and a reference frame. The type can be determined at run time.
//   Interesting value types might include int, real, string, direction, point,
//   frame, enumerations, etc. and can be added as needed.
//
// GENERAL COMPOUND
//
// A general Compound consists of
//   (1) a reference frame, called the Compound Frame C,
//   (2) a reference to its parent Compound (unless it is the top level Compound)
//   (3) zero or more sub-Compounds (sharing the same atom pool as the parent), 
//         positioned with a Transform X_CS from C to subcompound frame S. Subcompounds
//         can overlap, and their atoms may or may not be bonded to other atoms
//         of the compound.
//   (4) Local names for any of the subcompounds or sub-subcompounds, including atoms, bonds,
//       bond centers, and properties.
//
// The list of subcompounds is maintained in the order they
// were added to the Compound. Names are local to the Compound, that is, different
// compounds will have different names for the same objects. There may be multiple
// names for the same object, even within a single Compound. Subcompounds form a
// tree structure and their contained objects can be named using a pathname-like
// naming scheme.
//
// TODO: not sure if we need a special "Reference" subcompound for renaming
// sub-sub-sub-...-compounds. Maybe that is just a Property with value type
// "compound"? 
//
//

class CompoundRep {
public:
    CompoundRep() : myHandle(0), atomPool(0), isAtomPoolOwner(true)
    {
    }

    virtual ~CompoundRep() {myHandle=0; atomPool=0;}

    CompoundRep* clone() {
        CompoundRep* p = cloneImpl();
        myHandle=0;
    }

    virtual CompoundRep* cloneImpl() const {assert(false);} // TODO

    CompoundRep(Compound& parent, const String& nameInParent, const String& compoundType);
    CompoundRep(AtomPool&, const String& compoundType);
    explicit CompoundRep(const String& compoundType); // allocates its own pool

    CompoundAtomIndex getCompoundAtomIndex(AtomName key) const {
        std::map<AtomName,CompoundAtomIndex>::const_iterator abn = atomByName.find(key);
        return abn==atomByName.end() ? InvalidCompoundAtomIndex : abn->second;
    }

    CompoundBondIndex getCompoundBondIndex(BondName key) const {
        std::map<BondName,CompoundBondIndex>::const_iterator bbn = bondByName.find(key);
        return bbn==bondByName.end() ? InvalidCompoundBondIndex : bbn->second;
    }

    SubcompoundIndex getSubcompoundIndex(CompoundName key) const {
        std::map<CompoundName,CompoundId>::const_iterator sbn = subcompoundByName.find(key);
        return sbn==subcompoundByName.end() ? InvalidSubcompoundIndex : sbn->second;
    }

    void deleteAtom(CompoundAtomIndex id) {
        assert(0 <= id && id < atomInfo.size());
        atomPool->deleteAtom(atomInfo[id].poolAtomIndex);
    }
    void deleteBond(CompoundBondIndex);
    void deleteSubcompound(SubcompoundIndex);

    void deleteAtom(AtomName);
    void deleteBond(BondName);
    void deleteCompound(CompoundName);

private:
    friend class Compound;
    Compound* myHandle;

    // Back pointer for traversing the subcompound tree upwards.
    Compound*     parent;     // null if top level Compound
    SubcompoundIndex idInParent; // index into the parent's subcompounds array, if parent!=0

    bool isAtomPoolOwner;
    AtomPool* atomPool;

    String compoundType;

    // Here is the Compound-local data.
    Array<AtomInfo>        atomInfo;       // index by CompoundAtomIndex
    Array<BondInfo>        bondInfo;       // index by CompoundBondIndex
    Array<SubcompoundInfo> subcompounds;   // index by SubcompoundIndex; includes the sub-Compounds

    // This is our set of names for the local set of atoms and bonds, which refer to
    // atoms and bonds in the atom/bond pool, and the set of local subcompounds. 
    // There can be more than one name for an atom, bond, or subcompound.
    std::map<AtomName,     CompoundAtomIndex> atomByName;
    std::map<BondName,     CompoundBondIndex> bondByName;
    std::map<CompoundName, SubcompoundIndex>  subcompoundByName;
};

class Compound {
public:
    Compound() : rep(0) { }

    Compound(const Compound&);
    Compound& operator=(const Compound&);
    ~Compound() {delete rep; rep=0;}

    Compound& setCompoundType(String);
    const String& getCompoundType() const {return getRep().compoundType;}

    bool hasParentCompound() const {return getRep().parent != 0;}
    const Compound& getParentCompound() const {assert(hasParent()); return *getRep().parent;}
    Compound&       updParentCompound() {assert(hasParent()); return *updRep().parent;}

    bool isAtomPoolOwner() const {return getRep().isAtomPoolOwner;}
    const AtomPool& getAtomPool() const {assert(getRep().atomPool); return *getRep().atomPool;}
    AtomPool&       updAtomPool() {assert(getRep().atomPool); return *updRep().atomPool;}

    // This takes over ownership of the Compound::Rep if the original handle was the
    // owner. If not, an exception is thrown.
    SubcompoundIndex adoptSubcompound(Compound& c, const String& newName, const Transform& X_PC) {
        assert(c.isOwnerHandle());
        assert(usesSameAtomPool(c) && !c.isAtomPoolOwner() && !c.hasParentCompound());
        updRep().subcompounds.push_back(Compound()); // add an empty handle to reserve id
        updRep().subcompounds.back().rep = c.rep;
        c.updRep().myHandle = &updRep().subcompounds.back(); // change ownership
        c.updRep().parent = this;
        c.updRep().idInParent = SubcompoundIndex(subcompounds.size()-1);
    }

    // TODO
    SubcompoundIndex copyInSubcompound(const Compound&, const String& newName, const Transform& X_PC);

    // Add a new atom to the atom pool and a corresponding compound atom here.
    // We are given a "base" frame B attached to compound C, and then we position
    // the new atom relative to B.
    CompoundAtomIndex addAtom(const String& atomName, const Transform& X_CB, const Vec3& loc, Element, BioType=InvalidBioType);
    CompoundBondIndex addNamedBond(const String& bondName, const String& atom1, const String& atom2);
    CompoundBondIndex addBond(const String& atom1, const String& atom2);

    // Add a new atom to the atom pool and to this compound, bond the atom to an existing atom with
    // the bond direction being in the +x direction of the given coordinate frame. The existing atom
    // must have a location relative to the current compound.
    CompoundAtomIndex bondAtom(const String& newAtomName, const String& existingAtom, const Rotation& R_CB,
                            Real bondLength, Element, BioType=InvalidBioType)
    {
        addAtom(newAtomName, Transform(R_CB, getAtomLocationInCompound(existingAtom)), Vec3(bondLength,0,0),
                Element, BioType);
        addBond(existingAtom, newAtomName);
    }

    bool isOwnerHandle() const {
        return getRep().myHandle == this;
    }

    class Rep;
private:
    Compound::Rep* rep;
    const Compound::Rep& getRep() const {assert(rep); return *rep;}
    CompoundRep& updRep()               {assert(rep); return *rep;}
};

// A BondCenter is a small compound consisting of just its compound frame
// and a single atom located at the origin of that frame. The frame's x-axis
// points in the direction of the currently unknown "other" atom to which
// this one will eventually be bonded.
class BondCenter : public Compound {
public:
    BondCenter(Compound& parent, String bondCenterNameInParent,
               String bondAtom, const Rotation& R_PB)
        : Compound("bondCenter")
    {
        atomInfo.push_back(AtomInfo("bondAtom", parent.getAtomIndex(bondAtom), Vec3(0));
        parent.adoptSubcompound(*this, bondCenterNameInParent);
    }

private:
};

// CH4
//
//
//          H2     
//           \   y
//            \  |
//             . C --> x H1
//      (H4) .  /      
//         *   z    
//       H3
//
// We'll put the frame origin at the C, H1 is along x, H2 in the x-y plane, H3 a positive
// rotation around x by 120 degrees from H2, H4 a further 120 degrees.
//
class Methane : public Compound {
public:
    Methane(Compound& parent, String compoundNameInParent, const Transform& X_GM=Transform())
      : Compound(parent.atomPool, "methane")
    {
        construct();
        parent.adoptCompound(*this, compoundNameInParent, X_GM);
    }
    Methane() : Compound("methane") 
    {
        construct();
    }

    // Assign biotypes since methane is a complete molecule.
    void construct() {
        addAtom ("C",  Vec3(0), Carbon, MethaneC); // Cartesian location in frame C
        bondAtom("H1", "C", Rotation(), 0.1102, Hydrogen, MethaneH);
        for (int i=0; i<3; ++i) 
            bondAtom("H"+String(i+2), "C", 
                     Rotation::aboutZThenOldX(109.47*Deg2Rad, i*120*Deg2Rad),
                     0.1102, Hydrogen, MethaneH);
    }
};

class Methyl : public Compound {
public:
    Methyl(Compound& parent, String compoundNameInParent, const Transform& X_GM=Transform())
      : Compound(parent.atomPool, "methyl")
    {
        construct();
        parent.adoptCompound(*this, compoundNameInParent, X_GM);
    }
    Methyl() : Compound("methyl") 
    {
        construct();
    }

    // Assign tentative biotypes. These will likely change if this is incorporated
    // in a larger compound.
    void construct() {
        addAtom ("C",  Vec3(0), Carbon, MethylC); // Cartesian location in frame C
        for (int i=0; i<3; ++i) 
            bondAtom("H"+String(i+1), "C", 
                     Rotation::aboutZThenOldX(109.47*Deg2Rad, i*120*Deg2Rad),
                     0.1112, Hydrogen, MethylH);

        // A "reaction center" marks a place on this compound where it 
        // expects a bond to be formed to an as-yet-to-be-defined atom.
        // A reaction center is a small subcompound consisting of a frame and a 
        // single atom at its origin. The x-axis points in the direction
        // of the bond to be formed, that is from the reaction atom to its
        // bonded partner.
        BondCenter(*this, "methylate", "C", Rotation());
    }
};

class Ethane : public Compound {
public:
    Ethane() : Compound("ethane")
    {
        // Place two methyls here and give them unique compound names.
        Methyl m1(*this, "methyl1"); // carbon goes at Ethane's origin
        Methyl m2(*this, "methyl2"); // carbon goes at Ethane's origin

        // Forming this bond will cause methyl2 to reorient so that the
        // reaction center x axes are antiparallel, via rotation about y.
        connectBondCenters("CC", "m1/methylate", "m2/methylate", .15247);
        setBondIsRotatable("CC", true);

        // Now fill in the Ethane name table. This doesn't affect the
        // two methane name tables.
        setAtomName("methyl1/C",  "C1"); setAtomName("methyl2/C",  "C2");
        setAtomName("methyl1/H1", "H1"); setAtomName("methyl2/H1", "H4");
        setAtomName("methyl1/H2", "H2"); setAtomName("methyl2/H2", "H5");
        setAtomName("methyl1/H3", "H3"); setAtomName("methyl2/H3", "H6");

        setBioType("C1", EthaneC); setBioType("C2", EthaneC);
        for (int i=1; i<=6; ++i)
            setBioType("H" + String(i), EthaneH);
    }
};

class PolyAlanine : public Compound {
public:
    explicit PolyAlanine(int n) : Compound("polyalanine"), nResidues(n)
    {
        Compound prev = Alanine(*this, "A0");
        for (int i=1; i< nResidues; ++i) {
            Alanine alai(*this, "A"+String(i));
            addPeptideBond(prev.getCompound("carboxyl"), alai.getCompound("amino"));
        }
    }

private:
};

class TetrahedralCarbon : public Compound {
public:
    TetrahedralCarbon(Compound&, String atomName, Transform X_PC) {
        BondCenter(*this, "bond1", "C", Rotation());
        BondCenter(*this, "bond2", "C", Rotation(109.5*Deg2Rad,ZAxis));
        BondCenter(*this, "bond3", "C", Rotation::aboutZThenOldX(109.5*Deg2Rad,  120*Deg2Rad));
        BondCenter(*this, "bond4", "C", Rotation::aboutZThenOldX(109.5*Deg2Rad, -120*Deg2Rad));
    }
};

class Compound::BondCenter : public Compound {
public:
    static const int Rotatable = 1;
    static const int Restricted = 2;
    static const int Rigid = 3;
    explicit BondCenter(int rot=Rotatable);
private:
    int allowedRotation;
    bool isBonded;
    Compound::Bond* whichBond;
};


class Compound::Atom : public Compound {
public:
    Atom(Element e, int valence, int formalCharge=0) 
      : element(elt), formalCharge(formalCharge), bondCenters(valence)
    {
    }

    Element getElement()  const {return element;}
    int getFormalCharge() const {return formalCharge;}
    int getValence()      const {return (int)bondCenters.size();}

    int getNumBondCenters() const {return (int)bondCenters.size();} 
    const BondCenter& getBondCenter(int i) const {
        assert(0 <= i && i < getNumBondCenters());
        assert(bondCenters[i].isValid());
        return BondCenter::downcast(getSubcompound(bondCenters[i]);
    }

    // Add a bond center for a particular bond index (0..valence-1).
    void setBondCenter(int i, const String& name, const Rotation& orientation, 
                       const BondCenter& bc)
    {
        assert(0 <= i && i < getNumBondCenters());
        assert(!bondCenter[i].isValid());
        bondCenter[i] = addSubcompound(name, Transform(orientation,Vec3(0)), bc);
    }

    // Add a new bond center at the first unassigned slot. The bond center index (0..valence-1)
    // is returned.
    int addBondCenter(const String& name, const Rotation& orientation, const BondCenter& bc) {
        // Grab first unassigned bond center.
        for (int i=0; i < (int)bondCenters.size(); ++i)
            if (!bondCenters[i].isValid()) {
                setBondCenter(i,name,orientation,bc);
                return i;
            }

        assert(!"Tried to add too many bond centers");
    }

private:
    Element element;
    int formalCharge;
    std::vector<SubcompoundIndex> bondCenters; // indices into subcompound list
};

class TetrahedralCarbon : public Compound::Atom {
public:
    TetrahedralCarbon() : Atom(Carbon, 4)
    {
        // Defines bond center subcompounds, with +x pointing from this Atom to where
        // its bonded partner will be.
        setBondCenter(0, "bond1", Rotation()); // aligned with compound x
        setBondCenter(1, "bond2", Rotation(109.5*Deg2Rad,ZAxis)); // in compound xy plane
        setBondCenter(2, "bond3", Rotation::aboutZThenOldX(109.5*Deg2Rad,  120*Deg2Rad)); // on +z side
        setBondCenter(3, "bond4", Rotation::aboutZThenOldX(109.5*Deg2Rad, -120*Deg2Rad)); // on -z side   
    }
};

class Methyl : public Compound {
public:
    Methyl() {
        addNewAtom ("C",  TetrahedralCarbon(), Transform());
        bondNewAtom("H1", UnivalentAtom(Hydrogen), "bond", "C/bond2", CHbondLength);
        bondNewAtom("H2", UnivalentAtom(Hydrogen), "bond", "C/bond3", CHbondLength);
        bondNewAtom("H3", UnivalentAtom(Hydrogen), "bond", "C/bond4", CHbondLength);
        addSubcompoundName("bond", "C/bond1"); // new name for existing bond center
    }
};


class Ethane : public Compound {
public:
    Ethane() {
        addNewSubcompound ("methyl1", Methyl(), Transform());
        bondNewSubcompound("methyl2", Methyl(), "bond", "methyl1/bond");

        //OR

        addNewCompound("methyl2", Methyl());
        addNewBond("torsion", "methyl1/bond", "methyl2/bond");
    }
};



// TODO: should be subclass of AminoAcid
class Serine : public Compound::AminoAcid {
public:
    Serine() {
        addNewAtom ("Ca", TetrahedralCarbon(), Transform());
        bondNewAtom("N",  TrigonalAtom(Nitrogen), "bond1", "Ca/bond1",  NCbondLength);
        bondNewAtom("C",  TetrahedralCarbon(),    "bond1", "Ca/bond2",  CCbondLength);
        bondNewAtom("Ha", HydrogenAtom(),         "bond",  "Ca/bond3",  CHbondLength);
        bondNewAtom("Cb", TetrahedralCarbon(),    "bond1", "Ca/bond4",  CCbondLength);

        // At phi==0, the H on the nitrogen is in the same plane as C-CA-N but
        // in "trans". (TODO: check this)
        bondNewAtom("Hn", HydrogenAtom(), 
                    zTransform("N", "Ca", "C", NHbondLength, CaNHangle, DefaultPhiAngle));

    }
};


int main() {
try
  { 
    
  }
catch (const std::exception& e)
  {
    printf("EXCEPTION THROWN: %s\n", e.what());
  }
catch (...)
  {
    printf("UNKNOWN EXCEPTION THROWN\n");
  }    return 0;
}


// For now we'll allow only letters, digits, and underscore in names. Case matters.
bool Compound::isLegalName(const String& n) {
    if (n.size()==0) return false;
    for (int i=0; i<n.size(); ++i)
        if (!(isalnum(n[i]) || n[i]=='_'))
            return false;
    return true;
}

// Take pathname of the form xxx/yyy/zzz, check its validity and optionally
// return as a list of separate subsystem names. We return true if we're successful,
// false if the pathname is malformed in some way. In that case the last segment
// returned will be the one that caused trouble.
bool Compound::isLegalSubmodelPathname(const String& pathname, 
                                       std::vector<String>* segments)
{
    String t;
    const int end = pathname.size();
    int nxt = 0;
    if (segments) segments->clear();
    bool foundAtLeastOne = false;
    // for each segment
    while (nxt < end) {
        // for each character of a segment
        while (nxt < end && pathname[nxt] != '/')
            t += pathname[nxt++];
        foundAtLeastOne = true;
        if (segments) segments->push_back(t);
        if (!isLegalSubmodelName(t))
            return false;
        t.clear();
        ++nxt; // skip '/' (or harmless extra increment at end)
    }
    return foundAtLeastOne;
}

#endif
