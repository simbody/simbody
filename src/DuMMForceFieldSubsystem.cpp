/* Portions copyright (c) 2006 Stanford University and Michael Sherman.
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


/**@file
 *
 * Private implementation of DuMMForceFieldSubsystem.
 */

#include "Simbody.h"
#include "simbody/internal/ForceSubsystem.h"
#include "simbody/internal/DuMMForceFieldSubsystem.h"

#include "simbody/internal/DecorativeGeometry.h"

#include "ForceSubsystemRep.h"

#include <string>
#include <vector>
#include <cmath>
#include <cassert>
#include <set>
#include <map>
#include <algorithm>

namespace SimTK {

static const Real Pi = std::acos(Real(-1));
static const Real RadiansPerDegree = Pi/180;

// Convert energy from Kcal/mol to consistent units Da-A^2/ps^2.
static const Real EnergyUnitsPerKcal = 418.4; // exact 

// This is Coulomb's constant 1/(4*pi*e0) in units which convert
// e^2/A to kcal/mol, followed by conversion to consistent energy units
// This constant was calculated (by both me and Jay Ponder) from the
// NIST physical constants at http://physics.nist.gov/constants
// (2002 CODATA).
static const Real CoulombFac = 332.06371 * EnergyUnitsPerKcal;

class IntPair {
public:
    IntPair() {ints[0]=ints[1]=-1;}
    IntPair(int i1, int i2, bool canon=false) {
        ints[0]=i1; ints[1]=i2;
        if (canon) canonicalize();
    }
    int operator[](int i) const {assert(0<=i&&i<2); return ints[i];}
    bool isValid() const {return ints[0]>=0 && ints[1]>=0;}
    // canonical is low,high
    void canonicalize() {if(ints[0]>ints[1]) std::swap(ints[0],ints[1]);}
private:
    int ints[2];
};
inline bool operator<(const IntPair& i1, const IntPair& i2) {
    assert(i1.isValid() && i2.isValid());
    if (i1[0] < i2[0]) return true;
    if (i1[0] > i2[0]) return false;
    return i1[1] < i2[1];
}

class IntTriple {
public:
    IntTriple() {ints[0]=ints[1]=ints[2]=-1;}
    IntTriple(int i1, int i2, int i3, bool canon=false) {
        ints[0]= i1; ints[1]=i2; ints[2]=i3;
        if (canon) canonicalize();
    }
    int operator[](int i) const {assert(0<=i&&i<3); return ints[i];}
    bool isValid() const {return ints[0]>=0 && ints[1]>=0 && ints[2]>=0;}
    // canonical has 1st number <= last number; middle stays put
    void canonicalize() {if(ints[0]>ints[2]) std::swap(ints[0],ints[2]);}
private:
    int ints[3];
};
inline bool operator<(const IntTriple& i1, const IntTriple& i2) {
    assert(i1.isValid() && i2.isValid());
    if (i1[0] < i2[0]) return true;
    if (i1[0] > i2[0]) return false;
    if (i1[1] < i2[1]) return true;
    if (i1[1] > i2[1]) return false;
    return i1[2] < i2[2];
}

class IntQuad {
public:
    IntQuad() {ints[0]=ints[1]=ints[2]=ints[3]=-1;}
    IntQuad(int i1, int i2, int i3, int i4, bool canon=false) {
        ints[0]= i1; ints[1]=i2; ints[2]=i3; ints[3]=i4;
        if (canon) canonicalize();
    }
    int operator[](int i) const {assert(0<=i&&i<4); return ints[i];}
    bool isValid() const {return ints[0]>=0 && ints[1]>=0 && ints[2]>=0 && ints[3]>=0;}
    // canonical has 1st number <= last number; middle two must swap
    // if the outside ones do
    void canonicalize() {
        if(ints[0]>ints[3]) {
            std::swap(ints[0],ints[3]); 
            std::swap(ints[1],ints[2]);
        }
    }

private:
    int ints[4];
};
inline bool operator<(const IntQuad& i1, const IntQuad& i2) {
    assert(i1.isValid() && i2.isValid());
    if (i1[0] < i2[0]) return true;
    if (i1[0] > i2[0]) return false;
    if (i1[1] < i2[1]) return true;
    if (i1[1] > i2[1]) return false;
    if (i1[2] < i2[2]) return true;
    if (i1[2] > i2[2]) return false;
    return i1[3] < i2[3];
}

// Vdw combining functions
// -----------------------
// There are several in common use. The most common
// one, Lorentz-Berthelot is also the worst one!
// The pragmatically best seems to be the Waldman-Hagler rule, which
// we will use by default. In between is the Halgren-HHG
// rule. Another good rule is Tang-Toennies but it requires
// additional empirical data (the "sixth dispersion coefficient"
// C6) which we do not have available. An alternative to
// Tang-Toennies is Kong, which uses the Tang-Toennies radius
// formula, but Waldman-Hagler's well depth formula (and Kong
// came considerably before either of them).
//
// The Lennard-Jones 12-6 potential is specified as follows:
// Each atom type i has two parameters ri and ei, resp. the
// van der Waals radius and energy well depth. The radii are
// defined so that if two atoms of type i are separated by
// a distance dmin=2*ri, then the van der Waals energy is -ei.
// For a pair of atoms of types i and j we define an effective
// separation dmin_ij and well depth e_ij. Then if the vector
// from atom i to atom j is v, and d=|v| we have
//
//    Evdw(d) = e_ij * ( (dmin_ij/d)^12 - 2*(dmin_ij/d)^6 )
//
//    Fvdw_j(d) = -grad_j(Evdw) 
//              = 12 e_ij * ( (dmin_ij/d)^12 - (dmin_ij/d)^6 ) * v/d^2
//    Fvdw_i(d) = -Fvdw_j(d)
//
// Some cautions: it is common among force fields to specify
// the vdw size (1) either by radius or diameter, and (2) by
// minimum energy or zero crossing. In the latter case the
// symbol "sigma" is used instead of "r", with r=2^(1/6) * sigma
// (that is, sigma is smaller than r). We will be using the
// "radius at minimum energy" convention; note that that has to
// be doubled to produce the dmin used in the LJ formula.



static inline Real arithmeticMean(Real a, Real b) {
    return 0.5*(a+b);
}
static inline Real geometricMean(Real a, Real b) {
    return std::sqrt(a*b);
}
static inline Real harmonicMean(Real a, Real b) {
    return (2*a*b) / (a+b);
}


// cubicMean = (a^3+b^3)/(a^2+b^2)
static inline Real cubicMean(Real a, Real b) {
    return (a*a*a+b*b*b)/(a*a+b*b);
}

// Harmonic mean of harmonic & geometric means
// hhgMean = 4ab/(sqrt(a)+sqrt(b))^2
static inline Real hhgMean(Real a, Real b) {
    return harmonicMean(harmonicMean(a,b), geometricMean(a,b));
}

// Used in AMBER, CHARMM, and MM2/3 (but MMs don't use LJ)
static inline void vdwCombineLorentzBerthelot(
    Real ri, Real rj, Real ei, Real ej,
    Real& r, Real& e)
{
    r = arithmeticMean(ri,rj);
    e = geometricMean(ei,ej);
}

// Used in OPLS, DANG
static inline void vdwCombineJorgensen(
    Real ri, Real rj, Real ei, Real ej,
    Real& r, Real& e)
{
    r = geometricMean(ri,rj);
    e = geometricMean(ei,ej);
}

// Used in MMFF, AMOEBA (but with Buffered 14-7 rather than LJ)
static inline void vdwCombineHalgrenHHG(
    Real ri, Real rj, Real ei, Real ej,
    Real& r, Real& e)
{
    r = cubicMean(ri,rj);
    e = hhgMean(ei,ej);
}

static const Real oo6  = Real(1)/Real(6);
static const Real oo13 = Real(1)/Real(13);

// This doesn't seem to be used by anyone but it should be!
// Ref: Waldman, M. & Hagler, A.T. New combining rules for
// rare gas van der Waals parameters. 
// J. Comput. Chem. 14(9):1077 (1973).
static inline void vdwCombineWaldmanHagler(
    Real ri, Real rj, Real ei, Real ej,
    Real& r, Real& e)
{
    const Real ri3 = ri*ri*ri, ri6 = ri3*ri3;
    const Real rj3 = rj*rj*rj, rj6 = rj3*rj3;
    const Real er6 = geometricMean(ei*ri6, ej*rj6);
    const Real r6  = arithmeticMean(ri6, rj6);

    r = std::pow(r6, oo6);
    e = er6 / r6;
}

// This is a possible alternative to Waldman-Hagler. It uses 
// the same well depth combination term as WH, but with a different
// radius combination term which is the same as Tang-Toennies.
// Ref: Kong, C.L. Combining rules for intermolecular potential
// parameters. II. Rules for the Lennard-Jones (12-6) potential
// and the Morse potential. J. Chem. Phys. 59(5):2464 (1973).
// Comparison with WH: Delhommelle, J. & Millie, P. Inadequacy of 
// the Lorentz-Berthelot combining rules for accurate predictions
// of equilibrium properties by molecular simulation. Molecular
// Physics 99(8):619 (2001).

static inline void vdwCombineKong(
    Real ri, Real rj, Real ei, Real ej,
    Real& r, Real& e)
{
    const Real ri3 = ri*ri*ri, ri6 = ri3*ri3, ri12 = ri6*ri6;
    const Real rj3 = rj*rj*rj, rj6 = rj3*rj3, rj12 = rj6*rj6;
    const Real er6 = geometricMean(ei*ri6, ej*rj6);

    // calculate (ei*ri^12)^(1/13), etc.
    const Real eri12_13 = std::pow(ei*ri12, oo13);
    const Real erj12_13 = std::pow(ej*rj12, oo13);
    const Real er12_13  = arithmeticMean(eri12_13, erj12_13);
    const Real r6 =  std::pow(er12_13, 13) / er6;

    r = std::pow(r6, oo6);
    e = er6 / r6;
}

class Element {
public:
    Element() : atomicNumber(-1), mass(-1), defaultColor(Gray) { }
    Element(int anum, const char* sym, const char* nm, Real m)
        : atomicNumber(anum), mass(m), defaultColor(Gray), symbol(sym), name(nm)
    {
        assert(isValid());
    }
    bool isValid() const {return atomicNumber > 0 && mass > 0;}

    Element& setDefaultColor(const Vec3& c) {
        defaultColor = c;
        return *this;
    }

    // These are all Topological state variables, that is,
    // set during construction and constant thereafter.
    int atomicNumber;
    Real mass;         // in Daltons (Da, g/mol)
    Vec3 defaultColor;
    std::string symbol;
    std::string name;
};

class AtomClass {
public:
    AtomClass() : element(-1), valence(-1), vdwRadius(-1), vdwWellDepth(-1) { }
    AtomClass(int id, const char* nm, int e, int v, Real rad, Real wellKcal)
      : atomClassId(id), name(nm), element(e), valence(v), vdwRadius(rad), 
        vdwWellDepth(wellKcal*EnergyUnitsPerKcal)
    { 
        assert(isValid());
    }
    bool isValid() const {return atomClassId >= 0 && element > 0 && valence >= 0 
                                 && vdwRadius >= 0 && vdwWellDepth >= 0;}

    void invalidateTopologicalCache() {
        vdwDij.clear();
        vdwEij.clear();
    }

    void dump() const {
        printf("   %d(%s): element=%d, valence=%d vdwRad=%g, vdwDepth(Kcal)=%g\n",
            atomClassId, name.c_str(), element, valence, vdwRadius, vdwWellDepth/EnergyUnitsPerKcal);
        printf("    vdwDij:");
        for (int i=0; i< (int)vdwDij.size(); ++i)
            printf(" %g", vdwDij[i]);
        printf("\n    vdwEij:");
        for (int i=0; i< (int)vdwEij.size(); ++i)
            printf(" %g", vdwEij[i]/EnergyUnitsPerKcal);
        printf("\n");
    }

        // TOPOLOGICAL STATE VARIABLES
        //   Filled in during construction.

    int         atomClassId;
    std::string name;

    int     element;
    int     valence;       // # of direct bonds expected
    Real    vdwRadius;     // ri, Angstroms
    Real    vdwWellDepth;  // ei, Da-A^2/ps^2


        // TOPOLOGICAL CACHE ENTRIES
        //   These are calculated in realizeConstruction() from topological
        //   state variables (from here or others in the DuMM class).

    // After all types have been defined, we can calculate vdw 
    // combining rules for dmin and well depth energy. We only fill
    // in entries for pairings of this class with itself and with
    // higher-numbered atom types, so to find the entry for class c, 
    // index these arrays by c-atomClassId where atomClassId is the
    // class Id of the present AtomClass.
    // Note that different combining rules may be used but they
    // will always result in a pair of vdw parameters.
    std::vector<Real> vdwDij;   // A
    std::vector<Real> vdwEij;   // Da-A^2/ps^2
};

class ChargedAtomType {
public:
    ChargedAtomType() : chargedAtomTypeId(-1), atomClass(-1), partialCharge(CNT<Real>::getNaN()) { }
    ChargedAtomType(int id, const char* nm, int aclass, Real chg)
      : chargedAtomTypeId(id), name(nm), atomClass(aclass), partialCharge(chg) 
    { 
        assert(isValid());
    }
    bool isValid() const {return chargedAtomTypeId >= 0 && atomClass >= 0;}

    void dump() const {
        printf("    %d(%s): atomClass=%d, chg=%g\n", 
               chargedAtomTypeId, name.c_str(), atomClass, partialCharge);
    }

    // These are all Topological state variables, filled in during construction.
    // There are no calculations to be performed.
    int         chargedAtomTypeId;
    std::string name;

    int         atomClass;
    Real        partialCharge; // qi, in e (charge on proton)

};

// This represents bond-stretch information for a pair of atom types.
// Use an IntPair as a key.
class BondStretch {
public:
    BondStretch() : k(-1), d0(-1) { }
    BondStretch(Real stiffnessKcalPerASq, Real length) 
      : k(stiffnessKcalPerASq*EnergyUnitsPerKcal), d0(length) { 
        assert(isValid());
    }
    bool isValid() const {return k >= 0 && d0 >= 0; }
    Real k;  // in energy units per A^2, i.e. Da/ps^2
    Real d0; // distance at which force is 0 (in A)
};

class BondBend {
public:
    BondBend() : k(-1), theta0(-1) { }
    BondBend(Real stiffnessKcalPerRadSq, Real angleDeg) 
      : k(stiffnessKcalPerRadSq*EnergyUnitsPerKcal), 
        theta0(angleDeg*RadiansPerDegree) {
        assert(isValid());
    }
    bool isValid() const {return k >= 0 && (0 <= theta0 && theta0 <= Pi);}

    // Given a central atom location c bonded to atoms at r and s,
    // calculate the angle between them, the potential energy,
    // and forces on each of the three atoms.
    void harmonic(const Vec3& cG, const Vec3& rG, const Vec3& sG,
                  Real& theta, Real& pe, Vec3& cf, Vec3& rf, Vec3& sf) const;

    Real k;      // energy units per rad^2, i.e. Da-A^2/(ps^2-rad^2)
    Real theta0; // unstressed angle in radians
};

//
// Torsion term for atoms bonded r-x-y-s. Rotation occurs about
// the axis v=y-x, that is, a vector from x to y. We define a torsion
// angle theta using the "polymer convention" rather than the IUPAC
// one which is 180 degrees different. Ours is like this:
//             r                         r      s
//   theta=0    \             theta=180   \    / 
//               x--y                      x--y
//                   \
//                    s
// The sign convention is the same for IUPAC and polymer:
// A positive angle is defined by considering r-x fixed in space. Then
// using the right and rule around v (that is, thumb points from x to y)
// a positive rotation rotates y->s in the direction of your fingers.
//
// We use a periodic energy function like this:
//       E(theta) = sum E_n(1 + cos(n*theta - theta0_n))
// where n is the periodicity, E_n is the amplitude (kcal/mol) for
// term n, and theta0_n is the phase offset for term n. The torque
// term (applied about the v axis) is then
//       T(theta) = -[sum -n*E_n*sin(n*theta - theta0_n)]
// We have to translate this into forces on the four atoms.
// 
class TorsionTerm {
public:
    TorsionTerm() : periodicity(-1), amplitude(-1), theta0(-1) { }
    TorsionTerm(int n, Real amp, Real th0) 
      : periodicity(n), amplitude(amp*EnergyUnitsPerKcal), theta0(th0*RadiansPerDegree) {
        assert(isValid());
    }
    bool isValid() const {return periodicity > 0 && amplitude >= 0 
                                 && -Pi < theta0 && theta0 <= Pi;}
    Real energy(Real theta) const {
        return amplitude*(1 + std::cos(periodicity*theta-theta0));
    }
    Real torque(Real theta) const {
        return periodicity*amplitude*std::sin(periodicity*theta-theta0);
    }

    int  periodicity; // 1=360, 2=180, 3=120, etc.
    Real amplitude; // energy units (Da-A^2/ps^2)
    Real theta0;    // radians
};

class BondTorsion {
public:
    BondTorsion() { }
    void addTerm(const TorsionTerm& tt) {
        assert(!hasTerm(tt.periodicity));
        terms.push_back(tt);
    }
    bool isValid() const {return !terms.empty();}
    bool hasTerm(int n) const {
        for (int i=0; i<(int)terms.size(); ++i)
            if (terms[i].periodicity == n) return true;
        return false;
    }

    // Given atom locations r-x-y-s in the ground frame, calculate the
    // torsion angle, energy and a force on each atom so that the desired
    // pure torque is produced.
    void periodic(const Vec3& rG, const Vec3& xG, const Vec3& yG, const Vec3& sG,
                  Real& theta, Real& pe, 
                  Vec3& rf, Vec3& xf, Vec3& yf, Vec3& sf) const;
    
    std::vector<TorsionTerm> terms;
};


class AtomPlacement {
public:
    AtomPlacement() : atomId(-1) { }
    AtomPlacement(int a, const Vec3& s) : atomId(a), station(s) {
        assert(isValid());
    }
    bool isValid() const {return atomId >= 0;}

    int  atomId;
    Vec3 station;
};
inline bool operator<(const AtomPlacement& a1, const AtomPlacement& a2) {
    return a1.atomId < a2.atomId;
}
inline bool operator==(const AtomPlacement& a1, const AtomPlacement& a2) {
    return a1.atomId == a2.atomId;
}

class RigidGroupPlacement {
public:
    RigidGroupPlacement() : rigidGroupId(-1) { }
    RigidGroupPlacement(int g, const Transform& t) : rigidGroupId(g), placement(t) {
        assert(isValid());
    }
    bool isValid() const {return rigidGroupId >= 0;}

    int         rigidGroupId;
    Transform   placement;
};
inline bool operator<(const RigidGroupPlacement& r1, const RigidGroupPlacement& r2) {
    return r1.rigidGroupId < r2.rigidGroupId;
}
inline bool operator==(const RigidGroupPlacement& r1, const RigidGroupPlacement& r2) {
    return r1.rigidGroupId == r2.rigidGroupId;
}

typedef std::vector<int>                    AtomArray;
typedef std::vector<AtomPlacement>          AtomPlacementArray;
typedef std::set<AtomPlacement>             AtomPlacementSet;
typedef std::set<RigidGroupPlacement>       RigidGroupPlacementSet;

class Atom {
public:
    Atom() 
      : atomId(-1), chargedAtomTypeId(-1), bodyId(-1) {
    }
    Atom(int t, int aId) : atomId(aId), chargedAtomTypeId(t), bodyId(-1) {
        assert(isValid());
    }
    bool isValid() const {return atomId>=0 && chargedAtomTypeId>=0;}

    bool isBondedTo(int anum) const {
        for (int i=0; i<(int)bond12.size(); ++i)
            if (bond12[i] == anum) return true;
        return false;
    }

    void dump() const;


    void invalidateTopologicalCache() {
        bodyId = -1; station_B.setToNaN();
        bond13.clear(); bond14.clear(); bond15.clear();
        xbond12.clear(); xbond13.clear(); xbond14.clear(); xbond15.clear();
        stretch.clear(); bend.clear(); torsion.clear();
    }

public:
        // TOPOLOGICAL STATE VARIABLES
        //   Filled in during construction.

    int         atomId;
    int         chargedAtomTypeId;
    AtomArray   bond12;

        // TOPOLOGICAL CACHE ENTRIES
        //   These are calculated in realizeConstruction() from topological
        //   state variables (from here or others in the DuMM class).

    // After the atom or its containing group has been attached to a
    // body, we can fill these in.
    int  bodyId;
    Vec3 station_B; // atom's station fixed in body bodyId's frame

    // This is a set of lists which identify atoms nearby in the
    // molecules bond structure. The bond12 list above contains the directly
    // bonded (1-2) atoms; the 13 list below has the 1-(2)-3 bonded atoms,
    // etc. The current Atom is always "1" so it isn't stored.

    std::vector<IntPair>   bond13;
    std::vector<IntTriple> bond14;
    std::vector<IntQuad>   bond15;

    // These are shorter versions of the bond lists in which only those
    // bonds which include atoms from at least two bodies are included.
    // Note that each bond will appear twice in the overall data structure,
    // in the Atom entries for the atoms at either end. We avoid double
    // processing by only processing the instance in which the first atoms's
    // ID is the lower of the two. But we need to keep both copies because
    // these are used for scaling nearby interaction during non-bonded 
    // calculation.
    std::vector<int>       xbond12;
    std::vector<IntPair>   xbond13;
    std::vector<IntTriple> xbond14;
    std::vector<IntQuad>   xbond15;

    std::vector<BondStretch> stretch; // same length as cross-body 1-2 list
    std::vector<BondBend>    bend;    // same length as   " 1-3 list
    std::vector<BondTorsion> torsion; // same length as   " 1-4 list
};


class Bond {
public:
    Bond() { }
    Bond(int atom1, int atom2) : atoms(atom1,atom2) { 
        assert(isValid());
    }
    bool isValid() const {return atoms.isValid();}

    IntPair atoms;
};

class ChargeProperties {
public:
    Real     netCharge;
    Vec3     centerOfCharge;
    Vec3     dipoleMoment;
    SymMat33 quadrupoleMoment;
};

class GeometricProperties {
public:
    Transform obbFrame;
    Vec3      obbHalfLengths;
    Real      boundingSphereRadius;
    Vec3      boundingSphereCenter;
};


class RigidGroup {
public:
    RigidGroup() : rigidGroupId(-1), topologicalCacheValid(false) { }
    RigidGroup(const char* nm)
      : rigidGroupId(-1), name(nm), topologicalCacheValid(false) {
        assert(isValid());
    }

    bool isValid() const {return rigidGroupId >= 0;}

    bool isTopologicalCacheValid() const   {return topologicalCacheValid;}
    void invalidateTopologicalCache()      {topologicalCacheValid=false;}

    void placeAtom(int atomId, const Vec3& station) {
        std::pair<AtomPlacementSet::iterator, bool> ret =
            atomPlacements.insert(AtomPlacement(atomId,station));
        assert(ret.second); // must not have been there already
    }

    void placeRigidGroup(int groupId, const Transform& placement) {
        std::pair<RigidGroupPlacementSet::iterator, bool> ret =
            rigidGroupPlacements.insert(RigidGroupPlacement(groupId,placement));
        assert(ret.second); // must not have been there already

        //TODO: check for loops
    }

    // True if the atom has been placed in this group, or if any group that has
    // been placed here contains the atom, or a subgroup containing the atom
    // recursively.
    bool containsAtom      (int atomId,  const DuMMForceFieldSubsystemRep& mm) const;
    bool containsRigidGroup(int groupId, const DuMMForceFieldSubsystemRep& mm) const;

    // Transform all our atoms (including those in subgroups) by the indicated 
    // transformation and append the result to the supplied AtomPlacement container.
    // The results are not in any particular order. Note that for a subgroup we
    // compose the provided transform with that of the subgroup before applying
    // the results to the contained atoms and, recursively, to the contained sub-subgroups.
    void calculateAllAtomArray(DuMMForceFieldSubsystemRep& mm);

    // Transform the atoms on my all atom array into the parent frame, and append
    // them to the supplied placement array.
    void xformAllAtomArray(const Transform& X_PC, AtomPlacementArray& xformedAtoms) {
        assert(isTopologicalCacheValid());
        for (int i=0; i < (int)allAtoms.size(); ++i) {
            const AtomPlacement& ap = allAtoms[i];
            xformedAtoms.push_back(AtomPlacement(ap.atomId, X_PC*ap.station));
        }
    }

    // Recursively calculate composite properties for this group and all the
    // groups it contains. All groups were marked "invalid" at the beginning
    // of this step.
    void realizeTopologicalCache(DuMMForceFieldSubsystemRep& mm) {
        if (topologicalCacheValid)
            return;

        // Transform all attached atoms into this group's frame, then sort
        // the resulting list to make it easy to check for duplicates.
        allAtoms.clear();
        calculateAllAtomArray(mm);
        std::sort(allAtoms.begin(), allAtoms.end());
        // TODO: check for duplicates
        // TODO: mass properties
        // TODO: charge properties

        topologicalCacheValid = true;
    }


    void dump() const {
        printf("    rigidGroupId=%d(%s)\n", rigidGroupId, name);
        printf("      atom placements: ");
        AtomPlacementSet::const_iterator ap = atomPlacements.begin();
        while (ap != atomPlacements.end())
            std::cout << " " << ap->atomId << ":" << ap->station;
        printf("      group placements:\n");
        RigidGroupPlacementSet::const_iterator gp = rigidGroupPlacements.begin();
        while (gp != rigidGroupPlacements.end())
            std::cout << "      " << gp->rigidGroupId << ":" << gp->placement;
        std::cout << "      topological cache valid? " 
                  << isTopologicalCacheValid() << std::endl;
    }

    void clearAllCalculatedData() {
        topologicalCacheValid = false;
        allAtoms.clear();
        massProps      = MassProperties();
        chargeProps    = ChargeProperties();
        geometricProps = GeometricProperties();
    }

public:
        // TOPOLOGICAL STATE VARIABLES
        //   Filled in during construction.
    int                     rigidGroupId;
    std::string             name;
    AtomPlacementSet        atomPlacements;
    RigidGroupPlacementSet  rigidGroupPlacements;

        // TOPOLOGICAL CACHE ENTRIES
        //   These are calculated in realizeConstruction() from topological
        //   state variables (from here or others in the DuMM class).

    bool topologicalCacheValid;

    // This is an expansion of all the atom & group placements, with
    // all stations transformed to this group's frame.
    AtomPlacementArray  allAtoms; // sorted by atomId

    // These reflect composite properties built from the allAtoms list.
    MassProperties      massProps;
    ChargeProperties    chargeProps;
    GeometricProperties geometricProps;
};

// This is just a RigidGroup with some additional runtime data structures.
// In particular we want a complete list of all the atoms from all groups
// attached to this body, with their placements as measured in the body frame.
class Body : public RigidGroup {
public:
    Body() : bodyNum(-1) { }

    int bodyNum; // this the body number as known to the matter subsystem
    std::vector<int> shadowBodies; // if needed
};

// Assume units:
//    Ref: http://physics.nist.gov/constants (2002 CODATA)
//    charge  e=charge on proton=1.60217653e-19C
//    Avogadro's number N0=6.0221415e23 atoms/mole       
//    length  A=Angstroms=1e-10 m=0.1nm
//    mass    Da=g/mole
//    time    ps
//    That implies force = Da-A/ps^2
//    atomic mass unit = 1/12 mass(C)=1.66053886e-24 g
//      (specifically Carbon-12, unbound, in its rest state)
//    mass of 1 mole of Carbon-12 = 12g (exact), thus mass
//      of one Carbon-12 atom is 12 Da.
//    energy kcal/mole = 418.4 Da-A^2/ps^2
//    e0 in e^2/(A-kcal/mole)
//      = 8.854187817e-12 C^2/(m-J)
//          * (1/1.60217653e-19)^2 * 4184/6.0221415e23 * 1e-10
//      = 2.3964519142e-4
//    1/(4*pi*e0) = 332.06371
//    speed of light c=2.99792458e8 m/s (exact)
//    Joules(N-m)/Kcal = 4184 (exact)
//
// Note: we have to use consistent force units, meaning
//   Da-A/ps^2
//
// Jay Ponder's Tinker units, as of email 8/30/06:
// In any case, I've just updated all TINKER units to the following:
//       parameter (avogadro=6.0221415d+23)
//       parameter (boltzmann=0.8314472d0)
//       parameter (gasconst=1.9872065d-3)
//       parameter (lightspd=2.99792458d-2)
//       parameter (bohr=0.5291772108d0)
//       parameter (joule=4.184d0)
//       parameter (evolt=27.2113845d0)
//       parameter (hartree=627.509472d0)
//       parameter (electric=332.06371d0)
//       parameter (debye=4.8033324d0)
//       parameter (prescon=6.85695d+4)
//       parameter (convert=4.184d+2)


class DuMMForceFieldSubsystemRep : public ForceSubsystemRep {
public:
    DuMMForceFieldSubsystemRep()
      : ForceSubsystemRep("DuMMForceFieldSubsystem", "0.0.1")
    {
        topologicalCacheValid = false;

        vdwScale12=coulombScale12=vdwScale13=coulombScale13=0;
        vdwScale14=coulombScale14=vdwScale15=coulombScale15=1;
        loadElements();
        const int gid = addRigidGroup(RigidGroup("free atoms and groups"));
        assert(gid==0);
    }

    bool isValidElement(int atomicNumber) const {
        return 1 <= atomicNumber && atomicNumber < (int)elements.size() 
                && elements[atomicNumber].isValid();
    }

    bool isValidAtom(int atomNum) const {
        return 0 <= atomNum && atomNum < (int)atoms.size() && atoms[atomNum].isValid();
    }

    bool isValidBond(int bondNum) const {
        return 0 <= bondNum && bondNum < (int)bonds.size() && bonds[bondNum].isValid();
    }

    bool isValidRigidGroup(int rigidGroupId) const {
        return 0 <= rigidGroupId && rigidGroupId < (int)rigidGroups.size()
                && rigidGroups[rigidGroupId].isValid();
    }

    bool isValidBody(int bodyId) const {
        return 0 <= bodyId && bodyId < (int)bodies.size() && bodies[bodyId].isValid();
    }

    bool isValidChargedAtomType(int typeNum) const {
        return 0 <= typeNum && typeNum < (int)chargedAtomTypes.size() 
               && chargedAtomTypes[typeNum].isValid();
    }

    bool isValidAtomClass(int classNum) const {
        return 0 <= classNum && classNum < (int)atomClasses.size() 
               && atomClasses[classNum].isValid();
    }


    // We scale short range interactions but only for bonds which cross bodies.
    void scaleBondedAtoms(const Atom& a, Vector& vdwScale, Vector& coulombScale) const;
    void unscaleBondedAtoms(const Atom& a, Vector& vdwScale, Vector& coulombScale) const;

    void applyMixingRule(Real ri, Real rj, Real ei, Real ej, Real& dmin, Real& emin) const
    {
        Real rmin;
        vdwCombineWaldmanHagler(ri,rj,ei,ej,rmin,emin); // TODO: choices
        //vdwCombineJorgensen(ri,rj,ei,ej,rmin,emin);
        //vdwCombineHalgrenHHG(ri,rj,ei,ej,rmin,emin);
        //vdwCombineKong(ri,rj,ei,ej,rmin,emin);
        // NO NO NO!! :
        //vdwCombineLorentzBerthelot(ri,rj,ei,ej,rmin,emin);
        dmin = 2*rmin;
    }

    int addRigidGroup(const RigidGroup& group) {
        const int rigidGroupId = (int)rigidGroups.size();
        rigidGroups.push_back(group);
        rigidGroups[rigidGroupId].rigidGroupId = rigidGroupId;
        return rigidGroupId;
    }
    RigidGroup& updRigidGroup(int groupId) {
        assert(isValidRigidGroup(groupId));
        return rigidGroups[groupId];
    }
    const RigidGroup& getRigidGroup(int groupId) const {
        assert(isValidRigidGroup(groupId));
        return rigidGroups[groupId];
    }

    void placeAtomInRigidGroup(int atomId, int rigidGroupId, const Vec3& station) {
        assert(isValidAtom(atomId));
        assert(!getRigidGroup(rigidGroupId).containsAtom(atomId, *this));
        updRigidGroup(rigidGroupId).placeAtom(atomId, station);
    }

    void placeRigidGroupInRigidGroup
       (int childGroupId, int parentGroupId, const Transform& placement) 
    {
        assert(isValidRigidGroup(childGroupId) && isValidRigidGroup(parentGroupId));
        assert(!getRigidGroup(parentGroupId).containsRigidGroup(childGroupId, *this));
        updRigidGroup(parentGroupId).placeRigidGroup(childGroupId, placement);
    }

    int addAtom(int chargedAtomTypeId)
    {
        assert(isValidChargedAtomType(chargedAtomTypeId));
        const int atomId = (int)atoms.size();
        atoms.push_back(Atom(chargedAtomTypeId, atomId));
        return atomId;    
    }

    int addBond(int atom1, int atom2)
    {
        assert(isValidAtom(atom1) && isValidAtom(atom2));
        assert(atom1 != atom2);

        // Ensure that atom1 < atom2
        if (atom1 > atom2)
            std::swap(atom1,atom2);

        Atom& a1 = atoms[atom1];
        Atom& a2 = atoms[atom2];

        if (a1.isBondedTo(atom2)) {
            assert(a2.isBondedTo(atom1));
            for (int i=0; i < (int)bonds.size(); ++i)
                if (bonds[i].atoms[0]==atom1 && bonds[i].atoms[1]==atom2)
                    return i;
            assert(!"missing bond");
        }

        bonds.push_back(Bond(atom1,atom2));
        a1.bond12.push_back(atom2);
        a2.bond12.push_back(atom1);
        return (int)bonds.size() - 1;
    }

    int getNAtoms() const {return (int)atoms.size();}
    int getNBonds() const {return (int)bonds.size();}

    int getChargedAtomTypeNum(int atomId) const {
        assert(isValidAtom(atomId));
        return atoms[atomId].chargedAtomTypeId;
    }

    int getAtomClassNum(int atomId) const {
        assert(isValidAtom(atomId));
        const ChargedAtomType& type = chargedAtomTypes[getChargedAtomTypeNum(atomId)];
        return type.atomClass;
    }

    int getAtomElementNum(int atomId) const {
        assert(isValidAtom(atomId));
        const AtomClass& cl = atomClasses[getAtomClassNum(atomId)];
        return cl.element;
    }

    Real getAtomMass(int atomId) const {
        assert(isValidAtom(atomId));
        const Element& e = elements[getAtomElementNum(atomId)];
        return e.mass;
    }

    const Vec3& getAtomDefaultColor(int atomId) const {
        assert(isValidAtom(atomId));
        const Element& e = elements[getAtomElementNum(atomId)];
        return e.defaultColor;
    }

    Real getAtomRadius(int atomId) const {
        assert(isValidAtom(atomId));
        const AtomClass& cl = atomClasses[getAtomClassNum(atomId)];
        return cl.vdwRadius;
    }

    const Vec3& getAtomStation(int atomId) const {
        assert(isValidAtom(atomId));
        assert(topologicalCacheValid); // can't know station until after realizeConstruction()
        return atoms[atomId].station_B;
    }

    int getAtomBody(int atomId) const {
        assert(isValidAtom(atomId));
        assert(topologicalCacheValid); // can't know body until after realizeConstruction()
        return atoms[atomId].bodyId;
    }

    int getBondAtom(int b, int which) const {
        assert(isValidBond(b) && (which==0 || which==1));
        return bonds[b].atoms[which];
    }


    void addAtomClass(const AtomClass& atomClass);
    void addChargedAtomType(const ChargedAtomType& atomType);
    void addBondStretch(int class1, int class2, const BondStretch& bs);
    void addBondBend(int class1, int class2, int class3, const BondBend& bb);
    void addBondTorsion(int class1, int class2, int class3, int class4, 
                        const TorsionTerm& tt1, 
                        const TorsionTerm& tt2=TorsionTerm(),
                        const TorsionTerm& tt3=TorsionTerm());

    const BondStretch& getBondStretch(int class1, int class2) const;
    const BondBend&    getBondBend   (int class1, int class2, int class3) const;
    const BondTorsion& getBondTorsion(int class1, int class2, int class3, int class4) const;

    // Save the reciprocal so we can multiply instead of divide.
    void setVdw12ScaleFactor(Real fac) {assert(0<=fac&&fac<=1); vdwScale12=fac;}
    void setVdw13ScaleFactor(Real fac) {assert(0<=fac&&fac<=1); vdwScale13=fac;}
    void setVdw14ScaleFactor(Real fac) {assert(0<=fac&&fac<=1); vdwScale14=fac;}
    void setVdw15ScaleFactor(Real fac) {assert(0<=fac&&fac<=1); vdwScale15=fac;}

    void setCoulomb12ScaleFactor(Real fac) {assert(0<=fac&&fac<=1); coulombScale12=fac;}
    void setCoulomb13ScaleFactor(Real fac) {assert(0<=fac&&fac<=1); coulombScale13=fac;}
    void setCoulomb14ScaleFactor(Real fac) {assert(0<=fac&&fac<=1); coulombScale14=fac;}
    void setCoulomb15ScaleFactor(Real fac) {assert(0<=fac&&fac<=1); coulombScale15=fac;}


    void realizeConstruction(State& s) const;

    void realizeModeling(State& s) const {
        // Sorry, no choices available at the moment.
    }

    void realizeParameters(const State& s) const {
        // Nothing to compute here.
    }

    void realizeTime(const State& s) const {
        // Nothing to compute here.
    }

    void realizeConfiguration(const State& s) const {
        // Nothing to compute here.
    }

    void realizeMotion(const State& s) const {
        // Nothing to compute here.
    }


    void realizeDynamics(const State& s) const;

    void realizeReaction(const State& s) const {
        // Nothing to compute here.
    }

    void dump() const;

    DuMMForceFieldSubsystemRep* cloneSubsystemRep() const {
        return new DuMMForceFieldSubsystemRep(*this);
    }

private:
    void loadElements();

    void invalidateAllTopologicalCacheEntries() {
        topologicalCacheValid = false;

        // If any of these objects are invalid, the invalidateTopologicalCache()
        // call does nothing (i.e., it doesn't blow up!).

        // molecule
        for (int i=0; i < (int)atoms.size(); ++i)
            atoms[i].invalidateTopologicalCache();
        for (int i=0; i < (int)rigidGroups.size(); ++i)
            rigidGroups[i].invalidateTopologicalCache();
        for (int i=0; i < (int)bodies.size(); ++i)
            bodies[i].invalidateTopologicalCache();

        // force field
        for (int i=0; i < (int)atomClasses.size(); ++i)
            atomClasses[i].invalidateTopologicalCache();
    }

private:
        // TOPOLOGICAL STATE VARIABLES
        //   Filled in during construction.

    // molecule

    std::vector<Atom>       atoms;
    std::vector<Bond>       bonds;
    std::vector<RigidGroup> rigidGroups;
    // This defines the partitioning of atoms onto the matter subsystem's bodies.
    // The indices here correspond to the body numbers. Only entries for bodies on
    // which our atoms have been attached will be valid.
    std::vector<Body>       bodies;

    // force field

    // Force field description. These are not necessarily fully populated;
    // check the "isValid()" method to see if anything is there.
    std::vector<Element>             elements;
    std::vector<AtomClass>           atomClasses;
    std::vector<ChargedAtomType>     chargedAtomTypes;

    // These relate atom classes, not charged atom types.
    std::map<IntPair,   BondStretch> bondStretch;
    std::map<IntTriple, BondBend>    bondBend;
    std::map<IntQuad,   BondTorsion> bondTorsion;

    // Scale factors for nonbonded forces when applied to
    // atoms which are near in the graph formed by the bonds.
    Real vdwScale12, coulombScale12;    // default 0,0
    Real vdwScale13, coulombScale13;    // default 0,0
    Real vdwScale14, coulombScale14;    // default 1,1
    Real vdwScale15, coulombScale15;    // default 1,1

        // TOPOLOGICAL CACHE ENTRIES
        //   These are calculated in realizeConstruction() from topological
        //   state variables (from here or others in the DuMM class).
    bool topologicalCacheValid;
};


    /////////////////////////////
    // DuMMForceFieldSubsystem //
    /////////////////////////////

/*static*/ bool 
DuMMForceFieldSubsystem::isInstanceOf(const ForceSubsystem& s) {
    return DuMMForceFieldSubsystemRep::isA(s.getRep());
}
/*static*/ const DuMMForceFieldSubsystem&
DuMMForceFieldSubsystem::downcast(const ForceSubsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const DuMMForceFieldSubsystem&>(s);
}
/*static*/ DuMMForceFieldSubsystem&
DuMMForceFieldSubsystem::updDowncast(ForceSubsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<DuMMForceFieldSubsystem&>(s);
}

const DuMMForceFieldSubsystemRep& 
DuMMForceFieldSubsystem::getRep() const {
    return dynamic_cast<const DuMMForceFieldSubsystemRep&>(*rep);
}
DuMMForceFieldSubsystemRep&       
DuMMForceFieldSubsystem::updRep() {
    return dynamic_cast<DuMMForceFieldSubsystemRep&>(*rep);
}

DuMMForceFieldSubsystem::DuMMForceFieldSubsystem() {
    rep = new DuMMForceFieldSubsystemRep();
    rep->setMyHandle(*this);
}

void DuMMForceFieldSubsystem::defineAtomClass
   (int id, const char* className, int element, int valence, 
    Real vdwRadius, Real vdwWellDepth)
{
    updRep().addAtomClass(AtomClass(id, className, element, valence, 
                                    vdwRadius, vdwWellDepth));
}

void DuMMForceFieldSubsystem::defineChargedAtomType
   (int id, const char* typeName, int atomClass, Real partialCharge)
{
    updRep().addChargedAtomType(ChargedAtomType(id, typeName, atomClass, partialCharge));
}

void DuMMForceFieldSubsystem::defineBondStretch
   (int class1, int class2, Real stiffnessInKcalPerASq, Real nominalLengthInA)
{
    updRep().addBondStretch(class1, class2, BondStretch(stiffnessInKcalPerASq,nominalLengthInA));
}

void DuMMForceFieldSubsystem::defineBondBend
   (int class1, int class2, int class3, Real stiffnessInKcalPerRadSq, Real nominalAngleInDegrees)
{
    updRep().addBondBend(class1, class2, class3, BondBend(stiffnessInKcalPerRadSq,nominalAngleInDegrees));
}

void DuMMForceFieldSubsystem::defineBondTorsion
   (int class1, int class2, int class3, int class4, 
    int periodicity1, Real amp1InKcal, Real phase1InDegrees)
{
    updRep().addBondTorsion(class1, class2, class3, class4, 
                            TorsionTerm(periodicity1,amp1InKcal,phase1InDegrees));
}
void DuMMForceFieldSubsystem::defineBondTorsion
   (int class1, int class2, int class3, int class4, 
    int periodicity1, Real amp1InKcal, Real phase1InDegrees,
    int periodicity2, Real amp2InKcal, Real phase2InDegrees)
{
    updRep().addBondTorsion(class1, class2, class3, class4, 
                            TorsionTerm(periodicity1,amp1InKcal,phase1InDegrees),
                            TorsionTerm(periodicity2,amp2InKcal,phase2InDegrees));
}
void DuMMForceFieldSubsystem::defineBondTorsion
   (int class1, int class2, int class3, int class4, 
    int periodicity1, Real amp1InKcal, Real phase1InDegrees,
    int periodicity2, Real amp2InKcal, Real phase2InDegrees,
    int periodicity3, Real amp3InKcal, Real phase3InDegrees)
{
    updRep().addBondTorsion(class1, class2, class3, class4, 
                            TorsionTerm(periodicity1,amp1InKcal,phase1InDegrees),
                            TorsionTerm(periodicity2,amp2InKcal,phase2InDegrees),
                            TorsionTerm(periodicity3,amp3InKcal,phase3InDegrees));
}

void DuMMForceFieldSubsystem::setVdw12ScaleFactor(Real fac) {updRep().setVdw12ScaleFactor(fac);}
void DuMMForceFieldSubsystem::setVdw13ScaleFactor(Real fac) {updRep().setVdw13ScaleFactor(fac);}
void DuMMForceFieldSubsystem::setVdw14ScaleFactor(Real fac) {updRep().setVdw14ScaleFactor(fac);}
void DuMMForceFieldSubsystem::setVdw15ScaleFactor(Real fac) {updRep().setVdw15ScaleFactor(fac);}

void DuMMForceFieldSubsystem::setCoulomb12ScaleFactor(Real fac) {updRep().setCoulomb12ScaleFactor(fac);}
void DuMMForceFieldSubsystem::setCoulomb13ScaleFactor(Real fac) {updRep().setCoulomb13ScaleFactor(fac);}
void DuMMForceFieldSubsystem::setCoulomb14ScaleFactor(Real fac) {updRep().setCoulomb14ScaleFactor(fac);}
void DuMMForceFieldSubsystem::setCoulomb15ScaleFactor(Real fac) {updRep().setCoulomb15ScaleFactor(fac);}

int DuMMForceFieldSubsystem::addAtom(int chargedAtomType)
{
    return updRep().addAtom(chargedAtomType);
}

int DuMMForceFieldSubsystem::addBond(int atom1, int atom2)
{
    return updRep().addBond(atom1, atom2);
}

int DuMMForceFieldSubsystem::getNAtoms() const {
    return getRep().getNAtoms();
}
int DuMMForceFieldSubsystem::getNBonds() const {
    return getRep().getNBonds();
}
int DuMMForceFieldSubsystem::getBondAtom(int bond, int which) const {
    return getRep().getBondAtom(bond, which);
}

Real DuMMForceFieldSubsystem::getAtomMass(int atomNum) const {
    return getRep().getAtomMass(atomNum);
}
Vec3 DuMMForceFieldSubsystem::getAtomDefaultColor(int atomNum) const {
    return getRep().getAtomDefaultColor(atomNum);
}
Real DuMMForceFieldSubsystem::getAtomRadius(int atomNum) const {
    return getRep().getAtomRadius(atomNum);
}
Vec3 DuMMForceFieldSubsystem::getAtomStation(int atomNum) const {
    return getRep().getAtomStation(atomNum);
}
int DuMMForceFieldSubsystem::getAtomBody(int atomNum) const {
    return getRep().getAtomBody(atomNum);
}

void DuMMForceFieldSubsystem::dump() const {
    return getRep().dump();
}



    ////////////////////////////////
    // DuMMForceFieldSubsystemRep //
    ////////////////////////////////


void DuMMForceFieldSubsystemRep::addAtomClass(const AtomClass& atomClass) {
    const int id = atomClass.atomClassId;
    assert(id >= 0);
    if (id >= (int)atomClasses.size())
        atomClasses.resize(id+1);
    assert(!atomClasses[id].isValid());
    atomClasses[id] = atomClass;
}

void DuMMForceFieldSubsystemRep::addChargedAtomType(const ChargedAtomType& atomType) {
    const int id = atomType.chargedAtomTypeId;
    assert(id >= 0);
    if (id >= (int)chargedAtomTypes.size())
        chargedAtomTypes.resize(id+1);
    assert(!chargedAtomTypes[id].isValid());
    chargedAtomTypes[id] = atomType;
}

void DuMMForceFieldSubsystemRep::addBondStretch(int class1, int class2, const BondStretch& bs) {
    assert(isValidAtomClass(class1) && isValidAtomClass(class2));
    // Canonicalize the pair to have lowest class # first
    const IntPair key(class1,class2,true);
    std::pair<std::map<IntPair,BondStretch>::iterator, bool> ret = 
      bondStretch.insert(std::pair<IntPair,BondStretch>(key,bs));
    assert(ret.second); // must not have been there already
}

const BondStretch& 
DuMMForceFieldSubsystemRep::getBondStretch(int class1, int class2) const {
    const IntPair key(class1,class2,true);
    std::map<IntPair,BondStretch>::const_iterator bs = bondStretch.find(key);
    assert(bs != bondStretch.end());
    return bs->second;
}

void DuMMForceFieldSubsystemRep::addBondBend(int class1, int class2, int class3, const BondBend& bb) {
    assert(isValidAtomClass(class1) && isValidAtomClass(class2) && isValidAtomClass(class3));
    // Canonicalize the triple to have lowest type # first
    const IntTriple key(class1, class2, class3, true);
    std::pair<std::map<IntTriple,BondBend>::iterator, bool> ret = 
      bondBend.insert(std::pair<IntTriple,BondBend>(key,bb));
    assert(ret.second); // must not have been there already
}

const BondBend& 
DuMMForceFieldSubsystemRep::getBondBend(int class1, int class2, int class3) const {
    const IntTriple key(class1, class2, class3, true);
    std::map<IntTriple,BondBend>::const_iterator bb = bondBend.find(key);
    assert(bb != bondBend.end());
    return bb->second;
}

void DuMMForceFieldSubsystemRep::addBondTorsion
   (int class1, int class2, int class3, int class4, 
    const TorsionTerm& tt1, const TorsionTerm& tt2, const TorsionTerm& tt3) 
{
    assert(isValidAtomClass(class1) && isValidAtomClass(class2) 
            && isValidAtomClass(class3) && isValidAtomClass(class4));
    assert(tt1.isValid());

    // Canonicalize the quad to have lowest type # first
    const IntQuad key(class1, class2, class3, class4, true);
    BondTorsion bt;
    if (tt1.isValid()) bt.addTerm(tt1);
    if (tt2.isValid()) bt.addTerm(tt2);
    if (tt3.isValid()) bt.addTerm(tt3);

    std::pair<std::map<IntQuad,BondTorsion>::iterator, bool> ret = 
      bondTorsion.insert(std::pair<IntQuad,BondTorsion>(key,bt));
    assert(ret.second); // must not have been there already
}

const BondTorsion& 
DuMMForceFieldSubsystemRep::getBondTorsion
   (int class1, int class2, int class3, int class4) const
{
    const IntQuad key(class1, class2, class3, class4, true);
    std::map<IntQuad,BondTorsion>::const_iterator bt = bondTorsion.find(key);
    assert(bt != bondTorsion.end());
    return bt->second;
}

void DuMMForceFieldSubsystemRep::realizeConstruction(State& s) const {
    if (topologicalCacheValid)
        return; // already got this far

    // We need to write once onto the 'cache' portion of the object once
    // the topology is known.
    DuMMForceFieldSubsystemRep* mutableThis = 
        const_cast<DuMMForceFieldSubsystemRep*>(this);

    mutableThis->invalidateAllTopologicalCacheEntries();

        // force field

    // Calculate effective van der Waals parameters for all 
    // pairs of atom classes. We only fill in the diagonal
    // and upper triangle; that is, each class contains
    // parameters for like classes and all classes whose
    // (arbitrary) class number is higher.
    for (int i=0; i < (int)atomClasses.size(); ++i) {
        if (!atomClasses[i].isValid()) continue;

        AtomClass& iclass = mutableThis->atomClasses[i];
        iclass.vdwDij.resize((int)atomClasses.size()-i, CNT<Real>::getNaN());
        iclass.vdwEij.resize((int)atomClasses.size()-i, CNT<Real>::getNaN()); 
        for (int j=i; j < (int)atomClasses.size(); ++j) {
            const AtomClass& jclass = atomClasses[j];
            if (jclass.isValid())
                applyMixingRule(iclass.vdwRadius,    jclass.vdwRadius,
                                iclass.vdwWellDepth, jclass.vdwWellDepth,
                                iclass.vdwDij[j-i],  iclass.vdwEij[j-i]);

        }
    }

        // molecule

    // Process groups & bodies (bodies are treated as top-level groups)

    // We process groups recursively, so we need to allow the groups writable
    // access to the main DuMM object (i.e., "this").
    for (int gnum=0; gnum < (int)rigidGroups.size(); ++gnum) {
        RigidGroup& g = mutableThis->rigidGroups[gnum];
        assert(g.isValid()); // Shouldn't be any unused group numbers.
        g.realizeTopologicalCache(*mutableThis);
    }

    // Bodies, on the other hand, are always top level groups and the
    // calculation here assumes that all the groups have been processed.
    // Thus bodies need only read access to the main DuMM object, 
    // although we're passign the mutable one in so we can use the
    // same routine (TODO).
    for (int bnum=0; bnum < (int)bodies.size(); ++bnum) {
        Body& b = mutableThis->bodies[bnum];
        if (!b.isValid())
            continue; // OK for these to be unused.
        b.realizeTopologicalCache(*mutableThis);
    }

    // Assign body & station to every atom that has been assigned to a body.
    for (int anum=0; anum < (int)atoms.size(); ++anum) {
        Atom& a = mutableThis->atoms[anum];
        a.bodyId = -1;
    }
    for (int bnum=0; bnum < (int)bodies.size(); ++bnum) {
        const Body& b = bodies[bnum];
        if (!b.isValid())
            continue;   // Unused body numbers are OK.

        for (int i=0; i < (int)b.allAtoms.size(); ++i) {
            const AtomPlacement& ap = b.allAtoms[i];   assert(ap.isValid());
            Atom& a = mutableThis->atoms[ap.atomId]; assert(a.isValid());
            assert(a.bodyId == -1); // Can only be on one body!!
            a.bodyId    = bnum;
            a.station_B = ap.station;
        }
    }
    for (int anum=0; anum < (int)atoms.size(); ++anum) {
        const Atom& a = atoms[anum];
        assert(a.bodyId >= 0); // TODO catch unassigned atoms
    }

    // need to chase bonds to fill in the bonded data
    // Be sure only to find the *shortest* path between two atoms
    for (int anum=0; anum < (int)atoms.size(); ++anum) {
        Atom& a = mutableThis->atoms[anum];
        std::set<int> allBondedSoFar;   // to avoid duplicate paths

        // Only the bond12 list should be filled in at the moment. We'll sort
        // all the lists when they're done for good hygiene.
        std::sort(a.bond12.begin(), a.bond12.end());

        // Add this atom and its direct (1-2) bonds to the list of all bonded atoms.
        allBondedSoFar.insert(anum);
        for (int j=0; j < (int)a.bond12.size(); ++j)
            allBondedSoFar.insert(a.bond12[j]);

        // Find longer bond paths by building each list in turn from
        // the direct bonds of the atoms in the previous list.

        // build the bond13 list
        a.bond13.clear();
        for (int j=0; j < (int)a.bond12.size(); ++j) {
            const Atom& a12 = atoms[a.bond12[j]];
            const AtomArray& a12_12 = a12.bond12;
            for (int k=0; k < (int)a12_12.size(); ++k) {
                const int newAtom = a12_12[k];
                if (allBondedSoFar.find(newAtom) != allBondedSoFar.end())
                    continue; // there was already a shorter path
                allBondedSoFar.insert(newAtom);
                a.bond13.push_back(IntPair(a.bond12[j], newAtom));
            }
        }
        std::sort(a.bond13.begin(), a.bond13.end());

        // build the bond14 list
        a.bond14.clear();
        for (int j=0; j < (int)a.bond13.size(); ++j) {
            const Atom& a13 = atoms[a.bond13[j][1]];
            const AtomArray& a13_12 = a13.bond12;
            for (int k=0; k < (int)a13_12.size(); ++k) {
                const int newAtom = a13_12[k];
                if (allBondedSoFar.find(newAtom) != allBondedSoFar.end())
                    continue; // there was already a shorter path
                allBondedSoFar.insert(newAtom);
                a.bond14.push_back(IntTriple(a.bond13[j][0], a.bond13[j][1], newAtom));
            }
        }
        std::sort(a.bond14.begin(), a.bond14.end());

        // build the bond15 list
        a.bond15.clear();
        for (int j=0; j < (int)a.bond14.size(); ++j) {
            const Atom& a14 = atoms[a.bond14[j][1]];
            const AtomArray& a14_12 = a14.bond12;
            for (int k=0; k < (int)a14_12.size(); ++k) {
                const int newAtom = a14_12[k];
                if (allBondedSoFar.find(newAtom) != allBondedSoFar.end())
                    continue; // there was already a shorter path
                allBondedSoFar.insert(newAtom);
                a.bond15.push_back(IntQuad(a.bond14[j][0], a.bond14[j][1], a.bond14[j][2], newAtom));
            }
        }
        std::sort(a.bond15.begin(), a.bond15.end());

        // Fill in the cross-body bond lists. We only keep atoms which
        // are on a different body.
        a.xbond12.clear(); a.xbond13.clear(); a.xbond14.clear(); a.xbond15.clear();
        for (int j=0; j < (int)a.bond12.size(); ++j)
            if (atoms[a.bond12[j]].bodyId != a.bodyId)
                a.xbond12.push_back(a.bond12[j]);

        for (int j=0; j < (int)a.bond13.size(); ++j)
            if (   atoms[a.bond13[j][0]].bodyId != a.bodyId
                || atoms[a.bond13[j][1]].bodyId != a.bodyId)
                a.xbond13.push_back(a.bond13[j]);

        for (int j=0; j < (int)a.bond14.size(); ++j)
            if (   atoms[a.bond14[j][0]].bodyId != a.bodyId
                || atoms[a.bond14[j][1]].bodyId != a.bodyId
                || atoms[a.bond14[j][2]].bodyId != a.bodyId)
                a.xbond14.push_back(a.bond14[j]);

        for (int j=0; j < (int)a.bond15.size(); ++j)
            if (   atoms[a.bond15[j][0]].bodyId != a.bodyId
                || atoms[a.bond15[j][1]].bodyId != a.bodyId
                || atoms[a.bond15[j][2]].bodyId != a.bodyId
                || atoms[a.bond15[j][3]].bodyId != a.bodyId)
                a.xbond15.push_back(a.bond15[j]);

        // Save a BondStretch entry for each 1-2 bond
        a.stretch.resize(a.xbond12.size());
        for (int b12=0; b12 < (int)a.xbond12.size(); ++b12)
            a.stretch[b12] = getBondStretch(getAtomClassNum(anum), 
                                            getAtomClassNum(a.xbond12[b12]));

        // Save a BondBend entry for each 1-3 bond
        a.bend.resize(a.xbond13.size());
        for (int b13=0; b13 < (int)a.xbond13.size(); ++b13)
            a.bend[b13] = getBondBend(getAtomClassNum(anum), 
                                      getAtomClassNum(a.xbond13[b13][0]), 
                                      getAtomClassNum(a.xbond13[b13][1]));

        // Save a BondTorsion entry for each 1-4 bond
        a.torsion.resize(a.xbond14.size());
        for (int b14=0; b14 < (int)a.xbond14.size(); ++b14)
            a.torsion[b14] = getBondTorsion(getAtomClassNum(anum), 
                                            getAtomClassNum(a.xbond14[b14][0]), 
                                            getAtomClassNum(a.xbond14[b14][1]),
                                            getAtomClassNum(a.xbond14[b14][2]));
    }

    mutableThis->topologicalCacheValid = true;
}

// Cost of processing here (in flops): XXX
// Strategy:
//   for each body b we know about here
//     for each atom a on b
//          set scale factors on bonded atoms
//          for each body c > b
//            for each atom ac on c
//                 compute vector r=ac-a and distance d=|r|
//                 compute vdw forces
//                 compute charge forces
//                 add force contribution to body
//          reset scale factors on bonded atoms
//

void DuMMForceFieldSubsystemRep::realizeDynamics(const State& s) const 
{
    const MultibodySystem& mbs    = getMultibodySystem(); // my owner
    const MatterSubsystem& matter = mbs.getMatterSubsystem();

    // Temps for scale factors; initialize to 1
    Vector vdwScale((int)atoms.size(), Real(1)); 
    Vector coulombScale((int)atoms.size(), Real(1));

    // Get access to system-global cache entries.
    Real&                  pe              = mbs.updPotentialEnergy(s);
    Vector_<SpatialVec>&   rigidBodyForces = mbs.updRigidBodyForces(s);

    for (int b1=0; b1 < (int)bodies.size(); ++b1) {
        const Transform&          X_GB1  = matter.getBodyConfiguration(s,b1);
        const AtomPlacementArray& alist1 = bodies[b1].allAtoms;

        for (int i=0; i < (int)alist1.size(); ++i) {
            const int       a1num = alist1[i].atomId;
            const Atom&     a1 = atoms[a1num];
            const ChargedAtomType& a1type  = chargedAtomTypes[a1.chargedAtomTypeId];
            int                    a1cnum  = a1type.atomClass;
            const AtomClass&       a1class = atomClasses[a1cnum];
            const Vec3      a1Station_G = X_GB1.R()*a1.station_B;
            const Vec3      a1Pos_G     = X_GB1.T() + a1Station_G;
            const Real      q1Fac = CoulombFac*a1type.partialCharge;

            // Bonded. Note that each bond will appear twice so we only process
            // it the time when its 1st atom has a lower ID than its last.

            // Bond stretch (1-2)
            for (int b12=0; b12 < (int)a1.xbond12.size(); ++b12) {
                const int a2num = a1.xbond12[b12];
                assert(a2num != a1num);
                if (a2num < a1num)
                    continue; // don't process this bond this time

                const Atom& a2 = atoms[a2num];
                const int b2 = a2.bodyId;
                assert(b2 != b1);
                const Transform& X_GB2   = matter.getBodyConfiguration(s, a2.bodyId);
                const Vec3       a2Station_G = X_GB2.R()*a2.station_B;
                const Vec3       a2Pos_G     = X_GB2.T() + a2Station_G;
                const Vec3       r = a2Pos_G - a1Pos_G;
                const Real       d = r.norm();

                // TODO: come up with something for when d is 0; for relaxation
                // just needs to push away from zero; what to do for dynamics?

                const BondStretch& bs = a1.stretch[b12];
                const Real         x  = d - bs.d0;

                const Real eStretch =  bs.k*x*x; // no factor of 1/2!
                const Real fStretch = -2*bs.k*x; // sign is as would be applied to a2
                const Vec3 f2 = (fStretch/d) * r;
                pe += eStretch;
                rigidBodyForces[b2] += SpatialVec( a2Station_G % f2, f2);   // 15 flops
                rigidBodyForces[b1] -= SpatialVec( a1Station_G % f2, f2);   // 15 flops
            }

            // Bond bend (1-2-3)
            for (int b13=0; b13 < (int)a1.xbond13.size(); ++b13) {
                const int a2num = a1.xbond13[b13][0];
                const int a3num = a1.xbond13[b13][1];
                assert(a3num != a1num);
                if (a3num < a1num)
                    continue; // don't process this bond this time

                const Atom& a2 = atoms[a2num];
                const Atom& a3 = atoms[a3num];
                const int b2 = a2.bodyId;
                const int b3 = a3.bodyId;
                assert(!(b2==b1 && b3==b1)); // shouldn't be on the list if all on 1 body

                // TODO: These might be the same body but for now we don't care.
                const Transform& X_GB2   = matter.getBodyConfiguration(s, a2.bodyId);
                const Transform& X_GB3   = matter.getBodyConfiguration(s, a3.bodyId);
                const Vec3       a2Station_G = X_GB2.R()*a2.station_B;
                const Vec3       a3Station_G = X_GB3.R()*a3.station_B;
                const Vec3       a2Pos_G     = X_GB2.T() + a2Station_G;
                const Vec3       a3Pos_G     = X_GB3.T() + a3Station_G;

                Real angle, energy;
                Vec3 f1, f2, f3;
                const BondBend& bb = a1.bend[b13];
                // atom 2 is the central one
                bb.harmonic(a2Pos_G, a1Pos_G, a3Pos_G, angle, energy, f2, f1, f3);

                pe += energy;
                rigidBodyForces[b1] += SpatialVec( a1Station_G % f1, f1);   // 15 flops
                rigidBodyForces[b2] += SpatialVec( a2Station_G % f2, f2);   // 15 flops
                rigidBodyForces[b3] += SpatialVec( a3Station_G % f3, f3);   // 15 flops
            }

            // Bond torsion (1-2-3-4)
            for (int b14=0; b14 < (int)a1.xbond14.size(); ++b14) {
                const int a2num = a1.xbond14[b14][0];
                const int a3num = a1.xbond14[b14][1];
                const int a4num = a1.xbond14[b14][2];
                assert(a4num != a1num);
                if (a4num < a1num)
                    continue; // don't process this bond this time

                const Atom& a2 = atoms[a2num];
                const Atom& a3 = atoms[a3num];
                const Atom& a4 = atoms[a4num];
                const int b2 = a2.bodyId;
                const int b3 = a3.bodyId;
                const int b4 = a4.bodyId;
                assert(!(b2==b1 && b3==b1 && b4==b1)); // shouldn't be on the list if all on 1 body

                // TODO: These might be the same body but for now we don't care.
                const Transform& X_GB2   = matter.getBodyConfiguration(s, a2.bodyId);
                const Transform& X_GB3   = matter.getBodyConfiguration(s, a3.bodyId);
                const Transform& X_GB4   = matter.getBodyConfiguration(s, a4.bodyId);
                const Vec3       a2Station_G = X_GB2.R()*a2.station_B;
                const Vec3       a3Station_G = X_GB3.R()*a3.station_B;
                const Vec3       a4Station_G = X_GB4.R()*a4.station_B;
                const Vec3       a2Pos_G     = X_GB2.T() + a2Station_G;
                const Vec3       a3Pos_G     = X_GB3.T() + a3Station_G;
                const Vec3       a4Pos_G     = X_GB4.T() + a4Station_G;

                Real angle, energy;
                Vec3 f1, f2, f3, f4;
                const BondTorsion& bt = a1.torsion[b14];
                bt.periodic(a1Pos_G, a2Pos_G, a3Pos_G, a4Pos_G, 
                            angle, energy, f1, f2, f3, f4);

                pe += energy;
                rigidBodyForces[b1] += SpatialVec( a1Station_G % f1, f1);   // 15 flops
                rigidBodyForces[b2] += SpatialVec( a2Station_G % f2, f2);   // 15 flops
                rigidBodyForces[b3] += SpatialVec( a3Station_G % f3, f3);   // 15 flops
                rigidBodyForces[b4] += SpatialVec( a4Station_G % f4, f4);   // 15 flops
            }

            scaleBondedAtoms(a1,vdwScale,coulombScale);
            for (int b2=b1+1; b2 < (int)bodies.size(); ++b2) {
                const Transform&          X_GB2  = matter.getBodyConfiguration(s,b2);
                const AtomPlacementArray& alist2 = bodies[b2].allAtoms;

                for (int j=0; j < (int)alist2.size(); ++j) {
                    const int       a2num = alist2[j].atomId;
                    assert(a2num != a1num);
                    const Atom&     a2 = atoms[a2num];
                    const ChargedAtomType& a2type  = chargedAtomTypes[a2.chargedAtomTypeId];
                    int                    a2cnum  = a2type.atomClass;
                    const AtomClass&       a2class = atomClasses[a2cnum];
                    
                    const Vec3  a2Station_G = X_GB2.R()*a2.station_B; // 15 flops
                    const Vec3  a2Pos_G     = X_GB2.T() + a2Station_G;  // 3 flops
                    const Vec3  r = a2Pos_G - a1Pos_G; // from a1 to a2 (3 flops)
                    const Real  d2 = r.normSqr(); // 5 flops

                    // Check for cutoffs on d2?

                    const Real  ood = 1/std::sqrt(d2); // approx 40 flops
                    const Real  ood2 = ood*ood;

                    // Coulomb. This unfortunately needs the separation distance which
                    // is expensive. But if scale, q1, or q2 are zero we can skip that.

                    Real eCoulomb = 0, fCoulomb = 0;
                    const Real qq = coulombScale[a2num]*q1Fac*a2type.partialCharge; // 2 flops
                    eCoulomb = qq * ood; //  scale*(1/(4*pi*e0)) *  q1*q2/d       (1 flop)  
                    fCoulomb = eCoulomb; // -scale*(1/(4*pi*e0)) * -q1*q2/d^2 * d (factor of 1/d^2 missing)

                    // van der Waals.

                    // Get precomputed mixed dmin and emin. Must ask the lower-numbered atom class.
                    const Real dij = (a1cnum <= a2cnum ? a1class.vdwDij[a2cnum-a1cnum]
                                                       : a2class.vdwDij[a1cnum-a2cnum]);
                    const Real eij = (a1cnum <= a2cnum ? a1class.vdwEij[a2cnum-a1cnum]
                                                       : a2class.vdwEij[a1cnum-a2cnum]);

                    const Real ddij2  = dij*dij*ood2;        // (dmin_ij/d)^2 (2 flops)
                    const Real ddij6  = ddij2*ddij2*ddij2;   // 2 flops
                    const Real ddij12 = ddij6*ddij6;         // 1 flop

                    const Real eijScale = vdwScale[a2num]*eij;            // 1 flop
                    const Real eVdw =      eijScale * (ddij12 - 2*ddij6); // 3 flops
                    const Real fVdw = 12 * eijScale * (ddij12 - ddij6);   // factor of 1/d^2 missing (3 flops)
                    const Vec3 fj = ((fCoulomb+fVdw)*ood2) * r;      // to apply to atom j on b2 (5 flops)

                    pe += (eCoulomb + eVdw); // Da-A^2/ps^2  (2 flops)
                    rigidBodyForces[b2] += SpatialVec( a2Station_G % fj, fj);   // 15 flops
                    rigidBodyForces[b1] -= SpatialVec( a1Station_G % fj, fj);   // 15 flops
                }
            }
            unscaleBondedAtoms(a1,vdwScale,coulombScale);
        }
    }
}


// We scale short range interactions but only for bonds which cross bodies.
void DuMMForceFieldSubsystemRep::scaleBondedAtoms
   (const Atom& a, Vector& vdwScale, Vector& coulombScale) const 
{
    for (int i=0; i < (int)a.xbond12.size(); ++i) {
        const int ix = a.xbond12[i]; 
        vdwScale[ix]=vdwScale12; coulombScale[ix]=coulombScale12;
    }
    for (int i=0; i < (int)a.xbond13.size(); ++i) {
        const int ix = a.xbond13[i][1]; // the 2nd atom is the 1-3
        vdwScale[ix]=vdwScale13; coulombScale[ix]=coulombScale13;
    }
    if (vdwScale14 != 1 || coulombScale14 != 1)
        for (int i=0; i < (int)a.xbond14.size(); ++i) {
            const int ix = a.xbond14[i][2]; // the 3rd atom is the 1-4
            vdwScale[ix]=vdwScale14; coulombScale[ix]=coulombScale14;
        }
    if (vdwScale15 != 1 || coulombScale15 != 1)
        for (int i=0; i < (int)a.xbond15.size(); ++i) {
            const int ix = a.xbond15[i][3]; // the 4th atom is the 1-5
            vdwScale[ix]=vdwScale15; coulombScale[ix]=coulombScale15;
        }
}

void DuMMForceFieldSubsystemRep::unscaleBondedAtoms
   (const Atom& a, Vector& vdwScale, Vector& coulombScale) const 
{
    for (int i=0; i < (int)a.xbond12.size(); ++i) {
        const int ix = a.xbond12[i];    vdwScale[ix]=coulombScale[ix]=1;
    }
    for (int i=0; i < (int)a.xbond13.size(); ++i) {
        const int ix = a.xbond13[i][1]; vdwScale[ix]=coulombScale[ix]=1;
    }
    if (vdwScale14 != 1 || coulombScale14 != 1)
        for (int i=0; i < (int)a.xbond14.size(); ++i) {
            const int ix = a.xbond14[i][2]; vdwScale[ix]=coulombScale[ix]=1;
        }
    if (vdwScale15 != 1 || coulombScale15 != 1)
        for (int i=0; i < (int)a.xbond15.size(); ++i) {
            const int ix = a.xbond15[i][3]; vdwScale[ix]=coulombScale[ix]=1;
        }
}

void DuMMForceFieldSubsystemRep::loadElements() {
    elements.resize(93); // Room for 1-92. I guess that's a little ambitious!
    elements[1] =Element(1,  "H",  "Hydrogen",    1.008).setDefaultColor(Green);
    elements[2] =Element(2,  "He", "Helium",      4.003);
    elements[3] =Element(3,  "Li", "Lithium",     6.941);
    elements[6] =Element(6,  "C",  "Carbon",     12.011).setDefaultColor(Gray);
    elements[7] =Element(7,  "N",  "Nitrogen",   14.007).setDefaultColor(Blue);
    elements[8] =Element(8,  "O",  "Oxygen",     15.999).setDefaultColor(Red);
    elements[9] =Element(9,  "F",  "Fluorine",   18.998);
    elements[10]=Element(10, "Ne", "Neon",       20.180);
    elements[11]=Element(11, "Na", "Sodium",     22.990);
    elements[12]=Element(12, "Mg", "Magnesium",  24.305);
    elements[14]=Element(14, "Si", "Silicon",    28.086);
    elements[15]=Element(15, "P",  "Phosphorus", 30.974).setDefaultColor(Magenta);
    elements[16]=Element(16, "S",  "Sulphur",    32.066).setDefaultColor(Yellow);
    elements[17]=Element(17, "Cl", "Chlorine",   35.453);
    elements[18]=Element(18, "Ar", "Argon",      39.948);
    elements[19]=Element(19, "K",  "Potassium",  39.098);
    elements[20]=Element(20, "Ca", "Calcium",    40.078);
    elements[26]=Element(26, "Fe", "Iron",       55.845);
    elements[29]=Element(29, "Cu", "Copper",     63.546);
    elements[30]=Element(30, "Zn", "Zinc",       65.390);
    elements[36]=Element(36, "Kr", "Krypton",    83.800);
    elements[47]=Element(47, "Ag", "Silver",    107.868);
    elements[53]=Element(53, "I",  "Iodine",    126.904);
    elements[54]=Element(54, "Xe", "Xenon",     131.290);
    elements[79]=Element(79, "Au", "Gold",      196.967).setDefaultColor(Yellow);
    elements[92]=Element(92, "U",  "Uranium",   238.029);
}

void DuMMForceFieldSubsystemRep::dump() const 
{
    printf("Dump of DuMMForceFieldSubsystem:\n");
    printf("  NBodies=%d NAtoms=%d NAtomClasses=%d NChargedAtomTypes=%d NBonds=%d\n",
        bodies.size(), atoms.size(), atomClasses.size(), chargedAtomTypes.size(), bonds.size());
    for (int i=0; i < (int)bodies.size(); ++i) {
        printf("  Body %d:\n", i);
        bodies[i].dump();
    }
    for (int i=0; i < (int)atoms.size(); ++i) {
        printf("  Atom %d: ", i);
        atoms[i].dump();
    }
    for (int i=0; i < (int)atomClasses.size(); ++i) {
        if (!atomClasses[i].isValid()) continue;
        printf("  AtomClass %d:\n", i);
        atomClasses[i].dump();
    }
    for (int i=0; i < (int)chargedAtomTypes.size(); ++i) {
        if (!chargedAtomTypes[i].isValid()) continue;
        printf("  ChargedAtomType %d:\n", i);
        chargedAtomTypes[i].dump();
    }
}

    //////////////
    // BondBend //
    //////////////

// Given a central atom location c bonded to atoms at r and s,
// calculate the angle between them, the potential energy,
// and forces on each of the three atoms.
void BondBend::harmonic
   (const Vec3& cG, const Vec3& rG, const Vec3& sG,
    Real& theta, Real& pe, Vec3& cf, Vec3& rf, Vec3& sf) const
{
    const Vec3 r = rG - cG; //               3 flops
    const Vec3 s = sG - cG; //               3 flops
    const Real rr = ~r*r, ss = ~s*s;    // |r|^2, |s|^2 ( 10 flops)

    const Real rs = ~r * s; // r dot s      (5 flops)
    const Vec3 rxs = r % s; // r cross s    (9 flops)
    const Real rxslen = rxs.norm(); //      (~35 flops)
    theta = std::atan2(rxslen, rs); //       ~50 flops
    const Real bend = theta - theta0;   //   1 flop
    pe = k*bend*bend; // NOTE: no factor of 1/2 (2 flops)

    // p is unit vector perpendicular to r and s

    // TODO: come up with something for when rxslen is 0 (vectors r & s
    // aligned or opposite); for relaxation
    // just needs to push them apart; what to do for dynamics?
    // Here we'll just make up a direction perpendicular to both
    // vectors and use it.
    const UnitVec3 p = (rxslen != 0 ? UnitVec3(rxs/rxslen,true)  // ~11 flops
                                    : UnitVec3(r).perp()); 
    const Real ffac = -2*k*bend; // 2 flops
    rf = (ffac/rr)*(r % p);          // ~20 flops
    sf = (ffac/ss)*(p % s);          // ~20 flops
    cf = -(rf+sf); // makes the net force zero (6 flops)
}

    /////////////////
    // BondTorsion //
    /////////////////

// Given atom locations r-x-y-s in the ground frame, calculate the
// torsion angle, energy and a force on each atom so that the desired
// pure torque is produced.
// This code is modeled in part after Tinker's torsion code in
// etors1.f because I couldn't figure out how to do it myself
// (sherm 060905). Thanks, Jay!
void BondTorsion::periodic(const Vec3& rG, const Vec3& xG, const Vec3& yG, const Vec3& sG,
              Real& theta, Real& pe, 
              Vec3& rf, Vec3& xf, Vec3& yf, Vec3& sf) const
{
    // All vectors point along the r->x->y->s direction
    const Vec3 r  = xG - rG; //               3 flops
    const Vec3 s  = sG - yG; //               3 flops
    const Vec3 xy = yG - xG; //               3 flops

    // Create a unit vector v along the axis, using increasingly
    // desperate measures in case of overlapping atoms. If we
    // don't have a real axis (i.e., atoms x and y overlap)
    // we'll signal that with oov==0 (see below). We don't care
    // much what happens in that case, but we hope to do something
    // remotely plausible so a stuck minimization will have some
    // hope of getting unstuck.

    const Real vv = ~xy*xy;                     //   5 flops
    const Real oov = (vv==0 ? Real(0) 
                            : 1/std::sqrt(vv)); // ~40 flops
    const UnitVec3 v = 
        (oov != 0 ? UnitVec3(xy*oov,true)       //   4 flops
                   : ((r%s).norm() != 0 ? UnitVec3(r % s)
                                        : UnitVec3(r).perp()));

    // Calculate plane normals. Axis vector v serves as the "x" 
    // axis of both planes. Vectors r (r->x) and s (y->s) are in
    // the plane in a vaguely "y axis" way, so t=rXv is the "z" axis
    // (plane normal) for the first plane and u=vXs is the plane normal
    // for the second. When those normals are aligned theta is 0.
    const Vec3 t = r % v, u = v % s; // 18 flops

    // If either r or s are aligned with the axis, we can't generate
    // a torque so we're done.
    const Real tt = ~t*t, uu = ~u*u; // 10 flops
    if (tt == 0 || uu == 0) {
        pe = 0; rf=xf=yf=sf=Vec3(0);
        return;
    }

    const Vec3 txu = t % u;                 //   9 flops
    const Real ootu = 1/std::sqrt(tt*uu);   // ~40 flops
    const Real cth = (~t*u)*ootu;           //   6 flops
    const Real sth = (~v*txu)*ootu;         //   6 flops
    theta = std::atan2(sth,cth);            // ~50 flops

    Real torque = 0;
    pe = 0; 
    for (int i=0; i < (int)terms.size(); ++i) {
        pe += terms[i].energy(theta);
        torque += terms[i].torque(theta);
    }

    const Vec3 ry = yG-rG;    // from r->y        3 flops
    const Vec3 xs = sG-xG;    // from x->s        3 flops
    const Vec3 dedt =  (torque/tt)*(t % v);  // ~20 flops
    const Vec3 dedu = -(torque/uu)*(u % v);  // ~21 flops

    rf = dedt % v; // 9 flops
    sf = dedu % v; // 9 flops
    if (oov==0) {
        xf = -rf;   // No axis; this is just desperation.
        yf = -sf;   // At least it keeps the forces summing to 0.
    } else {
        xf = ((ry % dedt) + (dedu % s))*oov;
        yf = ((dedt % r) + (xs % dedu))*oov;
    }
}
    //////////
    // Atom //
    //////////

void Atom::dump() const {
    printf(" chargedAtomType=%d body=%d station=%g %g %g\n",
        chargedAtomTypeId, bodyId, station_B[0], station_B[1], station_B[2]);

    printf("    bond 1-2:");
    for (int i=0; i < (int)bond12.size(); ++i)
        printf(" %d", bond12[i]);
    printf("\n    bond 1-3:");
    for (int i=0; i < (int)bond13.size(); ++i)
        printf(" %d-%d", bond13[i][0], bond13[i][1]);
    printf("\n    bond 1-4:");
    for (int i=0; i < (int)bond14.size(); ++i)
        printf(" %d-%d-%d", bond14[i][0], bond14[i][1], bond14[i][2]);
    printf("\n    bond 1-5:");
    for (int i=0; i < (int)bond15.size(); ++i)
        printf(" %d-%d-%d-%d", bond15[i][0], bond15[i][1], bond15[i][2], bond15[i][3]);
    printf("\n");

    printf("    xbond 1-2:");
    for (int i=0; i < (int)xbond12.size(); ++i)
        printf(" %d", xbond12[i]);
    printf("\n    xbond 1-3:");
    for (int i=0; i < (int)xbond13.size(); ++i)
        printf(" %d-%d", xbond13[i][0], xbond13[i][1]);
    printf("\n    xbond 1-4:");
    for (int i=0; i < (int)xbond14.size(); ++i)
        printf(" %d-%d-%d", xbond14[i][0], xbond14[i][1], xbond14[i][2]);
    printf("\n    xbond 1-5:");
    for (int i=0; i < (int)xbond15.size(); ++i)
        printf(" %d-%d-%d-%d", xbond15[i][0], xbond15[i][1], xbond15[i][2], xbond15[i][3]);
    printf("\n");

    printf("    1-2 stretch:");
    for (int i=0; i < (int)stretch.size(); ++i)
        printf(" (%g,%g)", stretch[i].k, stretch[i].d0);
    printf("\n    1-3 bend:");
    for (int i=0; i < (int)bend.size(); ++i)
        printf(" (%g,%g)", bend[i].k, bend[i].theta0);
    printf("\n    1-4 torsion:\n");
    for (int i=0; i < (int)torsion.size(); ++i) {
        const BondTorsion& bt = torsion[i];
        printf("     ");
        for (int j=0; j<(int)bt.terms.size(); ++j) {
            const TorsionTerm& tt = bt.terms[j];
            printf(" (%d:%g,%g)", tt.periodicity, 
                                  tt.amplitude, tt.theta0);
        }
        printf("\n");
    }
    printf("\n");
}

    ////////////////
    // RigidGroup //
    ////////////////

// True if the atom has been placed in this group, or if any group that has
// been placed here contains the atom, or a subgroup containing the atom
// recursively.
bool RigidGroup::containsAtom(int atomId, const DuMMForceFieldSubsystemRep& mm) const {
    if (atomPlacements.find(AtomPlacement(atomId,Vec3(0))) != atomPlacements.end())
        return true;
    RigidGroupPlacementSet::const_iterator gp = rigidGroupPlacements.begin();
    while (gp != rigidGroupPlacements.end()) {
        const RigidGroup& g = mm.getRigidGroup(gp->rigidGroupId);
        if (g.containsAtom(atomId, mm))
            return true;
        ++gp;
    }
    return false;
}

bool RigidGroup::containsRigidGroup(int groupId, const DuMMForceFieldSubsystemRep& mm) const {
    if (rigidGroupPlacements.find(RigidGroupPlacement(groupId,Transform())) != rigidGroupPlacements.end())
        return true;
    RigidGroupPlacementSet::const_iterator gp = rigidGroupPlacements.begin();
    while (gp != rigidGroupPlacements.end()) {
        const RigidGroup& g = mm.getRigidGroup(gp->rigidGroupId);
        if (g.containsRigidGroup(groupId, mm))
            return true;
        ++gp;
    }
    return false;
}

// Append all our directly-placed atoms into the allAtoms array. Then transform
// all of our subgroups' atoms by the appropriate transform and append them to
// the allAtoms array. The results are not in any particular order.
void RigidGroup::calculateAllAtomArray(DuMMForceFieldSubsystemRep& mm)
{
    AtomPlacementSet::const_iterator ap = atomPlacements.begin();
    while(ap != atomPlacements.end())
        allAtoms.push_back(*ap);

    RigidGroupPlacementSet::iterator gp = rigidGroupPlacements.begin();
    while (gp != rigidGroupPlacements.end()) {
        RigidGroup& g = mm.updRigidGroup(gp->rigidGroupId);
        g.realizeTopologicalCache(mm);  // might already be done
        g.xformAllAtomArray(gp->placement, allAtoms);
        ++gp;
    }
}

} // namespace SimTK

