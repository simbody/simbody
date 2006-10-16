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
 * Private implementation of DuMMForceFieldSubsystem. Units here are uniformly
 * MD units: nanometers, daltons, picoseconds, with energy in kilojoules/mole.
 * We accept angles from users in degrees, but use only radians internally.
 */

#include "Simbody.h"
#include "simbody/internal/ForceSubsystem.h"
#include "simbody/internal/DuMMForceFieldSubsystem.h"

#include "simbody/internal/DecorativeGeometry.h"

#include "ForceSubsystemRep.h"

#include <string>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cassert>
#include <set>
#include <map>
#include <algorithm>

namespace SimTK {

// Conversion constants (multiply something in units on left to get
// equivalent in units on right).

/*static*/ const Real DuMMForceFieldSubsystem::Ang2Nm  = Real(0.1L);
/*static*/ const Real DuMMForceFieldSubsystem::Nm2Ang  = Real(10);
/*static*/ const Real DuMMForceFieldSubsystem::Kcal2KJ = (Real)SimTK_KCAL_TO_KJOULE;
/*static*/ const Real DuMMForceFieldSubsystem::KJ2Kcal = (Real)SimTK_KJOULE_TO_KCAL;
/*static*/ const Real DuMMForceFieldSubsystem::Deg2Rad = (Real)SimTK_DEGREE_TO_RADIAN;
/*static*/ const Real DuMMForceFieldSubsystem::Rad2Deg = (Real)SimTK_RADIAN_TO_DEGREE;
/*static*/ const Real DuMMForceFieldSubsystem::Sigma2Radius = 
                        (Real)std::pow(2.L,  1.L/6.L); // sigma < radius
/*static*/ const Real DuMMForceFieldSubsystem::Radius2Sigma = 
                        (Real)std::pow(2.L, -1.L/6.L);

// handy abbreviations
static const Real Pi      = (Real)SimTK_PI;
static const Real& Deg2Rad = DuMMForceFieldSubsystem::Deg2Rad;
static const Real& Rad2Deg = DuMMForceFieldSubsystem::Rad2Deg;
static const Real& KJ2Kcal = DuMMForceFieldSubsystem::KJ2Kcal;
static const Real& Kcal2KJ = DuMMForceFieldSubsystem::Kcal2KJ;

// This is Coulomb's constant 1/(4*pi*e0) in units which convert
// e^2/nm to kJ/mol.

static const Real CoulombFac = (Real)SimTK_COULOMB_CONSTANT_IN_MD;

//static const Real CoulombFac = 332.06371 * EnergyUnitsPerKcal;

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

static const Real oo6  = Real(1/6.L);
static const Real oo13 = Real(1/13.L);

// This doesn't seem to be used by anyone but it should be!
// Ref: Waldman, M. & Hagler, A.T. New combining rules for
// rare gas van der Waals parameters. 
// J. Comput. Chem. 14(9):1077 (1993).
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
    Real mass;         // in daltons (Da, g/mol, amu, u)
    Vec3 defaultColor;
    std::string symbol;
    std::string name;
};

class AtomClass {
public:
    AtomClass() : element(-1), valence(-1), vdwRadius(-1), vdwWellDepth(-1) { }
    AtomClass(int id, const char* nm, int e, int v, Real radInNm, Real wellDepthInKJ)
      : atomClassId(id), name(nm), element(e), valence(v), 
        vdwRadius(radInNm), vdwWellDepth(wellDepthInKJ)
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
        printf("   %d(%s): element=%d, valence=%d vdwRad=%g nm, vdwDepth(kJ)=%g (%g kcal)\n",
            atomClassId, name.c_str(), element, valence, vdwRadius, vdwWellDepth,
            vdwWellDepth*KJ2Kcal);
        printf("    vdwDij (nm):");
        for (int i=0; i< (int)vdwDij.size(); ++i)
            printf(" %g", vdwDij[i]);
        printf("\n    vdwEij (kJ):");
        for (int i=0; i< (int)vdwEij.size(); ++i)
            printf(" %g", vdwEij[i]);
        printf("\n");
    }

        // TOPOLOGICAL STATE VARIABLES
        //   Filled in during construction.

    int         atomClassId;
    std::string name;

    int     element;
    int     valence;       // # of direct bonds expected
    Real    vdwRadius;     // ri, nm
    Real    vdwWellDepth;  // ei, kJ=Da-nm^2/ps^2


        // TOPOLOGICAL CACHE ENTRIES
        //   These are calculated in realizeTopology() from topological
        //   state variables (from here or others in the DuMM class).

    // After all types have been defined, we can calculate vdw 
    // combining rules for dmin and well depth energy. We only fill
    // in entries for pairings of this class with itself and with
    // higher-numbered atom types, so to find the entry for class c, 
    // index these arrays by c-atomClassId where atomClassId is the
    // class Id of the present AtomClass.
    // Note that different combining rules may be used but they
    // will always result in a pair of vdw parameters.
    std::vector<Real> vdwDij;   // nm
    std::vector<Real> vdwEij;   // kJ=Da-A^2/ps^2
};

class ChargedAtomType {
public:
    ChargedAtomType() : chargedAtomTypeId(-1), atomClassId(-1), partialCharge(NTraits<Real>::NaN) { }
    ChargedAtomType(int id, const char* nm, int aclass, Real chg)
      : chargedAtomTypeId(id), name(nm), atomClassId(aclass), partialCharge(chg) 
    { 
        assert(isValid());
    }
    bool isValid() const {return chargedAtomTypeId >= 0 && atomClassId >= 0;}

    void dump() const {
        printf("    %d(%s): atomClassId=%d, chg=%g e\n", 
               chargedAtomTypeId, name.c_str(), atomClassId, partialCharge);
    }

    // These are all Topological state variables, filled in during construction.
    // There are no calculations to be performed.
    int         chargedAtomTypeId;
    std::string name;

    int         atomClassId;
    Real        partialCharge; // qi, in e (charge on proton)

};

// This represents bond-stretch information for a pair of atom types.
// Use an IntPair as a key.
class BondStretch {
public:
    BondStretch() : k(-1), d0(-1) { }
    BondStretch(Real stiffnessInKJperNmSq, Real lengthInNm) 
      : k(stiffnessInKJperNmSq), d0(lengthInNm) { 
        assert(isValid());
    }
    bool isValid() const {return k >= 0 && d0 >= 0; }
    Real k;  // in energy units (kJ=Da-nm^2/ps^2) per nm^2, i.e. Da/ps^2
    Real d0; // distance at which force is 0 (in nm)
};

class BondBend {
public:
    BondBend() : k(-1), theta0(-1) { }
    BondBend(Real stiffnessInKJPerRadSq, Real angleInDeg) 
      : k(stiffnessInKJPerRadSq), theta0(angleInDeg*Deg2Rad) {
        assert(isValid());
    }
    bool isValid() const {return k >= 0 && (0 <= theta0 && theta0 <= Pi);}

    // Given a central atom location c bonded to atoms at r and s,
    // calculate the angle between them, the potential energy,
    // and forces on each of the three atoms.
    void harmonic(const Vec3& cG, const Vec3& rG, const Vec3& sG, const Real& scale,
                  Real& theta, Real& pe, Vec3& cf, Vec3& rf, Vec3& sf) const;

    Real k;      // energy units kJ per rad^2, i.e. Da-nm^2/(ps^2-rad^2)
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
    TorsionTerm(int n, Real ampInKJ, Real th0InDeg) 
      : periodicity(n), amplitude(ampInKJ), theta0(th0InDeg*Deg2Rad) {
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
    Real amplitude; // energy units (kJ)
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
                  const Real& scale, Real& theta, Real& pe, 
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
    Vec3 station;   // in nm
};
inline bool operator<(const AtomPlacement& a1, const AtomPlacement& a2) {
    return a1.atomId < a2.atomId;
}
inline bool operator==(const AtomPlacement& a1, const AtomPlacement& a2) {
    return a1.atomId == a2.atomId;
}

class ClusterPlacement {
public:
    ClusterPlacement() : clusterId(-1) { }
    ClusterPlacement(int c, const Transform& t) : clusterId(c), placement(t) {
        assert(isValid());
    }
    bool isValid() const {return clusterId >= 0;}

    int         clusterId;
    Transform   placement;  // translation in nm
};
inline bool operator<(const ClusterPlacement& r1, const ClusterPlacement& r2) {
    return r1.clusterId < r2.clusterId;
}
inline bool operator==(const ClusterPlacement& r1, const ClusterPlacement& r2) {
    return r1.clusterId == r2.clusterId;
}

typedef std::vector<int>            AtomArray;
typedef std::vector<AtomPlacement>  AtomPlacementArray;
typedef std::set<AtomPlacement>     AtomPlacementSet;
typedef std::set<ClusterPlacement>  ClusterPlacementSet;

class Atom {
public:
    Atom() 
      : atomId(-1), chargedAtomTypeId(-1), bodyId(-1) {
    }
    Atom(int t, int aId) : atomId(aId), chargedAtomTypeId(t), bodyId(-1) {
        assert(isValid());
    }
    bool isValid() const {return atomId>=0 && chargedAtomTypeId>=0;}
    bool isAttachedToBody() const {return bodyId >= 0;}

    int getBodyId() const {assert(isAttachedToBody()); return bodyId;}

    void attachToBody(int bnum, const Vec3& s) {
        assert(!isAttachedToBody());
        bodyId = bnum;
        station_B = s;
    }

    bool isBondedTo(int anum) const {
        for (int i=0; i<(int)bond12.size(); ++i)
            if (bond12[i] == anum) return true;
        return false;
    }

    void dump() const;

    void invalidateTopologicalCache() {
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

    // After the atom or a containing cluster has been attached to a
    // body, we fill these in.
    int  bodyId;
    Vec3 station_B; // atom's station fixed in body bodyId's frame, in nm

        // TOPOLOGICAL CACHE ENTRIES
        //   These are calculated in realizeTopology() from topological
        //   state variables (from here or others in the DuMM class).

    // This is a group of lists which identify atoms nearby in the
    // molecule's bond structure. The bond12 list above contains the directly
    // bonded (1-2) atoms; the 13 list below has the 1-(2)-3 bonded atoms (that
    // is, it includes the path to the "3" atom), etc. The current Atom is
    // always atom "1" so it isn't stored.

    std::vector<IntPair>   bond13;
    std::vector<IntTriple> bond14;
    std::vector<IntQuad>   bond15;

    // These are shorter versions of the bond lists in which only those
    // bonds which include atoms from at least two bodies are included.
    // Note that each bond will appear twice in the overall data structure,
    // in the Atom entries for the atoms at either end. We avoid double
    // processing by only processing the instance in which the first atoms'
    // Id is the lower of the two. But we need to keep both copies because
    // these are also used for scaling nearby interaction during non-bonded 
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
    Real     netCharge;         // in proton charge units e
    Vec3     centerOfCharge;    // in nm
    Vec3     dipoleMoment;      // units?? TODO
    SymMat33 quadrupoleMoment;  // units?? TODO
};

class GeometricProperties {
public:
    Transform obbFrame;
    Vec3      obbHalfLengths;       // nm
    Real      boundingSphereRadius; // nm
    Vec3      boundingSphereCenter; // nm
};

//
// This class is a rigid grouping of atoms. It can contain atoms directly, and
// subclusters which can contain atoms or sub-subclusters, etc. As we build
// up a cluster, we keep a running "flat" view of all the atoms and all the clusters
// contained anywhere deep within, already transformed to this cluster's reference
// frame.
//
class Cluster {
public:
    Cluster() : clusterId(-1), topologicalCacheValid(false), bodyId(-1) { }
    Cluster(const char* nm)
      : clusterId(-1), name(nm), topologicalCacheValid(false), bodyId(-1) {
        // not valid yet -- still need Id assigned
    }

    bool isValid() const {return clusterId >= 0;}
    bool isAttachedToBody() const {return bodyId >= 0;}
    bool isTopLevelCluster() const {return parentClusters.empty();}

    bool isTopologicalCacheValid() const   {return topologicalCacheValid;}
    void invalidateTopologicalCache()      {topologicalCacheValid=false;}

    int getBodyId() const {assert(isAttachedToBody()); return bodyId;}

    const AtomPlacementSet& getDirectlyContainedAtoms() const {return directAtomPlacements;}
    const AtomPlacementSet& getAllContainedAtoms()      const {return allAtomPlacements;}
    AtomPlacementSet&       updAllContainedAtoms()            {return allAtomPlacements;}

    const ClusterPlacementSet& getDirectlyContainedClusters() const {return directClusterPlacements;}
    const ClusterPlacementSet& getAllContainedClusters()      const {return allClusterPlacements;}
    ClusterPlacementSet&       updAllContainedClusters()            {return allClusterPlacements;}

    bool containsAtom(int atomId) const {
        return allAtomPlacements.find(AtomPlacement(atomId,Vec3(0))) 
                != allAtomPlacements.end();
    }
    bool containsCluster(int clusterId) const {
        return allClusterPlacements.find(ClusterPlacement(clusterId,Transform())) 
                != allClusterPlacements.end();
    }

    // See if a cluster contains any atoms which are already in
    // any of the cluster trees to which this cluster is associated.
    // TODO: can only handle top-level cluster so we don't have to run up the
    //       ancestor branches.
    // If we find an atom common to both clusters we'll return it to permit
    // nice error messages, otherwise we return false and -1 for the atomId.
    bool overlapsWithCluster(const Cluster& test, int& anAtomIdInBothClusters) const {
        assert(isTopLevelCluster());

        const AtomPlacementSet& testAtoms = test.getAllContainedAtoms();
        const AtomPlacementSet& myAtoms   = getAllContainedAtoms();

        AtomPlacementSet::const_iterator ap = testAtoms.begin();
        while (ap != testAtoms.end()) {
            if (containsAtom(ap->atomId)) {
                anAtomIdInBothClusters = ap->atomId;
                return true;
            }
            ++ap;
        }
        anAtomIdInBothClusters = -1;
        return false;
    }

    // Return true if this cluster contains (directly or indirectly) any atom which has already
    // been attached to a body. If so return one of the attached atoms and its body, which can
    // be helpful in error messages.
    bool containsAnyAtomsAttachedToABody(int& atomId, int& bodyId, 
                                         const DuMMForceFieldSubsystemRep& mm) const;

    // Translation is in nm.
    void attachToBody(int bnum, const Transform& X_BR, DuMMForceFieldSubsystemRep& mm);

    // Place an atom in this cluster. To be valid, the atom must not
    // already be
    //   (a) in any of the trees of which this group is a part, or
    //   (b) attached to a body.
    // TODO: (c) at the moment we don't allow placing an atom in a group unless
    //           that group is a top-level group (i.e., it has no parents).
    // If this group is already attached to a body, then we will update
    // the atom entry to note that it is now attached to the body also.
    void placeAtom(int atomId, const Vec3& stationInNm, DuMMForceFieldSubsystemRep& mm);

    // Place a child cluster in this parent cluster. To be valid, the child 
    // must not 
    //   (a) already be contained in the parent group or one of the parent's subgroups, or
    //   (b) contain any atoms which are already present in the parent or any
    //       of the parent's subgroups, or
    //   (c) already be attached to a body.
    // TODO: (d) at the moment we don't allow adding a child group unless
    //           the parent (this) group is a top-level group (i.e., it has no parents).
    // If the parent is already attached to a body, then we will update
    // the child to note that it is now attached to the body also (and it
    // will update its contained atoms).
    // (translation is in nm)
    void placeCluster(int childClusterId, const Transform& placement, DuMMForceFieldSubsystemRep& mm);


    // Calculate the composite mass properties for this cluster, transformed into
    // the indicated frame. Translation part of the Transform is in nm, returned mass
    // proprties are in daltons and nm.
    MassProperties calcMassProperties
       (const Transform& tr, const DuMMForceFieldSubsystemRep& mm) const;


    // Recursively calculate composite properties for this group and all the
    // groups it contains. All groups were marked "invalid" at the beginning
    // of this step.
    void realizeTopologicalCache(DuMMForceFieldSubsystemRep& mm) {
        if (topologicalCacheValid)
            return;

        // TODO: charge & geometric properties

        topologicalCacheValid = true;
    }


    void dump() const {
        printf("    clusterId=%d(%s)\n", clusterId, name.c_str());
        printf("      direct atom placements (nm): ");
        AtomPlacementSet::const_iterator ap = directAtomPlacements.begin();
        while (ap != directAtomPlacements.end()) {
            std::cout << " " << ap->atomId << ":" << ap->station;
            ++ap;
        }
        printf("\n      all atom placements (nm): ");
        AtomPlacementSet::const_iterator aap = allAtomPlacements.begin();
        while (aap != allAtomPlacements.end()) {
            std::cout << " " << aap->atomId << ":" << aap->station;
            ++aap;
        }
        printf("\n      direct cluster placements (nm):\n");
        ClusterPlacementSet::const_iterator cp = directClusterPlacements.begin();
        while (cp != directClusterPlacements.end()) {
            std::cout << "      " << cp->clusterId << ":" << cp->placement;
            ++cp;
        }
        printf("\n      all cluster placements (nm):\n");
        ClusterPlacementSet::const_iterator acp = allClusterPlacements.begin();
        while (acp != allClusterPlacements.end()) {
            std::cout << "      " << acp->clusterId << ":" << acp->placement;
            ++acp;
        }
        printf("\n      parent cluster placements (nm):\n");
        ClusterPlacementSet::const_iterator pp = parentClusters.begin();
        while (pp != parentClusters.end()) {
            std::cout << "      " << pp->clusterId << ":" << pp->placement;
            ++pp;
        }

        std::cout     << "\n      topological cache valid? " 
                  << isTopologicalCacheValid() << std::endl;
        if (bodyId >= 0) 
            std::cout << "\n      attached to body " << bodyId << " at (nm) " << placement_B;
        else
            std::cout << "\n      NOT ATTACHED TO ANY BODY.";
        std::cout << std::endl;
    }

    void clearAllCalculatedData() {
        topologicalCacheValid = false;
        chargeProps    = ChargeProperties();
        geometricProps = GeometricProperties();
    }

private:
    // translation is in nm
    void noteNewChildCluster(int childClusterId, const Transform& X_PC) {
        std::pair<ClusterPlacementSet::iterator, bool> ret;
        ret = directClusterPlacements.insert(ClusterPlacement(childClusterId,X_PC));
        assert(ret.second); // must not have been there already

        ret = allClusterPlacements.insert(ClusterPlacement(childClusterId,X_PC));
        assert(ret.second); // must not have been there already
    }

    // translation is in nm
    void noteNewParentCluster(int parentClusterId, const Transform& X_PC) {
        std::pair<ClusterPlacementSet::iterator, bool> ret =
            parentClusters.insert(ClusterPlacement(parentClusterId,X_PC));
        assert(ret.second); // must not have been there already
    }

public:
        // TOPOLOGICAL STATE VARIABLES
        //   Filled in during construction.
    int                 clusterId;
    std::string         name;

    // These are the *directly* attached atoms and clusters.
    AtomPlacementSet    directAtomPlacements;
    ClusterPlacementSet directClusterPlacements;

    // These sets are kept up to date as we add atoms and clusters.
    // 'allAtomPlacements' contains *all* the atoms in this cluster
    // or its descendents, transformed into this cluster's frame.
    // 'allClusterPlacements' contains *all* the clusters in this
    // cluster or its subclusters, transformed into this cluster's frame.
    AtomPlacementSet    allAtomPlacements;
    ClusterPlacementSet allClusterPlacements;

    // This is a list of all the immediate parents of this cluster, if any.
    // This is updated whenever this cluster is placed in another one. The
    // body is *not* considered a parent cluster; it is handled separately
    // below. Note that whenever an atom or cluster is added to this cluster,
    // the atom or atoms involved [SHOULD BE: TODO] added to each ancestor.
    ClusterPlacementSet parentClusters;

    // After this cluster or a containing cluster has been attached to a
    // body, we can fill these in.
    int       bodyId;
    Transform placement_B; // cluster's placement fixed in body bodyId's frame (nm)

        // TOPOLOGICAL CACHE ENTRIES
        //   These are calculated in realizeTopology() from topological
        //   state variables (from here or others in the DuMM class).

    bool topologicalCacheValid;

    // These reflect composite properties built from the allAtoms list.
    ChargeProperties    chargeProps;
    GeometricProperties geometricProps;
};

// A Body has a reference to a top-level Cluster, plus some information used
// at runtime for fast body-by-body processing.
class Body {
public:
    Body() : clusterId(-1) { }
    explicit Body(int cId) : clusterId(cId) { 
        assert(isValid());
    }

    bool isValid() const {return clusterId >= 0;}

    void invalidateTopologicalCache() {allAtoms.clear();}
    void realizeTopologicalCache(const DuMMForceFieldSubsystemRep& mm);

    int getClusterId() const {assert(isValid()); return clusterId;}

    void dump() const {
        printf("    clusterId=%d\n", clusterId);
        printf("    shadowBodies=");
        for (int i=0; i < (int)shadowBodies.size(); ++i)
            printf(" %d", shadowBodies[i]);
        printf("\n");
        printf("    allAtoms=");
        for (int i=0; i < (int)allAtoms.size(); ++i) 
            printf(" %d(%g,%g,%g)(nm)", allAtoms[i].atomId,
                allAtoms[i].station[0], allAtoms[i].station[1], allAtoms[i].station[2]);
        printf("\n");
    }

    static std::string createClusterNameForBody(int bnum) {
        char buf[100];
        std::sprintf(buf, "Body %d", bnum);
        return std::string(buf);
    }

    int clusterId;
    std::vector<int> shadowBodies; // if needed

    // This is an expansion of all the atom & group placements, with
    // all stations transformed to this body's frame, sorted in order
    // of atomId, and built for speed!
    AtomPlacementArray  allAtoms;
};

class DuMMForceFieldSubsystemRep : public ForceSubsystemRep {
    friend class DuMMForceFieldSubsystem;
    static const char* ApiClassName; // "DuMMForceFieldSubsystem"
public:
    DuMMForceFieldSubsystemRep()
      : ForceSubsystemRep("DuMMForceFieldSubsystem", "0.0.1")
    {
        topologicalCacheValid = false;

        vdwMixingRule = DuMMForceFieldSubsystem::WaldmanHagler;
        vdwGlobalScaleFactor=coulombGlobalScaleFactor=bondStretchGlobalScaleFactor
            =bondBendGlobalScaleFactor=bondTorsionGlobalScaleFactor=1;
        vdwScale12=coulombScale12=vdwScale13=coulombScale13=0;
        vdwScale14=coulombScale14=vdwScale15=coulombScale15=1;
        loadElements();
        const int gid = addCluster(Cluster("free atoms and groups"));
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

    bool isValidCluster(int clusterId) const {
        return 0 <= clusterId && clusterId < (int)clusters.size()
                && clusters[clusterId].isValid();
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

    // Radii and returned diameter are given in nm, energies in kJ/mol.
    void applyMixingRule(Real ri, Real rj, Real ei, Real ej, Real& dmin, Real& emin) const
    {
        Real rmin;

        switch(vdwMixingRule) {
        case DuMMForceFieldSubsystem::WaldmanHagler:     
            vdwCombineWaldmanHagler(ri,rj,ei,ej,rmin,emin);     break;
        case DuMMForceFieldSubsystem::HalgrenHHG:         
            vdwCombineHalgrenHHG(ri,rj,ei,ej,rmin,emin);        break;
        case DuMMForceFieldSubsystem::Jorgensen:         
            vdwCombineJorgensen(ri,rj,ei,ej,rmin,emin);         break;
        case DuMMForceFieldSubsystem::LorentzBerthelot:  
            vdwCombineLorentzBerthelot(ri,rj,ei,ej,rmin,emin);  break;
        case DuMMForceFieldSubsystem::Kong:              
            vdwCombineKong(ri,rj,ei,ej,rmin,emin);              break;
        default: assert(!"unknown vdw mixing rule");
        };

        dmin = 2*rmin;
    }

    int addCluster(const Cluster& c) {
        const int clusterId = (int)clusters.size();
        clusters.push_back(c);
        clusters[clusterId].clusterId = clusterId;
        return clusterId;
    }
    Cluster& updCluster(int clusterId) {
        assert(isValidCluster(clusterId));
        return clusters[clusterId];
    }
    const Cluster& getCluster(int clusterId) const {
        assert(isValidCluster(clusterId));
        return clusters[clusterId];
    }
    Body& updBody(int bodyId) {
        assert(isValidBody(bodyId));
        return bodies[bodyId];
    }
    const Body& getBody(int bodyId) const {
        assert(isValidBody(bodyId));
        return bodies[bodyId];
    }


    int getNAtoms() const {return (int)atoms.size();}
    int getNBonds() const {return (int)bonds.size();}

    const Atom& getAtom(int atomId) const {
        assert(isValidAtom(atomId));
        return atoms[atomId];
    }
    Atom& updAtom(int atomId) {
        assert(isValidAtom(atomId));
        return atoms[atomId];
    }

    int getChargedAtomTypeId(int atomId) const {
        return getAtom(atomId).chargedAtomTypeId;
    }

    int getAtomClassId(int atomId) const {
        const ChargedAtomType& type = chargedAtomTypes[getChargedAtomTypeId(atomId)];
        return type.atomClassId;
    }

    int getAtomElementNum(int atomId) const {
        const AtomClass& cl = atomClasses[getAtomClassId(atomId)];
        return cl.element;
    }

    const Element& getElement(int element) const {
        assert(isValidElement(element));
        return elements[element];
    }


    const BondStretch& getBondStretch(int class1, int class2) const;
    const BondBend&    getBondBend   (int class1, int class2, int class3) const;
    const BondTorsion& getBondTorsion(int class1, int class2, int class3, int class4) const;

    void realizeTopology(State& s) const;

    void realizeModel(State& s) const {
        // Sorry, no choices available at the moment.
    }

    void realizeInstance(const State& s) const {
        // Nothing to compute here.
    }

    void realizeTime(const State& s) const {
        // Nothing to compute here.
    }

    void realizePosition(const State& s) const {
        // Nothing to compute here.
    }

    void realizeVelocity(const State& s) const {
        // Nothing to compute here.
    }


    void realizeDynamics(const State& s) const;

    void realizeAcceleration(const State& s) const {
        // Nothing to compute here.
    }

    void dump() const;

    DuMMForceFieldSubsystemRep* cloneSubsystemRep() const {
        return new DuMMForceFieldSubsystemRep(*this);
    }

private:
    void loadElements();

    void ensureBodyEntryExists(int bodyNum) {
        if (bodyNum >= (int)bodies.size())
            bodies.resize(bodyNum+1);
        if (!bodies[bodyNum].isValid()) {
            const int clusterId = 
                addCluster(Cluster(Body::createClusterNameForBody(bodyNum).c_str()));
            clusters[clusterId].attachToBody(bodyNum, Transform(), *this);
            bodies[bodyNum] = Body(clusterId);
        }
    }

    void invalidateAllTopologicalCacheEntries() {
        topologicalCacheValid = false;

        // If any of these objects are invalid, the invalidateTopologicalCache()
        // call does nothing (i.e., it doesn't blow up!).

        // molecule
        for (int i=0; i < (int)atoms.size(); ++i)
            atoms[i].invalidateTopologicalCache();
        for (int i=0; i < (int)clusters.size(); ++i)
            clusters[i].invalidateTopologicalCache();
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

    std::vector<Atom>    atoms;
    std::vector<Bond>    bonds;
    std::vector<Cluster> clusters;
    // This defines the partitioning of atoms onto the matter subsystem's bodies.
    // The indices here correspond to the body numbers. Only entries for bodies on
    // which our atoms have been attached will be valid.
    std::vector<Body>    bodies;

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

    // Which rule to use for combining van der Waals radii and energy well
    // depth for dissimilar atom classes.
    DuMMForceFieldSubsystem::VdwMixingRule  vdwMixingRule;

    // Scale factors for nonbonded forces when applied to
    // atoms which are near in the graph formed by the bonds.
    Real vdwScale12, coulombScale12;    // default 0,0
    Real vdwScale13, coulombScale13;    // default 0,0
    Real vdwScale14, coulombScale14;    // default 1,1
    Real vdwScale15, coulombScale15;    // default 1,1

    // Global scale factors for non-physical disabling or fiddling with
    // individual force field terms.
    Real vdwGlobalScaleFactor, coulombGlobalScaleFactor; 
    Real bondStretchGlobalScaleFactor, bondBendGlobalScaleFactor, 
         bondTorsionGlobalScaleFactor;

        // TOPOLOGICAL CACHE ENTRIES
        //   These are calculated in realizeTopology() from topological
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
   (int atomClassId, const char* atomClassName, int element, int valence, 
    Real vdwRadiusInNm, Real vdwWellDepthInKJPerMol)
{
    static const char* MethodName = "defineAtomClass";
    DuMMForceFieldSubsystemRep& mm = updRep();

        // Catch nonsense arguments.
    SimTK_APIARGCHECK1_ALWAYS(atomClassId >= 0, mm.ApiClassName, MethodName,
        "atom class Id %d invalid: must be nonnegative", atomClassId);
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidElement(element), mm.ApiClassName, MethodName,
        "element %d invalid: must be a valid atomic number and have an entry here",element);
    SimTK_APIARGCHECK1_ALWAYS(valence >= 0, mm.ApiClassName, MethodName, 
        "expected valence %d invalid: must be nonnegative", valence);
    SimTK_APIARGCHECK1_ALWAYS(vdwRadiusInNm >= 0, mm.ApiClassName, MethodName, 
        "van der Waals radius %g invalid: must be nonnegative", vdwRadiusInNm);
    SimTK_APIARGCHECK1_ALWAYS(vdwWellDepthInKJPerMol >= 0, mm.ApiClassName, MethodName, 
        "van der Waals energy well depth %g invalid: must be nonnegative", vdwWellDepthInKJPerMol);

        // Make sure there is a slot available for this atom class.
    if (atomClassId >= (int)mm.atomClasses.size())
        mm.atomClasses.resize(atomClassId+1);

        // Make sure this atom class hasn't already been defined.
    SimTK_APIARGCHECK2_ALWAYS(!mm.atomClasses[atomClassId].isValid(), mm.ApiClassName, MethodName, 
        "atom class Id %d is already in use for '%s'", atomClassId, 
        mm.atomClasses[atomClassId].name.c_str());

        // It's all good -- add the new atom class.
    mm.atomClasses[atomClassId] = AtomClass(atomClassId, atomClassName, element, valence, 
                                            vdwRadiusInNm, vdwWellDepthInKJPerMol);
}

void DuMMForceFieldSubsystem::defineChargedAtomType
   (int chargedAtomTypeId, const char* typeName, int atomClassId, Real partialChargeInE)
{
    static const char* MethodName = "defineChargedAtomType";
    DuMMForceFieldSubsystemRep& mm = updRep();

        // Check for nonsense arguments.
    SimTK_APIARGCHECK1_ALWAYS(chargedAtomTypeId >= 0, mm.ApiClassName, MethodName,
        "charged atom type Id %d invalid: must be nonnegative", chargedAtomTypeId);
    SimTK_APIARGCHECK1_ALWAYS(atomClassId >= 0, mm.ApiClassName, MethodName,
        "atom class Id %d invalid: must be nonnegative", atomClassId);
    // partialCharge is a signed quantity

        // Make sure the referenced atom class has already been defined.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtomClass(atomClassId), mm.ApiClassName, MethodName,
        "atom class %d is undefined", atomClassId);

        // Make sure there is a slot available for the new chargedAtomType.
    if (chargedAtomTypeId >= (int)mm.chargedAtomTypes.size())
        mm.chargedAtomTypes.resize(chargedAtomTypeId+1);

        // Check that this slot is not already in use.
    SimTK_APIARGCHECK2_ALWAYS(!mm.chargedAtomTypes[chargedAtomTypeId].isValid(), mm.ApiClassName, MethodName, 
        "charged atom type Id %d is already in use for '%s'", chargedAtomTypeId, 
        mm.chargedAtomTypes[chargedAtomTypeId].name.c_str());

        // Define the new charged atom type.
    mm.chargedAtomTypes[chargedAtomTypeId] = 
        ChargedAtomType(chargedAtomTypeId, typeName, atomClassId, partialChargeInE);
}

void DuMMForceFieldSubsystem::defineBondStretch
   (int class1, int class2, Real stiffnessInKJPerNmSq, Real nominalLengthInNm)
{
    static const char* MethodName = "defineBondStretch";
    DuMMForceFieldSubsystemRep& mm = updRep();

        // Watch for nonsense arguments.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtomClass(class1), mm.ApiClassName, MethodName, 
        "class1=%d which is not a valid atom class Id", class1);
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtomClass(class2), mm.ApiClassName, MethodName, 
        "class2=%d which is not a valid atom class Id", class2);
    SimTK_APIARGCHECK1_ALWAYS(stiffnessInKJPerNmSq >= 0, mm.ApiClassName, MethodName, 
        "stiffness %g is not valid: must be nonnegative", stiffnessInKJPerNmSq);
    SimTK_APIARGCHECK1_ALWAYS(nominalLengthInNm >= 0, mm.ApiClassName, MethodName, 
        "nominal length %g is not valid: must be nonnegative", nominalLengthInNm);

        // Attempt to insert the new bond stretch entry, canonicalizing first
        // so that the atom class pair has the lower class Id first.
    const IntPair key(class1,class2,true);
    std::pair<std::map<IntPair,BondStretch>::iterator, bool> ret = 
      mm.bondStretch.insert(std::pair<IntPair,BondStretch>
        (key, BondStretch(stiffnessInKJPerNmSq,nominalLengthInNm)));

        // Throw an exception if this bond stretch term was already defined. (std::map 
        // indicates that with a bool in the return value.)
    SimTK_APIARGCHECK2_ALWAYS(ret.second, mm.ApiClassName, MethodName, 
        "there was already a bond stretch term for atom class pair (%d,%d)", key[0], key[1]);
}

void DuMMForceFieldSubsystem::defineBondBend
   (int class1, int class2, int class3, Real stiffnessInKJPerRadSq, Real nominalAngleInDeg)
{
    static const char* MethodName = "defineBondBend";
    DuMMForceFieldSubsystemRep& mm = updRep();

        // Watch for nonsense arguments.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtomClass(class1), mm.ApiClassName, MethodName, 
        "class1=%d which is not a valid atom class Id", class1);
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtomClass(class2), mm.ApiClassName, MethodName, 
        "class2=%d which is not a valid atom class Id", class2);
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtomClass(class3), mm.ApiClassName, MethodName, 
        "class3=%d which is not a valid atom class Id", class3);
    SimTK_APIARGCHECK1_ALWAYS(stiffnessInKJPerRadSq >= 0, mm.ApiClassName, MethodName, 
        "stiffness %g is not valid: must be nonnegative", stiffnessInKJPerRadSq);
    SimTK_APIARGCHECK1_ALWAYS(0 <= nominalAngleInDeg && nominalAngleInDeg <= 180, 
        mm.ApiClassName, MethodName, 
        "nominal angle %g is not valid: must be between 0 and 180 degrees, inclusive", 
        nominalAngleInDeg);

        // Attempt to insert the new bond bend entry, canonicalizing first
        // by reversing the class Id triple if necessary so that the first 
        // classId is no larger than the third.
    const IntTriple key(class1, class2, class3, true);
    std::pair<std::map<IntTriple,BondBend>::iterator, bool> ret = 
      mm.bondBend.insert(std::pair<IntTriple,BondBend>
        (key, BondBend(stiffnessInKJPerRadSq,nominalAngleInDeg)));

        // Throw an exception if this bond bend term was already defined. (std::map 
        // indicates that with a bool in the return value.)
    SimTK_APIARGCHECK3_ALWAYS(ret.second, mm.ApiClassName, MethodName, 
        "there was already a bond bend term for atom class triple (%d,%d,%d)", 
        key[0], key[1], key[2]);
}

// 
// We allow up to 3 terms in a single torsion function, with three different
// periodicities. If any of these are unused, set the corresponding periodicity
// to -1.
//
void DuMMForceFieldSubsystem::defineBondTorsion
   (int class1, int class2, int class3, int class4, 
    int periodicity1, Real amp1InKJ, Real phase1InDegrees,
    int periodicity2, Real amp2InKJ, Real phase2InDegrees,
    int periodicity3, Real amp3InKJ, Real phase3InDegrees)
{
    static const char* MethodName = "defineBondTorsion";
    DuMMForceFieldSubsystemRep& mm = updRep();

        // Watch for nonsense arguments.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtomClass(class1), mm.ApiClassName, MethodName, 
        "class1=%d which is not a valid atom class Id", class1);
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtomClass(class2), mm.ApiClassName, MethodName, 
        "class2=%d which is not a valid atom class Id", class2);
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtomClass(class3), mm.ApiClassName, MethodName, 
        "class3=%d which is not a valid atom class Id", class3);
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtomClass(class4), mm.ApiClassName, MethodName, 
        "class4=%d which is not a valid atom class Id", class4);
    SimTK_APIARGCHECK_ALWAYS(periodicity1!=-1 || periodicity2!=-1 || periodicity3!=-1, 
        mm.ApiClassName, MethodName, "must be at least one torsion term supplied");

        // Canonicalize atom class quad by reversing order if necessary so that the
        // first class Id is numerically no larger than the fourth.
    const IntQuad key(class1, class2, class3, class4, true);

        // Now allocate an empty BondTorsion object and add terms to it as they are found.
    BondTorsion bt;
    if (periodicity1 != -1) {
            // No nonsense.
        SimTK_APIARGCHECK1_ALWAYS(1 <= periodicity1 && periodicity1 <= 6, mm.ApiClassName, MethodName, 
            "periodicity1(%d) is invalid: we require 1 <= periodicity <= 6", periodicity1);
        SimTK_APIARGCHECK1_ALWAYS(amp1InKJ >= 0, mm.ApiClassName, MethodName, 
            "amplitude1(%g) is not valid: must be nonnegative", amp1InKJ);
        SimTK_APIARGCHECK1_ALWAYS(0 <= phase1InDegrees && phase1InDegrees <= 180, mm.ApiClassName, MethodName, 
            "phaseAngle1(%g) is not valid: must be between 0 and 180 degrees, inclusive", phase1InDegrees);

            // No repeats.
        SimTK_APIARGCHECK1_ALWAYS((periodicity2 != periodicity1) && (periodicity3 != periodicity1), 
            mm.ApiClassName, MethodName,
            "only one term with a given periodicity may be specified (periodicity %d was repeated)",
            periodicity1);

            // Add the new term.
        bt.addTerm(TorsionTerm(periodicity1, amp1InKJ, phase1InDegrees));
    }
    if (periodicity2 != -1) {
            // No nonsense.
        SimTK_APIARGCHECK1_ALWAYS(1 <= periodicity2 && periodicity2 <= 6, mm.ApiClassName, MethodName, 
            "periodicity2(%d) is invalid: we require 1 <= periodicity <= 6", periodicity2);
        SimTK_APIARGCHECK1_ALWAYS(amp2InKJ >= 0, mm.ApiClassName, MethodName, 
            "amplitude2(%g) is not valid: must be nonnegative", amp2InKJ);
        SimTK_APIARGCHECK1_ALWAYS(0 <= phase2InDegrees && phase2InDegrees <= 180, mm.ApiClassName, MethodName, 
            "phaseAngle2(%g) is not valid: must be between 0 and 180 degrees, inclusive", phase2InDegrees);

            // No repeats.
        SimTK_APIARGCHECK1_ALWAYS(periodicity3 != periodicity2, mm.ApiClassName, MethodName,
            "only one term with a given periodicity may be specified (periodicity %d was repeated)",
            periodicity2);

            // Add the new term.
        bt.addTerm(TorsionTerm(periodicity2, amp2InKJ, phase2InDegrees));
    }
    if (periodicity3 != -1) {
            // No nonsense.
        SimTK_APIARGCHECK1_ALWAYS(1 <= periodicity3 && periodicity3 <= 6, mm.ApiClassName, MethodName, 
            "periodicity3(%d) is invalid: we require 1 <= periodicity <= 6", periodicity3);
        SimTK_APIARGCHECK1_ALWAYS(amp3InKJ >= 0, mm.ApiClassName, MethodName, 
            "amplitude3(%g) is not valid: must be nonnegative", amp3InKJ);
        SimTK_APIARGCHECK1_ALWAYS(0 <= phase3InDegrees && phase3InDegrees <= 180, mm.ApiClassName, MethodName, 
            "phaseAngle3(%g) is not valid: must be between 0 and 180 degrees, inclusive", phase3InDegrees);
            // (we've already checked for any possible repeats)

            // Add the new term.
        bt.addTerm(TorsionTerm(periodicity3, amp3InKJ, phase3InDegrees));
    }

        // Now try to insert the allegedly new BondTorsion specification into the bondTorsion map.
        // If it is already there the 2nd element in the returned pair will be 'false'.
    std::pair<std::map<IntQuad,BondTorsion>::iterator, bool> ret = 
      mm.bondTorsion.insert(std::pair<IntQuad,BondTorsion>(key,bt));

        // Throw an exception if terms for this bond torsion were already defined.
    SimTK_APIARGCHECK4_ALWAYS(ret.second, mm.ApiClassName, MethodName, 
        "bond torsion term(s) were already defined for atom class quad (%d,%d,%d,%d)", 
        key[0], key[1], key[2], key[3]);
}

// Convenient signature for a bond torsion with only one term.
void DuMMForceFieldSubsystem::defineBondTorsion
   (int class1, int class2, int class3, int class4, 
    int periodicity1, Real amp1InKJ, Real phase1InDegrees)
{
    defineBondTorsion(class1, class2, class3, class4, 
                      periodicity1,amp1InKJ,phase1InDegrees,
                      -1,0.,0., -1,0.,0.);
}

// Convenient signature for a bond torsion with two terms.
void DuMMForceFieldSubsystem::defineBondTorsion
   (int class1, int class2, int class3, int class4, 
    int periodicity1, Real amp1InKJ, Real phase1InDegrees,
    int periodicity2, Real amp2InKJ, Real phase2InDegrees)
{
    defineBondTorsion(class1, class2, class3, class4, 
                      periodicity1,amp1InKJ,phase1InDegrees,
                      periodicity2,amp2InKJ,phase2InDegrees,
                      -1,0.,0.);
}

void DuMMForceFieldSubsystem::setVdwMixingRule(VdwMixingRule rule) {
    static const char* MethodName = "setVdwMixingRule";
    DuMMForceFieldSubsystemRep& mm = updRep();
    mm.vdwMixingRule = rule; 
}

DuMMForceFieldSubsystem::VdwMixingRule 
DuMMForceFieldSubsystem::getVdwMixingRule() const {
    static const char* MethodName = "getVdwMixingRule";
    const DuMMForceFieldSubsystemRep& mm = getRep();
    return mm.vdwMixingRule; 
}

const char*
DuMMForceFieldSubsystem::getVdwMixingRuleName(VdwMixingRule rule) const {
    static const char* MethodName = "getVdwMixingRuleName";
    switch(rule) {
    case WaldmanHagler:     return "Waldman-Hagler";
    case HalgrenHHG:        return "Halgren-HHG";        
    case Jorgensen:         return "Jorgensen";        
    case LorentzBerthelot:  return "Lorentz-Berthelot"; 
    case Kong:              return "Kong";          
    default:
        SimTK_APIARGCHECK1_ALWAYS(false, "DuMMForceFieldSubsystem", MethodName,
        "Unknown van der Waals mixing rule %d", (int)rule);
    };
}

void DuMMForceFieldSubsystem::setVdw12ScaleFactor(Real fac) {
    static const char* MethodName = "setVdw12ScaleFactor";
    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac && fac <= 1, mm.ApiClassName, MethodName,
        "van der Waals energy scale factor (%g) for 1-2 bonded atoms was invalid: must be between 0 and 1, inclusive",
        fac);

    mm.vdwScale12=fac;
}
void DuMMForceFieldSubsystem::setVdw13ScaleFactor(Real fac) {
    static const char* MethodName = "setVdw13ScaleFactor";
    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac && fac <= 1, mm.ApiClassName, MethodName,
        "van der Waals energy scale factor (%g) for 1-3 bonded atoms was invalid: must be between 0 and 1, inclusive",
        fac);

    mm.vdwScale13=fac;
}
void DuMMForceFieldSubsystem::setVdw14ScaleFactor(Real fac) {
    static const char* MethodName = "setVdw14ScaleFactor";
    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac && fac <= 1, mm.ApiClassName, MethodName,
        "van der Waals energy scale factor (%g) for 1-4 bonded atoms was invalid: must be between 0 and 1, inclusive",
        fac);

    mm.vdwScale14=fac;
}
void DuMMForceFieldSubsystem::setVdw15ScaleFactor(Real fac) {
    static const char* MethodName = "setVdw15ScaleFactor";
    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac && fac <= 1, mm.ApiClassName, MethodName,
        "van der Waals energy scale factor (%g) for 1-5 bonded atoms was invalid: must be between 0 and 1, inclusive",
        fac);

    mm.vdwScale15=fac;
}

void DuMMForceFieldSubsystem::setCoulomb12ScaleFactor(Real fac) {
    static const char* MethodName = "setCoulomb12ScaleFactor";
    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac && fac <= 1, mm.ApiClassName, MethodName,
        "Coulomb scale factor (%g) for 1-2 bonded atoms was invalid: must be between 0 and 1, inclusive",
        fac);

    mm.coulombScale12=fac;
}

void DuMMForceFieldSubsystem::setCoulomb13ScaleFactor(Real fac) {
    static const char* MethodName = "setCoulomb13ScaleFactor";
    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac && fac <= 1, mm.ApiClassName, MethodName,
        "Coulomb scale factor (%g) for 1-3 bonded atoms was invalid: must be between 0 and 1, inclusive",
        fac);

    mm.coulombScale13=fac;
}
void DuMMForceFieldSubsystem::setCoulomb14ScaleFactor(Real fac) {
    static const char* MethodName = "setCoulomb14ScaleFactor";
    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac && fac <= 1, mm.ApiClassName, MethodName,
        "Coulomb scale factor (%g) for 1-4 bonded atoms was invalid: must be between 0 and 1, inclusive",
        fac);

    mm.coulombScale14=fac;
}
void DuMMForceFieldSubsystem::setCoulomb15ScaleFactor(Real fac) {
    static const char* MethodName = "setCoulomb15ScaleFactor";
    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac && fac <= 1, mm.ApiClassName, MethodName,
        "Coulomb scale factor (%g) for 1-5 bonded atoms was invalid: must be between 0 and 1, inclusive",
        fac);

    mm.coulombScale15=fac;
}

void DuMMForceFieldSubsystem::setVdwGlobalScaleFactor(Real fac) {
    static const char* MethodName = "setVdwScaleFactor";
    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac, mm.ApiClassName, MethodName,
        "Global van der Waals scale factor (%g) was invalid: must be nonnegative",
        fac);

    mm.vdwGlobalScaleFactor=fac;
}

void DuMMForceFieldSubsystem::setCoulombGlobalScaleFactor(Real fac) {
    static const char* MethodName = "setCoulombScaleFactor";
    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac, mm.ApiClassName, MethodName,
        "Global Coulomb scale factor (%g) was invalid: must be nonnegative",
        fac);

    mm.coulombGlobalScaleFactor=fac;
}
void DuMMForceFieldSubsystem::setBondStretchGlobalScaleFactor(Real fac) {
    static const char* MethodName = "setBondStretchScaleFactor";
    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac, mm.ApiClassName, MethodName,
        "Global bond stretch scale factor (%g) was invalid: must be nonnegative",
        fac);

    mm.bondStretchGlobalScaleFactor=fac;
}
void DuMMForceFieldSubsystem::setBondBendGlobalScaleFactor(Real fac) {
    static const char* MethodName = "setBondBendScaleFactor";
    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac, mm.ApiClassName, MethodName,
        "Global bond bend scale factor (%g) was invalid: must be nonnegative",
        fac);

    mm.bondBendGlobalScaleFactor=fac;
}
void DuMMForceFieldSubsystem::setBondTorsionGlobalScaleFactor(Real fac) {
    static const char* MethodName = "setBondTorsionScaleFactor";
    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac, mm.ApiClassName, MethodName,
        "Global bond torsion scale factor (%g) was invalid: must be nonnegative",
        fac);

    mm.bondTorsionGlobalScaleFactor=fac;
}

int DuMMForceFieldSubsystem::createCluster(const char* groupName)
{
    // Currently there is no error checking to do. We don't insist on unique group names.
    return updRep().addCluster(Cluster(groupName));
}

int DuMMForceFieldSubsystem::addAtom(int chargedAtomTypeId)
{
    static const char* MethodName = "addAtom";
    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(mm.isValidChargedAtomType(chargedAtomTypeId), mm.ApiClassName, MethodName, 
        "charged atom type %d is not valid", chargedAtomTypeId);

    const int atomId = (int)mm.atoms.size();
    mm.atoms.push_back(Atom(chargedAtomTypeId, atomId));
    return atomId;
}

void DuMMForceFieldSubsystem::placeAtomInCluster(int atomId, int clusterId, const Vec3& stationInNm)
{
    static const char* MethodName = "placeAtomInCluster";
    DuMMForceFieldSubsystemRep& mm = updRep();

        // Make sure that we've seen both the atomId and clusterId before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtom(atomId), mm.ApiClassName, MethodName,
        "atom Id %d is not valid", atomId);
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidCluster(clusterId), mm.ApiClassName, MethodName,
        "cluster Id %d is not valid", clusterId);

    Cluster& cluster = mm.updCluster(clusterId);

        // Make sure that this cluster doesn't already contain this atom, either directly
        // or recursively through its subclusters.
    SimTK_APIARGCHECK3_ALWAYS(!cluster.containsAtom(atomId), mm.ApiClassName, MethodName,
        "cluster %d('%s') already contains atom %d", clusterId, cluster.name.c_str(), atomId);

        // Add the atom to the cluster.
    cluster.placeAtom(atomId, stationInNm, mm);
}

void DuMMForceFieldSubsystem::placeClusterInCluster
   (int childClusterId, int parentClusterId, const Transform& placementInNm)
{
    static const char* MethodName = "placeClusterInCluster";
    DuMMForceFieldSubsystemRep& mm = updRep();

        // Make sure that we've seen both of these clusters before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidCluster(childClusterId), mm.ApiClassName, MethodName,
        "child cluster Id %d is not valid", childClusterId);
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidCluster(parentClusterId), mm.ApiClassName, MethodName,
        "parent cluster Id %d is not valid", parentClusterId);

    Cluster&       parent = mm.updCluster(parentClusterId);
    const Cluster& child  = mm.getCluster(childClusterId);

        // TODO: for now, make sure the parent is a top-level cluster, meaning that it does
        // not have any parent clusters (although it can be attached to a body). This restriction
        // should be relaxed but it is tricky to get all the parents' and ancestors' content
        // lists updated correctly so I'm deferring that for now (sherm 060928).
    SimTK_APIARGCHECK2_ALWAYS(parent.isTopLevelCluster(), mm.ApiClassName, MethodName,
        "parent cluster %d('%s') is not a top-level cluster so you cannot add a child cluster to it now",
        parentClusterId, parent.name.c_str());

        // Child must not already be attached to a body.
    SimTK_APIARGCHECK2_ALWAYS(!child.isAttachedToBody(), mm.ApiClassName, MethodName,
        "child cluster %d('%s') is already attached to a body so cannot now be placed in another cluster",
        childClusterId, child.name.c_str());

        // Make sure that parent cluster doesn't already contain child cluster, either directly
        // or recursively through its subclusters.
    SimTK_APIARGCHECK4_ALWAYS(!parent.containsCluster(childClusterId), mm.ApiClassName, MethodName,
        "parent cluster %d('%s') already contains child cluster %d('%s')", 
        parentClusterId, parent.name.c_str(), childClusterId, child.name.c_str());

        // Make sure the new child cluster doesn't contain any atoms which are already in
        // any of the trees to which the parent cluster is associated.
        // TODO: for now we need only look at the parent since we know it is top level.
    int atomId;
    SimTK_APIARGCHECK5_ALWAYS(!parent.overlapsWithCluster(child, atomId), mm.ApiClassName, MethodName,
        "parent cluster %d('%s') and would-be child cluster %d('%s') both contain atom %d"
        " so they cannot have a parent/child relationship",
        parentClusterId, parent.name.c_str(), childClusterId, child.name.c_str(), atomId);

        // Add the child cluster to the parent.
    parent.placeCluster(childClusterId, placementInNm, mm);
}

void DuMMForceFieldSubsystem::attachClusterToBody(int clusterId, int bodyNum, 
                                                  const Transform& placementInNm) 
{
    static const char* MethodName = "attachClusterToBody";
    DuMMForceFieldSubsystemRep& mm = updRep();

        // Make sure we've seen this cluster before, and that the body number is well formed.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidCluster(clusterId), mm.ApiClassName, MethodName,
        "cluster Id %d is not valid", clusterId);
    SimTK_APIARGCHECK1_ALWAYS(bodyNum >= 0, mm.ApiClassName, MethodName,
        "body number %d is not valid: must be nonnegative", bodyNum);

    const Cluster& child  = mm.getCluster(clusterId);

        // Child must not already be attached to a body.
    SimTK_APIARGCHECK3_ALWAYS(!child.isAttachedToBody(), mm.ApiClassName, MethodName,
        "cluster %d('%s') is already attached to body %d so cannot now be attached to a body",
        clusterId, child.name.c_str(), child.getBodyId());

        // None of the atoms in the child can be attached to any body.
    int atomId, bodyId;
    SimTK_APIARGCHECK4_ALWAYS(!child.containsAnyAtomsAttachedToABody(atomId,bodyId,mm), 
        mm.ApiClassName, MethodName,
        "cluster %d('%s') contains atom %d which is already attached to body %d"
        " so the cluster cannot now be attached to another body",
        clusterId, child.name.c_str(), atomId, bodyId);

        // Create an entry for the body if necessary, and its corresponding cluster.
    mm.ensureBodyEntryExists(bodyNum);
    Cluster& bodyCluster = mm.updCluster(mm.getBody(bodyNum).getClusterId());

        // Make sure that body cluster doesn't already contain child cluster, either directly
        // or recursively through its subclusters.
    SimTK_APIARGCHECK3_ALWAYS(!bodyCluster.containsCluster(clusterId), mm.ApiClassName, MethodName,
        "cluster %d('%s') is already attached (directly or indirectly) to body %d", 
        clusterId, child.name.c_str(), bodyNum);

        // OK, attach the cluster to the body's cluster.
    bodyCluster.placeCluster(clusterId, placementInNm, mm);
}

void DuMMForceFieldSubsystem::attachAtomToBody(int atomId, int bodyNum, const Vec3& stationInNm) 
{
    static const char* MethodName = "attachAtomToBody";
    DuMMForceFieldSubsystemRep& mm = updRep();

        // Make sure we've seen this atom before, and that the body number is well formed.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtom(atomId), mm.ApiClassName, MethodName,
        "atom Id %d is not valid", atomId);
    SimTK_APIARGCHECK1_ALWAYS(bodyNum >= 0, mm.ApiClassName, MethodName,
        "body number %d is not valid: must be nonnegative", bodyNum);

        // The atom must not already be attached to a body, even this one.
    SimTK_APIARGCHECK2_ALWAYS(!mm.getAtom(atomId).isAttachedToBody(), mm.ApiClassName, MethodName,
        "atom %d is already attached to body %d so cannot now be attached to a body",
        atomId, mm.getAtom(atomId).getBodyId());

        // Create an entry for the body if necessary, and its corresponding cluster.
    mm.ensureBodyEntryExists(bodyNum);
    Cluster& bodyCluster = mm.updCluster(mm.getBody(bodyNum).getClusterId());

        // Attach the atom to the body's cluster.
    bodyCluster.placeAtom(atomId, stationInNm, mm);
}

MassProperties DuMMForceFieldSubsystem::calcClusterMassProperties
   (int clusterId, const Transform& placementInNm) const
{
    static const char* MethodName = "calcClusterMassProperties";
    const DuMMForceFieldSubsystemRep& mm = getRep();

        // Make sure we've seen this cluster before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidCluster(clusterId), mm.ApiClassName, MethodName,
        "cluster Id %d is not valid", clusterId);

    return mm.getCluster(clusterId).calcMassProperties(placementInNm, mm);
}


int DuMMForceFieldSubsystem::addBond(int atom1Id, int atom2Id)
{
    static const char* MethodName = "addBond";
    DuMMForceFieldSubsystemRep& mm = updRep();

        // Make sure we've seen these atoms before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtom(atom1Id), mm.ApiClassName, MethodName,
        "atom1(%d) is not valid", atom1Id);
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtom(atom2Id), mm.ApiClassName, MethodName,
        "atom2(%d) is not valid", atom2Id);

        // An atom can't be bonded to itself.
    SimTK_APIARGCHECK1_ALWAYS(atom1Id != atom2Id, mm.ApiClassName, MethodName,
        "the same atom Id(%d) was given for both atoms, which makes no sense", atom1Id);

    // Ensure that atom1 < atom2
    if (atom1Id > atom2Id)
        std::swap(atom1Id,atom2Id);

    Atom& a1 = mm.updAtom(atom1Id);
    Atom& a2 = mm.updAtom(atom2Id);

    SimTK_APIARGCHECK2_ALWAYS(!a1.isBondedTo(atom2Id), mm.ApiClassName, MethodName,
        "atom %d is already bonded to atom %d; you can only do that once",
        atom1Id, atom2Id);

    mm.bonds.push_back(Bond(atom1Id,atom2Id));
    a1.bond12.push_back(atom2Id);
    a2.bond12.push_back(atom1Id);
    return (int)mm.bonds.size() - 1;
}

int DuMMForceFieldSubsystem::getNAtoms() const {
    return getRep().getNAtoms();
}
int DuMMForceFieldSubsystem::getNBonds() const {
    return getRep().getNBonds();
}

// 'which' is 0 or 1 to pick one of the two atoms whose Id we return.
int DuMMForceFieldSubsystem::getBondAtom(int bondId, int which) const {
    static const char* MethodName = "getBondAtom";
    const DuMMForceFieldSubsystemRep& mm = getRep();

        // Make sure we've seen this bond before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidBond(bondId), mm.ApiClassName, MethodName,
        "bond %d is not valid", bondId);

    SimTK_APIARGCHECK1_ALWAYS(which==0 || which==1, mm.ApiClassName, MethodName,
        "'which' was %d but must be 0 or 1 to choose one of the two atoms", which);

    return mm.bonds[bondId].atoms[which];
}

// Returned mass is in daltons (g/mol).
Real DuMMForceFieldSubsystem::getAtomMass(int atomId) const {
    static const char* MethodName = "getAtomMass";
    const DuMMForceFieldSubsystemRep& mm = getRep();

        // Make sure we've seen this atom before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtom(atomId), mm.ApiClassName, MethodName,
        "atom %d is not valid", atomId);

    const Element& e = mm.elements[mm.getAtomElementNum(atomId)];
    return e.mass;
}

// Returns the atomic number (number of protons in nucleus).
int DuMMForceFieldSubsystem::getAtomElement(int atomId) const {
    static const char* MethodName = "getAtomElement";
    const DuMMForceFieldSubsystemRep& mm = getRep();

        // Make sure we've seen this atom before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtom(atomId), mm.ApiClassName, MethodName,
        "atom %d is not valid", atomId);

    return mm.getAtomElementNum(atomId);
}

Vec3 DuMMForceFieldSubsystem::getAtomDefaultColor(int atomId) const {
    static const char* MethodName = "getAtomDefaultColor";
    const DuMMForceFieldSubsystemRep& mm = getRep();

        // Make sure we've seen this atom before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtom(atomId), mm.ApiClassName, MethodName,
        "atom %d is not valid", atomId);

    const Element& e = mm.elements[mm.getAtomElementNum(atomId)];
    return e.defaultColor;
}

// Returned radius is in nm.
Real DuMMForceFieldSubsystem::getAtomRadius(int atomId) const {
    static const char* MethodName = "getAtomRadius";
    const DuMMForceFieldSubsystemRep& mm = getRep();

        // Make sure we've seen this atom before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtom(atomId), mm.ApiClassName, MethodName,
        "atom %d is not valid", atomId);

    const AtomClass& cl = mm.atomClasses[mm.getAtomClassId(atomId)];
    return cl.vdwRadius;
}

// Returned station is in nm.
Vec3 DuMMForceFieldSubsystem::getAtomStationOnBody(int atomId) const {
    static const char* MethodName = "getAtomStationOnBody";
    const DuMMForceFieldSubsystemRep& mm = getRep();

        // Make sure we've seen this atom before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtom(atomId), mm.ApiClassName, MethodName,
        "atom %d is not valid", atomId);

    const Atom& a = mm.getAtom(atomId);

        // Atom must be attached to a body.
    SimTK_APIARGCHECK1_ALWAYS(a.isAttachedToBody(), mm.ApiClassName, MethodName,
        "atom %d is not attached to a body", atomId);

    return a.station_B;
}

// Returned placement is in nm.
Transform DuMMForceFieldSubsystem::getClusterPlacementOnBody(int clusterId) const {
    static const char* MethodName = "getClusterPlacementOnBody";
    const DuMMForceFieldSubsystemRep& mm = getRep();

        // Make sure we've seen this cluster before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidCluster(clusterId), mm.ApiClassName, MethodName,
        "cluster Id %d is not valid", clusterId);

    const Cluster& c = mm.getCluster(clusterId);

        // Cluster must be attached to a body.
    SimTK_APIARGCHECK2_ALWAYS(c.isAttachedToBody(), mm.ApiClassName, MethodName,
        "cluster %d('%s') is not attached to a body", clusterId, c.name.c_str());

    return c.placement_B;
}

// Returned station is in nm.
Vec3 DuMMForceFieldSubsystem::getAtomStationInCluster(int atomId, int clusterId) const {
    static const char* MethodName = "getAtomStationInCluster";
    const DuMMForceFieldSubsystemRep& mm = getRep();

        // Make sure that we've seen both the atomId and clusterId before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtom(atomId), mm.ApiClassName, MethodName,
        "atom Id %d is not valid", atomId);
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidCluster(clusterId), mm.ApiClassName, MethodName,
        "cluster Id %d is not valid", clusterId);

    const Cluster& c = mm.getCluster(clusterId);
    const AtomPlacementSet& atoms = c.getAllContainedAtoms();
    const AtomPlacementSet::const_iterator ap = 
        atoms.find(AtomPlacement(atomId,Vec3(0)));

        // We're going to be upset of this cluster doesn't contain this atom.
    SimTK_APIARGCHECK3_ALWAYS(ap != atoms.end(), mm.ApiClassName, MethodName,
        "cluster %d('%s') does not contain atom %d", clusterId, c.name.c_str(), atomId);

    return ap->station;
}

// Returned placement is in nm.
Transform DuMMForceFieldSubsystem::getClusterPlacementInCluster(int childClusterId, int parentClusterId) const {
    static const char* MethodName = "getClusterPlacementInCluster";
    const DuMMForceFieldSubsystemRep& mm = getRep();

        // Make sure that we've seen both of these clusters before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidCluster(childClusterId), mm.ApiClassName, MethodName,
        "child cluster Id %d is not valid", childClusterId);
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidCluster(parentClusterId), mm.ApiClassName, MethodName,
        "parent cluster Id %d is not valid", parentClusterId);

    const Cluster& parent = mm.getCluster(parentClusterId);
    const Cluster& child  = mm.getCluster(childClusterId);

    const ClusterPlacementSet& clusters = parent.getAllContainedClusters();
    const ClusterPlacementSet::const_iterator cp = 
        clusters.find(ClusterPlacement(childClusterId,Transform()));

        // We're going to be upset of the parent cluster doesn't contain the child.
    SimTK_APIARGCHECK4_ALWAYS(cp != clusters.end(), mm.ApiClassName, MethodName,
        "cluster %d('%s') does not contain cluster %d('%d')", 
        parentClusterId, parent.name.c_str(), childClusterId, child.name.c_str());

    return cp->placement;
}

int DuMMForceFieldSubsystem::getAtomBody(int atomId) const {
    static const char* MethodName = "getAtomBody";
    const DuMMForceFieldSubsystemRep& mm = getRep();

        // Make sure that we've seen this atomId before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtom(atomId), mm.ApiClassName, MethodName,
        "atom Id %d is not valid", atomId);

    const Atom& a = mm.getAtom(atomId);

        // Atom must be attached to a body.
    SimTK_APIARGCHECK1_ALWAYS(a.isAttachedToBody(), mm.ApiClassName, MethodName,
        "atom %d is not attached to a body", atomId);

    return a.getBodyId();
}


int DuMMForceFieldSubsystem::getClusterBody(int clusterId) const {
    static const char* MethodName = "getClusterBody";
    const DuMMForceFieldSubsystemRep& mm = getRep();

        // Make sure that we've seen this atomId before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidCluster(clusterId), mm.ApiClassName, MethodName,
        "cluster Id %d is not valid", clusterId);

    const Cluster& c = mm.getCluster(clusterId);

        // Cluster must be attached to a body.
    SimTK_APIARGCHECK2_ALWAYS(c.isAttachedToBody(), mm.ApiClassName, MethodName,
        "cluster %d('%s') is not attached to a body", clusterId, c.name.c_str());

    return c.getBodyId();
}

void DuMMForceFieldSubsystem::dump() const {
    return getRep().dump();
}



    ////////////////////////////////
    // DuMMForceFieldSubsystemRep //
    ////////////////////////////////

/*static*/ const char* DuMMForceFieldSubsystemRep::ApiClassName 
    = "DuMMForceFieldSubsystem";

const BondStretch& 
DuMMForceFieldSubsystemRep::getBondStretch(int class1, int class2) const {
    static const BondStretch dummy; // invalid
    const IntPair key(class1,class2,true);
    std::map<IntPair,BondStretch>::const_iterator bs = bondStretch.find(key);
    return (bs != bondStretch.end()) ? bs->second : dummy;
}

const BondBend& 
DuMMForceFieldSubsystemRep::getBondBend(int class1, int class2, int class3) const {
    static const BondBend dummy; // invalid
    const IntTriple key(class1, class2, class3, true);
    std::map<IntTriple,BondBend>::const_iterator bb = bondBend.find(key);
    return (bb != bondBend.end()) ? bb->second : dummy;
}

const BondTorsion& 
DuMMForceFieldSubsystemRep::getBondTorsion
   (int class1, int class2, int class3, int class4) const
{
    static const BondTorsion dummy; // invalid
    const IntQuad key(class1, class2, class3, class4, true);
    std::map<IntQuad,BondTorsion>::const_iterator bt = bondTorsion.find(key);
    return (bt != bondTorsion.end()) ? bt->second : dummy;
}

void DuMMForceFieldSubsystemRep::realizeTopology(State& s) const {
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

    // Process clusters & bodies (bodies are treated as top-level clusters)

    // We process clusters recursively, so we need to allow the clusters writable
    // access to the main DuMM object (i.e., "this").
    for (int cnum=0; cnum < (int)clusters.size(); ++cnum) {
        Cluster& c = mutableThis->clusters[cnum];
        assert(c.isValid()); // Shouldn't be any unused cluster numbers.
        c.realizeTopologicalCache(*mutableThis);
    }

    // Bodies, on the other hand, are always top level clusters and the
    // calculation here assumes that all the clusters have been processed.
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
            const AtomPlacement& ap = b.allAtoms[i]; assert(ap.isValid());
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
        allBondedSoFar.insert(a.bond12.begin(), a.bond12.end());

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

        const int c1 = getAtomClassId(anum);

        // Save a BondStretch entry for each cross-body 1-2 bond
        a.stretch.resize(a.xbond12.size());
        for (int b12=0; b12 < (int)a.xbond12.size(); ++b12) {
            const int c2 = getAtomClassId(a.xbond12[b12]);
            a.stretch[b12] = getBondStretch(c1, c2);

            SimTK_REALIZECHECK2_ALWAYS(a.stretch[b12].isValid(),
                Stage::Topology, getMySubsystemIndex(), getName(),
                "couldn't find bond stretch parameters for cross-body atom class pair (%d,%d)", 
                c1,c2);
        }

        // Save a BondBend entry for each cross-body 1-3 bond
        a.bend.resize(a.xbond13.size());
        for (int b13=0; b13 < (int)a.xbond13.size(); ++b13) {
            const int c2 = getAtomClassId(a.xbond13[b13][0]);
            const int c3 = getAtomClassId(a.xbond13[b13][1]);
            a.bend[b13] = getBondBend(c1, c2, c3);

            SimTK_REALIZECHECK3_ALWAYS(a.bend[b13].isValid(),
                Stage::Topology, getMySubsystemIndex(), getName(),
                "couldn't find bond bend parameters for cross-body atom class triple (%d,%d,%d)", 
                c1,c2,c3);
        }

        // Save a BondTorsion entry for each cross-body 1-4 bond
        a.torsion.resize(a.xbond14.size());
        for (int b14=0; b14 < (int)a.xbond14.size(); ++b14) {
            const int c2 = getAtomClassId(a.xbond14[b14][0]);
            const int c3 = getAtomClassId(a.xbond14[b14][1]);
            const int c4 = getAtomClassId(a.xbond14[b14][2]);
            a.torsion[b14] = getBondTorsion(c1, c2, c3, c4); 

            SimTK_REALIZECHECK4_ALWAYS(a.torsion[b14].isValid(),
                Stage::Topology, getMySubsystemIndex(), getName(),
                "couldn't find bond torsion parameters for cross-body atom class quad (%d,%d,%d,%d)", 
                c1,c2,c3,c4);
        }
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
    Real&                  pe              = mbs.updPotentialEnergy(s); // kJ
    Vector_<SpatialVec>&   rigidBodyForces = mbs.updRigidBodyForces(s); // kJ (torque), kJ/nm (force)

    for (int b1=0; b1 < (int)bodies.size(); ++b1) {
        const Transform&          X_GB1  = matter.getBodyPosition(s,b1);
        const AtomPlacementArray& alist1 = bodies[b1].allAtoms;

        for (int i=0; i < (int)alist1.size(); ++i) {
            const int       a1num = alist1[i].atomId;
            const Atom&     a1 = atoms[a1num];
            const ChargedAtomType& a1type  = chargedAtomTypes[a1.chargedAtomTypeId];
            int                    a1cnum  = a1type.atomClassId;
            const AtomClass&       a1class = atomClasses[a1cnum];
            const Vec3      a1Station_G = X_GB1.R()*a1.station_B;
            const Vec3      a1Pos_G     = X_GB1.T() + a1Station_G;
            const Real      q1Fac = coulombGlobalScaleFactor*CoulombFac*a1type.partialCharge;

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
                const Transform& X_GB2   = matter.getBodyPosition(s, a2.bodyId);
                const Vec3       a2Station_G = X_GB2.R()*a2.station_B;
                const Vec3       a2Pos_G     = X_GB2.T() + a2Station_G;
                const Vec3       r = a2Pos_G - a1Pos_G;
                const Real       d = r.norm();

                // TODO: come up with something for when d is 0; for relaxation
                // just needs to push away from zero; what to do for dynamics?

                const BondStretch& bs = a1.stretch[b12];
                const Real         x  = d - bs.d0;
                const Real         k  = bondStretchGlobalScaleFactor*bs.k;

                const Real eStretch =  k*x*x; // no factor of 1/2!
                const Real fStretch = -2*k*x; // sign is as would be applied to a2
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
                const Transform& X_GB2   = matter.getBodyPosition(s, a2.bodyId);
                const Transform& X_GB3   = matter.getBodyPosition(s, a3.bodyId);
                const Vec3       a2Station_G = X_GB2.R()*a2.station_B;
                const Vec3       a3Station_G = X_GB3.R()*a3.station_B;
                const Vec3       a2Pos_G     = X_GB2.T() + a2Station_G;
                const Vec3       a3Pos_G     = X_GB3.T() + a3Station_G;

                Real angle, energy;
                Vec3 f1, f2, f3;
                const BondBend& bb = a1.bend[b13];
                // atom 2 is the central one
                bb.harmonic(a2Pos_G, a1Pos_G, a3Pos_G, bondBendGlobalScaleFactor, angle, energy, f2, f1, f3);

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
                const Transform& X_GB2   = matter.getBodyPosition(s, a2.bodyId);
                const Transform& X_GB3   = matter.getBodyPosition(s, a3.bodyId);
                const Transform& X_GB4   = matter.getBodyPosition(s, a4.bodyId);
                const Vec3       a2Station_G = X_GB2.R()*a2.station_B;
                const Vec3       a3Station_G = X_GB3.R()*a3.station_B;
                const Vec3       a4Station_G = X_GB4.R()*a4.station_B;
                const Vec3       a2Pos_G     = X_GB2.T() + a2Station_G;
                const Vec3       a3Pos_G     = X_GB3.T() + a3Station_G;
                const Vec3       a4Pos_G     = X_GB4.T() + a4Station_G;

                Real angle, energy;
                Vec3 f1, f2, f3, f4;
                const BondTorsion& bt = a1.torsion[b14];
                bt.periodic(a1Pos_G, a2Pos_G, a3Pos_G, a4Pos_G, bondTorsionGlobalScaleFactor,
                            angle, energy, f1, f2, f3, f4);

                pe += energy;
                rigidBodyForces[b1] += SpatialVec( a1Station_G % f1, f1);   // 15 flops
                rigidBodyForces[b2] += SpatialVec( a2Station_G % f2, f2);   // 15 flops
                rigidBodyForces[b3] += SpatialVec( a3Station_G % f3, f3);   // 15 flops
                rigidBodyForces[b4] += SpatialVec( a4Station_G % f4, f4);   // 15 flops
            }

            scaleBondedAtoms(a1,vdwScale,coulombScale);
            for (int b2=b1+1; b2 < (int)bodies.size(); ++b2) {
                const Transform&          X_GB2  = matter.getBodyPosition(s,b2);
                const AtomPlacementArray& alist2 = bodies[b2].allAtoms;

                for (int j=0; j < (int)alist2.size(); ++j) {
                    const int       a2num = alist2[j].atomId;
                    assert(a2num != a1num);
                    const Atom&     a2 = atoms[a2num];
                    const ChargedAtomType& a2type  = chargedAtomTypes[a2.chargedAtomTypeId];
                    int                    a2cnum  = a2type.atomClassId;
                    const AtomClass&       a2class = atomClasses[a2cnum];
                    
                    const Vec3  a2Station_G = X_GB2.R()*a2.station_B; // 15 flops
                    const Vec3  a2Pos_G     = X_GB2.T() + a2Station_G;  // 3 flops
                    const Vec3  r = a2Pos_G - a1Pos_G; // from a1 to a2 (3 flops)
                    const Real  d2 = r.normSqr(); // 5 flops

                    // Check for cutoffs on d2?

                    const Real  ood = 1/std::sqrt(d2); // approx 40 flops
                    const Real  ood2 = ood*ood;        // 1 flop

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

                    const Real eijScale = vdwGlobalScaleFactor*vdwScale[a2num]*eij; // 2 flops
                    const Real eVdw =      eijScale * (ddij12 - 2*ddij6); // 3 flops
                    const Real fVdw = 12 * eijScale * (ddij12 - ddij6);   // factor of 1/d^2 missing (3 flops)
                    const Vec3 fj = ((fCoulomb+fVdw)*ood2) * r;      // to apply to atom j on b2 (5 flops)

                    pe += (eCoulomb + eVdw); // kJ (Da-nm^2/ps^2) (2 flops)
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

// Element masses are given in daltons (==g/mol==amu==u).
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
    printf("  NBodies=%d NClusters=%d NAtoms=%d NAtomClasses=%d NChargedAtomTypes=%d NBonds=%d\n",
        bodies.size(), clusters.size(), atoms.size(), 
        atomClasses.size(), chargedAtomTypes.size(), bonds.size());
    for (int i=0; i < (int)bodies.size(); ++i) {
        printf("  Body %d:\n", i);
        bodies[i].dump();
    }
    for (int i=0; i < (int)clusters.size(); ++i) {
        printf("  Cluster %d:\n", i);
        clusters[i].dump();
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
   (const Vec3& cG, const Vec3& rG, const Vec3& sG, const Real& scale,
    Real& theta, Real& pe, Vec3& cf, Vec3& rf, Vec3& sf) const
{
    const Real ks = scale*k; //              1 flop
    const Vec3 r = rG - cG; //               3 flops
    const Vec3 s = sG - cG; //               3 flops
    const Real rr = ~r*r, ss = ~s*s;    // |r|^2, |s|^2 ( 10 flops)

    const Real rs = ~r * s; // r dot s      (5 flops)
    const Vec3 rxs = r % s; // r cross s    (9 flops)
    const Real rxslen = rxs.norm(); //      (~35 flops)
    theta = std::atan2(rxslen, rs); //       ~50 flops
    const Real bend = theta - theta0;   //   1 flop
    pe = ks*bend*bend; // NOTE: no factor of 1/2 (2 flops)

    // p is unit vector perpendicular to r and s

    // TODO: come up with something for when rxslen is 0 (vectors r & s
    // aligned or opposite); for relaxation
    // just needs to push them apart; what to do for dynamics?
    // Here we'll just make up a direction perpendicular to both
    // vectors and use it.
    const UnitVec3 p = (rxslen != 0 ? UnitVec3(rxs/rxslen,true)  // ~11 flops
                                    : UnitVec3(r).perp()); 
    const Real ffac = -2*ks*bend; // 2 flops
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
              const Real& scale, Real& theta, Real& pe, 
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
        pe     += terms[i].energy(theta);
        torque += terms[i].torque(theta);
    }
    pe     *= scale;
    torque *= scale;

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

    /////////////
    // Cluster //
    /////////////


void Cluster::attachToBody(int bnum, const Transform& X_BR, DuMMForceFieldSubsystemRep& mm) {
    assert(!isAttachedToBody());
    bodyId = bnum;
    placement_B = X_BR;

    // Tell all the atoms directly contained in this cluster that they are
    // now attached to the body also. This will fail if any of the atoms are
    // alread attached -- no polygamy.
    AtomPlacementSet::const_iterator ap = directAtomPlacements.begin();
    while (ap != directAtomPlacements.end()) {
        Atom& a = mm.updAtom(ap->atomId);
        a.attachToBody(bnum, X_BR*ap->station);
        ++ap;
    }

    // Now do the same for our contained groups, who will in turn notify their
    // own atoms and subgroups.
    ClusterPlacementSet::const_iterator cp = directClusterPlacements.begin();
    while (cp != directClusterPlacements.end()) {
        Cluster& c = mm.updCluster(cp->clusterId);
        c.attachToBody(bnum, X_BR*cp->placement, mm);
        ++cp;
    }
}

// Return true if this cluster contains (directly or indirectly) any atom which has already
// been attached to a body. If so return one of the attached atoms and its body, which can
// be helpful in error messages.
bool Cluster::containsAnyAtomsAttachedToABody(int& atomId, int& bodyId, 
                                              const DuMMForceFieldSubsystemRep& mm) const 
{
    const AtomPlacementSet& myAtoms   = getAllContainedAtoms();
    AtomPlacementSet::const_iterator ap = myAtoms.begin();
    while (ap != myAtoms.end()) {
        const Atom& a = mm.getAtom(ap->atomId);
        if (a.isAttachedToBody()) {
            atomId = ap->atomId;
            bodyId = a.getBodyId();
            return true;
        }
        ++ap;
    }
    atomId = bodyId = -1;
    return false;
}

// Place an atom in this cluster. To be valid, the atom must not
// already be
//   (a) in any of the trees of which this group is apart, or
//   (b) attached to a body.
// TODO: (c) at the moment we don't allow placing an atom in a group unless
//           that group is a top-level group (i.e., it has no parents).
// If this group is already attached to a body, then we will update
// the atom entry to note that it is now attached to the body also.
void Cluster::placeAtom(int atomId, const Vec3& station, DuMMForceFieldSubsystemRep& mm) {
    assert(isTopLevelCluster()); // TODO
    assert(!mm.getAtom(atomId).isAttachedToBody());
    assert(!containsAtom(atomId));

    std::pair<AtomPlacementSet::iterator, bool> ret;
    ret = directAtomPlacements.insert(AtomPlacement(atomId,station));
    assert(ret.second); // must not have been there already

    ret = allAtomPlacements.insert(AtomPlacement(atomId,station));
    assert(ret.second); // must not have been there already

    if (isAttachedToBody())
        mm.updAtom(atomId).attachToBody(bodyId, placement_B*station);
}

// Place a child cluster in this parent cluster. To be valid, the child 
// must not 
//   (a) already be contained in the parent group or one of the parent's subgroups, or
//   (b) contain any atoms which are already present in the parent or any
//       of the parent's subgroups, or
//   (c) already be attached to a body.
// TODO: (d) at the moment we don't allow adding a child group unless
//           the parent (this) group is a top-level group (i.e., it has no parents).
// If the parent is already attached to a body, then we will update
// the child to note that it is now attached to the body also (and it
// will update its contained atoms).
void Cluster::placeCluster(int childClusterId, const Transform& placement, DuMMForceFieldSubsystemRep& mm) {
    assert(isTopLevelCluster()); // TODO

    Cluster& child = mm.updCluster(childClusterId);
    assert(!child.isAttachedToBody());
    assert(!containsCluster(childClusterId));

    // Make sure the new child cluster doesn't contain any atoms which are already in
    // any of the trees to which the parent cluster (this) is associated.
    // TODO: for now we need only look at the parent since we know it is top level.
    const AtomPlacementSet& childsAtoms  = child.getAllContainedAtoms();
    AtomPlacementSet&       parentsAtoms = updAllContainedAtoms();

    // Make sure none of the child's atoms are already in the parent.
    AtomPlacementSet::const_iterator ap = childsAtoms.begin();
    while (ap != childsAtoms.end()) {
        std::pair<AtomPlacementSet::iterator, bool> ret =
            parentsAtoms.insert(AtomPlacement(ap->atomId, placement*ap->station));
        assert(ret.second); // mustn't have been there already
        ++ap;
    }

    const ClusterPlacementSet& childsClusters  = child.getAllContainedClusters();
    ClusterPlacementSet&       parentsClusters = updAllContainedClusters();

    // Make sure none of the child's atoms are already in the parent.
    ClusterPlacementSet::const_iterator cp = childsClusters.begin();
    while (cp != childsClusters.end()) {
        std::pair<ClusterPlacementSet::iterator, bool> ret =
            parentsClusters.insert(ClusterPlacement(cp->clusterId, placement*cp->placement));
        assert(ret.second); // mustn't have been there already
        ++cp;
    }

    noteNewChildCluster(childClusterId, placement);
    child.noteNewParentCluster(clusterId, placement);

    if (isAttachedToBody())
        child.attachToBody(bodyId, placement_B*placement, mm);

    //TODO: check for loops
}



// Calculate the composite mass properties for this cluster, transformed into
// the indicated frame.
MassProperties Cluster::calcMassProperties
   (const Transform& tr, const DuMMForceFieldSubsystemRep& mm) const 
{
    Real       mass = 0;
    Vec3       com(0);
    Inertia inertia(0);

    // Calculate the mass properties in the local frame and transform last.
    AtomPlacementSet::const_iterator aap = allAtomPlacements.begin();
    while (aap != allAtomPlacements.end()) {
        const Real ma = mm.getElement(mm.getAtomElementNum(aap->atomId)).mass;
        mass += ma;
        com  += ma*aap->station;
        inertia += Inertia(aap->station, ma);
        ++aap;
    }
    com /= mass;
    return MassProperties(mass,com,inertia).calcTransformedMassProps(tr);
}

    //////////
    // Body //
    //////////

void Body::realizeTopologicalCache(const DuMMForceFieldSubsystemRep& mm) {
    allAtoms.clear();
    const Cluster& c = mm.getCluster(clusterId);
    AtomPlacementSet::const_iterator ap = c.getAllContainedAtoms().begin();
    while (ap != c.getAllContainedAtoms().end()) {
        allAtoms.push_back(*ap);
        ++ap;
    }
}

} // namespace SimTK

