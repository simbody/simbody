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
 * Contributors: Christopher Bruns, Randy Radmer                              *
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
 *
 * Private implementation of DuMMForceFieldSubsystem. Units here are uniformly
 * MD units: nanometers, daltons, picoseconds, with energy in kilojoules/mole.
 * We accept angles from users in degrees, but use only radians internally.
 */

#include "SimTKcommon.h"

#include "simbody/internal/common.h"
#include "simbody/internal/ForceSubsystem.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/DuMMForceFieldSubsystem.h"
#include "simbody/internal/MolecularMechanicsSystem.h"

#include "ForceSubsystemRep.h"

#include <string>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cassert>
#include <set>
#include <map>
#include <algorithm>

using namespace SimTK;

// This is Coulomb's constant 1/(4*pi*e0) in units which convert
// e^2/nm to kJ/mol.
static const Real CoulombFac = (Real)SimTK_COULOMB_CONSTANT_IN_MD;

template<class T>
class IdPair {
public:
    IdPair() {ids[0]=ids[1]=(T)DuMM::InvalidAtomId;}
    IdPair(T i1, T i2, bool canon=false) {
        ids[0]=i1; ids[1]=i2;
        if (canon) canonicalize();
    }
    T operator[](int i) const {assert(0<=i&&i<2); return ids[i];}
    bool isValid() const {return ids[0]>=0 && ids[1]>=0;}
    // canonical is low,high
    void canonicalize() {if(ids[0]>ids[1]) std::swap(ids[0],ids[1]);}
private:
    T ids[2];
};

template<class T>
inline bool operator<(const IdPair<T>& i1, const IdPair<T>& i2) {
    assert(i1.isValid() && i2.isValid());
    if (i1[0] < i2[0]) return true;
    if (i1[0] > i2[0]) return false;
    return i1[1] < i2[1];
}

typedef IdPair<DuMM::AtomId> AtomIdPair;
typedef IdPair<DuMM::AtomClassId> AtomClassIdPair;

template <class T>
class IdTriple {
public:
    IdTriple() {invalidate();}
    IdTriple(T i1, T i2, T i3, bool canon=false) {
        ids[0]= i1; ids[1]=i2; ids[2]=i3;
        if (canon) canonicalize();
    }
    T operator[](int i) const {assert(0<=i&&i<3); return ids[i];}
    bool isValid() const {return ids[0]>=0 && ids[1]>=0 && ids[2]>=0;}
    void invalidate() {ids[0]=ids[1]=ids[2]=(T)(int)(DuMM::InvalidAtomId);}
    // canonical has 1st number <= last number; middle stays put
    void canonicalize() {if(ids[0]>ids[2]) std::swap(ids[0],ids[2]);}
private:
    T ids[3];
};
template<class T>
inline bool operator<(const IdTriple<T>& i1, const IdTriple<T>& i2) {
    assert(i1.isValid() && i2.isValid());
    if (i1[0] < i2[0]) return true;
    if (i1[0] > i2[0]) return false;
    if (i1[1] < i2[1]) return true;
    if (i1[1] > i2[1]) return false;
    return i1[2] < i2[2];
}

typedef IdTriple<DuMM::AtomId> AtomIdTriple;
typedef IdTriple<DuMM::AtomClassId> AtomClassIdTriple;

template<class T>
class IdQuad {
public:
    IdQuad() {ids[0]=ids[1]=ids[2]=ids[3]=-1;}
    IdQuad(T i1, T i2, T i3, T i4, bool canon=false) {
        ids[0]= i1; ids[1]=i2; ids[2]=i3; ids[3]=i4;
        if (canon) canonicalize();
    }
    T operator[](int i) const {assert(0<=i&&i<4); return ids[i];}
    bool isValid() const {return ids[0]>=0 && ids[1]>=0 && ids[2]>=0 && ids[3]>=0;}
    // canonical has 1st number <= last number; middle two must swap
    // if the outside ones do
    void canonicalize() {
        // Id quad has additional case where 1 == 4 and 2 differs from 3
        if(    (ids[0]>ids[3]) 
            || ( (ids[0] == ids[3]) && (ids[1] > ids[2]) ) )
        {
            std::swap(ids[0],ids[3]); 
            std::swap(ids[1],ids[2]);
        }
    }

private:
    T ids[4];
};

template<class T>
inline bool operator<(const IdQuad<T>& i1, const IdQuad<T>& i2) {
    assert(i1.isValid() && i2.isValid());
    if (i1[0] < i2[0]) return true;
    if (i1[0] > i2[0]) return false;
    if (i1[1] < i2[1]) return true;
    if (i1[1] > i2[1]) return false;
    if (i1[2] < i2[2]) return true;
    if (i1[2] > i2[2]) return false;
    return i1[3] < i2[3];
}

typedef IdQuad<DuMM::AtomId> AtomIdQuad;
typedef IdQuad<DuMM::AtomClassId> AtomClassIdQuad;

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
        // Permit incomplete construction, i.e. radius and depth not yet set
        assert(isValid());
    }

    bool isValid() const {
        return atomClassId >= 0 
            && element > 0 
            && valence >= 0;
    }

    bool isComplete() const {
        return isValid() 
            && vdwRadius >= 0
            && vdwWellDepth >= 0;
    }

    void invalidateTopologicalCache() {
        vdwDij.clear();
        vdwEij.clear();
    }

    void dump() const {
        printf("   %d(%s): element=%d, valence=%d vdwRad=%g nm, vdwDepth(kJ)=%g (%g kcal)\n",
            (int) atomClassId, name.c_str(), element, valence, vdwRadius, vdwWellDepth,
            vdwWellDepth*DuMM::KJ2Kcal);
        printf("    vdwDij (nm):");
        for (int i=0; i< (int)vdwDij.size(); ++i)
            printf(" %g", vdwDij[i]);
        printf("\n    vdwEij (kJ):");
        for (int i=0; i< (int)vdwEij.size(); ++i)
            printf(" %g", vdwEij[i]);
        printf("\n");
    }

    std::ostream& generateSelfCode(std::ostream& os) const 
    {
        os << "defineAtomClass((DuMM::AtomClassId)";
        os << atomClassId << ", ";
        os << "\"" << name << "\", ";
        os << element << ", ";
        os << valence << ", ";
        os << vdwRadius << ", ";
        os << vdwWellDepth << ");";

        return os;
    }

        // TOPOLOGICAL STATE VARIABLES
        //   Filled in during construction.

    DuMM::AtomClassId         atomClassId;
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
    ChargedAtomType() : chargedAtomTypeId(DuMM::InvalidChargedAtomTypeId), atomClassId(DuMM::InvalidAtomId), partialCharge(NaN) { }
    ChargedAtomType(DuMM::ChargedAtomTypeId id, const char* nm, DuMM::AtomClassId aclass, Real chg)
      : chargedAtomTypeId(id), name(nm), atomClassId(aclass), partialCharge(chg) 
    { 
        assert(isValid());
    }
    bool isValid() const {return chargedAtomTypeId >= 0 && atomClassId >= 0;}

    void dump() const {
        printf("    %d(%s): atomClassId=%d, chg=%g e\n", 
              (int) chargedAtomTypeId, name.c_str(), (int) atomClassId, partialCharge);
    }

    std::ostream& generateSelfCode(std::ostream& os) const 
    {
        os << "defineChargedAtomType((DuMM::ChargedAtomTypeId)";
        os << chargedAtomTypeId << ", ";
        os << "\"" << name << "\", ";
        os << "(DuMM::AtomClassId)" << atomClassId << ", ";
        os << partialCharge << ");";

        return os;
    }

    // These are all Topological state variables, filled in during construction.
    // There are no calculations to be performed.
    DuMM::ChargedAtomTypeId         chargedAtomTypeId;
    std::string name;

    DuMM::AtomClassId         atomClassId;
    Real        partialCharge; // qi, in e (charge on proton)

};

// This represents bond-stretch information for a pair of atom types.
// Use an AtomIdPair as a key.
class BondStretch {
public:
    BondStretch() : k(-1), d0(-1), classes(DuMM::InvalidAtomClassId, DuMM::InvalidAtomClassId)
    {}
    BondStretch(AtomClassIdPair key, Real stiffnessInKJperNmSq, Real lengthInNm) 
      : classes(key), k(stiffnessInKJperNmSq), d0(lengthInNm) { 
        assert(isValid());
    }
    bool isValid() const {
        return (k >= 0 )
            && (d0 >= 0)
            && (classes[0] != DuMM::InvalidAtomClassId)
            && (classes[1] != DuMM::InvalidAtomClassId); 
    }

    std::ostream& generateSelfCode(std::ostream& os) const 
    {
        os << "defineBondStretch((DuMM::AtomClassId)";
        os << (int) classes[0] << ", (DuMM::AtomClassId)";
        os << (int) classes[1] << ", ";
        os << k << ", ";
        os << d0  << ");";

        return os;
    }

    AtomClassIdPair classes;
    Real k;  // in energy units (kJ=Da-nm^2/ps^2) per nm^2, i.e. Da/ps^2
    Real d0; // distance at which force is 0 (in nm)
};

class BondBend {
public:
    BondBend() : k(-1), theta0(-1), 
        classes(DuMM::InvalidAtomClassId, DuMM::InvalidAtomClassId, DuMM::InvalidAtomClassId) 
    { }
    BondBend(AtomClassIdTriple key, Real stiffnessInKJPerRadSq, Real angleInDeg) 
      : k(stiffnessInKJPerRadSq), theta0(angleInDeg*DuMM::Deg2Rad),
      classes(key)
    {
        assert(isValid());
    }
    bool isValid() const {return k >= 0 && (0 <= theta0 && theta0 <= Pi);}

    // Given a central atom location c bonded to atoms at r and s,
    // calculate the angle between them, the potential energy,
    // and forces on each of the three atoms.
    void harmonic(const Vec3& cG, const Vec3& rG, const Vec3& sG, const Real& scale,
                  Real& theta, Real& pe, Vec3& cf, Vec3& rf, Vec3& sf) const;

    std::ostream& generateSelfCode(std::ostream& os) const 
    {
        os << "defineBondBend((DuMM::AtomClassId)";
        os << (int) classes[0] << ", (DuMM::AtomClassId)";
        os << (int) classes[1] << ", (DuMM::AtomClassId)";
        os << (int) classes[2] << ", ";
        os << k << ", ";
        os << theta0 << ");";

        return os;
    }

    AtomClassIdTriple classes;
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
      : periodicity(n), amplitude(ampInKJ), theta0(th0InDeg*DuMM::Deg2Rad) {
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

    std::ostream& generateSelfCode(std::ostream& os) const 
    {
        os << ", " << periodicity;
        os << ", " << amplitude;
        os << ", " << theta0 * DuMM::Rad2Deg;
 
        return os;
    }

    int  periodicity; // 1=360, 2=180, 3=120, etc.
    Real amplitude; // energy units (kJ)
    Real theta0;    // radians
};

class BondTorsion {
public:
    BondTorsion() :
      classes(DuMM::InvalidAtomClassId, DuMM::InvalidAtomClassId, DuMM::InvalidAtomClassId, DuMM::InvalidAtomClassId)
    {}
    BondTorsion(AtomClassIdQuad key) : classes(key)
    { }
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

    // equality operator to help handle case where user innocently
    // attempts to add the same torsion a second time
    // WARNING: this is very inefficient
    bool operator==(const BondTorsion& other) const {
        if (terms.size() != other.terms.size()) return false;
        std::vector<TorsionTerm>::const_iterator iTerm;
        for (iTerm = terms.begin(); iTerm != terms.end(); ++iTerm) {
            const TorsionTerm& myTerm = *iTerm;
            if (! other.hasTerm(myTerm.periodicity) ) return false;

            std::vector<TorsionTerm>::const_iterator iOtherTerm;
            for (iOtherTerm = other.terms.begin(); iOtherTerm != other.terms.end(); ++iOtherTerm) {
                const TorsionTerm& otherTerm = *iOtherTerm;
                if (otherTerm.periodicity == myTerm.periodicity) {
                    if (myTerm.amplitude != otherTerm.amplitude) return false;
                    if (myTerm.theta0 != otherTerm.theta0) return false;
                }
            }
            
        }
        return true;
    }

    // Given atom locations r-x-y-s in the ground frame, calculate the
    // torsion angle, energy and a force on each atom so that the desired
    // pure torque is produced.
    void periodic(const Vec3& rG, const Vec3& xG, const Vec3& yG, const Vec3& sG,
                  const Real& scale, Real& theta, Real& pe, 
                  Vec3& rf, Vec3& xf, Vec3& yf, Vec3& sf) const;

    // Type 1 => normal torsion parameters
    // Type 2 => amber improper torsion parameters
    std::ostream& generateSelfCode(std::ostream& os, int torsionType = 1) const 
    {
        if (torsionType == 1)
            os << "defineBondTorsion((DuMM::AtomClassId)";
        else
            os << "defineAmberImproperTorsion((DuMM::AtomClassId)";

        os << classes[0];
        os << ", (DuMM::AtomClassId)" << classes[1];
        os << ", (DuMM::AtomClassId)" << classes[2];
        os << ", (DuMM::AtomClassId)" << classes[3];

        std::vector<TorsionTerm>::const_iterator term;
        for (term = terms.begin(); term != terms.end(); ++term)
            term->generateSelfCode(os);

        os << ");";
 
        return os;
    }

    AtomClassIdQuad classes;
    std::vector<TorsionTerm> terms;
};


class AtomPlacement {
public:
    AtomPlacement() : atomId(-1) { }
    AtomPlacement(DuMM::AtomId a, const Vec3& s) : atomId(a), station(s) {
        assert(isValid());
    }
    bool isValid() const {return atomId >= 0;}

    DuMM::AtomId  atomId;
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
    ClusterPlacement(DuMM::ClusterId c, const Transform& t) : clusterId(c), placement(t) {
        assert(isValid());
    }
    bool isValid() const {return clusterId >= 0;}

    DuMM::ClusterId         clusterId;
    Transform   placement;  // translation in nm
};
inline bool operator<(const ClusterPlacement& r1, const ClusterPlacement& r2) {
    return r1.clusterId < r2.clusterId;
}
inline bool operator==(const ClusterPlacement& r1, const ClusterPlacement& r2) {
    return r1.clusterId == r2.clusterId;
}

typedef std::vector<DuMM::AtomId>            AtomArray;
typedef std::vector<AtomPlacement>  AtomPlacementArray;
typedef std::set<AtomPlacement>     AtomPlacementSet;
typedef std::set<ClusterPlacement>  ClusterPlacementSet;

class Atom {
public:
    Atom() 
      : atomId(-1), chargedAtomTypeId(-1) {
    }
    Atom(DuMM::ChargedAtomTypeId t, DuMM::AtomId aId) : atomId(aId), chargedAtomTypeId(t) {
        assert(isValid());
    }

    bool isValid() const {return atomId>=0 && chargedAtomTypeId>=0;}
    // bool isValid() const {return atomId>=0;}

    bool isAttachedToBody() const {return bodyId >= 0;}

    MobilizedBodyId getBodyId() const {assert(isAttachedToBody()); return bodyId;}

    void attachToBody(MobilizedBodyId bnum, const Vec3& s) {
        assert(!isAttachedToBody());
        bodyId = bnum;
        station_B = s;
    }

    bool isBondedTo(DuMM::AtomId anum) const {
        for (int i=0; i<(int)bond12.size(); ++i)
            if (bond12[i] == anum) return true;
        return false;
    }

    void dump() const;

    void invalidateTopologicalCache() {
        bond13.clear(); bond14.clear(); bond15.clear();
        xbond12.clear(); xbond13.clear(); xbond14.clear(); xbond15.clear();
        shortPath13.clear(); shortPath14.clear(); shortPath15.clear();
        xshortPath13.clear(); xshortPath14.clear(); xshortPath15.clear();
        stretch.clear(); bend.clear(); torsion.clear();
        bonds3Atoms.invalidate();
        xbonds3Atoms.invalidate();
        aImproperTorsion14.clear(); aImproperTorsion.clear();
    }

public:
        // TOPOLOGICAL STATE VARIABLES
        //   Filled in during construction.

    DuMM::AtomId         atomId;
    DuMM::ChargedAtomTypeId         chargedAtomTypeId;
    AtomArray   bond12;

    // After the atom or a containing cluster has been attached to a
    // body, we fill these in.
    MobilizedBodyId bodyId;
    Vec3   station_B; // atom's station fixed in body bodyId's frame, in nm

        // TOPOLOGICAL CACHE ENTRIES
        //   These are calculated in realizeTopology() from topological
        //   state variables (from here or others in the DuMM class).

    // This is a group of lists which identify atoms nearby in the
    // molecule's bond structure. The bond12 list above contains the directly
    // bonded (1-2) atoms; the 13 list below has the 1-(2)-3 bonded atoms (that
    // is, it includes the path to the "3" atom), etc. The current Atom is
    // always atom "1" so it isn't stored.
    //
    // Note that the shortPath and xshortPath arrays give the shortest path between
    // two atoms, while the bond and xbond arrays give *all* connection paths,
    // with bonds3Atoms giving at most one.

    std::vector<AtomIdPair>   bond13;
    std::vector<AtomIdTriple> bond14;
    std::vector<AtomIdQuad>   bond15;
    std::vector<AtomIdPair>   shortPath13;
    std::vector<AtomIdTriple> shortPath14;
    std::vector<AtomIdQuad>   shortPath15;

    // This will be invalid unless we find that the current atom is directly
    // bonded to exactly three other atoms, in which case their atom Ids will
    // be stored here and isValid() will return true.
    AtomIdTriple bonds3Atoms;

    // These are shorter versions of the bond lists in which only those
    // bonds which include atoms from at least two bodies are included.
    // Note that each bond will appear twice in the overall data structure,
    // in the Atom entries for the atoms at either end. We avoid double
    // processing by only processing the instance in which the first atoms'
    // Id is the lower of the two. But we need to keep both copies because
    // these are also used for scaling nearby interaction during non-bonded 
    // calculation.
    // TODO: not sure the above comment about the need for both copies
    // is (a) right in the first place, and (b) in any case necessary for
    // the "bond" arrays since it would seem to apply only to the shortPath
    // arrays which are used for scaling.
    std::vector<DuMM::AtomId>       xbond12;
    std::vector<AtomIdPair>   xbond13;
    std::vector<AtomIdTriple> xbond14;
    std::vector<AtomIdQuad>   xbond15;
    std::vector<AtomIdPair>   xshortPath13;
    std::vector<AtomIdTriple> xshortPath14;
    std::vector<AtomIdQuad>   xshortPath15;

    // This is even less likely to be valid than bonds3Atoms above. It will
    // be valid iff (a) bonds3Atoms is valid, and (b) at least one of the
    // three atoms is on a different body from this one.
    AtomIdTriple xbonds3Atoms;
    std::vector<AtomIdTriple> aImproperTorsion14; // might have zero length
    std::vector<BondTorsion> aImproperTorsion; // might have zero length

    std::vector<BondStretch> stretch; // same length as cross-body 1-2 list
    std::vector<BondBend>    bend;    // same length as   " 1-3 list
    std::vector<BondTorsion> torsion; // same length as   " 1-4 list
};


class Bond {
public:
    Bond() { }
    Bond(DuMM::AtomId atom1, DuMM::AtomId atom2) : atoms(atom1,atom2) { 
        assert(isValid());
    }
    bool isValid() const {return atoms.isValid();}

    AtomIdPair atoms;
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
    Cluster() : clusterId(-1) { }
    Cluster(const char* nm)
        : clusterId(DuMM::InvalidClusterId), name(nm) {
        // not valid yet -- still need Id assigned
    }

    bool isValid() const {return clusterId >= 0;}
    bool isAttachedToBody() const {return bodyId >= 0;}
    bool isTopLevelCluster() const {return parentClusters.empty();}

    MobilizedBodyId getBodyId() const {assert(isAttachedToBody()); return bodyId;}

    const AtomPlacementSet& getDirectlyContainedAtoms() const {return directAtomPlacements;}
    const AtomPlacementSet& getAllContainedAtoms()      const {return allAtomPlacements;}
    AtomPlacementSet&       updAllContainedAtoms()            {return allAtomPlacements;}

    const ClusterPlacementSet& getDirectlyContainedClusters() const {return directClusterPlacements;}
    const ClusterPlacementSet& getAllContainedClusters()      const {return allClusterPlacements;}
    ClusterPlacementSet&       updAllContainedClusters()            {return allClusterPlacements;}

    bool containsAtom(DuMM::AtomId atomId) const {
        return allAtomPlacements.find(AtomPlacement(atomId,Vec3(0))) 
                != allAtomPlacements.end();
    }
    bool containsCluster(DuMM::ClusterId clusterId) const {
        return allClusterPlacements.find(ClusterPlacement(clusterId,Transform())) 
                != allClusterPlacements.end();
    }

    // See if a cluster contains any atoms which are already in
    // any of the cluster trees to which this cluster is associated.
    // TODO: can only handle top-level cluster so we don't have to run up the
    //       ancestor branches.
    // If we find an atom common to both clusters we'll return it to permit
    // nice error messages, otherwise we return false and -1 for the atomId.
    bool overlapsWithCluster(const Cluster& test, DuMM::AtomId& anAtomIdInBothClusters) const {
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
        anAtomIdInBothClusters = DuMM::InvalidAtomId;
        return false;
    }

    // Return true if this cluster contains (directly or indirectly) any atom which has already
    // been attached to a body. If so return one of the attached atoms and its body, which can
    // be helpful in error messages.
    bool containsAnyAtomsAttachedToABody(DuMM::AtomId& atomId, MobilizedBodyId& bodyId, 
                                         const DuMMForceFieldSubsystemRep& mm) const;

    // Translation is in nm.
    void attachToBody(MobilizedBodyId bnum, const Transform& X_BR, DuMMForceFieldSubsystemRep& mm);

    // Place an atom in this cluster. To be valid, the atom must not
    // already be
    //   (a) in any of the trees of which this group is a part, or
    //   (b) attached to a body.
    // TODO: (c) at the moment we don't allow placing an atom in a group unless
    //           that group is a top-level group (i.e., it has no parents).
    // If this group is already attached to a body, then we will update
    // the atom entry to note that it is now attached to the body also.
    void placeAtom(DuMM::AtomId atomId, const Vec3& stationInNm, DuMMForceFieldSubsystemRep& mm);

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
    void placeCluster(DuMM::ClusterId childClusterId, const Transform& placement, DuMMForceFieldSubsystemRep& mm);


    // Calculate the composite mass properties for this cluster, transformed into
    // the indicated frame. Translation part of the Transform is in nm, returned mass
    // proprties are in daltons and nm.
    MassProperties calcMassProperties
       (const Transform& tr, const DuMMForceFieldSubsystemRep& mm) const;


    // Recursively calculate composite properties for this group and all the
    // groups it contains. All groups were marked "invalid" at the beginning
    // of this step.
    void invalidateTopologicalCache() { // TODO
    }
    void realizeTopologicalCache(DuMMForceFieldSubsystemRep& mm) {
    }


    void dump() const {
        printf("    clusterId=%d(%s)\n", (int) clusterId, name.c_str());
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

        if (bodyId >= 0) 
            std::cout << "\n      attached to body " << bodyId << " at (nm) " << placement_B;
        else
            std::cout << "\n      NOT ATTACHED TO ANY BODY.";
        std::cout << std::endl;
    }

    void clearAllCalculatedData() {
        chargeProps    = ChargeProperties();
        geometricProps = GeometricProperties();
    }

private:
    // translation is in nm
    void noteNewChildCluster(DuMM::ClusterId childClusterId, const Transform& X_PC) {
        std::pair<ClusterPlacementSet::iterator, bool> ret;
        ret = directClusterPlacements.insert(ClusterPlacement(childClusterId,X_PC));
        assert(ret.second); // must not have been there already

        ret = allClusterPlacements.insert(ClusterPlacement(childClusterId,X_PC));
        assert(ret.second); // must not have been there already
    }

    // translation is in nm
    void noteNewParentCluster(DuMM::ClusterId parentClusterId, const Transform& X_PC) {
        std::pair<ClusterPlacementSet::iterator, bool> ret =
            parentClusters.insert(ClusterPlacement(parentClusterId,X_PC));
        assert(ret.second); // must not have been there already
    }

public:
        // TOPOLOGICAL STATE VARIABLES
        //   Filled in during construction.
    DuMM::ClusterId                 clusterId;
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
    MobilizedBodyId    bodyId;
    Transform placement_B; // cluster's placement fixed in body bodyId's frame (nm)

        // TOPOLOGICAL CACHE ENTRIES
        //   These are calculated in realizeTopology() from topological
        //   state variables (from here or others in the DuMM class).

    // These reflect composite properties built from the allAtoms list.
    ChargeProperties    chargeProps;
    GeometricProperties geometricProps;
};

// A DuMMBody has a reference to a top-level Cluster, plus some information used
// at runtime for fast body-by-body processing.
class DuMMBody {
public:
    DuMMBody() : clusterId(DuMM::InvalidClusterId) { }
    explicit DuMMBody(DuMM::ClusterId cId) : clusterId(cId) { 
        assert(isValid());
    }

    bool isValid() const {return clusterId >= 0;}

    void invalidateTopologicalCache() {allAtoms.clear();}
    void realizeTopologicalCache(const DuMMForceFieldSubsystemRep& mm);

    DuMM::ClusterId getClusterId() const {assert(isValid()); return clusterId;}

    void dump() const {
        printf("    clusterId=%d\n", (int) clusterId);
        printf("    shadowBodies=");
        for (int i=0; i < (int)shadowBodies.size(); ++i)
            printf(" %d", shadowBodies[i]);
        printf("\n");
        printf("    allAtoms=");
        for (int i=0; i < (int)allAtoms.size(); ++i) 
            printf(" %d(%g,%g,%g)(nm)", (int) allAtoms[i].atomId,
                allAtoms[i].station[0], allAtoms[i].station[1], allAtoms[i].station[2]);
        printf("\n");
    }

    static std::string createClusterNameForBody(int bnum) {
        char buf[100];
        std::sprintf(buf, "DuMMBody %d", bnum);
        return std::string(buf);
    }

    DuMM::ClusterId clusterId;
    std::vector<int> shadowBodies; // if needed

    // This is an expansion of all the atom & group placements, with
    // all stations transformed to this body's frame, sorted in order
    // of atomId, and built for speed!
    AtomPlacementArray  allAtoms;
};

class SimTK::DuMMForceFieldSubsystemRep : public ForceSubsystemRep {
    friend class DuMMForceFieldSubsystem;
    static const char* ApiClassName; // "DuMMForceFieldSubsystem"
public:
    DuMMForceFieldSubsystemRep()
      : ForceSubsystemRep("DuMMForceFieldSubsystem", "0.0.1")
    {
        vdwMixingRule = DuMMForceFieldSubsystem::WaldmanHagler;
        vdwGlobalScaleFactor=coulombGlobalScaleFactor=bondStretchGlobalScaleFactor
            =bondBendGlobalScaleFactor=bondTorsionGlobalScaleFactor
            =amberImproperTorsionGlobalScaleFactor=1;
        vdwScale12=coulombScale12=vdwScale13=coulombScale13=0;
        vdwScale14=coulombScale14=vdwScale15=coulombScale15=1;
        loadElements();
        const DuMM::ClusterId gid = addCluster(Cluster("free atoms and groups"));
        assert(gid==0);
    }

    // common checks when defining improper and proper torsions
    void DuMMForceFieldSubsystemRep::checkTorsion 
        (DuMM::AtomClassId class1, DuMM::AtomClassId class2, DuMM::AtomClassId class3, DuMM::AtomClassId class4, 
            int periodicity1, Real amp1InKJ, Real phase1InDegrees,
            int periodicity2, Real amp2InKJ, Real phase2InDegrees,
            int periodicity3, Real amp3InKJ, Real phase3InDegrees,
            const char* CallingMethodName) const;

    bool isValidElement(int atomicNumber) const {
        return 1 <= atomicNumber && atomicNumber < (int)elements.size() 
                && elements[atomicNumber].isValid();
    }

    bool isValidAtom(DuMM::AtomId atomNum) const {
        return 0 <= atomNum && atomNum < (DuMM::AtomId)atoms.size() && atoms[atomNum].isValid();
    }

    bool isValidBond(DuMM::BondId bondNum) const {
        return 0 <= bondNum && bondNum < (DuMM::BondId)bonds.size() && bonds[bondNum].isValid();
    }

    bool isValidCluster(DuMM::ClusterId clusterId) const {
        return 0 <= clusterId && clusterId < (DuMM::ClusterId)clusters.size()
                && clusters[clusterId].isValid();
    }

    bool isValidBody(int bodyId) const {
        return 0 <= bodyId && bodyId < (int)bodies.size() && bodies[bodyId].isValid();
    }

    bool isValidChargedAtomType(DuMM::ChargedAtomTypeId typeNum) const {
        return 0 <= typeNum && typeNum < (DuMM::ChargedAtomTypeId)chargedAtomTypes.size() 
               && chargedAtomTypes[typeNum].isValid();
    }

    bool isValidAtomClass(DuMM::AtomClassId classNum) const {
        return 0 <= classNum && classNum < (DuMM::AtomClassId)atomClasses.size() 
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

    DuMM::ClusterId addCluster(const Cluster& c) {

        invalidateSubsystemTopologyCache();

        const DuMM::ClusterId clusterId = (DuMM::ClusterId)clusters.size();
        clusters.push_back(c);
        clusters[clusterId].clusterId = clusterId;
        return clusterId;
    }
    Cluster& updCluster(DuMM::ClusterId clusterId) {
        assert(isValidCluster(clusterId));

        invalidateSubsystemTopologyCache();
        return clusters[clusterId];
    }
    const Cluster& getCluster(DuMM::ClusterId clusterId) const {
        assert(isValidCluster(clusterId));
        return clusters[clusterId];
    }
    DuMMBody& updBody(int bodyId) {
        assert(isValidBody(bodyId));

        invalidateSubsystemTopologyCache();
        return bodies[bodyId];
    }
    const DuMMBody& getBody(int bodyId) const {
        assert(isValidBody(bodyId));
        return bodies[bodyId];
    }


    int getNAtoms() const {return (int)atoms.size();}
    int getNBonds() const {return (int)bonds.size();}

    const Atom& getAtom(DuMM::AtomId atomId) const {
        assert(isValidAtom(atomId));
        return atoms[atomId];
    }
    Atom& updAtom(DuMM::AtomId atomId) {
        assert(isValidAtom(atomId));

        invalidateSubsystemTopologyCache();
        return atoms[atomId];
    }

    DuMM::ChargedAtomTypeId getChargedAtomTypeId(DuMM::AtomId atomId) const {
        return getAtom(atomId).chargedAtomTypeId;
    }

    DuMM::AtomClassId getAtomClassId(DuMM::AtomId atomId) const {
        const ChargedAtomType& type = chargedAtomTypes[getChargedAtomTypeId(atomId)];
        return type.atomClassId;
    }

    int getAtomElementNum(DuMM::AtomId atomId) const {
        const AtomClass& cl = atomClasses[getAtomClassId(atomId)];
        return cl.element;
    }

    const Element& getElement(int element) const {
        assert(isValidElement(element));
        return elements[element];
    }


    const BondStretch& getBondStretch(DuMM::AtomClassId class1, DuMM::AtomClassId class2) const;
    const BondBend&    getBondBend   (DuMM::AtomClassId class1, DuMM::AtomClassId class2, DuMM::AtomClassId class3) const;
    const BondTorsion& getBondTorsion(DuMM::AtomClassId class1, DuMM::AtomClassId class2, DuMM::AtomClassId class3, DuMM::AtomClassId class4) const;
    const BondTorsion& getAmberImproperTorsion(DuMM::AtomClassId class1, DuMM::AtomClassId class2, DuMM::AtomClassId class3, DuMM::AtomClassId class4) const;

    // Override virtual methods from Subsystem::Guts class.

    DuMMForceFieldSubsystemRep* cloneImpl() const {
        return new DuMMForceFieldSubsystemRep(*this);
    }

    int realizeSubsystemTopologyImpl(State& s) const;

    int realizeSubsystemModelImpl(State& s) const {
        // Sorry, no choices available at the moment.
        return 0;
    }

    int realizeSubsystemInstanceImpl(const State& s) const {
        // Nothing to compute here.
        return 0;
    }

    int realizeSubsystemTimeImpl(const State& s) const {
        // Nothing to compute here.
        return 0;
    }

    int realizeSubsystemPositionImpl(const State& s) const {
        // Nothing to compute here.
        return 0;
    }

    int realizeSubsystemVelocityImpl(const State& s) const {
        // Nothing to compute here.
        return 0;
    }


    int realizeSubsystemDynamicsImpl(const State& s) const;

    int realizeSubsystemAccelerationImpl(const State& s) const {
        // Nothing to compute here.
        return 0;
    }

    int realizeSubsystemReportImpl(const State& s) const {
        // Nothing to compute here.
        return 0;
    }

    void dump() const;


private:
    void loadElements();

    void ensureBodyEntryExists(MobilizedBodyId bodyNum) {
        if (bodyNum >= (int)bodies.size())
            bodies.resize(bodyNum+1);
        if (!bodies[bodyNum].isValid()) {
            const DuMM::ClusterId clusterId = 
                addCluster(Cluster(DuMMBody::createClusterNameForBody(bodyNum).c_str()));
            clusters[clusterId].attachToBody(bodyNum, Transform(), *this);
            bodies[bodyNum] = DuMMBody(clusterId);
        }
    }

    void invalidateAllTopologicalCacheEntries() {
        // If any of these objects are invalid, the invalidateTopologicalCache()
        // call does nothing (i.e., it doesn't blow up!).

        // molecule
        for (DuMM::AtomId i = (DuMM::AtomId)0; i < (DuMM::AtomId)atoms.size(); ++i)
            atoms[i].invalidateTopologicalCache();
        for (DuMM::ClusterId i = (DuMM::ClusterId)0; i < (DuMM::ClusterId)clusters.size(); ++i)
            clusters[i].invalidateTopologicalCache();
        for (int i=0; i < (int)bodies.size(); ++i)
            bodies[i].invalidateTopologicalCache();

        // force field
        for (DuMM::AtomClassId i = (DuMM::AtomClassId)0; i < (DuMM::AtomClassId)atomClasses.size(); ++i)
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
    std::vector<DuMMBody>    bodies;

    // force field

    // Force field description. These are not necessarily fully populated;
    // check the "isValid()" method to see if anything is there.
    std::vector<Element>             elements;
    std::vector<AtomClass>           atomClasses;
    std::vector<ChargedAtomType>     chargedAtomTypes;

    // These relate atom classes, not charged atom types.
    std::map<AtomClassIdPair,   BondStretch> bondStretch;
    std::map<AtomClassIdTriple, BondBend>    bondBend;
    std::map<AtomClassIdQuad,   BondTorsion> bondTorsion;
    std::map<AtomClassIdQuad,   BondTorsion> amberImproperTorsion;

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
         bondTorsionGlobalScaleFactor, amberImproperTorsionGlobalScaleFactor;

        // TOPOLOGICAL CACHE ENTRIES
        //   These are calculated in realizeTopology() from topological
        //   state variables (from here or others in the DuMM class).
};


    ////////////////////////////////
    // DUMM FORCE FIELD SUBSYSTEM //
    ////////////////////////////////

/*static*/ bool 
DuMMForceFieldSubsystem::isInstanceOf(const Subsystem& s) {
    return DuMMForceFieldSubsystemRep::isA(s.getSubsystemGuts());
}
/*static*/ const DuMMForceFieldSubsystem&
DuMMForceFieldSubsystem::downcast(const Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const DuMMForceFieldSubsystem&>(s);
}
/*static*/ DuMMForceFieldSubsystem&
DuMMForceFieldSubsystem::updDowncast(Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<DuMMForceFieldSubsystem&>(s);
}

const DuMMForceFieldSubsystemRep& 
DuMMForceFieldSubsystem::getRep() const {
    return dynamic_cast<const DuMMForceFieldSubsystemRep&>(ForceSubsystem::getRep());
}
DuMMForceFieldSubsystemRep&       
DuMMForceFieldSubsystem::updRep() {
    return dynamic_cast<DuMMForceFieldSubsystemRep&>(ForceSubsystem::updRep());
}

// Create Subsystem but don't associate it with any System. This isn't much use except
// for making std::vector's, which require a default constructor to be available.
DuMMForceFieldSubsystem::DuMMForceFieldSubsystem() 
  : ForceSubsystem()
{
    adoptSubsystemGuts(new DuMMForceFieldSubsystemRep());
}

DuMMForceFieldSubsystem::DuMMForceFieldSubsystem(MolecularMechanicsSystem& mms) 
  : ForceSubsystem()
{
    adoptSubsystemGuts(new DuMMForceFieldSubsystemRep());
    mms.setMolecularMechanicsForceSubsystem(*this); // steal ownership
}

void DuMMForceFieldSubsystem::dumpCForcefieldParameters(std::ostream& os, const String& methodName) const {
    const DuMMForceFieldSubsystemRep& mm = getRep();

    os << "void " << methodName << "(DuMMForceFieldSubsystem& dumm)" << std::endl;
    os << "{" << std::endl; // open method

    // 1) define atom classes
    for (int i=0; i < (int)mm.atomClasses.size(); ++i) {
        if (!mm.atomClasses[i].isValid()) continue;
        const AtomClass& atomClass = mm.atomClasses[i];

        os << "    dumm.";
        atomClass.generateSelfCode(os);
        os << std::endl;
    }

    os << std::endl;

    // 2) define charged atom types
    for (int i=0; i < (int)mm.chargedAtomTypes.size(); ++i) {
        if (!mm.chargedAtomTypes[i].isValid()) continue;

        const ChargedAtomType& chargedAtomType = mm.chargedAtomTypes[i];
        os << "    dumm.";
        chargedAtomType.generateSelfCode(os);
        os << std::endl;
    }

    os << std::endl;

    // 3) bond stretch parameters
    std::map<AtomClassIdPair, BondStretch>::const_iterator b;
    for (b = mm.bondStretch.begin(); b != mm.bondStretch.end(); ++b) {
        os << "    dumm.";
        b->second.generateSelfCode(os);
        os << std::endl;
    }

    os << std::endl;

    // 4) bond torsion parameters
    std::map<AtomClassIdQuad, BondTorsion>::const_iterator t;
    for (t = mm.bondTorsion.begin(); t != mm.bondTorsion.end(); ++t) {
        os << "    dumm.";
        t->second.generateSelfCode(os);
        os << std::endl;
    }

    os << std::endl;

    // 5) amber-style improper torsion parameters
    for (t = mm.amberImproperTorsion.begin(); t != mm.amberImproperTorsion.end(); ++t) {
        os << "    dumm.";
        t->second.generateSelfCode(os, 2);
        os << std::endl;
    }

    os << std::endl;

    // 6) global parameters

    // van der Waals mixing rule
    os << "    dumm.setVdwMixingRule(";
    switch (getVdwMixingRule()) {
        case WaldmanHagler:
            os << "DuMMForceFieldSubsystem::WaldmanHagler";
            break;
        case HalgrenHHG:
            os << "DuMMForceFieldSubsystem::HalgrenHHG";
            break;
        case Jorgensen:
            os << "DuMMForceFieldSubsystem::Jorgensen";
            break;
        case LorentzBerthelot:
            os << "DuMMForceFieldSubsystem::LorentzBerthelot";
            break;
        case Kong:
            os << "DuMMForceFieldSubsystem::Kong";
            break;
        default:
            assert(false);
            os << "DuMMForceFieldSubsystem::WaldmanHagler";
            break;
    }
    os << ");" << std::endl;

    os << "    dumm.setVdw12ScaleFactor(" << mm.vdwScale12 << ");" << std::endl;
    os << "    dumm.setVdw13ScaleFactor(" << mm.vdwScale13 << ");" << std::endl;
    os << "    dumm.setVdw14ScaleFactor(" << mm.vdwScale14 << ");" << std::endl;
    os << "    dumm.setVdw15ScaleFactor(" << mm.vdwScale15 << ");" << std::endl;

    os << "    dumm.setCoulomb12ScaleFactor(" << mm.coulombScale12 << ");" << std::endl;
    os << "    dumm.setCoulomb13ScaleFactor(" << mm.coulombScale13 << ");" << std::endl;
    os << "    dumm.setCoulomb14ScaleFactor(" << mm.coulombScale14 << ");" << std::endl;
    os << "    dumm.setCoulomb15ScaleFactor(" << mm.coulombScale15 << ");" << std::endl;

    os << "    dumm.setVdwGlobalScaleFactor(" << mm.vdwGlobalScaleFactor << ");" << std::endl;
    os << "    dumm.setCoulombGlobalScaleFactor(" << mm.coulombGlobalScaleFactor << ");" << std::endl;
    os << "    dumm.setBondStretchGlobalScaleFactor(" << mm.bondStretchGlobalScaleFactor << ");" << std::endl;
    os << "    dumm.setBondBendGlobalScaleFactor(" << mm.bondBendGlobalScaleFactor << ");" << std::endl;
    os << "    dumm.setBondTorsionGlobalScaleFactor(" << mm.bondTorsionGlobalScaleFactor << ");" << std::endl;
    os << "    dumm.setAmberImproperTorsionGlobalScaleFactor(" << mm.amberImproperTorsionGlobalScaleFactor << ");" << std::endl;

    os << "}" << std::endl; // end of method
}

void DuMMForceFieldSubsystem::defineIncompleteAtomClass
   (DuMM::AtomClassId atomClassId, const char* atomClassName, int element, int valence)
{
    static const char* MethodName = "defineIncompleteAtomClass";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

        // Catch nonsense arguments.
    SimTK_APIARGCHECK1_ALWAYS(atomClassId >= 0, mm.ApiClassName, MethodName,
        "atom class Id %d invalid: must be nonnegative", (int) atomClassId);
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidElement(element), mm.ApiClassName, MethodName,
        "element %d invalid: must be a valid atomic number and have an entry here",element);
    SimTK_APIARGCHECK1_ALWAYS(valence >= 0, mm.ApiClassName, MethodName, 
        "expected valence %d invalid: must be nonnegative", valence);

        // Make sure there is a slot available for this atom class.
    if (atomClassId >= (DuMM::AtomClassId)mm.atomClasses.size())
        mm.atomClasses.resize(atomClassId+1);

        // Make sure this atom class hasn't already been defined.
    SimTK_APIARGCHECK2_ALWAYS(!mm.atomClasses[atomClassId].isValid(), mm.ApiClassName, MethodName, 
        "atom class Id %d is already in use for '%s'", (int) atomClassId, 
        mm.atomClasses[atomClassId].name.c_str());

        // It's all good -- add the new atom class.
    mm.atomClasses[atomClassId] = AtomClass(atomClassId, atomClassName, element, valence, 
                                            NaN, NaN);
}

void DuMMForceFieldSubsystem::setAtomClassVdwParameters(DuMM::AtomClassId atomClassId, Real vdwRadiusInNm, Real vdwWellDepthInKJPerMol) 
{
    static const char* MethodName = "setAtomClsasVdwParameters";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(atomClassId >= 0, mm.ApiClassName, MethodName,
        "atom class Id %d invalid: must be nonnegative", (int) atomClassId);

    SimTK_APIARGCHECK1_ALWAYS(vdwRadiusInNm >= 0, mm.ApiClassName, MethodName, 
        "van der Waals radius %g invalid: must be nonnegative", vdwRadiusInNm);
    SimTK_APIARGCHECK1_ALWAYS(vdwWellDepthInKJPerMol >= 0, mm.ApiClassName, MethodName, 
        "van der Waals energy well depth %g invalid: must be nonnegative", vdwWellDepthInKJPerMol);

    AtomClass& atomClass = mm.atomClasses[atomClassId];
    atomClass.vdwRadius = vdwRadiusInNm;
    atomClass.vdwWellDepth = vdwWellDepthInKJPerMol;
}

bool DuMMForceFieldSubsystem::isValidAtomClass(DuMM::AtomClassId atomClassId) const {
    return getRep().isValidAtomClass(atomClassId);
}

void DuMMForceFieldSubsystem::defineIncompleteChargedAtomType
    (DuMM::ChargedAtomTypeId chargedAtomTypeId, const char* typeName, DuMM::AtomClassId atomClassId)
{
    static const char* MethodName = "defineChargedAtomType";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

        // Check for nonsense arguments.
    SimTK_APIARGCHECK1_ALWAYS(chargedAtomTypeId >= 0, mm.ApiClassName, MethodName,
        "charged atom type Id %d invalid: must be nonnegative", (int) chargedAtomTypeId);
    SimTK_APIARGCHECK1_ALWAYS(atomClassId >= 0, mm.ApiClassName, MethodName,
        "atom class Id %d invalid: must be nonnegative", (int) atomClassId);
    // partialCharge is a signed quantity

        // Make sure the referenced atom class has already been defined.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtomClass(atomClassId), mm.ApiClassName, MethodName,
        "atom class %d is undefined", (int) atomClassId);

        // Make sure there is a slot available for the new chargedAtomType.
    if (chargedAtomTypeId >= (int)mm.chargedAtomTypes.size())
        mm.chargedAtomTypes.resize(chargedAtomTypeId+1);

        // Check that this slot is not already in use.
    SimTK_APIARGCHECK2_ALWAYS(!mm.chargedAtomTypes[chargedAtomTypeId].isValid(), mm.ApiClassName, MethodName, 
        "charged atom type Id %d is already in use for '%s'", (int) chargedAtomTypeId, 
        mm.chargedAtomTypes[chargedAtomTypeId].name.c_str());

        // Define the new charged atom type.
    mm.chargedAtomTypes[chargedAtomTypeId] = 
        ChargedAtomType(chargedAtomTypeId, typeName, atomClassId, NaN);
}

void DuMMForceFieldSubsystem::setChargedAtomTypeCharge(DuMM::ChargedAtomTypeId chargedAtomTypeId, Real charge) {
    static const char* MethodName = "defineChargedAtomType";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

        // Check for nonsense arguments.
    SimTK_APIARGCHECK1_ALWAYS(chargedAtomTypeId >= 0, mm.ApiClassName, MethodName,
        "charged atom type Id %d invalid: must be nonnegative", (int) chargedAtomTypeId);

    ChargedAtomType& chargedAtomType = mm.chargedAtomTypes[chargedAtomTypeId];
    chargedAtomType.partialCharge = charge;
}

void DuMMForceFieldSubsystem::defineBondStretch
   (DuMM::AtomClassId class1, DuMM::AtomClassId class2, Real stiffnessInKJPerNmSq, Real nominalLengthInNm)
{
    static const char* MethodName = "defineBondStretch";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

        // Watch for nonsense arguments.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtomClass(class1), mm.ApiClassName, MethodName, 
        "class1=%d which is not a valid atom class Id", (int) class1);
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtomClass(class2), mm.ApiClassName, MethodName, 
        "class2=%d which is not a valid atom class Id", (int) class2);
    SimTK_APIARGCHECK1_ALWAYS(stiffnessInKJPerNmSq >= 0, mm.ApiClassName, MethodName, 
        "stiffness %g is not valid: must be nonnegative", stiffnessInKJPerNmSq);
    SimTK_APIARGCHECK1_ALWAYS(nominalLengthInNm >= 0, mm.ApiClassName, MethodName, 
        "nominal length %g is not valid: must be nonnegative", nominalLengthInNm);

        // Attempt to insert the new bond stretch entry, canonicalizing first
        // so that the atom class pair has the lower class Id first.
    const AtomClassIdPair key(class1,class2,true);
    std::pair<std::map<AtomClassIdPair,BondStretch>::iterator, bool> ret = 
      mm.bondStretch.insert(std::pair<AtomClassIdPair,BondStretch>
        (key, BondStretch(key,stiffnessInKJPerNmSq,nominalLengthInNm)));

        // Throw an exception if this bond stretch term was already defined. (std::map 
        // indicates that with a bool in the return value.)
    SimTK_APIARGCHECK2_ALWAYS(ret.second, mm.ApiClassName, MethodName, 
        "there was already a bond stretch term for atom class pair (%d,%d)", (int) key[0], (int) key[1]);
}

void DuMMForceFieldSubsystem::defineBondBend
   (DuMM::AtomClassId class1, DuMM::AtomClassId class2, DuMM::AtomClassId class3, Real stiffnessInKJPerRadSq, Real nominalAngleInDeg)
{
    static const char* MethodName = "defineBondBend";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

        // Watch for nonsense arguments.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtomClass(class1), mm.ApiClassName, MethodName, 
        "class1=%d which is not a valid atom class Id", (int) class1);
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtomClass(class2), mm.ApiClassName, MethodName, 
        "class2=%d which is not a valid atom class Id", (int) class2);
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtomClass(class3), mm.ApiClassName, MethodName, 
        "class3=%d which is not a valid atom class Id", (int) class3);
    SimTK_APIARGCHECK1_ALWAYS(stiffnessInKJPerRadSq >= 0, mm.ApiClassName, MethodName, 
        "stiffness %g is not valid: must be nonnegative", stiffnessInKJPerRadSq);
    SimTK_APIARGCHECK1_ALWAYS(0 <= nominalAngleInDeg && nominalAngleInDeg <= 180, 
        mm.ApiClassName, MethodName, 
        "nominal angle %g is not valid: must be between 0 and 180 degrees, inclusive", 
        nominalAngleInDeg);

        // Attempt to insert the new bond bend entry, canonicalizing first
        // by reversing the class Id triple if necessary so that the first 
        // classId is no larger than the third.
    const AtomClassIdTriple key(class1, class2, class3, true);
    std::pair<std::map<AtomClassIdTriple,BondBend>::iterator, bool> ret = 
      mm.bondBend.insert(std::pair<AtomClassIdTriple,BondBend>
        (key, BondBend(key, stiffnessInKJPerRadSq,nominalAngleInDeg)));

        // Throw an exception if this bond bend term was already defined. (std::map 
        // indicates that with a bool in the return value.)
    SimTK_APIARGCHECK3_ALWAYS(ret.second, mm.ApiClassName, MethodName, 
        "there was already a bond bend term for atom class triple (%d,%d,%d)", 
        (int) key[0], (int) key[1], (int) key[2]);
}


// 
// This is a utility method that checks for invalid inputs to the defineBondTorsion() and
// defineAmberImproperTorsion() functions.
//
void DuMMForceFieldSubsystemRep::checkTorsion 
(DuMM::AtomClassId class1, DuMM::AtomClassId class2, DuMM::AtomClassId class3, DuMM::AtomClassId class4, 
    int periodicity1, Real amp1InKJ, Real phase1InDegrees,
    int periodicity2, Real amp2InKJ, Real phase2InDegrees,
    int periodicity3, Real amp3InKJ, Real phase3InDegrees,
    const char* CallingMethodName) const
{
    // DuMMForceFieldSubsystemRep& mm = updRep();

        // Watch for nonsense arguments.
    SimTK_APIARGCHECK1_ALWAYS(isValidAtomClass(class1), ApiClassName, CallingMethodName,
        "class1=%d which is not a valid atom class Id", (int) class1);
    SimTK_APIARGCHECK1_ALWAYS(isValidAtomClass(class2), ApiClassName, CallingMethodName,
        "class2=%d which is not a valid atom class Id", (int) class2);
    SimTK_APIARGCHECK1_ALWAYS(isValidAtomClass(class3), ApiClassName, CallingMethodName,
        "class3=%d which is not a valid atom class Id", (int) class3);
    SimTK_APIARGCHECK1_ALWAYS(isValidAtomClass(class4), ApiClassName, CallingMethodName,
        "class4=%d which is not a valid atom class Id", (int) class4);
    SimTK_APIARGCHECK_ALWAYS(periodicity1!=-1 || periodicity2!=-1 || periodicity3!=-1,
        ApiClassName, CallingMethodName, "must be at least one torsion term supplied");


    if (periodicity1 != -1) {
            // No nonsense.
        SimTK_APIARGCHECK1_ALWAYS(1 <= periodicity1 && periodicity1 <= 6, ApiClassName, CallingMethodName,
            "periodicity1(%d) is invalid: we require 1 <= periodicity <= 6", periodicity1);
        SimTK_APIARGCHECK1_ALWAYS(amp1InKJ >= 0, ApiClassName, CallingMethodName,
            "amplitude1(%g) is not valid: must be nonnegative", amp1InKJ);
        SimTK_APIARGCHECK1_ALWAYS(0 <= phase1InDegrees && phase1InDegrees <= 180, ApiClassName, CallingMethodName,
            "phaseAngle1(%g) is not valid: must be between 0 and 180 degrees, inclusive", phase1InDegrees);

            // No repeats.
        SimTK_APIARGCHECK1_ALWAYS((periodicity2 != periodicity1) && (periodicity3 != periodicity1),
            ApiClassName, CallingMethodName,
            "only one term with a given periodicity may be specified (periodicity %d was repeated)",
            periodicity1);
    }
    if (periodicity2 != -1) {
            // No nonsense.
        SimTK_APIARGCHECK1_ALWAYS(1 <= periodicity2 && periodicity2 <= 6, ApiClassName, CallingMethodName,
            "periodicity2(%d) is invalid: we require 1 <= periodicity <= 6", periodicity2);
        SimTK_APIARGCHECK1_ALWAYS(amp2InKJ >= 0, ApiClassName, CallingMethodName,
            "amplitude2(%g) is not valid: must be nonnegative", amp2InKJ);
        SimTK_APIARGCHECK1_ALWAYS(0 <= phase2InDegrees && phase2InDegrees <= 180, ApiClassName, CallingMethodName,
            "phaseAngle2(%g) is not valid: must be between 0 and 180 degrees, inclusive", phase2InDegrees);

            // No repeats.
        SimTK_APIARGCHECK1_ALWAYS(periodicity3 != periodicity2, ApiClassName, CallingMethodName,
            "only one term with a given periodicity may be specified (periodicity %d was repeated)",
            periodicity2);
    }
    if (periodicity3 != -1) {
            // No nonsense.
        SimTK_APIARGCHECK1_ALWAYS(1 <= periodicity3 && periodicity3 <= 6, ApiClassName, CallingMethodName,
            "periodicity3(%d) is invalid: we require 1 <= periodicity <= 6", periodicity3);
        SimTK_APIARGCHECK1_ALWAYS(amp3InKJ >= 0, ApiClassName, CallingMethodName,
            "amplitude3(%g) is not valid: must be nonnegative", amp3InKJ);
        SimTK_APIARGCHECK1_ALWAYS(0 <= phase3InDegrees && phase3InDegrees <= 180, ApiClassName, CallingMethodName,
            "phaseAngle3(%g) is not valid: must be between 0 and 180 degrees, inclusive", phase3InDegrees);
            // (we've already checked for any possible repeats)
    }
}



// 
// We allow up to 3 terms in a single torsion function, with three different
// periodicities. If any of these are unused, set the corresponding periodicity
// to -1.
//
void DuMMForceFieldSubsystem::defineBondTorsion
   (DuMM::AtomClassId class1, DuMM::AtomClassId class2, DuMM::AtomClassId class3, DuMM::AtomClassId class4, 
    int periodicity1, Real amp1InKJ, Real phase1InDegrees,
    int periodicity2, Real amp2InKJ, Real phase2InDegrees,
    int periodicity3, Real amp3InKJ, Real phase3InDegrees)
{
    static const char* MethodName = "defineBondTorsion";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();
    mm.checkTorsion(class1, class2, class3, class4, 
                 periodicity1, amp1InKJ, phase1InDegrees,
                 periodicity2, amp2InKJ, phase2InDegrees,
                 periodicity3, amp3InKJ, phase3InDegrees,
                 MethodName);

        // Canonicalize atom class quad by reversing order if necessary so that the
        // first class Id is numerically no larger than the fourth.
    const AtomClassIdQuad key(class1, class2, class3, class4, true);

        // Now allocate an empty BondTorsion object and add terms to it as they are found.
    BondTorsion bt(key);
    if (periodicity1 != -1) {
         bt.addTerm(TorsionTerm(periodicity1, amp1InKJ, phase1InDegrees));
    }
    if (periodicity2 != -1) {
        bt.addTerm(TorsionTerm(periodicity2, amp2InKJ, phase2InDegrees));
    }
    if (periodicity3 != -1) {
        bt.addTerm(TorsionTerm(periodicity3, amp3InKJ, phase3InDegrees));
    }

    // If this torsion is already defined, this should ordinarily be an error
    // But, if the parameters are the same, let it slide
    if (mm.bondTorsion.find(key) != mm.bondTorsion.end()) {
        const BondTorsion& oldBondTorsion = mm.bondTorsion.find(key)->second;
        if (oldBondTorsion == bt) return;  // same, so let it slide
    }

        // Now try to insert the allegedly new BondTorsion specification into the bondTorsion map.
        // If it is already there the 2nd element in the returned pair will be 'false'.
    std::pair<std::map<AtomClassIdQuad,BondTorsion>::iterator, bool> ret = 
      mm.bondTorsion.insert(std::pair<AtomClassIdQuad,BondTorsion>(key,bt));

        // Throw an exception if terms for this bond torsion were already defined.
    SimTK_APIARGCHECK4_ALWAYS(ret.second, mm.ApiClassName, MethodName, 
        "bond torsion term(s) were already defined for atom class quad (%d,%d,%d,%d)", 
        (int)key[0], (int)key[1], (int)key[2], (int)key[3]);
}

// Convenient signature for a bond torsion with only one term.
void DuMMForceFieldSubsystem::defineBondTorsion
   (DuMM::AtomClassId class1, DuMM::AtomClassId class2, DuMM::AtomClassId class3, DuMM::AtomClassId class4, 
    int periodicity1, Real amp1InKJ, Real phase1InDegrees)
{
    defineBondTorsion(class1, class2, class3, class4, 
                      periodicity1,amp1InKJ,phase1InDegrees,
                      -1,0.,0., -1,0.,0.);
}

// Convenient signature for a bond torsion with two terms.
void DuMMForceFieldSubsystem::defineBondTorsion
   (DuMM::AtomClassId class1, DuMM::AtomClassId class2, DuMM::AtomClassId class3, DuMM::AtomClassId class4, 
    int periodicity1, Real amp1InKJ, Real phase1InDegrees,
    int periodicity2, Real amp2InKJ, Real phase2InDegrees)
{
    defineBondTorsion(class1, class2, class3, class4, 
                      periodicity1,amp1InKJ,phase1InDegrees,
                      periodicity2,amp2InKJ,phase2InDegrees,
                      -1,0.,0.);
}

// 
// This function is based on the defineTorsion function.
// As with the normal bond torsions, we allow up to 3 terms in a single torsion function,
// with three different periodicities. If any of these are unused, set the corresponding
// periodicity to -1.
//
void DuMMForceFieldSubsystem::defineAmberImproperTorsion
(DuMM::AtomClassId class1, DuMM::AtomClassId class2, DuMM::AtomClassId class3, DuMM::AtomClassId class4,
    int periodicity1, Real amp1InKJ, Real phase1InDegrees,
    int periodicity2, Real amp2InKJ, Real phase2InDegrees,
    int periodicity3, Real amp3InKJ, Real phase3InDegrees)
{
    static const char* MethodName = "defineAmberImproperTorsion";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();
    mm.checkTorsion(class1, class2, class3, class4, 
                 periodicity1, amp1InKJ, phase1InDegrees,
                 periodicity2, amp2InKJ, phase2InDegrees,
                 periodicity3, amp3InKJ, phase3InDegrees,
                 MethodName);

        // Unlike the normal bond torsions (see defineBondTorstion function) we do *not*
        // canonicalize atom class quad, because atom order does matter for amber improper torsions
    const AtomClassIdQuad key(class1, class2, class3, class4, false);

        // Now allocate an empty BondTorsion object and add terms to it as they are found.
    BondTorsion bt(key);
    // Add the new terms.
    if (periodicity1 != -1)
        bt.addTerm(TorsionTerm(periodicity1, amp1InKJ, phase1InDegrees));
    if (periodicity2 != -1)
        bt.addTerm(TorsionTerm(periodicity2, amp2InKJ, phase2InDegrees));
    if (periodicity3 != -1)
        bt.addTerm(TorsionTerm(periodicity3, amp3InKJ, phase3InDegrees));

        // Now try to insert the allegedly new BondTorsion specification into the
        // amberImproperTorsion map.  If it is already there the 2nd element in the
        // returned pair will be 'false'.
    std::pair<std::map<AtomClassIdQuad,BondTorsion>::iterator, bool> ret =
      mm.amberImproperTorsion.insert(std::pair<AtomClassIdQuad,BondTorsion>(key,bt));

        // Throw an exception if terms for this improper torsion were already defined.
    SimTK_APIARGCHECK4_ALWAYS(ret.second, mm.ApiClassName, MethodName,
        "amber improper torsion term(s) were already defined for atom class quad (%d,%d,%d,%d)",
        (int)key[0], (int)key[1], (int)key[2], (int)key[3]);
}

// Convenient signature for an amber improper torsion with only one term.
void DuMMForceFieldSubsystem::defineAmberImproperTorsion
   (DuMM::AtomClassId class1, DuMM::AtomClassId class2, DuMM::AtomClassId class3, DuMM::AtomClassId class4,
    int periodicity1, Real amp1InKJ, Real phase1InDegrees)
{
    defineAmberImproperTorsion(class1, class2, class3, class4,
                               periodicity1,amp1InKJ,phase1InDegrees,
                               -1,0.,0., -1,0.,0.);
}

// Convenient signature for an amber improper torsion with two terms.
void DuMMForceFieldSubsystem::defineAmberImproperTorsion
   (DuMM::AtomClassId class1, DuMM::AtomClassId class2, DuMM::AtomClassId class3, DuMM::AtomClassId class4,
    int periodicity1, Real amp1InKJ, Real phase1InDegrees,
    int periodicity2, Real amp2InKJ, Real phase2InDegrees)
{
    defineAmberImproperTorsion(class1, class2, class3, class4,
                               periodicity1,amp1InKJ,phase1InDegrees,
                               periodicity2,amp2InKJ,phase2InDegrees,
                               -1,0.,0.);
}

void DuMMForceFieldSubsystem::setVdwMixingRule(VdwMixingRule rule) {
    static const char* MethodName = "setVdwMixingRule";

    invalidateSubsystemTopologyCache();

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

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac && fac <= 1, mm.ApiClassName, MethodName,
        "van der Waals energy scale factor (%g) for 1-2 bonded atoms was invalid: must be between 0 and 1, inclusive",
        fac);

    mm.vdwScale12=fac;
}
void DuMMForceFieldSubsystem::setVdw13ScaleFactor(Real fac) {
    static const char* MethodName = "setVdw13ScaleFactor";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac && fac <= 1, mm.ApiClassName, MethodName,
        "van der Waals energy scale factor (%g) for 1-3 bonded atoms was invalid: must be between 0 and 1, inclusive",
        fac);

    mm.vdwScale13=fac;
}
void DuMMForceFieldSubsystem::setVdw14ScaleFactor(Real fac) {
    static const char* MethodName = "setVdw14ScaleFactor";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac && fac <= 1, mm.ApiClassName, MethodName,
        "van der Waals energy scale factor (%g) for 1-4 bonded atoms was invalid: must be between 0 and 1, inclusive",
        fac);

    mm.vdwScale14=fac;
}
void DuMMForceFieldSubsystem::setVdw15ScaleFactor(Real fac) {
    static const char* MethodName = "setVdw15ScaleFactor";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac && fac <= 1, mm.ApiClassName, MethodName,
        "van der Waals energy scale factor (%g) for 1-5 bonded atoms was invalid: must be between 0 and 1, inclusive",
        fac);

    mm.vdwScale15=fac;
}

void DuMMForceFieldSubsystem::setCoulomb12ScaleFactor(Real fac) {
    static const char* MethodName = "setCoulomb12ScaleFactor";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac && fac <= 1, mm.ApiClassName, MethodName,
        "Coulomb scale factor (%g) for 1-2 bonded atoms was invalid: must be between 0 and 1, inclusive",
        fac);

    mm.coulombScale12=fac;
}

void DuMMForceFieldSubsystem::setCoulomb13ScaleFactor(Real fac) {
    static const char* MethodName = "setCoulomb13ScaleFactor";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac && fac <= 1, mm.ApiClassName, MethodName,
        "Coulomb scale factor (%g) for 1-3 bonded atoms was invalid: must be between 0 and 1, inclusive",
        fac);

    mm.coulombScale13=fac;
}
void DuMMForceFieldSubsystem::setCoulomb14ScaleFactor(Real fac) {
    static const char* MethodName = "setCoulomb14ScaleFactor";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac && fac <= 1, mm.ApiClassName, MethodName,
        "Coulomb scale factor (%g) for 1-4 bonded atoms was invalid: must be between 0 and 1, inclusive",
        fac);

    mm.coulombScale14=fac;
}
void DuMMForceFieldSubsystem::setCoulomb15ScaleFactor(Real fac) {
    static const char* MethodName = "setCoulomb15ScaleFactor";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac && fac <= 1, mm.ApiClassName, MethodName,
        "Coulomb scale factor (%g) for 1-5 bonded atoms was invalid: must be between 0 and 1, inclusive",
        fac);

    mm.coulombScale15=fac;
}

void DuMMForceFieldSubsystem::setVdwGlobalScaleFactor(Real fac) {
    static const char* MethodName = "setVdwScaleFactor";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac, mm.ApiClassName, MethodName,
        "Global van der Waals scale factor (%g) was invalid: must be nonnegative",
        fac);

    mm.vdwGlobalScaleFactor=fac;
}

void DuMMForceFieldSubsystem::setCoulombGlobalScaleFactor(Real fac) {
    static const char* MethodName = "setCoulombScaleFactor";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac, mm.ApiClassName, MethodName,
        "Global Coulomb scale factor (%g) was invalid: must be nonnegative",
        fac);

    mm.coulombGlobalScaleFactor=fac;
}
void DuMMForceFieldSubsystem::setBondStretchGlobalScaleFactor(Real fac) {
    static const char* MethodName = "setBondStretchScaleFactor";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac, mm.ApiClassName, MethodName,
        "Global bond stretch scale factor (%g) was invalid: must be nonnegative",
        fac);

    mm.bondStretchGlobalScaleFactor=fac;
}
void DuMMForceFieldSubsystem::setBondBendGlobalScaleFactor(Real fac) {
    static const char* MethodName = "setBondBendScaleFactor";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac, mm.ApiClassName, MethodName,
        "Global bond bend scale factor (%g) was invalid: must be nonnegative",
        fac);

    mm.bondBendGlobalScaleFactor=fac;
}
void DuMMForceFieldSubsystem::setBondTorsionGlobalScaleFactor(Real fac) {
    static const char* MethodName = "setBondTorsionScaleFactor";
 
    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac, mm.ApiClassName, MethodName,
        "Global bond torsion scale factor (%g) was invalid: must be nonnegative",
        fac);

    mm.bondTorsionGlobalScaleFactor=fac;
}
void DuMMForceFieldSubsystem::setAmberImproperTorsionGlobalScaleFactor(Real fac) {
    static const char* MethodName = "setAmberImproperTorsionScaleFactor";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac, mm.ApiClassName, MethodName,
        "Global amber improper torsion scale factor (%g) was invalid: must be nonnegative",
        fac);

    mm.amberImproperTorsionGlobalScaleFactor=fac;
}

DuMM::ClusterId DuMMForceFieldSubsystem::createCluster(const char* groupName)
{
    // Currently there is no error checking to do. We don't insist on unique group names.
    return updRep().addCluster(Cluster(groupName));
}

DuMM::AtomId DuMMForceFieldSubsystem::addAtom(DuMM::ChargedAtomTypeId chargedAtomTypeId)
{
    static const char* MethodName = "addAtom";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(mm.isValidChargedAtomType(chargedAtomTypeId), mm.ApiClassName, MethodName, 
        "charged atom type %d is not valid", (int) chargedAtomTypeId);

    const DuMM::AtomId atomId = (const DuMM::AtomId)mm.atoms.size();
    mm.atoms.push_back(Atom(chargedAtomTypeId, atomId));
    return atomId;
}

void DuMMForceFieldSubsystem::placeAtomInCluster(DuMM::AtomId atomId, DuMM::ClusterId clusterId, const Vec3& stationInNm)
{
    static const char* MethodName = "placeAtomInCluster";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

        // Make sure that we've seen both the atomId and clusterId before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtom(atomId), mm.ApiClassName, MethodName,
        "atom Id %d is not valid", (int) atomId);
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidCluster(clusterId), mm.ApiClassName, MethodName,
        "cluster Id %d is not valid", (int) clusterId);

    Cluster& cluster = mm.updCluster(clusterId);

        // Make sure that this cluster doesn't already contain this atom, either directly
        // or recursively through its subclusters.
    SimTK_APIARGCHECK3_ALWAYS(!cluster.containsAtom(atomId), mm.ApiClassName, MethodName,
        "cluster %d('%s') already contains atom %d", (int) clusterId, cluster.name.c_str(), (int) atomId);

        // Add the atom to the cluster.
    cluster.placeAtom(atomId, stationInNm, mm);
}

void DuMMForceFieldSubsystem::placeClusterInCluster
   (DuMM::ClusterId childClusterId, DuMM::ClusterId parentClusterId, const Transform& placementInNm)
{
    static const char* MethodName = "placeClusterInCluster";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

        // Make sure that we've seen both of these clusters before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidCluster(childClusterId), mm.ApiClassName, MethodName,
        "child cluster Id %d is not valid", (int) childClusterId);
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidCluster(parentClusterId), mm.ApiClassName, MethodName,
        "parent cluster Id %d is not valid", (int) parentClusterId);

    Cluster&       parent = mm.updCluster(parentClusterId);
    const Cluster& child  = mm.getCluster(childClusterId);

        // TODO: for now, make sure the parent is a top-level cluster, meaning that it does
        // not have any parent clusters (although it can be attached to a body). This restriction
        // should be relaxed but it is tricky to get all the parents' and ancestors' content
        // lists updated correctly so I'm deferring that for now (sherm 060928).
    SimTK_APIARGCHECK2_ALWAYS(parent.isTopLevelCluster(), mm.ApiClassName, MethodName,
        "parent cluster %d('%s') is not a top-level cluster so you cannot add a child cluster to it now",
        (int) parentClusterId, parent.name.c_str());

        // Child must not already be attached to a body.
    SimTK_APIARGCHECK2_ALWAYS(!child.isAttachedToBody(), mm.ApiClassName, MethodName,
        "child cluster %d('%s') is already attached to a body so cannot now be placed in another cluster",
        (int) childClusterId, child.name.c_str());

        // Make sure that parent cluster doesn't already contain child cluster, either directly
        // or recursively through its subclusters.
    SimTK_APIARGCHECK4_ALWAYS(!parent.containsCluster(childClusterId), mm.ApiClassName, MethodName,
        "parent cluster %d('%s') already contains child cluster %d('%s')", 
        (int) parentClusterId, parent.name.c_str(), (int) childClusterId, child.name.c_str());

        // Make sure the new child cluster doesn't contain any atoms which are already in
        // any of the trees to which the parent cluster is associated.
        // TODO: for now we need only look at the parent since we know it is top level.
    DuMM::AtomId atomId;
    SimTK_APIARGCHECK5_ALWAYS(!parent.overlapsWithCluster(child, atomId), mm.ApiClassName, MethodName,
        "parent cluster %d('%s') and would-be child cluster %d('%s') both contain atom %d"
        " so they cannot have a parent/child relationship",
        (int) parentClusterId, parent.name.c_str(), (int) childClusterId, child.name.c_str(), (int) atomId);

        // Add the child cluster to the parent.
    parent.placeCluster(childClusterId, placementInNm, mm);
}

void DuMMForceFieldSubsystem::attachClusterToBody(DuMM::ClusterId clusterId, MobilizedBodyId bodyNum, 
                                                  const Transform& placementInNm) 
{
    static const char* MethodName = "attachClusterToBody";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

        // Make sure we've seen this cluster before, and that the body number is well formed.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidCluster(clusterId), mm.ApiClassName, MethodName,
        "cluster Id %d is not valid", (int) clusterId);
    SimTK_APIARGCHECK1_ALWAYS(bodyNum >= 0, mm.ApiClassName, MethodName,
        "body number %d is not valid: must be nonnegative", (int)bodyNum);

    const Cluster& child  = mm.getCluster(clusterId);

        // Child must not already be attached to a body.
    SimTK_APIARGCHECK3_ALWAYS(!child.isAttachedToBody(), mm.ApiClassName, MethodName,
        "cluster %d('%s') is already attached to body %d so cannot now be attached to a body",
        (int) clusterId, child.name.c_str(),  (int)child.getBodyId());

        // None of the atoms in the child can be attached to any body.
    DuMM::AtomId    atomId;
    MobilizedBodyId bodyId;
    SimTK_APIARGCHECK4_ALWAYS(!child.containsAnyAtomsAttachedToABody(atomId,bodyId,mm), 
        mm.ApiClassName, MethodName,
        "cluster %d('%s') contains atom %d which is already attached to body %d"
        " so the cluster cannot now be attached to another body",
        (int) clusterId, child.name.c_str(), (int) atomId, (int)bodyId);

        // Create an entry for the body if necessary, and its corresponding cluster.
    mm.ensureBodyEntryExists(bodyNum);
    Cluster& bodyCluster = mm.updCluster(mm.getBody(bodyNum).getClusterId());

        // Make sure that body cluster doesn't already contain child cluster, either directly
        // or recursively through its subclusters.
    SimTK_APIARGCHECK3_ALWAYS(!bodyCluster.containsCluster(clusterId), mm.ApiClassName, MethodName,
        "cluster %d('%s') is already attached (directly or indirectly) to body %d", 
        (int) clusterId, child.name.c_str(), (int)bodyNum);

        // OK, attach the cluster to the body's cluster.
    bodyCluster.placeCluster(clusterId, placementInNm, mm);
}

void DuMMForceFieldSubsystem::attachAtomToBody(DuMM::AtomId atomId, MobilizedBodyId bodyNum, const Vec3& stationInNm) 
{
    static const char* MethodName = "attachAtomToBody";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

        // Make sure we've seen this atom before, and that the body number is well formed.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtom(atomId), mm.ApiClassName, MethodName,
        "atom Id %d is not valid", (int) atomId);
    SimTK_APIARGCHECK1_ALWAYS(bodyNum >= 0, mm.ApiClassName, MethodName,
        "body number %d is not valid: must be nonnegative", (int)bodyNum);

        // The atom must not already be attached to a body, even this one.
    SimTK_APIARGCHECK2_ALWAYS(!mm.getAtom(atomId).isAttachedToBody(), mm.ApiClassName, MethodName,
        "atom %d is already attached to body %d so cannot now be attached to a body",
        (int) atomId, (int)mm.getAtom(atomId).getBodyId());

        // Create an entry for the body if necessary, and its corresponding cluster.
    mm.ensureBodyEntryExists(bodyNum);
    Cluster& bodyCluster = mm.updCluster(mm.getBody(bodyNum).getClusterId());

        // Attach the atom to the body's cluster.
    bodyCluster.placeAtom(atomId, stationInNm, mm);
}

MassProperties DuMMForceFieldSubsystem::calcClusterMassProperties
   (DuMM::ClusterId clusterId, const Transform& placementInNm) const
{
    static const char* MethodName = "calcClusterMassProperties";
    const DuMMForceFieldSubsystemRep& mm = getRep();

        // Make sure we've seen this cluster before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidCluster(clusterId), mm.ApiClassName, MethodName,
        "cluster Id %d is not valid", (int) clusterId);

    return mm.getCluster(clusterId).calcMassProperties(placementInNm, mm);
}


DuMM::BondId DuMMForceFieldSubsystem::addBond(DuMM::AtomId atom1Id, DuMM::AtomId atom2Id)
{
    static const char* MethodName = "addBond";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

        // Make sure we've seen these atoms before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtom(atom1Id), mm.ApiClassName, MethodName,
        "atom1(%d) is not valid", (int) atom1Id);
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtom(atom2Id), mm.ApiClassName, MethodName,
        "atom2(%d) is not valid", (int) atom2Id);

        // An atom can't be bonded to itself.
    SimTK_APIARGCHECK1_ALWAYS(atom1Id != atom2Id, mm.ApiClassName, MethodName,
        "the same atom Id(%d) was given for both atoms, which makes no sense", (int) atom1Id);

    // Ensure that atom1 < atom2
    if (atom1Id > atom2Id)
        std::swap(atom1Id,atom2Id);

    Atom& a1 = mm.updAtom(atom1Id);
    Atom& a2 = mm.updAtom(atom2Id);

    SimTK_APIARGCHECK2_ALWAYS(!a1.isBondedTo(atom2Id), mm.ApiClassName, MethodName,
        "atom %d is already bonded to atom %d; you can only do that once",
        (int) atom1Id, (int) atom2Id);

    mm.bonds.push_back(Bond(atom1Id,atom2Id));
    a1.bond12.push_back(atom2Id);
    a2.bond12.push_back(atom1Id);
    return (DuMM::BondId)(mm.bonds.size() - 1);
}

int DuMMForceFieldSubsystem::getNAtoms() const {
    return getRep().getNAtoms();
}
int DuMMForceFieldSubsystem::getNBonds() const {
    return getRep().getNBonds();
}

// 'which' is 0 or 1 to pick one of the two atoms whose Id we return.
DuMM::AtomId DuMMForceFieldSubsystem::getBondAtom(DuMM::BondId bondId, int which) const {
    static const char* MethodName = "getBondAtom";
    const DuMMForceFieldSubsystemRep& mm = getRep();

        // Make sure we've seen this bond before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidBond(bondId), mm.ApiClassName, MethodName,
        "bond %d is not valid", (int) bondId);

    SimTK_APIARGCHECK1_ALWAYS(which==0 || which==1, mm.ApiClassName, MethodName,
        "'which' was %d but must be 0 or 1 to choose one of the two atoms", which);

    return mm.bonds[bondId].atoms[which];
}

// Returned mass is in daltons (g/mol).
Real DuMMForceFieldSubsystem::getAtomMass(DuMM::AtomId atomId) const {
    static const char* MethodName = "getAtomMass";
    const DuMMForceFieldSubsystemRep& mm = getRep();

        // Make sure we've seen this atom before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtom(atomId), mm.ApiClassName, MethodName,
        "atom %d is not valid", (int) atomId);

    const Element& e = mm.elements[mm.getAtomElementNum(atomId)];
    return e.mass;
}

// Returns the atomic number (number of protons in nucleus).
int DuMMForceFieldSubsystem::getAtomElement(DuMM::AtomId atomId) const {
    static const char* MethodName = "getAtomElement";
    const DuMMForceFieldSubsystemRep& mm = getRep();

        // Make sure we've seen this atom before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtom(atomId), mm.ApiClassName, MethodName,
        "atom %d is not valid", (int) atomId);

    return mm.getAtomElementNum(atomId);
}

Vec3 DuMMForceFieldSubsystem::getAtomDefaultColor(DuMM::AtomId atomId) const {
    static const char* MethodName = "getAtomDefaultColor";
    const DuMMForceFieldSubsystemRep& mm = getRep();

        // Make sure we've seen this atom before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtom(atomId), mm.ApiClassName, MethodName,
        "atom %d is not valid", (int) atomId);

    const Element& e = mm.elements[mm.getAtomElementNum(atomId)];
    return e.defaultColor;
}

// Returned radius is in nm.
Real DuMMForceFieldSubsystem::getAtomRadius(DuMM::AtomId atomId) const {
    static const char* MethodName = "getAtomRadius";
    const DuMMForceFieldSubsystemRep& mm = getRep();

        // Make sure we've seen this atom before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtom(atomId), mm.ApiClassName, MethodName,
        "atom %d is not valid", (int) atomId);

    const AtomClass& cl = mm.atomClasses[mm.getAtomClassId(atomId)];
    return cl.vdwRadius;
}

// Returned station is in nm.
Vec3 DuMMForceFieldSubsystem::getAtomStationOnBody(DuMM::AtomId atomId) const {
    static const char* MethodName = "getAtomStationOnBody";
    const DuMMForceFieldSubsystemRep& mm = getRep();

        // Make sure we've seen this atom before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtom(atomId), mm.ApiClassName, MethodName,
        "atom %d is not valid", (int) atomId);

    const Atom& a = mm.getAtom(atomId);

        // Atom must be attached to a body.
    SimTK_APIARGCHECK1_ALWAYS(a.isAttachedToBody(), mm.ApiClassName, MethodName,
        "atom %d is not attached to a body", (int) atomId);

    return a.station_B;
}

// Returned placement is in nm.
Transform DuMMForceFieldSubsystem::getClusterPlacementOnBody(DuMM::ClusterId clusterId) const {
    static const char* MethodName = "getClusterPlacementOnBody";
    const DuMMForceFieldSubsystemRep& mm = getRep();

        // Make sure we've seen this cluster before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidCluster(clusterId), mm.ApiClassName, MethodName,
        "cluster Id %d is not valid", (int) clusterId);

    const Cluster& c = mm.getCluster(clusterId);

        // Cluster must be attached to a body.
    SimTK_APIARGCHECK2_ALWAYS(c.isAttachedToBody(), mm.ApiClassName, MethodName,
        "cluster %d('%s') is not attached to a body", (int) clusterId, c.name.c_str());

    return c.placement_B;
}

// Returned station is in nm.
Vec3 DuMMForceFieldSubsystem::getAtomStationInCluster(DuMM::AtomId atomId, DuMM::ClusterId clusterId) const {
    static const char* MethodName = "getAtomStationInCluster";
    const DuMMForceFieldSubsystemRep& mm = getRep();

        // Make sure that we've seen both the atomId and clusterId before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtom(atomId), mm.ApiClassName, MethodName,
        "atom Id %d is not valid", (int) atomId);
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidCluster(clusterId), mm.ApiClassName, MethodName,
        "cluster Id %d is not valid", (int) clusterId);

    const Cluster& c = mm.getCluster(clusterId);
    const AtomPlacementSet& atoms = c.getAllContainedAtoms();
    const AtomPlacementSet::const_iterator ap = 
        atoms.find(AtomPlacement(atomId,Vec3(0)));

        // We're going to be upset of this cluster doesn't contain this atom.
    SimTK_APIARGCHECK3_ALWAYS(ap != atoms.end(), mm.ApiClassName, MethodName,
        "cluster %d('%s') does not contain atom %d", (int) clusterId, c.name.c_str(), (int) atomId);

    return ap->station;
}

// Returned placement is in nm.
Transform DuMMForceFieldSubsystem::getClusterPlacementInCluster(DuMM::ClusterId childClusterId, DuMM::ClusterId parentClusterId) const {
    static const char* MethodName = "getClusterPlacementInCluster";
    const DuMMForceFieldSubsystemRep& mm = getRep();

        // Make sure that we've seen both of these clusters before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidCluster(childClusterId), mm.ApiClassName, MethodName,
        "child cluster Id %d is not valid", (int) childClusterId);
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidCluster(parentClusterId), mm.ApiClassName, MethodName,
        "parent cluster Id %d is not valid", (int) parentClusterId);

    const Cluster& parent = mm.getCluster(parentClusterId);
    const Cluster& child  = mm.getCluster(childClusterId);

    const ClusterPlacementSet& clusters = parent.getAllContainedClusters();
    const ClusterPlacementSet::const_iterator cp = 
        clusters.find(ClusterPlacement(childClusterId,Transform()));

        // We're going to be upset of the parent cluster doesn't contain the child.
    SimTK_APIARGCHECK4_ALWAYS(cp != clusters.end(), mm.ApiClassName, MethodName,
        "cluster %d('%s') does not contain cluster %d('%d')", 
        (int) parentClusterId, parent.name.c_str(), (int) childClusterId, child.name.c_str());

    return cp->placement;
}

MobilizedBodyId DuMMForceFieldSubsystem::getAtomBody(DuMM::AtomId atomId) const {
    static const char* MethodName = "getAtomBody";
    const DuMMForceFieldSubsystemRep& mm = getRep();

        // Make sure that we've seen this atomId before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtom(atomId), mm.ApiClassName, MethodName,
        "atom Id %d is not valid", (int) atomId);

    const Atom& a = mm.getAtom(atomId);

        // Atom must be attached to a body.
    SimTK_APIARGCHECK1_ALWAYS(a.isAttachedToBody(), mm.ApiClassName, MethodName,
        "atom %d is not attached to a body", (int) atomId);

    return a.getBodyId();
}


MobilizedBodyId DuMMForceFieldSubsystem::getClusterBody(DuMM::ClusterId clusterId) const {
    static const char* MethodName = "getClusterBody";
    const DuMMForceFieldSubsystemRep& mm = getRep();

        // Make sure that we've seen this atomId before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidCluster(clusterId), mm.ApiClassName, MethodName,
        "cluster Id %d is not valid", (int) clusterId);

    const Cluster& c = mm.getCluster(clusterId);

        // Cluster must be attached to a body.
    SimTK_APIARGCHECK2_ALWAYS(c.isAttachedToBody(), mm.ApiClassName, MethodName,
        "cluster %d('%s') is not attached to a body", (int) clusterId, c.name.c_str());

    return c.getBodyId();
}

void DuMMForceFieldSubsystem::dump() const {
    return getRep().dump();
}


    ////////////////////////////////////
    // DUMM FORCE FIELD SUBSYSTEM REP //
    ////////////////////////////////////

/*static*/ const char* DuMMForceFieldSubsystemRep::ApiClassName 
    = "DuMMForceFieldSubsystem";

const BondStretch& 
DuMMForceFieldSubsystemRep::getBondStretch(DuMM::AtomClassId class1, DuMM::AtomClassId class2) const {
    static const BondStretch dummy; // invalid
    const AtomClassIdPair key(class1,class2,true);
    std::map<AtomClassIdPair,BondStretch>::const_iterator bs = bondStretch.find(key);
    return (bs != bondStretch.end()) ? bs->second : dummy;
}

const BondBend& 
DuMMForceFieldSubsystemRep::getBondBend(DuMM::AtomClassId class1, DuMM::AtomClassId class2, DuMM::AtomClassId class3) const {
    static const BondBend dummy; // invalid
    const AtomClassIdTriple key(class1, class2, class3, true);
    std::map<AtomClassIdTriple,BondBend>::const_iterator bb = bondBend.find(key);
    return (bb != bondBend.end()) ? bb->second : dummy;
}

const BondTorsion& 
DuMMForceFieldSubsystemRep::getBondTorsion
   (DuMM::AtomClassId class1, DuMM::AtomClassId class2, DuMM::AtomClassId class3, DuMM::AtomClassId class4) const
{
    static const AtomClassIdQuad dummyKey(DuMM::InvalidAtomClassId, DuMM::InvalidAtomClassId, DuMM::InvalidAtomClassId, DuMM::InvalidAtomClassId);
    static const BondTorsion dummy(dummyKey); // invalid

    const AtomClassIdQuad key(class1, class2, class3, class4, true);
    std::map<AtomClassIdQuad,BondTorsion>::const_iterator bt = bondTorsion.find(key);
    return (bt != bondTorsion.end()) ? bt->second : dummy;
}

const BondTorsion& 
DuMMForceFieldSubsystemRep::getAmberImproperTorsion
   (DuMM::AtomClassId class1, DuMM::AtomClassId class2, DuMM::AtomClassId class3, DuMM::AtomClassId class4) const
{
//xxx -> Randy's warning flag
    printf("aImp--classes: %d-%d-%d-%d\n", (int)class1,
                                  (int)class2,
                                  (int)class3,
                                  (int)class4);
    std::map<AtomClassIdQuad,BondTorsion>::const_iterator i;
    for (i=amberImproperTorsion.begin(); i!=amberImproperTorsion.end(); i++) {
        printf( "aImp-matches: %d-%d-%d-%d\n", (int)i->first[0],
                                      (int)(i->first[1]),
                                      (int)(i->first[2]),
                                      (int)(i->first[3]) );
    }
    
    static const AtomClassIdQuad dummyKey(DuMM::InvalidAtomClassId, DuMM::InvalidAtomClassId, DuMM::InvalidAtomClassId, DuMM::InvalidAtomClassId);
    static const BondTorsion dummy(dummyKey); // invalid

    const AtomClassIdQuad key(class1, class2, class3, class4, false);
    std::map<AtomClassIdQuad,BondTorsion>::const_iterator bt = amberImproperTorsion.find(key);
    return (bt != amberImproperTorsion.end()) ? bt->second : dummy;
}

int DuMMForceFieldSubsystemRep::realizeSubsystemTopologyImpl(State& s) const 
{

    // At realization time, we need to verify that every atom has a valid atom class id
    for (DuMM::AtomId anum = (DuMM::AtomId)0; (DuMM::AtomId)anum < (DuMM::AtomId)atoms.size(); ++anum) {
        if ( ! isValidChargedAtomType(atoms[anum].chargedAtomTypeId) ) {
            throw Exception::Base("Atom must have valid charged atom type before realizing topology");
        }
    }

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
    for (DuMM::AtomClassId i = (DuMM::AtomClassId)0; i < (DuMM::AtomClassId)atomClasses.size(); ++i) {
        if (!atomClasses[i].isValid()) continue;
        if (!atomClasses[i].isComplete()) continue;

        AtomClass& iclass = mutableThis->atomClasses[i];
        iclass.vdwDij.resize((int)atomClasses.size()-i, NaN);
        iclass.vdwEij.resize((int)atomClasses.size()-i, NaN); 
        for (DuMM::AtomClassId j=i; j < (DuMM::AtomClassId)atomClasses.size(); ++j) {
            const AtomClass& jclass = atomClasses[j];
            if (jclass.isValid() && jclass.isComplete())
                applyMixingRule(iclass.vdwRadius,    jclass.vdwRadius,
                                iclass.vdwWellDepth, jclass.vdwWellDepth,
                                iclass.vdwDij[j-i],  iclass.vdwEij[j-i]);

        }
    }

        // molecule

    // Process clusters & bodies (bodies are treated as top-level clusters)

    // We process clusters recursively, so we need to allow the clusters writable
    // access to the main DuMM object (i.e., "this").
    for (DuMM::ClusterId cnum = (DuMM::ClusterId)0; cnum < (DuMM::ClusterId)clusters.size(); ++cnum) {
        Cluster& c = mutableThis->clusters[cnum];
        assert(c.isValid()); // Shouldn't be any unused cluster numbers.
        c.realizeTopologicalCache(*mutableThis);
    }

    // Bodies, on the other hand, are always top level clusters and the
    // calculation here assumes that all the clusters have been processed.
    // Thus bodies need only read access to the main DuMM object, 
    // although we're passign the mutable one in so we can use the
    // same routine (TODO).
    for (MobilizedBodyId bnum(0); bnum < (int)bodies.size(); ++bnum) {
        DuMMBody& b = mutableThis->bodies[bnum];
        if (!b.isValid())
            continue; // OK for these to be unused.
        b.realizeTopologicalCache(*mutableThis);
    }

    // Assign body & station to every atom that has been assigned to a body.
    for (DuMM::AtomId anum = (DuMM::AtomId)0; anum < (DuMM::AtomId)atoms.size(); ++anum) {
        Atom& a = mutableThis->atoms[anum];
        a.bodyId = InvalidMobilizedBodyId;
    }
    for (MobilizedBodyId bnum(0); bnum < (int)bodies.size(); ++bnum) {
        const DuMMBody& b = bodies[bnum];
        if (!b.isValid())
            continue;   // Unused body numbers are OK.

        for (int i=0; i < (int)b.allAtoms.size(); ++i) {
            const AtomPlacement& ap = b.allAtoms[i]; assert(ap.isValid());
            Atom& a = mutableThis->atoms[ap.atomId]; assert(a.isValid());
            assert(a.bodyId == InvalidMobilizedBodyId); // Can only be on one body!!
            a.bodyId    = bnum;
            a.station_B = ap.station;
        }
    }
    for (DuMM::AtomId anum = (DuMM::AtomId)0; anum < (DuMM::AtomId)atoms.size(); ++anum) {
        const Atom& a = atoms[anum];
        assert(a.bodyId >= 0); // TODO catch unassigned atoms
    }

    // need to chase bonds to fill in the bonded data
    // Be sure to distinguish the *shortest* path between two atoms from 
    // the set of all paths between atoms.
    for (DuMM::AtomId anum = (DuMM::AtomId)0; (DuMM::AtomId)anum < (DuMM::AtomId)atoms.size(); ++anum) {
        Atom& a = mutableThis->atoms[anum];

        // This set is used to avoid duplicate paths in the shortestPath calculation.
        std::set<DuMM::AtomId> allBondedSoFar;

        // Only the bond12 list should be filled in at the moment. We'll sort
        // all the lists when they're done for good hygiene.
        std::sort(a.bond12.begin(), a.bond12.end());

        // Add this atom and its direct (1-2) bonds to the list of all bonded atoms.
        allBondedSoFar.insert(anum);
        allBondedSoFar.insert(a.bond12.begin(), a.bond12.end());

        // Find longer bond paths by building each list in turn from
        // the direct bonds of the atoms in the previous list.

        // build the bond13 and shortPath13 lists
        // - bond1x list gives *all* paths between bonded atoms where all the
        // atoms are distinct (i.e., no fair retracing one of the bonds or
        // running around a short loop to get back to the first atom again).
        // - shortPath1x list gives *shortest* path between bonded atoms
        a.bond13.clear();
        a.shortPath13.clear();
        for (int j=0; j < (int)a.bond12.size(); ++j) {
            const Atom& a12 = atoms[a.bond12[j]];
            const AtomArray& a12_12 = a12.bond12;
            for (int k=0; k < (int)a12_12.size(); ++k) {
                const DuMM::AtomId newAtom = a12_12[k];
                assert(newAtom != a.bond12[j]);
                if (newAtom == anum)
                    continue; // no loop backs!
                a.bond13.push_back(AtomIdPair(a.bond12[j], newAtom));

                // if no shorter path, note this short route
                if (allBondedSoFar.find(newAtom) == allBondedSoFar.end()) {
                    allBondedSoFar.insert(newAtom);
                    a.shortPath13.push_back(AtomIdPair(a.bond12[j], newAtom));
                }
            }
        }
        std::sort(a.bond13.begin(), a.bond13.end());
        std::sort(a.shortPath13.begin(), a.shortPath13.end());

        // Randy was too big of a sissy to combine the bond14 and shortPath14 computations!
        // Or, discretion is sometimes the better part of valor.

        // build the bond14 list (all non-overlapping, non-looped paths)
        a.bond14.clear();
        for (int j=0; j < (int)a.bond13.size(); ++j) {
            const Atom& a13 = atoms[a.bond13[j][1]];
            const AtomArray& a13_12 = a13.bond12;
            for (int k=0; k < (int)a13_12.size(); ++k) {
                const DuMM::AtomId newAtom = a13_12[k];
                assert(newAtom != a.bond13[j][1]);
                // avoid repeated atoms (loop back)
                if (newAtom!=anum && newAtom!=a.bond13[j][0]) {
                    a.bond14.push_back(AtomIdTriple(a.bond13[j][0],
                                                 a.bond13[j][1], newAtom));
                }
            }
        }
        std::sort(a.bond14.begin(), a.bond14.end());

        // build the shortPath14 list
        a.shortPath14.clear();
        for (int j=0; j < (int)a.shortPath13.size(); ++j) {
            const Atom& a13 = atoms[a.shortPath13[j][1]];
            const AtomArray& a13_12 = a13.bond12;
            for (int k=0; k < (int)a13_12.size(); ++k) {
                const DuMM::AtomId newAtom = a13_12[k];

                 // check if there was already a shorter path
                if (allBondedSoFar.find(newAtom) == allBondedSoFar.end()) {
                    allBondedSoFar.insert(newAtom);
                    a.shortPath14.push_back(AtomIdTriple(a.shortPath13[j][0],
                                                      a.shortPath13[j][1], newAtom));
                }
            }
        }
        std::sort(a.shortPath14.begin(), a.shortPath14.end());


        // build the bond15 list
        a.bond15.clear();
        for (int j=0; j < (int)a.bond14.size(); ++j) {
            const Atom& a14 = atoms[a.bond14[j][2]];
            const AtomArray& a14_12 = a14.bond12;
            for (int k=0; k < (int)a14_12.size(); ++k) {
                const DuMM::AtomId newAtom = a14_12[k];
                assert(newAtom != a.bond14[j][2]);

                // avoid repeats and loop back
                if (newAtom!=anum && newAtom!=a.bond14[j][0] && newAtom!=a.bond14[j][1]) {
                    a.bond15.push_back(AtomIdQuad(a.bond14[j][0],
                                               a.bond14[j][1],
                                               a.bond14[j][2], newAtom));
                }
            }
        }
        std::sort(a.bond15.begin(), a.bond15.end());

        // build the shortPath15 list
        a.shortPath15.clear();
        for (int j=0; j < (int)a.shortPath14.size(); ++j) {
            const Atom& a14 = atoms[a.shortPath14[j][2]];
            const AtomArray& a14_12 = a14.bond12;
            for (int k=0; k < (int)a14_12.size(); ++k) {
                const DuMM::AtomId newAtom = a14_12[k];

                // check if there was already a shorter path
                if (allBondedSoFar.find(newAtom) == allBondedSoFar.end()) {
                    allBondedSoFar.insert(newAtom);
                    a.shortPath15.push_back(AtomIdQuad(a.shortPath14[j][0],
                                                    a.shortPath14[j][1],
                                                    a.shortPath14[j][2], newAtom));
                }
            }
        }
        std::sort(a.shortPath15.begin(), a.shortPath15.end());

        // Find all atom that are connected to three (and only three) other atoms
        // Then add all orderings of ths to the improper torsion list
        a.bonds3Atoms.invalidate();
        if (3 == (int)a.bond12.size()) {
            a.bonds3Atoms = AtomIdTriple(a.bond12[0], a.bond12[1], a.bond12[2]);
        }

        // Fill in the cross-body bond lists. We only keep atoms which
        // are on a different body. We do this both for the all-bond lists
        // and the shortest bond lists.
        a.xbond12.clear();
        for (int j=0; j < (int)a.bond12.size(); ++j)
            if (atoms[a.bond12[j]].bodyId != a.bodyId)
                a.xbond12.push_back(a.bond12[j]);

        a.xbond13.clear(); a.xshortPath13.clear();
        for (int j=0; j < (int)a.bond13.size(); ++j)
            if (   atoms[a.bond13[j][0]].bodyId != a.bodyId
                || atoms[a.bond13[j][1]].bodyId != a.bodyId)
                a.xbond13.push_back(a.bond13[j]);
        for (int j=0; j < (int)a.shortPath13.size(); ++j)
            if (   atoms[a.shortPath13[j][0]].bodyId != a.bodyId
                || atoms[a.shortPath13[j][1]].bodyId != a.bodyId)
                a.xshortPath13.push_back(a.shortPath13[j]);

        a.xbond14.clear(); a.xshortPath14.clear();
        for (int j=0; j < (int)a.bond14.size(); ++j)
            if (   atoms[a.bond14[j][0]].bodyId != a.bodyId
                || atoms[a.bond14[j][1]].bodyId != a.bodyId
                || atoms[a.bond14[j][2]].bodyId != a.bodyId)
                a.xbond14.push_back(a.bond14[j]);
        for (int j=0; j < (int)a.shortPath14.size(); ++j)
            if (   atoms[a.shortPath14[j][0]].bodyId != a.bodyId
                || atoms[a.shortPath14[j][1]].bodyId != a.bodyId
                || atoms[a.shortPath14[j][2]].bodyId != a.bodyId)
                a.xshortPath14.push_back(a.shortPath14[j]);

        a.xbond15.clear(); a.xshortPath15.clear();
        for (int j=0; j < (int)a.bond15.size(); ++j)
            if (   atoms[a.bond15[j][0]].bodyId != a.bodyId
                || atoms[a.bond15[j][1]].bodyId != a.bodyId
                || atoms[a.bond15[j][2]].bodyId != a.bodyId
                || atoms[a.bond15[j][3]].bodyId != a.bodyId)
                a.xbond15.push_back(a.bond15[j]);
        for (int j=0; j < (int)a.shortPath15.size(); ++j)
            if (   atoms[a.shortPath15[j][0]].bodyId != a.bodyId
                || atoms[a.shortPath15[j][1]].bodyId != a.bodyId
                || atoms[a.shortPath15[j][2]].bodyId != a.bodyId
                || atoms[a.shortPath15[j][3]].bodyId != a.bodyId)
                a.xshortPath15.push_back(a.shortPath15[j]);

        a.xbonds3Atoms.invalidate();
        // If there were 3 bonds, and at least one of them is
        // on a different body, then we win!
        if (a.bonds3Atoms.isValid() && 
            (   atoms[a.bonds3Atoms[0]].bodyId != a.bodyId
             || atoms[a.bonds3Atoms[1]].bodyId != a.bodyId
             || atoms[a.bonds3Atoms[2]].bodyId != a.bodyId))
            a.xbonds3Atoms = a.bonds3Atoms;

        const DuMM::AtomClassId c1 = getAtomClassId(anum);

        // Save a BondStretch entry for each cross-body 1-2 bond
        a.stretch.resize(a.xbond12.size());
        for (int b12=0; b12 < (int)a.xbond12.size(); ++b12) {
            const DuMM::AtomClassId c2 = getAtomClassId(a.xbond12[b12]);
            a.stretch[b12] = getBondStretch(c1, c2);

            SimTK_REALIZECHECK2_ALWAYS(a.stretch[b12].isValid(),
                Stage::Topology, getMySubsystemId(), getName(),
                "couldn't find bond stretch parameters for cross-body atom class pair (%d,%d)", 
                (int) c1,(int) c2);
        }

        // Save a BondBend entry for each cross-body 1-3 bond
        a.bend.resize(a.xbond13.size());
        for (int b13=0; b13 < (int)a.xbond13.size(); ++b13) {
            const DuMM::AtomClassId c2 = getAtomClassId(a.xbond13[b13][0]);
            const DuMM::AtomClassId c3 = getAtomClassId(a.xbond13[b13][1]);
            a.bend[b13] = getBondBend(c1, c2, c3);

            SimTK_REALIZECHECK3_ALWAYS(a.bend[b13].isValid(),
                Stage::Topology, getMySubsystemId(), getName(),
                "couldn't find bond bend parameters for cross-body atom class triple (%d,%d,%d)", 
                (int) c1,(int) c2,(int) c3);
        }

        // Save a BondTorsion entry for each cross-body 1-4 bond
        a.torsion.resize(a.xbond14.size());
        for (int b14=0; b14 < (int)a.xbond14.size(); ++b14) {
            const DuMM::AtomClassId c2 = getAtomClassId(a.xbond14[b14][0]);
            const DuMM::AtomClassId c3 = getAtomClassId(a.xbond14[b14][1]);
            const DuMM::AtomClassId c4 = getAtomClassId(a.xbond14[b14][2]);
            a.torsion[b14] = getBondTorsion(c1, c2, c3, c4); 

            SimTK_REALIZECHECK4_ALWAYS(a.torsion[b14].isValid(),
                Stage::Topology, getMySubsystemId(), getName(),
                "couldn't find bond torsion parameters for cross-body atom class quad (%d,%d,%d,%d)", 
                (int) c1,(int) c2,(int) c3,(int) c4);
        }

        // Save *all* Amber improper torsion entries if this atom is bonded to three, and only 
        // three ohter atoms, *and* a matching amber improper torsion term is found in the
        // amberImproperTorsion array
        // Note that by convention, the center atom is in the third position
        // Also note that unlike AMBER keeps only *one* match, but we keep *all*.
        // To correct for this we also scale my the total number of matches.  This
        // is how TINKER implements AMBER's improper torsions.
        a.aImproperTorsion.clear();
        a.aImproperTorsion14.clear();
        if (a.xbonds3Atoms.isValid()) {
            for (int i2=0; i2<3; i2++) {
                for (int i3=0; i3<3; i3++) {
                    if (i3==i2) continue;
                    for (int i4=0; i4<3; i4++) {
                        if (i4==i2 || i4==i3) continue;
                        static const BondTorsion bt = getAmberImproperTorsion(
                                                 getAtomClassId(a.xbonds3Atoms[i2]),
                                                 getAtomClassId(a.xbonds3Atoms[i3]),
                                                 c1,
                                                 getAtomClassId(a.xbonds3Atoms[i4]));
                        if (bt.isValid()) {
                            printf("anum=%d: i2=%d i3=%d i4=%d\n", (int) anum, i2, i3, i4);
                            a.aImproperTorsion14.push_back(AtomIdTriple(
                                                           a.xbonds3Atoms[i2],
                                                           a.xbonds3Atoms[i3],
                                                           a.xbonds3Atoms[i4]));
                            a.aImproperTorsion.push_back(bt);
                        }
                    }
                }
            }
        }

    }

    return 0;
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

int DuMMForceFieldSubsystemRep::realizeSubsystemDynamicsImpl(const State& s) const 
{
    const MultibodySystem&        mbs    = getMultibodySystem(); // my owner
    const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();

    // Temps for scale factors; initialize to 1
    Vector vdwScale((int)atoms.size(), Real(1)); 
    Vector coulombScale((int)atoms.size(), Real(1));

    // Get access to system-global cache entries.
    Real& pe = mbs.updPotentialEnergy(s, Stage::Dynamics); // kJ
    Vector_<SpatialVec>& rigidBodyForces = 
        mbs.updRigidBodyForces(s, Stage::Dynamics); // kJ (torque), kJ/nm (force)

    for (MobilizedBodyId b1(0); b1 < (int)bodies.size(); ++b1) {
        const Transform&          X_GB1  = matter.getMobilizedBody(b1).getBodyTransform(s);
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
                const MobilizedBodyId b2 = a2.bodyId;
                assert(b2 != b1);
                const Transform& X_GB2 = matter.getMobilizedBody(a2.bodyId).getBodyTransform(s);
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
                const MobilizedBodyId b2 = a2.bodyId;
                const MobilizedBodyId b3 = a3.bodyId;
                assert(!(b2==b1 && b3==b1)); // shouldn't be on the list if all on 1 body

                // TODO: These might be the same body but for now we don't care.
                const Transform& X_GB2   = matter.getMobilizedBody(a2.bodyId).getBodyTransform(s);
                const Transform& X_GB3   = matter.getMobilizedBody(a3.bodyId).getBodyTransform(s);
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
                const Transform& X_GB2   = matter.getMobilizedBody(a2.bodyId).getBodyTransform(s);
                const Transform& X_GB3   = matter.getMobilizedBody(a3.bodyId).getBodyTransform(s);
                const Transform& X_GB4   = matter.getMobilizedBody(a4.bodyId).getBodyTransform(s);
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

            // Amber improper torsion
            // Note that a1 is the *third* atom in the torsion
            for (int b14=0; b14 < (int)a1.aImproperTorsion14.size(); ++b14) {
                const int a2num = a1.aImproperTorsion14[b14][0];
                const int a3num = a1.aImproperTorsion14[b14][1];
                const int a4num = a1.aImproperTorsion14[b14][2];
                assert(a4num != a1num);
                //if (a4num < a1num)
                //    continue; // don't process this bond this time

                const Atom& a2 = atoms[a2num];
                const Atom& a3 = atoms[a3num];
                const Atom& a4 = atoms[a4num];
                const int b2 = a2.bodyId;
                const int b3 = a3.bodyId;
                const int b4 = a4.bodyId;
                assert(!(b2==b1 && b3==b1 && b4==b1)); // shouldn't be on the list if all on 1 body

                // TODO: These might be the same body but for now we don't care.
                const Transform& X_GB2   = matter.getMobilizedBody(a2.bodyId).getBodyTransform(s);
                const Transform& X_GB3   = matter.getMobilizedBody(a3.bodyId).getBodyTransform(s);
                const Transform& X_GB4   = matter.getMobilizedBody(a4.bodyId).getBodyTransform(s);
                const Vec3       a2Station_G = X_GB2.R()*a2.station_B;
                const Vec3       a3Station_G = X_GB3.R()*a3.station_B;
                const Vec3       a4Station_G = X_GB4.R()*a4.station_B;
                const Vec3       a2Pos_G     = X_GB2.T() + a2Station_G;
                const Vec3       a3Pos_G     = X_GB3.T() + a3Station_G;
                const Vec3       a4Pos_G     = X_GB4.T() + a4Station_G;

                Real angle, energy;
                Vec3 f1, f2, f3, f4;
                const BondTorsion& bt = a1.aImproperTorsion[b14];
                //bt.periodic(a1Pos_G, a2Pos_G, a3Pos_G, a4Pos_G, bondTorsionGlobalScaleFactor,
                //            angle, energy, f1, f2, f3, f4);
                bt.periodic(a2Pos_G, a3Pos_G, a1Pos_G, a4Pos_G,
                            amberImproperTorsionGlobalScaleFactor,
                            angle, energy, f2, f3, f1, f4);

                pe += energy;
                rigidBodyForces[b1] += SpatialVec( a1Station_G % f1, f1);   // 15 flops
                rigidBodyForces[b2] += SpatialVec( a2Station_G % f2, f2);   // 15 flops
                rigidBodyForces[b3] += SpatialVec( a3Station_G % f3, f3);   // 15 flops
                rigidBodyForces[b4] += SpatialVec( a4Station_G % f4, f4);   // 15 flops
            }

/*
            if (a1.aImproperTorsion14.isValid()) {
                const Atom& a2 = atoms[a1.aImproperTorsion14[0]];
                const Atom& a3 = atoms[a1.aImproperTorsion14[1]];
                const Atom& a4 = atoms[a1.aImproperTorsion14[2]];
                const int b2 = a2.bodyId;
                const int b3 = a3.bodyId;
                const int b4 = a4.bodyId;
                assert(!(b2==b1 && b3==b1 && b4==b1)); // shouldn't be on the list if all on 1 body

                // TODO: These might be the same body but for now we don't care.
                const Transform& X_GB2   = matter.getMobilizedBody(a2.bodyId).getBodyTransform(s);
                const Transform& X_GB3   = matter.getMobilizedBody(a3.bodyId).getBodyTransform(s);
                const Transform& X_GB4   = matter.getMobilizedBody(a4.bodyId).getBodyTransform(s);
                const Vec3       a2Station_G = X_GB2.R()*a2.station_B;
                const Vec3       a3Station_G = X_GB3.R()*a3.station_B;
                const Vec3       a4Station_G = X_GB4.R()*a4.station_B;
                const Vec3       a2Pos_G     = X_GB2.T() + a2Station_G;
                const Vec3       a3Pos_G     = X_GB3.T() + a3Station_G;
                const Vec3       a4Pos_G     = X_GB4.T() + a4Station_G;

                Real angle, energy;
                Vec3 f1, f2, f3, f4;
                const BondTorsion& bt = a1.aImproperTorsion;
//                bt.periodic(a1Pos_G, a2Pos_G, a3Pos_G, a4Pos_G, amberImproperTorsionGlobalScaleFactor,
//                            angle, energy, f1, f2, f3, f4);
                bt.periodic(a2Pos_G, a3Pos_G, a1Pos_G, a4Pos_G, amberImproperTorsionGlobalScaleFactor,
                            angle, energy, f2, f3, f1, f4);
                pe += energy;
                rigidBodyForces[b1] += SpatialVec( a1Station_G % f1, f1);   // 15 flops
                rigidBodyForces[b2] += SpatialVec( a2Station_G % f2, f2);   // 15 flops
                rigidBodyForces[b3] += SpatialVec( a3Station_G % f3, f3);   // 15 flops
                rigidBodyForces[b4] += SpatialVec( a4Station_G % f4, f4);   // 15 flops
            }
*/

            scaleBondedAtoms(a1,vdwScale,coulombScale);
            for (MobilizedBodyId b2(b1+1); b2 < (int)bodies.size(); ++b2) {
                const Transform&          X_GB2  = matter.getMobilizedBody(b2).getBodyTransform(s);
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

    return 0;
}


// We scale short range interactions but only when the shortest bonded path
// cross bodies.
void DuMMForceFieldSubsystemRep::scaleBondedAtoms
   (const Atom& a, Vector& vdwScale, Vector& coulombScale) const 
{
    for (int i=0; i < (int)a.xbond12.size(); ++i) {
        const int ix = a.xbond12[i]; // those are also the shortest paths!
        vdwScale[ix]=vdwScale12; coulombScale[ix]=coulombScale12;
    }
    for (int i=0; i < (int)a.xshortPath13.size(); ++i) {
        const int ix = a.xshortPath13[i][1]; // the 2nd atom is the 1-3
        vdwScale[ix]=vdwScale13; coulombScale[ix]=coulombScale13;
    }
    if (vdwScale14 != 1 || coulombScale14 != 1)
        for (int i=0; i < (int)a.xshortPath14.size(); ++i) {
            const int ix = a.xshortPath14[i][2]; // the 3rd atom is the 1-4
            vdwScale[ix]=vdwScale14; coulombScale[ix]=coulombScale14;
        }
    if (vdwScale15 != 1 || coulombScale15 != 1)
        for (int i=0; i < (int)a.xshortPath15.size(); ++i) {
            const int ix = a.xshortPath15[i][3]; // the 4th atom is the 1-5
            vdwScale[ix]=vdwScale15; coulombScale[ix]=coulombScale15;
        }
}

void DuMMForceFieldSubsystemRep::unscaleBondedAtoms
   (const Atom& a, Vector& vdwScale, Vector& coulombScale) const 
{
    for (int i=0; i < (int)a.xbond12.size(); ++i) {
        const int ix = a.xbond12[i];    vdwScale[ix]=coulombScale[ix]=1;
    }
    for (int i=0; i < (int)a.xshortPath13.size(); ++i) {
        const int ix = a.xshortPath13[i][1]; vdwScale[ix]=coulombScale[ix]=1;
    }
    if (vdwScale14 != 1 || coulombScale14 != 1)
        for (int i=0; i < (int)a.xshortPath14.size(); ++i) {
            const int ix = a.xshortPath14[i][2]; vdwScale[ix]=coulombScale[ix]=1;
        }
    if (vdwScale15 != 1 || coulombScale15 != 1)
        for (int i=0; i < (int)a.xshortPath15.size(); ++i) {
            const int ix = a.xshortPath15[i][3]; vdwScale[ix]=coulombScale[ix]=1;
        }
}

// Element masses are given in daltons (==g/mol==amu==u).
void DuMMForceFieldSubsystemRep::loadElements() {
    elements.resize(111); // Room for 1-110

    elements[1] = Element( 1, "H", "hydrogen", 1.007947 ).setDefaultColor(Green);
    elements[2] = Element( 2, "He", "helium", 4.003 );
    elements[3] = Element( 3, "Li", "lithium", 6.941 );
    elements[4] = Element( 4, "Be", "beryllium", 9.012 );
    elements[5] = Element( 5, "B", "boron", 10.811 );
    elements[6] = Element( 6, "C", "carbon", 12.01078 ).setDefaultColor(Gray);
    elements[7] = Element( 7, "N", "nitrogen", 14.00672 ).setDefaultColor(Blue);
    elements[8] = Element( 8, "O", "oxygen", 15.99943 ).setDefaultColor(Red);
    elements[9] = Element( 9, "F", "fluorine", 18.998 );
    elements[10] = Element( 10, "Ne", "neon", 20.180 );
    elements[11] = Element( 11, "Na", "sodium", 22.989769282 );
    elements[12] = Element( 12, "Mg", "magnesium", 24.30506 );
    elements[13] = Element( 13, "Al", "aluminum", 26.982 );
    elements[14] = Element( 14, "Si", "silicon", 28.086 );
    elements[15] = Element( 15, "P", "phosphorus", 30.9737622 ).setDefaultColor(Magenta);
    elements[16] = Element( 16, "S", "sulfur", 32.0655 ).setDefaultColor(Yellow);
    elements[17] = Element( 17, "Cl", "chlorine", 35.4532 );
    elements[18] = Element( 18, "Ar", "argon", 39.948 );
    elements[19] = Element( 19, "K", "potassium", 39.09831 );
    elements[20] = Element( 20, "Ca", "calcium", 40.0784 );
    elements[21] = Element( 21, "Sc", "scandium", 44.956 );
    elements[22] = Element( 22, "Ti", "titanium", 47.88 );
    elements[23] = Element( 23, "V", "vanadium", 50.942 );
    elements[24] = Element( 24, "Cr", "chromium", 51.996 );
    elements[25] = Element( 25, "Mn", "manganese", 54.9380455 );
    elements[26] = Element( 26, "Fe", "iron", 55.8452 );
    elements[27] = Element( 27, "Co", "cobalt", 58.9331955 );
    elements[28] = Element( 28, "Ni", "nickel", 58.69342 );
    elements[29] = Element( 29, "Cu", "copper", 63.5463 );
    elements[30] = Element( 30, "Zn", "zinc", 65.4094 );
    elements[31] = Element( 31, "Ga", "gallium", 69.723 );
    elements[32] = Element( 32, "Ge", "germanium", 72.61 );
    elements[33] = Element( 33, "As", "arsenic", 74.922 );
    elements[34] = Element( 34, "Se", "selenium", 78.963 );
    elements[35] = Element( 35, "Br", "bromine", 79.9041 );
    elements[36] = Element( 36, "Kr", "krypton", 83.80 );
    elements[37] = Element( 37, "Rb", "rubidium", 85.468 );
    elements[38] = Element( 38, "Sr", "strontium", 87.62 );
    elements[39] = Element( 39, "Y", "yttrium", 88.906 );
    elements[40] = Element( 40, "Zr", "zirconium", 91.224 );
    elements[41] = Element( 41, "Nb", "niobium", 92.906 );
    elements[42] = Element( 42, "Mo", "molybdenum", 95.94 );
    elements[43] = Element( 43, "Tc", "technetium", 97.907 );
    elements[44] = Element( 44, "Ru", "ruthenium", 101.07 );
    elements[45] = Element( 45, "Rh", "rhodium", 102.906 );
    elements[46] = Element( 46, "Pd", "palladium", 106.42 );
    elements[47] = Element( 47, "Ag", "silver", 107.868 );
    elements[48] = Element( 48, "Cd", "cadmium", 112.411 );
    elements[49] = Element( 49, "In", "indium", 114.82 );
    elements[50] = Element( 50, "Sn", "tin", 118.710 );
    elements[51] = Element( 51, "Sb", "antimony", 121.757 );
    elements[52] = Element( 52, "Te", "tellurium", 127.60 );
    elements[53] = Element( 53, "I", "iodine", 126.904 );
    elements[54] = Element( 54, "Xe", "xenon", 131.290 );
    elements[55] = Element( 55, "Cs", "cesium", 132.905 );
    elements[56] = Element( 56, "Ba", "barium", 137.327 );
    elements[57] = Element( 57, "La", "lanthanum", 138.906 );
    elements[58] = Element( 58, "Ce", "cerium", 140.115 );
    elements[59] = Element( 59, "Pr", "praseodymium", 140.908 );
    elements[60] = Element( 60, "Nd", "neodymium", 144.24 );
    elements[61] = Element( 61, "Pm", "promethium", 144.913 );
    elements[62] = Element( 62, "Sm", "samarium", 150.36 );
    elements[63] = Element( 63, "Eu", "europium", 151.965 );
    elements[64] = Element( 64, "Gd", "gadolinium", 157.25 );
    elements[65] = Element( 65, "Tb", "terbium", 158.925 );
    elements[66] = Element( 66, "Dy", "dysprosium", 162.50 );
    elements[67] = Element( 67, "Ho", "holmium", 164.930 );
    elements[68] = Element( 68, "Er", "erbium", 167.26 );
    elements[69] = Element( 69, "Tm", "thulium", 168.934 );
    elements[70] = Element( 70, "Yb", "ytterbium", 173.04 );
    elements[71] = Element( 71, "Lu", "lutetium", 174.967 );
    elements[72] = Element( 72, "Hf", "hafnium", 178.49 );
    elements[73] = Element( 73, "Ta", "tantalum", 180.948 );
    elements[74] = Element( 74, "W", "tungsten", 183.84 );
    elements[75] = Element( 75, "Re", "rhenium", 186.207 );
    elements[76] = Element( 76, "Os", "osmium", 190.2 );
    elements[77] = Element( 77, "Ir", "iridium", 192.22 );
    elements[78] = Element( 78, "Pt", "platinum", 195.08 );
    elements[79] = Element( 79, "Au", "gold", 196.967 ).setDefaultColor(Yellow);
    elements[80] = Element( 80, "Hg", "mercury", 200.59 );
    elements[81] = Element( 81, "Tl", "thallium", 204.383 );
    elements[82] = Element( 82, "Pb", "lead", 207.2 );
    elements[83] = Element( 83, "Bi", "bismuth", 208.980 );
    elements[84] = Element( 84, "Po", "polonium", 208.982 );
    elements[85] = Element( 85, "At", "astatine", 209.978 );
    elements[86] = Element( 86, "Rn", "radon", 222.018 );
    elements[87] = Element( 87, "Fr", "francium", 223.020 );
    elements[88] = Element( 88, "Ra", "radium", 226.025 );
    elements[89] = Element( 89, "Ac", "actinium", 227.028 );
    elements[90] = Element( 90, "Th", "thorium", 232.038 );
    elements[91] = Element( 91, "Pa", "protactinium", 231.038 );
    elements[92] = Element( 92, "U", "uranium", 238.028913 );
    elements[93] = Element( 93, "Np", "neptunium", 237.048 );
    elements[94] = Element( 94, "Pu", "plutonium", 244.064 );
    elements[95] = Element( 95, "Am", "americium", 243.061 );
    elements[96] = Element( 96, "Cm", "curium", 247.070 );
    elements[97] = Element( 97, "Bk", "berkelium", 247.070 );
    elements[98] = Element( 98, "Cf", "californium", 251.080 );
    elements[99] = Element( 99, "Es", "einsteinium", 252.083 );
    elements[100] = Element( 100, "Fm", "fermium", 257.095 );
    elements[101] = Element( 101, "Md", "mendelevium", 258.099 );
    elements[102] = Element( 102, "No", "nobelium", 259.101 );
    elements[103] = Element( 103, "Lr", "lawrencium", 260.105 );
    elements[104] = Element( 104, "Rf", "rutherfordium", 261 );
    elements[105] = Element( 105, "Db", "dubnium", 262 );
    elements[106] = Element( 106, "Sg", "seaborgium", 263 );
    elements[107] = Element( 107, "Bh", "bohrium", 262 );
    elements[108] = Element( 108, "Hs", "hassium", 265 );
    elements[109] = Element( 109, "Mt", "meitnerium", 266 );
    elements[110] = Element( 110, "Ds", "darmstadtium", 281 );
}

void DuMMForceFieldSubsystemRep::dump() const 
{
    printf("Dump of DuMMForceFieldSubsystem:\n");
    printf("  NBodies=%d NClusters=%d NAtoms=%d NAtomClasses=%d NChargedAtomTypes=%d NBonds=%d\n",
        bodies.size(), clusters.size(), atoms.size(), 
        atomClasses.size(), chargedAtomTypes.size(), bonds.size());
    for (int i=0; i < (int)bodies.size(); ++i) {
        printf("  DuMMBody %d:\n", i);
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

    ///////////////
    // BOND BEND //
    ///////////////

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

    //////////////////
    // BOND TORSION //
    //////////////////

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
    // ATOM //
    //////////

void Atom::dump() const {
    printf(" chargedAtomType=%d body=%d station=%g %g %g\n",
            (int) chargedAtomTypeId, (int)bodyId, station_B[0], station_B[1], station_B[2]);

    printf("          bond 1-2:");
    for (int i=0; i < (int)bond12.size(); ++i)
        printf(" %d", (int) bond12[i]);
    printf("\n          bond 1-3:");
    for (int i=0; i < (int)bond13.size(); ++i)
        printf(" %d-%d", (int) bond13[i][0], (int) bond13[i][1]);
    printf("\n          bond 1-4:");
    for (int i=0; i < (int)bond14.size(); ++i)
        printf(" %d-%d-%d", (int) bond14[i][0], (int) bond14[i][1], (int) bond14[i][2]);
    printf("\n          bond 1-5:");
    for (int i=0; i < (int)bond15.size(); ++i)
        printf(" %d-%d-%d-%d", (int) bond15[i][0], (int) bond15[i][1], (int) bond15[i][2], (int) bond15[i][3]);
    printf("\n     shortPath 1-3:");
    for (int i=0; i < (int)shortPath13.size(); ++i)
        printf(" %d-%d", (int) shortPath13[i][0], (int) shortPath13[i][1]);
    printf("\n     shortPath 1-4:");
    for (int i=0; i < (int)shortPath14.size(); ++i)
        printf(" %d-%d-%d", (int) shortPath14[i][0], (int) shortPath14[i][1], (int) shortPath14[i][2]);
    printf("\n     shortPath 1-5:");
    for (int i=0; i < (int)shortPath15.size(); ++i)
        printf(" %d-%d-%d-%d", (int) shortPath15[i][0], (int) shortPath15[i][1], (int) shortPath15[i][2], (int) shortPath15[i][3]);
    printf("\n       center of 3:");
    if (bonds3Atoms.isValid())
        printf(" %d-%d-%d", (int) bonds3Atoms[0], (int) bonds3Atoms[1], (int) bonds3Atoms[2]);
    printf("\n");

    printf("         xbond 1-2:");
    for (int i=0; i < (int)xbond12.size(); ++i)
        printf(" %d", (int) xbond12[i]);
    printf("\n         xbond 1-3:");
    for (int i=0; i < (int)xbond13.size(); ++i)
        printf(" %d-%d", (int) xbond13[i][0], (int) xbond13[i][1]);
    printf("\n         xbond 1-4:");
    for (int i=0; i < (int)xbond14.size(); ++i)
        printf(" %d-%d-%d", (int) xbond14[i][0], (int) xbond14[i][1], (int) xbond14[i][2]);
    printf("\n         xbond 1-5:");
    for (int i=0; i < (int)xbond15.size(); ++i)
        printf(" %d-%d-%d-%d", (int) xbond15[i][0], (int) xbond15[i][1], (int) xbond15[i][2], (int) xbond15[i][3]);
    printf("\n    xshortPath 1-3:");
    for (int i=0; i < (int)xshortPath13.size(); ++i)
        printf(" %d-%d", (int) xshortPath13[i][0], (int) xshortPath13[i][1]);
    printf("\n    xshortPath 1-4:");
    for (int i=0; i < (int)xshortPath14.size(); ++i)
        printf(" %d-%d-%d", (int) xshortPath14[i][0], (int) xshortPath14[i][1], (int) xshortPath14[i][2]);
    printf("\n    xshortPath 1-5:");
    for (int i=0; i < (int)xshortPath15.size(); ++i)
        printf(" %d-%d-%d-%d", (int) xshortPath15[i][0], (int) xshortPath15[i][1], (int) xshortPath15[i][2], (int) xshortPath15[i][3]);
    printf("\n      xcenter of 3:");
    if (xbonds3Atoms.isValid())
        printf(" %d-%d-%d", (int) xbonds3Atoms[0], (int) xbonds3Atoms[1], (int) xbonds3Atoms[2]);
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
    if (0<aImproperTorsion14.size()) {
        printf("\n    Amber improper torsion atoms:\n");
        for (int i=0; i < (int)aImproperTorsion14.size(); ++i) {
            const BondTorsion& bt = aImproperTorsion[i];
            printf("      %d-%d-x-%d:", (int) aImproperTorsion14[i][0],
                                        (int) aImproperTorsion14[i][1],
                                        (int) aImproperTorsion14[i][2]);
            for (int j=0; j<(int)bt.terms.size(); ++j) {
                const TorsionTerm& tt = bt.terms[j];
                printf(" (%d:%g,%g)", tt.periodicity, tt.amplitude, tt.theta0);
            }
            printf("\n");
        }
    }
    printf("\n");
}

    /////////////
    // CLUSTER //
    /////////////


void Cluster::attachToBody(MobilizedBodyId bnum, const Transform& X_BR, DuMMForceFieldSubsystemRep& mm) {
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
bool Cluster::containsAnyAtomsAttachedToABody(DuMM::AtomId& atomId, MobilizedBodyId& bodyId, 
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
    atomId = DuMM::InvalidAtomId;
    bodyId = InvalidMobilizedBodyId;
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
void Cluster::placeAtom(DuMM::AtomId atomId, const Vec3& station, DuMMForceFieldSubsystemRep& mm) {
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
void Cluster::placeCluster(DuMM::ClusterId childClusterId, const Transform& placement, DuMMForceFieldSubsystemRep& mm) {
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

    ///////////////
    // DUMM BODY //
    ///////////////

void DuMMBody::realizeTopologicalCache(const DuMMForceFieldSubsystemRep& mm) {
    allAtoms.clear();
    const Cluster& c = mm.getCluster(clusterId);
    AtomPlacementSet::const_iterator ap = c.getAllContainedAtoms().begin();
    while (ap != c.getAllContainedAtoms().end()) {
        allAtoms.push_back(*ap);
        ++ap;
    }
}

