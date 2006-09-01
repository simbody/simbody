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

#include "ForceSubsystemRep.h"

#include <string>
#include <vector>
#include <cmath>
#include <cassert>
#include <set>
#include <map>
#include <algorithm>

namespace SimTK {

class IntPair {
public:
    IntPair() {ints[0]=ints[1]=-1;}
    IntPair(int i1, int i2) {ints[0]=i1; ints[1]=i2;}
    int operator[](int i) const {assert(0<=i&&i<2); return ints[i];}
    bool isValid() const {return ints[0]>0 && ints[1]>0;}
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
    IntTriple(int i1, int i2, int i3) {ints[0]= i1; ints[1]=i2; ints[2]=i3;}
    int operator[](int i) const {assert(0<=i&&i<3); return ints[i];}
    bool isValid() const {return ints[0]>0 && ints[1]>0 && ints[2]>0;}
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
    IntQuad(int i1, int i2, int i3, int i4) {
        ints[0]= i1; ints[1]=i2; ints[2]=i3; ints[3]=i4;
    }
    int operator[](int i) const {assert(0<=i&&i<4); return ints[i];}
    bool isValid() const {return ints[0]>0 && ints[1]>0 && ints[2]>0 && ints[3]>0;}
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

class AtomType {
public:
    AtomType() : mass(-1), vdwRadius(-1), vdwWellDepth(-1), partialCharge(-1) { }
    AtomType(Real m, Real rad, Real well, Real chg)
      : mass(m), vdwRadius(rad), vdwWellDepth(well), partialCharge(chg) { }
    bool isValid() const {return mass >= 0;}

    void dump() const {
        printf("    mass=%g, vdwRad=%g, vdwDepth=%g, chg=%g\n",
            mass, vdwRadius, vdwWellDepth, partialCharge);
        printf("    vdwDij:");
        for (int i=0; i< (int)vdwDij.size(); ++i)
            printf(" %g", vdwDij[i]);
        printf("\n    vdwEij:");
        for (int i=0; i< (int)vdwEij.size(); ++i)
            printf(" %g", vdwEij[i]);
        printf("\n");
    }

    Real mass;
    Real vdwRadius;     // ri, length units
    Real vdwWellDepth;  // ei, energy units
    Real partialCharge; // qi, charge units

    // After all types have been defined, we can calculate vdw 
    // combining rules for dmin and well depth energy. We only fill
    // in entries for higher-numberd atom types, so to find the
    // entry for type t, index these arrays by t-t0 where t0 is the
    // type number of the present AtomType.
    // Note that different combining rules may be used but they
    // will always result in a pair of vdw parameters.
    std::vector<Real> vdwDij;
    std::vector<Real> vdwEij;
};

// This represents bond-stretch information for a pair of atom types.
// Use an IntPair as a key.
class BondStretch {
public:
    BondStretch() : k(-1), d0(-1) { }
    BondStretch(Real stiffness, Real length) : k(stiffness), d0(length) { 
        assert(isValid());
    }
    bool isValid() const {return k >= 0 && d0 >= 0; }
    Real k;
    Real d0; // distance at which force is 0
};

class BondAngle {
};

class BondTorsion {
};

typedef std::vector<int> AtomList;
class Atom {
public:
    Atom() : type(-1), bodyNum(-1), bodyAtomNum(-1) {
    }
    Atom(int t, int body, int bodyAtom, const Vec3& st)
      : type(t), bodyNum(body), bodyAtomNum(bodyAtom), station(st) {
    }
    bool isValid() const {return type>=0 && bodyNum>=0 && bodyAtomNum>=0;}

    bool isBondedTo(int anum) const {
        for (int i=0; i<(int)bond12.size(); ++i)
            if (bond12[i] == anum) return true;
        return false;
    }

    void dump() const {
        printf(" type=%d body=%d bodyAtomNum=%d station=%g %g %g\n",
            type, bodyNum, bodyAtomNum, station[0], station[1], station[2]);

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
        printf("\n");
    }

    int  type;
    int  bodyNum;
    int  bodyAtomNum;
    Vec3 station;

    // This is a set of lists which identify atoms nearby in the
    // molecules bond structure. The first list are the directly
    // bonded (1-2) atoms; the 2nd list has the 1-(2)-3 bonded atoms,
    // etc. The current Atom is always "1" so it isn't stored.

    std::vector<int>       bond12;
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
    std::vector<BondAngle>   bend;    // same length as   " 1-3 list
    std::vector<BondTorsion> torsion; // same length as   " 1-4 list
};


class Bond {
public:
    Bond() { }
    Bond(int atom1, int atom2) : atoms(atom1,atom2) { }

    IntPair atoms;
};

class Body {
public:
    Body() : bodyNum(-1) { }
    bool isValid() const {return bodyNum >= 0;}
    int nAtoms() const {return atoms.size();}
    int addAtom(int anum) {
        atoms.push_back(anum);
        return (int)atoms.size()-1;
    }

    void dump() const {
        printf("    bodyNum=%d\n      atoms:", bodyNum);
        for (int i=0; i < (int)atoms.size(); ++i)
            printf(" %d", atoms[i]);
        printf("\n");
    }

    int      bodyNum;
    AtomList atoms;
};

// Assume units:
//    Ref: http://physics.nist.gov/constants
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
//    1/(4*pi*e0) = 332.063711
//    speed of light c=2.99792458e8 m/s (exact)
//    Joules(N-m)/Kcal = 4184 (exact)
//
// Note: we have to use consistent force units, meaning
//   Da-A/ps^2

static const Real ForceUnitsPerKcal = 418.4; // convert energy to consistent units (Da-A^2/ps^2)
// This is 1/(4*pi*e0) in units which convert e^2/A to kcal/mole.
static const Real CoulombFac = 332.063711 * ForceUnitsPerKcal;

class DuMMForceFieldSubsystemRep : public ForceSubsystemRep {

    // Topological variables; set once.
    std::vector<Body>     bodies;
    std::vector<Atom>     atoms;
    std::vector<Bond>     bonds;

    // Force field description
    std::vector<AtomType>            types;
    std::map<IntPair,   BondStretch> bondStretch;
    std::map<IntTriple, BondAngle>   bondAngle;
    std::map<IntQuad,   BondTorsion> bondTorsion;

    // Scale factors for nonbonded forces when applied to
    // atoms which are near in the graph formed by the bonds.
    Real vdwScale12, coulombScale12;    // default 0,0
    Real vdwScale13, coulombScale13;    // default 0,0
    Real vdwScale14, coulombScale14;    // default 1,1
    Real vdwScale15, coulombScale15;    // default 1,1

    void ensureBodyExists(int bnum) {
        assert(bnum >= 0);
        if ((int)bodies.size() <= bnum)
            bodies.resize(bnum+1);
        if (!bodies[bnum].isValid()) {
            bodies[bnum].bodyNum = bnum;
            bodies[bnum].atoms.clear();
        }
    }

    bool isValidAtom(int atomNum) const {
        return 0 <= atomNum && atomNum < (int)atoms.size() && atoms[atomNum].isValid();
    }

    bool isValidAtomType(int typeNum) const {
        return 0 <= typeNum && typeNum < (int)types.size() && types[typeNum].isValid();
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

public:
    DuMMForceFieldSubsystemRep()
     : ForceSubsystemRep("DuMMForceFieldSubsystem", "0.0.1"), 
       built(false)
    {
        vdwScale12=coulombScale12=vdwScale13=coulombScale13=0;
        vdwScale14=coulombScale14=vdwScale15=coulombScale15=1;
    }

    int getNAtoms() const {return (int)atoms.size();}

    void addAtomType(int id, const AtomType& type) {
        assert(id >= 0);
        if (id >= (int)types.size())
            types.resize(id+1);
        assert(!types[id].isValid());
        types[id] = type;
    }

    void addBondStretch(int type1, int type2, const BondStretch& bs) {
        assert(isValidAtomType(type1) && isValidAtomType(type2));
        // Canonicalize the pair to have lowest type # first
        const IntPair key(std::min(type1,type2), std::max(type1,type2));
        std::pair<std::map<IntPair,BondStretch>::iterator, bool> ret = 
          bondStretch.insert(std::pair<IntPair,BondStretch>(key,bs));
        assert(ret.second); // must not have been there already
    }

    const BondStretch& getBondStretch(int type1, int type2) const {
        const IntPair key(std::min(type1,type2), std::max(type1,type2));
        std::map<IntPair,BondStretch>::const_iterator bs = bondStretch.find(key);
        assert(bs != bondStretch.end());
        return bs->second;
    }

    void setVdw12ScaleFactor(Real fac) {assert(0<=fac && fac<=1); vdwScale12=fac;}
    void setVdw13ScaleFactor(Real fac) {assert(0<=fac && fac<=1); vdwScale13=fac;}
    void setVdw14ScaleFactor(Real fac) {assert(0<=fac && fac<=1); vdwScale14=fac;}
    void setVdw15ScaleFactor(Real fac) {assert(0<=fac && fac<=1); vdwScale15=fac;}

    void setCoulomb12ScaleFactor(Real fac) {assert(0<=fac && fac<=1); coulombScale12=fac;}
    void setCoulomb13ScaleFactor(Real fac) {assert(0<=fac && fac<=1); coulombScale13=fac;}
    void setCoulomb14ScaleFactor(Real fac) {assert(0<=fac && fac<=1); coulombScale14=fac;}
    void setCoulomb15ScaleFactor(Real fac) {assert(0<=fac && fac<=1); coulombScale15=fac;}

    int addAtom(int body, int type, const Vec3& station)
    {
        ensureBodyExists(body); // create if necessary
        const int atomNum     = (int)atoms.size();
        const int bodyAtomNum = bodies[body].addAtom(atomNum);
        atoms.push_back(Atom(type,body,bodyAtomNum,station));
        return atomNum;    
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
    bool built;

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

void DuMMForceFieldSubsystem::defineAtomType
   (int id, Real mass, Real vdwRadius, Real vdwWellDepth, Real charge)
{
    updRep().addAtomType(id, AtomType(mass, vdwRadius, vdwWellDepth, charge));
}

void DuMMForceFieldSubsystem::defineBondStretch
   (int type1, int type2, Real stiffness, Real nominalLength)
{
    updRep().addBondStretch(type1, type2, BondStretch(stiffness,nominalLength));
}


int DuMMForceFieldSubsystem::addAtom(int body, int type, const Vec3& station)
{
    return updRep().addAtom(body, type, station);
}

int DuMMForceFieldSubsystem::addBond(int atom1, int atom2)
{
    return updRep().addBond(atom1, atom2);
}

int DuMMForceFieldSubsystem::getNAtoms() const {
    return getRep().getNAtoms();
}

void DuMMForceFieldSubsystem::dump() const {
    return getRep().dump();
}



    ////////////////////////////////
    // DuMMForceFieldSubsystemRep //
    ////////////////////////////////

void DuMMForceFieldSubsystemRep::realizeConstruction(State& s) const {
    // We need to write once onto the 'cache' portion of the object once
    // the topology is known.
    DuMMForceFieldSubsystemRep* mutableThis = 
        const_cast<DuMMForceFieldSubsystemRep*>(this);

    // Calculate effective van der Waals parameters for all 
    // pairs of atom types. We only fill in the diagonal
    // and upper triangle; that is, each type contains
    // parameters for like types and all types whose
    // (arbitrary) type number is higher.
    for (int i=0; i < (int)types.size(); ++i) {
        if (!types[i].isValid()) continue;

        AtomType& itype = mutableThis->types[i];
        itype.vdwDij.resize((int)types.size()-i, CNT<Real>::getNaN());
        itype.vdwEij.resize((int)types.size()-i, CNT<Real>::getNaN()); 
        for (int j=i; j < (int)types.size(); ++j) {
            const AtomType& jtype = types[j];
            if (jtype.isValid())
                applyMixingRule(itype.vdwRadius,    jtype.vdwRadius,
                                itype.vdwWellDepth, jtype.vdwWellDepth,
                                itype.vdwDij[j-i],  itype.vdwEij[j-i]);

        }
    }
    // need to chase bonds to fill in the bonded data
    // Be sure only to find the *shortest* path between two atoms
    for (int i=0; i < (int)atoms.size(); ++i) {
        Atom& a = mutableThis->atoms[i];
        std::set<int> allBondedSoFar;   // to avoid duplicate paths

        // Only the bond12 list should be filled in at the moment. We'll sort
        // all the lists when they're done for good hygiene.
        std::sort(a.bond12.begin(), a.bond12.end());

        // Add this atom and its direct (1-2) bonds to the list of all bonded atoms.
        allBondedSoFar.insert(i);
        for (int j=0; j < (int)a.bond12.size(); ++j)
            allBondedSoFar.insert(a.bond12[j]);

        // Find longer bond paths by building each list in turn from
        // the direct bonds of the atoms in the previous list.

        // build the bond13 list
        a.bond13.clear();
        for (int j=0; j < (int)a.bond12.size(); ++j) {
            const Atom& a12 = atoms[a.bond12[j]];
            const AtomList& a12_12 = a12.bond12;
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
            const AtomList& a13_12 = a13.bond12;
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
            const AtomList& a14_12 = a14.bond12;
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
            if (atoms[a.bond12[j]].bodyNum != a.bodyNum)
                a.xbond12.push_back(a.bond12[j]);

        for (int j=0; j < (int)a.bond13.size(); ++j)
            if (   atoms[a.bond13[j][0]].bodyNum != a.bodyNum
                || atoms[a.bond13[j][1]].bodyNum != a.bodyNum)
                a.xbond13.push_back(a.bond13[j]);

        for (int j=0; j < (int)a.bond14.size(); ++j)
            if (   atoms[a.bond14[j][0]].bodyNum != a.bodyNum
                || atoms[a.bond14[j][1]].bodyNum != a.bodyNum
                || atoms[a.bond14[j][2]].bodyNum != a.bodyNum)
                a.xbond14.push_back(a.bond14[j]);

        for (int j=0; j < (int)a.bond15.size(); ++j)
            if (   atoms[a.bond15[j][0]].bodyNum != a.bodyNum
                || atoms[a.bond15[j][1]].bodyNum != a.bodyNum
                || atoms[a.bond15[j][2]].bodyNum != a.bodyNum
                || atoms[a.bond15[j][3]].bodyNum != a.bodyNum)
                a.xbond15.push_back(a.bond15[j]);

        // Save a BondStretch entry for each 1-2 bond
        a.stretch.resize(a.xbond12.size());
        for (int b12=0; b12 < (int)a.xbond12.size(); ++b12)
            a.stretch[b12] = getBondStretch(a.type, atoms[a.xbond12[b12]].type);
    }

    mutableThis->built = true;
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
        const Transform&        X_GB1  = matter.getBodyConfiguration(s,b1);
        const std::vector<int>& alist1 = bodies[b1].atoms;

        for (int i=0; i < (int)alist1.size(); ++i) {
            const Atom&     a1 = atoms[alist1[i]];
            const AtomType& a1type = types[a1.type];
            const Vec3      a1Station_G = X_GB1.R()*a1.station;
            const Vec3      a1Pos_G     = X_GB1.T() + a1Station_G;
            const Real      q1Fac = CoulombFac*a1type.partialCharge;

            // Bonded. Note that each bond will appear twice so we only process
            // it the time when its 1st atom has a lower ID than its last.
            for (int b12=0; b12 < (int)a1.xbond12.size(); ++b12) {
                const int a2num = a1.xbond12[b12];
                assert(a2num != i);
                if (a2num > i)
                    continue; // don't process this bond this time

                const Atom& a2 = atoms[a2num];
                const int b2 = a2.bodyNum;
                assert(b2 != b1);
                const Transform& X_GB2   = matter.getBodyConfiguration(s, a2.bodyNum);
                const Vec3       a2Station_G = X_GB2.R()*a2.station;
                const Vec3       a2Pos_G     = X_GB2.T() + a2Station_G;
                const Vec3       r = a2Pos_G - a1Pos_G;
                const Real       d = r.norm();

                const BondStretch& bs = a1.stretch[b12];
                const Real         x  = d - bs.d0;

                const Real k = bs.k*ForceUnitsPerKcal;
                const Real eStretch = 0.5*k*x*x;
                const Real fStretch = -k*x; // sign is as would be applied to a2
                const Vec3 f2 = (fStretch/d) * r;
                pe += eStretch;
                rigidBodyForces[b2] += SpatialVec( a2Station_G % f2, f2);   // 15 flops
                rigidBodyForces[b1] -= SpatialVec( a1Station_G % f2, f2);   // 15 flops
            }

            scaleBondedAtoms(a1,vdwScale,coulombScale);
            for (int b2=b1+1; b2 < (int)bodies.size(); ++b2) {
                const Transform&  X_GB2 = matter.getBodyConfiguration(s,b2);
                const AtomList& alist2 = bodies[b2].atoms;

                for (int j=0; j < (int)alist2.size(); ++j) {
                    const Atom&     a2 = atoms[alist2[j]];
                    const AtomType& a2type = types[a2.type];
                    
                    const Vec3  a2Station_G = X_GB2.R()*a2.station; // 15 flops
                    const Vec3  a2Pos_G     = X_GB2.T() + a2Station_G;  // 3 flops
                    const Vec3  r = a2Pos_G - a1Pos_G; // from a1 to a2 (3 flops)
                    const Real  d2 = r.normSqr(); // 5 flops

                    // Check for cutoffs on d2?

                    const Real  ood2 = 1/d2; // approx 10 flops

                    // Coulomb. This unfortunately needs the separation distance which
                    // is expensive. But if scale, q1, or q2 are zero we can skip that.

                    Real eCoulomb = 0, fCoulomb = 0;
                    const Real qq = coulombScale[j]*q1Fac*a2type.partialCharge; // 2 flops

                    if (qq != 0.) {
                        const Real ood  = std::sqrt(ood2); // approx 30 flops
                        eCoulomb = qq * ood; //  scale*(1/(4*pi*e0)) *  q1*q2/d       (1 flop)  
                        fCoulomb = eCoulomb; // -scale*(1/(4*pi*e0)) * -q1*q2/d^2 * d (factor of 1/d^2 missing)
                    }

                    // van der Waals.

                    // Get precomputed mixed dmin and emin. Must ask the lower-numbered atom type.
                    const Real dij = (a1.type <= a2.type ? a1type.vdwDij[a2.type-a1.type]
                                                         : a2type.vdwDij[a1.type-a2.type]);
                    const Real eij = (a1.type <= a2.type ? a1type.vdwEij[a2.type-a1.type]
                                                         : a2type.vdwEij[a1.type-a2.type])
                                     * ForceUnitsPerKcal;

                    const Real ddij2  = dij*dij*ood2;        // (dmin_ij/d)^2 (2 flops)
                    const Real ddij6  = ddij2*ddij2*ddij2;   // 2 flops
                    const Real ddij12 = ddij6*ddij6;         // 1 flop

                    const Real eVdw =      eij * (ddij12 - 2*ddij6); // 3 flops
                    const Real fVdw = 12 * eij * (ddij12 - ddij6);   // factor of 1/d^2 missing (3 flops)
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

void DuMMForceFieldSubsystemRep::dump() const 
{
    printf("Dump of DuMMForceFieldSubsystem:\n");
    printf("  NBodies=%d NAtoms=%d NAtomTypes=%d NBonds=%d\n",
        bodies.size(), atoms.size(), types.size(), bonds.size());
    for (int i=0; i < (int)bodies.size(); ++i) {
        printf("  Body %d:\n", i);
        bodies[i].dump();
    }
    for (int i=0; i < (int)atoms.size(); ++i) {
        printf("  Atom %d: ", i);
        atoms[i].dump();
    }
    for (int i=0; i < (int)types.size(); ++i) {
        if (!types[i].isValid()) continue;
        printf("  AtomType %d:\n", i);
        types[i].dump();
    }
}

} // namespace SimTK

