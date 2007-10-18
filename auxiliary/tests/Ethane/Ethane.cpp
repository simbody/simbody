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
 * This is an outer block for simulating ethane in various ways with Simbody.
 * This is about testing Simbody, *not* studying ethane!
 */

#include "SimTKsimbody.h"

#include <string>
#include <iostream>
#include <exception>
#include <cmath>
using std::cout;
using std::endl;

using namespace SimTK;
using namespace DuMM; // for conversion constants

// mapping from atoms to Amber99 Charged Atom Types
enum {
    A99EthaneCarbon   = 13,    // closest I could find in Amber99; not really right
    A99EthaneHydrogen = 14,
    ShermDoubleBondedOxygen = 9999
};

class Molecule {
public:
    Molecule(MobilizedBodyId pId, const Transform& parentMobilizerFrame,
             const MolecularMechanicsSystem& mmSys)
        : parentId(pId), mobilizerFrameOnParent(parentMobilizerFrame), 
          mmSystem(mmSys)
    { }

    // Translate and rotation the molecule as a whole. Usually this just means move the base body,
    // but some molecules may not have a single base, so they can override this default.
    virtual void setMoleculeTransform(State& s, const Transform& pos) const
    {
        getMatter().getMobilizedBody(bodies[0]).setQToFitTransform(s, pos);
    }

    virtual void setMoleculeVelocity(State& s, const SpatialVec& vel) const
    {
        getMatter().getMobilizedBody(bodies[0]).setUToFitVelocity(s, vel);
    }

    // This routine must set the internal mobilities to their nominal values, both
    // for position and velocity states. This State must have already been realized
    // to at least Model stage.
    virtual void setDefaultInternalState(State& s) const = 0; 

    int getNAtoms()  const {return (int)atoms.size();}
    int getNBodies() const {return (int)bodies.size();}

    // return atomId of ith atom in Molecule
    DuMM::AtomId getAtom(int i) const {return atoms[i];}

    // return bodyNum of ith body; 0 is molecule's base body
    MobilizedBodyId getBodyId(int i) const {return bodies[i];}

    const SimbodyMatterSubsystem& getMatter() const {
        return mmSystem.getMatterSubsystem();
    }
    const DuMMForceFieldSubsystem& getDuMM() const {
        return mmSystem.getMolecularMechanicsForceSubsystem();
    }
protected:
    std::vector<DuMM::AtomId>       atoms;
    std::vector<MobilizedBodyId>    bodies;
    MobilizedBodyId                 parentId;
    Transform                       mobilizerFrameOnParent;
    const MolecularMechanicsSystem& mmSystem;
};


//
// Ribose sugar ring for showing "pucker" modes.
//
//             H5T--O5T   H51
//                    \  / 
//                H52--C5  
//                     |             
//                     |      O4
//                     |    /    \   O1--HO1
//                     |  /        \ |
//                 H4--C4            C1--H1
//                       \          /
//                    H3--C3------C2--H2
//                        |       |
//                   H3T--O3T     O2--HO2
//
//  Atom  Class  ChargedAtomType
//
// 0 O5T    22     1232 -0.6223
// 1 H5T    31     1233  0.4295
// 2 C5     1      1002  0.0558
// 3 H51    35     1003  0.0679
// 4 H52    35     1004  0.0679
//                      -------
//                      -0.0012
//
// 5 O4     23     1096 -0.3548
// 6 C4     1      1094  0.1065
// 7 H4     35     1095  0.1174
//                      -------
//                      -0.1309
//
// 8 O3T    22     1237 -0.6541
// 9 H3T    31     1238  0.4376
//10 C3     1      1010  0.2022
//11 H3     35     1011  0.0615
//                      -------
//                       0.0473
//
//12 O2     22     1237 -0.6541
//13 HO2    31     1238  0.4376
//14 C2     1      1010  0.2022
//15 H2     35     1011  0.0615
//                      -------
//                       0.0473
//
//16 O1     22     1237 -0.6541
//17 HO1    31     1238  0.4376
//18 C1     1      1010  0.2022
//19 H1     35     1011  0.0615
//                      -------
//                       0.0473
//
//                   ----------
//                       0.0098
//
// Charged atom type assignments are cobbled together from Amber99.
//
class Ribose : public Molecule {
public:
    Ribose(MobilizedBodyId pId, const Transform& parentMobilizerFrame,
           MolecularMechanicsSystem& mmSys) 
      : Molecule(pId, parentMobilizerFrame, mmSys)
    {
        DuMMForceFieldSubsystem& mm = mmSys.updMolecularMechanicsForceSubsystem();
        
        atoms.push_back(mm.addAtom((DuMM::ChargedAtomTypeId)1232)); // 0
        atoms.push_back(mm.addAtom((DuMM::ChargedAtomTypeId)1233)); // 1
        atoms.push_back(mm.addAtom((DuMM::ChargedAtomTypeId)1002)); // 2
        atoms.push_back(mm.addAtom((DuMM::ChargedAtomTypeId)1003)); // 3
        atoms.push_back(mm.addAtom((DuMM::ChargedAtomTypeId)1004)); // 4

        mm.addBond(atoms[0],atoms[1]); mm.addBond(atoms[2],atoms[3]);
        mm.addBond(atoms[2],atoms[4]); mm.addBond(atoms[0],atoms[2]);

        atoms.push_back(mm.addAtom((DuMM::ChargedAtomTypeId)1096)); // 5
        atoms.push_back(mm.addAtom((DuMM::ChargedAtomTypeId)1094)); // 6
        atoms.push_back(mm.addAtom((DuMM::ChargedAtomTypeId)1095)); // 7

        mm.addBond(atoms[5],atoms[6]); mm.addBond(atoms[6],atoms[7]);
        mm.addBond(atoms[2],atoms[6]);

        atoms.push_back(mm.addAtom((DuMM::ChargedAtomTypeId)1237)); // 8
        atoms.push_back(mm.addAtom((DuMM::ChargedAtomTypeId)1238)); // 9
        atoms.push_back(mm.addAtom((DuMM::ChargedAtomTypeId)1010)); //10
        atoms.push_back(mm.addAtom((DuMM::ChargedAtomTypeId)1011)); //11
        mm.addBond(atoms[8],atoms[9]); mm.addBond(atoms[10],atoms[11]);
        mm.addBond(atoms[8],atoms[10]); 
        mm.addBond(atoms[6],atoms[10]);

        atoms.push_back(mm.addAtom((DuMM::ChargedAtomTypeId)1237)); //12
        atoms.push_back(mm.addAtom((DuMM::ChargedAtomTypeId)1238)); //13
        atoms.push_back(mm.addAtom((DuMM::ChargedAtomTypeId)1010)); //14
        atoms.push_back(mm.addAtom((DuMM::ChargedAtomTypeId)1011)); //15
        mm.addBond(atoms[12],atoms[13]); mm.addBond(atoms[14],atoms[15]);
        mm.addBond(atoms[12],atoms[14]);
        mm.addBond(atoms[10],atoms[14]);

        atoms.push_back(mm.addAtom((DuMM::ChargedAtomTypeId)1237)); //16
        atoms.push_back(mm.addAtom((DuMM::ChargedAtomTypeId)1238)); //17
        atoms.push_back(mm.addAtom((DuMM::ChargedAtomTypeId)1010)); //18
        atoms.push_back(mm.addAtom((DuMM::ChargedAtomTypeId)1011)); //19
        mm.addBond(atoms[16],atoms[17]); mm.addBond(atoms[18],atoms[19]);
        mm.addBond(atoms[16],atoms[18]);
        mm.addBond(atoms[14],atoms[18]);
        mm.addBond(atoms[5],atoms[18]); // O4-C1

        // Define the clusers.
        //twoOxygens = mm.createCluster("two oxygens");
        //mm.placeAtomInCluster(getO(0), twoOxygens, Vec3(0));
        //mm.placeAtomInCluster(getO(1), twoOxygens, Vec3(0,0,getNominalOOBondLength()));
    }

    // Use body zero's frame as the "molecule frame" (call it "F" here). We'll move
    // that one to F' and record the relative motion. Then we move all the others
    // by the same relative transform.
    void setMoleculeTransform(State& s, const Transform& X_GFprime) const
    {
        mmSystem.realize(s,Stage::Position);
        const Transform X_GF = getMatter().getMobilizedBody(bodies[0]).getBodyTransform(s); // current
        const Transform X_FFprime = ~X_GF * X_GFprime; // relative transform

        for (int i=0; i < (int)bodies.size(); ++i) {
            mmSystem.realize(s,Stage::Position);
            const Transform& X_GB = getMatter().getMobilizedBody(bodies[i]).getBodyTransform(s);
            const Transform  X_FprimeBprime = ~X_GF * X_GB; // we want this to be the same as X_FB
            const Transform  X_FBprime = X_FFprime * X_FprimeBprime;
            const Transform  X_GBprime = X_GF*X_FBprime;
            // This only works since the mobilizers are all ground-attached Cartesian.
            // The Stage is reduced due to the change to the q's below.
            getMatter().getMobilizedBody(bodies[i]).setQToFitTransform(s, X_GBprime);
        }
    }

protected:
    //int twoOxygens; // cluster
};

class CartesianRibose : public Ribose {
public:
    CartesianRibose(MobilizedBodyId pId, MolecularMechanicsSystem& mmSys)
      : Ribose(pId,Transform(),mmSys)
    {
        SimbodyMatterSubsystem&  matter = 
            SimbodyMatterSubsystem::updDowncast(mmSys.updMatterSubsystem());
        DuMMForceFieldSubsystem& mm = mmSys.updMolecularMechanicsForceSubsystem();

        MobilizedBody& parent = matter.updMobilizedBody(parentId);

        //bodies.push_back(GroundId);
        //mm.attachAtomToBody(0,GroundId,Vec3(0));
        for (int i=0; i<20; ++i) {
            bodies.push_back(
                MobilizedBody::Cartesian(parent, Transform(),    // parent mobilizer frame
                                         Body::Rigid(MassProperties(mm.getAtomMass(atoms[i]),Vec3(0),Inertia(0))),
                                         Transform()));          // inboard mobilizer frame
            mm.attachAtomToBody(atoms[i],bodies.back(),Vec3(0));
        }
    }

    void setDefaultInternalState(State& s) const {
        //TODO: these are not right -- some of the pieces have the wrong chirality.
        Real q[]={
            0,0,0,
            0.0159955,0.0930929,0.0179526,
            0.0862463,-0.0788493,0.0804993,
            0.188797,-0.0643956,0.0460469,
            0.0793338,-0.048313,0.185042,
            -0.0734023,-0.252303,0.137272,
            0.0493894,-0.227588,0.0688814,
            0.0372042,-0.251271,-0.0368623,
            0.276019,-0.266582,0.172778-.1,
            0.32527,-0.335729,0.218463-.1,
            0.154191,-0.324668,0.12799-.1, //XXX
            0.175325,-0.401742,0.0536479-.1,
            0.121915,-0.521913,0.270613,
            0.0560667,-0.560277,0.329768,
            0.0800273,-0.388735,0.245242,
            0.0915603,-0.326613,0.334268,
            -0.156899,-0.38217,0.30165,
            -0.220948,-0.313277,0.281256,
            -0.0642466,-0.379147,0.196642,
            -0.0852995,-0.45703,0.122944};
        for (int i=0; i<20; ++i) {
            const MobilizedBody& b = getMatter().getMobilizedBody(bodies[i]);
            b.setQVector(s,Vector(3,q+3*i));
        }

        //for (int i=0; i<20; ++i) {
         //   if (bodies[i] != GroundId)
         //       getMatter().setMobilizerCoordsAsVec3(s,bodies[i],Vec3(i/10., i%5/10., i%2/2.));
       // }
    }


private:
};

// An oxygen (O2) molecule has just two atoms, with a maximum
// of 6 degrees of freedom. We provide the following models:
//   CartesianO2 -- 3 dofs referred to the Ground origin
//   InternalCartesianO2 -- the first oxygen is measured
//        from ground origin, 2nd w.r.t. first
//   InternalO2 -- same as Rigid, but adds a bond stretch
//                 mobility between the atoms.
//   RigidO2 -- a single body with both atoms attached
//              at a nominal bond length. The body then
//              has *5* dofs w.r.t. ground (can't rotate
//              around the line between the atoms).
// The first three models are equivalent; they just use
// different parameterization for the 6 dofs. The RigidO2
// model has one fewer degree of freedom. Specifically, the
// very high frequency O=O stretch term has been eliminated.

class OxygenMolecule : public Molecule {
public:
    OxygenMolecule(MobilizedBodyId pId, const Transform& parentMobilizerFrame,
                   MolecularMechanicsSystem& mmSys) 
      : Molecule(pId, parentMobilizerFrame, mmSys)
    {
        DuMMForceFieldSubsystem& mm = mmSys.updMolecularMechanicsForceSubsystem();
        
        // Create the atoms and bonds. Atom 0 is O0, atom 1 is O1. O0 will serve
        // as the base frame for the molecule.
        for (int i=0;i<2;++i) atoms.push_back(mm.addAtom((DuMM::ChargedAtomTypeId)ShermDoubleBondedOxygen));
        mm.addBond(getO(0),getO(1));

        // Define the clusers.
        twoOxygens = mm.createCluster("two oxygens");

        mm.placeAtomInCluster(getO(0), twoOxygens, Vec3(0));
        mm.placeAtomInCluster(getO(1), twoOxygens, Vec3(0,0,getNominalOOBondLength()));
    }

    // Get the atom number for each oxygen.
    DuMM::AtomId getO(int i) const {assert(i==0||i==1); return getAtom(i);}

    Real getNominalOOBondLength() const {
        return 1.21 * Ang2Nm;
    }
protected:
    DuMM::ClusterId twoOxygens; // cluster
};

class RigidO2 : public OxygenMolecule {
public:
    RigidO2(MobilizedBodyId pId, MolecularMechanicsSystem& mmSys)
      : OxygenMolecule(pId,Transform(),mmSys)
    {
        SimbodyMatterSubsystem&  matter = 
            SimbodyMatterSubsystem::updDowncast(mmSys.updMatterSubsystem());
        DuMMForceFieldSubsystem& mm     = mmSys.updMolecularMechanicsForceSubsystem();

        MobilizedBody& parent = matter.updMobilizedBody(parentId);

        // Align cluster reference frame with body's. (5 dofs!)
        // FreeLine prevents rotation about Z, so make sure the body has its
        // O=O axis arranged along Z (or rotate the frame here).
        /* This doesn't work: 
        bodies.push_back(
            matter.addRigidBody(
                MassProperties(0,Vec3(0),Inertia(0)),
                Transform(),            // inboard mobilizer frame
                parent, Transform(),    // parent mobilizer frame
                Mobilizer::Cartesian()));
        // y
        bodies.push_back(
            matter.addRigidBody(
                MassProperties(0,Vec3(0),Inertia(0)),
                Transform(Rotation::aboutX(-90*Deg2Rad)),            // inboard mobilizer frame
                bodies.back(), Transform(Rotation::aboutX(-90*Deg2Rad)),    // parent mobilizer frame
                Mobilizer::Pin()));
        // x
        MassProperties mprops = mm.calcClusterMassProperties(twoOxygens, Transform());
        cout << "Inertia:" << mprops.getInertia();
        cout << "inertia kludge:" << mprops.getInertia()+Inertia(0,0,.4);
        MassProperties mpropsKludge(mprops.getMass(), mprops.getMassCenter(), mprops.getInertia() + Inertia(0,0,.4));
        bodies.push_back(
            matter.addRigidBody(
                mpropsKludge,
                Transform(Rotation::aboutY(90*Deg2Rad)),            // inboard mobilizer frame
                bodies.back(), Transform(Rotation::aboutY(90*Deg2Rad)),    // parent mobilizer frame
                Mobilizer::Pin()));
        
        */
        MassProperties mprops = mm.calcClusterMassProperties(twoOxygens, Transform());
        MassProperties mpropsKludge(mprops.getMass(), mprops.getMassCenter(), 
                                    mprops.getInertia() + Inertia(0,0,.01));

        bodies.push_back(
            MobilizedBody::FreeLine(
                parent, Transform(0*Vec3(.1,.7,.19)),    // parent mobilizer frame
                Body::Rigid(mprops),
                Transform(0*Vec3(0,0,.3))));             // inboard mobilizer frame

        mm.attachClusterToBody(twoOxygens, bodies.back(), Transform()); 
    }

    void setDefaultInternalState(State& s) const { } // none
};

// ethane:
// atom 0 is carbon0
// atoms 2,3,4 are attached to carbon0
// atom 1 is carbon1
// atoms 5,6,7 are attached to carbon1
//
// pre-built rigid clusters:
//   the two carbons
//   methyl 1 (atom 0) and hydrogens 1,2,3
//   methyl 2 (atom 4) and hydrogens 5,6,7
// Any cluster or individual atom can be assigned to a body, provided
// the resulting set of assignments represents a partitioning of
// the atoms across the bodies.
//
// Ethane has 8 atoms, modeled as point masses. Consequently there
// can be a maximum of 24 dofs. Choosing one atom as a "base", the
// maximum number of internal coordinates is 21. If the base has
// 6 dofs instead of three, then the maximal internal set is 18.
// We provide a variety of models below as examples of what can
// be done with Simbody, not necessarily because these are good
// models!

class EthaneMolecule : public Molecule {
public:
    EthaneMolecule(MobilizedBodyId pId, const Transform& parentMobilizerFrame,
                   MolecularMechanicsSystem&);

    // find the atoms
    DuMM::AtomId getC(int i) const {assert(i==0||i==1); return getAtom(i);}
    DuMM::AtomId getH(int whichCarbon, int whichHydrogen) const {
        assert(0<=whichCarbon&&whichCarbon<=1);
        assert(0<=whichHydrogen&&whichHydrogen<=2);
        return getAtom(2+whichCarbon*3+whichHydrogen);
    }

    Real getNominalCCBondLength() const {
        return 1.53688 * Ang2Nm;
    }
    Real getNominalCHBondLength() const {
        return 1.09 * Ang2Nm;
    }
    Real getNominalHCCBondAngle() const {
        return 109.5 * Deg2Rad;
    }

protected:
    // Some pre-built atom clusters.
    DuMM::ClusterId twoCarbons;
    DuMM::ClusterId methyl[2];
};

class OneDofEthane : public EthaneMolecule {
public:
    OneDofEthane(bool allowStretch, MobilizedBodyId pId, MolecularMechanicsSystem&);

    void setDefaultInternalState(State& s) const {
        const MobilizedBody& b = getMatter().getMobilizedBody(getBodyId(1));
        const int ndof = b.getNumU(s);
        for (int i=0; i<ndof; ++i) {
            b.setOneQ(s, i, 0);
            b.setOneU(s, i, 0);
        }
    }

    // Set stretch around the nominal length.
    void setCCStretch(Real stretchInNm, State& s) const {
        const MobilizedBody& b1 = getMatter().getMobilizedBody(getBodyId(1));
        assert(b1.getNumU(s) == 2);    // must have been build with Cylinder mobilizer
        const MobilizedBodyId CBody = getDuMM().getAtomBody(getC(1));
        const MobilizedBody& b = getMatter().getMobilizedBody(CBody);
        b.setOneQ(s, 1, stretchInNm);
    }

    void setTorsionAngleDeg(Real angleInDeg, State& s) const {
        const MobilizedBodyId CBody = getDuMM().getAtomBody(getC(1));
        const MobilizedBody& b = getMatter().getMobilizedBody(CBody);
        b.setOneQ(s, 0, angleInDeg*Deg2Rad);
    }

    // Rate is rad/ps
    void setTorsionRate(Real rateInRadPerPs, State& s) const {
        const MobilizedBodyId CBody = getDuMM().getAtomBody(getC(1));
        const MobilizedBody& b = getMatter().getMobilizedBody(CBody);
        b.setOneU(s, 0, rateInRadPerPs);
    }
};

class RigidEthane : public EthaneMolecule {
public:
    RigidEthane(Real torsionAngleInDeg, MobilizedBodyId pId, MolecularMechanicsSystem&);
    void setDefaultInternalState(State& s) const { } // doesn't have any
};

// Here is a model of ethane with 14 internal coordinates: stretch & torsion
// between the carbons, and bend-stretch for the connection between the 
// hydrogens and their carbons. That is, for each hydrogen, we are
// modeling the CH stretch term, and the HCC bend term. However, we
// are *not* permitting HCH bending, so the hydrogens will always be
// found arranged exactly 120 degrees apart. 
class FloppyEthane : public EthaneMolecule {
public:
    FloppyEthane(MobilizedBodyId pId, MolecularMechanicsSystem&);

    void setDefaultInternalState(State& s) const;

    // Set stretch around the nominal length.
    void setCCStretch(Real stretchInNm, State& s) const {
        const MobilizedBodyId CBody = getDuMM().getAtomBody(getC(1));
        const MobilizedBody& b = getMatter().getMobilizedBody(CBody);
        b.setOneQ(s, 1, stretchInNm);
    }

    void setTorsionAngleDeg(Real angleInDeg, State& s) const {
        const MobilizedBodyId CBody = getDuMM().getAtomBody(getC(1));
        const MobilizedBody& b = getMatter().getMobilizedBody(CBody);
        b.setOneQ(s, 0, angleInDeg*Deg2Rad);
    }

    // Rate is rad/ps
    void setTorsionRate(Real rateInRadPerPs, State& s) const {
        const MobilizedBodyId CBody = getDuMM().getAtomBody(getC(1));
        const MobilizedBody& b = getMatter().getMobilizedBody(CBody);
        b.setOneU(s, 0, rateInRadPerPs);
    }
};

static const Transform BodyFrame;   // identity transform on any body

// How it actually looks now:
int main() {
try
  {
    MolecularMechanicsSystem mbs;

    SimbodyMatterSubsystem   matter(mbs);
    DuMMForceFieldSubsystem  mm(mbs);
    GeneralForceElements     forces(mbs);
    DecorationSubsystem      artwork(mbs);
    UniformGravitySubsystem  gravity(mbs, Vec3(0,0,0));

    Real accuracy = 1e-2;
    Real outputInterval = .01;
    Real simulationLength = 100;
    //Real outputInterval = .1;
    //Real simulationLength = 10;

    const Real torsControlGain = /*100000*/0;
    const Real desiredTorsAngle = /*Pi/3*/0;

    forces.addGlobalEnergyDrain(.01);


    // AMBER 99

    mm.setVdw14ScaleFactor(1/2.); // reduce energy by these factors
    mm.setCoulomb14ScaleFactor(1/1.2);

    mm.defineAtomClass_KA(1,  "Amber99 CT", 6, 4, 1.9080, 0.1094);
    mm.defineAtomClass_KA(2,  "Amber99 C",  6, 3, 1.9080, 0.0860);
    mm.defineAtomClass_KA(3,  "Amber99 CA", 6, 3, 1.9080, 0.0860);
    mm.defineAtomClass_KA(4,  "Amber99 CM", 6, 3, 1.9080, 0.0860);

    mm.defineAtomClass_KA(9,  "Amber99 CB", 6, 3, 1.9080, 0.0860);
    mm.defineAtomClass_KA(22, "Amber99 OH", 8, 2, 1.7210, 0.2104);
    mm.defineAtomClass_KA(23, "Amber99 OS", 8, 2, 1.6837, 0.1700);

    mm.defineAtomClass_KA(24, "Amber99 O",  8, 1, 1.6612, 0.2100);
    mm.defineAtomClass_KA(25, "Amber99 O2", 8, 1, 1.6612, 0.2100); 

    mm.defineAtomClass_KA(31, "Amber99 HO", 1, 1, 0.0001, 0.0000);
    //mm.defineAtomClass_KA(31, "Amber99 HO", 1, 1, 1., 0.1);//KLUDGE

    mm.defineAtomClass_KA(34, "Amber99 HC", 1, 1, 1.4870, 0.0157); 
    mm.defineAtomClass_KA(35, "Amber99 H1", 1, 1, 1.3870, 0.0157);

    mm.defineChargedAtomType_KA(13, "Amber99 Alanine CB", 1, -0.1825);
    mm.defineChargedAtomType_KA(14, "Amber99 Alanine HB", 34, 0.0603);

    mm.defineChargedAtomType_KA(1002, "Amber99 R-Adenosine C5'",   1,  0.0558);
    mm.defineChargedAtomType_KA(1003, "Amber99 R-Adenosine H5'1", 35,  0.0679);
    mm.defineChargedAtomType_KA(1004, "Amber99 R-Adenosine H5'2", 35,  0.0679);
    mm.defineChargedAtomType_KA(1006, "Amber99 R-Adenosine H4'",  35,  0.1174);
    mm.defineChargedAtomType_KA(1007, "Amber99 R-Adenosine O4'",  23, -0.3548);

    mm.defineChargedAtomType_KA(1010, "Amber99 R-Adenosine C3'",   1,  0.2022);
    mm.defineChargedAtomType_KA(1011, "Amber99 R-Adenosine H3'",  35,  0.0615);

    mm.defineChargedAtomType_KA(1094, "Amber99 R-Uracil C4'",   1,  0.1065);
    mm.defineChargedAtomType_KA(1095, "Amber99 R-Uracil H4'",  35,  0.1174);
    mm.defineChargedAtomType_KA(1096, "Amber99 R-Uracil O4'",  23, -0.3548);

    mm.defineChargedAtomType_KA(1101, "Amber99 R-Uracil C2'",   1,   0.0670);
    mm.defineChargedAtomType_KA(1102, "Amber99 R-Uracil H2'1",  35,  0.0972);
    mm.defineChargedAtomType_KA(1103, "Amber99 R-Uracil O2'",   22, -0.6139);
    mm.defineChargedAtomType_KA(1104, "Amber99 R-Uracil HO'2",  31,  0.4186);

    mm.defineChargedAtomType_KA(1232, "Amber99 R-5'-Hydroxyl O5'",  22, -0.6223);
    mm.defineChargedAtomType_KA(1233, "Amber99 R-5'-Hydroxyl H5T",  31,  0.4295);
    mm.defineChargedAtomType_KA(1237, "Amber99 R-5'-Hydroxyl O3'",  22, -0.6541);
    mm.defineChargedAtomType_KA(1238, "Amber99 R-5'-Hydroxyl H3T",  31,  0.4376);

    mm.defineBondStretch_KA( 1, 1, 310., 1.5260);
    mm.defineBondStretch_KA( 1,22, 320., 1.4100);
    mm.defineBondStretch_KA( 1,23, 320., 1.4100);
    mm.defineBondStretch_KA( 1,34, 340., 1.09);
    mm.defineBondStretch_KA( 1,35, 340., 1.09);
    mm.defineBondStretch_KA(22,31, 553., 0.9600);

    // I'm making this one up -- couldn't find O2 in Amber99
    mm.defineChargedAtomType_KA(9999, "Sherm's O2", 25, 0); // must be neutral by symmetry
    mm.defineBondStretch_KA(25,25, 570., 1.21); // bond length is right, stiffness is from C=O.

    mm.defineBondBend_KA( 1, 1, 1, 40., 109.5);
    mm.defineBondBend_KA( 1, 1,22, 50., 109.5);
    mm.defineBondBend_KA( 1, 1,23, 50., 109.5);
    mm.defineBondBend_KA( 1, 1,34, 50., 109.5);
    mm.defineBondBend_KA( 1, 1,35, 50., 109.5);
    mm.defineBondBend_KA( 1,22,31, 55., 108.5);
    mm.defineBondBend_KA( 1,23, 1, 60., 109.5);
    mm.defineBondBend_KA(22, 1,23, 50., 109.5); // made up (sherm)
    mm.defineBondBend_KA(22, 1,35, 50., 109.5);
    mm.defineBondBend_KA(22,31, 1, 50., 109.5);
    mm.defineBondBend_KA(23, 1,35, 50., 109.5);
    mm.defineBondBend_KA(34, 1,34, 35., 109.5);
    mm.defineBondBend_KA(35, 1,35, 35., 109.5);

    mm.defineBondTorsion_KA( 1, 1, 1, 1, 1, 0.2,  180.,
                                         2, 0.25, 180.,
                                         3, 0.18,   0.);
    mm.defineBondTorsion_KA( 1, 1, 1,22, 3, 0.156,  0.);
    mm.defineBondTorsion_KA( 1, 1, 1,23, 3, 0.156,  0.);
    mm.defineBondTorsion_KA( 1, 1, 1,34, 3, 0.160,  0.);
    mm.defineBondTorsion_KA( 1, 1, 1,35, 3, 0.156,  0.);
    mm.defineBondTorsion_KA( 1, 1,22,31, 1, 0.025,  0.,
                                         3, 0.160,  0.);
    mm.defineBondTorsion_KA( 1, 1,23, 1, 2, 0.100,180.,
                                         3, 0.383,  0.);
    mm.defineBondTorsion_KA( 1,23, 1, 1, 2, 0.850,180.,
                                         3, 0.100,  0.);
    mm.defineBondTorsion_KA( 1,23, 1,22, 1, 1.350,180.,
                                         2, 0.850,180.,
                                         3, 0.100,  0.);
    mm.defineBondTorsion_KA( 1,23, 1,23, 2, 0.850,180.,
                                         3, 0.100,  0.);
    mm.defineBondTorsion_KA( 1,23, 1,35, 3, 0.383,  0.);
    mm.defineBondTorsion_KA(22, 1, 1,23, 2, 1.175,  0.,
                                         3, 0.144,  0.);
    mm.defineBondTorsion_KA(22, 1, 1,22, 2, 1.175,  0.,
                                         3, 0.144,  0.);
    mm.defineBondTorsion_KA(22, 1, 1,35, 3, 0.156,  0.);
    mm.defineBondTorsion_KA(23, 1, 1,35, 3, 0.156,  0.);
    mm.defineBondTorsion_KA(23, 1,22,31, 3, 0.156,  0.); // made up
    mm.defineBondTorsion_KA(31,22, 1,35, 3, 0.167,  0.);

    mm.defineBondTorsion_KA(34, 1, 1,34, 3, 0.150,  0.);
    mm.defineBondTorsion_KA(35, 1, 1,35, 3, 0.156,  0.);

    mm.setVdwMixingRule( DuMMForceFieldSubsystem::LorentzBerthelot );

    // These are just for playing around with the force field terms.
    mm.setVdwGlobalScaleFactor(1);
    mm.setCoulombGlobalScaleFactor(1);
    mm.setBondStretchGlobalScaleFactor(1);
    mm.setBondBendGlobalScaleFactor(1);
    mm.setBondTorsionGlobalScaleFactor(1);


    //mbs.addForceSubsystem(gravity);

    /*

    const OneDofEthane ethane1(allowStretch, GroundId, mbs);
    const OneDofEthane ethane2(allowStretch, ethane1.getBodyId(0), mbs);
    const OneDofEthane ethane3(allowStretch, ethane2.getBodyId(0), mbs);
    const OneDofEthane ethane4(allowStretch, ethane3.getBodyId(0), mbs);
    const RigidEthane  rethane1(0, GroundId, mbs);
    const RigidEthane  rethane2(60, GroundId, mbs);
    */
    const bool allowStretch = false;
    const OneDofEthane ethane1(allowStretch, GroundId, mbs);
    //const RigidEthane  rethane1(0, GroundId, mbs);
    const RigidEthane  rethane2(60, GroundId, mbs);
    const FloppyEthane floppy1(GroundId, mbs);
    const RigidO2      rigidO2(GroundId, mbs);

    const CartesianRibose cribose(GroundId, mbs);

    /* Cartesian:  
    for (int i=0; i < mm.getNAtoms(); ++i) {
        MobilizedBodyId b = ethane.addRigidBody(
            MassProperties(mm.getAtomMass(i), Vec3(0), Inertia(0)), Transform(),
            GroundId, Transform(),
            Mobilizer::Cartesian());
        mm.attachAtomToBody(i, b, Vec3(0));
    }
    /**/

    //if (useCartesian && useRigid && wantConstraint) {
    //    int theConstraint =
    //           ethane.addConstantDistanceConstraint(firstCartesianBody-1, Vec3(0),
    //                                                firstCartesianBody, Vec3(0), 1.5);
    //    DecorativeLine purpleLine; purpleLine.setColor(Purple).setLineThickness(3);
    //    artwork.addRubberBandLine(firstCartesianBody-1, Vec3(0),
    //                              firstCartesianBody, Vec3(0), purpleLine);
    //}

    DecorativeLine crossBodyBond; crossBodyBond.setColor(Orange).setLineThickness(5);

    for (DuMM::BondId i = (DuMM::BondId)0; i < (DuMM::BondId)mm.getNBonds(); ++i) {
        const DuMM::AtomId    a1 = mm.getBondAtom(i,0), a2 = mm.getBondAtom(i,1);
        const MobilizedBodyId b1 = mm.getAtomBody(a1),  b2 = mm.getAtomBody(a2);
        if (b1==b2)
            artwork.addBodyFixedDecoration(b1, Transform(),
                                           DecorativeLine(mm.getAtomStationOnBody(a1), mm.getAtomStationOnBody(a2))
                                             .setColor(Gray).setLineThickness(3));
        else
            artwork.addRubberBandLine(b1, mm.getAtomStationOnBody(a1),
                                      b2, mm.getAtomStationOnBody(a2), crossBodyBond);
    }

    for (DuMM::AtomId anum = (DuMM::AtomId)0; anum < (DuMM::AtomId)mm.getNAtoms(); ++anum) {
        Real shrink = 0.25, opacity = mm.getAtomElement(anum)==1?0.5:1;
        Real r = mm.getAtomRadius(anum);
        if (r<.001) r=0.1; //nm
        //opacity=0.5;//XXX
        artwork.addBodyFixedDecoration(mm.getAtomBody(anum), mm.getAtomStationOnBody(anum),
            DecorativeSphere(shrink*r)
                .setColor(mm.getAtomDefaultColor(anum)).setOpacity(opacity).setResolution(3));
    }

    State s = mbs.realizeTopology();
    //matter.setUseEulerAngles(s,true);
    mbs.realizeModel(s);
   // gravity.setZeroHeight(s, -100);

    cribose.setDefaultInternalState(s);
    cribose.setMoleculeTransform(s, Transform(Rotation::aboutZ(Pi/2), Vec3(0,1,0)));

    floppy1.setDefaultInternalState(s);
    floppy1.setMoleculeTransform(s,Vec3(-1,0,0));
    //floppy1.setCCStretch(.1,s);
    //floppy1.setTorsionAngleDeg(80,s);
    //floppy1.setTorsionRate(10,s);

    ethane1.setDefaultInternalState(s);
    ethane1.setMoleculeTransform(s,Vec3(1,0,0));

    //rethane1.setDefaultInternalState(s);
   // rethane1.setMoleculeTransform(s,Vec3(0,0,-1));

    rethane2.setDefaultInternalState(s);
    rethane2.setMoleculeTransform(s,Vec3(-1,0,-1));

    rigidO2.setDefaultInternalState(s);

    const Transform o2pos( Rotation::aboutXThenNewY(0.5*Pi/2, 0.5*Pi/2),
                           Vec3(1,0,-1));
   // rigidO2.setMoleculeTransform(s,o2pos);
   // rigidO2.setMoleculeVelocity(s,SpatialVec(0*Vec3(1.1,1.2,3), Vec3(-.2,0,0)));

    /*

    if (allowStretch) ethane1.setCCStretch(0.03, s);
    ethane1.setTorsionAngleDeg(5, s);

    if (allowStretch) ethane2.setCCStretch(0.03, s);
    ethane2.setMoleculeTransform(s,Vec3(0,1,0));

    if (allowStretch) ethane3.setCCStretch(-0.03, s);
    ethane3.setMoleculeTransform(s,Transform(Rotation::aboutZ(Pi/2),Vec3(1,0,1)),);

    if (allowStretch) ethane4.setCCStretch(-0.03, s);
    ethane4.setMoleculeTransform(s,Vec3(-1,0,0));


   */

    /* Cartesian: 
    for (int i=0; i < mm.getNAtoms(); ++i) {
        int b = mm.getAtomBody(i);
        ethane.setMobilizerTransform(s, b, 
            Transform(mm.getAtomStationInCluster(i, wholeEthaneEclipsed)));
    }
    /**/


    mm.dump();

    VTKReporter display(mbs, 0.1);

    RungeKuttaMerson study(mbs, s);
    //CPodesIntegrator study(mbs,s);

    const Real h = outputInterval;
    const int interval = 1;
    const Real tstart = 0.;
    const Real tmax = simulationLength; //ps

    s.updTime() = tstart;
    display.report(s);

    study.setAccuracy(accuracy);
    study.initialize(); 

    std::vector<State> saveEm;
    saveEm.push_back(s);
    for (int i=0; i<100; ++i)
        saveEm.push_back(s);    // delay
    display.report(s);

    const Real Estart = mbs.getEnergy(s);

    int step = 0; bool flag=false;
    while (s.getTime() <= tmax) {

        mbs.realize(s);

        cout << s.getTime();
        cout << " deltaE=" << 100*(mbs.getEnergy(s)-Estart)
                                /(std::abs(Estart)+TinyReal) 
             << "% pe(kcal)=" << mbs.getPotentialEnergy(s)*KJ2Kcal
             << ", ke(kcal)=" << mbs.getKineticEnergy(s)*KJ2Kcal
             << " hNext(fs)=" << 1000*study.getPredictedNextStep();

        cout << "\n  System COM loc=" << matter.calcSystemMassCenterLocationInGround(s);
        cout << "\n  System COM vel=" << matter.calcSystemMassCenterVelocityInGround(s);
        cout << "\n  System COM acc=" << matter.calcSystemMassCenterAccelerationInGround(s);
        cout << endl;

        cout << "     q=" << matter.getQ(s) << endl;
        cout << "     u=" << matter.getU(s) << endl;
        cout << "  udot=" << matter.getUDot(s) << endl;

        cout << endl;

        if (!(step % interval)) {
            display.report(s);
            saveEm.push_back(s);
        }

        study.step(s.getTime() + h);
        ++step;
    }

    while(true) {
        for (int i=0; i < (int)saveEm.size(); ++i) {
            display.report(saveEm[i]);
            //display.report(saveEm[i]); // half speed
        }
        getchar();
    }

  }
catch (const std::exception& e)
  {
    printf("EXCEPTION THROWN: %s\n", e.what());
  }
    return 0;
}

EthaneMolecule::EthaneMolecule(MobilizedBodyId pId, const Transform& parentTransform,
                               MolecularMechanicsSystem& mmSys)
  : Molecule(pId,parentTransform,mmSys)
{
    SimbodyMatterSubsystem&  matter = 
        SimbodyMatterSubsystem::updDowncast(mmSys.updMatterSubsystem());
    DuMMForceFieldSubsystem& mm     = mmSys.updMolecularMechanicsForceSubsystem();

    twoCarbons = methyl[0] = methyl[1] = DuMM::InvalidClusterId;

    // Create the atoms and bonds. Atom 0 is C0, atom 1 is C1, 2-4 are the
    // hydrogens attached to C0, 5-7 are the hydrogens attached to C1.
    for (int i=0;i<2;++i) atoms.push_back(mm.addAtom((DuMM::ChargedAtomTypeId)A99EthaneCarbon));
    for (int i=0;i<6;++i) atoms.push_back(mm.addAtom((DuMM::ChargedAtomTypeId)A99EthaneHydrogen));

    mm.addBond(getC(0),getC(1));
    for (int c=0; c<2; ++c)
        for (int h=0; h<3; ++h)
            mm.addBond(getC(c),getH(c,h));

    // Define the clusers.
    twoCarbons = mm.createCluster("two carbons");
    methyl[0] =  mm.createCluster("methyl 0");
    methyl[1] =  mm.createCluster("methyl 1");
    
    // The "twoCarbons" cluster looks like this:        
    //          y
    //          |
    //          C0 --> ---- C1
    //         /     x
    //        z 
    // That is, the 1st carbon is at the origin, the 2nd is out along the +x
    // axis by the nominal C-C bond length.

    mm.placeAtomInCluster(getC(0), twoCarbons, Vec3(0));
    mm.placeAtomInCluster(getC(1), twoCarbons, Vec3(getNominalCCBondLength(),0,0));

    // Now build two identical methyl clusters. We'll worry about getting them
    // oriented properly when we place them onto bodies.
    // The methyl clusters should look like this:
    //
    //          H0     
    //           \   y
    //            \  |
    //             . C --> x
    //      (H2) .  /      
    //         *   z    
    //       H1
    //
    //
    // That is, H0 is in the (-x,+y) plane, tipped by the nominal
    // H-C-C bend angle. Then H1 is the H0 vector 
    // rotated +120 degrees about x (that is, out of the screen).
    // H2 is the H0 vector rotated 240 (=-120) degrees about x (into the
    // screen, not shown).

    const Vec3 H1pos = Rotation::aboutZ(getNominalHCCBondAngle())
                          * Vec3(getNominalCHBondLength(),0,0);

    for (int c=0; c<2; ++c) {
        mm.placeAtomInCluster(getC(c), methyl[c], Vec3(0));
        for (int h=0; h<3; ++h) {
            const Vec3 Hpos = Rotation::aboutX(h*120*Deg2Rad) * H1pos;
            mm.placeAtomInCluster(getH(c,h), methyl[c], Hpos);
        }
    }
}

OneDofEthane::OneDofEthane(bool allowStretch, MobilizedBodyId pId, MolecularMechanicsSystem& mmSys)
  : EthaneMolecule(pId,Transform(),mmSys)
{
    SimbodyMatterSubsystem&  matter = 
        SimbodyMatterSubsystem::updDowncast(mmSys.updMatterSubsystem());
    DuMMForceFieldSubsystem& mm     = mmSys.updMolecularMechanicsForceSubsystem();

    MobilizedBody& parent = matter.updMobilizedBody(parentId);

    const Rotation PinAboutX = Rotation::aboutY(90*Deg2Rad); // move z to +x

    // Mount the methyls onto bodies, methyl[0] first. Connect
    // them by either a pin or cylinder depending on allowStretch.

    const MassProperties m0mp(mm.calcClusterMassProperties(methyl[0], Transform()));
    const MassProperties m1mp(mm.calcClusterMassProperties(methyl[1], Transform()));

    bodies.push_back(
        MobilizedBody::Free( parent, Transform(),    // parent mobilizer frmae
            Body::Rigid(m0mp), Transform()));           // inboard mobilizer frame

    if (allowStretch) {
        bodies.push_back(
            MobilizedBody::Cylinder(
                matter.updMobilizedBody(getBodyId(0)), 
                Transform(PinAboutX, Vec3(getNominalCCBondLength(),0,0)),
                Body::Rigid(m1mp), Transform(PinAboutX, Vec3(0))));
    } else {
        bodies.push_back(
            MobilizedBody::Pin(
                matter.updMobilizedBody(getBodyId(0)),
                Transform(PinAboutX, Vec3(getNominalCCBondLength(),0,0)),
                Body::Rigid(m1mp), Transform(PinAboutX, Vec3(0))));
    }

    mm.attachClusterToBody(methyl[0], bodies[0], Transform());
    mm.attachClusterToBody(methyl[1], bodies[1], Transform(Rotation::aboutY(180*Deg2Rad)));
}

RigidEthane::RigidEthane(Real torsionAngleInDeg, MobilizedBodyId pId, MolecularMechanicsSystem& mmSys)
  : EthaneMolecule(pId,Transform(),mmSys)
{
    SimbodyMatterSubsystem&  matter = 
        SimbodyMatterSubsystem::updDowncast(mmSys.updMatterSubsystem());
    DuMMForceFieldSubsystem& mm     = mmSys.updMolecularMechanicsForceSubsystem();

    MobilizedBody& parent = matter.updMobilizedBody(parentId);

    const DuMM::ClusterId wholeEthaneCluster = mm.createCluster("rigid ethane");

    // If we choose to treat the entire ethane molecule as a rigid body, we'll align 
    // the 1st methyl group's reference frame with the body frame, and transform the
    // second by rotating it 180 degrees about y and shifting it by the nominal C-C
    // bond length in the +x direction. We'll then rotate about x to produce
    // a desired conformation.    
    //    H00                        H10
    //      \   y            y1 z1   /
    //       \  |             | /  /
    //          C0 --> -- <-- C1  
    //         /     x    x1
    //        z 
    mm.placeClusterInCluster(methyl[0], wholeEthaneCluster,  Transform());
    mm.placeClusterInCluster(methyl[1], wholeEthaneCluster, 
        Transform(Rotation::aboutYThenOldX(180*Deg2Rad, torsionAngleInDeg),
                  Vec3(getNominalCCBondLength(),0,0)));

    // Align cluster reference frame with body's.
    bodies.push_back(
        MobilizedBody::Free(parent, Transform(),    // parent mobilizer frame Mb
            Body::Rigid(mm.calcClusterMassProperties(wholeEthaneCluster, Transform())), 
            Transform()));                          // body mobilizer frame M

    mm.attachClusterToBody(wholeEthaneCluster, bodies[0], Transform()); 
}

// We will orient the carbon atom's body frames in opposite 
// directions, rotating 180 degrees about Y.
//    H00                        H10
//      \   y0           y1 z1   /
//       \  |             | /  /
//          C0 --> -- <-- C1  
//         /     x0   x1
//        z0 
// The cylinder joint between C0 & C1 will be defined along the +x0 direction (-x1).
// The bend-stretch joints connecting the H's to their C's will bend around their C's
// z axis, so that a bend of +109.5 yields the nominal H angle on either side. The
// stretch coordinate then operates along the outboard (child) body's NEW x axis
// direction.
//
// Warning: the reference configuration (all coordinates 0) has all the atoms
// on top of one another; don't realize past the Model stage until you have called
// setDefaultInternalState().
//
// The "molecule frame" is considered to be identical with C0's body frame.
//
FloppyEthane::FloppyEthane(MobilizedBodyId pId, MolecularMechanicsSystem& mmSys)
  : EthaneMolecule(pId,Transform(),mmSys)
{
    SimbodyMatterSubsystem&  matter = 
        SimbodyMatterSubsystem::updDowncast(mmSys.updMatterSubsystem());
    DuMMForceFieldSubsystem& mm     = mmSys.updMolecularMechanicsForceSubsystem();

    MobilizedBody& parent = matter.updMobilizedBody(parentId);

    // For C-C cylinder joint; rotation and translation are about the
    // Mobilizer frames' common Z axis.
    const Transform C0CylMobFrame(Rotation::aboutY( 90*Deg2Rad)); // move z to x0 direction
    const Transform C1CylMobFrame(Rotation::aboutY(-90*Deg2Rad)); // move z to -x1 direction
    const Transform HMobFrame;  // same as body frame for all H's

    // C0 is our base body, attached to parent by 6 dof joint
    bodies.push_back(
        MobilizedBody::Free(
            parent, Transform(),    // use parent body frame as reference
            Body::Rigid(MassProperties(mm.getAtomMass(getC(0)), Vec3(0), Inertia(0))),
            Transform()));          // use C0 body frame for Free mobilizer
    mm.attachAtomToBody(getC(0), bodies.back());

    // C1 is body 1, connected to C0 by a cylinder joint
    bodies.push_back(
        MobilizedBody::Cylinder(
            matter.updMobilizedBody(bodies[0]), C0CylMobFrame, 
            Body::Rigid(MassProperties(mm.getAtomMass(getC(1)), Vec3(0), Inertia(0))),
            C1CylMobFrame));
    mm.attachAtomToBody(getC(1), bodies.back());
           
    // Now attach 3 Hs to each C.
    for (int c=0; c<2; ++c) {
        const MobilizedBodyId Cbody = mm.getAtomBody(getC(c));
        for (int h=0; h<3; ++h) {
            const Transform CBendStretchMob(Rotation::aboutX(h*120*Deg2Rad));
            bodies.push_back(
                MobilizedBody::BendStretch(
                    matter.updMobilizedBody(Cbody), CBendStretchMob,
                    Body::Rigid(MassProperties(mm.getAtomMass(getH(c,h)), Vec3(0), Inertia(0))),
                    HMobFrame));
            mm.attachAtomToBody(getH(c,h), bodies.back());
        }
    }
}

void FloppyEthane::setDefaultInternalState(State& s) const {
    // All the bodies (1 per atom) except the 0th are internal. We'll set all internal q's
    // to the appropriate nominal values for generic EthaneMolecules, and
    // set all the internal u's to zero

    // C1
    const MobilizedBodyId CBody = getDuMM().getAtomBody(getC(1));
    const MobilizedBody& b = getMatter().getMobilizedBody(CBody);
    b.setOneQ(s, 0, 0); // torsion;
    b.setOneQ(s, 1, getNominalCCBondLength()); // stretch
    b.setOneU(s, 0, 0); // torsion rate
    b.setOneU(s, 1, 0); // stretch rate

    // H (bend,stretch)
    for (int c=0; c<2; ++c)
        for (int h=0; h<3; ++h) {
            const MobilizedBodyId HBody = getDuMM().getAtomBody(getH(c,h));
            const MobilizedBody& b = getMatter().getMobilizedBody(HBody);
            b.setOneQ(s, 0, getNominalHCCBondAngle());
            b.setOneQ(s, 1, getNominalCHBondLength());
            b.setOneU(s, 0, 0);
            b.setOneU(s, 1, 0);
        }
}
