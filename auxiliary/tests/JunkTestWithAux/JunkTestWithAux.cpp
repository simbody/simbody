/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-7 Stanford University and the Authors.         *
 * Authors: Randy Radmer, Michael Sherman                                     *
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
using std::cout;
using std::endl;

using namespace SimTK;
using namespace DuMM; // for conversion constants

//static const int NUM_WATERS = 23;
static const int NUM_WATERS = 20;

enum {
    ATOM_TYPE_OW = 1,
    ATOM_TYPE_HW = 2,
    ATOM_TYPE_CA_BASE = 11,
    ATOM_TYPE_HA_BASE = 21,
    CHARGED_ATOM_TYPE_OW = 1001,
    CHARGED_ATOM_TYPE_HW = 1002,
    CHARGED_ATOM_TYPE_CA_BASE = 1011,
    CHARGED_ATOM_TYPE_HA_BASE = 1021,
};

class Molecule {
public:
    Molecule(int parentBodyNum, const Transform& parentMobilizerFrame,
             const MolecularMechanicsSystem& mmSys)
        : parent(parentBodyNum), mobilizerFrameOnParent(parentMobilizerFrame), 
          mmSystem(&mmSys)
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
    MobilizedBodyId getBody(int i) const {return bodies[i];}

    const SimbodyMatterSubsystem& getMatter() const {
        return SimbodyMatterSubsystem::downcast(mmSystem->getMatterSubsystem());
    }
    const DuMMForceFieldSubsystem& getDuMM() const {
        return mmSystem->getMolecularMechanicsForceSubsystem();
    }
protected:
    std::vector<DuMM::AtomId>    atoms;
    std::vector<MobilizedBodyId> bodies;
    MobilizedBodyId       parent;
    Transform mobilizerFrameOnParent;
    const MolecularMechanicsSystem* mmSystem;
};


class Tip3p_water : public Molecule {
public:
    Tip3p_water(MobilizedBodyId parentBodyNum,
          MolecularMechanicsSystem& mmSys)
      : Molecule(parentBodyNum, Transform(), mmSys)
    {
        Real xH, yH;
        DuMMForceFieldSubsystem& mm = mmSys.updMolecularMechanicsForceSubsystem();

        // Create the atoms and bonds. Atom 0 is OW, atoms 1 and 2 are H1 and H2
        //  OW will serve as the base frame for the molecule.
        atoms.push_back(mm.addAtom((DuMM::ChargedAtomTypeId)CHARGED_ATOM_TYPE_OW));
        atoms.push_back(mm.addAtom((DuMM::ChargedAtomTypeId)CHARGED_ATOM_TYPE_HW));
        atoms.push_back(mm.addAtom((DuMM::ChargedAtomTypeId)CHARGED_ATOM_TYPE_HW));
        mm.addBond(getAtom(0),getAtom(1));
        mm.addBond(getAtom(0),getAtom(2));
        //mm.addBond(getAtom(1),getAtom(2));

        // Define the clusers.
        tip3p_water = mm.createCluster("tip3p_water");

        xH=getNominalHHBondLength()/2;
        yH=sqrt(getNominalOHBondLength()*getNominalOHBondLength()-xH*xH);
        //printf("2*xH = %f, sqrt(xH**2+yH**2) = %f\n", 2*xH, sqrt(xH*xH+yH*yH));
        mm.placeAtomInCluster(getAtom(0), tip3p_water, Vec3(0));
        mm.placeAtomInCluster(getAtom(1), tip3p_water, Vec3(-xH,yH,0));
        mm.placeAtomInCluster(getAtom(2), tip3p_water, Vec3( xH,yH,0));

        SimbodyMatterSubsystem&  matter =
            SimbodyMatterSubsystem::updDowncast(mmSys.updMatterSubsystem());

        MassProperties mprops = mm.calcClusterMassProperties(tip3p_water, Transform());

        bodies.push_back(
            MobilizedBody::Free(
                matter.updMobilizedBody(parentBodyNum), Transform(Vec3(0,0,0)),    // parent mobilizer frame
                Body::Rigid(mprops), Transform(Vec3(0,0,0))));    // inboard mobilizer frame

        mm.attachClusterToBody(tip3p_water, bodies.back(), Transform());
    }

    Real getNominalOHBondLength() const { return 0.9572 * Ang2Nm; }
    Real getNominalHHBondLength() const { return 1.5136 * Ang2Nm; }

    void setDefaultInternalState(State& s) const { } // none

protected:
    DuMM::ClusterId tip3p_water; // cluster
};


class Benzene : public Molecule {
public:
    Benzene(MobilizedBodyId parentBodyNum, const Transform& parentMobilizerFrame,
                  MolecularMechanicsSystem& mmSys)
      : Molecule(parentBodyNum, parentMobilizerFrame, mmSys)
    {
        Real a, x, y, rC, rH;
        DuMMForceFieldSubsystem& mm = mmSys.updMolecularMechanicsForceSubsystem();

        for(int i=0; i<6; i++) {
            atoms.push_back( mm.addAtom((DuMM::ChargedAtomTypeId)(CHARGED_ATOM_TYPE_CA_BASE+i)) );
        }
        for(int i=0; i<6; i++) {
            atoms.push_back( mm.addAtom((DuMM::ChargedAtomTypeId)(CHARGED_ATOM_TYPE_HA_BASE+i)) );
        }
        for( int i = 0; i < 6; i++ ) {
            mm.addBond( (DuMM::AtomId)i, (DuMM::AtomId)((i+1)%6) );
            mm.addBond( (DuMM::AtomId)i, (DuMM::AtomId)(i+6) );
            //mm.addBond(i,(i+2)%6);  // False bonds -- hack for improper torsions
        }

        // Define the clusers.
        benzene = mm.createCluster("benzene");

        rC=getNominalCCBondLength();
        rH=rC+getNominalCCBondLength();
        for (int i=0; i<6; i++) {
            a=(2*Pi/6.0)*i;
            x=cos(a);
            y=sin(a);
            mm.placeAtomInCluster((DuMM::AtomId)(i), benzene, Vec3(rC*x, rC*y, 0));
            mm.placeAtomInCluster((DuMM::AtomId)(i+6), benzene, Vec3(rH*x, rH*y, 0));
        }
    }

    // Get the atom number for each atom.

    Real getNominalCCBondLength() const { return 1.40 * Ang2Nm; }
    Real getNominalCHBondLength() const { return 1.08 * Ang2Nm; }

protected:
    DuMM::ClusterId benzene; // cluster
};


class RigidBenzene : public Benzene {
public:
    RigidBenzene(MobilizedBodyId parent, MolecularMechanicsSystem& mmSys)
                : Benzene(parent,Transform(),mmSys)
    {
        SimbodyMatterSubsystem&  matter =
            SimbodyMatterSubsystem::updDowncast(mmSys.updMatterSubsystem());
        DuMMForceFieldSubsystem& mm     = mmSys.updMolecularMechanicsForceSubsystem();

        MassProperties mprops = mm.calcClusterMassProperties(benzene, Transform());

        bodies.push_back(
            MobilizedBody::Free(
                matter.updMobilizedBody(parent), Transform(Vec3(0,0,0)),    // parent mobilizer frame
                Body::Rigid(mprops), Transform(Vec3(0,0,0))));    // inboard mobilizer frame

        mm.attachClusterToBody(benzene, bodies.back(), Transform());
    }

    void setDefaultInternalState(State& s) const { } // none
};


class FloppyBenzene : public Benzene {
public:
    FloppyBenzene(MobilizedBodyId parent, MolecularMechanicsSystem& mmSys)
                : Benzene(parent,Transform(),mmSys)
    { }

    void setDefaultInternalState(State& s) const { } // none
};



static const Transform BodyFrame;   // identity transform on any body

// How it actually looks now:
int main() {
try
  { MolecularMechanicsSystem mbs;
    SimbodyMatterSubsystem   matter(mbs);
    DuMMForceFieldSubsystem  mm(mbs);
    GeneralForceElements     forces(mbs);
    DecorationSubsystem      artwork(mbs);

    Real accuracy = 1e-2;
    Real outputInterval = .10;
    Real simulationLength = 100;
    //Real outputInterval = .1;
    //Real simulationLength = 10;

    const Real torsControlGain = /*100000*/0;
    const Real desiredTorsAngle = /*Pi/3*/0;

    forces.addGlobalEnergyDrain(.1);

    mm.setVdw14ScaleFactor(1/2.); // reduce energy by these factors
    mm.setCoulomb14ScaleFactor(1/1.2);

    mm.defineAtomClass_KA(ATOM_TYPE_OW, "TIP3P OW", 8, 1, 1.7683, 0.1520); 
    mm.defineAtomClass_KA(ATOM_TYPE_HW, "TIP3P HW", 1, 1, 1.0000, 0.0000); 

    mm.defineChargedAtomType_KA(CHARGED_ATOM_TYPE_OW, "TIP3P OW",
                                ATOM_TYPE_OW, -0.834);
    mm.defineChargedAtomType_KA(CHARGED_ATOM_TYPE_HW, "TIP3P HW",
                                ATOM_TYPE_HW,  0.417);

    mm.defineBondStretch_KA(ATOM_TYPE_OW, ATOM_TYPE_HW, 553., 0.9572);
    mm.defineBondStretch_KA(ATOM_TYPE_HW, ATOM_TYPE_HW, 553., 1.5136);

// CA
    mm.defineAtomClass_KA(ATOM_TYPE_CA_BASE+0, "BENZENE C1", 6, 1, 1.9080, 0.0860);
    mm.defineAtomClass_KA(ATOM_TYPE_CA_BASE+1, "BENZENE C2", 6, 1, 1.9080, 0.0860);
    mm.defineAtomClass_KA(ATOM_TYPE_CA_BASE+2, "BENZENE C3", 6, 1, 1.9080, 0.0860);
    mm.defineAtomClass_KA(ATOM_TYPE_CA_BASE+3, "BENZENE C4", 6, 1, 1.9080, 0.0860);
    mm.defineAtomClass_KA(ATOM_TYPE_CA_BASE+4, "BENZENE C5", 6, 1, 1.9080, 0.0860);
    mm.defineAtomClass_KA(ATOM_TYPE_CA_BASE+5, "BENZENE C6", 6, 1, 1.9080, 0.0860);
// HA
    mm.defineAtomClass_KA(ATOM_TYPE_HA_BASE+0, "BENZENE H1", 1, 1, 1.4590, 0.0150);
    mm.defineAtomClass_KA(ATOM_TYPE_HA_BASE+1, "BENZENE H2", 1, 1, 1.4590, 0.0150);
    mm.defineAtomClass_KA(ATOM_TYPE_HA_BASE+2, "BENZENE H3", 1, 1, 1.4590, 0.0150);
    mm.defineAtomClass_KA(ATOM_TYPE_HA_BASE+3, "BENZENE H4", 1, 1, 1.4590, 0.0150);
    mm.defineAtomClass_KA(ATOM_TYPE_HA_BASE+4, "BENZENE H5", 1, 1, 1.4590, 0.0150);
    mm.defineAtomClass_KA(ATOM_TYPE_HA_BASE+5, "BENZENE H6", 1, 1, 1.4590, 0.0150);

// CA
    mm.defineChargedAtomType_KA(CHARGED_ATOM_TYPE_CA_BASE+0, "BENZENE C1", ATOM_TYPE_CA_BASE+0, -0.115);
    mm.defineChargedAtomType_KA(CHARGED_ATOM_TYPE_CA_BASE+1, "BENZENE C2", ATOM_TYPE_CA_BASE+1, -0.115);
    mm.defineChargedAtomType_KA(CHARGED_ATOM_TYPE_CA_BASE+2, "BENZENE C3", ATOM_TYPE_CA_BASE+2, -0.115);
    mm.defineChargedAtomType_KA(CHARGED_ATOM_TYPE_CA_BASE+3, "BENZENE C4", ATOM_TYPE_CA_BASE+3, -0.115);
    mm.defineChargedAtomType_KA(CHARGED_ATOM_TYPE_CA_BASE+4, "BENZENE C5", ATOM_TYPE_CA_BASE+4, -0.115);
    mm.defineChargedAtomType_KA(CHARGED_ATOM_TYPE_CA_BASE+5, "BENZENE C6", ATOM_TYPE_CA_BASE+5, -0.115);
// HA
    mm.defineChargedAtomType_KA(CHARGED_ATOM_TYPE_HA_BASE+0, "BENZENE H1", ATOM_TYPE_HA_BASE+0,  0.115);
    mm.defineChargedAtomType_KA(CHARGED_ATOM_TYPE_HA_BASE+1, "BENZENE H2", ATOM_TYPE_HA_BASE+1,  0.115);
    mm.defineChargedAtomType_KA(CHARGED_ATOM_TYPE_HA_BASE+2, "BENZENE H3", ATOM_TYPE_HA_BASE+2,  0.115);
    mm.defineChargedAtomType_KA(CHARGED_ATOM_TYPE_HA_BASE+3, "BENZENE H4", ATOM_TYPE_HA_BASE+3,  0.115);
    mm.defineChargedAtomType_KA(CHARGED_ATOM_TYPE_HA_BASE+4, "BENZENE H5", ATOM_TYPE_HA_BASE+4,  0.115);
    mm.defineChargedAtomType_KA(CHARGED_ATOM_TYPE_HA_BASE+5, "BENZENE H6", ATOM_TYPE_HA_BASE+5,  0.115);

    for (int i=0; i<6; i++) {
//     CA-CA
        mm.defineBondStretch_KA( ATOM_TYPE_CA_BASE+i, ATOM_TYPE_CA_BASE+(1+i)%6, 469.0, 1.400);
//     CA-HA
        mm.defineBondStretch_KA( ATOM_TYPE_CA_BASE+i, ATOM_TYPE_HA_BASE+i, 367.0, 1.080);
//     CA-CA -- False bond, for Improper Torsions; note these are zero
        mm.defineBondStretch_KA( ATOM_TYPE_CA_BASE+i, ATOM_TYPE_CA_BASE+(2+i)%6, 0.0, 1.732);

//     CA-CA-CA
        mm.defineBondBend_KA( ATOM_TYPE_CA_BASE+i, ATOM_TYPE_CA_BASE+(1+i)%6, ATOM_TYPE_CA_BASE+(2+i)%6, 63.0, 120.00);
//     CA-CA-HA
        mm.defineBondBend_KA( ATOM_TYPE_CA_BASE+i, ATOM_TYPE_CA_BASE+(1+i)%6, ATOM_TYPE_HA_BASE+(1+i)%6, 50.0, 120.00);
//     HA-CA-CA
        mm.defineBondBend_KA( ATOM_TYPE_HA_BASE+i, ATOM_TYPE_CA_BASE+i, ATOM_TYPE_CA_BASE+(1+i)%6, 50.0, 120.00);

//     CA-CA-CA-CA
        mm.defineBondTorsion_KA( ATOM_TYPE_CA_BASE+i, ATOM_TYPE_CA_BASE+(1+i)%6, ATOM_TYPE_CA_BASE+(2+i)%6, ATOM_TYPE_CA_BASE+(3+i)%6, 2, 3.625, 180.);
//     CA-CA-CA-HA
        mm.defineBondTorsion_KA( ATOM_TYPE_CA_BASE+i, ATOM_TYPE_CA_BASE+(1+i)%6, ATOM_TYPE_CA_BASE+(2+i)%6, ATOM_TYPE_HA_BASE+(2+i)%6, 2, 3.625, 180.);
//     HA-CA-CA-CA
        mm.defineBondTorsion_KA( ATOM_TYPE_HA_BASE+i, ATOM_TYPE_CA_BASE+i, ATOM_TYPE_CA_BASE+(1+i)%6, ATOM_TYPE_CA_BASE+(2+i)%6, 2, 3.625, 180.);
//     HA-CA-CA-HA
        mm.defineBondTorsion_KA( ATOM_TYPE_HA_BASE+i, ATOM_TYPE_CA_BASE+i, ATOM_TYPE_CA_BASE+(1+i)%6, ATOM_TYPE_HA_BASE+(1+i)%6, 2, 3.625, 180.);
//     CA-CA-CA-HA -- Improper
        mm.defineBondTorsion_KA( ATOM_TYPE_CA_BASE+i, ATOM_TYPE_CA_BASE+(2+i)%6, ATOM_TYPE_CA_BASE+(1+i)%6, ATOM_TYPE_HA_BASE+(1+i)%6, 2, 1.1, 180.);
    }

    mm.setVdwMixingRule( DuMMForceFieldSubsystem::LorentzBerthelot );


    const RigidBenzene rBenzene(GroundId, mbs);

    std::vector<Tip3p_water> tip3p_waters;
    for (int i=0; i<NUM_WATERS; ++i) {
        tip3p_waters.push_back(Tip3p_water(GroundId, mbs));
    }

    DecorativeLine crossBodyBond; crossBodyBond.setColor(Orange).setLineThickness(5);

    for (DuMM::BondId i = (DuMM::BondId)0; i < (DuMM::BondId)mm.getNBonds(); ++i) {
        const DuMM::AtomId a1 = mm.getBondAtom(i, 0), a2 = mm.getBondAtom(i, 1);
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
    mbs.realizeModel(s);

    rBenzene.setDefaultInternalState(s);
    rBenzene.setMoleculeTransform(s,Vec3(0,0,0));

    int iSize;
    iSize=int(sqrt(NUM_WATERS/2.0));
    for (int i=0; i<NUM_WATERS/2; ++i) {
        Real x, y=-0.5, z;

        x=.5*(i%iSize-iSize/2);
        z=.5*(i/iSize-iSize/2);
        tip3p_waters[i].setDefaultInternalState(s);
        tip3p_waters[i].setMoleculeTransform(s,Vec3(x,y,z));
    }
    for (int i=NUM_WATERS/2; i<NUM_WATERS; ++i) {
        Real x, y=0.5, z;

        x=.5*(i%iSize-iSize/2);
        z=.5*(i/iSize-iSize/2);
        tip3p_waters[i].setDefaultInternalState(s);
        tip3p_waters[i].setMoleculeTransform(s,Vec3(x,y,z));
    }

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


