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

static const Real Pi      = (Real)SimTK_PI;

 // multiply to convert
static const Real& Deg2Rad = DuMMForceFieldSubsystem::Deg2Rad;
static const Real& Rad2Deg = DuMMForceFieldSubsystem::Rad2Deg;
static const Real& Ang2Nm  = DuMMForceFieldSubsystem::Ang2Nm;
static const Real& Nm2Ang  = DuMMForceFieldSubsystem::Nm2Ang;
static const Real& Kcal2KJ = DuMMForceFieldSubsystem::Kcal2KJ;
static const Real& KJ2Kcal = DuMMForceFieldSubsystem::KJ2Kcal;

// mapping from atoms to Amber99 Charged Atom Types
enum {
    A99EthaneCarbon   = 13,    // closest I could find in Amber99; not really right
    A99EthaneHydrogen = 14,
    ShermDoubleBondedOxygen = 9999
};

class Molecule {
public:
    Molecule(int parentBodyNum, const Transform& parentMobilizerFrame,
             const MolecularMechanicsSystem& mmSys)
        : parent(parentBodyNum), mobilizerFrameOnParent(parentMobilizerFrame), 
          mmSystem(mmSys)
    { }

    // Translate and rotation the molecule as a whole. Usually this just means move the base body,
    // but some molecules may not have a single base, so they can override this default.
    virtual void setMoleculeTransform(State& s, const Transform& pos) const
    {
        getMatter().setMobilizerTransform(s, bodies[0], pos);
    }

    virtual void setMoleculeVelocity(State& s, const SpatialVec& vel) const
    {
        getMatter().setMobilizerVelocity(s, bodies[0], vel);
    }

    // This routine must set the internal mobilities to their nominal values, both
    // for position and velocity states. This State must have already been realized
    // to at least Model stage.
    virtual void setDefaultInternalState(State& s) const = 0; 

    int getNAtoms()  const {return (int)atoms.size();}
    int getNBodies() const {return (int)bodies.size();}

    // return atomId of ith atom in Molecule
    int getAtom(int i) const {return atoms[i];}

    // return bodyNum of ith body; 0 is molecule's base body
    BodyId getBody(int i) const {return bodies[i];}

    const SimbodyMatterSubsystem& getMatter() const {
        return SimbodyMatterSubsystem::downcast(mmSystem.getMatterSubsystem());
    }
    const DuMMForceFieldSubsystem& getDuMM() const {
        return mmSystem.getMolecularMechanicsForceSubsystem();
    }
protected:
    std::vector<int>    atoms;
    std::vector<BodyId> bodies;
    BodyId       parent;
    Transform mobilizerFrameOnParent;
    const MolecularMechanicsSystem& mmSystem;
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
    OxygenMolecule(BodyId parentBodyNum, const Transform& parentMobilizerFrame,
                   MolecularMechanicsSystem& mmSys) 
      : Molecule(parentBodyNum, parentMobilizerFrame, mmSys)
    {
        DuMMForceFieldSubsystem& mm = mmSys.updMolecularMechanicsForceSubsystem();
        
        // Create the atoms and bonds. Atom 0 is O0, atom 1 is O1. O0 will serve
        // as the base frame for the molecule.
        for (int i=0;i<2;++i) atoms.push_back(mm.addAtom(ShermDoubleBondedOxygen));
        mm.addBond(getO(0),getO(1));

        // Define the clusers.
        twoOxygens = mm.createCluster("two oxygens");

        mm.placeAtomInCluster(getO(0), twoOxygens, Vec3(0));
        mm.placeAtomInCluster(getO(1), twoOxygens, Vec3(0,0,getNominalOOBondLength()));
    }

    // Get the atom number for each oxygen.
    int getO(int i) const {assert(i==0||i==1); return getAtom(i);}

    Real getNominalOOBondLength() const {
        return 1.21 * Ang2Nm;
    }
protected:
    int twoOxygens; // cluster
};

class RigidO2 : public OxygenMolecule {
public:
    RigidO2(BodyId parent, MolecularMechanicsSystem& mmSys)
      : OxygenMolecule(parent,Transform(),mmSys)
    {
        SimbodyMatterSubsystem&  matter = 
            SimbodyMatterSubsystem::updDowncast(mmSys.updMatterSubsystem());
        DuMMForceFieldSubsystem& mm     = mmSys.updMolecularMechanicsForceSubsystem();

        // Align cluster reference frame with body's. (5 dofs!)
        // FreeLine prevents rotation about Z, so make sure the body has its
        // O=O axis arranged along Z (or rotate the frame here).
        /* This doesn't work: 
        bodies.push_back(
            matter.addRigidBody(
                MassProperties(0,Vec3(0),Inertia(0)),
                Transform(),            // inboard mobilizer frame
                parent, Transform(),    // parent mobilizer frame
                Mobilizer::Cartesian));
        // y
        bodies.push_back(
            matter.addRigidBody(
                MassProperties(0,Vec3(0),Inertia(0)),
                Transform(Rotation::aboutX(-90*Deg2Rad)),            // inboard mobilizer frame
                bodies.back(), Transform(Rotation::aboutX(-90*Deg2Rad)),    // parent mobilizer frame
                Mobilizer::Pin));
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
                Mobilizer::Pin));
        
        */
        MassProperties mprops = mm.calcClusterMassProperties(twoOxygens, Transform());
        MassProperties mpropsKludge(mprops.getMass(), mprops.getMassCenter(), 
                                    mprops.getInertia() + Inertia(0,0,.01));

        bodies.push_back(
            matter.addRigidBody(
                mprops,
                Transform(0*Vec3(0,0,.3)),            // inboard mobilizer frame
                parent, Transform(0*Vec3(.1,.7,.19)),    // parent mobilizer frame
                Mobilizer::FreeLine));

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
    EthaneMolecule(BodyId parentBodyNum, const Transform& parentMobilizerFrame,
                   MolecularMechanicsSystem&);

    // find the atoms
    int getC(int i) const {assert(i==0||i==1); return getAtom(i);}
    int getH(int whichCarbon, int whichHydrogen) const {
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
    int twoCarbons;
    int methyl[2];
};

class OneDofEthane : public EthaneMolecule {
public:
    OneDofEthane(bool allowStretch, BodyId parent, MolecularMechanicsSystem&);

    void setDefaultInternalState(State& s) const {
        const int ndof = getMatter().getDOF(getBody(1));
        for (int i=0; i<ndof; ++i) {
            getMatter().setMobilizerQ(s, getBody(1), i, 0);
            getMatter().setMobilizerU(s, getBody(1), i, 0);
        }
    }

    // Set stretch around the nominal length.
    void setCCStretch(Real stretchInNm, State& s) const {
        assert(getMatter().getDOF(getBody(1)) == 2);    // must have been build with Cylinder mobilizer
        const BodyId CBody = getDuMM().getAtomBody(getC(1));
        getMatter().setMobilizerQ(s, CBody, 1, stretchInNm);
    }

    void setTorsionAngleDeg(Real angleInDeg, State& s) const {
        const BodyId CBody = getDuMM().getAtomBody(getC(1));
        getMatter().setMobilizerQ(s, CBody, 0, angleInDeg*Deg2Rad);
    }

    // Rate is rad/ps
    void setTorsionRate(Real rateInRadPerPs, State& s) const {
        const BodyId CBody = getDuMM().getAtomBody(getC(1));
        getMatter().setMobilizerU(s, CBody, 0, rateInRadPerPs);
    }
};

class RigidEthane : public EthaneMolecule {
public:
    RigidEthane(Real torsionAngleInDeg, BodyId parent, MolecularMechanicsSystem&);
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
    FloppyEthane(BodyId parent, MolecularMechanicsSystem&);

    void setDefaultInternalState(State& s) const;

    // Set stretch around the nominal length.
    void setCCStretch(Real stretchInNm, State& s) const {
        const BodyId CBody = getDuMM().getAtomBody(getC(1));
        getMatter().setMobilizerQ(s, CBody, 1, stretchInNm);
    }

    void setTorsionAngleDeg(Real angleInDeg, State& s) const {
        const BodyId CBody = getDuMM().getAtomBody(getC(1));
        getMatter().setMobilizerQ(s, CBody, 0, angleInDeg*Deg2Rad);
    }

    // Rate is rad/ps
    void setTorsionRate(Real rateInRadPerPs, State& s) const {
        const BodyId CBody = getDuMM().getAtomBody(getC(1));
        getMatter().setMobilizerU(s, CBody, 0, rateInRadPerPs);
    }
};

static const Transform BodyFrame;   // identity transform on any body

// How it actually looks now:
int main() {
try
  { SimbodyMatterSubsystem   matter;
    DuMMForceFieldSubsystem  mm;
    GeneralForceElements     forces;

    Real accuracy = 1e-2;
    Real outputInterval = .05;
    Real simulationLength = 30;
    //Real outputInterval = .000000001;
    //Real simulationLength = .0000001;

    const Real torsControlGain = /*100000*/0;
    const Real desiredTorsAngle = /*Pi/3*/0;

    //forces.addGlobalEnergyDrain(0.02);


    // AMBER 99

    mm.setVdw14ScaleFactor(1/2.); // reduce energy by these factors
    mm.setCoulomb14ScaleFactor(1/1.2);

    mm.defineAtomClass_KA(1,  "Amber99 CT", 6, 4, 1.9080, 0.1094);
    mm.defineAtomClass_KA(2,  "Amber99 C",  6, 3, 1.9080, 0.0860);
    mm.defineAtomClass_KA(3,  "Amber99 CA", 6, 3, 1.9080, 0.0860);
    mm.defineAtomClass_KA(4,  "Amber99 CM", 6, 3, 1.9080, 0.0860);
    mm.defineAtomClass_KA(24, "Amber99 O",  8, 1, 1.6612, 0.2100);
    mm.defineAtomClass_KA(25, "Amber99 O2", 8, 1, 1.6612, 0.2100); 
    mm.defineAtomClass_KA(34, "Amber99 HC", 1, 1, 1.4870, 0.0157); 

    mm.defineChargedAtomType_KA(13, "Amber99 Alanine CB", 1, -0.1825);
    mm.defineChargedAtomType_KA(14, "Amber99 Alanine HB", 34, 0.0603);


    mm.defineBondStretch_KA(1,1,  310., 1.5260);
    mm.defineBondStretch_KA(1,34, 340., 1.09);

    // I'm making this one up -- couldn't find O2 in Amber99
    mm.defineChargedAtomType_KA(9999, "Sherm's O2", 25, 0); // must be neutral by symmetry
    mm.defineBondStretch_KA(25,25, 570., 1.21); // bond length is right, stiffness is from C=O.

    mm.defineBondBend_KA(1, 1,34, 50, 109.5);
    mm.defineBondBend_KA(34,1,34, 35, 109.5);

    mm.defineBondTorsion_KA(34,1,1,34, 3, 0.150, 0);

    mm.setVdwMixingRule( DuMMForceFieldSubsystem::LorentzBerthelot );

    // These are just for playing around with the force field terms.
    mm.setVdwGlobalScaleFactor(1);
    mm.setCoulombGlobalScaleFactor(1);
    mm.setBondStretchGlobalScaleFactor(1);
    mm.setBondBendGlobalScaleFactor(1);
    mm.setBondTorsionGlobalScaleFactor(1);


    MolecularMechanicsSystem mbs;
    mbs.setMatterSubsystem(matter);
    mbs.setMolecularMechanicsForceSubsystem(mm);
    mbs.addForceSubsystem(forces);
    UniformGravitySubsystem gravity(Vec3(0,.01,0));

    //mbs.addForceSubsystem(gravity);

    /*

    const OneDofEthane ethane1(allowStretch, GroundId, mbs);
    const OneDofEthane ethane2(allowStretch, ethane1.getBody(0), mbs);
    const OneDofEthane ethane3(allowStretch, ethane2.getBody(0), mbs);
    const OneDofEthane ethane4(allowStretch, ethane3.getBody(0), mbs);
    const RigidEthane  rethane1(0, GroundId, mbs);
    const RigidEthane  rethane2(60, GroundId, mbs);
    */
    const bool allowStretch = false;
    const OneDofEthane ethane1(allowStretch, GroundId, mbs);
    const RigidEthane  rethane1(0, GroundId, mbs);
    //const RigidEthane  rethane2(60, GroundId, mbs);
    const FloppyEthane floppy1(GroundId, mbs);
    const RigidO2      rigidO2(GroundId, mbs);

    /* Cartesian:  
    for (int i=0; i < mm.getNAtoms(); ++i) {
        BodyId b = ethane.addRigidBody(
            MassProperties(mm.getAtomMass(i), Vec3(0), Inertia(0)), Transform(),
            GroundId, Transform(),
            Mobilizer::Cartesian);
        mm.attachAtomToBody(i, b, Vec3(0));
    }
    /**/

    State s;
    mbs.realize(s, Stage::Topology);
    //matter.setUseEulerAngles(s,true);
    mbs.realize(s, Stage::Model);
   // gravity.setZeroHeight(s, -100);

    floppy1.setDefaultInternalState(s);
    //floppy1.setMoleculeTransform(s,Vec3(-3,0,0));
    //floppy1.setCCStretch(.1,s);
    //floppy1.setTorsionAngleDeg(80,s);
    //floppy1.setTorsionRate(10,s);

    ethane1.setDefaultInternalState(s);
    ethane1.setMoleculeTransform(s,Vec3(1,0,0));

    rethane1.setDefaultInternalState(s);
    rethane1.setMoleculeTransform(s,Vec3(0,0,-1));

    //rethane2.setDefaultInternalState(s);
    //rethane2.setMoleculeTransform(s,Vec3(-1,0,-1));

    rigidO2.setDefaultInternalState(s);

    const Transform o2pos( Rotation::aboutXThenNewY(0.5*Pi/2, 0.5*Pi/2),
                           Vec3(1,0,-1));
    rigidO2.setMoleculeTransform(s,o2pos);
    rigidO2.setMoleculeVelocity(s,SpatialVec(0*Vec3(1.1,1.2,3), Vec3(-.2,0,0)));

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

    //if (useCartesian && useRigid && wantConstraint) {
    //    int theConstraint =
    //           ethane.addConstantDistanceConstraint(firstCartesianBody-1, Vec3(0),
    //                                                firstCartesianBody, Vec3(0), 1.5);
    //    DecorativeLine purpleLine; purpleLine.setColor(Purple).setLineThickness(3);
    //    display.addRubberBandLine(firstCartesianBody-1, Vec3(0),
    //                              firstCartesianBody, Vec3(0), purpleLine);
    //}

    DecorativeLine crossBodyBond; crossBodyBond.setColor(Orange).setLineThickness(5);

    for (int i=0; i<mm.getNBonds(); ++i) {
        const int    a1 = mm.getBondAtom(i,0), a2 = mm.getBondAtom(i,1);
        const BodyId b1 = mm.getAtomBody(a1),  b2 = mm.getAtomBody(a2);
        if (b1==b2)
            display.addDecoration(b1, Transform(),
                                  DecorativeLine(mm.getAtomStationOnBody(a1), mm.getAtomStationOnBody(a2))
                                    .setColor(Gray).setLineThickness(3));
        else
            display.addRubberBandLine(b1, mm.getAtomStationOnBody(a1),
                                      b2, mm.getAtomStationOnBody(a2), crossBodyBond);
    }

    for (int anum=0; anum < mm.getNAtoms(); ++anum) {
        Real shrink = 0.25, opacity = mm.getAtomElement(anum)==1?0.5:1;
        //opacity=0.5;//XXX
        display.addDecoration(mm.getAtomBody(anum), mm.getAtomStationOnBody(anum),
            DecorativeSphere(shrink*mm.getAtomRadius(anum))
                .setColor(mm.getAtomDefaultColor(anum)).setOpacity(opacity).setResolution(3));
    }


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

    int step = 0;
    while (s.getTime() <= tmax) {
        mbs.realize(s);

        cout << s.getTime();
        cout << " deltaE=" << 100*(mbs.getEnergy(s)-Estart)
                                /(std::abs(Estart)+NTraits<Real>::Tiny) 
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

EthaneMolecule::EthaneMolecule(BodyId parent, const Transform& parentTransform,
                               MolecularMechanicsSystem& mmSys)
  : Molecule(parent,parentTransform,mmSys)
{
    SimbodyMatterSubsystem&  matter = 
        SimbodyMatterSubsystem::updDowncast(mmSys.updMatterSubsystem());
    DuMMForceFieldSubsystem& mm     = mmSys.updMolecularMechanicsForceSubsystem();

   twoCarbons = methyl[0] = methyl[1] = -1;

    // Create the atoms and bonds. Atom 0 is C0, atom 1 is C1, 2-4 are the
    // hydrogens attached to C0, 5-7 are the hydrogens attached to C1.
    for (int i=0;i<2;++i) atoms.push_back(mm.addAtom(A99EthaneCarbon));
    for (int i=0;i<6;++i) atoms.push_back(mm.addAtom(A99EthaneHydrogen));

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

OneDofEthane::OneDofEthane(bool allowStretch, BodyId parent, MolecularMechanicsSystem& mmSys)
  : EthaneMolecule(parent,Transform(),mmSys)
{
    SimbodyMatterSubsystem&  matter = 
        SimbodyMatterSubsystem::updDowncast(mmSys.updMatterSubsystem());
    DuMMForceFieldSubsystem& mm     = mmSys.updMolecularMechanicsForceSubsystem();

    const Rotation PinAboutX = Rotation::aboutY(90*Deg2Rad); // move z to +x

    // Mount the methyls onto bodies, methyl[0] first. Connect
    // them by either a pin or cylinder depending on allowStretch.

    bodies.push_back(
        matter.addRigidBody(
            mm.calcClusterMassProperties(methyl[0], Transform()),
            Transform(),            // inboard mobilizer frame
            parent, Transform(),    // parent mobilizer frmae
            Mobilizer::Free));

    bodies.push_back(
        matter.addRigidBody(
            mm.calcClusterMassProperties(methyl[1], Transform()),      
            Transform(PinAboutX, Vec3(0)), 
            getBody(0), Transform(PinAboutX, Vec3(getNominalCCBondLength(),0,0)),
            allowStretch ? Mobilizer::Cylinder : Mobilizer::Pin));

    mm.attachClusterToBody(methyl[0], bodies[0], Transform());
    mm.attachClusterToBody(methyl[1], bodies[1], Transform(Rotation::aboutY(180*Deg2Rad)));
}

RigidEthane::RigidEthane(Real torsionAngleInDeg, BodyId parent, MolecularMechanicsSystem& mmSys)
  : EthaneMolecule(parent,Transform(),mmSys)
{
    SimbodyMatterSubsystem&  matter = 
        SimbodyMatterSubsystem::updDowncast(mmSys.updMatterSubsystem());
    DuMMForceFieldSubsystem& mm     = mmSys.updMolecularMechanicsForceSubsystem();

    const int wholeEthaneCluster = mm.createCluster("rigid ethane");

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
        matter.addRigidBody(
            mm.calcClusterMassProperties(wholeEthaneCluster, Transform()), 
            Transform(),            // inboard mobilizer frame
            parent, Transform(),    // parent mobilizer frmae
            Mobilizer::Free));

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
FloppyEthane::FloppyEthane(BodyId parent, MolecularMechanicsSystem& mmSys)
  : EthaneMolecule(parent,Transform(),mmSys)
{
    SimbodyMatterSubsystem&  matter = 
        SimbodyMatterSubsystem::updDowncast(mmSys.updMatterSubsystem());
    DuMMForceFieldSubsystem& mm     = mmSys.updMolecularMechanicsForceSubsystem();

    // For C-C cylinder joint; rotation and translation are about the
    // Mobilizer frames' common Z axis.
    const Transform C0CylMobFrame(Rotation::aboutY( 90*Deg2Rad)); // move z to x0 direction
    const Transform C1CylMobFrame(Rotation::aboutY(-90*Deg2Rad)); // move z to -x1 direction
    const Transform HMobFrame;  // same as body frame for all H's

    // C0 is our base body, attached to parent by 6 dof joint
    bodies.push_back(
        matter.addRigidBody(
            MassProperties(mm.getAtomMass(getC(0)), Vec3(0), Inertia(0)),
            Transform(),            // use C0 body frame for Free mobilizer
            parent, Transform(),    // use parent body frame as reference
            Mobilizer::Free));
    mm.attachAtomToBody(getC(0), bodies.back());

    // C1 is body 1, connected to C0 by a cylinder joint
    bodies.push_back(
        matter.addRigidBody(
            MassProperties(mm.getAtomMass(getC(1)), Vec3(0), Inertia(0)),
            C1CylMobFrame,
            bodies[0], C0CylMobFrame,
            Mobilizer::Cylinder));
    mm.attachAtomToBody(getC(1), bodies.back());
           
    // Now attach 3 Hs to each C.
    for (int c=0; c<2; ++c) {
        const BodyId Cbody = mm.getAtomBody(getC(c));
        for (int h=0; h<3; ++h) {
            const Transform CBendStretchMob(Rotation::aboutX(h*120*Deg2Rad));
            bodies.push_back(
                matter.addRigidBody(
                    MassProperties(mm.getAtomMass(getH(c,h)), Vec3(0), Inertia(0)),
                    HMobFrame,
                    Cbody, CBendStretchMob,
                    Mobilizer::BendStretch));
            mm.attachAtomToBody(getH(c,h), bodies.back());
        }
    }
}

void FloppyEthane::setDefaultInternalState(State& s) const {
    // All the bodies (1 per atom) except the 0th are internal. We'll set all internal q's
    // to the appropriate nominal values for generic EthaneMolecules, and
    // set all the internal u's to zero

    // C1
    const BodyId CBody = getDuMM().getAtomBody(getC(1));
    getMatter().setMobilizerQ(s, CBody, 0, 0); // torsion;
    getMatter().setMobilizerQ(s, CBody, 1, getNominalCCBondLength()); // stretch
    getMatter().setMobilizerU(s, CBody, 0, 0); // torsion rate
    getMatter().setMobilizerU(s, CBody, 1, 0); // stretch rate

    // H (bend,stretch)
    for (int c=0; c<2; ++c)
        for (int h=0; h<3; ++h) {
            const BodyId HBody = getDuMM().getAtomBody(getH(c,h));
            getMatter().setMobilizerQ(s, HBody, 0, getNominalHCCBondAngle());
            getMatter().setMobilizerQ(s, HBody, 1, getNominalCHBondLength());
            getMatter().setMobilizerU(s, HBody, 0, 0);
            getMatter().setMobilizerU(s, HBody, 1, 0);
        }
}
