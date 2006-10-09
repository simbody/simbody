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

#include "SimTKcommon.h"
#include "Simmatrix.h"
#include "Simbody.h"

#include "simbody/internal/DecorativeGeometry.h"
#include "simbody/internal/VTKReporter.h"
#include "simbody/internal/NumericalMethods.h"

#include "simbody/internal/DuMMForceFieldSubsystem.h"

#include <string>
#include <iostream>
#include <exception>
#include <cmath>
using std::cout;
using std::endl;

using namespace SimTK;


static const Real Pi = NTraits<Real>::Pi, RadiansPerDegree = Pi/180;
static const Real EnergyUnitsPerKcal = 418.4; // exact 
static const int  Ground = 0;       // ground is always body 0
static const Transform BodyFrame;   // identity transform on any body

// How it actually looks now:
int main() {
try {
    SimbodyMatterSubsystem   ethane;
    DuMMForceFieldSubsystem  mm;
    GeneralForceElements     forces;

    bool useRigid = true, useCartesian = true, wantConstraint=false;
    const Real vdwFac = 1;
    const Real chgFac = 1;
    const Real stretchFac =1;
    const Real bendFac = 1;
    const Real torsFac = 1;

    const Real torsControlGain = /*100000*/0;
    const Real desiredTorsAngle = /*Pi/3*/0;

    forces.addGlobalEnergyDrain(20);


    // AMBER 99

    mm.setVdw14ScaleFactor(1/2.); // reduce energy by these factors
    mm.setCoulomb14ScaleFactor(1/1.2);

    mm.defineAtomClass(1,  "Amber99 CT", 6, 4, 1.9080, vdwFac*0.1094);
    mm.defineAtomClass(34, "Amber99 HC", 1, 1, 1.4870, vdwFac*0.0157); 
    mm.defineChargedAtomType(13, "Amber99 Alanine CB", 1, -0.1825*chgFac);
    mm.defineChargedAtomType(14, "Amber99 Alanine HB", 34, 0.0603*chgFac);
    mm.defineBondStretch(1,1,  stretchFac*310., 1.5260);
    mm.defineBondStretch(1,34, stretchFac*340., 1.09);
    mm.defineBondBend(1,1,34, bendFac*50, 109.5);
    mm.defineBondBend(34,1,34, bendFac*35, 109.5);
    mm.defineBondTorsion(34,1,1,34, 3, torsFac*0.150, 0);

    mm.setVdwMixingRule( DuMMForceFieldSubsystem::LorentzBerthelot );


    MultibodySystem mbs;
    mbs.setMatterSubsystem(ethane);
    mbs.addForceSubsystem(mm);
    mbs.addForceSubsystem(forces);

    // ethane:
    // atom 0 is carbon1
    // atoms 1,2,3 are attached to carbon1
    // atom 4 is carbon2
    // atoms 5,6,7 are attached to carbon2

    // rigid clusters:
    //   group 1: the two carbons
    //   group 2: carbon 1 (atom 0) and hydrogens 1,2,3
    //   group 3: carbon 2 (atom 4) and hydrogens 5,6,7
    //   group 4: the entire molecule
    // Any cluster or individual atom can be assigned to a body, provided
    // the resulting set of assignments represents a partitioning of
    // the atoms across the bodies.

    const int twoCarbons           = mm.createCluster("two carbons");
    const int methyl1              = mm.createCluster("methyl 1");
    const int methyl2              = mm.createCluster("methyl 2");
    const int wholeEthaneEclipsed  = mm.createCluster("ethaneEclipsed");
    const int wholeEthaneStaggered = mm.createCluster("ethaneStaggered");

    const Real ccNominalBondLength = 1.53688; // A
    const Real chNominalBondLength = 1.09;    // A
    const Real hccNominalBondBend  = 109.5*RadiansPerDegree;

    // Create the atoms and bonds. H[0..2] are attached to C[0], the others to C[1].
    int C[2]; for (int i=0; i<2; ++i) C[i] = mm.addAtom(13);
    int H[6]; for (int i=0; i<6; ++i) H[i] = mm.addAtom(14);
    mm.addBond(C[0],C[1]);
    for (int i=0; i<3; ++i) {mm.addBond(H[i],C[0]); mm.addBond(H[i+3],C[1]);}

    // Now build clusters. The "twoCarbons" cluster looks like this:        
    //          y
    //          |
    //          C0 --> ---- C1
    //         /     x
    //        z 
    // That is, the 1st carbon is at the origin, the 2nd is out along the +x
    // axis by the nominal C-C bond length.

    mm.placeAtomInCluster(C[0], twoCarbons, Vec3(0));
    mm.placeAtomInCluster(C[1], twoCarbons, Vec3(ccNominalBondLength,0,0));

    // Now build two identical methyl clusters. We'll worry about getting them
    // oriented properly when we place them into larger clusters or onto bodies.
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
    mm.placeAtomInCluster(C[0], methyl1, Vec3(0));
    mm.placeAtomInCluster(C[1], methyl2, Vec3(0));

    const Vec3 H1pos = RotationMat::aboutZ(hccNominalBondBend)
                          * Vec3(chNominalBondLength,0,0);
    for (int i=0; i<3; ++i) {
        const Vec3 Hpos = RotationMat::aboutX(i*120*RadiansPerDegree) * H1pos;
        mm.placeAtomInCluster(H[i],   methyl1, Hpos);
        mm.placeAtomInCluster(H[3+i], methyl2, Hpos);
    }

    // If we choose to treat the entire ethane molecule as a rigid body, we'll align 
    // the 1st methyl group's reference frame with the body frame, and transform the
    // second by rotating it 180 degrees about y and shifting it by the nominal C-C
    // bond length in the +x direction. We'll then rotate about x to produce
    // a staggered conformation.    
    //    H00                        H10
    //      \   y            y1 z1   /
    //       \  |             | /  /
    //          C0 --> -- <-- C1  
    //         /     x    x1
    //        z 
    mm.placeClusterInCluster(methyl1, wholeEthaneEclipsed,  Transform());
    mm.placeClusterInCluster(methyl1, wholeEthaneStaggered, Transform());

    mm.placeClusterInCluster(methyl2, wholeEthaneEclipsed, 
        Transform(RotationMat::aboutY(180*RadiansPerDegree),
                  Vec3(ccNominalBondLength,0,0)));
    mm.placeClusterInCluster(methyl2, wholeEthaneStaggered, 
        Transform(RotationMat::aboutYThenOldX(180*RadiansPerDegree, 60*RadiansPerDegree),
                  Vec3(ccNominalBondLength,0,0)));

    cout << "mass props twoCarbons =" << mm.calcClusterMassProperties(twoCarbons, Transform(Vec3(.76844,1,0)));
    cout << "mass props methyl1    =" << mm.calcClusterMassProperties(methyl1);
    cout << "mass props methyl2    =" << mm.calcClusterMassProperties(methyl2);
    cout << "mass props methyl2(rot-45y) =" << mm.calcClusterMassProperties(methyl2,
        Transform(RotationMat::aboutY(-45*RadiansPerDegree)));
    cout << "mass props eclipsed   =" << mm.calcClusterMassProperties(wholeEthaneEclipsed);
    cout << "mass props staggered  =" << mm.calcClusterMassProperties(wholeEthaneStaggered);


    /* Whole ethane as a rigid body, eclipsed. 
    // Align cluster reference frame with body's.
    int b1 = ethane.addRigidBody(
        mm.calcClusterMassProperties(wholeEthaneStaggered, Transform()), 
        Transform(),            // inboard mobilizer frame
        Ground, Transform(),    // parent mobilizer frmae
        Mobilizer::Free);
    mm.attachClusterToBody(wholeEthaneStaggered, b1, Transform()); 
    /**/

    /* Methyls connected by a torsion/stretch (cylinder) mobilizer.
    int b1 = ethane.addRigidBody(
                mm.calcClusterMassProperties(methyl1, Transform()),
                Transform(),            // inboard mobilizer frame
                Ground, Transform(),    // parent mobilizer frmae
                Mobilizer::Free);
    int b2 = ethane.addRigidBody(
                mm.calcClusterMassProperties(methyl2, Transform()),      
                Transform(RotationMat::aboutY(90*RadiansPerDegree), Vec3(0)), // move z to +x
                b1, Transform(RotationMat::aboutY(90*RadiansPerDegree), // move z to +x
                              Vec3(ccNominalBondLength,0,0)),
                Mobilizer::Cylinder);
    mm.attachClusterToBody(methyl1, b1, Transform());
    mm.attachClusterToBody(methyl2, b2, Transform(RotationMat::aboutY(180*RadiansPerDegree)));
    /**/

    /* Cartesian:  */
    for (int i=0; i < mm.getNAtoms(); ++i) {
        int b = ethane.addRigidBody(
            MassProperties(mm.getAtomMass(i), Vec3(0), InertiaMat(0)), Transform(),
            Ground, Transform(),
            Mobilizer::Cartesian);
        mm.attachAtomToBody(i, b, Vec3(0));
    }
    /**/

    State s;
    mbs.realize(s, Stage::Built);
    mbs.realize(s, Stage::Modeled);

    /* Cartesian: */
    for (int i=0; i < mm.getNAtoms(); ++i) {
        int b = mm.getAtomBody(i);
        ethane.setMobilizerConfiguration(s, b, 
            Transform(mm.getAtomStationInCluster(i, wholeEthaneEclipsed)));
    }
    /**/



    mm.dump();

    //const Vec3 ccBond(1.53688, 0, 0);
    //Vec3 station[] = { Vec3(0), /*Vec3(-.3778,1.02422,0)*/Vec3(0), 
    //                            /*Vec3(-.3778,-0.514034,-0.885898)*/Vec3(0), 
    //                            /*Vec3(-.3778,-0.510199,0.888107)*/Vec3(0),
    //                   Vec3(0), Vec3(.3778,0.510199,0.888107), 
    //                            Vec3(.3778,0.514034,-0.885898),
    //                            Vec3(.3778,-1.02422,0)};


    VTKReporter display(mbs);

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
        const int a1 = mm.getBondAtom(i,0), a2 = mm.getBondAtom(i,1);
        const int b1 = mm.getAtomBody(a1),  b2 = mm.getAtomBody(a2);
        if (b1==b2)
            display.addDecoration(b1, Transform(),
                                  DecorativeLine(mm.getAtomStationOnBody(a1), mm.getAtomStationOnBody(a2))
                                    .setColor(Gray).setLineThickness(3));
        else
            display.addRubberBandLine(b1, mm.getAtomStationOnBody(a1),
                                      b2, mm.getAtomStationOnBody(a2), crossBodyBond);
    }

    for (int anum=0; anum < mm.getNAtoms(); ++anum) {
        display.addDecoration(mm.getAtomBody(anum), Transform(mm.getAtomStationOnBody(anum)),
            DecorativeSphere(0.25*mm.getAtomRadius(anum))
                .setColor(mm.getAtomDefaultColor(anum)).setOpacity(0.25).setResolution(3));
    }



    RungeKuttaMerson study(mbs, s);

    display.report(s);

    /*
    // Give the whole rigid body molecule an initial velocity.
    //ethane.setMobilizerVelocity(s, b1, SpatialVec(Vec3(0), Vec3(0,10,0)));
    // Apply position and velocity directly to the joint axis for the torsion
    // between the two carbons.
    if (useRigid) {
        ethane.setMobilizerQ(s, b2, 0,Pi/3);
        //ethane.setMobilizerU(s, b2, 0, 10);
        ethane.setMobilizerQ(s, bh1, 0, Pi/2); // 1st axis is bend
        ethane.setMobilizerQ(s, bh2, 0, Pi/2);
        ethane.setMobilizerQ(s, bh3, 0, Pi/2);
        ethane.setMobilizerQ(s, bh1, 1, 1.); // 2nd axis is slider
        ethane.setMobilizerQ(s, bh2, 1, 1.);
        ethane.setMobilizerQ(s, bh3, 1, 1.);

        ethane.setMobilizerQ(s, bh4, 0, -Pi/2); // 1st axis is bend
        ethane.setMobilizerQ(s, bh5, 0, -Pi/2);
        ethane.setMobilizerQ(s, bh6, 0, -Pi/2);
        ethane.setMobilizerQ(s, bh4, 1, 1.); // 2nd axis is slider
        ethane.setMobilizerQ(s, bh5, 1, 1.);
        ethane.setMobilizerQ(s, bh6, 1, 1.);
    }

    if (useCartesian) {
        // shift 2nd molecule up yoffs in y
        const Real yoffs = 4;

        ethane.setMobilizerConfiguration(s, firstCartesianBody+0, Transform(Vec3(0,yoffs,0)));
        //ethane.setMobilizerU(s, firstCartesianBody+0, 1, -10);

        // distort bond a little
        ethane.setMobilizerConfiguration(s, firstCartesianBody+1, Transform(Vec3(1.53688+.05, yoffs, 0)));

        ethane.setMobilizerConfiguration(s, firstCartesianBody+2, Transform(Vec3(-.3778, 1.02422 +yoffs, 0)));
        ethane.setMobilizerConfiguration(s, firstCartesianBody+3, Transform(Vec3(-.3778,-0.514034+yoffs,-0.885898)));
        ethane.setMobilizerConfiguration(s, firstCartesianBody+4, Transform(Vec3(-.3778,-0.510199+yoffs, 0.888107)));

        ethane.setMobilizerConfiguration(s, firstCartesianBody+5, Transform(Vec3( .3778+1.53688, 0.510199+yoffs, 0.888107)));
        ethane.setMobilizerConfiguration(s, firstCartesianBody+6, Transform(Vec3( .3778+1.53688, 0.514034+yoffs,-0.885898)));
        ethane.setMobilizerConfiguration(s, firstCartesianBody+7, Transform(Vec3( .3778+1.53688,-1.02422 +yoffs, 0)));
    }
    */

    display.report(s);

    const Real h = /*.0025*/0.01;
    const int interval = 1;
    const Real tstart = 0.;
    const Real tmax = 5; //ps

    study.setAccuracy(1e-7);
    study.initialize(); 

    std::vector<State> saveEm;
    saveEm.push_back(s);
    for (int i=0; i<100; ++i)
        saveEm.push_back(s);    // delay
    display.report(s);

    const Real Estart = mbs.getEnergy(s);

    s.updTime() = tstart;
    int step = 0;
    while (s.getTime() < tmax) {
        study.step(s.getTime() + h);

        cout << s.getTime();
        cout << " deltaE=" << (mbs.getEnergy(s)-Estart)/Estart
             << " (pe=" << mbs.getPotentialEnergy(s)/EnergyUnitsPerKcal
             << ", ke=" << mbs.getKineticEnergy(s)/EnergyUnitsPerKcal
             << ") hNext=" << study.getPredictedNextStep();
        //if (useRigid) {
        //    cout << " cctors=" << ethane.getMobilizerQ(s, b2, 0)/RadiansPerDegree
        //         << " ccstretch=" << ethane.getMobilizerQ(s, b2, 1)
        //         << " h1bend=" << ethane.getMobilizerQ(s, bh1, 0)/RadiansPerDegree
        //         << " h1stretch=" << ethane.getMobilizerQ(s, bh1, 1); // XXX
        //}
        cout << endl;

        if (!(step % interval)) {
            display.report(s);
            saveEm.push_back(s);
        }
        ++step;
    }
/*
    const Transform& c1X = molecule.getBodyConfiguration(s, 1);
    cout << "h11=" << ~c1X*molecule.getBodyConfiguration(s, 3) << endl;
    cout << "h12=" << ~c1X*molecule.getBodyConfiguration(s, 4) << endl;
    cout << "h13=" << ~c1X*molecule.getBodyConfiguration(s, 5) << endl;

    const Transform& c2X = molecule.getBodyConfiguration(s, 2);
    cout << "c2=" << ~c1X*c2X << endl;
    cout << "h21=" << ~c2X*molecule.getBodyConfiguration(s, 6) << endl;
    cout << "h22=" << ~c2X*molecule.getBodyConfiguration(s, 7) << endl;
    cout << "h23=" << ~c2X*molecule.getBodyConfiguration(s, 8) << endl;
*/
    while(true) {
        for (int i=0; i < (int)saveEm.size(); ++i) {
            display.report(saveEm[i]);
            //display.report(saveEm[i]); // half speed
        }
        getchar();
    }

}
catch (const std::exception& e) {
    printf("EXCEPTION THROWN: %s\n", e.what());
}
return 0;
}
