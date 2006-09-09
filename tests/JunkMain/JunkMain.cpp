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
 * This is just a disposable outer block for tests during development of
 * Simbody. Don't include this in the nightly test suite.
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


static const Real Pi = std::acos(-1.), RadiansPerDegree = Pi/180;
static const int  Ground = 0;       // ground is always body 0
static const Transform BodyFrame;   // identity transform on any body

class MySinusoid: public GeneralForceElements::UserForce {
public:
    MySinusoid(int b, int d, const Real& amp, const Real& w, const Real& ph=0) 
      : body(b), dof(d), amplitude(amp), period(w), phase(ph)
    {
    }

    // Implementation of pure virtual.
    void calc(const MatterSubsystem& matter, 
              const State&           state,
              Vector_<SpatialVec>&   bodyForces,
              Vector_<Vec3>&         particleForces,
              Vector&                mobilityForces,
              Real&                  pe) const 
    {
        matter.addInMobilityForce(state,body,dof,
            amplitude*std::sin(2*Pi*period*state.getTime() + phase),
            mobilityForces);
    }

    // Implementation of pure virtual;
    GeneralForceElements::UserForce* clone() const { 
        return new MySinusoid(*this); 
    }
private:
    int  body, dof;
    Real amplitude, period, phase;
};


static const Real g = 9.8;  // Earth gravity in m/s^2
static const Real d = 0.5;  // length of pendulum (m)
static const Real m = 1.;   // mass of pendulum (kg)
#ifdef NOTDEF
// How I want it to look:
int main() {
try {
    // Create the Subsystems and add them to the System.
    SimbodyMatterSubsystem  pend;
    UniformGravitySubsystem gravity(Vec3(0,-g,0));

    MultibodySystem system;
    mbs.setMatterSubsystem(pend);
    mbs.addForceSubsystem(gravity);

    // Build the multibody system.
    const Real g = 9.8;  // Earth gravity in m/s^2
    const Real d = 0.5;  // length of pendulum (m)
    const Real m = 1.;   // mass of pendulum (kg)

    const Transform mobilizerFrame(Vec3(0,d/2,0));
    const Transform mobilizerFrameOnGround = BodyFrame;

    const int pendBody =
        pend.addRigidBody(massProps, mobilizerFrame,  // the body
                          Ground, GroundFrame,        // its parent
                          Mobilizer::Pin);            // aligns z axes


    // Create visualization reporter.
    VTKReporter display(system);

    // Create study and get writable access to its State.
    MultibodyDynamicsStudy study(system);
    State& state = study.updState();

    // Study parameters.
    const Real startAngle       = 30;   // degrees
    const Real reportInterval   = .01;  // s
    const Real simulationLength = 100;  // s
    const Real tStart           = 0;
    const Real tEnd             = tStart + simulationLength;
    const Real accuracy         = 1e-4;

    // Set study options.
    study.setAccuracy(accuracy);

    // Set initial state.
    state.setTime(tStart);
    pend.setMobilizerQ(state, pendBody, 0, startAngle*RadiansPerDegree);

    // Evaluate the system at the current state without performing
    // any analysis.
    study.realize();
    display.report(state);    // Let's see what it looks like.

    // Perform initial condition analyses if any. This will fail if it can't
    // satisfy all constraints to tolerance.
    study.initialize();
    display.report(state);    // Let's see what it looks like.

    // Step until we pass the end time.
    while (state.getTime() < tEnd) {
        study.step(std::min(state.getTime() + reportInterval, tEnd));
        cout << " E=" << mbs.getEnergy(state)
             << " (pe=" << mbs.getPotentialEnergy(state)
             << ", ke=" << mbs.getKineticEnergy(state)
             << ") hNext=" << study.getPredictedNextStep() << endl;
            
        display.report(state);
    }
}
catch (const std::exception& e) {
    printf("EXCEPTION THROWN: %s\n", e.what());
}
return 0;
}
#endif

// How it actually looks now:
int main() {
try {
    SimbodyMatterSubsystem   molecule;
    SimbodyMatterSubsystem   ethane;
    DuMMForceFieldSubsystem  mm;
    GeneralForceElements     forces;

    bool useRigid = true;
    const Real vdwFac = 1;
    const Real chgFac = 1;
    const Real stretchFac = 1;
    const Real bendFac = 1;
    const Real torsFac = 1;

    forces.addGlobalEnergyDrain(10);

    // AMBER 99

    mm.setVdw14ScaleFactor(2.0);
    mm.setCoulomb14ScaleFactor(1.2);
    mm.defineAtomClass(1,  "Amber99 CT", 6, 4, 1.9080, vdwFac*0.1094);
    mm.defineAtomClass(34, "Amber99 HC", 1, 1, 1.4870, vdwFac*0.0157); 
    mm.defineChargedAtomType(13, "Amber99 Alanine CB", 1, -0.1825*chgFac);
    mm.defineChargedAtomType(14, "Amber99 Alanine HB", 34, 0.0603*chgFac);
    mm.defineBondStretch(1,1,  stretchFac*310., 1.5260);
    mm.defineBondStretch(1,34, stretchFac*340., 1.09);
    mm.defineBondBend(1,1,34, bendFac*50, 109.5);
    mm.defineBondBend(34,1,34, bendFac*35, 109.5);
    mm.defineBondTorsion(34,1,1,34, 3, torsFac*0.150, 0);
    
    MultibodySystem mbs;
    mbs.setMatterSubsystem(useRigid ? ethane : molecule);
    mbs.addForceSubsystem(mm);
    mbs.addForceSubsystem(forces);

    // ethane:
    // atom 0 is carbon1
    // atoms 1,2,3 are attached to carbon1
    // atom 4 is carbon2
    // atoms 5,6,7 are attached to carbon2

    const int nAtoms = 8;
    const Real massH = 1.008, massC = 12.011;
    // This is the description of the joint between the two bodies. We
    // want to rotate about the X axes, shifted by the length of a C-C bond.
    // Pin joints rotate around Z, so we need to rotate the joint frames
    // by +90 degrees around Y.
    const Vec3 ccBond(1.53688, 0, 0);
    RotationMat ccJointFrame; ccJointFrame.setToBodyFixed123(Vec3(0,Pi/2,0));

    int  type[] = {13,14,14,14, 13,14,14,14};
    int  body[] = {2,2,2,2, 3,3,3,3};
    Real mass[] = {massH, massC, massC, massC,
                   massH, massC, massC, massC};
    Vec3 station[] = { Vec3(0), Vec3(-.3778,1.02422,0), Vec3(-.3778,-0.514034,-0.885898), Vec3(-.3778,-0.510199,0.888107),
                       Vec3(0), Vec3(.3778,0.510199,0.888107), Vec3(.3778,0.514034,-0.885898),
                                Vec3(.3778,-1.02422,0) };

    // Collect mass, center of mass, and inertia
    Real mass1=0, mass2=0;
    Vec3 com1(0), com2(0);
    InertiaMat iner1(0), iner2(0);
    for (int i=0; i<4; ++i) {
        mass1 += mass[i]; mass2 += mass[i+4];
        com1 += mass[i]*station[i]; com2 += mass[i+4]*station[i+4];
        iner1 += InertiaMat(station[i], mass[i]);
        iner2 += InertiaMat(station[i+4], mass[i+4]);
    }
    com1 /= mass1; com2 /= mass2;

    MassProperties mprops1(mass1,com1,iner1);
    MassProperties mprops2(mass2,com2,iner2);
    int b0 = ethane.addRigidBody(MassProperties(0,Vec3(0),InertiaMat(0)), Transform(),
                                Ground, Transform(),
                                Mobilizer::Ball);
    int b1 = ethane.addRigidBody(mprops1, Transform(),
                                b0, Transform(),
                                Mobilizer::Cartesian);
    int b2 = ethane.addRigidBody(mprops2, Transform(ccJointFrame),
                                b1, Transform(ccJointFrame, ccBond),
                                Mobilizer::Pin);

    MassProperties mH(massH, Vec3(0), InertiaMat(0));
    MassProperties mC(massC, Vec3(0), InertiaMat(0));
    for (int i=0; i<2; ++i) 
        ethane.addRigidBody(mC, Transform(),
                              Ground, Transform(),
                              Mobilizer::Cartesian);
    for (int i=0; i<6; ++i) 
        ethane.addRigidBody(mH, Transform(),
                              Ground, Transform(),
                              Mobilizer::Cartesian);

    VTKReporter display(mbs);

// Rigid
    int c1 = mm.addAtom(body[0], 13, station[0]);
    int h11 = mm.addAtom(body[1], 14, station[1]);
    int h12 = mm.addAtom(body[2], 14, station[2]);
    int h13 = mm.addAtom(body[3], 14, station[3]);

    int c2 = mm.addAtom(body[4], 13, station[4]);
    int h21 = mm.addAtom(body[5], 14, station[5]);
    int h22 = mm.addAtom(body[6], 14, station[6]);
    int h23 = mm.addAtom(body[7], 14, station[7]);

    mm.addBond(c1,c2);
    mm.addBond(c1,h11); mm.addBond(c1,h12); mm.addBond(c1,h13);
    mm.addBond(c2,h21); mm.addBond(c2,h22); mm.addBond(c2,h23);

// Cartesian
    int cc1 = mm.addAtom(b2+1, 13, Vec3(0));
    int cc2 = mm.addAtom(b2+2, 13, Vec3(0));
    int ch11 = mm.addAtom(b2+3, 14, Vec3(0));
    int ch12 = mm.addAtom(b2+4, 14, Vec3(0));
    int ch13 = mm.addAtom(b2+5, 14, Vec3(0));
    int ch21 = mm.addAtom(b2+6, 14, Vec3(0));
    int ch22 = mm.addAtom(b2+7, 14, Vec3(0));
    int ch23 = mm.addAtom(b2+8, 14, Vec3(0));


    mm.addBond(cc1,cc2);
    mm.addBond(cc1,ch11); mm.addBond(cc1,ch12); mm.addBond(cc1,ch13);
    mm.addBond(cc2,ch21); mm.addBond(cc2,ch22); mm.addBond(cc2,ch23);

/*
    DecorativeLine ln; ln.setColor(Magenta).setLineThickness(3);
    display.addRubberBandLine(1, Vec3(0), 2, Vec3(0), ln);
    for (int i=3; i<=5; ++i)
        display.addRubberBandLine(1, Vec3(0), i, Vec3(0), ln);
    for (int i=6; i<=8; ++i)
        display.addRubberBandLine(2, Vec3(0), i, Vec3(0), ln);
*/
    for (int anum=0; anum < mm.getNAtoms(); ++anum) {
        display.addDecoration(mm.getAtomBody(anum), Transform(mm.getAtomStation(anum)),
            DecorativeSphere(0.25*mm.getAtomRadius(anum))
                .setColor(mm.getAtomDefaultColor(anum)).setOpacity(1).setResolution(3));
    }

    State s;
    mbs.realize(s, Stage::Built);
    //molecule.setUseEulerAngles(s, true);

    mm.dump();


    RungeKuttaMerson study(mbs, s);

    display.report(s);

    ethane.setMobilizerU(s, b1, 1, 10);

    ethane.setMobilizerQ(s, b2, 0, Pi/3);
    ethane.setMobilizerU(s, b2, 0, 100);

    const Real yoffs = 4;

    ethane.setMobilizerQ(s, b2+1, 1, yoffs); // shift 2nd molecule up yoffs in y
    //ethane.setMobilizerU(s, b2+1, 1, -10);

    ethane.setMobilizerQ(s, b2+2, 0, 1.53688+.05/*distort bond*/);
    ethane.setMobilizerQ(s, b2+2, 1, yoffs);

    ethane.setMobilizerQ(s, b2+3, 0, -.3778);
    ethane.setMobilizerQ(s, b2+3, 1,  1.02422+yoffs);
    ethane.setMobilizerQ(s, b2+3, 2,  0);

    ethane.setMobilizerQ(s, b2+4, 0, -.3778);
    ethane.setMobilizerQ(s, b2+4, 1, -0.514034+yoffs);
    ethane.setMobilizerQ(s, b2+4, 2, -0.885898);

    ethane.setMobilizerQ(s, b2+5, 0, -.3778);
    ethane.setMobilizerQ(s, b2+5, 1, -0.510199+yoffs);
    ethane.setMobilizerQ(s, b2+5, 2, 0.888107);


    ethane.setMobilizerQ(s, b2+6, 0, .3778+1.53688);
    ethane.setMobilizerQ(s, b2+6, 1,  0.510199+yoffs);
    ethane.setMobilizerQ(s, b2+6, 2,  0.888107);

    ethane.setMobilizerQ(s, b2+7, 0, .3778+1.53688);
    ethane.setMobilizerQ(s, b2+7, 1, 0.514034+yoffs);
    ethane.setMobilizerQ(s, b2+7, 2, -0.885898);

    ethane.setMobilizerQ(s, b2+8, 0, .3778+1.53688);
    ethane.setMobilizerQ(s, b2+8, 1, -1.02422+yoffs);
    ethane.setMobilizerQ(s, b2+8, 2, 0);

    display.report(s);

    const Real h = .01;
    const int interval = 1;
    const Real tstart = 0.;
    const Real tmax = 10; //ps

    study.setAccuracy(1e-2);
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
             << " (pe=" << mbs.getPotentialEnergy(s)
             << ", ke=" << mbs.getKineticEnergy(s)
             << ") hNext=" << study.getPredictedNextStep() << endl;

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
