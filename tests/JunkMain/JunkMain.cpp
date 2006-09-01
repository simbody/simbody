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
    DuMMForceFieldSubsystem  mm;
    GeneralForceElements     forces;

    const Real vdwRad5  = 2.2;
    const Real vdwRad10 = 2.0;
    const Real vdwRad20 = 1.5;
    const Real vdwFac = 1;
    const Real stretchFac = 0;

    mm.defineAtomType(5,  20., vdwRad5, vdwFac*0.2, 0);
    mm.defineAtomType(10, 14., vdwRad10, vdwFac*0.2, -1);
    mm.defineAtomType(20, 12., vdwRad20, vdwFac*0.3, 1);

    mm.defineBondStretch(5,5,   stretchFac*300.,2.5);
    mm.defineBondStretch(10,10, stretchFac*300.,2.0);
    mm.defineBondStretch(20,20, stretchFac*300.,1.5);
    mm.defineBondStretch(5,10,  stretchFac*300.,1.5);
    mm.defineBondStretch(5,20,  stretchFac*300.,1.5);
    mm.defineBondStretch(20,10, stretchFac*300.,5);

    
    MultibodySystem mbs;
    mbs.setMatterSubsystem(molecule);
    mbs.addForceSubsystem(mm);
    mbs.addForceSubsystem(forces);

    forces.addGlobalEnergyDrain(50);


    // Collect mass, center of mass, and inertia 
    const Real mass=1;
    MassProperties mprops(mass, Vec3(0), InertiaMat(0));

    // Create some point mass bodies
    for (int i=0; i<5; ++i) 
        molecule.addRigidBody(mprops, Transform(),
                              Ground, Transform(),
                              Mobilizer::Cartesian);

    VTKReporter display(mbs);
    display.setDefaultBodyColor(1, Red);
    display.setDefaultBodyColor(3, Red);
    display.setDefaultBodyColor(2, Green);
    display.setDefaultBodyColor(4, Green);
    display.setDefaultBodyColor(5, Yellow);

    int a1 = mm.addAtom(1, 10, Vec3(0));
    int a2 = mm.addAtom(2, 20, Vec3(0));
    int a3 = mm.addAtom(3, 10, Vec3(0));
    int a4 = mm.addAtom(4, 20, Vec3(0));
    int a5 = mm.addAtom(5,  5, Vec3(0));

    display.addDecoration(1,Transform(),
        DecorativeSphere(vdwRad10).setOpacity(0.3).setResolution(3));
    display.addDecoration(2,Transform(),
        DecorativeSphere(vdwRad20).setOpacity(1).setResolution(3));
    display.addDecoration(3,Transform(),
        DecorativeSphere(vdwRad10).setOpacity(0.3).setResolution(3));
    display.addDecoration(4,Transform(),
        DecorativeSphere(vdwRad20).setOpacity(1).setResolution(3));
    display.addDecoration(5,Transform(),
        DecorativeSphere(vdwRad5).setOpacity(.9).setResolution(3));

    mm.addBond(a1,a2); 
    mm.addBond(a1,a3); 
    mm.addBond(a2,a3);
    DecorativeLine ln; ln.setColor(Magenta).setLineThickness(3);
    display.addRubberBandLine(1, Vec3(0), 2, Vec3(0), ln);
    display.addRubberBandLine(1, Vec3(0), 3, Vec3(0), ln);
    display.addRubberBandLine(2, Vec3(0), 3, Vec3(0), ln);

    State s;
    mbs.realize(s, Stage::Built);

    mm.dump();


    RungeKuttaMerson study(mbs, s);

    display.report(s);

    const Real d = 0.6*(vdwRad10+vdwRad20);
    molecule.setMobilizerQ(s, 1, 0,  2*d);
    molecule.setMobilizerQ(s, 2, 0,  2*d);
    molecule.setMobilizerQ(s, 2, 1, 2*d);
    molecule.setMobilizerQ(s, 3, 0,  -2*d);
    molecule.setMobilizerQ(s, 4, 0,  -2*d);
    molecule.setMobilizerQ(s, 4, 1, -2*d);
    molecule.setMobilizerQ(s, 4, 2, 0.5*d);
    molecule.setMobilizerQ(s, 5, 0,  -d);
    molecule.setMobilizerQ(s, 5, 1,  2*d);
    molecule.setMobilizerQ(s, 5, 2,  -2*d);

    //molecule.setMobilizerU(s, 2, 0, -100);

    display.report(s);

    const Real h = .001;
    const int interval = 10;
    const Real tstart = 0.;
    const Real tmax = 5; //ps

    study.setAccuracy(1e-3);
    study.initialize(); 

    std::vector<State> saveEm;
    saveEm.push_back(s);
    for (int i=0; i<100; ++i)
        saveEm.push_back(s);    // delay
    display.report(s);

    s.updTime() = tstart;
    int step = 0;
    while (s.getTime() < tmax) {
        study.step(s.getTime() + h);

        cout << s.getTime();
        cout << " E=" << mbs.getEnergy(s)
             << " (pe=" << mbs.getPotentialEnergy(s)
             << ", ke=" << mbs.getKineticEnergy(s)
             << ") hNext=" << study.getPredictedNextStep() << endl;

        if (!(step % interval)) {
            display.report(s);
            saveEm.push_back(s);
        }
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
catch (const std::exception& e) {
    printf("EXCEPTION THROWN: %s\n", e.what());
}
return 0;
}
