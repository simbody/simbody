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
    SimbodyMatterSubsystem  pend;
    UniformGravitySubsystem gravity(Vec3(0,-g,0));
    GeneralForceElements    forces;
    
    MultibodySystem mbs;
    mbs.setMatterSubsystem(pend);
    mbs.addForceSubsystem(gravity);
    mbs.addForceSubsystem(forces);

    const Transform pinFrameOnGround = BodyFrame; // connect directly to ground frame z
    const Transform pinFrame(Vec3(0,d/2,0));      // joint frame in body frame (will use z)
    const Vec3      massStation(0,-d/2,0);        // where we'll put the point mass in the body frame

    // Collect mass, center of mass, and inertia 
    MassProperties mprops(m, massStation, InertiaMat(massStation, m));
    cout << "mprops about body frame: " << mprops.getMass() << ", " 
        << mprops.getCOM() << ", " << mprops.getInertia() << endl;

    const int pendBody =
      pend.addRigidBody(mprops, pinFrame,
                        Ground, pinFrameOnGround, 
                        Mobilizer::Pin   // rotation around z
                        );

    const Vec3 attachPt(1,-.5,0);
    const Real k = 10, theta0 = 90*RadiansPerDegree, c = 0.1;
    //forces.addMobilityLinearSpring(pendBody, 0, k, theta0);
    //forces.addMobilityLinearDamper(pendBody, 0, 10.);


    forces.addUserForce(new MySinusoid(pendBody, 0, 10000, 100, 0));

    forces.addTwoPointLinearSpring(Ground, attachPt, 
                                   pendBody, massStation,
                                   100000, 0);

    //forces.addTwoPointLinearDamper(Ground, Vec3(20,-20,0), 
    //                               pendBody, massStation,
    //                               1000);


    //forces.addTwoPointConstantForce(Ground, attachPt,
    //                                pendBody, massStation,
     //                               100000);

    //forces.addGlobalEnergyDrain(100);

    // want -g*sin(theta)*d = k*(theta-theta0)
    // so solve err(theta)=k*(theta-theta0)+g*d*sin(theta)=0
    // d err/d theta is k+g*d*cos(theta)
    Real guess = Pi/4;
    Real err = k*(guess-theta0)+g*d*std::sin(guess);
    while (std::abs(err) > 1e-10) {
        cout << "err=" << err << endl;
        Real deriv = k+g*d*std::cos(guess);
        guess -= err/deriv;
        err = k*(guess-theta0)+g*d*std::sin(guess);
    }

    cout << "Final angle should be: " << guess/RadiansPerDegree << " err=" << err << endl;

    VTKReporter display(mbs);
    DecorativeLine ln; ln.setColor(Magenta).setLineThickness(3);
    display.addRubberBandLine(Ground, attachPt, pendBody, massStation, ln);

    State s;
    mbs.realize(s, Stage::Built);

    RungeKuttaMerson study(mbs, s);

    display.report(s);

    const Real angleInDegrees = 30;
    pend.setMobilizerQ(s, pendBody, 0, angleInDegrees*RadiansPerDegree);

    display.report(s);

    const Real h = .0001;
    const Real tstart = 0.;
    const Real tmax = 100;

    study.setAccuracy(1e-4);

    study.initialize(); 
    display.report(s);
    s.updTime() = tstart;
    int step = 0;
    while (s.getTime() < tmax) {
        study.step(s.getTime() + h);
        printf("%6.2f %8.2f ", s.getTime(), 
            pend.getMobilizerQ(s,pendBody,0)/RadiansPerDegree);
        cout << " E=" << mbs.getEnergy(s)
             << " (pe=" << mbs.getPotentialEnergy(s)
             << ", ke=" << mbs.getKineticEnergy(s)
             << ") hNext=" << study.getPredictedNextStep() << endl;

        if (!(step % 10)) {
            display.report(s);
        }
        ++step;
    }




}
catch (const std::exception& e) {
    printf("EXCEPTION THROWN: %s\n", e.what());
}
return 0;
}
