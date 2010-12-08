/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
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
 * This is a simple demo of the need for constraint stabilization.
 */

#include "SimTKsimbody.h"

#include <cmath>
#include <cstdio>
#include <exception>
#include <vector>

using namespace std;
using namespace SimTK;

static const Real Deg2Rad = (Real)SimTK_DEGREE_TO_RADIAN,
                  Rad2Deg = (Real)SimTK_RADIAN_TO_DEGREE;

static const Transform GroundFrame;

static const Real m = 1;   // kg
static const Real g = 9.8; // meters/s^2; apply in –y direction
static const Real d = 0.5; // meters


void ff(Vector& v) {
    v = 23.;
}
    
int main(int argc, char** argv) {
  try { // If anything goes wrong, an exception will be thrown.

        // CREATE MULTIBODY SYSTEM AND ITS SUBSYSTEMS
    MultibodySystem         mbs;

    SimbodyMatterSubsystem  crankRocker(mbs);
    GeneralForceSubsystem   forces(mbs);
    DecorationSubsystem     viz(mbs);
    Force::UniformGravity   gravity(forces, crankRocker, Vec3(0, -g, 0));

        // ADD BODIES AND THEIR MOBILIZERS
    Body::Rigid crankBody  = Body::Rigid(MassProperties(.1, Vec3(0), 0.1*Gyration::brick(1,3,.5)))
                                .addDecoration(Transform(), 
                                               DecorativeEllipsoid(0.1*Vec3(1,3,.4))
                                                    .setResolution(10)
                                                    .setOpacity(.2));
    Body::Rigid sliderBody = Body::Rigid(MassProperties(.2, Vec3(0), 0.2*Gyration::brick(1,5,.5)))
                                .addDecoration(Transform(), 
                                               DecorativeEllipsoid(0.2*Vec3(1,5,.4))
                                                    .setColor(Blue)
                                                    .setResolution(10)
                                                    .setOpacity(.2));
    Body::Rigid longBar = Body::Rigid(MassProperties(0.01, Vec3(0), 0.01*Gyration::cylinderAlongX(.1, 5)))
        .addDecoration(Rotation(Pi/2,ZAxis), DecorativeCylinder(.01, 1));

    MobilizedBody::Pin
        crank(crankRocker.Ground(), Transform(), crankBody, Transform(Vec3(0, .15, 0)));
    MobilizedBody::Pin
        slider(crankRocker.Ground(), Transform(Vec3(2.5,0,0)), sliderBody, Transform(Vec3(0, 1, 0)));

    MobilizedBody::Universal
        rightConn(crank, Transform(Rotation(-Pi/2,YAxis),Vec3(0,-.3,0)), 
                  longBar, Transform(Rotation(-Pi/2,YAxis),Vec3(-1,0,0)));

    Constraint::Ball ball(slider, Vec3(0,-1, 0), rightConn, Vec3(1,0,0));
    ball.setDefaultRadius(0.01);

    //Constraint::Rod rodl(leftConn, Vec3(0), rightConn, Vec3(-2.5,0,0), 2.5);
    //Constraint::Rod rodr(leftConn, Vec3(2.5,-1,0), rightConn, Vec3(-1.5,0,0), std::sqrt(2.));
    //Constraint::Rod rodo(leftConn, Vec3(1.5,0,0), rightConn, Vec3(-2.5,1,0), std::sqrt(2.));
    //Constraint::Rod rodl(leftConn, Vec3(-2.5,0,0), rightConn, Vec3(-2.5,0,0), 5);
    //Constraint::Rod rodr(leftConn, Vec3(2.5,0,0), rightConn, Vec3(2.5,0,0), 5);
    //Constraint::Rod rodo(leftConn, Vec3(2.5,-1,0), rightConn, Vec3(-2.5,1,0), 2);

    // Add visualization line (orange=spring, black=constraint)
    //viz.addRubberBandLine(crank, Vec3(0,-3,0),
     //                     slider, Vec3(0,-5,0),
       //                   DecorativeLine().setColor(Black).setLineThickness(4));

    //forces.addMobilityConstantForce(leftPendulum, 0, 20);
    //forces.addCustomForce(ShermsForce(leftPendulum,rightPendulum));
    //forces.addGlobalEnergyDrain(1);

    Force::MobilityConstantForce(forces, crank, 0, 1);
    Force::MobilityLinearDamper(forces, crank, 0, 1.0);


    //Motion::Linear(crank, Vec3(a,b,c)); // crank(t)=at^2+bt+c
    //Motion::Linear lmot(rightConn, Vec3(a,b,c)); // both axes follow 
    //lmot.setAxis(1, Vec3(d,e,f));
    //Motion::Orientation(someBall, orientFuncOfT);
    //someBall.prescribeOrientation(orientFunc);
    //Motion::Relax(crank); // acc=vel=0, pos determined by some default relaxation solver

    // Like this, or should this be an Instance-stage mode of every mobilizer?
    //Motion::Lock(crank); // acc=vel=0; pos is discrete or fast
    //Motion::Steady(crank, Vec3(1,2,3)); // acc=0; vel=constant; pos integrated
    //Motion::Steady crankRate(crank, 1); // acc=0; vel=constant, same for each axis; pos integrated
    // or ...
    //crank.lock(state);
    //crank.setMotionType(state, Regular|Discrete|Fast|Prescribed, stage);

    // need a way to register a mobilizer with a particular relaxation solver,
    // switch between dynamic, continuous relaxed, end-of-step relaxed, discrete.
    // what about a "local" (explicit) relaxation, like q=(q1+q2)/2 ?


    State s = mbs.realizeTopology(); // returns a reference to the the default state
    mbs.realizeModel(s); // define appropriate states for this System

    //crankRate.setRate(s, 3);
    crank.setAngle(s, 5); //q 
    crank.setRate(s, 3);  //u

    Visualizer display(mbs);

    mbs.realize(s, Stage::Position);
    display.report(s);
    cout << "q=" << s.getQ() << endl;
    cout << "qErr=" << s.getQErr() << endl;

    crank.setAngle(s, -30*Deg2Rad);
    //crank.setQToFitRotation (s, Rotation::aboutZ(-.9*Pi/2));


    //TODO
    //rightPendulum.setUToFitLinearVelocity(s, Vec3(1.1,0,1.2));

    //crank.setUToFitAngularVelocity(s, 10*Vec3(.1,.2,.3));

    mbs.realize(s, Stage::Velocity);
    display.report(s);

    cout << "q=" << s.getQ() << endl;
    cout << "qErr=" << s.getQErr() << endl;


    // These are the SimTK Simmath integrators:
    RungeKuttaMersonIntegrator myStudy(mbs);
    //CPodesIntegrator myStudy(mbs, CPodes::BDF, CPodes::Newton);


    //myStudy.setMaximumStepSize(0.001);
    myStudy.setAccuracy(1e-3);
    //myStudy.setProjectEveryStep(true);
    //myStudy.setAllowInterpolation(false);
    //myStudy.setMaximumStepSize(.1);

    const Real dt = 1./30; // output intervals
    const Real finalTime = 10;

    myStudy.setFinalTime(finalTime);

    cout << "Hit ENTER in console to continue ...\n";
    getchar();
    display.report(s);

    cout << "Hit ENTER in console to continue ...\n";
    getchar();

    // Peforms assembly if constraints are violated.
    myStudy.initialize(s);
    myStudy.setProjectEveryStep(false);
    myStudy.setConstraintTolerance(.001);
    myStudy.initialize(myStudy.getState());

    cout << "Using Integrator " << std::string(myStudy.getMethodName()) << ":\n";
    cout << "ACCURACY IN USE=" << myStudy.getAccuracyInUse() << endl;
    cout << "CTOL IN USE=" << myStudy.getConstraintToleranceInUse() << endl;
    cout << "TIMESCALE=" << myStudy.getTimeScaleInUse() << endl;
    cout << "Y WEIGHTS=" << myStudy.getStateWeightsInUse() << endl;
    cout << "1/CTOLS=" << myStudy.getConstraintWeightsInUse() << endl;

    {
        const State& s = myStudy.getState();
        display.report(s);
        cout << "q=" << s.getQ() << endl;
        cout << "qErr=" << s.getQErr() << endl;
    }

    Integrator::SuccessfulStepStatus status;
    int nextReport = 0;
    while ((status=myStudy.stepTo(nextReport*dt, Infinity))
           != Integrator::EndOfSimulation) 
    {
        const State& s = myStudy.getState();
        mbs.realize(s);
        const Real crankAngle = crank.getBodyRotation(s).convertRotationToAngleAxis()[0] * Rad2Deg;
        printf("%5g %10.4g E=%10.8g h%3d=%g %s%s\n", s.getTime(), 
            crankAngle,
            mbs.calcEnergy(s), myStudy.getNumStepsTaken(),
            myStudy.getPreviousStepSizeTaken(),
            Integrator::successfulStepStatusString(status).c_str(),
            myStudy.isStateInterpolated()?" (INTERP)":"");
        printf("     qerr=%10.8g uerr=%10.8g uderr=%10.8g\n",
            crankRocker.getQErr(s).normRMS(),
            crankRocker.getUErr(s).normRMS(),
            crankRocker.getUDotErr(s).normRMS());

        display.report(s);


        if (status == Integrator::ReachedReportTime)
            ++nextReport;
    }

    for (int i=0; i<100; ++i)
        display.report(myStudy.getState());

    printf("Using Integrator %s:\n", myStudy.getMethodName());
    printf("# STEPS/ATTEMPTS = %d/%d\n", myStudy.getNumStepsTaken(), myStudy.getNumStepsAttempted());
    printf("# ERR TEST FAILS = %d\n", myStudy.getNumErrorTestFailures());
    printf("# REALIZE/PROJECT = %d/%d\n", myStudy.getNumRealizations(), myStudy.getNumProjections());

  } 
  catch (const exception& e) {
    printf("EXCEPTION THROWN: %s\n", e.what());
    exit(1);
  }
  catch (...) {
    printf("UNKNOWN EXCEPTION THROWN\n");
    exit(1);
  }

}

