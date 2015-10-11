/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): Free Water Test                        *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2009-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

/**@file
 * This is an outer block for simulating ??? in various ways with Simbody.
 * This is about testing Simbody, *not* studying ???!
 */

#include "SimTKsimbody.h"

#include <string>
#include <iostream>
#include <exception>

using std::cout;
using std::cin;
using std::endl;

using namespace SimTK;

static const Real Deg2Rad = (Real)SimTK_DEGREE_TO_RADIAN,
                  Rad2Deg = (Real)SimTK_RADIAN_TO_DEGREE;

static Real g = 9.8;
static Real m = 1;

int main(int argc, char** argv) {
    static const Transform GroundFrame;
    static const Rotation ZalongY(UnitVec3(XAxis), XAxis, UnitVec3(YAxis), ZAxis);
    static const Rotation ZalongX(UnitVec3(XAxis), ZAxis, UnitVec3(YAxis), YAxis);
    static const Rotation YalongZ(UnitVec3(ZAxis), XAxis, UnitVec3(XAxis), YAxis);
    static const Vec3 TestLoc(1,0,0);

  try { // If anything goes wrong, an exception will be thrown.

        // CREATE MULTIBODY SYSTEM AND ITS SUBSYSTEMS
    MultibodySystem             mbs;

    SimbodyMatterSubsystem      matter(mbs);
    GeneralForceSubsystem       forces(mbs);
    DecorationSubsystem         viz(mbs);
    //Force::UniformGravity       gravity(forces, matter, Vec3(0, -g, 0));

        // ADD BODIES AND THEIR MOBILIZERS
    Body::Rigid oxygen = Body::Rigid(MassProperties(m, Vec3(0), Inertia(0)));
    oxygen.addDecoration(DecorativeSphere(.1).setColor(Red).setOpacity(.3));
    Body::Rigid hydrogen = Body::Rigid(MassProperties(m, Vec3(0), Inertia(0)));
    hydrogen.addDecoration(DecorativeSphere(.05).setColor(Green).setOpacity(.3));

    MobilizedBody::Cartesian
        masslessFrame(matter.Ground(), Transform(ZalongY, Vec3(0)),
                      MassProperties(0,Vec3(0),Inertia(0)), Transform());

    MobilizedBody::Pin
        H1(masslessFrame, Transform(), hydrogen, Transform());

    MobilizedBody::Pin
        O(H1, ZalongY, oxygen, Transform(Vec3(0,1,0)));

    MobilizedBody::Universal
        H2(O, Transform(YalongZ), hydrogen, Transform(Vec3(0,0,1)));

    Force::MobilityLinearSpring(forces, H2, 1, 2, 60*Deg2Rad); // harmonic bend


    //MobilizedBody::SphericalCoords
    //    H2(matter.Ground(), Transform(ZUp, TestLoc),
    //           particle, Transform(),
    //           0*Deg2Rad,  false,   // azimuth offset, negated
    //           0,          false,   // zenith offset, negated
    //           ZAxis,      true);  // translation axis, negated

    //anAtom.setDefaultRadius(.1);
    //anAtom.setDefaultAngles(Vec2(0, 30*Deg2Rad));

    //viz.addRubberBandLine(matter.Ground(), TestLoc,
    //                      anAtom, Vec3(0),
    //                      DecorativeLine().setColor(Orange).setLineThickness(4));

    //Force::MobilityLinearSpring(forces, anAtom, 1, 2, 45*Deg2Rad); // harmonic bend
    //Force::MobilityLinearSpring(forces, anAtom, 2, 20, .5); // harmonic stretch

    //Force::MobilityLinearDamper(forces, anAtom, 0, .1); // harmonic bend
    //Force::MobilityLinearDamper(forces, anAtom, 1, .1); // harmonic bend
    //Force::MobilityLinearDamper(forces, anAtom, 2, .1); // harmonic stretch


    State s = mbs.realizeTopology(); // returns a reference to the the default state
    mbs.realizeModel(s); // define appropriate states for this System
    mbs.realize(s, Stage::Instance); // instantiate constraints if any

    Visualizer display(mbs);
    display.setBackgroundColor(White);
    display.setBackgroundType(Visualizer::SolidColor);

    mbs.realize(s, Stage::Velocity);
    display.report(s);

    cout << "q=" << s.getQ() << endl;
    cout << "u=" << s.getU() << endl;

    char c;
    //cout << "Default configuration shown. Ready?\n"; cin >> c;

    O.setOneQ(s, 0, 30*Deg2Rad);
    //H2.setOneQ(s, 0, 30*Deg2Rad);
    H2.setOneQ(s, 0, 30*Deg2Rad);
    H2.setOneQ(s, 1, 50*Deg2Rad);

    s.updU() = Vector(Vec7(0,0,0,1,1,1,0));

    //anAtom.setQToFitRotation(s, Rotation(-.9*Pi/2,YAxis));

    //while (true) {
    //    Real x;
    //    cout << "Torsion (deg)? "; cin >> x; if (x==1234) break;
    //    Vec2 a = anAtom.getAngles(s); a[0]=x*Deg2Rad; anAtom.setAngles(s, a);
    //    display.report(s);
    //    cout << "Bend (deg)? "; cin >> x; if (x==1234) break;
    //    a = anAtom.getAngles(s); a[1]=x*Deg2Rad; anAtom.setAngles(s, a);
    //    display.report(s);
    //    cout << "Radius? "; cin >> x; if (x==1234) break;
    //    anAtom.setRadius(s, x);
    //    display.report(s);
    //}
    //anAtom.setUToFitAngularVelocity(s, Vec3(.1,.2,.3));

    //anAtom.setAngle(s, 45*Deg2Rad);
    //anAtom.setTranslation(s, Vec2(.4, .1));

    mbs.realize(s, Stage::Dynamics);
    mbs.realize(s, Stage::Acceleration);

    cout << "q=" << s.getQ() << endl;
    cout << "u=" << s.getU() << endl;
    cout << "qdot=" << s.getQDot() << endl;
    cout << "udot=" << s.getUDot() << endl;
    cout << "qdotdot=" << s.getQDotDot() << endl;
    display.report(s);

    cout << "Initialized configuration shown. Ready? ";
    cin >> c;

    RungeKuttaMersonIntegrator myStudy(mbs);
    myStudy.setAccuracy(1e-4);

    const Real dt = .02; // output intervals
    const Real finalTime = 20;

    myStudy.setFinalTime(finalTime);

    // Peforms assembly if constraints are violated.
    myStudy.initialize(s);

    cout << "Using Integrator " << std::string(myStudy.getMethodName()) << ":\n";
    cout << "ACCURACY IN USE=" << myStudy.getAccuracyInUse() << endl;
    cout << "CTOL IN USE=" << myStudy.getConstraintToleranceInUse() << endl;
    cout << "TIMESCALE=" << mbs.getDefaultTimeScale() << endl;
    cout << "U WEIGHTS=" << s.getUWeights() << endl;
    cout << "Z WEIGHTS=" << s.getZWeights() << endl;
    cout << "1/QTOLS=" << s.getQErrWeights() << endl;
    cout << "1/UTOLS=" << s.getUErrWeights() << endl;

    Integrator::SuccessfulStepStatus status;
    int nextReport = 0;
    while ((status=myStudy.stepTo(nextReport*dt))
           != Integrator::EndOfSimulation)
    {
        const State& s = myStudy.getState();
        mbs.realize(s);
        printf("%5g %10.4g %10.4g %10.4g %10.4g E=%10.8g h%3d=%g %s%s\n", s.getTime(),
            H1.getAngle(s), O.getAngle(s), H2.getOneQ(s,0), H2.getOneQ(s,1),
            mbs.calcEnergy(s), myStudy.getNumStepsTaken(),
            myStudy.getPreviousStepSizeTaken(),
            Integrator::getSuccessfulStepStatusString(status).c_str(),
            myStudy.isStateInterpolated()?" (INTERP)":"");

        display.report(s);

        if (status == Integrator::ReachedReportTime)
            ++nextReport;
    }

    printf("Using Integrator %s:\n", myStudy.getMethodName());
    printf("# STEPS/ATTEMPTS = %d/%d\n", myStudy.getNumStepsTaken(), myStudy.getNumStepsAttempted());
    printf("# ERR TEST FAILS = %d\n", myStudy.getNumErrorTestFailures());
    printf("# REALIZE/PROJECT = %d/%d\n", myStudy.getNumRealizations(), myStudy.getNumProjections());

  }
  catch (const std::exception& e) {
    printf("EXCEPTION THROWN: %s\n", e.what());
    exit(1);
  }
  catch (...) {
    printf("UNKNOWN EXCEPTION THROWN\n");
    exit(1);
  }

}
