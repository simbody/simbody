/* -------------------------------------------------------------------------- *
 *                           SimTK Simbody(tm)                                *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-13 Stanford University and the Authors.        *
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

/* Test compliant contact. Attempts to run at 3X real time. */

#include "Simbody.h"

#include <string>
#include <iostream>
#include <exception>
#include <cmath>
#include <ctime>
using std::cout;
using std::endl;

using namespace SimTK;

#define ANIMATE // off to get more accurate CPU time (you can still playback)

namespace {

const Real mu_s = 0*.2;
const Real mu_d = 0*.1;
const Real mu_v = 0.01;

// material properties
// Steel
const Real steel_density = 8000.;  // kg/m^3
const Real steel_young   = 200e9;  // pascals (N/m)
const Real steel_poisson = 0.3;    // ratio
const Real steel_planestrain = 
    ContactMaterial::calcPlaneStrainStiffness(steel_young,steel_poisson);
const Real steel_dissipation = 0.001;

const ContactMaterial steel(steel_planestrain,steel_dissipation,0,0,0);

// Concrete
const Real concrete_density = 2300.;  // kg/m^3
const Real concrete_young   = 25e9;  // pascals (N/m)
const Real concrete_poisson = 0.15;    // ratio
const Real concrete_planestrain = 
    ContactMaterial::calcPlaneStrainStiffness(concrete_young,concrete_poisson);
const Real concrete_dissipation = 0.005;

const ContactMaterial concrete(concrete_planestrain,concrete_dissipation,
                               mu_s,mu_d,mu_v);

// Nylon
const Real nylon_density = 1100.;  // kg/m^3
const Real nylon_young   = 2.5e9;  // pascals (N/m)
const Real nylon_poisson = 0.4;    // ratio
const Real nylon_planestrain =
    ContactMaterial::calcPlaneStrainStiffness(nylon_young,nylon_poisson);
const Real nylon_dissipation = 0.005;

const ContactMaterial nylon(nylon_planestrain,nylon_dissipation,
                               mu_s,mu_d,mu_v);

// Rubber
const Real rubber_density = 1100.;  // kg/m^3
const Real rubber_young   = 0.01e9; // pascals (N/m)
const Real rubber_poisson = 0.5;    // ratio
const Real rubber_planestrain = 
    ContactMaterial::calcPlaneStrainStiffness(rubber_young,rubber_poisson);
const Real rubber_dissipation = 0.005;

const ContactMaterial rubber(rubber_planestrain,rubber_dissipation,
                               mu_s,mu_d,mu_v);

extern "C" void SimTK_version_SimTKlapack(int*,int*,int*);
extern "C" void SimTK_about_SimTKlapack(const char*, int, char*);

const int NColors =16;
Vec3 colors[NColors] = {
    Red,Green,Blue,Yellow,Orange,Magenta,Cyan,Purple,
    Red/1.5,Green/1.5,Blue/1.5,Yellow/1.5,Orange/1.5,Magenta/1.5,Cyan/1.5,Purple/1.5
};

const ContactMaterial wallMaterial = nylon;
const ContactMaterial softMaterial = rubber;
const ContactMaterial hardMaterial = nylon;

const int  NRubberBalls = 20, NHardBalls = 40;
const Real PendBallRadius = 3, RubberBallRadius = 2.5, HardBallRadius = 2;//m

const Real PendBallMass = rubber_density*(4./3.)*Pi*cube(PendBallRadius);
const Real RubberBallMass = rubber_density*(4./3.)*Pi*cube(RubberBallRadius);
const Real HardBallMass = steel_density*(4./3.)*Pi*cube(HardBallRadius);


// Write interesting integrator info to stdout.
void dumpIntegratorStats(const Integrator& integ);
}


//==============================================================================
//                                   MAIN
//==============================================================================
int main() {

    int major,minor,build;
    char out[100];
    const char* keylist[] = { "version", "library", "type", "debug", "authors", 
                              "copyright", "svn_revision", 0 };

    SimTK_version_SimTKcommon(&major,&minor,&build);
    std::printf("==> SimTKcommon library version: %d.%d.%d\n", major, minor, build);
    std::printf("    SimTK_about_SimTKcommon():\n");
    for (const char** p = keylist; *p; ++p) {
        SimTK_about_SimTKcommon(*p, 100, out);
        std::printf("      about(%s)='%s'\n", *p, out);
    }

    SimTK_version_simmath(&major,&minor,&build);
    std::printf("==> SimTKmath library version: %d.%d.%d\n", major, minor, build);
    std::printf("    SimTK_about_simmath():\n");
    for (const char** p = keylist; *p; ++p) {
        SimTK_about_simmath(*p, 100, out);
        std::printf("      about(%s)='%s'\n", *p, out);
    }

    SimTK_version_simbody(&major,&minor,&build);
    std::printf("==> simbody library version: %d.%d.%d\n", major, minor, build);
    std::printf("    SimTK_about_simbody():\n");
    for (const char** p = keylist; *p; ++p) {
        SimTK_about_simbody(*p, 100, out);
        std::printf("      about(%s)='%s'\n", *p, out);
    }


try
  { Real g = 9.8;   // m/s^2
    //Real g = 1.6;   // m/s^2 (moon)

    MultibodySystem         mbs;
    SimbodyMatterSubsystem  matter(mbs);

    GeneralForceSubsystem   forces(mbs);
    Force::Gravity gravityForces(forces, matter, -YAxis, g);

    ContactTrackerSubsystem  tracker(mbs);
    CompliantContactSubsystem contactForces(mbs, tracker);
    //contactForces.setTransitionVelocity(1e-3);
    contactForces.setTransitionVelocity(1e-2);

    //Force::Thermostat thermostat(forces, matter, 6000/*boltzmann??*/, 500, 1);
    //Force::Thermostat thermostat(forces, matter, 6000/*boltzmann??*/, 1000, 1, 0);
    //thermostat.setDefaultNumChains(3);
    //thermostat.setDisabledByDefault(true);

    // No, thank you.
    matter.setShowDefaultGeometry(false);

    MobilizedBody& Ground = matter.Ground();

    // Add the Ground contact geometry. Contact half spaces have -XAxis normals
    // (right hand wall) so we'll have to rotate them.
    const Rotation R_right; // identity
    const Rotation R_back(Pi/2,YAxis);
    const Rotation R_left(Pi,YAxis);
    const Rotation R_front(-Pi/2,YAxis);
    const Rotation R_ceiling(Pi/2,ZAxis);
    const Rotation R_floor(-Pi/2,ZAxis);


    Ground.updBody().addContactSurface(Transform(R_floor,Vec3(0)),
        ContactSurface(ContactGeometry::HalfSpace(),wallMaterial));
    Ground.updBody().addContactSurface(Transform(R_left,Vec3(-20,0,0)),
        ContactSurface(ContactGeometry::HalfSpace(),wallMaterial));
    Ground.updBody().addContactSurface(Transform(R_right,Vec3(20,0,0)),
        ContactSurface(ContactGeometry::HalfSpace(),wallMaterial));
    Ground.updBody().addContactSurface(Transform(R_front,Vec3(0,0,20)),
        ContactSurface(ContactGeometry::HalfSpace(),wallMaterial));
    Ground.updBody().addContactSurface(Transform(R_back,Vec3(0,0,-20)),
        ContactSurface(ContactGeometry::HalfSpace(),wallMaterial));

    DecorativeBrick smallWall(Vec3(.1,20,20)), tallWall(Vec3(.1,30,20));
    Ground.addBodyDecoration(Transform(R_floor,Vec3(0)),
        smallWall.setColor(Green).setOpacity(1));
    Ground.addBodyDecoration(Transform(R_left,Vec3(-20,20,0)), 
        tallWall.setColor(Yellow).setOpacity(.2));
    Ground.addBodyDecoration(Transform(R_right,Vec3(20,20,0)), 
        tallWall.setColor(Yellow).setOpacity(.2));
    Ground.addBodyDecoration(Transform(R_front,Vec3(0,20,20)), 
        smallWall.setColor(Gray).setOpacity(.05));
    Ground.addBodyDecoration(Transform(R_back, Vec3(0,20,-20)), 
        tallWall.setColor(Cyan).setOpacity(.2));

    // start with a 2body pendulum with big rubber balls
    Real linkLength = 20.; // m
    Vec3 pendGroundPt1 = Vec3(-10,50,-10);
    Vec3 pendGroundPt2 = Vec3(20,20,-10);
    Inertia pendBallInertia(PendBallMass*UnitInertia::sphere(PendBallRadius));
    const MassProperties pendMProps(PendBallMass, Vec3(0, -linkLength/2, 0), 
        pendBallInertia.shiftFromMassCenter(Vec3(0, -linkLength/2, 0), PendBallMass));

    Body::Rigid pendBodyInfo(pendMProps);
    pendBodyInfo.addContactSurface(Vec3(0, -linkLength/2, 0),
        ContactSurface(ContactGeometry::Sphere(PendBallRadius), softMaterial));
    pendBodyInfo.addDecoration(Vec3(0,-linkLength/2,0),
        DecorativeSphere(PendBallRadius).setColor(Gray).setOpacity(1));
    pendBodyInfo.addDecoration(Vec3(0),
        DecorativeLine(Vec3(0,linkLength/2,0),Vec3(0,-linkLength/2,0)));

    MobilizedBody::Ball pend1(matter.Ground(), Transform(pendGroundPt1),
                              pendBodyInfo, Transform(Vec3(0, linkLength/2, 0)));
    MobilizedBody::Ball pend2(pend1, Transform(Vec3(0,-linkLength/2,0)),
                              pendBodyInfo, Transform(Vec3(0, linkLength/2, 0)));

    // Add a "rod" constraint to make the two-link pendulum a loop.
    Constraint::Rod theConstraint(matter.Ground(), pendGroundPt2,
                                  pend2, Vec3(0, -linkLength/2, 0),
                                  20);


    Body::Rigid hardBallInfo(MassProperties(HardBallMass, Vec3(0), 
                                UnitInertia::sphere(HardBallRadius)));
    hardBallInfo.addContactSurface(Vec3(0),
        ContactSurface(ContactGeometry::Sphere(HardBallRadius),hardMaterial));
        //ContactSurface(ContactGeometry::Ellipsoid(Vec3(4,HardBallRadius,.7)),hardMaterial));
    DecorativeSphere hardSphere(HardBallRadius);
    //DecorativeEllipsoid hardSphere(Vec3(4,HardBallRadius,.7));

    Body::Rigid rubberBallInfo(MassProperties(RubberBallMass, Vec3(0), 
                                  UnitInertia::sphere(RubberBallRadius)));
    rubberBallInfo.addContactSurface(Vec3(0),
        ContactSurface(ContactGeometry::Sphere(RubberBallRadius),softMaterial));
        //ContactSurface(ContactGeometry::Ellipsoid(Vec3(.7,RubberBallRadius,5)),softMaterial));
    DecorativeSphere rubberSphere(RubberBallRadius);
    //DecorativeEllipsoid rubberSphere(Vec3(.7,RubberBallRadius,5));

    const Vec3 firstHardBallPos   = Vec3(-6,30,  0), 
               firstRubberBallPos = Vec3(18,30,-18);

    for (int i=0; i<NRubberBalls; ++i) {
        //MobilizedBody::Cartesian 
        MobilizedBody::Free 
            rubberBall(matter.Ground(),
            Transform(firstRubberBallPos+i*Vec3(0,2*RubberBallRadius+1,0)),
            rubberBallInfo, Transform());
        rubberBall.addBodyDecoration(Vec3(0), 
            rubberSphere.setColor((Real(i+1)/NRubberBalls)*Gray));
    }

    for (int i=0; i < NHardBalls; ++i) {
        //MobilizedBody::Cartesian 
        MobilizedBody::Free 
            hardBall(matter.Ground(),
            Transform(firstHardBallPos+i*Vec3(0,2*HardBallRadius+1,0)
                        + (i==NHardBalls-1)*Vec3(1e-14,0,1e-16)),
            hardBallInfo, Transform());
        hardBall.addBodyDecoration(Vec3(0), 
            hardSphere.setColor(colors[i % NColors]));
    }


    State s = mbs.realizeTopology();

    const Real TargetElapsedTime = 10;
    const Real FrameRate = 30;
    //const Real TimeScale = 1;
    //const Real FrameRate = 60;
    const Real TimeScale = 3; // i.e., 3X realtime
    Visualizer viz(mbs);
    viz.setShowFrameRate(true);
    viz.setShowSimTime(true);

    //viz.setMode(Visualizer::Sampling);
    //viz.setMode(Visualizer::PassThrough);
    viz.setMode(Visualizer::RealTime);
    viz.setDesiredFrameRate(FrameRate);
    viz.setRealTimeScale(TimeScale);

    // Illustrate the rod with a rubber band line.
    DecorativeLine rbProto; rbProto.setColor(Orange).setLineThickness(1);
    viz.addRubberBandLine(GroundIndex, pendGroundPt2,
                          pend2, Vec3(0,-linkLength/2,0), rbProto);

    mbs.realizeModel(s);
    bool suppressProjection = false;

    RungeKuttaMersonIntegrator ee(mbs); const Real Accuracy = .02;

    //RungeKutta3Integrator ee(mbs); const Real Accuracy = .05;

    //This runs the fastest:
    //SemiExplicitEuler2Integrator ee(mbs); const Real Accuracy = .3;
    //ee.setMaximumStepSize(.02); // 20ms

    ee.setAccuracy(Accuracy);
    ee.setConstraintTolerance(std::min(.001, Accuracy/10));


    viz.report(s);

    std::vector<State> saveEm;
    saveEm.reserve(10000);

    const Real h = TimeScale/FrameRate; // output every frame
    const Real tstart = 0.;
    const Real tmax = TimeScale*TargetElapsedTime;


    s.updTime() = tstart;
    mbs.realize(s, Stage::Acceleration);
    for (int i=0; i<25; ++i)
        saveEm.push_back(s);    // delay
    viz.report(s);

    ee.initialize(s);
    for (int i=0; i<25; ++i)
        saveEm.push_back(ee.getState());    // delay
    viz.report(ee.getState());

    cout << "Num mobilities=" << matter.getNumMobilities() << endl;
    cout << "Using Integrator " << std::string(ee.getMethodName()) << ":\n";
    cout << "ACCURACY IN USE=" << ee.getAccuracyInUse() << endl;
    cout << "CTOL IN USE=" << ee.getConstraintToleranceInUse() << endl;
    cout << "TIMESCALE=" << mbs.getDefaultTimeScale() << endl;
    cout << "U WEIGHTS=" << s.getUWeights() << endl;
    cout << "Z WEIGHTS=" << s.getZWeights() << endl;
    cout << "1/QTOLS=" << s.getQErrWeights() << endl;
    cout << "1/UTOLS=" << s.getUErrWeights() << endl;

    const double startCPU  = cpuTime();
    const double startThreadCPU  = threadCpuTime();
    const double startTime = realTime();

    int step = 0;
    while (ee.getTime() <= tmax) {
        const State& ss = ee.getState();
        if (!(step % 30)) {
            mbs.realize(ss);
            cout << ss.getTime() 
                //<< ": T=" << thermostat.getCurrentTemperature(ss)
             << " E=" << mbs.calcEnergy(ss)
             << " (pe=" << mbs.calcPotentialEnergy(ss)
             << ", ke=" << mbs.calcKineticEnergy(ss)
             << ") qerr=" << matter.getQErr(ss).normRMS()
             << " uerr=" << matter.getUErr(ss).normRMS()
             << " hNext=" << ee.getPredictedNextStepSize() << endl;
        }
        ++step;

        ee.stepTo(ss.getTime() + h);
        #ifdef ANIMATE
        viz.report(ss);
        #endif
        saveEm.push_back(ss);

        //if (std::abs(ss.getTime()-30) < .001)
        //    thermostat.setBathTemperature(ee.updAdvancedState(), 10);
        //if (std::abs(ss.getTime()-80) < .001)
        //    thermostat.setBathTemperature(ee.updAdvancedState(), 200);
    }

    std::cout << "Simulated " << ee.getTime() << " seconds in " <<
        realTime()-startTime << " elapsed s\n";
    std::cout << "Process CPU=" << cpuTime()-startCPU << std::endl;
    std::cout << "Thread CPU=" << threadCpuTime()-startThreadCPU << std::endl;


    //printCPodesStats(ee.getCPodes());
    dumpIntegratorStats(ee);
    viz.dumpStats(std::cout);

    while(true) {
        for (int i=0; i < (int)saveEm.size(); ++i) {
            viz.report(saveEm[i]);
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

namespace {
//==============================================================================
//                        DUMP INTEGRATOR STATS
//==============================================================================
void dumpIntegratorStats(const Integrator& integ) {
    const int evals = integ.getNumRealizations();
    std::cout << "\nDone -- simulated " << integ.getTime() << "s with " 
            << integ.getNumStepsTaken() << " steps, avg step=" 
        << (1000*integ.getTime())/integ.getNumStepsTaken() << "ms " 
        << (1000*integ.getTime())/evals << "ms/eval\n";

    printf("Used Integrator %s at accuracy %g:\n", 
        integ.getMethodName(), integ.getAccuracyInUse());
    printf("# STEPS/ATTEMPTS = %d/%d\n",  integ.getNumStepsTaken(), 
                                          integ.getNumStepsAttempted());
    printf("# ERR TEST FAILS = %d\n",     integ.getNumErrorTestFailures());
    printf("# REALIZE/PROJECT = %d/%d\n", integ.getNumRealizations(), 
                                          integ.getNumProjections());
}
}
