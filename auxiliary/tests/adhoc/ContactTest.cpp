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
 * Test Hunt & Crossley contact.
 */

#include "SimTKsimbody.h"
#include "SimTKsimbody_aux.h" // requires VTK

//#define USE_VTK
#ifdef USE_VTK
#define Visualizer VTKVisualizer
#endif

#include <string>
#include <iostream>
#include <exception>
#include <cmath>
#include <ctime>
using std::cout;
using std::endl;

using namespace SimTK;

static const Real RadiansPerDegree = Pi/180;

// material properties
// Steel
static const Real steel_density = 8000.;  // kg/m^3
static const Real steel_young   = 200e9;  // pascals (N/m)
static const Real steel_poisson = 0.3;    // ratio
static const Real steel_planestrain = 
                    steel_young/(1.-steel_poisson*steel_poisson);
static const Real steel_dissipation = 0.001;

// Concrete
static const Real concrete_density = 2300.;  // kg/m^3
static const Real concrete_young   = 25e9;  // pascals (N/m)
static const Real concrete_poisson = 0.15;    // ratio
static const Real concrete_planestrain = 
                    concrete_young/(1.-concrete_poisson*concrete_poisson);
static const Real concrete_dissipation = 0.005;

// Nylon
static const Real nylon_density = 1100.;  // kg/m^3
static const Real nylon_young   = 2.5e9;  // pascals (N/m)
static const Real nylon_poisson = 0.4;    // ratio
static const Real nylon_planestrain = 
                    nylon_young/(1.-nylon_poisson*nylon_poisson);
static const Real nylon_dissipation = 0.005;

// Rubber
static const Real rubber_density = 1100.;  // kg/m^3
static const Real rubber_young   = 0.01e9; // pascals (N/m)
static const Real rubber_poisson = 0.5;    // ratio
static const Real rubber_planestrain = 
                    rubber_young/(1.-rubber_poisson*rubber_poisson);
static const Real rubber_dissipation = 0.005;

extern "C" void SimTK_version_SimTKlapack(int*,int*,int*);
extern "C" void SimTK_about_SimTKlapack(const char*, int, char*);

static const int NColors =16;
static Vec3 colors[NColors] = {
    Red,Green,Blue,Yellow,Orange,Magenta,Cyan,Purple,
    Red/1.5,Green/1.5,Blue/1.5,Yellow/1.5,Orange/1.5,Magenta/1.5,Cyan/1.5,Purple/1.5
};

static void printFinalStats(const CPodes&);

int main() {

    int major,minor,build;
    char out[100];
    const char* keylist[] = { "version", "library", "type", "debug", "authors", "copyright", "svn_revision", 0 };

    SimTK_version_SimTKlapack(&major,&minor,&build);
    std::printf("==> SimTKlapack library version: %d.%d.%d\n", major, minor, build);
    std::printf("    SimTK_about_SimTKlapack():\n");
    for (const char** p = keylist; *p; ++p) {
        SimTK_about_SimTKlapack(*p, 100, out);
        std::printf("      about(%s)='%s'\n", *p, out);
    }

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
    Vec3 gravity(0.,-g,0.);

    MultibodySystem mbs;
    SimbodyMatterSubsystem  bouncers(mbs);
    HuntCrossleyContact     contact(mbs);
    DecorationSubsystem     artwork(mbs);
    GeneralForceSubsystem   forces(mbs);
   // Force::UniformGravity gravityForces(forces, bouncers, gravity);
    //Force::Gravity gravityForces(forces, bouncers, gravity);
    Force::Gravity gravityForces(forces, bouncers, -UnitVec3(YAxis), g);

    //Force::Thermostat thermostat(forces, bouncers, 6000/*boltzmann??*/,
    //    500, 1);
    //thermostat.setDefaultNumChains(3);
    //thermostat.setDisabledByDefault(true);

    // No, thank you.
    bouncers.setShowDefaultGeometry(false);

    // start with a 2body pendulum with big rubber balls
    Real linkLength = 20.; // m
    Real pendBallRadius = 3; // m
    Real pendMass = rubber_density*(4./3.)*Pi*std::pow(pendBallRadius,3);
    Vec3 pendGroundPt1 = Vec3(-10,50,-10);
    Vec3 pendGroundPt2 = Vec3(20,20,-10);
    Inertia pendBallInertia(pendMass*Gyration::sphere(pendBallRadius));
    const MassProperties pendMProps(pendMass, Vec3(0, -linkLength/2, 0), 
        pendBallInertia.shiftFromMassCenter(Vec3(0, -linkLength/2, 0), pendMass));

    MobilizedBody::Ball pend1(bouncers.Ground(), Transform(pendGroundPt1),
                              Body::Rigid(pendMProps), Transform(Vec3(0, linkLength/2, 0)));

    MobilizedBody::Ball pend2(pend1, Transform(Vec3(0,-linkLength/2,0)),
                              Body::Rigid(pendMProps), Transform(Vec3(0, linkLength/2, 0)));

    Constraint::Rod theConstraint(bouncers.Ground(), pendGroundPt2,
                                  pend2, Vec3(0, -linkLength/2, 0),
                                  20);


    const Real hardBallRadius = 1.25, rubberBallRadius = 1.5;  // m
    const Real hardBallMass = steel_density*(4./3.)*Pi*std::pow(hardBallRadius,3);
    const Real rubberBallMass = rubber_density*(4./3.)*Pi*std::pow(rubberBallRadius,3);

    const MassProperties hardBallMProps(hardBallMass, Vec3(0), 
        hardBallMass*Gyration::sphere(hardBallRadius));
    const MassProperties rubberBallMProps(rubberBallMass, Vec3(0), 
        rubberBallMass*Gyration::sphere(rubberBallRadius));
    const Vec3 firstHardBallPos = Vec3(-6,30,0), firstRubberBallPos = Vec3(18,30,-18);
    std::vector<MobilizedBodyIndex> balls;

    const int NRubberBalls = 6;
    for (int i=0; i<NRubberBalls; ++i)
        balls.push_back(
            MobilizedBody::Cartesian(bouncers.Ground(), Transform(firstRubberBallPos+i*Vec3(0,2*rubberBallRadius+1,0)),
                                     Body::Rigid(rubberBallMProps), Transform()));

    const int NHardBalls = 12;
    for (int i=0; i < NHardBalls; ++i)
        balls.push_back(
            MobilizedBody::Cartesian(bouncers.Ground(),
                Transform(firstHardBallPos+i*Vec3(0,2*hardBallRadius+1,0)
                          + (i==NHardBalls-1)*Vec3(1e-14,0,1e-16)),
                Body::Rigid(hardBallMProps), Transform()));

    // The k's here are the plane-strain moduli, that is, Y/(1-p^2) where Y is
    // Young's modulus and p is Poisson's ratio for the material. The c's are
    // the dissipation coefficients in units of 1/velocity.

    Real kwall=concrete_planestrain, khard=steel_planestrain, krubber=rubber_planestrain;
    Real cwall=concrete_dissipation, chard=steel_dissipation, crubber=rubber_dissipation;
    //const Real cwall = 0., chard = 0., crubber=0.;
    /*XXX*/kwall*=.001, khard*=.001;

    contact.addSphere(pend1, Vec3(0, -linkLength/2, 0), pendBallRadius, krubber, crubber);
    contact.addSphere(pend2, Vec3(0, -linkLength/2, 0), pendBallRadius, krubber, crubber);
    artwork.addBodyFixedDecoration(pend1, Transform(Vec3(0,-linkLength/2,0)), 
        DecorativeSphere(pendBallRadius).setColor(Gray).setOpacity(1));
    artwork.addBodyFixedDecoration(pend1, Transform(), DecorativeLine(Vec3(0,linkLength/2,0),Vec3(0,-linkLength/2,0)));
    artwork.addBodyFixedDecoration(pend2, Transform(Vec3(0,-linkLength/2,0)), 
        DecorativeSphere(pendBallRadius).setColor(Gray).setOpacity(1));
    artwork.addBodyFixedDecoration(pend2, Transform(), DecorativeLine(Vec3(0,linkLength/2,0),Vec3(0,-linkLength/2,0)));
 
    DecorativeLine rbProto; rbProto.setColor(Orange).setLineThickness(1);
    artwork.addRubberBandLine(GroundIndex, pendGroundPt2,pend2,Vec3(0,-linkLength/2,0), rbProto);

    contact.addHalfSpace(GroundIndex, UnitVec3(0,1,0), 0, kwall, cwall);
    contact.addHalfSpace(GroundIndex, UnitVec3(1,0,0), -20, kwall, cwall); // left
    contact.addHalfSpace(GroundIndex, UnitVec3(-1,0,0), -20, kwall, cwall); // right
    contact.addHalfSpace(GroundIndex, UnitVec3(0,0,-1), -20, kwall, cwall); // front
    contact.addHalfSpace(GroundIndex, UnitVec3(0,0,1), -20, kwall, cwall); // back

    artwork.addBodyFixedDecoration(GroundIndex, Transform(), 
        DecorativeBrick(Vec3(20,.1,20)).setColor(Green).setOpacity(1));
    artwork.addBodyFixedDecoration(GroundIndex, Transform(Vec3(-20,20,0)), 
        DecorativeBrick(Vec3(.1,30,20)).setColor(Yellow).setOpacity(.2));
    artwork.addBodyFixedDecoration(GroundIndex, Transform(Vec3(20,20,0)), 
        DecorativeBrick(Vec3(.1,30,20)).setColor(Yellow).setOpacity(.2));
    artwork.addBodyFixedDecoration(GroundIndex, Transform(Vec3(0,20,20)), 
        DecorativeBrick(Vec3(20,20,.1)).setColor(Gray).setOpacity(.05));
    artwork.addBodyFixedDecoration(GroundIndex, Transform(Vec3(0,20,-20)), 
        DecorativeBrick(Vec3(20,30,.1)).setColor(Cyan).setOpacity(.2));

    DecorativeSphere rubberSphere(rubberBallRadius);
    for (int i=0; i<NRubberBalls; ++i) {
        contact.addSphere(balls[i], Vec3(0), rubberBallRadius, krubber, crubber);
        artwork.addBodyFixedDecoration(balls[i], Transform(), 
                    rubberSphere.setColor((Real(i+1)/NRubberBalls)*Gray));
    }

    DecorativeSphere hardSphere(hardBallRadius);
    for (int i=0; i<NHardBalls; ++i) {
        contact.addSphere(balls[NRubberBalls+i], Vec3(0), hardBallRadius, khard, chard);
        artwork.addBodyFixedDecoration(balls[NRubberBalls+i], Transform(), 
                    hardSphere.setColor(colors[i % NColors]));

    }

    State s = mbs.realizeTopology();

    const Real FrameRate = 24;
    Visualizer vtk(mbs);
    //vtk.setCameraClippingPlanes(.5,200);
   // vtk.zoomCameraToShowAllGeometry();

    //vtk.setMode(Visualizer::Sampling);
    vtk.setDesiredFrameRate(FrameRate);
    //vtk.setMode(Visualizer::PassThrough);
    vtk.setMode(Visualizer::RealTime);
    vtk.setDesiredBufferLengthInSec(.15);
    vtk.setDesiredFrameRate(FrameRate);


    //bouncers.setUseEulerAngles(s, true);
    mbs.realizeModel(s);
    bool suppressProjection = false;

    RungeKuttaMersonIntegrator ee(mbs);
    //VerletIntegrator ee(mbs);
    //ee.setMaximumStepSize(.01);
    //CPodesIntegrator ee(mbs, CPodes::BDF, CPodes::Newton);
    //CPodesIntegrator ee(mbs, CPodes::BDF, CPodes::Functional);
    //ee.setOrderLimit(3); //CPodes only
    //ee.setProjectEveryStep(true);

    //ee.setAccuracy(1e-5);
    ee.setAccuracy(2e-2);
    //ee.setConstraintTolerance(1e-3);

    vtk.report(s);

    // set Modeling stuff (s)
    //   none
    //bouncers.setUseEulerAngles(s,true);

    std::vector<State> saveEm;
    saveEm.reserve(10000);

    const Real h = 1/FrameRate; // output every frame
    const Real tstart = 0.;
    const Real tmax = 30;


    s.updTime() = tstart;
    mbs.realize(s, Stage::Acceleration);
    for (int i=0; i<25; ++i)
        saveEm.push_back(s);    // delay
    vtk.report(s);

    ee.initialize(s);
    for (int i=0; i<25; ++i)
        saveEm.push_back(ee.getState());    // delay
    vtk.report(ee.getState());

    cout << "Using Integrator " << std::string(ee.getMethodName()) << ":\n";
    cout << "ACCURACY IN USE=" << ee.getAccuracyInUse() << endl;
    cout << "CTOL IN USE=" << ee.getConstraintToleranceInUse() << endl;
    cout << "TIMESCALE=" << ee.getTimeScaleInUse() << endl;
    cout << "Y WEIGHTS=" << ee.getStateWeightsInUse() << endl;
    cout << "1/CTOLS=" << ee.getConstraintWeightsInUse() << endl;

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
             << ") qerr=" << bouncers.getQErr(ss).normRMS()
             << " uerr=" << bouncers.getUErr(ss).normRMS()
             << " hNext=" << ee.getPredictedNextStepSize() << endl;
        }
        ++step;

        ee.stepTo(ss.getTime() + h);
        vtk.report(ss);
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


    //printFinalStats(ee.getCPodes());
    printf("Using Integrator %s:\n", ee.getMethodName());
    printf("# STEPS/ATTEMPTS = %d/%d\n", ee.getNumStepsTaken(), ee.getNumStepsAttempted());
    printf("# ERR TEST FAILS = %d\n", ee.getNumErrorTestFailures());
    printf("# REALIZE/PROJECT = %d/%d\n", ee.getNumRealizations(), ee.getNumProjections());

    vtk.dumpStats(std::cout);

    while(true) {
        for (int i=0; i < (int)saveEm.size(); ++i) {
            for (int j=0; j < 4; ++j) // 1/4 speed
                vtk.report(saveEm[i]);
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

static void printFinalStats(const CPodes& cpodes)
{
  Real h0u;
  int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  int nproj, nce, nsetupsP, nprf;
  int flag;

  h0u=NaN;
  nst=nfe=nsetups=nje=nfeLS=nni=ncfn=netf=nge=-1;
  nproj=nce=nsetupsP=nprf=-1;

  CPodes& cpode = const_cast<CPodes&>(cpodes);
  flag = cpode.getActualInitStep(&h0u);
  flag = cpode.getNumSteps(&nst);
  flag = cpode.getNumFctEvals(&nfe);
  flag = cpode.getNumLinSolvSetups(&nsetups);
  flag = cpode.getNumErrTestFails(&netf);
  flag = cpode.getNumNonlinSolvIters(&nni);
  flag = cpode.getNumNonlinSolvConvFails(&ncfn);
  flag = cpode.dlsGetNumJacEvals(&nje);
  flag = cpode.dlsGetNumFctEvals(&nfeLS);
  flag = cpode.getProjStats(&nproj, &nce, &nsetupsP, &nprf);
  flag = cpode.getNumGEvals(&nge);

 // h0u   = integ.getActualInitialStepSizeTaken();
 // nst   = integ.getNumStepsTaken();
 // nfe   = integ.getNumRealizations();
 // netf  = integ.getNumErrorTestFailures();
 // nproj = integ.getNumProjections();
 // nprf  = integ.getNumProjectionFailures();

  printf("\nFinal Statistics:\n");
  printf("h0u = %g\n",h0u);
  printf("nst = %-6d nfe  = %-6d nsetups = %-6d\n",
     nst, nfe, nsetups);
  printf("nfeLS = %-6d nje = %d\n",
     nfeLS, nje);
  printf("nni = %-6d ncfn = %-6d netf = %-6d \n",
     nni, ncfn, netf);
  printf("nproj = %-6d nce = %-6d nsetupsP = %-6d nprf = %-6d\n",
         nproj, nce, nsetupsP, nprf);
  printf("nge = %d\n", nge);

}
