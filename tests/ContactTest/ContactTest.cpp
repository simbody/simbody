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
 * Test Hunt & Crossley contact.
 */

#include "SimTKsimbody.h"

#include <string>
#include <iostream>
#include <exception>
#include <cmath>
using std::cout;
using std::endl;

using namespace SimTK;

static const Real Pi = std::acos(-1.), RadiansPerDegree = Pi/180;

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

    SimTK_version_SimTKcpodes(&major,&minor,&build);
    std::printf("==> SimTKcpodes library version: %d.%d.%d\n", major, minor, build);
    std::printf("    SimTK_about_SimTKcpodes():\n");
    for (const char** p = keylist; *p; ++p) {
        SimTK_about_SimTKcpodes(*p, 100, out);
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
    Vec3 gravity(0.,-g,0.);

    SimbodyMatterSubsystem bouncers;

    // start with a 2body pendulum with big rubber balls
    Real linkLength = 20.; // m
    Real pendBallRadius = 3; // m
    Real pendMass = rubber_density*(4./3.)*Pi*std::pow(pendBallRadius,3);
    Vec3 pendGroundPt1 = Vec3(-10,50,-10);
    Vec3 pendGroundPt2 = Vec3(20,20,-10);
    Inertia pendBallInertia(pendMass*Inertia::sphere(pendBallRadius));
    const MassProperties pendMProps(pendMass, Vec3(0, -linkLength/2, 0), 
        pendBallInertia.shiftFromCOM(Vec3(0, -linkLength/2, 0), pendMass));

    BodyId pend1 = bouncers.addRigidBody(pendMProps, Transform(Vec3(0, linkLength/2, 0)),
                          GroundId, Transform(pendGroundPt1),
                          Mobilizer(Mobilizer::Ball, false));

    BodyId pend2 = bouncers.addRigidBody(pendMProps, Transform(Vec3(0, linkLength/2, 0)),
                          pend1, Transform(Vec3(0,-linkLength/2,0)),
                          Mobilizer(Mobilizer::Ball, false));

    int theConstraint =
       bouncers.addConstantDistanceConstraint(GroundId, pendGroundPt2,
                                     pend2, Vec3(0, -linkLength/2, 0),20);


    const Real hardBallRadius = 1.25, rubberBallRadius = 1.5;  // m
    const Real hardBallMass = steel_density*(4./3.)*Pi*std::pow(hardBallRadius,3);
    const Real rubberBallMass = rubber_density*(4./3.)*Pi*std::pow(rubberBallRadius,3);

    const MassProperties hardBallMProps(hardBallMass, Vec3(0), 
        hardBallMass*Inertia::sphere(hardBallRadius));
    const MassProperties rubberBallMProps(rubberBallMass, Vec3(0), 
        rubberBallMass*Inertia::sphere(rubberBallRadius));
    const Vec3 firstHardBallPos = Vec3(-6,30,0), firstRubberBallPos = Vec3(18,30,-18);
    std::vector<BodyId> balls;

    const int NRubberBalls = 6;
    for (int i=0; i<NRubberBalls; ++i)
        balls.push_back( 
            bouncers.addRigidBody(rubberBallMProps, Transform(),
                                  GroundId, Transform(firstRubberBallPos+i*Vec3(0,2*rubberBallRadius+1,0)),
                                  Mobilizer(Mobilizer::Cartesian, false)));

    const int NHardBalls = 12;
    for (int i=0; i < NHardBalls; ++i)
        balls.push_back(
            bouncers.addRigidBody(hardBallMProps, Transform(),
                                  GroundId, Transform(firstHardBallPos+i*Vec3(0,2*hardBallRadius+1,0)
                                                + (i==NHardBalls-1)*Vec3(1e-14,0,1e-16)),
                                  Mobilizer(Mobilizer::Cartesian, false)));


    MultibodySystem mbs;
    mbs.setMatterSubsystem(bouncers);

    UniformGravitySubsystem gravityForces(gravity);
    mbs.addForceSubsystem(gravityForces);

    HuntCrossleyContact contact;
    mbs.addForceSubsystem(contact);

    // The k's here are the plane-strain moduli, that is, Y/(1-p^2) where Y is
    // Young's modulus and p is Poisson's ratio for the material. The c's are
    // the dissipation coefficients in units of 1/velocity.

    const Real kwall=concrete_planestrain, khard=steel_planestrain, krubber=rubber_planestrain;
    const Real cwall=concrete_dissipation, chard=steel_dissipation, crubber=rubber_dissipation;
    //const Real cwall = 0., chard = 0., crubber=0.;

    VTKReporter vtk(mbs, false); // suppress default geometry
    contact.addSphere(pend1, Vec3(0, -linkLength/2, 0), pendBallRadius, krubber, crubber);
    contact.addSphere(pend2, Vec3(0, -linkLength/2, 0), pendBallRadius, krubber, crubber);
    vtk.addDecoration(pend1, Transform(Vec3(0,-linkLength/2,0)), 
        DecorativeSphere(pendBallRadius).setColor(Gray).setOpacity(1));
    vtk.addDecoration(pend1, Transform(), DecorativeLine(Vec3(0,linkLength/2,0),Vec3(0,-linkLength/2,0)));
    vtk.addDecoration(pend2, Transform(Vec3(0,-linkLength/2,0)), 
        DecorativeSphere(pendBallRadius).setColor(Gray).setOpacity(1));
    vtk.addDecoration(pend2, Transform(), DecorativeLine(Vec3(0,linkLength/2,0),Vec3(0,-linkLength/2,0)));
 
    DecorativeLine rbProto; rbProto.setColor(Orange).setLineThickness(1);
    vtk.addRubberBandLine(GroundId, pendGroundPt2,pend2,Vec3(0,-linkLength/2,0), rbProto);

    contact.addHalfSpace(GroundId, UnitVec3(0,1,0), 0, kwall, cwall);
    contact.addHalfSpace(GroundId, UnitVec3(1,0,0), -20, kwall, cwall); // left
    contact.addHalfSpace(GroundId, UnitVec3(-1,0,0), -20, kwall, cwall); // right
    contact.addHalfSpace(GroundId, UnitVec3(0,0,-1), -20, kwall, cwall); // front
    contact.addHalfSpace(GroundId, UnitVec3(0,0,1), -20, kwall, cwall); // back

    vtk.addDecoration(GroundId, Transform(), 
        DecorativeBrick(Vec3(20,.1,20)).setColor(Green).setOpacity(1));
    vtk.addDecoration(GroundId, Transform(Vec3(-20,20,0)), 
        DecorativeBrick(Vec3(.1,30,20)).setColor(Yellow).setOpacity(.2));
    vtk.addDecoration(GroundId, Transform(Vec3(20,20,0)), 
        DecorativeBrick(Vec3(.1,30,20)).setColor(Yellow).setOpacity(.2));
    vtk.addDecoration(GroundId, Transform(Vec3(0,20,20)), 
        DecorativeBrick(Vec3(20,20,.1)).setColor(Gray).setOpacity(.05));
    vtk.addDecoration(GroundId, Transform(Vec3(0,20,-20)), 
        DecorativeBrick(Vec3(20,30,.1)).setColor(Cyan).setOpacity(.2));

    DecorativeSphere rubberSphere(rubberBallRadius);
    for (int i=0; i<NRubberBalls; ++i) {
        contact.addSphere(balls[i], Vec3(0), rubberBallRadius, krubber, crubber);
        vtk.addDecoration(balls[i], Transform(), 
            rubberSphere.setColor((Real(i+1)/NRubberBalls)*Gray));
    }

    DecorativeSphere hardSphere(hardBallRadius);
    for (int i=0; i<NHardBalls; ++i) {
        contact.addSphere(balls[NRubberBalls+i], Vec3(0), hardBallRadius, khard, chard);
        vtk.addDecoration(balls[NRubberBalls+i], Transform(), 
            hardSphere.setColor(colors[i % NColors]));

    }


    State s;
    mbs.realize(s, Stage::Topology);
    bool suppressProjection = false;
    RungeKuttaMerson ee(mbs, s, suppressProjection);
    //CPodesIntegrator ee(mbs, s);
    //ee.setProjectEveryStep(true);

    vtk.report(s);

    // set Modeling stuff (s)
    //   none
    //bouncers.setUseEulerAngles(s,true);

    std::vector<State> saveEm;

    const Real h = .1;
    const Real tstart = 0.;
    const Real tmax = 100;

    ee.setAccuracy(1e-2);
    ee.setConstraintTolerance(1e-3);

    s.updTime() = tstart;
    mbs.realize(s, Stage::Acceleration);
    for (int i=0; i<100; ++i)
        saveEm.push_back(s);    // delay
    vtk.report(s);

    ee.initialize();
    for (int i=0; i<100; ++i)
        saveEm.push_back(s);    // delay
    vtk.report(s);

    int step = 0;
    while (s.getTime() <= tmax) {
        if (!(step % 10)) {
            mbs.realize(s);
            cout << s.getTime() << ": E=" << mbs.getEnergy(s)
             << " (pe=" << mbs.getPotentialEnergy(s)
             << ", ke=" << mbs.getKineticEnergy(s)
             << ") qerr=" << bouncers.calcQConstraintNorm(s)
             << " uerr=" << bouncers.calcUConstraintNorm(s)
             << " hNext=" << ee.getPredictedNextStep() << endl;
        }
        ++step;

        ee.step(s.getTime() + h);
        vtk.report(s);
        saveEm.push_back(s);
    }

    while(true) {
        for (int i=0; i < (int)saveEm.size(); ++i) {
            vtk.report(saveEm[i]);
            //vtk.report(saveEm[i]); // half speed
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
