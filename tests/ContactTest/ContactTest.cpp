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
static const int  Ground = 0; // ground is always body 0

extern "C" void SimTK_version_SimTKlapack(int*,int*,int*);
extern "C" void SimTK_about_SimTKlapack(const char*, int, char*);

static const int NColors =8;
static Vec3 colors[NColors] = {Red,Green,Blue,Yellow,Orange,Magenta,Cyan,Purple};

int main() {

    int major,minor,build;
    char out[100];
    const char* keylist[] = { "version", "library", "type", "debug", "authors", "copyright", "svn_revision", 0 };

    //SimTK_version_SimTKlapack(&major,&minor,&build);
    //std::printf("==> SimTKlapack library version: %d.%d.%d\n", major, minor, build);
    //std::printf("    SimTK_about_SimTKlapack():\n");
    //for (const char** p = keylist; *p; ++p) {
    //    SimTK_about_SimTKlapack(*p, 100, out);
    //    std::printf("      about(%s)='%s'\n", *p, out);
    //}

    SimTK_version_SimTKcommon(&major,&minor,&build);
    std::printf("==> SimTKcommon library version: %d.%d.%d\n", major, minor, build);
    std::printf("    SimTK_about_SimTKcommon():\n");
    for (const char** p = keylist; *p; ++p) {
        SimTK_about_SimTKcommon(*p, 100, out);
        std::printf("      about(%s)='%s'\n", *p, out);
    }
    SimTK_version_simmatrix(&major,&minor,&build);
    std::printf("==> simmatrix library version: %d.%d.%d\n", major, minor, build);
    std::printf("    SimTK_about_simmatrix():\n");
    for (const char** p = keylist; *p; ++p) {
        SimTK_about_simmatrix(*p, 100, out);
        std::printf("      about(%s)='%s'\n", *p, out);
    }
    SimTK_version_simbody(&major,&minor,&build);
    std::printf("==> simbody library version: %d.%d.%d\n", major, minor, build);
    std::printf("    SimTK_about_simbody():\n");
    for (const char** p = keylist; *p; ++p) {
        SimTK_about_simbody(*p, 100, out);
        std::printf("      about(%s)='%s'\n", *p, out);
    }


try {
    Real g = 9.8;   // m/s^2
    Vec3 gravity(0.,-g,0.);

    SimbodyMatterSubsystem bouncers;

    // start with a 2body pendulum
    Real linkLength = 10.; // m
    Real pendMass   = 10.; // kg
    Real pendBallRadius = 2;
    InertiaMat pendBallInertia(pendMass*InertiaMat::sphere(pendBallRadius));
    const MassProperties pendMProps(pendMass, Vec3(0, -linkLength/2, 0), 
        pendBallInertia.shiftFromCOM(Vec3(0, -linkLength/2, 0), pendMass));

    int pend1 = bouncers.addRigidBody(pendMProps, Transform(Vec3(0, linkLength/2, 0)),
                          Ground, Transform(Vec3(0,40,-5)),
                          Mobilizer(Mobilizer::Ball, false));

    int pend2 = bouncers.addRigidBody(pendMProps, Transform(Vec3(0, linkLength/2, 0)),
                          pend1, Transform(Vec3(0,-linkLength/2,0)),
                          Mobilizer(Mobilizer::Ball, false));


    const Real hardBallMass = 10, rubberBallMass = 100;  // kg
    const Real hardBallRadius = 1, rubberBallRadius = 1.5;  // m
    const MassProperties hardBallMProps(hardBallMass, Vec3(0), 
        hardBallMass*InertiaMat::sphere(hardBallRadius));
    const MassProperties rubberBallMProps(rubberBallMass, Vec3(0), 
        rubberBallMass*InertiaMat::sphere(rubberBallRadius));
    const Vec3 firstHardBallPos = Vec3(-3,10,0), firstRubberBallPos = Vec3(3,10,0);
    std::vector<int> balls;

    const int NRubberBalls = 3;
    for (int i=0; i<NRubberBalls; ++i)
        balls.push_back( 
            bouncers.addRigidBody(rubberBallMProps, Transform(),
                              Ground, Transform(firstRubberBallPos+i*Vec3(0,2*rubberBallRadius+1,0)),
                              Mobilizer(Mobilizer::Cartesian, false)));

    const int NHardBalls = 12;
    for (int i=0; i < NHardBalls; ++i)
        balls.push_back(
            bouncers.addRigidBody(hardBallMProps, Transform(),
                              Ground, Transform(firstHardBallPos+i*Vec3(0,2*hardBallRadius+1,0)
                                                + (i==NHardBalls-1)*Vec3(1e-15,0,1e-15)),
                              Mobilizer(Mobilizer::Cartesian, false)));


    MultibodySystem mbs;
    mbs.setMatterSubsystem(bouncers);

    UniformGravitySubsystem gravityForces(gravity);
    mbs.addForceSubsystem(gravityForces);

    HuntCrossleyContact contact;
    mbs.addForceSubsystem(contact);

    
    const Real kwall = 100000, khard=200000, krubber=20000;
    const Real cwall = 0.001, chard = 1e-6, crubber=0.1;
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

    contact.addHalfSpace(Ground, UnitVec3(0,1,0), 0, kwall, cwall);
    contact.addHalfSpace(Ground, UnitVec3(1,0,0), -20, kwall, cwall); // left
    contact.addHalfSpace(Ground, UnitVec3(-1,0,0), -20, kwall, cwall); // right
    contact.addHalfSpace(Ground, UnitVec3(0,0,-1), -20, kwall, cwall); // front
    contact.addHalfSpace(Ground, UnitVec3(0,0,1), -20, kwall, cwall); // back

    vtk.addDecoration(Ground, Transform(), 
        DecorativeBrick(Vec3(20,.1,20)).setColor(Green).setOpacity(1));
    vtk.addDecoration(Ground, Transform(Vec3(-20,20,0)), 
        DecorativeBrick(Vec3(.1,30,20)).setColor(Yellow).setOpacity(.3));
    vtk.addDecoration(Ground, Transform(Vec3(20,20,0)), 
        DecorativeBrick(Vec3(.1,30,20)).setColor(Yellow).setOpacity(.3));
    vtk.addDecoration(Ground, Transform(Vec3(0,20,20)), 
        DecorativeBrick(Vec3(20,20,.1)).setColor(Gray).setOpacity(.05));
    vtk.addDecoration(Ground, Transform(Vec3(0,20,-20)), 
        DecorativeBrick(Vec3(20,30,.1)).setColor(Cyan).setOpacity(.3));

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
    mbs.realize(s, Stage::Built);
    bool suppressProjection = false;
    RungeKuttaMerson ee(mbs, s, suppressProjection);
    ee.setProjectEveryStep(false);

    vtk.report(s);

    // set Modeling stuff (s)
    //   none

    vtk.report(s);

    std::vector<State> saveEm;

    const Real h = .1;
    const Real tstart = 0.;
    const Real tmax = 100;

    ee.setAccuracy(1e-2);
    ee.setConstraintTolerance(1e-3);

    s.updTime() = tstart;
    ee.initialize();
    for (int i=0; i<100; ++i)
        saveEm.push_back(s);    // delay
    vtk.report(s);

    int step = 0;
    while (s.getTime() <= tmax) {
        if (!(step % 10)) {
            cout << s.getTime() << ": E=" << mbs.getEnergy(s)
             << " (pe=" << mbs.getPotentialEnergy(s)
             << ", ke=" << mbs.getKineticEnergy(s)
             << ") hNext=" << ee.getPredictedNextStep() << endl;
        }
        ++step;

        ee.step(s.getTime() + h);
        vtk.report(s);
        saveEm.push_back(s);
    }

    while(true)
    for (int i=0; i < (int)saveEm.size(); ++i)
        vtk.report(saveEm[i]);
}
catch (const std::exception& e) {
    printf("EXCEPTION THROWN: %s\n", e.what());
}
    return 0;
}
