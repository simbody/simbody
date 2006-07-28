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
static const int  Ground = 0; // ground is always body 0


void stateTest() {
  try {
    State s;
    s.setNSubsystems(1);
    s.advanceSubsystemToStage(0, Stage::Built);
    s.advanceSystemToStage(Stage::Built);

    Vector v3(3), v2(2);
    long q1 = s.allocateQ(0, v3);
    long q2 = s.allocateQ(0, v2);

    printf("q1,2=%d,%d\n", q1, q2);
    cout << s;

    long dv = s.allocateDiscreteVariable(0, Stage::Dynamics, new Value<int>(5));

    s.advanceSubsystemToStage(0, Stage::Modeled);
        //long dv2 = s.allocateDiscreteVariable(0, Stage::Configured, new Value<int>(5));

    Value<int>::downcast(s.updDiscreteVariable(0, dv)) = 71;
    cout << s.getDiscreteVariable(0, dv) << endl;

    s.advanceSystemToStage(Stage::Modeled);

    cout << s;

  }
  catch(const std::exception& e) {
      printf("*** STATE TEST EXCEPTION\n%s\n***\n", e.what());
  }

}

extern "C" void SimTK_version_SimTKlapack(int*,int*,int*);
extern "C" void SimTK_about_SimTKlapack(const char*, int, char*);

int main() {
    //stateTest();
    //exit(0);

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
    SimbodyMatterSubsystem test;


    //
    // G | *R |------ | A | ----- | B | ----- | C |
    // Bodies A,B,C are connected by sliding joints aligned with the x axis.
    // Body R is connected to ground origin by pin joint in z.
    // Ref config has everyone piled up at the origin.

    Real cubeSide = 1.; 
    Real m = 3.;

    // Sliders slide along the z axes of Mb and M, so we need a frame whose z axis
    // points in the x direction. (This one is 00-1,010,100.)
    RotationMat RzAlongX; 
    RzAlongX.setToBodyFixed123(Vec3(0,Pi/2,0));
    cout << "RzAlongX=" << RzAlongX;

    Transform jointFrame(RzAlongX,Vec3(0));

    MassProperties cube(m, Vec3(0), m*InertiaMat::brick(cubeSide/2,cubeSide/2,cubeSide/2));
    cout << "cube mprops about body frame: " << cube.getMass() << ", " 
        << cube.getCOM() << ", " << cube.getInertia() << endl;

    const int bodyR =
      test.addRigidBody(cube, Transform(),
                        Ground, Transform(), 
                        Mobilizer(Mobilizer::Pin, false)
                        );

    const int bodyA = 
      test.addRigidBody(cube, jointFrame,
                        bodyR, jointFrame, 
                        Mobilizer(Mobilizer::Sliding, false)
                        );
    const int bodyB = 
      test.addRigidBody(cube, jointFrame,
                        bodyA, jointFrame, 
                        Mobilizer(Mobilizer::Sliding, false)
                        );
    const int bodyC = 
      test.addRigidBody(cube, jointFrame,
                        bodyB, jointFrame, 
                        Mobilizer(Mobilizer::Sliding, false)
                        );


    MultibodySystem mbs;
    mbs.setMatterSubsystem(test);
    GeneralForceElements springs;
    mbs.addForceSubsystem(springs);

    State s;
    mbs.realize(s, Stage::Built);
    //cout << "mbs State as built: " << s;

    mbs.realize(s, Stage::Modeled);
    //cout << "mbs State as modeled: " << s;

    test.setJointQ(s,bodyR,0,0.);
    test.setJointQ(s,bodyA,0,1.);
    test.setJointQ(s,bodyB,0,1.);
    test.setJointQ(s,bodyC,0,1.);

    mbs.realize(s, Stage::Configured);

    const Real w = 1.; // rad/sec around z
    test.setJointU(s,bodyR,0,w); 
    mbs.realize(s, Stage::Moving);

    Transform bodyRConfig = test.getBodyConfiguration(s, bodyR);
    Transform bodyAConfig = test.getBodyConfiguration(s, bodyA);
    Transform bodyBConfig = test.getBodyConfiguration(s, bodyB);
    Transform bodyCConfig = test.getBodyConfiguration(s, bodyC);
    cout << "q=" << s.getQ() << endl;
    cout << "body R frame=" << bodyAConfig;
    cout << "body A frame=" << bodyAConfig;
    cout << "body B frame=" << bodyBConfig;
    cout << "body C frame=" << bodyCConfig;

    cout << "body velocities:\n";
    for (int i=0; i < test.getNBodies(); ++i)
        cout << "  V_GB[" << i << "]=" << test.getBodyVelocity(s, i) << endl;

    Vector_<SpatialVec> dEdR(test.getNBodies());
    dEdR.setToZero();
    dEdR[bodyC] = SpatialVec(Vec3(0), Vec3(-.5,0,0));
    Vector dEdQ;
    test.calcInternalGradientFromSpatial(s, dEdR, dEdQ);
    cout << "dEdR=" << dEdR << endl;
    cout << "dEdQ=" << dEdQ << endl;

    Vector_<SpatialVec> bodyForces;
    Vector_<Vec3>       particleForces;
    Vector              mobilityForces;

    test.resetForces(bodyForces, particleForces, mobilityForces);
    bodyForces[bodyC] = SpatialVec(Vec3(0), Vec3(-.5,0,0));

    mbs.realize(s, Stage::Dynamics);

    cout << "body centrifugal forces a,b,Ma+b:\n";
    for (int i=0; i < test.getNBodies(); ++i) {
        cout << "  body[" << i << "]=" << test.getCoriolisAcceleration(s,i) << ", "
            << test.getGyroscopicForce(s,i) << ", "
            << test.getCentrifugalForces(s, i) << endl;
    }


    Vector equivT;
    test.calcTreeEquivalentJointForces(s, bodyForces, equivT);
    cout << "body forces: " << bodyForces << endl;
    cout << "-> equiv joint forces: " << equivT << endl;


    Vector_<SpatialVec> bodyAccs;
    Vector udot;
    test.calcTreeUDot(s,mobilityForces,bodyForces,udot,bodyAccs);
    cout << "Udot from F=" << udot << endl;
    for (int i=0; i < test.getNBodies(); ++i)
        cout << "  A_GB[" << i << "]=" << bodyAccs[i] << endl;

    test.resetForces(bodyForces, particleForces, mobilityForces);
    test.calcTreeUDot(s,equivT,bodyForces,udot,bodyAccs);
    cout << "Udot from equiv=" << udot << endl;
    for (int i=0; i < test.getNBodies(); ++i)
        cout << "  A_GB[" << i << "]=" << bodyAccs[i] << endl;


    test.calcMInverseF(s,equivT,udot,bodyAccs);
    cout << "M^-1 f from equiv=" << udot << endl;
    for (int i=0; i < test.getNBodies(); ++i)
        cout << "  A_GB[" << i << "]=" << bodyAccs[i] << endl;

    const Real p[]={0,17.5,14.5,8.5};
    Vector sdfastFrcMat(4, p);
    cout << "sdfastFrcMat=" << sdfastFrcMat << endl;

    test.calcMInverseF(s,sdfastFrcMat,udot,bodyAccs);
    cout << "M^-1 f from sdfastFrcMat=" << udot << endl;
    for (int i=0; i < test.getNBodies(); ++i)
        cout << "  A_GB[" << i << "]=" << bodyAccs[i] << endl;

    test.calcTreeUDot(s,dEdQ,bodyForces,udot,bodyAccs);
    cout << "Udot from dEdQ=" << udot << endl;
    for (int i=0; i < test.getNBodies(); ++i)
        cout << "  A_GB[" << i << "]=" << bodyAccs[i] << endl;

    test.calcMInverseF(s,dEdQ,udot,bodyAccs);
    cout << "M^-1 f from dEdQ=" << udot << endl;
    for (int i=0; i < test.getNBodies(); ++i)
        cout << "  A_GB[" << i << "]=" << bodyAccs[i] << endl;

    mbs.realize(s, Stage::Reacting);

    SpatialVec bodyRAcc = test.getBodyAcceleration(s, bodyR);
    SpatialVec bodyAAcc = test.getBodyAcceleration(s, bodyA);
    SpatialVec bodyBAcc = test.getBodyAcceleration(s, bodyB);
    SpatialVec bodyCAcc = test.getBodyAcceleration(s, bodyC);
    cout << "body acc R,A,B,C: " << bodyRAcc << ", " 
         << bodyAAcc << ", " << bodyBAcc << ", " << bodyCAcc<< endl;

}
catch (const std::exception& e) {
    printf("EXCEPTION THROWN: %s\n", e.what());
}
    return 0;
}
