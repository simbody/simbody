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
 * A one-body pendulum, to test proper frame alignment and basic
 * functioning of Simbody.
 */

/* Sketch:
 *
 *     |           \           | g
 *     *--          *--        v
 *    / G          / Jb
 *
 *
 *   |           |
 *   *==---------*==---------W
 *  / J         / B         weight
 *   <--- L/2 ---|--- L/2 --->
 *
 *
 * The pendulum is a massless rod with origin frame
 * B, joint attachment frame J, and a point mass W.
 * The rod length is L, with the joint and mass
 * located in opposite directions along the B
 * frame X axis.
 *
 * There is a frame Jb on Ground which will connect
 * to J via a torsion joint around their mutual z axis.
 * Gravity is in the -y direction of the Ground frame.
 * Note that Jb may not be aligned with G, and J may
 * differ from B so the reference configuration may 
 * involve twisting the pendulum around somewhat.
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
static const int  GroundBodyNum = 0; // ground is always body 0


void stateTest() {
  try {
    State s;
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
    stateTest();
    //exit(0);

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
    SimbodyMatterSubsystem pend;

    Real L = 1.; 
    Real m = 3.;
    Real g = 9.8;
    Transform groundFrame;
    Transform baseFrame;

    //int baseBody =
     //   pend.addRigidBody(0, groundFrame, 
     //                     Mobilizer(Mobilizer::Pin, false),
      //                    Transform(), MassProperties(0.,Vec3(0.),InertiaMat(0.)));
    Transform jointFrame(Vec3(-L/2,0,0));
    MassProperties mprops(m, Vec3(L/2,0,0), InertiaMat(Vec3(L/2,0,0), m)+InertiaMat(1e-6,1e-6,1e-6));
    cout << "mprops about body frame: " << mprops.getMass() << ", " 
        << mprops.getCOM() << ", " << mprops.getInertia() << endl;

    Vec3 gravity(0.,-g,0.);
    cout << "period should be " << 2*std::acos(-1.)*std::sqrt(L/g) << " seconds." << endl;
    int theBody = 
      pend.addRigidBody(mprops, jointFrame,
                        0, Transform(), 
                        //Mobilizer(Mobilizer::Cartesian, false)
                        //Mobilizer(Mobilizer::Sliding, false)
                        //Mobilizer(Mobilizer::Pin, false)
                        //Mobilizer(Mobilizer::Ball, false)
                        Mobilizer(Mobilizer::Free, false)
                        );

/*
    int secondBody = 
      pend.addRigidBody(mprops, jointFrame,
                        theBody, Transform(Vec3(L/2,0,0)), 
                        //Mobilizer(Mobilizer::Cartesian, false),
                        //Mobilizer(Mobilizer::Sliding, false),
                        Mobilizer(Mobilizer::Pin, false),
                        //Mobilizer(Mobilizer::Ball, false),
                        //Mobilizer(Mobilizer::Free, false));
*/
    //int theConstraint =
    //    pend.addConstantDistanceConstraint(0, Vec3((L/2)*std::sqrt(2.)+1,1,0),
    //                                       theBody, Vec3(0,0,0),
    //                                       L/2+std::sqrt(2.));
 
    int ballConstraint =
        pend.addCoincidentStationsConstraint(0, Transform().T(),
                                             theBody, jointFrame.T()); 
/*
    Transform harderOne;
    harderOne.updR().setToBodyFixed123(Vec3(.1,.2,.3));
    harderOne.updT() = jointFrame.T()+Vec3(.1,.2,.3);
    int weldConstraint =
        pend.addWeldConstraint(0, Transform(),
                              theBody, harderOne);    
*/

    //pend.endConstruction();

    EmptyForcesSubsystem noForces;


    MultibodySystem mbs;
    mbs.addMatterSubsystem(pend);

    const Vec3 attachPt(1.5, 1, 0);
    TwoPointSpringSubsystem spring1(
        0, attachPt, 
        1, Vec3(L/2,0,0), 
        20, 1);
    mbs.addForceSubsystem(spring1);


    VTKReporter vtk(mbs);

    DecorativeLine rbProto; rbProto.setColor(Orange).setLineThickness(3);
    vtk.addRubberBandLine(0, attachPt, 1,Vec3(L/2,0,0), rbProto);

    DecorativeSphere sphere(0.25);
    sphere.setRepresentationToPoints();
    sphere.setResolution(2);
    vtk.addDecoration(0, Transform(Vec3(1,2,3)), sphere);
    sphere.setScale(0.5); sphere.setResolution(1);
    vtk.addDecoration(1, Transform(Vec3(0.1,0.2,0.3)), sphere);
    Quaternion qqq; qqq.setToAngleAxis(Pi/4, UnitVec3(1,0,0));
    vtk.addDecoration(1, Transform(RotationMat(qqq), Vec3(0,1,0)), DecorativeBrick(Vec3(.5,.1,.25)));
    DecorativeCylinder cyl(0.1); cyl.setOpacity(0.3);
    vtk.addDecoration(1, Transform(Vec3(-1,0,0)), 
        DecorativeCylinder(0.1).setOpacity(0.3));

    vtk.addDecoration(1, Transform(Vec3(3, 0, 0)), DecorativeSphere().setColor(Black));
    vtk.addDecoration(1, Transform(Vec3(3, 0.5, 0)), DecorativeSphere().setColor(Gray));
    vtk.addDecoration(1, Transform(Vec3(3, 1, 0)), DecorativeSphere().setColor(White));
    vtk.addDecoration(1, Transform(Vec3(3, 1.5, 0)), DecorativeSphere().setColor(Red));
    vtk.addDecoration(1, Transform(Vec3(3, 2, 0)), DecorativeSphere().setColor(Green));
    vtk.addDecoration(1, Transform(Vec3(3, 2.5, 0)), DecorativeSphere().setColor(Blue));
    vtk.addDecoration(1, Transform(Vec3(3, 3, 0)), DecorativeSphere().setColor(Yellow));
    vtk.addDecoration(1, Transform(Vec3(3, 3.5, 0)), DecorativeSphere().setColor(Orange));
    vtk.addDecoration(1, Transform(Vec3(3, 4, 0)), DecorativeSphere().setColor(Magenta));
    vtk.addDecoration(1, Transform(Vec3(3, 4.5, 0)), DecorativeSphere().setColor(Cyan));
    vtk.addDecoration(1, Transform(Vec3(3, 5, 0)), DecorativeSphere().setColor(Purple));
    State s;
    mbs.realize(s, Stage::Built);
    cout << "mbs State as built: " << s;

    //ExplicitEuler ee(mbs, s);
    bool suppressProjection = false;
    RungeKuttaMerson ee(mbs, s, suppressProjection);
    ee.setProjectEveryStep(true);

    vtk.report(s);

    // set Modeling stuff (s)
    pend.setUseEulerAngles(s, false); // this is the default
    //pend.setUseEulerAngles(s, true);
    mbs.realize(s, Stage::Modeled);

    cout << "mbs State as modeled: " << s;

    spring1.updGravity(s) = Vec3(0,-9.8,0);
    spring1.updStiffness(s) = 1000;


    //pend.setJointQ(s,1,0,0);
   // pend.setJointQ(s,1,3,-1.1);
   // pend.setJointQ(s,1,4,-2.2);
   // pend.setJointQ(s,1,5,-3.3);

    mbs.realize(s, Stage::Configured);

    Transform bodyConfig = pend.getBodyConfiguration(s, theBody);
    cout << "q=" << s.getQ() << endl;
    cout << "body frame: " << bodyConfig;

    pend.enforceConfigurationConstraints(s, 1e-3, 1e-10);
    mbs.realize(s, Stage::Configured);

    cout << "-------> STATE after realize(Configured):" << s;
    cout << "<------- STATE after realize(Configured)." << endl;

    cout << "after assembly body frame: " << pend.getBodyConfiguration(s,theBody); 

    Vector_<SpatialVec> dEdR(pend.getNBodies());
    dEdR[0] = 0;
    for (int i=1; i < pend.getNBodies(); ++i)
        dEdR[i] = SpatialVec(Vec3(0), Vec3(0.,2.,0.));
    Vector dEdQ;
    pend.calcInternalGradientFromSpatial(s, dEdR, dEdQ);
    cout << "dEdR=" << dEdR << endl;
    cout << "dEdQ=" << dEdQ << endl;

    pend.setJointU(s, 1, 0, 10.);

    Vector_<SpatialVec> bodyForces;
    Vector_<Vec3>       particleForces;
    Vector              mobilityForces;

    pend.resetForces(bodyForces, particleForces, mobilityForces);
    pend.addInGravity(s, gravity, bodyForces);
    pend.addInMobilityForce(s, 1, 0, 147, mobilityForces);

    mbs.realize(s, Stage::Moving);
    SpatialVec bodyVel = pend.getBodyVelocity(s, theBody);
    cout << "body vel: " << bodyVel << endl;

    cout << "wXwXr=" << bodyVel[0] % (bodyVel[0] % Vec3(2.5,0,0)) << endl;


    cout << "after applying gravity, body forces=" << bodyForces << endl;
    cout << "   joint forces=" << mobilityForces << endl;

    mbs.realize(s, Stage::Dynamics);
    Vector equivT;
    pend.calcTreeEquivalentJointForces(s, bodyForces, equivT);
    cout << "body forces -> equiv joint forces=" << equivT << endl;

    mbs.realize(s, Stage::Reacting);

    SpatialVec bodyAcc = pend.getBodyAcceleration(s, theBody);
    cout << "body acc: " << bodyAcc << endl;

    //pend.updQ(s) = Vector(4, &Vec4(1.,0.,0.,0.)[0]);
    //pend.updQ(s)[0] = -1.5; // almost hanging straight down
    pend.setJointU(s, 1, 0,   0.);
    //pend.setJointU(s, 1, 1,   0.);
    //pend.setJointU(s, 1, 2, -10.);
   // pend.setJointU(s, 1, 2,   0.);

    const Real angleInDegrees = 45;
    const Vec4 aa(angleInDegrees*RadiansPerDegree,0, 0, 1);
    Quaternion q; q.setToAngleAxis(aa);
    pend.setMobilizerConfiguration(s,1,Transform(RotationMat(q), Vec3(.1,.2,.3)));
    vtk.report(s);

    //pend.updQ(s)[2] = -.1;
    //pend.setJointQ(s, 1, 2, -0.999*std::acos(-1.)/2);

    const Real h = .01;
    const Real tstart = 0.;
    const Real tmax = 100;

    ee.setAccuracy(1e-6);
    ee.setConstraintTolerance(1e-9);

    ee.initialize(); 
    vtk.report(s);
    s.updTime() = tstart;
    int step = 0;
    while (s.getTime() < tmax) {
        ee.step(s.getTime() + h);
        cout << " E=" << mbs.getPotentialEnergy(s)+mbs.getKineticEnergy(s)
             << " (pe=" << mbs.getPotentialEnergy(s)
             << ", ke=" << mbs.getKineticEnergy(s)
             << ") hNext=" << ee.getPredictedNextStep() << endl;

        const Vector qdot = pend.getQDot(s);

        Transform  x = pend.getBodyConfiguration(s,theBody);
        SpatialVec v = pend.getBodyVelocity(s,theBody);

        //Vec3 err = x.T()-Vec3(2.5,0.,0.);
        //Real d = err.norm();
        //Real k = m*gravity.norm(); // stiffness, should balance at 1
        // Real c = 10.; // damping
        //Vec3 fk = -k*err;
        //Real fc = -c*pend.getU(s)[2];
        //pend.applyPointForce(s,theBody,Vec3(0,0,0),fk);
        //pend.applyJointForce(s,theBody,2,fc);

        if (!(step % 10)) {
            cout << s.getTime() << " " 
                 << s.getQ() << " " << s.getU() 
                 << " hNext=" << ee.getPredictedNextStep() << endl;
            cout << "body config=" << x;
            cout << "body velocity=" << v << endl;
            //cout << "err=" << err << " |err|=" << d << endl;
            //cout << "spring force=" << fk << endl;
            //cout << "damping joint forces=" << fc << endl;
            
            vtk.report(s);
        }
        ++step;


        mbs.realize(s, Stage::Reacting);
        /*
        cout << "CONSTRAINT ERRORS:\n";
        cout << "quat:" << Vec4::getAs(&s.getQ()[0]).norm()-1 << endl;
        cout << "   q:" << pend.getQConstraintErrors(s)
             << "(" << pend.calcQConstraintNorm(s) << ")\n";
        cout << "   u:" << pend.getUConstraintErrors(s)
             << "(" << pend.calcUConstraintNorm(s) << ")\n";
        cout << "udot: " << pend.getUDotConstraintErrors(s)
             << "(" << pend.calcUDotConstraintNorm(s) << ")\n\n";
        */

        const Vector udot = s.getUDot();
        Vector udot2;
        pend.calcTreeUDot(s, 
            pend.getAppliedMobilityForces(s),
            pend.getAppliedBodyForces(s),
            udot2);
        if (!(step % 100)) {
            cout << "udot = " << udot << endl;
            cout << "udot2= " << udot2 << endl;
        }
    }

}
catch (const std::exception& e) {
    printf("EXCEPTION THROWN: %s\n", e.what());
}
    return 0;
}
