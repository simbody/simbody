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
 * There is a frame Jb on GroundIndex which will connect
 * to J via a torsion joint around their mutual z axis.
 * Gravity is in the -y direction of the GroundIndex frame.
 * Note that Jb may not be aligned with G, and J may
 * differ from B so the reference configuration may 
 * involve twisting the pendulum around somewhat.
 */

#include "SimTKsimbody.h"

#include <string>
#include <iostream>
#include <exception>
#include <cmath>
using std::cout;
using std::endl;

using namespace SimTK;

static const Real RadiansPerDegree = Pi/180;

void stateTest() {
  try {
    State s;
    s.setNumSubsystems(1);
    s.advanceSubsystemToStage(SubsystemIndex(0), Stage::Topology);
    s.advanceSystemToStage(Stage::Topology);

    Vector v3(3), v2(2);
    int q1 = s.allocateQ(SubsystemIndex(0), v3);
    int q2 = s.allocateQ(SubsystemIndex(0), v2);

    printf("q1,2=%d,%d\n", q1, q2);
    cout << s;

    DiscreteVariableIndex dv = s.allocateDiscreteVariable(SubsystemIndex(0), Stage::Dynamics, new Value<int>(5));

    s.advanceSubsystemToStage(SubsystemIndex(0), Stage::Model);
        //long dv2 = s.allocateDiscreteVariable(SubsystemIndex(0), Stage::Position, new Value<int>(5));

    Value<int>::downcast(s.updDiscreteVariable(SubsystemIndex(0), dv)) = 71;
    cout << s.getDiscreteVariable(SubsystemIndex(0), dv) << endl;

    s.advanceSystemToStage(Stage::Model);

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

    SimTK_version_simbody(&major,&minor,&build);
    std::printf("==> simbody library version: %d.%d.%d\n", major, minor, build);
    std::printf("    SimTK_about_simbody():\n");
    for (const char** p = keylist; *p; ++p) {
        SimTK_about_simbody(*p, 100, out);
        std::printf("      about(%s)='%s'\n", *p, out);
    }


try {
    MultibodySystem mbs;
    SimbodyMatterSubsystem pend(mbs);
    GeneralForceSubsystem springs(mbs);
    HuntCrossleyContact contact(mbs);

    Real L = 1.; 
    Real m = 3.;
    Real g = 9.8;
    Transform groundFrame;
    Transform baseFrame;

    Transform jointFrame(Vec3(-L/2,0,0));
    MassProperties mprops(m, Vec3(L/2,0,0), Inertia(Vec3(L/2,0,0), m)+Inertia(1e-6,1e-6,1e-6));
    cout << "mprops about body frame: " << mprops.getMass() << ", " 
        << mprops.getMassCenter() << ", " << mprops.getInertia() << endl;

    Vec3 gravity(0.,-g,0.);
    Force::UniformGravity gravityForces(springs, pend, gravity);
    cout << "period should be " << 2*std::acos(-1.)*std::sqrt(L/g) << " seconds." << endl;

    MobilizedBody::Free aPendulum(pend.Ground(), Transform(), // ground, at origin
                                  Body::Rigid(mprops), jointFrame);

    const Real ballMass = 10;
    const Real ballRadius = 2;
    const MassProperties ballMProps(ballMass, Vec3(0), ballMass*UnitInertia::sphere(ballRadius));
    const Vec3 ballPos = Vec3(-3,5,0);

    MobilizedBody::Cartesian aBall(pend.Ground(), Transform(ballPos),
                                   Body::Rigid(ballMProps), Transform());

    MobilizedBody::Cartesian aBall2(pend.Ground(), Transform(ballPos+Vec3(0.1,10,0)),
                                    Body::Rigid(ballMProps), Transform());

    Constraint::Ball ballConstraint(pend.Ground(), Transform().p(),
                                    aPendulum, jointFrame.p());

    const Vec3 attachPt(1.5, 1, 0);
    Force::TwoPointLinearSpring(springs, pend.Ground(), attachPt, 
        aPendulum, Vec3(L/2,0,0), 
        100, 1);

    const Real k = 1000, c = 0.0;
    contact.addHalfSpace(GroundIndex, UnitVec3(0,1,0), 0, k, c); // h,k,c
    contact.addHalfSpace(GroundIndex, UnitVec3(1,0,0), -10, k, c); // h,k,c
    contact.addHalfSpace(GroundIndex, UnitVec3(-1,0,0), -10, k, c); // h,k,c

    contact.addSphere(aBall, Vec3(0), ballRadius, k, c); // r,k,c
    contact.addSphere(aBall2, Vec3(0), ballRadius, k, c); // r,k,c

    State s = mbs.realizeTopology();
    cout << "mbs State as built: " << s;

    Visualizer viz(mbs);
    viz.setBackgroundType(Visualizer::SolidColor);

    viz.addDecoration(GroundIndex, Transform(), DecorativeBrick(Vec3(20,.1,20)).setColor(1.5*Gray).setOpacity(.3));
    viz.addDecoration(GroundIndex, Transform(Vec3(-10,0,0)), DecorativeBrick(Vec3(.1,20,20)).setColor(Yellow).setOpacity(1));
    viz.addDecoration(GroundIndex, Transform(Vec3(10,0,0)), DecorativeBrick(Vec3(.1,20,20)).setColor(Yellow).setOpacity(1));

    DecorativeSphere bouncer(ballRadius);
    viz.addDecoration(aBall, Transform(), bouncer.setColor(Orange));
    viz.addDecoration(aBall2, Transform(), bouncer.setColor(Blue));

    DecorativeLine rbProto; rbProto.setColor(Orange).setLineThickness(3);
    viz.addRubberBandLine(GroundIndex, attachPt, aPendulum, Vec3(L/2,0,0), rbProto);

    DecorativeSphere sphere(0.25);
    sphere.setRepresentation(DecorativeGeometry::DrawPoints);
    sphere.setResolution(2);
    viz.addDecoration(GroundIndex, Transform(Vec3(1,2,3)), sphere);
    sphere.setScale(0.5); sphere.setResolution(1);
    viz.addDecoration(aPendulum, Transform(Vec3(0.1,0.2,0.3)), sphere);
    Quaternion qqq; qqq.setQuaternionFromAngleAxis(Pi/4, UnitVec3(1,0,0));
    viz.addDecoration(aPendulum, Transform(Rotation(qqq), Vec3(0,1,0)), DecorativeBrick(Vec3(.5,.1,.25)));
    DecorativeCylinder cyl(0.1); cyl.setOpacity(0.3);
    viz.addDecoration(aPendulum, Transform(Vec3(-1,0,0)), 
        DecorativeCylinder(0.1).setOpacity(0.3));

    viz.addDecoration(aPendulum, Transform(Vec3(3, 0, 0)), DecorativeSphere().setColor(Black));
    viz.addDecoration(aPendulum, Transform(Vec3(3, 0.5, 0)), DecorativeSphere().setColor(Gray));
    viz.addDecoration(aPendulum, Transform(Vec3(3, 1, 0)), DecorativeSphere().setColor(White));
    viz.addDecoration(aPendulum, Transform(Vec3(3, 1.5, 0)), DecorativeSphere().setColor(Red));
    viz.addDecoration(aPendulum, Transform(Vec3(3, 2, 0)), DecorativeSphere().setColor(Green));
    viz.addDecoration(aPendulum, Transform(Vec3(3, 2.5, 0)), DecorativeSphere().setColor(Blue));
    viz.addDecoration(aPendulum, Transform(Vec3(3, 3, 0)), DecorativeSphere().setColor(Yellow));
    viz.addDecoration(aPendulum, Transform(Vec3(3, 3.5, 0)), DecorativeSphere().setColor(Orange));
    viz.addDecoration(aPendulum, Transform(Vec3(3, 4, 0)), DecorativeSphere().setColor(Magenta));
    viz.addDecoration(aPendulum, Transform(Vec3(3, 4.5, 0)), DecorativeSphere().setColor(Cyan));
    viz.addDecoration(aPendulum, Transform(Vec3(3, 5, 0)), DecorativeSphere().setColor(Purple));


    viz.report(s);

    // set Modeling stuff (s)
    pend.setUseEulerAngles(s, false); // this is the default
    //pend.setUseEulerAngles(s, true);
    mbs.realizeModel(s);
    cout << "mbs State as modeled: " << s;

    printf("GLOBAL ny=%d q:y(%d,%d) u:y(%d,%d) z:y(%d,%d)\n",
        (int)s.getNY(), (int)s.getQStart(), (int)s.getNQ(), 
        (int)s.getUStart(), (int)s.getNU(), (int)s.getZStart(), (int)s.getNZ());

    mbs.realize(s, Stage::Instance);
    printf("  nyerr=%d qerr:yerr(%d,%d) uerr:yerr(%d,%d)\n",
        (int)s.getNYErr(), (int)s.getQErrStart(), (int)s.getNQErr(),
        (int)s.getUErrStart(), (int)s.getNUErr());
    printf("  nudoterr=%d\n", s.getNUDotErr());
    for (SubsystemIndex i(0); i<s.getNumSubsystems(); ++i) {
        printf("Subsys %d: q:y(%d,%d) u:y(%d,%d) z:y(%d,%d)\n",
            (int)i,(int)s.getQStart()+(int)s.getQStart(i),(int)s.getNQ(i),
              (int)s.getUStart()+(int)s.getUStart(i),(int)s.getNU(i),
              (int)s.getZStart()+(int)s.getZStart(i),(int)s.getNZ(i));
        printf("  qerr:yerr(%d,%d) uerr:yerr(%d,%d) uderr(%d,%d)\n",
            (int)s.getQErrStart()+(int)s.getQErrStart(i),(int)s.getNQErr(i),
            (int)s.getUErrStart()+(int)s.getUErrStart(i),(int)s.getNUErr(i),
            (int)s.getUDotErrStart(i),(int)s.getNUDotErr(i));
    }


    //pend.setJointQ(s,1,0,0);
   // pend.setJointQ(s,1,3,-1.1);
   // pend.setJointQ(s,1,4,-2.2);
   // pend.setJointQ(s,1,5,-3.3);

    mbs.realize(s, Stage::Position);
    Vector wts; mbs.calcYUnitWeights(s,wts);
    cout << "Y WEIGHTS: " << wts << endl;
    Vector tols; mbs.calcYErrUnitTolerances(s,tols);
    cout << "YERR TOLERANCES: " << tols << endl;

    Transform bodyConfig = aPendulum.getBodyTransform(s);
    cout << "q=" << s.getQ() << endl;
    cout << "body frame: " << bodyConfig;

    Vector dummy;
    pend.projectQConstraints(s, 1e-3, wts, tols, dummy, System::ProjectOptions::All);
    mbs.realize(s, Stage::Position);

    cout << "-------> STATE after realize(Position):" << s;
    cout << "<------- STATE after realize(Position)." << endl;

    cout << "after assembly body frame: " << aPendulum.getBodyTransform(s); 

    Vector_<SpatialVec> dEdR(pend.getNumBodies());
    dEdR[0] = 0;
    for (int i=1; i < pend.getNumBodies(); ++i)
        dEdR[i] = SpatialVec(Vec3(0), Vec3(0.,2.,0.));
    Vector dEdQ;
    pend.calcInternalGradientFromSpatial(s, dEdR, dEdQ);
    cout << "dEdR=" << dEdR << endl;
    cout << "dEdQ=" << dEdQ << endl;

    pend.getMobilizedBody(MobilizedBodyIndex(1)).setOneU(s,0,10.);

    Vector_<SpatialVec> bodyForces(pend.getNumBodies());
    Vector_<Vec3>       particleForces(pend.getNumParticles());
    Vector              mobilityForces(pend.getNumMobilities());

    bodyForces.setToZero();
    particleForces.setToZero();
    mobilityForces.setToZero();

    pend.addInMobilityForce(s, aPendulum, MobilizerUIndex(0), 147, mobilityForces);

    mbs.realize(s, Stage::Velocity);
    SpatialVec bodyVel = aPendulum.getBodyVelocity(s);
    cout << "body vel: " << bodyVel << endl;

    cout << "wXwXr=" << bodyVel[0] % (bodyVel[0] % Vec3(2.5,0,0)) << endl;


    cout << "after applying gravity, body forces=" << bodyForces << endl;
    cout << "   joint forces=" << mobilityForces << endl;

    mbs.realize(s, Stage::Dynamics);
    Vector equivT;
    pend.calcTreeEquivalentMobilityForces(s, bodyForces, equivT);
    cout << "body forces -> equiv joint forces=" << equivT << endl;

    mbs.realize(s, Stage::Acceleration);

    SpatialVec bodyAcc = aPendulum.getBodyAcceleration(s);
    cout << "body acc: " << bodyAcc << endl;

    aPendulum.setOneU(s, 0, 0.);

    const Real angleInDegrees = 45;
    const Vec4 aa(angleInDegrees*RadiansPerDegree,0, 0, 1);
    Quaternion q; q.setQuaternionFromAngleAxis(aa);
    aPendulum.setQToFitTransform(s,Transform(Rotation(q), Vec3(.1,.2,.3)));
    viz.report(s);

    //pend.updQ(s)[2] = -.1;
    //pend.setJointQ(s, 1, 2, -0.999*std::acos(-1.)/2);

    const Real h = 1./30;
    const Real tstart = 0.;
    const Real tmax = 20;

    RungeKuttaMersonIntegrator ee(mbs);
    ee.setFinalTime(tmax);

    s.updTime() = tstart;
    ee.initialize(s);   // assemble if needed
    s = ee.getState();
    viz.report(s);

    Integrator::SuccessfulStepStatus status;
    int step = 0;
    while ((status=ee.stepTo(step*h)) != Integrator::EndOfSimulation) {
        const State& s = ee.getState();

        if (!(step % 10)) {
            // This is so we can calculate potential energy (although logically
            // one should be able to do that at Stage::Position).
            mbs.realize(s, Stage::Dynamics);

            cout << " E=" << mbs.calcEnergy(s)
                 << " (pe=" << mbs.calcPotentialEnergy(s)
                 << ", ke=" << mbs.calcKineticEnergy(s)
                 << ") hNext=" << ee.getPredictedNextStepSize() << endl;

            const Vector qdot = pend.getQDot(s);

            Transform  x = aPendulum.getBodyTransform(s);
            SpatialVec v = aPendulum.getBodyVelocity(s);

            //Vec3 err = x.p()-Vec3(2.5,0.,0.);
            //Real d = err.norm();
            //Real k = m*gravity.norm(); // stiffness, should balance at 1
            // Real c = 10.; // damping
            //Vec3 fk = -k*err;
            //Real fc = -c*pend.getU(s)[2];
            //pend.applyPointForce(s,aPendulum,Vec3(0,0,0),fk);
            //pend.applyJointForce(s,aPendulum,2,fc);

            cout << s.getTime() << " " 
                 << s.getQ() << " " << s.getU() 
                 << " hNext=" << ee.getPredictedNextStepSize() << endl;
            cout << "body config=" << x;
            cout << "body velocity=" << v << endl;
            //cout << "err=" << err << " |err|=" << d << endl;
            //cout << "spring force=" << fk << endl;
            //cout << "damping joint forces=" << fc << endl;
            
        }

        viz.report(s);

        if (status == Integrator::ReachedReportTime)
            ++step;


        if (!(step % 100)) {
            mbs.realize(s, Stage::Acceleration);
            const Vector udot = s.getUDot();
            cout << "udot = " << udot << endl;
        }
    }

}
catch (const std::exception& e) {
    printf("EXCEPTION THROWN: %s\n", e.what());
}
    return 0;
}
