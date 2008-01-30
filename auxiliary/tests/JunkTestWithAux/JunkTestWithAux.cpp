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
 * This is a simple example of using a constraint.
 * Here we have two independent pendulums hanging from ground pins.
 * They can be connected by a spring or a distance constraint.
 */

#include "SimTKsimbody.h"
#include "SimTKsimbody_aux.h" // requires VTK

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

static const Vec3 hl(1, 0.5, 0.5); // body half lengths


// Create a free body and constrain it to a point on ground.
//
//                       * (1,2,0)
//                |
//           -----|-----
//          |     *     |   (0,0,0)
//           -----------

int main(int argc, char** argv) {
  try { // If anything goes wrong, an exception will be thrown.

        // CREATE MULTIBODY SYSTEM AND ITS SUBSYSTEMS
    MultibodySystem         mbs;

    SimbodyMatterSubsystem  matter(mbs);
    UniformGravitySubsystem gravity(mbs, Vec3(0, -g, 0));
    GeneralForceElements    forces(mbs);
    DecorationSubsystem     viz(mbs);

        // ADD BODIES AND THEIR MOBILIZERS
    const Body::Rigid body = Body::Rigid(MassProperties(m, Vec3(0), m*Inertia::brick(hl[0],hl[1],hl[2])))
                                  .addDecoration(Transform(), DecorativeBrick(hl).setOpacity(.5))
                                  .addDecoration(Transform(), DecorativeLine(Vec3(0), Vec3(0,1,0)).setColor(Green));

    MobilizedBody::Free mobilizedBody(matter.Ground(), Transform(), body, Transform());
    MobilizedBody::Free mobilizedBody0(mobilizedBody, Transform(Vec3(1,2,0)), body, Transform(Vec3(0,1,0)));
    //MobilizedBody::Ball mobilizedBody(mobilizedBody0, Transform(Vec3(1,2,0)), body, Transform(Vec3(0,1,0)));
    MobilizedBody::Free mobilizedBody2(mobilizedBody0, Vec3(-5,0,0), body, Transform());

#define HASC
    
    Constraint::Ball myc2(matter.Ground(), Vec3(-4,2,0),  mobilizedBody2, Vec3(0,1,0));
    Constraint::Ball myc(matter.Ground(), Vec3(1,2,0),  mobilizedBody, Vec3(0,1,0));
    Constraint::Ball ball(mobilizedBody0, Vec3(2,0,0), mobilizedBody2, Vec3(3,0,0));
    

    /*
    Constraint::Rod myc2(matter.Ground(), Vec3(-4,2,0),  mobilizedBody2, Vec3(0,1,0), 1);
    Constraint::Rod myc(matter.Ground(), Vec3(1,2,0),  mobilizedBody, Vec3(0,1,0), 1);
    Constraint::Rod ball(mobilizedBody0, Vec3(2,0,0), mobilizedBody2, Vec3(3,0,0), 5);
    */

    //Constraint::Rod rod(mobilizedBody, Vec3(0,0,1), mobilizedBody2, Vec3(-1,0,0), 5);
    //Constraint::Rod myc(matter.Ground(), Vec3(1,2,0),  mobilizedBody, Vec3(0,1,0), 1);

    //Constraint::Ball myc(mobilizedBody, Vec3(1,2,0),  mobilizedBody2, Vec3(0,1,0));
    //Constraint::PointInPlane myc(mobilizedBody, UnitVec3(0,1,0),  2, mobilizedBody2, Vec3(0,1,0));
    //Constraint::Rod myc(mobilizedBody, Vec3(0), mobilizedBody2, Vec3(0), 5);


    //Constraint::Rod myc(matter.Ground(), Vec3(1,2,0),  mobilizedBody, Vec3(0,1,0), 0.1);
    
    //Constraint::PointInPlane myc(matter.Ground(), UnitVec3(0,1,0), 2, // parallel to xz at height 2
    //                     mobilizedBody, Vec3(0,1,0));
    //Constraint::ConstantAngle myc(matter.Ground(), UnitVec3(1,1,1), mobilizedBody, UnitVec3(0,1,0), Pi/2);

    //Constraint::PointInPlane myc(matter.Ground(), UnitVec3(1,0,0), 1, 
    //                         mobilizedBody, Vec3(0,1,0));
     //Constraint::PointInPlane mycy(matter.Ground(), UnitVec3(0,1,0), 2, 
     //                        mobilizedBody, Vec3(0,1,0));
     //Constraint::PointInPlane mycz(matter.Ground(), UnitVec3(0,0,1), 0, 
     //                        mobilizedBody, Vec3(0,1,0));
    
    //Constraint::Ball redundant(matter.Ground(), Vec3(1,2,0), 
      //                   mobilizedBody, Vec3(0,1,0));
     //Constraint::Ball other(matter.Ground(), Vec3(0,0,0), 
      //                   mobilizedBody, Vec3(0,0,0));
    


    viz.addRubberBandLine(matter.Ground(), Vec3(1,2,0),
                          mobilizedBody, Vec3(0,1,0),
                          DecorativeLine().setColor(Orange).setLineThickness(4));
    viz.addBodyFixedDecoration(matter.Ground(), Vec3(1,2,0), DecorativeSphere(0.1).setOpacity(.2));
    viz.addBodyFixedDecoration(matter.Ground(), Transform(), DecorativeLine(Vec3(0),Vec3(1,1,1)).setColor(Black));

    State s = mbs.realizeTopology(); // returns a reference to the the default state
    //matter.setUseEulerAngles(s, true);
    mbs.realizeModel(s); // define appropriate states for this System

    //mobilizedBody0.setQ(s, .1);
    //mobilizedBody.setQ(s, .2);

    VTKVisualizer display(mbs);

    //{Real qq[] = {0.9719,0,0,-0.235396,0.542437,1.11082,0};
    // s.updQ()=Vector(4,qq);
    //}
    mbs.realize(s, Stage::Position);
    display.report(s);
    cout << "q=" << s.getQ() << endl;
    cout << "qErr=" << s.getQErr() << endl;

#ifdef HASC
    for (ConstraintIndex cid(0); cid < matter.getNConstraints(); ++cid) {
        const Constraint& c = matter.getConstraint(cid);
        int mp,mv,ma;
        c.getNumConstraintEquations(s, mp,mv,ma);

	    cout << "CONSTRAINT " << cid << " ancestor=" << c.getAncestorMobilizedBody().getMobilizedBodyIndex()
             << " " << c.getNumConstrainedBodies() << "constrained bodies, perr=" << c.getPositionError(s)
		     << endl;
        for (ConstrainedBodyIndex cid(0); cid < c.getNumConstrainedBodies(); ++cid)
            cout << "  constrained body: " << c.getConstrainedMobilizedBody(cid).getMobilizedBodyIndex() << endl;
        cout << c.getSubtree();
	    cout << "   d(perrdot)/du=" << c.calcPositionConstraintMatrixP(s);
        cout << "  ~d(Gt lambda)/dlambda=" << ~c.calcPositionConstraintMatrixPt(s);


	    cout << "   d(perr)/dq=" << c.calcPositionConstraintMatrixPQInverse(s);

        Matrix P = c.calcPositionConstraintMatrixP(s);
        Matrix PQ(mp,matter.getNQ(s));
        Vector out(matter.getNQ(s));
        Matrix Q(matter.getNQ(s), matter.getNU(s));
        Matrix Qinv(matter.getNU(s),matter.getNQ(s));
        for (int i=0; i<mp; ++i) {
            Vector in = ~P[i];
            matter.multiplyByQMatrixInverse(s, true, in, out);
            PQ[i] = ~out;
        }
        cout << " calculated d(perr)/dq=" << PQ;

        Vector u(matter.getNU(s)); u=0;
        for (int i=0; i < matter.getNU(s); ++i) {
            u[i]=1;
            matter.multiplyByQMatrix(s,false,u,Q(i));
            u[i]=0;
        }
        cout << " Q=" << Q;

    }
    const Constraint& c = matter.getConstraint(myc.getConstraintIndex());
    
#endif

    char ans;
    cout << "Default configuration shown. Ready? "; cin >> ans;

    mobilizedBody.setQToFitTransform (s, Transform(Rotation(.05,Vec3(1,1,1)),Vec3(.1,.2,.3)));
    mobilizedBody0.setQToFitTransform (s, Transform(Rotation(.05,Vec3(1,-1,1)),Vec3(.2,.2,.3)));
    mobilizedBody2.setQToFitTransform (s, Transform(Rotation(.05,Vec3(-1,1,1)),Vec3(.1,.2,.1)));
    mobilizedBody.setUToFitAngularVelocity(s, 100*Vec3(.1,.2,.3));

    mbs.realize(s, Stage::Velocity);
    display.report(s);

    cout << "q=" << s.getQ() << endl;
    cout << "qErr=" << s.getQErr() << endl;
    cout << "T_MbM=" << mobilizedBody.getMobilizerTransform(s).T() << endl;
    cout << "v_MbM=" << mobilizedBody.getMobilizerVelocity(s)[1] << endl;
    cout << "Unassembled configuration shown. Ready to assemble? "; cin >> ans;

    // These are the SimTK Simmath integrators:
    RungeKuttaMersonIntegrator myStudy(mbs);
    //CPodesIntegrator myStudy(mbs, CPodes::BDF, CPodes::Newton/*Functional*/);
    //myStudy.setOrderLimit(2); // cpodes only
    //VerletIntegrator myStudy(mbs);
    //ExplicitEulerIntegrator myStudy(mbs, .0005);


    //myStudy.setMaximumStepSize(0.001);
    myStudy.setAccuracy(1e-6);
    //myStudy.setProjectEveryStep(true);
    //myStudy.setProjectInterpolatedStates(false);
    myStudy.setConstraintTolerance(1e-7);
    //myStudy.setAllowInterpolation(false);
    //myStudy.setMaximumStepSize(.1);

    const Real dt = .02; // output intervals
    const Real finalTime = 2;

    myStudy.setFinalTime(finalTime);

    std::vector<State> saveEm;

    for (int i=0; i<50; ++i)
        saveEm.push_back(s);    // delay


    // Peforms assembly if constraints are violated.
    myStudy.initialize(s);

    for (int i=0; i<50; ++i)
        saveEm.push_back(s);    // delay

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
        cout << "T_MbM=" << mobilizedBody.getMobilizerTransform(s).T() << endl;
        cout << "PE=" << mbs.getPotentialEnergy(s) << " KE=" << mbs.getKineticEnergy(s) << " E=" << mbs.getEnergy(s) << endl;
        cout << "angle=" << std::acos(~mobilizedBody.expressBodyVectorInGround(s, Vec3(0,1,0)) * UnitVec3(1,1,1)) << endl;
        cout << "Assembled configuration shown. Ready to simulate? "; cin >> ans;
    }

    Integrator::SuccessfulStepStatus status;
    int nextReport = 0;

    mbs.resetAllCountersToZero();

    while ((status=myStudy.stepTo(nextReport*dt))
           != Integrator::EndOfSimulation) 
    {
        const State& s = myStudy.getState();
        mbs.realize(s, Stage::Acceleration);
        const Real angle = std::acos(~mobilizedBody.expressBodyVectorInGround(s, Vec3(0,1,0)) * UnitVec3(1,1,1));
        printf("%5g %10.4g E=%10.8g h%3d=%g %s%s\n", s.getTime(), 
            angle,
            mbs.getEnergy(s), myStudy.getNStepsTaken(),
            myStudy.getPreviousStepSizeTaken(),
            Integrator::successfulStepStatusString(status).c_str(),
            myStudy.isStateInterpolated()?" (INTERP)":"");
        printf("     qerr=%10.8g uerr=%10.8g uderr=%10.8g\n",
            matter.getQErr(s).normRMS(),
            matter.getUErr(s).normRMS(),
            s.getSystemStage() >= Stage::Acceleration ? matter.getUDotErr(s).normRMS() : Real(-1));
#ifdef HASC
		cout << "CONSTRAINT perr=" << c.getPositionError(s)
			 << " verr=" << c.getVelocityError(s)
			 << " aerr=" << c.getAccelerationError(s)
			 << endl;
#endif
		//cout << "   d(perrdot)/du=" << c.calcPositionConstraintMatrixP(s);
		//cout << "  ~d(f)/d lambda=" << c.calcPositionConstraintMatrixPT(s);
		//cout << "   d(perr)/dq=" << c.calcPositionConstraintMatrixPQInverse(s);
        cout << "Q=" << matter.getQ(s) << endl;

        Vector qdot;
        matter.calcQDot(s, s.getU(), qdot);
       // cout << "===> qdot =" << qdot << endl;

        Vector qdot2;
        matter.multiplyByQMatrix(s, false, s.getU(), qdot2);
       // cout << "===> qdot2=" << qdot2 << endl;

        Vector u1,u2;
        matter.multiplyByQMatrixInverse(s, false, qdot, u1);
        matter.multiplyByQMatrixInverse(s, false, qdot2, u2);
      //  cout << "===> u =" << s.getU() << endl;
      //  cout << "===> u1=" << u1 << endl;
      //  cout << "===> u2=" << u2 << endl;
       // cout << "     norm=" << (s.getU()-u2).normRMS() << endl;

        display.report(s);
        saveEm.push_back(s);

		if (status == Integrator::ReachedReportTime)
			++nextReport;
    }

    printf("Using Integrator %s:\n", myStudy.getMethodName());
    printf("# STEPS/ATTEMPTS = %d/%d\n", myStudy.getNStepsTaken(), myStudy.getNStepsAttempted());
    printf("# ERR TEST FAILS = %d\n", myStudy.getNErrorTestFailures());
    printf("# REALIZE/PROJECT = %d/%d\n", myStudy.getNRealizations(), myStudy.getNProjections());

    printf("System stats: realize %ldP %ldV %ldA, project %ld\n",
        mbs.getNumRealizationsOfThisStage(Stage::Position),
        mbs.getNumRealizationsOfThisStage(Stage::Velocity),
        mbs.getNumRealizationsOfThisStage(Stage::Acceleration),
        mbs.getNumProjectCalls());


    while(true) {
        for (int i=0; i < (int)saveEm.size(); ++i) {
            display.report(saveEm[i]);
            //display.report(saveEm[i]); // half speed
        }
        getchar();
    }
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

