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
 * Throwaway main program.
 */

#include "SimTKsimbody.h"
#include "SimTKsimbody_aux.h" // requires VTK

#include "simbody/internal/VisualizationReporter.h"
#define VTKVisualizer Visualizer

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
static const Real g = 9.8; // meters/s^2; apply in �y direction
static const Real d = 0.5; // meters

static const Vec3 hl(1, 0.5, 0.5); // body half lengths

class MyConstraintImplementation : public Constraint::Custom::Implementation {
public:
    MyConstraintImplementation(MobilizedBody& mobilizer, Real speed)
      : Implementation(mobilizer.updMatterSubsystem(), 0,1,0),  
        theMobilizer(), whichMobility(), prescribedSpeed(NaN)
    {
        theMobilizer = addConstrainedMobilizer(mobilizer);
        whichMobility = MobilizerUIndex(0);
        prescribedSpeed = speed;
    }
    MyConstraintImplementation* clone() const {return new MyConstraintImplementation(*this);}

    // Implementation of virtuals required for nonholonomic constraints.

    // One non-holonomic (well, velocity-level) constraint equation.
    //    verr = u - s
    //    aerr = udot
    // 
    void realizeVelocityErrors(const State& s, int mv,  Real* verr) const {
        assert(mv==1 && verr);
        *verr = getOneU(s, theMobilizer, whichMobility) - prescribedSpeed;
    }

    void realizeVelocityDotErrors(const State& s, int mv,  Real* vaerr) const {
        assert(mv==1 && vaerr);
        *vaerr = getOneUDot(s, theMobilizer, whichMobility, true);
    }

    // apply generalized force lambda to the mobility
    void applyVelocityConstraintForces
       (const State& s, int mv, const Real* multipliers,
        Vector_<SpatialVec>& bodyForcesInA,
        Vector&              mobilityForces) const
    {
        assert(mv==1 && multipliers);
        const Real lambda = *multipliers;
        addInOneMobilityForce(s, theMobilizer, whichMobility, lambda, mobilityForces);
    }

private:
    ConstrainedMobilizerIndex   theMobilizer;
    MobilizerUIndex             whichMobility;
    Real                        prescribedSpeed;
};

class MyConstraint : public Constraint::Custom {
public:
    explicit MyConstraint(MobilizedBody& mobilizer, Real speed) 
      : Custom(new MyConstraintImplementation(mobilizer, speed))
    {
    }
};


int main(int argc, char** argv) {
  try { // If anything goes wrong, an exception will be thrown.

        // CREATE MULTIBODY SYSTEM AND ITS SUBSYSTEMS
    MultibodySystem         mbs;

    SimbodyMatterSubsystem  matter(mbs);
    GeneralForceSubsystem    forces(mbs);
    DecorationSubsystem     viz(mbs);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -g, 0));

        // ADD BODIES AND THEIR MOBILIZERS
    const Body::Rigid body = Body::Rigid(MassProperties(m, Vec3(0), m*Gyration::brick(hl[0],hl[1],hl[2])))
                                  .addDecoration(Transform(), DecorativeBrick(hl).setOpacity(.5))
                                  .addDecoration(Transform(), DecorativeLine(Vec3(0), Vec3(0,1,0)).setColor(Green));

    MobilizedBody::Free mobilizedBody(matter.Ground(), Transform(), body, Transform());
    MobilizedBody::Free mobilizedBody0(mobilizedBody, Transform(Vec3(1,2,0)), body, Transform(Vec3(0,1,0)));
    MobilizedBody::Free mobilizedBody2(mobilizedBody0, Vec3(-5,0,0), body, Transform());

    const Body::Rigid gear1body = Body::Rigid(MassProperties(m, Vec3(0), m*Gyration::cylinderAlongZ(.5, .1)))
        .addDecoration(Transform(), DecorativeCircle(.5).setColor(Green).setOpacity(.7))
        .addDecoration(Transform(), DecorativeLine(Vec3(0), Vec3(.5,0,0)).setColor(Black).setLineThickness(4));
    const Body::Rigid gear2body = Body::Rigid(MassProperties(m, Vec3(0), m*Gyration::cylinderAlongZ(1.5, .1)))
        .addDecoration(Transform(), DecorativeCircle(1.5).setColor(Blue).setOpacity(.7))  
        .addDecoration(Transform(), DecorativeLine(Vec3(0), Vec3(1.5,0,0)).setColor(Black).setLineThickness(4));

    MobilizedBody::Pin gear1(mobilizedBody2, Vec3(-1,0,0), gear1body, Transform()); // along z
    MobilizedBody::Pin gear2(mobilizedBody2, Vec3(1,0,0), gear2body, Transform()); // along z
    Constraint::NoSlip1D(mobilizedBody2, Vec3(-.5,0,0), UnitVec3(0,1,0), gear1, gear2);

    Constraint::ConstantSpeed(gear1, 100.);
    
    //Constraint::Ball myc2(matter.Ground(), Vec3(-4,2,0),  mobilizedBody2, Vec3(0,1,0));
    Constraint::Weld myc(matter.Ground(), Vec3(1,2,0),  mobilizedBody, Vec3(0,1,0));
    Constraint::Ball ball1(mobilizedBody, Vec3(2,0,0), mobilizedBody0, Vec3(3,1,1));
    Constraint::Ball ball2(mobilizedBody0, Vec3(2,0,0), mobilizedBody2, Vec3(3,0,0));
    //Constraint::PointInPlane pip(mobilizedBody0, UnitVec3(0,-1,0), 3, mobilizedBody2, Vec3(3,0,0));

    //Constraint::ConstantOrientation ori(mobilizedBody, Rotation(), mobilizedBody0, Rotation());
    //Constraint::ConstantOrientation ori2(mobilizedBody2, Rotation(), mobilizedBody0, Rotation());
    //Constraint::Weld weld(mobilizedBody, Transform(Rotation(Pi/4, ZAxis), Vec3(1,1,0)),
      //                    mobilizedBody2, Transform(Rotation(-Pi/4, ZAxis), Vec3(-1,-1,0)));
    
    //MyConstraint xyz(gear1, 100.);

    viz.addBodyFixedDecoration(mobilizedBody, Transform(Vec3(1,2,3)), DecorativeText("hello world").setScale(.1));



/*
    class MyHandler : public ScheduledEventHandler {
    public:
        MyHandler(const Constraint& cons) : c(cons) { }
        Real getNextEventTime(const State&, bool includeCurrentTime) const {
            return .314;
        }
        void handleEvent(State& s, Real acc, const Vector& ywts, const Vector& cwts, Stage& modified,
                         bool& shouldTerminate) const 
        {
            cout << "<<<< TRIGGERED AT T=" << s.getTime() << endl;
            c.enable(s);
            modified = Stage::Model;
        }
    private:
        const Constraint& c;
    };

    mbs.updDefaultSubsystem().addEventHandler(new MyHandler(xyz));
*/


    State s = mbs.realizeTopology(); // returns a reference to the the default state

    //xyz.disable(s);

    //matter.setUseEulerAngles(s, true);
    mbs.realizeModel(s); // define appropriate states for this System

    //mobilizedBody0.setQ(s, .1);
    //mobilizedBody.setQ(s, .2);

    VTKVisualizer display(mbs);

    mbs.realize(s, Stage::Velocity);
    display.report(s);
    cout << "q=" << s.getQ() << endl;
    cout << "u=" << s.getU() << endl;
    cout << "qErr=" << s.getQErr() << endl;
    cout << "uErr=" << s.getUErr() << endl;

    for (ConstraintIndex cid(0); cid < matter.getNumConstraints(); ++cid) {
        const Constraint& c = matter.getConstraint(cid);
        int mp,mv,ma;
        c.getNumConstraintEquationsInUse(s, mp,mv,ma);

        cout << "CONSTRAINT " << cid << (c.isDisabled(s) ? "**DISABLED** " : "")
             << " constrained bodies=" << c.getNumConstrainedBodies();
        if (c.getNumConstrainedBodies()) cout << " ancestor=" << c.getAncestorMobilizedBody().getMobilizedBodyIndex();
        cout << " constrained mobilizers/nq/nu=" << c.getNumConstrainedMobilizers() 
                                           << "/" << c.getNumConstrainedQ(s) << "/" << c.getNumConstrainedU(s)
             << " mp,mv,ma=" << mp << "," << mv << "," << ma 
             << endl;
        for (ConstrainedBodyIndex cid(0); cid < c.getNumConstrainedBodies(); ++cid) {
            cout << "  constrained body: " << c.getMobilizedBodyFromConstrainedBody(cid).getMobilizedBodyIndex(); 
            cout << endl;
        }
        for (ConstrainedMobilizerIndex cmx(0); cmx < c.getNumConstrainedMobilizers(); ++cmx) {
            cout << "  constrained mobilizer " << c.getMobilizedBodyFromConstrainedMobilizer(cmx).getMobilizedBodyIndex() 
                  << ", q(" << c.getNumConstrainedQ(s, cmx) << ")="; 
            for (MobilizerQIndex i(0); i < c.getNumConstrainedQ(s, cmx); ++i)
                cout << " " << c.getConstrainedQIndex(s, cmx, i);                  
            cout << ", u(" << c.getNumConstrainedU(s, cmx) << ")=";
            for (MobilizerUIndex i(0); i < c.getNumConstrainedU(s, cmx); ++i)
                cout << " " << c.getConstrainedUIndex(s, cmx, i);
            cout << endl;
        }
        cout << c.getSubtree();
             
        if (mp) {
            cout << "perr=" << c.getPositionErrorsAsVector(s) << endl;
            cout << "   d(perrdot)/du=" << c.calcPositionConstraintMatrixP(s);
            cout << "  ~d(Pt lambda)/dlambda=" << ~c.calcPositionConstraintMatrixPt(s);
            cout << "   d(perr)/dq=" << c.calcPositionConstraintMatrixPNInv(s);

            Matrix P = c.calcPositionConstraintMatrixP(s);
            Matrix PQ(mp,matter.getNQ(s));
            Vector out(matter.getNQ(s));
            for (int i=0; i<mp; ++i) {
                Vector in = ~P[i];
                matter.multiplyByNInv(s, true, in, out);
                PQ[i] = ~out;
            }
            cout << " calculated d(perr)/dq=" << PQ;
        }


        if (mv) {
            cout << "verr=" << c.getVelocityErrorsAsVector(s) << endl;
            //cout << "   d(verrdot)/dudot=" << c.calcVelocityConstraintMatrixV(s);
            cout << "  ~d(Vt lambda)/dlambda=" << ~c.calcVelocityConstraintMatrixVt(s);
        }

    }
    const Constraint& c = matter.getConstraint(myc.getConstraintIndex());

    char ans;
    cout << "Default configuration shown. Ready? "; cin >> ans;

    mobilizedBody.setQToFitTransform (s, Transform(Rotation(.05,Vec3(1,1,1)),Vec3(.1,.2,.3)));
    mobilizedBody0.setQToFitTransform (s, Transform(Rotation(.05,Vec3(1,-1,1)),Vec3(.2,.2,.3)));
    mobilizedBody2.setQToFitTransform (s, Transform(Rotation(.05,Vec3(-1,1,1)),Vec3(.1,.2,.1)));
    mobilizedBody.setUToFitAngularVelocity(s, 10*Vec3(.1,.2,.3));
    mobilizedBody0.setUToFitAngularVelocity(s, 10*Vec3(.1,.2,.3));
    mobilizedBody2.setUToFitAngularVelocity(s, 10*Vec3(.1,.2,.3));

    //gear1.setUToFitAngularVelocity(s, Vec3(0,0,500)); // these should be opposite directions!
    //gear2.setUToFitAngularVelocity(s, Vec3(0,0,100));

    mbs.realize(s, Stage::Velocity);
    display.report(s);

    cout << "q=" << s.getQ() << endl;
    cout << "u=" << s.getU() << endl;
    cout << "qErr=" << s.getQErr() << endl;
    cout << "uErr=" << s.getUErr() << endl;
    cout << "p_MbM=" << mobilizedBody.getMobilizerTransform(s).p() << endl;
    cout << "v_MbM=" << mobilizedBody.getMobilizerVelocity(s)[1] << endl;
    cout << "Unassembled configuration shown. Ready to assemble? "; cin >> ans;

    // These are the SimTK Simmath integrators:
    RungeKuttaMersonIntegrator myStudy(mbs);
    //CPodesIntegrator myStudy(mbs, CPodes::BDF, CPodes::/*Newton*/Functional);
    //myStudy.setOrderLimit(2); // cpodes only
    //VerletIntegrator myStudy(mbs);
   // ExplicitEulerIntegrator myStudy(mbs, .0005); // fixed step
    //ExplicitEulerIntegrator myStudy(mbs); // variable step


    //myStudy.setMaximumStepSize(0.001);
    myStudy.setAccuracy(1e-6); myStudy.setAccuracy(1e-1);
    //myStudy.setProjectEveryStep(true);
    //myStudy.setProjectInterpolatedStates(false);
    myStudy.setConstraintTolerance(1e-7); myStudy.setConstraintTolerance(1e-2);
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
        cout << "u=" << s.getU() << endl;
        cout << "qErr=" << s.getQErr() << endl;
        cout << "uErr=" << s.getUErr() << endl;
        cout << "p_MbM=" << mobilizedBody.getMobilizerTransform(s).p() << endl;
        cout << "PE=" << mbs.calcPotentialEnergy(s) << " KE=" << mbs.calcKineticEnergy(s) << " E=" << mbs.calcEnergy(s) << endl;
        cout << "angle=" << std::acos(~mobilizedBody.expressVectorInGroundFrame(s, Vec3(0,1,0)) * UnitVec3(1,1,1)) << endl;
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
        const Real angle = std::acos(~mobilizedBody.expressVectorInGroundFrame(s, Vec3(0,1,0)) * UnitVec3(1,1,1));
        printf("%5g %10.4g E=%10.8g h%3d=%g %s%s\n", s.getTime(), 
            angle,
            mbs.calcEnergy(s), myStudy.getNumStepsTaken(),
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
        cout << "U=" << matter.getU(s) << endl;
        cout << "Multipliers=" << matter.getMultipliers(s) << endl;

        Vector qdot;
        matter.calcQDot(s, s.getU(), qdot);
       // cout << "===> qdot =" << qdot << endl;

        Vector qdot2;
        matter.multiplyByN(s, false, s.getU(), qdot2);
       // cout << "===> qdot2=" << qdot2 << endl;

        Vector u1,u2;
        matter.multiplyByNInv(s, false, qdot, u1);
        matter.multiplyByNInv(s, false, qdot2, u2);
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
    printf("# STEPS/ATTEMPTS = %d/%d\n", myStudy.getNumStepsTaken(), myStudy.getNumStepsAttempted());
    printf("# ERR TEST FAILS = %d\n", myStudy.getNumErrorTestFailures());
    printf("# REALIZE/PROJECT = %d/%d\n", myStudy.getNumRealizations(), myStudy.getNumProjections());

    printf("System stats: realize %dP %dV %dA, project %d\n",
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

