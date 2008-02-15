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

class MyHomeMadeWeldImplementation;

class MyHomeMadeWeld : public Constraint::Custom {
public:
    MyHomeMadeWeld(MobilizedBody& body1, MobilizedBody& body2);
    MyHomeMadeWeld(MobilizedBody& body1, const Transform& frame1,
                   MobilizedBody& body2, const Transform& frame2);
};




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
    GeneralForceSubsystem    forces(mbs);
    DecorationSubsystem     viz(mbs);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -g, 0));

        // ADD BODIES AND THEIR MOBILIZERS
    const Body::Rigid body = Body::Rigid(MassProperties(m, Vec3(0), m*Inertia::brick(hl[0],hl[1],hl[2])))
                                  .addDecoration(Transform(), DecorativeBrick(hl).setOpacity(.5))
                                  .addDecoration(Transform(), DecorativeLine(Vec3(0), Vec3(0,1,0)).setColor(Green));

    MobilizedBody::Free mobilizedBody(matter.Ground(), Transform(), body, Transform());
    MobilizedBody::Free mobilizedBody0(mobilizedBody, Transform(Vec3(1,2,0)), body, Transform(Vec3(0,1,0)));
    //MobilizedBody::Ball mobilizedBody(mobilizedBody0, Transform(Vec3(1,2,0)), body, Transform(Vec3(0,1,0)));
    MobilizedBody::Free mobilizedBody2(mobilizedBody0, Vec3(-5,0,0), body, Transform());

    const Body::Rigid gear1body = Body::Rigid(MassProperties(m, Vec3(0), m*Inertia::cylinderAlongZ(.5, .1)))
        .addDecoration(Transform(), DecorativeCircle(.5).setColor(Green).setOpacity(.7))
        .addDecoration(Transform(), DecorativeLine(Vec3(0), Vec3(.5,0,0)).setColor(Black).setLineThickness(4));
    const Body::Rigid gear2body = Body::Rigid(MassProperties(m, Vec3(0), m*Inertia::cylinderAlongZ(1.5, .1)))
        .addDecoration(Transform(), DecorativeCircle(1.5).setColor(Blue).setOpacity(.7))  
        .addDecoration(Transform(), DecorativeLine(Vec3(0), Vec3(1.5,0,0)).setColor(Black).setLineThickness(4));
    MobilizedBody::Pin gear1(mobilizedBody0, Vec3(-1,0,0), gear1body, Transform()); // along z
    MobilizedBody::Pin gear2(mobilizedBody0, Vec3(1,0,0), gear2body, Transform()); // along z
    Constraint::NoSlip1D(mobilizedBody0, Vec3(-.5,0,0), UnitVec3(0,1,0), gear1, gear2);

#define HASC
    
    //Constraint::Ball myc2(matter.Ground(), Vec3(-4,2,0),  mobilizedBody2, Vec3(0,1,0));
    Constraint::Ball myc(matter.Ground(), Vec3(1,2,0),  mobilizedBody, Vec3(0,1,0));
    Constraint::Ball ball(mobilizedBody0, Vec3(2,0,0), mobilizedBody2, Vec3(3,0,0));

    //Constraint::ConstantOrientation ori(mobilizedBody, Rotation(), mobilizedBody2, Rotation());


    //Constraint::Weld weld(mobilizedBody, Transform(Rotation(Pi/4, ZAxis), Vec3(1,1,0)),
     //                     mobilizedBody2, Transform(Rotation(-Pi/4, ZAxis), Vec3(-1,-1,0)));

    MyHomeMadeWeld hweld(mobilizedBody, Transform(Rotation(Pi/4, ZAxis), Vec3(1,1,0)),
                          mobilizedBody2, Transform(Rotation(-Pi/4, ZAxis), Vec3(-1,-1,0)));
    
    //Constraint::PointInPlane pip(mobilizedBody, UnitVec3(0,1,0),  2, mobilizedBody2, Vec3(0,1,0));
    //pip.setPlaneDisplayHalfWidth(40);

    //Constraint::PointOnLine pol(mobilizedBody, UnitVec3(1,0,0),  Vec3(5,0,0), mobilizedBody2, Vec3(0,1,0));
    //pol.setLineDisplayHalfLength(5);

    //Constraint::Rod rod(mobilizedBody, hl, mobilizedBody2, -hl, 5);
    //Constraint::ConstantAngle ang(mobilizedBody0, UnitVec3(0,1,0),
     //                             mobilizedBody2, UnitVec3(0,1,0), Pi/4);

    
    //Constraint::Rod myc2(matter.Ground(), Vec3(-4,2,0),  mobilizedBody2, Vec3(0,1,0), 1);
    //Constraint::Rod myc(matter.Ground(), Vec3(1,2,0),  mobilizedBody, Vec3(0,1,0), 1);
    //Constraint::Rod ball(mobilizedBody0, Vec3(2,0,0), mobilizedBody2, Vec3(3,0,0), 5);
    

    //Constraint::Rod rod(mobilizedBody, Vec3(0,0,1), mobilizedBody2, Vec3(-1,0,0), 5);
    //Constraint::Rod myc(matter.Ground(), Vec3(1,2,0),  mobilizedBody, Vec3(0,1,0), 1);

    //Constraint::Ball myc(mobilizedBody, Vec3(1,2,0),  mobilizedBody2, Vec3(0,1,0));
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

    viz.addBodyFixedDecoration(mobilizedBody, Transform(), DecorativeCircle());
    viz.addBodyFixedDecoration(mobilizedBody, Transform(Vec3(1,2,3)), DecorativeText("hello world").setScale(.1));
    State s = mbs.realizeTopology(); // returns a reference to the the default state
    //matter.setUseEulerAngles(s, true);
    mbs.realizeModel(s); // define appropriate states for this System

    //mobilizedBody0.setQ(s, .1);
    //mobilizedBody.setQ(s, .2);

    VTKVisualizer display(mbs);

    //{Real qq[] = {0.9719,0,0,-0.235396,0.542437,1.11082,0};
    // s.updQ()=Vector(4,qq);
    //}
    mbs.realize(s, Stage::Velocity);
    display.report(s);
    cout << "q=" << s.getQ() << endl;
    cout << "u=" << s.getU() << endl;
    cout << "qErr=" << s.getQErr() << endl;
    cout << "uErr=" << s.getUErr() << endl;
#ifdef HASC
    for (ConstraintIndex cid(0); cid < matter.getNConstraints(); ++cid) {
        const Constraint& c = matter.getConstraint(cid);
        int mp,mv,ma;
        c.getNumConstraintEquationsInUse(s, mp,mv,ma);

	    cout << "CONSTRAINT " << cid 
             << " constrained bodies=" << c.getNumConstrainedBodies() 
             << " ancestor=" << c.getAncestorMobilizedBody().getMobilizedBodyIndex()
             << " constrained mobilizers/nq/nu=" << c.getNumConstrainedMobilizers() 
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
            cout << "perr=" << c.getPositionError(s) << endl;
	        cout << "   d(perrdot)/du=" << c.calcPositionConstraintMatrixP(s);
            cout << "  ~d(Pt lambda)/dlambda=" << ~c.calcPositionConstraintMatrixPt(s);
	        cout << "   d(perr)/dq=" << c.calcPositionConstraintMatrixPQInverse(s);

            Matrix P = c.calcPositionConstraintMatrixP(s);
            Matrix PQ(mp,matter.getNQ(s));
            Vector out(matter.getNQ(s));
            for (int i=0; i<mp; ++i) {
                Vector in = ~P[i];
                matter.multiplyByQMatrixInverse(s, true, in, out);
                PQ[i] = ~out;
            }
            cout << " calculated d(perr)/dq=" << PQ;
        }


        if (mv) {
            cout << "verr=" << c.getVelocityError(s) << endl;
	        //cout << "   d(verrdot)/dudot=" << c.calcVelocityConstraintMatrixV(s);
            cout << "  ~d(Vt lambda)/dlambda=" << ~c.calcVelocityConstraintMatrixVt(s);
        }

    }
    const Constraint& c = matter.getConstraint(myc.getConstraintIndex());
    
#endif
    Matrix Q(matter.getNQ(s), matter.getNU(s));
    Matrix Qinv(matter.getNU(s),matter.getNQ(s));
    Vector u(matter.getNU(s)); u=0;
    for (int i=0; i < matter.getNU(s); ++i) {
        u[i]=1;
        matter.multiplyByQMatrix(s,false,u,Q(i));
        u[i]=0;
    }
    cout << " BigQ=" << Q;

    char ans;
    cout << "Default configuration shown. Ready? "; cin >> ans;

    mobilizedBody.setQToFitTransform (s, Transform(Rotation(.05,Vec3(1,1,1)),Vec3(.1,.2,.3)));
    mobilizedBody0.setQToFitTransform (s, Transform(Rotation(.05,Vec3(1,-1,1)),Vec3(.2,.2,.3)));
    mobilizedBody2.setQToFitTransform (s, Transform(Rotation(.05,Vec3(-1,1,1)),Vec3(.1,.2,.1)));
    mobilizedBody.setUToFitAngularVelocity(s, 10*Vec3(.1,.2,.3));

    gear1.setUToFitAngularVelocity(s, Vec3(0,0,100)); // these should be opposite directions!
    //gear2.setUToFitAngularVelocity(s, Vec3(0,0,100));

    mbs.realize(s, Stage::Velocity);
    display.report(s);

    cout << "q=" << s.getQ() << endl;
    cout << "u=" << s.getU() << endl;
    cout << "qErr=" << s.getQErr() << endl;
    cout << "uErr=" << s.getUErr() << endl;
    cout << "T_MbM=" << mobilizedBody.getMobilizerTransform(s).T() << endl;
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
    myStudy.setAccuracy(1e-6); //myStudy.setAccuracy(1e-1);
    //myStudy.setProjectEveryStep(true);
    //myStudy.setProjectInterpolatedStates(false);
    myStudy.setConstraintTolerance(1e-7); //myStudy.setConstraintTolerance(1e-2);
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
        cout << "U=" << matter.getU(s) << endl;

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



// This is copied from the built-in Weld constraint and should behave identically.
class MyHomeMadeWeldImplementation : public Constraint::Custom::Implementation {
    static Real getDefaultAxisDisplayLength() {return 1;}
    static Vec3 getDefaultFrameColor(int which) {
        return which==0 ? Blue : Purple;
    }
public:
    MyHomeMadeWeldImplementation(MobilizedBody& body1, const MobilizedBody& body2)
      : Implementation(body1.updMatterSubsystem(), 6,0,0),  
        axisDisplayLength(-1), // means "use default axis length"
        frameBColor(-1), frameFColor(-1) // means "use default colors"
    {   // default Transforms are identity, i.e. body frames
        B = addConstrainedBody(body1);
        F = addConstrainedBody(body2);
    }

    MyHomeMadeWeldImplementation(MobilizedBody& body1, const Transform& f1, 
                                 const MobilizedBody& body2, const Transform& f2)
      : Implementation(body1.updMatterSubsystem(), 6,0,0), defaultFrameB(f1), defaultFrameF(f2), 
        axisDisplayLength(-1), // means "use default axis length"
        frameBColor(-1), frameFColor(-1) // means "use default colors"
    {   // default Transforms are identity, i.e. body frames
        B = addConstrainedBody(body1);
        F = addConstrainedBody(body2);
    }

    MyHomeMadeWeldImplementation* cloneVirtual() const {return new MyHomeMadeWeldImplementation(*this);}

    // Draw the two frames.
    void calcDecorativeGeometryAndAppendVirtual
       (const State& s, Stage stage, std::vector<DecorativeGeometry>& geom) const;

    void setAxisDisplayLength(Real len) {
        // len == 0 means "don't display"
        // len < 0 means "use default"
        invalidateTopologyCache();
        axisDisplayLength = len >= 0 ? len : -1;
    }
    Real getAxisDisplayLength() const {
        return axisDisplayLength < 0 ? getDefaultAxisDisplayLength() : axisDisplayLength;
    }

    void setFrameColor(int which, const Vec3& color) {
        assert(which==0 || which==1);
        // color[0] < 0 means "use default color for this frame"
        invalidateTopologyCache();
        if (which==0) frameBColor = color[0] < 0 ? Vec3(-1) : color;
        else          frameFColor = color[0] < 0 ? Vec3(-1) : color;
    }
    Vec3 getFrameColor(int which) const {
        assert(which==0 || which==1);
        if (which==0) return frameBColor[0] < 0 ? getDefaultFrameColor(0) : frameBColor;
        else          return frameFColor[0] < 0 ? getDefaultFrameColor(1) : frameFColor;
    }

    // Implementation of virtuals required for holonomic constraints.

    // For theory, look at the ConstantOrientation (1st 3 equations) and 
    // Ball (last 3 equations) theory above. Otherwise just lay back and 
    // enjoy the ride.

    void realizePositionErrorsVirtual(const State& s, int mp,  Real* perr) const {
        assert(mp==6 && perr);

        const Rotation& R_AB = getBodyRotation(s, B, true);
        const Rotation& R_AF = getBodyRotation(s, F, true);
        const Rotation  RB = R_AB * defaultFrameB.R(); // now expressed in A
        const Rotation  RF = R_AF * defaultFrameF.R();

        // Orientation error
        Vec3::updAs(perr) = Vec3(~RF.x()*RB.y(),
                                 ~RF.y()*RB.z(),
                                 ~RF.z()*RB.x());

        const Vec3 p_AF1 = calcStationLocation(s, B, defaultFrameB.T(), true);
        const Vec3 p_AF2 = calcStationLocation(s, F, defaultFrameF.T(), true);

        // position error
        Vec3::updAs(perr+3) = p_AF2 - p_AF1;
    }

    void realizePositionDotErrorsVirtual(const State& s, int mp,  Real* pverr) const {
        assert(mp==6 && pverr);
        //TODO: should be able to get p info from State
        const Rotation& R_AB = getBodyRotation(s, B);
        const Rotation& R_AF = getBodyRotation(s, F);
        const Rotation  RB = R_AB * defaultFrameB.R(); // now expressed in A
        const Rotation  RF = R_AF * defaultFrameF.R();

        const Vec3&     w_AB = getBodyAngularVelocity(s, B, true);
        const Vec3&     w_AF = getBodyAngularVelocity(s, F, true);
        const Vec3      w_BF = w_AF-w_AB; // in A

        // orientation error
        Vec3::updAs(pverr) = Vec3( ~w_BF * (RF.x() % RB.y()),
                                   ~w_BF * (RF.y() % RB.z()),
                                   ~w_BF * (RF.z() % RB.x()) );

        //TODO: should be able to get p info from State
        const Transform&  X_AB   = getBodyTransform(s, B);
        const Vec3        p_AF2  = calcStationLocation(s, F, defaultFrameF.T());
        const Vec3        p_BC   = ~X_AB*p_AF2; // C is a material point of body B

        const Vec3        v_AF2   = calcStationVelocity(s, F, defaultFrameF.T(), true);
        const Vec3        v_AC    = calcStationVelocity(s, B, p_BC, true);
 
        // position error
        Vec3::updAs(pverr+3) = v_AF2 - v_AC;
    }

    void realizePositionDotDotErrorsVirtual(const State& s, int mp,  Real* paerr) const {
        assert(mp==6 && paerr);
        //TODO: should be able to get p and v info from State
        const Rotation& R_AB = getBodyRotation(s, B);
        const Rotation& R_AF = getBodyRotation(s, F);
        const Rotation  RB = R_AB * defaultFrameB.R(); // now expressed in A
        const Rotation  RF = R_AF * defaultFrameF.R();

        const Vec3&     w_AB = getBodyAngularVelocity(s, B);
        const Vec3&     w_AF = getBodyAngularVelocity(s, F);
        const Vec3      w_BF = w_AF-w_AB; // in A

        const Vec3&     b_AB = getBodyAngularAcceleration(s, B, true);
        const Vec3&     b_AF = getBodyAngularAcceleration(s, F, true);
        const Vec3      b_BF = b_AF-b_AB; // in A

        // orientation error
        Vec3::updAs(paerr) = 
             Vec3( dot( b_BF, RF.x() % RB.y() )
                        + dot( w_BF, (w_AF%RF.x()) % RB.y() - (w_AB%RB.y()) % RF.x()),
                   dot( b_BF, RF.y() % RB.z() )
                        + dot( w_BF, (w_AF%RF.y()) % RB.z() - (w_AB%RB.z()) % RF.y()),
                   dot( b_BF, RF.z() % RB.x() )
                        + dot( w_BF, (w_AF%RF.z()) % RB.x() - (w_AB%RB.x()) % RF.z()));

        const Transform&  X_AB   = getBodyTransform(s, B);
        const Vec3        p_AF2  = calcStationLocation(s, F, defaultFrameF.T());
        const Vec3        p_BC   = ~X_AB*p_AF2; // C is a material point of body B

        const Vec3        a_AF2  = calcStationAcceleration(s, F, defaultFrameF.T(), true);
        const Vec3        a_AC   = calcStationAcceleration(s, B, p_BC, true);

        // position error
        Vec3::updAs(paerr+3) = a_AF2 - a_AC;
    }

    void applyPositionConstraintForcesVirtual
       (const State& s, int mp, const Real* multipliers,
        Vector_<SpatialVec>& bodyForcesInA,
        Vector&              mobilityForces) const
    {
        assert(mp==6 && multipliers);

        const Vec3& torques = Vec3::getAs(multipliers);
        const Vec3& force_A = Vec3::getAs(multipliers+3);

        //TODO: should be able to get p info from State
        const Rotation& R_AB = getBodyRotation(s, B);
        const Rotation& R_AF = getBodyRotation(s, F);
        const Rotation  RB = R_AB * defaultFrameB.R(); // now expressed in A
        const Rotation  RF = R_AF * defaultFrameF.R();

        const Vec3 torque_F_A =   torques[0] * (RF.x() % RB.y())
                                + torques[1] * (RF.y() % RB.z())
                                + torques[2] * (RF.z() % RB.x());

        addInBodyTorque(s, F,  torque_F_A, bodyForcesInA);
        addInBodyTorque(s, B, -torque_F_A, bodyForcesInA);

        const Transform& X_AB  = getBodyTransform(s,B);
        const Vec3&      p_FF2 = defaultFrameF.T();
        const Vec3       p_AF2 = calcStationLocation(s, F, p_FF2);
        const Vec3       p_BC = ~X_AB * p_AF2;

        addInStationForce(s, F, p_FF2, force_A, bodyForcesInA);
        addInStationForce(s, B, p_BC, -force_A, bodyForcesInA);
    }

private:
    ConstrainedBodyIndex B; // aka "body 1"
    ConstrainedBodyIndex F; // aka "body 2"

    Transform       defaultFrameB; // on body 1, relative to B frame
    Transform       defaultFrameF; // on body 2, relative to F frame};

    // These are for visualization control only.
    Real axisDisplayLength; // for all 6 axes; <= 0 means "don't display"
    Vec3 frameBColor, frameFColor;
};

MyHomeMadeWeld::MyHomeMadeWeld(MobilizedBody& body1, MobilizedBody& body2) 
  : Custom(new MyHomeMadeWeldImplementation(body1, body2))
{
}
MyHomeMadeWeld::MyHomeMadeWeld(MobilizedBody& body1, const Transform& frame1,
               MobilizedBody& body2, const Transform& frame2) 
  : Custom(new MyHomeMadeWeldImplementation(body1, frame1, body2, frame2))
{
}

void MyHomeMadeWeldImplementation::calcDecorativeGeometryAndAppendVirtual
   (const State& s, Stage stage, std::vector<DecorativeGeometry>& geom) const
{
    // We can't generate the frames until we know the axis lengths to use, and we can't place
    // the geometry on the bodies until we know the body1 and body2 frame
    // placements, which might not be until Instance stage.
    if (stage == Stage::Instance && getAxisDisplayLength() != 0 
        && getMatterSubsystem().getShowDefaultGeometry()) 
    {
        geom.push_back(DecorativeFrame(getAxisDisplayLength())
                                            .setColor(getFrameColor(0))
                                            .setLineThickness(2)
                                            .setBodyId(getMobilizedBodyIndexOfConstrainedBody(B))
                                            .setTransform(defaultFrameB));

        // Draw connector line back to body origin.
        if (defaultFrameB.T().norm() >= SignificantReal)
            geom.push_back(DecorativeLine(Vec3(0), defaultFrameB.T())
                             .setColor(getFrameColor(0))
                             .setLineThickness(2)
                             .setBodyId(getMobilizedBodyIndexOfConstrainedBody(B)));
           

        geom.push_back(DecorativeFrame(0.67*getAxisDisplayLength())
                                            .setColor(getFrameColor(1))
                                            .setLineThickness(4)
                                            .setBodyId(getMobilizedBodyIndexOfConstrainedBody(F))
                                            .setTransform(defaultFrameF));

        if (defaultFrameF.T().norm() >= SignificantReal)
            geom.push_back(DecorativeLine(Vec3(0), defaultFrameF.T())
                             .setColor(getFrameColor(1))
                             .setLineThickness(4)
                             .setBodyId(getMobilizedBodyIndexOfConstrainedBody(F)));
    }
}

