/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-12 Stanford University and the Authors.        *
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

/**@file
 * This is a simple example of using a constraint.
 * Here we have two independent pendulums hanging from ground pins.
 * They can be connected by a spring or a distance constraint.
 */

#include "SimTKsimbody.h"

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

class ShermsForce : public Force::Custom::Implementation {
public:
    ShermsForce(const MobilizedBody& b1, const MobilizedBody& b2) : body1(b1), body2(b2) { }
    ShermsForce* clone() const {return new ShermsForce(*this);}

    void calcForce(const State& state,
              Vector_<SpatialVec>& bodyForces,
              Vector_<Vec3>&       particleForces,
              Vector&              mobilityForces) const override
    {
        const Vec3& pos1 = body1.getBodyTransform(state).p();
        const Vec3& pos2 = body2.getBodyTransform(state).p();
        const Real d = (pos2-pos1).norm();
        const Real k = 1000, d0 = 1;
        const Vec3 f = k*(d-d0)*(pos2-pos1)/d;
        body1.applyBodyForce(state, SpatialVec(Vec3(0),  f), bodyForces);
        body2.applyBodyForce(state, SpatialVec(Vec3(0), -f), bodyForces);
    }
    Real calcPotentialEnergy(const State& state) const override {
        return 0;
    }
private:
    const MobilizedBody& body1;
    const MobilizedBody& body2;
};
//template <class E> Vector_<E>
//operator*(const VectorView_<E>& l, const typename CNT<E>::StdNumber& r)
//  { return Vector_<E>(l)*=r; }

void ff(Vector& v) {
    v = 23.;
}

int main(int argc, char** argv) {
  try { // If anything goes wrong, an exception will be thrown.

        // CREATE MULTIBODY SYSTEM AND ITS SUBSYSTEMS
    MultibodySystem         mbs;

    SimbodyMatterSubsystem  twoPends(mbs);
    GeneralForceSubsystem    forces(mbs);
    DecorationSubsystem     viz(mbs);
    Force::UniformGravity gravity(forces, twoPends, Vec3(0, -g, 0));
    gravity.setDisabledByDefault(true);

        // ADD BODIES AND THEIR MOBILIZERS
    Body::Rigid pendulumBody = Body::Rigid(MassProperties(m, Vec3(0), Inertia(1)));
    pendulumBody.addDecoration(Transform(),
                               DecorativeBrick(Vec3(.1,.0667,.05)).setOpacity(.5));

    MobilizedBody:: Ball /*Gimbal*/ /*FreeLine*/ /*LineOrientation*/ /*Free*/
        leftPendulum(twoPends.Ground(),
                         Transform(Vec3(-1, 0, 0)),
                     pendulumBody,
                         Transform(Vec3(0, d, 0)));
/*
    MobilizedBody::Ball
        leftPendulum2(leftPendulum,
                         Transform(Vec3(0.5, 0, 0)),
                     pendulumBody,
                         Transform(Vec3(0, d, 0)));
*/

//    leftPendulum.setDefaultRadius(0.2); // for Ball artwork

    Vec3 radii(1.5/2.,1/3.,1/4.); radii*=.5; //radii=Vec3(.333,.5,1);
    MobilizedBody::Ellipsoid rightPendulum(twoPends.Ground(), pendulumBody);
    rightPendulum.setDefaultRadii(radii)
        .setDefaultInboardFrame(Transform(Rotation(),Vec3(1,0,0)))
        .setDefaultOutboardFrame(Transform( Rotation( SpaceRotationSequence, Pi/2, XAxis, -Pi/2, YAxis ), Vec3(0,d-radii[1],0)));

    //rightPendulum.setDefaultAngle(20*Deg2Rad);
   // rightPendulum.setDefaultRotation( Rotation(60*Deg2Rad, Vec3(0,0,1)) );

        // OPTIONALLY TIE TOGETHER WITH SPRING/DAMPER OR DISTANCE CONSTRAINT

    const Real distance = /*2*/1.5;      // nominal length for spring; length for constraint
    const Real stiffness = 100;   // only if spring is used
    const Real damping   = 10;     //          "

    char c;
    cout << "Constraint, spring, or nothing? c/s/n"; cin >> c;

    ConstraintIndex cid;
    const Vec3 leftAttachPt(.1,0.05,0);
    if (c == 'c') {

        cid =
        //Constraint::PointInPlane(twoPends.Ground(), UnitVec3(0,1,0), -2*d,
        //                         leftPendulum2, Vec3(0))
        Constraint::Rod(leftPendulum, leftAttachPt,
                        rightPendulum, Vec3(0),
                       distance)
        // Constraint::Ball(leftPendulum2, Vec3(.5,0,0),
        //                 twoPends.Ground(), Vec3(0,-d,0))
        .getConstraintIndex();

    } else if (c == 's') {
        Force::TwoPointLinearSpring(forces, leftPendulum, leftAttachPt,
                                       rightPendulum, Vec3(0),
                                       stiffness, distance);
        Force::TwoPointLinearDamper(forces, leftPendulum, leftAttachPt,
                                       rightPendulum, Vec3(0),
                                       damping);
    }

    // Add visualization line for spring (Rod constraint has one automatically)
    if (c=='s')
        viz.addRubberBandLine(leftPendulum, leftAttachPt,
                              rightPendulum, Vec3(0),
                              DecorativeLine().setColor(Orange).setLineThickness(4));

    //forces.addMobilityConstantForce(leftPendulum, 0, 20);
    //forces.addCustomForce(ShermsForce(leftPendulum,rightPendulum));
    //forces.addGlobalEnergyDrain(1);

    mbs.setHasTimeAdvancedEvents(false);

    cout << "HAS TIME ADVANCED EVENTS=" << mbs.hasTimeAdvancedEvents() << endl;

    Measure::Constant meas1(twoPends, 20);
    const Real amp = 3, freq = 100, phase = Pi/2;
    Measure::Sinusoid sint(twoPends, amp, freq, phase);

    Measure::Integrate twentyPlus10t(twoPends, Measure::Constant(twoPends, 10), meas1);

    State s = mbs.realizeTopology(); // returns a reference to the the default state
    //twoPends.setUseEulerAngles(s, true);
    mbs.realizeModel(s); // define appropriate states for this System

    gravity.enable(s);
    twentyPlus10t.setValue(s, 20);

    mbs.realize(s, Stage::Instance); // instantiate constraints

    cout << "meas1=" << meas1.getValue(s) << endl;


    if (cid.isValid()) {
        int mp, mv, ma;
        twoPends.getConstraint(cid).getNumConstraintEquationsInUse(s, mp, mv, ma);
        cout << "CONSTRAINT ID " << cid << " mp,v,a=" << mp << ", " << mv << ", " << ma << endl;
        cout << "CONSTRAINT -- " << twoPends.getConstraint(cid).getSubtree();
    }

    for (MobilizedBodyIndex i(0); i < twoPends.getNumBodies(); ++i) {
        const MobilizedBody& mb = twoPends.getMobilizedBody(i);
        cout << "Body " << i
             << " base=" << mb.getBaseMobilizedBody().getMobilizedBodyIndex()
             << endl;
    }


    SimbodyMatterSubtree sub(twoPends);
    sub.addTerminalBody(leftPendulum); sub.addTerminalBody(rightPendulum);
    sub.realizeTopology();
    cout << "SUB -- " << sub;

    SimbodyMatterSubtreeResults results;
    sub.initializeSubtreeResults(s, results);
    cout << "INIT RESULTS=" << results;

    Visualizer display(mbs);
    display.setBackgroundType(Visualizer::SolidColor);

    // gravity.disable(s);
    mbs.realize(s, Stage::Position);
    display.report(s);
    cout << "q=" << s.getQ() << endl;
    cout << "qErr=" << s.getQErr() << endl;
    cout << "p_MbM=" << rightPendulum.getMobilizerTransform(s).p() << endl;

    Vector_<SpatialVec> bodyForces;
    Vector_<Vec3> particleForces;
    Vector mobilityForces;
    gravity.calcForceContribution(s,bodyForces,particleForces,mobilityForces);
    cout << "Gravity forces: body:" << bodyForces << endl;
    cout << "                particle:" << particleForces << endl;
    cout << "                mobility:" << mobilityForces << endl;
    cout << "  PE=" << gravity.calcPotentialEnergyContribution(s) << endl;

    if (cid.isValid()) {
        const Constraint& c = twoPends.getConstraint(cid);
        cout << "CONSTRAINT perr=" << c.getPositionErrorsAsVector(s)
             << endl;
        cout << "   d(perrdot)/du=" << c.calcPositionConstraintMatrixP(s);
        cout << "   d(perr)/dq=" << c.calcPositionConstraintMatrixPNInv(s);
    }

    cout << "Default configuration shown. Ready? "; cin >> c;

    sub.copyPositionsFromState(s, results);
    cout << "POS RESULTS=" << results;

    //leftPendulum.setAngle(s, -60*Deg2Rad);
    //leftPendulum.setQToFitRotation(s, Rotation(-60*Deg2Rad,ZAxis));

    //rightPendulum.setQToFitTranslation(s, Vec3(0,1,0));
    leftPendulum.setQToFitRotation (s, Rotation(-.9*Pi/2,ZAxis));
    rightPendulum.setQToFitRotation(s, Rotation(-.9*Pi/2,YAxis));


    //TODO
    //rightPendulum.setUToFitLinearVelocity(s, Vec3(1.1,0,1.2));

    leftPendulum.setUToFitAngularVelocity(s, 10*Vec3(.1,.2,.3));
    rightPendulum.setUToFitAngularVelocity(s, 10*Vec3(.1,.2,.3));


    s.setTime(0);

    mbs.realize(s, Stage::Velocity);
    display.report(s);

    cout << "q=" << s.getQ() << endl;
    cout << "qErr=" << s.getQErr() << endl;
    cout << "p_MbM=" << rightPendulum.getMobilizerTransform(s).p() << endl;
    cout << "v_MbM=" << rightPendulum.getMobilizerVelocity(s)[1] << endl;
    cout << "Unassembled configuration shown. Ready to assemble? "; cin >> c;

    // These are the SimTK Simmath integrators:
    RungeKuttaMersonIntegrator myStudy(mbs);
    //CPodesIntegrator myStudy(mbs, CPodes::BDF, CPodes::Newton);
    //myStudy.setOrderLimit(2); // cpodes only
    //VerletIntegrator myStudy(mbs);


    //myStudy.setMaximumStepSize(0.001);
    myStudy.setAccuracy(1e-1);
    //myStudy.setProjectEveryStep(true);
    myStudy.setConstraintTolerance(1e-2);
    //myStudy.setAllowInterpolation(false);
    //myStudy.setMaximumStepSize(.1);

    const Real dt = 1./60; // output intervals
    const Real finalTime = 20;

    myStudy.setFinalTime(finalTime);

    // Peforms assembly if constraints are violated.
    myStudy.initialize(s);

    cout << "Using Integrator " << std::string(myStudy.getMethodName()) << ":\n";
    cout << "ACCURACY IN USE=" << myStudy.getAccuracyInUse() << endl;
    cout << "CTOL IN USE=" << myStudy.getConstraintToleranceInUse() << endl;
    cout << "TIMESCALE=" << mbs.getDefaultTimeScale() << endl;
    cout << "U WEIGHTS=" << s.getUWeights() << endl;
    cout << "Z WEIGHTS=" << s.getZWeights() << endl;
    cout << "1/QTOLS=" << s.getQErrWeights() << endl;
    cout << "1/UTOLS=" << s.getUErrWeights() << endl;

    {
        const State& s = myStudy.getState();
        display.report(s);
        cout << "q=" << s.getQ() << endl;
        cout << "qErr=" << s.getQErr() << endl;
        cout << "p_MbM=" << rightPendulum.getMobilizerTransform(s).p() << endl;
        cout << "Assembled configuration shown. Ready to simulate? "; cin >> c;
    }

    Integrator::SuccessfulStepStatus status;
    int nextReport = 0;
    int nextScheduledEvent = 0;
    Real schedule[] = {1.234, 3.1415, 3.14159, 4.5, 9.090909, 100.};
    while ((status=myStudy.stepTo(nextReport*dt,schedule[nextScheduledEvent]))
           != Integrator::EndOfSimulation)
    {
        const State& s = myStudy.getState();
        mbs.realize(s);
        const Real leftPendulumAngle = leftPendulum.getBodyRotation(s).convertRotationToAngleAxis()[0] * Rad2Deg;

        if (status == Integrator::ReachedScheduledEvent
            || std::abs(std::floor(s.getTime()+0.5)-s.getTime())<1e-4)
        {
            printf("%5g %10.4g E=%10.8g h%3d=%g %s%s\n", s.getTime(),
                leftPendulumAngle,
                mbs.calcEnergy(s), myStudy.getNumStepsTaken(),
                myStudy.getPreviousStepSizeTaken(),
                Integrator::getSuccessfulStepStatusString(status).c_str(),
                myStudy.isStateInterpolated()?" (INTERP)":"");
            printf("     qerr=%10.8g uerr=%10.8g uderr=%10.8g\n",
                twoPends.getQErr(s).normRMS(),
                twoPends.getUErr(s).normRMS(),
                twoPends.getUDotErr(s).normRMS());

            cout << "t=" << s.getTime() << "sint=" << sint.getValue(s) << "a*sin(wt+p)="
                << amp*std::sin(freq*s.getTime() + phase) << endl;

            cout << "20+10t=" << twentyPlus10t.getValue(s) << endl;

            if (cid.isValid()) {
                const Constraint& c = twoPends.getConstraint(cid);
                cout << "CONSTRAINT perr=" << c.getPositionErrorsAsVector(s)
                     << " verr=" << c.getVelocityErrorsAsVector(s)
                     << " aerr=" << c.getAccelerationErrorsAsVector(s)
                     << endl;
                //cout << "   d(perrdot)/du=" << c.calcPositionConstraintMatrixP(s);
                //cout << "  ~d(f)/d lambda=" << c.calcPositionConstraintMatrixPT(s);
                //cout << "   d(perr)/dq=" << c.calcPositionConstraintMatrixPQInverse(s);
            }
        }

        Vector qdot;
        twoPends.calcQDot(s, s.getU(), qdot);
       // cout << "===> qdot =" << qdot << endl;

        Vector qdot2;
        twoPends.multiplyByN(s, false, s.getU(), qdot2);
       // cout << "===> qdot2=" << qdot2 << endl;

        Vector u1,u2;
        twoPends.multiplyByNInv(s, false, qdot, u1);
        twoPends.multiplyByNInv(s, false, qdot2, u2);
      //  cout << "===> u =" << s.getU() << endl;
      //  cout << "===> u1=" << u1 << endl;
      //  cout << "===> u2=" << u2 << endl;
       // cout << "     norm=" << (s.getU()-u2).normRMS() << endl;

        //sub.copyPositionsFromState(s, results);
        //sub.copyVelocitiesFromState(s, results);
       // sub.copyAccelerationsFromState(s, results);
        //cout << results;

        display.report(s);
        //if (s.getTime() >= finalTime)
           // break;

        //status = myStudy.stepTo(s.getTime() + dt);

        //THIS CAN FAIL SOMETIMES
        //if (s.getTime() >= nextReport*dt)
        //    ++nextReport;

        if (status == Integrator::ReachedReportTime)
            ++nextReport;

        if (s.getTime() >= schedule[nextScheduledEvent])
            ++nextScheduledEvent;
    }

    printf("Using Integrator %s:\n", myStudy.getMethodName());
    printf("# STEPS/ATTEMPTS = %d/%d\n", myStudy.getNumStepsTaken(), myStudy.getNumStepsAttempted());
    printf("# ERR TEST FAILS = %d\n", myStudy.getNumErrorTestFailures());
    printf("# REALIZE/PROJECT = %d/%d\n", myStudy.getNumRealizations(), myStudy.getNumProjections());

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

