/* -------------------------------------------------------------------------- *
 *                           SimTK Simbody(tm)                                *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-12 Stanford University and the Authors.        *
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
 * A big, chunky fake RNA built of cylinders and ball joints.
 */

#include "SimTKsimbody.h"

#include <cmath>
#include <cstdio>
#include <exception>
#include <vector>
#include <ctime>

#ifdef _MSC_VER
#pragma warning(disable:4996) // don't warn about strerror, sprintf, etc.
#endif

using namespace std;
using namespace SimTK;

static const Real Deg2Rad = (Real)SimTK_DEGREE_TO_RADIAN;
static const int  GroundBodyNum = 0; // ground is always body 0

static const Real g = 9.8; // meters/s^2; apply in –y direction

static const Real DuplexRadius = 3; // A
static const Real HalfHeight = 10;  // A
static const Real CylinderSlop = 1; // A

static const int  NAtoms = 20;
static const Real AtomMass = 12;    // Daltons
static const Real AtomRadius = 1;   // A

static const Real ConnectorRadius     = 1;  // A
static const Real ConnectorHalfHeight = 3;  // A
static const Real ConnectorEndSlop    = 0.2;// A
static const Real ConnectorDensity    = 10;  // Dalton/A^3

static int NSegments = 3;

class MyRNAExample : public SimbodyMatterSubsystem {
    struct PerBodyInfo {
        PerBodyInfo(MobilizedBodyIndex b, bool d) : bnum(b), isDuplex(d) { }
        MobilizedBodyIndex bnum;
        bool isDuplex;
    };
    std::vector<PerBodyInfo> bodyInfo;
    MobilizedBodyIndex end1, end2;
public:
    MyRNAExample(MultibodySystem& mbs, int nsegs, bool shouldFlop) : SimbodyMatterSubsystem(mbs)
    {
        bodyInfo.push_back(PerBodyInfo(GroundIndex, false)); // placeholder for ground
        end1 = makeChain(GroundIndex, Vec3(0), nsegs, shouldFlop);
        end2 = makeChain(GroundIndex, Vec3(20,0,0), nsegs, shouldFlop);

        if (true) {
            Constraint::Rod theConstraint2(updMobilizedBody(end1), Vec3(0, -HalfHeight,0),
                                           updMobilizedBody(end2), Vec3(0, -HalfHeight,0), 10);
        }

    }

    void decorateBody(MobilizedBodyIndex bodyNum, Visualizer& display) const {
        assert(bodyInfo[bodyNum].bnum == bodyNum);
        if (bodyInfo[bodyNum].isDuplex)
            addDuplexDecorations(bodyNum, DuplexRadius, HalfHeight, CylinderSlop,
                                 NAtoms, AtomRadius, display);
        else
            addConnectorDecorations(bodyNum, ConnectorRadius, ConnectorHalfHeight,
                                    ConnectorEndSlop, display);
    }

    void decorateGlobal(Visualizer& display) const {
        // nothing
    }

private:

    MobilizedBodyIndex makeChain(MobilizedBodyIndex startBodyId, const Vec3& startOrigin, int nSegs, bool shouldFlop) {
        //MobilizedBodyIndex baseBodyIx = startBody;
        MobilizedBody baseBody = updMobilizedBody(startBodyId);
        Vec3 origin = startOrigin;
        MobilizedBody lastDup;
        for (int seg=0; seg < nSegs; ++seg) {

            MobilizedBody::Ball left1(
                baseBody, Transform(origin + Vec3(-DuplexRadius,-HalfHeight,0)),
                Body::Rigid(calcConnectorMassProps(ConnectorRadius, ConnectorHalfHeight, ConnectorDensity)),
                Transform(Vec3(0, ConnectorHalfHeight, 0)));
            left1.setDefaultRadius(1.5);
            bodyInfo.push_back(PerBodyInfo(left1, false));

            MobilizedBody::Ball left2(
                left1, Transform(Vec3(0, -ConnectorHalfHeight, 0)),
                Body::Rigid(calcConnectorMassProps(ConnectorRadius, ConnectorHalfHeight, ConnectorDensity)),
                Transform(Vec3(0, ConnectorHalfHeight, 0)));
            left2.setDefaultRadius(1.5);
            bodyInfo.push_back(PerBodyInfo(left2, false));

            MobilizedBody::Ball rt1(
                baseBody, Transform(origin + Vec3(DuplexRadius,-HalfHeight,0)),
                Body::Rigid(calcConnectorMassProps(ConnectorRadius, ConnectorHalfHeight, ConnectorDensity)),
                Transform(Vec3(0, ConnectorHalfHeight, 0)));
            rt1.setDefaultRadius(1.5);
            bodyInfo.push_back(PerBodyInfo(rt1, false));

            MobilizedBody::Ball rt2(
                rt1, Transform(Vec3(0, -ConnectorHalfHeight, 0)),
                Body::Rigid(calcConnectorMassProps(ConnectorRadius, ConnectorHalfHeight, ConnectorDensity)),
                Transform(Vec3(0, ConnectorHalfHeight, 0)));
            rt2.setDefaultRadius(1.5);
            bodyInfo.push_back(PerBodyInfo(rt2, false));

            MobilizedBody::Ball dup(
                rt2, Transform(Vec3(0, -ConnectorHalfHeight, 0)),
                Body::Rigid(calcDuplexMassProps(DuplexRadius, HalfHeight, NAtoms, AtomMass)),
                                Transform(Vec3(-DuplexRadius, HalfHeight, 0)));
            dup.setDefaultRadius(1.5);
            bodyInfo.push_back(PerBodyInfo(dup, true));

            if (!shouldFlop) {
                Constraint::Ball theConstraint(left2, Vec3(0, -ConnectorHalfHeight, 0),
                                               dup, Vec3(DuplexRadius, HalfHeight, 0));
                theConstraint.setDefaultRadius(1.5);
            }

            baseBody = lastDup = dup;
            origin = Vec3(0);
        }
        return lastDup;
    }

    MassProperties calcDuplexMassProps(
        Real halfHeight, Real r, int nAtoms, Real atomMass)
    {
        const Real pitch = 2*Pi/halfHeight;
        const Real trans = (2*halfHeight)/(nAtoms-1);
        const Real rot = pitch*trans;
        Inertia iner(0);
        Vec3 com(0);
        Real mass = 0;
        for (int i=0; i<nAtoms; ++i) {
            const Real h = halfHeight - i*trans;
            const Real th = i*rot;
            const Vec3 p1(-r*cos(th),h,r*sin(th)), p2(r*cos(th),h,-r*sin(th));
            mass += 2*atomMass;
            iner += Inertia(p1, atomMass) + Inertia(p2, atomMass);
            com += atomMass*p1 + atomMass*p2;
        }
        return MassProperties(mass,com/mass,iner);
    }

    MassProperties calcConnectorMassProps(Real r, Real halfHeight, Real density)
    {
        const Real volume = Pi*r*r*halfHeight;
        const Real mass = volume*density;
        const Vec3 com = Vec3(0);
        const Inertia iner = mass*UnitInertia::cylinderAlongY(r, halfHeight);

        return MassProperties(mass,com,iner);
    }

    void addDuplexDecorations(MobilizedBodyIndex bodyNum, Real r, Real halfHeight, Real slop, int nAtoms,
                              Real atomRadius, Visualizer& display) const
    {
        display.addDecoration(bodyNum, Transform(),
            DecorativeCylinder(r+atomRadius+slop, halfHeight).setColor(Cyan).setOpacity(0.4));

        const Real pitch = 2*Pi/halfHeight;
        const Real trans = (2*halfHeight)/(nAtoms-1);
        const Real rot = pitch*trans;
        for (int i=0; i<nAtoms; ++i) {
            const Real h = halfHeight - i*trans;
            const Real th = i*rot;
            const Vec3 p1(-r*cos(th),h,r*sin(th)), p2(r*cos(th),h,-r*sin(th));
            display.addDecoration(bodyNum, Transform(Vec3(p1)),
                DecorativeSphere(atomRadius).setColor(Red).setResolution(0.5));
            display.addDecoration(bodyNum, Transform(Vec3(p2)),
                DecorativeSphere(atomRadius).setColor(Green).setResolution(0.5));
        }
    }

    void addConnectorDecorations(MobilizedBodyIndex bodyNum, Real r, Real halfHeight, Real endSlop,
                                 Visualizer& display) const
    {
        display.addDecoration(bodyNum, Transform(),
            DecorativeCylinder(r, halfHeight-endSlop).setColor(Blue));
    }
};


int main(int argc, char** argv) {
    std::vector<State> saveEm;
    saveEm.reserve(1000);

try // If anything goes wrong, an exception will be thrown.
  { int nseg = NSegments;
    int shouldFlop = 0;
    if (argc > 1) sscanf(argv[1], "%d", &nseg);
    if (argc > 2) sscanf(argv[2], "%d", &shouldFlop);

    // Create a multibody system using Simbody.
    MultibodySystem mbs;
    MyRNAExample myRNA(mbs, nseg, shouldFlop != 0);
    GeneralForceSubsystem forces(mbs);
    Force::UniformGravity ugs(forces, myRNA, Vec3(0, -g, 0), -0.8);

    const Vec3 attachPt(150, -40, -50);

    Force::TwoPointLinearSpring(forces, myRNA.Ground(), attachPt,
                                   myRNA.getMobilizedBody(MobilizedBodyIndex(myRNA.getNumBodies()-1)), Vec3(0),
                                   1000.,  // stiffness
                                   1.);    // natural length

    Force::GlobalDamper(forces, myRNA, 1000);


    State s = mbs.realizeTopology();


    //myRNA.getConstraint(ConstraintIndex(myRNA.getNConstraints()-4)).disable(s);


    //myRNA.setUseEulerAngles(s,true);
    mbs.realizeModel(s);
    mbs.realize(s, Stage::Position);

    for (ConstraintIndex cid(0); cid < myRNA.getNumConstraints(); ++cid) {
        const Constraint& c = myRNA.getConstraint(cid);
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

        cout << "   d(perrdot)/du=" << c.calcPositionConstraintMatrixP(s);
        cout << "   d(perrdot)/du=" << ~c.calcPositionConstraintMatrixPt(s);

        cout << "   d(perr)/dq=" << c.calcPositionConstraintMatrixPNInv(s);
    }





    SimbodyMatterSubtree sub(myRNA);
    sub.addTerminalBody(myRNA.getMobilizedBody(MobilizedBodyIndex(7)));
    sub.addTerminalBody(myRNA.getMobilizedBody(MobilizedBodyIndex(10)));
    //sub.addTerminalBody(myRNA.getMobilizedBody(MobilizedBodyIndex(20)));
    sub.realizeTopology();
    cout << "sub.ancestor=" << sub.getAncestorMobilizedBodyIndex();
//    cout << "  sub.terminalBodies=" << sub.getTerminalBodies() << endl;
//    cout << "sub.allBodies=" << sub.getAllBodies() << endl;
    for (SubtreeBodyIndex i(0); i < (int)sub.getAllBodies().size(); ++i) {
       cout << "sub.parent[" << i << "]=" << sub.getParentSubtreeBodyIndex(i);
//       cout << "  sub.children[" << i << "]=" << sub.getChildSubtreeBodyIndexs(i) << endl;
    }


    printf("# quaternions in use = %d\n", myRNA.getNumQuaternionsInUse(s));
    for (MobilizedBodyIndex i(0); i<myRNA.getNumBodies(); ++i) {
        printf("body %2d: using quat? %s; quat index=%d\n",
            (int)i, myRNA.isUsingQuaternion(s,i) ? "true":"false",
            (int)myRNA.getQuaternionPoolIndex(s,i));
    }

    // And a study using the Runge Kutta Merson integrator
    bool suppressProject = false;

    RungeKuttaMersonIntegrator myStudy(mbs);
    //RungeKuttaFeldbergIntegrator myStudy(mbs);
    //RungeKutta3Integrator myStudy(mbs);
    //CPodesIntegrator  myStudy(mbs);
    //VerletIntegrator myStudy(mbs);
    //ExplicitEulerIntegrator myStudy(mbs);

    myStudy.setAccuracy(1e-2);
    myStudy.setConstraintTolerance(1e-3);
    myStudy.setProjectEveryStep(false);

    Visualizer display(mbs);
    display.setBackgroundColor(White);
    display.setBackgroundType(Visualizer::SolidColor);
    display.setMode(Visualizer::RealTime);

    for (MobilizedBodyIndex i(1); i<myRNA.getNumBodies(); ++i)
        myRNA.decorateBody(i, display);
    myRNA.decorateGlobal(display);

    DecorativeLine rbProto; rbProto.setColor(Orange).setLineThickness(3);
    display.addRubberBandLine(GroundIndex, attachPt,MobilizedBodyIndex(myRNA.getNumBodies()-1),Vec3(0), rbProto);
    //display.addRubberBandLine(GroundIndex, -attachPt,myRNA.getNumBodies()-1,Vec3(0), rbProto);

    const Real dt = 1./30; // output intervals

    printf("time  nextStepSize\n");

    s.updTime() = 0;
    for (int i=0; i<50; ++i)
        saveEm.push_back(s);    // delay
    display.report(s);

    myStudy.initialize(s);
    cout << "Using Integrator " << std::string(myStudy.getMethodName()) << ":\n";
    cout << "ACCURACY IN USE=" << myStudy.getAccuracyInUse() << endl;
    cout << "CTOL IN USE=" << myStudy.getConstraintToleranceInUse() << endl;
    cout << "TIMESCALE=" << mbs.getDefaultTimeScale() << endl;
    cout << "U WEIGHTS=" << s.getUWeights() << endl;
    cout << "Z WEIGHTS=" << s.getZWeights() << endl;
    cout << "1/QTOLS=" << s.getQErrWeights() << endl;
    cout << "1/UTOLS=" << s.getUErrWeights() << endl;

    saveEm.push_back(myStudy.getState());
    for (int i=0; i<50; ++i)
        saveEm.push_back(myStudy.getState());    // delay
    display.report(myStudy.getState());


    const double startReal = realTime(), startCPU = cpuTime();
    int stepNum = 0;
    for (;;) {
        const State& ss = myStudy.getState();

        mbs.realize(ss);

        if ((stepNum++%100)==0) {
            printf("%5g qerr=%10.4g uerr=%10.4g hNext=%g\n", ss.getTime(),
                myRNA.getQErr(ss).normRMS(), myRNA.getUErr(ss).normRMS(),
                myStudy.getPredictedNextStepSize());
            printf("      E=%14.8g (pe=%10.4g ke=%10.4g)\n",
                mbs.calcEnergy(ss), mbs.calcPotentialEnergy(ss), mbs.calcKineticEnergy(ss));

            cout << "QERR=" << ss.getQErr() << endl;
            cout << "UERR=" << ss.getUErr() << endl;
        }

        //if (s.getTime() - std::floor(s.getTime()) < 0.2)
        //    display.addEphemeralDecoration( DecorativeSphere(10).setColor(Green) );

        display.report(ss);
        saveEm.push_back(ss);

        if (ss.getTime() >= 10)
            break;

        // TODO: should check for errors or have or teach RKM to throw.
        myStudy.stepTo(ss.getTime() + dt, Infinity);
    }

    printf("CPU time=%gs, REAL time=%gs\n", cpuTime()-startCPU, realTime()-startReal);
    printf("Using Integrator %s:\n", myStudy.getMethodName());
    printf("# STEPS/ATTEMPTS = %d/%d\n", myStudy.getNumStepsTaken(), myStudy.getNumStepsAttempted());
    printf("# ERR TEST FAILS = %d\n", myStudy.getNumErrorTestFailures());
    printf("# CONVERGENCE FAILS = %d\n", myStudy.getNumConvergenceTestFailures());
    printf("# REALIZE/PROJECT = %d/%d\n", myStudy.getNumRealizations(), myStudy.getNumProjections());
    printf("# PROJECTION FAILS = %d\n", myStudy.getNumProjectionFailures());

    display.dumpStats(std::cout);

    while(true) {
        for (int i=0; i < (int)saveEm.size(); ++i) {
            display.report(saveEm[i]);
            //display.report(saveEm[i]); // half speed
        }
        getchar();
    }
  }
catch (const exception& e)
  {
    printf("EXCEPTION THROWN: %s\n", e.what());
    exit(1);
  }
}

