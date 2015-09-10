/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2009-12 Stanford University and the Authors.        *
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
 * This is an outer block for simulating ??? in various ways with Simbody.
 * This is about testing Simbody, *not* studying ???!
 */

#include "Simbody.h"

#include <string>
#include <iostream>
#include <exception>

using std::cout;
using std::endl;

using namespace SimTK;

const Real ReportInterval=0.033;
const Real RunTime=20;

class ShowLocking : public DecorationGenerator {
public:
    ShowLocking(const MobilizedBody& mobod, 
        const Constraint::ConstantSpeed& lock)
        : m_mobod(mobod), m_lock(lock)
    {
    }
    void generateDecorations(const State& state, 
                             Array_<DecorativeGeometry>& geometry) override
    {
        if (!m_lock.isDisabled(state)) {
            const Transform& X_BM = m_mobod.getDefaultOutboardFrame();
            geometry.push_back(DecorativeSphere(.5)
                .setTransform(X_BM)
                .setColor(Red).setOpacity(.25)
                .setBodyId(m_mobod.getMobilizedBodyIndex()));
            geometry.push_back(DecorativeText("LOCKED")
                .setTransform(X_BM+Vec3(.5,0,0))
                .setBodyId(m_mobod.getMobilizedBodyIndex()));
        }
    }
private:
    const MobilizedBody&                    m_mobod;
    const Constraint::ConstantSpeed&        m_lock;
};

class StateSaver : public PeriodicEventReporter {
public:
    StateSaver(const MultibodySystem&           system,
               const Constraint::ConstantSpeed& lock,
               const Integrator&                integ,
               Real reportInterval)
    :   PeriodicEventReporter(reportInterval), 
        m_system(system), m_lock(lock), m_integ(integ) 
    {   m_states.reserve(2000); }

    ~StateSaver() {}

    void clear() {m_states.clear();}
    int getNumSavedStates() const {return m_states.size();}
    const State& getState(int n) const {return m_states[n];}

    void handleEvent(const State& s) const override {
        const SimbodyMatterSubsystem& matter=m_system.getMatterSubsystem();
        const SpatialVec PG = matter.calcSystemMomentumAboutGroundOrigin(s);

        const bool isLocked = !m_lock.isDisabled(s);

        printf("%3d: %5g mom=%g,%g E=%g %s", m_integ.getNumStepsTaken(),
            s.getTime(),
            PG[0].norm(), PG[1].norm(), m_system.calcEnergy(s),
            isLocked?"LOCKED":"FREE");

        if (isLocked) {
            m_system.realize(s, Stage::Acceleration);
            printf(" lambda=%g", m_lock.getMultiplier(s));
        }
        cout << " Triggers=" << s.getEventTriggers() << endl;

        m_states.push_back(s);
    }
private:
    const MultibodySystem&              m_system;
    const Constraint::ConstantSpeed&    m_lock;
    const Integrator&                   m_integ;
    mutable Array_<State,int>           m_states;
};

class LockOn: public TriggeredEventHandler {
public:
    LockOn(const MultibodySystem& system,
        const MobilizedBody& mobod, Real lockangle, // must be 1dof
        const Constraint::ConstantSpeed& lock,
        Real low, Real high) 
    :   TriggeredEventHandler(Stage::Position), 
        m_mbs(system), m_mobod(mobod), m_lockangle(lockangle),
        m_lock(lock), m_low(low), m_high(high)
    { 
        //getTriggerInfo().setTriggerOnRisingSignTransition(false);
    }

    const Array_<Real>& getOnTimes() const {return m_onTimes;}

    Real getValue(const State& state) const override {
        if (!m_lock.isDisabled(state)) 
            return 0; // already locked
        const Real qdist = m_mobod.getOneQ(state, 0) - m_lockangle;
        return qdist;
    }

    void handleEvent
       (State& s, Real accuracy, bool& shouldTerminate) const override 
    {
        const SimbodyMatterSubsystem& matter = m_mbs.getMatterSubsystem();
        assert(m_lock.isDisabled(s));

        const Vector uin = s.getU();
        cout << "BEFORE u=" << uin << endl;
        cout << "before uerr=" << s.getUErr() << endl;

        // This is a ground-connected system -- system "center of mass" doesn't
        // really mean anything since ground's mass is infinite.
        SpatialVec PG = matter.calcSystemMomentumAboutGroundOrigin(s);
        SpatialVec PC = matter.calcSystemCentralMomentum(s);

        printf("Locking: BEFORE q=%.15g u=%.15g\n",
            m_mobod.getOneQ(s,0), m_mobod.getOneU(s,0));
        printf("  %5g G mom=%g,%g C mom=%g,%g E=%g\n", s.getTime(),
            PG[0].norm(), PG[1].norm(), PC[0].norm(), PC[1].norm(), m_mbs.calcEnergy(s));

        const Real CoefRest = 0.0, MinVel=.05;

        const UIndex ux = m_mobod.getFirstUIndex(s);
        const Real saveQ = m_mobod.getOneQ(s,0);

        m_lock.enable(s);
        m_mbs.realize(s, Stage::Dynamics);

        // We're using Poisson's definition of the coefficient of 
        // restitution, relating impulses, rather than Newton's, 
        // relating velocities, since Newton's can produce non-physical 
        // results for a multibody system. For Poisson, calculate the impulse
        // that would bring the velocity to zero, multiply by the coefficient
        // of restitution to calculate the rest of the impulse, then apply
        // both impulses to produce changes in velocity. In most cases this
        // will produce the same rebound velocity as Newton, but not always.

        Vector desiredDeltaV; // in constraint space
        Vector myDeltaV(1); // just one scalar for this constraint
        myDeltaV[0] = -uin[ux]; // dump the current velocity first
        m_lock.setMyPartInConstraintSpaceVector(s, myDeltaV, desiredDeltaV);

        Vector lambda, f, deltaU;


        cout << "USING POISSON COEF " << CoefRest << endl;

        Vector lambda0; // impulse that gets this velocity to zero
        matter.solveForConstraintImpulses(s, desiredDeltaV, lambda0);

        // This is the total impulse we want, in constraint space.
        lambda = (1+CoefRest)*lambda0;
        // Convert constraint impulse to generalized impulse.
        matter.multiplyByGTranspose(s,lambda,f);
        // Convert generalized impulse to generalized speed change.
        matter.multiplyByMInv(s,f,deltaU);
        // If the new speed is slow enough, we'll declare that it "stuck" and
        // drive it to exactly zero instead (using lambda0 instead of lambda).
        Real achievedU = uin[ux]+deltaU[ux];
        if (CoefRest > 0 && std::abs(achievedU) < MinVel) {
            cout << "  rebound " << achievedU 
                 << " too slow; drive to zero instead\n";
            cout << "  deltaU was=" << deltaU << endl;
            lambda = lambda0; f /= 1+CoefRest; deltaU /= 1+CoefRest;
            cout << "  new deltaU=" << deltaU << endl;
        }

        // Now update all the generalized speeds.
        s.updU() = uin + deltaU;
        cout << "lambda=" << lambda << endl;
        cout << "f=" << f << endl;
        cout << "du=" << deltaU << endl;       

        cout << "AFTER u=" << s.getU() << endl;

        // If this leaves us rebounding then don't turn on the constraint.
        if (std::abs(s.getU()[ux]) > SignificantReal) {
            printf("REBOUND vel=%g -- not locking.\n", s.getU()[1]);
            m_lock.disable(s);
            // Move q a small distance in the direction of u so that the
            // witness function will be non-zero and we'll consider an 
            // immediate reverse as a new event (transitions away from zero
            // are not events).
            m_mobod.setOneQ(s,0,m_lockangle+sign(s.getU()[ux])*SignificantReal); 

            m_mbs.realize(s, Stage::Acceleration);
            PG = matter.calcSystemMomentumAboutGroundOrigin(s);
            PC = matter.calcSystemCentralMomentum(s);
            printf("  %5g G mom=%g,%g C mom=%g,%g E=%g\n", s.getTime(),
                PG[0].norm(), PG[1].norm(), PC[0].norm(), PC[1].norm(), 
                m_mbs.calcEnergy(s));
            cout << "  uerr=" << s.getUErr() << endl;
            return;
        }

        // Tidy up the q to be exactly zero when we lock to avoid
        // accidental retriggering fo this event. We'll put this back
        // if we decide below not to hold the lock.
        m_mobod.setOneQ(s,0,m_lockangle);

        // Can we really hold the lock?
        m_mbs.realize(s, Stage::Acceleration);
        const Real lockForce = m_lock.getMultiplier(s);

        if (lockForce < m_low || lockForce > m_high) {
            m_lock.disable(s); // oops can't lock
            s.updU() = uin;
            m_mobod.setOneQ(s,0, saveQ); 
            printf("CAN'T LOCK: force would have been %g\n", lockForce);
            return;
        }

        printf("LOCKED: reaction force is now %g\n", lockForce);

        m_onTimes.push_back(s.getTime());

        printf("  after q=%.15g\n", m_mobod.getOneQ(s,0));

        PG = matter.calcSystemMomentumAboutGroundOrigin(s);
        PC = matter.calcSystemCentralMomentum(s);
        printf("  %5g G mom=%g,%g C mom=%g,%g E=%g\n", s.getTime(),
            PG[0].norm(), PG[1].norm(), PC[0].norm(), PC[1].norm(), m_mbs.calcEnergy(s));
        cout << "  after uerr=" << s.getUErr() << endl;
    }

private:
    const MultibodySystem&                  m_mbs; 
    const MobilizedBody&                    m_mobod;
    const Real                              m_lockangle;
    const Constraint::ConstantSpeed&        m_lock;
    const Real                              m_low, m_high;

    mutable Array_<Real> m_onTimes;
};

class LockOff: public TriggeredEventHandler {
public:
    LockOff(const MultibodySystem& system,
        MobilizedBody& mobod, Real lockangle, 
        Constraint::ConstantSpeed& lock,
        Real low, Real high) 
    :   TriggeredEventHandler(Stage::Acceleration), 
        m_system(system), m_mobod(mobod), m_lockangle(lockangle), m_lock(lock), 
        m_low(low), m_high(high)
    { 
        getTriggerInfo().setTriggerOnRisingSignTransition(false);
    }

    const Array_<Real>& getOffTimes() const {return m_offTimes;}

    Real getValue(const State& state) const override {
        if (m_lock.isDisabled(state)) return 0;
        const Real f = m_lock.getMultiplier(state);
        const Real mid = (m_high+m_low)/2;
        return f > mid ? m_high - f : f - m_low;
    }

    void handleEvent
       (State& s, Real accuracy, bool& shouldTerminate) const override 
    {
        assert(!m_lock.isDisabled(s));

        m_system.realize(s, Stage::Acceleration);
        printf("\nUNLOCK: disabling at t=%g q=%g lambda=%g",
            s.getTime(), s.getQ()[1], m_lock.getMultiplier(s));
        cout << " Triggers=" << s.getEventTriggers() << "\n\n";

        m_mobod.setOneQ(s, 0, m_lockangle); // avoid retriggering lock

        m_lock.disable(s);

        m_offTimes.push_back(s.getTime());
    }

private:
    const MultibodySystem&                  m_system; 
    const MobilizedBody&                    m_mobod;
    const Real                              m_lockangle;
    const Constraint::ConstantSpeed&        m_lock;
    const Real                              m_low;
    const Real                              m_high;

    mutable Array_<Real> m_offTimes;
};

static const Real Deg2Rad = (Real)SimTK_DEGREE_TO_RADIAN,
                  Rad2Deg = (Real)SimTK_RADIAN_TO_DEGREE;



static Real g = 9.8;

int main(int argc, char** argv) {
    static const Transform GroundFrame;
    static const Rotation ZUp(UnitVec3(XAxis), XAxis, UnitVec3(YAxis), ZAxis);
    static const Vec3 TestLoc(1,0,0);

  try { // If anything goes wrong, an exception will be thrown.

        // CREATE MULTIBODY SYSTEM AND ITS SUBSYSTEMS
    MultibodySystem             mbs;

    SimbodyMatterSubsystem      matter(mbs);
    GeneralForceSubsystem       forces(mbs);
    Force::Gravity              gravity(forces, matter, Vec3(0, -g, 0));

        // ADD BODIES AND THEIR MOBILIZERS
    const Vec3 thighHDim(.5,2,.25); 
    const Real thighVol=8*thighHDim[0]*thighHDim[1]*thighHDim[2];
    const Vec3 calfHDim(.25,2,.125); 
    const Real calfVol=8*calfHDim[0]*calfHDim[1]*calfHDim[2];
    const Real density = 1000; // water
    const Real thighMass = density*thighVol, calfMass = density*calfVol;
    Body::Rigid thighBody = 
        Body::Rigid(MassProperties(10*thighMass, Vec3(0), 
                        10*thighMass*UnitInertia::brick(thighHDim)));
    thighBody.addDecoration(Transform(), DecorativeBrick(thighHDim)
                                             .setColor(Red).setOpacity(.3));
    Body::Rigid calfBody = 
        Body::Rigid(MassProperties(calfMass, Vec3(0), 
                        calfMass*UnitInertia::brick(calfHDim)));
    calfBody.addDecoration(Transform(), DecorativeBrick(calfHDim)
                                            .setColor(Blue).setOpacity(.3));
    Body::Rigid footBody = 
        Body::Rigid(MassProperties(10*calfMass, Vec3(0), 
                        10*calfMass*UnitInertia::brick(calfHDim)));
    footBody.addDecoration(Transform(), DecorativeBrick(calfHDim)
                                            .setColor(Black).setOpacity(.3));
    MobilizedBody::Pin thigh(matter.Ground(), Vec3(0),
                             thighBody, Vec3(0,thighHDim[1],0));
    MobilizedBody::Pin calf(thigh, Vec3(0,-thighHDim[1],0),
                             calfBody, Vec3(0,calfHDim[1],0));
    MobilizedBody::Pin foot(calf, Vec3(0,-calfHDim[1],0),
                             footBody, Vec3(0,calfHDim[1],0));
    //Constraint::PrescribedMotion pres(matter, 
    //    new Function::Constant(Pi/4,1), foot, MobilizerQIndex(0));
    Constraint::PrescribedMotion pres(matter, 
       new Function::Sinusoid(Pi/4,2*Pi,-Pi/4), foot, MobilizerQIndex(0));

    Constraint::ConstantSpeed lock(calf,0);
    lock.setDisabledByDefault(true);

    Visualizer viz(mbs);
    viz.addDecorationGenerator(new ShowLocking(calf,lock));
    mbs.addEventReporter(new Visualizer::Reporter(viz, ReportInterval));

    //ExplicitEulerIntegrator integ(mbs);
    //CPodesIntegrator integ(mbs,CPodes::BDF,CPodes::Newton);
    //RungeKuttaFeldbergIntegrator integ(mbs);
    //RungeKuttaMersonIntegrator integ(mbs);
    RungeKutta3Integrator integ(mbs);
    //VerletIntegrator integ(mbs);
    integ.setAccuracy(1e-3);
    //integ.setAllowInterpolation(false);

    StateSaver* stateSaver = new StateSaver(mbs,lock,integ,ReportInterval);
    mbs.addEventReporter(stateSaver);

    const Real low=-110000*4, high=110000*4;
    const Real lockAngle = 0;
    LockOn* lockOn = new LockOn(mbs,calf,lockAngle,lock,low,high);
    mbs.addEventHandler(lockOn);

    LockOff* lockOff = new LockOff(mbs,calf,lockAngle,lock,low,high);
    mbs.addEventHandler(lockOff);
  
    State s = mbs.realizeTopology(); // returns a reference to the the default state
    mbs.realizeModel(s); // define appropriate states for this System
    mbs.realize(s, Stage::Instance); // instantiate constraints if any

    thigh.setAngle(s, 20*Deg2Rad);
    calf.setAngle(s, 90*Deg2Rad);
    //calf.setRate(s, -10);

    mbs.realize(s, Stage::Velocity);
    viz.report(s);

    mbs.realize(s, Stage::Acceleration);


    cout << "q=" << s.getQ() << endl;
    cout << "u=" << s.getU() << endl;
    cout << "qerr=" << s.getQErr() << endl;
    cout << "uerr=" << s.getUErr() << endl;
    cout << "udoterr=" << s.getUDotErr() << endl;
    cout << "mults=" << s.getMultipliers() << endl;
    cout << "qdot=" << s.getQDot() << endl;
    cout << "udot=" << s.getUDot() << endl;
    cout << "qdotdot=" << s.getQDotDot() << endl;
    viz.report(s);

    cout << "Initial configuration shown. Next? ";
    getchar();

    Assembler(mbs).assemble(s);
    viz.report(s);
    cout << "Assembled configuration shown. Ready? ";
    getchar();
    
    // Simulate it.
    const double start = realTime();

    // TODO: misses some transitions if interpolating
    //integ.setAllowInterpolation(false);
    TimeStepper ts(mbs, integ);
    ts.initialize(s);
    ts.stepTo(RunTime);

    const double timeInSec = realTime()-start;
    const int evals = integ.getNumRealizations();
    cout << "Done -- took " << integ.getNumStepsTaken() << " steps in " <<
        timeInSec << "s for " << ts.getTime() << "s sim (avg step=" 
        << (1000*ts.getTime())/integ.getNumStepsTaken() << "ms) " 
        << (1000*ts.getTime())/evals << "ms/eval\n";
    cout << "On times: " << lockOn->getOnTimes() << endl;
    cout << "Off times: " << lockOff->getOffTimes() << endl;

    printf("Used Integrator %s at accuracy %g:\n", 
        integ.getMethodName(), integ.getAccuracyInUse());
    printf("# STEPS/ATTEMPTS = %d/%d\n", integ.getNumStepsTaken(), integ.getNumStepsAttempted());
    printf("# ERR TEST FAILS = %d\n", integ.getNumErrorTestFailures());
    printf("# REALIZE/PROJECT = %d/%d\n", integ.getNumRealizations(), integ.getNumProjections());

    while(true) {
        for (int i=0; i < stateSaver->getNumSavedStates(); ++i) {
            viz.report(stateSaver->getState(i));
        }
        getchar();
    }

  } 
  catch (const std::exception& e) {
    printf("EXCEPTION THROWN: %s\n", e.what());
    exit(1);
  }
  catch (...) {
    printf("UNKNOWN EXCEPTION THROWN\n");
    exit(1);
  }

}
