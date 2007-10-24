/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmath(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-7 Stanford University and the Authors.         *
 * Authors: Michael Sherman, Peter Eastman                                    *
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

#include "SimTKcommon.h"

#include "SimTKmath.h"
#include "simmath/CPodesIntegrator.h" // will be Simmath
#include "simmath/TimeStepper.h" // will be Simmath

#include "PendulumSystem.h"

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

using namespace SimTK;

using std::printf;
using std::cout;
using std::endl;

class PeriodicEventHandler : public ScheduledEventHandler {
public:
    static int eventCount;
    static Real lastEventTime;
    static Real interval;
    Real getNextEventTime(const State&) const {
        return lastEventTime+interval;
    }
    void handleEvent(State& state, Real accuracy, const Vector& yWeights, const Vector& ooConstraintTols, Stage& lowestModified, bool& shouldTerminate) {
        
        // This should be triggered every (interval) time units.
        
        ASSERT(state.getTime() == lastEventTime+interval);
        eventCount++;
        lastEventTime = state.getTime();
    }
};

class ZeroPositionHandler : public TriggeredEventHandler {
public:
    static int eventCount;
    static Real lastEventTime;
    static bool hasAccelerated;
    ZeroPositionHandler(PendulumSystem& pendulum) : TriggeredEventHandler(Stage::Velocity), pendulum(pendulum) {
    }
    Real getValue(const State& state) const {
        return state.getQ(pendulum.getGuts().getSubsysIndex())[0];
    }
    void handleEvent(State& state, Real accuracy, const Vector& yWeights, const Vector& ooConstraintTols, Stage& lowestModified, bool& shouldTerminate) {
        
        // This should be triggered when the pendulum crosses x == 0.
        
        Real x = state.getQ(pendulum.getGuts().getSubsysIndex())[0];
        ASSERT(std::abs(x) < 0.01);
        ASSERT(state.getTime() > lastEventTime);
        eventCount++;
        lastEventTime = state.getTime();
        if (state.getTime() > 7 && !hasAccelerated) {

            // Multiply the pendulum's velocity by sqrt(1.5), which should multiply its total energy by 1.5 (since
            // at x=0, all of its energy is kinetic).
            
            hasAccelerated = true;
            SubsystemId subsys = pendulum.getGuts().getSubsysIndex();
            state.updU(subsys) *= std::sqrt(1.5);
            lowestModified = Stage::Velocity;
        }
    }
private:
    PendulumSystem& pendulum;
};

class ZeroVelocityHandler : public TriggeredEventHandler {
public:
    static int eventCount;
    static Real lastEventTime;
    ZeroVelocityHandler(PendulumSystem& pendulum) : TriggeredEventHandler(Stage::Velocity), pendulum(pendulum) {
        getTriggerInfo().setTriggerOnFallingSignTransition(false);
    }
    Real getValue(const State& state) const {
        return state.getU(pendulum.getGuts().getSubsysIndex())[0];
    }
    void handleEvent(State& state, Real accuracy, const Vector& yWeights, const Vector& ooConstraintTols, Stage& lowestModified, bool& shouldTerminate) {
        
        // This should be triggered when the pendulum reaches its farthest point in the
        // negative direction: q[0] == -1, u[0] == 0.
        
        Vector u = state.getU(pendulum.getGuts().getSubsysIndex());
        ASSERT(std::abs(u[0]) < 0.01);
        ASSERT(state.getTime() > lastEventTime);
        eventCount++;
        lastEventTime = state.getTime();
    }
private:
    PendulumSystem& pendulum;
};

class PeriodicEventReporter : public ScheduledEventReporter {
public:
    static int eventCount;
    static Real lastEventTime;
    static Real interval;
    PeriodicEventReporter(PendulumSystem& pendulum) : pendulum(pendulum) {
    }
    Real getNextEventTime(const State&) const {
        return lastEventTime+interval;
    }
    void handleEvent(const State& state) {
        
        // This should be triggered every (interval) time units.
        
        ASSERT(state.getTime() == lastEventTime+interval);
        eventCount++;
        lastEventTime = state.getTime();
        
        // Verify conservation of energy.
        
        const Vector q = state.getQ(pendulum.getGuts().getSubsysIndex());
        const Vector u = state.getU(pendulum.getGuts().getSubsysIndex());
        Real energy = pendulum.getMass(state)*(0.5*(u[0]*u[0]+u[1]*u[1])+pendulum.getGravity(state)*(1.0+q[1]));
        Real expectedEnergy = pendulum.getMass(state)*pendulum.getGravity(state);
        if (ZeroPositionHandler::hasAccelerated)
            expectedEnergy *= 1.5;
        ASSERT(std::abs(1.0-energy/expectedEnergy) < 0.05);
    }
private:
    PendulumSystem& pendulum;
};

int ZeroVelocityHandler::eventCount = 0;
Real ZeroVelocityHandler::lastEventTime = 0.0;
int PeriodicEventHandler::eventCount = 0;
Real PeriodicEventHandler::lastEventTime = 0.0;
Real PeriodicEventHandler::interval = 0.0;
int ZeroPositionHandler::eventCount = 0;
Real ZeroPositionHandler::lastEventTime = 0.0;
bool ZeroPositionHandler::hasAccelerated = false;
int PeriodicEventReporter::eventCount = 0;
Real PeriodicEventReporter::lastEventTime = 0.0;
Real PeriodicEventReporter::interval = 0.0;

void testIntegrator (Integrator& integ, PendulumSystem& sys) {
    ZeroVelocityHandler::eventCount = 0;
    ZeroVelocityHandler::lastEventTime = 0.0;
    PeriodicEventHandler::eventCount = 0;
    PeriodicEventHandler::lastEventTime = 0.0;
    ZeroPositionHandler::eventCount = 0;
    ZeroPositionHandler::lastEventTime = 0.0;
    ZeroPositionHandler::hasAccelerated = false;
    PeriodicEventReporter::eventCount = 0;
    PeriodicEventReporter::lastEventTime = 0.0;

    const Real t0=0;
    const Real tFinal = 20.003;
    const Real qi[] = {1,0}; // (x,y)=(1,0)
    const Real ui[] = {0,0}; // v=0
    const Vector q0(2, qi);
    const Vector u0(2, ui);

    sys.setDefaultMass(10);
    sys.setDefaultTimeAndState(t0, q0, u0);
    integ.setAccuracy(1e-4);
    integ.setConstraintTolerance(1e-4);
    integ.setFinalTime(tFinal);
    
    TimeStepper ts(sys);
    ts.setIntegrator(integ);
    ts.initialize(sys.getDefaultState());
    
    // Try taking a series of steps of constant size.
    
    ASSERT(ts.getTime() == 0.0);
    Real time = 1.0;
    for (; time < 5.0; time += 1.0) {
        ts.stepTo(time);
        ASSERT(ts.getTime() == time);
    }
    ASSERT(!ZeroPositionHandler::hasAccelerated);
    
    // Try some steps of random sizes.
    
    static Random::Uniform random(0.0, 1.0);
    for (; time < 10.0; time += random.getValue()) {
        ts.stepTo(time);
        ASSERT(ts.getTime() == time);
    }
    ASSERT(ZeroPositionHandler::hasAccelerated);
    
    // Try one large step that goes beyond tFinal.
    
    ts.stepTo(50.0);
    ASSERT(ts.getTime() == tFinal);
    ASSERT(integ.getTerminationReason() == Integrator::ReachedFinalTime);
    ASSERT(ZeroVelocityHandler::eventCount > 10);
    ASSERT(PeriodicEventHandler::eventCount == (int) (ts.getTime()/PeriodicEventHandler::interval));
    ASSERT(ZeroPositionHandler::eventCount > 10);
    ASSERT(PeriodicEventReporter::eventCount == (int) (ts.getTime()/PeriodicEventReporter::interval));
}

int main () {
  try {
    PendulumSystem sys;
    sys.updDefaultSubsystem().addEventHandler(new ZeroVelocityHandler(sys));
    sys.updDefaultSubsystem().addEventHandler(new PeriodicEventHandler());
    sys.updDefaultSubsystem().addEventHandler(new ZeroPositionHandler(sys));
    sys.updDefaultSubsystem().addEventReporter(new PeriodicEventReporter(sys));
    sys.realizeTopology();

    // Test with various intervals for the event handler and event reporter, ones that are either
    // large or small compared to the expected internal step size of the integrator.

    for (int i = 0; i < 4; ++i) {
        PeriodicEventHandler::interval = (i == 0 || i == 1 ? 0.01 : 2.0);
        PeriodicEventReporter::interval = (i == 0 || i == 2 ? 0.015 : 1.5);
        
        // Test the BDF integrator in both normal and single step modes.
        
        CPodesIntegrator bdfInteg(sys, CPodes::BDF);
        testIntegrator(bdfInteg, sys);
        bdfInteg.setReturnEveryInternalStep(true);
        testIntegrator(bdfInteg, sys);
        
        // Test the Adams integrator in both normal and single step modes.
        
        CPodesIntegrator adamsInteg(sys, CPodes::Adams);
        testIntegrator(adamsInteg, sys);
        adamsInteg.setReturnEveryInternalStep(true);
        testIntegrator(adamsInteg, sys);
        
        // Calling setUseCPodesProjection() after the integrator has been initialized should
        // produce an exception.
        
        try {
            adamsInteg.setUseCPodesProjection();
            assert(false);
        }
        catch (...) {
        }
        
        // Try having CPODES do the projection instead of the System.
        
        CPodesIntegrator projInteg(sys, CPodes::BDF);
        projInteg.setUseCPodesProjection();
        testIntegrator(projInteg, sys);
    }
    cout << "Done" << endl;
    return 0;
  }
  catch (std::exception& e) {
    std::printf("FAILED: %s\n", e.what());
    return 1;
  }
}
