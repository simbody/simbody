#ifndef SimTK_SIMMATH_INTEGRATOR_TEST_FRAMEWORK_H_
#define SimTK_SIMMATH_INTEGRATOR_TEST_FRAMEWORK_H_

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

/**
 * This file contains code which is used for testing various integrators.
 */

#include "SimTKcommon.h"

#include "SimTKmath.h"
#include "simmath/TimeStepper.h"

#include "PendulumSystem.h"

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

using namespace SimTK;

using std::printf;
using std::cout;
using std::endl;

class PeriodicHandler : public PeriodicEventHandler {
public:
    static int eventCount;
    static PeriodicHandler* handler;
    PeriodicHandler() : PeriodicEventHandler(1.0) {
    }
    void handleEvent(State& state, Real accuracy, bool& shouldTerminate) const {
        
        // This should be triggered every (interval) time units.
        
        ASSERT(state.getTime() == getNextEventTime(state, true));
        eventCount++;
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
    void handleEvent(State& state, Real accuracy, bool& shouldTerminate) const {
        
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
            SubsystemIndex subsys = pendulum.getGuts().getSubsysIndex();
            state.updU(subsys) *= std::sqrt(1.5);
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
    void handleEvent(State& state, Real accuracy, bool& shouldTerminate) const {
        
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

class PeriodicReporter : public PeriodicEventReporter {
public:
    static int eventCount;
    static PeriodicReporter* reporter;
    PeriodicReporter(PendulumSystem& pendulum) : PeriodicEventReporter(1.0), pendulum(pendulum) {
    }
    void handleEvent(const State& state) const {
        
        // This should be triggered every (interval) time units.
        
        ASSERT(state.getTime() == getNextEventTime(state, true));
        eventCount++;
        
        // Verify conservation of energy.
        
        const Vector q = state.getQ(pendulum.getGuts().getSubsysIndex());
        const Vector u = state.getU(pendulum.getGuts().getSubsysIndex());
        Real energy =   pendulum.getMass(state)*(0.5*(u[0]*u[0]+u[1]*u[1])
                      + pendulum.getGravity(state)*(1.0+q[1]));
        Real expectedEnergy = pendulum.getMass(state)
                              * pendulum.getGravity(state);
        if (ZeroPositionHandler::hasAccelerated)
            expectedEnergy *= 1.5;
        ASSERT(std::abs(1.0-energy/expectedEnergy) < 0.05);
    }
private:
    PendulumSystem& pendulum;
};

class OnceOnlyEventReporter : public ScheduledEventReporter {
public:
    static bool hasOccurred;
    OnceOnlyEventReporter() {
    }
    Real getNextEventTime(const State&, bool includeCurrentTime) const {
        return 5.0;
    }
    void handleEvent(const State& state) const {
        ASSERT(!hasOccurred);
        hasOccurred = true;
    }
};

class DiscontinuousReporter : public TriggeredEventReporter {
public:
    static int eventCount;
    DiscontinuousReporter() : TriggeredEventReporter(Stage::Time) {
    }
    Real getValue(const State& state) const {
        Real step = std::floor(state.getTime());
        step = std::fmod(step, 4.0);
        if (step == 0.0)
            return 1.0;
        if (step == 2.0)
            return -1.0;
        return 0.0;
    }
    void handleEvent(const State& state) const {
        
        // This should be triggered when the value goes to 0, but not when it leaves 0.
        
        Real t = state.getTime();
        Real phase = std::fmod(t, 2.0);
        ASSERT(std::abs(phase-1.0) < 0.01);
        eventCount++;
    }
};


int ZeroVelocityHandler::eventCount = 0;
Real ZeroVelocityHandler::lastEventTime = 0.0;
int PeriodicHandler::eventCount = 0;
PeriodicHandler* PeriodicHandler::handler = 0;
int ZeroPositionHandler::eventCount = 0;
Real ZeroPositionHandler::lastEventTime = 0.0;
bool ZeroPositionHandler::hasAccelerated = false;
int PeriodicReporter::eventCount = 0;
PeriodicReporter* PeriodicReporter::reporter = 0;
bool OnceOnlyEventReporter::hasOccurred = false;
int DiscontinuousReporter::eventCount = 0;

void testIntegrator (Integrator& integ, PendulumSystem& sys, Real accuracy=1e-4) {
    ZeroVelocityHandler::eventCount = 0;
    ZeroVelocityHandler::lastEventTime = 0.0;
    PeriodicHandler::eventCount = 0;
    ZeroPositionHandler::eventCount = 0;
    ZeroPositionHandler::lastEventTime = 0.0;
    ZeroPositionHandler::hasAccelerated = false;
    PeriodicReporter::eventCount = 0;
    OnceOnlyEventReporter::hasOccurred = false;
    DiscontinuousReporter::eventCount = 0;

    const Real t0=0;
    const Real tFinal = 20.003;
    const Real qi[] = {1,0}; // (x,y)=(1,0)
    const Real ui[] = {0,0}; // v=0
    const Vector q0(2, qi);
    const Vector u0(2, ui);

    sys.setDefaultMass(10);
    sys.setDefaultTimeAndState(t0, q0, u0);
    integ.setAccuracy(accuracy);
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
    ASSERT(!OnceOnlyEventReporter::hasOccurred);
    ASSERT(!ZeroPositionHandler::hasAccelerated);
    
    // Try some steps of random sizes.
    
    static Random::Uniform random(0.0, 1.0);
    for (; time < 10.0; time += random.getValue()) {
        ASSERT(OnceOnlyEventReporter::hasOccurred == (ts.getTime() >= 5.0));
        ts.stepTo(time);
        ASSERT(ts.getTime() == time);
    }
    ASSERT(ZeroPositionHandler::hasAccelerated);
    
    // Try one large step that goes beyond tFinal.
    
    ts.stepTo(50.0);
    ASSERT(ts.getTime() == tFinal);
    ASSERT(integ.getTerminationReason() == Integrator::ReachedFinalTime);
    ASSERT(ZeroVelocityHandler::eventCount > 10);
    ASSERT(PeriodicHandler::eventCount == (int) (ts.getTime()/PeriodicHandler::handler->getEventInterval())+1);
    ASSERT(ZeroPositionHandler::eventCount > 10);
    ASSERT(PeriodicReporter::eventCount == (int) (ts.getTime()/PeriodicReporter::reporter->getEventInterval())+1);
    ASSERT(DiscontinuousReporter::eventCount == (int) (ts.getTime()/2.0));
}

#endif /*SimTK_SIMMATH_INTEGRATOR_TEST_FRAMEWORK_H_*/
