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

#include "SimTKmath.h"

#include "PendulumSystem.h"

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

using namespace SimTK;

using std::printf;
using std::cout;
using std::endl;

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
    void handleEvent(State& state, Real accuracy, const Vector& yWeights, const Vector& ooConstraintTols, Stage& lowestModified, bool& shouldTerminate) const {
        
        // This should be triggered when the pendulum reaches its farthest point in the
        // negative direction: q == -1, u == 0.
        
        Real q = state.getQ(pendulum.getGuts().getSubsysIndex())[0];
        Real u = state.getU(pendulum.getGuts().getSubsysIndex())[0];
        ASSERT(std::abs(q+1.0) < 0.05);
        ASSERT(std::abs(u) < 0.01);
        ASSERT(state.getTime() > lastEventTime);
        eventCount++;
        lastEventTime = state.getTime();
    }
private:
    PendulumSystem& pendulum;
};

class PeriodicHandler : public ScheduledEventHandler {
public:
    static int eventCount;
    static Real lastEventTime;
    Real getNextEventTime(const State&, bool includeCurrentTime) const {
        return lastEventTime+1.5;
    }
    void handleEvent(State& state, Real accuracy, const Vector& yWeights, const Vector& ooConstraintTols, Stage& lowestModified, bool& shouldTerminate) const {
        
        // This should be triggered every 1.5 time units.
        
        ASSERT(state.getTime() == lastEventTime+1.5);
        eventCount++;
        lastEventTime = state.getTime();
    }
};

class ZeroPositionReporter : public TriggeredEventReporter {
public:
    static int eventCount;
    static Real lastEventTime;
    ZeroPositionReporter(PendulumSystem& pendulum) : TriggeredEventReporter(Stage::Velocity), pendulum(pendulum) {
    }
    Real getValue(const State& state) const {
        return state.getQ(pendulum.getGuts().getSubsysIndex())[0];
    }
    void handleEvent(const State& state) const {
        
        // This should be triggered when the pendulum crosses q == 0.
        
        Real q = state.getQ(pendulum.getGuts().getSubsysIndex())[0];
        ASSERT(std::abs(q) < 0.01);
        ASSERT(state.getTime() > lastEventTime);
        eventCount++;
        lastEventTime = state.getTime();
    }
private:
    PendulumSystem& pendulum;
};

class PeriodicReporter : public ScheduledEventReporter {
public:
    static int eventCount;
    static Real lastEventTime;
    Real getNextEventTime(const State&, bool includeCurrentTime) const {
        return lastEventTime*2;
    }
    void handleEvent(const State& state) const {
        
        // This should be triggered every 1.5 time units.
        
        ASSERT(state.getTime() == lastEventTime*2);
        eventCount++;
        lastEventTime = state.getTime();
    }
};

int ZeroVelocityHandler::eventCount = 0;
Real ZeroVelocityHandler::lastEventTime = 0.0;
int PeriodicHandler::eventCount = 0;
Real PeriodicHandler::lastEventTime = 0.0;
int ZeroPositionReporter::eventCount = 0;
Real ZeroPositionReporter::lastEventTime = 0.0;
int PeriodicReporter::eventCount = 0;
Real PeriodicReporter::lastEventTime = 0.5;

int main () {
  try {
    PendulumSystem sys;
    sys.addEventHandler(new ZeroVelocityHandler(sys));
    sys.addEventHandler(new PeriodicHandler());
    sys.addEventReporter(new ZeroPositionReporter(sys));
    sys.addEventReporter(new PeriodicReporter());
    sys.realizeTopology();

    RungeKuttaMersonIntegrator integ(sys);

    const Real t0=0;
    const Real qi[] = {1,0}; // (x,y)=(1,0)
    const Real ui[] = {0,0}; // v=0
    const Vector q0(2, qi);
    const Vector u0(2, ui);

    sys.setDefaultMass(10);
    sys.setDefaultTimeAndState(t0, q0, u0);

    integ.setAccuracy(1e-2);
    integ.setConstraintTolerance(1e-4);

    const Real tFinal = 20.003;
    const Real hReport = 1.;

    integ.setFinalTime(tFinal);
    
    TimeStepper ts(sys);
    ts.setIntegrator(integ);
    ts.initialize(sys.getDefaultState());
    ASSERT(ts.getTime() == 0.0);
    ts.stepTo(10.0);
    ASSERT(ts.getTime() == 10.0);
    ts.stepTo(50.0);
    ASSERT(ts.getTime() == tFinal);
    ASSERT(integ.getTerminationReason() == Integrator::ReachedFinalTime);
    ASSERT(ZeroVelocityHandler::eventCount >= 10);
    ASSERT(PeriodicHandler::eventCount == (int) (ts.getTime()/1.5));
    ASSERT(ZeroPositionReporter::eventCount > 10);
    ASSERT(PeriodicReporter::eventCount == (int) (std::log(ts.getTime())/std::log(2.0))+1);
    cout << "Done" << endl;
    return 0;
  }
  catch (std::exception& e) {
    std::printf("FAILED: %s\n", e.what());
    return 1;
  }
}
