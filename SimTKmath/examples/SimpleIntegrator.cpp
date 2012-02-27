/* -------------------------------------------------------------------------- *
 *                                Simbody(tm)                                 *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2012 Stanford University and the Authors.           *
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


/* This example shows how to make use of the Simbody numerical integrators
to integrate an ordinary ODE. These integrators are designed to integrate the
complicated DAEs that arise in multibody dynamics so there is more overhead
in setting up the problem than you would expect for just an ODE. However, you
can use the framework supplied here for any simple ODE. 

The approach here is to provide stubs for unwanted features such as constraint
projection and event handling so that we can ignore them. */

// We need only the SimTKcommon and SimTKmath declarations since we're not
// doing any multibody dynamics here.
#include "SimTKmath.h"
#include "SimTKcommon/internal/SystemGuts.h"

#include <cstdio>
#include <iostream>

using namespace SimTK;

// Fill in your function here. This one is
//      sin( w t )
const Real w = (2 * Pi) / 4; // period
Real myFunction(Real t) {
    return std::sin(w*t);
}

// For testing, put the analytic result here. This one is
//     (1-cos(w t))/w
Real analyticIntegral(Real t) {
    return (1-std::cos(w*t)) / w;
}

// ---------- NOTHING TO SEE HERE, BUDDY. KEEP WALKING. ----------
// Simbody integrators act on System objects to advance their States, so we'll
// have to make a minimal System and give it a State that corresponds to the
// function of interest.
//
// A System is actually two classes: System::Guts does the work while System
// provides a pleasant veneer to make usage easier.
class MySystem;
class MySystemGuts : public System::Guts {
    friend class MySystem;

    // Implement required System::Guts virtuals.
    MySystemGuts* cloneImpl() const {return new MySystemGuts(*this);}

    // During realizeTopology() we allocate the needed State.
    int realizeTopologyImpl(State& state) const {
        // HERE'S WHERE THE IC GETS SET
        Vector zInit(1, 0.); // initial value for z
        state.allocateZ(SubsystemIndex(0), zInit);
        return 0;
    }

    // During realizeAcceleration() we calculate the State derivative.
    int realizeAccelerationImpl(const State& state) const {
        const Real t = state.getTime();
        // HERE'S THE CALL TO YOUR FUNCTION
        state.updZDot()[0] = myFunction(t);
        return 0;
    }

    // Disable prescribe and project since we have no constraints or
    // prescribed state variables to worry about.
    int prescribeImpl(State&, Stage) const {return 0;}
    int projectImpl(State&, Real, const Vector&, const Vector&, 
                    Vector&, System::ProjectOptions) const {return 0;}
};

class MySystem : public System {
public:
    MySystem() {
        adoptSystemGuts(new MySystemGuts());
        DefaultSystemSubsystem defsub(*this);
    }
};
// ---------- OKAY ----------

const Real startTime = 0;
const Real finalTime = 1;
const Real reportInterval = 0.1;
int main () {
    try {
        MySystem sys;
        State initState = sys.realizeTopology();
        initState.setTime(startTime);

        RungeKuttaMersonIntegrator integ(sys);
        //RungeKutta3Integrator integ(sys);
        integ.setAccuracy(1e-6);
        integ.setFinalTime(finalTime);
        integ.setReturnEveryInternalStep(true);
        
        Real timer = realTime();
        const int NReps = 1;
        for (int i=0; i < NReps; ++i) {
            initState.setTime(startTime);
            integ.initialize(initState);
            while (true) {
                const Real t = integ.getTime();
                // Use this for even report intervals:
                //Integrator::SuccessfulStepStatus 
                //    status = integ.stepTo(t + reportInterval);

                // Use this for variable step output.
                Integrator::SuccessfulStepStatus 
                    status = integ.stepTo(Infinity);

                if (status == Integrator::EndOfSimulation)
                    break;

                const State& state = integ.getState();
                const Real numerical = state.getZ()[0];
                const Real analytic = analyticIntegral(state.getTime());

                printf("%10g %10g\n", t, numerical-analytic);
            }
        }
        printf("Runtime=%g\n", (realTime()-timer)/NReps);

    } 
    catch (std::exception& e) {
        std::printf("FAILED: %s\n", e.what());
        return 1;
    }

    return 0;
}
