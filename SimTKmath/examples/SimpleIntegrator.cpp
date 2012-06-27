/* -------------------------------------------------------------------------- *
 *                  Simbody(tm) Example: Simple Integrator                    *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2012 Stanford University and the Authors.           *
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

/* This example shows how to make use of the Simbody numerical integrators
to integrate an ordinary ODE. These integrators are designed to integrate the
complicated DAEs that arise in multibody dynamics so there is more overhead
in setting up the problem than you would expect for just an ODE. However, you
can use the framework supplied here for any simple ODE, and you can easily
extend it to integrate DAEs (ODEs + constraints) if you need to do so. */

// We need only the SimTKmath (and implicitly, SimTKcommon) declarations since 
// we're not doing any multibody dynamics here.
#include "SimTKmath.h"

#include <cstdio>
#include <iostream>

using namespace SimTK;

// We'll use this class to communicate the desired period to the
// function we're integrating. This is an example to show how you can
// pass any constant stuff you want to your function.
struct SomeStuff {
    Real w;
};

// Fill in your function here. This one is
//      sin( w t )
Real myFunction(const SomeStuff& stuff, Real t) {
    const Real w = stuff.w; // get the period
    return std::sin(w*t);
}

// For testing, put the analytic result here. Normally you wouldn't know
// this -- it is the function we're going to calculate by numerical
// integration. However, in the case of the above myFunction() the integral is:
//     (1-cos(w t))/w
Real analyticIntegral(const SomeStuff& stuff, Real t) {
    const Real w = stuff.w; // get the period
    return (1-std::cos(w*t)) / w;
}

// This class has the form required by the Simbody integrators; that is, it
// is a SimTK::System. It calls myFunction() to calculate its state derivative. 
// It is defined after the main program below.
class MySystem : public System {
public:
    explicit MySystem(const SomeStuff&);
};


const Real startTime = 0;
const Real finalTime = 1;
const Real reportInterval = 0.1;
int main () {
    try {
        SomeStuff myStuff; // constant data to pass to myFunction()
        myStuff.w = Pi/2; 
        MySystem sys(myStuff);
        State initState = sys.realizeTopology();
        initState.setTime(startTime);

        RungeKuttaMersonIntegrator integ(sys);
        //RungeKutta3Integrator integ(sys);
        //RungeKuttaFeldbergIntegrator integ(sys);
        //ExplicitEulerIntegrator integ(sys);
        //VerletIntegrator integ(sys);
        //CPodesIntegrator integ(sys);   // implicit integrator

        integ.setAccuracy(1e-6);
        integ.setFinalTime(finalTime);

        // Use this to get output that is spaced according to the local
        // complexity of the function being integrated. By default you will
        // just get output at specified reporting times, with the integrator
        // taking as many steps as necessary to get there with the specified
        // accuracy.
        //integ.setReturnEveryInternalStep(true);
        
        initState.setTime(startTime);
        integ.initialize(initState);

        Integrator::SuccessfulStepStatus status;
        std::printf("%6s %12s\n", "time", "error");
        // Keep stepping until the integrator returns EndOfSimulation.
        while (true) {
            const Real prevT = integ.getTime();

            // To get evenly-spaced reporting intervals, disable 
            // returnEveryInternalStep above and use this line.
            status = integ.stepTo(prevT + reportInterval);

            // Use this for variable step output, along with setting 
            // returnEveryInternalStep above.
            // status = integ.stepTo(Infinity);

            //std::printf("step to %g, status=%s\n", integ.getTime(),
            //    Integrator::getSuccessfulStepStatusString(status).c_str());

            if (status == Integrator::EndOfSimulation)
                break;

            // Note that the integrator maintains its own internal State.
            const State& state = integ.getState();
            const Real t = state.getTime();
            const Real numerical = state.getZ()[0];
            const Real analytic = analyticIntegral(myStuff, state.getTime());

            std::printf("%6g %12g\n", t, numerical-analytic);
        }

        std::printf("\nDone. Used %s with %d function calls.\n",
            integ.getMethodName(), integ.getNumRealizations());
        std::printf("  %d steps taken out of %d attempted.\n", 
            integ.getNumStepsTaken(), integ.getNumStepsAttempted());
    } 
    catch (const std::exception& e) {
        std::printf("FAILED: %s\n", e.what());
        return 1;
    }

    return 0;
}


//--------------- SYSTEM CLASS NEEDED FOR INTEGRATOR BOOKKEEPING ---------------
// Simbody integrators act on System objects to advance their States, so we'll
// have to make a minimal System and give it a State that corresponds to the
// function of interest.
//
// A System is actually two classes: System::Guts does the work while System
// provides a pleasant veneer to make usage easier.
class MySystemGuts : public System::Guts {
    friend class MySystem;

    MySystemGuts(const SomeStuff& stuff) : stuff(stuff) {}

    // Implement required System::Guts virtuals.
    MySystemGuts* cloneImpl() const {return new MySystemGuts(*this);}

    // During realizeTopology() we allocate the needed State.
    int realizeTopologyImpl(State& state) const {
        // HERE'S WHERE THE IC GETS SET
        const Vector zInit(1, 0.); // initial value for the one z we want
        state.allocateZ(SubsystemIndex(0), zInit);
        return 0;
    }

    // During realizeAcceleration() we calculate the State derivative.
    int realizeAccelerationImpl(const State& state) const {
        const Real t = state.getTime();
        // HERE'S THE CALL TO YOUR FUNCTION
        state.updZDot()[0] = myFunction(stuff, t);
        return 0;
    }

private:
    SomeStuff stuff;
};

// Define MySystem's constructor so that it creates a MySystemGuts object and
// a default Subsystem.
MySystem::MySystem(const SomeStuff& stuff) {
    adoptSystemGuts(new MySystemGuts(stuff));
    DefaultSystemSubsystem defsub(*this);
}
//------------------------ END OF OBSCURE BOOKKEEPING --------------------------

