#ifndef SimTK_SIMMATH_TIMESTEPPER_H_
#define SimTK_SIMMATH_TIMESTEPPER_H_

/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-15 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
 * Contributors: Michael Sherman                                              *
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

#include "SimTKcommon.h"
#include "simmath/internal/common.h"
#include "simmath/Integrator.h"

namespace SimTK {
class Integrator;

/** A %TimeStepper advances the State of a System through time, using an 
Integrator for continuous parts of the trajectory and handling discrete events
and reporting. The purpose of a %TimeStepper is to simplify typical usage of
an Integrator; you can use an Integrator directly if you want more control.

Typical usage:
@code
    // Create a suitable Integrator object for your System.
    RungeKuttaMersonIntegrator integrator(system);
    TimeStepper stepper(integrator);
    stepper.initialize(initialState);
    stepper.stepTo(finalTime);
@endcode

The stepTo() method invokes the Integrator repeatedly to advance time. It 
detects events which may occur, calls event handlers and reporters as
appropriate, then continues advancing time until the target time is reached.
**/
class SimTK_SIMMATH_EXPORT TimeStepper {
public:
    class EventChangeActionFailed;

    /** Create a %TimeStepper to advance a System using an Integrator. Ownership
    of the Integrator object does *not* pass to the %TimeStepper; you must not
    destruct the Integrator while a TimeStepper is using it. **/
    explicit TimeStepper(Integrator& integrator);

    /** Create a %TimeStepper leaving the Integrator unspecified. You must call 
    setIntegrator() before calling initialize(). **/
    TimeStepper();

    /** Set the Integrator this %TimeStepper will use to advance the System.
    This is required if you didn't supply an Integrator at construction. 
    Ownership of the Integrator object does *not* pass to the TimeStepper; you 
    must not destruct the Integrator while a TimeStepper is using it.  **/
    void setIntegrator(Integrator& integrator);

    /** Supply the %TimeStepper with a starting state, which is *copied* to 
    initialize the Integrator. This must be called after the Integrator has
    been set, and before the first call to stepTo(). The specified state is 
    *copied* into the Integrator's internally maintained state; subsequent 
    changes to the State object passed in here will not affect the 
    simulation. **/
    void initialize(const State&);

    /** Use the Integrator to advance the System's State from the current time
    to the specified time. This method will repeatedly invoke the Integrator as 
    necessary, handling any events which occur. The %TimeStepper must be
    initialized by calling initialize() before this method may be called.
    
    When this method returns, the System will usually have been advanced all the
    way to the specified final time, regardless of how many integration steps
    that requires. There are situations where it may return sooner, however: 
    if setReportAllSignificantStates() was set to true, and a "significant"
    state occurred; if {@link Integrator#setFinalTime setFinalTime()} was 
    invoked on the Integrator, and the final time was reached; or if an event 
    handler requested that the simulation terminate immediately. **/
    Integrator::SuccessfulStepStatus stepTo(double time);

    /** Get the current State of the System being integrated. Usually this will 
    correspond to the time specified in the most recent call to stepTo(). This
    is the State maintained internally by the Integrator; this call just 
    forwards to Integrator::getState().
    @see Integrator::getState() **/
    const State& getState() const;

    /** Get the current time of the System being integrated. This is just 
    shorthand for getState().getTime(). **/
    double getTime() const {return getState().getTime();}

    /** Get a const reference to the Integrator being used. **/
    const Integrator& getIntegrator() const;

    /** Get a writable reference to the Integrator being used. **/
    Integrator& updIntegrator();

    /** Note that the Integrator is not destructed with the %TimeStepper; it is
    an independent object. **/
    ~TimeStepper();

    /** @name             Advanced/Debugging/Deprecated 
    You probably won't need to use these. **/
    /**@{**/
    /** (Advanced) Set whether the %TimeStepper's stepTo() method should return 
    after any "significant" state is returned by the Integrator. Normally 
    stepTo() will only return when the specified time has been reached or when 
    the simulation is terminated, with output obtained only through Event 
    reporters and handlers. However, if this flag is set `true`, stepTo() will 
    return whenever the Integrator reports that something interesting has 
    happened, such as when an event occurs or at the start of a new continuous 
    interval after an event has been handled. This will not return every time
    a step is taken, unless the Integrator has been told separately to return
    every step. 
    @see Integrator::setReturnEveryInternalStep() **/
    void setReportAllSignificantStates(bool b);

    /** (Advanced) Return the current value of the "report all significant 
    states" flag.
    @see setReportAllSignificantState() **/
    bool getReportAllSignificantStates() const;

    /** <b>(Deprecated)</b> There was no need to supply a System to a 
    %TimeStepper since there is always a System associated with the Integrator.
    Use the default constructor `TimeStepper()` instead of this one. For 
    temporary backwards compatibility we're still allowing this signature 
    but the given System is ignored. **/
    DEPRECATED_14("use TimeStepper()")
    explicit TimeStepper(const System&) : TimeStepper() {}
    /** <b>(Deprecated)</b> There was no need to supply a System to a 
    %TimeStepper since there is always a System associated with the Integrator.
    Use the constructor `TimeStepper(Integrator)` instead of this one.
    For temporary backwards compatibility we're still allowing this signature 
    but the given System is ignored. **/
    DEPRECATED_14("use TimeStepper(Integrator)") 
    TimeStepper(const System&, Integrator& integ) 
    :   TimeStepper(integ) {}
    /**@}**/

private:
friend class TimeStepperRep;

    class TimeStepperRep* rep;
};

//==============================================================================
//                         TIME STEPPER EXCEPTIONS
//==============================================================================
// A TimeStepper may also return with exceptions thrown by its Integrator.

/** The %TimeStepper invoked an Event handler ("change" action) that failed. **/
class TimeStepper::EventChangeActionFailed : public Exception::Base {
public:
    EventChangeActionFailed(const char* fn, int ln, Real t, const char* msg) 
    :   Base(fn,ln) {
        setMessage("TimeStepper halted at time=" + String(t) + 
                   " because an EventAction failed, with message:\n" +
                   "  '" + msg + "'.");
    }
};

} // namespace SimTK

#endif // SimTK_SIMMATH_TIMESTEPPER_H_
