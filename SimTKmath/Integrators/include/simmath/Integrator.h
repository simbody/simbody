#ifndef SimTK_SIMMATH_INTEGRATOR_H_
#define SimTK_SIMMATH_INTEGRATOR_H_

/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-15 Stanford University and the Authors.        *
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

/** @file
 * This is the header file that user code should include to pick up the 
 * SimTK Simmath numerical integration tools.
 */

#include "SimTKcommon.h"
#include "simmath/internal/common.h"

namespace SimTK {

    
class IntegratorRep;

//==============================================================================
//                               INTEGRATOR
//==============================================================================
/** An %Integrator is a Study that can advance the State of a System as it
evolves continuously through time, and detect discrete events. This is an 
abstract handle class. Subclasses implement a variety of different integration 
methods, which vary in their speed, accuracy, stability, and various other 
properties.

<h3>Usage</h3>

An %Integrator can be used alone, but is most often used in combination with a 
TimeStepper. A TimeStepper automates much of the work of using an %Integrator: 
invoking it repeatedly to advance time, calling event handlers and reporters, 
reintializing the %Integrator as necessary, and so on. A typical use of an 
%Integrator generally resembles the following:

@code
    // Instantiate an Integrator subclass which is appropriate for your problem.
    RungeKuttaMersonIntegrator integ(system);
    // Set configuration options on the Integrator.
    integ.setAccuracy(1e-3);
    // Create a TimeStepper and use it to run the simulation.
    TimeStepper stepper(integ);
    stepper.initialize(initialState);
    stepper.stepTo(finalTime);
@endcode

<h3>Mathematical Overview</h3>

Given a continuous system of differential equations for state variables y, and
optionally a manifold (set of algebraic equations) on which the solution must
lie, an %Integrator object will advance that system through time. If the full 
system is continuous, this is sufficient. However, most interesting systems 
consist of both continuous and discrete equations, in which case the overall 
time advancement is handled by a TimeStepper object which will use an 
%Integrator as an "inner loop" for advancing the system across continuous 
intervals between discrete updates. In that case, in addition to continuous 
state variables y there will be a set of discrete variables d which are held 
constant during an integration interval.

The continuous part of the system is an ODE-on-a-manifold system suitable for 
solution via coordinate projection[1], structured like this: <pre>
    (1)  y' = f(d;t,y)        differential equations
    (2)  0  = c(d;t,y)        algebraic equations (manifold is c=0)
    (3)  e  = e(d;t,y)        event witnesses (watch for zero crossings)
</pre> with initial conditions t0,y0 such that c=0. The "d;" is a reminder that 
the overall system is dependent on discrete variables d as well as y, but d 
cannot change during a continuous interval.

By "ODE on a manifold" we mean that the ODE (1) automatically satisfies the 
condition that IF c==0, THEN c'=0, where <pre>
    c' = dc/dt = partial(c)/partial(t) + [partial(c)/partial(y)]*y'.
</pre>
This is a less stringent condition than an ODE with "first integral" invariant,
in which c'=0 regardless of whether c=0.

[1] Hairer, Lubich, Wanner, "Geometric Numerical Integration: 
Structure-Preserving Algorithms for Ordinary Differential Equations", 2nd ed., 
section IV.4, pg 109ff, Springer, 2006.

The discrete variables d are updated by the time stepper upon occurence of 
specific events, which terminate a continuous interval. The %Integrator detects
these using a set of scalar-valued event witness functions shown as equation 
(3) above. An event witness function for a particular event must be designed so
that it has a zero crossing when the event occurs. The integrator can thus 
watch for sign changes in event triggers and terminate the current step when a
zero crossing occurs, notifying the System and giving it a chance to handle the
event; that is, update its state variables discontinuously.

The zero crossings of continuous event witness functions will be isolated 
quickly; discontinuous ones have to be "binary chopped" which is more expensive
but they will still work.

We are given a set of weights W for the y's, a set of tolerances T for the
constraint errors, and localization requirements for the event witnesses. Given
an accuracy specification (like 0.1%), an %Integrator is expected to 
solve for y(t) such that the local error |W*yerr|_RMS <= accuracy, 
and |T*c(t,y)|_RMS <= accuracy at all times. They are also required to localize
witness function zero crossings to satisfy witness localization requirements.
**/
class SimTK_SIMMATH_EXPORT Integrator : public Study {
public:
    /** Default constructor creates an empty handle. **/
    Integrator() : rep(0) { }

    /** Destructor destroys the internally-maintained State, but not the
    System which is an independent object. **/
    ~Integrator();

    // These are the exceptions that can be thrown by this class.
    class InitializationFailed;
    class StepSizeTooSmall;
    class StepFailed;
    class TriedToAdvancePastFinalTime;
    class CantAskForEventInfoWhenNoEventTriggered;

    /** Get the name of this integration method. **/
    const char* getMethodName() const;
    /** Get the minimum integration order this %Integrator may use. **/
    int getMethodMinOrder() const;
    /** Get the maximum integration order this %Integrator may use. **/
    int getMethodMaxOrder() const;

    /** Get whether this %Integrator provides error control. An error controlled
    %Integrator will dynamically adjust its step size to maintain the level of 
    accuracy specified with setAccuracy(). An %Integrator which does not provide
    error control cannot do this, and will usually ignore the value specified 
    with setAccuracy(). **/
    bool methodHasErrorControl() const;

    /** Supply the integrator with a starting state. This must be called before
    the first call to stepBy() or stepTo(). The specified state is *copied* into
    the %Integrator's internally maintained current State; subsequent changes to
    the State object passed in here will not affect the integration. The 
    internal State will be realized, the System's Initialization event will be
    triggered and its actions performed on the internal state, statistics will
    be reset, and the %Integrator will be made ready to take its first step.

    @throws Integrator::InitializationFailed if something goes wrong, with
            a message providing more detail. **/
    void initialize(const State& initialState);

    /** Return a const reference to the System to which this %Integrator is
    being applied. **/
    const System& getSystem() const;

    /** Return a const reference to the State corresponding to the "current" 
    time at the end of the last call to stepTo() or stepBy(). This may be an 
    interpolated value earlier than getAdvancedState(), in which case the State
    object referenced here is *not* the same object as would be returned by
    getAdvancedState(). You should not expect this to return the same State
    object from call to call during a simulation. **/
    const State& getState() const;

    /** Get the time of the current State. This is equivalent to calling 
    getState().getTime(). **/
    double getTime() const {return getState().getTime();}

    /** Get whether getState() will return an interpolated state or just the 
    same thing as getAdvancedState() does. In most cases, you should not have 
    reason to care whether the state is interpolated or not. **/
    bool isStateInterpolated() const;

    /** Return a const reference to the state representing the trajectory point
    to which the integrator has irreversibly advanced. This may be later than 
    the State return by getState(), which may be an interpolation. Note that 
    this will not necessarily return a reference the same object as 
    getState(). **/
    const State& getAdvancedState() const;

    /** Get the time of the advanced State. This is equivalent to calling 
    getAdvancedState().getTime(). **/
    double getAdvancedTime() const {return getAdvancedState().getTime();}

    /** Get a writable reference to the advanced state. This is only available
    if you have write access to the %Integrator. **/
    State& updAdvancedState();

    /** Get the accuracy which is being used for error control. Usually this is
    the default accuracy or the value that was last specified with
    setAccuracy(). **/
    Real getAccuracyInUse() const;

    /** Get the constraint tolerance which is being used for error control.
    Usually this is the default tolerance based on the accuracy setting, or the
    value that was last specified with setConstraintTolerance(). **/
    Real getConstraintToleranceInUse() const;

    /** When a step is successful, it will return an indication of what caused 
    it to stop where it did. When unsuccessful it will throw an exception so you
    won't see any return value.
    
    When return of control is due ONLY to reaching a report time, (status is 
    ReachedReportTime) the integrator's getState() method may return an
    interpolated value at an earlier time than its getAdvancedState() method 
    would return. For the other returns, and whenever the report time coincides 
    with the end of an internal step, getState() and getAdvancedState() will be
    identical.
    
    Note: we ensure algorithmically that no report time, scheduled time, or
    final time t can occur *within* an event witness localization window, that 
    is, we will never have t_low < t < t_high for any interesting t. Further, 
    t_report, t_scheduled and t_final can coincide with t_high but only t_report
    can be at t_low. The interior of the interval t_low:t_high is a "no man's 
    land" where we don't understand the solution, so must be avoided. **/
    enum SuccessfulStepStatus {
        /** Stopped only to report; might be interpolated. **/
        ReachedReportTime    =1,
        /** Localized an event; this is the *before* state (interpolated). **/
        ReachedEventTrigger  =2,
        /** Reached the limit provided in stepTo() (scheduled event). **/
        ReachedScheduledEvent=3,
        /** User requested control whenever an internal step is successful. **/
        TimeHasAdvanced      =4,
        /** Took a lot of internal steps but didn't return control yet. **/
        ReachedStepLimit     =5,
        /** Termination; don't call again. **/
        EndOfSimulation      =6,
        /** The beginning of a continuous interval: either the start of the 
        simulation, or t_high after an event handler has modified the state. **/
        StartOfContinuousInterval=7,
        /** This will never be returned as a status. **/
        InvalidSuccessfulStepStatus = -1
    };
    /** Get a human readable description of the reason a step returned. **/
    static String getSuccessfulStepStatusString(SuccessfulStepStatus);

    /** Integrate the System until something happens which requires outside 
    processing, and return a status code describing what happened.
    @param reportTime           the time of the next scheduled report
    @param scheduledEventTime   the time of the next scheduled event **/
    SuccessfulStepStatus stepTo(double reportTime, 
                                double scheduledEventTime=Infinity);
    
    /** Integrate the System until something happens which requires outside 
    processing, and return a status code describing what happened.
    @param reportInterval       
        The interval from the current time (as returned by getTime()) until the
        next scheduled report.
    @param scheduledEventInterval   
        The interval from the current time (as returned by getTime() until the 
        next scheduled event. 

    This is equivalent to @code
        stepTo(getTime()+reportInterval, getTime()+scheduledEventInterval);
    @endcode
    **/
    SuccessfulStepStatus stepBy(double reportInterval, 
                                double scheduledEventInterval=Infinity);

        // TERMINATION //

    /** Once the simulation has ended, getTerminationReason() may be called to 
    find out what caused it to end. **/
    enum TerminationReason {
        /** The simulation reached the time specified by setFinalTime(). **/
        ReachedFinalTime                 = 1,
        /** An error occurred which the %Integrator was unable to handle. **/
        AnUnrecoverableErrorOccurred     = 2,
        /** An event handler function requested that the simulation terminate
        immediately. **/
        EventHandlerRequestedTermination = 3,
        /** This will be returned if getTerminationReason() is called before 
        the simulation has ended. **/
        InvalidTerminationReason         = -1
    };

    /** (Advanced) %Force the current integration Study to terminate. The 
    System's Termination event will be triggered and its actions performed on 
    the %Integrator's internal State, which may be examined after termination
    along with statistics. However, subsequent calls to stepTo() will fail 
    unless initialize() is called again. This will typically be invoked due
    to the %Integrator reaching a user-specified final time, or due to 
    a termination request or failure from an event handler. After this call,
    isSimulationOver() will return `true`. **/
    void terminate(TerminationReason reason);

    /** Get whether the simulation has terminated. If this returns true, you 
    should not call stepTo() or stepBy() again. **/
    bool isSimulationOver() const;

    /** Get the reason the simulation terminated. This should only be called if
    isSimulationOver() returns true. **/
    TerminationReason getTerminationReason() const;

    /** Get a human readable description of the termination reason. **/
    static String getTerminationReasonString(TerminationReason);

    /** Reset all statistics to zero. **/
    void resetAllStatistics();

    /** Get the size of the first successful step after the last 
    initialize() call. **/
    double getActualInitialStepSizeTaken() const;

    /** Get the size of the most recent successful step. **/
    double getPreviousStepSizeTaken() const;

    /** Return which of the continuous state variables y={q,u,z} had the worst 
    error and hence limited the size of the previous step taken. The result
    will be an invalid index value if there is no previous step or if the 
    integrator in use does not support this feature; you can check 
    with `isValid()`. **/
    SystemYIndex getPreviousStepWorstState() const;

    /** Get the step size that will be attempted first on the next call to 
    stepTo() or stepBy(). **/
    double getPredictedNextStepSize() const;

    /** Get the total number of steps that have been attempted (successfully or
    unsuccessfully) since the last call to resetAllStatistics(). **/
    int getNumStepsAttempted() const;
    /** Get the total number of steps that have been successfully taken since 
    the last call to resetAllStatistics(). **/
    int getNumStepsTaken() const; 
    /** Get the total number of state realizations that have been performed 
    since the last call to resetAllStatistics(). **/
    int getNumRealizations() const;
    /** Get the total number of times a state positions Q have been projected
    since the last call to resetAllStatistics(). **/
    int getNumQProjections() const;
    /** Get the total number of times a state velocities U have been projected
    since the last call to resetAllStatistics(). **/
    int getNumUProjections() const;
    /** Get the total number of times a state has been projected (counting 
    both Q and U projections) since the last call to resetAllStatistics(). **/
    int getNumProjections() const;
    /** Get the number of attempted steps that have failed due to the error 
    being unacceptably high since the last call to resetAllStatistics(). **/
    int getNumErrorTestFailures() const;
    /** Get the number of attempted steps that failed due to non-convergence of
    internal step iterations. This is most common with iterative methods 
    but can occur if for some reason a step can't be completed. It is reset
    to zero by resetAllStatistics. **/
    int getNumConvergenceTestFailures() const;
    /** Get the number of attempted steps that have failed due to an error when 
    realizing the state since the last call to resetAllStatistics(). **/
    int getNumRealizationFailures() const;
    /** Get the number of attempted steps that have failed due to an error 
    when projecting the state positions (Q) since the last call to 
    resetAllStatistics(). **/
    int getNumQProjectionFailures() const;
    /** Get the number of attempted steps that have failed due to an error 
    when projecting the state velocities (U) since the last call to 
    resetAllStatistics(). **/
    int getNumUProjectionFailures() const;
    /** Get the number of attempted steps that have failed due to an error 
    when projecting the state (either a Q- or U-projection) since
    the last call to resetAllStatistics(). **/
    int getNumProjectionFailures() const;
    /** For iterative methods, get the number of internal step iterations in 
    steps that led to convergence (not necessarily successful steps). Reset to 
    zero by resetAllStatistics(). **/
    int getNumConvergentIterations() const;
    /** For iterative methods, get the number of internal step iterations in 
    steps that did not lead to convergence. Reset to zero by 
    resetAllStatistics(). **/
    int getNumDivergentIterations() const;
    /** For iterative methods, this is the total number of internal step 
    iterations taken regardless of whether those iterations led to convergence 
    or to successful steps. This is the sum of the number of convergent and 
    divergent iterations which are available separately. **/
    int getNumIterations() const;

    /** Set the time at which the simulation should end. The default is 
    Infinity. Some integrators may not support this option. **/
    void setFinalTime(Real tFinal);


    /** Set the minimum step size that should ever be used. The default depends 
    on the integration method. Some integrators may not support this option. **/
    void setMinimumStepSize(Real hmin);

    /** Set the maximum step size that should ever be used. The default depends
    on the integration method. Some integrators may not support this option. **/
    void setMaximumStepSize(Real hmax);
    
    /** Set the integrator to use a single fixed step size for all steps. This 
    is exactly equivalent to calling setInitialStepSize(), setMinimumStepSize(),
    and setMaximumStepSize(), passing the same value to each one. This will 
    therefore not work correctly if the integrator does not support minimum 
    and/or maximum step sizes. **/
    void setFixedStepSize(Real stepSize);

    /** Set the overall accuracy that should be used for integration. If the 
    %Integrator does not support error control (methodHasErrorControl() 
    returns `false`), this may have no effect. **/
    void setAccuracy(Real accuracy);

    /** Set the tolerance within which constraints must be satisfied. **/
    void setConstraintTolerance(Real consTol);

    /** Set the maximum number of steps that may be taken within a single call
    to stepTo() or stepBy(). If this many internal steps occur before 
    reaching the report time, it will return control with a ReachedStepLimit
    status. If nSteps <= 0, the number of steps will be unlimited. **/
    void setInternalStepLimit(int nSteps);

    /** Set whether the %Integrator should return from stepTo() or stepBy() 
    after every internal step, even if no event has occurred and the report time
    has not been reached. The default is `false`. If you are using a TimeStepper
    and want the TimeStepper also to return every step, you must *also* tell it
    to do so via TimeStepper::setReportAllSignificantStates(). **/
    void setReturnEveryInternalStep(bool shouldReturn);

    /** @name                  Event Handling 
    These methods are needed for dealing with events that the %Integrator
    detected. Normally these will be used only by a TimeStepper but you can
    use them directly if you want. **/
    /**@{**/
    /** Get the window `(tLow, tHigh]` within which one or more event witness
    zero crossings have been localized. This may be called only when 
    stepTo() or stepBy() has returned status ReachedEventTrigger.
    @returns Vec2(tLow,tHigh) **/
    Vec2 getEventWindow() const;

    /** Get pointers to the event witnesses whose zero crossings have been
    localized to within the event window returned by getEventWindow(). This may 
    be called only when stepTo() or stepBy() has returned status 
    ReachedEventTrigger. **/
    const Array_<const EventTrigger::Witness*>& getTriggeredWitnesses() const;

    /** Get the estimated times for the zero crossings of each witness that
    is returned by getTriggeredWitnesses(), with corresponding indexing. This 
    may be called only when stepTo() or stepBy() has returned status 
    ReachedEventTrigger. **/
    const Array_<Real>& getEstimatedTriggerTimes() const;

    /** Get the observed zero crossing direction (rising or falling) for each 
    witness that is returned by getTriggeredWitnesses(), with corresponding 
    indexing. This may only be called when stepTo() or stepBy() has returned 
    status ReachedEventTrigger. **/
    const Array_<Event::TriggerDirection>& getWitnessTransitionsSeen() const;

    /** After an event handler has made a discontinuous change to the 
    %Integrator's "advanced state", this method must be called to reinitialize 
    the %Integrator. Event handlers can do varying amounts of damage and some 
    events will require no reinitialization, or minimal reinitialization, 
    depending on details of the particular integration method. So after the 
    handler has mangled our State, we tell the %Integrator the lowest Stage 
    which was changed and allow the %Integrator to figure out how much 
    reinitialization to do.
    
    If "shouldTerminate" is passed in true, the %Integrator will wrap
    things up and report that the end of the simulation has been reached. **/
    void reinitialize(Stage stage, bool shouldTerminate);
    /**@}**/

    /** @name               Advanced/Debugging
    You probably won't need to use these. **/
    /**@{**/
    /** (Advanced) Set the initial step size that should be attempted. Normally
    an %Integrator will choose a reasonable default, try it, and then adjust 
    quickly to a good step size based on error control. The default depends on 
    the integration method and the System default time scale. Some integrators 
    may not support this option. **/
    void setInitialStepSize(Real hinit);

    /** (Advanced) Use infinity norm (maximum absolute value) instead of 
    default RMS norm to evaluate whether accuracy has been achieved for 
    states and for constraint tolerance. The infinity norm is more strict 
    but may permit use of a looser accuracy request. **/
    void setUseInfinityNorm(bool useInfinityNorm);

    /** (Advanced) Are we currently using the infinity norm? **/
    bool isInfinityNormInUse() const;

    /** (Advanced) Set whether the system should be projected back to the 
    constraint manifold after every step. If this is false, projection will only
    be performed when the constraint error exceeds the allowed tolerance. **/
    void setProjectEveryStep(bool forceProject);   

    /** (Advanced) Set whether the %Integrator is permitted to return 
    interpolated states for reporting purposes which may be less accurate than 
    the "real" states that form the trajectory. Setting this to `false` may 
    significantly affect performance, since the %Integrator will be forced 
    to decrease its step size at every scheduled reporting time.
    
    This option is generally only meaningful if interpolated states are 
    less accurate than other states on the trajectory. If an %Integrator can
    produce interpolated states that have the same accuracy as the rest of 
    the trajectory, it may ignore this option. **/
    void setAllowInterpolation(bool shouldInterpolate);

    /** (Advanced) Constraint projection may use an out-of-date iteration
    matrix for efficiency. You can force strict use of a current iteration
    matrix recomputed at each iteration if you want. **/
    void setForceFullNewton(bool forceFullNewton);

    /** (Advanced) Set whether interpolated states should be projected back to 
    the constraint manifold after interpolation is performed. The default is
    `true`. **/
    void setProjectInterpolatedStates(bool shouldProject);
    /**@}**/

    /** Return `true` if the given Study is an %Integrator. **/
    static bool isA(const Study& study)
    {   return dynamic_cast<const Integrator*>(&study) != nullptr; }

    /** Downcast the given const Study to a const %Integrator, assuming that it 
    actually is an %Integrator. At least in Debug builds, an exception will be 
    thrown if the Study is not an %Integrator. **/
    static const Integrator& downcast(const Study& study)
    {   return SimTK_DYNAMIC_CAST_DEBUG<const Integrator&>(study); }

    /** Downcast the given writable Study to a writable reference to 
    %Integrator, assuming that it actually is an %Integrator. At least in Debug 
    builds, an exception will be thrown if the Study is not an %Integrator. **/
    static Integrator& updDowncast(Study& study)
    {   return SimTK_DYNAMIC_CAST_DEBUG<Integrator&>(study); }

    /** @cond **/ // hide from doxygen
    /** (Deprecated) Use getSuccessfulStepStatusString() instead. **/
    static String successfulStepStatusString(SuccessfulStepStatus stat)
    {   return getSuccessfulStepStatusString(stat); }
    /** @endcond **/

protected:
    const IntegratorRep& getRep() const {assert(rep); return *rep;}
    IntegratorRep&       updRep()       {assert(rep); return *rep;}

    // opaque implementation
    IntegratorRep* rep;

private:
    // Implement the Study interface.
    const System& getSystemVirtual() const override
    {   return getSystem(); }
    const State& getCurrentStateVirtual() const override
    {   return getState(); } 
    const State& getInternalStateVirtual() const override
    {   return getAdvancedState(); }  
    State& updInternalStateVirtual() override
    {   return updAdvancedState(); } 
    Real getAccuracyInUseVirtual() const override
    {   return getAccuracyInUse(); }
    Real getConstraintToleranceInUseVirtual() const override
    {   return getConstraintToleranceInUse(); }

};

//==============================================================================
//                         INTEGRATOR EXCEPTIONS
//==============================================================================

/** Something went wrong during the initialize() call. **/
class Integrator::InitializationFailed : public Exception::Base {
public:
    InitializationFailed(const char* fn, int ln, const char* msg) 
    :   Base(fn,ln) {
        setMessage("Integrator initialization failed apparently because:\n"
                    + String(msg));
    }
};

/** The %Integrator was unable to achieve some required goal, such as the
specified accuracy, after reducing the step size repeatedly until it went
below the allowed minimum step size. **/
class Integrator::StepSizeTooSmall : public Exception::Base {
public:
    StepSizeTooSmall(const char* fn, int ln, Real t, Real h, Real hmin) 
    :   Base(fn,ln) {
        setMessage("At time=" + String(t) + 
        " the integrator failed to take a step with step size "
        + String(h) + " which was already at or below the minimum allowed size " 
        + String(hmin));
    }
};

/** The current step failed for some reason other than being reduced to
below the minimum. **/
class Integrator::StepFailed : public Exception::Base {
public:
    StepFailed(const char* fn, int ln, Real t, const char* msg) 
    :   Base(fn,ln) {
        setMessage("Integrator step failed at time " + String(t) 
                    + " apparently because:\n" + String(msg));
    }
};
/** Once the final time has been reached, a subsequent call to stepTo()
will thrown this exception. initialize() must be called again before the
%Integrator can be used. **/
class Integrator::TriedToAdvancePastFinalTime : public Exception::Base {
public:
    TriedToAdvancePastFinalTime
        (const char* fn, int ln, Real t, const char* msg) 
    :   Base(fn,ln) {
        setMessage("stepTo() called at time " + String(t) 
                    + " but the final step has already been reported."
                    + " You must call initialize() again to restart.");
    }
};

/** Some methods provide information about events that occurred in the
current step; those cannot be called if no event occurred. **/
class Integrator::CantAskForEventInfoWhenNoEventTriggered 
:   public Exception::Base {
public:
    CantAskForEventInfoWhenNoEventTriggered
        (const char* fn, int ln, const char* methodname, Real t) 
    :   Base(fn,ln) {
        setMessage("Method Integrator::" + String(methodname)
                    + "() was called at time " + String(t) 
                    + " but no event had triggered.");
    }
};


} // namespace SimTK

#endif // SimTK_SIMMATH_INTEGRATOR_H_
