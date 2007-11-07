#ifndef SimTK_SIMMATH_INTEGRATOR_H_
#define SimTK_SIMMATH_INTEGRATOR_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmath(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-7 Stanford University and the Authors.         *
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

/** @file
 * This is the header file that user code should include to pick up the 
 * SimTK Simmath numerical integration tools.
 */

#include "SimTKcommon.h"
#include "simmath/internal/common.h"

namespace SimTK {

    
class IntegratorRep;

/**
 * Given a continuous system of differential equations for state variables y, and
 * optionally a manifold (set of algebraic equations) on which the solution must lie,
 * an Integrator object will advance that system through time. If the full system
 * is continuous, this is sufficient. However, most interesting systems consist
 * of both continuous and discrete equations, in
 * which case the overall time advancement is handled by a TimeStepper object which
 * will use an Integrator as an "inner loop" for advancing the system across 
 * continuous intervals between discrete updates. In that case, in addition to 
 * continuous state variables y there will be a set of discrete variables d which
 * are held constant during an integration interval.
 *
 * The continuous part of the system is an ODE-on-a-manifold
 * system suitable for solution via coordinate projection[1], structured like
 * this:
 *         (1)  y' = f[d](t,y)        differential equations
 *         (2)  0  = c[d](t,y)        algebraic equations (manifold is c=0)
 *         (3)  e  = e[d](t,y)        event triggers (watch for zero crossings)
 * with initial conditions t0,y0 such that c=0. The [d] is a reminder that
 * the overall system is dependent on discrete variables d as well as y, but
 * d cannot change during a continuous interval.
 * 
 * By "ODE on a manifold" we mean that the ODE (1) automatically satisfies
 * the condition that IF c==0, THEN c'=0, where
 *     c'=partial(c)/partial(t) + [partial(c)/partial(y)]*y'.
 * This is a less stringent condition than an ODE with "first integral" invariant,
 * in which c'=0 regardless of whether c=0.
 *
 * [1] Hairer, Lubich, Wanner, "Geometric Numerical Integration: Structure-Preserving
 * Algorithms for Ordinary Differential Equations", 2nd ed., section IV.4, pg 109ff,
 * Springer, 2006.
 *
 * The discrete variables d are updated by the time stepper upon occurence of specific
 * events, which terminate a continuous interval. The Integrator detects these
 * using a set of scalar-valued event trigger functions shown as equation (3) above. 
 * An event trigger function for a particular event must be designed so that it has a
 * zero crossing when the event occurs. The integrator can thus watch for sign changes
 * in event triggers and terminate the current step when a zero crossing occurs,
 * notifying the system and giving it a chance to handle the event; that is, 
 * update its state variables discontinuously.
 *
 * The zero crossings of continuous event trigger functions will
 * be isolated quickly; discontinuous ones have to be "binary chopped" which
 * is more expensive but they will still work.
 *
 * We are given a set of weights W for the y's, and a set of tolerances T
 * for the constraint errors. Given an accuracy specification (like 0.1%),
 * the integrators here are expected to solve for y(t) such that the
 * local error |W*y|_RMS <= accuracy, and |T*c(t,y)|_RMS <= accuracy at
 * all times.
 *
 * TODO: isolation tolerances for witnesses; dealing with simultaneity.
 *
 */
class SimTK_SIMMATH_EXPORT Integrator {
public:
    Integrator() : rep(0) { }

    // These are the exceptions that can be thrown by this class.

    class InitializationFailed;
    class StepSizeTooSmall;
    class StepFailed;
    class TriedToAdvancePastFinalTime;
    class CantAskForEventInfoWhenNoEventTriggered;

    const char* getMethodName();
    int         getMethodMinOrder();
    int         getMethodMaxOrder();
    bool        methodHasErrorControl();

    // Supply the integrator with a starting state. This is *copied*
    // into the integrator's internally maintained current state; 
    // subsequent changes to the State object passed in here will not
    // affect the integration.
    void initialize(const State&);

    // After the TimeStepper has made a discontinuous change to the Integrator's
    // "advanced state" we need to reinitialize the Integrator. However, event handlers
    // can do varying amounts of damage and some events will require no
    // reinitialization, or minimal reinitialization, depending on details
    // of the particular integration method. So after the handler has 
    // mangled our State, we tell the Integrator the lowest Stage which was
    // changed and allow the Integrator to figure out how much reinitialization
    // to do. For example, if Stage::Model has not been altered then
    // the Integrator will not have to reallocate any of its internal
    // data structures since sizes will be unchanged.
    // TODO: actually this should not be allowed if Model or Topology
    // stage changes have been made since the integrator would have no
    // way to create a reasonable state. Those changes require that
    // initialize() be called again, providing a workable Model-stage
    // state.
    // If 'shouldTerminate" is passed in true, the integrator will wrap
    // things up and report that the end of the simulation has been reached.
    void reinitialize(Stage, bool shouldTerminate);

    // Return a state which satisfies the caller's step request. This 
    // may be an interpolated value earlier than getAdvancedState().
    const State& getState() const;
    Real         getTime() const {return getState().getTime();}

    // This is mostly for prurient interest; it shouldn't matter. This
    // says whether getState() above will return an interpolated state
    // or just the same thing as getAdvancedState() does.
    bool isStateInterpolated() const;

    // Return the state representing the trajectory point to which the
    // integrator has irreversibly advanced. This may be later than the
    // state return by getState().
    const State& getAdvancedState() const;
    Real         getAdvancedTime() const {return getAdvancedState().getTime();}

    State& updAdvancedState();

    // Report the internally-maintained quantities used for accuracy control.
    // These may or may not be the same as what the DynamicSystem would return
    // at the current state. 
    Real getAccuracyInUse() const;
    Real getConstraintToleranceInUse() const;
    Real getTimeScaleInUse() const;
    const Vector& getStateWeightsInUse() const;
    const Vector& getConstraintWeightsInUse() const;

    // When a step is successful, it will return an indication of what
    // caused it to stop where it did. When unsuccessful it will throw
    // an exception so you won't see any return value.
    //
    // When return of control is due ONLY to reaching a report time,
    // (status is ReachedReportTime) the integrator's getState() method may return an
    // interpolated value at an earlier time than its getAdvancedState()
    // method would return. For the other returns, and whenever the report
    // time coincides with the end of an internal step, getState() and
    // getAdvancedState() will be identical.
    //
    // Note: we ensure algorithmically that no report time, scheduled time,
    // or final time t can occur *within* an event window, that is, we will
    // never have t_low < t < t_high for any interesting t. Further, t_report,
    // t_scheduled and t_final can coincide with t_high but only t_report can
    // be at t_low. The interior of t_low:t_high is a "no mans land" where we
    // don't understand the solution, so must be avoided.

    enum SuccessfulStepStatus {
        ReachedReportTime    =1, // stopped only to report; might be interpolated
        ReachedEventTrigger  =2, // localized an event; this is the *before* state (interpolated)
        ReachedScheduledEvent=3, // reached the limit provided in stepTo() (scheduled event)
        TimeHasAdvanced      =4, // user requested control whenever an internal step is successful
        ReachedStepLimit     =5, // took a lot of internal steps but didn't return control yet
        EndOfSimulation      =6, // termination; don't call again
        StartOfContinuousInterval=7, // the beginning of a continuous interval: either the start of the simulation, or t_high
                                     // after an event handler has modified the state.

        InvalidSuccessfulStepStatus = -1
    };
    static String successfulStepStatusString(SuccessfulStepStatus);

    SuccessfulStepStatus stepTo(Real reportTime, Real timeLimit=Infinity);
    SuccessfulStepStatus stepBy(Real interval, Real timeLimit=Infinity);


    // The following methods are callable only when stepTo() or stepBy() returns
    // "ReachedEventTrigger".
    Vec2               getEventWindow()         const; // w==(getTime(),getAdvancedTime()]
    const Array<int>&  getTriggeredEvents()     const; // indices corresponding to event trigger functions
    const Array<Real>& getEstimatedEventTimes() const; // all in w==(tLow,tHigh]
    const Array<EventStatus::EventTrigger>&
                       getEventTransitionsSeen() const;


        // TERMINATION //

    enum TerminationReason {
        ReachedFinalTime                 = 1,
        AnUnrecoverableErrorOccurred     = 2,
        EventHandlerRequestedTermination = 3,

        InvalidTerminationReason         = -1
    };

    // Don't call stepTo() or stepBy() again if this returns true.
    bool isSimulationOver() const;

    // Callable only when EndOfSimulation has been reached (isSimulationOver()==true).
    TerminationReason getTerminationReason() const;

    // Statistics (mutable)
    void resetAllStatistics();                // reset all stats to zero

    // What was the size of the first successful step after the last initialize() call?
    Real getActualInitialStepSizeTaken() const;

    // What was the size of the most recent successful step?
    Real getPreviousStepSizeTaken() const;

    // What step size will be attempted first on the next step() call?
    Real getPredictedNextStepSize() const;

    long getNStepsAttempted() const;
    long getNStepsTaken() const; 
    long getNRealizations() const;
    long getNProjections() const;
    long getNErrorTestFailures() const;
    long getNRealizationFailures() const;
    long getNProjectionFailures() const;

    // Control over integrator behavior. Not all integration methods
    // respond to all the controls. Defaults are shown.

    void setFinalTime(Real tFinal);         // Infinity
    void setInitialStepSize(Real hinit);    // default depends on method
    void setMinimumStepSize(Real hmin);     // default depends on method
    void setMaximumStepSize(Real hmax);     // Infinity

    void setAccuracy(Real accuracy);        
    void setRelativeTolerance(Real relTol);
    void setAbsoluteTolerance(Real absTol);
    void setConstraintTolerance(Real consTol);

    // If this many internal steps occur during stepTo() before reaching
    // the report time, we'll return control with a ReachedStepLimit status.
    // nSteps <= 0 means no limit.
    void setInternalStepLimit(int nSteps);  // 500

    void setReturnEveryInternalStep(bool shouldReturn); // false

    void setProjectEveryStep(bool forceProject);
    void setAllowInterpolation(bool shouldInterpolate);
    void setProjectInterpolatedStates(bool shouldProject);

protected:
    const IntegratorRep& getRep() const {assert(rep); return *rep;}
    IntegratorRep&       updRep()       {assert(rep); return *rep;}

    // opaque implementation for binary compatibility
    IntegratorRep* rep;
    friend class IntegratorRep;
};

} // namespace SimTK

#endif // SimTK_SIMMATH_INTEGRATOR_H_
