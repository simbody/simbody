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
 * An Integrator is an object that can simulate the behavior of a System through
 * time.  This is an abstract class.  Subclasses implement a variety of different
 * integration methods, which vary in their speed, accuracy, stability, and various
 * other properties.
 * 
 * <h3>Usage</h3>
 * 
 * An Integrator is most often used in combination with a TimeStepper.  The
 * TimeStepper automates much of the work of using an Integrator: invoking it repeatedly to
 * advance time, calling event handlers, reintializing the Integrator as necessary,
 * and so on.  A typical use of an Integrator generally resembles the following:
 * 
 * <pre>
 * // Instantiate an Integrator subclass which is appropriate for your problem.
 * VerletIntegrator integ(system);
 * // Set configuration options on the Integrator.
 * integ.setAccuracy(1e-3);
 * // Create a TimeStepper and use it to run the simulation.
 * TimeStepper stepper(system, integ);
 * stepper.initialize(initialState);
 * stepper.stepTo(finalTime);
 * </pre>
 * 
 * <h3>Mathematical Overview</h3>
 * 
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

    /// Get the name of this integration method
    const char* getMethodName();
    /// Get the minimum order this Integrator may use
    int         getMethodMinOrder();
    /// Get the maximum order this Integrator may use
    int         getMethodMaxOrder();
    /// Get whether this Integrator provides error control.  An error controlled Integrator will dynamically adjust its
    /// step size to maintain the level of accuracy specified with setAccuracy().  An Integrator which does not provide
    /// error control cannot do this, and will usually ignore the value specified with setAccuracy().
    bool        methodHasErrorControl();

    /// Supply the integrator with a starting state.  This must be called before the first call to stepBy() or stepTo().
    /// The specified state is copied into the Integrator's internally maintained current state; subsequent changes to
    /// the State object passed in here will not affect the integration.
    void initialize(const State& state);

    /// After an event handler has made a discontinuous change to the Integrator's
    /// "advanced state", this method must be called to reinitialize the Integrator.
    /// Event handlers can do varying amounts of damage and some events will require no
    /// reinitialization, or minimal reinitialization, depending on details
    /// of the particular integration method. So after the handler has 
    /// mangled our State, we tell the Integrator the lowest Stage which was
    /// changed and allow the Integrator to figure out how much reinitialization
    /// to do.
    ///
    /// If 'shouldTerminate" is passed in true, the Integrator will wrap
    /// things up and report that the end of the simulation has been reached.
    void reinitialize(Stage stage, bool shouldTerminate);

    /// Return a State corresponding to the "current" time at the end of the last call to
    /// stepTo() or stepBy(). This may be an interpolated value earlier than getAdvancedState().
    const State& getState() const;
    /// Get the time of the current State.  This is equivalent to calling getState().getTime().
    Real         getTime() const {return getState().getTime();}

    /// Get whether getState() will return an interpolated state or just the same thing as getAdvancedState() does.
    /// In most cases, you should not have reason to care whether the state is interpolated or not.
    bool isStateInterpolated() const;

    /// Return the state representing the trajectory point to which the
    /// integrator has irreversibly advanced. This may be later than the
    /// state return by getState().
    const State& getAdvancedState() const;
    /// Get the time of the advanced State.  This is equivalent to calling getAdvancedState().getTime().
    Real         getAdvancedTime() const {return getAdvancedState().getTime();}

    /// Get a non-const reference to the advanced state.
    State& updAdvancedState();

    /// Get the accuracy which is being used for error control.  Usually this is the same value that was
    /// specified to setAccuracy().
    Real getAccuracyInUse() const;
    /// Get the constraint tolerance which is being used for error control.  Usually this is the same value that was
    /// specified to setConstraintTolerance().
    Real getConstraintToleranceInUse() const;
    /// Get the characteristic timescale of the System being integrated.  This is calculated by calling
    /// calcTimescale() on the System.
    Real getTimeScaleInUse() const;
    /// Get the vector of weights being used to normalize errors in the state variables.  This is calculated by calling
    /// calcYUnitWeights() on the System.
    const Vector& getStateWeightsInUse() const;
    /// Get the vector of weights being used to normalize constraint errors.  This is calculated by calling
    /// calcYErrUnitTolerances() on the System.
    const Vector& getConstraintWeightsInUse() const;

    /// When a step is successful, it will return an indication of what
    /// caused it to stop where it did. When unsuccessful it will throw
    /// an exception so you won't see any return value.
    ///
    /// When return of control is due ONLY to reaching a report time,
    /// (status is ReachedReportTime) the integrator's getState() method may return an
    /// interpolated value at an earlier time than its getAdvancedState()
    /// method would return. For the other returns, and whenever the report
    /// time coincides with the end of an internal step, getState() and
    /// getAdvancedState() will be identical.
    ///
    /// Note: we ensure algorithmically that no report time, scheduled time,
    /// or final time t can occur *within* an event window, that is, we will
    /// never have t_low < t < t_high for any interesting t. Further, t_report,
    /// t_scheduled and t_final can coincide with t_high but only t_report can
    /// be at t_low. The interior of t_low:t_high is a "no man's land" where we
    /// don't understand the solution, so must be avoided.
    enum SuccessfulStepStatus {
        /// stopped only to report; might be interpolated
        ReachedReportTime    =1,
        /// localized an event; this is the *before* state (interpolated)
        ReachedEventTrigger  =2,
        /// reached the limit provided in stepTo() (scheduled event)
        ReachedScheduledEvent=3,
        /// user requested control whenever an internal step is successful
        TimeHasAdvanced      =4,
        /// took a lot of internal steps but didn't return control yet
        ReachedStepLimit     =5,
        /// termination; don't call again
        EndOfSimulation      =6,
        /// the beginning of a continuous interval: either the start of the simulation, or t_high after an event handler has modified the state.
        StartOfContinuousInterval=7,
        InvalidSuccessfulStepStatus = -1
    };
    /// Get a human readable description of the reason a step returned.
    static String successfulStepStatusString(SuccessfulStepStatus);

    /// Integrate the System until something happens which requires outside processing, and return a status code describing what happened.
    /// @param reportTime           the time of the next scheduled report
    /// @param scheduledEventTime   the time of the next scheduled event
    SuccessfulStepStatus stepTo(Real reportTime, Real scheduledEventTime=Infinity);
    /// Integrate the System until something happens which requires outside processing, and return a status code describing what happened.
    /// @param interval             the interval from the current time (as returned by getTime()) until the next scheduled report
    /// @param scheduledEventTime   the time of the next scheduled event
    SuccessfulStepStatus stepBy(Real interval, Real scheduledEventTime=Infinity);


    /// Get the window (tLow, tHigh] within which one or more events have been localized.  This may only be called
    /// when stepTo() or stepBy() has returned ReachedEventTrigger.
    Vec2 getEventWindow() const;
    /// Get the IDs of all events which have been localized within the event window.  This may only be called
    /// when stepTo() or stepBy() has returned ReachedEventTrigger.
    const std::vector<int>& getTriggeredEvents() const;
    /// Get the estimated times of all events which have been localized within the event window.  This may only be called
    /// when stepTo() or stepBy() has returned ReachedEventTrigger.
    const std::vector<Real>& getEstimatedEventTimes() const;
    /// Get EventTriggers describing the events which have been localized within the event window.  This may only be called
    /// when stepTo() or stepBy() has returned ReachedEventTrigger.
    const std::vector<EventStatus::EventTrigger>& getEventTransitionsSeen() const;


        // TERMINATION //

    /// Once the simulation has ended, getTerminationReason() may be called to find out what caused it to end.
    enum TerminationReason {
        /// The simulation reached the time specified by setFinalTime().
        ReachedFinalTime                 = 1,
        /// An error occurred which the Integrator was unable to handle.
        AnUnrecoverableErrorOccurred     = 2,
        /// An event handler function requested that the simulation terminate immediately.
        EventHandlerRequestedTermination = 3,
        /// This will be returned if getTerminationReason() is called before the simulation has ended.
        InvalidTerminationReason         = -1
    };

    /// Get whether the simulation has terminated.  If this returns true, you should not
    /// call stepTo() or stepBy() again.
    bool isSimulationOver() const;

    /// Get the reason the simulation terminated.  This should only be called if
    /// isSimulationOver() returns true.
    TerminationReason getTerminationReason() const;

    /// Reset all statistics to zero.
    void resetAllStatistics();

    /// Get the size of the first successful step after the last initialize() call.
    Real getActualInitialStepSizeTaken() const;

    /// Get the size of the most recent successful step.
    Real getPreviousStepSizeTaken() const;

    /// Get the step size that will be attempted first on the next call to stepTo() or stepBy().
    Real getPredictedNextStepSize() const;

    /// Get the total number of steps that have been attempted (successfully or unsuccessfully) since
    /// the last call to resetAllStatistics().
    long getNStepsAttempted() const;
    /// Get the total number of steps that have been successfully taken since the last call to resetAllStatistics().
    long getNStepsTaken() const; 
    /// Get the total number of state realizations that have been performed since the last call to resetAllStatistics().
    long getNRealizations() const;
    /// Get the total number of times a state has been projected since the last call to resetAllStatistics().
    long getNProjections() const;
    /// Get the number of attempted steps that have failed due to the error being unacceptably high since
    /// the last call to resetAllStatistics().
    long getNErrorTestFailures() const;
    /// Get the number of attempted steps that have failed due to an error when realizing the state since
    /// the last call to resetAllStatistics().
    long getNRealizationFailures() const;
    /// Get the number of attempted steps that have failed due to an error when projecting the state since
    /// the last call to resetAllStatistics().
    long getNProjectionFailures() const;

    /// Set the time at which the simulation should end.  The default is infinity.  Some integrators may
    /// not support this option.
    void setFinalTime(Real tFinal);
    /// Set the initial step size that should be attempted.  The default depends on the integration method.
    /// Some integrators may not support this option.
    void setInitialStepSize(Real hinit);
    /// Set the minimum step size that should ever be used.  The default depends on the integration method.
    /// Some integrators may not support this option.
    void setMinimumStepSize(Real hmin);
    /// Set the maximum step size that should ever be used.  The default depends on the integration method.
    /// Some integrators may not support this option.
    void setMaximumStepSize(Real hmax);
    
    /// Set the integrator to use a single fixed step size for all steps.  This is exactly equivalent
    /// to calling setInitialStepSize(), setMinimumStepSize(), and setMaximumStepSize(), passing the
    /// same value to each one.  This will therefore not work correctly if the integrator does not
    /// support minimum and/or maximum step sizes.
    void setFixedStepSize(Real stepSize);

    /// Set the overall accuracy that should be used for integration.  If the Integrator does not support
    /// error control (methodHasErrorControl() returns false), this may have no effect.
    void setAccuracy(Real accuracy);
    /// Set the relative error tolerance to be used for integration.  Some Integrators do not allow
    /// relative and absolute tolerances to be specified independently, in which case this will have
    /// no effect.  Even for Integrators which do support it, it is usually sufficient to call only
    /// setAccuracy(), and allow the Integrator to select relative and absolute tolerances accordingly.
    void setRelativeTolerance(Real relTol);
    /// Set the absolute error tolerance to be used for integration.  Some Integrators do not allow
    /// relative and absolute tolerances to be specified independently, in which case this will have
    /// no effect.  Even for Integrators which do support it, it is usually sufficient to call only
    /// setAccuracy(), and allow the Integrator to select relative and absolute tolerances accordingly.
    void setAbsoluteTolerance(Real absTol);
    /// Set the tolerance within which constraints must be satisfied.
    void setConstraintTolerance(Real consTol);

    /// Set the maximum number of steps that may be taken within a single call to stepTo() or stepBy().
    /// If this many internal steps occur before reaching the report time, it will return control with
    /// a ReachedStepLimit status.  If nSteps <= 0, the number of steps will be unlimited.
    void setInternalStepLimit(int nSteps);

    /// Set whether the Integrator should return from stepTo() or stepBy() after every internal step,
    /// even if no event has occurred and the report time has not been reached.  The default is false.
    void setReturnEveryInternalStep(bool shouldReturn);

    /// Set whether the system should be projected back to the constraint manifold after every step.
    /// If this is false, projection will only be performed when the constraint error exceeds the
    /// allowed tolerance.
    void setProjectEveryStep(bool forceProject);    
    /// Set whether the Integrator is permitted to return interpolated states for reporting purposes
    /// which may be less accurate than the "real" states that form the trajectory.  Setting this to
    /// true may significantly affect performance, since the Integrator will be forced to decrease its
    /// step size at every scheduled reporting time.
    ///
    /// This option is generally only meaningful if interpolated states are less accurate than other
    /// states on the trajectory.  If an Integrator can produce interpolated states that have the same
    /// accuracy as the rest of the trajectory, it may ignore this option.
    void setAllowInterpolation(bool shouldInterpolate);
    /// Set whether interpolated states should be projected back to the constraint manifold after
    /// interpolation is performed.
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
