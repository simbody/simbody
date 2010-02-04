#ifndef SimTK_SIMMATH_INTEGRATOR_REP_H_
#define SimTK_SIMMATH_INTEGRATOR_REP_H_

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
 * This is the declaration of the abstract IntegratorRep class which
 * represents the implementation of the Integrator class, and DAESystemRep
 * which implements the Integrator::DAESystem class.
 */

#include "SimTKcommon.h"
#include "simmath/Integrator.h"
//#include "SimTKcpodes/DynamicSystem.h"

//#include "DynamicSystemRep.h"

#include <exception>
#include <limits>
#include <algorithm>
#include <vector>

namespace SimTK {
//class DynamicSystemRep;

    ///////////////////////////
    // INTEGRATOR EXCEPTIONS //
    ///////////////////////////

class Integrator::InitializationFailed : public Exception::Base {
public:
    InitializationFailed(const char* fn, int ln, const char* msg) : Base(fn,ln) {
        setMessage("Integrator initialization failed apparently because:\n"
                   + String(msg));
    }
};


class Integrator::StepSizeTooSmall : public Exception::Base {
public:
    StepSizeTooSmall(const char* fn, int ln, Real t, Real h, Real hmin) : Base(fn,ln) {
        setMessage("At time=" + String(t) + 
                   " the integrator failed to take a step with step size "
                  + String(h) + " which was already at or below the minimum allowed size " 
                  + String(hmin));
    }
};

class Integrator::StepFailed : public Exception::Base {
public:
    StepFailed(const char* fn, int ln, Real t, const char* msg) : Base(fn,ln) {
        setMessage("Integrator step failed at time " + String(t) + " apparently because:\n"
                   + String(msg));
    }
};

class Integrator::CantAskForEventInfoWhenNoEventTriggered : public Exception::Base {
public:
    CantAskForEventInfoWhenNoEventTriggered
       (const char* fn, int ln, const char* methodname, Real t) : Base(fn,ln)
    {
        setMessage("Method Integrator::" + String(methodname)
                   + "() was called at time " + String(t) 
                   + " but no event had triggered.");
    }
};

// 
// Specification of the stepTo() method below:
//
// On entry we expect advancedState to be valid. Time may have already
// been advanced past the indicated reportTime in which case we'll
// be able to return an interpolated result immediately. We will never
// advance past scheduledEventTime or finalTime.
//
// The Integrator has a simple state machine to deal with the fact that the
// computation is interruptable due to the need for returning control
// to the caller (usually the time stepper) for various reasons. We enter
// this computation either due to a call to the routine stepTo(), or from
// within stepTo() after taking an internal step that didn't require a return.
// Here's a sketch:
//
//          advance (or report interpolated value)
//           ------
//          |      ^   
//          v      |no  
//   (TimeAdvanced)--return?->(StepReturned)--final?->(ReturnedFinalTime)
//        ^                                 |no     
//        <-------advance-------------------<
//
// In words, we take a step, we report that step if necessary (returning
// control), then take another step unless we just reported one at the
// final time. The state transition table below gives more detail about
// what exactly gets reported. "->Xyz" below means return control
// to caller of stepTo() with status Xyz, "INTERP t" means "interpolate
// back to time t".
//
// There are additional states below used to distinguish between triggered
// event conditions (which take some special handling) and the others.
// The only difference between state StepReturned and StepReturnedEvent
// is that some event-info query methods are allowed in the latter state.
//
// Current State     Condition         Action           Next State
// ----------------  ----------------  --------------  ----------------
// Uninitialized          any          throw error            X
// ReturnedFinalTime      any          throw error            X
//
// StepReturned or   t_adv>=t_final    ->Done         ReturnedFinalTime
// StepReturnedEvent    OTHERWISE      ADVANCE       TimeAdvanced[Event]
//
// TimeAdvancedEvent t_report<t_low    INTERP t_report  (unchanged)
//                                     ->Report
//
//                      OTHERWISE      INTERP t_low
//                  (event triggered)  ->Triggered    StepReturnedEvent
//
// TimeAdvanced      t_report<t_adv    INTERP t_report 
//                                     ->Report          (unchanged)
//
//                         else
//                   t_adv==t_sched    ->Scheduled      StepReturned
//
//                         else
//                   user wants control
//                     at end of step  ->EndOfStep      StepReturned
//
//                         else
//                   t_adv==t_report   ->Report         StepReturned
//
//                         else
//                   t_adv>=t_final    ->Final          StepReturned       
//
//                         else
//                   hit step limit    ->StepLimit      StepReturned
//                 
//                      OTHERWISE      ADVANCE        TimeAdvanced[Event]
// ----------------  ----------------  --------------  ----------------
//
// Notes:
// * for the above conditions to be sufficient, the integrator must
//   ensure that the interior of the interval (t_low,t_high) to
//   which an event trigger has been localized DOES NOT contain
//   t_report, t_sched, or t_final. It is OK if t_report is exactly
//   t_low, and it is OK if t_report, t_sched, or t_final is exactly
//   t_high.
// * the intent of this state diagram is to make sure that the
//   most significant end-of-step reason is reported, but that each
//   step is reported at most once. Note that when the return status
//   is "Done", the final trajectory point has *already* been reported.
//   Simultaneous occurrences must
//   be anticipated and handled by the caller (time stepper).
// * note that an interpolated report doesn't count as reporting the
//   step since we haven't reached t_advanced. But the interpolation
//   to t_low does since t_advanced==t_high in that case and there
//   is nothing in between t_low and t_high. Also, a report that
//   coincides with the end of an internal step counts as a step report.

class IntegratorRep {
public:
    IntegratorRep(Integrator* handle,
                  const System&);
    virtual ~IntegratorRep() { }
    // no default constructor, no copy or copy assign

    // The DynamicSystem must be successfully realized to Stage::Model before this
    // call. At this point the integrator can query the DynamicSystem about the
    // problem size and allocate appropriately-sized internal data structures.
    void initialize(const State&);
    void reinitialize(Stage, bool shouldTerminate);
    // We also give the concrete integration method a chance to initialize itself
    // after the generic initialization is done.
    virtual void methodInitialize(const State&) { }
    virtual void methodReinitialize(Stage stage, bool shouldTerminate) { }

    // The integrator has already been initialized. Take as many internal steps
    // as needed to get up to or past reportTime, but don't advance past the
    // next scheduled event time or the simulation final time.
    virtual Integrator::SuccessfulStepStatus 
        stepTo(Real reportTime, Real scheduledEventTime) = 0;

    bool isSimulationOver() const {
        return stepCommunicationStatus == FinalTimeHasBeenReturned;
    }
    
    // Get the reason the simulation ended.  This should only be invoked If isSimulationOver() returns true.
    Integrator::TerminationReason getTerminationReason() const {
        assert (isSimulationOver());
        return terminationReason;
    }

    // This represents the Integrator's finite state machine for dealing with 
    // communication of internal steps to the caller. After initialization, the
    // computation will be in state CompletedInternalStepNoEvent, but with t_prev=t_advanced=t_initial.
    // We will process this as though we had just taken a step so that various conditions will
    // result in the first stepTo() call returning immediately.
    // Note: event handling and subsequent partial reinitialization of the integrator MUST NOT
    // change the integrator's communication state.
    enum StepCommunicationStatus {
        CompletedInternalStepNoEvent,   // Time has advanced from t_prev to t_advanced; no event
                                        //   triggered. We have not yet returned to the caller with
                                        //   time t_advanced, although we may have returned at interpolated
                                        //   report times t_prev < t_report < t_advanced.
        CompletedInternalStepWithEvent, // Time has advanced from t_prev to t_advanced with an
                                        //   event trigger localized to the interval (t_low,t_high]
                                        //   where t_high==t_advanced. We have not yet returned to
                                        //   the caller with time t_low, although we may have
                                        //   returned at interpolated report times 
                                        //   t_prev < t_report < t_low.
        StepHasBeenReturnedNoEvent,     // We have already returned control to the caller at t_advanced
                                        //   with the strongest stopping reason given;
                                        //   no more returns are allowed for this step.
        StepHasBeenReturnedWithEvent,   // We have already returned control to the caller at t_low
                                        //   with "EventTriggered" as the stopping reason;
                                        //   no more returns are allowed for this step.
        FinalTimeHasBeenReturned,       // Any call to stepTo() when the integrator is already in this
                                        //   state is a fatal error.

        InvalidStepCommunicationStatus = -1
    };

    const State& getAdvancedState() const {
        return advancedState;
    }
    Real getAdvancedTime() const {return advancedState.getTime();}

    const State& getState() const {
        return useInterpolatedState ? interpolatedState : advancedState;
    }
    bool isStateInterpolated() const {return useInterpolatedState;}

    Real getAccuracyInUse() const {
        return accuracyInUse;
    }

    Real getConstraintToleranceInUse() const {
        return consTol;
    }

    Real getTimeScaleInUse() const {
        return timeScaleInUse;
    }

    const Vector& getStateWeightsInUse() const {
        return stateWeightsInUse;
    }

    const Vector& getConstraintWeightsInUse() const {
        return constraintWeightsInUse;
    }

    // What was the size of the first successful step after the last initialize() call?
    virtual Real getActualInitialStepSizeTaken() const = 0;

    // What was the size of the most recent successful step?
    virtual Real getPreviousStepSizeTaken() const = 0;

    // What step size will be attempted first on the next step() call?
    virtual Real getPredictedNextStepSize() const = 0;

    virtual long getNumStepsAttempted() const = 0;
    virtual long getNumStepsTaken() const = 0; 
    virtual long getNumErrorTestFailures() const = 0;
    virtual long getNumConvergentIterations() const = 0;
    virtual long getNumDivergentIterations() const = 0;
    virtual long getNumIterations() const = 0;

    virtual void resetMethodStatistics() {
    }

    // Cubic Hermite interpolation. See Hairer, et al. Solving
    // ODEs I, 2nd rev. ed., pg 190. Given (t0,y0,y0'),(t1,y1,y1')
    // with y0 and y1 at least 3rd order accurate,
    // we can obtain a 3rd order accurate interpolation yt for
    // a time t0 < t < t1 using Hermite interpolation. Let
    // f0=y0', f1=y1', h=t1-t0, d=(t-t0)/h (0<=d<=1).
    //   u(d)=y(t0+dh)=(1-d)y0 + dy1 + d(d-1)[(1-2d)(y1-y0) + (d-1)hf0 + dhf1]
    // Rearrange so we only have to through each array once:
    //   u(d)=cy0*y0 + cy1*y1 + cf0*f0 + cf1*f1
    //   cy0= 1 - d^2(3 - 2*d)
    //   cy1=     d^2(3 - 2*d) = (1-cy0)
    //   cf0=d*(d-1)^2*h = h*d*(d-1) * (d-1)
    //   cf1=d^2*(d-1)*h = h*d*(d-1) * d
    //   
    // Note: you can't get a 3rd order estimate of the derivative ft, only
    // the state yt.
    //
    // It is OK if yt is the same object as y0 or y1.
    static void interpolateOrder3(const Real& t0, const Vector& y0, const Vector& f0,
                                  const Real& t1, const Vector& y1, const Vector& f1,
                                  const Real& t, Vector& yt) {
        assert(t0 < t1);
        assert(t0 <= t && t <= t1);
        assert(f0.size()==y0.size() && y1.size()==y0.size() && f1.size()==y0.size());

        const Real h = t1-t0, d=(t-t0)/h;
        const Real cy1 = d*d*(3-2*d), cy0 = 1-cy1;
        const Real hdd1 = h*d*(d-1), cf1=hdd1*d, cf0=cf1-hdd1;

        yt = cy0*y0 + cy1*y1 + cf0*f0 + cf1*f1; // + O(h^4)
    }

    // We have bracketed a zero crossing for some function f(t)
    // between (tLow,fLow) and (tHigh,fHigh), not including
    // the end points. That means tHigh > tLow, and
    // sign(fLow) != sign(fHigh) (sign(x) returns -1, 0, 1).
    // We want to estimate time tRoot with
    // tLow<tRoot<tHigh such that f(tRoot) is zero.
    // For a nicely behaved continuous function this is just
    // the secant method:
    //      x = fHigh/(fHigh-fLow)  (0 <= x <= 1)
    //      tRoot = tHigh - x*(tHigh-tLow)
    // However, if the function appears to be discontinuous we'll
    // simply bisect the interval (x==0.5 above). We decide it is
    // discontinuous if either end point is exactly zero, which
    // would occur with a boolean function, for example. Also, if
    // the time interval is already at or below the smallest allowable
    // localization window, we'll bisect.
    //
    // One further twist, taken from CVODES, is to allow the caller
    // to provide a bias (> 0) which will bias our returned tRoot
    // further into the lower (bias<1) or upper (bias>1) half-interval.
    // If bias==1 this is the pure secant method. Bias has no effect
    // on functions deemed to be discontinuous; we'll always bisect
    // the interval in that case.
    //
    // Finally, we won't return tRoot very close to either end of
    // the interval. Instead we define a buffer zone at either end
    // of width 10% of the interval, and push tRoot away from the
    // edges if it gets any closer. And in any case we'll require
    // at least 1/2 minWindow from either end, even if 10% of the
    // interval is smaller than that.
    // 
    // Note that "minWindow" here is *not* the desired localization window
    // based on user accuracy requirements, but the (typically much smaller)
    // smallest allowable localization window based on numerical roundoff
    // considerations.

    static Real estimateRootTime(Real tLow, Real fLow, Real tHigh, Real fHigh,
                                 Real bias, Real minWindow)
    {
        assert(tLow < tHigh);
        assert(sign(fLow) != sign(fHigh));
        assert(bias > 0);
        assert(minWindow > 0);

        const Real h = tHigh-tLow;

        if (fLow==0 || fHigh==0 || h <= minWindow) {
            // bisecting
            return tLow + h/2;
        }

        // Use secant method.

        const Real x = fHigh/(fHigh-bias*fLow);
        Real tRoot = tHigh - x*h;

        // If tRoot is too close to either end point we'll assume bad behavior
        // and guess a value 10% of the interval away from the end.

        const Real BufferZone = std::max(0.1*h, minWindow/2);
        tRoot = std::max(tRoot, tLow  + BufferZone);
        tRoot = std::min(tRoot, tHigh - BufferZone);

        return tRoot;
    }

    // Here we look at a pair of event trigger function values across an interval
    // and decide if there are any events triggering. If so we return that event
    // as an "event candidate". Optionally, pass in the current list of event candidates
    // and we'll only look at those (that is, the list can only be narrowed). Don't
    // use the same array for the current list and new list.
    //
    // For purposes of this method, events are specified by their indices in the array
    // of trigger functions, NOT by their event IDs.
    void findEventCandidates(int nEvents, 
                             const Array_<SystemEventTriggerIndex>* viableCandidates,
                             const Array_<Event::Trigger>*          viableCandidateTransitions,
                             Real tLow,  const Vector& eLow, 
                             Real tHigh, const Vector& eHigh,
                             Real bias,  Real          minWindow,
                             Array_<SystemEventTriggerIndex>&   candidates,
                             Array_<Real>&                      timeEstimates,
                             Array_<Event::Trigger>&            transitions,
                             Real&                              earliestTimeEst, 
                             Real&                              narrowestWindow) const
    {
        int nCandidates;
        if (viableCandidates) {
            nCandidates = (int)viableCandidates->size();
            assert(viableCandidateTransitions && (int)viableCandidateTransitions->size()==nCandidates);
        } else {
            assert(!viableCandidateTransitions);
            nCandidates = nEvents;
        }

        candidates.clear();
        timeEstimates.clear();
        transitions.clear();
        earliestTimeEst = narrowestWindow = Infinity;
        for (int i=0; i<nCandidates; ++i) {
            const SystemEventTriggerIndex e = viableCandidates ? (*viableCandidates)[i] 
                                                                 : SystemEventTriggerIndex(i);
            Event::Trigger transitionSeen =
                Event::maskTransition(
                    Event::classifyTransition(sign(eLow[e]), sign(eHigh[e])),
                    eventTriggerInfo[e].calcTransitionMask());

            if (transitionSeen != Event::NoEventTrigger) {
                // Replace the transition we just saw with the appropriate one for
                // reporting purposes. For example, if the event trigger only wants
                // negative-to-positive transitions but we just saw negative-to-zero
                // we'll report that as negative-to-positive.
                transitionSeen = eventTriggerInfo[e].calcTransitionToReport(transitionSeen);
                candidates.push_back(e);
				narrowestWindow = std::max(
					std::min(narrowestWindow, 
					         accuracyInUse*timeScaleInUse*eventTriggerInfo[e].getRequiredLocalizationTimeWindow()),
				    minWindow);

                // Set estimated event trigger time for the viable candidates.
                timeEstimates.push_back(estimateRootTime(tLow, eLow[e], tHigh, eHigh[e],
                                                         bias, minWindow));
                transitions.push_back(transitionSeen);
                earliestTimeEst = std::min(earliestTimeEst, timeEstimates.back());
            }
        }
    }
    
    /// Given a list of events, specified by their indices in the list of trigger functions,
    /// convert them to the corresponding event IDs.
    void findEventIds(const Array_<SystemEventTriggerIndex>& indices, Array_<EventId>& ids) {
        for (int i = 0; i < (int)indices.size(); ++i)
            ids.push_back(eventTriggerInfo[indices[i]].getEventId());
    }

    // TODO: these utilities don't really belong here
    static Real calcWeightedRMSNorm(const Vector& values, const Vector& weights) {
        assert(weights.size() == values.size());
        if (values.size()==0) return 0;
        Real sumsq = 0;
        for (int i=0; i<values.size(); ++i) {
            const Real wv = values[i]*weights[i];
            sumsq += wv*wv;
        }
        return std::sqrt(sumsq/weights.size());
    }

    static Real calcWeightedInfinityNorm(const Vector& values, const Vector& weights) {
        assert(weights.size() == values.size());
        if (values.size()==0) return 0;
        Real maxval = 0;
        for (int i=0; i<values.size(); ++i) {
            const Real wv = std::abs(values[i]*weights[i]);
            if (wv > maxval) maxval=wv;
        }
        return maxval;
    }
    virtual const char* getMethodName() const = 0;
    virtual int getMethodMinOrder() const = 0;
    virtual int getMethodMaxOrder() const = 0;
    virtual bool methodHasErrorControl() const = 0;

protected:
    const System& getSystem() const {return sys;}

    // This is information we extract from the DynamicSystem during initialization or
    // reinitialization and save here. Only the state variable weights can be updated during
    // simulation (since they are configuration dependent). However, the weights must remain
    // constant across a time step, and usually across many time steps.
    bool getDynamicSystemHasTimeAdvancedEvents()      const {return systemHasTimeAdvancedEvents;}
    Real getDynamicSystemTimescale()                  const {return timeScaleInUse;}
    const Array_<System::EventTriggerInfo>&
        getDynamicSystemEventTriggerInfo()            const {return eventTriggerInfo;}
    const Vector& getDynamicSystemOneOverTolerances() const {return constraintWeightsInUse;}
    const Vector& getDynamicSystemWeights()           const {return stateWeightsInUse;}

    const State& getInterpolatedState() const {return interpolatedState;}
    State&       updInterpolatedState()       {return interpolatedState;}

    State& updAdvancedState() {return advancedState;}

    void setAdvancedState(const Real& t, const Vector& y) {
        advancedState.updY() = y;
        advancedState.updTime() = t;
    }

    void setTriggeredEvents(Real tlo, Real thi,
                            const Array_<EventId>&  eventIds,
                            const Array_<Real>& estEventTimes,
                            const Array_<Event::Trigger>& transitionsSeen)
    {
        assert(tPrev <= tlo && tlo < thi && thi <= advancedState.getTime());
        tLow = tlo;
        tHigh = thi;

        const int n = eventIds.size();
        assert(n > 0 && estEventTimes.size()==n && transitionsSeen.size()==n);
        triggeredEvents.resize(n); estimatedEventTimes.resize(n); eventTransitionsSeen.resize(n);
        Array_<int> eventOrder; // will be a permutation of 0:n-1
        calcEventOrder(eventIds, estEventTimes, eventOrder);
        for (int i=0; i<(int)eventOrder.size(); ++i) {
            const int ipos = eventOrder[i];
            triggeredEvents[i] = eventIds[ipos];

            assert(tlo < estEventTimes[ipos] && estEventTimes[ipos] <= thi);
            estimatedEventTimes[i] = estEventTimes[ipos];

            // TODO: assert that the transition is one of the allowed ones for this event
            eventTransitionsSeen[i] = transitionsSeen[ipos];
        }
    }

    Real getEventWindowLow()  const {return tLow;}
    Real getEventWindowHigh() const {return tHigh;}

    const Array_<EventId>&  getTriggeredEvents()  const {return triggeredEvents;}
    const Array_<Real>& getEstimatedEventTimes()  const {return estimatedEventTimes;}
    const Array_<Event::Trigger>&
                       getEventTransitionsSeen() const {return eventTransitionsSeen;}

    // This determines which state will be returned by getState().
    void setUseInterpolatedState(bool shouldUse) {
        useInterpolatedState = shouldUse;
    }

    void setStepCommunicationStatus(StepCommunicationStatus scs) {
        stepCommunicationStatus = scs;
    }

    StepCommunicationStatus getStepCommunicationStatus() const {
        return stepCommunicationStatus;
    }

    const Real&   getPreviousTime()   const {return tPrev;}
    const Vector& getPreviousY()      const {return yPrev;}
    const Vector& getPreviousYDot()   const {return ydotPrev;}
    const Vector& getPreviousEventTriggers() const {return triggersPrev;}
    Real&   updPreviousTime()   {return tPrev;}
    Vector& updPreviousY()      {return yPrev;}
    Vector& updPreviousYDot()   {return ydotPrev;}
    Vector& updPreviousEventTriggers() {return triggersPrev;}
    Array_<System::EventTriggerInfo>& updEventTriggerInfo() {return eventTriggerInfo;}
    Vector& updConstraintWeightsInUse() {return constraintWeightsInUse;}

    // State must already have been evaluated through Stage::Acceleration
    // or this will throw a stage violation.
    void saveStateAsPrevious(const State& s) {
        tPrev        = s.getTime();
        yPrev        = s.getY();
        ydotPrev     = s.getYDot();
        triggersPrev = s.getEventTriggers();
    }

    // collect user requests
    Real userInitStepSize, userMinStepSize, userMaxStepSize;
    Real userAccuracy; // use for relTol, absTol, constraintTol
    Real userRelTol, userAbsTol, userConsTol; // for fussy people
    Real userFinalTime; // never go past this
    long userInternalStepLimit; // that is, in a single call to step(); 0=no limit

    // three-state booleans
    int  userReturnEveryInternalStep;   // -1 (not supplied), 0(false), 1(true)
    int  userProjectEveryStep;          //      "
    int  userAllowInterpolation;        //      "
    int  userProjectInterpolatedStates; //      "

    // Mark all user-supplied options "not supplied by user".
    void initializeUserStuff() {
        userInitStepSize = userMinStepSize = userMaxStepSize = -1.;
        userAccuracy = userRelTol = userAbsTol = userConsTol = -1.;
        userFinalTime = -1.;
        userInternalStepLimit = -1;

        // booleans
        userReturnEveryInternalStep = userProjectEveryStep = userAllowInterpolation
            = userProjectInterpolatedStates = -1;

        accuracyInUse = NaN;
        consTol  = NaN;
        relTol   = NaN;
        absTol   = NaN;
    }

    // Required accuracy and constraint tolerance.
    // These are obtained during initialization from the user requests above, and
    // given default values if the user was agnostic.

    Real consTol;   // fraction of the constraint unit tolerance to be applied to each constraint

    // These are here just to accommodate integration methods which use these concepts.
    Real relTol;    // relative tolerance to be used for state variable (default==accuracy)
    Real absTol;    // absolute tolerance to be used for state variables near zero 

    void setAccuracyAndTolerancesFromUserRequests() {
        accuracyInUse = (userAccuracy != -1. ? userAccuracy : 1e-3);
        consTol       = (userConsTol  != -1. ? userConsTol  : 0.1*accuracyInUse); 
        relTol        = (userRelTol   != -1. ? userRelTol   : accuracyInUse); 
        absTol        = (userAbsTol   != -1. ? userAbsTol   : 0.1*accuracyInUse); 
    }
    
    // If this is set to true, the next call to stepTo() will return immediately
    // with result code StartOfContinuousInterval.
    bool startOfContinuousInterval;
    
    // The reason the simulation ended.  If isSimulationOver() returns false,
    // the value of this field is meaningless.
    Integrator::TerminationReason terminationReason;

    // Realize the supplied state through Stage::Acceleration
    // and bump statistics appropriately. Throws
    // an exception if it fails, with failure statistics bumped.
    void realizeStateDerivatives(const State& s) const {
        if (s.getSystemStage() < Stage::Acceleration) {
            ++statsRealizations; ++statsRealizationFailures;
            getSystem().realize(s, Stage::Acceleration);
            --statsRealizationFailures;
        }
    }

    // Project the supplied state onto the constraint manifold and
    // remove the corresponding errors from errEst. Throws an exception
    // if it fails. Updates stats.
    void projectStateAndErrorEstimate(State& s, Vector& errEst, System::ProjectOptions opts=System::ProjectOptions::All) {
        ++statsProjections; ++statsProjectionFailures;
        getSystem().project(s,
            consTol,getDynamicSystemWeights(),getDynamicSystemOneOverTolerances(),errEst,opts);
        --statsProjectionFailures;
    }

    // Set the advanced state and then evaluate state derivatives. Throws an
    // exception if it fails. Updates stats.
    void setAdvancedStateAndRealizeDerivatives(const Real& t, const Vector& y) 
    {
        setAdvancedState(t,y);
        realizeStateDerivatives(getAdvancedState());
    }

    // Set the advanced state and then evaluate constraint errors. Throws an
    // exception if it fails. Never counts as a realization because
    // we only need to realize kinematics.
    void setAdvancedStateAndRealizeKinematics(const Real& t, const Vector& y)
    {
        setAdvancedState(t,y);
        getSystem().realize(getAdvancedState(), Stage::Velocity);
    }


    void resetIntegratorStatistics() {
        statsRealizations = 0;
        statsProjections = 0;
        statsRealizationFailures = 0;
        statsProjectionFailures = 0;
    }

    long getNumRealizations() const {
        return statsRealizations;
    } 
    long getNumProjections() const {
        return statsProjections;
    } 
    long getNumRealizationFailures() const {
        return statsRealizationFailures;
    } 
    long getNumProjectionFailures() const {
        return statsProjectionFailures;
    } 

private:
    class EventSorter {
    public:
        EventSorter() : ordinal(-1), id(-1), estTime(NaN) { }
        EventSorter(int ord, int eventId, Real estEventTime) 
          : ordinal(ord), id(eventId), estTime(estEventTime) { }
        bool operator<(const EventSorter& r) const {
            if (estTime < r.estTime) return true;
            if (estTime > r.estTime) return false;
            return id < r.id;
        }

        int  ordinal; // original position in the eventIds array
        int  id;
        Real estTime;
    };

    void calcEventOrder(const Array_<EventId>&  eventIds,
                        const Array_<Real>& estEventTimes,
                        Array_<int>&        eventOrder)
    {
        const size_t n = eventIds.size();
        assert(estEventTimes.size()==n);
        eventOrder.resize(n);
        if (n==0) 
            return;

        if (n==1) {
            eventOrder[0] = 0;
            return;
        }

        // otherwise sort
        Array_<EventSorter> events(n);
        for (unsigned i=0; i<n; ++i)
            events[i] = EventSorter(i, eventIds[i], estEventTimes[i]);
        std::sort(events.begin(), events.end());
        for (size_t i=0; i<n; ++i) 
            eventOrder[i] = events[i].ordinal;
    }

private:
    Integrator* myHandle;
    friend class Integrator;

    const System& sys;


protected:
    // realization and projection stats are shared by all integrators;
    // others are left to the individual integration methods
    mutable long statsProjectionFailures;
    mutable long statsProjections;
    mutable long statsRealizations;
    mutable long statsRealizationFailures;
private:

        // DYNAMIC SYSTEM INFORMATION
        // Information extracted from the DynamicSystem describing properties we need
        // to know about BEFORE taking a time step. These properties are always frozen
        // across an integration step, but may be updated as discrete updates during
        // time stepping.

    // Some DynamicSystems need to get control whenever time is advanced
    // successfully (and irreversibly) so that they can do discrete updates.
    // This is Stage::Model information.
    bool systemHasTimeAdvancedEvents;


    // Event trigger functions specify which sign transitions should be
    // monitored: rising, falling, to/from zero, or combinations of those.
    // They also provide a localization window width.
    // This is Stage::Instance information.

    Array_<System::EventTriggerInfo> eventTriggerInfo;

    // A unitless fraction.
    Real accuracyInUse;

    // The time scale is an estimate as to what is the smallest length of time
    // on which we can expect any significant changes to occur. This can be
    // used to estimate an initial step size, and it also establishes the units
    // in which event trigger localization windows are defined (that is, the 
    // window is given as some fraction of the time scale).
    // This is Stage::Instance information.
    Real timeScaleInUse;

    // These are the constraint unit tolerances (actually 1/tol), given in 
    // physical units. The actual tolerances enforced during integration
    // will be these numbers scaled by the user-requested accuracy.
    // This is Stage::Instance information.
    Vector constraintWeightsInUse;

    // These are the weights for each continuous state variable. These are
    // configuration-dependent (i.e., Stage::Position) but never change during
    // an integration interval so are considered to be constants. The integrator
    // is expected to ask the DynamicSystem to update them from time to time.
    Vector stateWeightsInUse;


        // INTEGRATOR INTERNAL STATE 
        // Persists between stepTo() calls.
        // A concrete integrator is free to have more state.

    StepCommunicationStatus stepCommunicationStatus;

    // We save both the step size we have to take next and the one
    // we wish we could take. Once we're done with whatever's causing
    // us to come up short (like event isolation), we may be able to
    // jump right back up to the ideal one.
    Real    nextStepSizeToTry;  // This is what we'll actually try first.
    Real    idealNextStepSize;  // But if accuracy were the only concern
                                //   we would try this (>=nextStepSizeToTry).

    // This is how far the integrator has advanced the system
    // trajectory. We might allow an interpolated state a
    // little earlier than this, but otherwise there is no
    // going back from this point.
    State   advancedState;

    // When stepCommunicationStatus indicates that an event has
    // triggered, the following arrays are valid. Events are sorted
    // by estimated occurrence time within the window; identical-time
    // events are sorted in ascending order of event Id.

    // These are the events that the integrator has algorithmically
    // determined are now triggered. You may not get the same result
    // simply comparing the trigger function values at tLow and tHigh.
    Array_<EventId>  triggeredEvents;

    // These are the estimated times corresponding to the triggeredEvents.
    // They are in ascending order although there may be duplicates.
    Array_<Real> estimatedEventTimes;
    
    // Which transition was seen for each triggered event (this is 
    // only a single transition, not an OR-ed together set). This is
    // the integrator's algorithmic determination of the transition to
    // be reported -- you might not get the same answer just looking
    // at the event trigger functions at tLow and tHigh.
    Array_<Event::Trigger> eventTransitionsSeen;

    // When we have successfully localized a triggering event into
    // the time interval (tLow,tHigh] we record the bounds here.
    Real    tLow, tHigh;

    State   interpolatedState;    // might be unused
    bool    useInterpolatedState;

    // Use these to record the continuous part of the previous
    // accepted state. We use these in combination with the 
    // continuous contents of advancedState to fill in the
    // interpolated state.
    Real    tPrev;
    Vector  yPrev;

    // Save continuous derivatives and event function values
    // from previous accepted state also so we don't have to
    // recalculate for restarts.
    Vector  ydotPrev;
    Vector  triggersPrev;

    // We'll leave the various arrays above sized as they are and full
    // of garbage. They'll be resized when first assigned to something
    // meaningful.
	void invalidateIntegratorInternalState() {
		stepCommunicationStatus = InvalidStepCommunicationStatus;
		nextStepSizeToTry       = NaN;
		idealNextStepSize       = NaN;
        tLow = tHigh            = NaN;
        useInterpolatedState    = false;
        tPrev                   = NaN;
	}

    // suppress
    IntegratorRep(const IntegratorRep&);
    IntegratorRep& operator=(const IntegratorRep&);
};

} // namespace SimTK

#endif // SimTK_SIMMATH_INTEGRATOR_REP_H_


