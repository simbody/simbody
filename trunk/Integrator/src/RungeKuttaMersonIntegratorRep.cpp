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
 * This is the private (library side) implementation of the 
 * RungeKuttaMersonIntegratorRep class which is a concrete class
 * implementing the abstract IntegratorRep.
 */

#include "SimTKcommon.h"
#include "simmath/Integrator.h"

#include "RungeKuttaMersonIntegratorRep.h"

#include <exception>
#include <limits>

using namespace SimTK;

void RungeKuttaMersonIntegratorRep::methodInitialize(const State& initState) {
    reconstructForNewModel();

    initializeIntegrationParameters();

    haveTakenAStep = false;
    initialized = true;
}

Integrator::SuccessfulStepStatus 
RungeKuttaMersonIntegratorRep::stepTo(Real reportTime, 
                                      Real scheduledEventTime)
{
  try
  { assert(initialized);
    assert(reportTime >= getState().getTime());
    assert(scheduledEventTime >= getState().getTime());
    
    // If this is the start of a continuous interval, return immediately so
    // the current state will be seen as part of the trajectory.
    if (startOfContinuousInterval) {
        startOfContinuousInterval = false;
        return Integrator::StartOfContinuousInterval;
    }
    
    // tMax is the time beyond which we cannot advance internally.
    Real tMax = std::min(scheduledEventTime, finalTime);
    if (!allowInterpolation) tMax = std::min(tMax, reportTime);

    // tReturn is the next scheduled time to return control to the user.
    // Events may cause an earlier return. If interpolation is not allowed
    // then tReturn and tMax are the same thing.
    const Real tReturn = std::min(reportTime, tMax);

    // Count the number of internal steps taken during this call to stepTo().
    // Normally this is limited by maxNumInternalSteps, meaning even if
    // we don't reach reportTime we'll return after that many steps.
    int internalStepsTaken = 0;

    // hNeededForAccuracy tracks the step size we would like to have taken,
    // based only on accuracy considerations and the maximum limit, if any,
    // on individual steps. The actual step we take may be limited by the
    // occurrence of some unfortunate event such as a time limit being
    // reached or event triggering. In those cases we want to resume with
    // hNeededForAccuracy after the crisis has passed.
    Real hNeededForAccuracy = std::min(predictedNextStep, maxStepSize);

    bool wasLastStep=false;
    for(;;) { // MAIN STEPPING LOOP

        // At this point the system's advancedState is the one to which
        // we last advanced successfully. It has been realized through the
        // Acceleration Stage. Behavior depends on the current state of
        // the step communication state machine.

        switch (getStepCommunicationStatus()) {
          case FinalTimeHasBeenReturned:
            assert(!"can't call after final time reported"); //TODO: throw
            break;

          case StepHasBeenReturnedWithEvent:
            setUseInterpolatedState(false);
            setStepCommunicationStatus(StepHasBeenReturnedNoEvent);
            break;

          case StepHasBeenReturnedNoEvent:
            if (getAdvancedTime() >= finalTime) {
                setUseInterpolatedState(false);
                setStepCommunicationStatus(FinalTimeHasBeenReturned);
                terminationReason = Integrator::ReachedFinalTime;
                return Integrator::EndOfSimulation;
            }
            // need to advance
            break;

          case CompletedInternalStepWithEvent: {
            // Time has been advanced, but the step ended by triggering an event
            // at tAdvanced (==tHigh). We will return the last-good state at tLow<tHigh,
            // but first we need to dispatch any pending reports that are supposed
            // to occur before tLow in which case we interpolate back and
            // return control before moving on.
            if (reportTime <= getEventWindowLow()) {
                // Report time reached: take a brief time-out to make an interpolated report.
                // After that we'll come right back here with the user having advanced
                // the reportTime (hopefully).
                if (reportTime < getAdvancedTime()) {
                    createInterpolatedState(reportTime);
                    setUseInterpolatedState(true);
                }
                else
                    setUseInterpolatedState(false);
                // No change to step communication status -- this doesn't count as
                // reporting the current step since it is earlier than advancedTime.
                return Integrator::ReachedReportTime;
            }

            // The advanced state is at tHigh, but needs event handling.

            // Report the last "pre-event" state.
            createInterpolatedState(getEventWindowLow()); 
            setUseInterpolatedState(true);
            setStepCommunicationStatus(StepHasBeenReturnedWithEvent);
            return Integrator::ReachedEventTrigger;
          }

          case CompletedInternalStepNoEvent: {
            // Time has been advanced. If there is a report due before tAdvanced,
            // then we interpolate back and return control before moving on.
            if (reportTime <= getAdvancedTime()) {
                // Report time reached: take a brief time-out to make an interpolated report.
                // After that we'll come right back here with the user having advanced
                // the reportTime (hopefully).
                if (reportTime < getAdvancedTime()) {
                    createInterpolatedState(reportTime);
                    setUseInterpolatedState(true);
                }
                else
                    setUseInterpolatedState(false);
                // No change to step communication status -- this doesn't count as
                // reporting the current step since it is earlier than advancedTime.
                return Integrator::ReachedReportTime;
            }

            Integrator::SuccessfulStepStatus reportReason = Integrator::InvalidSuccessfulStepStatus;
            setUseInterpolatedState(false);

            if (getAdvancedTime() >= scheduledEventTime) {
                reportReason = Integrator::ReachedScheduledEvent;
            } else if (returnEveryInternalStep) {
                reportReason = Integrator::TimeHasAdvanced;
            } else if (getAdvancedTime() >= reportTime) {
                reportReason = Integrator::ReachedReportTime;
            } else if (getAdvancedTime() >= finalTime) {
                reportReason = Integrator::ReachedReportTime;
            } else if (maxNumInternalSteps > 0 && internalStepsTaken >= maxNumInternalSteps) {
                // Last-ditch excuse: too much work.
                reportReason = Integrator::ReachedStepLimit; // too much work
            }

            if (reportReason != Integrator::InvalidSuccessfulStepStatus) {
                setStepCommunicationStatus(StepHasBeenReturnedNoEvent);
                return reportReason;
            }

            // no return required; need to advance time
            break;
          }

          default:
            assert(!"unrecognized stepCommunicationStatus");
        }

        // If a report or event is requested at the current time, return immediately.
        if (getState().getTime() == reportTime)
            return Integrator::ReachedReportTime;
        if (getState().getTime() == scheduledEventTime)
            return Integrator::ReachedScheduledEvent;

            // NO EXCUSE LEFT: MUST ADVANCE TIME
    
        // We will now try to take the biggest step
        // we can consistent with accuracy requirements and tLimit,
        // and will continue to take steps until we pass a report time or
        // encounter a limit or event trigger.

        // Keep track of whether we had to use an excessively small step size
        // in order to stop at tMax. In that case we don't want to use the
        // shrunken step to predict the next one.
        bool trialStepSizeWasLimitedByTMax = false;

        wasLastStep = false;

        // First guess at the next step size. Note that we carry around the 
        // *time* we want at the end of the step rather than the step size itself
        // because we want to be able to depend on specified times being reached
        // exactly, rather than being off by a bit or two.
        Real trialTime = getAdvancedTime() + hNeededForAccuracy; 

        // If this step would fall very near the next return time (which might
        // be due to an interpolation time needed for reporting or a hard time
        // limit), we'll stretch it or shrink it by up to 5% to get there. 
        if (std::abs(trialTime - tReturn) <= 0.05*(trialTime - getAdvancedTime()))
            trialTime = tReturn;
        
        // Under no circumstances can we pass tMax. If that causes us to shrink
        // substantially we'll make note of that and not attempt to 
        // grow the step size from here if we succeed with this tiny thing.
        // Also, if trialTime falls just a little shy of tMax we'll stretch it.
        // Note that the above may have stretched it to the next report time;
        // if that now brings us to within a sliver of tMax we'll stretch it
        // a little more to avoid having to take that sliver step later.
        if ((getAdvancedTime() + 1.05*(trialTime - getAdvancedTime())) >= tMax) {
            const Real originalStepSize = trialTime - getAdvancedTime();
            const Real adjustedStepSize = tMax      - getAdvancedTime();
            trialStepSizeWasLimitedByTMax = adjustedStepSize < 0.8*originalStepSize;
            trialTime = tMax;
        }

        // At this point we know we are going to have to advance the state, meaning
        // that the current state will become the previous one. We need to remember
        // the continuous pieces of the current state for restarts and interpolation.
        // These will be updated below only when we make irreversible progress. Otherwise
        // we'll use them to put things back the way we found them after any failures.
        realizeStateDerivatives(getAdvancedState());
        saveStateAsPrevious(getAdvancedState());

        // Now we're going to attempt to take as big a step as we can
        // take. We'll keep shrinking trialH until something
        // works or we die.
        Real lastTrialStepAttempted = NaN;
        AttemptedStepResult stepAttemptResult;
        for (;;) {  // STEP ATTEMPT LOOP
            Real scalarError, tLow, tHigh;
            stepAttemptResult = 
                attemptAStep(getPreviousTime(), getPreviousY(), getPreviousYDot(),
                             trialTime, reportTime, scalarError, tLow, tHigh);

            // We'll reset this if we have to stop early because of an event.
            lastTrialStepAttempted = trialTime - getPreviousTime();

            switch (stepAttemptResult) {
              case SuccessWithoutEvent:
                adjustStepSize(lastTrialStepAttempted, trialStepSizeWasLimitedByTMax, 
                               scalarError, hNeededForAccuracy);
                break;


              case SuccessWithEvent:
                lastTrialStepAttempted = getAdvancedTime() - getPreviousTime();
                // don't change the step size
                break;

              case ErrorTestFailure:
                ++statsErrorTestFailures;
                if (lastTrialStepAttempted <= minStepSize ) {
                    // Restore the state and bail out. No need to realize it.
                    setAdvancedState(getPreviousTime(), getPreviousY());
                    SimTK_THROW3(Integrator::StepSizeTooSmall, getAdvancedTime(),
                        lastTrialStepAttempted, minStepSize);
                }
                //Error estimate is available.
                adjustStepSize(lastTrialStepAttempted, trialStepSizeWasLimitedByTMax, scalarError,
                               hNeededForAccuracy);
                break;

              case EvaluationFailure:
                if (lastTrialStepAttempted <= minStepSize ) {
                    // Restore the state and bail out. No need to realize it.
                    setAdvancedState(getPreviousTime(), getPreviousY());
                    SimTK_THROW3(Integrator::StepSizeTooSmall, getAdvancedTime(),
                        lastTrialStepAttempted, minStepSize);
                }
                // No error estimate available.
                adjustStepSize(lastTrialStepAttempted, trialStepSizeWasLimitedByTMax, -1,
                               hNeededForAccuracy);
                break;

              default: assert(!"unrecognized attempAStep() return");
            }

                
            if (hNeededForAccuracy != lastTrialStepAttempted)
                ++statsStepSizeChanges;
            predictedNextStep = hNeededForAccuracy;

            if (stepAttemptResult == SuccessWithoutEvent || stepAttemptResult == SuccessWithEvent)
                break;

            // Step failure; back to previous time for another try.

            trialTime = getPreviousTime() + hNeededForAccuracy;
            trialStepSizeWasLimitedByTMax = false;
        }

        // SUCCESS (attempted step succeeded)

        // Completed a step at h=lastTrialStepAttempted, hNeededForAccuracy is set properly
        ++internalStepsTaken;
        ++statsStepsTaken;
        previousSuccessfulStepSizeTaken = lastTrialStepAttempted;
        if (!haveTakenAStep) {
            actualInitialStepSizeTaken = previousSuccessfulStepSizeTaken;
            haveTakenAStep = true;
        }

        setStepCommunicationStatus(stepAttemptResult==SuccessWithEvent
                                     ? CompletedInternalStepWithEvent
                                     : CompletedInternalStepNoEvent);

    } // END OF MAIN STEP LOOP

    //NOTREACHED
    assert(!"can't get here!!");
  
  } catch (const std::exception& e) {
    setAdvancedState(getPreviousTime(),getPreviousY()); // restore
    SimTK_THROW2(Integrator::StepFailed, getAdvancedState().getTime(), e.what());
  } catch (...) {
    setAdvancedState(getPreviousTime(),getPreviousY()); // restore
    SimTK_THROW2(Integrator::StepFailed, getAdvancedState().getTime(), "UNKNOWN EXCEPTION");
  }

    // can't happen
    return Integrator::InvalidSuccessfulStepStatus;
}

// Advance from (t0,y0) to (t1,y1) where yd0=f(t0,y0)
// is already known. We are passed in the next reporting
// time tReport (=infinity if none) just in case we need
// to localize an event and tReport turns out to fall within
// the event localization window. In that case we adjust the
// window boundaries to get tReport at one end or the other,
// since we don't know anything about the state in the no-mans-land
// of the event window interior.
//
// If return is SuccessWithoutEvent, then scalarErrorEstimate
// will be valid and <=1. If return is ErrorTestFailure then
// scalarErrorEstimate will be valid and >1.
// If return is SuccessWithEvent, then scalarErrorEstimate will
// be <1 but refers to the full step, not the reduced step ending
// at the event occurrence. It should not be used to change the
// step size. In addition, tLow and tHigh will be the bracketing
// times for the event triggers: 
//     tPrev <= tLow < tEvent <= tHigh <= tAdvanced,
// For other returns we'll set scalarErrorEstimate=-1,
// tLow=tHigh=Infinity.
RungeKuttaMersonIntegratorRep::AttemptedStepResult 
RungeKuttaMersonIntegratorRep::attemptAStep
   (Real t0, const Vector& y0, const Vector& yd0,
    Real t1, Real tReport, Real& scalarErrorEstimate,
    Real& tLow, Real& tHigh)
{
    assert(t1 > t0);
    ++statsStepsAttempted;

    scalarErrorEstimate = -1;
    tLow = tHigh = Infinity;

    // Takes the RK step with error estimate, projects back to
    // manifold and adjusts error estimate. If anything goes wrong suppresses
    // exceptions and returns false.
    const bool RKstepOK = takeAnRK4MStep(t0, y0, yd0, t1, errEst);
    if (!RKstepOK) 
        return EvaluationFailure;

    // See if we achieved the required accuracy.
    scalarErrorEstimate = IntegratorRep::calcWeightedRMSNorm(errEst,getDynamicSystemWeights()) 
                          / getAccuracyInUse();
    if (scalarErrorEstimate > 1)
        return ErrorTestFailure;

    // Passed error test; see if we can realize the new state.
    bool finalRealizeOK = true;
    try {realizeStateDerivatives(getAdvancedState());}
    catch (...) {finalRealizeOK=false;}

    if (!finalRealizeOK) {
        scalarErrorEstimate = -1; // meaningless now
        return EvaluationFailure;
    }

    // The step succeeded. Check for event triggers. If there aren't
    // any, we're done with the step. Otherwise our goal will be to
    // find the "exact" time at which the first event occurs within
    // the current interval. Then we will advance the system to that
    // time only and forget about the rest of the interval.
    //
    // Here's how this is done:
    // First, determine which events trigger across the *whole* interval.
    // At least one of those events, and no other events, are eventually
    // going to be reported as having triggered after localization, and the
    // trigger type will be the original one. That is,
    // we don't care *what* we see during localization (which could be viewed
    // as more accurate); we have already decided which events are candidates.
    // Any trigger that came and went during the current interval has missed
    // the boat permanently; that's a user error (bad design of event trigger
    // or maybe excessively loose accuracy).
    // Second, with the list of candidates and their transitions in hand we want to find the
    // time at which the first one(s) trigger, to within a localization window
    // (tlo,thi]; that is, we do not see the event triggering at tlo but we
    // do see it triggering at thi, and (thi-tlo)<=w for some time window w.
    // We then chop back the advancedState by interpolation to thi and
    // return to the caller who will deal with interpolations needed
    // for reporting, then a final interpolation at tlo.
    // Note that we should never actually return the state at thi as part of the
    // trajectory, since it is invalid until the event handler is called to fix it.
    //
    // TODO: if the next report time would occur in the interval (tlo,thi) we narrow
    // the window as needed to get either tlo=treport or thi=treport.


    const Vector& e0 = getPreviousEvents();
    const Vector& e1 = getAdvancedState().getEvents();
    assert(e0.size() == e1.size() && e0.size() == getAdvancedState().getNEvents());

    const Real MinWindow = SignificantReal*getAdvancedTime();
    std::vector<int> eventCandidates, newEventCandidates;
    std::vector<EventStatus::EventTrigger> 
        eventCandidateTransitions, newEventCandidateTransitions;
    std::vector<Real> eventTimeEstimates, newEventTimeEstimates;

    Real earliestTimeEst, narrowestWindow;

    findEventCandidates(e0.size(), 0, 0, t0, e0, t1, e1, 1., MinWindow,
                        eventCandidates, eventTimeEstimates, eventCandidateTransitions,
                        earliestTimeEst, narrowestWindow);

    if (eventCandidates.empty()) {
        // This is the normal return.
        return SuccessWithoutEvent;
    }

    tLow = t0; tHigh = t1;

    if ((tHigh-tLow) <= narrowestWindow && !(tLow < tReport && tReport < tHigh)) {
        findEventIds(eventCandidates);
        setTriggeredEvents(tLow, tHigh, eventCandidates, eventTimeEstimates, eventCandidateTransitions);
        return SuccessWithEvent;     // localized already; advanced state is right (tHigh==tAdvanced)
    }


    // We now have a list of candidates in the (tLow,tHigh] interval, but that
    // interval is too wide. We have to narrow the interval until all the
    // triggered events in the interval are happy with the interval width.
    // From above we have earliestTimeEst which is the time at which we
    // think the first event is triggering.

    Vector eLow = e0, eHigh = e1;
    Real bias = 1; // neutral

    // There is an event in (tLow,tHigh], with the eariest occurrence
    // estimated at tMid=earliestTimeEst, tLow<tMid<tHigh. 
    // Decide whether the earliest occurence is actually in the
    // (tLow,tMid] interval or (tMid,tHigh].

    // Remember which side of the interval the root estimate was in over
    // the last two iterations. -1 => (tLow,tMid], 1 => (tMid,tHigh], 0 => not
    // valid yet.
    int sideTwoItersAgo=0, sidePrevIter=0;
    do {
        if (sideTwoItersAgo != 0 && sidePrevIter != 0) {
            if (sideTwoItersAgo != sidePrevIter)
                bias = 1; // this is good; alternating intervals
            else 
                bias = sidePrevIter < 0 ? bias/2 : bias*2;
        }

        const Real tMid = (tLow < tReport && tReport < tHigh) 
                          ? tReport : earliestTimeEst;

        createInterpolatedState(tMid);

        // Failure to evaluate at the interpolated state is a disaster of some kind,
        // not something we expect to be able to recover from, so this will throw
        // an exception if it fails.
        realizeStateDerivatives(getInterpolatedState());

        const Vector& eMid = getInterpolatedState().getEvents();

        // TODO: should search in the wider interval first

        // First guess: it is in (tLow,tMid].
        findEventCandidates(e0.size(), &eventCandidates, &eventCandidateTransitions,
                            tLow, eLow, tMid, eMid, bias, MinWindow,
                            newEventCandidates, newEventTimeEstimates, newEventCandidateTransitions,
                            earliestTimeEst, narrowestWindow);

        if (!newEventCandidates.empty()) {
            sideTwoItersAgo = sidePrevIter;
            sidePrevIter = -1;

            // We guessed right -- tMid is our new tHigh, earliestTimeEst is
            // our new tMid, narrowestWindow is our localization requirement.
            tHigh = tMid; eHigh = eMid;
            eventCandidates = newEventCandidates;
            eventTimeEstimates = newEventTimeEstimates;
            // these will still be the original transitions, but only the ones
            // which are still candidates are retained
            eventCandidateTransitions = newEventCandidateTransitions;
            continue;
        }

        // Nope. It must be in the upper part of the interval (tMid,tHigh].
        findEventCandidates(e0.size(), &eventCandidates,  &eventCandidateTransitions,
                            tMid, eMid, tHigh, eHigh, bias, MinWindow,
                            newEventCandidates, newEventTimeEstimates, newEventCandidateTransitions,
                            earliestTimeEst, narrowestWindow);

        assert(!newEventCandidates.empty()); // TODO: I think this can happen if
                                             // we land exactly on a zero in eMid.

        sideTwoItersAgo = sidePrevIter;
        sidePrevIter = 1;

        tLow = tMid; eLow = eMid;
        eventCandidates = newEventCandidates;
        eventTimeEstimates = newEventTimeEstimates;
        // these will still be the original transitions, but only the ones
        // which are still candidates are retained
        eventCandidateTransitions = newEventCandidateTransitions;

    } while ((tHigh-tLow) > narrowestWindow);

    findEventIds(eventCandidates);
    setTriggeredEvents(tLow, tHigh, eventCandidates, eventTimeEstimates, eventCandidateTransitions);


    // We have to throw away all of the advancedState that occurred after
    // tHigh, by interpolating back to tHigh.
    // TODO: I think a smarter algorithm could avoid this last realize() half
    // the time since we may have just interpolated and realized this end
    // of the interval above.
    if (tHigh < getAdvancedTime()) {
        backUpAdvancedStateByInterpolation(tHigh);
        // Failure to realize here should never happen since we already
        // succeeded all the way to advancedTime. So if this fails it
        // will throw an exception that will kill the simulation.
        realizeStateDerivatives(getAdvancedState());
    }

    return SuccessWithEvent;
}


// Adjust step size after a successful or failed step. 
// Input errEst tells us what happened:
//    errEst < 0,       step failed, no error estimate available
//    0 <= errEst <= 1, step was successful at current step size
//    errEst > 1,       error test failure at current step size
// In the "no error estimate" case we just cut the step size in half.
// Otherwise, we try to choose the largest step size that would
// get the error estimate below 1. For that we must know the order
// of the integration method. Runge Kutta Merson is a 4th order
// method so 1/4 is the correct exponent.
// A number of heuristics are applied to improve practical behavior:
//    * In a single adjustment, growth is restricted to 5x,
//      shrinking to 0.1x.
//    * We won't grow at all if the step we just executed was
//      artificially shrunk due, e.g., to us being near tMax or on
//      the hunt for an event trigger. (We'll still allow shrinkage
//      in that case, although that is unlikely to be required.)
//    * We apply a safety factor (0.9) to the step size; that is,
//      we always underpredict the step size by 10% to avoid
//      thrashing.
//    * If the previous step was successful, and the outcome here
//      would produce only a small change in step size, we'll leave
//      the step size alone to provide some step size stability.
//    * If we have to shrink the step due to error test failure, we
//      will do it by at least 10%.
//    * The returned step size is limited to the range
//      [minStepSize,maxStepSize].
//
void RungeKuttaMersonIntegratorRep::adjustStepSize
   (const Real& hTaken, bool hWasArtificiallyLimited, 
    const Real& errEst, Real& hNeededForAccuracy)
{
    const Real Safety=0.9, MinShrink=0.1, MaxGrow=5;
    const Real HysteresisLow=0.9, HysteresisHigh=1.2;
    const Real Order=4; // of the error *estimate*

    if (errEst < 0) {
        // Step failed without giving us an error estimate.
        hNeededForAccuracy = std::max(0.5*hTaken, minStepSize);
        return;
    }

    // We have a valid err estimate.

    Real optimalStepScale = 
        errEst > 0 ? Safety/std::pow(errEst, 1/Order)
                   : MaxGrow;
    optimalStepScale = 
        std::min(std::max(optimalStepScale, MinShrink), MaxGrow);

    // If the previous step failed, we require the step to shrink
    // by at least HysteresisLow. If it succeeded we won't change
    // the step at all unless the calculated change is outside
    // the hysteresis interval.

    if (errEst > 1) 
        optimalStepScale = std::min(optimalStepScale, HysteresisLow);
    else if (HysteresisLow < optimalStepScale && optimalStepScale < HysteresisHigh)
        optimalStepScale = 1;

    // Prevent growth if the current step's size was artificially limited.
    if (hWasArtificiallyLimited && optimalStepScale > 1)
        optimalStepScale = 1;

    hNeededForAccuracy = optimalStepScale * hTaken;
    hNeededForAccuracy = std::min(std::max(hNeededForAccuracy, 
                                           minStepSize),
                                  maxStepSize);
}


// For a discussion of the Runge-Kutta-Merson method, see Hairer,
// Norsett & Wanner, Solving ODEs I, 2nd rev. ed. pp. 166-8. This is
// a 5-stage, first-same-as-last (FSAL) 4th order method which
// gives us an embedded 3rd order method as well, so we can extract
// a 4th-order error estimate for the 3rd-order result, which error
// estimate can then be used for step size control, since it will
// behave as h^4. We then propagate the 4th order result (whose error
// is unknown), which Hairer calls "local extrapolation".
// We call the initial state (t0,y0) and want (t0+h,y1). We are
// given the initial derivative f0=f(t0,y0), which most likely
// is left over from an evaluation at the end of the last step.
// 
// We will call the derivatives at stage f1,f2,f3,f4 but these
// are done with only two temporaries fa and fb. (What we're calling
// "f" Hairer calls "k".)


bool RungeKuttaMersonIntegratorRep::takeAnRK4MStep
   (Real t0, const Vector& y0, const Vector& f0,
    Real t1, Vector& y1err)
{
    assert(t1 > t0);

    Vector& ysave = ytmp[0]; // rename temps
    Vector& fa    = ytmp[1];
    Vector& fb    = ytmp[2];

    const Real h = t1-t0;

  try
  { setAdvancedStateAndRealizeDerivatives(t0+h/3, y0 + (h/3)*f0);
    fa = getAdvancedState().getYDot(); // fa=f1

    setAdvancedStateAndRealizeDerivatives(t0+h/3, y0 + (h/6)*(f0+fa)); // f0+f1
    fa = getAdvancedState().getYDot(); // fa=f2

    setAdvancedStateAndRealizeDerivatives(t0+h/2, y0 + (h/8)*(f0 + 3*fa)); // f0+3f2
    fb = getAdvancedState().getYDot(); // fb=f3

    // We'll need this for error estimation.
    ysave = y0 + (h/2)*(f0 - 3*fa + 4*fb); // f0-3f2+4f3
    setAdvancedStateAndRealizeDerivatives(t1, ysave);
    fa = getAdvancedState().getYDot(); // fa=f4

    // Final value. This is the 4th order accurate estimate for y1=y(t0+h)+O(h^5):
    // y1 = y0 + (h/6)*(f0 + 4 f3 + f4). Evaluate through kinematics only.
    setAdvancedStateAndRealizeKinematics(t1, y0 + (h/6)*(f0 + 4*fb + fa));
    // YErr is valid now
  } catch (...) { 
    return false; 
  }

    // This is an embedded 3rd-order estimate y1hat=y(t0+h)+O(h^4). (Apparently
    // Merson thought it was 5th order, but that is only true if
    // the function is linear w/constant coefficients; not bloody likely!)
    //     y1hat = y0 + (h/10)*(f0 + 3 f2 + 4 f3 + 2 f4)
    //
    // We don't actually have any need for y1hat, just its 4th-order
    // error estimate y1hat-y1=(1/5)(y1-ysave) (easily verified from the above).

    const Vector& y1 = getAdvancedState().getY();
    for (int i=0; i<y1.size(); ++i)
        y1err[i] = 0.2*std::abs(y1[i]-ysave[i]);

    if (!projectEveryStep) {
        const Real constraintError = 
            IntegratorRep::calcWeightedInfinityNorm(getAdvancedState().getYErr(),
                                                    getDynamicSystemOneOverTolerances());
        if (constraintError <= consTol)
            return true; // no need to project
    }

    // Project back to manifold and reduce error estimate appropriately. This
    // requires only kinematic evaluations, so doesn't count as a stage!
    try {projectStateAndErrorEstimate(updAdvancedState(), y1err);} 
    catch (...) {return false;}

    // Don't do final evaluation yet because we won't need it if
    // we fail the error test.

    return true;
}

// Interpolate the advanced state back to an earlier part of the interval,
// forgetting about the rest of the interval. This is necessary, for
// example after we have localized an event trigger to an interval tLow:tHigh
// where tHigh < tAdvanced.
void RungeKuttaMersonIntegratorRep::backUpAdvancedStateByInterpolation(Real t) {
    State& advanced = updAdvancedState();
    Vector& yinterp = ytmp[0]; // rename

    assert(getPreviousTime() <= t && t <= advanced.getTime());
    interpolateOrder3(getPreviousTime(),  getPreviousY(),  getPreviousYDot(),
                      advanced.getTime(), advanced.getY(), advanced.getYDot(),
                      t, yinterp);
    advanced.updY() = yinterp;
    advanced.updTime() = t;
    getSystem().realize(advanced, Stage::Velocity); // cheap 

    // Ignore any user request not to project interpolated states here -- this
    // is the actual advanced state which will be propagated through the
    // rest of the trajectory so we can't allow it not to satisfy the 
    // constraints!

    // But it is OK if it just *barely* satisfies the constraints so we
    // won't get carried away if the user isn't being finicky about it.
    if (!projectEveryStep) {
        const Real constraintError = 
            IntegratorRep::calcWeightedInfinityNorm(advanced.getYErr(),
                                                    getDynamicSystemOneOverTolerances());
        if (constraintError <= consTol)
            return; // no need to project
    }

    // no error estimate to project here; just pass an empty Vector
    projectStateAndErrorEstimate(advanced, Vector());
}

// Create an interpolated state at time t, which is between tPrev and tCurrent.
// If we haven't yet delivered an interpolated state in this interval, we have
// to initialize its discrete part from the advanced state.
void RungeKuttaMersonIntegratorRep::createInterpolatedState(Real t) {
    const State& advanced = getAdvancedState();
    State&       interp   = updInterpolatedState();
    interp = advanced; // pick up discrete stuff. TODO: no need to do this twice in the same interval
    interpolateOrder3(getPreviousTime(),  getPreviousY(),  getPreviousYDot(),
                      advanced.getTime(), advanced.getY(), advanced.getYDot(),
                      t, interp.updY());
    interp.updTime() = t;
    getSystem().realize(interp, Stage::Velocity); // cheap  

    if (!projectInterpolatedStates)
        return; // leave 'em in "as is" condition

    if (!projectEveryStep) {
        const Real constraintError = 
            IntegratorRep::calcWeightedInfinityNorm(interp.getYErr(),
                                                    getDynamicSystemOneOverTolerances());
        if (constraintError <= consTol)
            return; // no need to project
    }

    // no error estimate to project here; just pass an empty Vector
    projectStateAndErrorEstimate(interp, Vector());
}

const char* RungeKuttaMersonIntegratorRep::getMethodName() const {
    return "RungeKuttaMerson";
}

int RungeKuttaMersonIntegratorRep::getMethodMinOrder() const {
    return 4;
}

int RungeKuttaMersonIntegratorRep::getMethodMaxOrder() const {
    return 4;
}

bool RungeKuttaMersonIntegratorRep::methodHasErrorControl() const {
    return true;
}
