/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmath(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007-2008 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
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

#include "SimTKcommon.h"
#include "AbstractIntegratorRep.h"

using namespace SimTK;

AbstractIntegratorRep::AbstractIntegratorRep(Integrator* handle, const System& sys, int minOrder, int maxOrder, const std::string& methodName, bool hasErrorControl) :
    IntegratorRep(handle, sys), minOrder(minOrder), maxOrder(maxOrder), methodName(methodName), hasErrorControl(hasErrorControl) {
}

void AbstractIntegratorRep::methodInitialize(const State& state) {
    initialized = true;
    if (userInitStepSize != -1)
        currentStepSize = userInitStepSize;
    else
        currentStepSize = 0.1*getDynamicSystemTimescale();
    if (userMinStepSize != -1)
        currentStepSize = std::max(currentStepSize, userMinStepSize);
    if (userMaxStepSize != -1)
        currentStepSize = std::min(currentStepSize, userMaxStepSize);
    lastStepSize = currentStepSize;
    actualInitialStepSizeTaken = (hasErrorControl ? NaN : currentStepSize);
    resetMethodStatistics();
 }

// Create an interpolated state at time t, which is between tPrev and tCurrent.
// If we haven't yet delivered an interpolated state in this interval, we have
// to initialize its discrete part from the advanced state.
void AbstractIntegratorRep::createInterpolatedState(Real t) {
    const State& advanced = getAdvancedState();
    State&       interp   = updInterpolatedState();
    interp = advanced; // pick up discrete stuff.
    interpolateOrder3(getPreviousTime(),  getPreviousY(),  getPreviousYDot(),
                      advanced.getTime(), advanced.getY(), advanced.getYDot(),
                      t, interp.updY());
    interp.updTime() = t;
    getSystem().realize(interp, Stage::Velocity); // cheap  
    if (userProjectInterpolatedStates == 0) // default is yes
        return; // leave 'em in "as is" condition
    if (userProjectEveryStep != 1) {
        const Real constraintError =  IntegratorRep::calcWeightedRMSNorm(interp.getYErr(), getDynamicSystemOneOverTolerances());
        if (constraintError <= consTol)
            return; // no need to project
    }

    // no error estimate to project here; just pass an empty Vector
    projectStateAndErrorEstimate(interp, Vector());
}

Integrator::SuccessfulStepStatus AbstractIntegratorRep::stepTo(Real reportTime, Real scheduledEventTime) {
    try {
      assert(initialized);
      assert(reportTime >= getState().getTime());
      assert(scheduledEventTime >= getState().getTime());
      
      // If this is the start of a continuous interval, return immediately so
      // the current state will be seen as part of the trajectory.
      if (startOfContinuousInterval) {
          // The set of constraints or event triggers might have changed.
          getSystem().calcEventTriggerInfo(getAdvancedState(), updEventTriggerInfo());
          getSystem().calcYErrUnitTolerances(getAdvancedState(), updConstraintWeightsInUse());
          startOfContinuousInterval = false;
          return Integrator::StartOfContinuousInterval;
      }
      
      // tMax is the time beyond which we cannot advance internally.
      const Real finalTime = (userFinalTime == -1.0 ? Infinity : userFinalTime);
      const Real tMax = std::min(scheduledEventTime, finalTime);

      // tReturn is the next scheduled time to return control to the user.
      // Events may cause an earlier return.
      const Real tReturn = std::min(reportTime, tMax);

      // Count the number of internal steps taken during this call to stepTo().
      // Normally this is limited by maxNumInternalSteps, meaning even if
      // we don't reach reportTime we'll return after that many steps.
      int internalStepsTaken = 0;

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

            case StepHasBeenReturnedWithEvent:
              setUseInterpolatedState(false);
              // Fall through to the next case.

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
              } else if (userReturnEveryInternalStep == 1) {
                  reportReason = Integrator::TimeHasAdvanced;
              } else if (getAdvancedTime() >= reportTime) {
                  reportReason = Integrator::ReachedReportTime;
              } else if (getAdvancedTime() >= finalTime) {
                  reportReason = Integrator::ReachedReportTime;
              } else if (userInternalStepLimit > 0 && internalStepsTaken >= userInternalStepLimit) {
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

          // At this point we know we are going to have to advance the state, meaning
          // that the current state will become the previous one. We need to remember
          // the continuous pieces of the current state for restarts and interpolation.
          // These will be updated below only when we make irreversible progress. Otherwise
          // we'll use them to put things back the way we found them after any failures.
          realizeStateDerivatives(getAdvancedState());
          saveStateAsPrevious(getAdvancedState());
          
          // Now take a step and see whether an event occurred.
          bool eventOccurred = takeOneStep(tMax, reportTime);
          ++internalStepsTaken;
          ++statsStepsTaken;
          setStepCommunicationStatus(eventOccurred ? CompletedInternalStepWithEvent : CompletedInternalStepNoEvent);
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

bool AbstractIntegratorRep::adjustStepSize(Real err, int errOrder, bool hWasArtificiallyLimited) {
    const Real Safety = 0.9, MinShrink = 0.1, MaxGrow = 5;
    const Real HysteresisLow = 0.9, HysteresisHigh = 1.2;
    
    Real newStepSize = Safety*currentStepSize*std::pow(getAccuracyInUse()/err, 1.0/errOrder);
    if (newStepSize > currentStepSize) {
        if (hWasArtificiallyLimited || newStepSize < HysteresisHigh*currentStepSize)
            newStepSize = currentStepSize;
    }
    if (newStepSize < currentStepSize && err <= getAccuracyInUse())
        newStepSize = currentStepSize;
    newStepSize = std::min(newStepSize, MaxGrow*currentStepSize);
    newStepSize = std::max(newStepSize, MinShrink*currentStepSize);
    if (newStepSize < currentStepSize)
        newStepSize = std::min(newStepSize, HysteresisLow*currentStepSize);
    if (userMinStepSize != -1)
        newStepSize = std::max(newStepSize, userMinStepSize);
    if (userMaxStepSize != -1)
        newStepSize = std::min(newStepSize, userMaxStepSize);
    bool success = (newStepSize >= currentStepSize);
    currentStepSize = newStepSize;
    return success;
}

bool AbstractIntegratorRep::takeOneStep(Real tMax, Real tReport)
{
    Real t1;
    State& advanced = updAdvancedState();

    // Make sure we're realized through Stage::Acceleration. This won't
    // do anything if the state is already realized that far.
    realizeStateDerivatives(advanced);

    // Record the starting values for the continuous variables and their
    // derivatives in case we have to back up.
    const Real t0 = getPreviousTime();
    const Vector q0         = advanced.getQ();
    const Vector qdot0      = advanced.getQDot();
    const Vector qdotdot0   = advanced.getQDotDot();
    const Vector u0         = advanced.getU();
    const Vector udot0      = advanced.getUDot();
    const Vector z0         = advanced.getZ();
    const Vector zdot0      = advanced.getZDot();
    
    Vector err(advanced.getNY());
    bool stepSucceeded = false;
    do {
        bool hWasArtificiallyLimited = (tMax < t0+currentStepSize);
        t1 = std::min(tMax, t0+currentStepSize);
        int errOrder;
        bool converged = attemptAStep(t0, t1, q0, qdot0, qdotdot0, u0, udot0, z0, zdot0, err, errOrder);
        Real rmsErr = (converged ? calcWeightedRMSNorm(err, getDynamicSystemWeights()) : Infinity);
        if (hasErrorControl)
            stepSucceeded = adjustStepSize(rmsErr, errOrder, hWasArtificiallyLimited);
        else
            stepSucceeded = true;
        if (!stepSucceeded)
            statsErrorTestFailures++;
		else { // step succeeded
			lastStepSize = t1-t0;
			if (isNaN(actualInitialStepSizeTaken))
				actualInitialStepSizeTaken = lastStepSize;
		}
    } while (!stepSucceeded);
    
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

    realizeStateDerivatives(getAdvancedState());
    const Vector& e0 = getPreviousEventTriggers();
    const Vector& e1 = getAdvancedState().getEventTriggers();
    assert(e0.size() == e1.size() && e0.size() == getAdvancedState().getNEventTriggers());

	const Real MinWindow = std::max(SignificantReal, SignificantReal*getAdvancedTime());
    std::vector<SystemEventTriggerIndex> eventCandidates, newEventCandidates;
    std::vector<Event::Trigger> 
        eventCandidateTransitions, newEventCandidateTransitions;
    std::vector<Real> eventTimeEstimates, newEventTimeEstimates;

    Real earliestTimeEst, narrowestWindow;

    findEventCandidates(e0.size(), 0, 0, t0, e0, t1, e1, 1., MinWindow,
                        eventCandidates, eventTimeEstimates, eventCandidateTransitions,
                        earliestTimeEst, narrowestWindow);

    if (eventCandidates.empty()) {
        // This is the normal return.
        return false;
    }

    Real tLow = t0;
    Real tHigh = t1;

    if ((tHigh-tLow) <= narrowestWindow && !(tLow < tReport && tReport < tHigh)) {
        std::vector<EventId> ids;
        findEventIds(eventCandidates, ids);
        setTriggeredEvents(tLow, tHigh, ids, eventTimeEstimates, eventCandidateTransitions);
        return true;     // localized already; advanced state is right (tHigh==tAdvanced)
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

        const Vector& eMid = getInterpolatedState().getEventTriggers();

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

    std::vector<EventId> ids;
    findEventIds(eventCandidates, ids);
    setTriggeredEvents(tLow, tHigh, ids, eventTimeEstimates, eventCandidateTransitions);


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

    return true;
}

// Interpolate the advanced state back to an earlier part of the interval,
// forgetting about the rest of the interval. This is necessary, for
// example after we have localized an event trigger to an interval tLow:tHigh
// where tHigh < tAdvanced.
void AbstractIntegratorRep::backUpAdvancedStateByInterpolation(Real t) {
    State& advanced = updAdvancedState();
    Vector yinterp;

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
    if (userProjectEveryStep != 1) {
        const Real constraintError = 
            IntegratorRep::calcWeightedRMSNorm(advanced.getYErr(),
                                               getDynamicSystemOneOverTolerances());
        if (constraintError <= consTol)
            return; // no need to project
    }

    // no error estimate to project here; just pass an empty Vector
    projectStateAndErrorEstimate(advanced, Vector());
}

Real AbstractIntegratorRep::getActualInitialStepSizeTaken() const {
    return actualInitialStepSizeTaken;
}

Real AbstractIntegratorRep::getPreviousStepSizeTaken() const {
    return lastStepSize;
}

Real AbstractIntegratorRep::getPredictedNextStepSize() const {
    return currentStepSize;
}

long AbstractIntegratorRep::getNStepsAttempted() const {
    assert(initialized);
    return statsStepsAttempted;
}

long AbstractIntegratorRep::getNStepsTaken() const {
    assert(initialized);
    return statsStepsTaken;
}

long AbstractIntegratorRep::getNErrorTestFailures() const {
    return statsErrorTestFailures;
}

void AbstractIntegratorRep::resetMethodStatistics() {
    statsStepsTaken = 0;
    statsStepsAttempted = 0;
    statsErrorTestFailures = 0;
}

const char* AbstractIntegratorRep::getMethodName() const {
    return methodName.c_str();
}

int AbstractIntegratorRep::getMethodMinOrder() const {
    return minOrder;
}

int AbstractIntegratorRep::getMethodMaxOrder() const {
    return maxOrder;
}

bool AbstractIntegratorRep::methodHasErrorControl() const {
    return hasErrorControl;
}
