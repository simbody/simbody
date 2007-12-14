/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmath(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
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
#include "ExplicitEulerIntegratorRep.h"

using namespace SimTK;

ExplicitEulerIntegratorRep::ExplicitEulerIntegratorRep(Integrator* handle, const System& sys) : IntegratorRep(handle, sys) {
}

void ExplicitEulerIntegratorRep::methodInitialize(const State& state) {
    initialized = true;
    resetMethodStatistics();
 }

// Create an interpolated state at time t, which is between tPrev and tCurrent.
// If we haven't yet delivered an interpolated state in this interval, we have
// to initialize its discrete part from the advanced state.
void ExplicitEulerIntegratorRep::createInterpolatedState(Real t) {
    const State& advanced = getAdvancedState();
    State&       interp   = updInterpolatedState();
    interp = advanced; // pick up discrete stuff.
    const Real weight1 = (getAdvancedTime()-t)/(getAdvancedTime()-getPreviousTime());
    const Real weight2 = 1.0-weight1;
    interp.updY() = weight1*getPreviousY()+weight2*getAdvancedState().getY();
    interp.updTime() = t;
    getSystem().realize(interp, Stage::Velocity); // cheap  
    if (userProjectInterpolatedStates != 1)
        return; // leave 'em in "as is" condition
    if (userProjectEveryStep != 1) {
        const Real constraintError =  IntegratorRep::calcWeightedInfinityNorm(interp.getYErr(), getDynamicSystemOneOverTolerances());
        if (constraintError <= consTol)
            return; // no need to project
    }

    // no error estimate to project here; just pass an empty Vector
    projectStateAndErrorEstimate(interp, Vector());
}

Integrator::SuccessfulStepStatus ExplicitEulerIntegratorRep::stepTo(Real reportTime, Real scheduledEventTime) {
    try {
      assert(initialized);
      assert(reportTime >= getState().getTime());
      assert(scheduledEventTime >= getState().getTime());
      
      // If this is the start of a continuous interval, return immediately so
      // the current state will be seen as part of the trajectory.
      if (startOfContinuousInterval) {
          startOfContinuousInterval = false;
          return Integrator::StartOfContinuousInterval;
      }
      
      // tMax is the time beyond which we cannot advance internally.
      const Real finalTime = (userFinalTime == -1.0 ? Infinity : userFinalTime);
      Real tMax = std::min(scheduledEventTime, finalTime);

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
          Real tNext = std::min(tMax, getAdvancedTime()+userInitStepSize);
          bool eventOccurred = takeOneStep(getPreviousTime(), getPreviousY(), getPreviousYDot(), tNext, reportTime);
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

bool ExplicitEulerIntegratorRep::takeOneStep(Real t0, const Vector& y0, const Vector& f0, Real t1, Real tReport)
{
    assert(t1 > t0);
    setAdvancedStateAndRealizeDerivatives(t1, y0 + (t1-t0)*f0);
    projectStateAndErrorEstimate(updAdvancedState(), Vector());
    realizeStateDerivatives(getAdvancedState());
    
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
        return false;
    }

    Real tLow = t0;
    Real tHigh = t1;

    if ((tHigh-tLow) <= narrowestWindow && !(tLow < tReport && tReport < tHigh)) {
        findEventIds(eventCandidates);
        setTriggeredEvents(tLow, tHigh, eventCandidates, eventTimeEstimates, eventCandidateTransitions);
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
        setAdvancedStateAndRealizeDerivatives(tHigh, y0 + (tHigh-t0)*f0);
        // Failure to realize here should never happen since we already
        // succeeded all the way to advancedTime. So if this fails it
        // will throw an exception that will kill the simulation.
        realizeStateDerivatives(getAdvancedState());
    }

    return true;
}

Real ExplicitEulerIntegratorRep::getActualInitialStepSizeTaken() const {
    return userInitStepSize;
}

Real ExplicitEulerIntegratorRep::getPreviousStepSizeTaken() const {
    return userInitStepSize;
}

Real ExplicitEulerIntegratorRep::getPredictedNextStepSize() const {
    return userInitStepSize;
}

long ExplicitEulerIntegratorRep::getNStepsAttempted() const {
    assert(initialized);
    return statsStepsTaken;
}

long ExplicitEulerIntegratorRep::getNStepsTaken() const {
    assert(initialized);
    return statsStepsTaken;
}

long ExplicitEulerIntegratorRep::getNErrorTestFailures() const {
    return 0;
}

void ExplicitEulerIntegratorRep::resetMethodStatistics() {
    statsStepsTaken = 0;
}

const char* ExplicitEulerIntegratorRep::getMethodName() const {
    return "ExplicitEuler";
}

int ExplicitEulerIntegratorRep::getMethodMinOrder() const {
    return 1;
}

int ExplicitEulerIntegratorRep::getMethodMaxOrder() const {
    return 1;
}

bool ExplicitEulerIntegratorRep::methodHasErrorControl() const {
    return false;
}
