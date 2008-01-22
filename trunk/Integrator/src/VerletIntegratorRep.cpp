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
#include "VerletIntegratorRep.h"

using namespace SimTK;

VerletIntegratorRep::VerletIntegratorRep(Integrator* handle, const System& sys) : IntegratorRep(handle, sys) {
}

void VerletIntegratorRep::methodInitialize(const State& state) {
    initialized = true;
    if (userInitStepSize != -1)
        currentStepSize = userInitStepSize;
    else
        currentStepSize = 0.1*getDynamicSystemTimescale();
    if (userMinStepSize != -1)
        currentStepSize = std::max(currentStepSize, userMinStepSize);
    if (userMaxStepSize != -1)
        currentStepSize = std::min(currentStepSize, userMaxStepSize);
    resetMethodStatistics();
 }

// Create an interpolated state at time t, which is between tPrev and tCurrent.
// If we haven't yet delivered an interpolated state in this interval, we have
// to initialize its discrete part from the advanced state.
void VerletIntegratorRep::createInterpolatedState(Real t) {
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

Integrator::SuccessfulStepStatus VerletIntegratorRep::stepTo(Real reportTime, Real scheduledEventTime) {
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
          bool eventOccurred = takeOneStep(getPreviousTime(), tMax, reportTime);
          ++internalStepsTaken;
          ++statsStepsTaken;
          setStepCommunicationStatus(eventOccurred 
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

/**
 * Given initial values for all the continuous variables y=(q,u,z) and their derivatives
 * (not necessarily what's in advancedState currently), take a trial step of size 
 * h=(t1-t0), optimistically storing the result in advancedState.
 * Also estimate the absolute error in each element of y, and store them in yErrEst.
 */

bool VerletIntegratorRep::attemptAStep(Real t0, Real t1, 
                                       const Vector& q0, const Vector& qdot0, const Vector& qdotdot0, 
                                       const Vector& u0, const Vector& udot0, 
                                       const Vector& z0, const Vector& zdot0, 
                                       Vector& yErrEst)
{
    // We will catch any exceptions thrown by realize() or project() and simply treat that
    // as a failure to take a step due to the step size being too big. The idea is that the
    // caller should reduce the step size and try again, giving up only when the step size
    // goes below the allowed minimum.

  try
  { statsStepsAttempted++;
    const Real h = t1-t0;
    State& advanced = updAdvancedState();

    const int nq = advanced.getNQ();
    const int nu = advanced.getNU();
    const int nz = advanced.getNZ();

    VectorView qErrEst = yErrEst(    0, nq);
    VectorView uErrEst = yErrEst(   nq, nu);
    VectorView zErrEst = yErrEst(nq+nu, nz);

    Vector dummyErrEst; // use this when we don't want the error estimate projected
    
    // Calculate the new positions q (3rd order) and initial (1st order) 
    // estimate for the velocities u and auxiliary variables z.
    
    // These are final values (the q's will get projected, though).
    advanced.updTime() = t1;
    advanced.updQ()    = q0 + qdot0*h + 0.5*qdotdot0*h*h;

    const Vector u1_est = u0 + udot0*h;
    const Vector z1_est = z0 + zdot0*h;

    advanced.updU()    = u1_est; // these will change
    advanced.updZ()    = z1_est;

    // Here we'll project the q's to their final constrained values, and project
    // our estimated u's also. TODO: does this mess up the error estimation?

    getSystem().realize(advanced, Stage::Velocity);
    if (userProjectEveryStep == 1 || IntegratorRep::calcWeightedRMSNorm(advanced.getYErr(), getDynamicSystemOneOverTolerances()) > consTol)
        projectStateAndErrorEstimate(advanced, dummyErrEst);

    // Get new values for the derivatives.
    realizeStateDerivatives(advanced);
    
    // We're going to integrate the u's and z's with the 2nd order implicit
    // trapezoid rule: u(t+h) = u(t) + h*(f(u(t))+f(u(t+h)))/2. Unfortunately this is an
    // implicit method so we have to iterate to refine u(t+h) until this equation
    // is acceptably satisfied. We're using functional iteration here which has a
    // very limited radius of convergence.
    
    const Real tol = std::min(1e-4, 0.1*getAccuracyInUse());
    Vector usave(nu), zsave(nz); // temporaries
    bool converged = false;
    for (int i = 0; !converged && i < 10; ++i) {
        // At this point we know that the advanced state has been realized
        // through the Acceleration level, so its uDot and zDot reflect
        // the u and z state values it contains.
        usave = advanced.getU();
        zsave = advanced.getZ();

        // Get these references now -- as soon as we change u or z they
        // will be invalid.
        const Vector& udot1 = advanced.getUDot();
        const Vector& zdot1 = advanced.getZDot();
        
        // Refine u and z estimates.
        advanced.setU(u0 + 0.5*(udot0 + udot1)*h);
        advanced.setZ(z0 + 0.5*(zdot0 + zdot1)*h);

        // Project only the velocities; we're done with the q's for good.
        getSystem().realize(advanced, Stage::Velocity);
        if (userProjectEveryStep == 1 || IntegratorRep::calcWeightedRMSNorm(advanced.getYErr(), getDynamicSystemOneOverTolerances()) > consTol)
            projectStateAndErrorEstimate(advanced, dummyErrEst, System::ProjectOptions::VelocityOnly);

        // Calculate fresh derivatives UDot and ZDot.
        realizeStateDerivatives(advanced);

        // Calculate convergence as the ratio of the norm of the last delta to
        // the norm of the values prior to the last change. We're using the 2-norm
        // but this ratio would be the same if we used the RMS norm. TinyReal is
        // there to keep us out of trouble if we started at zero.
        
        const Real convergenceU = (advanced.getU()-usave).norm()/(usave.norm()+TinyReal);
        const Real convergenceZ = (advanced.getZ()-zsave).norm()/(zsave.norm()+TinyReal);
        converged = std::max(convergenceU,convergenceZ) <= tol;
    }

    // Now that we have achieved 2nd order estimates of u and z, we can use them to calculate a 3rd order
    // error estimate for q and 2nd order error estimates for u and z. Note that we have already
    // realized the state with the new values, so QDot reflects the new u's.

    qErrEst = q0 + 0.5*(qdot0+advanced.getQDot())*h // implicit trapezoid rule integral
              - advanced.getQ();                    // Verlet integral

    uErrEst = u1_est                    // explicit Euler integral
              - advanced.getU();        // implicit trapezoid rule integral
    zErrEst = z1_est - advanced.getZ(); // ditto for z's

    // TODO: because we're only projecting velocities here, we aren't going to get our
    // position errors reduced here, which is a shame.
    if (userProjectEveryStep == 1 || IntegratorRep::calcWeightedRMSNorm(advanced.getYErr(), getDynamicSystemOneOverTolerances()) > consTol)
        projectStateAndErrorEstimate(advanced, yErrEst, System::ProjectOptions::VelocityOnly);
    
    // Two different integrators were used to estimate errors: trapezoidal for Q, and
    // explicit Euler for U and Z.  This means that the U and Z errors are of a different
    // order than the Q errors.  We therefore multiply them by h so everything will be
    // of the same order.
    
    uErrEst *= h; zErrEst *= h; // everything is 3rd order in h now
    return converged;
  }
  catch (std::exception ex) {
    return false;
  }
}

/**
 * Evaluate the error that occurred in the step we just attempted, and select a new
 * step size accordingly.
 * 
 * @param err     the error estimate from the step that was just attempted
 * @param hWasArtificiallyLimited   tells whether the step size was artificially
 * reduced due to a scheduled event time.  If this is true, we will never attempt
 * to increase the step size.
 * @return true if the step should be accepted, false if it should be rejected and retried
 * with a smaller step size
 */

bool VerletIntegratorRep::adjustStepSize(Real err, bool hWasArtificiallyLimited) {
    const Real Safety = 0.9, MinShrink = 0.1, MaxGrow = 5;
    const Real HysteresisLow = 0.9, HysteresisHigh = 1.5;
    const Real Order = 3; // of the error *estimate*
    
    Real newStepSize = Safety*currentStepSize*std::pow(getAccuracyInUse()/err, 1.0/Order);
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

bool VerletIntegratorRep::takeOneStep(Real t0, Real tMax, Real tReport)
{
    Real t1;
    State& advanced = updAdvancedState();

    // Make sure we're realized through Stage::Acceleration. This won't
    // do anything if the state is already realized that far.
    realizeStateDerivatives(advanced);

    // Record the starting values for the continuous variables and their
    // derivatives in case we have to back up.
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
        bool converged = attemptAStep(t0, t1, q0, qdot0, qdotdot0, u0, udot0, z0, zdot0, err);
        Real rmsErr = (converged ? calcWeightedRMSNorm(err, getDynamicSystemWeights()) : Infinity);
        stepSucceeded = adjustStepSize(rmsErr, hWasArtificiallyLimited);
        if (!stepSucceeded)
            statsErrorTestFailures++;
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
void VerletIntegratorRep::backUpAdvancedStateByInterpolation(Real t) {
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

Real VerletIntegratorRep::getActualInitialStepSizeTaken() const {
    return userInitStepSize;
}

Real VerletIntegratorRep::getPreviousStepSizeTaken() const {
    return userInitStepSize;
}

Real VerletIntegratorRep::getPredictedNextStepSize() const {
    return userInitStepSize;
}

long VerletIntegratorRep::getNStepsAttempted() const {
    assert(initialized);
    return statsStepsAttempted;
}

long VerletIntegratorRep::getNStepsTaken() const {
    assert(initialized);
    return statsStepsTaken;
}

long VerletIntegratorRep::getNErrorTestFailures() const {
    return statsErrorTestFailures;
}

void VerletIntegratorRep::resetMethodStatistics() {
    statsStepsTaken = 0;
    statsStepsAttempted = 0;
    statsErrorTestFailures = 0;
}

const char* VerletIntegratorRep::getMethodName() const {
    return "Verlet";
}

int VerletIntegratorRep::getMethodMinOrder() const {
    return 2;
}

int VerletIntegratorRep::getMethodMaxOrder() const {
    return 3;
}

bool VerletIntegratorRep::methodHasErrorControl() const {
    return false;
}
