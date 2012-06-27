#ifndef SimTK_SIMMATH_INTEGRATOR_REP_H_
#define SimTK_SIMMATH_INTEGRATOR_REP_H_

/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-12 Stanford University and the Authors.        *
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
 * This is the declaration of the abstract IntegratorRep class which
 * represents the implementation of the Integrator class, and DAESystemRep
 * which implements the Integrator::DAESystem class.
 */

#include "SimTKcommon.h"
#include "simmath/Integrator.h"

#include <exception>
#include <limits>
#include <algorithm>
#include <vector>

namespace SimTK {

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
//   Simultaneous occurrences must be anticipated and handled by the caller
//   (time stepper).
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

    // The System must be successfully realized to Stage::Instance before this
    // call. At this point the integrator can query the System about the
    // problem size and allocate appropriately-sized internal data structures.
    void initialize(const State&);
    // This is for use after return from an event handler.
    void reinitialize(Stage, bool shouldTerminate);
    // We also give the concrete integration method a chance to initialize 
    // itself after the generic initialization is done.
    virtual void methodInitialize(const State&) { }
    // This is for use after return from an event handler.
    virtual void methodReinitialize(Stage stage, bool shouldTerminate) { }

    // The integrator has already been initialized. Take as many internal steps
    // as needed to get up to or past reportTime, but don't advance past the
    // next scheduled event time or the simulation final time.
    virtual Integrator::SuccessfulStepStatus 
        stepTo(Real reportTime, Real scheduledEventTime) = 0;

    bool isSimulationOver() const {
        return stepCommunicationStatus == FinalTimeHasBeenReturned;
    }
    
    // Get the reason the simulation ended. This should only be invoked if 
    // isSimulationOver() returns true.
    Integrator::TerminationReason getTerminationReason() const {
        assert (isSimulationOver());
        return terminationReason;
    }

    // This represents the Integrator's finite state machine for dealing with 
    // communication of internal steps to the caller. After initialization, the
    // computation will be in state CompletedInternalStepNoEvent, but with 
    // t_prev=t_advanced=t_initial. We will process this as though we had just 
    // taken a step so that various conditions will result in the first 
    // stepTo() call returning immediately. Note: event handling and subsequent
    // partial reinitialization of the integrator MUST NOT change the 
    // integrator's communication state.
    enum StepCommunicationStatus {
        // Time has advanced from t_prev to t_advanced; no event triggered. We
        // have not yet returned to the caller with time t_advanced, although 
        // we may have returned at interpolated report times 
        // t_prev < t_report < t_advanced.
        CompletedInternalStepNoEvent,   

        // Time has advanced from t_prev to t_advanced with an event trigger 
        // localized to the interval (t_low,t_high] where t_high==t_advanced. 
        // We have not yet returned to the caller with time t_low, although we 
        // may have returned at interpolated report times 
        // t_prev < t_report < t_low.
        CompletedInternalStepWithEvent, 

        // We have already returned control to the caller at t_advanced with 
        // the strongest stopping reason given; no more returns are allowed for
        // this step.
        StepHasBeenReturnedNoEvent,     

        // We have already returned control to the caller at t_low with 
        // "EventTriggered" as the stopping reason; no more returns are allowed
        // for this step.
        StepHasBeenReturnedWithEvent,   

        // Any call to stepTo() when the integrator is already in this state is
        // a fatal error.
        FinalTimeHasBeenReturned,       

        InvalidStepCommunicationStatus = -1
    };

    const State& getAdvancedState() const {return advancedState;}
    Real         getAdvancedTime()  const {return advancedState.getTime();}

    const State& getState() const 
    {   return useInterpolatedState ? interpolatedState : advancedState; }
    bool isStateInterpolated() const {return useInterpolatedState;}

    Real getAccuracyInUse() const {return accuracyInUse;}
    Real getConstraintToleranceInUse() const {return consTol;}
    Real getTimeScaleInUse() const {return timeScaleInUse;}

    // What was the size of the first successful step after the last initialize() call?
    virtual Real getActualInitialStepSizeTaken() const = 0;

    // What was the size of the most recent successful step?
    virtual Real getPreviousStepSizeTaken() const = 0;

    // What step size will be attempted first on the next step() call?
    virtual Real getPredictedNextStepSize() const = 0;

    virtual int getNumStepsAttempted() const = 0;
    virtual int getNumStepsTaken() const = 0; 
    virtual int getNumErrorTestFailures() const = 0;
    virtual int getNumConvergenceTestFailures() const = 0;
    virtual int getNumConvergentIterations() const = 0;
    virtual int getNumDivergentIterations() const = 0;
    virtual int getNumIterations() const = 0;

    virtual void resetMethodStatistics() {
    }

    // Cubic Hermite interpolation. See Hairer, et al. Solving ODEs I, 2nd rev.
    // ed., pg 190. Given (t0,y0,y0'),(t1,y1,y1') with y0 and y1 at least 3rd 
    // order accurate, we can obtain a 3rd order accurate interpolation yt for
    // a time t0 < t < t1 using Hermite interpolation. Let
    // f0=y0', f1=y1', h=t1-t0, d=(t-t0)/h (0<=d<=1). Then the interpolating
    // function is:
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
    static void interpolateOrder3
       (const Real& t0, const Vector& y0, const Vector& f0,
        const Real& t1, const Vector& y1, const Vector& f1,
        const Real& t, Vector& yt) {
        assert(t0 < t1);
        assert(t0 <= t && t <= t1);
        assert(f0.size()==y0.size() && y1.size()==y0.size() 
            && f1.size()==y0.size());

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

    // Here we look at a pair of event trigger function values across an 
    // interval and decide if there are any events triggering. If so we return
    // that event as an "event candidate". Optionally, pass in the current list
    // of event candidates and we'll only look at those (that is, the list can
    // only be narrowed). Don't use the same array for the current list and new
    // list.
    //
    // For purposes of this method, events are specified by their indices in 
    // the array of trigger functions, NOT by their event IDs.
    void findEventCandidates
       (int nEvents, 
        const Array_<SystemEventTriggerIndex>*  viableCandidates,
        const Array_<Event::Trigger>*           viableCandidateTransitions,
        Real    tLow,   const Vector&   eLow, 
        Real    tHigh,  const Vector&   eHigh,
        Real    bias,   Real            minWindow,
        Array_<SystemEventTriggerIndex>&        candidates,
        Array_<Real>&                           timeEstimates,
        Array_<Event::Trigger>&                 transitions,
        Real&                                   earliestTimeEst, 
        Real&                                   narrowestWindow) const
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

    // Calculate the error norm using RMS or Inf norm, and report which y
    // was dominant.
    Real calcErrorNorm(const State& s, const Vector& yErrEst, 
                       int& worstY) const {
        const int nq=s.getNQ(), nu=s.getNU(), nz=s.getNZ();
        int worstQ, worstU, worstZ;
        Real qNorm, uNorm, zNorm, maxNorm;
        if (userUseInfinityNorm == 1) {
            qNorm = calcWeightedInfNormQ(s, s.getUWeights(), yErrEst(0,nq),
                                         worstQ);
            uNorm = calcWeightedInfNorm(getPreviousUScale(), yErrEst(nq,nu),
                                        worstU);
            zNorm = calcWeightedInfNorm(getPreviousZScale(), yErrEst(nq+nu,nz),
                                        worstZ);
        } else {
            qNorm = calcWeightedRMSNormQ(s, s.getUWeights(), yErrEst(0,nq),
                                         worstQ);
            uNorm = calcWeightedRMSNorm(getPreviousUScale(), yErrEst(nq,nu),
                                        worstU);
            zNorm = calcWeightedRMSNorm(getPreviousZScale(), yErrEst(nq+nu,nz),
                                        worstZ);
        }

        // Find the largest of the three norms and report the corresponding
        // worst offender within q, u, or z.
        if (qNorm >= uNorm) {
            if (qNorm >= zNorm) 
                 {maxNorm = qNorm; worstY = worstQ;}        // q>=u && q>=z
            else {maxNorm = zNorm; worstY = nq+nu+worstZ;}  // z>q>=u
        } else { // qNorm < uNorm
            if (uNorm >= zNorm) 
                 {maxNorm = uNorm; worstY = nq + worstU;}   // u>q && u>=z
            else {maxNorm = zNorm; worstY = nq+nu+worstZ;}  // z>u>q
        }

        return maxNorm;
    }


    // Given a proposed absolute change dq to generalized coordinates q, scale 
    // it to produce fq, such that fq_i is the fraction of q_i's "unit change"
    // represented by dq_i. A suitable norm of fq can then be compared directly
    // with the relative accuracy requirement. 
    //
    // Because q's are not independent of u's (qdot=N(q)*u), the "unit change"
    // of q is related to the unit change of u. We want fq=Wq*dq, but we 
    // determine Wq from Wu via Wq = N*Wu*pinv(N). Wq is block diagonal while 
    // Wu is diagonal. State must be realized to Position stage.
    void scaleDQ(const State& state, const Vector& Wu,
                 const Vector& dq, Vector& dqw) const // in/out 
    {
        const System& system = getSystem();
        const int nq = state.getNQ();
        const int nu = state.getNU();
        assert(dq.size() == nq);
        assert(Wu.size() == nu);
        dqw.resize(nq);
        if (nq==0) return;
        Vector du(nu);
        system.multiplyByNPInv(state, dq, du);
        du.rowScaleInPlace(Wu);
        system.multiplyByN(state, du, dqw);
    }
    // Calculate |Wq*dq|_RMS=|N*Wu*pinv(N)*dq|_RMS
    Real calcWeightedRMSNormQ(const State& state, const Vector& Wu,
                              const Vector& dq, int& worstQ) const
    {
        const int nq = state.getNQ();
        Vector dqw(nq);
        scaleDQ(state, Wu, dq, dqw);
        return dqw.normRMS(&worstQ);
    }
    // Calculate |Wq*dq|_Inf=|N*Wu*pinv(N)*dq|_Inf
    Real calcWeightedInfNormQ(const State& state, const Vector& Wu,
                              const Vector& dq, int& worstQ) const
    {
        const int nq = state.getNQ();
        Vector dqw(nq);
        scaleDQ(state, Wu, dq, dqw);
        return dqw.normInf(&worstQ);
    }

    // TODO: these utilities don't really belong here
    static Real calcWeightedRMSNorm(const Vector& weights, const Vector& values, 
                                    int& worstOne) {
        return values.weightedNormRMS(weights, &worstOne);
    }

    static Real calcWeightedInfNorm(const Vector& weights, const Vector& values,
                                    int& worstOne) {
        return values.weightedNormInf(weights, &worstOne);
    }

    virtual const char* getMethodName() const = 0;
    virtual int getMethodMinOrder() const = 0;
    virtual int getMethodMaxOrder() const = 0;
    virtual bool methodHasErrorControl() const = 0;

protected:
    const System& getSystem() const {return sys;}

    // This is information we extract from the System during initialization or
    // reinitialization and save here.
    bool getDynamicSystemHasTimeAdvancedEvents()      const {return systemHasTimeAdvancedEvents;}
    Real getDynamicSystemTimescale()                  const {return timeScaleInUse;}
    const Array_<EventTriggerInfo>&
        getDynamicSystemEventTriggerInfo()            const {return eventTriggerInfo;}

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

    const Real&   getPreviousTime()          const {return tPrev;}

    // q,u,z are views into y.
    const Vector& getPreviousY()            const {return yPrev;}
    const Vector& getPreviousQ()            const {return qPrev;}
    const Vector& getPreviousU()            const {return uPrev;}
    const Vector& getPreviousZ()            const {return zPrev;}

    const Vector& getPreviousUScale()       const {return uScalePrev;}
    const Vector& getPreviousZScale()       const {return zScalePrev;}

    // qdot,udot,zdot are views into ydot.
    const Vector& getPreviousYDot()         const {return ydotPrev;}
    const Vector& getPreviousQDot()         const {return qdotPrev;}
    const Vector& getPreviousUDot()         const {return udotPrev;}
    const Vector& getPreviousZDot()         const {return zdotPrev;}

    const Vector& getPreviousQDotDot()      const {return qdotdotPrev;}

    const Vector& getPreviousEventTriggers() const {return triggersPrev;}

    Array_<EventTriggerInfo>& updEventTriggerInfo() {return eventTriggerInfo;}

    // Given an array of state variables v (either u or z) and corresponding
    // weights w (1/absolute scale), return relative scale min(wi,1/|vi|) for
    // each variable i. That is, we choose the current value of vi as its
    // scale when it is large enough, otherwise use absolute scale.
    static void calcRelativeScaling(const Vector& v, const Vector& w, 
                                    Vector& vScale)
    {
        const int nv = v.size();
        assert(w.size() == nv);
        vScale.resize(nv);
        for (int i=0; i<nv; ++i) {
            const Real vi = std::abs(v[i]);
            const Real wi = w[i];
            vScale[i] = vi*wi > 1 ? 1/vi : wi;
        }
    }


    // State must already have been evaluated through Stage::Acceleration
    // or this will throw a stage violation. We calculate the scaling for
    // u and z here which may include relative scaling based on their current
    // values. This scaling is frozen during a step attempt.
    void saveStateAsPrevious(const State& s) {
        const int nq = s.getNQ(), nu = s.getNU(), nz = s.getNZ();

        tPrev        = s.getTime();

        yPrev        = s.getY();
        qPrev.viewAssign(yPrev(0,     nq));
        uPrev.viewAssign(yPrev(nq,    nu));
        zPrev.viewAssign(yPrev(nq+nu, nz));

        calcRelativeScaling(s.getU(), s.getUWeights(), uScalePrev); 
        calcRelativeScaling(s.getZ(), s.getZWeights(), zScalePrev);

        ydotPrev     = s.getYDot();
        qdotPrev.viewAssign(ydotPrev(0,     nq));
        udotPrev.viewAssign(ydotPrev(nq,    nu));
        zdotPrev.viewAssign(ydotPrev(nq+nu, nz));

        qdotdotPrev  = s.getQDotDot();
        triggersPrev = s.getEventTriggers();
    }

    // collect user requests
    Real userInitStepSize, userMinStepSize, userMaxStepSize;
    Real userAccuracy; // also use for constraintTol
    Real userConsTol; // for fussy people
    Real userFinalTime; // never go past this
    int  userInternalStepLimit; // that is, in a single call to step(); 0=no limit

    // three-state booleans
    int  userUseInfinityNorm;           // -1 (not supplied), 0(false), 1(true)
    int  userReturnEveryInternalStep;   //      "
    int  userProjectEveryStep;          //      "
    int  userAllowInterpolation;        //      "
    int  userProjectInterpolatedStates; //      "
    int  userForceFullNewton;           //      "

    // Mark all user-supplied options "not supplied by user".
    void initializeUserStuff() {
        userInitStepSize = userMinStepSize = userMaxStepSize = -1.;
        userAccuracy = userConsTol = -1.;
        userFinalTime = -1.;
        userInternalStepLimit = -1;

        // booleans
        userUseInfinityNorm = userReturnEveryInternalStep = 
            userProjectEveryStep = userAllowInterpolation = 
            userProjectInterpolatedStates = userForceFullNewton = -1;

        accuracyInUse = NaN;
        consTol  = NaN;
    }

    // Required accuracy and constraint tolerance.
    // These are obtained during initialization from the user requests above, and
    // given default values if the user was agnostic.

    Real consTol;   // fraction of the constraint unit tolerance to be applied to each constraint

    void setAccuracyAndTolerancesFromUserRequests() {
        accuracyInUse = (userAccuracy != -1. ? userAccuracy : 1e-3);
        consTol       = (userConsTol  != -1. ? userConsTol  : 0.1*accuracyInUse); 
    }
    
    // If this is set to true, the next call to stepTo() will return immediately
    // with result code StartOfContinuousInterval. This is set by initialize()
    // and by reinitialize() after an event handler call that modified the
    // state.
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

    // State should have had its q's prescribed and realized through Position
    // stage. This will attempt to project q's and the q part of the yErrEst
    // (if yErrEst is not length zero). Returns false if we fail which you
    // can consider a convergence failure for the step. This is intended for
    // use during integration.
    // Stats are properly updated. State is realized through Position stage
    // on successful return.
    bool localProjectQAndQErrEstNoThrow(State& s, Vector& yErrEst,
        bool& anyChanges, Real projectionLimit=Infinity) 
    {
        ProjectOptions options;
        options.setRequiredAccuracy(getConstraintToleranceInUse());
        options.setProjectionLimit(projectionLimit);
        options.setOption(ProjectOptions::LocalOnly);
        options.setOption(ProjectOptions::DontThrow);
        if (userProjectEveryStep==1) 
            options.setOption(ProjectOptions::ForceProjection);
        if (userUseInfinityNorm==1)
            options.setOption(ProjectOptions::UseInfinityNorm);
        if (userForceFullNewton==1)
            options.setOption(ProjectOptions::ForceFullNewton);

        anyChanges = false;
        ProjectResults results;
        // Nothing happens here if position constraints were already satisfied
        // unless we set the ForceProjection option above.
        if (yErrEst.size()) {
            VectorView qErrEst = yErrEst(0, s.getNQ());
            getSystem().projectQ(s, qErrEst, options, results);
        } else {
            getSystem().projectQ(s, yErrEst, options, results);
        }
        if (results.getExitStatus() != ProjectResults::Succeeded) {
            ++statsQProjectionFailures;
            return false;
        }
        anyChanges = results.getAnyChangeMade();
        if (anyChanges)
            ++statsQProjections;

        return true;
    }

    // State should have had its q's and u's prescribed and realized through 
    // Velocity stage. This will attempt to project u's and the u part of the 
    // yErrEst (if yErrEst is not length zero). Returns false if we fail which 
    // you can consider a convergence failure for the step. This is intended for
    // use during integration.
    // Stats are properly updated. State is realized through Velocity stage
    // on successful return.
    bool localProjectUAndUErrEstNoThrow(State& s, Vector& yErrEst,
        bool& anyChanges, Real projectionLimit=Infinity) 
    {
        ProjectOptions options;
        options.setRequiredAccuracy(getConstraintToleranceInUse());
        options.setProjectionLimit(projectionLimit);
        options.setOption(ProjectOptions::LocalOnly);
        options.setOption(ProjectOptions::DontThrow);
        if (userProjectEveryStep==1) 
            options.setOption(ProjectOptions::ForceProjection);
        if (userUseInfinityNorm==1)
            options.setOption(ProjectOptions::UseInfinityNorm);
        if (userForceFullNewton==1)
            options.setOption(ProjectOptions::ForceFullNewton);

        anyChanges = false;
        ProjectResults results;
        // Nothing happens here if velocity constraints were already satisfied
        // unless we set the ForceProjection option above.
        if (yErrEst.size()) {
            VectorView uErrEst = yErrEst(s.getNQ(), s.getNU());
            getSystem().projectU(s, uErrEst, options, results);
        } else {
            getSystem().projectU(s, yErrEst, options, results);
        }
        if (results.getExitStatus() != ProjectResults::Succeeded) {
            ++statsUProjectionFailures;
            return false;
        }
        anyChanges = results.getAnyChangeMade();
        if (anyChanges)
            ++statsUProjections;

        return true;
    }

    // Given a state which has just had t and y updated, realize it through
    // velocity stage taking care of prescribed motion, projection, and 
    // throwing an exception if anything goes wrong. This is not for use 
    // during normal integration but good for initialization and generation
    // of interpolated states.
    // Extra options to consider are whether to restrict projection to the
    // local neighborhood, and whether to unconditionally force projection.
    // Stats are properly updated. State is realized through Velocity stage
    // on return.
    void realizeAndProjectKinematicsWithThrow(State& s,
        ProjectOptions::Option xtraOption1=ProjectOptions::None,
        ProjectOptions::Option xtraOption2=ProjectOptions::None,
        ProjectOptions::Option xtraOption3=ProjectOptions::None) 
    {
        const System& system = getSystem();
        ProjectOptions options;
        options.setRequiredAccuracy(getConstraintToleranceInUse());
        options.setOption(xtraOption1);
        options.setOption(xtraOption2);
        options.setOption(xtraOption3);
        if (userProjectEveryStep==1) 
            options.setOption(ProjectOptions::ForceProjection);
        if (userUseInfinityNorm==1)
            options.setOption(ProjectOptions::UseInfinityNorm);
        if (userForceFullNewton==1)
            options.setOption(ProjectOptions::ForceFullNewton);

        system.realize(s, Stage::Time);
        system.prescribeQ(s);
        system.realize(s, Stage::Position);

        Vector dummy; // no error estimate to project
        ProjectResults results;
        ++statsQProjectionFailures; // assume failure, then fix if no throw
        system.projectQ(s, dummy, options, results);
        --statsQProjectionFailures; // false alarm -- it succeeded
        if (results.getAnyChangeMade())
            ++statsQProjections;

        system.prescribeU(s);
        system.realize(s, Stage::Velocity);

        results.clear();
        ++statsUProjectionFailures; // assume failure, then fix if no throw
        system.projectU(s, dummy, options, results);
        --statsUProjectionFailures; // false alarm -- it succeeded
        if (results.getAnyChangeMade())
            ++statsUProjections;
    }

    // Set the advanced state and then evaluate state derivatives. Throws an
    // exception if it fails. Updates stats.
    void setAdvancedStateAndRealizeDerivatives(const Real& t, const Vector& y) 
    {
        const System& system = getSystem();
        State& advanced = updAdvancedState();

        setAdvancedState(t,y);

        system.realize(advanced, Stage::Time);
        system.prescribeQ(advanced); // set q_p
        system.realize(advanced, Stage::Position);
        system.prescribeU(advanced); // set u_p

        // Now realize Velocity, Dynamics, and Acceleration stages.
        realizeStateDerivatives(getAdvancedState());
    }

    // Set the advanced state and then evaluate constraint errors. Throws an
    // exception if it fails. Never counts as a realization because
    // we only need to realize kinematics.
    void setAdvancedStateAndRealizeKinematics(const Real& t, const Vector& y)
    {
        const System& system = getSystem();
        State& advanced = updAdvancedState();

        setAdvancedState(t,y);

        system.realize(advanced, Stage::Time);
        system.prescribeQ(advanced); // set q_p
        system.realize(advanced, Stage::Position);
        system.prescribeU(advanced); // set u_p

        // Now realize remaining kinematics.
        system.realize(advanced, Stage::Velocity);
    }


    void resetIntegratorStatistics() {
        statsRealizations = 0;
        statsQProjections = statsUProjections = 0;
        statsRealizationFailures = 0;
        statsQProjectionFailures = statsUProjectionFailures = 0;
    }

    int getNumRealizations() const {return statsRealizations;} 
    int getNumQProjections() const {return statsQProjections;}
    int getNumUProjections() const {return statsUProjections;}

    int getNumRealizationFailures() const {return statsRealizationFailures;} 
    int getNumQProjectionFailures() const {return statsQProjectionFailures;} 
    int getNumUProjectionFailures() const {return statsUProjectionFailures;} 

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
        const unsigned n = eventIds.size();
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
        for (unsigned i=0; i<n; ++i) 
            eventOrder[i] = events[i].ordinal;
    }

private:
    Integrator* myHandle;
    friend class Integrator;

    const System& sys;


protected:
    // realization and projection stats are shared by all integrators;
    // others are left to the individual integration methods
    mutable int statsQProjectionFailures, statsUProjectionFailures;
    mutable int statsQProjections, statsUProjections;
    mutable int statsRealizations;
    mutable int statsRealizationFailures;
private:

        // SYSTEM INFORMATION
        // Information extracted from the System describing properties we need
        // to know about BEFORE taking a time step. These properties are always
        // frozen across an integration step, but may be updated as discrete 
        // updates during time stepping.

    // Some Systems need to get control whenever time is advanced
    // successfully (and irreversibly) so that they can do discrete updates.
    // This is Stage::Model information.
    bool systemHasTimeAdvancedEvents;


    // Event trigger functions specify which sign transitions should be
    // monitored: rising, falling, to/from zero, or combinations of those.
    // They also provide a localization window width.
    // This is Stage::Instance information.

    Array_<EventTriggerInfo> eventTriggerInfo;

    // A unitless fraction.
    Real accuracyInUse;

    // The time scale is an estimate as to what is the smallest length of time
    // on which we can expect any significant changes to occur. This can be
    // used to estimate an initial step size, and it also establishes the units
    // in which event trigger localization windows are defined (that is, the 
    // window is given as some fraction of the time scale).
    // This is Stage::Instance information.
    Real timeScaleInUse;

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

    // These combine weightings from the Prev state with the possible
    // requirement of relative accuracy using the current values of u and z.
    // The result is min( wxi, 1/|xi|) for the relative ones where wxi is
    // the weight (1/absolute scale) and xi is either ui or zi.
    Vector  uScalePrev;
    Vector  zScalePrev;

    // Save continuous derivatives and event function values
    // from previous accepted state also so we don't have to
    // recalculate for restarts.
    Vector  ydotPrev;
    Vector  qdotdotPrev;
    Vector  triggersPrev;

    // These are views into yPrev and ydotPrev.
    Vector qPrev, uPrev, zPrev;
    Vector qdotPrev, udotPrev, zdotPrev;

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


