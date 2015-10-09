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
#include <iostream>
using std::cout; using std::endl;

namespace SimTK {

//==============================================================================
//                             INTEGRATOR REP
//==============================================================================
 
/* Specification of the stepTo() method below:

On entry we expect advancedState to be valid. Time may have already been 
advanced past the indicated reportTime in which case we'll be able to return an
interpolated result immediately. We will never advance past scheduledEventTime 
or finalTime.

The Integrator has a simple state machine to deal with the fact that the
computation is interruptable due to the need for returning control
to the caller (usually the TimeStepper) for various reasons. We enter
this computation either due to a call to the routine stepTo(), or from
within stepTo() after taking an internal step that didn't require a return.
Here's a sketch:

         advance (or report interpolated value)
          ------
         |      ^   
         v      |no  
  (TimeAdvanced)--return?->(StepReturned)--final?->(ReturnedFinalTime)
       ^                                 |no     
       <-------advance-------------------<

In words, we take a step, we report that step if necessary (returning
control), then take another step unless we just reported one at the
final time. The state transition table below gives more detail about
what exactly gets reported. "->Xyz" below means return control
to caller of stepTo() with status Xyz, "INTERP t" means "interpolate
back to time t".

There are additional states below used to distinguish between witnessed event
triggers (which take some special handling) and the others. The only difference
between state StepReturned and StepReturnedEvent is that some event-info query 
methods are allowed in the latter state.

Current State     Condition         Action           Next State
----------------  ----------------  --------------  ----------------
Uninitialized          any          throw error            X
ReturnedFinalTime      any          throw error            X

StepReturned or   t_adv>=t_final    ->Done         ReturnedFinalTime
StepReturnedEvent    OTHERWISE      ADVANCE       TimeAdvanced[Event]

TimeAdvancedEvent t_report<t_low    INTERP t_report  (unchanged)
                                    ->Report

                     OTHERWISE      INTERP t_low
                 (event triggered)  ->Triggered    StepReturnedEvent

TimeAdvanced      t_report<t_adv    INTERP t_report 
                                    ->Report          (unchanged)

                        else
                  t_adv==t_sched    ->Scheduled      StepReturned

                        else
                  user wants control
                    at end of step  ->EndOfStep      StepReturned

                        else
                  t_adv==t_report   ->Report         StepReturned

                        else
                  t_adv>=t_final    ->Final          StepReturned       

                        else
                  hit step limit    ->StepLimit      StepReturned
                
                     OTHERWISE      ADVANCE        TimeAdvanced[Event]
----------------  ----------------  --------------  ----------------

Notes:
* for the above conditions to be sufficient, the integrator must
  ensure that the interior of the interval [t_low,t_high] to
  which an event trigger has been localized DOES NOT contain
  t_report, t_sched, or t_final. It is OK if t_report is exactly
  t_low, and it is OK if t_report, t_sched, or t_final is exactly
  t_high.
* the intent of this state diagram is to make sure that the
  most significant end-of-step reason is reported, but that each
  step is reported at most once. Note that when the return status
  is "Done", the final trajectory point has *already* been reported.
  Simultaneous occurrences must be anticipated and handled by the caller
  (TimeStepper).
* note that an interpolated report doesn't count as reporting the
  step since we haven't reached t_advanced. But the interpolation
  to t_low does since t_advanced==t_high in that case and there
  is nothing in between t_low and t_high. Also, a report that
  coincides with the end of an internal step counts as a step report.
*/
class IntegratorRep {
public:
    /** This is an abbreviation for the type of an array of currently-active 
    event witnesses. **/
    using ActiveWitnessList = Array_<const EventWitness*, ActiveWitnessIndex>;
    /** This is an abbreviation for the type of an array used to represent a 
    subset of active event witnesses by indexing into an ActiveWitnessList. **/
    using ActiveWitnessSubset = Array_<ActiveWitnessIndex>;

    IntegratorRep(Integrator* handle, const System& system);
    virtual ~IntegratorRep() {}

    // no default constructor, no copy or copy assign
    IntegratorRep() = delete;
    IntegratorRep(const IntegratorRep&) = delete;
    IntegratorRep& operator=(const IntegratorRep&) = delete;

    // The System must be successfully realized to Stage::Instance before this
    // call. At this point the integrator can query the System about the
    // problem size and allocate appropriately-sized internal data structures.
    void initialize(const State&);

    // This is for use after return from an event handler.
    void reinitialize(Stage, bool shouldTerminate);

    // Make the integrator stop and invoke Termination event actions.
    void terminate(Integrator::TerminationReason reason);

    bool isSimulationOver() const 
    {   return m_stepCommunicationStatus == FinalTimeHasBeenReturned; }

    bool isStartOfContinuousInterval() const
    {   return m_startOfContinuousInterval; }
    
    // Get the reason the simulation ended. This should only be invoked if 
    // isSimulationOver() returns true.
    Integrator::TerminationReason getTerminationReason() const {
        assert(isSimulationOver());
        return m_terminationReason;
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

    const State& getAdvancedState() const {return m_advancedState;}
    double       getAdvancedTime()  const {return m_advancedState.getTime();}

    const State& getState() const 
    {   return m_useInterpolatedState ? m_interpolatedState : m_advancedState; }
    bool isStateInterpolated() const {return m_useInterpolatedState;}

    Real getAccuracyInUse() const {return m_accuracyInUse;}
    Real getConstraintToleranceInUse() const {return m_consTolInUse;}
    double getTimeScaleInUse() const {return m_timeScaleInUse;}

    // We also give the concrete integration method a chance to initialize 
    // itself after the generic initialization is done.
    virtual void methodInitialize(const State&) {}

    // This is for use after return from an event handler.
    virtual void methodReinitialize(Stage, bool) {}

    // In case there is any method-specific termination to do.
    virtual void methodTerminate(Integrator::TerminationReason) {}

    // The integrator has already been initialized. Take as many internal steps
    // as needed to get up to or past reportTime, but don't advance past the
    // next scheduled event time or the simulation final time.
    virtual Integrator::SuccessfulStepStatus 
        stepTo(double reportTime, double scheduledEventTime) = 0;

    // Size of the first successful step after the last initialize() call?
    virtual double getActualInitialStepSizeTaken() const = 0;

    // What was the size of the most recent successful step?
    virtual double getPreviousStepSizeTaken() const = 0;

    // What step size will be attempted first on the next step() call?
    virtual double getPredictedNextStepSize() const = 0;

    virtual SystemYIndex getPreviousStepWorstState() const
    {   return SystemYIndex(); } // invalid value

    virtual int getNumStepsAttempted() const = 0;
    virtual int getNumStepsTaken() const = 0; 
    virtual int getNumErrorTestFailures() const = 0;
    virtual int getNumConvergenceTestFailures() const = 0;
    virtual int getNumConvergentIterations() const = 0;
    virtual int getNumDivergentIterations() const = 0;
    virtual int getNumIterations() const = 0;

    virtual void resetMethodStatistics() {}

    virtual const char* getMethodName() const = 0;
    virtual int getMethodMinOrder() const = 0;
    virtual int getMethodMaxOrder() const = 0;
    virtual bool methodHasErrorControl() const = 0;

    // Cubic Hermite interpolation. See Hairer, et al. Solving ODEs I, 2nd rev.
    // ed., pg 190. Given (t0,y0,y0'),(t1,y1,y1') with y0 and y1 at least 3rd 
    // order accurate, obtain a 3rd order accurate interpolation yt for a time
    // t0 < t < t1 using Hermite interpolation.
    //
    // It is OK if yt is the same object as y0 or y1.
    static void interpolateOrder3
       (const double& t0, const Vector& y0, const Vector& f0,
        const double& t1, const Vector& y1, const Vector& f1,
        const double& t, Vector& yt);

    // We have bracketed a zero crossing for some function f(t) between 
    // (tLow,fLow) and (tHigh,fHigh), not including the end points. That means 
    // tHigh > tLow, and sign(fLow) != sign(fHigh) (sign(x) returns -1, 0, 1).
    // Estimate time tRoot with tLow<tRoot<tHigh such that f(tRoot) is zero.
    //
    // You can provide a bias (> 0) which will bias the returned tRoot into the
    // lower (bias<1) or upper (bias>1) half-interval. If bias==1 this is the 
    // pure secant method. 
    // 
    // Note that "minWindow" here is *not* the desired localization window
    // based on user accuracy requirements, but the (typically much smaller)
    // smallest allowable localization window based on numerical roundoff
    // considerations.
    static double estimateRootTime
       (double tLow,  const EventWitness::Value& fLow, 
        double tHigh, const EventWitness::Value& fHigh,
        Real bias, double minWindow);

    // Here we look at a pair of event witness function values across an 
    // interval and decide if there are any events triggering. If so we return
    // that witness as an "event candidate". Optionally, pass in the current 
    // list of event candidates and we'll only look at those (that is, the list
    // can only be narrowed). Don't use the same array for the current list and
    // new list.
    void findEventCandidates
       (const ActiveWitnessSubset*                  viableCandidates,
        const Array_<EventWitness::TransitionMask>* viableCandidateTransitions,
        double tLow,   const Array_<EventWitness::Value>&   eLow, 
        double tHigh,  const Array_<EventWitness::Value>&   eHigh,
        Real   bias,   double minWindow,
        ActiveWitnessSubset&                    candidates,
        Array_<double>&                         timeEstimates,
        Array_<EventWitness::TransitionMask>&   transitions,
        double&                                 earliestTimeEst, 
        double&                                 narrowestWindow) const;

    // Calculate the error norm using RMS or Inf norm, and report which y
    // was dominant.
    Real calcErrorNorm(const State& s, const Vector& yErrEst, 
                       int& worstY) const;


    // Given a proposed absolute change dq to generalized coordinates q, scale 
    // it to produce fq, such that fq_i is the fraction of q_i's "unit change"
    // represented by dq_i. A suitable norm of fq can then be compared directly
    // with the relative accuracy requirement. State must already be realized 
    // to Position stage. Note that dq scaling is based on u weighting Wu,
    // because u's and dq's are intimately related.
    void scaleDQ(const State& state, const Vector& Wu,
                 const Vector& dq, Vector& dqw) const; 

    // Calculate |Wq*dq|_RMS=|N*Wu*pinv(N)*dq|_RMS
    Real calcWeightedRMSNormQ(const State& state, const Vector& Wu,
                              const Vector& dq, int& worstQ) const {
        const int nq = state.getNQ();
        Vector dqw(nq);
        scaleDQ(state, Wu, dq, dqw);
        return dqw.normRMS(&worstQ);
    }
    // Calculate |Wq*dq|_Inf=|N*Wu*pinv(N)*dq|_Inf
    Real calcWeightedInfNormQ(const State& state, const Vector& Wu,
                              const Vector& dq, int& worstQ) const {
        const int nq = state.getNQ();
        Vector dqw(nq);
        scaleDQ(state, Wu, dq, dqw);
        return dqw.normInf(&worstQ);
    }


protected:
    const System& getSystem() const {return m_system;}
    const Integrator& getOwnerHandle() const {return *m_myHandle;}
    Integrator& updOwnerHandle() {return *m_myHandle;}

    // This is information we extract from the System during initialization or
    // reinitialization and save here.
    bool getSystemHasTimeAdvancedEvent() const {
        return m_systemHasTimeAdvancedEvent; 
    }
    double getSystemTimescale() const 
    {   return m_timeScaleInUse; }

    const State& getInterpolatedState() const {return m_interpolatedState;}
    State&       updInterpolatedState()       {return m_interpolatedState;}

    State& updAdvancedState() {return m_advancedState;}

    void setAdvancedState(const double& t, const Vector& y) {
        m_advancedState.updY() = y;
        m_advancedState.updTime() = t;
    }

    const ActiveWitnessList& getWitnesses() const {return m_witnesses;}
    ActiveWitnessList& updWitnesses() {return m_witnesses;}

    void setTriggeredEvents
       (double tlo, double thi,
        const Array_<ActiveWitnessIndex>&           triggeringWitnesses,
        const Array_<double>&                       estEventTimes,
        const Array_<EventWitness::TransitionMask>& transitionsSeen);

    Real getEventWindowLow()  const {return tLow;}
    Real getEventWindowHigh() const {return tHigh;}

    const Array_<const EventWitness*>&  getTriggeredWitnesses() const 
    {   return m_triggeredWitnesses; }
    const Array_<double>& getEstimatedTriggerTimes() const 
    {   return m_estimatedTriggerTimes; }
    const Array_<EventWitness::TransitionMask>& getWitnessTransitionsSeen() const 
    {   return m_witnessTransitionsSeen; }

    // This determines which state will be returned by getState().
    void setUseInterpolatedState(bool shouldUse) 
    {   m_useInterpolatedState = shouldUse; }
    void setStepCommunicationStatus(StepCommunicationStatus scs) 
    {   m_stepCommunicationStatus = scs; }
    StepCommunicationStatus getStepCommunicationStatus() const 
    {   return m_stepCommunicationStatus; }

    void setTerminationReason(Integrator::TerminationReason reason)
    {   m_terminationReason = reason; }

    void setIsStartOfContinuousInterval(bool isStart)
    {   m_startOfContinuousInterval = isStart; }

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

    const Array_<EventWitness::Value>& getPreviousWitnessValues() const 
    {   return m_witnessValuesPrev; }
    const Array_<int>& getPreviousWitnessDerivs() const 
    {   return m_witnessDerivsPrev; }
    const Array_<EventWitness::Value>& getAdvancedWitnessValues() const 
    {   return m_advancedWitnessValues; }
    Array_<EventWitness::Value>& updAdvancedWitnessValues() 
    {   return m_advancedWitnessValues; }

    // Given an array of state variables v (either u or z) and corresponding
    // weights w (1/absolute scale), return relative scale min(wi,1/|vi|) for
    // each variable i. That is, we choose the current value of vi as its
    // scale when it is large enough, otherwise use absolute scale.
    static void calcRelativeScaling
       (const Vector& v, const Vector& w, Vector& vScale) {
        const int nv = v.size();
        assert(w.size() == nv);
        vScale.resize(nv);
        for (int i=0; i<nv; ++i) {
            const Real vi = std::abs(v[i]);
            const Real wi = w[i];
            vScale[i] = vi*wi > 1 ? 1/vi : wi;
        }
    }

    // Calculate the values of all witnesses currently in m_witnesses, using
    // the value if it is non-zero, non-forbidden, otherwise the lowest 
    // non-zero, non-forbidden derivative.
    // Record which derivative was used to obtain value and sign.
    void calcStartStepWitnessValues(const State&                 state, 
                                    Array_<EventWitness::Value>& values,
                                    Array_<int>&                 whichDeriv) {
        values.resize(m_witnesses.size());
        whichDeriv.resize(m_witnesses.size());

        for (ActiveWitnessIndex awx(0); awx < m_witnesses.size(); ++awx) {
            const EventWitness& w = *m_witnesses[awx];
            const int nDerivs = w.getNumTimeDerivatives();
            for (int d=0; d <= nDerivs; ++d) {
                whichDeriv[awx] = d;
                values[awx] = w.calcWitnessValue(getOwnerHandle(), state, d);
                const int sign = values[awx].getSign();
                if (sign && !w.isSignForbidden(sign))
                    break;
            }
            // Highest deriv is accepted even if zero.
        }
    }

    // Calculate the end-of-step values of all witnesses currently in 
    // m_witnesses, using the same derivative as was used at step start.
    void calcEndStepWitnessValues(const State&                  state, 
                                  const Array_<int>&            whichDeriv,
                                  Array_<EventWitness::Value>&  values) {
        values.resize(m_witnesses.size());
        assert(whichDeriv.size()==m_witnesses.size());

        for (ActiveWitnessIndex awx(0); awx < m_witnesses.size(); ++awx) {
            const EventWitness& w = *m_witnesses[awx];
            const int d = whichDeriv[awx];
            values[awx] = w.calcWitnessValue(getOwnerHandle(), state, d);
        }
    }

    // Calculate the values of all witnesses currently in 
    // m_witnesses, ignoring derivatives.
    void calcWitnessValues(const State&                  state, 
                           Array_<EventWitness::Value>&  values) {
        values.resize(m_witnesses.size());

        for (ActiveWitnessIndex awx(0); awx < m_witnesses.size(); ++awx) {
            const EventWitness& w = *m_witnesses[awx];
            values[awx] = w.calcWitnessValue(getOwnerHandle(), state);
        }
    }

    // We're about to take a step. Record the current time and state from the
    // advanced state as the previous ones. We calculate the scaling for
    // u and z here which may include relative scaling based on their current
    // values. This scaling is frozen during a step attempt.
    void saveAdvancedTimeAndStateAsPrevious() {
        const State& s = getAdvancedState();
        const int nq = s.getNQ(), nu = s.getNU(), nz = s.getNZ();

        tPrev        = s.getTime();

        yPrev        = s.getY();
        qPrev.viewAssign(yPrev(0,     nq));
        uPrev.viewAssign(yPrev(nq,    nu));
        zPrev.viewAssign(yPrev(nq+nu, nz));

        calcRelativeScaling(s.getU(), s.getUWeights(), uScalePrev); 
        calcRelativeScaling(s.getZ(), s.getZWeights(), zScalePrev);
    }

    // We have evaluated state derivatives through Stage::Acceleration at the
    // start of a new step; record as previous values to make restarts easy. 
    // Must have called saveTimeAndStateAsPrevious() on this same state already.
    // We can now calculate the values of all witness functions and record them.
    void saveAdvancedStateDerivsAsPrevious() {
        const State& s = getAdvancedState();
        const int nq = s.getNQ(), nu = s.getNU(), nz = s.getNZ();

        ydotPrev     = s.getYDot();
        qdotPrev.viewAssign(ydotPrev(0,     nq));
        udotPrev.viewAssign(ydotPrev(nq,    nu));
        zdotPrev.viewAssign(ydotPrev(nq+nu, nz));

        qdotdotPrev  = s.getQDotDot();

        getSystem().findActiveEventWitnesses(*m_myHandle, m_witnesses);
        calcStartStepWitnessValues(s, m_witnessValuesPrev, m_witnessDerivsPrev);

#ifndef NDEBUG
        showWitnessValues("START OF STEP", s);
#endif
    }

    // For debugging.
    void showWitnessValues(const String& msg, const State& s) const {
        printf("%s: witness values at t=@%.15g\n", msg.c_str(), s.getTime());
        cout << "  q   =" << s.getQ() << endl;
        cout << "  qdot=" << s.getQDot() << endl;
        for (ActiveWitnessIndex awx(0); awx < m_witnesses.size(); ++awx) {
            const EventWitness& w = *m_witnesses[awx];
            const int whichD = m_witnessDerivsPrev[awx];
            for (int d=0; d <= w.getNumTimeDerivatives(); ++d) {
                const EventWitness::Value v = 
                    w.calcWitnessValue(getOwnerHandle(), s, d);
                const Real value = v.getValue();
                const int sign = v.getSign();
                if (isNaN(value)) continue;
                printf("  %2d %svalue=%g sign=%d",
                       (int)w.getEventTriggerId(),
                       d==0?"":(d==1?"d":"dd"),
                       value, sign);
                if (w.isSignForbidden(sign))
                    printf(" (FORBIDDEN)");
                if (d == whichD) printf(" <-- USING THIS ONE");
                printf("\n");
                if (d>=whichD && sign != 0 && !w.isSignForbidden(sign)) 
                    break;
            }
        }
    }

    // Advanced state must already have been evaluated through 
    // Stage::Acceleration or this will throw a stage violation.
    void saveAdvancedStateAndDerivsAsPrevious() {
        saveAdvancedTimeAndStateAsPrevious();
        saveAdvancedStateDerivsAsPrevious();
    }


    // Required accuracy and constraint tolerance.
    // These are obtained during initialization from the user requests, 
    // and given default values if the user was agnostic.
    void setAccuracyAndTolerancesFromUserRequests() {
        m_accuracyInUse = (userAccuracy != -1? userAccuracy:Real(1e-3));
        m_consTolInUse  = (userConsTol  != -1? userConsTol :m_accuracyInUse/10); 
    }


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
        bool& anyChanges, Real projectionLimit=Infinity);

    // State should have had its q's and u's prescribed and realized through 
    // Velocity stage. This will attempt to project u's and the u part of the 
    // yErrEst (if yErrEst is not length zero). Returns false if we fail which 
    // you can consider a convergence failure for the step. This is intended for
    // use during integration. Stats are properly updated. State is realized 
    // through Velocity stage on successful return.
    bool localProjectUAndUErrEstNoThrow(State& s, Vector& yErrEst,
        bool& anyChanges, Real projectionLimit=Infinity);

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
        ProjectOptions::Option xtraOption3=ProjectOptions::None);

    // Set the advanced state and then evaluate state derivatives. Throws an
    // exception if it fails. Updates stats.
    void setAdvancedStateAndRealizeDerivatives(const double& t, const Vector& y);

    // Set the advanced state and then evaluate constraint errors. Throws an
    // exception if it fails. Never counts as a realization because
    // we only need to realize kinematics.
    void setAdvancedStateAndRealizeKinematics(const double& t, const Vector& y);


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

protected:
    // collect user requests
    Real userInitStepSize, userMinStepSize, userMaxStepSize;
    Real userAccuracy; // also use for constraintTol
    Real userConsTol; // for fussy people
    Real userFinalTime; // never go past this
    int  userInternalStepLimit; // max in a single call to step(); 0=no limit

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

        m_accuracyInUse = NaN;
        m_consTolInUse  = NaN;
    }

    // realization and projection stats are shared by all integrators;
    // others are left to the individual integration methods
    mutable int statsQProjectionFailures, statsUProjectionFailures;
    mutable int statsQProjections, statsUProjections;
    mutable int statsRealizations;
    mutable int statsRealizationFailures;


private:
friend class Integrator;


    static void calcEventOrder
       (const Array_<ActiveWitnessIndex>&    triggeringWitnesses,
        const Array_<double>&                estEventTimes,
        Array_<unsigned>&                    witnessOrder);

    Integrator*     m_myHandle;
    const System&   m_system;


        // SYSTEM INFORMATION
        // Information extracted from the System describing properties we need
        // to know about BEFORE taking a time step. These properties are always
        // frozen across an integration step, but may be updated as discrete 
        // updates during time stepping.

    // Some Systems need to get control whenever time is advanced
    // successfully (and irreversibly) so that they can do discrete updates.
    // This is Stage::Model information.
    bool m_systemHasTimeAdvancedEvent;


    // Event trigger functions specify which sign transitions should be
    // monitored: rising, falling, to/from zero, or combinations of those.
    // They also provide a localization window width.
    // This is Stage::Instance information.

    //Array_<EventTriggerInfo> eventTriggerInfo;

    // A unitless fraction.
    Real m_accuracyInUse;

    // Fraction of constraint unit tolerance to be applied to each constraint.
    Real m_consTolInUse; 

    // The time scale is an estimate as to what is the smallest length of time
    // on which we can expect any significant changes to occur. This can be
    // used to estimate an initial step size, and it also establishes the units
    // in which event trigger localization windows are defined (that is, the 
    // window is given as some fraction of the time scale).
    // This is Stage::Instance information.
    double                          m_timeScaleInUse;

        // INTEGRATOR INTERNAL STATE 
        // Persists between stepTo() calls.
        // A concrete integrator is free to have more state.

    StepCommunicationStatus         m_stepCommunicationStatus;
    
    // If this is set to true, the next call to stepTo() will return immediately
    // with result code StartOfContinuousInterval. This is set by initialize()
    // and by reinitialize() after an event handler call that modified the
    // state.
    bool                            m_startOfContinuousInterval;
    
    // The reason the simulation ended.  If isSimulationOver() returns false,
    // the value of this field is meaningless.
    Integrator::TerminationReason   m_terminationReason;

    // We save both the step size we have to take next and the one
    // we wish we could take. Once we're done with whatever's causing
    // us to come up short (like event isolation), we may be able to
    // jump right back up to the ideal one.
    double  m_nextStepSizeToTry;  // This is what we'll actually try first.
    double  m_idealNextStepSize;  // But if accuracy were the only concern
                                  //   we would try this (>=nextStepSizeToTry).

    // This is how far the integrator has advanced the system
    // trajectory. We might allow an interpolated state a
    // little earlier than this, but otherwise there is no
    // going back from this point.
    State   m_advancedState;

    // This is a temporary in which to evaluate witness functions. The 
    // start-of-step value was saved in m_witnessValuesPrev.
    Array_<EventWitness::Value>             m_advancedWitnessValues;

    // When stepCommunicationStatus indicates that an event has triggered due to
    // a witness function zero crossing, the following arrays are valid.
    // Witnesses are sorted by estimated occurrence time within the window; 
    // identical-time witnesses are sorted in ascending order 
    // of ActiveWitnessIndex.

    // These are the witnesses that the integrator has algorithmically
    // determined have triggered in the current step. You may not get the same 
    // result simply comparing the witness function values at tLow and tHigh.
    Array_<const EventWitness*>             m_triggeredWitnesses;

    // These are the estimated times corresponding to the triggeredEvents.
    // They are in ascending order although there may be duplicates.
    Array_<double>                          m_estimatedTriggerTimes;
    
    // Which transition (rising or falling) was seen for each triggered witness
    // (this is only a single transition, not an OR-ed together set). This is
    // the integrator's algorithmic determination of the transition to
    // be reported -- you might not get the same answer just looking
    // at the event trigger functions at tLow and tHigh.
    Array_<EventWitness::TransitionMask>    m_witnessTransitionsSeen;

    // When we have successfully localized a triggering event into
    // the time interval (tLow,tHigh] we record the bounds here.
    double  tLow, tHigh;

    State   m_interpolatedState;    // might be unused
    bool    m_useInterpolatedState;

    // Use these to record the continuous part of the previous
    // accepted state. We use these in combination with the 
    // continuous contents of advancedState to fill in the
    // interpolated state.
    double  tPrev;
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

    // These are the witnesses we will consider during the current continuous
    // step. We'll record the values at the beginning of the step and then
    // compare those with recalculated end-of-step values to see if an event
    // has occurred. We'll reference these by "active witness index" 
    // [0..m_witnesses.size()-1].
    Array_<const EventWitness*, 
           ActiveWitnessIndex>      m_witnesses;
    Array_<EventWitness::Value>     m_witnessValuesPrev; // same length
    Array_<int>                     m_witnessDerivsPrev; // which deriv used?

    // These are views into yPrev and ydotPrev.
    Vector qPrev, uPrev, zPrev;
    Vector qdotPrev, udotPrev, zdotPrev;

    // We'll leave the various arrays above sized as they are and full
    // of garbage. They'll be resized when first assigned to something
    // meaningful.
    void invalidateIntegratorInternalState() {
        m_stepCommunicationStatus = InvalidStepCommunicationStatus;
        m_nextStepSizeToTry       = dNaN; // times are always double
        m_idealNextStepSize       = dNaN;
        tLow = tHigh              = dNaN;
        m_useInterpolatedState    = false;
        tPrev                     = dNaN;
    }
};



} // namespace SimTK

#endif // SimTK_SIMMATH_INTEGRATOR_REP_H_


