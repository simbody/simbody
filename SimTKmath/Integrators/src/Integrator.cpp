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
 * This is the private (library side) implementation of the Simmath
 * Integrator family of classes.
 */

#include "SimTKcommon.h"
#include "simmath/Integrator.h"

#include "IntegratorRep.h"

#include <exception>
#include <stdexcept>
#include <limits>
#include <iostream>
#include <algorithm>

using std::cout; using std::endl;

using namespace SimTK;

//==============================================================================
//                                INTEGRATOR
//==============================================================================


Integrator::~Integrator() {
    delete rep; rep=0;
}

void Integrator::resetAllStatistics() {
    updRep().resetIntegratorStatistics();
    updRep().resetMethodStatistics();
}

void Integrator::initialize(const State& initState) {
    updRep().initialize(initState);
}

void Integrator::reinitialize(Stage g, bool shouldTerminate) {
    updRep().reinitialize(g,shouldTerminate);
}

void Integrator::terminate(TerminationReason reason) {
    updRep().terminate(reason);
}


Integrator::SuccessfulStepStatus 
Integrator::stepTo(double reportTime, double advanceLimit) {
    return updRep().stepTo(reportTime, advanceLimit);
}

Integrator::SuccessfulStepStatus 
Integrator::stepBy(double interval, double advanceIntervalLimit) {
    const double t = getRep().getState().getTime();
    return updRep().stepTo(t + interval, t + advanceIntervalLimit);
}

bool Integrator::isSimulationOver() const {
    return getRep().isSimulationOver();
}

Integrator::TerminationReason Integrator::getTerminationReason() const {
    return getRep().getTerminationReason();
}

/*static*/ String Integrator::
getSuccessfulStepStatusString(SuccessfulStepStatus stat) {
    switch(stat) { 
        case ReachedReportTime: return "ReachedReportTime";
        case ReachedEventTrigger: return "ReachedEventTrigger";
        case ReachedScheduledEvent: return "ReachedScheduledEvent";
        case TimeHasAdvanced: return "TimeHasAdvanced";
        case ReachedStepLimit: return "ReachedStepLimit";
        case EndOfSimulation: return "EndOfSimulation";
        case StartOfContinuousInterval: return "StartOfContinuousInterval";
        case InvalidSuccessfulStepStatus: return "InvalidSuccessfulStepStatus";
        default: return String("UNKNOWN SUCCESSFUL STEP STATUS ") 
                      + String((int)stat);
    }
}

/*static*/ String Integrator::
getTerminationReasonString(TerminationReason reason) {
    switch(reason) { 
        case ReachedFinalTime: return "ReachedFinalTime";
        case AnUnrecoverableErrorOccurred: return "AnUnrecoverableErrorOccurred";
        case EventHandlerRequestedTermination: return "EventHandlerRequestedTermination";
        case InvalidTerminationReason: return "InvalidTerminationReason";
        default: return String("UNKNOWN TERMINATION REASON ") 
                      + String((int)reason);
    }
}

// The following methods are callable only when the Integrator has just
// returned a step that ended with an event triggering.

Vec2 Integrator::getEventWindow() const {
    if (getRep().getStepCommunicationStatus() != IntegratorRep::StepHasBeenReturnedWithEvent) {
        SimTK_THROW2(CantAskForEventInfoWhenNoEventTriggered, "getEventWindow",
                     getRep().getState().getTime());
        //NOTREACHED
    }
    assert(getRep().getEventWindowLow()==getRep().getState().getTime());
    assert(getRep().getEventWindowHigh()==getRep().getAdvancedTime());
    return Vec2(getRep().getEventWindowLow(), getRep().getEventWindowHigh());
}

const Array_<const EventTrigger::Witness*>& 
Integrator::getTriggeredWitnesses() const {
    if (getRep().getStepCommunicationStatus() 
        != IntegratorRep::StepHasBeenReturnedWithEvent) {
        SimTK_THROW2(CantAskForEventInfoWhenNoEventTriggered, 
                     "getTriggeredEvents", getRep().getState().getTime());
        //NOTREACHED
    }
    return getRep().getTriggeredWitnesses();
}

const Array_<Real>&
Integrator::getEstimatedTriggerTimes() const {
    if (getRep().getStepCommunicationStatus() 
        != IntegratorRep::StepHasBeenReturnedWithEvent) {
        SimTK_THROW2(CantAskForEventInfoWhenNoEventTriggered, 
                     "getEstimatedEventTimes", getRep().getState().getTime());
        //NOTREACHED
    }
    return getRep().getEstimatedTriggerTimes();
}

const Array_<Event::TriggerDirection>&
Integrator::getWitnessTransitionsSeen() const {
    if (getRep().getStepCommunicationStatus() 
        != IntegratorRep::StepHasBeenReturnedWithEvent) {
        SimTK_THROW2(CantAskForEventInfoWhenNoEventTriggered, 
                     "getEventTransitionsSeen", getRep().getState().getTime());
        //NOTREACHED
    }
    return getRep().getWitnessTransitionsSeen();
}

const State& Integrator::getState() const {
    return getRep().getState();
}

const System& Integrator::getSystem() const {
    return getRep().getSystem();
}

bool Integrator::isStateInterpolated() const {
    return getRep().isStateInterpolated();
}

const State& Integrator::getAdvancedState() const {
    return getRep().getAdvancedState();
}

State& Integrator::updAdvancedState() {
    return updRep().updAdvancedState();
}


Real Integrator::getAccuracyInUse() const {
    return getRep().getAccuracyInUse();
}

Real Integrator::getConstraintToleranceInUse() const {
    return getRep().getConstraintToleranceInUse();
}

double Integrator::getActualInitialStepSizeTaken() const {
    return getRep().getActualInitialStepSizeTaken();
}
double Integrator::getPreviousStepSizeTaken() const {
    return getRep().getPreviousStepSizeTaken();
}
SystemYIndex Integrator::getPreviousStepWorstState() const {
    return getRep().getPreviousStepWorstState();
}
double Integrator::getPredictedNextStepSize() const {
    return getRep().getPredictedNextStepSize();
}
int Integrator::getNumStepsAttempted() const {
    return getRep().getNumStepsAttempted();
}
int Integrator::getNumStepsTaken() const {
    return getRep().getNumStepsTaken();
}
int Integrator::getNumRealizations() const {
    return getRep().getNumRealizations();
}
int Integrator::getNumQProjections() const {
    return getRep().getNumQProjections();
}
int Integrator::getNumUProjections() const {
    return getRep().getNumUProjections();
}
int Integrator::getNumProjections() const {
    return getRep().getNumQProjections()
         + getRep().getNumUProjections();
}
int Integrator::getNumErrorTestFailures() const {
    return getRep().getNumErrorTestFailures();
}
int Integrator::getNumConvergenceTestFailures() const {
    return getRep().getNumConvergenceTestFailures();
}
int Integrator::getNumRealizationFailures() const {
    return getRep().getNumRealizationFailures();
}
int Integrator::getNumQProjectionFailures() const {
    return getRep().getNumQProjectionFailures();
}
int Integrator::getNumUProjectionFailures() const {
    return getRep().getNumUProjectionFailures();
}
int Integrator::getNumProjectionFailures() const {
    return getRep().getNumQProjectionFailures()
         + getRep().getNumUProjectionFailures();
}
int Integrator::getNumConvergentIterations() const {
    return getRep().getNumConvergentIterations();
}
int Integrator::getNumDivergentIterations() const {
    return getRep().getNumDivergentIterations();
}
int Integrator::getNumIterations() const {
    return getRep().getNumIterations();
}

void Integrator::setFinalTime(Real tFinal) {
    assert(tFinal == -1. || (0. <= tFinal));
    updRep().userFinalTime = tFinal;
}

void Integrator::setInternalStepLimit(int nSteps) {
    updRep().userInternalStepLimit = nSteps > 0 ? nSteps : -1;
}

void Integrator::setInitialStepSize(Real z) {
    assert(z == -1. || z > 0.);
    assert(getRep().userMinStepSize==-1. || z >= getRep().userMinStepSize);
    assert(getRep().userMaxStepSize==-1. || z <= getRep().userMaxStepSize);
    updRep().userInitStepSize = z;
}
void Integrator::setMinimumStepSize(Real z) { 
    assert(z == -1. || z > 0.);
    assert(getRep().userInitStepSize==-1. || z <= getRep().userInitStepSize);
    assert(getRep().userMaxStepSize ==-1. || z <= getRep().userMaxStepSize);
    updRep().userMinStepSize = z;
}
void Integrator::setMaximumStepSize(Real z) {
    assert(z == -1. || z > 0.);
    assert(getRep().userInitStepSize==-1. || z >= getRep().userInitStepSize);
    assert(getRep().userMinStepSize ==-1. || z >= getRep().userMinStepSize);
    updRep().userMaxStepSize = z;
}
void Integrator::setFixedStepSize(Real stepSize) {
    updRep().userInitStepSize = stepSize;
    updRep().userMinStepSize = stepSize;
    updRep().userMaxStepSize = stepSize;
}
void Integrator::setAccuracy(Real accuracy) {
    assert(accuracy == -1. || (0. < accuracy && accuracy < 1.));
    updRep().userAccuracy = accuracy;
}
void Integrator::setConstraintTolerance(Real consTol) {
    assert(consTol == -1. || (0. < consTol && consTol <= 1.));
    updRep().userConsTol=consTol;
}
void Integrator::setUseInfinityNorm(bool useInfinityNorm) {
    updRep().userUseInfinityNorm = useInfinityNorm ? 1 : 0;
}
bool Integrator::isInfinityNormInUse() const
{   return getRep().userUseInfinityNorm == 1; }

void Integrator::setForceFullNewton(bool forceFullNewton) {
    updRep().userForceFullNewton = forceFullNewton ? 1 : 0;
}
void Integrator::setReturnEveryInternalStep(bool shouldReturn) {
    updRep().userReturnEveryInternalStep = shouldReturn ? 1 : 0;
}
void Integrator::setProjectEveryStep(bool forceProject) {
    updRep().userProjectEveryStep = forceProject ? 1 : 0;
}
void Integrator::setAllowInterpolation(bool shouldInterpolate) {
    updRep().userAllowInterpolation = shouldInterpolate ? 1 : 0;
}
void Integrator::setProjectInterpolatedStates(bool shouldProject) {
    updRep().userProjectInterpolatedStates = shouldProject ? 1 : 0;
}

bool Integrator::methodHasErrorControl() const {
    return getRep().methodHasErrorControl();
}

const char* Integrator::getMethodName() const {
    return getRep().getMethodName();
}

int Integrator::getMethodMinOrder() const {
    return getRep().getMethodMinOrder();
}

int Integrator::getMethodMaxOrder() const {
    return getRep().getMethodMaxOrder();
}

//==============================================================================
//                              INTEGRATOR REP
//==============================================================================

IntegratorRep::IntegratorRep
       (Integrator*               handle,
        const System&             system)
:   m_myHandle(handle), m_system(system)
{
    invalidateIntegratorInternalState();
    initializeUserStuff();
    resetIntegratorStatistics();
    resetMethodStatistics();
}

//------------------------------------------------------------------------------
//                               INITIALIZE
//------------------------------------------------------------------------------

void IntegratorRep::initialize(const State& initState) {
  try
  { invalidateIntegratorInternalState();

    // Copy the supplied initial state into the integrator.
    updAdvancedState() = initState;

    // Freeze the number and kinds of state variables.
    getSystem().realizeModel(updAdvancedState());
     
    // Freeze problem dimensions; at this point the state represents an
    // instance of a physically realizable system.
    getSystem().realize(getAdvancedState(), Stage::Instance);

    // Some Systems need to get control whenever time is advanced
    // successfully (and irreversibly) so that they can do discrete updates.
    m_systemHasTimeAdvancedEvent = getSystem().hasTimeAdvancedEvents();

    // Allocate integrator-local data structures now that we know the sizes.
    const int ny = getAdvancedState().getNY();
    const int nc = getAdvancedState().getNYErr();
    m_timeScaleInUse = getSystem().getDefaultTimeScale();

    // Set accuracy and consTol to their user-requested values or
    // to the appropriate defaults.
    setAccuracyAndTolerancesFromUserRequests();

    // This will be set to something meaningful when the simulation ends.
    setTerminationReason(Integrator::InvalidTerminationReason);

    // Realize Time stage, apply prescribed motions, and attempt to project q's 
    // and u's onto the
    // position and velocity constraint manifolds to drive constraint errors
    // below accuracy*unitTolerance for each constraint. We'll allow project()
    // to throw an exception if it fails since we can't recover from here.
    // However, we won't set the LocalOnly option which means project() is
    // allowed to thrash around wildly to attempt to find *some* solution.
    // Also force repeated updates of iteration matrix to maximize chances of 
    // finding a solution; we're not in a hurry here.
    realizeAndProjectKinematicsWithThrow(updAdvancedState(),
        ProjectOptions::ForceProjection, ProjectOptions::ForceFullNewton);

    // Inform any state-using System elements (especially Measures) that we 
    // are starting a simulation and give them a chance to initialize their own
    // continuous (z) variables and discrete variables.

    // If an Initialization event action throws an exception we're dead since
    // we have no way to recover here.
    EventsAndCauses triggeredEvents;
    Array_<EventId> ignoredEventIds;
    getSystem().noteEventOccurrence({&getSystem().getInitializationTrigger()},
                                    triggeredEvents, ignoredEventIds);

    EventChangeResult result;
    getSystem().performEventChangeActions
       (updOwnerHandle(), triggeredEvents, result);

    // A failure during initialization isn't recoverable, and a request to
    // terminate at initialization is, well, odd.
    if (result.getExitStatus() != EventChangeResult::Succeeded) {
        String reason("  An initialization event handler ");

        switch (result.getExitStatus()) {
        case EventChangeResult::ShouldTerminate:
            reason += "requested termination ";
            reason += result.getMessage().empty() ? String("(no message).")
                : "with message '" + result.getMessage() + "'.";
            break;
        case EventChangeResult::Failed:
            reason += "failed ";
            reason += result.getMessage().empty() 
                ? String("(no message).")
                : "with message '" + result.getMessage() + "'.";
            break;
        default:
            reason += "returned an unknown exit status " 
                      + String((int)result.getExitStatus()) + ".";
        }

        throw std::runtime_error(reason);
    }

    // Now evaluate the state through the Acceleration stage to calculate
    // the initial state derivatives.
    realizeStateDerivatives(getAdvancedState());

    // Now that we have valid update values for auto-update discrete variables,
    // use them to reinitialize the discrete state variables. This one time
    // only, the value swapped in from the discrete variable here is not yet
    // suitable for use as an update value, during time integration since it
    // has never been realized to Instance stage. So we'll force a 
    // re-evaluation.
    updAdvancedState().autoUpdateDiscreteVariables();
    getAdvancedState().invalidateAllCacheAtOrAbove(Stage::Instance);
    // Re-realize to fill in the swapped-in update values.
    realizeStateDerivatives(getAdvancedState());

    // If there are initialization report actions, perform them now.
    getSystem().performEventReportActions(getOwnerHandle(), triggeredEvents);

    // Record the continuous parts of this now-realized initial state as the 
    // previous state as well (previous state is used when we have to back up
    // from a failed step attempt).
    saveAdvancedStateAndDerivsAsPrevious();

    // The initial state is set so it looks like we just *completed* a step to 
    // get here. That way if the first reportTime is zero, this will get 
    // reported.
    setStepCommunicationStatus(CompletedInternalStepNoEvent);
    setIsStartOfContinuousInterval(true);

    // Call this virtual method in case the concrete integrator class has 
    // additional initialization needs. Note that it gets the State as we have
    // adjusted it, NOT the one originally passed to initialize().
    methodInitialize(getAdvancedState());

  } catch (const std::exception& e) {
    SimTK_THROW1(Integrator::InitializationFailed, e.what());
  } /* catch (...) {
    SimTK_THROW1(Integrator::InitializationFailed, "UNKNOWN EXCEPTION");
  } */
}

//------------------------------------------------------------------------------
//                               REINITIALIZE
//------------------------------------------------------------------------------
void IntegratorRep::reinitialize(Stage stage, bool shouldTerminate) {
    if (stage < Stage::Report) {
        setIsStartOfContinuousInterval(true);
        setUseInterpolatedState(false);
    }

    if (shouldTerminate)
        terminate(Integrator::EventHandlerRequestedTermination);

    methodReinitialize(stage,shouldTerminate);
}

//------------------------------------------------------------------------------
//                               TERMINATE
//------------------------------------------------------------------------------
void IntegratorRep::terminate(Integrator::TerminationReason reason) {
    const System& sys = getSystem();

    setStepCommunicationStatus(FinalTimeHasBeenReturned);
    setUseInterpolatedState(false);
    setTerminationReason(reason);

    EventsAndCauses triggeredEvents;
    Array_<EventId> ignoredEventIds;
    getSystem().noteEventOccurrence({&getSystem().getTerminationTrigger()},
                                    triggeredEvents, ignoredEventIds);

    EventChangeResult result;
    sys.performEventChangeActions
        (updOwnerHandle(), triggeredEvents, result);
    // Ignore handler failure and requests for termination here.

    // Make sure the final state is fully realized.
    sys.realize(getAdvancedState(), Stage::Report);
    sys.performEventReportActions(getOwnerHandle(), triggeredEvents);

    methodTerminate(reason);
}

//------------------------------------------------------------------------------
//                          FIND EVENT CANDIDATES
//------------------------------------------------------------------------------
/* Here we look at pairs of event witness function values across an interval 
and decide if there are any events triggering. If so we return those witness
triggers as "event candidates". Optionally, pass in the current list of event 
candidates and we'll only look at those (that is, the list can only be 
narrowed). Don't use the same array for the current list and new list.

For purposes of this method, events are specified by their indices in 
the array of trigger functions, NOT by their event IDs. */
void IntegratorRep::findEventCandidates
   (const ActiveWitnessSubset*              viableCandidates,
    const Array_<Event::TriggerDirection>*  viableCandidateTransitions,
    double tLow,   const Vector&   eLow,  // [nWitnesses]
    double tHigh,  const Vector&   eHigh, // [nWitnesses]
    Real   bias,   double          minWindow,
    ActiveWitnessSubset&                    candidates,
    Array_<double>&                         timeEstimates,
    Array_<Event::TriggerDirection>&        transitions,
    double&                                 earliestTimeEst, 
    double&                                 narrowestWindow) const
{
    const int nWitnesses = (int)m_witnesses.size();

    int nCandidates;
    if (viableCandidates) {
        nCandidates = (int)viableCandidates->size();
        assert(   viableCandidateTransitions 
                && (int)viableCandidateTransitions->size()==nCandidates);
    } else {
        assert(!viableCandidateTransitions);
        nCandidates = nWitnesses;
    }

    candidates.clear();
    timeEstimates.clear();
    transitions.clear();
    earliestTimeEst = narrowestWindow = dInfinity;
    for (int i=0; i<nCandidates; ++i) {
        const ActiveWitnessIndex awx = viableCandidates ?
                        (*viableCandidates)[i] : ActiveWitnessIndex(i);
        const EventTrigger::Witness& aw = *m_witnesses[awx];

        Event::TriggerDirection transitionSeen =
            Event::maskTransition(
                Event::classifyTransition(sign(eLow[awx]), sign(eHigh[awx])),
                aw.calcTransitionMask());

        if (transitionSeen != Event::NoEventTrigger) {
            candidates.push_back(awx);
            const Real relWindow = // unitless scale factor
                aw.getAccuracyRelativeTimeLocalizationWindow();
            const double absWindow = 
                (m_accuracyInUse*relWindow)*m_timeScaleInUse;
            narrowestWindow = 
                std::max(std::min(narrowestWindow, absWindow),
                         minWindow);

            // Set estimated event trigger time for the viable candidates.
            timeEstimates.push_back
                        (estimateRootTime(tLow, eLow[awx], tHigh, eHigh[awx],
                                          bias, minWindow));
            transitions.push_back(transitionSeen);
            earliestTimeEst = std::min(earliestTimeEst, timeEstimates.back());
        }
    }
}


//------------------------------------------------------------------------------
//                           SET TRIGGERED EVENTS
//------------------------------------------------------------------------------
void IntegratorRep::
setTriggeredEvents(double tlo, double thi,
                   const Array_<ActiveWitnessIndex>&      triggeringWitnesses,
                   const Array_<double>&                  estEventTimes,
                   const Array_<Event::TriggerDirection>& transitionsSeen)
{
    assert(tPrev <= tlo && tlo < thi && thi <= m_advancedState.getTime());
    tLow = tlo;
    tHigh = thi;

    const int n = triggeringWitnesses.size();
    assert(n > 0 && estEventTimes.size()==n && transitionsSeen.size()==n);
    m_triggeredWitnesses.resize(n); m_estimatedTriggerTimes.resize(n); 
    m_witnessTransitionsSeen.resize(n);
    Array_<unsigned> triggerOrder; // will be a permutation of 0:n-1
    calcEventOrder(triggeringWitnesses, estEventTimes, triggerOrder);
    for (unsigned i=0; i < triggerOrder.size(); ++i) {
        const unsigned ipos = triggerOrder[i];
        m_triggeredWitnesses[i] = m_witnesses[triggeringWitnesses[ipos]];

        assert(tlo < estEventTimes[ipos] && estEventTimes[ipos] <= thi);
        m_estimatedTriggerTimes[i] = estEventTimes[ipos];

        // TODO: assert that the transition is one of the allowed ones for
        // this event.
        m_witnessTransitionsSeen[i] = transitionsSeen[ipos];
    }
}


//------------------------------------------------------------------------------
//                             CALC ERROR NORM
//------------------------------------------------------------------------------
/* Calculate error norm using RMS or Inf norm, and report which y was the worst
offender. */
Real IntegratorRep::
calcErrorNorm(const State& s, const Vector& yErrEst, int& worstY) const {
    const int nq=s.getNQ(), nu=s.getNU(), nz=s.getNZ();
    int worstQ, worstU, worstZ;
    Real qNorm, uNorm, zNorm, maxNorm;

    if (userUseInfinityNorm == 1) {
        qNorm = calcWeightedInfNormQ(s, s.getUWeights(), yErrEst(0,nq), worstQ);
        uNorm = yErrEst(nq,   nu).weightedNormInf(getPreviousUScale(), &worstU);
        zNorm = yErrEst(nq+nu,nz).weightedNormInf(getPreviousZScale(), &worstZ);
    } else {
        qNorm = calcWeightedRMSNormQ(s, s.getUWeights(), yErrEst(0,nq), worstQ);
        uNorm = yErrEst(nq,   nu).weightedNormRMS(getPreviousUScale(), &worstU);
        zNorm = yErrEst(nq+nu,nz).weightedNormRMS(getPreviousZScale(), &worstZ);
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


//------------------------------------------------------------------------------
//                                SCALE DQ
//------------------------------------------------------------------------------
/* Given a proposed absolute change dq to generalized coordinates q, scale it to
produce fq, such that fq_i is the fraction of q_i's "unit change" represented by
dq_i. A suitable norm of fq can then be compared directly with the relative 
accuracy requirement. 

Because dq's are not independent of u's (qdot=N(q)*u), the "unit change" of q is
related to the unit change of u. We want fq=Wq*dq, but we determine Wq from Wu 
via Wq = N*Wu*pinv(N). Wq is block diagonal while Wu is diagonal. State must 
already be realized to Position stage. */
void IntegratorRep::
scaleDQ(const State& state, const Vector& Wu,
        const Vector& dq, Vector& dqw) const
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


//------------------------------------------------------------------------------
//                   LOCAL PROJECT Q AND QERREST NO THROW
//------------------------------------------------------------------------------
/* State should have had its q's prescribed and realized through Position
stage. This will attempt to project q's and the q part of the yErrEst (if 
yErrEst is not length zero). Returns false if we fail which you can consider a 
convergence failure for the step. This is intended for use during integration.
Stats are properly updated. State is realized through Position stage on 
successful return. */
bool IntegratorRep::
localProjectQAndQErrEstNoThrow(State& s, Vector& yErrEst,
                               bool& anyChanges, Real projectionLimit) 
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


//------------------------------------------------------------------------------
//                   LOCAL PROJECT U AND UERREST NO THROW
//------------------------------------------------------------------------------
/* State should have had its q's and u's prescribed and realized through 
Velocity stage. This will attempt to project u's and the u part of the yErrEst
(if yErrEst is not length zero). Returns false if we fail which you can consider
a convergence failure for the step. This is intended for use during integration.
Stats are properly updated. State is realized through Velocity stage on 
successful return. */
bool IntegratorRep::
localProjectUAndUErrEstNoThrow(State& s, Vector& yErrEst,
                               bool& anyChanges, Real projectionLimit) 
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


//------------------------------------------------------------------------------
//                  REALIZE AND PROJECT KINEMATICS WITH THROW
//------------------------------------------------------------------------------
/* Given a state which has just had t and y updated, realize it through velocity
stage taking care of prescribed motion, projection, and throwing an exception if
anything goes wrong. This is not for use during normal integration but good for
initialization and generation of interpolated states. Extra options to consider
are whether to restrict projection to the local neighborhood, and whether to 
unconditionally force projection. Stats are properly updated. State is realized
through Velocity stage on return. */
void IntegratorRep::
realizeAndProjectKinematicsWithThrow(State&                 s,
                                     ProjectOptions::Option xtraOption1,
                                     ProjectOptions::Option xtraOption2,
                                     ProjectOptions::Option xtraOption3) 
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


//------------------------------------------------------------------------------
//                  SET ADVANCED STATE AND REALIZE DERIVATIVES
//------------------------------------------------------------------------------
/* Set the advanced state and then evaluate state derivatives. Throws an
exception if it fails. Updates stats. */
void IntegratorRep::
setAdvancedStateAndRealizeDerivatives(const double& t, const Vector& y) {
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


//------------------------------------------------------------------------------
//                  SET ADVANCED STATE AND REALIZE KINEMATICS
//------------------------------------------------------------------------------
/* Set the advanced state and then evaluate constraint errors. Throws an
exception if it fails. Never counts as a realization because we only need to 
realize kinematics. */
void IntegratorRep::
setAdvancedStateAndRealizeKinematics(const double& t, const Vector& y) {
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


//------------------------------------------------------------------------------
//                            INTERPOLATE ORDER 3
//------------------------------------------------------------------------------
/* Cubic Hermite interpolation. See Hairer, et al. Solving ODEs I, 2nd rev. ed.,
pg 190. Given (t0,y0,y0'),(t1,y1,y1') with y0 and y1 at least 3rd order 
accurate, we can obtain a 3rd order accurate interpolation yt for a time 
t0 < t < t1 using Hermite interpolation. Let f0=y0', f1=y1', h=t1-t0, 
d=(t-t0)/h (0<=d<=1). Then the interpolating function is:

  u(d)=y(t0+dh)=(1-d)y0 + dy1 + d(d-1)[(1-2d)(y1-y0) + (d-1)hf0 + dhf1]

Rearrange so we only have to through each array once:
  u(d)=cy0*y0 + cy1*y1 + cf0*f0 + cf1*f1
  cy0= 1 - d^2(3 - 2*d)
  cy1=     d^2(3 - 2*d) = (1-cy0)
  cf0=d*(d-1)^2*h = h*d*(d-1) * (d-1)
  cf1=d^2*(d-1)*h = h*d*(d-1) * d
  
Note: you can't get a 3rd order estimate of the derivative ft, only the 
state yt.

Cost is about 20 + 7*n flops (n is the vector length).

It is OK if yt is the same object as y0 or y1. */
/*static*/ void IntegratorRep::interpolateOrder3
   (const double& t0, const Vector& y0, const Vector& f0,
    const double& t1, const Vector& y1, const Vector& f1,
    const double& t, Vector& yt) {
    assert(t0 < t1);
    assert(t0 <= t && t <= t1);
    assert(f0.size()==y0.size() && y1.size()==y0.size() 
        && f1.size()==y0.size());

    const double h=t1-t0, d=(t-t0)/h, d1=d-1, hdd1=h*d*d1, dd32=d*d*(3-2*d);
    const Real cy1 = Real(dd32),   cy0 = Real(1-dd32);
    const Real cf1 = Real(hdd1*d), cf0 = Real(hdd1*d1);

    yt = cy0*y0 + cy1*y1 + cf0*f0 + cf1*f1; // + O(h^4)
}



//------------------------------------------------------------------------------
//                          ESTIMATE ROOT TIME
//------------------------------------------------------------------------------
/* We have bracketed a zero crossing for some function f(t) between (tLow,fLow)
and (tHigh,fHigh), not including the end points. That means tHigh > tLow, and
sign(fLow) != sign(fHigh) (sign(x) returns -1, 0, 1). We want to estimate 
time tRoot with tLow<tRoot<tHigh such that f(tRoot) is zero. For a nicely 
behaved continuous function this is just the secant method:
     x = fHigh/(fHigh-fLow)  (0 <= x <= 1)
     tRoot = tHigh - x*(tHigh-tLow)
However, if the function appears to be discontinuous we'll simply bisect the
interval (x==0.5 above). We decide it is discontinuous if either end point is
exactly zero, which would occur with a boolean function, for example. Also, 
if the time interval is already at or below the smallest allowable 
localization window, we'll bisect.

One further twist, taken from CVODES, is to allow the caller to provide a 
bias (> 0) which will bias our returned tRoot further into the lower (bias<1)
or upper (bias>1) half-interval. If bias==1 this is the pure secant method. 
Bias has no effect on functions deemed to be discontinuous; we'll always 
bisect the interval in that case.

Finally, we won't return tRoot very close to either end of the interval. 
Instead we define a buffer zone at either end of width 10% of the interval, 
and push tRoot away from the edges if it gets any closer. And in any case 
we'll require at least 1/2 minWindow from either end, even if 10% of the
interval is smaller than that.

Note that "minWindow" here is *not* the desired localization window based on
user accuracy requirements, but the (typically much smaller) smallest 
allowable localization window based on numerical roundoff considerations. */
/*static*/ double IntegratorRep::
estimateRootTime(double tLow, Real fLow, double tHigh, Real fHigh,
                 Real bias, double minWindow)
{
    assert(tLow < tHigh);
    assert(sign(fLow) != sign(fHigh));
    assert(bias > 0);
    assert(minWindow > 0);

    const double h = tHigh-tLow;

    if (fLow==0 || fHigh==0 || h <= minWindow) {
        // bisecting
        return tLow + h/2;
    }

    // Use secant method.

    const Real x = fHigh/(fHigh-bias*fLow);
    double tRoot = tHigh - double(x)*h;

    // If tRoot is too close to either end point we'll assume bad behavior
    // and guess a value 10% of the interval away from the end.

    const double BufferZone = std::max(0.1*h, minWindow/2);
    tRoot = std::max(tRoot, tLow  + BufferZone);
    tRoot = std::min(tRoot, tHigh - BufferZone);

    return tRoot;
}


//------------------------------------------------------------------------------
//                           CALC EVENT ORDER
//------------------------------------------------------------------------------

namespace {
// Local helper class to facilitate sort by estimated occurrence time.
class EventSorter {
public:
    EventSorter(int ord, ActiveWitnessIndex awx, double estEventTime) 
        : ordinal(ord), index(awx), estTime(estEventTime) { }

    // Order by estimated time, or in case of tie by active witness index.
    bool operator<(const EventSorter& r) const {
        if (estTime < r.estTime) return true;
        if (estTime > r.estTime) return false;
        return index < r.index;
    }

    unsigned            ordinal; // position in the triggeringWitnesses array
    ActiveWitnessIndex  index;   // position in the m_witnesses array
    double              estTime; // when we think it triggered
};
}

/*static*/ void IntegratorRep::
calcEventOrder(const Array_<ActiveWitnessIndex>&    triggeringWitnesses,
               const Array_<double>&                estEventTimes,
               Array_<unsigned>&                    witnessOrder)
{
    const unsigned n = triggeringWitnesses.size();
    assert(estEventTimes.size()==n);
    witnessOrder.resize(n);
    if (n==0) 
        return;

    if (n==1) {
        witnessOrder[0] = 0;
        return;
    }

    // otherwise sort
    Array_<EventSorter> events;
    for (unsigned i=0; i<n; ++i)
        events.emplace_back(i, triggeringWitnesses[i], estEventTimes[i]);
    std::sort(events.begin(), events.end());
    for (unsigned i=0; i<n; ++i) 
        witnessOrder[i] = events[i].ordinal;
}