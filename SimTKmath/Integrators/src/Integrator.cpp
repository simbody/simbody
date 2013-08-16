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
 * This is the private (library side) implementation of the Simmath
 * Integrator family of classes.
 */

#include "SimTKcommon.h"
#include "simmath/Integrator.h"

#include "IntegratorRep.h"

#include <exception>
#include <limits>
#include <iostream>

using std::cout; using std::endl;

namespace SimTK {

    //////////////////////////////////
    // IMPLEMENTATION OF INTEGRATOR //
    //////////////////////////////////

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

Integrator::SuccessfulStepStatus 
Integrator::stepTo(Real reportTime, Real advanceLimit) {
    return updRep().stepTo(reportTime, advanceLimit);
}

Integrator::SuccessfulStepStatus 
Integrator::stepBy(Real interval, Real advanceIntervalLimit) {
    const Real t = getRep().getState().getTime();
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

const Array_<EventId>& 
Integrator::getTriggeredEvents() const {
    if (getRep().getStepCommunicationStatus() != IntegratorRep::StepHasBeenReturnedWithEvent) {
        SimTK_THROW2(CantAskForEventInfoWhenNoEventTriggered, "getTriggeredEvents",
                     getRep().getState().getTime());
        //NOTREACHED
    }
    return getRep().getTriggeredEvents();
}

const Array_<Real>&
Integrator::getEstimatedEventTimes() const {
    if (getRep().getStepCommunicationStatus() != IntegratorRep::StepHasBeenReturnedWithEvent) {
        SimTK_THROW2(CantAskForEventInfoWhenNoEventTriggered, "getEstimatedEventTimes",
                     getRep().getState().getTime());
        //NOTREACHED
    }
    return getRep().getEstimatedEventTimes();
}

const Array_<Event::Trigger>&
Integrator::getEventTransitionsSeen() const {
    if (getRep().getStepCommunicationStatus() != IntegratorRep::StepHasBeenReturnedWithEvent) {
        SimTK_THROW2(CantAskForEventInfoWhenNoEventTriggered, "getEventTransitionsSeen",
                     getRep().getState().getTime());
        //NOTREACHED
    }
    return getRep().getEventTransitionsSeen();
}

const State& Integrator::getState() const {
    return getRep().getState();
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

Real Integrator::getActualInitialStepSizeTaken() const {
    return getRep().getActualInitialStepSizeTaken();
}
Real Integrator::getPreviousStepSizeTaken() const {
    return getRep().getPreviousStepSizeTaken();
}
Real Integrator::getPredictedNextStepSize() const {
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

    //////////////////////////////////////
    // IMPLEMENTATION OF INTEGRATOR REP //
    //////////////////////////////////////

IntegratorRep::IntegratorRep
       (Integrator*               handle,
        const System&             system)
  : myHandle(handle), sys(system)
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
    systemHasTimeAdvancedEvents = getSystem().hasTimeAdvancedEvents();

    // Allocate integrator-local data structures now that we know the sizes.
    const int ny = getAdvancedState().getNY();
    const int nc = getAdvancedState().getNYErr();
    const int ne = getAdvancedState().getNEventTriggers();
    timeScaleInUse = getSystem().getDefaultTimeScale();

    // Set accuracy and consTol to their user-requested values or
    // to the appropriate defaults.
    setAccuracyAndTolerancesFromUserRequests();

    // Now we can ask the System for information about its event triggers.
    // (Event trigger info is Instance stage information.) Note that the event
    // localization windows here are expressed in terms of the system timescale
    // -- be sure to multiply by timescale before using.
    getSystem().calcEventTriggerInfo(getAdvancedState(), eventTriggerInfo);

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

    // Handler is allowed to throw an exception if it fails since we don't
    // have a way to recover.
    HandleEventsOptions handleOpts(getConstraintToleranceInUse());
    HandleEventsResults results;
    getSystem().handleEvents(updAdvancedState(),
                             Event::Cause::Initialization,
                             Array_<EventId>(),
                             handleOpts, results);
    SimTK_ERRCHK_ALWAYS(
        results.getExitStatus()!=HandleEventsResults::ShouldTerminate,
        "Integrator::initialize()", 
        "An initialization event handler requested termination.");

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

    // Record the continuous parts of this now-realized initial state as the 
    // previous state as well (previous state is used when we have to back up
    // from a failed step attempt).
    saveStateAndDerivsAsPrevious(getAdvancedState());

    // The initial state is set so it looks like we just *completed* a step to 
    // get here. That way if the first reportTime is zero, this will get 
    // reported.
    setStepCommunicationStatus(CompletedInternalStepNoEvent);
    startOfContinuousInterval = true;

    // This will be set to something meaningful when the simulation ends.
    terminationReason = Integrator::InvalidTerminationReason;

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

void IntegratorRep::reinitialize(Stage stage, bool shouldTerminate) {
    if (stage < Stage::Report) {
        startOfContinuousInterval = true;
        setUseInterpolatedState(false);
    }
    if (shouldTerminate) {
        setStepCommunicationStatus(FinalTimeHasBeenReturned);
        terminationReason = Integrator::EventHandlerRequestedTermination;
    }
    methodReinitialize(stage,shouldTerminate);
}

} // namespace SimTK


