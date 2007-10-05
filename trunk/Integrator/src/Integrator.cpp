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
 * This is the private (library side) implementation of the Simmath
 * Integrator family of classes.
 */

#include "SimTKcommon.h"
#include "simmath/Integrator.h"

#include "simmath/IntegratorRep.h"

#include <exception>
#include <limits>
#include <iostream>

using std::cout; using std::endl;

namespace SimTK {

    //////////////////////////////////
    // IMPLEMENTATION OF INTEGRATOR //
    //////////////////////////////////

void Integrator::resetAllStatistics() {
    updRep().resetIntegratorStatistics();
    updRep().resetMethodStatistics();
}

void Integrator::initialize(const State& initState) {
    updRep().updAdvancedState() = initState;
    getRep().getSystem().realize(getRep().getAdvancedState(), Stage::Model);
    updRep().initialize(getRep().getAdvancedState());
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
    const Real t = getRep().getAdvancedState().getTime();
    return updRep().stepTo(t + interval, t + advanceIntervalLimit);
}

bool Integrator::isSimulationOver() const {
    return getRep().isSimulationOver();
}

Integrator::TerminationReason Integrator::getTerminationReason() const {
    return getRep().getTerminationReason();
}

/*static*/ String Integrator::successfulStepStatusString(SuccessfulStepStatus stat) {
    switch(stat) { 
        case ReachedReportTime: return "ReachedReportTime";
        case ReachedEventTrigger: return "ReachedEventTrigger";
        case ReachedScheduledEvent: return "ReachedScheduledEvent";
        case TimeHasAdvanced: return "TimeHasAdvanced";
        case ReachedStepLimit: return "ReachedStepLimit";
        case EndOfSimulation: return "EndOfSimulation";
        case StartOfContinuousInterval: return "StartOfContinuousInterval";
        case InvalidSuccessfulStepStatus: return "InvalidSuccessfulStepStatus";
        default: return String("UNKNOWN SUCCESSFUL STEP STATUS ") + String((int)stat);
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

const Array<int>& 
Integrator::getTriggeredEvents() const {
    if (getRep().getStepCommunicationStatus() != IntegratorRep::StepHasBeenReturnedWithEvent) {
        SimTK_THROW2(CantAskForEventInfoWhenNoEventTriggered, "getTriggeredEvents",
                     getRep().getState().getTime());
        //NOTREACHED
    }
    return getRep().getTriggeredEvents();
}

const Array<Real>&
Integrator::getEstimatedEventTimes() const {
    if (getRep().getStepCommunicationStatus() != IntegratorRep::StepHasBeenReturnedWithEvent) {
        SimTK_THROW2(CantAskForEventInfoWhenNoEventTriggered, "getEstimatedEventTimes",
                     getRep().getState().getTime());
        //NOTREACHED
    }
    return getRep().getEstimatedEventTimes();
}

const Array<EventStatus::EventTrigger>&
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

Real Integrator::getTimeScaleInUse() const {
    return getRep().getTimeScaleInUse();
}

const Vector& Integrator::getStateWeightsInUse() const {
    return getRep().getStateWeightsInUse();
}

const Vector& Integrator::getConstraintWeightsInUse() const {
    return getRep().getConstraintWeightsInUse();
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
long Integrator::getNStepsAttempted() const {
    return getRep().getNStepsAttempted();
}
long Integrator::getNStepsTaken() const {
    return getRep().getNStepsTaken();
}
long Integrator::getNRealizations() const {
    return getRep().getNRealizations();
}
long Integrator::getNProjections() const {
    return getRep().getNProjections();
}
long Integrator::getNErrorTestFailures() const {
    return getRep().getNErrorTestFailures();
}
long Integrator::getNRealizationFailures() const {
    return getRep().getNRealizationFailures();
}
long Integrator::getNProjectionFailures() const {
    return getRep().getNProjectionFailures();
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
    assert(rep.userMinStepSize==-1. || z >= rep.userMinStepSize);
    assert(rep.userMaxStepSize==-1. || z <= rep.userMaxStepSize);
    updRep().userInitStepSize = z;
}
void Integrator::setMinimumStepSize(Real z) { 
    assert(z == -1. || z > 0.);
    assert(rep.userInitStepSize==-1. || z <= rep.userInitStepSize);
    assert(rep.userMaxStepSize ==-1. || z <= rep.userMaxStepSize);
    updRep().userMinStepSize = z;
}
void Integrator::setMaximumStepSize(Real z) {
    assert(z == -1. || z > 0.);
    assert(rep.userInitStepSize==-1. || z >= rep.userInitStepSize);
    assert(rep.userMinStepSize ==-1. || z >= rep.userMinStepSize);
    updRep().userMaxStepSize = z;
}
void Integrator::setAccuracy(Real accuracy) {
    assert(accuracy == -1. || (0. < accuracy && accuracy < 1.));
    updRep().userAccuracy = accuracy;
}
void Integrator::setRelativeTolerance(Real relTol) {
    assert(relTol  == -1. || (0. < relTol  && relTol  <= 1.));
    updRep().userRelTol=relTol;
}
void Integrator::setAbsoluteTolerance(Real absTol) {
    assert(absTol  == -1. || (0. < absTol  && absTol  <= 1.));
    updRep().userAbsTol=absTol;
}
void Integrator::setConstraintTolerance(Real consTol) {
    assert(consTol == -1. || (0. < consTol && consTol <= 1.));
    updRep().userConsTol=consTol;
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

bool Integrator::methodHasErrorControl() {
    return getRep().methodHasErrorControl();
}

const char* Integrator::getMethodName() {
    return getRep().getMethodName();
}

int Integrator::getMethodMinOrder() {
    return getRep().getMethodMinOrder();
}

int Integrator::getMethodMaxOrder() {
    return getRep().getMethodMaxOrder();
}

    //////////////////////////////////////
    // IMPLEMENTATION OF INTEGRATOR REP //
    //////////////////////////////////////

IntegratorRep::IntegratorRep
       (Integrator*               handle,
        const System&             system)
  : myHandle(handle), sys(system),
    stepCommunicationStatus(InvalidStepCommunicationStatus),
    nextStepSizeToTry(NaN), idealNextStepSize(NaN),
    tLow(NaN), tHigh(NaN),
    useInterpolatedState(false),
    tPrev(NaN)
{
    initializeUserStuff();
    resetIntegratorStatistics();
    resetMethodStatistics();
}


void IntegratorRep::initialize(const State& initState) {
  try
  { updAdvancedState() = initState;
     
    // Freeze problem dimensions.
    getSystem().realize(getAdvancedState(), Stage::Model);

    // Some DynamicSystems need to get control whenever time is advanced
    // successfully (and irreversibly) so that they can do discrete updates.
    systemHasTimeAdvancedEvents = getSystem().hasTimeAdvancedEvents();

    // Allocate data structures now that we know the sizes.
    const int ny = getSystem().getDefaultState().getNY();
    const int nc = getSystem().getDefaultState().getNYErr();
    const int ne = getSystem().getDefaultState().getNEvents();
    stateWeightsInUse.resize(ny);
    constraintWeightsInUse.resize(nc);

    // Freeze instance variables (parameters).
    getSystem().realize(getAdvancedState(), Stage::Instance);
    timeScaleInUse = getSystem().calcTimescale(getAdvancedState());

    // Set accuracy, consTol, relTol, absTol to their user-requested values or to
    // the appropriate defaults.
    setAccuracyAndTolerancesFromUserRequests();

    // Now we can ask the DynamicSystem for information about its event triggers.
    // (Event trigger info is Instance stage information.) Note that the event
    // localization windows here are expressed in terms of the system timescale -- 
    // be sure to multiply by timescale before using.
    getSystem().calcEventTriggerInfo(getAdvancedState(), eventTriggerInfo);

    // Obtain the constraint error tolerance units (actually 1/tolerance) for each constraint.
    // Tolerance units are Instance stage information; they cannot change with time. The
    // actual tolerances we'll use will be these unit tolerances scaled by the user's
    // accuracy request.
    getSystem().calcYErrUnitTolerances(getAdvancedState(), constraintWeightsInUse);

    // Obtain the state variable weights. Weights are Position-stage information but
    // are expected to remain constant over substantial intervals, so that we can
    // consider them to be constant during a time step. These should be recalculated
    // from time to time during a simulation, whenever "substantial" changes to the
    // configuration have been made.
    getSystem().realize(getAdvancedState(), Stage::Position);
    getSystem().calcYUnitWeights(getAdvancedState(), stateWeightsInUse);

    // Using the constraint unit tolerances and state weights we just calculated, project the
    // states to drive constraint errors below accuracy*unitTolerance for each constraint. 
    getSystem().realize(getAdvancedState(), Stage::Velocity); // all kinematics
    projectStateAndErrorEstimate(updAdvancedState(), Vector()); // no err est here

    // Refresh the weights in case anything significant was done by
    // the initial projection.
    getSystem().calcYUnitWeights(getAdvancedState(), stateWeightsInUse);

    realizeStateDerivatives(getAdvancedState());

    // Record the continuous parts of this now-realized initial state
    // as the previous state as well.
    saveStateAsPrevious(getAdvancedState());

    // The initial state is set so it looks like we just *completed* a step to 
    // get here. That way if the first reportTime is zero, this will get reported.
    setStepCommunicationStatus(CompletedInternalStepNoEvent);
    startOfContinuousInterval = true;
    terminationReason = Integrator::InvalidTerminationReason;

    // Call this virtual method in case the concrete class has additional 
    // integration needs. Note that it gets the State as we have adjusted it,
    // NOT the one originally passed to initialize().
    methodInitialize(getAdvancedState());

  } catch (const std::exception& e) {
    SimTK_THROW1(Integrator::InitializationFailed, e.what());
  } catch (...) {
    SimTK_THROW1(Integrator::InitializationFailed, "UNKNOWN EXCEPTION");
  }
}

void IntegratorRep::reinitialize(Stage stage, bool shouldTerminate) {
    if (stage < Stage::Report)
        startOfContinuousInterval = true;
    if (shouldTerminate) {
        setStepCommunicationStatus(FinalTimeHasBeenReturned);
        terminationReason = Integrator::EventHandlerRequestedTermination;
    }
    methodReinitialize(stage,shouldTerminate);
}

} // namespace SimTK


