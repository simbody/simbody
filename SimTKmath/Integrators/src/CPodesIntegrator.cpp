/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
 * Contributors: Michael Sherman                                              *
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

#include "simmath/CPodesIntegrator.h"

#include "IntegratorRep.h"
#include "CPodesIntegratorRep.h"

using namespace SimTK;


//------------------------------------------------------------------------------
//                            CPODES INTEGRATOR
//------------------------------------------------------------------------------

CPodesIntegrator::CPodesIntegrator(const System& sys, CPodes::LinearMultistepMethod method) {
    rep = new CPodesIntegratorRep(this, sys, method);
}

CPodesIntegrator::CPodesIntegrator(const System& sys, CPodes::LinearMultistepMethod method, CPodes::NonlinearSystemIterationType iterationType)  {
    rep = new CPodesIntegratorRep(this, sys, method, iterationType);
}

void CPodesIntegrator::setUseCPodesProjection() {
    CPodesIntegratorRep& cprep = dynamic_cast<CPodesIntegratorRep&>(*rep);
    cprep.setUseCPodesProjection();
}

void CPodesIntegrator::setOrderLimit(int order) {
    CPodesIntegratorRep& cprep = dynamic_cast<CPodesIntegratorRep&>(*rep);
    cprep.setOrderLimit(order);
}



//------------------------------------------------------------------------------
//                          CPODES INTEGRATOR REP
//------------------------------------------------------------------------------
// This class implements the abstract CPodesSystem interface understood by our 
// C++ interface to CPodes.

class CPodesIntegratorRep::CPodesSystemImpl : public CPodesSystem {
public:
    CPodesSystemImpl(CPodesIntegratorRep& integ, const System& system) 
    :   integ(integ), system(system) {}

    // Calculate ydot = f(t,y).
    int explicitODE(Real t, const Vector& y, Vector& ydot) const {
        try { 
            integ.setAdvancedStateAndRealizeDerivatives(t,y);
        }
        catch(...) { return CPodes::RecoverableError; } // assume recoverable
        ydot = integ.getAdvancedState().getYDot();
        return CPodes::Success;
    }

    // Calculate yerr = c(t,y).
    int constraint(Real t, const Vector& y, Vector& yerr) const {
        try { 
            integ.setAdvancedStateAndRealizeKinematics(t,y);
        }
        catch(...) { return CPodes::RecoverableError; } // assume recoverable
        yerr = integ.getAdvancedState().getYErr();
        return CPodes::Success;
    }

    // Given a state (t,y) not on the constraint manifold, return ycorr
    // such that (t,y+ycorr+eps) is on the manifold, with 
    // ||eps||_wrms <= epsProj. 'err' passed in as the integrator's current 
    // error estimate for state y; optionally project it to eliminate the 
    // portion normal to the manifold.
    int project(Real t, const Vector& y, 
                Vector& ycorr, Real epsProj, Vector& err) const {
        integ.setAdvancedState(t,y);
        State& advanced = integ.updAdvancedState();
       
        try {
            system.realize(advanced, Stage::Time);
            system.prescribeQ(advanced); // set q_p
            system.realize(advanced, Stage::Position);
            bool anyChanges;
            if (!integ.localProjectQAndQErrEstNoThrow
                    (advanced, err, anyChanges, Infinity)) //TODO: proj limit?
                return CPodes::RecoverableError;

            system.prescribeU(advanced); // set u_p
            system.realize(advanced, Stage::Velocity);
            if (!integ.localProjectUAndUErrEstNoThrow
                    (advanced, err, anyChanges, Infinity)) //TODO: proj limit?
                return CPodes::RecoverableError;
        }
        catch (...) { return CPodes::RecoverableError; } // assume recoverable
        ycorr = advanced.getY()-y;
        return CPodes::Success;
    }
    
    /**
     * Calculate the event trigger functions.
     */
    int root(Real t, const Vector& y, const Vector& yp, Vector& gout) const {
        try { 
            integ.setAdvancedStateAndRealizeDerivatives(t,y);
        }
        catch(...) { return CPodes::RecoverableError; } // assume recoverable
        integ.calcWitnessValues(integ.getAdvancedState(), gout);
        return CPodes::Success;
    }
private:
    CPodesIntegratorRep& integ;
    const System& system;
};

void CPodesIntegratorRep::init
   (CPodes::LinearMultistepMethod method, 
    CPodes::NonlinearSystemIterationType iterationType) 
{
    cpodes = new CPodes(CPodes::ExplicitODE, method, iterationType);
    cps = new CPodesSystemImpl(*this, getSystem());
    initialized = false;
    useCpodesProjection = false;
}

CPodesIntegratorRep::CPodesIntegratorRep
   (Integrator* handle, const System& sys, 
    CPodes::LinearMultistepMethod method) 
:   IntegratorRep(handle, sys), method(method) {
    init(method, method == CPodes::Adams ? CPodes::Functional : CPodes::Newton);
}

CPodesIntegratorRep::CPodesIntegratorRep
   (Integrator* handle, const System& sys, 
    CPodes::LinearMultistepMethod method, 
    CPodes::NonlinearSystemIterationType iterationType)
:   IntegratorRep(handle, sys), method(method) {
    init(method, iterationType);
}

CPodesIntegratorRep::~CPodesIntegratorRep() {
    delete cpodes;
    delete cps;
}

void CPodesIntegratorRep::methodInitialize(const State& state) {
    if (state.getSystemStage() < Stage::Model)
        reconstructForNewModel();
    initializeIntegrationParameters();
    initialized = true;
    pendingReturnCode = -1;
    previousStartTime = 0.0;
    resetMethodStatistics();
    getSystem().realize(state, Stage::Velocity);
    const int ny = state.getY().size();
    const int nc = state.getNYErr();
    Vector ydot(ny);
    if (cps->explicitODE(state.getTime(), Vector(state.getY()), ydot) 
        != CPodes::Success) 
    {
        SimTK_THROW1(Integrator::InitializationFailed, 
                     "Failed to calculate ydot");
    }
    int retval;
    //TODO: change this to do abstol only for q, reltol for u&z
    Real relTol = getAccuracyInUse();
    Real absTol = relTol/10; //TODO: base on weights
    if ((retval=cpodes->init(*cps, state.getTime(), 
                             Vector(state.getY()), ydot, 
                             CPodes::ScalarScalar, relTol, &absTol)) 
        != CPodes::Success) 
    {
        printf("init() returned %d\n", retval);
        SimTK_THROW1(Integrator::InitializationFailed, "init() failed");
    }
    cpodes->lapackDense(ny);
    cpodes->setNonlinConvCoef(Real(0.01)); // TODO (default is 0.1)
    if (useCpodesProjection) {
        const int nqerr = state.getNQErr(), nuerr = state.getNUErr();
        const Real tol = getConstraintToleranceInUse();
        Vector constraintTols(nqerr+nuerr);
        constraintTols(0,nqerr) = tol*state.getQErrWeights();
        constraintTols(nqerr,nuerr) = tol*state.getUErrWeights();
        cpodes->projInit(CPodes::L2Norm, CPodes::Nonlinear, 
                         constraintTols);
        cpodes->lapackDenseProj(nc, ny, CPodes::ProjectWithQRPivot);
    }
    else {
        cpodes->projDefine();
    }

    const int nWitnesses = (int)getWitnesses().size();
    cpodes->rootInit(nWitnesses);
    if (nWitnesses) {
        Array_<int> rootDir(nWitnesses);
        for (ActiveWitnessIndex awx(0); awx < nWitnesses; ++awx) {
            auto& w = *getWitnesses()[awx];
            if (w.getTriggerOnFallingSignTransition()) {
                if (w.getTriggerOnRisingSignTransition())
                    rootDir[awx] = 0; // All transitions
                else 
                    rootDir[awx] = -1; // Falling transitions only
            } else
                rootDir[awx] = 1; // Rising transitions only
        }
        cpodes->setRootDirection(rootDir);
    }
}

void CPodesIntegratorRep::methodReinitialize
   (Stage stage, bool shouldTerminate) {
    if (stage < Stage::Report) {
        pendingReturnCode = -1;
        State state = getAdvancedState();
        getSystem().realize(state, Stage::Acceleration);
        //TODO: change this to do abstol only for q, reltol for u&z
        Real relTol = getAccuracyInUse();
        Real absTol = relTol/10; //TODO: base on weights
        cpodes->reInit(*cps, state.getTime(), 
                       Vector(state.getY()), Vector(state.getYDot()), 
                       CPodes::ScalarScalar, relTol, &absTol);
    }
}

void CPodesIntegratorRep::initializeIntegrationParameters() {
    if (userInitStepSize != -1)
        cpodes->setInitStep(userInitStepSize);
    if (userMinStepSize != -1)
        cpodes->setMinStep(userMinStepSize);
    if (userMaxStepSize != -1)
        cpodes->setMaxStep(userMaxStepSize);
    if (userFinalTime != -1.) 
        cpodes->setStopTime(userFinalTime);
    if (userInternalStepLimit != -1) 
        cpodes->setMaxNumSteps(userInternalStepLimit);
    if (userProjectEveryStep != -1)
        if (userProjectEveryStep==1)
            cpodes->setProjFrequency(1); // every step
}

void CPodesIntegratorRep::reconstructForNewModel() {
    initialized = false;
    delete cpodes;
    cpodes = new CPodes(CPodes::ExplicitODE, CPodes::BDF, CPodes::Newton);
}

// Create an interpolated state at time t, which is between tPrev and tCurrent.
// If we haven't yet delivered an interpolated state in this interval, we have
// to initialize its discrete part from the advanced state.
void CPodesIntegratorRep::createInterpolatedState(Real t) {
    const System& system  = getSystem();
    const State& advanced = getAdvancedState();
    State&       interp   = updInterpolatedState();
    interp = advanced; // pick up discrete stuff.
    Vector yout(advanced.getY().size());
    cpodes->getDky(t, 0, yout);
    interp.updY() = yout;
    interp.updTime() = t;

    if (userProjectInterpolatedStates == 0) {
        system.realize(interp, Stage::Time);
        system.prescribeQ(interp);
        system.realize(interp, Stage::Position);
        system.prescribeU(interp);
        system.realize(interp, Stage::Velocity);
        return;
    }

    // We may need to project onto constraint manifold. Allow project()
    // to throw an exception if it fails since there is no way to recover here.
    realizeAndProjectKinematicsWithThrow(interp, ProjectOptions::LocalOnly);
}

// Take a step. See AbstractIntegratorRep::stepTo() for how this is supposed
// to behave. We have to go through some contortions to squeeze CPodes into
// that mold.
Integrator::SuccessfulStepStatus CPodesIntegratorRep::
stepTo(Real reportTime, Real scheduledEventTime) {
    assert(initialized);
    assert(reportTime >= getState().getTime());
    assert(scheduledEventTime >= getState().getTime());

    if (getStepCommunicationStatus() == FinalTimeHasBeenReturned) {
        SimTK_ERRCHK2_ALWAYS(!"EndOfSimulation already returned",
            "Integrator::stepTo()",
            "Attempted stepTo(t=%g) but final time %g had already been "
            "reached and returned."
            "\nCheck for Integrator::EndOfSimulation status, or use the "
            "Integrator::initialize() method to restart.",
            reportTime, userFinalTime);
    }
    
    // If this is the start of a continuous interval, return immediately so
    // the current state will be seen as part of the trajectory.
    if (isStartOfContinuousInterval()) {
        // The set of event witnesses might have changed.
        getSystem().findActiveEventWitnesses(getOwnerHandle(), updWitnesses());
        pendingReturnCode = -1; // forget post-event state
        setIsStartOfContinuousInterval(false);
        return Integrator::StartOfContinuousInterval;
    }

    CPodes::StepMode mode;
    if (userFinalTime != -1 || userAllowInterpolation==0)
        mode = (userReturnEveryInternalStep == 1) ? CPodes::OneStepTstop
                                                  : CPodes::NormalTstop;
    else
        mode = (userReturnEveryInternalStep == 1) ? CPodes::OneStep
                                                  : CPodes::Normal;

    // Keep taking steps until something interesting happens at tMax or
    // earlier.
    Real tMax = std::min(reportTime, scheduledEventTime);

    // We might have to create a fake tstop to prevent interpolation.
    bool isFakeTstop = false;
    if (userAllowInterpolation==0) {
        if (userFinalTime == -1 || userFinalTime > tMax) {
            isFakeTstop = true;
            cpodes->setStopTime(tMax);
        }
    }

    // Assume we'll return at the advanced state; we'll change this below
    // if necessary.
    setUseInterpolatedState(false);
    while (true) {
        Real tret;
        int res;
        const bool usePendingReturnCode = (pendingReturnCode != -1);
        if (usePendingReturnCode) {            
            // The last time returned was an event or report time. The 
            // integrator has already gone beyond that time, so reset 
            // everything to how it was after the last call to cpodes->step()
            // and then process the step.
            
            res = pendingReturnCode;
            tret = previousTimeReturned;
            if (savedY.size() > 0) { 
                setAdvancedStateAndRealizeKinematics(tret, savedY);
            } else {
                updAdvancedState().updTime() = tret;
            }
            pendingReturnCode = -1;
        }
        else if (tMax == getState().getTime()) {
            
            // A report or event is scheduled for the current time, so return 
            // immediately.
            
            res = CPodes::Success;
            tret = tMax;
            previousTimeReturned = tret;
            updAdvancedState().updTime() = tret;
        }
        else {
            // We're going to advance time now.
    
            // Auto-update discrete variables. This update is not allowed to 
            // affect any computations performed at the current state value so
            // does not invalidate any stage. Swap the discrete state update 
            // cache entries with the state variables.
            updAdvancedState().autoUpdateDiscreteVariables();

            previousStartTime = getAdvancedTime();
            Vector yout(getAdvancedState().getY().size());
            Vector ypout(getAdvancedState().getY().size()); // ignored
            int oldSteps=0, oldTestFailures=0, oldNonlinIterations=0, 
                oldNonlinConvFailures=0;
            cpodes->getNumSteps(&oldSteps);
            cpodes->getNumErrTestFails(&oldTestFailures);
            cpodes->getNumNonlinSolvIters(&oldNonlinIterations);
            cpodes->getNumNonlinSolvConvFails(&oldNonlinConvFailures);

            //---------------------step------------------------
            res = cpodes->step(tMax, &tret, yout, ypout, mode);
            //-------------------------------------------------
            if (res == CPodes::TstopReturn && isFakeTstop)
                res = CPodes::Success;

            if (res == CPodes::TooClose) {              
                // This happens when the user asked the integrator to advance 
                // time by a tiny amount, comparable to numerical precision.
                // Since CPODES cannot advance time by such small increments, 
                // and the state would not change significantly in that time 
                // anyway, just set the time while leaving the rest of the 
                // state unchanged.             
                tret = tMax;
                yout = getAdvancedState().getY();
                res = CPodes::Success;
            }

            int newSteps=0, newTestFailures=0, newNonlinIterations=0, 
                newNonlinConvFailures=0;
            cpodes->getNumSteps(&newSteps);
            cpodes->getNumErrTestFails(&newTestFailures);
            cpodes->getNumNonlinSolvIters(&newNonlinIterations);
            cpodes->getNumNonlinSolvConvFails(&newNonlinConvFailures);
            statsStepsTaken += newSteps-oldSteps;
            statsErrorTestFailures += newTestFailures-oldTestFailures;
            // Project stats were already updated in project() above.
            statsIterations += newNonlinIterations-oldNonlinIterations;
            statsConvergenceTestFailures += newNonlinConvFailures-oldNonlinConvFailures;
 
            // This takes care of prescribed motion.
            setAdvancedStateAndRealizeKinematics(tret, yout);
            previousTimeReturned = tret;
        }

        realizeStateDerivatives(getAdvancedState());
        
        // Check for integration errors.        
        if (res == CPodes::TooMuchWork) {         
            // The maximum number of steps was reached.          
            setStepCommunicationStatus(IntegratorRep::StepHasBeenReturnedNoEvent);
            return Integrator::ReachedStepLimit;
        }
        if (res < 0) {           
            // An error of some sort occurred.           
            SimTK_THROW2(Integrator::StepFailed, getAdvancedState().getTime(), 
                         "CPodes::step() returned an error");
        }

        // No error occurred.      

        // If a triggered event was isolated to the window (tLo,tHi], 
        // CPodes will have returned with tret==tHi, which is where the
        // advancedState is now. We need to return an interpolated, "last good"
        // state at tLo, with return status ReachedEventTrigger. The calling
        // TimeStepper will then invoke the event
        // handler to fix up the advanced state at tHi without reporting it.
        // (If the event handler modified the state, startOfContinuousInterval
        // will be set.) Then we will get called here again and need to report
        // the now-fixed advanced state at tHi as part of the trajectory,
        // either as an ordinary time advanced state or as a 
        // StartOfContinuousInterval state if something happened.
        // But there might be some reports due prior to tLo first.

        if (res == CPodes::RootReturn) {
            Real tLo, tHi;
            cpodes->getRootWindow(&tLo, &tHi);
            tret = tLo;
        }
        
        // Determine the correct return code.
        
        if (tret >= reportTime && reportTime <= scheduledEventTime) {          
            // We reached the scheduled report time.  
            // If necessary, generate an interpolated state.      
            if (tret > tMax) {
                setUseInterpolatedState(true);
                createInterpolatedState(tMax);
                realizeStateDerivatives(getInterpolatedState());
            }
            savedY.resize(0);
            pendingReturnCode = res;
            setStepCommunicationStatus(IntegratorRep::StepHasBeenReturnedNoEvent);
            return Integrator::ReachedReportTime;
        }

        if (tret >= scheduledEventTime) {           
            // We reached a scheduled event time.            
            savedY.resize(0);
            if (tret > scheduledEventTime) {              
                // Back up the advanced state to the event time.               
                savedY = getAdvancedState().getY();
                createInterpolatedState(scheduledEventTime);
                setAdvancedStateAndRealizeDerivatives(scheduledEventTime,
                                              getInterpolatedState().getY());
            }
            pendingReturnCode = res;
            setStepCommunicationStatus(IntegratorRep::StepHasBeenReturnedWithEvent);
            return Integrator::ReachedScheduledEvent;
        }


        if (res == CPodes::RootReturn) {           
            // An event was triggered in the interval (tLo,tHi]. We're 
            // going to return with an interpolated state at tLo; CPodes has
            // already capped the advanced state at tHi.
            
            ActiveWitnessSubset witnessIndices;
            Array_<double> eventTimes;
            Array_<Event::TriggerDirection> eventTransitions;
            const unsigned nWitnesses = getWitnesses().size();
            int* eventFlags = new int[nWitnesses];
            cpodes->getRootInfo(eventFlags);
            for (ActiveWitnessIndex i(0); i < nWitnesses; ++i)
                if (eventFlags[i] != 0) {
                    witnessIndices.push_back(i);
                    eventTimes.push_back(previousTimeReturned);
                    eventTransitions.push_back(eventFlags[i] == 1 
                        ? Event::Rising : Event::Falling);
                }
            delete[] eventFlags;

            // Generate an interpolated state at tLo.      
            setUseInterpolatedState(true);
            createInterpolatedState(tret);
            realizeStateDerivatives(getInterpolatedState());

            setTriggeredEvents(tret, previousTimeReturned, 
                               witnessIndices, eventTimes, eventTransitions);

            // For next time, we'll treat the state at tHi as an ordinary
            // trajectory step, but this will only get used if the event
            // handler makes no changes. Otherwise, we'll be called with
            // startOfContinuousInterval==true and the handler-modified
            // state will get reported above instead.
            pendingReturnCode = CPodes::Success;
            setStepCommunicationStatus(IntegratorRep::StepHasBeenReturnedWithEvent);
            return Integrator::ReachedEventTrigger;
        }

        if (res == CPodes::TstopReturn) {
            // The specified final time was reached.  
            if (usePendingReturnCode) {
                // The final step result was already reported; now just return
                // the same state but with an "end of simulation" status. No
                // further calls should be made to stepTo().
                setStepCommunicationStatus(IntegratorRep::FinalTimeHasBeenReturned);
                setTerminationReason(Integrator::ReachedFinalTime);
                return Integrator::EndOfSimulation;
            }
            // This is our first encounter with the final time. Report this as
            // a TimeHasAdvanced state if the user has asked for those, 
            // otherwise pretend there was a report scheduled there. This should
            // match the behavior of AbstractIntegratorRep.
            savedY.resize(0);
            pendingReturnCode = res;
            setStepCommunicationStatus(IntegratorRep::StepHasBeenReturnedNoEvent);
            return (userReturnEveryInternalStep == 1)
                ? Integrator::TimeHasAdvanced
                : Integrator::ReachedReportTime;
        }

        if (userReturnEveryInternalStep == 1) {           
            // The user asked to be notified of every internal step.           
            setStepCommunicationStatus(IntegratorRep::StepHasBeenReturnedNoEvent);
            return Integrator::TimeHasAdvanced;
        }

        // Otherwise keep going -- no one wants to see this step.
    }
}

Real CPodesIntegratorRep::getActualInitialStepSizeTaken() const {
    assert(initialized);
    Real size;
    cpodes->getActualInitStep(&size);
    return size;
}

Real CPodesIntegratorRep::getPreviousStepSizeTaken() const {
    assert(initialized);
    Real size;
    cpodes->getLastStep(&size);
    return size;
}

Real CPodesIntegratorRep::getPredictedNextStepSize() const {
    assert(initialized);
    Real size;
    cpodes->getCurrentStep(&size);
    return size;
}

int CPodesIntegratorRep::getNumStepsAttempted() const {
    assert(initialized);
    return statsStepsTaken+statsErrorTestFailures+statsConvergenceTestFailures;
}

int CPodesIntegratorRep::getNumStepsTaken() const {
    assert(initialized);
    return statsStepsTaken;
}

int CPodesIntegratorRep::getNumErrorTestFailures() const {
    assert(initialized);
    return statsErrorTestFailures;
}

int CPodesIntegratorRep::getNumConvergenceTestFailures() const {
    assert(initialized);
    return statsConvergenceTestFailures;
}

int CPodesIntegratorRep::getNumIterations() const {
    assert(initialized);
    return statsIterations;
}

void CPodesIntegratorRep::resetMethodStatistics() {
    statsStepsTaken = 0;
    statsErrorTestFailures = 0;
    statsConvergenceTestFailures = 0;
    statsIterations = 0;
}

const char* CPodesIntegratorRep::getMethodName() const {
    return (method == CPodes::BDF ? "CPodesBDF" : "CPodesAdams");
}

int CPodesIntegratorRep::getMethodMinOrder() const {
    return 1;
}

int CPodesIntegratorRep::getMethodMaxOrder() const {
    return (method == CPodes::BDF ? 5 : 12);
}

bool CPodesIntegratorRep::methodHasErrorControl() const {
    return true;
}

void CPodesIntegratorRep::setUseCPodesProjection() {
    SimTK_APIARGCHECK_ALWAYS(!initialized, "CPodesIntegrator", 
        "setUseCPodesProjection",
        "This method may not be invoked after the integrator has been initialized.");
    useCpodesProjection = true;
}

void CPodesIntegratorRep::setOrderLimit(int order) {
    cpodes->setMaxOrd(order);
}


