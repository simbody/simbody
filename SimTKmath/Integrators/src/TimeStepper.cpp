/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman, Peter Eastman                                    *
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
 * TimeStepper family of classes.
 */

#include "SimTKcommon.h"
#include "simmath/TimeStepper.h"

#include "TimeStepperRep.h"

#include <exception>
#include <limits>

namespace SimTK {

    ////////////////////////////////////
    // IMPLEMENTATION OF TIME STEPPER //
    ////////////////////////////////////

TimeStepper::TimeStepper(const System& sys) {
    rep = new TimeStepperRep(this, sys);
}

TimeStepper::TimeStepper(const System& sys, Integrator& integrator) {
    rep = new TimeStepperRep(this, sys);
    setIntegrator(integrator);
}

TimeStepper::~TimeStepper() {
    if (rep && rep->myHandle==this)
        delete rep;
    rep = 0;
}

void TimeStepper::setIntegrator(Integrator& integrator) {
    rep->setIntegrator(integrator);
}

const Integrator& TimeStepper::getIntegrator() const {
    return rep->getIntegrator();
}

Integrator& TimeStepper::updIntegrator() {
    return rep->updIntegrator();
}

const State& TimeStepper::getState() const {
    return rep->getState();
}

void TimeStepper::initialize(const State& initState) {
    updIntegrator().initialize(initState);
    rep->lastEventTime = -Infinity;
    rep->lastReportTime = -Infinity;
}

Integrator::SuccessfulStepStatus TimeStepper::stepTo(Real reportTime) {
    return( rep->stepTo(reportTime) );
}

bool TimeStepper::getReportAllSignificantStates() const {
    return rep->getReportAllSignificantStates();
}

void TimeStepper::setReportAllSignificantStates(bool b) {
    rep->setReportAllSignificantStates(b);
}


    ////////////////////////////////////////
    // IMPLEMENTATION OF TIME STEPPER REP //
    ////////////////////////////////////////

TimeStepperRep::TimeStepperRep(TimeStepper* handle, const System& system) 
:   myHandle(handle), system(system), integ(0), 
    reportAllSignificantStates(false) {}

Integrator::SuccessfulStepStatus TimeStepperRep::stepTo(Real time) {
    // Handler is allowed to throw an exception if it fails since we don't
    // have a way to recover.
    HandleEventsOptions handleOpts(integ->getConstraintToleranceInUse());

    if (integ->isInfinityNormInUse())
        handleOpts.setOption(HandleEventsOptions::UseInfinityNorm);

    Array_<EventId> scheduledEventIds, scheduledReportIds;
    while (!integ->isSimulationOver()) {
        Real nextScheduledEvent  = Infinity;
        Real nextScheduledReport = Infinity;
        Real currentTime         = integ->getTime();
        system.realize(integ->getState(), Stage::Time);
        system.realize(integ->getAdvancedState(), Stage::Time);
        system.calcTimeOfNextScheduledEvent 
           (integ->getState(), nextScheduledEvent,  scheduledEventIds,  
            lastEventTime != currentTime);  // whether to allow now as an answer
        system.calcTimeOfNextScheduledReport
           (integ->getState(), nextScheduledReport, scheduledReportIds, 
            lastReportTime != currentTime); // whether to allow now as an answer

        Real reportTime = std::min(nextScheduledReport, time);
        Real eventTime  = std::min(nextScheduledEvent,  time);

        //---------------- take continuous step ----------------
        Integrator::SuccessfulStepStatus status = 
            integ->stepTo(reportTime, eventTime);
        //------------------------------------------------------

        Stage lowestModified = Stage::Report;
        bool shouldTerminate;
        switch (status) {
            case Integrator::ReachedStepLimit: {
                if (reportAllSignificantStates)
                    return status;
                continue;
            }
            case Integrator::StartOfContinuousInterval: {
                if (reportAllSignificantStates)
                    return status;
                continue;
            }
            case Integrator::ReachedReportTime: {
                if (integ->getTime() >= nextScheduledReport) {
                    system.reportEvents(integ->getState(),
                        Event::Cause::Scheduled,
                        scheduledReportIds);
                    lastReportTime = integ->getTime();
                }
                if (integ->getTime() >= time || reportAllSignificantStates)
                    return status;
                continue;
            }
            case Integrator::ReachedScheduledEvent: {
                HandleEventsResults results;
                system.handleEvents(integ->updAdvancedState(),
                                    Event::Cause::Scheduled,
                                    scheduledEventIds,
                                    handleOpts, results);
                lowestModified = results.getLowestModifiedStage();
                shouldTerminate = 
                    results.getExitStatus()==HandleEventsResults::ShouldTerminate;
                lastEventTime = integ->getTime();
                break;
            }
            case Integrator::TimeHasAdvanced: {
                HandleEventsResults results;
                system.handleEvents(integ->updAdvancedState(),
                                    Event::Cause::TimeAdvanced,
                                    Array_<EventId>(),
                                    handleOpts, results);
                lowestModified = results.getLowestModifiedStage();
                shouldTerminate = 
                    results.getExitStatus()==HandleEventsResults::ShouldTerminate;
                break;
            }
            case Integrator::ReachedEventTrigger: {
                HandleEventsResults results;
                system.handleEvents(integ->updAdvancedState(),
                                    Event::Cause::Triggered,
                                    integ->getTriggeredEvents(),
                                    handleOpts, results);
                lowestModified = results.getLowestModifiedStage();
                shouldTerminate = 
                    results.getExitStatus()==HandleEventsResults::ShouldTerminate;
                break;
            }
            case Integrator::EndOfSimulation: {
                HandleEventsResults results;
                system.handleEvents(integ->updAdvancedState(),
                                    Event::Cause::Termination,
                                    Array_<EventId>(),
                                    handleOpts, results);
                lowestModified = results.getLowestModifiedStage();
                shouldTerminate = 
                    results.getExitStatus()==HandleEventsResults::ShouldTerminate;
                break;
            }
            default: assert(!"Unrecognized return from stepTo()");
        }
        integ->reinitialize(lowestModified, shouldTerminate);
        if (reportAllSignificantStates)
            return status;
    }
    return Integrator::EndOfSimulation;
}

} // namespace SimTK



