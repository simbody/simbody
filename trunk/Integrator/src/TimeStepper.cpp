/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmath(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
 * Authors: Michael Sherman, Peter Eastman                                    *
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

void TimeStepper::stepTo(Real reportTime) {
    rep->stepTo(reportTime);
}

bool TimeStepper::getReportAllSignificantStates() {
    return rep->getReportAllSignificantStates();
}

void TimeStepper::setReportAllSignificantStates(bool b) {
    rep->setReportAllSignificantStates(b);
}


    ////////////////////////////////////////
    // IMPLEMENTATION OF TIME STEPPER REP //
    ////////////////////////////////////////

TimeStepperRep::TimeStepperRep(TimeStepper* handle, const System& system) : myHandle(handle), system(system), integ(0), reportAllSignificantStates(false) {
}

void TimeStepperRep::stepTo(Real time) {
    while (!integ->isSimulationOver()) {
        std::vector<EventId> scheduledEventIds;
        std::vector<EventId> scheduledReportIds;
        Real nextScheduledEvent = NTraits<Real>::getInfinity();
        Real nextScheduledReport = NTraits<Real>::getInfinity();
        Real currentTime = integ->getTime();
        system.calcTimeOfNextScheduledEvent(integ->getState(), nextScheduledEvent, scheduledEventIds, lastEventTime != currentTime);
        system.calcTimeOfNextScheduledReport(integ->getState(), nextScheduledReport, scheduledReportIds, lastReportTime != currentTime);
        Real reportTime = std::min(nextScheduledReport, time);
        Real eventTime = std::min(nextScheduledEvent, time);
        Integrator::SuccessfulStepStatus status = integ->stepTo(reportTime, eventTime);
        Stage lowestModified = Stage::Report;
        bool shouldTerminate;
        switch (status) {
            case Integrator::ReachedStepLimit: {
                if (reportAllSignificantStates)
                    return;
                continue;
            }
            case Integrator::StartOfContinuousInterval: {
                if (reportAllSignificantStates)
                    return;
                continue;
            }
            case Integrator::ReachedReportTime: {
                if (integ->getTime() >= nextScheduledReport) {
                    system.reportEvents(integ->getState(),
                        System::ScheduledEvents,
                        scheduledReportIds);
                    lastReportTime = integ->getTime();
                }
                if (integ->getTime() >= time || reportAllSignificantStates)
                    return;
                continue;
            }
            case Integrator::ReachedScheduledEvent: {
                system.handleEvents(integ->updAdvancedState(),
                    System::ScheduledEvents,
                    scheduledEventIds,
                    integ->getAccuracyInUse(),
                    integ->getStateWeightsInUse(),
                    integ->getConstraintWeightsInUse(),
                    lowestModified, shouldTerminate);
                lastEventTime = integ->getTime();
                break;
            }
            case Integrator::TimeHasAdvanced: {
                system.handleEvents(integ->updAdvancedState(),
                    System::TimeAdvancedEvent,
                    std::vector<EventId>(),
                    integ->getAccuracyInUse(),
                    integ->getStateWeightsInUse(),
                    integ->getConstraintWeightsInUse(),
                    lowestModified, shouldTerminate);
                break;
            }
            case Integrator::ReachedEventTrigger: {
                system.handleEvents(integ->updAdvancedState(),
                    System::TriggeredEvents,
                    integ->getTriggeredEvents(),
                    integ->getAccuracyInUse(),
                    integ->getStateWeightsInUse(),
                    integ->getConstraintWeightsInUse(),
                    lowestModified, shouldTerminate);
                break;
            }
            case Integrator::EndOfSimulation: {
                system.handleEvents(integ->updAdvancedState(),
                    System::TerminationEvent,
                    std::vector<EventId>(),
                    integ->getAccuracyInUse(),
                    integ->getStateWeightsInUse(),
                    integ->getConstraintWeightsInUse(),
                    lowestModified, shouldTerminate);
                break;
            }
            default: assert(!"Unrecognized return from stepTo()");
        }
        integ->reinitialize(lowestModified, shouldTerminate);
        if (reportAllSignificantStates)
            return;
    }
}

} // namespace SimTK



