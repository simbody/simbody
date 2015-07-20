/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-15 Stanford University and the Authors.        *
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

using namespace SimTK;

//==============================================================================
//                               TIME STEPPER
//==============================================================================

TimeStepper::TimeStepper() {
    rep = new TimeStepperRep(this);
}

TimeStepper::TimeStepper(Integrator& integrator) {
    rep = new TimeStepperRep(this);
    setIntegrator(integrator);
}

TimeStepper::~TimeStepper() {
    if (rep && rep->m_myHandle==this)
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
    rep->m_lastEventTime = -Infinity;
    rep->m_lastReportTime = -Infinity;
}

Integrator::SuccessfulStepStatus TimeStepper::stepTo(double reportTime) {
    return( rep->stepTo(reportTime) );
}

bool TimeStepper::getReportAllSignificantStates() const {
    return rep->getReportAllSignificantStates();
}

void TimeStepper::setReportAllSignificantStates(bool b) {
    rep->setReportAllSignificantStates(b);
}


//==============================================================================
//                             TIME STEPPER REP
//==============================================================================

Integrator::SuccessfulStepStatus TimeStepperRep::
stepTo(double requestedTime) {
    const System& sys = getSystem();

    // These will be re-used for each step. We're allocating them here to
    // avoid repeated heap allocations.
    EventTriggers   scheduledReportTimers, scheduledChangeTimers;
    EventsAndCauses triggeredEvents;
    Array_<EventId> ignoredEventIds;

    // The TimeAdvanced trigger list is always the same.
    const EventTriggers timeAdvancedTrigger{&sys.getTimeAdvancedTrigger()};

    while (!m_integ->isSimulationOver()) {
        double timeOfNextReport=Infinity, timeOfNextChange=Infinity;

        // Clear the lists of triggered events (otherwise they will get
        // appended to by noteEventOccurrences()).
        triggeredEvents.clear(); ignoredEventIds.clear();

        sys.realize(m_integ->getState(), Stage::Time);
        sys.realize(m_integ->getAdvancedState(), Stage::Time);

        /* The same Timer may appear both on the report and change lists since
        an Event may have both report and change actions. In that case, the
        integrator will return ReachedReportTime first, at which point 
        report actions should be performed, then on reentry will immediately
        return ReachedScheduledEvent at the same time.
        TimeAdvanced events behave similarly, but witness-triggered events
        get reported *before* the event, then the state that actually caused
        the event must not appear as part of the trajectory, so its report
        actions must be invoked after its change actions have fixed whatever
        is wrong. */

        //TODO: shouldn't we guarantee realization to Acceleration stage here?
        sys.findNextScheduledEventTimes
           (*m_integ, m_lastReportTime, m_lastEventTime,
            timeOfNextReport, scheduledReportTimers,
            timeOfNextChange, scheduledChangeTimers);

        const double reportTime = std::min(timeOfNextReport, requestedTime);
        const double eventTime  = std::min(timeOfNextChange, requestedTime);

        //---------------- take continuous step ----------------
        Integrator::SuccessfulStepStatus status = 
            m_integ->stepTo(reportTime, eventTime);
        //------------------------------------------------------

        EventChangeResult result;
        switch (status) {
            case Integrator::ReachedStepLimit: // not an event
                if (m_reportAllSignificantStates)
                    return status;
                continue;

            case Integrator::StartOfContinuousInterval: // not an event
                if (m_reportAllSignificantStates)
                    return status;
                continue;

            case Integrator::EndOfSimulation:
                continue; // we'll fall out of the loop next

            case Integrator::ReachedReportTime: {
                if (m_integ->getTime() >= timeOfNextReport) {
                    m_lastReportTime = m_integ->getTime();
                    sys.noteEventOccurrence(scheduledReportTimers,
                                            triggeredEvents, ignoredEventIds);
                    sys.performEventReportActions(*m_integ, triggeredEvents);
                }
                if (   m_integ->getTime() >= requestedTime 
                    || m_reportAllSignificantStates)
                    return status;
                continue;
            }

            // The remaining cases may cause state changes.

            case Integrator::ReachedScheduledEvent: {
                m_lastEventTime = m_integ->getTime();
                sys.noteEventOccurrence(scheduledChangeTimers,
                                        triggeredEvents, ignoredEventIds);
                sys.performEventChangeActions(*m_integ, triggeredEvents,
                                              result);
                break;
            }
            case Integrator::TimeHasAdvanced: {
                sys.noteEventOccurrence(timeAdvancedTrigger,
                                        triggeredEvents, ignoredEventIds);
                // Report the state we reached, just before performing any
                // TimeAdvanced state-change actions.
                sys.performEventReportActions(*m_integ, triggeredEvents);
                // Now modify the state, with time unchanged.
                sys.performEventChangeActions(*m_integ, triggeredEvents,
                                              result);
                break;
            }
            case Integrator::ReachedEventTrigger: {
                sys.noteEventOccurrence(m_integ->getTriggeredWitnesses(),
                                        triggeredEvents, ignoredEventIds);
                // First fix the bad state containing an unresolved issue.
                sys.performEventChangeActions(*m_integ, triggeredEvents,
                                              result);
                // Then report the repaired state.
                sys.performEventReportActions(*m_integ, triggeredEvents);
                break;
            }

            default: SimTK_ASSERT1_ALWAYS(!"Unrecognized return from stepTo()",
            "TimeStepper::stepTo(): integrator returned unknown status %d.\n",
            (int)status);
        }

        // If we get here we just returned from an event change action with
        // status returned in "results".

        // If the change action did nothing, just continue with the next step.
        if (   result.getExitStatus() == EventChangeResult::Succeeded
            && result.getLowestModifiedStage() == Stage::Infinity)
            continue;

        // If an Event handler fails, we're dead.
        if (result.getExitStatus() == EventChangeResult::Failed) {
            m_integ->terminate(Integrator::AnUnrecoverableErrorOccurred);
            SimTK_THROW2(TimeStepper::EventChangeActionFailed,
                         m_integ->getAdvancedTime(),
                         result.getMessage().c_str());
        }

        // Reinitialize the integrator as needed due to the state change.
        const bool shouldTerminate = 
            result.getExitStatus() == EventChangeResult::ShouldTerminate;
        const Stage lowestModified = result.getLowestModifiedStage();
        m_integ->reinitialize(lowestModified, shouldTerminate);

        if (m_reportAllSignificantStates)
            return status;
    }

    return Integrator::EndOfSimulation;
}





