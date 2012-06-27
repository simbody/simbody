/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
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

#include "SimTKcommon/internal/EventReporter.h"

// Workaround for a Microsoft compiler bug
#pragma optimize("g", off)

namespace SimTK {

EventReporter::~EventReporter() {
}

class TriggeredEventReporter::TriggeredEventReporterImpl {
public:
    TriggeredEventReporterImpl(Stage requiredStage) : requiredStage(requiredStage) {
    }
    EventTriggerInfo triggerInfo;
    Stage requiredStage;
};

TriggeredEventReporter::TriggeredEventReporter(Stage requiredStage) {
    impl = new TriggeredEventReporterImpl(requiredStage);
}

TriggeredEventReporter::TriggeredEventReporter(const TriggeredEventReporter& clone) {
    impl = new TriggeredEventReporterImpl(*clone.impl);
}

TriggeredEventReporter& TriggeredEventReporter::operator=(const TriggeredEventReporter& clone) {
    impl = new TriggeredEventReporterImpl(*clone.impl);
    return *this;
}

TriggeredEventReporter::~TriggeredEventReporter() {
    delete impl;
}

EventTriggerInfo& TriggeredEventReporter::getTriggerInfo() {
    return impl->triggerInfo;
}

Stage TriggeredEventReporter::getRequiredStage() const {
    return impl->requiredStage;
}

ScheduledEventReporter::~ScheduledEventReporter() {
}

class PeriodicEventReporter::PeriodicEventReporterImpl {
public:
    PeriodicEventReporterImpl(Real eventInterval) : eventInterval(eventInterval) {
        SimTK_APIARGCHECK1_ALWAYS(eventInterval > 0.0, "PeriodicEventReporterImpl", "PeriodicEventReporterImpl", "The interval was %d.  It must be > 0", eventInterval);
    }
    Real eventInterval;
};

PeriodicEventReporter::PeriodicEventReporter(Real eventInterval) {
    impl = new PeriodicEventReporterImpl(eventInterval);
}

PeriodicEventReporter::~PeriodicEventReporter() {
    delete impl;
}

Real PeriodicEventReporter::getNextEventTime(const State& state, bool includeCurrentTime) const {
    Real currentTime = state.getTime();
    long count = (long)std::floor(currentTime/impl->eventInterval);
    volatile Real eventTime = count*impl->eventInterval;
    while (eventTime < currentTime || (eventTime == currentTime && !includeCurrentTime)) {
        count++;
        eventTime = count*impl->eventInterval;
    }
    return eventTime;
}

Real PeriodicEventReporter::getEventInterval() const {
    return impl->eventInterval;
}

void PeriodicEventReporter::setEventInterval(Real eventInterval) {
    SimTK_APIARGCHECK1_ALWAYS(eventInterval > 0.0, "PeriodicEventReporter", "setEventInterval", "The interval was %d.  It must be > 0", eventInterval);
    impl->eventInterval = eventInterval;
}

} // namespace SimTK
