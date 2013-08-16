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

#include "SimTKcommon/internal/EventHandler.h"

// Workaround for a Microsoft compiler bug
#pragma optimize("g", off)

namespace SimTK {

EventHandler::~EventHandler() {
}

class TriggeredEventHandler::TriggeredEventHandlerImpl {
public:
    TriggeredEventHandlerImpl(Stage requiredStage) : requiredStage(requiredStage) {
    }
    EventTriggerInfo triggerInfo;
    Stage requiredStage;
};

TriggeredEventHandler::TriggeredEventHandler(Stage requiredStage) {
    impl = new TriggeredEventHandlerImpl(requiredStage);
}

TriggeredEventHandler::TriggeredEventHandler(const TriggeredEventHandler& clone) {
    impl = new TriggeredEventHandlerImpl(*clone.impl);
}

TriggeredEventHandler& TriggeredEventHandler::operator=(const TriggeredEventHandler& clone) {
    impl = new TriggeredEventHandlerImpl(*clone.impl);
    return *this;
}

TriggeredEventHandler::~TriggeredEventHandler() {
    delete impl;
}

EventTriggerInfo& TriggeredEventHandler::getTriggerInfo() {
    return impl->triggerInfo;
}

Stage TriggeredEventHandler::getRequiredStage() const {
    return impl->requiredStage;
}

ScheduledEventHandler::~ScheduledEventHandler() {
}

class PeriodicEventHandler::PeriodicEventHandlerImpl {
public:
    PeriodicEventHandlerImpl(Real eventInterval) : eventInterval(eventInterval) {
        SimTK_APIARGCHECK1_ALWAYS(eventInterval > 0.0, "PeriodicEventHandlerImpl", "PeriodicEventHandlerImpl", "The interval was %d.  It must be > 0", eventInterval);
    }
    Real eventInterval;
};

PeriodicEventHandler::PeriodicEventHandler(Real eventInterval) {
    impl = new PeriodicEventHandlerImpl(eventInterval);
}

PeriodicEventHandler::~PeriodicEventHandler() {
    delete impl;
}

Real PeriodicEventHandler::getNextEventTime(const State& state, bool includeCurrentTime) const {
    Real currentTime = state.getTime();
    long count = (long)std::floor(currentTime/impl->eventInterval);
    volatile Real eventTime = count*impl->eventInterval;
    while (eventTime < currentTime || (eventTime == currentTime && !includeCurrentTime)) {
        count++;
        eventTime = count*impl->eventInterval;
    }
    return eventTime;
}

Real PeriodicEventHandler::getEventInterval() const {
    return impl->eventInterval;
}

void PeriodicEventHandler::setEventInterval(Real eventInterval) {
    SimTK_APIARGCHECK1_ALWAYS(eventInterval > 0.0, "PeriodicEventHandler", "setEventInterval", "The interval was %d.  It must be > 0", eventInterval);
    impl->eventInterval = eventInterval;
}

} // namespace SimTK
