/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
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

#include "SimTKcommon/internal/EventReporter.h"

namespace SimTK {

class TriggeredEventReporter::TriggeredEventReporterImpl {
public:
    TriggeredEventReporterImpl(Stage requiredStage) : requiredStage(requiredStage) {
    }
    System::EventTriggerInfo triggerInfo;
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

System::EventTriggerInfo& TriggeredEventReporter::getTriggerInfo() {
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
    int count = std::floor(currentTime/impl->eventInterval);
    Real eventTime = count*impl->eventInterval;
    if (eventTime < currentTime || (eventTime == currentTime && !includeCurrentTime)) {
        count++;
        eventTime = count*impl->eventInterval;
    }
    return eventTime;
}

Real PeriodicEventReporter::getEventInterval() {
    return impl->eventInterval;
}

void PeriodicEventReporter::setEventInterval(Real eventInterval) {
    SimTK_APIARGCHECK1_ALWAYS(eventInterval > 0.0, "PeriodicEventReporter", "setEventInterval", "The interval was %d.  It must be > 0", eventInterval);
    impl->eventInterval = eventInterval;
}

} // namespace SimTK
