/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2015 Stanford University and the Authors.           *
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


/**@file
 *
 * Implementation of non-inline methods from the EventTrigger::Timer classes.
 */

#include "SimTKcommon/basics.h"
#include "SimTKcommon/internal/Event.h"
#include "SimTKcommon/internal/EventTrigger.h"
#include "SimTKcommon/internal/EventTrigger_Timer.h"

#include <cassert>
#include <algorithm>
#include <initializer_list>
#include <cmath>

using namespace SimTK;

//==============================================================================
//                         DESIGNATED EVENT TIMER
//==============================================================================

EventTrigger::Timer::Designated::Designated(const std::string& description,
                                            double t) 
:   Super(description) {
    checkTime("Designated", t);
    m_triggerTimes.push_back(t);
}

void EventTrigger::Timer::Designated::
insertDesignatedTime(double t) {
    checkTime("insertDesignatedTime", t);
    if (m_triggerTimes.empty() || t > m_triggerTimes.back()) {
        m_triggerTimes.push_back(t);
        return;
    }
    // Find the first time that is >= t; cost is log(n). There must be such
    // an element since we already checked if t is the largest.
    auto p = std::lower_bound(m_triggerTimes.begin(),m_triggerTimes.end(),t);
    assert(p != m_triggerTimes.end()); 
    if (t != *p) // ignore if t==*p
        m_triggerTimes.insert(p, t);
}

void EventTrigger::Timer::Designated::
checkTime(const char* methodName, double t) const {
    SimTK_APIARGCHECK2_ALWAYS(t >= 0, 
        "EventTrigger::Timer::Designated", methodName,
        "Illegal trigger time %g (%s).", t,
        getTriggerDescription().c_str());
}

void EventTrigger::Timer::Designated::
constructHelper() {
    for (auto t : m_triggerTimes) checkTime("Designated", t);
    std::sort(m_triggerTimes.begin(), m_triggerTimes.end()); 
    auto newEnd = std::unique(m_triggerTimes.begin(), m_triggerTimes.end());
    m_triggerTimes.resize((unsigned)std::distance(m_triggerTimes.begin(), 
                                                  newEnd));
}

void EventTrigger::Timer::Designated::
insertHelper(Array_<double>& times) {
    for (auto t : times) checkTime("insertDesignatedTimes", t);
    std::sort(times.begin(), times.end());
    Array_<double> merged(times.size() + m_triggerTimes.size());
    std::merge(times.begin(), times.end(),
                m_triggerTimes.begin(), m_triggerTimes.end(),
                merged.begin());
    auto newEnd = std::unique(merged.begin(), merged.end());
    merged.resize((unsigned)std::distance(merged.begin(), newEnd));
    m_triggerTimes.swap(merged);
}

// The System is not used in this calculation.
double EventTrigger::Timer::Designated::
calcTimeOfNextTriggerVirtual(const System&, const State& state, 
                             double timeOfLastTrigger) const 
{
    const double currentTime = state.getTime();
    // Find the first time that is >= currentTime (O(log n) time).
    auto p = std::lower_bound(m_triggerTimes.begin(), m_triggerTimes.end(),
                              currentTime);
    if (p == m_triggerTimes.end())
        return Infinity; // ran off the end

    // Selected event time is at least currentTime, but we still don't want
    // one earlier than timeOfLastTrigger (which might have been at the
    // currentTime).
    while (*p <= timeOfLastTrigger) {
        if (++p == m_triggerTimes.end())
            return Infinity; // event time was the last entry
    }
    return *p;
}


//==============================================================================
//                            PERIODIC EVENT TIMER
//==============================================================================

// System is not used for this calculation.
double EventTrigger::Timer::Periodic::
calcTimeOfNextTriggerVirtual(const System&, const State& state, 
                             double timeOfLastTrigger) const 
{
    const double currentTime = state.getTime();
    long long count = (long long)std::floor(currentTime/m_period);
    if (currentTime == 0 && !getShouldTriggerAtZero())
        count = 1;
    double eventTime = count*m_period;
    
    while (eventTime < currentTime || eventTime <= timeOfLastTrigger) {
        ++count; // skip this one
        eventTime = count * m_period;
    }
    return eventTime;
}


