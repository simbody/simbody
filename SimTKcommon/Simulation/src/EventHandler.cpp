/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-15 Stanford University and the Authors.        *
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

#include "SimTKcommon/internal/EventHandler.h"

using namespace SimTK;

Real PeriodicEventHandler::
getNextEventTime(const State& state, bool includeCurrentTime) const {
    Real currentTime = state.getTime();
    long long count = (long long)std::floor(currentTime/m_eventInterval);
    Real eventTime = count*m_eventInterval;
    while (    eventTime < currentTime 
           || (eventTime == currentTime && !includeCurrentTime)) 
    {
        count++;
        eventTime = count*m_eventInterval;
    }
    return eventTime;
}

void PeriodicEventHandler::
setEventInterval(Real eventInterval) {
    SimTK_APIARGCHECK1_ALWAYS(eventInterval > 0, 
        "PeriodicEventHandler", "setEventInterval", 
        "The interval was %g but must be greater than zero.", eventInterval);
    m_eventInterval = eventInterval;
}

