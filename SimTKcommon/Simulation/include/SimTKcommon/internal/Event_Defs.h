#ifndef SimTK_SimTKCOMMON_EVENT_DEFS_H_
#define SimTK_SimTKCOMMON_EVENT_DEFS_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-15 Stanford University and the Authors.        *
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

/** @file
 *
 * This file declares just the basic definitions needed by Simbody's various
 * Event header files so that they don't all have to include one another.
 */

#include "SimTKcommon/basics.h"

#include <memory>
#include <string>
#include <utility>

namespace SimTK {
class System;
class State;
class Event;
class EventAction;
class EventTrigger;
class EventChangeResult;
class HandleEventsOptions;
class HandleEventsResults;
    
/** @class SimTK::EventId
This class provides type-safe, System-unique integer IDs for Event objects
in a given System. **/
SimTK_DEFINE_UNIQUE_INDEX_TYPE(EventId);

/** @class SimTK::EventTriggerId
This class provides type-safe, System-unique integer IDs for EventTrigger 
objects in a System or one of its States. **/
SimTK_DEFINE_UNIQUE_INDEX_TYPE(EventTriggerId);

/** @class SimTK::EventActionIndex
Type-safe index for EventAction objects within an Event. These are
only unique within a particular Event. **/
SimTK_DEFINE_UNIQUE_INDEX_TYPE(EventActionIndex);

/** A container for const pointers to EventTrigger objects. **/
using EventTriggers = Array_<const EventTrigger*>;

/** These are all the possible causes for events.

@par Initialization
A Study has performed its own initialization and is about to start.

@par Triggered
An event trigger function underwent a monitored sign transition.

@par Scheduled
An integrator reached a previously-scheduled time for the Event to occur.

@par TimeAdvanced
An integrator completed an internal step, meaning that it has reached
a point where time has advanced irreversibly.

@par Signaled
A flag in the State has been explicitly set, meaning that a particular
Event has occurred. Anyone with write access to a State can set these,
but typically they are set in event handlers associated with one of
the other kinds of events.

@par Termination
The Study has finished. If a Termination event handler signals more
Events, those signaled events are not processed by the Study; that is,
the signals remain set in the final State.

In case several of these causes are detected in a single step, they are 
sequentialized in the order shown, like this:

1. The occurrence of triggered events is reported and the triggering state 
    and a list of triggered events are passed to the event handler for 
    processing (meaning the state, but not the time, is modified). Note 
    that simultaneity *within* the set of triggered events may also require 
    special handling; we're not talking about that here, just simultaneity 
    of *causes*.
2. Next, using the state resulting from step 1, the time is checked to see 
    if scheduled events have occurred. If so, a list of those events is 
    passed to the event handler for processing.
3. Next, if this system has requested time-advanced events, the event 
    handler is called with the state that resulted from step 2 and the "time 
    advanced" cause noted. No event list is passed in that case. The state 
    may be modified.
4. Last, if the final time has been reached or if any of the event handlers
    asked for termination, we pass the state to the event handler again 
    noting that we have reached termination. The state may be modified and 
    the result will be the final state of the simulation.
**/
enum class EventCause {
    Initialization  = 1,
    Triggered       = 2,
    Scheduled       = 3,
    TimeAdvanced    = 4,
    Signaled        = 5,
    Termination     = 6,
    Invalid         = -1
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_EVENT_DEFS_H_
