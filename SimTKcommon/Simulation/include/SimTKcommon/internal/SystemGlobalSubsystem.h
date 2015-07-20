#ifndef SimTK_SimTKCOMMON_SYSTEM_GLOBAL_SUBSYSTEM_H_
#define SimTK_SimTKCOMMON_SYSTEM_GLOBAL_SUBSYSTEM_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-15 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors: Peter Eastman                                                *
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

#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"
#include "SimTKcommon/internal/State.h"
#include "SimTKcommon/internal/System.h"
#include "SimTKcommon/internal/Subsystem.h"

#include <cassert>

namespace SimTK {
class Study;
class ScheduledEventHandler;
class ScheduledEventReporter;
class TriggeredEventHandler;
class TriggeredEventReporter;

//==============================================================================
//                         SYSTEM GLOBAL SUBSYSTEM
//==============================================================================
/** This Subsystem manages resources that are global to the System rather than
associated with any particular Subsystem. For example, we maintain lists of
event handlers and reporters here. There is one of these automatically defined
in every System after it is constructed.

To obtain the global subsystem for a System, call getSystemGlobalSubsystem() or 
getSystemGlobalSubsystem() on it. Also, a System can be implicitly converted
to a Subsystem, in which case it actually returns a reference to this subsystem. 

@par Realization order
For the Stage::Topology and Stage::Model realizations, the 
%SystemGlobalSubsystem is realized *first* so that resources it allocates can 
be used during the same realization Stage by other Subsystems. For the later 
stages Stage::Instance through Stage::Report, this subsystem is realized *last*,
after all other subsystems. That way if EventTrigger functions realized here 
depend on computations done in other subsystems, those computations will have 
been completed by the time they are needed. **/
class SimTK_SimTKCOMMON_EXPORT SystemGlobalSubsystem : public Subsystem {
public:
    /** Don't call this constructor; it will be invoked automatically during 
    System construction to add this %Subsystem to the System. **/
    explicit SystemGlobalSubsystem(System& sys);

    /** Add a ScheduledEventHandler to the System. This must be called before 
    the Model stage is realized. The subsystem assumes ownership of the object 
    passed to this method, and will delete it when the subsystem is deleted. **/
    void adoptEventHandler(ScheduledEventHandler* handler);

    /** Add a TriggeredEventHandler to the System. This must be called before 
    the Model stage is realized. The subsystem assumes ownership of the object 
    passed to this method, and will delete it when the subsystem is deleted. **/
    void adoptEventHandler(TriggeredEventHandler* handler);

    /** Add a ScheduledEventReporter to the System. This must be called before 
    the Model stage is realized. The subsystem assumes ownership of the object 
    passed to this method, and will delete it when the subsystem is deleted. **/
    void adoptEventReporter(ScheduledEventReporter* reporter);

    /** Add a TriggeredEventReporter to the System. This must be called before 
    the Model stage is realized. The subsystem assumes ownership of the object 
    passed to this method, and will delete it when the subsystem is deleted. **/
    void adoptEventReporter(TriggeredEventReporter* reporter);

    /** Add a new Event to this System and obtain a unique, small-integer
    EventId for it. This subsystem takes over ownership of the heap-allocated 
    Event object. **/
    EventId adoptEvent(Event* eventp);

    /** Add a new EventTrigger to this System and obtain a unique, 
    small-integer EventTriggerId for it. This subsystem takes over ownership 
    of the heap-allocated Event object. **/
    EventTriggerId adoptEventTrigger(EventTrigger* triggerp);

    /** Return the number n of Event objects contained in this System. The
    corresponding EventIds are EventId(0) through EventId(n-1). Don't 
    confuse this with the number of event occurrences. **/
    int getNumEvents() const;

    /** Return a const reference to an Event specified by EventId. Throws an 
    exception if the EventId is invalid or if there is no Event corresponding to
    it in this %System. **/
    const Event& getEvent(EventId id) const;

    /** Return a writable reference to an Event specified by EventId. Throws an 
    exception if the EventId is invalid or if there is no Event corresponding to
    it in this System. **/
    Event& updEvent(EventId id);

    /** Check whether this %System has an Event with the given EventId. **/
    bool hasEvent(EventId id) const;

    /** Return the number n of EventTrigger objects contained in this %System.
    The corresponding EventTriggerIds are EventTriggerId(0) through 
    EventTriggerId(n-1). **/
    int getNumEventTriggers() const;

    /** Return a const reference to an EventTrigger specified by EventTriggerId. 
    Throws an exception if the EventTriggerId is invalid or if there is no 
    EventTrigger corresponding to it in this %System. **/
    const EventTrigger& getEventTrigger(EventTriggerId id) const;

    /** Return a writable reference to an EventTrigger specified by EventTriggerId. 
    Throws an exception if the EventTriggerId is invalid or if there is no 
    EventTrigger corresponding to it in this %System. **/
    EventTrigger& updEventTrigger(EventTriggerId id);

    /** Check whether this %System has an EventTrigger with the given 
    EventTriggerId. **/
    bool hasEventTrigger(EventTriggerId id) const;

    /** Return the EventId for the built-in Initialization Event. This Event
    is triggered once at the start of a simulation. **/
    EventId getInitializationEventId() const;
    /** Return the EventId for the built-in TimeAdvanced Event. This Event is
    triggered whenever an integrator successfully completes a step. **/
    EventId getTimeAdvancedEventId() const;
    /** Return the EventId for the built-in Termination Event. This Event is
    triggered once at the end of a simulation. **/
    EventId getTerminationEventId() const;
    /** Return the EventId for the built-in ExtremeValueIsolated Event. This
    Event is triggered whenever a watched value reaches and extremum of 
    interest (could be a maximum, minimum, or either). **/
    EventId getExtremeValueIsolatedEventId() const;

    /** Return the EventTriggerId for the built-in InitializationTrigger/ **/
    EventTriggerId getInitializationTriggerId() const;
    /** Return the EventTriggerId for the built-in TimeAdvancedTrigger. **/
    EventTriggerId getTimeAdvancedTriggerId() const;
    /** Return the EventTriggerId for the built-in TerminationTrigger. **/
    EventTriggerId getTerminationTriggerId() const;

    /** Return a list of event witnesses currently present for the %System and
    internal %State from the given `study`. **/   
    void findActiveEventWitnesses
      (const Study&                         study, 
       Array_<const EventTrigger::Witness*,
              ActiveWitnessIndex>&          witnesses) const;

    void findActiveEventTimers
       (const Study&                            study, 
        Array_<const EventTrigger::Timer*,
               ActiveTimerIndex>&               timers) const;

    void findNextScheduledEventTimes
       (const Study&        study,
        double              timeOfLastReport,
        double              timeOfLastChange,
        double&             timeOfNextReport,
        EventTriggers&      reportTimers,
        double&             timeOfNextChange,
        EventTriggers&      changeTimers) const;

    void noteEventOccurrence
       (const EventTriggers&        triggers,
        EventsAndCauses&            appendTriggeredEvents,
        Array_<EventId>&            appendIgnoredEvents) const;

    void performEventReportActions
       (const Study&                study,
        const EventsAndCauses&      triggeredEvents) const;

    void performEventChangeActions
       (Study&                      study,
        const EventsAndCauses&      triggeredEvents,
        EventChangeResult&          result) const;

    /** @name                    Deprecated 
    Don't use these methods; they are here temporarily for backwards 
    compatibility but should be replaced by the new methods indicated. **/
    /**@{**/
    /** (Deprecated) Use `adoptEventHandler()` instead; changed in 
    Simbody 4.0. **/
    void addEventHandler(ScheduledEventHandler* handler)
    {   adoptEventHandler(handler); }
    /** (Deprecated) Use `adoptEventHandler()` instead; changed in 
    Simbody 4.0. **/
    void addEventHandler(TriggeredEventHandler* handler)
    {   adoptEventHandler(handler); }
    /** (Deprecated) Use `adoptEventReporter()` instead; changed in 
    Simbody 4.0. **/
    void addEventReporter(ScheduledEventReporter* reporter)
    {   adoptEventReporter(reporter); }
    /** (Deprecated) Use `adoptEventReporter()` instead; changed in 
    Simbody 4.0. **/
    void addEventReporter(TriggeredEventReporter* reporter)
    {   adoptEventReporter(reporter); }
    /**@}**/

    /** @cond **/  // don't let doxygen see this private class
    class Guts;
    /** @endcond **/
private:
    const Guts& getGuts() const;
    Guts& updGuts();
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_SYSTEM_GLOBAL_SUBSYSTEM_H_
