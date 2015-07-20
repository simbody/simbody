#ifndef SimTK_SimTKCOMMON_EVENT_TRIGGER_H_
#define SimTK_SimTKCOMMON_EVENT_TRIGGER_H_

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

#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"
#include "SimTKcommon/internal/State.h"

#include <set>

namespace SimTK {

//==============================================================================
//                               EVENT TRIGGER
//==============================================================================
/** This is the base class for the various types of triggers that can cause the
occurrence of an Event. An %EventTrigger object checks for some condition that
occurs during a simulation; when that condition occurs the associated Event
objects are triggered.

Derived %EventTrigger classes represent particular conditions. For example, an
EventTrigger::Timer is a trigger that occurs at a particular time or times known
in advance, while an EventTrigger::Witness is a trigger that causes an Event 
occurrence when a particular state-dependent transition is observed.

Triggers can be built into a System or can be added and removed dynamically 
during simulation, in which case they reside in a State rather than
in the System. States are intended to be separate from any particular instance
of a System, which is what allows them to be serialized and then read back in
to be used with a new instance of the same System. So if you are implementing 
your own dynamic %EventTrigger, be sure to avoid use of pointers or references 
into the System. Rather, you should reference System objects by name, index, or 
id number. **/
class SimTK_SimTKCOMMON_EXPORT EventTrigger {
public:
    class Timer;
    class Witness;

    virtual ~EventTrigger() {}

    /** Make a copy of this Trigger, including the EventId list identifying the
    Event objects it triggers. **/
    EventTrigger* clone() const {return cloneVirtual();}

    /** Return the description of this trigger that was supplied during
    construction or via setTriggerDescription(). **/
    const std::string& getTriggerDescription() const 
    {   return m_triggerDescription; }

    EventTrigger& setTriggerDescription(const std::string& description) {
        m_triggerDescription = description;
        return *this;
    }

    /** Add an Event to this %Trigger by providing its EventId. We do not check
    whether an Event with this Id currently exists but simply record the given
    number. The order of added Events is preserved so when this %Trigger is seen
    the Event actions will occur in this order. If this EventId has already
    been added we'll quietly do nothing. **/
    EventTrigger& addEvent(EventId eventId) {
        SimTK_APIARGCHECK_ALWAYS(eventId.isValid(),
            "EventTrigger", "addEvent", "Expected a valid EventId.");
        auto p = std::find(m_eventIds.begin(), m_eventIds.end(), eventId);
        if (p == m_eventIds.end())
            m_eventIds.push_back(eventId);
        return *this;
    }

    /** Remove the Event given by `eventId` from this %Trigger. The order of any
    remaining events is preserved. If this is not on the list of triggered
    events, this method quietly does nothing. **/
    EventTrigger& removeEvent(EventId eventId) {
        SimTK_APIARGCHECK_ALWAYS(eventId.isValid(),
            "EventTrigger", "removeEvent", "Expected a valid EventId.");
        auto p = std::find(m_eventIds.begin(), m_eventIds.end(), eventId);
        if (p != m_eventIds.end())
            m_eventIds.erase(p);
        return *this;
    }

    /** Remove all events triggered by this %EventTrigger. Simulations will 
    still watch for this trigger, but nothing will happen as a result except 
    that the occurrence count will increase. **/
    void clearEvents() {m_eventIds.clear();}

    /** Returns true if the Event given by `eventId` is triggered by this 
    %EventTrigger. **/
    bool hasEvent(EventId eventId) const {
        SimTK_APIARGCHECK_ALWAYS(eventId.isValid(),
            "EventTrigger", "hasEvent", "Expected a valid EventId.");
        auto p = std::find(m_eventIds.begin(), m_eventIds.end(), eventId);
        return p != m_eventIds.end();
    }

    /** Return the number of events associated with this %Trigger. **/
    int getNumEvents() const {return m_eventIds.size();}

    /** Return the i'th EventId associated with this %Trigger. The index must
    be in range 0 to getNumEvents()-1. EventId ordering is the same order as
    they were added with addEvent(). **/
    EventId getEventId(int i) const {return m_eventIds[i];}

    /** Return the number of times this %Trigger condition has occurred since 
    the count was last reset. @see resetStatistics() **/
    long long getNumOccurrences() const {return m_count;}

    /** Increment the number of occurrences by 1. This is for use by time
    steppers when they determine that this trigger has been pulled. This 
    method is const because the counter is mutable. **/
    void noteOccurrence() const {++m_count;}

    /** Reset to zero any statistics being kept by this trigger. **/
    void resetStatistics() {
        m_count=0;  // reset base class statistics
        resetStatisticsVirtual();
    }

    /** The EventTriggerId is set when this %EventTrigger is adopted by a System
    or a State and will be invalid before then. **/
    EventTriggerId getEventTriggerId() const {return m_triggerId;}

protected:
    /** Construct an %EventTrigger base class object with no associated 
    events. A description can be useful for reporting and debugging. **/
    explicit EventTrigger(const std::string& triggerDescription)
    :   m_triggerDescription(triggerDescription), m_count(0) {}

    /** Invoke the concrete class copy constructor to return a new copy of
    the concrete %EventTrigger object. Make sure your copy constructor invokes
    the superclass copy constructor in its initializer list; the best way
    to do that is to design your %EventTrigger so that it can use the 
    compiler-generated default copy constructor. As discussed in the class
    documentation, your %EventTrigger must not contain pointers or references 
    into a System; otherwise the copy will be tied to the same instance of that
    System. **/
    virtual EventTrigger* cloneVirtual() const = 0;

    /** If your derived %EventTrigger keeps statistics, reset them to zero. **/
    virtual void resetStatisticsVirtual() {}

private:
friend class SystemGlobalSubsystem;
    std::string         m_triggerDescription;

    Array_<EventId>     m_eventIds;
    mutable long long   m_count;

    // Set when this EventTrigger is adopted by a System or State.
    EventTriggerId      m_triggerId; 
};



//-------------------------- INITIALIZATION TRIGGER ----------------------------
/** This is the EventTrigger that "causes" an Initialization event; it will be 
the listed trigger whenever a report or change action is invoked on the built-in
Initialization event.  Every System includes one of these triggers; you
will not need to create one yourself. **/
class InitializationTrigger: public EventTrigger {
public:
    /** Don't create one of these; they are built in to every System. **/
    explicit InitializationTrigger(EventId initEventId)
    :   EventTrigger("Study initialization") {addEvent(initEventId);}

    /** Return true if the concrete type of the given EventTrigger is 
    InitializationTrigger. **/
    static bool isA(const EventTrigger& et)
    {   return dynamic_cast<const InitializationTrigger*>(&et) != nullptr; }

    /** Downcast a const reference to an EventTrigger to a const reference to 
    this type InitializationTrigger. A std::bad_cast exception is thrown if the
    type is wrong, at least in Debug builds. **/
    static const InitializationTrigger& downcast(const EventTrigger& et)
    {   return SimTK_DYNAMIC_CAST_DEBUG<const InitializationTrigger&>(et); }

    /** Downcast a writable reference to an EventTrigger to a writable reference
    to this type InitializationTrigger. A std::bad_cast exception is thrown if 
    the type is wrong, at least in Debug builds. **/
    static InitializationTrigger& updDowncast(EventTrigger& et)
    {   return SimTK_DYNAMIC_CAST_DEBUG<InitializationTrigger&>(et); }

private:
    InitializationTrigger* cloneVirtual() const override
    {   return new InitializationTrigger(*this); }
};

//--------------------------- TIME ADVANCED TRIGGER ----------------------------
/** This is the EventTrigger that "causes" a TimeAdvanced event; it will be the
listed trigger whenever a report or change action is invoked on the built-in
TimeAdvanced event. Every System includes one of these triggers; you
will not need to create one yourself. **/
class TimeAdvancedTrigger: public EventTrigger {
public:
    /** Don't create one of these; they are built in to every System. **/
    explicit TimeAdvancedTrigger(EventId timeAdvancedEventId)
    :   EventTrigger("time advanced irreversibly") 
    {   addEvent(timeAdvancedEventId); }

    /** Return true if the concrete type of the given EventTrigger is 
    TimeAdvancedTrigger. **/
    static bool isA(const EventTrigger& et)
    {   return dynamic_cast<const TimeAdvancedTrigger*>(&et) != nullptr; }

    /** Downcast a const reference to an EventTrigger to a const reference to 
    this type TimeAdvancedTrigger. A std::bad_cast exception is thrown if the
    type is wrong, at least in Debug builds. **/
    static const TimeAdvancedTrigger& downcast(const EventTrigger& et)
    {   return SimTK_DYNAMIC_CAST_DEBUG<const TimeAdvancedTrigger&>(et); }

    /** Downcast a writable reference to an EventTrigger to a writable reference
    to this type TimeAdvancedTrigger. A std::bad_cast exception is thrown if 
    the type is wrong, at least in Debug builds. **/
    static TimeAdvancedTrigger& updDowncast(EventTrigger& et)
    {   return SimTK_DYNAMIC_CAST_DEBUG<TimeAdvancedTrigger&>(et); }

private:
    TimeAdvancedTrigger* cloneVirtual() const override
    {   return new TimeAdvancedTrigger(*this); }
};

//---------------------------- TERMINATION TRIGGER -----------------------------
/** This is the EventTrigger that "causes" a Termination event; it will be the
listed trigger whenever a report or change action is invoked on the built-in
Termination event. Every System includes one of these triggers; you
will not need to create one yourself. **/
class TerminationTrigger: public EventTrigger {
public:
    /** Don't create one of these; they are built in to every System. **/
    explicit TerminationTrigger(EventId termEventId)
    :   EventTrigger("Study terminated") {addEvent(termEventId);}

    /** Return true if the concrete type of the given EventTrigger is 
    TerminationTrigger. **/
    static bool isA(const EventTrigger& et)
    {   return dynamic_cast<const TerminationTrigger*>(&et) != nullptr; }

    /** Downcast a const reference to an EventTrigger to a const reference to 
    this type TerminationTrigger. A std::bad_cast exception is thrown if the
    type is wrong, at least in Debug builds. **/
    static const TerminationTrigger& downcast(const EventTrigger& et)
    {   return SimTK_DYNAMIC_CAST_DEBUG<const TerminationTrigger&>(et); }

    /** Downcast a writable reference to an EventTrigger to a writable reference
    to this type TerminationTrigger. A std::bad_cast exception is thrown if 
    the type is wrong, at least in Debug builds. **/
    static TerminationTrigger& updDowncast(EventTrigger& et)
    {   return SimTK_DYNAMIC_CAST_DEBUG<TerminationTrigger&>(et); }

private:
    TerminationTrigger* cloneVirtual() const override
    {   return new TerminationTrigger(*this); }
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_EVENT_TRIGGER_H_
