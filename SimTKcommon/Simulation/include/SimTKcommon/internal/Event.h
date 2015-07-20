#ifndef SimTK_SimTKCOMMON_EVENT_H_
#define SimTK_SimTKCOMMON_EVENT_H_

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
 * This file declares the types needed for Simbody's support for Events.
 */

#include "SimTKcommon/basics.h"
#include "SimTKcommon/internal/Event_Defs.h"
#include "SimTKcommon/internal/EventAction.h"


#include <memory>
#include <string>
#include <utility>

namespace SimTK {
    
//==============================================================================
//                                  EVENT
//==============================================================================
/** An event is "something that happens" during a Study that is advancing
through time. Its occurrence interrupts the normal flow of computation, allowing
an event handler to adjust the System state prior to resuming the Study.

The occurrence of an event is signaled by EventTrigger objects that are 
constructed to watch time and state for specified conditions. Many triggers may
cause the same Event, and EventTrigger objects may come and go during the
simulation. The two most common kinds of triggers are timers and witnesses,
represented by EventTrigger::Timer and EventTrigger::Witness objects. Timers 
trigger Events at regular intervals or specific times; witnesses are scalar 
functions designed to cross zero upon occurrence of an event.

An Event is associated with a list of actions to be taken when the Event is
triggered. Actions are represented by EventAction objects. Actions can be 
reports, which leave the trajectory unchanged, or state changes that do affect 
the subsequent trajectory.

Events are System-level resources and each Event is assigned a System-unique
small integer EventId. Some events are predefined by the System, others
are added during extended construction of the System, by Subsystems and by 
elements within Subystems. Actions can be added to existing Events during
extended construction, and new Trigger objects can be added to trigger existing
Events at any time. **/
class SimTK_SimTKCOMMON_EXPORT Event {
public:
    class Initialization; // predefined Events
    class TimeAdvanced;
    class Termination;
    class ExtremeValueIsolated;

    /** Destructor is public and virtual for good cleanup behavior. **/
    virtual ~Event() = default;

    /** Create a new copy of this %Event on the heap, including copies of all
    the EventAction objects it contains. **/
    Event* clone() const {
        Event* p = cloneVirtual();
        return p;
    }

    /** Return the description of this event that was supplied during
    construction or via setEventDescription(). **/
    const std::string& getEventDescription() const 
    {   return m_eventDescription; }

    Event& setEventDescription(const std::string& description) {
        m_eventDescription = description;
        return *this;
    }

    /** The EventId is set when this %Event is adopted by a System and will be
    invalid before then. **/
    EventId getEventId() const {return m_eventId;}

    /** Return the number of times this %Event has been triggered since 
    the count was last reset. Multiple triggers simultaneously causing the
    same %Event to occur only counts as one occurrence.
    @see resetStatistics() **/
    long long getNumOccurrences() const {return m_numOccurrences;}

    /** Increment the number of occurrences by 1. This is for use by time
    steppers when they determine that %Event has occurred. This 
    method is const because the counter is mutable. Be sure to call this just
    once if multiple triggers simultaneously cause the same %Event. **/
    void noteOccurrence() const {++m_numOccurrences;}

    /** Reset to zero any statistics being kept by this %Event. **/
    void resetStatistics() {
        m_numOccurrences=0;  // reset base class statistics
        resetStatisticsVirtual();
    }

    /** Add an EventAction to be taken when this %Event is triggered. The 
    %Event takes over ownership of the heap-allocated EventAction and will 
    destruct it when the %Event is destructed. The order in which actions are 
    added is maintained and the actions will be invoked in that order, except 
    that *all* actions that make state changes will be executed before *any* 
    actions that are just reports. 
    @returns The index assigned to this EventAction object within 
    this %Event. **/
    EventActionIndex adoptEventAction(EventAction* actionp);

    /** Return the number of EventAction objects associated with this
    %Event. **/
    int getNumActions() const {return (int)m_actions.size();}

    /** Obtain a const reference to one of this %Event's EventAction objects, 
    given the index of that EventAction as returned from adoptEventAction(). **/
    const EventAction& getAction(EventActionIndex index) const
    {   return *m_actions[index]; }

    /** Obtain a writable reference to one of this %Event's EventAction objects, 
    given the index of that EventAction as returned from adoptEventAction(). **/
    EventAction& updAction(EventActionIndex index)
    {   return *m_actions[index]; }

    /** Does this Event have any actions that can change the state of a study?
    If not it has only reporting actions. **/
    bool hasChangeAction() const {return m_hasChangeAction;}

    /** Does this Event have any reporting actions? If so these actions will
    be taken after all change actions have been taken. **/
    bool hasReportAction() const {return m_hasReportAction;}

    /** This method is invoked when this %Event has been triggered but performs
    only Report EventActions. **/
    void performReportActions
       (const Study&            study,
        const EventTriggers&    triggers) const;

    /** This method is invoked when this %Event has been triggered but performs
    only Change Actions. **/
    void performChangeActions
       (Study&                  study,
        const EventTriggers&    triggers,
        EventChangeResult&      result) const;



    /** This is useful for debugging; it translates an EventCause into a 
    readable string. **/
    static const char* getCauseName(EventCause);

    /** EventTrigger::Witness triggers respond to zero crossings of their 
    associated witness function. This enum defines constants for use in 
    specifying which kind of zero crossing has been seen, or which kinds are 
    considered interesting. For the latter purpose, these can be or'ed 
    together to make a mask. **/
    enum TriggerDirection {
        NoEventTrigger          =0x0000,    // must be 0

        PositiveToNegative      =0x0001,    // 1
        NegativeToPositive      =0x0002,    // 2

        Falling                 =(PositiveToNegative), // 1
        Rising                  =(NegativeToPositive), // 2
        AnySignChange           =(PositiveToNegative|NegativeToPositive)    // 3
    };

    /** This is useful for debugging; it translates an EventTrigger or a mask
    formed by a union of EventTriggers, into a readable string. **/
    static std::string 
    eventTriggerDirectionString(TriggerDirection direction);


    /** Classify a before/after sign transition. Before and after must both
    be -1,0, or 1 as returned by the SimTK::sign() function applied to
    the trigger function value at the beginning and end of a step. **/
    static TriggerDirection classifyTransition(int before, int after) {
        if (before==after)
            return NoEventTrigger;
        if (before==0)
            return NoEventTrigger; // Do not report transitions away from zero.
        if (before==1)
            return PositiveToNegative;
        // before==-1
        return NegativeToPositive;
    }

    /** Given an observed transition, weed out ignorable ones using the supplied
    mask. That is, the return will indicate NoEventTrigger unless the original
    Trigger was present in the mask. **/
    static TriggerDirection maskTransition(TriggerDirection transition, 
                                           TriggerDirection mask) 
    {
        // we're depending on NoEventTrigger==0
        return TriggerDirection(transition & mask); 
    }

protected:
    explicit Event(const std::string& eventDescription);
    explicit Event(std::string&& eventDescription);

    Event() = delete;
    Event(const Event&) = default;
    //Event(const Event&&) = default;
    Event& operator=(const Event&) = delete;
    Event& operator=(const Event&&) = delete;

    /** Derived classes must implement this method, like this:
    @code
    class MyEvent : public Event {
        MyEvent* cloneVirtual() const override {
            return new MyEvent(*this);
        }
    };
    @endcode **/
    virtual Event* cloneVirtual() const = 0;

    /** If your derived %Event keeps statistics, reset them to zero. **/
    virtual void resetStatisticsVirtual() {}
private:
friend class SystemGlobalSubsystem;
    // TOPOLOGY STATE
    std::string         m_eventDescription;
    EventId             m_eventId;  // Set when Event is adopted by a System

    Array_<ClonePtr<EventAction>,EventActionIndex>  m_actions;
    bool                m_hasReportAction;
    bool                m_hasChangeAction;
    bool                m_allowsInterpolationForReport;

    // STATISTICS
    mutable long long   m_numOccurrences;

    void initialize() {
        // don't clear the event description here
        m_eventId.invalidate();
        m_hasReportAction = m_hasChangeAction = false;
        m_allowsInterpolationForReport = true;
        m_numOccurrences = 0;
    }
};

using EventAndCausesPair = std::pair<const Event*, EventTriggers>;

class EventsAndCauses : public Array_<EventAndCausesPair> {
public:
    using Array_::Array_;
};


//==============================================================================
//                           EVENT :: INITIALIZATION
//==============================================================================
/** This event occurs at the start of any Study to permit one-time initial
condition actions to be taken. Every System includes one of these events; you
will not need to create one yourself. **/
class Event::Initialization : public Event {
public:
    /** Don't create one of these; they are built in to every System. **/
    Initialization() : Event("Initialization") {}

    /** Return true if the concrete type of the given Event is 
    Event::Initialization. **/
    static bool isA(const Event& e)
    {   return dynamic_cast<const Initialization*>(&e) != nullptr; }

    /** Downcast a const reference to an Event to a const reference to this
    type Event::Initialization. A std::bad_cast exception is thrown if the type
    is wrong, at least in Debug builds. **/
    static const Initialization& downcast(const Event& e)
    {   return SimTK_DYNAMIC_CAST_DEBUG<const Initialization&>(e); }

    /** Downcast a writable reference to an Event to a writable reference to 
    this type Event::Initialization. A std::bad_cast exception is thrown if 
    the type is wrong, at least in Debug builds. **/
    static Initialization& updDowncast(Event& e)
    {   return SimTK_DYNAMIC_CAST_DEBUG<Initialization&>(e); }
private:
    Initialization* cloneVirtual() const override 
    {   return new Initialization(*this); }
};


//==============================================================================
//                           EVENT :: TIME ADVANCED
//==============================================================================
/** This event occurs during a time-stepping study whenever time advances
irreversably. Every System includes one of these events; you will not need to 
create one yourself. **/
class Event::TimeAdvanced : public Event {
public:
    /** Don't create one of these; they are built in to every System. **/
    TimeAdvanced() : Event("Time advanced") {}

    /** Return true if the concrete type of the given Event is 
    Event::TimeAdvanced. **/
    static bool isA(const Event& e)
    {   return dynamic_cast<const TimeAdvanced*>(&e) != nullptr; }

    /** Downcast a const reference to an Event to a const reference to this
    type Event::TimeAdvanced. A std::bad_cast exception is thrown if the type
    is wrong, at least in Debug builds. **/
    static const TimeAdvanced& downcast(const Event& e) {
        return SimTK_DYNAMIC_CAST_DEBUG<const TimeAdvanced&>(e);
    }

    /** Downcast a writable reference to an Event to a writable reference to 
    this type Event::TimeAdvanced. A std::bad_cast exception is thrown if 
    the type is wrong, at least in Debug builds. **/
    static TimeAdvanced& updDowncast(Event& e)
    {   return SimTK_DYNAMIC_CAST_DEBUG<TimeAdvanced&>(e); }
private:
    TimeAdvanced* cloneVirtual() const override 
    {   return new TimeAdvanced(*this); }
};


//==============================================================================
//                            EVENT :: TERMINATION
//==============================================================================
/** This event occurs at the end of any Study to permit one-time final
actions to be taken. Every System includes one of these events; you
will not need to create one yourself. **/
class Event::Termination : public Event {
public:
    /** Don't create one of these; they are built in to every System. **/
    Termination() : Event("Termination") {}

    /** Return true if the concrete type of the given Event is 
    Event::Termination. **/
    static bool isA(const Event& e)
    {   return dynamic_cast<const Termination*>(&e) != nullptr; }

    /** Downcast a const reference to an Event to a const reference to this
    type Event::Termination. A std::bad_cast exception is thrown if the type
    is wrong, at least in Debug builds. **/
    static const Termination& downcast(const Event& e)
    {   return SimTK_DYNAMIC_CAST_DEBUG<const Termination&>(e); }

    /** Downcast a writable reference to an Event to a writable reference to 
    this type Event::Termination. A std::bad_cast exception is thrown if 
    the type is wrong, at least in Debug builds. **/
    static Termination& updDowncast(Event& e)
    {   return SimTK_DYNAMIC_CAST_DEBUG<Termination&>(e); }

private:
    Termination* cloneVirtual() const override 
    {   return new Termination(*this); }
};


//==============================================================================
//                       EVENT :: EXTREME VALUE ISOLATED
//==============================================================================
/** This Event marks the occurrence of a minimum or maximum value of a function
that is being monitored, typically a Measure. The function's time derivative is
used as an EventTrigger::Witness; when the derivative changes sign we have seen
an extreme value. Normally no further action is required when this event occurs
since it will already have localized the extreme value. However, for debugging
it can be useful to report these occurrences, and perhaps which trigger(s)
caused the Event to occur. **/
class Event::ExtremeValueIsolated : public Event {
public:
    /** Don't create one of these; they are built in to every System. **/
    ExtremeValueIsolated() : Event("Extreme value isolated") {}

    /** Return true if the concrete type of the given Event is 
    Event::ExtremeValueIsolated. **/
    static bool isA(const Event& e)
    {   return dynamic_cast<const ExtremeValueIsolated*>(&e) != nullptr; }

    /** Downcast a const reference to an Event to a const reference to this
    type Event::Termination. A std::bad_cast exception is thrown if the type
    is wrong, at least in Debug builds. **/
    static const ExtremeValueIsolated& downcast(const Event& e)
    {   return SimTK_DYNAMIC_CAST_DEBUG<const ExtremeValueIsolated&>(e); }

    /** Downcast a writable reference to an Event to a writable reference to 
    this type Event::ExtremeValueIsolated. A std::bad_cast exception is thrown if 
    the type is wrong, at least in Debug builds. **/
    static ExtremeValueIsolated& updDowncast(Event& e)
    {   return SimTK_DYNAMIC_CAST_DEBUG<ExtremeValueIsolated&>(e); }

private:
    ExtremeValueIsolated* cloneVirtual() const override 
    {   return new ExtremeValueIsolated(*this); }
};

//==============================================================================
//                           EVENT CHANGE RESULT
//==============================================================================
/** TODO: This is for accumulating the results from a chain of Event Change
Actions. For example, if any one of them wants to terminate, that should be
the final result. This will replace HandleEventsResults.
**/
class EventChangeResult {
public:
    EventChangeResult() {clear();}

    /** These are ordered; we keep the highest (worst) Status value we encounter
    in processing a series of event actions. A message is recorded if an action
    returns ShouldTerminate or Failed status. **/
    enum Status {
        /** All attempted event change actions were successful and time stepping
        may continue. This is the status after construction or a call to 
        clear(); that is, if we have yet to perform any event actions. **/
        Succeeded               = 0,
        /** The event change actions were successful but the event requires 
        time stepping to terminate. An explanation may have been placed in
        the message argument; this will be from the first event action that
        returned ShouldTerminate status. **/
        ShouldTerminate         = 1,
        /** An event change action was unable to successfully handle the
        event. This is likely to be a fatal error. A human-readable 
        explanation is in the message argument. Execution of event actions is
        terminated by the first one to return Failed status. **/
        Failed                  = 2    
    };

    /** Restore this object to its default-constructed state, with the return
    status set to Succeeded. **/
    void clear() {
        m_exitStatus = Succeeded;
        m_message.clear();
        m_lowestModifiedStage = Stage::Infinity;
    }

    /** Return the status resulting from all the event actions executed so 
    far; the worst one encountered is the overally status. **/
    Status getExitStatus() const {return m_exitStatus;}

    /** Return a human-readable name for the given exit status. **/
    static const char *getExitStatusName(Status status) {
        switch (status) {
        case Succeeded: return "Succeeded";
        case ShouldTerminate: return "ShouldTerminate";
        case Failed: return "Failed";
        }
        return "UNKNOWN EventChangeResult::Status";
    }

    /** Return the lowest (earliest) Stage state variable modified by any
    executed event action. **/
    Stage getLowestModifiedStage() const {return m_lowestModifiedStage;}

    /** Return a human-readable message supplied by whichever event action
    caused us to set the exit status as it is now. **/
    const std::string& getMessage() const {return m_message;}

    /** Report the exit status of an event action. If it is "more severe" than
    the worst exit status we've seen so far, then we'll update the exit status
    and record the message. Otherwise both will be ignored. **/
    void reportExitStatus(Status status, const std::string& message) {
        if (status > m_exitStatus) {
            m_exitStatus = status;
            m_message = message;
        }
    }
    /** Alternate signature that does not take a message; use this for 
    successful returns. **/
    void reportExitStatus(Status status) {
        if (status > m_exitStatus) {
            m_exitStatus = status;
            m_message.clear();
        }
    }

    /** Report the lowest (earliest) Stage modified by any executed event
    action. Individual event actions need not set this; it will be determined
    at the end by examing stage version numbers. It can't be determined by
    looking just at the final realization stage of the modified State because
    event actions may realize the state themselves after changing it. **/
    void setLowestModifiedStage(Stage stage) {m_lowestModifiedStage=stage;}

private:
    Status      m_exitStatus;
    std::string m_message;
    Stage       m_lowestModifiedStage;
};

//==============================================================================
//              HANDLE EVENTS OPTIONS and HANDLE EVENTS RESULTS
//==============================================================================
/** Options for the handleEvent() method. Accuracy should be be set by the
caller, but if not the default is 1e-4. **/
class HandleEventsOptions {
public:
    enum Option {
        /** Take all defaults. **/
        None            = 0x0000,
        /** Normally failure to meet the accuracy requirements throws an
        exception. This will force the handleEvent() method to quietly return bad 
        status instead. **/
        DontThrow       = 0x0001,
        /** Use the stricter infinity (max absolute value) norm rather than
        the default RMS norm to determine when accuracy has been achieved. **/
        UseInfinityNorm = 0x0002
    };


    HandleEventsOptions() {clear();}
    explicit HandleEventsOptions(Real accuracy) 
    {   clear(); setAccuracy(accuracy); }
    explicit HandleEventsOptions(Option opt)
    {   clear(); setOption(opt); }

    /** Restore this object to its default-constructed state (no options
    selected, default accuracy). A reference to the
    newly-cleared object is returned. **/
    HandleEventsOptions& clear() 
    {   optionSet=0; setAccuracyDefaults(); return *this; }

    /** The norm of the constraint errors must be driven to below this value
    for a project() to be considered successful. Normally an RMS norm is used
    but you can override that to use an infinity norm instead. **/
    HandleEventsOptions& setAccuracy(Real accuracy) {
        assert(accuracy > 0);
        requiredAccuracy = accuracy;
        return *this;
    }

    /** Remove a given option from the set. Nothing happens if the option wasn't
    already set. **/
    HandleEventsOptions& clearOption(Option opt) 
    {   optionSet &= ~(unsigned)opt; return *this; }
    /** Select a given option from the set. Nothing happens if the option wasn't
    already set. **/
    HandleEventsOptions& setOption  (Option opt) 
    {   optionSet |= (unsigned)opt; return *this; }

    /** Return the current value for the accuracy option. **/
    Real getAccuracy()       const {return requiredAccuracy;}

    bool isOptionSet(Option opt) const {return (optionSet&(unsigned)opt) != 0;}

    static Real getDefaultAccuracy() {return Real(1e-4);}

    // Set operators: not, or, and, set difference
    HandleEventsOptions& operator|=(const HandleEventsOptions& opts) 
    {   optionSet |= opts.optionSet; return *this; }
    HandleEventsOptions& operator&=(const HandleEventsOptions& opts) 
    {   optionSet &= opts.optionSet; return *this; }
    HandleEventsOptions& operator-=(const HandleEventsOptions& opts) 
    {   optionSet &= ~opts.optionSet; return *this; }

    HandleEventsOptions& operator|=(Option opt) {setOption(opt); return *this;}
    HandleEventsOptions& operator-=(Option opt) {clearOption(opt); return *this;}

private:
    Real     requiredAccuracy;
    unsigned optionSet;

    void setAccuracyDefaults() {
        requiredAccuracy = getDefaultAccuracy();
    }
};

/** Results returned by the handleEvent() method. In addition to return 
status, this records the lowest stage in the state that was modified by
the handler. The caller can use this to determine how much reinitialization
is required before time stepping can proceed. **/
class HandleEventsResults {
public:
    HandleEventsResults() {clear();}

    enum Status {
        /** This object has not been filled in yet and holds no results. **/
        Invalid                 = -1,
        /** The handleEvent() operation was successful and time stepping
        may continue. **/
        Succeeded               = 0,
        /** The handleEvent() call was successful but the event requires 
        time stepping to terminate. An explanation may have been placed in
        the message argument. **/
        ShouldTerminate         = 1,
        /** The handleEvent() call was unable to successfully handle the
        event. This is likely to be a fatal error. A human-readable 
        explanation is in the message argument. **/
        Failed                  = 2    
    };

    /** Restore this object to its default-constructed state, with the return
    status set to Invalid. **/
    HandleEventsResults& clear() {
        m_exitStatus = Invalid;
        m_anyChangeMade = false;
        m_lowestModifiedStage = Stage::Infinity; // i.e., nothing modified
        m_message.clear();
        return *this;
    }
    bool    isValid()           const {return m_exitStatus != Invalid;}
    Status  getExitStatus()     const {return m_exitStatus;}

    bool getAnyChangeMade()     const 
    {   assert(isValid()); return m_anyChangeMade; }
    Stage getLowestModifiedStage() const 
    {   assert(isValid()); return m_lowestModifiedStage; }
    const String& getMessage() const
    {   assert(isValid()); return m_message; }

    HandleEventsResults& setExitStatus(Status status) 
    {   m_exitStatus=status; return *this; }
    HandleEventsResults& setAnyChangeMade(bool changeMade) 
    {   m_anyChangeMade=changeMade; return *this; }
    HandleEventsResults& setLowestModifiedStage(Stage stage) 
    {   m_lowestModifiedStage=stage; return *this; }
    HandleEventsResults& setMessage(const String& message) 
    {   m_message=message; return *this; }
private:
    Status  m_exitStatus;
    bool    m_anyChangeMade;
    Stage   m_lowestModifiedStage;
    String  m_message;
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_EVENT_H_
