#ifndef SimTK_SimTKCOMMON_EVENT_REPORTER_H_
#define SimTK_SimTKCOMMON_EVENT_REPORTER_H_

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

#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"
#include "SimTKcommon/internal/State.h"
#include "SimTKcommon/internal/System.h"

namespace SimTK {
class EventTriggerInfo; // defined below

//==============================================================================
//                             EVENT REPORTER
//==============================================================================
/** An %EventReporter is an object that defines an event that can occur within a
system. It is an abstract class. Subclasses define how to determine when the 
event occurs, and what happens when it does. You will not generally subclass 
%EventReporter directly. Instead, subclass ScheduledEventReporter (for events 
that occur at a particular time that you can calculate in advance) or 
TriggeredEventReporter (for events that occur when some condition is satisfied 
within the system state). ScheduledEventReporter also has another subclass, 
PeriodicEventReporter, for the common situation of events that occur at regular
intervals.

%EventReporter is very similar to EventHandler, but differs in that it is not 
permitted to modify the state of the system when the event occurs. It can 
observe the system's behavior, but not alter it.

Once you have created an %EventReporter, you can add it to a System by calling
System::adoptEventReporter(). **/
class SimTK_SimTKCOMMON_EXPORT EventReporter {
public:
    virtual ~EventReporter() {}

    /** Create a new, unowned copy of this %EventReporter on the heap. The
    System-related information will not be copied. **/
    EventReporter* clone() const {
        EventReporter* p = cloneVirtual();
        p->clearCache();
        return p;
    }

    /** Return the the description given at construction (if any) to be used 
    for the Event created by this %EventReporter. **/
    const std::string& getEventDescription() const {return m_eventDescription;}
    
    /** This method is invoked to report the event.  
    @param state
        The state of the system when the event occurred. Note that `state` is
        const so cannot be changed here. **/
    virtual void handleEvent(const State& state) const = 0;

    /** Return `true` if this %EventHandler has been adopted by a System. **/
    bool isInSystem() const {return !m_system.empty();}

    /** Get a reference to the System that owns this %EventReporter, if any.
    Otherwise throws an exception; use isInSystem() to check. **/
    const System& getSystem() const {
        SimTK_ASSERT_ALWAYS(isInSystem(), 
            "EventReporter::getSystem(): "
            "This EventReporter hasn't been adopted by a System.");
        return *m_system;
    }

    /** Return the Event created for this %EventReporter after it is adopted by
    a System. Throws an exception if this %EventReporter is unowned; use
    isInSystem() to check. **/
    const Event& getEvent() const 
    {   return getSystem().getEvent(m_eventId); }

    /** Return the EventTrigger created for this %EventReporter after it is 
    adopted by a System. This may be a Timer or Witness depending on the
    concrete type of %EventReporter. **/
    const EventTrigger& getEventTrigger() const 
    {   return getSystem().getEventTrigger(m_triggerId);}
protected:
    /** Constructor for use by concrete derived classes. **/
    explicit EventReporter(const std::string& eventDescription) 
    :   m_eventDescription(eventDescription) {}

    /** A concrete %EventReporter should implement this method to create an
    identical copy of itself on the heap and return the pointer. **/
    virtual EventReporter* cloneVirtual() const {
        SimTK_THROW2(Exception::UnimplementedVirtualMethod, "EventReporter", 
                     "cloneVirtual");
    }

    /** Get a writable reference to the System that owns this %EventReporter, 
    if any. Otherwise throws an exception; use isInSystem() to check. **/
    System& updSystem() {
        SimTK_ASSERT_ALWAYS(isInSystem(), 
            "EventReporter::updSystem(): "
            "This EventReporter hasn't been adopted by a System.");
        return *m_system;
    }

    /** Return a writable reference to the Event created for this %EventReporter 
    after it is adopted by a System. Throws an exception if this %EventReporter 
    is unowned; use isInSystem() to check. **/
    Event& updEvent() 
    {   return updSystem().updEvent(m_eventId); }

    /** Return a writable reference to the EventTrigger created for this 
    %EventReporter after it is adopted by a System. This may be a Timer or 
    Witness depending on the concrete type of %EventReporter. **/
    EventTrigger& updEventTrigger() 
    {   return updSystem().updEventTrigger(m_triggerId);}
private:
friend class SystemGlobalSubsystem;

    void clearCache() {
        m_system.reset();
        m_eventId.invalidate(); 
        m_triggerId.invalidate();
    }

    // Set during construction.
    std::string     m_eventDescription;

    // These are set when the EventReporter is adopted by a System. These 
    // should not be copied.
    ReferencePtr<System>    m_system;
    EventId                 m_eventId;
    EventTriggerId          m_triggerId;
};

//==============================================================================
//                         SCHEDULED EVENT REPORTER
//==============================================================================
/** %ScheduledEventReporter is a still-abstract subclass of EventHandler for 
events that occur at particular times that are known in advance. This includes 
events that occur just once, at multiple times, or periodically. The only 
requirement is that it must be able to calculate from the current time and 
state the next time at which the event should occur. **/
class SimTK_SimTKCOMMON_EXPORT ScheduledEventReporter : public EventReporter {
    using Super = EventReporter;
public:
    /** Invokes the base class clone() method but makes the return type more
    specific. **/
    ScheduledEventReporter* clone() const 
    {   return static_cast<ScheduledEventReporter*>(Super::clone()); }

    /** Determine the next time at which this event is scheduled to occur.    
    @param      state                
        The current time and state of the system.
    @param      includeCurrentTime    
        If true, return the next event whose time is >= the current time.
        If false, only return the next event *after* (not at) the current time.  
    @returns The time. 
    @bug Should be called `calcNextEventTime()`. **/    
    virtual Real getNextEventTime(const State& state, 
                                  bool includeCurrentTime) const = 0;

protected:
    explicit ScheduledEventReporter
       (const std::string& eventDescription = "ScheduledEventReporter Event")
    :   EventReporter(eventDescription) {}
};



//==============================================================================
//                           EVENT TRIGGER INFO
//==============================================================================
/** This class is used by TriggeredEventHandler and TriggeredEventReporter to
specify the properties of their event witness functions. This will be used to
construct an appropriate EventWitness object for the handler or reporter. Not 
all properties of an EventWitness can be set here; anything unspecified will 
have its default value. See EventWitness for more information.

The properties you can set here are:
  - Whether to watch for rising sign transitions, falling, or both. [BOTH]
  - The localization window in units of the System's timescale. [10%]
    (That is then the "unit" window which is reduced by the accuracy setting.)
    
The default values are shown in brackets above.
@see TriggeredEventHandler, TriggeredEventReporter **/
class EventTriggerInfo {
public:
    /** Construct default object; see class description for default values. **/
    EventTriggerInfo() {setDefaults();}

    /** Construct default object with a given EventId. **/
    explicit EventTriggerInfo(EventId eid) : m_eventId(eid) {setDefaults();}

    EventId getEventId() const {return m_eventId;}
    void setEventId(EventId eid) {m_eventId=eid;}

    bool shouldTriggerOnRisingSignTransition()  const
    {   return m_triggerOnRising; }

    bool shouldTriggerOnFallingSignTransition() const
    {   return m_triggerOnFalling; }

    Real getRequiredLocalizationTimeWindow() const
    {   return m_localizationWindow; }

    // These return the modified 'this', like assignment operators.
    EventTriggerInfo& setTriggerOnRisingSignTransition(bool shouldTrigger)
    {   m_triggerOnRising = shouldTrigger; return *this; }

    EventTriggerInfo& setTriggerOnFallingSignTransition(bool shouldTrigger)
    {   m_triggerOnFalling = shouldTrigger; return *this; }

    EventTriggerInfo& 
    setRequiredLocalizationTimeWindow(double fractionOfTimeScale) {
        SimTK_APIARGCHECK_ALWAYS(fractionOfTimeScale > 0, "EventTriggerInfo",
            "setRequiredLocalizationTimeWindow", 
            "Localization window (fraction of time scale) must be "
            "greater than zero");  
        m_localizationWindow = fractionOfTimeScale; 
        return *this;
    }


private:
    void setDefaults() {
        m_triggerOnRising = m_triggerOnFalling = true;
        m_localizationWindow = 0.1; // 10% of System timescale
    }

    EventId     m_eventId;
    bool        m_triggerOnRising;
    bool        m_triggerOnFalling;
    double      m_localizationWindow;
};


//==============================================================================
//                         TRIGGERED EVENT REPORTER
//==============================================================================
/** %TriggeredEventReporter is a still-abstract subclass of EventHandler for 
events that occur when some specific state-dependent condition occurs during
system time stepping. This is implemented by means of an "event witness 
function". For any State, the time stepper must be able to calculate the value 
of the function. When the witness passes through 0, that indicates the event has
occurred so the corresponding Event is triggered. You can also customize when 
the event is triggered, for example to specify that it is only triggered when 
the function goes from negative to positive, but not when it goes from positive
to negative. **/
class SimTK_SimTKCOMMON_EXPORT TriggeredEventReporter : public EventReporter {
    using Super = EventReporter;
public:  
    /** Invokes the base class clone() method but makes the return type more
    specific. **/
    TriggeredEventReporter* clone() const 
    {   return static_cast<TriggeredEventReporter*>(Super::clone()); }

    /** Calculate the value of the event witness function using the given
    `state`. This will be called during realization of the Stage that was 
    specified in the constructor, *after* all subsystems have been realized. 
    @bug Should be called `calcWitnessValue()`. **/  
    virtual Real getValue(const State& state) const = 0;
    
    /** Get a writable EventTriggerInfo object which can be used to customize 
    the treatment of the witness function. **/
    EventTriggerInfo& updTriggerInfo() {return m_triggerInfo;}

    /** Get a reference to the currently-set EventTriggerInfo for this
    triggered event. **/
    const EventTriggerInfo& getTriggerInfo() const {return m_triggerInfo;}
    
    /** Get the Stage at which the witness function will be evaluated (as 
    provided in the constructor). **/
    Stage getRequiredStage() const {return m_requiredStage;}

protected:
    /** Construct a new TriggeredEventReporter and provide the Stage at which
    its witness function can be evaluated.
    @param  requiredStage    
          The stage at which the witness function will be evaluated.
    @param  eventDescription
        An optional human-readable description of the Event created for this
        EventReporter. If you don't provide one, a generic one will be used. **/    
    explicit TriggeredEventReporter(Stage requiredStage,
        const std::string& eventDescription = "TriggeredEventReporter::Event") 
    :   EventReporter(eventDescription), m_requiredStage(requiredStage) {} 

private:
    EventTriggerInfo    m_triggerInfo;
    Stage               m_requiredStage;
};


//==============================================================================
//                         PERIODIC EVENT REPORTER
//==============================================================================
/** PeriodicEventReporter is a still-abstract subclass of ScheduledEventReporter
which generates a series of uniformly spaced events at regular intervals. This 
allows you to very easily create event reporters with this behavior. **/
class SimTK_SimTKCOMMON_EXPORT PeriodicEventReporter 
:   public ScheduledEventReporter {
    using Super = ScheduledEventReporter;
public:
    /** Invokes the base class clone() method but makes the return type more
    specific. **/
    PeriodicEventReporter* clone() const 
    {   return static_cast<PeriodicEventReporter*>(Super::clone()); }

    /** Create a PeriodicEventReporter and specify the time interval between
    event occurrences. A generic event description will be provided including
    "PeriodicEventReporter" and the period; see the other constructor if you 
    want to supply your own.
    @param      eventInterval
        The time interval at which events should occur. Must be strictly 
        greater than zero. **/
    explicit PeriodicEventReporter
       (Real               eventInterval) 
    :   ScheduledEventReporter("PeriodicEventReporter event, period=" 
                                                        + String(eventInterval))
    {   setEventInterval(eventInterval); }

    /** Create a PeriodicEventReporter, specify the time interval between
    event occurrences, and provide your own description.    
    @param      eventInterval
        The time interval at which events should occur. Must be strictly 
        greater than zero. 
    @param      eventDescription
        A human-readable description of the Event to be created for this 
        EventReporter. **/
   PeriodicEventReporter(Real               eventInterval, 
                         const std::string& eventDescription) 
    :   ScheduledEventReporter(eventDescription)
    {   setEventInterval(eventInterval); }

    /** Determine the next time at which this periodic event will occur.    
    @param      state                
        The current time and state of the system.
    @param      includeCurrentTime    
        If true, return the time of next occurrence for this event *at or after*
        the current time. If false, return the time *after* (not at) the current
        time.  
    @returns The time.  
    @bug Should be called `calcNextEventTime()`. **/    
    Real getNextEventTime(const State& state, 
                          bool includeCurrentTime) const override;
    
    /** Get the time interval at which events are to occur. This was specified
    in the constructor or setEventInterval(). **/   
    Real getEventInterval() const {return m_eventInterval;}
    
    /** Change the time interval at which events are to occur.   
    @param      eventInterval
        The time interval at which events should occur. Must be strictly 
        greater than zero. **/
    void setEventInterval(Real eventInterval);

private:
    Real                m_eventInterval;
};


} // namespace SimTK

#endif // SimTK_SimTKCOMMON_EVENT_REPORTER_H_
