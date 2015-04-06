#ifndef SimTK_SimTKCOMMON_EVENT_HANDLER_H_
#define SimTK_SimTKCOMMON_EVENT_HANDLER_H_

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

//==============================================================================
//                             EVENT HANDLER
//==============================================================================
/** An %EventHandler is an object that defines an event that can occur within a 
system. It is an abstract class. Subclasses define how to determine when the 
event occurs, and what happens when it does. You will not generally subclass 
%EventHandler directly. Instead, subclass ScheduledEventHandler (for events 
that occur at a particular time that you can calculate in advance) or 
TriggeredEventHandler (for events that occur when some condition is satisfied 
within the system state). ScheduledEventHandler also has another subclass, 
PeriodicEventHandler, for the common situation of events that occur at regular 
intervals.

An %EventHandler should be thought of as an integral part of the system it 
belongs to, and may alter the physical properties or behavior of the system.  
If you merely want to observe the system upon even occurrence, and not alter 
its behavior in any way, you should use a EventReporter instead.

Once you have created an %EventHandler, you can add it to a System by calling
System::addEventHandler(). **/
class SimTK_SimTKCOMMON_EXPORT EventHandler {
public:
    virtual ~EventHandler() {}
    
    /** This method is invoked to handle the event. It is given a State which 
    describes the system at the time when the event occurs, and it is permitted
    to modify any aspect of the state except the time. In doing so, it should 
    respect the specified accuracy requirements for the continuous variables 
    and constraints.
     
    @param state
        The state of the system when the event occurred. This method should
        modify \a state to reflect the changes caused by the event.
    @param accuracy          
        The accuracy to which this simulation is being computed. If your 
        handler performs any approximate operation it should do so consistent
        with the simulation accuracy.
    @param shouldTerminate   
        If the event handler sets this to true, it will cause the simulation to
        terminate immediately. This does not necessarily indicate an error
        condition.
    **/   
    virtual void handleEvent(State& state, Real accuracy, 
                             bool& shouldTerminate) const = 0;
};


//==============================================================================
//                         SCHEDULED EVENT HANDLER
//==============================================================================
/** %ScheduledEventHandler is a subclass of EventHandler for events that occur
at a particular time that is known in advance. This includes events that occur 
multiple times. The only requirement is that it must be able to calculate from 
the current time and state the next time at which the event should occur. **/
class SimTK_SimTKCOMMON_EXPORT ScheduledEventHandler : public EventHandler {
public:
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
};


//==============================================================================
//                         TRIGGERED EVENT HANDLER
//==============================================================================
/** %TriggeredEventHandler is a subclass of EventHandler for events that occur
when some condition is satisfied within the system. This is implemented by means
of an "event witness function". For any State, the handler must be able to 
calculate the value of the function. When it passes through 0, that indicates 
the event has occurred. You can also customize when the event is triggered, for
example to specify that it is only triggered when the function goes from 
negative to positive, not when it goes from positive to negative. **/
class SimTK_SimTKCOMMON_EXPORT TriggeredEventHandler : public EventHandler {
public:
    /** Construct a new TriggeredEventHandler and provide the Stage at which
    its witness function can be evaluated.
    @param  requiredStage    
          The stage at which the witness function will be evaluated. **/    
    explicit TriggeredEventHandler(Stage requiredStage) 
    :   m_requiredStage(requiredStage) {}
    
    /** Calculate the value of the event witness function using the given
    `state`. This will be called during realization of the Stage that was 
    specified in the constructor, *after* all subsystems have been realized. 
    @bug Should be called `calcValue()`. **/  
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

private:
    EventTriggerInfo    m_triggerInfo;
    Stage               m_requiredStage;
};


//==============================================================================
//                         PERIODIC EVENT HANDLER
//==============================================================================
/** PeriodicEventHandler is a subclass of ScheduledEventHandler which generates
a series of uniformly spaced events at regular intervals. This allows you to 
very easily create event handlers with this behavior. **/
class SimTK_SimTKCOMMON_EXPORT PeriodicEventHandler 
:   public ScheduledEventHandler {
public:
    /** Create a PeriodicEventHandler and specify the time interval between
    event occurrences.    
    @param      eventInterval
        The time interval at which events should occur. Must be strictly 
        greater than zero. **/
    explicit PeriodicEventHandler(Real eventInterval) 
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

#endif // SimTK_SimTKCOMMON_EVENT_HANDLER_H_
