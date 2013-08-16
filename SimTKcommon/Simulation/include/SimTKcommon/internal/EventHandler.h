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
 * Portions copyright (c) 2007-12 Stanford University and the Authors.        *
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

/** An EventHandler is an object that defines an event that can occur within a 
system. It is an abstract class.  Subclasses define how to determine when the 
event occurs, and what happens when it does.  You will not generally subclass 
EventHandler directly.  Instead, subclass ScheduledEventHandler (for events 
that occur at a particular time that is know in advance) or 
TriggeredEventHandler (for events that occur when some condition is satisfied 
within the system).  ScheduledEventHandler also has another subclass, 
PeriodicEventHandler, for the common situation of events that occur at regular 
intervals.

An EventHandler should be thought of as an integral part of the system it 
belongs to, and may alter the physical properties or behavior of the system.  
If you merely want to observe the system but not to alter it, you should 
generally use a EventReporter instead.

Once you have created an EventHandler, you can add it to a System by calling
addEventHandler() on the System. **/
class SimTK_SimTKCOMMON_EXPORT EventHandler {
public:
    virtual ~EventHandler();
    
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

/**
 * ScheduledEventHandler is a subclass of EventHandler for events that occur at a particular
 * time that is known in advance.  This includes events that occur multiple times.  The only
 * requirement is that, at any time, it must be able to report the next time at which an
 * event will occur.
 */

class SimTK_SimTKCOMMON_EXPORT ScheduledEventHandler : public EventHandler {
public:
    virtual ~ScheduledEventHandler();
    
    /**
     * Get the next time at which an event will occur.
     * 
     * @param state                 the current state of the system
     * @param includeCurrentTime    if true, return the next event whose time is >= the current time.
     *                              If false, only return events after (not at) the current time.
     */
    
    virtual Real getNextEventTime(const State& state, 
                                  bool includeCurrentTime) const = 0;
};

/**
 * TriggeredEventHandler is a subclass of EventHandler for events that occur when some condition
 * is satisfied within the system.  This is implemented by means of an "event trigger function".
 * For any State, the handler must be able to calculate the value of the function.  When it
 * passes through 0, that indicates the event has occurred.  You can also customize when the
 * event is triggered, for example to specify that it is only triggered when the function
 * goes from negative to positive, not when it goes from positive to negative.
 */

class SimTK_SimTKCOMMON_EXPORT TriggeredEventHandler : public EventHandler {
public:
    class TriggeredEventHandlerImpl;
    TriggeredEventHandler(const TriggeredEventHandler& clone);
    TriggeredEventHandler& operator=(const TriggeredEventHandler& clone);
    virtual ~TriggeredEventHandler();
    
    /**
     * Construct a new TriggeredEventHandler.
     * 
     * @param   requiredStage    
     *      The stage at which the trigger function will be evaluated.
     */    
    TriggeredEventHandler(Stage requiredStage);
    
    /**
     * Get the value of the event trigger function for a State.
     */  
    virtual Real getValue(const State&) const = 0;
    
    /**
     * Get an EventTriggerInfo object which can be used to customize when the
     * event occurs.
     */  
    EventTriggerInfo& getTriggerInfo();
    
    /**
     * Get the stage at which the trigger function will be evaluated.
     */  
    Stage getRequiredStage() const;
private:
    TriggeredEventHandlerImpl* impl;
};

/** PeriodicEventHandler is a subclass of ScheduledEventHandler which generates
a series of uniformly spaced events at regular intervals. This allows you to 
very easily create event handlers with this behavior. **/
class SimTK_SimTKCOMMON_EXPORT PeriodicEventHandler 
:   public ScheduledEventHandler {
public:
    class PeriodicEventHandlerImpl;
    ~PeriodicEventHandler();
    Real getNextEventTime(const State& state, bool includeCurrentTime) const;
    
    /**
     * Create a PeriodicEventHandler.
     * 
     * @param   eventInterval
     *     The time interval at which events should occur.
     */
    PeriodicEventHandler(Real eventInterval);
    
    /**
     * Get the time interval at which events occur.
     */   
    Real getEventInterval() const;
    
    /**
     * Set the time interval at which events occur.
     */   
    void setEventInterval(Real eventInterval);
private:
    PeriodicEventHandlerImpl* impl;
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_EVENT_HANDLER_H_
