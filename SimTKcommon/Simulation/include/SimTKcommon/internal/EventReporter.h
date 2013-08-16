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
 * Portions copyright (c) 2007-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
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
#include "SimTKcommon/internal/System.h"

namespace SimTK {

/**
 * An EventReporter is an object that defines an event that can occur within a system.
 * It is an abstract class.  Subclasses define how to determine when the event occurs,
 * and what happens when it does.  You will not generally subclass EventReporter
 * directly.  Instead, subclass ScheduledEventReporter (for events that occur at a
 * particular time that is know in advance) or TriggeredEventReporter (for events that
 * occur when some condition is satisfied within the system).  ScheduledEventReporter
 * also has another subclass, PeriodicEventReporter, for the common situation of
 * events that occur at regular intervals.
 * 
 * EventReporter is very similar to EventHandler, but differs in that it is not permitted
 * to modify the state of the system.  It can observe the system's behavior, but not
 * alter it.  This means that adding an EventReporter to a System is not considered a
 * change to the physical system it represents. 
 * 
 * Once you have created an EventReporter, you can add it to a System by calling
 * addEventReporter() on the System.
 */

class SimTK_SimTKCOMMON_EXPORT EventReporter {
public:
    virtual ~EventReporter();
    
    /**
     * This method is invoked to handle the event.  It is given a State which describes the
     * system at the time when the event occurs.
     */
    
    virtual void handleEvent(const State& state) const = 0;
};

/**
 * ScheduledEventReporter is a subclass of EventReporter for events that occur at a particular
 * time that is known in advance.  This includes events that occur multiple times.  The only
 * requirement is that, at any time, it must be able to report the next time at which an
 * event will occur.
 */

class SimTK_SimTKCOMMON_EXPORT ScheduledEventReporter : public EventReporter {
public:
    virtual ~ScheduledEventReporter();
    
    /**
     * Get the next time at which an event will occur.
     * 
     * @param state                 the current state of the system
     * @param includeCurrentTime    if true, return the next event whose time is >= the current time.
     *                              If false, only return events after (not at) the current time.
     */
    
    virtual Real getNextEventTime(const State& state, bool includeCurrentTime) const = 0;
};

/**
 * TriggeredEventReporter is a subclass of EventReporter for events that occur when some condition
 * is satisfied within the system.  This is implemented by means of an "event trigger function".
 * For any State, the handler must be able to calculate the value of the function.  When it
 * passes through 0, that indicates the event has occurred.  You can also customize when the
 * event is triggered, for example to specify that it is only triggered when the function
 * goes from negative to positive, not when it goes from positive to negative.
 */

class SimTK_SimTKCOMMON_EXPORT TriggeredEventReporter : public EventReporter {
public:
    class TriggeredEventReporterImpl;
    TriggeredEventReporter(const TriggeredEventReporter& clone);
    TriggeredEventReporter& operator=(const TriggeredEventReporter& clone);
    virtual ~TriggeredEventReporter();
    
    /**
     * Construct a new TriggeredEventReporter.
     * 
     * @param requiredStage    the stage at which the trigger function will be evaluated
     */
    
    TriggeredEventReporter(Stage requiredStage);
    
    /**
     * Get the value of the event trigger function for a State.
     */
    
    virtual Real getValue(const State&) const = 0;
    
    /**
     * Get an EventTriggerInfo object which can be used to customize when the event occurs.
     */
    
    EventTriggerInfo& getTriggerInfo();
    
    /**
     * Get the stage at which the trigger function will be evaluated.
     */
    
    Stage getRequiredStage() const;
private:
    TriggeredEventReporterImpl* impl;
};

/**
 * PeriodicEventReporter is a subclass of ScheduledEventReporter which generates a series
 * of uniformly spaced events at regular intervals.  This allows you to very easily create
 * event reporters with this behavior.
 */

class SimTK_SimTKCOMMON_EXPORT PeriodicEventReporter : public ScheduledEventReporter {
public:
    class PeriodicEventReporterImpl;
    ~PeriodicEventReporter();
    Real getNextEventTime(const State& state, bool includeCurrentTime) const;
    
    /**
     * Create a PeriodicEventReporter.
     * 
     * @param eventInterval       the time interval at which events should occur.
     */

    PeriodicEventReporter(Real eventInterval);
    
    /**
     * Get the time interval at which events occur.
     */
    
    Real getEventInterval() const;
    
    /**
     * Set the time interval at which events occur.
     */
    
    void setEventInterval(Real eventInterval);
private:
    PeriodicEventReporterImpl* impl;
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_EVENT_REPORTER_H_
