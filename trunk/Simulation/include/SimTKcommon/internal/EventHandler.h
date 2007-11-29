#ifndef SimTK_SimTKCOMMON_EVENT_HANDLER_H_
#define SimTK_SimTKCOMMON_EVENT_HANDLER_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"
#include "SimTKcommon/internal/State.h"
#include "SimTKcommon/internal/System.h"

namespace SimTK {

/**
 * An EventHandler is an object that defines an event that can occur within a system.
 * It is an abstract class.  Subclasses define how to determine when the event occurs,
 * and what happens when it does.  You will not generally subclass EventHandler
 * directly.  Instead, subclass ScheduledEventHandler (for events that occur at a
 * particular time that is know in advance) or TriggeredEventHandler (for events that
 * occur when some condition is satisfied within the system).  ScheduledEventHandler
 * also has another subclass, PeriodicEventHandler, for the common situation of
 * events that occur at regular intervals.
 * 
 * An EventHandler should be thought of as an integral part of the system it belongs
 * to, and may alter the physical properties or behavior of the system.  If you merely
 * want to observe the system but not to alter it, you should generally use a
 * EventReporter instead.
 * 
 * Once you have created an EventHandler, you can add it to a System by calling
 * updDefaultSubsystem().addEventHandler() on the System.
 */

class SimTK_SimTKCOMMON_EXPORT EventHandler {
public:
    
    /**
     * This method is invoked to handle the event.  It is given a State which describes the
     * system at the time when the event occurs, and it is permitted to modify any aspect of
     * the state except the time.  In doing so, it should respect the specified accuracy
     * requirements for the continuous variables and constraints.
     * 
     * @param state             the state of the system when the event occurred.  This method should
     *                          modify it to reflect the changes caused by the event.
     * @param accuracy          the overall accuracy for the simulation.  This acts as a multiplier for
     *                          the weights and constraint tolerances.
     * @param yWeights          a vector of weights for the continuous state variables.  If this method
     *                          modifies a variable, its new value should be accurate to within accuracy*yWeight.
     * @param ooConstraintTols  a vector of constraint tolerances.  If this method modifies the
     *                          state, each constraint should be satisfied to within accuracy*ooConstraintTol.
     * @param lowestModified    if this method modifies the state, it should set this to the lowest stage
     *                          which it modified, so that the realization cache can be updated.
     * @param shouldTerminate   if the event handler sets this to true, it will cause the simulation to terminate
     *                          immediately.
     */
    
    virtual void handleEvent(State& state, Real accuracy, const Vector& yWeights, const Vector& ooConstraintTols, Stage& lowestModified, bool& shouldTerminate) = 0;
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
    
    virtual Real getNextEventTime(const State&, bool includeCurrentTime) const = 0;
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
     * @param requiredStage    the stage at which the trigger function will be evaluated
     */
    
    TriggeredEventHandler(Stage requiredStage);
    
    /**
     * Get the value of the event trigger function for a State.
     */
    
    virtual Real getValue(const State&) const = 0;
    
    /**
     * Get an EventTriggerInfo object which can be used to customize when the event occurs.
     */
    
    System::EventTriggerInfo& getTriggerInfo();
    
    /**
     * Get the stage at which the trigger function will be evaluated.
     */
    
    Stage getRequiredStage() const;
private:
    TriggeredEventHandlerImpl* impl;
};

/**
 * PeriodicEventHandler is a subclass of ScheduledEventHandler which generates a series
 * of uniformly spaced events at regular intervals.  This allows you to very easily create
 * event handlers with this behavior.
 */

class SimTK_SimTKCOMMON_EXPORT PeriodicEventHandler : public ScheduledEventHandler {
public:
    class PeriodicEventHandlerImpl;
    ~PeriodicEventHandler();
    Real getNextEventTime(const State& state, bool includeCurrentTime) const;
    
    /**
     * Create a PeriodicEventHandler.
     * 
     * @param eventInterval       the time interval at which events should occur.
     */

    PeriodicEventHandler(Real eventInterval);
    
    /**
     * Get the time interval at which events occur.
     */
    
    Real getEventInterval();
    
    /**
     * Set the time interval at which events occur.
     */
    
    void setEventInterval(Real eventInterval);
private:
    PeriodicEventHandlerImpl* impl;
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_EVENT_HANDLER_H_
