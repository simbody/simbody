#ifndef SimTK_SimTKCOMMON_EVENT_REPORTER_H_
#define SimTK_SimTKCOMMON_EVENT_REPORTER_H_

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

#include "SimTKcommon.h"

namespace SimTK {

/**
 * An EventReporter is an object that defines an event that can occur within a system.
 * It is an abstract class.  Subclasses define how to determine when the event occurs,
 * and what happens when it does.  You will not generally subclass EventReporter
 * directly.  Instead, subclass ScheduledEventReporter (for events that occur at a
 * particular time that is know in advance) or TriggeredEventReporter (for events that
 * occur when some condition is satisfied within the system).
 * 
 * EventReporter is very similar to EventHandler, but differs in that it is not permitted
 * to modify the state of the system.  It can observe the system's behavior, but not
 * alter it.  This means that adding an EventReporter to a System is not considered a
 * change to the physical system it represents. 
 * 
 * Once you have created an EventReporter, you can add it to a System by calling
 * getDefaultSubsystem().addEventReporter() on the System.
 */

class SimTK_SimTKCOMMON_EXPORT EventReporter {
public:
    
    /**
     * This method is invoked to handle the event.  It is given a State which describes the
     * system at the time when the event occurs.
     */
    
    virtual void handleEvent(const State& state) = 0;
};

/**
 * ScheduledEventReporter is a subclass of EventReporter for events that occur at a particular
 * time that is known in advance.  This includes events that occur multiple times.  The only
 * requirement is that, at any time, it must be able to report the next time at which an
 * event will occur.
 */

class SimTK_SimTKCOMMON_EXPORT ScheduledEventReporter : public EventReporter {
public:
    
    /**
     * Get the next time at which an event will occur.
     */
    
    virtual Real getNextEventTime(const State&) const = 0;
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
    class TriggeredEventReporterRep;
    
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
    
    System::EventTriggerInfo& getTriggerInfo();
    
    /**
     * Get the stage at which the trigger function will be evaluated.
     */
    
    Stage getRequiredStage() const;
private:
    TriggeredEventReporterRep* rep;
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_EVENT_REPORTER_H_
