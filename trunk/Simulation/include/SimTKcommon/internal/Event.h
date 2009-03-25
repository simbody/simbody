#ifndef SimTK_SimTKCOMMON_EVENT_H_
#define SimTK_SimTKCOMMON_EVENT_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-9 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
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

/** @file
 *
 * This file declares the types needed for SimTK's support for Events.
 */

#include "SimTKcommon/basics.h"

namespace SimTK {

/// This unique integer type is for identifying an event in the full System-level
/// view of the State. This is a global resource and applies to all System-level
/// events, not just triggered events requiring other State resources.
/// @see EventIndex for Subsystem-local event indexing
SimTK_DEFINE_UNIQUE_INDEX_TYPE(SystemEventIndex);
/// Unique integer type for Subsystem-local event indexing
/// @see SystemEventIndex
SimTK_DEFINE_UNIQUE_INDEX_TYPE(EventIndex);

/// This unique integer type is for identifying a triggered event in the full System-level
/// view of the State. More precisely, this is the index of the slot in
/// the global array in the cache allocated to hold the value of that event's
/// trigger function.
/// @see EventTriggerIndex for Subsystem-local event indexing
SimTK_DEFINE_UNIQUE_INDEX_TYPE(SystemEventTriggerIndex);
/// Unique integer type for Subsystem-local event indexing
/// @see SystemEventTriggerIndex
SimTK_DEFINE_UNIQUE_INDEX_TYPE(EventTriggerIndex);

/// This unique integer type is for identifying a triggered event within a particular Stage
/// of the full System-level view of the State. (Event triggers for a particular Stage
/// are stored consecutively within the full collection of event triggers.) That is,
/// the EventTriggerByStageIndex will be 0 for the first event trigger at that stage.
/// @see SystemEventTriggerIndex for System-global event indexing
/// @see EventTriggerByStageIndex for Subsystem-local, per-stage event indexing
SimTK_DEFINE_UNIQUE_INDEX_TYPE(SystemEventTriggerByStageIndex);
/// Unique integer type for Subsystem-local, per-stage event indexing
/// @see SystemEventTriggerIndex
SimTK_DEFINE_UNIQUE_INDEX_TYPE(EventTriggerByStageIndex);

/**
 * An Event is "something that happens" during a Study that is advancing
 * through time. Its occurence interrupts the normal flow of computation, allowing
 * an event Handler to adjust the State prior to resuming the Study.
 *
 * Events are allocated by Subsystems, but require some System global 
 * resources. All Events are given a unique SystemEventId. Some Events
 * require other State resources, such as slots for the values of trigger
 * functions in the case of Triggered events.
 *
 * Events can be allocated at Topology, Model, and Instance Stages. All
 * Event resources are assigned when the Instance stage is realized. However, if
 * an Event requires state variables, then it must be allocated by Model
 * stage.
 */
class Event {
public:

    /**
     * These are all the possible causes for events.
     *
     * @par Initialization
     * A Study has performed its own initialization and is about to start.
     * @par Triggered
     * An event trigger function underwent a monitored sign transition.
     * @par Scheduled
     * An integrator reached a previously-scheduled time for the Event to occur.
     * @par TimeAdvanced
     * An integrator completed an internal step, meaning that it has reached
     * a point where time has advanced irreversibly.
     * @par Signaled
     * A flag in the State has been explicitly set, meaning that a particular
     * Event has occurred. Anyone with write access to a State can set these,
     * but typically they are set in event handlers associated with one of
     * the other kinds of events.
     * @par Termination
     * The Study has finished. If a Termination event handler signals more
     * Events, those signaled events are not processed by the Study; that is,
     * the signals remain set in the final State.
     *
     * In case several of these causes are detected in 
     * a single step, they are sequentialized in the order shown,
     * like this:
     *    1. The occurrence of triggered events is reported and
     *       the triggering state and a list of triggered events
     *       are passed to the event handler for processing (meaning
     *       the state, but not the time, is modified). [Note that
     *       simultaneity *within* the set of triggered events may
     *       also require special handling; we're not talking about
     *       that here, just simultaneity of *causes*.]
     *    2. Next, using the state resulting from step 1, the time is checked
     *       to see if scheduled events have occurred. If so, a list of
     *       those events is passed to the event handler for processing.
     *    3. Next, if this system has requested time-advanced events,
     *       the event handler is called with the state that resulted
     *       from step 2 and the "time advanced" cause noted. No event
     *       list is passed in that case. The state may be modified.
     *    4. Last, if the final time has been reached or if any of
     *       the event handlers asked for termination, we pass the
     *       state to the event handler again noting that we have
     *       reached termination. The state may be modified and the
     *       result will be the final state of the simulation.
     */
    class Cause {
    public:
        enum Num {
            Initialization  = 1,
            Triggered       = 2,
            Scheduled       = 3,
            TimeAdvanced    = 4,
            Signaled        = 5,
            Termination     = 6,
            Invalid         = -1
        };

        Cause() : value(Invalid) {}
        Cause(Num n) : value(n) {}           // implicit conversion
        operator Num() const {return value;} // implicit conversion
        Cause& operator=(Num n) {value=n; return *this;}

        bool isValid() const {Initialization<=value && value<=Termination;}

    private:
        Num value;
    };

    /// This is useful for debugging; it translates an Event::Cause
    /// into a readable string.
    SimTK_SimTKCOMMON_EXPORT static const char* getCauseName(Cause);


    /**
     * Triggered Events respond to zero crossings of their associated trigger
     * function. This enum defines constants for use in specifying which kind
     * of zero crossing has been seen, or which kinds are considered interesting.
     * For the latter purpose, these can be or'ed together to make a mask.
     */
    enum Trigger {
        NoEventTrigger          =0x0000,    // must be 0

        PositiveToNegative      =0x0001,    // 1
        NegativeToPositive      =0x0002,    // 2

        Falling                 =(PositiveToNegative), // 1
        Rising                  =(NegativeToPositive), // 2
        AnySignChange           =(PositiveToNegative|NegativeToPositive)    // 3
    };

    /// This is useful for debugging; it translates an Event::Trigger or a
    /// mask formed by a union of Event::Triggers, into a readable string.
    SimTK_SimTKCOMMON_EXPORT static std::string eventTriggerString(Trigger);


    /**
     * Classify a before/after sign transition. Before and after must both
     * be -1,0, or 1 as returned by the SimTK::sign() function applied to
     * the trigger function value at the beginning and end of a step.
     */
    static Trigger classifyTransition(int before, int after) {
        if (before==after)
            return NoEventTrigger;
        if (before==0)
            return NoEventTrigger; // Do not report transitions away from zero.
        if (before==1)
            return PositiveToNegative;
        // before==-1
        return NegativeToPositive;
    }

    /**
     * Given an observed transition, weed out ignorable ones using the supplied mask.
     * That is, the return will indicate NoEventTrigger unless the original Trigger
     * was present in the mask.
     */
    static Trigger maskTransition(Trigger transition, Trigger mask) {
        return Trigger(transition & mask); // we're depending on NoEventTrigger==0
    }

private:
};

/*
class EventStatus {
public:
    EventStatus() { initialize(); }
    // default destructor, copy constructor, copy assignment

    // Event trigger (which zero crossings cause triggering). Can be
    // OR'ed together to make a mask.
    enum EventTrigger {
        NoEventTrigger          =0x0000,    // must be 0

        PositiveToNegative      =0x0001,    // 1
        NegativeToPositive      =0x0002,    // 2

        Falling                 =(PositiveToNegative), // 1
        Rising                  =(NegativeToPositive), // 2
        AnySignChange           =(PositiveToNegative|NegativeToPositive)    // 3
    };

    bool isEventPending() const {return transitionSeen != NoEventTrigger;}
    EventTrigger getEventTrigger() const {return transitionSeen;}
    Real getLastTriggerTime() const {return lastTriggerTime;}
    Real getLastTriggerTimeBestGuess() const {return lastTriggerTimeBestGuess;}
    Real getBeforeValue() const {return beforeValue;}
    Real getAfterValue() const {return afterValue;}
    Real getLocalizationWindow() const {return localizationWindow;}

    void setEventTriggered(EventTrigger transition, Real triggerTime,
                           Real actualTimeEst, Real window,
                           Real before, Real after)
    {
        assert(transition != NoEventTrigger);
        assert(triggerTime >= 0 && actualTimeEst >= 0 
               && triggerTime >= actualTimeEst);

        transitionSeen = transition;
        lastTriggerTime = triggerTime;
        lastTriggerTimeBestGuess = actualTimeEst;
        localizationWindow = window;
        beforeValue = before;
        afterValue  = after;
    }

    void clearEventTrigger() {
        transitionSeen = NoEventTrigger;
    }

    // Classify a before/after sign transition.
    static EventTrigger classifyTransition(int before, int after) {
        if (before==after)
            return NoEventTrigger;
        if (before==0)
            return NoEventTrigger; // Do not report transitions away from zero.
        if (before==1)
            return PositiveToNegative;
        // before==-1
        return NegativeToPositive;
    }

    static EventTrigger maskTransition(EventTrigger transition, EventTrigger mask) {
        return EventTrigger(transition & mask); // we're depending on NoEventTrigger==0
    }

    SimTK_SimTKCOMMON_EXPORT static String eventTriggerString(EventTrigger e);
private:
    void initialize() {
        transitionSeen = NoEventTrigger;
        lastTriggerTime = lastTriggerTimeBestGuess = localizationWindow
            = beforeValue = afterValue = NaN;
    }

    EventTrigger transitionSeen;
    Real         lastTriggerTime; // digital
    Real         lastTriggerTimeBestGuess; // analog, <=lastTriggerTime
    Real         localizationWindow;
    Real         beforeValue, afterValue;
};
*/

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_EVENT_H_
