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
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
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

namespace SimTK {
    
/** @class SimTK::EventId
This is a class to represent unique IDs for events in a type-safe way. These
are created and managed by the System, not the State. **/
SimTK_DEFINE_UNIQUE_INDEX_TYPE(EventId);

/** @class SimTK::SystemEventTriggerIndex
This unique integer type is for identifying a triggered event in the full 
System-level view of the State. More precisely, this is the index of the slot 
in the global array in the cache allocated to hold the value of that event's
trigger function.
@see EventTriggerIndex **/
SimTK_DEFINE_UNIQUE_INDEX_TYPE(SystemEventTriggerIndex);

/** @class SimTK::SystemEventTriggerByStageIndex
This unique integer type is for identifying a triggered event within a 
particular Stage of the full System-level view of the State. (Event triggers 
for a particular Stage are stored consecutively within the full collection of 
event triggers.) That is, the EventTriggerByStageIndex will be 0 for the first 
event trigger at that stage.
@see EventTriggerByStageIndex
**/
SimTK_DEFINE_UNIQUE_INDEX_TYPE(SystemEventTriggerByStageIndex);

/** @class SimTK::EventTriggerByStageIndex
Unique integer type for Subsystem-local, per-stage event indexing.
@see SystemEventTriggerByStageIndex **/
SimTK_DEFINE_UNIQUE_INDEX_TYPE(EventTriggerByStageIndex);

/** An Event is "something that happens" during a Study that is advancing
through time. Its occurence interrupts the normal flow of computation, allowing
an event Handler to adjust the State prior to resuming the Study.

Events are allocated by Subsystems, but require some System global resources. 
All Events are given a unique SystemEventId. Some Events require other State 
resources, such as slots for the values of trigger functions in the case of 
Triggered events.

Events can be allocated at Topology, Model, and Instance Stages. All Event 
resources are assigned when the Instance stage is realized. However, if an 
Event requires state variables, then it must be allocated by Model stage. **/
class Event {
public:

    /** These are all the possible causes for events.

    @par Initialization
    A Study has performed its own initialization and is about to start.
    @par Triggered
    An event trigger function underwent a monitored sign transition.
    @par Scheduled
    An integrator reached a previously-scheduled time for the Event to occur.
    @par TimeAdvanced
    An integrator completed an internal step, meaning that it has reached
    a point where time has advanced irreversibly.
    @par Signaled
    A flag in the State has been explicitly set, meaning that a particular
    Event has occurred. Anyone with write access to a State can set these,
    but typically they are set in event handlers associated with one of
    the other kinds of events.
    @par Termination
    The Study has finished. If a Termination event handler signals more
    Events, those signaled events are not processed by the Study; that is,
    the signals remain set in the final State.

    In case several of these causes are detected in 
    a single step, they are sequentialized in the order shown,
    like this:
       1. The occurrence of triggered events is reported and
          the triggering state and a list of triggered events
          are passed to the event handler for processing (meaning
          the state, but not the time, is modified). [Note that
          simultaneity *within* the set of triggered events may
          also require special handling; we're not talking about
          that here, just simultaneity of *causes*.]
       2. Next, using the state resulting from step 1, the time is checked
          to see if scheduled events have occurred. If so, a list of
          those events is passed to the event handler for processing.
       3. Next, if this system has requested time-advanced events,
          the event handler is called with the state that resulted
          from step 2 and the "time advanced" cause noted. No event
          list is passed in that case. The state may be modified.
       4. Last, if the final time has been reached or if any of
          the event handlers asked for termination, we pass the
          state to the event handler again noting that we have
          reached termination. The state may be modified and the
          result will be the final state of the simulation.
    **/
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

        bool isValid() const {return Initialization<=value && value<=Termination;}

    private:
        Num value;
    };

    /** This is useful for debugging; it translates an Event::Cause into a 
    readable string. **/
    SimTK_SimTKCOMMON_EXPORT static const char* getCauseName(Cause);


    /** Triggered Events respond to zero crossings of their associated trigger
    function. This enum defines constants for use in specifying which kind
    of zero crossing has been seen, or which kinds are considered interesting.
    For the latter purpose, these can be or'ed together to make a mask. **/
    enum Trigger {
        NoEventTrigger          =0x0000,    // must be 0

        PositiveToNegative      =0x0001,    // 1
        NegativeToPositive      =0x0002,    // 2

        Falling                 =(PositiveToNegative), // 1
        Rising                  =(NegativeToPositive), // 2
        AnySignChange           =(PositiveToNegative|NegativeToPositive)    // 3
    };

    /** This is useful for debugging; it translates an Event::Trigger or a mask
    formed by a union of Event::Triggers, into a readable string. **/
    SimTK_SimTKCOMMON_EXPORT static std::string eventTriggerString(Trigger);


    /** Classify a before/after sign transition. Before and after must both
    be -1,0, or 1 as returned by the SimTK::sign() function applied to
    the trigger function value at the beginning and end of a step. **/
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

    /** Given an observed transition, weed out ignorable ones using the supplied
    mask. That is, the return will indicate NoEventTrigger unless the original
    Trigger was present in the mask. **/
    static Trigger maskTransition(Trigger transition, Trigger mask) {
        // we're depending on NoEventTrigger==0
        return Trigger(transition & mask); 
    }

private:
};


/** This class is used to communicate between the System and an Integrator 
regarding the properties of a particular event trigger function. Currently 
these are:
  - Whether to watch for rising sign transitions, falling, or both. [BOTH]
  - Whether to watch for transitions to and from zero. [NO]
  - The localization window in units of the System's timescale. [10%]
    (That is then the "unit" window which is reduced by the accuracy setting.)
    
The default values are shown in brackets above. **/
class SimTK_SimTKCOMMON_EXPORT EventTriggerInfo {
public:
    EventTriggerInfo();
    explicit EventTriggerInfo(EventId eventId);
    ~EventTriggerInfo();
    EventTriggerInfo(const EventTriggerInfo&);
    EventTriggerInfo& operator=(const EventTriggerInfo&);

    EventId getEventId() const; // returns -1 if not set
    bool shouldTriggerOnRisingSignTransition()  const; // default=true
    bool shouldTriggerOnFallingSignTransition() const; // default=true
    Real getRequiredLocalizationTimeWindow()    const; // default=0.1

    // These return the modified 'this', like assignment operators.
    EventTriggerInfo& setEventId(EventId);
    EventTriggerInfo& setTriggerOnRisingSignTransition(bool);
    EventTriggerInfo& setTriggerOnFallingSignTransition(bool);
    EventTriggerInfo& setRequiredLocalizationTimeWindow(Real);

    Event::Trigger calcTransitionMask() const {
        unsigned mask = 0;
        if (shouldTriggerOnRisingSignTransition()) {
            mask |= Event::NegativeToPositive;
        }
        if (shouldTriggerOnFallingSignTransition()) {
            mask |= Event::PositiveToNegative;
        }
        return Event::Trigger(mask);
    }

    Event::Trigger calcTransitionToReport
       (Event::Trigger transitionSeen) const
    {
        // report -1 to 1 or 1 to -1 as appropriate
        if (transitionSeen & Event::Rising)
            return Event::NegativeToPositive;
        if (transitionSeen & Event::Falling)
            return Event::PositiveToNegative;
        assert(!"impossible event transition situation");
        return Event::NoEventTrigger;
    }

private:
    class EventTriggerInfoRep;

    // opaque implementation for binary compatibility
    EventTriggerInfoRep* rep;

    const EventTriggerInfoRep& getRep() const {assert(rep); return *rep;}
    EventTriggerInfoRep&       updRep()       {assert(rep); return *rep;}
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
    HandleEventsResults() : m_lowestModifiedStage(Stage::Infinity) {clear();}

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
