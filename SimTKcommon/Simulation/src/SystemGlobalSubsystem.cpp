/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-15 Stanford University and the Authors.        *
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


/**@file
 *
 * Implementation of SystemGlobalSubsystem.
 */

#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"
#include "SimTKcommon/internal/Subsystem.h"
#include "SimTKcommon/internal/SubsystemGuts.h"
#include "SimTKcommon/internal/System.h"
#include "SimTKcommon/internal/Study.h"
#include "SimTKcommon/internal/SystemGlobalSubsystem.h"
#include "SimTKcommon/internal/Event.h"
#include "SimTKcommon/internal/EventAction.h"
#include "SimTKcommon/internal/EventTrigger.h"
#include "SimTKcommon/internal/EventTrigger_Timer.h"
#include "SimTKcommon/internal/EventTrigger_Witness.h"

#include "SimTKcommon/internal/EventHandler.h"
#include "SimTKcommon/internal/EventReporter.h"

#include <cassert>
#include <map>
#include <set>
#include <algorithm>
#include <utility>
#include <typeinfo>

using namespace SimTK;

//==============================================================================
//                       EVENT TRIGGER COLLECTION
//==============================================================================
// This local class holds a set of EventTrigger objects. It is intended for
// two purposes: once as a member of the System, for triggers that are always
// present, and once as the value type of a discrete state variable, for
// triggers that come and go at run time.
//
// The discrete variable is allocated at the start of realizeTopology() and
// contains no run time triggers at that point. After that triggers can be
// added and removed.
//
// Accessing the discrete variable for update doesn't invalidate any stage.
// However, adding or removing a trigger invalidates the results cache entry
// that holds the value for that trigger.
namespace {

class EventTriggerCollection {
public:
    unsigned adoptTimer(EventTrigger::Timer* timer) {
        return adoptThing(timer, m_timers, m_freeTimers);
    }
    void removeTimer(unsigned timerIndex) {
        removeThing(timerIndex, m_timers, m_freeTimers);
    }
    unsigned adoptWitness(EventTrigger::Witness* witness) {
        return adoptThing(witness, m_witnesses, m_freeWitnesses);
    }
    void removeWitness(unsigned witnessIndex) {
        removeThing(witnessIndex, m_witnesses, m_freeWitnesses);
    }
private:
    template <class T>
    static unsigned adoptThing(T* thing, 
                               Array_<std::unique_ptr<T>>& things,
                               Array_<unsigned>&           freeSlots) {
        const unsigned thingIndex = findFreeSlot(things,freeSlots);
        assert(!things[thingIndex]); // slot must be empty!
        things[thingIndex].reset(thing);
        return thingIndex;
    }

    template <class T>
    static void removeThing(unsigned thingIndex, 
                            Array_<std::unique_ptr<T>>& things,
                            Array_<unsigned>&           freeSlots) {
        assert(things[thingIndex]);
        if (thingIndex == things.size()-1)
            things.pop_back();
        else {
            things[thingIndex].reset(nullptr);
            freeSlots.push_back(thingIndex);
        }
    }

    template <class T>
    static unsigned findFreeSlot(Array_<std::unique_ptr<T>>& things,
                                 Array_<unsigned>&           freeSlots) {
        if (!freeSlots.empty()) {
            const unsigned slot = freeSlots.back();
            freeSlots.pop_back();
            return slot;
        }
        things.emplace_back(nullptr); // make room
        return things.size()-1;
    }
private:
    Array_<std::unique_ptr<EventTrigger::Timer>>    m_timers;
    Array_<unsigned>                                m_freeTimers;
    Array_<std::unique_ptr<EventTrigger::Witness>>  m_witnesses;
    Array_<unsigned>                                m_freeWitnesses;
};

}   // anonymous file-scope namespace


//==============================================================================
//                     SYSTEM GLOBAL SUBSYSTEM :: GUTS
//==============================================================================
// This is the implementation class for SystemGlobalSubsystem.

class SystemGlobalSubsystem::Guts : public Subsystem::Guts {
public:

    Guts() : Subsystem::Guts("SystemGlobalSubsystem", "0.0.1") {

    }
    
    ~Guts() = default;
    
    Guts* cloneImpl() const override {
        return new Guts(*this);
    }
        
    int realizeSubsystemTopologyImpl(State& s) const override;
    
    int realizeSubsystemModelImpl(State& s) const override {
        return 0;
    }

    int realizeSubsystemInstanceImpl(const State& s) const override {
        return 0;        
    }
    int realizeSubsystemTimeImpl(const State& s) const override {
        return realizeEventWitnesses(s, Stage::Time);
    }
    int realizeSubsystemPositionImpl(const State& s) const override {
        return realizeEventWitnesses(s, Stage::Position);
    }
    int realizeSubsystemVelocityImpl(const State& s) const override {
        return realizeEventWitnesses(s, Stage::Velocity);
    }
    int realizeSubsystemDynamicsImpl(const State& s) const override {
        return realizeEventWitnesses(s, Stage::Dynamics);
    }
    int realizeSubsystemAccelerationImpl(const State& s) const override {
        return realizeEventWitnesses(s, Stage::Acceleration);
    }
    int realizeSubsystemReportImpl(const State& s) const override {
        return realizeEventWitnesses(s, Stage::Report);
    }

    const Event::Initialization& getInitializationEvent() const {
        return Event::Initialization::downcast
                                        (*m_events[m_initializationEventId]);
    }

    const Event::TimeAdvanced& getTimeAdvancedEvent() const {
        return Event::TimeAdvanced::downcast
                                        (*m_events[m_timeAdvancedEventId]);
    }

    const Event::Termination& getTerminationEvent() const {
        return Event::Termination::downcast
                                        (*m_events[m_terminationEventId]);
    }

    const InitializationTrigger& getInitializationTrigger() const {
        return InitializationTrigger::downcast
                                        (*m_triggers[m_initializationTriggerId]);
    }

    const TimeAdvancedTrigger& getTimeAdvancedTrigger() const {
        return TimeAdvancedTrigger::downcast
                                        (*m_triggers[m_timeAdvancedTriggerId]);
    }

    const TerminationTrigger& getTerminationTrigger() const {
        return TerminationTrigger::downcast
                                        (*m_triggers[m_terminationTriggerId]);
    }

    static const int MaxDeriv = EventTrigger::Witness::MaxDeriv;

private:
    // Evaluate all witness values and derivatives for which g is the depends-on
    // stage. The result goes into the corresponding cache entry.
    int realizeEventWitnesses(const State& s, Stage g) const;


private:
friend class SystemGlobalSubsystem;

    //  TOPOLOGY STATE VARIABLES
    Array_<ClonePtr<Event>,EventId>                  m_events;
    Array_<ClonePtr<EventTrigger>,EventTriggerId>    m_triggers;

    // These store their assigned EventIds and EventTriggerIds.
    Array_<ClonePtr<ScheduledEventHandler>>  m_scheduledEventHandlers;
    Array_<ClonePtr<TriggeredEventHandler>>  m_triggeredEventHandlers;
    Array_<ClonePtr<ScheduledEventReporter>> m_scheduledEventReporters;
    Array_<ClonePtr<TriggeredEventReporter>> m_triggeredEventReporters;

    EventId         m_initializationEventId;    // for predefined Events
    EventId         m_timeAdvancedEventId;
    EventId         m_terminationEventId;
    EventId         m_extremeValueIsolatedEventId;

    EventTriggerId  m_initializationTriggerId;  // for predefined Triggers
    EventTriggerId  m_timeAdvancedTriggerId;
    EventTriggerId  m_terminationTriggerId;

    // TOPOLOGY CACHE VARIABLES

    // These reference objects that are in m_triggers.

    // Timers are always evaluated at the beginning of a step when the state
    // has been realized to Stage::Acceleration.
    Array_<EventTrigger::Timer*,EventTimerIndex>       m_timers;

    // Witness values and derivatives are partitioned by depends-on stage. This
    // array allows us to access values and derivatives by witness.
    Array_<EventTrigger::Witness*,EventWitnessIndex>   m_witnesses;

    // These arrays provides access to witness values and derivatives by stage.
    static const int NWitnessStages = Stage::NValid;
    Array_<EventWitnessIndex,int>  
        m_witnessesByStage[NWitnessStages][1+MaxDeriv];

    // Each of these (var,cache) pairs holds an Array of real values of
    // witness functions or their derivatives (mixed together).
    using   ValArray = Array_<Real,int>;
    DiscreteVariableIndex           m_witnessVars[NWitnessStages];
    CacheEntryIndex                 m_witnessValues[NWitnessStages];

    void clearCache() {
        m_timers.clear(); m_witnesses.clear();
        for (int g=0; g <NWitnessStages; ++g) {
            for (int d=0; d <= MaxDeriv; ++d)
                m_witnessesByStage[g][d].clear();
            m_witnessVars[g].invalidate();
            m_witnessValues[g].invalidate();
        }
    }
};

//------------------------------------------------------------------------------
//                     REALIZE SUBSYSTEM TOPOLOGY IMPL
//------------------------------------------------------------------------------
int SystemGlobalSubsystem::Guts::
realizeSubsystemTopologyImpl(State& s) const {
    auto mThis = const_cast<SystemGlobalSubsystem::Guts*>(this);
    mThis->clearCache();

    // Find all the Timers and Witnesses and dynamic_cast them once here.
    // Partition witness values and derivatives by depends-on stage of their 
    // functions.
    for (auto& trigger : mThis->m_triggers) {
        auto p = trigger.upd();

        // Deal with Timers.
        if (auto tp = dynamic_cast<EventTrigger::Timer*>(p)) {
            const EventTimerIndex timerIndex(m_timers.size());
            tp->m_timerIndex = timerIndex;
            mThis->m_timers.emplace_back(tp);
            continue;
        }
        
        // Deal with Witnesses.
        if (auto wp = dynamic_cast<EventTrigger::Witness*>(p)) {
            const EventWitnessIndex witnessIndex(m_witnesses.size());
            wp->m_witnessIndex = witnessIndex;
            mThis->m_witnesses.emplace_back(wp);
            // We'll calculate only up to MaxDeriv derivatives.
            const int nDerivs = std::min(MaxDeriv, wp->getNumTimeDerivatives());
            for (int deriv=0; deriv <= nDerivs; ++deriv) {
                const Stage g = wp->getDependsOnStage(deriv);
                mThis->m_witnessesByStage[g][deriv].push_back(witnessIndex);
            }
            continue;
        }

        // Nothing to do for other types of Triggers.
    }

    // Allocate auto-update state variables and corresponding cache entries
    // for witness values (including derivatives). The number of witnesses may
    // change at run time so we'll allocate variables even if there are 
    // currently no witness functions at some stages.
    for (Stage g(Stage::Topology); g <= Stage::Acceleration; ++g) {
        int nAtThisStage = 0;
        for (int deriv=0; deriv <= MaxDeriv; ++deriv) {
            for (auto wx : m_witnessesByStage[g][deriv]) {
                auto wp = m_witnesses[wx];
                wp->m_valueIndex[deriv] = nAtThisStage++;
            }
        }
        mThis->m_witnessVars[g] = 
            allocateAutoUpdateDiscreteVariable(s, Stage::Report,
                new Value<ValArray>(ValArray(nAtThisStage,NaN)), g);
        mThis->m_witnessValues[g] =
            getDiscreteVarUpdateIndex(s, m_witnessVars[g]);
    }

    // Allocate state variables and cache entries for the witness results.
    return 0;
}

//------------------------------------------------------------------------------
//                        REALIZE EVENT WITNESSES
//------------------------------------------------------------------------------
// For a given Stage, realize any witness value or derivative functions that
// have that stage as their depends-on stage. The results go into the cache
// entry for that stage, which is marked valid here.
int SystemGlobalSubsystem::Guts::
realizeEventWitnesses(const State& state, Stage g) const {
    //assert(Stage::Topology <= g && g <= Stage::Acceleration);
    //const System& system = this->getSystem();

    //ValArray& cache = Value<ValArray>::updDowncast
    //                        (updCacheEntry(state, m_witnessValues[g])).upd();
    //for (int deriv=0; deriv <= MaxDeriv; ++deriv) {
    //    const auto& witnesses = m_witnessesByStage[g][deriv];

    //    for (auto wx : witnesses) {
    //        const EventTrigger::Witness* const wp = m_witnesses[wx];
    //        const Real value = wp->calcWitnessValue(system, state, deriv);
    //        const int ix = wp->getValueIndex(deriv);
    //        cache[ix] = value;
    //    }
    //}
    //markCacheValueRealized(state, m_witnessValues[g]);

    return 0;
}

//==============================================================================
//               LOCAL CLASSES FOR EVENTHANDLER/REPORTER SUPPORT
//==============================================================================
// The EventHandler/EventReport facility preceded the current Event
// implementation. These classes are reimplemented here in terms of the new
// facility, with the aid of these file-local classes.
namespace {
//---------------------------- EVENT HANDLER EVENT -----------------------------
/* This is the System-wide Event that is generated by an EventHandler. */
class EventHandler_Event : public Event {
public:
    explicit EventHandler_Event(const EventHandler& eh) 
    :   Event(eh.getEventDescription().empty() 
                ? std::string("EventHandler Event") 
                : eh.getEventDescription()) {}
private:
    EventHandler_Event* cloneVirtual() const override 
    {   return new EventHandler_Event(*this); }
};

//--------------------------- EVENT REPORTER EVENT -----------------------------
/* This is the System-wide Event that is generated by an EventReporter. */
class EventReporter_Event : public Event {
public:
    explicit EventReporter_Event(const EventReporter& er) 
    :   Event(er.getEventDescription().empty() 
                ? std::string("EventReporter Event") 
                : er.getEventDescription()) {}
private:
    EventReporter_Event* cloneVirtual() const override 
    {   return new EventReporter_Event(*this); }
};


//---------------------------- EVENT HANDLER ACTION ----------------------------
/* This is the EventAction to be taken when an EventHandler-defined event
occurs. We just need to call the EventHandler's handleEvent() method. */
class EventHandler_Action : public EventAction {
public:
    EventHandler_Action(const EventHandler& handler) 
    :   EventAction(Change), m_handler(handler) {}

private:
    EventHandler_Action* cloneVirtual() const override 
    {   return new EventHandler_Action(*this); }

    void changeVirtual
       (Study&                  study,
        const Event&            event,
        const EventTriggers&    triggers,
        EventChangeResult&      result) const override
    {
        bool shouldTerminate = false;
        m_handler.handleEvent(study.updInternalState(), 
                              study.getAccuracyInUse(), shouldTerminate);
        result.reportExitStatus(shouldTerminate 
                              ? EventChangeResult::ShouldTerminate
                              : EventChangeResult::Succeeded);
    }

    const EventHandler& m_handler;
};

//--------------------------- EVENT REPORTER ACTION ----------------------------
/* This is the EventAction to be taken when an EventReporter-defined event
occurs. We just need to call the EventReporter's handleEvent() method. */
class EventReporter_Action : public EventAction {
public:
    EventReporter_Action(const EventReporter& reporter) 
    :   EventAction(Report), m_reporter(reporter) {}

private:
    EventReporter_Action* cloneVirtual() const override 
    {   return new EventReporter_Action(*this); }

    void reportVirtual
       (const Study&            study,
        const Event&            event,
        const EventTriggers&    triggers) const override
    {
        m_reporter.handleEvent(study.getCurrentState());
    }

    const EventReporter& m_reporter;
};

//--------------------- SCHEDULED EVENT HANDLER TIMER --------------------------
/* This is the EventTrigger::Timer generated by a ScheduledEventHandler. */
class ScheduledEventHandler_Timer : public EventTrigger::Timer {
    using Super = EventTrigger::Timer;
public:
    ScheduledEventHandler_Timer(const ScheduledEventHandler& handler) 
    :   Super("ScheduledEventHandler timer"), m_handler(handler) {}

private:
    ScheduledEventHandler_Timer* cloneVirtual() const override 
    {   return new ScheduledEventHandler_Timer(*this); }

    double calcTimeOfNextTriggerVirtual
       (const System& system, const State& state, 
        double timeOfLastTrigger) const override 
    {
        return (double)m_handler.getNextEventTime
                                (state, state.getTime() > timeOfLastTrigger);
    }


    const ScheduledEventHandler& m_handler;
};

//--------------------- SCHEDULED EVENT REPORTER TIMER -------------------------
/* This is the EventTrigger::Timer generated by a ScheduledEventReporter. */
class ScheduledEventReporter_Timer : public EventTrigger::Timer {
    using Super = EventTrigger::Timer;
public:
    ScheduledEventReporter_Timer(const ScheduledEventReporter& reporter) 
    :   Super("ScheduledEventReporter timer"), m_reporter(reporter) {}

private:
    ScheduledEventReporter_Timer* cloneVirtual() const override 
    {   return new ScheduledEventReporter_Timer(*this); }

    double calcTimeOfNextTriggerVirtual
       (const System& system, const State& state, 
        double timeOfLastTrigger) const override 
    {
        return (double)m_reporter.getNextEventTime
            (state, state.getTime() > timeOfLastTrigger);
    }

    const ScheduledEventReporter& m_reporter;
};

//--------------------- TRIGGERED EVENT HANDLER WITNESS ------------------------
/* This is the EventTrigger::Witness generated by a TriggeredEventHandler. */
class TriggeredEventHandler_Witness : public EventTrigger::Witness {
    using Super = EventTrigger::Witness;
public:
    TriggeredEventHandler_Witness(const TriggeredEventHandler& handler) 
    :   Super("TriggeredEventHandler witness"), m_handler(handler) {}

private:
    TriggeredEventHandler_Witness* cloneVirtual() const override 
    {   return new TriggeredEventHandler_Witness(*this); }

    Real calcWitnessValueVirtual(const System& /*system*/,
                                 const State&  state, 
                                 int           /*derivOrder*/) const override 
    {   return m_handler.getValue(state); }

    Stage getDependsOnStageVirtual(int /*derivOrder*/) const override
    {   return m_handler.getRequiredStage(); }

    int getNumTimeDerivativesVirtual() const override 
    {   return 0; }

private:
    const TriggeredEventHandler& m_handler;
};

//--------------------- TRIGGERED EVENT REPORTER WITNESS -----------------------
/* This is the EventTrigger::Witness generated by a TriggeredEventReporter. */
class TriggeredEventReporter_Witness : public EventTrigger::Witness {
    using Super = EventTrigger::Witness;
public:
    TriggeredEventReporter_Witness(const TriggeredEventReporter& reporter) 
    :   Super("TriggeredEventReporter witness"), m_reporter(reporter) {}

private:
    TriggeredEventReporter_Witness* cloneVirtual() const override 
    {   return new TriggeredEventReporter_Witness(*this); }

    Real calcWitnessValueVirtual(const System& /*system*/,
                                 const State&  state, 
                                 int           /*derivOrder*/) const override 
    {   return m_reporter.getValue(state); }

    Stage getDependsOnStageVirtual(int /*derivOrder*/) const override
    {   return m_reporter.getRequiredStage(); }

    int getNumTimeDerivativesVirtual() const override 
    {   return 0; }

private:
    const TriggeredEventReporter& m_reporter;
};

} // anonymous file-scope namespace


//==============================================================================
//                       SYSTEM GLOBAL SUBSYSTEM
//==============================================================================

SystemGlobalSubsystem::SystemGlobalSubsystem(System& sys) {
    adoptSubsystemGuts(new SystemGlobalSubsystem::Guts());

    SystemGlobalSubsystem::Guts& guts = updGuts();
    guts.m_initializationEventId = adoptEvent(new Event::Initialization());
    guts.m_timeAdvancedEventId   = adoptEvent(new Event::TimeAdvanced());
    guts.m_terminationEventId    = adoptEvent(new Event::Termination());
    guts.m_extremeValueIsolatedEventId = 
                                adoptEvent(new Event::ExtremeValueIsolated());

    guts.m_initializationTriggerId = 
        adoptEventTrigger(new InitializationTrigger(guts.m_initializationEventId));
    guts.m_timeAdvancedTriggerId = 
        adoptEventTrigger(new TimeAdvancedTrigger(guts.m_timeAdvancedEventId));
    guts.m_terminationTriggerId = 
        adoptEventTrigger(new TerminationTrigger(guts.m_terminationEventId));

    sys.adoptSubsystem(*this);
}

const SystemGlobalSubsystem::Guts& SystemGlobalSubsystem::getGuts() const {
    return SimTK_DYNAMIC_CAST_DEBUG<const SystemGlobalSubsystem::Guts&>
                                                        (getSubsystemGuts());
}

SystemGlobalSubsystem::Guts& SystemGlobalSubsystem::updGuts() {
    return SimTK_DYNAMIC_CAST_DEBUG<SystemGlobalSubsystem::Guts&>
                                                        (updSubsystemGuts());
}

//------------------------------------------------------------------------------
//                       ADOPT EVENT HANDLER (SCHEDULED)
//------------------------------------------------------------------------------
void SystemGlobalSubsystem::
adoptEventHandler(ScheduledEventHandler* handler) {
    SystemGlobalSubsystem::Guts& guts = updGuts();
    auto evnt = new EventHandler_Event(*handler);
    auto action = new EventHandler_Action(*handler);
    const EventActionIndex eax = evnt->adoptEventAction(action);
    const EventId eid = adoptEvent(evnt);
    auto timer = new ScheduledEventHandler_Timer(*handler);
    timer->addEvent(eid);
    const EventTriggerId tid = adoptEventTrigger(timer);
    handler->m_system    = &updSystem(); // This class is a friend.
    handler->m_eventId   = eid;
    handler->m_triggerId = tid;
    guts.m_scheduledEventHandlers.emplace_back(handler);
}

//------------------------------------------------------------------------------
//                       ADOPT EVENT HANDLER (TRIGGERED)
//------------------------------------------------------------------------------
void SystemGlobalSubsystem::
adoptEventHandler(TriggeredEventHandler* handler) {
    SystemGlobalSubsystem::Guts& guts = updGuts();
    auto evnt = new EventHandler_Event(*handler);
    auto action = new EventHandler_Action(*handler);
    const EventActionIndex eax = evnt->adoptEventAction(action);
    const EventId eid = adoptEvent(evnt);
    auto witness = new TriggeredEventHandler_Witness(*handler);
    witness->addEvent(eid);

    // Apply trigger info from the TriggeredEventHandler interface to the
    // witness we're creating here.
    const EventTriggerInfo& etInfo = handler->getTriggerInfo();
    witness->setTriggerOnRisingSignTransition
       (etInfo.shouldTriggerOnRisingSignTransition());
    witness->setTriggerOnFallingSignTransition
       (etInfo.shouldTriggerOnFallingSignTransition());
    witness->setAccuracyRelativeTimeLocalizationWindow
       (etInfo.getRequiredLocalizationTimeWindow());

    const EventTriggerId tid = adoptEventTrigger(witness);
    handler->m_system    = &updSystem(); // This class is a friend.
    handler->m_eventId   = eid;
    handler->m_triggerId = tid;
    handler->updTriggerInfo().setEventId(eid); //TODO: get rid of this
    guts.m_triggeredEventHandlers.emplace_back(handler);
}


//------------------------------------------------------------------------------
//                       ADOPT EVENT REPORTER (SCHEDULED)
//------------------------------------------------------------------------------
void SystemGlobalSubsystem::
adoptEventReporter(ScheduledEventReporter* reporter) {
    SystemGlobalSubsystem::Guts& guts = updGuts();
    auto evnt = new EventReporter_Event(*reporter);
    auto action = new EventReporter_Action(*reporter);
    const EventActionIndex eax = evnt->adoptEventAction(action);
    const EventId eid = adoptEvent(evnt);
    auto timer = new ScheduledEventReporter_Timer(*reporter);
    timer->addEvent(eid);
    const EventTriggerId tid = adoptEventTrigger(timer);
    reporter->m_system    = &updSystem(); // This class is a friend.
    reporter->m_eventId   = eid;
    reporter->m_triggerId = tid;
    guts.m_scheduledEventReporters.emplace_back(reporter);
}

//------------------------------------------------------------------------------
//                       ADOPT EVENT REPORTER (TRIGGERED)
//------------------------------------------------------------------------------
void SystemGlobalSubsystem::
adoptEventReporter(TriggeredEventReporter* reporter) {
    SystemGlobalSubsystem::Guts& guts = updGuts();
    auto evnt = new EventReporter_Event(*reporter);
    auto action = new EventReporter_Action(*reporter);
    const EventActionIndex eax = evnt->adoptEventAction(action);
    const EventId eid = adoptEvent(evnt);
    auto witness = new TriggeredEventReporter_Witness(*reporter);
    witness->addEvent(eid);

    // Apply trigger info from the TriggeredEventReporter interface to the
    // witness we're creating here.
    const EventTriggerInfo& etInfo = reporter->getTriggerInfo();
    witness->setTriggerOnRisingSignTransition
       (etInfo.shouldTriggerOnRisingSignTransition());
    witness->setTriggerOnFallingSignTransition
       (etInfo.shouldTriggerOnFallingSignTransition());
    witness->setAccuracyRelativeTimeLocalizationWindow
       (etInfo.getRequiredLocalizationTimeWindow());

    const EventTriggerId tid = adoptEventTrigger(witness);
    reporter->m_system    = &updSystem(); // This class is a friend.
    reporter->m_eventId   = eid;
    reporter->m_triggerId = tid;
    reporter->updTriggerInfo().setEventId(eid); //TODO: get rid of this
    guts.m_triggeredEventReporters.emplace_back(reporter);
}

//------------------------------------------------------------------------------
//                                 ADOPT EVENT
//------------------------------------------------------------------------------
EventId SystemGlobalSubsystem::
adoptEvent(Event* eventp) {
    SimTK_APIARGCHECK_ALWAYS(eventp != nullptr, "SystemGlobalSubsystem",
                             "adoptEvent", "Event pointer can't be null.");
    eventp->m_eventId = EventId(getGuts().m_events.size());
    updGuts().m_events.emplace_back(eventp);
    return eventp->m_eventId;
}


//------------------------------------------------------------------------------
//                             ADOPT EVENT TRIGGER
//------------------------------------------------------------------------------
EventTriggerId SystemGlobalSubsystem::
adoptEventTrigger(EventTrigger* triggerp) {
    SimTK_APIARGCHECK_ALWAYS(triggerp != nullptr, "SystemGlobalSubsystem",
        "adoptEventTrigger", "EventTrigger pointer can't be null.");
    triggerp->m_triggerId = EventTriggerId(getGuts().m_triggers.size());
    updGuts().m_triggers.emplace_back(triggerp);
    return triggerp->m_triggerId;
}


int SystemGlobalSubsystem::getNumEvents() const 
{   return (int)getGuts().m_events.size(); }

const Event& SystemGlobalSubsystem::getEvent(EventId id) const {
    SimTK_APIARGCHECK(id.isValid(), "SystemGlobalSubsystem",
                      "getEvent", "Uninitialized (invalid) EventId.");
    const auto& guts = getGuts();
    SimTK_INDEXCHECK_ALWAYS((int)id, getNumEvents(),
                            "SystemGlobalSubsystem::getEvent()");
    const Event* eventp = guts.m_events[id].get();
    SimTK_ERRCHK1_ALWAYS(eventp != nullptr, "SystemGlobalSubsystem::getEvent()",
                  "There is no Event associated with EventId(%d).", (int)id);
    return *eventp;
}

Event& SystemGlobalSubsystem::updEvent(EventId id) {
    SimTK_APIARGCHECK(id.isValid(), "SystemGlobalSubsystem",
                      "updEvent", "Uninitialized (invalid) EventId.");
    auto& guts = updGuts();
    SimTK_INDEXCHECK_ALWAYS((int)id, getNumEvents(),
                            "SystemGlobalSubsystem::updEvent()");
    Event* eventp = guts.m_events[id].upd();
    SimTK_ERRCHK1_ALWAYS(eventp != nullptr, "SystemGlobalSubsystem::updEvent()",
                  "There is no Event associated with EventId(%d).", (int)id);
    return *eventp;
}

bool SystemGlobalSubsystem::hasEvent(EventId id) const {
    SimTK_APIARGCHECK(id.isValid(), "SystemGlobalSubsystem",
                      "hasEvent", "Uninitialized (invalid) EventId.");
    return    SimTK::isIndexInRange((int)id, getNumEvents()) 
           && getGuts().m_events[id].get() != nullptr;
}

int SystemGlobalSubsystem::getNumEventTriggers() const 
{   return (int)getGuts().m_triggers.size(); }

const EventTrigger& SystemGlobalSubsystem::
getEventTrigger(EventTriggerId id) const {
    SimTK_APIARGCHECK(id.isValid(), "SystemGlobalSubsystem", "getEventTrigger", 
                      "Uninitialized (invalid) EventTriggerId.");
    const auto& guts = getGuts();
    SimTK_INDEXCHECK_ALWAYS((int)id, getNumEventTriggers(),
                            "SystemGlobalSubsystem::getEventTrigger()");
    const EventTrigger* eventp = guts.m_triggers[id].get();
    SimTK_ERRCHK1_ALWAYS(eventp != nullptr, 
        "SystemGlobalSubsystem::getEventTrigger()",
        "There is no EventTrigger associated with EventTriggerId(%d).", 
        (int)id);
    return *eventp;
}

EventTrigger& SystemGlobalSubsystem::
updEventTrigger(EventTriggerId id) {
    SimTK_APIARGCHECK(id.isValid(), "SystemGlobalSubsystem",
                      "updEventTrigger", "Uninitialized (invalid) EventId.");
    auto& guts = updGuts();
    SimTK_INDEXCHECK_ALWAYS((int)id, getNumEventTriggers(),
                            "SystemGlobalSubsystem::updEventTrigger()");
    EventTrigger* eventp = guts.m_triggers[id].upd();
    SimTK_ERRCHK1_ALWAYS(eventp != nullptr, 
        "SystemGlobalSubsystem::updEventTrigger()",
        "There is no EventTrigger associated with EventTriggerId(%d).", 
        (int)id);
    return *eventp;
}

bool SystemGlobalSubsystem::
hasEventTrigger(EventTriggerId id) const {
    SimTK_APIARGCHECK(id.isValid(), "SystemGlobalSubsystem",
        "hasEventTrigger", "Uninitialized (invalid) EventTriggerId.");
    return    SimTK::isIndexInRange((int)id, getNumEventTriggers()) 
           && getGuts().m_triggers[id].get() != nullptr;
}

EventId SystemGlobalSubsystem::getInitializationEventId() const
{   return getGuts().m_initializationEventId; }
EventId SystemGlobalSubsystem::getTimeAdvancedEventId() const
{   return getGuts().m_timeAdvancedEventId; }
EventId SystemGlobalSubsystem::getTerminationEventId() const
{   return getGuts().m_terminationEventId; }
EventId SystemGlobalSubsystem::getExtremeValueIsolatedEventId() const
{   return getGuts().m_extremeValueIsolatedEventId; }

EventTriggerId SystemGlobalSubsystem::getInitializationTriggerId() const
{   return getGuts().m_initializationTriggerId; }
EventTriggerId SystemGlobalSubsystem::getTimeAdvancedTriggerId() const
{   return getGuts().m_timeAdvancedTriggerId; }
EventTriggerId SystemGlobalSubsystem::getTerminationTriggerId() const
{   return getGuts().m_terminationTriggerId; }

void SystemGlobalSubsystem::
findActiveEventWitnesses(const Study&                          study, 
                         Array_<const EventTrigger::Witness*,
                                ActiveWitnessIndex>&           witnesses) const
{
    auto& guts = getGuts();

    //TODO: collect dynamic witnesses from state.
    witnesses.assign(guts.m_witnesses.begin(), guts.m_witnesses.end());
}


void SystemGlobalSubsystem::
findActiveEventTimers(const Study&                          study, 
                      Array_<const EventTrigger::Timer*,
                             ActiveTimerIndex>&             timers) const {
    auto& guts = getGuts();

    //TODO: collect dynamic timers from state.
    timers.assign(guts.m_timers.begin(), guts.m_timers.end());
}

void SystemGlobalSubsystem::
findNextScheduledEventTimes
   (const Study&        study,
    double              timeOfLastReport,
    double              timeOfLastChange,
    double&             timeOfNextReport,
    EventTriggers&      reportTimers,
    double&             timeOfNextChange,
    EventTriggers&      changeTimers) const 
{
    auto& system = study.getSystem();
    auto& guts = getGuts();

    timeOfNextReport = timeOfNextChange = dInfinity;
    reportTimers.clear(); changeTimers.clear();

    for (auto timerp : guts.m_timers) {
        bool hasChangeAction = false;
        for (int i = 0; i < timerp->getNumEvents(); ++i) {
            const EventId eid = timerp->getEventId(i);
            const Event& e = system.getEvent(eid);
            if ((hasChangeAction=e.hasChangeAction()))
                break;
        }
        if (hasChangeAction) {
            const double t = timerp->calcTimeOfNextTrigger
                            (system, study.getCurrentState(), timeOfLastChange);
            if (t > timeOfNextChange)
                continue; // This one is not interesting.

            if (t < timeOfNextChange) {
                changeTimers.clear(); // forget previous earliest
                timeOfNextChange = t;
            }
            // Add to list if new winner or same as previous winner.
            changeTimers.push_back(timerp);
        } else { // timer triggers only report actions
            const double t = timerp->calcTimeOfNextTrigger
                            (system, study.getCurrentState(), timeOfLastReport);
            if (t > timeOfNextReport)
                continue; // This one is not interesting.

            if (t < timeOfNextReport) {
                reportTimers.clear(); // forget previous earliest
                timeOfNextReport = t;
            }

            reportTimers.push_back(timerp);
        }
    }

}


//------------------------------------------------------------------------------
//                           NOTE EVENT OCCURENCE
//------------------------------------------------------------------------------
/* We're given a list of event triggers that a time stepper declares have 
occurred simultaneously. Each of those contains a list of EventIds that are
caused by that trigger. We assume that the triggers are unique, but several
triggers may cause the same event. However, each caused event should occur only
once, and for each unique event we need to know which triggers caused it. We
map each EventId to its corresponding Event object, ignoring any EventIds that
are not recognized.

We're assuming these are *very* short lists (typically one trigger causing
a single event) so are using the least-overhead algorithms possible. A better 
algorithm would be required if many events could be triggered at once. But this
one is good even if many triggers cause the same event. 

Mutable occurrence counters are bumped here, once per trigger and once per 
unique event caused. */
void SystemGlobalSubsystem::
noteEventOccurrence(const EventTriggers&    triggers,
                    EventsAndCauses&        appendTriggeredEvents,
                    Array_<EventId>&        appendIgnoredEvents) const {
    // We expect there to be very few events (typically, 1) so this linear 
    // search should be fastest, despite its apparent O(n^2) complexity.
    for (auto trigger : triggers) {
        trigger->noteOccurrence(); // bump mutable counter
        for (int i = 0; i < trigger->getNumEvents(); ++i) {
            const EventId eid = trigger->getEventId(i);

            if (!hasEvent(eid)) { // ignore unrecognized EventId
                auto p = std::find(appendIgnoredEvents.cbegin(), //linear search
                                   appendIgnoredEvents.cend(), eid);
                if (p == appendIgnoredEvents.cend())
                    appendIgnoredEvents.push_back(eid);
                continue;
            }

            // Find or insert output entry e for this event.
            const Event& evnt = getEvent(eid);
            auto e = std::find_if // linear search
               (appendTriggeredEvents.begin(), appendTriggeredEvents.end(),
                [&evnt](const EventAndCausesPair& et){return et.first==&evnt;});
            if (e == appendTriggeredEvents.end()) {
                evnt.noteOccurrence(); // This is a new event; bump counter
                appendTriggeredEvents.emplace_back
                                            (&evnt, EventTriggers{trigger});
            } else {
                // Just add this trigger as a cause for event at e.
                e->second.push_back(trigger);
            }
        }
    }
}

//------------------------------------------------------------------------------
//                        PERFORM EVENT REPORT ACTIONS
//------------------------------------------------------------------------------
void SystemGlobalSubsystem::
performEventReportActions
   (const Study&            study,
    const EventsAndCauses&  triggeredEvents) const {
    assert(!triggeredEvents.empty());

    for (auto&& et : triggeredEvents) {
        et.first->performReportActions(study, et.second);
    }
}

//------------------------------------------------------------------------------
//                        PERFORM EVENT CHANGE ACTIONS
//------------------------------------------------------------------------------
void SystemGlobalSubsystem::
performEventChangeActions
   (Study&                  study,
    const EventsAndCauses&  triggeredEvents,
    EventChangeResult&      result) const {
    assert(!triggeredEvents.empty());

    State& state  = study.updInternalState();

    // Save the stage version numbers so we can look for changes.
    Array_<StageVersion> stageVersions;
    state.getSystemStageVersions(stageVersions);

    // Results are accumulated by the actions. Start empty.
    result.clear();

    for (auto&& et : triggeredEvents) {
        et.first->performChangeActions(study, et.second, result);
    }

    // Note the lowest stage whose version was changed by the actions.
    const Stage lowestModified = 
        state.getLowestSystemStageDifference(stageVersions);
    result.setLowestModifiedStage(lowestModified);
}

