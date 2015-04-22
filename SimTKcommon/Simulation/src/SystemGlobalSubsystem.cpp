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
#include "SimTKcommon/internal/SystemGlobalSubsystem.h"
#include "SimTKcommon/internal/EventHandler.h"
#include "SimTKcommon/internal/EventReporter.h"

#include <cassert>
#include <map>
#include <set>

namespace SimTK {


//==============================================================================
//                         SYSTEM GLOBAL SUBSYSTEM
//==============================================================================

class SystemGlobalSubsystem::Guts : public Subsystem::Guts {
public:

    Guts() : Subsystem::Guts("SystemGlobalSubsystem::Guts", "0.0.1") { }
    
    ~Guts() {
        for (unsigned i=0; i < m_scheduledEventHandlers.size(); ++i)
            delete m_scheduledEventHandlers[i];
        for (unsigned i=0; i < m_triggeredEventHandlers.size(); ++i)
            delete m_triggeredEventHandlers[i];
        for (unsigned i=0; i < m_scheduledEventReporters.size(); ++i)
            delete m_scheduledEventReporters[i];
        for (unsigned i=0; i < m_triggeredEventReporters.size(); ++i)
            delete m_triggeredEventReporters[i];
    }
    
    Guts* cloneImpl() const override {
        return new Guts(*this);
    }
        
    const Array_<ScheduledEventHandler*>& getScheduledEventHandlers() const {
        return m_scheduledEventHandlers;
    }
    
    Array_<ScheduledEventHandler*>& updScheduledEventHandlers() {
        invalidateSubsystemTopologyCache();
        return m_scheduledEventHandlers;
    }
    
    const Array_<TriggeredEventHandler*>& getTriggeredEventHandlers() const {
        return m_triggeredEventHandlers;
    }
    
    Array_<TriggeredEventHandler*>& updTriggeredEventHandlers() {
        invalidateSubsystemTopologyCache();
        return m_triggeredEventHandlers;
    }
    
    const Array_<ScheduledEventReporter*>& getScheduledEventReporters() const {
        return m_scheduledEventReporters;
    }
    
    Array_<ScheduledEventReporter*>& updScheduledEventReporters() const {
        invalidateSubsystemTopologyCache();
        return m_scheduledEventReporters;
    }
    
    const Array_<TriggeredEventReporter*>& getTriggeredEventReporters() const {
        return m_triggeredEventReporters;
    }
    
    Array_<TriggeredEventReporter*>& updTriggeredEventReporters() const {
        invalidateSubsystemTopologyCache();
        return m_triggeredEventReporters;
    }

    int realizeSubsystemTopologyImpl(State& s) const override {
        auto mThis = const_cast<SystemGlobalSubsystem::Guts*>(this);
        mThis->clearCache();
        const System& sys = getSystem();
        const SubsystemIndex myIx = getMySubsystemIndex();

        // Allocate EventIds for scheduled events and reports.
        for (auto seh : m_scheduledEventHandlers)
            mThis->m_scheduledEventIds.push_back(sys.createNewEventId(myIx));
        for (auto ser : m_scheduledEventReporters)
            mThis->m_scheduledReportIds.push_back(sys.createNewEventId(myIx));

        // Allocate EventIds and witness function slots for triggered 
        // events & reports.
        for (auto teh : m_triggeredEventHandlers) {
            mThis->m_triggeredEventIds.push_back(sys.createNewEventId(myIx));
            const EventTriggerByStageIndex index = 
                s.allocateEventTrigger(myIx, teh->getRequiredStage(), 1);
            mThis->m_triggeredEventIndices.push_back(index);
        }
        for (auto ter : m_triggeredEventReporters) {
            mThis->m_triggeredReportIds.push_back(sys.createNewEventId(myIx));
            const EventTriggerByStageIndex index = 
                s.allocateEventTrigger(myIx, ter->getRequiredStage(), 1);
            mThis->m_triggeredReportIndices.push_back(index);
        }
        return 0;
    }
    
    int realizeSubsystemModelImpl(State& s) const override {
        return 0;
    }

    int realizeEventTriggers(const State& s, Stage g) const {
        Vector& triggers = s.updEventTriggersByStage(getMySubsystemIndex(), g);
        for (unsigned i=0; i < m_triggeredEventHandlers.size(); ++i) {
            if (g == m_triggeredEventHandlers[i]->getRequiredStage())
                triggers[m_triggeredEventIndices[i]] = 
                    m_triggeredEventHandlers[i]->getValue(s);
        }
        for (unsigned i=0; i < m_triggeredEventReporters.size(); ++i) {
            if (g == m_triggeredEventReporters[i]->getRequiredStage())
                triggers[m_triggeredReportIndices[i]] = 
                    m_triggeredEventReporters[i]->getValue(s);
        }
        return 0;
    }
    
    int realizeSubsystemInstanceImpl(const State& s) const override {
        return 0;        
    }
    int realizeSubsystemTimeImpl(const State& s) const override {
        return realizeEventTriggers(s, Stage::Time);
    }
    int realizeSubsystemPositionImpl(const State& s) const override {
        return realizeEventTriggers(s, Stage::Position);
    }
    int realizeSubsystemVelocityImpl(const State& s) const override {
        return realizeEventTriggers(s, Stage::Velocity);
    }
    int realizeSubsystemDynamicsImpl(const State& s) const override {
        return realizeEventTriggers(s, Stage::Dynamics);
    }
    int realizeSubsystemAccelerationImpl(const State& s) const override {
        return realizeEventTriggers(s, Stage::Acceleration);
    }
    int realizeSubsystemReportImpl(const State& s) const override {
        return realizeEventTriggers(s, Stage::Report);
    }

    void calcEventTriggerInfoImpl
       (const State& s, Array_<EventTriggerInfo>& trigInfo) const override 
    {       
        // Loop over all registered TriggeredEventHandlers and 
        // TriggeredEventReporters, and ask each one for its EventTriggerInfo.
        trigInfo.resize(  m_triggeredEventHandlers.size()
                        + m_triggeredEventReporters.size());

        for (unsigned i=0; i < m_triggeredEventHandlers.size(); ++i) {
            const Stage stage = m_triggeredEventHandlers[i]->getRequiredStage();
            const SystemEventTriggerIndex index
               (  m_triggeredEventIndices[i]
                + s.getEventTriggerStartByStage(stage)
                + s.getEventTriggerStartByStage(getMySubsystemIndex(), stage));
            trigInfo[index] = m_triggeredEventHandlers[i]->getTriggerInfo();
            trigInfo[index].setEventId(m_triggeredEventIds[i]);
        }

        for (unsigned i=0; i < m_triggeredEventReporters.size(); ++i) {
            const Stage stage = m_triggeredEventReporters[i]->getRequiredStage();
            const SystemEventTriggerIndex index
               (  m_triggeredReportIndices[i]
                + s.getEventTriggerStartByStage(stage)
                + s.getEventTriggerStartByStage(getMySubsystemIndex(), stage));
            trigInfo[index] = m_triggeredEventReporters[i]->getTriggerInfo();
            trigInfo[index].setEventId(m_triggeredReportIds[i]);
        }
    }
    void calcTimeOfNextScheduledEventImpl(const State& s, Real& tNextEvent, 
        Array_<EventId>& eventIds, bool includeCurrentTime) const override {      
        // Loop over all registered ScheduledEventHandlers, and ask each one 
        // when its next event occurs.
        tNextEvent = Infinity;
        for (unsigned i = 0; i < m_scheduledEventHandlers.size(); ++i) {
            Real time = m_scheduledEventHandlers[i]->getNextEventTime
                                                        (s, includeCurrentTime);
            if (time <= tNextEvent 
                && (time > s.getTime() 
                    || (includeCurrentTime && time == s.getTime()))) 
            {
                if (time < tNextEvent)
                    eventIds.clear();
                tNextEvent = time;
                eventIds.push_back(m_scheduledEventIds[i]);
            }
        }
    }
    void calcTimeOfNextScheduledReportImpl(const State& s, Real& tNextEvent, 
        Array_<EventId>& eventIds, bool includeCurrentTime) const override {      
        // Loop over all registered ScheduledEventReporters, and ask each one 
        // when its next event occurs.      
        tNextEvent = Infinity;
        for (unsigned i = 0; i < m_scheduledEventReporters.size(); ++i) {
            Real time = m_scheduledEventReporters[i]->getNextEventTime
                                                        (s, includeCurrentTime);
            if (time <= tNextEvent 
                && (time > s.getTime() 
                    || (includeCurrentTime && time == s.getTime()))) 
            {
                if (time < tNextEvent)
                    eventIds.clear();
                tNextEvent = time;
                eventIds.push_back(m_scheduledReportIds[i]);
            }
        }
    }
    void handleEventsImpl(State& s, Event::Cause cause, 
                      const Array_<EventId>& eventIds, 
                      const HandleEventsOptions& options, 
                      HandleEventsResults& results) const override 
    {
        const Real accuracy = options.getAccuracy();
        bool shouldTerminate = false;
        
        // Build a set of the ids for quick lookup.      
        std::set<EventId> idSet;
        for (unsigned i=0; i < eventIds.size(); ++i)
            idSet.insert(eventIds[i]);
        
        // Process triggered events and reports.     
        if (cause == Event::Cause::Triggered) {
            for (unsigned i = 0; i < m_triggeredEventHandlers.size(); ++i) {
                if (idSet.find(m_triggeredEventIds[i]) != idSet.end()) {
                    bool eventShouldTerminate = false;
                    m_triggeredEventHandlers[i]->handleEvent
                                            (s, accuracy, eventShouldTerminate);
                    if (eventShouldTerminate)
                        shouldTerminate = true;
                }
            }
            for (unsigned i=0; i < m_triggeredEventReporters.size(); ++i) {
                if (idSet.find(m_triggeredReportIds[i]) != idSet.end())
                    m_triggeredEventReporters[i]->handleEvent(s);
            }
        }
        
        // Process scheduled events and reports.       
        if (cause == Event::Cause::Scheduled) {
            for (unsigned i=0; i < m_scheduledEventHandlers.size(); ++i) {
                if (idSet.find(m_scheduledEventIds[i]) != idSet.end()) {
                    bool eventShouldTerminate = false;
                    m_scheduledEventHandlers[i]->handleEvent
                                            (s, accuracy, eventShouldTerminate);
                    if (eventShouldTerminate)
                        shouldTerminate = true;
                }
            }
            for (unsigned i=0; i < m_scheduledEventReporters.size(); ++i) {
                if (idSet.find(m_scheduledReportIds[i]) != idSet.end())
                    m_scheduledEventReporters[i]->handleEvent(s);
            }
        }

        // Assume some change was made.
        results.setAnyChangeMade(true);

        results.setExitStatus(shouldTerminate 
            ? HandleEventsResults::ShouldTerminate
            : HandleEventsResults::Succeeded);
    }

    void reportEventsImpl(const State& s, Event::Cause cause, 
                          const Array_<EventId>& eventIds) const override
    {
        // Build a set of the ids for quick lookup.        
        std::set<EventId> idSet;
        for (unsigned i=0; i < eventIds.size(); ++i)
            idSet.insert(eventIds[i]);
        
        // Process triggered reports.       
        if (cause == Event::Cause::Triggered) {
            for (unsigned i=0; i < m_triggeredEventReporters.size(); ++i) {
                if (idSet.find(m_triggeredReportIds[i]) != idSet.end())
                    m_triggeredEventReporters[i]->handleEvent(s);
            }
        }
        
        // Process scheduled reports.      
        if (cause == Event::Cause::Scheduled) {
            for (unsigned i=0; i < m_scheduledEventReporters.size(); ++i) {
                if (idSet.find(m_scheduledReportIds[i]) != idSet.end())
                    m_scheduledEventReporters[i]->handleEvent(s);
            }
        }
    }

private:
friend class SystemGlobalSubsystem;

    //  TOPOLOGY STATE VARIABLES
    Array_<ScheduledEventHandler*>          m_scheduledEventHandlers;
    Array_<TriggeredEventHandler*>          m_triggeredEventHandlers;

    // We allow these to be added later since they are harmless.
    mutable Array_<ScheduledEventReporter*> m_scheduledEventReporters;
    mutable Array_<TriggeredEventReporter*> m_triggeredEventReporters;

    // TOPOLOGY CACHE VARIABLES
    Array_<EventId>                         m_scheduledEventIds;
    Array_<EventId>                         m_scheduledReportIds;

    Array_<EventId>                         m_triggeredEventIds;
    Array_<EventTriggerByStageIndex>        m_triggeredEventIndices;

    Array_<EventId>                         m_triggeredReportIds;
    Array_<EventTriggerByStageIndex>        m_triggeredReportIndices;

    void clearCache() {
        m_scheduledEventIds.clear();
        m_scheduledReportIds.clear();       
        m_triggeredEventIds.clear();
        m_triggeredEventIndices.clear();       
        m_triggeredReportIds.clear();
        m_triggeredReportIndices.clear();
    }
};

SystemGlobalSubsystem::SystemGlobalSubsystem(System& sys) {
    adoptSubsystemGuts(new SystemGlobalSubsystem::Guts());
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

/*
 * Add a ScheduledEventHandler to the System.  This must be called before the 
 * Model stage is realized.
 * 
 * The System assumes ownership of the object passed to this method, and will 
 * delete it when the System is deleted.
 */
void SystemGlobalSubsystem::addEventHandler(ScheduledEventHandler* handler) {
    updGuts().updScheduledEventHandlers().push_back(handler);
}

/*
 * Add a TriggeredEventHandler to the System.  This must be called before the 
 * Model stage is realized.
 * 
 * The System assumes ownership of the object passed to this method, and will 
 * delete it when the System is deleted.
 */
void SystemGlobalSubsystem::
addEventHandler(TriggeredEventHandler* handler) {
    updGuts().updTriggeredEventHandlers().push_back(handler);
}

/*
 * Add a ScheduledEventReporter to the System.  This must be called before the 
 * Model stage is realized.
 * 
 * The System assumes ownership of the object passed to this method, and will 
 * delete it when the System is deleted.
 * 
 * Note that this method is const.  Because an EventReporter cannot affect the
 * behavior of the system being simulated, it is permitted to add one to a 
 * const System.
 */
void SystemGlobalSubsystem::
addEventReporter(ScheduledEventReporter* handler) const {
    getGuts().updScheduledEventReporters().push_back(handler);
}

/*
 * Add a TriggeredEventReporter to the System.  This must be called before the 
 * Model stage is realized.
 * 
 * The System assumes ownership of the object passed to this method, and will 
 * delete it when the System is deleted.
 * 
 * Note that this method is const.  Because an EventReporter cannot affect the 
 * behavior of the system being simulated, it is permitted to add one to a 
 * const System.
 */
void SystemGlobalSubsystem::
addEventReporter(TriggeredEventReporter* handler) const {
    getGuts().updTriggeredEventReporters().push_back(handler);
}

} // namespace SimTK

