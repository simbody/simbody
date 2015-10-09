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
 * Implementation of System and System::Guts.
 */

#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"
#include "SimTKcommon/internal/System.h"
#include "SimTKcommon/internal/SystemGlobalSubsystem.h"
#include "SimTKcommon/internal/SystemGuts.h"

#include "SimTKcommon/internal/EventReporter.h"

#include <cassert>
#include <map>
#include <set>

using namespace SimTK;



//==============================================================================
//                                   SYSTEM
//==============================================================================
// Most of the System:: methods just forward to System::Guts.

// Put new unowned Guts into this handle and take over ownership. Add the
// SystemGlobalSubsystem as the 0th subsystem entry. Allocate needed resources.
System::System(System::Guts* guts) : m_guts(nullptr) {
    SimTK_ASSERT_ALWAYS(guts, "System::System(): can't adopt null Guts.");
    SimTK_ASSERT_ALWAYS(guts->getNumSubsystems() == 0,
        "System::System(): can't adopt Guts that already has Subsystems.");
    m_guts = guts;
    m_guts->setOwnerHandle(*this);

    // Now add the Default Subsystem, which will be SubsystemIndex(0).
    SystemGlobalSubsystem sysglobal(*this); // Invokes adoptSubsystem().

}

System::~System() {
    if (m_guts && isOwnerHandle())
        delete m_guts;
    m_guts = nullptr; // just to be tidy
}

bool System::isOwnerHandle() const 
{   return isEmptyHandle() || &m_guts->getSystem()==this; }
bool System::isSameSystem(const System& otherSystem) const 
{   return m_guts && (m_guts == otherSystem.m_guts); }

const String& System::getName()    const {return getSystemGuts().getName();}
const String& System::getVersion() const {return getSystemGuts().getVersion();}

System& System::setDefaultTimeScale(double tc) 
{   updSystemGuts().setDefaultTimeScale(tc); return *this; }
double System::getDefaultTimeScale() const
{   return getSystemGuts().getDefaultTimeScale(); }
System& System::setDefaultLengthScale(Real lc)
{   updSystemGuts().setDefaultLengthScale(lc); return *this; }
Real System::getDefaultLengthScale() const
{   return getSystemGuts().getDefaultLengthScale(); }

System& System::setUpDirection(const CoordinateDirection& up) 
{   updSystemGuts().setUpDirection(up); return *this; }
CoordinateDirection System::getUpDirection() const
{   return getSystemGuts().getUpDirection(); }
System& System::setUseUniformBackground(bool useUniformBackground)
{   updSystemGuts().setUseUniformBackground(useUniformBackground);
    return *this; }
bool System::getUseUniformBackground() const
{   return getSystemGuts().getUseUniformBackground(); }

void System::resetAllCountersToZero() {updSystemGuts().resetAllCounters();}
int System::getNumRealizationsOfThisStage(Stage g) const 
{return getSystemGuts().nRealizationsOfStage[g];}
int System::getNumRealizeCalls() const {return getSystemGuts().nRealizeCalls;}

int System::getNumPrescribeQCalls() const 
{return getSystemGuts().nPrescribeQCalls;}
int System::getNumPrescribeUCalls() const 
{return getSystemGuts().nPrescribeUCalls;}

int System::getNumProjectQCalls() const {return getSystemGuts().nProjectQCalls;}
int System::getNumProjectUCalls() const {return getSystemGuts().nProjectUCalls;}
int System::getNumFailedProjectQCalls() const 
{   return getSystemGuts().nFailedProjectQCalls; }
int System::getNumFailedProjectUCalls() const 
{   return getSystemGuts().nFailedProjectUCalls; }
int System::getNumQProjections() const {return getSystemGuts().nQProjections;}
int System::getNumUProjections() const {return getSystemGuts().nUProjections;} 
int System::getNumQErrorEstimateProjections() const 
{   return getSystemGuts().nQErrEstProjections; }
int System::getNumUErrorEstimateProjections() const 
{   return getSystemGuts().nUErrEstProjections; }

const State& System::getDefaultState() const 
{   return getSystemGuts().getDefaultState(); }

SubsystemIndex System::adoptSubsystem(Subsystem& child) 
{   return updSystemGuts().adoptSubsystem(child); }

int System::getNumSubsystems() const 
{   return getSystemGuts().getNumSubsystems(); }

const Subsystem& System::getSubsystem(SubsystemIndex i) const 
{   return getSystemGuts().getSubsystem(i); }
Subsystem& System::updSubsystem(SubsystemIndex i) 
{   return updSystemGuts().updSubsystem(i); }

const SystemGlobalSubsystem& System::getSystemGlobalSubsystem() const {
    return static_cast<const SystemGlobalSubsystem&>
                            (getSystemGuts().getSubsystem(SubsystemIndex(0)));
}
SystemGlobalSubsystem& System::updSystemGlobalSubsystem() {
    return static_cast<SystemGlobalSubsystem&>
                            (updSystemGuts().updSubsystem(SubsystemIndex(0)));
}

// TODO: this should be a Model stage variable allocated by the base class.
// Currently it is just a Topology stage variable stored in the base class.
void System::setHasTimeAdvancedEvents(bool hasEm) {
    updSystemGuts().setHasTimeAdvancedEvents(hasEm);
    getSystemGuts().invalidateSystemTopologyCache();
}
bool System::hasTimeAdvancedEvents() const {
    return getSystemGuts().hasTimeAdvancedEvents();
}
bool System::systemTopologyHasBeenRealized() const 
{   return getSystemGuts().systemTopologyHasBeenRealized(); }
StageVersion System::getSystemTopologyCacheVersion() const
{   return getSystemGuts().getSystemTopologyCacheVersion(); }
void System::setSystemTopologyCacheVersion(StageVersion topoVersion) const
{   getSystemGuts().setSystemTopologyCacheVersion(topoVersion); }
void System::invalidateSystemTopologyCache() const
{   getSystemGuts().invalidateSystemTopologyCache(); }



const State& System::realizeTopology() const 
{   return getSystemGuts().realizeTopology(); }
void System::realizeModel(State& s) const 
{   getSystemGuts().realizeModel(s); }

void System::realize(const State& s, Stage g) const 
{   getSystemGuts().realize(s,g); }

void System::calcDecorativeGeometryAndAppend
   (const State& s, Stage g, Array_<DecorativeGeometry>& geom) const 
{   getSystemGuts().calcDecorativeGeometryAndAppend(s,g,geom); }

void System::
multiplyByN(const State& s, const Vector& u, Vector& dq) const
{   getSystemGuts().multiplyByN(s,u,dq); }
void System::
multiplyByNTranspose(const State& s, const Vector& fq, Vector& fu) const
{   getSystemGuts().multiplyByNTranspose(s,fq,fu); }
void System::
multiplyByNPInv(const State& s, const Vector& dq, Vector& u) const
{   getSystemGuts().multiplyByNPInv(s,dq,u); }
void System::
multiplyByNPInvTranspose(const State& s, const Vector& fu, Vector& fq) const
{   getSystemGuts().multiplyByNPInvTranspose(s,fu,fq); }

bool System::prescribeQ(State& s) const
{   return getSystemGuts().prescribeQ(s); }
bool System::prescribeU(State& s) const
{   return getSystemGuts().prescribeU(s); }
void System::getFreeQIndex(const State& s, Array_<SystemQIndex>& freeQs) const
{   return getSystemGuts().getFreeQIndex(s,freeQs); }
void System::getFreeUIndex(const State& s, Array_<SystemUIndex>& freeUs) const
{   return getSystemGuts().getFreeUIndex(s,freeUs); }

void System::project(State& state, Real accuracy) const {
    const ProjectOptions projOptions(accuracy);
    ProjectResults projResults;
    Vector noErrEst;    // empty

    realize(state, Stage::Time);
    prescribeQ(state);
    realize(state, Stage::Position);
    projectQ(state, noErrEst, projOptions, projResults);
    prescribeU(state);
    realize(state, Stage::Velocity);
    projectU(state, noErrEst, projOptions, projResults);
}

void System::projectQ(State& state, Real accuracy) const {
    const ProjectOptions projOptions(accuracy);
    ProjectResults projResults;
    Vector noErrEst;    // empty

    realize(state, Stage::Time);
    prescribeQ(state);
    realize(state, Stage::Position);
    projectQ(state, noErrEst, projOptions, projResults);
}

void System::projectU(State& state, Real accuracy) const {
    const ProjectOptions projOptions(accuracy);
    ProjectResults projResults;
    Vector noErrEst;    // empty

    realize(state, Stage::Position);
    prescribeU(state);
    realize(state, Stage::Velocity);
    projectU(state, noErrEst, projOptions, projResults);
}

void System::projectQ(State& s, Vector& qErrEst, 
             const ProjectOptions& options, ProjectResults& results) const
{   getSystemGuts().projectQ(s,qErrEst,options,results); }
void System::projectU(State& s, Vector& uErrEst, 
             const ProjectOptions& options, ProjectResults& results) const
{   getSystemGuts().projectU(s,uErrEst,options,results); }


// These methods can't be defined until SystemGlobalSubsystem has been
// declared -- they just forward to it.
// TODO: consider moving event stuff to System.

void System::adoptEventHandler(ScheduledEventHandler* handler)
{   updSystemGlobalSubsystem().adoptEventHandler(handler); }
void System::adoptEventHandler(TriggeredEventHandler* handler)
{   updSystemGlobalSubsystem().adoptEventHandler(handler); }
void System::adoptEventReporter(ScheduledEventReporter* reporter)
{   updSystemGlobalSubsystem().adoptEventReporter(reporter); }
void System::adoptEventReporter(TriggeredEventReporter* reporter)
{   updSystemGlobalSubsystem().adoptEventReporter(reporter); }

EventId System::adoptEvent(Event* eventp)
{   return updSystemGlobalSubsystem().adoptEvent(eventp); }

int System::getNumEvents() const
{   return getSystemGlobalSubsystem().getNumEvents(); }

const Event& System::getEvent(EventId id) const
{   return getSystemGlobalSubsystem().getEvent(id); }

Event& System::updEvent(EventId id)
{   return updSystemGlobalSubsystem().updEvent(id); }

bool System::hasEvent(EventId id) const
{   return getSystemGlobalSubsystem().hasEvent(id); }

const Event::Initialization& System::getInitializationEvent() const {
    const EventId eid = getSystemGlobalSubsystem().getInitializationEventId(); 
    return Event::Initialization::downcast(getEvent(eid));
}
const Event::TimeAdvanced& System::getTimeAdvancedEvent() const {
    const EventId eid = getSystemGlobalSubsystem().getTimeAdvancedEventId(); 
    return Event::TimeAdvanced::downcast(getEvent(eid));
}
const Event::Termination& System::getTerminationEvent() const {
    const EventId eid = getSystemGlobalSubsystem().getTerminationEventId(); 
    return Event::Termination::downcast(getEvent(eid));
}
const Event::ExtremeValueIsolated& System::
getExtremeValueIsolatedEvent() const {
    const EventId eid = getSystemGlobalSubsystem()
                                            .getExtremeValueIsolatedEventId(); 
    return Event::ExtremeValueIsolated::downcast(getEvent(eid));
}

Event::Initialization& System::updInitializationEvent() {
    const EventId eid = getSystemGlobalSubsystem().getInitializationEventId(); 
    return Event::Initialization::updDowncast(updEvent(eid));
}
Event::TimeAdvanced& System::updTimeAdvancedEvent() {
    const EventId eid = getSystemGlobalSubsystem().getTimeAdvancedEventId(); 
    return Event::TimeAdvanced::updDowncast(updEvent(eid));
}
Event::Termination& System::updTerminationEvent() {
    const EventId eid = getSystemGlobalSubsystem().getTerminationEventId(); 
    return Event::Termination::updDowncast(updEvent(eid));
}
Event::ExtremeValueIsolated& System::updExtremeValueIsolatedEvent() {
    const EventId eid = getSystemGlobalSubsystem()
                                            .getExtremeValueIsolatedEventId(); 
    return Event::ExtremeValueIsolated::updDowncast(updEvent(eid));
}

EventTriggerId System::adoptEventTrigger(EventTrigger* triggerp)
{   return updSystemGlobalSubsystem().adoptEventTrigger(triggerp); }


int System::getNumEventTriggers() const
{   return getSystemGlobalSubsystem().getNumEventTriggers(); }

const EventTrigger& System::getEventTrigger(EventTriggerId id) const
{   return getSystemGlobalSubsystem().getEventTrigger(id); }

EventTrigger& System::updEventTrigger(EventTriggerId id)
{   return updSystemGlobalSubsystem().updEventTrigger(id); }

bool System::hasEventTrigger(EventTriggerId id) const
{   return getSystemGlobalSubsystem().hasEventTrigger(id); }

const InitializationTrigger& System::getInitializationTrigger() const {
    auto& sgs = getSystemGlobalSubsystem();
    const EventTriggerId tid = sgs.getInitializationTriggerId(); 
    return InitializationTrigger::downcast(sgs.getEventTrigger(tid));
}
const TimeAdvancedTrigger& System::getTimeAdvancedTrigger() const {
    auto& sgs = getSystemGlobalSubsystem();
    const EventTriggerId tid = sgs.getTimeAdvancedTriggerId(); 
    return TimeAdvancedTrigger::downcast(sgs.getEventTrigger(tid));
}
const TerminationTrigger& System::getTerminationTrigger() const {
    auto& sgs = getSystemGlobalSubsystem();
    const EventTriggerId tid = sgs.getTerminationTriggerId(); 
    return TerminationTrigger::downcast(sgs.getEventTrigger(tid));
}

void System::findActiveEventWitnesses
   (const Study&                                        study, 
    Array_<const EventWitness*, ActiveWitnessIndex>&    witnesses) const 
{
    return getSystemGlobalSubsystem()
            .findActiveEventWitnesses(study, witnesses); 
}

void System::findActiveEventTimers
   (const Study&                                    study, 
    Array_<const EventTimer*, ActiveTimerIndex>&    timers) const 
{
    return getSystemGlobalSubsystem()
            .findActiveEventTimers(study, timers); 
}

void System::findNextScheduledEventTimes
   (const Study&    study,
    double          timeOfLastReport,
    double          timeOfLastChange,
    double&         timeOfNextReport,
    EventTriggers&  reportTimers,
    double&         timeOfNextChange,
    EventTriggers&  changeTimers) const
{
    return getSystemGlobalSubsystem().findNextScheduledEventTimes
       (study,timeOfLastReport,timeOfLastChange,
        timeOfNextReport, reportTimers,
        timeOfNextChange, changeTimers);
}

void System::noteEventOccurrence
   (const EventTriggers&    triggers,
    EventsAndCauses&        appendTriggeredEvents,
    Array_<EventId>&        appendIgnoredEvents) const {
    getSystemGlobalSubsystem().noteEventOccurrence
       (triggers,appendTriggeredEvents,appendIgnoredEvents); 
}

void System::performEventReportActions
   (const Study&            study,
    const EventsAndCauses&  triggeredEvents) const {
    getSystemGlobalSubsystem().performEventReportActions 
       (study, triggeredEvents);
}

void System::performEventChangeActions
   (Study&                  study,
    const EventsAndCauses&  triggeredEvents,
    EventChangeResult&      result) const {
    getSystemGlobalSubsystem().performEventChangeActions 
       (study, triggeredEvents, result);
}

System::operator const Subsystem&() const 
{   return getSystemGlobalSubsystem(); }
System::operator Subsystem&() 
{   return updSystemGlobalSubsystem(); }


//==============================================================================
//                             SYSTEM :: GUTS
//==============================================================================
// Implementation of non-inline System::Guts methods.

//------------------------------------------------------------------------------
//                             ADOPT SUBSYSTEM
//------------------------------------------------------------------------------
SubsystemIndex System::Guts::adoptSubsystem(Subsystem& child) {
    assert(child.hasGuts() && !child.isInSystem()); // TODO
    assert(child.isOwnerHandle());

    // This is a topology change.
    invalidateSystemTopologyCache();

    const SubsystemIndex ix = SubsystemIndex((int)m_subsystems.size());
    m_subsystems.resize(ix+1); // grow
    Subsystem& subsys = m_subsystems.back(); // refer to the empty handle

    // Take over ownership of the child's guts, leaving the child
    // as a non-owner reference to the same guts.
    subsys.adoptSubsystemGuts(&child.updSubsystemGuts());
    subsys.setSystem(*m_myHandle, ix);

    return ix;
}


//------------------------------------------------------------------------------
//                            REALIZE TOPOLOGY
//------------------------------------------------------------------------------
// The SystemGlobalSubsystem gets realized first at Topology stage so that
// the resources it allocates are immediately available to the other Subsystems
// in their realizeTopology() methods.
const State& System::Guts::realizeTopology() const {
    if (systemTopologyHasBeenRealized())
        return m_defaultState;

    auto mThis = const_cast<System::Guts*>(this); // mutable temporarily
    State& writableState = mThis->m_defaultState;

    // Initialize the default state.
    writableState.clear();
    writableState.setNumSubsystems(getNumSubsystems());
    for (SubsystemIndex ssx(0); ssx<getNumSubsystems(); ++ssx) 
        writableState.initializeSubsystem(ssx, getSubsystem(ssx).getName(), 
                                          getSubsystem(ssx).getVersion());

    // Realize the SystemGlobalSubsystem (SubsystemIndex 0) first.
    if (writableState.getSubsystemStage(SubsystemIndex(0)) < Stage::Topology)
        getSubsystem(SubsystemIndex(0)).getSubsystemGuts()
                                       .realizeSubsystemTopology(writableState);

    // Allow the concrete System subclass to do its processing.
    realizeTopologyImpl(writableState); // defaultState is mutable

    // Realize any subsystems that the subclass didn't already take care of.
    // (only a subset are typically handled explicitly).
    for (SubsystemIndex ssx(1); ssx<getNumSubsystems(); ++ssx) {
        if (writableState.getSubsystemStage(ssx) < Stage::Topology)
            getSubsystem(ssx).getSubsystemGuts()
                                    .realizeSubsystemTopology(writableState);
    }
    
    // Force the defaultState's Topology stage version number to match the
    // Topology cache version in this System.
    writableState.setSystemTopologyStageVersion
                                    (getSystemTopologyCacheVersion());
    writableState.advanceSystemToStage(Stage::Topology);

    mThis->m_systemTopologyRealized = true;
    ++nRealizationsOfStage[Stage::Topology]; // mutable counter

    // Realize the model using the default settings of the Model variables.
    // This allocates all the later-stage State variables. But the Model stage
    // realization will get re-done if the any Model-stage state variables are
    // changed from these default values.
    realizeModel(writableState);

    return m_defaultState;
}



//------------------------------------------------------------------------------
//                              REALIZE MODEL
//------------------------------------------------------------------------------
// The SystemGlobalSubsystem gets realized first at Model stage so that
// any resources it allocates are immediately available to the other Subsystems
// in their realizeModel() methods.
void System::Guts::realizeModel(State& s) const {
    SimTK_STAGECHECK_TOPOLOGY_REALIZED_ALWAYS(systemTopologyHasBeenRealized(),
        "System", getName(), "System::Guts::realizeModel()");
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Topology, 
        "System::Guts::realizeModel()");
    SimTK_STAGECHECK_TOPOLOGY_VERSION_ALWAYS(
        getSystemTopologyCacheVersion(), s.getSystemTopologyStageVersion(),
        "System", getName(), "System::Guts::realizeModel()");

    if (s.getSystemStage() < Stage::Model) {
        // Realize the SystemGlobalSubsystem (SubsystemIndex 0) first.
        if (s.getSubsystemStage(SubsystemIndex(0)) < Stage::Model)
            getSubsystem(SubsystemIndex(0)).getSubsystemGuts()
                                                    .realizeSubsystemModel(s);

        realizeModelImpl(s);    // Allow the subclass to do its processing.

        // Realize any subsystems that the subclass didn't already take care of.
        for (SubsystemIndex ssx(1); ssx<getNumSubsystems(); ++ssx) {
            if (s.getSubsystemStage(ssx) < Stage::Model)
                getSubsystem(ssx).getSubsystemGuts().realizeSubsystemModel(s);
        }

        s.advanceSystemToStage(Stage::Model);
        ++nRealizationsOfStage[Stage::Model]; // mutable counter
    }
}



//------------------------------------------------------------------------------
//                             REALIZE INSTANCE
//------------------------------------------------------------------------------
void System::Guts::realizeInstance(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Instance).prev(), 
        "System::Guts::realizeInstance()");

    if (s.getSystemStage() < Stage::Instance) {
        realizeInstanceImpl(s); // Allow the subclass to do its processing.

        // Realize any subsystems that the subclass didn't already take care of.
        for (SubsystemIndex ssx(1); ssx<getNumSubsystems(); ++ssx) {
            if (s.getSubsystemStage(ssx) < Stage::Instance)
                getSubsystem(ssx).getSubsystemGuts()
                                                .realizeSubsystemInstance(s);
        }

        // The Default Subsystem (SubsystemIndex 0) gets realized last.
        if (s.getSubsystemStage(SubsystemIndex(0)) < Stage::Instance)
            getSubsystem(SubsystemIndex(0)).getSubsystemGuts()
                                                .realizeSubsystemInstance(s);

        s.advanceSystemToStage(Stage::Instance);
        ++nRealizationsOfStage[Stage::Instance]; // mutable counter
    }
}



//------------------------------------------------------------------------------
//                              REALIZE TIME
//------------------------------------------------------------------------------
void System::Guts::realizeTime(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Time).prev(), 
        "System::Guts::realizeTime()");

    if (s.getSystemStage() < Stage::Time) {
        realizeTimeImpl(s);     // Allow the subclass to do processing.

        // Realize any subsystems that the subclass didn't already take care of.
        for (SubsystemIndex ssx(1); ssx<getNumSubsystems(); ++ssx) {
            if (s.getSubsystemStage(ssx) < Stage::Time)
                getSubsystem(ssx).getSubsystemGuts().realizeSubsystemTime(s);
        }

        // The Default Subsystem (SubsystemIndex 0) gets realized last.
        if (s.getSubsystemStage(SubsystemIndex(0)) < Stage::Time)
            getSubsystem(SubsystemIndex(0)).getSubsystemGuts()
                                                    .realizeSubsystemTime(s);
        
        s.advanceSystemToStage(Stage::Time);
        ++nRealizationsOfStage[Stage::Time]; // mutable counter
    }
}



//------------------------------------------------------------------------------
//                            REALIZE POSITION
//------------------------------------------------------------------------------
void System::Guts::realizePosition(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Position).prev(), 
        "System::Guts::realizePosition()");

    if (s.getSystemStage() < Stage::Position) {
        realizePositionImpl(s);     // Allow the subclass to do processing.

        // Realize any subsystems that the subclass didn't already take care of.
        for (SubsystemIndex ssx(1); ssx<getNumSubsystems(); ++ssx) {
            if (s.getSubsystemStage(ssx) < Stage::Position)
                getSubsystem(ssx).getSubsystemGuts()
                                                .realizeSubsystemPosition(s);
        }

        // The Default Subsystem (SubsystemIndex 0) gets realized last.
        if (s.getSubsystemStage(SubsystemIndex(0)) < Stage::Position)
            getSubsystem(SubsystemIndex(0)).getSubsystemGuts()
                                                .realizeSubsystemPosition(s);

        s.advanceSystemToStage(Stage::Position);
        ++nRealizationsOfStage[Stage::Position]; // mutable counter
    }
}



//------------------------------------------------------------------------------
//                            REALIZE VELOCITY
//------------------------------------------------------------------------------
void System::Guts::realizeVelocity(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Velocity).prev(), 
        "System::Guts::realizeVelocity()");

    if (s.getSystemStage() < Stage::Velocity) {
        realizeVelocityImpl(s);     // Allow the subclass to do processing.

        // Realize any subsystems that the subclass didn't already take care of.
        for (SubsystemIndex ssx(1); ssx<getNumSubsystems(); ++ssx) {
            if (s.getSubsystemStage(ssx) < Stage::Velocity)
                getSubsystem(ssx).getSubsystemGuts()
                                                .realizeSubsystemVelocity(s);
        }

        // The Default Subsystem (SubsystemIndex 0) gets realized last.
        if (s.getSubsystemStage(SubsystemIndex(0)) < Stage::Velocity)
            getSubsystem(SubsystemIndex(0)).getSubsystemGuts()
                                                .realizeSubsystemVelocity(s);

        s.advanceSystemToStage(Stage::Velocity);
        ++nRealizationsOfStage[Stage::Velocity]; // mutable counter
    }
}



//------------------------------------------------------------------------------
//                            REALIZE DYNAMICS
//------------------------------------------------------------------------------
void System::Guts::realizeDynamics(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Dynamics).prev(), 
        "System::Guts::realizeDynamics()");

    if (s.getSystemStage() < Stage::Dynamics) {
        realizeDynamicsImpl(s);     // Allow the subclass to do processing.

        // Realize any subsystems that the subclass didn't already take care of.
        for (SubsystemIndex ssx(1); ssx<getNumSubsystems(); ++ssx) {
            if (s.getSubsystemStage(ssx) < Stage::Dynamics)
                getSubsystem(ssx).getSubsystemGuts()
                                                .realizeSubsystemDynamics(s);
        }

        // The Default Subsystem (SubsystemIndex 0) gets realized last.
        if (s.getSubsystemStage(SubsystemIndex(0)) < Stage::Dynamics)
            getSubsystem(SubsystemIndex(0)).getSubsystemGuts()
                                                .realizeSubsystemDynamics(s);

        s.advanceSystemToStage(Stage::Dynamics);
        ++nRealizationsOfStage[Stage::Dynamics]; // mutable counter
    }
}



//------------------------------------------------------------------------------
//                           REALIZE ACCELERATION
//------------------------------------------------------------------------------
void System::Guts::realizeAcceleration(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), 
                               Stage(Stage::Acceleration).prev(), 
                               "System::Guts::realizeAcceleration()");

    if (s.getSystemStage() < Stage::Acceleration) {
        realizeAccelerationImpl(s); // Allow the subclass to do processing.

        // Realize any subsystems that the subclass didn't already take care of.
        for (SubsystemIndex ssx(1); ssx<getNumSubsystems(); ++ssx) {
            if (s.getSubsystemStage(ssx) < Stage::Acceleration)
                getSubsystem(ssx).getSubsystemGuts()
                                            .realizeSubsystemAcceleration(s);
        }

        // The Default Subsystem (SubsystemIndex 0) gets realized last.
        if (s.getSubsystemStage(SubsystemIndex(0)) < Stage::Acceleration)
            getSubsystem(SubsystemIndex(0)).getSubsystemGuts()
                                            .realizeSubsystemAcceleration(s);

        s.advanceSystemToStage(Stage::Acceleration);
        ++nRealizationsOfStage[Stage::Acceleration]; // mutable counter
    }
}



//------------------------------------------------------------------------------
//                              REALIZE REPORT
//------------------------------------------------------------------------------
void System::Guts::realizeReport(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Report).prev(), 
        "System::Guts::realizeReport()");

    if (s.getSystemStage() < Stage::Report) {
        realizeReportImpl(s);       // Allow the subclass to do processing.

        // Realize any subsystems that the subclass didn't already take care of.
        for (SubsystemIndex ssx(1); ssx<getNumSubsystems(); ++ssx) {
            if (s.getSubsystemStage(ssx) < Stage::Report)
                getSubsystem(ssx).getSubsystemGuts().realizeSubsystemReport(s);
        }

        // The Default Subsystem (SubsystemIndex 0) gets realized last.
        if (s.getSubsystemStage(SubsystemIndex(0)) < Stage::Report)
            getSubsystem(SubsystemIndex(0)).getSubsystemGuts()
                                                .realizeSubsystemReport(s);

        s.advanceSystemToStage(Stage::Report);
        ++nRealizationsOfStage[Stage::Report]; // mutable counter
    }
}



//------------------------------------------------------------------------------
//                    MULTIPLY BY N, ~N, N+, ~(N+)
//------------------------------------------------------------------------------

void System::Guts::multiplyByN(const State& s, const Vector& u, 
                               Vector& dq) const {
    SimTK_STAGECHECK_GE(s.getSystemStage(), Stage::Position,
        "System::Guts::multiplyByN()");
    return multiplyByNImpl(s,u,dq);
}
void System::Guts::multiplyByNTranspose(const State& s, const Vector& fq, 
                                        Vector& fu) const {
    SimTK_STAGECHECK_GE(s.getSystemStage(), Stage::Position,
        "System::Guts::multiplyByNTranspose()");
    return multiplyByNTransposeImpl(s,fq,fu);
}
void System::Guts::multiplyByNPInv(const State& s, const Vector& dq, 
                                   Vector& u) const {
    SimTK_STAGECHECK_GE(s.getSystemStage(), Stage::Position,
        "System::Guts::multiplyByNPInv()");
    return multiplyByNPInvImpl(s,dq,u);
}
void System::Guts::multiplyByNPInvTranspose(const State& s, const Vector& fu, 
                                            Vector& fq) const {
    SimTK_STAGECHECK_GE(s.getSystemStage(), Stage::Position,
        "System::Guts::multiplyByNPInvTranspose()");
    return multiplyByNPInvTransposeImpl(s,fu,fq);
}



//------------------------------------------------------------------------------
//                              PRESCRIBE Q
//------------------------------------------------------------------------------
bool System::Guts::prescribeQ(State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Time,
                               "System::Guts::prescribeQ()");
    ++nPrescribeQCalls; // mutable counter
    return prescribeQImpl(s);
}


//------------------------------------------------------------------------------
//                              PRESCRIBE U
//------------------------------------------------------------------------------
bool System::Guts::prescribeU(State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Position,
                               "System::Guts::prescribeU()");
    ++nPrescribeUCalls; // mutable counter
    return prescribeUImpl(s);
}


//------------------------------------------------------------------------------
//                          GET FREE Q(U) INDEX
//------------------------------------------------------------------------------
void System::Guts::
getFreeQIndex(const State& s, Array_<SystemQIndex>& freeQs) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Instance,
                               "System::Guts::getFreeQIndex()");
    return getFreeQIndexImpl(s, freeQs);
}
void System::Guts::
getFreeUIndex(const State& s, Array_<SystemUIndex>& freeUs) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Instance,
                               "System::Guts::getFreeUIndex()");
    return getFreeUIndexImpl(s, freeUs);
}



//------------------------------------------------------------------------------
//                                PROJECT Q
//------------------------------------------------------------------------------
void System::Guts::projectQ(State& s, Vector& qErrEst, 
                            const ProjectOptions& options, 
                            ProjectResults& results) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Position,
                               "System::Guts::projectQ()");

    nProjectQCalls++;       // counters are mutable
    nFailedProjectQCalls++; // assume this will throw an exception
    //---------------------------------------------------------
    projectQImpl(s,qErrEst,options,results);
    //---------------------------------------------------------
    if (results.getExitStatus()==ProjectResults::Succeeded) {
        nFailedProjectQCalls--; // never mind!
        if (results.getAnyChangeMade()) ++nQProjections;
        if (qErrEst.size()) ++nQErrEstProjections;
    }
}



//------------------------------------------------------------------------------
//                                PROJECT U
//------------------------------------------------------------------------------
void System::Guts::projectU(State& s, Vector& uErrEst, 
                            const ProjectOptions& options, 
                            ProjectResults& results) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Velocity,
                               "System::Guts::projectU()");

    nProjectUCalls++;       // counters are mutable
    nFailedProjectUCalls++; // assume this will throw an exception
    //---------------------------------------------------------
    projectUImpl(s,uErrEst,options,results);
    //---------------------------------------------------------
    if (results.getExitStatus()==ProjectResults::Succeeded) {
        nFailedProjectUCalls--; // never mind!
        if (results.getAnyChangeMade()) ++nUProjections;
        if (uErrEst.size()) ++nUErrEstProjections;
    }
}



//------------------------------------------------------------------------------
//                                 REALIZE
//------------------------------------------------------------------------------
void System::Guts::realize(const State& s, Stage g) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Model, 
        "System::Guts::realize()");

    Stage stageNow = Stage::Empty;
    while ((stageNow=s.getSystemStage()) < g) {
        switch (stageNow) {
        case Stage::Model:        realizeInstance(s);     break;
        case Stage::Instance:     realizeTime(s);         break;
        case Stage::Time:         realizePosition(s);     break;
        case Stage::Position:     realizeVelocity(s);     break;
        case Stage::Velocity:     realizeDynamics(s);     break;
        case Stage::Dynamics:     realizeAcceleration(s); break;
        case Stage::Acceleration: realizeReport(s);       break;
        default: assert(!"System::Guts::realize(): bad stage");
        }
    }
}



//------------------------------------------------------------------------------
//                   CALC DECORATIVE GEOMETRY AND APPEND
//------------------------------------------------------------------------------
void System::Guts::calcDecorativeGeometryAndAppend
   (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const 
{
    assert(stage==Stage::Topology || s.getSystemStage() >= stage);
    for (SubsystemIndex i(0); i<getNumSubsystems(); ++i)
        m_subsystems[i].getSubsystemGuts()
                              .calcDecorativeGeometryAndAppend(s, stage, geom);
}


//------------------------------------------------------------------------------
// Default implementations of virtual methods. Most quietly do nothing. A
// concrete System will typically override these.
//------------------------------------------------------------------------------


// These multiplyByN methods must be implemented if used.
void System::Guts::multiplyByNImpl
   (const State& state, const Vector& u, Vector& dq) const
{   SimTK_THROW2(Exception::UnimplementedVirtualMethod, "System::Guts", 
                 "multiplyByNImpl"); }
void System::Guts::multiplyByNTransposeImpl
   (const State& state, const Vector& fq, Vector& fu) const
{   SimTK_THROW2(Exception::UnimplementedVirtualMethod, "System::Guts", 
                 "multiplyByNTransposeImpl"); }
void System::Guts::multiplyByNPInvImpl
   (const State& state, const Vector& dq, Vector& u) const 
{   SimTK_THROW2(Exception::UnimplementedVirtualMethod, "System::Guts", 
                 "multiplyByNPInvImpl"); }
void System::Guts::multiplyByNPInvTransposeImpl
   (const State& state, const Vector& fu, Vector& fq) const 
{   SimTK_THROW2(Exception::UnimplementedVirtualMethod, "System::Guts", 
                 "multiplyByNPInvTransposeImpl"); }


