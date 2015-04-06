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
 * Implementation of System, System::Guts, and System::GutsRep, and also
 * EventTriggerInfo and EventTriggerInfoRep.
 */

#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"
#include "SimTKcommon/internal/Subsystem.h"
#include "SimTKcommon/internal/SubsystemGuts.h"
#include "SimTKcommon/internal/System.h"
#include "SimTKcommon/internal/SystemGuts.h"
#include "SimTKcommon/internal/EventHandler.h"
#include "SimTKcommon/internal/EventReporter.h"

#include "SystemGutsRep.h"

#include <cassert>
#include <map>
#include <set>

namespace SimTK {

//==============================================================================
//                                   SYSTEM
//==============================================================================
// Most of the System:: methods just forward to System::Guts.

bool System::isEmptyHandle() const {return guts==0;}
bool System::isOwnerHandle() const {return guts==0 || &guts->getSystem()==this;}
bool System::isSameSystem(const System& otherSystem) const {
    return guts && (guts==otherSystem.guts);
}

System::System(const System& src) : guts(0) {
    if (src.guts) {
        guts = src.guts->clone();
        guts->setOwnerHandle(*this);
    }
}

System& System::operator=(const System& src) {
    if (!isSameSystem(src)) {
        if (isOwnerHandle())
            delete guts;
        guts=0;
        if (src.guts) {
            guts = src.guts->clone();
            guts->setOwnerHandle(*this);
        }
    }
    return *this;
}

System::~System() {
    if (guts && isOwnerHandle())
        delete guts;
    guts=0;
}

void System::adoptSystemGuts(System::Guts* g) {
    SimTK_ASSERT_ALWAYS(g, "System::adoptSystemGuts(): can't adopt null Guts.");
    SimTK_ASSERT_ALWAYS(!guts,
    "System::adoptSystemGuts(): this System handle is already in use.");
    SimTK_ASSERT_ALWAYS(g->getNumSubsystems() == 0,
    "System::adoptSystemGuts(): can't adopt Guts that already has Subsystems.");
    guts = g;
    guts->setOwnerHandle(*this);

    // Now add the Default Subsystem, which will be SubsystemIndex(0).
    DefaultSystemSubsystem defsub(*this); // invokes adoptSubsystem().

}

const String& System::getName()    const {return getSystemGuts().getName();}
const String& System::getVersion() const {return getSystemGuts().getVersion();}

System& System::setDefaultTimeScale(Real tc) 
{   updSystemGuts().updRep().setDefaultTimeScale(tc); return *this; }
Real System::getDefaultTimeScale() const
{   return getSystemGuts().getRep().getDefaultTimeScale(); }
System& System::setDefaultLengthScale(Real lc)
{   updSystemGuts().updRep().setDefaultLengthScale(lc); return *this; }
Real System::getDefaultLengthScale() const
{   return getSystemGuts().getRep().getDefaultLengthScale(); }

System& System::setUpDirection(const CoordinateDirection& up) 
{   updSystemGuts().updRep().setUpDirection(up); return *this; }
CoordinateDirection System::getUpDirection() const
{   return getSystemGuts().getRep().getUpDirection(); }
System& System::setUseUniformBackground(bool useUniformBackground)
{   updSystemGuts().updRep().setUseUniformBackground(useUniformBackground);
    return *this; }
bool System::getUseUniformBackground() const
{   return getSystemGuts().getRep().getUseUniformBackground(); }

void System::resetAllCountersToZero() {updSystemGuts().updRep().resetAllCounters();}
int System::getNumRealizationsOfThisStage(Stage g) const {return getSystemGuts().getRep().nRealizationsOfStage[g];}
int System::getNumRealizeCalls() const {return getSystemGuts().getRep().nRealizeCalls;}

int System::getNumPrescribeQCalls() const {return getSystemGuts().getRep().nPrescribeQCalls;}
int System::getNumPrescribeUCalls() const {return getSystemGuts().getRep().nPrescribeUCalls;}

int System::getNumProjectQCalls() const {return getSystemGuts().getRep().nProjectQCalls;}
int System::getNumProjectUCalls() const {return getSystemGuts().getRep().nProjectUCalls;}
int System::getNumFailedProjectQCalls() const {return getSystemGuts().getRep().nFailedProjectQCalls;}
int System::getNumFailedProjectUCalls() const {return getSystemGuts().getRep().nFailedProjectUCalls;}
int System::getNumQProjections() const {return getSystemGuts().getRep().nQProjections;}
int System::getNumUProjections() const {return getSystemGuts().getRep().nUProjections;} 
int System::getNumQErrorEstimateProjections() const {return getSystemGuts().getRep().nQErrEstProjections;}
int System::getNumUErrorEstimateProjections() const {return getSystemGuts().getRep().nUErrEstProjections;}

int System::getNumHandlerCallsThatChangedStage(Stage g) const {return getSystemGuts().getRep().nHandlerCallsThatChangedStage[g];}
int System::getNumHandleEventCalls() const {return getSystemGuts().getRep().nHandleEventsCalls;}
int System::getNumReportEventCalls() const {return getSystemGuts().getRep().nReportEventsCalls;}

const State& System::getDefaultState() const {return getSystemGuts().getDefaultState();}
State& System::updDefaultState() {return updSystemGuts().updDefaultState();}

SubsystemIndex System::adoptSubsystem(Subsystem& child) 
{   return updSystemGuts().adoptSubsystem(child); }

int System::getNumSubsystems() const 
{   return getSystemGuts().getNumSubsystems(); }

const Subsystem& System::getSubsystem(SubsystemIndex i) const 
{   return getSystemGuts().getSubsystem(i); }
Subsystem& System::updSubsystem(SubsystemIndex i) 
{   return updSystemGuts().updSubsystem(i); }

const DefaultSystemSubsystem& System::getDefaultSubsystem() const {
    return static_cast<const DefaultSystemSubsystem&>
                            (getSystemGuts().getSubsystem(SubsystemIndex(0)));
}
DefaultSystemSubsystem& System::updDefaultSubsystem() {
    return static_cast<DefaultSystemSubsystem&>
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

void System::handleEvents
   (State& s, Event::Cause cause, const Array_<EventId>& eventIds,
    const HandleEventsOptions& options, HandleEventsResults& results) const
{   getSystemGuts().handleEvents(s,cause,eventIds,options,results); }

void System::reportEvents
   (const State& s, Event::Cause cause, const Array_<EventId>& eventIds) const
{   getSystemGuts().reportEvents(s,cause,eventIds); }

void System::
calcEventTriggerInfo(const State& s, Array_<EventTriggerInfo>& info) const
{   getSystemGuts().calcEventTriggerInfo(s,info); }

void System::calcTimeOfNextScheduledEvent
   (const State& s, Real& tNextEvent,
    Array_<EventId>& eventIds, bool includeCurrentTime) const
{   getSystemGuts().calcTimeOfNextScheduledEvent
                                (s,tNextEvent,eventIds,includeCurrentTime); }

void System::calcTimeOfNextScheduledReport
   (const State& s, Real& tNextEvent,
    Array_<EventId>& eventIds, bool includeCurrentTime) const
{   getSystemGuts().calcTimeOfNextScheduledReport
                                (s,tNextEvent,eventIds,includeCurrentTime); }

EventId System::createNewEventId(SubsystemIndex ssx) const {
    return getSystemGuts().createNewEventId(ssx);
}

SubsystemIndex System::findEventIdOwner(EventId id) const {
    return getSystemGuts().findEventIdOwner(id);
}


//==============================================================================
//                             SYSTEM :: GUTS
//==============================================================================

// This is also the default constructor.
System::Guts::Guts(const String& name, const String& version) {
    rep = new GutsRep(name,version);
    // note that the GutsRep object currently has no owner handle
}

System::Guts::~Guts() {
    delete rep;
    rep=0;
}

// Copy constructor
System::Guts::Guts(const Guts& src) : rep(0) {
    if (src.rep) {
        rep = new GutsRep(*src.rep);
        // note that the GutsRep object currently has no owner handle
    }
}

// Copy assignment is suppressed
    

const System& System::Guts::getSystem() const {
    assert(rep->myHandle);
    return *rep->myHandle;
}

System& System::Guts::updSystem() {
    assert(rep->myHandle);
    return *rep->myHandle;
}

void System::Guts::setOwnerHandle(System& sys) {
    assert(!rep->myHandle);
    rep->myHandle = &sys;
}

bool System::Guts::hasOwnerHandle() const {
    return rep->myHandle != 0;
}

void System::Guts::setHasTimeAdvancedEvents(bool hasEm) {
    updRep().hasTimeAdvancedEventsFlag = hasEm;
}
bool System::Guts::hasTimeAdvancedEvents() const {
    return getRep().hasTimeAdvancedEventsFlag;
}

const String& System::Guts::getName()    const {return getRep().getName();}
const String& System::Guts::getVersion() const {return getRep().getVersion();}

bool System::Guts::systemTopologyHasBeenRealized() const 
{   return getRep().systemTopologyHasBeenRealized(); }
StageVersion System::Guts::getSystemTopologyCacheVersion() const
{   return getRep().getSystemTopologyCacheVersion(); }
void System::Guts::setSystemTopologyCacheVersion(StageVersion topoVersion) const
{   getRep().setSystemTopologyCacheVersion(topoVersion); }
void System::Guts::invalidateSystemTopologyCache() const 
{   return getRep().invalidateSystemTopologyCache(); } // mutable

const State& System::Guts::getDefaultState() const {
    SimTK_STAGECHECK_TOPOLOGY_REALIZED_ALWAYS(systemTopologyHasBeenRealized(),
        "System", getName(), "System::Guts::getDefaultState()");

    return getRep().getDefaultState();
}

State& System::Guts::updDefaultState() {
    SimTK_STAGECHECK_TOPOLOGY_REALIZED_ALWAYS(systemTopologyHasBeenRealized(),
        "System", getName(), "System::Guts::updDefaultState()");

    return updRep().updDefaultState();
}

int System::Guts::getNumSubsystems() const {return getRep().getNumSubsystems();}
const Subsystem& System::Guts::getSubsystem(SubsystemIndex i) const 
{   return getRep().getSubsystem(i); }
Subsystem& System::Guts::updSubsystem(SubsystemIndex i) 
{   return updRep().updSubsystem(i); }

System::Guts* System::Guts::clone() const 
{   return cloneImpl(); }



//------------------------------------------------------------------------------
//                            REALIZE TOPOLOGY
//------------------------------------------------------------------------------
const State& System::Guts::realizeTopology() const {
    if (getRep().systemTopologyHasBeenRealized())
        return getRep().defaultState;

    auto mThis = const_cast<System::Guts*>(this);
    System::Guts::GutsRep& rep = mThis->updRep(); // mutable temporarily

    // Initialize EventIds and map.
    rep.nextAvailableEventId = EventId(1);
    rep.eventOwnerMap.clear();

    // Initialize the default state.
    State& defaultState = rep.defaultState;
    defaultState.clear();
    defaultState.setNumSubsystems(getNumSubsystems());
    for (SubsystemIndex ssx(0); ssx<getNumSubsystems(); ++ssx) 
        defaultState.initializeSubsystem(ssx, getSubsystem(ssx).getName(), 
                                            getSubsystem(ssx).getVersion());

    // Allow the concrete System subclass to do its processing.
    realizeTopologyImpl(defaultState); // defaultState is mutable

    // Realize any subsystems that the subclass didn't already take care of.
    // (only a subset are typically handled explicitly).
    for (SubsystemIndex ssx(1); ssx<getNumSubsystems(); ++ssx) {
        if (defaultState.getSubsystemStage(ssx) < Stage::Topology)
            getSubsystem(ssx).getSubsystemGuts()
                                    .realizeSubsystemTopology(defaultState);
    }

    // The Default Subsystem (SubsystemIndex 0) always gets realized last.
    if (defaultState.getSubsystemStage(SubsystemIndex(0)) < Stage::Topology)
        getSubsystem(SubsystemIndex(0)).getSubsystemGuts()
                                    .realizeSubsystemTopology(defaultState);
       
    // Force the defaultState's Topology stage version number to match the
    // Topology cache version in this System.
    defaultState.setSystemTopologyStageVersion
                                    (getRep().getSystemTopologyCacheVersion());
    defaultState.advanceSystemToStage(Stage::Topology);

    rep.systemTopologyRealized = true;
    ++rep.nRealizationsOfStage[Stage::Topology]; // mutable counter

    // Realize the model using the default settings of the Model variables.
    // This allocates all the later-stage State variables.
    realizeModel(defaultState);

    return defaultState;
}



//------------------------------------------------------------------------------
//                              REALIZE MODEL
//------------------------------------------------------------------------------
void System::Guts::realizeModel(State& s) const {
    SimTK_STAGECHECK_TOPOLOGY_REALIZED_ALWAYS(systemTopologyHasBeenRealized(),
        "System", getName(), "System::Guts::realizeModel()");
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Topology, 
        "System::Guts::realizeModel()");
    SimTK_STAGECHECK_TOPOLOGY_VERSION_ALWAYS(
        getSystemTopologyCacheVersion(), s.getSystemTopologyStageVersion(),
        "System", getName(), "System::Guts::realizeModel()");

    if (s.getSystemStage() < Stage::Model) {
        realizeModelImpl(s);    // Allow the subclass to do its processing.

        // Realize any subsystems that the subclass didn't already take care of.
        for (SubsystemIndex ssx(1); ssx<getNumSubsystems(); ++ssx) {
            if (s.getSubsystemStage(ssx) < Stage::Model)
                getSubsystem(ssx).getSubsystemGuts().realizeSubsystemModel(s);
        }

        // The Default Subsystem (SubsystemIndex 0) always gets realized last.
        if (s.getSubsystemStage(SubsystemIndex(0)) < Stage::Model)
            getSubsystem(SubsystemIndex(0)).getSubsystemGuts()
                                                    .realizeSubsystemModel(s);

        s.advanceSystemToStage(Stage::Model);
        ++getRep().nRealizationsOfStage[Stage::Model]; // mutable counter
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

        // The Default Subsystem (SubsystemIndex 0) always gets realized last.
        if (s.getSubsystemStage(SubsystemIndex(0)) < Stage::Instance)
            getSubsystem(SubsystemIndex(0)).getSubsystemGuts()
                                                .realizeSubsystemInstance(s);

        s.advanceSystemToStage(Stage::Instance);
        ++getRep().nRealizationsOfStage[Stage::Instance]; // mutable counter
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

        // The Default Subsystem (SubsystemIndex 0) always gets realized last.
        if (s.getSubsystemStage(SubsystemIndex(0)) < Stage::Time)
            getSubsystem(SubsystemIndex(0)).getSubsystemGuts()
                                                    .realizeSubsystemTime(s);
        
        s.advanceSystemToStage(Stage::Time);
        ++getRep().nRealizationsOfStage[Stage::Time]; // mutable counter
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

        // The Default Subsystem (SubsystemIndex 0) always gets realized last.
        if (s.getSubsystemStage(SubsystemIndex(0)) < Stage::Position)
            getSubsystem(SubsystemIndex(0)).getSubsystemGuts()
                                                .realizeSubsystemPosition(s);

        s.advanceSystemToStage(Stage::Position);
        ++getRep().nRealizationsOfStage[Stage::Position]; // mutable counter
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

        // The Default Subsystem (SubsystemIndex 0) always gets realized last.
        if (s.getSubsystemStage(SubsystemIndex(0)) < Stage::Velocity)
            getSubsystem(SubsystemIndex(0)).getSubsystemGuts()
                                                .realizeSubsystemVelocity(s);

        s.advanceSystemToStage(Stage::Velocity);
        ++getRep().nRealizationsOfStage[Stage::Velocity]; // mutable counter
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

        // The Default Subsystem (SubsystemIndex 0) always gets realized last.
        if (s.getSubsystemStage(SubsystemIndex(0)) < Stage::Dynamics)
            getSubsystem(SubsystemIndex(0)).getSubsystemGuts()
                                                .realizeSubsystemDynamics(s);

        s.advanceSystemToStage(Stage::Dynamics);
        ++getRep().nRealizationsOfStage[Stage::Dynamics]; // mutable counter
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

        // The Default Subsystem (SubsystemIndex 0) always gets realized last.
        if (s.getSubsystemStage(SubsystemIndex(0)) < Stage::Acceleration)
            getSubsystem(SubsystemIndex(0)).getSubsystemGuts()
                                            .realizeSubsystemAcceleration(s);

        s.advanceSystemToStage(Stage::Acceleration);
        ++getRep().nRealizationsOfStage[Stage::Acceleration]; // mutable counter
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

        // The Default Subsystem (SubsystemIndex 0) always gets realized last.
        if (s.getSubsystemStage(SubsystemIndex(0)) < Stage::Report)
            getSubsystem(SubsystemIndex(0)).getSubsystemGuts()
                                                .realizeSubsystemReport(s);

        s.advanceSystemToStage(Stage::Report);
        ++getRep().nRealizationsOfStage[Stage::Report]; // mutable counter
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
    getRep().nPrescribeQCalls++; // mutable counter
    return prescribeQImpl(s);
}


//------------------------------------------------------------------------------
//                              PRESCRIBE U
//------------------------------------------------------------------------------
bool System::Guts::prescribeU(State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Position,
                               "System::Guts::prescribeU()");
    getRep().nPrescribeUCalls++; // mutable counter
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
    const System::Guts::GutsRep& rep = getRep();

    rep.nProjectQCalls++;       // counters are mutable
    rep.nFailedProjectQCalls++; // assume this will throw an exception
    //---------------------------------------------------------
    projectQImpl(s,qErrEst,options,results);
    //---------------------------------------------------------
    if (results.getExitStatus()==ProjectResults::Succeeded) {
        rep.nFailedProjectQCalls--; // never mind!
        if (results.getAnyChangeMade()) rep.nQProjections++;
        if (qErrEst.size()) rep.nQErrEstProjections++;
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
    const System::Guts::GutsRep& rep = getRep();

    rep.nProjectUCalls++;       // counters are mutable
    rep.nFailedProjectUCalls++; // assume this will throw an exception
    //---------------------------------------------------------
    projectUImpl(s,uErrEst,options,results);
    //---------------------------------------------------------
    if (results.getExitStatus()==ProjectResults::Succeeded) {
        rep.nFailedProjectUCalls--; // never mind!
        if (results.getAnyChangeMade()) rep.nUProjections++;
        if (uErrEst.size()) rep.nUErrEstProjections++;
    }
}



//------------------------------------------------------------------------------
//                             HANDLE EVENTS
//------------------------------------------------------------------------------
void System::Guts::handleEvents
   (State& s, Event::Cause cause, const Array_<EventId>& eventIds,
    const HandleEventsOptions& options, HandleEventsResults& results) const
{
         // TODO: is Model the right stage?
   SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Model, 
        "System::Guts::handleEvents()");
    const Real savedTime = s.getTime();

    // Save the stage version numbers so we can look for changes.
    Array_<StageVersion> stageVersions;
    s.getSystemStageVersions(stageVersions);

    handleEventsImpl(s, cause, eventIds, options, results);

    // Note the lowest stage whose version was changed by the handler.
    const Stage lowestModified = s.getLowestSystemStageDifference(stageVersions);
    results.setLowestModifiedStage(lowestModified);

    getRep().nHandleEventsCalls++; // mutable counters
    getRep().nHandlerCallsThatChangedStage[lowestModified]++;

    SimTK_ASSERT_ALWAYS(s.getTime() == savedTime,
        "System::Guts::handleEvents(): handleEventsImpl() tried to change the time");
}



//------------------------------------------------------------------------------
//                              REPORT EVENTS
//------------------------------------------------------------------------------
void System::Guts::reportEvents
   (const State& s, Event::Cause cause, const Array_<EventId>& eventIds) const
{
        // TODO: is Model this the right stage?
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Model, 
        "System::Guts::reportEvents()");
    reportEventsImpl(s,cause,eventIds);

    getRep().nReportEventsCalls++; // mutable counter
    getRep().nHandlerCallsThatChangedStage[Stage::Report]++;
}



//------------------------------------------------------------------------------
//                     CALC TIME OF NEXT SCHEDULED EVENT
//------------------------------------------------------------------------------
void System::Guts::calcTimeOfNextScheduledEvent
    (const State& s, Real& tNextEvent, Array_<EventId>& eventIds, 
     bool includeCurrentTime) const
{
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Time,
        "System::Guts::calcTimeOfNextScheduledEvent()");
    tNextEvent = Infinity;
    calcTimeOfNextScheduledEventImpl(s,tNextEvent,eventIds,includeCurrentTime);
}



//------------------------------------------------------------------------------
//                    CALC TIME OF NEXT SCHEDULED REPORT
//------------------------------------------------------------------------------
void System::Guts::calcTimeOfNextScheduledReport
   (const State& s, Real& tNextEvent, Array_<EventId>& eventIds, 
    bool includeCurrentTime) const
{
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Time,
        "System::Guts::calcTimeOfNextScheduledReport()");
    tNextEvent = Infinity;
    calcTimeOfNextScheduledReportImpl(s,tNextEvent,eventIds,includeCurrentTime);
}



//------------------------------------------------------------------------------
//                          CALC EVENT TRIGGER INFO
//------------------------------------------------------------------------------
void System::Guts::calcEventTriggerInfo
   (const State& s, Array_<EventTriggerInfo>& info) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Instance,
        "System::Guts::calcEventTriggerInfo()");
    calcEventTriggerInfoImpl(s,info);
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
        getRep().subsystems[i].getSubsystemGuts()
                              .calcDecorativeGeometryAndAppend(s, stage, geom);
}


SubsystemIndex System::Guts::adoptSubsystem(Subsystem& src) {
    return updRep().adoptSubsystem(src);
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


//------------------------------------------------------------------------------
//                          HANDLE EVENTS IMPL
//------------------------------------------------------------------------------
// This is the default implementation but shouldn't need to be overridden. It
// distributes handleEvents() calls to each relevant Subsystem, in order of
// SubsystemIndex.
void System::Guts::handleEventsImpl
   (State& state, Event::Cause cause, const Array_<EventId>& eventIds,
    const HandleEventsOptions& options, HandleEventsResults& results) const
{
    // Event handlers should not be able to modify the time. We'll save time 
    // here and restore it at the end if it might have been changed.
    const Real savedT = state.getTime();

    // Note: caller will take care of lowest modified stage.
    results.setExitStatus(HandleEventsResults::Succeeded); // assume success

    // If we have been supplied a list of eventIds, group them by Subsystem
    // ownership and only invoke handleEvents() on Subsystems that own one or
    // more of the eventIds. If there are no eventIds then this is a generic
    // event like Initialization and all the Subsystems get a call.
    
    const bool isGenericEvent = eventIds.empty();
    Array_<Array_<EventId>, SubsystemIndex> eventsForSubsystem;
    if (!isGenericEvent) {
        eventsForSubsystem.resize(getNumSubsystems());
        for (EventId id : eventIds) {
            const SubsystemIndex owner = findEventIdOwner(id);
            if (!owner.isValid()) continue; // ignore unrecognized events
            eventsForSubsystem[owner].push_back(id);
        }
    }

    for (SubsystemIndex sx(0); sx < getNumSubsystems(); ++sx) {
        if (!isGenericEvent && eventsForSubsystem[sx].empty())
            continue; // No events for this Subsystem

        const Subsystem::Guts& sub = getRep().subsystems[sx].getSubsystemGuts();
        HandleEventsResults subsysResults;
        //---------------------------------------------------------------------
        sub.handleEvents(state, cause, 
                         isGenericEvent ? eventIds /*0-length*/
                                        : eventsForSubsystem[sx],
                         options, subsysResults);
        //---------------------------------------------------------------------

        if (subsysResults.getAnyChangeMade())
            results.setAnyChangeMade(true);

        if (subsysResults.getExitStatus() == HandleEventsResults::Failed) {
            results.setExitStatus(HandleEventsResults::Failed);
            results.setMessage(subsysResults.getMessage());
            break; // must quit at first failure
        }

        if (subsysResults.getExitStatus()==HandleEventsResults::ShouldTerminate)
            results.setExitStatus(HandleEventsResults::ShouldTerminate);
            // termination required, but keep going
    }

    // Event handlers must not change the time.
    if (state.getTime() != savedT)
        state.setTime(savedT);
}



//------------------------------------------------------------------------------
//                            REPORT EVENTS IMPL
//------------------------------------------------------------------------------
// This is the default implementation but shouldn't need to be overridden. It
// distributes reportEvents() calls to each relevant Subsystem, in order of
// SubsystemIndex.
int System::Guts::reportEventsImpl
   (const State& s, Event::Cause cause, const Array_<EventId>& eventIds) const
{
    // If we have been supplied a list of eventIds, group them by Subsystem
    // ownership and only invoke reportEvents() on Subsystems that own one or
    // more of the eventIds. If there are no eventIds then this is a generic
    // reporting event and all the Subsystems get a call.
    
    const bool isGenericEvent = eventIds.empty();
    Array_<Array_<EventId>, SubsystemIndex> eventsForSubsystem;
    if (!isGenericEvent) {
        eventsForSubsystem.resize(getNumSubsystems());
        for (EventId id : eventIds) {
            const SubsystemIndex owner = findEventIdOwner(id);
            if (!owner.isValid()) continue; // ignore unrecognized events
            eventsForSubsystem[owner].push_back(id);
        }
    }

    for (SubsystemIndex sx(0); sx < getNumSubsystems(); ++sx) {
        if (!isGenericEvent && eventsForSubsystem[sx].empty())
            continue; // No events for this Subsystem

        const Subsystem::Guts& sub = getRep().subsystems[sx].getSubsystemGuts();

        //---------------------------------------------------------------------
        sub.reportEvents(s, cause, 
                         isGenericEvent ? eventIds /*0-length*/
                                        : eventsForSubsystem[sx]);
        //---------------------------------------------------------------------
    }
    return 0;
}



//------------------------------------------------------------------------------
//                       CALC EVENT TRIGGER INFO IMPL
//------------------------------------------------------------------------------
// This is the default implementation but shouldn't need to be overridden.
int System::Guts::calcEventTriggerInfoImpl
   (const State& s, Array_<EventTriggerInfo>& info) const {

    // Loop over each subsystem, get its EventTriggerInfos, and combine all of 
    // them into a single list.
    
    info.clear();
    for (SubsystemIndex sx(0); sx < getNumSubsystems(); ++sx) {
        const Subsystem::Guts& sub = getRep().subsystems[sx].getSubsystemGuts();
        Array_<EventTriggerInfo> subinfo;
        sub.calcEventTriggerInfo(s, subinfo);
        for (auto e = subinfo.begin(); e != subinfo.end(); ++e)
            info.push_back(*e);
    }
    return 0;
}



//------------------------------------------------------------------------------
//                    CALC TIME OF NEXT SCHEDULED EVENT IMPL
//------------------------------------------------------------------------------
// This is the default implementation but shouldn't need to be overridden.
int System::Guts::calcTimeOfNextScheduledEventImpl
   (const State& s, Real& tNextEvent, Array_<EventId>& eventIds, 
    bool includeCurrentTime) const
{
    tNextEvent = Infinity;
    eventIds.clear();
    Array_<EventId> ids;
    for (SubsystemIndex sx(0); sx < getNumSubsystems(); ++sx) {
        const Subsystem::Guts& sub = getRep().subsystems[sx].getSubsystemGuts();
        Real time;
        ids.clear();
        sub.calcTimeOfNextScheduledEvent(s, time, ids, includeCurrentTime);
        if (time <= tNextEvent) {
            tNextEvent = time;
            if (time < tNextEvent) 
                eventIds.clear(); // otherwise just accumulate
            for (int i = 0; i < (int)ids.size(); ++i)
                eventIds.push_back(ids[i]);
        }
    }
    return 0;
}



//------------------------------------------------------------------------------
//                   CALC TIME OF NEXT SCHEDULED REPORT IMPL
//------------------------------------------------------------------------------
// This is the default implementation but shouldn't need to be overridden.
int System::Guts::calcTimeOfNextScheduledReportImpl
   (const State& s, Real& tNextEvent, Array_<EventId>& eventIds, 
    bool includeCurrentTime) const
{
    tNextEvent = Infinity;
    eventIds.clear();
    Array_<EventId> ids;
    for (SubsystemIndex sx(0); sx < getNumSubsystems(); ++sx) {
        const Subsystem::Guts& sub = getRep().subsystems[sx].getSubsystemGuts();
        Real time;
        ids.clear();
        sub.calcTimeOfNextScheduledReport(s, time, ids, includeCurrentTime);
        if (time <= tNextEvent) {
            tNextEvent = time;
            if (time < tNextEvent) 
                eventIds.clear(); // otherwise just accumulate
            for (int i = 0; i < (int)ids.size(); ++i)
                eventIds.push_back(ids[i]);
        }
    }
    return 0;
}

EventId System::Guts::createNewEventId(SubsystemIndex ssx) const {
    return getRep().createNewEventId(ssx);
}

SubsystemIndex System::Guts::findEventIdOwner(EventId id) const {
    return getRep().findEventIdOwner(id);
}

//==============================================================================
//                         SYSTEM :: GUTS :: GUTS REP
//==============================================================================
// All inline currently.



//==============================================================================
//                         DEFAULT SYSTEM SUBSYSTEM
//==============================================================================

class DefaultSystemSubsystem::Guts : public Subsystem::Guts {
public:

    Guts() : Subsystem::Guts("DefaultSystemSubsystem::Guts", "0.0.1") { }
    
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
        auto mThis = const_cast<DefaultSystemSubsystem::Guts*>(this);
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
friend class DefaultSystemSubsystem;

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

DefaultSystemSubsystem::DefaultSystemSubsystem(System& sys) {
    adoptSubsystemGuts(new DefaultSystemSubsystem::Guts());
    sys.adoptSubsystem(*this);
}

const DefaultSystemSubsystem::Guts& DefaultSystemSubsystem::getGuts() const {
    return SimTK_DYNAMIC_CAST_DEBUG<const DefaultSystemSubsystem::Guts&>
                                                        (getSubsystemGuts());
}

DefaultSystemSubsystem::Guts& DefaultSystemSubsystem::updGuts() {
    return SimTK_DYNAMIC_CAST_DEBUG<DefaultSystemSubsystem::Guts&>
                                                        (updSubsystemGuts());
}

/*
 * Add a ScheduledEventHandler to the System.  This must be called before the 
 * Model stage is realized.
 * 
 * The System assumes ownership of the object passed to this method, and will 
 * delete it when the System is deleted.
 */
void DefaultSystemSubsystem::addEventHandler(ScheduledEventHandler* handler) {
    updGuts().updScheduledEventHandlers().push_back(handler);
}

/*
 * Add a TriggeredEventHandler to the System.  This must be called before the 
 * Model stage is realized.
 * 
 * The System assumes ownership of the object passed to this method, and will 
 * delete it when the System is deleted.
 */
void DefaultSystemSubsystem::
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
void DefaultSystemSubsystem::
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
void DefaultSystemSubsystem::
addEventReporter(TriggeredEventReporter* handler) const {
    getGuts().updTriggeredEventReporters().push_back(handler);
}

} // namespace SimTK

