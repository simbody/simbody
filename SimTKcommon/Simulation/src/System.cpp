/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-13 Stanford University and the Authors.        *
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
    SimTK_ASSERT_ALWAYS(g, "System::adoptSystemGuts(): can't adopt null Guts");
    SimTK_ASSERT_ALWAYS(!guts,
        "System::adoptSystemGuts(): this System handle is already in use");
    guts = g;
    guts->setOwnerHandle(*this);
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

SubsystemIndex System::adoptSubsystem(Subsystem& child) {return updSystemGuts().adoptSubsystem(child);}
int System::getNumSubsystems() const {return getSystemGuts().getNumSubsystems();}
const Subsystem& System::getSubsystem(SubsystemIndex i) const {return getSystemGuts().getSubsystem(i);}
Subsystem& System::updSubsystem(SubsystemIndex i) {return updSystemGuts().updSubsystem(i);}
const DefaultSystemSubsystem& System::getDefaultSubsystem() const {
    return static_cast<const DefaultSystemSubsystem&>(getSystemGuts().getSubsystem(SubsystemIndex(0)));
}
DefaultSystemSubsystem& System::updDefaultSubsystem() {
    return static_cast<DefaultSystemSubsystem&>(updSystemGuts().updSubsystem(SubsystemIndex(0)));
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



const State& System::realizeTopology() const {return getSystemGuts().realizeTopology();}
void System::realizeModel(State& s) const {getSystemGuts().realizeModel(s);}
void System::realize(const State& s, Stage g) const {getSystemGuts().realize(s,g);}
void System::calcDecorativeGeometryAndAppend
   (const State& s, Stage g, Array_<DecorativeGeometry>& geom) const
{   getSystemGuts().calcDecorativeGeometryAndAppend(s,g,geom); }

void System::multiplyByN(const State& s, const Vector& u, Vector& dq) const
{   getSystemGuts().multiplyByN(s,u,dq); }
void System::multiplyByNTranspose(const State& s, const Vector& fq, Vector& fu) const
{   getSystemGuts().multiplyByNTranspose(s,fq,fu); }
void System::multiplyByNPInv(const State& s, const Vector& dq, Vector& u) const
{   getSystemGuts().multiplyByNPInv(s,dq,u); }
void System::multiplyByNPInvTranspose(const State& s, const Vector& fu, Vector& fq) const
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
void System::reportEvents(const State& s, Event::Cause cause, const Array_<EventId>& eventIds) const
{   getSystemGuts().reportEvents(s,cause,eventIds); }
void System::calcEventTriggerInfo(const State& s, Array_<EventTriggerInfo>& info) const
{   getSystemGuts().calcEventTriggerInfo(s,info); }
void System::calcTimeOfNextScheduledEvent(const State& s, Real& tNextEvent,
                                          Array_<EventId>& eventIds, bool includeCurrentTime) const
{   getSystemGuts().calcTimeOfNextScheduledEvent(s,tNextEvent,eventIds,includeCurrentTime); }
void System::calcTimeOfNextScheduledReport(const State& s, Real& tNextEvent,
                                          Array_<EventId>& eventIds, bool includeCurrentTime) const
{   getSystemGuts().calcTimeOfNextScheduledReport(s,tNextEvent,eventIds,includeCurrentTime); }




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
const Subsystem& System::Guts::getSubsystem(SubsystemIndex i) const {return getRep().getSubsystem(i);}
Subsystem& System::Guts::updSubsystem(SubsystemIndex i) {return updRep().updSubsystem(i);}

System::Guts* System::Guts::clone() const {
    return cloneImpl();
}



//------------------------------------------------------------------------------
//                            REALIZE TOPOLOGY
//------------------------------------------------------------------------------
const State& System::Guts::realizeTopology() const {
    State& defaultState = getRep().defaultState; // mutable
    if (getRep().systemTopologyHasBeenRealized())
        return defaultState;

    defaultState.clear();
    defaultState.setNumSubsystems(getNumSubsystems());
    for (SubsystemIndex i(0); i<getNumSubsystems(); ++i)
        defaultState.initializeSubsystem(i, getRep().subsystems[i].getName(),
                                            getRep().subsystems[i].getVersion());

    // Allow the concrete System subclass to do its processing.
    realizeTopologyImpl(defaultState); // defaultState is mutable

    // Realize any subsystems that the subclass didn't already take care of.
    for (SubsystemIndex i(0); i<getNumSubsystems(); ++i)
        if (getRep().subsystems[i].getStage(defaultState) < Stage::Topology)
            getRep().subsystems[i].getSubsystemGuts()
                                  .realizeSubsystemTopology(defaultState);

    // Force the defaultState's Topology stage version number to match the
    // Topology cache version in this System.
    defaultState.setSystemTopologyStageVersion
       (getRep().getSystemTopologyCacheVersion());
    defaultState.advanceSystemToStage(Stage::Topology);

    getRep().systemTopologyRealized = true; // mutable
    getRep().nRealizationsOfStage[Stage::Topology]++; // mutable counter

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
        // Allow the subclass to do its processing.
        realizeModelImpl(s);
        // Realize any subsystems that the subclass didn't already take care of.
        for (SubsystemIndex i(0); i<getNumSubsystems(); ++i)
            if (getRep().subsystems[i].getStage(s) < Stage::Model)
                getRep().subsystems[i].getSubsystemGuts()
                                      .realizeSubsystemModel(s);
        s.advanceSystemToStage(Stage::Model);

        getRep().nRealizationsOfStage[Stage::Model]++; // mutable counter
    }
}



//------------------------------------------------------------------------------
//                             REALIZE INSTANCE
//------------------------------------------------------------------------------
void System::Guts::realizeInstance(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Instance).prev(),
        "System::Guts::realizeInstance()");
    if (s.getSystemStage() < Stage::Instance) {
        realizeInstanceImpl(s);    // take care of the Subsystems
        // Realize any subsystems that the subclass didn't already take care of.
        for (SubsystemIndex i(0); i<getNumSubsystems(); ++i)
            if (getRep().subsystems[i].getStage(s) < Stage::Instance)
                getRep().subsystems[i].getSubsystemGuts()
                                      .realizeSubsystemInstance(s);
        s.advanceSystemToStage(Stage::Instance);

        getRep().nRealizationsOfStage[Stage::Instance]++; // mutable counter
    }
}



//------------------------------------------------------------------------------
//                              REALIZE TIME
//------------------------------------------------------------------------------
void System::Guts::realizeTime(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Time).prev(),
        "System::Guts::realizeTime()");
    if (s.getSystemStage() < Stage::Time) {
        // Allow the subclass to do processing.
        realizeTimeImpl(s);
        // Realize any subsystems that the subclass didn't already take care of.
        for (SubsystemIndex i(0); i<getNumSubsystems(); ++i)
            if (getRep().subsystems[i].getStage(s) < Stage::Time)
                getRep().subsystems[i].getSubsystemGuts()
                                      .realizeSubsystemTime(s);
        s.advanceSystemToStage(Stage::Time);

        getRep().nRealizationsOfStage[Stage::Time]++; // mutable counter
    }
}



//------------------------------------------------------------------------------
//                            REALIZE POSITION
//------------------------------------------------------------------------------
void System::Guts::realizePosition(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Position).prev(),
        "System::Guts::realizePosition()");
    if (s.getSystemStage() < Stage::Position) {
        // Allow the subclass to do processing.
        realizePositionImpl(s);
        // Realize any subsystems that the subclass didn't already take care of.
        for (SubsystemIndex i(0); i<getNumSubsystems(); ++i)
            if (getRep().subsystems[i].getStage(s) < Stage::Position)
                getRep().subsystems[i].getSubsystemGuts()
                                      .realizeSubsystemPosition(s);
        s.advanceSystemToStage(Stage::Position);

        getRep().nRealizationsOfStage[Stage::Position]++; // mutable counter
    }
}



//------------------------------------------------------------------------------
//                            REALIZE VELOCITY
//------------------------------------------------------------------------------
void System::Guts::realizeVelocity(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Velocity).prev(),
        "System::Guts::realizeVelocity()");
    if (s.getSystemStage() < Stage::Velocity) {
        // Allow the subclass to do processing.
        realizeVelocityImpl(s);
        // Realize any subsystems that the subclass didn't already take care of.
        for (SubsystemIndex i(0); i<getNumSubsystems(); ++i)
            if (getRep().subsystems[i].getStage(s) < Stage::Velocity)
                getRep().subsystems[i].getSubsystemGuts()
                                      .realizeSubsystemVelocity(s);
        s.advanceSystemToStage(Stage::Velocity);

        getRep().nRealizationsOfStage[Stage::Velocity]++; // mutable counter
    }
}



//------------------------------------------------------------------------------
//                            REALIZE DYNAMICS
//------------------------------------------------------------------------------
void System::Guts::realizeDynamics(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Dynamics).prev(),
        "System::Guts::realizeDynamics()");
    if (s.getSystemStage() < Stage::Dynamics) {
        // Allow the subclass to do processing.
        realizeDynamicsImpl(s);
        // Realize any subsystems that the subclass didn't already take care of.
        for (SubsystemIndex i(0); i<getNumSubsystems(); ++i)
            if (getRep().subsystems[i].getStage(s) < Stage::Dynamics)
                getRep().subsystems[i].getSubsystemGuts()
                                      .realizeSubsystemDynamics(s);
        s.advanceSystemToStage(Stage::Dynamics);

        getRep().nRealizationsOfStage[Stage::Dynamics]++; // mutable counter
    }
}



//------------------------------------------------------------------------------
//                           REALIZE ACCELERATION
//------------------------------------------------------------------------------
void System::Guts::realizeAcceleration(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Acceleration).prev(),
        "System::Guts::realizeAcceleration()");
    if (s.getSystemStage() < Stage::Acceleration) {
        // Allow the subclass to do processing.
        realizeAccelerationImpl(s);
        // Realize any subsystems that the subclass didn't already take care of.
        for (SubsystemIndex i(0); i<getNumSubsystems(); ++i)
            if (getRep().subsystems[i].getStage(s) < Stage::Acceleration)
                getRep().subsystems[i].getSubsystemGuts()
                                      .realizeSubsystemAcceleration(s);
        s.advanceSystemToStage(Stage::Acceleration);

        getRep().nRealizationsOfStage[Stage::Acceleration]++; // mutable counter
    }
}



//------------------------------------------------------------------------------
//                              REALIZE REPORT
//------------------------------------------------------------------------------
void System::Guts::realizeReport(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Report).prev(),
        "System::Guts::realizeReport()");
    if (s.getSystemStage() < Stage::Report) {
        // Allow the subclass to do processing.
        realizeReportImpl(s);
        // Realize any subsystems that the subclass didn't already take care of.
        for (SubsystemIndex i(0); i<getNumSubsystems(); ++i)
            if (getRep().subsystems[i].getStage(s) < Stage::Report)
                getRep().subsystems[i].getSubsystemGuts().realizeSubsystemReport(s);
        s.advanceSystemToStage(Stage::Report);

        getRep().nRealizationsOfStage[Stage::Report]++; // mutable counter
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

    Array_<EventId> eventsForSubsystem;
    for (SubsystemIndex sx(0); sx < getNumSubsystems(); ++sx) {
        const Subsystem::Guts& sub = getRep().subsystems[sx].getSubsystemGuts();

        if (!eventIds.empty()) {
            getSystem().getDefaultSubsystem().findSubsystemEventIds
               (sub.getMySubsystemIndex(), state, eventIds, eventsForSubsystem);
            if (eventsForSubsystem.empty())
                continue; // no events for this Subsystem
        }

        //---------------------------------------------------------------------
        HandleEventsResults subsysResults;
        sub.handleEvents(state,cause,eventsForSubsystem,options,subsysResults);
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

    Array_<EventId> eventsForSubsystem;
    for (SubsystemIndex sx(0); sx < getNumSubsystems(); ++sx) {
        const Subsystem::Guts& sub = getRep().subsystems[sx].getSubsystemGuts();

        if (!eventIds.empty()) {
            getSystem().getDefaultSubsystem().findSubsystemEventIds
               (sub.getMySubsystemIndex(), s, eventIds, eventsForSubsystem);
            if (eventsForSubsystem.empty())
                continue; // no reporting events for this Subsystem
        }

        //---------------------------------------------------------------------
        sub.reportEvents(s, cause, eventsForSubsystem);
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
        for (Array_<EventTriggerInfo>::const_iterator e = subinfo.begin();
             e != subinfo.end(); ++e)
        {
            info.push_back(*e);
        }
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



//==============================================================================
//                         SYSTEM :: GUTS :: GUTS REP
//==============================================================================
// All inline currently.



//==============================================================================
//                         DEFAULT SYSTEM SUBSYSTEM
//==============================================================================

// Instance-stage cache entry that records the motion method for every q, u,
// and udot.
// TODO: not used yet
class MotionInfo {
public:
    // Allocated at end of realize(Model) when we know nq and nu and initialized
    // to "free". Subsystems fill these in during realize(Instance).
    Array_<unsigned char> qMotion;      // nq
    Array_<unsigned char> uMotion;      // nu
    Array_<unsigned char> udotMotion;   // nu

    // We'll calculate these index arrays once all Subsystems have had a
    // chance to set the motion info for their variables.
    Array_<int>& freeQIndex;    // nfq indices of each q still marked "free"
    Array_<int>& freeUIndex;    // nfu indices of each u still marked "free"
    Array_<int>& freeUDotIndex; // nfudot indices of udot still marked "free"
};

// Whether to use relative accuracy for u and z variables.
// TODO: not used yet
class StateVarUseRelAccuracyUZ {
public:
    // Whether to scale by 1/u and 1/z when u and z are large. We never
    // use relative accuracy for q.
    Array_<bool> useRelativeAccuracyU; // nu
    Array_<bool> useRelativeAccuracyZ; // nz
};

// This class stores various information used by the default subsystem
// that needs to be stored in the state.
class CachedEventInfo {
public:
    CachedEventInfo() : eventIdCounter(0) {}
    int eventIdCounter;
    std::map<int, SubsystemIndex> eventOwnerMap;
    Array_<EventId> scheduledEventIds;
    Array_<EventTriggerByStageIndex> triggeredEventIndices;
    Array_<EventId> triggeredEventIds;
    Array_<EventId> scheduledReportIds;
    Array_<EventTriggerByStageIndex> triggeredReportIndices;
    Array_<EventId> triggeredReportIds;
};

class DefaultSystemSubsystem::Guts : public Subsystem::Guts {
public:

    Guts() : Subsystem::Guts("DefaultSystemSubsystem::Guts", "0.0.1") { }

    ~Guts() {
        for (int i = 0; i < (int)scheduledEventHandlers.size(); ++i)
            delete scheduledEventHandlers[i];
        for (int i = 0; i < (int)triggeredEventHandlers.size(); ++i)
            delete triggeredEventHandlers[i];
        for (int i = 0; i < (int)scheduledEventReporters.size(); ++i)
            delete scheduledEventReporters[i];
        for (int i = 0; i < (int)triggeredEventReporters.size(); ++i)
            delete triggeredEventReporters[i];
    }

    Guts* cloneImpl() const override {
        return new Guts(*this);
    }

    const Array_<ScheduledEventHandler*>& getScheduledEventHandlers() const {
        return scheduledEventHandlers;
    }

    Array_<ScheduledEventHandler*>& updScheduledEventHandlers() {
        invalidateSubsystemTopologyCache();
        return scheduledEventHandlers;
    }

    const Array_<TriggeredEventHandler*>& getTriggeredEventHandlers() const {
        return triggeredEventHandlers;
    }

    Array_<TriggeredEventHandler*>& updTriggeredEventHandlers() {
        invalidateSubsystemTopologyCache();
        return triggeredEventHandlers;
    }

    const Array_<ScheduledEventReporter*>& getScheduledEventReporters() const {
        return scheduledEventReporters;
    }

    Array_<ScheduledEventReporter*>& updScheduledEventReporters() const {
        invalidateSubsystemTopologyCache();
        return scheduledEventReporters;
    }

    const Array_<TriggeredEventReporter*>& getTriggeredEventReporters() const {
        return triggeredEventReporters;
    }

    Array_<TriggeredEventReporter*>& updTriggeredEventReporters() const {
        invalidateSubsystemTopologyCache();
        return triggeredEventReporters;
    }

    const CachedEventInfo& getCachedEventInfo(const State& s) const {
        return Value<CachedEventInfo>::downcast
           (getCacheEntry(s, cachedEventInfoIndex)).get();
    }

    CachedEventInfo& updCachedEventInfo(const State& s) const {
        return Value<CachedEventInfo>::downcast
           (updCacheEntry(s, cachedEventInfoIndex)).upd();
    }

    int realizeSubsystemTopologyImpl(State& s) const override {
        cachedEventInfoIndex = s.allocateCacheEntry(getMySubsystemIndex(),
                                                    Stage::Topology,
                                                    new Value<CachedEventInfo>());
        CachedEventInfo& info = updCachedEventInfo(s);
        info.scheduledEventIds.clear();
        info.triggeredEventIndices.clear();
        info.triggeredEventIds.clear();
        info.scheduledReportIds.clear();
        info.triggeredReportIndices.clear();
        info.triggeredReportIds.clear();
        info.eventIdCounter = 0;
        for (Array_<ScheduledEventHandler*>::const_iterator
                 e = scheduledEventHandlers.begin();
                 e != scheduledEventHandlers.end(); ++e) {
            EventId id;
            createScheduledEvent(s, id);
            info.scheduledEventIds.push_back(id);
        }
        for (Array_<TriggeredEventHandler*>::const_iterator
                 e = triggeredEventHandlers.begin();
                 e != triggeredEventHandlers.end(); ++e) {
            EventId id;
            EventTriggerByStageIndex index;
            createTriggeredEvent(s, id, index, (*e)->getRequiredStage());
            info.triggeredEventIds.push_back(id);
            info.triggeredEventIndices.push_back(index);
        }
        for (Array_<ScheduledEventReporter*>::const_iterator
                 e = scheduledEventReporters.begin();
                 e != scheduledEventReporters.end(); ++e) {
            EventId id;
            createScheduledEvent(s, id);
            info.scheduledReportIds.push_back(id);
        }
        for (Array_<TriggeredEventReporter*>::const_iterator
                 e = triggeredEventReporters.begin();
                 e != triggeredEventReporters.end(); ++e) {
            EventId id;
            EventTriggerByStageIndex index;
            createTriggeredEvent(s, id, index, (*e)->getRequiredStage());
            info.triggeredReportIds.push_back(id);
            info.triggeredReportIndices.push_back(index);
        }
        return 0;
    }

    int realizeSubsystemModelImpl(State& s) const override {
        return 0;
    }

    int realizeEventTriggers(const State& s, Stage g) const {
        const CachedEventInfo& info = getCachedEventInfo(s);
        Vector& triggers = s.updEventTriggersByStage(getMySubsystemIndex(), g);
        for (int i = 0; i < (int)triggeredEventHandlers.size(); ++i) {
            if (g == triggeredEventHandlers[i]->getRequiredStage())
                triggers[info.triggeredEventIndices[i]] =
                    triggeredEventHandlers[i]->getValue(s);
        }
        for (int i = 0; i < (int)triggeredEventReporters.size(); ++i) {
            if (g == triggeredEventReporters[i]->getRequiredStage())
                triggers[info.triggeredReportIndices[i]] =
                    triggeredEventReporters[i]->getValue(s);
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

        const CachedEventInfo& info = getCachedEventInfo(s);
        trigInfo.resize(  triggeredEventHandlers.size()
                        + triggeredEventReporters.size());

        for (int i = 0; i < (int)triggeredEventHandlers.size(); ++i) {
            const Stage stage = triggeredEventHandlers[i]->getRequiredStage();
            const SystemEventTriggerIndex index
               (  info.triggeredEventIndices[i]
                + s.getEventTriggerStartByStage(stage)
                + s.getEventTriggerStartByStage(getMySubsystemIndex(), stage));
            trigInfo[index] = triggeredEventHandlers[i]->getTriggerInfo();
            trigInfo[index].setEventId(info.triggeredEventIds[i]);
        }

        for (int i = 0; i < (int)triggeredEventReporters.size(); ++i) {
            const Stage stage = triggeredEventReporters[i]->getRequiredStage();
            const SystemEventTriggerIndex index
               (  info.triggeredReportIndices[i]
                + s.getEventTriggerStartByStage(stage)
                + s.getEventTriggerStartByStage(getMySubsystemIndex(), stage));
            trigInfo[index] = triggeredEventReporters[i]->getTriggerInfo();
            trigInfo[index].setEventId(info.triggeredReportIds[i]);
        }
    }
    void calcTimeOfNextScheduledEventImpl(const State& s, Real& tNextEvent,
        Array_<EventId>& eventIds, bool includeCurrentTime) const override {
        // Loop over all registered ScheduledEventHandlers, and ask each one
        // when its next event occurs.

        const CachedEventInfo& info = getCachedEventInfo(s);
        tNextEvent = Infinity;
        for (int i = 0; i < (int)scheduledEventHandlers.size(); ++i) {
            Real time = scheduledEventHandlers[i]->getNextEventTime
                                                        (s, includeCurrentTime);
            if (time <= tNextEvent
                && (time > s.getTime()
                    || (includeCurrentTime && time == s.getTime())))
            {
                if (time < tNextEvent)
                    eventIds.clear();
                tNextEvent = time;
                eventIds.push_back(info.scheduledEventIds[i]);
            }
        }
    }
    void calcTimeOfNextScheduledReportImpl(const State& s, Real& tNextEvent,
        Array_<EventId>& eventIds, bool includeCurrentTime) const override {
        // Loop over all registered ScheduledEventReporters, and ask each one
        // when its next event occurs.

        const CachedEventInfo& info = getCachedEventInfo(s);
        tNextEvent = Infinity;
        for (int i = 0; i < (int)scheduledEventReporters.size(); ++i) {
            Real time = scheduledEventReporters[i]->getNextEventTime
                                                        (s, includeCurrentTime);
            if (time <= tNextEvent
                && (time > s.getTime()
                    || (includeCurrentTime && time == s.getTime())))
            {
                if (time < tNextEvent)
                    eventIds.clear();
                tNextEvent = time;
                eventIds.push_back(info.scheduledReportIds[i]);
            }
        }
    }
    void handleEventsImpl(State& s, Event::Cause cause,
                      const Array_<EventId>& eventIds,
                      const HandleEventsOptions& options,
                      HandleEventsResults& results) const override
    {
        const CachedEventInfo& info = getCachedEventInfo(s);
        const Real accuracy = options.getAccuracy();
        bool shouldTerminate = false;

        // Build a set of the ids for quick lookup.
        std::set<EventId> idSet;
        for (int i = 0; i < (int)eventIds.size(); ++i)
            idSet.insert(eventIds[i]);

        // Process triggered events.

        if (cause == Event::Cause::Triggered) {
            for (int i = 0; i < (int)triggeredEventHandlers.size(); ++i) {
                if (idSet.find(info.triggeredEventIds[i]) != idSet.end()) {
                    bool eventShouldTerminate = false;
                    triggeredEventHandlers[i]->handleEvent
                                            (s, accuracy, eventShouldTerminate);
                    if (eventShouldTerminate)
                        shouldTerminate = true;
                }
            }
            for (int i = 0; i < (int)triggeredEventReporters.size(); ++i) {
                if (idSet.find(info.triggeredReportIds[i]) != idSet.end())
                    triggeredEventReporters[i]->handleEvent(s);
            }
        }

        // Process scheduled events.

        if (cause == Event::Cause::Scheduled) {
            for (int i = 0; i < (int)scheduledEventHandlers.size(); ++i) {
                if (idSet.find(info.scheduledEventIds[i]) != idSet.end()) {
                    bool eventShouldTerminate = false;
                    scheduledEventHandlers[i]->handleEvent
                                            (s, accuracy, eventShouldTerminate);
                    if (eventShouldTerminate)
                        shouldTerminate = true;
                }
            }
            for (int i = 0; i < (int)scheduledEventReporters.size(); ++i) {
                if (idSet.find(info.scheduledReportIds[i]) != idSet.end())
                    scheduledEventReporters[i]->handleEvent(s);
            }
        }

        // Assume some change was made.
        results.setAnyChangeMade(true);

        results.setExitStatus(shouldTerminate
            ? HandleEventsResults::ShouldTerminate
            : HandleEventsResults::Succeeded);
    }

    void reportEventsImpl(const State& s, Event::Cause cause,
                      const Array_<EventId>& eventIds) const override {
        const CachedEventInfo& info = getCachedEventInfo(s);

        // Build a set of the ids for quick lookup.

        std::set<EventId> idSet;
        for (int i = 0; i < (int)eventIds.size(); ++i)
            idSet.insert(eventIds[i]);

        // Process triggered events.

        if (cause == Event::Cause::Triggered) {
            for (int i = 0; i < (int)triggeredEventReporters.size(); ++i) {
                if (idSet.find(info.triggeredReportIds[i]) != idSet.end())
                    triggeredEventReporters[i]->handleEvent(s);
            }
        }

        // Process scheduled events.

        if (cause == Event::Cause::Scheduled) {
            for (int i = 0; i < (int)scheduledEventReporters.size(); ++i) {
                if (idSet.find(info.scheduledReportIds[i]) != idSet.end())
                    scheduledEventReporters[i]->handleEvent(s);
            }
        }
    }

private:
    mutable CacheEntryIndex                 cachedEventInfoIndex;
    mutable Array_<ScheduledEventHandler*>  scheduledEventHandlers;
    mutable Array_<TriggeredEventHandler*>  triggeredEventHandlers;
    mutable Array_<ScheduledEventReporter*> scheduledEventReporters;
    mutable Array_<TriggeredEventReporter*> triggeredEventReporters;
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

/**
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

/**
 * Generate a new, globally unique event ID.  Typically you will not call this
 * directly.  When a Subsystem needs to obtain an event ID for an event it
 * defines, it should do so by calling Subsystem::Guts::createScheduledEvent()
 * or Subsystem::Guts::createTriggeredEvent().
 */
EventId DefaultSystemSubsystem::
createEventId(SubsystemIndex subsys, const State& state) const {
    // Must use "upd" here because this is called from realize()
    // while we're still filling in the CachedEventInfo.
    CachedEventInfo& info = getGuts().updCachedEventInfo(state);
    int id = info.eventIdCounter++;
    info.eventOwnerMap[id] = subsys;
    return EventId(id);
}

/**
 * Given a list of event IDs, filter it to produce a list of those events
 * belonging to a particular Subsystem.
 *
 * @param subsys       the Subsystem for which to find events
 * @param state        the State which produced the events
 * @param allEvents    a list of event IDs to filter
 * @param eventsForSubsystem  on exit, this Array_ contains the filtered list
 *                            of event IDs belonging to the specified Subsystem.
 */
void DefaultSystemSubsystem::findSubsystemEventIds
   (SubsystemIndex subsys, const State& state, const Array_<EventId>& allEvents,
    Array_<EventId>& eventsForSubsystem) const
{
    const CachedEventInfo& info = getGuts().getCachedEventInfo(state);
    eventsForSubsystem.clear();
    for (int i = 0; i < (int)allEvents.size(); ++i) {
        std::map<int, SubsystemIndex>::const_iterator p =
            info.eventOwnerMap.find(allEvents[i]);
        assert(p != info.eventOwnerMap.end());
        if (p->second == subsys)
            eventsForSubsystem.push_back(allEvents[i]);
    }
}


} // namespace SimTK

