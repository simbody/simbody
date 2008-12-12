/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-7 Stanford University and the Authors.         *
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

#include "SystemGutsRep.h"

#include <cassert>

namespace SimTK {

    ////////////
    // SYSTEM //
    ////////////

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

// Don't use ordinary delete, assignment, or copy here. Must go
// through the library-side VFT to get access to the correct client-side
// virtual destructor and clone method.
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

void System::resetAllCountersToZero() {updSystemGuts().updRep().resetAllCounters();}
long System::getNumRealizationsOfThisStage(Stage g) const {return getSystemGuts().getRep().nRealizationsOfStage[g];}
long System::getNumRealizeCalls() const {return getSystemGuts().getRep().nRealizeCalls;}
long System::getNumQProjections() const {return getSystemGuts().getRep().nQProjections;}
long System::getNumUProjections() const {return getSystemGuts().getRep().nUProjections;} 
long System::getNumQErrorEstimateProjections() const {return getSystemGuts().getRep().nQErrEstProjections;}
long System::getNumUErrorEstimateProjections() const {return getSystemGuts().getRep().nUErrEstProjections;}
long System::getNumProjectCalls() const {return getSystemGuts().getRep().nProjectCalls;}
long System::getNumHandlerCallsThatChangedStage(Stage g) const {return getSystemGuts().getRep().nHandlerCallsThatChangedStage[g];}
long System::getNumHandleEventCalls() const {return getSystemGuts().getRep().nHandleEventsCalls;}
long System::getNumReportEventCalls() const {return getSystemGuts().getRep().nReportEventsCalls;}

const State& System::getDefaultState() const {return getSystemGuts().getDefaultState();}
State& System::updDefaultState() {return updSystemGuts().updDefaultState();}

SubsystemIndex System::adoptSubsystem(Subsystem& child) {return updSystemGuts().adoptSubsystem(child);}
int System::getNSubsystems() const {return getSystemGuts().getNSubsystems();}
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
bool System::systemTopologyHasBeenRealized() const {
    return getSystemGuts().systemTopologyHasBeenRealized();
}
const State& System::realizeTopology() const {return getSystemGuts().realizeTopology();}
void System::realizeModel(State& s) const {getSystemGuts().realizeModel(s);}
void System::realize(const State& s, Stage g) const {getSystemGuts().realize(s,g);}
void System::calcDecorativeGeometryAndAppend(const State& s, Stage g, std::vector<DecorativeGeometry>& geom) const {
    getSystemGuts().calcDecorativeGeometryAndAppend(s,g,geom);
}

Real System::calcTimescale(const State& s) const {return getSystemGuts().calcTimescale(s);}
void System::calcYUnitWeights(const State& s, Vector& weights) const
  { getSystemGuts().calcYUnitWeights(s,weights); }
void System::project(State& s, Real consAccuracy, const Vector& yweights,
                     const Vector& ootols, Vector& yerrest, System::ProjectOptions opts) const
  { getSystemGuts().project(s,consAccuracy,yweights,ootols,yerrest, opts); }
void System::calcYErrUnitTolerances(const State& s, Vector& tolerances) const
  { getSystemGuts().calcYErrUnitTolerances(s,tolerances); }
void System::handleEvents(State& s, EventCause cause, const std::vector<EventId>& eventIds,
                          Real accuracy, const Vector& yWeights, const Vector& ooConstraintTols,
                          Stage& lowestModified, bool& shouldTerminate) const
  { getSystemGuts().handleEvents(s,cause,eventIds,accuracy,yWeights,ooConstraintTols,
                               lowestModified,shouldTerminate); }
void System::reportEvents(const State& s, EventCause cause, const std::vector<EventId>& eventIds) const
  { getSystemGuts().reportEvents(s,cause,eventIds); }
void System::calcEventTriggerInfo(const State& s, std::vector<EventTriggerInfo>& info) const
  { getSystemGuts().calcEventTriggerInfo(s,info); }
void System::calcTimeOfNextScheduledEvent(const State& s, Real& tNextEvent,
                                          std::vector<EventId>& eventIds, bool includeCurrentTime) const
  { getSystemGuts().calcTimeOfNextScheduledEvent(s,tNextEvent,eventIds,includeCurrentTime); }
void System::calcTimeOfNextScheduledReport(const State& s, Real& tNextEvent,
                                          std::vector<EventId>& eventIds, bool includeCurrentTime) const
  { getSystemGuts().calcTimeOfNextScheduledReport(s,tNextEvent,eventIds,includeCurrentTime); }

const char* System::getEventCauseName(System::EventCause cause) {
    switch(cause) {
    case TriggeredEvents:   return "TriggeredEvents";
    case ScheduledEvents:   return "ScheduledEvents";
    case TimeAdvancedEvent: return "TimeAdvancedEvent";
    case TerminationEvent:  return "TerminationEvent";
    case InvalidEventCause: return "InvalidEventCause";
    }
    return "UNRECOGNIZED EVENT CAUSE";
}


    //////////////////
    // SYSTEM::GUTS //
    //////////////////

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

bool System::Guts::systemTopologyHasBeenRealized() const {
    return getRep().systemTopologyHasBeenRealized();
}

void System::Guts::invalidateSystemTopologyCache() const {
    return getRep().invalidateSystemTopologyCache(); // mutable
}

const State& System::Guts::getDefaultState() const {
    return getRep().getDefaultState();
}

State& System::Guts::updDefaultState() {
    return updRep().updDefaultState();
}

int System::Guts::getNSubsystems() const {return getRep().getNSubsystems();}
const Subsystem& System::Guts::getSubsystem(SubsystemIndex i) const {return getRep().getSubsystem(i);}
Subsystem& System::Guts::updSubsystem(SubsystemIndex i) {return updRep().updSubsystem(i);}

System::Guts* System::Guts::clone() const {
    return cloneImpl();
}

const State& System::Guts::realizeTopology() const {
    State& defaultState = getRep().defaultState; // mutable
    if (!getRep().systemTopologyHasBeenRealized()) {
        defaultState.clear();
        defaultState.setNSubsystems(getNSubsystems());
        for (SubsystemIndex i(0); i<getNSubsystems(); ++i) 
            defaultState.initializeSubsystem(i, getRep().subsystems[i].getName(), 
                                                getRep().subsystems[i].getVersion());
        
        // Allow the subclass to do processing.
        realizeTopologyImpl(defaultState); // defaultState is mutable
        // Realize any subsystems that the subclass didn't already take care of.
        for (SubsystemIndex i(0); i<getNSubsystems(); ++i)
            if (getRep().subsystems[i].getStage(defaultState) < Stage::Topology)
                getRep().subsystems[i].getSubsystemGuts().realizeSubsystemTopology(defaultState);
        getRep().systemTopologyRealized = true; // mutable
        defaultState.advanceSystemToStage(Stage::Topology);

        // Realize the model using the default settings of the Model variables.
        // This allocates all the later-stage State variables.
        realizeModel(defaultState);

        // Now realize the default state to the highest Stage.
        // TODO: this is problematic if the default state is not a valid one.
        //realize(defaultState, Stage::HighestValid);

        getRep().nRealizationsOfStage[Stage::Topology]++; // mutable counter
    }
    return defaultState;
}

void System::Guts::realizeModel(State& s) const {
    SimTK_STAGECHECK_TOPOLOGY_REALIZED_ALWAYS(systemTopologyHasBeenRealized(),
        "System", getName(), "System::Guts::realizeModel()");
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Topology, 
        "System::Guts::realizeModel()");
    if (s.getSystemStage() < Stage::Model) {
        // Allow the subclass to do processing.
        realizeModelImpl(s);
        // Realize any subsystems that the subclass didn't already take care of.
        for (SubsystemIndex i(0); i<getNSubsystems(); ++i)
            if (getRep().subsystems[i].getStage(s) < Stage::Model)
                getRep().subsystems[i].getSubsystemGuts().realizeSubsystemModel(s);
        s.advanceSystemToStage(Stage::Model);

        getRep().nRealizationsOfStage[Stage::Model]++; // mutable counter
    }
}
void System::Guts::realizeInstance(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Instance).prev(), 
        "System::Guts::realizeInstance()");
    if (s.getSystemStage() < Stage::Instance) {
        realizeInstanceImpl(s);    // take care of the Subsystems
        // Realize any subsystems that the subclass didn't already take care of.
        for (SubsystemIndex i(0); i<getNSubsystems(); ++i)
            if (getRep().subsystems[i].getStage(s) < Stage::Instance)
                getRep().subsystems[i].getSubsystemGuts().realizeSubsystemInstance(s);
        s.advanceSystemToStage(Stage::Instance);

        getRep().nRealizationsOfStage[Stage::Instance]++; // mutable counter
    }
}
void System::Guts::realizeTime(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Time).prev(), 
        "System::Guts::realizeTime()");
    if (s.getSystemStage() < Stage::Time) {
        // Allow the subclass to do processing.
        realizeTimeImpl(s);
        // Realize any subsystems that the subclass didn't already take care of.
        for (SubsystemIndex i(0); i<getNSubsystems(); ++i)
            if (getRep().subsystems[i].getStage(s) < Stage::Time)
                getRep().subsystems[i].getSubsystemGuts().realizeSubsystemTime(s);
        s.advanceSystemToStage(Stage::Time);

        getRep().nRealizationsOfStage[Stage::Time]++; // mutable counter
    }
}
void System::Guts::realizePosition(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Position).prev(), 
        "System::Guts::realizePosition()");
    if (s.getSystemStage() < Stage::Position) {
        // Allow the subclass to do processing.
        realizePositionImpl(s);
        // Realize any subsystems that the subclass didn't already take care of.
        for (SubsystemIndex i(0); i<getNSubsystems(); ++i)
            if (getRep().subsystems[i].getStage(s) < Stage::Position)
                getRep().subsystems[i].getSubsystemGuts().realizeSubsystemPosition(s);
        s.advanceSystemToStage(Stage::Position);

        getRep().nRealizationsOfStage[Stage::Position]++; // mutable counter
    }
}
void System::Guts::realizeVelocity(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Velocity).prev(), 
        "System::Guts::realizeVelocity()");
    if (s.getSystemStage() < Stage::Velocity) {
        // Allow the subclass to do processing.
        realizeVelocityImpl(s);
        // Realize any subsystems that the subclass didn't already take care of.
        for (SubsystemIndex i(0); i<getNSubsystems(); ++i)
            if (getRep().subsystems[i].getStage(s) < Stage::Velocity)
                getRep().subsystems[i].getSubsystemGuts().realizeSubsystemVelocity(s);
        s.advanceSystemToStage(Stage::Velocity);

        getRep().nRealizationsOfStage[Stage::Velocity]++; // mutable counter
    }
}
void System::Guts::realizeDynamics(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Dynamics).prev(), 
        "System::Guts::realizeDynamics()");
    if (s.getSystemStage() < Stage::Dynamics) {
        // Allow the subclass to do processing.
        realizeDynamicsImpl(s);
        // Realize any subsystems that the subclass didn't already take care of.
        for (SubsystemIndex i(0); i<getNSubsystems(); ++i)
            if (getRep().subsystems[i].getStage(s) < Stage::Dynamics)
                getRep().subsystems[i].getSubsystemGuts().realizeSubsystemDynamics(s);
        s.advanceSystemToStage(Stage::Dynamics);

        getRep().nRealizationsOfStage[Stage::Dynamics]++; // mutable counter
    }
}
void System::Guts::realizeAcceleration(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Acceleration).prev(), 
        "System::Guts::realizeAcceleration()");
    if (s.getSystemStage() < Stage::Acceleration) {
        // Allow the subclass to do processing.
        realizeAccelerationImpl(s);
        // Realize any subsystems that the subclass didn't already take care of.
        for (SubsystemIndex i(0); i<getNSubsystems(); ++i)
            if (getRep().subsystems[i].getStage(s) < Stage::Acceleration)
                getRep().subsystems[i].getSubsystemGuts().realizeSubsystemAcceleration(s);
        s.advanceSystemToStage(Stage::Acceleration);

        getRep().nRealizationsOfStage[Stage::Acceleration]++; // mutable counter
    }
}
void System::Guts::realizeReport(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Report).prev(), 
        "System::Guts::realizeReport()");
    if (s.getSystemStage() < Stage::Report) {
        // Allow the subclass to do processing.
        realizeReportImpl(s);
        // Realize any subsystems that the subclass didn't already take care of.
        for (SubsystemIndex i(0); i<getNSubsystems(); ++i)
            if (getRep().subsystems[i].getStage(s) < Stage::Report)
                getRep().subsystems[i].getSubsystemGuts().realizeSubsystemReport(s);
        s.advanceSystemToStage(Stage::Report);

        getRep().nRealizationsOfStage[Stage::Report]++; // mutable counter
    }
}

Real System::Guts::calcTimescale(const State& s) const {
    SimTK_STAGECHECK_GE(s.getSystemStage(), Stage::Instance,
        "System::Guts::calcTimescale()");
    return calcTimescaleImpl(s);
}

void System::Guts::calcYUnitWeights(const State& s, Vector& weights) const {
    SimTK_STAGECHECK_GE(s.getSystemStage(), Stage::Position,
        "System::Guts::calcYUnitWeights()");
    calcYUnitWeightsImpl(s,weights);
}

void System::Guts::calcYErrUnitTolerances(const State& s, Vector& tolerances) const {
    SimTK_STAGECHECK_GE(s.getSystemStage(), Stage::Instance,
        "System::Guts::calcYErrUnitTolerances()");
    calcYErrUnitTolerancesImpl(s,tolerances);
}

void System::Guts::project(State& s, Real consAccuracy, const Vector& yweights,
                           const Vector& ootols, Vector& yerrest, System::ProjectOptions opts) const
{
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Instance, // TODO: is this the right stage?
        "System::Guts::project()");
    projectImpl(s,consAccuracy,yweights,ootols,yerrest,opts);

    getRep().nProjectCalls++; // mutable counter
    if (opts&ProjectOptions::Q) getRep().nQProjections++;
    if (opts&ProjectOptions::U) getRep().nUProjections++;
    if (opts&ProjectOptions::QError) getRep().nQErrEstProjections++;
    if (opts&ProjectOptions::UError) getRep().nUErrEstProjections++;
}

void System::Guts::handleEvents
   (State& s, EventCause cause, const std::vector<EventId>& eventIds,
    Real accuracy, const Vector& yWeights, const Vector& ooConstraintTols,
    Stage& lowestModified, bool& shouldTerminate) const
{
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Model, // TODO: is this the right stage?
        "System::Guts::handleEvents()");
    const Real savedTime = s.getTime();
    handleEventsImpl(s,cause,eventIds,accuracy,yWeights,ooConstraintTols,
                     lowestModified, shouldTerminate);

    getRep().nHandleEventsCalls++; // mutable counters
    getRep().nHandlerCallsThatChangedStage[lowestModified]++;

    SimTK_ASSERT_ALWAYS(s.getTime() == savedTime,
        "System::Guts::handleEvents(): handleEventsImpl() tried to change the time");
}

void System::Guts::reportEvents
   (const State& s, EventCause cause, const std::vector<EventId>& eventIds) const
{
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Model, // TODO: is this the right stage?
        "System::Guts::reportEvents()");
    const Real savedTime = s.getTime();
    reportEventsImpl(s,cause,eventIds);

    getRep().nReportEventsCalls++; // mutable counter
    getRep().nHandlerCallsThatChangedStage[Stage::Report]++;
}

void System::Guts::calcTimeOfNextScheduledEvent
    (const State& s, Real& tNextEvent, std::vector<EventId>& eventIds, bool includeCurrentTime) const
{
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Time,
        "System::Guts::calcTimeOfNextScheduledEvent()");
    tNextEvent = Infinity;
    calcTimeOfNextScheduledEventImpl(s,tNextEvent,eventIds,includeCurrentTime);
}

void System::Guts::calcTimeOfNextScheduledReport
    (const State& s, Real& tNextEvent, std::vector<EventId>& eventIds, bool includeCurrentTime) const
{
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Time,
        "System::Guts::calcTimeOfNextScheduledReport()");
    tNextEvent = Infinity;
    calcTimeOfNextScheduledReportImpl(s,tNextEvent,eventIds,includeCurrentTime);
}

void System::Guts::calcEventTriggerInfo(const State& s, std::vector<EventTriggerInfo>& info) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Instance,
        "System::Guts::calcEventTriggerInfo()");
    calcEventTriggerInfoImpl(s,info);
}

void System::Guts::realize(const State& s, Stage g) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Model, 
        "System::Guts::realize()");

    Stage stageNow = Stage::Empty;
    while ((stageNow=s.getSystemStage()) < g) {
        switch (stageNow) {
        case Stage::ModelIndex:        realizeInstance(s);     break;
        case Stage::InstanceIndex:     realizeTime(s);         break;
        case Stage::TimeIndex:         realizePosition(s);     break;
        case Stage::PositionIndex:     realizeVelocity(s);     break;
        case Stage::VelocityIndex:     realizeDynamics(s);     break;
        case Stage::DynamicsIndex:     realizeAcceleration(s); break;
        case Stage::AccelerationIndex: realizeReport(s);       break;
        default: assert(!"System::Guts::realize(): bad stage");
        }
    }
}

void System::Guts::calcDecorativeGeometryAndAppend(const State& s, Stage stage, std::vector<DecorativeGeometry>& geom) const {
    assert(stage==Stage::Topology || s.getSystemStage() >= stage);
    for (SubsystemIndex i(0); i<getNSubsystems(); ++i)
        getRep().subsystems[i].getSubsystemGuts().calcDecorativeGeometryAndAppend(s, stage, geom);
}


SubsystemIndex System::Guts::adoptSubsystem(Subsystem& src) {
    return updRep().adoptSubsystem(src);
}


// These are the default implementations for the System virtual functions.
// Note that this DOES NOT cause binary compatibility problems. The addresses of
// these functions will be supplied from the library side, but these addresses will
// get filled in to the default virtual function table on the *client* side which
// knows where to put each function by name.

int System::Guts::realizeTopologyImpl(State& s) const { 
    return 0;
}
int System::Guts::realizeModelImpl(State& s) const {
    return 0;
}
int System::Guts::realizeInstanceImpl(const State& s) const { 
    return 0;
}
int System::Guts::realizeTimeImpl(const State& s) const { 
    return 0;
}
int System::Guts::realizePositionImpl(const State& s) const { 
    return 0;
}
int System::Guts::realizeVelocityImpl(const State& s) const { 
    return 0;
}
int System::Guts::realizeDynamicsImpl(const State& s) const { 
    return 0;
}
int System::Guts::realizeAccelerationImpl(const State& s) const { 
    return 0;
}
int System::Guts::realizeReportImpl(const State& s) const { 
    return 0;
}

Real System::Guts::calcTimescaleImpl(const State&) const {
    return 0.1; // TODO!!!
}

int System::Guts::calcYUnitWeightsImpl(const State& s, Vector& weights) const {
    weights.resize(s.getNY());
    VectorView qwts = weights(s.getQStart(), s.getNQ());   // writable views
    VectorView uwts = weights(s.getUStart(), s.getNU());
    VectorView zwts = weights(s.getZStart(), s.getNZ());

    for (SubsystemIndex i(0); i<getNSubsystems(); ++i) {
        const Subsystem::Guts& sub = getRep().subsystems[i].getSubsystemGuts();
        sub.calcQUnitWeights(s, qwts(s.getQStart(i), s.getNQ(i)));
        sub.calcUUnitWeights(s, uwts(s.getUStart(i), s.getNU(i)));
        sub.calcZUnitWeights(s, zwts(s.getZStart(i), s.getNZ(i)));
    }
    return 0;
}

int System::Guts::projectImpl(State&, Real consAccuracy, const Vector& yweights,
                              const Vector& ootols, Vector& yerrest, System::ProjectOptions) const
{
    SimTK_THROW2(Exception::UnimplementedVirtualMethod, "System", "project"); 
    return std::numeric_limits<int>::min();
}


int System::Guts::calcYErrUnitTolerancesImpl(const State& s, Vector& ootols) const {
    ootols.resize(s.getNYErr());
    VectorView qtols = ootols(s.getQErrStart(), s.getNQErr()); // writable views
    VectorView utols = ootols(s.getUErrStart(), s.getNUErr());

    for (SubsystemIndex i(0); i<getNSubsystems(); ++i) {
        const Subsystem::Guts& sub = getRep().subsystems[i].getSubsystemGuts();
        sub.calcQErrUnitTolerances(s, qtols(s.getQErrStart(i), s.getNQErr(i)));
        sub.calcUErrUnitTolerances(s, utols(s.getUErrStart(i), s.getNUErr(i)));
    }
    return 0;
}

int System::Guts::handleEventsImpl
   (State& s, EventCause cause, const std::vector<EventId>& eventIds,
    Real accuracy, const Vector& yWeights, const Vector& ooConstraintTols,
    Stage& lowestModified, bool& shouldTerminate) const
{
    // Event handlers should not be able to modify the time.
    
    std::set<SubsystemIndex> temp;
    State restricted;
    s.createRestrictedState(restricted, Stage::Time, temp);
    
    // Loop over each subsystem, see which events belong to it, and allow it to handle those events.
    
    lowestModified = Stage::HighestValid;
    shouldTerminate = false;
    std::vector<EventId> eventsForSubsystem;
    for (SubsystemIndex i(0); i<getNSubsystems(); ++i) {
        Stage subsysLowestModified = Stage::HighestValid;
        bool subsysShouldTerminate = false;
        getSystem().getDefaultSubsystem().findSubsystemEventIds(getRep().subsystems[i].getMySubsystemIndex(), restricted, eventIds, eventsForSubsystem);
        if (eventsForSubsystem.size() > 0) {
            getRep().subsystems[i].getSubsystemGuts().handleEvents(restricted, cause, eventsForSubsystem, accuracy, yWeights, ooConstraintTols, subsysLowestModified, subsysShouldTerminate);
            if (subsysLowestModified < lowestModified)
                lowestModified = subsysLowestModified;
            if (subsysShouldTerminate)
                shouldTerminate = true;
        }
    }
    return 0;
}

int System::Guts::reportEventsImpl(const State& s, EventCause cause, const std::vector<EventId>& eventIds) const
{
    // Loop over each subsystem, see which events belong to it, and allow it to handle those events.
    
    std::vector<EventId> eventsForSubsystem;
    for (SubsystemIndex i(0); i<getNSubsystems(); ++i) {
        getSystem().getDefaultSubsystem().findSubsystemEventIds(getRep().subsystems[i].getMySubsystemIndex(), s, eventIds, eventsForSubsystem);
        if (eventsForSubsystem.size() > 0) {
            getRep().subsystems[i].getSubsystemGuts().reportEvents(s, cause, eventsForSubsystem);
        }
    }
    return 0;
}

int System::Guts::calcEventTriggerInfoImpl(const State& s, std::vector<System::EventTriggerInfo>& info) const {

    // Loop over each subsystem, get its EventTriggerInfos, and combine all of them into a single list.
    
    info.clear();
    for (SubsystemIndex i(0); i<getNSubsystems(); ++i) {
        std::vector<System::EventTriggerInfo> subinfo;
        getRep().subsystems[i].getSubsystemGuts().calcEventTriggerInfo(s, subinfo);
        for (std::vector<EventTriggerInfo>::const_iterator e = subinfo.begin(); e != subinfo.end(); e++) {
            info.push_back(*e);
        }
    }
    return 0;
}

int System::Guts::calcTimeOfNextScheduledEventImpl(const State& s, Real& tNextEvent, std::vector<EventId>& eventIds, bool includeCurrentTime) const
{
    tNextEvent = Infinity;
    eventIds.clear();
    for (SubsystemIndex i(0); i<getNSubsystems(); ++i) {
        Real time;
        std::vector<EventId> ids;
        getRep().subsystems[i].getSubsystemGuts().calcTimeOfNextScheduledEvent(s, time, ids, includeCurrentTime);
        if (time < tNextEvent) {
            tNextEvent = time;
            eventIds.clear();
            for (int i = 0; i < (int)ids.size(); ++i)
                eventIds.push_back(ids[i]);
        }
    }
    return 0;
}

int System::Guts::calcTimeOfNextScheduledReportImpl(const State& s, Real& tNextEvent, std::vector<EventId>& eventIds, bool includeCurrentTime) const
{
    tNextEvent = Infinity;
    eventIds.clear();
    for (SubsystemIndex i(0); i<getNSubsystems(); ++i) {
        Real time;
        std::vector<EventId> ids;
        getRep().subsystems[i].getSubsystemGuts().calcTimeOfNextScheduledReport(s, time, ids, includeCurrentTime);
        if (time < tNextEvent) {
            tNextEvent = time;
            eventIds.clear();
            for (int i = 0; i < (int)ids.size(); ++i)
                eventIds.push_back(ids[i]);
        }
    }
    return 0;
}

    ///////////////////////////
    // SYSTEM::GUTS::GUTSREP //
    ///////////////////////////

// All inline currently.



    ////////////////////////
    // EVENT TRIGGER INFO //
    ////////////////////////

System::EventTriggerInfo::EventTriggerInfo() : rep(0) {
    rep = new System::EventTriggerInfoRep(this);
}
System::EventTriggerInfo::~EventTriggerInfo() {
    if (getRep().myHandle == this)
        delete rep;
    rep = 0;
}

System::EventTriggerInfo::EventTriggerInfo(EventId eventId) : rep(0) {
    rep = new System::EventTriggerInfoRep(this);
    rep->eventId = eventId;
}

System::EventTriggerInfo::EventTriggerInfo(const System::EventTriggerInfo& src) : rep(0) {
    rep = new System::EventTriggerInfoRep(src.getRep());
    rep->myHandle = this;
}

System::EventTriggerInfo& 
System::EventTriggerInfo::operator=(const System::EventTriggerInfo& src) {
    if (&src != this) {
        if (getRep().myHandle == this)
            delete rep;
        rep = new System::EventTriggerInfoRep(src.getRep());
        rep->myHandle = this;
    }
    return *this;
}

EventId System::EventTriggerInfo::getEventId() const {
    return getRep().eventId;
}
bool System::EventTriggerInfo::shouldTriggerOnRisingSignTransition() const {
    return getRep().triggerOnRising;
}
bool System::EventTriggerInfo::shouldTriggerOnFallingSignTransition() const {
    return getRep().triggerOnFalling;
}
Real System::EventTriggerInfo::getRequiredLocalizationTimeWindow()    const {
    return getRep().localizationWindow;
}

System::EventTriggerInfo& 
System::EventTriggerInfo::setEventId(EventId id) {
    updRep().eventId = id; 
    return *this;
}
System::EventTriggerInfo& 
System::EventTriggerInfo::setTriggerOnRisingSignTransition(bool shouldTrigger) {
    updRep().triggerOnRising = shouldTrigger; 
    return *this;
}
System::EventTriggerInfo& 
System::EventTriggerInfo::setTriggerOnFallingSignTransition(bool shouldTrigger) {
    updRep().triggerOnFalling = shouldTrigger; 
    return *this;
}
System::EventTriggerInfo& 
System::EventTriggerInfo::setRequiredLocalizationTimeWindow(Real w) {
    assert(w > 0);
    updRep().localizationWindow = w; 
    return *this;
}
    ////////////////////////////
    // EVENT TRIGGER INFO REP //
    ////////////////////////////

// All inline currently.


} // namespace SimTK

