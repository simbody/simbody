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
            System::Guts::destruct(guts);
        guts=0;
        if (src.guts) {
            guts = src.guts->clone();
            guts->setOwnerHandle(*this);
        }
    }
    return *this;
}

System::~System() {
    // Must delete using the library-side VFT, so that we can get access
    // to the client side virtual destructor to destruct this client-side
    // System::Guts object.
    if (guts && isOwnerHandle())
        System::Guts::destruct(guts);
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

const State& System::getDefaultState() const {return getSystemGuts().getDefaultState();}
State& System::updDefaultState() {return updSystemGuts().updDefaultState();}

SubsystemId System::adoptSubsystem(Subsystem& child) {return updSystemGuts().adoptSubsystem(child);}
int System::getNSubsystems() const {return getSystemGuts().getNSubsystems();}
const Subsystem& System::getSubsystem(SubsystemId i) const {return getSystemGuts().getSubsystem(i);}
Subsystem& System::updSubsystem(SubsystemId i) {return updSystemGuts().updSubsystem(i);}
const DefaultSystemSubsystem& System::getDefaultSubsystem() const {
    return static_cast<const DefaultSystemSubsystem&>(getSystemGuts().getSubsystem(SubsystemId(0)));
}
DefaultSystemSubsystem& System::updDefaultSubsystem() {
    return static_cast<DefaultSystemSubsystem&>(updSystemGuts().updSubsystem(SubsystemId(0)));
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
void System::calcDecorativeGeometryAndAppend(const State& s, Stage g, Array<DecorativeGeometry>& geom) const {
    getSystemGuts().calcDecorativeGeometryAndAppend(s,g,geom);
}

Real System::calcTimescale(const State& s) const {return getSystemGuts().calcTimescale(s);}
void System::calcYUnitWeights(const State& s, Vector& weights) const
  { getSystemGuts().calcYUnitWeights(s,weights); }
void System::project(State& s, Real consAccuracy, const Vector& yweights,
                     const Vector& ootols, Vector& yerrest) const
  { getSystemGuts().project(s,consAccuracy,yweights,ootols,yerrest); }
void System::calcYErrUnitTolerances(const State& s, Vector& tolerances) const
  { getSystemGuts().calcYErrUnitTolerances(s,tolerances); }
void System::handleEvents(State& s, EventCause cause, const Array<int>& eventIds,
                          Real accuracy, const Vector& yWeights, const Vector& ooConstraintTols,
                          Stage& lowestModified, bool& shouldTerminate) const
  { getSystemGuts().handleEvents(s,cause,eventIds,accuracy,yWeights,ooConstraintTols,
                               lowestModified,shouldTerminate); }
void System::calcEventTriggerInfo(const State& s, Array<EventTriggerInfo>& info) const
  { getSystemGuts().calcEventTriggerInfo(s,info); }
void System::calcTimeOfNextScheduledEvent(const State& s, Real& tNextEvent,
                                          Array<int>& eventIds) const
  { getSystemGuts().calcTimeOfNextScheduledEvent(s,tNextEvent,eventIds); }

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

// Default constructor is inline, but calls librarySideConstuction() here.
void System::Guts::librarySideConstruction(const String& name, const String& version) {
    rep = new GutsRep(name,version);
    // note that the GutsRep object currently has no owner handle
}

// Destructor is inline, but calls librarySideDestruction() here.
void System::Guts::librarySideDestruction() {
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
const Subsystem& System::Guts::getSubsystem(SubsystemId i) const {return getRep().getSubsystem(i);}
Subsystem& System::Guts::updSubsystem(SubsystemId i) {return updRep().updSubsystem(i);}

void System::Guts::registerDestructImpl(DestructImplLocator f) {
    updRep().destructp = f;
}
void System::Guts::registerCloneImpl(CloneImplLocator f) {
    updRep().clonep = f;
}

void System::Guts::registerRealizeTopologyImpl(RealizeWritableStateImplLocator f) {
    updRep().realizeTopologyp = f;
}
void System::Guts::registerRealizeModelImpl(RealizeWritableStateImplLocator f) {
    updRep().realizeModelp = f;
}
void System::Guts::registerRealizeInstanceImpl(RealizeConstStateImplLocator f) {
    updRep().realizeInstancep = f;
}
void System::Guts::registerRealizeTimeImpl(RealizeConstStateImplLocator f) {
    updRep().realizeTimep = f;
}
void System::Guts::registerRealizePositionImpl(RealizeConstStateImplLocator f) {
    updRep().realizePositionp = f;
}
void System::Guts::registerRealizeVelocityImpl(RealizeConstStateImplLocator f) {
    updRep().realizeVelocityp = f;
}
void System::Guts::registerRealizeDynamicsImpl(RealizeConstStateImplLocator f) {
    updRep().realizeDynamicsp = f;
}
void System::Guts::registerRealizeAccelerationImpl(RealizeConstStateImplLocator f) {
    updRep().realizeAccelerationp = f;
}
void System::Guts::registerRealizeReportImpl(RealizeConstStateImplLocator f) {
    updRep().realizeReportp = f;
}



void System::Guts::registerCalcTimescaleImpl(CalcTimescaleImplLocator f) {
    updRep().calcTimescalep = f;
}
void System::Guts::registerCalcYUnitWeightsImplLocator(CalcUnitWeightsImplLocator f) {
    updRep().calcYUnitWeightsp = f;
}
void System::Guts::registerProjectImpl(ProjectImplLocator f) {
    updRep().projectp = f;
}
void System::Guts::registerCalcYErrUnitTolerancesImplLocator(CalcUnitWeightsImplLocator f) {
    updRep().calcYErrUnitTolerancesp = f;
}
void System::Guts::registerHandleEventsImpl(HandleEventsImplLocator f) {
    updRep().handleEventsp = f;
}
void System::Guts::registerCalcEventTriggerInfoImpl(CalcEventTriggerInfoImplLocator f) {
    updRep().calcEventTriggerInfop = f;
}
void System::Guts::registerCalcTimeOfNextScheduledEventImpl(CalcTimeOfNextScheduledEventImplLocator f) {
    updRep().calcTimeOfNextScheduledEventp = f;
}

System::Guts* System::Guts::clone() const {
    return getRep().clonep(*this);
}


/*static*/void System::Guts::destruct(System::Guts* gutsp) {
    if (gutsp)
        gutsp->getRep().destructp(gutsp);
}

const State& System::Guts::realizeTopology() const {
    State& defaultState = getRep().defaultState; // mutable
    if (!getRep().systemTopologyHasBeenRealized()) {
        defaultState.clear();
        defaultState.setNSubsystems(getNSubsystems());
        for (SubsystemId i(0); i<getNSubsystems(); ++i) 
            defaultState.initializeSubsystem(i, getRep().subsystems[i].getName(), 
                                                getRep().subsystems[i].getVersion());
        
        getRep().realizeTopologyp(*this,defaultState); // defaultState is mutable
        getRep().systemTopologyRealized = true; // mutable
        defaultState.advanceSystemToStage(Stage::Topology);

        // Realize the model using the default settings of the Model variables.
        // This allocates all the later-stage State variables.
        realizeModel(defaultState);

        // Now realize the default state to the highest Stage.
        // TODO: this is problematic if the default state is not a valid one.
        //realize(defaultState, Stage::HighestValid);
    }
    return defaultState;
}

void System::Guts::realizeModel(State& s) const {
    SimTK_STAGECHECK_TOPOLOGY_REALIZED_ALWAYS(systemTopologyHasBeenRealized(),
        "System", getName(), "System::Guts::realizeModel()");
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Topology, 
        "System::Guts::realizeModel()");
    if (s.getSystemStage() < Stage::Model) {
        getRep().realizeModelp(*this,s);    // take care of the Subsystems
        s.advanceSystemToStage(Stage::Model);
    }
}
void System::Guts::realizeInstance(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Instance).prev(), 
        "System::Guts::realizeInstance()");
    if (s.getSystemStage() < Stage::Instance) {
        getRep().realizeInstancep(*this,s);    // take care of the Subsystems
        s.advanceSystemToStage(Stage::Instance);
    }
}
void System::Guts::realizeTime(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Time).prev(), 
        "System::Guts::realizeTime()");
    if (s.getSystemStage() < Stage::Time) {
        getRep().realizeTimep(*this,s);    // take care of the Subsystems
        s.advanceSystemToStage(Stage::Time);
    }
}
void System::Guts::realizePosition(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Position).prev(), 
        "System::Guts::realizePosition()");
    if (s.getSystemStage() < Stage::Position) {
        getRep().realizePositionp(*this,s);    // take care of the Subsystems
        s.advanceSystemToStage(Stage::Position);
    }
}
void System::Guts::realizeVelocity(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Velocity).prev(), 
        "System::Guts::realizeVelocity()");
    if (s.getSystemStage() < Stage::Velocity) {
        getRep().realizeVelocityp(*this,s);    // take care of the Subsystems
        s.advanceSystemToStage(Stage::Velocity);
    }
}
void System::Guts::realizeDynamics(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Dynamics).prev(), 
        "System::Guts::realizeDynamics()");
    if (s.getSystemStage() < Stage::Dynamics) {
        getRep().realizeDynamicsp(*this,s);    // take care of the Subsystems
        s.advanceSystemToStage(Stage::Dynamics);
    }
}
void System::Guts::realizeAcceleration(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Acceleration).prev(), 
        "System::Guts::realizeAcceleration()");
    if (s.getSystemStage() < Stage::Acceleration) {
        getRep().realizeAccelerationp(*this,s);    // take care of the Subsystems
        s.advanceSystemToStage(Stage::Acceleration);
    }
}
void System::Guts::realizeReport(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Report).prev(), 
        "System::Guts::realizeReport()");
    if (s.getSystemStage() < Stage::Report) {
        getRep().realizeReportp(*this,s);    // take care of the Subsystems
        s.advanceSystemToStage(Stage::Report);
    }
}

Real System::Guts::calcTimescale(const State& s) const {
    SimTK_STAGECHECK_GE(s.getSystemStage(), Stage::Instance,
        "System::Guts::calcTimescale()");
    return getRep().calcTimescalep(*this,s);
}

void System::Guts::calcYUnitWeights(const State& s, Vector& weights) const {
    SimTK_STAGECHECK_GE(s.getSystemStage(), Stage::Position,
        "System::Guts::calcYUnitWeights()");
    getRep().calcYUnitWeightsp(*this,s,weights);
}

void System::Guts::calcYErrUnitTolerances(const State& s, Vector& tolerances) const {
    SimTK_STAGECHECK_GE(s.getSystemStage(), Stage::Instance,
        "System::Guts::calcYErrUnitTolerances()");
    getRep().calcYErrUnitTolerancesp(*this,s,tolerances);
}

void System::Guts::project(State& s, Real consAccuracy, const Vector& yweights,
                     const Vector& ootols, Vector& yerrest) const
{
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Instance, // TODO: is this the right stage?
        "System::Guts::project()");
    getRep().projectp(*this,s,consAccuracy,yweights,ootols,yerrest);
}

void System::Guts::handleEvents
   (State& s, EventCause cause, const Array<int>& eventIds,
    Real accuracy, const Vector& yWeights, const Vector& ooConstraintTols,
    Stage& lowestModified, bool& shouldTerminate) const
{
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Model, // TODO: is this the right stage?
        "System::Guts::handleEvents()");
    const Real savedTime = s.getTime();
    getRep().handleEventsp(*this,s,cause,eventIds,accuracy,yWeights,ooConstraintTols,
                     lowestModified, shouldTerminate);
    SimTK_ASSERT_ALWAYS(s.getTime() == savedTime,
        "System::Guts::handleEvents(): handleEventsImpl() tried to change the time");
}

void System::Guts::calcTimeOfNextScheduledEvent
    (const State& s, Real& tNextEvent, Array<int>& eventIds) const
{
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Time,
        "System::Guts::calcTimeOfNextScheduledEvent()");
    tNextEvent = CNT<Real>::getInfinity();
    getRep().calcTimeOfNextScheduledEventp(*this,s,tNextEvent,eventIds);
}

void System::Guts::calcEventTriggerInfo(const State& s, Array<EventTriggerInfo>& info) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Instance,
        "System::Guts::calcEventTriggerInfo()");
    getRep().calcEventTriggerInfop(*this,s,info);
}

void System::Guts::realize(const State& s, Stage g) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Model, 
        "System::Guts::realize()");

    Stage stageNow;
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

void System::Guts::calcDecorativeGeometryAndAppend(const State& s, Stage stage, Array<DecorativeGeometry>& geom) const {
    assert(stage==Stage::Topology || s.getSystemStage() >= stage);
    for (SubsystemId i(0); i<getNSubsystems(); ++i)
        getRep().subsystems[i].getSubsystemGuts().calcDecorativeGeometryAndAppend(s, stage, geom);
}


SubsystemId System::Guts::adoptSubsystem(Subsystem& src) {
    return updRep().adoptSubsystem(src);
}


// These are the default implementations for the System virtual functions.
// Note that this DOES NOT cause binary compatibility problems. The addresses of
// these functions will be supplied from the library side, but these addresses will
// get filled in to the default virtual function table on the *client* side which
// knows where to put each function by name.

int System::Guts::realizeTopologyImpl(State& s) const { 
    for (SubsystemId i(0); i<getNSubsystems(); ++i)
        getRep().subsystems[i].getSubsystemGuts().realizeSubsystemTopology(s);
    return 0;
}
int System::Guts::realizeModelImpl(State& s) const {
    for (SubsystemId i(0); i<getNSubsystems(); ++i)
        getRep().subsystems[i].getSubsystemGuts().realizeSubsystemModel(s);
    return 0;
}
int System::Guts::realizeInstanceImpl(const State& s) const { 
    for (SubsystemId i(0); i<getNSubsystems(); ++i)
        getRep().subsystems[i].getSubsystemGuts().realizeSubsystemInstance(s);
    return 0;
}
int System::Guts::realizeTimeImpl(const State& s) const { 
    for (SubsystemId i(0); i<getNSubsystems(); ++i)
        getRep().subsystems[i].getSubsystemGuts().realizeSubsystemTime(s);
    return 0;
}
int System::Guts::realizePositionImpl(const State& s) const { 
    for (SubsystemId i(0); i<getNSubsystems(); ++i)
        getRep().subsystems[i].getSubsystemGuts().realizeSubsystemPosition(s);
    return 0;
}
int System::Guts::realizeVelocityImpl(const State& s) const { 
    for (SubsystemId i(0); i<getNSubsystems(); ++i)
        getRep().subsystems[i].getSubsystemGuts().realizeSubsystemVelocity(s);
    return 0;
}
int System::Guts::realizeDynamicsImpl(const State& s) const { 
    for (SubsystemId i(0); i<getNSubsystems(); ++i)
        getRep().subsystems[i].getSubsystemGuts().realizeSubsystemDynamics(s);
    return 0;
}
int System::Guts::realizeAccelerationImpl(const State& s) const { 
    for (SubsystemId i(0); i<getNSubsystems(); ++i)
        getRep().subsystems[i].getSubsystemGuts().realizeSubsystemAcceleration(s);
    return 0;
}
int System::Guts::realizeReportImpl(const State& s) const { 
    for (SubsystemId i(0); i<getNSubsystems(); ++i)
        getRep().subsystems[i].getSubsystemGuts().realizeSubsystemReport(s);
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

    for (SubsystemId i(0); i<getNSubsystems(); ++i) {
        const Subsystem::Guts& sub = getRep().subsystems[i].getSubsystemGuts();
        sub.calcQUnitWeights(s, qwts(s.getQStart(i), s.getNQ(i)));
        sub.calcUUnitWeights(s, uwts(s.getUStart(i), s.getNU(i)));
        sub.calcZUnitWeights(s, zwts(s.getZStart(i), s.getNZ(i)));
    }
    return 0;
}

int System::Guts::projectImpl(State&, Real consAccuracy, const Vector& yweights,
                        const Vector& ootols, Vector& yerrest) const
{
    SimTK_THROW2(Exception::UnimplementedVirtualMethod, "System", "project"); 
    return std::numeric_limits<int>::min();
}


int System::Guts::calcYErrUnitTolerancesImpl(const State& s, Vector& ootols) const {
    ootols.resize(s.getNYErr());
    VectorView qtols = ootols(s.getQErrStart(), s.getNQErr()); // writable views
    VectorView utols = ootols(s.getUErrStart(), s.getNUErr());

    for (SubsystemId i(0); i<getNSubsystems(); ++i) {
        const Subsystem::Guts& sub = getRep().subsystems[i].getSubsystemGuts();
        sub.calcQErrUnitTolerances(s, qtols(s.getQErrStart(i), s.getNQErr(i)));
        sub.calcUErrUnitTolerances(s, utols(s.getUErrStart(i), s.getNUErr(i)));
    }
    return 0;
}

int System::Guts::handleEventsImpl
   (State& s, EventCause cause, const Array<int>& eventIds,
    Real accuracy, const Vector& yWeights, const Vector& ooConstraintTols,
    Stage& lowestModified, bool& shouldTerminate) const
{
    // Loop over each subsystem, see which events belong to it, and allow it to handle those events.
    
    lowestModified = Stage::HighestValid;
    shouldTerminate = false;
    Array<int> eventsForSubsystem;
    for (SubsystemId i(0); i<getNSubsystems(); ++i) {
        Stage subsysLowestModified = Stage::HighestValid;
        bool subsysShouldTerminate = false;
        getSystem().getDefaultSubsystem().findSubsystemEventIds(getRep().subsystems[i].getMySubsystemId(), s, eventIds, eventsForSubsystem);
        if (eventsForSubsystem.size() > 0) {
            getRep().subsystems[i].getSubsystemGuts().handleEvents(s, cause, eventsForSubsystem, accuracy, yWeights, ooConstraintTols, subsysLowestModified, subsysShouldTerminate);
            if (subsysLowestModified < lowestModified)
                lowestModified = subsysLowestModified;
            if (subsysShouldTerminate)
                shouldTerminate = true;
        }
    }
    return 0;
}

int System::Guts::calcEventTriggerInfoImpl(const State& s, Array<System::EventTriggerInfo>& info) const {

    // Loop over each subsystem, get its EventTriggerInfos, and combine all of them into a single list.
    
    info.clear();
    for (SubsystemId i(0); i<getNSubsystems(); ++i) {
        Array<System::EventTriggerInfo> subinfo;
        getRep().subsystems[i].getSubsystemGuts().calcEventTriggerInfo(s, subinfo);
        for (Array<EventTriggerInfo>::const_iterator e = subinfo.begin(); e != subinfo.end(); e++) {
            info += *e;
        }
    }
    return 0;
}

int System::Guts::calcTimeOfNextScheduledEventImpl(const State& s, Real& tNextEvent, Array<int>& eventIds) const
{
    tNextEvent = Infinity;
    eventIds.clear();
    for (SubsystemId i(0); i<getNSubsystems(); ++i) {
        Real time;
        Array<int> ids;
        getRep().subsystems[i].getSubsystemGuts().calcTimeOfNextScheduledEvent(s, time, ids);
        if (time < tNextEvent) {
            tNextEvent = time;
            eventIds.clear();
            for (int i = 0; i < ids.size(); ++i)
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

System::EventTriggerInfo::EventTriggerInfo(int eventId) : rep(0) {
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

int System::EventTriggerInfo::getEventId() const {
    return getRep().eventId;
}
bool System::EventTriggerInfo::shouldTriggerOnRisingSignTransition() const {
    return getRep().triggerOnRising;
}
bool System::EventTriggerInfo::shouldTriggerOnFallingSignTransition() const {
    return getRep().triggerOnFalling;
}
bool System::EventTriggerInfo::shouldTriggerOnZeroTransitions()       const {
    return getRep().triggerOnZero;
}
Real System::EventTriggerInfo::getRequiredLocalizationTimeWindow()    const {
    return getRep().localizationWindow;
}

System::EventTriggerInfo& 
System::EventTriggerInfo::setEventId(int id) {
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
System::EventTriggerInfo::setTriggerOnZeroTransitions(bool shouldTrigger) {
    updRep().triggerOnZero = shouldTrigger; 
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

