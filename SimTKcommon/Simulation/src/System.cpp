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
#include "SimTKcommon/internal/EventHandler.h"
#include "SimTKcommon/internal/EventReporter.h"

#include "SystemGutsRep.h"

#include <cassert>
#include <map>
#include <set>

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
int System::getNumQProjections() const {return getSystemGuts().getRep().nQProjections;}
int System::getNumUProjections() const {return getSystemGuts().getRep().nUProjections;} 
int System::getNumQErrorEstimateProjections() const {return getSystemGuts().getRep().nQErrEstProjections;}
int System::getNumUErrorEstimateProjections() const {return getSystemGuts().getRep().nUErrEstProjections;}
int System::getNumPrescribeCalls() const {return getSystemGuts().getRep().nPrescribeCalls;}
int System::getNumProjectCalls() const {return getSystemGuts().getRep().nProjectCalls;}
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
bool System::systemTopologyHasBeenRealized() const {
    return getSystemGuts().systemTopologyHasBeenRealized();
}
const State& System::realizeTopology() const {return getSystemGuts().realizeTopology();}
void System::realizeModel(State& s) const {getSystemGuts().realizeModel(s);}
void System::realize(const State& s, Stage g) const {getSystemGuts().realize(s,g);}
void System::calcDecorativeGeometryAndAppend(const State& s, Stage g, Array_<DecorativeGeometry>& geom) const {
    getSystemGuts().calcDecorativeGeometryAndAppend(s,g,geom);
}

Real System::calcTimescale(const State& s) const {return getSystemGuts().calcTimescale(s);}
void System::calcYUnitWeights(const State& s, Vector& weights) const
{   getSystemGuts().calcYUnitWeights(s,weights); }
void System::prescribe(State& s, Stage g) const
{   getSystemGuts().prescribe(s,g); }
void System::project(State& s, Real consAccuracy, const Vector& yweights,
                     const Vector& ootols, Vector& yerrest, System::ProjectOptions opts) const
{   getSystemGuts().project(s,consAccuracy,yweights,ootols,yerrest, opts); }
void System::calcYErrUnitTolerances(const State& s, Vector& tolerances) const
{   getSystemGuts().calcYErrUnitTolerances(s,tolerances); }
void System::handleEvents(State& s, Event::Cause cause, const Array_<EventId>& eventIds,
                          Real accuracy, const Vector& yWeights, const Vector& ooConstraintTols,
                          Stage& lowestModified, bool& shouldTerminate) const
{   getSystemGuts().handleEvents(s,cause,eventIds,accuracy,yWeights,ooConstraintTols,
                               lowestModified,shouldTerminate); }
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



    //////////////////
    // SYSTEM::GUTS //
    //////////////////

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

int System::Guts::getNumSubsystems() const {return getRep().getNumSubsystems();}
const Subsystem& System::Guts::getSubsystem(SubsystemIndex i) const {return getRep().getSubsystem(i);}
Subsystem& System::Guts::updSubsystem(SubsystemIndex i) {return updRep().updSubsystem(i);}

System::Guts* System::Guts::clone() const {
    return cloneImpl();
}

const State& System::Guts::realizeTopology() const {
    State& defaultState = getRep().defaultState; // mutable
    if (!getRep().systemTopologyHasBeenRealized()) {
        defaultState.clear();
        defaultState.setNumSubsystems(getNumSubsystems());
        for (SubsystemIndex i(0); i<getNumSubsystems(); ++i) 
            defaultState.initializeSubsystem(i, getRep().subsystems[i].getName(), 
                                                getRep().subsystems[i].getVersion());
        
        // Allow the subclass to do processing.
        realizeTopologyImpl(defaultState); // defaultState is mutable
        // Realize any subsystems that the subclass didn't already take care of.
        for (SubsystemIndex i(0); i<getNumSubsystems(); ++i)
            if (getRep().subsystems[i].getStage(defaultState) < Stage::Topology)
                getRep().subsystems[i].getSubsystemGuts().realizeSubsystemTopology(defaultState);
        getRep().systemTopologyRealized = true; // mutable
        defaultState.advanceSystemToStage(Stage::Topology);
        getRep().nRealizationsOfStage[Stage::Topology]++; // mutable counter

        // Realize the model using the default settings of the Model variables.
        // This allocates all the later-stage State variables.
        realizeModel(defaultState);
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
        for (SubsystemIndex i(0); i<getNumSubsystems(); ++i)
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
        for (SubsystemIndex i(0); i<getNumSubsystems(); ++i)
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
        for (SubsystemIndex i(0); i<getNumSubsystems(); ++i)
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
        for (SubsystemIndex i(0); i<getNumSubsystems(); ++i)
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
        for (SubsystemIndex i(0); i<getNumSubsystems(); ++i)
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
        for (SubsystemIndex i(0); i<getNumSubsystems(); ++i)
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
        for (SubsystemIndex i(0); i<getNumSubsystems(); ++i)
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
        for (SubsystemIndex i(0); i<getNumSubsystems(); ++i)
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

void System::Guts::prescribe(State& s, Stage g) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), g.prev(),
                               "System::Guts::prescribe()");
    prescribeImpl(s,g);
    getRep().nPrescribeCalls++; // mutable counter
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
   (State& s, Event::Cause cause, const Array_<EventId>& eventIds,
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
   (const State& s, Event::Cause cause, const Array_<EventId>& eventIds) const
{
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Model, // TODO: is this the right stage?
        "System::Guts::reportEvents()");
    const Real savedTime = s.getTime();
    reportEventsImpl(s,cause,eventIds);

    getRep().nReportEventsCalls++; // mutable counter
    getRep().nHandlerCallsThatChangedStage[Stage::Report]++;
}

void System::Guts::calcTimeOfNextScheduledEvent
    (const State& s, Real& tNextEvent, Array_<EventId>& eventIds, bool includeCurrentTime) const
{
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Time,
        "System::Guts::calcTimeOfNextScheduledEvent()");
    tNextEvent = Infinity;
    calcTimeOfNextScheduledEventImpl(s,tNextEvent,eventIds,includeCurrentTime);
}

void System::Guts::calcTimeOfNextScheduledReport
    (const State& s, Real& tNextEvent, Array_<EventId>& eventIds, bool includeCurrentTime) const
{
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Time,
        "System::Guts::calcTimeOfNextScheduledReport()");
    tNextEvent = Infinity;
    calcTimeOfNextScheduledReportImpl(s,tNextEvent,eventIds,includeCurrentTime);
}

void System::Guts::calcEventTriggerInfo(const State& s, Array_<EventTriggerInfo>& info) const {
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

void System::Guts::calcDecorativeGeometryAndAppend(const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const {
    assert(stage==Stage::Topology || s.getSystemStage() >= stage);
    for (SubsystemIndex i(0); i<getNumSubsystems(); ++i)
        getRep().subsystems[i].getSubsystemGuts().calcDecorativeGeometryAndAppend(s, stage, geom);
}


SubsystemIndex System::Guts::adoptSubsystem(Subsystem& src) {
    return updRep().adoptSubsystem(src);
}

// Default implementations of virtual methods. These quietly succeed
// at doing nothing.

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

    for (SubsystemIndex i(0); i<getNumSubsystems(); ++i) {
        const Subsystem::Guts& sub = getRep().subsystems[i].getSubsystemGuts();
        sub.calcQUnitWeights(s, qwts(s.getQStart(i), s.getNQ(i)));
        sub.calcUUnitWeights(s, uwts(s.getUStart(i), s.getNU(i)));
        sub.calcZUnitWeights(s, zwts(s.getZStart(i), s.getNZ(i)));
    }
    return 0;
}

int System::Guts::prescribeImpl(State&, Stage) const
{
    SimTK_THROW2(Exception::UnimplementedVirtualMethod, "System", "prescribe"); 
    return std::numeric_limits<int>::min();
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

    for (SubsystemIndex i(0); i<getNumSubsystems(); ++i) {
        const Subsystem::Guts& sub = getRep().subsystems[i].getSubsystemGuts();
        sub.calcQErrUnitTolerances(s, qtols(s.getQErrStart(i), s.getNQErr(i)));
        sub.calcUErrUnitTolerances(s, utols(s.getUErrStart(i), s.getNUErr(i)));
    }
    return 0;
}

int System::Guts::handleEventsImpl
   (State& state, Event::Cause cause, const Array_<EventId>& eventIds,
    Real accuracy, const Vector& yWeights, const Vector& ooConstraintTols,
    Stage& lowestModified, bool& shouldTerminate) const
{
    // Event handlers should not be able to modify the time. We'll
    // save time here and restore it at the end if it might have been
    // changed.
    const Real savedT = state.getTime();

    // Loop over each subsystem, see which events belong to it, and 
    // allow it to handle those events.
    
    lowestModified = Stage::Infinity;
    shouldTerminate = false;
    Array_<EventId> eventsForSubsystem;
    for (SubsystemIndex i(0); i<getNumSubsystems(); ++i) {
        Stage subsysLowestModified = Stage::Infinity;
        bool subsysShouldTerminate = false;
        getSystem().getDefaultSubsystem().findSubsystemEventIds
           (getRep().subsystems[i].getMySubsystemIndex(), state, eventIds, 
            eventsForSubsystem);
        if (eventsForSubsystem.size() > 0) {
            getRep().subsystems[i].getSubsystemGuts().handleEvents
               (state, cause, eventsForSubsystem, accuracy, yWeights, 
                ooConstraintTols, subsysLowestModified, subsysShouldTerminate);
            if (subsysLowestModified < lowestModified)
                lowestModified = subsysLowestModified;
            if (subsysShouldTerminate)
                shouldTerminate = true;
        }
    }

    if (lowestModified <= Stage::Time)
        state.setTime(savedT);

    return 0;
}

int System::Guts::reportEventsImpl(const State& s, Event::Cause cause, const Array_<EventId>& eventIds) const
{
    // Loop over each subsystem, see which events belong to it, and allow it to handle those events.
    
    Array_<EventId> eventsForSubsystem;
    for (SubsystemIndex i(0); i<getNumSubsystems(); ++i) {
        getSystem().getDefaultSubsystem().findSubsystemEventIds(getRep().subsystems[i].getMySubsystemIndex(), s, eventIds, eventsForSubsystem);
        if (eventsForSubsystem.size() > 0) {
            getRep().subsystems[i].getSubsystemGuts().reportEvents(s, cause, eventsForSubsystem);
        }
    }
    return 0;
}

int System::Guts::calcEventTriggerInfoImpl(const State& s, Array_<EventTriggerInfo>& info) const {

    // Loop over each subsystem, get its EventTriggerInfos, and combine all of them into a single list.
    
    info.clear();
    for (SubsystemIndex i(0); i<getNumSubsystems(); ++i) {
        Array_<EventTriggerInfo> subinfo;
        getRep().subsystems[i].getSubsystemGuts().calcEventTriggerInfo(s, subinfo);
        for (Array_<EventTriggerInfo>::const_iterator e = subinfo.begin(); e != subinfo.end(); e++) {
            info.push_back(*e);
        }
    }
    return 0;
}

int System::Guts::calcTimeOfNextScheduledEventImpl(const State& s, Real& tNextEvent, Array_<EventId>& eventIds, bool includeCurrentTime) const
{
    tNextEvent = Infinity;
    eventIds.clear();
    for (SubsystemIndex i(0); i<getNumSubsystems(); ++i) {
        Real time;
        Array_<EventId> ids;
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

int System::Guts::calcTimeOfNextScheduledReportImpl(const State& s, Real& tNextEvent, Array_<EventId>& eventIds, bool includeCurrentTime) const
{
    tNextEvent = Infinity;
    eventIds.clear();
    for (SubsystemIndex i(0); i<getNumSubsystems(); ++i) {
        Real time;
        Array_<EventId> ids;
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


    //////////////////////////////
    // DEFAULT SYSTEM SUBSYSTEM //
    //////////////////////////////

class DefaultSystemSubsystem::Guts : public Subsystem::Guts {
public:
    /**
     * This class stores various information used by the default subsystem that needs to be stored in the state.
     */

    class CacheInfo {
    public:
        CacheInfo() : eventIdCounter(0) {}
        mutable int eventIdCounter;
        mutable std::map<int, SubsystemIndex> eventOwnerMap;
        Array_<EventId> scheduledEventIds;
        Array_<EventTriggerByStageIndex> triggeredEventIndices;
        Array_<EventId> triggeredEventIds;
        Array_<EventId> scheduledReportIds;
        Array_<EventTriggerByStageIndex> triggeredReportIndices;
        Array_<EventId> triggeredReportIds;
    };
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
    
    Guts* cloneImpl() const {
        return new Guts(*this);
    }
        
    const Array_<ScheduledEventHandler*>& getScheduledEventHandlers() const {
        return scheduledEventHandlers;
    }
    
    Array_<ScheduledEventHandler*>& updScheduledEventHandlers() {
        return scheduledEventHandlers;
    }
    
    const Array_<TriggeredEventHandler*>& getTriggeredEventHandlers() const {
        return triggeredEventHandlers;
    }
    
    Array_<TriggeredEventHandler*>& updTriggeredEventHandlers() {
        return triggeredEventHandlers;
    }
    
    const Array_<ScheduledEventReporter*>& getScheduledEventReporters() const {
        return scheduledEventReporters;
    }
    
    Array_<ScheduledEventReporter*>& updScheduledEventReporters() const {
        return scheduledEventReporters;
    }
    
    const Array_<TriggeredEventReporter*>& getTriggeredEventReporters() const {
        return triggeredEventReporters;
    }
    
    Array_<TriggeredEventReporter*>& updTriggeredEventReporters() const {
        return triggeredEventReporters;
    }
    
    const CacheInfo& getCacheInfo(const State& s) const {
        return Value<CacheInfo>::downcast(getCacheEntry(s, cacheInfoIndex)).get();
    }
    
    CacheInfo& updCacheInfo(const State& s) const {
        return Value<CacheInfo>::downcast(updCacheEntry(s, cacheInfoIndex)).upd();
    }

    int realizeSubsystemTopologyImpl(State& s) const {
        cacheInfoIndex = s.allocateCacheEntry(getMySubsystemIndex(), Stage::Topology, new Value<CacheInfo>());
        return 0;
    }
    
    int realizeSubsystemModelImpl(State& s) const {
        return 0;
    }

    int realizeEventTriggers(const State& s, Stage g) const {
        const CacheInfo& info = getCacheInfo(s);
        Vector& triggers = s.updEventTriggersByStage(getMySubsystemIndex(), g);
        for (int i = 0; i < (int)triggeredEventHandlers.size(); ++i) {
            if (g == triggeredEventHandlers[i]->getRequiredStage())
                triggers[info.triggeredEventIndices[i]] = triggeredEventHandlers[i]->getValue(s);
        }
        for (int i = 0; i < (int)triggeredEventReporters.size(); ++i) {
            if (g == triggeredEventReporters[i]->getRequiredStage())
                triggers[info.triggeredReportIndices[i]] = triggeredEventReporters[i]->getValue(s);
        }
        return 0;
    }
    
    int realizeSubsystemInstanceImpl(const State& s) const {
        CacheInfo& info = updCacheInfo(s);
        info.scheduledEventIds.clear();
        info.triggeredEventIndices.clear();
        info.triggeredEventIds.clear();
        info.scheduledReportIds.clear();
        info.triggeredReportIndices.clear();
        info.triggeredReportIds.clear();
        info.eventIdCounter = 0;
        if (scheduledEventHandlers.size() > 0)
            for (Array_<ScheduledEventHandler*>::const_iterator e = scheduledEventHandlers.begin(); e != scheduledEventHandlers.end(); e++) {
                EventId id;
                createScheduledEvent(s, id);
                info.scheduledEventIds.push_back(id);
            }
        if (triggeredEventHandlers.size() > 0)
            for (Array_<TriggeredEventHandler*>::const_iterator e = triggeredEventHandlers.begin(); e != triggeredEventHandlers.end(); e++) {
                EventId id;
                EventTriggerByStageIndex index;
                createTriggeredEvent(s, id, index, (*e)->getRequiredStage());
                info.triggeredEventIds.push_back(id);
                info.triggeredEventIndices.push_back(index);
            }
        if (scheduledEventReporters.size() > 0)
            for (Array_<ScheduledEventReporter*>::const_iterator e = scheduledEventReporters.begin(); e != scheduledEventReporters.end(); e++) {
                EventId id;
                createScheduledEvent(s, id);
                info.scheduledReportIds.push_back(id);
            }
        if (triggeredEventReporters.size() > 0)
            for (Array_<TriggeredEventReporter*>::const_iterator e = triggeredEventReporters.begin(); e != triggeredEventReporters.end(); e++) {
                EventId id;
                EventTriggerByStageIndex index;
                createTriggeredEvent(s, id, index, (*e)->getRequiredStage());
                info.triggeredReportIds.push_back(id);
                info.triggeredReportIndices.push_back(index);
            }
        return 0;
    }
    int realizeSubsystemTimeImpl(const State& s) const {
        return realizeEventTriggers(s, Stage::Time);
    }
    int realizeSubsystemPositionImpl(const State& s) const {
        return realizeEventTriggers(s, Stage::Position);
    }
    int realizeSubsystemVelocityImpl(const State& s) const {
        return realizeEventTriggers(s, Stage::Velocity);
    }
    int realizeSubsystemDynamicsImpl(const State& s) const {
        return realizeEventTriggers(s, Stage::Dynamics);
    }
    int realizeSubsystemAccelerationImpl(const State& s) const {
        return realizeEventTriggers(s, Stage::Acceleration);
    }
    int realizeSubsystemReportImpl(const State& s) const {
        return realizeEventTriggers(s, Stage::Report);
    }
    void calcEventTriggerInfo(const State& s, Array_<EventTriggerInfo>& trigInfo) const {
        
        // Loop over all registered TriggeredEventHandlers and TriggeredEventReporters, and ask
        // each one for its EventTriggerInfo.
        
        const CacheInfo& info = getCacheInfo(s);
        trigInfo.resize(triggeredEventHandlers.size()+triggeredEventReporters.size());
        for (int i = 0; i < (int)triggeredEventHandlers.size(); ++i) {
            Stage stage = triggeredEventHandlers[i]->getRequiredStage();
            SystemEventTriggerIndex index = SystemEventTriggerIndex
                                                     (info.triggeredEventIndices[i]
                                                      + s.getEventTriggerStartByStage(stage)
                                                      + s.getEventTriggerStartByStage(getMySubsystemIndex(), stage));
            trigInfo[index] = triggeredEventHandlers[i]->getTriggerInfo();
            trigInfo[index].setEventId(info.triggeredEventIds[i]);
        }
        for (int i = 0; i < (int)triggeredEventReporters.size(); ++i) {
            Stage stage = triggeredEventReporters[i]->getRequiredStage();
            SystemEventTriggerIndex index = SystemEventTriggerIndex
                                                     (info.triggeredReportIndices[i]
                                                      + s.getEventTriggerStartByStage(stage)
                                                      + s.getEventTriggerStartByStage(getMySubsystemIndex(), stage));
            trigInfo[index] = triggeredEventReporters[i]->getTriggerInfo();
            trigInfo[index].setEventId(info.triggeredReportIds[i]);
        }
    }
    void calcTimeOfNextScheduledEvent(const State& s, Real& tNextEvent, Array_<EventId>& eventIds, bool includeCurrentTime) const {
        
        // Loop over all registered ScheduledEventHandlers, and ask each one when its next event occurs.
        
        const CacheInfo& info = getCacheInfo(s);
        tNextEvent = Infinity;
        for (int i = 0; i < (int)scheduledEventHandlers.size(); ++i) {
            Real time = scheduledEventHandlers[i]->getNextEventTime(s, includeCurrentTime);
            if (time <= tNextEvent && (time > s.getTime() || (includeCurrentTime && time == s.getTime()))) {
                if (time < tNextEvent)
                    eventIds.clear();
                tNextEvent = time;
                eventIds.push_back(info.scheduledEventIds[i]);
            }
        }
    }
    void calcTimeOfNextScheduledReport(const State& s, Real& tNextEvent, Array_<EventId>& eventIds, bool includeCurrentTime) const {
        
        // Loop over all registered ScheduledEventReporters, and ask each one when its next event occurs.
        
        const CacheInfo& info = getCacheInfo(s);
        tNextEvent = Infinity;
        for (int i = 0; i < (int)scheduledEventReporters.size(); ++i) {
            Real time = scheduledEventReporters[i]->getNextEventTime(s, includeCurrentTime);
            if (time <= tNextEvent && (time > s.getTime() || (includeCurrentTime && time == s.getTime()))) {
                if (time < tNextEvent)
                    eventIds.clear();
                tNextEvent = time;
                eventIds.push_back(info.scheduledReportIds[i]);
            }
        }
    }
    void handleEvents(State& s, Event::Cause cause, const Array_<EventId>& eventIds, Real accuracy, const Vector& yWeights,
            const Vector& ooConstraintTols, Stage& lowestModified, bool& shouldTerminate) const {
        const CacheInfo& info = getCacheInfo(s);
        lowestModified = Stage::Infinity;
        shouldTerminate = false;
        
        // Build a set of the ids for quick lookup.
        
        std::set<EventId> idSet;
        for (int i = 0; i < (int)eventIds.size(); ++i)
            idSet.insert(eventIds[i]);
        
        // Process triggered events.
        
        if (cause == Event::Cause::Triggered) {
            for (int i = 0; i < (int)triggeredEventHandlers.size(); ++i) {
                if (idSet.find(info.triggeredEventIds[i]) != idSet.end()) {
                    Stage eventLowestModified = Stage::Infinity;
                    bool eventShouldTerminate = false;
                    triggeredEventHandlers[i]->handleEvent(s, accuracy, yWeights, ooConstraintTols, eventLowestModified, eventShouldTerminate);
                    if (eventLowestModified < lowestModified)
                        lowestModified = eventLowestModified;
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
                    Stage eventLowestModified = Stage::Infinity;
                    bool eventShouldTerminate = false;
                    scheduledEventHandlers[i]->handleEvent(s, accuracy, yWeights, ooConstraintTols, eventLowestModified, eventShouldTerminate);
                    if (eventLowestModified < lowestModified)
                        lowestModified = eventLowestModified;
                    if (eventShouldTerminate)
                        shouldTerminate = true;
                }
            }
            for (int i = 0; i < (int)scheduledEventReporters.size(); ++i) {
                if (idSet.find(info.scheduledReportIds[i]) != idSet.end())
                    scheduledEventReporters[i]->handleEvent(s);
            }
        }
    }

    void reportEvents(const State& s, Event::Cause cause, const Array_<EventId>& eventIds) const {
        const CacheInfo& info = getCacheInfo(s);
        
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
    mutable CacheEntryIndex                 cacheInfoIndex;
    mutable Array_<ScheduledEventHandler*>  scheduledEventHandlers;
    mutable Array_<TriggeredEventHandler*>  triggeredEventHandlers;
    mutable Array_<ScheduledEventReporter*> scheduledEventReporters;
    mutable Array_<TriggeredEventReporter*> triggeredEventReporters;
};

std::ostream& operator<<(std::ostream& o, const DefaultSystemSubsystem::Guts::CacheInfo& info) {
    o << "DefaultSystemSubsystemGuts::CacheInfo";
    return o;
}

DefaultSystemSubsystem::DefaultSystemSubsystem(System& sys) {
    adoptSubsystemGuts(new DefaultSystemSubsystem::Guts());
    sys.adoptSubsystem(*this);
}

const DefaultSystemSubsystem::Guts& DefaultSystemSubsystem::getGuts() const {
    return dynamic_cast<const DefaultSystemSubsystem::Guts&>(getSubsystemGuts());
}

DefaultSystemSubsystem::Guts& DefaultSystemSubsystem::updGuts() {
    return dynamic_cast<DefaultSystemSubsystem::Guts&>(updSubsystemGuts());
}

/**
 * Add a ScheduledEventHandler to the System.  This must be called before the Model stage is realized.
 * 
 * The System assumes ownership of the object passed to this method, and will delete it when the
 * System is deleted.
 */

void DefaultSystemSubsystem::addEventHandler(ScheduledEventHandler* handler) {
    updGuts().updScheduledEventHandlers().push_back(handler);
}

/**
 * Add a TriggeredEventHandler to the System.  This must be called before the Model stage is realized.
 * 
 * The System assumes ownership of the object passed to this method, and will delete it when the
 * System is deleted.
 */

void DefaultSystemSubsystem::addEventHandler(TriggeredEventHandler* handler) {
    updGuts().updTriggeredEventHandlers().push_back(handler);
}

/**
 * Add a ScheduledEventReporter to the System.  This must be called before the Model stage is realized.
 * 
 * The System assumes ownership of the object passed to this method, and will delete it when the
 * System is deleted.
 * 
 * Note that this method is const.  Because an EventReporter cannot affect the behavior of the system
 * being simulated, it is permitted to add one to a const System.
 */

void DefaultSystemSubsystem::addEventReporter(ScheduledEventReporter* handler) const {
    getGuts().updScheduledEventReporters().push_back(handler);
}

/**
 * Add a TriggeredEventReporter to the System.  This must be called before the Model stage is realized.
 * 
 * The System assumes ownership of the object passed to this method, and will delete it when the
 * System is deleted.
 * 
 * Note that this method is const.  Because an EventReporter cannot affect the behavior of the system
 * being simulated, it is permitted to add one to a const System.
 */

void DefaultSystemSubsystem::addEventReporter(TriggeredEventReporter* handler) const {
    getGuts().updTriggeredEventReporters().push_back(handler);
}

/**
 * Generate a new, globally unique event ID.  Typically you will not call this directly.  When a Subsystem
 * needs to obtain an event ID for an event it defines, it should do so by calling
 * Subsystem::Guts::createScheduledEvent() or Subsystem::Guts::createTriggeredEvent().
 */

EventId DefaultSystemSubsystem::createEventId(SubsystemIndex subsys, const State& state) const {
    const Guts::CacheInfo& info = getGuts().getCacheInfo(state);
    int id = info.eventIdCounter++;
    info.eventOwnerMap[id] = subsys;
    return EventId(id);
}

/**
 * Given a list of event IDs, filter it to produce a list of those events belonging to a particular Subsystem.
 * 
 * @param subsys       the Subsystem for which to find events
 * @param state        the State which produced the events
 * @param allEvents    a list of event IDs to filter
 * @param eventsForSubsystem    on exit, this Array_ contains the filtered list of event IDs belonging to the
 *                              specified Subsystem.
 */

void DefaultSystemSubsystem::findSubsystemEventIds
   (SubsystemIndex subsys, const State& state, const Array_<EventId>& allEvents, 
    Array_<EventId>& eventsForSubsystem) const 
{
    const Guts::CacheInfo& info = getGuts().getCacheInfo(state);
    eventsForSubsystem.clear();
    for (int i = 0; i < (int)allEvents.size(); ++i) {
        if (info.eventOwnerMap[allEvents[i]] == subsys)
            eventsForSubsystem.push_back(allEvents[i]);
    }
}


} // namespace SimTK

