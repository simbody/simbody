/* Portions copyright (c) 2006-7 Stanford University and Michael Sherman.
 * Contributors:
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */


/**@file
 *
 * Implementation of System and SystemRep.
 */

#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"
#include "SimTKcommon/internal/System.h"
#include "SimTKcommon/internal/Subsystem.h"

#include "SystemRep.h"

#include <cassert>

namespace SimTK {

    ////////////
    // System //
    ////////////

bool System::isEmptyHandle() const {return rep==0;}
bool System::isOwnerHandle() const {return rep==0 || rep->myHandle==this;}
bool System::isSameSystem(const System& otherSystem) const {
    return rep && (rep==otherSystem.rep);
}


void System::librarySideConstruction(const String& name, const String& version) {
    rep = new SystemRep(name,version);
    rep->setMyHandle(*this);
}

void System::librarySideDestruction() {
    if (getRep().myHandle==this)
        delete rep; 
    rep=0;
}

void System::registerRealizeTopologyImpl(RealizeWritableStateImplLocator f) {
    updRep().realizeTopologyp = f;
}
void System::registerRealizeModelImpl(RealizeWritableStateImplLocator f) {
    updRep().realizeModelp = f;
}
void System::registerRealizeInstanceImpl(RealizeConstStateImplLocator f) {
    updRep().realizeInstancep = f;
}
void System::registerRealizeTimeImpl(RealizeConstStateImplLocator f) {
    updRep().realizeTimep = f;
}
void System::registerRealizePositionImpl(RealizeConstStateImplLocator f) {
    updRep().realizePositionp = f;
}
void System::registerRealizeVelocityImpl(RealizeConstStateImplLocator f) {
    updRep().realizeVelocityp = f;
}
void System::registerRealizeDynamicsImpl(RealizeConstStateImplLocator f) {
    updRep().realizeDynamicsp = f;
}
void System::registerRealizeAccelerationImpl(RealizeConstStateImplLocator f) {
    updRep().realizeAccelerationp = f;
}
void System::registerRealizeReportImpl(RealizeConstStateImplLocator f) {
    updRep().realizeReportp = f;
}

void System::registerCalcDecorativeGeometryAndAppendImpl(CalcDecorativeGeometryAndAppendImplLocator f) {
    updRep().calcDecorativeGeometryAndAppendp = f;
}
void System::registerCloneImpl(CloneImplLocator f) {
    updRep().clonep = f;
}

void System::registerCalcTimescaleImpl(CalcTimescaleImplLocator f) {
    updRep().calcTimescalep = f;
}
void System::registerCalcYUnitWeightsImplLocator(CalcUnitWeightsImplLocator f) {
    updRep().calcYUnitWeightsp = f;
}
void System::registerProjectImpl(ProjectImplLocator f) {
    updRep().projectp = f;
}
void System::registerCalcYErrUnitTolerancesImplLocator(CalcUnitWeightsImplLocator f) {
    updRep().calcYErrUnitTolerancesp = f;
}
void System::registerHandleEventsImpl(HandleEventsImplLocator f) {
    updRep().handleEventsp = f;
}
void System::registerCalcEventTriggerInfoImpl(CalcEventTriggerInfoImplLocator f) {
    updRep().calcEventTriggerInfop = f;
}
void System::registerCalcTimeOfNextScheduledEventImpl(CalcTimeOfNextScheduledEventImplLocator f) {
    updRep().calcTimeOfNextScheduledEventp = f;
}

System::System(const System& src) : rep(0) {
    if (src.rep) {
        rep = new SystemRep(*src.rep);
        rep->setMyHandle(*this);
    }
}

System& System::operator=(const System& src) {
    if (&src != this) {
        if (isOwnerHandle()) delete rep; 
        rep=0;
        if (src.rep) {
            rep = new SystemRep(*src.rep);
            rep->setMyHandle(*this);
        }
    }
    return *this;
}


const String& System::getName()    const {return getRep().getName();}
const String& System::getVersion() const {return getRep().getVersion();}

int System::getNSubsystems() const {return getRep().getNSubsystems();}
const Subsystem& System::getSubsystem(SubsystemId i) const {return getRep().getSubsystem(i);}
Subsystem& System::updSubsystem(SubsystemId i) {return updRep().updSubsystem(i);}

// TODO: this should be a Model stage variable allocated by the base class.
// Currently it is just a Topology stage variable stored in the base class.
void System::setHasTimeAdvancedEvents(State&, bool hasEm) const {
    getRep().hasTimeAdvancedEventsFlag = hasEm; // mutable
    getRep().invalidateSystemTopologyCache();
}
bool System::hasTimeAdvancedEvents(const State&) const {
    return getRep().hasTimeAdvancedEventsFlag;
}

const State& System::realizeTopology() const {
    State& defaultState = getRep().defaultState; // mutable
    if (!getRep().systemTopologyHasBeenRealized()) {
        defaultState.clear();
        defaultState.setNSubsystems(getNSubsystems());
        for (SubsystemId i(0); i<getNSubsystems(); ++i) 
            defaultState.initializeSubsystem(i, getRep().subsystems[i].getName(), 
                                                getRep().subsystems[i].getVersion());
        
        realizeTopologyImpl(defaultState); // defaultState is mutable
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

void System::realizeModel(State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Topology, 
        "System::realizeModel()");
    if (s.getSystemStage() < Stage::Model) {
        realizeModelImpl(s);    // take care of the Subsystems
        s.advanceSystemToStage(Stage::Model);
    }
}
void System::realizeInstance(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Instance).prev(), 
        "System::realizeInstance()");
    if (s.getSystemStage() < Stage::Instance) {
        realizeInstanceImpl(s);    // take care of the Subsystems
        s.advanceSystemToStage(Stage::Instance);
    }
}
void System::realizeTime(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Time).prev(), 
        "System::realizeTime()");
    if (s.getSystemStage() < Stage::Time) {
        realizeTimeImpl(s);    // take care of the Subsystems
        s.advanceSystemToStage(Stage::Time);
    }
}
void System::realizePosition(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Position).prev(), 
        "System::realizePosition()");
    if (s.getSystemStage() < Stage::Position) {
        realizePositionImpl(s);    // take care of the Subsystems
        s.advanceSystemToStage(Stage::Position);
    }
}
void System::realizeVelocity(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Velocity).prev(), 
        "System::realizeVelocity()");
    if (s.getSystemStage() < Stage::Velocity) {
        realizeVelocityImpl(s);    // take care of the Subsystems
        s.advanceSystemToStage(Stage::Velocity);
    }
}
void System::realizeDynamics(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Dynamics).prev(), 
        "System::realizeDynamics()");
    if (s.getSystemStage() < Stage::Dynamics) {
        realizeDynamicsImpl(s);    // take care of the Subsystems
        s.advanceSystemToStage(Stage::Dynamics);
    }
}
void System::realizeAcceleration(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Acceleration).prev(), 
        "System::realizeAcceleration()");
    if (s.getSystemStage() < Stage::Acceleration) {
        realizeAccelerationImpl(s);    // take care of the Subsystems
        s.advanceSystemToStage(Stage::Acceleration);
    }
}
void System::realizeReport(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Report).prev(), 
        "System::realizeReport()");
    if (s.getSystemStage() < Stage::Report) {
        realizeReportImpl(s);    // take care of the Subsystems
        s.advanceSystemToStage(Stage::Report);
    }
}

Real System::calcTimescale(const State& s) const {
    SimTK_STAGECHECK_GE(s.getSystemStage(), Stage::Instance,
        "System::calcTimescale()");
    return calcTimescaleImpl(s);
}

void System::calcYUnitWeights(const State& s, Vector& weights) const {
    SimTK_STAGECHECK_GE(s.getSystemStage(), Stage::Position,
        "System::calcYUnitWeights()");
    calcYUnitWeightsImpl(s,weights);
}

void System::calcYErrUnitTolerances(const State& s, Vector& tolerances) const {
    SimTK_STAGECHECK_GE(s.getSystemStage(), Stage::Instance,
        "System::calcYErrUnitTolerances()");
    calcYErrUnitTolerancesImpl(s,tolerances);
}

void System::project(State& s, Real consAccuracy, const Vector& yweights,
                     const Vector& ootols, Vector& yerrest) const
{
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Instance, // TODO: is this the right stage?
        "System::project()");
    projectImpl(s,consAccuracy,yweights,ootols,yerrest);
}

void System::handleEvents
   (State& s, EventCause cause, const Array<int>& eventIds,
    Real accuracy, const Vector& yWeights, const Vector& ooConstraintTols,
    Stage& lowestModified, bool& shouldTerminate) const
{
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Model, // TODO: is this the right stage?
        "System::handleEvents()");
    const Real savedTime = s.getTime();
    handleEventsImpl(s,cause,eventIds,accuracy,yWeights,ooConstraintTols,
                     lowestModified, shouldTerminate);
    SimTK_ASSERT_ALWAYS(s.getTime() == savedTime,
        "System::handleEvents(): handleEventsImpl() tried to change the time");
}

void System::calcTimeOfNextScheduledEvent
    (const State& s, Real& tNextEvent, Array<int>& eventIds) const
{
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Time,
        "System::calcTimeOfNextScheduledEvent()");
    tNextEvent = CNT<Real>::getInfinity();
    calcTimeOfNextScheduledEventImpl(s,tNextEvent,eventIds);
}

void System::calcEventTriggerInfo(const State& s, Array<EventTriggerInfo>& info) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Instance,
        "System::calcEventTriggerInfo()");
    calcEventTriggerInfoImpl(s,info);
}

void System::realize(const State& s, Stage g) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Model, 
        "System::realize()");

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
        default: assert(!"System::realize(): bad stage");
        }
    }
}

void System::calcDecorativeGeometryAndAppend(const State& s, Stage stage, Array<DecorativeGeometry>& geom) const {
    assert(stage==Stage::Topology || s.getSystemStage() >= stage);
    for (SubsystemId i(0); i<getNSubsystems(); ++i)
        getRep().subsystems[i].calcDecorativeGeometryAndAppend(s, stage, geom);
}

bool System::topologyHasBeenRealized() const {
    return getRep().systemTopologyHasBeenRealized();
}

const State& System::getDefaultState() const {
    return getRep().getDefaultState();
}

State& System::updDefaultState() {
    return updRep().updDefaultState();
}


Real System::calcYErrorNorm(const State& s, const Vector& y_err) const {
    return getRep().calcYErrorNorm(s,y_err);
}

SubsystemId System::adoptSubsystem(Subsystem& src) {
    return updRep().adoptSubsystem(src);
}


// These are the default implementations for the System virtual functions.
// Note that this DOES NOT cause binary compatibility problems. The addresses of
// these functions will be supplied from the library side, but these addresses will
// get filled in to the default virtual function table on the *client* side which
// knows where to put each function by name.

int System::realizeTopologyImpl(State& s) const { 
    for (SubsystemId i(0); i<getNSubsystems(); ++i)
        getRep().subsystems[i].realizeSubsystemTopology(s);
    return 0;
}
int System::realizeModelImpl(State& s) const {
    for (SubsystemId i(0); i<getNSubsystems(); ++i)
        getRep().subsystems[i].realizeSubsystemModel(s);
    return 0;
}
int System::realizeInstanceImpl(const State& s) const { 
    for (SubsystemId i(0); i<getNSubsystems(); ++i)
        getRep().subsystems[i].realizeSubsystemInstance(s);
    return 0;
}
int System::realizeTimeImpl(const State& s) const { 
    for (SubsystemId i(0); i<getNSubsystems(); ++i)
        getRep().subsystems[i].realizeSubsystemTime(s);
    return 0;
}
int System::realizePositionImpl(const State& s) const { 
    for (SubsystemId i(0); i<getNSubsystems(); ++i)
        getRep().subsystems[i].realizeSubsystemPosition(s);
    return 0;
}
int System::realizeVelocityImpl(const State& s) const { 
    for (SubsystemId i(0); i<getNSubsystems(); ++i)
        getRep().subsystems[i].realizeSubsystemVelocity(s);
    return 0;
}
int System::realizeDynamicsImpl(const State& s) const { 
    for (SubsystemId i(0); i<getNSubsystems(); ++i)
        getRep().subsystems[i].realizeSubsystemDynamics(s);
    return 0;
}
int System::realizeAccelerationImpl(const State& s) const { 
    for (SubsystemId i(0); i<getNSubsystems(); ++i)
        getRep().subsystems[i].realizeSubsystemAcceleration(s);
    return 0;
}
int System::realizeReportImpl(const State& s) const { 
    for (SubsystemId i(0); i<getNSubsystems(); ++i)
        getRep().subsystems[i].realizeSubsystemReport(s);
    return 0;
}

int System::calcDecorativeGeometryAndAppendImpl
   (const State&, Stage, Array<DecorativeGeometry>&) const
{
    return 0;
}

System* System::cloneImpl() const {
    return new System(*this);
}

Real System::calcTimescaleImpl(const State&) const {
    return 0.1; // TODO!!!
}

int System::calcYUnitWeightsImpl(const State& s, Vector& weights) const {
    weights.resize(s.getNY());
    VectorView qwts = weights(s.getQStart(), s.getNQ());   // writable views
    VectorView uwts = weights(s.getUStart(), s.getNU());
    VectorView zwts = weights(s.getZStart(), s.getNZ());

    for (SubsystemId i(0); i<getNSubsystems(); ++i) {
        const Subsystem& sub = getRep().subsystems[i];
        sub.calcQUnitWeights(s, qwts(s.getQStart(i), s.getNQ(i)));
        sub.calcUUnitWeights(s, uwts(s.getUStart(i), s.getNU(i)));
        sub.calcZUnitWeights(s, zwts(s.getZStart(i), s.getNZ(i)));
    }
    return 0;
}

int System::projectImpl(State&, Real consAccuracy, const Vector& yweights,
                        const Vector& ootols, Vector& yerrest) const
{
    SimTK_THROW2(Exception::UnimplementedVirtualMethod, "System", "project"); 
    return std::numeric_limits<int>::min();
}


int System::calcYErrUnitTolerancesImpl(const State& s, Vector& ootols) const {
    ootols.resize(s.getNYErr());
    VectorView qtols = ootols(s.getQErrStart(), s.getNQErr()); // writable views
    VectorView utols = ootols(s.getUErrStart(), s.getNUErr());

    for (SubsystemId i(0); i<getNSubsystems(); ++i) {
        const Subsystem& sub = getRep().subsystems[i];
        sub.calcQErrUnitTolerances(s, qtols(s.getQErrStart(i), s.getNQErr(i)));
        sub.calcUErrUnitTolerances(s, utols(s.getUErrStart(i), s.getNUErr(i)));
    }
    return 0;
}

int System::handleEventsImpl
   (State&, EventCause, const Array<int>& eventIds,
    Real accuracy, const Vector& yWeights, const Vector& ooConstraintTols,
    Stage& lowestModified, bool& shouldTerminate) const
{
    SimTK_THROW2(Exception::UnimplementedVirtualMethod, "System", "handleEvents"); 
    return std::numeric_limits<int>::min();
}

int System::calcEventTriggerInfoImpl(const State& s, Array<System::EventTriggerInfo>& info) const {
    info.resize(s.getNEvents());
    for (int i=0; i<info.size(); ++i) 
        info[i] = EventTriggerInfo(i);
    return 0;
}

int System::calcTimeOfNextScheduledEventImpl
    (const State&, Real& tNextEvent, Array<int>& eventIds) const
{
    tNextEvent = NTraits<Real>::Infinity;
    eventIds.clear();
    return 0;
}


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


void System::adoptPrivateImplementation(System::PrivateImplementation* p,
                                        System::ClonePrivateImplementation clone,
                                        System::DestructPrivateImplementation destruct)
{
    updRep().adoptPrivateImplementation(p,clone,destruct);
}

const System::PrivateImplementation& System::getPrivateImplementation() const {
    return getRep().getPrivateImplementation();
}

System::PrivateImplementation& System::updPrivateImplementation() {
    return updRep().updPrivateImplementation();
}

    // IMPLEMENTATION OF EVENT TRIGGER INFO

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

    ////////////////
    // SYSTEM REP //
    ////////////////



} // namespace SimTK

