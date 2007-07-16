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

System::~System() {
    if (isOwnerHandle()) delete rep; 
    rep=0;
}

System::System(const System& src) : rep(0) {
    if (src.rep) {
        rep = src.rep->clone();
        rep->setMyHandle(*this);
    }
}

System& System::operator=(const System& src) {
    if (&src != this) {
        if (isOwnerHandle()) delete rep; 
        rep=0;
        if (src.rep) {
            rep = src.rep->clone();
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

const State& System::realizeTopology() const {
    return getRep().realizeTopology();
}

bool System::topologyHasBeenRealized() const {
    return getRep().systemTopologyHasBeenRealized();
}

const State& System::getDefaultState() const {
    return getRep().getDefaultState();
}

void System::realizeModel(State& s) const {
    return getRep().realizeModel(s);
}

void System::realize(const State& s, Stage g) const {
    getRep().realize(s,g);
}

void System::calcDecorativeGeometryAndAppend(const State& s, Stage g, Array<DecorativeGeometry>& dg) const {
    getRep().calcDecorativeGeometryAndAppend(s,g,dg);
}

Real System::calcTimescale(const State& s) const {
    return getRep().calcTimescale(s);
}

void System::calcYUnitWeights(const State& s, Vector& weights) const {
    return getRep().calcYUnitWeights(s,weights);
}

void System::calcYErrUnitTolerances(const State& s, Vector& tolerances) const {
    return getRep().calcYErrUnitTolerances(s,tolerances);
}

Real System::calcYErrorNorm(const State& s, const Vector& y_err) const {
    return getRep().calcYErrorNorm(s,y_err);
}

SubsystemId System::adoptSubsystem(Subsystem& src) {
    return updRep().adoptSubsystem(src);
}

    ////////////////
    // SYSTEM REP //
    ////////////////

const State& SystemRep::realizeTopology() const {
    if (!systemTopologyHasBeenRealized()) {
        defaultState.clear();
        defaultState.setNSubsystems(getNSubsystems());
        for (SubsystemId i(0); i<getNSubsystems(); ++i) 
            defaultState.initializeSubsystem(i, subsystems[i].getName(), 
                                                subsystems[i].getVersion());
        
        realizeTopologyImpl(defaultState); // defaultState is mutable
        systemTopologyRealized = true;
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

void SystemRep::realizeModel(State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage::Topology, 
        "System::realizeModel()");
    if (s.getSystemStage() < Stage::Model) {
        realizeModelImpl(s);    // take care of the Subsystems
        s.advanceSystemToStage(Stage::Model);
    }
}
void SystemRep::realizeInstance(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Instance).prev(), 
        "System::realizeInstance()");
    if (s.getSystemStage() < Stage::Instance) {
        realizeInstanceImpl(s);    // take care of the Subsystems
        s.advanceSystemToStage(Stage::Instance);
    }
}
void SystemRep::realizeTime(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Time).prev(), 
        "System::realizeTime()");
    if (s.getSystemStage() < Stage::Time) {
        realizeTimeImpl(s);    // take care of the Subsystems
        s.advanceSystemToStage(Stage::Time);
    }
}
void SystemRep::realizePosition(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Position).prev(), 
        "System::realizePosition()");
    if (s.getSystemStage() < Stage::Position) {
        realizePositionImpl(s);    // take care of the Subsystems
        s.advanceSystemToStage(Stage::Position);
    }
}
void SystemRep::realizeVelocity(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Velocity).prev(), 
        "System::realizeVelocity()");
    if (s.getSystemStage() < Stage::Velocity) {
        realizeVelocityImpl(s);    // take care of the Subsystems
        s.advanceSystemToStage(Stage::Velocity);
    }
}
void SystemRep::realizeDynamics(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Dynamics).prev(), 
        "System::realizeDynamics()");
    if (s.getSystemStage() < Stage::Dynamics) {
        realizeDynamicsImpl(s);    // take care of the Subsystems
        s.advanceSystemToStage(Stage::Dynamics);
    }
}
void SystemRep::realizeAcceleration(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Acceleration).prev(), 
        "System::realizeAcceleration()");
    if (s.getSystemStage() < Stage::Acceleration) {
        realizeAccelerationImpl(s);    // take care of the Subsystems
        s.advanceSystemToStage(Stage::Acceleration);
    }
}
void SystemRep::realizeReport(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getSystemStage(), Stage(Stage::Report).prev(), 
        "System::realizeReport()");
    if (s.getSystemStage() < Stage::Report) {
        realizeReportImpl(s);    // take care of the Subsystems
        s.advanceSystemToStage(Stage::Report);
    }
}

void SystemRep::realize(const State& s, Stage g) const {
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

void SystemRep::calcDecorativeGeometryAndAppend(const State& s, Stage stage, Array<DecorativeGeometry>& geom) const {
    assert(stage==Stage::Topology || s.getSystemStage() >= stage);
    for (int i=0; i<getNSubsystems(); ++i)
        subsystems[i].calcDecorativeGeometryAndAppend(s, stage, geom);
}


} // namespace SimTK

