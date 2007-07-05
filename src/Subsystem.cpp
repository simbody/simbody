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
 * Implementation of Subsystem, SubsystemRep and DefaultSystemSubsystem.
 */

#include "simbody/internal/common.h"
#include "simbody/internal/System.h"
#include "simbody/internal/Subsystem.h"

#include "SystemRep.h"
#include "SubsystemRep.h"

#include <cassert>

namespace SimTK {

    ///////////////
    // SUBSYSTEM //
    ///////////////

bool Subsystem::isEmptyHandle() const {return rep==0;}
bool Subsystem::isOwnerHandle() const {return rep==0 || rep->myHandle==this;}
bool Subsystem::isSameSubsystem(const Subsystem& otherSubsystem) const {
    return rep && (rep==otherSubsystem.rep);
}


Subsystem::~Subsystem() {
    if (isOwnerHandle()) delete rep; 
    rep=0;
}

Subsystem::Subsystem(const Subsystem& src) : rep(0) {
    if (src.rep) {
		assert(!src.rep->isInSystem()); // TODO
        rep = src.rep->clone();	// create a new object
        rep->setMyHandle(*this);
    }
}

Subsystem& Subsystem::operator=(const Subsystem& src) {
    if (&src != this) {
        if (isOwnerHandle()) delete rep; 
        rep=0;
        if (src.rep) {
			assert(!src.rep->isInSystem()); // TODO
			rep = src.rep->clone();	// create a new object
			rep->setMyHandle(*this);
        }
    }
    return *this;
}



const String& Subsystem::getName()    const {return getRep().getName();}
const String& Subsystem::getVersion() const {return getRep().getVersion();}

//void Subsystem::realize(const State& s, Stage g) const {
//    getRep().realize(s,g);
//}

void Subsystem::calcDecorativeGeometryAndAppend(const State& s, Stage stage, Array<DecorativeGeometry>& geom) const {
    getRep().calcDecorativeGeometryAndAppend(s,stage,geom);
}

void Subsystem::endConstruction() {updRep().endConstruction();}

bool Subsystem::isInSystem() const {return getRep().isInSystem();}
bool Subsystem::isInSameSystem(const Subsystem& otherSubsystem) const {
	return getRep().isInSameSystem(otherSubsystem);
}
const System& Subsystem::getSystem() const {return getRep().getSystem();}
System&       Subsystem::updSystem()	   {return updRep().updSystem();}
SubsystemId   Subsystem::getMySubsystemId() const {return getRep().getMySubsystemId();}

const Vector& Subsystem::getQ(const State& s) const {return getRep().getQ(s);}
const Vector& Subsystem::getU(const State& s) const {return getRep().getU(s);}
const Vector& Subsystem::getZ(const State& s) const {return getRep().getZ(s);}

const Vector& Subsystem::getQErr(const State& s)    const {return getRep().getQErr(s);}
const Vector& Subsystem::getUErr(const State& s)    const {return getRep().getUErr(s);}
const Vector& Subsystem::getUDotErr(const State& s) const {return getRep().getUDotErr(s);}

void Subsystem::calcQUnitWeights(const State& s, Vector& weights) const {
    getRep().calcQUnitWeights(s,weights);
}
void Subsystem::calcUUnitWeights(const State& s, Vector& weights) const {
    getRep().calcUUnitWeights(s,weights);
}
void Subsystem::calcZUnitWeights(const State& s, Vector& weights) const {
    getRep().calcZUnitWeights(s,weights);
}
void Subsystem::calcQErrUnitTolerances(const State& s, Vector& tolerances) const {
    getRep().calcQErrUnitTolerances(s,tolerances);
}
void Subsystem::calcUErrUnitTolerances(const State& s, Vector& tolerances) const {
    getRep().calcUErrUnitTolerances(s,tolerances);
}

    ///////////////////
    // SUBSYSTEM REP //
    ///////////////////

void SubsystemRep::realizeSubsystemTopology(State& s) const {
    SimTK_STAGECHECK_EQ_ALWAYS(getStage(s), Stage::Empty, 
        "Subsystem::realizeTopology()");
    realizeSubsystemTopologyImpl(s);
    subsystemTopologyRealized = true;   // mark the subsystem itself (mutable)
    advanceToStage(s, Stage::Topology); // mark the State as well
}
void SubsystemRep::realizeSubsystemModel(State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage::Topology, 
        "Subsystem::realizeModel()");
    if (getStage(s) < Stage::Model) {
        realizeSubsystemModelImpl(s);
        advanceToStage(s, Stage::Model);
    }
}
void SubsystemRep::realizeSubsystemInstance(const State& s) const { 
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Instance).prev(), 
        "Subsystem::realizeInstance()");
    if (getStage(s) < Stage::Instance) {
        realizeSubsystemInstanceImpl(s);
        advanceToStage(s, Stage::Instance);
    }
}
void SubsystemRep::realizeSubsystemTime(const State& s) const { 
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Time).prev(), 
        "Subsystem::realizeTime()");
    if (getStage(s) < Stage::Time) {
        realizeSubsystemTimeImpl(s);
        advanceToStage(s, Stage::Time);
    }
}
void SubsystemRep::realizeSubsystemPosition(const State& s) const { 
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Position).prev(), 
        "Subsystem::realizePosition()");
    if (getStage(s) < Stage::Position) {
        realizeSubsystemPositionImpl(s);
        advanceToStage(s, Stage::Position);
    }
}
void SubsystemRep::realizeSubsystemVelocity(const State& s) const { 
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Velocity).prev(), 
        "Subsystem::realizeVelocity()");
    if (getStage(s) < Stage::Velocity) {
        realizeSubsystemVelocityImpl(s);
        advanceToStage(s, Stage::Velocity);
    }
}
void SubsystemRep::realizeSubsystemDynamics(const State& s) const { 
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Dynamics).prev(), 
        "Subsystem::realizeDynamics()");
    if (getStage(s) < Stage::Dynamics) {
        realizeSubsystemDynamicsImpl(s);
        advanceToStage(s, Stage::Dynamics);
    }
}
void SubsystemRep::realizeSubsystemAcceleration(const State& s) const { 
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Acceleration).prev(), 
        "Subsystem::realizeAcceleration()");
    if (getStage(s) < Stage::Acceleration) {
        realizeSubsystemAccelerationImpl(s);
        advanceToStage(s, Stage::Acceleration);
    }
}
void SubsystemRep::realizeSubsystemReport(const State& s) const { 
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Report).prev(), 
        "Subsystem::realizeReport()");
    if (getStage(s) < Stage::Report) {
        realizeSubsystemReportImpl(s);
        advanceToStage(s, Stage::Report);
    }
}

void SubsystemRep::invalidateSubsystemTopologyCache() {
    subsystemTopologyRealized = false;
    if (isInSystem()) 
        updSystem().updRep().invalidateSystemTopologyCache();
}

    //////////////////////////////
    // DEFAULT SYSTEM SUBSYSTEM //
    //////////////////////////////


DefaultSystemSubsystem::DefaultSystemSubsystem() {
    rep = new DefaultSystemSubsystemRep();
    rep->setMyHandle(*this);
}

DefaultSystemSubsystem::DefaultSystemSubsystem(const String& sysName, const String& sysVersion) {
    rep = new DefaultSystemSubsystemRep();
    rep->setMyHandle(*this);
}

/*static*/ bool
DefaultSystemSubsystem::isInstanceOf(const Subsystem& s) {
    return DefaultSystemSubsystemRep::isA(s.getRep());
}

/*static*/ const DefaultSystemSubsystem&
DefaultSystemSubsystem::downcast(const Subsystem& s) {
    assert(DefaultSystemSubsystemRep::isA(s.getRep()));
    return reinterpret_cast<const DefaultSystemSubsystem&>(s);
}

/*static*/ DefaultSystemSubsystem&
DefaultSystemSubsystem::updDowncast(Subsystem& s) {
    assert(DefaultSystemSubsystemRep::isA(s.getRep()));
    return reinterpret_cast<DefaultSystemSubsystem&>(s);
}

DefaultSystemSubsystemRep& DefaultSystemSubsystem::updRep() {
    return DefaultSystemSubsystemRep::downcast(*rep);
}

const DefaultSystemSubsystemRep& DefaultSystemSubsystem::getRep() const {
    return DefaultSystemSubsystemRep::downcast(*rep);
}

} // namespace SimTK

