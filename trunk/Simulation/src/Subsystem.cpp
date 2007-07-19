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

#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"
#include "SimTKcommon/internal/System.h"
#include "SimTKcommon/internal/Subsystem.h"

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

void Subsystem::librarySideConstruction(const String& name, const String& version) {
    rep = new SubsystemRep(name,version);
    rep->setMyHandle(*this);
}

void Subsystem::librarySideDestruction() {
    if (rep && rep->myHandle==this)
        delete rep; 
    rep=0;
}

void Subsystem::registerRealizeTopologyImpl(RealizeWritableStateImplLocator f) {
    updRep().realizeTopologyp = f;
}
void Subsystem::registerRealizeModelImpl(RealizeWritableStateImplLocator f) {
    updRep().realizeModelp = f;
}
void Subsystem::registerRealizeInstanceImpl(RealizeConstStateImplLocator f) {
    updRep().realizeInstancep = f;
}
void Subsystem::registerRealizeTimeImpl(RealizeConstStateImplLocator f) {
    updRep().realizeTimep = f;
}
void Subsystem::registerRealizePositionImpl(RealizeConstStateImplLocator f) {
    updRep().realizePositionp = f;
}
void Subsystem::registerRealizeVelocityImpl(RealizeConstStateImplLocator f) {
    updRep().realizeVelocityp = f;
}
void Subsystem::registerRealizeDynamicsImpl(RealizeConstStateImplLocator f) {
    updRep().realizeDynamicsp = f;
}
void Subsystem::registerRealizeAccelerationImpl(RealizeConstStateImplLocator f) {
    updRep().realizeAccelerationp = f;
}
void Subsystem::registerRealizeReportImpl(RealizeConstStateImplLocator f) {
    updRep().realizeReportp = f;
}

void Subsystem::registerCalcQUnitWeightsImpl(CalcUnitWeightsImplLocator f) {
    updRep().calcQUnitWeightsp = f;
}
void Subsystem::registerCalcUUnitWeightsImpl(CalcUnitWeightsImplLocator f) {
    updRep().calcUUnitWeightsp = f;
}
void Subsystem::registerCalcZUnitWeightsImpl(CalcUnitWeightsImplLocator f) {
    updRep().calcZUnitWeightsp = f;
}
void Subsystem::registerCalcQErrUnitTolerancesImpl(CalcUnitWeightsImplLocator f) {
    updRep().calcQErrUnitTolerancesp = f;
}
void Subsystem::registerCalcUErrUnitTolerancesImpl(CalcUnitWeightsImplLocator f) {
    updRep().calcUErrUnitTolerancesp = f;
}
void Subsystem::registerCalcDecorativeGeometryAndAppendImpl(CalcDecorativeGeometryAndAppendImplLocator f) {
    updRep().calcDecorativeGeometryAndAppendp = f;
}
void Subsystem::registerCloneImpl(CloneImplLocator f) {
    updRep().clonep = f;
}

const String& Subsystem::getName()    const {return getRep().getName();}
const String& Subsystem::getVersion() const {return getRep().getVersion();}


bool Subsystem::isInSystem() const {return getRep().isInSystem();}
bool Subsystem::isInSameSystem(const Subsystem& otherSubsystem) const {
	return getRep().isInSameSystem(otherSubsystem);
}
const System& Subsystem::getSystem() const {return getRep().getSystem();}
System&       Subsystem::updSystem()	   {return updRep().updSystem();}
SubsystemId   Subsystem::getMySubsystemId() const {return getRep().getMySubsystemId();}

int Subsystem::allocateQ(State& s, const Vector& qInit) const {
    return s.allocateQ(getRep().getMySubsystemId(), qInit);
}

int Subsystem::allocateU(State& s, const Vector& uInit) const {
    return s.allocateU(getRep().getMySubsystemId(), uInit);
}

int Subsystem::allocateZ(State& s, const Vector& zInit) const {
    return s.allocateZ(getRep().getMySubsystemId(), zInit);
}

int Subsystem::allocateQErr(State& s, int nqerr) const {
    return s.allocateQErr(getRep().getMySubsystemId(), nqerr);
}

int Subsystem::allocateUErr(State& s, int nuerr) const {
    return s.allocateUErr(getRep().getMySubsystemId(), nuerr);
}

int Subsystem::allocateUDotErr(State& s, int nudoterr) const {
    return s.allocateUDotErr(getRep().getMySubsystemId(), nudoterr);
}

int Subsystem::allocateDiscreteVariable(State& s, Stage g, AbstractValue* v) const {
    return s.allocateDiscreteVariable(getRep().getMySubsystemId(), g, v);
}

int Subsystem::allocateCacheEntry(State& s, Stage g, AbstractValue* v) const {
    return s.allocateCacheEntry(getRep().getMySubsystemId(), g, v);
}

void Subsystem::advanceToStage(const State& s, Stage g) const {
    s.advanceSubsystemToStage(getRep().getMySubsystemId(), g);
}

Stage Subsystem::getStage(const State& s) const {
    return s.getSubsystemStage(getRep().getMySubsystemId());
}
const AbstractValue& Subsystem::getDiscreteVariable(const State& s, int index) const {
    return s.getDiscreteVariable(getRep().getMySubsystemId(), index);
}

AbstractValue& Subsystem::updDiscreteVariable(State& s, int index) const {
    return s.updDiscreteVariable(getRep().getMySubsystemId(), index);
}

const AbstractValue& Subsystem::getCacheEntry(const State& s, int index) const {
    return s.getCacheEntry(getRep().getMySubsystemId(), index);
}

AbstractValue& Subsystem::updCacheEntry(const State& s, int index) const {
    return s.updCacheEntry(getRep().getMySubsystemId(), index);
}

const Vector& Subsystem::getQ(const State& s) const {return s.getQ(getRep().getMySubsystemId());}
const Vector& Subsystem::getU(const State& s) const {return s.getU(getRep().getMySubsystemId());}
const Vector& Subsystem::getZ(const State& s) const {return s.getZ(getRep().getMySubsystemId());}

Vector& Subsystem::updQ(State& s) const {return s.updQ(getRep().getMySubsystemId());}
Vector& Subsystem::updU(State& s) const {return s.updU(getRep().getMySubsystemId());}
Vector& Subsystem::updZ(State& s) const {return s.updZ(getRep().getMySubsystemId());}

const Vector& Subsystem::getQDot   (const State& s) const {return s.getQDot(getRep().getMySubsystemId());}
const Vector& Subsystem::getUDot   (const State& s) const {return s.getUDot(getRep().getMySubsystemId());}
const Vector& Subsystem::getZDot   (const State& s) const {return s.getZDot(getRep().getMySubsystemId());}
const Vector& Subsystem::getQDotDot(const State& s) const {return s.getQDotDot(getRep().getMySubsystemId());}

Vector& Subsystem::updQDot   (const State& s) const {return s.updQDot(getRep().getMySubsystemId());}
Vector& Subsystem::updUDot   (const State& s) const {return s.updUDot(getRep().getMySubsystemId());}
Vector& Subsystem::updZDot   (const State& s) const {return s.updZDot(getRep().getMySubsystemId());}
Vector& Subsystem::updQDotDot(const State& s) const {return s.updQDotDot(getRep().getMySubsystemId());}

const Vector& Subsystem::getQErr(const State& s) const {return s.getQErr(getRep().getMySubsystemId());}
const Vector& Subsystem::getUErr(const State& s) const {return s.getUErr(getRep().getMySubsystemId());}
const Vector& Subsystem::getUDotErr(const State& s) const {return s.getUDotErr(getRep().getMySubsystemId());}

Vector& Subsystem::updQErr(const State& s) const {return s.updQErr(getRep().getMySubsystemId());}
Vector& Subsystem::updUErr(const State& s) const {return s.updUErr(getRep().getMySubsystemId());}
Vector& Subsystem::updUDotErr(const State& s) const {return s.updUDotErr(getRep().getMySubsystemId());}

int Subsystem::getQStart(const State& s) const {return s.getQStart(getRep().getMySubsystemId());}
int Subsystem::getNQ(const State& s)     const {return s.getNQ(getRep().getMySubsystemId());}

int Subsystem::getUStart(const State& s) const {return s.getUStart(getRep().getMySubsystemId());}
int Subsystem::getNU(const State& s)     const {return s.getNU(getRep().getMySubsystemId());}

int Subsystem::getZStart(const State& s) const {return s.getZStart(getRep().getMySubsystemId());}
int Subsystem::getNZ(const State& s)     const {return s.getNZ(getRep().getMySubsystemId());}

int Subsystem::getQErrStart(const State& s) const {return s.getQErrStart(getRep().getMySubsystemId());}
int Subsystem::getNQErr(const State& s)     const {return s.getNQErr(getRep().getMySubsystemId());}

int Subsystem::getUErrStart(const State& s) const {return s.getUErrStart(getRep().getMySubsystemId());}
int Subsystem::getNUErr(const State& s)     const {return s.getNUErr(getRep().getMySubsystemId());}

int Subsystem::getUDotErrStart(const State& s) const {return s.getUDotErrStart(getRep().getMySubsystemId());}
int Subsystem::getNUDotErr(const State& s)     const {return s.getNUDotErr(getRep().getMySubsystemId());}


void Subsystem::invalidateSubsystemTopologyCache() const {
    getRep().invalidateSubsystemTopologyCache();
}

bool Subsystem::subsystemTopologyCacheHasBeenRealized() const {
    return getRep().subsystemTopologyCacheHasBeenRealized();
}

    // wrappers for Subsystem virtuals

void Subsystem::realizeSubsystemTopology(State& s) const {
    SimTK_STAGECHECK_EQ_ALWAYS(getStage(s), Stage::Empty, 
        "Subsystem::realizeSubsystemTopology()");
    getRep().realizeTopologyp(*this,s);
    getRep().subsystemTopologyRealized = true; // mark the subsystem itself (mutable)
    advanceToStage(s, Stage::Topology);  // mark the State as well
}
void Subsystem::realizeSubsystemModel(State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage::Topology, 
        "Subsystem::realizeSubsystemModel()");
    if (getStage(s) < Stage::Model) {
        getRep().realizeModelp(*this,s);
        advanceToStage(s, Stage::Model);
    }
}
void Subsystem::realizeSubsystemInstance(const State& s) const { 
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Instance).prev(), 
        "Subsystem::realizeSubsystemInstance()");
    if (getStage(s) < Stage::Instance) {
        getRep().realizeInstancep(*this,s);
        advanceToStage(s, Stage::Instance);
    }
}
void Subsystem::realizeSubsystemTime(const State& s) const { 
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Time).prev(), 
        "Subsystem::realizeTime()");
    if (getStage(s) < Stage::Time) {
        getRep().realizeTimep(*this,s);
        advanceToStage(s, Stage::Time);
    }
}
void Subsystem::realizeSubsystemPosition(const State& s) const { 
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Position).prev(), 
        "Subsystem::realizeSubsystemPosition()");
    if (getStage(s) < Stage::Position) {
        getRep().realizePositionp(*this,s);
        advanceToStage(s, Stage::Position);
    }
}
void Subsystem::realizeSubsystemVelocity(const State& s) const { 
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Velocity).prev(), 
        "Subsystem::realizeSubsystemVelocity()");
    if (getStage(s) < Stage::Velocity) {
        getRep().realizeVelocityp(*this,s);
        advanceToStage(s, Stage::Velocity);
    }
}
void Subsystem::realizeSubsystemDynamics(const State& s) const { 
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Dynamics).prev(), 
        "Subsystem::realizeSubsystemDynamics()");
    if (getStage(s) < Stage::Dynamics) {
        getRep().realizeDynamicsp(*this,s);
        advanceToStage(s, Stage::Dynamics);
    }
}
void Subsystem::realizeSubsystemAcceleration(const State& s) const { 
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Acceleration).prev(), 
        "Subsystem::realizeSubsystemAcceleration()");
    if (getStage(s) < Stage::Acceleration) {
        getRep().realizeAccelerationp(*this,s);
        advanceToStage(s, Stage::Acceleration);
    }
}
void Subsystem::realizeSubsystemReport(const State& s) const { 
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Report).prev(), 
        "Subsystem::realizeSubsystemReport()");
    if (getStage(s) < Stage::Report) {
        getRep().realizeReportp(*this,s);
        advanceToStage(s, Stage::Report);
    }
}


void Subsystem::calcQUnitWeights(const State& s, Vector& weights) const {
    getRep().calcQUnitWeightsp(*this,s,weights);
}
void Subsystem::calcUUnitWeights(const State& s, Vector& weights) const {
    getRep().calcUUnitWeightsp(*this,s,weights);
}
void Subsystem::calcZUnitWeights(const State& s, Vector& weights) const {
    getRep().calcZUnitWeightsp(*this,s,weights);
}
void Subsystem::calcQErrUnitTolerances(const State& s, Vector& tolerances) const {
    getRep().calcQErrUnitTolerancesp(*this,s,tolerances);
}
void Subsystem::calcUErrUnitTolerances(const State& s, Vector& tolerances) const {
    getRep().calcUErrUnitTolerancesp(*this,s,tolerances);
}

void Subsystem::calcDecorativeGeometryAndAppend(const State& s, Stage stage, Array<DecorativeGeometry>& geom) const {
    getRep().calcDecorativeGeometryAndAppendp(*this,s,stage,geom);
}

Subsystem* Subsystem::clone() const {
    return getRep().clonep(*this);
}


    // default implementations for Subsystem /*virtual*/s
/*virtual*/ int Subsystem::realizeSubsystemTopologyImpl(State& s) const {
    return 0;
}
/*virtual*/ int Subsystem::realizeSubsystemModelImpl(State& s) const {
    return 0; 
}
/*virtual*/ int Subsystem::realizeSubsystemInstanceImpl(const State& s) const {
    return 0; 
}
/*virtual*/ int Subsystem::realizeSubsystemTimeImpl(const State& s) const {
    return 0; 
}
/*virtual*/ int Subsystem::realizeSubsystemPositionImpl(const State& s) const {
    return 0; 
}
/*virtual*/ int Subsystem::realizeSubsystemVelocityImpl(const State& s) const {
    return 0; 
}
/*virtual*/ int Subsystem::realizeSubsystemDynamicsImpl(const State& s) const {
    return 0; 
}
/*virtual*/ int Subsystem::realizeSubsystemAccelerationImpl(const State& s) const {
    return 0; 
}
/*virtual*/ int Subsystem::realizeSubsystemReportImpl(const State& s) const {
    return 0; 
}

/*virtual*/ int Subsystem::calcQUnitWeightsImpl(const State& s, Vector& weights) const {
    weights.resize(getNQ(s));
    weights = 1; // default says everyone's opinion is just as valid
    return 0;
}
/*virtual*/ int Subsystem::calcUUnitWeightsImpl(const State& s, Vector& weights) const {
    weights.resize(getNU(s));
    weights = 1;
    return 0;
}
/*virtual*/ int Subsystem::calcZUnitWeightsImpl(const State& s, Vector& weights) const {
    weights.resize(getNZ(s));
    weights = 1;
    return 0;
}
/*virtual*/ int Subsystem::calcQErrUnitTolerancesImpl(const State& s, Vector& tolerances) const {
    tolerances.resize(getNQErr(s));
    tolerances = 1;
    return 0;
}
/*virtual*/ int Subsystem::calcUErrUnitTolerancesImpl(const State& s, Vector& tolerances) const {
    tolerances.resize(getNUErr(s));
    tolerances = 1;
    return 0;
}
/*virtual*/ int Subsystem::calcDecorativeGeometryAndAppendImpl
                                (const State&, Stage, Array<DecorativeGeometry>&) const
{
    return 0;
}

/*virtual*/ Subsystem* Subsystem::cloneImpl() const {
    return new Subsystem(*this);
}

void Subsystem::adoptPrivateImplementation(Subsystem::PrivateImplementation* p,
                                           Subsystem::ClonePrivateImplementation clone,
                                           Subsystem::DestructPrivateImplementation destruct)
{
    updRep().adoptPrivateImplementation(p,clone,destruct);
}

const Subsystem::PrivateImplementation& Subsystem::getPrivateImplementation() const {
    return getRep().getPrivateImplementation();
}

Subsystem::PrivateImplementation& Subsystem::updPrivateImplementation() {
    return updRep().updPrivateImplementation();
}
    ///////////////////
    // SUBSYSTEM REP //
    ///////////////////

void SubsystemRep::invalidateSubsystemTopologyCache() const {
    subsystemTopologyRealized = false;
    if (isInSystem()) 
        getSystem().getRep().invalidateSystemTopologyCache();
}

    //////////////////////////////
    // DEFAULT SYSTEM SUBSYSTEM //
    //////////////////////////////

} // namespace SimTK

