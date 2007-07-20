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

bool Subsystem::isEmptyHandle() const {return guts==0;}
bool Subsystem::isOwnerHandle() const {return guts==0 || &guts->getOwnerSubsystemHandle()==this;}
bool Subsystem::isSameSubsystem(const Subsystem& otherSubsystem) const {
    return guts && (guts==otherSubsystem.guts);
}


Subsystem::Subsystem(const Subsystem& src) : guts(0) {
    if (src.guts) {
        guts = src.guts->clone();
        guts->setOwnerSubsystemHandle(*this);
    }
}

Subsystem& Subsystem::operator=(const Subsystem& src) {
    if (!isSameSubsystem(src)) {
        if (isOwnerHandle()) delete guts; 
        guts=0;
        if (src.guts) {
            guts = src.guts->clone();
            guts->setOwnerSubsystemHandle(*this);
        }
    }
    return *this;
}

Subsystem::~Subsystem() {
    //TODO: delete should probably be called from library side VFT
    if (guts && isOwnerHandle())
        delete guts;
    guts=0;
}

void Subsystem::adoptSubsystemGuts(Subsystem::Guts* g) {
    SimTK_ASSERT_ALWAYS(g, "Subsystem::adoptSubsystemGuts(): can't adopt null Guts");
    SimTK_ASSERT_ALWAYS(!guts,
        "Subsystem::adoptSubsystemGuts(): this Subsystem handle is already in use");
    guts = g;
    guts->setOwnerSubsystemHandle(*this);
}

void Subsystem::setSystem(System& sys, SubsystemId id) {
    updSubsystemGuts().setSystem(sys,id);
}

const String& Subsystem::getName()    const {return getSubsystemGuts().getName();}
const String& Subsystem::getVersion() const {return getSubsystemGuts().getVersion();}

bool Subsystem::subsystemTopologyHasBeenRealized() const {
    return getSubsystemGuts().subsystemTopologyHasBeenRealized();
}

void Subsystem::invalidateSubsystemTopologyCache() const {
    return getSubsystemGuts().invalidateSubsystemTopologyCache(); // mutable
}

bool Subsystem::isInSystem() const {return getSubsystemGuts().isInSystem();}
bool Subsystem::isInSameSystem(const Subsystem& otherSubsystem) const {
    return getSubsystemGuts().isInSameSystem(otherSubsystem);
}

const System& Subsystem::getSystem() const {return getSubsystemGuts().getSystem();}
System&       Subsystem::updSystem()       {return updSubsystemGuts().updSystem();}

SubsystemId Subsystem::getMySubsystemId() const {
    return getSubsystemGuts().getMySubsystemId();
}

const Vector& Subsystem::getQ(const State& s) const {return getSubsystemGuts().getQ(s);}
const Vector& Subsystem::getU(const State& s) const {return getSubsystemGuts().getU(s);}
const Vector& Subsystem::getZ(const State& s) const {return getSubsystemGuts().getZ(s);}
const Vector& Subsystem::getQDot(const State& s) const {return getSubsystemGuts().getQDot(s);}
const Vector& Subsystem::getUDot(const State& s) const {return getSubsystemGuts().getUDot(s);}
const Vector& Subsystem::getZDot(const State& s) const {return getSubsystemGuts().getZDot(s);}
const Vector& Subsystem::getQDotDot(const State& s) const {return getSubsystemGuts().getQDotDot(s);}
const Vector& Subsystem::getQErr(const State& s) const {return getSubsystemGuts().getQErr(s);}
const Vector& Subsystem::getUErr(const State& s) const {return getSubsystemGuts().getUErr(s);}
const Vector& Subsystem::getUDotErr(const State& s) const {return getSubsystemGuts().getUDotErr(s);}

Vector& Subsystem::updQ(State& s) const {return getSubsystemGuts().updQ(s);}
Vector& Subsystem::updU(State& s) const {return getSubsystemGuts().updU(s);}
Vector& Subsystem::updZ(State& s) const {return getSubsystemGuts().updZ(s);}

Vector& Subsystem::updQDot(const State& s) const {return getSubsystemGuts().updQDot(s);}
Vector& Subsystem::updUDot(const State& s) const {return getSubsystemGuts().updUDot(s);}
Vector& Subsystem::updZDot(const State& s) const {return getSubsystemGuts().updZDot(s);}
Vector& Subsystem::updQDotDot(const State& s) const {return getSubsystemGuts().updQDotDot(s);}
Vector& Subsystem::updQErr(const State& s) const {return getSubsystemGuts().updQErr(s);}
Vector& Subsystem::updUErr(const State& s) const {return getSubsystemGuts().updUErr(s);}
Vector& Subsystem::updUDotErr(const State& s) const {return getSubsystemGuts().updUDotErr(s);}

int Subsystem::getQStart      (const State& s) const {return getSubsystemGuts().getQStart(s);}
int Subsystem::getNQ          (const State& s) const {return getSubsystemGuts().getNQ(s);}
int Subsystem::getUStart      (const State& s) const {return getSubsystemGuts().getUStart(s);}
int Subsystem::getNU          (const State& s) const {return getSubsystemGuts().getNU(s);}
int Subsystem::getZStart      (const State& s) const {return getSubsystemGuts().getZStart(s);}
int Subsystem::getNZ          (const State& s) const {return getSubsystemGuts().getNZ(s);}
int Subsystem::getQErrStart   (const State& s) const {return getSubsystemGuts().getQErrStart(s);}
int Subsystem::getNQErr       (const State& s) const {return getSubsystemGuts().getNQErr(s);}
int Subsystem::getUErrStart   (const State& s) const {return getSubsystemGuts().getUErrStart(s);}
int Subsystem::getNUErr       (const State& s) const {return getSubsystemGuts().getNUErr(s);}
int Subsystem::getUDotErrStart(const State& s) const {return getSubsystemGuts().getUDotErrStart(s);}
int Subsystem::getNUDotErr    (const State& s) const {return getSubsystemGuts().getNUDotErr(s);}

    /////////////////////
    // SUBSYSTEM::GUTS //
    /////////////////////

// Default constructor is inline, but calls librarySideConstuction() here.
void Subsystem::Guts::librarySideConstruction(const String& name, const String& version) {
    rep = new GutsRep(name,version);
    // note that the GutsRep object currently has no owner handle
}

// Destructor is inline, but calls librarySideDestruction() here.
void Subsystem::Guts::librarySideDestruction() {
    delete rep; 
    rep=0;
}


// Copy constructor
Subsystem::Guts::Guts(const Guts& src) : rep(0) {
    if (src.rep) {
        rep = new Subsystem::Guts::GutsRep(*src.rep);
        // note that the GutsRep object currently has no owner handle
    }
}

// Copy assignment is suppressed
    

const Subsystem& Subsystem::Guts::getOwnerSubsystemHandle() const {
    assert(rep->myHandle);
    return *rep->myHandle;
}

Subsystem& Subsystem::Guts::updOwnerSubsystemHandle() {
    assert(rep->myHandle);
    return *rep->myHandle;
}

void Subsystem::Guts::setOwnerSubsystemHandle(Subsystem& sys) {
    // might be the first owner or a replacement
    rep->myHandle = &sys;
}

bool Subsystem::Guts::hasOwnerSubsystemHandle() const {
    return rep->myHandle != 0;
}

void Subsystem::Guts::setSystem(System& sys, SubsystemId id) {
    updRep().setSystem(sys,id);
}

const String& Subsystem::Guts::getName()    const {return getRep().getName();}
const String& Subsystem::Guts::getVersion() const {return getRep().getVersion();}


void Subsystem::Guts::registerCloneImpl(CloneImplLocator f) {
    updRep().clonep = f;
}

void Subsystem::Guts::registerRealizeTopologyImpl(RealizeWritableStateImplLocator f) {
    updRep().realizeTopologyp = f;
}
void Subsystem::Guts::registerRealizeModelImpl(RealizeWritableStateImplLocator f) {
    updRep().realizeModelp = f;
}
void Subsystem::Guts::registerRealizeInstanceImpl(RealizeConstStateImplLocator f) {
    updRep().realizeInstancep = f;
}
void Subsystem::Guts::registerRealizeTimeImpl(RealizeConstStateImplLocator f) {
    updRep().realizeTimep = f;
}
void Subsystem::Guts::registerRealizePositionImpl(RealizeConstStateImplLocator f) {
    updRep().realizePositionp = f;
}
void Subsystem::Guts::registerRealizeVelocityImpl(RealizeConstStateImplLocator f) {
    updRep().realizeVelocityp = f;
}
void Subsystem::Guts::registerRealizeDynamicsImpl(RealizeConstStateImplLocator f) {
    updRep().realizeDynamicsp = f;
}
void Subsystem::Guts::registerRealizeAccelerationImpl(RealizeConstStateImplLocator f) {
    updRep().realizeAccelerationp = f;
}
void Subsystem::Guts::registerRealizeReportImpl(RealizeConstStateImplLocator f) {
    updRep().realizeReportp = f;
}

void Subsystem::Guts::registerCalcQUnitWeightsImpl(CalcUnitWeightsImplLocator f) {
    updRep().calcQUnitWeightsp = f;
}
void Subsystem::Guts::registerCalcUUnitWeightsImpl(CalcUnitWeightsImplLocator f) {
    updRep().calcUUnitWeightsp = f;
}
void Subsystem::Guts::registerCalcZUnitWeightsImpl(CalcUnitWeightsImplLocator f) {
    updRep().calcZUnitWeightsp = f;
}
void Subsystem::Guts::registerCalcQErrUnitTolerancesImpl(CalcUnitWeightsImplLocator f) {
    updRep().calcQErrUnitTolerancesp = f;
}
void Subsystem::Guts::registerCalcUErrUnitTolerancesImpl(CalcUnitWeightsImplLocator f) {
    updRep().calcUErrUnitTolerancesp = f;
}
void Subsystem::Guts::registerCalcDecorativeGeometryAndAppendImpl(CalcDecorativeGeometryAndAppendImplLocator f) {
    updRep().calcDecorativeGeometryAndAppendp = f;
}

bool Subsystem::Guts::isInSystem() const {return getRep().isInSystem();}
bool Subsystem::Guts::isInSameSystem(const Subsystem& otherSubsystem) const {
	return getRep().isInSameSystem(otherSubsystem);
}
const System& Subsystem::Guts::getSystem() const {return getRep().getSystem();}
System&       Subsystem::Guts::updSystem()	     {return updRep().updSystem();}
SubsystemId   Subsystem::Guts::getMySubsystemId() const {return getRep().getMySubsystemId();}

int Subsystem::Guts::allocateQ(State& s, const Vector& qInit) const {
    return s.allocateQ(getRep().getMySubsystemId(), qInit);
}

int Subsystem::Guts::allocateU(State& s, const Vector& uInit) const {
    return s.allocateU(getRep().getMySubsystemId(), uInit);
}

int Subsystem::Guts::allocateZ(State& s, const Vector& zInit) const {
    return s.allocateZ(getRep().getMySubsystemId(), zInit);
}

int Subsystem::Guts::allocateQErr(State& s, int nqerr) const {
    return s.allocateQErr(getRep().getMySubsystemId(), nqerr);
}

int Subsystem::Guts::allocateUErr(State& s, int nuerr) const {
    return s.allocateUErr(getRep().getMySubsystemId(), nuerr);
}

int Subsystem::Guts::allocateUDotErr(State& s, int nudoterr) const {
    return s.allocateUDotErr(getRep().getMySubsystemId(), nudoterr);
}

int Subsystem::Guts::allocateDiscreteVariable(State& s, Stage g, AbstractValue* v) const {
    return s.allocateDiscreteVariable(getRep().getMySubsystemId(), g, v);
}

int Subsystem::Guts::allocateCacheEntry(State& s, Stage g, AbstractValue* v) const {
    return s.allocateCacheEntry(getRep().getMySubsystemId(), g, v);
}

void Subsystem::Guts::advanceToStage(const State& s, Stage g) const {
    s.advanceSubsystemToStage(getRep().getMySubsystemId(), g);
}

Stage Subsystem::Guts::getStage(const State& s) const {
    return s.getSubsystemStage(getRep().getMySubsystemId());
}
const AbstractValue& Subsystem::Guts::getDiscreteVariable(const State& s, int index) const {
    return s.getDiscreteVariable(getRep().getMySubsystemId(), index);
}

AbstractValue& Subsystem::Guts::updDiscreteVariable(State& s, int index) const {
    return s.updDiscreteVariable(getRep().getMySubsystemId(), index);
}

const AbstractValue& Subsystem::Guts::getCacheEntry(const State& s, int index) const {
    return s.getCacheEntry(getRep().getMySubsystemId(), index);
}

AbstractValue& Subsystem::Guts::updCacheEntry(const State& s, int index) const {
    return s.updCacheEntry(getRep().getMySubsystemId(), index);
}

const Vector& Subsystem::Guts::getQ(const State& s) const {return s.getQ(getRep().getMySubsystemId());}
const Vector& Subsystem::Guts::getU(const State& s) const {return s.getU(getRep().getMySubsystemId());}
const Vector& Subsystem::Guts::getZ(const State& s) const {return s.getZ(getRep().getMySubsystemId());}

Vector& Subsystem::Guts::updQ(State& s) const {return s.updQ(getRep().getMySubsystemId());}
Vector& Subsystem::Guts::updU(State& s) const {return s.updU(getRep().getMySubsystemId());}
Vector& Subsystem::Guts::updZ(State& s) const {return s.updZ(getRep().getMySubsystemId());}

const Vector& Subsystem::Guts::getQDot   (const State& s) const {return s.getQDot(getRep().getMySubsystemId());}
const Vector& Subsystem::Guts::getUDot   (const State& s) const {return s.getUDot(getRep().getMySubsystemId());}
const Vector& Subsystem::Guts::getZDot   (const State& s) const {return s.getZDot(getRep().getMySubsystemId());}
const Vector& Subsystem::Guts::getQDotDot(const State& s) const {return s.getQDotDot(getRep().getMySubsystemId());}

Vector& Subsystem::Guts::updQDot   (const State& s) const {return s.updQDot(getRep().getMySubsystemId());}
Vector& Subsystem::Guts::updUDot   (const State& s) const {return s.updUDot(getRep().getMySubsystemId());}
Vector& Subsystem::Guts::updZDot   (const State& s) const {return s.updZDot(getRep().getMySubsystemId());}
Vector& Subsystem::Guts::updQDotDot(const State& s) const {return s.updQDotDot(getRep().getMySubsystemId());}

const Vector& Subsystem::Guts::getQErr(const State& s) const {return s.getQErr(getRep().getMySubsystemId());}
const Vector& Subsystem::Guts::getUErr(const State& s) const {return s.getUErr(getRep().getMySubsystemId());}
const Vector& Subsystem::Guts::getUDotErr(const State& s) const {return s.getUDotErr(getRep().getMySubsystemId());}

Vector& Subsystem::Guts::updQErr(const State& s) const {return s.updQErr(getRep().getMySubsystemId());}
Vector& Subsystem::Guts::updUErr(const State& s) const {return s.updUErr(getRep().getMySubsystemId());}
Vector& Subsystem::Guts::updUDotErr(const State& s) const {return s.updUDotErr(getRep().getMySubsystemId());}

int Subsystem::Guts::getQStart(const State& s) const {return s.getQStart(getRep().getMySubsystemId());}
int Subsystem::Guts::getNQ(const State& s)     const {return s.getNQ(getRep().getMySubsystemId());}

int Subsystem::Guts::getUStart(const State& s) const {return s.getUStart(getRep().getMySubsystemId());}
int Subsystem::Guts::getNU(const State& s)     const {return s.getNU(getRep().getMySubsystemId());}

int Subsystem::Guts::getZStart(const State& s) const {return s.getZStart(getRep().getMySubsystemId());}
int Subsystem::Guts::getNZ(const State& s)     const {return s.getNZ(getRep().getMySubsystemId());}

int Subsystem::Guts::getQErrStart(const State& s) const {return s.getQErrStart(getRep().getMySubsystemId());}
int Subsystem::Guts::getNQErr(const State& s)     const {return s.getNQErr(getRep().getMySubsystemId());}

int Subsystem::Guts::getUErrStart(const State& s) const {return s.getUErrStart(getRep().getMySubsystemId());}
int Subsystem::Guts::getNUErr(const State& s)     const {return s.getNUErr(getRep().getMySubsystemId());}

int Subsystem::Guts::getUDotErrStart(const State& s) const {return s.getUDotErrStart(getRep().getMySubsystemId());}
int Subsystem::Guts::getNUDotErr(const State& s)     const {return s.getNUDotErr(getRep().getMySubsystemId());}


void Subsystem::Guts::invalidateSubsystemTopologyCache() const {
    getRep().invalidateSubsystemTopologyCache();
}

bool Subsystem::Guts::subsystemTopologyHasBeenRealized() const {
    return getRep().subsystemTopologyHasBeenRealized();
}

    // wrappers for Subsystem::Guts virtuals


Subsystem::Guts* Subsystem::Guts::clone() const {
    return getRep().clonep(*this);
}

void Subsystem::Guts::realizeSubsystemTopology(State& s) const {
    SimTK_STAGECHECK_EQ_ALWAYS(getStage(s), Stage::Empty, 
        "Subsystem::Guts::realizeSubsystemTopology()");
    getRep().realizeTopologyp(*this,s);
    getRep().subsystemTopologyRealized = true; // mark the subsystem itself (mutable)
    advanceToStage(s, Stage::Topology);  // mark the State as well
}
void Subsystem::Guts::realizeSubsystemModel(State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage::Topology, 
        "Subsystem::Guts::realizeSubsystemModel()");
    if (getStage(s) < Stage::Model) {
        getRep().realizeModelp(*this,s);
        advanceToStage(s, Stage::Model);
    }
}
void Subsystem::Guts::realizeSubsystemInstance(const State& s) const { 
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Instance).prev(), 
        "Subsystem::Guts::realizeSubsystemInstance()");
    if (getStage(s) < Stage::Instance) {
        getRep().realizeInstancep(*this,s);
        advanceToStage(s, Stage::Instance);
    }
}
void Subsystem::Guts::realizeSubsystemTime(const State& s) const { 
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Time).prev(), 
        "Subsystem::Guts::realizeTime()");
    if (getStage(s) < Stage::Time) {
        getRep().realizeTimep(*this,s);
        advanceToStage(s, Stage::Time);
    }
}
void Subsystem::Guts::realizeSubsystemPosition(const State& s) const { 
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Position).prev(), 
        "Subsystem::Guts::realizeSubsystemPosition()");
    if (getStage(s) < Stage::Position) {
        getRep().realizePositionp(*this,s);
        advanceToStage(s, Stage::Position);
    }
}
void Subsystem::Guts::realizeSubsystemVelocity(const State& s) const { 
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Velocity).prev(), 
        "Subsystem::Guts::realizeSubsystemVelocity()");
    if (getStage(s) < Stage::Velocity) {
        getRep().realizeVelocityp(*this,s);
        advanceToStage(s, Stage::Velocity);
    }
}
void Subsystem::Guts::realizeSubsystemDynamics(const State& s) const { 
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Dynamics).prev(), 
        "Subsystem::Guts::realizeSubsystemDynamics()");
    if (getStage(s) < Stage::Dynamics) {
        getRep().realizeDynamicsp(*this,s);
        advanceToStage(s, Stage::Dynamics);
    }
}
void Subsystem::Guts::realizeSubsystemAcceleration(const State& s) const { 
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Acceleration).prev(), 
        "Subsystem::Guts::realizeSubsystemAcceleration()");
    if (getStage(s) < Stage::Acceleration) {
        getRep().realizeAccelerationp(*this,s);
        advanceToStage(s, Stage::Acceleration);
    }
}
void Subsystem::Guts::realizeSubsystemReport(const State& s) const { 
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Report).prev(), 
        "Subsystem::Guts::realizeSubsystemReport()");
    if (getStage(s) < Stage::Report) {
        getRep().realizeReportp(*this,s);
        advanceToStage(s, Stage::Report);
    }
}


void Subsystem::Guts::calcQUnitWeights(const State& s, Vector& weights) const {
    getRep().calcQUnitWeightsp(*this,s,weights);
}
void Subsystem::Guts::calcUUnitWeights(const State& s, Vector& weights) const {
    getRep().calcUUnitWeightsp(*this,s,weights);
}
void Subsystem::Guts::calcZUnitWeights(const State& s, Vector& weights) const {
    getRep().calcZUnitWeightsp(*this,s,weights);
}
void Subsystem::Guts::calcQErrUnitTolerances(const State& s, Vector& tolerances) const {
    getRep().calcQErrUnitTolerancesp(*this,s,tolerances);
}
void Subsystem::Guts::calcUErrUnitTolerances(const State& s, Vector& tolerances) const {
    getRep().calcUErrUnitTolerancesp(*this,s,tolerances);
}

void Subsystem::Guts::calcDecorativeGeometryAndAppend(const State& s, Stage stage, Array<DecorativeGeometry>& geom) const {
    getRep().calcDecorativeGeometryAndAppendp(*this,s,stage,geom);
}



    // default implementations for Subsystem::Guts virtuals
/*virtual*/ int Subsystem::Guts::realizeSubsystemTopologyImpl(State& s) const {
    return 0;
}
/*virtual*/ int Subsystem::Guts::realizeSubsystemModelImpl(State& s) const {
    return 0; 
}
/*virtual*/ int Subsystem::Guts::realizeSubsystemInstanceImpl(const State& s) const {
    return 0; 
}
/*virtual*/ int Subsystem::Guts::realizeSubsystemTimeImpl(const State& s) const {
    return 0; 
}
/*virtual*/ int Subsystem::Guts::realizeSubsystemPositionImpl(const State& s) const {
    return 0; 
}
/*virtual*/ int Subsystem::Guts::realizeSubsystemVelocityImpl(const State& s) const {
    return 0; 
}
/*virtual*/ int Subsystem::Guts::realizeSubsystemDynamicsImpl(const State& s) const {
    return 0; 
}
/*virtual*/ int Subsystem::Guts::realizeSubsystemAccelerationImpl(const State& s) const {
    return 0; 
}
/*virtual*/ int Subsystem::Guts::realizeSubsystemReportImpl(const State& s) const {
    return 0; 
}

/*virtual*/ int Subsystem::Guts::calcQUnitWeightsImpl(const State& s, Vector& weights) const {
    weights.resize(getNQ(s));
    weights = 1; // default says everyone's opinion is just as valid
    return 0;
}
/*virtual*/ int Subsystem::Guts::calcUUnitWeightsImpl(const State& s, Vector& weights) const {
    weights.resize(getNU(s));
    weights = 1;
    return 0;
}
/*virtual*/ int Subsystem::Guts::calcZUnitWeightsImpl(const State& s, Vector& weights) const {
    weights.resize(getNZ(s));
    weights = 1;
    return 0;
}
/*virtual*/ int Subsystem::Guts::calcQErrUnitTolerancesImpl(const State& s, Vector& tolerances) const {
    tolerances.resize(getNQErr(s));
    tolerances = 1;
    return 0;
}
/*virtual*/ int Subsystem::Guts::calcUErrUnitTolerancesImpl(const State& s, Vector& tolerances) const {
    tolerances.resize(getNUErr(s));
    tolerances = 1;
    return 0;
}
/*virtual*/ int Subsystem::Guts::calcDecorativeGeometryAndAppendImpl
                                (const State&, Stage, Array<DecorativeGeometry>&) const
{
    return 0;
}

    ///////////////////
    // SUBSYSTEM REP //
    ///////////////////

void Subsystem::Guts::GutsRep::invalidateSubsystemTopologyCache() const {
    subsystemTopologyRealized = false;
    if (isInSystem()) 
        getSystem().getSystemGuts().invalidateSystemTopologyCache();
}

    //////////////////////////////
    // DEFAULT SYSTEM SUBSYSTEM //
    //////////////////////////////

class DefaultSystemSubsystemGuts : public Subsystem::Guts {
public:
    DefaultSystemSubsystemGuts() : Guts("DefaultSystemSubsystemGuts", "0.0.1") { }
    DefaultSystemSubsystemGuts* cloneImpl() const {
        return new DefaultSystemSubsystemGuts(*this);
    }

};

DefaultSystemSubsystem::DefaultSystemSubsystem(System& sys) {
    adoptSubsystemGuts(new DefaultSystemSubsystemGuts());
    sys.adoptSubsystem(*this);
}

} // namespace SimTK

