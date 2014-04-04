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
 * Contributors: Peter Eastman                                                *
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
 * Implementation of Subsystem, Subsystem::Guts and DefaultSystemSubsystem.
 */

#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"
#include "SimTKcommon/internal/Measure.h"
#include "SimTKcommon/internal/State.h"
#include "SimTKcommon/internal/EventHandler.h"
#include "SimTKcommon/internal/EventReporter.h"
#include "SimTKcommon/internal/System.h"
#include "SimTKcommon/internal/Subsystem.h"

#include "SimTKcommon/internal/MeasureImplementation.h"

#include "SystemGutsRep.h"

#include <cassert>

namespace SimTK {

//==============================================================================
//                                 SUBSYSTEM
//==============================================================================

bool Subsystem::isEmptyHandle() const {return guts==0;}
bool Subsystem::isOwnerHandle() const 
{   return guts==0 || &guts->getOwnerSubsystemHandle()==this; }
bool Subsystem::isSameSubsystem(const Subsystem& otherSubsystem) const 
{   return guts && (guts==otherSubsystem.guts); }

Subsystem::Subsystem(const Subsystem& src) : guts(0) {
    if (src.guts) {
        guts = src.guts->clone();
        guts->setOwnerSubsystemHandle(*this);
    }
}

Subsystem& Subsystem::operator=(const Subsystem& src) {
    if (!isSameSubsystem(src)) {
        if (isOwnerHandle())
            delete guts;
        guts=0;
        if (src.guts) {
            guts = src.guts->clone();
            guts->setOwnerSubsystemHandle(*this);
        }
    }
    return *this;
}

Subsystem::~Subsystem() {
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

void Subsystem::setSystem(System& sys, SubsystemIndex id) {
    updSubsystemGuts().setSystem(sys,id);
}

const String& Subsystem::getName()    const {return getSubsystemGuts().getName();}
const String& Subsystem::getVersion() const {return getSubsystemGuts().getVersion();}

bool Subsystem::subsystemTopologyHasBeenRealized() const {
    return getSubsystemGuts().subsystemTopologyHasBeenRealized();
}

void Subsystem::invalidateSubsystemTopologyCache() const {
    getSubsystemGuts().invalidateSubsystemTopologyCache(); // mutable
}

MeasureIndex Subsystem::adoptMeasure(AbstractMeasure& m)
{   return updSubsystemGuts().adoptMeasure(m); }
AbstractMeasure Subsystem::getMeasure(MeasureIndex mx) const
{   return getSubsystemGuts().getMeasure(mx); }


bool Subsystem::isInSystem() const {return getSubsystemGuts().isInSystem();}
bool Subsystem::isInSameSystem(const Subsystem& otherSubsystem) const {
    return getSubsystemGuts().isInSameSystem(otherSubsystem);
}

const System& Subsystem::getSystem() const {return getSubsystemGuts().getSystem();}
System&       Subsystem::updSystem()       {return updSubsystemGuts().updSystem();}

SubsystemIndex Subsystem::getMySubsystemIndex() const {
    return getSubsystemGuts().getMySubsystemIndex();
}

QIndex Subsystem::allocateQ(State& s, const Vector& qInit)const {return getSubsystemGuts().allocateQ(s,qInit);}
UIndex Subsystem::allocateU(State& s, const Vector& uInit)const {return getSubsystemGuts().allocateU(s,uInit);}
ZIndex Subsystem::allocateZ(State& s, const Vector& zInit)const {return getSubsystemGuts().allocateZ(s,zInit);}

DiscreteVariableIndex Subsystem::allocateDiscreteVariable
   (State& s, Stage g, AbstractValue* v) const 
{   return getSubsystemGuts().allocateDiscreteVariable(s,g,v); }
DiscreteVariableIndex Subsystem::allocateAutoUpdateDiscreteVariable
   (State& s, Stage invalidates, AbstractValue* v, Stage updateDependsOn) const
{   return getSubsystemGuts().allocateAutoUpdateDiscreteVariable(s,invalidates,v,updateDependsOn); }

CacheEntryIndex Subsystem::allocateCacheEntry   
   (const State& s, Stage dependsOn, Stage computedBy, AbstractValue* v) const 
{   return getSubsystemGuts().allocateCacheEntry(s,dependsOn,computedBy,v); }

QErrIndex       Subsystem::allocateQErr         (const State& s, int nqerr)    const {return getSubsystemGuts().allocateQErr(s,nqerr);}
UErrIndex       Subsystem::allocateUErr         (const State& s, int nuerr)    const {return getSubsystemGuts().allocateUErr(s,nuerr);}
UDotErrIndex    Subsystem::allocateUDotErr      (const State& s, int nudoterr) const {return getSubsystemGuts().allocateUDotErr(s,nudoterr);}
EventTriggerByStageIndex Subsystem::allocateEventTriggersByStage(const State& s, Stage g, int ntriggers) const {return getSubsystemGuts().allocateEventTriggersByStage(s,g,ntriggers);}

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
const Vector& Subsystem::getMultipliers(const State& s) const {return getSubsystemGuts().getMultipliers(s);}
const Vector& Subsystem::getEventTriggersByStage(const State& s, Stage g) const {return getSubsystemGuts().getEventTriggersByStage(s,g);}


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
Vector& Subsystem::updMultipliers(const State& s) const {return getSubsystemGuts().updMultipliers(s);}
Vector& Subsystem::updEventTriggersByStage(const State& s, Stage g) const {return getSubsystemGuts().updEventTriggersByStage(s,g);}

Stage Subsystem::getStage(const State& s) const {return getSubsystemGuts().getStage(s);}
const AbstractValue& Subsystem::getDiscreteVariable(const State& s, DiscreteVariableIndex index) const {
    return getSubsystemGuts().getDiscreteVariable(s, index);
}
AbstractValue& Subsystem::updDiscreteVariable(State& s, DiscreteVariableIndex index) const {
    return getSubsystemGuts().updDiscreteVariable(s, index);
}
const AbstractValue& Subsystem::getCacheEntry(const State& s, CacheEntryIndex index) const {
    return getSubsystemGuts().getCacheEntry(s, index);
}
AbstractValue& Subsystem::updCacheEntry(const State& s, CacheEntryIndex index) const {
    return getSubsystemGuts().updCacheEntry(s, index);
}
bool Subsystem::isCacheValueRealized(const State& s, CacheEntryIndex cx) const {
    return getSubsystemGuts().isCacheValueRealized(s, cx);
}
void Subsystem::markCacheValueRealized(const State& s, CacheEntryIndex cx) const {
    getSubsystemGuts().markCacheValueRealized(s, cx);
}
void Subsystem::markCacheValueNotRealized(const State& s, CacheEntryIndex cx) const {
    getSubsystemGuts().markCacheValueNotRealized(s, cx);
}

SystemQIndex Subsystem::getQStart      (const State& s) const {return getSubsystemGuts().getQStart(s);}
int Subsystem::getNQ          (const State& s) const {return getSubsystemGuts().getNQ(s);}
SystemUIndex Subsystem::getUStart      (const State& s) const {return getSubsystemGuts().getUStart(s);}
int Subsystem::getNU          (const State& s) const {return getSubsystemGuts().getNU(s);}
SystemZIndex Subsystem::getZStart      (const State& s) const {return getSubsystemGuts().getZStart(s);}
int Subsystem::getNZ          (const State& s) const {return getSubsystemGuts().getNZ(s);}
SystemQErrIndex Subsystem::getQErrStart   (const State& s) const {return getSubsystemGuts().getQErrStart(s);}
int Subsystem::getNQErr       (const State& s) const {return getSubsystemGuts().getNQErr(s);}
SystemUErrIndex Subsystem::getUErrStart   (const State& s) const {return getSubsystemGuts().getUErrStart(s);}
int Subsystem::getNUErr       (const State& s) const {return getSubsystemGuts().getNUErr(s);}
SystemUDotErrIndex Subsystem::getUDotErrStart(const State& s) const {return getSubsystemGuts().getUDotErrStart(s);}
int Subsystem::getNUDotErr    (const State& s) const {return getSubsystemGuts().getNUDotErr(s);}
SystemMultiplierIndex Subsystem::getMultipliersStart(const State& s) const {return getSubsystemGuts().getMultipliersStart(s);}
int Subsystem::getNMultipliers    (const State& s) const {return getSubsystemGuts().getNMultipliers(s);}
SystemEventTriggerByStageIndex Subsystem::getEventTriggerStartByStage(const State& s, Stage g) const {return getSubsystemGuts().getEventTriggerStartByStage(s,g);}
int Subsystem::getNEventTriggersByStage   (const State& s, Stage g) const {return getSubsystemGuts().getNEventTriggersByStage(s,g);}



//==============================================================================
//                           SUBSYSTEM :: GUTS
//==============================================================================

// This serves as default constructor.
Subsystem::Guts::Guts(const String& name, const String& version)
:   m_subsystemName(name), m_subsystemVersion(version),
    m_mySystem(0), m_mySubsystemIndex(InvalidSubsystemIndex), m_myHandle(0),
    m_subsystemTopologyRealized(false)
{ 
}

// Copy constructor isn't very useful. Note that it doesn't copy Measures.
Subsystem::Guts::Guts(const Subsystem::Guts& src) 
:   m_subsystemName(src.m_subsystemName), 
    m_subsystemVersion(src.m_subsystemVersion),
    m_mySystem(0), m_mySubsystemIndex(InvalidSubsystemIndex), m_myHandle(0),
    m_subsystemTopologyRealized(false)
{
}

// Destructor must unreference and possibly delete measures.
Subsystem::Guts::~Guts() {
    for (MeasureIndex mx(0); mx < m_measures.size(); ++mx)
        if (m_measures[mx]->decrRefCount()==0) delete m_measures[mx];
    m_myHandle = 0;
    invalidateSubsystemTopologyCache();
}

  

MeasureIndex Subsystem::Guts::adoptMeasure(AbstractMeasure& m) {
    SimTK_ASSERT(m.hasImpl(), "Subsystem::Guts::adoptMeasure()");

    // In Debug mode check that this measure hasn't already been adopted.
    // This is an expensive check if there are lots of measures.
    SimTK_ASSERT(std::find(m_measures.begin(), m_measures.end(), &m.getImpl())
                 == m_measures.end(), "Subsystem::Guts::adoptMeasure()");

    invalidateSubsystemTopologyCache();
    const MeasureIndex mx(m_measures.size());
    m_measures.push_back(&m.updImpl());
    m_measures.back()->incrRefCount();
    m_measures.back()->setSubsystem(updOwnerSubsystemHandle(), mx);
    return mx;
}

bool Subsystem::Guts::isInSameSystem(const Subsystem& otherSubsystem) const {
	return isInSystem() && otherSubsystem.isInSystem()
        && getSystem().isSameSystem(otherSubsystem.getSystem());
}

QIndex Subsystem::Guts::allocateQ(State& s, const Vector& qInit) const {
    return s.allocateQ(getMySubsystemIndex(), qInit);
}

UIndex Subsystem::Guts::allocateU(State& s, const Vector& uInit) const {
    return s.allocateU(getMySubsystemIndex(), uInit);
}

ZIndex Subsystem::Guts::allocateZ(State& s, const Vector& zInit) const {
    return s.allocateZ(getMySubsystemIndex(), zInit);
}

DiscreteVariableIndex Subsystem::Guts::allocateDiscreteVariable
   (State& s, Stage g, AbstractValue* v) const 
{   return s.allocateDiscreteVariable(getMySubsystemIndex(), g, v); }
DiscreteVariableIndex Subsystem::Guts::allocateAutoUpdateDiscreteVariable
   (State& s, Stage invalidates, AbstractValue* v, Stage updateDependsOn) const
{   return s.allocateAutoUpdateDiscreteVariable
               (getMySubsystemIndex(),invalidates,v,updateDependsOn); }

CacheEntryIndex Subsystem::Guts::allocateCacheEntry
   (const State& s, Stage dependsOn, Stage computedBy, AbstractValue* v) const 
{   return s.allocateCacheEntry
               (getMySubsystemIndex(), dependsOn, computedBy, v); }

QErrIndex Subsystem::Guts::allocateQErr(const State& s, int nqerr) const {
    return s.allocateQErr(getMySubsystemIndex(), nqerr);
}

UErrIndex Subsystem::Guts::allocateUErr(const State& s, int nuerr) const {
    return s.allocateUErr(getMySubsystemIndex(), nuerr);
}

UDotErrIndex Subsystem::Guts::
allocateUDotErr(const State& s, int nudoterr) const {
    return s.allocateUDotErr(getMySubsystemIndex(), nudoterr);
}

EventTriggerByStageIndex Subsystem::Guts::
allocateEventTriggersByStage(const State& s, Stage g, int ntriggers) const {
    return s.allocateEventTrigger(getMySubsystemIndex(),g,ntriggers);
}


Stage Subsystem::Guts::getStage(const State& s) const {
    return s.getSubsystemStage(getMySubsystemIndex());
}
const AbstractValue& Subsystem::Guts::
getDiscreteVariable(const State& s, DiscreteVariableIndex index) const {
    return s.getDiscreteVariable(getMySubsystemIndex(), index);
}

AbstractValue& Subsystem::Guts::
updDiscreteVariable(State& s, DiscreteVariableIndex index) const {
    return s.updDiscreteVariable(getMySubsystemIndex(), index);
}

const AbstractValue& Subsystem::Guts::
getCacheEntry(const State& s, CacheEntryIndex index) const {
    return s.getCacheEntry(getMySubsystemIndex(), index);
}

AbstractValue& Subsystem::Guts::
updCacheEntry(const State& s, CacheEntryIndex index) const {
    return s.updCacheEntry(getMySubsystemIndex(), index);
}

bool Subsystem::Guts::
isCacheValueRealized(const State& s, CacheEntryIndex cx) const {
    return s.isCacheValueRealized(getMySubsystemIndex(), cx);
}
void Subsystem::Guts::
markCacheValueRealized(const State& s, CacheEntryIndex cx) const {
    s.markCacheValueRealized(getMySubsystemIndex(), cx);
}
void Subsystem::Guts::
markCacheValueNotRealized(const State& s, CacheEntryIndex cx) const {
    s.markCacheValueNotRealized(getMySubsystemIndex(), cx);
}

const Vector& Subsystem::Guts::getQ(const State& s) const {return s.getQ(getMySubsystemIndex());}
const Vector& Subsystem::Guts::getU(const State& s) const {return s.getU(getMySubsystemIndex());}
const Vector& Subsystem::Guts::getZ(const State& s) const {return s.getZ(getMySubsystemIndex());}
const Vector& Subsystem::Guts::getUWeights(const State& s) const {return s.getUWeights(getMySubsystemIndex());}
const Vector& Subsystem::Guts::getZWeights(const State& s) const {return s.getZWeights(getMySubsystemIndex());}

Vector& Subsystem::Guts::updQ(State& s) const {return s.updQ(getMySubsystemIndex());}
Vector& Subsystem::Guts::updU(State& s) const {return s.updU(getMySubsystemIndex());}
Vector& Subsystem::Guts::updZ(State& s) const {return s.updZ(getMySubsystemIndex());}

const Vector& Subsystem::Guts::getQDot   (const State& s) const {return s.getQDot(getMySubsystemIndex());}
const Vector& Subsystem::Guts::getUDot   (const State& s) const {return s.getUDot(getMySubsystemIndex());}
const Vector& Subsystem::Guts::getZDot   (const State& s) const {return s.getZDot(getMySubsystemIndex());}
const Vector& Subsystem::Guts::getQDotDot(const State& s) const {return s.getQDotDot(getMySubsystemIndex());}

Vector& Subsystem::Guts::updQDot   (const State& s) const {return s.updQDot(getMySubsystemIndex());}
Vector& Subsystem::Guts::updUDot   (const State& s) const {return s.updUDot(getMySubsystemIndex());}
Vector& Subsystem::Guts::updZDot   (const State& s) const {return s.updZDot(getMySubsystemIndex());}
Vector& Subsystem::Guts::updQDotDot(const State& s) const {return s.updQDotDot(getMySubsystemIndex());}

const Vector& Subsystem::Guts::getQErr(const State& s) const {return s.getQErr(getMySubsystemIndex());}
const Vector& Subsystem::Guts::getUErr(const State& s) const {return s.getUErr(getMySubsystemIndex());}
const Vector& Subsystem::Guts::getQErrWeights(const State& s) const {return s.getQErrWeights(getMySubsystemIndex());}
const Vector& Subsystem::Guts::getUErrWeights(const State& s) const {return s.getUErrWeights(getMySubsystemIndex());}

const Vector& Subsystem::Guts::getUDotErr(const State& s) const {return s.getUDotErr(getMySubsystemIndex());}
const Vector& Subsystem::Guts::getMultipliers(const State& s) const {return s.getMultipliers(getMySubsystemIndex());}
const Vector& Subsystem::Guts::getEventTriggersByStage(const State& s, Stage g) const
{   return s.getEventTriggersByStage(getMySubsystemIndex(),g); }

Vector& Subsystem::Guts::updQErr(const State& s) const {return s.updQErr(getMySubsystemIndex());}
Vector& Subsystem::Guts::updUErr(const State& s) const {return s.updUErr(getMySubsystemIndex());}
Vector& Subsystem::Guts::updUDotErr(const State& s) const {return s.updUDotErr(getMySubsystemIndex());}
Vector& Subsystem::Guts::updMultipliers(const State& s) const {return s.updMultipliers(getMySubsystemIndex());}
Vector& Subsystem::Guts::updEventTriggersByStage(const State& s, Stage g) const
{   return s.updEventTriggersByStage(getMySubsystemIndex(),g); }

SystemQIndex Subsystem::Guts::getQStart(const State& s) const {return s.getQStart(getMySubsystemIndex());}
int Subsystem::Guts::getNQ(const State& s)     const {return s.getNQ(getMySubsystemIndex());}

SystemUIndex Subsystem::Guts::getUStart(const State& s) const {return s.getUStart(getMySubsystemIndex());}
int Subsystem::Guts::getNU(const State& s)     const {return s.getNU(getMySubsystemIndex());}

SystemZIndex Subsystem::Guts::getZStart(const State& s) const {return s.getZStart(getMySubsystemIndex());}
int Subsystem::Guts::getNZ(const State& s)     const {return s.getNZ(getMySubsystemIndex());}

SystemQErrIndex Subsystem::Guts::getQErrStart(const State& s) const {return s.getQErrStart(getMySubsystemIndex());}
int Subsystem::Guts::getNQErr(const State& s)     const {return s.getNQErr(getMySubsystemIndex());}

SystemUErrIndex Subsystem::Guts::getUErrStart(const State& s) const {return s.getUErrStart(getMySubsystemIndex());}
int Subsystem::Guts::getNUErr(const State& s)     const {return s.getNUErr(getMySubsystemIndex());}

SystemUDotErrIndex Subsystem::Guts::getUDotErrStart(const State& s) const {return s.getUDotErrStart(getMySubsystemIndex());}
int Subsystem::Guts::getNUDotErr(const State& s)     const {return s.getNUDotErr(getMySubsystemIndex());}

SystemMultiplierIndex Subsystem::Guts::getMultipliersStart(const State& s) const {return s.getMultipliersStart(getMySubsystemIndex());}
int Subsystem::Guts::getNMultipliers(const State& s)     const {return s.getNMultipliers(getMySubsystemIndex());}

SystemEventTriggerByStageIndex Subsystem::Guts::getEventTriggerStartByStage(const State& s, Stage g) const {return s.getEventTriggerStartByStage(getMySubsystemIndex(),g);}
int Subsystem::Guts::getNEventTriggersByStage   (const State& s, Stage g) const {return s.getNEventTriggersByStage(getMySubsystemIndex(),g);}

// Invalidating a Subsystem's topology cache forces invalidation of the
// whole System's topology cache, which will in turn invalidate all the other
// Subsystem's topology caches.
void Subsystem::Guts::invalidateSubsystemTopologyCache() const {
    if (m_subsystemTopologyRealized) {
        m_subsystemTopologyRealized = false;
        if (isInSystem()) 
            getSystem().getSystemGuts().invalidateSystemTopologyCache();
    }
}


//------------------------------------------------------------------------------
//                         CREATE SCHEDULED EVENT
//------------------------------------------------------------------------------
/*
 * A Subsystem should invoke this method during Instance stage for each 
 * scheduled event it defines.
 * It allocates a global event ID for the event, and registers that ID as 
 * belonging to this Subsystem.
 * 
 * @param state     the State which is being realized
 * @param eventId   on exit, the newly allocated event ID is stored here
 */
void Subsystem::Guts::
createScheduledEvent(const State& state, EventId& eventId) const {
    eventId = getSystem().getDefaultSubsystem()
                         .createEventId(getMySubsystemIndex(), state);
}

//------------------------------------------------------------------------------
//                         CREATE TRIGGERED EVENT
//------------------------------------------------------------------------------
/*
 * A Subsystem should invoke this method during Instance stage for each 
 * triggered event it defines. It allocates a global event ID for the event, 
 * registers that ID as belonging to this Subsystem, and allocates space in the
 * State for the event trigger function.
 * 
 * @param state     the State which is being realized
 * @param eventId   on exit, the newly allocated event ID is stored here
 * @param triggerFunctionIndex  
 *      on exit, the index corresponding to the event's trigger function
 *      is stored here (this is a local, per-Subsystem, per-Stage index)
 * @param stage     the Stage at which the event will be evaluated
 */
void Subsystem::Guts::
createTriggeredEvent(const State& state, EventId& eventId, 
                     EventTriggerByStageIndex& triggerFunctionIndex, 
                     Stage stage) const 
{
    eventId = getSystem().getDefaultSubsystem()
                         .createEventId(getMySubsystemIndex(), state);
    triggerFunctionIndex = 
        state.allocateEventTrigger(getMySubsystemIndex(), stage, 1);
}


    // wrappers for Subsystem::Guts virtuals

//------------------------------------------------------------------------------
//                                  CLONE
//------------------------------------------------------------------------------
Subsystem::Guts* Subsystem::Guts::clone() const {
    return cloneImpl();
}

//------------------------------------------------------------------------------
//                     REALIZE SUBSYSTEM TOPOLOGY
//------------------------------------------------------------------------------
void Subsystem::Guts::realizeSubsystemTopology(State& s) const {
    SimTK_STAGECHECK_EQ_ALWAYS(getStage(s), Stage::Empty, 
        "Subsystem::Guts::realizeSubsystemTopology()");
    realizeSubsystemTopologyImpl(s);

    // Realize this Subsystem's Measures.
    for (MeasureIndex mx(0); mx < m_measures.size(); ++mx)
        m_measures[mx]->realizeTopology(s);

    m_subsystemTopologyRealized = true; // mark subsys itself (mutable)
    advanceToStage(s, Stage::Topology);  // mark the State as well
}

//------------------------------------------------------------------------------
//                        REALIZE SUBSYSTEM MODEL
//------------------------------------------------------------------------------
void Subsystem::Guts::realizeSubsystemModel(State& s) const {
    SimTK_STAGECHECK_TOPOLOGY_REALIZED_ALWAYS(subsystemTopologyHasBeenRealized(),
        "Subsystem", getName(), "Subsystem::Guts::realizeSubsystemModel()");

    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage::Topology, 
        "Subsystem::Guts::realizeSubsystemModel()");
    if (getStage(s) < Stage::Model) {
        realizeSubsystemModelImpl(s);

        // Realize this Subsystem's Measures.
        for (MeasureIndex mx(0); mx < m_measures.size(); ++mx)
            m_measures[mx]->realizeModel(s);

        advanceToStage(s, Stage::Model);
    }
}

//------------------------------------------------------------------------------
//                     REALIZE SUBSYSTEM INSTANCE
//------------------------------------------------------------------------------
void Subsystem::Guts::realizeSubsystemInstance(const State& s) const { 
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Instance).prev(), 
        "Subsystem::Guts::realizeSubsystemInstance()");
    if (getStage(s) < Stage::Instance) {
        realizeSubsystemInstanceImpl(s);

        // Realize this Subsystem's Measures.
        for (MeasureIndex mx(0); mx < m_measures.size(); ++mx)
            m_measures[mx]->realizeInstance(s);

        advanceToStage(s, Stage::Instance);
    }
}

//------------------------------------------------------------------------------
//                         REALIZE SUBSYSTEM TIME
//------------------------------------------------------------------------------
void Subsystem::Guts::realizeSubsystemTime(const State& s) const { 
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Time).prev(), 
        "Subsystem::Guts::realizeTime()");
    if (getStage(s) < Stage::Time) {
        realizeSubsystemTimeImpl(s);

        // Realize this Subsystem's Measures.
        for (MeasureIndex mx(0); mx < m_measures.size(); ++mx)
            m_measures[mx]->realizeTime(s);

        advanceToStage(s, Stage::Time);
    }
}

//------------------------------------------------------------------------------
//                     REALIZE SUBSYSTEM POSITION
//------------------------------------------------------------------------------
void Subsystem::Guts::realizeSubsystemPosition(const State& s) const { 
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Position).prev(), 
        "Subsystem::Guts::realizeSubsystemPosition()");
    if (getStage(s) < Stage::Position) {
        realizeSubsystemPositionImpl(s);

        // Realize this Subsystem's Measures.
        for (MeasureIndex mx(0); mx < m_measures.size(); ++mx)
            m_measures[mx]->realizePosition(s);

        advanceToStage(s, Stage::Position);
    }
}

//------------------------------------------------------------------------------
//                       REALIZE SUBSYSTEM VELOCITY
//------------------------------------------------------------------------------
void Subsystem::Guts::realizeSubsystemVelocity(const State& s) const { 
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Velocity).prev(), 
        "Subsystem::Guts::realizeSubsystemVelocity()");
    if (getStage(s) < Stage::Velocity) {
        realizeSubsystemVelocityImpl(s);

        // Realize this Subsystem's Measures.
        for (MeasureIndex mx(0); mx < m_measures.size(); ++mx)
            m_measures[mx]->realizeVelocity(s);

        advanceToStage(s, Stage::Velocity);
    }
}

//------------------------------------------------------------------------------
//                       REALIZE SUBSYSTEM DYNAMICS
//------------------------------------------------------------------------------
void Subsystem::Guts::realizeSubsystemDynamics(const State& s) const { 
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Dynamics).prev(), 
        "Subsystem::Guts::realizeSubsystemDynamics()");
    if (getStage(s) < Stage::Dynamics) {
        realizeSubsystemDynamicsImpl(s);

        // Realize this Subsystem's Measures.
        for (MeasureIndex mx(0); mx < m_measures.size(); ++mx)
            m_measures[mx]->realizeDynamics(s);

        advanceToStage(s, Stage::Dynamics);
    }
}

//------------------------------------------------------------------------------
//                     REALIZE SUBSYSTEM ACCELERATION
//------------------------------------------------------------------------------
void Subsystem::Guts::realizeSubsystemAcceleration(const State& s) const { 
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Acceleration).prev(), 
        "Subsystem::Guts::realizeSubsystemAcceleration()");
    if (getStage(s) < Stage::Acceleration) {
        realizeSubsystemAccelerationImpl(s);

        // Realize this Subsystem's Measures.
        for (MeasureIndex mx(0); mx < m_measures.size(); ++mx)
            m_measures[mx]->realizeAcceleration(s);

        advanceToStage(s, Stage::Acceleration);
    }
}

//------------------------------------------------------------------------------
//                         REALIZE SUBSYSTEM REPORT
//------------------------------------------------------------------------------
void Subsystem::Guts::realizeSubsystemReport(const State& s) const { 
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Report).prev(), 
        "Subsystem::Guts::realizeSubsystemReport()");
    if (getStage(s) < Stage::Report) {
        realizeSubsystemReportImpl(s);

        // Realize this Subsystem's Measures.
        for (MeasureIndex mx(0); mx < m_measures.size(); ++mx)
            m_measures[mx]->realizeReport(s);

        advanceToStage(s, Stage::Report);
    }
}

//------------------------------------------------------------------------------
//                  CALC DECORATIVE GEOMETRY AND APPEND
//------------------------------------------------------------------------------
void Subsystem::Guts::calcDecorativeGeometryAndAppend
   (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const 
{
    calcDecorativeGeometryAndAppendImpl(s,stage,geom);
}

//------------------------------------------------------------------------------
//                              HANDLE EVENTS
//------------------------------------------------------------------------------
void Subsystem::Guts::handleEvents(State& state, Event::Cause cause, 
    const Array_<EventId>& eventIds,
    const HandleEventsOptions& options, HandleEventsResults& results) const
{
    // Invoke Measure handlers where appropriate (TODO: just initialization
    // so far). Initialize measures first in case the Subsystem initialization
    // handler references measures.
    if (cause == Event::Cause::Initialization) {
        for (MeasureIndex mx(0); mx < m_measures.size(); ++mx)
            m_measures[mx]->initialize(state);
    }

    // assume success
    results.clear(); results.setExitStatus(HandleEventsResults::Succeeded); 
    handleEventsImpl(state, cause, eventIds, options, results);
}

//------------------------------------------------------------------------------
//                               REPORT EVENTS
//------------------------------------------------------------------------------
void Subsystem::Guts::reportEvents(const State& state, Event::Cause cause, 
                                   const Array_<EventId>& eventIds) const
{
    reportEventsImpl(state, cause, eventIds);
}

//------------------------------------------------------------------------------
//                         CALC EVENT TRIGGER INFO
//------------------------------------------------------------------------------
void Subsystem::Guts::
calcEventTriggerInfo(const State& state, Array_<EventTriggerInfo>& info) const {
    calcEventTriggerInfoImpl(state, info);
}

//------------------------------------------------------------------------------
//                     CALC TIME OF NEXT SCHEDULED EVENT
//------------------------------------------------------------------------------
void Subsystem::Guts::
calcTimeOfNextScheduledEvent(const State& state, Real& tNextEvent, 
                             Array_<EventId>& eventIds, 
                             bool includeCurrentTime) const 
{
    tNextEvent = Infinity;
    eventIds.clear();
    calcTimeOfNextScheduledEventImpl(state, tNextEvent, eventIds, 
                                     includeCurrentTime);
}

//------------------------------------------------------------------------------
//                    CALC TIME OF NEXT SCHEDULED REPORT
//------------------------------------------------------------------------------
void Subsystem::Guts::
calcTimeOfNextScheduledReport(const State& state, Real& tNextReport, 
                              Array_<EventId>& eventIds, 
                              bool includeCurrentTime) const 
{
    tNextReport = Infinity;
    eventIds.clear();
    calcTimeOfNextScheduledReportImpl(state, tNextReport, eventIds, 
                                      includeCurrentTime);
}



} // namespace SimTK

