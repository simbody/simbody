#ifndef SimTK_SimTKCOMMON_SUBSYSTEM_GUTS_H_
#define SimTK_SimTKCOMMON_SUBSYSTEM_GUTS_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-15 Stanford University and the Authors.        *
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

#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"
#include "SimTKcommon/internal/State.h"
#include "SimTKcommon/internal/Event_Defs.h"

#include <cassert>

namespace SimTK {

class System;
class DecorativeGeometry;

//==============================================================================
//                           SUBSYSTEM :: GUTS
//==============================================================================
/** The abstract parent of all Subsystem implementation classes. You must 
extend this class if you implement a new concrete Subsystem class. The 
Simbody user's API for Subsystems is defined exclusively in the Subsystem 
handle class declaration; the API declared here is for use by Subsystem 
implementors. **/
class SimTK_SimTKCOMMON_EXPORT Subsystem::Guts {
public:

/** Destructor is virtual to permit cleanup of derived classes. **/
virtual ~Guts();

/** Copy constructor results in a default-constructed object though with
the name and version string copied. Copy assignment is suppressed. **/
Guts(const Guts&);

/** Report back the name supplied on construction; this is not interpreted
in any way by Simbody. **/
const String& getName()    const {return m_subsystemName;}
/** Report back the version string supplied on construction; this i not
interpreted in any way by Simbody. **/
const String& getVersion() const {return m_subsystemVersion;}

/** @name                   State access methods
These convenience methods are inline pass-throughs to the State methods of the 
same name but insert this %Subsystem's SubsystemIndex as the first argument. 
That is the value returned by the getMySubsystemIndex() method. An exception 
will be thrown if this %Subsystem is not contained in a System.  

See the SimTK::State documentation for the meaning of these methods; the 
behavior is identical here. An identical set of methods is present in 
the %Subsystem handle class; in both cases they were added because it got
annoying to have to dig up the subsystem index for every call to a State
method. **/
/**@{**/

QIndex allocateQ(State& s, const Vector& qInit) const 
{   return s.allocateQ(getMySubsystemIndex(), qInit); }
UIndex allocateU(State& s, const Vector& uInit) const 
{   return s.allocateU(getMySubsystemIndex(), uInit); }
ZIndex allocateZ(State& s, const Vector& zInit) const 
{   return s.allocateZ(getMySubsystemIndex(), zInit); }

DiscreteVariableIndex 
allocateDiscreteVariable(State& s, Stage g, AbstractValue* v) const 
{   return s.allocateDiscreteVariable(getMySubsystemIndex(), g, v); }
DiscreteVariableIndex allocateAutoUpdateDiscreteVariable
   (State& s, Stage invalidates, AbstractValue* v, Stage updateDependsOn) const
{   return s.allocateAutoUpdateDiscreteVariable
               (getMySubsystemIndex(),invalidates,v,updateDependsOn); }

CacheEntryIndex allocateCacheEntry
   (const State& s, Stage dependsOn, Stage computedBy, AbstractValue* v,
    unsigned numExtras=0) const 
{   return s.allocateCacheEntry
               (getMySubsystemIndex(), dependsOn, computedBy, v, numExtras); }
CacheEntryIndex allocateCacheEntry
   (const State& state, Stage g, AbstractValue* v, unsigned numExtras=0) const 
{   return allocateCacheEntry(state, g, g, v, numExtras); }
CacheEntryIndex allocateLazyCacheEntry   
   (const State& state, Stage earliest, AbstractValue* v, 
    unsigned numExtras=0) const 
{   return allocateCacheEntry(state, earliest, Stage::Infinity, v, numExtras); }

QErrIndex allocateQErr(const State& s, int nqerr) const 
{   return s.allocateQErr(getMySubsystemIndex(), nqerr); }
UErrIndex allocateUErr(const State& s, int nuerr) const 
{   return s.allocateUErr(getMySubsystemIndex(), nuerr); }
UDotErrIndex allocateUDotErr(const State& s, int nudoterr) const 
{   return s.allocateUDotErr(getMySubsystemIndex(), nudoterr); }

const Vector& getQ(const State& s) const 
{   return s.getQ(getMySubsystemIndex()); }
const Vector& getU(const State& s) const 
{   return s.getU(getMySubsystemIndex()); }
const Vector& getZ(const State& s) const 
{   return s.getZ(getMySubsystemIndex()); }
const Vector& getUWeights(const State& s) const 
{   return s.getUWeights(getMySubsystemIndex()); }
const Vector& getZWeights(const State& s) const 
{   return s.getZWeights(getMySubsystemIndex()); }

Vector& updQ(State& s) const {return s.updQ(getMySubsystemIndex());}
Vector& updU(State& s) const {return s.updU(getMySubsystemIndex());}
Vector& updZ(State& s) const {return s.updZ(getMySubsystemIndex());}

const Vector& getQDot   (const State& s) const 
{   return s.getQDot(getMySubsystemIndex()); }
const Vector& getUDot   (const State& s) const 
{   return s.getUDot(getMySubsystemIndex()); }
const Vector& getZDot   (const State& s) const 
{   return s.getZDot(getMySubsystemIndex()); }
const Vector& getQDotDot(const State& s) const 
{   return s.getQDotDot(getMySubsystemIndex()); }

Vector& updQDot   (const State& s) const 
{   return s.updQDot(getMySubsystemIndex()); }
Vector& updUDot   (const State& s) const 
{   return s.updUDot(getMySubsystemIndex()); }
Vector& updZDot   (const State& s) const 
{   return s.updZDot(getMySubsystemIndex()); }
Vector& updQDotDot(const State& s) const 
{   return s.updQDotDot(getMySubsystemIndex()); }

const Vector& getQErr(const State& s) const 
{   return s.getQErr(getMySubsystemIndex()); }
const Vector& getUErr(const State& s) const 
{   return s.getUErr(getMySubsystemIndex()); }
const Vector& getQErrWeights(const State& s) const 
{   return s.getQErrWeights(getMySubsystemIndex()); }
const Vector& getUErrWeights(const State& s) const 
{   return s.getUErrWeights(getMySubsystemIndex()); }

const Vector& getUDotErr(const State& s) const 
{   return s.getUDotErr(getMySubsystemIndex()); }
const Vector& getMultipliers(const State& s) const 
{   return s.getMultipliers(getMySubsystemIndex()); }

Vector& updQErr(const State& s) const 
{   return s.updQErr(getMySubsystemIndex()); }
Vector& updUErr(const State& s) const 
{   return s.updUErr(getMySubsystemIndex()); }
Vector& updUDotErr(const State& s) const 
{   return s.updUDotErr(getMySubsystemIndex()); }
Vector& updMultipliers(const State& s) const 
{   return s.updMultipliers(getMySubsystemIndex()); }

SystemQIndex getQStart(const State& s) const 
{   return s.getQStart(getMySubsystemIndex()); }
int getNQ(const State& s)     const 
{   return s.getNQ(getMySubsystemIndex()); }

SystemUIndex getUStart(const State& s) const 
{   return s.getUStart(getMySubsystemIndex()); }
int getNU(const State& s)     const 
{   return s.getNU(getMySubsystemIndex()); }

SystemZIndex getZStart(const State& s) const 
{   return s.getZStart(getMySubsystemIndex()); }
int getNZ(const State& s)     const 
{   return s.getNZ(getMySubsystemIndex()); }

SystemQErrIndex getQErrStart(const State& s) const 
{   return s.getQErrStart(getMySubsystemIndex()); }
int getNQErr(const State& s) const 
{   return s.getNQErr(getMySubsystemIndex()); }

SystemUErrIndex getUErrStart(const State& s) const 
{   return s.getUErrStart(getMySubsystemIndex()); }
int getNUErr(const State& s)     const 
{   return s.getNUErr(getMySubsystemIndex()); }

SystemUDotErrIndex getUDotErrStart(const State& s) const 
{   return s.getUDotErrStart(getMySubsystemIndex()); }
int getNUDotErr(const State& s)     const 
{   return s.getNUDotErr(getMySubsystemIndex()); }

SystemMultiplierIndex getMultipliersStart(const State& s) const 
{   return s.getMultipliersStart(getMySubsystemIndex()); }
int getNMultipliers(const State& s)     const 
{   return s.getNMultipliers(getMySubsystemIndex()); }


// For convenience.
void setQ(State& s, const Vector& q) const {
    SimTK_ASSERT(q.size() == getNQ(s), "Subsystem::Guts::setQ()");
    updQ(s) = q;
}
void setU(State& s, const Vector& u) const {
    SimTK_ASSERT(u.size() == getNU(s), "Subsystem::Guts::setU()");
    updU(s) = u;
}
void setZ(State& s, const Vector& z) const {
    SimTK_ASSERT(z.size() == getNZ(s), "Subsystem::Guts::setZ()");
    updZ(s) = z;
}

Stage getStage(const State& s) const 
{   return s.getSubsystemStage(getMySubsystemIndex()); }
void advanceToStage(const State& s, Stage g) const 
{   s.advanceSubsystemToStage(getMySubsystemIndex(), g); }

const AbstractValue& 
getDiscreteVariable(const State& s, DiscreteVariableIndex index) const 
{   return s.getDiscreteVariable(getMySubsystemIndex(), index); }
AbstractValue& updDiscreteVariable(State& s, DiscreteVariableIndex index) const 
{   return s.updDiscreteVariable(getMySubsystemIndex(), index); }
const AbstractValue& getCacheEntry(const State& s, CacheEntryIndex index) const 
{   return s.getCacheEntry(getMySubsystemIndex(), index); }
AbstractValue& updCacheEntry(const State& s, CacheEntryIndex index) const 
{   return s.updCacheEntry(getMySubsystemIndex(), index); }
Real getDiscreteVarLastUpdateTime(const State& s, DiscreteVariableIndex dx) const
{   return s.getDiscreteVarLastUpdateTime(getMySubsystemIndex(),dx); }
CacheEntryIndex 
getDiscreteVarUpdateIndex(const State& s, DiscreteVariableIndex dx) const
{   return s.getDiscreteVarUpdateIndex(getMySubsystemIndex(),dx); }
const AbstractValue& 
getDiscreteVarUpdateValue(const State& s, DiscreteVariableIndex dx) const
{   return s.getDiscreteVarUpdateValue(getMySubsystemIndex(),dx); }
AbstractValue& 
updDiscreteVarUpdateValue(const State& s, DiscreteVariableIndex dx) const
{   return s.updDiscreteVarUpdateValue(getMySubsystemIndex(),dx); }
bool isDiscreteVarUpdateValueRealized
   (const State& s, DiscreteVariableIndex dx) const
{   return s.isDiscreteVarUpdateValueRealized(getMySubsystemIndex(),dx); }
void markDiscreteVarUpdateValueRealized
   (const State& s, DiscreteVariableIndex dx) const
{   return s.markDiscreteVarUpdateValueRealized(getMySubsystemIndex(),dx); }

bool isCacheValueRealized(const State& s, CacheEntryIndex cx) const 
{   return s.isCacheValueRealized(getMySubsystemIndex(), cx); }
bool isCacheValueRealizedWithExtras
   (const State& s, CacheEntryIndex cx,
    unsigned numExtras, const StageVersion* extraVersions,
    StageVersion& cacheValueVersion) const 
{   return s.isCacheValueRealizedWithExtras(getMySubsystemIndex(), cx,
                                            numExtras, extraVersions,
                                            cacheValueVersion); }

void markCacheValueRealized(const State& s, CacheEntryIndex cx) const 
{   s.markCacheValueRealized(getMySubsystemIndex(), cx); }
void markCacheValueRealizedWithExtras
   (const State& s, CacheEntryIndex cx,
    unsigned numExtras, const StageVersion* extraVersions) const 
{   s.markCacheValueRealizedWithExtras(getMySubsystemIndex(), cx,
                                       numExtras, extraVersions); }

void markCacheValueNotRealized(const State& s, CacheEntryIndex cx) const 
{   s.markCacheValueNotRealized(getMySubsystemIndex(), cx); }
/**@}**/

/** Add a new Measure to this Subsystem. The returned MeasureIndex is local
to this Subsystem and can be used to access the Measure later. **/
MeasureIndex adoptMeasure(AbstractMeasure& m);

/** Return the Measure whose index within this Subsystem is given. The 
index should be as it was returned by adoptMeasure(). This is templatized 
for use when you know the value type of the Measure; if you don't then you
should call getMeasure() which will return the Measure as an AbstractMeasure
instead. 
@see getMeasure() **/
template <class T> Measure_<T> getMeasure_(MeasureIndex mx) const
{   return Measure_<T>::getAs(getMeasure(mx)); }

/** Return the Measure whose index within this Subsystem is given, as an
AbstractMeasure (that is, its value type is not specified). The 
index should be as it was returned by adoptMeasure(). 
@see getMeasure_ **/
AbstractMeasure getMeasure(MeasureIndex mx) const {
    SimTK_ASSERT(0 <= mx && mx < m_measures.size(), 
                 "Subsystem::Guts::getMeasure()");
    //TODO: the returned object is writable; that doesn't seem like a good idea.
    return AbstractMeasure(m_measures[mx]);
}

/** Call the initialize() method of each Measure Implementation object. This
is normally called at the start of a simulation, as an action taken by the
Initialization event. **/
void initializeMeasures(State& state) const;

bool isInSystem() const {
    return m_mySystem != nullptr;
    std::cout << m_mySystem;
}

bool isInSameSystem(const Subsystem& otherSubsystem) const;

const System& getSystem() const {
    SimTK_ASSERT(isInSystem(), "Subsystem::getSystem()");
    return *m_mySystem;
}
System& updSystem() {
    SimTK_ASSERT(isInSystem(), "Subsystem::updSystem()");
    return *m_mySystem;
}

/** This Subsystem has just been added to the indicated System. Record the
System and this Subsystem's index within it, and acquire any needed 
System-global resources. **/
void setSystem(System& sys, SubsystemIndex index);

SubsystemIndex getMySubsystemIndex() const {
    SimTK_ASSERT(isInSystem(), "Subsystem::getMySubsystemIndex()");
    return m_mySubsystemIndex;
}


/** Returns \c true if this subsystem's realizeTopology() method has been
called since the last topological change or call to 
invalidateSubsystemTopologyCache(). **/
bool subsystemTopologyHasBeenRealized() const
{   return m_subsystemTopologyRealized; }

/** Always call this method when a topological change is made to this 
%Subsystem to indicate that any Stage::Topology cache values may need
recomputation. If the %Subsystem belongs to a System, the %System's overall 
topology will also be invalidated, since its Stage::Topology cannot be valid
if any of its %Subsystem topology stages are not valid. However, the stages
of other %Subsystems in the same %System are not affected. A subsequent call
to realizeTopology() is required before any computed values may be 
obtained. **/
void invalidateSubsystemTopologyCache() const;

// These are wrappers for the virtual methods defined below. They
// are used to ensure good behavior. Most of them deal automatically with
// the Subsystem's Measures, as well as invoking the corresponding virtual
// for the Subsystem's own processing.

Subsystem::Guts* clone() const;

// Realize this subsystem's part of the State from Stage-1 to Stage
// for the indicated stage. After doing some checking, these routines
// call the concrete subsystem's corresponding virtual method, and
// on return they make sure the stage has been properly updated.
// Note that these will do nothing if the Subsystem stage is already
// at or greater than the indicated stage.
void realizeSubsystemTopology    (State&) const;
void realizeSubsystemModel       (State&) const;
void realizeSubsystemInstance    (const State&) const;
void realizeSubsystemTime        (const State&) const;
void realizeSubsystemPosition    (const State&) const;
void realizeSubsystemVelocity    (const State&) const;
void realizeSubsystemDynamics    (const State&) const;
void realizeSubsystemAcceleration(const State&) const;
void realizeSubsystemReport      (const State&) const;

// Generate decorative geometry computable at a specific stage. This will
// throw an exception if this subsystem's state hasn't already been realized
// to that stage. Note that the list is not inclusive -- you have to
// request geometry from each stage to get all of it.
// The generated geometry will be *appended* to the supplied output vector.
void calcDecorativeGeometryAndAppend
    (const State&, Stage, Array_<DecorativeGeometry>&) const;

// At the time a Subsystem is added to a System (see setSystem()) it is given
// a chance to request System-level resources for its own use. This usually
// means Events and EventActions. Concrete Subsystems can implement
// acquireSystemResourcesImpl() to get more System goodies. This will throw
// an exception if the Subsystem is not already part of a System.
void acquireSystemResources();

/** Return a const reference to the Subsystem handle object that is the unique 
owner of this Subsystem::Guts object. **/
const Subsystem& getOwnerSubsystemHandle() const {
    SimTK_ASSERT(m_myHandle, "Subsystem::getOwnerSubsystemHandle()");
    return *m_myHandle;
}
/** Return a writable reference to the Subsystem handle object that is the
unique owner of this Subsystem::Guts object. **/
Subsystem& updOwnerSubsystemHandle() {
    SimTK_ASSERT(m_myHandle, "Subsystem::getOwnerSubsystemHandle()");
    return *m_myHandle;
}

/** Provide a reference to the Subsystem handle object that is the unique 
owner of this Subsystem::Guts object. The owner Subsystem is responsible for
deleting this object at destruction. Subsystem::Guts objects are not 
reference counted. **/
void setOwnerSubsystemHandle(Subsystem& subsys) {m_myHandle=&subsys;}

/** Check whether this Subsystem::Guts object is currently owned by some
Subsystem handle object. **/
bool hasOwnerSubsystemHandle() const {return m_myHandle != nullptr;}


//------------------------------------------------------------------------------
                                  protected:

/** This constructor is for use in the constructors of derived Subsystems. 
This serves as a default constructor since both arguments have defaults. 
The name and version strings are not interpreted by Simbody in any way; they
are simply stored and returned as given. **/
explicit Guts(const String& name="<NONAME>", const String& version="0.0.0");

Guts& operator=(const Guts&) = delete;

// These virtual methods should be overridden in concrete Subsystems as
// necessary. They should never be called directly; instead call the
// wrapper routines above, which have the same name but without the "Impl"
// (implementation) at the end.
    
// The "realize..." wrappers will call the "realize...Impl" methods below
// only when the current stage for the Subsystem is the one just prior
// to the stage being realized. For example, realizeSubsystemVelocityImpl()
// is called by realizeSubsystemVelocity() only when the passed-in State
// shows this subsystem's stage to be exactly Stage::Position.
//
// The default implementations provided here do nothing. That means the
// wrappers will simply check that the current stage is correct and
// advance it if necessary.

// The destructor is already virtual; see above.

virtual Subsystem::Guts* cloneImpl() const = 0;

virtual void acquireSystemResourcesImpl() {}

virtual int realizeSubsystemTopologyImpl(State& s)           const {return 0;}
virtual int realizeSubsystemModelImpl   (State& s)           const {return 0;}
virtual int realizeSubsystemInstanceImpl(const State& s)     const {return 0;}
virtual int realizeSubsystemTimeImpl    (const State& s)     const {return 0;}
virtual int realizeSubsystemPositionImpl(const State& s)     const {return 0;}
virtual int realizeSubsystemVelocityImpl(const State& s)     const {return 0;}
virtual int realizeSubsystemDynamicsImpl(const State& s)     const {return 0;}
virtual int realizeSubsystemAccelerationImpl(const State& s) const {return 0;}
virtual int realizeSubsystemReportImpl  (const State& s)     const {return 0;}

virtual int calcDecorativeGeometryAndAppendImpl
    (const State&, Stage, Array_<DecorativeGeometry>&) const {return 0;}


//------------------------------------------------------------------------------
                                    private:
friend class Subsystem;

    // TOPOLOGY STAGE STATE //
ReferencePtr<Subsystem>     m_myHandle; // the owner handle of this Guts object

String                      m_subsystemName, m_subsystemVersion;


// This is the list of Measures belonging to this Subsystem. 
// TODO: Oddly, these have shallow copy semantics so can't be ClonePtrs.
Array_<AbstractMeasure::Implementation*> 
                            m_measures;

// These are set when the Subsystem is adopted by a System.
ReferencePtr<System>        m_mySystem; // System to which this Subsys belongs
SubsystemIndex              m_mySubsystemIndex; // Subsystem # within System

    // TOPOLOGY STAGE CACHE //
bool                        m_subsystemTopologyRealized;
};


//==============================================================================
//                           SUBSYSTEM INLINES
//==============================================================================
// These had to wait for Subsystem::Guts to be defined.

inline SubsystemIndex Subsystem::getMySubsystemIndex() const
{   return getSubsystemGuts().getMySubsystemIndex(); }

inline const String& Subsystem::getName()    const 
{   return getSubsystemGuts().getName(); }
inline const String& Subsystem::getVersion() const 
{   return getSubsystemGuts().getVersion(); }

inline bool Subsystem::subsystemTopologyHasBeenRealized() const {
    return getSubsystemGuts().subsystemTopologyHasBeenRealized();
}

inline void Subsystem::invalidateSubsystemTopologyCache() const {
    getSubsystemGuts().invalidateSubsystemTopologyCache(); // mutable
}

inline MeasureIndex Subsystem::adoptMeasure(AbstractMeasure& m)
{   return updSubsystemGuts().adoptMeasure(m); }
inline AbstractMeasure Subsystem::getMeasure(MeasureIndex mx) const
{   return getSubsystemGuts().getMeasure(mx); }


inline bool Subsystem::isInSystem() const 
{   return getSubsystemGuts().isInSystem(); }
inline bool Subsystem::isInSameSystem(const Subsystem& otherSubsystem) const 
{   return getSubsystemGuts().isInSameSystem(otherSubsystem); }

inline const System& Subsystem::getSystem() const 
{   return getSubsystemGuts().getSystem(); }
inline System& Subsystem::updSystem()       
{   return updSubsystemGuts().updSystem(); }
inline void Subsystem::setSystem(System& sys, SubsystemIndex id) 
{   updSubsystemGuts().setSystem(sys,id); }

inline bool Subsystem::isOwnerHandle() const 
{   return guts==nullptr || &guts->getOwnerSubsystemHandle()==this; }

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_SUBSYSTEM_GUTS_H_
