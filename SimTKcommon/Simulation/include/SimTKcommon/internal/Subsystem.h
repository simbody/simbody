#ifndef SimTK_SimTKCOMMON_SUBSYSTEM_H_
#define SimTK_SimTKCOMMON_SUBSYSTEM_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-14 Stanford University and the Authors.        *
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
#include "SimTKcommon/internal/Measure.h"

#include <cassert>

namespace SimTK {

class System;

/** A %Subsystem is expected to be part of a larger System and to have
interdependencies with other subsystems of that same %System. It must NOT have 
dependencies on objects which are outside the %System. Consequently construction
of any concrete %Subsystem requires specification of a %System at that time.
Subsystems go through an extended construction phase in which their contents and
interdependencies are created. Thus all of a System's Subsystems generally need
to be available simultaneously during construction, so that they can reference 
each other.

There are three distinct users of this class:
   - the System class
   - the concrete Subsystems derived from this class
   - the end user of a concrete Subsystem

Only end user methods are public here. Methods intended for use by the concrete
Subsystem class implementation can be found in the Subsystem::Guts class which 
is defined in a separate header file. End users need not look over there. **/
class SimTK_SimTKCOMMON_EXPORT Subsystem {
public:
class Guts; // local; name is Subsystem::Guts
friend class Guts;

/** Default constructor creates and empty handle with a null Subsystem::Guts
pointer. **/
Subsystem() : guts(0) {}

/** Copy constructor clones the Subsystem::Guts object if there is one and
makes this the owner handle of the new clone. This is typically not very 
useful. **/
Subsystem(const Subsystem&);
/** Copy assignment deletes the Subsystem::Guts object if there is one and then
behaves like the copy constructor. Probably not useful in most cases. **/
Subsystem& operator=(const Subsystem&);

/** Destructor deletes the referenced Subsystem::Guts object if this is the
owner handle of that object, otherwise does nothing. Note that Subsystem::Guts
objects are not reference counted so any other handles pointing to the same
object will be invalid after the owner handle is destructed. **/
~Subsystem();


/** @name                   State access methods
These convenience methods are inline pass-throughs to the State methods of the 
same name but insert this %Subsystem's SubsystemIndex as the first argument. 
That is the value returned by the getMySubsystemIndex() method. An exception 
will be thrown if this %Subsystem is not contained in a System.  

See the SimTK::State documentation for the meaning of these methods; the 
behavior is identical here. **/
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
   (const State& s, Stage dependsOn, Stage computedBy, AbstractValue* v) const 
{   return s.allocateCacheEntry
               (getMySubsystemIndex(), dependsOn, computedBy, v); }

CacheEntryIndex allocateCacheEntry
    (const State& state, Stage g, AbstractValue* v) const 
{   return allocateCacheEntry(state, g, g, v); }
CacheEntryIndex allocateLazyCacheEntry   
    (const State& state, Stage earliest, AbstractValue* v) const 
{   return allocateCacheEntry(state, earliest, Stage::Infinity, v); }

QErrIndex allocateQErr(const State& s, int nqerr) const 
{   return s.allocateQErr(getMySubsystemIndex(), nqerr); }
UErrIndex allocateUErr(const State& s, int nuerr) const 
{   return s.allocateUErr(getMySubsystemIndex(), nuerr); }
UDotErrIndex allocateUDotErr(const State& s, int nudoterr) const 
{   return s.allocateUDotErr(getMySubsystemIndex(), nudoterr); }
EventTriggerByStageIndex 
allocateEventTriggersByStage(const State& s, Stage g, int ntriggers) const 
{   return s.allocateEventTrigger(getMySubsystemIndex(),g,ntriggers); }

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
const Vector& getEventTriggersByStage(const State& s, Stage g) const
{   return s.getEventTriggersByStage(getMySubsystemIndex(),g); }

Vector& updQErr(const State& s) const 
{   return s.updQErr(getMySubsystemIndex()); }
Vector& updUErr(const State& s) const 
{   return s.updUErr(getMySubsystemIndex()); }
Vector& updUDotErr(const State& s) const 
{   return s.updUDotErr(getMySubsystemIndex()); }
Vector& updMultipliers(const State& s) const 
{   return s.updMultipliers(getMySubsystemIndex()); }
Vector& updEventTriggersByStage(const State& s, Stage g) const
{   return s.updEventTriggersByStage(getMySubsystemIndex(),g); }

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

SystemEventTriggerByStageIndex getEventTriggerStartByStage(const State& s, Stage g) const 
{   return s.getEventTriggerStartByStage(getMySubsystemIndex(),g); }
int getNEventTriggersByStage   (const State& s, Stage g) const 
{   return s.getNEventTriggersByStage(getMySubsystemIndex(),g); }


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
void markCacheValueRealized(const State& s, CacheEntryIndex cx) const 
{   s.markCacheValueRealized(getMySubsystemIndex(), cx); }
void markCacheValueNotRealized(const State& s, CacheEntryIndex cx) const 
{   s.markCacheValueNotRealized(getMySubsystemIndex(), cx); }
/**@}**/

/** @name                  Miscellaneous bookkeeping
These methods are not commonly used by end users. They are concerned with
the mechanics of creating and managing a %Subsystem, typically a concern of
%Subsystem developers. These are mostly inline pass-throughs to the 
Subsystem::Guts object pointed to by this handle. **/
/**@{**/

/** Obtain the Subsystem name if one was given on construction of the concrete
Subsystem. Simbody does not interpret this in any way. **/
inline const String& getName()    const;
/** Obtain the Subsystem version string if one was given on construction. 
Simbody does not interpret this in any way. **/
inline const String& getVersion() const;

/** Return \c true if this %Subsystem is contained in a System. **/
inline bool isInSystem() const;
/** Return \c true if this %Subsystem is contained in the same System as 
contains the given \a otherSubsystem. Returns \c false if either %Subsystem is 
not contained in a %System, or if the %Systems don't match. **/
inline bool isInSameSystem(const Subsystem& otherSubsystem) const;

/** Return a const reference to the System that contains this %Subsystem. Throws
an exception if this is not contained in a System; call isInSystem() first if 
you aren't sure. **/
inline const System& getSystem() const;
/** Return a writable reference to the System that contains this %Subsystem. 
Throws an exception if this is not contained in a System; call isInSystem() 
first if you aren't sure. **/
inline System& updSystem();
/** Inform this %Subsystem of the System that contains it, as well as the
SubsystemIndex which the System has assigned to it. **/
inline void setSystem(System& system, SubsystemIndex subx);

/** Return the SubsystemIndex within the containing System. An exception will
be thrown (in Debug mode at least) if this Subsystem is not in a System. **/
inline SubsystemIndex getMySubsystemIndex() const;

/** Return true if this handle has a null Subsystem::Guts pointer. **/
inline bool isEmptyHandle() const {return guts==0;}

/** Determine if \c this Subsystem handle refers to the same Subsystem::Guts
object as handle \a otherSubsystem. There can be multiple handles on the same
underlying Subsystem::Guts object but at most one can be the owner. **/
inline bool isSameSubsystem(const Subsystem& otherSubsystem) const
{   return guts && (guts==otherSubsystem.guts); }

/** Is this Subsystem handle the owner of the Subsystem::Guts object it points
to? This is \c true if the handle is empty or if its Guts object points back 
to this handle. **/
inline bool isOwnerHandle() const; 

/** Returns \c true if this %Subsystem's realizeTopology() method has been
called since the last topological change or call to 
invalidateSubsystemTopologyCache(). **/
inline bool subsystemTopologyHasBeenRealized() const;
/** Always call this method when a topological change is made to this 
%Subsystem to indicate that any Stage::Topology cache values may need
recomputation. If the %Subsystem belongs to a System, the %System's overall 
topology will also be invalidated, since its Stage::Topology cannot be valid
if any of its %Subsystem topology stages are not valid. However, the stages
of other %Subsystems in the same %System are not affected. A subsequent call
to realizeTopology() is required before any computed values may be 
obtained. **/
inline void invalidateSubsystemTopologyCache() const;

// Add a new Measure to this Subsystem. This method is generally used by Measure
// constructors to install a newly-constructed Measure into its Subsystem.
inline MeasureIndex adoptMeasure(AbstractMeasure&);
inline AbstractMeasure getMeasure(MeasureIndex) const;
template <class T> Measure_<T> getMeasure_(MeasureIndex mx) const
{   return Measure_<T>::getAs(getMeasure(mx));}

// dynamic_cast the returned reference to a reference to your concrete Guts
// class.
const Subsystem::Guts& getSubsystemGuts() const {assert(guts); return *guts;}
Subsystem::Guts&       updSubsystemGuts()       {assert(guts); return *guts;}

// Put new Guts into this *empty* handle and take over ownership.
// If this handle is already in use, this routine will throw
// an exception.
void adoptSubsystemGuts(Subsystem::Guts* g);

explicit Subsystem(Subsystem::Guts* g) : guts(g) { }
bool hasGuts() const {return guts!=0;}
/**@}**/

private:
// This is the only data member in this class. Also, any class derived from
// Subsystem must have *NO* data members at all (data goes in the Guts 
// class).
Guts* guts;
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_SUBSYSTEM_H_
