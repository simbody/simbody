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
Subsystem class can be found in the Subsystem::Guts class which is defined in a
separate header file. End users need not look over there. **/
class SimTK_SimTKCOMMON_EXPORT Subsystem {
public:
    class Guts; // local; name is Subsystem::Guts
    friend class Guts;

public:
    Subsystem() : guts(0) { } // an empty handle
    Subsystem(const Subsystem&);
    Subsystem& operator=(const Subsystem&);
    ~Subsystem();

    const String& getName()    const;
    const String& getVersion() const;

    // These call the corresponding State method, supplying this Subsystem's
    // SubsystemIndex. The returned indices are local to this Subsystem.
    QIndex allocateQ(State&, const Vector& qInit) const;
    UIndex allocateU(State&, const Vector& uInit) const;
    ZIndex allocateZ(State&, const Vector& zInit) const;

    DiscreteVariableIndex allocateDiscreteVariable
       (State&, Stage invalidates, AbstractValue* v) const;
    DiscreteVariableIndex allocateAutoUpdateDiscreteVariable
       (State&, Stage invalidates, AbstractValue* v, Stage updateDependsOn) const; 

    CacheEntryIndex allocateCacheEntry
       (const State&, Stage dependsOn, Stage computedBy, AbstractValue* v) const;
    CacheEntryIndex allocateCacheEntry   
       (const State& state, Stage g, AbstractValue* v) const 
    {   return allocateCacheEntry(state, g, g, v); }
    CacheEntryIndex allocateLazyCacheEntry   
       (const State& state, Stage earliest, AbstractValue* v) const 
    {   return allocateCacheEntry(state, earliest, Stage::Infinity, v); }

    QErrIndex allocateQErr         (const State&, int nqerr) const;
    UErrIndex allocateUErr         (const State&, int nuerr) const;
    UDotErrIndex allocateUDotErr      (const State&, int nudoterr) const;
    EventTriggerByStageIndex allocateEventTriggersByStage
       (const State&, Stage, int ntriggers) const;

    // These return views on State shared global resources. The views
    // are private to this subsystem, but the global resources themselves
    // are not allocated until the *System* advances to stage Model.
    // Note that there is no subsystem equivalent of the State's "y"
    // vector because in general a subsystem's state variables will
    // not be contiguous. However, a subsystem's q's, u's, and z's
    // will all be contiguous within those arrays.
    const Vector& getQ(const State&) const;
    const Vector& getU(const State&) const;
    const Vector& getZ(const State&) const;
    const Vector& getQDot(const State&) const;
    const Vector& getUDot(const State&) const;
    const Vector& getZDot(const State&) const;
    const Vector& getQDotDot(const State&) const;
    const Vector& getQErr(const State&) const;
    const Vector& getUErr(const State&) const;
    const Vector& getUDotErr(const State&) const;
    const Vector& getMultipliers(const State&) const;
    const Vector& getEventTriggersByStage(const State&, Stage) const;

    // These return writable access to this subsystem's partition in the
    // State pool of continuous variables. These can be called at Stage::Model
    // or higher, and if necesary they invalidate the Position (q), Velocity (u),
    // or Dynamics (z) stage respectively.
    Vector& updQ(State&) const; // invalidates Stage::Position
    Vector& updU(State&) const; // invalidates Stage::Velocity
    Vector& updZ(State&) const; // invalidates Stage::Dynamics

    // For convenience.
    void setQ(State& s, const Vector& q) const {
        assert(q.size() == getNQ(s));
        updQ(s) = q;
    }
    void setU(State& s, const Vector& u) const {
        assert(u.size() == getNU(s));
        updU(s) = u;
    }
    void setZ(State& s, const Vector& z) const {
        assert(z.size() == getNZ(s));
        updZ(s) = z;
    }

    // These update the State cache which is mutable; hence, const State. They
    // can be called only if the previous stage has already been realized, e.g.,
    // updQDot() is allowed only while realizing the Velocity stage, requiring
    // that Position stage has already been realized.
    Vector& updQDot(const State&) const;
    Vector& updUDot(const State&) const;
    Vector& updZDot(const State&) const;
    Vector& updQDotDot(const State&) const;
    Vector& updQErr(const State&) const;
    Vector& updUErr(const State&) const;
    Vector& updUDotErr(const State&) const;
    Vector& updMultipliers(const State&) const;
    Vector& updEventTriggersByStage(const State&, Stage) const;

    // These pull out the State entries which belong exclusively to
    // this Subsystem. These variables and cache entries are available
    // as soon as this subsystem is at stage Model.
    Stage getStage(const State&) const;
    const AbstractValue& getDiscreteVariable(const State& s, DiscreteVariableIndex dx) const;

    Real getDiscreteVarLastUpdateTime(const State& s, DiscreteVariableIndex dx) const
    {   return s.getDiscreteVarLastUpdateTime(getMySubsystemIndex(),dx); }
    CacheEntryIndex getDiscreteVarUpdateIndex(const State& s, DiscreteVariableIndex dx) const
    {   return s.getDiscreteVarUpdateIndex(getMySubsystemIndex(),dx); }
    const AbstractValue& getDiscreteVarUpdateValue(const State& s, DiscreteVariableIndex dx) const
    {   return s.getDiscreteVarUpdateValue(getMySubsystemIndex(),dx); }
    AbstractValue& updDiscreteVarUpdateValue(const State& s, DiscreteVariableIndex dx) const
    {   return s.updDiscreteVarUpdateValue(getMySubsystemIndex(),dx); }
    bool isDiscreteVarUpdateValueRealized(const State& s, DiscreteVariableIndex dx) const
    {   return s.isDiscreteVarUpdateValueRealized(getMySubsystemIndex(),dx); }
    void markDiscreteVarUpdateValueRealized(const State& s, DiscreteVariableIndex dx) const
    {   return s.markDiscreteVarUpdateValueRealized(getMySubsystemIndex(),dx); }

    // State is *not* mutable here -- must have write access to change state variables.
    AbstractValue& updDiscreteVariable(State&, DiscreteVariableIndex) const;

    const AbstractValue& getCacheEntry(const State&, CacheEntryIndex) const;
    // State is mutable here.
    AbstractValue& updCacheEntry(const State&, CacheEntryIndex) const;

    bool isCacheValueRealized(const State&, CacheEntryIndex) const;
    void markCacheValueRealized(const State&, CacheEntryIndex) const;
    void markCacheValueNotRealized(const State&, CacheEntryIndex) const;

    // Dimensions. These are valid at System Stage::Model while access to the 
    // various arrays may have stricter requirements. Hence it is better to use
    // these routines than to get a reference to a Vector above and ask for 
    // its size().

    SystemQIndex getQStart      (const State&) const;
    int getNQ          (const State&) const;
    SystemUIndex getUStart      (const State&) const;
    int getNU          (const State&) const;
    SystemZIndex getZStart      (const State&) const;
    int getNZ          (const State&) const;
    SystemQErrIndex getQErrStart   (const State&) const;
    int getNQErr       (const State&) const;
    SystemUErrIndex getUErrStart   (const State&) const;
    int getNUErr       (const State&) const;
    SystemUDotErrIndex getUDotErrStart(const State&) const;
    int getNUDotErr    (const State&) const;
    SystemMultiplierIndex getMultipliersStart (const State&) const;
    int getNMultipliers     (const State&) const;
    SystemEventTriggerByStageIndex getEventTriggerStartByStage(const State&, Stage) const;
    int getNEventTriggersByStage   (const State&, Stage) const;

	bool isInSystem() const;
	bool isInSameSystem(const Subsystem& otherSubsystem) const;

	const System& getSystem() const;
	System&       updSystem();

	SubsystemIndex getMySubsystemIndex() const;

    // Is this handle the owner of this rep? This is true if the
    // handle is empty or if its rep points back here.
    bool isOwnerHandle() const;
    bool isEmptyHandle() const;

    // There can be multiple handles on the same Subsystem.
    bool isSameSubsystem(const Subsystem& otherSubsystem) const;

    bool subsystemTopologyHasBeenRealized() const;
    void invalidateSubsystemTopologyCache() const;

    // Add a new Measure to this Subsystem. This method is generally used by Measure
    // constructors to install a newly-constructed Measure into its Subsystem.
    MeasureIndex adoptMeasure(AbstractMeasure&);

    AbstractMeasure getMeasure(MeasureIndex) const;
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
    void setSystem(System&, SubsystemIndex);

    explicit Subsystem(Subsystem::Guts* g) : guts(g) { }
    bool hasGuts() const {return guts!=0;}

private:
    // This is the only data member in this class. Also, any class derived from
    // Subsystem must have *NO* data members at all (data goes in the Guts 
    // class).
    Guts* guts;
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_SUBSYSTEM_H_
