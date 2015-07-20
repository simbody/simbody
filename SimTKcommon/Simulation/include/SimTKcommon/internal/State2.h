#ifndef SimTK_SimTKCOMMON_STATE2_H_
#define SimTK_SimTKCOMMON_STATE2_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2015 Stanford University and the Authors.           *
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

/** @file
TODO: This is currently just a sketch for an alternate State implementation.
Declares the user-visible part of a SimTK::State, the implementation is
done in a separate internal class. **/

#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"

#include <ostream>
#include <cassert>
#include <algorithm>
#include <memory>
#include <cstdint>

namespace SimTK {

SimTK_DEFINE_UNIQUE_INDEX_TYPE(VariableIndex);
SimTK_DEFINE_UNIQUE_INDEX_TYPE(CacheIndex);

class State2;


enum class ContinuousVarPoolIndex {
    Q  = 0,
    U  = 1,
    Z  = 2
};


//==============================================================================
//                                VARIABLE
//==============================================================================
/** An abstract state variable may be discrete, continuous, or a pool of 
continuous variables.

A variable's lifetime is determined either by an "allocation stage" below which
the variable must be deleted (and that must be Topology or Model stage), or by 
a "parent variable" whose deletion requires deletion of its children. Ultimately
some parent or ancestor has an allocation stage instead of a parent.

Every variable has a value. For continuous variables that value is a Vector of
Reals. For discrete variables the value is of an arbitrary type.
When a variable's value is modified (or more accurately, when its value is
returned with write access), several things happen:
- the "time of last modification" is set to the current time
- an integer "value version" is incremented 
- the "invalidated stage" (if any) and all higher stages are invalidated, for 
  the State as a whole and for each of its subsystems
- a specific set of cache entry dependents is invalidated

Modifying a pool's value directly is also a modification of each of its members.
Modifying a variable that is a pool member is also considered a modification
to the pool variable, but not to the other members of the pool. This allows
cache entries to depend either on a whole pool, such as any q, or to just
a particular member of that pool.
**/
class Variable {
public:
    Variable() 
    :   m_timeLastModified(dNaN), m_valueVersion(-1) {}


    double getTimeLastModified() const {return m_timeLastModified;}
    int64_t getValueVersion() const {return m_valueVersion;}

    virtual ~Variable() {}
    virtual Variable* clone() const = 0;

    /** Note that this variable has been modified, and invalidate its
    dependents. This is invoked automatically whenever write access to this
    variable's value is obtained. **/
    inline void markVariableModified(State2& state, const double& t);

private:
friend class State2;

    void deleteVariable(State2& state, VariableIndex);

    //void completeAdoption(State2& state, VariableIndex myIndex) {
    //    SimTK_ASSERT1_ALWAYS(!m_myIndex.isValid(),
    //        "Variable::completeAdoption(): Variable %d already adopted.",
    //        (int)m_myIndex);
    //    SimTK_ASSERT1_ALWAYS(&state.getVariable(myIndex) == this,
    //        "Variable::completeAdoption(): Variable not found at index %d.",
    //        (int)myIndex);
    //    m_myIndex = myIndex;

    //    // Register as child with parent variable, if any.
    //    if (m_parent.isValid()) {
    //        state.updVariable(m_parent).addChild(m_myIndex);
    //    }
    //}




    //virtual void deleteVariableVirtual(State2& state) = 0;
    //virtual void completeAdoptionVirtual(State2& state) = 0;


    VariableIndex           m_parent;
    Array_<VariableIndex>   m_children;

    Stage                   m_allocationStage;
    Stage                   m_invalidatedStage;
    Array_<CacheIndex>      m_dependents;

    VariableIndex           m_myIndex; // assigned by State

    // Run time information.
    double                  m_timeLastModified;
    int64_t                 m_valueVersion;


};

//==============================================================================
//                              TIME VAR
//==============================================================================
/** This variable holds the current time for a State. Time is always
a double even when Real is float. **/
class TimeVar : public Variable {
public:
    TimeVar() : m_time(dNaN) {}
    TimeVar* clone() const override {return new TimeVar(*this);}

    double getTime() const {return m_time;}
    double& updTime(State2& state) {
        //assert(&state.getTimeVariable() == this);
        markVariableModified(state, m_time);
        return m_time;
    }
private:
    double  m_time;
};

//==============================================================================
//                              CONTINUOUS VAR
//==============================================================================
/** A continuous variable has a scalar or vector value that is maintained at
runtime in a pool shared with other similar continuous variables. An initial
value is stored here, and whenever the pool is invalidated the current value
is transferred back from the pool and becomes the new initial value.

All continuous variables are either position kinematic variables q, velocity
kinematic variables u, or general first order variables z. Consequently the
State maintains three pools, Q, U, and Z that gather all the values for the
corresponding variables. The pools are also considered variables and accessing
a pool for write is considered a change to the pool and to *all* its members. 
Accessing any member for write is a change to that variable and to its pool but
not to any other members of the pool.

The invalidated stage is the same for a pool and all its members, but cache 
entries may be dependencies of just a variable or the whole pool. **/
class ContinuousVar : public Variable {
public:
    explicit ContinuousVar(const Vector& initialValue)
    :   m_savedValue(initialValue) {}

    ContinuousVar* clone() const override 
    {   return new ContinuousVar(*this); }

private:
    Vector          m_savedValue;
    VariableIndex   m_pool;
    Vector          m_pooledValue;
};

//==============================================================================
//                           CONTINUOUS VAR POOL
//==============================================================================
/** Once allocated, this variable contains the current values for a collection
of related continuous variables. For example, all the z variables from every
subsystem are gathered in a Z pool. The pool is invalidated by the addition
or removal of any ContinuousVar that says it is a member of this pool. In that
case values are first copied from the pool back to the constituent variables,
then the pool is marked invalid. It is resized, initialized, and validated
when referenced for the first time.
**/
class ContinuousVarPool : public Variable {
public:
    ContinuousVarPool* clone() const override 
    {   return new ContinuousVarPool(*this); }

private:
    Array_<VariableIndex> m_members;
    Vector                m_value;
};

//==============================================================================
//                               DISCRETE VAR 
//==============================================================================
class DiscreteVar : public Variable {
public:
    DiscreteVar* clone() const override 
    {   return new DiscreteVar(*this); }

    const AbstractValue& getAbstractValue() const {return *m_value;}

    template <class T>
    const T& getValue() const {return Value<T>::downcast(*m_value).get();}

    template <class T>
    inline T& updValue(State2& state);

private:
    Array_<VariableIndex>   m_childVariables;
    Array_<CacheIndex>      m_childCacheEntries;

    CacheIndex              m_autoUpdateEntry; // if any
    ClonePtr<AbstractValue> m_value;
};

//==============================================================================
//                                CACHE ENTRY
//==============================================================================
/** An entry in the State's cache of stored values whose validity is managed
automatically so that a stale value cannot be returned. 

If a cache entry has a computedBy stage, it is assumed to be valid if the
subsystem that owns it has reached that stage. No other checks will be 
performed in that case.

If there is no computedBy stage (lazy cache entry), or computedBy has not
yet been reached by the owner subsystem, the value *may* still be valid,
provided all three of these conditions have been met:
- the subsystem has reached the cache entry's dependsOn stage, AND
- the subsystem's dependsOn stage version number has not changed since
  the cache entry's value was last computed, AND
- the prerequisiteHasChanged flag is false.

The prerequisiteHasChanged flag is cleared, and the dependsOn stage version
recorded, whenever the cache entry is explicitly marked as up to date.
The prerequisiteHasChanged flag is set whenever one of a designated set of
prerequisite state variables and other cache entries is modified. It is *not* 
affected by changes to the subsystem's stage, so the value can be invalid even
if no prerequisite has changed.

Note that this is optimized for fast repeated access once evaluated; validity
can be determined in constant time no matter how many prerequisites there are.
Instead, the cost of invalidation is mostly borne at the time any prerequisite 
variable or cache entry is modified.
**/
class Cache {
public:
    Cache() 
    :   m_timeLastUpdated(dNaN), m_valueVersion(-1),
        m_lastDependsOnStageVersion(-1) {}

    /** See class documentation for conditions under which a cache entry's
    value is considered to be up to date. **/
    inline bool isValueUpToDate(const State2& state) const;

    /** Mark this cache entry invalid if it is currently valid. This increments
    the value version number and sets the invalid flag. It does not change the 
    stored value or timeLastUpdated, and will not affect access to the cache 
    entry if it has a computedBy stage and its subsystem is already past that 
    stage. **/
    inline void noteThatPrerequisiteHasChanged(State2& state);

    /** After the cache value has been computed, call this method to note that
    the value is current as of time `t`. The dependsOn stage version of the
    owner subsystem is recorded, and the prerequisiteHasChanged flag is 
    cleared. **/
    inline void markValueUpToDate(const State2& state, const double& t);

    // Not abstract.
    Cache* clone() const {return new Cache(*this);}

private:
friend class State2;
    SubsystemIndex          m_owner;
    VariableIndex           m_parent;
    Stage                   m_allocationStage;
    Stage                   m_dependsOnStage;   // Empty if not stage-dependent
    Stage                   m_computedByStage;  // Infinity if never

    VariableIndex           m_autoUpdateVariable; // if any
    Array_<VariableIndex>   m_prerequisites;

    Array_<CacheIndex>      m_dependents;

    CacheIndex              m_myIndex; // set when added to State

    // Run time information.
    double                  m_timeLastUpdated;
    int64_t                 m_valueVersion;
    ClonePtr<AbstractValue> m_value;
    int64_t                 m_lastDependsOnStageVersion; // when last updated
    bool                    m_prerequisiteHasChanged;    // since last updated
};

//==============================================================================
//                             SUBSYSTEM INFO
//==============================================================================
class SubsystemInfo {
public:
    SubsystemInfo() {
        for (Stage s = Stage::LowestValid; s <= Stage::HighestValid; ++s)
            m_stageVersions[s] = -1;
    }
    Stage getCurrentStage() const {return m_currentStage;}
    int64_t getStageVersion(Stage stage) const {
        return m_stageVersions[stage];
    }

private:
    mutable Stage        m_currentStage;
    mutable int64_t      m_stageVersions[Stage::NValid];
};


//==============================================================================
//                                STATE
//==============================================================================
/** Container for all of a Simbody System's state variables and cached 
results that depend on those variables.

Variables come in three flavors: time t, continuous variables y, and discrete
variables d. Continuous variables include the kinematic variables q and u, and
arbitrary variables z.

The State understands that a System is composed of Subsystems and maintains 
per-subsystem information. Each variable is owned by a particular subsystem,
although continuous variables can be pooled across subsystems for convenient
access.


**/
class /*SimTK_SimTKCOMMON_EXPORT*/ State2 {
public:
    State2() {
        m_variables.push_back(CloneOnWritePtr<TimeVar>(new TimeVar()));
        auto time  = new TimeVar();
        auto qpool = new ContinuousVarPool();
        auto upool = new ContinuousVarPool();
        auto zpool = new ContinuousVarPool();
        m_time  = adoptVariable(time);
        m_qPool = adoptVariable(qpool);
        m_uPool = adoptVariable(upool);
        m_zPool = adoptVariable(zpool);
    }


    ~State2() {}

    void invalidateStage(Stage stage);

    const SubsystemInfo& getSubsystemInfo(SubsystemIndex sx) const {
        return m_subsystemInfo[sx];
    }


    VariableIndex adoptVariable(Variable* var) { 
        const VariableIndex vx(takeFreeVarSlot());
        m_variables.push_back(CloneOnWritePtr<Variable>(var)); // take ownership
        //m_variables.back()->completeAdoption(*this, vx);
        return vx;
    }

    CacheIndex adoptCacheEntry(Cache* entry) {
        const CacheIndex cx(takeFreeCacheSlot());
        m_cacheEntries.push_back(CloneOnWritePtr<Cache>(entry)); // own
        //m_cacheEntries.back()->completeAdoption(*this, cx);
        return cx;
    }

    const Variable& getVariable(VariableIndex i) const 
    {   return *m_variables[i]; }
    Variable& updVariable(VariableIndex i) 
    {   return *m_variables[i]; }
    const Cache& getCacheEntry(CacheIndex i) const 
    {   return *m_cacheEntries[i]; }
    Cache& updCacheEntry(CacheIndex i) const 
    {   return *m_cacheEntries[i]; }

    /** Return current value of time; will be NaN in an uninitialized 
    %State. **/ 
    double getTime() const {
        return static_cast<const TimeVar*>
                    (m_variables[m_time].get())->getTime();
    }

    double& updTime() {
        return static_cast<TimeVar*>
                    (m_variables[m_time].upd())->updTime(*this);
    }

    void setTime(double t) {
        updTime() = t;
    }

    VariableIndex getTimeIndex()  const {return m_time;}
    VariableIndex getQPoolIndex() const {return m_qPool;}
    VariableIndex getUPoolIndex() const {return m_uPool;}
    VariableIndex getZPoolIndex() const {return m_zPool;}

    /** Get the pooled q variables from all Subsystems. Valid after the pool
    has been constructed at Model stage. **/
    const Vector& getQ() const;

    /** Get the pooled u variables from all Subsystems. Valid after the pool
    has been constructed at Model stage. **/
    const Vector& getU() const;

    /** Get the pooled z variables from all Subsystems, regardless of 
    their invalidates stage (that is, z=zp zv zd za zr packed together). **/
    const Vector& getZ() const;

    /** Get the pooled variable for one of zp zv zd za zr, including 
    contributions from all the subsystems packed together. **/
    const Vector& getZ(Stage stage) const;

    /** Access pooled Q for update; invalidates Position stage. **/
    Vector& updQ();

    /** Access pooled U for update; invalidates Velocity stage. **/
    Vector& updU();

    /** Access pooled Z for update; invalidates lowest Z stage. **/
    Vector& updZ();

private:   
    // Find an available slot in the Variables array and return its index.
    // The contents of that slot will be null.
    VariableIndex takeFreeVarSlot() {
        if (!m_freeVarSlots.empty()) {
            const VariableIndex vx = m_freeVarSlots.back();
            m_freeVarSlots.pop_back();
            return vx;
        }
        const VariableIndex vx(m_variables.size());
        m_variables.push_back(); // create a null entry at the end
        return vx;
    }

    // Find an available slot in the Cache array and return its index.
    // The contents of that slot will be null.
    CacheIndex takeFreeCacheSlot() {
        if (!m_freeCacheSlots.empty()) {
            const CacheIndex cx = m_freeCacheSlots.back();
            m_freeCacheSlots.pop_back();
            return cx;
        }
        const CacheIndex cx(m_cacheEntries.size());
        m_cacheEntries.push_back(); // create a null entry at the end
        return cx;
    }


private:
    Array_<SubsystemInfo, SubsystemIndex>               m_subsystemInfo;

    // All variables, with free slots null and tracked.
    Array_<CloneOnWritePtr<Variable>, VariableIndex>    m_variables;
    Array_<VariableIndex>                               m_freeVarSlots;

    VariableIndex                                       m_time;
    VariableIndex                                       m_qPool;
    VariableIndex                                       m_uPool;
    VariableIndex                                       m_zPool;

    Array_<VariableIndex>                               m_autoUpdateVars;


    // All cache entries, with free slots null and tracked. 
    mutable Array_<CloneOnWritePtr<Cache>, CacheIndex>  m_cacheEntries;
    mutable Array_<CacheIndex>                          m_freeCacheSlots;


    mutable Stage       m_currentSystemStage;
    mutable int64_t     m_systemStageVersions[Stage::NValid];
};


//------------------------------------------------------------------------------
//                        VARIABLE IMPLEMENTATION
//------------------------------------------------------------------------------
inline void Variable::markVariableModified(State2& state, const double& t) {
    assert(state.getTime() == t);
    m_timeLastModified = t;
    ++m_valueVersion;
    state.invalidateStage(m_invalidatedStage);
    for (CacheIndex cx : m_dependents)
        state.updCacheEntry(cx).noteThatPrerequisiteHasChanged(state);
}

//------------------------------------------------------------------------------
//                        DISCRETE VAR IMPLEMENTATION
//------------------------------------------------------------------------------

template <class T>
inline T& DiscreteVar::updValue(State2& state) {
    markVariableModified(state, state.getTime());
    return Value<T>::updDowncast(*m_value).upd();
}

//------------------------------------------------------------------------------
//                          CACHE IMPLEMENTATION
//------------------------------------------------------------------------------
// This gets called *a lot*!
SimTK_FORCE_INLINE bool Cache::isValueUpToDate(const State2& state) const {
    const SubsystemInfo& ss = state.getSubsystemInfo(m_owner);
    const Stage currentStage = ss.getCurrentStage();
    if (currentStage >= m_computedByStage)
        return true;
    if (currentStage < m_dependsOnStage)
        return false;
    // We're between dependsOn and computedBy; lazy evaluation territory.
    if (m_prerequisiteHasChanged)
        return false;
    const int64_t dependsOnVersion = ss.getStageVersion(m_dependsOnStage);
    return m_lastDependsOnStageVersion == dependsOnVersion;
}

inline void Cache::noteThatPrerequisiteHasChanged(State2& state) {
    if (!m_prerequisiteHasChanged) {
        m_prerequisiteHasChanged = true;
        ++m_valueVersion; // the version it will have when next updated
        for (CacheIndex cx : m_dependents)
            state.updCacheEntry(cx).noteThatPrerequisiteHasChanged(state);
    }
}

inline void Cache::markValueUpToDate(const State2& state, const double& t) {
    const SubsystemInfo& ss = state.getSubsystemInfo(m_owner);
    m_timeLastUpdated = t;
    m_lastDependsOnStageVersion = ss.getStageVersion(m_dependsOnStage);
    m_prerequisiteHasChanged = false;
}

} // namespace SimTK


#endif // SimTK_SimTKCOMMON_STATE2_H_
