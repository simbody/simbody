#ifndef SimTK_SimTKCOMMON_STATE_IMPL_H_
#define SimTK_SimTKCOMMON_STATE_IMPL_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-14 Stanford University and the Authors.        *
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

/** @file
This is part of the internal implementation of SimTK::State and does not
contain any user-visible objects. **/

/** @cond **/   // Hide from Doxygen

// This header is logically an extension of State.h and is intended to be
// included from within State.h and nowhere else.

namespace SimTK {

// These local classes
//      DiscreteVarInfo
//      CacheVarInfo
//      EventInfo
// contain the information needed for each discrete variable and cache entry,
// including a reference to its current value and everything needed to 
// understand its dependencies, as well as information for each Event allocated
// by a Subsystem.
//
// These are intended as elements in an allocation stack as described above,
// so it is expected that they will get reaallocated, copied, and destructed 
// during allocation as the Array_ gets resized from time to time. However, the
// discrete variable and cache entry *values* must remain in a fixed location in 
// memory once allocated, because callers are permitted to retain references to
// these values once they have been allocated. So pointers to the AbstractValues
// are kept in these objects, and only the pointers are copied around. That 
// means the actual value object will not be deleted by the destructor; be sure
// to do that explicitly in the higher-level destructor or you'll have a nasty
// leak.

//==============================================================================
//                           DISCRETE VAR INFO
//==============================================================================
class DiscreteVarInfo {
public:
    DiscreteVarInfo()
    :   allocationStage(Stage::Empty), invalidatedStage(Stage::Empty),
        value(0), timeLastUpdated(NaN) {}

    DiscreteVarInfo(Stage allocation, Stage invalidated, AbstractValue* v)
    :   allocationStage(allocation), invalidatedStage(invalidated),
        autoUpdateEntry(), value(v), timeLastUpdated(NaN) 
    {   assert(isReasonable()); }

    // Default copy constructor, copy assignment, destructor are shallow.

    // Use this to make this entry contain a *copy* of the source value.
    // If the destination already has a value, the new value must be
    // assignment compatible.
    DiscreteVarInfo& deepAssign(const DiscreteVarInfo& src) {
        assert(src.isReasonable());
         
        allocationStage   = src.allocationStage;
        invalidatedStage  = src.invalidatedStage;
        autoUpdateEntry   = src.autoUpdateEntry;
        if (value) *value = *src.value;
        else        value = src.value->clone();
        timeLastUpdated   = src.timeLastUpdated;
        return *this;
    }

    // For use in the containing class's destructor.
    void deepDestruct() {delete value; value=0;}
    const Stage& getAllocationStage()  const {return allocationStage;}

    // Exchange value pointers (should be from this dv's update cache entry).
    void swapValue(Real updTime, AbstractValue*& other) 
    {   std::swap(value, other); timeLastUpdated=updTime; }

    const AbstractValue& getValue() const {assert(value); return *value;}
    Real                 getTimeLastUpdated() const 
    {   assert(value); return timeLastUpdated; }
    AbstractValue&       updValue(Real updTime)
    {   assert(value); timeLastUpdated=updTime; return *value; }

    const Stage&    getInvalidatedStage() const {return invalidatedStage;}
    CacheEntryIndex getAutoUpdateEntry()  const {return autoUpdateEntry;}
    void setAutoUpdateEntry(CacheEntryIndex cx) {autoUpdateEntry = cx;}

private:
    // These are fixed at construction.
    Stage           allocationStage;
    Stage           invalidatedStage;
    CacheEntryIndex autoUpdateEntry;

    // These change at run time.
    AbstractValue*  value;
    Real            timeLastUpdated;

    bool isReasonable() const
    {    return (allocationStage==Stage::Topology 
                 || allocationStage==Stage::Model)
             && (invalidatedStage > allocationStage)
             && (value != 0); }
};



//==============================================================================
//                            CACHE ENTRY INFO
//==============================================================================
class CacheEntryInfo {
public:
    CacheEntryInfo()
    :   allocationStage(Stage::Empty), dependsOnStage(Stage::Empty), 
        computedByStage(Stage::Empty), value(0), versionWhenLastComputed(-1) {}

    CacheEntryInfo
       (Stage allocation, Stage dependsOn, Stage computedBy, AbstractValue* v)
    :   allocationStage(allocation), dependsOnStage(dependsOn), 
        computedByStage(computedBy), value(v), versionWhenLastComputed(0) 
    {   assert(isReasonable()); }

    bool isCurrent(const Stage& current, const StageVersion versions[]) const 
    {   if (current >= computedByStage) return true;
        if (current <  dependsOnStage)  return false;
        return versions[dependsOnStage] == versionWhenLastComputed;}

    StageVersion getVersionWhenLastComputed() const 
    {   return versionWhenLastComputed; }

    // These affect only the explicit "last computed" flag which does not fully
    // determine whether the value is current; see isCurrent() above.
    void invalidate() {versionWhenLastComputed = 0;}
    void markAsComputed(const StageVersion versions[])
    {   versionWhenLastComputed = versions[dependsOnStage];}

    // Default copy constructor, copy assignment, destructor are shallow.

    // Use this to make this entry contain a *copy* of the source value.
    CacheEntryInfo& deepAssign(const CacheEntryInfo& src) {
        assert(src.isReasonable());

        allocationStage   = src.allocationStage;
        dependsOnStage    = src.dependsOnStage;
        computedByStage   = src.computedByStage;
        associatedVar     = src.associatedVar;
        if (value) *value = *src.value;
        else        value = src.value->clone();
        versionWhenLastComputed = src.versionWhenLastComputed;
        return *this;
    }

    // For use in the containing class's destructor.
    void deepDestruct() {delete value; value=0;}
    const Stage& getAllocationStage() const {return allocationStage;}

    // Exchange values with a discrete variable (presumably this
    // cache entry has been determined to be that variable's update
    // entry but we're not checking here).
    void swapValue(Real updTime, DiscreteVarInfo& dv) 
    {   dv.swapValue(updTime, value); }
    const AbstractValue& getValue() const {assert(value); return *value;}
    AbstractValue&       updValue()       {assert(value); return *value;}

    const Stage&          getDependsOnStage()  const {return dependsOnStage;}
    const Stage&          getComputedByStage() const {return computedByStage;}
    DiscreteVariableIndex getAssociatedVar()   const {return associatedVar;}
    void setAssociatedVar(DiscreteVariableIndex dx)  {associatedVar=dx;}
private:
    // These are fixed at construction.
    Stage                   allocationStage;
    Stage                   dependsOnStage;
    Stage                   computedByStage;
    DiscreteVariableIndex   associatedVar;  // if this is an auto-update entry

    // These change at run time.
    AbstractValue*          value;
    StageVersion            versionWhenLastComputed;//version of Stage dependsOn

    bool isReasonable() const
    {    return (   allocationStage==Stage::Topology
                 || allocationStage==Stage::Model
                 || allocationStage==Stage::Instance)
             && (computedByStage >= dependsOnStage)
             && (value != 0)
             && (versionWhenLastComputed >= 0); }
};



//==============================================================================
//                               TRIGGER INFO
//==============================================================================
class TriggerInfo {
public:
    TriggerInfo() 
    :   allocationStage(Stage::Empty), firstIndex(-1), nslots(0) {}

    TriggerInfo(Stage allocation, int index, int n)
    :   allocationStage(allocation), firstIndex(index), nslots(n)
    {   assert(isReasonable()); assert(n>0);}

    // Default copy constructor, copy assignment, destructor are fine since 
    // there is no heap object owned here.

    int getFirstIndex() const {return firstIndex;}
    int getNumSlots() const {return nslots;}

    // These the the "virtual" methods required by template methods elsewhere.
    TriggerInfo& deepAssign(const TriggerInfo& src) 
    {   return operator=(src); }
    void         deepDestruct() {}
    const Stage& getAllocationStage() const {return allocationStage;}
private:
    // These are fixed at construction.
    Stage                   allocationStage;
    int                     firstIndex;
    int                     nslots;

    bool isReasonable() const {
        return (    allocationStage==Stage::Topology
                 || allocationStage==Stage::Model
                 || allocationStage==Stage::Instance);
    }
};



//==============================================================================
//                           CONTINUOUS VAR INFO
//==============================================================================
// Used for q, u, and z (in separate stacks).
// These accumulate default values for this subsystem's use of shared
// global state variables. After the System is advanced to Stage::Model,
// the state will allocate those globals and copy the initial
// values stored here into them. Some of these are allocated at Topology
// stage, and some at Model stage. If Model stage is invalidated, variables
// allocated then are forgotten while the ones allocated at Topology stage
// remain.
class ContinuousVarInfo {
public:
    ContinuousVarInfo() : allocationStage(Stage::Empty), firstIndex(-1) {}

    ContinuousVarInfo(Stage         allocation, 
                      int           index,  // QIndex, UIndex, or ZIndex
                      const Vector& initVals, 
                      const Vector& varWeights=Vector())
    :   allocationStage(allocation), firstIndex(index), initialValues(initVals)
    {   assert(isReasonable());
        assert(varWeights.size()==0 || varWeights.size()==initVals.size());
        assert(weightsAreOK(varWeights));
        if (varWeights.size()) weights=varWeights;
        else weights=Vector(initVals.size(), Real(1));
    }

    int getFirstIndex() const {return firstIndex;}
    int getNumVars() const {return initialValues.size();}
    const Vector& getInitialValues() const {return initialValues;}
    const Vector& getWeights() const {return weights;}

    // Default copy constructor, copy assignment, destructor are fine since 
    // there is no heap object owned here.

    // These the the "virtual" methods required by template methods elsewhere.
    ContinuousVarInfo& deepAssign(const ContinuousVarInfo& src) 
    {   return operator=(src); }
    void               deepDestruct() {}
    const Stage&       getAllocationStage() const {return allocationStage;}
private:
    // These are fixed at construction.
    Stage     allocationStage;
    int       firstIndex;   // a QIndex, UIndex, or ZIndex
    Vector    initialValues;
    Vector    weights; // only used for u and z

    static bool weightsAreOK(const Vector& wts) {
        for (int i=0; i<wts.size(); ++i)
            if (wts[i] <= 0) return false;
        return true;
    }

    bool isReasonable() const {
        return (    allocationStage==Stage::Topology
                 || allocationStage==Stage::Model);
    }
};

//==============================================================================
//                           CONSTRAINT ERR INFO
//==============================================================================
// Used for qerr, uerr, and udoterr.
class ConstraintErrInfo {
public:
    ConstraintErrInfo() : allocationStage(Stage::Empty), firstIndex(-1) {}

    ConstraintErrInfo(Stage         allocation,
                      int           index, // QErr, UErr, or UDotErrIndex
                      int           nerr,
                      const Vector& varWeights=Vector())
    :   allocationStage(allocation), firstIndex(index)
    {   assert(isReasonable());
        assert(varWeights.size()==0 || varWeights.size()==nerr);
        assert(weightsAreOK(varWeights));
        if (varWeights.size()) weights=varWeights;
        else weights=Vector(nerr, Real(1));
    }

    int getFirstIndex() const {return firstIndex;}
    int getNumErrs() const {return weights.size();}
    const Vector& getWeights() const {return weights;}

    // Default copy constructor, copy assignment, destructor are fine since 
    // there is no heap object owned here.

    // These the the "virtual" methods required by template methods elsewhere.
    ConstraintErrInfo& deepAssign(const ConstraintErrInfo& src) 
    {   return operator=(src); }
    void               deepDestruct() {}
    const Stage&       getAllocationStage() const {return allocationStage;}
private:
    // These are fixed at construction.
    Stage     allocationStage;
    int       firstIndex;   // a QErrIndex, UErrIndex, or UDotErrIndex
    Vector    weights;      // only used for u and z

    static bool weightsAreOK(const Vector& wts) {
        for (int i=0; i<wts.size(); ++i)
            if (wts[i] <= 0) return false;
        return true;
    }

    bool isReasonable() const {
        return (    allocationStage==Stage::Topology
                 || allocationStage==Stage::Model
                 || allocationStage==Stage::Instance);
    }
};



//==============================================================================
//                           PER SUBSYSTEM INFO
//==============================================================================
// This internal utility class is used to capture all the information needed for
// a single subsystem within the StateImpl.
class PerSubsystemInfo {
public:
    PerSubsystemInfo() : currentStage(Stage::Empty)     {initialize();}
    PerSubsystemInfo(const String& n, const String& v) 
      : name(n), version(v), currentStage(Stage::Empty) {initialize();}

    // Everything will properly clean itself up except for the AbstractValues
    // stored in the discrete variable and cache entry arrays. Be sure to
    // delete those prior to allowing the arrays themselves to destruct.
    ~PerSubsystemInfo() {
        clearAllStacks();
    }

    // Copy constructor copies all variables but cache only through
    // Instance stage. Note that this must be done in conjunction with
    // copying the whole state or our global resource indices will
    // be nonsense.
    PerSubsystemInfo(const PerSubsystemInfo& src) : currentStage(Stage::Empty) {
        initialize();
        copyFrom(src, Stage::Instance);
    }

    PerSubsystemInfo& operator=(const PerSubsystemInfo& src) {
        // destination is already initialized; copyFrom() will try
        // to reuse space and will properly clean up unused stuff
        if (&src != this)
            copyFrom(src, Stage::Instance);
        return *this;
    }

    // Back up to the stage just before g if this subsystem thinks
    // it is already at g or beyond. Note that we may be backing up
    // over many stages here. Careful: invalidating the stage
    // for a subsystem must also invalidate the same stage for all
    // the other subsystems and the system as a whole but we don't
    // take care of that here. Also, you can't invalidate Stage::Empty.
    void invalidateStageJustThisSubsystem(Stage g) {
        assert(g > Stage::Empty);
        restoreToStage(g.prev());
    }

    // Advance from stage g-1 to stage g. This is called at the end
    // of realize(g). You can't use this to "advance" to Stage::Empty.
    // It is a fatal error if the current stage isn't g-1.
    void advanceToStage(Stage g) const {
        assert(g > Stage::Empty);
        assert(currentStage == g.prev());

        // This validates whatever the current version number is of Stage g.
        currentStage = g;
    }


    void clearReferencesToModelStageGlobals() {
        qstart.invalidate(); ustart.invalidate(); 
        zstart.invalidate();
        q.clear(); u.clear(); z.clear();
        uWeights.clear(); zWeights.clear();
        qdot.clear(); udot.clear(); zdot.clear(); qdotdot.clear();
    }

    void clearReferencesToInstanceStageGlobals() {
        // These are late-allocated state variables.
        qerrWeights.clear(); uerrWeights.clear();

        // These are all mutable cache entries.
        qerrstart.invalidate();uerrstart.invalidate();udoterrstart.invalidate();
        qerr.clear();uerr.clear();udoterr.clear();multipliers.clear();

        for (int j=0; j<Stage::NValid; ++j) {
            triggerstart[j].invalidate();
            triggers[j].clear();
        }
    }

    QIndex getNextQIndex() const {
        if (qInfo.empty()) return QIndex(0);
        const ContinuousVarInfo& last = qInfo.back();
        return QIndex(last.getFirstIndex()+last.getNumVars());
    }
    UIndex getNextUIndex() const {
        if (uInfo.empty()) return UIndex(0);
        const ContinuousVarInfo& last = uInfo.back();
        return UIndex(last.getFirstIndex()+last.getNumVars());
    }
    ZIndex getNextZIndex() const {
        if (zInfo.empty()) return ZIndex(0);
        const ContinuousVarInfo& last = zInfo.back();
        return ZIndex(last.getFirstIndex()+last.getNumVars());
    }

    QErrIndex getNextQErrIndex() const {
        if (qerrInfo.empty()) return QErrIndex(0);
        const ConstraintErrInfo& last = qerrInfo.back();
        return QErrIndex(last.getFirstIndex()+last.getNumErrs());
    }
    UErrIndex getNextUErrIndex() const {
        if (uerrInfo.empty()) return UErrIndex(0);
        const ConstraintErrInfo& last = uerrInfo.back();
        return UErrIndex(last.getFirstIndex()+last.getNumErrs());
    }
    UDotErrIndex getNextUDotErrIndex() const {
        if (udoterrInfo.empty()) return UDotErrIndex(0);
        const ConstraintErrInfo& last = udoterrInfo.back();
        return UDotErrIndex(last.getFirstIndex()+last.getNumErrs());
    }
    DiscreteVariableIndex getNextDiscreteVariableIndex() const {
        return DiscreteVariableIndex(discreteInfo.size());
    }
    CacheEntryIndex getNextCacheEntryIndex() const {
        return CacheEntryIndex(cacheInfo.size());
    }
    EventTriggerByStageIndex getNextEventTriggerByStageIndex(Stage g) const {
        if (triggerInfo[g].empty()) return EventTriggerByStageIndex(0);
        const TriggerInfo& last = triggerInfo[g].back();
        return EventTriggerByStageIndex
            (last.getFirstIndex()+last.getNumSlots());
    }


    String name;
    String version;

        // DEFINITIONS //

    // State variables (continuous or discrete) can be defined (allocated) 
    // during realization of Topology or Model stages. Cache entries,
    // constraint error slots, and event trigger slots can be defined during 
    // realization of Topology, Model, or Instance stages. No further 
    // allocations are allowed. Then, when one of these stages is invalidated, 
    // all the definitions that occurred during realization of that stage must 
    // be forgotten.
    // 
    // To do that the allocation entries are stored in arrays which are really 
    // stacks, with definitions pushed onto the ends as the stage is advanced 
    // and popped off the ends as the stage is reduced.

    // Topology and Model stage definitions. 
    Array_<ContinuousVarInfo>      qInfo, uInfo, zInfo;
    Array_<DiscreteVarInfo>        discreteInfo;

    // Topology, Model, and Instance stage definitions.
    mutable Array_<ConstraintErrInfo>   qerrInfo, uerrInfo, udoterrInfo;
    mutable Array_<TriggerInfo>         triggerInfo[Stage::NValid];
    mutable Array_<CacheEntryInfo>      cacheInfo;
   
        // GLOBAL RESOURCE ALLOCATIONS //

    // These are our own private views into partitions of the global
    // state and cache entries of the same names. The State will assign
    // contiguous blocks to this subsystem when the *System* stage is raised
    // to Model or Instance stage, and they are invalidated whenever that 
    // stage is invalidated. The starting indices are filled in here at 
    // the time the views are built.

    // Model stage global resources and views into them.
    SystemQIndex    qstart;
    SystemUIndex    ustart;
    SystemZIndex    zstart;
    Vector          q, u, z;
    Vector          uWeights, zWeights;

    mutable Vector  qdot, udot, zdot, qdotdot;

    // Instance stage global resources and views into them.
    Vector          qerrWeights, uerrWeights;

    // Note that multipliers just use the same indices as udoterr.
    mutable SystemQErrIndex                 qerrstart;
    mutable SystemUErrIndex                 uerrstart;
    mutable SystemUDotErrIndex              udoterrstart;
    mutable SystemEventTriggerByStageIndex  triggerstart[Stage::NValid];

    mutable Vector  qerr, uerr;

    mutable Vector  udoterr, multipliers; // same size and partioning
    mutable Vector  triggers[Stage::NValid];

    // The currentStage is the highest stage of this subsystem that is valid,
    // meaning that it has been realized since the last change to any variable
    // that might affect it. All stages less than currentStage are also valid,
    // and all higher stages are invalid.
    // Each stage has a "stage version" which is like a serial number that is
    // bumped every time the stage is invalidated by a variable change. Cache
    // entries that are calculated from a particular stage version can record
    // the version number to allow a quick check later -- if the current version
    // of a cache entry's dependsOn stage is different than the one stored with
    // the cache entry, then that cache entry cannot be valid.
    mutable Stage        currentStage;
    mutable StageVersion stageVersions[Stage::NValid];

private:
    // This is for use in constructors and for resetting an existing State into
    // its just-constructed condition.
    void initialize() {
        clearAllStacks();
        qstart.invalidate();ustart.invalidate();zstart.invalidate();
        qerrstart.invalidate();uerrstart.invalidate();udoterrstart.invalidate();
        for (int j=0; j<Stage::NValid; ++j) {
            triggerstart[j].invalidate();
            stageVersions[j] = 1; // never 0
        }
        currentStage = Stage::Empty;
    }

    // Manage allocation stacks.

    void clearContinuousVars();
    void clearConstraintErrs();
    void clearDiscreteVars();
    void clearEventTriggers(int g);
    void clearCache();

    void clearAllStacks();

    void popContinuousVarsBackToStage(const Stage& g);
    void popDiscreteVarsBackToStage(const Stage& g);
    void popConstraintErrsBackToStage(const Stage& g);
    void popCacheBackToStage(const Stage& g);
    void popEventTriggersBackToStage(const Stage& g);

    void popAllStacksBackToStage(const Stage& g);

    // Call once each for qInfo, uInfo, zInfo.
    void copyContinuousVarInfoThroughStage
       (const Array_<ContinuousVarInfo>& src, const Stage& g,
        Array_<ContinuousVarInfo>& dest);

    void copyDiscreteVarsThroughStage
       (const Array_<DiscreteVarInfo>& src, const Stage& g);

    // Call once each for qerrInfo, uerrInfo, udoterrInfo.
    void copyConstraintErrInfoThroughStage
       (const Array_<ConstraintErrInfo>& src, const Stage& g,
        Array_<ConstraintErrInfo>& dest);

    void copyCacheThroughStage
       (const Array_<CacheEntryInfo>& src, const Stage& g);

    void copyEventsThroughStage
       (const Array_<TriggerInfo>& src, const Stage& g,
        Array_<TriggerInfo>& dest);

    void copyAllStacksThroughStage(const PerSubsystemInfo& src, const Stage& g);

    // Restore this subsystem to the way it last was at realize(g) for a given 
    // Stage g; that is, invalidate all stages > g. Allocations will be 
    // forgotten as Instance, Model, and Topology stages are invalidated.
    void restoreToStage(Stage g);

    // Utility which makes "this" a copy of the source subsystem exactly as it
    // was after being realized to stage maxStage. If maxStage >= Model then
    // all the subsystem-private state variables will be copied, but only
    // cached computations up through maxStage come through. We clear
    // our references to global variables regardless -- those will have to
    // be repaired at the System (State global) level.
    void copyFrom(const PerSubsystemInfo& src, Stage maxStage);
};


//==============================================================================
//                                 STATE IMPL
//==============================================================================

class StateImpl {
public:
    StateImpl() 
    :   t(NaN), currentSystemStage(Stage::Empty) {initializeStageVersions();} 

    // We'll do the copy constructor and assignment explicitly here
    // to get tight control over what's allowed.
    StateImpl(const StateImpl& src);

    StateImpl& operator=(const StateImpl& src);

    ~StateImpl() {}   // default destructor

    // Copies all the variables but not the cache.
    StateImpl* clone() const {return new StateImpl(*this);}

    const Stage& getSystemStage() const {return currentSystemStage;}
    Stage&       updSystemStage() const {return currentSystemStage;} // mutable


    const PerSubsystemInfo& getSubsystem(int subsystem) const {
        SimTK_INDEXCHECK(subsystem, (int)subsystems.size(), 
                         "StateImpl::getSubsystem()");
        return subsystems[subsystem];
    }

    PerSubsystemInfo& updSubsystem(int subsystem) {
        SimTK_INDEXCHECK(subsystem, (int)subsystems.size(), 
                         "StateImpl::updSubsystem()");
        return subsystems[subsystem];
    }

    const Stage& getSubsystemStage(int subsystem) const {
        return subsystems[subsystem].currentStage;
    }
    Stage& updSubsystemStage(int subsystem) const {
        return subsystems[subsystem].currentStage; // mutable
    }

    const StageVersion* getSubsystemStageVersions(int subsystem) const {
        return subsystems[subsystem].stageVersions;
    }

    // Back up the System stage just before stg if it thinks
    // it is already at stg or beyond. Note that we may be backing up
    // over many stages here. Careful: invalidating the stage
    // for the system must also invalidate the same stage for all
    // the subsystems (because we trash the shared resource pool
    // here if we back up earlier than Stage::Model) but we don't
    // take care of that here. Also, you can't invalidate Stage::Empty.
    void invalidateJustSystemStage(Stage stg);

    // Advance the System stage from stg-1 to stg. It is a fatal error if
    // we're not already at stg-1, and you can't advance to Stage::Empty.
    // Also, you can't advance the system to stg unless ALL subsystems have
    // already gotten there.
    void advanceSystemToStage(Stage stg) const;
    
    void setNumSubsystems(int i) {
        assert(i >= 0);
        subsystems.clear();
        subsystems.resize(i);
    }
    
    void initializeSubsystem
       (SubsystemIndex i, const String& name, const String& version) {
        updSubsystem(i).name = name;
        updSubsystem(i).version = version;
    }
      
    SubsystemIndex addSubsystem(const String& name, const String& version) {
        const SubsystemIndex sx(subsystems.size());
        subsystems.push_back(PerSubsystemInfo(name,version));
        return sx;
    }
    
    int getNumSubsystems() const {return (int)subsystems.size();}
    
    const String& getSubsystemName(SubsystemIndex subsys) const {
        return subsystems[subsys].name;
    }
    const String& getSubsystemVersion(SubsystemIndex subsys) const {
        return subsystems[subsys].version;
    }

    // Make sure the stage is no higher than g-1 for *any* subsystem and
    // hence for the system stage also. TODO: this should be more selective.
    void invalidateAll(Stage g) {
        invalidateJustSystemStage(g);
        for (SubsystemIndex i(0); i<(int)subsystems.size(); ++i)
            subsystems[i].invalidateStageJustThisSubsystem(g);
    }

    // Make sure the stage is no higher than g-1 for *any* subsystem and
    // hence for the system stage also. Same as invalidateAll() except this
    // requires only const access and can't be used for g below Instance.
    void invalidateAllCacheAtOrAbove(Stage g) const {
        SimTK_STAGECHECK_GE_ALWAYS(g, Stage::Instance, 
            "StateImpl::invalidateAllCacheAtOrAbove()");

        // We promise not to hurt this State; get non-const access just so
        // we can call these methods.
        StateImpl* mthis = const_cast<StateImpl*>(this);
        mthis->invalidateJustSystemStage(g);
        for (SubsystemIndex i(0); i<(int)subsystems.size(); ++i)
            mthis->subsystems[i].invalidateStageJustThisSubsystem(g);
    }
    
    // Move the stage for a particular subsystem from g-1 to g. No other
    // subsystems are affected, nor the global system stage.
    void advanceSubsystemToStage(SubsystemIndex subsys, Stage g) const {
        subsystems[subsys].advanceToStage(g);
        // We don't automatically advance the System stage even if this brings
        // ALL the subsystems up to stage g.
    }
     
    // We don't expect State entry allocations to be performance critical so
    // we'll keep error checking on even in Release mode.
    
    QIndex allocateQ(SubsystemIndex subsys, const Vector& qInit) {
        SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), Stage::Model, 
                                   "StateImpl::allocateQ()");
        // We are currently realizing the next stage.
        const Stage allocStage = getSubsystemStage(subsys).next();
        PerSubsystemInfo& ss = subsystems[subsys];
        const QIndex nxt(ss.getNextQIndex());
        ss.qInfo.push_back(ContinuousVarInfo(allocStage,nxt,qInit,Vector())); 
        return nxt;
    }
    
    UIndex allocateU(SubsystemIndex subsys, const Vector& uInit) {
        SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), Stage::Model, 
                                   "StateImpl::allocateU()");
        const Stage allocStage = getSubsystemStage(subsys).next();
        PerSubsystemInfo& ss = subsystems[subsys];
        const UIndex nxt(ss.getNextUIndex());
        ss.uInfo.push_back(ContinuousVarInfo(allocStage,nxt,uInit,Vector())); 
        return nxt;
    }
    ZIndex allocateZ(SubsystemIndex subsys, const Vector& zInit) {
        SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), Stage::Model, 
                                   "StateImpl::allocateZ()");
        const Stage allocStage = getSubsystemStage(subsys).next();
        PerSubsystemInfo& ss = subsystems[subsys];
        const ZIndex nxt(ss.getNextZIndex());
        ss.zInfo.push_back(ContinuousVarInfo(allocStage,nxt,zInit,Vector())); 
        return nxt;
    }
    
    QErrIndex allocateQErr(SubsystemIndex subsys, int nqerr) const {
        SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), Stage::Instance, 
                                   "StateImpl::allocateQErr()");
        const Stage allocStage = getSubsystemStage(subsys).next();
        const PerSubsystemInfo& ss = subsystems[subsys];
        const QErrIndex nxt(ss.getNextQErrIndex());
        ss.qerrInfo.push_back(ConstraintErrInfo(allocStage,nxt,nqerr,Vector())); 
        return nxt;
    }
    UErrIndex allocateUErr(SubsystemIndex subsys, int nuerr) const {
        SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), 
                                   Stage::Instance,"StateImpl::allocateUErr()");
        const Stage allocStage = getSubsystemStage(subsys).next();
        const PerSubsystemInfo& ss = subsystems[subsys];
        const UErrIndex nxt(ss.getNextUErrIndex());
        ss.uerrInfo.push_back(ConstraintErrInfo(allocStage,nxt,nuerr,Vector())); 
        return nxt;
    }
    UDotErrIndex allocateUDotErr(SubsystemIndex subsys, int nudoterr) const {
        SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), Stage::Instance,
                                   "StateImpl::allocateUDotErr()");
        const Stage allocStage = getSubsystemStage(subsys).next();
        const PerSubsystemInfo& ss = subsystems[subsys];
        const UDotErrIndex nxt(ss.getNextUDotErrIndex());
        ss.udoterrInfo.push_back
           (ConstraintErrInfo(allocStage,nxt,nudoterr,Vector())); 
        return nxt;
    }
    EventTriggerByStageIndex allocateEventTrigger
       (SubsystemIndex subsys, Stage g, int nt) const {
        SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), Stage::Instance, 
                                   "StateImpl::allocateEventTrigger()");
        const Stage allocStage = getSubsystemStage(subsys).next();
        const PerSubsystemInfo& ss = subsystems[subsys];
        const EventTriggerByStageIndex 
            nxt(ss.getNextEventTriggerByStageIndex(g));
        ss.triggerInfo[g].push_back(TriggerInfo(allocStage,nxt,nt)); 
        return nxt;
    }
    
    // Topology- and Model-stage State variables can only be added during 
    // construction; that is, while stage <= Topology. Other entries can be 
    // added while stage < Model.
    DiscreteVariableIndex allocateDiscreteVariable
       (SubsystemIndex subsys, Stage invalidates, AbstractValue* vp) 
    {
        SimTK_STAGECHECK_RANGE_ALWAYS(Stage(Stage::LowestRuntime).prev(),
                                      invalidates, Stage::HighestRuntime, 
            "StateImpl::allocateDiscreteVariable()");
    
        const Stage maxAcceptable = (invalidates <= Stage::Model 
                                     ? Stage::Empty : Stage::Topology);
        SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), 
            maxAcceptable.next(), "StateImpl::allocateDiscreteVariable()");

        const Stage allocStage = getSubsystemStage(subsys).next();    
        PerSubsystemInfo& ss = subsystems[subsys];
        const DiscreteVariableIndex nxt(ss.getNextDiscreteVariableIndex());
        ss.discreteInfo.push_back
           (DiscreteVarInfo(allocStage,invalidates,vp));
        return nxt;
    }
    
    // Cache entries can be allocated while stage < Instance.
    CacheEntryIndex allocateCacheEntry
       (SubsystemIndex subsys, Stage dependsOn, Stage computedBy,
        AbstractValue* vp) const 
    {
        SimTK_STAGECHECK_RANGE_ALWAYS(Stage(Stage::LowestRuntime).prev(), 
                                      dependsOn, Stage::HighestRuntime, 
            "StateImpl::allocateCacheEntry()");
        SimTK_STAGECHECK_RANGE_ALWAYS(dependsOn, computedBy, Stage::Infinity, 
            "StateImpl::allocateCacheEntry()");
        SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), 
            Stage::Instance, "StateImpl::allocateCacheEntry()");

        const Stage allocStage = getSubsystemStage(subsys).next();    
        const PerSubsystemInfo& ss = subsystems[subsys];
        const CacheEntryIndex nxt(ss.getNextCacheEntryIndex());
        ss.cacheInfo.push_back
           (CacheEntryInfo(allocStage, dependsOn, computedBy, vp)); //mutable
        return nxt;
    }

    // Allocate a discrete variable and a corresponding cache entry for
    // updating it, and connect them together.
    DiscreteVariableIndex allocateAutoUpdateDiscreteVariable
       (SubsystemIndex subsys, Stage invalidates, AbstractValue* vp,
        Stage updateDependsOn)
    {
        const DiscreteVariableIndex dx = 
            allocateDiscreteVariable(subsys,invalidates,vp->clone());
        const CacheEntryIndex       cx = 
            allocateCacheEntry(subsys,updateDependsOn,Stage::Infinity,vp);

        PerSubsystemInfo& ss = subsystems[subsys];
        DiscreteVarInfo& dvinfo = ss.discreteInfo[dx];
        CacheEntryInfo&  ceinfo = ss.cacheInfo[cx];
        dvinfo.setAutoUpdateEntry(cx);
        ceinfo.setAssociatedVar(dx);
        return dx;
    }
    
        // State dimensions for shared continuous variables.
    
    int getNY() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::getNY()");
        return y.size();
    }
    
    SystemYIndex getQStart() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::getQStart()");
        return SystemYIndex(0); // q's come first
    }
    int getNQ() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::getNQ()");
        return q.size();
    }
    
    SystemYIndex getUStart() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::getUStart()");
        return SystemYIndex(q.size()); // u's come right after q's
    }
    int getNU() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::getNU()");
        return u.size();
    }
    
    SystemYIndex getZStart() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::getZStart()");
        return SystemYIndex(q.size() + u.size()); // q,u, then z
    }
    int getNZ() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::getNZ()");
        return z.size();
    }
    
    int getNYErr() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
                            "StateImpl::getNYErr()");
        return yerr.size();
    }
    
    SystemYErrIndex getQErrStart() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
                            "StateImpl::getQErrStart()");
        return SystemYErrIndex(0); // qerr's come first
    }
    int getNQErr() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
                            "StateImpl::getNQErr()");
        return qerr.size();
    }
    
    SystemYErrIndex getUErrStart() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
                            "StateImpl::getUErrStart()");
        return SystemYErrIndex(qerr.size()); // uerr's follow qerrs
    }
    int getNUErr() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
                            "StateImpl::getNUErr()");
        return uerr.size();
    }
    
    // UDot errors are independent of qerr & uerr.
    // This is used for multipliers also.
    int getNUDotErr() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
                            "StateImpl::getNUDotErr()");
        return udoterr.size();
    }
    
    int getNEventTriggers() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
                            "StateImpl::getNEventTriggers()");
        return allTriggers.size();
    }
    
    SystemEventTriggerIndex getEventTriggerStartByStage(Stage g) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
                            "StateImpl::getEventTriggerStartByStage()");
        int nxt = 0;
        for (int j=0; j<g; ++j)
            nxt += triggers[j].size();
        return SystemEventTriggerIndex(nxt); // g starts where g-1 leaves off
    }
    
    int getNEventTriggersByStage(Stage g) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
                            "StateImpl::getNEventTriggersByStage()");
        return triggers[g].size();
    }
    
        // Subsystem dimensions.
    
    SystemQIndex getQStart(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::getQStart(subsys)");
        return getSubsystem(subsys).qstart;
    }
    int getNQ(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::getNQ(subsys)");
        return getSubsystem(subsys).q.size();
    }
    
    SystemUIndex getUStart(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::getUStart(subsys)");
        return getSubsystem(subsys).ustart;
    }
    int getNU(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::getNU(subsys)");
        return getSubsystem(subsys).u.size();
    }
    
    SystemZIndex getZStart(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::getZStart(subsys)");
        return getSubsystem(subsys).zstart;
    }
    int getNZ(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::getNZ(subsys)");
        return getSubsystem(subsys).z.size();
    }
    
    SystemQErrIndex getQErrStart(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
                            "StateImpl::getQErrStart(subsys)");
        return getSubsystem(subsys).qerrstart;
    }
    int getNQErr(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
                            "StateImpl::getNQErr(subsys)");
        return getSubsystem(subsys).qerr.size();
    }
    
    SystemUErrIndex getUErrStart(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
                            "StateImpl::getUErrStart(subsys)");
        return getSubsystem(subsys).uerrstart;
    }
    int getNUErr(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
                            "StateImpl::getNUErr(subsys)");
        return getSubsystem(subsys).uerr.size();
    }
    
    // These are used for multipliers also.
    SystemUDotErrIndex getUDotErrStart(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
                            "StateImpl::getUDotErrStart(subsys)");
        return getSubsystem(subsys).udoterrstart;
    }
    int getNUDotErr(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
                            "StateImpl::getNUDotErr(subsys)");
        return getSubsystem(subsys).udoterr.size();
    }
    
    SystemEventTriggerByStageIndex getEventTriggerStartByStage
       (SubsystemIndex subsys, Stage g) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
                            "StateImpl::getEventTriggerStartByStage(subsys)");
        return getSubsystem(subsys).triggerstart[g];
    }
    
    int getNEventTriggersByStage(SubsystemIndex subsys, Stage g) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
                            "StateImpl::getNEventTriggersByStage(subsys)");
        return getSubsystem(subsys).triggers[g].size();
    }
    
        // Per-subsystem access to the global shared variables.
    
    const Vector& getQ(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::getQ(subsys)");
        return getSubsystem(subsys).q;
    }
    const Vector& getU(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::getU(subsys)");
        return getSubsystem(subsys).u;
    }
    const Vector& getZ(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::getZ(subsys)");
        return getSubsystem(subsys).z;
    }

    const Vector& getUWeights(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
            "StateImpl::getUWeights(subsys)");
        return getSubsystem(subsys).uWeights;
    }
    const Vector& getZWeights(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
            "StateImpl::getZWeights(subsys)");
        return getSubsystem(subsys).zWeights;
    }

    const Vector& getQDot(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::getQDot(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Velocity, 
                            "StateImpl::getQDot(subsys)");
        return getSubsystem(subsys).qdot;
    }
    const Vector& getUDot(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::getUDot(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Acceleration, 
                            "StateImpl::getUDot(subsys)");
        return getSubsystem(subsys).udot;
    }
    const Vector& getZDot(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::getZDot(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Dynamics, 
                            "StateImpl::getZDot(subsys)");
        return getSubsystem(subsys).zdot;
    }
    const Vector& getQDotDot(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::getQDotDot(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Acceleration, 
                            "StateImpl::getQDotDot(subsys)");
        return getSubsystem(subsys).qdotdot;
    }
    
    Vector& updQ(SubsystemIndex subsys) {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::updQ(subsys)");
        invalidateAll(Stage::Position);
        return updSubsystem(subsys).q;
    }
    Vector& updU(SubsystemIndex subsys) {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::updU(subsys)");
        invalidateAll(Stage::Velocity);
        return updSubsystem(subsys).u;
    }
    Vector& updZ(SubsystemIndex subsys) {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::updZ(subsys)");
        invalidateAll(Stage::Dynamics);
        return updSubsystem(subsys).z;
    }

    Vector& updUWeights(SubsystemIndex subsys) {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
            "StateImpl::updUWeights(subsys)");
        invalidateAll(Stage::Report);
        return updSubsystem(subsys).uWeights;
    }
    Vector& updZWeights(SubsystemIndex subsys) {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
            "StateImpl::updZWeights(subsys)");
        invalidateAll(Stage::Report);
        return updSubsystem(subsys).zWeights;
    }
    
        // These are mutable so the routines are const.
    
    Vector& updQDot(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::updQDot(subsys)");
        return getSubsystem(subsys).qdot;
    }
    Vector& updUDot(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::updUDot(subsys)");
        return getSubsystem(subsys).udot;
    }
    Vector& updZDot(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::updZDot(subsys)");
        return getSubsystem(subsys).zdot;
    }
    Vector& updQDotDot(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::updQDotDot(subsys)");
        return getSubsystem(subsys).qdotdot;
    }
    
    
    const Vector& getQErr(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
                            "StateImpl::getQErr(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Position, 
                            "StateImpl::getQErr(subsys)");
        return getSubsystem(subsys).qerr;
    }
    const Vector& getUErr(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
                            "StateImpl::getUErr(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Velocity, 
                            "StateImpl::getUErr(subsys)");
        return getSubsystem(subsys).uerr;
    }

    const Vector& getQErrWeights(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
            "StateImpl::getQErrWeights(subsys)");
        return getSubsystem(subsys).qerrWeights;
    }
    const Vector& getUErrWeights(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
            "StateImpl::getUErrWeights(subsys)");
        return getSubsystem(subsys).uerrWeights;
    }

    const Vector& getUDotErr(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
                            "StateImpl::getUDotErr(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Acceleration, 
                            "StateImpl::getUDotErr(subsys)");
        return getSubsystem(subsys).udoterr;
    }
    const Vector& getMultipliers(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
                            "StateImpl::getMultipliers(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Acceleration, 
                            "StateImpl::getMultipliers(subsys)");
        return getSubsystem(subsys).multipliers;
    }
    
    const Vector& getEventTriggersByStage(SubsystemIndex subsys, Stage g) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
                            "StateImpl::getEventTriggersByStage(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), g, 
                            "StateImpl::getEventTriggersByStage(subsys)");
        return getSubsystem(subsys).triggers[g];
    }
    
    Vector& updQErr(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
                            "StateImpl::updQErr(subsys)");
        return getSubsystem(subsys).qerr;
    }
    Vector& updUErr(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
                            "StateImpl::updUErr(subsys)");
        return getSubsystem(subsys).uerr;
    }
    
    Vector& updQErrWeights(SubsystemIndex subsys) {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
            "StateImpl::updQErrWeights(subsys)");
        invalidateAll(Stage::Position);
        return updSubsystem(subsys).qerrWeights;
    }
    Vector& updUErrWeights(SubsystemIndex subsys) {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
            "StateImpl::updUErrWeights(subsys)");
        invalidateAll(Stage::Velocity);
        return updSubsystem(subsys).uerrWeights;
    }

    Vector& updUDotErr(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
                            "StateImpl::updUDotErr(subsys)");
        return getSubsystem(subsys).udoterr;
    }
    Vector& updMultipliers(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
                            "StateImpl::updMultipliers(subsys)");
        return getSubsystem(subsys).multipliers;
    }
    Vector& updEventTriggersByStage(SubsystemIndex subsys, Stage g) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
                            "StateImpl::updEventTriggersByStage(subsys)");
        return getSubsystem(subsys).triggers[g];
    }
    
        // Direct access to the global shared state and cache entries.
        // Time is allocated in Stage::Topology, State in Stage::Model, and
        // Cache in Stage::Instance.
    
    const Real& getTime() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Topology, 
                            "StateImpl::getTime()");
        return t;
    }
    
    const Vector& getY() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::getY()");
        return y;
    }
    
    const Vector& getQ() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::getQ()");
        return q;
    }
    
    const Vector& getU() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::getU()");
        return u;
    }
    
    const Vector& getZ() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::getZ()");
        return z;
    }
        
    const Vector& getUWeights() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::getUWeights()");
        return uWeights;
    }
    
    const Vector& getZWeights() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::getZWeights()");
        return zWeights;
    }
       
    // You can call these as long as stage >= allocation stage, but the
    // stage will be backed up if necessary to one stage prior to the invalidated stage.
    Real& updTime() {  // Back to Stage::Time-1
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Topology, 
                            "StateImpl::updTime()");
        invalidateAll(Stage::Time);
        return t;
    }
    
    Vector& updY() {    // Back to Stage::Position-1
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::updY()");
        invalidateAll(Stage::Position);
        return y;
    }
    
    Vector& updQ() {    // Stage::Position-1
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::updQ()");
        invalidateAll(Stage::Position);
        return q;
    }
    
    Vector& updU() {     // Stage::Velocity-1
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::updU()");
        invalidateAll(Stage::Velocity);
        return u;
    }
    
    Vector& updZ() {     // Stage::Dynamics-1
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::updZ()");
        invalidateAll(Stage::Dynamics);
        return z;
    }
     
    Vector& updUWeights() { // Invalidates Report stage only
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
            "StateImpl::updUWeights()");
        invalidateAll(Stage::Report);
        return uWeights;
    }
    
    Vector& updZWeights() { // Invalidates Report stage only
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
            "StateImpl::updZWeights()");
        invalidateAll(Stage::Dynamics);
        return zWeights;
    }   

    // These cache entries you can get at their "computedBy" stages.
    const Vector& getYDot() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Acceleration, 
                            "StateImpl::getYDot()");
        return ydot;
    }
    
    const Vector& getQDot() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Velocity, 
                            "StateImpl::getQDot()");
        return qdot;
    }
    
    const Vector& getZDot() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Dynamics, 
                            "StateImpl::getZDot()");
        return zdot;
    }
    
    const Vector& getUDot() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Acceleration, 
                            "StateImpl::getUDot()");
        return udot;
    }
    
    const Vector& getQDotDot() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Acceleration, 
                            "StateImpl::getQDotDot()");
        return qdotdot;
    }
    
    // Cache updates are allowed any time after they have been allocated.
    Vector& updYDot() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::updYDot()");
        return ydot;
    }
    
    Vector& updQDot() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::updQDot()");
        return qdot;
    }
    
    Vector& updUDot() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::updUDot()");
        return udot;
    }
    
    Vector& updZDot() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::updZDot()");
        return zdot;
    }
    
    Vector& updQDotDot() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, 
                            "StateImpl::updQDotDot()");
        return qdotdot;
    }
    
    
    const Vector& getYErr() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Velocity, 
                            "StateImpl::getYErr()");
        return yerr;
    }
    
    const Vector& getQErr() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Position, 
                            "StateImpl::getQErr()");
        return qerr;
    }
    const Vector& getUErr() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Velocity, 
                            "StateImpl::getUErr()");
        return uerr;
    }
    
    const Vector& getQErrWeights() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
            "StateImpl::getQErrWeights()");
        return qerrWeights;
    }
    const Vector& getUErrWeights() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
            "StateImpl::getUErrWeights()");
        return uerrWeights;
    }

    const Vector& getUDotErr() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Acceleration, 
                            "StateImpl::getUDotErr()");
        return udoterr;
    }
    const Vector& getMultipliers() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Acceleration, 
                            "StateImpl::getMultipliers()");
        return multipliers;
    }
    
    Vector& updYErr() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
                            "StateImpl::updYErr()");
        return yerr;
    }
    Vector& updQErr() const{
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
                            "StateImpl::updQErr()");
        return qerr;
    }
    Vector& updUErr() const{
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
                            "StateImpl::updUErr()");
        return uerr;
    }

    Vector& updQErrWeights() {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
            "StateImpl::updQErrWeights()");
        invalidateAll(Stage::Position);
        return qerrWeights;
    }
    Vector& updUErrWeights() {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
            "StateImpl::updUErrWeights()");
        invalidateAll(Stage::Velocity);
        return uerrWeights;
    }

    Vector& updUDotErr() const{
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
                            "StateImpl::updUDotErr()");
        return udoterr;
    }
    Vector& updMultipliers() const{
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
                            "StateImpl::updMultipliers()");
        return multipliers;
    }
    
    const Vector& getEventTriggers() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Acceleration, 
                            "StateImpl::getEventTriggers()");
        return allTriggers;
    }
    const Vector& getEventTriggersByStage(Stage g) const {
        SimTK_STAGECHECK_GE(getSystemStage(), g, 
                            "StateImpl::getEventTriggersByStage()");
        return triggers[g];
    }
    
    // These are mutable; hence 'const'.
    Vector& updEventTriggers() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
                            "StateImpl::updEventTriggers()");
        return allTriggers;
    }
    Vector& updEventTriggersByStage(Stage g) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, 
                            "StateImpl::updEventTriggersByStage()");
        return triggers[g];
    }

    CacheEntryIndex getDiscreteVarUpdateIndex
       (SubsystemIndex subsys, DiscreteVariableIndex index) const {
        const PerSubsystemInfo& ss = subsystems[subsys];
        SimTK_INDEXCHECK(index,(int)ss.discreteInfo.size(),
            "StateImpl::getDiscreteVarUpdateIndex()");
        const DiscreteVarInfo& dv = ss.discreteInfo[index];
        return dv.getAutoUpdateEntry();
    } 

    Stage getDiscreteVarAllocationStage
       (SubsystemIndex subsys, DiscreteVariableIndex index) const {
        const PerSubsystemInfo& ss = subsystems[subsys];
        SimTK_INDEXCHECK(index,(int)ss.discreteInfo.size(),
            "StateImpl::getDiscreteVarAllocationStage()");
        const DiscreteVarInfo& dv = ss.discreteInfo[index];
        return dv.getAllocationStage();
    } 

    Stage getDiscreteVarInvalidatesStage
       (SubsystemIndex subsys, DiscreteVariableIndex index) const {
        const PerSubsystemInfo& ss = subsystems[subsys];
        SimTK_INDEXCHECK(index,(int)ss.discreteInfo.size(),
            "StateImpl::getDiscreteVarInvalidatesStage()");
        const DiscreteVarInfo& dv = ss.discreteInfo[index];
        return dv.getInvalidatedStage();
    } 

    // You can access at any time a variable that was allocated during 
    // realizeTopology(), but don't access others until you have done 
    // realizeModel().
    const AbstractValue& 
    getDiscreteVariable(SubsystemIndex subsys, DiscreteVariableIndex index) const {
        const PerSubsystemInfo& ss = subsystems[subsys];
    
        SimTK_INDEXCHECK(index,(int)ss.discreteInfo.size(),
                         "StateImpl::getDiscreteVariable()");
        const DiscreteVarInfo& dv = ss.discreteInfo[index];
    
        if (dv.getAllocationStage() > Stage::Topology) {
            SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Model, 
                                "StateImpl::getDiscreteVariable()");
        }
    
        return dv.getValue();
    }
    Real getDiscreteVarLastUpdateTime(SubsystemIndex subsys, DiscreteVariableIndex index) const {
        const PerSubsystemInfo& ss = subsystems[subsys];
        SimTK_INDEXCHECK(index,(int)ss.discreteInfo.size(),
                         "StateImpl::getDiscreteVarLastUpdateTime()");
        const DiscreteVarInfo& dv = ss.discreteInfo[index];
        return dv.getTimeLastUpdated();
    }

    const AbstractValue& getDiscreteVarUpdateValue
       (SubsystemIndex subsys, DiscreteVariableIndex index) const {
        const CacheEntryIndex cx = getDiscreteVarUpdateIndex(subsys,index);
        SimTK_ERRCHK2(cx.isValid(), "StateImpl::getDiscreteVarUpdateValue()", 
            "Subsystem %d has a discrete variable %d but it does not have an"
            " associated update cache variable.", (int)subsys, (int)index);
        return getCacheEntry(subsys, cx);
    }
    AbstractValue& updDiscreteVarUpdateValue
       (SubsystemIndex subsys, DiscreteVariableIndex index) const {
        const CacheEntryIndex cx = getDiscreteVarUpdateIndex(subsys,index);
        SimTK_ERRCHK2(cx.isValid(), "StateImpl::updDiscreteVarUpdateValue()", 
            "Subsystem %d has a discrete variable %d but it does not have an"
            " associated update cache variable.", (int)subsys, (int)index);
        return updCacheEntry(subsys, cx);
    }
    bool isDiscreteVarUpdateValueRealized
       (SubsystemIndex subsys, DiscreteVariableIndex index) const {
        const CacheEntryIndex cx = getDiscreteVarUpdateIndex(subsys,index);
        SimTK_ERRCHK2(cx.isValid(), "StateImpl::isDiscreteVarUpdateValueRealized()", 
            "Subsystem %d has a discrete variable %d but it does not have an"
            " associated update cache variable.", (int)subsys, (int)index);
        return isCacheValueRealized(subsys, cx);
    }
    void markDiscreteVarUpdateValueRealized
       (SubsystemIndex subsys, DiscreteVariableIndex index) const {
        const CacheEntryIndex cx = getDiscreteVarUpdateIndex(subsys,index);
        SimTK_ERRCHK2(cx.isValid(), 
            "StateImpl::markDiscreteVarUpdateValueRealized()", 
            "Subsystem %d has a discrete variable %d but it does not have an"
            " associated update cache variable.", (int)subsys, (int)index);
        markCacheValueRealized(subsys, cx);
    }

    // You can update at any time a variable that was allocated during 
    // realizeTopology(), but later variables must wait until you have done 
    // realizeModel(). This always backs the stage up to one earlier than the 
    // variable's "invalidates" stage, and if there is an auto-update cache 
    // entry it is also invalidated, regardless of its "depends on" stage.
    AbstractValue& 
    updDiscreteVariable(SubsystemIndex subsys, DiscreteVariableIndex index) {
        PerSubsystemInfo& ss = subsystems[subsys];
    
        SimTK_INDEXCHECK(index,(int)ss.discreteInfo.size(),
                         "StateImpl::updDiscreteVariable()");
        DiscreteVarInfo& dv = ss.discreteInfo[index];

        if (dv.getAllocationStage() > Stage::Topology) {
            SimTK_STAGECHECK_GE(getSubsystemStage(subsys), 
                Stage::Model, "StateImpl::updDiscreteVariable()");
        }
    
        invalidateAll(dv.getInvalidatedStage());

        // Invalidate the auto-update entry, if any.
        const CacheEntryIndex cx = dv.getAutoUpdateEntry();
        if (cx.isValid()) {
            CacheEntryInfo& ce = ss.cacheInfo[cx];
            ce.invalidate();
        }
    
        // We're now marking this variable as having been updated at the 
        // current time.
        return dv.updValue(t);
    }
    
    Stage getCacheEntryAllocationStage
       (SubsystemIndex subsys, CacheEntryIndex index) const {
        const PerSubsystemInfo& ss = subsystems[subsys];
        SimTK_INDEXCHECK(index,(int)ss.cacheInfo.size(),
            "StateImpl::getCacheEntryAllocationStage()");
        const CacheEntryInfo& ce = ss.cacheInfo[index];
        return ce.getAllocationStage();
    } 

    // Stage >= ce.stage
    // This method gets called a lot, so make it fast in Release mode.
    const AbstractValue& 
    getCacheEntry(SubsystemIndex subsys, CacheEntryIndex index) const {
        const PerSubsystemInfo& ss = subsystems[subsys];
    
        SimTK_INDEXCHECK(index,(int)ss.cacheInfo.size(),
                         "StateImpl::getCacheEntry()");
        const CacheEntryInfo& ce = ss.cacheInfo[index];

        // These two unconditional checks are somewhat expensive; I'm leaving
        // them in though because they catch so many user errors.
        // (sherm 20130222).
        SimTK_STAGECHECK_GE_ALWAYS(ss.currentStage, 
            ce.getDependsOnStage(), "StateImpl::getCacheEntry()");

        if (ss.currentStage < ce.getComputedByStage()) {
            const StageVersion currentDependsOnVersion = 
                ss.stageVersions[ce.getDependsOnStage()];
            const StageVersion lastCacheVersion = 
                ce.getVersionWhenLastComputed();

            if (lastCacheVersion != currentDependsOnVersion) {
                SimTK_THROW4(Exception::CacheEntryOutOfDate,
                    ss.currentStage, ce.getDependsOnStage(), 
                    currentDependsOnVersion, lastCacheVersion);
            }
        }

        // If we get here then we're either past the "computed by" stage, or
        // we're past "depends on" with an explicit validation having been made 
        // at the "depends on" stage's current version.

        return ce.getValue();
    }
    
    // You can access a cache entry for update any time after it has been
    // allocated. This does not affect the stage.
    AbstractValue& 
    updCacheEntry(SubsystemIndex subsys, CacheEntryIndex index) const {
        const PerSubsystemInfo& ss = subsystems[subsys];
    
        SimTK_INDEXCHECK(index,(int)ss.cacheInfo.size(),
                         "StateImpl::updCacheEntry()");
        CacheEntryInfo& ce = ss.cacheInfo[index];
    
        return ce.updValue();
    }

    bool isCacheValueRealized(SubsystemIndex subx, CacheEntryIndex cx) const {
        const PerSubsystemInfo& ss = subsystems[subx];
        SimTK_INDEXCHECK(cx,(int)ss.cacheInfo.size(),
                         "StateImpl::isCacheValueRealized()");
        const CacheEntryInfo& ce = ss.cacheInfo[cx];
        return ce.isCurrent(getSubsystemStage(subx), 
                            getSubsystemStageVersions(subx));
    }
    void markCacheValueRealized(SubsystemIndex subx, CacheEntryIndex cx) const {
        const PerSubsystemInfo& ss = subsystems[subx];
        SimTK_INDEXCHECK(cx,(int)ss.cacheInfo.size(),
                         "StateImpl::markCacheValueRealized()");
        CacheEntryInfo& ce = ss.cacheInfo[cx];
    
        // If this cache entry depends on anything, it can't be valid unless 
        // we're at least *working* on its depends-on stage, meaning the current
        // stage would have to be the one before that. The depends-on stage is 
        // required to be at least Stage::Topology, so its prev() stage exists.
        SimTK_STAGECHECK_GE(getSubsystemStage(subx), 
                            ce.getDependsOnStage().prev(), 
                            "StateImpl::markCacheValueRealized()");

        ce.markAsComputed(getSubsystemStageVersions(subx));
    }

    void markCacheValueNotRealized
       (SubsystemIndex subx, CacheEntryIndex cx) const {
        const PerSubsystemInfo& ss = subsystems[subx];
        SimTK_INDEXCHECK(cx,(int)ss.cacheInfo.size(),
            "StateImpl::markCacheValueNotRealized()");
        CacheEntryInfo& ce = ss.cacheInfo[cx];

        ce.invalidate();
    }

    StageVersion getSystemTopologyStageVersion() const
    {   return systemStageVersions[Stage::Topology]; }

    void setSystemTopologyStageVersion(StageVersion topoVersion)
    {   assert(topoVersion>0); 
        systemStageVersions[Stage::Topology]=topoVersion; }

    // Capture the stage versions only for currently-realized stages.
    void getSystemStageVersions(Array_<StageVersion>& versions) const {
        versions.resize(currentSystemStage+1);
        for (int i=0; i <= currentSystemStage; ++i)
            versions[i] = systemStageVersions[i];
    }

    // If the current state is realized at least as high as the previous one, 
    // then report Stage::Infinity if all of those stage versions match.
    // Otherwise report either the first mismatch or the first now-invalid
    // stage if lower stages match.
    Stage getLowestSystemStageDifference
       (const Array_<StageVersion>& prevVersions) const {
        const int nRealizedBefore = (int)prevVersions.size();
        const int nRealizedNow    = (int)currentSystemStage+1; // count from 0
        const int nRealizedBoth   = std::min(nRealizedBefore,nRealizedNow);

        // First check the stages both had in common.
        Stage g=Stage::Topology; 
        for (; g < nRealizedBoth; ++g)
            if (systemStageVersions[g] != prevVersions[g])
                return g;

        // All stages that were valid before and now have identical versions.
        // If that's all there was before then nothing has changed.
        return nRealizedNow >= nRealizedBefore ? Stage::Infinity
                                               : g; // 1st unrealized stage
    }

    void autoUpdateDiscreteVariables();
    
    String toString() const;    
    String cacheToString() const;

private:
    void initializeStageVersions() {
        for (int i=0; i < Stage::NValid; ++i)
            systemStageVersions[i] = 1; // never 0
    }

private:
        // Shared global resource State variables //

    // We consider time t to be a state variable allocated at Topology stage,
    // with its "invalidated" stage Stage::Time. The value of t is NaN in an 
    // Empty State, and is initialized to zero when the System stage advances
    // to Stage::Topology (i.e., when the System is realized to stage Topology).
    Real            t;

    // The continuous state variables are allocated at Model stage, and given
    // their specified initial values when the System stage advances to
    // Stage::Model (i.e., when the System is realized to Model stage).
    Vector          y; // All the continuous state variables together {q,u,z}

        // These are views into y.
    Vector          q; // Stage::Position continuous variables
    Vector          u; // Stage::Velocity continuous variables
    Vector          z; // Stage::Dynamics continuous variables

        // These are not views; there are no qWeights.
    Vector          uWeights; // scaling for u
    Vector          zWeights; // scaling for z

        // These are Instance stage state variables.
    Vector          qerrWeights; // Scaling for perrs
    Vector          uerrWeights; // Scaling for pdoterrs and verrs

        // Shared global resource Cache entries //

    // This is the System's currently highest-valid Stage.
    mutable Stage        currentSystemStage;
    // This contains a counter for each system stage which is bumped each
    // time that stage is invalidated. This allows detection of a state
    // that has been changed even after a subsequent realization. The
    // Topology stage entry should match the System's Topology version.
    mutable StageVersion systemStageVersions[Stage::NValid];

        // DIFFERENTIAL EQUATIONS

    // All the state derivatives taken together (qdot,udot,zdot)
    mutable Vector  ydot; 

    // These are views into ydot.
    mutable Vector  qdot;       // Stage::Velocity
    mutable Vector  udot;       // Stage::Acceleration
    mutable Vector  zdot;       // Stage::Acceleration

    // This is an independent cache entry.
    mutable Vector  qdotdot;    // Stage::Acceleration

        // ALGEBRAIC EQUATIONS

    mutable Vector  yerr;        // All constraint errors together (qerr,uerr)
    mutable Vector  udoterr;     // Stage::Acceleration (Index 1 constraints)
    mutable Vector  multipliers; // Stage::Acceleration (Index 1 algebraic vars)

    // These are views into yerr.
    mutable Vector  qerr;       // Stage::Position (Index 3 constraints)
    mutable Vector  uerr;       // Stage::Velocity (Index 2 constraints)

        // DISCRETE EQUATIONS

    // All the event triggers together, ordered by stage.
    mutable Vector  allTriggers;

    // These are views into allTriggers.
    mutable Vector  triggers[Stage::NValid];
    

        // Subsystem support //

    // Subsystem 0 (always present) is for the System as a whole. Its name
    // and version are the System name and version.
    Array_<PerSubsystemInfo> subsystems;

    // Return true only if all subsystems are realized to at least Stage g.
    bool allSubsystemsAtLeastAtStage(Stage g) const {
        for (SubsystemIndex i(0); i < (int)subsystems.size(); ++i)
            if (subsystems[i].currentStage < g)
                return false;
        return true;
    }

};


//==============================================================================
//             INLINE IMPLEMENTATIONS OF STATE METHODS
//==============================================================================
// These mostly just forward to the StateImpl object.

inline void State::setNumSubsystems(int i) {
    updImpl().setNumSubsystems(i);
}
inline void State::initializeSubsystem
   (SubsystemIndex subsys, const String& name, const String& version) {
    updImpl().initializeSubsystem(subsys, name, version);
}

inline SubsystemIndex State::addSubsystem
   (const String& name, const String& version) {
    return updImpl().addSubsystem(name, version);
}
inline int State::getNumSubsystems() const {
    return getImpl().getNumSubsystems();
}
inline const String& State::getSubsystemName(SubsystemIndex subsys) const {
    return getImpl().getSubsystemName(subsys);
}
inline const String& State::getSubsystemVersion(SubsystemIndex subsys) const {
    return getImpl().getSubsystemVersion(subsys);
}
inline const Stage& State::getSubsystemStage(SubsystemIndex subsys) const {
    return getImpl().getSubsystemStage(subsys);
}
inline const Stage& State::getSystemStage() const {
    return getImpl().getSystemStage();
}
inline void State::invalidateAll(Stage stage) {
    updImpl().invalidateAll(stage);
}
inline void State::invalidateAllCacheAtOrAbove(Stage stage) const {
    getImpl().invalidateAllCacheAtOrAbove(stage);
}
inline void State::advanceSubsystemToStage(SubsystemIndex subsys, Stage stage) const {
    getImpl().advanceSubsystemToStage(subsys, stage);
}
inline void State::advanceSystemToStage(Stage stage) const {
    getImpl().advanceSystemToStage(stage);
}

inline StageVersion State::getSystemTopologyStageVersion() const 
{   return getImpl().getSystemTopologyStageVersion(); }

// Continuous state variables
inline QIndex State::allocateQ(SubsystemIndex subsys, const Vector& qInit) {
    return updImpl().allocateQ(subsys, qInit);
}
inline UIndex State::allocateU(SubsystemIndex subsys, const Vector& uInit) {
    return updImpl().allocateU(subsys, uInit);
}
inline ZIndex State::allocateZ(SubsystemIndex subsys, const Vector& zInit) {
    return updImpl().allocateZ(subsys, zInit);
}

// Constraint errors and multipliers
inline QErrIndex State::allocateQErr(SubsystemIndex subsys, int nqerr) const {
    return getImpl().allocateQErr(subsys, nqerr);
}
inline UErrIndex State::allocateUErr(SubsystemIndex subsys, int nuerr) const {
    return getImpl().allocateUErr(subsys, nuerr);
}
inline UDotErrIndex State::
allocateUDotErr(SubsystemIndex subsys, int nudoterr) const {
    return getImpl().allocateUDotErr(subsys, nudoterr);
}

// Event witness functions
inline EventTriggerByStageIndex State::
allocateEventTrigger(SubsystemIndex subsys, Stage stage, int nevent) const {
    return getImpl().allocateEventTrigger(subsys, stage, nevent);
}

// Discrete Variables
inline DiscreteVariableIndex State::
allocateDiscreteVariable(SubsystemIndex subsys, Stage stage, AbstractValue* v) {
    return updImpl().allocateDiscreteVariable(subsys, stage, v);
}
inline DiscreteVariableIndex State::
allocateAutoUpdateDiscreteVariable
   (SubsystemIndex subsys, Stage invalidates, AbstractValue* v,
    Stage updateDependsOn) {
    return updImpl().allocateAutoUpdateDiscreteVariable
       (subsys, invalidates, v, updateDependsOn); 
}

inline CacheEntryIndex State::
getDiscreteVarUpdateIndex
   (SubsystemIndex subsys, DiscreteVariableIndex index) const {
    return getImpl().getDiscreteVarUpdateIndex(subsys, index);
}
inline Stage State::
getDiscreteVarAllocationStage
   (SubsystemIndex subsys, DiscreteVariableIndex index) const {
    return getImpl().getDiscreteVarAllocationStage(subsys, index);
}
inline Stage State::
getDiscreteVarInvalidatesStage
   (SubsystemIndex subsys, DiscreteVariableIndex index) const {
    return getImpl().getDiscreteVarInvalidatesStage(subsys, index);
}

inline const AbstractValue& State::
getDiscreteVariable(SubsystemIndex subsys, DiscreteVariableIndex index) const {
    return getImpl().getDiscreteVariable(subsys, index);
}
inline Real State::
getDiscreteVarLastUpdateTime
   (SubsystemIndex subsys, DiscreteVariableIndex index) const {
    return getImpl().getDiscreteVarLastUpdateTime(subsys, index);
}
inline const AbstractValue& State::
getDiscreteVarUpdateValue
   (SubsystemIndex subsys, DiscreteVariableIndex index) const {
    return getImpl().getDiscreteVarUpdateValue(subsys, index);
}
inline AbstractValue& State::
updDiscreteVarUpdateValue
   (SubsystemIndex subsys, DiscreteVariableIndex index) const {
    return getImpl().updDiscreteVarUpdateValue(subsys, index);
}
inline bool State::
isDiscreteVarUpdateValueRealized
   (SubsystemIndex subsys, DiscreteVariableIndex index) const {
    return getImpl().isDiscreteVarUpdateValueRealized(subsys, index);
}
inline void State::
markDiscreteVarUpdateValueRealized
   (SubsystemIndex subsys, DiscreteVariableIndex index) const {
    getImpl().markDiscreteVarUpdateValueRealized(subsys, index);
}

inline AbstractValue& State::
updDiscreteVariable
   (SubsystemIndex subsys, DiscreteVariableIndex index) {
    return updImpl().updDiscreteVariable(subsys, index);
}
inline void State::
setDiscreteVariable
   (SubsystemIndex i, DiscreteVariableIndex index, const AbstractValue& v) {
    updDiscreteVariable(i,index) = v;
}

// Cache Entries
inline CacheEntryIndex State::
allocateCacheEntry(SubsystemIndex subsys, Stage dependsOn, Stage computedBy, 
                   AbstractValue* v) const {
    return getImpl().allocateCacheEntry(subsys, dependsOn, computedBy, v);
}

inline Stage State::
getCacheEntryAllocationStage(SubsystemIndex subsys, CacheEntryIndex index) const {
    return getImpl().getCacheEntryAllocationStage(subsys, index);
}
inline const AbstractValue& State::
getCacheEntry(SubsystemIndex subsys, CacheEntryIndex index) const {
    return getImpl().getCacheEntry(subsys, index);
}
inline AbstractValue& State::
updCacheEntry(SubsystemIndex subsys, CacheEntryIndex index) const {
    return getImpl().updCacheEntry(subsys, index);
}

inline bool State::
isCacheValueRealized(SubsystemIndex subx, CacheEntryIndex cx) const {
    return getImpl().isCacheValueRealized(subx, cx); 
}
inline void State::
markCacheValueRealized(SubsystemIndex subx, CacheEntryIndex cx) const {
    getImpl().markCacheValueRealized(subx, cx); 
}
inline void State::
markCacheValueNotRealized(SubsystemIndex subx, CacheEntryIndex cx) const {
    getImpl().markCacheValueNotRealized(subx, cx); 
}

// Global Resource Dimensions

inline int State::getNY() const {
    return getImpl().getNY();
}
inline int State::getNQ() const {
    return getImpl().getNQ();
}
inline SystemYIndex State::getQStart() const {
    return getImpl().getQStart();
}
inline int State::getNU() const {
    return getImpl().getNU();
}
inline SystemYIndex State::getUStart() const {
    return getImpl().getUStart();
}
inline int State::getNZ() const {
    return getImpl().getNZ();
}
inline SystemYIndex State::getZStart() const {
    return getImpl().getZStart();
}
inline int State::getNYErr() const {
    return getImpl().getNYErr();
}
inline int State::getNQErr() const {
    return getImpl().getNQErr();
}
inline SystemYErrIndex State::getQErrStart() const {
    return getImpl().getQErrStart();
}
inline int State::getNUErr() const {
    return getImpl().getNUErr();
}
inline SystemYErrIndex State::getUErrStart() const {
    return getImpl().getUErrStart();
}
inline int State::getNUDotErr() const {
    return getImpl().getNUDotErr();
}
inline int State::getNMultipliers() const {
    return getNUDotErr();
}
inline int State::getNEventTriggers() const {
    return getImpl().getNEventTriggers();
}
inline int State::getNEventTriggersByStage(Stage stage) const {
    return getImpl().getNEventTriggersByStage(stage);
}



// Per-Subsystem Dimensions
inline SystemQIndex State::getQStart(SubsystemIndex subsys) const {
    return getImpl().getQStart(subsys);
}
inline int State::getNQ(SubsystemIndex subsys) const {
    return getImpl().getNQ(subsys);
}
inline SystemUIndex State::getUStart(SubsystemIndex subsys) const {
    return getImpl().getUStart(subsys);
}
inline int State::getNU(SubsystemIndex subsys) const {
    return getImpl().getNU(subsys);
}
inline SystemZIndex State::getZStart(SubsystemIndex subsys) const {
    return getImpl().getZStart(subsys);
}
inline int State::getNZ(SubsystemIndex subsys) const {
    return getImpl().getNZ(subsys);
}
inline SystemQErrIndex State::getQErrStart(SubsystemIndex subsys) const {
    return getImpl().getQErrStart(subsys);
}
inline int State::getNQErr(SubsystemIndex subsys) const {
    return getImpl().getNQErr(subsys);
}
inline SystemUErrIndex State::getUErrStart(SubsystemIndex subsys) const {
    return getImpl().getUErrStart(subsys);
}
inline int State::getNUErr(SubsystemIndex subsys) const {
    return getImpl().getNUErr(subsys);
}
inline SystemUDotErrIndex State::getUDotErrStart(SubsystemIndex subsys) const {
    return getImpl().getUDotErrStart(subsys);
}
inline int State::getNUDotErr(SubsystemIndex subsys) const {
    return getImpl().getNUDotErr(subsys);
}
inline SystemMultiplierIndex State::getMultipliersStart(SubsystemIndex subsys) const {
    return SystemMultiplierIndex(getUDotErrStart(subsys));
}
inline int State::getNMultipliers(SubsystemIndex subsys) const {
    return getNUDotErr(subsys);
}
inline SystemEventTriggerByStageIndex State::
getEventTriggerStartByStage(SubsystemIndex subsys, Stage stage) const {
    return getImpl().getEventTriggerStartByStage(subsys, stage);
}
inline int State::
getNEventTriggersByStage(SubsystemIndex subsys, Stage stage) const {
    return getImpl().getNEventTriggersByStage(subsys, stage);
}



inline const Vector& State::
getEventTriggersByStage(SubsystemIndex subsys, Stage stage) const {
    return getImpl().getEventTriggersByStage(subsys, stage);
}
inline Vector& State::
updEventTriggersByStage(SubsystemIndex subsys, Stage stage) const {
    return getImpl().updEventTriggersByStage(subsys, stage);
}
inline const Vector& State::getQ(SubsystemIndex subsys) const {
    return getImpl().getQ(subsys);
}
inline const Vector& State::getU(SubsystemIndex subsys) const {
    return getImpl().getU(subsys);
}
inline const Vector& State::getZ(SubsystemIndex subsys) const {
    return getImpl().getZ(subsys);
}
inline const Vector& State::getUWeights(SubsystemIndex subsys) const {
    return getImpl().getUWeights(subsys);
}
inline const Vector& State::getZWeights(SubsystemIndex subsys) const {
    return getImpl().getZWeights(subsys);
}
inline Vector& State::updQ(SubsystemIndex subsys) {
    return updImpl().updQ(subsys);
}
inline Vector& State::updU(SubsystemIndex subsys) {
    return updImpl().updU(subsys);
}
inline Vector& State::updZ(SubsystemIndex subsys) {
    return updImpl().updZ(subsys);
}
inline Vector& State::updUWeights(SubsystemIndex subsys) {
    return updImpl().updUWeights(subsys);
}
inline Vector& State::updZWeights(SubsystemIndex subsys) {
    return updImpl().updZWeights(subsys);
}
inline const Vector& State::getQDot(SubsystemIndex subsys) const {
    return getImpl().getQDot(subsys);
}
inline const Vector& State::getUDot(SubsystemIndex subsys) const {
    return getImpl().getUDot(subsys);
}
inline const Vector& State::getZDot(SubsystemIndex subsys) const {
    return getImpl().getZDot(subsys);
}
inline const Vector& State::getQDotDot(SubsystemIndex subsys) const {
    return getImpl().getQDotDot(subsys);
}
inline Vector& State::updQDot(SubsystemIndex subsys) const {
    return getImpl().updQDot(subsys);
}
inline Vector& State::updUDot(SubsystemIndex subsys) const {
    return getImpl().updUDot(subsys);
}
inline Vector& State::updZDot(SubsystemIndex subsys) const {
    return getImpl().updZDot(subsys);
}
inline Vector& State::updQDotDot(SubsystemIndex subsys) const {
    return getImpl().updQDotDot(subsys);
}
inline const Vector& State::getQErr(SubsystemIndex subsys) const {
    return getImpl().getQErr(subsys);
}
inline const Vector& State::getUErr(SubsystemIndex subsys) const {
    return getImpl().getUErr(subsys);
}
inline const Vector& State::getQErrWeights(SubsystemIndex subsys) const {
    return getImpl().getQErrWeights(subsys);
}
inline const Vector& State::getUErrWeights(SubsystemIndex subsys) const {
    return getImpl().getUErrWeights(subsys);
}
inline const Vector& State::getUDotErr(SubsystemIndex subsys) const {
    return getImpl().getUDotErr(subsys);
}
inline const Vector& State::getMultipliers(SubsystemIndex subsys) const {
    return getImpl().getMultipliers(subsys);
}
inline Vector& State::updQErr(SubsystemIndex subsys) const {
    return getImpl().updQErr(subsys);
}
inline Vector& State::updUErr(SubsystemIndex subsys) const {
    return getImpl().updUErr(subsys);
}
inline Vector& State::updQErrWeights(SubsystemIndex subsys) {
    return updImpl().updQErrWeights(subsys);
}
inline Vector& State::updUErrWeights(SubsystemIndex subsys) {
    return updImpl().updUErrWeights(subsys);
}
inline Vector& State::updUDotErr(SubsystemIndex subsys) const {
    return getImpl().updUDotErr(subsys);
}
inline Vector& State::updMultipliers(SubsystemIndex subsys) const {
    return getImpl().updMultipliers(subsys);
}

inline SystemEventTriggerIndex State::
getEventTriggerStartByStage(Stage stage) const {
    return getImpl().getEventTriggerStartByStage(stage);
}

inline const Vector& State::getEventTriggers() const {
    return getImpl().getEventTriggers();
}
inline const Vector& State::getEventTriggersByStage(Stage stage) const {
    return getImpl().getEventTriggersByStage(stage);
}

inline Vector& State::updEventTriggers() const {
    return getImpl().updEventTriggers();
}
inline Vector& State::updEventTriggersByStage(Stage stage) const {
    return getImpl().updEventTriggersByStage(stage);
}

inline const Real& State::getTime() const {
    return getImpl().getTime();
}
inline const Vector& State::getY() const {
    return getImpl().getY();
}
inline const Vector& State::getQ() const {
    return getImpl().getQ();
}
inline const Vector& State::getU() const {
    return getImpl().getU();
}
inline const Vector& State::getZ() const {
    return getImpl().getZ();
}
inline const Vector& State::getUWeights() const {
    return getImpl().getUWeights();
}
inline const Vector& State::getZWeights() const {
    return getImpl().getZWeights();
}
Real& State::updTime() {
    return updImpl().updTime();
}
inline Vector& State::updY() {
    return updImpl().updY();
}
inline void State::setTime(Real t) {
    updTime() = t;
}
inline void State::setY(const Vector& y) {
    updY() = y;
}
inline Vector& State::updQ() {
    return updImpl().updQ();
}
inline Vector& State::updU() {
    return updImpl().updU();
}
inline Vector& State::updZ() {
    return updImpl().updZ();
}
inline Vector& State::updUWeights() {
    return updImpl().updUWeights();
}
inline Vector& State::updZWeights() {
    return updImpl().updZWeights();
}
inline void State::setQ(const Vector& q) {
    updQ() = q;
}
inline void State::setU(const Vector& u) {
    updU() = u;
}
inline void State::setZ(const Vector& z) {
    updZ() = z;
}
inline const Vector& State::getYDot() const {
    return getImpl().getYDot();
}
inline const Vector& State::getQDot() const {
    return getImpl().getQDot();
}
inline const Vector& State::getZDot() const {
    return getImpl().getZDot();
}
inline const Vector& State::getUDot() const {
    return getImpl().getUDot();
}
inline const Vector& State::getQDotDot() const {
    return getImpl().getQDotDot();
}
inline Vector& State::updYDot() const {
    return getImpl().updYDot();
}
inline Vector& State::updQDot() const {
    return getImpl().updQDot();
}
inline Vector& State::updZDot() const {
    return getImpl().updZDot();
}
inline Vector& State::updUDot() const {
    return getImpl().updUDot();
}
inline Vector& State::updQDotDot() const {
    return getImpl().updQDotDot();
}
inline const Vector& State::getYErr() const {
    return getImpl().getYErr();
}
inline const Vector& State::getQErr() const {
    return getImpl().getQErr();
}
inline const Vector& State::getUErr() const {
    return getImpl().getUErr();
}
inline const Vector& State::getQErrWeights() const {
    return getImpl().getQErrWeights();
}
inline const Vector& State::getUErrWeights() const {
    return getImpl().getUErrWeights();
}
inline const Vector& State::getUDotErr() const {
    return getImpl().getUDotErr();
}
inline const Vector& State::getMultipliers() const {
    return getImpl().getMultipliers();
}
inline Vector& State::updYErr() const {
    return getImpl().updYErr();
}
inline Vector& State::updQErr() const {
    return getImpl().updQErr();
}
inline Vector& State::updUErr() const {
    return getImpl().updUErr();
}
inline Vector& State::updQErrWeights() {
    return updImpl().updQErrWeights();
}
inline Vector& State::updUErrWeights() {
    return updImpl().updUErrWeights();
}
inline Vector& State::updUDotErr() const {
    return getImpl().updUDotErr();
}
inline Vector& State::updMultipliers() const {
    return getImpl().updMultipliers();
}

inline void State::
setSystemTopologyStageVersion(StageVersion topoVersion)
{   return updImpl().setSystemTopologyStageVersion(topoVersion); }

inline void State::
getSystemStageVersions(Array_<StageVersion>& versions) const {
    return getImpl().getSystemStageVersions(versions); 
}
inline Stage State::
getLowestSystemStageDifference(const Array_<StageVersion>& prev) const {
    return getImpl().getLowestSystemStageDifference(prev); 
}
inline void State::autoUpdateDiscreteVariables() {
    updImpl().autoUpdateDiscreteVariables(); 
}

inline String State::toString() const {
    return getImpl().toString();
}
inline String State::cacheToString() const {
    return getImpl().cacheToString();
}


} // namespace SimTK

/** @endcond **/   // End of hiding from Doxygen

#endif // SimTK_SimTKCOMMON_STATE_IMPL_H_


