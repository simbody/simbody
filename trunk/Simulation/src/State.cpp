/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-9 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
 * Contributors: Peter Eastman                                                *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"
#include "SimTKcommon/internal/Event.h"
#include "SimTKcommon/internal/State.h"

#include <cassert>
#include <algorithm>
#include <ostream>
#include <set>

namespace SimTK {

// These static methods implement a specialized stacking mechanism for State
// resources that can be allocated a different stages, which we'll call an
// "allocation stack". A resource that was
// allocated at a later stage must be forgotten again when that stage is
// subsequently invalidated, and keeping their allocations stacked by
// stage allows that to be done efficiently.
//
// The method are templatized and expect the stacks to be in std::vectors
// of the same template. The template value must be a type that supports
// three methods (the tempate analog to virtual functions):
//      deepAssign()            a non-shallow assignment, i.e. clone the value
//      deepDestruct()          destroy any owned heap space
//      getAllocationStage()    return the stage being worked on when this was allocated
// The template type must otherwise support shallow copy semantics so that
// the std::vector can move them around without causing any heap activity.

// Clear the contents of an allocation stack, freeing up all associated heap space.
template <class T>
static void clearAllocationStack(std::vector<T>& stack) {
    for (int i=stack.size()-1; i >= 0; --i)
        stack[i].deepDestruct();
    stack.clear();
}

// Resize the given allocation stack, taking care to free the heap space if the size is reduced.
template <class T>
static void resizeAllocationStack(std::vector<T>& stack, int newSize) {
    assert(newSize >= 0);
    for (int i = stack.size()-1; i >= newSize; --i)
        stack[i].deepDestruct();
    stack.resize(newSize);
}

// Keep only those stack entries whose allocation stage is <= the supplied one.
template <class T>
static void popAllocationStackBackToStage(std::vector<T>& stack, const Stage& g) {
    unsigned newSize = stack.size();
    while (newSize > 0 && stack[newSize-1].getAllocationStage() > g)
        stack[--newSize].deepDestruct();
    stack.resize(newSize); 
}

// Make this allocation stack the same as the source, copying only through the given stage.
template <class T>
static void copyAllocationStackThroughStage(std::vector<T>& stack, const std::vector<T>& src, const Stage& g) {
    unsigned nVarsToCopy = src.size(); // assume we'll copy all
    while (nVarsToCopy && src[nVarsToCopy-1].getAllocationStage() > g)
        --nVarsToCopy;
    resizeAllocationStack(stack, nVarsToCopy);
    for (unsigned i=0; i < nVarsToCopy; ++i)
        stack[i].deepAssign(src[i]);
}

// These local classes
//      DiscreteVarInfo
//      CacheVarInfo
//      EventInfo
// contain the information needed for each discrete variable and cache entry,
// including a reference to its current value and everything needed to understand
// its dependencies, as well as information for each Event allocated by a
// Subsystem.
//
// These are intended as elements in an allocation stack as described above,
// so it is expected that they will get reaallocated, copied, and destructed during 
// allocation as the std::vector gets resized from time to time. However, the
// discrete variable and cache entry *values* must remain in a fixed location in 
// memory once allocated, because callers are permitted to retain references to these 
// values once they have been allocated. So pointers to the AbstractValues are kept
// in these objects, and only the pointers are copied around. That means the actual 
// value object will not be deleted by the destructor; be sure to do that explicitly
// in the higher-level destructor or you'll have a nasty leak.

//------------------------------------------------------------------------------
//                              MechanicalStateInfo
//
// Here we have a block of continuous, real-valued, second order state variables
// q (generalized positions) paired with first order state variables u 
// (generalized speeds) such that qdot = N(q)*u for some invertible block 
// diagonal matrix N. This block is allocated during realization of Topology
// or Model stage. The actual number nq of q's and nu <= nq of u's in all such 
// blocks must be known by the end of this Subsystem's realization of its Model 
// stage. Then:
//
//  System-level resources allocated upon realize(Model):
//      - nq slots in the global q, qdot, and qdotdot Vectors
//      - nu slots in the global u and udot Vectors
//
// The fact that these are continuous variables with time derivatives does not
// necessarily mean we have to use numerical integration to find their values.
// In some cases we may have closed forms for the integrals q and u also. In
// particular, prescribed motion (which provides udots as a function of t,q,
// and u) may also provide a formula for u=u(t,q), and if it does it may 
// additionally provide a formula for q=q(t). So we may need zero, one, or
// two integrations for any mobility u. By the end of the Subsystem's realization 
// of Instance stage, we must know the number nq_integ <= nq of q's and 
// nu_integ <= nu of u's that are defined by differential equations and thus 
// require integration of qdots and udots respectively. (A second order integrator 
// like Verlet may find q from qdotdot rather than qdot.) Then:
//
//  System-level resources allocated upon realize(Instance):
//      - nqInteg slots in the q partition of the yInteg Vector
//      - nuInteg slots in the u partition of the yInteg Vector
//
//------------------------------------------------------------------------------
class MechanicalStateInfo {
public:
    MechanicalStateInfo() : allocationStage(Stage::Empty) {}

    MechanicalStateInfo(Stage allocation)
    :   allocationStage(allocation)
    {   assert(isReasonable()); }

    // Default copy constructor, copy assignment, destructor are fine since there
    // is no heap object owned here.

    // These the the "virtual" methods required by template methods elsewhere.
    MechanicalStateInfo&   deepAssign(const MechanicalStateInfo& src) {return operator=(src);}
    void         deepDestruct() {}
    const Stage& getAllocationStage() const {return allocationStage;}
private:
    // This is fixed at construction.
    Stage allocationStage;

    // These are initialized at construction but can be changed during realization
    // of allocationStage+1.
    Vector  qInit; // generalized position (2nd order variables)
    Vector  uInit; // generalized speeds (1st order variables); nu <= nq

    // These are assigned at realization of allocationStage+1.
    SystemQIndex    globalQIndex;
    SystemUIndex    globalUIndex;

    // These are initialized to nq,nu during realization of allocationStage+1 but
    // can be changed during realization of allocationStage+2.
    int     nqInteg;    // number of q's to be calculated by integration
    int     nuInteg;    // number of u's to be calculated by integration

    //SystemQIntegIndex   globalQIntegIndex;
    //SystemUIntegIndex   globalUIntegIndex;


    bool isReasonable() const {
        return (    allocationStage==Stage::Topology
                 || allocationStage==Stage::Model);
    }
};


//------------------------------------------------------------------------------
//                              AuxiliaryStateInfo
//
// Here we have a block of continuous real-valued state variables z, allocated
// during realization of Topology or Model stage. The actual number nz of z's in 
// all such blocks must be known by the end of this Subsystem's realization of
// its Model stage. Then:
//
//  System-level resources allocated upon realize(Model):
//      - nz slots in the global z Vector
//      - nz slots in the global zdot Vector
//
// By the end of the Subsystem's realization of Instance stage, we must know
// the number nz_integ of the z's are defined by differential equations and
// thus require integration of zdots. Then:
//
//  System-level resources allocated upon realize(Instance):
//      - nz_integ slots in the z partition of the y_integ Vector
//
//------------------------------------------------------------------------------
class AuxiliaryStateInfo {
public:
    AuxiliaryStateInfo() : allocationStage(Stage::Empty) {}

    AuxiliaryStateInfo(Stage allocation)
    :   allocationStage(allocation)
    {   assert(isReasonable()); }

    // Default copy constructor, copy assignment, destructor are fine since there
    // is no heap object owned here.

    // These the the "virtual" methods required by template methods elsewhere.
    AuxiliaryStateInfo&   deepAssign(const AuxiliaryStateInfo& src) {return operator=(src);}
    void         deepDestruct() {}
    const Stage& getAllocationStage() const {return allocationStage;}
private:
    // This is fixed at construction.
    Stage allocationStage;

    // This is initialized at construction but can be changed during realization
    // of allocationStage+1.
    Vector  zInit; // nz and initial values for z's

    // This is assigned at realization of allocationStage+1.
    SystemZIndex    globalZIndex;

    // This is initialized to nz during realization of allocationStage+1 but
    // can be changed during realization of allocationStage+2.
    int     nzInteg;    // number of z's to be calculated by integration

    //SystemZIntegIndex   globalZIntegIndex;

    bool isReasonable() const {
        return (    allocationStage==Stage::Topology
                 || allocationStage==Stage::Model);
    }
};

//------------------------------------------------------------------------------
//                              DiscreteVarInfo
//------------------------------------------------------------------------------
class DiscreteVarInfo {
public:
    DiscreteVarInfo()
    :   allocationStage(Stage::Empty), invalidatedStage(Stage::Empty),
        value(0), timeLastUpdated(NaN) {}

    DiscreteVarInfo(Stage allocation, Stage invalidated, AbstractValue* v,
                    CacheEntryIndex autoUpdate = CacheEntryIndex())
    :   allocationStage(allocation), invalidatedStage(invalidated), value(v),
        autoUpdateEntry(autoUpdate), timeLastUpdated(NaN) 
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

    // Exchange value pointers.
    void swapValue(AbstractValue*& other) {std::swap(value, other);}

    const AbstractValue& getValue() const {assert(value); return *value;}
    Real                 getTimeLastUpdated() const {assert(value); return timeLastUpdated;}
    AbstractValue&       updValue(Real updTime)
    {   assert(value); assert(updTime >= 0); 
        timeLastUpdated=updTime; return *value; }

    const Stage&    getInvalidatedStage() const {return invalidatedStage;}
    CacheEntryIndex getAutoUpdateEntry()  const {return autoUpdateEntry;}

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
             && (value != 0)
             && (timeLastUpdated >= 0 || isNaN(timeLastUpdated)); }
};

//------------------------------------------------------------------------------
//                              CacheEntryInfo
//------------------------------------------------------------------------------
class CacheEntryInfo {
public:
    CacheEntryInfo()
    :   allocationStage(Stage::Empty), dependsOnStage(Stage::Empty), computedByStage(Stage::Empty),
        value(0), versionWhenLastComputed(-1) {}

    CacheEntryInfo(Stage allocation, Stage dependsOn, Stage computedBy, AbstractValue* v)
    :   allocationStage(allocation), dependsOnStage(dependsOn), computedByStage(computedBy),
        value(v), versionWhenLastComputed(0) 
    {   assert(isReasonable()); }

    bool isValid(const Stage& current, const StageVersion versions[]) const 
    {   if (current >= computedByStage) return true;
        if (current <  dependsOnStage)  return false;
        return versions[dependsOnStage] == versionWhenLastComputed;}

    // These affect only the explicit validity flag which does not fully
    // determine validity; see isValid() above.
    void invalidate() {versionWhenLastComputed = -1;}
    void markAsValid(const StageVersion versions[])
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

    // Exchange value pointers.
    void swapValue(AbstractValue*& other) {std::swap(value, other);}
    const AbstractValue& getValue() const {assert(value); return *value;}
    AbstractValue&       updValue()       {assert(value); return *value;}

    const Stage&          getDependsOnStage()  const {return dependsOnStage;}
    const Stage&          getComputedByStage() const {return computedByStage;}
    DiscreteVariableIndex getAssociatedVar()   const {return associatedVar;}
private:
    // These are fixed at construction.
    Stage                   allocationStage;
    Stage                   dependsOnStage;
    Stage                   computedByStage;
    DiscreteVariableIndex   associatedVar;  // if this is an auto-update entry

    // These change at run time.
    AbstractValue*          value;
    StageVersion            versionWhenLastComputed;   // version of Stage dependsOn

    bool isReasonable() const
    {    return (   allocationStage==Stage::Topology
                 || allocationStage==Stage::Model
                 || allocationStage==Stage::Instance)
             && (computedByStage >= dependsOnStage)
             && (value != 0)
             && (versionWhenLastComputed >= 0); }
};

//--------------------------------- EventInfo ----------------------------------
//
//------------------------------------------------------------------------------
class EventInfo {
public:
    EventInfo() : allocationStage(Stage::Empty), triggerStage(Stage::Empty) {}

    EventInfo(Stage allocation, Event::Cause cause, Stage triggered=Stage::Empty)
    :   allocationStage(allocation), cause(cause), triggerStage(triggered)
    {   assert(isReasonable()); }

    // Default copy constructor, copy assignment, destructor are fine since there
    // is no heap object owned here.

    // These the the "virtual" methods required by template methods elsewhere.
    EventInfo&   deepAssign(const EventInfo& src) {return operator=(src);}
    void         deepDestruct() {}
    const Stage& getAllocationStage() const {return allocationStage;}
private:
    // These are fixed at construction.
    Stage                   allocationStage;
    Event::Cause            cause;
    Stage                   triggerStage;   // only if EventCause==Triggered

    // These are assigned and filled in when the System Stage is advanced to Instance.
    SystemEventIndex               globalEventId;          // unique Id for this Event
    SystemEventTriggerByStageIndex globalByStageTriggerId; // only if EventCause==Triggered

    bool isReasonable() const {
        return (    allocationStage==Stage::Topology
                 || allocationStage==Stage::Model
                 || allocationStage==Stage::Instance)
            && cause.isValid()
            && (cause != Event::Cause::Triggered || triggerStage >= Stage::Time);
    }
};

//--------------------------- PerSubsystemInfo ---------------------------------
//
// This internal utility class is used to capture all the information needed for
// a single subsystem within the StateData.
//------------------------------------------------------------------------------
class PerSubsystemInfo {
public:
    PerSubsystemInfo() : currentStage(Stage::Empty)     {initialize();}
    PerSubsystemInfo(const String& n, const String& v) 
      : name(n), version(v), currentStage(Stage::Empty) {initialize();}

    // Everything will properly clean itself up except for the AbstractValues
    // stored in the discrete variable and cache entry arrays. Be sure to
    // delete those prior to allowing the arrays themselves to destruct.
    ~PerSubsystemInfo() {
        clearDiscreteVars();
        clearCache();
    }

    // Copy constructor copies all variables but cache only through
    // modeled stage. Note that this must be done in conjunction with
    // copying the whole state or our global resource indices will
    // be nonsense.
    PerSubsystemInfo(const PerSubsystemInfo& src) : currentStage(Stage::Empty) {
        initialize();
        copyFrom(src, Stage::Model);
    }

    PerSubsystemInfo& operator=(const PerSubsystemInfo& src) {
        if (&src != this)
            copyFrom(src, Stage::Model);
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
    void advanceToStage(Stage g) {
        assert(g > Stage::Empty);
        assert(currentStage == g.prev());

        // This validates whatever the current version number is of Stage g.
        currentStage = g;
    }


    void clearReferencesToModelStageGlobals() {
        qstart.invalidate(); ustart.invalidate(); zstart.invalidate();
        q.clear(); u.clear(); z.clear();
        qdot.clear(); udot.clear(); zdot.clear(); qdotdot.clear();
    }

    void clearReferencesToInstanceStageGlobals() {
        qerrstart.invalidate(); uerrstart.invalidate(); udoterrstart.invalidate();
        qerr.clear(); uerr.clear(); udoterr.clear(); multipliers.clear();

        for (int j=0; j<Stage::NValid; ++j) {
            triggerstart[j].invalidate();
            triggers[j].clear();
        }
    }


    String name;
    String version;

        // DEFINITIONS //

    // Discrete variables can be defined (allocated) during realization
    // of Topology or Model stages. Cache entries and Events can be defined during
    // realization of Topology, Model, or Instance stages. No further
    // allocations are allowed. Then, when one of these stages is invalidated,
    // all the definitions that occurred during realization of that stage
    // must be forgotten.
    // 
    // To do that the discrete variables and cache entries are stored
    // in arrays which are really stacks, with definitions pushed onto the ends
    // as the stage is advanced and popped off the ends as the stage is reduced.

    // Topology and Model stage definitions.
    std::vector<MechanicalStateInfo>    quInfo;
    std::vector<AuxiliaryStateInfo>     zInfo;
    std::vector<DiscreteVarInfo>        discreteInfo;

    // Topology, Model, and Instance stage definitions.
    //mutable std::vector<ConstraintInfo> constraintInfo;
    mutable std::vector<EventInfo>      eventInfo;
    mutable std::vector<CacheEntryInfo> cacheInfo;

        // AGGREGATE GLOBAL RESOURCE NEEDS //
    
    // These accumulate default values for this subsystem's use of shared
    // global state variables. After the System is advanced to Stage::Model,
    // the state will allocate those globals and copy these initial
    // values into them. The lengths of these Vectors define the 
    // needs of this Subsystem.
    Vector qInit, uInit, zInit;

    // For constraints we need just lengths (nmultipliers==nudoterr).
    int nqerr, nuerr, nudoterr;

    // For event trigger functions we need just lengths, for each stage.
    int ntriggers[Stage::NValid];

        // GLOBAL RESOURCE ALLOCATIONS //

    // These are our own private views into partitions of the global
    // state and cache entries of the same names. The State will assign
    // contiguous blocks to this subsystem when the *System* stage is raised
    // to Model or Instance stage, and they are invalidated whenever that 
    // stage is invalidated. The starting indices are filled in here at 
    // the time the views are built.

    // Model stage global resources and views into them.
    SystemQIndex                    qstart;
    SystemUIndex                    ustart;
    SystemZIndex                    zstart;
    Vector q, u, z;
    mutable Vector qdot, udot, zdot, qdotdot;

    // Instance stage global resources and views into them.
    // Note that multipliers just use the same indices as udoterr.
    SystemQErrIndex                 qerrstart;
    SystemUErrIndex                 uerrstart;
    SystemUDotErrIndex              udoterrstart;
    SystemEventTriggerByStageIndex  triggerstart[Stage::NValid];
    mutable Vector qerr, uerr;
    mutable Vector udoterr, multipliers; // same size and partioning
    mutable Vector triggers[Stage::NValid];

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
    // This is for use in constructors and for resetting an existing State into its
    // post-construction condition.
    void initialize() {
        nqerr = nuerr = nudoterr= 0;
        qstart.invalidate(); ustart.invalidate(); zstart.invalidate();
        qerrstart.invalidate(); uerrstart.invalidate(); udoterrstart.invalidate();
        for (int j=0; j<Stage::NValid; ++j) {
            ntriggers[j] = 0;
            triggerstart[j].invalidate();
            stageVersions[j] = 1;
        }
        currentStage = Stage::Empty;
    }

    // Manage allocation stacks for DiscreteVars, CacheEntries, and Events.

    void clearDiscreteVars()                        {clearAllocationStack(discreteInfo);}
    void clearCache()                               {clearAllocationStack(cacheInfo);}
    void clearEvents()                              {clearAllocationStack(eventInfo);}

    void resizeDiscreteVars        (int newSize)    {resizeAllocationStack(discreteInfo,newSize);}
    void resizeCache               (int newSize)    {resizeAllocationStack(cacheInfo,newSize);}
    void resizeEvents              (int newSize)    {resizeAllocationStack(eventInfo,newSize);}

    void popDiscreteVarsBackToStage(const Stage& g) {popAllocationStackBackToStage(discreteInfo,g);}
    void popCacheBackToStage       (const Stage& g) {popAllocationStackBackToStage(cacheInfo,g);}
    void popEventsBackToStage      (const Stage& g) {popAllocationStackBackToStage(eventInfo,g);}

    void copyDiscreteVarsThroughStage(const std::vector<DiscreteVarInfo>& src, const Stage& g)
    {   copyAllocationStackThroughStage(discreteInfo, src, g); }
    void copyCacheThroughStage(const std::vector<CacheEntryInfo>& src, const Stage& g)
    {   copyAllocationStackThroughStage(cacheInfo, src, g); }
    void copyEventsThroughStage(const std::vector<EventInfo>& src, const Stage& g)
    {   copyAllocationStackThroughStage(eventInfo, src, g); }

    // Restore this subsystem to the way it last was at realize(g) for a given Stage g; 
    // that is, invalidate all stages > g. Allocations will be forgotten as Instance, 
    // Model, and Topology stages are invalidated.
    void restoreToStage(Stage g) {
        if (currentStage <= g)
			return;

        if (g < Stage::Instance) {
			clearReferencesToInstanceStageGlobals();

            // TODO: this assumes that all constraint error and event trigger
            // allocations are performed at Stage::Instance. Should allow
            // them to be done earlier also.
            nqerr = nuerr = nudoterr = 0;
            for (int i = 0; i < Stage::NValid; i++)
                ntriggers[i] = 0;
        }

        if (g < Stage::Model) {
            clearReferencesToModelStageGlobals();

            // TODO: this assumes that all continuous variable
            // allocations are performed at Stage::Model. Should allow
            // them to be done earlier also.
            qInit.clear(); uInit.clear(); zInit.clear();
        }

        if (g == Stage::Empty) {
            // Throw out everything, reset stage versions to 1. Leave
            // name and version alone.
            clearDiscreteVars(); clearCache(); clearEvents();
            initialize();
            return;
        }

        // Backup all the allocation stacks.
        popDiscreteVarsBackToStage(g);
        popCacheBackToStage(g);
        popEventsBackToStage(g);

        // Raise the version number for every stage that we're invalidating.
        for (int i=currentStage; i > g; --i)
            stageVersions[i]++;
        currentStage = g;
    }

    // Utility which makes "this" a copy of the source subsystem exactly as it
    // was after being realized to stage maxStage. If maxStage >= Model then
    // all the subsystem-private state variables will be copied, but only
    // cached computations up through maxStage come through. We clear
    // our references to global variables regardless -- those will have to
    // be repaired at the System (State global) level.
    void copyFrom(const PerSubsystemInfo& src, Stage maxStage) {
        const Stage targetStage = std::min<Stage>(src.currentStage, maxStage);

        // Forget any references to global resources.
        clearReferencesToInstanceStageGlobals();
        clearReferencesToModelStageGlobals();

        // Make sure the destination state doesn't have anything past targetStage.
        restoreToStage(targetStage);

        name     = src.name;
        version  = src.version;
        copyDiscreteVarsThroughStage(src.discreteInfo, targetStage);
        copyCacheThroughStage(src.cacheInfo, targetStage);
        copyEventsThroughStage(src.eventInfo, targetStage);

        if (targetStage >= Stage::Model) {
            // Get the specifications for global state variables.
            qInit = src.qInit;
            uInit = src.uInit;
            zInit = src.zInit;
            if (targetStage >= Stage::Instance) {
                // Get the specifications for global cache entries.
                nqerr = src.nqerr;
                nuerr = src.nuerr;
                nudoterr = src.nudoterr;
                for (int j=0; j<Stage::NValid; ++j)
                    ntriggers[j] = src.ntriggers[j];
            }
        }

        // Set stage versions so that any cache entries we copied can still
        // be valid if they were valid in the source.
        for (int i=0; i<=targetStage; ++i)
            stageVersions[i] = src.stageVersions[i];

        // Subsystem stage should now match what we copied.
        currentStage = targetStage;
    }
};

class StateData {
public:
    StateData() 
    :   t(NaN), currentSystemStage(Stage::Empty), lowestModifiedSystemStage(Stage::Topology),
        myHandle(0) {} 

    // We'll do the copy constructor and assignment explicitly here
    // to get tight control over what's allowed, and to make sure
    // we don't copy the handle pointer.
    StateData(const StateData& src)
    :   myHandle(0), currentSystemStage(Stage::Empty), lowestModifiedSystemStage(Stage::Topology)
    {
        subsystems = src.subsystems;
        if (src.currentSystemStage >= Stage::Topology) {
            advanceSystemToStage(Stage::Topology);
            t = src.t;
            if (src.currentSystemStage >= Stage::Model) {
                advanceSystemToStage(Stage::Model);
                // careful -- don't allow reallocation
                y = src.y;
            }
        }
    }

    StateData& operator=(const StateData& src) {
        if (&src == this) return *this;
        invalidateJustSystemStage(Stage::Topology);
        for (SubsystemIndex i(0); i<(int)subsystems.size(); ++i)
            subsystems[i].invalidateStageJustThisSubsystem(Stage::Topology);
        subsystems = src.subsystems;
        if (src.currentSystemStage >= Stage::Topology) {
            advanceSystemToStage(Stage::Topology);
            t = src.t;
            if (src.currentSystemStage >= Stage::Model) {
                advanceSystemToStage(Stage::Model);
                // careful -- don't allow reallocation
                y = src.y;
            }
        }
        // don't mess with the handle pointer!
        return *this;
    }

    ~StateData() {   // default destructor
    }

    // Copies all the variables but not the cache.
    StateData* clone() const {return new StateData(*this);}

    const Stage& getSystemStage() const {return currentSystemStage;}
    Stage&       updSystemStage() const {return currentSystemStage;} // mutable


    const PerSubsystemInfo& getSubsystem(int subsystem) const {
        SimTK_INDEXCHECK(0, subsystem, (int)subsystems.size(), "StateData::getSubsystem()");
        return subsystems[subsystem];
    }

    PerSubsystemInfo& updSubsystem(int subsystem) {
        SimTK_INDEXCHECK(0, subsystem, (int)subsystems.size(), "StateData::updSubsystem()");
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


    // Back up the System stage just before g if it thinks
    // it is already at g or beyond. Note that we may be backing up
    // over many stages here. Careful: invalidating the stage
    // for the system must also invalidate the same stage for all
    // the subsystems (because we trash the shared resource pool
    // here if we back up earlier than Stage::Model) but we don't
    // take care of that here. Also, you can't invalidate Stage::Empty.
    void invalidateJustSystemStage(Stage g) {
        assert(g > Stage::Empty);
        if (currentSystemStage < g)
            return;

        if (currentSystemStage >= Stage::Instance && Stage::Instance >= g) {
            // We are "uninstancing" this State. Trash all the shared
            // cache entries that are allocated at Instance stage.

            // First make sure no subsystem is looking at the
            // shared cache entries any more.
            for (SubsystemIndex i(0); i < (int)subsystems.size(); ++i)
                subsystems[i].clearReferencesToInstanceStageGlobals();

            // Next get rid of the global views of these cache entries.
            qerr.clear(); uerr.clear();             // yerr views
            for (int j=0; j<Stage::NValid; ++j)
                triggers[j].clear();                // event trigger views

            // Finally nuke the actual cache data.
            yerr.unlockShape();        yerr.clear();
            udoterr.unlockShape();     udoterr.clear();
            multipliers.unlockShape(); multipliers.clear();
            allTriggers.unlockShape(); allTriggers.clear();
        }
        if (currentSystemStage >= Stage::Model && Stage::Model >= g) {
            // We are "unmodeling" this State. Trash all the global
            // shared states & corresponding cache entries.

            // First make sure no subsystem is looking at the
            // global shared state any more.
            for (SubsystemIndex i(0); i < (int)subsystems.size(); ++i)
                subsystems[i].clearReferencesToModelStageGlobals();

            // Next get rid of the global views of these state variables
            // and corresponding cache entries.
            q.clear(); u.clear(); z.clear();            // y views
            qdot.clear(); udot.clear(); zdot.clear();   // ydot views

            // Finally nuke the actual state data.
            y.unlockShape();           y.clear(); 
            ydot.unlockShape();        ydot.clear(); 
            qdotdot.unlockShape();     qdotdot.clear();
        }
        if (currentSystemStage >= Stage::Topology && Stage::Topology >= g) {
            // We're invalidating the topology stage. Time is considered
            // a topology stage variable so needs to be invalidated here.
            t = NaN;
        }

        // Now record the lowest System Stage we have invalidated, and
        // set the current System Stage one lower than the one being invalidated.
        lowestModifiedSystemStage = std::min(lowestModifiedSystemStage, g);
        currentSystemStage        = g.prev();
    }

    // Advance the System stage from g-1 to g. It is a fatal error if
    // we're not already at g-1, and you can't advance to Stage::Empty.
    // Also, you can't advance the system to g unless ALL subsystems have
    // already gotten there.
    void advanceSystemToStage(Stage g) {
        assert(g > Stage::Empty);
        assert(currentSystemStage == g.prev());
        assert(allSubsystemsAtLeastAtStage(g));

        if (g == Stage::Topology) {
            // As the final "Topology" step, initialize time to 0 (it's NaN before this).
            t = 0;
        }
        else if (g == Stage::Model) {
            // We know the shared state pool sizes now. Allocate the
            // states and matching shared cache pools.
            int nq=0, nu=0, nz=0;

            // Count up all 
            for (SubsystemIndex i(0); i<(int)subsystems.size(); ++i) {
                nq += subsystems[i].qInit.size();
                nu += subsystems[i].uInit.size();
                nz += subsystems[i].zInit.size();
            }

            // Allocate the actual shared state variables & cache 
            // entries and make sure no one can accidentally change the size.
            y.resize(nq+nu+nz);             y.lockShape();
            ydot.resize(nq+nu+nz);          ydot.lockShape();
            qdotdot.resize(nq);             qdotdot.lockShape();

            // Allocate subviews of the shared state & cache entries.
            q.viewAssign(y(0,nq));
            u.viewAssign(y(nq,nu));
            z.viewAssign(y(nq+nu,nz));

            qdot.viewAssign(ydot(0,     nq));
            udot.viewAssign(ydot(nq,    nu));
            zdot.viewAssign(ydot(nq+nu, nz));

            // Now partition the global resources among the subsystems and copy
            // in the initial values for the state variables.
            SystemQIndex nxtq(0);
            SystemUIndex nxtu(0);
            SystemZIndex nxtz(0);

            for (SubsystemIndex i(0); i<(int)subsystems.size(); ++i) {
                PerSubsystemInfo& ss = subsystems[i];
                const int nq=ss.qInit.size(), nu=ss.uInit.size(), nz=ss.zInit.size();

                // Assign the starting indices.
                ss.qstart=nxtq; ss.ustart=nxtu; ss.zstart=nxtz;
 
                // Build the views.
                ss.q.viewAssign(q(nxtq, nq)); ss.q = ss.qInit;
                ss.qdot.viewAssign(qdot(nxtq, nq));
                ss.qdotdot.viewAssign(qdotdot(nxtq, nq));
                ss.u.viewAssign(u(nxtu, nu)); ss.u = ss.uInit;
                ss.udot.viewAssign(udot(nxtu, nu));
                ss.z.viewAssign(z(nxtz, nz)); ss.z = ss.zInit;
                ss.zdot.viewAssign(zdot(nxtz, nz));

                // Consume the slots.
                nxtq += nq; nxtu += nu; nxtz += nz;
            }
        }
        else if (g == Stage::Instance) {
            // We know the shared state pool sizes now. Allocate the
            // states and matching shared cache pools.
            int nqerr=0, nuerr=0, nudoterr=0, nAllTriggers=0;
            int ntriggers[Stage::NValid];
            for (int j=0; j<Stage::NValid; ++j)
                ntriggers[j] = 0;

            // Count up all 
            for (SubsystemIndex i(0); i<(int)subsystems.size(); ++i) {
                nqerr    += subsystems[i].nqerr;
                nuerr    += subsystems[i].nuerr;
                nudoterr += subsystems[i].nudoterr;

                for (int j=0; j<Stage::NValid; ++j)
                    ntriggers[j] += subsystems[i].ntriggers[j];
            }
            for (int j=0; j<Stage::NValid; ++j)
                nAllTriggers += ntriggers[j];

            // Allocate the actual shared state variables & cache 
            // entries and make sure no one can accidentally change the size.
            yerr.resize(nqerr+nuerr);         yerr.lockShape();
            udoterr.resize(nudoterr);         udoterr.lockShape();
            multipliers.resize(nudoterr);     multipliers.lockShape(); // same size as udoterr
            allTriggers.resize(nAllTriggers); allTriggers.lockShape();

            // Allocate subviews of the shared state & cache entries.

            qerr.viewAssign(yerr(0,     nqerr));
            uerr.viewAssign(yerr(nqerr, nuerr));

            int stageStart=0;
            for (int j=0; j<Stage::NValid; ++j) {
                triggers[j].viewAssign(allTriggers(stageStart, ntriggers[j]));
                stageStart += ntriggers[j];
            }

            // Now partition the global resources among the subsystems and copy
            // in the initial values for the state variables.
            SystemQErrIndex nxtqerr(0);
            SystemUErrIndex nxtuerr(0);
            SystemUDotErrIndex nxtudoterr(0);
            SystemEventTriggerByStageIndex nxttrigger[Stage::NValid];
            for (int j=0; j<Stage::NValid; ++j)
                nxttrigger[j] = SystemEventTriggerByStageIndex(0);

            for (SubsystemIndex i(0); i<(int)subsystems.size(); ++i) {
                PerSubsystemInfo& ss = subsystems[i];

                // Assign the starting indices.
                ss.qerrstart=nxtqerr; ss.uerrstart=nxtuerr; ss.udoterrstart=nxtudoterr;

                // Build the views.
                ss.qerr.viewAssign(qerr(nxtqerr, ss.nqerr));
                ss.uerr.viewAssign(uerr(nxtuerr, ss.nuerr));
                ss.udoterr.viewAssign(udoterr(nxtudoterr, ss.nudoterr));
                // multipliers have same partitioning as udoterr
                ss.multipliers.viewAssign(multipliers(nxtudoterr, ss.nudoterr));

                // Consume the slots.
                nxtqerr += ss.nqerr; nxtuerr += ss.nuerr; nxtudoterr += ss.nudoterr;

                // Same thing for event trigger slots, but by stage.
                for (int j=0; j<Stage::NValid; ++j) {
                    ss.triggerstart[j] = nxttrigger[j];
                    ss.triggers[j].viewAssign(triggers[j](nxttrigger[j], ss.ntriggers[j]));
                    nxttrigger[j] += ss.ntriggers[j];
                }

            }
        }

        // All cases fall through to here.

        currentSystemStage = g;
    }

    void            setMyHandle(StateRep& s) {myHandle = &s;}
    const StateRep& getMyHandle() const   {assert(myHandle); return *myHandle;}
    StateRep&       updMyHandle()         {assert(myHandle); return *myHandle;}
private:
    friend class StateRep;

        // Shared global resource State variables //

    // We consider time t to be a state variable allocated at Topology stage,
    // with its "invalidated" stage Stage::Time. The value of t is NaN in an Empty
    // State, and is initialized to zero when the System stage advances
    // to Stage::Topology (i.e., when the System is realized to stage Topology).
    Real            t;

    // The continuous state variables are allocated at Model stage, and given
    // their specified initial values when the System stage advances to
    // Stage::Model (i.e., when the System is realized to Model stage).
    Vector          y; // All the continuous state variables taken together {q,u,z}

        // These are views into y.
    Vector          q; // Stage::Position continuous variables
    Vector          u; // Stage::Velocity continuous variables
    Vector          z; // Stage::Dynamics continuous variables


        // Shared global resource Cache entries //

    // This is the System's currently highest-valid Stage.
    mutable Stage   currentSystemStage;

    // This tracks the lowest invalidated System Stage since the last time
    // it was reset, which is done by an external call. It is *never* higher
    // than currentSystemStage+1.
    mutable Stage   lowestModifiedSystemStage;

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

    mutable Vector  yerr;        // All constraint errors taken together (qerr,uerr)
    mutable Vector  udoterr;     // Stage::Acceleration (Index 1 constraints)
    mutable Vector  multipliers; // Stage::Acceleration (Index 1 algebraic variables)

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
    std::vector<PerSubsystemInfo> subsystems;

    // Return true only if all subsystems have been realized to at least Stage g.
    bool allSubsystemsAtLeastAtStage(Stage g) const {
        for (SubsystemIndex i(0); i < (int)subsystems.size(); ++i)
            if (subsystems[i].currentStage < g)
                return false;
        return true;
    }

private:
    StateRep* myHandle;
};


class StateRep {
public:
    StateRep() : data(new StateData()) {
        data->setMyHandle(*this);
    }
    
    
    // Restore state to default-constructed condition
    void clear() {
        delete data;
        data = new StateData();
        data->setMyHandle(*this);
    }
    
    ~StateRep() {
        if (data->myHandle == this)
            delete data;
        data=0;
    }
    
    // copy constructor
    StateRep(const StateRep& src) 
      : data(src.data->clone()) {
        data->setMyHandle(*this);
    }
    
    // This constructor creates restricted states.
    StateRep(const StateRep& src, const EnumerationSet<Stage>& restrictedStages, const std::set<SubsystemIndex>& restrictedSubsystems) 
      : data(src.data), restrictedStages(restrictedStages), restrictedSubsystems(restrictedSubsystems) {
    }
    
    // copy assignment
    StateRep& operator=(const StateRep& src) {
        if (&src == this) return *this;
        if (!data) {
            // we're defining this state here (if src is not empty)
            if (src.data) {
                data = src.data->clone();
                data->setMyHandle(*this);
            }
            return *this;
        }
    
        // Assignment or redefinition
        if (src.data) *data = *src.data;
        else {delete data; data=0;}
        return *this;
    }
    
    void setNSubsystems(int i) {
        assert(i >= 0);
        data->subsystems.clear();
        data->subsystems.resize(i);
    }
    
    void initializeSubsystem(SubsystemIndex i, const String& name, const String& version) {
        data->updSubsystem(i).name = name;
        data->updSubsystem(i).version = version;
    }
    
    
    SubsystemIndex addSubsystem(const String& name, const String& version) {
        const SubsystemIndex sx(data->subsystems.size());
        data->subsystems.push_back(PerSubsystemInfo(name,version));
        return sx;
    }
    
    int getNSubsystems() const {return (int)data->subsystems.size();}
    
    const String& getSubsystemName(SubsystemIndex subsys) const {
        return data->subsystems[subsys].name;
    }
    const String& getSubsystemVersion(SubsystemIndex subsys) const {
        return data->subsystems[subsys].version;
    }
    
    const Stage& getSystemStage() const {
        return data->getSystemStage();
    }
    
    const Stage& getSubsystemStage(SubsystemIndex subx) const {
        SimTK_ASSERT(data, "StateRep::getSubsystemStage(): no data"); // can't happen(?)
        return data->getSubsystemStage(subx);
    }
     
    const StageVersion* getSubsystemStageVersions(SubsystemIndex subx) const {
        SimTK_ASSERT(data, "StateRep::getSubsystemStageVersions(): no data"); // can't happen(?)
        return data->getSubsystemStageVersions(subx);
    }  

    // Make sure the stage is no higher than g-1 for *any* subsystem and
    // hence for the system stage also. TODO: this should be more selective.
    void invalidateAll(Stage g) const {
        SimTK_ASSERT(data, "StateRep::invalidateAll(): no data");
    
        data->invalidateJustSystemStage(g);
        for (SubsystemIndex i(0); i<(int)data->subsystems.size(); ++i)
            data->subsystems[i].invalidateStageJustThisSubsystem(g);
    }
    
    // Move the stage for a particular subsystem from g-1 to g. No other subsystems
    // are affected, nor the global system stage.
    void advanceSubsystemToStage(SubsystemIndex subsys, Stage g) const {
        SimTK_ASSERT(data, "StateRep::advanceSubsystemToStage(): no data");
    
        data->subsystems[subsys].advanceToStage(g);
        // We don't automatically advance the System stage even if this brings
        // ALL the subsystems up to stage g.
    }
    
    // Move the system stage from g-1 to g. Don't call this until ALL 
    // subsystem have been advanced to at least stage g.
    void advanceSystemToStage(Stage g) const {
        SimTK_ASSERT(data, "StateRep::advanceToStage(): no data");
    
        // Terrible things will happen if either of these conditions is not met:
        //   (1) the system is at stage g-1 now, AND
        //   (2) ALL subsystems have already been advanced to stage g.
        data->advanceSystemToStage(g);
    }
    
    // We don't expect State entry allocations to be performance critical so
    // we'll keep error checking on even in Release mode.
    
    QIndex allocateQ(SubsystemIndex subsys, const Vector& qInit) {
        SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), Stage::Model, "StateRep::allocateQ()");
        const QIndex nxt(data->subsystems[subsys].qInit.size());
        data->subsystems[subsys].qInit.resizeKeep(nxt + qInit.size());
        data->subsystems[subsys].qInit(nxt, qInit.size()) = qInit;
        return nxt;
    }
    
    UIndex allocateU(SubsystemIndex subsys, const Vector& uInit) {
        SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), Stage::Model, "StateRep::allocateU()");
        const UIndex nxt(data->subsystems[subsys].uInit.size());
        data->subsystems[subsys].uInit.resizeKeep(nxt + uInit.size());
        data->subsystems[subsys].uInit(nxt, uInit.size()) = uInit;
        return nxt;
    }
    ZIndex allocateZ(SubsystemIndex subsys, const Vector& zInit) {
        SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), Stage::Model, "StateRep::allocateZ()");
        const ZIndex nxt(data->subsystems[subsys].zInit.size());
        data->subsystems[subsys].zInit.resizeKeep(nxt + zInit.size());
        data->subsystems[subsys].zInit(nxt, zInit.size()) = zInit;
        return nxt;
    }
    
    QErrIndex allocateQErr(SubsystemIndex subsys, int nqerr) const {
        SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), Stage::Instance, "StateRep::allocateQErr()");
        const QErrIndex nxt(data->subsystems[subsys].nqerr);
        data->subsystems[subsys].nqerr += nqerr;
        return nxt;
    }
    UErrIndex allocateUErr(SubsystemIndex subsys, int nuerr) const {
        SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), Stage::Instance, "StateRep::al()");
        const UErrIndex nxt(data->subsystems[subsys].nuerr);
        data->subsystems[subsys].nuerr += nuerr;
        return nxt;
    }
    UDotErrIndex allocateUDotErr(SubsystemIndex subsys, int nudoterr) const {
        SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), Stage::Instance, "StateRep::allocateUDotErr()");
        const UDotErrIndex nxt(data->subsystems[subsys].nudoterr);
        data->subsystems[subsys].nudoterr += nudoterr;
        return nxt;
    }
    EventTriggerByStageIndex allocateEventTrigger(SubsystemIndex subsys, Stage g, int nt) const {
        SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), Stage::Instance, "StateRep::allocateEvent()");
        const EventTriggerByStageIndex nxt(data->subsystems[subsys].ntriggers[g]);
        data->subsystems[subsys].ntriggers[g] += nt;
        return nxt;
    }
    
    // Topology- and Model-stage State variables can only be added during construction; that is,
    // while stage <= Topology. Other entries can be added while stage < Model.
    DiscreteVariableIndex allocateDiscreteVariable(SubsystemIndex subsys, Stage invalidates, AbstractValue* vp) {
        SimTK_STAGECHECK_RANGE_ALWAYS(Stage(Stage::LowestRuntime).prev(), invalidates, Stage::HighestRuntime, 
            "StateRep::allocateDiscreteVariable()");
    
        const Stage maxAcceptable = (invalidates <= Stage::Model ? Stage::Empty : Stage::Topology);
        SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), 
            maxAcceptable.next(), "StateRep::allocateDiscreteVariable()");
    
        PerSubsystemInfo& ss = data->subsystems[subsys];
        const DiscreteVariableIndex nxt(ss.discreteInfo.size());
        ss.discreteInfo.push_back(DiscreteVarInfo(getSubsystemStage(subsys).next(),invalidates,vp));
        return nxt;
    }
    
    // Cache entries can be allocated while stage < Instance.
    CacheEntryIndex allocateCacheEntry(SubsystemIndex subsys, Stage g, AbstractValue* vp) const {
        SimTK_STAGECHECK_RANGE_ALWAYS(Stage(Stage::LowestRuntime).prev(), g, Stage::HighestRuntime, 
            "StateRep::allocateCacheEntry()");
        SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), 
            Stage::Instance, "StateRep::allocateCacheEntry()");

        const PerSubsystemInfo& ss = data->subsystems[subsys];
        const CacheEntryIndex nxt(ss.cacheInfo.size());
        ss.cacheInfo.push_back(CacheEntryInfo(getSubsystemStage(subsys).next(),g,g,vp)); // mutable
        return nxt;
    }
    
        // State dimensions for shared continuous variables.
    
    int getNY() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getNY()");
        return data->y.size();
    }
    
    SystemYIndex getQStart() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getQStart()");
        return SystemYIndex(0); // q's come first
    }
    int getNQ() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getNQ()");
        return data->q.size();
    }
    
    SystemYIndex getUStart() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getUStart()");
        return SystemYIndex(data->q.size()); // u's come right after q's
    }
    int getNU() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getNU()");
        return data->u.size();
    }
    
    SystemYIndex getZStart() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getZStart()");
        return SystemYIndex(data->q.size() + data->u.size()); // q,u, then z
    }
    int getNZ() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getNZ()");
        return data->z.size();
    }
    
    int getNYErr() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getNYErr()");
        return data->yerr.size();
    }
    
    SystemYErrIndex getQErrStart() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getQErrStart()");
        return SystemYErrIndex(0); // qerr's come first
    }
    int getNQErr() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getNQErr()");
        return data->qerr.size();
    }
    
    SystemYErrIndex getUErrStart() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getUErrStart()");
        return SystemYErrIndex(data->qerr.size()); // uerr's follow qerrs
    }
    int getNUErr() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getNUErr()");
        return data->uerr.size();
    }
    
    // UDot errors are independent of qerr & uerr.
    // This is used for multipliers also.
    int getNUDotErr() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getNUDotErr()");
        return data->udoterr.size();
    }
    
    int getNEventTriggers() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getNEventTriggers()");
        return data->allTriggers.size();
    }
    
    SystemEventTriggerIndex getEventTriggerStartByStage(Stage g) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getEventTriggerStartByStage()");
        int nxt = 0;
        for (int j=0; j<g; ++j)
            nxt += data->triggers[j].size();
        return SystemEventTriggerIndex(nxt); // g starts where g-1 leaves off
    }
    
    int getNEventTriggersByStage(Stage g) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getNEventTriggersByStage()");
        return data->triggers[g].size();
    }
    
        // Subsystem dimensions.
    
    SystemQIndex getQStart(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getQStart(subsys)");
        return data->getSubsystem(subsys).qstart;
    }
    int getNQ(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getNQ(subsys)");
        return data->getSubsystem(subsys).q.size();
    }
    
    SystemUIndex getUStart(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getUStart(subsys)");
        return data->getSubsystem(subsys).ustart;
    }
    int getNU(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getNU(subsys)");
        return data->getSubsystem(subsys).u.size();
    }
    
    SystemZIndex getZStart(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getZStart(subsys)");
        return data->getSubsystem(subsys).zstart;
    }
    int getNZ(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getNZ(subsys)");
        return data->getSubsystem(subsys).z.size();
    }
    
    SystemQErrIndex getQErrStart(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getQErrStart(subsys)");
        return data->getSubsystem(subsys).qerrstart;
    }
    int getNQErr(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getNQErr(subsys)");
        return data->getSubsystem(subsys).qerr.size();
    }
    
    SystemUErrIndex getUErrStart(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getUErrStart(subsys)");
        return data->getSubsystem(subsys).uerrstart;
    }
    int getNUErr(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getNUErr(subsys)");
        return data->getSubsystem(subsys).uerr.size();
    }
    
    // These are used for multipliers also.
    SystemUDotErrIndex getUDotErrStart(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getUDotErrStart(subsys)");
        return data->getSubsystem(subsys).udoterrstart;
    }
    int getNUDotErr(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getNUDotErr(subsys)");
        return data->getSubsystem(subsys).udoterr.size();
    }
    
    SystemEventTriggerByStageIndex getEventTriggerStartByStage(SubsystemIndex subsys, Stage g) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getEventTriggerStartByStage(subsys)");
        return data->getSubsystem(subsys).triggerstart[g];
    }
    
    int getNEventTriggersByStage(SubsystemIndex subsys, Stage g) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getNEventTriggersByStage(subsys)");
        return data->getSubsystem(subsys).triggers[g].size();
    }
    
        // Per-subsystem access to the global shared variables.
    
    const Vector& getQ(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getQ(subsys)");
        return data->getSubsystem(subsys).q;
    }
    const Vector& getU(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getU(subsys)");
        return data->getSubsystem(subsys).u;
    }
    const Vector& getZ(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getZ(subsys)");
        return data->getSubsystem(subsys).z;
    }
    
    const Vector& getQDot(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getQDot(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Velocity, "StateRep::getQDot(subsys)");
        return data->getSubsystem(subsys).qdot;
    }
    const Vector& getUDot(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getUDot(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Acceleration, "StateRep::getUDot(subsys)");
        return data->getSubsystem(subsys).udot;
    }
    const Vector& getZDot(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getZDot(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Dynamics, "StateRep::getZDot(subsys)");
        return data->getSubsystem(subsys).zdot;
    }
    const Vector& getQDotDot(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getQDotDot(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Acceleration, "StateRep::getQDotDot(subsys)");
        return data->getSubsystem(subsys).qdotdot;
    }
    
    Vector& updQ(SubsystemIndex subsys) {
        assert(data);
        checkCanModify(Stage::Position);
        checkCanModify(subsys);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::updQ(subsys)");
        invalidateAll(Stage::Position);
        return data->updSubsystem(subsys).q;
    }
    Vector& updU(SubsystemIndex subsys) {
        assert(data);
        checkCanModify(Stage::Velocity);
        checkCanModify(subsys);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::updU(subsys)");
        invalidateAll(Stage::Velocity);
        return data->updSubsystem(subsys).u;
    }
    Vector& updZ(SubsystemIndex subsys) {
        assert(data);
        checkCanModify(Stage::Dynamics);
        checkCanModify(subsys);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::updZ(subsys)");
        invalidateAll(Stage::Dynamics);
        return data->updSubsystem(subsys).z;
    }
    
        // These are mutable so the routines are const.
    
    Vector& updQDot(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::updQDot(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage(Stage::Velocity).prev(), "StateRep::updQDot(subsys)");
        return data->getSubsystem(subsys).qdot;
    }
    Vector& updUDot(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::updUDot(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage(Stage::Acceleration).prev(), "StateRep::updUDot(subsys)");
        return data->getSubsystem(subsys).udot;
    }
    Vector& updZDot(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::updZDot(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage(Stage::Dynamics).prev(), "StateRep::updZDot(subsys)");
        return data->getSubsystem(subsys).zdot;
    }
    Vector& updQDotDot(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::updQDotDot(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage(Stage::Acceleration).prev(), "StateRep::updQDotDot(subsys)");
        return data->getSubsystem(subsys).qdotdot;
    }
    
    
    const Vector& getQErr(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getQErr(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Position, "StateRep::getQErr(subsys)");
        return data->getSubsystem(subsys).qerr;
    }
    const Vector& getUErr(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getUErr(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Velocity, "StateRep::getUErr(subsys)");
        return data->getSubsystem(subsys).uerr;
    }
    const Vector& getUDotErr(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getUDotErr(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Acceleration, "StateRep::getUDotErr(subsys)");
        return data->getSubsystem(subsys).udoterr;
    }
    const Vector& getMultipliers(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getMultipliers(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Acceleration, "StateRep::getMultipliers(subsys)");
        return data->getSubsystem(subsys).multipliers;
    }
    
    const Vector& getEventTriggersByStage(SubsystemIndex subsys, Stage g) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getEventTriggersByStage(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), g, "StateRep::getEventTriggersByStage(subsys)");
        return data->getSubsystem(subsys).triggers[g];
    }
    
    Vector& updQErr(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::updQErr(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage(Stage::Position).prev(), "StateRep::updQErr(subsys)");
        return data->getSubsystem(subsys).qerr;
    }
    Vector& updUErr(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::updUErr(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage(Stage::Velocity).prev(), "StateRep::updUErr(subsys)");
        return data->getSubsystem(subsys).uerr;
    }
    Vector& updUDotErr(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::updUDotErr(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage(Stage::Acceleration).prev(), 
                            "StateRep::updUDotErr(subsys)");
        return data->getSubsystem(subsys).udoterr;
    }
    Vector& updMultipliers(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::updMultipliers(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage(Stage::Acceleration).prev(), 
                            "StateRep::updMultipliers(subsys)");
        return data->getSubsystem(subsys).multipliers;
    }
    Vector& updEventTriggersByStage(SubsystemIndex subsys, Stage g) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::updEventTriggersByStage(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), g.prev(), "StateRep::updEventTriggersByStage(subsys)");
        return data->getSubsystem(subsys).triggers[g];
    }
    
        // Direct access to the global shared state and cache entries.
        // Time is allocated in Stage::Topology, State in Stage::Model, and
        // Cache in Stage::Instance.
    
    const Real& getTime() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Topology, "StateRep::getTime()");
        return data->t;
    }
    
    const Vector& getY() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getY()");
        return data->y;
    }
    
    const Vector& getQ() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getQ()");
        return data->q;
    }
    
    const Vector& getU() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getU()");
        return data->u;
    }
    
    const Vector& getZ() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getZ()");
        return data->z;
    }
    
    
    // You can call these as long as stage >= allocation stage, but the
    // stage will be backed up if necessary to one stage prior to the invalidated stage.
    Real& updTime() {  // Back to Stage::Time-1
        assert(data);
        checkCanModify(Stage::Time);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Topology, "StateRep::updTime()");
        invalidateAll(Stage::Time);
        return data->t;
    }
    
    Vector& updY() {    // Back to Stage::Position-1
        assert(data);
        checkCanModifyY();
        checkCanModifyAnySubsystem();
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::updY()");
        invalidateAll(Stage::Position);
        return data->y;
    }
    
    Vector& updQ() {    // Stage::Position-1
        assert(data);
        checkCanModify(Stage::Position);
        checkCanModifyAnySubsystem();
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::updQ()");
        invalidateAll(Stage::Position);
        return data->q;
    }
    
    Vector& updU() {     // Stage::Velocity-1
        assert(data);
        checkCanModify(Stage::Velocity);
        checkCanModifyAnySubsystem();
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::updU()");
        invalidateAll(Stage::Velocity);
        return data->u;
    }
    
    Vector& updZ() {     // Stage::Dynamics-1
        assert(data);
        checkCanModify(Stage::Dynamics);
        checkCanModifyAnySubsystem();
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::updZ()");
        invalidateAll(Stage::Dynamics);
        return data->z;
    }
    
    // These cache entries you can get at their "computedBy" stages.
    const Vector& getYDot() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Acceleration, "StateRep::getYDot()");
        return data->ydot;
    }
    
    const Vector& getQDot() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Velocity, "StateRep::getQDot()");
        return data->qdot;
    }
    
    const Vector& getZDot() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Dynamics, "StateRep::getZDot()");
        return data->zdot;
    }
    
    const Vector& getUDot() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Acceleration, "StateRep::getUDot()");
        return data->udot;
    }
    
    const Vector& getQDotDot() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Acceleration, "StateRep::getQDotDot()");
        return data->qdotdot;
    }
    
    // Cache updates are allowed while realizing their "dependsOn" stages.
    Vector& updYDot() const {
        assert(data);
        checkCanModifyY();
        SimTK_STAGECHECK_GE(getSystemStage(), Stage(Stage::Acceleration).prev(), "StateRep::updYDot()");
        return data->ydot;
    }
    
    Vector& updQDot() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage(Stage::Velocity).prev(), "StateRep::updQDot()");
        return data->qdot;
    }
    
    Vector& updUDot() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage(Stage::Acceleration).prev(), "StateRep::updUDot()");
        return data->udot;
    }
    
    Vector& updZDot() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage(Stage::Dynamics).prev(), "StateRep::updZDot()");
        return data->zdot;
    }
    
    Vector& updQDotDot() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage(Stage::Acceleration).prev(), "StateRep::updQDotDot()");
        return data->qdotdot;
    }
    
    
    const Vector& getYErr() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Velocity, "StateRep::getYErr()");
        return data->yerr;
    }
    
    const Vector& getQErr() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Position, "StateRep::getQErr()");
        return data->qerr;
    }
    const Vector& getUErr() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Velocity, "StateRep::getUErr()");
        return data->uerr;
    }
    const Vector& getUDotErr() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Acceleration, "StateRep::getUDotErr()");
        return data->udoterr;
    }
    const Vector& getMultipliers() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Acceleration, "StateRep::getMultipliers()");
        return data->multipliers;
    }
    
    Vector& updYErr() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage(Stage::Velocity).prev(), "StateRep::updYErr()");
        return data->yerr;
    }
    Vector& updQErr() const{
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage(Stage::Position).prev(), "StateRep::updQErr()");
        return data->qerr;
    }
    Vector& updUErr() const{
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage(Stage::Velocity).prev(), "StateRep::updUErr()");
        return data->uerr;
    }
    Vector& updUDotErr() const{
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage(Stage::Acceleration).prev(), 
                            "StateRep::updUDotErr()");
        return data->udoterr;
    }
    Vector& updMultipliers() const{
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage(Stage::Acceleration).prev(), 
                            "StateRep::updMultipliers()");
        return data->multipliers;
    }
    
    const Vector& getEventTriggers() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Acceleration, "StateRep::getEventTriggers()");
        return data->allTriggers;
    }
    const Vector& getEventTriggersByStage(Stage g) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), g, "StateRep::getEventTriggersByStage()");
        return data->triggers[g];
    }
    
    // These are mutable; hence 'const'.
    Vector& updEventTriggers() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage(Stage::Acceleration).prev(), "StateRep::updEventTriggers()");
        return data->allTriggers;
    }
    Vector& updEventTriggersByStage(Stage g) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), g.prev(), "StateRep::updEventTriggersByStage()");
        return data->triggers[g];
    }
    
    // You can access a Model stage variable any time, but don't access others
    // until you have realized the Model stage.
    const AbstractValue& 
    getDiscreteVariable(SubsystemIndex subsys, int index) const {
        const PerSubsystemInfo& ss = data->subsystems[subsys];
    
        SimTK_INDEXCHECK(0,index,(int)ss.discreteInfo.size(),"StateRep::getDiscreteVariable()");
        const DiscreteVarInfo& dv = ss.discreteInfo[index];
    
        if (dv.getInvalidatedStage() > Stage::Model) {
            SimTK_STAGECHECK_GE(getSubsystemStage(subsys), 
                Stage::Model, "StateRep::getDiscreteVariable()");
        }
    
        return dv.getValue();
    }
    
    // You can update a Model stage variable from Topology stage, but higher variables 
    // must wait until you have realized the Model stage. This always backs the 
    // stage up to one earlier than the variable's stage.
    AbstractValue& 
    updDiscreteVariable(SubsystemIndex subsys, int index) {
        checkCanModify(subsys);
        PerSubsystemInfo& ss = data->subsystems[subsys];
    
        SimTK_INDEXCHECK(0,index,(int)ss.discreteInfo.size(),"StateRep::updDiscreteVariable()");
        DiscreteVarInfo& dv = ss.discreteInfo[index];
    
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), 
            std::min(dv.getInvalidatedStage().prev(), Stage(Stage::Model)), 
            "StateRep::updDiscreteVariable()");
    
        invalidateAll(dv.getInvalidatedStage());
    
        // We're now marking this variable as having been updated at the current time.
        return dv.updValue(data->t);
    }
    
    // Stage >= ce.stage
    const AbstractValue& 
    getCacheEntry(SubsystemIndex subsys, int index) const {
        const PerSubsystemInfo& ss = data->subsystems[subsys];
    
        SimTK_INDEXCHECK(0,index,(int)ss.cacheInfo.size(),"StateRep::getCacheEntry()");
        const CacheEntryInfo& ce = ss.cacheInfo[index];
    
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), 
            ce.getComputedByStage(), "StateRep::getCacheEntry()");
    
        return ce.getValue();
    }
    
    // Stage >= ce.stage-1; does not change stage
    AbstractValue& 
    updCacheEntry(SubsystemIndex subsys, int index) const {
        const PerSubsystemInfo& ss = data->subsystems[subsys];
    
        SimTK_INDEXCHECK(0,index,(int)ss.cacheInfo.size(),"StateRep::updCacheEntry()");
        CacheEntryInfo& ce = ss.cacheInfo[index];
    
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), 
            ce.getComputedByStage().prev(), "StateRep::updCacheEntry()");
    
        return ce.updValue();
    }

    bool isCacheValueValid(SubsystemIndex subx, CacheEntryIndex cx) const {
        const PerSubsystemInfo& ss = data->subsystems[subx];
        SimTK_INDEXCHECK(0,cx,(int)ss.cacheInfo.size(),"StateRep::isCacheValueValid()");
        const CacheEntryInfo& ce = ss.cacheInfo[cx];
        return ce.isValid(getSubsystemStage(subx), getSubsystemStageVersions(subx));
    }
    void markCacheValueValid(SubsystemIndex subx, CacheEntryIndex cx) const {
        const PerSubsystemInfo& ss = data->subsystems[subx];
        SimTK_INDEXCHECK(0,cx,(int)ss.cacheInfo.size(),"StateRep::markCacheValueValid()");
        CacheEntryInfo& ce = ss.cacheInfo[cx];
    
        // If this cache entry depends on anything, it can't be valid unless we're
        // at least *working* on its depends-on stage, meaning the current stage would
        // have to be the one before that. The depends-on stage is required to be at
        // least Stage::Topology, so its prev() stage exists.
        SimTK_STAGECHECK_GE(getSubsystemStage(subx), 
            ce.getDependsOnStage().prev(), "StateRep::markCacheValueValid()");

        ce.markAsValid(getSubsystemStageVersions(subx));
    }
    const Stage& getLowestStageModified() const {
        return data->lowestModifiedSystemStage;
    }
    void resetLowestStageModified() const {
        data->lowestModifiedSystemStage = std::min(Stage::Infinity, data->currentSystemStage.next());
    }
    
    const EnumerationSet<Stage>& getRestrictedStages() const {
        return restrictedStages;
    }

    const std::set<SubsystemIndex>& getRestrictedSubsystems() const {
        return restrictedSubsystems;
    }
    
    // Verify that a particular stage may be modified.
    void checkCanModify(Stage stage) const {
        SimTK_ASSERT1_ALWAYS(!restrictedStages.contains(stage),
                "Modification of state data for stage %s has been restricted.", stage.getName().c_str());
    }
    
    // Verify that a particular subsystem may be modified.
    void checkCanModify(SubsystemIndex subsystem) const {
        SimTK_ASSERT1_ALWAYS(restrictedSubsystems.find(subsystem) == restrictedSubsystems.end(),
                "Modification of state data for subsystem %d has been restricted.", (int) subsystem);
    }
    
    // Verify that this State permits all state variables to be modified.
    void checkCanModifyY() const {
        checkCanModify(Stage::Position);
        checkCanModify(Stage::Velocity);
        checkCanModify(Stage::Dynamics);
    }
    
    // Verify that this State permits unrestricted modifications.
    void checkCanModifyAnySubsystem() const {
        SimTK_ASSERT_ALWAYS(restrictedSubsystems.empty(),
                "Modification of state data has been restricted.");
    }
    
    String toString() const {
        String out;
        out += "<State>\n";
    
        out += "<Real name=time>" + String(data->t) + "</Real>\n";
    
        out += "<Vector name=q size=" + String(data->q.size()) + ">";
        if (data->q.size()) out += "\n";
        for (long i=0; i<data->q.size(); ++i)
            out += String(data->q[i]) + "\n";
        out += "</Vector>\n";
    
        out += "<Vector name=u size=" + String(data->u.size()) + ">";
        if (data->u.size()) out += "\n";
        for (long i=0; i<data->u.size(); ++i)
            out += String(data->u[i]) + "\n";
        out += "</Vector>\n";
    
        out += "<Vector name=z size=" + String(data->z.size()) + ">";
        if (data->z.size()) out += "\n";
        for (long i=0; i<data->z.size(); ++i)
            out += String(data->z[i]) + "\n";
        out += "</Vector>\n";
    
    
        for (SubsystemIndex ss(0); ss < (int)data->subsystems.size(); ++ss) {
            const PerSubsystemInfo& info = data->subsystems[ss];
            out += "<Subsystem index=" + String(ss) + " name=" + info.name 
                + " version=" + info.version + ">\n";
    
            out += "  <DISCRETE VARS TODO>\n";
    
            out += "  <Vector name=qInit size=" + String(info.qInit.size()) + ">\n";
            out += "  <Vector name=uInit size=" + String(info.uInit.size()) + ">\n";
            out += "  <Vector name=zInit size=" + String(info.zInit.size()) + ">\n";
    
            out += "  <Vector name=q size=" + String(info.q.size()) + ">\n";
            out += "  <Vector name=u size=" + String(info.u.size()) + ">\n";
            out += "  <Vector name=z size=" + String(info.z.size()) + ">\n";
    
            out += "</Subsystem>\n";
        }
    
        out += "</State>\n";
        return out;
    }
    
    String cacheToString() const {
        String out;
        out += "<Cache>\n";
        out += "<Stage>" + getSystemStage().getName() + "</Stage>\n";
    
        for (SubsystemIndex ss(0); ss < (int)data->subsystems.size(); ++ss) {
            const PerSubsystemInfo& info = data->subsystems[ss];
            out += "<Subsystem index=" + String(ss) + " name=" + info.name 
                + " version=" + info.version + ">\n";
            out += "  <Stage>" + info.currentStage.getName() + "</Stage>\n";
    
            out += "  <DISCRETE CACHE TODO>\n";
    
            out += "  <Vector name=qdot size=" + String(info.qdot.size()) + ">\n";
            out += "  <Vector name=udot size=" + String(info.udot.size()) + ">\n";
            out += "  <Vector name=zdot size=" + String(info.zdot.size()) + ">\n";
            out += "  <Vector name=qdotdot size=" + String(info.qdotdot.size()) + ">\n";
    
            out += "  <Vector name=qerr size=" + String(info.qerr.size()) + ">\n";
            out += "  <Vector name=uerr size=" + String(info.uerr.size()) + ">\n";
            out += "  <Vector name=udoterr size=" + String(info.udoterr.size()) + ">\n";
            out += "  <Vector name=multipliers size=" + String(info.multipliers.size()) + ">\n";
    
            for (int j=0; j<Stage::NValid; ++j) {
                out += "  <Vector name=triggers[";
                out += Stage::getValue(j).getName();
                out += "] size=" + String(info.triggers[j].size()) + ">\n";
            }
    
            out += "</Subsystem>\n";
        }
    
        out += "<Vector name=qdot size=" + String(data->qdot.size()) + ">";
        if (data->qdot.size()) out += "\n";
        for (long i=0; i<data->qdot.size(); ++i)
            out += String(data->qdot[i]) + "\n";
        out += "</Vector>\n";
    
        out += "<Vector name=udot size=" + String(data->udot.size()) + ">";
        if (data->udot.size()) out += "\n";
        for (long i=0; i<data->udot.size(); ++i)
            out += String(data->udot[i]) + "\n";
        out += "</Vector>\n";
    
        out += "<Vector name=zdot size=" + String(data->zdot.size()) + ">";
        if (data->zdot.size()) out += "\n";
        for (long i=0; i<data->zdot.size(); ++i)
            out += String(data->zdot[i]) + "\n";
        out += "</Vector>\n";
    
        out += "<Vector name=qdotdot size=" + String(data->qdotdot.size()) + ">";
        if (data->qdotdot.size()) out += "\n";
        for (long i=0; i<data->qdotdot.size(); ++i)
            out += String(data->qdotdot[i]) + "\n";
        out += "</Vector>\n";
    
        out += "<Vector name=qerr size=" + String(data->qerr.size()) + ">";
        if (data->qerr.size()) out += "\n";
        for (long i=0; i<data->qerr.size(); ++i)
            out += String(data->qerr[i]) + "\n";
        out += "</Vector>\n";
    
        out += "<Vector name=uerr size=" + String(data->uerr.size()) + ">";
        if (data->uerr.size()) out += "\n";
        for (long i=0; i<data->uerr.size(); ++i)
            out += String(data->uerr[i]) + "\n";
        out += "</Vector>\n";
    
        out += "<Vector name=udoterr size=" + String(data->udoterr.size()) + ">";
        if (data->udoterr.size()) out += "\n";
        for (long i=0; i<data->udoterr.size(); ++i)
            out += String(data->udoterr[i]) + "\n";
        out += "</Vector>\n";
    
        out += "<Vector name=multipliers size=" + String(data->multipliers.size()) + ">";
        if (data->multipliers.size()) out += "\n";
        for (long i=0; i<data->multipliers.size(); ++i)
            out += String(data->multipliers[i]) + "\n";
        out += "</Vector>\n";
    
        out += "</Cache>\n";
        return out;
    }
    StateData*                  data;
    EnumerationSet<Stage>       restrictedStages;
    std::set<SubsystemIndex>    restrictedSubsystems;
};

State::State() {
    rep = new StateRep();
}
State::~State() {
    delete rep;
}
State::State(const State& state) {
    rep = new StateRep(*state.rep);
}
void State::clear() {
    rep->clear();
}
void State::setNSubsystems(int i) {
    rep->setNSubsystems(i);
}
void State::initializeSubsystem(SubsystemIndex subsys, const String& name, const String& version) {
    rep->initializeSubsystem(subsys, name, version);
}
State& State::operator=(const State& state) {
    *rep = *state.rep;
    return *this;
}
SubsystemIndex State::addSubsystem(const String& name, const String& version) {
    return rep->addSubsystem(name, version);
}
int State::getNSubsystems() const {
    return rep->getNSubsystems();
}
const String& State::getSubsystemName(SubsystemIndex subsys) const {
    return rep->getSubsystemName(subsys);
}
const String& State::getSubsystemVersion(SubsystemIndex subsys) const {
    return rep->getSubsystemVersion(subsys);
}
const Stage& State::getSubsystemStage(SubsystemIndex subsys) const {
    return rep->getSubsystemStage(subsys);
}
const Stage& State::getSystemStage() const {
    return rep->getSystemStage();
}
void State::invalidateAll(Stage stage) const {
    rep->invalidateAll(stage);
}
void State::advanceSubsystemToStage(SubsystemIndex subsys, Stage stage) const {
    rep->advanceSubsystemToStage(subsys, stage);
}
void State::advanceSystemToStage(Stage stage) const {
    rep->advanceSystemToStage(stage);
}
QIndex State::allocateQ(SubsystemIndex subsys, const Vector& qInit) {
    return rep->allocateQ(subsys, qInit);
}
UIndex State::allocateU(SubsystemIndex subsys, const Vector& uInit) {
    return rep->allocateU(subsys, uInit);
}
ZIndex State::allocateZ(SubsystemIndex subsys, const Vector& zInit) {
    return rep->allocateZ(subsys, zInit);
}
QErrIndex State::allocateQErr(SubsystemIndex subsys, int nqerr) const {
    return rep->allocateQErr(subsys, nqerr);
}
UErrIndex State::allocateUErr(SubsystemIndex subsys, int nuerr) const {
    return rep->allocateUErr(subsys, nuerr);
}
UDotErrIndex State::allocateUDotErr(SubsystemIndex subsys, int nudoterr) const {
    return rep->allocateUDotErr(subsys, nudoterr);
}
EventTriggerByStageIndex State::allocateEventTrigger(SubsystemIndex subsys, Stage stage, int nevent) const {
    return rep->allocateEventTrigger(subsys, stage, nevent);
}
DiscreteVariableIndex State::allocateDiscreteVariable(SubsystemIndex subsys, Stage stage, AbstractValue* v) {
    return rep->allocateDiscreteVariable(subsys, stage, v);
}
CacheEntryIndex State::allocateCacheEntry(SubsystemIndex subsys, Stage stage, AbstractValue* v) const {
    return rep->allocateCacheEntry(subsys, stage, v);
}
int State::getNY() const {
    return rep->getNY();
}
SystemYIndex State::getQStart() const {
    return rep->getQStart();
}
int State::getNQ() const {
    return rep->getNQ();
}
SystemYIndex State::getUStart() const {
    return rep->getUStart();
}
int State::getNU() const {
    return rep->getNU();
}
SystemYIndex State::getZStart() const {
    return rep->getZStart();
}
int State::getNZ() const {
    return rep->getNZ();
}
int State::getNYErr() const {
    return rep->getNYErr();
}
SystemYErrIndex State::getQErrStart() const {
    return rep->getQErrStart();
}
int State::getNQErr() const {
    return rep->getNQErr();
}
SystemYErrIndex State::getUErrStart() const {
    return rep->getUErrStart();
}
int State::getNUErr() const {
    return rep->getNUErr();
}
int State::getNUDotErr() const {
    return rep->getNUDotErr();
}
int State::getNMultipliers() const {
    return getNUDotErr();
}
SystemQIndex State::getQStart(SubsystemIndex subsys) const {
    return rep->getQStart(subsys);
}
int State::getNQ(SubsystemIndex subsys) const {
    return rep->getNQ(subsys);
}
SystemUIndex State::getUStart(SubsystemIndex subsys) const {
    return rep->getUStart(subsys);
}
int State::getNU(SubsystemIndex subsys) const {
    return rep->getNU(subsys);
}
SystemZIndex State::getZStart(SubsystemIndex subsys) const {
    return rep->getZStart(subsys);
}
int State::getNZ(SubsystemIndex subsys) const {
    return rep->getNZ(subsys);
}
SystemQErrIndex State::getQErrStart(SubsystemIndex subsys) const {
    return rep->getQErrStart(subsys);
}
int State::getNQErr(SubsystemIndex subsys) const {
    return rep->getNQErr(subsys);
}
SystemUErrIndex State::getUErrStart(SubsystemIndex subsys) const {
    return rep->getUErrStart(subsys);
}
int State::getNUErr(SubsystemIndex subsys) const {
    return rep->getNUErr(subsys);
}
SystemUDotErrIndex State::getUDotErrStart(SubsystemIndex subsys) const {
    return rep->getUDotErrStart(subsys);
}
int State::getNUDotErr(SubsystemIndex subsys) const {
    return rep->getNUDotErr(subsys);
}
SystemMultiplierIndex State::getMultipliersStart(SubsystemIndex i) const {
    return SystemMultiplierIndex(getUDotErrStart(i));
}
int State::getNMultipliers(SubsystemIndex i) const {
    return getNUDotErr(i);
}
int State::getNEventTriggers() const {
    return rep->getNEventTriggers();
}
SystemEventTriggerIndex State::getEventTriggerStartByStage(Stage stage) const {
    return rep->getEventTriggerStartByStage(stage);
}
int State::getNEventTriggersByStage(Stage stage) const {
    return rep->getNEventTriggersByStage(stage);
}
SystemEventTriggerByStageIndex State::getEventTriggerStartByStage(SubsystemIndex subsys, Stage stage) const {
    return rep->getEventTriggerStartByStage(subsys, stage);
}
int State::getNEventTriggersByStage(SubsystemIndex subsys, Stage stage) const {
    return rep->getNEventTriggersByStage(subsys, stage);
}
const Vector& State::getEventTriggers() const {
    return rep->getEventTriggers();
}
const Vector& State::getEventTriggersByStage(Stage stage) const {
    return rep->getEventTriggersByStage(stage);
}
const Vector& State::getEventTriggersByStage(SubsystemIndex subsys, Stage stage) const {
    return rep->getEventTriggersByStage(subsys, stage);
}
Vector& State::updEventTriggers() const {
    return rep->updEventTriggers();
}
Vector& State::updEventTriggersByStage(Stage stage) const {
    return rep->updEventTriggersByStage(stage);
}
Vector& State::updEventTriggersByStage(SubsystemIndex subsys, Stage stage) const {
    return rep->updEventTriggersByStage(subsys, stage);
}
const Vector& State::getQ(SubsystemIndex subsys) const {
    return rep->getQ(subsys);
}
const Vector& State::getU(SubsystemIndex subsys) const {
    return rep->getU(subsys);
}
const Vector& State::getZ(SubsystemIndex subsys) const {
    return rep->getZ(subsys);
}
Vector& State::updQ(SubsystemIndex subsys) {
    return rep->updQ(subsys);
}
Vector& State::updU(SubsystemIndex subsys) {
    return rep->updU(subsys);
}
Vector& State::updZ(SubsystemIndex subsys) {
    return rep->updZ(subsys);
}
const Vector& State::getQDot(SubsystemIndex subsys) const {
    return rep->getQDot(subsys);
}
const Vector& State::getUDot(SubsystemIndex subsys) const {
    return rep->getUDot(subsys);
}
const Vector& State::getZDot(SubsystemIndex subsys) const {
    return rep->getZDot(subsys);
}
const Vector& State::getQDotDot(SubsystemIndex subsys) const {
    return rep->getQDotDot(subsys);
}
Vector& State::updQDot(SubsystemIndex subsys) const {
    return rep->updQDot(subsys);
}
Vector& State::updUDot(SubsystemIndex subsys) const {
    return rep->updUDot(subsys);
}
Vector& State::updZDot(SubsystemIndex subsys) const {
    return rep->updZDot(subsys);
}
Vector& State::updQDotDot(SubsystemIndex subsys) const {
    return rep->updQDotDot(subsys);
}
const Vector& State::getQErr(SubsystemIndex subsys) const {
    return rep->getQErr(subsys);
}
const Vector& State::getUErr(SubsystemIndex subsys) const {
    return rep->getUErr(subsys);
}
const Vector& State::getUDotErr(SubsystemIndex subsys) const {
    return rep->getUDotErr(subsys);
}
const Vector& State::getMultipliers(SubsystemIndex subsys) const {
    return rep->getMultipliers(subsys);
}
Vector& State::updQErr(SubsystemIndex subsys) const {
    return rep->updQErr(subsys);
}
Vector& State::updUErr(SubsystemIndex subsys) const {
    return rep->updUErr(subsys);
}
Vector& State::updUDotErr(SubsystemIndex subsys) const {
    return rep->updUDotErr(subsys);
}
Vector& State::updMultipliers(SubsystemIndex subsys) const {
    return rep->updMultipliers(subsys);
}
const Real& State::getTime() const {
    return rep->getTime();
}
const Vector& State::getY() const {
    return rep->getY();
}
const Vector& State::getQ() const {
    return rep->getQ();
}
const Vector& State::getU() const {
    return rep->getU();
}
const Vector& State::getZ() const {
    return rep->getZ();
}
Real& State::updTime() {
    return rep->updTime();
}
Vector& State::updY() {
    return rep->updY();
}
void State::setTime(Real t) {
    updTime() = t;
}
void State::setY(const Vector& y) {
    updY() = y;
}
Vector& State::updQ() {
    return rep->updQ();
}
Vector& State::updU() {
    return rep->updU();
}
Vector& State::updZ() {
    return rep->updZ();
}
void State::setQ(const Vector& q) {
    updQ() = q;
}
void State::setU(const Vector& u) {
    updU() = u;
}
void State::setZ(const Vector& z) {
    updZ() = z;
}
const Vector& State::getYDot() const {
    return rep->getYDot();
}
const Vector& State::getQDot() const {
    return rep->getQDot();
}
const Vector& State::getZDot() const {
    return rep->getZDot();
}
const Vector& State::getUDot() const {
    return rep->getUDot();
}
const Vector& State::getQDotDot() const {
    return rep->getQDotDot();
}
Vector& State::updYDot() const {
    return rep->updYDot();
}
Vector& State::updQDot() const {
    return rep->updQDot();
}
Vector& State::updZDot() const {
    return rep->updZDot();
}
Vector& State::updUDot() const {
    return rep->updUDot();
}
Vector& State::updQDotDot() const {
    return rep->updQDotDot();
}
const Vector& State::getYErr() const {
    return rep->getYErr();
}
const Vector& State::getQErr() const {
    return rep->getQErr();
}
const Vector& State::getUErr() const {
    return rep->getUErr();
}
const Vector& State::getUDotErr() const {
    return rep->getUDotErr();
}
const Vector& State::getMultipliers() const {
    return rep->getMultipliers();
}
Vector& State::updYErr() const {
    return rep->updYErr();
}
Vector& State::updQErr() const {
    return rep->updQErr();
}
Vector& State::updUErr() const {
    return rep->updUErr();
}
Vector& State::updUDotErr() const {
    return rep->updUDotErr();
}
Vector& State::updMultipliers() const {
    return rep->updMultipliers();
}
const AbstractValue& State::getDiscreteVariable(SubsystemIndex subsys, DiscreteVariableIndex index) const {
    return rep->getDiscreteVariable(subsys, index);
}
AbstractValue& State::updDiscreteVariable(SubsystemIndex subsys, DiscreteVariableIndex index) {
    return rep->updDiscreteVariable(subsys, index);
}
void State::setDiscreteVariable(SubsystemIndex i, DiscreteVariableIndex index, const AbstractValue& v) {
    updDiscreteVariable(i,index) = v;
}
const AbstractValue& State::getCacheEntry(SubsystemIndex subsys, CacheEntryIndex index) const {
    return rep->getCacheEntry(subsys, index);
}
AbstractValue& State::updCacheEntry(SubsystemIndex subsys, CacheEntryIndex index) const {
    return rep->updCacheEntry(subsys, index);
}


bool State::isCacheValueValid(SubsystemIndex subx, CacheEntryIndex cx) const {
    return rep->isCacheValueValid(subx, cx); 
}
void State::markCacheValueValid(SubsystemIndex subx, CacheEntryIndex cx) const {
    rep->markCacheValueValid(subx, cx); 
}
const Stage& State::getLowestStageModified() const {
    return rep->getLowestStageModified(); 
}
void State::resetLowestStageModified() const {
    rep->resetLowestStageModified(); 
}


void State::createRestrictedState
   (State&                   restrictedState,
    EnumerationSet<Stage>    restrictedStages, 
    std::set<SubsystemIndex> restrictedSubsystems)
{
    restrictedStages |= getRestrictedStages();
    const std::set<SubsystemIndex>& currentSubsystems = getRestrictedSubsystems();
    restrictedSubsystems.insert(currentSubsystems.begin(), currentSubsystems.end());
    delete restrictedState.rep;
    restrictedState.rep = new StateRep(*rep, restrictedStages, restrictedSubsystems);
}
const EnumerationSet<Stage>& State::getRestrictedStages() const {
    return rep->getRestrictedStages();
}
const std::set<SubsystemIndex>& State::getRestrictedSubsystems() const {
    return rep->getRestrictedSubsystems();
}
String State::toString() const {
    return rep->toString();
}
String State::cacheToString() const {
    return rep->cacheToString();
}

std::ostream& 
operator<<(std::ostream& o, const State& s) {
    o << "STATE:" << std::endl;
    o << s.toString() << std::endl;
    o << "CACHE:" << std::endl;
    return o << s.cacheToString() << std::endl;
}

} // namespace SimTK

