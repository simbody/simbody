/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-10 Stanford University and the Authors.        *
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
// resources that can be allocated at different stages, which we'll call an
// "allocation stack". A resource that was
// allocated at a later stage must be forgotten again when that stage is
// subsequently invalidated, and keeping their allocations stacked by
// stage allows that to be done efficiently.
//
// The method are templatized and expect the stacks to be in std::vectors
// of the same template. The template value must be a type that supports
// three methods (the template analog to virtual functions):
//      deepAssign()            a non-shallow assignment, i.e. clone the value
//      deepDestruct()          destroy any owned heap space
//      getAllocationStage()    return the stage being worked on when this was allocated
// The template type must otherwise support shallow copy semantics so that
// the Array_ can move them around without causing any heap activity.

// Clear the contents of an allocation stack, freeing up all associated heap space.
template <class T>
static void clearAllocationStack(Array_<T>& stack) {
    for (int i=stack.size()-1; i >= 0; --i)
        stack[i].deepDestruct();
    stack.clear();
}

// Resize the given allocation stack, taking care to free the heap space if the size is reduced.
template <class T>
static void resizeAllocationStack(Array_<T>& stack, int newSize) {
    assert(newSize >= 0);
    for (int i = stack.size()-1; i >= newSize; --i)
        stack[i].deepDestruct();
    stack.resize(newSize);
}

// Keep only those stack entries whose allocation stage is <= the supplied one.
template <class T>
static void popAllocationStackBackToStage(Array_<T>& stack, const Stage& g) {
    unsigned newSize = stack.size();
    while (newSize > 0 && stack[newSize-1].getAllocationStage() > g)
        stack[--newSize].deepDestruct();
    stack.resize(newSize); 
}

// Make this allocation stack the same as the source, copying only through the given stage.
template <class T>
static void copyAllocationStackThroughStage
   (Array_<T>& stack, const Array_<T>& src, const Stage& g) 
{
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
// allocation as the Array_ gets resized from time to time. However, the
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

    // These are the "virtual" methods required by template methods elsewhere.
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
// the number nz_integ of the z's that are defined by differential equations and
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

    DiscreteVarInfo(Stage allocation, Stage invalidated, AbstractValue* v)
    :   allocationStage(allocation), invalidatedStage(invalidated), value(v),
        autoUpdateEntry(), timeLastUpdated(NaN) 
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
    Real                 getTimeLastUpdated() const {assert(value); return timeLastUpdated;}
    AbstractValue&       updValue(Real updTime)
    {   assert(value); assert(updTime >= 0); 
        timeLastUpdated=updTime; return *value; }

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

    bool isCurrent(const Stage& current, const StageVersion versions[]) const 
    {   if (current >= computedByStage) return true;
        if (current <  dependsOnStage)  return false;
        return versions[dependsOnStage] == versionWhenLastComputed;}

    StageVersion getVersionWhenLastComputed() const {return versionWhenLastComputed;}

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
        clearDiscreteVars();
        clearCache();
        clearEvents();
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
        // destination is already initialized; copyFrom() will try
        // to reuse space and will properly clean up unused stuff
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
    void advanceToStage(Stage g) const {
        assert(g > Stage::Empty);
        assert(currentStage == g.prev());

        if (g == Stage::Topology) {
            qTopoSize = qInit.size();   // record now in case we have to
            uTopoSize = uInit.size();   //  back up later
            zTopoSize = zInit.size();
        }

        // This validates whatever the current version number is of Stage g.
        currentStage = g;
    }


    void clearReferencesToModelStageGlobals() {
        qstart.invalidate(); ustart.invalidate(); 
        zstart.invalidate();
        q.clear(); u.clear(); z.clear();
        qdot.clear(); udot.clear(); zdot.clear(); qdotdot.clear();
    }

    void clearReferencesToInstanceStageGlobals() const {
        // These are all mutable cache entries.
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

    // Topology and Model stage definitions. TODO: some of these aren't used yet
    Array_<MechanicalStateInfo>    quInfo;
    Array_<AuxiliaryStateInfo>     zInfo;
    Array_<DiscreteVarInfo>        discreteInfo;

    // Topology, Model, and Instance stage definitions.
    //mutable Array_<ConstraintInfo> constraintInfo;
    mutable Array_<EventInfo>      eventInfo;
    mutable Array_<CacheEntryInfo> cacheInfo;

        // AGGREGATE GLOBAL RESOURCE NEEDS //
    
    // These accumulate default values for this subsystem's use of shared
    // global state variables. After the System is advanced to Stage::Model,
    // the state will allocate those globals and copy these initial
    // values into them. The lengths of these Vectors define the 
    // needs of this Subsystem. Some of these are allocated at Topology
    // stage, and some at Model stage so we record the size upon realizeTopology()
    // in case we need to wipe out just Model stage later.
    Vector qInit, uInit, zInit;
    mutable int qTopoSize, uTopoSize, zTopoSize;

    // For constraints we need just lengths (nmultipliers==nudoterr).
    mutable int nqerr, nuerr, nudoterr;

    // For event trigger functions we need just lengths, for each stage.
    mutable int ntriggers[Stage::NValid];

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
    mutable SystemQErrIndex                 qerrstart;
    mutable SystemUErrIndex                 uerrstart;
    mutable SystemUDotErrIndex              udoterrstart;
    mutable SystemEventTriggerByStageIndex  triggerstart[Stage::NValid];
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
        clearDiscreteVars(); clearCache(); clearEvents();
        qInit.clear(); uInit.clear(); zInit.clear();
        qTopoSize = uTopoSize = zTopoSize = 0;
        nqerr = nuerr = nudoterr= 0;
        qstart.invalidate(); ustart.invalidate(); zstart.invalidate();
        qerrstart.invalidate(); uerrstart.invalidate(); udoterrstart.invalidate();
        for (int j=0; j<Stage::NValid; ++j) {
            ntriggers[j] = 0;
            triggerstart[j].invalidate();
            stageVersions[j] = 1; // never 0
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

    void copyDiscreteVarsThroughStage
       (const Array_<DiscreteVarInfo>& src, const Stage& g)
    {   copyAllocationStackThroughStage(discreteInfo, src, g); }
    void copyCacheThroughStage(const Array_<CacheEntryInfo>& src, const Stage& g)
    {   copyAllocationStackThroughStage(cacheInfo, src, g); }
    void copyEventsThroughStage(const Array_<EventInfo>& src, const Stage& g)
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

            // Some of these may have been allocated at Topology stage so
            // we shrink them back just to that size.
            qInit.resizeKeep(qTopoSize); 
            uInit.resizeKeep(uTopoSize); 
            zInit.resizeKeep(zTopoSize);
        }

        if (g == Stage::Empty) {
            // Throw out everything, reset stage versions to 1. Leave
            // name and version alone.
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

        if (targetStage >= Stage::Topology) {
            qTopoSize=src.qTopoSize; uTopoSize=src.uTopoSize; zTopoSize=src.zTopoSize;
            if (targetStage == Stage::Topology) {
                qInit = src.qInit(0, qTopoSize); 
                uInit = src.uInit(0, uTopoSize); 
                zInit = src.zInit(0, zTopoSize); 
            } else { // targetStage is at least Model
                qInit = src.qInit;
                uInit = src.uInit;
                zInit = src.zInit;
            }

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
        // be valid if they were valid in the source and depended only on
        // things we copied.
        for (int i=0; i<=targetStage; ++i)
            stageVersions[i] = src.stageVersions[i];
        // The rest of the stages need to be invalidated in the destination
        // since we didn't copy any state information from those stages.
        for (int i=targetStage+1; i<=src.currentStage; ++i)
            stageVersions[i] = src.stageVersions[i] + 1;

        // Subsystem stage should now match what we copied.
        currentStage = targetStage;
    }
};


//==============================================================================
//                                 STATE IMPL
//==============================================================================

class StateImpl {
public:
    StateImpl() 
    :   t(NaN), currentSystemStage(Stage::Empty), 
        lowestModifiedSystemStage(Stage::Topology) {} 

    // We'll do the copy constructor and assignment explicitly here
    // to get tight control over what's allowed.
    StateImpl(const StateImpl& src)
    :   currentSystemStage(Stage::Empty), lowestModifiedSystemStage(Stage::Topology)
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

    StateImpl& operator=(const StateImpl& src) {
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
        return *this;
    }

    ~StateImpl() {   // default destructor
    }

    // Copies all the variables but not the cache.
    StateImpl* clone() const {return new StateImpl(*this);}

    const Stage& getSystemStage() const {return currentSystemStage;}
    Stage&       updSystemStage() const {return currentSystemStage;} // mutable


    const PerSubsystemInfo& getSubsystem(int subsystem) const {
        SimTK_INDEXCHECK(subsystem, (int)subsystems.size(), "StateImpl::getSubsystem()");
        return subsystems[subsystem];
    }

    PerSubsystemInfo& updSubsystem(int subsystem) {
        SimTK_INDEXCHECK(subsystem, (int)subsystems.size(), "StateImpl::updSubsystem()");
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
            q.clear(); u.clear(); z.clear(); // y views
            // Finally nuke the actual y data.
            y.unlockShape(); y.clear(); 

            qdot.clear(); udot.clear(); zdot.clear();   // ydot views
            ydot.unlockShape();        ydot.clear();    // ydot data
            qdotdot.unlockShape();     qdotdot.clear(); // qdotdot data (no views)
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
    void advanceSystemToStage(Stage g) const {
        assert(g > Stage::Empty);
        assert(currentSystemStage == g.prev());
        assert(allSubsystemsAtLeastAtStage(g));

        if (g == Stage::Topology) {
            // As the final "Topology" step, initialize time to 0 (it's NaN before this).
            const_cast<StateImpl*>(this)->t = 0;
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
            // We need write access temporarily to set up the state.
            StateImpl* wThis = const_cast<StateImpl*>(this);
            wThis->y.resize(nq+nu+nz);      wThis->y.lockShape();

            ydot.resize(nq+nu+nz);          ydot.lockShape();
            qdotdot.resize(nq);             qdotdot.lockShape();

            // Allocate subviews of the shared state & cache entries.
            wThis->q.viewAssign(wThis->y(0,nq));
            wThis->u.viewAssign(wThis->y(nq,nu));
            wThis->z.viewAssign(wThis->y(nq+nu,nz));

            qdot.viewAssign(ydot(0,     nq));
            udot.viewAssign(ydot(nq,    nu));
            zdot.viewAssign(ydot(nq+nu, nz));

            // Now partition the global resources among the subsystems and copy
            // in the initial values for the state variables.
            SystemQIndex nxtq(0);
            SystemUIndex nxtu(0);
            SystemZIndex nxtz(0);

            for (SubsystemIndex i(0); i<(int)subsystems.size(); ++i) {
                PerSubsystemInfo& ss = const_cast<PerSubsystemInfo&>(subsystems[i]);
                const int nq=ss.qInit.size(), nu=ss.uInit.size(), nz=ss.zInit.size();

                // Assign the starting indices.
                ss.qstart=nxtq; ss.ustart=nxtu; ss.zstart=nxtz;
 
                // Build the views.
                ss.q.viewAssign(wThis->q(nxtq, nq)); 
                ss.q = ss.qInit;
                ss.u.viewAssign(wThis->u(nxtu, nu)); 
                ss.u = ss.uInit;
                ss.z.viewAssign(wThis->z(nxtz, nz)); 
                ss.z = ss.zInit;

                ss.qdot.viewAssign(qdot(nxtq, nq));
                ss.qdotdot.viewAssign(qdotdot(nxtq, nq));
                ss.udot.viewAssign(udot(nxtu, nu));
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
                const PerSubsystemInfo& ss = subsystems[i];

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

    
    void setNumSubsystems(int i) {
        assert(i >= 0);
        subsystems.clear();
        subsystems.resize(i);
    }
    
    void initializeSubsystem(SubsystemIndex i, const String& name, const String& version) {
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
    
    // Move the stage for a particular subsystem from g-1 to g. No other subsystems
    // are affected, nor the global system stage.
    void advanceSubsystemToStage(SubsystemIndex subsys, Stage g) const {
        subsystems[subsys].advanceToStage(g);
        // We don't automatically advance the System stage even if this brings
        // ALL the subsystems up to stage g.
    }
    

    
    // We don't expect State entry allocations to be performance critical so
    // we'll keep error checking on even in Release mode.
    
    QIndex allocateQ(SubsystemIndex subsys, const Vector& qInit) {
        SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), Stage::Model, "StateImpl::allocateQ()");
        const QIndex nxt(subsystems[subsys].qInit.size());
        subsystems[subsys].qInit.resizeKeep(nxt + qInit.size());
        subsystems[subsys].qInit(nxt, qInit.size()) = qInit;
        return nxt;
    }
    
    UIndex allocateU(SubsystemIndex subsys, const Vector& uInit) {
        SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), Stage::Model, "StateImpl::allocateU()");
        const UIndex nxt(subsystems[subsys].uInit.size());
        subsystems[subsys].uInit.resizeKeep(nxt + uInit.size());
        subsystems[subsys].uInit(nxt, uInit.size()) = uInit;
        return nxt;
    }
    ZIndex allocateZ(SubsystemIndex subsys, const Vector& zInit) {
        SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), Stage::Model, "StateImpl::allocateZ()");
        const ZIndex nxt(subsystems[subsys].zInit.size());
        subsystems[subsys].zInit.resizeKeep(nxt + zInit.size());
        subsystems[subsys].zInit(nxt, zInit.size()) = zInit;
        return nxt;
    }
    
    QErrIndex allocateQErr(SubsystemIndex subsys, int nqerr) const {
        SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), Stage::Instance, "StateImpl::allocateQErr()");
        const QErrIndex nxt(subsystems[subsys].nqerr);
        subsystems[subsys].nqerr += nqerr;
        return nxt;
    }
    UErrIndex allocateUErr(SubsystemIndex subsys, int nuerr) const {
        SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), Stage::Instance, "StateImpl::al()");
        const UErrIndex nxt(subsystems[subsys].nuerr);
        subsystems[subsys].nuerr += nuerr;
        return nxt;
    }
    UDotErrIndex allocateUDotErr(SubsystemIndex subsys, int nudoterr) const {
        SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), Stage::Instance, "StateImpl::allocateUDotErr()");
        const UDotErrIndex nxt(subsystems[subsys].nudoterr);
        subsystems[subsys].nudoterr += nudoterr;
        return nxt;
    }
    EventTriggerByStageIndex allocateEventTrigger(SubsystemIndex subsys, Stage g, int nt) const {
        SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), Stage::Instance, "StateImpl::allocateEvent()");
        const EventTriggerByStageIndex nxt(subsystems[subsys].ntriggers[g]);
        subsystems[subsys].ntriggers[g] += nt;
        return nxt;
    }
    
    // Topology- and Model-stage State variables can only be added during construction; that is,
    // while stage <= Topology. Other entries can be added while stage < Model.
    DiscreteVariableIndex allocateDiscreteVariable(SubsystemIndex subsys, Stage invalidates, AbstractValue* vp) {
        SimTK_STAGECHECK_RANGE_ALWAYS(Stage(Stage::LowestRuntime).prev(), invalidates, Stage::HighestRuntime, 
            "StateImpl::allocateDiscreteVariable()");
    
        const Stage maxAcceptable = (invalidates <= Stage::Model ? Stage::Empty : Stage::Topology);
        SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), 
            maxAcceptable.next(), "StateImpl::allocateDiscreteVariable()");
    
        PerSubsystemInfo& ss = subsystems[subsys];
        const DiscreteVariableIndex nxt(ss.discreteInfo.size());
        ss.discreteInfo.push_back(DiscreteVarInfo(getSubsystemStage(subsys).next(),invalidates,vp));
        return nxt;
    }
    
    // Cache entries can be allocated while stage < Instance.
    CacheEntryIndex allocateCacheEntry(SubsystemIndex subsys, Stage dependsOn, Stage computedBy, AbstractValue* vp) const {
        SimTK_STAGECHECK_RANGE_ALWAYS(Stage(Stage::LowestRuntime).prev(), dependsOn, Stage::HighestRuntime, 
            "StateImpl::allocateCacheEntry()");
        SimTK_STAGECHECK_RANGE_ALWAYS(dependsOn, computedBy, Stage::Infinity, 
            "StateImpl::allocateCacheEntry()");
        SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), 
            Stage::Instance, "StateImpl::allocateCacheEntry()");

        const PerSubsystemInfo& ss = subsystems[subsys];
        const CacheEntryIndex nxt(ss.cacheInfo.size());
        ss.cacheInfo.push_back(CacheEntryInfo(getSubsystemStage(subsys).next(),dependsOn,computedBy,vp)); // mutable
        return nxt;
    }

    // Allocate a discrete variable and a corresponding cache entry for
    // updating it, and connect them together.
    DiscreteVariableIndex
    allocateAutoUpdateDiscreteVariable(SubsystemIndex subsys, Stage invalidates, AbstractValue* vp,
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
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::getNY()");
        return y.size();
    }
    
    SystemYIndex getQStart() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::getQStart()");
        return SystemYIndex(0); // q's come first
    }
    int getNQ() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::getNQ()");
        return q.size();
    }
    
    SystemYIndex getUStart() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::getUStart()");
        return SystemYIndex(q.size()); // u's come right after q's
    }
    int getNU() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::getNU()");
        return u.size();
    }
    
    SystemYIndex getZStart() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::getZStart()");
        return SystemYIndex(q.size() + u.size()); // q,u, then z
    }
    int getNZ() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::getNZ()");
        return z.size();
    }
    
    int getNYErr() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, "StateImpl::getNYErr()");
        return yerr.size();
    }
    
    SystemYErrIndex getQErrStart() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, "StateImpl::getQErrStart()");
        return SystemYErrIndex(0); // qerr's come first
    }
    int getNQErr() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, "StateImpl::getNQErr()");
        return qerr.size();
    }
    
    SystemYErrIndex getUErrStart() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, "StateImpl::getUErrStart()");
        return SystemYErrIndex(qerr.size()); // uerr's follow qerrs
    }
    int getNUErr() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, "StateImpl::getNUErr()");
        return uerr.size();
    }
    
    // UDot errors are independent of qerr & uerr.
    // This is used for multipliers also.
    int getNUDotErr() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, "StateImpl::getNUDotErr()");
        return udoterr.size();
    }
    
    int getNEventTriggers() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, "StateImpl::getNEventTriggers()");
        return allTriggers.size();
    }
    
    SystemEventTriggerIndex getEventTriggerStartByStage(Stage g) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, "StateImpl::getEventTriggerStartByStage()");
        int nxt = 0;
        for (int j=0; j<g; ++j)
            nxt += triggers[j].size();
        return SystemEventTriggerIndex(nxt); // g starts where g-1 leaves off
    }
    
    int getNEventTriggersByStage(Stage g) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, "StateImpl::getNEventTriggersByStage()");
        return triggers[g].size();
    }
    
        // Subsystem dimensions.
    
    SystemQIndex getQStart(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::getQStart(subsys)");
        return getSubsystem(subsys).qstart;
    }
    int getNQ(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::getNQ(subsys)");
        return getSubsystem(subsys).q.size();
    }
    
    SystemUIndex getUStart(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::getUStart(subsys)");
        return getSubsystem(subsys).ustart;
    }
    int getNU(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::getNU(subsys)");
        return getSubsystem(subsys).u.size();
    }
    
    SystemZIndex getZStart(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::getZStart(subsys)");
        return getSubsystem(subsys).zstart;
    }
    int getNZ(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::getNZ(subsys)");
        return getSubsystem(subsys).z.size();
    }
    
    SystemQErrIndex getQErrStart(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, "StateImpl::getQErrStart(subsys)");
        return getSubsystem(subsys).qerrstart;
    }
    int getNQErr(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, "StateImpl::getNQErr(subsys)");
        return getSubsystem(subsys).qerr.size();
    }
    
    SystemUErrIndex getUErrStart(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, "StateImpl::getUErrStart(subsys)");
        return getSubsystem(subsys).uerrstart;
    }
    int getNUErr(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, "StateImpl::getNUErr(subsys)");
        return getSubsystem(subsys).uerr.size();
    }
    
    // These are used for multipliers also.
    SystemUDotErrIndex getUDotErrStart(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, "StateImpl::getUDotErrStart(subsys)");
        return getSubsystem(subsys).udoterrstart;
    }
    int getNUDotErr(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, "StateImpl::getNUDotErr(subsys)");
        return getSubsystem(subsys).udoterr.size();
    }
    
    SystemEventTriggerByStageIndex getEventTriggerStartByStage(SubsystemIndex subsys, Stage g) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, "StateImpl::getEventTriggerStartByStage(subsys)");
        return getSubsystem(subsys).triggerstart[g];
    }
    
    int getNEventTriggersByStage(SubsystemIndex subsys, Stage g) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, "StateImpl::getNEventTriggersByStage(subsys)");
        return getSubsystem(subsys).triggers[g].size();
    }
    
        // Per-subsystem access to the global shared variables.
    
    const Vector& getQ(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::getQ(subsys)");
        return getSubsystem(subsys).q;
    }
    const Vector& getU(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::getU(subsys)");
        return getSubsystem(subsys).u;
    }
    const Vector& getZ(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::getZ(subsys)");
        return getSubsystem(subsys).z;
    }
    
    const Vector& getQDot(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::getQDot(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Velocity, "StateImpl::getQDot(subsys)");
        return getSubsystem(subsys).qdot;
    }
    const Vector& getUDot(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::getUDot(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Acceleration, "StateImpl::getUDot(subsys)");
        return getSubsystem(subsys).udot;
    }
    const Vector& getZDot(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::getZDot(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Dynamics, "StateImpl::getZDot(subsys)");
        return getSubsystem(subsys).zdot;
    }
    const Vector& getQDotDot(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::getQDotDot(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Acceleration, "StateImpl::getQDotDot(subsys)");
        return getSubsystem(subsys).qdotdot;
    }
    
    Vector& updQ(SubsystemIndex subsys) {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::updQ(subsys)");
        invalidateAll(Stage::Position);
        return updSubsystem(subsys).q;
    }
    Vector& updU(SubsystemIndex subsys) {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::updU(subsys)");
        invalidateAll(Stage::Velocity);
        return updSubsystem(subsys).u;
    }
    Vector& updZ(SubsystemIndex subsys) {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::updZ(subsys)");
        invalidateAll(Stage::Dynamics);
        return updSubsystem(subsys).z;
    }
    
        // These are mutable so the routines are const.
    
    Vector& updQDot(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::updQDot(subsys)");
        return getSubsystem(subsys).qdot;
    }
    Vector& updUDot(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::updUDot(subsys)");
        return getSubsystem(subsys).udot;
    }
    Vector& updZDot(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::updZDot(subsys)");
        return getSubsystem(subsys).zdot;
    }
    Vector& updQDotDot(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::updQDotDot(subsys)");
        return getSubsystem(subsys).qdotdot;
    }
    
    
    const Vector& getQErr(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, "StateImpl::getQErr(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Position, "StateImpl::getQErr(subsys)");
        return getSubsystem(subsys).qerr;
    }
    const Vector& getUErr(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, "StateImpl::getUErr(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Velocity, "StateImpl::getUErr(subsys)");
        return getSubsystem(subsys).uerr;
    }
    const Vector& getUDotErr(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, "StateImpl::getUDotErr(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Acceleration, "StateImpl::getUDotErr(subsys)");
        return getSubsystem(subsys).udoterr;
    }
    const Vector& getMultipliers(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, "StateImpl::getMultipliers(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Acceleration, "StateImpl::getMultipliers(subsys)");
        return getSubsystem(subsys).multipliers;
    }
    
    const Vector& getEventTriggersByStage(SubsystemIndex subsys, Stage g) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, "StateImpl::getEventTriggersByStage(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), g, "StateImpl::getEventTriggersByStage(subsys)");
        return getSubsystem(subsys).triggers[g];
    }
    
    Vector& updQErr(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, "StateImpl::updQErr(subsys)");
        return getSubsystem(subsys).qerr;
    }
    Vector& updUErr(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, "StateImpl::updUErr(subsys)");
        return getSubsystem(subsys).uerr;
    }
    Vector& updUDotErr(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, "StateImpl::updUDotErr(subsys)");
        return getSubsystem(subsys).udoterr;
    }
    Vector& updMultipliers(SubsystemIndex subsys) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, "StateImpl::updMultipliers(subsys)");
        return getSubsystem(subsys).multipliers;
    }
    Vector& updEventTriggersByStage(SubsystemIndex subsys, Stage g) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, "StateImpl::updEventTriggersByStage(subsys)");
        return getSubsystem(subsys).triggers[g];
    }
    
        // Direct access to the global shared state and cache entries.
        // Time is allocated in Stage::Topology, State in Stage::Model, and
        // Cache in Stage::Instance.
    
    const Real& getTime() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Topology, "StateImpl::getTime()");
        return t;
    }
    
    const Vector& getY() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::getY()");
        return y;
    }
    
    const Vector& getQ() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::getQ()");
        return q;
    }
    
    const Vector& getU() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::getU()");
        return u;
    }
    
    const Vector& getZ() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::getZ()");
        return z;
    }
    
    
    // You can call these as long as stage >= allocation stage, but the
    // stage will be backed up if necessary to one stage prior to the invalidated stage.
    Real& updTime() {  // Back to Stage::Time-1
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Topology, "StateImpl::updTime()");
        invalidateAll(Stage::Time);
        return t;
    }
    
    Vector& updY() {    // Back to Stage::Position-1
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::updY()");
        invalidateAll(Stage::Position);
        return y;
    }
    
    Vector& updQ() {    // Stage::Position-1
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::updQ()");
        invalidateAll(Stage::Position);
        return q;
    }
    
    Vector& updU() {     // Stage::Velocity-1
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::updU()");
        invalidateAll(Stage::Velocity);
        return u;
    }
    
    Vector& updZ() {     // Stage::Dynamics-1
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::updZ()");
        invalidateAll(Stage::Dynamics);
        return z;
    }
    
    // These cache entries you can get at their "computedBy" stages.
    const Vector& getYDot() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Acceleration, "StateImpl::getYDot()");
        return ydot;
    }
    
    const Vector& getQDot() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Velocity, "StateImpl::getQDot()");
        return qdot;
    }
    
    const Vector& getZDot() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Dynamics, "StateImpl::getZDot()");
        return zdot;
    }
    
    const Vector& getUDot() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Acceleration, "StateImpl::getUDot()");
        return udot;
    }
    
    const Vector& getQDotDot() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Acceleration, "StateImpl::getQDotDot()");
        return qdotdot;
    }
    
    // Cache updates are allowed any time after they have been allocated.
    Vector& updYDot() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::updYDot()");
        return ydot;
    }
    
    Vector& updQDot() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::updQDot()");
        return qdot;
    }
    
    Vector& updUDot() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::updUDot()");
        return udot;
    }
    
    Vector& updZDot() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::updZDot()");
        return zdot;
    }
    
    Vector& updQDotDot() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateImpl::updQDotDot()");
        return qdotdot;
    }
    
    
    const Vector& getYErr() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Velocity, "StateImpl::getYErr()");
        return yerr;
    }
    
    const Vector& getQErr() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Position, "StateImpl::getQErr()");
        return qerr;
    }
    const Vector& getUErr() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Velocity, "StateImpl::getUErr()");
        return uerr;
    }
    const Vector& getUDotErr() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Acceleration, "StateImpl::getUDotErr()");
        return udoterr;
    }
    const Vector& getMultipliers() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Acceleration, "StateImpl::getMultipliers()");
        return multipliers;
    }
    
    Vector& updYErr() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, "StateImpl::updYErr()");
        return yerr;
    }
    Vector& updQErr() const{
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, "StateImpl::updQErr()");
        return qerr;
    }
    Vector& updUErr() const{
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, "StateImpl::updUErr()");
        return uerr;
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
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Acceleration, "StateImpl::getEventTriggers()");
        return allTriggers;
    }
    const Vector& getEventTriggersByStage(Stage g) const {
        SimTK_STAGECHECK_GE(getSystemStage(), g, "StateImpl::getEventTriggersByStage()");
        return triggers[g];
    }
    
    // These are mutable; hence 'const'.
    Vector& updEventTriggers() const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, "StateImpl::updEventTriggers()");
        return allTriggers;
    }
    Vector& updEventTriggersByStage(Stage g) const {
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Instance, "StateImpl::updEventTriggersByStage()");
        return triggers[g];
    }

    CacheEntryIndex getDiscreteVarUpdateIndex(SubsystemIndex subsys, DiscreteVariableIndex index) const {
        const PerSubsystemInfo& ss = subsystems[subsys];
        SimTK_INDEXCHECK(index,(int)ss.discreteInfo.size(),
            "StateImpl::getDiscreteVarUpdateIndex()");
        const DiscreteVarInfo& dv = ss.discreteInfo[index];
        return dv.getAutoUpdateEntry();
    } 

    Stage getDiscreteVarAllocationStage(SubsystemIndex subsys, DiscreteVariableIndex index) const {
        const PerSubsystemInfo& ss = subsystems[subsys];
        SimTK_INDEXCHECK(index,(int)ss.discreteInfo.size(),
            "StateImpl::getDiscreteVarAllocationStage()");
        const DiscreteVarInfo& dv = ss.discreteInfo[index];
        return dv.getAllocationStage();
    } 

    Stage getDiscreteVarInvalidatesStage(SubsystemIndex subsys, DiscreteVariableIndex index) const {
        const PerSubsystemInfo& ss = subsystems[subsys];
        SimTK_INDEXCHECK(index,(int)ss.discreteInfo.size(),
            "StateImpl::getDiscreteVarInvalidatesStage()");
        const DiscreteVarInfo& dv = ss.discreteInfo[index];
        return dv.getInvalidatedStage();
    } 

    // You can access at any time a variable that was allocated during realizeTopology(), 
    // but don't access others until you have done realizeModel().
    const AbstractValue& 
    getDiscreteVariable(SubsystemIndex subsys, DiscreteVariableIndex index) const {
        const PerSubsystemInfo& ss = subsystems[subsys];
    
        SimTK_INDEXCHECK(index,(int)ss.discreteInfo.size(),"StateImpl::getDiscreteVariable()");
        const DiscreteVarInfo& dv = ss.discreteInfo[index];
    
        if (dv.getAllocationStage() > Stage::Topology) {
            SimTK_STAGECHECK_GE(getSubsystemStage(subsys), 
                Stage::Model, "StateImpl::getDiscreteVariable()");
        }
    
        return dv.getValue();
    }
    Real getDiscreteVarLastUpdateTime(SubsystemIndex subsys, DiscreteVariableIndex index) const {
        const PerSubsystemInfo& ss = subsystems[subsys];
        SimTK_INDEXCHECK(index,(int)ss.discreteInfo.size(),"StateImpl::getDiscreteVarLastUpdateTime()");
        const DiscreteVarInfo& dv = ss.discreteInfo[index];
        return dv.getTimeLastUpdated();
    }

    const AbstractValue& getDiscreteVarUpdateValue(SubsystemIndex subsys, DiscreteVariableIndex index) const {
        const CacheEntryIndex cx = getDiscreteVarUpdateIndex(subsys,index);
        SimTK_ERRCHK2(cx.isValid(), "StateImpl::getDiscreteVarUpdateValue()", 
            "Subsystem %d has a discrete variable %d but it does not have an"
            " associated update cache variable.", (int)subsys, (int)index);
        return getCacheEntry(subsys, cx);
    }
    AbstractValue& updDiscreteVarUpdateValue(SubsystemIndex subsys, DiscreteVariableIndex index) const {
        const CacheEntryIndex cx = getDiscreteVarUpdateIndex(subsys,index);
        SimTK_ERRCHK2(cx.isValid(), "StateImpl::updDiscreteVarUpdateValue()", 
            "Subsystem %d has a discrete variable %d but it does not have an"
            " associated update cache variable.", (int)subsys, (int)index);
        return updCacheEntry(subsys, cx);
    }
    bool isDiscreteVarUpdateValueRealized(SubsystemIndex subsys, DiscreteVariableIndex index) const {
        const CacheEntryIndex cx = getDiscreteVarUpdateIndex(subsys,index);
        SimTK_ERRCHK2(cx.isValid(), "StateImpl::isDiscreteVarUpdateValueRealized()", 
            "Subsystem %d has a discrete variable %d but it does not have an"
            " associated update cache variable.", (int)subsys, (int)index);
        return isCacheValueRealized(subsys, cx);
    }
    void markDiscreteVarUpdateValueRealized(SubsystemIndex subsys, DiscreteVariableIndex index) const {
        const CacheEntryIndex cx = getDiscreteVarUpdateIndex(subsys,index);
        SimTK_ERRCHK2(cx.isValid(), "StateImpl::markDiscreteVarUpdateValueRealized()", 
            "Subsystem %d has a discrete variable %d but it does not have an"
            " associated update cache variable.", (int)subsys, (int)index);
        markCacheValueRealized(subsys, cx);
    }

    // You can update at any time a variable that was allocated during realizeTopology(), 
    // but later variables must wait until you have done realizeModel(). This always 
    // backs the stage up to one earlier than the variable's stage.
    AbstractValue& 
    updDiscreteVariable(SubsystemIndex subsys, DiscreteVariableIndex index) {
        PerSubsystemInfo& ss = subsystems[subsys];
    
        SimTK_INDEXCHECK(index,(int)ss.discreteInfo.size(),"StateImpl::updDiscreteVariable()");
        DiscreteVarInfo& dv = ss.discreteInfo[index];

        if (dv.getAllocationStage() > Stage::Topology) {
            SimTK_STAGECHECK_GE(getSubsystemStage(subsys), 
                Stage::Model, "StateImpl::updDiscreteVariable()");
        }
    
        invalidateAll(dv.getInvalidatedStage());
    
        // We're now marking this variable as having been updated at the current time.
        return dv.updValue(t);
    }
    
    Stage getCacheEntryAllocationStage(SubsystemIndex subsys, CacheEntryIndex index) const {
        const PerSubsystemInfo& ss = subsystems[subsys];
        SimTK_INDEXCHECK(index,(int)ss.cacheInfo.size(),
            "StateImpl::getCacheEntryAllocationStage()");
        const CacheEntryInfo& ce = ss.cacheInfo[index];
        return ce.getAllocationStage();
    } 

    // Stage >= ce.stage
    const AbstractValue& 
    getCacheEntry(SubsystemIndex subsys, CacheEntryIndex index) const {
        const PerSubsystemInfo& ss = subsystems[subsys];
    
        SimTK_INDEXCHECK(index,(int)ss.cacheInfo.size(),"StateImpl::getCacheEntry()");
        const CacheEntryInfo& ce = ss.cacheInfo[index];
    
        const Stage stageNow = getSubsystemStage(subsys);
        SimTK_STAGECHECK_GE_ALWAYS(stageNow, 
            ce.getDependsOnStage(), "StateImpl::getCacheEntry()");

        if (stageNow < ce.getComputedByStage()) {
            const StageVersion currentDependsOnVersion = 
                getSubsystemStageVersions(subsys)[ce.getDependsOnStage()];
            const StageVersion lastCacheVersion = 
                ce.getVersionWhenLastComputed();

            if (lastCacheVersion != currentDependsOnVersion) {
                SimTK_THROW4(Exception::CacheEntryOutOfDate,
                    stageNow, ce.getDependsOnStage(), currentDependsOnVersion, lastCacheVersion);
            }
        }

        // If we get here then we're either past the "computed by" stage, or we're
        // past "depends on" with an explicit validation having been made at the
        // "depends on" stage's current version.

        return ce.getValue();
    }
    
    // You can access a cache entry for update any time after it has been allocated.
    // This does not affect the stage.
    AbstractValue& 
    updCacheEntry(SubsystemIndex subsys, CacheEntryIndex index) const {
        const PerSubsystemInfo& ss = subsystems[subsys];
    
        SimTK_INDEXCHECK(index,(int)ss.cacheInfo.size(),"StateImpl::updCacheEntry()");
        CacheEntryInfo& ce = ss.cacheInfo[index];
    
        return ce.updValue();
    }

    bool isCacheValueRealized(SubsystemIndex subx, CacheEntryIndex cx) const {
        const PerSubsystemInfo& ss = subsystems[subx];
        SimTK_INDEXCHECK(cx,(int)ss.cacheInfo.size(),"StateImpl::isCacheValueRealized()");
        const CacheEntryInfo& ce = ss.cacheInfo[cx];
        return ce.isCurrent(getSubsystemStage(subx), getSubsystemStageVersions(subx));
    }
    void markCacheValueRealized(SubsystemIndex subx, CacheEntryIndex cx) const {
        const PerSubsystemInfo& ss = subsystems[subx];
        SimTK_INDEXCHECK(cx,(int)ss.cacheInfo.size(),"StateImpl::markCacheValueRealized()");
        CacheEntryInfo& ce = ss.cacheInfo[cx];
    
        // If this cache entry depends on anything, it can't be valid unless we're
        // at least *working* on its depends-on stage, meaning the current stage would
        // have to be the one before that. The depends-on stage is required to be at
        // least Stage::Topology, so its prev() stage exists.
        SimTK_STAGECHECK_GE(getSubsystemStage(subx), 
            ce.getDependsOnStage().prev(), "StateImpl::markCacheValueRealized()");

        ce.markAsComputed(getSubsystemStageVersions(subx));
    }

    void markCacheValueNotRealized(SubsystemIndex subx, CacheEntryIndex cx) const {
        const PerSubsystemInfo& ss = subsystems[subx];
        SimTK_INDEXCHECK(cx,(int)ss.cacheInfo.size(),
            "StateImpl::markCacheValueNotRealized()");
        CacheEntryInfo& ce = ss.cacheInfo[cx];

        ce.invalidate();
    }

    const Stage& getLowestStageModified() const {
        return lowestModifiedSystemStage;
    }
    void resetLowestStageModified() const {
        lowestModifiedSystemStage = std::min(Stage::Infinity, currentSystemStage.next());
    }

    void autoUpdateDiscreteVariables() {
        // TODO: make this more efficient
        for (SubsystemIndex subx(0); subx < subsystems.size(); ++subx) {
            PerSubsystemInfo& ss = subsystems[subx];
            Array_<DiscreteVarInfo>& dvars = ss.discreteInfo;
            for (DiscreteVariableIndex dx(0); dx < dvars.size(); ++dx) {
                DiscreteVarInfo& dinfo = dvars[dx];
                const CacheEntryIndex cx = dinfo.getAutoUpdateEntry();
                if (!cx.isValid()) continue; // not an auto-update variable
                CacheEntryInfo& cinfo = ss.cacheInfo[cx];
                if (cinfo.isCurrent(getSubsystemStage(subx), getSubsystemStageVersions(subx)))
                    cinfo.swapValue(getTime(), dinfo);
                cinfo.invalidate();
            }
        }
    }
    
    String toString() const {
        String out;
        out += "<State>\n";
    
        out += "<Real name=time>" + String(t) + "</Real>\n";
    
        out += "<Vector name=q size=" + String(q.size()) + ">";
        if (q.size()) out += "\n";
        for (int i=0; i<q.size(); ++i)
            out += String(q[i]) + "\n";
        out += "</Vector>\n";
    
        out += "<Vector name=u size=" + String(u.size()) + ">";
        if (u.size()) out += "\n";
        for (int i=0; i<u.size(); ++i)
            out += String(u[i]) + "\n";
        out += "</Vector>\n";
    
        out += "<Vector name=z size=" + String(z.size()) + ">";
        if (z.size()) out += "\n";
        for (int i=0; i<z.size(); ++i)
            out += String(z[i]) + "\n";
        out += "</Vector>\n";
    
    
        for (SubsystemIndex ss(0); ss < (int)subsystems.size(); ++ss) {
            const PerSubsystemInfo& info = subsystems[ss];
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
    
        for (SubsystemIndex ss(0); ss < (int)subsystems.size(); ++ss) {
            const PerSubsystemInfo& info = subsystems[ss];
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
    
        out += "<Vector name=qdot size=" + String(qdot.size()) + ">";
        if (qdot.size()) out += "\n";
        for (int i=0; i<qdot.size(); ++i)
            out += String(qdot[i]) + "\n";
        out += "</Vector>\n";
    
        out += "<Vector name=udot size=" + String(udot.size()) + ">";
        if (udot.size()) out += "\n";
        for (int i=0; i<udot.size(); ++i)
            out += String(udot[i]) + "\n";
        out += "</Vector>\n";
    
        out += "<Vector name=zdot size=" + String(zdot.size()) + ">";
        if (zdot.size()) out += "\n";
        for (int i=0; i<zdot.size(); ++i)
            out += String(zdot[i]) + "\n";
        out += "</Vector>\n";
    
        out += "<Vector name=qdotdot size=" + String(qdotdot.size()) + ">";
        if (qdotdot.size()) out += "\n";
        for (int i=0; i<qdotdot.size(); ++i)
            out += String(qdotdot[i]) + "\n";
        out += "</Vector>\n";
    
        out += "<Vector name=qerr size=" + String(qerr.size()) + ">";
        if (qerr.size()) out += "\n";
        for (int i=0; i<qerr.size(); ++i)
            out += String(qerr[i]) + "\n";
        out += "</Vector>\n";
    
        out += "<Vector name=uerr size=" + String(uerr.size()) + ">";
        if (uerr.size()) out += "\n";
        for (int i=0; i<uerr.size(); ++i)
            out += String(uerr[i]) + "\n";
        out += "</Vector>\n";
    
        out += "<Vector name=udoterr size=" + String(udoterr.size()) + ">";
        if (udoterr.size()) out += "\n";
        for (int i=0; i<udoterr.size(); ++i)
            out += String(udoterr[i]) + "\n";
        out += "</Vector>\n";
    
        out += "<Vector name=multipliers size=" + String(multipliers.size()) + ">";
        if (multipliers.size()) out += "\n";
        for (int i=0; i<multipliers.size(); ++i)
            out += String(multipliers[i]) + "\n";
        out += "</Vector>\n";
    
        out += "</Cache>\n";
        return out;
    }

private:
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
    Array_<PerSubsystemInfo> subsystems;

    // Return true only if all subsystems have been realized to at least Stage g.
    bool allSubsystemsAtLeastAtStage(Stage g) const {
        for (SubsystemIndex i(0); i < (int)subsystems.size(); ++i)
            if (subsystems[i].currentStage < g)
                return false;
        return true;
    }

};




//==============================================================================
//                                   STATE
//==============================================================================

State::State() {
    impl = new StateImpl();
}
// Restore state to default-constructed condition
void State::clear() {
    delete impl;
    impl = new StateImpl();
}
State::~State() {
    delete impl; impl=0;
}
// copy constructor
State::State(const State& state) {
    impl = new StateImpl(*state.impl);
}
    
// copy assignment
State& State::operator=(const State& src) {
    if (&src == this) return *this;
    if (!impl) {
        // we're defining this state here (if src is not empty)
        if (src.impl)
            impl = src.impl->clone();
        return *this;
    }

    // Assignment or redefinition
    if (src.impl) *impl = *src.impl;
    else {delete impl; impl=0;}
    return *this;
}


void State::setNumSubsystems(int i) {
    updImpl().setNumSubsystems(i);
}
void State::initializeSubsystem(SubsystemIndex subsys, const String& name, const String& version) {
    updImpl().initializeSubsystem(subsys, name, version);
}

SubsystemIndex State::addSubsystem(const String& name, const String& version) {
    return updImpl().addSubsystem(name, version);
}
int State::getNumSubsystems() const {
    return getImpl().getNumSubsystems();
}
const String& State::getSubsystemName(SubsystemIndex subsys) const {
    return getImpl().getSubsystemName(subsys);
}
const String& State::getSubsystemVersion(SubsystemIndex subsys) const {
    return getImpl().getSubsystemVersion(subsys);
}
const Stage& State::getSubsystemStage(SubsystemIndex subsys) const {
    return getImpl().getSubsystemStage(subsys);
}
const Stage& State::getSystemStage() const {
    return getImpl().getSystemStage();
}
void State::invalidateAll(Stage stage) {
    updImpl().invalidateAll(stage);
}
void State::invalidateAllCacheAtOrAbove(Stage stage) const {
    getImpl().invalidateAllCacheAtOrAbove(stage);
}
void State::advanceSubsystemToStage(SubsystemIndex subsys, Stage stage) const {
    getImpl().advanceSubsystemToStage(subsys, stage);
}
void State::advanceSystemToStage(Stage stage) const {
    getImpl().advanceSystemToStage(stage);
}
QIndex State::allocateQ(SubsystemIndex subsys, const Vector& qInit) {
    return updImpl().allocateQ(subsys, qInit);
}
UIndex State::allocateU(SubsystemIndex subsys, const Vector& uInit) {
    return updImpl().allocateU(subsys, uInit);
}
ZIndex State::allocateZ(SubsystemIndex subsys, const Vector& zInit) {
    return updImpl().allocateZ(subsys, zInit);
}
QErrIndex State::allocateQErr(SubsystemIndex subsys, int nqerr) const {
    return getImpl().allocateQErr(subsys, nqerr);
}
UErrIndex State::allocateUErr(SubsystemIndex subsys, int nuerr) const {
    return getImpl().allocateUErr(subsys, nuerr);
}
UDotErrIndex State::allocateUDotErr(SubsystemIndex subsys, int nudoterr) const {
    return getImpl().allocateUDotErr(subsys, nudoterr);
}
EventTriggerByStageIndex State::allocateEventTrigger(SubsystemIndex subsys, Stage stage, int nevent) const {
    return getImpl().allocateEventTrigger(subsys, stage, nevent);
}
DiscreteVariableIndex State::allocateDiscreteVariable(SubsystemIndex subsys, Stage stage, AbstractValue* v) {
    return updImpl().allocateDiscreteVariable(subsys, stage, v);
}
DiscreteVariableIndex
State::allocateAutoUpdateDiscreteVariable(SubsystemIndex subsys, Stage invalidates, AbstractValue* v,
                                          Stage updateDependsOn) {
    return updImpl().allocateAutoUpdateDiscreteVariable(subsys, invalidates, v, updateDependsOn); 
}
CacheEntryIndex State::allocateCacheEntry(SubsystemIndex subsys, Stage dependsOn, Stage computedBy, AbstractValue* v) const {
    return getImpl().allocateCacheEntry(subsys, dependsOn, computedBy, v);
}
int State::getNY() const {
    return getImpl().getNY();
}
SystemYIndex State::getQStart() const {
    return getImpl().getQStart();
}
int State::getNQ() const {
    return getImpl().getNQ();
}
SystemYIndex State::getUStart() const {
    return getImpl().getUStart();
}
int State::getNU() const {
    return getImpl().getNU();
}
SystemYIndex State::getZStart() const {
    return getImpl().getZStart();
}
int State::getNZ() const {
    return getImpl().getNZ();
}
int State::getNYErr() const {
    return getImpl().getNYErr();
}
SystemYErrIndex State::getQErrStart() const {
    return getImpl().getQErrStart();
}
int State::getNQErr() const {
    return getImpl().getNQErr();
}
SystemYErrIndex State::getUErrStart() const {
    return getImpl().getUErrStart();
}
int State::getNUErr() const {
    return getImpl().getNUErr();
}
int State::getNUDotErr() const {
    return getImpl().getNUDotErr();
}
int State::getNMultipliers() const {
    return getNUDotErr();
}
SystemQIndex State::getQStart(SubsystemIndex subsys) const {
    return getImpl().getQStart(subsys);
}
int State::getNQ(SubsystemIndex subsys) const {
    return getImpl().getNQ(subsys);
}
SystemUIndex State::getUStart(SubsystemIndex subsys) const {
    return getImpl().getUStart(subsys);
}
int State::getNU(SubsystemIndex subsys) const {
    return getImpl().getNU(subsys);
}
SystemZIndex State::getZStart(SubsystemIndex subsys) const {
    return getImpl().getZStart(subsys);
}
int State::getNZ(SubsystemIndex subsys) const {
    return getImpl().getNZ(subsys);
}
SystemQErrIndex State::getQErrStart(SubsystemIndex subsys) const {
    return getImpl().getQErrStart(subsys);
}
int State::getNQErr(SubsystemIndex subsys) const {
    return getImpl().getNQErr(subsys);
}
SystemUErrIndex State::getUErrStart(SubsystemIndex subsys) const {
    return getImpl().getUErrStart(subsys);
}
int State::getNUErr(SubsystemIndex subsys) const {
    return getImpl().getNUErr(subsys);
}
SystemUDotErrIndex State::getUDotErrStart(SubsystemIndex subsys) const {
    return getImpl().getUDotErrStart(subsys);
}
int State::getNUDotErr(SubsystemIndex subsys) const {
    return getImpl().getNUDotErr(subsys);
}
SystemMultiplierIndex State::getMultipliersStart(SubsystemIndex i) const {
    return SystemMultiplierIndex(getUDotErrStart(i));
}
int State::getNMultipliers(SubsystemIndex i) const {
    return getNUDotErr(i);
}
int State::getNEventTriggers() const {
    return getImpl().getNEventTriggers();
}
SystemEventTriggerIndex State::getEventTriggerStartByStage(Stage stage) const {
    return getImpl().getEventTriggerStartByStage(stage);
}
int State::getNEventTriggersByStage(Stage stage) const {
    return getImpl().getNEventTriggersByStage(stage);
}
SystemEventTriggerByStageIndex State::getEventTriggerStartByStage(SubsystemIndex subsys, Stage stage) const {
    return getImpl().getEventTriggerStartByStage(subsys, stage);
}
int State::getNEventTriggersByStage(SubsystemIndex subsys, Stage stage) const {
    return getImpl().getNEventTriggersByStage(subsys, stage);
}
const Vector& State::getEventTriggers() const {
    return getImpl().getEventTriggers();
}
const Vector& State::getEventTriggersByStage(Stage stage) const {
    return getImpl().getEventTriggersByStage(stage);
}
const Vector& State::getEventTriggersByStage(SubsystemIndex subsys, Stage stage) const {
    return getImpl().getEventTriggersByStage(subsys, stage);
}
Vector& State::updEventTriggers() const {
    return getImpl().updEventTriggers();
}
Vector& State::updEventTriggersByStage(Stage stage) const {
    return getImpl().updEventTriggersByStage(stage);
}
Vector& State::updEventTriggersByStage(SubsystemIndex subsys, Stage stage) const {
    return getImpl().updEventTriggersByStage(subsys, stage);
}
const Vector& State::getQ(SubsystemIndex subsys) const {
    return getImpl().getQ(subsys);
}
const Vector& State::getU(SubsystemIndex subsys) const {
    return getImpl().getU(subsys);
}
const Vector& State::getZ(SubsystemIndex subsys) const {
    return getImpl().getZ(subsys);
}
Vector& State::updQ(SubsystemIndex subsys) {
    return updImpl().updQ(subsys);
}
Vector& State::updU(SubsystemIndex subsys) {
    return updImpl().updU(subsys);
}
Vector& State::updZ(SubsystemIndex subsys) {
    return updImpl().updZ(subsys);
}
const Vector& State::getQDot(SubsystemIndex subsys) const {
    return getImpl().getQDot(subsys);
}
const Vector& State::getUDot(SubsystemIndex subsys) const {
    return getImpl().getUDot(subsys);
}
const Vector& State::getZDot(SubsystemIndex subsys) const {
    return getImpl().getZDot(subsys);
}
const Vector& State::getQDotDot(SubsystemIndex subsys) const {
    return getImpl().getQDotDot(subsys);
}
Vector& State::updQDot(SubsystemIndex subsys) const {
    return getImpl().updQDot(subsys);
}
Vector& State::updUDot(SubsystemIndex subsys) const {
    return getImpl().updUDot(subsys);
}
Vector& State::updZDot(SubsystemIndex subsys) const {
    return getImpl().updZDot(subsys);
}
Vector& State::updQDotDot(SubsystemIndex subsys) const {
    return getImpl().updQDotDot(subsys);
}
const Vector& State::getQErr(SubsystemIndex subsys) const {
    return getImpl().getQErr(subsys);
}
const Vector& State::getUErr(SubsystemIndex subsys) const {
    return getImpl().getUErr(subsys);
}
const Vector& State::getUDotErr(SubsystemIndex subsys) const {
    return getImpl().getUDotErr(subsys);
}
const Vector& State::getMultipliers(SubsystemIndex subsys) const {
    return getImpl().getMultipliers(subsys);
}
Vector& State::updQErr(SubsystemIndex subsys) const {
    return getImpl().updQErr(subsys);
}
Vector& State::updUErr(SubsystemIndex subsys) const {
    return getImpl().updUErr(subsys);
}
Vector& State::updUDotErr(SubsystemIndex subsys) const {
    return getImpl().updUDotErr(subsys);
}
Vector& State::updMultipliers(SubsystemIndex subsys) const {
    return getImpl().updMultipliers(subsys);
}
const Real& State::getTime() const {
    return getImpl().getTime();
}
const Vector& State::getY() const {
    return getImpl().getY();
}
const Vector& State::getQ() const {
    return getImpl().getQ();
}
const Vector& State::getU() const {
    return getImpl().getU();
}
const Vector& State::getZ() const {
    return getImpl().getZ();
}
Real& State::updTime() {
    return updImpl().updTime();
}
Vector& State::updY() {
    return updImpl().updY();
}
void State::setTime(Real t) {
    updTime() = t;
}
void State::setY(const Vector& y) {
    updY() = y;
}
Vector& State::updQ() {
    return updImpl().updQ();
}
Vector& State::updU() {
    return updImpl().updU();
}
Vector& State::updZ() {
    return updImpl().updZ();
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
    return getImpl().getYDot();
}
const Vector& State::getQDot() const {
    return getImpl().getQDot();
}
const Vector& State::getZDot() const {
    return getImpl().getZDot();
}
const Vector& State::getUDot() const {
    return getImpl().getUDot();
}
const Vector& State::getQDotDot() const {
    return getImpl().getQDotDot();
}
Vector& State::updYDot() const {
    return getImpl().updYDot();
}
Vector& State::updQDot() const {
    return getImpl().updQDot();
}
Vector& State::updZDot() const {
    return getImpl().updZDot();
}
Vector& State::updUDot() const {
    return getImpl().updUDot();
}
Vector& State::updQDotDot() const {
    return getImpl().updQDotDot();
}
const Vector& State::getYErr() const {
    return getImpl().getYErr();
}
const Vector& State::getQErr() const {
    return getImpl().getQErr();
}
const Vector& State::getUErr() const {
    return getImpl().getUErr();
}
const Vector& State::getUDotErr() const {
    return getImpl().getUDotErr();
}
const Vector& State::getMultipliers() const {
    return getImpl().getMultipliers();
}
Vector& State::updYErr() const {
    return getImpl().updYErr();
}
Vector& State::updQErr() const {
    return getImpl().updQErr();
}
Vector& State::updUErr() const {
    return getImpl().updUErr();
}
Vector& State::updUDotErr() const {
    return getImpl().updUDotErr();
}
Vector& State::updMultipliers() const {
    return getImpl().updMultipliers();
}



CacheEntryIndex State::getDiscreteVarUpdateIndex(SubsystemIndex subsys, DiscreteVariableIndex index) const {
    return getImpl().getDiscreteVarUpdateIndex(subsys, index);
}
Stage State::getDiscreteVarAllocationStage(SubsystemIndex subsys, DiscreteVariableIndex index) const {
    return getImpl().getDiscreteVarAllocationStage(subsys, index);
}
Stage State::getDiscreteVarInvalidatesStage(SubsystemIndex subsys, DiscreteVariableIndex index) const {
    return getImpl().getDiscreteVarInvalidatesStage(subsys, index);
}
const AbstractValue& State::getDiscreteVariable(SubsystemIndex subsys, DiscreteVariableIndex index) const {
    return getImpl().getDiscreteVariable(subsys, index);
}
Real State::getDiscreteVarLastUpdateTime(SubsystemIndex subsys, DiscreteVariableIndex index) const {
    return getImpl().getDiscreteVarLastUpdateTime(subsys, index);
}
const AbstractValue& State::getDiscreteVarUpdateValue(SubsystemIndex subsys, DiscreteVariableIndex index) const {
    return getImpl().getDiscreteVarUpdateValue(subsys, index);
}
AbstractValue& State::updDiscreteVarUpdateValue(SubsystemIndex subsys, DiscreteVariableIndex index) const {
    return getImpl().updDiscreteVarUpdateValue(subsys, index);
}
bool State::isDiscreteVarUpdateValueRealized(SubsystemIndex subsys, DiscreteVariableIndex index) const {
    return getImpl().isDiscreteVarUpdateValueRealized(subsys, index);
}
void State::markDiscreteVarUpdateValueRealized(SubsystemIndex subsys, DiscreteVariableIndex index) const {
    getImpl().markDiscreteVarUpdateValueRealized(subsys, index);
}


AbstractValue& State::updDiscreteVariable(SubsystemIndex subsys, DiscreteVariableIndex index) {
    return updImpl().updDiscreteVariable(subsys, index);
}
void State::setDiscreteVariable(SubsystemIndex i, DiscreteVariableIndex index, const AbstractValue& v) {
    updDiscreteVariable(i,index) = v;
}

Stage State::getCacheEntryAllocationStage(SubsystemIndex subsys, CacheEntryIndex index) const {
    return getImpl().getCacheEntryAllocationStage(subsys, index);
}
const AbstractValue& State::getCacheEntry(SubsystemIndex subsys, CacheEntryIndex index) const {
    return getImpl().getCacheEntry(subsys, index);
}
AbstractValue& State::updCacheEntry(SubsystemIndex subsys, CacheEntryIndex index) const {
    return getImpl().updCacheEntry(subsys, index);
}


bool State::isCacheValueRealized(SubsystemIndex subx, CacheEntryIndex cx) const {
    return getImpl().isCacheValueRealized(subx, cx); 
}
void State::markCacheValueRealized(SubsystemIndex subx, CacheEntryIndex cx) const {
    getImpl().markCacheValueRealized(subx, cx); 
}
void State::markCacheValueNotRealized(SubsystemIndex subx, CacheEntryIndex cx) const {
    getImpl().markCacheValueNotRealized(subx, cx); 
}
const Stage& State::getLowestStageModified() const {
    return getImpl().getLowestStageModified(); 
}
void State::resetLowestStageModified() const {
    getImpl().resetLowestStageModified(); 
}
void State::autoUpdateDiscreteVariables() {
    updImpl().autoUpdateDiscreteVariables(); 
}

String State::toString() const {
    return getImpl().toString();
}
String State::cacheToString() const {
    return getImpl().cacheToString();
}

std::ostream& 
operator<<(std::ostream& o, const State& s) {
    o << "STATE:" << std::endl;
    o << s.toString() << std::endl;
    o << "CACHE:" << std::endl;
    return o << s.cacheToString() << std::endl;
}

} // namespace SimTK

