/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-12 Stanford University and the Authors.        *
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

#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"
#include "SimTKcommon/internal/Event.h"
#include "SimTKcommon/internal/State.h"

#include <cassert>
#include <algorithm>
#include <ostream>
#include <set>
#include <mutex>
#include <thread>

namespace SimTK {


//==============================================================================
//                                   STATE
//==============================================================================

State::State() {
    impl = new StateImpl();
}

// Restore state to default-constructed condition
void State::clear() {
    delete impl; impl = new StateImpl();
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


// See StateImpl.h for inline method implementations.


std::ostream& 
operator<<(std::ostream& o, const State& s) {
    o << "STATE:" << std::endl;
    o << s.toString() << std::endl;
    o << "CACHE:" << std::endl;
    return o << s.cacheToString() << std::endl;
}


//==============================================================================
//                          PER SUBSYSTEM INFO
//==============================================================================

// These static methods implement a specialized stacking mechanism for State
// resources that can be allocated at different stages, which we'll call an
// "allocation stack". A resource that was
// allocated at a later stage must be forgotten again when that stage is
// subsequently invalidated, and keeping their allocations stacked by
// stage allows that to be done efficiently.
//
// The method are templatized and expect the stacks to be in Arrays
// of the same template. The template value must be a type that supports
// three methods (the template analog to virtual functions):
//      deepAssign()            a non-shallow assignment, i.e. clone the value
//      deepDestruct()          destroy any owned heap space
//      getAllocationStage()    return the stage being worked on when this was 
//                              allocated
// The template type must otherwise support shallow copy semantics so that
// the Array_ can move them around without causing any heap activity.

// Clear the contents of an allocation stack, freeing up all associated heap 
// space.
template <class T>
static void clearAllocationStack(Array_<T>& stack) {
    for (int i=stack.size()-1; i >= 0; --i)
        stack[i].deepDestruct();
    stack.clear();
}

// Resize the given allocation stack, taking care to free the heap space if the
// size is reduced.
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

// Make this allocation stack the same as the source, copying only through the 
// given stage.
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

void PerSubsystemInfo::clearContinuousVars() {
    clearAllocationStack(qInfo); 
    clearAllocationStack(uInfo);                                
    clearAllocationStack(zInfo); 
}

void PerSubsystemInfo::clearConstraintErrs() {
    clearAllocationStack(qerrInfo);
    clearAllocationStack(uerrInfo);
    clearAllocationStack(udoterrInfo);
}

void PerSubsystemInfo::clearDiscreteVars()
{   clearAllocationStack(discreteInfo); }
void PerSubsystemInfo::clearEventTriggers(int g)
{   clearAllocationStack(triggerInfo[g]); }
void PerSubsystemInfo::clearCache()
{   clearAllocationStack(cacheInfo); }

void PerSubsystemInfo::clearAllStacks() {
    clearContinuousVars(); clearDiscreteVars();
    clearConstraintErrs(); clearCache();
    for (int i=0; i < Stage::NValid; ++i)
        clearEventTriggers(i);
}

void PerSubsystemInfo::popContinuousVarsBackToStage(const Stage& g) { 
    popAllocationStackBackToStage(qInfo,g);
    popAllocationStackBackToStage(uInfo,g);
    popAllocationStackBackToStage(zInfo,g);
}
void PerSubsystemInfo::popDiscreteVarsBackToStage(const Stage& g) 
{   popAllocationStackBackToStage(discreteInfo,g); }

void PerSubsystemInfo::popConstraintErrsBackToStage(const Stage& g) {
    popAllocationStackBackToStage(qerrInfo,g);
    popAllocationStackBackToStage(uerrInfo,g);
    popAllocationStackBackToStage(udoterrInfo,g);
}
void PerSubsystemInfo::popCacheBackToStage(const Stage& g) 
{   popAllocationStackBackToStage(cacheInfo,g); }

void PerSubsystemInfo::popEventTriggersBackToStage(const Stage& g) {
    for (int i=0; i < Stage::NValid; ++i)
        popAllocationStackBackToStage(triggerInfo[i],g); 
}

void PerSubsystemInfo::popAllStacksBackToStage(const Stage& g) {
    popContinuousVarsBackToStage(g);
    popDiscreteVarsBackToStage(g);
    popConstraintErrsBackToStage(g);
    popCacheBackToStage(g);
    popEventTriggersBackToStage(g);
}

void PerSubsystemInfo::copyContinuousVarInfoThroughStage
   (const Array_<ContinuousVarInfo>& src, const Stage& g,
    Array_<ContinuousVarInfo>& dest)
{   copyAllocationStackThroughStage(dest, src, g); }

void PerSubsystemInfo::copyDiscreteVarsThroughStage
   (const Array_<DiscreteVarInfo>& src, const Stage& g)
{   copyAllocationStackThroughStage(discreteInfo, src, g); }

void PerSubsystemInfo::copyConstraintErrInfoThroughStage
   (const Array_<ConstraintErrInfo>& src, const Stage& g,
    Array_<ConstraintErrInfo>& dest)
{   copyAllocationStackThroughStage(dest, src, g); }

void PerSubsystemInfo::copyCacheThroughStage
   (const Array_<CacheEntryInfo>& src, const Stage& g)
{   copyAllocationStackThroughStage(cacheInfo, src, g); }

void PerSubsystemInfo::copyEventsThroughStage
   (const Array_<TriggerInfo>& src, const Stage& g,
    Array_<TriggerInfo>& dest)
{   copyAllocationStackThroughStage(dest, src, g); }

void PerSubsystemInfo::copyAllStacksThroughStage
   (const PerSubsystemInfo& src, const Stage& g)
{
    copyContinuousVarInfoThroughStage(src.qInfo, g, qInfo);
    copyContinuousVarInfoThroughStage(src.uInfo, g, uInfo);
    copyContinuousVarInfoThroughStage(src.zInfo, g, zInfo);

    copyDiscreteVarsThroughStage(src.discreteInfo, g);

    copyConstraintErrInfoThroughStage(src.qerrInfo,    g, qerrInfo);
    copyConstraintErrInfoThroughStage(src.uerrInfo,    g, uerrInfo);
    copyConstraintErrInfoThroughStage(src.udoterrInfo, g, udoterrInfo);

    copyCacheThroughStage(src.cacheInfo, g);
    for (int i=0; i < Stage::NValid; ++i)
        copyEventsThroughStage(src.triggerInfo[i], g, triggerInfo[i]);
}

void PerSubsystemInfo::restoreToStage(Stage g) {
    if (currentStage <= g)
        return;

    if (g < Stage::Instance) {
        clearReferencesToInstanceStageGlobals();
    }

    if (g < Stage::Model) {
        clearReferencesToModelStageGlobals();
    }

    if (g == Stage::Empty) {
        // Throw out everything, reset stage versions to 1. Leave
        // name and version alone.
        initialize();
        return;
    }

    // Backup all the allocation stacks.
    popAllStacksBackToStage(g);

    // Raise the version number for every stage that we're invalidating.
    for (int i=currentStage; i > g; --i)
        stageVersions[i]++;
    currentStage = g;
}

void PerSubsystemInfo::copyFrom(const PerSubsystemInfo& src, Stage maxStage) {
    const Stage targetStage = std::min<Stage>(src.currentStage, maxStage);

    // Forget any references to global resources.
    clearReferencesToInstanceStageGlobals();
    clearReferencesToModelStageGlobals();

    // Make sure destination state doesn't have anything past targetStage.
    restoreToStage(targetStage);

    name     = src.name;
    version  = src.version;
    copyAllStacksThroughStage(src, targetStage);

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


//==============================================================================
//                              STATE IMPL
//==============================================================================


//------------------------------------------------------------------------------
//                           COPY CONSTRUCTOR
//------------------------------------------------------------------------------
StateImpl::StateImpl(const StateImpl& src)
:   currentSystemStage(Stage::Empty)
{
    initializeStageVersions();
    //TODO: When you copy a state, you don't get the same mutex, is that okay? 
    // If not, how should we fix it? "mutex objects cannot be copied/moved (both 
    //the copy constructor and assignment operator are deleted for this type)."
    
    // Make sure that no copied cache entry could accidentally think
    // it was up to date. We'll change some of these below if appropriate.
    // (We're skipping the Empty stage 0.)
    for (int i=1; i <= src.currentSystemStage; ++i)
        systemStageVersions[i] = src.systemStageVersions[i]+1;

    subsystems = src.subsystems;
    if (src.currentSystemStage >= Stage::Topology) {
        advanceSystemToStage(Stage::Topology);
        systemStageVersions[Stage::Topology] = 
            src.systemStageVersions[Stage::Topology];
        t = src.t;
        if (src.currentSystemStage >= Stage::Model) {
            advanceSystemToStage(Stage::Model);
            systemStageVersions[Stage::Model] = 
                src.systemStageVersions[Stage::Model];
            // careful -- don't allow reallocation
            y = src.y;
            uWeights = src.uWeights;
            zWeights = src.zWeights;
        }
        if (src.currentSystemStage >= Stage::Instance) {
            advanceSystemToStage(Stage::Instance);
            systemStageVersions[Stage::Instance] = 
                src.systemStageVersions[Stage::Instance];
            // careful -- don't allow reallocation
            qerrWeights = src.qerrWeights;
            uerrWeights = src.uerrWeights;
        }
    }
}

//------------------------------------------------------------------------------
//                           COPY ASSIGNMENT
//------------------------------------------------------------------------------
StateImpl& StateImpl::operator=(const StateImpl& src) {
    if (&src == this) return *this;
    invalidateJustSystemStage(Stage::Topology);
    for (SubsystemIndex i(0); i<(int)subsystems.size(); ++i)
        subsystems[i].invalidateStageJustThisSubsystem(Stage::Topology);

    // Make sure that no copied cache entry could accidentally think
    // it was up to date. We'll change some of these below if appropriate.
    // (We're skipping the Empty stage 0.)
    for (int i=1; i <= src.currentSystemStage; ++i)
        systemStageVersions[i] = src.systemStageVersions[i]+1;

    subsystems = src.subsystems;
    if (src.currentSystemStage >= Stage::Topology) {
        advanceSystemToStage(Stage::Topology);
        systemStageVersions[Stage::Topology] = 
            src.systemStageVersions[Stage::Topology];
        t = src.t;
        if (src.currentSystemStage >= Stage::Model) {
            advanceSystemToStage(Stage::Model);
            systemStageVersions[Stage::Model] = 
                src.systemStageVersions[Stage::Model];
            // careful -- don't allow reallocation
            y = src.y;
            uWeights = src.uWeights;
            zWeights = src.zWeights;
        }
        if (src.currentSystemStage >= Stage::Instance) {
            advanceSystemToStage(Stage::Instance);
            systemStageVersions[Stage::Instance] = 
                src.systemStageVersions[Stage::Instance];
            // careful -- don't allow reallocation
            qerrWeights = src.qerrWeights;
            uerrWeights = src.uerrWeights;
        }
    }
    return *this;
}

//------------------------------------------------------------------------------
//                     INVALIDATE JUST SYSTEM STAGE
//------------------------------------------------------------------------------
// Back up the System stage just before stg if it thinks
// it is already at stg or beyond. Note that we may be backing up
// over many stages here. Careful: invalidating the stage
// for the system must also invalidate the same stage for all
// the subsystems (because we trash the shared resource pool
// here if we back up earlier than Stage::Model) but we don't
// take care of that here. Also, you can't invalidate Stage::Empty.
void StateImpl::invalidateJustSystemStage(Stage stg) {
    assert(stg > Stage::Empty);
    if (currentSystemStage < stg)
        return;

    if (currentSystemStage >= Stage::Instance && Stage::Instance >= stg) {
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
        qerrWeights.unlockShape(); qerrWeights.clear();
        uerrWeights.unlockShape(); uerrWeights.clear();
        udoterr.unlockShape();     udoterr.clear();
        multipliers.unlockShape(); multipliers.clear();
        allTriggers.unlockShape(); allTriggers.clear();
    }
    if (currentSystemStage >= Stage::Model && Stage::Model >= stg) {
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
        uWeights.unlockShape(); uWeights.clear();
        zWeights.unlockShape(); zWeights.clear();

        qdot.clear(); udot.clear(); zdot.clear();   // ydot views
        ydot.unlockShape();        ydot.clear();    // ydot data
        qdotdot.unlockShape();     qdotdot.clear(); // qdotdot data (no views)
    }
    if (currentSystemStage >= Stage::Topology && Stage::Topology >= stg) {
        // We're invalidating the topology stage. Time is considered
        // a topology stage variable so needs to be invalidated here.
        t = NaN;
    }

    // Raise the version number for every stage that we're invalidating and
    // set the current System Stage one lower than the one being invalidated.
    for (int i=currentSystemStage; i >= stg; --i)
        systemStageVersions[i]++;
    currentSystemStage = stg.prev();
}

//------------------------------------------------------------------------------
//                         ADVANCE SYSTEM TO STAGE
//------------------------------------------------------------------------------
// Advance the System stage from stg-1 to stg. It is a fatal error if
// we're not already at stg-1, and you can't advance to Stage::Empty.
// Also, you can't advance the system to stg unless ALL subsystems have
// already gotten there.
void StateImpl::advanceSystemToStage(Stage stg) const {
    assert(stg > Stage::Empty);
    assert(currentSystemStage == stg.prev());
    assert(allSubsystemsAtLeastAtStage(stg));

    if (stg == Stage::Topology) {
        // As the final "Topology" step, initialize time to 0 (it's NaN 
        // before this).
        const_cast<StateImpl*>(this)->t = 0;
    }
    else if (stg == Stage::Model) {
        // We know the shared state pool sizes now. Allocate the
        // states and matching shared cache pools.
        int nq=0, nu=0, nz=0; // total sizes
        Array_<int> ssnq(subsystems.size(), 0); // per subsystem sizes
        Array_<int> ssnu(subsystems.size(), 0);
        Array_<int> ssnz(subsystems.size(), 0);

        // Count up all 
        for (SubsystemIndex i(0); i<subsystems.size(); ++i) {
            const PerSubsystemInfo& ss = subsystems[i];
            for (unsigned j=0; j<ss.qInfo.size(); ++j)
                ssnq[i] += ss.qInfo[j].getNumVars();
            nq += ssnq[i];
            for (unsigned j=0; j<ss.uInfo.size(); ++j)
                ssnu[i] += ss.uInfo[j].getNumVars();
            nu += ssnu[i];
            for (unsigned j=0; j<ss.zInfo.size(); ++j)
                ssnz[i] += ss.zInfo[j].getNumVars();
            nz += ssnz[i];
        }

        // Allocate the actual shared state variables & cache 
        // entries and make sure no one can accidentally change the size.
        // We need write access temporarily to set up the state.
        StateImpl* wThis = const_cast<StateImpl*>(this);
        wThis->y.resize(nq+nu+nz);      wThis->y.lockShape();
        wThis->uWeights.resize(nu);     wThis->uWeights.lockShape();
        wThis->zWeights.resize(nz);     wThis->zWeights.lockShape();

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
            PerSubsystemInfo& ss = 
                const_cast<PerSubsystemInfo&>(subsystems[i]);
            const int nq=ssnq[i], nu=ssnu[i], nz=ssnz[i];

            // Assign the starting indices.
            ss.qstart=nxtq; ss.ustart=nxtu; ss.zstart=nxtz;
 
            // Build the views.
            ss.q.viewAssign(wThis->q(nxtq, nq));
            int nxt=0;
            for (unsigned j=0; j<ss.qInfo.size(); ++j) {
                const int nv = ss.qInfo[j].getNumVars();
                ss.q(nxt, nv) = ss.qInfo[j].getInitialValues();
                nxt += nv;
            }

            ss.u.viewAssign(wThis->u(nxtu, nu)); 
            ss.uWeights.viewAssign(wThis->uWeights(nxtu, nu));
            nxt=0;
            for (unsigned j=0; j<ss.uInfo.size(); ++j) {
                const int nv = ss.uInfo[j].getNumVars();
                ss.u(nxt, nv)        = ss.uInfo[j].getInitialValues();
                ss.uWeights(nxt, nv) = ss.uInfo[j].getWeights();
                nxt += nv;
            }

            ss.z.viewAssign(wThis->z(nxtz, nz)); 
            ss.zWeights.viewAssign(wThis->zWeights(nxtz, nz));
            nxt=0;
            for (unsigned j=0; j<ss.zInfo.size(); ++j) {
                const int nv = ss.zInfo[j].getNumVars();
                ss.z(nxt, nv)        = ss.zInfo[j].getInitialValues();
                ss.zWeights(nxt, nv) = ss.zInfo[j].getWeights();
                nxt += nv;
            }

            ss.qdot.viewAssign(qdot(nxtq, nq));
            ss.qdotdot.viewAssign(qdotdot(nxtq, nq));
            ss.udot.viewAssign(udot(nxtu, nu));
            ss.zdot.viewAssign(zdot(nxtz, nz));

            // Consume the slots.
            nxtq += nq; nxtu += nu; nxtz += nz;
        }
    }
    else if (stg == Stage::Instance) {
        // We know the shared cache pool sizes now. Allocate them.

        // Global sizes.
        int nqerr=0, nuerr=0, nudoterr=0, nAllTriggers=0;
        Array_<int> ntriggers(Stage::NValid, 0);

        // Per-subsystem sizes.
        const unsigned nss = subsystems.size();
        Array_<int> ssnqerr(nss, 0), ssnuerr(nss, 0), ssnudoterr(nss, 0);
        Array_< Array_<int> > ssntriggers(nss);
        for (unsigned i=0; i<nss; ++i)
            ssntriggers[i].resize(Stage::NValid, 0);

        // Count up all 
        for (SubsystemIndex i(0); i<nss; ++i) {
            const PerSubsystemInfo& ss = subsystems[i];
            for (unsigned j=0; j<ss.qerrInfo.size(); ++j)
                ssnqerr[i] += ss.qerrInfo[j].getNumErrs();
            nqerr += ssnqerr[i];
            for (unsigned j=0; j<ss.uerrInfo.size(); ++j)
                ssnuerr[i] += ss.uerrInfo[j].getNumErrs();
            nuerr += ssnuerr[i];
            for (unsigned j=0; j<ss.udoterrInfo.size(); ++j)
                ssnudoterr[i] += ss.udoterrInfo[j].getNumErrs();
            nudoterr += ssnudoterr[i];

            Array_<int>& ssntrigs = ssntriggers[i];
            for (int g=0; g<Stage::NValid; ++g)
                for (unsigned j=0; j<ss.triggerInfo[g].size(); ++j)
                    ssntrigs[g] += ss.triggerInfo[g][j].getNumSlots();

            for (int g=0; g<Stage::NValid; ++g)
                ntriggers[g] += ssntrigs[g];
        }
        for (int g=0; g<Stage::NValid; ++g)
            nAllTriggers += ntriggers[g];

        // We need write access temporarily to set up the state.
        StateImpl* wThis = const_cast<StateImpl*>(this);
        wThis->qerrWeights.resize(nqerr); wThis->qerrWeights.lockShape();
        wThis->uerrWeights.resize(nuerr); wThis->uerrWeights.lockShape();

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
        for (int g=0; g<Stage::NValid; ++g)
            nxttrigger[g] = SystemEventTriggerByStageIndex(0);

        for (SubsystemIndex i(0); i<(int)subsystems.size(); ++i) {
            PerSubsystemInfo& ss = 
                const_cast<PerSubsystemInfo&>(subsystems[i]);
            const int nqerr=ssnqerr[i], nuerr=ssnuerr[i], 
                        nudoterr=ssnudoterr[i];
            const Array_<int>& ssntrigs = ssntriggers[i];

            // Build the views. Only weights need initialization.
            ss.qerr.viewAssign(qerr(nxtqerr, nqerr));
            ss.qerrWeights.viewAssign(wThis->qerrWeights(nxtqerr, nqerr));
            int nxt=0;
            for (unsigned j=0; j<ss.qerrInfo.size(); ++j) {
                const int nerr = ss.qerrInfo[j].getNumErrs();
                ss.qerrWeights(nxt, nerr) = ss.qerrInfo[j].getWeights();
                nxt += nerr;
            }
            ss.uerr.viewAssign(uerr(nxtuerr, nuerr));
            ss.uerrWeights.viewAssign(wThis->uerrWeights(nxtuerr, nuerr));
            nxt=0;
            for (unsigned j=0; j<ss.uerrInfo.size(); ++j) {
                const int nerr = ss.uerrInfo[j].getNumErrs();
                ss.uerrWeights(nxt, nerr) = ss.uerrInfo[j].getWeights();
                nxt += nerr;
            }

            ss.udoterr.viewAssign(udoterr(nxtudoterr, nudoterr));
            // multipliers have same partitioning as udoterr
            ss.multipliers.viewAssign(multipliers(nxtudoterr, nudoterr));

            // Assign the starting indices.
            ss.qerrstart=nxtqerr; ss.uerrstart=nxtuerr; 
            ss.udoterrstart=nxtudoterr;

            // Consume the slots.
            nxtqerr += nqerr; nxtuerr += nuerr; nxtudoterr += nudoterr;

            // Same thing for event trigger slots, but by stage.
            for (int g=0; g<Stage::NValid; ++g) {
                ss.triggerstart[g] = nxttrigger[g];
                ss.triggers[g].viewAssign
                    (triggers[g](nxttrigger[g], ssntrigs[g]));
                nxttrigger[g] += ssntrigs[g];
            }

        }
    }

    // All cases fall through to here.
    currentSystemStage = stg;
}

//------------------------------------------------------------------------------
//                     AUTO UPDATE DISCRETE VARIABLES
//------------------------------------------------------------------------------
void StateImpl::autoUpdateDiscreteVariables() {
    // TODO: make this more efficient
    for (SubsystemIndex subx(0); subx < subsystems.size(); ++subx) {
        PerSubsystemInfo& ss = subsystems[subx];
        Array_<DiscreteVarInfo>& dvars = ss.discreteInfo;
        for (DiscreteVariableIndex dx(0); dx < dvars.size(); ++dx) {
            DiscreteVarInfo& dinfo = dvars[dx];
            const CacheEntryIndex cx = dinfo.getAutoUpdateEntry();
            if (!cx.isValid()) continue; // not an auto-update variable
            CacheEntryInfo& cinfo = ss.cacheInfo[cx];
            if (cinfo.isCurrent(getSubsystemStage(subx), 
                                getSubsystemStageVersions(subx)))
                cinfo.swapValue(getTime(), dinfo);
            cinfo.invalidate();
        }
    }
}

//------------------------------------------------------------------------------
//                              TO STRING
//------------------------------------------------------------------------------
// TODO: this is imcomplete. We really need a full serialization capability
// for States; this is at most useful for debugging.
String StateImpl::toString() const {
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
    
    out += "<Vector name=uWeights size=" + String(uWeights.size()) + ">";
    if (uWeights.size()) out += "\n";
    for (int i=0; i<uWeights.size(); ++i)
        out += String(uWeights[i]) + "\n";
    out += "</Vector>\n";
    
    out += "<Vector name=zWeights size=" + String(zWeights.size()) + ">";
    if (zWeights.size()) out += "\n";
    for (int i=0; i<zWeights.size(); ++i)
        out += String(zWeights[i]) + "\n";
    out += "</Vector>\n";    
    
    for (SubsystemIndex ss(0); ss < (int)subsystems.size(); ++ss) {
        const PerSubsystemInfo& info = subsystems[ss];
        out += "<Subsystem index=" + String(ss) + " name=" + info.name 
            + " version=" + info.version + ">\n";
    
        out += "  <DISCRETE VARS TODO>\n";
        
        out += "  <Vector name=q size=" + String(info.q.size()) + ">\n";
        out += "  <Vector name=u size=" + String(info.u.size()) + ">\n";
        out += "  <Vector name=z size=" + String(info.z.size()) + ">\n";

        out += "  <Vector name=uWeights size=" + String(info.uWeights.size()) + ">\n";
        out += "  <Vector name=zWeights size=" + String(info.zWeights.size()) + ">\n";    
        out += "</Subsystem>\n";
    }
    
    out += "</State>\n";
    return out;
}
    
//------------------------------------------------------------------------------
//                             CACHE TO STRING
//------------------------------------------------------------------------------
// This is just for debugging -- you wouldn't expect to serialize the cache
// since it can always be regenerated.
String StateImpl::cacheToString() const {
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
        out += "  <Vector name=qerrWeights size=" + String(info.qerrWeights.size()) + ">\n";
        out += "  <Vector name=uerr size=" + String(info.uerr.size()) + ">\n";
        out += "  <Vector name=uerrWeights size=" + String(info.uerrWeights.size()) + ">\n";

        out += "  <Vector name=udoterr size=" + String(info.udoterr.size()) + ">\n";
        out += "  <Vector name=multipliers size=" + String(info.multipliers.size()) + ">\n";

    
        for (int j=0; j<Stage::NValid; ++j) {
            out += "  <Vector name=triggers[";
            out += Stage(j).getName();
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
    
    out += "<Vector name=qerrWeights size=" + String(qerrWeights.size()) + ">";
    if (qerrWeights.size()) out += "\n";
    for (int i=0; i<qerrWeights.size(); ++i)
        out += String(qerrWeights[i]) + "\n";
    out += "</Vector>\n";

    out += "<Vector name=uerr size=" + String(uerr.size()) + ">";
    if (uerr.size()) out += "\n";
    for (int i=0; i<uerr.size(); ++i)
        out += String(uerr[i]) + "\n";
    out += "</Vector>\n";
        
    out += "<Vector name=uerrWeights size=" + String(uerrWeights.size()) + ">";
    if (uerrWeights.size()) out += "\n";
    for (int i=0; i<uerrWeights.size(); ++i)
        out += String(uerrWeights[i]) + "\n";
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


} // namespace SimTK

