#ifndef SimTK_SimTKCOMMON_SYSTEM_GUTSREP_H_
#define SimTK_SimTKCOMMON_SYSTEM_GUTSREP_H_

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


// This header is internal source code and is not part of the SimTKcommon
// API or distribution. This is the private, opaque implementation of
// the System::Guts class, which contains just a pointer to the
// object declared here.

#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"
#include "SimTKcommon/internal/State.h"
#include "SimTKcommon/internal/Subsystem.h"

#include "SimTKcommon/internal/System.h"
#include "SimTKcommon/internal/SystemGuts.h"

#include <unordered_map>

namespace SimTK {

class System::Guts::GutsRep {
public:

    GutsRep(const String& name, const String& version) 
      : myHandle(nullptr), 
        systemName(name), systemVersion(version),
        defaultTimeScale(Real(0.1)), 
        defaultLengthScale(Real(1)),
        defaultUpDirection(YAxis), 
        useUniformBackground(false),
        hasTimeAdvancedEventsFlag(false),
        topologyCacheVersion(1) // zero isn't allowed
    {
        clearCache();
        resetAllCounters();
    }

    // Default constructor invokes the one above.
    GutsRep() : GutsRep("<NONAME>", "2.2.0") {}

    GutsRep(const GutsRep& src)
    :   myHandle(nullptr),
        systemName(src.systemName), systemVersion(src.systemVersion),
        subsystems(src.subsystems),
        defaultTimeScale(src.defaultTimeScale),
        defaultLengthScale(src.defaultLengthScale),
        defaultUpDirection(src.defaultUpDirection), 
        useUniformBackground(src.useUniformBackground),
        hasTimeAdvancedEventsFlag(src.hasTimeAdvancedEventsFlag),
        topologyCacheVersion(src.topologyCacheVersion)
    {
        clearCache();
        resetAllCounters();
    }

    ~GutsRep() {
        clearMyHandle();
        subsystems.clear();
        invalidateSystemTopologyCache();
    }

    const String& getName()    const {return systemName;}
    const String& getVersion() const {return systemVersion;}

    void setDefaultTimeScale(Real tc) 
    {   defaultTimeScale = tc; }
    Real getDefaultTimeScale() const {return defaultTimeScale;}
    void setDefaultLengthScale(Real lc)
    {   defaultLengthScale = lc; }
    Real getDefaultLengthScale() const {return defaultLengthScale;}

    void setUpDirection(const CoordinateDirection& up) 
    {   defaultUpDirection = up; }
    CoordinateDirection getUpDirection() const {return defaultUpDirection;}
    void setUseUniformBackground(bool useUniform)
    {   useUniformBackground = useUniform; }
    bool getUseUniformBackground() const {return useUniformBackground;}

    const State& getDefaultState() const {return defaultState;}
    State&       updDefaultState()       {return defaultState;}

    int              getNumSubsystems()             const {return (int)subsystems.size();}
    const Subsystem& getSubsystem(SubsystemIndex i) const {return subsystems[i];}
    Subsystem&       updSubsystem(SubsystemIndex i)       {return subsystems[i];}


    // Take over ownership from the Subsystem handle, allocate a new
    // subsystem slot for it, and return the slot number. This is only 
    // allowed if the supplied Subsystem already has a rep, but is
    // NOT part of some other System.
    SubsystemIndex adoptSubsystem(Subsystem& child) {
        assert(child.hasGuts() && !child.isInSystem()); // TODO
        assert(child.isOwnerHandle());

        // This is a topology change.
        invalidateSystemTopologyCache();

        const SubsystemIndex id = SubsystemIndex((int)subsystems.size());
        subsystems.resize(id+1); // grow
        Subsystem& s = subsystems.back(); // refer to the empty handle

        // Take over ownership of the child's guts, leaving the child
        // as a non-owner reference to the same guts.
        s.adoptSubsystemGuts(&child.updSubsystemGuts());
        s.setSystem(*myHandle, id);

        return id;
    }

    void setMyHandle(System& h) {myHandle = &h;}
    void clearMyHandle() {myHandle=nullptr;}

    bool systemTopologyHasBeenRealized() const 
    {   return systemTopologyRealized; }

    StageVersion getSystemTopologyCacheVersion() const
    {   return topologyCacheVersion; }

    // Use this cautiously if at all!
    void setSystemTopologyCacheVersion(StageVersion version) const {
        assert(version>0); 
        auto mThis = const_cast<System::Guts::GutsRep*>(this);
        mThis->topologyCacheVersion = version; 
    }

    // Invalidating the System topology cache requires invalidating all
    // Subsystem topology caches also so that the next realizeTopology()
    // will cause them to request their State resources again so we can
    // build up the defaultState.
    void invalidateSystemTopologyCache() const {
        if (systemTopologyRealized) {
            auto mThis = const_cast<System::Guts::GutsRep*>(this);
            mThis->invalidateCache();
        }
    }

    // This is a Topology-stage computation to be invoked only from
    // realizeTopology() calls.
    EventId createNewEventId(SubsystemIndex ssx) const {
        auto mThis = const_cast<System::Guts::GutsRep*>(this);
        mThis->eventOwnerMap[nextAvailableEventId] = ssx;
        return mThis->nextAvailableEventId++;
    }

    SubsystemIndex findEventIdOwner(EventId id) const {
        auto owner = eventOwnerMap.find(id);
        if (owner == eventOwnerMap.end())
            return SubsystemIndex(); // invalid
        return owner->second;
    }

private:
friend class System;
friend class System::Guts;

    void invalidateCache() {
        // Mark system topology invalid *first* so that the invalidate
        // subsystem calls below don't recurse back here!
        clearCache();
        for (SubsystemIndex i(0); i < (int)subsystems.size(); ++i)
            subsystems[i].invalidateSubsystemTopologyCache();
    }

    void clearCache() {
        // Mark system topology invalid *first* so that the invalidate
        // subsystem calls below don't recurse back here!
        systemTopologyRealized = false;
        topologyCacheVersion++;
        defaultState.clear();
        nextAvailableEventId = EventId(1);
        eventOwnerMap.clear();
    }
    
    System*                 myHandle;     // the owner handle of these guts

    String                  systemName, systemVersion;
    StableArray<Subsystem>  subsystems;

        // TOPOLOGY STAGE STATE //

    Real                    defaultTimeScale;       // scaling hints
    Real                    defaultLengthScale;

    CoordinateDirection     defaultUpDirection;     // visualization hint
    bool                    useUniformBackground;   // visualization hint

    bool hasTimeAdvancedEventsFlag; //TODO: should be in State as a Model variable
       
    
    // TOPOLOGY STAGE CACHE //

    // This should only be true when *all* subsystems have successfully
    // completed realizeTopology(). Anything which invalidates topology for
    // one of the contained subsystem must invalidate topology for the system
    // as a whole also.
    bool                    systemTopologyRealized;
    StageVersion            topologyCacheVersion;

    // This is only meaningful if systemTopologyRealized==true. In that case
    // its Topology stage version will match the above. A State with a different
    // Topology version cannot be used with this Subsystem.
    State                   defaultState;

    // A hash function class for EventIds.
    struct EventIdHash {
        std::hash<int>::result_type operator()(EventId eid) const 
        {   return std::hash<int>()((int)eid); }
    };
    // This is initialized to EventId(1) whenever Topology is invalidated.
    EventId                 nextAvailableEventId;
    std::unordered_map<EventId, SubsystemIndex, EventIdHash>  
                            eventOwnerMap;

        // STATISTICS //
    mutable int nRealizationsOfStage[Stage::NValid];
    mutable int nRealizeCalls; // counts realizeTopology(), realizeModel(), realize()

    mutable int nPrescribeQCalls, nPrescribeUCalls;

    mutable int nProjectQCalls, nProjectUCalls;
    mutable int nFailedProjectQCalls, nFailedProjectUCalls;
    mutable int nQProjections, nUProjections; // the ones that did something
    mutable int nQErrEstProjections, nUErrEstProjections;

    mutable int nHandlerCallsThatChangedStage[Stage::NValid];
    mutable int nHandleEventsCalls;
    mutable int nReportEventsCalls;

    void resetAllCounters() {
        for (int i=0; i<Stage::NValid; ++i)
            nRealizationsOfStage[i] = nHandlerCallsThatChangedStage[i] = 0;
        nRealizeCalls = nPrescribeQCalls = nPrescribeUCalls = 0;
        nProjectQCalls = nProjectUCalls = 0;
        nFailedProjectQCalls = nFailedProjectUCalls = 0;
        nQProjections = nUProjections = 0;
        nQErrEstProjections = nUErrEstProjections = 0;
        nHandleEventsCalls = nReportEventsCalls = 0;
    }

};


} // namespace SimTK

#endif // SimTK_SimTKCOMMON_SYSTEM_GUTSREP_H_
