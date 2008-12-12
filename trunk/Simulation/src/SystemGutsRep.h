#ifndef SimTK_SimTKCOMMON_SYSTEM_GUTSREP_H_
#define SimTK_SimTKCOMMON_SYSTEM_GUTSREP_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-7 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
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

#include "SubsystemGutsRep.h"

namespace SimTK {

class System::Guts::GutsRep {
public:
    GutsRep() {new (this) GutsRep("<NONAME>", "0.0.0");}

    GutsRep(const String& name, const String& version) 
      : systemName(name), systemVersion(version), 
        myHandle(0),
        systemTopologyRealized(false), hasTimeAdvancedEventsFlag(false)
    {
        resetAllCounters();
    }

    GutsRep(const GutsRep& src)
    :   systemName(src.systemName),
        systemVersion(src.systemVersion),
        myHandle(0),
        subsystems(src.subsystems),
        hasTimeAdvancedEventsFlag(src.hasTimeAdvancedEventsFlag),
        systemTopologyRealized(false)
    {
        resetAllCounters();
    }


    ~GutsRep() {
        clearMyHandle();
        subsystems.clear();
        invalidateSystemTopologyCache();
    }

    const String& getName()    const {return systemName;}
    const String& getVersion() const {return systemVersion;}

    const State& getDefaultState() const {return defaultState;}
    State&       updDefaultState()       {return defaultState;}

    int              getNSubsystems()               const {return (int)subsystems.size();}
    const Subsystem& getSubsystem(SubsystemIndex i) const {return subsystems[i];}
    Subsystem&       updSubsystem(SubsystemIndex i)       {return subsystems[i];}


	// Take over ownership from the Subsystem handle, allocate a new
    // subsystem slot for it, and return the slot number. This is only 
    // allowed if the supplied Subsytem already has a rep, but is
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
    void clearMyHandle() {myHandle=0;}


    bool systemTopologyHasBeenRealized() const {
        return systemTopologyRealized;
    }

    void invalidateSystemTopologyCache() const {
        systemTopologyRealized = false;
        defaultState.clear();
    }

protected:
    String systemName;
    String systemVersion;
	StableArray<Subsystem> subsystems;

private:
    friend class System;
    friend class System::Guts;
    System* myHandle;     // the owner handle of these guts

        // TOPOLOGY STAGE STATE //

    bool hasTimeAdvancedEventsFlag; //TODO: should be in State as a Model variable

        // TOPOLOGY STAGE CACHE //

    // This should only be true when *all* subsystems have successfully
    // completed realizeTopology(). Anything which invalidates topology for
    // one of the contained subsystem must invalidate topology for the system
    // as a whole also.
    mutable bool systemTopologyRealized;

    // This is only meaningful if systemTopologyRealized==true.
    mutable State defaultState;

        // STATISTICS //
    mutable long nRealizationsOfStage[Stage::NValid];
    mutable long nRealizeCalls; // counts realizeTopology(), realizeModel(), realize()

    mutable long nQProjections, nUProjections;
    mutable long nQErrEstProjections, nUErrEstProjections;
    mutable long nProjectCalls;

    mutable long nHandlerCallsThatChangedStage[Stage::NValid];
    mutable long nHandleEventsCalls;
    mutable long nReportEventsCalls;

    void resetAllCounters() {
        for (int i=0; i<Stage::NValid; ++i)
            nRealizationsOfStage[i] = nHandlerCallsThatChangedStage[i] = 0;
        nRealizeCalls = nProjectCalls = nHandleEventsCalls = nReportEventsCalls = 0;
        nQProjections = nUProjections = 0;
        nQErrEstProjections = nUErrEstProjections = 0;
    }

};


////////////////////////////
// EVENT TRIGGER INFO REP //
////////////////////////////

class System::EventTriggerInfoRep {
public:
    explicit EventTriggerInfoRep(System::EventTriggerInfo* h)
      : myHandle(h), eventId(EventId(InvalidIndex)), triggerOnRising(true), triggerOnFalling(true), localizationWindow(0.1)
    {
        assert(h);
    }

private:
    System::EventTriggerInfo* myHandle;
    friend class System::EventTriggerInfo;

    EventId  eventId;
    bool triggerOnRising;
    bool triggerOnFalling;
    Real localizationWindow;
};


} // namespace SimTK

#endif // SimTK_SimTKCOMMON_SYSTEM_GUTSREP_H_
