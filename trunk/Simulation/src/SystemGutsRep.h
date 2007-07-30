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
    GutsRep() 
      : systemName("<NONAME>"), systemVersion("0.0.0"), 
        myHandle(0),
        systemTopologyRealized(false), hasTimeAdvancedEventsFlag(false)
    {
        clearAllFunctionPointers();
    }
    GutsRep(const String& name, const String& version) 
      : systemName(name), systemVersion(version), 
        myHandle(0),
        systemTopologyRealized(false), hasTimeAdvancedEventsFlag(false)
    {
    }

    GutsRep(const GutsRep& src) {
        systemName = src.systemName;
        systemVersion = src.systemVersion;
        myHandle = 0;
        subsystems = src.subsystems;
        copyAllFunctionPointers(src);
        hasTimeAdvancedEventsFlag = src.hasTimeAdvancedEventsFlag;
        systemTopologyRealized = false;
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

    int              getNSubsystems()            const {return subsystems.size();}
    const Subsystem& getSubsystem(SubsystemId i) const {return subsystems[i];}
    Subsystem&       updSubsystem(SubsystemId i)       {return subsystems[i];}


	// Take over ownership from the Subsystem handle, allocate a new
    // subsystem slot for it, and return the slot number. This is only 
    // allowed if the supplied Subsytem already has a rep, but is
    // NOT part of some other System.
	SubsystemId adoptSubsystem(Subsystem& child) {
		assert(child.hasGuts() && !child.isInSystem()); // TODO
        assert(child.isOwnerHandle());

        // This is a topology change.
        invalidateSystemTopologyCache();

        const SubsystemId id = SubsystemId((int)subsystems.size());
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

        // POINTERS TO CLIENT-SIDE FUNCTION LOCATORS

        // This is a virtual function table, but the addresses are
        // determined at run time so that we don't have to depend on a
        // particular ordering in the client side virtual function table.

    System::Guts::DestructImplLocator                     destructp;
    System::Guts::CloneImplLocator                        clonep;

    System::Guts::RealizeWritableStateImplLocator         realizeTopologyp;
    System::Guts::RealizeWritableStateImplLocator         realizeModelp;
    System::Guts::RealizeConstStateImplLocator            realizeInstancep;
    System::Guts::RealizeConstStateImplLocator            realizeTimep;
    System::Guts::RealizeConstStateImplLocator            realizePositionp;
    System::Guts::RealizeConstStateImplLocator            realizeVelocityp;
    System::Guts::RealizeConstStateImplLocator            realizeDynamicsp;
    System::Guts::RealizeConstStateImplLocator            realizeAccelerationp;
    System::Guts::RealizeConstStateImplLocator            realizeReportp;

    System::Guts::CalcTimescaleImplLocator                calcTimescalep;
    System::Guts::CalcUnitWeightsImplLocator              calcYUnitWeightsp;
    System::Guts::ProjectImplLocator                      projectp;
    System::Guts::CalcUnitWeightsImplLocator              calcYErrUnitTolerancesp;
    System::Guts::HandleEventsImplLocator                 handleEventsp;
    System::Guts::CalcEventTriggerInfoImplLocator         calcEventTriggerInfop;
    System::Guts::CalcTimeOfNextScheduledEventImplLocator calcTimeOfNextScheduledEventp;

    void clearAllFunctionPointers() {
        destructp = 0;
        clonep    = 0;

        realizeTopologyp = 0;
        realizeModelp = 0;
        realizeInstancep = 0;
        realizeTimep = 0;
        realizePositionp = 0;
        realizeVelocityp = 0;
        realizeDynamicsp = 0;
        realizeAccelerationp = 0;
        realizeReportp = 0;

        calcTimescalep = 0;
        calcYUnitWeightsp = 0;
        projectp = 0;
        calcYErrUnitTolerancesp = 0;
        handleEventsp = 0;
        calcEventTriggerInfop = 0;
        calcTimeOfNextScheduledEventp = 0;
    }

    void copyAllFunctionPointers(const GutsRep& src) {
        destructp = src.destructp;
        clonep    = src.clonep;

        realizeTopologyp = src.realizeTopologyp;
        realizeModelp    = src.realizeModelp;
        realizeInstancep = src.realizeInstancep;
        realizeTimep     = src.realizeTimep;
        realizePositionp = src.realizePositionp;
        realizeVelocityp = src.realizeVelocityp;
        realizeDynamicsp = src.realizeDynamicsp;
        realizeAccelerationp = src.realizeAccelerationp;
        realizeReportp   = src.realizeReportp;

        calcTimescalep                  = src.calcTimescalep;
        calcYUnitWeightsp               = src.calcYUnitWeightsp;
        projectp                        = src.projectp;
        calcYErrUnitTolerancesp         = src.calcYErrUnitTolerancesp;
        handleEventsp                   = src.handleEventsp;
        calcEventTriggerInfop           = src.calcEventTriggerInfop;
        calcTimeOfNextScheduledEventp   = src.calcTimeOfNextScheduledEventp;
    }


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

};


////////////////////////////
// EVENT TRIGGER INFO REP //
////////////////////////////

class System::EventTriggerInfoRep {
public:
    explicit EventTriggerInfoRep(System::EventTriggerInfo* h)
      : myHandle(h), eventId(-1), triggerOnRising(true), triggerOnFalling(true),
        triggerOnZero(false), localizationWindow(0.1)
    {
        assert(h);
    }

private:
    System::EventTriggerInfo* myHandle;
    friend class System::EventTriggerInfo;

    int  eventId;
    bool triggerOnRising;
    bool triggerOnFalling;
    bool triggerOnZero;
    Real localizationWindow;
};


} // namespace SimTK

#endif // SimTK_SimTKCOMMON_SYSTEM_GUTSREP_H_
