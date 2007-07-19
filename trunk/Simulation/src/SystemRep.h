#ifndef SimTK_SYSTEM_GUTSREP_H_
#define SimTK_SYSTEM_GUTSREP_H_

/* Portions copyright (c) 2006-7 Stanford University and Michael Sherman.
 * Contributors:
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"
#include "SimTKcommon/internal/State.h"
#include "SimTKcommon/internal/System.h"
#include "SimTKcommon/internal/Subsystem.h"

#include "SubsystemRep.h"

namespace SimTK {

class System::GutsRep {
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
		assert(child.hasRep() && !child.isInSystem()); // TODO
        assert(child.isOwnerHandle());

        // This is a topology change.
        invalidateSystemTopologyCache();

        const SubsystemId id = SubsystemId((int)subsystems.size());
        subsystems.resize(id+1); // grow
		Subsystem& s = subsystems.back(); // refer to the empty handle

		s.setRep(child.updRep());		 // reference the passed-in rep
		s.updRep().setMyHandle(s);	     // steal ownership
        s.updRep().setSystem(*myHandle, id);
		return id;
	}

    // Default treats all state variable identically. Should be asking the 
    // subsystems. TODO
    Real calcYErrorNorm(const State& s, const Vector& y_err) const {
        assert(y_err.size() == s.getY().size());
        SimTK_STAGECHECK_GE(s.getSystemStage(), Stage::Position,
            "System::calcYErrorNorm()");
        return y_err.size()==0 ? 0 : std::sqrt( y_err.normSqr()/y_err.size() );
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
        clonep = 0;

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
        clonep = src.clonep;

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

        // TOPOLOGY STAGE CACHE //

    // This should only be true when *all* subsystems have successfully
    // completed realizeTopology(). Anything which invalidates topology for
    // one of the contained subsystem must invalidate topology for the system
    // as a whole also.
    mutable bool systemTopologyRealized;

    // This is only meaningful if topologyRealized==true.
    mutable State defaultState;

    mutable bool hasTimeAdvancedEventsFlag; //TODO: should be in State as a Model variable

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

// TODO
class StudyRep {
public:
    StudyRep(const System& sys)
      : myHandle(0), system(new System(sys))
    {
        system->realize(state, Stage::Topology);
    }

    virtual ~StudyRep() {
        delete system;
    }

    StudyRep* clone() const {
        StudyRep* dup = cloneStudyRep();
        dup->myHandle = 0;
        return dup;
    }
    virtual StudyRep* cloneStudyRep() const = 0;

    const System& getSystem() const {return *system;}
    const State&  getState()  const {return state;}
    State&        updState()        {return state;}

    void setMyHandle(Study& h) {myHandle = &h;}
    void clearMyHandle() {myHandle=0;}
private:
    friend class Study;
    Study* myHandle;     // the owner of this rep

    System* system;
    State   state;
};

} // namespace SimTK

#endif // SimTK_SYSTEM_GUTSREP_H_
