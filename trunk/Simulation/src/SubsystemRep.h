#ifndef SimTK_SUBSYSTEM_REP_H_
#define SimTK_SUBSYSTEM_REP_H_

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

namespace SimTK {

class SubsystemRep {
public:
	SubsystemRep(const String& name, const String& version) 
      : subsystemName(name), subsystemVersion(version),
        mySystem(0), mySubsystemId(InvalidSubsystemId), myHandle(0),
        privateImplementation(0), subsystemTopologyRealized(false)
    { 
        clearAllFunctionPointers();
    }

    SubsystemRep(const SubsystemRep& src) {
        subsystemName = src.subsystemName;
        subsystemVersion = src.subsystemVersion;
        mySystem = 0;
        mySubsystemId = InvalidSubsystemId;
        myHandle = 0;
        privateImplementation = 0;
        if (src.privateImplementation && src.clonePrivateImplementationp) {
            privateImplementation = 
                src.clonePrivateImplementationp(src.privateImplementation);
        }
        subsystemTopologyRealized = false;
        copyAllFunctionPointers(src);
    }

    ~SubsystemRep() { 
        if (privateImplementation && destructPrivateImplementationp) {
            destructPrivateImplementationp(privateImplementation);
            privateImplementation=0;
        }
    }

    void adoptPrivateImplementation
       (Subsystem::PrivateImplementation* p,
        Subsystem::ClonePrivateImplementation clone,
        Subsystem::DestructPrivateImplementation destruct)
    {
        SimTK_ASSERT_ALWAYS(p && clone && destruct, 
            "Subsystem::adoptPrivateImplementation(): incomplete specification");
        privateImplementation = p;
        clonePrivateImplementationp = clone;
        destructPrivateImplementationp = destruct;
    }

    const Subsystem::PrivateImplementation& getPrivateImplementation() const {
        SimTK_ASSERT(privateImplementation,
            "Subsystem::getPrivateImplementation()");
        return *privateImplementation;
    }

    Subsystem::PrivateImplementation& updPrivateImplementation() {
        SimTK_ASSERT(privateImplementation,
            "Subsystem::updPrivateImplementation()");
        return *privateImplementation;
    }

    const String& getName()    const {return subsystemName;}
    const String& getVersion() const {return subsystemVersion;}

    void invalidateSubsystemTopologyCache() const;

    bool subsystemTopologyCacheHasBeenRealized() const {
        return subsystemTopologyRealized;
    }

	bool isInSystem() const {return mySystem != 0;}
	bool isInSameSystem(const Subsystem& otherSubsystem) const {
		return isInSystem() && otherSubsystem.isInSystem()
            && getSystem().isSameSystem(otherSubsystem.getSystem());
	}

	const System& getSystem() const {
        SimTK_ASSERT(isInSystem(), "Subsystem::getSystem()");
		return *mySystem;
	}
	System& updSystem() {
        SimTK_ASSERT(isInSystem(), "Subsystem::updSystem()");
		return *mySystem;
	}
	void setSystem(System& sys, SubsystemId id) {
        SimTK_ASSERT(!isInSystem(), "Subsystem::setSystem()");
        SimTK_ASSERT(id.isValid(), "Subsystem::setSystem()");
		mySystem = &sys;
		mySubsystemId = id;
	}
	SubsystemId getMySubsystemId() const {
		SimTK_ASSERT(isInSystem(), "Subsystem::getMySubsystemId()");
		return mySubsystemId;
	}

    void setMyHandle(Subsystem& h) {myHandle = &h;}
    const Subsystem& getMyHandle() const {assert(myHandle); return *myHandle;}
    Subsystem& updMyHandle() {assert(myHandle); return *myHandle;}
    void clearMyHandle() {myHandle=0;}

private:
    String      subsystemName;
    String      subsystemVersion;
	System*     mySystem;       // the System to which this Subsystem belongs
	SubsystemId mySubsystemId;  // Subsystem # within System

    friend class Subsystem;
    Subsystem* myHandle;	// the owner handle of this rep

    // Private implementation of concrete subsystem, if any.
    Subsystem::PrivateImplementation* privateImplementation;

        // POINTERS TO CLIENT-SIDE FUNCTION LOCATORS

        // This is a virtual function table, but the addresses are
        // determined at run time so that we don't have to depend on a
        // particular ordering in the client side virtual function table.

    Subsystem::RealizeWritableStateImplLocator         realizeTopologyp;
    Subsystem::RealizeWritableStateImplLocator         realizeModelp;
    Subsystem::RealizeConstStateImplLocator            realizeInstancep;
    Subsystem::RealizeConstStateImplLocator            realizeTimep;
    Subsystem::RealizeConstStateImplLocator            realizePositionp;
    Subsystem::RealizeConstStateImplLocator            realizeVelocityp;
    Subsystem::RealizeConstStateImplLocator            realizeDynamicsp;
    Subsystem::RealizeConstStateImplLocator            realizeAccelerationp;
    Subsystem::RealizeConstStateImplLocator            realizeReportp;

    Subsystem::CalcUnitWeightsImplLocator                   calcQUnitWeightsp;
    Subsystem::CalcUnitWeightsImplLocator                   calcUUnitWeightsp;
    Subsystem::CalcUnitWeightsImplLocator                   calcZUnitWeightsp;
    Subsystem::CalcUnitWeightsImplLocator                   calcQErrUnitTolerancesp;
    Subsystem::CalcUnitWeightsImplLocator                   calcUErrUnitTolerancesp;
    Subsystem::CalcDecorativeGeometryAndAppendImplLocator   calcDecorativeGeometryAndAppendp;
    Subsystem::CloneImplLocator                             clonep;

        // These routines allow us to manipulate the concrete subsystem's
        // private implementation which we store here but otherwise ignore.
    Subsystem::ClonePrivateImplementation    clonePrivateImplementationp;
    Subsystem::DestructPrivateImplementation destructPrivateImplementationp;

    void clearAllFunctionPointers() {
        realizeTopologyp = 0;
        realizeModelp = 0;
        realizeInstancep = 0;
        realizeTimep = 0;
        realizePositionp = 0;
        realizeVelocityp = 0;
        realizeDynamicsp = 0;
        realizeAccelerationp = 0;
        realizeReportp = 0;

        calcQUnitWeightsp = 0;
        calcUUnitWeightsp = 0;
        calcZUnitWeightsp = 0;
        calcQErrUnitTolerancesp = 0;
        calcUErrUnitTolerancesp = 0;
        calcDecorativeGeometryAndAppendp = 0;
        clonep = 0;

        clonePrivateImplementationp = 0;
        destructPrivateImplementationp = 0;
    }

    void copyAllFunctionPointers(const SubsystemRep& src) {
        realizeTopologyp = src.realizeTopologyp;
        realizeModelp = src.realizeModelp;
        realizeInstancep = src.realizeInstancep;
        realizeTimep = src.realizeTimep;
        realizePositionp = src.realizePositionp;
        realizeVelocityp = src.realizeVelocityp;
        realizeDynamicsp = src.realizeDynamicsp;
        realizeAccelerationp = src.realizeAccelerationp;
        realizeReportp = src.realizeReportp;

        calcQUnitWeightsp = src.calcQUnitWeightsp;
        calcUUnitWeightsp = src.calcUUnitWeightsp;
        calcZUnitWeightsp = src.calcZUnitWeightsp;
        calcQErrUnitTolerancesp = src.calcQErrUnitTolerancesp;
        calcUErrUnitTolerancesp = src.calcUErrUnitTolerancesp;
        calcDecorativeGeometryAndAppendp = src.calcDecorativeGeometryAndAppendp;
        clonep = src.clonep;

        clonePrivateImplementationp = src.clonePrivateImplementationp;
        destructPrivateImplementationp = src.destructPrivateImplementationp;
    }

        // TOPOLOGY CACHE

    mutable bool subsystemTopologyRealized;

private:
    // suppress automatic copy assignment operator
    SubsystemRep& operator=(const SubsystemRep&);
};

} // namespace SimTK

#endif // SimTK_SUBSYSTEM_REP_H_
