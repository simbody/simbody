#ifndef SimTK_SimTKCOMMON_SUBSYSTEM_GUTSREP_H_
#define SimTK_SimTKCOMMON_SUBSYSTEM_GUTSREP_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-12 Stanford University and the Authors.        *
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
// the Subsystem::Guts class, which contains just a pointer to the
// object declared here.

#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"
#include "SimTKcommon/internal/State.h"
#include "SimTKcommon/internal/System.h"
#include "SimTKcommon/internal/Subsystem.h"
#include "SimTKcommon/internal/SubsystemGuts.h"
#include "SimTKcommon/internal/Measure.h"
#include "SimTKcommon/internal/MeasureImplementation.h"

#include <algorithm>

namespace SimTK {

class Subsystem::Guts::GutsRep {
    friend class Subsystem::Guts;
public:
	GutsRep(const String& name, const String& version) 
      : subsystemName(name), subsystemVersion(version),
        mySystem(0), mySubsystemIndex(InvalidSubsystemIndex), myHandle(0),
        subsystemTopologyRealized(false)
    { 
    }

    GutsRep(const GutsRep& src)
    :   subsystemName(src.subsystemName),
        subsystemVersion(src.subsystemVersion),
        mySystem(0),
        mySubsystemIndex(InvalidSubsystemIndex),
        myHandle(0),
        subsystemTopologyRealized(false)
    {
    }

    ~GutsRep() {
        for (MeasureIndex mx(0); mx < measures.size(); ++mx)
            if (measures[mx]->decrRefCount()==0) delete measures[mx];
        clearMyHandle();
        invalidateSubsystemTopologyCache();
    }


    const String& getName()    const {return subsystemName;}
    const String& getVersion() const {return subsystemVersion;}

    void invalidateSubsystemTopologyCache() const;

    bool subsystemTopologyHasBeenRealized() const {
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
	void setSystem(System& sys, SubsystemIndex id) {
        SimTK_ASSERT(!isInSystem(), "Subsystem::setSystem()");
        SimTK_ASSERT(id.isValid(), "Subsystem::setSystem()");
		mySystem = &sys;
		mySubsystemIndex = id;
	}
	SubsystemIndex getMySubsystemIndex() const {
		SimTK_ASSERT(isInSystem(), "Subsystem::getMySubsystemIndex()");
		return mySubsystemIndex;
	}

    void setMyHandle(Subsystem& h) {myHandle = &h;}
    const Subsystem& getMyHandle() const {assert(myHandle); return *myHandle;}
    Subsystem& updMyHandle() {assert(myHandle); return *myHandle;}
    void clearMyHandle() {myHandle=0;}

    AbstractMeasure getMeasure(MeasureIndex mx) const {
        assert(0 <= mx && mx < measures.size());
        return AbstractMeasure(measures[mx]);
    }

    MeasureIndex adoptMeasure(AbstractMeasure& m) {
        assert(m.hasImpl());
        // This is an expensive check if there are lots of measures.
        assert(std::find(measures.begin(), measures.end(), &m.getImpl())
                == measures.end());

        invalidateSubsystemTopologyCache();
        const MeasureIndex mx(measures.size());
        measures.push_back(&m.updImpl());
        measures.back()->incrRefCount();
        measures.back()->setSubsystem(updMyHandle(), mx);
        return mx;
    }

private:
    String      subsystemName;
    String      subsystemVersion;
	System*     mySystem;       // the System to which this Subsystem belongs
	SubsystemIndex mySubsystemIndex;  // Subsystem # within System

    friend class Subsystem;
    Subsystem* myHandle;	// the owner handle of this rep

    Array_<AbstractMeasure::Implementation*> measures;

        // TOPOLOGY CACHE

    mutable bool subsystemTopologyRealized;

private:
    // suppress automatic copy assignment operator
    GutsRep& operator=(const GutsRep&);
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_SUBSYSTEM_GUTSREP_H_
