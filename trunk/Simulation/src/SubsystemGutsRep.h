#ifndef SimTK_SimTKCOMMON_SUBSYSTEM_GUTSREP_H_
#define SimTK_SimTKCOMMON_SUBSYSTEM_GUTSREP_H_

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

    std::vector<AbstractMeasure::Implementation*> measures;

        // TOPOLOGY CACHE

    mutable bool subsystemTopologyRealized;

private:
    // suppress automatic copy assignment operator
    GutsRep& operator=(const GutsRep&);
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_SUBSYSTEM_GUTSREP_H_
