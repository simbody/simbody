#ifndef SimTK_SIMBODY_SUBSYSTEM_REP_H_
#define SimTK_SIMBODY_SUBSYSTEM_REP_H_

/* Portions copyright (c) 2006 Stanford University and Michael Sherman.
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
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/State.h"
#include "simbody/internal/System.h"
#include "simbody/internal/Subsystem.h"

namespace SimTK {

class SubsystemRep {
public:
	SubsystemRep() 
      : subsystemName("<NONAME>"), subsystemVersion("0.0.0"),
        mySystem(0), mySubsystemIndex(-1), myHandle(0) 
    { 
    }

	SubsystemRep(const String& name, const String& version) 
      : subsystemName(name), subsystemVersion(version),
        mySystem(0), mySubsystemIndex(-1), myHandle(0) 
    { 
    }

    virtual ~SubsystemRep() { 
    }

    const String& getName()    const {return subsystemName;}
    const String& getVersion() const {return subsystemVersion;}

    void advanceToStage(const State& s, Stage g) const {
        s.advanceSubsystemToStage(getMySubsystemIndex(), g);
    }

    // These pull out the State entries which belong exclusively to
    // this Subsystem. These variables and cache entries are available
    // as soon as this subsystem is at stage Modeled.
    Stage getStage(const State& s) const {
        return s.getSubsystemStage(getMySubsystemIndex());
    }
    const AbstractValue& getDiscreteVariable(const State& s, int index) const {
        return s.getDiscreteVariable(getMySubsystemIndex(), index);
    }
    // State is *not* mutable here -- must have write access to change state variables.
    AbstractValue& updDiscreteVariable(State& s, int index) const {
        return s.updDiscreteVariable(getMySubsystemIndex(), index);
    }
    const AbstractValue& getCacheEntry(const State& s, int index) const {
        return s.getCacheEntry(getMySubsystemIndex(), index);
    }
    // State is mutable here.
    AbstractValue& updCacheEntry(const State& s, int index) const {
        return s.updCacheEntry(getMySubsystemIndex(), index);
    }

    // These return views on State shared global resources. The views
    // are private to this subsystem, but the global resources themselves
    // are not allocated until the *System* advances to stage Modeled.
    // Note that there is no subsystem equivalent of the State's "Y"
    // vector because in general a subsystem's state variables will
    // not be contiguous. However, a subsystem's Q's, U's, and Z's
    // will all be contiguous within those arrays.

    const Vector& getQ(const State& s) const {return s.getQ(getMySubsystemIndex());}
    const Vector& getU(const State& s) const {return s.getU(getMySubsystemIndex());}
    const Vector& getZ(const State& s) const {return s.getZ(getMySubsystemIndex());}

    // Not mutable: must have a writable state.
    Vector& updQ(State& s) const {return s.updQ(getMySubsystemIndex());}
    Vector& updU(State& s) const {return s.updU(getMySubsystemIndex());}
    Vector& updZ(State& s) const {return s.updZ(getMySubsystemIndex());}

    const Vector& getQDot   (const State& s) const {return s.getQDot(getMySubsystemIndex());}
    const Vector& getUDot   (const State& s) const {return s.getUDot(getMySubsystemIndex());}
    const Vector& getZDot   (const State& s) const {return s.getZDot(getMySubsystemIndex());}
    const Vector& getQDotDot(const State& s) const {return s.getQDotDot(getMySubsystemIndex());}

    // These are mutable
    Vector& updQDot   (const State& s) const {return s.updQDot(getMySubsystemIndex());}
    Vector& updUDot   (const State& s) const {return s.updUDot(getMySubsystemIndex());}
    Vector& updZDot   (const State& s) const {return s.updZDot(getMySubsystemIndex());}
    Vector& updQDotDot(const State& s) const {return s.updQDotDot(getMySubsystemIndex());}

    SubsystemRep* clone() const {
		assert(!isInSystem()); // TODO
        SubsystemRep* dup = cloneSubsystemRep();
        dup->myHandle = 0;
        return dup;
    }

    void realize(const State&, Stage) const;

    virtual SubsystemRep* cloneSubsystemRep() const = 0;
    virtual void endConstruction() { }

    virtual void realizeTopology(State& s) const { 
    }
    virtual void realizeModel(State& s) const { 
    }
    virtual void realizeInstance(const State& s) const { 
    }
    virtual void realizeTime(const State& s) const { 
    }
    virtual void realizePosition(const State& s) const { 
    }
    virtual void realizeVelocity(const State& s) const { 
    }
    virtual void realizeDynamics(const State& s) const { 
    }
    virtual void realizeAcceleration(const State& s) const { 
    }
    virtual void realizeReport(const State& s) const { 
    }

	bool isInSystem() const {return mySystem != 0;}
	bool isInSameSystem(const System& sys) const {
		return mySystem && mySystem==&sys;
	}
	const System& getSystem() const {
		assert(isInSystem()); // TODO
		return *mySystem;
	}
	System& updSystem() {
		assert(isInSystem()); // TODO
		return *mySystem;
	}
	void setSystem(System& sys, int ix) {
		assert(!isInSystem()); // TODO
		assert(ix >= 0);
		mySystem = &sys;
		mySubsystemIndex = ix;
	}
	int getMySubsystemIndex() const {
		SimTK_ASSERT(isInSystem(), "getMySubsystemIndex()");
		return mySubsystemIndex;
	}

    void setMyHandle(Subsystem& h) {myHandle = &h;}
    const Subsystem& getMyHandle() const {assert(myHandle); return *myHandle;}
    void clearMyHandle() {myHandle=0;}

protected:
    String  subsystemName;
    String  subsystemVersion;
	System* mySystem;		    // the System to which this Subsystem belongs
	int     mySubsystemIndex;	// Subsystem # within System

private:
    friend class Subsystem;
    Subsystem* myHandle;	// the owner handle of this rep
};

// This is the concrete Rep class implementing DefaultSystemSubsystem, i.e.,
// subsystem 0.
class DefaultSystemSubsystemRep : public SubsystemRep {
public:
    DefaultSystemSubsystemRep() : SubsystemRep() { }
    DefaultSystemSubsystemRep(const String& sysName, const String& sysVersion)
        : SubsystemRep(sysName, sysVersion)
    {
    }
    // The one pure virtual.
    DefaultSystemSubsystemRep* cloneSubsystemRep() const {
        return new DefaultSystemSubsystemRep(*this);
    }

    SimTK_DOWNCAST(DefaultSystemSubsystemRep, SubsystemRep);
};

} // namespace SimTK

#endif // SimTK_SIMBODY_SYSTEM_REP_H_
