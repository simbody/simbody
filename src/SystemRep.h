#ifndef SimTK_SIMBODY_SYSTEM_REP_H_
#define SimTK_SIMBODY_SYSTEM_REP_H_

/* Copyright (c) 2006 Stanford University and Michael Sherman.
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

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/State.h"
#include "simbody/internal/System.h"


namespace SimTK {

class AnalyticGeometry;
class DecorativeGeometry;
class System;
class Subsystem;
class Study;

// An object of this class is stored as a construction-stage discrete variable
// in a System's State, one for each Subsystem. The appropriate object can
// then be found by each Subsystem, from which it can find the rest of
// its private state entries. The 0th object is reserved for the System
// itself.
struct SubsystemDescriptor {
    SubsystemDescriptor(int ix, const String& nm, const String& vers)
      : subsystemIndex(ix), name(nm), version(vers),
        constructionVarsIndex(-1), modelingVarsIndex(-1)
    { 
        SimTK_ASSERT1_ALWAYS(subsystemIndex >= 0,
            "SubsystemDescriptor(int): index was negative (%d)",
            subsystemIndex);
    }

    int     subsystemIndex;
    String  name;
    String  version;

    // These are DiscreteVariable indices within the enclosing State.
    // The first is expected to be Construction stage, the second is
    // Modeling stage. The construction stage variable (if present)
    // is typically there only for sanity checking the state.
    // The modeling variable contains real information, plus the
    // indices of all the remaining variables and cache entries
    // belonging to this subsystem.
    int     constructionVarsIndex;
    int     modelingVarsIndex;
};

std::ostream& 
operator<<(std::ostream& o, const SubsystemDescriptor&);

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

    Stage getStage(const State& s) const {
        return s.getSubsystemStage(getMySubsystemIndex());
    }

    SubsystemRep* clone() const {
		assert(!isInSystem()); // TODO
        SubsystemRep* dup = cloneSubsystemRep();
        dup->myHandle = 0;
        return dup;
    }

    virtual SubsystemRep* cloneSubsystemRep() const = 0;
    virtual void endConstruction() { }
    virtual void realizeConstruction(State&) const = 0;
    virtual void realizeModeling(State&) const = 0;

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

class SystemRep {
public:
    SystemRep() 
      : systemName("<NONAME>"), systemVersion("0.0.0"), subsystems(1), myHandle(0)
    {
    }
    SystemRep(int nSubsystems, const String& name, const String& version) 
      : systemName(name), systemVersion(version), subsystems(nSubsystems), myHandle(0)
    {
    }
    virtual ~SystemRep() {
        clearMyHandle();
        subsystems.clear();
    }

    const String& getName()    const {return systemName;}
    const String& getVersion() const {return systemVersion;}

    int              getNSubsystems()    const {return subsystems.size();}
    const Subsystem& getSubsystem(int i) const {return subsystems[i];}
    Subsystem&       updSubsystem(int i)       {return subsystems[i];}

    SystemRep* clone() const {
        SystemRep* dup = cloneSystemRep();
        dup->myHandle = 0;
        return dup;
    }

	// Take over ownership from the Subsystem handle, and return the
	// new handle. This is only allowed if the supplied Subsytem
	// already has a rep, but is NOT part of some other System.
	Subsystem& takeOverSubsystem(Subsystem& src) {
		assert(src.hasRep() && !src.isInSystem()); // TODO
		// Push an empty handle on the subsystem list
		subsystems.resize(subsystems.size()+1);
		Subsystem& s = subsystems.back();
		s.setRep(src.updRep());			 // reference the passed-in rep
		s.updRep().setMyHandle(s);	     // steal ownership
		return s;
	}

    virtual SystemRep* cloneSystemRep() const = 0;

    virtual void realizeConstruction (State& s)       const = 0;
    virtual void realizeModeling     (State& s)       const = 0;
    virtual void realizeParameters   (const State& s) const { }
    virtual void realizeTime         (const State& s) const { }
    virtual void realizeConfiguration(const State& s) const { }
    virtual void realizeMotion       (const State& s) const { }
    virtual void realizeDynamics     (const State& s) const { }
    virtual void realizeReaction     (const State& s) const { }

    void realize(const State& s, Stage g) const;

    void setMyHandle(System& h) {myHandle = &h;}
    void clearMyHandle() {myHandle=0;}

protected:
    String systemName;
    String systemVersion;
	StableArray<Subsystem> subsystems;

private:
    friend class System;
    System* myHandle;     // the owner of this rep
};


class StudyRep {
public:
    StudyRep(const System& sys)
      : myHandle(0), system(new System(sys))
    {
        system->realize(state, Stage::Built);
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

#endif // SimTK_SIMBODY_SYSTEM_REP_H_
