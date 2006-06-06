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

    virtual void realizeConstruction(State& s) const { 
        advanceToStage(s, Stage::Built);
    }
    virtual void realizeModeling(State& s) const { 
        advanceToStage(s, Stage::Modeled);
    }
    virtual void realizeParameters(const State& s) const { 
        advanceToStage(s, Stage::Parametrized);
    }
    virtual void realizeTime(const State& s) const { 
        advanceToStage(s, Stage::Timed);
    }
    virtual void realizeConfiguration(const State& s) const { 
        advanceToStage(s, Stage::Configured);
    }
    virtual void realizeMotion(const State& s) const { 
        advanceToStage(s, Stage::Moving);
    }
    virtual void realizeDynamics(const State& s) const { 
        advanceToStage(s, Stage::Dynamics);
    }
    virtual void realizeReaction(const State& s) const { 
        advanceToStage(s, Stage::Reacting);
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

class SystemRep {
public:
    SystemRep() 
      : systemName("<NONAME>"), systemVersion("0.0.0"), subsystems(1), myHandle(0)
    {
        takeOverSubsystem(0, DefaultSystemSubsystem());
    }
    SystemRep(int nSubsystems, const String& name, const String& version) 
      : systemName(name), systemVersion(version), subsystems(nSubsystems), myHandle(0)
    {
        assert(nSubsystems >= 1);
        takeOverSubsystem(0, DefaultSystemSubsystem());
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

	// Take over ownership from the Subsystem handle, install
    // it into a particular subsystem slot, and return the
	// new handle. This is only allowed if (a) the subsystem
    // slot is empty, and (b) the supplied Subsytem
	// already has a rep, but is NOT part of some other System.
	Subsystem& takeOverSubsystem(int subsys, Subsystem& src) {
        assert(0 <= subsys && subsys < getNSubsystems());
        assert(!subsystems[subsys].hasRep());
		assert(src.hasRep() && !src.isInSystem()); // TODO
		Subsystem& s = subsystems[subsys];
		s.setRep(src.updRep());			 // reference the passed-in rep
		s.updRep().setMyHandle(s);	     // steal ownership
        s.updRep().setSystem(*myHandle, subsys);
		return s;
	}

    virtual SystemRep* cloneSystemRep() const = 0;


    virtual void realizeConstruction(State& s) const { 
        for (int i=0; i<getNSubsystems(); ++i)
            subsystems[i].realize(s, Stage::Built);
    }
    virtual void realizeModeling(State& s) const {
        for (int i=0; i<getNSubsystems(); ++i)
            subsystems[i].realize(s, Stage::Modeled);
    }
    virtual void realizeParameters(const State& s) const { 
        for (int i=0; i<getNSubsystems(); ++i)
            subsystems[i].realize(s, Stage::Parametrized);
    }
    virtual void realizeTime(const State& s) const { 
        for (int i=0; i<getNSubsystems(); ++i)
            subsystems[i].realize(s, Stage::Timed);
    }
    virtual void realizeConfiguration(const State& s) const { 
        for (int i=0; i<getNSubsystems(); ++i)
            subsystems[i].realize(s, Stage::Configured);
    }
    virtual void realizeMotion(const State& s) const { 
        for (int i=0; i<getNSubsystems(); ++i)
            subsystems[i].realize(s, Stage::Moving);
    }
    virtual void realizeDynamics(const State& s) const { 
        for (int i=0; i<getNSubsystems(); ++i)
            subsystems[i].realize(s, Stage::Dynamics);
    }
    virtual void realizeReaction(const State& s) const { 
        for (int i=0; i<getNSubsystems(); ++i)
            subsystems[i].realize(s, Stage::Reacting);
    }

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
