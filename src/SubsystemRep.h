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

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/System.h"
#include "simbody/internal/Subsystem.h"

namespace SimTK {

class SubsystemRep {
public:
	SubsystemRep() 
      : subsystemName("<NONAME>"), subsystemVersion("0.0.0"),
        mySystem(0), mySubsystemId(InvalidSubsystemId), myHandle(0), 
        subsystemTopologyRealized(false) 
    { 
    }

	SubsystemRep(const String& name, const String& version) 
      : subsystemName(name), subsystemVersion(version),
        mySystem(0), mySubsystemId(InvalidSubsystemId), myHandle(0),
        subsystemTopologyRealized(false)
    { 
    }

    virtual ~SubsystemRep() { 
    }

    const String& getName()    const {return subsystemName;}
    const String& getVersion() const {return subsystemVersion;}

    // Use these to allocate state variables and cache entries that are owned
    // by this Subsystem.

    // qdot, qdotdot also allocated in cache
    int allocateQ(State& s, const Vector& qInit) const {
        return s.allocateQ(getMySubsystemId(), qInit);
    }

    // udot is also allocated in the cache
    int allocateU(State& s, const Vector& uInit) const {
        return s.allocateU(getMySubsystemId(), uInit);
    }

    // zdot is also allocated in the cache
    int allocateZ(State& s, const Vector& zInit) const {
        return s.allocateZ(getMySubsystemId(), zInit);
    }

    // qerr, uerr, udoterr are all cache entries, not variables
    int allocateQErr(State& s, int nqerr) const {
        return s.allocateQErr(getMySubsystemId(), nqerr);
    }

    int allocateUErr(State& s, int nuerr) const {
        return s.allocateUErr(getMySubsystemId(), nuerr);
    }

    int allocateUDotErr(State& s, int nudoterr) const {
        return s.allocateUDotErr(getMySubsystemId(), nudoterr);
    }

    int allocateDiscreteVariable(State& s, Stage g, AbstractValue* v) const {
        return s.allocateDiscreteVariable(getMySubsystemId(), g, v);
    }

    int allocateCacheEntry(State& s, Stage g, AbstractValue* v) const {
        return s.allocateCacheEntry(getMySubsystemId(), g, v);
    }

    void advanceToStage(const State& s, Stage g) const {
        s.advanceSubsystemToStage(getMySubsystemId(), g);
    }

    // These pull out the State entries which belong exclusively to
    // this Subsystem. These variables and cache entries are available
    // as soon as this subsystem is at stage Model.
    Stage getStage(const State& s) const {
        return s.getSubsystemStage(getMySubsystemId());
    }
    const AbstractValue& getDiscreteVariable(const State& s, int index) const {
        return s.getDiscreteVariable(getMySubsystemId(), index);
    }
    // State is *not* mutable here -- must have write access to change state variables.
    AbstractValue& updDiscreteVariable(State& s, int index) const {
        return s.updDiscreteVariable(getMySubsystemId(), index);
    }
    const AbstractValue& getCacheEntry(const State& s, int index) const {
        return s.getCacheEntry(getMySubsystemId(), index);
    }
    // State is mutable here.
    AbstractValue& updCacheEntry(const State& s, int index) const {
        return s.updCacheEntry(getMySubsystemId(), index);
    }

    // These return views on State shared global resources. The views
    // are private to this subsystem, but the global resources themselves
    // are not allocated until the *System* advances to stage Model.
    // Note that there is no subsystem equivalent of the State's "Y"
    // vector because in general a subsystem's state variables will
    // not be contiguous. However, a subsystem's Q's, U's, and Z's
    // will all be contiguous within those arrays.

    const Vector& getQ(const State& s) const {return s.getQ(getMySubsystemId());}
    const Vector& getU(const State& s) const {return s.getU(getMySubsystemId());}
    const Vector& getZ(const State& s) const {return s.getZ(getMySubsystemId());}

    // Not mutable: must have a writable state.
    Vector& updQ(State& s) const {return s.updQ(getMySubsystemId());}
    Vector& updU(State& s) const {return s.updU(getMySubsystemId());}
    Vector& updZ(State& s) const {return s.updZ(getMySubsystemId());}

    const Vector& getQDot   (const State& s) const {return s.getQDot(getMySubsystemId());}
    const Vector& getUDot   (const State& s) const {return s.getUDot(getMySubsystemId());}
    const Vector& getZDot   (const State& s) const {return s.getZDot(getMySubsystemId());}
    const Vector& getQDotDot(const State& s) const {return s.getQDotDot(getMySubsystemId());}

    // These are mutable
    Vector& updQDot   (const State& s) const {return s.updQDot(getMySubsystemId());}
    Vector& updUDot   (const State& s) const {return s.updUDot(getMySubsystemId());}
    Vector& updZDot   (const State& s) const {return s.updZDot(getMySubsystemId());}
    Vector& updQDotDot(const State& s) const {return s.updQDotDot(getMySubsystemId());}

    const Vector& getQErr(const State& s) const {return s.getQErr(getMySubsystemId());}
    const Vector& getUErr(const State& s) const {return s.getUErr(getMySubsystemId());}
    const Vector& getUDotErr(const State& s) const {return s.getUDotErr(getMySubsystemId());}

    // These are mutable
    Vector& updQErr(const State& s) const {return s.updQErr(getMySubsystemId());}
    Vector& updUErr(const State& s) const {return s.updUErr(getMySubsystemId());}
    Vector& updUDotErr(const State& s) const {return s.updUDotErr(getMySubsystemId());}

    // Dimensions. These are valid at Stage::Model while access to the various
    // arrays may have stricter requirements. Hence it is better to use these
    // routines than to get a reference to a Vector and ask for its size().

    int getQStart(const State& s) const {return s.getQStart(getMySubsystemId());}
    int getNQ(const State& s)     const {return s.getNQ(getMySubsystemId());}

    int getUStart(const State& s) const {return s.getUStart(getMySubsystemId());}
    int getNU(const State& s)     const {return s.getNU(getMySubsystemId());}

    int getZStart(const State& s) const {return s.getZStart(getMySubsystemId());}
    int getNZ(const State& s)     const {return s.getNZ(getMySubsystemId());}

    int getQErrStart(const State& s) const {return s.getQErrStart(getMySubsystemId());}
    int getNQErr(const State& s)     const {return s.getNQErr(getMySubsystemId());}

    int getUErrStart(const State& s) const {return s.getUErrStart(getMySubsystemId());}
    int getNUErr(const State& s)     const {return s.getNUErr(getMySubsystemId());}

    int getUDotErrStart(const State& s) const {return s.getUDotErrStart(getMySubsystemId());}
    int getNUDotErr(const State& s)     const {return s.getNUDotErr(getMySubsystemId());}

    SubsystemRep* clone() const {
		assert(!isInSystem()); // TODO
        SubsystemRep* dup = cloneSubsystemRep();
        dup->myHandle = 0;
        return dup;
    }

    //void realize(const State&, Stage) const;

    // Default implementation does nothing but blow up if the State hasn't been realized to at
    // least the requested stage.
    virtual void calcDecorativeGeometryAndAppend
       (const State& s, Stage stage, Array<DecorativeGeometry>& geom) const
    {
        assert(stage==Stage::Topology || getStage(s) >= stage);
    }
        

    virtual SubsystemRep* cloneSubsystemRep() const = 0;

    // These routines wrap the virtual realizeSubsystem...Impl() methods to ensure
    // good behavior such as checking that stage requirements are met and
    // updating the stage at the end. Note that these will do nothing if
    // the Subsystem stage is already at or greater than the indicated stage.
    void realizeSubsystemTopology(State& s) const;
    void realizeSubsystemModel(State& s) const;
    void realizeSubsystemInstance(const State& s) const;
    void realizeSubsystemTime(const State& s) const;
    void realizeSubsystemPosition(const State& s) const;
    void realizeSubsystemVelocity(const State& s) const;
    void realizeSubsystemDynamics(const State& s) const;
    void realizeSubsystemAcceleration(const State& s) const;
    void realizeSubsystemReport(const State& s) const;

    // These virtual methods should be overridden in concrete Subsystems as
    // necessary. They should never be called directly; instead call the
    // wrapper routines above. The wrappers will call these only when
    // the current stage for the
    // Subsystem is the one just prior to the stage being realized. For example,
    // realizeSubsystemVelocityImpl() is called by realizeSubsystemVelocity()
    // only when the passed-in State shows this subsystem's stage to be 
    // exactly Stage::Position.
    //
    // The default implementations provided here do nothing. That means the
    // wrappers will simply check that the current stage is correct and
    // advance it if necessary.

    virtual void realizeSubsystemTopologyImpl(State& s) const { 
    }
    virtual void realizeSubsystemModelImpl(State& s) const { 
    }
    virtual void realizeSubsystemInstanceImpl(const State& s) const { 
    }
    virtual void realizeSubsystemTimeImpl(const State& s) const { 
    }
    virtual void realizeSubsystemPositionImpl(const State& s) const { 
    }
    virtual void realizeSubsystemVelocityImpl(const State& s) const { 
    }
    virtual void realizeSubsystemDynamicsImpl(const State& s) const { 
    }
    virtual void realizeSubsystemAccelerationImpl(const State& s) const { 
    }
    virtual void realizeSubsystemReportImpl(const State& s) const { 
    }

    virtual void calcQUnitWeights(const State& s, Vector& weights) const {
        weights.resize(getNQ(s));
        weights = 1; // default says everyone's opinion is just as valid
    }
    virtual void calcUUnitWeights(const State& s, Vector& weights) const {
        weights.resize(getNU(s));
        weights = 1;
    }
    virtual void calcZUnitWeights(const State& s, Vector& weights) const {
        weights.resize(getNZ(s));
        weights = 1;
    }
    virtual void calcQErrUnitTolerances(const State& s, Vector& tolerances) const {
        tolerances.resize(getNQErr(s));
        tolerances = 1;
    }
    virtual void calcUErrUnitTolerances(const State& s, Vector& tolerances) const {
        tolerances.resize(getNUErr(s));
        tolerances = 1;
    }

	bool isInSystem() const {return mySystem != 0;}
	bool isInSameSystem(const Subsystem& otherSubsystem) const {
		return isInSystem() && otherSubsystem.isInSystem()
            && getSystem().isSameSystem(otherSubsystem.getSystem());
	}

	const System& getSystem() const {
		assert(isInSystem()); // TODO
		return *mySystem;
	}
	System& updSystem() {
		assert(isInSystem()); // TODO
		return *mySystem;
	}
	void setSystem(System& sys, SubsystemId id) {
		assert(!isInSystem()); // TODO
		assert(id.isValid());
		mySystem = &sys;
		mySubsystemId = id;
	}
	SubsystemId getMySubsystemId() const {
		SimTK_ASSERT(isInSystem(), "getMySubsystemId()");
		return mySubsystemId;
	}

    bool subsystemTopologyHasBeenRealized() const {
        return subsystemTopologyRealized;
    }

    // This must invalidate the topology of the containing System as well.
    void invalidateSubsystemTopologyCache();

    void setMyHandle(Subsystem& h) {myHandle = &h;}
    const Subsystem& getMyHandle() const {assert(myHandle); return *myHandle;}
    Subsystem& updMyHandle() {assert(myHandle); return *myHandle;}
    void clearMyHandle() {myHandle=0;}

protected:
    String      subsystemName;
    String      subsystemVersion;
	System*     mySystem;       // the System to which this Subsystem belongs
	SubsystemId mySubsystemId;  // Subsystem # within System

private:
    friend class Subsystem;
    Subsystem* myHandle;	// the owner handle of this rep

    mutable bool subsystemTopologyRealized;
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

#endif // SimTK_SUBSYSTEM_REP_H_
