#ifndef SimTK_SYSTEM_REP_H_
#define SimTK_SYSTEM_REP_H_

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

#include "SubsystemRep.h"

namespace SimTK {

class SystemRep {
public:
    SystemRep() 
      : systemName("<NONAME>"), systemVersion("0.0.0"), 
        myHandle(0), systemTopologyRealized(false)
    {
    }
    SystemRep(const String& name, const String& version) 
      : systemName(name), systemVersion(version), 
        myHandle(0), systemTopologyRealized(false)
    {
    }
    virtual ~SystemRep() {
        clearMyHandle();
        subsystems.clear();
        invalidateSystemTopologyCache();
    }

    const String& getName()    const {return systemName;}
    const String& getVersion() const {return systemVersion;}

    // This is available after realizeTopology().
    const State& getDefaultState() const {
        SimTK_ASSERT_ALWAYS(systemTopologyHasBeenRealized(),
            "System::getDefaultState(): realizeTopology() must be called first.");
        return defaultState;
    }

    int              getNSubsystems()            const {return subsystems.size();}
    const Subsystem& getSubsystem(SubsystemId i) const {return subsystems[i];}
    Subsystem&       updSubsystem(SubsystemId i)       {return subsystems[i];}

    SystemRep* clone() const {
        SystemRep* dup = cloneSystemRep();
        dup->myHandle = 0;
        return dup;
    }

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

    virtual SystemRep* cloneSystemRep() const = 0;

    // These routines wrap the virtual realize...Impl() methods to ensure
    // good behavior such as checking that stage requirements are met and
    // updating the stage at the end. Note that these will do nothing if
    // the System stage is already at or greater than the indicated stage.
    const State& realizeTopology()           const;
    void realizeModel       (State &s)       const;
    void realizeInstance    (const State &s) const;
    void realizeTime        (const State &s) const;
    void realizePosition    (const State &s) const;
    void realizeVelocity    (const State &s) const;
    void realizeDynamics    (const State &s) const;
    void realizeAcceleration(const State &s) const;
    void realizeReport      (const State &s) const;

    // For a State that has already been realized to Stage::Model or higher,
    // this routine conveniently realizes the State up to the indicated Stage
    // (if needed), one Stage at a time, using the above routines. This will
    // throw an exception if the State hasn't already been realized to at
    // least Stage::Model.
    void realize(const State& s, Stage g) const;

    // Override these to change the evaluation order of the Subsystems.
    // The default is to evaluate them in increasing order of SubsystemId.
    // These methods should not be called directly; they are invoked by the
    // above wrapper methods. Note: the wrappers *will not* call these
    // routines if the system stage has already met the indicated stage level.
    // If fact these routines will be called only when the system stage
    // is at the level just prior to the one indicated here. For example,
    // realizeVelocityImpl() will be called only if the passed-in State
    // has been determined to have its system stage exactly Stage::Position.
    virtual void realizeTopologyImpl(State& s) const { 
        for (int i=0; i<getNSubsystems(); ++i)
            subsystems[i].getRep().realizeSubsystemTopology(s);
    }
    virtual void realizeModelImpl(State& s) const {
        for (int i=0; i<getNSubsystems(); ++i)
            subsystems[i].getRep().realizeSubsystemModel(s);
    }
    virtual void realizeInstanceImpl(const State& s) const { 
        for (int i=0; i<getNSubsystems(); ++i)
            subsystems[i].getRep().realizeSubsystemInstance(s);
    }
    virtual void realizeTimeImpl(const State& s) const { 
        for (int i=0; i<getNSubsystems(); ++i)
            subsystems[i].getRep().realizeSubsystemTime(s);
    }
    virtual void realizePositionImpl(const State& s) const { 
        for (int i=0; i<getNSubsystems(); ++i)
            subsystems[i].getRep().realizeSubsystemPosition(s);
    }
    virtual void realizeVelocityImpl(const State& s) const { 
        for (int i=0; i<getNSubsystems(); ++i)
            subsystems[i].getRep().realizeSubsystemVelocity(s);
    }
    virtual void realizeDynamicsImpl(const State& s) const { 
        for (int i=0; i<getNSubsystems(); ++i)
            subsystems[i].getRep().realizeSubsystemDynamics(s);
    }
    virtual void realizeAccelerationImpl(const State& s) const { 
        for (int i=0; i<getNSubsystems(); ++i)
            subsystems[i].getRep().realizeSubsystemAcceleration(s);
    }
    virtual void realizeReportImpl(const State& s) const { 
        for (int i=0; i<getNSubsystems(); ++i)
            subsystems[i].getRep().realizeSubsystemReport(s);
    }

    void calcDecorativeGeometryAndAppend(const State&, Stage, Array<DecorativeGeometry>&) const;

    void calcYUnitWeights(const State& s, Vector& weights) const {
        weights.resize(s.getNY());
        VectorView qwts = weights(s.getQStart(), s.getNQ());   // writable views
        VectorView uwts = weights(s.getUStart(), s.getNU());
        VectorView zwts = weights(s.getZStart(), s.getNZ());

        for (int i=0; i<getNSubsystems(); ++i) {
            const Subsystem& sub = subsystems[i];
            sub.calcQUnitWeights(s, qwts(s.getQStart(i), s.getNQ(i)));
            sub.calcUUnitWeights(s, uwts(s.getUStart(i), s.getNU(i)));
            sub.calcZUnitWeights(s, zwts(s.getZStart(i), s.getNZ(i)));
        }
    }

    void calcYErrUnitTolerances(const State& s, Vector& tolerances) const {
        tolerances.resize(s.getNYErr());
        VectorView qtols = tolerances(s.getQErrStart(), s.getNQErr()); // writable views
        VectorView utols = tolerances(s.getUErrStart(), s.getNUErr());

        for (int i=0; i<getNSubsystems(); ++i) {
            const Subsystem& sub = subsystems[i];
            sub.calcQErrUnitTolerances(s, qtols(s.getQErrStart(i), s.getNQErr(i)));
            sub.calcUErrUnitTolerances(s, utols(s.getUErrStart(i), s.getNUErr(i)));
        }
    }

    virtual Real calcTimescale(const State& s) const {
        SimTK_STAGECHECK_GE(s.getSystemStage(), Stage::Instance,
            "System::calcTimescale()");
        return 0.1; // TODO!!!
    }

    // Default treats all state variable identically. Should be asking the 
    // subsystems. TODO
    virtual Real calcYErrorNorm(const State& s, const Vector& y_err) const {
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

    void invalidateSystemTopologyCache() {
        systemTopologyRealized = false;
        defaultState.clear();
    }

protected:
    String systemName;
    String systemVersion;
	StableArray<Subsystem> subsystems;

private:
    friend class System;
    System* myHandle;     // the owner of this rep

        // TOPOLOGY STAGE CACHE //

    // This should only be true when *all* subsystems have successfully
    // completed realizeTopology(). Anything which invalidates topology for
    // one of the contained subsystem must invalidate topology for the system
    // as a whole also.
    mutable bool systemTopologyRealized;

    // This is only meaningful if topologyRealized==true.
    mutable State defaultState;
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

#endif // SimTK_SYSTEM_REP_H_
