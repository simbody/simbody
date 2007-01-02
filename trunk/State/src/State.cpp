/* Portions copyright (c) 2005-6 Stanford University and Michael Sherman.
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

#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"
#include "SimTKcommon/State.h"

#include <cassert>
#include <ostream>

namespace SimTK {

class DiscreteVariableRep {
public:
    DiscreteVariableRep() : value(0), myHandle(0) { }
    DiscreteVariableRep(const DiscreteVariableRep& src) {
        assert(src.isValid());
        stage = src.stage;
        value = src.getValue().clone();
        myHandle = 0;
    }
    DiscreteVariableRep& operator=(const DiscreteVariableRep& src) {
        if (&src != this) {
            assert(src.isValid());
            assert(stage == src.stage);
            if (value) *value = src.getValue();
            else value = src.getValue().clone();
        }
        return *this;
    }
    ~DiscreteVariableRep() {
        clear(); myHandle = 0;
    }

    DiscreteVariableRep(Stage g, AbstractValue* vp) 
      : stage(g), value(vp), myHandle(0)
    {
        assert(g == Stage::Topology || Stage::isInRuntimeRange(g));
        assert(vp);
    }

    DiscreteVariableRep* clone() const {return new DiscreteVariableRep(*this);}

    bool isValid() const {
        return (stage==Stage::Topology || Stage::isInRuntimeRange(stage)) && value && myHandle; 
    }

    Stage getStage() const {return stage;}

    const AbstractValue& getValue() const {assert(value); return *value;}
    AbstractValue&       updValue()       {assert(value); return *value;}

    void                    setMyHandle(DiscreteVariable& dv) {myHandle = &dv;}
    const DiscreteVariable& getMyHandle() const   {assert(myHandle); return *myHandle;}
    DiscreteVariable&       updMyHandle()         {assert(myHandle); return *myHandle;}
private:
    void clear() {stage=Stage::Invalid; delete value;}

    Stage stage;
    AbstractValue* value;

    DiscreteVariable* myHandle;
};

// This internal utility class is used to capture all the information needed for
// a single subsystem within the StateRep.
class PerSubsystemInfo {
public:
    PerSubsystemInfo()      {initialize();}
    PerSubsystemInfo(const String& n, const String& v) 
      : name(n), version(v) {initialize();}

    ~PerSubsystemInfo() {   // default destructor
    }

    // Copy constructor copies all variables but cache only through
    // modeled stage. Note that this must be done in conjunction with
    // copying the whole state or our global resource indices will
    // be nonsense.
    PerSubsystemInfo(const PerSubsystemInfo& src) {
        initialize();
        copyFrom(src, Stage::Model);
    }

    PerSubsystemInfo& operator=(const PerSubsystemInfo& src) {
        if (&src != this) {
            copyFrom(src, Stage::Model);
        }
        return *this;
    }

    // Back up to the stage just before g if this subsystem thinks
    // it is already at g or beyond. Note that we may be backing up
    // over many stages here. Careful: invalidating the stage
    // for a subsystem must also invalidate the same stage for all
    // the other subsystems and the system as a whole but we don't
    // take care of that here. Also, you can't invalidate Stage::Empty.
    void invalidateStageJustThisSubsystem(Stage g) {
        assert(g > Stage::Empty);
        restoreToStage(g.prev());
    }

    // Advance from stage g-1 to stage g. This is called at the end
    // of realize(g). You can't use this to "advance" to Stage::Empty.
    // It is a fatal error if the current stage isn't g-1.
    void advanceToStage(Stage g) {
        assert(g > Stage::Empty);
        assert(currentStage == g.prev());
        // Record data needed to back up to this stage later.
        if (g == Stage::Topology)
            nDiscreteWhenBuilt = discrete.size();
        cacheSize[g] = cache.size();
        currentStage = g;
    }


    void clearReferencesToStateGlobals() {
        qstart=ustart=zstart=qerrstart=uerrstart=udoterrstart = -1;
        q.clear(); u.clear(); z.clear();
        qdot.clear(); udot.clear(); zdot.clear(); qdotdot.clear();
        qerr.clear(); uerr.clear(); udoterr.clear();
    }

    String name;
    String version;

    // These accumulate default values for this subsystem's use of shared
    // global state variables. After the System is advanced to Stage::Model,
    // the state will allocate those globals and copy these initial
    // values into them. The lengths of these Vectors define the 
    // needs of this Subsystem.
    Vector qInit, uInit, zInit;

    // For constraints we need just lengths.
    int nqerr, nuerr, nudoterr;

    // These are our own private views into partitions of the global
    // state and cache entries of the same names. These are valid
    // only after the *System* stage is raised to Model, and they
    // are invalidated whenever the System's Model stage is invalidated.

    // The State will assign contiguous blocks to this subsystem. The
    // starting indices are filled in here at the time the views are built.
    int qstart, ustart, zstart, qerrstart, uerrstart, udoterrstart;
    Vector q, u, z;
    mutable Vector qdot, udot, zdot, qdotdot;
    mutable Vector qerr, uerr, udoterr;

    // Each of these discrete variable & cache entries has its own Stage.
    StableArray<DiscreteVariable> discrete;
    int nDiscreteWhenBuilt;

    mutable StableArray<CacheEntry> cache;
    int cacheSize[Stage::NValid];   // cache size at the end of each stage

    mutable Stage currentStage;

private:
    // This is for use in constructors. Everything that doesn't have
    // a reasonable default constructor gets set here.
    void initialize() {
        nqerr = nuerr = nudoterr= 0;
        qstart=ustart=zstart=qerrstart=uerrstart=udoterrstart = -1;
        nDiscreteWhenBuilt = -1;
        for (int i=0; i < Stage::NValid; ++i)
            cacheSize[i] = -1;
        cacheSize[Stage::Empty] = 0;
        currentStage = Stage::Empty;
    }


    // Set all the allocation sizes to zero, leaving name and version alone.
    // The stage is set back to "Empty".
    void restoreToEmptyStage() {
        if (currentStage > Stage::Empty) {
            clearReferencesToStateGlobals();
            discrete.clear(); cache.clear();
            qInit.clear(); uInit.clear(); zInit.clear();
            initialize();
        }
    }

    // Put this Subsystem back to the way it was just after realize(Topology).
    // That means: no shared global state vars, discrete vars and cache
    // only as they were after realize(Topology).
    // The stage is set back to "Topology". Nothing happens if the stage
    // is already at Topology or below.
    void restoreToTopologyStage() {
        if (currentStage <= Stage::Topology)
            return;
        clearReferencesToStateGlobals();
        discrete.resize(nDiscreteWhenBuilt);
        restoreCacheToStage(Stage::Topology);
        qInit.clear(); uInit.clear(); zInit.clear();
        nqerr = nuerr = nudoterr = 0;
        currentStage = Stage::Topology;
    }

    // Restore this subsystem to the way it last was at realize(Stage).
    void restoreToStage(Stage g) {
        if (g==Stage::Empty)    {restoreToEmptyStage();    return;}
        if (g==Stage::Topology) {restoreToTopologyStage(); return;}
        if (currentStage <= g) return;

        // State variables remain unchanged since they are all allocated
        // after realize(Model). Cache gets shrunk to the length it
        // had after realize(g) in case some entries were late additions.
        restoreCacheToStage(g);
        currentStage = g;
    }

    // Don't call this unless you know the current stage is > g.
    // We shrink the subsystem's cache to the size it had the last
    // time realize(g) was performed, and we forget about all the
    // later sizes we used to know about. CurrentStage is not changed,
    // so the subsystem will be a mess after this call!
    void restoreCacheToStage(Stage g) {
        assert(currentStage > g);
        assert(cacheSize[g] >= 0);
        cache.resize(cacheSize[g]);
        for (int i=((int)g+1); i<Stage::NValid; ++i)
            cacheSize[i] = -1;
    }


    // Utility which makes "this" a copy of the source subsystem exactly as it
    // was after being realized to stage maxStage. If maxStage >= Model then
    // all the subsystem-private state variables will be copied, but only
    // cached computations up through maxStage come through. We clear
    // our references to global variables regardless -- those will have to
    // be repaired at the System (State global) level.
    void copyFrom(const PerSubsystemInfo& src, Stage maxStage) {
        const Stage targetStage = std::min<Stage>(src.currentStage, maxStage);
        name     = src.name;
        version  = src.version;

        if (targetStage < Stage::Topology) {
            restoreToEmptyStage();  // don't copy anything
            return;
        }

        // At "Topology" stage we need to copy all the private state
        // variables that were present at realize(Topology) (those will
        // mostly be modeling variables). There can't
        // be any global ones since those aren't allocated until
        // Model stage. We also need to copy any cache results
        // that were available after realize(Topology).
        if (targetStage == Stage::Topology) {
            restoreToEmptyStage();
            discrete.resize(src.nDiscreteWhenBuilt);
            for (int i=0; i<src.nDiscreteWhenBuilt; ++i)
                discrete[i] = src.discrete[i];
            nDiscreteWhenBuilt = discrete.size();
            // don't copy any global shared resources

            cache.resize(src.cacheSize[Stage::Topology]);
            for (int i=0; i<(int)cache.size(); ++i)
                cache[i] = src.cache[i];
            cacheSize[Stage::Topology] = cache.size();
            currentStage = Stage::Topology;
            return;
        }

        // This is the general case where Stage > Topology.
        clearReferencesToStateGlobals();
        qInit = src.qInit;
        uInit = src.uInit;
        zInit = src.zInit;
        nqerr = src.nqerr;
        nuerr = src.nuerr;
        nudoterr = src.nudoterr;

        // Copy *all* state variables since no more can be allocated
        // after Model stage.
        discrete = src.discrete;

        // Copy only the cache as it was at the end of targetStage since
        // more might have been added later.
        cache.resize(src.cacheSize[targetStage]);
        for (int i=0; i<(int)cache.size(); ++i)
            cache[i] = src.cache[i];
        for (int i=0; i<=targetStage; ++i)
            cacheSize[i] = src.cacheSize[i];
        for (int i=targetStage+1; i<Stage::NValid; ++i)
            cacheSize[i] = -1;
        currentStage = targetStage;
    }
};

class StateRep {
public:
    StateRep() 
      : t(CNT<Real>::getNaN()), systemStage(Stage::Empty), 
        myHandle(0) 
    { 
    }

    ~StateRep() {   // default destructor
    }

    const Stage& getSystemStage() const {
        return systemStage;
    }
    Stage& updSystemStage() const {
        return systemStage; // mutable
    }

    const PerSubsystemInfo& getSubsystem(int subsystem) const {
        SimTK_INDEXCHECK(0, subsystem, (int)subsystems.size(), "StateRep::getSubsystem()");
        return subsystems[subsystem];
    }

    PerSubsystemInfo& updSubsystem(int subsystem) {
        SimTK_INDEXCHECK(0, subsystem, (int)subsystems.size(), "StateRep::updSubsystem()");
        return subsystems[subsystem];
    }

    const Stage& getSubsystemStage(int subsystem) const {
        return subsystems[subsystem].currentStage;
    }
    Stage& updSubsystemStage(int subsystem) const {
        return subsystems[subsystem].currentStage; // mutable
    }

    // We'll do the copy constructor and assignment explicitly here
    // to get tight control over what's allowed, and to make sure
    // we don't copy the handle pointer.
    StateRep(const StateRep& src) : myHandle(0), systemStage(Stage::Empty) {
        subsystems = src.subsystems;
        if (src.systemStage >= Stage::Topology) {
            advanceSystemToStage(Stage::Topology);
            if (src.systemStage >= Stage::Model) {
                advanceSystemToStage(Stage::Model);
                t = src.t;
                // careful -- don't allow reallocation
                y = src.y;
            }
        }
    }

    StateRep& operator=(const StateRep& src) {
        if (&src == this) return *this;
        invalidateJustSystemStage(Stage::Topology);
        for (int i=0; i<(int)subsystems.size(); ++i)
            subsystems[i].invalidateStageJustThisSubsystem(Stage::Topology);
        subsystems = src.subsystems;
        if (src.systemStage >= Stage::Topology) {
            advanceSystemToStage(Stage::Topology);
            if (src.systemStage >= Stage::Model) {
                advanceSystemToStage(Stage::Model);
                t = src.t;
                y = src.y;
            }
        }
        // don't mess with the handle pointer!
        return *this;
    }

    // Copies all the variables but not the cache.
    StateRep* clone() const {return new StateRep(*this);}

    // Back up the System stage just before g if it thinks
    // it is already at g or beyond. Note that we may be backing up
    // over many stages here. Careful: invalidating the stage
    // for the system must also invalidate the same stage for all
    // the subsystems (because we trash the shared resource pool
    // here if we back up earlier than Stage::Model) but we don't
    // take care of that here. Also, you can't invalidate Stage::Empty.
    void invalidateJustSystemStage(Stage g) {
        assert(g > Stage::Empty);
        if (systemStage >= g) {
            if (systemStage >= Stage::Model && Stage::Model >= g) {
                // We are "unmodeling" this State. Trash all the global
                // shared states & corresponding cache entries.

                // First make sure no subsystem is looking at the
                // global shared state any more.
                for (int i=0; i < (int)subsystems.size(); ++i)
                    subsystems[i].clearReferencesToStateGlobals();

                t = CNT<Real>::getNaN();
                // Nuke all the global views.
                q.clear(); u.clear(); z.clear();
                qdot.clear(); udot.clear(); zdot.clear();
                qerr.clear(); uerr.clear();
                // Nuke the actual data.
                y.unlockShape();       y.clear(); 
                ydot.unlockShape();    ydot.clear(); 
                qdotdot.unlockShape(); qdotdot.clear();
                yerr.unlockShape();    yerr.clear();
                udoterr.unlockShape(); udoterr.clear();
            }
            systemStage = g.prev();
        }
    }

    // Advance the System stage from g-1 to g. It is a fatal error if
    // we're not already at g-1, and you can't advance to Stage::Empty.
    // Also, you can't advance the system to g unless ALL subsystems have
    // already gotten there.
    void advanceSystemToStage(Stage g) {
        assert(g > Stage::Empty);
        assert(systemStage == g.prev());
        assert(allSubsystemsAtLeastAtStage(g));

        if (g == Stage::Model) {
            // We know the shared state pool sizes now. Allocate the
            // states and matching shared cache pools.
            int nq=0, nu=0, nz=0, nqerr=0, nuerr=0, nudoterr=0;
            for (int i=0; i<(int)subsystems.size(); ++i) {
                nq += subsystems[i].qInit.size();
                nu += subsystems[i].uInit.size();
                nz += subsystems[i].zInit.size();
                nqerr    += subsystems[i].nqerr;
                nuerr    += subsystems[i].nuerr;
                nudoterr += subsystems[i].nudoterr;
            }
            // Allocate the actual shared state variables & cache 
            // entries and make sure no one can accidentally change the size.
            y.resize(nq+nu+nz);             y.lockShape();
            ydot.resize(nq+nu+nz);          ydot.lockShape();
            qdotdot.resize(nq);             qdotdot.lockShape();
            yerr.resize(nqerr+nuerr);       yerr.lockShape();
            udoterr.resize(nudoterr);       udoterr.lockShape();

            // Allocate subviews of the shared state & cache entries.
            q.viewAssign(y(0,nq));
            u.viewAssign(y(nq,nu));
            z.viewAssign(y(nq+nu,nz));

            qdot.viewAssign(ydot(0,     nq));
            udot.viewAssign(ydot(nq,    nu));
            zdot.viewAssign(ydot(nq+nu, nz));

            qerr.viewAssign(yerr(0,     nqerr));
            uerr.viewAssign(yerr(nqerr, nuerr));

            // Now partition the global resources among the subsystems.
            int nxtq=0, nxtu=0, nxtz=0, nxtqerr=0, nxtuerr=0, nxtudoterr=0;
            for (int i=0; i<(int)subsystems.size(); ++i) {
                PerSubsystemInfo& ss = subsystems[i];
                const int nq=ss.qInit.size(), nu=ss.uInit.size(), nz=ss.zInit.size();

                // Assign the starting indices.
                ss.qstart=nxtq; ss.ustart=nxtu; ss.zstart=nxtz;
                ss.qerrstart=nxtqerr; ss.uerrstart=nxtuerr; ss.udoterrstart=nxtudoterr;

                // Build the views.
                ss.q.viewAssign(q(nxtq, nq)); ss.q = ss.qInit;
                ss.qdot.viewAssign(qdot(nxtq, nq));
                ss.qdotdot.viewAssign(qdotdot(nxtq, nq));
                ss.u.viewAssign(u(nxtu, nu)); ss.u = ss.uInit;
                ss.udot.viewAssign(udot(nxtu, nu));
                ss.z.viewAssign(z(nxtz, nz)); ss.z = ss.zInit;
                ss.zdot.viewAssign(zdot(nxtz, nz));

                ss.qerr.viewAssign(qerr(nxtqerr, ss.nqerr));
                ss.uerr.viewAssign(uerr(nxtuerr, ss.nuerr));
                ss.udoterr.viewAssign(udoterr(nxtudoterr, ss.nudoterr));

                // Consume the slots.
                nxtq += nq; nxtu += nu; nxtz += nz;
                nxtqerr += ss.nqerr; nxtuerr += ss.nuerr; nxtudoterr += ss.nudoterr;
            }
        }

        systemStage = g;
    }

    void         setMyHandle(State& s) {myHandle = &s;}
    const State& getMyHandle() const   {assert(myHandle); return *myHandle;}
    State&       updMyHandle()         {assert(myHandle); return *myHandle;}
private:
    friend class State;

        // Shared global resource State variables //

    Real            t; // Stage::Time (time)
    Vector          y; // All the continuous state variables taken together {q,u,z}

        // These are views into y.
    Vector          q; // Stage::Position continuous variables
    Vector          u; // Stage::Velocity continuous variables
    Vector          z; // Stage::Dynamics continuous variables


        // Shared global resource Cache entries //
    mutable Stage   systemStage;

        // DIFFERENTIAL EQUATIONS
    mutable Vector  ydot; // All the state derivatives taken together (qdot,udot,zdot)

        // These are views into ydot.
    mutable Vector  qdot;       // Stage::Velocity
    mutable Vector  udot;       // Stage::Acceleration
    mutable Vector  zdot;       // Stage::Acceleration

        // This is an independent cache entry.
    mutable Vector  qdotdot;    // Stage::Acceleration

        // ALGEBRAIC EQUATIONS
    mutable Vector  yerr;       // All constraint errors taken together (qerr,uerr)
    mutable Vector  udoterr;    // Stage::Acceleration (Index 1 constraints)

        // These are views into yerr.
    mutable Vector  qerr;       // Stage::Position (Index 3 constraints)
    mutable Vector  uerr;       // Stage::Velocity (Index 2 constraints)

        // Subsystem support //

    // Subsystem 0 (always present) is for the System as a whole. Its name
    // and version are the System name and version.
    std::vector<PerSubsystemInfo> subsystems;

    // Return true only if all subsystems have been realized to at least Stage g.
    bool allSubsystemsAtLeastAtStage(Stage g) const {
        for (int i=0; i < (int)subsystems.size(); ++i)
            if (subsystems[i].currentStage < g)
                return false;
        return true;
    }

private:
    State* myHandle;
};

    // DISCRETE VARIABLE

DiscreteVariable::~DiscreteVariable() {
    delete rep; rep=0;
}

DiscreteVariable::DiscreteVariable(Stage g, AbstractValue* v) {
    rep = new DiscreteVariableRep(g, v);
    rep->setMyHandle(*this);
}

// copy constructor
DiscreteVariable::DiscreteVariable(const DiscreteVariable& src) : rep(0) {
    if (src.rep) {
        rep = src.rep->clone();
        rep->setMyHandle(*this);
    }
}

// copy assignment
// If "this" has a rep already then the stages must match and the
// value types must match. Otherwise, we just duplicate the src
// and treat this as the definition of "this".
DiscreteVariable& 
DiscreteVariable::operator=(const DiscreteVariable& src) {
    if (&src == this) return *this;

    if (rep) { // "this" has already been defined
        assert(src.rep && src.getStage() == getStage());
        updValue() = src.getValue();
    } else { // "this" is not defined yet
        if (src.rep) {
            // defining
            rep = src.rep->clone();
            rep->setMyHandle(*this);
        }
    }
    return *this;
}

Stage DiscreteVariable::getStage() const {
    return rep->getStage();
}

const AbstractValue& 
DiscreteVariable::getValue() const {
    return rep->getValue();
}

AbstractValue&
DiscreteVariable::updValue() {
    return rep->updValue();
}



    // STATE
State::State() : rep(new StateRep()) {
    rep->setMyHandle(*this);
}

State::~State() {
    delete rep; rep=0;
}

// copy constructor
State::State(const State& src) 
  : rep(src.rep->clone()) {
    rep->setMyHandle(*this);
}

// copy assignment
State& State::operator=(const State& src) {
    if (&src == this) return *this;
    if (!rep) {
        // we're defining this state here (if src is not empty)
        if (src.rep) {
            rep = src.rep->clone();
            rep->setMyHandle(*this);
        }
        return *this;
    }

    // Assignment or redefinition
    if (src.rep) *rep = *src.rep;
    else {delete rep; rep=0;}
    return *this;
}

void State::setNSubsystems(int i) {
    assert(i >= 0);
    updRep().subsystems.clear();
    updRep().subsystems.resize(i);
}

void State::initializeSubsystem(int i, const String& name, const String& version) {
    updRep().updSubsystem(i).name = name;
    updRep().updSubsystem(i).version = version;
}


int State::addSubsystem(const String& name, const String& version) {
    rep->subsystems.push_back(
        PerSubsystemInfo(name,version));
    return (int)rep->subsystems.size() - 1;
}

int State::getNSubsystems() const {return (int)rep->subsystems.size();}

const String& State::getSubsystemName(int subsys) const {
    return rep->subsystems[subsys].name;
}
const String& State::getSubsystemVersion(int subsys) const {
    return rep->subsystems[subsys].version;
}

const Stage& State::getSystemStage() const {
    return rep->getSystemStage();
}

const Stage& State::getSubsystemStage(int subsys) const {
    SimTK_ASSERT(rep, "State::getStage(): no rep"); // can't happen(?)
    return rep->getSubsystemStage(subsys);
}

// Make sure the stage is no higher than g-1 for *any* subsystem and
// hence for the system stage also. TODO: this should be more selective.
void State::invalidateAll(Stage g) const {
    SimTK_ASSERT(rep, "State::invalidateAll(): no rep");

    rep->invalidateJustSystemStage(g);
    for (int i=0; i<(int)rep->subsystems.size(); ++i)
        rep->subsystems[i].invalidateStageJustThisSubsystem(g);
}

// Move the stage for a particular subsystem from g-1 to g. No other subsystems
// are affected, nor the global system stage.
void State::advanceSubsystemToStage(int subsys, Stage g) const {
    SimTK_ASSERT(rep, "State::advanceSubsystemToStage(): no rep");

    rep->subsystems[subsys].advanceToStage(g);
    // We don't automatically advance the System stage even if this brings
    // ALL the subsystems up to stage g.
}

// Move the system stage from g-1 to g. Don't call this until ALL 
// subsystem have been advanced to at least stage g.
void State::advanceSystemToStage(Stage g) const {
    SimTK_ASSERT(rep, "State::advanceToStage(): no rep");

    // Terrible things will happen if either of these conditions is not met:
    //   (1) the system is at stage g-1 now, AND
    //   (2) ALL subsystems have already been advanced to stage g.
    rep->advanceSystemToStage(g);
}

// We don't expect State entry allocations to be performance critical so
// we'll keep error checking on even in Release mode.

int State::allocateQ(int subsys, const Vector& qInit) {
    SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), Stage::Model, "State::allocateQ()");
    const int nxt = rep->subsystems[subsys].qInit.size();
    rep->subsystems[subsys].qInit.resizeKeep(nxt + qInit.size());
    rep->subsystems[subsys].qInit(nxt, qInit.size()) = qInit;
    return nxt;
}

int State::allocateU(int subsys, const Vector& uInit) {
    SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), Stage::Model, "State::allocateU()");
    const int nxt = rep->subsystems[subsys].uInit.size();
    rep->subsystems[subsys].uInit.resizeKeep(nxt + uInit.size());
    rep->subsystems[subsys].uInit(nxt, uInit.size()) = uInit;
    return nxt;
}
int State::allocateZ(int subsys, const Vector& zInit) {
    SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), Stage::Model, "State::allocateZ()");
    const int nxt = rep->subsystems[subsys].zInit.size();
    rep->subsystems[subsys].zInit.resizeKeep(nxt + zInit.size());
    rep->subsystems[subsys].zInit(nxt, zInit.size()) = zInit;
    return nxt;
}

int State::allocateQErr(int subsys, int nqerr) {
    SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), Stage::Model, "State::allocateQErr()");
    const int nxt = rep->subsystems[subsys].nqerr;
    rep->subsystems[subsys].nqerr += nqerr;
    return nxt;
}
int State::allocateUErr(int subsys, int nuerr) {
    SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), Stage::Model, "State::allocateUErr()");
    const int nxt = rep->subsystems[subsys].nuerr;
    rep->subsystems[subsys].nuerr += nuerr;
    return nxt;
}
int State::allocateUDotErr(int subsys, int nudoterr) {
    SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), Stage::Model, "State::allocateUDotErr()");
    const int nxt = rep->subsystems[subsys].nudoterr;
    rep->subsystems[subsys].nudoterr += nudoterr;
    return nxt;
}


// Topology- and Model-stage State variables can only be added during construction; that is,
// while stage <= Topology. Other entries can be added while stage < Model.
int State::allocateDiscreteVariable(int subsys, Stage g, AbstractValue* vp) {
    SimTK_STAGECHECK_RANGE_ALWAYS(Stage(Stage::LowestRuntime).prev(), g, Stage::HighestRuntime, 
        "State::allocateDiscreteVariable()");

    const Stage maxAcceptable = (g <= Stage::Model ? Stage::Empty : Stage::Topology);
    SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), 
        maxAcceptable.next(), "State::allocateDiscreteVariable()");

    PerSubsystemInfo& ss = rep->subsystems[subsys];
    const int nxt = ss.discrete.size();
    ss.discrete.push_back(DiscreteVariable(g,vp));
    return nxt;
}

// Cache entries can be allocated while stage < Model, even if they are Model-stage entries.
int State::allocateCacheEntry(int subsys, Stage g, AbstractValue* vp) {
    SimTK_STAGECHECK_RANGE_ALWAYS(Stage(Stage::LowestRuntime).prev(), g, Stage::HighestRuntime, 
        "State::allocateCacheEntry()");
    SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), 
        Stage::Model, "State::allocateCacheEntry()");

    PerSubsystemInfo& ss = rep->subsystems[subsys];
    const int nxt = ss.cache.size();
    ss.cache.push_back(CacheEntry(g,vp));
    return nxt;
}

    // State dimensions for shared continuous variables.

int State::getNY() const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getNY()");
    return rep->y.size();
}

int State::getQStart() const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getQStart()");
    return 0; // q's come first
}
int State::getNQ() const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getNQ()");
    return rep->q.size();
}

int State::getUStart() const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getUStart()");
    return rep->q.size(); // u's come right after q's
}
int State::getNU() const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getNU()");
    return rep->u.size();
}

int State::getZStart() const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getZStart()");
    return rep->q.size() + rep->u.size(); // q,u, then z
}
int State::getNZ() const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getNZ()");
    return rep->z.size();
}

int State::getNYErr() const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getNYErr()");
    return rep->yerr.size();
}

int State::getQErrStart() const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getQErrStart()");
    return 0; // qerr's come first
}
int State::getNQErr() const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getNQErr()");
    return rep->qerr.size();
}

int State::getUErrStart() const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getUErrStart()");
    return rep->qerr.size(); // uerr's follow qerrs
}
int State::getNUErr() const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getNUErr()");
    return rep->uerr.size();
}

// UDot errors are independent of qerr & uerr.
int State::getNUDotErr() const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getNUDotErr()");
    return rep->udoterr.size();
}

    // Subsystem dimensions.

int State::getQStart(int subsys) const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getQStart(subsys)");
    return rep->getSubsystem(subsys).qstart;
}
int State::getNQ(int subsys) const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getNQ(subsys)");
    return rep->getSubsystem(subsys).q.size();
}

int State::getUStart(int subsys) const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getUStart(subsys)");
    return rep->getSubsystem(subsys).ustart;
}
int State::getNU(int subsys) const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getNU(subsys)");
    return rep->getSubsystem(subsys).u.size();
}

int State::getZStart(int subsys) const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getZStart(subsys)");
    return rep->getSubsystem(subsys).zstart;
}
int State::getNZ(int subsys) const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getNZ(subsys)");
    return rep->getSubsystem(subsys).z.size();
}

int State::getQErrStart(int subsys) const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getQErrStart(subsys)");
    return rep->getSubsystem(subsys).qerrstart;
}
int State::getNQErr(int subsys) const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getNQErr(subsys)");
    return rep->getSubsystem(subsys).qerr.size();
}

int State::getUErrStart(int subsys) const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getUErrStart(subsys)");
    return rep->getSubsystem(subsys).uerrstart;
}
int State::getNUErr(int subsys) const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getNUErr(subsys)");
    return rep->getSubsystem(subsys).uerr.size();
}

int State::getUDotErrStart(int subsys) const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getUDotErrStart(subsys)");
    return rep->getSubsystem(subsys).udoterrstart;
}
int State::getNUDotErr(int subsys) const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getNUDotErr(subsys)");
    return rep->getSubsystem(subsys).udoterr.size();
}

    // Per-subsystem access to the global shared variables.

const Vector& State::getQ(int subsys) const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getQ(subsys)");
    return rep->getSubsystem(subsys).q;
}
const Vector& State::getU(int subsys) const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getU(subsys)");
    return rep->getSubsystem(subsys).u;
}
const Vector& State::getZ(int subsys) const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getZ(subsys)");
    return rep->getSubsystem(subsys).z;
}

const Vector& State::getQDot(int subsys) const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getQDot(subsys)");
    SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Velocity, "State::getQDot(subsys)");
    return rep->getSubsystem(subsys).qdot;
}
const Vector& State::getUDot(int subsys) const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getUDot(subsys)");
    SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Acceleration, "State::getUDot(subsys)");
    return rep->getSubsystem(subsys).udot;
}
const Vector& State::getZDot(int subsys) const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getZDot(subsys)");
    SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Dynamics, "State::getZDot(subsys)");
    return rep->getSubsystem(subsys).zdot;
}
const Vector& State::getQDotDot(int subsys) const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getQDotDot(subsys)");
    SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Acceleration, "State::getQDotDot(subsys)");
    return rep->getSubsystem(subsys).qdotdot;
}

Vector& State::updQ(int subsys) {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::updQ(subsys)");
    invalidateAll(Stage::Position);
    return rep->updSubsystem(subsys).q;
}
Vector& State::updU(int subsys) {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::updU(subsys)");
    invalidateAll(Stage::Velocity);
    return rep->updSubsystem(subsys).u;
}
Vector& State::updZ(int subsys) {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::updZ(subsys)");
    invalidateAll(Stage::Dynamics);
    return rep->updSubsystem(subsys).z;
}

    // These are mutable so the routines are const.

Vector& State::updQDot(int subsys) const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::updQDot(subsys)");
    SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage(Stage::Velocity).prev(), "State::updQDot(subsys)");
    return rep->getSubsystem(subsys).qdot;
}
Vector& State::updUDot(int subsys) const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::updUDot(subsys)");
    SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage(Stage::Acceleration).prev(), "State::updUDot(subsys)");
    return rep->getSubsystem(subsys).udot;
}
Vector& State::updZDot(int subsys) const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::updZDot(subsys)");
    SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage(Stage::Dynamics).prev(), "State::updZDot(subsys)");
    return rep->getSubsystem(subsys).zdot;
}
Vector& State::updQDotDot(int subsys) const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::updQDotDot(subsys)");
    SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage(Stage::Acceleration).prev(), "State::updQDotDot(subsys)");
    return rep->getSubsystem(subsys).qdotdot;
}


const Vector& State::getQErr(int subsys) const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getQErr(subsys)");
    SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Position, "State::getQErr(subsys)");
    return rep->getSubsystem(subsys).qerr;
}
const Vector& State::getUErr(int subsys) const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getUErr(subsys)");
    SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Velocity, "State::getUErr(subsys)");
    return rep->getSubsystem(subsys).uerr;
}
const Vector& State::getUDotErr(int subsys) const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getUDotErr(subsys)");
    SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Acceleration, "State::getUDotErr(subsys)");
    return rep->getSubsystem(subsys).udoterr;
}

Vector& State::updQErr(int subsys) const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::updQErr(subsys)");
    SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage(Stage::Position).prev(), "State::updQErr(subsys)");
    return rep->getSubsystem(subsys).qerr;
}
Vector& State::updUErr(int subsys) const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::updUErr(subsys)");
    SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage(Stage::Velocity).prev(), "State::updUErr(subsys)");
    return rep->getSubsystem(subsys).uerr;
}
Vector& State::updUDotErr(int subsys) const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::updUDotErr(subsys)");
    SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage(Stage::Acceleration).prev(), 
                        "State::updUDotErr(subsys)");
    return rep->getSubsystem(subsys).udoterr;
}

    // Direct access to the global shared state and cache entries.
    // These are allocated once the System Stage is Stage::Model.

const Real& State::getTime() const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getTime()");
    return rep->t;
}

const Vector& State::getY() const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getY()");
    return rep->y;
}

const Vector& State::getQ() const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getQ()");
    return rep->q;
}

const Vector& State::getU() const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getU()");
    return rep->u;
}

const Vector& State::getZ() const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::getZ()");
    return rep->z;
}


// You can call these as long as stage >= Model, but the
// stage will be backed up if necessary to the indicated stage.
Real& State::updTime() {  // Stage::Time-1
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::updTime()");
    invalidateAll(Stage::Time);
    return rep->t;
}

Vector& State::updY() {    // Back to Stage::Position-1
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::updY()");
    invalidateAll(Stage::Position);
    return rep->y;
}

Vector& State::updQ() {    // Stage::Position-1
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::updQ()");
    invalidateAll(Stage::Position);
    return rep->q;
}

Vector& State::updU() {     // Stage::Velocity-1
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::updU()");
    invalidateAll(Stage::Velocity);
    return rep->u;
}

Vector& State::updZ() {     // Stage::Dynamics-1
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "State::updZ()");
    invalidateAll(Stage::Dynamics);
    return rep->z;
}

const Vector& State::getYDot() const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Acceleration, "State::getYDot()");
    return rep->ydot;
}

const Vector& State::getQDot() const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Velocity, "State::getQDot()");
    return rep->qdot;
}

const Vector& State::getZDot() const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Dynamics, "State::getZDot()");
    return rep->zdot;
}

const Vector& State::getUDot() const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Acceleration, "State::getUDot()");
    return rep->udot;
}

const Vector& State::getQDotDot() const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Acceleration, "State::getQDotDot()");
    return rep->qdotdot;
}


Vector& State::updYDot() const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage(Stage::Acceleration).prev(), "State::updYDot()");
    return rep->ydot;
}

Vector& State::updQDot() const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage(Stage::Velocity).prev(), "State::updQDot()");
    return rep->qdot;
}

Vector& State::updUDot() const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage(Stage::Acceleration).prev(), "State::updUDot()");
    return rep->udot;
}

Vector& State::updZDot() const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage(Stage::Dynamics).prev(), "State::updZDot()");
    return rep->zdot;
}

Vector& State::updQDotDot() const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage(Stage::Acceleration).prev(), "State::updQDotDot()");
    return rep->qdotdot;
}


const Vector& State::getYErr() const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Velocity, "State::getYErr()");
    return rep->yerr;
}

const Vector& State::getQErr() const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Position, "State::getQErr()");
    return rep->qerr;
}
const Vector& State::getUErr() const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Velocity, "State::getUErr()");
    return rep->uerr;
}
const Vector& State::getUDotErr() const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Acceleration, "State::getUDotErr()");
    return rep->udoterr;
}

Vector& State::updYErr() const {
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage(Stage::Velocity).prev(), "State::updYErr()");
    return rep->yerr;
}
Vector& State::updQErr() const{
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage(Stage::Position).prev(), "State::updQErr()");
    return rep->qerr;
}
Vector& State::updUErr() const{
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage(Stage::Velocity).prev(), "State::updUErr()");
    return rep->uerr;
}
Vector& State::updUDotErr() const{
    assert(rep);
    SimTK_STAGECHECK_GE(getSystemStage(), Stage(Stage::Acceleration).prev(), 
                        "State::updUDotErr()");
    return rep->udoterr;
}

// You can access a Model stage variable any time, but don't access others
// until you have realized the Model stage.
const AbstractValue& 
State::getDiscreteVariable(int subsys, int index) const {
    const PerSubsystemInfo& ss = rep->subsystems[subsys];

    SimTK_INDEXCHECK(0,index,(int)ss.discrete.size(),"State::getDiscreteVariable()");
    const DiscreteVariable& dv = ss.discrete[index];

    if (dv.getStage() > Stage::Model) {
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), 
            Stage::Model, "State::getDiscreteVariable()");
    }

    return dv.getValue();
}

// You can update a Model stage variable from Topology stage, but higher variables 
// must wait until you have realized the Model stage. This always backs the 
// stage up to one earlier than the variable's stage.
AbstractValue& 
State::updDiscreteVariable(int subsys, int index) {
    PerSubsystemInfo& ss = rep->subsystems[subsys];

    SimTK_INDEXCHECK(0,index,(int)ss.discrete.size(),"State::updDiscreteVariable()");
    DiscreteVariable& dv = ss.discrete[index];

    SimTK_STAGECHECK_GE(getSubsystemStage(subsys), 
        std::min(dv.getStage().prev(), Stage(Stage::Model)), 
        "State::updDiscreteVariable()");

    invalidateAll(dv.getStage());

    return dv.updValue();
}

// Stage >= ce.stage
const AbstractValue& 
State::getCacheEntry(int subsys, int index) const {
    const PerSubsystemInfo& ss = rep->subsystems[subsys];

    SimTK_INDEXCHECK(0,index,(int)ss.cache.size(),"State::getCacheEntry()");
    const CacheEntry& ce = ss.cache[index];

    SimTK_STAGECHECK_GE(getSubsystemStage(subsys), 
        ce.getStage(), "State::getCacheEntry()");

    return ce.getValue();
}

// Stage >= ce.stage-1; does not change stage
AbstractValue& 
State::updCacheEntry(int subsys, int index) const {
    const PerSubsystemInfo& ss = rep->subsystems[subsys];

    SimTK_INDEXCHECK(0,index,(int)ss.cache.size(),"State::updCacheEntry()");
    CacheEntry& ce = ss.cache[index];

    SimTK_STAGECHECK_GE(getSubsystemStage(subsys), 
        ce.getStage().prev(), "State::updCacheEntry()");

    return ce.updValue();
}

String State::toString() const {
    String out;
    out += "<State>\n";

    out += "<Real name=time>" + String(rep->t) + "</Real>\n";

    out += "<Vector name=q size=" + String(rep->q.size()) + ">";
    if (rep->q.size()) out += "\n";
    for (long i=0; i<rep->q.size(); ++i)
        out += String(rep->q[i]) + "\n";
    out += "</Vector>\n";

    out += "<Vector name=u size=" + String(rep->u.size()) + ">";
    if (rep->u.size()) out += "\n";
    for (long i=0; i<rep->u.size(); ++i)
        out += String(rep->u[i]) + "\n";
    out += "</Vector>\n";

    out += "<Vector name=z size=" + String(rep->z.size()) + ">";
    if (rep->z.size()) out += "\n";
    for (long i=0; i<rep->z.size(); ++i)
        out += String(rep->z[i]) + "\n";
    out += "</Vector>\n";


    for (int ss=0; ss < (int)rep->subsystems.size(); ++ss) {
        const PerSubsystemInfo& info = rep->subsystems[ss];
        out += "<Subsystem index=" + String(ss) + " name=" + info.name 
            + " version=" + info.version + ">\n";

        out += "  <DISCRETE VARS TODO>\n";

        out += "  <Vector name=qInit size=" + String(info.qInit.size()) + ">\n";
        out += "  <Vector name=uInit size=" + String(info.uInit.size()) + ">\n";
        out += "  <Vector name=zInit size=" + String(info.zInit.size()) + ">\n";

        out += "  <Vector name=q size=" + String(info.q.size()) + ">\n";
        out += "  <Vector name=u size=" + String(info.u.size()) + ">\n";
        out += "  <Vector name=z size=" + String(info.z.size()) + ">\n";

        out += "</Subsystem>\n";
    }

    out += "</State>\n";
    return out;
}

String State::cacheToString() const {
    String out;
    out += "<Cache>\n";
    out += "<Stage>" + getSystemStage().name() + "</Stage>\n";

    for (int ss=0; ss < (int)rep->subsystems.size(); ++ss) {
        const PerSubsystemInfo& info = rep->subsystems[ss];
        out += "<Subsystem index=" + String(ss) + " name=" + info.name 
            + " version=" + info.version + ">\n";
        out += "  <Stage>" + info.currentStage.name() + "</Stage>\n";

        out += "  <DISCRETE CACHE TODO>\n";

        out += "  <Vector name=qdot size=" + String(info.qdot.size()) + ">\n";
        out += "  <Vector name=udot size=" + String(info.udot.size()) + ">\n";
        out += "  <Vector name=zdot size=" + String(info.zdot.size()) + ">\n";
        out += "  <Vector name=qdotdot size=" + String(info.qdotdot.size()) + ">\n";

        out += "  <Vector name=qerr size=" + String(info.qerr.size()) + ">\n";
        out += "  <Vector name=uerr size=" + String(info.uerr.size()) + ">\n";
        out += "  <Vector name=udoterr size=" + String(info.udoterr.size()) + ">\n";

        out += "</Subsystem>\n";
    }

    out += "<Vector name=qdot size=" + String(rep->qdot.size()) + ">";
    if (rep->qdot.size()) out += "\n";
    for (long i=0; i<rep->qdot.size(); ++i)
        out += String(rep->qdot[i]) + "\n";
    out += "</Vector>\n";

    out += "<Vector name=udot size=" + String(rep->udot.size()) + ">";
    if (rep->udot.size()) out += "\n";
    for (long i=0; i<rep->udot.size(); ++i)
        out += String(rep->udot[i]) + "\n";
    out += "</Vector>\n";

    out += "<Vector name=zdot size=" + String(rep->zdot.size()) + ">";
    if (rep->zdot.size()) out += "\n";
    for (long i=0; i<rep->zdot.size(); ++i)
        out += String(rep->zdot[i]) + "\n";
    out += "</Vector>\n";

    out += "<Vector name=qdotdot size=" + String(rep->qdotdot.size()) + ">";
    if (rep->qdotdot.size()) out += "\n";
    for (long i=0; i<rep->qdotdot.size(); ++i)
        out += String(rep->qdotdot[i]) + "\n";
    out += "</Vector>\n";

    out += "<Vector name=qerr size=" + String(rep->qerr.size()) + ">";
    if (rep->qerr.size()) out += "\n";
    for (long i=0; i<rep->qerr.size(); ++i)
        out += String(rep->qerr[i]) + "\n";
    out += "</Vector>\n";

    out += "<Vector name=uerr size=" + String(rep->uerr.size()) + ">";
    if (rep->uerr.size()) out += "\n";
    for (long i=0; i<rep->uerr.size(); ++i)
        out += String(rep->uerr[i]) + "\n";
    out += "</Vector>\n";

    out += "<Vector name=udoterr size=" + String(rep->udoterr.size()) + ">";
    if (rep->udoterr.size()) out += "\n";
    for (long i=0; i<rep->udoterr.size(); ++i)
        out += String(rep->udoterr[i]) + "\n";
    out += "</Vector>\n";

    out += "</Cache>\n";
    return out;
}


std::ostream& 
operator<<(std::ostream& o, const State& s) {
    o << "STATE:" << std::endl;
    o << s.toString() << std::endl;
    o << "CACHE:" << std::endl;
    return o << s.cacheToString() << std::endl;
}

} // namespace SimTK

