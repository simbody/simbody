/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-7 Stanford University and the Authors.         *
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

#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"
#include "SimTKcommon/internal/State.h"

#include <cassert>
#include <ostream>
#include <set>

using std::set;

namespace SimTK {

class DiscreteVariableRep {
public:
    DiscreteVariableRep() : value(0), myHandle(0), stage(Stage::Empty) { }
    DiscreteVariableRep(const DiscreteVariableRep& src) : stage(src.stage) {
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
        assert(g == Stage::Topology || g.isInRuntimeRange());
        assert(vp);
    }

    DiscreteVariableRep* clone() const {return new DiscreteVariableRep(*this);}

    bool isValid() const {
        return (stage==Stage::Topology || stage.isInRuntimeRange()) && value && myHandle; 
    }

    Stage getStage() const {return stage;}

    const AbstractValue& getValue() const {assert(value); return *value;}
    AbstractValue&       updValue()       {assert(value); return *value;}

    void                    setMyHandle(DiscreteVariable& dv) {myHandle = &dv;}
    const DiscreteVariable& getMyHandle() const   {assert(myHandle); return *myHandle;}
    DiscreteVariable&       updMyHandle()         {assert(myHandle); return *myHandle;}
private:
    void clear() {stage=Stage::Empty; delete value;}

    Stage stage;
    AbstractValue* value;

    DiscreteVariable* myHandle;
};


/*static*/ String 
EventStatus::eventTriggerString(EventTrigger e) {
    // Catch special combos first
    if (e==NoEventTrigger)        return "NoEventTrigger";
    if (e==Falling)               return "Falling";
    if (e==Rising)                return "Rising";
    if (e==AnySignChange)         return "AnySignChange";

    // Not a special combo; unmask one at a time.
    const EventTrigger triggers[] =
     { PositiveToNegative,NegativeToPositive,NoEventTrigger };
    const char *triggerNames[] =
     { "PositiveToNegative","NegativeToPositive" };

    String s;
    for (int i=0; triggers[i] != NoEventTrigger; ++i)
        if (e & triggers[i]) {
            if (s.size()) s += "|";
            s += triggerNames[i];
            e = EventTrigger((unsigned)e & ~((unsigned)triggers[i])); 
        }

    // should have accounted for everything by now
    if (e != NoEventTrigger) {
        char buf[128];
        std::sprintf(buf, "0x%x", (unsigned)e);
        if (s.size()) s += " + ";
        s += "UNRECOGNIZED EVENT TRIGGER GARBAGE ";
        s += buf;
    }
    return s;
}

// This internal utility class is used to capture all the information needed for
// a single subsystem within the StateData.
class PerSubsystemInfo {
public:
    PerSubsystemInfo() : currentStage(Stage::Empty)     {initialize();}
    PerSubsystemInfo(const String& n, const String& v) 
      : name(n), version(v), currentStage(Stage::Empty) {initialize();}

    ~PerSubsystemInfo() {   // default destructor
    }

    // Copy constructor copies all variables but cache only through
    // modeled stage. Note that this must be done in conjunction with
    // copying the whole state or our global resource indices will
    // be nonsense.
    PerSubsystemInfo(const PerSubsystemInfo& src) : currentStage(Stage::Empty) {
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

        for (int j=0; j<Stage::NValid; ++j) {
            eventstart[j] = -1;
            events[j].clear();
        }
    }

    String name;
    String version;

    // These accumulate default values for this subsystem's use of shared
    // global state variables. After the System is advanced to Stage::Model,
    // the state will allocate those globals and copy these initial
    // values into them. The lengths of these Vectors define the 
    // needs of this Subsystem.
    Vector qInit, uInit, zInit;

    // For constraints we need just lengths (nmultipliers==nudoterr).
    int nqerr, nuerr, nudoterr;

    // For events we need just lengths, for each stage.
    int nevents[Stage::NValid];

    // These are our own private views into partitions of the global
    // state and cache entries of the same names. These are valid
    // only after the *System* stage is raised to Model, and they
    // are invalidated whenever the System's Model stage is invalidated.

    // The State will assign contiguous blocks to this subsystem. The
    // starting indices are filled in here at the time the views are built.
    // Note that multipliers just use the same indices as udoterr.
    int qstart, ustart, zstart, qerrstart, uerrstart, udoterrstart;
    int eventstart[Stage::NValid];
    Vector q, u, z;
    mutable Vector qdot, udot, zdot, qdotdot;
    mutable Vector qerr, uerr;
    mutable Vector udoterr, multipliers; // same size and partioning
    mutable Vector events[Stage::NValid];

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
        for (int j=0; j<Stage::NValid; ++j) {
            nevents[j] = 0;
            eventstart[j] = -1;
        }
        nDiscreteWhenBuilt = -1;
        for (int j=0; j<Stage::NValid; ++j)
            cacheSize[j] = -1;
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
        for (int j=0; j<Stage::NValid; ++j)
            nevents[j] = 0;
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
        for (int j=0; j<Stage::NValid; ++j)
            nevents[j] = src.nevents[j];

        // Copy *all* discrete state variables since no more can be allocated
        // after Model stage. Also, make sure we know how to back up to Stage::Topology
        // later if necessary.
        discrete = src.discrete;
        nDiscreteWhenBuilt = src.nDiscreteWhenBuilt;

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

class StateData {
public:
    StateData() 
      : t(CNT<Real>::getNaN()), systemStage(Stage::Empty), 
        myHandle(0) 
    { 
    }

    ~StateData() {   // default destructor
    }

    const Stage& getSystemStage() const {
        return systemStage;
    }
    Stage& updSystemStage() const {
        return systemStage; // mutable
    }

    const PerSubsystemInfo& getSubsystem(int subsystem) const {
        SimTK_INDEXCHECK(0, subsystem, (int)subsystems.size(), "StateData::getSubsystem()");
        return subsystems[subsystem];
    }

    PerSubsystemInfo& updSubsystem(int subsystem) {
        SimTK_INDEXCHECK(0, subsystem, (int)subsystems.size(), "StateData::updSubsystem()");
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
    StateData(const StateData& src) : myHandle(0), systemStage(Stage::Empty) {
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

    StateData& operator=(const StateData& src) {
        if (&src == this) return *this;
        invalidateJustSystemStage(Stage::Topology);
        for (SubsystemIndex i(0); i<(int)subsystems.size(); ++i)
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
    StateData* clone() const {return new StateData(*this);}

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
                for (SubsystemIndex i(0); i < (int)subsystems.size(); ++i)
                    subsystems[i].clearReferencesToStateGlobals();

                t = CNT<Real>::getNaN();
                // Nuke all the global views.
                q.clear(); u.clear(); z.clear();
                qdot.clear(); udot.clear(); zdot.clear();
                qerr.clear(); uerr.clear();
                for (int j=0; j<Stage::NValid; ++j)
                    events[j].clear();

                // Nuke the actual data.
                y.unlockShape();           y.clear(); 
                ydot.unlockShape();        ydot.clear(); 
                qdotdot.unlockShape();     qdotdot.clear();
                yerr.unlockShape();        yerr.clear();
                udoterr.unlockShape();     udoterr.clear();
                multipliers.unlockShape(); multipliers.clear();
                allEvents.unlockShape();   allEvents.clear();
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
            int nq=0, nu=0, nz=0, nqerr=0, nuerr=0, nudoterr=0, nAllEvents=0;
            int nevents[Stage::NValid];
            for (int j=0; j<Stage::NValid; ++j)
                nevents[j] = 0;

            // Count up all 
            for (SubsystemIndex i(0); i<(int)subsystems.size(); ++i) {
                nq += subsystems[i].qInit.size();
                nu += subsystems[i].uInit.size();
                nz += subsystems[i].zInit.size();
                nqerr    += subsystems[i].nqerr;
                nuerr    += subsystems[i].nuerr;
                nudoterr += subsystems[i].nudoterr;

                for (int j=0; j<Stage::NValid; ++j)
                    nevents[j] += subsystems[i].nevents[j];
            }
            for (int j=0; j<Stage::NValid; ++j)
                nAllEvents += nevents[j];

            // Allocate the actual shared state variables & cache 
            // entries and make sure no one can accidentally change the size.
            y.resize(nq+nu+nz);             y.lockShape();
            ydot.resize(nq+nu+nz);          ydot.lockShape();
            qdotdot.resize(nq);             qdotdot.lockShape();
            yerr.resize(nqerr+nuerr);       yerr.lockShape();
            udoterr.resize(nudoterr);       udoterr.lockShape();
            multipliers.resize(nudoterr);   multipliers.lockShape(); // same size as udoterr
            allEvents.resize(nAllEvents);   allEvents.lockShape();

            // Allocate subviews of the shared state & cache entries.
            q.viewAssign(y(0,nq));
            u.viewAssign(y(nq,nu));
            z.viewAssign(y(nq+nu,nz));

            qdot.viewAssign(ydot(0,     nq));
            udot.viewAssign(ydot(nq,    nu));
            zdot.viewAssign(ydot(nq+nu, nz));

            qerr.viewAssign(yerr(0,     nqerr));
            uerr.viewAssign(yerr(nqerr, nuerr));

            int stageStart=0;
            for (int j=0; j<Stage::NValid; ++j) {
                events[j].viewAssign(allEvents(stageStart, nevents[j]));
                stageStart += nevents[j];
            }

            // Now partition the global resources among the subsystems and copy
            // in the initial values for the state variables.
            int nxtq=0, nxtu=0, nxtz=0, nxtqerr=0, nxtuerr=0, nxtudoterr=0;
            int nxtevent[Stage::NValid];
            for (int j=0; j<Stage::NValid; ++j)
                nxtevent[j] = 0;

            for (SubsystemIndex i(0); i<(int)subsystems.size(); ++i) {
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
                // multipliers have same partitioning as udoterr
                ss.multipliers.viewAssign(multipliers(nxtudoterr, ss.nudoterr));

                // Consume the slots.
                nxtq += nq; nxtu += nu; nxtz += nz;
                nxtqerr += ss.nqerr; nxtuerr += ss.nuerr; nxtudoterr += ss.nudoterr;

                // Same thing for event slots, but by stage.
                for (int j=0; j<Stage::NValid; ++j) {
                    ss.eventstart[j] = nxtevent[j];
                    ss.events[j].viewAssign(events[j](nxtevent[j], ss.nevents[j]));
                    nxtevent[j] += ss.nevents[j];
                }

            }

            // As the final "modeling" step, initialize time to 0 (it's NaN before this).
            t = 0;
        }

        systemStage = g;
    }

    void            setMyHandle(StateRep& s) {myHandle = &s;}
    const StateRep& getMyHandle() const   {assert(myHandle); return *myHandle;}
    StateRep&       updMyHandle()         {assert(myHandle); return *myHandle;}
private:
    friend class StateRep;

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

    // All the state derivatives taken together (qdot,udot,zdot)
    mutable Vector  ydot; 

    // These are views into ydot.
    mutable Vector  qdot;       // Stage::Velocity
    mutable Vector  udot;       // Stage::Acceleration
    mutable Vector  zdot;       // Stage::Acceleration

    // This is an independent cache entry.
    mutable Vector  qdotdot;    // Stage::Acceleration

        // ALGEBRAIC EQUATIONS

    mutable Vector  yerr;        // All constraint errors taken together (qerr,uerr)
    mutable Vector  udoterr;     // Stage::Acceleration (Index 1 constraints)
    mutable Vector  multipliers; // Stage::Acceleration (Index 1 algebraic variables)

    // These are views into yerr.
    mutable Vector  qerr;       // Stage::Position (Index 3 constraints)
    mutable Vector  uerr;       // Stage::Velocity (Index 2 constraints)

        // DISCRETE EQUATIONS

    // All the events together, ordered by stage.
    mutable Vector  allEvents;

    // These are views into allEvents.
    mutable Vector  events[Stage::NValid];
    

        // Subsystem support //

    // Subsystem 0 (always present) is for the System as a whole. Its name
    // and version are the System name and version.
    std::vector<PerSubsystemInfo> subsystems;

    // Return true only if all subsystems have been realized to at least Stage g.
    bool allSubsystemsAtLeastAtStage(Stage g) const {
        for (SubsystemIndex i(0); i < (int)subsystems.size(); ++i)
            if (subsystems[i].currentStage < g)
                return false;
        return true;
    }

private:
    StateRep* myHandle;
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



class StateRep {
public:
    StateRep() : data(new StateData()) {
        data->setMyHandle(*this);
    }
    
    
    // Restore state to default-constructed condition
    void clear() {
        delete data;
        data = new StateData();
        data->setMyHandle(*this);
    }
    
    ~StateRep() {
        if (data->myHandle == this)
            delete data;
        data=0;
    }
    
    // copy constructor
    StateRep(const StateRep& src) 
      : data(src.data->clone()) {
        data->setMyHandle(*this);
    }
    
    // This constructor creates restricted states.
    StateRep(const StateRep& src, EnumerationSet<Stage>& restrictedStages, set<SubsystemIndex>& restrictedSubsystems) 
      : data(src.data), restrictedStages(restrictedStages), restrictedSubsystems(restrictedSubsystems) {
    }
    
    // copy assignment
    StateRep& operator=(const StateRep& src) {
        if (&src == this) return *this;
        if (!data) {
            // we're defining this state here (if src is not empty)
            if (src.data) {
                data = src.data->clone();
                data->setMyHandle(*this);
            }
            return *this;
        }
    
        // Assignment or redefinition
        if (src.data) *data = *src.data;
        else {delete data; data=0;}
        return *this;
    }
    
    void setNSubsystems(int i) {
        assert(i >= 0);
        data->subsystems.clear();
        data->subsystems.resize(i);
    }
    
    void initializeSubsystem(SubsystemIndex i, const String& name, const String& version) {
        data->updSubsystem(i).name = name;
        data->updSubsystem(i).version = version;
    }
    
    
    int addSubsystem(const String& name, const String& version) {
        data->subsystems.push_back(
            PerSubsystemInfo(name,version));
        return (int)data->subsystems.size() - 1;
    }
    
    int getNSubsystems() const {return (int)data->subsystems.size();}
    
    const String& getSubsystemName(SubsystemIndex subsys) const {
        return data->subsystems[subsys].name;
    }
    const String& getSubsystemVersion(SubsystemIndex subsys) const {
        return data->subsystems[subsys].version;
    }
    
    const Stage& getSystemStage() const {
        return data->getSystemStage();
    }
    
    const Stage& getSubsystemStage(SubsystemIndex subsys) const {
        SimTK_ASSERT(data, "StateRep::getStage(): no data"); // can't happen(?)
        return data->getSubsystemStage(subsys);
    }
    
    // Make sure the stage is no higher than g-1 for *any* subsystem and
    // hence for the system stage also. TODO: this should be more selective.
    void invalidateAll(Stage g) const {
        SimTK_ASSERT(data, "StateRep::invalidateAll(): no data");
    
        data->invalidateJustSystemStage(g);
        for (SubsystemIndex i(0); i<(int)data->subsystems.size(); ++i)
            data->subsystems[i].invalidateStageJustThisSubsystem(g);
    }
    
    // Move the stage for a particular subsystem from g-1 to g. No other subsystems
    // are affected, nor the global system stage.
    void advanceSubsystemToStage(SubsystemIndex subsys, Stage g) const {
        SimTK_ASSERT(data, "StateRep::advanceSubsystemToStage(): no data");
    
        data->subsystems[subsys].advanceToStage(g);
        // We don't automatically advance the System stage even if this brings
        // ALL the subsystems up to stage g.
    }
    
    // Move the system stage from g-1 to g. Don't call this until ALL 
    // subsystem have been advanced to at least stage g.
    void advanceSystemToStage(Stage g) const {
        SimTK_ASSERT(data, "StateRep::advanceToStage(): no data");
    
        // Terrible things will happen if either of these conditions is not met:
        //   (1) the system is at stage g-1 now, AND
        //   (2) ALL subsystems have already been advanced to stage g.
        data->advanceSystemToStage(g);
    }
    
    // We don't expect State entry allocations to be performance critical so
    // we'll keep error checking on even in Release mode.
    
    int allocateQ(SubsystemIndex subsys, const Vector& qInit) {
        SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), Stage::Model, "StateRep::allocateQ()");
        checkCanModify(subsys);
        const int nxt = data->subsystems[subsys].qInit.size();
        data->subsystems[subsys].qInit.resizeKeep(nxt + qInit.size());
        data->subsystems[subsys].qInit(nxt, qInit.size()) = qInit;
        return nxt;
    }
    
    int allocateU(SubsystemIndex subsys, const Vector& uInit) {
        SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), Stage::Model, "StateRep::allocateU()");
        checkCanModify(subsys);
        const int nxt = data->subsystems[subsys].uInit.size();
        data->subsystems[subsys].uInit.resizeKeep(nxt + uInit.size());
        data->subsystems[subsys].uInit(nxt, uInit.size()) = uInit;
        return nxt;
    }
    int allocateZ(SubsystemIndex subsys, const Vector& zInit) {
        SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), Stage::Model, "StateRep::allocateZ()");
        checkCanModify(subsys);
        const int nxt = data->subsystems[subsys].zInit.size();
        data->subsystems[subsys].zInit.resizeKeep(nxt + zInit.size());
        data->subsystems[subsys].zInit(nxt, zInit.size()) = zInit;
        return nxt;
    }
    
    int allocateQErr(SubsystemIndex subsys, int nqerr) {
        SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), Stage::Model, "StateRep::allocateQErr()");
        checkCanModify(subsys);
        const int nxt = data->subsystems[subsys].nqerr;
        data->subsystems[subsys].nqerr += nqerr;
        return nxt;
    }
    int allocateUErr(SubsystemIndex subsys, int nuerr) {
        SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), Stage::Model, "StateRep::al()");
        checkCanModify(subsys);
        const int nxt = data->subsystems[subsys].nuerr;
        data->subsystems[subsys].nuerr += nuerr;
        return nxt;
    }
    int allocateUDotErr(SubsystemIndex subsys, int nudoterr) {
        SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), Stage::Model, "StateRep::allocateUDotErr()");
        checkCanModify(subsys);
        const int nxt = data->subsystems[subsys].nudoterr;
        data->subsystems[subsys].nudoterr += nudoterr;
        return nxt;
    }
    int allocateEvent(SubsystemIndex subsys, Stage g, int ne) {
        SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), Stage::Model, "StateRep::allocateEvent()");
        checkCanModify(g);
        checkCanModify(subsys);
        const int nxt = data->subsystems[subsys].nevents[g];
        data->subsystems[subsys].nevents[g] += ne;
        return nxt;
    }
    
    // Topology- and Model-stage State variables can only be added during construction; that is,
    // while stage <= Topology. Other entries can be added while stage < Model.
    int allocateDiscreteVariable(SubsystemIndex subsys, Stage g, AbstractValue* vp) {
        SimTK_STAGECHECK_RANGE_ALWAYS(Stage(Stage::LowestRuntime).prev(), g, Stage::HighestRuntime, 
            "StateRep::allocateDiscreteVariable()");
    
        const Stage maxAcceptable = (g <= Stage::Model ? Stage::Empty : Stage::Topology);
        SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), 
            maxAcceptable.next(), "StateRep::allocateDiscreteVariable()");
        checkCanModify(g);
        checkCanModify(subsys);
    
        PerSubsystemInfo& ss = data->subsystems[subsys];
        const int nxt = ss.discrete.size();
        ss.discrete.push_back(DiscreteVariable(g,vp));
        return nxt;
    }
    
    // Cache entries can be allocated while stage < Model, even if they are Model-stage entries.
    int allocateCacheEntry(SubsystemIndex subsys, Stage g, AbstractValue* vp) {
        SimTK_STAGECHECK_RANGE_ALWAYS(Stage(Stage::LowestRuntime).prev(), g, Stage::HighestRuntime, 
            "StateRep::allocateCacheEntry()");
        SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), 
            Stage::Model, "StateRep::allocateCacheEntry()");
        checkCanModify(g);
        checkCanModify(subsys);

        PerSubsystemInfo& ss = data->subsystems[subsys];
        const int nxt = ss.cache.size();
        ss.cache.push_back(CacheEntry(g,vp));
        return nxt;
    }
    
        // State dimensions for shared continuous variables.
    
    int getNY() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getNY()");
        return data->y.size();
    }
    
    int getQStart() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getQStart()");
        return 0; // q's come first
    }
    int getNQ() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getNQ()");
        return data->q.size();
    }
    
    int getUStart() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getUStart()");
        return data->q.size(); // u's come right after q's
    }
    int getNU() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getNU()");
        return data->u.size();
    }
    
    int getZStart() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getZStart()");
        return data->q.size() + data->u.size(); // q,u, then z
    }
    int getNZ() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getNZ()");
        return data->z.size();
    }
    
    int getNYErr() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getNYErr()");
        return data->yerr.size();
    }
    
    int getQErrStart() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getQErrStart()");
        return 0; // qerr's come first
    }
    int getNQErr() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getNQErr()");
        return data->qerr.size();
    }
    
    int getUErrStart() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getUErrStart()");
        return data->qerr.size(); // uerr's follow qerrs
    }
    int getNUErr() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getNUErr()");
        return data->uerr.size();
    }
    
    // UDot errors are independent of qerr & uerr.
    // This is used for multipliers also.
    int getNUDotErr() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getNUDotErr()");
        return data->udoterr.size();
    }
    
    int getNEvents() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getNEvents()");
        return data->allEvents.size();
    }
    
    int getEventStartByStage(Stage g) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getEventStartByStage()");
        int nxt = 0;
        for (int j=0; j<g; ++j)
            nxt += data->events[j].size();
        return nxt; // g starts where g-1 leaves off
    }
    
    int getNEventsByStage(Stage g) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getNEventsByStage()");
        return data->events[g].size();
    }
    
        // Subsystem dimensions.
    
    int getQStart(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getQStart(subsys)");
        return data->getSubsystem(subsys).qstart;
    }
    int getNQ(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getNQ(subsys)");
        return data->getSubsystem(subsys).q.size();
    }
    
    int getUStart(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getUStart(subsys)");
        return data->getSubsystem(subsys).ustart;
    }
    int getNU(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getNU(subsys)");
        return data->getSubsystem(subsys).u.size();
    }
    
    int getZStart(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getZStart(subsys)");
        return data->getSubsystem(subsys).zstart;
    }
    int getNZ(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getNZ(subsys)");
        return data->getSubsystem(subsys).z.size();
    }
    
    int getQErrStart(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getQErrStart(subsys)");
        return data->getSubsystem(subsys).qerrstart;
    }
    int getNQErr(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getNQErr(subsys)");
        return data->getSubsystem(subsys).qerr.size();
    }
    
    int getUErrStart(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getUErrStart(subsys)");
        return data->getSubsystem(subsys).uerrstart;
    }
    int getNUErr(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getNUErr(subsys)");
        return data->getSubsystem(subsys).uerr.size();
    }
    
    // These are used for multipliers also.
    int getUDotErrStart(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getUDotErrStart(subsys)");
        return data->getSubsystem(subsys).udoterrstart;
    }
    int getNUDotErr(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getNUDotErr(subsys)");
        return data->getSubsystem(subsys).udoterr.size();
    }
    
    int getEventStartByStage(SubsystemIndex subsys, Stage g) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getEventStartByStage(subsys)");
        return data->getSubsystem(subsys).eventstart[g];
    }
    
    int getNEventsByStage(SubsystemIndex subsys, Stage g) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getNEventsByStage(subsys)");
        return data->getSubsystem(subsys).events[g].size();
    }
    
        // Per-subsystem access to the global shared variables.
    
    const Vector& getQ(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getQ(subsys)");
        return data->getSubsystem(subsys).q;
    }
    const Vector& getU(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getU(subsys)");
        return data->getSubsystem(subsys).u;
    }
    const Vector& getZ(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getZ(subsys)");
        return data->getSubsystem(subsys).z;
    }
    
    const Vector& getQDot(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getQDot(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Velocity, "StateRep::getQDot(subsys)");
        return data->getSubsystem(subsys).qdot;
    }
    const Vector& getUDot(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getUDot(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Acceleration, "StateRep::getUDot(subsys)");
        return data->getSubsystem(subsys).udot;
    }
    const Vector& getZDot(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getZDot(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Dynamics, "StateRep::getZDot(subsys)");
        return data->getSubsystem(subsys).zdot;
    }
    const Vector& getQDotDot(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getQDotDot(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Acceleration, "StateRep::getQDotDot(subsys)");
        return data->getSubsystem(subsys).qdotdot;
    }
    
    Vector& updQ(SubsystemIndex subsys) {
        assert(data);
        checkCanModify(Stage::Position);
        checkCanModify(subsys);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::updQ(subsys)");
        invalidateAll(Stage::Position);
        return data->updSubsystem(subsys).q;
    }
    Vector& updU(SubsystemIndex subsys) {
        assert(data);
        checkCanModify(Stage::Velocity);
        checkCanModify(subsys);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::updU(subsys)");
        invalidateAll(Stage::Velocity);
        return data->updSubsystem(subsys).u;
    }
    Vector& updZ(SubsystemIndex subsys) {
        assert(data);
        checkCanModify(Stage::Dynamics);
        checkCanModify(subsys);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::updZ(subsys)");
        invalidateAll(Stage::Dynamics);
        return data->updSubsystem(subsys).z;
    }
    
        // These are mutable so the routines are const.
    
    Vector& updQDot(SubsystemIndex subsys) const {
        assert(data);
        checkCanModify(Stage::Position);
        checkCanModify(subsys);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::updQDot(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage(Stage::Velocity).prev(), "StateRep::updQDot(subsys)");
        return data->getSubsystem(subsys).qdot;
    }
    Vector& updUDot(SubsystemIndex subsys) const {
        assert(data);
        checkCanModify(Stage::Velocity);
        checkCanModify(subsys);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::updUDot(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage(Stage::Acceleration).prev(), "StateRep::updUDot(subsys)");
        return data->getSubsystem(subsys).udot;
    }
    Vector& updZDot(SubsystemIndex subsys) const {
        assert(data);
        checkCanModify(Stage::Dynamics);
        checkCanModify(subsys);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::updZDot(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage(Stage::Dynamics).prev(), "StateRep::updZDot(subsys)");
        return data->getSubsystem(subsys).zdot;
    }
    Vector& updQDotDot(SubsystemIndex subsys) const {
        assert(data);
        checkCanModify(Stage::Position);
        checkCanModify(subsys);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::updQDotDot(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage(Stage::Acceleration).prev(), "StateRep::updQDotDot(subsys)");
        return data->getSubsystem(subsys).qdotdot;
    }
    
    
    const Vector& getQErr(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getQErr(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Position, "StateRep::getQErr(subsys)");
        return data->getSubsystem(subsys).qerr;
    }
    const Vector& getUErr(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getUErr(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Velocity, "StateRep::getUErr(subsys)");
        return data->getSubsystem(subsys).uerr;
    }
    const Vector& getUDotErr(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getUDotErr(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Acceleration, "StateRep::getUDotErr(subsys)");
        return data->getSubsystem(subsys).udoterr;
    }
    const Vector& getMultipliers(SubsystemIndex subsys) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getMultipliers(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage::Acceleration, "StateRep::getMultipliers(subsys)");
        return data->getSubsystem(subsys).multipliers;
    }
    
    const Vector& getEventsByStage(SubsystemIndex subsys, Stage g) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getEventsByStage(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), g, "StateRep::getEventsByStage(subsys)");
        return data->getSubsystem(subsys).events[g];
    }
    
    Vector& updQErr(SubsystemIndex subsys) const {
        assert(data);
        checkCanModify(Stage::Position);
        checkCanModify(subsys);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::updQErr(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage(Stage::Position).prev(), "StateRep::updQErr(subsys)");
        return data->getSubsystem(subsys).qerr;
    }
    Vector& updUErr(SubsystemIndex subsys) const {
        assert(data);
        checkCanModify(Stage::Velocity);
        checkCanModify(subsys);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::updUErr(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage(Stage::Velocity).prev(), "StateRep::updUErr(subsys)");
        return data->getSubsystem(subsys).uerr;
    }
    Vector& updUDotErr(SubsystemIndex subsys) const {
        assert(data);
        checkCanModify(Stage::Velocity);
        checkCanModify(subsys);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::updUDotErr(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage(Stage::Acceleration).prev(), 
                            "StateRep::updUDotErr(subsys)");
        return data->getSubsystem(subsys).udoterr;
    }
    Vector& updMultipliers(SubsystemIndex subsys) const {
        assert(data);
        checkCanModify(subsys);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::updMultipliers(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), Stage(Stage::Acceleration).prev(), 
                            "StateRep::updMultipliers(subsys)");
        return data->getSubsystem(subsys).multipliers;
    }
    Vector& updEventsByStage(SubsystemIndex subsys, Stage g) const {
        assert(data);
        checkCanModify(subsys);
        checkCanModify(g);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::updEventsByStage(subsys)");
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), g.prev(), "StateRep::updEventsByStage(subsys)");
        return data->getSubsystem(subsys).events[g];
    }
    
        // Direct access to the global shared state and cache entries.
        // These are allocated once the System Stage is Stage::Model.
    
    const Real& getTime() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getTime()");
        return data->t;
    }
    
    const Vector& getY() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getY()");
        return data->y;
    }
    
    const Vector& getQ() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getQ()");
        return data->q;
    }
    
    const Vector& getU() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getU()");
        return data->u;
    }
    
    const Vector& getZ() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::getZ()");
        return data->z;
    }
    
    
    // You can call these as long as stage >= Model, but the
    // stage will be backed up if necessary to the indicated stage.
    Real& updTime() {  // Stage::Time-1
        assert(data);
        checkCanModify(Stage::Time);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::updTime()");
        invalidateAll(Stage::Time);
        return data->t;
    }
    
    Vector& updY() {    // Back to Stage::Position-1
        assert(data);
        checkCanModifyY();
        checkCanModifyAnySubsystem();
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::updY()");
        invalidateAll(Stage::Position);
        return data->y;
    }
    
    Vector& updQ() {    // Stage::Position-1
        assert(data);
        checkCanModify(Stage::Position);
        checkCanModifyAnySubsystem();
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::updQ()");
        invalidateAll(Stage::Position);
        return data->q;
    }
    
    Vector& updU() {     // Stage::Velocity-1
        assert(data);
        checkCanModify(Stage::Velocity);
        checkCanModifyAnySubsystem();
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::updU()");
        invalidateAll(Stage::Velocity);
        return data->u;
    }
    
    Vector& updZ() {     // Stage::Dynamics-1
        assert(data);
        checkCanModify(Stage::Dynamics);
        checkCanModifyAnySubsystem();
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Model, "StateRep::updZ()");
        invalidateAll(Stage::Dynamics);
        return data->z;
    }
    
    const Vector& getYDot() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Acceleration, "StateRep::getYDot()");
        return data->ydot;
    }
    
    const Vector& getQDot() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Velocity, "StateRep::getQDot()");
        return data->qdot;
    }
    
    const Vector& getZDot() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Dynamics, "StateRep::getZDot()");
        return data->zdot;
    }
    
    const Vector& getUDot() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Acceleration, "StateRep::getUDot()");
        return data->udot;
    }
    
    const Vector& getQDotDot() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Acceleration, "StateRep::getQDotDot()");
        return data->qdotdot;
    }
    
    
    Vector& updYDot() const {
        assert(data);
        checkCanModifyY();
        SimTK_STAGECHECK_GE(getSystemStage(), Stage(Stage::Acceleration).prev(), "StateRep::updYDot()");
        return data->ydot;
    }
    
    Vector& updQDot() const {
        assert(data);
        checkCanModify(Stage::Position);
        checkCanModifyAnySubsystem();
        SimTK_STAGECHECK_GE(getSystemStage(), Stage(Stage::Velocity).prev(), "StateRep::updQDot()");
        return data->qdot;
    }
    
    Vector& updUDot() const {
        assert(data);
        checkCanModify(Stage::Velocity);
        checkCanModifyAnySubsystem();
        SimTK_STAGECHECK_GE(getSystemStage(), Stage(Stage::Acceleration).prev(), "StateRep::updUDot()");
        return data->udot;
    }
    
    Vector& updZDot() const {
        assert(data);
        checkCanModify(Stage::Dynamics);
        checkCanModifyAnySubsystem();
        SimTK_STAGECHECK_GE(getSystemStage(), Stage(Stage::Dynamics).prev(), "StateRep::updZDot()");
        return data->zdot;
    }
    
    Vector& updQDotDot() const {
        assert(data);
        checkCanModify(Stage::Position);
        checkCanModifyAnySubsystem();
        SimTK_STAGECHECK_GE(getSystemStage(), Stage(Stage::Acceleration).prev(), "StateRep::updQDotDot()");
        return data->qdotdot;
    }
    
    
    const Vector& getYErr() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Velocity, "StateRep::getYErr()");
        return data->yerr;
    }
    
    const Vector& getQErr() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Position, "StateRep::getQErr()");
        return data->qerr;
    }
    const Vector& getUErr() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Velocity, "StateRep::getUErr()");
        return data->uerr;
    }
    const Vector& getUDotErr() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Acceleration, "StateRep::getUDotErr()");
        return data->udoterr;
    }
    const Vector& getMultipliers() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Acceleration, "StateRep::getMultipliers()");
        return data->multipliers;
    }
    
    Vector& updYErr() const {
        assert(data);
        checkCanModifyY();
        checkCanModifyAnySubsystem();
        SimTK_STAGECHECK_GE(getSystemStage(), Stage(Stage::Velocity).prev(), "StateRep::updYErr()");
        return data->yerr;
    }
    Vector& updQErr() const{
        assert(data);
        checkCanModify(Stage::Position);
        checkCanModifyAnySubsystem();
        SimTK_STAGECHECK_GE(getSystemStage(), Stage(Stage::Position).prev(), "StateRep::updQErr()");
        return data->qerr;
    }
    Vector& updUErr() const{
        assert(data);
        checkCanModify(Stage::Velocity);
        checkCanModifyAnySubsystem();
        SimTK_STAGECHECK_GE(getSystemStage(), Stage(Stage::Velocity).prev(), "StateRep::updUErr()");
        return data->uerr;
    }
    Vector& updUDotErr() const{
        assert(data);
        checkCanModify(Stage::Velocity);
        checkCanModifyAnySubsystem();
        SimTK_STAGECHECK_GE(getSystemStage(), Stage(Stage::Acceleration).prev(), 
                            "StateRep::updUDotErr()");
        return data->udoterr;
    }
    Vector& updMultipliers() const{
        assert(data);
        checkCanModifyAnySubsystem();
        SimTK_STAGECHECK_GE(getSystemStage(), Stage(Stage::Acceleration).prev(), 
                            "StateRep::updMultipliers()");
        return data->multipliers;
    }
    
    const Vector& getEvents() const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), Stage::Acceleration, "StateRep::getEvents()");
        return data->allEvents;
    }
    const Vector& getEventsByStage(Stage g) const {
        assert(data);
        SimTK_STAGECHECK_GE(getSystemStage(), g, "StateRep::getEventsByStage()");
        return data->events[g];
    }
    
    // These are mutable; hence 'const'.
    Vector& updEvents() const {
        assert(data);
        checkCanModifyAnySubsystem();
        SimTK_STAGECHECK_GE(getSystemStage(), Stage(Stage::Acceleration).prev(), "StateRep::updEvents()");
        return data->allEvents;
    }
    Vector& updEventsByStage(Stage g) const {
        assert(data);
        checkCanModify(g);
        SimTK_STAGECHECK_GE(getSystemStage(), g.prev(), "StateRep::updEventsByStage()");
        return data->events[g];
    }
    
    // You can access a Model stage variable any time, but don't access others
    // until you have realized the Model stage.
    const AbstractValue& 
    getDiscreteVariable(SubsystemIndex subsys, int index) const {
        const PerSubsystemInfo& ss = data->subsystems[subsys];
    
        SimTK_INDEXCHECK(0,index,(int)ss.discrete.size(),"StateRep::getDiscreteVariable()");
        const DiscreteVariable& dv = ss.discrete[index];
    
        if (dv.getStage() > Stage::Model) {
            SimTK_STAGECHECK_GE(getSubsystemStage(subsys), 
                Stage::Model, "StateRep::getDiscreteVariable()");
        }
    
        return dv.getValue();
    }
    
    // You can update a Model stage variable from Topology stage, but higher variables 
    // must wait until you have realized the Model stage. This always backs the 
    // stage up to one earlier than the variable's stage.
    AbstractValue& 
    updDiscreteVariable(SubsystemIndex subsys, int index) {
        checkCanModify(subsys);
        PerSubsystemInfo& ss = data->subsystems[subsys];
    
        SimTK_INDEXCHECK(0,index,(int)ss.discrete.size(),"StateRep::updDiscreteVariable()");
        DiscreteVariable& dv = ss.discrete[index];
    
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), 
            std::min(dv.getStage().prev(), Stage(Stage::Model)), 
            "StateRep::updDiscreteVariable()");
    
        invalidateAll(dv.getStage());
    
        return dv.updValue();
    }
    
    // Stage >= ce.stage
    const AbstractValue& 
    getCacheEntry(SubsystemIndex subsys, int index) const {
        const PerSubsystemInfo& ss = data->subsystems[subsys];
    
        SimTK_INDEXCHECK(0,index,(int)ss.cache.size(),"StateRep::getCacheEntry()");
        const CacheEntry& ce = ss.cache[index];
    
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), 
            ce.getStage(), "StateRep::getCacheEntry()");
    
        return ce.getValue();
    }
    
    // Stage >= ce.stage-1; does not change stage
    AbstractValue& 
    updCacheEntry(SubsystemIndex subsys, int index) const {
        checkCanModify(subsys);
        const PerSubsystemInfo& ss = data->subsystems[subsys];
    
        SimTK_INDEXCHECK(0,index,(int)ss.cache.size(),"StateRep::updCacheEntry()");
        CacheEntry& ce = ss.cache[index];
    
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), 
            ce.getStage().prev(), "StateRep::updCacheEntry()");
    
        return ce.updValue();
    }
    
    EnumerationSet<Stage> getRestrictedStages() const {
        return restrictedStages;
    }

    set<SubsystemIndex> getRestrictedSubsystems() const {
        return restrictedSubsystems;
    }
    
    // Verify that a particular stage may be modified.
    void checkCanModify(Stage stage) const {
        SimTK_ASSERT1_ALWAYS(!restrictedStages.contains(stage),
                "Modification of state data for stage %s has been restricted.", stage.getName().c_str());
    }
    
    // Verify that a particular subsystem may be modified.
    void checkCanModify(SubsystemIndex subsystem) const {
        SimTK_ASSERT1_ALWAYS(restrictedSubsystems.find(subsystem) == restrictedSubsystems.end(),
                "Modification of state data for subsystem %d has been restricted.", (int) subsystem);
    }
    
    // Verify that this State permits all state variables to be modified.
    void checkCanModifyY() const {
        checkCanModify(Stage::Position);
        checkCanModify(Stage::Velocity);
        checkCanModify(Stage::Dynamics);
    }
    
    // Verify that this State permits unrestricted modifications.
    void checkCanModifyAnySubsystem() const {
        SimTK_ASSERT_ALWAYS(restrictedSubsystems.empty(),
                "Modification of state data has been restricted.");
    }
    
    String toString() const {
        String out;
        out += "<State>\n";
    
        out += "<Real name=time>" + String(data->t) + "</Real>\n";
    
        out += "<Vector name=q size=" + String(data->q.size()) + ">";
        if (data->q.size()) out += "\n";
        for (long i=0; i<data->q.size(); ++i)
            out += String(data->q[i]) + "\n";
        out += "</Vector>\n";
    
        out += "<Vector name=u size=" + String(data->u.size()) + ">";
        if (data->u.size()) out += "\n";
        for (long i=0; i<data->u.size(); ++i)
            out += String(data->u[i]) + "\n";
        out += "</Vector>\n";
    
        out += "<Vector name=z size=" + String(data->z.size()) + ">";
        if (data->z.size()) out += "\n";
        for (long i=0; i<data->z.size(); ++i)
            out += String(data->z[i]) + "\n";
        out += "</Vector>\n";
    
    
        for (SubsystemIndex ss(0); ss < (int)data->subsystems.size(); ++ss) {
            const PerSubsystemInfo& info = data->subsystems[ss];
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
    
    String cacheToString() const {
        String out;
        out += "<Cache>\n";
        out += "<Stage>" + getSystemStage().getName() + "</Stage>\n";
    
        for (SubsystemIndex ss(0); ss < (int)data->subsystems.size(); ++ss) {
            const PerSubsystemInfo& info = data->subsystems[ss];
            out += "<Subsystem index=" + String(ss) + " name=" + info.name 
                + " version=" + info.version + ">\n";
            out += "  <Stage>" + info.currentStage.getName() + "</Stage>\n";
    
            out += "  <DISCRETE CACHE TODO>\n";
    
            out += "  <Vector name=qdot size=" + String(info.qdot.size()) + ">\n";
            out += "  <Vector name=udot size=" + String(info.udot.size()) + ">\n";
            out += "  <Vector name=zdot size=" + String(info.zdot.size()) + ">\n";
            out += "  <Vector name=qdotdot size=" + String(info.qdotdot.size()) + ">\n";
    
            out += "  <Vector name=qerr size=" + String(info.qerr.size()) + ">\n";
            out += "  <Vector name=uerr size=" + String(info.uerr.size()) + ">\n";
            out += "  <Vector name=udoterr size=" + String(info.udoterr.size()) + ">\n";
            out += "  <Vector name=multipliers size=" + String(info.multipliers.size()) + ">\n";
    
            for (int j=0; j<Stage::NValid; ++j) {
                out += "  <Vector name=events[";
                out += Stage::getValue(j).getName();
                out += "] size=" + String(info.events[j].size()) + ">\n";
            }
    
            out += "</Subsystem>\n";
        }
    
        out += "<Vector name=qdot size=" + String(data->qdot.size()) + ">";
        if (data->qdot.size()) out += "\n";
        for (long i=0; i<data->qdot.size(); ++i)
            out += String(data->qdot[i]) + "\n";
        out += "</Vector>\n";
    
        out += "<Vector name=udot size=" + String(data->udot.size()) + ">";
        if (data->udot.size()) out += "\n";
        for (long i=0; i<data->udot.size(); ++i)
            out += String(data->udot[i]) + "\n";
        out += "</Vector>\n";
    
        out += "<Vector name=zdot size=" + String(data->zdot.size()) + ">";
        if (data->zdot.size()) out += "\n";
        for (long i=0; i<data->zdot.size(); ++i)
            out += String(data->zdot[i]) + "\n";
        out += "</Vector>\n";
    
        out += "<Vector name=qdotdot size=" + String(data->qdotdot.size()) + ">";
        if (data->qdotdot.size()) out += "\n";
        for (long i=0; i<data->qdotdot.size(); ++i)
            out += String(data->qdotdot[i]) + "\n";
        out += "</Vector>\n";
    
        out += "<Vector name=qerr size=" + String(data->qerr.size()) + ">";
        if (data->qerr.size()) out += "\n";
        for (long i=0; i<data->qerr.size(); ++i)
            out += String(data->qerr[i]) + "\n";
        out += "</Vector>\n";
    
        out += "<Vector name=uerr size=" + String(data->uerr.size()) + ">";
        if (data->uerr.size()) out += "\n";
        for (long i=0; i<data->uerr.size(); ++i)
            out += String(data->uerr[i]) + "\n";
        out += "</Vector>\n";
    
        out += "<Vector name=udoterr size=" + String(data->udoterr.size()) + ">";
        if (data->udoterr.size()) out += "\n";
        for (long i=0; i<data->udoterr.size(); ++i)
            out += String(data->udoterr[i]) + "\n";
        out += "</Vector>\n";
    
        out += "<Vector name=multipliers size=" + String(data->multipliers.size()) + ">";
        if (data->multipliers.size()) out += "\n";
        for (long i=0; i<data->multipliers.size(); ++i)
            out += String(data->multipliers[i]) + "\n";
        out += "</Vector>\n";
    
        out += "</Cache>\n";
        return out;
    }
    StateData* data;
    EnumerationSet<Stage> restrictedStages;
    set<SubsystemIndex> restrictedSubsystems;
};

State::State() {
    rep = new StateRep();
}
State::~State() {
    delete rep;
}
State::State(const State& state) {
    rep = new StateRep(*state.rep);
}
void State::clear() {
    rep->clear();
}
void State::setNSubsystems(int i) {
    rep->setNSubsystems(i);
}
void State::initializeSubsystem(SubsystemIndex subsys, const String& name, const String& version) {
    rep->initializeSubsystem(subsys, name, version);
}
State& State::operator=(const State& state) {
    *rep = *state.rep;
    return *this;
}
int State::addSubsystem(const String& name, const String& version) {
    return rep->addSubsystem(name, version);
}
int State::getNSubsystems() const {
    return rep->getNSubsystems();
}
const String& State::getSubsystemName(SubsystemIndex subsys) const {
    return rep->getSubsystemName(subsys);
}
const String& State::getSubsystemVersion(SubsystemIndex subsys) const {
    return rep->getSubsystemVersion(subsys);
}
const Stage& State::getSubsystemStage(SubsystemIndex subsys) const {
    return rep->getSubsystemStage(subsys);
}
const Stage& State::getSystemStage() const {
    return rep->getSystemStage();
}
void State::invalidateAll(Stage stage) const {
    rep->invalidateAll(stage);
}
void State::advanceSubsystemToStage(SubsystemIndex subsys, Stage stage) const {
    rep->advanceSubsystemToStage(subsys, stage);
}
void State::advanceSystemToStage(Stage stage) const {
    rep->advanceSystemToStage(stage);
}
int State::allocateQ(SubsystemIndex subsys, const Vector& qInit) {
    return rep->allocateQ(subsys, qInit);
}
int State::allocateU(SubsystemIndex subsys, const Vector& uInit) {
    return rep->allocateU(subsys, uInit);
}
int State::allocateZ(SubsystemIndex subsys, const Vector& zInit) {
    return rep->allocateZ(subsys, zInit);
}
int State::allocateQErr(SubsystemIndex subsys, int nqerr) {
    return rep->allocateQErr(subsys, nqerr);
}
int State::allocateUErr(SubsystemIndex subsys, int nuerr) {
    return rep->allocateUErr(subsys, nuerr);
}
int State::allocateUDotErr(SubsystemIndex subsys, int nudoterr) {
    return rep->allocateUDotErr(subsys, nudoterr);
}
int State::allocateEvent(SubsystemIndex subsys, Stage stage, int nevent) {
    return rep->allocateEvent(subsys, stage, nevent);
}
int State::allocateDiscreteVariable(SubsystemIndex subsys, Stage stage, AbstractValue* v) {
    return rep->allocateDiscreteVariable(subsys, stage, v);
}
int State::allocateCacheEntry(SubsystemIndex subsys, Stage stage, AbstractValue* v) {
    return rep->allocateCacheEntry(subsys, stage, v);
}
int State::getNY() const {
    return rep->getNY();
}
int State::getQStart() const {
    return rep->getQStart();
}
int State::getNQ() const {
    return rep->getNQ();
}
int State::getUStart() const {
    return rep->getUStart();
}
int State::getNU() const {
    return rep->getNU();
}
int State::getZStart() const {
    return rep->getZStart();
}
int State::getNZ() const {
    return rep->getNZ();
}
int State::getNYErr() const {
    return rep->getNYErr();
}
int State::getQErrStart() const {
    return rep->getQErrStart();
}
int State::getNQErr() const {
    return rep->getNQErr();
}
int State::getUErrStart() const {
    return rep->getUErrStart();
}
int State::getNUErr() const {
    return rep->getNUErr();
}
int State::getNUDotErr() const {
    return rep->getNUDotErr();
}
int State::getNMultipliers() const {
    return getNUDotErr();
}
int State::getQStart(SubsystemIndex subsys) const {
    return rep->getQStart(subsys);
}
int State::getNQ(SubsystemIndex subsys) const {
    return rep->getNQ(subsys);
}
int State::getUStart(SubsystemIndex subsys) const {
    return rep->getUStart(subsys);
}
int State::getNU(SubsystemIndex subsys) const {
    return rep->getNU(subsys);
}
int State::getZStart(SubsystemIndex subsys) const {
    return rep->getZStart(subsys);
}
int State::getNZ(SubsystemIndex subsys) const {
    return rep->getNZ(subsys);
}
int State::getQErrStart(SubsystemIndex subsys) const {
    return rep->getQErrStart(subsys);
}
int State::getNQErr(SubsystemIndex subsys) const {
    return rep->getNQErr(subsys);
}
int State::getUErrStart(SubsystemIndex subsys) const {
    return rep->getUErrStart(subsys);
}
int State::getNUErr(SubsystemIndex subsys) const {
    return rep->getNUErr(subsys);
}
int State::getUDotErrStart(SubsystemIndex subsys) const {
    return rep->getUDotErrStart(subsys);
}
int State::getNUDotErr(SubsystemIndex subsys) const {
    return rep->getNUDotErr(subsys);
}
int State::getMultipliersStart(SubsystemIndex i) const {
    return getUDotErrStart(i);
}
int State::getNMultipliers(SubsystemIndex i) const {
    return getNUDotErr(i);
}
int State::getNEvents() const {
    return rep->getNEvents();
}
int State::getEventStartByStage(Stage stage) const {
    return rep->getEventStartByStage(stage);
}
int State::getNEventsByStage(Stage stage) const {
    return rep->getNEventsByStage(stage);
}
int State::getEventStartByStage(SubsystemIndex subsys, Stage stage) const {
    return rep->getEventStartByStage(subsys, stage);
}
int State::getNEventsByStage(SubsystemIndex subsys, Stage stage) const {
    return rep->getNEventsByStage(subsys, stage);
}
const Vector& State::getEvents() const {
    return rep->getEvents();
}
const Vector& State::getEventsByStage(Stage stage) const {
    return rep->getEventsByStage(stage);
}
const Vector& State::getEventsByStage(SubsystemIndex subsys, Stage stage) const {
    return rep->getEventsByStage(subsys, stage);
}
Vector& State::updEvents() const {
    return rep->updEvents();
}
Vector& State::updEventsByStage(Stage stage) const {
    return rep->updEventsByStage(stage);
}
Vector& State::updEventsByStage(SubsystemIndex subsys, Stage stage) const {
    return rep->updEventsByStage(subsys, stage);
}
const Vector& State::getQ(SubsystemIndex subsys) const {
    return rep->getQ(subsys);
}
const Vector& State::getU(SubsystemIndex subsys) const {
    return rep->getU(subsys);
}
const Vector& State::getZ(SubsystemIndex subsys) const {
    return rep->getZ(subsys);
}
Vector& State::updQ(SubsystemIndex subsys) {
    return rep->updQ(subsys);
}
Vector& State::updU(SubsystemIndex subsys) {
    return rep->updU(subsys);
}
Vector& State::updZ(SubsystemIndex subsys) {
    return rep->updZ(subsys);
}
const Vector& State::getQDot(SubsystemIndex subsys) const {
    return rep->getQDot(subsys);
}
const Vector& State::getUDot(SubsystemIndex subsys) const {
    return rep->getUDot(subsys);
}
const Vector& State::getZDot(SubsystemIndex subsys) const {
    return rep->getZDot(subsys);
}
const Vector& State::getQDotDot(SubsystemIndex subsys) const {
    return rep->getQDotDot(subsys);
}
Vector& State::updQDot(SubsystemIndex subsys) const {
    return rep->updQDot(subsys);
}
Vector& State::updUDot(SubsystemIndex subsys) const {
    return rep->updUDot(subsys);
}
Vector& State::updZDot(SubsystemIndex subsys) const {
    return rep->updZDot(subsys);
}
Vector& State::updQDotDot(SubsystemIndex subsys) const {
    return rep->updQDotDot(subsys);
}
const Vector& State::getQErr(SubsystemIndex subsys) const {
    return rep->getQErr(subsys);
}
const Vector& State::getUErr(SubsystemIndex subsys) const {
    return rep->getUErr(subsys);
}
const Vector& State::getUDotErr(SubsystemIndex subsys) const {
    return rep->getUDotErr(subsys);
}
const Vector& State::getMultipliers(SubsystemIndex subsys) const {
    return rep->getMultipliers(subsys);
}
Vector& State::updQErr(SubsystemIndex subsys) const {
    return rep->updQErr(subsys);
}
Vector& State::updUErr(SubsystemIndex subsys) const {
    return rep->updUErr(subsys);
}
Vector& State::updUDotErr(SubsystemIndex subsys) const {
    return rep->updUDotErr(subsys);
}
Vector& State::updMultipliers(SubsystemIndex subsys) const {
    return rep->updMultipliers(subsys);
}
const Real& State::getTime() const {
    return rep->getTime();
}
const Vector& State::getY() const {
    return rep->getY();
}
const Vector& State::getQ() const {
    return rep->getQ();
}
const Vector& State::getU() const {
    return rep->getU();
}
const Vector& State::getZ() const {
    return rep->getZ();
}
Real& State::updTime() {
    return rep->updTime();
}
Vector& State::updY() {
    return rep->updY();
}
void State::setTime(Real t) {
    updTime() = t;
}
void State::setY(const Vector& y) {
    updY() = y;
}
Vector& State::updQ() {
    return rep->updQ();
}
Vector& State::updU() {
    return rep->updU();
}
Vector& State::updZ() {
    return rep->updZ();
}
void State::setQ(const Vector& q) {
    updQ() = q;
}
void State::setU(const Vector& u) {
    updU() = u;
}
void State::setZ(const Vector& z) {
    updZ() = z;
}
const Vector& State::getYDot() const {
    return rep->getYDot();
}
const Vector& State::getQDot() const {
    return rep->getQDot();
}
const Vector& State::getZDot() const {
    return rep->getZDot();
}
const Vector& State::getUDot() const {
    return rep->getUDot();
}
const Vector& State::getQDotDot() const {
    return rep->getQDotDot();
}
Vector& State::updYDot() const {
    return rep->updYDot();
}
Vector& State::updQDot() const {
    return rep->updQDot();
}
Vector& State::updZDot() const {
    return rep->updZDot();
}
Vector& State::updUDot() const {
    return rep->updUDot();
}
Vector& State::updQDotDot() const {
    return rep->updQDotDot();
}
const Vector& State::getYErr() const {
    return rep->getYErr();
}
const Vector& State::getQErr() const {
    return rep->getQErr();
}
const Vector& State::getUErr() const {
    return rep->getUErr();
}
const Vector& State::getUDotErr() const {
    return rep->getUDotErr();
}
const Vector& State::getMultipliers() const {
    return rep->getMultipliers();
}
Vector& State::updYErr() const {
    return rep->updYErr();
}
Vector& State::updQErr() const {
    return rep->updQErr();
}
Vector& State::updUErr() const {
    return rep->updUErr();
}
Vector& State::updUDotErr() const {
    return rep->updUDotErr();
}
Vector& State::updMultipliers() const {
    return rep->updMultipliers();
}
const AbstractValue& State::getDiscreteVariable(SubsystemIndex subsys, int index) const {
    return rep->getDiscreteVariable(subsys, index);
}
AbstractValue& State::updDiscreteVariable(SubsystemIndex subsys, int index) {
    return rep->updDiscreteVariable(subsys, index);
}
void State::setDiscreteVariable(SubsystemIndex i, int index, const AbstractValue& v) {
    updDiscreteVariable(i,index) = v;
}
const AbstractValue& State::getCacheEntry(SubsystemIndex subsys, int index) const {
    return rep->getCacheEntry(subsys, index);
}
AbstractValue& State::updCacheEntry(SubsystemIndex subsys, int index) const {
    return rep->updCacheEntry(subsys, index);
}
void State::createRestrictedState(State& restrictedState, EnumerationSet<Stage> restrictedStages, std::set<SubsystemIndex> restrictedSubsystems) {
    restrictedStages |= getRestrictedStages();
    std::set<SubsystemIndex> currentSubsystems = getRestrictedSubsystems();
    restrictedSubsystems.insert(currentSubsystems.begin(), currentSubsystems.end());
    delete restrictedState.rep;
    restrictedState.rep = new StateRep(*rep, restrictedStages, restrictedSubsystems);
}
EnumerationSet<Stage> State::getRestrictedStages() const {
    return rep->getRestrictedStages();
}
set<SubsystemIndex> State::getRestrictedSubsystems() const {
    return rep->getRestrictedSubsystems();
}
String State::toString() const {
    return rep->toString();
}
String State::cacheToString() const {
    return rep->cacheToString();
}

std::ostream& 
operator<<(std::ostream& o, const State& s) {
    o << "STATE:" << std::endl;
    o << s.toString() << std::endl;
    o << "CACHE:" << std::endl;
    return o << s.cacheToString() << std::endl;
}

} // namespace SimTK

