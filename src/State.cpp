/* Copyright (c) 2005-6 Stanford University and Michael Sherman.
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
#include "Simmatrix.h"
#include "simbody/internal/common.h"
#include "simbody/internal/State.h"

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
        assert(g == Stage::Built || Stage::isInRuntimeRange(g));
        assert(vp);
    }

    DiscreteVariableRep* clone() const {return new DiscreteVariableRep(*this);}

    bool isValid() const {
        return (stage==Stage::Built || Stage::isInRuntimeRange(stage)) && value && myHandle; 
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
        copyFrom(src, Stage::Modeled);
    }

    PerSubsystemInfo& operator=(const PerSubsystemInfo& src) {
        if (&src != this) {
            copyFrom(src, Stage::Modeled);
        }
        return *this;
    }

    // Back up to the stage just before g if this subsystem thinks
    // it is already at g or beyond. Note that we may be backing up
    // over many stages here. Careful: invalidating the stage
    // for a subsystem must also invalidate the same stage for all
    // the other subsystems and the system as a whole but we don't
    // take care of that here. Also, you can't invalidate Stage::Allocated.
    void invalidateStageJustThisSubsystem(Stage g) {
        assert(g > Stage::Allocated);
        restoreToStage(g.prev());
    }

    // Advance from stage g-1 to stage g. This is called at the end
    // of realize(g). You can't use this to "advance" to Stage::Allocated.
    // It is a fatal error if the current stage isn't g-1.
    void advanceToStage(Stage g) {
        assert(g > Stage::Allocated);
        assert(currentStage == g.prev());
        // Record data needed to back up to this stage later.
        if (g == Stage::Built)
            nDiscreteWhenBuilt = discrete.size();
        cacheSize[g] = cache.size();
        currentStage = g;
    }

    String name;
    String version;
    std::vector<int> qs, us, zs;    // indices into q,u,z and qdot, etc.

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
        nDiscreteWhenBuilt = -1;
        for (int i=0; i < Stage::NValid; ++i)
            cacheSize[i] = -1;
        cacheSize[Stage::Allocated] = 0;
        currentStage = Stage::Allocated;
    }

    // Set all the allocation sizes to zero, leaving name and version alone.
    // The stage is set back to "Allocated".
    void restoreToAllocatedStage() {
        if (currentStage > Stage::Allocated) {
            qs.clear(); us.clear(); zs.clear();
            discrete.clear(); cache.clear();
            initialize();
        }
    }

    // Put this Subsystem back to the way it was just after realize(Built).
    // That means: no shared global state vars, discrete vars and cache
    // only as they were after realize(Built).
    // The stage is set back to "Built". Nothing happens if the stage
    // is already at Built or below.
    void restoreToBuiltStage() {
        if (currentStage <= Stage::Built)
            return;
        qs.clear(); us.clear(); zs.clear(); // nuke the globals
        discrete.resize(nDiscreteWhenBuilt);
        restoreCacheToStage(Stage::Built);
        currentStage = Stage::Built;
    }

    // Restore this subsystem to the way it last was at realize(Stage).
    void restoreToStage(Stage g) {
        if (g==Stage::Allocated) {restoreToAllocatedStage(); return;}
        if (g==Stage::Built)     {restoreToBuiltStage();     return;}
        if (currentStage <= g) return;

        // State variables remain unchanged since they are all allocated
        // after realize(Modeled). Cache gets shrunk to the length it
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
    // was after being realized to stage maxStage. If maxStage >= Modeled then
    // all the state variables we be copied, but only cached computations
    // up through maxStage come through.
    void copyFrom(const PerSubsystemInfo& src, Stage maxStage) {
        const Stage targetStage = std::min<Stage>(src.currentStage, maxStage);
        name     = src.name;
        version  = src.version;

        if (targetStage < Stage::Built) {
            restoreToAllocatedStage();  // don't copy anything
            return;
        }

        // At "Built" stage we need to copy all the private state
        // variables that were present at realize(Built) (those will
        // mostly be modeling variables). There can't
        // be any global ones since those aren't allocated until
        // Modeled stage. We also need to copy any cache results
        // that were available after realize(Built).
        if (targetStage == Stage::Built) {
            restoreToAllocatedStage();
            discrete.resize(src.nDiscreteWhenBuilt);
            for (int i=0; i<src.nDiscreteWhenBuilt; ++i)
                discrete[i] = src.discrete[i];
            nDiscreteWhenBuilt = discrete.size();
            // don't copy any global shared resources

            cache.resize(src.cacheSize[Stage::Built]);
            for (int i=0; i<(int)cache.size(); ++i)
                cache[i] = src.cache[i];
            cacheSize[Stage::Built] = cache.size();
            currentStage = Stage::Built;
            return;
        }

        // This is the general case where Stage > Built.

        // Copy *all* state variables since no more can be allocated
        // after Modeled stage.
        discrete = src.discrete;
        qs = src.qs; us = src.us; zs = src.zs;

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
      : t(CNT<Real>::getNaN()), systemStage(Stage::Allocated), 
        subsystems(1), myHandle(0) 
    { 
    }

    const Stage& getSystemStage() const {
        return systemStage;
    }
    Stage& updSystemStage() const {
        return systemStage; // mutable
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
    StateRep(const StateRep& src) : myHandle(0) {
        t = src.t; q = src.q; u = src.u; z = src.z;
        systemStage = std::min<Stage>(src.systemStage, Stage::Modeled);
        subsystems = src.subsystems;
    }

    StateRep& operator=(const StateRep& src) {
        if (&src == this) return *this;
        t = src.t; q = src.q; u = src.u; z = src.z;
        systemStage = std::min<Stage>(src.systemStage, Stage::Modeled);
        subsystems = src.subsystems;
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
    // here if we back up earlier than Modeled) but we don't
    // take care of that here. Also, you can't invalidate Stage::Allocated.
    void invalidateJustSystemStage(Stage g) {
        assert(g > Stage::Allocated);
        if (systemStage >= g) {
            if (g <= Stage::Modeled) {
                t = CNT<Real>::getNaN();
                q.resize(0); u.resize(0); z.resize(0);
            }
            systemStage = g.prev();
        }
    }

    // Advance the System stage from g-1 to g. It is a fatal error if
    // we're not already at g-1, and you can't advance to Stage::Allocated.
    // Also, you can't advance the system to g unless ALL subsystems have
    // already gotten there.
    void advanceSystemToStage(Stage g) {
        assert(g > Stage::Allocated);
        assert(systemStage == g.prev());
        assert(allSystemsAtLeastAtStage(g));

        if (g == Stage::Modeled) {
            // We know the shared state pool sizes now. Allocate the
            // matching shared cache pools.
            // TODO: q.setLockRows(true);
            //u.setLockRows(true);
            //z.setLockRows(true);
            qdot.resize(q.size());    //qdot.setLockRows(true);
            qdotdot.resize(q.size()); //qdotdot.setLockRows(true);
            udot.resize(u.size());    //udot.setLockRows(true);
            zdot.resize(z.size());    //zdot.setLockRows(true);
        }

        systemStage = g;
    }

    void         setMyHandle(State& s) {myHandle = &s;}
    const State& getMyHandle() const   {assert(myHandle); return *myHandle;}
    State&       updMyHandle()         {assert(myHandle); return *myHandle;}
private:
    friend class State;

        // Shared global resource State variables //

    Real            t; // Stage::Timed (time)
    Vector          q; // Stage::Configured continuous variables
    Vector          u; // Stage::Moving continuous variables
    Vector          z; // Stage::Dynamics continuous variables


        // Shared global resource Cache entries //
    mutable Stage   systemStage;

    mutable Vector  qdot;       // Stage::Moving
    mutable Vector  qdotdot;    // Stage::Reacting
    mutable Vector  udot;       // Stage::Reacting
    mutable Vector  zdot;       // Stage::Reacting

        // Subsystem support //

    // Subsystem 0 (always present) is for the System as a whole. Its name
    // and version are the System name and version.
    std::vector<PerSubsystemInfo> subsystems;

    // Return true only if all subsystems have been realized to at least Stage g.
    bool allSystemsAtLeastAtStage(Stage g) const {
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
State::State()
  : rep(new StateRep()) {
    rep->setMyHandle(*this);
}

State::State(const String& name, const String& version)
  : rep(new StateRep()) {
    rep->setMyHandle(*this);
    rep->subsystems[0].name = name;
    rep->subsystems[0].version = version;
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
    SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), Stage::Modeled, "State::allocateQ()");

    // Allocate global Q
    const int nxt = rep->q.size();
    rep->q.resize(nxt+qInit.size());
    rep->q(nxt, qInit.size()) = qInit;

    // Map to local subsystem Q
    rep->subsystems[subsys].qs.push_back(nxt);
    return (int)rep->subsystems[subsys].qs.size() - 1;
}

int State::allocateU(int subsys, const Vector& uInit) {
    SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), Stage::Modeled, "State::allocateU()");

    // Allocate global U
    const int nxt = rep->u.size();
    rep->u.resize(nxt+uInit.size());
    rep->u(nxt, uInit.size()) = uInit;

    // Map to local subsystem U
    rep->subsystems[subsys].us.push_back(nxt);
    return (int)rep->subsystems[subsys].us.size() - 1;
}
int State::allocateZ(int subsys, const Vector& zInit) {
    SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), Stage::Modeled, "State::allocateZ()");

    // Allocate global Z
    const int nxt = rep->z.size();
    rep->z.resize(nxt+zInit.size());
    rep->z(nxt, zInit.size()) = zInit;

    // Map to local subsystem Z
    rep->subsystems[subsys].zs.push_back(nxt);
    return (int)rep->subsystems[subsys].zs.size() - 1;
}

// Construction- and Modeling-stage State variables can only be added during construction; that is,
// while stage <= Built. Other entries can be added while stage < Modeled.
int State::allocateDiscreteVariable(int subsys, Stage g, AbstractValue* vp) {
    SimTK_STAGECHECK_RANGE_ALWAYS(Stage(Stage::LowestRuntime).prev(), g, Stage::HighestRuntime, 
        "State::allocateDiscreteVariable()");

    const Stage maxAcceptable = (g <= Stage::Modeled ? Stage::Allocated : Stage::Built);
    SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), 
        maxAcceptable.next(), "State::allocateDiscreteVariable()");

    PerSubsystemInfo& ss = rep->subsystems[subsys];
    const int nxt = ss.discrete.size();
    ss.discrete.push_back(DiscreteVariable(g,vp));
    return nxt;
}

// Cache entries can be allocated while stage < Modeled, even if they are Modeled-stage entries.
int State::allocateCacheEntry(int subsys, Stage g, AbstractValue* vp) {
    SimTK_STAGECHECK_RANGE_ALWAYS(Stage(Stage::LowestRuntime).prev(), g, Stage::HighestRuntime, 
        "State::allocateCacheEntry()");
    SimTK_STAGECHECK_LT_ALWAYS(getSubsystemStage(subsys), 
        Stage::Modeled, "State::allocateCacheEntry()");

    PerSubsystemInfo& ss = rep->subsystems[subsys];
    const int nxt = ss.cache.size();
    ss.cache.push_back(CacheEntry(g,vp));
    return nxt;
}

// You can call these as long as stage >= Modeled.
const Real&
State::getTime() const {
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Modeled, "State::getTime()");
    return rep->t;
}

const Vector&
State::getQ() const {
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Modeled, "State::getQ()");
    return rep->q;
}

const Vector&
State::getU() const {
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Modeled, "State::getU()");
    return rep->u;
}

const Vector&
State::getZ() const {
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Modeled, "State::getZ()");
    return rep->z;
}


// You can call these as long as stage >= Modeled, but the
// stage will be backed up if necessary to the indicated stage.
Real&
State::updTime() {  // Stage::Timed-1
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Modeled, "State::updTime()");
    invalidateAll(Stage::Timed);
    return rep->t;
}

Vector&
State::updQ() {    // Stage::Configured-1
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Modeled, "State::updQ()");
    invalidateAll(Stage::Configured);
    return rep->q;
}

Vector&
State::updU() {     // Stage::Moving-1
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Modeled, "State::updU()");
    invalidateAll(Stage::Moving);
    return rep->u;
}

Vector&
State::updZ() {     // Stage::Dynamics-1
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Modeled, "State::updZ()");
    invalidateAll(Stage::Dynamics);
    return rep->z;
}


const Vector&
State::getQDot() const {
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Moving, "State::getQDot()");
    return rep->qdot;
}

const Vector&
State::getZDot() const {
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Dynamics, "State::getZDot()");
    return rep->zdot;
}

const Vector&
State::getUDot() const {
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Reacting, "State::getUDot()");
    return rep->udot;
}

const Vector&
State::getQDotDot() const {
    SimTK_STAGECHECK_GE(getSystemStage(), Stage::Reacting, "State::getQDotDot()");
    return rep->qdotdot;
}

Vector&
State::updQDot() const {
    SimTK_STAGECHECK_GE(getSystemStage(), Stage(Stage::Moving).prev(), "State::updQDot()");
    return rep->qdot;
}

Vector&
State::updZDot() const {
    SimTK_STAGECHECK_GE(getSystemStage(), Stage(Stage::Dynamics).prev(), "State::updZDot()");
    return rep->zdot;
}

Vector&
State::updUDot() const {
    SimTK_STAGECHECK_GE(getSystemStage(), Stage(Stage::Reacting).prev(), "State::updUDot()");
    return rep->udot;
}

Vector&
State::updQDotDot() const {
    SimTK_STAGECHECK_GE(getSystemStage(), Stage(Stage::Reacting).prev(), "State::updQDotDot()");
    return rep->qdotdot;
}

// You can access a Modeling variable any time, but don't access others
// until you have realized the Modeled stage.
const AbstractValue& 
State::getDiscreteVariable(int subsys, int index) const {
    const PerSubsystemInfo& ss = rep->subsystems[subsys];

    SimTK_INDEXCHECK(0,index,(int)ss.discrete.size(),"State::getDiscreteVariable()");
    const DiscreteVariable& dv = ss.discrete[index];

    if (dv.getStage() > Stage::Modeled) {
        SimTK_STAGECHECK_GE(getSubsystemStage(subsys), 
            Stage::Modeled, "State::getDiscreteVariable()");
    }

    return dv.getValue();
}

// You can update a Modeling variable from Built stage, but higher variables 
// must wait until you have realized the Modeled stage. This always backs the 
// stage up to one earlier than the variable's stage.
AbstractValue& 
State::updDiscreteVariable(int subsys, int index) {
    PerSubsystemInfo& ss = rep->subsystems[subsys];

    SimTK_INDEXCHECK(0,index,(int)ss.discrete.size(),"State::updDiscreteVariable()");
    DiscreteVariable& dv = ss.discrete[index];

    SimTK_STAGECHECK_GE(getSubsystemStage(subsys), 
        std::min(dv.getStage().prev(), Stage(Stage::Modeled)), 
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

static String indexToString(const std::vector<int>& ix, const String& name) {
    String out;
    out += "  <Indices name=qs size=" + String(ix.size()) + ">";
    if (ix.size()) out += "\n";

    for (int i=0; i<(int)ix.size(); ++i)
        out += "    " + String(ix[i]) + "\n";
    out += "  </Indices>\n";
    return out;
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

        out += indexToString(info.qs, "qs");
        out += indexToString(info.us, "us");
        out += indexToString(info.zs, "zs");

        out += "</Subsystem>\n";
    }

    out += "</State>\n";
    return out;
}

String State::cacheToString() const {
    String out;
    out += "<Cache>\n";

    for (int ss=0; ss < (int)rep->subsystems.size(); ++ss) {
        const PerSubsystemInfo& info = rep->subsystems[ss];
        out += "<Subsystem index=" + String(ss) + " name=" + info.name 
            + " version=" + info.version + ">\n";
        out += "  <Stage>" + info.currentStage.name() + "</Stage>\n";


        out += "  <DISCRETE CACHE TODO>\n";
        out += "</Subsystem>\n";
    }

    out += "<Vector name=qdot size=" + String(rep->qdot.size()) + ">";
    if (rep->qdot.size()) out += "\n";
    for (long i=0; i<rep->qdot.size(); ++i)
        out += String(rep->qdot[i]) + "\n";
    out += "</Vector>\n";

    out += "<Vector name=qdotdot size=" + String(rep->qdotdot.size()) + ">";
    if (rep->qdotdot.size()) out += "\n";
    for (long i=0; i<rep->qdotdot.size(); ++i)
        out += String(rep->qdotdot[i]) + "\n";
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


    out += "</Cache>\n";
    return out;
}

} // namespace SimTK

