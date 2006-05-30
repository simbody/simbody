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

class StateRep {
public:
    StateRep() : t(CNT<Real>::getNaN()), subsystems(1), myHandle(0) { 
    }

    const Stage& getCurrentStage(int subsystem) const {
        return subsystems[subsystem].currentStage;
    }
    Stage& updCurrentStage(int subsystem) const {
        return subsystems[subsystem].currentStage;
    }

    // We'll do the copy constructor and assignment explicitly here
    // to get tight control over what's allowed, and to make sure
    // we don't copy the handle pointer.
    StateRep(const StateRep& src) : myHandle(0) {
        t = src.t; q = src.q; u = src.u; z = src.z;
        discrete   = src.discrete;
        cache      = src.cache; // TODO: shouldn't copy cache(?)
        subsystems = src.subsystems;
    }

    StateRep& operator=(const StateRep& src) {
        if (&src == this) return *this;
        t = src.t; q = src.q; u = src.u; z = src.z;
        discrete   = src.discrete;
        qdot       = src.qdot; qdotdot = src.qdotdot; udot = src.udot; zdot = src.zdot;
        cache      = src.cache;
        subsystems = src.subsystems;
        // don't mess with the handle pointer!
        return *this;
    }

    // Copies everything but the handle pointer.
    StateRep* clone() const {return new StateRep(*this);}

    // Done realizing the model. All state variable and cache allocations
    // are now done.
    void finishModeling() {

        // TODO: q.setLockRows(true);
        //u.setLockRows(true);
        //z.setLockRows(true);
        qdot.resize(q.size());    //qdot.setLockRows(true);
        qdotdot.resize(q.size()); //qdotdot.setLockRows(true);
        udot.resize(u.size());    //udot.setLockRows(true);
        zdot.resize(z.size());    //zdot.setLockRows(true);
    }


    void         setMyHandle(State& s) {myHandle = &s;}
    const State& getMyHandle() const   {assert(myHandle); return *myHandle;}
    State&       updMyHandle()         {assert(myHandle); return *myHandle;}

private:
    friend class State;

        // State variables //

    Real            t; // Stage::Timed (time)
    Vector          q; // Stage::Configured continuous variables
    Vector          u; // Stage::Moving continuous variables
    Vector          z; // Stage::Dynamics continuous variables

    // Each of these discrete entries has its own Stage.
    StableArray<DiscreteVariable> discrete;

        // Cache //

    mutable Vector  qdot;       // Stage::Moving
    mutable Vector  qdotdot;    // Stage::Reacting
    mutable Vector  udot;       // Stage::Reacting
    mutable Vector  zdot;       // Stage::Reacting

    // Each of these discrete entries has its own Stage.
    mutable StableArray<CacheEntry> cache;

        // Subsystem support //

    struct PerSubsystemInfo {
        PerSubsystemInfo() : currentStage(Stage::Allocated) { }
        PerSubsystemInfo(const String& n, const String& v) 
            : name(n), version(v), currentStage(Stage::Allocated) { }
        String name;
        String version;
        std::vector<int> qs, us, zs;    // indices into q,u,z and qdot, etc.
        std::vector<int> discreteVars;
        std::vector<int> cacheEntries;
        mutable Stage    currentStage;
    };

    // Subsystem 0 (always present) is the System as a whole.
    std::vector<PerSubsystemInfo> subsystems;

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
        StateRep::PerSubsystemInfo(name,version));
    return (int)rep->subsystems.size() - 1;
}

int State::getNSubsystems() const {return (int)rep->subsystems.size();}

const String& State::getSubsystemName(int subsys) const {
    return rep->subsystems[subsys].name;
}
const String& State::getSubsystemVersion(int subsys) const {
    return rep->subsystems[subsys].version;
}

const Stage& State::getStage(int subsys) const {
    SimTK_ASSERT(rep, "State::getStage(): no rep"); // can't happen(?)
    return rep->getCurrentStage(subsys);
}

// Make sure the stage is no higher than g. Currently we have to reduce
// the currentState for *all* subsystems to g, regardless of which one
// is registering the change. TODO: this should be more selective.
void State::invalidateStage(int /*subsys*/, Stage g) const {
    SimTK_ASSERT(rep, "State::invalidateStage(): no rep");

    for (int i=0; i<(int)rep->subsystems.size(); ++i)
        if (rep->getCurrentStage(i) >= g)
            rep->updCurrentStage(i) = g.prev();    // mutable

    //TODO: special handling should be done when invalidating "Modeled" or "Built".
}

// Move the stage for a particular subsystem from g-1 to g. No other subsystems
// are affected.
void State::advanceToStage(int subsys, Stage g) const {
    SimTK_ASSERT(rep, "State::advanceToStage(): no rep");

    // TODO: need a real error here, not an assert
    SimTK_ASSERT2(getStage(subsys) == g.prev(),
        "State::advanceToStage(%s) called but current stage is %s",
        g.name().c_str(), getStage(subsys).name().c_str());

    if (g == Stage::Modeled)
        rep->finishModeling();

    rep->updCurrentStage(subsys) = g; // mutable
}

// We don't expect State entry allocations to be performance critical so
// we'll keep error checking on even in Release mode.

int State::allocateQ(int subsys, const Vector& qInit) {
    SimTK_STAGECHECK_LT_ALWAYS(getStage(subsys), Stage::Modeled, "State::allocateQ()");

    // Allocate global Q
    const int nxt = rep->q.size();
    rep->q.resize(nxt+qInit.size());
    rep->q(nxt, qInit.size()) = qInit;

    // Map to local subsystem Q
    rep->subsystems[subsys].qs.push_back(nxt);
    return (int)rep->subsystems[subsys].qs.size() - 1;
}

int State::allocateU(int subsys, const Vector& uInit) {
    SimTK_STAGECHECK_LT_ALWAYS(getStage(subsys), Stage::Modeled, "State::allocateU()");

    // Allocate global U
    const int nxt = rep->u.size();
    rep->u.resize(nxt+uInit.size());
    rep->u(nxt, uInit.size()) = uInit;

    // Map to local subsystem U
    rep->subsystems[subsys].us.push_back(nxt);
    return (int)rep->subsystems[subsys].us.size() - 1;
}
int State::allocateZ(int subsys, const Vector& zInit) {
    SimTK_STAGECHECK_LT_ALWAYS(getStage(subsys), Stage::Modeled, "State::allocateZ()");

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
    SimTK_STAGECHECK_LT_ALWAYS(getStage(subsys), maxAcceptable.next(), "State::allocateDiscreteVariable()");

    // Allocate global DiscreteVariable
    const int nxt = rep->discrete.size();
    rep->discrete.push_back(DiscreteVariable(g,vp));

    // Map to local subsystem DiscreteVariable
    rep->subsystems[subsys].discreteVars.push_back(nxt);
    return (int)rep->subsystems[subsys].discreteVars.size() - 1;
}

// Cache entries can be allocated while stage < Modeled, even if they are Modeled-stage entries.
int State::allocateCacheEntry(int subsys, Stage g, AbstractValue* vp) {
    SimTK_STAGECHECK_RANGE_ALWAYS(Stage(Stage::LowestRuntime).prev(), g, Stage::HighestRuntime, 
        "State::allocateCacheEntry()");
    SimTK_STAGECHECK_LT_ALWAYS(getStage(subsys), Stage::Modeled, "State::allocateCacheEntry()");

    // Allocate global CacheEntry
    const int nxt = rep->cache.size();
    rep->cache.push_back(CacheEntry(g,vp));

    // Map to local subsystem CacheEntry
    rep->subsystems[subsys].cacheEntries.push_back(nxt);
    return (int)rep->subsystems[subsys].cacheEntries.size() - 1;
}

// You can call these as long as stage >= Modeled.
const Real&
State::getTime() const {
    SimTK_STAGECHECK_GE(getStage(), Stage::Modeled, "State::getTime()");
    return rep->t;
}

const Vector&
State::getQ() const {
    SimTK_STAGECHECK_GE(getStage(), Stage::Modeled, "State::getQ()");
    return rep->q;
}

const Vector&
State::getU() const {
    SimTK_STAGECHECK_GE(getStage(), Stage::Modeled, "State::getU()");
    return rep->u;
}

const Vector&
State::getZ() const {
    SimTK_STAGECHECK_GE(getStage(), Stage::Modeled, "State::getZ()");
    return rep->z;
}


// You can call these as long as stage >= Modeled, but the
// stage will be backed up if necessary to the indicated stage.
Real&
State::updTime() {  // Stage::Timed-1
    SimTK_STAGECHECK_GE(getStage(), Stage::Modeled, "State::updTime()");
    invalidateStage(Stage::Timed);
    return rep->t;
}

Vector&
State::updQ() {    // Stage::Configured-1
    SimTK_STAGECHECK_GE(getStage(), Stage::Modeled, "State::updQ()");
    invalidateStage(Stage::Configured);
    return rep->q;
}

Vector&
State::updU() {     // Stage::Moving-1
    SimTK_STAGECHECK_GE(getStage(), Stage::Modeled, "State::updU()");
    invalidateStage(Stage::Moving);
    return rep->u;
}

Vector&
State::updZ() {     // Stage::Dynamics-1
    SimTK_STAGECHECK_GE(getStage(), Stage::Modeled, "State::updZ()");
    invalidateStage(Stage::Dynamics);
    return rep->z;
}


const Vector&
State::getQDot() const {
    SimTK_STAGECHECK_GE(getStage(), Stage::Moving, "State::getQDot()");
    return rep->qdot;
}

const Vector&
State::getZDot() const {
    SimTK_STAGECHECK_GE(getStage(), Stage::Dynamics, "State::getZDot()");
    return rep->zdot;
}

const Vector&
State::getUDot() const {
    SimTK_STAGECHECK_GE(getStage(), Stage::Reacting, "State::getUDot()");
    return rep->udot;
}

const Vector&
State::getQDotDot() const {
    SimTK_STAGECHECK_GE(getStage(), Stage::Reacting, "State::getQDotDot()");
    return rep->qdotdot;
}

Vector&
State::updQDot() const {
    SimTK_STAGECHECK_GE(getStage(), Stage(Stage::Moving).prev(), "State::updQDot()");
    return rep->qdot;
}

Vector&
State::updZDot() const {
    SimTK_STAGECHECK_GE(getStage(), Stage(Stage::Dynamics).prev(), "State::updZDot()");
    return rep->zdot;
}

Vector&
State::updUDot() const {
    SimTK_STAGECHECK_GE(getStage(), Stage(Stage::Reacting).prev(), "State::updUDot()");
    return rep->udot;
}

Vector&
State::updQDotDot() const {
    SimTK_STAGECHECK_GE(getStage(), Stage(Stage::Reacting).prev(), "State::updQDotDot()");
    return rep->qdotdot;
}

// You can access a Modeling variable any time, but don't access others
// until you have realized the Modeled stage.
const AbstractValue& 
State::getDiscreteVariable(int subsys, int index) const {
    const StateRep::PerSubsystemInfo& info = rep->subsystems[subsys];

    SimTK_INDEXCHECK(0,index,(int)info.discreteVars.size(),"State::getDiscreteVariable()");
    const DiscreteVariable& dv = rep->discrete[info.discreteVars[index]];

    if (dv.getStage() > Stage::Modeled) {
        SimTK_STAGECHECK_GE(getStage(subsys), Stage::Modeled, "State::getDiscreteVariable()");
    }

    return dv.getValue();
}

// You can update a Modeling variable from Built stage, but higher variables 
// must wait until you have realized the Modeled stage. This always backs the 
// stage up to one earlier than the variable's stage.
AbstractValue& 
State::updDiscreteVariable(int subsys, int index) {
    const StateRep::PerSubsystemInfo& info = rep->subsystems[subsys];

    SimTK_INDEXCHECK(0,index,(int)info.discreteVars.size(),"State::updDiscreteVariable()");
    DiscreteVariable& dv = rep->discrete[info.discreteVars[index]];

    SimTK_STAGECHECK_GE(getStage(subsys), std::min(dv.getStage().prev(), Stage(Stage::Modeled)), 
        "State::updDiscreteVariable()");

    invalidateStage(subsys, dv.getStage());

    return dv.updValue();
}

// Stage >= ce.stage
const AbstractValue& 
State::getCacheEntry(int subsys, int index) const {
    const StateRep::PerSubsystemInfo& info = rep->subsystems[subsys];

    SimTK_INDEXCHECK(0,index,(int)info.cacheEntries.size(),"State::getCacheEntry()");
    const CacheEntry& ce = rep->cache[info.cacheEntries[index]];

    SimTK_STAGECHECK_GE(getStage(subsys), ce.getStage(), "State::getCacheEntry()");

    return ce.getValue();
}

// Stage >= ce.stage-1; does not change stage
AbstractValue& 
State::updCacheEntry(int subsys, int index) const {
    const StateRep::PerSubsystemInfo& info = rep->subsystems[subsys];

    SimTK_INDEXCHECK(0,index,(int)info.cacheEntries.size(),"State::updCacheEntry()");
    CacheEntry& ce = rep->cache[info.cacheEntries[index]];

    SimTK_STAGECHECK_GE(getStage(subsys), ce.getStage().prev(), "State::updCacheEntry()");

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

    out += "<DISCRETE TODO>\n";

    for (int ss=0; ss < (int)rep->subsystems.size(); ++ss) {
        const StateRep::PerSubsystemInfo& info = rep->subsystems[ss];
        out += "<Subsystem index=" + String(ss) + " name=" + info.name 
            + " version=" + info.version + ">\n";

        out += indexToString(info.qs, "qs");
        out += indexToString(info.us, "us");
        out += indexToString(info.zs, "zs");
        out += indexToString(info.discreteVars, "discreteVars");
        out += indexToString(info.cacheEntries, "cacheEntries");

        out += "</Subsystem>\n";
    }

    out += "</State>\n";
    return out;
}

String State::cacheToString() const {
    String out;
    out += "<Cache>\n";

    for (int ss=0; ss < (int)rep->subsystems.size(); ++ss) {
        const StateRep::PerSubsystemInfo& info = rep->subsystems[ss];
        out += "<Subsystem index=" + String(ss) + " name=" + info.name 
            + " version=" + info.version + ">\n";
        out += "  <Stage>" + info.currentStage.name() + "</Stage>\n";
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

    out += "<DISCRETE CACHE TODO>\n";

    out += "</Cache>\n";
    return out;
}

} // namespace SimTK

