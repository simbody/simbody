#ifndef SimTK_SimTKCOMMON_STATE_H_
#define SimTK_SimTKCOMMON_STATE_H_

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

#include <ostream>
#include <cassert>
#include <set>

namespace SimTK {

SimTK_DEFINE_UNIQUE_ID_TYPE(SubsystemId)

// TODO: these need an option to have associated "update" variables in the cache,
// analogous to the derivative variables qdot,udot,zdot that we create
// for the continuous variables. Consider whether "discrete variable" should
// be reserved for those that are updated in time, with something else like
// "parameter variable" for those that just hold externally set data.

class SimTK_SimTKCOMMON_EXPORT DiscreteVariable {
public:
    DiscreteVariable() : rep(0) { }
    DiscreteVariable(const DiscreteVariable&);
    DiscreteVariable& operator=(const DiscreteVariable&);
    ~DiscreteVariable();

    // This takes ownership of the AbstractValue pointer.
    DiscreteVariable(Stage, AbstractValue* vp);

    Stage getStage() const;
    const AbstractValue& getValue() const;
    AbstractValue&       updValue();

private:
    class DiscreteVariableRep* rep;
};

class SimTK_SimTKCOMMON_EXPORT CacheEntry : public DiscreteVariable {
public:
    CacheEntry() : DiscreteVariable() { }

    // This takes ownership of the AbstractValue pointer.
    CacheEntry(Stage g, AbstractValue* vp)
        : DiscreteVariable(g,vp) { }

    CacheEntry(const CacheEntry& ce) : DiscreteVariable(ce) { }
    CacheEntry& operator=(const CacheEntry& ce) {
        DiscreteVariable::operator=(ce);
        return *this;
    }
};

class EventStatus {
public:
    EventStatus() { initialize(); }
    // default destructor, copy constructor, copy assignment

    // Event trigger (which zero crossings cause triggering). Can be
    // OR'ed together to make a mask.
    enum EventTrigger {
        NoEventTrigger          =0x0000,    // must be 0

        ZeroToPositive          =0x0001,    //  1
        ZeroToNegative          =0x0002,    //  2
        PositiveToZero          =0x0004,    //  4
        NegativeToZero          =0x0008,    //  8
        PositiveToNegative      =0x0010,    // 16
        NegativeToPositive      =0x0020,    // 32

        ZeroToNonzero           =(ZeroToPositive|ZeroToNegative),        // 3
        NonzeroToZero           =(NegativeToZero|PositiveToZero),        // 12
        PositiveToNonpositive   =(PositiveToZero|PositiveToNegative),    // 20 
        Falling                 =(PositiveToNonpositive|ZeroToNegative), // 22
        NegativeToNonnegative   =(NegativeToZero|NegativeToPositive),    // 40
        Rising                  =(NegativeToNonnegative|ZeroToPositive), // 41
        AnySignChange           =(ZeroToNonzero|NegativeToNonnegative    // 63
                                  |PositiveToNonpositive)
    };

    bool isEventPending() const {return transitionSeen != NoEventTrigger;}
    EventTrigger getEventTrigger() const {return transitionSeen;}
    Real getLastTriggerTime() const {return lastTriggerTime;}
    Real getLastTriggerTimeBestGuess() const {return lastTriggerTimeBestGuess;}
    Real getBeforeValue() const {return beforeValue;}
    Real getAfterValue() const {return afterValue;}
    Real getLocalizationWindow() const {return localizationWindow;}

    void setEventTriggered(EventTrigger transition, Real triggerTime,
                           Real actualTimeEst, Real window,
                           Real before, Real after)
    {
        assert(transition != NoEventTrigger);
        assert(triggerTime >= 0 && actualTimeEst >= 0 
               && triggerTime >= actualTimeEst);

        transitionSeen = transition;
        lastTriggerTime = triggerTime;
        lastTriggerTimeBestGuess = actualTimeEst;
        localizationWindow = window;
        beforeValue = before;
        afterValue  = after;
    }

    void clearEventTrigger() {
        transitionSeen = NoEventTrigger;
    }

    // Classify a before/after sign transition.
    static EventTrigger classifyTransition(int before, int after) {
        if (before==after) return NoEventTrigger;
        if (before==0)
            return after==1 ? ZeroToPositive : ZeroToNegative;
        if (before==1)
            return after==0 ? PositiveToZero : PositiveToNegative;
        // before==-1
        return after==0 ? NegativeToZero : NegativeToPositive;
    }

    static EventTrigger maskTransition(EventTrigger transition, EventTrigger mask) {
        return EventTrigger(transition & mask); // we're depending on NoEventTrigger==0
    }

    SimTK_SimTKCOMMON_EXPORT static String eventTriggerString(EventTrigger e);
private:
    void initialize() {
        transitionSeen = NoEventTrigger;
        lastTriggerTime = lastTriggerTimeBestGuess = localizationWindow
            = beforeValue = afterValue = NaN;
    }

    EventTrigger transitionSeen;
    Real         lastTriggerTime; // digital
    Real         lastTriggerTimeBestGuess; // analog, <=lastTriggerTime
    Real         localizationWindow;
    Real         beforeValue, afterValue;
};

/**
 * This is the handle class for the hidden State implementation.
 *
 * This object is intended to contain all State information for a
 * SimTK::System, except topological information which is stored
 * in the System itself. A System is "const" after its topology
 * has been constructed and realized.
 *
 * Systems are composed of Subsystems, and the State supports that
 * concept by allowing per-subsystem partitioning of the total System
 * state. This allows subsytems to have their own private state
 * variables, while permitting the system to allow shared access
 * to state among the subsystems when necessary.
 *
 * The State provides services reflecting the structure of the equations
 * it expects to find in the System. Three different views of the
 * same state information are supported to accommodate three 
 * different users:
 *    - the System as a whole
 *    - Subsystems contained in the System
 *    - numerical methods operating on the state
 * Typically numerical methods have a much less nuanced view of the state
 * than do the System or Subsystems.
 *
 * The system is expected to be a "hybrid DAE", that is, a mixture
 * of continuous and discrete dynamic equations, and algebraic constraints.
 * There is an independent variable t, continuous state variables y, and
 * discrete state variables d.
 *
 * The continuous part is an ODE-on-a-manifold system suitable for solution
 * via coordinate projection, structured like this for the view taken
 * by numerical methods:
 *         (1)  y' = f(d;t,y)         differential equations
 *         (2)  c  = c(d;t,y)         algebraic equations (manifold is c=0)
 *         (3)  e  = e(d;t,y)         event triggers (watch for zero crossings)
 * with initial conditions t0,y0,d0 such that c=0. The discrete variables d
 * are updated upon occurence of specific events, which are
 * detected using the set of scalar-valued event trigger functions e (3).
 *
 * In the more detailed view as seen from the System, we consider y={q,u,z}
 * to be partitioned into position variables q, velocity variables u, and
 * auxiliary variables z. There will be algebraic constraints involving q, u,
 * and u's time derivatives udot. The system is now assumed to look like this:
 *      (4) qdot    = Q(q) u
 *      (5) zdot    = zdot(d;t,q,u,z)
 *
 *      (6)           M(q) udot + ~G(q) mult = f(d;t,q,u,z)
 *      (7) udotErr = G(q) udot - b(d;t,q,u) = 0
 *      (8) uErr    = v(d;t,q,u)             = 0
 *      (9) qErr    = p(d;t,q)               = 0
 *
 * Here G = [A;V] with A(q) being the coefficient matrix for constraints
 * appearing only at the acceleration level, and V(q)=partial(v)/partial(u).
 * (Note that v in Eq 8 is assumed to include equations resulting from
 * differentiation of Eq 9, as well as ones first introduced at the
 * velocity level (nonholonomic constraints).
 *
 * If a system allocates nq q's, nu u's, and nz z's the State will also
 * allocate matching cache variables qdot, qdotdot, udot, and zdot. If
 * np position constraints (9), nv velocity constraints (8) and 
 * na acceleration constraints (7) are allocated, the state creates
 * cache entries of like sizes qErr, uErr, udotErr. In addition room
 * for the na Lagrange multipliers 'mult' is allocated in the cache.
 *
 * In the final view, the Subsystem view, the same variables and cache
 * entries exist, but only the ones allocated by that Subsystem are
 * visible. All of a Subsystem's q's are conscutive in memory, as are
 * its u's, uErr's, etc., but the q's are not adjacent to the u's as
 * they are for the System's view.
 *
 * The default constructor creates a State containing no state variables
 * and with its realization cache stage set to Stage::Empty.
 * During Subsystem construction, variables and cache entries for any
 * stage can be allocated, however *all* Model stage variables
 * must be allocated during this time. At the end of construction,
 * call advanceSubsystemToStage(Topology) which will put the Subsystem
 * at Stage::Topology. Then the Subsystems realize their Model stages, during which 
 * variables at any stage > Model, and cache entries at any stage
 * >= Model can be allocated. After that call advanceSubsystemToStage(Model)
 * which sets the stage to Stage::Model and disallows further allocation.
 *
 * Note that there is a global Stage for the state as a whole, and individual
 * Stages for each subsystem. The global stage can never be higher than
 * the lowest subsystem stage. Global resources are allocated when the
 * global Stage advances to "Model" and tossed out if that stage is
 * invalidated. Note that subsystems will "register" their use of the
 * global variable pools during their own modeling stages, but that the
 * actual global resources won't exist until the *system* has been
 * advanced to Model stage.
 */
class SimTK_SimTKCOMMON_EXPORT State {
public:
    /// Create an empty State.
    State();
    ~State();

    /// Restore State to default-constructed condition.
    void clear();

    /// Set the number of subsystems in this state. This is done during
    /// initialization of the State by a System; it completely wipes out
    /// anything that used to be in the state so use cautiously!
    void setNSubsystems(int i);

    /// Set the name and version for a given subsystem, which must already
    /// have a slot allocated.
    void initializeSubsystem(SubsystemId, const String& name, const String& version);

    /// Make the current State a copy of the source state, copying only
    /// state variables and not the cache. If the source state hasn't
    /// been realized to Model stage, then we don't copy its state
    /// variables either, except those associated with the Topology stage.
    State(const State&);

    /// Make the current State a copy of the source state, copying only
    /// state variables and not the cache. If the source state hasn't
    /// been realized to Model stage, then we don't copy its state
    /// variables either, except those associated with the Topology stage.
    State& operator=(const State&);

    /// Register a new subsystem as a client of this State. The
    /// supplied strings are stored with the State but are not
    /// interpreted by it. The intent is that they can be used to
    /// perform "sanity checks" on deserialized States to make
    /// sure they match the currently instantiated System.
    /// The subsystem index (a small integer) is returned.
    int addSubsystem(const String& name, const String& version);

    int getNSubsystems() const;
    const String& getSubsystemName   (SubsystemId) const;
    const String& getSubsystemVersion(SubsystemId) const;
    const Stage&  getSubsystemStage  (SubsystemId) const;

    /// This returns the *global* stage for this State.
    const Stage& getSystemStage() const;

    /// If any subsystem or the system stage is currently at or
    /// higher than the passed-in one, back up to the stage just prior.
    /// Otherwise do nothing.
    void invalidateAll(Stage) const;  // cache is mutable

    /// Advance the current stage by one to the indicated stage.
    /// The stage is passed in just to give us a chance to verify
    /// that all is as expected. You can only advance one stage at
    /// a time. Advancing to "Topology" and "Model" stages affect
    /// what you can do later.
    void advanceSubsystemToStage(SubsystemId, Stage) const;
    void advanceSystemToStage(Stage) const;

    /// These are shared among all the subsystems and are not allocated until
    /// the *System* is advanced to Stage::Model. The returned index is
    /// local to each subsystem. After the System is modeled, we guarantee that
    /// all the q's for a subsystem will be contiguous, and similarly for u's
    /// and z's. However, q,u,z will *not* be contiguous with each other.
    /// The *global* y is contiguous, and global q,u,z are contiguous within
    /// y, in that order.

    int allocateQ(SubsystemId, const Vector& qInit); // qdot, qdotdot also allocated in cache
    int allocateU(SubsystemId, const Vector& uInit); // udot                    "
    int allocateZ(SubsystemId, const Vector& zInit); // zdot                    "

    /// Slots for constraint errors are handled similarly, although these are
    /// just cache entries not state variables. Q errors and U errors
    /// will each be contiguous for a given subsystem, but *not* with each other.
    /// However, yerr={qerr,uerr} *is* a single contiguous vector.
    /// UDotErr is a separate quantity, not part of yerr. Again the UDotErr's for
    /// each subsystem will be contiguous within the larger UDotErr Vector.
    /// Allocating a UDotErr has the side effect of allocating another Vector
    /// of the same size in the cache for the corresponding Lagrange multipliers,
    /// and these are partitioned identically to UDotErrs.

    int allocateQErr   (SubsystemId, int nqerr);    // these are cache entries
    int allocateUErr   (SubsystemId, int nuerr);
    int allocateUDotErr(SubsystemId, int nudoterr);

    /// Slots for event witness values are similar to constraint errors.
    /// However, this also allocates a discrete state variable to hold
    /// the "triggered" indication. The Stage here is the stage at which
    /// the event witness function can first be examined.
    int allocateEvent(SubsystemId, Stage, int nevent);

    /// These are private to each subsystem and are allocated immediately.
    /// TODO: do discrete variables need an "update" variable in the cache?
    int allocateDiscreteVariable(SubsystemId, Stage, AbstractValue* v);
    int allocateCacheEntry      (SubsystemId, Stage, AbstractValue* v);
    
    /// Dimensions. These are valid at Stage::Model while access to the various
    /// arrays may have stricter requirements. Hence it is better to use these
    /// routines than to get a reference to a Vector and ask for its size().

    int getNY() const; // = nq+nu+nz
    int getQStart() const; int getNQ() const;
    int getUStart() const; int getNU() const;
    int getZStart() const; int getNZ() const;

    int getNYErr() const; // = nqerr+nuerr
    int getQErrStart() const; int getNQErr() const;
    int getUErrStart() const; int getNUErr() const;

    int getNUDotErr() const;
    int getNMultipliers() const;

    int getQStart(SubsystemId)       const; int getNQ(SubsystemId)       const;
    int getUStart(SubsystemId)       const; int getNU(SubsystemId)       const;
    int getZStart(SubsystemId)       const; int getNZ(SubsystemId)       const;

    int getQErrStart(SubsystemId)    const; int getNQErr(SubsystemId)    const;
    int getUErrStart(SubsystemId)    const; int getNUErr(SubsystemId)    const;
    int getUDotErrStart(SubsystemId) const; int getNUDotErr(SubsystemId) const;
    int getMultipliersStart(SubsystemId i) const;
    int getNMultipliers(SubsystemId i)     const;

        // Event handling
    int getNEvents() const; // total
    int getEventStartByStage(Stage) const; // per-stage
    int getNEventsByStage(Stage) const;
    int getEventStartByStage(SubsystemId, Stage) const;
    int getNEventsByStage(SubsystemId, Stage) const;

    const Vector& getEvents() const;
    const Vector& getEventsByStage(Stage) const;
    const Vector& getEventsByStage(SubsystemId, Stage) const;

    Vector& updEvents() const; // mutable
    Vector& updEventsByStage(Stage) const;
    Vector& updEventsByStage(SubsystemId, Stage) const;

    /// Per-subsystem access to the global shared variables.
    const Vector& getQ(SubsystemId) const;
    const Vector& getU(SubsystemId) const;
    const Vector& getZ(SubsystemId) const;

    Vector& updQ(SubsystemId);
    Vector& updU(SubsystemId);
    Vector& updZ(SubsystemId);

    /// Per-subsystem access to the shared cache entries.
    const Vector& getQDot(SubsystemId) const;
    const Vector& getUDot(SubsystemId) const;
    const Vector& getZDot(SubsystemId) const;
    const Vector& getQDotDot(SubsystemId) const;

    Vector& updQDot(SubsystemId) const;    // these are mutable
    Vector& updUDot(SubsystemId) const;
    Vector& updZDot(SubsystemId) const;
    Vector& updQDotDot(SubsystemId) const;

    const Vector& getQErr(SubsystemId) const;
    const Vector& getUErr(SubsystemId) const;
    const Vector& getUDotErr(SubsystemId) const;
    const Vector& getMultipliers(SubsystemId) const;
    Vector& updQErr(SubsystemId) const;    // these are mutable
    Vector& updUErr(SubsystemId) const;
    Vector& updUDotErr(SubsystemId) const;
    Vector& updMultipliers(SubsystemId) const;

    /// You can call these as long as *system* stage >= Model.
    const Real&   getTime() const;
    const Vector& getY() const; // {Q,U,Z} packed and in that order

    /// These are just views into Y.
    const Vector& getQ() const;
    const Vector& getU() const;
    const Vector& getZ() const;

    /// You can call these as long as stage >= Model, but the
    /// stage will be backed up if necessary to the indicated stage.
    Real&   updTime();  // Back up to Stage::Time-1
    Vector& updY();     // Back up to Stage::Congfigured-1

    /// An alternate syntax equivalent to updTime() and updY().
    void setTime(Real t);
    void setY(const Vector& y);

    /// These are just views into Y.
    Vector& updQ();     // Back up to Stage::Position-1
    Vector& updU();     // Back up to Stage::Velocity-1
    Vector& updZ();     // Back up to Stage::Dynamics-1

    /// Alternate interface.
    void setQ(const Vector& q);
    void setU(const Vector& u);
    void setZ(const Vector& z);

    const Vector& getYDot()    const; // Stage::Acceleration

    /// These are just views into YDot.
    const Vector& getQDot()    const; // Stage::Velocity
    const Vector& getZDot()    const; // Stage::Dynamics
    const Vector& getUDot()    const; // Stage::Acceleration

    /// This has its own space, not a view.
    const Vector& getQDotDot() const; // Stage::Acceleration

    /// These are mutable
    Vector& updYDot() const;    // Stage::Acceleration-1
    Vector& updQDot() const;    // Stage::Velocity-1     (view into YDot)
    Vector& updZDot() const;    // Stage::Dynamics-1            "
    Vector& updUDot() const;    // Stage::Acceleration-1        "

    /// This is a separate shared cache entry, not part of YDot. If you
    /// have a direct 2nd order integrator you can integrate QDotDot
    /// (twice) to get Q.
    Vector& updQDotDot() const; // Stage::Acceleration-1

    /// Return the current constraint errors for all constraints.
    const Vector& getYErr() const; // {QErr,UErr} packed and in that order

    /// These are just views into YErr.
    const Vector& getQErr() const;  // Stage::Position (index 3 constraints)
    const Vector& getUErr() const;  // Stage::Velocity (index 2 constraints)

    /// These have their own space, the are not views.
    const Vector& getUDotErr()     const; // Stage::Acceleration (index 1 constraints)
    const Vector& getMultipliers() const; // Stage::Acceleration

    /// These are mutable
    Vector& updYErr() const; // Stage::Dynamics-1
    Vector& updQErr() const; // Stage::Position-1 (view into YErr)
    Vector& updUErr() const; // Stage::Velocity-1        "

    Vector& updUDotErr()     const; // Stage::Acceleration-1 (not a view)
    Vector& updMultipliers() const; // Stage::Acceleration-1 (not a view)

    /// OK if dv.stage==Model or stage >= Model
    const AbstractValue& getDiscreteVariable(SubsystemId, int index) const;

    /// OK if dv.stage==Model or stage >= Model; set stage to dv.stage-1
    AbstractValue&       updDiscreteVariable(SubsystemId, int index);

    /// Alternate interface to updDiscreteVariable.
    void setDiscreteVariable(SubsystemId i, int index, const AbstractValue& v);

    /// Stage >= ce.stage
    const AbstractValue& getCacheEntry(SubsystemId, int index) const;

    /// Stage >= ce.stage-1; does not change stage
    AbstractValue& updCacheEntry(SubsystemId, int index) const; // mutable
    
    /// Create a new State object which shares all the same data as this one,
    /// such that modifying either one will modify both of them.  The new object
    /// restricts which stages and subsystems may be modified.  Any
    /// attempt to modify restricted data through that object will produce an
    /// exception.
    ///
    /// This method can only add restrictions, not remove them.  If this State was
    /// itself created by createRestrictedState(), the new state will inherit
    /// all of the restrictions from this one, in addition to any that are specified
    /// in the arguments.
    State createRestrictedState(EnumerationSet<Stage> restrictedStages, std::set<SubsystemId> restrictedSubsystems);

    /// Get the set of stages which cannot be modified in this State.  Attempting
    /// to modify any of these stages will produce an exception.
    EnumerationSet<Stage> getRestrictedStages() const;

    /// Get the set of subsystems which cannot be modified in this State.  Attempting
    /// to modify any of these subsystems will produce an exception.
    std::set<SubsystemId> getRestrictedSubsystems() const;

    String toString() const;
    String cacheToString() const;

private:
    class StateRep* rep;
    const StateRep& getRep() const {assert(rep); return *rep;}
    StateRep&       updRep()       {assert(rep); return *rep;}
    State(State& state, EnumerationSet<Stage>& restrictedStages, std::set<SubsystemId>& restrictedSubsystems);
};

SimTK_SimTKCOMMON_EXPORT std::ostream& 
operator<<(std::ostream& o, const State& s);

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_STATE_H_
