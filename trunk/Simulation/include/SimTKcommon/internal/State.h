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
 * Portions copyright (c) 2005-9 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
 * Contributors: Peter Eastman                                                *
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
#include "SimTKcommon/internal/Event.h"

#include <ostream>
#include <cassert>
#include <set>

namespace SimTK {


/// Provide a unique integer type for identifying Subsystems.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(SubsystemIndex);

/// This unique integer type is for indexing the global, System-level "y-like"
/// arrays, that is, the arrays in which all of the various Subsystems' continuous
/// state variables q, u, and z have been collected into contiguous memory.
/// This type should be used for the global y and its global time derivative yDot.
/// Note that there is no Subsystem-local equivalent of the y array.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(SystemYIndex);

/// This unique integer type is for indexing global "q-like" arrays, that is, arrays
/// that inherently have the same dimension as the total number of second
/// order state variables (generalized coordinates) in the full System-level
/// view of the State. This type should be used for the global q and its global
/// time derivatives qDot and qDotDot.
/// @see QIndex for Subsystem-local q indexing
SimTK_DEFINE_UNIQUE_INDEX_TYPE(SystemQIndex);
/// Unique integer type for Subsystem-local q indexing
/// @see SystemQIndex
SimTK_DEFINE_UNIQUE_INDEX_TYPE(QIndex);

/// This unique integer type is for indexing global "u-like" arrays, that is, arrays
/// that inherently have the same dimension as the total number of mobilities
/// (generalized speeds) in the full System-level view of the State. This type
/// should be used for the global u and its global time derivative uDot.
/// @see UIndex for Subsystem-local u indexing
SimTK_DEFINE_UNIQUE_INDEX_TYPE(SystemUIndex);
/// Unique integer type for Subsystem-local u indexing
/// @see SystemUIndex
SimTK_DEFINE_UNIQUE_INDEX_TYPE(UIndex);

/// This unique integer type is for indexing global "z-like" arrays, that is, arrays
/// that inherently have the same dimension as the total number of auxiliary
/// state variables in the full System-level view of the State. This type
/// should be used for the global z and its global time derivative zDot.
/// @see ZIndex for Subsystem-local z indexing
SimTK_DEFINE_UNIQUE_INDEX_TYPE(SystemZIndex);
/// Unique integer type for Subsystem-local u indexing
/// @see SystemZIndex
SimTK_DEFINE_UNIQUE_INDEX_TYPE(ZIndex);

/// This unique integer type is for selecting discrete variables. These indices
/// are always Subsystem-local, that is, the first discrete variable belonging
/// to each Subsystem has index 0.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(DiscreteVariableIndex);

/// This unique integer type is for selecting non-shared cache entries. These indices
/// are always Subsystem-local, that is, the first explicitly-allocated cache
/// entry belonging to each Subsystem has index 0.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(CacheEntryIndex);

/// This unique integer type is for indexing the global, System-level "yErr-like"
/// arrays, that is, the arrays in which all of the various Subsystems' qErr and
/// uErr constraint equation slots have been collected together.
/// Note that there is no Subsystem-local equivalent of the yErr array.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(SystemYErrIndex);

/// This unique integer type is for indexing global "qErr-like" arrays, that is, arrays
/// that inherently have the same dimension as the total number of position-level
/// constraint equations in the full System-level view of the State.
/// @see QErrIndex for Subsystem-local qErr indexing
SimTK_DEFINE_UNIQUE_INDEX_TYPE(SystemQErrIndex);
/// Unique integer type for Subsystem-local qErr indexing
/// @see SystemQErrIndex
SimTK_DEFINE_UNIQUE_INDEX_TYPE(QErrIndex);

/// This unique integer type is for indexing global "uErr-like" arrays, that is, arrays
/// that inherently have the same dimension as the total number of velocity-level
/// constraint equations in the full System-level view of the State.
/// @see UErrIndex for Subsystem-local uErr indexing
SimTK_DEFINE_UNIQUE_INDEX_TYPE(SystemUErrIndex);
/// Unique integer type for Subsystem-local uErr indexing
/// @see SystemUErrIndex
SimTK_DEFINE_UNIQUE_INDEX_TYPE(UErrIndex);

/// This unique integer type is for indexing global "uDotErr-like" arrays, that is, arrays
/// that inherently have the same dimension as the total number of acceleration-level
/// constraint equations in the full System-level view of the State.
/// @see UDotErrIndex for Subsystem-local uDotErr indexing
SimTK_DEFINE_UNIQUE_INDEX_TYPE(SystemUDotErrIndex);
/// Unique integer type for Subsystem-local uDotErr indexing
/// @see SystemUDotErrIndex
SimTK_DEFINE_UNIQUE_INDEX_TYPE(UDotErrIndex);

/// This unique integer type is for indexing global "multiplier-like" arrays, that is, arrays
/// that inherently have the same dimension as the total number of Lagrange multipliers
/// in the full System-level view of the State.
/// @see MultiplierIndex for Subsystem-local multiplier indexing
SimTK_DEFINE_UNIQUE_INDEX_TYPE(SystemMultiplierIndex);
/// Unique integer type for Subsystem-local multiplier indexing
/// @see SystemMultiplierIndex
SimTK_DEFINE_UNIQUE_INDEX_TYPE(MultiplierIndex);

/// This is the type to use for Stage version numbers. Whenever any state variable is
/// modified, we increment the stage version for the stage(s) that depend on it.
/// -1 means "unintialized". 0 is never used as a stage version, but is allowed as
/// a cache value which is guaranteed never to look valid. 
typedef int StageVersion;


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
 *
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
 * <pre>
 *      (1)  y' = f(d;t,y)         differential equations
 *      (2)  c  = c(d;t,y)         algebraic equations (manifold is c=0)
 *      (3)  e  = e(d;t,y)         event triggers (watch for zero crossings)
 * </pre>
 * with initial conditions t0,y0,d0 such that c=0. The discrete variables d
 * are updated upon occurence of specific events, which are
 * detected using the set of scalar-valued event trigger functions e (3).
 *
 * In the more detailed view as seen from the System, we consider y={q,u,z}
 * to be partitioned into position variables q, velocity variables u, and
 * auxiliary variables z. There will be algebraic constraints involving q, u,
 * and u's time derivatives udot. The system is now assumed to look like this:
 * <pre>
 *      (4) qdot    = N(q) u
 *      (5) zdot    = zdot(d;t,q,u,z)
 *
 *      (6) M(q) udot + ~G(q) mult = f(d;t,q,u,z)
 *          G(q) udot              = b(d;t,q,u)
 *
 *      (7) udotErr = [ pdotdot(d;t,q,u,udot) ] = 0
 *                    [ vdot(d;t,q,u,udot)    ]
 *                    [ a(d;t,q,u,udot)       ] 
 *
 *      (8) uErr    = [ pdot(d;t,q,u) ]
 *                    [ v(d;t,q,u)    ]         = 0
 *
 *      (9) qErr    = [ p(d;t,q) ]              = 0
 *                    [ n(q)     ]
 * </pre>
 * The q's can also be dealt with directly as second order variables via
 * <pre>
 *     (10) qdotdot = Ndot(q,qdot) u + N(q) udot
 * </pre>
 *
 * Here G = [P;V;A] with A(q) being the coefficient matrix for constraints
 * appearing only at the acceleration level, and V(q)=partial(v)/partial(u)
 * the coefficient matrix for the velocity (nonholonomic) constraints, and
 * P(q)=partial(pdot)/partial(u) is the coefficient matrix of the first
 * time derivatives of the position (holonomic) constraints.
 * Note that uErr in Eq 8 is assumed to include equations resulting from
 * differentiation of p() in Eq 9, as well as ones first introduced at the
 * velocity level (nonholonomic constraints), and udotErr is similarly 
 * built from acceleration-only constraints a() and derivatives of higher-level
 * constraints.
 *
 * If a system allocates nq q's, nu u's, and nz z's the State will also
 * allocate matching cache variables qdot, qdotdot, udot, and zdot. If
 * mp position (holonomic) constraints (9), mpv velocity constraints (8) and 
 * mpva acceleration constraints (7) are allocated, the state creates
 * cache entries of like sizes qErr, uErr, udotErr. In addition room
 * for the mpva Lagrange multipliers 'mult' is allocated in the cache.
 *
 * In the final view, the Subsystem view, the same variables and cache
 * entries exist, but only the ones allocated by that Subsystem are
 * visible. All of a Subsystem's q's are consecutive in memory, as are
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
 * the lowest subsystem stage. Global state resources are allocated when the
 * global Stage advances to "Model" and tossed out if that stage is
 * invalidated. Similarly cache resources are allocated at stage Instance
 * and forgotten when Instance is invalidated. Note that subsystems will
 * "register" their use of the global variable pools during their own modeling
 * stages, but that the actual global resources won't exist until the *System* 
 * has been advanced to Model or Instance stage.
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
    /// anything that used to be in the State so use cautiously!
    void setNSubsystems(int i);

    /// Set the name and version for a given subsystem, which must already
    /// have a slot allocated.
    void initializeSubsystem(SubsystemIndex, const String& name, const String& version);

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
    SubsystemIndex addSubsystem(const String& name, const String& version);

    int getNSubsystems() const;
    const String& getSubsystemName   (SubsystemIndex) const;
    const String& getSubsystemVersion(SubsystemIndex) const;
    const Stage&  getSubsystemStage  (SubsystemIndex) const;

    /// This returns the *global* stage for this State.
    const Stage& getSystemStage() const;

    /// If any subsystem or the system stage is currently at or
    /// higher than the passed-in one, back up to the stage just prior.
    /// Otherwise do nothing.
    void invalidateAll(Stage) const;  // cache is mutable

    /// Advance a particular Subsystem's current stage by one to
    /// the indicated stage. The stage is passed in just to give us a
    /// chance to verify that all is as expected. You can only advance one stage at
    /// a time. Advancing to Topology, Model, or Instance stage affects
    /// what you can do later.
    /// @see advanceSystemToStage()
    void advanceSubsystemToStage(SubsystemIndex, Stage) const;
    /// Advance the System-level current stage by one to the indicated stage.
    /// This can only be done if <em>all</em> Subsystem have already been
    /// advanced to this Stage.
    /// @see advanceSubsystemToStage()
    void advanceSystemToStage(Stage) const;

    /// These continuous state variables are shared among all the Subsystems
    /// and are not allocated until the *System* is advanced to Stage::Model.
    /// The returned index is local to each Subsystem. After the System is modeled,
    /// we guarantee that all the q's for a Subsystem will be contiguous, and
    /// similarly for u's and z's. However, q,u,z will *not* be contiguous with
    /// each other. The *global* y={q,u,z} is contiguous, and global q,u,z are
    /// contiguous within y, in that order. Corresponding cache entries for
    /// the derivatives of these variables are allocated at Model stage also.

    QIndex allocateQ(SubsystemIndex, const Vector& qInit); // qdot, qdotdot also allocated in cache
    UIndex allocateU(SubsystemIndex, const Vector& uInit); // udot                    "
    ZIndex allocateZ(SubsystemIndex, const Vector& zInit); // zdot                    "

    /// These constraint error cache entries are shared among all the subsystems
    /// and are not allocated until the *System* is advanced to Stage::Instance.
    /// The returned index is local to each Subsystem. Q errors and U errors will
    /// each be contiguous for a given subsystem, but *not* with each other. However,
    /// the System-level yerr={qerr,uerr} *is* a single contiguous vector.
    /// UDotErr is a separate quantity, not part of yerr. Again the UDotErr's for
    /// each subsystem will be contiguous within the larger UDotErr Vector.
    /// Allocating a UDotErr has the side effect of allocating another Vector
    /// of the same size in the cache for the corresponding Lagrange multipliers,
    /// and these are partitioned identically to UDotErrs.

    QErrIndex    allocateQErr   (SubsystemIndex, int nqerr) const;    // these are cache entries
    UErrIndex    allocateUErr   (SubsystemIndex, int nuerr) const;
    UDotErrIndex allocateUDotErr(SubsystemIndex, int nudoterr) const;

    /// This method is used to obtain one or more consecutive EventIndex's which are
    /// unique within the indicated Subsystem. These can be allocated in a State that
    /// has not yet been advanced to Instance stage, and the stage at which they
    /// are allocated is remembered so that reducing the stage below that causes
    /// the more-recently-allocated EventIndex's to be forgotten. When Instance stage is
    /// realized, global SystemEventIndex's will be allocated such that each
    /// (Subsystem,EventIndex) pair has a unique SystemEventIndex. Also, SystemEventIndex's
    /// are doled out consecutively within each Subsystem so finding the first one allows
    /// you to compute the rest for a given Subsystem.
    EventIndex allocateEventIndex(SubsystemIndex, int nevent=1) const;

    /// Some Events require a slot in the State cache to hold the current value
    /// of the event trigger function (a.k.a. event "witness" function).
    /// The Stage here is the stage at which the trigger function's value
    /// should be examined by a TimeStepper. The returned index is local to
    /// the Subsystem and also to the Stage. These can be allocated in a State
    /// that has not yet been realized to Instance stage, and will be forgotten
    /// appropriately if the State is later backed up to an earlier Stage.
    /// When this State is realized to Instance stage, global event trigger slots
    /// will be allocated, collecting all same-Stage event triggers together
    /// consecutively for the convenience of the TimeStepper. Within a stage,
    /// a given Subsystem's event trigger slots for that stage are consecutive.
    EventTriggerByStageIndex allocateEventTrigger(SubsystemIndex, Stage, EventIndex, int nevent=1) const;

    //OBSOLETE
    EventTriggerByStageIndex allocateEventTrigger(SubsystemIndex, Stage, int nevent) const;

    /// You can allocate a new DiscreteVariable in any State whose stage has not yet been
    /// advanced to Model stage. The stage at allocation (Empty or Topology) is remembered
    /// so that the appropriate discrete variables can be forgotten if the State's stage
    /// is reduced back to that stage later after advancing past it. DiscreteVariables are
    /// private to each Subsystem and allocated immediately. The returned
    /// index is unique within the Subsystem but there is no corresponding global index.
    ///
    /// The Stage supplied here in the call is the lowest subsystem stage which is invalidated
    /// by a change made to this discrete variable. You may access the value of the discrete
    /// variable for reading (via getDiscreteVariable()) or writing (via updDiscreteVariable())
    /// any time after it has been allocated. Access for writing has the side effect of
    /// reducing the subsystem and system stages for this State to one stage below the one
    /// supplied here, that is, the stage supplied here is invalidated. Note that you must
    /// have write access to the State in order to change the value of any state variable.
    ///
    /// Ownership of the AbstractValue object supplied here is taken over by the State --
    /// don't delete the object after this call!
    /// @see getDiscreteVariable()
    /// @see updDiscreteVariable()
    DiscreteVariableIndex allocateDiscreteVariable(SubsystemIndex, Stage invalidates, AbstractValue*);


    /// This method allocates a DiscreteVariable and a CacheEntry of the same value type.
    /// This can be done if the State hasn't yet been advanced to Model stage.
    /// The discrete variable is allocated as described for allocateDiscreteVariable(),
    /// except that the \a invalidates stage must be higher than Stage::Time.
    /// The cache entry is allocated as described for allocateCacheEntry() without an
    /// automatic calculation (\a latest) stage. The cache entry is then considered to be the
    /// "update" value for the discrete variable. Update values play a similar role for
    /// discrete variables as derivatives play for continuous variables. That is, they
    /// define how the variable is to be updated when a TimeStepper accepts a step.
    ///
    /// Whenever the value of the discrete variable is accessed, the value of the 
    /// update cache entry is returned instead if it is valid. In that way the results
    /// are always calculated using the value as it will be \e after an update. That
    /// means that no results will change when the swap occurs, so no stage needs
    /// to be invalidated upon updating. However, any \e explicit change to the discrete variable
    /// will invalidate the \a invalidates stage just as for a non-updating discrete
    /// variable. Note that this is entirely analogous to the treatment of continuous
    /// variables like q: the integrator ensures that only updated values of q are
    /// seen when evaluations are made at intermediate or trial steps.
    ///
    /// Update occurs as follows: at the start of every continuous interval, after all
    /// other pending events have been handled, a time stepper should call the State
    /// method updateDiscreteVariables(). That method looks at all the
    /// updating discrete variables to see which ones have valid update values. For
    /// each valid value, the discrete variable and its update value are swapped, and
    /// the new cache value is marked invalid. No stage is invalidated by the swap; see
    /// above for why that isn't necessary.
    ///
    /// Ownership of the AbstractValue object supplied here is taken over by the State --
    /// don't delete the object after this call!
    /// @see allocateDiscreteVariable()
    /// @see allocateCacheEntry()
    std::pair<DiscreteVariableIndex, CacheEntryIndex>
        allocateUpdatingDiscreteVariable(SubsystemIndex, Stage invalidates, AbstractValue*,
                                         Stage earliestUpdate); 

    /// You can allocate a new CacheEntry in any State whose stage has not yet been
    /// advanced to Instance stage. The stage at allocation (Empty, Topology, or Model)
    /// is remembered so that the appropriate cache entries can be forgotten if the
    /// State's stage is reduced back to that stage later after advancing past it. CacheEntries
    /// are private to each Subsystem and allocated immediately. The returned index is
    /// unique within the Subsystem and there is no corresponding global index.
    ///
    /// There are two Stages supplied explicitly as arguments to this method: \a earliest
    /// and \a latest. The \a earliest Stage is the stage at which the cache entry
    /// \e could be calculated. Hence if the Subsystem stage is reduced below \a earliest
    /// the cache entry is known to be invalid. The \a latest Stage, if any, is the stage at
    /// which the cache entry is \e guaranteed to have been calculated. For stages
    /// \a earliest through \a latest-1, the cache entry \e may be valid, if it has 
    /// already been calculated. In that case an explicit validity indicator will have
    /// been set at the time it was computed. That indicator is cleared automatically
    /// whenever the Subsystem stage is reduced below \a earliest. The validity indicator
    /// need not have been set in order for the cache entry to be deemed valid at
    /// \a latest stage.
    ///
    /// If \a latest is given as Stage::Empty then there is no guarantee that this Subsystem
    /// will automatically calculate a value for this cache entry. In that case the only
    /// way the cache entry can become valid is if the calculation is performed and the
    /// validity indicator is set explicitly.
    ///
    /// Prior to the Subsystem advancing to \a earliest stage, and prior to \a latest 
    /// stage unless the validity indicator is set, attempts to look at the value via
    /// getCacheEntry() will throw an exception. However, you may access the cache entry
    /// for writing via updCacheEntry() any time after stage \a earliest-1. If you 
    /// evaluate it prior to \a latest, be sure to explicitly mark it valid.
    /// Note that cache entries are mutable so you do not need write
    /// access to the State in order to access a cache entry for writing.
    ///
    /// Ownership of the AbstractValue object supplied here is taken over by the State --
    /// don't delete the object after this call! 
    /// @see getCacheEntry()
    /// @see updCacheEntry()
    CacheEntryIndex allocateCacheEntry(SubsystemIndex, Stage earliest, Stage latest,
                                       AbstractValue*) const;

    /// This is an abbreviation for allocation of a cache entry whose earliest
    /// and latest Stages are the same. That is, this cache entry is guaranteed to
    /// be valid if its Subsystem has advanced to the supplied Stage or later, and
    /// is guaranteed to be invalid below that Stage.
    CacheEntryIndex allocateCacheEntry(SubsystemIndex sx, Stage g, AbstractValue* v) const;
    /* should be ->  {   return allocateCacheEntry(sx,g,g,v); } */

    
    /// @name Global Resource Dimensions
    ///
    /// These are the dimensions of the global shared state and cache resources,
    /// as well as the dimensions of the per-Subsystem partioning of those
    /// resources. State resource dimensions (including cache resources directly
    /// related to state variables) are known after the System has been
    /// realized to Model stage. Other cache resource dimensions are known after
    /// the System has been realized to Instance stage. Access to the actual data arrays
    /// may have stricter requirements (for example, you can't ask to look at UErr
    /// arrays until Velocity stage). Hence it is better to use these explicit 
    /// dimension-providing methods than to get a reference to a Vector and ask
    /// for its size().
    ///
    /// @see Per-Subsystem Dimensions group
    /// @see Global-to-Subsystem Maps group
    /// @{

    /// Get the total number ny=nq+nu+nz of shared continuous state variables.
    /// This is also the number of state derivatives in the cache entry ydot.
    /// Callable at Model stage.
    int getNY() const;
    /// Get total number of shared q's (generalized coordinates; second order
    /// state variables). This is also the number of first and second q time
    /// derivatives in the cache entries qdot and qdotdot.
    /// Callable at Model stage.
    int getNQ() const;
    /// Returns the y index at which the q's begin. Callable at Model stage.
    SystemYIndex getQStart() const;
    /// Get total number of shared u's (generalized speeds; mobilities). 
    /// This is also the number of u time derivatives in the cache entry udot.
    /// Callable at Model stage.
    int getNU() const;
    /// Returns the y index at which the u's begin. Callable at Model stage.
    SystemYIndex getUStart() const;
    /// Get total number of shared z's (auxiliary state variables). 
    /// This is also the number of z time derivatives in the cache entry zdot.
    /// Callable at Model stage.
    int getNZ() const;
    /// Returns the y index at which the z's begin. Callable at Model stage.
    SystemYIndex getZStart() const;
    /// Get the total number nyerr=nqerr+nuerr of shared cache entries for
    /// position-level and velocity-level constraint errors.
    /// Callable at Instance stage.
    int getNYErr() const;
    /// Return the total number nqerr=mp+nQuaternions of cache entries for
    /// position-level constraint errors. Callable at Instance stage.
    int getNQErr() const;
    /// Returns the yErr index at which the qErr's begin. Callable at Instance stage.
    SystemYErrIndex getQErrStart() const; 
    /// Return the total number nuerr=mp+mv of cache entries for
    /// velocity-level constraint errors (including also errors in the 
    /// time derivatives of position-level constraints). Callable at Instance stage.
    int getNUErr() const;
    /// Returns the yErr index at which the uErr's begin. Callable at Instance stage.
    SystemYErrIndex getUErrStart() const; 
    /// Return the total number nudotErr=mp+mv+ma of cache entries for
    /// acceleration-level constraint errors (including also errors in the 
    /// second time derivatives of position-level constraints and the first
    /// time derivatives of velocity-level constraints). Callable at Instance stage.
    int getNUDotErr() const;
    /// Return the total number of constraint multipliers; necessarily the same
    /// as the number of acceleration-level constraint errors nUDotErr. Callable
    /// at Instance stage.
    /// @see getNUDotErr()
    int getNMultipliers() const; // =mp+mv+ma, necessarily the same as NUDotErr
    /// Return the total number of event trigger function slots in the cache.
    /// Callable at Instance stage.
    int getNEventTriggers() const;
    /// Return the size of the partition of event trigger functions which are 
    /// evaluated at a given Stage. Callable at Instance stage.
    int getNEventTriggersByStage(Stage) const;
    /// Return the index within the global event trigger array at which the
    /// first of the event triggers associated with a particular Stage are stored;
    /// the rest follow contiguously. Callable at Instance stage.
    SystemEventTriggerIndex getEventTriggerStartByStage(Stage) const; // per-stage

    /// @}

    /// @name Per-Subsystem Dimensions
    ///
    /// These are the dimensions and locations within the global resource arrays
    /// of state and cache resources allocated to a particular Subsystem. Note
    /// that a Subsystem has contiguous q's, contiguous u's, and contiguous z's
    /// but that the q-, u-, and z-partitions are not contiguous. Hence there is
    /// no Subsystem equivalent of the global y vector.
    ///
    /// These serve as a mapping from Subsystem-local indices for the various
    /// shared resources to their global resource indices.
    ///
    /// @see Global Resource Dimensions 
    /// @{

    SystemQIndex getQStart(SubsystemIndex)       const; int getNQ(SubsystemIndex)       const;
    SystemUIndex getUStart(SubsystemIndex)       const; int getNU(SubsystemIndex)       const;
    SystemZIndex getZStart(SubsystemIndex)       const; int getNZ(SubsystemIndex)       const;

    SystemQErrIndex       getQErrStart(SubsystemIndex)    const; int getNQErr(SubsystemIndex)    const;
    SystemUErrIndex       getUErrStart(SubsystemIndex)    const; int getNUErr(SubsystemIndex)    const;
    SystemUDotErrIndex    getUDotErrStart(SubsystemIndex) const; int getNUDotErr(SubsystemIndex) const;
    SystemMultiplierIndex getMultipliersStart(SubsystemIndex) const;
    int getNMultipliers(SubsystemIndex)     const;

    SystemEventTriggerByStageIndex getEventTriggerStartByStage(SubsystemIndex, Stage) const;
    int getNEventTriggersByStage(SubsystemIndex, Stage) const;

    /// @}

    /// @name Global-to-Subsystem Maps
    ///
    /// Once the dimensions and allocations of the global shared resources
    /// are known, you can call these methods to map a global resource index
    /// to Subsystem to which it belongs and the index by which that resource
    /// is known locally to the Subsystem.
    ///
    /// @see Global Resource Dimensions 
    /// @see Per-Subsystem Dimensions group
    /// @{

    /// For a given global q, return the Subsystem that allocated it and the Subsystem-local
    /// index by which it is known; callable at Model stage.
    void mapQToSubsystem(SystemQIndex, SubsystemIndex&, QIndex&) const;
    /// For a given global u, return the Subsystem that allocated it and the Subsystem-local
    /// index by which it is known; callable at Model stage.
    void mapUToSubsystem(SystemUIndex, SubsystemIndex&, UIndex&) const;
    /// For a given global z, return the Subsystem that allocated it and the Subsystem-local
    /// index by which it is known; callable at Model stage.
    void mapZToSubsystem(SystemZIndex, SubsystemIndex&, ZIndex&) const;
    /// For a given global qErr index, return the Subsystem that allocated it and the Subsystem-local
    /// index by which it is known; callable at Instance stage.
    void mapQErrToSubsystem(SystemQErrIndex, SubsystemIndex&, QErrIndex&) const;
    /// For a given global uErr index, return the Subsystem that allocated it and the Subsystem-local
    /// index by which it is known; callable at Instance stage.
    void mapUErrToSubsystem(SystemUErrIndex, SubsystemIndex&, UErrIndex&) const;
    /// For a given global uDotErr index, return the Subsystem that allocated it and the Subsystem-local
    /// index by which it is known; callable at Instance stage.
    void mapUDotErrToSubsystem(SystemUDotErrIndex, SubsystemIndex&, UDotErrIndex&) const;
    /// For a given global multiplier index, return the Subsystem that allocated it and the Subsystem-local
    /// index by which it is known; callable at Instance stage. This is necessarily the same Subsystme
    /// and index as for the corresponding global uDotErr.
    void mapMultiplierToSubsystem(SystemMultiplierIndex, SubsystemIndex&, MultiplierIndex&) const;
    /// For a given global event trigger function index, return the Subsystem that allocated it and
    /// the Subsystem-local index by which it is known; callable at Instance stage.
    void mapEventTriggerToSubsystem(SystemEventTriggerIndex, SubsystemIndex&, EventTriggerIndex&) const;
    /// For a given global event trigger function index, return the Stage at which that
    /// trigger function should be evaluated; callable at Instance stage.
    void mapEventTriggerToStage(SystemEventTriggerIndex, Stage&, SystemEventTriggerByStageIndex&) const;
    /// Given a Subsystem-wide event index, map that to a particular Stage and an index
    /// within that Stage.
    void mapSubsystemEventTriggerToStage(EventTriggerIndex, Stage&, EventTriggerByStageIndex&) const;

    /// @}

    const Vector& getEventTriggers() const;
    const Vector& getEventTriggersByStage(Stage) const;
    const Vector& getEventTriggersByStage(SubsystemIndex, Stage) const;

    Vector& updEventTriggers() const; // mutable
    Vector& updEventTriggersByStage(Stage) const;
    Vector& updEventTriggersByStage(SubsystemIndex, Stage) const;

    /// Per-subsystem access to the global shared variables.
    const Vector& getQ(SubsystemIndex) const;
    const Vector& getU(SubsystemIndex) const;
    const Vector& getZ(SubsystemIndex) const;

    Vector& updQ(SubsystemIndex);
    Vector& updU(SubsystemIndex);
    Vector& updZ(SubsystemIndex);

    /// Per-subsystem access to the shared cache entries.
    const Vector& getQDot(SubsystemIndex) const;
    const Vector& getUDot(SubsystemIndex) const;
    const Vector& getZDot(SubsystemIndex) const;
    const Vector& getQDotDot(SubsystemIndex) const;

    Vector& updQDot(SubsystemIndex) const;    // these are mutable
    Vector& updUDot(SubsystemIndex) const;
    Vector& updZDot(SubsystemIndex) const;
    Vector& updQDotDot(SubsystemIndex) const;

    const Vector& getQErr(SubsystemIndex) const;
    const Vector& getUErr(SubsystemIndex) const;
    const Vector& getUDotErr(SubsystemIndex) const;
    const Vector& getMultipliers(SubsystemIndex) const;
    Vector& updQErr(SubsystemIndex) const;    // these are mutable
    Vector& updUErr(SubsystemIndex) const;
    Vector& updUDotErr(SubsystemIndex) const;
    Vector& updMultipliers(SubsystemIndex) const;

    /// You can call these as long as *system* stage >= Model.
    const Real&   getTime() const;
    const Vector& getY() const; // {Q,U,Z} packed and in that order

    /// These are just views into Y.
    const Vector& getQ() const;
    const Vector& getU() const;
    const Vector& getZ() const;

    /// You can call these as long as System stage >= Model, but the
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
    const AbstractValue& getDiscreteVariable         (SubsystemIndex, DiscreteVariableIndex) const;
    Real                 getDiscreteVarLastUpdateTime(SubsystemIndex, DiscreteVariableIndex) const;
    CacheEntryIndex      getDiscreteVarUpdateEntry   (SubsystemIndex, DiscreteVariableIndex) const;
    const AbstractValue& getDiscreteVarPrevValue     (SubsystemIndex, DiscreteVariableIndex) const;

    /// OK if dv.stage==Model or stage >= Model; set stage to dv.stage-1
    AbstractValue&       updDiscreteVariable(SubsystemIndex, DiscreteVariableIndex);

    /// Alternate interface to updDiscreteVariable.
    void setDiscreteVariable(SubsystemIndex, DiscreteVariableIndex, const AbstractValue&);

    /// Stage >= ce.stage
    const AbstractValue& getCacheEntry(SubsystemIndex, CacheEntryIndex) const;

    /// Stage >= ce.stage-1; does not change stage
    AbstractValue& updCacheEntry(SubsystemIndex, CacheEntryIndex) const; // mutable

    bool isCacheValueValid(SubsystemIndex, CacheEntryIndex) const;
    void markCacheValueValid(SubsystemIndex, CacheEntryIndex) const;

    /// Return the lowest System Stage that was invalidated since the last time this
    /// "low water mark" was reset. The returned value is never higher than the
    /// one just following the State's current System Stage, but can be lower.
    const Stage& getLowestStageModified() const;

    /// Reset the invalid System Stage "low water mark" to the one just above the
    /// State's current System Stage.
    void resetLowestStageModified() const;
    
    /// Transform a State into one which shares all the same data as this one,
    /// such that modifying either one will modify both of them.  The new State
    /// restricts which stages and subsystems may be modified.  Any
    /// attempt to modify restricted data through that object will produce an
    /// exception.
    ///
    /// This method can only add restrictions, not remove them.  If this State was
    /// itself created by createRestrictedState(), the new state will inherit
    /// all of the restrictions from this one, in addition to any that are specified
    /// in the arguments.
    void createRestrictedState(State&                       restrictedState, 
                               EnumerationSet<Stage>        restrictedStages,
                               std::set<SubsystemIndex>     restrictedSubsystems);

    /// Get the set of stages which cannot be modified in this State.  Attempting
    /// to modify any of these stages will produce an exception.
    const EnumerationSet<Stage>& getRestrictedStages() const;

    /// Get the set of subsystems which cannot be modified in this State.  Attempting
    /// to modify any of these subsystems will produce an exception.
    const std::set<SubsystemIndex>& getRestrictedSubsystems() const;

    String toString() const;
    String cacheToString() const;

private:
    class StateRep* rep;
    const StateRep& getRep() const {assert(rep); return *rep;}
    StateRep&       updRep()       {assert(rep); return *rep;}
};

SimTK_SimTKCOMMON_EXPORT std::ostream& 
operator<<(std::ostream& o, const State& s);

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_STATE_H_
