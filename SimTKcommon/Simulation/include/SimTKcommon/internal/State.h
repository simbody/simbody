#ifndef SimTK_SimTKCOMMON_STATE_H_
#define SimTK_SimTKCOMMON_STATE_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-15 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors: Peter Eastman                                                *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

/** @file
Declares the user-visible part of a SimTK::State, the implementation is
done in a separate internal class. **/

// Note: the StateImpl.h header is included at the end of this file to 
// complete the inline part of the State definition.
#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"

#include <ostream>
#include <cassert>
#include <algorithm>
#include <mutex>
#include <array>

namespace SimTK {


/// @class SimTK::SubsystemIndex
/// Provide a unique integer type for identifying Subsystems.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(SubsystemIndex);

/// @class SimTK::SystemYIndex
/// This unique integer type is for indexing the global, System-level "y-like"
/// arrays, that is, the arrays in which all of the various Subsystems' continuous
/// state variables q, u, and z have been collected into contiguous memory.
/// This type should be used for the global y and its global time derivative yDot.
/// Note that there is no Subsystem-local equivalent of the y array.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(SystemYIndex);

/// @class SimTK::SystemQIndex
/// This unique integer type is for indexing global "q-like" arrays, that is, arrays
/// that inherently have the same dimension as the total number of second
/// order state variables (generalized coordinates) in the full System-level
/// view of the State. This type should be used for the global q and its global
/// time derivatives qDot and qDotDot.
/// @see QIndex for Subsystem-local q indexing
SimTK_DEFINE_UNIQUE_INDEX_TYPE(SystemQIndex);
/// @class SimTK::QIndex
/// Unique integer type for Subsystem-local q indexing
/// @see SystemQIndex
SimTK_DEFINE_UNIQUE_INDEX_TYPE(QIndex);

/// @class SimTK::SystemUIndex
/// This unique integer type is for indexing global "u-like" arrays, that is, arrays
/// that inherently have the same dimension as the total number of mobilities
/// (generalized speeds) in the full System-level view of the State. This type
/// should be used for the global u and its global time derivative uDot.
/// @see UIndex for Subsystem-local u indexing
SimTK_DEFINE_UNIQUE_INDEX_TYPE(SystemUIndex);
/// @class SimTK::UIndex
/// Unique integer type for Subsystem-local u indexing
/// @see SystemUIndex
SimTK_DEFINE_UNIQUE_INDEX_TYPE(UIndex);

/// @class SimTK::SystemZIndex
/// This unique integer type is for indexing global "z-like" arrays, that is, arrays
/// that inherently have the same dimension as the total number of auxiliary
/// state variables in the full System-level view of the State. This type
/// should be used for the global z and its global time derivative zDot.
/// @see ZIndex for Subsystem-local z indexing
SimTK_DEFINE_UNIQUE_INDEX_TYPE(SystemZIndex);
/// @class SimTK::ZIndex
/// Unique integer type for Subsystem-local z indexing
/// @see SystemZIndex
SimTK_DEFINE_UNIQUE_INDEX_TYPE(ZIndex);

/// @class SimTK::DiscreteVariableIndex
/// This unique integer type is for selecting discrete variables. These indices
/// are always Subsystem-local, that is, the first discrete variable belonging
/// to each Subsystem has index 0.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(DiscreteVariableIndex);

/// @class SimTK::CacheEntryIndex
/// This unique integer type is for selecting non-shared cache entries. These indices
/// are always Subsystem-local, that is, the first explicitly-allocated cache
/// entry belonging to each Subsystem has index 0.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(CacheEntryIndex);

/// @class SimTK::SystemYErrIndex
/// This unique integer type is for indexing the global, System-level "yErr-like"
/// arrays, that is, the arrays in which all of the various Subsystems' qErr and
/// uErr constraint equation slots have been collected together.
/// Note that there is no Subsystem-local equivalent of the yErr array.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(SystemYErrIndex);

/// @class SimTK::SystemQErrIndex
/// This unique integer type is for indexing global "qErr-like" arrays, that is, arrays
/// that inherently have the same dimension as the total number of position-level
/// constraint equations in the full System-level view of the State.
/// @see QErrIndex for Subsystem-local qErr indexing
SimTK_DEFINE_UNIQUE_INDEX_TYPE(SystemQErrIndex);
/// @class SimTK::QErrIndex
/// Unique integer type for Subsystem-local qErr indexing
/// @see SystemQErrIndex
SimTK_DEFINE_UNIQUE_INDEX_TYPE(QErrIndex);

/// @class SimTK::SystemUErrIndex
/// This unique integer type is for indexing global "uErr-like" arrays, that is, arrays
/// that inherently have the same dimension as the total number of velocity-level
/// constraint equations in the full System-level view of the State.
/// @see UErrIndex for Subsystem-local uErr indexing
SimTK_DEFINE_UNIQUE_INDEX_TYPE(SystemUErrIndex);
/// @class SimTK::UErrIndex
/// Unique integer type for Subsystem-local uErr indexing
/// @see SystemUErrIndex
SimTK_DEFINE_UNIQUE_INDEX_TYPE(UErrIndex);

/// @class SimTK::SystemUDotErrIndex
/// This unique integer type is for indexing global "uDotErr-like" arrays, that is, arrays
/// that inherently have the same dimension as the total number of acceleration-level
/// constraint equations in the full System-level view of the State.
/// @see UDotErrIndex for Subsystem-local uDotErr indexing
SimTK_DEFINE_UNIQUE_INDEX_TYPE(SystemUDotErrIndex);
/// @class SimTK::UDotErrIndex
/// Unique integer type for Subsystem-local uDotErr indexing
/// @see SystemUDotErrIndex
SimTK_DEFINE_UNIQUE_INDEX_TYPE(UDotErrIndex);

/// @class SimTK::SystemMultiplierIndex
/// This unique integer type is for indexing global "multiplier-like" arrays, that is, arrays
/// that inherently have the same dimension as the total number of Lagrange multipliers
/// in the full System-level view of the State.
/// @see MultiplierIndex for Subsystem-local multiplier indexing
SimTK_DEFINE_UNIQUE_INDEX_TYPE(SystemMultiplierIndex);
/// @class SimTK::MultiplierIndex
/// Unique integer type for Subsystem-local multiplier indexing
/// @see SystemMultiplierIndex
SimTK_DEFINE_UNIQUE_INDEX_TYPE(MultiplierIndex);

class StateImpl;
class PerSubsystemInfo;
class DiscreteVarInfo;
class CacheEntryInfo;

/** This object is intended to contain all state information for a 
SimTK::System, except topological information which is stored in the %System 
itself. 

A SimTK::System is "const" after its topology has been constructed and realized.
Anything changeable is instead stored in a %State. Systems contain a set of 
Subsystem objects, and %State supports that concept by allowing per-subsystem 
partitioning of the total system state. This allows subsystems to have their own
private state variables, while permitting the system to allow shared access to 
state among the subsystems when necessary.

The %State provides services reflecting the structure of the equations it 
expects to find in the System. Three different views of the same state 
information are supported to accommodate three different users:
   - the system as a whole
   - Subsystems contained in the system
   - numerical methods operating on the state

Typically numerical methods have a much less nuanced view of the state
than do the system or subsystems.

The system is expected to be a "hybrid DAE", that is, a mixture of continuous 
and discrete dynamic equations, and algebraic constraints. There is an 
independent variable t, continuous state variables y, and discrete state 
variables d.

The continuous part is an ODE-on-a-manifold system suitable for solution via 
coordinate projection, structured like this for the view taken by numerical 
methods:
<pre>
     (1)  y' = f(d;t,y)         differential equations
     (2)  c  = c(d;t,y)         algebraic equations (manifold is c=0)
     (3)  e  = e(d;t,y)         event witnesses (watch for zero crossings)
</pre>
with initial conditions t0,y0,d0 such that c=0. The discrete variables d are 
updated upon occurrence of specific events. When those events are functions of 
time or state, they are detected using the set of scalar-valued event witness 
functions e (3).

In the more detailed view as seen from the System, we consider y={q,u,z} to 
be partitioned into twice-differentiable position variables q, and 
singly-differentiable velocity variables u and auxiliary variables z. There will
be algebraic constraints involving q, u, and u's time derivatives udot. The 
system is now assumed to look like this:
<pre>
     (4) qdot    = N(q) u
     (5) zdot    = zdot(d;t,q,u,z)

     (6) M(q) udot + ~G(q) mult = f(d;t,q,u,z)
         G(q) udot              = b(d;t,q,u)

                   [ pdotdot(d;t,q,u,udot) ]
     (7) udotErr = [ vdot(d;t,q,u,udot)    ] = 0
                   [ a(d;t,q,u,udot)       ] 

     (8) uErr    = [ pdot(d;t,q,u) ]         = 0
                   [ v(d;t,q,u)    ]

     (9) qErr    = [ p(d;t,q) ]              = 0
                   [ n(q)     ]
</pre>
The q's can also be dealt with directly as second order variables via
<pre>
    (10) qdotdot = Ndot(q,qdot) u + N(q) udot
</pre>

Here G = [P;V;A] with A(q) being the coefficient matrix for constraints
appearing only at the acceleration level, and V(q)=partial(v)/partial(u)
the coefficient matrix for the velocity (nonholonomic) constraints, and
P(q)=partial(pdot)/partial(u) is the coefficient matrix of the first
time derivatives of the position (holonomic) constraints.
Note that uErr in Eq 8 is assumed to include equations resulting from
differentiation of p() in Eq 9, as well as ones first introduced at the
velocity level (nonholonomic constraints), and udotErr is similarly 
built from acceleration-only constraints a() and derivatives of higher-level
constraints.

If a system allocates nq q's, nu u's, and nz z's the State will also allocate 
matching cache variables qdot, qdotdot, udot, and zdot. If mp position 
(holonomic) constraints (9), mpv velocity constraints (8) and mpva acceleration
constraints (7) are allocated, the state creates cache entries of like sizes 
qErr, uErr, udotErr. In addition room for the mpva Lagrange multipliers 'mult' 
is allocated in the cache.

In the final view, the Subsystem view, the same variables and cache entries 
exist, but only the ones allocated by that Subsystem are visible. All of a 
Subsystem's q's are consecutive in memory, as are its u's, uErr's, etc., but 
the q's are not adjacent to the u's as they are for the System's view.

The default constructor creates a %State containing no state variables and with 
its realization cache stage set to Stage::Empty. During subsystem construction,
variables and cache entries for any stage can be allocated, however *all* 
Model stage variables must be allocated during this time. At the end of 
construction, call advanceSubsystemToStage(Topology) which will put the 
subsystem at Stage::Topology. Then the subsystems realize their Model stages, 
during which variables at any stage > Model, and cache entries at any stage
>= Model can be allocated. After that call advanceSubsystemToStage(Model)
which sets the stage to Stage::Model and disallows further state allocation.

Note that there is a global Stage for the state as a whole, and individual
Stages for each subsystem. The global stage can never be higher than
the lowest subsystem stage. Global state resources are allocated when the
global Stage advances to "Model" and tossed out if that stage is
invalidated. Similarly, cache resources are allocated at stage Instance
and forgotten when Instance is invalidated. Note that subsystems will
"register" their use of the global variable pools during their own modeling
stages, but that the actual global resources won't exist until the *system* 
has been advanced to Model or Instance stage. 

<h3>Implementation note</h3>
This class is actually a "handle" that contains only a pointer to the 
hidden implementation class, StateImpl. Most of the methods here are inline
and simply forward to the implementation object; the implementations of
these must wait for the declaration of StateImpl in StateImpl.h. **/
class SimTK_SimTKCOMMON_EXPORT State {
public:
/// Create an empty State; this allocates an implementation object on the heap.
State();

/// The copy constructor has deep copy semantics; that is, this creates a new
/// *copy* of the source object, *not* a reference to the original object.
/// This makes the new %State contain a copy of the state information in the
/// source %State, copying only state variables and not the cache. If the source
/// state hasn't been realized to at least Stage::Model, then we don't copy its
/// state variables either, except those associated with the Topology stage.
State(const State& source);

/// The move constructor is very fast. The source object is left empty.
State(State&& source);

/// Copy assignment has deep copy semantics; that is, `this` %State will 
/// contain a *copy* of the source, *not* a reference into it. This makes the 
/// current %State contain a copy of the state information in the source %State,
/// copying only state variables and not the cache. If the source state hasn't
/// been realized to at least Stage::Model, then we don't copy its state
/// variables either, except those associated with the Topology stage.
State& operator=(const State&);

/// Move assignment is very fast. The source object is left in a valid but
/// undefined condition.
State& operator=(State&& source);

/// Destruct this %State object and free up the heap space it is using.
~State();

/// Restore State to default-constructed condition.
void clear();

/// Set the number of subsystems in this state. This is done during
/// initialization of the State by a System; it completely wipes out
/// anything that used to be in the State so use cautiously!
inline void setNumSubsystems(int i);

/// Set the name and version for a given subsystem, which must already
/// have a slot allocated.
inline void initializeSubsystem(SubsystemIndex, const String& name, 
                                const String& version);

/// Register a new subsystem as a client of this State. The
/// supplied strings are stored with the State but are not
/// interpreted by it. The intent is that they can be used to
/// perform "sanity checks" on deserialized States to make
/// sure they match the currently instantiated System.
/// The subsystem index (a small integer) is returned.
inline SubsystemIndex addSubsystem(const String& name, const String& version);

/// Return the number of %Subsystems known to this %State.
inline int getNumSubsystems() const;
inline const String& getSubsystemName   (SubsystemIndex) const;
inline const String& getSubsystemVersion(SubsystemIndex) const;
inline const Stage&  getSubsystemStage  (SubsystemIndex) const;

/// This returns the *global* stage for this State.
inline const Stage& getSystemStage() const;

/// If any subsystem or the system stage is currently at or
/// higher than the passed-in one, back up to the stage just prior;
/// otherwise do nothing. This is for use if you have write
/// access to the state and can invalidate even Topology and Model
/// stages which may destroy state variables. "All" here refers to
/// all Subysystems.
inline void invalidateAll(Stage);

/// If any subsystem or the system stage is currently at or
/// higher than the passed-in one, back up to the stage just prior;
/// otherwise do nothing. This const method can only be used to
/// invalidate Stage::Instance or higher. To invalidate Model or 
/// Topology stage you must have write access to the state because
/// invalidating those stages can destroy state variables in addition 
/// to cache entries. "All" here refers to all Subsystems.
inline void invalidateAllCacheAtOrAbove(Stage) const;

/// Advance a particular Subsystem's current stage by one to
/// the indicated stage. The stage is passed in just to give us a
/// chance to verify that all is as expected. You can only advance one stage at
/// a time. Advancing to Topology, Model, or Instance stage affects
/// what you can do later.
/// @see advanceSystemToStage()
inline void advanceSubsystemToStage(SubsystemIndex, Stage) const;
/// Advance the System-level current stage by one to the indicated stage.
/// This can only be done if <em>all</em> Subsystem have already been
/// advanced to this Stage.
/// @see advanceSubsystemToStage()
inline void advanceSystemToStage(Stage) const;

/** The Topology stage version number (an integer) stored in this %State must 
match the topology cache version number stored in the System for which it is 
allegedly a state. 
@see System::getSystemTopologyCacheVersion() **/
inline StageVersion getSystemTopologyStageVersion() const;

/** @name                  Continuous state variables
These continuous state variables are shared among all the subsystems
and are not allocated until the \e system is advanced to Stage::Model.
The returned index is local to each subsystem. After the system is modeled,
we guarantee that all the q's for a subsystem will be contiguous, and
similarly for u's and z's. However, q,u,z will \e not be contiguous with
each other. The \e global y={q,u,z} is contiguous, and global q,u,z are
contiguous within y, in that order. Corresponding cache entries for
the derivatives of these variables are allocated at Model stage also. **/
/**@{**/
/** Allocate generalized coordinates \e q, which are second order continuous
state variables. Matching cache entries \e qdot and \e qdotdot are allocated to 
hold the first and second time derivatives of \e q. The supplied vector \a qInit
is used to specify the number of \e q's to be allocated and their initial 
values. The Subsystem-local QIndex of the first allocated \e q is returned; the 
others follow consecutively. **/
inline QIndex allocateQ(SubsystemIndex, const Vector& qInit);
/** Allocate generalized speeds \e u, which are first order continuous
state variables related to the derivatives of the second order \e q's by
qdot=N(q)*u, for a %System-defined coupling matrix N. A matching cache entry
\e udot is allocated to hold the time derivative of \e u. The supplied vector 
\a uInit is used to specify the number of \e u's to be allocated and their 
initial values. The Subsystem-local UIndex of the first allocated \e u is 
returned; the others follow consecutively. **/
inline UIndex allocateU(SubsystemIndex, const Vector& uInit); 
/** Allocate auxiliary first order continuous state variables \e z. A matching 
cache entry \e zdot is allocated to hold the time derivative of \e z. The 
supplied vector \a zInit is used to specify the number of \e z's to be allocated
and their initial values. The Subsystem-local ZIndex of the first allocated 
\e z is returned; the others follow consecutively. **/
inline ZIndex allocateZ(SubsystemIndex, const Vector& zInit);
/**@}**/

/** @name               Constraint errors and multipliers
These constraint error cache entries are shared among all the subsystems
and are not allocated until the \e system is advanced to Stage::Instance.
The returned index is local to each subsystem. Q errors and U errors will
each be contiguous for a given subsystem, but \e not with each other. However,
the System-level yerr={qerr,uerr} \e is a single contiguous vector.
UDotErr is a separate quantity, not part of yerr. Again the UDotErr's for
each subsystem will be contiguous within the larger UDotErr Vector.
Allocating a UDotErr has the side effect of allocating another Vector
of the same size in the cache for the corresponding Lagrange multipliers,
and these are partitioned identically to UDotErrs. **/
/**@{**/
/** Allocate \a nqerr cache slots to hold the current error for position-level
(holonomic) constraint equations. **/
inline QErrIndex    allocateQErr   (SubsystemIndex, int nqerr) const;   
/** Allocate \a nuerr cache slots to hold the current error for velocity-level
(nonholonomic and holonomic first derivative) constraint equations. **/
inline UErrIndex    allocateUErr   (SubsystemIndex, int nuerr) const;
/** Allocate \a nudoterr cache slots to hold the current error for 
acceleration-level (acceleration-only, nonholonomic first derivative, and 
holonomic second derivative) constraint equations. This also allocates the
same number of slots in the constraint multipliers vector. **/
inline UDotErrIndex allocateUDotErr(SubsystemIndex, int nudoterr) const;
/**@}**/

/** @name                      Discrete Variables

You can allocate a new DiscreteVariable in any State whose stage has not yet 
been advanced to Model stage. The stage at allocation (Empty or Topology) is 
remembered so that the appropriate discrete variables can be forgotten if the 
%State's stage is reduced back to that stage later after advancing past it. 
DiscreteVariables are private to each Subsystem and allocated immediately. The 
returned index is unique within the Subsystem and there is no corresponding 
global index. **/
/**@{**/
/** The Stage supplied here in the call is the earliest subsystem stage which is
invalidated by a change made to this discrete variable. You may access the value
of the discrete variable for reading (via getDiscreteVariable()) or writing (via
updDiscreteVariable()) any time after it has been allocated. Access for writing 
has the side effect of reducing the subsystem and system stages for this State 
to one stage below the one supplied here, that is, the stage supplied here is 
invalidated. Note that you must have write access to the State in order to 
change the value of any state variable.

Ownership of the AbstractValue object supplied here is taken over by the 
State -- don't delete the object after this call!
@see getDiscreteVariable()
@see updDiscreteVariable() **/
inline DiscreteVariableIndex 
allocateDiscreteVariable(SubsystemIndex, Stage invalidates, AbstractValue*);

/** This method allocates a DiscreteVariable whose value should be updated
automatically after each time step. A CacheEntry of the same value type as
the variable is allocated to hold the update value. The discrete variable is 
allocated as described for allocateDiscreteVariable(), except that the 
\a invalidates stage must be higher than Stage::Time. The cache entry is 
allocated as described for allocateCacheEntry() without an automatic calculation
(\a latest) stage. The cache entry is then considered to be the "update" value
for the discrete variable. Update values play a similar role for
discrete variables as derivatives play for continuous variables. That is, they
define how the variable is to be updated when a TimeStepper accepts a step.

Update occurs as follows: at the start of every continuous interval, after all
other pending events have been handled, a time stepper should call the State
method autoUpdateDiscreteVariables(). That method looks at all the
auto-update discrete variables to see which ones have valid update values. For
each valid value, the discrete variable and its update value are swapped, and
the new cache value is marked invalid. 

@note No stage is invalidated by the swap even though this is clearly modifying
the state variable. It is up to the user of this variable to make sure that is
reasonable, by using the <em>update value</em>, not the <em>variable value</em>
for computations during realize(). In that way the results are always calculated 
using the value as it will be \e after an update. That means that no results 
will change when the swap occurs, so no stage needs to be invalidated upon 
updating. If you do use both values, make sure that all computed results remain
unchanged from the end of one step to the beginning of the next. 

The above behavior is entirely analogous to the treatment of continuous
variables like q: the integrator ensures that only updated values of q are
seen when evaluations are made at intermediate or trial steps; you should
do the same. In contrast to this auto-update behavior, any \e explicit change 
to the discrete variable will invalidate the variable's \a invalidates stage 
just as for a non-auto-updating discrete variable. The auto-update cache entry
is always invalidated by an explicit change to the variable, as well as by
the \a updateDependsOn stage being invalidated.

Ownership of the AbstractValue object supplied here is taken over by the 
State -- don't delete the object after this call! A clone() of this value will 
be used in the auto-update cache entry so there will be two objects of this type
around at run time that get swapped back and forth between the state variable 
and the cache entry.

You can allocate discrete variables in a State at Topology stage or Model
stage but not later. That is, you allocate the variable while the State
is in Stage::Empty, and then it appears when you do realizeTopology(); or,
you allocate the variable when the State is in Stage::Topology and it
appears when you do realizeModel().

@see allocateDiscreteVariable()
@see allocateCacheEntry() **/
inline DiscreteVariableIndex
allocateAutoUpdateDiscreteVariable(SubsystemIndex, Stage invalidates, 
                                   AbstractValue*, Stage updateDependsOn); 
/** For an auto-updating discrete variable, return the CacheEntryIndex for 
its associated update cache entry, otherwise return an invalid index. **/
inline CacheEntryIndex 
getDiscreteVarUpdateIndex(SubsystemIndex, DiscreteVariableIndex) const;
/** At what stage was this State when this discrete variable was allocated? The 
answer must be Stage::Empty or Stage::Topology. **/
inline Stage 
getDiscreteVarAllocationStage(SubsystemIndex, DiscreteVariableIndex) const;
/** What is the earliest stage that is invalidated when this discrete variable
is modified? All later stages are also invalidated. This stage was set
when the discrete variable was allocated and can't be changed with
unallocating it first. **/
inline Stage 
getDiscreteVarInvalidatesStage(SubsystemIndex, DiscreteVariableIndex) const;


/** Get the current value of the indicated discrete variable. This requires
only that the variable has already been allocated and will fail otherwise. **/
inline const AbstractValue& 
getDiscreteVariable(SubsystemIndex, DiscreteVariableIndex) const;
/** Return the time of last update for this discrete variable. **/
inline double 
getDiscreteVarLastUpdateTime(SubsystemIndex, DiscreteVariableIndex) const;
/** For an auto-updating discrete variable, return the current value of its 
associated update cache entry; this is the value the discrete variable will have
the next time it is updated. This will fail if the value is not valid or if this
is not an auto-update discrete variable. **/ 
inline const AbstractValue& 
getDiscreteVarUpdateValue(SubsystemIndex, DiscreteVariableIndex) const;
/** For an auto-updating discrete variable, return a writable reference to
the value of its associated update cache entry. This will be the value that this
discrete variable will have when it is next updated. Don't forget to mark
the cache entry valid after you have updated it. This will fail if this is
not an auto-update discrete variable. **/
inline AbstractValue& 
updDiscreteVarUpdateValue(SubsystemIndex, DiscreteVariableIndex) const;
/** Check whether the update value for this auto-update discrete variable has
already been computed since the last change to state variables it depends on. 
**/
inline bool 
isDiscreteVarUpdateValueRealized(SubsystemIndex, DiscreteVariableIndex) const;
/** Mark the update value for this auto-update discrete variable as up-to-date
with respect to the state variables it depends on. **/
inline void 
markDiscreteVarUpdateValueRealized(SubsystemIndex, DiscreteVariableIndex) const;

/** Get a writable reference to the value stored in the indicated discrete
state variable dv, and invalidate stage dv.invalidates and all higher stages.
The current time is recorded as the variable's "last update time". **/
inline AbstractValue& 
updDiscreteVariable(SubsystemIndex, DiscreteVariableIndex);
/** Alternate interface to updDiscreteVariable. **/
inline void 
setDiscreteVariable
   (SubsystemIndex, DiscreteVariableIndex, const AbstractValue&);
/**@}**/

/** @name                      Cache Entries

You can allocate a new CacheEntry in any State whose stage has not yet been
advanced to Instance stage. The stage at allocation (Empty, Topology, or 
Model) is remembered so that the appropriate cache entries can be forgotten
if the State's stage is reduced back to that stage later after advancing 
past it. CacheEntries are private to each Subsystem and allocated 
immediately. The returned index is unique within the Subsystem and there 
is no corresponding global index. **/
/**@{**/

/** There are two Stages supplied explicitly as arguments to this method: 
\a earliest and \a latest. The \a earliest Stage is the stage at which the 
cache entry \e could be calculated. Hence if the Subsystem stage is reduced
below \a earliest the cache entry is known to be invalid. The \a latest 
Stage, if any, is the stage at which the cache entry is \e guaranteed to 
have been calculated (typically as the result of a System-wide realize()
call to that stage). For stages \a earliest through \a latest-1, the 
cache entry \e may be valid, if it has already been calculated. In that 
case an explicit validity indicator will have been set at the time it was 
computed, via markCacheValueRealized(). That indicator is cleared 
automatically whenever the Subsystem stage is reduced below \a earliest. 
The validity indicator need not have been set in order for the cache entry 
to be deemed valid at \a latest stage.

If \a latest is given as Stage::Infinity then there is no guarantee that 
this Subsystem will automatically calculate a value for this cache entry,
which makes it a "lazy" evaluation that is done only if requested. In that 
case the only way the cache entry can become valid is if the calculation 
is performed and the validity indicator is set explicitly with
markCacheValueRealized(). Here is how we suggest you structure lazy
evaluation of a cache entry CE of type CEType and CacheEntryIndex CEIndex
(this is pseudocode):

(1) Allocate your lazy cache entry something like this:
\code
    CEIndex = s.allocateLazyCacheEntry(subsys,stage,new Value<CEType>());
\endcode
(2) Write a realizeCE() method structured like this:
\code
    void realizeCE(const State& s) const {
        if (s.isCacheValueRealized(subsys,CEIndex)) 
            return;
        // calculate the cache entry, update with updCacheEntry()
        s.markCacheValueRealized(subsys,CEIndex);
    }
\endcode
(3) Write a getCE() method structured like this:
\code
    const CEType& getCE(const State& s) const {
        realizeCE(s); // make sure CE has been calculated
        return Value<CEType>::downcast(s.getCacheEntry(subsys,CEIndex));
    }
\endcode
(4) Write an updCE() method like this:
\code
    CEType& updCE(const State& s) const {
        return Value<CEType>::updDowncast(s.updCacheEntry(subsys,CEIndex));
    }
\endcode

Then access CE \e only through your getCE() method. There
should be only one place in your code where isCacheValueRealized() and
markCacheValueRealized() are called for a particular cache entry. If
you do this from multiple locations there is a high probability of a bug
being introduced, especially due to later modification of the code.

Prior to the Subsystem advancing to \a earliest stage, and prior to \a latest 
stage unless the validity indicator is set, attempts to look at the value via
getCacheEntry() will throw an exception. However, you may access the cache entry
for writing via updCacheEntry() any time after stage \a earliest-1. If you 
evaluate it prior to \a latest, be sure to explicitly mark it valid.
Note that cache entries are mutable so you do not need write
access to the State in order to access a cache entry for writing.

Ownership of the AbstractValue object supplied here is taken over by the 
State -- don't delete the object after this call! 
@see getCacheEntry(), updCacheEntry()
@see allocateLazyCacheEntry(), isCacheValueRealized(), markCacheValueRealized() **/
inline CacheEntryIndex 
allocateCacheEntry(SubsystemIndex, Stage earliest, Stage latest,
                   AbstractValue* value, unsigned numExtras=0) const;

/** This is an abbreviation for allocation of a cache entry whose earliest
and latest Stages are the same. That is, this cache entry is guaranteed to
be valid if its Subsystem has advanced to the supplied Stage or later, and
is guaranteed to be invalid below that Stage. **/
inline CacheEntryIndex 
allocateCacheEntry(SubsystemIndex sx, Stage stage, 
                   AbstractValue* value, unsigned numExtras=0) const
{   return allocateCacheEntry(sx, stage, stage, value, numExtras); }

/** This is an abbreviation for allocation of a lazy cache entry. The 
\a earliest stage at which this \e can be evaluated is provided; but there is no
stage at which the cache entry will automatically be evaluated. Instead you have
to evaluate it explicitly when someone asks for it, and then call
markCacheValueRealized() to indicate that the value is available. The value
is automatically invalidated when the indicated stage \a earliest is
invalidated in the State.
@see allocateCacheEntry(), isCacheValueRealized(), markCacheValueRealized() **/
inline CacheEntryIndex 
allocateLazyCacheEntry(SubsystemIndex sx, Stage earliest, 
                       AbstractValue* value, unsigned numExtras=0) const {
    return allocateCacheEntry(sx, earliest, Stage::Infinity, 
                              value, numExtras);
}

/** At what stage was this State when this cache entry was allocated?
The answer must be Stage::Empty, Stage::Topology, or Stage::Model. **/
inline Stage 
getCacheEntryAllocationStage(SubsystemIndex, CacheEntryIndex) const;

/** Retrieve a const reference to the value contained in a particular cache 
entry. The value must be up to date with respect to the state variables it 
depends on or this will throw an exception. No calculation will be 
performed here.
@see updCacheEntry()
@see allocateCacheEntry(), isCacheValueRealized(), markCacheValueRealized() **/
inline const AbstractValue& 
getCacheEntry(SubsystemIndex, CacheEntryIndex) const;

/** Retrieve a writable reference to the value contained in a particular cache 
entry. You can access a cache entry for writing any time after it has been
allocated. This does not affect the current stage. The cache entry will
neither be invalidated nor marked valid by accessing it here.
@see getCacheEntry()
@see allocateCacheEntry(), isCacheValueRealized(), markCacheValueRealized() **/
inline AbstractValue& 
updCacheEntry(SubsystemIndex, CacheEntryIndex) const; // mutable 

/** Check whether the value in a particular cache entry has been recalculated
since the last change to the state variables it depends on. Validity can
result either from an explicit call to markCacheValueRealized() or by
this %State's stage reaching the \a latest stage specified when the cache
entry was allocated, after which the value is \e presumed valid. If this
method returns true, then you can access the value with getCacheEntry()
without getting an exception thrown.
@see allocateCacheEntry(), markCacheValueRealized(), getCacheEntry() **/
inline bool isCacheValueRealized(SubsystemIndex, CacheEntryIndex) const;

/** (Advanced) Check whether cache entry value is up to date, dealing with
extra dependencies. The number of extras must match the number provided
when the cache entry was allocated. If the cache value is up to date, its
stage version is returned for the convenience of higher-level computations
that might have a dependency on this one. Otherwise the stage version is
returned with an invalid value.
@see markCacheValueRealizedWithExtras() **/
inline bool isCacheValueRealizedWithExtras
   (SubsystemIndex, CacheEntryIndex,
    unsigned numExtras, const StageVersion* extraVersions,
    StageVersion& cacheValueVersion) const;

/** Mark the value of a particular cache entry as up to date after it has
been recalculated. This %State's current stage must be at least the
\a earliest stage as supplied when this cache entry was allocated, and
it is unnecessary to call this method if the stage has reached the
specified \a latest stage since after that we'll \e presume that the
cache entry's value has been realized. Note that if the \a latest stage
was given as Stage::Infinity then it is always necessary to call this
method prior to accessing the cache entry's value. After a cache entry
has been marked valid here, isCacheValueRealized() will return true. The
cache entry is marked invalid automatically whenever a change occurs to
a state variable on which it depends.
@see allocateCacheEntry(), isCacheValueRealized(), getCacheEntry() **/
inline void markCacheValueRealized(SubsystemIndex, CacheEntryIndex) const;

/** (Advanced) Mark a cache entry up to date after recalculating it and
include version information for extra dependencies. The number of extras must 
match the number provided when the cache entry was allocated. 
@see isCacheValueRealizedWithExtras() **/
inline void markCacheValueRealizedWithExtras
   (SubsystemIndex, CacheEntryIndex,
    unsigned numExtras, const StageVersion* extraVersions) const;

/** (Advanced) Normally cache entries are invalidated automatically, however 
this method allows manual invalidation of the value of a particular cache entry.
After a cache entry has been marked invalid here, isCacheValueRealized() or
isCacheValueRealizedWithExtras() will return false.
@see isCacheValueRealized(), markCacheValueRealized() **/
inline void markCacheValueNotRealized(SubsystemIndex, CacheEntryIndex) const;
/**@}**/

/// @name                Global Resource Dimensions
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
inline int getNY() const;
/// Get total number of shared q's (generalized coordinates; second order
/// state variables). This is also the number of first and second q time
/// derivatives in the cache entries qdot and qdotdot.
/// Callable at Model stage.
inline int getNQ() const;
/// Returns the y index at which the q's begin. Callable at Model stage.
inline SystemYIndex getQStart() const;
/// Get total number of shared u's (generalized speeds; mobilities). 
/// This is also the number of u time derivatives in the cache entry udot.
/// Callable at Model stage.
inline int getNU() const;
/// Returns the y index at which the u's begin. Callable at Model stage.
inline SystemYIndex getUStart() const;
/// Get total number of shared z's (auxiliary state variables). 
/// This is also the number of z time derivatives in the cache entry zdot.
/// Callable at Model stage.
inline int getNZ() const;
/// Returns the y index at which the z's begin. Callable at Model stage.
inline SystemYIndex getZStart() const;
/// Get the total number nyerr=nqerr+nuerr of shared cache entries for
/// position-level and velocity-level constraint errors.
/// Callable at Instance stage.
inline int getNYErr() const;
/// Return the total number nqerr=mp+nQuaternions of cache entries for
/// position-level constraint errors. Callable at Instance stage.
inline int getNQErr() const;
/// Returns the yErr index at which the qErr's begin. Callable at Instance stage.
inline SystemYErrIndex getQErrStart() const; 
/// Return the total number nuerr=mp+mv of cache entries for
/// velocity-level constraint errors (including also errors in the 
/// time derivatives of position-level constraints). Callable at Instance stage.
inline int getNUErr() const;
/// Returns the yErr index at which the uErr's begin. Callable at Instance stage.
inline SystemYErrIndex getUErrStart() const; 
/// Return the total number nudotErr=mp+mv+ma of cache entries for
/// acceleration-level constraint errors (including also errors in the 
/// second time derivatives of position-level constraints and the first
/// time derivatives of velocity-level constraints). Callable at Instance stage.
inline int getNUDotErr() const;
/// Return the total number of constraint multipliers; necessarily the same
/// as the number of acceleration-level constraint errors nUDotErr. Callable
/// at Instance stage.
/// @see getNUDotErr()
inline int getNMultipliers() const; // =mp+mv+ma, necessarily the same as NUDotErr
/// Returns a mutex that should be used to lock the state whenever multiple 
/// threads are asynchronously writing/updating a common state cache. A lock
/// should always be used when thread-safe state is not guarenteed. If multiple
/// threads are simply reading from the cache, locking the state may not be
/// necessary.
inline std::mutex& getStateLock() const;

/// @}

/// @name                   Per-Subsystem Dimensions
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

inline SystemQIndex getQStart(SubsystemIndex) const; 
inline int getNQ(SubsystemIndex) const;
inline SystemUIndex getUStart(SubsystemIndex) const; 
inline int getNU(SubsystemIndex) const;
inline SystemZIndex getZStart(SubsystemIndex) const; 
inline int getNZ(SubsystemIndex) const;

inline SystemQErrIndex getQErrStart(SubsystemIndex) const; 
inline int getNQErr(SubsystemIndex) const;
inline SystemUErrIndex getUErrStart(SubsystemIndex) const; 
inline int getNUErr(SubsystemIndex) const;
inline SystemUDotErrIndex getUDotErrStart(SubsystemIndex) const; 
inline int getNUDotErr(SubsystemIndex) const;
inline SystemMultiplierIndex getMultipliersStart(SubsystemIndex) const;
inline int getNMultipliers(SubsystemIndex) const;

/// Per-subsystem access to the global shared variables.
inline const Vector& getQ(SubsystemIndex) const;
inline const Vector& getU(SubsystemIndex) const;
inline const Vector& getZ(SubsystemIndex) const;

inline const Vector& getUWeights(SubsystemIndex) const;
inline const Vector& getZWeights(SubsystemIndex) const;

inline Vector& updQ(SubsystemIndex);
inline Vector& updU(SubsystemIndex);
inline Vector& updZ(SubsystemIndex);

inline Vector& updUWeights(SubsystemIndex);
inline Vector& updZWeights(SubsystemIndex);

/// Per-subsystem access to the shared cache entries.
inline const Vector& getQDot(SubsystemIndex) const;
inline const Vector& getUDot(SubsystemIndex) const;
inline const Vector& getZDot(SubsystemIndex) const;
inline const Vector& getQDotDot(SubsystemIndex) const;

inline Vector& updQDot(SubsystemIndex) const;    // these are mutable
inline Vector& updUDot(SubsystemIndex) const;
inline Vector& updZDot(SubsystemIndex) const;
inline Vector& updQDotDot(SubsystemIndex) const;

inline const Vector& getQErr(SubsystemIndex) const;
inline const Vector& getUErr(SubsystemIndex) const;
inline const Vector& getUDotErr(SubsystemIndex) const;
inline const Vector& getMultipliers(SubsystemIndex) const;

inline const Vector& getQErrWeights(SubsystemIndex) const;
inline const Vector& getUErrWeights(SubsystemIndex) const;

inline Vector& updQErr(SubsystemIndex) const;    // these are mutable
inline Vector& updUErr(SubsystemIndex) const;
inline Vector& updUDotErr(SubsystemIndex) const;
inline Vector& updMultipliers(SubsystemIndex) const;

inline Vector& updQErrWeights(SubsystemIndex);
inline Vector& updUErrWeights(SubsystemIndex);

/// You can call these as long as *system* stage >= Model.
inline const double& getTime() const;
inline const Vector& getY() const; // {Q,U,Z} packed and in that order

/// These are just views into Y.
inline const Vector& getQ() const;
inline const Vector& getU() const;
inline const Vector& getZ() const;


/** Get a unit weighting (1/unit change) for each u that can be used to 
weight a vector du so that the disparate elements are comparable in physical
effect. This permits mixing of generalized speeds
that have different units, and scaling of generalized speeds that have
differing amounts of leverage due to their positions in the multibody tree.
This can be used to create a scaled norm that represents the overall
significance of a change du to u.

Define a unit change di for each ui such that a change
ui+eps*di to each generalized speed in turn produces a physical velocity change
of roughly equal significance. Then a diagonal matrix Wu=diag(1/di) is 
a weighting matrix such that wdu=Wu*du is a vector in which each element wdu_i
has units of "unit change" for its corresponding ui. This method returns a
vector which is the diagonal of Wu.

These same weights on u also determine the scaling of the generalized
coordinates q, because q and u are related via qdot=N*u. For cases where
qdot_i=u_i, the numerical value of the unit change to q_i is just di because 
dP/dq_i == dV/du_i. Otherwise, they are related by Wq = N*Wu*pinv(N) where 
Wq is the weighting matrix for dq (block diagonal), and pinv() is the
pseudoinverse.

For example, say you define unit scaling for an angle coordinate to be 1 radian
(about 57 degrees), meaning that a 1 radian change of coordinate produces
(roughly) one length unit of meaningful position change. Then if a generalized
coordinate is measured in radians, its unit scale would be 1. If instead you 
created a generalized coordinate with units of degrees, its unit scale would 
be 57 degrees. That would allow mixing of such coordinates in the same system 
by bringing the coordinates into a physically-meaningful basis.
Scaling is defined in the u basis where each variable is independent;
the N matrix couples variables in the q basis. So here the units would actually
be 1 radian/time unit and 57 degrees/time unit (numerically identical).

This is allocated and set to 1 at the end of realize(Model). **/
inline const Vector& getUWeights() const;    // diag(Wu)

/** Get a unit weighting (1/unit change) for each z that can be used to 
weight a vector dz so that the disparate elements are comparable in physical
effect. This defines a weighting matrix Wz=diag(1/unitchange_zi) such
that wdz=Wz*dz is a vector in which each element wdz_i has units of
"unit change" for its corresponding zi.  This method returns a
vector which is the diagonal of Wz. **/
inline const Vector& getZWeights() const;

/** Set u weights (and q weights indirectly). You can call this after Model 
stage has been realized. This will invalidate just Report stage because it is 
not used in calculating udots. **/
inline Vector& updUWeights();

/** Set z weights. You can call this after Model stage has been realized. This
will invalidate just Report stage because it is not used in calculating 
zdots. **/
inline Vector& updZWeights();

/// You can call these as long as System stage >= Model, but the
/// stage will be backed up if necessary to the indicated stage.
inline double& updTime();  // Back up to Stage::Time-1
inline Vector& updY();     // Back up to Stage::Dynamics-1

/// An alternate syntax equivalent to updTime() and updY().
inline void setTime(double t);
inline void setY(const Vector& y);

/// These are just views into Y.
inline Vector& updQ();     // Back up to Stage::Position-1
inline Vector& updU();     // Back up to Stage::Velocity-1
inline Vector& updZ();     // Back up to Stage::Dynamics-1

/// Alternate interface.
inline void setQ(const Vector& q);
inline void setU(const Vector& u);
inline void setZ(const Vector& z);

inline const Vector& getYDot()    const; // Stage::Acceleration

/// These are just views into YDot.
inline const Vector& getQDot()    const; // Stage::Velocity
inline const Vector& getZDot()    const; // Stage::Dynamics
inline const Vector& getUDot()    const; // Stage::Acceleration

/// This has its own space, not a view.
inline const Vector& getQDotDot() const; // Stage::Acceleration

/// These are mutable
inline Vector& updYDot() const;    // Stage::Acceleration-1
inline Vector& updQDot() const;    // Stage::Velocity-1     (view into YDot)
inline Vector& updZDot() const;    // Stage::Dynamics-1            "
inline Vector& updUDot() const;    // Stage::Acceleration-1        "

/// This is a separate shared cache entry, not part of YDot. If you
/// have a direct 2nd order integrator you can integrate QDotDot
/// (twice) to get Q.
inline Vector& updQDotDot() const; // Stage::Acceleration-1

/// Return the current constraint errors for all constraints. This
/// is {QErr,UErr} packed and in that order.
inline const Vector& getYErr() const;  // Stage::Velocity

/// These are just views into YErr.
inline const Vector& getQErr() const;  // Stage::Position (index 3 constraints)
inline const Vector& getUErr() const;  // Stage::Velocity (index 2 constraints)

/// These have their own space, they are not views.
inline const Vector& getUDotErr()     const; // Stage::Acceleration (index 1 constraints)
inline const Vector& getMultipliers() const; // Stage::Acceleration

/** Get the unit weighting (1/unit error) for each of the mp+mquat position 
inline constraints equations. Allocated and initialized to 1 on realize(Instance). **/
inline const Vector& getQErrWeights() const;

/** Get the unit weighting (1/unit error) for each of the mp+mv velocity-level 
inline constraint equations, meaning mp time derivatives of position (holonomic) 
inline constraint equations followed by mv velocity (nonholonomic) constraints.
Typically the weight of position constraint derivatives is just the
position constraint weight times the System's characteristic time scale. 

There is no entry corresponding to quaternions here since they do not 
produce velocity-level constraints in Simbody's forumulation.

This is allocated and initialized to 1 on realize(Instance). **/
inline const Vector& getUErrWeights() const;

/** Set the unit weighting (1/unit error) for each of the mp+mquat position 
inline constraint equations. You can call this after the weight variable is allocated 
at the end of Instance stage. Position stage is invalidated to force 
recalculation of weighted position constraint errors. **/
inline Vector& updQErrWeights();

/** Set the unit weighting (1/unit error) for each of the mp+mv velocity-level
inline constraints. You can call this after the weight variable is allocated at the 
end of Instance stage. Velocity stage is invalidated to force recalculation of 
weighted velocity-level constraint errors. **/
inline Vector& updUErrWeights();

/// These are mutable
inline Vector& updYErr() const; // Stage::Velocity-1
inline Vector& updQErr() const; // Stage::Position-1 (view into YErr)
inline Vector& updUErr() const; // Stage::Velocity-1        "

inline Vector& updUDotErr()     const; // Stage::Acceleration-1 (not a view)
inline Vector& updMultipliers() const; // Stage::Acceleration-1 (not a view)

/** @name                 Advanced/Obscure/Debugging
Don't call these methods unless you know what you're doing. **/
/**@{**/

/** (Advanced) Record the current version numbers of each valid System-level 
stage. This can be used to unambiguously determine what stages have been 
changed by some opaque operation, even if that operation realized the stages 
after modifying them. This is particularly useful for event handlers as a way 
for a time stepper to know how much damage may have been done by a handler, and
thus how much reinitialization is required before continuing on.
@see getLowestSystemStageDifference() **/
inline void getSystemStageVersions(Array_<StageVersion>& versions) const;

/** (Advanced) Given a list of per-stage version numbers extracted by an 
earlier call to getSystemStageVersions(), note the lowest system stage in the 
current State whose version number differs from the corresponding previous 
version number. Returns Stage::Infinity if all the stages present in 
\a prevVersions are valid and have identical versions now, even if there are 
additional valid stages now, since nothing the caller cared about before has 
been changed. If the current State is not realized as far as the previous one, 
then the first unrealized stage is returned if all the lower versions match.
@see getSystemStageVersions() **/
inline Stage getLowestSystemStageDifference
   (const Array_<StageVersion>& prevVersions) const;

/** (Advanced) This explicitly modifies the Topology stage version; don't
use this method unless you know what you're doing! This can be used to force
compatibility with a System that has had Topology changes since this %State
was created. This has no effect on the realization level.
@see getSystemTopologyStageVersion(), System::getSystemTopologyCacheVersion()
**/
inline void setSystemTopologyStageVersion(StageVersion topoVersion);

/** (Advanced) Return a StageVersion for time, meaning an integer that is
incremented whenever time is changed. Will be 1 or greater if time has been
initialized; otherwise, time is meaningless and may be NaN. **/
inline StageVersion getTimeStageVersion() const;
/** (Advanced) Return a StageVersion for q, meaning an integer that is
incremented whenever any q is changed (or more precisely whenever q or y is
returned with writable access). Will be 1 or greater if q has been
allocated. **/
inline StageVersion getQStageVersion() const;
/** (Advanced) Return a StageVersion for u, meaning an integer that is
incremented whenever any u is changed (or more precisely whenever u or y is
returned with writable access). Will be 1 or greater if u has been
allocated. **/
inline StageVersion getUStageVersion() const;
/** (Advanced) Return a StageVersion for z, meaning an integer that is
incremented whenever any z is changed (or more precisely whenever z or y is
returned with writable access). Will be 1 or greater if z has been
allocated. **/
inline StageVersion getZStageVersion() const;

/** (Advanced) Return a reference to the cache entry information for a 
particular cache entry. No validity checking is performed. **/
SimTK_FORCE_INLINE const CacheEntryInfo& 
getCacheEntryInfo(SubsystemIndex, CacheEntryIndex) const;

/** (Advanced) Return a reference to the discrete variable information for a 
particular discrete variable. **/
SimTK_FORCE_INLINE const DiscreteVarInfo& 
getDiscreteVarInfo(SubsystemIndex, DiscreteVariableIndex) const;

/** (Advanced) Return a reference to the per-subsystem information in the
state. **/
SimTK_FORCE_INLINE const PerSubsystemInfo& 
getPerSubsystemInfo(SubsystemIndex) const;

/** (Advanced) This is called at the beginning of every integration step to set
the values of auto-update discrete variables from the values stored in their 
associated cache entries. **/
inline void autoUpdateDiscreteVariables();

/** (Debugging) Not suitable for serialization. **/
inline String toString() const;
/** (Debugging) Not suitable for serialization. **/
inline String cacheToString() const;
/**@}**/

//------------------------------------------------------------------------------
// The implementation class and associated inline methods are defined in a
// separate header file included below.
                                private:
StateImpl* impl;

const StateImpl& getImpl() const {assert(impl); return *impl;}
StateImpl&       updImpl()       {assert(impl); return *impl;}
};

// Dump state and cache to a stream for debugging; this is not serialization.
SimTK_SimTKCOMMON_EXPORT std::ostream& 
operator<<(std::ostream& o, const State& s);

} // namespace SimTK

// This completes the inline definition of State.
#include "SimTKcommon/internal/StateImpl.h"

#endif // SimTK_SimTKCOMMON_STATE_H_
