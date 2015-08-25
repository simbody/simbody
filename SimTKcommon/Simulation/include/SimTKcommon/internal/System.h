#ifndef SimTK_SimTKCOMMON_SYSTEM_H_
#define SimTK_SimTKCOMMON_SYSTEM_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-15 Stanford University and the Authors.        *
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

#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"
#include "SimTKcommon/internal/State.h"
#include "SimTKcommon/internal/Subsystem.h"
#include "SimTKcommon/internal/SubsystemGuts.h"
#include "SimTKcommon/internal/Event.h"
#include "SimTKcommon/internal/EventTrigger.h"
#include "SimTKcommon/internal/EventTrigger_Timer.h"
#include "SimTKcommon/internal/EventTrigger_Witness.h"

#include <cassert>

namespace SimTK {

class DecorativeGeometry;
class SystemGlobalSubsystem;
class ScheduledEventHandler;
class ScheduledEventReporter;
class TriggeredEventHandler;
class TriggeredEventReporter;
class RealizeOptions;
class RealizeResults;
class ProjectOptions;
class ProjectResults;

/** @class SimTK::ActiveWitnessIndex
This unique integer type is for indexing the array of active event witnesses
returned by System::findActiveEventWitnesses() for a particular State. This 
numbering is valid only for that state and will change meaning every time 
witnesses are added or removed. **/
SimTK_DEFINE_UNIQUE_INDEX_TYPE(ActiveWitnessIndex);

/** @class SimTK::ActiveTimerIndex
This unique integer type is for indexing the array of active event timers
returned by System::findActiveEventTimers() for a particular State. This 
numbering is valid only for that state and will change meaning every time timers
are added or removed. **/
SimTK_DEFINE_UNIQUE_INDEX_TYPE(ActiveTimerIndex);


//==============================================================================
//                                 SYSTEM
//==============================================================================
/** This is the base class that serves as the parent of all SimTK %System 
objects; most commonly Simbody's MultibodySystem class. %System is also the
object type on which the Simbody numerical integrators and other solvers 
operate, and many of the methods provided here are primarily for use by 
those solvers.

A %System serves as a mediator for a group of interacting Subsystem objects. 
All will share a single system State, and typically subsystems will need access
to content in the state which is managed by other subsystems. A %System provides 
a unique SubsystemIndex (a small positive integer) for each of its subsystems, 
and the subsystems are constructed knowing their indices. The indices are used 
subsequently by the subsystems to find their own entries in the system state, 
and by each subsystem to refer to others within the same system. Index 0 is 
reserved for use by the %System itself, e.g. for system-global state variables 
and events. Those are maintained in the SystemGlobalSubsystem, which has 
SubsystemIndex 0. 

Concrete %Systems understand the kinds of subsystems they contain. For example, 
a MultibodySystem might contain a mechanical subsystem, some force subsystems, 
and a geometry subsystem. At each computation stage, a subsystem is realized 
in a single operation. That operation can refer to computations from 
already-realized subsystems, but cannot initiate computation in other 
subsystems. The %System must know the proper order in which to realize the 
subsystems at each stage, and that ordering is likely to vary with stage. For 
example, at Position stage the mechanical positions must be realized before 
the configuration-dependent force elements. However, at Acceleration stage, 
the force elements must be realized before the mechanical accelerations can 
be calculated.

There are three distinct users of this class:
  - %System Users: people who are making use of a concrete %System (which will
    inherit methods from this class).
  - Numerical solvers: TimeStepper and Integrator objects require a 
    %System object whose State can be advanced through time.
  - %System Developers: people who are writing concrete System classes.

Note that %System Users may include people who are writing studies, reporters, 
modelers and so on as well as end users who are accessing the %System directly.

Methods intended for %System Users, numerical solvers, and a few bookkeeping 
methods are in the main %System class, which is a SimTK Handle class. A Handle
class consists only of a single pointer, which in this case points to a 
System::Guts object. The System::Guts class is abstract, and virtual methods to
be implemented by %System Developers in the concrete %System are defined there, 
along with other utilities of use to the concrete %System Developer but not to 
the end user or numerical solvers. The System::Guts class is declared 
in a separate header file, and only people who are writing their own %System
classes need look there. **/
class SimTK_SimTKCOMMON_EXPORT System {
public:
class Guts; // local; name is System::Guts


//------------------------------------------------------------------------------
/**@name                    System-wide properties

These methods set parameters that provide %System-level defaults that can be
used by various algorithms. This includes nominal time and length scales, and
hints that can be useful when visualizing this %System (for example, does
it make sense to show a ground plane?). There are defaults for all of these
parameters so you don't need to set them unless the defaults are not 
satisfactory. **/
/**@{**/

/** This is a hint to other software as to which way this System's designer 
considers to be "up". For visualization, this is the best direction to use as 
the default up direction for the camera, and the opposite direction is likely to
be a good direction in which to apply gravitational forces. The default up 
direction is +YAxis, which is the same as the OpenGL convention for the camera 
up direction. You can set this to any of the coordinate axes in the positive or 
negative direction. For example, use setUpDirection(ZAxis) for the "virtual 
world" convention where ground is the x-y plane and setUpDirection(-ZAxis) 
might be appropriate for an aviation convention in which Z is directed downward.
A visualizer that is showing a ground plane should make the ground plane normal 
be this up direction.
@see setUseUniformBackground() **/
System& setUpDirection(const CoordinateDirection& up);

/** This is a hint to visualization software that this System is best viewed 
against a uniform background (e.g.\ all white) rather than against a ground 
plane. A molecular system will typically set this flag so that the visualizer
will not attempt to place the molecule on the ground. The default is to 
consider this system best viewed with a ground plane displayed, perpendicular 
to the "up" direction and located at a height of zero.
@see setUpDirection() **/
System& setUseUniformBackground(bool useUniformBackground);

/** (Advanced) This is a hint used for some default behaviors, such as 
determining an initial step size for an integrator, the default time window for
event localization, and the default unit tolerance for a constraint error 
function derivative from the original constraint's unit tolerance. Most users 
can ignore this and just take the default.

This should be set to roughly the time scale at which you expect to see 
interesting things happen, that is the scale at which you might choose to view 
reporting output. An orbital simulation using seconds as time units might set 
this to 10 or 100s, for example, while a biomechanical simulation could use 
0.1s. This will affect the time scale on which velocity constraints are
stabilized, with longer time scales being more demanding since there is more 
time to drift. By default this is 0.1 time units, so 100ms for systems measuring
time in seconds and 100fs for systems measuring time in ps. **/
System& setDefaultTimeScale(double timeScale);

/** (Advanced) This is a hint that can be used to get a sense of what a "unit 
length" looks like for this System in the units it uses. Most users can ignore
this and just take the default.

For example, a model of a small toy car expressed in MKS units might set this to 
0.01 (1 cm). The default for this is 1 length unit, meaning 1 meter in MKS and 1 
nm in MD units. **/
System& setDefaultLengthScale(Real unitLength);

/** This determines whether this System wants to be notified whenever time
advances irreversibly. If set true, time advancement is treated as an
event, and the handleEvents() method is invoked with its \a cause argument
set to indicate a time-advanced event occurred. By default, time advancement 
proceeds silently. **/
void setHasTimeAdvancedEvents(bool); // default=false

/** Get the current setting of the "up" direction hint. **/
CoordinateDirection getUpDirection() const;
/** Get the current setting of the "use uniform background" visualization
hint. **/
bool getUseUniformBackground() const;
/** Get the currently-set value for the default time scale, 0.1 time units
if it has never been set. **/
double getDefaultTimeScale() const;
/** Get the currently-set value for the default length scale, 1 length unit
if it has never been set. **/
Real getDefaultLengthScale() const;
/** Return the current value of the flag indicating whether this %System wants
an event generated whenever time advances irreversibly. **/
bool hasTimeAdvancedEvents() const;
/**@}**/



//------------------------------------------------------------------------------
/**@name                         Realization

These methods provide the ability to compute the consequences that follow from 
the values of state variables, leaving the results in the appropriate cache. 
Usually this refers to values in a given State object, but we also use the same
concept for values assigned to %System member variables during the %System's
(extended) construction, which we call *Topological* variables. The 
realizeTopology() method is called at the end of construction, allocating
resources, performing computations, and creating a default State, with all
results stored in a cache that is maintained by the %System itself. The other 
realize() methods are given a State object to work with and store their results 
in the cache of that State object without making any changes to the %System. **/
/**@{**/

/** The following call must be made after any topological change has been made 
to this %System, before the %System can be used to perform any computations. 
Perhaps surprisingly, the method is const. That's because the topology cannot be
changed by this method. Various mutable "cache" entries will get calculated, 
including the default State, a reference to which is returned. The returned 
State has a Topology stage version number that matches the %System, and will 
have already been realized through the Model Stage, using the defaults for the 
Model-stage variables, meaning that all later stage variables have been 
allocated and set to their default values as well. You can access this same 
default State again using getDefaultState(). If the current topology has already
been realized, this call does nothing but return a reference to the 
already-built default State. 
@see getDefaultState(), realizeModel(), realize() **/
const State& realizeTopology() const;

/** This is available after realizeTopology(), and will throw an exception if 
realizeTopology() has not been called since the most recent topological change
to this %System. This method returns the same reference returned by 
realizeTopology(). The State to which a reference is returned was created by the
most recent realizeTopology() call. It has a Topology version number that 
matches the one currently in this %System, and has already been realized through
the Model Stage, using default values for all the Model-stage variables. All 
later-stage variables have been allocated and set to their default values. You 
can use this state directly to obtain information about the %System in its 
default state or you can use this state to initialize other States to which you
have write access. Those States are then suitable for further computation with 
this %System. 
@see realizeTopology() **/
const State& getDefaultState() const;

/** Realize the model to be used for subsequent computations with the given
`state`. This call is required if Model-stage variables are changed from their 
default values. The %System topology must already have been realized (that is, 
realizeTopology() must have been called since the last topological change made 
to the %System). Also, the supplied `state` must already have been initialized
to work with this %System either by copying the default state or some other 
State of this %System. If it has already been realized to Stage::Model or 
higher, nothing happens here. Otherwise, all Model-stage State resources are
acquired and initialized. This method modifies its `state` argument, but makes 
no changes at all to the %System itself (not even to the %System's mutable 
cache) and is hence const.

@note Any Model-stage state resources, including both state variables and cache
entries, are destroyed when a Model-stage state variable is modified. Then, when
realizeModel() is called they are acquired again. The number, types and memory 
locations of those resources will change, so any existing references or pointers
to them will be invalid.  
@see realizeTopology(), realize() **/
void realizeModel(State& state) const;

/** Realize the given `state` to the indicated `stage`. The passed-in State 
must have been initialized to work with this %System, and it must already have 
been realized through Stage::Model, since this realize() method doesn't have 
write access to the State. If the state has already been realized to the 
requested stage or higher, nothing happens here. Otherwise, the state is 
realized one stage at a time until it reaches the requested stage. 
@see realizeTopology(), realizeModel() **/
void realize(const State& state, Stage stage = Stage::HighestRuntime) const;
/**@}**/


//------------------------------------------------------------------------------
/**@name                    The Constrained System

These methods deal with prescribed motion and constraints on the allowable 
values of q and u. They are primarily for use by numerical integrators. 
Prescribed motion is always resolved exactly because it involves explicit
computation of q=q(t) and u=u(t,q). Constraints are given by implicit
equations like g(q)=0 and are solved to a specified tolerance. Also, prescribed
motion has precedence over constraints in the sense that no changes will be
made to prescribed variables in order to satisfy constraints; only the free
variables are modified by constraint solvers. **/
/**@{**/


/** Apply prescribed motion and attempt to satisfy all position and velocity 
constraints by projecting the generalized coordinates q and generalized speeds 
u onto the position and velocity constraint manifolds. This method expects to 
be given a \a state that is close to satisfying the constraints already. It 
will apply prescribed motion and then attempt a least-squares projection to 
satisfy constraints, working first on q's and then on u's. If you are trying to
assemble a system from an arbitrary configuration, would like control over
which states get modified, or are attempting to track markers, use an 
Assembler solver instead.

@param[in,out]  state   
    A State whose q's and u's are to be modified. Must already have been 
    realized to at least Stage::Model. On return will be realized to 
    Stage::Velocity.
@param[in]      accuracy    
    If provided, sets the required accuracy to which the constraints must
    be satisfied; an exception will be thrown if that accuracy cannot be 
    achieved. If accuracy is left off or is nonpositive, the accuracy used is
    the default as provided by the ProjectOptions class.

This method is provided for convenience and involves several calls to realize(),
prescribeQ(), prescribeU(), projectQ(), and projectU(); call those directly if 
you want more control.

On return the given \a state's q and u variables will have been modified and
the \a state will be realized through Velocity stage.


<h3>Theory</h3>

This solver projects the given state back on to the position or velocity
constraint manifold, by the shortest path possible in the scaled norm defined in
that state, and modifying only free (non-prescribed) variables. Constraint 
errors are scaled by their unit error weightings, then satisfied to a given 
accuracy using an RMS norm, or optionally using the stricter infinity norm.
You can specify accuracy here explicitly; otherwise it uses the default as
provided by the ProjectOptions class.

Position constraints are satisfied as follows:
<pre>
    solve |Tp*perr(t;q+dq)|_n <= accuracy for dq (n==rms or inf)
    such that |Wq*dq|_2 is minimized
</pre>
Here Tp=diag(1./unit_p) scales each position constraint error to a fraction of
its unit error. Wq=N*Wu*pinv(N) weights dq to include both the "unit change"
weightings on u and the artifactual configuration-dependent weightings on q
generated by choice of orientation coordinates such as quaternions or rotation
angles. (N as in qdot=N*u; pinv() is pseudoinverse.) We do not allow relative 
weighting on dq based on the current values of q; Simbody always solves q to 
absolute accuracy since arbitrary translations and rotations by 2pi should not 
affect physically-significant results.

Velocity constraints are satisfied as follows:
<pre>
    solve |Tpv*pverr(t,q;u+du)|_n <= accuracy for du (n==rms or inf)
    such that |Eu*du|_2 is minimized
</pre> 
Here Tpv=diag(ts./unit_p 1./unit_v) where ts is the system's time scale used to 
scale the holonomic constraint's unit errors, and unit_v is the unit error for 
the holonomic constraints. The error weighting matrix Eu combines relative
and absolute accuracy requirements as follows:
<pre>
    Eu_i={ min(Wu_i, 1/u_i), relative accuracy OK for u_i
         {       Wu_i,       otherwise
</pre>

@see Assembler, projectQ(), projectU(), prescribeQ(), prescribeU() **/
void project(State& state, Real accuracy=-1) const;

/** Set prescribed positions q=q(t) and attempt to satisfy position constraints
by modifying the remaining free q's. The given \a state will be realized as
necessary and on return will have been realized through Position stage. See 
project() for more information.

@param[in,out]  state   
    A State whose q's are to be modified. Must already have been realized to at
    least Stage::Model. On return will be realized to Stage::Position.
@param[in]      accuracy    
    If provided, sets the required accuracy to which the constraints must
    be satisfied; an exception will be thrown if that accuracy cannot be 
    achieved. If accuracy is left off or is nonpositive, the accuracy used is
    the default as provided by the ProjectOptions class.
      
@see project(), projectU(), ProjectOptions **/
void projectQ(State& state, Real accuracy=-1) const;

/** Set prescribed velocities u=u(t,q) and attempt to satisfy velocity
constraints by modifying the remaining free u's. The given \a state will be 
realized as necessary and on return will have been realized through Velocity 
stage. Note that no q's are modified by this method; you should already have
ensured that they satisfy position constraints. See project() for more 
information.

@param[in,out]  state   
    A State whose u's are to be modified. Must already have been realized to at
    least Stage::Model. On return will be realized to Stage::Velocity.
@param[in]      accuracy    
    If provided, sets the required accuracy to which the constraints must
    be satisfied; an exception will be thrown if that accuracy cannot be 
    achieved. If accuracy is left off or is nonpositive, the accuracy used is
    the default as provided by the ProjectOptions class.
        
@see project(), projectQ(), ProjectOptions **/
void projectU(State& state, Real accuracy=-1) const;

/** (Advanced) Project free q's so that position constraints are satisfied and 
remove the corresponding error from the supplied error estimate. This is 
primarily intended for use by numerical integration algorithms. You must 
already have set prescribed q's; this method will not modify them but may 
depend on their current values. State must be realized to Position stage on 
entry and will still be realized through Position stage on return. 

If the norm of the position constraint errors is already less than or equal to
the accuracy requested in \a options on entry, nothing will happen unless you 
have selected the "force projection" option. You can find out what actually 
happened by looking in the returned \a results.

Optionally, this method can also project out the constraint-normal component of
the passed-in error estimate vector \a qErrEst. This is part of the integration
of the continuous DAE system and thus should never require an integrator 
restart. Just pass a 0-length (default constructed) Vector in place of the
error estimate if you don't want error projection.

Some options are available, see ProjectOptions for complete list. For example,
- use infinity norm
- local projection only
- force projection
- don't throw an exception
- force full Newton

This method and projectU() are not for satisfying acceleration constraints, 
which does not involve modifications to state variables. Acceleration 
constraints are satisfied automatically when you realize a state to 
Acceleration stage using realize(state); the resulting udots will satisfy the 
acceleration constraints (if possible), even if the position and velocity 
constraints are not satisfied.

See the Theory section in the documentation for project() for more information.

@see ProjectOptions, ProjectResults, projectU(), project() **/
void projectQ(State& state, Vector& qErrEst, 
             const ProjectOptions& options, ProjectResults& results) const;

/** (Advanced) Project free u's so that velocity constraints are satisfied and 
remove the corresponding error from the supplied error estimate. This is 
primarily intended for use by numerical integration algorithms. You must 
already have dealt with q's and already set prescribed u's; this method will 
not modify them but may depend on their current values. State must be realized 
to Velocity stage on entry and will still be realized through Velocity stage on 
return. 

If the norm of the velocity constraint errors is already less than or equal to 
the accuracy requested in \a options on entry, nothing will happen unless you 
have selected the "force projection" option. You can find out what actually 
happened by looking in the returned \a results.

Optionally, this method can also project out the constraint-normal component of
the passed-in error estimate vector \a uErrEst. This is part of the integration
of the continuous DAE system and thus should never require an integrator 
restart. Just pass a 0-length (default constructed) Vector in place of the
error estimate if you don't want error projection.

See the Theory section in the documentation for project(), and the comments
for projectQ(), for more information.

@see ProjectOptions, ProjectResults, projectQ(), project() **/
void projectU(State& state, Vector& uErrEst, 
             const ProjectOptions& options, ProjectResults& results) const;


/** Set values for prescribed positions q and velocities u.
Prescribed positions are functions of time q(t) and prescribed velocities are 
functions of time and position u(t,q). Both can also depend on earlier-stage 
discrete variables such as modeling and instance parameters.

@param[in,out]  state
    The State to be modified. Time and the values of non-prescribed q's are 
    obtained from \a state and prescribed q's and u's are modified on return.
    The \a state will be realized as needed and on return will have been
    realized through Position stage. The prescribed velocities will have been
    set but not yet realized.

This is a convenience method that performs several realize() calls and uses
prescribeQ() and prescribeU(). If you want more control, use those methods
directly instead of this one.

@see prescribeQ(), prescribeU(), project(), realize() **/
void prescribe(State& state) const {
    realize(state, Stage::Time);
    prescribeQ(state);
    realize(state, Stage::Position);
    prescribeU(state);
}

/** This solver sets prescribed q's to their known values q=q(t) as a function
of time (and possibly earlier-stage discrete state variables).

We expect the supplied State already to have been realized to Stage::Time. 
Note that the \e derivatives qdot of prescribed q's (which are
of necessity also prescribed but are not themselves state variables) are set in
the subsequent realize() call. For example, prescribeU sets the 
prescribed u's, then the next realize(Velocity) call will use them to calculate
the prescribed qdots=N*u. realize(Dynamics) calculates known forces and the 
prescribed udoti=udoti(t,q,u). realize(Acceleration) calculates the
remaining udots, lambdas, taus, and all the zdots. 

Note that this method is \e not used to set prescribed udots, because those are
not state variables. Instead, prescribed udots (which depend on time, positions,
and velocities) are set as part of realize(Dynamics). 

@return \c true if any change was made to \a state. 
@see prescribe(), prescribeU() **/
bool prescribeQ(State& state) const;

/** This solver sets prescribed u's to their known values u=u(t,q) as a 
function of time and position variables q (and possibly earlier-stage discrete 
state variables).

We expect the supplied State already to have been realized to Stage::Position. 
Note that the \e derivatives of prescribed variables (which are
of necessity also prescribed but are not themselves state variables) are set in
the subsequent realize() calls. For example, prescribeU() sets the 
prescribed u's, then the next realize(Velocity) call will use them to calculate
the prescribed qdots=N*u. realize(Dynamics) calculates known forces and the 
prescribed udoti=udoti(t,q,u). realize(Acceleration) calculates the
remaining udots, lambdas, taus, and all the zdots. 

@return \c true if any change was made to \a state.
@see prescribe(), prescribeQ() **/
bool prescribeU(State& state) const;

/** Return the indices of the q's in the given \a state that are free; that is,
they are not affected by prescribeQ(). The indices are copied into the return
argument \a freeQs in ascending order; missing ones are prescribed. \a freeQs is 
resized if necessary to the number of free q's, nfq. \a state must be realized
through Stage::Instance. **/
void getFreeQIndex(const State& state, Array_<SystemQIndex>& freeQs) const;

/** Return the indices of the u's in the given \a state that are free; that is,
they are not affected by prescribeU(). The indices are copied into the return
argument \a freeUs in ascending order; missing ones are prescribed. \a freeUs is
resized if necessary to the number of free u's, nfu. \a state must be realized
through Stage::Instance. **/
void getFreeUIndex(const State& state, Array_<SystemUIndex>& freeUs) const;
/**@}**/


//------------------------------------------------------------------------------
/**@name              Event handlers and reporters

These methods provide a simplified interface to the Simbody event system. Each
instance creates one Event, one EventTrigger, and one EventAction, combined
into an EventHandler or EventReporter. The methods below are used to attach 
EventHandler and EventReporter objects to this %System. These are actually 
added to the SystemGlobalSubsystem that is contained in every System. **/
/**@{**/
/** Add a ScheduledEventHandler to this %System, which takes over ownership
of the event handler object. **/
void adoptEventHandler(ScheduledEventHandler* handler);
/** Add a TriggeredEventHandler to this %System, which takes over ownership
of the event handler object. **/
void adoptEventHandler(TriggeredEventHandler* handler);
/** Add a ScheduledEventReporter to this %System, which takes over ownership
of the event reporter object. **/
void adoptEventReporter(ScheduledEventReporter* reporter);
/** Add a TriggeredEventReporter to this %System, which takes over ownership
of the event reporter object. **/
void adoptEventReporter(TriggeredEventReporter* reporter);
/**@}**/


//------------------------------------------------------------------------------
/**@name                    The Discrete System

These methods deal with the discrete (event-driven) aspects of this %System. 
These methods are primarily for the use of Integrator objects that detect
events and TimeStepper objects that deal with them by calling appropriate
handlers. In general, events interrupt the continuous evolution of a system,
dividing a simulation into continuous segments interrupted by discontinuous
changes of state. There are also discrete variables that can be updated without
disrupting the continuous flow. These are called "auto-update variables" and are
updated at the beginning of every integration step, right *after* the initial 
state derivatives have been recorded. Thus the step proceeds with those 
derivatives even if the update could cause them to change. These updates are
initiated by Integrator objects at the start of every continuous step, while
trajectory-affecting discrete variable updates are initiated only by TimeStepper
objects via calls to event handlers.

@see State::allocateAutoUpdateDiscreteVariable()
**/
/**@{**/

/** Add a new Event to this %System and obtain a unique, small-integer EventId 
for it. The %System takes over ownership of the heap-allocated Event object. **/
EventId adoptEvent(Event* eventp);

/** Return the number n of Event objects contained in this %System. The
corresponding EventIds are EventId(0) through EventId(n-1). Don't 
confuse this with the number of event occurrences at runtime. **/
int getNumEvents() const;

/** Return a const reference to an Event specified by EventId. Throws an 
exception if the EventId is invalid or if there is no Event corresponding to
it in this %System. **/
const Event& getEvent(EventId id) const;

/** Return a writable reference to an Event specified by EventId. Throws an 
exception if the EventId is invalid or if there is no Event corresponding to
it in this %System. **/
Event& updEvent(EventId id);

/** Check whether this %System has an Event with the given EventId. **/
bool hasEvent(EventId id) const;

/** Return a const reference to the built-in Initialization Event. This Event
is triggered once at the start of a simulation. **/
const Event::Initialization& getInitializationEvent() const;
/** Return a writable reference to the built-in Initialization Event. This Event
is triggered once at the start of a simulation. **/
Event::Initialization& updInitializationEvent();

/** Return a const reference to the built-in TimeAdvanced Event. This Event is
triggered whenever an integrator successfully completes a step. **/
const Event::TimeAdvanced& getTimeAdvancedEvent() const;
/** Return a const reference to the built-in TimeAdvanced Event. This Event is
triggered whenever an integrator successfully completes a step. **/
Event::TimeAdvanced& updTimeAdvancedEvent();

/** Return a const reference to the built-in Termination Event. This Event is
triggered once at the end of a simulation. **/
const Event::Termination& getTerminationEvent() const;
/** Return a writable reference to the built-in Termination Event. This Event is
triggered once at the end of a simulation. **/
Event::Termination& updTerminationEvent();

/** Return a const reference to the built-in ExtremeValueIsolated Event. This
Event is triggered whenever a watched value reaches and extremum of 
interest (could be a maximum, minimum, or either). **/
const Event::ExtremeValueIsolated& getExtremeValueIsolatedEvent() const;
/** Return a writable reference to the built-in ExtremeValueIsolated Event. This
Event is triggered whenever a watched value reaches and extremum of 
interest (could be a maximum, minimum, or either). **/
Event::ExtremeValueIsolated& updExtremeValueIsolatedEvent();

/** Add a new EventTrigger to this %System and obtain a unique, small-integer 
EventTriggerId for it. The %System takes over ownership of the heap-allocated 
EventTrigger object. **/
EventTriggerId adoptEventTrigger(EventTrigger* triggerp);

/** Return the number n of EventTrigger objects contained in this %System.
The corresponding EventTriggerIds are EventTriggerId(0) through 
EventTriggerId(n-1). **/
int getNumEventTriggers() const;

/** Return a const reference to an EventTrigger specified by EventTriggerId. 
Throws an exception if the EventTriggerId is invalid or if there is no 
EventTrigger corresponding to it in this %System. **/
const EventTrigger& getEventTrigger(EventTriggerId id) const;

/** Return a writable reference to an EventTrigger specified by EventTriggerId. 
Throws an exception if the EventTriggerId is invalid or if there is no 
EventTrigger corresponding to it in this %System. **/
EventTrigger& updEventTrigger(EventTriggerId id);

/** Check whether this %System has an EventTrigger with the given 
EventTriggerId. **/
bool hasEventTrigger(EventTriggerId id) const;

/** Return a reference to the built-in EventTrigger that accompanies
initialization of a time stepping study (the Initialization event). **/
const InitializationTrigger& getInitializationTrigger() const;
/** Return a reference to the built-in EventTrigger that accompanies
the TimeAdvanced event. **/
const TimeAdvancedTrigger& getTimeAdvancedTrigger() const;
/** Return a reference to the built-in EventTrigger that accompanies
termination of a time stepping study (the Termination event). **/
const TerminationTrigger& getTerminationTrigger() const;

/** Return a list of event witnesses currently present for this %System in the
given Study's internal state. These are the witnesses that must be 
watched during the next time step to determine if an event has occurred. **/
void findActiveEventWitnesses
   (const Study&                            study, 
    Array_<const EventTrigger::Witness*, 
           ActiveWitnessIndex>&             witnesses) const; 

/** Determine if we are approaching a significant zero crossing of some
event witness and if so provide an estimate of when it will occur. We only
consider the interval between `study`'s internal time and the given
`timeLimit`; zero crossings outside that interval don't count. We report the 
earliest zero crossing we anticipate, and the required time localization 
window within which it must be isolated. If several witnesses have this zero
crossing, then the window is the narrowest of those. We return a list of the
witnesses we think will trigger. If there are no anticipated zeroes in the
interval, the next time and window will be Infinity and the list of witnesses
will be empty. **/
void predictNextWitnessTriggerTime
   (const Study&                            study,
    double                                  timeLimit,
    double&                                 timeOfNextZeroCrossing,
    double&                                 timeWindow,
    Array_<const EventTrigger::Witness*>&   witnesses) const; 

    

/** Return a list of event timers currently present for this %System in the
given study's internal state. These are the timers that must be consulted
to see whether the upcoming time step should be shortened because of a scheduled
report or change event. **/
void findActiveEventTimers
   (const Study&                            study, 
    Array_<const EventTrigger::Timer*,
           ActiveTimerIndex>&               timers) const;


/** Consult the Event Timers to find the next time at which a scheduled report 
Event occurs, and the next time for a scheduled change Event. Return the 
Timer(s) that trigger each of those. We are given the time of last report and 
time of last change and won't return Triggers that are at exactly those times 
(to avoid duplicates). Note that all the Timers in one of the returned lists had
simultaneous next trigger times (to within timer noise). **/
void findNextScheduledEventTimes
   (const Study&        study,
    double              timeOfLastReport,
    double              timeOfLastChange,
    double&             timeOfNextReport,
    EventTriggers&      reportTimers,
    double&             timeOfNextChange,
    EventTriggers&      changeTimers) const;

/** Given a list of EventTrigger objects that have occurred simultaneously,
determine the unique set of Event objects whose actions should be performed. As
a side effect, the mutable occurrence counter for each trigger and each unique 
event is incremented by one.

@param[in]          triggers    
    List of triggers that occurred simultaneously.
@param[in,out]      appendTriggeredEvents
    List of (event,causes) pairs where "causes" is a subset of `triggers`. The
    list will be *appended* to; clear it first if you don't want it to grow.
@param[in,out]      appendIgnoredEvents
    List of unrecognized EventId values that were in `triggers`. The list will
    be *appended* to; clear it first if you don't want it to grow.

Note that when multiple triggers cause the same event, the Event actions still
must be performed only once. We map the EventId values stored in the triggers
to Event objects, and append a list pairing each Event with the triggers that
caused it. If there is no Event corresponding to an EventId, that EventId is
ignored except to append it (once) to the `appendIgnoredEvents` list, which
may be useful for issuing diagnostics or for debugging.

The order of the given triggers is maintained in the output arrays; all events
caused by the first trigger come first, then all events from the second trigger,
and so on, except that duplicate events will not appear.

Be sure to clear the output arrays before calling this method if you want them
to start fresh; we're just appending here. Typically these will be very short 
lists -- most commonly there will be one trigger, one event, and no unrecognized
events. **/
void noteEventOccurrence
   (const EventTriggers&    triggers,
    EventsAndCauses&        appendTriggeredEvents,
    Array_<EventId>&        appendIgnoredEvents) const;
    

/** Execute any report() actions associated with Events whose occurrence has
been triggered. The supplied Study should have invoked noteEventOccurrence()
and the `triggeredEvents` parameter here should be the output from that method.
The state being reported is that returned by `study.getCurrentState()` and may 
be an interpolated state. 
@see noteEventOccurrence(), performEventChangeActions() **/
void performEventReportActions
   (const Study&            study,
    const EventsAndCauses&  triggeredEvents) const;

/** Execute any change() actions associated with Events whose occurrence has
been triggered. The supplied Study should have invoked noteEventOccurrence()
and the `triggeredEvents` parameter here should be the output from that method.
The state that triggered the Events is obtained from `study.updInternalState()`
and may be modified here. 
@see noteEventOccurrence(), performEventReportActions() **/
void performEventChangeActions
   (Study&                  study,
    const EventsAndCauses&  triggeredEvents,
    EventChangeResult&      result) const;

/**@}**/


//------------------------------------------------------------------------------
/**@name                      The Fast System

(NOT IMPLEMENTED YET) These methods are for state variables whose dynamics
are considered very fast compared to the rest of the system, meaning that
their values can be determined by solving for algebraic equilibrium 
conditions rather than integrating differential equations. **/
/**@{**/

/** This optional method should modify fast variables at the given stage 
until they satisfy some relaxation criteria. The criteria may involve 
anything in the State *except* fast variables at higher stages. Anything
that can be calculated from the state is also fair game provided that 
those calculations do not depend on higher-stage fast variables. 
"Relaxation" criteria may require that fast variables satisfy some 
implicit or explicit algebraic conditions (constraints), or reach some 
minimization or maximization condition. A common criterion is that fast
q's are moved to minimize potential energy; that can be achieved by 
driving the fast mobilities' lambdas (calculated joint torques) to zero
since they are the potential energy gradient. That may require repeated 
realization to Acceleration stage.

Note that when q's are fast, the corresponding u's and udots are 
\e prescribed (to zero), they are not \e fast. And when u's are fast, their 
udots are also zero (their q's are regular integrated variables). Fast 
z's have zero zdots. Any other variables (that is, x-partition variables) 
can also be fast but don't have derivatives.

TODO: should take options and return results. **/
void relax(State& state, Stage stage, Real accuracy=-1) const;
/**@}**/


//------------------------------------------------------------------------------
/**@name             Kinematic differential equations

These methods are primarily for use by numerical integration methods.

Generalized coordinates q are not independent of generalized speeds u;
they are related by the kinematic differential equation qdot=N(q)*u, where
N is an nqXnu block diagonal, invertible matrix in the sense that 
u=pinv(N)*qdot, where pinv(N) is the pseudoinverse of N. N has full column rank,
so pinv(N)*N = I, but N*pinv(N) != I.

Just as N provides the relation between velocities expressed in u-space and
the equivalent in q-space, its transpose ~N relates forces in q-space to
their equivalent in u-space: fu=~N*fq, and fq=~pinv(N)*fu. (Note that
~pinv(X)==pinv(~X).) This satisfies power=~fq*qdot==~fu*u as it must.

We provide fast O(n) operators for multiplication by N, pinv(N), and their 
transposes. N is often mostly an identity matrix so very little computation
is required. **/
/**@{**/
/** Calculate dq=N*u in O(n) time (very fast). **/
void multiplyByN(const State& state, const Vector& u, 
                 Vector& dq) const;
/** Calculate fu=~N*fq in O(n) time (very fast). **/
void multiplyByNTranspose(const State& state, const Vector& fq, 
                          Vector& fu) const;
/** Calculate u=pinv(N)*dq in O(n) time (very fast). **/
void multiplyByNPInv(const State& state, const Vector& dq, 
                     Vector& u) const;
/** Calculate fq=~pinv(N)*fu in O(n) time (very fast). **/
void multiplyByNPInvTranspose(const State& state, const Vector& fu, 
                              Vector& fq) const;
/**@}**/


//------------------------------------------------------------------------------
/**@name                         Statistics

The %System keeps mutable statistics internally, initialized to zero at 
construction. These <em>must not</em> affect results in any way. **/
/**@{**/

/** Zero out the statistics for this System. Although the statistics are 
mutable, we only allow them to be reset by a caller who has write access since
setting the stats to zero is associated with construction. **/
void resetAllCountersToZero();

    // Realization

/** Whenever the system was realized from Stage-1 to the indicated Stage,
this counter is bumped. Note that a single call to realize() can cause 
several counters to get bumped. **/
int getNumRealizationsOfThisStage(Stage) const;

/** Return the total number of calls to realizeTopology(), realizeModel(),
or realize(), regardless of whether these routines actually did
anything when called. **/
int getNumRealizeCalls() const;

    // Prescribed motion

/** Return the total number of calls to the System's prescribeQ() method. **/
int getNumPrescribeQCalls() const;
/** Return the total number of calls to the System's prescribeU() method. **/
int getNumPrescribeUCalls() const;

    // Projection

/** Return the total number of calls to projectQ(), regardless of
whether the call did anything. **/
int getNumProjectQCalls() const;
/** Return the total number of calls to projectQ() that failed. **/
int getNumFailedProjectQCalls() const;
/** How many of the successful projectQ() calls actually did a constraint 
projection, rather than returning quickly? **/
int getNumQProjections() const;
/** How many of the projectQ() calls that did a constraint projection also
projected an error estimate? **/
int getNumQErrorEstimateProjections() const;

/** Return the total number of calls to projectU(), regardless of
whether the call did anything. **/
int getNumProjectUCalls() const;
/** Return the total number of calls to projectU() that failed. **/
int getNumFailedProjectUCalls() const;
/** How many of the successful projectU() calls actually did a constraint 
projection, rather than returning quickly? **/
int getNumUProjections() const;
/** How many of the projectU() calls that did a constraint projection also
projected an error estimate? **/
int getNumUErrorEstimateProjections() const;

/**@}**/


//------------------------------------------------------------------------------
/**@name                Construction and bookkeeping

These methods are mostly for use by concrete Systems and will not typically
be of interest to users. System does not currently support copy construction
or assignment. **/
/**@{**/
// No default constructor; must supply a System::Guts object.
System() = delete;
System(const System&) = delete;
System(const System&&) = delete;
System& operator=(const System&) = delete;

/** Destructor here will invoke the virtual destructor in the System::Guts
class and thus clean up all unneeded detritus owned by this %System. **/
~System();

/** Return the name assigned to this %System (arbitrary). **/
const String& getName()    const;
/** Return the version string assigned to this %System (arbitrary). **/
const String& getVersion() const;

/** Take over ownership of the supplied subsystem and install it into 
the next free subsystem slot. The new slot index is returned. **/
SubsystemIndex adoptSubsystem(Subsystem& child);

/** How may Subsystems are in here? **/
int getNumSubsystems() const;
/** Obtain read-only access to a particular subsystem by its index. **/
const Subsystem& getSubsystem(SubsystemIndex)   const;
/** Obtain writable access to a particular subsystem by its index. **/
Subsystem& updSubsystem(SubsystemIndex);

/** Get read-only access to the default subsystem which is present in every 
system. **/
const SystemGlobalSubsystem& getSystemGlobalSubsystem() const;
/** Get writable access to the default subsystem which is present in every 
system. **/
SystemGlobalSubsystem& updSystemGlobalSubsystem();

/** Implicitly convert this System into a const Subsystem reference; this 
actually returns a reference to the SystemGlobalSubsystem contained in this 
System. **/
operator const Subsystem&() const;
/** Implicitly convert this System into a writable Subsystem reference; this 
actually returns a reference to the SystemGlobalSubsystem contained in this 
System. **/
operator Subsystem&();

/** There can be multiple handles referring to the same System::Guts object; 
they are considered to be the same %System. Two empty handles are not 
considered to be the same %System. **/
bool isSameSystem(const System& otherSystem) const;

/**@}**/

//------------------------------------------------------------------------------
/**@name              Advanced/obscure/debugging/obsolete

You probably don't want to use these methods. **/
/**@{**/

/** (Advanced) You can check whether realizeTopology() has been called since the
last topological change to this System. If you don't check and just plunge
ahead you are likely to encounter an exception since very few things will work 
without topology having been realized. **/
bool systemTopologyHasBeenRealized() const;

/** (Advanced) Return the current version number of this system's Topology
cache information. This is a counter that is incremented each time the Topology
is invalidated. Any State to be used with this System must have the same
Topology version number as the System does. The version number is returned
regardless of whether topology has been realized; you can check that with
systemTopologyHasBeenRealized(). 
@see State::getSystemTopologyStageVersion() **/
StageVersion getSystemTopologyCacheVersion() const;

/** (Really advanced) Set the current version number of this System's 
Topology cache information. Don't use this method unless you really know what
you're doing! This has no effect on realization status; if Topology has not
yet been realized this is the version number it will have after the next 
realizeTopology() call. **/
void setSystemTopologyCacheVersion(StageVersion topoVersion) const;

/** (Advanced) Mark the Topology stage of this system and all its subsystems
"not realized." This is normally handled automatically whenever you make a 
Topology-stage change to any Subsystem. Occasionally you may want to force 
recomputation of the Topology-stage cache, for example during testing. After 
this call the method systemTopologyHasBeenRealized() will return false and you 
will not be able to call getDefaultState(). A subsequent call to
realizeTopology() will invoke all the subsystems' realizeTopology() methods.
The Topology stage version number will have changed, so all previously-created
State objects will be invalid. **/
void invalidateSystemTopologyCache() const;

/** (Advanced) Generate all decorative geometry computable at a specific stage;
this is useful for visualizers. This will throw an exception if the state hasn't
already been realized to the given stage. Note that the list is not inclusive --
you have to request geometry from each stage to get all of it. This routine 
asks each Subsystem in succession to generate its decorative geometry and
append it to the end of the array. If the stage is Stage::Topology, 
realizeTopology() must already have been called but the State is ignored. **/
void calcDecorativeGeometryAndAppend(const State& state, Stage stage, 
                                     Array_<DecorativeGeometry>& geom) const;

/** Obtain a const reference to the System::Guts object to which this %System
handle refers. You should then `dynamic_cast` the returned reference to a 
reference to your concrete Guts class. **/
const Guts& getSystemGuts() const {assert(m_guts); return *m_guts;}

/** Obtain a writable reference to the System::Guts object to which this %System
handle refers. You should then `dynamic_cast` the returned reference to a 
reference to your concrete Guts class. **/
Guts& updSystemGuts() {assert(m_guts); return *m_guts;}

/** Return true if this %System handle is not empty. **/
bool hasGuts() const {return m_guts != nullptr;}

/** Internal use only. **/
bool isOwnerHandle() const;
/** Internal use only. **/
bool isEmptyHandle() const {return m_guts == nullptr;}

/** <b>(Deprecated)</b> Use getSystemGlobalSubsystem() instead. The name changed
in Simbody 4.0. **/
DEPRECATED_14("use getSystemGlobalSubsystem() instead")
const SystemGlobalSubsystem& getDefaultSubsystem() const
{   return getSystemGlobalSubsystem(); }
/** <b>(Deprecated)</b> Use updSystemGlobalSubsystem() instead. The name changed
in Simbody 4.0. **/
DEPRECATED_14("use updSystemGlobalSubsystem() instead")
SystemGlobalSubsystem& updDefaultSubsystem()
{   return updSystemGlobalSubsystem(); }

/** <b>(Deprecated)</b> Use `adoptEventHandler()` instead. The name changed 
in Simbody 4.0. **/
DEPRECATED_14("use adoptEventHandler() instead")
void addEventHandler(ScheduledEventHandler* handler)
{   adoptEventHandler(handler); }
/** <b>(Deprecated)</b> Use `adoptEventHandler()` instead. The name changed 
in Simbody 4.0. **/
DEPRECATED_14("use adoptEventHandler() instead")
void addEventHandler(TriggeredEventHandler* handler)
{   adoptEventHandler(handler); }
/** <b>(Deprecated)</b> Use `adoptEventReporter()` instead. The name changed 
in Simbody 4.0. **/
DEPRECATED_14("use adoptEventReporter() instead")
void addEventReporter(ScheduledEventReporter* reporter)
{   adoptEventReporter(reporter); }
/** <b>(Deprecated)</b> Use `adoptEventReporter()` instead. The name changed 
in Simbody 4.0. **/
DEPRECATED_14("use adoptEventReporter() instead")
void addEventReporter(TriggeredEventReporter* reporter)
{   adoptEventReporter(reporter); }
/**@}**/

protected:
/** Put new *unowned* System::Guts object into a new %System handle which takes
over ownership. If Guts is already owned by some other %System handle this
routine will throw an exception. Also, Guts must not yet contain any 
subsystems because we need to put the SystemGlobalSubsystem in the 0th slot. 
Any resources that are required for services provided by the base %System
class will be allocated here and will thus be present when you start putting
together your concrete %System. **/
explicit System(System::Guts* guts);

private:
// This is the only data member in this class. Also, any class derived from
// System must have *NO* data members at all (data goes in the Guts class).
Guts*   m_guts;
};


//==============================================================================
//                     PROJECT OPTIONS and PROJECT RESULTS
//==============================================================================
/** Options for the advanced project() methods. The default is to require
project() methods to reduce constraint errors to an RMS norm of 1e-4, while
asking them to attempt 10X tighter accuracy if possible. **/
class ProjectOptions {
public:
    enum Option {
        /** Take all defaults. **/
        None            = 0x0000,
        /** This option says we expect the state to be close to a solution 
        already and restricts projection to move downhill in the local 
        vicinity. This should be used during any continuous integration to 
        prevent erroneous jumps in the state. **/
        LocalOnly       = 0x0001,
        /** Normally failure to meet the accuracy requirements throws an
        exception. This will force the project() method to quietly return bad 
        status instead. **/
        DontThrow       = 0x0002,
        /** Use the stricter infinity (max absolute value) norm rather than
        the default RMS norm to determine when accuracy has been achieved. **/
        UseInfinityNorm = 0x0004,
        /** Normally a project() method will return immediately after 
        evaluating the norm if it is already at or below the required accuracy.
        This option forces it to make at least one iteration. **/
        ForceProjection = 0x0008,
        /** A project() method is free to use an out-of-date Jacobian when
        solving the nonlinear system. This option forces recalculation of
        the Jacobian at the start of each iteration. **/
        ForceFullNewton = 0x0010
    };

    /** Default constructor sets options to their default values. **/
    ProjectOptions() {clear();}
    /** This constructor allows the default accuracy to be overridden while
    leaving all other options at their default values. If \a accuracy is 
    nonpositive, the default accuracy will be used. **/
    explicit ProjectOptions(Real accuracy) 
    {   clear(); setRequiredAccuracy(accuracy); }
    /** This constructor creates default options except one setting one
    non-default Option. You can add more with the setOption() method. **/
    explicit ProjectOptions(Option opt)
    {   clear(); setOption(opt); }

    /** Restore this object to its default-constructed state (no options
    selected, default accuracy and overshoot). A reference to the
    newly-cleared object is returned. **/
    ProjectOptions& clear() 
    {   optionSet=0; setAccuracyDefaults(); return *this; }

    /** The norm of the constraint errors must be driven to below this value
    for a project() to be considered successful. Normally an RMS norm is used
    but you can override that to use an infinity norm instead. If accuracy
    is nonpositive we'll use the default; see getDefaultRequiredAccuracy(). **/
    ProjectOptions& setRequiredAccuracy(Real accuracy) {
        requiredAccuracy = accuracy > 0 ? accuracy 
                                        : getDefaultRequiredAccuracy();
        return *this;
    }

    /** Project will attempt to reach accuracy*overshoot but settle for 
    just accuracy. **/ 
    ProjectOptions& setOvershootFactor(Real overshoot) {
        assert(0 < overshoot && overshoot <= 1);
        desiredOvershoot = overshoot;
        return *this;
    }

    /** Project will fail immediately if the initial norm is greater than
    the projection limit, with status FailureToConverge. **/ 
    ProjectOptions& setProjectionLimit(Real limit) {
        assert(limit > 0);
        projectionLimit = limit;
        return *this;
    }

    /** Remove a given option from the set. Nothing happens if the option wasn't
    already set. **/
    ProjectOptions& clearOption(Option opt) 
    {   optionSet &= ~(unsigned)opt; return *this; }
    /** Set a particular option. Multiple calls are or'ed together. **/
    ProjectOptions& setOption  (Option opt) 
    {   optionSet |= (unsigned)opt; return *this; }

    /** Return the current value for the required accuracy option. **/
    Real getRequiredAccuracy()       const {return requiredAccuracy;}
    /** Return the factor by which a project() method should try to do better
    than the required accuracy. **/
    Real getOvershootFactor() const {return desiredOvershoot;}
    /** Return the maximum norm we're allowed to attempt to correct. **/
    Real getProjectionLimit() const {return projectionLimit;}

    bool isOptionSet(Option opt) const {return (optionSet&(unsigned)opt) != 0;}

    static Real getDefaultRequiredAccuracy() {return Real(1e-4);}
    static Real getDefaultOvershootFactor()  {return Real(0.1);} //i.e., 1e-5

    // Set operators: not, or, and, set difference
    ProjectOptions& operator|=(const ProjectOptions& opts) 
    {   optionSet |= opts.optionSet; return *this; }
    ProjectOptions& operator&=(const ProjectOptions& opts) 
    {   optionSet &= opts.optionSet; return *this; }
    ProjectOptions& operator-=(const ProjectOptions& opts) 
    {   optionSet &= ~opts.optionSet; return *this; }

    ProjectOptions& operator|=(Option opt) {setOption(opt); return *this;}
    ProjectOptions& operator-=(Option opt) {clearOption(opt); return *this;}

private:
    Real     requiredAccuracy;
    Real     desiredOvershoot; // try for accuracy*overshoot
    Real     projectionLimit;  // abort if initial norm is worse than this
    unsigned optionSet;

    void setAccuracyDefaults() {
        requiredAccuracy = getDefaultRequiredAccuracy();
        desiredOvershoot = getDefaultOvershootFactor(); 
        projectionLimit  = Infinity; // we'll try from however far away
    }
};

/** Results for advanced users of project() methods. **/
class ProjectResults {
public:
    ProjectResults() {clear();}

    enum Status {
        /** This object has not been filled in yet and holds no results. **/
        Invalid                 = -1,
        /** The project() was successful either because no projection was
        necessary or projection was able to achieve the required accuracy. **/
        Succeeded               = 0,
        /** Projection converged but was unable to achieve the required
        accuracy. **/
        FailedToAchieveAccuracy = 1,
        /** The Newton iterations were diverging. This is especially common
        when the LocalOnly option is set since project() will quit at the 
        first sign of divergence in that case. This is also the return if
        a projection limit has been set and the initial norm is larger. **/
        FailedToConverge        = 2    
    };

    /** Restore this object to its default-constructed state, with the return
    status set to Invalid. **/
    ProjectResults& clear() {
        m_exitStatus = Invalid;
        m_anyChangeMade = m_projectionLimitExceeded = false;
        m_numIterations = 0;
        m_worstError = -1;
        m_normOnEntrance = m_normOnExit = NaN;
        return *this;
    }
    bool    isValid()        const {return m_exitStatus != Invalid;}
    Status  getExitStatus()  const {return m_exitStatus;}

    bool getAnyChangeMade()  const {assert(isValid());return m_anyChangeMade;}
    int  getNumIterations()  const {assert(isValid());return m_numIterations;}
    Real getNormOnEntrance() const {assert(isValid());return m_normOnEntrance;}
    Real getNormOnExit()     const {assert(isValid());return m_normOnExit;}
    int  getWorstErrorOnEntrance()    const 
    {   assert(isValid());return m_worstError; }
    bool getProjectionLimitExceeded() const 
    {   assert(isValid());return m_projectionLimitExceeded; }

    ProjectResults& setExitStatus(Status status) 
    {   m_exitStatus=status; return *this; }
    ProjectResults& setAnyChangeMade(bool changeMade) 
    {   m_anyChangeMade=changeMade; return *this; }
    ProjectResults& setProjectionLimitExceeded(bool limitExceeded) 
    {   m_projectionLimitExceeded=limitExceeded; return *this; }
    ProjectResults& setNumIterations(int numIterations) 
    {   m_numIterations=numIterations; return *this; }
    ProjectResults& setNormOnEntrance(Real norm, int worstError) 
    {   m_normOnEntrance=norm; m_worstError=worstError; return *this; }
    ProjectResults& setNormOnExit(Real norm) 
    {   m_normOnExit=norm; return *this; }
private:
    Status  m_exitStatus;
    bool    m_anyChangeMade;
    bool    m_projectionLimitExceeded;
    int     m_numIterations;
    int     m_worstError;       // index of worst error on entrance
    Real    m_normOnEntrance;   // in selected rms or infinity norm
    Real    m_normOnExit;
};



//==============================================================================
//                    REALIZE OPTIONS and REALIZE RESULTS 
//==============================================================================
/** (NOT USED YET) Options for the advanced realize() methods. **/
class RealizeOptions {
    unsigned int optionSet;
    explicit RealizeOptions(unsigned o) : optionSet(o) { }
public:

    enum Option {
        None      = 0x00,
        DontThrow = 0x01
    };


    RealizeOptions() : optionSet(0) { }

    // This is an implicit conversion
    RealizeOptions(Option opt) : optionSet((unsigned)opt) { }

    // Implicit conversion to bool when needed
    operator bool() const {return optionSet != 0;}
    bool isEmpty() const {return optionSet==0;}

    bool isOptionSet(Option opt) const {return (optionSet&(unsigned)opt) != 0;}
    void clear() {optionSet=0;}
    void clearOption(Option opt) {optionSet &= ~(unsigned)opt;}
    void setOption  (Option opt) {optionSet |= (unsigned)opt;}

    // Set operators: or, and
    RealizeOptions& operator|=(RealizeOptions opts) {optionSet |= opts.optionSet; return *this;}
    RealizeOptions& operator&=(RealizeOptions opts) {optionSet &= opts.optionSet; return *this;}

    RealizeOptions& operator|=(Option opt) {setOption(opt); return *this;}
    RealizeOptions& operator-=(Option opt) {clearOption(opt); return *this;}
};

/** (NOT USED YET) Results for advanced users of realize() methods. **/
class RealizeResults {
};




} // namespace SimTK

#endif // SimTK_SimTKCOMMON_SYSTEM_H_
