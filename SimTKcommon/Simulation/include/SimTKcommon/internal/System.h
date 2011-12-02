#ifndef SimTK_SimTKCOMMON_SYSTEM_H_
#define SimTK_SimTKCOMMON_SYSTEM_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Simbody: SimTKcommon                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-11 Stanford University and the Authors.        *
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
#include "SimTKcommon/internal/State.h"
#include "SimTKcommon/internal/Subsystem.h"

#include <cassert>

namespace SimTK {

class DecorativeGeometry;
class DefaultSystemSubsystem;
class ScheduledEventHandler;
class ScheduledEventReporter;
class TriggeredEventHandler;
class TriggeredEventReporter;

//==============================================================================
//                    REALIZE OPTIONS and REALIZE RESULTS 
//==============================================================================
/** Options for the advanced realize() methods. **/
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

/** Results for advanced users of realize() methods. **/
class RealizeResults {
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

    ProjectOptions() {clear();}
    explicit ProjectOptions(Real accuracy) 
    {   clear(); setRequiredAccuracy(accuracy); }
    explicit ProjectOptions(Option opt)
    {   clear(); setOption(opt); }

    /** Restore this object to its default-constructed state (no options
    selected, default accuracy and overshoot). A reference to the
    newly-cleared object is returned. **/
    ProjectOptions& clear() 
    {   optionSet=0; setAccuracyDefaults(); return *this; }

    /** The norm of the constraint errors must be driven to below this value
    for a project() to be considered successful. Normally an RMS norm is used
    but you can override that to use an infinity norm instead. **/
    ProjectOptions& setRequiredAccuracy(Real accuracy) {
        assert(accuracy > 0);
        requiredAccuracy = accuracy;
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
    /** Select a given option from the set. Nothing happens if the option wasn't
    already set. **/
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
    bool    isValid()           const {return m_exitStatus != Invalid;}
    Status  getExitStatus()     const {return m_exitStatus;}

    bool getAnyChangeMade()     const {assert(isValid());return m_anyChangeMade;}
    int  getNumIterations()     const {assert(isValid());return m_numIterations;}
    Real getNormOnEntrance()    const {assert(isValid());return m_normOnEntrance;}
    Real getNormOnExit()        const {assert(isValid());return m_normOnExit;}
    int  getWorstErrorOnEntrance()    const {assert(isValid());return m_worstError;}
    bool getProjectionLimitExceeded() const {assert(isValid());return m_projectionLimitExceeded;}

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
//                                 SYSTEM
//==============================================================================
/** This is the handle class that serves as the abstract parent of all System
handles.

A System serves as a mediator for a group of interacting Subsystems. All will 
share a single system State, and typically subsystems will need access to 
content in the state which is produced by other subsystems. 

A System provides a unique SubsystemIndex (a small positive integer) for each
of its subsystems, and the subsystems are constructed knowing their indices. 
The indices are used subsequently by the subsystems to find their own entries 
in the system state, and by each subsystem to refer to others within the same
system. Index 0 is reserved for use by the System itself, e.g. for 
system-global state variables.

Concrete Systems understand the kinds of subsystems they contain. For example, 
a MultibodySystem might contain a mechanical subsystem, some force subsystems, 
and a geometry subsystem. At each computation stage, a subsystem is realized 
in a single operation. That operation can refer to computations from 
already-realized subsystems, but cannot initiate computation in other 
subsystems. The System must know the proper order with which to realize the 
subsystems at each stage, and that ordering is likely to vary with stage. For 
example, at Position stage the mechanical positions must be realized before 
the configuration-dependent force elements. However, at Acceleration stage, 
the force elements must be realized before the mechanical accelerations can 
be calculated.

There are two distinct users of this class:
  - System Users: people who are making use of a concrete System (which will
    inherit methods from this class)
  - System Developers: people who are writing concrete System classes
Note that System Users include people who are writing Studies, Reporters, 
Modelers and so on as well as end users who are accessing the System directly.

Only methods intended for System Users and a few bookkeeping methods are in 
the main System class, which is a SimTK Handle class, meaning that it consists 
only of a single pointer, which points to a System::Guts class. The Guts class 
is abstract, and virtual methods to be implemented by System Developers in the 
concrete System are defined there, along with other utilities of use to the 
concrete System Developer but not to the end user. The Guts class is declared 
in a separate header file, and only people who are writing their own System
classes need look there. **/
class SimTK_SimTKCOMMON_EXPORT System {
public:
class Guts; // local; name is System::Guts
friend class Guts;
private:
// This is the only data member in this class. Also, any class derived from
// System must have *NO* data members at all (data goes in the Guts class).
Guts* guts;
public:

/** Add a ScheduledEventHandler to this System, which takes over ownership
of the event handler object. The handler is actually
added to the DefaultSystemSubsystem that is contained in this System. **/
inline void addEventHandler(ScheduledEventHandler* handler);
/** Add a TriggeredEventHandler to this System, which takes over ownership
of the event handler object. The handler is actually
added to the DefaultSystemSubsystem that is contained in this System. **/
inline void addEventHandler(TriggeredEventHandler* handler);
/** Add a ScheduledEventReporter to this System, which takes over ownership
of the event reporter object. The handler is actually
added to the DefaultSystemSubsystem that is contained in this System. **/
inline void addEventReporter(ScheduledEventReporter* handler) const;
/** Add a TriggeredEventReporter to this System, which takes over ownership
of the event reporter object. The handler is actually
added to the DefaultSystemSubsystem that is contained in this System. **/
inline void addEventReporter(TriggeredEventReporter* handler) const;

/** (Advanced) This is a hint used for some default behaviors, such as 
determining an initial step size for an integrator, or the default unit error 
for a constraint error derivative from the original constraint. Most users can
ignore this and just take the default.

This should be set to roughly the time scale at which you expect to see 
interesting things happen, that is the scale at which you might choose to
view reporting output. An orbital simulation using seconds as time units might
set this to 10 or 100s, for example, while a biomechanical simulation could
use 0.1s. This will affect the time scale on which velocity constraints are
stabilized, with longer time scales being more demanding since there is more time to
drift. By default this is 0.1 time units, so 100ms for systems measuring time in 
seconds and 100fs for systems measuring time in ps. **/
System& setDefaultTimeScale(Real tc);
/** Get the currently-set value for the default time scale, 0.1 time units
if it has never been set. **/
Real getDefaultTimeScale() const;
/** (Advanced) This is a hint that can be used to get a sense of what a "unit 
length" looks like for this System in the units it uses. Most users can ignore
this and just take the default.

For example, a model of a small toy car expressed in MKS units might set this to 
0.01 (1 cm). The default for this is 1 length unit, meaning 1 meter in MKS and 1 
nm in MD units. **/
System& setDefaultLengthScale(Real lc);
/** Get the currently-set value for the default length scale, 1 length unit
if it has never been set. **/
Real getDefaultLengthScale() const;

/** This is a hint to visualization software as to which way this System's
designer considers to be "up".\ This is the best direction to use as the 
default up direction for the camera. The default up direction is  +YAxis, 
which is the same as the OpenGL convention for the camera up direction. You 
can set this to any of the coordinate axes in the positive or negative 
direction. For example, use setUpDirection(ZAxis) for the "virtual world" 
convention where ground is the x-y plane, or use setUpDirection(-ZAxis) for 
the aviation convention where +z points towards the ground. A visualizer that 
is showing a ground plane should make the ground plane normal be this up 
direction.
@see setUseUniformBackground() **/
System& setUpDirection(const CoordinateDirection& up);
/** Get the current setting of the "up" direction hint. **/
CoordinateDirection getUpDirection() const;

/** This is a hint to visualization software that this System is best viewed 
against a uniform background (e.g.\ all white) rather than against a ground 
plane.\ A molecular system will typically set this flag so that the visualizer
will not attempt to place the molecule on the ground. The default is to 
consider this system best viewed with a ground plane displayed, perpendicular 
to the "up" direction and located at a height of zero.
@see setUpDirection() **/
System& setUseUniformBackground(bool useUniformBackground);
/** Get the current setting of the "use uniform background" visualization
hint. **/
bool getUseUniformBackground() const;




    /////////////////
    // REALIZATION //
    /////////////////

/**@name                         Realization
These methods provide the ability to compute the consequences that follow
from the current values of state variables, leaving the results in the
appropriate cache. **/
/**@{**/

/** The following call must be made after any topological change has been 
made to this System, before the System can be used to perform any 
computations. Perhaps surprisingly, the method is const. That's because 
the topology cannot be changed by this method. Various mutable "cache" 
entries will get calculated, including the default State, a reference 
to which is returned. The returned State has a Topology stage version number
that matches the System, and will have already been realized 
through the Model Stage, using the defaults for the Model-stage 
variables, meaning that all later stage variables have been allocated 
and set to their default values as well. You can access this same 
default State again using getDefaultState(). If the current topology 
has already been realized, this call does nothing but return a reference 
to the already-built default State. **/
const State& realizeTopology() const;

/** This is available after realizeTopology(), and will throw an
exception if realizeTopology() has not been called since the
most recent topological change to this System. This method returns the
same reference returned by realizeTopology(). The State to which
a reference is returned was created by the most recent
realizeTopology() call. It has a Topology version number that matches the one
currently in this System, and has already been realized through the
Model Stage, using default values for all the Model-stage variables.
All later-stage variables have been allocated and set to their
default values. You can use this state directly to obtain information
about the System in its default state or you can use this state
to initialize other States to which you have write access. Those
States are then suitable for further computation with this System. **/
const State& getDefaultState() const;
State&       updDefaultState();


/** This call is required if Model-stage variables are changed from
their default values. The System topology must already have been
realized (that is, realizeTopology() must have been called since
the last topological change made to the System). Also, the supplied
State must already have been initialized to work with this System
either by copying the default state or some other State of this System.
If it has already been realized to Stage::Model or higher, nothing
happens here. Otherwise, all the state variables at Stage::Instance or
higher are allocated or reallocated (if necessary), and reinitialized
to their default values. NOTE: any State information at Stage::Instance
or higher in the passed-in State is *destroyed* here. The number, types
and memory locations of those state variables will change, so any
existing references or pointers to them are invalid after this call.
Note that this routine modifies its State argument, but makes no changes
at all to the System itself and is hence const. **/
void realizeModel(State&) const;

/** Realize the entire System to the indicated Stage. The passed-in
State must have been initialized to work with this System, and
it must already have been realized through Stage::Model, since
the realize() method doesn't have write access to the State.
If the state has already been realized to the requested stage
or higher, nothing happens here. Otherwise, the state is realized
one stage at a time until it reaches the requested stage. **/
void realize(const State& s, Stage g = Stage::HighestRuntime) const;

/** (Advanced) You can check whether realizeTopology() has been called since the
last topological change to this Syatem. If you don't check and just plunge
ahead you are likely to encounter an exception since very few things
will work without topology having been realized. **/
bool systemTopologyHasBeenRealized() const;

/** (Advanced) Return the current version number of this system's Topology
cache information. This is a counter that is incremented each time the Topology
is invalidated. Any State to be used with this System must have the same
Topology version number as the System does. The version number is returned
regardless of whether topology has been realized; you can check that with
systemTopologyHasBeenRealized(). 
@see State::getSystemTopologyStageVersion() **/
StageVersion getSystemTopologyCacheVersion() const;

/** (Really advanced) Set the current version number of this system's 
Topology cache information. Don't use this method unless you really know what
you're doing! This has no effect on realization status; if topology has not
yet been realized this is the version number it will have after the next 
realizeTopology() call. **/
void setSystemTopologyCacheVersion(StageVersion topoVersion) const;

/** (Advanced) Mark the Topology stage of this system and all its subsystems
"not realized." This is normally handled automatically by whenever you make a 
Topology-stage change to any subsystem. Occasionally you may want to force 
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
asks each subsystem in succession to generate its decorative geometry and
append it to the end of the array. If the stage is Stage::Topology, 
realizeTopology() must already have been called but the State is ignored. **/
void calcDecorativeGeometryAndAppend(const State&, Stage, 
                                     Array_<DecorativeGeometry>&) const;
/**@}**/

    ///////////////////////////
    // THE CONTINUOUS SYSTEM //
    ///////////////////////////

/**@name                    The Continuous System
These methods deal with the continuous (smoothly-evolving) aspects of this 
%System. They are primarly for use by numerical integrators. 
TODO: this is currenly just a vague design **/
/**@{**/

/** Free q's are the subset of the given state's q's that need to be 
integrated from corresponding qdots and optionally qdotdots. You can call this 
at Instance stage after we know which q's are prescribed. The result is in
a form suitable for constructing a view of a full q array that shows only
the free q's. **/
const Array_<int>& getFreeQIndex(const State& state) const;

/** Free u's are the subset of the given state's u's that need to be 
integrated from corresponding udots. You can call this at Instance stage after 
we know which u's are prescribed. The result is in a form suitable for 
constructing a view of a full u array that shows only the free u's. **/
const Array_<int>& getFreeUIndex(const State& state) const;

/** Free udots are the subset of the given state's udots that are determined
from forces. You can call this at Instance stage after we know which udots are 
prescribed. The result is in a form suitable for constructing a view of a full 
udot array that shows only the free udots. **/
const Array_<int>& getFreeUDotIndex(const State& state) const;

const Array_<bool>& getUseRelativeAccuracyU(const State& state) const;
Array_<bool>& updUseRelativeAccuracyU(State& state) const;

const Array_<bool>& getUseRelativeAccuracyZ(const State& state) const;
Array_<bool>& updUseRelativeAccuracyZ(State& state) const;


/** Given a qdot-like vector dq, weight it to produce dqw=Wq*dq, where
Wq=N*Wu*pinv(N). Wu is the set of u weights. **/
void calcWeightedDQ(const State& state, const Vector& dq, Vector& dqw) const;

/** Given a udot-like vector du, weight it to produce duw=Eu*du, where
Eu_i=min(Wu_i, 1/u_i) if we're looking for relative accuracy on u_i, otherwise
Eu_i=Wu_i. Wu is the set of u weights. **/
void calcWeightedDU(const State& state, const Vector& du, Vector& duw) const;

/** Given a zdot-like vector dz, weight it to produce dzw=Ez*dz, where
Ez_i=min(Wz_i, 1/z_i) if we're looking for relative accuracy on z_i, otherwise
Ez_i=Wz_i. Wz is the set of z weights. **/
void calcWeightedDZ(const State& state, const Vector& dz, Vector& dzw) const;

/**@name           Kinematic differential equations

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

        // UNCONSTRAINED



        // CONSTRAINED

/** This optional solver should set state variables q and u to known values
as a function of time and earlier-stage state variables.
    - prescribeQ sets each prescribed qi=qi(t).
    - prescribeU sets each prescribed ui=ui(t,q).

In each case we expect the supplied State already to have been realized to the 
previous stage. Note that the \e derivatives of prescribed variables (which are
of necessity also prescribed but are not themselves state variables) are set in
the subsequent realize() call. For example, prescribeU sets the 
prescribed u's, then the next realize(Velocity) call will use them to calculate
the prescribed qdots=N*u. realize(Dynamics) calculates known forces and the 
prescribed udoti=udoti(t,q,u). realize(Acceleration) calculates the
remaining udots, lambdas, taus, and all the zdots. 

Note that this method is \e not used to set prescribed udots, because those are
not state variables. Instead, prescribed udots (which depend on time, positions,
and velocities) are set as part of realize(Dynamics). 

@return \c true if any change was made to \a state **/
bool prescribeQ(State& state) const;
bool prescribeU(State& state) const;

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

@see prescribeQ(), prescribeU(), project(), realize() **/
void prescribe(State& state) const {
    realize(state, Stage::Time);
    prescribeQ(state);
    realize(state, Stage::Position);
    prescribeU(state);
}

/** This solver projects the given state back on to the position or velocity
constraint manifold, by the shortest path possible in the scaled norm defined in
that state, and modifying only free (non-prescribed) variables. Constraint 
errors are scaled by their unit error weightings, then satisfied to a given 
accuracy using an RMS norm, or optionally using the stricter infinity norm. 

Optionally, this method can also project out the constraint-normal component of
the passed-in error estimate vector yerrest.
This is part of the integration of the continuous DAE system and thus 
should never require an integrator restart. 

Options
- use infinity norm
- local projection only
- force projection
- don't throw an exception
- force full Newton

This method is not for satisfying acceleration constraints, which does not
involve modifications to state variables. Acceleration constraints are 
satisfied automatically when you realize a state to Acceleration stage using
realize(state); the resulting udots will satisfy the acceleration constraints
(if possible), even if the position and velocity constraints are not satisfied.

<h3>Theory</h3>
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
 **/


void project(State& state, Real accuracy=1e-4) const {
    ProjectResults projResults;
    ProjectOptions projOptions;
    projOptions.setRequiredAccuracy(accuracy);

    realize(state, Stage::Time);
    prescribeQ(state);
    realize(state, Stage::Position);
    projectQ(state, Vector(), projOptions, projResults);
    prescribeU(state);
    realize(state, Stage::Velocity);
    projectU(state, Vector(), projOptions, projResults);
}

/** Advanced: project free q's so that position constraints are satisfied and 
remove
the corresponding error from the supplied error estimate. This is primarily
intended for use by numerical integration algorithms. You must already have
set prescribed q's; this method will not modify them but may depend on their
current values. State must be realized to Position stage on entry and will 
still be realized through Position stage on return. 

If the norm of perr is already less than or equal to accuracy on entry, nothing
will happen unless you have selected the "ForceProjection" option. You can 
find out what actually happened by looking in the returned \a results.
@see ProjectOptions, ProjectResults, project()
**/
void projectQ(State&, Vector& qErrEst, 
             const ProjectOptions& options, ProjectResults& results) const;
void projectU(State&, Vector& uErrEst, 
             const ProjectOptions& options, ProjectResults& results) const;

        // FAST VARIABLES

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
void relax(State&, Stage, Real accuracy) const;
/**@}**/

    ////////////////////////////////
    // THE DISCRETE (SLOW) SYSTEM //
    ////////////////////////////////

/**@name                    The Discrete System
These methods deal with the discrete (event-driven) aspects of this %System. **/
/**@{**/

/** This determines whether this System wants to be notified whenever time
advances irreversibly. If set true, time advancement is treated as an
event, and the handleEvents() method is invoked with its \a cause argument
set to indicate a time-advanced event occurred. By default, time advancement 
proceeds silently. **/
void setHasTimeAdvancedEvents(bool); // default=false
/** Return the current value of the flag indicating whether this %System wants
an event generated whenever time advances irreversibly. **/
bool hasTimeAdvancedEvents() const;

/** This solver handles a set of events which a TimeStepper has denoted as 
having occurred at the given time and state. The event handler may make 
discontinuous changes in the State, in general both to discrete and continuous
variables, but \e not to time or topological information. If changes are made 
to continuous variables, the handler is required to make sure the returned 
state satisfies the constraints to the accuracy level specified in \a options.

On return, the handleEvents() method will set the output variable \a results
to indicate what happened. If any invoked handler determines that the 
occurrence of some event requires that the simulation be terminated, that 
information is returned in \a results and a well-behaved TimeStepper will stop
when it sees that.

Simbody will automatically set a field in \a results that says how much of
the \a state was changed by the handler so that the calling TimeStepper will
be able to determine how much reinitialization is required.

@see HandleEventsOptions, HandleEventsResults **/
void handleEvents(State&                        state, 
                  Event::Cause                  cause, 
                  const Array_<EventId>&        eventIds,
                  const HandleEventsOptions&    options,
                  HandleEventsResults&          results) const;
    
/** This method is similar to handleEvents(), but does not allow the State 
to be modified.  It is used for scheduled events that were marked as 
being reports. **/
void reportEvents(const State&                  state, 
                  Event::Cause                  cause, 
                  const Array_<EventId>&        eventIds) const;

/** This routine provides the Integrator with information it needs about the
individual event trigger functions, such as which sign transitions are
relevant and how tightly we need to localize. This is considered 
Instance stage information so cannot change during a continuous integration
interval (so an Integrator can process it upon restart(Instance)), 
however it can be updated whenever a discrete update is made to the 
State. A default implementation is provided which returns default 
EventTriggerInfo for each event trigger in \a state. The \a state must already be 
realized to Stage::Instance. **/
void calcEventTriggerInfo(const State&              state,
                          Array_<EventTriggerInfo>& triggerInfo) const;

/** This routine should be called to determine if and when there is an event
scheduled to occur at a particular time. This is a \e lot cheaper than
making the Integrator hunt these down like ordinary state-dependent events.
The returned time can be passed to the Integrator's stepping function as
the advance time limit. **/
void calcTimeOfNextScheduledEvent(const State&      state, 
                                  Real&             tNextEvent, 
                                  Array_<EventId>&  eventIds, 
                                  bool              includeCurrentTime) const;

/** This routine is similar to calcTimeOfNextScheduledEvent(), but is used for
"reporting events" which do not modify the state. Events returned by this
method should be handled by invoking reportEvents() instead of 
handleEvents(). **/
void calcTimeOfNextScheduledReport(const State&     state, 
                                   Real&            tNextEvent, 
                                   Array_<EventId>& eventIds, 
                                   bool             includeCurrentTime) const;
/**@}**/

    ////////////////
    // STATISTICS //
    ////////////////

/**@name                         Statistics
The System keeps mutable statistics internally, initialized to zero at 
construction. These *must not* affect results in any way. **/
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

    // Event handling and reporting

/** handleEvents() reports the lowest Stage it modified and we bump
the counter for that Stage. We also count reportEvents() calls here
as having "changed" Stage::Report. **/
int getNumHandlerCallsThatChangedStage(Stage) const;

/** This is the total number of calls to handleEvents() regardless
of the outcome. **/
int getNumHandleEventCalls() const;

/** This is the total number of calls to reportEvents() regardless
of the outcome. **/
int getNumReportEventCalls() const;
/**@}**/

/**@name                Construction and bookkeeping
These methods are mostly for use by concrete Systems and will not typically
be of interest to users. **/
/**@{**/
System() : guts(0) { }
System(const System&);
System& operator=(const System&);
~System();

const String& getName()    const;
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
const DefaultSystemSubsystem& getDefaultSubsystem() const;
/** Get writable access to the default subsystem which is present in every 
system. **/
DefaultSystemSubsystem& updDefaultSubsystem();

/** Implicitly convert this System into a const Subsystem reference; this 
actually returns a reference to the DefaultSystemSubsystem contained in this 
System. **/
inline operator const Subsystem&() const; // implemented below
/** Implicitly convert this System into a writable Subsystem reference; this 
actually returns a reference to the DefaultSystemSubsystem contained in this 
System. **/
inline operator Subsystem&();

// Internal use only
bool isOwnerHandle() const;
bool isEmptyHandle() const;

/** There can be multiple handles referring to the same System::Guts object; 
they are considered to be the same System. **/
bool isSameSystem(const System& otherSystem) const;


/** Obtain a const reference to the System::Guts object to which this handle
refers. You should then dynamic_cast the returned reference to a reference to 
your concrete Guts class. **/
const System::Guts& getSystemGuts() const {assert(guts); return *guts;}
/** Obtain a writable reference to the System::Guts object to which this handle
refers. You should then dynamic_cast the returned reference to a reference to 
your concrete Guts class. **/
System::Guts&       updSystemGuts()       {assert(guts); return *guts;}

/** Put new *unowned* Guts into this *empty* handle and take over ownership.
If this handle is already in use, or if Guts is already owned this
routine will throw an exception. **/
void adoptSystemGuts(System::Guts* g);

explicit System(System::Guts* g) : guts(g) { }
bool hasGuts() const {return guts!=0;}
/**@}**/
};


/** This is a concrete Subsystem that is part of every System.\ It provides a 
variety of services for the System, such as maintaining lists of event handlers
and reporters, and acting as a source of globally unique event IDs. 

To obtain the default subsystem for a System, call getDefaultSubsystem() or 
updDefaultSubsystem() on it. Also, a System can be implicitly converted
to a Subsystem, in which case it actually returns a reference to
this Subsystem. **/
class SimTK_SimTKCOMMON_EXPORT DefaultSystemSubsystem : public Subsystem {
public:
    explicit DefaultSystemSubsystem(System& sys);
    void addEventHandler(ScheduledEventHandler* handler);
    void addEventHandler(TriggeredEventHandler* handler);
    void addEventReporter(ScheduledEventReporter* handler) const;
    void addEventReporter(TriggeredEventReporter* handler) const;
    EventId createEventId(SubsystemIndex subsys, const State& state) const;
    void findSubsystemEventIds
       (SubsystemIndex subsys, const State& state, 
        const Array_<EventId>& allEvents, 
        Array_<EventId>& eventsForSubsystem) const;

    /** @cond **/  // don't let doxygen see this private class
    class Guts;
    /** @endcond **/
private:
    const Guts& getGuts() const;
    Guts& updGuts();
};

inline void System::addEventHandler(ScheduledEventHandler* handler)
{   updDefaultSubsystem().addEventHandler(handler); }
inline void System::addEventHandler(TriggeredEventHandler* handler)
{   updDefaultSubsystem().addEventHandler(handler); }
inline void System::addEventReporter(ScheduledEventReporter* handler) const
{   getDefaultSubsystem().addEventReporter(handler); }
inline void System::addEventReporter(TriggeredEventReporter* handler) const
{   getDefaultSubsystem().addEventReporter(handler); }

inline System::operator const Subsystem&() const {return getDefaultSubsystem();}
inline System::operator Subsystem&() {return updDefaultSubsystem();}

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_SYSTEM_H_
