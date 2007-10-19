#ifndef SimTK_SimTKCOMMON_SYSTEM_H_
#define SimTK_SimTKCOMMON_SYSTEM_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-7 Stanford University and the Authors.         *
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

namespace SimTK {

class Subsystem;
class DecorativeGeometry;
class DefaultSystemSubsystem;

/**
 * The handle class which serves as the abstract parent of all Systems.
 *
 * A System serves as a mediator for a group of interacting Subsystems.
 * All will share a single system State, and typically subsystems will
 * need access to content in the state which is produced by other
 * subsystems. 
 *
 * A System provides a unique SubsystemId (a small positive integer)
 * for each of its subsystems, and the subsystems are constructed
 * knowing their indices. The indices are used subsequently by the subsystems
 * to find their own entries in the system state, and by each subsystem
 * to refer to others within the same system. Index 0 is reserved for 
 * use by the System itself, e.g. for system-global state variables.
 *
 * Concrete Systems understand the kinds of subsystems they contain. 
 * For example, a MultibodySystem might contain a mechanical subsystem,
 * a force subsystem, and a geometry subsystem. At each computation
 * stage, a subsystem is realized in a single operation. That operation
 * can refer to computations from already-realized subsystems, but
 * cannot initiate computation in other subsystems. The System must
 * know the proper order with which to realize the subsystems at each
 * stage, and that ordering is likely to vary with stage. For example,
 * at Position stage the mechanical positions must be realized
 * before the configuration-dependent force elements. However, at
 * Acceleration stage, the force elements must be realized before the
 * mechanical accelerations can be calculated.
 *
 * There are two distinct users of this class:
 *   - System Users: people who are making use of a concrete System (which will
 *     inherit methods from this class)
 *   - System Developers: people who are writing concrete System classes
 * Note that System Users include people who are writing
 * Studies, Reporters, Modelers and so on as well as end users who
 * are accessing the System directly.
 *
 * Only methods intended for System Users and a few bookkeeping methods
 * are in the main System class, which is a SimTK Handle class, meaning
 * that it consists only of a single pointer, which points to a 
 * System::Guts class. The Guts class is abstract, and virtual methods
 * to be implemented by System Developers in the concrete
 * System are defined there, along
 * with other utilities of use to the concrete System Developer but
 * not to the end user. The Guts class is declared in a separate
 * header file, and only people who are writing their own System
 * classes need look there.
 */
class SimTK_SimTKCOMMON_EXPORT System {
public:
    class Guts; // local; name is System::Guts
    friend class Guts;
private:
    // This is the only data member in this class. Also, any class derived from
    // System must have *NO* data members at all (data goes in the Guts class).
    Guts* guts;
public:
    System() : guts(0) { }
    System(const System&);
    System& operator=(const System&);
    ~System();

    const String& getName()    const;
    const String& getVersion() const;


        /////////////////
        // REALIZATION //
        /////////////////

    /// The following call must be made after any topological change
    /// has been made to this System, before the System can be used
    /// to perform any computations. Perhaps surprisingly, the method
    /// is const. That's because the topology cannot be changed by this method.
    /// Various mutable "cache" entries will get calculated, including the
    /// default State, a reference to which is returned.
    /// The returned State has already been realized through the Model
    /// Stage, using the defaults for the Model-stage variables, meaning
    /// that all later stage variables have been allocated and set to
    /// their default values as well. You can access
    /// this same default State again using getDefaultState().
    /// If the current topology has already been realized, this call does 
    /// nothing but return a reference to the already-built default State.
    const State& realizeTopology() const;


    /// This is available after realizeTopology(), and will throw an
    /// exception if realizeTopology() has not been called since the
    /// most recent topological change to this System. This method returns the
    /// same reference returned by realizeTopology(). The State to which
    /// a reference is returned was created by the most recent
    /// realizeTopology() call. It has already been realized through the
    /// Model Stage, using default values for all the Model-stage variables.
    /// All later-stage variables have been allocated and set to their
    /// default values. You can use this state directly to obtain information
    /// about the System in its default state or you can use this state
    /// to initialize other States to which you have write access. Those
    /// States are then suitable for further computation with this System.
    const State& getDefaultState() const;
    State&       updDefaultState();

    /// This call is required if Model-stage variables are changed from
    /// their default values. The System topology must already have been
    /// realized (that is, realizeTopology() must have been called since
    /// the last topological change made to the System). Also, the supplied
    /// State must already have been initialized to work with this System
    /// either by copying the default state or some other State of this System.
    /// If it has already been realized to Stage::Model or higher, nothing
    /// happens here. Otherwise, all the state variables at Stage::Instance or
    /// higher are allocated or reallocated (if necessary), and reinitialized
    /// to their default values. NOTE: any State information at Stage::Instance
    /// or higher in the passed-in State is *destroyed* here. The number, types
    /// and memory locations of those state variables will change, so any
    /// existing references or pointers to them are invalid after this call.
    /// Note that this routine modifies its State argument, but makes no changes
    /// at all to the System itself and is hence const.
    void realizeModel(State&) const;

    /// Realize the entire System to the indicated Stage. The passed-in
    /// State must have been initialized to work with this System, and
    /// it must already have been realized through Stage::Model, since
    /// the realize() method doesn't have write access to the State.
    /// If the state has already been realized to the requested stage
    /// or higher, nothing happens here. Otherwise, the state is realized
    /// one stage at a time until it reaches the requested stage. 
    void realize(const State& s, Stage g = Stage::HighestValid) const;

    // Generate all decorative geometry computable at a specific stage. This will
    // throw an exception if the state hasn't already been realized
    // to that stage. Note that the list is not inclusive -- you have to
    // request geometry from each stage to get all of it.
    // This routine asks each subsystem in succession to generate its decorative geometry
    // and append it to the end of the Array.
    // If the stage is Stage::Topology, realizeTopology() must already have
    // been called but the State is ignored.
    void calcDecorativeGeometryAndAppend(const State&, Stage, Array<DecorativeGeometry>&) const;

        ///////////////////////////
        // THE CONTINUOUS SYSTEM //
        ///////////////////////////

            // UNCONSTRAINED

    /// This operator can be called at Stage::Instance or higher and 
    /// returns a rough estimate of a length of time we consider significant
    /// for this system. For example, this could be the period of the highest-frequency
    /// oscillation that we care about. This can be used as a hint by 
    /// numerical integrators in choosing their initial step size, and 
    /// suggests how velocity variables should be scaled relative to their
    /// corresponding position variables.
    Real calcTimescale(const State&) const;

    /// This operator can be called at Stage::Position to calculate a weighting
    /// vector w, with one entry for each state variable y={q,u,z}, ordered
    /// the same as Y in the State and calculated specifically for the current
    /// values of Y in the State. Weight wi is proportional to the "importance"
    /// of state variable yi with respect to some criteria determined by
    /// the System, such that wi*dyi=1 indicates that a change dyi in state
    /// yi produces approximately a unit change in the weighting criteria. This
    /// is intended for use by numerical integration methods for step size control.
    /// The idea is to allow creation of a weighted RMS norm which returns 1 just when
    /// all the state variable changes have a unit effect. The norm is RMS(W*dy) where
    /// W=diag(w). A value of 1 for this norm would typically be a huge error. 
    /// For example, if your accuracy requirement is 0.1%, you would test that
    /// the weighted RMS norm is <= .001.
    /// We expect this operation to be fairly expensive and thus the integrator
    /// is expected to invoke it only occasionally.
    void calcYUnitWeights(const State&, Vector& weights) const;

            // CONSTRAINED

    /// This optional solver projects the given State back on to the constraint
    /// manifold, by the shortest path possible in the weighted norm
    /// given by the supplied weights, satisfying the constraints by
    /// reducing the supplied tolerance norm to below consAccuracy.
    /// May also project out the constraint-normal component of the
    /// passed-in error estimate vector yerrest.
    /// This is part of the integration of the continuous DAE system and
    /// thus should never require an integrator restart. The System
    /// author must ensure that only position and velocity stage, continuous
    /// variables are updated by this call.
    /// On return the state will be realized to at least Stage::Velocity.
    void project(State&, Real consAccuracy, const Vector& yweights,
                 const Vector& ootols, Vector& yerrest) const;

    /// This provides scaling information for each of the position and velocity
    /// constraints (YErr) in the State. The tolerance is the absolute error in the
    /// constraint which is considered a "unit violation" of that state. Then
    /// if T=diag(tol) and c the vector of constraint errors, we can use a
    /// weighted RMS norm condition like RMS(T*c) <= accuracy to define
    /// when constraints have been adequately met.
    /// This is expected to be a cheap operation but not to change during
    /// a study. State must be realized to Stage::Model.
    void calcYErrUnitTolerances(const State&, Vector& tolerances) const;


        /////////////////////////
        // THE DISCRETE SYSTEM //
        /////////////////////////

    class EventTriggerInfo;

    /// This determines whether this System wants to be notified whenever time
    /// advances irreversibly. If set true, time advancement is treated as an
    /// event. Otherwise, time advancement proceeds silently.
    /// TODO: currently not using State so this is a Topology stage variable,
    /// but should probably be Model stage.
    void setHasTimeAdvancedEvents(bool); // default=false
    bool hasTimeAdvancedEvents() const;

    /// These are all the possible causes for events. (1) An event trigger
    /// function may have undergone a monitored sign transition, (2) we
    /// may have reached the time for a scheduled event, (3) the integrator
    /// may have completed an internal step (assuming we registered our
    /// interest in those), (4) or a termination condition may have been
    /// detected. In case several of these causes are detected in 
    /// a single step, they are sequentialized in the order shown,
    /// like this:
    ///    1. The occurrence of triggered events is reported and
    ///       the triggering state and a list of triggered events
    ///       are passed to the event handler for processing (meaning
    ///       the state, but not the time, is modified). [Note that
    ///       simultaneity *within* the set of triggered events may
    ///       also require special handling; we're not talking about
    ///       that here, just simultaneity of *causes*.]
    ///    2. Next, using the state resulting from step 1, the time is checked
    ///       to see if scheduled events have occurred. If so, a list of
    ///       those events is passed to the event handler for processing.
    ///    3. Next, if this system has requested time-advanced events,
    ///       the event handler is called with the state that resulted
    ///       from step 2 and the "time advanced" cause noted. No event
    ///       list is passed in that case. The state may be modified.
    ///    4. Last, if the final time has been reached or if any of
    ///       the event handlers asked for termination, we pass the
    ///       state to the event handler again noting that we have
    ///       reached termination. The state may be modified and the
    ///       result will be the final state of the simulation.

    enum EventCause {
        TriggeredEvents    =1,
        ScheduledEvents    =2,
        TimeAdvancedEvent  =3,
        TerminationEvent   =4,

        InvalidEventCause  = -1
    };
    static const char* getEventCauseName(EventCause);


    /// This solver handles a set of events which a TimeStepper has
    /// denoted as having occurred. The event handler may make discontinuous
    /// changes in the State, in general both to discrete and continuous variables,
    /// but NOT to time. It cannot change topological information.
    /// If changes are made to continuous variables, the handler
    /// is required to make sure the returned state satisfies the constraints to the
    /// indicated accuracy level.
    /// On return, the handleEvents routine should set the output variable
    /// lowestModified to the Stage level of the lowest-stage variable
    /// it modified. This information tells the time stepper how much of 
    /// a restart it must perform on the underlying numerical integrator. 
    /// When in doubt, set lowestModified to Stage::Model, which will cause
    /// a complete restart.
    /// Finally, if the handler determines that the occurrence of some event
    /// requires that the simulation be terminated it should set 'shouldTerminate'
    /// to true before returning.

    void handleEvents
       (State&, EventCause, const Array<int>& eventIds,
        Real accuracy, const Vector& yWeights, const Vector& ooConstraintTols,
        Stage& lowestModified, bool& shouldTerminate) const;
    
    /// This method is similar to handleEvents(), but does not allow the State to be
    /// modified.  It is used for scheduled events that were marked as being reports.
    
    void reportEvents(const State& s, EventCause cause, const Array<int>& eventIds) const;

    /// This routine provides the Integrator with information it needs about the
    /// individual event trigger functions, such as which sign transitions are
    /// relevant and how tightly we need to localize. This is considered 
    /// Instance stage information so cannot change during a continuous integration
    /// interval (so an Integrator can process it upon restart(Instance)), however it can be
    /// updated whenever a discrete update is made to the State.
    /// A default implementation is provided which returns default EventTriggerInfo for
    /// each event trigger in State.
    /// State must already be realized to Stage::Instance.
    void calcEventTriggerInfo(const State&, Array<EventTriggerInfo>&) const;

    /// This routine should be called to determine if and when there is an event
    /// scheduled to occur at a particular time. This is a *lot* cheaper than
    /// making the Integrator hunt these down like ordinary state-dependent events.
    /// The returned time can be passed to the Integrator's stepping function as
    /// the advance time limit.
    /// If this routine sets isReport to true, the event should be handled by
    /// calling reportEvents() instead of handleEvents().
    void calcTimeOfNextScheduledEvent
        (const State&, Real& tNextEvent, Array<int>& eventIds, bool& isReport) const;


    //TODO: these operators should be provided by the Vector class where they
    //can be perfomed more efficiently.

    static Real calcWeightedRMSNorm(const Vector& values, const Vector& weights) {
        assert(weights.size() == values.size());
        if (values.size()==0) return 0;
        Real sumsq = 0;
        for (int i=0; i<values.size(); ++i) {
            const Real wv = weights[i]*values[i];
            sumsq += wv*wv;
        }
        return std::sqrt(sumsq/weights.size());
    }

    static Real calcWeightedInfinityNorm(const Vector& values, const Vector& weights) {
        assert(weights.size() == values.size());
        if (values.size()==0) return 0;
        Real maxval = 0;
        for (int i=0; i<values.size(); ++i) {
            const Real wv = std::abs(weights[i]*values[i]);
            if (wv > maxval) maxval=wv;
        }
        return maxval;
    }

    /// Take over ownership of the supplied subsystem and install it into 
    /// the next free subsystem slot. The new slot index is returned.
    SubsystemId adoptSubsystem(Subsystem& child);

    /// How may Subsystems are in here?
    int getNSubsystems() const;
    /// Obtain read-only access to a particular subsystem by its index.
    const Subsystem& getSubsystem(SubsystemId)   const;
    /// Obtain writable access to a particular subsystem by its index.
    Subsystem&       updSubsystem(SubsystemId);
    /// Get read-only access to the default subsystem which is present in every system.
    const DefaultSystemSubsystem& getDefaultSubsystem() const;
    /// Get writable access to the default subsystem which is present in every system.
    DefaultSystemSubsystem& updDefaultSubsystem();

    // Internal use only
    bool isOwnerHandle() const;
    bool isEmptyHandle() const;

    // There can be multiple handles on the same System.
    bool isSameSystem(const System& otherSystem) const;

    /// You can check whether realizeTopology() has been called since the last
    /// topological change to this Syatem. If you don't check and just plunge
    /// ahead you are likely to encounter an exception since very few things
    /// will work without topology having been realized.
    bool systemTopologyHasBeenRealized() const;

    // dynamic_cast the returned reference to a reference to your concrete Guts
    // class.
    const System::Guts& getSystemGuts() const {assert(guts); return *guts;}
    System::Guts&       updSystemGuts()       {assert(guts); return *guts;}

    // Put new *unowned* Guts into this *empty* handle and take over ownership.
    // If this handle is already in use, or if Guts is already owned this
    // routine will throw an exception.
    void adoptSystemGuts(System::Guts* g);

    explicit System(System::Guts* g) : guts(g) { }
    bool hasGuts() const {return guts!=0;}

private:
    class EventTriggerInfoRep;
};


/// This class is used to communicate between the System and an 
/// Integrator regarding the properties of a particular event trigger
/// function. Currently these are:
///   - Whether to watch for rising sign transitions, falling, or both. [BOTH]
///   - Whether to watch for transitions to and from zero. [NO]
///   - The localization window in units of the System's timescale. [10%]
///     (That is then the "unit" window which is reduced by the
///      accuracy setting.)
/// The default values are shown in brackets above.
///
class SimTK_SimTKCOMMON_EXPORT System::EventTriggerInfo {
public:
    EventTriggerInfo();
    explicit EventTriggerInfo(int eventId);
    ~EventTriggerInfo();
    EventTriggerInfo(const EventTriggerInfo&);
    EventTriggerInfo& operator=(const EventTriggerInfo&);

    int getEventId() const; // returns -1 if not set
    bool shouldTriggerOnRisingSignTransition()  const; // default=true
    bool shouldTriggerOnFallingSignTransition() const; // default=true
    bool shouldTriggerOnZeroTransitions()       const; // default=false
    Real getRequiredLocalizationTimeWindow()    const; // default=0.1

    // These return the modified 'this', like assignment operators.
    EventTriggerInfo& setEventId(int);
    EventTriggerInfo& setTriggerOnRisingSignTransition(bool);
    EventTriggerInfo& setTriggerOnFallingSignTransition(bool);
    EventTriggerInfo& setTriggerOnZeroTransitions(bool);
    EventTriggerInfo& setRequiredLocalizationTimeWindow(Real);

    // TODO: switch to SignTransitionSet and trash EventStatus altogether.
    EventStatus::EventTrigger calcTransitionMask() const {
        unsigned mask = 0;
        if (shouldTriggerOnRisingSignTransition()) {
            mask |= (EventStatus::NegativeToPositive|EventStatus::NegativeToZero);
            if (shouldTriggerOnZeroTransitions())
                mask |= EventStatus::ZeroToPositive;
        }
        if (shouldTriggerOnFallingSignTransition()) {
            mask |= (EventStatus::PositiveToNegative|EventStatus::PositiveToZero);
            if (shouldTriggerOnZeroTransitions())
                mask |= EventStatus::ZeroToNegative;
        }
        return EventStatus::EventTrigger(mask);
    }

    EventStatus::EventTrigger calcTransitionToReport
       (EventStatus::EventTrigger transitionSeen) const
    {
        if (shouldTriggerOnZeroTransitions())
            return transitionSeen; // what we saw is what we report
        // otherwise report -1 to 1 or 1 to -1 as appropriate
        if (transitionSeen & EventStatus::Rising)
            return EventStatus::NegativeToPositive;
        if (transitionSeen & EventStatus::Falling)
            return EventStatus::PositiveToNegative;
        assert(!"impossible event transition situation");
        return EventStatus::NoEventTrigger;
    }

private:
    // opaque implementation for future binary compatibility
    System::EventTriggerInfoRep* rep;

    const System::EventTriggerInfoRep& getRep() const {assert(rep); return *rep;}
    System::EventTriggerInfoRep&       updRep()       {assert(rep); return *rep;}
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_SYSTEM_H_
