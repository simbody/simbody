#ifndef SimTK_SimTKCOMMON_SYSTEM_GUTS_H_
#define SimTK_SimTKCOMMON_SYSTEM_GUTS_H_

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
class System;

// See below for definitions.
static void systemDestructImplLocator(System::Guts*);
static System::Guts* systemCloneImplLocator(const System::Guts&);
static int systemRealizeTopologyImplLocator(const System::Guts&, State&);
static int systemRealizeModelImplLocator   (const System::Guts&, State&);
static int systemRealizeInstanceImplLocator(const System::Guts&, const State&);
static int systemRealizeTimeImplLocator    (const System::Guts&, const State&);
static int systemRealizePositionImplLocator(const System::Guts&, const State&);
static int systemRealizeVelocityImplLocator(const System::Guts&, const State&);
static int systemRealizeDynamicsImplLocator(const System::Guts&, const State&);
static int systemRealizeAccelerationImplLocator(const System::Guts&, const State&);
static int systemRealizeReportImplLocator      (const System::Guts&, const State&);
static Real systemCalcTimescaleImplLocator(const System::Guts&, const State&);
static int  systemCalcYUnitWeightsImplLocator(const System::Guts&, const State&, Vector& weights);
static int  systemProjectImplLocator(const System::Guts&, State&, Real, const Vector&, const Vector&,
                                         Vector&);
static int  systemCalcYErrUnitTolerancesImplLocator(const System::Guts&, const State&, Vector& ootols);
static int  systemHandleEventsImplLocator(const System::Guts&, State&, System::EventCause, const Array<int>&,
                                              Real, const Vector&, const Vector&, Stage&, bool&);
static int  systemCalcEventTriggerInfoImplLocator(const System::Guts&, const State&, Array<System::EventTriggerInfo>&);
static int  systemCalcTimeOfNextScheduledEventImplLocator(const System::Guts&, const State&, Real&, Array<int>&);

/**
 * This is the declaration for the System::Guts class, the abstract object to
 * which a System handle points. This is in a separate header file from System
 * because only people who are extending the System class to make their own
 * Systems need to be aware of the details. End users access only methods from
 * the System class and classes derived from System, never anything from
 * System::Guts or its derived classes.
 *
 * Below is the physical layout of memory for a System, and which
 * portions are allocated by the client program and which by the
 * binary library code. For binary compatiblity, only the side
 * which allocated a piece of memory can access it. So for example,
 * the client code can use the C++ System::Guts Virtual Function Table (VFT)
 * to call the concrete Guts methods. But the library side, when
 * calling those same methods, must go through its own explicitly-
 * managed VFT since it can't know what ordering was used for the
 * methods in the VFT on the client side.
 *
 *               CLIENT SIDE                    .  LIBRARY SIDE
 *                                              .
 *       System              System::Guts       . System::Guts::GutsRep
 *   ---------------       ------------------   .   -------------
 *  | System::Guts* | --> | System::GutsRep* | --> |   GutsRep   |
 *   ---------------       ------------------   .  |             |
 *          ^             | Concrete Guts    |  .  |  Position   |
 *          |             | class data and   |  .  | independent |
 *   ===============      | client-side VFT  |  .  |  Guts VFT   |
 *   Concrete System       ------------------   .  |             |
 *     adds no data                             .  | Other opaque|
 *       members                                .  |   stuff     |
 *                                              .   -------------
 *
 * If the concrete System::Guts class also has an opaque implementation,
 * as it will for concrete Systems provided by the SimTK Core, then
 * the System author should expose only the data-free handle class 
 * derived from System.
 */
class SimTK_SimTKCOMMON_EXPORT System::Guts {
    class GutsRep;
    friend class GutsRep;

    // This is the only data member in this class.
    GutsRep* rep; // opaque implementation of System::Guts base class.
public:
    // Constructor must be inline for binary compatibility. Note that this
    // serves as a default constructor since both arguments have defaults.
    inline explicit Guts(const String& name="<NONAME>", 
                         const String& version="0.0.0");

    // This won't be called directly from library-side code. Instead,
    // a method from the explicit virtual function table will be invoked
    // which will know where to find this on in the C++ VFT on the client side.
    virtual ~Guts() {librarySideDestruction();}

    const String& getName()    const;
    const String& getVersion() const;

    void setHasTimeAdvancedEvents(State&, bool hasEm) const;
    bool hasTimeAdvancedEvents(const State&) const;

        //////////////////////////////
        // EVALUATION (REALIZATION) //
        //////////////////////////////

    // These are the routines to which the System class forwards requests.

    const State& getDefaultState() const;
    State&       updDefaultState();

    void realize(const State& s, Stage g = Stage::HighestValid) const;

    SubsystemId adoptSubsystem(Subsystem& child);

    int getNSubsystems() const;
    const Subsystem& getSubsystem(SubsystemId)   const;
    Subsystem&       updSubsystem(SubsystemId);

    // Obtain the owner handle for this System::Guts object.
    const System& getSystem() const;
    System& updSystem();

    void setOwnerHandle(System&);
    bool hasOwnerHandle() const;

    explicit Guts(class GutsRep* r) : rep(r) { }
    bool           hasRep() const {return rep!=0;}
    const GutsRep& getRep() const {assert(rep); return *rep;}
    GutsRep&       updRep() const {assert(rep); return *rep;}

    bool systemTopologyHasBeenRealized() const;
    void invalidateSystemTopologyCache() const;

    // Call this routine to invoke the client-side virtual destructor,
    // by going through the library-side explicit virtual function table.
    static void destruct(System::Guts*);

    // Wrap the cloneImpl virtual method.
    System::Guts* clone() const;

    // These routines wrap the virtual realize...Impl() methods to ensure
    // good behavior such as checking that stage requirements are met and
    // updating the stage at the end. Note that these will do nothing if
    // the System stage is already at or greater than the indicated stage.

    const State& realizeTopology() const;
    void realizeModel(State&) const;
    void realizeInstance    (const State& s) const;
    void realizeTime        (const State& s) const;
    void realizePosition    (const State& s) const;
    void realizeVelocity    (const State& s) const;
    void realizeDynamics    (const State& s) const;
    void realizeAcceleration(const State& s) const;
    void realizeReport      (const State& s) const;

    // These wrap the other virtual methods.
    Real calcTimescale(const State&) const;
    void calcYUnitWeights(const State&, Vector& weights) const;
    void project(State&, Real consAccuracy, const Vector& yweights,
                 const Vector& ootols, Vector& yerrest) const;
    void calcYErrUnitTolerances(const State&, Vector& tolerances) const;
    void handleEvents
       (State&, EventCause, const Array<int>& eventIds,
        Real accuracy, const Vector& yWeights, const Vector& ooConstraintTols,
        Stage& lowestModified, bool& shouldTerminate) const;
    void calcEventTriggerInfo(const State&, Array<EventTriggerInfo>&) const;
    void calcTimeOfNextScheduledEvent
        (const State&, Real& tNextEvent, Array<int>& eventIds) const;

    void calcDecorativeGeometryAndAppend(const State&, Stage, 
                                         Array<DecorativeGeometry>&) const;


protected:
    Guts(const Guts&);  // copies the base class; for use from derived class copy constructors
    
    // The destructor is already virtual; see above.

    virtual System::Guts* cloneImpl() const = 0;

    // Override these to change the evaluation order of the Subsystems.
    // The default is to evaluate them in increasing order of SubsystemId.
    // These methods should not be called directly; they are invoked by the
    // above wrapper methods. Note: the wrappers *will not* call these
    // routines if the system stage has already met the indicated stage level.
    // If fact these routines will be called only when the system stage
    // is at the level just prior to the one indicated here. For example,
    // realizeVelocityImpl() will be called only if the passed-in State
    // has been determined to have its system stage exactly Stage::Position.
    virtual int realizeTopologyImpl(State&) const;
    virtual int realizeModelImpl(State&) const;
    virtual int realizeInstanceImpl(const State&) const;
    virtual int realizeTimeImpl(const State&) const;
    virtual int realizePositionImpl(const State&) const;
    virtual int realizeVelocityImpl(const State&) const;
    virtual int realizeDynamicsImpl(const State&) const;
    virtual int realizeAccelerationImpl(const State&) const;
    virtual int realizeReportImpl(const State&) const;

    virtual Real calcTimescaleImpl(const State&) const;

    virtual int calcYUnitWeightsImpl(const State&, Vector& weights) const;

    virtual int projectImpl(State&, Real consAccuracy, const Vector& yweights,
                            const Vector& ootols, Vector& yerrest) const;
    virtual int calcYErrUnitTolerancesImpl(const State&, Vector& tolerances) const;

    virtual int handleEventsImpl
       (State&, EventCause, const Array<int>& eventIds,
        Real accuracy, const Vector& yWeights, const Vector& ooConstraintTols,
        Stage& lowestModified, bool& shouldTerminate) const;

    virtual int calcEventTriggerInfoImpl(const State&, Array<EventTriggerInfo>&) const;

    virtual int calcTimeOfNextScheduledEventImpl
        (const State&, Real& tNextEvent, Array<int>& eventIds) const;
private:
    Guts& operator=(const Guts&); // suppress default copy assignment operator

    // These typedefs are used internally to manage the binary-compatible
    // handling of the virtual function table.

    // This first entry calls the virtual destructor above to delete the
    // heap-allocated object pointed to by the passed-in pointer.
    typedef void (*DestructImplLocator)(System::Guts*);
    typedef System::Guts* (*CloneImplLocator)(const System::Guts&);

    typedef int (*RealizeWritableStateImplLocator)(const System::Guts&, State&);
    typedef int (*RealizeConstStateImplLocator)(const System::Guts&, const State&);
    typedef Real (*CalcTimescaleImplLocator)(const System::Guts&, const State&);
    typedef int (*CalcUnitWeightsImplLocator)(const System::Guts&, const State&, Vector& weights);

    typedef int (*ProjectImplLocator)(const System::Guts&, State&, Real, const Vector&, const Vector&,
                                             Vector&);

    typedef int (*HandleEventsImplLocator)
       (const System::Guts&, State&, EventCause, const Array<int>&,
        Real, const Vector&, const Vector&, Stage&, bool&);
    typedef int (*CalcEventTriggerInfoImplLocator)
       (const System::Guts&, const State&, Array<EventTriggerInfo>&);
    typedef int (*CalcTimeOfNextScheduledEventImplLocator)
       (const System::Guts&, const State&, Real&, Array<int>&);

    class EventTriggerInfoRep;

    void librarySideConstruction(const String& name, const String& version);
    void librarySideDestruction();

    void registerDestructImpl(DestructImplLocator);
    void registerCloneImpl(CloneImplLocator);

    void registerRealizeTopologyImpl    (RealizeWritableStateImplLocator);
    void registerRealizeModelImpl       (RealizeWritableStateImplLocator);
    void registerRealizeInstanceImpl    (RealizeConstStateImplLocator);
    void registerRealizeTimeImpl        (RealizeConstStateImplLocator);
    void registerRealizePositionImpl    (RealizeConstStateImplLocator);
    void registerRealizeVelocityImpl    (RealizeConstStateImplLocator);
    void registerRealizeDynamicsImpl    (RealizeConstStateImplLocator);
    void registerRealizeAccelerationImpl(RealizeConstStateImplLocator);
    void registerRealizeReportImpl      (RealizeConstStateImplLocator);

    void registerCalcTimescaleImpl(CalcTimescaleImplLocator);
    void registerCalcYUnitWeightsImplLocator(CalcUnitWeightsImplLocator);
    void registerProjectImpl(ProjectImplLocator);
    void registerCalcYErrUnitTolerancesImplLocator(CalcUnitWeightsImplLocator);
    void registerHandleEventsImpl(HandleEventsImplLocator);
    void registerCalcEventTriggerInfoImpl(CalcEventTriggerInfoImplLocator);
    void registerCalcTimeOfNextScheduledEventImpl(CalcTimeOfNextScheduledEventImplLocator);

    // We want the locator functions to have access to the protected "Impl"
    // virtual methods, so we make them friends.

    friend void systemDestructImplLocator(System::Guts*);
    friend System::Guts* systemCloneImplLocator(const System::Guts&);

    friend int systemRealizeTopologyImplLocator(const System::Guts&, State&);
    friend int systemRealizeModelImplLocator   (const System::Guts&, State&);
    friend int systemRealizeInstanceImplLocator(const System::Guts&, const State&);
    friend int systemRealizeTimeImplLocator    (const System::Guts&, const State&);
    friend int systemRealizePositionImplLocator(const System::Guts&, const State&);
    friend int systemRealizeVelocityImplLocator(const System::Guts&, const State&);
    friend int systemRealizeDynamicsImplLocator(const System::Guts&, const State&);
    friend int systemRealizeAccelerationImplLocator(const System::Guts&, const State&);
    friend int systemRealizeReportImplLocator      (const System::Guts&, const State&);

    friend Real systemCalcTimescaleImplLocator(const System::Guts&, const State&);
    friend int  systemCalcYUnitWeightsImplLocator(const System::Guts&, const State&, Vector& weights);
    friend int  systemProjectImplLocator(const System::Guts&, State&, Real, const Vector&, const Vector&,
                                         Vector&);
    friend int  systemCalcYErrUnitTolerancesImplLocator(const System::Guts&, const State&, Vector& ootols);
    friend int  systemHandleEventsImplLocator(const System::Guts&, State&, EventCause, const Array<int>&,
                                              Real, const Vector&, const Vector&, Stage&, bool&);
    friend int  systemCalcEventTriggerInfoImplLocator(const System::Guts&, const State&, Array<EventTriggerInfo>&);
    friend int  systemCalcTimeOfNextScheduledEventImplLocator(const System::Guts&, const State&, Real&, Array<int>&);
};


// These are used to supply the client-side virtual function to the library, without
// the client and library having to agree on the layout of the virtual function tables.

static void systemDestructImplLocator(System::Guts* sysp)
  { delete sysp; } // invokes virtual destructor
static System::Guts* systemCloneImplLocator(const System::Guts& sys)
  { return sys.cloneImpl(); }

static int systemRealizeTopologyImplLocator(const System::Guts& sys, State& state)
  { return sys.realizeTopologyImpl(state); }
static int systemRealizeModelImplLocator(const System::Guts& sys, State& state)
  { return sys.realizeModelImpl(state); }
static int systemRealizeInstanceImplLocator(const System::Guts& sys, const State& state)
  { return sys.realizeInstanceImpl(state); }
static int systemRealizeTimeImplLocator(const System::Guts& sys, const State& state)
  { return sys.realizeTimeImpl(state); }
static int systemRealizePositionImplLocator(const System::Guts& sys, const State& state)
  { return sys.realizePositionImpl(state); }
static int systemRealizeVelocityImplLocator(const System::Guts& sys, const State& state)
  { return sys.realizeVelocityImpl(state); }
static int systemRealizeDynamicsImplLocator(const System::Guts& sys, const State& state)
  { return sys.realizeDynamicsImpl(state); }
static int systemRealizeAccelerationImplLocator(const System::Guts& sys, const State& state)
  { return sys.realizeAccelerationImpl(state); }
static int systemRealizeReportImplLocator(const System::Guts& sys, const State& state)
  { return sys.realizeReportImpl(state); }


static Real systemCalcTimescaleImplLocator(const System::Guts& sys, const State& state)
  { return sys.calcTimescaleImpl(state); }

static int  systemCalcYUnitWeightsImplLocator(const System::Guts& sys, const State& state, Vector& weights)
  { return sys.calcYUnitWeightsImpl(state, weights); }

static int  systemProjectImplLocator
   (const System::Guts& sys, State& state, Real consAccuracy,
    const Vector& yWeights, const Vector& ooConstraintTols, Vector& yErrest)
  { return sys.projectImpl(state, consAccuracy, yWeights, ooConstraintTols, yErrest); }

static int  systemCalcYErrUnitTolerancesImplLocator(const System::Guts& sys, const State& state, Vector& ootols)
  { return sys.calcYErrUnitTolerancesImpl(state, ootols); }

static int  systemHandleEventsImplLocator
   (const System::Guts& sys, State& state, System::EventCause cause, const Array<int>& eventIds,
    Real accuracy, const Vector& yWeights, const Vector& ooConstraintTols, Stage& lowestModified, bool& shouldTerminate)
  { return sys.handleEventsImpl(state, cause, eventIds, accuracy, yWeights, ooConstraintTols, lowestModified, shouldTerminate); }

static int  systemCalcEventTriggerInfoImplLocator(const System::Guts& sys, const State& state, Array<System::EventTriggerInfo>& info)
  { return sys.calcEventTriggerInfoImpl(state,info); }

static int  systemCalcTimeOfNextScheduledEventImplLocator
   (const System::Guts& sys, const State& state, Real& tNextEvent, Array<int>& eventIds)
  { return sys.calcTimeOfNextScheduledEventImpl(state,tNextEvent,eventIds); }


// Constructor must be inline so that it has access to the above static
// functions which are private to the client-side compilation unit in which the
// client-side virtual function table is understood.
inline System::Guts::Guts(const String& name, const String& version) : rep(0)
{
    librarySideConstruction(name, version);

    // Teach the library code how to call client side virtual functions by
    // calling through the client side compilation unit's private static
    // locator functions.
    registerDestructImpl(systemDestructImplLocator);
    registerCloneImpl(systemCloneImplLocator);

    registerRealizeTopologyImpl    (systemRealizeTopologyImplLocator);
    registerRealizeModelImpl       (systemRealizeModelImplLocator);
    registerRealizeInstanceImpl    (systemRealizeInstanceImplLocator);
    registerRealizeTimeImpl        (systemRealizeTimeImplLocator);
    registerRealizePositionImpl    (systemRealizePositionImplLocator);
    registerRealizeVelocityImpl    (systemRealizeVelocityImplLocator);
    registerRealizeDynamicsImpl    (systemRealizeDynamicsImplLocator);
    registerRealizeAccelerationImpl(systemRealizeAccelerationImplLocator);
    registerRealizeReportImpl      (systemRealizeReportImplLocator);


    registerCalcTimescaleImpl(systemCalcTimescaleImplLocator);
    registerCalcYUnitWeightsImplLocator(systemCalcYUnitWeightsImplLocator);
    registerProjectImpl(systemProjectImplLocator);
    registerCalcYErrUnitTolerancesImplLocator(systemCalcYErrUnitTolerancesImplLocator);
    registerHandleEventsImpl(systemHandleEventsImplLocator);
    registerCalcEventTriggerInfoImpl(systemCalcEventTriggerInfoImplLocator);
    registerCalcTimeOfNextScheduledEventImpl(systemCalcTimeOfNextScheduledEventImplLocator);
}

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_SYSTEM_GUTS_H_
