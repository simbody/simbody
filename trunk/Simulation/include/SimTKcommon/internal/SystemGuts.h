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
#include "SimTKcommon/internal/System.h"

namespace SimTK {

class Subsystem;
class DecorativeGeometry;

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
 * which allocated a piece of memory can access it. Exception: both
 * the client and library side must agree on the virtual function
 * table (VFT) ordering of the client's virtual functions.
 * @verbatim
 *               CLIENT SIDE                    .  LIBRARY SIDE
 *                                              .
 *       System              System::Guts       . System::Guts::GutsRep
 *   ---------------       ------------------   .   -------------
 *  | System::Guts* | --> | System::GutsRep* | --> |   GutsRep   |
 *   ---------------       ------------------   .  |             |
 *          ^             | Concrete Guts    |  .  |  Opaque     |
 *          |             | class data and   |  .  |  stuff      |
 *   ===============      | virt func table  |  .  |             |
 *   Concrete System       ------------------   .  |             |
 *     adds no data                             .   -------------
 *       members
 * @endverbatim
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
    // Note that this serves as a default constructor since both arguments have defaults.
    explicit Guts(const String& name="<NONAME>", 
                  const String& version="0.0.0");
    virtual ~Guts();

    const String& getName()    const;
    const String& getVersion() const;

    void setHasTimeAdvancedEvents(bool hasEm);
    bool hasTimeAdvancedEvents() const;

        //////////////////////////////
        // EVALUATION (REALIZATION) //
        //////////////////////////////

    // These are the routines to which the System class forwards requests.

    const State& getDefaultState() const;
    State&       updDefaultState();

    void realize(const State& s, Stage g = Stage::HighestRuntime) const;

    SubsystemIndex adoptSubsystem(Subsystem& child);

    int getNSubsystems() const;
    const Subsystem& getSubsystem(SubsystemIndex)   const;
    Subsystem&       updSubsystem(SubsystemIndex);

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
    void prescribe(State&, Stage) const;
    void project(State&, Real consAccuracy, const Vector& yweights,
                 const Vector& ootols, Vector& yerrest, System::ProjectOptions) const;
    void calcYErrUnitTolerances(const State&, Vector& tolerances) const;
    void handleEvents
        (State&, Event::Cause, const std::vector<EventId>& eventIds,
        Real accuracy, const Vector& yWeights, const Vector& ooConstraintTols,
        Stage& lowestModified, bool& shouldTerminate) const;
    void reportEvents(const State&, Event::Cause, const std::vector<EventId>& eventIds) const;
    void calcEventTriggerInfo(const State&, std::vector<EventTriggerInfo>&) const;
    void calcTimeOfNextScheduledEvent(const State&, Real& tNextEvent, std::vector<EventId>& eventIds, bool includeCurrentTime) const;
    void calcTimeOfNextScheduledReport(const State&, Real& tNextEvent, std::vector<EventId>& eventIds, bool includeCurrentTime) const;

    void calcDecorativeGeometryAndAppend(const State&, Stage, 
                                         std::vector<DecorativeGeometry>&) const;


protected:
    Guts(const Guts&);  // copies the base class; for use from derived class copy constructors
    
    // The destructor is already virtual; see above.

    virtual System::Guts* cloneImpl() const = 0;

    // Override these to change the evaluation order of the Subsystems.
    // The default is to evaluate them in increasing order of SubsystemIndex.
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

    virtual int prescribeImpl(State&, Stage) const;
    virtual int projectImpl(State&, Real consAccuracy, const Vector& yweights,
                            const Vector& ootols, Vector& yerrest, System::ProjectOptions) const;
    virtual int calcYErrUnitTolerancesImpl(const State&, Vector& tolerances) const;

    virtual int handleEventsImpl
        (State&, Event::Cause, const std::vector<EventId>& eventIds,
        Real accuracy, const Vector& yWeights, const Vector& ooConstraintTols,
        Stage& lowestModified, bool& shouldTerminate) const;

    virtual int reportEventsImpl(const State&, Event::Cause, const std::vector<EventId>& eventIds) const;

    virtual int calcEventTriggerInfoImpl(const State&, std::vector<EventTriggerInfo>&) const;

    virtual int calcTimeOfNextScheduledEventImpl(const State&, Real& tNextEvent, std::vector<EventId>& eventIds, bool includeCurrentTime) const;
    virtual int calcTimeOfNextScheduledReportImpl(const State&, Real& tNextEvent, std::vector<EventId>& eventIds, bool includeCurrentTime) const;
private:
    Guts& operator=(const Guts&); // suppress default copy assignment operator

    class EventTriggerInfoRep;

};


} // namespace SimTK

#endif // SimTK_SimTKCOMMON_SYSTEM_GUTS_H_
