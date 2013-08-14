#ifndef SimTK_SimTKCOMMON_SYSTEM_GUTS_H_
#define SimTK_SimTKCOMMON_SYSTEM_GUTS_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
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
 * <pre>
 *               CLIENT SIDE                    |  LIBRARY SIDE
 *                                              |
 *       System              System::Guts       | System::Guts::GutsRep
 *   ---------------       ------------------   |   -------------
 *  | System::Guts* | --> | System::GutsRep* | --> |   GutsRep   |
 *   ---------------       ------------------   |  |             |
 *          ^             | Concrete Guts    |  |  |  Opaque     |
 *          |             | class data and   |  |  |  stuff      |
 *   ===============      | virt func table  |  |  |             |
 *   Concrete System       ------------------   |  |             |
 *     adds no data                             |   -------------
 *       members
 * </pre>
 *
 * If the concrete System::Guts class also has an opaque implementation,
 * as it will for concrete Systems provided by Simbody, then
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

    int getNumSubsystems() const;
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
    StageVersion getSystemTopologyCacheVersion() const;
    void setSystemTopologyCacheVersion(StageVersion topoVersion) const;
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
    void multiplyByN(const State& state, const Vector& u, 
                     Vector& dq) const;
    void multiplyByNTranspose(const State& state, const Vector& fq, 
                              Vector& fu) const;
    void multiplyByNPInv(const State& state, const Vector& dq, 
                         Vector& u) const;
    void multiplyByNPInvTranspose(const State& state, const Vector& fu, 
                                  Vector& fq) const;

    bool prescribeQ(State&) const;
    bool prescribeU(State&) const;

    void projectQ(State&, Vector& qErrEst, 
                  const ProjectOptions& options, ProjectResults& results) const;
    void projectU(State&, Vector& uErrEst, 
                  const ProjectOptions& options, ProjectResults& results) const;

    void handleEvents
        (State&, Event::Cause, const Array_<EventId>& eventIds,
         const HandleEventsOptions& options,
         HandleEventsResults& results) const;
    void reportEvents(const State&, Event::Cause, const Array_<EventId>& eventIds) const;
    void calcEventTriggerInfo(const State&, Array_<EventTriggerInfo>&) const;
    void calcTimeOfNextScheduledEvent(const State&, Real& tNextEvent, Array_<EventId>& eventIds, bool includeCurrentTime) const;
    void calcTimeOfNextScheduledReport(const State&, Real& tNextEvent, Array_<EventId>& eventIds, bool includeCurrentTime) const;

    void calcDecorativeGeometryAndAppend(const State&, Stage, 
                                         Array_<DecorativeGeometry>&) const;


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
    virtual int realizeTopologyImpl(State& state)       const {return 0;}
    virtual int realizeModelImpl   (State& state)       const {return 0;}
    virtual int realizeInstanceImpl(const State& state) const {return 0;}
    virtual int realizeTimeImpl    (const State& state) const {return 0;}
    virtual int realizePositionImpl(const State& state) const {return 0;}
    virtual int realizeVelocityImpl(const State& state) const {return 0;}
    virtual int realizeDynamicsImpl(const State& state) const {return 0;}
    virtual int realizeAccelerationImpl(const State& state) const {return 0;}
    virtual int realizeReportImpl  (const State& state) const {return 0;}

    virtual void multiplyByNImpl(const State& state, const Vector& u, 
                                 Vector& dq) const;
    virtual void multiplyByNTransposeImpl(const State& state, const Vector& fq, 
                                          Vector& fu) const;
    virtual void multiplyByNPInvImpl(const State& state, const Vector& dq, 
                                     Vector& u) const;
    virtual void multiplyByNPInvTransposeImpl(const State& state, const Vector& fu, 
                                              Vector& fq) const;

    // Defaults assume no prescribed motion; hence, no change made.
    virtual bool prescribeQImpl(State&) const {return false;}
    virtual bool prescribeUImpl(State&) const {return false;}

    // Defaults assume no constraints and return success meaning "all 
    // constraints satisfied".
    virtual void projectQImpl(State& state, Vector& qErrEst, 
             const ProjectOptions& options, ProjectResults& results) const
    {   results.clear(); results.setExitStatus(ProjectResults::Succeeded); }
    virtual void projectUImpl(State& state, Vector& uErrEst, 
             const ProjectOptions& options, ProjectResults& results) const
    {   results.clear(); results.setExitStatus(ProjectResults::Succeeded); }

    virtual void handleEventsImpl
       (State& state, Event::Cause cause, const Array_<EventId>& eventIds,
        const HandleEventsOptions& options, HandleEventsResults& results) const;

    virtual int reportEventsImpl(const State& state, Event::Cause cause, 
                                 const Array_<EventId>& eventIds) const;

    virtual int calcEventTriggerInfoImpl(const State& state, 
                                         Array_<EventTriggerInfo>& info) const;

    virtual int calcTimeOfNextScheduledEventImpl
       (const State& state, Real& tNextEvent, Array_<EventId>& eventIds, 
        bool includeCurrentTime) const;
    virtual int calcTimeOfNextScheduledReportImpl
       (const State& state, Real& tNextEvent, Array_<EventId>& eventIds, 
        bool includeCurrentTime) const;

private:
    Guts& operator=(const Guts&); // suppress default copy assignment operator

    class EventTriggerInfoRep;

};


} // namespace SimTK

#endif // SimTK_SimTKCOMMON_SYSTEM_GUTS_H_
