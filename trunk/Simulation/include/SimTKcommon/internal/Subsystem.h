#ifndef SimTK_SimTKCOMMON_SUBSYSTEM_H_
#define SimTK_SimTKCOMMON_SUBSYSTEM_H_

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
#include "SimTKcommon/internal/Measure.h"
#include "SimTKcommon/internal/System.h"

#include <cassert>

namespace SimTK {

class System;
class DecorativeGeometry;
class DefaultSystemSubsystemGuts;
class ScheduledEventHandler;
class ScheduledEventReporter;
class TriggeredEventHandler;
class TriggeredEventReporter;

/**
 * The abstract parent of all Subsystems.
 *
 * A Subsystem is expected to be part of a larger System and to have
 * interdependencies with other subsystems of that same system. It
 * must NOT have dependencies on objects which are outside the System.
 * Consequently construction of any concrete subsystem requires
 * specification of a system at that time.
 * Subsystems go through an extended construction phase in which
 * their contents and interdependencies are created. Thus all
 * of a System's Subsystems generally need to be available simultaneously 
 * during construction, so that they can reference each other.
 *
 * There are three distinct users of this class:
 *    - the System class
 *    - the concrete Subsystems derived from this class
 *    - the end user of a concrete Subsystem
 * Only end user methods are public here. Methods intended for
 * use by the concrete Subsystem class can be found in the Subsystem::Guts
 * class which is defined in a separate header file. End users need not
 * look over there -- trust me, you'll find it disturbing if you do!
 */
class SimTK_SimTKCOMMON_EXPORT Subsystem {
public:
    class Guts; // local; name is Subsystem::Guts
    friend class Guts;
private:
    // This is the only data member in this class. Also, any class derived from
    // Subsystem must have *NO* data members at all (data goes in the Guts class).
    Guts* guts;
public:
    Subsystem() : guts(0) { } // an empty handle
    Subsystem(const Subsystem&);
    Subsystem& operator=(const Subsystem&);
    ~Subsystem();

    const String& getName()    const;
    const String& getVersion() const;

    // These call the corresponding State method, supplying this Subsystem's
    // SubsystemIndex. The returned indices are local to this Subsystem.
    QIndex allocateQ(State&, const Vector& qInit) const;
    UIndex allocateU(State&, const Vector& uInit) const;
    ZIndex allocateZ(State&, const Vector& zInit) const;
    DiscreteVariableIndex allocateDiscreteVariable(State&, Stage, AbstractValue* v) const;

    CacheEntryIndex allocateCacheEntry   (const State&, Stage dependsOn, Stage computedBy, AbstractValue* v) const;
    CacheEntryIndex allocateCacheEntry   (const State& state, Stage g, AbstractValue* v) const 
    {   return allocateCacheEntry(state, g, g, v); }
    QErrIndex allocateQErr         (const State&, int nqerr) const;
    UErrIndex allocateUErr         (const State&, int nuerr) const;
    UDotErrIndex allocateUDotErr      (const State&, int nudoterr) const;
    EventTriggerByStageIndex allocateEventTriggersByStage(const State&, Stage, int ntriggers) const;

    // These return views on State shared global resources. The views
    // are private to this subsystem, but the global resources themselves
    // are not allocated until the *System* advances to stage Model.
    // Note that there is no subsystem equivalent of the State's "y"
    // vector because in general a subsystem's state variables will
    // not be contiguous. However, a subsystem's q's, u's, and z's
    // will all be contiguous within those arrays.
    const Vector& getQ(const State&) const;
    const Vector& getU(const State&) const;
    const Vector& getZ(const State&) const;
    const Vector& getQDot(const State&) const;
    const Vector& getUDot(const State&) const;
    const Vector& getZDot(const State&) const;
    const Vector& getQDotDot(const State&) const;
    const Vector& getQErr(const State&) const;
    const Vector& getUErr(const State&) const;
    const Vector& getUDotErr(const State&) const;
    const Vector& getMultipliers(const State&) const;
    const Vector& getEventTriggersByStage(const State&, Stage) const;

    // These return writable access to this subsystem's partition in the
    // State pool of continuous variables. These can be called at Stage::Model
    // or higher, and if necesary they invalidate the Position (q), Velocity (u),
    // or Dynamics (z) stage respectively.
    Vector& updQ(State&) const; // invalidates Stage::Position
    Vector& updU(State&) const; // invalidates Stage::Velocity
    Vector& updZ(State&) const; // invalidates Stage::Dynamics

    // For convenience.
    void setQ(State& s, const Vector& q) const {
        assert(q.size() == getNQ(s));
        updQ(s) = q;
    }
    void setU(State& s, const Vector& u) const {
        assert(u.size() == getNU(s));
        updU(s) = u;
    }
    void setZ(State& s, const Vector& z) const {
        assert(z.size() == getNZ(s));
        updZ(s) = z;
    }

    // These update the State cache which is mutable; hence, const State. They
    // can be called only if the previous stage has already been realized, e.g.,
    // updQDot() is allowed only while realizing the Velocity stage, requiring
    // that Position stage has already been realized.
    Vector& updQDot(const State&) const;
    Vector& updUDot(const State&) const;
    Vector& updZDot(const State&) const;
    Vector& updQDotDot(const State&) const;
    Vector& updQErr(const State&) const;
    Vector& updUErr(const State&) const;
    Vector& updUDotErr(const State&) const;
    Vector& updMultipliers(const State&) const;
    Vector& updEventTriggersByStage(const State&, Stage) const;

    // These pull out the State entries which belong exclusively to
    // this Subsystem. These variables and cache entries are available
    // as soon as this subsystem is at stage Model.
    Stage getStage(const State&) const;
    const AbstractValue& getDiscreteVariable(const State&, DiscreteVariableIndex) const;
    // State is *not* mutable here -- must have write access to change state variables.
    AbstractValue& updDiscreteVariable(State&, DiscreteVariableIndex) const;
    const AbstractValue& getCacheEntry(const State&, CacheEntryIndex) const;
    // State is mutable here.
    AbstractValue& updCacheEntry(const State&, CacheEntryIndex) const;

    bool isCacheValueCurrent(const State&, CacheEntryIndex) const;
    void markCacheValueRealized(const State&, CacheEntryIndex) const;

    // Dimensions. These are valid at System Stage::Model while access to the various
    // arrays may have stricter requirements. Hence it is better to use these
    // routines than to get a reference to a Vector above and ask for its size().

    SystemQIndex getQStart      (const State&) const;
    int getNQ          (const State&) const;
    SystemUIndex getUStart      (const State&) const;
    int getNU          (const State&) const;
    SystemZIndex getZStart      (const State&) const;
    int getNZ          (const State&) const;
    SystemQErrIndex getQErrStart   (const State&) const;
    int getNQErr       (const State&) const;
    SystemUErrIndex getUErrStart   (const State&) const;
    int getNUErr       (const State&) const;
    SystemUDotErrIndex getUDotErrStart(const State&) const;
    int getNUDotErr    (const State&) const;
    SystemMultiplierIndex getMultipliersStart (const State&) const;
    int getNMultipliers     (const State&) const;
    SystemEventTriggerByStageIndex getEventTriggerStartByStage(const State&, Stage) const;
    int getNEventTriggersByStage   (const State&, Stage) const;

	bool isInSystem() const;
	bool isInSameSystem(const Subsystem& otherSubsystem) const;

	const System& getSystem() const;
	System&       updSystem();

	SubsystemIndex getMySubsystemIndex() const;

    // Is this handle the owner of this rep? This is true if the
    // handle is empty or if its rep points back here.
    bool isOwnerHandle() const;
    bool isEmptyHandle() const;

    // There can be multiple handles on the same Subsystem.
    bool isSameSubsystem(const Subsystem& otherSubsystem) const;

    bool subsystemTopologyHasBeenRealized() const;
    void invalidateSubsystemTopologyCache() const;

    // Add a new Measure to this Subsystem. This method is generally used by Measure
    // constructors to install a newly-constructed Measure into its Subsystem.
    MeasureIndex adoptMeasure(AbstractMeasure&);

    AbstractMeasure getMeasure(MeasureIndex) const;
    template <class T> Measure_<T> getMeasure_(MeasureIndex mx) const
    {   return Measure_<T>::getAs(getMeasure(mx));}

    // dynamic_cast the returned reference to a reference to your concrete Guts
    // class.
    const Subsystem::Guts& getSubsystemGuts() const {assert(guts); return *guts;}
    Subsystem::Guts&       updSubsystemGuts()       {assert(guts); return *guts;}

    // Put new Guts into this *empty* handle and take over ownership.
    // If this handle is already in use, this routine will throw
    // an exception.
    void adoptSubsystemGuts(Subsystem::Guts* g);
    void setSystem(System&, SubsystemIndex);

    explicit Subsystem(Subsystem::Guts* g) : guts(g) { }
    bool hasGuts() const {return guts!=0;}
};


/**
 * This is a concrete Subsystem that is part of every System.  It provides a variety of services
 * for the System, such as maintaining lists of event handlers and reporters, and acting as a
 * source of globally unique event IDs.  To obtain the default subsystem for a System, call
 * getDefaultSubsystem() or updDefaultSubsystem() on it.
 */
class SimTK_SimTKCOMMON_EXPORT DefaultSystemSubsystem : public Subsystem {
public:
    DefaultSystemSubsystem(System& sys);
    void addEventHandler(ScheduledEventHandler* handler);
    void addEventHandler(TriggeredEventHandler* handler);
    void addEventReporter(ScheduledEventReporter* handler) const;
    void addEventReporter(TriggeredEventReporter* handler) const;
    EventId createEventId(SubsystemIndex subsys, const State& state) const;
    void findSubsystemEventIds(SubsystemIndex subsys, const State& state, const std::vector<EventId>& allEvents, std::vector<EventId>& eventsForSubsystem) const;
private:
    const DefaultSystemSubsystemGuts& getGuts() const;
    DefaultSystemSubsystemGuts& updGuts();
};


} // namespace SimTK

#endif // SimTK_SimTKCOMMON_SUBSYSTEM_H_
