#ifndef SimTK_SimTKCOMMON_SUBSYSTEM_GUTS_H_
#define SimTK_SimTKCOMMON_SUBSYSTEM_GUTS_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-10 Stanford University and the Authors.        *
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

#include <cassert>

namespace SimTK {

class System;
class DecorativeGeometry;


/**
 * The abstract parent of all Subsystem "Guts" implementation classes.
 */
class SimTK_SimTKCOMMON_EXPORT Subsystem::Guts {
public:
    Guts(const Guts&);
    Guts& operator=(const Guts&);

    // This constructor is for use by concrete Subsystems. Note that this
    // serves as a default constructor since both arguments have defaults.
    explicit Guts(const String& name="<NONAME>", 
                  const String& version="0.0.0");
    virtual ~Guts();

    const String& getName()    const;
    const String& getVersion() const;

    // Use these to allocate state variables and cache entries that are owned
    // by this Subsystem.

    // qdot, qdotdot also allocated in cache
    QIndex allocateQ(State& s, const Vector& qInit) const;
    // udot is also allocated in the cache
    UIndex allocateU(State& s, const Vector& uInit) const;
    // zdot is also allocated in the cache
    ZIndex allocateZ(State& s, const Vector& zInit) const;

    DiscreteVariableIndex allocateDiscreteVariable(State& s, Stage g, AbstractValue* v) const;
    DiscreteVariableIndex allocateAutoUpdateDiscreteVariable
       (State&, Stage invalidates, AbstractValue* v, Stage updateDependsOn) const; 

    // Cache entries
    CacheEntryIndex allocateCacheEntry
       (const State&, Stage dependsOn, Stage computedBy, AbstractValue* v) const;
    CacheEntryIndex allocateCacheEntry
       (const State& state, Stage g, AbstractValue* v) const 
    {   return allocateCacheEntry(state, g, g, v); }
    CacheEntryIndex allocateLazyCacheEntry   
       (const State& state, Stage earliest, AbstractValue* v) const 
    {   return allocateCacheEntry(state, earliest, Stage::Infinity, v); }

    // qerr, uerr, udoterr are all cache entries, not variables
    // allocating udoterr also allocates matching multipliers
    QErrIndex allocateQErr(const State& s, int nqerr) const;
    UErrIndex allocateUErr(const State& s, int nuerr) const;
    UDotErrIndex allocateUDotErr(const State& s, int nudoterr) const;
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

    Real getDiscreteVarLastUpdateTime(const State& s, DiscreteVariableIndex dx) const
    {   return s.getDiscreteVarLastUpdateTime(getMySubsystemIndex(),dx); }
    CacheEntryIndex getDiscreteVarUpdateIndex(const State& s, DiscreteVariableIndex dx) const
    {   return s.getDiscreteVarUpdateIndex(getMySubsystemIndex(),dx); }
    const AbstractValue& getDiscreteVarUpdateValue(const State& s, DiscreteVariableIndex dx) const
    {   return s.getDiscreteVarUpdateValue(getMySubsystemIndex(),dx); }
    AbstractValue& updDiscreteVarUpdateValue(const State& s, DiscreteVariableIndex dx) const
    {   return s.updDiscreteVarUpdateValue(getMySubsystemIndex(),dx); }
    bool isDiscreteVarUpdateValueRealized(const State& s, DiscreteVariableIndex dx) const
    {   return s.isDiscreteVarUpdateValueRealized(getMySubsystemIndex(),dx); }
    void markDiscreteVarUpdateValueRealized(const State& s, DiscreteVariableIndex dx) const
    {   return s.markDiscreteVarUpdateValueRealized(getMySubsystemIndex(),dx); }

    // State is *not* mutable here -- must have write access to change state variables.
    AbstractValue& updDiscreteVariable(State&, DiscreteVariableIndex) const;
    const AbstractValue& getCacheEntry(const State&, CacheEntryIndex) const;
    // State is mutable here.
    AbstractValue& updCacheEntry(const State&, CacheEntryIndex) const;

    bool isCacheValueRealized(const State&, CacheEntryIndex) const;
    void markCacheValueRealized(const State&, CacheEntryIndex) const;
    void markCacheValueNotRealized(const State&, CacheEntryIndex) const;

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
    SystemMultiplierIndex getMultipliersStart(const State&) const;
    int getNMultipliers(const State&) const;
    SystemEventTriggerByStageIndex   getEventTriggerStartByStage(const State&, Stage) const;
    int getNEventTriggersByStage(const State&, Stage) const;

    MeasureIndex adoptMeasure(AbstractMeasure& m);
    AbstractMeasure getMeasure(MeasureIndex) const;
    template <class T> Measure_<T> getMeasure_(MeasureIndex mx) const
    {   return Measure_<T>::getAs(getMeasure(mx));}

	bool isInSystem() const;
	bool isInSameSystem(const Subsystem& otherSubsystem) const;

	const System& getSystem() const;
	System&       updSystem();

	SubsystemIndex getMySubsystemIndex() const;

    // Internal use only
    const Subsystem& getOwnerSubsystemHandle() const;
    Subsystem& updOwnerSubsystemHandle();
    void setOwnerSubsystemHandle(Subsystem&);
    bool hasOwnerSubsystemHandle() const;

    void setSystem(System&, SubsystemIndex);

    class GutsRep;
    explicit Guts(GutsRep* r) : rep(r) { }
    bool                hasRep() const {return rep!=0;}
    const GutsRep& getRep() const {assert(rep); return *rep;}
    GutsRep&       updRep() const {assert(rep); return *rep;}
	void setRep(GutsRep& r) {assert(!rep); rep = &r;}

    bool subsystemTopologyHasBeenRealized() const;
    void invalidateSubsystemTopologyCache() const;

    // These are wrappers for the virtual methods defined below. They
    // are used to ensure good behavior.

    Subsystem::Guts* clone() const;

    // Realize this subsystem's part of the State from Stage-1 to Stage
    // for the indicated stage. After doing some checking, these routines
    // call the concrete subsystem's corresponding virtual method, and
    // on return they make sure the stage has been properly updated.
    // Note that these will do nothing if the Subsystem stage is already
    // at or greater than the indicated stage.
    void realizeSubsystemTopology    (State&) const;
    void realizeSubsystemModel       (State&) const;
    void realizeSubsystemInstance    (const State&) const;
    void realizeSubsystemTime        (const State&) const;
    void realizeSubsystemPosition    (const State&) const;
    void realizeSubsystemVelocity    (const State&) const;
    void realizeSubsystemDynamics    (const State&) const;
    void realizeSubsystemAcceleration(const State&) const;
    void realizeSubsystemReport      (const State&) const;

    // Calculate weights and tolerances.

    // Given the current values of our own position variables,
    // calculate a "unit weighting" for changes to these position
    // variables. 
    // @param[out] weights must be resizable or already the right size
    // @pre State must be realized to >= Stage::Position (this subsystem).
    void calcQUnitWeights(const State&, Vector& weights) const;

    // Given the current values of our own position variables,
    // calculate a "unit weighting" for changes to our velocity
    // variables.
    // @param[out] weights must be resizable or already the right size
    // @pre State must be realized to >= Stage::Position (this subsystem).
    void calcUUnitWeights(const State&, Vector& weights) const;

    // Given the current values of our own position variables,
    // calculate a "unit weighting" for changes to our auxiliary
    // variables.
    // @param[out] weights must be resizable or already the right size
    // @pre State must be realized to >= Stage::Position (this subsystem).
    void calcZUnitWeights(const State&, Vector& weights) const;

    // Calculate a "unit tolerance" for errors in each of this subsystem's
    // position-level constraints.
    // @param[out] tolerances must be resizable or already the right size
    // @pre State must be realized to >= Stage::Model (this subsystem).
    // @remark Tolerances are expected to be constant during a study; typically
    //         they just reflect the units in which the contraint equations
    //         are calculated, e.g. angles or lengths.
    void calcQErrUnitTolerances(const State&, Vector& tolerances) const;

    // Calculate a "unit tolerance" for errors in each of this subsystem's
    // velocity-level constraints.
    // @param[out] tolerances must be resizable or already the right size
    // @pre State must be realized to >= Stage::Model (this subsystem).
    // @remark Tolerances are expected to be constant during a study; typically
    //         they just reflect the units in which the contraint equations
    //         are calculated, e.g. angles/time or lengths/time.
    void calcUErrUnitTolerances(const State&, Vector& tolerances) const;

    // Generate decorative geometry computable at a specific stage. This will
    // throw an exception if this subsystem's state hasn't already been realized
    // to that stage. Note that the list is not inclusive -- you have to
    // request geometry from each stage to get all of it.
    // The generated geometry will be *appended* to the supplied output vector.
    void calcDecorativeGeometryAndAppend
       (const State&, Stage, Array_<DecorativeGeometry>&) const;
    
    void createScheduledEvent(const State& state, EventId& eventId) const;
    void createTriggeredEvent(const State& state, EventId& eventId, 
                              EventTriggerByStageIndex& triggerFunctionIndex,
                              Stage stage) const;

    // These methods are called by the corresponding methods of System.
    // Each subsystem is responsible for defining its own events, and
    // System then combines the information from them, and dispatches events
    // to the appropriate subsystems for handling when they occur.
    virtual void calcEventTriggerInfo
       (const State&, Array_<EventTriggerInfo>&) const;
    virtual void calcTimeOfNextScheduledEvent
       (const State&, Real& tNextEvent, Array_<EventId>& eventIds, 
        bool includeCurrentTime) const;
    virtual void calcTimeOfNextScheduledReport
       (const State&, Real& tNextEvent, Array_<EventId>& eventIds, 
        bool includeCurrentTime) const;
    virtual void handleEvents
       (State&, Event::Cause, const Array_<EventId>& eventIds,
        Real accuracy, const Vector& yWeights, const Vector& ooConstraintTols,
        Stage& lowestModified, bool& shouldTerminate) const;
    virtual void reportEvents
       (const State&, Event::Cause, const Array_<EventId>& eventIds) const;

protected:
    // These virtual methods should be overridden in concrete Subsystems as
    // necessary. They should never be called directly; instead call the
    // wrapper routines above, which have the same name but without the "Impl"
    // (implementation) at the end.
    
    // The "realize..." wrappers will call the "realize...Impl" methods below
    // only when the current stage for the Subsystem is the one just prior
    // to the stage being realized. For example, realizeSubsystemVelocityImpl()
    // is called by realizeSubsystemVelocity() only when the passed-in State
    // shows this subsystem's stage to be exactly Stage::Position.
    //
    // The default implementations provided here do nothing. That means the
    // wrappers will simply check that the current stage is correct and
    // advance it if necessary.

    // The destructor is already virtual; see above.

    virtual Subsystem::Guts* cloneImpl() const = 0;

    virtual int realizeSubsystemTopologyImpl(State& s) const;
    virtual int realizeSubsystemModelImpl(State& s) const;
    virtual int realizeSubsystemInstanceImpl(const State& s) const;
    virtual int realizeSubsystemTimeImpl(const State& s) const;
    virtual int realizeSubsystemPositionImpl(const State& s) const;
    virtual int realizeSubsystemVelocityImpl(const State& s) const;
    virtual int realizeSubsystemDynamicsImpl(const State& s) const;
    virtual int realizeSubsystemAccelerationImpl(const State& s) const;
    virtual int realizeSubsystemReportImpl(const State& s) const;

    virtual int calcQUnitWeightsImpl(const State& s, Vector& weights) const;
    virtual int calcUUnitWeightsImpl(const State& s, Vector& weights) const;
    virtual int calcZUnitWeightsImpl(const State& s, Vector& weights) const;
    virtual int calcQErrUnitTolerancesImpl(const State& s, Vector& tolerances) const;
    virtual int calcUErrUnitTolerancesImpl(const State& s, Vector& tolerances) const;
    virtual int calcDecorativeGeometryAndAppendImpl
       (const State&, Stage, Array_<DecorativeGeometry>&) const;

    void advanceToStage(const State& s, Stage g) const;

private:
    // this is the only data member in the base class
    GutsRep* rep; // opaque implementation of Subsystem::Guts base class.

friend class GutsRep;
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_SUBSYSTEM_GUTS_H_
