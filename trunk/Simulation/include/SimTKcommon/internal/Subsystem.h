#ifndef SimTK_SimTKCOMMON_SUBSYSTEM_H_
#define SimTK_SimTKCOMMON_SUBSYSTEM_H_

/* Portions copyright (c) 2006-7 Stanford University and Michael Sherman.
 * Contributors:
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"
#include "SimTKcommon/internal/State.h"

namespace SimTK {

class System;
class DecorativeGeometry;

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
 * use by the concrete Subsystem class are protected; anything
 * else is private.
 */
class SimTK_SimTKCOMMON_EXPORT Subsystem {
public:
    Subsystem() : rep(0) { } // an empty handle

    // These are allowed only for empty handles and are here just
    // so we can make std::vectors of Subsystems.
    Subsystem(const Subsystem& empty) {
        SimTK_ASSERT_ALWAYS(empty.isEmptyHandle(),
            "Subsystem const copy constructor only allowed if source is an empty handle.");
        rep = 0;
    }
    Subsystem& operator=(const Subsystem& empty) {
        SimTK_ASSERT_ALWAYS(isEmptyHandle() && empty.isEmptyHandle(),
            "Subsystem const copy assignment only allowed if source and destination are both empty handles.");
        return *this;
    }


    // This constructor is for use by concrete Subsystems. The indicated System adopts
    // this handle's rep; that is, some Subsystem handle inside the System will become
    // the owner of our rep, while the current handle will continue to reference it.
    inline explicit Subsystem(System&, const String& name="<NONAME>", 
                                       const String& version="0.0.0");

    //TODO: is this right? May need a private static function to pass
    //to the library which can find the virtual destructor entry.
    inline virtual ~Subsystem() {librarySideDestruction();}


    const String& getName()    const;
    const String& getVersion() const;


    // If a concrete Subsystem needs a private implementation (which it 
    // should if you hope to achieve binary compatibility), derive that
    // private implementation class from this base class. That permits
    // it to be stored with the Subsystem base class's private implementation,
    // which is called SubsystemRep.
    class PrivateImplementation {
    };

    // If a concrete Subsystem class has a private implementation, store it in
    // the SubsystemRep class by passing in a heap pointer and static functions
    // for clone and destruct. The SubsystemRep will take over ownership of
    // the heap space and manage it with the supplied routines.
    typedef void (*DestructPrivateImplementation)(PrivateImplementation*);
    typedef PrivateImplementation* (*ClonePrivateImplementation)(const PrivateImplementation*);
    void adoptPrivateImplementation(PrivateImplementation*,
                                    ClonePrivateImplementation,
                                    DestructPrivateImplementation);

    // Downcast the returned reference to a reference to your private implementation
    // class. Note that the PrivateImplementation class is not abstract, so this is
    // not a dynamic_cast and there is no automated way to check that you are looking
    // at the right type of object. So be careful!
    const PrivateImplementation& getPrivateImplementation() const;
    PrivateImplementation&       updPrivateImplementation();



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

    // These return writable access to this subsystem's partition in the
    // State pool of continuous variables. These can be called at Stage::Model
    // or higher, and if necesary they invalidate the Position (q), Velocity (u),
    // or Dynamics (z) stage respectively.
    Vector& updQ(State&) const; // invalidates Stage::Position
    Vector& updU(State&) const; // invalidates Stage::Velocity
    Vector& updZ(State&) const; // invalidates Stage::Dynamics

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

    // Dimensions. These are valid at System Stage::Model while access to the various
    // arrays may have stricter requirements. Hence it is better to use these
    // routines than to get a reference to a Vector above and ask for its size().

    int getQStart      (const State&) const;
    int getNQ          (const State&) const;
    int getUStart      (const State&) const;
    int getNU          (const State&) const;
    int getZStart      (const State&) const;
    int getNZ          (const State&) const;
    int getQErrStart   (const State&) const;
    int getNQErr       (const State&) const;
    int getUErrStart   (const State&) const;
    int getNUErr       (const State&) const;
    int getUDotErrStart(const State&) const;
    int getNUDotErr    (const State&) const;

	bool isInSystem() const;
	bool isInSameSystem(const Subsystem& otherSubsystem) const;

	const System& getSystem() const;
	System&       updSystem();

	SubsystemId getMySubsystemId() const;

    // Is this handle the owner of this rep? This is true if the
    // handle is empty or if its rep points back here.
    bool isOwnerHandle() const;
    bool isEmptyHandle() const;

    // There can be multiple handles on the same Subsystem.
    bool isSameSubsystem(const Subsystem& otherSubsystem) const;

    // Internal use only
    explicit Subsystem(class SubsystemRep* r) : rep(r) { }
    bool                hasRep() const {return rep!=0;}
    const SubsystemRep& getRep() const {assert(rep); return *rep;}
    SubsystemRep&       updRep() const {assert(rep); return *rep;}
	void setRep(SubsystemRep& r) {assert(!rep); rep = &r;}

    bool subsystemTopologyCacheHasBeenRealized() const;
    void invalidateSubsystemTopologyCache() const;

    /// @name
    /// Realize this subsystem's part of the State from Stage-1 to Stage
    /// for the indicated stage. After doing some checking, these routines
    /// call the concrete subsystem's corresponding virtual method, and
    /// on return they make sure the stage has been properly updated.
    /// Note that these will do nothing if the Subsystem stage is already
    /// at or greater than the indicated stage.
    //@{
    void realizeSubsystemTopology    (State&) const;
    void realizeSubsystemModel       (State&) const;
    void realizeSubsystemInstance    (const State&) const;
    void realizeSubsystemTime        (const State&) const;
    void realizeSubsystemPosition    (const State&) const;
    void realizeSubsystemVelocity    (const State&) const;
    void realizeSubsystemDynamics    (const State&) const;
    void realizeSubsystemAcceleration(const State&) const;
    void realizeSubsystemReport      (const State&) const;
    //@}

    /// @name
    /// Calculate weights and tolerances.
    //@{
    /// Given the current values of our own position variables,
    /// calculate a "unit weighting" for changes to these position
    /// variables. 
    /// @param[out] weights must be resizable or already the right size
    /// @pre State must be realized to >= Stage::Position (this subsystem).
    void calcQUnitWeights(const State&, Vector& weights) const;

    /// Given the current values of our own position variables,
    /// calculate a "unit weighting" for changes to our velocity
    /// variables.
    /// @param[out] weights must be resizable or already the right size
    /// @pre State must be realized to >= Stage::Position (this subsystem).
    void calcUUnitWeights(const State&, Vector& weights) const;

    /// Given the current values of our own position variables,
    /// calculate a "unit weighting" for changes to our auxiliary
    /// variables.
    /// @param[out] weights must be resizable or already the right size
    /// @pre State must be realized to >= Stage::Position (this subsystem).
    void calcZUnitWeights(const State&, Vector& weights) const;

    /// Calculate a "unit tolerance" for errors in each of this subsystem's
    /// position-level constraints.
    /// @param[out] tolerances must be resizable or already the right size
    /// @pre State must be realized to >= Stage::Model (this subsystem).
    /// @remark Tolerances are expected to be constant during a study; typically
    ///         they just reflect the units in which the contraint equations
    ///         are calculated, e.g. angles or lengths.
    void calcQErrUnitTolerances(const State&, Vector& tolerances) const;

    /// Calculate a "unit tolerance" for errors in each of this subsystem's
    /// velocity-level constraints.
    /// @param[out] tolerances must be resizable or already the right size
    /// @pre State must be realized to >= Stage::Model (this subsystem).
    /// @remark Tolerances are expected to be constant during a study; typically
    ///         they just reflect the units in which the contraint equations
    ///         are calculated, e.g. angles/time or lengths/time.
    void calcUErrUnitTolerances(const State&, Vector& tolerances) const;
    //@}

    // Generate decorative geometry computable at a specific stage. This will
    // throw an exception if this subsystem's state hasn't already been realized
    // to that stage. Note that the list is not inclusive -- you have to
    // request geometry from each stage to get all of it.
    // The generated geometry will be *appended* to the supplied output Array.
    void calcDecorativeGeometryAndAppend(const State&, Stage, Array<DecorativeGeometry>&) const;

    Subsystem* clone() const;

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
       (const State&, Stage, Array<DecorativeGeometry>&) const;
    virtual Subsystem* cloneImpl() const;

    // Use these to allocate state variables and cache entries that are owned
    // by this Subsystem.

    // qdot, qdotdot also allocated in cache
    int allocateQ(State& s, const Vector& qInit) const;
    // udot is also allocated in the cache
    int allocateU(State& s, const Vector& uInit) const;
    // zdot is also allocated in the cache
    int allocateZ(State& s, const Vector& zInit) const;
    // qerr, uerr, udoterr are all cache entries, not variables
    int allocateQErr(State& s, int nqerr) const;
    int allocateUErr(State& s, int nuerr) const;
    int allocateUDotErr(State& s, int nudoterr) const;
    int allocateDiscreteVariable(State& s, Stage g, AbstractValue* v) const;
    int allocateCacheEntry(State& s, Stage g, AbstractValue* v) const;
    void advanceToStage(const State& s, Stage g) const;

    // These pull out the State entries which belong exclusively to
    // this Subsystem. These variables and cache entries are available
    // as soon as this subsystem is at stage Model.
    Stage getStage(const State& s) const;
    const AbstractValue& getDiscreteVariable(const State& s, int index) const;
    // State is *not* mutable here -- must have write access to change state variables.
    AbstractValue& updDiscreteVariable(State& s, int index) const;
    const AbstractValue& getCacheEntry(const State& s, int index) const;
    // State is mutable here.
    AbstractValue& updCacheEntry(const State& s, int index) const;

private:
    class SubsystemRep* rep;
    friend class SubsystemRep;
    friend class System;

    // These typedefs are used internally to manage the binary-compatible
    // handling of the virtual function table.

    typedef int (*RealizeWritableStateImplLocator)(const Subsystem&, State&);
    typedef int (*RealizeConstStateImplLocator)(const Subsystem&, const State&);
    typedef int (*CalcUnitWeightsImplLocator)(const Subsystem&, const State&, Vector& weights);
    typedef int (*CalcDecorativeGeometryAndAppendImplLocator)
       (const Subsystem&, const State&, Stage, Array<DecorativeGeometry>&);
    typedef Subsystem* (*CloneImplLocator)(const Subsystem&);

    void librarySideConstruction(System& sys, const String& name, const String& version);
    void librarySideDestruction();

    void registerRealizeTopologyImpl    (RealizeWritableStateImplLocator);
    void registerRealizeModelImpl       (RealizeWritableStateImplLocator);
    void registerRealizeInstanceImpl    (RealizeConstStateImplLocator);
    void registerRealizeTimeImpl        (RealizeConstStateImplLocator);
    void registerRealizePositionImpl    (RealizeConstStateImplLocator);
    void registerRealizeVelocityImpl    (RealizeConstStateImplLocator);
    void registerRealizeDynamicsImpl    (RealizeConstStateImplLocator);
    void registerRealizeAccelerationImpl(RealizeConstStateImplLocator);
    void registerRealizeReportImpl      (RealizeConstStateImplLocator);

    void registerCalcQUnitWeightsImpl(CalcUnitWeightsImplLocator);
    void registerCalcUUnitWeightsImpl(CalcUnitWeightsImplLocator);
    void registerCalcZUnitWeightsImpl(CalcUnitWeightsImplLocator);
    void registerCalcQErrUnitTolerancesImpl(CalcUnitWeightsImplLocator);
    void registerCalcUErrUnitTolerancesImpl(CalcUnitWeightsImplLocator);
    void registerCalcDecorativeGeometryAndAppendImpl(CalcDecorativeGeometryAndAppendImplLocator);
    void registerCloneImpl(CloneImplLocator);

    // We want the locator functions to have access to the protected "Impl"
    // virtual methods, so we make them friends.

    friend int subsystemRealizeTopologyImplLocator(const Subsystem&, State&);
    friend int subsystemRealizeModelImplLocator(const Subsystem&, State&);
    friend int subsystemRealizeInstanceImplLocator(const Subsystem&, const State&);
    friend int subsystemRealizeTimeImplLocator(const Subsystem&, const State&);
    friend int subsystemRealizePositionImplLocator(const Subsystem&, const State&);
    friend int subsystemRealizeVelocityImplLocator(const Subsystem&, const State&);
    friend int subsystemRealizeDynamicsImplLocator(const Subsystem&, const State&);
    friend int subsystemRealizeAccelerationImplLocator(const Subsystem&, const State&);
    friend int subsystemRealizeReportImplLocator(const Subsystem&, const State&);

    friend int subsystemCalcQUnitWeightsImplLocator(const Subsystem&, const State&, Vector&);
    friend int subsystemCalcUUnitWeightsImplLocator(const Subsystem&, const State&, Vector&);
    friend int subsystemCalcZUnitWeightsImplLocator(const Subsystem&, const State&, Vector&);
    friend int subsystemCalcQErrUnitTolerancesImplLocator(const Subsystem&, const State&, Vector&);
    friend int subsystemCalcUErrUnitTolerancesImplLocator(const Subsystem&, const State&, Vector&);
    friend int subsystemCalcDecorativeGeometryAndAppendImplLocator
                (const Subsystem&, const State&, Stage, Array<DecorativeGeometry>&);
    friend Subsystem* subsystemCloneImplLocator(const Subsystem&);
};


// These are used to supply the client-side virtual function to the library, without
// the client and library having to agree on the layout of the virtual function tables.

static int subsystemRealizeTopologyImplLocator(const Subsystem& sys, State& state)
  { return sys.realizeSubsystemTopologyImpl(state); }
static int subsystemRealizeModelImplLocator(const Subsystem& sys, State& state)
  { return sys.realizeSubsystemModelImpl(state); }
static int subsystemRealizeInstanceImplLocator(const Subsystem& sys, const State& state)
  { return sys.realizeSubsystemInstanceImpl(state); }
static int subsystemRealizeTimeImplLocator(const Subsystem& sys, const State& state)
  { return sys.realizeSubsystemTimeImpl(state); }
static int subsystemRealizePositionImplLocator(const Subsystem& sys, const State& state)
  { return sys.realizeSubsystemPositionImpl(state); }
static int subsystemRealizeVelocityImplLocator(const Subsystem& sys, const State& state)
  { return sys.realizeSubsystemVelocityImpl(state); }
static int subsystemRealizeDynamicsImplLocator(const Subsystem& sys, const State& state)
  { return sys.realizeSubsystemDynamicsImpl(state); }
static int subsystemRealizeAccelerationImplLocator(const Subsystem& sys, const State& state)
  { return sys.realizeSubsystemAccelerationImpl(state); }
static int subsystemRealizeReportImplLocator(const Subsystem& sys, const State& state)
  { return sys.realizeSubsystemReportImpl(state); }

static int subsystemCalcQUnitWeightsImplLocator(const Subsystem& sys, const State& s, Vector& weights)
  { return sys.calcQUnitWeightsImpl(s,weights); }
static int subsystemCalcUUnitWeightsImplLocator(const Subsystem& sys, const State& s, Vector& weights)
  { return sys.calcUUnitWeightsImpl(s,weights); }
static int subsystemCalcZUnitWeightsImplLocator(const Subsystem& sys, const State& s, Vector& weights)
  { return sys.calcZUnitWeightsImpl(s,weights); }
static int subsystemCalcQErrUnitTolerancesImplLocator(const Subsystem& sys, const State& s, Vector& tolerances)
  { return sys.calcQErrUnitTolerancesImpl(s,tolerances); }
static int subsystemCalcUErrUnitTolerancesImplLocator(const Subsystem& sys, const State& s, Vector& tolerances)
  { return sys.calcUErrUnitTolerancesImpl(s,tolerances); }
static int subsystemCalcDecorativeGeometryAndAppendImplLocator
   (const Subsystem& sys, const State& s, Stage g, Array<DecorativeGeometry>& geom)
  { return sys.calcDecorativeGeometryAndAppendImpl(s,g,geom); }
static Subsystem* subsystemCloneImplLocator(const Subsystem& sys)
  { return sys.cloneImpl(); }

// Default constructor must be inline so that it has access to the above static
// functions which are private to the client-side compilation unit in which the
// client-side virtual function table is understood.
inline Subsystem::Subsystem(System& sys, const String& name, const String& version) : rep(0)
{
    librarySideConstruction(sys, name, version);

    // Teach the library code how to call client side virtual functions by
    // calling through the client side compilation unit's private static
    // locator functions.
    registerRealizeTopologyImpl    (subsystemRealizeTopologyImplLocator);
    registerRealizeModelImpl       (subsystemRealizeModelImplLocator);
    registerRealizeInstanceImpl    (subsystemRealizeInstanceImplLocator);
    registerRealizeTimeImpl        (subsystemRealizeTimeImplLocator);
    registerRealizePositionImpl    (subsystemRealizePositionImplLocator);
    registerRealizeVelocityImpl    (subsystemRealizeVelocityImplLocator);
    registerRealizeDynamicsImpl    (subsystemRealizeDynamicsImplLocator);
    registerRealizeAccelerationImpl(subsystemRealizeAccelerationImplLocator);
    registerRealizeReportImpl      (subsystemRealizeReportImplLocator);

    registerCalcQUnitWeightsImpl(subsystemCalcQUnitWeightsImplLocator);
    registerCalcUUnitWeightsImpl(subsystemCalcUUnitWeightsImplLocator);
    registerCalcZUnitWeightsImpl(subsystemCalcZUnitWeightsImplLocator);
    registerCalcQErrUnitTolerancesImpl(subsystemCalcQErrUnitTolerancesImplLocator);
    registerCalcUErrUnitTolerancesImpl(subsystemCalcUErrUnitTolerancesImplLocator);
    registerCalcDecorativeGeometryAndAppendImpl(subsystemCalcDecorativeGeometryAndAppendImplLocator);
    registerCloneImpl(subsystemCloneImplLocator);
}



/**
 * This is a concrete Subsystem used by default as the 0th Subsystem of
 * every System. Feel free to replace it with something useful!
 */
class SimTK_SimTKCOMMON_EXPORT DefaultSystemSubsystem : public Subsystem {
public:
    DefaultSystemSubsystem(System& sys)
      : Subsystem(sys, "DefaultSystemSubsystem")
    {
    }
    DefaultSystemSubsystem(System& sys, const String& sysName, const String& sysVersion)
        : Subsystem(sys, sysName, sysVersion)
    {
    }

    /*virtual*/ DefaultSystemSubsystem* cloneImpl() const {
        return new DefaultSystemSubsystem(*this);
    }

    SimTK_DOWNCAST(DefaultSystemSubsystem, Subsystem);
};


} // namespace SimTK

#endif // SimTK_SimTKCOMMON_SUBSYSTEM_H_
