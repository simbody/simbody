#ifndef SimTK_SimTKCOMMON_MEASURE_H_
#define SimTK_SimTKCOMMON_MEASURE_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-13 Stanford University and the Authors.        *
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

/** @file
 *
 * This file declares the base class AbstractMeasure for all derived Measure
 * handle classes, and the handle classes for built-in Measures. Measure
 * handles provide the end user API, while the implementations of Measures
 * derive from the abstract Measure::Implementation class defined in
 * MeasureImplementation.h. Measure Implementation classes provide the Measure
 * developer's API.
 */

#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"

#include <cassert>

/**
 * Every concrete measure handle class "MH" derived from a parent handle "PH"
 * that derives directly or indirectly from the AbstractMeasure handle class
 * should include the macro SimTK_MEASURE_HANDLE_PREAMBLE(MH,PH) at the
 * beginning of the public section of its declaration. It performs the following
 * declarations:
 * <pre>
 *   MH::Implementation         the handle's local implementation class
 *   MH::MH()                   default constructor creates an empty handle
 *   MH::MH(Implementation*)    create a handle referencing an existing object
 *   MH::MH(Subsystem&, Implementation*)
 *                              create a handle referencing an existing but
 *                                unowned object, then installs that in the
 *                                given Subsystem which becomes the owner
 * </pre>
 *
 * MH::Implementation must be defined elsewhere as a class derived
 * from Measure::Implementation.
 */

// Helper macro shared by SimTK_MEASURE_HANDLE_PREAMBLE and
// SimTK_MEASURE_HANDLE_PREAMBLE_ABSTRACT.
#define SimTK_MEASURE_HANDLE_PREAMBLE_BASE(MH,PH) \
    class Implementation;                                           \
    explicit MH(Implementation* imp) : PH(imp) {}                   \
    MH(SimTK::Subsystem& sub, Implementation* imp,                  \
       const SimTK::AbstractMeasure::SetHandle& sh)                 \
    :   PH(sub,imp,sh) {}                                           \
    MH& operator=(const MH& src) {PH::operator=(src); return *this;}\
    MH& shallowAssign(const MH& src) {PH::shallowAssign(src); return *this;}\
    MH& deepAssign(const MH& src) {PH::deepAssign(src); return *this;}


// The default constructor for concrete classes should instantiate
// a default-constructed Implementation object if no Implementation object
// is provided.
#define SimTK_MEASURE_HANDLE_PREAMBLE(MH,PH)    \
    SimTK_MEASURE_HANDLE_PREAMBLE_BASE(MH,PH)   \
    MH() : PH(new Implementation()) {}          \
    explicit MH(SimTK::Subsystem& sub)          \
    : PH(sub,new Implementation(), typename PH::SetHandle()) {}



// The default constructor for a still-abstract derived class can't
// instantiate an Implementation.
#define SimTK_MEASURE_HANDLE_PREAMBLE_ABSTRACT(MH,PH)   \
    SimTK_MEASURE_HANDLE_PREAMBLE_BASE(MH,PH)           \
    MH() : PH() {}

/**
 * Every measure handle class "MH" derived directly or indirectly from the
 * abstract measure handle class "Measure" should include this macro at
 * the end of the "public" section of its declaration. The macro expects
 * there to be a local class, MH::Implementation,
 * already declared. (MH::Implementation is the type of MH's
 * implementation object to be derived from Measure::Implementation
 * and defined elsewhere.) Then the following
 * type-safe downcast methods will be added to MH's definition:
 * <pre>
 *   MH::getAs(const AbstractMeasure&)  generic handle to const MH (static)
 *   MH::updAs(AbstractMeasure&)        generic handle to writable MH (static)
 *   MH::isA(AbstractMeasure&)          test if generic handle is of type MH
 *   getImpl()  generic implementation to const MH::Implementation
 *   updImpl()  generic implementation to writable MH::Implementation
 * </pre>
 * Type checking for the handle class conversions is done only in Debug
 * builds.
 */
#define SimTK_MEASURE_HANDLE_POSTSCRIPT(MH,PH) \
    static bool isA(const SimTK::AbstractMeasure& m)                        \
    {   return dynamic_cast<const Implementation*>(&m.getImpl()) != 0; }    \
    static const MH& getAs(const SimTK::AbstractMeasure& m)                 \
    { assert(isA(m)); return static_cast<const MH&>(m); }                 \
    static MH& updAs(SimTK::AbstractMeasure& m)                             \
    {   assert(isA(m)); return static_cast<MH&>(m); }                       \
    const Implementation& getImpl() const                                   \
    {   return SimTK_DYNAMIC_CAST_DEBUG<const Implementation&>              \
                    (SimTK::AbstractMeasure::getImpl());}                   \
    Implementation& updImpl()                                               \
    {   return SimTK_DYNAMIC_CAST_DEBUG<Implementation&>                    \
                    (SimTK::AbstractMeasure::updImpl());}

namespace SimTK {

class State;
class Subsystem;
class System;
class EventId;

/// Define a unique integral type for safe indexing of Measures.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(MeasureIndex);

//==============================================================================
//                            ABSTRACT MEASURE
//==============================================================================
/** This is the base class for all Measure handle classes. This class is not
templatized, and represents a Measure generically without knowledge of its
value type. This is useful for most of the basic operations that are
performed on measures, such as realization, adoption by a Subsystem, and
other bookkeeping tasks. For most user purposes, the still-generic derived
class Measure_<T> is a more useful handle since its value, of known type T,
can be obtained. All the built-in concrete Measure types derive from
Measure_<T> rather than AbstractMeasure. The various concrete Measures,
built in or otherwise, may set restrictions on the kinds of types that
are allowable as the template argument.

Note that handles always consist of exactly one pointer, and handle classes
are always concrete (meaning they have no virtual methods). **/
class SimTK_SimTKCOMMON_EXPORT AbstractMeasure {
protected:
    /// An object of this type is used as a dummy argument to make sure the
    /// automatically-generated handle constructor's signature doesn't conflict
    /// with an explicitly-defined one.
    class SetHandle {};

public:
    class Implementation; // local; name is AbstractMeasure::Implementation

    /// Provide an Implementation for this AbstractMeasure and bump its
    /// reference count. This is also the default constructor for the
    /// base class producing an empty handle.
    explicit AbstractMeasure(Implementation* g=0);

    /// Construct this handle with a given Implementation object (whose
    /// reference count will be bumped) and then let the given Subsystem
    /// adopt this Measure (which will again bump the Implementation's
    /// reference count, leaving us with two new handles).
    AbstractMeasure(Subsystem&, Implementation* g, const SetHandle&);

    /// Shallow copy constructor copies the pointer from the source
    /// Implementation object and bumps its reference count.
    AbstractMeasure(const AbstractMeasure&);

    /// Shallow assignment operator results in this handle referencing
    /// the same Implementation object as does the source.
    /// @see shallowAssign for more details
    AbstractMeasure& operator=(const AbstractMeasure& source)
    {   return shallowAssign(source); }

    /// Destructor decrements the Implementation's reference count and
    /// deletes the object if the count goes to zero.
    ~AbstractMeasure();

    /// Shallow assignment operator destructs the current Implementation
    /// object (meaning its reference count is decremented and the object
    /// actually deleted only if the count goes to zero), then copies the
    /// Implementation pointer from the source and bumps its reference count.
    /// This is what the copy assignment operator= does for this class.
    AbstractMeasure& shallowAssign(const AbstractMeasure&);

    /// Deep assignment clones the Implementation object pointed to by
    /// the source handle, so that this handle ends up pointing to a
    /// new Measure object similar to the original but not yet contained
    /// in any Subsystem.
    AbstractMeasure& deepAssign(const AbstractMeasure& source);

    /// Every Measure can produce a value, and some can provide one or
    /// more total derivatives with respect to time of that value. This
    /// method reports how many are available: 1 -> first derivative
    /// d/dt is available, 2 -> first and second derivative d^2/dt^2
    /// are available, etc. We interpret the "0th derivative" to be the
    /// Measure's value.
    /// @return The maximum available derivative order.
    int getNumTimeDerivatives() const;

    /// At what Stage can we expect the value of this AbstractMeasure or
    /// one of its time derivatives to be available? Users of Measures will
    /// typically impose restrictions on the levels they will accept.
    /// @param[in] derivOrder
    ///     Which derivative level is to be checked: 0 -> the value,
    ///     1 -> the 1st time derivative, etc. Must not be higher than the
    ///     value returned by getNumTimeDerivatives().
    /// @return The Stage after which this value is available.
    Stage getDependsOnStage(int derivOrder=0) const;


    /// There can be multiple handles on the same Measure.
    bool isSameMeasure(const AbstractMeasure& other) const
    {   return impl && impl==other.impl;}

    bool isEmptyHandle() const {return !hasImpl();}

    /// Test whether this Measure object has been adopted by a Subsystem.
    bool isInSubsystem() const;
    /// Return a reference to the Subsystem that owns this Measure. Will
    /// throw an exception if the Measure is not currently owned by any
    /// Subsystem.
    const Subsystem& getSubsystem()  const;
    /// Return the MeasureIndex by which this Measure is known to the Subsystem
    /// that owns it. Will throw an exception if the Measure is not currently
    /// owned by any Subsystem.
    MeasureIndex getSubsystemMeasureIndex() const;

    // Internal use only

    // dynamic_cast the returned reference to a reference to your concrete
    // Implementation class.
    const Implementation& getImpl() const {assert(impl); return *impl;}
    Implementation&       updImpl()       {assert(impl); return *impl;}
    bool                  hasImpl() const {return impl!=0;}

    int getRefCount() const;
private:
    // This is the only data member in this class. Also, any class derived
    // from AbstractMeasure must have *NO* data members at all (data goes
    // in the Implementation class).
    Implementation* impl;

friend class Implementation;
};


//==============================================================================
//                               MEASURE <T>
//==============================================================================
/** This is the base handle class for all Measures whose value type is known,
including all the Simbody built-in %Measure types. **/
template <class T>
class Measure_ : public AbstractMeasure {
public:
    /** This class is still abstract so we don't want it to allocate an
    Implementation object in its default constructor. **/
    SimTK_MEASURE_HANDLE_PREAMBLE_ABSTRACT(Measure_, AbstractMeasure);

    /** Retrieve the Value of this Measure or one of its time derivatives,
    assuming the supplied State has been realized to at least the required
    stage for the selected value or derivative, as reported by
    getDependsOnStage(). If the stage is sufficient but the corresponding
    value has not yet been computed, it will be computed first with its
    value going into this State's cache so that subsequent calls do not
    require further computation. **/
    const T& getValue(const State& s, int derivOrder=0) const
    {   return getImpl().getValue(s,derivOrder); }

    /** Change the default value associated with this %Measure. How this is
    used varies with the %Measure type but generally it is the value that
    the %Measure will have before any calculations are performed. Note
    that this does not require a State since it is a Topology-stage
    change; you have to call realizeTopology() again if you call this
    method. **/
    Measure_& setDefaultValue(const T& defaultValue)
    {   updImpl().setDefaultValue(defaultValue); return *this; }

    /** Obtain a reference to the default value associated with
    this %Measure. **/
    const T& getDefaultValue() const
    {   return getImpl().getDefaultValue(); }

    // These are built-in Measures with local class names.

    // Templatized measures may have restrictions on the allowable template
    // type and may be specialized for particular types.
    class Zero;         // T is any numerical type
    class One;          // T is any numerical type
    class Constant;     // T is any assignable type
    class Time;         // T is any type for which T(t) makes sense.
    class Variable;     // T is any assignable type (state)
    class Result;       // T is any assignable type (cache)
    class SampleAndHold;// T is any assignable type
    class Delay;        // T is any assignable type

    // This requires any numerical type.
    class Plus;
    class Minus;
    class Scale;
    class Differentiate;

    // These find extreme values *in time*, not among inputs at the same
    // time. They perform elementwise on aggregate types.
    class Extreme;  // base class for min/max/minabs/maxabs
    class Minimum;  // most positive value
    class Maximum;  // most negative value
    class MinAbs;   // the signed quantity whose absolute value was min
    class MaxAbs;   // the signed quantity whose absolute value was max

    // These accept floating point numerical template arguments only.
    class Integrate;
    class Sinusoid;

    SimTK_MEASURE_HANDLE_POSTSCRIPT(Measure_, AbstractMeasure);
};

/** This typedef is a convenient abbreviation for the most common kind of
%Measure -- one that returns a single Real result; the underlying class is
Measure_; look there for documentation. **/
typedef Measure_<Real> Measure;


//==============================================================================
//                                CONSTANT
//==============================================================================
/** This creates a Measure whose value is a Topology-stage constant of any
type T. This does not have to be part of a Subsystem, but if it is then changing
the constant's value invalidates the containing Subsystem's Topology.
@tparam T   This can be any type that supports copy construction. **/
template <class T>
class Measure_<T>::Constant : public Measure_<T> {
public:
    SimTK_MEASURE_HANDLE_PREAMBLE(Constant, Measure_<T>);

    /** Create a constant measure that is not part of any Subsystem, and
    provide the constant value. **/
    explicit Constant(const T& value)
    :   Measure_<T>(new Implementation(value)) {}

    /** Create a constant measure with the given \a value and install it into
    the given Subsystem. **/
    Constant(Subsystem& sub, const T& value)
    :   Measure_<T>(sub, new Implementation(value),
                    AbstractMeasure::SetHandle()) {}

    /** Change the value returned by this %Measure. Note that this does not
    require a State since it is a Topology-stage change. **/
    Constant& setValue(const T& value)
    {   updImpl().setValue(value); return *this; }

    SimTK_MEASURE_HANDLE_POSTSCRIPT(Constant, Measure_<T>);
};

//==============================================================================
//                                   ZERO
//==============================================================================
/** This creates a Measure::Constant whose value is always T(0) and can't be
changed. This class is specialized for Vector so that you must provide a
size at construction. **/
template <class T>
class Measure_<T>::Zero : public Measure_<T>::Constant {
public:
    Zero();
    explicit Zero(Subsystem& sub);
};

template <>
class Measure_< Vector >::Zero : public Measure_< Vector >::Constant {
public:
    explicit Zero(int size);
    Zero(Subsystem& sub, int size);
};

//==============================================================================
//                                    ONE
//==============================================================================
/** This creates a Measure::Constant whose value is always T(1) and can't be
changed. This class is specialized for Vector so that you must provide a
size at construction. **/
template <class T>
class Measure_<T>::One : public Measure_<T>::Constant {
public:
    One();
    explicit One(Subsystem& sub);
};

template <>
class Measure_< Vector >::One : public Measure_< Vector >::Constant {
public:
    explicit One(int size);
    One(Subsystem& sub, int size);
};

//==============================================================================
//                                   TIME
//==============================================================================
/** This creates a Measure::Time whose value is always T(time). **/
template <class T>
class Measure_<T>::Time : public Measure_<T> {
public:
    SimTK_MEASURE_HANDLE_PREAMBLE(Time, Measure_<T>);

    SimTK_MEASURE_HANDLE_POSTSCRIPT(Time, Measure_<T>);
};

//==============================================================================
//                                 VARIABLE
//==============================================================================
/** This creates a Measure whose value is a discrete State variable of any
type T. **/
template <class T>
class Measure_<T>::Variable : public Measure_<T> {
public:
    SimTK_MEASURE_HANDLE_PREAMBLE(Variable, Measure_<T>);

    // TODO: should not require invalidated Stage here. Instead,
    // should have a unique "generation" counter for this variable
    // and allow subsequent users to check it.
    Variable(Subsystem& sub, Stage invalidates, const T& defaultValue)
    :   Measure_<T>(sub, new Implementation(invalidates, defaultValue),
                    AbstractMeasure::SetHandle()) {}


    void setValue(State& state, const T& value) const
    {   getImpl().setValue(state, value); }

    SimTK_MEASURE_HANDLE_POSTSCRIPT(Variable, Measure_<T>);
};

//==============================================================================
//                                 RESULT
//==============================================================================
/** This Measure holds the result of some externally-determined computation,
and helps to coordinate the validity of that computation with respect
to the state variables. The value must be set manually and explicitly
marked valid when it is complete. The value will be automatically
invalidated after a state change at or below a specified Stage, and changing
the value here will automatically invalidate a specified Stage and above.

In constrast to Measure::Variable, Measure::Result is not a state variable;
it is a cache variable meaning that it works with a const State. It is
expected that the result can be recalculated from the state variables when
needed, or contains ephemeral information that can be discarded.

No provision for derivatives is made; there is only the one result. **/
template <class T>
class Measure_<T>::Result : public Measure_<T> {
public:
    SimTK_MEASURE_HANDLE_PREAMBLE(Result, Measure_<T>);

    // TODO: should not require invalidated Stage here. Instead,
    // should have a unique "generation" counter for this cache entry
    // and allow subsequent users of the value to check it.

    /// Create a new Result measure and add it to the indicated
    /// subsystem. The Result measure's value depends on earlier
    /// values calculated at \a dependsOn stage, so whenever that
    /// stage or earlier is invalidated, so is the value here.
    /// In addition, any change made to the value here will invalidate
    /// stage \a invalidated and all subsequent stages. The
    /// \a invalidated stage must be later than the \a dependsOn stage,
    /// and you cannot update the value here in a state that hasn't
    /// already been realized to the dependsOn stage.
    /// Set \a dependsOn = Stage::Empty if your value computation doesn't
    /// depend on anything else (that will be interpreted as Stage::Topology,
    /// however, since space for the value gets allocated then).
    /// Set \a invalidated = Stage::Infinity if you don't want anything
    /// invalidated when this value is changed.
    Result(Subsystem& sub, Stage dependsOn, Stage invalidated)
    :   Measure_<T>(sub, new Implementation(dependsOn, invalidated),
                    AbstractMeasure::SetHandle()) {}

    /// Get the \a dependsOn stage for this measure's value.
    Stage getDependsOnStage() const {return getImpl().getDependsOnStage();}
    /// Get the \a invalidated stage for this measure's value.
    Stage getInvalidatedStage() const {return getImpl().getInvalidatedStage();}
    /// Change the \a dependsOn stage for this measure's value, which must
    /// be strictly less than the current setting for the \a invalidated
    /// stage. If you set the dependsOn stage to Stage::Empty it will be
    /// interpreted as Stage::Topology since the value must always depend
    /// on at least topology. Setting the dependsOn stage is itself a
    /// topological change requiring reallocation of the value if it has
    /// already been allocated; you must call realizeTopology() again.
    Result& setDependsOnStage(Stage dependsOn)
    {   updImpl().setDependsOnStage(dependsOn); return *this; }
    /// Change the \a invalidated stage for this measure's value, which must
    /// be strictly greater than the current setting for the \a dependsOn
    /// stage. This is a topological change requiring reallocation of the
    /// value if it has already been allocated; you must call
    /// realizeTopology() again.
    Result& setInvalidatedStage(Stage invalidated)
    {   updImpl().setInvalidatedStage(invalidated); return *this; }

    /// Normally a Result measure's value is not considered valid unless
    /// we are notified explicitly that it is, via a call to markAsValid()
    /// or setValue(); this method allows the value to be assumed valid
    /// after the \a dependsOn stage has been realized. The reason this
    /// exists is that in some cases it is difficult to find a place from
    /// which to call markAsValid(). That means you must set the value during
    /// realization of the dependsOn stage, but there is no way for this to
    /// be checked automatically. Thus use of this feature can lead
    /// to very difficult-to-find bugs; you should try hard to find a place
    /// to call markAsValid() before resorting to this method. This is a
    /// topological change requiring reallocation of the value if it has
    /// already been allocated; you must call realizeTopology() again.
    Result& setIsPresumedValidAtDependsOnStage(bool presume)
    {   updImpl().setIsPresumedValidAtDependsOnStage(presume); return *this; }

    /// Return the value of the "presumed valid at dependsOn stage" flag.
    /// @see setIsPresumedValidAtDependsOnStage() for a discussion.
    bool getIsPresumedValidAtDependsOnStage() const
    {   return getImpl().getIsPresumedValidAtDependsOnStage(); }


    /// Obtain write access to the Measure's value in order to modify it.
    /// Calling this method marks the result invalid; you must explicitly
    /// validate it when you're done. Also, if this Measure was created with
    /// an \a invalidates stage, then that stage and all later stages in
    /// the given state are marked invalid also.
    T& updValue(const State& state) const
    {   return getImpl().updValue(state); }

    /// Mark the current value as valid. This is done automatically if you
    /// call setValue() but must be done manually if you use updValue() to
    /// access the value. Note that you cannot mark this valid if
    /// \a state hasn't been realized already to at least the stage prior
    /// to this measure's \a dependsOn stage; that is, you must at least
    /// be working on the dependsOn stage at the time this is called.
    void markAsValid(const State& state) const {getImpl().markAsValid(state);}

    /// Check whether the value contained in this Measure is currently
    /// valid. If true, you can call getValue() to obtain it. If false,
    /// calling getValue() will throw an exception.
    bool isValid(const State& state) const {
      return getImpl().isValid(state);}

    /// Manually mark the contained value as invalid. This will also
    /// invalidate any stages in \a state that depend on this value. This
    /// is called automatically whenever the updValue() method is invoked,
    /// and the value starts out invalid. The value becomes valid either
    /// by calling markAsValid() explicitly or calling setValue().
    /// @warning If you have set this Result measure to be presumed valid
    /// at dependsOn stage, this method will have no effect.
    void markAsNotValid(const State& state) const
    {   getImpl().markAsNotValid(state); }

    /// Set a new value and mark it as valid. For more flexibility, use the
    /// updValue() method and markAsValid() manually when you are done with
    /// the new value.
    void setValue(const State& state, const T& value) const
    {   updValue(state) = value; markAsValid(state); }

    SimTK_MEASURE_HANDLE_POSTSCRIPT(Result, Measure_<T>);
};

//==============================================================================
//                                 SINUSOID
//==============================================================================
/** This measure produces a sinusoidal function of time:
<pre>
     m(t) = a*sin(w*t+p)
</pre>
where a=amplitude in arbitrary units, w=frequency in rad/unit time, p=phase
in radians. **/
template <class T>
class Measure_<T>::Sinusoid : public Measure_<T> {
public:
    SimTK_MEASURE_HANDLE_PREAMBLE(Sinusoid, Measure_<T>);

    Sinusoid(Subsystem& sub,
             const T& amplitude,
             const T& frequency,
             const T& phase=T(0))
    :   Measure_<T>(sub, new Implementation(amplitude,frequency,phase),
                    AbstractMeasure::SetHandle()) {}

    SimTK_MEASURE_HANDLE_POSTSCRIPT(Sinusoid, Measure_<T>);
};

//==============================================================================
//                                   PLUS
//==============================================================================
/** This Measure is the sum of two Measures of the same type T.
@tparam T    Any type that supports a plus operator that returns a sum
             as another object of type T. In particular, Real, Vec<N>,
             and Vector will work. **/
template <class T>
class Measure_<T>::Plus : public Measure_<T> {
public:
    SimTK_MEASURE_HANDLE_PREAMBLE(Plus, Measure_<T>);

    Plus(Subsystem& sub, const Measure_<T>& left, const Measure_<T>& right)
    :   Measure_<T>(sub, new Implementation(left, right),
                    AbstractMeasure::SetHandle())
    {   SimTK_ERRCHK_ALWAYS
           (   this->getSubsystem().isSameSubsystem(left.getSubsystem())
            && this->getSubsystem().isSameSubsystem(right.getSubsystem()),
            "Measure_<T>::Plus::ctor()",
            "Arguments must be in the same Subsystem as this Measure.");
    }

    SimTK_MEASURE_HANDLE_POSTSCRIPT(Plus, Measure_<T>);
};

//==============================================================================
//                                   MINUS
//==============================================================================
/** This Measure is the difference of two Measures of the same type T.
@tparam T    Any type that supports a subtract operator that returns the
             difference as another object of type T. In particular, Real,
             Vec<N>, and Vector will work. **/
template <class T>
class Measure_<T>::Minus : public Measure_<T> {
public:
    SimTK_MEASURE_HANDLE_PREAMBLE(Minus, Measure_<T>);

    Minus(Subsystem& sub, const Measure_<T>& left, const Measure_<T>& right)
    :   Measure_<T>(sub, new Implementation(left, right),
                    AbstractMeasure::SetHandle())
    {   SimTK_ERRCHK_ALWAYS
           (   this->getSubsystem().isSameSubsystem(left.getSubsystem())
            && this->getSubsystem().isSameSubsystem(right.getSubsystem()),
            "Measure_<T>::Minus::ctor()",
            "Arguments must be in the same Subsystem as this Measure.");
    }

    SimTK_MEASURE_HANDLE_POSTSCRIPT(Minus, Measure_<T>);
};

//==============================================================================
//                                   SCALE
//==============================================================================
/** This Measure multiplies some other Measure by a Real scale factor.
@tparam T    Any type that supports a scalar multiply, with the scalar on
             the left, that returns the product as another object of type T.
             In particular, Real, Vec<N>, and Vector will work. **/
template <class T>
class Measure_<T>::Scale : public Measure_<T> {
public:
    SimTK_MEASURE_HANDLE_PREAMBLE(Scale, Measure_<T>);

    Scale(Subsystem& sub, Real factor, const Measure_<T>& operand)
    :   Measure_<T>(sub, new Implementation(factor, operand),
                    AbstractMeasure::SetHandle())
    {   SimTK_ERRCHK_ALWAYS
           (this->getSubsystem().isSameSubsystem(operand.getSubsystem()),
            "Measure_<T>::Scale::ctor()",
            "Argument must be in the same Subsystem as this Measure.");
    }

    /** Get the operand (thing being scaled) measure for this measure. **/
    const Measure_<T>& getOperandMeasure() const
    { return getImpl().getOperandMeasure(); }

    SimTK_MEASURE_HANDLE_POSTSCRIPT(Scale, Measure_<T>);
};

//==============================================================================
//                                 INTEGRATE
//==============================================================================
/** This measure yields the time integral of a given derivative measure,
initializing with an initial condition measure of the same type T. The type
T can be Real or a fixed size Vec type like Vec<3>, or a variable-size
Vector. But in the case of a Vector you must say at construction what size
the Vector will be during the simulation, so that an appropriate number
of state variables can be allocated. **/
template <class T>
class Measure_<T>::Integrate : public Measure_<T> {
public:
    SimTK_MEASURE_HANDLE_PREAMBLE(Integrate, Measure_<T>);

    /** Create a new measure that will use %Measure \a ic's value for initial
    conditions, and then integrate the \a deriv %Measure. If type T is a
    variable-length Vector you must also provide a const Vector to be used
    to size and initialized the number of state variables at allocation. In
    that case both the \a ic and \a deriv measures must return Vector values
    of the same size as the allocation constant. **/
    Integrate(Subsystem&            subsystem,
              const Measure_<T>&    deriv,
              const Measure_<T>&    ic,
              const T&              initAlloc=T(0))
    :   Measure_<T>(subsystem, new Implementation(deriv,ic,initAlloc),
                    AbstractMeasure::SetHandle()) {}

    /** Set the current value of this measure by modifying the state variables
    that hold the integral. **/
    void setValue(State& s, const T& value) const
    {   return getImpl().setValue(s, value); }

    /** Get the integrand (derivative) measure for this integral. **/
    const Measure_<T>& getDerivativeMeasure() const

    {   return getImpl().getDerivativeMeasure(); }
    /** Get the measure whose value is used as an initial condition for the
    integral at the start of an integration. **/
    const Measure_<T>& getInitialConditionMeasure() const
    {   return getImpl().getInitialConditionMeasure(); }

    Integrate& setDerivativeMeasure(const Measure_<T>& d)
    {   updImpl().setDerivativeMeasure(d); return *this; }
    Integrate& setInitialConditionMeasure(const Measure_<T>& ic)
    {   updImpl().setInitialConditionMeasure(ic); return *this; }

    SimTK_MEASURE_HANDLE_POSTSCRIPT(Integrate, Measure_<T>);
};

//==============================================================================
//                               DIFFERENTIATE
//==============================================================================
/** This Measure operator returns the time derivative of its operand measure,
or a numerical approximation of the time derivative if an analytic one is not
available.

If the operand measure provides its own derivative measure, then the value of
the Differentiate operator is just the value of the operand's derivative
measure, and this measure will have one fewer available derivatives than does
the operand. If the operand does not have a derivative, then we will estimate
it by the following method:
    - retrieve the previous value f0 and previous derivative fdot0 of the
      operand measure, and their sample time t0
    - obtain the current value f(t) of the operand
    - estimate fdot(t)=2(f-f0)/(t-t0) - fdot0  (fit a quadratic)
    - record new samples f(t), fdot(t) with timestamp t

Special cases:
    - if t==t0 then fdot(t)=fdot0 (if available) else fdot(t)=0
    - if fdot0 not available, fdot(t)=(f-f0)/(t-t0) (first order estimate)

At initialization of a timestepping study beginning at t=t0, we sample the
operand and record its initial value f0 at t0, and set fdot0=NaN. This
ensures that we'll return zero as the initial derivative (for lack of anything
better) and then use the first order method for the first step's derivative.
**/
template <class T>
class Measure_<T>::Differentiate : public Measure_<T> {
public:
    SimTK_MEASURE_HANDLE_PREAMBLE(Differentiate, Measure_<T>);

    /** Create a measure whose value is the time derivative of the given
    \a operand measure.
    @param  subsystem   The Subsystem into which this measure will be placed.
    @param  operand     The Measure to be differentiated. **/
    Differentiate(Subsystem& subsystem, const Measure_<T>& operand)
    :   Measure_<T>(subsystem, new Implementation(operand),
                    AbstractMeasure::SetHandle()) {}

    /** Test whether the derivative returned as the value of this measure is
    being estimated numerically, either because the operand measure is unable
    to supply its derivative or because setForceUseApproximation(true) has
    been called. **/
    bool isUsingApproximation() const
    {   return getImpl().isUsingApproximation(); }

    /** Get a reference to the measure that is being differentiated by this
    measure. **/
    const Measure_<T>& getOperandMeasure() const
    {   return getImpl().getOperandMeasure(); }

    /** Set the measure that is to be differentiated by this measure. This is
    a topology-stage change so you'll have to call realizeTopology() again on
    the enclosing System before using it. **/
    Differentiate& setOperandMeasure(const Measure_<T>& operand)
    {   updImpl().setOperandMeasure(operand); return *this; }

    /** Force use of numerical approximation for the derivative, even if the
    operand measure can supply its own derivative. This is not recommended!
    This is a Topology-stage change. **/
    void setForceUseApproximation(bool mustApproximate)
    {   updImpl().setForceUseApproximation(mustApproximate); }

    /** Check the current value of the flag which forces this measure to use
    numerical approximation regardless of whether the operand can supply its
    own derivative. Note that even if the flag is currently false (the default)
    we may still have to use approximation; see isUsingApproximation(). **/
    bool getForceUseApproximation() const
    {   return getImpl().getForceUseApproximation(); }

    SimTK_MEASURE_HANDLE_POSTSCRIPT(Differentiate, Measure_<T>);
};

//==============================================================================
//                 EXTREME, MINIMUM, MAXIMUM, MINABS, MAXABS
//==============================================================================
/** This Measure tracks extreme values attained by the elements of its source
operand since the last initialize() call or explicit call to setValue(). The
extreme is either minimum or maximum and may be determined by the actual or
absolute value of the operand. In any case the value of the %Extreme measure is
the actual extreme value of the operand, not its absolute value.

The template type T must be the same as the template type of the operand
measure. If T is a Vec or Vector type, each element is treated separately.

Normally %Extreme is not used directly; it is the common implementation
underlying the Minimum, Maximum, MinAbs, and MaxAbs measures.

Information available from this Measure:
    - the current extreme value (at the source's stage)
    - the value at the start of the current step (Time stage)
    - time of last change
    - the DiscreteVariableIndex holding the previous value
    - the CacheEntryIndex designated for the current/next value
    - a reference to the source Measure

The time derivative fdot of f(t)=min_t0_t(s(t')) where s is the
source measure and t0 <= t' <= t is
<pre>
     fdot(t) = s(t) < f(ti) && sdot(t) < 0 ? sdot(t) : 0
</pre>
where ti is the time at the start of the current step.

At the start of a continuous interval, the updated value (if any)
replaces the stored value.

<h3>Extreme isolation (not implemented yet)</h3>
If the time derivative of the
source operand is available, the measure will arrange to ensure
precise isolation of the minimum values by defining a triggered
event that watches for negative-to-positive sign changes of the derivative.
Then if you output reporting data upon the occurrence of triggered
events (as well as your regularly scheduled output) your data will
include the precise minimum (to within a specifiable isolation time window).

Additional information available in this case:
    - the EventId of the extreme-trapping event if there is one
    - the time derivative of the Extreme measure
**/
template <class T>
class Measure_<T>::Extreme : public Measure_<T> {
public:
    SimTK_MEASURE_HANDLE_PREAMBLE(Extreme, Measure_<T>);

    enum Operation {
        MaxAbs, // default
        Maximum,
        MinAbs,
        Minimum
    };

    /** Default behavior for the %Extreme measure is to find the operand's
    value that is of maximum absolute value. You can change that to minimum
    and/or actual value. **/
    Extreme(Subsystem& sub, const Measure_<T>& operand, Operation op=MaxAbs)
    :   Measure_<T>(sub, new Implementation(operand, op),
                    AbstractMeasure::SetHandle()) {}

    /** Set the operation to be performed. The default operation is MaxAbs. **/
    Extreme& setOperation(Operation op)
    {   updImpl().setOperation(op); return *this; }

    /** Return the operation currently being performed by this measure. **/
    Operation getOperation() const {return getImpl().getOperation();}

    /** Return the time at which the reported extreme value first occurred.
    This is the current time if the operand is at its extreme value now,
    otherwise it is the time that the extreme value first occurred during a
    time stepping study. The \a state must be realized to the level required
    to evaluate the operand measure. **/
    Real getTimeOfExtremeValue(const State& state) const
    {   return getImpl().getTimeOfExtremeValue(state); }

    void setValue(State& s, const T& value) const
    {   return getImpl().setValue(s, value); }

    const Measure_<T>& getOperandMeasure() const
    {   return getImpl().getOperandMeasure(); }

    Extreme& setOperandMeasure(const Measure_<T>& s)
    {   updImpl().setOperandMeasure(s); return *this; }

    SimTK_MEASURE_HANDLE_POSTSCRIPT(Extreme, Measure_<T>);
};

/** Track the minimum value of the operand (signed).
@see Measure_::Extreme **/
template <class T>
class Measure_<T>::Minimum : public Measure_<T>::Extreme {
    typedef typename Measure_<T>::Extreme Super;
public:
    Minimum(Subsystem& sub, const Measure_<T>& operand)
    :   Super(sub, operand, Super::Minimum) {}
};

/** Track the maximum value of the operand (signed).
@see Measure_::Extreme **/
template <class T>
class Measure_<T>::Maximum : public Measure_<T>::Extreme {
    typedef typename Measure_<T>::Extreme Super;
public:
    Maximum(Subsystem& sub, const Measure_<T>& operand)
    :   Super(sub, operand, Super::Maximum) {}
};

/** Track the value of the operand that is of maximum absolute value.
@see Measure_::Extreme **/
template <class T>
class Measure_<T>::MaxAbs : public Measure_<T>::Extreme {
    typedef typename Measure_<T>::Extreme Super;
public:
    MaxAbs(Subsystem& sub, const Measure_<T>& operand)
    :   Super(sub, operand, Super::MaxAbs) {}
};

/** Track the value of the operand that is of minimum absolute value (not very
useful).
@see Measure_::Extreme **/
template <class T>
class Measure_<T>::MinAbs : public Measure_<T>::Extreme {
    typedef typename Measure_<T>::Extreme Super;
public:
    MinAbs(Subsystem& sub, const Measure_<T>& operand)
    :   Super(sub, operand, Super::MinAbs) {}
};

//==============================================================================
//                                 DELAY
//==============================================================================
/** (CAUTION: still under development)
This is a %Measure whose value at time t is the value that its
\a source operand had at time t-delay for a specified \a delay. For times prior
to the start of a simulation this %Measure behaves as though the source value
had been constant at its initial value.

When the \a source %Measure can provide a time derivative dvalue we use saved
(t,value,dvalue) triples surrounding the required time to construct a cubic
Hermite interpolant giving a third-order accurate estimate of the delayed
value. Otherwise we use more data points to construct the cubic interpolant but
the accuracy cannot be guaranteed. If there aren't enough data points, then
linear interpolation is used instead. There is an option to force use of linear
interpolation if you prefer.

In the case where the delayed time is within the current step, we would need
the current \a source value in order to interpolate. We assume that is not
available (commonly the current value depends on the delayed value) so have to
extrapolate beyond the last buffered value in that case. Extrapolation is
considerably less accurate than interpolation, so when step sizes are large
compared to delay times the accuracy of the delayed value is reduced. In cases
where the \a source does not depend on its delayed value, you can request that
the current value be used if necessary, ensuring consistent accuracy.
Alternatively, you can set the maximum integrator step size to be just less
than the minimum delay time, guaranteeing that there will always be an entry
already in the buffer that is later than any requested delayed time. That
could have a substantial performance penalty if steps much larger than the
delay would otherwise have been taken.

<h3>%Implementation</h3>
This %Measure maintains a variable-sized buffer holding values that the
\a source measure and its time derivative (if available) had at each time step
starting just prior to t-delay until just before the current time t. New values
are added to the end of the buffer as integrator steps are completed, and old
values that are no longer needed are removed from the beginning. When a value is
requested at current time t, the %Measure interpolates using values from just
prior to t-delay and just afterwards to approximate the value at t-delay.

@bug Only linear interpolation implemented so far.
@bug There should be an option for the measure to specify a sampling interval
that would force the integrator to provide interpolated states at least that
often.
@bug The current implementation involves a lot of unnecessary copying because
it uses an auto-update state variable and cache entry. It should be using
an end-of-step event handler to update the buffer directly in the state.

@see Measure_::Integrate, Measure_::Differentiate **/
template <class T>
class Measure_<T>::Delay : public Measure_<T> {
public:
    /** @cond **/
    SimTK_MEASURE_HANDLE_PREAMBLE(Delay, Measure_<T>);
    /** @endcond **/

    /** Create a %Measure whose output is the same as the given \a source
    measure but delayed by a time \a delay. **/
    Delay(Subsystem& sub, const Measure_<T>& source, Real delay)
    :   Measure_<T>(sub, new Implementation(source, delay),
                    AbstractMeasure::SetHandle()) {}

    /** (Advanced) Restrict the %Delay measure to use only linear interpolation
    to estimate delayed values. By default it uses cubic interpolation whenever
    possible. Cubic interpolation will almost always be better but can be
    unstable in some circumstances. Despite its name this flag also applies
    to extrapolation if we have to do any. This is a topological change. **/
    Delay& setUseLinearInterpolationOnly(bool linearOnly)
    {   updImpl().setUseLinearInterpolationOnly(linearOnly); return *this; }

    /** (Advanced) Allow the %Delay measure to refer to the current
    value when estimating the delayed value. Normally we expect that the current
    value might depend on the delayed value so is not available at the time
    we ask for the delayed value. That means that if the delayed time is
    between the current time and the last saved time (that is, it is a time
    during the current integration step), the measure will have to
    extrapolate from the last-saved values to avoid requiring the current
    value to be available. With this approach the "depends on" time for a %Delay
    measure is just Time stage since it does not depend on any current
    calculations. However, \e extrapolation is much less accurate than
    \e interpolation so if you don't mind the "depends on" stage for a %Delay
    measure being the same stage as for its \a source measure, then you can get
    nicer interpolated values. This is a topological change. **/
    Delay& setCanUseCurrentValue(bool canUseCurrentValue)
    {   updImpl().setCanUseCurrentValue(canUseCurrentValue); return *this; }

    /** Replace the \a source measure. This is a topological change. **/
    Delay& setSourceMeasure(const Measure_<T>& source)
    {   updImpl().setSourceMeasure(source); return *this; }

    /** Change the \a delay time. This is a topological change. **/
    Delay& setDelay(Real delay)
    {   updImpl().setDelay(delay); return *this; }

    /** Return the value of the "use linear interpolation only" flag. **/
    bool getUseLinearInterpolationOnly() const
    {   return getImpl().getUseLinearInterpolationOnly(); }

    /** Return the value of the "can use current value" flag. **/
    bool getCanUseCurrentValue() const
    {   return getImpl().getCanUseCurrentValue(); }

    /** Obtain a reference to the \a source %Measure. **/
    const Measure_<T>& getSourceMeasure() const
    {   return getImpl().getSourceMeasure(); }

    /** Get the amount of time by which this %Measure is delaying its
    \a source %Measure. **/
    Real getDelay() const
    {   return getImpl().getDelay(); }

    /** @cond **/
    SimTK_MEASURE_HANDLE_POSTSCRIPT(Delay, Measure_<T>);
    /** @endcond **/
};

//==============================================================================
//                              SAMPLE AND HOLD
//==============================================================================
/** NOT IMPLEMENTED YET --
This is a Measure operator which, upon occurrence of a designated event,
samples its source Measure and then holds its value in a discrete state
variable until the next occurrence of the event. Any type of data can be
sampled this way.

Information available from this Measure:
 - the held value (Time stage)
 - time of last sample
 - the DiscreteVariableIndex holding the sampled value
 - a reference to the operand Measure

Study initialization is always considered a sampling event.
This measure has no time derivative. **/
template <class T>
class Measure_<T>::SampleAndHold : public Measure_<T> {
public:
    SimTK_MEASURE_HANDLE_PREAMBLE(SampleAndHold, Measure_<T>);

    SampleAndHold(Subsystem& sub, const Measure_<T>& source, EventId e);

    /// Set the held value to a particular value, unrelated to the source.
    /// The time stamp will be taken from the supplied State.
    void setValue(State& s, const T& value) const;

    /// %Force this Measure to sample its input at the current time.
    void sample(State& s) const;

    const Measure_<T>& getSource() const;
    EventId            getEventId() const;

    SampleAndHold& setSource(const Measure_<T>& s);
    SampleAndHold& setEventId(EventId);

    SimTK_MEASURE_HANDLE_POSTSCRIPT(SampleAndHold, Measure_<T>);
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_MEASURE_H_
