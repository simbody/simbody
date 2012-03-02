#ifndef SimTK_SimTKCOMMON_MEASURE_H_
#define SimTK_SimTKCOMMON_MEASURE_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-10 Stanford University and the Authors.        *
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
 * Every measure handle class "MH" derived directly or indirectly from the
 * AbstractMeasure handle class should include this macro at
 * the beginning of the "public" section of its declaration. It performs
 * the following declarations:
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
#define SimTK_MEASURE_HANDLE_PREAMBLE_BASE(MH,PH) \
    class Implementation;                                           \
    explicit MH(Implementation* imp) : PH(imp) {}                   \
    MH(Subsystem& sub, Implementation* imp,                         \
       const AbstractMeasure::SetHandle& sh)                        \
    :   PH(sub,imp,sh) {}                                           \
    MH& operator=(const MH& src) {PH::operator=(src); return *this;}\
    MH& shallowAssign(const MH& src) {PH::shallowAssign(src); return *this;}\
    MH& deepAssign(const MH& src) {PH::deepAssign(src); return *this;}

// The default constructor for concrete classes should instantiate
// a default-constructed Implementation object if no Implementation
// is provided.
#define SimTK_MEASURE_HANDLE_PREAMBLE(MH,PH)    \
    SimTK_MEASURE_HANDLE_PREAMBLE_BASE(MH,PH)   \
    MH() : PH(new Implementation()) {}          \
    explicit MH(Subsystem& sub)                 \
    : PH(sub,new Implementation(), typename PH::SetHandle()) {}

// The default constructor for a still-abstract derived class can't
// instantiate an Implementation.
#define SimTK_MEASURE_HANDLE_PREAMBLE_ABSTRACT(MH,PH)   \
    SimTK_MEASURE_HANDLE_PREAMBLE_BASE(MH,PH)           \
    MH() : PH() {}                                      \
    explicit MH(Subsystem& sub) : PH(sub) {}    

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
 *   getImpl(const AbstractMeasure::Implementation&)  
 *      generic implementation to const MH::Implementation
 *   updImpl(AbstractMeasure::Implementation&)
 *      generic implementation to writable MH::Implementation
 * </pre>
 * Type checking for the handle class conversions is done only in Debug
 * builds.
 */
#define SimTK_MEASURE_HANDLE_POSTSCRIPT(MH,PH) \
    static bool isA(const AbstractMeasure& m)                               \
    {   return dynamic_cast<const Implementation*>(&m.getImpl()) != 0; }    \
    static const MH& getAs(const AbstractMeasure& m)                        \
    {   assert(isA(m)); return static_cast<const MH&>(m); }                 \
    static MH& updAs(AbstractMeasure& m)                                    \
    {   assert(isA(m)); return static_cast<MH&>(m); }                       \
    const Implementation& getImpl() const                                   \
    {   return dynamic_cast<const Implementation&>                          \
                    (AbstractMeasure::getImpl());}                          \
    Implementation& updImpl()                                               \
    {   return dynamic_cast<Implementation&>(AbstractMeasure::updImpl());} 

namespace SimTK {

class State;
class Subsystem;
class System;
class EventId;

/// Define a unique integral type for safe indexing of Measures. 
SimTK_DEFINE_UNIQUE_INDEX_TYPE(MeasureIndex);

    //////////////////////
    // ABSTRACT MEASURE //
    //////////////////////

/**
 * This is the base class for all Measure handle classes. This class is not
 * templatized, and represents a Measure generically without knowledge of its
 * value type. This is useful for most of the basic operations that are 
 * performed on measures, such as realization, adoption by a Subsystem, and 
 * other bookkeeping tasks. For most user purposes, the still-generic derived 
 * class Measure_<T> is a more useful handle since its value, of known type T, 
 * can be obtained. All the built-in concrete Measure types derive from 
 * Measure_<T> rather than AbstractMeasure. The various concrete Measures, 
 * built in or otherwise, may set restrictions on the kinds of types that 
 * are allowable as the template argument.
 *
 * Note that handles always consist of exactly one pointer, and handle classes
 * are always concrete (meaning they have no virtual methods).
 */
class SimTK_SimTKCOMMON_EXPORT AbstractMeasure {
protected:
	// An object of this type is used as a dummy argument to make sure the 
    // automatically-generated handle constructor's signature doesn't conflict 
    // with an explicitly-defined one.
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


    /////////////////
    // MEASURE_<T> //
    /////////////////

/**
 * This is the base handle class for all Measures whose value type is known.
 */
template <class T>
class Measure_ : public AbstractMeasure {
public:
    /// This class is still abstract so we don't want it to allocate an
    /// Implementation object in its default constructor.
    SimTK_MEASURE_HANDLE_PREAMBLE_ABSTRACT(Measure_, AbstractMeasure);

    /// Retrieve the Value of this Measure or one of its time
    /// derivatives, assuming the State has been realized to at least the
    /// required stage for the desired value, as reported by
    /// getValueDependence(). If the stage is sufficient but the corresponding 
    /// value has not yet been computed, it will be computed first with its 
    /// value going into this State's cache so that subsequent calls do not 
    /// require further computation.
    inline const T& getValue(const State& s, int derivOrder=0) const 
    {   return getImpl().getValue(s,derivOrder); }

    // These are built-in Measures with local class names. 

    // Templatized measures may have restrictions on the allowable template
    // type and may be specialized for particular types.
    class Zero;         // T is any numerical type
    class One;          // T is any numerical type
    class Constant;     // T is any assignable type
    class Time;         // T is any type for which T(t) makes sense.
    class Variable;     // T is any assignable type (state)
    class Result;       // T is any assignable type (cache)
    class SampleAndHold;//    "

    // This requires any numerical type.
    class Plus;
    class Minus;
    class Scale;
    class Differentiate;

    // These accept any type that supports operator<, elementwise for 
    // vectors and matrices. TODO
    class Minimum;
    class Maximum;

    // These accept floating point numerical template arguments only.
    class Integrate;
    class Sinusoid;

    SimTK_MEASURE_HANDLE_POSTSCRIPT(Measure_, AbstractMeasure);
};

/// A convenient abbreviation for the most common kind of Measure --
/// one which returns a single Real result.
typedef Measure_<Real> Measure;


    //////////////
    // CONSTANT //
    //////////////

/**
 * This creates a Measure whose value is a Topology-stage constant of any
 * type T. This does not have to be part of a Subsystem, but if it is then
 * changing the constant's value invalidates the containing Subsystem's 
 * Topology.
 */
template <class T>
class Measure_<T>::Constant : public Measure_<T> {
public:
    SimTK_MEASURE_HANDLE_PREAMBLE(Constant, Measure_<T>);

    explicit Constant(const T& value)
    :   Measure_<T>(new Implementation(value)) {}

    Constant(Subsystem& sub, const T& value)
    :   Measure_<T>(sub, new Implementation(value), 
                    AbstractMeasure::SetHandle()) {}

    /// Note that this does not require a State since it is a Topology-stage 
    /// change.
    Constant& setValue(const T& value) 
    {   updImpl().setValue(value); return *this; }

    SimTK_MEASURE_HANDLE_POSTSCRIPT(Constant, Measure_<T>);
};

    //////////
    // ZERO //
    //////////

/**
 * This creates a Measure::Constant whose value is always T(0) and can't
 * be changed.
 */
template <class T>
class Measure_<T>::Zero : public Measure_<T>::Constant {
public:
    SimTK_MEASURE_HANDLE_PREAMBLE(Zero, Measure_<T>::Constant);

    SimTK_MEASURE_HANDLE_POSTSCRIPT(Zero, Measure_<T>::Constant);
private:
    Zero& setValue(const T&) 
    {   return *reinterpret_cast<Zero*>(0);} // suppressed!
};

    /////////
    // ONE //
    /////////

/**
 * This creates a Measure::Constant whose value is always T(1) and can't
 * be changed.
 */
template <class T>
class Measure_<T>::One : public Measure_<T>::Constant {
public:
    SimTK_MEASURE_HANDLE_PREAMBLE(One, Measure_<T>::Constant);

    SimTK_MEASURE_HANDLE_POSTSCRIPT(One, Measure_<T>::Constant);
private:
    One& setValue(const T&)     
    {   return *reinterpret_cast<One*>(0);} // suppressed!
};

    //////////
    // TIME //
    //////////

/**
 * This creates a Measure::Time whose value is always T(time).
 */
template <class T>
class Measure_<T>::Time : public Measure_<T> {
public:
    SimTK_MEASURE_HANDLE_PREAMBLE(Time, Measure_<T>);

    SimTK_MEASURE_HANDLE_POSTSCRIPT(Time, Measure_<T>);
};

    //////////////
    // VARIABLE //
    //////////////

/**
 * This creates a Measure whose value is a discrete State variable
 * of any type T.
 */
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

    /// Note that this does not require a State since it is a Topology-stage 
    /// change.
    Variable& setDefaultValue(const T& defaultValue) 
    {   updImpl().setDefaultValue(defaultValue); return *this; }

    const T& getDefaultValue() const
    {   return getImpl().getDefaultValue(); }

    void setValue(State& state, const T& value) const
    {   getImpl().setValue(state, value); }

    SimTK_MEASURE_HANDLE_POSTSCRIPT(Variable, Measure_<T>);
};

    ////////////
    // RESULT //
    ////////////

/**
 * This Measure holds the result of some externally-determined computation,
 * and helps to coordinate the validity of that computation with respect
 * to the state variables. The value must be set manually and explicitly
 * marked valid when it is complete. The value will be automatically 
 * invalidated after a state change at or below a specified Stage, and changing
 * the value here will automatically invalidate a specified Stage and above.
 *
 * In constrast to Measure::Variable, Measure::Result is not a state variable;
 * it is a cache variable meaning that it works with a const State. It is
 * expected that the result can be recalculated from the state variables when
 * needed, or contains ephemeral information that can be discarded.
 *
 * No provision for derivatives is made; there is only the one result.
 */
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
    bool isValid(const State& state) const {return getImpl().isValid(state);}

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

    //////////////
    // SINUSOID //
    //////////////

/**
 * This measure produces a sinusoidal function of time:
 * <pre>
 *      m(t) = a*sin(w*t+p)
 * </pre>
 * where a=amplitude in arbitrary units, w=frequency in rad/unit time, 
 * p=phase in rad.
 */
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

    //////////
    // PLUS //
    //////////

/**
 * This Measure is the sum of two Measures of the same type T.
 */
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

    ///////////
    // MINUS //
    ///////////

/**
 * This Measure is the difference of two Measures of the same type T.
 */
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

    ///////////
    // SCALE //
    ///////////

/**
 * This Measure multiplies some other Measure by a Real scale factor.
 */
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

    SimTK_MEASURE_HANDLE_POSTSCRIPT(Scale, Measure_<T>);
};

    ///////////////
    // INTEGRATE //
    ///////////////

template <class T>
class Measure_<T>::Integrate : public Measure_<T> {
public:
    SimTK_MEASURE_HANDLE_PREAMBLE(Integrate, Measure_<T>);

    Integrate(Subsystem& sub, const Measure_<T>& deriv, const Measure_<T>& ic)
    :   Measure_<T>(sub, new Implementation(deriv,ic), 
                    AbstractMeasure::SetHandle()) {}

    void setValue(State& s, const T& value) const 
    {   return getImpl().setValue(s, value); }
    const Measure_<T>& getDerivativeMeasure() const 
    {   return getImpl().getDerivativeMeasure(); }
    const Measure_<T>& getInitialConditionMeasure() const 
    {   return getImpl().getInitialConditionMeasure(); }

    Integrate& setDerivativeMeasure(const Measure_<T>& d)
    {   updImpl().setDerivativeMeasure(d); return *this; }
    Integrate& setInitialConditionMeasure(const Measure_<T>& ic)
    {   updImpl().setInitialConditionMeasure(ic); return *this; }

    SimTK_MEASURE_HANDLE_POSTSCRIPT(Integrate, Measure_<T>);
};

    ///////////////////
    // DIFFERENTIATE //
    ///////////////////

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

    /////////////
    // MINIMUM //
    /////////////

/** NOT IMPLEMENTED YET --
 * This Measure tracks the minimum value attained by its source operand
 * since the last initialize() call. If the time derivative of the 
 * source operand is available, the measure will arrange to ensure
 * precise isolation of the minimum values by defining a triggered 
 * event that watches for negative-to-positive sign changes of the derivative.
 * Then if you output reporting data upon the occurrence of triggered
 * events (as well as your regularly scheduled output) your data will
 * include the precise minimum (to within a specifiable isolation time window).
 *
 * Information available from this Measure:
 *  - the current value (at the source's stage)
 *  - the value at the start of the current step (Time stage)
 *  - time of last change
 *  - the DiscreteVariableIndex holding the previous value
 *  - the CacheEntryIndex designated for the current/next value
 *  - a reference to the source Measure
 *
 * And if a source derivative is available
 *  - the EventId of the minimum-trapping event if there is one
 *  - the time derivative of the Minimum measure
 *  - a reference to the source derivative Measure
 *
 * The time derivative fdot of f(t)=min_t0_t(s(t')) where s is the
 * source measure and t0 <= t' <= t is
 * <pre>
 *      fdot(t) = s(t) < f(ti) && sdot(t) < 0 ? sdot(t) : 0
 * </pre>
 * where ti is the time at the start of the current step.
 *
 * At the start of a continuous interval, the updated value (if any)
 * replaces the stored value.
 */

template <class T>
class Measure_<T>::Minimum : public Measure_<T> {
public:
    SimTK_MEASURE_HANDLE_PREAMBLE(Minimum, Measure_<T>);

    Minimum(Subsystem& sub, const Measure_<T>& source, const Measure_<T>& sourceDot)
    :   Measure_<T>(sub, new Implementation(source, sourceDot), 
                    AbstractMeasure::SetHandle()) {}

    void setValue(State& s, const T& value) const 
    {   return getImpl().setValue(s, value); }
    const Measure_<T>& getSourceMeasure() const 
    {   return getImpl().getSourceMeasure(); }
    const Measure_<T>& getSourceDerivativeMeasure() const 
    {   return getImpl().getSourceDerivativeMeasure(); }

    Minimum& setSourceMeasure(const Measure_<T>& s)
    {   updImpl().setSourceMeasure(s); return *this; }
    Minimum& setSourceDerivativeMeasure(const Measure_<T>& sdot)
    {   updImpl().setSourceDerivativeMeasure(sdot); return *this; }

    SimTK_MEASURE_HANDLE_POSTSCRIPT(Minimum, Measure_<T>);
};

/** NOT IMPLEMENTED YET --
 * This is a Measure operator which, upon occurrence of a designated event, 
 * samples its source Measure and then holds its value in a discrete state 
 * variable until the next occurrence of the event. Any type of data can be 
 * sampled this way.
 *
 * Information available from this Measure:
 *  - the held value (Time stage)
 *  - time of last sample
 *  - the DiscreteVariableIndex holding the sampled value
 *  - a reference to the operand Measure
 *
 * Study initialization is always considered a sampling event.
 * This measure has no time derivative.
 */
template <class T>
class Measure_<T>::SampleAndHold : public Measure_<T> {
public:
    SimTK_MEASURE_HANDLE_PREAMBLE(SampleAndHold, Measure_<T>);

    SampleAndHold(Subsystem& sub, const Measure_<T>& source, EventId e);

    /// Set the held value to a particular value, unrelated to the source.
    /// The time stamp will be taken from the supplied State.
    void setValue(State& s, const T& value) const;

    /// Force this Measure to sample its input at the current time.
    void sample(State& s) const;

    const Measure_<T>& getSource() const;
    EventId            getEventId() const;

    SampleAndHold& setSource(const Measure_<T>& s);
    SampleAndHold& setEventId(EventId);

    SimTK_MEASURE_HANDLE_POSTSCRIPT(SampleAndHold, Measure_<T>);
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_MEASURE_H_
