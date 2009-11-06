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
 * Portions copyright (c) 2008-9 Stanford University and the Authors.         *
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
#define SimTK_MEASURE_HANDLE_PREAMBLE(MH,PH) \
    class Implementation;                                           \
    explicit MH() : PH(new Implementation()) {}                     \
    explicit MH(Implementation* imp) : PH(imp) {}                   \
    MH(Subsystem& sub, Implementation* imp, const SetHandle& sh)    \
    :   PH(sub,imp,sh) {}

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
class EventIndex;

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
	// This is used to make sure the automatically-generated handle
	// constructor doesn't conflict with an explicitly-defined one.
	class SetHandle {};

public:
    class Implementation; // local; name is AbstractMeasure::Implementation

    explicit AbstractMeasure(Implementation* g=0);
    AbstractMeasure(Subsystem&, Implementation* g, const SetHandle&);
    AbstractMeasure(const AbstractMeasure&);
    AbstractMeasure& operator=(const AbstractMeasure&);
    ~AbstractMeasure();

    /// Every Measure can produce a value, and some can provide one or
    /// more total derivatives with respect to time of that value. This
    /// method reports how many are available: 1 -> first derivative
    /// d/dt is available, 2 -> first and second derivative d^2/dt^2
    /// are available, etc. We interpret the "0th derivative" to be the
    /// Measure's value.
    /// @return The maximum available derivative order.
    inline int getNumTimeDerivatives() const;

    /// At what Stage can we expect the value of this AbstractMeasure or
    /// one of its time derivatives to be available? Users of Measures will 
    /// typically impose restrictions on the levels they will accept.
    /// @param[in] state 
    ///     The State to be consulted to determine the maximum Stage on
    ///     which the requested value may depend, in case that Stage is
    ///     State-dependendent (e.g., on a Model-stage variable).
    /// @param[in] derivOrder
    ///     Which derivative level is to be checked: 0 -> the value,
    ///     1 -> the 1st time derivative, etc. Must not be higher than the
    ///     value returned by getNumTimeDerivatives().
    /// @return The Stage after which this value is available.                  
    inline Stage 
    getDependsOnStage(const State& state, int derivOrder=0) const;


    /// There can be multiple handles on the same Measure.
    bool isSameMeasure(const AbstractMeasure& other) const
    {   return impl && impl==other.impl;}

    /// Test whether this Measure object has been adopted by a Subsystem.
    inline bool isInSubsystem() const;
    /// Return a reference to the Subsystem that owns this Measure. Will
    /// throw an exception if the Measure is not currently owned by any
    /// Subsystem.
    inline const Subsystem& getSubsystem()  const;
    /// Return the MeasureIndex by which this Measure is know to the Subsystem 
    /// that owns it. Will throw an exception if the Measure is not currently 
    /// owned by any Subsystem.
    inline MeasureIndex getSubsystemMeasureIndex() const;

    // Internal use only

    // dynamic_cast the returned reference to a reference to your concrete
    // Implementation class.
    const Implementation& getImpl() const {assert(impl); return *impl;}
    Implementation&       updImpl()       {assert(impl); return *impl;}
    bool                  hasImpl() const {return impl!=0;}


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
    SimTK_MEASURE_HANDLE_PREAMBLE(Measure_, AbstractMeasure);

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
    class Variable;     // T is any assignable type
    class SampleAndHold;//    "

    // This requires any numerical type.
    class Plus;

    // These accept any type that supports operator<, elementwise for 
    // vectors and matrices.
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
    :   Measure_<T>(sub, new Implementation(value), SetHandle()) {}

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

    explicit Zero(Subsystem& sub)
    :   Measure_<T>::Constant(sub, new Implementation(), SetHandle()) {}

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

    explicit One(Subsystem& sub)
    :   Measure_<T>::Constant(sub, new Implementation()) {}

    SimTK_MEASURE_HANDLE_POSTSCRIPT(One, Measure_<T>::Constant);
private:
    One& setValue(const T&)     
    {   return *reinterpret_cast<Zero*>(0);} // suppressed!
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
                    SetHandle()) {}

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
    :   Measure_<T>(sub, new Implementation(amplitude,frequency,phase), SetHandle()) {}

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
    :   Measure_<T>(sub,
                    new Implementation(left, right), SetHandle())
    {   SimTK_ERRCHK_ALWAYS
           (   this->getSubsystem().isSameSubsystem(left.getSubsystem())
            && this->getSubsystem().isSameSubsystem(right.getSubsystem()),
            "Measure_<T>::Plus::ctor()",
            "Arguments must be in the same Subsystem as this Measure.");
    }


    SimTK_MEASURE_HANDLE_POSTSCRIPT(Plus, Measure_<T>);
};

    ///////////////
    // INTEGRATE //
    ///////////////

template <class T>
class Measure_<T>::Integrate : public Measure_<T> {
public:
    SimTK_MEASURE_HANDLE_PREAMBLE(Integrate, Measure_<T>);

    Integrate(Subsystem& sub, const Measure_<T>& deriv, const Measure_<T>& ic)
    :   Measure_<T>(sub, new Implementation(deriv,ic), SetHandle()) {}

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


    /////////////
    // MINIMUM //
    /////////////

/**
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
    :   Measure_<T>(sub, new Implementation(source, sourceDot), SetHandle()) {}

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

/**
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

    SampleAndHold(Subsystem& sub, const Measure_<T>& source, EventIndex e);

    /// Set the held value to a particular value, unrelated to the source.
    /// The time stamp will be taken from the supplied State.
    void setValue(State& s, const T& value) const;

    /// Force this Measure to sample its input at the current time.
    void sample(State& s) const;

    const Measure_<T>& getSource() const;
    EventIndex         getEvent() const;

    SampleAndHold& setSource(const Measure_<T>& s);
    SampleAndHold& setEvent(EventIndex);

    SimTK_MEASURE_HANDLE_POSTSCRIPT(SampleAndHold, Measure_<T>);
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_MEASURE_H_
