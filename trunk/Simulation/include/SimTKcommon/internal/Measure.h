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
 * This file declares the base class Measure for all derived Measure 
 * handle classes, and the handle classes for built-in Measures. Measure handles
 * provide the end user API, while the implementations ("Guts") of Measures
 * derive from the abstract Measure::Guts class defined in MeasureGuts.h.
 * Measure Guts classes provide the Measure developer's API.
 */

#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"

/**
 * Every measure handle class "MH" derived directly or indirectly from the
 * abstract measure handle class "Measure" should include this macro at
 * the beginning of the "public" section of its declaration. It performs
 * the following declarations:
 * <pre>
 *   MH::Guts                   the handle's local implementation class
 *   MH::MH()                   default constructor creates an empty handle
 *   MH::MH(Guts*)              create a handle referencing an existing object
 *   MH::MH(Subsystem&, Guts*)  create a handle referencing an existing but
 *                                unowned object, then installs that in the
 *                                given Subsystem which becomes the owner
 * </pre>
 *
 * MH::Guts must be defined elsewhere as a class derived from Measure::Guts.
 */
#define SimTK_MEASURE_HANDLE_PREAMBLE(MH,PH)    \
    class Guts;                                 \
    explicit MH(Guts* g=0) : PH(g) {}           \
    MH(Subsystem& sub, Guts* g) : PH(sub,g) {}

/**
 * Every measure handle class "MH" derived directly or indirectly from the
 * abstract measure handle class "Measure" should include this macro at
 * the end of the "public" section of its declaration. The macro expects
 * there to be a local class, MH::Guts,
 * already declared. (MH::Guts is the type of MH's implementation object to be
 * derived from Measure::Guts and defined elsewhere.) Then the following
 * type-safe downcast methods will be added to MH's definition:
 * <pre>
 *   MH::getAs(const Measure&)      generic handle to const MH (static)
 *   MH::updAs(Measure&)            generic handle to writable MH (static)
 *   MH::isA(Measure&)              test if generic handle is of type MH
 *   getGuts(const Measure::Guts&)  generic implementation to const MH::Guts
 *   updGuts(Measure::Guts&)        generic implementation to writable MH::Guts
 * </pre>
 * Type checking for the handle class conversions is done only in Debug
 * builds.
 */
#define SimTK_MEASURE_HANDLE_POSTSCRIPT(MH,PH)                  \
    static bool isA(const Measure& m)                           \
    {   return dynamic_cast<const Guts*>(&m.getGuts()) != 0; }  \
    static const MH& getAs(const Measure& m)                    \
    {   assert(isA(m)); return static_cast<const MH&>(m); }     \
    static MH& updAs(Measure& m)                                \
    {   assert(isA(m)); return static_cast<MH&>(m); }           \
    const Guts& getGuts() const                                 \
    {   return dynamic_cast<const Guts&>(Measure::getGuts());}  \
    Guts& updGuts()                                             \
    {   return dynamic_cast<Guts&>(Measure::updGuts());} 

namespace SimTK {

class State;
class Subsystem;
class System;
class EventIndex;

/// Define a unique integral type for safe indexing of Measures. 
SimTK_DEFINE_UNIQUE_INDEX_TYPE(MeasureIndex);

    /////////////
    // MEASURE //
    /////////////

/**
 * This is the base class for all Measure handle classes. This class is not
 * templatized, and represents a Measure generically without knowledge of its
 * value type. This is useful for most of the basic operations that are performed
 * on measures, such as realization, adoption by a Subsystem, and other bookkeeping
 * tasks. For most user purposes, the still-generic derived class Measure_<T> is
 * a more useful handle since its value, of known type T, can be obtained. All
 * the built-in concrete Measure types derive from Measure_<T> rather than Measure.
 * The various concrete Measures, built in or otherwise, may set restrictions on
 * the kinds of types that are allowable as the template argument.
 *
 * Note that handles always consist of exactly one pointer, and handle classes are
 * always concrete (meaning they have no virtual methods).
 */
class SimTK_SimTKCOMMON_EXPORT Measure {
public:
    class Guts; // local; name is Measure::Guts

    explicit Measure(Guts* g=0);
    Measure(Subsystem&, Guts* g);
    Measure(const Measure&);
    Measure& operator=(const Measure&);
    ~Measure();

    // At what Stage can we expect the value of this Measure to be available?
    Stage getValueDependence(const State&) const;

    // There can be multiple handles on the same Measure_.
    bool isSameMeasure(const Measure& other) const
    {   return guts && guts==other.guts;}

    bool isInSubsystem() const;
    const Subsystem& getSubsystem() const;
    MeasureIndex     getSubsystemMeasureIndex() const;

    // Internal use only

    // dynamic_cast the returned reference to a reference to your concrete Guts
    // class.
    const Guts& getGuts() const {assert(guts); return *guts;}
    Guts&       updGuts()       {assert(guts); return *guts;}

    bool hasGuts() const {return guts!=0;}

public:
    // These are built-in Measures with local class names. 

    // Templatized measures may have restrictions on the allowable template
    // type and may be specialized for particular types.
    template <class T> class Constant_;     // T is any assignable type
    template <class T> class SampleAndHold_;//    "

    // These accept any type that supports operator<, elementwise for 
    // vectors and matrices.
    template <class T> class Minimum_;
    template <class T> class Maximum_;

    // These accept floating point numerical template arguments only.
    template <class T> class Integrate_;
    template <class T> class Sinusoid_;



    // Default template arguments for the built-in Measures.
    typedef Constant_<Real> Constant;
    typedef Integrate_<Real> Integrate;
    typedef Sinusoid_<Real> Sinusoid;


private:
    // This is the only data member in this class. Also, any class derived from
    // Measure must have *NO* data members at all (data goes in the Guts class).
    Guts* guts;

friend class Guts;
};


    /////////////////
    // MEASURE_<T> //
    /////////////////

/**
 * This is the base handle class for all Measures whose value type is known.
 */
template <class T>
class Measure_ : public Measure {
public:
    SimTK_MEASURE_HANDLE_PREAMBLE(Measure_, Measure);

    const T& getValue(const State& s) const {return getGuts().getValue(s);}

    SimTK_MEASURE_HANDLE_POSTSCRIPT(Measure_, Measure);
};


    //////////////////
    // CONSTANT_<T> //
    //////////////////

/**
 * This creates a Measure whose value is a Topology-stage constant of any
 * type T. Changing the constant's value invalidates the containing Subsystem's
 * Topology.
 */
template <class T>
class Measure::Constant_ : public Measure_<T> {
public:
    SimTK_MEASURE_HANDLE_PREAMBLE(Constant_, Measure_<T>);

    Constant_(Subsystem& sub, const T& value)
    :   Measure_<T>(sub, new Guts(value)) {}

    Constant_(System& sys, const T& value)
    {   new (this) Constant_(sys.updDefaultSubsystem(), value); }

    /// Note that this does not require a State since it is a Topology-stage change.
    Constant_& setValue(const T& value) {updGuts().setValue(value); return *this;}

    SimTK_MEASURE_HANDLE_POSTSCRIPT(Constant_, Measure_<T>);
};


    //////////////////
    // SINUSOID_<T> //
    //////////////////

/**
 * This measure produces a sinusoidal function of time:
 *  m(t) = a*sin(w*t+p)
 * where a=amplitude in arbitrary units, w=frequency in rad/s, p=phase in rad.
 */
template <class T>
class Measure::Sinusoid_ : public Measure_<T> {
public:
    SimTK_MEASURE_HANDLE_PREAMBLE(Sinusoid_, Measure_<T>);

    Sinusoid_(Subsystem& sub,
              const T& amplitude, 
              const T& frequency,
              const T& phase=T(0))
    :   Measure_<T>(sub, new Guts(amplitude,frequency,phase)) {}

    Sinusoid_(System& sys,
              const T& amplitude, 
              const T& frequency,
              const T& phase=T(0))
    {   new (this) Sinusoid(sys.updDefaultSubsystem(), amplitude, frequency, phase); }

    SimTK_MEASURE_HANDLE_POSTSCRIPT(Sinusoid_, Measure_<T>);
};


    ///////////////////
    // INTEGRATE_<T> //
    ///////////////////

template <class T>
class Measure::Integrate_ : public Measure_<T> {
public:
    SimTK_MEASURE_HANDLE_PREAMBLE(Integrate_, Measure_<T>);

    Integrate_(Subsystem& sub, const Measure_<T>& deriv, const Measure_<T>& ic)
    :   Measure_<T>(sub, new Guts(deriv,ic)) {}
    Integrate_(System& sys, const Measure_<T>& deriv, const Measure_<T>& ic)
    {   new (this) Integrate_(sys.updDefaultSubsystem(), deriv, ic);}

    void setValue(State& s, const T& value) const {return getGuts().setValue(s, value);}
    const Measure_<T>& getDerivativeMeasure() const {return getGuts().getDerivativeMeasure();}
    const Measure_<T>& getInitialConditionMeasure() const {return getGuts().getInitialConditionMeasure();}

    Integrate_& setDerivativeMeasure(const Measure_<T>& d)
    {   updGuts().setDerivativeMeasure(d); return *this; }
    Integrate_& setInitialConditionMeasure(const Measure_<T>& ic)
    {   updGuts().setInitialConditionMeasure(ic); return *this; }

    SimTK_MEASURE_HANDLE_POSTSCRIPT(Integrate_, Measure_<T>);
};


    /////////////////
    // MINIMUM_<T> //
    /////////////////

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
class Measure::Minimum_ : public Measure_<T> {
public:
    SimTK_MEASURE_HANDLE_PREAMBLE(Minimum_, Measure_<T>);

    Minimum_(Subsystem& sub, const Measure_<T>& source, const Measure_<T>& sourceDot)
    :   Measure_<T>(sub, new Guts(source, sourceDot)) {}
    Minimum_(System& sys, const Measure_<T>& source, const Measure_<T>& sourceDot)
    {   new (this) Minimum_(sys.updDefaultSubsystem(), source, sourceDot);}

    void setValue(State& s, const T& value) const {return getGuts().setValue(s, value);}
    const Measure_<T>& getSourceMeasure() const {return getGuts().getSourceMeasure();}
    const Measure_<T>& getSourceDerivativeMeasure() const {return getGuts().getSourceDerivativeMeasure();}

    Minimum_& setSourceMeasure(const Measure_<T>& s)
    {   updGuts().setSourceMeasure(s); return *this; }
    Minimum_& setSourceDerivativeMeasure(const Measure_<T>& sdot)
    {   updGuts().setSourceDerivativeMeasure(sdot); return *this; }

    SimTK_MEASURE_HANDLE_POSTSCRIPT(Minimum_, Measure_<T>);
};

/**
 * This is a Measure operator which, upon occurrence of a designated event, samples
 * its source Measure and then holds its value in a discrete state variable
 * until the next occurrence of the event. Any type of data can be sampled this way.
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
class Measure::SampleAndHold_ : public Measure_<T> {
public:
    SimTK_MEASURE_HANDLE_PREAMBLE(SampleAndHold_, Measure_<T>);

    SampleAndHold_(Subsystem& sub, const Measure_<T>& source, EventIndex e);
    SampleAndHold_(System&    sys, const Measure_<T>& source, EventIndex e);

    /// Set the held value to a particular value, unrelated to the source.
    /// The time stamp will be taken from the supplied State.
    void setValue(State& s, const T& value) const;

    /// Force this Measure to sample its input at the current time.
    void sample(State& s) const;

    const Measure_<T>& getSource() const;
    EventIndex         getEvent() const;

    SampleAndHold_& setSource(const Measure_<T>& s);
    SampleAndHold_& setEvent(EventIndex);

    SimTK_MEASURE_HANDLE_POSTSCRIPT(SampleAndHold_, Measure_<T>);
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_MEASURE_H_
