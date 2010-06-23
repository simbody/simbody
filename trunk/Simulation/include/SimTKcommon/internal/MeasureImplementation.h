#ifndef SimTK_SimTKCOMMON_MEASURE_IMPLEMENTATION_H_
#define SimTK_SimTKCOMMON_MEASURE_IMPLEMENTATION_H_

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

#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"
#include "SimTKcommon/internal/State.h"
#include "SimTKcommon/internal/Measure.h"
#include "SimTKcommon/internal/Subsystem.h"
#include "SimTKcommon/internal/System.h"
#include "SimTKcommon/internal/SubsystemGuts.h"

#include <cmath>


namespace SimTK {


    /////////////////////////////
    // MEASURE::IMPLEMENTATION //
    /////////////////////////////

/**
 * The abstract parent of all Measure Implementation classes.
 */
class SimTK_SimTKCOMMON_EXPORT AbstractMeasure::Implementation {
protected:
    // This constructor is for use by concrete Measure::Implementations. Note 
    // that this serves as a default constructor since the argument has a 
    // default.
    explicit Implementation(const std::string& name="<NONAME>")
    :   measureName(name), copyNumber(0), mySubsystem(0), refCount(0) {}

    // Base class copy constructor copies the name, removes the Subsystem
    // and sets the reference count to zero. This gets used by the clone()
    // methods in the concrete classes.
    Implementation(const Implementation& src)
    :   measureName(src.measureName), copyNumber(src.copyNumber+1),
        mySubsystem(0), refCount(0) {}
    
    // Base class copy assignment operator copies the name, removes the
    // Subsystem, and sets the reference count to zero. This is probably
    // not used.
    Implementation& operator=(const Implementation& src) {
        if (&src != this)
        {   measureName=src.measureName; copyNumber=src.copyNumber+1;
            refCount=0; mySubsystem=0; }
        return *this; 
    }

    // destructor is virtual

    // Increment the reference count and return its new value.
    int incrRefCount() const {return ++refCount;}

    // Decrement the reference count and return its new value.
    int decrRefCount() const {return --refCount;}

    // Get the current value of the reference counter.
    int getRefCount() const {return refCount;}

    const std::string& getName()        const {return measureName;}
    int                getCopyNumber()  const {return copyNumber;}

    // This is a deep copy of the concrete Implementation object, except the
    // Subsystem will have been removed. The reference count on the new object
    // will be zero; be sure to increment it if you put it in a handle.
    Implementation* clone() const {return cloneVirtual();}

    // realizeTopology() is pure virtual below for Measure_<T> to supply.
    void realizeModel       (State& s)       const {realizeMeasureModelVirtual(s);}
    void realizeInstance    (const State& s) const {realizeMeasureInstanceVirtual(s);}
    void realizeTime        (const State& s) const {realizeMeasureTimeVirtual(s);}
    void realizePosition    (const State& s) const {realizeMeasurePositionVirtual(s);}
    void realizeVelocity    (const State& s) const {realizeMeasureVelocityVirtual(s);}
    void realizeDynamics    (const State& s) const {realizeMeasureDynamicsVirtual(s);}
    void realizeAcceleration(const State& s) const {realizeMeasureAccelerationVirtual(s);}
    void realizeReport      (const State& s) const {realizeMeasureReportVirtual(s);}

    // This should be called at the start of a numerical integration study to
    // cause initial conditions to get set for any Measures that have integrated
    // state variables.
    void initialize(State& s) const {initializeVirtual(s);}

    int getNumTimeDerivatives() const {return getNumTimeDerivativesVirtual();}

    Stage getDependsOnStage(int derivOrder) const {
        SimTK_ERRCHK2(0 <= derivOrder && derivOrder <= getNumTimeDerivatives(),
            "Measure::getDependsOnStage()",
            "derivOrder %d was out of range; this Measure allows 0-%d.",
            derivOrder, getNumTimeDerivatives()); 
        return getDependsOnStageVirtual(derivOrder); 
    }


    void setSubsystem(Subsystem& sub, MeasureIndex mx) 
    {   assert(!mySubsystem && mx.isValid()); 
        mySubsystem = &sub; myIndex = mx; }

    bool             isInSubsystem() const {return mySubsystem != 0;}
    const Subsystem& getSubsystem() const {assert(mySubsystem); return *mySubsystem;}
    Subsystem&       updSubsystem() {assert(mySubsystem); return *mySubsystem;}
    MeasureIndex     getSubsystemMeasureIndex() const {assert(mySubsystem); return myIndex;}
    SubsystemIndex   getSubsystemIndex() const
    {   return getSubsystem().getMySubsystemIndex(); }

    void invalidateTopologyCache() const
    {   if (isInSubsystem()) getSubsystem().invalidateSubsystemTopologyCache(); }

    Stage getStage(const State& s) const {return getSubsystem().getStage(s);}

    // VIRTUALS //
    // Ordinals must retain the same meaning from release to release
    // to preserve binary compatibility.

    /* 0*/virtual ~Implementation() {}
    /* 1*/virtual Implementation* cloneVirtual() const = 0;

    /* 2*/virtual void realizeTopology(State&)const = 0;

    /* 3*/virtual void realizeMeasureModelVirtual(State&) const {}
    /* 4*/virtual void realizeMeasureInstanceVirtual(const State&) const {}
    /* 5*/virtual void realizeMeasureTimeVirtual(const State&) const {}
    /* 6*/virtual void realizeMeasurePositionVirtual(const State&) const {}
    /* 7*/virtual void realizeMeasureVelocityVirtual(const State&) const {}
    /* 8*/virtual void realizeMeasureDynamicsVirtual(const State&) const {}
    /* 9*/virtual void realizeMeasureAccelerationVirtual(const State&) const {}
    /*10*/virtual void realizeMeasureReportVirtual(const State&) const {}

    /*11*/virtual void initializeVirtual(State&) const {}
    /*12*/virtual int  
          getNumTimeDerivativesVirtual() const {return 0;}
    /*13*/virtual Stage 
          getDependsOnStageVirtual(int order) const = 0;

private:
    std::string     measureName;
    int             copyNumber; // bumped each time we do a deep copy

    // These are set when this Measure is adopted by a Subsystem.
    Subsystem*      mySubsystem;
    MeasureIndex    myIndex;

    // Measures have shallow copy semantics so they share the Implementation 
    // objects, which are only deleted when the refCount goes to zero.
    mutable int     refCount;

friend class AbstractMeasure;
friend class Subsystem::Guts;
friend class Subsystem::Guts::GutsRep;
};

    //////////////////////////////////
    // ABSTRACT MEASURE DEFINITIONS //
    //////////////////////////////////

// These had to wait for Implementation to be defined.

inline AbstractMeasure::
AbstractMeasure(Implementation* g) 
:   impl(g)
{   if (impl) impl->incrRefCount(); }

inline AbstractMeasure::
AbstractMeasure(Subsystem& sub, Implementation* g, const SetHandle&) 
:   impl(g) {
    SimTK_ERRCHK(hasImpl(), "AbstractMeasure::ctor()",
        "An empty Measure handle can't be put in a Subsystem.");
    impl->incrRefCount();
    sub.adoptMeasure(*this);
}

// Shallow copy constructor.
inline AbstractMeasure::AbstractMeasure(const AbstractMeasure& src) 
:   impl(0) {
    if (src.impl) {
        impl = src.impl;
        impl->incrRefCount();
    }
}

// Shallow assignment.
inline AbstractMeasure& AbstractMeasure::
shallowAssign(const AbstractMeasure& src) {
    if (impl != src.impl) {
        if (impl && impl->decrRefCount()==0) delete impl;
        impl = src.impl;
        impl->incrRefCount();
    }
    return *this;
}

// Note that even if the source and destination are currently pointing
// to the same Implementation, we still have to make a new copy so that
// afterwards the destination has its own, refcount==1 copy.
inline AbstractMeasure& AbstractMeasure::
deepAssign(const AbstractMeasure& src) {
    if (&src != this) {
        if (impl && impl->decrRefCount()==0) delete impl;
        if (src.impl) {
            impl = src.impl->clone();
            impl->incrRefCount();
        } else
            impl = 0;
    }
    return *this;
}

inline AbstractMeasure::
~AbstractMeasure()
{   if (impl && impl->decrRefCount()==0) delete impl;}

inline bool AbstractMeasure::
isInSubsystem() const
{   return hasImpl() && getImpl().isInSubsystem(); }

inline const Subsystem& AbstractMeasure::
getSubsystem() const
{   return getImpl().getSubsystem(); }

inline MeasureIndex AbstractMeasure::
getSubsystemMeasureIndex() const
{   return getImpl().getSubsystemMeasureIndex();}

inline int AbstractMeasure::
getNumTimeDerivatives() const
{   return getImpl().getNumTimeDerivatives(); }

inline Stage AbstractMeasure::
getDependsOnStage(int derivOrder) const
{   return getImpl().getDependsOnStage(derivOrder); }

inline int AbstractMeasure::
getRefCount() const
{   return getImpl().getRefCount(); }


    /////////////////////////////////
    // MEASURE_<T>::IMPLEMENTATION //
    /////////////////////////////////

/**
 * This is the base Implementation class for all Measures whose value type 
 * is known. This class is still abstract but provides many services 
 * related to the values of the derived Measure and its derivatives, 
 * all of which require cache entries of type T.
 *
 * The constructor needs to be told how many type-T cache entries
 * to allocate. 
 */
template <class T>
class Measure_<T>::Implementation : public AbstractMeasure::Implementation {
public:
    const T& getValue(const State& s, int derivOrder) const {
        SimTK_ERRCHK2(0 <= derivOrder && derivOrder <= getNumTimeDerivatives(),
            "Measure_<T>::getValue()",
            "derivOrder %d was out of range; this Measure allows 0-%d.",
            derivOrder, getNumTimeDerivatives()); 

        SimTK_ERRCHK2
            (   getDependsOnStage(derivOrder)==Stage::Empty
             || (isInSubsystem() 
                 && getStage(s)>=getDependsOnStage(derivOrder))
             || (!isInSubsystem() 
                 && s.getSystemStage()>=getDependsOnStage(derivOrder)),
            "Measure_<T>::getValue()",
            "Expected State to have been realized to at least stage "
            "%s but stage was %s.", 
            getDependsOnStage(derivOrder).getName().c_str(), 
            (isInSubsystem() ? getStage(s) : s.getSystemStage())
                .getName().c_str());

        if (derivOrder < getNumCacheEntries()) {
            if (!isCacheValueRealized(s,derivOrder)) {
                T& value = updCacheEntry(s,derivOrder);
                calcCachedValueVirtual(s, derivOrder, value);
                markCacheValueRealized(s,derivOrder);
                return value;
            }
            return getCacheEntry(s,derivOrder);
        }

        // We can't handle it here -- punt to the concrete Measure
        // for higher order derivatives.
        return getUncachedValueVirtual(s,derivOrder); 
    }

protected:
    // numValues is one greater than the number of derivatives; i.e., there
    // is room for the value ("0th" derivative) also. The default is to
    // allocate just room for the value.
    explicit Implementation(int numValues=1) : derivIx(numValues) {}

    // Satisfy the realizeTopology() pure virtual here now that we know
    // the data type T.
    void realizeTopology(State& s) const {
        Implementation* mutableThis = const_cast<Implementation*>(this);
        // Allocate cache entries.
        for (int i=0; i < getNumCacheEntries(); ++i) {
            const Stage dependsOn = getDependsOnStage(i);
            mutableThis->derivIx[i] =
               this->getSubsystem().allocateLazyCacheEntry
                    (s, dependsOn, new Value<T>());
        }

        // Call the concrete class virtual if any.
        realizeMeasureTopologyVirtual(s);
    }

    int getNumCacheEntries() const {return (int)derivIx.size();}

    const T& getCacheEntry(const State& s, int derivOrder) const {
        SimTK_ERRCHK2(0 <= derivOrder && derivOrder < getNumCacheEntries(),
            "Measure_<T>::Implementation::getCacheEntry()",
            "Derivative order %d is out of range; only %d cache entries"
            " were allocated.", derivOrder, getNumCacheEntries());

        return Value<T>::downcast(
            this->getSubsystem().getCacheEntry(s, derivIx[derivOrder]));
    }

    T& updCacheEntry(const State& s, int derivOrder) const {
        SimTK_ERRCHK2(0 <= derivOrder && derivOrder < getNumCacheEntries(),
            "Measure_<T>::Implementation::updCacheEntry()",
            "Derivative order %d is out of range; only %d cache entries"
            " were allocated.", derivOrder, getNumCacheEntries());

        return Value<T>::updDowncast(
            this->getSubsystem().updCacheEntry(s, derivIx[derivOrder]));
    }

    bool isCacheValueRealized(const State& s, int derivOrder) const {
        SimTK_ERRCHK2(0 <= derivOrder && derivOrder < getNumCacheEntries(),
            "Measure_<T>::Implementation::isCacheValueRealized()",
            "Derivative order %d is out of range; only %d cache entries"
            " were allocated.", derivOrder, getNumCacheEntries());

        return this->getSubsystem().isCacheValueRealized(s, derivIx[derivOrder]);
    }

    void markCacheValueRealized(const State& s, int derivOrder) const {
        SimTK_ERRCHK2(0 <= derivOrder && derivOrder < getNumCacheEntries(),
            "Measure_<T>::Implementation::markCacheValueRealized()",
            "Derivative order %d is out of range; only %d cache entries"
            " were allocated.", derivOrder, getNumCacheEntries());

        this->getSubsystem().markCacheValueRealized(s, derivIx[derivOrder]);
    }

    // VIRTUALS //
    // Ordinals must retain the same meaning from release to release
    // to preserve binary compatibility.

    /* 0*/virtual void realizeMeasureTopologyVirtual(State&) const {}
    /* 1*/virtual void 
    calcCachedValueVirtual(const State&, int derivOrder, T& value) const
    {   SimTK_ERRCHK1_ALWAYS(!"implemented", 
        "Measure_<T>::Implementation::calcCachedValueVirtual()",
        "This method should have been overridden by the derived"
        " Measure but was not. It is needed to calculate the"
        " cached value for derivOrder=%d.", derivOrder); }

    // This is only called when derivOrder >= the number of cache 
    // entries we have, but still <= the number of derivatives the
    // Measure says it can deliver.
    /* 2*/virtual const T& 
    getUncachedValueVirtual(const State&, int derivOrder) const
    {   SimTK_ERRCHK1_ALWAYS(!"implemented", 
            "Measure_<T>::Implementation::getUncachedValueVirtual()",
            "This method should have been overridden by the derived"
            " Measure but was not. It is needed to return the uncached"
            " value at derivOrder=%d.", derivOrder);
        return *reinterpret_cast<T*>(0);
    }

    // STATICS //
    static const T& getValueZero() {
        static T zero(0);
        return zero;
    }

    static const T& getValueOne() {
        static T one(1);
        return one;
    }

private:
    Array_<CacheEntryIndex> derivIx;
};



    //////////////////////////////
    // CONSTANT::IMPLEMENTATION //
    //////////////////////////////

template <class T>
class Measure_<T>::Constant::Implementation 
:   public Measure_<T>::Implementation 
{
public:
    // We don't want the base class to allocate *any* cache entries.
    Implementation() : Measure_<T>::Implementation(0) {}
    explicit Implementation(const T& value) 
    :   Measure_<T>::Implementation(0), value(value) {}

    // Allow this to be overridden by derived classes so they can
    // refuse to allow their values to change.
    virtual void setValue(const T& v) {
        value = v;
        this->invalidateTopologyCache();
    }

    // Implementations of virtual methods.
    // Measure_<T> virtuals:
    // No cached values.

    const T& getUncachedValueVirtual(const State&, int derivOrder) const 
    {   return derivOrder>0 ? this->getValueZero() : value; }

    // AbstractMeasure virtuals:
    Implementation* cloneVirtual() const {return new Implementation(*this);}
    Stage getDependsOnStageVirtual(int derivOrder) const 
    {   return derivOrder>0 ? Stage::Empty : Stage::Topology; }
    int getNumTimeDerivativesVirtual() const 
    {   return std::numeric_limits<int>::max(); }

private:
    T value;
};

    //////////////////////////
    // ZERO::IMPLEMENTATION //
    //////////////////////////

template <class T>
class Measure_<T>::Zero::Implementation 
:   public Constant::Implementation 
{
public:
    Implementation() : Constant::Implementation(T(0)) {} 

    // VIRTUALS
    // From Constant::Implementation:
    void setValue(const T&) {
        SimTK_ERRCHK_ALWAYS(!"invalid", "Measure_<T>::Zero::setValue()",
            "You can't change the value of Zero!");
    }

    // From AbstractMeasure:
    Implementation* cloneVirtual() const {return new Implementation(*this);}
    Stage getDependsOnStageVirtual(int) const 
    {   return Stage::Empty; }
};

    /////////////////////////
    // ONE::IMPLEMENTATION //
    /////////////////////////


template <class T>
class Measure_<T>::One::Implementation 
:   public Constant::Implementation 
{
public:
    Implementation() : Constant::Implementation(T(1)) {} 

    // VIRTUALS
    // From Constant::Implementation:
    void setValue(const T&) {
        SimTK_ERRCHK_ALWAYS(!"invalid", "Measure_<T>::One::setValue()",
            "You can't change the value of One!");
    }
    // From AbstractMeasure:
    Implementation* cloneVirtual() const {return new Implementation(*this);}
    Stage getDependsOnStageVirtual(int) const 
    {   return Stage::Empty; }
};


    //////////////////////////
    // TIME::IMPLEMENTATION //
    //////////////////////////

template <class T>
class Measure_<T>::Time::Implementation {};

template <>
class Measure_<Real>::Time::Implementation 
:   public Measure_<Real>::Implementation 
{
public:
    // We don't want the base class to allocate *any* cache entries.
    Implementation() : Measure_<Real>::Implementation(0) {}

    // Implementations of virtual methods.
    // Measure_<Real> virtuals:
    // No cached values.

    const Real& getUncachedValueVirtual(const State& s, int derivOrder) const 
    {   return derivOrder==0 ? s.getTime()
            : (derivOrder==1 ? this->getValueOne() 
                             : this->getValueZero()); }

    // AbstractMeasure virtuals:
    Implementation* cloneVirtual() const {return new Implementation(*this);}
    Stage getDependsOnStageVirtual(int derivOrder) const 
    {   return derivOrder>0 ? Stage::Empty : Stage::Time; }

    // Value is t, 1st derivative is 1, the rest are 0.
    int getNumTimeDerivativesVirtual() const 
    {   return std::numeric_limits<int>::max(); }
};

    //////////////////////////////
    // VARIABLE::IMPLEMENTATION //
    //////////////////////////////

template <class T>
class Measure_<T>::Variable::Implementation 
:   public Measure_<T>::Implementation 
{
public:
    // We don't want the base class to allocate *any* cache entries;
    // we'll use the variable as its own value and zeroes for all
    // the derivatives.
    Implementation() : Measure_<T>::Implementation(0) {}

    Implementation(Stage invalidates, const T& defaultValue) 
    :   Measure_<T>::Implementation(0),
        invalidatedStage(invalidates), defaultValue(defaultValue) {}

    // Copy constructor should not copy the variable.
    Implementation(const Implementation& source)
    :   Measure_<T>::Implementation(0),
        invalidatedStage(source.invalidatedStage), 
        defaultValue(source.defaultValue) {}

    void setDefaultValue(const T& v) {
        defaultValue = v;
        this->invalidateTopologyCache();
    }

    void setInvalidatedStage(Stage invalidates) {
        invalidatedStage = invalidates;
        this->invalidateTopologyCache();
    }

    const T& getDefaultValue()     const {return defaultValue;}
    Stage    getInvalidatedStage() const {return invalidatedStage;}

    void setValue(State& state, const T& value) const 
    {   updVarValue(state) = value; }

    // Implementations of virtual methods.
    Implementation* cloneVirtual() const {return new Implementation(*this);}

    int getNumTimeDerivativesVirtual() const 
    {   return std::numeric_limits<int>::max(); }

    // Discrete variable is available after Model stage; but all its 
    // derivatives are zero so are always available.
    Stage getDependsOnStageVirtual(int derivOrder) const 
    {   return derivOrder>0 ? Stage::Empty : Stage::Model;}

    const T& getUncachedValueVirtual(const State& s, int derivOrder) const 
    {   return derivOrder>0 ? this->getValueZero() : getVarValue(s); }

    // No cached values.

    void realizeMeasureTopologyVirtual(State& s) const {
        discreteVarIndex = this->getSubsystem().allocateDiscreteVariable
            (s, invalidatedStage, new Value<T>(defaultValue));
    }
private:
    const T& getVarValue(const State& s) const {
        assert(discreteVarIndex.isValid());
        return Value<T>::downcast(
            this->getSubsystem().getDiscreteVariable(s, discreteVarIndex));
    }
    T& updVarValue(State& s) const {
        assert(discreteVarIndex.isValid());
        return Value<T>::downcast(
            this->getSubsystem().updDiscreteVariable(s, discreteVarIndex));
    }

    // TOPOLOGY STATE
    Stage   invalidatedStage; // TODO this shouldn't be needed
    T       defaultValue;

    // TOPOLOGY CACHE
    mutable DiscreteVariableIndex discreteVarIndex;
};

    //////////////////////////////
    // SINUSOID::IMPLEMENTATION //
    //////////////////////////////

template <class T>
class Measure_<T>::Sinusoid::Implementation
:   public Measure_<T>::Implementation 
{
    static const int NumDerivs = 3;
public:
    Implementation() 
    :   Measure_<T>::Implementation(NumDerivs+1),
        a(CNT<T>::getNaN()), w(CNT<T>::getNaN()), p(CNT<T>::getNaN()) {}

    Implementation(const T& amplitude, 
                   const T& frequency, 
                   const T& phase=T(0))
    :   Measure_<T>::Implementation(NumDerivs+1),
        a(amplitude), w(frequency), p(phase) {}

    // Default copy constructor is fine.

    // Implementations of virtual methods.
    Implementation* cloneVirtual() const {return new Implementation(*this);}

    int getNumTimeDerivativesVirtual() const {return NumDerivs;}

    Stage getDependsOnStageVirtual(int order) const 
    {   return Stage::Time; }

    void calcCachedValueVirtual(const State& s, int derivOrder, T& value) const {
        // We need to allow the compiler to select std::sin or SimTK::sin
        // based on the argument type.
        using std::sin; using std::cos;

        assert(NumDerivs == 3);
        const Real t = s.getTime();
        const T arg = w*t + p;

        switch (derivOrder) {
        case 0: value =        a*sin(arg); break;
        case 1: value =      w*a*cos(arg); break;
        case 2: value =   -w*w*a*sin(arg); break;
        case 3: value = -w*w*w*a*cos(arg); break;
        default: SimTK_ASSERT1_ALWAYS(!"out of range",
                     "Measure::Sinusoid::Implementation::calcCachedValueVirtual():"
                     " derivOrder %d is out of range 0-3.", derivOrder);
        }
    }

    // There are no uncached values.

private:
    // TOPOLOGY STATE
    T a, w, p;

    // TOPOLOGY CACHE
    // nothing
};

    //////////////////////////
    // PLUS::IMPLEMENTATION //
    //////////////////////////

template <class T>
class Measure_<T>::Plus::Implementation
:   public Measure_<T>::Implementation 
{
public:
    // TODO: Currently allocates just one cache entry.
    // left and right will be empty handles.
    Implementation() {}

    Implementation(const Measure_<T>& left, 
                   const Measure_<T>& right)
    :   left(left), right(right) {}

    // Default copy constructor gives us a new Implementation object,
    // but with references to the *same* operand measures.

    // Implementations of virtual methods.

    // This uses the default copy constructor.
    Implementation* cloneVirtual() const 
    {   return new Implementation(*this); }

    // TODO: Let this be settable up to the min number of derivatives 
    // provided by the arguments.
    int getNumTimeDerivativesVirtual() const {return 0;} 
    //{   return std::min(left.getNumTimeDerivatives(), 
    //                    right.getNumTimeDerivatives()); }

    Stage getDependsOnStageVirtual(int order) const 
    {   return Stage(std::max(left.getDependsOnStage(order),
                              right.getDependsOnStage(order))); }


    void calcCachedValueVirtual(const State& s, int derivOrder, T& value) const {
        value = left.getValue(s,derivOrder) + right.getValue(s,derivOrder);
    }

    // There are no uncached values.

private:
    // TOPOLOGY STATE
    const Measure_<T> left;
    const Measure_<T> right;

    // TOPOLOGY CACHE
    // nothing
};

    ///////////////////////////
    // MINUS::IMPLEMENTATION //
    ///////////////////////////

template <class T>
class Measure_<T>::Minus::Implementation
:   public Measure_<T>::Implementation 
{
public:
    // TODO: Currently allocates just one cache entry.
    // left and right will be empty handles.
    Implementation() {}

    Implementation(const Measure_<T>& left, 
                   const Measure_<T>& right)
    :   left(left), right(right) {}

    // Default copy constructor gives us a new Implementation object,
    // but with references to the *same* operand measures.

    // Implementations of virtual methods.

    // This uses the default copy constructor.
    Implementation* cloneVirtual() const 
    {   return new Implementation(*this); }

    // TODO: Let this be settable up to the min number of derivatives 
    // provided by the arguments.
    int getNumTimeDerivativesVirtual() const {return 0;} 
    //{   return std::min(left.getNumTimeDerivatives(), 
    //                    right.getNumTimeDerivatives()); }

    Stage getDependsOnStageVirtual(int order) const 
    {   return Stage(std::max(left.getDependsOnStage(order),
                              right.getDependsOnStage(order))); }


    void calcCachedValueVirtual(const State& s, int derivOrder, T& value) const {
        value = left.getValue(s,derivOrder) - right.getValue(s,derivOrder);
    }

    // There are no uncached values.

private:
    // TOPOLOGY STATE
    const Measure_<T> left;
    const Measure_<T> right;

    // TOPOLOGY CACHE
    // nothing
};


    ///////////////////////////
    // SCALE::IMPLEMENTATION //
    ///////////////////////////

template <class T>
class Measure_<T>::Scale::Implementation
:   public Measure_<T>::Implementation 
{
public:
    // TODO: Currently allocates just one cache entry.
    // scale will be uninitialized, operand will be empty handle.
    Implementation() : factor(NaN) {}

    Implementation(Real factor, const Measure_<T>& operand)
    :   factor(factor), operand(operand) {}

    // Default copy constructor gives us a new Implementation object,
    // but with references to the *same* operand measure.

    void setScaleFactor(Real sf) {
        factor = sf;
        this->invalidateTopologyCache();
    }

    // Implementations of virtual methods.

    // This uses the default copy constructor.
    Implementation* cloneVirtual() const 
    {   return new Implementation(*this); }

    // TODO: Let this be settable up to the min number of derivatives 
    // provided by the arguments.
    int getNumTimeDerivativesVirtual() const {return 0;} 
    //{   return std::min(left.getNumTimeDerivatives(), 
    //                    right.getNumTimeDerivatives()); }

    Stage getDependsOnStageVirtual(int order) const 
    {   return operand.getDependsOnStage(order); }


    void calcCachedValueVirtual(const State& s, int derivOrder, T& value) const {
        value = factor * operand.getValue(s,derivOrder);
    }

    // There are no uncached values.

private:
    // TOPOLOGY STATE
    const Real        factor;
    const Measure_<T> operand;

    // TOPOLOGY CACHE
    // nothing
};

    ///////////////////////////////
    // INTEGRATE::IMPLEMENTATION //
    ///////////////////////////////

template <class T>
class Measure_<T>::Integrate::Implementation 
:   public Measure_<T>::Implementation {
public:
    // We don't want any cache entries allocated -- we'll use the
    // State z variable as its own value and let the derivative 
    // Measure supply its value.

    // The derivative and initialConditions Measures will be
    // empty handles if this is default constructed.
    Implementation() : Measure_<T>::Implementation(0) {}

    // Here we're shallow-copying the Measure handles so we'll
    // be referring to the original Measures.
    Implementation(const Measure_<T>& deriv, const Measure_<T>& ic)
    :   Measure_<T>::Implementation(0), 
        derivMeasure(deriv), icMeasure(ic) {}

    // Copy constructor shallow-copies the referenced measures, but
    // we don't want to share our state variable.
    Implementation(const Implementation& source)
    :   Measure_<T>::Implementation(0), 
        derivMeasure(source.derivMeasure), icMeasure(source.icMeasure) {}

    void setValue(State& s, const T& value) const
    {   assert(zIndex >= 0); 
        this->getSubsystem().updZ(s)[zIndex] = value; }
    
    const Measure_<T>& getDerivativeMeasure() const
    {   SimTK_ERRCHK(!derivMeasure.isEmptyHandle(), 
            "Measure_<T>::Integrate::getDerivativeMeasure()",
            "No derivative measure is available for this integrated measure."); 
        return derivMeasure; }

    const Measure_<T>& getInitialConditionMeasure() const
    {   SimTK_ERRCHK(!icMeasure.isEmptyHandle(), 
            "Measure_<T>::Integrate::getInitialConditionMeasure()",
            "No initial condition measure is available for this "
            "integrated measure."); 
        return icMeasure; }

    void setDerivativeMeasure(const Measure_<T>& d)
    {   derivMeasure = d; this->invalidateTopologyCache(); }
    void setInitialConditionMeasure(const Measure_<T>& ic)
    {   icMeasure = ic; this->invalidateTopologyCache(); }

    // Implementations of virtuals.

    // This uses the copy constructor defined above.
    Implementation* cloneVirtual() const 
    {   return new Implementation(*this); }

    int getNumTimeDerivativesVirtual() const {return 1;}

    // There are no cached values.

    const T& getUncachedValueVirtual(const State& s, int derivOrder) const
    {   if (derivOrder>0) 
            return getDerivativeMeasure().getValue(s);
        else { 
            assert(zIndex.isValid()); 
            return this->getSubsystem().getZ(s)[zIndex];
        }
    }
    Stage getDependsOnStageVirtual(int derivOrder) const 
    {   return derivOrder>0 ? getDerivativeMeasure().getDependsOnStage(0)
                            : Stage::Time; }

    void initializeVirtual(State& s) const {
        assert(zIndex.isValid());
        Real& z = this->getSubsystem().updZ(s)[zIndex];
        if (!icMeasure.isEmptyHandle()) 
             z = icMeasure.getValue(s);
        else z = 0;
    }

    void realizeMeasureTopologyVirtual(State& s) const {
        static const Vector zero(1, Real(0));
        zIndex = this->getSubsystem().allocateZ(s, zero);
    }

    void realizeMeasureAccelerationVirtual(const State& s) const {
        assert(zIndex.isValid());
        Real& zdot = this->getSubsystem().updZDot(s)[zIndex];
        if (!derivMeasure.isEmptyHandle()) 
             zdot = derivMeasure.getValue(s);
        else zdot = 0;
    }

private:
    // TOPOLOGY STATE
    const Measure_<T> derivMeasure; // just handles
    const Measure_<T> icMeasure;

    // TOPOLOGY CACHE
    mutable ZIndex zIndex;
};

    ///////////////////////////////////
    // DIFFERENTIATE::IMPLEMENTATION //
    ///////////////////////////////////

// This helper class is the contents of the discrete state variable and corresponding
// cache entry maintained by this measure. The variable is auto-update,
// meaning the value of the cache entry replaces the state variable at the
// start of each step.
// TODO: This was a local class in Measure_<T>::Differentiate::Implementation
// but VC++ 8 (2005) failed to properly instantiate the templatized operator<<()
// in that case; doing it this way is a workaround.
template <class T>
class Measure_Differentiate_Result {
public:
    Measure_Differentiate_Result() : derivIsGood(false) {}
    T       operand;    // previous value of operand
    T       operandDot; // previous value of derivative
    bool    derivIsGood; // do we think the deriv is a good one?
};


// Dummy for Value<Measure_Differentiate_Result>.
template <class T> inline std::ostream& 
operator<<(std::ostream& o, 
           const Measure_Differentiate_Result<T>&)
{   assert(!"not implemented"); return o; }

template <class T>
class Measure_<T>::Differentiate::Implementation
:   public Measure_<T>::Implementation 
{
    typedef Measure_Differentiate_Result<T> Result;
public:
    // Don't allocate any cache entries in the base class.
    Implementation() : Measure_<T>::Implementation(0) {}

    Implementation(const Measure_<T>& operand)
    :   Measure_<T>::Implementation(0),
        operand(operand), forceUseApprox(false), isApproxInUse(false) {}

    // Default copy constructor gives us a new Implementation object,
    // but with reference to the *same* operand measure.

    void setForceUseApproximation(bool mustApproximate) {
        forceUseApprox = mustApproximate;
        this->invalidateTopologyCache();
    }

    void setOperandMeasure(const Measure_<T>& operand) {
        this->operand = operand;
        this->invalidateTopologyCache();
    }

    bool getForceUseApproximation() const {return forceUseApprox;}
    bool isUsingApproximation() const {return isApproxInUse;}
    const Measure_<T>& getOperandMeasure() const {return operand;}

    // Implementations of virtual methods.

    // This uses the default copy constructor.
    Implementation* cloneVirtual() const 
    {   return new Implementation(*this); }

    // This has one fewer than the operand.
    int getNumTimeDerivativesVirtual() const
    {   if (!isApproxInUse) return operand.getNumTimeDerivatives()-1;
        else return 0; }

    Stage getDependsOnStageVirtual(int order) const 
    {   if (!isApproxInUse) return operand.getDependsOnStage(order+1);
        else return operand.getDependsOnStage(order); }


    // We're not using the Measure_<T> base class cache services, but
    // we do have one of our own. It looks uncached from the base class
    // point of view which is why we're implementing it here.
    const T& getUncachedValueVirtual(const State& s, int derivOrder) const
    {   if (!isApproxInUse) 
            return operand.getValue(s, derivOrder+1);

        ensureDerivativeIsRealized(s);
        const Subsystem& subsys = this->getSubsystem();
        const Result& result = Value<Result>::downcast
                                (subsys.getDiscreteVarUpdateValue(s,resultIx));
        return result.operandDot; // has a value but might not be a good one
    }

    void initializeVirtual(State& s) const {
        if (!isApproxInUse) return;

        assert(resultIx.isValid());
        const Subsystem& subsys = this->getSubsystem();
        Result& result = Value<Result>::updDowncast
                            (subsys.updDiscreteVariable(s,resultIx));
        result.operand = operand.getValue(s);
        result.operandDot = this->getValueZero();
        result.derivIsGood = false;
    }

    void realizeMeasureTopologyVirtual(State& s) const {
        isApproxInUse = (forceUseApprox || operand.getNumTimeDerivatives()==0);
        if (!isApproxInUse)
            return;

        resultIx = this->getSubsystem()
            .allocateAutoUpdateDiscreteVariable(s, operand.getDependsOnStage(0),
                new Value<Result>(), operand.getDependsOnStage(0));
    }

    void ensureDerivativeIsRealized(const State& s) const {
        assert(resultIx.isValid());
        const Subsystem& subsys = this->getSubsystem();
        if (subsys.isDiscreteVarUpdateValueRealized(s,resultIx))
            return;

        const Real t0 = subsys.getDiscreteVarLastUpdateTime(s,resultIx);
        const Result& prevResult = Value<Result>::downcast
           (subsys.getDiscreteVariable(s,resultIx));
        const T&   f0         = prevResult.operand;
        const T&   fdot0      = prevResult.operandDot;   // may be invalid
        const bool good0     = prevResult.derivIsGood;

        const Real t  = s.getTime();
        Result& result = Value<Result>::updDowncast
           (subsys.updDiscreteVarUpdateValue(s,resultIx));
        T&         f          = result.operand;          // renaming
        T&         fdot       = result.operandDot;
        bool&      good       = result.derivIsGood;

        f = operand.getValue(s);
        good = false;
        if (!isFinite(t0))
            fdot = this->getValueZero(); 
        else if (t == t0) {
            fdot = fdot0;
            good = good0;
        } else {
            fdot = (f-f0)/(t-t0); // 1st order
            if (good0)
                fdot = Real(2)*fdot - fdot0; // now 2nd order
            good = true; // either 1st or 2nd order estimate
        }
        subsys.markDiscreteVarUpdateValueRealized(s,resultIx);
    }
private:
    // TOPOLOGY STATE
    Measure_<T>     operand;
    bool            forceUseApprox;

    // TOPOLOGY CACHE
    mutable bool                    isApproxInUse;
    mutable DiscreteVariableIndex   resultIx;    // auto-update
};


} // namespace SimTK




#endif // SimTK_SimTKCOMMON_MEASURE_IMPLEMENTATION_H_
