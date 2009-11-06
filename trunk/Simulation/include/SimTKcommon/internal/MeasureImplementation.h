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
    :   measureName(name), mySubsystem(0), refCount(0) {}

    Implementation(const Implementation& src)
    :   measureName(src.measureName), mySubsystem(0), refCount(0) {}
        
    Implementation& operator=(const Implementation& src) {
        if (&src != this)
        {   measureName=src.measureName; refCount=0; mySubsystem=0; }
        return *this; 
    }

    // destructor is virtual

    // Increment the reference count and return its new value.
    int incrRefCount() const {return ++refCount;}

    // Decrement the reference count and return its new value.
    int decrRefCount() const {return --refCount;}

    const std::string& getName() const {return measureName;}

    Implementation* clone() const {return cloneVirtual();}

    void realizeTopology    (State& s)       const {realizeMeasureTopologyVirtual(s);}
    void realizeModel       (State& s)       const {realizeMeasureModelVirtual(s);}
    void realizeInstance    (const State& s) const {realizeMeasureInstanceVirtual(s);}
    void realizeTime        (const State& s) const {realizeMeasureTimeVirtual(s);}
    void realizePosition    (const State& s) const {realizeMeasurePositionVirtual(s);}
    void realizeVelocity    (const State& s) const {realizeMeasureVelocityVirtual(s);}
    void realizeDynamics    (const State& s) const {realizeMeasureDynamicsVirtual(s);}
    void realizeAcceleration(const State& s) const {realizeMeasureAccelerationVirtual(s);}
    void realizeReport      (const State& s) const {realizeMeasureReportVirtual(s);}

    void  initialize(State& s)               const {initializeVirtual(s);}

    int getNumTimeDerivatives() const {return getNumTimeDerivativesVirtual();}

    Stage getDependsOnStage(const State& s, int derivOrder) const {
        SimTK_ERRCHK2(0 <= derivOrder && derivOrder <= getNumTimeDerivatives(),
            "Measure::getDependsOnStage()",
            "derivOrder %d was out of range; this Measure allows 0-%d.",
            derivOrder, getNumTimeDerivatives()); 
        return getDependsOnStageVirtual(s,derivOrder); 
    }


    void setSubsystem(Subsystem& sub, MeasureIndex mx) 
    {   assert(!mySubsystem && mx.isValid()); 
        mySubsystem = &sub; myIndex = mx; }

    bool isInSubsystem() const {return mySubsystem != 0;}
    const Subsystem& getSubsystem() const {assert(mySubsystem); return *mySubsystem;}
    Subsystem& updSubsystem() {assert(mySubsystem); return *mySubsystem;}
    MeasureIndex getSubsystemMeasureIndex() const {assert(mySubsystem); return myIndex;}
    SubsystemIndex getSubsystemIndex() const
    {   return getSubsystem().getMySubsystemIndex(); }
    void invalidateTopologyCache() const
    {   if (isInSubsystem()) getSubsystem().invalidateSubsystemTopologyCache(); }

    Stage getStage(const State& s) const {return getSubsystem().getStage(s);}

    // VIRTUALS //
    // Ordinals must retain the same meaning from release to release
    // to preserve binary compatibility.

    /* 0*/virtual ~Implementation() {}
    /* 1*/virtual Implementation* cloneVirtual() const = 0;

    /* 2*/virtual void realizeMeasureTopologyVirtual(State&) const {}
    /* 3*/virtual void realizeMeasureModelVirtual(State&) const {}
    /* 4*/virtual void realizeMeasureInstanceVirtual(const State&) const {}
    /* 5*/virtual void realizeMeasureTimeVirtual(const State&) const {}
    /* 6*/virtual void realizeMeasurePositionVirtual(const State&) const {}
    /* 7*/virtual void realizeMeasureVelocityVirtual(const State&) const {}
    /* 8*/virtual void realizeMeasureDynamicsVirtual(const State&) const {}
    /* 9*/virtual void realizeMeasureAccelerationVirtual(const State&) const {}
    /*10*/virtual void realizeMeasureReportVirtual(const State&) const {}

    /*11*/virtual void initializeVirtual(State&) const {}
    /*12*/virtual int  getNumTimeDerivativesVirtual() const {return 0;}

    /*13*/virtual Stage 
            getDependsOnStageVirtual(const State&, int order) const = 0;

private:
    std::string     measureName;

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

inline AbstractMeasure::AbstractMeasure(const AbstractMeasure& src) 
:   impl(0) {
    if (src.impl) {
        impl = src.impl->mySubsystem ? src.impl : src.impl->clone();
        impl->incrRefCount();
    }
}

inline AbstractMeasure& AbstractMeasure::
operator=(const AbstractMeasure& src) {
    if (&src != this) {
        if (impl && impl->decrRefCount()==0) delete impl;
        impl = src.impl->mySubsystem ? src.impl : src.impl->clone();
        impl->incrRefCount();
    }
    return *this;
}

inline AbstractMeasure::
~AbstractMeasure()
{   if (impl && getImpl().decrRefCount()==0) delete impl;}

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
getDependsOnStage(const State& s, int derivOrder) const
{   return getImpl().getDependsOnStage(s,derivOrder); }



    /////////////////////////////////
    // MEASURE_<T>::IMPLEMENTATION //
    /////////////////////////////////

/**
 * This is the base Implementation class for all Measures whose value type 
 * is known. Logically this is an abstract base class. However, you can 
 * create an empty Measure<T> handle that can be assigned to any derived 
 * concrete Measure.
 */
template <class T>
class Measure_<T>::Implementation : public AbstractMeasure::Implementation {
public:
    // Has implicit default constructor.

    const T& getValue(const State& s, int derivOrder) const {
        SimTK_ERRCHK2(0 <= derivOrder && derivOrder <= getNumTimeDerivatives(),
            "Measure_<T>::getValue()",
            "derivOrder %d was out of range; this Measure allows 0-%d.",
            derivOrder, getNumTimeDerivatives()); 

        SimTK_ERRCHK2
            (   getDependsOnStage(s,derivOrder)==Stage::Empty
             || (getStage(s)>=getDependsOnStage(s,derivOrder)),
            "Measure_<T>::getValue()",
            "Expected State to have been realized to at least stage "
            "%s but stage was %s.", 
            getDependsOnStage(s,derivOrder).getName().c_str(), 
            getStage(s).getName().c_str());

        return getValueVirtual(s,derivOrder); 
    }

    // VIRTUALS //
    // Ordinals must retain the same meaning from release to release
    // to preserve binary compatibility.

    /* 0*/virtual const T& 
    getValueVirtual(const State&, int derivOrder) const = 0;

    // STATICS //
    static const T& getValueZero() {
        static T zero(0);
        return zero;
    }

    static const Value<T>& getValueOne() {
        static T one(1);
        return one;
    }
};

    //////////////////////////////
    // CONSTANT::IMPLEMENTATION //
    //////////////////////////////

template <class T>
class Measure_<T>::Constant::Implementation 
:   public Measure_<T>::Implementation 
{
public:
    Implementation() {}
    explicit Implementation(const T& value) : value(value) {}

    // Allow this to be overridden by derived classes so they can
    // refuse to allow their values to change.
    virtual void setValue(const T& v) {
        value = v;
        this->invalidateTopologyCache();
    }

    // Implementations of virtual methods.
    // Measure_<T> virtuals:
    const T& getValueVirtual(const State&, int derivOrder) const 
    {   return derivOrder ? this->getValueZero() : value; }

    // AbstractMeasure virtuals:
    Implementation* cloneVirtual() const {return new Implementation(*this);}
    Stage getDependsOnStageVirtual(const State&, int) const 
    {   return Stage::Topology; }

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
    Stage getDependsOnStageVirtual(const State&, int) const 
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
    Stage getDependsOnStageVirtual(const State&, int) const 
    {   return Stage::Empty; }
};

    //////////////////////////////
    // VARIABLE::IMPLEMENTATION //
    //////////////////////////////

template <class T>
class Measure_<T>::Variable::Implementation 
:   public Measure_<T>::Implementation 
{
public:
    Implementation(Stage invalidates, const T& defaultValue) 
    : invalidatedStage(invalidates), defaultValue(defaultValue) {}

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
    Stage getDependsOnStageVirtual(const State&, int order) const 
    {   return order ? Stage::Topology : Stage::Model;}

    const T& getValueVirtual(const State& s, int order) const 
    {   return order ? this->getValueZero() : getVarValue(s); }

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
    typedef Vec<NumDerivs+1,T>  ValType;     // cache entry type
    typedef Value<ValType>      AbsValType;

public:
    Implementation(const T& amplitude, 
                   const T& frequency, 
                   const T& phase=T(0))
    :   a(amplitude), w(frequency), p(phase) {}

    // Implementations of virtual methods.
    Implementation* cloneVirtual() const {return new Implementation(*this);}

    int getNumTimeDerivativesVirtual() const {return NumDerivs;}

    const T& getValueVirtual(const State& s, int order) const
    {   ensureRealized(s);
        return getValueEntry(s)[order]; }

    Stage getDependsOnStageVirtual(const State&, int order) const 
    {   return Stage::Time; }

    void realizeMeasureTopologyVirtual(State& s) const
    {   cacheIndex = this->getSubsystem().allocateCacheEntry
           (s, Stage::Time, Stage::Infinity, 
            new AbsValType(ValType::getNaN())); }


private:
    const ValType& getValueEntry(const State& s) const {
        assert(cacheIndex.isValid());
        return AbsValType::downcast(
            this->getSubsystem().getCacheEntry(s, cacheIndex));
    }

    ValType& updValueEntry(const State& s) const {
        assert(cacheIndex.isValid());
        return AbsValType::updDowncast(
            this->getSubsystem().updCacheEntry(s, cacheIndex));
    }

    void ensureRealized(const State& s) const {
        // We need to allow the compiler to select std::sin or SimTK::sin
        // based on the argument type.
        using std::sin; using std::cos;
        assert(cacheIndex.isValid());

        // If already calculated since last time change we don't need
        // to do anything.
        if (this->getSubsystem().isCacheValueCurrent(s, cacheIndex))
            return;

        ValType& entry = updValueEntry(s);

        assert(NumDerivs == 3);
        const Real t = s.getTime();
        const T arg = w*t + p;
        const T as = a*sin(arg);
        const T ac = a*cos(arg);
        entry[0] =  as;        //      a sin wt+p
        entry[1] =  w*ac;      //    w a cos wt+p
        entry[2] = -w*w*as;    // -w^2 a sin wt+p
        entry[3] = -w*w*w*ac;  // -w^3 a cos wt+p

        // Mark this as having been realized.
        this->getSubsystem().markCacheValueRealized(s, cacheIndex);
    }

    // TOPOLOGY STATE
    T a, w, p;

    // TOPOLOGY CACHE
    mutable CacheEntryIndex cacheIndex;
};

    //////////////////////////
    // PLUS::IMPLEMENTATION //
    //////////////////////////

template <class T>
class Measure_<T>::Plus::Implementation
:   public Measure_<T>::Implementation 
{
public:
    Implementation(const Measure_<T>& left, 
                   const Measure_<T>& right)
    :   left(left), right(right) {}

    // Implementations of virtual methods.
    Implementation* cloneVirtual() const {return new Implementation(*this);}

    // TODO: Let this be settable up to the min number of derivatives 
    // provided by the arguments.
    int getNumTimeDerivativesVirtual() const {return 0;} 
    //{   return std::min(left.getNumTimeDerivatives(), 
    //                    right.getNumTimeDerivatives()); }

    const T& getValueVirtual(const State& s, int order) const
    {   ensureRealized(s);
        return getValueEntry(s); }

    Stage getDependsOnStageVirtual(const State& s, int order) const 
    {   return Stage(std::max(left.getDependsOnStage(s,order),
                              right.getDependsOnStage(s,order))); }

    void realizeMeasureTopologyVirtual(State& s) const
    {   cacheIndex = this->getSubsystem().allocateCacheEntry
           (s, getDependsOnStage(s,0), Stage::Infinity, 
            new Value<T>()); }

private:
    const T& getValueEntry(const State& s) const {
        assert(cacheIndex.isValid());
        return Value<T>::downcast(
            this->getSubsystem().getCacheEntry(s, cacheIndex));
    }

    T& updValueEntry(const State& s) const {
        assert(cacheIndex.isValid());
        return Value<T>::updDowncast(
            this->getSubsystem().updCacheEntry(s, cacheIndex));
    }

    void ensureRealized(const State& s) const {
        assert(cacheIndex.isValid());

        // If already calculated since last time change we don't need
        // to do anything.
        if (this->getSubsystem().isCacheValueCurrent(s, cacheIndex))
            return;

        T& entry = updValueEntry(s);

        entry = left.getValue(s) + right.getValue(s);

        // Mark this as having been realized.
        this->getSubsystem().markCacheValueRealized(s, cacheIndex);
    }

    // TOPOLOGY STATE
    const Measure_<T>&   left;
    const Measure_<T>&   right;

    // TOPOLOGY CACHE
    mutable CacheEntryIndex cacheIndex;
};

    ///////////////////////////////
    // INTEGRATE::IMPLEMENTATION //
    ///////////////////////////////

template <class T>
class Measure_<T>::Integrate::Implementation 
:   public Measure_<T>::Implementation {
public:
    Implementation() : derivMeasure(0), icMeasure(0) {}
    Implementation(const Measure_<T>& deriv, const Measure_<T>& ic)
    :   derivMeasure(&deriv), icMeasure(&ic) {}

    void setValue(State& s, const T& value) const
    {    assert(zIndex >= 0); this->getSubsystem().updZ(s)[zIndex] = value; }
    
    const Measure_<T>& getDerivativeMeasure() const
    {   SimTK_ERRCHK(derivMeasure, 
            "Measure_<T>::Integrate::getDerivativeMeasure()",
            "No derivative measure is available for this integrated measure."); 
        return *derivMeasure; }

    const Measure_<T>& getInitialConditionMeasure() const
    {   SimTK_ERRCHK(icMeasure, 
            "Measure_<T>::Integrate::getInitialConditionMeasure()",
            "No initial condition measure is available for this "
            "integrated measure."); 
        return *icMeasure; }

    void setDerivativeMeasure(const Measure_<T>& d)
    {   derivMeasure = &d; this->invalidateTopologyCache(); }
    void setInitialConditionMeasure(const Measure_<T>& ic)
    {   icMeasure = &ic; this->invalidateTopologyCache(); }

    // Implementations of virtuals.
    Implementation* cloneVirtual() const {return new Implementation(*this);}

    int getNumTimeDerivativesVirtual() const {return 1;}

    const T& getValueVirtual(const State& s, int order) const
    {   if (order) 
            return getDerivativeMeasure().getValue(s);
        else { 
            assert(zIndex.isValid()); 
            return this->getSubsystem().getZ(s)[zIndex];
        }
    }
    Stage getDependsOnStageVirtual(const State& s, int order) const 
    {   return order ? getDerivativeMeasure().getDependsOnStage(s,0)
                     : Stage::Time; }

    void initializeVirtual(State& s) const {
        assert(zIndex >= 0);
        Real& z = this->getSubsystem().updZ(s)[zIndex];
        if (icMeasure) z = icMeasure->getValue(s);
        else z = 0;
    }

    void realizeMeasureTopologyVirtual(State& s) const {
        static const Vector zero(1, Real(0));
        zIndex = this->getSubsystem().allocateZ(s, zero);
    }

    void realizeMeasureAccelerationVirtual(const State& s) const {
        assert(zIndex.isValid());
        Real& zdot = this->getSubsystem().updZDot(s)[zIndex];
        if (derivMeasure) zdot = derivMeasure->getValue(s);
        else zdot = 0;
    }

private:
    // TOPOLOGY STATE
    const Measure_<T>*   derivMeasure;
    const Measure_<T>*   icMeasure;

    // TOPOLOGY CACHE
    mutable ZIndex zIndex;
};


} // namespace SimTK

#endif // SimTK_SimTKCOMMON_MEASURE_IMPLEMENTATION_H_
