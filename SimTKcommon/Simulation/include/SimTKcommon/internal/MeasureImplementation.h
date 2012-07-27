#ifndef SimTK_SimTKCOMMON_MEASURE_IMPLEMENTATION_H_
#define SimTK_SimTKCOMMON_MEASURE_IMPLEMENTATION_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
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
#include "SimTKcommon/internal/Measure.h"
#include "SimTKcommon/internal/Subsystem.h"
#include "SimTKcommon/internal/System.h"
#include "SimTKcommon/internal/SubsystemGuts.h"

#include <cmath>


namespace SimTK {


//==============================================================================
//                    ABSTRACT MEASURE :: IMPLEMENTATION
//==============================================================================

/**
 * The abstract parent of all Measure Implementation classes.
 */
class SimTK_SimTKCOMMON_EXPORT AbstractMeasure::Implementation {
protected:
    /** This default constructor is for use by concrete measure implementation
    classes. **/
    Implementation() : copyNumber(0), mySubsystem(0), refCount(0) {}

    /** Base class copy constructor removes the Subsystem
    and sets the reference count to zero. This gets used by the clone()
    methods in the concrete classes. **/
    Implementation(const Implementation& src)
    :   copyNumber(src.copyNumber+1), mySubsystem(0), refCount(0) {}
    
    /** Base class copy assignment operator removes the
    Subsystem, and sets the reference count to zero. This is probably
    not used. **/
    Implementation& operator=(const Implementation& src) {
        if (&src != this)
        {   copyNumber=src.copyNumber+1;
            refCount=0; mySubsystem=0; }
        return *this; 
    }

    // destructor is virtual; below

    // Increment the reference count and return its new value.
    int incrRefCount() const {return ++refCount;}

    // Decrement the reference count and return its new value.
    int decrRefCount() const {return --refCount;}

    // Get the current value of the reference counter.
    int getRefCount() const {return refCount;}

    int           getCopyNumber()  const {return copyNumber;}

    /** This is a deep copy of the concrete Implementation object, except the
    Subsystem will have been removed. The reference count on the new object
    will be zero; be sure to increment it if you put it in a handle. **/
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

    /** This should be called at the start of a time stepping study to
    cause this %Measure to set its state variables (if any) in the supplied
    state to their initial conditions. **/
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

//==============================================================================
//                      ABSTRACT MEASURE DEFINITIONS
//==============================================================================
// These had to wait for AbstractMeasure::Implementation to be defined.

inline AbstractMeasure::
AbstractMeasure(Implementation* g) 
:   impl(g)
{   if (impl) impl->incrRefCount(); }

inline AbstractMeasure::
AbstractMeasure(Subsystem& sub, Implementation* g, const SetHandle&) 
:   impl(g) {
    SimTK_ERRCHK(hasImpl(), "AbstractMeasure::AbstractMeasure()",
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


/** @cond **/ // Hide from Doxygen.
// This is a helper class that makes it possible to treat Real, Vec, and 
// Vector objects uniformly.
template <class T> class Measure_Num {
};

template <> class Measure_Num<float> {
public:
    typedef float Element;
    static int size(const float&) {return 1;}
    static const float& get(const float& v, int i) {assert(i==0); return v;}
    static float& upd(float& v, int i) {assert(i==0); return v;}
    static void makeNaNLike(const float&, float& nanValue) 
    {   nanValue = CNT<float>::getNaN();}
    static void makeZeroLike(const float&, float& zeroValue) {zeroValue=0.f;}
};

template <> class Measure_Num<double> {
public:
    typedef double Element;
    static int size(const double&) {return 1;}
    static const double& get(const double& v, int i) {assert(i==0); return v;}
    static double& upd(double& v, int i) {assert(i==0); return v;}
    static void makeNaNLike(const double&, double& nanValue) 
    {  nanValue = CNT<double>::getNaN(); }
    static void makeZeroLike(const double&, double& zeroValue) {zeroValue=0.;}
};

// We only support stride 1 (densely packed) Vec types.
template <int M, class E>
class Measure_Num< Vec<M,E,1> > {
    typedef Vec<M,E,1> T;
public:
    typedef E Element;
    static int size(const T&) {return M;}
    static const E& get(const T& v, int i) {return v[i];}
    static E& upd(T& v, int i) {return v[i];}
    static void makeNaNLike (const T&, T& nanValue)  {nanValue.setToNaN();}
    static void makeZeroLike(const T&, T& zeroValue) {zeroValue.setToZero();}
};

// We only support column major (densely packed) Mat types.
template <int M, int N, class E>
class Measure_Num< Mat<M,N,E> > {
    typedef Mat<M,N,E> T;
public:
    typedef E Element;
    static int size(const T&) {return N;} // number of columns
    static const typename T::TCol& get(const T& m, int j) {return m.col(j);}
    static typename T::TCol& upd(T& m, int j) {return m.col(j);}
    static void makeNaNLike (const T&, T& nanValue)  {nanValue.setToNaN();}
    static void makeZeroLike(const T&, T& zeroValue) {zeroValue.setToZero();}
};


template <class E>
class Measure_Num< Vector_<E> > {
    typedef Vector_<E> T;
public:
    typedef E Element;
    static int size(const T& v) {return v.size();}
    static const E& get(const T& v, int i) {return v[i];}
    static E& upd(T& v, int i) {return v[i];}
    static void makeNaNLike(const T& v, T& nanValue) 
    {   nanValue.resize(v.size()); nanValue.setToNaN(); }
    static void makeZeroLike(const T& v, T& zeroValue)
    {   zeroValue.resize(v.size()); zeroValue.setToZero(); }

};

/** @endcond **/

//==============================================================================
//                       MEASURE_<T> :: IMPLEMENTATION
//==============================================================================
/** This is the base Implementation class for all Measures whose value type is 
known. This class is still abstract but provides many services related to the 
values of the derived Measure and its derivatives, all of which require cache 
entries of type T.

The constructor needs to be told how many type-T cache entries to allocate. **/ 
template <class T>
class Measure_<T>::Implementation : public AbstractMeasure::Implementation {
public:
    const T& getValue(const State& s, int derivOrder) const {
        SimTK_ERRCHK2(0 <= derivOrder && derivOrder <= getNumTimeDerivatives(),
            "Measure_<T>::getValue()",
            "derivOrder %d was out of range; this Measure allows 0-%d.",
            derivOrder, getNumTimeDerivatives()); 

        // We require the stage to have been advanced to at least the one
        // before this measure's depends-on stage since this will get called
        // towards the end of the depends-on stage realization.
        if (getDependsOnStage(derivOrder) != Stage::Empty) {
            Stage prevStage = getDependsOnStage(derivOrder).prev();

            SimTK_ERRCHK2
                (   ( isInSubsystem() && getStage(s)>=prevStage)
                 || (!isInSubsystem() && s.getSystemStage()>=prevStage),
                "Measure_<T>::getValue()",
                "Expected State to have been realized to at least stage "
                "%s but stage was %s.", 
                prevStage.getName().c_str(), 
                (isInSubsystem() ? getStage(s) : s.getSystemStage())
                    .getName().c_str());
        }

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

    /** Set a new default value for this %Measure. This is a topological 
    change. **/
    void setDefaultValue(const T& defaultValue) {
        this->defaultValue = defaultValue;
        Measure_Num<T>::makeZeroLike(defaultValue, zeroValue);
        this->invalidateTopologyCache(); 
    }

    /** Return a reference to the value that this %Measure will use to 
    initialize its value-level state resource (state variable or cache entry) 
    during the next call to realizeTopology(). **/
    const T& getDefaultValue() const {return defaultValue;}

    void setIsPresumedValidAtDependsOnStage(bool presume) 
    {   presumeValidAtDependsOnStage = presume;
        this->invalidateTopologyCache(); }

    bool getIsPresumedValidAtDependsOnStage() const 
    {   return presumeValidAtDependsOnStage; }

protected:
    explicit Implementation(const T& defaultValue, int numCacheEntries=1)
    :   presumeValidAtDependsOnStage(false),
        defaultValue(defaultValue),
        derivIx(numCacheEntries) 
    {
        Measure_Num<T>::makeZeroLike(defaultValue, zeroValue);
    }

    /** Argument \a numCacheEntries should be one greater than the number of 
    derivatives; that is, there is room for the value ("0th" derivative) also. 
    The default is to allocate just room for the value. **/
    explicit Implementation(int numCacheEntries=1) 
    :   presumeValidAtDependsOnStage(false),
        defaultValue(),
        derivIx(numCacheEntries) 
    {
        Measure_Num<T>::makeZeroLike(defaultValue, zeroValue);
    }

    /** Copy constructor copies the \e number of cache entries from the source,
    but not the cache indices themselves as those must be allocated uniquely 
    for the copy. **/
    Implementation(const Implementation& source) 
    :   presumeValidAtDependsOnStage(source.presumeValidAtDependsOnStage),
        defaultValue(source.defaultValue),
        derivIx(source.derivIx.size()) 
    {
        Measure_Num<T>::makeZeroLike(defaultValue, zeroValue);
    }


    /** Return the number of elements in the data type of this %Measure; for
    Vector measures this is determined by the size of the default value. **/
    int size() const {return Measure_Num<T>::size(defaultValue);}

    /** Return the number of cache entries allocated for the value and
    derivatives of this %Measure. **/
    int getNumCacheEntries() const {return (int)derivIx.size();}

    /** Get a const reference to the value stored in one of this %Measure's
    cache entries, indexed by the derivative order (with the value treated as
    the 0th derivative). **/
    const T& getCacheEntry(const State& s, int derivOrder) const {
        SimTK_ERRCHK2(0 <= derivOrder && derivOrder < getNumCacheEntries(),
            "Measure_<T>::Implementation::getCacheEntry()",
            "Derivative order %d is out of range; only %d cache entries"
            " were allocated.", derivOrder, getNumCacheEntries());

        return Value<T>::downcast(
            this->getSubsystem().getCacheEntry(s, derivIx[derivOrder]));
    }

    /** Get a writable reference to the value stored in one of this %Measure's
    cache entries, indexed by the derivative order (with the value treated as
    the 0th derivative). **/
    T& updCacheEntry(const State& s, int derivOrder) const {
        SimTK_ERRCHK2(0 <= derivOrder && derivOrder < getNumCacheEntries(),
            "Measure_<T>::Implementation::updCacheEntry()",
            "Derivative order %d is out of range; only %d cache entries"
            " were allocated.", derivOrder, getNumCacheEntries());

        return Value<T>::updDowncast(
            this->getSubsystem().updCacheEntry(s, derivIx[derivOrder]));
    }

    /** Determine whether a particular one of this %Measure's cache entries has
    already been realized since the given state was modified. **/
    bool isCacheValueRealized(const State& s, int derivOrder) const {
        SimTK_ERRCHK2(0 <= derivOrder && derivOrder < getNumCacheEntries(),
            "Measure_<T>::Implementation::isCacheValueRealized()",
            "Derivative order %d is out of range; only %d cache entries"
            " were allocated.", derivOrder, getNumCacheEntries());

        return this->getSubsystem().isCacheValueRealized(s, derivIx[derivOrder]);
    }

    /** Mark one of this %Measure's cache entries up to date; call this after
    you have calculated a value or derivative and stored it in the
    corresponding cache entry. **/
    void markCacheValueRealized(const State& s, int derivOrder) const {
        SimTK_ERRCHK2(0 <= derivOrder && derivOrder < getNumCacheEntries(),
            "Measure_<T>::Implementation::markCacheValueRealized()",
            "Derivative order %d is out of range; only %d cache entries"
            " were allocated.", derivOrder, getNumCacheEntries());

        this->getSubsystem().markCacheValueRealized(s, derivIx[derivOrder]);
    }

    /** Invalidate one of this %Measure's cache entries. This is not normally
    necessary since the cache entries will be invalidated automatically when
    state variables they depend on change. However, this can be useful in
    some cases, particularly during debugging and testing. **/
    void markCacheValueNotRealized(const State& s, int derivOrder) const {
        SimTK_ERRCHK2(0 <= derivOrder && derivOrder < getNumCacheEntries(),
            "Measure_<T>::Implementation::markCacheValueNotRealized()",
            "Derivative order %d is out of range; only %d cache entries"
            " were allocated.", derivOrder, getNumCacheEntries());

        this->getSubsystem().markCacheValueNotRealized(s, derivIx[derivOrder]);
    }

    // VIRTUALS //
    // Ordinals must retain the same meaning from release to release
    // to preserve binary compatibility.

    /** Concrete measures can override this to allocate Topology-stage
    resources. **/
    /* 0*/virtual void realizeMeasureTopologyVirtual(State&) const {}

    /** Concrete measures must override this if the state cache is used for
    precalculated values or derivatives. **/
    /* 1*/virtual void 
    calcCachedValueVirtual(const State&, int derivOrder, T& value) const
    {   SimTK_ERRCHK1_ALWAYS(!"implemented", 
        "Measure_<T>::Implementation::calcCachedValueVirtual()",
        "This method should have been overridden by the derived"
        " Measure but was not. It is needed to calculate the"
        " cached value for derivOrder=%d.", derivOrder); }

    /** This is only called when derivOrder >= the number of cache 
    entries we have, but still <= the number of derivatives the
    %Measure says it can deliver. You don't need to override this if that
    condition can't occur. This is commonly used for functions whose derivatives
    above a certain order are zero. **/
    /* 2*/virtual const T& 
    getUncachedValueVirtual(const State&, int derivOrder) const
    {   SimTK_ERRCHK1_ALWAYS(!"implemented", 
            "Measure_<T>::Implementation::getUncachedValueVirtual()",
            "This method should have been overridden by the derived"
            " Measure but was not. It is needed to return the uncached"
            " value at derivOrder=%d.", derivOrder);
        return *reinterpret_cast<T*>(0);
    }

    /** Return a reference to a zero of the same type and size as this 
    %Measure's value. **/
    const T& getValueZero() const {return zeroValue;}

private:
    // Satisfy the realizeTopology() pure virtual here now that we know the 
    // data type T. Allocate lazy- or auto-validated- cache entries depending 
    // on the setting of presumeValidAtDependsOnStage.
    void realizeTopology(State& s) const FINAL_11 {
        // Allocate cache entries. Initialize the value cache entry to
        // the given defaultValue; all the derivative cache entries should be
        // initialized to a NaN of the same size.
        if (getNumCacheEntries()) {
            derivIx[0] = presumeValidAtDependsOnStage 
                ? this->getSubsystem().allocateCacheEntry
                        (s, getDependsOnStage(0), new Value<T>(defaultValue))
                : this->getSubsystem().allocateLazyCacheEntry
                        (s, getDependsOnStage(0), new Value<T>(defaultValue));

            if (getNumCacheEntries() > 1) {
                T nanValue; Measure_Num<T>::makeNaNLike(defaultValue, nanValue);
                for (int i=1; i < getNumCacheEntries(); ++i) {
                    derivIx[i] = presumeValidAtDependsOnStage 
                        ? this->getSubsystem().allocateCacheEntry
                            (s, getDependsOnStage(i), new Value<T>(nanValue))
                        : this->getSubsystem().allocateLazyCacheEntry
                            (s, getDependsOnStage(i), new Value<T>(nanValue));
                }
            }
        }

        // Call the concrete class virtual if any.
        realizeMeasureTopologyVirtual(s);
    }

//------------------------------------------------------------------------------
private:
    // TOPOLOGY STATE
    bool    presumeValidAtDependsOnStage;
    T       defaultValue;
    T       zeroValue;

    // TOPOLOGY CACHE
    mutable Array_<CacheEntryIndex> derivIx;
};



//==============================================================================
//                       CONSTANT :: IMPLEMENTATION
//==============================================================================
template <class T>
class Measure_<T>::Constant::Implementation 
:   public Measure_<T>::Implementation 
{
public:
    // We don't want the base class to allocate *any* cache entries.
    Implementation() : Measure_<T>::Implementation(0) {}
    explicit Implementation(const T& value) 
    :   Measure_<T>::Implementation(value,0) {}

    /** Changing the value of a %Constant measure is a topological change;
    if this is a Vector measure you can change the size here too. **/
    void setValue(const T& v) {this->setDefaultValue(v);}

    // Implementations of virtual methods.
    // Measure_<T> virtuals:
    // No cached values.

    const T& getUncachedValueVirtual(const State&, int derivOrder) const 
        OVERRIDE_11
    {   return derivOrder>0 ? this->getValueZero() : this->getDefaultValue(); }

    // AbstractMeasure virtuals:
    Implementation* cloneVirtual() const OVERRIDE_11
    {   return new Implementation(*this); }
    Stage getDependsOnStageVirtual(int derivOrder) const OVERRIDE_11 
    {   return derivOrder>0 ? Stage::Empty : Stage::Topology; }
    int getNumTimeDerivativesVirtual() const OVERRIDE_11 
    {   return std::numeric_limits<int>::max(); }
};



//==============================================================================
//                        MEASURE ZERO and ONE
//==============================================================================
// These had to wait for Constant::Implementation to be declared.

template <class T> inline
Measure_<T>::Zero::Zero() : Constant(T(0)) {}
template <class T> inline
Measure_<T>::Zero::Zero(Subsystem& sub) : Constant(sub, T(0)) {}

inline Measure_< Vector >::Zero::Zero(int size) 
:   Constant(Vector(size, Real(0))) {}
inline Measure_< Vector >::Zero::Zero(Subsystem& sub, int size) 
:  Constant(sub, Vector(size, Real(0))) {}

template <class T> inline
Measure_<T>::One::One() : Constant(T(1)) {}
template <class T> inline
Measure_<T>::One::One(Subsystem& sub) : Constant(sub, T(1)) {}

inline Measure_< Vector >::One::One(int size) 
:   Constant(Vector(size, Real(1))) {}
inline Measure_< Vector >::One::One(Subsystem& sub, int size) 
:  Constant(sub, Vector(size, Real(1))) {}



//==============================================================================
//                         TIME :: IMPLEMENTATION
//==============================================================================
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
        OVERRIDE_11
    {   return derivOrder==0 ? s.getTime()
            : (derivOrder==1 ? SimTK::One 
                             : SimTK::Zero); } 

    // AbstractMeasure virtuals:
    Implementation* cloneVirtual() const OVERRIDE_11
    {   return new Implementation(*this); }
    Stage getDependsOnStageVirtual(int derivOrder) const OVERRIDE_11
    {   return derivOrder>0 ? Stage::Empty : Stage::Time; }

    // Value is t, 1st derivative is 1, the rest are 0.
    int getNumTimeDerivativesVirtual() const OVERRIDE_11 
    {   return std::numeric_limits<int>::max(); }
};



//==============================================================================
//                        VARIABLE :: IMPLEMENTATION
//==============================================================================
template <class T>
class Measure_<T>::Variable::Implementation 
:   public Measure_<T>::Implementation 
{
public:
    // We don't want the base class to allocate *any* cache entries;
    // we'll use the variable as its own value and zeroes for all
    // the derivatives.
    Implementation() 
    :   Measure_<T>::Implementation(0),
        invalidatedStage(Stage::Empty) {}

    Implementation(Stage invalidated, const T& defaultValue) 
    :   Measure_<T>::Implementation(defaultValue, 0),
        invalidatedStage(invalidated) {}

    // Copy constructor should not copy the variable.
    Implementation(const Implementation& source)
    :   Measure_<T>::Implementation(source.getDefaultValue(), 0),
        invalidatedStage(source.invalidatedStage) {}

    void setInvalidatedStage(Stage invalidates) {
        invalidatedStage = invalidates;
        this->invalidateTopologyCache();
    }

    Stage getInvalidatedStage() const {return invalidatedStage;}

    /** Change the value of this %Measure in the given \a state. Invalidates
    cache entries in that \a state for any stage at or above the "invalidates" 
    stage that was set when this %Measure was constructed. **/
    void setValue(State& state, const T& value) const 
    {   updVarValue(state) = value; }

    // Implementations of virtual methods.
    Implementation* cloneVirtual() const OVERRIDE_11
    {   return new Implementation(*this); }

    int getNumTimeDerivativesVirtual() const OVERRIDE_11 
    {   return std::numeric_limits<int>::max(); }

    // Discrete variable is available after Model stage; but all its 
    // derivatives are zero so are always available.
    Stage getDependsOnStageVirtual(int derivOrder) const OVERRIDE_11 
    {   return derivOrder>0 ? Stage::Empty : Stage::Model;}

    const T& getUncachedValueVirtual(const State& s, int derivOrder) const
        OVERRIDE_11
    {   return derivOrder>0 ? this->getValueZero() : getVarValue(s); }

    // No cached values.

    void realizeMeasureTopologyVirtual(State& s) const OVERRIDE_11 {
        discreteVarIndex = this->getSubsystem().allocateDiscreteVariable
            (s, invalidatedStage, new Value<T>(getDefaultValue()));
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

    // TOPOLOGY CACHE
    mutable DiscreteVariableIndex discreteVarIndex;
};



//==============================================================================
//                         RESULT :: IMPLEMENTATION
//==============================================================================
template <class T>
class Measure_<T>::Result::Implementation 
:   public Measure_<T>::Implementation 
{
public:
    // We want the base class to allocate a single cache entry of type T.
    Implementation() 
    :   Measure_<T>::Implementation(1), 
        dependsOnStage(Stage::Topology), invalidatedStage(Stage::Infinity) {}

    Implementation(Stage dependsOn, Stage invalidated) 
    :   Measure_<T>::Implementation(1), 
        dependsOnStage(dependsOn==Stage::Empty ? Stage::Topology : dependsOn), 
        invalidatedStage(invalidated)
    {   SimTK_ERRCHK2_ALWAYS(invalidated > dependsOn,"Measure::Result::ctor()",
            "Got invalidated stage %s and dependsOn stage %s which is illegal "
            "because the invalidated stage must be later than dependsOn.",
            invalidated.getName().c_str(), dependsOn.getName().c_str());
    }

    // Copy constructor will not copy the cache entry index.
    Implementation(const Implementation& source)
    :   Measure_<T>::Implementation(source),
        dependsOnStage(source.dependsOnStage), 
        invalidatedStage(source.invalidatedStage) {} 

    void setDependsOnStage(Stage dependsOn) {
        if (dependsOn == Stage::Empty) dependsOn = Stage::Topology;
        SimTK_ERRCHK2_ALWAYS(dependsOn < getInvalidatedStage(),
            "Measure::Result::setDependsOnStage()",
            "The provided dependsOn stage %s is illegal because it is not "
            "less than the current invalidated stage %s. Change the "
            "invalidated stage first with setInvalidatedStage().",
            dependsOn.getName().c_str(), 
            getInvalidatedStage().getName().c_str());

        dependsOnStage = dependsOn;
        this->invalidateTopologyCache();
    }

    void setInvalidatedStage(Stage invalidated) {
        SimTK_ERRCHK2_ALWAYS(invalidated > getDependsOnStage(),
            "Measure::Result::setInvalidatedStage()",
            "The provided invalidated stage %s is illegal because it is not "
            "greater than the current dependsOn stage %s. Change the "
            "dependsOn stage first with setDependsOnStage().",
            invalidated.getName().c_str(),
            getDependsOnStage().getName().c_str());

        invalidatedStage = invalidated;
        this->invalidateTopologyCache();
    }


    Stage getDependsOnStage()   const {return dependsOnStage;}
    Stage getInvalidatedStage() const {return invalidatedStage;}


    void markAsValid(const State& state) const
    {   const Stage subsystemStage = this->getSubsystem().getStage(state);
        SimTK_ERRCHK3_ALWAYS(subsystemStage >= getDependsOnStage().prev(),
            "Measure::Result::markAsValid()",
            "This Result Measure cannot be marked valid in a State where this "
            "measure's Subsystem has been realized only to stage %s, because "
            "its value was declared to depend on stage %s. To mark it valid, "
            "we require that the State have been realized at least to the "
            "previous stage (%s in this case); that is, you must at least be "
            "*working on* the dependsOn stage in order to claim this result is "
            "available.",
            subsystemStage.getName().c_str(),
            getDependsOnStage().getName().c_str(),
            getDependsOnStage().prev().getName().c_str());
        this->markCacheValueRealized(state, 0); }

    bool isValid(const State& state) const
    {   return this->isCacheValueRealized(state, 0); }
    
    void markAsNotValid(const State& state) const
    {   this->markCacheValueNotRealized(state, 0); 
        state.invalidateAllCacheAtOrAbove(invalidatedStage); }

    T& updValue(const State& state) const 
    {   markAsNotValid(state); return this->updCacheEntry(state, 0); }


    // Implementations of virtual methods.
    Implementation* cloneVirtual() const OVERRIDE_11
    {   return new Implementation(*this); }

    int getNumTimeDerivativesVirtual() const OVERRIDE_11 {return 0;} 

    /** Cache value is available after its "depends on" stage has been 
    realized; but all its derivatives are zero so are always available. **/
    Stage getDependsOnStageVirtual(int derivOrder) const OVERRIDE_11 
    {   return derivOrder>0 ? Stage::Empty : dependsOnStage;}

    void calcCachedValueVirtual(const State&, int derivOrder, T& value) const
        OVERRIDE_11
    {   SimTK_ERRCHK_ALWAYS(!"calcCachedValueVirtual() implemented",
        "Measure_<T>::Result::getValue()",
        "Measure_<T>::Result::getValue() was called when the value was not "
        "yet valid. For most Measure types, this would have initiated "
        "computation of the value, but Result measures must have their values "
        "calculated and set externally, and then marked valid."); }

private:
    // TOPOLOGY STATE
    Stage   dependsOnStage;
    Stage   invalidatedStage;
};



//==============================================================================
//                        SINUSOID :: IMPLEMENTATION
//==============================================================================
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
    Implementation* cloneVirtual() const OVERRIDE_11
    {   return new Implementation(*this); }

    int getNumTimeDerivativesVirtual() const OVERRIDE_11 {return NumDerivs;}

    Stage getDependsOnStageVirtual(int order) const OVERRIDE_11 
    {   return Stage::Time; }

    void calcCachedValueVirtual(const State& s, int derivOrder, T& value) const
        OVERRIDE_11
    {
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



//==============================================================================
//                          PLUS :: IMPLEMENTATION
//==============================================================================
template <class T>
class Measure_<T>::Plus::Implementation : public Measure_<T>::Implementation {
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
    Implementation* cloneVirtual() const OVERRIDE_11 
    {   return new Implementation(*this); }

    // TODO: Let this be settable up to the min number of derivatives 
    // provided by the arguments.
    int getNumTimeDerivativesVirtual() const OVERRIDE_11 {return 0;} 
    //{   return std::min(left.getNumTimeDerivatives(), 
    //                    right.getNumTimeDerivatives()); }

    Stage getDependsOnStageVirtual(int order) const OVERRIDE_11 
    {   return Stage(std::max(left.getDependsOnStage(order),
                              right.getDependsOnStage(order))); }


    void calcCachedValueVirtual(const State& s, int derivOrder, T& value) const
        OVERRIDE_11
    {
        value = left.getValue(s,derivOrder) + right.getValue(s,derivOrder);
    }

    // There are no uncached values.

private:
    // TOPOLOGY STATE
    Measure_<T> left;
    Measure_<T> right;

    // TOPOLOGY CACHE
    // nothing
};



//==============================================================================
//                          MINUS :: IMPLEMENTATION
//==============================================================================
template <class T>
class Measure_<T>::Minus::Implementation : public Measure_<T>::Implementation {
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
    Implementation* cloneVirtual() const OVERRIDE_11 
    {   return new Implementation(*this); }

    // TODO: Let this be settable up to the min number of derivatives 
    // provided by the arguments.
    int getNumTimeDerivativesVirtual() const OVERRIDE_11 {return 0;} 
    //{   return std::min(left.getNumTimeDerivatives(), 
    //                    right.getNumTimeDerivatives()); }

    Stage getDependsOnStageVirtual(int order) const OVERRIDE_11 
    {   return Stage(std::max(left.getDependsOnStage(order),
                              right.getDependsOnStage(order))); }


    void calcCachedValueVirtual(const State& s, int derivOrder, T& value) const
        OVERRIDE_11
    {
        value = left.getValue(s,derivOrder) - right.getValue(s,derivOrder);
    }

    // There are no uncached values.

private:
    // TOPOLOGY STATE
    Measure_<T> left;
    Measure_<T> right;

    // TOPOLOGY CACHE
    // nothing
};



//==============================================================================
//                          SCALE :: IMPLEMENTATION
//==============================================================================
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
    Implementation* cloneVirtual() const OVERRIDE_11 
    {   return new Implementation(*this); }

    // TODO: Let this be settable up to the min number of derivatives 
    // provided by the arguments.
    int getNumTimeDerivativesVirtual() const OVERRIDE_11 {return 0;} 
    //{   return std::min(left.getNumTimeDerivatives(), 
    //                    right.getNumTimeDerivatives()); }

    Stage getDependsOnStageVirtual(int order) const OVERRIDE_11
    {   return operand.getDependsOnStage(order); }


    void calcCachedValueVirtual(const State& s, int derivOrder, T& value) const
        OVERRIDE_11
    {
        value = factor * operand.getValue(s,derivOrder);
    }

    // There are no uncached values.

private:
    // TOPOLOGY STATE
    Real        factor;
    Measure_<T> operand;

    // TOPOLOGY CACHE
    // nothing
};



//==============================================================================
//                        INTEGRATE :: IMPLEMENTATION
//==============================================================================
/** The implementation for %Integrate measures allocates a continuous state
variable or variables from the State's z pool and generates zdot values to be 
integrated into those z variables. The z's are then copied into a type T,
Time-stage cache entry so that we can return the value as a type T reference.
Derivative requests are passed through to the integrand so only one cache
entry is required here. **/
template <class T>
class Measure_<T>::Integrate::Implementation 
:   public Measure_<T>::Implementation {
public:
    /** The derivative and initialConditions Measures will be empty handles if
    this is default constructed. **/
    Implementation() : Measure_<T>::Implementation(1) {}

    /** Here we're shallow-copying the Measure handles so we'll be referring to
    the original Measures. **/
    Implementation(const Measure_<T>& deriv, const Measure_<T>& ic,
                   const T& defaultValue)
    :   Measure_<T>::Implementation(defaultValue, 1), 
        derivMeasure(deriv), icMeasure(ic) {}

    /** Copy constructor shallow-copies the referenced measures, but we don't 
    want to share our state variables. **/
    Implementation(const Implementation& source)
    :   Measure_<T>::Implementation(source.getDefaultValue(), 1), 
        derivMeasure(source.derivMeasure), icMeasure(source.icMeasure) {}

    /** Set the value of the state variables(s) that hold the integral. This
    cannot be used to change the size if the type T is a Vector; the supplied
    \a value must be the same length as the default value of this %Measure. **/
    void setValue(State& s, const T& value) const
    {   assert(zIndex >= 0);
        for (int i=0; i < this->size(); ++i)
            this->getSubsystem().updZ(s)[zIndex+i] = 
                Measure_Num<T>::get(value, i); }
    
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
    Implementation* cloneVirtual() const OVERRIDE_11 
    {   return new Implementation(*this); }

    /** This measure has one more time derivative than the integrand. **/
    int getNumTimeDerivativesVirtual() const OVERRIDE_11 
    {   int integralDerivs = getDerivativeMeasure().getNumTimeDerivatives();
        // Careful - can't add 1 to max int and stay an int.
        if (integralDerivs < std::numeric_limits<int>::max())
            ++integralDerivs;        
        return integralDerivs; }

    void calcCachedValueVirtual(const State& s, int derivOrder, T& value) const
        OVERRIDE_11
    {   assert(derivOrder == 0); // only one cache entry
        assert(Measure_Num<T>::size(value) == this->size());
        assert(zIndex.isValid());
        const Vector& allZ = this->getSubsystem().getZ(s);
        for (int i=0; i < this->size(); ++i)
            Measure_Num<T>::upd(value,i) = allZ[i];
    }

    const T& getUncachedValueVirtual(const State& s, int derivOrder) const
        OVERRIDE_11
    {   assert(derivOrder > 0); // 0th entry is cached
        return getDerivativeMeasure().getValue(s, derivOrder-1);
    }

    Stage getDependsOnStageVirtual(int derivOrder) const OVERRIDE_11 
    {   return derivOrder>0 
            ? getDerivativeMeasure().getDependsOnStage(derivOrder-1)
            : Stage::Time; }

    /** Initialize the state to the current value of the initial condition
    measure, if there is one, otherwise to the default value. **/
    void initializeVirtual(State& s) const OVERRIDE_11 {
        assert(zIndex.isValid());
        Vector& allZ = this->getSubsystem().updZ(s);
        if (!icMeasure.isEmptyHandle()) {
             const T& ic = icMeasure.getValue(s);
             for (int i=0; i < this->size(); ++i)
                 allZ[zIndex+i] = Measure_Num<T>::get(ic,i);
        } else {
             for (int i=0; i < this->size(); ++i)
                 allZ[zIndex+i] = Measure_Num<T>::get(getDefaultValue(),i);
        }
    }

    /** Allocate one Real continuous state variable z per element of this 
    %Measure's data type T, using the default value to determine 
    how many are needed (if that's not part of the type T), and initialize them
    to the corresponding element from the default value. **/
    void realizeMeasureTopologyVirtual(State& s) const OVERRIDE_11 {
        Vector init(this->size());
        for (int i=0; i < this->size(); ++i) 
            init[i] = Measure_Num<T>::get(getDefaultValue(),i);
        zIndex = this->getSubsystem().allocateZ(s, init);
    }

    /** Set the zdots to the integrand (derivative measure) value. If no
    integrand was provided it is treated as though it were zero. **/
    void realizeMeasureAccelerationVirtual(const State& s) const OVERRIDE_11 {
        assert(zIndex.isValid());
        Vector& allZDot = this->getSubsystem().updZDot(s);
        if (!derivMeasure.isEmptyHandle()) {
            const T& deriv = derivMeasure.getValue(s);
             for (int i=0; i < this->size(); ++i)
                 allZDot[zIndex+i] = Measure_Num<T>::get(deriv,i);
        } else {
            allZDot(zIndex,this->size()) = 0; // derivative is zero
        }
    }

private:
    // TOPOLOGY STATE
    Measure_<T> derivMeasure; // just handles
    Measure_<T> icMeasure;

    // TOPOLOGY CACHE
    mutable ZIndex zIndex;  // This is the first index if more than one z.
};



//==============================================================================
//                      DIFFERENTIATE :: IMPLEMENTATION
//==============================================================================

/** @cond **/ // Hide from Doxygen.
// This helper class is the contents of the discrete state variable and 
// corresponding cache entry maintained by this measure. The variable is 
// auto-update, meaning the value of the cache entry replaces the state 
// variable at the start of each step.
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
/** @endcond **/


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
        const bool good0      = prevResult.derivIsGood;

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



//==============================================================================
//                         EXTREME :: IMPLEMENTATION
//==============================================================================
template <class T>
class Measure_<T>::Extreme::Implementation : public Measure_<T>::Implementation 
{
    typedef typename Measure_<T>::Extreme Extreme;
    typedef typename Extreme::Operation   Operation;
public:
    /** Default constructor leaves the operand measure unspecified; no base
    class cache entries are allocated. **/
    Implementation() 
    :   Measure_<T>::Implementation(0), operation(Extreme::MaxAbs) {}

    /** Construct a measure that returns the extreme value taken on by the
    operand measure during a time stepping study. **/
    Implementation(const Measure_<T>& operand, Operation op)
    :   Measure_<T>::Implementation(0), operand(operand), operation(op) {}

    // Default copy constructor gives us a new Implementation object,
    // but with reference to the *same* operand measure.

    /** Set the operand measure for this %Extreme measure; this is a Topology
    stage change so you'll have to call realizeTopology() again if you call
    this. **/
    void setOperandMeasure(const Measure_<T>& operand) {
        this->operand = operand;
        this->invalidateTopologyCache();
    }

    /** Set the particular operation to be performed by this %Extreme measure;
    this is a Topology stage change so you'll have to call realizeTopology() 
    again if you call this. **/
    void setOperation(Operation op) {
        this->operation = op;
        this->invalidateTopologyCache();
    }

    /** Return a reference to the operand measure for this %Extreme measure. **/
    const Measure_<T>& getOperandMeasure() const {return operand;}

    /** Return the particular operation being performed by this %Extreme
    measure. **/
    Operation getOperation() const {return operation;}

    /** Set the current extreme value stored in this %Extreme measure's state
    variable. **/
    void setValue(State& s, const T& value) const {
        assert(extremeIx.isValid());
        const Subsystem& subsys = this->getSubsystem();
        T& prevMin = Value<T>::updDowncast
                            (subsys.updDiscreteVariable(s,extremeIx));
        prevMin = value;
    }

    /** Return the time at which the extreme was last updated. This will be
    the current time if the operand is currently at its lowest value, otherwise
    it will be sometime in the past. **/
    Real getTimeOfExtremeValue(const State& s) const {
        const Subsystem& subsys = this->getSubsystem();
        const bool hasNewExtreme = ensureExtremeHasBeenUpdated(s);
        Real tUpdate;
        if (hasNewExtreme)
            tUpdate = s.getTime(); // i.e., now
        else
            tUpdate = subsys.getDiscreteVarLastUpdateTime(s,extremeIx);
        return tUpdate;
    }

    // Implementations of virtual methods.

    // This uses the default copy constructor.
    Implementation* cloneVirtual() const OVERRIDE_11 
    {   return new Implementation(*this); }

    /** Extreme(f(t)) has the same number of derivatives as f except that
    they are all zero unless f(t) is a new extreme. **/
    int getNumTimeDerivativesVirtual() const OVERRIDE_11
    {   return operand.getNumTimeDerivatives(); }

    /** The depends-on stage for this measure is the same as for its 
    operand. **/
    Stage getDependsOnStageVirtual(int order) const OVERRIDE_11
    {   return operand.getDependsOnStage(order); }


    /** We're not using the Measure_<T> base class cache services, but
    we do have one of our own. It looks uncached from the base class
    point of view which is why we're implementing it here. **/
    const T& getUncachedValueVirtual(const State& s, int derivOrder) const
        OVERRIDE_11
    {   
        const Subsystem& subsys = this->getSubsystem();
        const bool hasNewExtreme = ensureExtremeHasBeenUpdated(s);
        if (derivOrder > 0) {
            // TODO: should be handled elementwise and zero unless the
            // derivative is acting in the direction that changes the 
            // extreme.
            return hasNewExtreme ? operand.getValue(s, derivOrder)
                                 : this->getValueZero();
        }
        if (hasNewExtreme) {
            const T& newExt = Value<T>::downcast
                                (subsys.getDiscreteVarUpdateValue(s,extremeIx));
            return newExt;
        } else {
            const T& currentExt = Value<T>::downcast
                                (subsys.getDiscreteVariable(s,extremeIx));
            return currentExt;
        }
    }

    /** At start of a time stepping study, this should be called to set the
    current extreme value to the current value of the operand. **/
    void initializeVirtual(State& s) const OVERRIDE_11 {
        setValue(s, operand.getValue(s));
    }

    /** Allocate the auto-updated state variable that holds the extreme seen
    so far. We'll assume that changes to this variable invalidate Dynamics
    (force) stage so that any forces that depend on it will be recomputed if
    it changes. **/
    void realizeMeasureTopologyVirtual(State& s) const OVERRIDE_11 {
        // TODO: this should be NaN once initialization is working properly.
        T initVal = getDefaultValue();
        switch(operation) {
        case Minimum: initVal = Infinity; break;
        case Maximum: initVal = -Infinity; break;
        case MinAbs:  initVal = Infinity; break;
        case MaxAbs:  initVal = 0; break;
        };
        extremeIx = this->getSubsystem()
            .allocateAutoUpdateDiscreteVariable(s, Stage::Dynamics,
                new Value<T>(initVal), operand.getDependsOnStage(0));
    }

    /** Here we make sure that the cache entry is updated if the current value
    of the operand is less than the previous one, and return a bool indicating
    whether we have a new minimum. We don't want to create an update entry
    unless the minimum has changed, because we would like the state 
    variable's last update value to reflect the last actual change. **/
    bool ensureExtremeHasBeenUpdated(const State& s) const {
        assert(extremeIx.isValid());
        const Subsystem& subsys = this->getSubsystem();
        if (subsys.isDiscreteVarUpdateValueRealized(s,extremeIx))
            return true; // already updated to new minimum

        const T& prevExtreme = Value<T>::downcast
                                    (subsys.getDiscreteVariable(s,extremeIx));
        const T& f = operand.getValue(s);
        // Search to see if any element has reached a new extreme.
        bool foundNewExt = false;
        for (int i=0; i < this->size() && !foundNewExt; ++i) 
            foundNewExt = isNewExtreme(Measure_Num<T>::get(f,i), 
                                       Measure_Num<T>::get(prevExtreme,i));
        if (!foundNewExt)
            return false; // no change to output

        // We have encountered a new record holder for minimum.
        T& newExtreme = Value<T>::updDowncast
                                (subsys.updDiscreteVarUpdateValue(s,extremeIx));

        for (int i=0; i < this->size(); ++i)
            Measure_Num<T>::upd(newExtreme,i) =
                extremeOf(Measure_Num<T>::get(f,i),
                          Measure_Num<T>::get(prevExtreme,i));

        subsys.markDiscreteVarUpdateValueRealized(s,extremeIx);
        return true;
    }
private:
    // Return true if newVal is "more extreme" than oldExtreme, according
    // to the operation we're performing.
    bool isNewExtreme(const typename Measure_Num<T>::Element& newVal,
                      const typename Measure_Num<T>::Element& oldExtreme) const
    {
        switch (operation) {
        case Extreme::Maximum: return newVal > oldExtreme;
        case Extreme::Minimum: return newVal < oldExtreme;
        case Extreme::MaxAbs: return std::abs(newVal) > std::abs(oldExtreme);
        case Extreme::MinAbs: return std::abs(newVal) < std::abs(oldExtreme);
        };
        SimTK_ASSERT1_ALWAYS(!"recognized", 
            "Measure::Extreme::Implementation::isNewExtreme(): "
            "unrecognized operation %d", (int)operation);
        return false; /*NOTREACHED*/
    }

    // Given the value of one element of the operand, and that value's time
    // derivative, determine whether the derivative is moving the element
    // into the "more extreme" direction, according to the operation.
    bool isExtremeDir(const typename Measure_Num<T>::Element& value,
                      const typename Measure_Num<T>::Element& deriv) const 
    {
        const int sv = sign(value), sd = sign(deriv);
        if (sd == 0) return false; // derivative is zero; not changing
        switch (operation) {
        case Extreme::Maximum: return sd ==  1; // getting larger
        case Extreme::Minimum: return sd == -1; // getting smaller
        case Extreme::MaxAbs: return sv==0 || sd==sv; // abs is growing
        case Extreme::MinAbs: return sd == -sv;
        };
        SimTK_ASSERT1_ALWAYS(!"recognized", 
            "Measure::Extreme::Implementation::isExtremeDir(): "
            "unrecognized operation %d", (int)operation);
        return false; /*NOTREACHED*/
    }

    typename Measure_Num<T>::Element 
    extremeOf(const typename Measure_Num<T>::Element& newVal,
              const typename Measure_Num<T>::Element& oldExtreme) const
    {
        return isNewExtreme(newVal,oldExtreme) ? newVal : oldExtreme;
    }

    // TOPOLOGY STATE
    Measure_<T>                     operand;
    Operation                       operation;

    // TOPOLOGY CACHE
    mutable DiscreteVariableIndex   extremeIx;    // min so far; auto-update
};


} // namespace SimTK




#endif // SimTK_SimTKCOMMON_MEASURE_IMPLEMENTATION_H_
