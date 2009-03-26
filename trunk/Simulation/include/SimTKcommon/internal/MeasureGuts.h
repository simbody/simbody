#ifndef SimTK_SimTKCOMMON_MEASURE_GUTS_H_
#define SimTK_SimTKCOMMON_MEASURE_GUTS_H_

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

namespace SimTK {

    ///////////////////
    // MEASURE::GUTS //
    ///////////////////

/**
 * The abstract parent of all Measure_ Guts classes.
 */
class SimTK_SimTKCOMMON_EXPORT Measure::Guts {
protected:
    // This constructor is for use by concrete Subsystems. Note that this
    // serves as a default constructor since the argument has a default.
    explicit Guts(const std::string& name="<NONAME>")
    :   measureName(name), refCount(0), mySubsystem(0) {}

    Guts(const Guts& src)
    :   measureName(src.measureName), refCount(0), mySubsystem(0) {}
        
    Guts& operator=(const Guts& src)
    {   if (&src != this) {measureName=src.measureName; refCount=0; mySubsystem=0;}
        return *this; }

    // destructor is virtual

    // Increment the reference count and return its new value.
    int incrRefCount() const {return ++refCount;}

    // Decrement the reference count and return its new value.
    int decrRefCount() const {return --refCount;}

    const std::string& getName() const {return measureName;}
    Guts* clone() const {return cloneImpl();}

    void realizeTopology    (State& s)       const {realizeMeasureTopologyImpl(s);}
    void realizeModel       (State& s)       const {realizeMeasureModelImpl(s);}
    void realizeInstance    (const State& s) const {realizeMeasureInstanceImpl(s);}
    void realizeTime        (const State& s) const {realizeMeasureTimeImpl(s);}
    void realizePosition    (const State& s) const {realizeMeasurePositionImpl(s);}
    void realizeVelocity    (const State& s) const {realizeMeasureVelocityImpl(s);}
    void realizeDynamics    (const State& s) const {realizeMeasureDynamicsImpl(s);}
    void realizeAcceleration(const State& s) const {realizeMeasureAccelerationImpl(s);}
    void realizeReport      (const State& s) const {realizeMeasureReportImpl(s);}

    void  initialize(State& s)               const {initializeImpl(s);}
    Stage getValueDependence(const State& s) const {return getValueDependenceImpl(s);}

    void setSubsystem(Subsystem& sub, MeasureIndex mx) 
    {   assert(!mySubsystem && mx.isValid()); 
        mySubsystem = &sub; myIndex = mx; }

    bool isInSubsystem() const {return mySubsystem != 0;}
    const Subsystem& getSubsystem() const {assert(mySubsystem); return *mySubsystem;}
    Subsystem& updSubsystem() {assert(mySubsystem); return *mySubsystem;}
    MeasureIndex getSubsystemMeasureIndex() const {assert(mySubsystem); return myIndex;}

    void invalidateTopologyCache() const
    {   if (isInSubsystem()) getSubsystem().invalidateSubsystemTopologyCache(); }

    // VIRTUALS //
    // Ordinals must retain the same meaning from release to release
    // to preserve binary compatibility.

    /* 0*/virtual ~Guts() {}
    /* 1*/virtual Guts* cloneImpl() const = 0;

    /* 2*/virtual void realizeMeasureTopologyImpl(State&) const {}
    /* 3*/virtual void realizeMeasureModelImpl(State&) const {}
    /* 4*/virtual void realizeMeasureInstanceImpl(const State&) const {}
    /* 5*/virtual void realizeMeasureTimeImpl(const State&) const {}
    /* 6*/virtual void realizeMeasurePositionImpl(const State&) const {}
    /* 7*/virtual void realizeMeasureVelocityImpl(const State&) const {}
    /* 8*/virtual void realizeMeasureDynamicsImpl(const State&) const {}
    /* 9*/virtual void realizeMeasureAccelerationImpl(const State&) const {}
    /*10*/virtual void realizeMeasureReportImpl(const State&) const {}

    /*11*/virtual void initializeImpl(State&) const {}
    /*12*/virtual Stage getValueDependenceImpl(const State&) const = 0;

private:
    std::string     measureName;

    // These are set when this Measure is adopted by a Subsystem.
    Subsystem*      mySubsystem;
    MeasureIndex    myIndex;

    // Measures have shallow copy semantics so they share the Guts objects,
    // which are only deleted when the refCount goes to zero.
    mutable int     refCount;

friend class Measure;
friend class Subsystem::Guts;
friend class Subsystem::Guts::GutsRep;
};

    ////////////////////////////
    // MEASURE IMPLEMENTATION //
    ////////////////////////////

inline Measure::Measure(Guts* g) : guts(g)
{   if (guts) guts->incrRefCount();}
inline Measure::Measure(Subsystem& sub, Guts* g, const SetHandle&) : guts(g) {
    SimTK_ASSERT_ALWAYS(guts, "An empty Measure handle can't be put in a Subsystem.");
    guts->incrRefCount();
    sub.adoptMeasure(*this);
}
inline Measure::Measure(const Measure& src) : guts(0) {
    if (src.guts) {
        guts = src.guts->mySubsystem ? src.guts : src.guts->clone();
        guts->incrRefCount();
    }
}
inline Measure& Measure::operator=(const Measure& src) {
    if (&src != this) {
        if (guts && guts->decrRefCount()==0) delete guts;
        guts = src.guts->mySubsystem ? src.guts : src.guts->clone();
        guts->incrRefCount();
    }
    return *this;
}
inline Measure::~Measure()
{   if (guts && guts->decrRefCount()==0) delete guts;}
inline bool Measure::isInSubsystem() const
{   return guts && guts->isInSubsystem();}
inline const Subsystem& Measure::getSubsystem() const
{   assert(guts); return guts->getSubsystem();}
inline MeasureIndex Measure::getSubsystemMeasureIndex() const
{   assert(guts); return guts->getSubsystemMeasureIndex();}

inline Stage Measure::getValueDependence(const State& s) const
{   assert(guts); return guts->getValueDependence(s);}

    ///////////////////////
    // MEASURE_<T>::GUTS //
    ///////////////////////

/**
 * This is the base Guts class for all Measures whose value type is known. Logically
 * this is an abstract base class. However, you can create an empty Measure<T> handle
 * that can be assigned to any derived concrete Measure.
 */
template <class T>
class Measure_<T>::Guts : public Measure::Guts {
public:
    const T& getValue(const State& s) const {return getValueImpl(s);}

    virtual const T& getValueImpl(const State& s) const = 0;
};

    ///////////////////////
    // CONSTANT<T>::GUTS //
    ///////////////////////

template <class T>
class Measure::Constant_<T>::Guts : public Measure_<T>::Guts {
public:
    Guts(const T& value) : value(value) {}

    void setValue(const T& v) {
        value = v;
        this->invalidateTopologyCache();
    }

    // Implementations of virtual methods.
    Guts* cloneImpl() const {return new Guts(*this);}
    const T& getValueImpl(const State&) const {return value;}
    Stage   getValueDependenceImpl(const State&) const {return Stage::Topology;}

private:
    T value;
};


    ///////////////////////
    // SINUSOID<T>::GUTS //
    ///////////////////////

template <class T>
class Measure::Sinusoid_<T>::Guts : public Measure_<T>::Guts {
public:
    Guts(const T& amplitude, 
         const T& frequency, 
         const T& phase=T(0))
    :   a(amplitude), w(frequency), p(phase) {}

    // Implementations of virtual methods.
    Guts* cloneImpl() const {return new Guts(*this);}
    const T& getValueImpl(const State& s) const
    {   return getValueEntry(s);}

    Stage getValueDependenceImpl(const State&) const {return Stage::Time;}

    void realizeMeasureTopologyImpl(State& s) const {
        cacheEntryIndex = getSubsystem().allocateCacheEntry
            (s, Stage::Time, new Value<T>(CNT<T>::getNaN()));
    }
    void realizeMeasureTimeImpl(const State& s) const {
        assert(cacheEntryIndex >= 0);
        updValueEntry(s) = a*std::sin(w*s.getTime() + p);
    }

private:
    const T& getValueEntry(const State& s) const
    {   assert(cacheEntryIndex >= 0);
        return Value<T>::downcast(getSubsystem().getCacheEntry(s, cacheEntryIndex));}
    T& updValueEntry(const State& s) const
    {   assert(cacheEntryIndex >= 0);
        return Value<T>::downcast(getSubsystem().updCacheEntry(s, cacheEntryIndex));}

    // TOPOLOGY STATE
    T a, w, p;

    // TOPOLOGY CACHE
    mutable CacheEntryIndex cacheEntryIndex;
};

    ////////////////////////
    // INTEGRATE<T>::GUTS //
    ////////////////////////

template <class T>
class Measure::Integrate_<T>::Guts : public Measure_<T>::Guts {
public:
    Guts() : derivMeasure(0), icMeasure(0) {}
    Guts(const Measure_<T>& deriv, const Measure_<T>& ic)
    :   derivMeasure(&deriv), icMeasure(&ic) {}

    void setValue(State& s, const T& value) const
    {    assert(zIndex >= 0); getSubsystem().updZ(s)[zIndex] = value; }
    
    const Measure_<T>& getDerivativeMeasure() const
    {   assert(derivMeasure); return *derivMeasure; }
    const Measure_<T>& getInitialConditionMeasure() const
    {   assert(icMeasure); return *icMeasure; }

    void setDerivativeMeasure(const Measure_<T>& d)
    {   derivMeasure = &d; this->invalidateTopologyCache(); }
    void setInitialConditionMeasure(const Measure_<T>& ic)
    {   icMeasure = &ic; this->invalidateTopologyCache(); }

    // Implementations of virtuals.
    Guts* cloneImpl() const {return new Guts(*this);}
    const T& getValueImpl(const State& s) const
    {   assert(zIndex.isValid()); return getSubsystem().getZ(s)[zIndex]; }
    Stage getValueDependenceImpl(const State&) const {return Stage::Time;}

    void initializeImpl(State& s) const {
        assert(zIndex >= 0);
        Real& z = getSubsystem().updZ(s)[zIndex];
        if (icMeasure) z = icMeasure->getValue(s);
        else z = 0;
    }

    void realizeMeasureTopologyImpl(State& s) const {
        static const Vector zero(1, Real(0));
        zIndex = getSubsystem().allocateZ(s, zero);
    }

    void realizeMeasureAccelerationImpl(const State& s) const {
        assert(zIndex.isValid());
        Real& zdot = getSubsystem().updZDot(s)[zIndex];
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

#endif // SimTK_SimTKCOMMON_MEASURE_GUTS_H_
