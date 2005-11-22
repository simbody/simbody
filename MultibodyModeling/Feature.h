#ifndef SIMTK_FEATURE_H_
#define SIMTK_FEATURE_H_

/* Copyright (c) 2005-6 Stanford University and Michael Sherman.
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
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/** @file
 * User-visible definitions for the objects that go into building a multibody system.
 * This is not the data structure used at run time, so the emphasis is on 
 * nice behavior for the caller. We'll have plenty of time for speed later.
 *
 * Feature: Station, Direction, Frame, MassElement, ...
 * Placement: constant, expression or feature
 * Body: is a Frame, has (Feature,Placement) pairs
 *
 */


#include "SimbodyCommon.h"
#include "Placement.h"

#include "simtk/SimTK.h"
#include "simmatrix/SmallMatrix.h"

#include <iostream>
#include <cassert>

namespace simtk {

// Declared below. Indenting shows inheritance structure.
class Feature;
class   Station;
class   Direction;
class   Orientation;
class   Frame;
class   RealMeasure;
class     RealParameter;
class   StationMeasure;
class     StationParameter;
class   DirectionMeasure;
class   OrientationMeasure;

/**
 * Features form a tree because many Features have sub-Features (children).
 * Parent features own their children; destructing the parent destructs
 * all the children.
 *
 * Most features require placement in order to be useful (e.g.,
 * a Station has to have a location). A Feature's Placement may
 * have a constant value, or can involve sibling, parent, or 
 * ancestor Features but not its own children. 
 *
 * Every feature has a name and an index by which it is known
 * to its parent. Parentless features can exist but they can't
 * be placed. Parentless features still have a name but they do
 * not have an index.
 */
class Feature {
public:
    Feature() : rep(0) { }
    Feature(const Feature&);    // external placements are not copied or assigned
    Feature& operator=(const Feature&);
    ~Feature();


    // Read-only access to subfeatures.
    const Feature&          getSubfeature      (const String&) const; // generic
    const RealParameter&    getRealParameter   (const String&) const; // type checked
    const StationParameter& getStationParameter(const String&) const; //   "
    const RealMeasure&      getRealMeasure     (const String&) const;
    const StationMeasure&   getStationMeasure  (const String&) const;
    const Station&          getStation         (const String&) const;
    const Direction&        getDirection       (const String&) const;
    const Orientation&      getOrientation     (const String&) const;
    const Frame&            getFrame           (const String&) const;

    // Writable access to subfeatures, e.g. allowing placement.
    Feature&                updSubfeature      (const String&);   // generic
    RealParameter&          updRealParameter   (const String&);   // type checked
    StationParameter&       updStationParameter(const String&);   //   "
    RealMeasure&            updRealMeasure     (const String&);
    StationMeasure&         updStationMeasure  (const String&);
    Station&                updStation         (const String&);
    Direction&              updDirection       (const String&);
    Orientation&            updOrientation     (const String&);
    Frame&                  updFrame           (const String&);

    // Create a new subfeature on this feature with a given name and type, and
    // optionally create a placement for it using the prototype placement supplied.
    RealParameter&    addRealParameter
                        (const String&, const Placement& = Placement());
    RealMeasure&      addRealMeasure
                        (const String&, const Placement& = Placement());
    StationParameter& addStationParameter
                        (const String&, const Placement& = Placement());
    StationMeasure&   addStationMeasure
                        (const String&, const Placement& = Placement());
    Station&          addStation
                        (const String&, const Placement& = Placement());
    Direction&        addDirection
                        (const String&, const Placement& = Placement());
    Orientation&      addOrientation
                        (const String&, const Placement& = Placement());
    Frame&            addFrame
                        (const String&, const Placement& = Placement());

    // This is similar to the "add" routines above, except that the newly created
    // feature is modeled on the prototype feature supplied here. Again a placement
    // may be supplied or not.

    // This generic routine will create a subfeature of the specific type supplied
    // in the first argument. If a Placement is provided, it must be appropriate
    // for that kind of feature or a runtime error will occur.
    Feature&          addSubfeatureLike
                        (const Feature&,
                         const String&, const Placement& = Placement());
   
    bool           hasParentFeature() const;
    int            getIndexInParent() const; // -1 if no parent
    const Feature& getParentFeature() const;

    bool isSameFeature(const Feature& f) const {return &f == this;}

    // True if this is the same feature as f or if the feature's placement
    // depends on f's placement.
    bool dependsOn(const Feature& f) const;

    String getName()            const;
    String getFullName()        const; // "ancestors.parent.name"
    String getFeatureTypeName() const; // "Station", "Frame", etc. (for messages only)

    // Subfeatures of this feature
    int            getNSubfeatures()  const;
    const Feature& getSubfeature(int) const;
    Feature&       updSubfeature(int);

    const Feature& operator[](int i) const           {return getSubfeature(i);}
    Feature&       operator[](int i)                 {return updSubfeature(i);}
    const Feature& operator[](const String& s) const {return getSubfeature(s);}
    Feature&       operator[](const String& s)       {return updSubfeature(s);}

    // True if this feature has been placed; its children may still be unplaced.
    bool hasPlacement() const;
    const Placement& getPlacement() const;

    // Place this Feature using the supplied Placement expression
    // as a prototype. We will choose an owner Feature for the new
    // Placement, and then add a copy of the supplied one to that
    // owner. Then this Feature will refer to that copy as its Placement.
    void place(const Placement&);

    String toString(const String& linePrefix="") const;

    // For internal use only.
    bool                    hasRep() const {return rep != 0;}
    const class FeatureRep& getRep() const {assert(rep); return *rep;}
    class FeatureRep&       updRep()       {assert(rep); return *rep;}
    void                    setRep(FeatureRep* fp) {assert(!rep); rep=fp;}
protected:
    class FeatureRep* rep;
    friend class FeatureRep;
};

// Global operators involving Features.
std::ostream& operator<<(std::ostream& o, const Feature&);

// Although these operators appear to act on Features, they actually
// create a Placement referring to the Features and then perform
// the operations on the Placement.

// unary(feature)
inline Placement operator+(const Feature& f) {return  Placement(f);}
inline Placement operator-(const Feature& f) {return -Placement(f);}
inline RealPlacement      length(const Feature& f)    {return length(Placement(f));}
inline DirectionPlacement normalize(const Feature& f) {return normalize(Placement(f));}

// binary(feature,feature)
inline Placement operator+(const Feature& l, const Feature& r) 
  { return Placement(l)+Placement(r); }
inline Placement operator-(const Feature& l, const Feature& r) 
  { return Placement(l)-Placement(r); }
inline Placement operator*(const Feature& l, const Feature& r) 
  { return Placement(l)*Placement(r); }
inline Placement operator/(const Feature& l, const Feature& r) 
  { return Placement(l)/Placement(r); }

// binary(placement,feature)
inline Placement operator+(const Placement& l, const Feature& r) 
  { return l+Placement(r); }
inline Placement operator-(const Placement& l, const Feature& r) 
  { return l-Placement(r); }
inline Placement operator*(const Placement& l, const Feature& r) 
  { return l*Placement(r); }
inline Placement operator/(const Placement& l, const Feature& r) 
  { return l/Placement(r); }

// binary(feature,placement)
inline Placement operator+(const Feature& l, const Placement& r) 
  { return Placement(l)+r; }
inline Placement operator-(const Feature& l, const Placement& r) 
  { return Placement(l)-r; }
inline Placement operator*(const Feature& l, const Placement& r) 
  { return Placement(l)*r; }
inline Placement operator/(const Feature& l, const Placement& r) 
  { return Placement(l)/r; }

// binary(real,feature)
inline Placement operator+(const Real& l, const Feature& r) 
  { return RealPlacement(l)+Placement(r); }
inline Placement operator-(const Real& l, const Feature& r) 
  { return RealPlacement(l)-Placement(r); }
inline Placement operator*(const Real& l, const Feature& r) 
  { return RealPlacement(l)*Placement(r); }
inline Placement operator/(const Real& l, const Feature& r) 
  { return RealPlacement(l)/Placement(r); }

// binary(feature,real)
inline Placement operator+(const Feature& l, const Real& r) 
  { return Placement(l)+RealPlacement(r); }
inline Placement operator-(const Feature& l, const Real& r) 
  { return Placement(l)-RealPlacement(r); }
inline Placement operator*(const Feature& l, const Real& r) 
  { return Placement(l)*RealPlacement(r); }
inline Placement operator/(const Feature& l, const Real& r) 
  { return Placement(l)/RealPlacement(r); }

/**
 * This is an expression yielding a value suitable for use
 * as a placement. The only child features it can have are
 * other Measures (of any type), or Parameters, which 
 * are a kind of Measure anyway.
 *
 * Like other features, measures have a name and are owned
 * by a parent feature. A Measure needs a Placement of a
 * particular type, which may be resolved internally or
 * may require an external placement. Similarly its child
 * measures require placements which may be resolved 
 * internally (meaning up to the parent measure) or 
 * externally (meaning the parent's parent or higher).
 * A fully placed Measure's value is an expression
 * that can be evaluated at run time.
 */
class RealMeasure : public Feature {
public:
    explicit RealMeasure(const String& name);
    RealMeasure(const RealMeasure&);
    RealMeasure& operator=(const RealMeasure&);
    ~RealMeasure();

    static bool               isInstanceOf(const Feature&);
    static const RealMeasure& downcast(const Feature&);
    static RealMeasure&       downcast(Feature&);
protected:
    RealMeasure() { }
};


/**
 * A parameter is a RealMeasure with constant value. The value
 * must be a constant with respect to its parent Feature. That means it
 * is literally a constant, or it is a calculated value (a measure)
 * owned higher up the feature tree, i.e., from its parent's parent
 * or ancestors.
 *
 * Parameters cannot have subfeatures.
 */
class RealParameter : public RealMeasure {
public:
    explicit RealParameter(const String& name);
    RealParameter(const RealParameter&);
    RealParameter& operator=(const RealParameter&);
    ~RealParameter();

    static bool                 isInstanceOf(const Feature&);
    static const RealParameter& downcast(const Feature&);
    static RealParameter&       downcast(Feature&);
protected:
    RealParameter() { }
};

class Vec3Measure : public Feature {
public:
    explicit Vec3Measure(const String& name);
    Vec3Measure(const Vec3Measure&);
    Vec3Measure& operator=(const Vec3Measure&);
    ~Vec3Measure();

    static bool               isInstanceOf(const Feature&);
    static const Vec3Measure& downcast(const Feature&);
    static Vec3Measure&       downcast(Feature&);
protected:
    Vec3Measure() { }
};

class Vec3Parameter : public Vec3Measure {
public:
    explicit Vec3Parameter(const String& name);
    Vec3Parameter(const Vec3Parameter&);
    Vec3Parameter& operator=(const Vec3Parameter&);
    ~Vec3Parameter();

    static bool                 isInstanceOf(const Feature&);
    static const Vec3Parameter& downcast(const Feature&);
    static Vec3Parameter&       downcast(Feature&);
protected:
    Vec3Parameter() { }
};

class StationMeasure : public Feature {
public:
    explicit StationMeasure(const String& name);
    StationMeasure(const StationMeasure&);
    StationMeasure& operator=(const StationMeasure&);
    ~StationMeasure();

    static bool                  isInstanceOf(const Feature&);
    static const StationMeasure& downcast(const Feature&);
    static StationMeasure&       downcast(Feature&);
protected:
    StationMeasure() { }
};

class StationParameter : public StationMeasure {
public:
    explicit StationParameter(const String& name);
    StationParameter(const StationParameter&);
    StationParameter& operator=(const StationParameter&);
    ~StationParameter();

    static bool                    isInstanceOf(const Feature&);
    static const StationParameter& downcast(const Feature&);
    static StationParameter&       downcast(Feature&);
protected:
    StationParameter() { }
};

class Station : public Feature {
public:
    explicit Station(const String& name);
    Station(const String& name, const Vec3& defaultValue);
    Station(const Station&);
    Station& operator=(const Station&);
    ~Station();

    static bool           isInstanceOf(const Feature&);
    static const Station& downcast(const Feature&);
    static Station&       downcast(Feature&);
protected:
    Station() { }
};


class DirectionMeasure : public Feature {
public:
    explicit DirectionMeasure(const String& name);
    DirectionMeasure(const DirectionMeasure&);
    DirectionMeasure& operator=(const DirectionMeasure&);
    ~DirectionMeasure();

    static bool                  isInstanceOf(const Feature&);
    static const DirectionMeasure& downcast(const Feature&);
    static DirectionMeasure&       downcast(Feature&);
protected:
    DirectionMeasure() { }
};

// Sorry, no DirectionParameter (doesn't make sense to have one).

class Direction : public Feature {
public:
    explicit Direction(const String& name);
    Direction(const String& name, const Vec3& defaultValue);
    Direction(const Direction&);
    Direction& operator=(const Direction&);
    ~Direction();

    static bool             isInstanceOf(const Feature&);
    static const Direction& downcast(const Feature&);
    static Direction&       downcast(Feature&);
protected:
    Direction() { }
};

class OrientationMeasure : public Feature {
public:
    explicit OrientationMeasure(const String& name);
    OrientationMeasure(const OrientationMeasure&);
    OrientationMeasure& operator=(const OrientationMeasure&);
    ~OrientationMeasure();

    static bool                  isInstanceOf(const Feature&);
    static const OrientationMeasure& downcast(const Feature&);
    static OrientationMeasure&       downcast(Feature&);
protected:
    OrientationMeasure() { }
};

// Sorry, no OrientationParameter (doesn't make sense to have one).

class Orientation : public Feature {
public:
    explicit Orientation(const String& name);
    Orientation(const String& name, const Mat33& defaultValue);
    Orientation(const Orientation&);
    Orientation& operator=(const Orientation&);
    ~Orientation();

    const Direction& getAxis(int) const;
    const Direction&   x()        const {return getAxis(0);}
    const Direction&   y()        const {return getAxis(1);}
    const Direction&   z()        const {return getAxis(2);}

    static bool               isInstanceOf(const Feature&);
    static const Orientation& downcast(const Feature&);
    static Orientation&       downcast(Feature&);
protected:
    Orientation() { }
};

class Frame : public Feature {
public:
    explicit Frame(const String& name);
    Frame(const Frame&);
    Frame& operator=(const Frame&);
    ~Frame();

    const Station&     getOrigin() const;
    const Orientation& getOrientation() const;
    const Direction&   getAxis(int i) const {return getOrientation().getAxis(i);}
    const Direction&   x()            const {return getOrientation().x();}
    const Direction&   y()            const {return getOrientation().y();}
    const Direction&   z()            const {return getOrientation().z();}

    static bool         isInstanceOf(const Feature&);
    static const Frame& downcast(const Feature&);
    static Frame&       downcast(Feature&);
protected:
    Frame() { }
};

} // namespace simtk

#endif // SIMTK_FEATURE_H_
