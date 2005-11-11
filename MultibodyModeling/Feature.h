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

namespace simtk {

// Declared below.
class Feature;
class RealParameter;
class StationParameter;
class RealMeasure;
class StationMeasure;
class Station;
class Direction;
class Orientation;
class Frame;


/**
 * Features form a tree because many Features have sub-Features (children).
 * Parent features own their children; destructing the parent destructs
 * all the children).
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

    bool           hasParentFeature() const;
    int            getIndexInParent() const; // -1 if no parent
    const Feature& getParentFeature() const;

    String getName()            const;
    String getFullName()        const; // "ancestors.parent.name"
    String getFeatureTypeName() const; // "Station", "Frame", etc.

    const Feature&          getFeature         (const String&) const; // generic
    const RealParameter&    getRealParameter   (const String&) const; // type checked
    const StationParameter& getStationParameter(const String&) const; //   "
    const RealMeasure&      getRealMeasure     (const String&) const;
    const StationMeasure&   getStationMeasure  (const String&) const;
    const Station&          getStation         (const String&) const;
    const Direction&        getDirection       (const String&) const;
    const Orientation&      getOrientation     (const String&) const;
    const Frame&            getFrame           (const String&) const;

    Feature&          updFeature         (const String&);   // generic
    RealParameter&    updRealParameter   (const String&);   // type checked
    StationParameter& updStationParameter(const String&);   //   "
    RealMeasure&      updRealMeasure     (const String&);
    StationMeasure&   updStationMeasure  (const String&);
    Station&          updStation         (const String&);
    Direction&        updDirection       (const String&);
    Orientation&      updOrientation     (const String&);
    Frame&            updFrame           (const String&);

    RealParameter&    addRealParameter   (const String&);
    StationParameter& addStationParameter(const String&);
    RealMeasure&      addRealMeasure     (const String&);
    StationMeasure&   addStationMeasure  (const String&);
    Station&          addStation         (const String&);
    Direction&        addDirection       (const String&);
    Orientation&      addOrientation     (const String&);
    Frame&            addFrame           (const String&);

    Feature&          addFeatureLike         (const Feature&,          const String&);
    RealParameter&    addRealParameterLike   (const RealParameter&,    const String&);
    StationParameter& addStationParameterLike(const StationParameter&, const String&);
    RealMeasure&      addRealMeasureLike     (const RealMeasure&,      const String&);
    StationMeasure&   addStationMeasureLike  (const StationMeasure&,   const String&);
    Station&          addStationLike         (const Station&,          const String&);
    Direction&        addDirectionLike       (const Direction&,        const String&);
    Orientation&      addOrientationLike     (const Orientation&,      const String&);
    Frame&            addFrameLike           (const Frame&,            const String&);

    // Subfeatures of this feature
    int            getNFeatures() const;
    const Feature& getFeature(int) const;
    Feature&       updFeature(int);
    void           placeFeature(const String&, const Placement&);
    void           placeFeature(int, const Placement&);

    // true if Feature & all children have been placed
    bool hasPlacement() const;
    void setPlacement(const Placement&);
    const Placement& getPlacement() const;

    String toString(const String& linePrefix="") const;

protected:
    class FeatureRep* rep;
    friend class FeatureRep;
};
std::ostream& operator<<(std::ostream& o, const Feature&);

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

    void set(const Placement&);

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

    void set(const Real&);

    static bool                 isInstanceOf(const Feature&);
    static const RealParameter& downcast(const Feature&);
    static RealParameter&       downcast(Feature&);
};

class StationMeasure : public Feature {
public:
    explicit StationMeasure(const String& name);
    StationMeasure(const StationMeasure&);
    StationMeasure& operator=(const StationMeasure&);
    ~StationMeasure();

    void set(const Placement&);

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

    void set(const Vec3&);

    static bool                    isInstanceOf(const Feature&);
    static const StationParameter& downcast(const Feature&);
    static StationParameter&       downcast(Feature&);
};

class Station : public Feature {
public:
    explicit Station(const String& name);
    Station(const String& name, const Vec3& defaultValue);
    Station(const Station&);
    Station& operator=(const Station&);
    ~Station();
    
    void placeStation(const StationPlacement&);

    static bool           isInstanceOf(const Feature&);
    static const Station& downcast(const Feature&);
    static Station&       downcast(Feature&);
};

class Direction : public Feature {
public:
    explicit Direction(const String& name);
    Direction(const String& name, const Vec3& defaultValue);
    Direction(const Direction&);
    Direction& operator=(const Direction&);
    ~Direction();

    void placeDirection(const DirectionPlacement&);

    static bool             isInstanceOf(const Feature&);
    static const Direction& downcast(const Feature&);
    static Direction&       downcast(Feature&);
};

class Orientation : public Feature {
public:
    explicit Orientation(const String& name);
    Orientation(const String& name, const Mat33& defaultValue);
    Orientation(const Orientation&);
    Orientation& operator=(const Orientation&);
    ~Orientation();

    const Direction& getAxis(int) const;

    void placeOrientation(const OrientationPlacement&);

    static bool               isInstanceOf(const Feature&);
    static const Orientation& downcast(const Feature&);
    static Orientation&       downcast(Feature&);
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
    const Direction&   x()            const {return getOrientation().getAxis(0);}
    const Direction&   y()            const {return getOrientation().getAxis(1);}
    const Direction&   z()            const {return getOrientation().getAxis(2);}

    void placeFrame(const FramePlacement&);

    static bool         isInstanceOf(const Feature&);
    static const Frame& downcast(const Feature&);
    static Frame&       downcast(Feature&);
protected:
    Frame() { }
};

} // namespace simtk

#endif // SIMTK_FEATURE_H_
