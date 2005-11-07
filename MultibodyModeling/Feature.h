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
class Parameter;
class Measure;
class Station;
class Direction;
class Orientation;
class Frame;


/**
 * Features form a tree because many Features have child Features.
 * Parent features own their children (destructing the parent 
 * destructs all the children).
 *
 * Most features require placement in order to be useful (e.g.,
 * a Station has to have a location). Placement expressions can
 * be constants or can involve parent or ancestor placement
 * expressions, but not children or siblings.
 *
 * Every feature has a name and an index by which it is known
 * to its parent. Parentless features can exist but they can't
 * be placed. Parentless features still have a name but they do
 * not have an index.
 */
class Feature {
public:
    Feature() : rep(0) { }
    Feature(const Feature&);    // external placements are not copied
    Feature& operator=(const Feature&);
    ~Feature();

    String         getName()          const; // name and index as known to parent
    int            getIndexInParent() const; // -1 if no parent
    const Feature& getParentFeature() const;

    String getFullName() const; // "ancestors/parent/name"
    String getFeatureTypeName() const; // "Station", "Frame", etc.

    const Feature& getFeature    (const String&) const;
    Feature&       updFeature    (const String&);
    Feature&       addFeatureLike(const Feature&, const String&);

    const Parameter&   getParameter  (const String&) const;
    const Measure&     getMeasure    (const String&) const;
    const Station&     getStation    (const String&) const;
    const Direction&   getDirection  (const String&) const;
    const Orientation& getOrientation(const String&) const;
    const Frame&       getFrame      (const String&) const;

    Parameter&   updParameter  (const String&);
    Measure&     updMeasure    (const String&);
    Station&     updStation    (const String&);
    Direction&   updDirection  (const String&);
    Orientation& updOrientation(const String&);
    Frame&       updFrame      (const String&);

    Parameter&   addParameter  (const String&);
    Measure&     addMeasure    (const String&);
    Station&     addStation    (const String&);
    Direction&   addDirection  (const String&);
    Orientation& addOrientation(const String&);
    Frame&       addFrame      (const String&);

    Parameter&   addParameterLike  (const Parameter&,   const String&);
    Measure&     addMeasureLike    (const Measure&,     const String&);
    Station&     addStationLike    (const Station&,     const String&);
    Direction&   addDirectionLike  (const Direction&,   const String&);
    Orientation& addOrientationLike(const Orientation&, const String&);
    Frame&       addFrameLike      (const Frame&,       const String&);

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
 * A parameter is a feature with a numerical value. The value
 * must be a constant with respect to its owner. That means it
 * is literally a constant, or it is a calculated value (a measure)
 * owned higher up the feature tree, i.e., from its owner's parent
 * or ancestors.
 *
 * Parameters cannot have subfeatures.
 */
class Parameter : public Feature {
public:
    explicit Parameter(const String& name);
    Parameter(const Parameter&);
    Parameter& operator=(const Parameter&);
    ~Parameter();

    void set(const Real&);

    static bool             isInstanceOf(const Feature&);
    static const Parameter& downcast(const Feature&);
    static Parameter&       downcast(Feature&);
};

/**
 * This is an expression yielding a value suitable for use
 * as a placement. The only child features it can have are
 * parameters.
 *
 * Like other features, measures have a name and are owned
 * by some parent feature. A Measure's value is an expression
 * that can be evaluated at run time.
 */
class Measure : public Feature {
public:
    explicit Measure(const String& name);
    Measure(const Measure&);
    Measure& operator=(const Measure&);
    ~Measure();

    void set(const Placement&);

    static bool           isInstanceOf(const Feature&);
    static const Measure& downcast(const Feature&);
    static Measure&       downcast(Feature&);
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
std::ostream& operator<<(std::ostream& o, const Station&);

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
std::ostream& operator<<(std::ostream& o, const Direction&);

class Orientation : public Feature {
public:
    explicit Orientation(const String& name);
    Orientation(const String& name, const Mat33& defaultValue);
    Orientation(const Orientation&);
    Orientation& operator=(const Orientation&);
    ~Orientation();

    void placeOrientation(const OrientationPlacement&);

    static bool               isInstanceOf(const Feature&);
    static const Orientation& downcast(const Feature&);
    static Orientation&       downcast(Feature&);
};
std::ostream& operator<<(std::ostream& o, const Orientation&);

class Frame : public Feature {
public:
    explicit Frame(const String& name);
    Frame(const Frame&);
    Frame& operator=(const Frame&);
    ~Frame();

    const Station& getOrigin() const;
    const Direction& getAxis(int) const;
    const Direction& x() const;
    const Direction& y() const;
    const Direction& z() const;

    void placeFrame(const FramePlacement&);

    static bool         isInstanceOf(const Feature&);
    static const Frame& downcast(const Feature&);
    static Frame&       downcast(Feature&);
};
//std::ostream& operator<<(std::ostream& o, const Frame&);



} // namespace simtk

#endif // SIMTK_FEATURE_H_
