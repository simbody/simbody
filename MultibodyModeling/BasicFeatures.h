#ifndef SIMTK_BASIC_FEATURES_H_
#define SIMTK_BASIC_FEATURES_H_

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
 * Client side declarations of run-of-the-mill basic Features.
 */


#include "simbody/SimbodyCommon.h"
#include "Feature.h"

#include "simtk/SimTK.h"
#include "simmatrix/SmallMatrix.h"

#include <iostream>
#include <cassert>

namespace simtk {

class Placement;

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

    const RealPlacement& getPlacement() const;
    const Real&          getValue()     const;

    static bool               isInstanceOf(const Subsystem&);
    static const RealMeasure& downcast(const Subsystem&);
    static RealMeasure&       downcast(Subsystem&);
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

    const RealPlacement& getPlacement() const;
    const Real&          getValue()     const;

    static bool                 isInstanceOf(const Subsystem&);
    static const RealParameter& downcast(const Subsystem&);
    static RealParameter&       downcast(Subsystem&);
protected:
    RealParameter() { }
};

class Vec3Measure : public Feature {
public:
    explicit Vec3Measure(const String& name);
    Vec3Measure(const Vec3Measure&);
    Vec3Measure& operator=(const Vec3Measure&);
    ~Vec3Measure();

    const Vec3Placement& getPlacement() const;
    const Vec3&          getValue()     const;

    static bool               isInstanceOf(const Subsystem&);
    static const Vec3Measure& downcast(const Subsystem&);
    static Vec3Measure&       downcast(Subsystem&);
protected:
    Vec3Measure() { }
};

class Vec3Parameter : public Vec3Measure {
public:
    explicit Vec3Parameter(const String& name);
    Vec3Parameter(const Vec3Parameter&);
    Vec3Parameter& operator=(const Vec3Parameter&);
    ~Vec3Parameter();

    const Vec3Placement& getPlacement() const;
    const Vec3&          getValue()     const;

    static bool                 isInstanceOf(const Subsystem&);
    static const Vec3Parameter& downcast(const Subsystem&);
    static Vec3Parameter&       downcast(Subsystem&);
protected:
    Vec3Parameter() { }
};

class StationMeasure : public Feature {
public:
    explicit StationMeasure(const String& name);
    StationMeasure(const StationMeasure&);
    StationMeasure& operator=(const StationMeasure&);
    ~StationMeasure();

    const StationPlacement& getPlacement() const;
    const Vec3&             getValue()     const;

    static bool                  isInstanceOf(const Subsystem&);
    static const StationMeasure& downcast(const Subsystem&);
    static StationMeasure&       downcast(Subsystem&);
protected:
    StationMeasure() { }
};

class StationParameter : public StationMeasure {
public:
    explicit StationParameter(const String& name);
    StationParameter(const StationParameter&);
    StationParameter& operator=(const StationParameter&);
    ~StationParameter();

    const StationPlacement& getPlacement() const;
    const Vec3&             getValue()     const;

    static bool                    isInstanceOf(const Subsystem&);
    static const StationParameter& downcast(const Subsystem&);
    static StationParameter&       downcast(Subsystem&);
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

    const StationPlacement& getPlacement() const;
    const Vec3&             getValue()     const;

    static bool           isInstanceOf(const Subsystem&);
    static const Station& downcast(const Subsystem&);
    static Station&       downcast(Subsystem&);
protected:
    Station() { }
};


class DirectionMeasure : public Feature {
public:
    explicit DirectionMeasure(const String& name);
    DirectionMeasure(const DirectionMeasure&);
    DirectionMeasure& operator=(const DirectionMeasure&);
    ~DirectionMeasure();

    const DirectionPlacement& getPlacement() const;
    const UnitVec3&           getValue()     const;

    static bool                    isInstanceOf(const Subsystem&);
    static const DirectionMeasure& downcast(const Subsystem&);
    static DirectionMeasure&       downcast(Subsystem&);
protected:
    DirectionMeasure() { }
};

// Sorry, no DirectionParameter (doesn't make sense to have one).
// If you want one, it should be achieved by rotating a DirectionFeature
// using a RealParameter representing an angle.

class Direction : public Feature {
public:
    explicit Direction(const String& name);
    Direction(const String& name, const Vec3& defaultValue);    // normalizes automatically
    Direction(const Direction&);
    Direction& operator=(const Direction&);
    ~Direction();

    const DirectionPlacement& getPlacement() const;
    const UnitVec3&           getValue()     const;

    static bool             isInstanceOf(const Subsystem&);
    static const Direction& downcast(const Subsystem&);
    static Direction&       downcast(Subsystem&);
protected:
    Direction() { }
};

class OrientationMeasure : public Feature {
public:
    explicit OrientationMeasure(const String& name);
    OrientationMeasure(const OrientationMeasure&);
    OrientationMeasure& operator=(const OrientationMeasure&);
    ~OrientationMeasure();

    const OrientationPlacement& getPlacement() const;
    const MatRotation&          getValue()     const;

    static bool                      isInstanceOf(const Subsystem&);
    static const OrientationMeasure& downcast(const Subsystem&);
    static OrientationMeasure&       downcast(Subsystem&);
protected:
    OrientationMeasure() { }
};

// Sorry, no OrientationParameter (doesn't make sense to have one).

class Orientation : public Feature {
public:
    explicit Orientation(const String& name);
    Orientation(const String& name, const MatRotation& defaultValue);
    Orientation(const Orientation&);
    Orientation& operator=(const Orientation&);
    ~Orientation();

    const OrientationPlacement& getPlacement() const;
    const MatRotation&          getValue()     const;

    const Direction& getAxis(int) const;
    const Direction& x()          const {return getAxis(0);}
    const Direction& y()          const {return getAxis(1);}
    const Direction& z()          const {return getAxis(2);}

    static bool               isInstanceOf(const Subsystem&);
    static const Orientation& downcast(const Subsystem&);
    static Orientation&       downcast(Subsystem&);
protected:
    Orientation() { }
};

class InertiaMeasure : public Feature {
public:
    explicit InertiaMeasure(const String& name);
    InertiaMeasure(const InertiaMeasure&);
    InertiaMeasure& operator=(const InertiaMeasure&);
    ~InertiaMeasure();

    const InertiaPlacement& getPlacement() const;
    const MatInertia&       getValue()     const;

    static bool                  isInstanceOf(const Subsystem&);
    static const InertiaMeasure& downcast(const Subsystem&);
    static InertiaMeasure&       downcast(Subsystem&);
protected:
    InertiaMeasure() { }
};

class FrameFeature : public Feature {
public:
    explicit FrameFeature(const String& name);
    FrameFeature(const FrameFeature&);
    FrameFeature& operator=(const FrameFeature&);
    ~FrameFeature();

    const FramePlacement& getPlacement() const;
    const Mat34&          getValue()     const;

    const Station&     getOrigin()      const;
    const Orientation& getOrientation() const;
    const Direction&   getAxis(int i)   const {return getOrientation().getAxis(i);}
    const Direction&   x()              const {return getOrientation().x();}
    const Direction&   y()              const {return getOrientation().y();}
    const Direction&   z()              const {return getOrientation().z();}

    static bool                isInstanceOf(const Subsystem&);
    static const FrameFeature& downcast(const Subsystem&);
    static FrameFeature&       downcast(Subsystem&);
protected:
    FrameFeature() { }
};

} // namespace simtk

#endif // SIMTK_BASIC_FEATURES_H_
