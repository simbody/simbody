#ifndef SIMTK_BASIC_PLACEMENTS_H_
#define SIMTK_BASIC_PLACEMENTS_H_

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
 * Definitions for mundane Placement types.
 */

#include "simbody/internal/SimbodyCommon.h"
#include "simbody/internal/Geometry.h"
#include "simbody/internal/Mechanics.h"
#include "simbody/internal/Placement.h"

#include <iostream>

namespace simtk {

class Subsystem;
class Feature;
class Station;
class Direction;
class Orientation;
class FrameFeature;
class RealMeasure;
class   RealParameter;
class Vec3Measure;
class   Vec3Parameter;
class StationMeasure;
class   StationParameter;
class DirectionMeasure;
class OrientationMeasure;
class InertiaMeasure;

// Declared below.
class Placement;
class   RealPlacement;
class   Vec2Placement;
class   Vec3Placement;
class   StationPlacement;
class   DirectionPlacement;
class   OrientationPlacement;
class   FramePlacement;
class   PlacementList;

class PlacementValue;

class RealPlacement : public Placement {
public:
    RealPlacement() { }
    RealPlacement(const Real&);
    RealPlacement(const RealParameter&);
    RealPlacement(const RealMeasure&);
    explicit RealPlacement(const Feature&);
    explicit RealPlacement(const Placement&);

    static bool          canConvert(const Placement&);
    static RealPlacement convert(const Placement&);

    static bool                 isInstanceOf(const Placement&);

    static const RealPlacement& downcast(const Placement&);
    static RealPlacement&       downcast(Placement&);

    // For internal use only.
    explicit RealPlacement(class RealPlacementRep* r) 
      : Placement(reinterpret_cast<PlacementRep*>(r)) { }
    const RealPlacementRep& getRep() const
      { return *reinterpret_cast<const RealPlacementRep*>(rep); }
};

class Vec3Placement : public Placement {
public:
    Vec3Placement() { }
    Vec3Placement(const Vec3&);
    Vec3Placement(const Vec3Measure&);
    Vec3Placement(const Vec3Parameter&);
    explicit Vec3Placement(const Feature&);
    explicit Vec3Placement(const Placement&);

    static bool          canConvert(const Placement&);
    static Vec3Placement convert(const Placement&);

    // For internal use only.
    explicit Vec3Placement(class Vec3PlacementRep* r) 
      : Placement(reinterpret_cast<PlacementRep*>(r)) { }
    const Vec3PlacementRep& getRep() const
      { return *reinterpret_cast<const Vec3PlacementRep*>(rep); }    
    static bool                 isInstanceOf(const Placement&);
    static const Vec3Placement& downcast(const Placement&);
    static Vec3Placement&       downcast(Placement&);
};

class StationPlacement : public Placement {
public:
    StationPlacement() { }
    StationPlacement(const Station&);          // implicit conversion
    StationPlacement(const StationMeasure&);   // implicit conversion
    StationPlacement(const StationParameter&); // implicit conversion
    explicit StationPlacement(const Vec3&);
    explicit StationPlacement(const FrameFeature&);   // use the origin
    explicit StationPlacement(const Feature&);
    explicit StationPlacement(const Placement&);

    // For internal use only.
    explicit StationPlacement(class StationPlacementRep* r) 
      : Placement(reinterpret_cast<PlacementRep*>(r)) { }
    const StationPlacementRep& getRep() const
      { return *reinterpret_cast<const StationPlacementRep*>(rep); }    
    static bool                    isInstanceOf(const Placement&);
    static const StationPlacement& downcast(const Placement&);
    static StationPlacement&       downcast(Placement&);
};

class DirectionPlacement : public Placement {
public:
    DirectionPlacement() { }
    DirectionPlacement(const Direction&);           // implicit conversion
    DirectionPlacement(const DirectionMeasure&);    // implicit conversion

    explicit DirectionPlacement(const UnitVec3&);
    explicit DirectionPlacement(const Vec3&);
    explicit DirectionPlacement(const Feature&);
    explicit DirectionPlacement(const Placement&);

    DirectionPlacement(const Orientation&, int i);  // use the i'th axis
    DirectionPlacement(const FrameFeature&, int i);        // use the i'th axis

    // For internal use only.
    explicit DirectionPlacement(class DirectionPlacementRep* r) 
      : Placement(reinterpret_cast<PlacementRep*>(r)) { }
    const DirectionPlacementRep& getRep() const
      { return *reinterpret_cast<const DirectionPlacementRep*>(rep); } 
    static bool                      isInstanceOf(const Placement&);
    static const DirectionPlacement& downcast(const Placement&);
    static DirectionPlacement&       downcast(Placement&);
};

// Three, mutually orthogonal, right handed directions.
class OrientationPlacement : public Placement {
public:
    OrientationPlacement() { }
    OrientationPlacement(const Orientation&);           // implicit conversion
    OrientationPlacement(const OrientationMeasure&);    // implicit conversion

    explicit OrientationPlacement(const RotationMat&);
    explicit OrientationPlacement(const FrameFeature&);        // use the orientation
    explicit OrientationPlacement(const Feature&);
    explicit OrientationPlacement(const Placement&);

    // Returns the transpose of this rotation, which is the inverse rotation matrix.
    OrientationPlacement invert() const;

    // Create a set of axes with the indicated Z axis, but arbitrary X & Y.
    static OrientationPlacement createFromZAxis(const DirectionPlacement& z);

    // For internal use only.
    explicit OrientationPlacement(class OrientationPlacementRep* r) 
      : Placement(reinterpret_cast<PlacementRep*>(r)) { }
    const OrientationPlacementRep& getRep() const
      { return *reinterpret_cast<const OrientationPlacementRep*>(rep); } 
    static bool                        isInstanceOf(const Placement&);
    static const OrientationPlacement& downcast(const Placement&);
    static OrientationPlacement&       downcast(Placement&);
};

// A symmetric matrix.
class InertiaPlacement : public Placement {
public:
    InertiaPlacement() { }
    InertiaPlacement(const InertiaMeasure&);    // implicit conversion

    // Construct the inertia of a point mass located at p having mass m.
    InertiaPlacement(const StationPlacement& p, const RealPlacement& m);

    // Construct from principal moments.
    InertiaPlacement(const RealPlacement& Ixx, const RealPlacement& Iyy, const RealPlacement& Izz);

    // Construct from general inertia terms. Note the order.
    InertiaPlacement(const RealPlacement& Ixx, const RealPlacement& Iyy, const RealPlacement& Izz,
                     const RealPlacement& Ixy, const RealPlacement& Ixz, const RealPlacement& Iyz);

    explicit InertiaPlacement(const InertiaMat&);
    explicit InertiaPlacement(const Feature&);
    explicit InertiaPlacement(const Placement&);

    InertiaPlacement changeAxes(const OrientationPlacement&) const;
    InertiaPlacement shiftFromCOM(const StationPlacement& to,
                                  const RealPlacement&    totalMass) const;
    InertiaPlacement shiftToCOM(const StationPlacement& to,
                                const RealPlacement&    totalMass) const;

    // For internal use only.
    explicit InertiaPlacement(class InertiaPlacementRep* r) 
      : Placement(reinterpret_cast<PlacementRep*>(r)) { }
    const InertiaPlacementRep& getRep() const
      { return *reinterpret_cast<const InertiaPlacementRep*>(rep); } 
    static bool                    isInstanceOf(const Placement&);
    static const InertiaPlacement& downcast(const Placement&);
    static InertiaPlacement&       downcast(Placement&);
};

class FramePlacement : public Placement {
public:
    FramePlacement() { }
    FramePlacement(const FrameFeature&);                       // implicit conversion

    explicit FramePlacement(const Frame&);

    // Inherit the orientation from the Station feature's parent.
    explicit FramePlacement(const Station&);

    // Inherit the origin from the Orientation feature's parent.
    explicit FramePlacement(const Orientation&);

    explicit FramePlacement(const Feature&);
    explicit FramePlacement(const Placement&);

    FramePlacement(const OrientationPlacement&, const StationPlacement&);

    // For internal use only.
    explicit FramePlacement(class FramePlacementRep* r) 
      : Placement(reinterpret_cast<PlacementRep*>(r)) { }
    const FramePlacementRep& getRep() const
      { return *reinterpret_cast<const FramePlacementRep*>(rep); } 
    static bool                  isInstanceOf(const Placement&);
    static const FramePlacement& downcast(const Placement&);
    static FramePlacement&       downcast(Placement&);
};

} // namespace simtk

#endif // SIMTK_BASIC_PLACEMENTS_H_
