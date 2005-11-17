#ifndef SIMTK_PLACEMENT_H_
#define SIMTK_PLACEMENT_H_

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
 * Placements definitions.
 */

#include "SimbodyCommon.h"

#include <iostream>

namespace simtk {

class Feature;
class Station;
class Direction;
class Orientation;
class Frame;
class RealMeasure;
class   RealParameter;
class Vec3Measure;
class   Vec3Parameter;
class StationMeasure;
class   StationParameter;

// Declared below.
class Placement;
class FeaturePlacement;
class RealPlacement;
class Vec2Placement;
class Vec3Placement;
class StationPlacement;
class DirectionPlacement;
class OrientationPlacement;
class FramePlacement;


// 
// <placement>     ::= constant | parameter | feature.placement
//                     | func(<placement> ...)
//
// Features provide a list of the placement types they need. Each of
// these list elements must be associated with a <placementExpr> eventually.
// Then these are evaluated at the appropriate runtime stage.
//
// <placement>'s have value types scalar, station, direction,
// orientation, or frame.
//                     

class Placement {
public:
    Placement() : rep(0) { }
    Placement(const Placement&);
    Placement& operator=(const Placement&);
    ~Placement();

    Placement(const Feature&);  // implicit conversion to a FeaturePlacement
    Placement(const Real&);     // implicit conversion to a RealPlacement
    Placement(const Vec3&);     // implicit conversion to a Vec3Placement

    bool           hasOwner() const;
    const Feature& getOwner() const;
    int            getIndexInOwner() const;

    bool isConstant() const;  // a plain old numerical value?
    bool dependsOn(const Feature&) const; // recursive dependency check

    String toString(const String& linePrefix="") const;

    // For internal use only.
    bool                      hasRep() const {return rep != 0;}
    const class PlacementRep& getRep() const {assert(rep); return *rep;}
    class PlacementRep&       updRep()       {assert(rep); return *rep;}
    void                      setRep(PlacementRep* pp) {assert(!rep); rep=pp;}
protected:
    class PlacementRep* rep;
    friend class PlacementRep;
};

// Global operators involving Placements.
std::ostream& operator<<(std::ostream& o, const Placement&);

// unary
Placement operator+(const Placement& f);
Placement operator-(const Placement& f);
RealPlacement      length(const Placement& f);
DirectionPlacement normalize(const Placement& f);

// binary
Placement operator+(const Placement& l, const Placement& r); 
Placement operator-(const Placement& l, const Placement& r); 
Placement operator*(const Placement& l, const Placement& r); 
Placement operator/(const Placement& l, const Placement& r);

/**
 * Create a Placement which is evaluated by returning the
 * value of the indicated Feature's Placement, with an
 * optional index selecting some subcomponent of the
 * Placement (e.g., selecting a particular axis from
 * a Frame). Resolution of this reference is
 * deferred until runtime; the Feature doesn't even need
 * to *have* a placement at the time we create this
 * reference.
 */
class FeaturePlacement : public Placement {
public:
    FeaturePlacement() { }
    explicit FeaturePlacement(const Feature&);
    FeaturePlacement(const Feature&, int index);

    static bool                    isInstanceOf(const Placement&);
    static const FeaturePlacement& downcast(const Placement&);
    static FeaturePlacement&       downcast(Placement&);
};

class RealPlacement : public Placement {
public:
    RealPlacement() { }
    RealPlacement(const Real&);
    RealPlacement(const RealParameter&);
    RealPlacement(const RealMeasure&);
    RealPlacement(const Feature&);

    static RealPlacement negate(const RealPlacement&);
    static RealPlacement plus  (const RealPlacement& l,
                                const RealPlacement& r);
    static RealPlacement minus (const RealPlacement& l,
                                const RealPlacement& r);
    static RealPlacement times (const RealPlacement& l,
                                const RealPlacement& r);
    static RealPlacement divide(const RealPlacement& l,
                                const RealPlacement& r);
    static RealPlacement length(const Vec3Placement&);

    static bool                 isInstanceOf(const Placement&);
    static const RealPlacement& downcast(const Placement&);
    static RealPlacement&       downcast(Placement&);
};

class Vec3Placement : public Placement {
public:
    Vec3Placement() { }
    Vec3Placement(const Vec3&);
    Vec3Placement(const Vec3Parameter&);
    Vec3Placement(const Vec3Measure&);
    //Vec3Placement(const Feature&);

    static Vec3Placement plus  (const Placement& l,
                                const Placement& r);
    static Vec3Placement minus (const Placement& l,
                                const Placement& r);
    static Vec3Placement scale (const Placement& l,
                                const Placement& r);
    static Vec3Placement cast  (const Placement& s);

    static bool                 isInstanceOf(const Placement&);
    static const Vec3Placement& downcast(const Placement&);
    static Vec3Placement&       downcast(Placement&);
};

class StationPlacement : public Placement {
public:
    StationPlacement() { }
    StationPlacement(const Vec3&);    // implicit conversion
    StationPlacement(const Station&); //   "
    StationPlacement(const Frame&);   //   "
    StationPlacement(const Feature&);

    static StationPlacement   plus (const StationPlacement&,
                                    const Vec3Placement&);
    static StationPlacement   minus(const StationPlacement&,
                                    const Vec3Placement&);
    static StationPlacement   cast(const Vec3Placement&);

    static bool                    isInstanceOf(const Placement&);
    static const StationPlacement& downcast(const Placement&);
    static StationPlacement&       downcast(Placement&);
};

class DirectionPlacement : public Placement {
public:
    DirectionPlacement() { }
    DirectionPlacement(const Vec3&);      // implicit conversion
    DirectionPlacement(const Direction&); //   "
    DirectionPlacement(const Feature&);

    static DirectionPlacement normalize(const Vec3Placement& v);
    static DirectionPlacement normalize(const StationPlacement& v);

    static bool                      isInstanceOf(const Placement&);
    static const DirectionPlacement& downcast(const Placement&);
    static DirectionPlacement&       downcast(Placement&);
};

// Three, mutually orthogonal, right handed directions.
class OrientationPlacement : public Placement {
public:
    OrientationPlacement() { }
    OrientationPlacement(const Mat33&);      // implicit conversion
    OrientationPlacement(const Feature&);

    static bool                        isInstanceOf(const Placement&);
    static const OrientationPlacement& downcast(const Placement&);
    static OrientationPlacement&       downcast(Placement&);
};

class FramePlacement : public Placement {
public:
    FramePlacement() { }
    FramePlacement(const Frame&);   // implicit conversions
    FramePlacement(const Station&); //   orientation inherited from 
                                    //   Station's parent
    FramePlacement(const Orientation&); // origin inherited from
                                        // Orientation's parent
    FramePlacement(const Orientation&, const Station&);
    FramePlacement(const Feature&);

    static bool                  isInstanceOf(const Placement&);
    static const FramePlacement& downcast(const Placement&);
    static FramePlacement&       downcast(Placement&);
};

} // namespace simtk

#endif // SIMTK_PLACEMENT_H_
