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

class Subsystem;
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
class DirectionMeasure;
class OrientationMeasure;

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


// 
// Features provide a list of the placement types they need. Each of
// these list elements must be associated with a <placementExpr> eventually.
// Then these are evaluated at the appropriate runtime stage.
//
// <placement>'s have value types real, vec2, vec3, station, direction,
// orientation, frame, list.
//                     

class Placement {
public:
    Placement() : rep(0) { }
    Placement(const Placement&);
    Placement& operator=(const Placement&);
    ~Placement();

    // Implicit conversions to a Placement of appropriate concrete type.
    Placement(const Real&);
    Placement(const Vec3&);
    Placement(const Mat33&);

    // These create a "feature placement" (reference to the placement
    // of a feature) of a particular type, but starting with a generic
    // Subsystem, which must turn out to be a Feature. Note that
    // it is not necessary for the Feature to have a placement at the time
    // we create this reference; it need only have one by the time we
    // require a numerical value.
    // The first constructor below is also an implicit
    // conversion, so any Feature may be used as a Placement. The second
    // still references the indicated Feature (not one of its subfeatures)
    // but then selects a subcomponent of the resulting Placement as 
    // indicated by the index.
    Placement(const Subsystem&);
    Placement(const Subsystem&, int index);

    void realize(/*State,*/ Stage) const;

    // This will throw an exception if this Placement's value has not 
    // already been realize()'d.
    const PlacementValue& getValue() const;

    bool             hasOwner() const;
    const Subsystem& getOwner() const;
    int              getIndexInOwner() const;

    bool isConstant() const;  // a plain old numerical value?
    bool dependsOn(const Feature&) const; // recursive dependency check

    String toString(const String& linePrefix="") const;

    // For internal use only.
    explicit Placement(class PlacementRep*);
    bool                      hasRep() const            {return rep != 0;}
    const class PlacementRep& getRep() const            {assert(rep); return *rep;}
    class PlacementRep&       updRep()                  {assert(rep); return *rep;}
    void                      setRep(PlacementRep* pp)  {assert(!rep); rep=pp;}
    void checkPlacementConsistency(const Subsystem* expOwner, 
                                   int              expIndexInOwner,
                                   const Subsystem& expRoot) const;
protected:
    class PlacementRep* rep;
    friend class PlacementRep;
};

// Global operators involving Placements. Note that these actually
// represent families of operators overloaded based on their argument types.

std::ostream& operator<<(std::ostream& o, const Placement&);

// unary
Placement negate    (const Placement&);
Placement abs       (const Placement&);
Placement sqrt      (const Placement&);
Placement exp       (const Placement&);
Placement log       (const Placement&);
Placement sin       (const Placement&);
Placement cos       (const Placement&);
Placement asin      (const Placement&);
Placement acos      (const Placement&);
Placement length    (const Placement&);
Placement normalize (const Placement&);

// binary
Placement add       (const Placement& l, const Placement& r); 
Placement subtract  (const Placement& l, const Placement& r); 
Placement multiply  (const Placement& l, const Placement& r); 
Placement divide    (const Placement& l, const Placement& r);
Placement distance  (const Placement& l, const Placement& r);
Placement angle     (const Placement& l, const Placement& r);
Placement dot       (const Placement& l, const Placement& r);
Placement cross     (const Placement& l, const Placement& r);

// alternative access to some of above routines via overloaded operators
inline Placement operator- (const Placement& p)                     {return negate(p);}
inline Placement operator+ (const Placement& l, const Placement& r) {return add(l,r);} 
inline Placement operator- (const Placement& l, const Placement& r) {return subtract(l,r);}
inline Placement operator* (const Placement& l, const Placement& r) {return multiply(l,r);}
inline Placement operator/ (const Placement& l, const Placement& r) {return divide(l,r);}

class RealPlacement : public Placement {
public:
    RealPlacement() { }
    RealPlacement(const Real&);
    RealPlacement(const RealParameter&);
    RealPlacement(const RealMeasure&);
    explicit RealPlacement(const Feature&);

    // For internal use only.
    explicit RealPlacement(class RealPlacementRep* r) 
      : Placement(reinterpret_cast<PlacementRep*>(r)) { }
    const RealPlacementRep& getRep() const
      { return *reinterpret_cast<const RealPlacementRep*>(rep); }
    static bool                 isInstanceOf(const Placement&);
    static const RealPlacement& downcast(const Placement&);
    static RealPlacement&       downcast(Placement&);
};

class Vec3Placement : public Placement {
public:
    Vec3Placement() { }
    Vec3Placement(const Vec3&);
    Vec3Placement(const Vec3Measure&);
    Vec3Placement(const Vec3Parameter&);
    explicit Vec3Placement(const Feature&);

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
    explicit StationPlacement(const Frame&);   // use the origin
    explicit StationPlacement(const Feature&);

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

    explicit DirectionPlacement(const Vec3&);
    explicit DirectionPlacement(const Feature&);

    DirectionPlacement(const Orientation&, int i);  // use the i'th axis
    DirectionPlacement(const Frame&, int i);        // use the i'th axis

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

    explicit OrientationPlacement(const Mat33&);
    explicit OrientationPlacement(const Frame&);        // use the orientation
    explicit OrientationPlacement(const Feature&);

    // For internal use only.
    explicit OrientationPlacement(class OrientationPlacementRep* r) 
      : Placement(reinterpret_cast<PlacementRep*>(r)) { }
    const OrientationPlacementRep& getRep() const
      { return *reinterpret_cast<const OrientationPlacementRep*>(rep); } 
    static bool                        isInstanceOf(const Placement&);
    static const OrientationPlacement& downcast(const Placement&);
    static OrientationPlacement&       downcast(Placement&);
};

class FramePlacement : public Placement {
public:
    FramePlacement() { }
    FramePlacement(const Frame&);                       // implicit conversion

    // Inherit the orientation from the Station feature's parent.
    explicit FramePlacement(const Station&);

    // Inherit the origin from the Orientation feature's parent.
    explicit FramePlacement(const Orientation&);

    FramePlacement(const Orientation&, const Station&);

    //explicit FramePlacement(const Placement&);

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

#endif // SIMTK_PLACEMENT_H_
