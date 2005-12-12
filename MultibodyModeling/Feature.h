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

// Declared below and elsewhere. Indenting shows inheritance structure.
class Subsystem;
class   Feature;
class       Station;
class       Direction;
class       Orientation;
class       Frame;
class           Body;
class               RigidBody;
class               DeformableBody;
class       RealMeasure;
class           RealParameter;
class       StationMeasure;
class           StationParameter;
class       DirectionMeasure;
class       OrientationMeasure;
class       MassElement;
class           PointMassElement;
class           SphereMassElement;
class           CylinderMassElement;
class           BrickMassElement;
class   Joint;
class       PinJoint;
class       BallJoint;
class   Multibody;

/**
 * Subsystems form a tree because many of them have sub-Subsystems (children).
 * Parent subsystems own their children; i.e., destructing the parent destructs
 * all the children.
 *
 * Low-level "elemental" subsystems are called Features. Features require
 * a Placement in order to be useful (e.g., a Station has to have a location).
 * A Feature's Placement may have a constant value, or can involve sibling,
 * parent, or ancestor Features but not its own children. 
 *
 * Every Subsystem has a name and an index by which it is known
 * to its parent. Parentless subsystems can exist; they still have a
 * name but they do not have an index. We call these "root" or "prototype"
 * subsystems depending on intended use.
 */
class Subsystem {
public:
    Subsystem() : rep(0) { }
    Subsystem(const Subsystem&);    // external placements are not copied or assigned
    Subsystem& operator=(const Subsystem&);
    ~Subsystem();

    // calculate values for all fully-resolved placements
    void realize(/*State,*/ Stage) const;

    const Placement&      getPlacement()  const; // fails if Subsystem is not a Feature
    const PlacementValue& getValue()      const; //   "
    void place(const Placement&);                //   "

    // Read-only access to child subsystems.
    const Subsystem&        getSubsystem       (const String&) const; // generic
    const Feature&          getFeature         (const String&) const; // slightly less generic
    const RealParameter&    getRealParameter   (const String&) const; // type checked
    const Vec3Parameter&    getVec3Parameter   (const String&) const; //   "
    const StationParameter& getStationParameter(const String&) const;
    const RealMeasure&      getRealMeasure     (const String&) const;
    const Vec3Measure&      getVec3Measure     (const String&) const;
    const StationMeasure&   getStationMeasure  (const String&) const;
    const Station&          getStation         (const String&) const;
    const Direction&        getDirection       (const String&) const;
    const Orientation&      getOrientation     (const String&) const;
    const Frame&            getFrame           (const String&) const;

    // Writable access to subsystems, e.g. allowing feature placement.
    Subsystem&              updSubsystem       (const String&);   // generic
    Feature&                updFeature         (const String&);   // less generic
    RealParameter&          updRealParameter   (const String&);   // type checked
    Vec3Parameter&          updVec3Parameter   (const String&);   //   "
    StationParameter&       updStationParameter(const String&);
    RealMeasure&            updRealMeasure     (const String&);
    Vec3Measure&            updVec3Measure     (const String&);
    StationMeasure&         updStationMeasure  (const String&);
    Station&                updStation         (const String&);
    Direction&              updDirection       (const String&);
    Orientation&            updOrientation     (const String&);
    Frame&                  updFrame           (const String&);

    // Create a new feature on this subsystem with a given name and type, and
    // optionally create a placement for it using the prototype placement supplied.
    RealParameter&    addRealParameter
                        (const String&, const Placement& = Placement());
    Vec3Parameter&    addVec3Parameter
                        (const String&, const Placement& = Placement());
    RealMeasure&      addRealMeasure
                        (const String&, const Placement& = Placement());
    Vec3Measure&      addVec3Measure
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
    // Subsystem is modeled on the prototype Subsystem supplied here.
    Subsystem&        addSubsystemLike(const Subsystem&, const String&);

    // This generic routine will create a Feature of the specific type supplied
    // in the Subsystem argument, which must turn out to be a Feature. If a 
    // Placement is provided, it must be appropriate for that kind of Feature.
    Feature&          addFeatureLike
                        (const Subsystem&,
                         const String&, const Placement& = Placement());
   
    bool              hasParentSubsystem() const;
    int               getIndexInParent()   const; // -1 if no parent
    const Subsystem&  getParentSubsystem() const;

    bool isSameSubsystem(const Subsystem& s) const;

    String getName()              const;
    String getFullName()          const; // "ancestors/parent/name"


    // Subsystems which are children of this subsystem
    int              getNSubsystems()  const;
    const Subsystem& getSubsystem(int) const;
    Subsystem&       updSubsystem(int);

    const Subsystem& operator[](int i) const           {return getSubsystem(i);}
    Subsystem&       operator[](int i)                 {return updSubsystem(i);}
    const Subsystem& operator[](const String& s) const {return getSubsystem(s);}
    Subsystem&       operator[](const String& s)       {return updSubsystem(s);}

    String toString(const String& linePrefix="") const;

    // For internal use only.
    bool                      hasRep() const            {return rep != 0;}
    const class SubsystemRep& getRep() const            {assert(rep); return *rep;}
    class SubsystemRep&       updRep()                  {assert(rep); return *rep;}
    void                      setRep(SubsystemRep* fp)  {assert(!rep); rep=fp;}
    void checkSubsystemConsistency(const Subsystem* expParent,
                                   int              expIndexInParent,
                                   const Subsystem& expRoot) const;
protected:
    class SubsystemRep* rep;
    friend class SubsystemRep;
};

class Feature : public Subsystem {
public:
    Feature(const Feature&);    // external placements are not copied or assigned
    Feature& operator=(const Feature&);
    ~Feature();

    String getFeatureTypeName() const; // "Station", "Frame", etc. (for messages only)

    // Return value of this feature's placement. Requires (a) that there is
    // a placement, and (b) that its value is available due to a prior realize()
    // call.
    const PlacementValue& getValue() const;


    // True if this feature has been placed; its children may still be unplaced.
    bool hasPlacement() const;
    const Placement& getPlacement() const;

    // Place this Feature using the supplied Placement expression
    // as a prototype. We will choose an owner Subsystem for the new
    // Placement, and then add a copy of the prototype to that
    // owner. Then this Feature will refer to that copy as its Placement.
    void place(const Placement&);


    // True if this is the same feature as f or if the feature's placement
    // depends on f's placement.
    bool dependsOn(const Feature& f) const;

    const class FeatureRep& getRep() const
      { return *reinterpret_cast<const FeatureRep*>(rep); }    
    class FeatureRep& updRep()
      { return *reinterpret_cast<FeatureRep*>(rep); } 

    static bool           isInstanceOf(const Subsystem&);
    static const Feature& downcast(const Subsystem&);
    static Feature&       downcast(Subsystem&);
protected:
    Feature() { }
};

// Global operators involving Subsystems and Features.
std::ostream& operator<<(std::ostream& o, const Subsystem&);

// Although these operators appear to act on Features, they actually
// create a Placement referring to the Features and then perform
// the operations on the Placement.

// unary
inline Placement negate    (const Feature& f) {return negate   (Placement(f));}
inline Placement abs       (const Feature& f) {return abs      (Placement(f));}
inline Placement sqrt      (const Feature& f) {return sqrt     (Placement(f));}
inline Placement exp       (const Feature& f) {return exp      (Placement(f));}
inline Placement log       (const Feature& f) {return log      (Placement(f));}
inline Placement sin       (const Feature& f) {return sin      (Placement(f));}
inline Placement cos       (const Feature& f) {return cos      (Placement(f));}
inline Placement asin      (const Feature& f) {return asin     (Placement(f));}
inline Placement acos      (const Feature& f) {return acos     (Placement(f));}
inline Placement length    (const Feature& f) {return length   (Placement(f));}
inline Placement normalize (const Feature& f) {return normalize(Placement(f));}

// binary (feature,feature)
inline Placement add       (const Feature& l, const Feature& r) {return add(Placement(l),Placement(r));} 
inline Placement subtract  (const Feature& l, const Feature& r) {return subtract(Placement(l),Placement(r));} 
inline Placement multiply  (const Feature& l, const Feature& r) {return multiply(Placement(l),Placement(r));} 
inline Placement divide    (const Feature& l, const Feature& r) {return divide(Placement(l),Placement(r));}
inline Placement distance  (const Feature& l, const Feature& r) {return distance(Placement(l),Placement(r));}
inline Placement angle     (const Feature& l, const Feature& r) {return angle(Placement(l),Placement(r));}
inline Placement dot       (const Feature& l, const Feature& r) {return dot(Placement(l),Placement(r));}
inline Placement cross     (const Feature& l, const Feature& r) {return cross(Placement(l),Placement(r));}

// binary (feature,placement)
inline Placement add       (const Feature& l, const Placement& r) {return add(Placement(l),r);} 
inline Placement subtract  (const Feature& l, const Placement& r) {return subtract(Placement(l),r);} 
inline Placement multiply  (const Feature& l, const Placement& r) {return multiply(Placement(l),r);} 
inline Placement divide    (const Feature& l, const Placement& r) {return divide(Placement(l),r);}
inline Placement distance  (const Feature& l, const Placement& r) {return distance(Placement(l),r);}
inline Placement angle     (const Feature& l, const Placement& r) {return angle(Placement(l),r);}
inline Placement dot       (const Feature& l, const Placement& r) {return dot(Placement(l),r);}
inline Placement cross     (const Feature& l, const Placement& r) {return cross(Placement(l),r);}

// binary (placement,feature)
inline Placement add       (const Placement& l, const Feature& r) {return add(l,Placement(r));} 
inline Placement subtract  (const Placement& l, const Feature& r) {return subtract(l,Placement(r));} 
inline Placement multiply  (const Placement& l, const Feature& r) {return multiply(l,Placement(r));} 
inline Placement divide    (const Placement& l, const Feature& r) {return divide(l,Placement(r));}
inline Placement distance  (const Placement& l, const Feature& r) {return distance(l,Placement(r));}
inline Placement angle     (const Placement& l, const Feature& r) {return angle(l,Placement(r));}
inline Placement dot       (const Placement& l, const Feature& r) {return dot(l,Placement(r));}
inline Placement cross     (const Placement& l, const Feature& r) {return cross(l,Placement(r));}

// Operator alternates for some of the above

// unary
inline Placement operator-(const Feature& f) {return negate(f);}

// binary (ff,fp,pf)
inline Placement operator+(const Feature&   l, const Feature&   r) {return add(l,r);}
inline Placement operator+(const Feature&   l, const Placement& r) {return add(l,r);}
inline Placement operator+(const Placement& l, const Feature&   r) {return add(l,r);}
inline Placement operator-(const Feature&   l, const Feature&   r) {return subtract(l,r);} 
inline Placement operator-(const Feature&   l, const Placement& r) {return subtract(l,r);} 
inline Placement operator-(const Placement& l, const Feature&   r) {return subtract(l,r);} 
inline Placement operator*(const Feature&   l, const Feature&   r) {return multiply(l,r);} 
inline Placement operator*(const Feature&   l, const Placement& r) {return multiply(l,r);} 
inline Placement operator*(const Placement& l, const Feature&   r) {return multiply(l,r);} 
inline Placement operator/(const Feature&   l, const Feature&   r) {return divide(l,r);} 
inline Placement operator/(const Feature&   l, const Placement& r) {return divide(l,r);} 
inline Placement operator/(const Placement& l, const Feature&   r) {return divide(l,r);} 

// binary (feature,Real; Real,feature)
inline Placement operator+(const Feature& l, const Real&    r) {return add(l,RealPlacement(r));}
inline Placement operator+(const Real&    l, const Feature& r) {return add(RealPlacement(l),r);} 
inline Placement operator-(const Feature& l, const Real&    r) {return subtract(l,RealPlacement(r));} 
inline Placement operator-(const Real&    l, const Feature& r) {return subtract(RealPlacement(l),r);}  
inline Placement operator*(const Feature& l, const Real&    r) {return multiply(l,RealPlacement(r));} 
inline Placement operator*(const Real&    l, const Feature& r) {return multiply(RealPlacement(l),r);} 
inline Placement operator/(const Feature& l, const Real&    r) {return divide(l,RealPlacement(r));} 
inline Placement operator/(const Real&    l, const Feature& r) {return divide(RealPlacement(l),r);}

// binary (feature,Vec3; Vec3,feature)
inline Placement operator+(const Feature& l, const Vec3&    r) {return add(l,Vec3Placement(r));}
inline Placement operator+(const Vec3&    l, const Feature& r) {return add(Vec3Placement(l),r);} 
inline Placement operator-(const Feature& l, const Vec3&    r) {return subtract(l,Vec3Placement(r));} 
inline Placement operator-(const Vec3&    l, const Feature& r) {return subtract(Vec3Placement(l),r);}  
inline Placement operator*(const Feature& l, const Vec3&    r) {return multiply(l,Vec3Placement(r));} 
inline Placement operator*(const Vec3&    l, const Feature& r) {return multiply(Vec3Placement(l),r);} 
inline Placement operator/(const Feature& l, const Vec3&    r) {return divide(l,Vec3Placement(r));} 
inline Placement operator/(const Vec3&    l, const Feature& r) {return divide(Vec3Placement(l),r);} 

// binary (feature,Mat33; Mat33,feature) (TODO: Mat33Placement)
inline Placement operator+(const Feature& l, const Mat33&    r) {return add(l,Placement(r));}
inline Placement operator+(const Mat33&    l, const Feature& r) {return add(Placement(l),r);} 
inline Placement operator-(const Feature& l, const Mat33&    r) {return subtract(l,Placement(r));} 
inline Placement operator-(const Mat33&    l, const Feature& r) {return subtract(Placement(l),r);}  
inline Placement operator*(const Feature& l, const Mat33&    r) {return multiply(l,Placement(r));} 
inline Placement operator*(const Mat33&    l, const Feature& r) {return multiply(Placement(l),r);} 
inline Placement operator/(const Feature& l, const Mat33&    r) {return divide(l,Placement(r));} 
inline Placement operator/(const Mat33&    l, const Feature& r) {return divide(Placement(l),r);} 

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
    const Real& getValue() const;

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
    const Real& getValue() const;

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
    const Vec3& getValue() const;

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
    const Vec3& getValue() const;

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
    const Vec3& getValue() const;

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
    const Vec3& getValue() const;

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
    const Vec3& getValue() const;

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
    const Vec3& getValue() const;

    static bool                    isInstanceOf(const Subsystem&);
    static const DirectionMeasure& downcast(const Subsystem&);
    static DirectionMeasure&       downcast(Subsystem&);
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

    const DirectionPlacement& getPlacement() const;
    const Vec3& getValue() const;

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
    const Mat33& getValue() const;

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
    Orientation(const String& name, const Mat33& defaultValue);
    Orientation(const Orientation&);
    Orientation& operator=(const Orientation&);
    ~Orientation();

    const OrientationPlacement& getPlacement() const;
    const Mat33& getValue() const;

    const Direction& getAxis(int) const;
    const Direction&   x()        const {return getAxis(0);}
    const Direction&   y()        const {return getAxis(1);}
    const Direction&   z()        const {return getAxis(2);}

    static bool               isInstanceOf(const Subsystem&);
    static const Orientation& downcast(const Subsystem&);
    static Orientation&       downcast(Subsystem&);
protected:
    Orientation() { }
};

class Frame : public Feature {
public:
    explicit Frame(const String& name);
    Frame(const Frame&);
    Frame& operator=(const Frame&);
    ~Frame();

    const FramePlacement& getPlacement() const;
    const Mat34& getValue() const;

    const Station&     getOrigin() const;
    const Orientation& getOrientation() const;
    const Direction&   getAxis(int i) const {return getOrientation().getAxis(i);}
    const Direction&   x()            const {return getOrientation().x();}
    const Direction&   y()            const {return getOrientation().y();}
    const Direction&   z()            const {return getOrientation().z();}

    static bool         isInstanceOf(const Subsystem&);
    static const Frame& downcast(const Subsystem&);
    static Frame&       downcast(Subsystem&);
protected:
    Frame() { }
};

} // namespace simtk

#endif // SIMTK_FEATURE_H_
