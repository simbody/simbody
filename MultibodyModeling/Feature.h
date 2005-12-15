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
#include "Subsystem.h"
#include "Placement.h"
#include "BasicPlacements.h"

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

inline Placement add       (const Feature& l, const Real& r) {return add(l,RealPlacement(r));} 
inline Placement subtract  (const Feature& l, const Real& r) {return subtract(l,RealPlacement(r));} 
inline Placement multiply  (const Feature& l, const Real& r) {return multiply(l,RealPlacement(r));} 
inline Placement divide    (const Feature& l, const Real& r) {return divide(l,RealPlacement(r));}
inline Placement distance  (const Feature& l, const Real& r) {return distance(l,RealPlacement(r));}
inline Placement angle     (const Feature& l, const Real& r) {return angle(l,RealPlacement(r));}
inline Placement dot       (const Feature& l, const Real& r) {return dot(l,RealPlacement(r));}
inline Placement cross     (const Feature& l, const Real& r) {return cross(l,RealPlacement(r));}

inline Placement add       (const Feature& l, const Vec3& r) {return add(l,Vec3Placement(r));} 
inline Placement subtract  (const Feature& l, const Vec3& r) {return subtract(l,Vec3Placement(r));} 
inline Placement multiply  (const Feature& l, const Vec3& r) {return multiply(l,Vec3Placement(r));} 
inline Placement distance  (const Feature& l, const Vec3& r) {return distance(l,Vec3Placement(r));}
inline Placement angle     (const Feature& l, const Vec3& r) {return angle(l,Vec3Placement(r));}
inline Placement dot       (const Feature& l, const Vec3& r) {return dot(l,Vec3Placement(r));}
inline Placement cross     (const Feature& l, const Vec3& r) {return cross(l,Vec3Placement(r));}

// TODO: Mat33Placement
inline Placement add       (const Feature& l, const Mat33& r) {return add(l,Placement(r));} 
inline Placement subtract  (const Feature& l, const Mat33& r) {return subtract(l,Placement(r));} 
inline Placement multiply  (const Feature& l, const Mat33& r) {return multiply(l,Placement(r));} 


// binary (placement,feature)
inline Placement add       (const Placement& l, const Feature& r) {return add(l,Placement(r));} 
inline Placement subtract  (const Placement& l, const Feature& r) {return subtract(l,Placement(r));} 
inline Placement multiply  (const Placement& l, const Feature& r) {return multiply(l,Placement(r));} 
inline Placement divide    (const Placement& l, const Feature& r) {return divide(l,Placement(r));}
inline Placement distance  (const Placement& l, const Feature& r) {return distance(l,Placement(r));}
inline Placement angle     (const Placement& l, const Feature& r) {return angle(l,Placement(r));}
inline Placement dot       (const Placement& l, const Feature& r) {return dot(l,Placement(r));}
inline Placement cross     (const Placement& l, const Feature& r) {return cross(l,Placement(r));}

inline Placement add       (const Real& l, const Feature& r) {return add(RealPlacement(l),r);} 
inline Placement subtract  (const Real& l, const Feature& r) {return subtract(RealPlacement(l),r);} 
inline Placement multiply  (const Real& l, const Feature& r) {return multiply(RealPlacement(l),r);} 
inline Placement divide    (const Real& l, const Feature& r) {return divide(RealPlacement(l),r);}
inline Placement distance  (const Real& l, const Feature& r) {return distance(RealPlacement(l),r);}
inline Placement angle     (const Real& l, const Feature& r) {return angle(RealPlacement(l),r);}
inline Placement dot       (const Real& l, const Feature& r) {return dot(RealPlacement(l),r);}
inline Placement cross     (const Real& l, const Feature& r) {return cross(RealPlacement(l),r);}

inline Placement add       (const Vec3& l, const Feature& r) {return add(Vec3Placement(l),r);} 
inline Placement subtract  (const Vec3& l, const Feature& r) {return subtract(Vec3Placement(l),r);} 
inline Placement multiply  (const Vec3& l, const Feature& r) {return multiply(Vec3Placement(l),r);} 
inline Placement divide    (const Vec3& l, const Feature& r) {return divide(Vec3Placement(l),r);}
inline Placement distance  (const Vec3& l, const Feature& r) {return distance(Vec3Placement(l),r);}
inline Placement angle     (const Vec3& l, const Feature& r) {return angle(Vec3Placement(l),r);}
inline Placement dot       (const Vec3& l, const Feature& r) {return dot(Vec3Placement(l),r);}
inline Placement cross     (const Vec3& l, const Feature& r) {return cross(Vec3Placement(l),r);}

// TODO: Mat33Placement
inline Placement add       (const Mat33& l, const Feature& r) {return add(Placement(l),r);} 
inline Placement subtract  (const Mat33& l, const Feature& r) {return subtract(Placement(l),r);} 
inline Placement multiply  (const Mat33& l, const Feature& r) {return multiply(Placement(l),r);} 
inline Placement divide    (const Mat33& l, const Feature& r) {return divide(Placement(l),r);}

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
inline Placement operator+(const Feature& l, const Real&    r) {return add(l,r);}
inline Placement operator+(const Real&    l, const Feature& r) {return add(l,r);} 
inline Placement operator-(const Feature& l, const Real&    r) {return subtract(l,r);} 
inline Placement operator-(const Real&    l, const Feature& r) {return subtract(l,r);}  
inline Placement operator*(const Feature& l, const Real&    r) {return multiply(l,r);} 
inline Placement operator*(const Real&    l, const Feature& r) {return multiply(l,r);} 
inline Placement operator/(const Feature& l, const Real&    r) {return divide(l,r);} 
inline Placement operator/(const Real&    l, const Feature& r) {return divide(l,r);}

// binary (feature,Vec3; Vec3,feature)
inline Placement operator+(const Feature& l, const Vec3&    r) {return add(l,r);}
inline Placement operator+(const Vec3&    l, const Feature& r) {return add(l,r);} 
inline Placement operator-(const Feature& l, const Vec3&    r) {return subtract(l,r);} 
inline Placement operator-(const Vec3&    l, const Feature& r) {return subtract(l,r);}  
inline Placement operator*(const Feature& l, const Vec3&    r) {return multiply(l,r);} 
inline Placement operator*(const Vec3&    l, const Feature& r) {return multiply(l,r);}  
inline Placement operator/(const Vec3&    l, const Feature& r) {return divide(l,r);} 

// binary (feature,Mat33; Mat33,feature) (TODO: Mat33Placement)
inline Placement operator+(const Feature& l, const Mat33&    r) {return add(l,r);}
inline Placement operator+(const Mat33&    l, const Feature& r) {return add(l,r);} 
inline Placement operator-(const Feature& l, const Mat33&    r) {return subtract(l,r);} 
inline Placement operator-(const Mat33&    l, const Feature& r) {return subtract(l,r);}  
inline Placement operator*(const Feature& l, const Mat33&    r) {return multiply(l,r);} 
inline Placement operator*(const Mat33&    l, const Feature& r) {return multiply(l,r);} 
inline Placement operator/(const Mat33&    l, const Feature& r) {return divide(l,r);} 

} // namespace simtk

#endif // SIMTK_FEATURE_H_
