#ifndef SimTK_SIMMATH_GEO_H_
#define SimTK_SIMMATH_GEO_H_

/* -------------------------------------------------------------------------- *
 *                        SimTK Simbody: SimTKmath                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2011 Stanford University and the Authors.           *
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

/** @file
Defines geometric primitive shapes and algorthms. **/

#include "SimTKcommon.h"
#include "simmath/internal/common.h"

#include <cassert>

namespace SimTK {

//==============================================================================
//                                    GEO
//==============================================================================
/** The Geo class collects geometric primitives intended to deal with raw, 
fixed-size geometric shapes occupying minimal memory and providing
maximum performance through small inline methods and larger high performance
algorithms. Subclasses collect algorithms relevant to particular shapes. 
There are no virtual methods or class hierarchies here; each subclass is a
"POD" (plain old data) class.

The Geo class itself is dataless and provides only static methods. **/
class SimTK_SIMMATH_EXPORT Geo {
public:
template <class P> class Point_;
template <class P> class Line_;
template <class P> class Plane_;
template <class P> class Circle_;
template <class P> class Box_;
template <class P> class Sphere_;
template <class P> class Triangle_;
template <class P> class BicubicPatch_;

typedef Point_<Real>        Point;
typedef Line_<Real>         Line;
typedef Plane_<Real>        Plane;
typedef Circle_<Real>       Circle;
typedef Box_<Real>          Box;
typedef Sphere_<Real>       Sphere;
typedef Triangle_<Real>     Triangle;
typedef BicubicPatch_<Real> BicubicPatch;

};



//==============================================================================
//                               GEO POINT
//==============================================================================
/** A 3d point primitive represented by a Vec3 from the origin of an unspecified 
frame, and a collection of point-related utility methods. **/
template <class P>
class SimTK_SIMMATH_EXPORT Geo::Point_ {
typedef P               RealP;
typedef Vec<3,RealP>    Vec3P;

public:
/** Construct an uninitialized Point object; the location will be garbage. **/
Point_() {}
/** Construct a Point with the given location.\ Also serves as implicit 
conversion from Vec3 to Geo::Point. **/
Point_(const Vec3P& location) : p(location) {} 

/** Change the location of this point. **/
Point_& setLocation(const Vec3P& location) {p=location; return *this;}

/** Get the location of this Point. **/
const Vec3P& getLocation() const {return p;}

/** Find the distance between this point and another one whose location is
expressed in the same frame. Cost is about 35 flops. **/
RealP findDistance(const Point_& p2) const 
{   return findDistance(*this, p2); }
/** Find the square of the distance between this point and another one whose
location is expressed in the same frame. Cost is 5 flops. **/
RealP findDistanceSqr(const Point_& p2) const 
{   return findDistanceSqr(*this, p2); }

/**@name                 Point-related utilities
These static methods work with points or collections of points. **/
/**@{**/
static RealP findDistance(const Point_<P>& p1, const Point_<P>& p2)
{   return (p2.p-p1.p).norm(); }
static RealP findDistanceSqr(const Point_<P>& p1, const Point_<P>& p2)
{   return (p2.p-p1.p).normSqr(); }
/**@}**/

private:
Vec3P   p;
};



//==============================================================================
//                              GEO SPHERE
//==============================================================================
/** A geometric primitive representing a sphere by its radius and center 
point, and a collection of sphere-related utility methods. **/
template <class RealP>
class SimTK_SIMMATH_EXPORT Geo::Sphere_ {
typedef Vec<3,RealP>        Vec3P;
typedef Vec<4,RealP>        Vec4P;
public:
/** Construct an uninitialized Point object; the center point and radius 
will be garbage. **/
Sphere_() {}
/** Construct a sphere from its radius and center location. **/
Sphere_(RealP radius, const Vec3P& location)
:   rc(radius, location[0], location[1], location[2]) {}
/** Change the radius of this sphere. **/
Sphere_& setRadius(RealP radius) 
{   rc[0]=radius; return *this; }
/** Change the location of this sphere. **/
Sphere_& setLocation(const Vec3P& location) 
{   Vec3P::updAs(&rc[1])=location; return *this; }

const Vec3P& getCenter() const {return Vec3P::getAs(&rc[1]);}
RealP getRadius() const {return rc[0];}

/**@name                 Sphere-related utilities
These static methods work with spheres or collections of spheres. **/
/**@{**/

/**@}**/

private:
// Store together to make sure the compiler doesn't introduce any padding.
// rc[0] is the radius, rc[1..3] is the center point
Vec4P   rc;    
};

} // namespace SimTK

#endif // SimTK_SIMMATH_CONTACT_GEOMETRY_H_
