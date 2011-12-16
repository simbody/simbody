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
#include <cmath>
#include <algorithm>

namespace SimTK {

//==============================================================================
//                                    GEO
//==============================================================================
/** The Geo class collects geometric primitives intended to deal with raw, 
fixed-size geometric shapes occupying minimal memory and providing
maximum performance through small inline methods and larger high performance
algorithms. Subclasses collect algorithms relevant to particular shapes. 
There are no virtual methods or class hierarchies here; each subclass is a
"POD" (plain old data) class. The general idea is to make it so that these
common methods are implemented in only one place in Simbody.

The Geo class itself is dataless and provides only static methods. **/
class SimTK_SIMMATH_EXPORT Geo {
public:
template <class P> class Point_;
template <class P> class Sphere_;
template <class P> class LineSeg_;
template <class P> class Line_;
template <class P> class Plane_;
template <class P> class Circle_;
template <class P> class Box_;
template <class P> class Triangle_;
template <class P> class BicubicPatch_;

typedef Point_<Real>        Point;
typedef Sphere_<Real>       Sphere;
typedef LineSeg_<Real>      LineSeg;
typedef Line_<Real>         Line;
typedef Plane_<Real>        Plane;
typedef Circle_<Real>       Circle;
typedef Box_<Real>          Box;
typedef Triangle_<Real>     Triangle;
typedef BicubicPatch_<Real> BicubicPatch;

/** Return the default tolerance to use for degeneracy tests and other tests
for "too small" or "near enough" that arise in dealing with geometry primitives.
The value depends on the precision being used; we use the SimTK constant
SignificantReal which is eps^(7/8) where eps is the resolution of P. That
makes this tolerance around 2e-14 in double precision and 9e-7 in float. **/
template <class P> static P getDefaultTol() 
{   return NTraits<P>::getSignificant(); }
template <class P> static P getEps() 
{   return NTraits<P>::getEps(); }
template <class P> static P getNaN() 
{   return NTraits<P>::getNaN(); }
template <class P> static P getInfinity() 
{   return NTraits<P>::getInfinity(); }

/** Stretch a dimension by a given tolerance amount. The result is the
given \a length increased by at least an absolute amount \a tol, or by a
relative amount length*tol if length > 1. Cost is 3 flops. **/
template <class P> static P stretchBy(P length, P tol)
{   assert(tol >= getEps<P>()); 
    return length + std::max(length*tol, tol); }

/** Stretch a dimension using the default tolerance for this precision as
the tolerance in stretchBy(). Cost is 3 flops. **/
template <class P> static P stretch(P length)
{   return stretchBy(length, getDefaultTol<P>()); }

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

/** Calculate the distance between this point and another one whose location is
expressed in the same frame (expensive). Cost is about 30 flops. **/
RealP calcDistance(const Vec3P& p2) const 
{   return calcDistance(p, p2); }
/** Find the square of the distance between this point and another one whose
location is expressed in the same frame (cheap). Cost is 8 flops. **/
RealP findDistanceSqr(const Vec3P& p2) const 
{   return findDistanceSqr(p, p2); }

/**@name                 Point-related utilities
These static methods work with points or collections of points. **/
/**@{**/
/** Calculate the distance between two points (expensive). Cost is about 
30 flops. **/
static RealP calcDistance(const Vec3P& p1, const Vec3P& p2)
{   return std::sqrt(findDistanceSqr(p1,p2)); }
/** Find the square of the distance between two points (cheap). Cost is 
8 flops. **/
static RealP findDistanceSqr(const Vec3P& p1, const Vec3P& p2)
{   return (p2-p1).normSqr(); }
/** Find the point midway between two points. Cost is 4 flops. **/
static Vec3P findMidpoint(const Vec3P& p1, const Vec3P& p2)
{   return (p1+p2)/2; }
/**@}**/

private:
Vec3P   p;
};



//==============================================================================
//                              GEO SPHERE
//==============================================================================
/** A geometric primitive representing a sphere by its radius and center 
point, and a collection of sphere-related utility methods. **/
template <class P>
class Geo::Sphere_ {
typedef P           RealP;
typedef Vec<3,P>    Vec3P;
typedef Vec<4,P>    Vec4P;
public:
/** Construct an uninitialized Sphere object; the center point and radius 
will be garbage. **/
Sphere_() {}
/** Construct a sphere from its center location and radius. **/
Sphere_(const Vec3P& location, RealP radius)
:   cr(location[0], location[1], location[2], radius) {}
/** Change the radius of this sphere. **/
Sphere_& setRadius(RealP radius) 
{   cr[3]=radius; return *this; }
/** Change the location of this sphere. **/
Sphere_& setLocation(const Vec3P& location) 
{   Vec3P::updAs(&cr[0])=location; return *this; }

/** Modify this sphere to scale its radius by a fractional amount f,
that is we set radius to f*radius.
@return A reference to this now-resized sphere. **/
Sphere_& scaleBy(RealP f)
{   setRadius(f*getRadius()); return *this; }

/** Stretch this sphere by a small amount to ensure that there will be no 
roundoff problems if this is used as a bounding sphere. The amount to stretch
depends on the default tolerance for this precision, the radius, and the
position of the sphere in space. A very large sphere, or a sphere that is very
far from the origin, must be stretched more than a small one at the origin. 
@see Geo class for tolerance information. **/
Sphere_& stretchBoundary() {
    const RealP tol = Geo::getDefaultTol<P>();
    const RealP maxdim = max(getCenter().abs());
    const RealP scale = std::max(maxdim, getRadius());
    updRadius() += std::max(scale*Geo::getEps<P>(), tol);
    return *this; 
}
/** Return true if a given point is strictly outside this sphere. Just touching
the sphere does not qualify. **/
bool isPointOutside(const Vec3P& p) const {
    const RealP r2 = Geo::Point_<P>::findDistanceSqr(p, getCenter());
    return r2 > square(getRadius());
}
const Vec3P& getCenter() const {return Vec3P::getAs(&cr[0]);}
Vec3P& updCenter() {return Vec3P::updAs(&cr[0]);}
RealP getRadius() const {return cr[3];}
RealP& updRadius() {return cr[3];}

/**@name                 Sphere-related utilities
These static methods work with spheres or collections of spheres. 
  - Minimum sphere methods calculate the smallest sphere around a given set of
    points such that no point is outside the sphere, although some may be on 
    its surface. How many and specifically which points were actually used to 
    define the sphere is returned; there will never be more than 4. Roundoff
    errors are expected so defining points may not be exactly on the surface
    and some points may be slightly outside.
  - Bounding sphere methods make use of the minimum sphere calculations but
    address roundoff by stretching the sphere enough to guarantee that all
    points are strictly inside the sphere. That avoids trouble later when these
    are used to look for possible intersections -- you could miss one using
    a minimum sphere but you won't if you use a bounding sphere.
    
**/
/**@{**/

/** Create a tiny bounding sphere around a single point. The center is the
point and the radius is tiny but non-zero. **/
static Sphere_ calcBoundingSphere(const Vec3P& p)
{   return Sphere_(p, 0).stretchBoundary(); }

/** Create a minimal bounding sphere around two points. The center is the 
midpoint, and the radius is half the distance between the points, plus a small 
amount to avoid roundoff problems later. Cost is about 35 flops. **/
static Sphere_ calcBoundingSphere(const Vec3P& p0, const Vec3P& p1)
{   const RealP r = Point_<P>::calcDistance(p0,p1)/2;
    return Sphere_(Point_<P>::findMidpoint(p0,p1), r).stretchBoundary(); }

/** Create a minimal bounding sphere around three points. **/
SimTK_SIMMATH_EXPORT static Sphere_ calcBoundingSphere
   (const Vec3P& p0, const Vec3P& p1, const Vec3P& p2) {
    Array_<int> which;
    Sphere_ minSphere = calcMinimumSphere(p0,p1,p2,which);
    return minSphere.stretchBoundary();
}

/** Create a minimal bounding sphere around four points. **/
SimTK_SIMMATH_EXPORT static Sphere_ calcBoundingSphere
   (const Vec3P& p0, const Vec3P& p1, const Vec3P& p2, const Vec3P& p3) {
    Array_<int> which;
    Sphere_ minSphere = calcMinimumSphere(p0,p1,p2,p3,which);
    return minSphere.stretchBoundary();
}

/** Create a minimal bounding sphere around a collection of n points. 
This has expected O(n) performance and yields a perfect bounding sphere. **/
static Sphere_ calcBoundingSphere(const Array_<Vec3P>& points) {
    Array_<int> which; 
    Sphere_ minSphere = calcMinimumSphere(points, which);
    return minSphere.stretchBoundary();
}

/** Create a minimum sphere around a single point. The center is the
point and the radius is zero. There is always 1 support point. **/
static Sphere_ calcMinimumSphere(const Vec3P& p0, Array_<int>& which) 
{   which.clear(); which.push_back(0); return Sphere_(p0,0); }

/** Create a minimum sphere around two points. The center is the 
midpoint, and the radius is half the distance between the points. There will
be two support points for the circle unless the given points are within
twice machine epsilon of each other. In that case, we treat these as a single 
point and report that only 1 point was used to define the sphere. Points
far from the origin will produce a larger sphere because of roundoff.
Cost is about 40 flops. **/
static Sphere_ calcMinimumSphere(const Vec3P& p0, const Vec3P& p1,
                                 Array_<int>& which);

/** Create a minimum sphere around three points. There can be 1, 2, or 3
support points returned. **/
SimTK_SIMMATH_EXPORT static Sphere_ 
calcMinimumSphere(const Vec3P& p0, const Vec3P& p1, const Vec3P& p2, 
                  Array_<int>& which);

/** Create a minimum sphere around four points. There can be 1, 2, 3, or 4
support points returned.  **/
SimTK_SIMMATH_EXPORT static Sphere_ 
calcMinimumSphere(const Vec3P& p0, const Vec3P& p1, const Vec3P& p2, 
                  const Vec3P& p3, Array_<int>& which);

/** Create an optimal minimum sphere around a collection of n points. This has 
expected O(n) performance and yields a perfect minimum sphere. There can be
1, 2, 3, or 4 support points used to define the sphere and \a which reports
which of the input points were used. **/
SimTK_SIMMATH_EXPORT static Sphere_ 
calcMinimumSphere(const Array_<Vec3P>& points, Array_<int>& which);

/**@}**/

private:
// Store together to make sure the compiler doesn't introduce any padding.
// cr[0..2] is the center point, cr[3] is the radius
Vec4P   cr;    
};


//==============================================================================
//                               GEO LINESEG
//==============================================================================
/** A 3d line segment primitive represented by its end points
in an unspecified frame, and a collection of line segment-related utility 
methods. We support a one-parameter t representation of the line segment, where
t=0 gives the first end point and t=1 gives the second end point. **/
template <class P>
class SimTK_SIMMATH_EXPORT Geo::LineSeg_ {
typedef P               RealP;
typedef Vec<3,RealP>    Vec3P;

public:
/** Construct an uninitialized LineSeg object; the end points will 
be garbage (NaN in Debug builds). **/
LineSeg_() {}

/** Construct a LineSeg with the given end points. When an orientation is
needed we define the line segment to go from \a e0 to \a e1. **/
LineSeg_(const Vec3P& e0, const Vec3P& e1) 
{   setEndpoints(e0,e1); }

/** Change the end points of this line segment. **/
LineSeg_& setEndpoints(const Vec3P& e0, const Vec3P& e1) 
{   e[0]=e0; e[1]=e1; return *this; }

/** Change one end point of this line segment. **/
LineSeg_& setEndpoint(int which, const Vec3P& p)
{   assert(which==0 || which==1); e[which] = p; return *this; }

/** Determine whether this line segment is degenerate to a given tolerance,
meaning that its length is less than or equal to the tolerance. Use a tolerance
of zero if you want to check only for exact degeneracy. By default the
tolerance is SignificantReal for the precision in use, about 2e-14 in double, 
9e-7 in float. Cost is 10 flops. **/
bool isDegenerate(RealP tol = Geo::getDefaultTol<P>()) const
{   assert(tol >= 0); return calcLengthSqr() <= square(tol); }

/** Get the location of an end point.\ Order is the same as construction. You
can use operator[] instead for a more compact notation. **/
const Vec3P& getEndpoint(int which) const 
{   assert(which==0 || which==1); return e[which]; }

/** Get a writable reference to the location of an end point.\ Order is the 
same as construction. You can use operator[] instead for a more compact 
notation. **/
Vec3P& updEndpoint(int which) 
{   assert(which==0 || which==1); return e[which]; }

/** Access an end point by indexing the line segment. **/
const Vec3P& operator[](int which) const {return getEndpoint(which);}
/** Get writable access to an end point by indexing the line segment. **/
Vec3P& operator[](int which) {return updEndpoint(which);}

/** Calculate the length of this line segment (expensive). Cost is about
40 flops. **/
RealP calcLength() const {return (e[1]-e[0]).norm(); }

/** Calculate the square of the length of this line segment (cheap). Cost 
is 8 flops. **/
RealP calcLengthSqr() const {return (e[1]-e[0]).normSqr(); }

/** Return a point along the line segment given by its coordinate t, where
t==0 gives end point 0, t==1 gives end point 1 and other values interpolate
or extrapolate linearly. Note that values outside the range [0,1] are not
on the line segment. Cost is 10 flops. **/
Vec3P findPoint(RealP t) const
{   return t*e[0] + (1-t)*e[1]; }

/** Return the center point of the line segment; this is the same as
findPoint(1/2) but faster. Cost is 4 flops. **/
Vec3P findMidpoint() const
{   return (e[0]+e[1]) / RealP(2); }

/** Calculate a minimal bounding sphere for this line segment. This is the
sphere whose center is the segment midpoint and whose radius is half the 
segment length, plus a little slop. Cost is about 45 flops. **/
Sphere_<P> calcBoundingSphere() const 
{   return Geo::Sphere_<P>::calcBoundingSphere(e[0],e[1]); }

/** Find the distance between this line segment and a point expressed in the 
same frame. Cost is XXX flops. **/
RealP findDistanceToPoint(const Vec3P& p2) const
{SimTK_ASSERT_ALWAYS(!"implemented", 
"Geo::LineSeg_::findDistanceToPoint(): Not implemented yet.");
return Geo::getNaN<P>();}

/** Find the square of the distance between this line segment and a point 
expressed in the same frame. Cost is XXX flops. **/
RealP findDistanceToPointSqr(const Vec3P& p2) const
{SimTK_ASSERT_ALWAYS(!"implemented", 
"Geo::LineSeg_::findDistanceToPointSqr(): Not implemented yet.");
return Geo::getNaN<P>();}

/**@name                 Line segment-related utilities
These static methods work with points or collections of points. **/
/**@{**/

/**@}**/

private:
Vec3P   e[2];
};

//==============================================================================
//                              GEO TRIANGLE
//==============================================================================
/** A geometric primitive representing a triangle by its vertices as points
in some unspecified frame, and a collection of triangle-related utility 
methods. We support a u-v parameterization for the triangle. **/
template <class P>
class SimTK_SIMMATH_EXPORT Geo::Triangle_ {
typedef P               RealP;
typedef Vec<2,P>        Vec2P;
typedef Vec<3,P>        Vec3P;
typedef UnitVec<P,1>    UnitVec3P;
public:
/** Construct an uninitialized Triangle object; the vertices will 
be garbage. **/
Triangle_() {}
/** Construct a triangle from its vertices. **/
Triangle_(const Vec3P& v0, const Vec3P& v1, const Vec3P& v2)
{  setVertices(v0,v1,v2); }
/** Construct a triangle from vertices stored in array which is presumed
to contain at least three points. **/
explicit Triangle_(const Vec3P* vertices)
{   setVertices(vertices); }
/** Construct a triangle from an indirect list of vertices stored in array 
which is presumed to contain at least three pointers to points. **/
explicit Triangle_(const Vec3P** vertexPointers)
{   setVertices(*vertexPointers[0], *vertexPointers[1], *vertexPointers[2]); }

/** Check whether this triangle is degenerate to within some tolerance. By
default the tolerance is the SignificantReal value for the precision being
used. We define the triangle to be degenerate if any of the line segments
of its boundary are degenerate. Use a tolerance of zero if you only want to 
check for exact degeneracy. Cost is 30 flops. **/
bool isDegenerate(RealP tol = Geo::getDefaultTol<P>) const {
    return LineSeg_<P>(v[0],v[1]).isDegenerate(tol)
        || LineSeg_<P>(v[0],v[2]).isDegenerate(tol)
        || LineSeg_<P>(v[1],v[2]).isDegenerate(tol);
}



/** Change one vertex of this triangle. **/
Triangle_& setVertex(int i, const Vec3P& p)
{   assert(0<=i && i<3); v[i] = p; return *this; }
/** Change the vertices of this triangle. **/
Triangle_& setVertices(const Vec3P& v0, const Vec3P& v1, const Vec3P& v2)
{   v[0]=v0; v[1]=v1; v[2]=v2; return *this; }
/** Change the vertices of this triangle, taking the new ones from an array 
which is presumed to contain at least three points. **/
Triangle_& setVertices(const Vec3P* vertices)
{   v[0]=vertices[0]; v[1]=vertices[1]; v[2]=vertices[2]; return *this; }

/** Access a vertex by index.\ Order is the same as construction. You can use
operator[] instead for a more compact notation. **/
const Vec3P& getVertex(int i) const
{   SimTK_INDEXCHECK(i,3,"Geo::Triangle_::getVertex()"); 
    return v[i]; }
/** Get writable access to a vertex by index.\ Order is the same as 
construction. You can use operator[] instead for a more compact notation.**/
Vec3P& updVertex(int i)
{   SimTK_INDEXCHECK(i,3,"Geo::Triangle_::updVertex()");  
    return v[i]; }

/** Access a vertex by indexing the triangle. **/
const Vec3P& operator[](int i) const {return getVertex(i);}
/** Get writable access to a vertex by indexing the triangle. **/
Vec3P& operator[](int i) {return updVertex(i);}

/** Return a LineSeg_ containing an edge of this triangle, with edges numbered
in a counterclockwise direction so that edge0 is v0v1, edge1 is v1v2, and 
edge2 is v2v0. **/
LineSeg_<P> getEdge(int i) const
{   SimTK_INDEXCHECK(i,3,"Geo::Triangle_::getEdge()");
    return LineSeg_<P>(v[i],v[(i+1)%3]); }

/** Return a point on the triangle's face given by its (u,v) coordinates. 
Cost is 13 flops. **/
Vec3P findPoint(const Vec2P& uv) const
{   return uv[0]*v[0] + uv[1]*v[1] + (1-uv[0]-uv[1])*v[2]; }

/** Return the centroid point on the triangle's face; this is the same as
findPoint(1/3,1/3) but faster. Cost is 7 flops. **/
Vec3P findCentroid() const
{   return (v[0]+v[1]+v[2]) / RealP(3); }

/** Calculate the unit normal to the triangle face taking the vertices in
counterclockwise order. Cost is about 50 flops. **/
UnitVec3P calcUnitNormal() const 
{   return UnitVec3P((v[1]-v[0]) % (v[2]-v[0])); }

/** Calculate the smallest bounding sphere enclosing the three vertices of 
this triangle. We guarantee that no vertex is outside the sphere, but for
numerical reasons it is possible for the sphere to be a little too big. **/
Sphere_<P> calcBoundingSphere() const
{   return Geo::Sphere_<P>::calcBoundingSphere(v[0],v[1],v[2]); }

/** Return the area of this triangle. Cost is about 50 flops. **/
RealP calcArea() const 
{   return ((v[1]-v[0]) % (v[2]-v[0])).norm() / 2; }

/** Given a location in space, find the point of this triangular face that
is closest to that location. If the answer is not unique then one of the 
equidistant points is returned. **/
Vec3P findNearestPoint(const Vec3P& position, Vec2P& uv) const
{SimTK_ASSERT_ALWAYS(!"implemented", 
"Geo::Triangle_::findNearestPoint(): Not implemented yet.");
return Vec3P();}

/** Determine whether a given ray intersects this triangle. TODO: is a 
perfect hit on the boundary an intersection? **/
bool intersectsRay(const Vec3P& origin, const UnitVec3P& direction, 
                   RealP& distance, Vec2P& uv) const
{SimTK_ASSERT_ALWAYS(!"implemented", 
"Geo::Triangle_::intersectsRay(): Not implemented yet."); return false;}

/** Determine yes/no whether this triangle overlaps another one. Note that
exactly touching is not overlapping. **/
bool overlapsTriangle(const Triangle_<P>& other) const;
/** Determine whether this triangle intersects another one, and if so then
if they are not coplanar return the line segment representing their 
intersection. If the triangles are coplanar and intersect we'll return true
and set \a isCoplanar true, but not return a line segment. Note that the 
triangles may meet at a point so the line segment may be degenerate in any
case. **/
bool intersectsTriangle(const Triangle_<P>& other, LineSeg_<P>& seg,
                        bool& isCoplanar) const;

/**@name                 Triangle-related utilities
These static methods work with triangles or collections of triangles. **/
/**@{**/

/**@}**/

private:
// Vertices in the order they were supplied in the constructor.
Vec3P   v[3];   
};

} // namespace SimTK

#endif // SimTK_SIMMATH_CONTACT_GEOMETRY_H_
