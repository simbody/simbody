#ifndef SimTK_SIMMATH_GEO_TRIANGLE_H_
#define SimTK_SIMMATH_GEO_TRIANGLE_H_

/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2011-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

/** @file
Defines primitive operations on triangles. **/

#include "SimTKcommon.h"
#include "simmath/internal/common.h"
#include "simmath/internal/Geo.h"
#include "simmath/internal/Geo_Point.h"
#include "simmath/internal/Geo_LineSeg.h"

#include <cassert>
#include <cmath>
#include <algorithm>

namespace SimTK {

//==============================================================================
//                              GEO TRIANGLE
//==============================================================================
/** A geometric primitive representing a triangle by its vertices as points
in some unspecified frame, and a collection of triangle-related utility
methods. We support a u-v barycentric parameterization for the triangle. **/
template <class P>
class Geo::Triangle_ {
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
Cost is 17 flops. **/
Vec3P findPoint(const Vec2P& uv) const
{   return uv[0]*v[0] + uv[1]*v[1] + (1-uv[0]-uv[1])*v[2]; }

/** Return the centroid point on the triangle's face; this is the same as
findPoint(1/3,1/3) but faster. Cost is 9 flops. **/
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
{   return Geo::Point_<P>::calcBoundingSphere(v[0],v[1],v[2]); }

/** Return the area of this triangle. Cost is about 40 flops. **/
RealP calcArea() const
{   return ((v[1]-v[0]) % (v[2]-v[0])).norm() / 2; }

/** Return the square of the area of this triangle. Cost is 23 flops. **/
RealP calcAreaSqr() const
{   return ((v[1]-v[0]) % (v[2]-v[0])).normSqr() / 4; }

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
SimTK_SIMMATH_EXPORT bool overlapsTriangle(const Triangle_<P>& other) const;
/** Determine whether this triangle intersects another one, and if so then
if they are not coplanar return the line segment representing their
intersection. If the triangles are coplanar and intersect we'll return true
and set \a isCoplanar true, but not return a line segment. Note that the
triangles may meet at a point so the line segment may be degenerate in any
case. **/
SimTK_SIMMATH_EXPORT bool intersectsTriangle(const Triangle_<P>& other, LineSeg_<P>& seg,
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

#endif // SimTK_SIMMATH_GEO_TRIANGLE_H_
