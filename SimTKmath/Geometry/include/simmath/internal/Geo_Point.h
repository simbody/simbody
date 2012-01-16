#ifndef SimTK_SIMMATH_GEO_POINT_H_
#define SimTK_SIMMATH_GEO_POINT_H_

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
Defines primitive computations involving points. **/

#include "SimTKcommon.h"
#include "simmath/internal/common.h"
#include "simmath/internal/Geo.h"

#include <cassert>
#include <cmath>
#include <algorithm>

namespace SimTK {

//==============================================================================
//                               GEO POINT
//==============================================================================
/** A 3d point primitive represented by a Vec3 from the origin of an unspecified 
frame, and a collection of point-related utility methods. **/
template <class P>
class Geo::Point_ {
typedef P               RealP;
typedef Vec<3,P>        Vec3P;
typedef Mat<3,3,P>      Mat33P;
typedef SymMat<3,P>     SymMat33P;
typedef UnitVec<P,1>    UnitVec3P;
typedef Rotation_<P>    RotationP;
typedef Transform_<P>   TransformP;

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

/**@name            Miscellaneous point-related utilities
These static methods work with points or collections of points. Collections
of points are represented either as an Array of point locations or as
an indirect Array of pointers to point locations, which can save a lot of
copying for large point sets. **/
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

/** Given a set of points, find the one that is the furthest in a given
direction, and return its index and location along that direction. There must 
be at least one point in the set. **/
SimTK_SIMMATH_EXPORT static void
findSupportPoint(const Array_<Vec3P>& points, const UnitVec3P& direction,
                 int& most, RealP& mostCoord);

/** Alternate signature taking an array of pointers to points rather than the
points themselves. **/
SimTK_SIMMATH_EXPORT static void
findSupportPointIndirect(const Array_<const Vec3P*>& points, 
                         const UnitVec3P& direction,
                         int& most, RealP& mostCoord);

/** Given a set of points, find the two points that are the most extreme along
a given direction (not necessarily distinct), and return their indices and
locations along the given direction. There must be at least one point
in the set. **/
SimTK_SIMMATH_EXPORT static void
findExtremePoints(const Array_<Vec3P>& points, const UnitVec3P& direction,
                  int& least, int& most,
                  RealP& leastCoord, RealP& mostCoord);

/** Alternate signature taking an array of pointers to points rather than the
points themselves. **/
SimTK_SIMMATH_EXPORT static void
findExtremePointsIndirect(const Array_<const Vec3P*>& points, 
                          const UnitVec3P& direction,
                          int& least, int& most,
                          RealP& leastCoord, RealP& mostCoord);

/** Given a set of points, calculate the centroid (average location) of those 
points. Cost is about 3*n+10 flops for n points. **/
SimTK_SIMMATH_EXPORT static Vec3P
calcCentroid(const Array_<Vec3P>& points_F);

/** Alternate signature taking an array of pointers to points rather than the
points themselves. **/
SimTK_SIMMATH_EXPORT static Vec3P
calcCentroidIndirect(const Array_<const Vec3P*>& points_F);

/** Given a set of points, calculate the centroid (average location) and
covariance matrix of those points. **/
SimTK_SIMMATH_EXPORT static void
calcCovariance(const Array_<Vec3P>& points_F,
               Vec3P& centroid, SymMat33P& covariance);

/** Alternate signature taking an array of pointers to points rather than the
points themselves. **/
SimTK_SIMMATH_EXPORT static void
calcCovarianceIndirect(const Array_<const Vec3P*>& points_F,
                       Vec3P& centroid, SymMat33P& covariance);

/** Given a set of points in an unspecified frame F, find the principal
component directions describing the distribution of the points in space. The
result is a frame P with origin at the centroid, x axis along the direction
of maximum dispersion, y axis along the direction of minimum dispersion, and
z=x X y. Note that clustering of points affects the directions. **/
SimTK_SIMMATH_EXPORT static void
calcPrincipalComponents(const Array_<Vec3P>& points_F,
                        TransformP&          X_FP);

/** Alternate signature taking an array of pointers to points rather than the
points themselves. **/
SimTK_SIMMATH_EXPORT static void
calcPrincipalComponentsIndirect(const Array_<const Vec3P*>& points_F,
                                TransformP&                 X_FP);

/**@}**/

/**@name              Axis-aligned bounding box creation
These static methods create a minimal axis-aligned box that includes all
of a set of given points. **/
/**@{**/

/** Given a set of points, find the six points that are the most extreme along
the axial directions (not necessarily distinct points). Return the indices of
the extreme points and the locations of the box corners. Note that the corners
do not necessarily correspond to any points in the set. There must be at least 
one point in the set. **/
SimTK_SIMMATH_EXPORT static void
findAxisAlignedExtremePoints(const Array_<Vec3P>& points,
                             int least[3], int most[3],
                             Vec3P& low, Vec3P& high);

/** Alternate signature taking an array of pointers to points rather than the
points themselves. **/
SimTK_SIMMATH_EXPORT static void
findAxisAlignedExtremePointsIndirect(const Array_<const Vec3P*>& points,
                                     int least[3], int most[3],
                                     Vec3P& low, Vec3P& high);

/** Calculate the smallest axis-aligned bounding box including all n given
points. Cost is O(n). **/
SimTK_SIMMATH_EXPORT static Geo::AlignedBox_<P> 
calcAxisAlignedBoundingBox(const Array_<Vec3P>& points,
                           Array_<int>&         support);

/** Alternate signature doesn't return support points. **/
static Geo::AlignedBox_<P> 
calcAxisAlignedBoundingBox(const Array_<Vec3P>& points)
{   Array_<int> support; 
    return calcAxisAlignedBoundingBox(points,support); }

/** Alternate signature taking an array of pointers to points rather than the
points themselves. **/
SimTK_SIMMATH_EXPORT static Geo::AlignedBox_<P> 
calcAxisAlignedBoundingBoxIndirect(const Array_<const Vec3P*>& points,
                                   Array_<int>&                support);

/** Alternate signature doesn't return support points. **/
static Geo::AlignedBox_<P> 
calcAxisAlignedBoundingBoxIndirect(const Array_<const Vec3P*>& points)
{   Array_<int> support; 
    return calcAxisAlignedBoundingBoxIndirect(points,support); }

/**@}**/

/**@name                Oriented bounding box creation
These static methods create a tight-fitting oriented bounding box (OBB) that 
includes all of a set of given points. The OBB is not guaranteed to be minimal
but will usually be very good. You can optionally obtain the set of support
points that determined the size of the box. **/
/**@{**/

/** Given a set of points, find the six points that are the most extreme along
specified orientation directions (not necessarily distinct points). The points
are given in an arbitrary frame F. We have an oriented "box" frame B given
by its orientation in F, R_FB. The origin of the B frame is coincident with
the F frame. We'll find the points that are the most extreme along the B frame 
axis directions, and we'll also return the corner points <em>in B</em> (that 
is, the points having minimum and maximum x,y,z values in B). Note that the
corners do not necessarily correspond to any points in the set. If you want to
know where the corners are in F, just compute R_FB*low_B and R_FB*high_B on 
return. There must be at least one point in the given set. **/
SimTK_SIMMATH_EXPORT static void
findOrientedExtremePoints(const Array_<Vec3P>&  points_F, 
                          const RotationP&      R_FB,
                          int least[3], int most[3],
                          Vec3P& low_B, Vec3P& high_B);

/** Alternate signature taking an array of pointers to points rather than the
points themselves. **/
SimTK_SIMMATH_EXPORT static void
findOrientedExtremePointsIndirect(const Array_<const Vec3P*>&  points_F, 
                                  const RotationP&             R_FB,
                                  int least[3], int most[3],
                                  Vec3P& low_B, Vec3P& high_B);

/** Calculate a tight-fitting oriented bounding box (OBB) that includes all
n given points. The OBB is not guaranteed to be minimal but will usually be
very good unless you suppress optimization to save runtime. Cost is O(n). **/
SimTK_SIMMATH_EXPORT static Geo::OrientedBox_<P> 
calcOrientedBoundingBox(const Array_<Vec3P>& points,
                        Array_<int>&         support,
                        bool                 optimize=true);

/** Alternate signature doesn't return support points. **/
static Geo::OrientedBox_<P> 
calcOrientedBoundingBox(const Array_<Vec3P>& points)
{   Array_<int> support; 
    return calcOrientedBoundingBox(points,support); }

/** Alternate signature taking an array of pointers to points rather than the
points themselves. **/
SimTK_SIMMATH_EXPORT static Geo::OrientedBox_<P> 
calcOrientedBoundingBoxIndirect(const Array_<const Vec3P*>& points,
                                Array_<int>&                support,
                                bool                        optimize=true);

/** Alternate signature doesn't return support points. **/
static Geo::OrientedBox_<P> 
calcOrientedBoundingBoxIndirect(const Array_<const Vec3P*>& points,
                                bool                        optimize=true)
{   Array_<int> support; 
    return calcOrientedBoundingBoxIndirect(points,support,optimize); }
/**@}**/

/**@name                 Sphere-related utilities
These static methods work with spheres or collections of spheres.

<h3>Bounding spheres</h3>
Bounding sphere methods calculate the smallest sphere around a given set of
points such that no point is outside the sphere, although some may be on 
its surface. How many and specifically which points were actually used to 
define the sphere can be returned; there will never be more than 4. This 
information is primarily used to construct bounding sphere algorithms; users
normally just need the sphere so can use the simpler signatures.

Bounding sphere methods address roundoff by stretching the sphere enough to 
guarantee that all points are strictly inside the sphere and that later tests
can produce only false positives not false negatives which might cause a
contact to be missed. To do that we have to account not just for machine
precision, but for relative errors caused by spheres of large radius or
spheres that are located far from the origin. These adjustments ensure that
if a test point appears numerically to be outside the sphere, it really cannot 
contact anything that is inside the sphere. 

We use a bounding sphere method due originally to Emo Welzl that computes a 
near-perfect minimal bounding sphere around a set of points with expected O(n) 
run time. Our implementation has been extensively modified to deal with
singular cases so you do not have to precondition the points before asking
for their bounding sphere. 

We also provide a conventional fast and dumb approximate bounding sphere using
Ritter's method as described in Christer Ericson's book. This is mostly 
useful for testing the Welzl method's accuracy and performance and should
not generally be used. A Welzl bounding sphere should never be larger than
a Ritter sphere and should normally be substantially smaller. **/
/**@{**/

/** Create a tiny bounding sphere around a single point. The center is the
point and the radius is tiny but non-zero. **/
static Sphere_<P> calcBoundingSphere(const Vec3P& p)
{   return Sphere_<P>(p, 0).stretchBoundary(); }


/** Create a minimal bounding sphere around two points. Some care is taken
to avoid roundoff problems if the points are far from the origin or very
close together. **/
static Sphere_<P> calcBoundingSphere(const Vec3P& p0, const Vec3P& p1) {
    Array_<int> which;
    return calcBoundingSphere(p0,p1,which);
}


/** Create a minimal bounding sphere around three points. **/
static Sphere_<P> calcBoundingSphere
   (const Vec3P& p0, const Vec3P& p1, const Vec3P& p2) {
    Array_<int> which;
    return calcBoundingSphere(p0,p1,p2,false,which);
}

/** Create a minimal bounding sphere around four points. **/
static Sphere_<P> calcBoundingSphere
   (const Vec3P& p0, const Vec3P& p1, const Vec3P& p2, const Vec3P& p3) {
    Array_<int> which;
    return calcBoundingSphere(p0,p1,p2,p3,false,which);
}

/** Create a minimal bounding sphere around a collection of n points. 
This has expected O(n) performance and usually yields a near-perfect 
bounding sphere. **/
static Sphere_<P> calcBoundingSphere(const Array_<Vec3P>& points) {
    Array_<int> which; 
    return calcBoundingSphere(points, which);
}

/** This signature takes an std::vector rather than a SimTK::Array_; no
extra copying is required. **/
static Sphere_<P> 
calcBoundingSphere(const std::vector<Vec3P>& points) {
    return calcBoundingSphere // no copy done here
                (ArrayViewConst_<Vec3P>(points));
}

/** Create a minimal bounding sphere around a collection of n points, given
indirectly as an array of pointers. This has expected O(n) performance and 
yields a perfect bounding sphere. **/
static Sphere_<P> calcBoundingSphereIndirect(const Array_<const Vec3P*>& points) {
    Array_<int> which; 
    return calcBoundingSphereIndirect(points, which);
}

/** This signature takes an std::vector rather than a SimTK::Array_; no
extra copying is required. **/
static Sphere_<P> 
calcBoundingSphere(const std::vector<const Vec3P*>& points) {
    return calcBoundingSphereIndirect // no copy done here
                (ArrayViewConst_<const Vec3P*>(points));
}

/** Create one-point bounding sphere and return the (trivial) support 
point, of which there is always one. **/
static Sphere_<P> calcBoundingSphere(const Vec3P& p0, Array_<int>& which) 
{   which.clear(); which.push_back(0); 
    return Sphere_<P>(p0,0).stretchBoundary(); }

/** Create a minimum sphere around two points. The center is the 
midpoint, and the radius is roughly half the distance between the points,
possibly expanded in the face of roundoff to ensure that neither point tests
outside. There will be two support points for the circle unless the given 
points are very close to one another. In that case, we treat these as a single 
point and report in \a which that only 1 point was used to define the sphere. 
Points far from the origin will produce a larger sphere because of roundoff.
Cost is about 45 flops. **/
SimTK_SIMMATH_EXPORT static Sphere_<P> 
calcBoundingSphere(const Vec3P& p0, const Vec3P& p1, Array_<int>& which);

/** Create a minimum sphere around three points. There can be 1, 2, or 3
support points returned in \a which. You can optionally force use of the 
3-point circumsphere, which will not always be minimal. Even if
\a forceCircumsphere is set \c true, if the points are
singular (coincident, collinear, coplanar) then it may not be
possible to generate a circumsphere and fewer support points will be used. **/
SimTK_SIMMATH_EXPORT static Sphere_<P> 
calcBoundingSphere(const Vec3P& p0, const Vec3P& p1, const Vec3P& p2, 
                   bool forceCircumsphere, Array_<int>& which);

/** Create a minimum sphere around four points. There can be 1, 2, 3, or 4
support points returned in \a which. You can optionally force use of the 
4-point circumsphere, which will not always be minimal. Even if
\a forceCircumsphere is set \c true, if the points are
singular (coincident, collinear, coplanar, cospherical) then it may not be
possible to generate a circumsphere and fewer support points will be used. **/
SimTK_SIMMATH_EXPORT static Sphere_<P> 
calcBoundingSphere(const Vec3P& p0, const Vec3P& p1, const Vec3P& p2, 
                   const Vec3P& p3, bool forceCircumsphere, Array_<int>& which);

/** Create an optimal minimum sphere around a collection of n points. This has 
expected O(n) performance and yields a near-perfect minimum sphere. There can be
1, 2, 3, or 4 support points used to define the sphere and \a which reports
which of the input points were used. **/
SimTK_SIMMATH_EXPORT static Sphere_<P> 
calcBoundingSphere(const Array_<Vec3P>& points, Array_<int>& which);

/** Alternate signature works with an array of pointers to points. **/
SimTK_SIMMATH_EXPORT static Sphere_<P> 
calcBoundingSphereIndirect(const Array_<const Vec3P*>& points, 
                           Array_<int>& which);

/** Calculate an approximate bounding sphere.\ You should normally use
calcBoundingSphere() which will give a smaller sphere. **/
SimTK_SIMMATH_EXPORT static Sphere_<P> 
calcApproxBoundingSphere(const Array_<Vec3P>& points);

/** This signature takes an std::vector rather than a SimTK::Array_; no
extra copying is required. **/
static Sphere_<P> 
calcApproxBoundingSphere(const std::vector<Vec3P>& points) {
    return calcApproxBoundingSphere // no copy done here
                (ArrayViewConst_<Vec3P>(points));
}

/** Alternate signature works with an array of pointers to points. **/
SimTK_SIMMATH_EXPORT static Sphere_<P> 
calcApproxBoundingSphereIndirect(const Array_<const Vec3P*>& points);

/** This signature takes an std::vector rather than a SimTK::Array_; no
extra copying is required. **/
static Sphere_<P> 
calcApproxBoundingSphereIndirect(const std::vector<const Vec3P*>& points) {
    return calcApproxBoundingSphereIndirect // no copy done here
                (ArrayViewConst_<const Vec3P*>(points));
}
/**@}**/


private:
Vec3P   p;
};

} // namespace SimTK

#endif // SimTK_SIMMATH_GEO_POINT_H_
