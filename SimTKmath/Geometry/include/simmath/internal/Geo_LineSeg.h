#ifndef SimTK_SIMMATH_GEO_LINESEG_H_
#define SimTK_SIMMATH_GEO_LINESEG_H_

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
Collects primitive operations involving line segments. **/

#include "SimTKcommon.h"
#include "simmath/internal/common.h"
#include "simmath/internal/Geo.h"
#include "simmath/internal/Geo_Sphere.h"

#include <cassert>
#include <cmath>
#include <algorithm>

namespace SimTK {


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
{   return Geo::Point_<P>::calcBoundingSphere(e[0],e[1]); }

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


} // namespace SimTK

#endif // SimTK_SIMMATH_GEO_LINESEG_H_
