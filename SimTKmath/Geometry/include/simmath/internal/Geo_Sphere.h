#ifndef SimTK_SIMMATH_GEO_SPHERE_H_
#define SimTK_SIMMATH_GEO_SPHERE_H_

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
Defines primitive operations on spheres. **/

#include "SimTKcommon.h"
#include "simmath/internal/common.h"
#include "simmath/internal/Geo.h"

#include <cassert>
#include <cmath>
#include <algorithm>

namespace SimTK {


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
Sphere_(const Vec3P& center, RealP radius)
:   cr(center[0], center[1], center[2], radius) {assert(radius>=0);}
/** Change the radius of this sphere. **/
Sphere_& setRadius(RealP radius) 
{   assert(radius>=0); cr[3]=radius; return *this; }
/** Change the center location of this sphere. **/
Sphere_& setCenter(const Vec3P& center) 
{   Vec3P::updAs(&cr[0])=center; return *this; }

/** Modify this sphere to scale its radius by a fractional amount f,
that is we set radius to f*radius.
@return A reference to this now-resized sphere. **/
Sphere_& scaleBy(RealP f)
{   setRadius(f*getRadius()); return *this; }

/** Stretch this sphere in place by a small amount to ensure that there will 
be no roundoff problems if this is used as a bounding sphere. The amount to 
stretch depends on the default tolerance for this precision, the radius, and 
the position of the sphere in space. A very large sphere, or a sphere that is 
very far from the origin, must be stretched more than a small one at the 
origin. Cost is 6 flops.
@see Geo class for tolerance information. **/
Sphere_& stretchBoundary() {
    const RealP tol = Geo::getDefaultTol<P>();
    const RealP maxdim = max(getCenter().abs());
    const RealP scale = std::max(maxdim, getRadius());
    updRadius() += std::max(scale*Geo::getEps<P>(), tol);
    return *this; 
}

/** Return the volume of this sphere (4/3 pi r^3). **/
RealP findVolume() const 
{   return (RealP(4)/3) * NTraits<P>::getPi() * cube(getRadius()); }
/** Return the surface area of this sphere (4 pi r^2). **/
RealP findArea() const 
{   return 4 * NTraits<P>::getPi() * square(getRadius()); }

/** Return true if a given point is strictly outside this sphere. Just touching
the sphere does not qualify. **/
bool isPointOutside(const Vec3P& p) const {
    const RealP r2 = Geo::Point_<P>::findDistanceSqr(p, getCenter());
    return r2 > square(getRadius());
}
/** Return true if a given point is more than a given tolerance outside the
sphere. **/
bool isPointOutside(const Vec3P& p, RealP tol) const {
    assert(tol >= 0);
    const RealP r2 = Geo::Point_<P>::findDistanceSqr(p, getCenter());
    return r2 > square(getRadius()+tol);
}
/** Get the location of the sphere's center point. **/
const Vec3P& getCenter() const {return Vec3P::getAs(&cr[0]);}
/** Get a writable reference to the sphere's center point. **/
Vec3P& updCenter() {return Vec3P::updAs(&cr[0]);}
/** Get the sphere's radius. **/
RealP getRadius() const {return cr[3];}
/** Get a writable reference to the sphere's radius. **/
RealP& updRadius() {return cr[3];}


private:
// Store together to make sure the compiler doesn't introduce any padding.
// cr[0..2] is the center point, cr[3] is the radius
Vec4P   cr;    
};


} // namespace SimTK

#endif // SimTK_SIMMATH_GEO_SPHERE_H_
