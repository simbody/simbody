#ifndef SimTK_SIMMATH_CONTACT_GEOMETRY_H_
#define SimTK_SIMMATH_CONTACT_GEOMETRY_H_

/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman, Michael Sherman                                    *
 * Contributors: Ian Stavness, Andreas Scholz                                                              *
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
Defines the ContactGeometry class and its API-visible local subclasses for
individual contact shapes. **/

#include "SimTKcommon.h"
#include "simmath/internal/common.h"
#include "simmath/internal/OrientedBoundingBox.h"
#include "simmath/internal/Geodesic.h"

#include <cassert>

namespace SimTK {

/** @class SimTK::ContactGeometryTypeId
This is a unique integer type for quickly identifying specific types of 
contact geometry for fast lookup purposes. **/
SimTK_DEFINE_UNIQUE_INDEX_TYPE(ContactGeometryTypeId);

class ContactGeometryImpl;
class OBBTreeNodeImpl;
class OBBTree;
class Plane;



//==============================================================================
//                             CONTACT GEOMETRY
//==============================================================================
/** A ContactGeometry object describes the shape of all or part of the boundary
of a solid object, for the purpose of modeling with Simbody physical 
effects that occur at the surface of that object, such as contact and 
wrapping forces. Surfaces may be finite or infinite (e.g. a halfspace). 
Surfaces may be smooth or discrete (polyhedral). Smooth surfaces are defined
implicitly as f(P)=0 (P=[px,py,pz]), and optionally may provide a surface
parameterization P=f(u,v). An implicit representation is valid for any P;
parametric representations may have limited validity, singular points, or may
be defined only in a local neighborhood.

A variety of operators are implemented by each specific surface type. Some of
these are designed to support efficient implementation of higher-level 
algorithms that deal in pairs of interacting objects, such as broad- and
narrow-phase contact and minimimum-distance calculations.

The idea here is to collect all the important knowledge about a particular
kind of geometric shape in one place, adding operators as needed to support
new algorithms from time to time. 

All surfaces provide these operations:
  - find closest point to a given point
  - find intersection with a given ray
  - find most extreme point in a given direction
  - return the outward-facing surface normal at a point
  - generate a polygonal mesh that approximates the surface
  - return a unique integer id that may be used to quickly determine the 
    concrete type of a generic surface

Finite surfaces provide
  - a bounding sphere that encloses the entire surface
  - a bounding volume hierarchy with tight-fitting leaf nodes containing
    only simple primitives

Smooth surfaces provide
  - Min/max curvatures and directions
  - Calculate a geodesic between two points on the surface, or 
    starting at a point for a given direction and length

Individual surface types generally support additional operations that may
be used by specialized algorithms that know they are working with that 
particular kind of surface. For example, an algorithm for determining 
ellipsoid-halfspace contact is likely to take advantage of special properties
of both surfaces.

We do not require
detailed solid geometry, but neither can the surface be treated without some
information about the solid it bounds. For example, for contact we must know
which side of the surface is the "inside". However, we don't need a fully
consistent treatment of the solid; for ease of modeling we require only that
the surface behave properly in those locations at which it is evaluated at run
time. The required behavior may vary depending on the algorithm using it.

This is the base class for surface handles; user code will typically 
reference one of the local classes it defines instead for specific shapes. **/
class SimTK_SIMMATH_EXPORT ContactGeometry {
public:
class HalfSpace;
class Cylinder;
class Sphere;
class Ellipsoid;
class SmoothHeightMap;
class TriangleMesh;

// TODO
class Cone;
class Torus;

/** Base class default constructor creates an empty handle. **/
ContactGeometry() : impl(0) {}
/** Copy constructor makes a deep copy. **/
ContactGeometry(const ContactGeometry& src);
/** Copy assignment makes a deep copy. **/
ContactGeometry& operator=(const ContactGeometry& src);
/** Base class destructor deletes the implementation object.\ Note that this
is not virtual; handles should consist of just a pointer to the 
implementation. **/
~ContactGeometry();

/** Generate a DecorativeGeometry that matches the shape of this ContactGeometry **/
DecorativeGeometry createDecorativeGeometry() const;

/** Given a point, find the nearest point on the surface of this object. If 
multiple points on the surface are equally close to the specified point, this 
may return any of them.
@param[in]  position    The point in question.
@param[out] inside      On exit, this is set to true if the specified point is 
                        inside this object, false otherwise.
@param[out] normal      On exit, this contains the surface normal at the 
                        returned point.
@return The point on the surface of the object which is closest to the 
specified point. **/
Vec3 findNearestPoint(const Vec3& position, bool& inside, UnitVec3& normal) const;

/** Determine whether this object intersects a ray, and if so, find the 
intersection point.
@param[in]  origin      The position at which the ray begins.
@param[in]  direction   The ray direction.
@param[out] distance    If an intersection is found, the distance from the ray 
                        origin to the intersection point is stored in this. 
                        Otherwise, it is left unchanged.
@param[out] normal      If an intersection is found, the surface normal of the
                        intersection point is stored in this. Otherwise, it is 
                        left unchanged.
@return \c true if an intersection is found, \c false otherwise. **/
bool intersectsRay(const Vec3& origin, const UnitVec3& direction, 
                   Real& distance, UnitVec3& normal) const;

/** Get a bounding sphere which completely encloses this object.
@param[out] center  On exit, this contains the location of the center of the 
                    bounding sphere.
@param[out] radius  On exit, this contains the radius of the bounding 
                    sphere. **/
void getBoundingSphere(Vec3& center, Real& radius) const;

/** Returns \c true if this is a smooth surface, meaning that it can provide
meaningful curvature information and continuous derivatives with respect to its
parameterization. **/
bool isSmooth() const;

/** Compute the principal curvatures and their directions, and the surface 
normal, at a given point on a smooth surface.
@param[in]      point        
    A point at which to compute the curvature.
@param[out]     curvature    
    On return, this will contain the maximum (curvature[0]) and minimum 
    (curvature[1]) curvatures of the surface at the point.
@param[out]     orientation  
    On return, this will contain the orientation of the surface at the given
    point as follows: the x axis along the direction of maximum curvature, the 
    y axis along the direction of minimum curvature, and the z axis along the 
    surface normal. These vectors are expressed in the surface's coordinate 
    frame.

Non-smooth surfaces will not implement this method and will throw an exception
if you call it. **/
void calcCurvature(const Vec3& point, Vec2& curvature, 
                   Rotation& orientation) const;

/** Our smooth surfaces define a function f(P)=0 that provides an implicit 
representation of the surface. P=(x,y,z) is any point in space expressed in 
the surface's coordinate frame S (that is, given by a vector P-So, expressed in
S). The function is positive inside the object, 0 on the surface, and negative 
outside the object. The returned Function object supports first and second 
partial derivatives with respect to the three function arguments x, y, and z.
Evaluation of the function and its derivatives is cheap. 

Non-smooth surfaces will not implement this method and will throw an exception
if you call it. **/
const Function& getImplicitFunction() const;

/** Calculate the value of the implicit surface function, at a given point.
@param[in]      point
    A point at which to compute the surface value.
@return
    The value of the implicit surface function at the point. **/
Real calcSurfaceValue(const Vec3& point) const;

/** Calculate the implicit surface outward facing unit normal at the given
point. This is determined using the implicit surface function gradient
so is undefined if the point is at a flat spot on the surface where the
gradient is zero. **/
UnitVec3 calcSurfaceUnitNormal(const Vec3& point) const;

/** Calculate the gradient of the implicit surface function, at a given point.
@param[in]      point
    A point at which to compute the surface gradient.
@return
    The gradient of the implicit surface function at the point. **/
Vec3 calcSurfaceGradient(const Vec3& point) const;

/** Calculate the hessian of the implicit surface function, at a given point.
@param[in]      point
    A point at which to compute the surface hessian.
@return
    The hessian of the implicit surface function at the point. **/
Mat33 calcSurfaceHessian(const Vec3& point) const;

/** For an implicit surface, return the Gaussian curvature at the given
point (which might not be on the surface). Here is the formula:
<pre>
        ~grad(f) * Adjoint(H) * grad(f)
   Kg = --------------------------------
                |grad(f)|^4
</pre>
where grad(f) is Df/Dx, Hessian H is D grad(f)/Dx and Adjoint is a 3x3
matrix A where A(i,j)=determinant(H with row i and column j removed). 
Ref: Goldman, R. "Curvature formulas for implicit curves and surfaces",
Comp. Aided Geometric Design 22 632-658 (2005).

Gaussian curvature is the product of the two principal curvatures, Kg=k1*k2.
So for example, the Gaussian curvature anywhere on a sphere is 1/r^2. Note
that despite the name, Gaussian curvature has units of 1/length^2 rather than
curvature units of 1/length.

Here is what the (symmetric) adjoint matrix looks like:
<pre>
adjH  =  [ fyy*fzz - fyz^2, fxz*fyz - fxy*fzz, fxy*fyz - fxz*fyy  ]
         [      (1,2),      fxx*fzz - fxz^2,   fxy*fxz - fxx*fyz  ]
         [      (1,3),           (2,3),        fxx*fyy - fxy^2    ]
</pre>
**/
Real calcGaussianCurvature(const Vec3& point) const;

/** For an implicit surface, return the curvature k of the surface at a given
point p in a given direction tp. Make sure the point is on the surface and the 
direction vector lies in the tangent plane and has unit length |tp| = 1. Then
</pre>
k = ~tp * H * tp ,
<pre>
where H is the Hessian matrix evaluated at p. 
**/
Real calcSurfaceCurvatureInDirection(const Vec3& point, const UnitVec3& direction) const;

/** Returns \c true if this surface is known to be convex. This can be true
for smooth or polygonal surfaces. **/
bool isConvex() const;

/** Given a direction expressed in the surface's frame S, return the point P on 
the surface that is the furthest in that direction (or one of those points if
there is more than one). This will be the point such that dot(P-So, direction)
is maximal for the surface (where So is the origin of the surface). This is 
particularly useful for convex surfaces and should be very fast for them. **/
Vec3 calcSupportPoint(UnitVec3 direction) const;

/** ContactTrackerSubsystem uses this id for fast identification of specific
surface shapes. **/
ContactGeometryTypeId getTypeId() const;

/** Calculate surface curvature at a point using differential geometry as 
suggested by Harris 2006, "Curvature of ellipsoids and other surfaces" Ophthal.
Physiol. Opt. 26:497-501, although the equations here come directly from 
Harris' reference Struik 1961, Lectures on Classical Differential Geometry, 
2nd ed. republished by Dover 1988. Equation and page numbers below are from 
Struik.

This method works for any smooth surface for which there is a local (u,v) 
surface parameterization; it is not restricted to ellipsoids or convex shapes, 
and (u,v) must be distinct but do not have to be perpendicular. Both must be 
perpendicular to the surface normal.

First fundamental form:  I  = E du^2 + 2F dudv + G dv^2
<pre>   E = dxdu^2, F = ~dxdu * dxdv, G=dxdv^2  </pre>

Second fundamental form: II = e du^2 + 2f dudv + g dv^2
<pre>   e = - ~d2xdu2 * nn, f = - ~d2xdudv * nn, g = - ~d2xdv2 * nn </pre>

Given a direction t=dv/du, curvature k is
<pre>
         II   e + 2f t + g t^2
     k = -- = ----------------   (eq. 6-3)
         I    E + 2F t + G t^2
</pre>

We want minimum and maximum values for k to get principal curvatures. We can 
find those as the solutions to dk/dt=0.
<pre>   dk/dt = (E + 2Ft + Gt^2)(f+gt) - (e + 2ft + gt^2)(F+Gt) </pre>

When dk/dt=0, k =(f+gt)/(F+Gt) = (e+ft)/(E+Ft) (eq. 6-4). That provides a 
quadratic equation for the two values of t:
<pre>   A t^2 + B t + C = 0     </pre>
where A=Fg-Gf, B=Eg-Ge, C=Ef-Fe  (eq. 6-5a).

In case the u and v tangent directions are the min and max curvature directions
(on a sphere, for example), they must be perpendicular so F=f=0 (eq. 6-6). Then
the curvatures are ku = e/E and kv = g/G (eq. 6-8).

We're going to return principal curvatures kmax and kmin such that kmax >= kmin,
along with the perpendicular tangent unit directions dmax,dmin that are the 
corresponding principal curvature directions, oriented so that (dmax,dmin,nn) 
form a right-handed coordinate frame.

Cost: given a point P, normalized normal nn, unnormalized u,v tangents and 
second derivatives <pre>
    curvatures: ~115
    directions:  ~50
                ----
                ~165
</pre>  **/
static Vec2 evalParametricCurvature(const Vec3& P, const UnitVec3& nn,
                                    const Vec3& dPdu, const Vec3& dPdv,
                                    const Vec3& d2Pdu2, const Vec3& d2Pdv2, 
                                    const Vec3& d2Pdudv,
                                    Transform& X_EP);

/** This utility method is useful for characterizing the relative geometry of
two locally-smooth surfaces in contact, in a way that is useful for later
application of Hertz compliant contact theory for generating forces. We assume
that contact points Q1 on surface1 and Q2 on surface2 have been determined with
the following properties:
    - the surface normals are aligned but opposite
    - points Q1 and Q2 are separated only along the normal (no tangential 
      separation)

Then the local regions near Q1 and Q2 may be fit with paraboloids P1 and P2
that have their origins at Q1 and Q2, and have the same normals and curvatures
at the origins as do the original surfaces. We will behave here as though
Q1 and Q2 are coincident in space at a point Q; imagine sliding them along
the normal until that happens. Now we define the equations of P1 and P2 in
terms of the maximum and minimum curvatures of surface1 and surface2 at Q:<pre>
    P1: -2z = kmax1 x1^2 + kmin1 y1^2
    P2:  2z = kmax2 x2^2 + kmin2 y2^2  
</pre>
Although the origin Q and z direction are shared, the x,y directions for the 
two paraboloids, though in the same plane z=0, are relatively rotated. Note
that the kmins might be negative; the surfaces do not have to be convex. 

For Hertz contact, we need to know the difference (relative) surface
between the two paraboloids. The difference is a paraboloid P with equation
<pre>
    P: -2z = kmax x^2 + kmin y^2
</pre>
It shares the origin Q and z direction (oriented as for P1), but has its
own principal directions x,y which are coplanar with x1,y1 and x2,y2 but
rotated into some unknown composite orientation. The purpose of this method
is to calculate kmax and kmin, and optionally (depending which signature you
call), x and y, the directions of maximum and minimum curvature (resp.). The
curvature directions are also the principal axes of the contact ellipse formed
by the deformed surfaces, so are necessary (for example) if you want to draw
that ellipse.

Cost is about 220 flops. If you don't need the curvature directions, call the
other overloaded signature which returns only kmax and kmin and takes only 
about 1/3 as long. 

@param[in]          R_SP1  
    The orientation of the P1 paraboloid's frame, expressed in some frame S 
    (typically the frame of the surface to which P1 is fixed). R_SP1.x() is 
    the direction of maximum curvature; y() is minimum curvature; z is the
    contact normal pointing away from surface 1.
@param[in]          k1      
    The maximum (k1[0]) and minimum (k1[1]) curvatures for P1 at the contact 
    point Q1 on surface1. Negative curvatures are handled correctly here but
    may cause trouble for your force model if the resulting contact is 
    conforming.
@param[in]          x2      
    The direction of maximum curvature for paraboloid P2. \a x2 must be in the 
    x1,y1 plane provided in \a R_SP1 and expressed in the S frame.
@param[in]          k2      
    The maximum (k2[0]) and minimum (k2[1]) curvatures for P2 at the contact 
    point Q2 on surface2. Negative curvatures are handled correctly here but
    may cause trouble for your force model if the resulting contact is 
    conforming.
@param[out]         R_SP    
    The orientation of the difference paraboloid P's frame, expressed in the 
    same S frame as was used for P1. R_SP.x() is the direction of maximum 
    curvature of P at the contact point; y() is the minimum curvature 
    direction; z() is the unchanged contact normal pointing away from surface1.
@param[out]         k       
    The maximum (k[0]) and minimum(k[1]) curvatures for the difference 
    paraboloid P at the contact point Q. If either of these is negative or
    zero then the surfaces are conforming and you can't use a point contact
    force model. Note that if k1>0 and k2>0 (i.e. surfaces are convex at Q)
    then k>0 too. If some of the surface curvatures are concave, it is still
    possible that k>0, depending on the relative curvatures.

@see The other signature for combineParaboloids() that is much cheaper if
you just need the curvatures \a k but not the directions \a R_SP. **/
static void combineParaboloids(const Rotation& R_SP1, const Vec2& k1,
                               const UnitVec3& x2, const Vec2& k2,
                               Rotation& R_SP, Vec2& k);

/** This is a much faster version of combineParaboloids() for when you just
need the curvatures of the difference paraboloid, but not the directions 
of those curvatures. Cost is about 70 flops. See the other overload of
this method for details. **/
static void combineParaboloids(const Rotation& R_SP1, const Vec2& k1,
                               const UnitVec3& x2, const Vec2& k2,
                               Vec2& k);


/** @name                  Geodesic Evaluators **/
/**@{**/

/** Given two points, find a geodesic curve connecting them.
If a preferred starting point is provided, find the geodesic curve that
is closest to that point. Otherwise, find the shortest length geodesic.

@param[in] xP            Coordinates of the first point.
@param[in] xQ            Coordinates of the second point.
@param[in] xSP           (Optional) Coordinates of a preferred point for the geodesic to be near
@param[in] options       Parameters related to geodesic calculation
@param[out] geod         On exit, this contains a geodesic between P and Q.
**/
void initGeodesic(const Vec3& xP, const Vec3& xQ, const Vec3& xSP,
        const GeodesicOptions& options, Geodesic& geod) const;


/** Given two points and previous geodesic curve close to the points, find
a geodesic curve connecting the points that is close to the previous geodesic.

@param[in] xP            Coordinates of the first point.
@param[in] xQ            Coordinates of the second point.
@param[in] prevGeod      A previous geodesic that should be near the new one.
@param[in] options       Parameters related to geodesic calculation
@param[out] geod         On exit, this contains a geodesic between P and Q.
**/
// XXX if xP and xQ are the exact end-points of prevGeod; then geod = prevGeod;
void continueGeodesic(const Vec3& xP, const Vec3& xQ, const Geodesic& prevGeod,
        const GeodesicOptions& options, Geodesic& geod);


/** Compute a geodesic curve starting at the given point, starting in the
 * given direction, and terminating at the given length.

@param[in] xP            Coordinates of the starting point for the geodesic.
@param[in] tP            The starting tangent direction for the geodesic.
@param[in] terminatingLength   The length that the resulting geodesic should have.
@param[in] options       Parameters related to geodesic calculation
@param[out] geod         On exit, this contains the calculated geodesic
**/
// XXX what to do if tP is not in the tangent plane at P -- project it?
void shootGeodesicInDirectionUntilLengthReached
   (const Vec3& xP, const UnitVec3& tP, const Real& terminatingLength, 
    const GeodesicOptions& options, Geodesic& geod) const;

/** Given an already-calculated geodesic on this surface connecting points
P and Q, fill in the sensitivity of point P with respect to a change of
tangent direction at Q. If there are interior points stored with the geodesic,
then we'll calculate the interior sensitivities also.

@param[in,out]  geodesic        An already-calculated geodesic.
@param[in]      initSensitivity
    Initial conditions for the Jacobi field calculation. If this is the whole
    geodesic then the initial conditions are (0,1) for the sensitivity and
    its arc length derivative. However, if we are continuing from another
    geodesic, then the end sensitivity for that geodesic is the initial 
    conditions for this one. 
**/
void calcGeodesicReverseSensitivity
   (Geodesic& geodesic,
    const Vec2& initSensitivity = Vec2(0,1)) const; // j, jdot at end point


/** Compute a geodesic curve starting at the given point, starting in the
 * given direction, and terminating when it hits the given plane.

@param[in] xP            Coordinates of the starting point for the geodesic.
@param[in] tP            The starting tangent direction for the geodesic.
@param[in] terminatingPlane   The plane in which the end point of the resulting geodesic should lie.
@param[in] options       Parameters related to geodesic calculation
@param[out] geod         On exit, this contains the calculated geodesic
**/
// XXX what to do if tP is not in the tangent plane at P -- project it?
// XXX what to do if we don't hit the plane
void shootGeodesicInDirectionUntilPlaneHit(const Vec3& xP, const UnitVec3& tP,
        const Plane& terminatingPlane, const GeodesicOptions& options,
        Geodesic& geod) const;


/** Utility method to find geodesic between P and Q using split geodesic 
method with initial shooting directions tPhint and -tQhint. **/
void calcGeodesic(const Vec3& xP, const Vec3& xQ,
        const Vec3& tPhint, const Vec3& tQhint, Geodesic& geod) const;

/** Utility method to find geodesic between P and Q using the orthogonal
method, with initial direction tPhint and initial length lengthHint. **/
void calcGeodesicUsingOrthogonalMethod(const Vec3& xP, const Vec3& xQ,
        const Vec3& tPhint, Real lengthHint, Geodesic& geod) const;

/** This signature makes a guess at the initial direction and length
and then calls the other signature. **/
void calcGeodesicUsingOrthogonalMethod(const Vec3& xP, const Vec3& xQ,
        Geodesic& geod) const
{
    const Vec3 r_PQ = xQ - xP;
    const Real lengthHint = r_PQ.norm();
    const UnitVec3 n = calcSurfaceUnitNormal(xP);
    // Project r_PQ into the tangent plane.
    const Vec3 t_PQ = r_PQ - (~r_PQ*n)*n;
    const Real tLength = t_PQ.norm();
    const UnitVec3 tPhint =
        tLength != 0 ? UnitVec3(t_PQ/tLength, true)
                     : n.perp(); // some arbitrary perpendicular to n
    calcGeodesicUsingOrthogonalMethod(xP, xQ, Vec3(tPhint), lengthHint, geod);           
}


/**
 * Utility method to calculate the "geodesic error" between one geodesic
 * shot from P in the direction tP and another geodesic shot from Q in the
 * direction tQ. We optionally return the resulting "kinked" geodesic in
 * case anyone wants it; if the returned error is below tolerance then that
 * geodesic is the good one.
 **/
Vec2 calcSplitGeodError(const Vec3& P, const Vec3& Q,
                   const UnitVec3& tP, const UnitVec3& tQ,
                   Geodesic* geod=0) const;



/** Analytically compute a geodesic curve starting at the given point, starting in the
 * given direction, and terminating at the given length. Only possible for a few simple
 * shapes, such as spheres and cylinders.

@param[in] xP            Coordinates of the starting point for the geodesic.
@param[in] tP            The starting tangent direction for the geodesic.
@param[in] terminatingLength   The length that the resulting geodesic should have.
@param[in] options       Parameters related to geodesic calculation
@param[out] geod         On exit, this contains the calculated geodesic
**/
// XXX what to do if tP is not in the tangent plane at P -- project it?
void shootGeodesicInDirectionUntilLengthReachedAnalytical
   (const Vec3& xP, const UnitVec3& tP, const Real& terminatingLength,
    const GeodesicOptions& options, Geodesic& geod) const;


/** Analytically compute a geodesic curve starting at the given point, starting in the
 * given direction, and terminating when it hits the given plane. Only possible
 * for a few simple shapes, such as spheres and cylinders.

@param[in] xP            Coordinates of the starting point for the geodesic.
@param[in] tP            The starting tangent direction for the geodesic.
@param[in] terminatingPlane   The plane in which the end point of the resulting geodesic should lie.
@param[in] options       Parameters related to geodesic calculation
@param[out] geod         On exit, this contains the calculated geodesic
**/
// XXX what to do if tP is not in the tangent plane at P -- project it?
// XXX what to do if we don't hit the plane
void shootGeodesicInDirectionUntilPlaneHitAnalytical(const Vec3& xP, const UnitVec3& tP,
        const Plane& terminatingPlane, const GeodesicOptions& options,
        Geodesic& geod) const;


/** Utility method to analytically find geodesic between P and Q with initial shooting
 directions tPhint and tQhint. Only possible for a few simple shapes, such as spheres
 and cylinders.
 **/
void calcGeodesicAnalytical(const Vec3& xP, const Vec3& xQ,
        const Vec3& tPhint, const Vec3& tQhint, Geodesic& geod) const;

/**
 * Utility method to analytically calculate the "geodesic error" between one geodesic
 * shot from P in the direction tP and another geodesic shot from Q in the
 * direction tQ. We optionally return the resulting "kinked" geodesic in
 * case anyone wants it; if the returned error is below tolerance then that
 * geodesic is the good one. Only possible for a few simple shapes, such as spheres
 and cylinders.
 **/
Vec2 calcSplitGeodErrorAnalytical(const Vec3& P, const Vec3& Q,
                   const UnitVec3& tP, const UnitVec3& tQ,
                   Geodesic* geod=0) const;

/**@}**/


/** @name Geodesic-related Debugging **/
/**@{**/


/**
 * Compute rotation matrix using the normal at the given point and the
 * given direction.
 **/
Rotation calcTangentBasis(const Vec3& point, const Vec3& dir) {
    const UnitVec3 n = calcSurfaceUnitNormal(point);
    Rotation R_GS;
    R_GS.setRotationFromTwoAxes(n, ZAxis, dir, XAxis);
    return R_GS;
}

/** Get the plane associated with the
    geodesic hit plane event handler  **/
const Plane& getPlane() const;
/** Set the plane associated with the
    geodesic hit plane event handler  **/
void setPlane(const Plane& plane) const;
/** Get the geodesic for access by visualizer **/
const Geodesic& getGeodP() const;
/** Get the geodesic for access by visualizer **/
const Geodesic& getGeodQ() const;
const int getNumGeodesicsShot() const;
void addVizReporter(ScheduledEventReporter* reporter) const;
/**@}**/



explicit ContactGeometry(ContactGeometryImpl* impl); /**< Internal use only. **/
bool isOwnerHandle() const;                          /**< Internal use only. **/
bool isEmptyHandle() const;                          /**< Internal use only. **/
bool hasImpl() const {return impl != 0;}             /**< Internal use only. **/
/** Internal use only. **/
const ContactGeometryImpl& getImpl() const {assert(impl); return *impl;}
/** Internal use only. **/
ContactGeometryImpl& updImpl() {assert(impl); return *impl; }

protected:
ContactGeometryImpl* impl; /**< Internal use only. **/
};



//==============================================================================
//                                 HALF SPACE
//==============================================================================
/** This ContactGeometry subclass represents an object that occupies the 
entire half-space x>0. This is useful for representing walls and floors. **/
class SimTK_SIMMATH_EXPORT ContactGeometry::HalfSpace : public ContactGeometry {
public:
HalfSpace();

/** Return true if the supplied ContactGeometry object is a halfspace. **/
static bool isInstance(const ContactGeometry& geo)
{   return geo.getTypeId()==classTypeId(); }
/** Cast the supplied ContactGeometry object to a const halfspace. **/
static const HalfSpace& getAs(const ContactGeometry& geo)
{   assert(isInstance(geo)); return static_cast<const HalfSpace&>(geo); }
/** Cast the supplied ContactGeometry object to a writable halfspace. **/
static HalfSpace& updAs(ContactGeometry& geo)
{   assert(isInstance(geo)); return static_cast<HalfSpace&>(geo); }

/** Obtain the unique id for HalfSpace contact geometry. **/
static ContactGeometryTypeId classTypeId();

class Impl; /**< Internal use only. **/
const Impl& getImpl() const; /**< Internal use only. **/
Impl& updImpl(); /**< Internal use only. **/
};



//==============================================================================
//                                CYLINDER
//==============================================================================
/** This ContactGeometry subclass represents a cylinder centered at the
origin, with radius r in the x-y plane, and infinite length along z.
TODO: should allow finite length to be specified. **/
class SimTK_SIMMATH_EXPORT ContactGeometry::Cylinder : public ContactGeometry {
public:
explicit Cylinder(Real radius);
Real getRadius() const;
void setRadius(Real radius);

/** Return true if the supplied ContactGeometry object is a sphere. **/
static bool isInstance(const ContactGeometry& geo)
{   return geo.getTypeId()==classTypeId(); }
/** Cast the supplied ContactGeometry object to a const sphere. **/
static const Cylinder& getAs(const ContactGeometry& geo)
{   assert(isInstance(geo)); return static_cast<const Cylinder&>(geo); }
/** Cast the supplied ContactGeometry object to a writable sphere. **/
static Cylinder& updAs(ContactGeometry& geo)
{   assert(isInstance(geo)); return static_cast<Cylinder&>(geo); }

/** Obtain the unique id for Cylinder contact geometry. **/
static ContactGeometryTypeId classTypeId();

class Impl; /**< Internal use only. **/
const Impl& getImpl() const; /**< Internal use only. **/
Impl& updImpl(); /**< Internal use only. **/
};



//==============================================================================
//                                  SPHERE
//==============================================================================
/** This ContactGeometry subclass represents a sphere centered at the 
origin. **/
class SimTK_SIMMATH_EXPORT ContactGeometry::Sphere : public ContactGeometry {
public:
explicit Sphere(Real radius);
Real getRadius() const;
void setRadius(Real radius);

/** Return true if the supplied ContactGeometry object is a sphere. **/
static bool isInstance(const ContactGeometry& geo)
{   return geo.getTypeId()==classTypeId(); }
/** Cast the supplied ContactGeometry object to a const sphere. **/
static const Sphere& getAs(const ContactGeometry& geo)
{   assert(isInstance(geo)); return static_cast<const Sphere&>(geo); }
/** Cast the supplied ContactGeometry object to a writable sphere. **/
static Sphere& updAs(ContactGeometry& geo)
{   assert(isInstance(geo)); return static_cast<Sphere&>(geo); }

/** Obtain the unique id for Sphere contact geometry. **/
static ContactGeometryTypeId classTypeId();

class Impl; /**< Internal use only. **/
const Impl& getImpl() const; /**< Internal use only. **/
Impl& updImpl(); /**< Internal use only. **/
};



//==============================================================================
//                                  ELLIPSOID
//==============================================================================
/** This ContactGeometry subclass represents an ellipsoid centered at the 
origin, with its principal axes pointing along the x, y, and z axes and
half dimensions a,b, and c (all > 0) along those axes, respectively. The
implicit equation f(x,y,z)=0 of the ellipsoid surface is <pre>
    f(x,y,z) = Ax^2+By^2+Cz^2 - 1
    where A=1/a^2, B=1/b^2, C=1/c^2     
</pre>
A,B, and C are the squares of the principal curvatures ka=1/a, kb=1/b, and
kc=1/c.

The interior of the ellipsoid consists of all points such that f(x,y,z)<0 and
points exterior satisfy f(x,y,z)>0. The region around any point (x,y,z) on an 
ellipsoid surface is locally an elliptic paraboloid with equation <pre>
    -2 z' = kmax x'^2 + kmin y'^2   
</pre>
where z' is measured along the the outward unit normal n at (x,y,z), x' is
measured along the the unit direction u of maximum curvature, and y' is
measured along the unit direction v of minimum curvature. kmax,kmin are the 
curvatures with kmax >= kmin > 0. The signs of the mutually perpendicular
vectors u and v are chosen so that (u,v,n) forms a right-handed coordinate 
system for the paraboloid. **/
class SimTK_SIMMATH_EXPORT ContactGeometry::Ellipsoid : public ContactGeometry {
public:
/** Construct an Ellipsoid given its three principal half-axis dimensions a,b,c
(all positive) along the local x,y,z directions respectively. The curvatures 
(reciprocals of radii) are precalculated here at a cost of about 30 flops. **/
explicit Ellipsoid(const Vec3& radii);
/** Obtain the three half-axis dimensions a,b,c used to define this
ellipsoid. **/
const Vec3& getRadii() const;
/** Set the three half-axis dimensions a,b,c (all positive) used to define this
ellipsoid, overriding the current radii and recalculating the principal 
curvatures at a cost of about 30 flops. 
@param[in] radii    The three half-dimensions of the ellipsoid, in the 
                    ellipsoid's local x, y, and z directions respectively. **/
void setRadii(const Vec3& radii);

/** For efficiency we precalculate the principal curvatures whenever the 
ellipsoid radii are set; this avoids having to repeatedly perform these three
expensive divisions at runtime. The curvatures are ka=1/a, kb=1/b, and kc=1/c 
so that the ellipsoid's implicit equation can be written Ax^2+By^2+Cz^2=1, 
with A=ka^2, etc. **/
const Vec3& getCurvatures() const;

/** Given a point \a P =(x,y,z) on the ellipsoid surface, return the unique unit
outward normal to the ellipsoid at that point. If \a P is not on the surface, 
the result is the same as for the point obtained by scaling the vector 
\a P - O until it just touches the surface. That is, we compute 
P'=findPointInThisDirection(P) and then return the normal at P'. Cost is about
40 flops regardless of whether P was initially on the surface. 
@param[in] P    A point on the ellipsoid surface, measured and expressed in the
                ellipsoid's local frame. See text for what happens if \a P is
                not actually on the ellipsoid surface.
@return The outward-facing unit normal at point \a P (or at the surface point
pointed to by \a P).
@see findPointInSameDirection() **/
UnitVec3 findUnitNormalAtPoint(const Vec3& P) const;

/** Given a unit direction \a n, find the unique point P on the ellipsoid 
surface at which the outward-facing normal is \a n. Cost is about 40 flops. 
@param[in] n    The unit vector for which we want to find a match on the
                ellipsoid surface, expressed in the ellipsoid's local frame. 
@return The point on the ellipsoid's surface at which the outward-facing
normal is the same as \a n. The point is measured and expressed in the 
ellipsoid's local frame. **/
Vec3 findPointWithThisUnitNormal(const UnitVec3& n) const;

/** Given a direction d defined by the vector Q-O for an arbitrary point in 
space Q=(x,y,z)!=O, find the unique point P on the ellipsoid surface that is 
in direction d from the ellipsoid origin O. That is, P=s*d for some scalar 
s > 0 such that f(P)=0. Cost is about 40 flops. 
@param[in] Q    A point in space measured from the ellipsoid origin but not the
                origin.
@return P, the intersection of the ray in the direction Q-O with the ellipsoid
surface **/
Vec3 findPointInSameDirection(const Vec3& Q) const;

/** Given a point Q on the surface of the ellipsoid, find the approximating
paraboloid at Q in a frame P where OP=Q, Pz is the outward-facing unit
normal to the ellipsoid at Q, Px is the direction of maximum curvature
and Py is the direction of minimum curvature. k=(kmax,kmin) are the returned
curvatures with kmax >= kmin > 0. The equation of the resulting paraboloid 
in the P frame is -2z = kmax*x^2 + kmin*y^2. Cost is about 260 flops; you can
save a little time if you already know the normal at Q by using the other
overloaded signature for this method.
 
@warning It is up to you to make sure that Q is actually on the ellipsoid
surface. If it is not you will quietly get a meaningless result.

@param[in]  Q       A point on the surface of this ellipsoid, measured and
                    expressed in the ellipsoid's local frame.
@param[out] X_EP    The frame of the paraboloid P, measured and expressed in
                    the ellipsoid local frame E. X_EP.p() is \a Q, X_EP.x() 
                    is the calculated direction of maximum curvature kmax; y() 
                    is the direction of minimum curvature kmin; z is the 
                    outward facing normal at \a Q.
@param[out] k       The maximum (k[0]) and minimum (k[1]) curvatures of the
                    ellipsoid (and paraboloid P) at point \a Q.
@see findParaboloidAtPointWithNormal() **/
void findParaboloidAtPoint(const Vec3& Q, Transform& X_EP, Vec2& k) const;

/** If you already have both a point and the unit normal at that point, this 
will save about 40 flops by trusting that you have provided the correct normal;
be careful -- no one is going to check that you got this right. The results are
meaningless if the point and normal are not consistent. Cost is about 220 flops.
@see findParaboloidAtPoint() for details **/
void findParaboloidAtPointWithNormal(const Vec3& Q, const UnitVec3& n,
    Transform& X_EP, Vec2& k) const;

/** Return true if the supplied ContactGeometry object is an Ellipsoid. **/
static bool isInstance(const ContactGeometry& geo)
{   return geo.getTypeId()==classTypeId(); }
/** Cast the supplied ContactGeometry object to a const Ellipsoid. **/
static const Ellipsoid& getAs(const ContactGeometry& geo)
{   assert(isInstance(geo)); return static_cast<const Ellipsoid&>(geo); }
/** Cast the supplied ContactGeometry object to a writable Ellipsoid. **/
static Ellipsoid& updAs(ContactGeometry& geo)
{   assert(isInstance(geo)); return static_cast<Ellipsoid&>(geo); }

/** Obtain the unique id for Ellipsoid contact geometry. **/
static ContactGeometryTypeId classTypeId();

class Impl; /**< Internal use only. **/
const Impl& getImpl() const; /**< Internal use only. **/
Impl& updImpl(); /**< Internal use only. **/
};



//==============================================================================
//                            SMOOTH HEIGHT MAP
//==============================================================================
/** This ContactGeometry subclass represents a smooth surface fit through a
set of sampled points using bicubic patches to provide C2 continuity. It is
particularly useful as a bounded terrain. The boundary is an axis-aligned
rectangle in the local x-y plane. Within the boundary, every (x,y) location has
a unique height z, so caves and overhangs cannot be represented.

The surface is parameterized as z=f(x,y) where x,y,z are measured in 
the surface's local coordinate frame. This can also be described as the 
implicit function F(x,y,z)=f(x,y)-z=0, as though this were an infinitely thick
slab in the -z direction below the surface. **/
class SimTK_SIMMATH_EXPORT 
ContactGeometry::SmoothHeightMap : public ContactGeometry {
public:
/** Create a SmoothHeightMap surface using an already-existing BicubicSurface.
The BicubicSurface object is referenced, not copied, so it may be shared with
other users. **/
explicit SmoothHeightMap(const BicubicSurface& surface);

/** Return a reference to the BicubicSurface object being used by this
SmoothHeightMap. **/
const BicubicSurface& getBicubicSurface() const;

/** (Advanced) Return a reference to the oriented bounding box tree for this
surface. **/
const OBBTree& getOBBTree() const;

/** Return true if the supplied ContactGeometry object is a SmoothHeightMap. **/
static bool isInstance(const ContactGeometry& geo)
{   return geo.getTypeId()==classTypeId(); }
/** Cast the supplied ContactGeometry object to a const SmoothHeightMap. **/
static const SmoothHeightMap& getAs(const ContactGeometry& geo)
{   assert(isInstance(geo)); return static_cast<const SmoothHeightMap&>(geo); }
/** Cast the supplied ContactGeometry object to a writable SmoothHeightMap. **/
static SmoothHeightMap& updAs(ContactGeometry& geo)
{   assert(isInstance(geo)); return static_cast<SmoothHeightMap&>(geo); }

/** Obtain the unique id for SmoothHeightMap contact geometry. **/
static ContactGeometryTypeId classTypeId();

class Impl; /**< Internal use only. **/
const Impl& getImpl() const; /**< Internal use only. **/
Impl& updImpl(); /**< Internal use only. **/
};



//==============================================================================
//                              TRIANGLE MESH
//==============================================================================
/** This ContactGeometry subclass represents an arbitrary shape described by a 
mesh of triangular faces. The mesh surface must satisfy the following 
requirements:
  - It must be closed, so that any point can unambiguously be classified as 
    either inside or outside.
  - It may not intersect itself anywhere, even at a single point.
  - It must be an oriented manifold.
  - The vertices for each face must be ordered counter-clockwise when viewed
    from the outside. That is, if v0, v1, and v2 are the locations of the 
    three vertices for a face, the cross product (v1-v0)%(v2-v0) must point 
    outward.
  - The length of every edge must be non-zero.

It is your responsibility to ensure that any mesh you create meets these 
requirements. The constructor will detect many incorrect meshes and signal them
by throwing an exception, but it is not guaranteed to detect all possible 
problems.  If a mesh fails to satisfy any of these requirements, the results of
calculations performed with it are undefined. For example, collisions involving
it might fail to be detected, or contact forces on it might be calculated 
incorrectly. **/
class SimTK_SIMMATH_EXPORT ContactGeometry::TriangleMesh 
:   public ContactGeometry {
public:
class OBBTreeNode;
/** Create a TriangleMesh.
@param vertices     The positions of all vertices in the mesh.
@param faceIndices  The indices of the vertices that make up each face. The 
                    first three elements are the vertices in the first face, 
                    the next three elements are the vertices in the second 
                    face, etc.
@param smooth       If true, the mesh will be treated as a smooth surface, and 
                    normal vectors will be smoothly interpolated between 
                    vertices. If false, it will be treated as a faceted mesh 
                    with a constant normal vector over each face. **/
TriangleMesh(const ArrayViewConst_<Vec3>& vertices, const ArrayViewConst_<int>& faceIndices, bool smooth=false);
/** Create a TriangleMesh based on a PolygonalMesh object. If any faces of the 
PolygonalMesh have more than three vertices, they are automatically 
triangulated.
@param mesh      The PolygonalMesh from which to construct a triangle mesh.
@param smooth    If true, the mesh will be treated as a smooth surface, and 
                 normal vectors will be smoothly interpolated between vertices.
                 If false, it will be treated as a faceted mesh with a constant
                 normal vector over each face. **/
explicit TriangleMesh(const PolygonalMesh& mesh, bool smooth=false);
/** Get the number of edges in the mesh. **/
int getNumEdges() const;
/** Get the number of faces in the mesh. **/
int getNumFaces() const;
/** Get the number of vertices in the mesh. **/
int getNumVertices() const;
/** Get the position of a vertex in the mesh.
@param index  The index of the vertex to get.
@return The position of the specified vertex. **/
const Vec3& getVertexPosition(int index) const;
/** Get the index of one of the edges of a face. Edge 0 connects vertices 0 
and 1. Edge 1 connects vertices 1 and 2. Edge 2 connects vertices 0 and 2.
@param face    The index of the face.
@param edge    The index of the edge within the face (0, 1, or 2).
@return The index of the specified edge. **/
int getFaceEdge(int face, int edge) const;
/** Get the index of one of the vertices of a face.
@param face    The index of the face.
@param vertex  The index of the vertex within the face (0, 1, or 2).
@return The index of the specified vertex. **/
int getFaceVertex(int face, int vertex) const;
/** Get the index of one of the faces shared by an edge
@param edge    The index of the edge.
@param face    The index of the face within the edge (0 or 1).
@return The index of the specified face. **/
int getEdgeFace(int edge, int face) const;
/** Get the index of one of the vertices shared by an edge.
@param edge    The index of the edge.
@param vertex  The index of the vertex within the edge (0 or 1).
@return The index of the specified vertex. **/
int getEdgeVertex(int edge, int vertex) const;
/** Find all edges that intersect a vertex.
@param vertex  The index of the vertex.
@param edges   The indices of all edges intersecting the vertex will be added
               to this. **/
void findVertexEdges(int vertex, Array_<int>& edges) const;
/** Get the normal vector for a face. This points outward from the mesh.
@param face    The index of the face. **/
const UnitVec3& getFaceNormal(int face) const;
/** Get the area of a face.
@param face    The index of the face. **/
Real getFaceArea(int face) const;
/** Calculate the location of a point on the surface, in the local frame of
the TriangleMesh. Cost is 11 flops. 
@param face    The index of the face containing the point.
@param uv      The point within the face, specified by its barycentric uv 
               coordinates. **/
Vec3 findPoint(int face, const Vec2& uv) const;
/** Calculate the location of a face's centroid, that is, the point uv=(1/3,1/3)
which is the average of the three vertex locations. This is a common special 
case of findPoint() that can be calculated more quickly (7 flops).
@param face    The index of the face whose centroid is of interest. **/
Vec3 findCentroid(int face) const;
/** Calculate the normal vector at a point on the surface.
@param face    The index of the face containing the point.
@param uv      The point within the face, specified by its barycentric uv 
               coordinates. **/
UnitVec3 findNormalAtPoint(int face, const Vec2& uv) const;
/** Given a point, find the nearest point on the surface of this object. If 
multiple points on the surface are equally close to the specified point, this 
may return any of them.
@param position    The point in question.
@param inside      On exit, this is set to true if the specified point is 
                   inside this object, false otherwise.
@param normal      On exit, this contains the surface normal at the returned 
                   point.
@return The point on the surface of the object which is closest to the 
specified point. **/
Vec3 findNearestPoint(const Vec3& position, bool& inside, UnitVec3& normal) const;
/** Given a point, find the nearest point on the surface of this object. If 
multiple points on the surface are equally close to the specified point, this 
may return any of them.
@param position    The point in question.
@param inside      On exit, this is set to true if the specified point is 
                   inside this object, false otherwise.
@param face        On exit, this contains the index of the face containing the 
                   returned point.
@param uv          On exit, this contains the barycentric coordinates (u and v)
                   of the returned point within its face.
@return The point on the surface of the object which is closest to the 
specified point. **/
Vec3 findNearestPoint(const Vec3& position, bool& inside, int& face, Vec2& uv) const;

/** Given a point and a face of this object, find the point of the face that is
nearest the given point. If multiple points on the face are equally close to 
the specified point, this may return any of them.
@param position    The point in question.
@param face        The face to be examined.
@param uv          On exit, this contains the barycentric coordinates (u and v)
                   of the returned point within the face.
@return The face point, in the surface's frame, that is closest to the 
specified point. **/
Vec3 findNearestPointToFace(const Vec3& position, int face, Vec2& uv) const;


/** Determine whether this mesh intersects a ray, and if so, find the 
intersection point.
@param origin     The position at which the ray begins.
@param direction  The ray direction.
@param distance   If an intersection is found, the distance from the ray origin
                  to the intersection point is stored in this. Otherwise, it is
                  left unchanged.
@param normal     If an intersection is found, the surface normal of the 
                  intersection point is stored in this. Otherwise, it is left 
                  unchanged.
@return \c true if an intersection is found, \c false otherwise. **/
bool intersectsRay(const Vec3& origin, const UnitVec3& direction, Real& distance, UnitVec3& normal) const;
/** Determine whether this mesh intersects a ray, and if so, find what face it 
hit.
@param origin     The position at which the ray begins.
@param direction  The ray direction.
@param distance   If an intersection is found, the distance from the ray origin
                  to the intersection point is stored in this. Otherwise, it is
                  left unchanged.
@param face       If an intersection is found, the index of the face hit by the
                  ray is stored in this. Otherwise, it is left unchanged.
@param uv         If an intersection is found, the barycentric coordinates (u 
                  and v) of the intersection point within the hit face are 
                  stored in this. Otherwise, it is left unchanged.
@return \c true if an intersection is found, \c false otherwise. **/
bool intersectsRay(const Vec3& origin, const UnitVec3& direction, Real& distance, int& face, Vec2& uv) const;
/** Get the OBBTreeNode which forms the root of this mesh's Oriented Bounding 
Box Tree. **/
OBBTreeNode getOBBTreeNode() const;

/** Generate a PolygonalMesh from this TriangleMesh; useful mostly for debugging
because you can create a DecorativeMesh from this and then look at it. **/
PolygonalMesh createPolygonalMesh() const;

/** Return true if the supplied ContactGeometry object is a triangle mesh. **/
static bool isInstance(const ContactGeometry& geo)
{   return geo.getTypeId()==classTypeId(); }
/** Cast the supplied ContactGeometry object to a const triangle mesh. **/
static const TriangleMesh& getAs(const ContactGeometry& geo)
{   assert(isInstance(geo)); return static_cast<const TriangleMesh&>(geo); }
/** Cast the supplied ContactGeometry object to a writable triangle mesh. **/
static TriangleMesh& updAs(ContactGeometry& geo)
{   assert(isInstance(geo)); return static_cast<TriangleMesh&>(geo); }

/** Obtain the unique id for TriangleMesh contact geometry. **/
static ContactGeometryTypeId classTypeId();

class Impl; /**< Internal use only. **/
const Impl& getImpl() const; /**< Internal use only. **/
Impl& updImpl(); /**< Internal use only. **/
};



//==============================================================================
//                       TRIANGLE MESH :: OBB TREE NODE
//==============================================================================
/** This class represents a node in the Oriented Bounding Box Tree for a 
TriangleMesh. Each node has an OrientedBoundingBox that fully encloses all 
triangles contained within it or its  children. This is a binary tree: each 
non-leaf node has two children. Triangles are stored only in the leaf nodes. **/
class SimTK_SIMMATH_EXPORT ContactGeometry::TriangleMesh::OBBTreeNode {
public:
OBBTreeNode(const OBBTreeNodeImpl& impl);
/** Get the OrientedBoundingBox which encloses all triangles in this node or 
its children. **/
const OrientedBoundingBox& getBounds() const;
/** Get whether this is a leaf node. **/
bool isLeafNode() const;
/** Get the first child node. Calling this on a leaf node will produce an 
exception. **/
const OBBTreeNode getFirstChildNode() const;
/** Get the second child node. Calling this on a leaf node will produce an 
exception. **/
const OBBTreeNode getSecondChildNode() const;
/** Get the indices of all triangles contained in this node. Calling this on a
non-leaf node will produce an exception. **/
const Array_<int>& getTriangles() const;
/** Get the number of triangles inside this node. If this is not a leaf node,
this is the total number of triangles contained by all children of this
node. **/
int getNumTriangles() const;

private:
const OBBTreeNodeImpl* impl;
};


//==============================================================================
//                     GEODESIC EVALUATOR helper classes
//==============================================================================


/**
 * A simple plane class
 **/
class Plane {
public:
    Plane() : m_normal(1,0,0), m_offset(0) { }
    Plane(const Vec3& normal, const Real& offset)
    :   m_normal(normal), m_offset(offset) { }

    Real getDistance(const Vec3& pt) const {
        return ~m_normal*pt - m_offset;
    }

    Vec3 getNormal() const {
        return m_normal;
    }

    Real getOffset() const {
        return m_offset;
    }

private:
    Vec3 m_normal;
    Real m_offset;
}; // class Plane


/**
 * A event handler to terminate integration when geodesic hits the plane.
 * For use with a ParticleOnSurfaceSystem
 **/
class GeodHitPlaneEvent : public TriggeredEventHandler {
public:
    GeodHitPlaneEvent()
    :   TriggeredEventHandler(Stage::Position) { }

    explicit GeodHitPlaneEvent(const Plane& aplane)
    :   TriggeredEventHandler(Stage::Position) {
        plane = aplane;
    }

    // event is triggered if distance of geodesic endpoint to plane is zero
    Real getValue(const State& state) const {
        if (!enabled) {
            return 1;
        }
        Vec3 endpt(&state.getQ()[0]);
        Real dist =  plane.getDistance(endpt);
//        std::cout << "dist = " << dist << std::endl;
        return dist;
    }

    // This method is called whenever this event occurs.
    void handleEvent(State& state, Real accuracy, bool& shouldTerminate) const {
        if (!enabled) {
            return;
        }

        // This should be triggered when geodesic endpoint to plane is zero.
        Vec3 endpt;
        const Vector& q = state.getQ();
        endpt[0] = q[0]; endpt[1] = q[1]; endpt[2] = q[2];
        Real dist = plane.getDistance(endpt);

//        ASSERT(std::abs(dist) < 0.01 );
        shouldTerminate = true;
//        std::cout << "hit plane!" << std::endl;
    }

    void setPlane(const Plane& aplane) const {
        plane = aplane;
    }

    const Plane& getPlane() const {
        return plane;
    }

    const void setEnabled(bool enabledFlag) {
        enabled = enabledFlag;
    }

    const bool isEnabled() {
        return enabled;
    }

private:
    mutable Plane plane;
    bool enabled;

}; // class GeodHitPlaneEvent

/**
 * This class generates decoration for contact points and straight line path segments
 **/
class PathDecorator : public DecorationGenerator {
public:
    PathDecorator(const Vector& x, const Vec3& O, const Vec3& I, const Vec3& color) :
            m_x(x), m_O(O), m_I(I), m_color(color) { }

    virtual void generateDecorations(const State& state,
            Array_<DecorativeGeometry>& geometry) {
//        m_system.realize(state, Stage::Position);

        Vec3 P, Q;
        P[0] = m_x[0]; P[1] = m_x[1]; P[2] = m_x[2];
        Q[0] = m_x[3]; Q[1] = m_x[4]; Q[2] = m_x[5];

        geometry.push_back(DecorativeSphere(0.05).setColor(Black).setTransform(m_O));
        geometry.push_back(DecorativeSphere(0.05).setColor(Black).setTransform(P));
        geometry.push_back(DecorativeSphere(0.05).setColor(Black).setTransform(Q));
        geometry.push_back(DecorativeSphere(0.05).setColor(Black).setTransform(m_I));

        geometry.push_back(DecorativeLine(m_O,P)
                .setColor(m_color)
                .setLineThickness(2));
        geometry.push_back(DecorativeLine(Q,m_I)
                .setColor(m_color)
                .setLineThickness(2));

    }

private:
    const Vector& m_x; // x = ~[P Q]
    const Vec3& m_O;
    const Vec3& m_I;
    const Vec3& m_color;
    Rotation R_plane;
    Vec3 offset;
}; // class DecorationGenerator


/**
 * This class generates decoration for a plane
 **/
class PlaneDecorator : public DecorationGenerator {
public:
    PlaneDecorator(const Plane& plane, const Vec3& color) :
            m_plane(plane), m_color(color) { }

    virtual void generateDecorations(const State& state,
            Array_<DecorativeGeometry>& geometry) {
//        m_system.realize(state, Stage::Position);

        // draw plane
        R_plane.setRotationFromOneAxis(UnitVec3(m_plane.getNormal()),
                CoordinateAxis::XCoordinateAxis());
        offset = 0;
        offset[0] = m_plane.getOffset();
        geometry.push_back(
                DecorativeBrick(Vec3(0.01,1,1))
                .setTransform(Transform(R_plane, R_plane*offset))
                .setColor(m_color)
                .setOpacity(.2));
    }

private:
    const Plane& m_plane;
    const Vec3& m_color;
    Rotation R_plane;
    Vec3 offset;
}; // class DecorationGenerator


} // namespace SimTK

#endif // SimTK_SIMMATH_CONTACT_GEOMETRY_H_
