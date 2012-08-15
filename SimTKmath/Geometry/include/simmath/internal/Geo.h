#ifndef SimTK_SIMMATH_GEO_H_
#define SimTK_SIMMATH_GEO_H_

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

The Geo class itself is dataless and provides only static methods. It is also
used as a namespace for geometric primitives to allow these names to be used
elsewhere for more significant objects. **/
class SimTK_SIMMATH_EXPORT Geo {
public:
template <class P> class Point_;
template <class P> class Sphere_;
template <class P> class LineSeg_;
template <class P> class Line_;
template <class P> class Plane_;
template <class P> class Circle_;
template <class P> class Box_;
template <class P> class AlignedBox_;
template <class P> class OrientedBox_;
template <class P> class Triangle_;
template <class P> class CubicHermiteCurve_;
template <class P> class BicubicHermitePatch_;
template <class P> class CubicBezierCurve_;
template <class P> class BicubicBezierPatch_;

typedef Point_<Real>        Point;
typedef Sphere_<Real>       Sphere;
typedef LineSeg_<Real>      LineSeg;
typedef Line_<Real>         Line;
typedef Plane_<Real>        Plane;
typedef Circle_<Real>       Circle;
typedef Box_<Real>          Box;
typedef AlignedBox_<Real>   AlignedBox;
typedef OrientedBox_<Real>  OrientedBox;
typedef Triangle_<Real>     Triangle;
typedef CubicHermiteCurve_<Real>    CubicHermiteCurve;
typedef BicubicHermitePatch_<Real>  BicubicHermitePatch;
typedef CubicBezierCurve_<Real>     CubicBezierCurve;
typedef BicubicBezierPatch_<Real>   BicubicBezierPatch;

/**@name          Differential geometry of curve segments

These methods calculate geometric quantities from given parametric ones, in an
arbitrary parameter u, such that points on the curve segment are given by 
P(u), with 0<=u<=1. We consider the direction of increasing arc length to be
the same as the direction of u. These are utility methods that can be used
with any parametric curve. The idea is that you use the curve evaluators
to determine the arguments here, using the same u value for each of them. If
you don't use the same u value you'll get meaningless results.

The primary reference for this material is the book "Lectures on Classical
Differential Geometry, 2nd ed." (chapter 1) by Dirk Struik, 1961, republished
by Dover in 1988. Notation and equation numbers are from that reference.
Two exceptions: (1) we use c for the curvature vector rather than Struik's bold
k, so we can use k for the scalar value of curvature, and (2) we define the 
curve normal n to point \e away (outward) from the center of curvature, while 
Struik defined it to point inward. Using our definition, if you have a
parametric circle the normal points towards the outside, which is more 
conventional and analogous to surface normals. Our right-handed curve frame 
is thus x,y,z=n,t,b rather than Struik's "moving trihedron" frame t,n,b. For 
us the binormal b=n X t, while with Struik's definition it is t X n. Note that 
only the normal is reversed from Struik's, the tangent and binormal vectors 
are the same. **/
/**@{**/

/** Given the parametric derivative Pu(u)=dP/du, determine whether the point
P(u) is at a cusp, that is, a place where the arc length s does not change when
parameter u does, so ds/du=0 (within tolerance). Note that an inflection point,
where the curvature is zero, is not a cusp; see isInflectionPoint(). Cost is 6 
flops. **/
template <class RealP, int S> static bool
isCusp(const Vec<3,RealP,S>& Pu) 
{   return Pu.normSqr() < getDefaultTolSqr<RealP>(); }

/** Given the parametric derivatives Pu(u)=dP/du, and Puu(u)=d2P/du2 determine 
whether point P(u) is at an inflection point on the curve, that is, a flat 
place where there is no curvature (within tolerance). We will also return true
if this point is a cusp, where curvature is not defined. Cost is 15 flops.

<h3>Theory</h3>
If Puu is zero or parallel to Pu (meaning it changes the parametric tangent's
length but not its direction) then we are at an inflection point. We consider 
a cusp to be an inflection point even though curvature is undefined there, 
so we define an inflection point to be any place where |Pu X Puu|==0 to 
within tolerance.

@see isCusp() **/
template <class RealP, int S> static bool
isInflectionPoint(const Vec<3,RealP,S>& Pu, const Vec<3,RealP,S>& Puu) 
{ return (Pu % Puu).normSqr() < getDefaultTolSqr<RealP>(); }

/** Calculate the unit tangent vector t=dP/ds, given Pu=dP/du. This is
undefined at a cusp (Pu==0). See Struik, eq. 2-2.
Cost is about 40 flops. **/
template <class RealP, int S> static UnitVec<RealP,1> 
calcUnitTangent(const Vec<3,RealP,S>& Pu) {
    const RealP dsdu = Pu.norm();               // ~25 flops
    SimTK_ERRCHK(dsdu >= getDefaultTol<RealP>(), 
        "Geo::calcUnitTangent()", "Unit tangent undefined at a cusp.");

    return UnitVec<RealP,1>(Pu/dsdu, true);     // ~13 flops
}

/** Return the curvature vector c=dt/ds=d2P/ds2, given Pu=dP/du and 
Puu=d2P/du2. This vector points \e towards the center of curvature, with
length equal to the magnitude of the curvature (it's not a unit vector). 
Curvature is undefined at a cusp (where Pu==0).
Since we define the curve unit normal n to point \e away
from the center of curvature (see above), we have c=-k*n where k is the 
(scalar) curvature. Cost is about 30 flops. 

<h3>Theory</h3>
See Struik, eqn. 4-3, 4-4. Let prime denote differentiation with respect to 
arclength: <pre>
    u' = du/ds = 1/|Pu|
    u'' = -(~Pu Puu)/Pu^4 = -(~Pu Puu) * u'^4
    t = P' = Pu u' = Pu/|Pu|
    c = t' = P'' = Puu u'^2 + Pu u''
      = Puu/Pu^2 - Pu (~Pu Puu)/Pu^4
      = -k * n, k is signed curvature, n is unit normal
</pre> 
Note that c can be used to determine the magnitude |k|, but you can't get
the sign until the normal n has been defined. We're going to define n from c,
so that n=-c/|c| so dot(c,n) is always negative meaning that k is always 
positive for us. **/
template <class RealP, int S> static Vec<3,RealP> 
calcCurvatureVector(const Vec<3,RealP,S>& Pu, const Vec<3,RealP,S>& Puu) {
    const RealP Pu2 = Pu.normSqr();                 // (ds/du)^2, 5 flops
    SimTK_ERRCHK(Pu2 >= getDefaultTolSqr<RealP>(), 
        "Geo::calcCurvatureVector()", "Curvature undefined at a cusp.");
    const RealP PuPuu     = dot(Pu,Puu);
    const RealP uPrimeSqr = 1/Pu2;                          // ~10 flops
    const RealP u2Prime   = -PuPuu * square(uPrimeSqr);     //   8 flops
    return uPrimeSqr*Puu + u2Prime*Pu;                      //   9 flops
}

/** In our definition, the unit normal vector n points in the "outward" 
direction, that is, it points away from the center of curvature (opposite
the curvature vector c so n=-c/|c|). This convention is the opposite of 
Struik's, where he has the normal point in the same direction
as the curvature vector. The normal is undefined at a cusp (Pu(u)==0), and is
an arbitrary perpendicular to the tangent at an inflection point (a flat 
place where there is no curvature, i.e. |Pu(u) X Puu(u)|==0). If the curve is 
a straight line then every point is an inflection point, so the normal is 
arbitrary everywhere. Cost is about 80 flops. 
@see calcCurvatureVector() for theory. **/
template <class RealP, int S> static UnitVec<RealP,1> 
calcUnitNormal(const Vec<3,RealP,S>& Pu, const Vec<3,RealP,S>& Puu) {
    typedef UnitVec<RealP,1> UnitVec3P;

    const RealP Pu2 = Pu.normSqr();                 // (ds/du)^2, 5 flops
    SimTK_ERRCHK(Pu2 >= getDefaultTolSqr<RealP>(), 
        "Geo::calcUnitNormal()", "The normal is undefined at a cusp.");

    // Now check if we're at an inflection point, meaning |Pu X Puu|==0.
    // Use this handy identity from Struik eq. 3-9 to avoid calculating
    // the cross product: (Pu X Puu)^2 = Pu^2Puu^2 - (~Pu*Puu)^2
    const RealP Puu2 = Puu.normSqr();                   // 5 flops
    const RealP PuPuu = dot(Pu,Puu);                    // 5 flops
    const RealP PuXPuu2 = Pu2*Puu2 - square(PuPuu);     // 3 flops
    if (PuXPuu2 < getDefaultTolSqr<RealP>())            // 1 flop
        return UnitVec3P(Pu).perp();
    // Calculate the curvature vector, negate, and normalize.
    const RealP uPrimeSqr = 1/Pu2;                      // ~10 flops
    const RealP u2Prime   = -PuPuu * square(uPrimeSqr); //   3 flops
    const Vec<3,RealP> c = uPrimeSqr*Puu + u2Prime*Pu;  //   9 flops
    return UnitVec3P(-c);                               // ~40 flops
}

/** Return the the curvature k (always positive), and a frame
whose origin is a point along the curve, x axis is the outward unit normal n,
y is the unit tangent t, and z=x X y is the binormal b, which is a normal
to the osculating plane. So the vectors n,t,b form a right-handed set; this
convention is different from Struik's since he has n pointing the opposite
direction. This frame is undefined at a cusp (Pu==0), and the normal is
arbitrary at an inflection point (|Pu(u) X Puu(u)|==0) or if the 
curve is a line. Cost is about 115 flops.
**/
template <class RealP, int S> static RealP 
calcCurveFrame(const Vec<3,RealP,S>&    P, 
               const Vec<3,RealP,S>&    Pu, 
               const Vec<3,RealP,S>&    Puu,
               Transform_<RealP>&       X_FP) {
    typedef UnitVec<RealP,1> UnitVec3P;

    const RealP Pu2 = Pu.normSqr();                 // (ds/du)^2, 5 flops
    SimTK_ERRCHK(Pu2 >= getDefaultTolSqr<RealP>(), 
        "Geo::calcCurveFrame()", "Curve frame is undefined at a cusp.");

    // Set the point P(u) as the frame origin.
    X_FP.updP() = P;

    // Calculate the unit tangent t, our y axis.
    const RealP uPrimeSqr = 1/Pu2;                          // ~10 flops
    const RealP uPrime    = std::sqrt(uPrimeSqr);           // ~20 flops
    const UnitVec3P t(uPrime*Pu, true);                     //   3 flops 

    // Next calculate unit normal n, our x axis. See calcUnitNormal() above
    // for theory.
    const RealP Puu2 = Puu.normSqr();                       // 5 flops
    const RealP PuPuu = dot(Pu,Puu);                        // 5 flops
    const RealP PuXPuu2 = Pu2*Puu2 - square(PuPuu);         // 3 flops
    UnitVec3P n; // unit normal
    RealP     k; // curvature magnitude
    if (PuXPuu2 < getDefaultTolSqr<RealP>()) {              // 1 flop
        k = 0;
        n = t.perp(); // arbitrary
    } else {
        // Calculate the curvature vector, negate, and normalize.
        const RealP u2Prime = -PuPuu * square(uPrimeSqr);   //   8 flops
        const Vec<3,RealP> c = uPrimeSqr*Puu + u2Prime*Pu;  //   9 flops
        k = c.norm();                       // curvature >= 0, ~25 flops
        n = UnitVec3P((-1/k)*c, true);                      // ~13 flops
    }

    // Finally calculate the binormal, our z axis. No need to normalize
    // here because n and t are perpendicular unit vectors.
    const UnitVec3P b(n % t, true);                      // 9 flops

    // Construct the coordinate frame without normalizing.
    X_FP.updR().setRotationFromUnitVecsTrustMe(n,t,b);

    return k;
}

/** Return k^2, the square of the scalar curvature k, given Pu=dP/du and
Puu=d2P/du2. Using our definition for the curve normal n (see above), k is 
always positive so the curvature is the positive square root of the value 
returned here. Curvature is undefined at a cusp (where Pu==0) and is zero
at an inflection point (|Pu X Puu|==0). Cost is about 30 flops. 

<h3>Theory</h3>
<pre>
                  |Pu X Puu|^2
    k^2 = |c|^2 = ------------
                     |Pu|^6
</pre> See Struik, pg. 17, eq. 5-5a. **/
template <class RealP, int S> static RealP
calcCurvatureSqr(const Vec<3,RealP,S>& Pu, const Vec<3,RealP,S>& Puu) {
    const RealP Pu2 = Pu.normSqr();   // (ds/du)^2, 5 flops
    SimTK_ERRCHK(Pu2 >= getDefaultTolSqr<RealP>(), 
        "Geo::calcCurvatureSqr()", "Curvature is undefined at a cusp.");
    const RealP num = cross(Pu,Puu).normSqr();      //  14 flops
    const RealP den = cube(Pu2);                    //   2 flops
    return num/den;                                 // ~10 flops
}

/** Return tau, the torsion or "second curvature" given Pu=dP/du, 
Puu=d2P/du2, Puuu=d3P/du3. Torsion is a signed quantity related to the
rate of change of the osculating plane binormal b, with db/ds=tau*n where
n is the "outward" unit normal (see above). Torsion is undefined at either a 
cusp (where Pu==0) or an inflection point (where |Pu X Puu|==0). Cost is about
30 flops. 

<h3>Theory</h3>
<pre>
          ~(Pu X Puu) * Puuu
    tau = ------------------
             |Pu X Puu|^2
</pre> See Struik, pg. 17, eq. 5-5b and discussion on page 16. **/
template <class RealP, int S> static RealP
calcTorsion(const Vec<3,RealP,S>& Pu, const Vec<3,RealP,S>& Puu, 
            const Vec<3,RealP,S>& Puuu) {
    const Vec<3,RealP> PuXPuu = cross(Pu,Puu);              //   9 flops
    const RealP PuXPuu2 = PuXPuu.normSqr();                 //   5 flops  
    SimTK_ERRCHK(PuXPuu2 >= getDefaultTolSqr<RealP>(), "Geo::calcTorsion()", 
        "Torsion is undefined at a cusp or inflection point.");
    const RealP num = dot(PuXPuu, Puuu);                    //   5 flops
    return num/PuXPuu2;                                     // ~10 flops
}
/**@}**/

/**@name                 Lines **/

/** Find the points of closest approach on two lines L0 and L1, each 
represented by an origin point and a direction. Points and vectors must be in 
a common frame. We return point x0 on L0 and x1 on L1 such that the distance
|x1-x0| is the smallest for any points on the two lines. If the lines are 
parallel or nearly so (all points same distance) we'll pick the point
on each line closest to midway between the origins as the closest points and 
return and indication that the returned points weren't unique. 

@param[in]      p0      The origin point of line L0, that is, any point
                            through which line L0 passes.
@param[in]      d0      A unit vector giving the direction of L0.
@param[in]      p1      The origin point of line L1.
@param[in]      d1      A unit vector giving the direction of L1.
@param[out]     x0      The point of L0 that is closest to L1.
@param[out]     x1      The point of L1 that is closest to L2.
@param[out]     linesAreParallel 
                        True if the lines were treated as effectively parallel.

Cost is about 65 flops. **/
template <class RealP> static void
findClosestPointsOfTwoLines
   (const Vec<3,RealP>& p0, const UnitVec<RealP,1>& d0,
    const Vec<3,RealP>& p1, const UnitVec<RealP,1>& d1,
    Vec<3,RealP>& x0, Vec<3,RealP>& x1, bool& linesAreParallel)
{
    const Vec<3,RealP> w = p1-p0; // vector from p0 to p1; 3 flops
    const RealP s2Theta = (d0 % d1).normSqr(); // sin^2(angle); 14 flops
    const RealP d =  dot(w,d0);     // 5 flops
    const RealP e = -dot(w,d1);     // 6 flops

    RealP t0, t1; // line parameters of closest points
    if (s2Theta < square(NTraits<RealP>::getSignificant())) { // 3 flops
        // Lines are parallel. Return a point on each line midway between
        // the origin points.
        linesAreParallel = true; // parallel
        t0 = d/2; t1 = e/2;      // 2 flops
    } else {
        linesAreParallel = false;
        const RealP cTheta = dot(d0,d1); // cos(angle between lines); 5 flops
        const RealP oos2Theta = 1/s2Theta; // about 10 flops
        t0 = (e*cTheta + d) * oos2Theta;   // 3 flops
        t1 = (d*cTheta + e) * oos2Theta;   // 3 flops
    }

    x0 = p0 + t0*d0;    // 6 flops
    x1 = p1 + t1*d1;    // 6 flops
}

/**@name                 Miscellaneous utilities **/
/**@{**/

/** Return the default tolerance to use for degeneracy tests and other tests
for "too small" or "near enough" that arise in dealing with geometry primitives.
The value depends on the precision being used; we use the SimTK constant
SignificantReal which is eps^(7/8) where eps is the resolution of the
template argument RealP (which must be \c float or \c double) That
makes this tolerance around 2e-14 in double precision and 9e-7 in float. **/
template <class RealP> static RealP getDefaultTol() 
{   return NTraits<RealP>::getSignificant(); }
/** Returns the square of the default tolerance. 
@see getDefaultTol() **/
template <class RealP> static RealP getDefaultTolSqr()
{   return square(getDefaultTol<RealP>()); }

/** Return machine precision for floating point calculations at precision
RealP. **/
template <class RealP> static RealP getEps() 
{   return NTraits<RealP>::getEps(); }
/** Return a NaN (not a number) at precision RealP. **/
template <class RealP> static RealP getNaN() 
{   return NTraits<RealP>::getNaN(); }
/** Return Infinity at precision RealP.\ You can negate this for -Infinity. **/
template <class RealP> static RealP getInfinity() 
{   return NTraits<RealP>::getInfinity(); }

/** Stretch a dimension by a given tolerance amount. The result is the
given \a length increased by at least an absolute amount \a tol, or by a
relative amount length*tol if length > 1. Don't call this with \a tol less
than machine precision or an exception will be thrown. Cost is 3 flops. **/
template <class RealP> static RealP stretchBy(RealP length, RealP tol) {
    SimTK_ERRCHK2(tol >= getEps<RealP>(), "Geo::stretchBy()",
        "The supplied tolerance %g is too small; must be at least %g"
        " for significance at this precision.", 
        (double)tol, (double)getEps<RealP>());

    return length + std::max(length*tol, tol);
}

/** Stretch a dimension using the default tolerance for this precision as
the tolerance in stretchBy(). Cost is 3 flops. **/
template <class RealP> static RealP stretch(RealP length)
{   return stretchBy(length, getDefaultTol<RealP>()); }

/**@}**/

};


} // namespace SimTK

#endif // SimTK_SIMMATH_GEO_H_
