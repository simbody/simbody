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
Non-inline static methods from the Geo::Sphere_ class. **/

#include "SimTKcommon.h"
#include "simmath/internal/common.h"
#include "simmath/internal/Geo.h"

namespace SimTK {

//==============================================================================
//                            GEO :: SPHERE
//==============================================================================


//==============================================================================
//                       CALC BOUNDING SPHERE
//==============================================================================

// (Note that the 3-point method from Nicolas Capens computes only the
// sphere that passes through all 3 points, which can be *way* too big.)
// This method produces the minimum bounding sphere for all cases except
// some very, very small triangles. It is modified from Christer Ericson's blog:
// "Minimum bounding circle (sphere) for a triangle (tetrahedron)"
// July 27, 2007. http://realtimecollisiondetection.net/blog/?p=20
// Accessed 12/12/2011 and implemented by Sherm.
// Cost is about 100 flops.
template <class P> /*static*/
Geo::Sphere_<P> Geo::Sphere_<P>::
calcBoundingSphere(const Vec3P& a, const Vec3P& b, const Vec3P& c) {
    const RealP tol = Geo::getDefaultTol<P>();
    const Vec3P ab = b - a, ac = c - a;                     //  6 flops
    const RealP ab2 = ab.normSqr(), ac2 = ac.normSqr();     // 10 flops
    const RealP abac = dot(ab, ac); // == |ab||ac|cos theta   ( 5 flops)

    // The expression below is equivalent to 
    //      d = 2 * (|ab||ac|)^2 * (1-cos^2 theta)
    // where theta is the angle between side ab and ac. That quantity can't
    // be negative except for roundoff error. Note that the area of this 
    // triangle is area = 1/2 sqrt(ab2*ac2-square(abac)) so d = 8*area^2
    const RealP d = 2*(ab2*ac2 - square(abac));             // 4 flops

    RealP rad;
    Vec3P ctr;

    // If d is very small, either the points were nearly collinear or 
    // the triangle is very small. In that case we'll use an approximate
    // method that puts the center at the midpoint of the longest edge.
    // That's perfect for the collinear case and still OK for the very-
    // teeny triangle case.
    if (d <= tol) {                                         // 1 flop
        // Triangle is near singular.
        const Vec3P bc = c - b;
        const RealP bc2 = bc.normSqr();
        RealP d2;    // length squared of longest edge (diameter^2=4*rad^2)
        Vec3P other; // the other possible radial point
        Vec3P ctr2;  // twice the center point
        if (ab2 >= ac2) {
            if (ab2 >= bc2) // ab is longest
                ctr2 = a+b, d2 = ab2, other=c;
            else // bc is longest
                ctr2 = b+c, d2 = bc2, other=a;
        } else { // ac > ab
            if (ac2 >= bc2) // ac is longest
                ctr2 = a+c, d2 = ac2, other=b;
            else // bc is longest
                ctr2 = b+c, d2 = bc2, other=a;
        }
        ctr = ctr2/2;
        rad = std::sqrt(std::max(d2/4, (other-ctr).normSqr()));
    } else {
        // Triangle is non-singular. 
        const RealP ood = 1/d; // (> 0)  ~10 flops
        // s controls height over ac, t over ab, 1-s-t over bc
        const RealP s = (ab2*ac2 - ac2*abac) * ood; // 4 flops
        if (s <= 0) // 1 flop
            return calcBoundingSphere(a,c);
        const RealP t = (ab2*ac2 - ab2*abac) * ood; // 4 flops
        if (t <= 0) // 1 flop
            return calcBoundingSphere(a,b);
        if (s+t >= 1) // 2 flops
            return calcBoundingSphere(b,c); // 1-s-t <= 0

        // Must calculate circumsphere.
        const Vec3P aToCenter = s*ab + t*ac; // 9 flops
        ctr = a + aToCenter;                 // 3 flops
        rad = aToCenter.norm();              // ~35 flops
    }

    // Stretch the sphere out a little to avoid roundoff troubles later.
    return Sphere_<P>(ctr, (1+tol)*rad);    // 2 flops
}

// Create a minimum bounding sphere around four points (a.k.a. a tetrahedron). 
// Method due to Nicolas Capens at 
// http://www.flipcode.com/archives/Smallest_Enclosing_Spheres.shtml
// modified by Peter Eastman and Michael Sherman.
template <class P> /*static*/
Geo::Sphere_<P> Geo::Sphere_<P>::
calcBoundingSphere(const Vec3P& p0, const Vec3P& p1, 
                   const Vec3P& p2, const Vec3P& p3) {
    const RealP tol = Geo::getDefaultTol<P>();
    const Vec3P a = p1-p0;
    const Vec3P b = p2-p0;
    const Vec3P c = p3-p0;
    RealP denom = 2*(a[0]*b[1]*c[2] + a[1]*b[2]*c[0] + a[2]*b[0]*c[1] 
                    - a[2]*b[1]*c[0] - a[0]*b[2]*c[1] - a[1]*b[0]*c[2]);

    RealP radius;
    Vec3P center;

    // If the triangle is very small use an approximate method that just
    // puts a circle about the triangle's center that passes through the
    // furthest-out vertex.
    if (denom <= tol) {
        center = (p0+p1+p2+p3)/4;
        const RealP rsq = std::max((p0-center).normSqr(), 
                          std::max((p1-center).normSqr(), 
                          std::max((p2-center).normSqr(),
                                   (p3-center).normSqr())));
        radius = std::sqrt(rsq);
    } else {
        const Vec3P o = (c.normSqr()*(a % b) +
                         b.normSqr()*(c % a) +
                         a.normSqr()*(b % c)) / denom;
        radius = o.norm();
        center = p0+o;
    }

    // Stretch the sphere out a little to avoid roundoff troubles later.
    return Sphere_<P>(center, (1+tol)*radius);
}

// Explicit instantiations for float and double.
template class Geo::Sphere_<float>;
template class Geo::Sphere_<double>;


}  // End of namespace SimTK