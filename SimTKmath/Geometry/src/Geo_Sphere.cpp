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

// These local helpers find the minimum or maximum of several values and 
// return which one was the extreme value.
template <class P>
inline static void minOf(P a, P b, P& minVal, int& which) {
    minVal=a; which=0;
    if (b<minVal) minVal=b, which=1;
}

template <class P>
inline static void minOf(P a, P b, P c, P& minVal, int& which) {
    minVal=a; which=0;
    if (b<minVal) minVal=b, which=1;
    if (c<minVal) minVal=c, which=2;
}

template <class P>
inline static void minOf(P a, P b, P c, P d, P& minVal, int& which) {
    minVal=a; which=0;
    if (b<minVal) minVal=b, which=1;
    if (c<minVal) minVal=c, which=2;
    if (d<minVal) minVal=d, which=3;
}
template <class P>
inline static void maxOf(P a, P b, P& maxVal, int& which) {
    maxVal=a; which=0;
    if (b>maxVal) maxVal=b, which=1;
}

template <class P>
inline static void maxOf(P a, P b, P c, P& maxVal, int& which) {
    maxVal=a; which=0;
    if (b>maxVal) maxVal=b, which=1;
    if (c>maxVal) maxVal=c, which=2;
}

template <class P>
inline static void maxOf(P a, P b, P c, P d, P& maxVal, int& which) {
    maxVal=a; which=0;
    if (b>maxVal) maxVal=b, which=1;
    if (c>maxVal) maxVal=c, which=2;
    if (d>maxVal) maxVal=d, which=3;
}

// This helper replaces the 0-based indices in "which" with corresponding
// indices in "map". We assume map is long enough.
inline static void fixWhich(const int* map, Array_<int>& which) {
    for (unsigned i=0; i<which.size(); ++i)
        which[i] = map[which[i]];
}

//==============================================================================
//                       CALC MINIMUM SPHERE - 2 POINTS
//==============================================================================
/* I bet you think this is easy! Unfortunately if the points in question are
far from the origin (at 100, say) but close together, we can't compute the 
separation very accurately (around 100*eps for this example). So the resulting
sphere can be too big (annoying but OK) or too small (very bad). Also, if the
line length is <= tol, we want to consider these a single point rather than
two supporting points for the sphere; in that case we pick one of the points
and center the sphere around that with radius tol. About 45 flops. */
template <class P> /*static*/
Geo::Sphere_<P> Geo::Sphere_<P>::
calcMinimumSphere(const Vec3P& p0, const Vec3P& p1, Array_<int>& which) {
    const RealP tol = Geo::getDefaultTol<P>();

    // Choose the tentative center, subject to scaled roundoff.
    const Vec3P ctr = (p0 + p1)/2; // 6 flops

    // Rather than using (p0-p1)/2 as the radius which can result in a sphere
    // that doesn't include the points by a large margin (I tried it), measure 
    // the actual distances from the center to each point and use the larger 
    // as the radius. If that radius is <= tol/2 (meaning line length <= tol) 
    // then we'll just move the center to the closer point, set the radius to 
    // tol and return with 1 support point (i.e., we consider this a point not 
    // a line).
    const RealP p0d2 = (p0-ctr).normSqr();  // squared dist from center
    const RealP p1d2 = (p1-ctr).normSqr();  // (8 flops each)

    RealP rad;
    which.clear();
    if (p0d2 >= p1d2) {
        rad = std::sqrt(p0d2);          // 20 flops
        if (rad <= tol/2) {             //  2 flops
            which.push_back(1);
            return Sphere_(p1, tol); 
        }
    } else { // p1d2 > p0d2
        rad = std::sqrt(p1d2);          // same cost here
        if (rad <= tol/2) {
            which.push_back(0);
            return Sphere_(p0, tol); 
        }
    }
    // use both points
    which.push_back(0); which.push_back(1);
    return Sphere_(ctr, rad); 
}

//==============================================================================
//                       CALC BOUNDING SPHERE - 3 POINTS
//==============================================================================
/* (Note that the 3- and 4-point methods from Nicolas Capens compute only the
sphere that passes through all points (the circumsphere), which can be *way* 
too big.)

This method produces the minimum bounding sphere around three points (i.e.,
a triangle) for all cases except some very, very small triangles where it will
give a small, but perhaps not minimal, bounding sphere. It is 
modified from Christer Ericson's blog: "Minimum bounding circle (sphere) for 
a triangle (tetrahedron)" July 27, 2007. 
http://realtimecollisiondetection.net/blog/?p=20
Accessed 12/12/2011 and reimplemented by Sherm. The main change is handling of
singular triangles and roundoff problems.

Cost is about 110 flops for a typical case.

Implementation
--------------
We have a triangle with vertices a,b,c with no restrictions on their placement 
(i.e., they can be singular in several different ways). Algorithm:
    - Calculate the area of the triangle.
    - If the area is very small:
         - Find the longest edge from ab, ac, bc.
         - Drop the vertex that is not used in the longest edge.
         - Calculate the bounding sphere of the longest edge.
         - If the dropped vertex is outside, stretch the sphere to include it.
         - Return that sphere as the result (might not be perfect but will
           be very good).
    - Else (area is OK):
         - Calculate barycentric coordinates s,t, u=1-s-t of the 
           circumsphere center P (the point equidistant to all vertices).
         - If s,t, or u < 0, P lies outside the triangle and a 2-point sphere 
           will be smaller. Drop the corresponding vertex and return the 
           bounding sphere for the remaining edge (2 supporting points).
         - Otherwise return the circumsphere (3 supporting points).

How to calculate the circumsphere:
    We express center P as P=a + s*(b-a) + t*(c-a). This point is to be 
    equidistant from all vertices, giving us two equations for s and t:
           (P-b)^2 = (P-a)^2, (P-c)^2 = (P-a)^2
    Substituting the definition of P gives us the following system of equations
    for s and t:
           [ab^2  abac] [s]   [ab^2/2]
           [abac  ac^2] [t] = [ac^2/2]   or   M x = b
    We'll use Cramer's rule to solve the equation, which requires computing 
    determinants of three matrices: M, Ms, Mt where the latter are M with its 
    1st or 2nd column replaced by b. Let dm=det(M), ds=det(Ms), dt=det(Mt). 
    Then s=ds/dm, t=dt/dm. Determinant dm is related to the area of the 
    triangle by 2*a = sqrt(dm), so dm = 4 a^2. Note that dm can't be negative
    except possibly due to roundoff error. Below we'll actually calculate
    2*dm, 2*ds, etc. since that saves some multiplies.
*/
template <class P> /*static*/
Geo::Sphere_<P> Geo::Sphere_<P>::
calcMinimumSphere(const Vec3P& a, const Vec3P& b, const Vec3P& c,
                  Array_<int>& which) {
    const RealP tol = Geo::getDefaultTol<P>();

    const Vec3P ab = b - a, ac = c - a;                     //  6 flops
    const RealP ab2 = ab.normSqr(), ac2 = ac.normSqr();     // 10 flops
    const RealP abac = dot(ab, ac); // == |ab||ac|cos theta     5 flops
    const RealP ab2ac2 = ab2*ac2;                           //  1 flop

    // The expression below is equivalent to 
    //      dm2 = 2*dm = 2 * (|ab||ac|)^2 * (1-cos^2 theta)
    // where theta is the angle between side ab and ac. That quantity can't
    // be negative except for roundoff error. Note that the area of this 
    // triangle is area = 1/2 sqrt(ab2*ac2-square(abac)) so dm2 = 8*area^2
    const RealP dm2 = 2*(ab2ac2 - square(abac));             //  3 flops

    // We'll use one of three ways to find the center point and determine
    // the set of supporting vertices (in which):
    //  1) singular (2 or 1 support vertices)
    //  2) non-singular with 2 support vertices
    //  3) non-singular using circumsphere (3 support vertices)

    Vec3P ctr;

    // Consider the triangle singular if its area is less than tol.
    if (dm2 <= 8*square(tol)) {                               // 3 flops
        // Triangle is near singular or very small.
        const Vec3P bc = c - b;
        const RealP bc2 = bc.normSqr();
        // Find the longest edge.
        RealP maxLen2; int maxEdge;
        maxOf(ab2,ac2,bc2,maxLen2,maxEdge);
        // Now create a sphere around that edge using 2 or 1 support points,
        // but we're only going to keep the center.
        Sphere_<P> edgeSphere; int map[2];
        switch(maxEdge) {
        case 0: map[0]=0; map[1]=1; // (a,b) drop c
            edgeSphere=calcMinimumSphere(a,b,which); break;
        case 1: map[0]=0; map[1]=2; // (a,c) drop b
            edgeSphere=calcMinimumSphere(a,c,which); break;
        case 2: map[0]=1; map[1]=2; // (b,c) drop a
            edgeSphere=calcMinimumSphere(b,c,which); break;
        };
        fixWhich(map, which); // fix the indices
        ctr = edgeSphere.getCenter();
    } else {
        // Triangle is non-singular. It is still possible that we won't need
        // all 3 points to be on the sphere surface. If one or more of the
        // barycentric coordinates is negative, it means that the circumsphere
        // center is outside the triangle. Pick the most negative, then drop
        // the opposite vertex. The circumsphere around the remaining edge will 
        // include the dropped one for free. 

        // s controls P's height over ac, t over ab, u=1-s-t over bc.
        const RealP ds2 = ab2ac2 - ac2*abac;                 // 2 flops
        const RealP dt2 = ab2ac2 - ab2*abac;                 // 2 flops
        const RealP du2 = dm2-ds2-dt2;                       // 2 flops
        const RealP oodm2 = 1/dm2; // (> 0)                   ~10 flops
        const RealP s = ds2*oodm2, t = dt2*oodm2, u = du2*oodm2; //  3 flops

        RealP minBary; int minBaryIx;
        minOf(s,t,u,minBary,minBaryIx); // 2 flops

        if (minBary <= 0) { // 1 flop
            Sphere_<P> edgeSphere; int map[2];
            switch(minBaryIx) {
            case 0: // s is the most negative
                map[0]=0; map[1]=2; // sphere around ac includes b
                edgeSphere=calcMinimumSphere(a,c,which); 
                assert(!edgeSphere.isPointOutside(b));
                break;
            case 1: // t is the most negative
                map[0]=0; map[1]=1; // sphere around ab includes c
                edgeSphere=calcMinimumSphere(a,b,which); 
                assert(!edgeSphere.isPointOutside(c));
                break;
            case 2: // u=1-s-t is the most negative
                map[0]=1; map[1]=2; // sphere around bc includes a
                edgeSphere=calcMinimumSphere(b,c,which); 
                assert(!edgeSphere.isPointOutside(a));
                break;
            };
            fixWhich(map, which);
            ctr = edgeSphere.getCenter();
        } else { // s,t,u > 0
            // All barycentric coordinates are positive. The circumsphere's 
            // center will be inside the triangle and thus of minimal size.
            which.clear(); 
            which.push_back(0); which.push_back(1); which.push_back(2);
            const Vec3P aToCenter = s*ab + t*ac;    //   9 flops
            ctr = a + aToCenter;                    //   3 flops
        }
    }

    // All methods lead here. At this point we have the support points
    // and center point chosen but still need to pick the radius.

    // This cleans up ugly roundoff errors that cause the center not to be
    // exactly equidistant from the points. This makes sure that the radius
    // touches the outermost point with the rest inside. This is 
    // expensive but worth it to ensure a trouble-free sphere.
    RealP rmax2; int rmaxIx;
    maxOf((a-ctr).normSqr(),(b-ctr).normSqr(),(c-ctr).normSqr(),
            rmax2, rmaxIx);                          // 27 flops
    const RealP rad = std::sqrt(rmax2);              // 20 flops

    return Sphere_<P>(ctr, rad);
}



//==============================================================================
//                       CALC BOUNDING SPHERE - 4 POINTS
//==============================================================================
/* This is a 4 point bounding sphere calculator based on Crister Ericson's
suggestion on his blog page -- see the 3 point routine above. The derivation
here was done by Sherm since Ericson didn't work out the 4-point case.

Cost is about 190 flops in a typical case where all 4 points are used.

Implementation
--------------
We have a tetrahedron with vertices a,b,c,d with no restrictions on their 
placement (i.e., they can be singular in a bunch of different ways). Algorithm:
    - Calculate the volume of the tetrahedron.
    - If the volume is very small:
         - Find the face of largest area from abc, abd, acd, bcd.
         - Drop the vertex that is not used in the largest face.
         - Calculate the bounding sphere of the largest face.
         - If the dropped vertex is outside, stretch the sphere to include it.
         - Return that sphere as the result (might not be perfect but will
           be very good).
    - Else (volume is OK):
         - Calculate barycentric coordinates s,t,u, v=1-s-t-u of the 
           circumsphere center P (the point equidistant to all vertices).
         - If s,t,u, or v < 0, P lies outside the tetrahedron and a 3-point 
           sphere will be smaller. Drop the corresponding vertex and return the
           bounding sphere for the remaining face. (3 support points).
         - Otherwise return the circumsphere (4 support points).

How to calculate the circumsphere:
    We express center P as P=a + s*(b-a) + t*(c-a) + u*(d-a). This point is to 
    be equidistant from all vertices, giving us 3 equations for s, t, and u:
           (P-b)^2 = (P-a)^2, (P-c)^2 = (P-a)^2, (P-d)^2 = (P-a)^2
    Substituting the definition of P gives us the following system of equations
    for s, t, and u:
         [ab^2  abac  abad] [s]   [ab^2/2]
         [abac  ac^2  acad] [t] = [ac^2/2]   or   M x = b
         [abad  acad  ad^2] [u] = [ad^2/2]
    We'll use Cramer's rule to solve the equation, which requires computing 
    determinants of four matrices: M, Ms, Mt, Mu where the latter are M with 
    its 1st, 2nd, or 3rd column replaced by b. Let dm=det(M), ds=det(Ms), 
    dt=det(Mt), du=det(Mu). Then s=ds/dm, t=dt/dm, u=du/dm. There are going to 
    be lots of common subexpressions, and we can pick up face areas along the 
    way. Also, determinant dm is related to the volume of the tetrahedron by 
    6*v = sqrt(dm), so dm = 36 v^2. Note that dm can't be negative except
    possibly due to roundoff error. Below we'll work with 2*dm, 2*ds, etc.
    because it saves some multiplies.
*/
template <class P> /*static*/
Geo::Sphere_<P> Geo::Sphere_<P>::
calcMinimumSphere(const Vec3P& a, const Vec3P& b, 
                  const Vec3P& c, const Vec3P& d,
                  Array_<int>& which) {
    const RealP tol = Geo::getDefaultTol<P>();

    const Vec3P ab = b-a, ac = c-a, ad = d-a;               //  9 flops
    const RealP abac = dot(ab, ac); // == |ab||ac|cos theta   ( 5 flops)
    const RealP abad = dot(ab, ad); // == |ab||ad|cos theta   ( 5 flops)
    const RealP acad = dot(ac, ad); // == |ac||ad|cos theta   ( 5 flops)
    // 25 flops in the next two lines.
    const RealP ab2=ab.normSqr(), ac2=ac.normSqr(), ad2=ad.normSqr(); 
    const RealP abac2=square(abac), abad2=square(abad), acad2=square(acad);
    const RealP ab2ac2=ab2*ac2, ab2ad2=ab2*ad2, ac2ad2=ac2*ad2;
    const RealP ab2ac2ad2=ab2ac2*ad2;
   
    const RealP dm2 = 2*(  ab2ac2ad2          // 10 flops
                         - ab2*acad2
                         - ad2*abac2
                         - ac2*abad2
                         + 2*abac*abad*acad); // 72 v^2

    // We'll use one of three ways to find the center point and determine
    // the set of supporting vertices (in which):
    //  1) singular (3, 2, or 1 support vertices)
    //  2) non-singular with 3 support vertices
    //  3) non-singular using circumsphere (4 support vertices)

    Vec3P ctr;

    // Consider the tetrahedron singular if its volume is less than tol.
    if (dm2 <= 72*square(tol)) {                             // 3 flops
        // Tetrahedron is near singular or very small.
        const Vec3P bc = c-b, cd = d-c;
        const RealP bccd = dot(bc, cd), bccd2 = square(bccd);
        const RealP bc2=bc.normSqr(), cd2=cd.normSqr(), bc2cd2=bc2*cd2;
        // Find the largest face.The expressions below are equivalent to, for 
        // example: (|ab||ac|)^2 * (1-cos^2 theta)
        // where theta is the angle between side ab and ac. Note that the area
        // of this triangle is area = 1/2 sqrt(ab2*ac2-square(abac)).
        const RealP abc4Area2 = ab2ac2-abac2; // 4*area(abc)^2
        const RealP abd4Area2 = ab2ad2-abad2; // 4*area(abd)^2
        const RealP acd4Area2 = ac2ad2-acad2; // 4*area(acd)^2
        const RealP bcd4Area2 = bc2cd2-bccd2; // 4*area(bcd)^2
        RealP maxArea2; int maxFace=0;
        maxOf(abc4Area2,abd4Area2,acd4Area2,bcd4Area2,maxArea2,maxFace);
        // Now create a sphere around that edge using 3,2, or 1 support points,
        // but we're only going to keep the center.
        Sphere_<P> faceSphere; int map[3];
        switch(maxFace) {
        case 0: map[0]=0; map[1]=1; map[2]=2; map[3]=3; //abc, drop d
            faceSphere=calcMinimumSphere(a,b,c,which); break;
        case 1: map[0]=0; map[1]=1; map[2]=3; map[3]=2; //abd, drop c
            faceSphere=calcMinimumSphere(a,b,d,which); break;
        case 2: map[0]=0; map[1]=2; map[2]=3; map[3]=1; //acd, drop b
            faceSphere=calcMinimumSphere(a,c,d,which); break;
        case 3: map[0]=1; map[1]=2; map[2]=3; map[3]=0; //bcd, drop a
            faceSphere=calcMinimumSphere(b,c,d,which); break;
        };
        fixWhich(map, which); // fix the indices
        ctr = faceSphere.getCenter();
    } else {
        // Tetrahedron is non-singular. It is still possible that we won't need
        // all four points to be on the sphere surface. If one or more of the
        // barycentric coordinates is negative, it means that the circumsphere
        // center is outside the tetrahedron. Pick the most negative, then drop
        // the opposite vertex. The circumsphere around the remaining face will 
        // include the dropped one for free. 

        // s controls height over acd, t over abd, u over abc, 1-s-t-u over bcd.
        const RealP ds2 =  ab2ac2ad2            // 12 flops
                         + ac2*abad*acad
                         + ad2*abac*acad
                         - ab2*acad2
                         - ac2ad2*abac
                         - ac2ad2*abad;
        const RealP dt2 =  ab2ac2ad2            // 12 flops
                         + ab2*abad*acad
                         + ad2*abac*abad
                         - ac2*abad2
                         - ab2ad2*acad
                         - ab2ad2*abac;
        const RealP du2 =  ab2ac2ad2            // 12 flops
                         + ab2*abac*acad
                         + ac2*abac*abad
                         - ad2*abac2
                         - ab2ac2*acad
                         - ab2ac2*abad;    
        const RealP dv2 = dm2-ds2-dt2-du2;      //  3 flops  
        const RealP oodm2 = 1/dm2; // (> 0)       ~10 flops
        const RealP s=ds2*oodm2, t=dt2*oodm2, u=du2*oodm2, v=dv2*oodm2; 
                                                //  4 flops

        RealP minBary; int minBaryIx;
        minOf(s,t,u,v,minBary,minBaryIx);   // 3 flops

        if (minBary <= 0) { // 1 flop
            Sphere_<P> faceSphere; int map[3];
            switch(minBaryIx) {
            case 0: // s is the most negative
                map[0]=0; map[1]=2; map[2]=3; // sphere around acd includes b
                faceSphere=calcMinimumSphere(a,c,d,which); 
                assert(!faceSphere.isPointOutside(b));
                break;
            case 1: // t is the most negative
                map[0]=0; map[1]=1; map[2]=3; // sphere around abd includes c
                faceSphere=calcMinimumSphere(a,b,d,which); 
                assert(!faceSphere.isPointOutside(c));
                break;
            case 2: // u is the most negative
                map[0]=0; map[1]=1; map[2]=2; // sphere around abc includes d
                faceSphere=calcMinimumSphere(a,b,c,which); 
                assert(!faceSphere.isPointOutside(d));
                break;
            case 3: // v is the most negative
                map[0]=1; map[1]=2; map[2]=3; // sphere around bcd includes a
                faceSphere=calcMinimumSphere(b,c,d,which); 
                assert(!faceSphere.isPointOutside(a));
                break;
            };
            fixWhich(map, which);
            ctr = faceSphere.getCenter();
        } else { // s,t,u,v > 0
            // All barycentric coordinates are positive. The circumsphere's 
            // center will be inside the tetrahedron and thus of minimal size.
            which.clear(); 
            which.push_back(0); which.push_back(1);  
            which.push_back(2); which.push_back(3);
            const Vec3P aToCenter = s*ab + t*ac + u*ad;     //  12 flops
            ctr = a + aToCenter;                            //   3 flops
        }
    }

    // All methods lead here. At this point we have the support points
    // and center point chosen but still need to pick the radius.

    // This cleans up ugly roundoff errors that cause the center not to be
    // exactly equidistant from the points. This makes sure that the radius
    // touches the outermost point with the rest inside. This is 
    // expensive but worth it to ensure a trouble-free sphere.
    RealP rmax2; int rmaxIx;
    maxOf((a-ctr).normSqr(),(b-ctr).normSqr(),
          (c-ctr).normSqr(),(d-ctr).normSqr(),
          rmax2, rmaxIx);                            // 35 flops
    const RealP rad = std::sqrt(rmax2);              // 20 flops

    return Sphere_<P>(ctr, rad);
}



//==============================================================================
//                       CALC BOUNDING SPHERE - N POINTS
//==============================================================================
// This is called recursively to calculate the minimum bounding sphere for a
// set of points (like the vertices of a mesh). It uses an algorithm developed 
// by Emo Welzl which, despite appearances, has O(n) expected running time.
// The implementation here is based on a description by Nicolas Capens at 
// http://www.flipcode.com/archives/Smallest_Enclosing_Spheres.shtml.
// As described there, the algorithm is highly susceptible to numerical 
// instabilities. Bernd Gartner describes an improved version in "Fast and 
// robust smallest enclosing balls", Proc. 7th Annual ACM European Symposium on 
// Algorithms, v. 1643 Lecture Notes in Computer Science, pp. 325-338, 1999.
// The implementation here follows Capens but the primitives have been 
// reworked so that they deal nicely with singularities and roundoff.

template <class P> static
Geo::Sphere_<P>
findWelzlSphere(const Array_<const Vec<3,P>*>& p, Array_<int>& ix,
                int bIn, Array_<int>& which, int recursionLevel) {

    // The first bIn points should be an independent support set, and we're 
    // hoping to calculate their circumsphere here. Although in theory the 
    // algorithm wouldn't have gotten here if all bIn points weren't 
    // independent, roundoff issues may cause the underlying primitive to 
    // disagree. So we might use only a subset of the available points, and
    // "which" will list the ones we used. If necessary we'll rearrange the
    // first few entries in "ix" to make sure that the indices of the new
    // support set come first.
    
    Geo::Sphere_<P> minSphere;

    switch(bIn) {
    // Create a bounding sphere for 0, 1, 2, 3, or 4 points. Note which
    // points were actually used, and how many. If we manage to use the
    // maximum of 4 support points, we return; otherwise, we'll fall through
    // and hunt for another support point.
    case 0: 
        minSphere = Geo::Sphere_<P>(Vec<3,P>(0),0);
        which.clear();
        break;
    case 1:
        minSphere = Geo::Sphere_<P>::calcMinimumSphere(*p[ix[0]],which);
        break;
    case 2:
        minSphere = Geo::Sphere_<P>::calcMinimumSphere
                            (*p[ix[0]],*p[ix[1]],which);
        break;
    case 3:
        minSphere = Geo::Sphere_<P>::calcMinimumSphere
                            (*p[ix[0]],*p[ix[1]],*p[ix[2]],which);
        break;
    case 4:
        // Never need more than 4 points.
        minSphere = Geo::Sphere_<P>::calcMinimumSphere
                            (*p[ix[0]],*p[ix[1]],*p[ix[2]],*p[ix[3]],which);
        if (which.size() == 4) {
            // We can return now after fixing up the indices in which.
            fixWhich(ix.begin(), which);
            return minSphere;
        }
        break;
    }

    // We're here because we used fewer than 4 support points. 
    
    // First, fix the indices in which (the primitives number from 0).
    fixWhich(ix.begin(), which);
    const int bActual = (int)which.size(); // <= 3

    // It is possible that we didn't use all the bIn points we were given, but
    // instead used a smaller number, bActual (probably bIn-1). In that case we
    // have to make sure that the *first* bActual are the support points. So 
    // we'll look at the entries [bActual..bIn-1]; if one is now part of the 
    // support set we'll swap it with a now-unused point in [0..bActual-1].
    for (int i=bActual; i < bIn; ++i) {
        const int* w = std::find(which.begin(), which.end(), ix[i]);
        if (w == which.end()) continue; // not being used
        // ix[i] is part of the support set
        bool swapped=false; // There has to be one to swap with!
        for (int j=0; j < bActual; ++j) {
            const int* ww = std::find(which.begin(), which.end(), ix[j]);
            if (ww == which.end()) {
                // ix[j] is not part of the support set; swap with i
                std::swap(ix[i], ix[j]); swapped=true;
                break;
            }
        }
        assert(swapped); // can't happen
    }

    // The indices of the support points are the first bActual entries in ix,
    // and there may be unused points after that that we already processed
    // but didn't need. Now run through all subsequent points and update the
    // sphere to include them.
  
    for (int i = bIn; i < (int)ix.size(); ++i) {
        if (minSphere.isPointOutside(*p[ix[i]])) {
            // This point is outside the current bounding sphere.  
            // Move it to the start of the list. (Without reordering; I *think*
            // that is necessary to avoid messing up other recursions, but I'm
            // not sure -- Sherm 111216.)
            for (int j = i; j > 0; --j)
                std::swap(ix[j], ix[j-1]);
            
            // Update the bounding sphere, taking the new point into account.
            ArrayView_<int> toBoundIx(ix.begin(), &ix[i]+1);
            minSphere = findWelzlSphere<P>(p, toBoundIx, bActual+1, which,
                                           recursionLevel+1);
        }
    }

    return minSphere;
}

// This signature takes an array of points, creates an array of pointers to 
// those points and calls the other signature.
template <class P> /*static*/
Geo::Sphere_<P> Geo::Sphere_<P>::
calcMinimumSphere(const Array_<Vec3P>& points, Array_<int>& which) {
    Array_<const Vec3P*> p(points.size());
    for (unsigned i=0; i<points.size(); ++i)
        p[i] = &points[i];
    return calcMinimumSphere(p, which);
}

template <class P> /*static*/
Geo::Sphere_<P> Geo::Sphere_<P>::
calcMinimumSphere(const Array_<const Vec3P*>& points, Array_<int>& which) {
    const unsigned npoints = points.size();

    // Allocate and initialize an array of point indices. These will get 
    // moved around during the computation.
    Array_<int> ix(npoints); for (unsigned i=0; i<npoints; ++i) ix[i] = i; 

    if (npoints < 10) {
        // Not worth rearranging.
        return findWelzlSphere<P>(points, ix, 1, which, 0);
    }

    // There are enough points that we'll try to improve the ordering so that
    // the bounding sphere gets large quickly. This optimization helps *a lot* 
    // for large numbers of points.

    // Find the six points that have the most extreme
    // x,y, and z coordinates (not necessarily six unique points) and move
    // them to the front so they get processed first.
    Vec3P lo=*points[0], hi=*points[0]; // initialize extremes
    int   ilo[3], ihi[3]; for (int i=0; i<3; ++i) ilo[i]=ihi[i]=0;
    for (unsigned i=0; i<points.size(); ++i) {
        const Vec3P& p = *points[i];
        if (p[0] > hi[0]) hi[0]=p[0], ihi[0]=i;
        if (p[0] < lo[0]) lo[0]=p[0], ilo[0]=i;
        if (p[1] > hi[1]) hi[1]=p[1], ihi[1]=i;
        if (p[1] < lo[1]) lo[1]=p[1], ilo[1]=i;
        if (p[2] > hi[2]) hi[2]=p[2], ihi[2]=i;
        if (p[2] < lo[2]) lo[2]=p[2], ilo[2]=i;
    }
    // Find the nx <= 6 unique extreme points.
    std::set<int> pending;
    pending.insert(ilo, ilo+3); pending.insert(ihi, ihi+3);
    const int nx = pending.size();
    // Go through the first n points. If the point is already an extreme,
    // remove it from the pending list. If not, swap it with one of the
    // extreme points.
    for (int i=0; i < nx; ++i) {
        std::set<int>::iterator p = pending.find(ix[i]);
        if (p != pending.end()) {
            pending.erase(p);
            continue;
        }
        p = pending.begin(); // first unmoved extreme
        const int extremeIx = *p;
        pending.erase(p);
        std::swap(ix[i], ix[extremeIx]);
    }

    return findWelzlSphere<P>(points, ix, 1, which, 0);
}


// Explicit instantiations for float and double.
template class Geo::Sphere_<float>;
template class Geo::Sphere_<double>;


}  // End of namespace SimTK