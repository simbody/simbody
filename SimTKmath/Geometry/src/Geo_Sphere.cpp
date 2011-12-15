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
//                       CALC BOUNDING SPHERE - 3 POINTS
//==============================================================================
/* (Note that the 3- and 4-point methods from Nicolas Capens compute only the
sphere that passes through all points, which can be *way* too big.)

This method produces the minimum bounding sphere around three points (i.e.,
a triangle) for all cases except some very, very small triangles. It is 
modified from Christer Ericson's blog: "Minimum bounding circle (sphere) for 
a triangle (tetrahedron)" July 27, 2007. 
http://realtimecollisiondetection.net/blog/?p=20
Accessed 12/12/2011 and reimplemented by Sherm. The main change is handling of
singular triangles.

Cost is about 90 flops for a typical case.

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
           [ab^2  abac] [s]   [ab^2]
       2 * [abac  ac^2] [t] = [ac^2]   or   M x = b
    We'll use Cramer's rule to solve the equation, which requires computing 
    determinants of three matrices: M, Ms, Mt where the latter are M with its 
    1st or 2nd column replaced by b. Let dm=det(M), ds=det(Ms), dt=det(Mt). 
    Then s=ds/dm, t=dt/dm. Determinant dm is related to the area of the 
    triangle by 2*a = sqrt(dm/2), so dm = 8 a^2. Note that dm can't be negative
    except possibly due to roundoff error.
*/
template <class P> /*static*/
Geo::Sphere_<P> Geo::Sphere_<P>::
calcMinimumSphere(const Vec3P& a, const Vec3P& b, const Vec3P& c,
                  Array_<int>& which) {
    const RealP tol = Geo::getDefaultTol<P>();
    const Vec3P* pts[3] = {&a,&b,&c}; // to support indexing
    int map[3]; // to support rearranging

    const Vec3P ab = b - a, ac = c - a;                     //  6 flops
    const RealP ab2 = ab.normSqr(), ac2 = ac.normSqr();     // 10 flops
    const RealP abac = dot(ab, ac); // == |ab||ac|cos theta   ( 5 flops)
    const RealP ab2ac2 = ab2*ac2; // 1 flop

    // The expression below is equivalent to 
    //      dm = 2 * (|ab||ac|)^2 * (1-cos^2 theta)
    // where theta is the angle between side ab and ac. That quantity can't
    // be negative except for roundoff error. Note that the area of this 
    // triangle is area = 1/2 sqrt(ab2*ac2-square(abac)) so dm = 8*area^2
    const RealP dm = 2*(ab2ac2 - square(abac));              // 3 flops

    // Consider the triangle singular if its area is less than tol.
    if (dm <= 8*square(tol)) {                               // 3 flops
        // Triangle is near singular or very small.
        const Vec3P bc = c - b;
        const RealP bc2 = bc.normSqr();
        // Find the longest edge.
        RealP maxLen2=ab2; int maxEdge=0;
        if (ac2 > maxLen2) maxLen2=ac2, maxEdge=1;
        if (bc2 > maxLen2) maxLen2=bc2, maxEdge=2;
        // Fill map[0..1] with longest edge; map[2] is the dropped vertex.
        switch(maxEdge) {
        case 0: map[0]=0; map[1]=1; map[2]=2; break; // use (a,b) drop c
        case 1: map[0]=0; map[1]=2; map[2]=1; break; // use (a,c) drop b
        case 2: map[0]=1; map[1]=2; map[2]=0; break; // use (b,c) drop a
        };
        const Vec3P& p0=*pts[map[0]]; const Vec3P& p1=*pts[map[1]];
        const Vec3P& p2=*pts[map[2]]; // this is the dropped one
        Sphere_<P> edgeSphere=calcMinimumSphere(p0,p1,which);
        for (unsigned i=0; i<which.size(); ++i)
            which[i] = map[which[i]]; // fix indices
        const Vec3P& ctr = edgeSphere.getCenter();
        RealP        rad = edgeSphere.getRadius();
        // If the point we dropped isn't inside, it can't be far from one
        // of the vertices we used. We'll just replace that vertex.
        const RealP  r2  = (p2-ctr).normSqr();
        if (which.size()==2 && r2 > square(rad)) {
            rad=std::sqrt(r2); // grow sphere & replace nearer support point
            const RealP d02 = Point_<P>::findDistanceSqr(*pts[which[0]], p2);
            const RealP d12 = Point_<P>::findDistanceSqr(*pts[which[1]], p2);
            if (d02 <= d12) which[0]=map[2]; // which[0] was closer
            else which[1]=map[2];            // which[1] was closer
        }
        // We used 2 or fewer support vertices.
        return Sphere_<P>(ctr, rad);
    } 

    // Triangle is non-singular. It is still very likely that we won't need
    // all three points to be on the sphere surface.

    // s controls P's height over ac, t over ab, 1-s-t over bc.
    const RealP ds = ab2ac2 - ac2*abac;     // 2 flops
    if (ds <= 0) {      // implies s<=0       (1 flop)
        map[0]=0; map[1]=2; // sphere around ac includes b
        const Sphere_<P> sph=calcMinimumSphere(a,c,which); 
        for (unsigned i=0; i<which.size(); ++i) which[i]=map[which[i]];
        return sph;
    }

    const RealP dt = ab2ac2 - ab2*abac;     // 2 flops
    if (dt <= 0) {      // implies t<=0       (1 flop)
        map[0]=0; map[1]=1; // sphere around ab includes c
        const Sphere_<P> sph=calcMinimumSphere(a,b,which); 
        for (unsigned i=0; i<which.size(); ++i) which[i]=map[which[i]];
        return sph;
    }

    if (ds+dt >= dm) {  // implies u=1-s-t<=0 (2 flops)
        map[0]=1; map[1]=2; // sphere around bc includes a
        const Sphere_<P> sph=calcMinimumSphere(b,c,which); 
        for (unsigned i=0; i<which.size(); ++i) which[i]=map[which[i]];
        return sph;
    }

    // Must calculate circumsphere and use all the points.
    const RealP oodm = 1/dm; // (> 0)        ~10 flops
    const RealP s = ds*oodm, t = dt*oodm; //   2 flops
    const Vec3P aToCenter = s*ab + t*ac;  //   9 flops
    const Vec3P ctr = a + aToCenter;      //   3 flops
    const RealP rad = aToCenter.norm();   // ~25 flops

    // We're using all the points!
    which.clear(); which.push_back(0); which.push_back(1); which.push_back(2);
    return Sphere_<P>(ctr, rad);    // 3 flops
}



//==============================================================================
//                       CALC BOUNDING SPHERE - 4 POINTS
//==============================================================================
/* This is a 4 point bounding sphere calculator based on Crister Ericson's
suggestion on his blog page -- see the 3 point routine above. The derivation
here was done by Sherm since Ericson didn't work out the 4-point case.

Cost is about 160 flops in a typical case.

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
         [ab^2  abac  abad] [s]   [ab^2]
       2*[abac  ac^2  acad] [t] = [ac^2]   or   M x = b
         [abad  acad  ad^2] [u] = [ad^2]
    We'll use Cramer's rule to solve the equation, which requires computing 
    determinants of four matrices: M, Ms, Mt, Mu where the latter are M with 
    its 1st, 2nd, or 3rd column replaced by b. Let dm=det(M), ds=det(Ms), 
    dt=det(Mt), du=det(Mu). Then s=ds/dm, t=dt/dm, u=du/dm. There are going to 
    be lots of common subexpressions, and we can pick up face areas along the 
    way. Also, determinant dm is related to the volume of the tetrahedron by 
    6*v = sqrt(dm/2), so dm = 72 v^2. Note that dm can't be negative except
    possibly due to roundoff error.
*/
template <class P> /*static*/
Geo::Sphere_<P> Geo::Sphere_<P>::
calcMinimumSphere(const Vec3P& a, const Vec3P& b, 
                  const Vec3P& c, const Vec3P& d,
                  Array_<int>& which) {
    const RealP tol = Geo::getDefaultTol<P>();
    const Vec3P* pts[4] = {&a,&b,&c,&d}; // to support indexing
    int map[4]; // to support rearranging

    const Vec3P ab = b-a, ac = c-a, ad = d-a;               //  9 flops
    const RealP abac = dot(ab, ac); // == |ab||ac|cos theta   ( 5 flops)
    const RealP abad = dot(ab, ad); // == |ab||ad|cos theta   ( 5 flops)
    const RealP acad = dot(ac, ad); // == |ac||ad|cos theta   ( 5 flops)
    // 25 flops in the next two lines.
    const RealP ab2=ab.normSqr(), ac2=ac.normSqr(), ad2=ad.normSqr(); 
    const RealP abac2=square(abac), abad2=square(abad), acad2=square(acad);
    const RealP ab2ac2=ab2*ac2, ab2ad2=ab2*ad2, ac2ad2=ac2*ad2;
    const RealP ab2ac2ad2=ab2ac2*ad2;
   
    const RealP dm = 2*(  ab2ac2ad2             // 10 flops
                         - ab2*acad2
                         - ad2*abac2
                         - ac2*abad2
                         + 2*abac*abad*acad); // 72 v^2

    // Consider the tetrahedron singular if its volume is less than tol.
    if (dm <= 72*square(tol)) {                             // 3 flops
        // Tetrahedron is near singular or very small.
        const Vec3P bc = c-b, cd = d-c;
        const RealP bccd = dot(bc, cd), bccd2 = square(bccd);
        const RealP bc2=bc.normSqr(), cd2=cd.normSqr(), bc2cd2=bc2*cd2;
        // The expressions below are equivalent to, for example 
        //      (|ab||ac|)^2 * (1-cos^2 theta)
        // where theta is the angle between side ab and ac. Note that the area
        // of this triangle is area = 1/2 sqrt(ab2*ac2-square(abac)).
        const RealP abc4Area2 = ab2ac2-abac2; // 4*area(abc)^2
        const RealP abd4Area2 = ab2ad2-abad2; // 4*area(abd)^2
        const RealP acd4Area2 = ac2ad2-acad2; // 4*area(acd)^2
        const RealP bcd4Area2 = bc2cd2-bccd2; // 4*area(bcd)^2
        // Find the largest face.
        RealP maxArea=abc4Area2;  int   maxFace=0;
        if (abd4Area2 > maxArea) maxArea=abd4Area2, maxFace=1;
        if (acd4Area2 > maxArea) maxArea=acd4Area2, maxFace=2;
        if (bcd4Area2 > maxArea) maxArea=bcd4Area2, maxFace=3;
        // Fill map[0..2] with largest face; map[3] is the dropped vertex.
        switch(maxFace) {
        case 0: map[0]=0; map[1]=1; map[2]=2; map[3]=3; break; //abc, drop d
        case 1: map[0]=0; map[1]=1; map[2]=3; map[3]=2; break; //abd, drop c
        case 2: map[0]=0; map[1]=2; map[2]=3; map[3]=1; break; //acd, drop b
        case 3: map[0]=1; map[1]=2; map[2]=3; map[3]=0; break; //bcd, drop a
        };
        const Vec3P& p0=*pts[map[0]]; const Vec3P& p1=*pts[map[1]];
        const Vec3P& p2=*pts[map[2]];
        const Vec3P& p3=*pts[map[3]]; // this is the dropped one
        Sphere_<P> faceSphere=calcMinimumSphere(p0,p1,p2,which);
        for (unsigned i=0; i<which.size(); ++i)
            which[i] = map[which[i]]; // fix indices
        const Vec3P& ctr = faceSphere.getCenter();
        RealP        rad = faceSphere.getRadius();
        // If the point we dropped isn't inside, it can't be far from one
        // of the vertices we used. We'll just replace that vertex.
        const RealP  r2  = (p3-ctr).normSqr();
        if (which.size()==3 && r2 > square(rad)) {
            rad=std::sqrt(r2); // grow sphere & replace nearest support point
            const RealP d03 = Point_<P>::findDistanceSqr(*pts[which[0]], p3);
            const RealP d13 = Point_<P>::findDistanceSqr(*pts[which[1]], p3);
            const RealP d23 = Point_<P>::findDistanceSqr(*pts[which[2]], p3);
            if (d03 <= d13 && d03 <= d23) which[0]=map[3]; // which[0] closest
            else if (d13 <= d03 && d13 <= d23) which[1]=map[3]; // [1] closest
            else which[2]=map[3];   // which[2] was closest
        }
        // We used 3 or fewer support vertices.
        return Sphere_<P>(ctr, rad);
    }

    // Tetrahedron is non-singular. It is still very likely that we won't need
    // all four points to be on the sphere surface.

    // s controls height over acd, t over abd, u over abc, 1-s-t-u over bcd.
    const RealP ds =   ab2ac2ad2            // 12 flops
                     + ac2*abad*acad
                     + ad2*abac*acad
                     - ab2*acad2
                     - ac2ad2*abac
                     - ac2ad2*abad;
    if (ds <= 0) { // implies s<=0 (1 flop)
        map[0]=0; map[1]=2; map[2]=3; // sphere around acd includes b
        const Sphere_<P> sph=calcMinimumSphere(a,c,d,which); 
        for (unsigned i=0; i<which.size(); ++i) which[i]=map[which[i]];
        return sph;
    }
    const RealP dt =   ab2ac2ad2            // 12 flops
                     + ab2*abad*acad
                     + ad2*abac*abad
                     - ac2*abad2
                     - ab2ad2*acad
                     - ab2ad2*abac;
    if (dt <= 0) { // implies t<=0 (1 flop)
        map[0]=0; map[1]=1; map[2]=3; // sphere around abd includes c
        const Sphere_<P> sph=calcMinimumSphere(a,b,d,which); 
        for (unsigned i=0; i<which.size(); ++i) which[i]=map[which[i]];
        return sph;
    }

    const RealP du =   ab2ac2ad2            // 12 flops
                     + ab2*abac*acad
                     + ac2*abac*abad
                     - ad2*abac2
                     - ab2ac2*acad
                     - ab2ac2*abad;
    if (du <= 0) { // implies u<=0 (1 flop)
        map[0]=0; map[1]=1; map[2]=2; // sphere around abc includes d
        const Sphere_<P> sph=calcMinimumSphere(a,b,c,which); 
        for (unsigned i=0; i<which.size(); ++i) which[i]=map[which[i]];
        return sph;
    }

    if (ds+dt+du >= dm) { // implies v=1-s-t-u<=0 (3 flops)
        map[0]=1; map[1]=2; map[2]=3; // sphere around bcd includes a
        const Sphere_<P> sph=calcMinimumSphere(b,c,d,which); 
        for (unsigned i=0; i<which.size(); ++i) which[i]=map[which[i]];
        return sph;
    }

    // Must calculate circumsphere.
    const RealP oodm = 1/dm; // (> 0)                     ~10 flops
    const RealP s = ds*oodm, t = dt*oodm, u = du*oodm; //   3 flops
    const Vec3P aToCenter = s*ab + t*ac + u*ad;        //  12 flops
    const Vec3P ctr = a + aToCenter;                   //   3 flops
    const RealP rad = aToCenter.norm();                // ~25 flops

    // We're using all the points!
    which.clear(); which.push_back(0); which.push_back(1); 
    which.push_back(2); which.push_back(3);
    return Sphere_<P>(ctr, rad);    // 3 flops
}



//==============================================================================
//                       CALC BOUNDING SPHERE - N POINTS
//==============================================================================
template <class P> /*static*/
Geo::Sphere_<P> Geo::Sphere_<P>::
calcMinimumSphere(const Array_<Vec3P>& points, Array_<int>& which) {
    switch(points.size()) {
    case 0: 
        which.clear();
        return Sphere_(Vec3P(0),0);
    case 1:
        return calcMinimumSphere(points[0],which);
    case 2:
        return calcMinimumSphere(points[0],points[1],which);
    case 3:
        return calcMinimumSphere(points[0],points[1],points[2],which);
    case 4:
        return calcMinimumSphere(points[0],points[1],points[2],points[3],which);
    default: 
        SimTK_ASSERT_ALWAYS(!"implemented",
            "Sphere_<P>::calcMinimumSphere(): only up to 4 so far.");
    }
    return Sphere_();
}



// Explicit instantiations for float and double.
template class Geo::Sphere_<float>;
template class Geo::Sphere_<double>;


}  // End of namespace SimTK