/* -------------------------------------------------------------------------- *
 *                        SimTK Simbody: SimTKmath                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-11 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
 * Contributors: Michael Sherman                                              *
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

#include "SimTKcommon.h"
#include "simmath/internal/common.h"
#include "simmath/internal/Geo.h"
#include "simmath/internal/Geo_Point.h"
#include "simmath/internal/Geo_Sphere.h"
#include "simmath/internal/BicubicSurface.h"
#include "simmath/internal/ContactGeometry.h"

#include "ContactGeometryImpl.h"

#include <cmath>
#include <map>
#include <set>

using namespace SimTK;
using std::map;
using std::pair;
using std::set;
using std::string;

//==============================================================================
//                            CONTACT GEOMETRY
//==============================================================================

ContactGeometry::ContactGeometry(ContactGeometryImpl* impl) : impl(impl) {
    assert(impl);
    impl->setMyHandle(*this);
}

ContactGeometry::~ContactGeometry() {
    if (isOwnerHandle())
        delete impl;
    impl = 0;
}

bool ContactGeometry::isOwnerHandle() const {
    return (impl == 0 || impl->getMyHandle() == this);
}

bool ContactGeometry::isEmptyHandle() const {
    return (impl == 0);
}

ContactGeometry::ContactGeometry(const ContactGeometry& src) : impl(0) {
    if (src.impl) {
        impl = src.impl->clone();
        impl->setMyHandle(*this);
    }
}

ContactGeometry& ContactGeometry::operator=(const ContactGeometry& src) {
    if (&src != this) {
        if (isOwnerHandle())
            delete impl;
        impl = 0;
        if (src.impl) {
            impl = src.impl->clone();
            impl->setMyHandle(*this);
        }
    }
    return *this;
}

ContactGeometryTypeId ContactGeometry::
getTypeId() const {return getImpl().getTypeId();}

bool ContactGeometry::intersectsRay(const Vec3& origin, const UnitVec3& direction, Real& distance, UnitVec3& normal) const {
    return getImpl().intersectsRay(origin, direction, distance, normal);
}

Vec3 ContactGeometry::findNearestPoint(const Vec3& position, bool& inside, UnitVec3& normal) const {
    return getImpl().findNearestPoint(position, inside, normal);
}

void ContactGeometry::getBoundingSphere(Vec3& center, Real& radius) const {
    getImpl().getBoundingSphere(center, radius);
}

bool ContactGeometry::isSmooth() const {return getImpl().isSmooth();}
bool ContactGeometry::isConvex() const {return getImpl().isConvex();}

void ContactGeometry::calcCurvature(const Vec3& point, Vec2& curvature, 
                                    Rotation& orientation) const 
{   getImpl().calcCurvature(point, curvature, orientation); }

const Function& ContactGeometry::getImplicitFunction() const 
{   return getImpl().getImplicitFunction(); }

Vec3 ContactGeometry::calcSupportPoint(UnitVec3 direction) const 
{   return getImpl().calcSupportPoint(direction); }

/*static*/Vec2 ContactGeometry::
evalParametricCurvature(const Vec3& P, const UnitVec3& nn,
                        const Vec3& dPdu, const Vec3& dPdv,
                        const Vec3& d2Pdu2, const Vec3& d2Pdv2, 
                        const Vec3& d2Pdudv,
                        Transform& X_EP)
{
    // All this is 42 flops
    Real E =  ~dPdu*dPdu,  F =  ~dPdu*dPdv,   G =  ~dPdv*dPdv;
    Real e =-(~d2Pdu2*nn), f =-(~d2Pdudv*nn), g =-(~d2Pdv2*nn);
    Real A = F*g-G*f, B = E*g-G*e, C = E*f-F*e;

    Real kmax, kmin;
    UnitVec3 dmax;
    if (std::abs(F) < SignificantReal) {
        Real ku = e/E, kv = g/G; // two divides ~20 flops
        if (ku < kv) {
            kmax=kv, kmin=ku;
            dmax=UnitVec3(dPdv); // normalizing, ~35 flops
        } else {
            kmax=ku, kmin=kv;
            dmax=UnitVec3(dPdu); // normalizing, ~35 flops
        }
    } else {
        // ~40 flops
        // t = (-b +/- sqrt(b^2-4ac)) / 2a
        // Discriminant must be nonnegative for real surfaces
        // but could be slightly negative due to numerical noise.
        Real sqrtd = std::sqrt(std::max(B*B - 4*A*C, Real(0)));
        Vec2 t = Vec2(sqrtd - B, -sqrtd - B) / (2*A);

        // Two divides + misc: ~30 flops
        Real kr = (e + f*t[0])/(E+F*t[0]); // Struik, eq. 6-4, pg 80
        Real ks = (e + f*t[1])/(E+F*t[1]); // (works only because these are extremes)
                                           // otherwise use eq. 6-3.

        if (kr < ks) {
            kmax=ks, kmin=kr;
            dmax = UnitVec3(t[1]*dPdv + dPdu); // Sdir, normalizing, ~40 flops
        } else {
            kmax=kr, kmin=ks;
            dmax = UnitVec3(t[0]*dPdv + dPdu); // Rdir, normalizing, ~40 flops
        }
    }

    // y=z%x ensures right handed; already unit vec (9 flops)
    UnitVec3 dmin = UnitVec3(nn % dmax, true);
    X_EP.updR().setRotationFromUnitVecsTrustMe(dmax, dmin, nn);
    X_EP.updP() = P; // the origin point

    return Vec2(kmax, kmin);
}

// See the documentation in the header file for a complete description of
// what's being calculated here. This comment adds implementation information
// that isn't relevant to the API user. 
// 
// Given two paraboloids P1 and P2 sharing a common normal z and origin,
// compute the paraboloid that represents their difference, and express that 
// paraboloid in a frame that has been rotated around z so that x and y 
// coincide with the principal curvature directions. P1 and P2 may be
// elliptic (kmax>=kmin>=0) or hyperbolic (kmax>=0>kmin). If the surfaces are 
// non-conforming, their difference will be elliptic with kmax>=kmin>0. 
//
// We assume the paraboloids represent surfaces and that each has its z axis 
// oriented away from the surface, pointing outside the "bowl" of the elliptic 
// paraboloid or away from the convex direction of a hyperbolic paraboloid. 
// That's the opposite sense from a standard paraboloid parameterization.
// The z axes are antiparallel. We will return the resulting difference 
// paraboloid in a frame whose z axis is coincident with P1's z axis, and thus 
// antiparallel to P2's z axis.
//
//     P1: z = -(kmax1/2 x1^2 + kmin1/2 y1^2)
//     P2: z =   kmax2/2 x2^2 + kmin2/2 y2^2
//      P: z = -( kmax/2  x^2 +  kmin/2  y^2)
// Thus the right-handed coordinate frames are:
//     P1: (x1,y1, z)
//     P2: (x2,y2,-z)
//      P: ( x, y, z)
// The above distinctions don't matter a whole lot for the implementation here,
// but still, I thought you might like to know anyway.
//
// Cost is about 70 flops to get the curvatures kmax,kmin. Then if you want
// the curvature directions too it costs another 150 flops. 

// This local static helper method calculates the curvatures and returns 
// intermediates necessary for calculating the directions, but doesn't actually
// calculate them. So we use only about 70 flops here.
static void combineParaboloidsHelper
   (const Rotation& R_SP1, const Vec2& k1,
    const UnitVec3& x2, const Vec2& k2,
    Real& cos2w, Real& sin2w, Real& kdiff1, Real& kdiff2, Vec2& k)
{
    const UnitVec3& x1 = R_SP1.x(); // P1 kmax direction
    const UnitVec3& y1 = R_SP1.y(); // P1 kmin direction
    const UnitVec3& z  = R_SP1.z(); // P1, P, -P2 normal

    const Real ksum1  = k1[0]+k1[1], ksum2  = k2[0]+k2[1]; // 4 flops
    kdiff1 = k1[0]-k1[1], kdiff2 = k2[0]-k2[1];

    // w is angle between x1, x2 max curvature directions defined
    // using right hand rule rotation of x1 about z until it is
    // coincident with x2. But ... we want -90 <= w <= 90, meaning
    // cos(w) >= 0. If necessary we flip x2 180 degrees around z, 
    // since -x2 is an equally good max curvature direction.
    const Real dotx1x2 = dot(x1,x2);         // 5 flops
    const UnitVec3 x2p = dotx1x2 < 0 ? -x2 : x2;
    const Real cosw = std::abs(dotx1x2);
    const Real sinw = dot(cross(x1,x2p), z); // signed, 14 flops

    // We'll need cos(2w), sin(2w); luckily these are easy to get (5 flops).
    cos2w = 2*square(cosw) - 1; // double angle formulas
    sin2w = 2*sinw*cosw;

    // Compute min/max curvatures of the difference surface.
    // See KL Johnson 1987 Ch. 4 and Appendix 2, and J-F Antoine, et al. 
    // 2006 pg 661. ~35 flops
    const Real ksum = ksum1 + ksum2;
    const Real kdiff = std::sqrt(square(kdiff1) + square(kdiff2)
                                 + 2*kdiff1*kdiff2*cos2w);
    k = Vec2(ksum + kdiff, ksum - kdiff)/2; // kmax, kmin (4 flops)  
}

// This is the full version that calculates the curvatures for 70 flops
// then spends another 150 to get the curvature directions.
/*static*/ void ContactGeometry::combineParaboloids
   (const Rotation& R_SP1, const Vec2& k1,
    const UnitVec3& x2, const Vec2& k2,
    Rotation& R_SP, Vec2& k)
{
    Real cos2w, sin2w, kdiff1, kdiff2;
    combineParaboloidsHelper(R_SP1, k1, x2, k2,
                             cos2w, sin2w, kdiff1, kdiff2, k);

    // Now find the rotated coordinate system by solving for the
    // angle -90 <= alpha <= 90 by which we need to rotate the x1 axis
    // about z to align it with the x axis. See KL Johnson Appendix 2
    // again, noting that beta = theta-alpha, then solving for tan(2 alpha).
    // This is about 130 flops.
    const Real yy = kdiff2*sin2w, xx = kdiff2*cos2w + kdiff1;
    Real a = std::atan2(yy,xx) / 2; // yy==xx==0 -> a=0
    Real cosa = std::cos(a), sina = std::sin(a);

    // Perform the actual rotations of x1,y1 to get x,y (18 flops)
    const UnitVec3& x1 = R_SP1.x(); // P1 kmax direction
    const UnitVec3& y1 = R_SP1.y(); // P1 kmin direction
    const UnitVec3& z  = R_SP1.z(); // P1, P, -P2 normal
    R_SP.setRotationFromUnitVecsTrustMe(UnitVec3(cosa*x1 + sina*y1, true),
                                        UnitVec3(cosa*y1 - sina*x1, true),
                                        z);
}

// This is the abridged version that costs only 70 flops but gives you just
// the curvatures without the curvature directions.
/*static*/ void ContactGeometry::combineParaboloids
   (const Rotation& R_SP1, const Vec2& k1,
    const UnitVec3& x2, const Vec2& k2,
    Vec2& k)
{
    Real cos2w, sin2w, kdiff1, kdiff2; // unneeded
    combineParaboloidsHelper(R_SP1, k1, x2, k2,
                             cos2w, sin2w, kdiff1, kdiff2, k);     
}



//==============================================================================
//                             HALF SPACE & IMPL
//==============================================================================
ContactGeometry::HalfSpace::HalfSpace()
:   ContactGeometry(new HalfSpace::Impl()) {}

/*static*/ ContactGeometryTypeId ContactGeometry::HalfSpace::classTypeId() 
{   return ContactGeometry::HalfSpace::Impl::classTypeId(); }

// Point position is given in the half space frame.
Vec3 ContactGeometry::HalfSpace::Impl::findNearestPoint
   (const Vec3& position, bool& inside, UnitVec3& normal) const {
    inside = (position[0] >= 0);
    normal = -UnitVec3(XAxis); // this does not require normalization
    return Vec3(0, position[1], position[2]);
}

bool ContactGeometry::HalfSpace::Impl::intersectsRay
   (const Vec3& origin, const UnitVec3& direction, 
    Real& distance, UnitVec3& normal) const 
{
    if (std::abs(direction[0]) < SignificantReal)
        return false; // ray is parallel to halfspace surface

    const Real t = origin[0]/direction[0];
    if (t > 0)
        return false; // ray points away from surface

    distance = -t;
    normal = -UnitVec3(XAxis); // cheap; no normalization required
    return true;
}

void ContactGeometry::HalfSpace::Impl::getBoundingSphere
   (Vec3& center, Real& radius) const 
{   center = Vec3(0);
    radius = Infinity; }

const ContactGeometry::HalfSpace::Impl& ContactGeometry::HalfSpace::
getImpl() const {
    assert(impl);
    return static_cast<const HalfSpace::Impl&>(*impl);
}

ContactGeometry::HalfSpace::Impl& ContactGeometry::HalfSpace::
updImpl() {
    assert(impl);
    return static_cast<HalfSpace::Impl&>(*impl);
}



//==============================================================================
//                               SPHERE & IMPL
//==============================================================================

ContactGeometry::Sphere::Sphere(Real radius) 
:   ContactGeometry(new Sphere::Impl(radius)) {}

/*static*/ ContactGeometryTypeId ContactGeometry::Sphere::classTypeId() 
{   return ContactGeometry::Sphere::Impl::classTypeId(); }

Real ContactGeometry::Sphere::getRadius() const {
    return getImpl().getRadius();
}

void ContactGeometry::Sphere::setRadius(Real radius) {
    updImpl().setRadius(radius);
}

const ContactGeometry::Sphere::Impl& ContactGeometry::Sphere::getImpl() const {
    assert(impl);
    return static_cast<const Sphere::Impl&>(*impl);
}

ContactGeometry::Sphere::Impl& ContactGeometry::Sphere::updImpl() {
    assert(impl);
    return static_cast<Sphere::Impl&>(*impl);
}

Vec3 ContactGeometry::Sphere::Impl::findNearestPoint(const Vec3& position, bool& inside, UnitVec3& normal) const {
    inside = (position.normSqr() <= radius*radius);
    normal = UnitVec3(position); // expensive -- normalizing
    return normal*radius;
}

bool ContactGeometry::Sphere::Impl::intersectsRay
   (const Vec3& origin, const UnitVec3& direction, 
    Real& distance, UnitVec3& normal) const 
{
    Real b = -~direction*origin;;
    Real c = origin.normSqr() - radius*radius;
    if (c > 0) {
        // Ray origin is outside sphere.

        if (b <= 0)
          return false;  // Ray points away from center of sphere.
        Real d = b*b - c;
        if (d < 0)
          return false;
        Real root = std::sqrt(d);
        distance = b - root;
      }
    else {
        // Ray origin is inside sphere.

        Real d = b*b - c;
        if (d < 0)
          return false;
        distance = b + std::sqrt(d);
      }
    normal = UnitVec3(origin+distance*direction);
    return true;
}

void ContactGeometry::Sphere::Impl::
getBoundingSphere(Vec3& center, Real& radius) const {
    center = Vec3(0);
    radius = this->radius;
}

void ContactGeometry::Sphere::Impl::
calcCurvature(const Vec3& point, Vec2& curvature, Rotation& orientation) const {
    orientation = Rotation(UnitVec3(point), ZAxis, fabs(point[0]) > 0.5 ? Vec3(0, 1, 0) : Vec3(1, 0, 0), XAxis);
    curvature = 1/radius;
}

Real SphereImplicitFunction::
calcValue(const Vector& x) const {
    return 1-(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])/square(ownerp->getRadius());
}

Real SphereImplicitFunction::
calcDerivative(const Array_<int>& derivComponents, const Vector& x) const {
    if (derivComponents.size() == 1)
        return 2*x[derivComponents[0]]/square(ownerp->getRadius());
    if (derivComponents[0] == derivComponents[1])
        return 2/square(ownerp->getRadius());
    return 0;
}



//==============================================================================
//                               ELLIPSOID & IMPL
//==============================================================================

ContactGeometry::Ellipsoid::Ellipsoid(const Vec3& radii)
:   ContactGeometry(new Ellipsoid::Impl(radii)) {}

void ContactGeometry::Ellipsoid::setRadii(const Vec3& radii) 
{   updImpl().setRadii(radii); }

/*static*/ ContactGeometryTypeId ContactGeometry::Ellipsoid::classTypeId()
{   return ContactGeometry::Ellipsoid::Impl::classTypeId(); }

const Vec3& ContactGeometry::Ellipsoid::getRadii() const 
{   return getImpl().getRadii(); }

const Vec3& ContactGeometry::Ellipsoid::getCurvatures() const 
{   return getImpl().getCurvatures(); }

UnitVec3 ContactGeometry::Ellipsoid::
findUnitNormalAtPoint(const Vec3& Q) const
{   return getImpl().findUnitNormalAtPoint(Q); }

Vec3 ContactGeometry::Ellipsoid::
findPointWithThisUnitNormal(const UnitVec3& nn) const
{   return getImpl().findPointWithThisUnitNormal(nn); }

Vec3 ContactGeometry::Ellipsoid::
findPointInSameDirection(const Vec3& Q) const
{   return getImpl().findPointInSameDirection(Q); }

void ContactGeometry::Ellipsoid::
findParaboloidAtPoint(const Vec3& Q, Transform& X_EP, Vec2& k) const
{   return getImpl().findParaboloidAtPoint(Q,X_EP,k); }

void ContactGeometry::Ellipsoid::
findParaboloidAtPointWithNormal(const Vec3& Q, const UnitVec3& nn,
                                Transform& X_EP, Vec2& k) const
{   return getImpl().findParaboloidAtPointWithNormal(Q,nn,X_EP,k); }


const ContactGeometry::Ellipsoid::Impl& ContactGeometry::Ellipsoid::
getImpl() const {
    assert(impl);
    return static_cast<const Ellipsoid::Impl&>(*impl);
}

ContactGeometry::Ellipsoid::Impl& ContactGeometry::Ellipsoid::
updImpl() {
    assert(impl);
    return static_cast<Ellipsoid::Impl&>(*impl);
}

// Given a point Q on an ellipsoid, with outward unit normal nn at Q: find the 
// principal curvatures at the point and their directions. The result is a 
// coordinate frame with origin Q, z axis the ellipsoid normal nn at Q, x axis 
// is the direction dmax of maximum curvature kmax, y axis the direction dmin 
// of minimum curvature kmin, such that [dmax dmin n] forms a right-handed set.
// This is equivalent to fitting an elliptic paraboloid 
// z = -kmax/2 x^2 -kmin/2 y^2 to the ellipsoid at point Q. Note that for
// an ellipsoid we have kmax>=kmin>0.
//
// We'll find the ellipse on the central plane perpendicular to the normal by 
// intersecting the plane equation with the ellipsoid equation but working in 
// the plane frame P=[u v n], where u and v are arbitrary axes in the plane.
// Our goal is to obtain an equation for the ellipse in P and then rotate the 
// P frame about its normal until we get the ellipse in standard form 
// Ru^2+Sv^2=1 in which case d/R and d/S are the ellipsoid curvatures (d is the
// distance from the point on the ellipsoid to the plane).
// ref: McArthur, Neil. "Principal radii of curvature at a point on an 
// ellipsoid", Mathematical Notes 24 pp. xvi-xvii, 1929.
//
// In its own frame E=[x y z] the ellipsoid surface is the set of points such 
// that
//    ~e * diag(A,B,C) * e = 1
// where e is a vector expressed in E. The plane is the set of points 
// satisfying ~e * n = 0. We can write rotation matrix R_EP=[u v n] where 
// u,v,n are expressed in E. Now we can put the ellipsoid in P:
//   ~(R_EP*p) * diag(A,B,C) * (R_EP*p) = 1
// We can intersect that with the plane just by dropping the n coordinate of 
// p so p=[u v 0] (u,v scalars here), and the intersection equation is
//    A(u*ux + v*vx)^2 + B(u*uy+v*vy)^2 + C(u*uz + v*vz)^2 = 1
// which is
//    R u^2 + S v^2 + T u*v = 1
// with
//    R =   A ux^2  + B uy^2  + C uz^2
//    S =   A vx^2  + B vy^2  + C vz^2
//    T = 2(A ux*vx + B uy*vy + C uz*vz)
//
// We want to find a rotation about n that eliminates the cross term Tuv, 
// leaving us with
//    R' u'^2 + S' v'^2 = 1
// for new constants R' and S' and new basis u' and v'.
//
// Method
// ------
// We'll calculate an angle theta where theta=0 would be along u and 
// theta=pi/2 would be along v. Then theta+pi/2 is a perpendicular direction 
// that has the other curvature extreme. Per "Dr Rob" at Mathforum.org 2000:
//   t2t = tan(2*theta) = T/(R-S)
//   theta = atan(t2t)/2, c = cos(theta), s = sin(theta)
//   R' = Rc^2 + Tsc + Ss^2   (theta direction)
//   S' = Rs^2 - Tsc + Sc^2   (theta+pi/2 direction)
// Directions are u' = c*u + s*v, v' = c*v - s*u; these are automatically unit
// vectors.
//
// Optimization
// ------------
// The above requires an atan() to get 2*theta then sin & cos(theta) at
// a cost of about 120 flops. We can use half angle formulas to work
// exclusively with 2*theta, but then we'll have to normalize u' and v' 
// at the end:
//   t2t = tan(2*theta) = T/(R-S)
//   c2t = cos(2*theta) = 1/sqrt(1 + t2t^2)
//   s2t = sin(2*theta) = t2t*cos2t;
//   2*R' = R+S + Rc2t - Sc2t + Ts2t
//   2*S' = R+S - Rc2t + Sc2t - Ts2t
// By multiplying the u',v' formulas above by 2*c we change the lengths
// but get expressions that are easily converted to double angles:
//   u' = normalize((1+c2t)*u + s2t*v)
//   v' = normalize((1+c2t)*v - s2t*u)
// (but actually v' is n X u' which is cheap). This saves about 30 
// flops over the straightforward method above.
//
// Cost: given a point and normalized normal
//    curvatures ~160 flops
//    directions ~ 60 flops more
//               ----
//               ~220 flops
//
// So: Given an ellipsoid in its own frame E, with equation Ax^2+By^2+Cz^2=1, a 
// point Q=(x,y,z) on its surface, and the unit outward normal vector nn at Q,
// return (kmax,kmin) the principal curvatures at Q, and a Transform with 
// x=dmax, y=dmin, z=nn, O=Q that gives the principal curvature directions. 
// (Note: A=1/a^2, B=1/b^2, C=1/c^2 where a,b,c are the ellipsoid radii.)
void ContactGeometry::Ellipsoid::Impl::
findParaboloidAtPointWithNormal(const Vec3& Q, const UnitVec3& nn,
                                Transform& X_EP, Vec2& k) const
{
    const Real A = square(curvatures[0]), B = square(curvatures[1]), 
               C = square(curvatures[2]);

    // Sanity checks in debug.
    SimTK_ERRCHK(std::abs(A*Q[0]*Q[0]+B*Q[1]*Q[1]+C*Q[2]*Q[2]-1) < SqrtEps,
        "ContactGeometry::Ellipsoid::findParaboloidAtPointWithNormal()",
        "The given point was not on the surface of the ellipsoid.");
    SimTK_ERRCHK((nn-findUnitNormalAtPoint(Q)).normSqr() < SqrtEps,
        "ContactGeometry::Ellipsoid::findParaboloidAtPointWithNormal()",
        "The given normal was not consistent with the given point.");

    UnitVec3 tu = nn.perp();    // ~40 flops
    UnitVec3 tv(nn % tu, true); // y = z X x for plane, already normalized (9 flops)
    
    // 27 flops to get R,S,T
    Real R=   A*square(tu[0]) + B*square(tu[1]) + C*square(tu[2]);
    Real S=   A*square(tv[0]) + B*square(tv[1]) + C*square(tv[2]);
    Real T=2*(A*tu[0]*tv[0]   + B*tu[1]*tv[1]   + C*tu[2]*tv[2]);

    // T will be zero for spheres (A=B=C) and for various "clean" points
    // on the ellipsoid where tu[i]*tv[i]==0, i=0,1,2. In that case we
    // already have the ellipse we're looking for with R,S.
    // R==S means curvature is the same in every direction (that's called
    // an "umbilic" point). In that case tu and tv are good directions.
    // I *believe* R==S -> T==0 but I don't have a proof.
    Real kmax2, kmin2; // squared curvatures of ellipse
    UnitVec3 dmax;
    if (std::abs(R-S) < SignificantReal*std::max(R,S)) {
        kmax2 = kmin2 = (R+S)/2;
        dmax = tu;
    } else if (std::abs(T) < SignificantReal) {
        if (R < S) kmax2=S, dmax=tv, kmin2=R;
        else       kmax2=R, dmax=tu, kmin2=S;
    } else { // T,R-S both nonzero
        Real tan2t = T/(R-S);       // ~20 flops
        Real cos2t = 1/std::sqrt(1 + square(tan2t)); // ~40 flops
        Real sin2t = tan2t*cos2t;   //   1 flop
        // 11 flops here
        Real term = R*cos2t-S*cos2t+T*sin2t;
        Real Rp = (R+S + term)/2;
        Real Sp = (R+S - term)/2;

        // Sort into kmax, kmin; at most one normalization done below
        if (Rp < Sp) {
            kmax2=Sp, kmin2=Rp;
            dmax = UnitVec3((1+cos2t)*tv - sin2t*tu); // Sdir, must normalize, ~50 flops
        } else {
            kmax2=Rp,kmin2=Sp;
            dmax = UnitVec3((1+cos2t)*tu + sin2t*tv); // Rdir, must normalize, ~50 flops
        }
    }

    Real d = ~Q * nn; // distance along normal from center to point on ellipsoid (5 flops)
    Real kmax = d * kmax2, kmin = d * kmin2; // surface curvatures (2 flops)

    X_EP.updP() = Q; // the origin point
    Rotation& R_EP = X_EP.updR();
    // 9 flops
    UnitVec3 dmin = UnitVec3(nn % dmax, true); // y=z%x ensures right handedness (already unit vector too)
    R_EP.setRotationFromUnitVecsTrustMe(dmax, dmin, nn);

    k = Vec2(kmax, kmin);
}


// Peter E. says he implemented this from David Eberly's web site
// http://www.geometrictools.com/Documentation/DistancePointToEllipsoid.pdf
// Eberly says he got it from John Hart's article in Graphics Gems 4, page
// 113 "Distance to an Ellipsoid". Both Eberly and Hart recommend using a
// Newton iteration to solve this problem because the largest root is directly
// downhill given appropriate starting points, which they provide. However,
// the implementation here uses a direct solution of the 6th-order polynomial
// then searches for the largest real root. That is likely to be *much* slower
// than the recommended approach, although that should be measured.
//
// I asked Peter and he said he did not try and reject the Newton approach;
// he just took the direct approach. I believe the Newton method would be
// *much* faster, but Eberly hints that there are special cases that can
// cause convergence troubles and must be dealt with carefully. If the
// existing routine turns out to be a bottleneck, it would be worth revisiting
// this implementation. -- Sherm 20110203.
//
// TODO: use faster method?
Vec3 ContactGeometry::Ellipsoid::Impl::
findNearestPoint(const Vec3& position, bool& inside, UnitVec3& normal) const {
    Real a2 = radii[0]*radii[0];
    Real b2 = radii[1]*radii[1];
    Real c2 = radii[2]*radii[2];
    Real a4 = a2*a2;
    Real b4 = b2*b2;
    Real c4 = c2*c2;
    Real px2 = position[0]*position[0];
    Real py2 = position[1]*position[1];
    Real pz2 = position[2]*position[2];
    Real a2b2 = a2*b2;
    Real b2c2 = b2*c2;
    Real a2c2 = a2*c2;
    Real a2b2c2 = a2b2*c2;
    Vector coeff(7);
    coeff[0] = 1;
    coeff[1] = 2*(a2+b2+c2);
    coeff[2] = -(a2*px2+b2*py2+c2*pz2) + a4+b4+c4 + 4*(a2b2+b2c2+a2c2);
    coeff[3] = -2*((a2b2+a2c2)*px2+(a2b2+b2c2)*py2+(b2c2+a2c2)*pz2) + 2*(a4*(b2+c2)+b4*(a2+c2)+c4*(a2+b2)) + 8*a2b2c2;
    coeff[4] = -a2*(b4+4*b2c2+c4)*px2-b2*(a4+4*a2c2+c4)*py2-c2*(a4+4*a2b2+b4)*pz2 + 4*(a2+b2+c2)*a2b2c2 + a4*b4+a4*c4+b4*c4;
    coeff[5] = 2*a2b2c2*(-(b2+c2)*px2-(a2+c2)*py2-(a2+b2)*pz2 + a2b2+b2c2+a2c2);
    coeff[6] = a2b2c2*(-b2c2*px2-a2c2*py2-a2b2*pz2+a2b2c2);
    Vector_<complex<Real> > roots(6);
    PolynomialRootFinder::findRoots(coeff, roots);
    Real root = NTraits<Real>::getMostNegative();
    for (int i = 0; i < 6; i++)
        if (fabs(roots[i].imag()) < 1e-10 && (roots[i].real()) > (root))
            root = roots[i].real();
    Vec3 result(position[0]*a2/(root+a2), position[1]*b2/(root+b2), position[2]*c2/(root+c2));
    Vec3 ri2(1/a2, 1/b2, 1/c2);
    inside = (position[0]*position[0]*ri2[0] + position[1]*position[1]*ri2[1] + position[2]*position[2]*ri2[2] < 1.0);
    normal = UnitVec3(result[0]*ri2[0], result[1]*ri2[1], result[2]*ri2[2]);
    return result;
}

// Peter says he took this algorithm from Art of Illusion but can't remember
// where it came from. It is similar to an algorithm presented in this thread:
// http://www.ogre3d.org/forums/viewtopic.php?f=2&t=26442&start=0
// and is most likely a special case of the general ray-quadric intersection
// method presented by Cychosz and Waggenspack in Graphics Gems III, pg. 275,
// "Intersecting a ray with a quadric surface."
bool ContactGeometry::Ellipsoid::Impl::intersectsRay
   (const Vec3& origin, const UnitVec3& direction,
    Real& distance, UnitVec3& normal) const
{
    Real rx2 = radii[0]*radii[0];
    Real sy = rx2/(radii[1]*radii[1]);
    Real sz = rx2/(radii[2]*radii[2]);
    Vec3 scaledDir(direction[0], sy*direction[1], sz*direction[2]);
    Real b = -(~scaledDir*origin);
    Real c = origin[0]*origin[0] + sy*origin[1]*origin[1] + sz*origin[2]*origin[2] - rx2;
    if (c > 0) {
        // Ray origin is outside ellipsoid.

        if (b <= 0)
          return false;  // Ray points away from the ellipsoid.
        Real a = ~scaledDir*direction;;
        Real d = b*b - a*c;
        if (d < 0)
          return false;
        distance = (b - std::sqrt(d))/a;
    }
    else {
        // Ray origin is inside ellipsoid.

        Real a = ~scaledDir*direction;;
        Real d = b*b - a*c;
        if (d < 0)
          return false;
        distance = (b + std::sqrt(d))/a;
    }
    Vec3 pos = origin+distance*direction;
    normal = UnitVec3(pos[0], pos[1]*sy, pos[2]*sz);
    return true;
}

void ContactGeometry::Ellipsoid::Impl::
getBoundingSphere(Vec3& center, Real& radius) const {
    center = Vec3(0);
    radius = max(radii);
}

void ContactGeometry::Ellipsoid::Impl::
calcCurvature(const Vec3& point, Vec2& curvature, Rotation& orientation) const {
    Transform transform;
    findParaboloidAtPoint(point, transform, curvature);
    orientation = transform.R();
}

Real EllipsoidImplicitFunction::
calcValue(const Vector& x) const {
    const Vec3& radii = ownerp->getRadii();
    return 1-x[0]*x[0]/(radii[0]*radii[0])-x[1]*x[1]/(radii[1]*radii[1])-x[2]*x[2]/(radii[2]*radii[2]);
}

Real EllipsoidImplicitFunction::
calcDerivative(const Array_<int>& derivComponents, const Vector& x) const {
    const Vec3& radii = ownerp->getRadii();
    if (derivComponents.size() == 1) {
        int c = derivComponents[0];
        return 2*x[c]/(radii[c]*radii[c]);
    }
    if (derivComponents[0] == derivComponents[1]) {
        int c = derivComponents[0];
        return 2/(radii[c]*radii[c]);
    }
    return 0;
}



//==============================================================================
//                          SMOOTH HEIGHT MAP & IMPL
//==============================================================================

ContactGeometry::SmoothHeightMap::
SmoothHeightMap(const BicubicSurface& surface) 
:   ContactGeometry(new SmoothHeightMap::Impl(surface)) {}

/*static*/ ContactGeometryTypeId ContactGeometry::SmoothHeightMap::
classTypeId() 
{   return ContactGeometry::SmoothHeightMap::Impl::classTypeId(); }

const BicubicSurface& ContactGeometry::SmoothHeightMap::
getBicubicSurface() const {return getImpl().getBicubicSurface();}

const ContactGeometry::SmoothHeightMap::Impl& ContactGeometry::SmoothHeightMap::
getImpl() const {
    assert(impl);
    return static_cast<const SmoothHeightMap::Impl&>(*impl);
}

ContactGeometry::SmoothHeightMap::Impl& ContactGeometry::SmoothHeightMap::
updImpl() {
    assert(impl);
    return static_cast<SmoothHeightMap::Impl&>(*impl);
}

// This is the main constructor.
ContactGeometry::SmoothHeightMap::Impl::
Impl(const BicubicSurface& surface) 
:   surface(surface) { 
    implicitFunction.setOwner(*this); 

    // Create bounding sphere.
    // TODO: fake this using mesh; this needs to be done correctly instead
    // by the BicubicSurface itself. Using 5 subdivisions per patch.
    PolygonalMesh mesh = surface.createPolygonalMesh(5);

    // Collect all the vertices.
    const int n = mesh.getNumVertices();
    Array_<const Vec3*> points(n);
    for (int i=0; i<n; ++i)
        points[i] = &mesh.getVertexPosition(i);
    boundingSphere = Geo::Point::calcBoundingSphere(points);
    // Add 10% as a hack to make it less likely we'll miss part of the surface.
    boundingSphere.updRadius() *= 1.1;
}

// This constructor is used by clone() to avoid recalculating the bounding
// sphere.
ContactGeometry::SmoothHeightMap::Impl::
Impl(const BicubicSurface& surface, const Geo::Sphere& boundingSphere) 
:   surface(surface), boundingSphere(boundingSphere) 
{   implicitFunction.setOwner(*this); }

Vec3 ContactGeometry::SmoothHeightMap::Impl::
findNearestPoint(const Vec3& position, bool& inside, UnitVec3& normal) const {
    assert(false);
    return Vec3(NaN);
}

bool ContactGeometry::SmoothHeightMap::Impl::intersectsRay
   (const Vec3& origin, const UnitVec3& direction, 
    Real& distance, UnitVec3& normal) const 
{
    assert(false);
    return true;
}

Real SmoothHeightMapImplicitFunction::
calcValue(const Vector& p) const {
    const BicubicSurface&      surf = ownerp->getBicubicSurface();
    BicubicSurface::PatchHint& hint = ownerp->updHint();
    const Real z = surf.calcValue(Vec2(p[0],p[1]), hint);
    //TODO: this is negated from convention
    return z - p[2]; // negative outside, positive inside
}

// First deriv with respect to z (component 2) is -1, all higher derivs are 0.
// Higher partials involving z are zero.
Real SmoothHeightMapImplicitFunction::
calcDerivative(const Array_<int>& derivComponents, const Vector& p) const {
    if (derivComponents.empty()) return calcValue(p);
    if (derivComponents.size() == 1 && derivComponents[0]==2)
        return -1;
    for (unsigned i=0; i<derivComponents.size(); ++i)
        if (derivComponents[i]==2) return 0;

    // We're asking only for derivatives in x and y.
    const BicubicSurface&      surf = ownerp->getBicubicSurface();
    BicubicSurface::PatchHint& hint = ownerp->updHint();
    const Real d = surf.calcDerivative(derivComponents, Vec2(p[0],p[1]), hint);
    return d;
}




//==============================================================================
//                              TRIANGLE MESH
//==============================================================================

ContactGeometry::TriangleMesh::TriangleMesh
   (const ArrayViewConst_<Vec3>& vertices, 
    const ArrayViewConst_<int>& faceIndices, bool smooth) 
:   ContactGeometry(new TriangleMesh::Impl(vertices, faceIndices, smooth)) {}

ContactGeometry::TriangleMesh::TriangleMesh
   (const PolygonalMesh& mesh, bool smooth) 
:   ContactGeometry(new TriangleMesh::Impl(mesh, smooth)) {}

/*static*/ ContactGeometryTypeId ContactGeometry::TriangleMesh::classTypeId() 
{   return ContactGeometry::TriangleMesh::Impl::classTypeId(); }


int ContactGeometry::TriangleMesh::getNumEdges() const {
    return getImpl().edges.size();
}

int ContactGeometry::TriangleMesh::getNumFaces() const {
    return getImpl().faces.size();
}

int ContactGeometry::TriangleMesh::getNumVertices() const {
    return getImpl().vertices.size();
}

const Vec3& ContactGeometry::TriangleMesh::getVertexPosition(int index) const {
    assert(index >= 0 && index < getNumVertices());
    return getImpl().vertices[index].pos;
}

int ContactGeometry::TriangleMesh::getFaceEdge(int face, int edge) const {
    assert(face >= 0 && face < getNumFaces());
    assert(edge >= 0 && edge < 3);
    return getImpl().faces[face].edges[edge];
}

int ContactGeometry::TriangleMesh::getFaceVertex(int face, int vertex) const {
    assert(face >= 0 && face < getNumFaces());
    assert(vertex >= 0 && vertex < 3);
    return getImpl().faces[face].vertices[vertex];
}

int ContactGeometry::TriangleMesh::getEdgeFace(int edge, int face) const {
    assert(edge >= 0 && edge < getNumEdges());
    assert(face >= 0 && face < 2);
    return getImpl().edges[edge].faces[face];
}

int ContactGeometry::TriangleMesh::getEdgeVertex(int edge, int vertex) const {
    assert(edge >= 0 && edge < getNumEdges());
    assert(vertex >= 0 && vertex < 2);
    return getImpl().edges[edge].vertices[vertex];
}

const UnitVec3& ContactGeometry::TriangleMesh::getFaceNormal(int face) const {
    assert(face >= 0 && face < getNumFaces());
    return getImpl().faces[face].normal;
}

Real ContactGeometry::TriangleMesh::getFaceArea(int face) const {
    assert(face >= 0 && face < getNumFaces());
    return getImpl().faces[face].area;
}

void ContactGeometry::TriangleMesh::
findVertexEdges(int vertex, Array_<int>& edges) const {
    // Begin at an arbitrary edge which intersects the vertex.
    
    int firstEdge = getImpl().vertices[vertex].firstEdge;
    int previousEdge = firstEdge;
    int previousFace = getImpl().edges[firstEdge].faces[0];
    
    // Walk around the vertex, using each edge to find the next face and each 
    // face to find the next edge.
    
    do {
        edges.push_back(previousEdge);
        const ContactGeometry::TriangleMesh::Impl::Edge& 
            edge = getImpl().edges[previousEdge];
        int nextFace = (edge.faces[0] == previousFace ? edge.faces[1] 
                                                      : edge.faces[0]);
        const ContactGeometry::TriangleMesh::Impl::Face& 
            face = getImpl().faces[nextFace];
        int nextEdge;
        if (    face.edges[0] != previousEdge
            && (face.vertices[0] == vertex || face.vertices[1] == vertex))
            nextEdge = face.edges[0];
        else if (   face.edges[1] != previousEdge 
                 && (face.vertices[1] == vertex || face.vertices[2] == vertex))
            nextEdge = face.edges[1];
        else
            nextEdge = face.edges[2];
        previousEdge = nextEdge;
        previousFace = nextFace;
    } while (previousEdge != firstEdge);
}

Vec3 ContactGeometry::TriangleMesh::findPoint(int face, const Vec2& uv) const {
    return getImpl().findPoint(face, uv);
}

Vec3 ContactGeometry::TriangleMesh::findCentroid(int face) const {
    return getImpl().findCentroid(face);
}

UnitVec3 ContactGeometry::TriangleMesh::
findNormalAtPoint(int face, const Vec2& uv) const {
    return getImpl().findNormalAtPoint(face, uv);
}

Vec3 ContactGeometry::TriangleMesh::findNearestPoint
   (const Vec3& position, bool& inside, UnitVec3& normal) const {
    return getImpl().findNearestPoint(position, inside, normal);
}

Vec3 ContactGeometry::TriangleMesh::findNearestPoint
   (const Vec3& position, bool& inside, int& face, Vec2& uv) const {
    return getImpl().findNearestPoint(position, inside, face, uv);
}

Vec3 ContactGeometry::TriangleMesh::findNearestPointToFace
   (const Vec3& position, int face, Vec2& uv) const {
    return getImpl().findNearestPointToFace(position, face, uv);
}

bool ContactGeometry::TriangleMesh::intersectsRay
   (const Vec3& origin, const UnitVec3& direction, Real& distance, 
    UnitVec3& normal) const {
    return getImpl().intersectsRay(origin, direction, distance, normal);
}

bool ContactGeometry::TriangleMesh::intersectsRay
   (const Vec3& origin, const UnitVec3& direction, Real& distance, int& face, 
    Vec2& uv) const {
    return getImpl().intersectsRay(origin, direction, distance, face, uv);
}

ContactGeometry::TriangleMesh::OBBTreeNode 
ContactGeometry::TriangleMesh::getOBBTreeNode() const {
    return OBBTreeNode(getImpl().obb);
}

PolygonalMesh ContactGeometry::TriangleMesh::createPolygonalMesh() const {
    PolygonalMesh mesh;
    getImpl().createPolygonalMesh(mesh);
    return mesh;
}

const ContactGeometry::TriangleMesh::Impl& 
ContactGeometry::TriangleMesh::getImpl() const {
    assert(impl);
    return static_cast<const TriangleMesh::Impl&>(*impl);
}

ContactGeometry::TriangleMesh::Impl& 
ContactGeometry::TriangleMesh::updImpl() {
    assert(impl);
    return static_cast<TriangleMesh::Impl&>(*impl);
}



//==============================================================================
//                            TRIANGLE MESH IMPL
//==============================================================================

Vec3 ContactGeometry::TriangleMesh::Impl::findPoint
   (int face, const Vec2& uv) const {
    const Face& f = faces[face];
    return             uv[0] * vertices[f.vertices[0]].pos
           +           uv[1] * vertices[f.vertices[1]].pos
           +  (1-uv[0]-uv[1])* vertices[f.vertices[2]].pos;
}

// same as findPoint(face, (1/3,1/3)) but faster
Vec3 ContactGeometry::TriangleMesh::Impl::findCentroid(int face) const {
    const Face& f = faces[face];
    return (  vertices[f.vertices[0]].pos
            + vertices[f.vertices[1]].pos
            + vertices[f.vertices[2]].pos) / 3;
}

UnitVec3 ContactGeometry::TriangleMesh::Impl::findNormalAtPoint
   (int face, const Vec2& uv) const {
    const Face& f = faces[face];
    if (smooth)
        return UnitVec3(            uv[0] * vertices[f.vertices[0]].normal
                        +           uv[1] * vertices[f.vertices[1]].normal
                        +  (1-uv[0]-uv[1])* vertices[f.vertices[2]].normal);
    return f.normal;
}

Vec3 ContactGeometry::TriangleMesh::Impl::
findNearestPoint(const Vec3& position, bool& inside, UnitVec3& normal) const {
    int face;
    Vec2 uv;
    Vec3 nearestPoint = findNearestPoint(position, inside, face, uv);
    normal = findNormalAtPoint(face, uv);
    return nearestPoint;
}

Vec3 ContactGeometry::TriangleMesh::Impl::
findNearestPoint(const Vec3& position, bool& inside, int& face, Vec2& uv) const 
{
    Real distance2;
    Vec3 nearestPoint = obb.findNearestPoint(*this, position, MostPositiveReal, distance2, face, uv);
    Vec3 delta = position-nearestPoint;
    inside = (~delta*faces[face].normal < 0);
    return nearestPoint;
}

bool ContactGeometry::TriangleMesh::Impl::
intersectsRay(const Vec3& origin, const UnitVec3& direction, Real& distance, 
              UnitVec3& normal) const {
    int face;
    Vec2 uv;
    if (!intersectsRay(origin, direction, distance, face, uv))
        return false;
    normal = findNormalAtPoint(face, uv);
    return true;
}

bool ContactGeometry::TriangleMesh::Impl::
intersectsRay(const Vec3& origin, const UnitVec3& direction, Real& distance, 
              int& face, Vec2& uv) const {
    Real boundsDistance;
    if (!obb.bounds.intersectsRay(origin, direction, boundsDistance))
        return false;
    return obb.intersectsRay(*this, origin, direction, distance, face, uv);
}

void ContactGeometry::TriangleMesh::Impl::
getBoundingSphere(Vec3& center, Real& radius) const {
    center = boundingSphereCenter;
    radius = boundingSphereRadius;
}

void ContactGeometry::TriangleMesh::Impl::
createPolygonalMesh(PolygonalMesh& mesh) const {
    for (unsigned vx=0; vx < vertices.size(); ++vx)
        mesh.addVertex(vertices[vx].pos);
    for (unsigned fx=0; fx < faces.size(); ++fx) {
        const Face& face = faces[fx];
        const ArrayViewConst_<int> verts(face.vertices, face.vertices+3);
        mesh.addFace(verts);
    }
}

ContactGeometry::TriangleMesh::Impl::Impl
   (const ArrayViewConst_<Vec3>& vertexPositions, 
    const ArrayViewConst_<int>& faceIndices, bool smooth) 
:   ContactGeometryImpl(), smooth(smooth) {
    init(vertexPositions, faceIndices);
}

ContactGeometry::TriangleMesh::Impl::Impl
   (const PolygonalMesh& mesh, bool smooth) 
:   ContactGeometryImpl(), smooth(smooth) 
{   // Create the mesh, triangulating faces as necessary.
    Array_<Vec3>    vertexPositions;
    Array_<int>     faceIndices;
    for (int i = 0; i < mesh.getNumVertices(); i++)
        vertexPositions.push_back(mesh.getVertexPosition(i));
    for (int i = 0; i < mesh.getNumFaces(); i++) {
        int numVert = mesh.getNumVerticesForFace(i);
        if (numVert < 3)
            continue; // Ignore it.
        if (numVert == 3) {
            faceIndices.push_back(mesh.getFaceVertex(i, 0));
            faceIndices.push_back(mesh.getFaceVertex(i, 1));
            faceIndices.push_back(mesh.getFaceVertex(i, 2));
        }
        else if (numVert == 4) {
            // Split it into two triangles.
            
            faceIndices.push_back(mesh.getFaceVertex(i, 0));
            faceIndices.push_back(mesh.getFaceVertex(i, 1));
            faceIndices.push_back(mesh.getFaceVertex(i, 2));
            faceIndices.push_back(mesh.getFaceVertex(i, 2));
            faceIndices.push_back(mesh.getFaceVertex(i, 3));
            faceIndices.push_back(mesh.getFaceVertex(i, 0));
        }
        else {
            // Add a vertex at the center, then split it into triangles.
            
            Vec3 center(0);
            for (int j = 0; j < numVert; j++)
                center += vertexPositions[mesh.getFaceVertex(i, j)];
            center /= numVert;
            vertexPositions.push_back(center);
            int newIndex = vertexPositions.size()-1;
            for (int j = 0; j < numVert-1; j++) {
                faceIndices.push_back(mesh.getFaceVertex(i, j));
                faceIndices.push_back(mesh.getFaceVertex(i, j+1));
                faceIndices.push_back(newIndex);
            }
        }
    }
    init(vertexPositions, faceIndices);
    
    // Make sure the mesh normals are oriented correctly.
    
    Vec3 origin(0);
    for (int i = 0; i < 3; i++)
        origin += vertices[faces[0].vertices[i]].pos;
    origin /= 3; // this is the face centroid

    const UnitVec3 direction = -faces[0].normal;
    // Calculate a ray origin that is guaranteed to be outside the
    // mesh. If the topology is right (face 0 normal points outward), we'll be
    // outside on the side containing face 0. If it is wrong, we'll be outside
    // on the opposite side of the mesh. Then we'll shoot a ray back along the
    // direction we came from (that is, towards the interior of the mesh from
    // outside). We'll hit *some* face. If the topology is right, the hit 
    // face's normal will be pointing back at us. If it is wrong, the face 
    // normal will also be pointing inwards, in roughly the same direction as 
    // the ray.
    origin -= max(obb.bounds.getSize())*direction;
    Real distance;
    int face;
    Vec2 uv;
    bool intersects = intersectsRay(origin, direction, distance, face, uv);
    assert(intersects);
    // Now dot the hit face normal with the ray direction; correct topology
    // will have them pointing in more-or-less opposite directions.
    if (dot(faces[face].normal, direction) > 0) {
        // We need to invert the mesh topology.
        
        for (int i = 0; i < (int) faces.size(); i++) {
            Face& f = faces[i];
            int temp = f.vertices[0];
            f.vertices[0] = f.vertices[1];
            f.vertices[1] = temp;
            temp = f.edges[1];
            f.edges[1] = f.edges[2];
            f.edges[2] = temp;
            f.normal *= -1;
        }
        for (int i = 0; i < (int) vertices.size(); i++)
            vertices[i].normal *= -1;
    }
}

void ContactGeometry::TriangleMesh::Impl::init
   (const Array_<Vec3>& vertexPositions, const Array_<int>& faceIndices) 
{   SimTK_APIARGCHECK_ALWAYS(faceIndices.size()%3 == 0, 
        "ContactGeometry::TriangleMesh::Impl", "TriangleMesh::Impl", 
        "The number of indices must be a multiple of 3.");
    int numFaces = faceIndices.size()/3;
    
    // Create the vertices.
    
    for (int i = 0; i < (int) vertexPositions.size(); i++)
        vertices.push_back(Vertex(vertexPositions[i]));
    
    // Create the faces and build lists of all the edges.
    
    map<pair<int, int>, int> forwardEdges;
    map<pair<int, int>, int> backwardEdges;
    for (int i = 0; i < numFaces; i++) {
        int start = i*3;
        int v1 = faceIndices[start], v2 = faceIndices[start+1], 
            v3 = faceIndices[start+2];
        SimTK_APIARGCHECK1_ALWAYS
           (   v1 >= 0 && v1 < (int) vertices.size() 
            && v2 >= 0 && v2 < (int) vertices.size() 
            && v3 >= 0 && v3 < (int) vertices.size(),
            "ContactGeometry::TriangleMesh::Impl", "TriangleMesh::Impl",
            "Face %d contains a vertex with an illegal index.", i);
        Vec3 cross =   (vertexPositions[v2]-vertexPositions[v1])
                     % (vertexPositions[v3]-vertexPositions[v1]);
        Real norm = cross.norm();
        cross *= 1.0/norm;
        SimTK_APIARGCHECK1_ALWAYS(norm > 0, 
            "ContactGeometry::TriangleMesh::Impl", "TriangleMesh::Impl",
            "Face %d is degenerate.", i);
        faces.push_back(Face(v1, v2, v3, cross, 0.5*norm));
        int edges[3][2] = {{v1, v2}, {v2, v3}, {v3, v1}};
        for (int j = 0; j < 3; j++) {
            SimTK_APIARGCHECK1_ALWAYS(edges[j][0] != edges[j][1], 
                "ContactGeometry::TriangleMesh::Impl", "TriangleMesh::Impl",
                "Vertices %d appears twice in a single face.", edges[j][0]);
            if (edges[j][0] < edges[j][1]) {
                SimTK_APIARGCHECK2_ALWAYS
                   (forwardEdges.find(pair<int, int>(edges[j][0], edges[j][1])) 
                    == forwardEdges.end(),
                    "ContactGeometry::TriangleMesh::Impl", "TriangleMesh::Impl",
                    "Multiple faces have an edge between vertices %d and %d"
                    " in the same order.", edges[j][0], edges[j][1]);
                forwardEdges[pair<int, int>(edges[j][0], edges[j][1])] = i;
            }
            else {
                SimTK_APIARGCHECK2_ALWAYS
                   (backwardEdges.find(pair<int, int>(edges[j][1], edges[j][0]))
                    == backwardEdges.end(),
                    "ContactGeometry::TriangleMesh::Impl", "TriangleMesh::Impl",
                    "Multiple faces have an edge between vertices %d and %d"
                    " in the same order.", edges[j][1], edges[j][0]);
                backwardEdges[pair<int, int>(edges[j][1], edges[j][0])] = i;
            }
        }
    }
    
    // Create the edges.
    
    SimTK_APIARGCHECK_ALWAYS(forwardEdges.size() == backwardEdges.size(), 
        "ContactGeometry::TriangleMesh::Impl", "TriangleMesh::Impl",
        "Each edge must be shared by exactly two faces.");
    for (map<pair<int, int>, int>::iterator iter = forwardEdges.begin(); 
         iter != forwardEdges.end(); ++iter) {
        int vert1 = iter->first.first;
        int vert2 = iter->first.second;
        int face1 = iter->second;
        map<pair<int, int>, int>::iterator iter2 = 
            backwardEdges.find(pair<int, int>(vert1, vert2));
        SimTK_APIARGCHECK_ALWAYS(iter2 != backwardEdges.end(), 
            "ContactGeometry::TriangleMesh::Impl", "TriangleMesh::Impl",
            "Each edge must be shared by exactly two faces.");
        int face2 = iter2->second;
        edges.push_back(Edge(vert1, vert2, face1, face2));
    }
    
    // Record the edges for each face.
    
    for (int i = 0; i < (int) edges.size(); i++) {
        Edge& edge = edges[i];
        int f[2] = {edge.faces[0], edge.faces[1]};
        for (int j = 0; j < 2; j++) {
            Face& face = faces[f[j]];
            if ((edge.vertices[0] == face.vertices[0] || edge.vertices[0] == face.vertices[1]) &&
                    (edge.vertices[1] == face.vertices[0] || edge.vertices[1] == face.vertices[1]))
                face.edges[0] = i;
            else if ((edge.vertices[0] == face.vertices[1] || edge.vertices[0] == face.vertices[2]) &&
                    (edge.vertices[1] == face.vertices[1] || edge.vertices[1] == face.vertices[2]))
                face.edges[1] = i;
            else if ((edge.vertices[0] == face.vertices[2] || edge.vertices[0] == face.vertices[0]) &&
                    (edge.vertices[1] == face.vertices[2] || edge.vertices[1] == face.vertices[0]))
                face.edges[2] = i;
            else
                SimTK_ASSERT_ALWAYS(false, 
                    "Face and edge vertices are inconsistent.");
        }
    }
    
    // Record a single edge for each vertex.
    
    for (int i = 0; i < (int) edges.size(); i++) {
        vertices[edges[i].vertices[0]].firstEdge = i;
        vertices[edges[i].vertices[1]].firstEdge = i;
    }
    for (int i = 0; i < (int) vertices.size(); i++)
        SimTK_APIARGCHECK1_ALWAYS(vertices[i].firstEdge >= 0, 
            "ContactGeometry::TriangleMesh::Impl", "TriangleMesh::Impl",
            "Vertex %d is not part of any face.", i);
    
    // Calculate a normal for each vertex.
    
    Vector_<Vec3> vertNorm(vertices.size(), Vec3(0));
    for (int i = 0; i < (int) faces.size(); i++) {
        const Face& f = faces[i];
        UnitVec3 edgeDir[3];
        for (int j = 0; j < 3; j++) {
            edgeDir[j] = UnitVec3(  vertices[f.vertices[(j+1)%3]].pos
                                  - vertices[f.vertices[j]].pos);
        }
        for (int j = 0; j < 3; j++) {
            Real angle = std::acos(~edgeDir[j]*edgeDir[(j+2)%3]);
            vertNorm[f.vertices[j]] += f.normal*angle;
        }
    }
    for (int i = 0; i < (int) vertices.size(); i++)
        vertices[i].normal = UnitVec3(vertNorm[i]);
    
    // Create the OBBTree.
    
    Array_<int> allFaces(faces.size());
    for (int i = 0; i < (int) allFaces.size(); i++)
        allFaces[i] = i;
    createObbTree(obb, allFaces);
    
    // Find the bounding sphere.
    Array_<const Vec3*> points(vertices.size());
    for (int i = 0; i < (int) vertices.size(); i++)
        points[i] = &vertices[i].pos;
    const Geo::Sphere bnd = Geo::Point::calcBoundingSphere(points);
    boundingSphereCenter = bnd.getCenter();
    boundingSphereRadius = bnd.getRadius();
}

void ContactGeometry::TriangleMesh::Impl::createObbTree
   (OBBTreeNodeImpl& node, const Array_<int>& faceIndices) 
{   // Find all vertices in the node and build the OrientedBoundingBox.
    node.numTriangles = faceIndices.size();
    set<int> vertexIndices;
    for (int i = 0; i < (int) faceIndices.size(); i++) 
        for (int j = 0; j < 3; j++)
            vertexIndices.insert(faces[faceIndices[i]].vertices[j]);
    Vector_<Vec3> points((int)vertexIndices.size());
    int index = 0;
    for (set<int>::iterator iter = vertexIndices.begin(); 
                            iter != vertexIndices.end(); ++iter)
        points[index++] = vertices[*iter].pos;
    node.bounds = OrientedBoundingBox(points);
    if (faceIndices.size() > 3) {

        // Order the axes by size.

        int axisOrder[3];
        const Vec3& size = node.bounds.getSize();
        if (size[0] > size[1]) {
            if (size[0] > size[2]) {
                axisOrder[0] = 0;
                if (size[1] > size[2]) {
                    axisOrder[1] = 1;
                    axisOrder[2] = 2;
                }
                else {
                    axisOrder[1] = 2;
                    axisOrder[2] = 1;
                }
            }
            else {
                axisOrder[0] = 2;
                axisOrder[1] = 0;
                axisOrder[2] = 1;
            }
        }
        else if (size[0] > size[2]) {
            axisOrder[0] = 1;
            axisOrder[1] = 0;
            axisOrder[2] = 2;
        }
        else {
            if (size[1] > size[2]) {
                axisOrder[0] = 1;
                axisOrder[1] = 2;
            }
            else {
                axisOrder[0] = 2;
                axisOrder[1] = 1;
            }
            axisOrder[2] = 0;
        }

        // Try splitting along each axis.

        for (int i = 0; i < 3; i++) {
            Array_<int> child1Indices, child2Indices;
            splitObbAxis(faceIndices, child1Indices, child2Indices, 
                         axisOrder[i]);
            if (child1Indices.size() > 0 && child2Indices.size() > 0) {
                // It was successfully split, so create the child nodes.

                node.child1 = new OBBTreeNodeImpl();
                node.child2 = new OBBTreeNodeImpl();
                createObbTree(*node.child1, child1Indices);
                createObbTree(*node.child2, child2Indices);
                return;
            }
        }
    }
    
    // This is a leaf node.
    
    node.triangles.insert(node.triangles.begin(), faceIndices.begin(), 
                          faceIndices.end());
}

void ContactGeometry::TriangleMesh::Impl::splitObbAxis
   (const Array_<int>& parentIndices, Array_<int>& child1Indices, 
    Array_<int>& child2Indices, int axis) 
{   // For each face, find its minimum and maximum extent along the axis.
    Vector minExtent(parentIndices.size());
    Vector maxExtent(parentIndices.size());
    for (int i = 0; i < (int) parentIndices.size(); i++) {
        int* vertexIndices = faces[parentIndices[i]].vertices;
        Real minVal = vertices[vertexIndices[0]].pos[axis];
        Real maxVal = vertices[vertexIndices[0]].pos[axis];
        minVal = std::min(minVal, vertices[vertexIndices[1]].pos[axis]);
        maxVal = std::max(maxVal, vertices[vertexIndices[1]].pos[axis]);
        minExtent[i] = std::min(minVal, vertices[vertexIndices[2]].pos[axis]);
        maxExtent[i] = std::max(maxVal, vertices[vertexIndices[2]].pos[axis]);
    }
    
    // Select a split point that tries to put as many faces as possible 
    // entirely on one side or the other.
    
    Real split = 0.5*(median(minExtent)+median(maxExtent));
    
    // Choose a side for each face.
    
    for (int i = 0; i < (int) parentIndices.size(); i++) {
        if (maxExtent[i] <= split)
            child1Indices.push_back(parentIndices[i]);
        else if (minExtent[i] >= split)
            child2Indices.push_back(parentIndices[i]);
        else if (0.5*(minExtent[i]+maxExtent[i]) <= split)
            child1Indices.push_back(parentIndices[i]);
        else
            child2Indices.push_back(parentIndices[i]);
    }
}

Vec3 ContactGeometry::TriangleMesh::Impl::findNearestPointToFace
   (const Vec3& position, int face, Vec2& uv) const {
    // Calculate the distance between a point in space and a face of the mesh.
    // This algorithm is based on a description by David Eberly found at 
    // http://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf.
    
    const ContactGeometry::TriangleMesh::Impl::Face& fc = faces[face];
    const Vec3& vert1 = vertices[fc.vertices[0]].pos;
    const Vec3& vert2 = vertices[fc.vertices[1]].pos;
    const Vec3& vert3 = vertices[fc.vertices[2]].pos;
    const Vec3 e0 = vert2-vert1;
    const Vec3 e1 = vert3-vert1;
    const Vec3 delta = vert1-position;
    const Real a = e0.normSqr();
    const Real b = ~e0*e1;
    const Real c = e1.normSqr();
    const Real d = ~e0*delta;
    const Real e = ~e1*delta;
    const Real f = delta.normSqr();
    const Real det = a*c-b*b;
    Real s = b*e-c*d;
    Real t = b*d-a*e;
    if (s+t <= det) {
        if (s < 0) {
            if (t < 0) {
                // Region 4

                if (d < 0) {
                    s = (-d >= a ? 1 : -d/a);
                    t = 0;
                }
                else {
                    s = 0;
                    t = (e >= 0 ? 0 : (-e >= c ? 1 : -e/c));
                }
            }
            else {
                // Region 3

                s = 0;
                t = (e >= 0 ? 0 : (-e >= c ? 1 : -e/c));
            }
        }
        else if (t < 0) {
            // Region 5

            s = (d >= 0 ? 0 : (-d >= a ? 1 : -d/a));
            t = 0;
        }
        else {
            // Region 0

            const Real invDet = 1.0/det;
            s *= invDet;
            t *= invDet;
        }
    }
    else {
        if (s < 0) {
            // Region 2

            Real temp0 = b+d;
            Real temp1 = c+e;
            if (temp1 > temp0) {
                Real numer = temp1-temp0;
                Real denom = a-2*b+c;
                s = (numer >= denom ? 1 : numer/denom);
                t = 1-s;
            }
            else {
                s = 0;
                t = (temp1 <= 0 ? 1 : (e >= 0 ? 0 : -e/c));
            }
        }
        else if (t < 0) {
            // Region 6

            Real temp0 = b+e;
            Real temp1 = a+d;
            if (temp1 > temp0) {
                Real numer = temp1-temp0;
                Real denom = a-2*b+c;
                t = (numer >= denom ? 1 : numer/denom);
                s = 1-t;
            }
            else {
                s = (temp1 <= 0 ? 1 : (e >= 0 ? 0 : -d/a));
                t = 0;
            }
        }
        else {
            // Region 1

            const Real numer = c+e-b-d;
            if (numer <= 0)
                s = 0;
            else {
                const Real denom = a-2*b+c;
                s = (numer >= denom ? 1 : numer/denom);
            }
            t = 1-s;
        }
    }
    uv = Vec2(1-s-t, s);
    return vert1 + s*e0 + t*e1;
}


//==============================================================================
//                            OBB TREE NODE IMPL
//==============================================================================

OBBTreeNodeImpl::OBBTreeNodeImpl(const OBBTreeNodeImpl& copy) 
:   bounds(copy.bounds), triangles(copy.triangles), 
    numTriangles(copy.numTriangles) {
    if (copy.child1 == NULL) {
        child1 = NULL;
        child2 = NULL;
    }
    else {
        child1 = new OBBTreeNodeImpl(*copy.child1);
        child2 = new OBBTreeNodeImpl(*copy.child2);
    }
}

OBBTreeNodeImpl::~OBBTreeNodeImpl() {
    if (child1 != NULL)
        delete child1;
    if (child2 != NULL)
        delete child2;
}

Vec3 OBBTreeNodeImpl::findNearestPoint
   (const ContactGeometry::TriangleMesh::Impl& mesh, 
    const Vec3& position, Real cutoff2, 
    Real& distance2, int& face, Vec2& uv) const 
{
    Real tol = 100*Eps;
    if (child1 != NULL) {
        // Recursively check the child nodes.
        
        Real child1distance2 = MostPositiveReal, 
             child2distance2 = MostPositiveReal;
        int child1face, child2face;
        Vec2 child1uv, child2uv;
        Vec3 child1point, child2point;
        Real child1BoundsDist2 = 
            (child1->bounds.findNearestPoint(position)-position).normSqr();
        Real child2BoundsDist2 = 
            (child2->bounds.findNearestPoint(position)-position).normSqr();
        if (child1BoundsDist2 < child2BoundsDist2) {
            if (child1BoundsDist2 < cutoff2) {
                child1point = child1->findNearestPoint(mesh, position, cutoff2, child1distance2, child1face, child1uv);
                if (child2BoundsDist2 < child1distance2 && child2BoundsDist2 < cutoff2)
                    child2point = child2->findNearestPoint(mesh, position, cutoff2, child2distance2, child2face, child2uv);
            }
        }
        else {
            if (child2BoundsDist2 < cutoff2) {
                child2point = child2->findNearestPoint(mesh, position, cutoff2, child2distance2, child2face, child2uv);
                if (child1BoundsDist2 < child2distance2 && child1BoundsDist2 < cutoff2)
                    child1point = child1->findNearestPoint(mesh, position, cutoff2, child1distance2, child1face, child1uv);
            }
        }
        if (   child1distance2 <= child2distance2*(1+tol) 
            && child2distance2 <= child1distance2*(1+tol)) {
            // Decide based on angle which one to use.
            
            if (  std::abs(~(child1point-position)*mesh.faces[child1face].normal) 
                > std::abs(~(child2point-position)*mesh.faces[child2face].normal))
                child2distance2 = MostPositiveReal;
            else
                child1distance2 = MostPositiveReal;
        }
        if (child1distance2 < child2distance2) {
            distance2 = child1distance2;
            face = child1face;
            uv = child1uv;
            return child1point;
        }
        else {
            distance2 = child2distance2;
            face = child2face;
            uv = child2uv;
            return child2point;
        }
    }    
    // This is a leaf node, so check each triangle for its distance to the point.
    
    distance2 = MostPositiveReal;
    Vec3 nearestPoint;
    for (int i = 0; i < (int) triangles.size(); i++) {
        Vec2 triangleUV;
        Vec3 p = mesh.findNearestPointToFace(position, triangles[i], triangleUV);
        Vec3 offset = p-position;
        // TODO: volatile to work around compiler bug
        volatile Real d2 = offset.normSqr(); 
        if (d2 < distance2 || (d2 < distance2*(1+tol) && std::abs(~offset*mesh.faces[triangles[i]].normal) > std::abs(~offset*mesh.faces[face].normal))) {
            nearestPoint = p;
            distance2 = d2;
            face = triangles[i];
            uv = triangleUV;
        }
    }
    return nearestPoint;
}

bool OBBTreeNodeImpl::
intersectsRay(const ContactGeometry::TriangleMesh::Impl& mesh,
              const Vec3& origin, const UnitVec3& direction, Real& distance, 
              int& face, Vec2& uv) const {
    if (child1 != NULL) {
        // Recursively check the child nodes.
        
        Real child1distance, child2distance;
        int child1face, child2face;
        Vec2 child1uv, child2uv;
        bool child1intersects = child1->bounds.intersectsRay(origin, direction, child1distance);
        bool child2intersects = child2->bounds.intersectsRay(origin, direction, child2distance);
        if (child1intersects) {
            if (child2intersects) {
                // The ray intersects both child nodes.  First check the closer one.
                
                if (child1distance < child2distance) {
                    child1intersects = child1->intersectsRay(mesh, origin,  direction, child1distance, child1face, child1uv);
                    if (!child1intersects || child2distance < child1distance)
                        child2intersects = child2->intersectsRay(mesh, origin,  direction, child2distance, child2face, child2uv);
                }
                else {
                    child2intersects = child2->intersectsRay(mesh, origin,  direction, child2distance, child2face, child2uv);
                    if (!child2intersects || child1distance < child2distance)
                        child1intersects = child1->intersectsRay(mesh, origin,  direction, child1distance, child1face, child1uv);
                }
            }
            else
                child1intersects = child1->intersectsRay(mesh, origin,  direction, child1distance, child1face, child1uv);
        }
        else if (child2intersects)
            child2intersects = child2->intersectsRay(mesh, origin,  direction, child2distance, child2face, child2uv);
        
        // If either one had an intersection, return the closer one.
        
        if (child1intersects && (!child2intersects || child1distance < child2distance)) {
            distance = child1distance;
            face = child1face;
            uv = child1uv;
            return true;
        }
        if (child2intersects) {
            distance = child2distance;
            face = child2face;
            uv = child2uv;
            return true;
        }
        return false;
    }
    
    // This is a leaf node, so check each triangle for an intersection with the 
    // ray.
    
    bool foundIntersection = false;
    for (int i = 0; i < (int) triangles.size(); i++) {
        const UnitVec3& faceNormal = mesh.faces[triangles[i]].normal;
        double vd = ~faceNormal*direction;
        if (vd == 0.0)
            continue; // The ray is parallel to the plane.
        const Vec3& vert1 = mesh.vertices[mesh.faces[triangles[i]].vertices[0]].pos;
        double v0 = ~faceNormal*(vert1-origin);
        double t = v0/vd;
        if (t < 0.0)
            continue; // Ray points away from plane of triangle.
        if (foundIntersection && t >= distance)
            continue; // We already have a closer intersection.

        // Determine whether the intersection point is inside the triangle by projecting onto
        // a plane and computing the barycentric coordinates.

        Vec3 ri = origin+direction*t;
        const Vec3& vert2 = mesh.vertices[mesh.faces[triangles[i]].vertices[1]].pos;
        const Vec3& vert3 = mesh.vertices[mesh.faces[triangles[i]].vertices[2]].pos;
        int axis1, axis2;
        if (std::abs(faceNormal[1]) > std::abs(faceNormal[0])) {
            if (std::abs(faceNormal[2]) > std::abs(faceNormal[1])) {
                axis1 = 0;
                axis2 = 1;
            }
            else {
                axis1 = 0;
                axis2 = 2;
            }
        }
        else {
            if (std::abs(faceNormal[2]) > std::abs(faceNormal[0])) {
                axis1 = 0;
                axis2 = 1;
            }
            else {
                axis1 = 1;
                axis2 = 2;
            }
        }
        Vec2 pos(ri[axis1]-vert1[axis1], ri[axis2]-vert1[axis2]);
        Vec2 edge1(vert1[axis1]-vert2[axis1], vert1[axis2]-vert2[axis2]);
        Vec2 edge2(vert1[axis1]-vert3[axis1], vert1[axis2]-vert3[axis2]);
        double denom = 1.0/(edge1%edge2);
        edge2 *= denom;
        double v = edge2%pos;
        if (v < 0.0 || v > 1.0)
            continue;
        edge1 *= denom;
        double w = pos%edge1;
        if (w < 0.0 || w > 1.0)
            continue;
        double u = 1.0-v-w;
        if (u < 0.0 || u > 1.0)
            continue;
        
        // It intersects.
        
        distance = t;
        face = triangles[i];
        uv = Vec2(u, v);
        foundIntersection = true;
    }
    return foundIntersection;
}




//==============================================================================
//                               OBB TREE NODE
//==============================================================================

ContactGeometry::TriangleMesh::OBBTreeNode::
OBBTreeNode(const OBBTreeNodeImpl& impl) : impl(&impl) {}

const OrientedBoundingBox& 
ContactGeometry::TriangleMesh::OBBTreeNode::getBounds() const {
    return impl->bounds;
}

bool ContactGeometry::TriangleMesh::OBBTreeNode::isLeafNode() const {
    return (impl->child1 == NULL);
}

const ContactGeometry::TriangleMesh::OBBTreeNode 
ContactGeometry::TriangleMesh::OBBTreeNode::getFirstChildNode() const {
    SimTK_ASSERT_ALWAYS(impl->child1, 
        "Called getFirstChildNode() on a leaf node");
    return OBBTreeNode(*impl->child1);
}

const ContactGeometry::TriangleMesh::OBBTreeNode 
ContactGeometry::TriangleMesh::OBBTreeNode::getSecondChildNode() const {
    SimTK_ASSERT_ALWAYS(impl->child2, 
        "Called getFirstChildNode() on a leaf node");
    return OBBTreeNode(*impl->child2);
}

const Array_<int>& ContactGeometry::TriangleMesh::OBBTreeNode::
getTriangles() const {
    SimTK_ASSERT_ALWAYS(impl->child2 == NULL, 
        "Called getTriangles() on a non-leaf node");
    return impl->triangles;
}

int ContactGeometry::TriangleMesh::OBBTreeNode::getNumTriangles() const {
    return impl->numTriangles;
}
