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

#include "SimTKcommon.h"
#include "simmath/internal/common.h"
#include "simmath/internal/Geo.h"
#include "simmath/internal/Geo_Point.h"
#include "simmath/internal/Geo_Sphere.h"
#include "simmath/internal/ContactGeometry.h"

#include "ContactGeometryImpl.h"

#include <iostream>
#include <cmath>
#include <map>
#include <set>

using namespace SimTK;
using std::map;
using std::pair;
using std::set;
using std::string;
using std::cout; using std::endl;


//==============================================================================
//                  CONTACT GEOMETRY :: ELLIPSOID & IMPL
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

DecorativeGeometry ContactGeometry::Ellipsoid::Impl::createDecorativeGeometry() const {
    return DecorativeEllipsoid(radii);
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


//TODO: just an axis-aligned leaf box for now
void ContactGeometry::Ellipsoid::Impl::createOBBTree() {
    OBBNode& root = obbTree.updRoot();
    root.box.setHalfLengths(radii);
    root.normal = UnitVec3(XAxis);  // doesn't matter
    root.coneHalfAngle = Pi;        // has all possible normals
    root.pointOnSurface = Vec3(radii[0],0,0); // doesn't matter
    root.children.clear(); // This is a leaf

    // Leaf contents.
    root.centerUW = Vec2(0,0);
    root.dims = Vec2(Pi, Pi/2); // u in [-Pi,Pi], v in [-Pi/2,Pi/2]
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
        return -2*x[c]/(radii[c]*radii[c]);
    }
    if (   derivComponents.size() == 2
        && derivComponents[0] == derivComponents[1]) {
        int c = derivComponents[0];
        return -2/(radii[c]*radii[c]);
    }
    // A mixed second derivative, or any higher derivative is zero.
    return 0;
}
