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
 * Contributors: Andreas Scholz, Ian Stavness                                 *
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
//                    CONTACT GEOMETRY :: SPHERE & IMPL
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

DecorativeGeometry ContactGeometry::Sphere::Impl::createDecorativeGeometry() const {
    return DecorativeSphere(radius);
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

//TODO: just an axis-aligned leaf box for now
void ContactGeometry::Sphere::Impl::createOBBTree() {
    OBBNode& root = obbTree.updRoot();
    root.box.setHalfLengths(Vec3(radius));
    root.normal = UnitVec3(XAxis);  // doesn't matter
    root.coneHalfAngle = Pi;        // has all possible normals
    root.pointOnSurface = Vec3(radius,0,0); // doesn't matter
    root.children.clear(); // This is a leaf

    // Leaf contents.
    root.centerUW = Vec2(0,0);
    root.dims = Vec2(Pi, Pi/2); // u in [-Pi,Pi], v in [-Pi/2,Pi/2]
}


// Compute geodesic between two points P and Q on a sphere analytically. Since a geodesic on a
// sphere is a great circle it is parameterized by
//
//        p(phi) = R * (e1*cos(phi) + e2*sin(phi)) ,
//
// where R is the radius of the sphere and the angle phi parameterizes the great circle with
// respect to an orthonormal basis {e1, e2}. By definition P = p(0) and the geodesic goes from
// P to Q, where Q = p(angle). Make sure e1 . e2 = 0 and |e1| = |e2| = 1.
static void setGeodesicToArc(const UnitVec3& e1, const UnitVec3& e2,
                             Real R, Real angle, Geodesic& geod)
{
    // Check if e1 and e2 are orthogonal.
    assert(std::abs(~e1*e2) <= SignificantReal);

    // Clear current geodesic.
    geod.clear();

    // TODO: Make this generic, so long geodesics are sampled more than short ones.
    const int numGeodesicSamples = 12;

    // Total arc length and orientation.
    const Real orientation = Real(sign(angle));
    const Real L = R*angle*orientation;

    // Increment of phi in loop.
    const Real deltaPhi = std::abs(angle / Real(numGeodesicSamples-1));

    const Real k = 1/R; // curvature
    for (int i = 0; i < numGeodesicSamples; ++i){
        Real phi = Real(i)*angle / Real(numGeodesicSamples-1);
        const Real sphi = sin(phi), cphi = cos(phi);

        // Trust me, this is already normalized by definition of the input.
        UnitVec3 n(e1*cphi + e2*sphi, true);

        Vec3 p = R*n;

        // t = dp/dphi, hence pointing into direction of increasing phi. 
        Vec3 t = (-e1*sphi + e2*cphi)*orientation;

        // Though not needed, we use an orthogonalizing constructor for the rotation.
        geod.addFrenetFrame(Transform(Rotation(n, ZAxis, t, YAxis), p));

        // Current arc length s.
        Real s = R*Real(i)*deltaPhi;
        geod.addArcLength(s);

        // Solve the scalar Jacobi equation
        //
        //        j''(s) + K(s)*j(s) = 0 ,                                     (1)
        //
        // where K is the Gaussian curvature and (.)' := d(.)/ds denotes differentiation
        // with respect to the arc length s. Then, j is the directional sensitivity and
        // we obtain the corresponding variational vector field by multiplying b*j. For
        // a sphere, K = R^(-2) and the solution of equation (1) becomes
        //
        //        j  = R * sin(1/R * s)                                        (2)
        //          j' =     cos(1/R * s) ,                                      (3)
        //
        // where equation (2) is the standard solution of a non-damped oscillator. Its
        // period is 2*pi*R and its amplitude is R.

        // Forward directional sensitivity from P to Q
        Vec2 jPQ(R*sin(k * s), cos(k * s));
        geod.addDirectionalSensitivityPtoQ(jPQ);

        // Backwards directional sensitivity from Q to P
        Vec2 jQP(R*sin(k * (L-s)), cos(k * (L-s)));
        geod.addDirectionalSensitivityQtoP(jQP);


        // TODO: positional sensitivity
        geod.addPositionalSensitivityPtoQ(Vec2(NaN));
        geod.addPositionalSensitivityQtoP(Vec2(NaN));

        geod.addCurvature(k);
    }
    geod.setTorsionAtP(0); geod.setTorsionAtQ(0);
    geod.setBinormalCurvatureAtP(k); geod.setBinormalCurvatureAtQ(k);

    geod.setIsConvex(true); // Curve on sphere is always convex.
    geod.setIsShortest(false); // TODO
    geod.setAchievedAccuracy(SignificantReal); // TODO: accuracy of length?
//    geod.setInitialStepSizeHint(integ.getActualInitialStepSizeTaken()); // TODO
}


void ContactGeometry::Sphere::Impl::
calcGeodesicAnalytical(const Vec3& xP, const Vec3& xQ,
                       const Vec3& tPhint, const Vec3& tQhint,
                       Geodesic& geod) const
{
    // Build an orthonormal basis {e1, e2, e3}.
    const UnitVec3 e1(xP), e_OQ(xQ);
    const Vec3 arcAxis = e1 % e_OQ;

    const Real sinAngle = arcAxis.norm();
    const Real cosAngle = ~e1*e_OQ;

    UnitVec3 e3(arcAxis/sinAngle, true);

    // Tangent vectors tP and tQ at P and Q corresponding to a positive rotation
    // of the arc around e3.
    UnitVec3 tP(e3 % e1,   true);
    UnitVec3 tQ(e3 % e_OQ, true);

    // Average moment of of hint vectors applied e3.
    Real MP = ~(e1   % tPhint)*e3;
    Real MQ = ~(e_OQ % tQhint)*e3;
    Real M  =  (MP + MQ) / 2;

    // Small angle between e_OP and e_OQ corresponding to a short geodesic.
    Real temp = atan2(sinAngle, cosAngle);
    Real angleRightHanded, angleLeftHanded;

    // Left-handed angle will always be negative, right-handed angle will be positive.
    if (temp >= 0) {
        angleRightHanded = temp;
        angleLeftHanded  = temp - 2*Pi;
    }

    else {
        angleLeftHanded  = temp;
        angleRightHanded = temp + 2*Pi;
    }

    // Orientation of arc. A negative angle means a left-handed rotation around e3. 
    Real angle;
    if (M >= 0)    {
        angle = angleRightHanded;
    }

    else {
        angle = angleLeftHanded;
    }

    // Create the last unit vector to form the orthonormal basis to describe the arc.
    UnitVec3 e2(e3 % e1, true);

    setGeodesicToArc(e1, e2, radius, angle, geod);
}

void ContactGeometry::Sphere::Impl::shootGeodesicInDirectionUntilLengthReachedAnalytical(const Vec3& xP, const UnitVec3& tP,
        const Real& terminatingLength, const GeodesicOptions& options, Geodesic& geod) const {

    UnitVec3 e_OP(xP);
    Real angle = terminatingLength/radius;

    setGeodesicToArc(e_OP, tP, radius, angle, geod);
}

void ContactGeometry::Sphere::Impl::shootGeodesicInDirectionUntilPlaneHitAnalytical(const Vec3& xP, const UnitVec3& tP,
        const Plane& terminatingPlane, const GeodesicOptions& options,
        Geodesic& geod) const {

    UnitVec3 e_OP(xP);

    // solve ~( e_OP * cos(t) + tP * sin(t) - pt_on_plane )*plane_normal = 0
    // for sphere plane offset is zero, therefore pt_on_plane = 0
    Real a = ~e_OP*terminatingPlane.getNormal();
    Real b = ~tP*terminatingPlane.getNormal();
    Real alpha = std::atan2(a,b);
    Real angle = (alpha > 0 ? Pi-alpha : -alpha);
//    std::cout << "a=" << a << ", b=" << b << ", alpha = " << alpha << std::endl;

    setGeodesicToArc(e_OP, tP, radius, angle, geod);
}


Real SphereImplicitFunction::
calcValue(const Vector& x) const {
    return 1-(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])/square(ownerp->getRadius());
}

Real SphereImplicitFunction::
calcDerivative(const Array_<int>& derivComponents, const Vector& x) const {
    if (derivComponents.size() == 1)
        return -2*x[derivComponents[0]]/square(ownerp->getRadius());
    if (   derivComponents.size() == 2 
        && derivComponents[0] == derivComponents[1])
        return -2/square(ownerp->getRadius());
    // A mixed second derivative, or any higher derivative is zero.
    return 0;
}

