/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2012 Stanford University and the Authors.           *
 * Authors: Ian Stavness                                                      *
 * Contributors: Michael Sherman, Andreas Scholz                              *
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
//                   CONTACT GEOMETRY :: CYLINDER & IMPL
//==============================================================================

ContactGeometry::Cylinder::Cylinder(Real radius)
:   ContactGeometry(new Cylinder::Impl(radius)) {}

/*static*/ ContactGeometryTypeId ContactGeometry::Cylinder::classTypeId()
{   return ContactGeometry::Cylinder::Impl::classTypeId(); }

Real ContactGeometry::Cylinder::getRadius() const {
    return getImpl().getRadius();
}

void ContactGeometry::Cylinder::setRadius(Real radius) {
    updImpl().setRadius(radius);
}

const ContactGeometry::Cylinder::Impl& ContactGeometry::Cylinder::getImpl() const {
    assert(impl);
    return static_cast<const Cylinder::Impl&>(*impl);
}

ContactGeometry::Cylinder::Impl& ContactGeometry::Cylinder::updImpl() {
    assert(impl);
    return static_cast<Cylinder::Impl&>(*impl);
}

DecorativeGeometry ContactGeometry::Cylinder::Impl::createDecorativeGeometry() const {
    DecorativeCylinder cyl(radius, radius*2);
    // DecorativeCylinder's axis is defined as the y-axis,
    // whereas ContactGeometry::Cylinder axis is defined as the z-axis
    cyl.setTransform(Rotation(Pi/2, XAxis));
    return cyl;
}

Vec3 ContactGeometry::Cylinder::Impl::findNearestPoint(const Vec3& position, bool& inside, UnitVec3& normal) const {

    normal = calcSurfaceUnitNormal(position);

    // long axis is z-axis, project to x-y plane
    Vec2 xy_position(position(0), position(1));
    inside = (xy_position.normSqr() <= radius*radius);

    // nearestPoint = point_on_surface_in_xy_plane + height_in_z
    Vec3 nearestPoint = normal*radius + Vec3(0,0,position(2));

    return nearestPoint;
}

bool ContactGeometry::Cylinder::Impl::intersectsRay
   (const Vec3& origin, const UnitVec3& direction,
    Real& distance, UnitVec3& normal) const
{
    // cylinder axis is z-axis, project to x-y plane
    const Vec3 xy_vec(direction(0),direction(1),0);
    const Real xy_vec_norm = xy_vec.norm();
    const UnitVec3 xy_direction(xy_vec/xy_vec_norm, true); // don't renormalize
    const Vec3 xy_origin(origin(0),origin(1),0);
    Real xy_distance;
    Real b = -~xy_direction*xy_origin;
    Real c = xy_origin.normSqr() - radius*radius;
    if (c > 0) {
        // Ray origin is outside cylinder.

        if (b <= 0)
          return false;  // Ray points away from axis of cylinder.
        Real d = b*b - c;
        if (d < 0)
          return false;
        Real root = std::sqrt(d);
        xy_distance = b - root;
      }
    else {
        // Ray origin is inside cylinder.

        Real d = b*b - c;
        if (d < 0)
          return false;
        xy_distance = b + std::sqrt(d);
      }
    distance = xy_distance/xy_vec_norm;
    normal = UnitVec3(xy_origin+xy_distance*xy_direction);
    return true;
}

void ContactGeometry::Cylinder::Impl::getBoundingSphere
    (Vec3& center, Real& radius) const {
    center = Vec3(0);
    radius = Infinity;
}

void ContactGeometry::Cylinder::Impl::
calcCurvature(const Vec3& point, Vec2& curvature, Rotation& orientation) const {

    // long axis (direction of min curvature) points in <0,0,1>
    orientation = Rotation(calcSurfaceUnitNormal(point), ZAxis, Vec3(0, 0, 1), YAxis);
    curvature[0] = 1/radius;
    curvature[1] = 0;

}

// Sample geodesic between two points P and Q on a cylinder analytically.
static void setGeodesicToHelicalArc(Real R, Real phiP, Real angle, Real m, Real c, Geodesic& geod)
{
   	// Clear current geodesic.
	geod.clear();

    const Real sqrt1m2 = sqrt(1+m*m);   // Avoid repeated calculation.
    const Real kappa = 1 / (R*(1+m*m)); // Curvature in tangent direction.
    const Real kb = 1 / (R*(1+1/(m*m))); // Curvature in binormal direction
                                         //   (slope is 1/m).
	const Real tau   = m*kappa;         // Torsion (signed).

	// Arc length of the helix. Always
	const Real L = R * sqrt1m2 * std::abs(angle);

	// Orientation of helix. 
	const Real orientation = angle < 0 ? Real(-1) : Real(1);

	// TODO: Make this generic, so long geodesics are sampled more than short ones.
    const int numGeodesicSamples = 12;
	const Real deltaPhi = std::abs(angle / Real(numGeodesicSamples-1));

    for (int i = 0; i < numGeodesicSamples; ++i)
	{
		// Watch out: Angle phi has an offset phiP
        Real phi = Real(i)*angle/Real(numGeodesicSamples-1) + phiP;
        const Real sphi = sin(phi), cphi = cos(phi);

		// Evaluate helix.
        Vec3	 p( R*cphi, R*sphi, R*m*(phi - phiP) + c);

        // We'll normalize so UnitVec3 doesn't have to do it.
		UnitVec3 t((orientation/sqrt1m2)*Vec3(-sphi, cphi, m), true);
		UnitVec3 n(Vec3(cphi, sphi, 0), true);

        // Though not needed, we use an orthogonalizing constructor for the rotation.
        geod.addFrenetFrame(Transform(Rotation(n, ZAxis, t, YAxis), p));

		// Current arc length s.
		Real s = R * sqrt1m2 * (Real(i)*deltaPhi);
        geod.addArcLength(s);
		geod.addCurvature(kappa);		

		// Solve the scalar Jacobi equation
		//
		//        j''(s) + K(s)*j(s) = 0 ,                                     (1)
		//
		// where K is the Gaussian curvature and (.)' := d(.)/ds denotes differentiation
		// with respect to the arc length s. Then, j is the directional sensitivity and
		// we obtain the corresponding variational vector field by multiplying b*j. For
		// a cylinder, K = 0 and the solution of equation (1) becomes
		//
		//        j  = s				                                       (2)
		//		  j' = 1 ,							                           (3)
		//
		// so the Jacobi field increases linearly in s.

		// Forward directional sensitivity from P to Q
		Vec2 jPQ(s, 1);
		geod.addDirectionalSensitivityPtoQ(jPQ);

		// Backwards directional sensitivity from Q to P
		Vec2 jQP(L-s, 1);
		geod.addDirectionalSensitivityQtoP(jQP);


        // TODO: positional sensitivity
        geod.addPositionalSensitivityPtoQ(Vec2(NaN));
        geod.addPositionalSensitivityQtoP(Vec2(NaN));
    }

	// Only compute torsion and binormal curvature at the end points.
	geod.setTorsionAtP(tau); geod.setTorsionAtQ(tau);
    geod.setBinormalCurvatureAtP(kb); geod.setBinormalCurvatureAtQ(kb);

    geod.setIsConvex(true); // Curve on cylinder is always convex.

    geod.setIsShortest(false); // TODO
    geod.setAchievedAccuracy(SignificantReal); // TODO: accuracy of length?
//    geod.setInitialStepSizeHint(integ.getActualInitialStepSizeTaken()); // TODO
}

// Compute geodesic between two points P and Q on a cylinder analytically. Since a geodesic on a
// cylinder is a helix it is parameterized by
//
//			       [ R * cos(phi)    ]
//        p(phi) = [ R * sin(phi)    ]
//				   [ R * m * phi + c ]
//
// where R is the radius of the cylinder, phi parameterizes the opening angle of the helix, m is  
// the slope and c is an offset. We define the geodesic from P to Q, hence c = Pz.
void ContactGeometry::Cylinder::Impl::
calcGeodesicAnalytical(const Vec3& xP, const Vec3& xQ,
                       const Vec3& tPhint, const Vec3& tQhint,
                       Geodesic& geod) const
{
	// Compute angle between P and Q. Save both the positive (right handed) and the negative
	// (left handed) angle.
	Real phiP = atan2(xP[1], xP[0]);
	Real phiQ = atan2(xQ[1], xQ[0]);

	Real temp = phiQ - phiP;
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

	// Compute "moment" of tPhint at P and tQhint at Q around z-Axis.
	// Make sure tPhint and tQhint are unit vectors, otherwise moments are scaled.
	Real MP = xP[0]*tPhint[1] - xP[1]*tPhint[0];
	Real MQ = xQ[0]*tQhint[1] - xQ[1]*tQhint[0];

	// Average moment.
	Real M = (MP + MQ) / 2;

	// Decide whether helix is right or left handed. The sign of angle stores the 
	// information about the orientation (right handed, if positive)
	Real angle;
	if (M >= 0)	{
		angle = angleRightHanded;
	}

	else {
		angle = angleLeftHanded;
	}

	// Offset and slope.
	Real c =  xP[2];
	Real m = (xQ[2] - xP[2]) / (angle * radius);

	setGeodesicToHelicalArc(radius, phiP, angle, m, c, geod);
}

void ContactGeometry::Cylinder::Impl::shootGeodesicInDirectionUntilLengthReachedAnalytical(const Vec3& xP, const UnitVec3& tP,
        const Real& terminatingLength, const GeodesicOptions& options, Geodesic& geod) const {

    //TODO for Andreas :)
}

void ContactGeometry::Cylinder::Impl::shootGeodesicInDirectionUntilPlaneHitAnalytical(const Vec3& xP, const UnitVec3& tP,
        const Plane& terminatingPlane, const GeodesicOptions& options,
        Geodesic& geod) const {

    //TODO for Andreas :)
}

Real CylinderImplicitFunction::
calcValue(const Vector& x) const {
    return 1-(x[0]*x[0]+x[1]*x[1])/square(ownerp->getRadius());
}

Real CylinderImplicitFunction::
calcDerivative(const Array_<int>& derivComponents, const Vector& x) const {
    if (derivComponents.size() == 1 && derivComponents[0] < 2)
        return -2*x[derivComponents[0]]/square(ownerp->getRadius());
    if (derivComponents.size() == 2 &&
        derivComponents[0] == derivComponents[1] &&
        derivComponents[0] < 2 )
        return -2/square(ownerp->getRadius());
    return 0;
}
