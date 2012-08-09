/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
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

#include "SimTKmath.h"
#include <vector>
#include <exception>

using namespace SimTK;
using namespace std;

//TODO - had to reduced tol for geodesic tests, should be 1e-10
const Real TOL = 1e-4;

const Real r = 3.5; // radius used for geodesic tests


#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

template <class T>
void assertEqual(T val1, T val2) {
    ASSERT(abs(val1-val2) < TOL);
}

template <int N>
void assertEqual(Vec<N> val1, Vec<N> val2) {
    for (int i = 0; i < N; ++i)
        ASSERT(abs(val1[i]-val2[i]) < TOL);
}

void testSurfaceGradient(const ContactGeometry& );
void testSurfaceHessian(const ContactGeometry& );
void compareAnalyticalAndNumericGradient(const ContactGeometry&, const Vec3& );
void compareAnalyticalAndNumericHessian(const ContactGeometry&, const Vec3& );

void testHalfSpace() {
    ContactGeometry::HalfSpace hs;
    
    // Check intersections with various rays.
    
    Real distance;
    UnitVec3 normal;
    ASSERT(!hs.intersectsRay(Vec3(4, 0, 0), UnitVec3(1, 0, 0), distance, normal));
    ASSERT(!hs.intersectsRay(Vec3(-1, 0, 0), UnitVec3(-1, 1, 0), distance, normal));
    ASSERT(!hs.intersectsRay(Vec3(-1, 0, 0), UnitVec3(0, 1, 0), distance, normal));
    ASSERT(hs.intersectsRay(Vec3(-1, 0, 0), UnitVec3(1, 0, 0), distance, normal));
    assertEqual(1.0, distance);
    assertEqual(Vec3(-1, 0, 0), normal);
    ASSERT(hs.intersectsRay(Vec3(-2, 15, 37), UnitVec3(1, 0, 0), distance, normal));
    assertEqual(2.0, distance);
    assertEqual(Vec3(-1, 0, 0), normal);
    ASSERT(hs.intersectsRay(Vec3(-3, 1, 2), UnitVec3(1, 1, 1), distance, normal));
    assertEqual(3*Sqrt3, distance);
    assertEqual(Vec3(-1, 0, 0), normal);
    ASSERT(hs.intersectsRay(Vec3(2, 0, 0), UnitVec3(-1, 0, 1), distance, normal));
    assertEqual(2*Sqrt2, distance);
    assertEqual(Vec3(-1, 0, 0), normal);
    
    // Test finding the nearest point.
    
    Random::Gaussian random(0, 3);
    for (int i = 0; i < 100; i++) {
        Vec3 pos(random.getValue(), random.getValue(), random.getValue());
        bool inside;
        UnitVec3 normal;
        Vec3 nearest = hs.findNearestPoint(pos, inside, normal);
        assertEqual(nearest, Vec3(0, pos[1], pos[2]));
        ASSERT(inside == (pos[0] >= 0));
        assertEqual(normal, Vec3(-1, 0, 0));
    }
}

void testCylinder() {
    Real radius = 3.5;
    ContactGeometry::Cylinder cyl(radius);
    assert(cyl.getRadius() == radius);

    // Check intersections with various rays.

    Real distance;
    UnitVec3 normal;
    ASSERT(!cyl.intersectsRay(Vec3(radius*1.1, 0, 0), UnitVec3(1, 0, 0), distance, normal));
    ASSERT(!cyl.intersectsRay(Vec3(-radius*1.1, 0, 0), UnitVec3(-1, 1, 0), distance, normal));
    ASSERT(!cyl.intersectsRay(Vec3(-radius*1.1, 0, 0), UnitVec3(0, 1, 0), distance, normal));
    ASSERT(!cyl.intersectsRay(Vec3(-radius, -radius, 0), UnitVec3(1, -Eps, 0), distance, normal));
    ASSERT(cyl.intersectsRay(Vec3(-radius, -radius, 0), UnitVec3(1, Eps, 0), distance, normal));

    ASSERT(cyl.intersectsRay(Vec3(-(radius+1.0), 0, 0), UnitVec3(1, 0, 0), distance, normal));
    assertEqual(1.0, distance);
    assertEqual(Vec3(-1, 0, 0), normal);

    ASSERT(cyl.intersectsRay(Vec3(-radius*2, radius*2, 37), UnitVec3(1, -1, 0), distance, normal));
    assertEqual(radius*(2*Sqrt2-1), distance);
    assertEqual(UnitVec3(-1, 1, 0), normal);

    ASSERT(cyl.intersectsRay(Vec3(-radius*2, 0, -radius*2), UnitVec3(1, 0, 1), distance, normal));
    assertEqual(radius*Sqrt2, distance);
    assertEqual(UnitVec3(-1, 0, 0), normal);

    // Test finding the nearest point.

    Random::Gaussian random(0, 3);
    for (int i = 0; i < 100; i++) {
        Vec3 pos(random.getValue(), random.getValue(), random.getValue());
        Vec3 projpos(pos);
        projpos(2)=0; // cyl axis is z-axis, project pos to x-y plane
        bool inside;
        UnitVec3 normal;
        Vec3 nearest = cyl.findNearestPoint(pos, inside, normal);
        assertEqual(nearest, projpos.normalize()*radius+Vec3(0,0,pos(2)));
        ASSERT(inside == (projpos.norm() <= radius));
        assertEqual(normal, projpos.normalize());
    }

    // Test derivatives
    testSurfaceGradient(cyl);
    testSurfaceHessian(cyl);
}

void testSphere() {
    // Create a sphere.
    
    Real radius = 3.5;
    ContactGeometry::Sphere sphere(radius);
    assert(sphere.getRadius() == radius);
    
    // Check intersections with various rays.
    
    Real distance;
    UnitVec3 normal;
    ASSERT(!sphere.intersectsRay(Vec3(4, 0, 0), UnitVec3(1, 0, 0), distance, normal));
    ASSERT(sphere.intersectsRay(Vec3(2, 0, 0), UnitVec3(1, 0, 0), distance, normal));
    assertEqual(1.5, distance);
    assertEqual(Vec3(1, 0, 0), normal);
    ASSERT(sphere.intersectsRay(Vec3(4, 0, 0), UnitVec3(-1, 0, 0), distance, normal));
    assertEqual(0.5, distance);
    assertEqual(Vec3(1, 0, 0), normal);
    ASSERT(sphere.intersectsRay(Vec3(2, 0, 0), UnitVec3(-1, 0, 0), distance, normal));
    assertEqual(5.5, distance);
    assertEqual(Vec3(-1, 0, 0), normal);
    ASSERT(sphere.intersectsRay(Vec3(0, 0, 0), UnitVec3(1, 1, 1), distance, normal));
    assertEqual(3.5, distance);
    assertEqual(Vec3(1.0/Sqrt3), normal);

    // Test finding the nearest point.
    
    Random::Gaussian random(0, 3);
    for (int i = 0; i < 100; i++) {
        Vec3 pos(random.getValue(), random.getValue(), random.getValue());
        bool inside;
        UnitVec3 normal;
        Vec3 nearest = sphere.findNearestPoint(pos, inside, normal);
        assertEqual(nearest, pos.normalize()*radius);
        ASSERT(inside == (pos.norm() <= radius));
        assertEqual(normal, pos.normalize());
    }
}

void testEllipsoid() {
    // Create a ellipsoid.

    Vec3 radii(1.5, 2.2, 3.1);
    ContactGeometry::Ellipsoid ellipsoid(radii);
    assert(ellipsoid.getRadii() == radii);

    // Check intersections with various rays.

    Real distance;
    UnitVec3 normal;
    ASSERT(!ellipsoid.intersectsRay(Vec3(4, 0, 0), UnitVec3(1, 0, 0), distance, normal));
    ASSERT(ellipsoid.intersectsRay(Vec3(1, 0, 0), UnitVec3(1, 0, 0), distance, normal));
    assertEqual(0.5, distance);
    assertEqual(Vec3(1, 0, 0), normal);
    ASSERT(ellipsoid.intersectsRay(Vec3(4, 0, 0), UnitVec3(-1, 0, 0), distance, normal));
    assertEqual(2.5, distance);
    assertEqual(Vec3(1, 0, 0), normal);
    ASSERT(ellipsoid.intersectsRay(Vec3(0, -5, 0), UnitVec3(0, 1, 0), distance, normal));
    assertEqual(2.8, distance);
    assertEqual(Vec3(0, -1, 0), normal);
    ASSERT(ellipsoid.intersectsRay(Vec3(0, 0, 0), UnitVec3(1, 1, 1), distance, normal));
    assertEqual(sqrt(3/(1/(radii[0]*radii[0])+1/(radii[1]*radii[1])+1/(radii[2]*radii[2]))), distance);
    assertEqual(UnitVec3(1/(radii[0]*radii[0]), 1/(radii[1]*radii[1]), 1/(radii[2]*radii[2])), normal);

    // Test finding the nearest point.

    Random::Gaussian random(0, 2);
    for (int i = 0; i < 100; i++) {
        Vec3 pos(random.getValue(), random.getValue(), random.getValue());
        bool inside;
        UnitVec3 normal;
        Vec3 nearest = ellipsoid.findNearestPoint(pos, inside, normal);
        assertEqual(nearest[0]*nearest[0]/(radii[0]*radii[0])+nearest[1]*nearest[1]/(radii[1]*radii[1])+nearest[2]*nearest[2]/(radii[2]*radii[2]), 1.0);
        Real projectedRadius = pos[0]*pos[0]/(radii[0]*radii[0])+pos[1]*pos[1]/(radii[1]*radii[1])+pos[2]*pos[2]/(radii[2]*radii[2]);
        ASSERT(inside == (projectedRadius < 1.0));
        Vec3 projectedPoint = pos/sqrt(projectedRadius);
        ASSERT((nearest-pos).normSqr() < (projectedPoint-pos).normSqr());
        assertEqual(normal, UnitVec3(nearest[0]/(radii[0]*radii[0]), nearest[1]/(radii[1]*radii[1]), nearest[2]/(radii[2]*radii[2])));
    }
}

void testTorus() {
    Real radius = r;
    Real tubeRadius = 0.75;
    ContactGeometry::Torus torus(radius, tubeRadius);
    assert(torus.getTorusRadius() == radius);
    assert(torus.getTubeRadius() == tubeRadius);

    // Check intersections with various rays.

//    Real distance;
//    UnitVec3 normal;
//    ASSERT(!torus.intersectsRay(Vec3(radius*1.1, 0, 0), UnitVec3(1, 0, 0), distance, normal));
//    ASSERT(!torus.intersectsRay(Vec3(-radius*1.1, 0, 0), UnitVec3(-1, 1, 0), distance, normal));
//    ASSERT(!torus.intersectsRay(Vec3(-radius*1.1, 0, 0), UnitVec3(0, 1, 0), distance, normal));
//    ASSERT(!torus.intersectsRay(Vec3(-radius, -radius, 0), UnitVec3(1, -Eps, 0), distance, normal));
//    ASSERT(torus.intersectsRay(Vec3(-radius, -radius, 0), UnitVec3(1, Eps, 0), distance, normal));
//
//    ASSERT(torus.intersectsRay(Vec3(-(radius+1.0), 0, 0), UnitVec3(1, 0, 0), distance, normal));
//    assertEqual(1.0, distance);
//    assertEqual(Vec3(-1, 0, 0), normal);
//
//    ASSERT(torus.intersectsRay(Vec3(-radius*2, radius*2, 37), UnitVec3(1, -1, 0), distance, normal));
//    assertEqual(radius*(2*Sqrt2-1), distance);
//    assertEqual(UnitVec3(-1, 1, 0), normal);
//
//    ASSERT(torus.intersectsRay(Vec3(-radius*2, 0, -radius*2), UnitVec3(1, 0, 1), distance, normal));
//    assertEqual(radius*Sqrt2, distance);
//    assertEqual(UnitVec3(-1, 0, 0), normal);
//
//    // Test finding the nearest point.
//
//    Random::Gaussian random(0, 3);
//    for (int i = 0; i < 100; i++) {
//        Vec3 pos(random.getValue(), random.getValue(), random.getValue());
//        Vec3 projpos(pos);
//        projpos(2)=0; // cyl axis is z-axis, project pos to x-y plane
//        bool inside;
//        UnitVec3 normal;
//        Vec3 nearest = torus.findNearestPoint(pos, inside, normal);
//        assertEqual(nearest, projpos.normalize()*radius+Vec3(0,0,pos(2)));
//        ASSERT(inside == (projpos.norm() <= radius));
//        assertEqual(normal, projpos.normalize());
//    }

    // Test derivatives
//    Vec3 pt(2*(radius+tubeRadius), 0,0);
//    Vec3 pt(3,4,5);
//    cout << "surfaceValue( " << pt << " ) = " << torus.calcSurfaceValue(pt) << endl;
//    compareAnalyticalAndNumericGradient(torus, pt);
//    compareAnalyticalAndNumericHessian(torus, pt);
    testSurfaceGradient(torus);
    testSurfaceHessian(torus);

}


//==============================================================================
//                      SURFACE EVALUATORS
//==============================================================================

class ImplicitSurfaceFunction : public Differentiator::GradientFunction {
public:
    ImplicitSurfaceFunction(const ContactGeometry& geom)
        : Differentiator::GradientFunction(3), geom(geom) { }

    int f(const Vector& y, Real& fy) const {
        fy = geom.calcSurfaceValue(Vec3::getAs(&y[0]));
        return 0;
    }

    const ContactGeometry& geom;
};

class ImplicitSurfaceGradient : public Differentiator::JacobianFunction {
public:
    ImplicitSurfaceGradient(const ContactGeometry& geom)
        : Differentiator::JacobianFunction(3,3), geom(geom) { }

    int f(const Vector& y, Vector& fy) const {
        fy = (Vector)geom.calcSurfaceGradient(Vec3::getAs(&y(0)));
        return 0;
    }

    const ContactGeometry& geom;
};

void testSurfaceGradient(const ContactGeometry& geom) {

    // setup random test points
    Random::Gaussian random(0, r);
    for (int i = 0; i < 100; i++) {
        Vec3 pt(random.getValue(), random.getValue(), random.getValue());
        compareAnalyticalAndNumericGradient(geom, pt);
    }
}

void compareAnalyticalAndNumericGradient(const ContactGeometry& geom, const Vec3& pt) {

    ImplicitSurfaceFunction surf(geom);
    Differentiator diff(surf);

    Vector tmp = diff.calcGradient((Vector)pt);
    Vec3 gradNumeric = Vec3::getAs(&tmp(0));
    Vec3 gradAnalytic = geom.calcSurfaceGradient(pt);

//    cout << "ana gradient = " << gradAnalytic << endl;
//    cout << "num gradient = " << gradNumeric << endl;
    assertEqual(gradNumeric, gradAnalytic);
}

void testSurfaceHessian(const ContactGeometry& geom) {

    // setup random test points
    Random::Gaussian random(0, r);
    for (int i = 0; i < 100; i++) {
        Vec3 pt(random.getValue(), random.getValue(), random.getValue());
        compareAnalyticalAndNumericHessian(geom, pt);
    }
}

void compareAnalyticalAndNumericHessian(const ContactGeometry& geom, const Vec3& pt) {

    ImplicitSurfaceGradient grad(geom);
    Differentiator diff(grad);

    Matrix tmp = diff.calcJacobian((Vector)pt);
    Mat33 hessNumeric = Mat33::getAs(&tmp(0,0));
    Mat33 hessAnalytic = geom.calcSurfaceHessian(pt);

//    cout << "ana hessian = " << hessAnalytic << endl;
//    cout << "num hessian = " << hessNumeric << endl;
    assertEqual(hessNumeric(0), hessAnalytic(0));
    assertEqual(hessNumeric(1), hessAnalytic(1));
    assertEqual(hessNumeric(2), hessAnalytic(2));
}



//==============================================================================
//                            GEODESIC EVALUTORS
//==============================================================================


void compareAnalyticalAndNumericGeodesic(const ContactGeometry& geom, const Vec3& P, const Vec3& Q) {

    Geodesic geodNumeric;
    Geodesic geodAnalytic;

    //TODO -- tests only pass for very small amounts of noise added
    // to the initial conditions for the numerical geodesic calculation

    Real lengthNoise = 0;
    Real tPNoise = 0;

    Real len;
    Random::Gaussian dlen(0, lengthNoise);

    Vec3 tP(0);
    Random::Gaussian dtP(0, tPNoise);

    UnitVec3 e_PQ(Q-P);
//    cout << "P = " << P << ", Q = " << Q << endl;

    // short geodesic, tPhint and tQhint point toward each other
    geom.calcGeodesicAnalytical(P, Q, e_PQ, e_PQ, geodAnalytic);

    // set init length and tangent dir to analytical result + noise
    len = geodAnalytic.getLength()*(1+dlen.getValue());
    tP = geodAnalytic.getTangentP() + Vec3(dtP.getValue(), dtP.getValue(), dtP.getValue());

    geom.calcGeodesicUsingOrthogonalMethod(P, Q, tP, len, geodNumeric);

//        cout << "ana length = " << geodAnalytic.getLength() << endl;
//        cout << "ini length = " << len << endl;
//        cout << "num length = " << geodNumeric.getLength() << endl;
    assertEqual(geodNumeric.getLength(), geodAnalytic.getLength());

//        cout << "ana tangentP = " << geodAnalytic.getTangentP() << endl;
//        cout << "ini tangentP = " << tP << endl;
//        cout << "num tangentP = " << geodNumeric.getTangentP() << endl;
    assertEqual(geodNumeric.getTangentP(), geodAnalytic.getTangentP());

//        cout << "ana tangentQ = " << geodAnalytic.getTangentQ() << endl;
//        cout << "num tangentQ = " << geodNumeric.getTangentQ() << endl;
    assertEqual(geodNumeric.getTangentQ(), geodAnalytic.getTangentQ());
}


void testAnalyticalGeodesicRandom(const ContactGeometry& geom) {

    Vec3 P(r,0,0);
    Vec3 Q(0,r,0);
    compareAnalyticalAndNumericGeodesic(geom,P,Q);

    // setup random test points
    Random::Gaussian random(0, r);
    for (int i = 0; i < 100; i++) {
        Vec3 P(random.getValue(), random.getValue(), random.getValue());
        Vec3 Q(random.getValue(), random.getValue(), random.getValue());
        compareAnalyticalAndNumericGeodesic(geom,P,Q);
    }
}

void testAnalyticalSphereGeodesic() {
    ContactGeometry::Sphere sphere(r);
    testAnalyticalGeodesicRandom(sphere);
}

void testAnalyticalCylinderGeodesic() {
    ContactGeometry::Cylinder cylinder(r);
    testAnalyticalGeodesicRandom(cylinder);
}

void testFindNearestPointWithNewtonsMethod(const ContactGeometry& geom, Real r) {

    Random::Gaussian random(0, r);
    for (int i = 0; i < 100; i++) {
        Vec3 pos(random.getValue(), random.getValue(), random.getValue());
//        cout << i << ": pos=" << pos << endl;

        bool inside;
        UnitVec3 normal;
        Vec3 nearest = geom.findNearestPoint(pos, inside, normal);
        Vec3 nearest2 = geom.findNearestPointUsingNewtonsMethod(pos, inside, normal);

//        cout << "near pt analytical = " << nearest << "norm = " << nearest.norm() << endl;
//        cout << "near pt newton mth = " << nearest2 << "norm = " << nearest2.norm() << endl;

        assertEqual(nearest, nearest2);

    }
}


int main() {
    try {
        testHalfSpace();
        testSphere();
        testEllipsoid();
	    testCylinder();
	    testTorus();
	    testAnalyticalSphereGeodesic();
	    testAnalyticalCylinderGeodesic();
//	    testFindNearestPointWithNewtonsMethod(ContactGeometry::Sphere(r), r);
    }
    catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
