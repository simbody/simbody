/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
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

#include "SimTKsimbody.h"

using namespace SimTK;
using namespace std;

const Real TOL = 1e-10;

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

void testContainsPoint() {
    OrientedBoundingBox box(Vec3(0), Vec3(1, 2, 3));
    ASSERT(box.containsPoint(Vec3(0)));
    ASSERT(box.containsPoint(Vec3(1, 2, 3)));
    ASSERT(box.containsPoint(Vec3(0.5, 1.0, 1.5)));
    ASSERT(!box.containsPoint(Vec3(-0.5, 1.0, 1.5)));
    ASSERT(!box.containsPoint(Vec3(0.5, -1.0, 1.5)));
    ASSERT(!box.containsPoint(Vec3(0.5, 1.0, -1.5)));
    ASSERT(!box.containsPoint(Vec3(1.01, 1.0, 1.5)));
    ASSERT(!box.containsPoint(Vec3(0.5, 2.01, 1.5)));
    ASSERT(!box.containsPoint(Vec3(0.5, 1.0, 3.01)));
    box = OrientedBoundingBox(Vec3(3, 2, 1), Vec3(1, 2, 3));
    ASSERT(box.containsPoint(Vec3(3, 2, 1)));
    ASSERT(box.containsPoint(Vec3(4, 4, 4)));
    ASSERT(box.containsPoint(Vec3(3.5, 3.0, 2.5)));
    ASSERT(!box.containsPoint(Vec3(2.5, 3.0, 2.5)));
    ASSERT(!box.containsPoint(Vec3(4.5, 1.0, 2.5)));
    ASSERT(!box.containsPoint(Vec3(4.5, 3.0, -0.5)));
    ASSERT(!box.containsPoint(Vec3(5.01, 3.0, 2.5)));
    ASSERT(!box.containsPoint(Vec3(3.5, 4.01, 2.5)));
    ASSERT(!box.containsPoint(Vec3(3.5, 3.0, 4.01)));
    box = OrientedBoundingBox(Rotation(0.25*Pi, ZAxis), Vec3(1, 1, 1));
    ASSERT(box.containsPoint(Vec3(0)));
    ASSERT(box.containsPoint(Vec3(0, Sqrt2-1e-10, 0)));
    ASSERT(!box.containsPoint(Vec3(0, Sqrt2+1e-10, 0)));
    ASSERT(box.containsPoint(Vec3(0, Sqrt2-1e-10, 1e-10)));
    ASSERT(!box.containsPoint(Vec3(0, Sqrt2-1e-10, -1e-10)));
}

void verifyCorners(Vec3 expected[8], Vec3 found[8]) {
    for (int i = 0; i < 8; i++) {
        bool match = false;
        for (int j = 0; j < 8 && !match; j++) {
            if (abs(expected[i][0]-found[j][0]) < TOL && abs(expected[i][1]-found[j][1]) < TOL && abs(expected[i][2]-found[j][2]) < TOL) {
                match = true;
                break;
            }
        }
        ASSERT(match);
    }
}

void testGetCorners() {
    Vec3 corners[8];
    OrientedBoundingBox(Vec3(0), Vec3(1, 2, 3)).getCorners(corners);
    {
        Vec3 expected[] = {Vec3(0), Vec3(1, 0, 0), Vec3(0, 2, 0), Vec3(1, 2, 0), Vec3(0, 0, 3), Vec3(1, 0, 3), Vec3(0, 2, 3), Vec3(1, 2, 3)};
        verifyCorners(expected, corners);
    }
    OrientedBoundingBox(Vec3(3, 2, 1), Vec3(1, 2, 3)).getCorners(corners);
    {
        Vec3 expected[] = {Vec3(3, 2, 1), Vec3(4, 2, 1), Vec3(3, 4, 1), Vec3(4, 4, 1), Vec3(3, 2, 4), Vec3(4, 2, 4), Vec3(3, 4, 4), Vec3(4, 4, 4)};
        verifyCorners(expected, corners);
    }
    OrientedBoundingBox(Rotation(0.25*Pi, ZAxis), Vec3(1, 1, 1)).getCorners(corners);
    {
        Real d = 0.5*Sqrt2;
        Vec3 expected[] = {Vec3(0), Vec3(-d, d, 0), Vec3(d, d, 0), Vec3(0, Sqrt2, 0), Vec3(0, 0, 1), Vec3(-d, d, 1), Vec3(d, d, 1), Vec3(0, Sqrt2, 1)};
        verifyCorners(expected, corners);
    }
}

void verifyBoxIntersection(bool shouldIntersect, OrientedBoundingBox box1, OrientedBoundingBox box2) {
    if (box1.intersectsBox(box2) != shouldIntersect)
        std::cout << "a"<< std::endl;
    if (box2.intersectsBox(box1) != shouldIntersect)
        std::cout << "b"<< std::endl;
    ASSERT(box1.intersectsBox(box2) == shouldIntersect);
    ASSERT(box2.intersectsBox(box1) == shouldIntersect);
}

void testIntersectsBox() {
    // Try boxes with identical orientations.
    
    verifyBoxIntersection(true, OrientedBoundingBox(Vec3(0), Vec3(1, 2, 3)), OrientedBoundingBox(Vec3(0), Vec3(3, 2, 1)));
    verifyBoxIntersection(true, OrientedBoundingBox(Vec3(0), Vec3(1, 2, 3)), OrientedBoundingBox(Vec3(1-1e-10, 2-1e-10, 3-1e-10), Vec3(3, 2, 1)));
    verifyBoxIntersection(false, OrientedBoundingBox(Vec3(0), Vec3(1, 2, 3)), OrientedBoundingBox(Vec3(1+1e-10, 2-1e-10, 3-1e-10), Vec3(3, 2, 1)));
    verifyBoxIntersection(false, OrientedBoundingBox(Vec3(0), Vec3(1, 2, 3)), OrientedBoundingBox(Vec3(1-1e-10, 2+1e-10, 3-1e-10), Vec3(3, 2, 1)));
    verifyBoxIntersection(false, OrientedBoundingBox(Vec3(0), Vec3(1, 2, 3)), OrientedBoundingBox(Vec3(1-1e-10, 2-1e-10, 3+1e-10), Vec3(3, 2, 1)));
    verifyBoxIntersection(true, OrientedBoundingBox(Vec3(0), Vec3(1, 2, 3)), OrientedBoundingBox(Vec3(-3+1e-10, -2+1e-10, -1+1e-10), Vec3(3, 2, 1)));
    verifyBoxIntersection(false, OrientedBoundingBox(Vec3(0), Vec3(1, 2, 3)), OrientedBoundingBox(Vec3(-3-1e-10, -2+1e-10, -1+1e-10), Vec3(3, 2, 1)));
    verifyBoxIntersection(false, OrientedBoundingBox(Vec3(0), Vec3(1, 2, 3)), OrientedBoundingBox(Vec3(-3+1e-10, -2-1e-10, -1+1e-10), Vec3(3, 2, 1)));
    verifyBoxIntersection(false, OrientedBoundingBox(Vec3(0), Vec3(1, 2, 3)), OrientedBoundingBox(Vec3(-3+1e-10, -2+1e-10, -1-1e-10), Vec3(3, 2, 1)));
    
    // Try some rotations by 90 degrees.
    
    verifyBoxIntersection(true, OrientedBoundingBox(Vec3(-0.5, 0, 0), Vec3(1, 2, 3)), OrientedBoundingBox(Rotation(0.5*Pi, ZAxis), Vec3(3, 2, 1)));
    verifyBoxIntersection(true, OrientedBoundingBox(Vec3(-0.5, -2+1e-10, 0), Vec3(1, 2, 3)), OrientedBoundingBox(Rotation(0.5*Pi, ZAxis), Vec3(3, 2, 1)));
    verifyBoxIntersection(false, OrientedBoundingBox(Vec3(-0.5, -2-1e-10, 0), Vec3(1, 2, 3)), OrientedBoundingBox(Rotation(0.5*Pi, ZAxis), Vec3(3, 2, 1)));
    verifyBoxIntersection(true, OrientedBoundingBox(Vec3(-0.5, 3-1e-10, 0), Vec3(1, 2, 3)), OrientedBoundingBox(Rotation(0.5*Pi, ZAxis), Vec3(3, 2, 1)));
    verifyBoxIntersection(false, OrientedBoundingBox(Vec3(-0.5, 3+1e-10, 0), Vec3(1, 2, 3)), OrientedBoundingBox(Rotation(0.5*Pi, ZAxis), Vec3(3, 2, 1)));
    verifyBoxIntersection(true, OrientedBoundingBox(Rotation(0.5*Pi, ZAxis), Vec3(1, 2, 3)), OrientedBoundingBox(Rotation(0.5*Pi, ZAxis), Vec3(3, 2, 1)));
    
    // Try rotations by 45 degrees.

    const Real d = 0.5*Sqrt2;
    verifyBoxIntersection(true, OrientedBoundingBox(Vec3(-0.5, 0, 0), Vec3(1, 1, 1)), OrientedBoundingBox(Rotation(0.25*Pi, ZAxis), Vec3(1, 1, 1)));
    verifyBoxIntersection(true, OrientedBoundingBox(Vec3(-0.5, -1+1e-10, 0), Vec3(1, 1, 1)), OrientedBoundingBox(Rotation(0.25*Pi, ZAxis), Vec3(1, 1, 1)));
    verifyBoxIntersection(false, OrientedBoundingBox(Vec3(-0.5, -1-1e-10, 0), Vec3(1, 1, 1)), OrientedBoundingBox(Rotation(0.25*Pi, ZAxis), Vec3(1, 1, 1)));
    verifyBoxIntersection(true, OrientedBoundingBox(Vec3(-0.5, Sqrt2-1e-10, 0), Vec3(1, 1, 1)), OrientedBoundingBox(Rotation(0.25*Pi, ZAxis), Vec3(1, 1, 1)));
    verifyBoxIntersection(false, OrientedBoundingBox(Vec3(-0.5, Sqrt2+1e-10, 0), Vec3(1, 1, 1)), OrientedBoundingBox(Rotation(0.25*Pi, ZAxis), Vec3(1, 1, 1)));
    verifyBoxIntersection(true, OrientedBoundingBox(Vec3(d-1e-10, 0, 0), Vec3(1, 1, 1)), OrientedBoundingBox(Rotation(0.25*Pi, ZAxis), Vec3(1, 1, 1)));
    verifyBoxIntersection(false, OrientedBoundingBox(Vec3(d+1e-10, 0, 0), Vec3(1, 1, 1)), OrientedBoundingBox(Rotation(0.25*Pi, ZAxis), Vec3(1, 1, 1)));
    verifyBoxIntersection(true, OrientedBoundingBox(Vec3(-1-d+1e-10, 0, 0), Vec3(1, 1, 1)), OrientedBoundingBox(Rotation(0.25*Pi, ZAxis), Vec3(1, 1, 1)));
    verifyBoxIntersection(false, OrientedBoundingBox(Vec3(-1-d-1e-10, 0, 0), Vec3(1, 1, 1)), OrientedBoundingBox(Rotation(0.25*Pi, ZAxis), Vec3(1, 1, 1)));
    verifyBoxIntersection(true, OrientedBoundingBox(Vec3(0, 0, d-1e-10), Vec3(1, 1, 1)), OrientedBoundingBox(Rotation(-0.25*Pi, XAxis), Vec3(1, 1, 1)));
    verifyBoxIntersection(false, OrientedBoundingBox(Vec3(0, 0, d+1e-10), Vec3(1, 1, 1)), OrientedBoundingBox(Rotation(-0.25*Pi, XAxis), Vec3(1, 1, 1)));
    verifyBoxIntersection(true, OrientedBoundingBox(Vec3(0, 0, -1-d+1e-10), Vec3(1, 1, 1)), OrientedBoundingBox(Rotation(-0.25*Pi, XAxis), Vec3(1, 1, 1)));
    verifyBoxIntersection(false, OrientedBoundingBox(Vec3(0, 0, -1-d-1e-10), Vec3(1, 1, 1)), OrientedBoundingBox(Rotation(-0.25*Pi, XAxis), Vec3(1, 1, 1)));
    
    // The following case detects a bug that came up in a different test case.
    
    Rotation r;
    r.setRotationToBodyFixedXYZ(Vec3(Pi/2, 0, -1.1));;
    verifyBoxIntersection(true, OrientedBoundingBox(Transform(r, Vec3(-0.95, 1.1, 2.5)), Vec3(2e-10, 1.118, 1)), OrientedBoundingBox(Vec3(0, -50, -50), Vec3(100, 100, 100)));
}

void verifyRayIntersection(const OrientedBoundingBox& box, const Vec3& origin, const UnitVec3& direction, bool shouldIntersect, Real expectedDistance) {
    Real distance;
    ASSERT(shouldIntersect == box.intersectsRay(origin, direction, distance));
    assertEqual(expectedDistance, distance);
}

void testIntersectsRay() {
    // Try rays starting inside the box.
    
    verifyRayIntersection(OrientedBoundingBox(Transform(), Vec3(1, 2, 3)), Vec3(0.5, 0.5, 0.5), UnitVec3(1, 0, 0), true, 0);
    verifyRayIntersection(OrientedBoundingBox(Transform(), Vec3(1, 2, 3)), Vec3(0.5, 0.5, 0.5), UnitVec3(0, 1, 0), true, 0);
    verifyRayIntersection(OrientedBoundingBox(Transform(), Vec3(1, 2, 3)), Vec3(0.5, 0.5, 0.5), UnitVec3(0, 0, 1), true, 0);
    
    // Try rays that hit it on various sides.
    
    verifyRayIntersection(OrientedBoundingBox(Transform(), Vec3(1, 2, 3)), Vec3(-1.5, 0.5, 0.5), UnitVec3(1, 0, 0), true, 1.5);
    verifyRayIntersection(OrientedBoundingBox(Transform(), Vec3(1, 2, 3)), Vec3(2.5, 0.5, 0.5), UnitVec3(-1, 0, 0), true, 1.5);
    verifyRayIntersection(OrientedBoundingBox(Transform(), Vec3(1, 2, 3)), Vec3(0.5, -0.5, 0.5), UnitVec3(0, 1, 0), true, 0.5);
    verifyRayIntersection(OrientedBoundingBox(Transform(), Vec3(1, 2, 3)), Vec3(0.5, 2.5, 0.5), UnitVec3(0, -1, 0), true, 0.5);
    verifyRayIntersection(OrientedBoundingBox(Transform(), Vec3(1, 2, 3)), Vec3(0.5, 0.5, -1.0), UnitVec3(0, 0, 1), true, 1.0);
    verifyRayIntersection(OrientedBoundingBox(Transform(), Vec3(1, 2, 3)), Vec3(0.5, 0.5, 4.0), UnitVec3(0, 0, -1), true, 1.0);
    
    // Try a box at an angle.

    verifyRayIntersection(OrientedBoundingBox(Rotation(-0.25*Pi, ZAxis), Vec3(2, 2, 2)), Vec3(0, 0, 0.5), UnitVec3(1, 0, 0), true, 0);
    verifyRayIntersection(OrientedBoundingBox(Rotation(-0.25*Pi, ZAxis), Vec3(2, 2, 2)), Vec3(-1, 0, 0.5), UnitVec3(1, 0, 0), true, 1.0);
}

void testCreateFromPoints() {
    Random::Uniform random(0, 1);
    for (int trial = 0; trial < 100; trial++) {
        // Select a volume in which to generate points.
        
        Vec3 size(10*random.getValue(), 10*random.getValue(), 10*random.getValue());
        Rotation rotation;
        rotation.setRotationToBodyFixedXYZ(Vec3(random.getValue(), random.getValue(), random.getValue()));
        Transform transform(rotation, Vec3(10*random.getValue(), 10*random.getValue(), 10*random.getValue()));
        
        // Create a set of points inside it.
        
        int numPoints = 50*random.getValue()+1;
        Vector_<Vec3> points(numPoints);
        for (int i = 0; i < numPoints; i++)
            points[i] = transform*Vec3(size[0]*random.getValue(), size[1]*random.getValue(), size[2]*random.getValue());
        
        // Create a bounding box from them.
        
        OrientedBoundingBox box(points);
        
        // Verify that it contains all the points.
        
        for (int i = 0; i < numPoints; i++) {
            ASSERT(box.containsPoint(points[i]));
        }
        
        // Verify that it gives a reasonably tight fit to them.
        
        Real expectedVolume = size[0]*size[1]*size[2];
        Real volume = box.getSize()[0]*box.getSize()[1]*box.getSize()[2];
        ASSERT(volume < 2*expectedVolume);
    }
}

void testFindNearestPoint() {
    Vec3 size(1, 1.5, 3);
    Transform trans(Rotation(0.3, XAxis), Vec3(1, 2, 0.5));
    OrientedBoundingBox box(trans, size);
    
    // First test some points inside the box.
    
    Random::Uniform random(0, 1);
    for (int i = 0; i < 100; i++) {
        Vec3 p(random.getValue()*size[0], random.getValue()*size[1], random.getValue()*size[2]);
        p = trans*p;
        assertEqual(p, box.findNearestPoint(p));
    }
    
    // Try some points outside the box.
    
    assertEqual(size, ~trans*box.findNearestPoint(trans*(size+Vec3(1, 2, 3))));
    assertEqual(Vec3(1, 1.5, 0.25), ~trans*box.findNearestPoint(trans*Vec3(2, 3, 0.25)));
    assertEqual(Vec3(0, 0, 0), ~trans*box.findNearestPoint(trans*Vec3(-1, -1, -2)));
    assertEqual(Vec3(0, 0, 0.5), ~trans*box.findNearestPoint(trans*Vec3(-1, -1, 0.5)));
}

int main() {
    try {
        testContainsPoint();
        testGetCorners();
        testIntersectsBox();
        testIntersectsRay();
        testCreateFromPoints();
        testFindNearestPoint();
    }
    catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
