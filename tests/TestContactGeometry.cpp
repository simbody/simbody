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
#include <vector>
#include <exception>

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

int main() {
    try {
        testHalfSpace();
        testSphere();
    }
    catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
