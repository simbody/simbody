/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
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

#include "SimTKcommon.h"

#include <iostream>

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

using std::cout;
using std::endl;
using namespace SimTK;
using namespace std;

const Real TOL = 1e-10;

template <class T>
void assertEqual(T val1, T val2, double tol = TOL) {
    ASSERT(abs(val1-val2) < tol);
}

template <int N>
void assertEqual(Vec<N> val1, Vec<N> val2, double tol) {
    for (int i = 0; i < N; ++i)
        ASSERT(abs(val1[i]-val2[i]) < tol);
}

template <int N, class T>
void assertEqualMat(Mat<N,N,T> val1, Mat<N,N> val2, double tol=TOL) {
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            ASSERT(abs(val1(i, j)-val2(i, j)) < tol);
}

template <int N>
void testInverse() {
    Random::Gaussian random;
    Mat<N,N> identity(1);
    Mat<N,N> mat;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            mat(i, j) = random.getValue();
    assertEqualMat(mat*mat.invert(), identity);
    mat = ~mat;
    assertEqualMat(mat*mat.invert(), identity);
    assertEqualMat((-mat)*(-mat).invert(), identity);
}

void testDotProducts() {
    Vec3 v1(1, 2, 3);
    Vec3 v2(-1, -2, -3);
    Row3 r1(0.1, 0.2, 0.3);
    Row3 r2(-0.1, -0.2, -0.3);
    assertEqual(dot(v1, v2), -14.0);
    assertEqual(dot(r1, r2), -0.14);
    assertEqual(dot(v1, r2), -1.4);
    assertEqual(dot(r1, v2), -1.4);
    assertEqual(r1*v2, -1.4);
    SpatialVec sv(Vec3(1, 2, 3), Vec3(4, 5, 6));
    SpatialRow sr(Row3(1, 2, 3), Row3(4, 5, 6));
    assertEqual(sr*sv, 91.0);
}

int main() {
    try {
        testInverse<1>();
        testInverse<2>();
        testInverse<3>();
        testInverse<5>();
        testInverse<10>();
        testDotProducts();
    } catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
