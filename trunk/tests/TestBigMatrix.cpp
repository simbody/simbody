/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007-8 Stanford University and the Authors.           *
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

template <class T, int N>
void testVector(const T& value, const Vec<N>& expected) {
    ASSERT(value.size() == N);
    for (int i = 0; i < N; ++i) {
        if (isNaN(expected[i])) {
            ASSERT(isNaN(value[i]));
        }
        else {
            ASSERT(value[i] == expected[i]);
        }
    }
}

template <class T, int M, int N>
void testMatrix(const T& value, const Mat<M, N>& expected) {
    ASSERT(value.nrow() == M);
    ASSERT(value.ncol() == N);
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < N; ++j) {
            if (isNaN(expected(i, j))) {
                ASSERT(isNaN(value(i, j)));
            }
            else {
                ASSERT(value(i, j) == expected(i, j));
            }
        }
}


void testMatDivision() {
    // scalar / matrix seems to be causing trouble
    SpatialMat m1;
    Mat33 m2;
    // Each of the following causes compile errors
    1.0/m2;
    1.0/m1;
}

int main() {
    try {
        // Currently, this only tests a small number of operations that were recently added.
        // It should be expanded into a more comprehensive test of the big matrix classes.

        testMatDivision();
        
        Matrix m(Mat22(1, 2, 3, 4));
        testMatrix<Matrix,2,2>(m, Mat22(1, 2, 3, 4));
        m += 3;
        testMatrix<Matrix,2,2>(m, Mat22(4, 2, 3, 7));
        m -= 3;
        testMatrix<Matrix,2,2>(m, Mat22(1, 2, 3, 4));
        testMatrix<Matrix,2,2>(m-1, Mat22(0, 2, 3, 3));
        testMatrix<Matrix,2,2>(m+1, Mat22(2, 2, 3, 5));
        testMatrix<Matrix,2,2>(1-m, Mat22(0, -2, -3, -3));
        testMatrix<Matrix,2,2>(1+m, Mat22(2, 2, 3, 5));
        Vector v(Vec3(1, 2, 3));
        testVector(v, Vec3(1, 2, 3));
        v += 2;
        testVector(v, Vec3(3, 4, 5));
        v -= 2;
        testVector(v, Vec3(1, 2, 3));
        testVector(v-1, Vec3(0, 1, 2));
        testVector(v+1, Vec3(2, 3, 4));
        testVector(1-v, Vec3(0, -1, -2));
        testVector(1+v, Vec3(2, 3, 4));
        RowVector r(Row3(1, 2, 3));
        testVector(r, Vec3(1, 2, 3));
        r += 2;
        testVector(r, Vec3(3, 4, 5));
        r -= 2;
        testVector(r, Vec3(1, 2, 3));
        testVector(r-1, Vec3(0, 1, 2));
        testVector(r+1, Vec3(2, 3, 4));
        testVector(1-r, Vec3(0, -1, -2));
        testVector(1+r, Vec3(2, 3, 4));

        Matrix mm( Mat23( 1, 2, 3,
                          7, 8, 9 ) );
        testMatrix<Matrix,2,3>(mm, Mat23(1,2,3,7,8,9));

            // Test copying a column or row of a Matrix into
            // a Vector or RowVector.

        // Test assignment constructor
        Vector vv = mm(1); testVector(vv, Vec2(2,8));
        // Test copy assignment
        vv = mm(0); testVector(vv, Vec2(1,7));
        // Test assignment constructor
        RowVector rr = mm[1]; testVector(rr, Vec3(7,8,9));
        // Test copy assignment
        rr = mm[0]; testVector(rr, Vec3(1,2,3));

            // Test copying a row into a Vector and column into RowVector.

        // Test assignment (copy) constructor
        RowVector rrr = ~mm(1); 
        testVector(rrr, Vec2(2,8));
        // Test copy assignment
        rrr = ~mm(0); testVector(rrr, Vec2(1,7));

        // Test assignment (copy) constructor
        Vector vvv = ~mm[1]; 
        testVector(vvv, Vec3(7,8,9));
        // Test copy assignment
        vvv = ~mm[0]; testVector(vvv, Vec3(1,2,3));

    } catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
