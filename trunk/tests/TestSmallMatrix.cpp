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

template <int N, class T>
void assertEqual(Mat<N,N,T> val1, Mat<N,N> val2, double tol=TOL) {
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
    assertEqual(mat*mat.invert(), identity);
    mat = ~mat;
    assertEqual(mat*mat.invert(), identity);
    assertEqual((-mat)*(-mat).invert(), identity);
}

int main() {
    try {
        testInverse<1>();
        testInverse<2>();
        testInverse<3>();
        testInverse<5>();
        testInverse<10>();
    } catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
