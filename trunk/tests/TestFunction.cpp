/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmath(tm)                         *
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

#include "SimTKmath.h"

using namespace SimTK;
using namespace std;

const Real TOL = 1e-10;

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

void assertEqual(Real val1, Real val2, double tol=TOL) {
    ASSERT(abs(val1-val2) < tol);
}

template <int N>
void assertEqual(Vec<N> val1, Vec<N> val2, double tol=TOL) {
    for (int i = 0; i < N; ++i)
        ASSERT(abs(val1[i]-val2[i]) < tol);
}

template <class T>
void assertEqual(Vector_<T> val1, Vector_<T> val2, double tol=TOL) {
    ASSERT(val1.size() == val2.size());
    for (int i = 0; i < val1.size(); ++i)
        assertEqual(val1[i], val2[i], tol);
}

void testConstant() {
    Function_<Vec3>::Constant f(Vec3(1, 2, 3), 2);
    ASSERT(f.getArgumentSize() == 2);
    Vector x(2);
    assertEqual(Vec3(1, 2, 3), f.calcValue(x));
    Array_<int> derivComponents(1);
    assertEqual(Vec3(0), f.calcDerivative(derivComponents, x));
}

void testLinear() {
    Vector_<Vec3> coeff(3);
    coeff[0] = Vec3(1, 2, 3);
    coeff[1] = Vec3(4, 3, 2);
    coeff[2] = Vec3(-1, -2, -3);
    Function_<Vec3>::Linear f(coeff);
    ASSERT(f.getArgumentSize() == 2);
    assertEqual(Vec3(-1, -2, -3), f.calcValue(Vector(Vec2(0, 0))));
    assertEqual(Vec3(0, 0, 0), f.calcValue(Vector(Vec2(1, 0))));
    assertEqual(Vec3(-2.5, -2.5, -2.5), f.calcValue(Vector(Vec2(0.5, -0.5))));
    Array_<int> derivComponents(1);
    derivComponents[0] = 1;
    assertEqual(Vec3(4, 3, 2), f.calcDerivative(derivComponents, Vector(Vec2(1, 0))));
    Array_<int> derivComponents2(2);
    assertEqual(Vec3(0, 0, 0), f.calcDerivative(derivComponents2, Vector(Vec2(1, 0))));
}

void testPolynomial() {
    Vector_<Vec3> coeff(3);
    coeff[0] = Vec3(1, 2, 3);
    coeff[1] = Vec3(4, 3, 2);
    coeff[2] = Vec3(-1, -2, -3);
    Function_<Vec3>::Polynomial f(coeff);
    ASSERT(f.getArgumentSize() == 1);
    assertEqual(Vec3(-1, -2, -3), f.calcValue(Vector(Vec1(0))));
    assertEqual(Vec3(4, 3, 2), f.calcValue(Vector(Vec1(1))));
    assertEqual(Vec3(11, 12, 13), f.calcValue(Vector(Vec1(2))));
    Array_<int> derivComponents(1);
    assertEqual(Vec3(4, 3, 2), f.calcDerivative(derivComponents, Vector(Vec1(0))));
    assertEqual(Vec3(6, 7, 8), f.calcDerivative(derivComponents, Vector(Vec1(1))));
    assertEqual(Vec3(8, 11, 14), f.calcDerivative(derivComponents, Vector(Vec1(2))));
    Array_<int> derivComponents2(2);
    assertEqual(Vec3(2, 4, 6), f.calcDerivative(derivComponents2, Vector(Vec1(0))));
    assertEqual(Vec3(2, 4, 6), f.calcDerivative(derivComponents2, Vector(Vec1(1))));
    Array_<int> derivComponents3(3);
    assertEqual(Vec3(0, 0, 0), f.calcDerivative(derivComponents3, Vector(Vec1(1))));
}

void testRealFunction() {
    Vector coeff(3);
    coeff[0] = 1.0;
    coeff[1] = 4.0;
    coeff[2] = -1.0;
    Function::Linear f(coeff);
    ASSERT(f.getArgumentSize() == 2);
    assertEqual(-1, f.calcValue(Vector(Vec2(0, 0))));
    assertEqual(0, f.calcValue(Vector(Vec2(1, 0))));
    assertEqual(-2.5, f.calcValue(Vector(Vec2(0.5, -0.5))));
    Array_<int> derivComponents(1);
    derivComponents[0] = 1;
    assertEqual(4, f.calcDerivative(derivComponents, Vector(Vec2(1, 0))));
    Array_<int> derivComponents2(2);
    assertEqual(0, f.calcDerivative(derivComponents2, Vector(Vec2(1, 0))));
}

int main () {
    try {
        testConstant();
        testLinear();
        testPolynomial();
        testRealFunction();
        cout << "Done" << endl;
        return 0;
    }
    catch (std::exception& e) {
        std::printf("FAILED: %s\n", e.what());
        return 1;
    }
}
