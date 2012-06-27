/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
 * Contributors: Michael Sherman                                              *
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

// Note: this file was moved from Simmath to SimTKcommon 20100601; see the
// Simmath repository for earlier history.

#include "SimTKcommon.h"
#include "SimTKcommon/Testing.h"


// We'll use some std::vectors to check interoperability between Array_<T>
// and std::vector<T>.
#include <vector>

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
    const Vec3 df = f.calcDerivative(derivComponents, x);
    assertEqual(Vec3(0), df);
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
    std::vector<int> derivComponents(1);
    derivComponents[0] = 1;
    assertEqual(Vec3(4, 3, 2), f.calcDerivative(derivComponents, Vector(Vec2(1, 0))));
    std::vector<int> derivComponents2(2);
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
    std::vector<int> derivComponents(1);
    assertEqual(Vec3(4, 3, 2), f.calcDerivative(derivComponents, Vector(Vec1(0))));
    assertEqual(Vec3(6, 7, 8), f.calcDerivative(derivComponents, Vector(Vec1(1))));
    assertEqual(Vec3(8, 11, 14), f.calcDerivative(derivComponents, Vector(Vec1(2))));
    std::vector<int> derivComponents2(2);
    assertEqual(Vec3(2, 4, 6), f.calcDerivative(derivComponents2, Vector(Vec1(0))));
    assertEqual(Vec3(2, 4, 6), f.calcDerivative(derivComponents2, Vector(Vec1(1))));
    std::vector<int> derivComponents3(3);
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

void testSinusoid() { 
    Real a=11.23, w=1.1, p=Pi/4;
    Vector t1(1,.23), t2(1,-3.2), t3(1,14.1);
    Function::Sinusoid s1(a,w,p);
    SimTK_TEST_EQ(s1.calcValue(Vector(1,0.)), a*std::sin(p));
    SimTK_TEST_EQ(s1.calcValue(t1), a*std::sin(w*t1[0]+p));
    SimTK_TEST_EQ(s1.calcValue(t2), a*std::sin(w*t2[0]+p));
    SimTK_TEST_EQ(s1.calcValue(t3), a*std::sin(w*t3[0]+p));

    // Do enough of these to make sure we reach the general forumula.
    Array_<int> deriv; // 0th derivative is function
    SimTK_TEST_EQ(s1.calcDerivative(deriv, t1), s1.calcValue(t1)); // 0th
    deriv.push_back(0); // 1st deriv
    SimTK_TEST_EQ(s1.calcDerivative(deriv, t2), a*w*std::cos(w*t2[0]+p));
    deriv.push_back(0); // 2nd deriv
    SimTK_TEST_EQ(s1.calcDerivative(deriv, t3), -a*w*w*std::sin(w*t3[0]+p));
    deriv.push_back(0); // 3rd deriv
    SimTK_TEST_EQ(s1.calcDerivative(deriv, t1), -a*w*w*w*std::cos(w*t1[0]+p));
    deriv.push_back(0); // 4th deriv
    SimTK_TEST_EQ(s1.calcDerivative(deriv, t1), a*w*w*w*w*std::sin(w*t1[0]+p));
    deriv.push_back(0); // 5th deriv
    SimTK_TEST_EQ(s1.calcDerivative(deriv, t2), a*w*w*w*w*w*std::cos(w*t2[0]+p));
    deriv.push_back(0); // 6th deriv
    SimTK_TEST_EQ(s1.calcDerivative(deriv, t2), -a*w*w*w*w*w*w*std::sin(w*t2[0]+p));
    deriv.push_back(0); // 7th deriv
    SimTK_TEST_EQ(s1.calcDerivative(deriv, t2), -a*w*w*w*w*w*w*w*std::cos(w*t2[0]+p));
}

void testStep() {
    Function::Step s1(-1,1,0,1); // y in [-1,1] as x in [0,1]
    SimTK_TEST(s1.calcValue(Vector(1,Zero)) == -1);    // x0 -> y0
    SimTK_TEST(s1.calcValue(Vector(1,One)) ==  1);    // x1 -> y1
    SimTK_TEST(s1.calcValue(Vector(1,OneHalf)) == 0); // 1/2 -> (y1+y0)/2
    SimTK_TEST(s1.calcValue(Vector(1,Real(-29))) == -1);
    SimTK_TEST(s1.calcValue(Vector(1,Real(234.3))) == 1);

    // First & second derivs should be zero at either end.
    Array_<int> derivOrder1(1);
    SimTK_TEST(s1.calcDerivative(derivOrder1, Vector(1,Zero)) == 0);
    SimTK_TEST(s1.calcDerivative(derivOrder1, Vector(1,One)) == 0);
    SimTK_TEST(s1.calcDerivative(derivOrder1, Vector(1,Real(-29))) == 0);
    SimTK_TEST(s1.calcDerivative(derivOrder1, Vector(1,Real(234.3))) == 0);
    Array_<int> derivOrder2(2);
    SimTK_TEST(s1.calcDerivative(derivOrder2, Vector(1,Zero)) == 0);
    SimTK_TEST(s1.calcDerivative(derivOrder2, Vector(1,One)) == 0);
    SimTK_TEST(s1.calcDerivative(derivOrder2, Vector(1,Real(-29))) == 0);
    SimTK_TEST(s1.calcDerivative(derivOrder2, Vector(1,Real(234.3))) == 0);
    Array_<int> derivOrder3(3); // don't know much about 3rd derivative
    SimTK_TEST(s1.calcDerivative(derivOrder3, Vector(1,Real(-29))) == 0);
    SimTK_TEST(s1.calcDerivative(derivOrder3, Vector(1,Real(234.3))) == 0);

    // Try a more general step with x0,x1 reversed also.
    // Here y goes from -221.3 to 47.9 as x goes from 1000 down to -333.
    const Real y0=-221.3, y1=47.9, x0=1000, x1=-333;
    Function::Step s2(y0,y1,x0,x1);
    SimTK_TEST(s2.calcValue(Vector(1,x0)) == y0);    // x0 -> y0
    SimTK_TEST(s2.calcValue(Vector(1,x1)) == y1);    // x1 -> y1
    SimTK_TEST_EQ(s2.calcValue(Vector(1,(x1+x0)/2)), (y1+y0)/2); // (x1+x0)/2 -> (y1+y0)/2
    SimTK_TEST(s2.calcValue(Vector(1,x0+100)) == y0); // note sign
    SimTK_TEST(s2.calcValue(Vector(1,x1-100)) == y1);

    // Calculate 3rd deriv by differencing 2nd
    const Real x = -22.701, dx = 1e-6;
    const Real d2m=s2.calcDerivative(derivOrder2, Vector(1, x-dx));
    const Real d2p=s2.calcDerivative(derivOrder2, Vector(1, x+dx));
    const Real d3approx = (d2p-d2m)/(2*dx); // approx 10 digits
    SimTK_TEST_EQ_TOL(s2.calcDerivative(derivOrder3, Vector(1,x)), 
                      d3approx, 1e-8);

    // Try interpolating a Vec3
    Function_<Vec3>::Step sv(Vec3(1,2,3), Vec3(4,5,6), 0, 1);
    SimTK_TEST(sv.calcValue(Vector(1,OneHalf)) == Vec3(2.5,3.5,4.5));
    SimTK_TEST(sv.calcDerivative(derivOrder2, Vector(1, -29.3)) == Vec3(0));
}

int main () {
    SimTK_START_TEST("TestFunction");

        SimTK_SUBTEST(testConstant);
        SimTK_SUBTEST(testLinear);
        SimTK_SUBTEST(testPolynomial);
        SimTK_SUBTEST(testSinusoid);
        SimTK_SUBTEST(testRealFunction);
        SimTK_SUBTEST(testStep);

    SimTK_END_TEST();
}
