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
#include <vector>

using namespace SimTK;
using namespace std;

const Real TOL = 1e-10;

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

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

void testSpline() {
    Vector_<Vec3> coeff(5);
    coeff[0] = Vec3(0, 1, 2);
    coeff[1] = Vec3(1, 4, 1);
    coeff[2] = Vec3(2, 2, 20);
    coeff[3] = Vec3(1, -1, 2);
    coeff[4] = Vec3(0, 0, 1);
    Vector x(Vec5(0, 1, 2, 5, 10));
    
    // Create a linear spline, and verify that it interpolates linearly between the control points.
    
    Spline<3> spline(1, x, coeff);
    for (int i = 0; i < x.size(); ++i)
        assertEqual(coeff[i], spline.calcValue(Vector(1, x[i])));
    vector<int> deriv;
    deriv.push_back(0);
    for (int i = 0; i < x.size()-1; ++i) {
        for (int j = 0; j < 10; ++j) {
            Real fract = (i+1.0)/12.0;
            Real t = x[i]+fract*(x[i+1]-x[i]);
            assertEqual(spline.calcValue(Vector(1, t)), coeff[i]+fract*(coeff[i+1]-coeff[i]));
            assertEqual(spline.calcDerivative(deriv, Vector(1, t)), (coeff[i+1]-coeff[i])/(x[i+1]-x[i]));
        }
    }
    
    // Create a cubic spline and verify the derivative calculations.
    
    spline = Spline<3>(3, x, coeff);
    Real delta = 1e-10;
    for (int i = 0; i < x.size()-1; ++i) {
        for (int j = 0; j < 10; ++j) {
            Real fract = (i+1.0)/12.0;
            Real t = x[i]+fract*(x[i+1]-x[i]);
            Vec3 value1 = spline.calcValue(Vector(1, t-delta));
            Vec3 value2 = spline.calcValue(Vector(1, t+delta));
            assertEqual(spline.calcDerivative(deriv, Vector(1, t)), (value2-value1)/(2*delta), 1e-4);
        }
    }
}

void testSplineFitter() {
    Real stddev = 0.5;
    int n = 100;
    Random::Gaussian random(0.0, stddev);
    Vector x(n);
    Vector_<Vec3> truey(n);
    Vector_<Vec3> y(n);
    for (int i = 0; i < x.size(); ++i) {
        x[i] = i*0.1;
        truey[i] = Vec3(sin(x[i]), 3.0*sin(2*x[i]), cos(x[i]));
        y[i] = truey[i] + Vec3(random.getValue(), random.getValue(), random.getValue());
    }
    SplineFitter<3> fitter = SplineFitter<3>::fitFromGCV(1, x, y);
    Spline<3> spline1 = fitter.getSpline();
    
    // The fitting should have reduced the error.
    
    Vec3 originalError = mean(abs(y-truey));
    Vec3 fittedError = mean(abs(spline1.getControlPointValues()-truey));
    ASSERT(fittedError[0] < originalError[0]);
    ASSERT(fittedError[1] < originalError[1]);
    ASSERT(fittedError[2] < originalError[2]);
    
    // If we perform the fitting again, explicitly specifying the same value for the smoothing parameter,
    // it should produce identical results.
    
    assertEqual(SplineFitter<3>::fitForSmoothingParameter(1, x, y, fitter.getSmoothingParameter()).getSpline().getControlPointValues(), spline1.getControlPointValues());
    
    // If we specify a smoothing parameter of 0, it should exactly reproduce the original data.

    assertEqual(SplineFitter<3>::fitForSmoothingParameter(1, x, y, 0.0).getSpline().getControlPointValues(), y);
}

int main () {
    try {
        testSpline();
        testSplineFitter();
        cout << "Done" << endl;
        return 0;
    }
    catch (std::exception& e) {
        std::printf("FAILED: %s\n", e.what());
        return 1;
    }
}
