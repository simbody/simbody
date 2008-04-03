#ifndef SimTK_SIMMATH_SPLINE_H_
#define SimTK_SIMMATH_SPLINE_H_

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

#include "simmath/internal/Function.h"
#include "simmath/internal/GCVSPLUtil.h"

#include <limits>

namespace SimTK {

/**
 * This class implements a non-uniform B-spline curve.  It requires the spline degree to be odd (linear, cubic, quintic, etc.), but
 * supports arbitrarily high degrees.  Only spline curves are supported, not surfaces or higher dimensional objects, but the curve
 * may be defined in an arbitrary dimensional space.  That is, a Spline is a Function that calculates an arbitrary number of output
 * values based on a single input value.
 */

template <int N>
class Spline : public Function<N> {
public:
    /**
     * Create a Spline object based on a set of control points.
     * 
     * @param degree the degree of the spline to create.  This must be a positive odd value.
     * @param x      the values of the independent variable for each control point
     * @param y      the values of the dependent variables for each control point
     */
    Spline(int degree, const Vector& x, const Vector_<Vec<N> >& y) : impl(new SplineImpl(degree, x, y)) {
    }
    Spline(const Spline& copy) : impl(copy.impl) {
        impl->referenceCount++;
    }
    Spline operator=(const Spline& copy) {
        return Spline(copy);
    }
    ~Spline() {
        impl->referenceCount--;
        if (impl->referenceCount == 0)
            delete impl;
    }
    /**
     * Calculate the values of the dependent variables at a particular value of the independent variable.
     * 
     * @param x a Vector of length 1 containing the value of the independent variable.
     */
    Vec<N> calcValue(const Vector& x) const {
        assert(x.size() == 1);
        return impl->getValue(x[0]);
    }
    /**
     * Calculate a derivative of the spline function.  See the Function class for details.  Because Spline
     * only allows a single independent variable, all elements of derivComponents should be 0, and its
     * length determines the order of the derivative to calculate.
     */
    Vec<N> calcDerivative(const std::vector<int>& derivComponents, const Vector& x) const {
        assert(x.size() == 1);
        assert(derivComponents.size() > 0);
        return impl->getDerivative(derivComponents.size(), x[0]);
    }
    int getArgumentSize() const {
        return 1;
    }
    int getMaxDerivativeOrder() const {
        return std::numeric_limits<int>::max();
    }
    /**
     * Get the locations (that is, the values of the independent variable) of the control points.
     */
    const Vector& getControlPointLocations() const {
        return impl->x;
    }
    /**
     * Get the values of the dependent variables at the control points.
     */
    const Vector_<Vec<N> >& getControlPointValues() const {
        return impl->y;
    }
    /**
     * Get the degree of the spline.
     */
    int getSplineDegree() const {
        return impl->degree;
    }
private:
    class SplineImpl;
    Spline();
    SplineImpl* impl;
};

template <int N>
class Spline<N>::SplineImpl {
public:
    SplineImpl(int degree, const Vector& x, const Vector_<Vec<N> >& y) : degree(degree), x(x), y(y), referenceCount(1) {
    }
    Vec<N> getValue(Real t) const {
        return GCVSPLUtil::splder(0, degree, t, x, y);
    }
    Vec<N> getDerivative(int derivOrder, Real t) const {
        return GCVSPLUtil::splder(derivOrder, degree, t, x, y);
    }
    int referenceCount;
    int degree;
    Vector x;
    Vector_<Vec<N> > y;
};

} // namespace SimTK

#endif // SimTK_SIMMATH_SPLINE_H_


