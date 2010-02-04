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
 * Portions copyright (c) 2008-2009 Stanford University and the Authors.      *
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
 * may be defined in an arbitrary dimensional space.  That is, a Spline_ is a Function_ that calculates an arbitrary number of output
 * values based on a single input value.  The template argument must be either Real or Vec<N> from some integer N.  The name "Spline"
 * (with no trailing _) may be used as a synonym for Spline_<Real>.
 */

template <class T>
class Spline_ : public Function_<T> {
public:
    /**
     * Create a Spline_ object based on a set of control points.
     * 
     * @param degree the degree of the spline to create.  This must be a positive odd value.
     * @param x      the values of the independent variable for each control point
     * @param y      the values of the dependent variables for each control point
     */
    Spline_(int degree, const Vector& x, const Vector_<T>& y) : impl(new SplineImpl(degree, x, y)) {
    }
    Spline_(const Spline_& copy) : impl(copy.impl) {
        impl->referenceCount++;
    }
    Spline_() : impl(NULL) {
    }
    Spline_& operator=(const Spline_& copy) {
        if (impl) {
            impl->referenceCount--;
            if (impl->referenceCount == 0)
                delete impl;
        }
        impl = copy.impl;
        assert(impl);
        impl->referenceCount++;
        return *this;
    }
    ~Spline_() {
        if (impl) {
            impl->referenceCount--;
            if (impl->referenceCount == 0)
                delete impl;
        }
    }
    /**
     * Calculate the values of the dependent variables at a particular value of the independent variable.
     * 
     * @param x a Vector of length 1 containing the value of the independent variable.
     */
    T calcValue(const Vector& x) const {
        assert(impl);
        assert(x.size() == 1);
        return impl->getValue(x[0]);
    }
    /**
     * Calculate a derivative of the spline function.  See the Function_ class for details.  Because Spline_
     * only allows a single independent variable, all elements of derivComponents should be 0, and its
     * length determines the order of the derivative to calculate.
     */
    T calcDerivative(const Array_<int>& derivComponents, const Vector& x) const {
        assert(impl);
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
        assert(impl);
        return impl->x;
    }
    /**
     * Get the values of the dependent variables at the control points.
     */
    const Vector_<T>& getControlPointValues() const {
        assert(impl);
        return impl->y;
    }
    /**
     * Get the degree of the spline.
     */
    int getSplineDegree() const {
        assert(impl);
        return impl->degree;
    }
private:
    class SplineImpl;
    SplineImpl* impl;
};

typedef Spline_<Real> Spline;

template <class T>
class Spline_<T>::SplineImpl {
public:
    SplineImpl(int degree, const Vector& x, const Vector_<T>& y) : referenceCount(1), degree(degree), x(x), y(y) {
    }
    ~SplineImpl() {
        assert(referenceCount == 0);
    }
    T getValue(Real t) const {
        return GCVSPLUtil::splder(0, degree, t, x, y);
    }
    T getDerivative(int derivOrder, Real t) const {
        return GCVSPLUtil::splder(derivOrder, degree, t, x, y);
    }
    int referenceCount;
    int degree;
    Vector x;
    Vector_<T> y;
};

} // namespace SimTK

#endif // SimTK_SIMMATH_SPLINE_H_


