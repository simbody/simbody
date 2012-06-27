#ifndef SimTK_SIMMATH_SPLINE_H_
#define SimTK_SIMMATH_SPLINE_H_

/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
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

#include "SimTKcommon.h"
#include "simmath/internal/common.h"
#include "simmath/internal/GCVSPLUtil.h"

#include <limits>

namespace SimTK {

/** This class implements a non-uniform Bezier curve. It requires the spline 
degree to be odd (linear, cubic, quintic, etc.), but supports arbitrarily high
degrees. Only spline curves are supported, not surfaces or higher dimensional 
objects, but the curve may be defined in an arbitrary dimensional space.
That is, a Spline_ is a Function_ that calculates an arbitrary number of 
output values based on a single input value. The template argument must be 
either Real or Vec<N> from some integer N.  The name "Spline" (with no trailing
"_") may be used as a synonym for Spline_<Real>.

Most users should generate Spline_ objects using SplineFitter which can be
used to generate spline curves through data points rather than requiring
Bezier control points.

@see SplineFitter for best-fitting a spline through sampled data.
@see BicubicSurface for fitting a smooth surface to 2D sample data.
**/
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
    Spline_(int degree, const Vector& x, const Vector_<T>& y) 
    :   impl(new SplineImpl(degree, x, y)) {}

    Spline_(const Spline_& copy) : impl(copy.impl) {
        if (impl) impl->referenceCount++;
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
        if (impl) impl->referenceCount++;
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

    /** This provides compatibility with std::vector without requiring any copying. **/
    T calcDerivative(const std::vector<int>& derivComponents, const Vector& x) const 
    {   return calcDerivative(ArrayViewConst_<int>(derivComponents),x); }
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


