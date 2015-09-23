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
 * Portions copyright (c) 2008-13 Stanford University and the Authors.        *
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
    /** Create a Spline_ object based on a set of control points.\ See
    SplineFitter for a nicer way to create these objects.

    @param[in]  degree
        The degree of the spline to create. This must be a positive odd value.
    @param[in]  x
        Values of the independent variable for each Bezier control point.
    @param[in]  y
        Values of the dependent variables for each Bezier control point.
    **/
    Spline_(int degree, const Vector& x, const Vector_<T>& y)
    :   impl(new SplineImpl(degree, x, y)) {}

    /** Default constructor creates an empty %Spline_ handle; not very
    useful. **/
    Spline_() : impl(nullptr) {}

    /** Copy constructor is shallow and reference-counted; that is, the new
    %Spline_ refers to the same object as does \a source. **/
    Spline_(const Spline_& source) : impl(source.impl)
    {   if (impl) impl->referenceCount++; }

    /** Copy assignment is shallow and reference-counted; that is, after the
    assignment this %Spline_ refers to the same object as does \a source. **/
    Spline_& operator=(const Spline_& source) {
        if (impl) {
            impl->referenceCount--;
            if (impl->referenceCount == 0)
                delete impl;
        }
        impl = source.impl;
        if (impl) impl->referenceCount++;
        return *this;
    }

    /** Destructor decrements the reference count and frees the heap space if
    this is the last reference. **/
    ~Spline_() {
        if (impl) {
            impl->referenceCount--;
            if (impl->referenceCount == 0)
                delete impl;
        }
    }

    /** Calculate the values of the dependent variables at a particular value
    of the independent variable.
    @param[in]  x    The value of the independent variable.
    @returns         The corresponding values of the dependent variables. **/
    T calcValue(Real x) const {
        assert(impl);
        return impl->getValue(x);
    }

    /** Calculate a derivative of the spline function with respect to its
    independent variable, at the given value.
    @param[in]  order
        Which derivative? order==1 returns first derivative, order==2 is second,
        etc. Must be >= 1.
    @param[in]  x
        The value of the independent variable about which the derivative is
        taken.
    @returns
        The \a order'th derivative of the dependent variables at \a x. **/
    T calcDerivative(int order, Real x) const {
        assert(impl);
        assert(order > 0);
        return impl->getDerivative(order, x);
    }

    /** Get the locations (that is, the values of the independent variable) for
    each of the Bezier control points. **/
    const Vector& getControlPointLocations() const {
        assert(impl);
        return impl->x;
    }
    /** Get the values of the dependent variables at each of the Bezier control
    points. **/
    const Vector_<T>& getControlPointValues() const {
        assert(impl);
        return impl->y;
    }

    /** Get the degree of the spline. **/
    int getSplineDegree() const {
        assert(impl);
        return impl->degree;
    }

    /** Alternate signature provided to implement the generic Function_
    interface expects a one-element Vector for the independent variable. **/
    T calcValue(const Vector& x) const override {
        assert(x.size() == 1);
        return calcValue(x[0]);
    }
    /** Alternate signature provided to implement the generic Function_
    interface expects an awkward \a derivComponents argument, and
    takes a one-element Vector for the independent variable. Because
    Spline_ only allows a single independent variable, all elements of
    derivComponents should be 0; only its length determines the order of the
    derivative to calculate. **/
    T calcDerivative(const Array_<int>& derivComponents, const Vector& x) const
        override
    {   assert(x.size() == 1);
        return calcDerivative((int)derivComponents.size(), x[0]); }
    /** For the Function_ style interface, this provides compatibility
    with std::vector. No copying or heap allocation is required. **/
    T calcDerivative(const std::vector<int>& derivComponents,
                     const Vector& x) const
    {   assert(x.size() == 1);
        return calcDerivative((int)derivComponents.size(), x[0]); }

    /** Required by the Function_ interface. **/
    int getArgumentSize() const override {return 1;}
    /** Required by the Function_ interface. **/
    int getMaxDerivativeOrder() const override
    {   return std::numeric_limits<int>::max(); }
    /** Required by the Function_ interface. **/
    Spline_* clone() const override {return new Spline_(*this);}

private:
    class SplineImpl;
    SplineImpl* impl;
};

/** Provide a convenient name for a scalar-valued Spline_. **/
typedef Spline_<Real> Spline;

/** This is the implementation class that supports the Spline_ interface. **/
template <class T>
class Spline_<T>::SplineImpl {
public:
    SplineImpl(int degree, const Vector& x, const Vector_<T>& y)
    :   referenceCount(1), degree(degree), x(x), y(y) {}
    ~SplineImpl() {
        assert(referenceCount == 0);
    }
    T getValue(Real t) const {
        return GCVSPLUtil::splder(0, degree, t, x, y);
    }
    T getDerivative(int derivOrder, Real t) const {
        return GCVSPLUtil::splder(derivOrder, degree, t, x, y);
    }
    int         referenceCount;
    int         degree;
    Vector      x;
    Vector_<T>  y;
};

} // namespace SimTK

#endif // SimTK_SIMMATH_SPLINE_H_


