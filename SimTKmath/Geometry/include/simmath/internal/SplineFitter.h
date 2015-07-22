#ifndef SimTK_SIMMATH_SPLINE_FITTER_H_
#define SimTK_SIMMATH_SPLINE_FITTER_H_

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
#include "simmath/internal/Spline.h"
#include "simmath/internal/GCVSPLUtil.h"

namespace SimTK {

/** Given a set of data points, this class creates a Spline_ which interpolates
or approximates them. The data points are assumed to represent a smooth curve
plus uncorrelated additive noise. It attempts to separate these from each other
and return a Spline_ which represents the original curve without noise.

The fitting is done based on a <i>smoothing parameter</i>. When the parameter
is 0, the spline will exactly interpolate the data points. Larger values of the
smoothing parameter produce smoother curves that may vary more from the original
data. Since you generally do not know in advance what value for the smoothing
parameter is "best", several different methods are provided for selecting it
automatically.

If you have no prior information about the structure of the input data, call
fitFromGCV(): @code
    SplineFitter<Vec3> fitter = SplineFitter::fitFromGCV(degree, x, y);
    Spline_<Vec3>      spline = fitter.getSpline();
@endcode

This chooses a value of the smoothing parameter to minimize the <i>Generalized
Cross Validation</i> function. It also estimates the true mean squared error of
the data the the number of degrees of freedom of the residual (that is, the
number of degrees of freedom not explained by the spline), which can be queried
by calling getMeanSquaredError() and getDegreesOfFreedom(). Alternatively, if
you have prior knowledge of the error variance or residual degrees of freedom,
you can call fitFromErrorVariance() or fitFromDOF() instead.  Finally, you can
explicitly specify the smoothing parameter to use by calling
fitForSmoothingParameter().

For more information on the GCVSPL algorithm, see Woltring, H.J. (1986),
"A FORTRAN package for generalized, cross-validatory spline smoothing and
differentiation." Advances in Engineering Software 8(2):104-113.  Also, while
this class provides access to the most important features of the algorithm,
there are a few advanced options which it does not expose directly. If you need
those options, you can access them using the GCVSPLUtil class. **/
template <class T>
class SplineFitter {
public:
    SplineFitter(const SplineFitter& copy) : impl(copy.impl) {
        impl->referenceCount++;
    }
    SplineFitter operator=(const SplineFitter& copy) {
        impl = copy.impl;
        impl->referenceCount++;
        return *this;
    }
    ~SplineFitter() {
        impl->referenceCount--;
        if (impl->referenceCount == 0)
            delete impl;
    }

    /** Perform a fit, choosing a value of the smoothing parameter that
    minimizes the Generalized Cross Validation function.
    @param[in] degree The degree of the spline to create. This must be a
                        positive, odd integer.
    @param[in] x      The value of the independent variable for each data point.
    @param[in] y      The values of the dependent variables for each data point.
    @return A SplineFitter object containing the desired Spline_. **/
    static SplineFitter fitFromGCV
       (int degree, const Vector& x, const Vector_<T>& y) {
        Vector_<T> coeff;
        Vector wk;
        int ier;
        GCVSPLUtil::gcvspl(x, y, Vector(x.size(), 1.0), T(1),
                           degree, 2, 0, coeff, wk, ier);
        return SplineFitter<T>(new SplineFitterImpl
           (degree, Spline_<T>(degree, x, coeff), wk[3], wk[4], wk[2]));
    }

    /** Perform a fit, choosing a value of the smoothing parameter based on the
    known error variance in the data.
    @param[in] degree The degree of the spline to create. This must be a
                        positive, odd integer.
    @param[in] x      The value of the independent variable for each data point.
    @param[in] y      The values of the dependent variables for each data point.
    @param[in] error  The variance of the error in the data.
    @return A SplineFitter object containing the desired Spline_. **/
    static SplineFitter fitFromErrorVariance
       (int degree, const Vector& x, const Vector_<T>& y, Real error) {
        Vector_<T> coeff;
        Vector wk;
        int ier;
        GCVSPLUtil::gcvspl(x, y, Vector(x.size(), 1.0), T(1),
                           degree, 3, error, coeff, wk, ier);
        return SplineFitter<T>(new SplineFitterImpl
           (degree, Spline_<T>(degree, x, coeff), wk[3], wk[4], wk[2]));
    }

    /** Perform a fit, choosing a value of the smoothing parameter based on the
    expect number of degrees of freedom of the residual.
    @param[in] degree The degree of the spline to create. This must be a
                        positive, odd value.
    @param[in] x      The value of the independent variable for each data point.
    @param[in] y      The values of the dependent variables for each data point.
    @param[in] dof    The expected number of degrees of freedom.
    @return A SplineFitter object containing the desired Spline_. **/
    static SplineFitter fitFromDOF
       (int degree, const Vector& x, const Vector_<T>& y, Real dof) {
        Vector_<T> coeff;
        Vector wk;
        int ier;
        GCVSPLUtil::gcvspl(x, y, Vector(x.size(), 1.0), T(1),
                           degree, 4, dof, coeff, wk, ier);
        return SplineFitter<T>(new SplineFitterImpl
           (degree, Spline_<T>(degree, x, coeff), wk[3], wk[4], wk[2]));
    }

    /** Perform a fit, using a specified fixed value for the smoothing
    parameter.
    @param[in] degree The degree of the spline to create. This must be a
                        positive, odd value.
    @param[in] x      The value of the independent variable for each data point.
    @param[in] y      The values of the dependent variables for each data point.
    @param[in] p      The value of the smoothing parameter.
    @return A SplineFitter object containing the desired Spline_. **/
    static SplineFitter fitForSmoothingParameter
       (int degree, const Vector& x, const Vector_<T>& y, Real p) {
        Vector_<T> coeff;
        Vector wk;
        int ier;
        GCVSPLUtil::gcvspl(x, y, Vector(x.size(), 1.0), T(1),
                           degree, 1, p, coeff, wk, ier);
        return SplineFitter<T>(new SplineFitterImpl
           (degree, Spline_<T>(degree, x, coeff), wk[3], wk[4], wk[2]));
    }
    /**
     * Get the Spline_ that was generated by the fitting.
     */
    const Spline_<T>& getSpline() {
        return impl->spline;
    }
    /**
     * Get the smoothing parameter that was used for the fitting.
     */
    Real getSmoothingParameter() {
        return impl->p;
    }
    /**
     * Get the estimate of the true mean squared error in the data that was determined by the fitting.
     */
    Real getMeanSquaredError() {
        return impl->error;
    }
    /**
     * Get the estimate of the number of degrees of freedom of the residual that was determined by the fitting.
     */
    Real getDegreesOfFreedom() {
        return impl->dof;
    }
private:
    class SplineFitterImpl;
    SplineFitter(SplineFitterImpl *impl) : impl(impl) {
    }
    SplineFitterImpl* impl;
};

template <class T>
class SplineFitter<T>::SplineFitterImpl {
public:
    SplineFitterImpl(int degree, const Spline_<T>& spline, Real p, Real error, Real dof) : referenceCount(1), degree(degree), spline(spline), p(p), error(error), dof(dof) {
    }
    ~SplineFitterImpl() {
        assert(referenceCount == 0);
    }
    int referenceCount;
    int degree;
    Spline_<T> spline;
    Real p, error, dof;
};

} // namespace SimTK

#endif // SimTK_SIMMATH_SPLINE_FITTER_H_


