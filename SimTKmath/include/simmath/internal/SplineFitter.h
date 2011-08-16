#ifndef SimTK_SIMMATH_SPLINE_FITTER_H_
#define SimTK_SIMMATH_SPLINE_FITTER_H_

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

#include "SimTKcommon.h"
#include "simmath/internal/common.h"
#include "simmath/internal/Spline.h"
#include "simmath/internal/GCVSPLUtil.h"

namespace SimTK {

/**
 * Given a set of data points, this class creates a Spline_ which interpolates or approximates them.  The
 * data points are assumed to represent a smooth curve plus uncorrelated additive noise.  It attempts to
 * separate these from each other and return a Spline_ which represents the original curve without noise.
 * 
 * The fitting is done based on a <i>smoothing parameter</i>.  When the parameter is 0, the spline will exactly
 * interpolate the data points.  Larger values of the smoothing parameter produce smoother curves that may vary
 * more from the original data.  Since you generally do not know in advance what value for the smoothing
 * parameter is "best", several different methods are provided for selecting it automatically.
 * 
 * If you have no prior information about the structure of the input data, call fitFromGCV():
 * 
 * <pre>
 * SplineFitter<Vec3> fitter = SplineFitter::fitFromGCV(degree, x, y);
 * Spline_<Vec3> spline = fitter.getSpline();
 * </pre>
 * 
 * This chooses a value of the smoothing parameter to minimize the <i>Generalized Cross Validation</i> function.
 * It also estimates the true mean squared error of the data the the number of degrees of freedom of the residual
 * (that is, the number of degrees of freedom not explained by the spline), which can
 * be queried by calling getMeanSquaredError() and getDegreesOfFreedom().  Alternatively, if you have prior
 * knowledge of the error variance or residual degrees of freedom, you can call fitFromErrorVariance()
 * or fitFromDOF() instead.  Finally, you can explicitly specify the smoothing parameter to use by calling
 * fitForSmoothingParameter().
 * 
 * For more information on the GCVSPL algorithm, see Woltring, H.J. (1986), A FORTRAN package for generalized,
 * cross-validatory spline smoothing and differentiation. Advances in Engineering Software 8(2):104-113.  Also,
 * while this class provides access to the most important features of the algorithm, there are a few advanced
 * options which it does not expose directly.  If you need those options, you can access them using the
 * GCVSPLUtil class.
 */

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
    /**
     * Perform a fit, choosing a value of the smoothing parameter that minimizes the Generalized Cross Validation function.
     * 
     * @param degree the degree of the spline to create.  This must be a positive odd value.
     * @param x      the values of the independent variable for each data point
     * @param y      the values of the dependent variables for each data point
     */
    static SplineFitter fitFromGCV(int degree, const Vector& x, const Vector_<T>& y) {
        Vector_<T> coeff;
        Vector wk;
        int ier;
        GCVSPLUtil::gcvspl(x, y, Vector(x.size(), 1.0), static_cast<T>(1), degree, 2, 0, coeff, wk, ier);
        return SplineFitter<T>(new SplineFitterImpl(degree, Spline_<T>(degree, x, coeff), wk[3], wk[4], wk[2]));
    }
    /**
     * Perform a fit, choosing a value of the smoothing parameter based on the known error variance in the data.
     * 
     * @param degree the degree of the spline to create.  This must be a positive odd value.
     * @param x      the values of the independent variable for each data point
     * @param y      the values of the dependent variables for each data point
     * @param error  the variance of the error in the data
     */
    static SplineFitter fitFromErrorVariance(int degree, const Vector& x, const Vector_<T>& y, Real error) {
        Vector_<T> coeff;
        Vector wk;
        int ier;
        GCVSPLUtil::gcvspl(x, y, Vector(x.size(), 1.0), 1, degree, 3, error, coeff, wk, ier);
        return SplineFitter<T>(new SplineFitterImpl(degree, Spline_<T>(degree, x, coeff), wk[3], wk[4], wk[2]));
    }
    /**
     * Perform a fit, choosing a value of the smoothing parameter based on the expect number of degrees of
     * freedom of the residual.
     * 
     * @param degree the degree of the spline to create.  This must be a positive odd value.
     * @param x      the values of the independent variable for each data point
     * @param y      the values of the dependent variables for each data point
     * @param dof    the expected number of degrees of freedom
     */
    static SplineFitter fitFromDOF(int degree, const Vector& x, const Vector_<T>& y, Real dof) {
        Vector_<T> coeff;
        Vector wk;
        int ier;
        GCVSPLUtil::gcvspl(x, y, Vector(x.size(), 1.0), static_cast<T>(1), degree, 4, dof, coeff, wk, ier);
        return SplineFitter<T>(new SplineFitterImpl(degree, Spline_<T>(degree, x, coeff), wk[3], wk[4], wk[2]));
    }
    /**
     * Perform a fit, using a specified fixed value for the smoothing parameter based.
     * 
     * @param degree the degree of the spline to create.  This must be a positive odd value.
     * @param x      the values of the independent variable for each data point
     * @param y      the values of the dependent variables for each data point
     * @param p      the value of the smoothing parameter
     */
    static SplineFitter fitForSmoothingParameter(int degree, const Vector& x, const Vector_<T>& y, Real p) {
        Vector_<T> coeff;
        Vector wk;
        int ier;
        GCVSPLUtil::gcvspl(x, y, Vector(x.size(), 1.0), static_cast<T>(1), degree, 1, p, coeff, wk, ier);
        return SplineFitter<T>(new SplineFitterImpl(degree, Spline_<T>(degree, x, coeff), wk[3], wk[4], wk[2]));
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


