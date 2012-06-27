#ifndef SimTK_SimTKCOMMON_POLYNOMIALROOTFINDER_H_
#define SimTK_SimTKCOMMON_POLYNOMIALROOTFINDER_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-12 Stanford University and the Authors.        *
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

#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"

namespace SimTK {

/**
 * This class provides static methods for finding the roots of polynomials.  There are
 * specialized methods for quadratic and cubic polynomials, as well as general methods
 * for polynomials of arbitrary degree.  In each case, there are methods for polynomials
 * with both real and complex coefficients.
 * 
 * There are two different algorithms used by this class.  The specialized methods for
 * quadratic polynomials calculate the roots by explicit evaluation of the quadratic
 * formula.  They use the evaluation method described in section 5.6 of "Numerical
 * Recipes in C++, Second Edition", by Press, Teukolsky, Vetterling, and Flannery.
 * In addition, the method for quadratic polynomials with real coefficients performs
 * an extra check to detect when the discriminant is zero to within machine precision.
 * This helps to prevent round-off error from producing a tiny imaginary part in a
 * multiple root.
 * 
 * The methods for cubic and arbitrary degree polynomials use the Jenkins-Traub method,
 * as implemented in the classic RPOLY and CPOLY functions:
 * 
 * Jenkins, M. A. and Traub, J. F. (1972), Algorithm 419: Zeros of a Complex Polynomial, Comm. ACM, 15, 97-99.
 * 
 * Jenkins, M. A. (1975), Algorithm 493: Zeros of a Real Polynomial, ACM TOMS, 1, 178-189.
 * 
 * This is an iterative method that provides rapid convergence and high accuracy in most cases.
 */

class SimTK_SimTKCOMMON_EXPORT PolynomialRootFinder {
public:
    class ZeroLeadingCoefficient;
    /**
     * Find the roots of a quadratic polynomial with real coefficients.
     * 
     * @param coefficients     The polynomial coefficients in order of decreasing powers
     * @param roots            On exit, the roots of the polynomial are stored in this
     */
    template <class T>
    static void findRoots(const Vec<3,T>& coefficients, Vec<2,complex<T> >& roots);
    /**
     * Find the roots of a quadratic polynomial with complex coefficients.
     * 
     * @param coefficients     The polynomial coefficients in order of decreasing powers
     * @param roots            On exit, the roots of the polynomial are stored in this
     */
    template <class T>
    static void findRoots(const Vec<3,complex<T> >& coefficients, Vec<2,complex<T> >& roots);
    /**
     * Find the roots of a cubic polynomial with real coefficients.
     * 
     * @param coefficients     The polynomial coefficients in order of decreasing powers
     * @param roots            On exit, the roots of the polynomial are stored in this
     */
    template <class T>
    static void findRoots(const Vec<4,T>& coefficients, Vec<3,complex<T> >& roots);
    /**
     * Find the roots of a cubic polynomial with complex coefficients.
     * 
     * @param coefficients     The polynomial coefficients in order of decreasing powers
     * @param roots            On exit, the roots of the polynomial are stored in this
     */
    template <class T>
    static void findRoots(const Vec<4,complex<T> >& coefficients, Vec<3,complex<T> >& roots);
    /**
     * Find the roots of a polynomial of arbitrary degree with real coefficients.
     * 
     * @param coefficients     The polynomial coefficients in order of decreasing powers
     * @param roots            On exit, the roots of the polynomial are stored in this
     */
    template <class T>
    static void findRoots(const Vector_<T>& coefficients, Vector_<complex<T> >& roots);
    /**
     * Find the roots of a polynomial of arbitrary degree with complex coefficients.
     * 
     * @param coefficients     The polynomial coefficients in order of decreasing powers
     * @param roots            On exit, the roots of the polynomial are stored in this
     */
    template <class T>
    static void findRoots(const Vector_<complex<T> >& coefficients, Vector_<complex<T> >& roots);
};

/**
 * This is an exception which is thrown by all of the PolynomialRootFinder::findRoots() methods.
 * It indicates that the leading polynomial coefficient is zero.  This means that the polynomial
 * is not really of the stated degree, and does not have the expected number of roots.
 */

class PolynomialRootFinder::ZeroLeadingCoefficient : public Exception::Base {
public:
    ZeroLeadingCoefficient(const char* fn, int ln) : Base(fn,ln) {
        setMessage("Attempting to find roots of a polynomial whose leading coefficient is 0.");
    }
    virtual ~ZeroLeadingCoefficient() throw() { }
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_POLYNOMIALROOTFINDER_H_
