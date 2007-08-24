/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
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

#include "SimTKcommon/PolynomialRootFinder.h"
#include "SimTKcommon/internal/rpoly.h"
#include "SimTKcommon/internal/cpoly.h"

namespace SimTK {

/**
 * Find the roots of a quadratic polynomial with real coefficients.
 * 
 * @param coefficients     The polynomial coefficients in order of decreasing powers
 * @param roots            On exit, the roots of the polynomial are stored in this
 */

template <class T>
void PolynomialRootFinder::findRoots(const Vec<3,T>& coefficients, Vec<2,complex<T> >& roots) {
    T a = coefficients[0], b = coefficients[1], c = coefficients[2];
    if (a == 0.0)
        SimTK_THROW(ZeroLeadingCoefficient);
    T b2 = b*b;
    T discriminant = b2 - 4.0*a*c;
    T tol = 2.0*NTraits<T>::getEps()*b2;
    if (discriminant < tol && discriminant > -tol) {

        // b^2 == 4ac to within machine precision, so make the roots identical.
        
        T root = -b/(2.0*a);
        roots[0] = complex<T>(root, 0.0);
        roots[1] = complex<T>(root, 0.0);
        return;
    }
    if (b == 0.0) {
        
        // The coefficient of the linear term is zero, which makes the formula simpler.
        
        if (discriminant >= 0.0) {
            T root = std::sqrt(discriminant)/2.0*a;
            roots[0] = root;
            roots[1] = -root;
        }
        else {
            T root = sqrt(-discriminant)/2.0*a;
            roots[0] = Complex(0.0, root);
            roots[1] = Complex(0.0, -root);
        }
        return;
    }
    complex<T> q = ((T) -0.5)*(b+(b > 0.0 ? std::sqrt(complex<T>(discriminant)) : -std::sqrt(complex<T>(discriminant))));
    roots[0] = q/a;
    roots[1] = c/q;
}

/**
 * Find the roots of a quadratic polynomial with complex coefficients.
 * 
 * @param coefficients     The polynomial coefficients in order of decreasing powers
 * @param roots            On exit, the roots of the polynomial are stored in this
 */

template <class T>
void PolynomialRootFinder::findRoots(const Vec<3,complex<T> >& coefficients, Vec<2,complex<T> >& roots) {
    complex<T> a = coefficients[0], b = coefficients[1], c = coefficients[2];
    if (a == (T) 0.0)
        SimTK_THROW(ZeroLeadingCoefficient);
    complex<T> b2 = b*b;
    complex<T> discriminant = b2 - ((T) 4.0)*a*c;
    if (b == (T) 0.0) {
        
        // The coefficient of the linear term is zero, which makes the formula simpler.
        
        complex<T> root = std::sqrt(discriminant)/((T) 2.0)*a;
        roots[0] = root;
        roots[1] = -root;
        return;
    }
    T temp = (conj(b)*sqrt(discriminant)).real();
    complex<T> q = ((T) -0.5)*(b+(temp > 0.0 ? std::sqrt(discriminant) : -std::sqrt(discriminant)));
    roots[0] = q/a;
    roots[1] = c/q;
}

/**
 * Find the roots of a cubic polynomial with real coefficients.
 * 
 * @param coefficients     The polynomial coefficients in order of decreasing powers
 * @param roots            On exit, the roots of the polynomial are stored in this
 */

template <class T>
void PolynomialRootFinder::findRoots(const Vec<4,T>& coefficients, Vec<3,complex<T> >& roots) {
    if (coefficients[0] == 0.0)
        SimTK_THROW(ZeroLeadingCoefficient);
    T coeff[4] = {coefficients[0], coefficients[1], coefficients[2], coefficients[3]};
    T rootr[3];
    T rooti[3];
    RPoly<T>().findRoots(coeff, 3, rootr, rooti);
    roots[0] = Complex(rootr[0], rooti[0]);
    roots[1] = Complex(rootr[1], rooti[1]);
    roots[2] = Complex(rootr[2], rooti[2]);
    
    
    
//    T scale = 1.0/coefficients[0];
//    T a = scale*coefficients[1], b = scale*coefficients[2], c = scale*coefficients[3];
//    T q = (a*a-3.0*b)/9.0;
//    T r = (2.0*a*a*a - 9.0*a*b + 27.0*c)/54.0;
//    T r2 = r*r;
//    T q3 = q*q*q;
//    T diff = r2-q3;
//    T tol = 2.0*NTraits<T>::getEps()*r2;
//    if (diff < tol && diff > -tol) {
//
//        // r^2 == q^3 to within machine precision, so set theta=0 and simplify the formulas.
//        
//        T mult = -2.0*std::sqrt(q);
//        T sub = a/3.0;
//        T twopi = 2.0*NTraits<T>::getPi();
//        roots[0] = mult - sub;
//        T root = mult*std::cos(twopi/3.0) - sub;
//        roots[1] = root;
//        roots[2] = root;
//        return;
//    }
//    if (r2 < q3) {
//        
//        // There are three real roots.
//        
//        T theta = std::acos(r/std::sqrt(q3));
//        T mult = -2.0*std::sqrt(q);
//        T sub = a/3.0;
//        T twopi = 2.0*NTraits<T>::getPi();
//        roots[0] = mult*std::cos(theta/3.0) - sub;
//        roots[1] = mult*std::cos((theta+twopi)/3.0) - sub;
//        roots[2] = mult*std::cos((theta-twopi)/3.0) - sub;
//        return;
//    }
//    T aa = std::pow(std::abs(r)+std::sqrt(diff), 1.0/3.0);
//    if (r >= 0.0)
//        aa = -aa;
//    T bb = (aa == 0.0 ? 0.0 : q/aa);
//    roots[0] = aa+bb-a/3.0;
//    T rootr = -0.5*(aa+bb)-a/3.0;
//    T rooti = 0.5*NTraits<T>::getSqrt3()*(aa-bb);
//    roots[1] = complex<T>(rootr, rooti);
//    roots[2] = complex<T>(rootr, -rooti);
}


/**
 * Find the roots of a cubic polynomial with complex coefficients.
 * 
 * @param coefficients     The polynomial coefficients in order of decreasing powers
 * @param roots            On exit, the roots of the polynomial are stored in this
 */

template <class T>
void PolynomialRootFinder::findRoots(const Vec<4,complex<T> >& coefficients, Vec<3,complex<T> >& roots) {
    if (coefficients[0] == (T) 0.0)
        SimTK_THROW(ZeroLeadingCoefficient);
    T coeffr[4] = {coefficients[0].real(), coefficients[1].real(), coefficients[2].real(), coefficients[3].real()};
    T coeffi[4] = {coefficients[0].imag(), coefficients[1].imag(), coefficients[2].imag(), coefficients[3].imag()};
    T rootr[3];
    T rooti[3];
    CPoly<T>().findRoots(coeffr, coeffi, 3, rootr, rooti);
    roots[0] = Complex(rootr[0], rooti[0]);
    roots[1] = Complex(rootr[1], rooti[1]);
    roots[2] = Complex(rootr[2], rooti[2]);
}

/**
 * Find the roots of a polynomial of arbitrary degree with real coefficients.
 * 
 * @param coefficients     The polynomial coefficients in order of decreasing powers
 * @param roots            On exit, the roots of the polynomial are stored in this
 */

template <class T>
void PolynomialRootFinder::findRoots(const Vector_<T>& coefficients, Vector_<complex<T> >& roots) {
    if (coefficients[0] == 0.0)
        SimTK_THROW(ZeroLeadingCoefficient);
    int n = roots.size();
    assert(coefficients.size() == n+1);
    T *coeff = new T[n+1];
    T *rootr = new T[n];
    T *rooti = new T[n];
    try {
        for (int i = 0; i < n+1; ++i)
            coeff[i] = coefficients[i];
        RPoly<T>().findRoots(coeff, n, rootr, rooti);
        for (int i = 0; i < n; ++i)
            roots[i] = Complex(rootr[i], rooti[i]);
    }
    catch (...) {
        delete[] coeff;
        delete[] rootr;
        delete[] rooti;
        throw;
    }
}

/**
 * Find the roots of a polynomial of arbitrary degree with complex coefficients.
 * 
 * @param coefficients     The polynomial coefficients in order of decreasing powers
 * @param roots            On exit, the roots of the polynomial are stored in this
 */

template <class T>
void PolynomialRootFinder::findRoots(const Vector_<complex<T> >& coefficients, Vector_<complex<T> >& roots) {
    if (coefficients[0] == (T) 0.0)
        SimTK_THROW(ZeroLeadingCoefficient);
    int n = roots.size();
    assert(coefficients.size() == n+1);
    T *coeffr = new T[n+1];
    T *coeffi = new T[n+1];
    T *rootr = new T[n];
    T *rooti = new T[n];
    try {
        for (int i = 0; i < n+1; ++i) {
            coeffr[i] = coefficients[i].real();
            coeffi[i] = coefficients[i].imag();
        }
        CPoly<T>().findRoots(coeffr, coeffi, n, rootr, rooti);
        for (int i = 0; i < n; ++i)
            roots[i] = Complex(rootr[i], rooti[i]);
    }
    catch (...) {
        delete[] coeffr;
        delete[] coeffi;
        delete[] rootr;
        delete[] rooti;
        throw;
    }
}

template void PolynomialRootFinder::findRoots(const Vec<3,float>& coefficients, Vec<2,complex<float> >& roots);
template void PolynomialRootFinder::findRoots(const Vec<3,complex<float> >& coefficients, Vec<2,complex<float> >& roots);
template void PolynomialRootFinder::findRoots(const Vec<4,float>& coefficients, Vec<3,complex<float> >& roots);
template void PolynomialRootFinder::findRoots(const Vec<4,complex<float> >& coefficients, Vec<3,complex<float> >& roots);
template void PolynomialRootFinder::findRoots(const Vector_<float>& coefficients, Vector_<complex<float> >& roots);
template void PolynomialRootFinder::findRoots(const Vector_<complex<float> >& coefficients, Vector_<complex<float> >& roots);

template void PolynomialRootFinder::findRoots(const Vec<3,double>& coefficients, Vec<2,complex<double> >& roots);
template void PolynomialRootFinder::findRoots(const Vec<3,complex<double> >& coefficients, Vec<2,complex<double> >& roots);
template void PolynomialRootFinder::findRoots(const Vec<4,double>& coefficients, Vec<3,complex<double> >& roots);
template void PolynomialRootFinder::findRoots(const Vec<4,complex<double> >& coefficients, Vec<3,complex<double> >& roots);
template void PolynomialRootFinder::findRoots(const Vector_<double>& coefficients, Vector_<complex<double> >& roots);
template void PolynomialRootFinder::findRoots(const Vector_<complex<double> >& coefficients, Vector_<complex<double> >& roots);

template void PolynomialRootFinder::findRoots(const Vec<3,long double>& coefficients, Vec<2,complex<long double> >& roots);
template void PolynomialRootFinder::findRoots(const Vec<3,complex<long double> >& coefficients, Vec<2,complex<long double> >& roots);
template void PolynomialRootFinder::findRoots(const Vec<4,long double>& coefficients, Vec<3,complex<long double> >& roots);
template void PolynomialRootFinder::findRoots(const Vec<4,complex<long double> >& coefficients, Vec<3,complex<long double> >& roots);
template void PolynomialRootFinder::findRoots(const Vector_<long double>& coefficients, Vector_<complex<long double> >& roots);
template void PolynomialRootFinder::findRoots(const Vector_<complex<long double> >& coefficients, Vector_<complex<long double> >& roots);

} // namespace SimTK
