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

#include "SimTKcommon/internal/PolynomialRootFinder.h"
#include "rpoly.h"
#include "cpoly.h"

namespace SimTK {

template <class T>
void PolynomialRootFinder::findRoots(const Vec<3,T>& coefficients, Vec<2,complex<T> >& roots) {
    T a = coefficients[0], b = coefficients[1], c = coefficients[2];
    if (a == 0.0)
        SimTK_THROW(ZeroLeadingCoefficient);
    T b2 = b*b;
    T discriminant = b2 - (T) 4.0*a*c;
    T tol = (T) 2.0*NTraits<T>::getEps()*b2;
    if (discriminant < tol && discriminant > -tol) {

        // b^2 == 4ac to within machine precision, so make the roots identical.
        
        T root = -b/((T) 2.0*a);
        roots[0] = complex<T>(root, 0.0);
        roots[1] = complex<T>(root, 0.0);
        return;
    }
    if (b == 0.0) {
        
        // The coefficient of the linear term is zero, which makes the formula simpler.
        
        if (discriminant >= 0.0) {
            T root = std::sqrt(discriminant)/(T) 2.0*a;
            roots[0] = root;
            roots[1] = -root;
        }
        else {
            T root = std::sqrt(-discriminant)/(T) 2.0*a;
            roots[0] = Complex(0.0, root);
            roots[1] = Complex(0.0, -root);
        }
        return;
    }
    complex<T> q = ((T) -0.5)*(b+(b > 0.0 ? std::sqrt(complex<T>(discriminant)) : -std::sqrt(complex<T>(discriminant))));
    roots[0] = q/a;
    roots[1] = c/q;
}

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
}

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
        for (int i = 0; i < n; ++i) // in case these don't get filled in
            rootr[i] = rooti[i] = NTraits<T>::getNaN(); 
        const int nrootsFound = RPoly<T>().findRoots(coeff, n, rootr, rooti);
        for (int i = 0; i < n; ++i)
            roots[i] = Complex(rootr[i], rooti[i]);
        SimTK_ERRCHK_ALWAYS(nrootsFound != -1,
            "PolynomialRootFinder::findRoots()",
            "Leading coefficient is zero; can't solve.");
        SimTK_ERRCHK1_ALWAYS(nrootsFound > 0,
            "PolynomialRootFinder::findRoots()",
            "Failure to find any roots for polynomial of order %d.", n);
    }
    catch (...) {
        delete[] coeff;
        delete[] rootr;
        delete[] rooti;
        throw;
    }
    delete[] coeff;
    delete[] rootr;
    delete[] rooti;
}

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
    delete[] coeffr;
    delete[] coeffi;
    delete[] rootr;
    delete[] rooti;
}

template SimTK_SimTKCOMMON_EXPORT void PolynomialRootFinder::findRoots<float>(const Vec<3,float>& coefficients, Vec<2,complex<float> >& roots);
template SimTK_SimTKCOMMON_EXPORT void PolynomialRootFinder::findRoots<float>(const Vec<3,complex<float> >& coefficients, Vec<2,complex<float> >& roots);
template SimTK_SimTKCOMMON_EXPORT void PolynomialRootFinder::findRoots<float>(const Vec<4,float>& coefficients, Vec<3,complex<float> >& roots);
template SimTK_SimTKCOMMON_EXPORT void PolynomialRootFinder::findRoots<float>(const Vec<4,complex<float> >& coefficients, Vec<3,complex<float> >& roots);
template SimTK_SimTKCOMMON_EXPORT void PolynomialRootFinder::findRoots<float>(const Vector_<float>& coefficients, Vector_<complex<float> >& roots);
template SimTK_SimTKCOMMON_EXPORT void PolynomialRootFinder::findRoots<float>(const Vector_<complex<float> >& coefficients, Vector_<complex<float> >& roots);

template SimTK_SimTKCOMMON_EXPORT void PolynomialRootFinder::findRoots<double>(const Vec<3,double>& coefficients, Vec<2,complex<double> >& roots);
template SimTK_SimTKCOMMON_EXPORT void PolynomialRootFinder::findRoots<double>(const Vec<3,complex<double> >& coefficients, Vec<2,complex<double> >& roots);
template SimTK_SimTKCOMMON_EXPORT void PolynomialRootFinder::findRoots<double>(const Vec<4,double>& coefficients, Vec<3,complex<double> >& roots);
template SimTK_SimTKCOMMON_EXPORT void PolynomialRootFinder::findRoots<double>(const Vec<4,complex<double> >& coefficients, Vec<3,complex<double> >& roots);
template SimTK_SimTKCOMMON_EXPORT void PolynomialRootFinder::findRoots<double>(const Vector_<double>& coefficients, Vector_<complex<double> >& roots);
template SimTK_SimTKCOMMON_EXPORT void PolynomialRootFinder::findRoots<double>(const Vector_<complex<double> >& coefficients, Vector_<complex<double> >& roots);

} // namespace SimTK
