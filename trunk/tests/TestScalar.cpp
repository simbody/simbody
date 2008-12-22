/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Michael Sherman                                                   *
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

#include <iostream>
using std::cout;
using std::endl;

using namespace SimTK;

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS((cond), "Assertion failed");}

// These are only valid for numbers near 1.
void assertEqual(float v1, float v2) {
    ASSERT(std::abs(v1-v2) < NTraits<float>::getSignificant());
}
void assertEqual(double v1, double v2) {
    ASSERT(std::abs(v1-v2) < NTraits<double>::getSignificant());
}
template <class P>
void assertEqual(const std::complex<P>& v1, const std::complex<P>& v2) {
    assertEqual(v1.real(), v2.real());
    assertEqual(v1.imag(), v2.imag());
}
template <class P>
void assertEqual(const conjugate<P>& v1, const conjugate<P>& v2) {
    assertEqual(v1.real(), v2.real());
    assertEqual(v1.imag(), v2.imag());
}
template <class P>
void assertEqual(const std::complex<P>& v1, const conjugate<P>& v2) {
    assertEqual(v1.real(), v2.real());
    assertEqual(v1.imag(), v2.imag());
}
template <class P>
void assertEqual(const conjugate<P>& v1, const std::complex<P>& v2) {
    assertEqual(v1.real(), v2.real());
    assertEqual(v1.imag(), v2.imag());
}
template <class P>
void assertEqual(const negator<P>& v1, const negator<P>& v2) {
    assertEqual(-v1, -v2);  // P, P
}
template <class P>
void assertEqual(const P& v1, const negator<P>& v2) {
    assertEqual(-v1, -v2);  // P, P
}
template <class P>
void assertEqual(const negator<P>& v1, const P& v2) {
    assertEqual(-v1, -v2);  // P, P
}
template <class P>
void assertEqual(const negator<std::complex<P> >& v1, const conjugate<P>& v2) {
    assertEqual(-v1, -v2);  // complex, conjugate
}
template <class P>
void assertEqual(const negator<conjugate<P> >& v1, const std::complex<P>& v2) {
    assertEqual(-v1, -v2);  // conjugate, complex
}
template <class P>
void assertEqual(const std::complex<P>& v1, const negator<conjugate<P> >& v2) {
    assertEqual(-v1, -v2); // complex, conjugate
}
template <class P>
void assertEqual(const conjugate<P>& v1, const negator<std::complex<P> >& v2) {
    assertEqual(-v1, -v2); // conjugate, complex
}

void testIsNaN() {
    const float  fltRegular = -12.34f;
    const double dblRegular = -12.34;
    const float fltNaN = NTraits<float>::getNaN();
    const double dblNaN = NTraits<double>::getNaN();
    const float nfltNaN = -fltNaN;
    const double ndblNaN = -dblNaN;

    ASSERT(isNaN(fltNaN) && isNaN(dblNaN));
    ASSERT(isNaN(nfltNaN) && isNaN(ndblNaN));
    ASSERT(!isNaN(fltRegular) && !isNaN(dblRegular));

    std::complex<float> cflt(fltRegular, -2*fltRegular);
    std::complex<double> cdbl(dblRegular, -2*dblRegular);
    conjugate<float> cjflt(fltRegular, -2*fltRegular);
    conjugate<double> cjdbl(dblRegular, -2*dblRegular);

    ASSERT(!isNaN(cflt) && !isNaN(cdbl));
    ASSERT(!isNaN(cjflt) && !isNaN(cjdbl));

    // Reference the same memory as a negator of its contents.
    negator<std::complex<float>>&   ncflt  = reinterpret_cast<negator<std::complex<float>>&> (cflt);
    negator<std::complex<double>>&  ncdbl  = reinterpret_cast<negator<std::complex<double>>&>(cdbl);
    negator<conjugate<float>>&      ncjflt = reinterpret_cast<negator<conjugate<float>>&>    (cjflt);
    negator<conjugate<double>>&     ncjdbl = reinterpret_cast<negator<conjugate<double>>&>   (cjdbl);

    // Test that negators are working properly.
    assertEqual(ncflt, -cflt);
    assertEqual(-ncflt, cflt);
    assertEqual(ncjflt, -cjflt);
    assertEqual(-ncjflt, cjflt);

    ASSERT(!isNaN(ncflt) && !isNaN(ncdbl));
    ASSERT(!isNaN(ncjflt) && !isNaN(ncjdbl));

    // Should be NaN if either or both parts are NaN.
    cflt.imag(fltNaN);  cdbl.imag(dblNaN);
    cjflt.imag() = fltNaN; cjdbl.imag() = dblNaN;

    // Imaginary only is NaN.
    ASSERT(isNaN(cflt) && isNaN(cdbl));
    ASSERT(isNaN(cjflt) && isNaN(cjdbl));
    ASSERT(isNaN(ncflt) && isNaN(ncdbl));
    ASSERT(isNaN(ncjflt) && isNaN(ncjdbl));

    cflt.real(fltNaN);  cdbl.real(dblNaN);
    cjflt.real() = fltNaN; cjdbl.real() = dblNaN;

    // Both parts are NaN.
    ASSERT(isNaN(cflt) && isNaN(cdbl));
    ASSERT(isNaN(cjflt) && isNaN(cjdbl));
    ASSERT(isNaN(ncflt) && isNaN(ncdbl));
    ASSERT(isNaN(ncjflt) && isNaN(ncjdbl));

    // Restore imaginary part to normal.
    cflt.imag(fltRegular); cdbl.imag(dblRegular);
    cjflt.imag() = fltRegular; cjdbl.imag() = dblRegular;

    // Real part only is NaN;
    ASSERT(isNaN(cflt) && isNaN(cdbl));
    ASSERT(isNaN(cjflt) && isNaN(cjdbl));
    ASSERT(isNaN(ncflt) && isNaN(ncdbl));
    ASSERT(isNaN(ncjflt) && isNaN(ncjdbl));
}

void testIsFinite() {
    const float  fltRegular = -12.34f;
    const double dblRegular = -12.34;
    const float fltNaN = NTraits<float>::getNaN();
    const double dblNaN = NTraits<double>::getNaN();
    const float nfltNaN = -fltNaN;
    const double ndblNaN = -dblNaN;
    const float fltInf = NTraits<float>::getInfinity();
    const double dblInf = NTraits<double>::getInfinity();
    const float nfltInf = -fltInf;
    const double ndblInf = -dblInf;

    ASSERT(isFinite(fltRegular) && isFinite(dblRegular));
    ASSERT(!isFinite(fltNaN) && !isFinite(dblNaN));
    ASSERT(!isFinite(fltInf) && !isFinite(dblInf));
    ASSERT(!isFinite(nfltInf) && !isFinite(ndblInf));

    std::complex<float> cflt(fltRegular, -2*fltRegular);
    std::complex<double> cdbl(dblRegular, -2*dblRegular);
    conjugate<float> cjflt(fltRegular, -2*fltRegular);
    conjugate<double> cjdbl(dblRegular, -2*dblRegular);

    ASSERT(isFinite(cflt) && isFinite(cdbl));
    ASSERT(isFinite(cjflt) && isFinite(cjdbl));

    // Reference the same memory as a negator of its contents.
    negator<std::complex<float>>&   ncflt  = reinterpret_cast<negator<std::complex<float>>&> (cflt);
    negator<std::complex<double>>&  ncdbl  = reinterpret_cast<negator<std::complex<double>>&>(cdbl);
    negator<conjugate<float>>&      ncjflt = reinterpret_cast<negator<conjugate<float>>&>    (cjflt);
    negator<conjugate<double>>&     ncjdbl = reinterpret_cast<negator<conjugate<double>>&>   (cjdbl);

    // Test that negators are working properly.
    assertEqual(ncflt, -cflt);
    assertEqual(-ncflt, cflt);
    assertEqual(ncjflt, -cjflt);
    assertEqual(-ncjflt, cjflt);

    ASSERT(isFinite(ncflt) && isFinite(ncdbl));
    ASSERT(isFinite(ncjflt) && isFinite(ncjdbl));

    // Should be finite only if both parts are finite.
    cflt.imag(fltInf);  cdbl.imag(nfltInf);
    cjflt.imag() = fltNaN; cjdbl.imag() = dblInf;

    // Imaginary only is NaN.
    ASSERT(!isFinite(cflt) && !isFinite(cdbl));
    ASSERT(!isFinite(cjflt) && !isFinite(cjdbl));
    ASSERT(!isFinite(ncflt) && !isFinite(ncdbl));
    ASSERT(!isFinite(ncjflt) && !isFinite(ncjdbl));

    cflt.real(fltInf);  cdbl.real(nfltInf);
    cjflt.real() = fltNaN; cjdbl.real() = dblInf;

    // Both parts are non-finite.
    ASSERT(!isFinite(cflt) && !isFinite(cdbl));
    ASSERT(!isFinite(cjflt) && !isFinite(cjdbl));
    ASSERT(!isFinite(ncflt) && !isFinite(ncdbl));
    ASSERT(!isFinite(ncjflt) && !isFinite(ncjdbl));

    // Restore imaginary part to normal.
    cflt.imag(fltRegular); cdbl.imag(dblRegular);
    cjflt.imag() = fltRegular; cjdbl.imag() = dblRegular;

    // Real part only is non-finite;
    ASSERT(!isFinite(cflt) && !isFinite(cdbl));
    ASSERT(!isFinite(cjflt) && !isFinite(cjdbl));
    ASSERT(!isFinite(ncflt) && !isFinite(ncdbl));
    ASSERT(!isFinite(ncjflt) && !isFinite(ncjdbl));
}


int main() {
    try {
        testIsNaN();
        testIsFinite();

    } catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}

