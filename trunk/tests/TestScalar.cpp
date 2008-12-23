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
    const float scale = std::max(std::max(std::abs(v1), std::abs(v2)), 0.1f);
    ASSERT(std::abs(v1-v2) < scale*NTraits<float>::getSignificant());
}
void assertEqual(double v1, double v2) {
    const double scale = std::max(std::max(std::abs(v1), std::abs(v2)), 0.1);
    ASSERT(std::abs(v1-v2) < scale*NTraits<double>::getSignificant());
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
    const negator<float>&           nflt   = reinterpret_cast<const negator<float>&>(fltRegular);
    const negator<double>&          ndbl   = reinterpret_cast<const negator<double>&>(dblRegular);
    negator<std::complex<float> >&   ncflt  = reinterpret_cast<negator<std::complex<float> >&> (cflt);
    negator<std::complex<double> >&  ncdbl  = reinterpret_cast<negator<std::complex<double> >&>(cdbl);
    negator<conjugate<float> >&      ncjflt = reinterpret_cast<negator<conjugate<float> >&>    (cjflt);
    negator<conjugate<double> >&     ncjdbl = reinterpret_cast<negator<conjugate<double> >&>   (cjdbl);

    // Test that negators are working properly.
    assertEqual(nflt, -fltRegular);
    assertEqual(ndbl, -dblRegular);
    assertEqual(ncflt, -cflt);
    assertEqual(-ncflt, cflt);
    assertEqual(ncjflt, -cjflt);
    assertEqual(-ncjflt, cjflt);

    ASSERT(!isNaN(nflt) && !isNaN(ndbl));
    ASSERT(!isNaN(ncflt) && !isNaN(ncdbl));
    ASSERT(!isNaN(ncjflt) && !isNaN(ncjdbl));

    // Should be NaN if either or both parts are NaN.
    cflt = std::complex<float>(cflt.real(), fltNaN);
    cdbl = std::complex<double>(cdbl.real(), dblNaN);
    cjflt = conjugate<float>(cjflt.real(), fltNaN);
    cjdbl = conjugate<double>(cjdbl.real(), dblNaN);

    // Imaginary only is NaN.
    ASSERT(isNaN(cflt) && isNaN(cdbl));
    ASSERT(isNaN(cjflt) && isNaN(cjdbl));
    ASSERT(isNaN(ncflt) && isNaN(ncdbl));
    ASSERT(isNaN(ncjflt) && isNaN(ncjdbl));

    cflt = std::complex<float>(fltNaN, cflt.imag());
    cdbl = std::complex<double>(dblNaN, cdbl.imag());
    cjflt = conjugate<float>(fltNaN, cjflt.imag());
    cjdbl = conjugate<double>(dblNaN, cjdbl.imag());

    // Both parts are NaN.
    ASSERT(isNaN(cflt) && isNaN(cdbl));
    ASSERT(isNaN(cjflt) && isNaN(cjdbl));
    ASSERT(isNaN(ncflt) && isNaN(ncdbl));
    ASSERT(isNaN(ncjflt) && isNaN(ncjdbl));

    // Restore imaginary part to normal.
    cflt = std::complex<float>(cflt.real(), fltRegular);
    cdbl = std::complex<double>(cdbl.real(), dblRegular);
    cjflt = conjugate<float>(cjflt.real(), fltRegular);
    cjdbl = conjugate<double>(cjdbl.real(), dblRegular);

    // Real part only is NaN;
    ASSERT(isNaN(cflt) && isNaN(cdbl));
    ASSERT(isNaN(cjflt) && isNaN(cjdbl));
    ASSERT(isNaN(ncflt) && isNaN(ncdbl));
    ASSERT(isNaN(ncjflt) && isNaN(ncjdbl));
}

void testIsInf() {
    const float  fltRegular = -12.34f;
    const double dblRegular = -12.34;
    const float fltInf = NTraits<float>::getInfinity();
    const double dblInf = NTraits<double>::getInfinity();
    const float mfltInf = -fltInf;
    const double mdblInf = -dblInf;
    const negator<float>& nfltInf = reinterpret_cast<const negator<float>&>(fltInf);
    const negator<double>& ndblInf = reinterpret_cast<const negator<double>&>(dblInf);

    ASSERT(nfltInf == -fltInf);
    ASSERT(ndblInf == -dblInf);

    ASSERT(isInf(fltInf) && isInf(dblInf));
    ASSERT(isInf(mfltInf) && isInf(mdblInf));
    ASSERT(isInf(nfltInf) && isInf(ndblInf));
    ASSERT(!isInf(fltRegular) && !isInf(dblRegular));

    std::complex<float> cflt(fltRegular, -2*fltRegular);
    std::complex<double> cdbl(dblRegular, -2*dblRegular);
    conjugate<float> cjflt(fltRegular, -2*fltRegular);
    conjugate<double> cjdbl(dblRegular, -2*dblRegular);

    ASSERT(!isInf(cflt) && !isInf(cdbl));
    ASSERT(!isInf(cjflt) && !isInf(cjdbl));

    // Reference the same memory as a negator of its contents.
    const negator<float>&           nflt   = reinterpret_cast<const negator<float>&>(fltRegular);
    const negator<double>&          ndbl   = reinterpret_cast<const negator<double>&>(dblRegular);
    negator<std::complex<float> >&   ncflt  = reinterpret_cast<negator<std::complex<float> >&> (cflt);
    negator<std::complex<double> >&  ncdbl  = reinterpret_cast<negator<std::complex<double> >&>(cdbl);
    negator<conjugate<float> >&      ncjflt = reinterpret_cast<negator<conjugate<float> >&>    (cjflt);
    negator<conjugate<double> >&     ncjdbl = reinterpret_cast<negator<conjugate<double> >&>   (cjdbl);

    // Test that negators are working properly.
    assertEqual(nflt, -fltRegular);
    assertEqual(ndbl, -dblRegular);
    assertEqual(ncflt, -cflt);
    assertEqual(-ncflt, cflt);
    assertEqual(ncjflt, -cjflt);
    assertEqual(-ncjflt, cjflt);

    ASSERT(!isInf(nflt) && !isInf(ndbl));
    ASSERT(!isInf(ncflt) && !isInf(ncdbl));
    ASSERT(!isInf(ncjflt) && !isInf(ncjdbl));

    // Should be Inf if either or both parts are Inf, as long as neither
    // part is NaN.
    cflt = std::complex<float>(cflt.real(), fltInf);
    cdbl = std::complex<double>(cdbl.real(), dblInf);
    cjflt = conjugate<float>(cjflt.real(), fltInf);
    cjdbl = conjugate<double>(cjdbl.real(), dblInf);

    // Imaginary only is Inf.
    ASSERT(isInf(cflt) && isInf(cdbl));
    ASSERT(isInf(cjflt) && isInf(cjdbl));
    ASSERT(isInf(ncflt) && isInf(ncdbl));
    ASSERT(isInf(ncjflt) && isInf(ncjdbl));

    cflt = std::complex<float>(fltInf, cflt.imag());
    cdbl = std::complex<double>(dblInf, cdbl.imag());
    cjflt = conjugate<float>(fltInf, cjflt.imag());
    cjdbl = conjugate<double>(dblInf, cjdbl.imag());

    // Both parts are Inf.
    ASSERT(isInf(cflt) && isInf(cdbl));
    ASSERT(isInf(cjflt) && isInf(cjdbl));
    ASSERT(isInf(ncflt) && isInf(ncdbl));
    ASSERT(isInf(ncjflt) && isInf(ncjdbl));

    // Restore imaginary part to normal.
    cflt = std::complex<float>(cflt.real(), fltRegular);
    cdbl = std::complex<double>(cdbl.real(), dblRegular);
    cjflt = conjugate<float>(cjflt.real(), fltRegular);
    cjdbl = conjugate<double>(cjdbl.real(), dblRegular);

    // Real part only is Inf;
    ASSERT(isInf(cflt) && isInf(cdbl));
    ASSERT(isInf(cjflt) && isInf(cjdbl));
    ASSERT(isInf(ncflt) && isInf(ncdbl));
    ASSERT(isInf(ncjflt) && isInf(ncjdbl));

    // Set real part to minus infinity.
    cflt = std::complex<float>(mfltInf, cflt.imag());
    cdbl = std::complex<double>(mdblInf, cdbl.imag());
    cjflt = conjugate<float>(mfltInf, cjflt.imag());
    cjdbl = conjugate<double>(mdblInf, cjdbl.imag());

    ASSERT(isInf(cflt) && isInf(cdbl));
    ASSERT(isInf(cjflt) && isInf(cjdbl));
    ASSERT(isInf(ncflt) && isInf(ncdbl));
    ASSERT(isInf(ncjflt) && isInf(ncjdbl));

    // Set real part to NaN.
    const float fltNaN = NTraits<float>::getNaN();
    const double dblNaN = NTraits<double>::getNaN();
    cflt = std::complex<float>(fltNaN, cflt.imag());
    cdbl = std::complex<double>(dblNaN, cdbl.imag());
    cjflt = conjugate<float>(fltNaN, cjflt.imag());
    cjdbl = conjugate<double>(dblNaN, cjdbl.imag());

    ASSERT(!isInf(cflt) && !isInf(cdbl));
    ASSERT(!isInf(cjflt) && !isInf(cjdbl));
    ASSERT(!isInf(ncflt) && !isInf(ncdbl));
    ASSERT(!isInf(ncjflt) && !isInf(ncjdbl));
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
    const float mfltInf = -fltInf;
    const double mdblInf = -dblInf;

    ASSERT(isFinite(fltRegular) && isFinite(dblRegular));
    ASSERT(!isFinite(fltNaN) && !isFinite(dblNaN));
    ASSERT(!isFinite(fltInf) && !isFinite(dblInf));
    ASSERT(!isFinite(mfltInf) && !isFinite(mdblInf));

    std::complex<float> cflt(fltRegular, -2*fltRegular);
    std::complex<double> cdbl(dblRegular, -2*dblRegular);
    conjugate<float> cjflt(fltRegular, -2*fltRegular);
    conjugate<double> cjdbl(dblRegular, -2*dblRegular);

    ASSERT(isFinite(cflt) && isFinite(cdbl));
    ASSERT(isFinite(cjflt) && isFinite(cjdbl));

    // Reference the same memory as a negator of its contents.
    const negator<float>&           nflt   = reinterpret_cast<const negator<float>&>(fltRegular);
    const negator<double>&          ndbl   = reinterpret_cast<const negator<double>&>(dblRegular);
    negator<std::complex<float> >&   ncflt  = reinterpret_cast<negator<std::complex<float> >&> (cflt);
    negator<std::complex<double> >&  ncdbl  = reinterpret_cast<negator<std::complex<double> >&>(cdbl);
    negator<conjugate<float> >&      ncjflt = reinterpret_cast<negator<conjugate<float> >&>    (cjflt);
    negator<conjugate<double> >&     ncjdbl = reinterpret_cast<negator<conjugate<double> >&>   (cjdbl);

    // Test that negators are working properly.
    assertEqual(nflt, -fltRegular);
    assertEqual(ndbl, -dblRegular);
    assertEqual(ncflt, -cflt);
    assertEqual(-ncflt, cflt);
    assertEqual(ncjflt, -cjflt);
    assertEqual(-ncjflt, cjflt);

    ASSERT(isFinite(nflt) && isFinite(ndbl));
    ASSERT(isFinite(ncflt) && isFinite(ncdbl));
    ASSERT(isFinite(ncjflt) && isFinite(ncjdbl));

    // Should be finite only if both parts are finite.
    cflt = std::complex<float>(cflt.real(),  fltInf);
    cdbl = std::complex<double>(cdbl.real(), mdblInf);
    cjflt = conjugate<float>(cjflt.real(),   fltNaN);
    cjdbl = conjugate<double>(cjdbl.real(),  dblInf);

    // Imaginary only is NaN.
    ASSERT(!isFinite(cflt) && !isFinite(cdbl));
    ASSERT(!isFinite(cjflt) && !isFinite(cjdbl));
    ASSERT(!isFinite(ncflt) && !isFinite(ncdbl));
    ASSERT(!isFinite(ncjflt) && !isFinite(ncjdbl));

    cflt = std::complex<float> (fltInf, cflt.imag());
    cdbl = std::complex<double>(mdblInf, cdbl.imag());
    cjflt = conjugate<float>   (fltNaN, cjflt.imag());
    cjdbl = conjugate<double>  (dblInf, cjdbl.imag());

    // Both parts are non-finite.
    ASSERT(!isFinite(cflt) && !isFinite(cdbl));
    ASSERT(!isFinite(cjflt) && !isFinite(cjdbl));
    ASSERT(!isFinite(ncflt) && !isFinite(ncdbl));
    ASSERT(!isFinite(ncjflt) && !isFinite(ncjdbl));

    // Restore imaginary part to normal.
    cflt = std::complex<float>(cflt.real(), fltRegular);
    cdbl = std::complex<double>(cdbl.real(), dblRegular);
    cjflt = conjugate<float>(cjflt.real(), fltRegular);
    cjdbl = conjugate<double>(cjdbl.real(), dblRegular);

    // Real part only is non-finite;
    ASSERT(!isFinite(cflt) && !isFinite(cdbl));
    ASSERT(!isFinite(cjflt) && !isFinite(cjdbl));
    ASSERT(!isFinite(ncflt) && !isFinite(ncdbl));
    ASSERT(!isFinite(ncjflt) && !isFinite(ncjdbl));
}

void testSignBit() {
    const unsigned char ucm=0xff, ucz=0, ucp=27;
    const unsigned short usm=0xffff, usz=0, usp=2342;
    const unsigned int   uim=0xffffffff, uiz=0, uip=2342344;
    const unsigned long ulm=(unsigned long)-23423L, ulz=0, ulp=234234UL;
    const unsigned long long ullm=(unsigned long long)-234234234LL, ullz=0, ullp=234234234ULL;

    ASSERT(!(signBit(ucm)||signBit(ucz)||signBit(ucp)));
    ASSERT(!(signBit(usm)||signBit(usz)||signBit(usp)));
    ASSERT(!(signBit(uim)||signBit(uiz)||signBit(uip)));
    ASSERT(!(signBit(ulm)||signBit(ulz)||signBit(ulp)));
    ASSERT(!(signBit(ullm)||signBit(ullz)||signBit(ullp)));
    
    // Note that signBit(char) doesn't exist.

    const signed char cm=-23, cz=0, cp=99;
    const short sm=-1234, sz=0, sp=23423;
    const int im=-2342343, iz=0, ip=29472383;
    const long lm=-43488, lz=0, lp=3454545;
    const long long llm=-2342342343433LL, llz=0, llp=874578478478574LL;

    ASSERT(signBit(cm) && !(signBit(cz)||signBit(cp)));
    ASSERT(signBit(sm) && !(signBit(sz)||signBit(sp)));
    ASSERT(signBit(im) && !(signBit(iz)||signBit(ip)));
    ASSERT(signBit(lm) && !(signBit(lz)||signBit(lp)));
    ASSERT(signBit(llm) && !(signBit(llz)||signBit(llp)));

    const float fm=-12398.34f, fz=0, fp=4354.331f;
    const double dm=-234234.454, dz=0, dp=345345.2342;
    float mfz=-fz; double mdz=-dz;// -0

    ASSERT(signBit(fm) && !(signBit(fz)||signBit(fp)));
    ASSERT(signBit(dm) && !(signBit(dz)||signBit(dp)));
    ASSERT(signBit(mfz) && signBit(mdz));

    // Note: signBit of negated float or double should be the
    // *same* as the underlying float or double; it is the
    // interpretation of that bit that is supposed to be 
    // different.
    const negator<float>& nfm=reinterpret_cast<const negator<float>&>(fm);
    const negator<float>& nfz=reinterpret_cast<const negator<float>&>(fz);
    const negator<float>& nfp=reinterpret_cast<const negator<float>&>(fp);
    const negator<float>& nmfz=reinterpret_cast<const negator<float>&>(mfz);
    const negator<double>& ndm=reinterpret_cast<const negator<double>&>(dm);
    const negator<double>& ndz=reinterpret_cast<const negator<double>&>(dz);
    const negator<double>& ndp=reinterpret_cast<const negator<double>&>(dp);
    const negator<double>& nmdz=reinterpret_cast<const negator<double>&>(mdz);

    ASSERT(signBit(nfm) && !(signBit(nfz)||signBit(nfp)));
    ASSERT(signBit(ndm) && !(signBit(ndz)||signBit(ndp)));
    ASSERT(signBit(nmfz) && signBit(nmdz));

    const float fltInf = NTraits<float>::getInfinity();
    const double dblInf = NTraits<double>::getInfinity();
    const float mfltInf = -fltInf;
    const double mdblInf = -dblInf;

    ASSERT(!signBit(fltInf) && !signBit(dblInf));
    ASSERT(signBit(mfltInf) && signBit(mdblInf));
}


void testSign() {
    const unsigned char ucm=0xff, ucz=0, ucp=27;
    const unsigned short usm=0xffff, usz=0, usp=2342;
    const unsigned int   uim=0xffffffff, uiz=0, uip=2342344;
    const unsigned long ulm=(unsigned long)-23423L, ulz=0, ulp=234234UL;
    const unsigned long long ullm=(unsigned long long)-234234234LL, ullz=0, ullp=234234234ULL;

    ASSERT(sign(ucm)==1 && sign(ucz)==0 && sign(ucp)==1);
    ASSERT(sign(usm)==1 && sign(usz)==0 && sign(usp)==1);
    ASSERT(sign(uim)==1 && sign(uiz)==0 && sign(uip)==1);
    ASSERT(sign(ulm)==1 && sign(ulz)==0 && sign(ulp)==1);
    ASSERT(sign(ullm)==1 && sign(ullz)==0 && sign(ullp)==1);

    // Note that sign(char) doesn't exist.

    const signed char cm=-23, cz=0, cp=99;
    const short sm=-1234, sz=0, sp=23423;
    const int im=-2342343, iz=0, ip=29472383;
    const long lm=-43488, lz=0, lp=3454545;
    const long long llm=-2342342343433LL, llz=0, llp=874578478478574LL;

    ASSERT(sign(cm)==-1 && sign(cz)==0 && sign(cp)==1);
    ASSERT(sign(sm)==-1 && sign(sz)==0 && sign(sp)==1);
    ASSERT(sign(im)==-1 && sign(iz)==0 && sign(ip)==1);
    ASSERT(sign(lm)==-1 && sign(lz)==0 && sign(lp)==1);
    ASSERT(sign(llm)==-1 && sign(llz)==0 && sign(llp)==1);

    const float fm=-12398.34f, fz=0, fp=4354.331f;
    const double dm=-234234.454, dz=0, dp=345345.2342;
    float mfz=-fz; double mdz=-dz;// -0

    ASSERT(sign(fm)==-1 && sign(fz)==0 && sign(fp)==1);
    ASSERT(sign(dm)==-1 && sign(dz)==0 && sign(dp)==1);
    ASSERT(sign(mfz)==0 && sign(mdz)==0); // doesn't matter if it's -0

    // Note: sign of negated float or double should be the
    // *opposite* as the underlying float or double.
    const negator<float>& nfm=reinterpret_cast<const negator<float>&>(fm);
    const negator<float>& nfz=reinterpret_cast<const negator<float>&>(fz);
    const negator<float>& nfp=reinterpret_cast<const negator<float>&>(fp);
    const negator<float>& nmfz=reinterpret_cast<const negator<float>&>(mfz);
    const negator<double>& ndm=reinterpret_cast<const negator<double>&>(dm);
    const negator<double>& ndz=reinterpret_cast<const negator<double>&>(dz);
    const negator<double>& ndp=reinterpret_cast<const negator<double>&>(dp);
    const negator<double>& nmdz=reinterpret_cast<const negator<double>&>(mdz);

    ASSERT(sign(nfm)==1 && sign(nfz)==0 && sign(nfp)==-1);
    ASSERT(sign(ndm)==1 && sign(ndz)==0 && sign(ndp)==-1);
    ASSERT(sign(nmfz)==0 && sign(nmdz)==0); // doesn't matter if it's -0

    const float fltInf = NTraits<float>::getInfinity();
    const double dblInf = NTraits<double>::getInfinity();
    const float mfltInf = -fltInf;
    const double mdblInf = -dblInf;
    const negator<float>& nfltInf = reinterpret_cast<const negator<float>&>(fltInf);
    const negator<double>& ndblInf = reinterpret_cast<const negator<double>&>(dblInf);

    ASSERT(sign(fltInf)==1 && sign(dblInf)==1);
    ASSERT(sign(mfltInf)==-1 && sign(mdblInf)==-1);
    ASSERT(sign(nfltInf)==-1 && sign(ndblInf)==-1);
}

void testSquareAndCube() {
    const float fval = -23.33f;
    const double dval = -234443.441;
    const negator<float>& nfval = reinterpret_cast<const negator<float>&>(fval);
    const negator<double>& ndval = reinterpret_cast<const negator<double>&>(dval);

    // Basic test.
    assertEqual(square(fval), fval*fval);
    assertEqual(square(dval), dval*dval);
    assertEqual(cube(fval), fval*fval*fval);
    assertEqual(cube(dval), dval*dval*dval);

    // Test scalar negators.
    assertEqual(square(nfval), nfval*nfval);
    assertEqual(square(nfval), fval*fval);
    assertEqual(square(ndval), ndval*ndval);
    assertEqual(square(ndval), dval*dval);
    assertEqual(cube(nfval), nfval*nfval*nfval);
    assertEqual(cube(nfval), -fval*fval*fval);
    assertEqual(cube(ndval), ndval*ndval*ndval);
    assertEqual(cube(ndval), -dval*dval*dval);

    // Create complex and conjugate values.

    std::complex<float> fc(-234.343f, 45345e7f);
    std::complex<double> dc(-234.343, 45345e7);
    conjugate<float> fcj(-19.1e3f, -454.234f);
    conjugate<double> dcj(-19.1e3, -454.234);

    // Manual conjugates
    std::complex<float>  fcmj(fcj.real(), fcj.imag());
    std::complex<double> dcmj(dcj.real(), dcj.imag());
    ASSERT(fcj == fcmj);
    ASSERT(dcj == dcmj);
    ASSERT(fcj*fcj == fcmj*fcmj);
    ASSERT(dcj*dcj == dcmj*dcmj);
    ASSERT(fcj*fcj*fcj == fcmj*fcmj*fcmj);
    ASSERT(dcj*dcj*dcj == dcmj*dcmj*dcmj);

    // Negators of complex an conjugate.
    negator<std::complex<float> >&   nfc  = reinterpret_cast<negator<std::complex<float> >&> (fc);
    negator<std::complex<double> >&  ndc  = reinterpret_cast<negator<std::complex<double> >&>(dc);
    negator<conjugate<float> >&      nfcj = reinterpret_cast<negator<conjugate<float> >&>    (fcj);
    negator<conjugate<double> >&     ndcj = reinterpret_cast<negator<conjugate<double> >&>   (dcj);

    ASSERT(nfc == -fc);
    ASSERT(ndc == -dc);
    ASSERT(nfcj == -fcj);
    ASSERT(ndcj == -dcj);


    // Basic complex and conjugate tests.
    assertEqual(square(fc), fc*fc);
    assertEqual(cube(fc), fc*fc*fc);
    assertEqual(square(dc), dc*dc);
    assertEqual(cube(dc), dc*dc*dc);
    assertEqual(square(fcj), fcj*fcj);
    assertEqual(cube(fcj), fcj*fcj*fcj);
    assertEqual(square(dcj), dcj*dcj);
    assertEqual(cube(dcj), dcj*dcj*dcj);

    // Tests involving negators of complex and conjugate.
    assertEqual(square(nfc), nfc*nfc); 
    assertEqual(square(nfc), fc*fc);
    assertEqual(square(ndc), ndc*ndc);
    assertEqual(square(ndc), dc*dc);

    assertEqual(cube(nfc), nfc*nfc*nfc); 
    assertEqual(cube(nfc), -fc*fc*fc);
    assertEqual(cube(ndc), ndc*ndc*ndc);
    assertEqual(cube(ndc), -dc*dc*dc);

    assertEqual(square(nfcj), nfcj*nfcj); 
    assertEqual(square(nfcj), fcj*fcj);
    assertEqual(square(ndcj), ndcj*ndcj);
    assertEqual(square(ndcj), dcj*dcj);

    assertEqual(cube(nfcj), nfcj*nfcj*nfcj); 
    assertEqual(cube(nfcj), -fcj*fcj*fcj);
    assertEqual(cube(ndcj), ndcj*ndcj*ndcj);
    assertEqual(cube(ndcj), -dcj*dcj*dcj);
}

int main() {
    try {
        testIsNaN();
        testIsInf();
        testIsFinite();
        testSignBit();
        testSign();
        testSquareAndCube();
    } catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}

