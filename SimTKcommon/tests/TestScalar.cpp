/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-15 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
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
#include "SimTKcommon/Testing.h"

#include <iostream>
#include <bitset>
#include <set>
using std::cout; using std::endl;


using namespace SimTK;


void testIsNaN() {
    const float  fltRegular = -12.34f;
    const double dblRegular = -12.34;
    const float fltNaN = NTraits<float>::getNaN();
    const double dblNaN = NTraits<double>::getNaN();
    const float nfltNaN = -fltNaN;
    const double ndblNaN = -dblNaN;

    SimTK_TEST(isNaN(fltNaN) && isNaN(dblNaN));
    SimTK_TEST(isNaN(nfltNaN) && isNaN(ndblNaN));
    SimTK_TEST(!isNaN(fltRegular) && !isNaN(dblRegular));

    std::complex<float> cflt(fltRegular, -2*fltRegular);
    std::complex<double> cdbl(dblRegular, -2*dblRegular);
    conjugate<float> cjflt(fltRegular, -2*fltRegular);
    conjugate<double> cjdbl(dblRegular, -2*dblRegular);

    SimTK_TEST(!isNaN(cflt) && !isNaN(cdbl));
    SimTK_TEST(!isNaN(cjflt) && !isNaN(cjdbl));

    // Reference the same memory as a negator of its contents.
    const negator<float>&           nflt   = reinterpret_cast<const negator<float>&>(fltRegular);
    const negator<double>&          ndbl   = reinterpret_cast<const negator<double>&>(dblRegular);
    negator<std::complex<float> >&   ncflt  = reinterpret_cast<negator<std::complex<float> >&> (cflt);
    negator<std::complex<double> >&  ncdbl  = reinterpret_cast<negator<std::complex<double> >&>(cdbl);
    negator<conjugate<float> >&      ncjflt = reinterpret_cast<negator<conjugate<float> >&>    (cjflt);
    negator<conjugate<double> >&     ncjdbl = reinterpret_cast<negator<conjugate<double> >&>   (cjdbl);

    // Test that negators are working properly.
    SimTK_TEST_EQ(nflt, -fltRegular);
    SimTK_TEST_EQ(ndbl, -dblRegular);
    SimTK_TEST_EQ(ncflt, -cflt);
    SimTK_TEST_EQ(-ncflt, cflt);
    SimTK_TEST_EQ(ncjflt, -cjflt);
    SimTK_TEST_EQ(-ncjflt, cjflt);

    SimTK_TEST(!isNaN(nflt) && !isNaN(ndbl));
    SimTK_TEST(!isNaN(ncflt) && !isNaN(ncdbl));
    SimTK_TEST(!isNaN(ncjflt) && !isNaN(ncjdbl));

    // Should be NaN if either or both parts are NaN.
    cflt = std::complex<float>(cflt.real(), fltNaN);
    cdbl = std::complex<double>(cdbl.real(), dblNaN);
    cjflt = conjugate<float>(cjflt.real(), fltNaN);
    cjdbl = conjugate<double>(cjdbl.real(), dblNaN);

    // Imaginary only is NaN.
    SimTK_TEST(isNaN(cflt) && isNaN(cdbl));
    SimTK_TEST(isNaN(cjflt) && isNaN(cjdbl));
    SimTK_TEST(isNaN(ncflt) && isNaN(ncdbl));
    SimTK_TEST(isNaN(ncjflt) && isNaN(ncjdbl));

    cflt = std::complex<float>(fltNaN, cflt.imag());
    cdbl = std::complex<double>(dblNaN, cdbl.imag());
    cjflt = conjugate<float>(fltNaN, cjflt.imag());
    cjdbl = conjugate<double>(dblNaN, cjdbl.imag());

    // Both parts are NaN.
    SimTK_TEST(isNaN(cflt) && isNaN(cdbl));
    SimTK_TEST(isNaN(cjflt) && isNaN(cjdbl));
    SimTK_TEST(isNaN(ncflt) && isNaN(ncdbl));
    SimTK_TEST(isNaN(ncjflt) && isNaN(ncjdbl));

    // Restore imaginary part to normal.
    cflt = std::complex<float>(cflt.real(), fltRegular);
    cdbl = std::complex<double>(cdbl.real(), dblRegular);
    cjflt = conjugate<float>(cjflt.real(), fltRegular);
    cjdbl = conjugate<double>(cjdbl.real(), dblRegular);

    // Real part only is NaN;
    SimTK_TEST(isNaN(cflt) && isNaN(cdbl));
    SimTK_TEST(isNaN(cjflt) && isNaN(cjdbl));
    SimTK_TEST(isNaN(ncflt) && isNaN(ncdbl));
    SimTK_TEST(isNaN(ncjflt) && isNaN(ncjdbl));
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

    SimTK_TEST(nfltInf == -fltInf);
    SimTK_TEST(ndblInf == -dblInf);

    SimTK_TEST(isInf(fltInf) && isInf(dblInf));
    SimTK_TEST(isInf(mfltInf) && isInf(mdblInf));
    SimTK_TEST(isInf(nfltInf) && isInf(ndblInf));
    SimTK_TEST(!isInf(fltRegular) && !isInf(dblRegular));

    std::complex<float> cflt(fltRegular, -2*fltRegular);
    std::complex<double> cdbl(dblRegular, -2*dblRegular);
    conjugate<float> cjflt(fltRegular, -2*fltRegular);
    conjugate<double> cjdbl(dblRegular, -2*dblRegular);

    SimTK_TEST(!isInf(cflt) && !isInf(cdbl));
    SimTK_TEST(!isInf(cjflt) && !isInf(cjdbl));

    // Reference the same memory as a negator of its contents.
    const negator<float>&           nflt   = reinterpret_cast<const negator<float>&>(fltRegular);
    const negator<double>&          ndbl   = reinterpret_cast<const negator<double>&>(dblRegular);
    negator<std::complex<float> >&   ncflt  = reinterpret_cast<negator<std::complex<float> >&> (cflt);
    negator<std::complex<double> >&  ncdbl  = reinterpret_cast<negator<std::complex<double> >&>(cdbl);
    negator<conjugate<float> >&      ncjflt = reinterpret_cast<negator<conjugate<float> >&>    (cjflt);
    negator<conjugate<double> >&     ncjdbl = reinterpret_cast<negator<conjugate<double> >&>   (cjdbl);

    // Test that negators are working properly.
    SimTK_TEST_EQ(nflt, -fltRegular);
    SimTK_TEST_EQ(ndbl, -dblRegular);
    SimTK_TEST_EQ(ncflt, -cflt);
    SimTK_TEST_EQ(-ncflt, cflt);
    SimTK_TEST_EQ(ncjflt, -cjflt);
    SimTK_TEST_EQ(-ncjflt, cjflt);

    SimTK_TEST(!isInf(nflt) && !isInf(ndbl));
    SimTK_TEST(!isInf(ncflt) && !isInf(ncdbl));
    SimTK_TEST(!isInf(ncjflt) && !isInf(ncjdbl));

    // Should be Inf if either or both parts are Inf, as long as neither
    // part is NaN.
    cflt = std::complex<float>(cflt.real(), fltInf);
    cdbl = std::complex<double>(cdbl.real(), dblInf);
    cjflt = conjugate<float>(cjflt.real(), fltInf);
    cjdbl = conjugate<double>(cjdbl.real(), dblInf);

    // Imaginary only is Inf.
    SimTK_TEST(isInf(cflt) && isInf(cdbl));
    SimTK_TEST(isInf(cjflt) && isInf(cjdbl));
    SimTK_TEST(isInf(ncflt) && isInf(ncdbl));
    SimTK_TEST(isInf(ncjflt) && isInf(ncjdbl));

    cflt = std::complex<float>(fltInf, cflt.imag());
    cdbl = std::complex<double>(dblInf, cdbl.imag());
    cjflt = conjugate<float>(fltInf, cjflt.imag());
    cjdbl = conjugate<double>(dblInf, cjdbl.imag());

    // Both parts are Inf.
    SimTK_TEST(isInf(cflt) && isInf(cdbl));
    SimTK_TEST(isInf(cjflt) && isInf(cjdbl));
    SimTK_TEST(isInf(ncflt) && isInf(ncdbl));
    SimTK_TEST(isInf(ncjflt) && isInf(ncjdbl));

    // Restore imaginary part to normal.
    cflt = std::complex<float>(cflt.real(), fltRegular);
    cdbl = std::complex<double>(cdbl.real(), dblRegular);
    cjflt = conjugate<float>(cjflt.real(), fltRegular);
    cjdbl = conjugate<double>(cjdbl.real(), dblRegular);

    // Real part only is Inf;
    SimTK_TEST(isInf(cflt) && isInf(cdbl));
    SimTK_TEST(isInf(cjflt) && isInf(cjdbl));
    SimTK_TEST(isInf(ncflt) && isInf(ncdbl));
    SimTK_TEST(isInf(ncjflt) && isInf(ncjdbl));

    // Set real part to minus infinity.
    cflt = std::complex<float>(mfltInf, cflt.imag());
    cdbl = std::complex<double>(mdblInf, cdbl.imag());
    cjflt = conjugate<float>(mfltInf, cjflt.imag());
    cjdbl = conjugate<double>(mdblInf, cjdbl.imag());

    SimTK_TEST(isInf(cflt) && isInf(cdbl));
    SimTK_TEST(isInf(cjflt) && isInf(cjdbl));
    SimTK_TEST(isInf(ncflt) && isInf(ncdbl));
    SimTK_TEST(isInf(ncjflt) && isInf(ncjdbl));

    // Set real part to NaN.
    const float fltNaN = NTraits<float>::getNaN();
    const double dblNaN = NTraits<double>::getNaN();
    cflt = std::complex<float>(fltNaN, cflt.imag());
    cdbl = std::complex<double>(dblNaN, cdbl.imag());
    cjflt = conjugate<float>(fltNaN, cjflt.imag());
    cjdbl = conjugate<double>(dblNaN, cjdbl.imag());

    SimTK_TEST(!isInf(cflt) && !isInf(cdbl));
    SimTK_TEST(!isInf(cjflt) && !isInf(cjdbl));
    SimTK_TEST(!isInf(ncflt) && !isInf(ncdbl));
    SimTK_TEST(!isInf(ncjflt) && !isInf(ncjdbl));
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

    SimTK_TEST(isFinite(fltRegular) && isFinite(dblRegular));
    SimTK_TEST(!isFinite(fltNaN) && !isFinite(dblNaN));
    SimTK_TEST(!isFinite(fltInf) && !isFinite(dblInf));
    SimTK_TEST(!isFinite(mfltInf) && !isFinite(mdblInf));

    std::complex<float> cflt(fltRegular, -2*fltRegular);
    std::complex<double> cdbl(dblRegular, -2*dblRegular);
    conjugate<float> cjflt(fltRegular, -2*fltRegular);
    conjugate<double> cjdbl(dblRegular, -2*dblRegular);

    SimTK_TEST(isFinite(cflt) && isFinite(cdbl));
    SimTK_TEST(isFinite(cjflt) && isFinite(cjdbl));

    // Reference the same memory as a negator of its contents.
    const negator<float>&           nflt   = reinterpret_cast<const negator<float>&>(fltRegular);
    const negator<double>&          ndbl   = reinterpret_cast<const negator<double>&>(dblRegular);
    negator<std::complex<float> >&   ncflt  = reinterpret_cast<negator<std::complex<float> >&> (cflt);
    negator<std::complex<double> >&  ncdbl  = reinterpret_cast<negator<std::complex<double> >&>(cdbl);
    negator<conjugate<float> >&      ncjflt = reinterpret_cast<negator<conjugate<float> >&>    (cjflt);
    negator<conjugate<double> >&     ncjdbl = reinterpret_cast<negator<conjugate<double> >&>   (cjdbl);

    // Test that negators are working properly.
    SimTK_TEST_EQ(nflt, -fltRegular);
    SimTK_TEST_EQ(ndbl, -dblRegular);
    SimTK_TEST_EQ(ncflt, -cflt);
    SimTK_TEST_EQ(-ncflt, cflt);
    SimTK_TEST_EQ(ncjflt, -cjflt);
    SimTK_TEST_EQ(-ncjflt, cjflt);

    SimTK_TEST(isFinite(nflt) && isFinite(ndbl));
    SimTK_TEST(isFinite(ncflt) && isFinite(ncdbl));
    SimTK_TEST(isFinite(ncjflt) && isFinite(ncjdbl));

    // Should be finite only if both parts are finite.
    cflt = std::complex<float>(cflt.real(),  fltInf);
    cdbl = std::complex<double>(cdbl.real(), mdblInf);
    cjflt = conjugate<float>(cjflt.real(),   fltNaN);
    cjdbl = conjugate<double>(cjdbl.real(),  dblInf);

    // Imaginary only is NaN.
    SimTK_TEST(!isFinite(cflt) && !isFinite(cdbl));
    SimTK_TEST(!isFinite(cjflt) && !isFinite(cjdbl));
    SimTK_TEST(!isFinite(ncflt) && !isFinite(ncdbl));
    SimTK_TEST(!isFinite(ncjflt) && !isFinite(ncjdbl));

    cflt = std::complex<float> (fltInf, cflt.imag());
    cdbl = std::complex<double>(mdblInf, cdbl.imag());
    cjflt = conjugate<float>   (fltNaN, cjflt.imag());
    cjdbl = conjugate<double>  (dblInf, cjdbl.imag());

    // Both parts are non-finite.
    SimTK_TEST(!isFinite(cflt) && !isFinite(cdbl));
    SimTK_TEST(!isFinite(cjflt) && !isFinite(cjdbl));
    SimTK_TEST(!isFinite(ncflt) && !isFinite(ncdbl));
    SimTK_TEST(!isFinite(ncjflt) && !isFinite(ncjdbl));

    // Restore imaginary part to normal.
    cflt = std::complex<float>(cflt.real(), fltRegular);
    cdbl = std::complex<double>(cdbl.real(), dblRegular);
    cjflt = conjugate<float>(cjflt.real(), fltRegular);
    cjdbl = conjugate<double>(cjdbl.real(), dblRegular);

    // Real part only is non-finite;
    SimTK_TEST(!isFinite(cflt) && !isFinite(cdbl));
    SimTK_TEST(!isFinite(cjflt) && !isFinite(cjdbl));
    SimTK_TEST(!isFinite(ncflt) && !isFinite(ncdbl));
    SimTK_TEST(!isFinite(ncjflt) && !isFinite(ncjdbl));
}

void testSignBit() {
    const unsigned char ucm=0xff, ucz=0, ucp=27;
    const unsigned short usm=0xffff, usz=0, usp=2342;
    const unsigned int   uim=0xffffffff, uiz=0, uip=2342344;
    const unsigned long ulm=(unsigned long)-23423L, ulz=0, ulp=234234UL;
    const unsigned long long ullm=(unsigned long long)-234234234LL, ullz=0, ullp=234234234ULL;

    SimTK_TEST(!(signBit(ucm)||signBit(ucz)||signBit(ucp)));
    SimTK_TEST(!(signBit(usm)||signBit(usz)||signBit(usp)));
    SimTK_TEST(!(signBit(uim)||signBit(uiz)||signBit(uip)));
    SimTK_TEST(!(signBit(ulm)||signBit(ulz)||signBit(ulp)));
    SimTK_TEST(!(signBit(ullm)||signBit(ullz)||signBit(ullp)));
    
    // Note that signBit(char) doesn't exist.

    const signed char cm=-23, cz=0, cp=99;
    const short sm=-1234, sz=0, sp=23423;
    const int im=-2342343, iz=0, ip=29472383;
    const long lm=-43488, lz=0, lp=3454545;
    const long long llm=-2342342343433LL, llz=0, llp=874578478478574LL;

    SimTK_TEST(signBit(cm) && !(signBit(cz)||signBit(cp)));
    SimTK_TEST(signBit(sm) && !(signBit(sz)||signBit(sp)));
    SimTK_TEST(signBit(im) && !(signBit(iz)||signBit(ip)));
    SimTK_TEST(signBit(lm) && !(signBit(lz)||signBit(lp)));
    SimTK_TEST(signBit(llm) && !(signBit(llz)||signBit(llp)));

    const float fm=-12398.34f, fz=0, fp=4354.331f;
    const double dm=-234234.454, dz=0, dp=345345.2342;
    float mfz=-fz; double mdz=-dz;// -0 for some compilers

    SimTK_TEST(signBit(fm) && !(signBit(fz)||signBit(fp)));
    SimTK_TEST(signBit(dm) && !(signBit(dz)||signBit(dp)));

    // Can't be sure whether the compiler will actually have produced
    // a minus zero here.
    // SimTK_TEST(signBit(mfz) && signBit(mdz));

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

    SimTK_TEST(signBit(nfm) && !(signBit(nfz)||signBit(nfp)));
    SimTK_TEST(signBit(ndm) && !(signBit(ndz)||signBit(ndp)));
    SimTK_TEST(signBit(nmfz)==signBit(mfz) && signBit(nmdz)==signBit(mdz));

    const float fltInf = NTraits<float>::getInfinity();
    const double dblInf = NTraits<double>::getInfinity();
    const float mfltInf = -fltInf;
    const double mdblInf = -dblInf;

    SimTK_TEST(!signBit(fltInf) && !signBit(dblInf));
    SimTK_TEST(signBit(mfltInf) && signBit(mdblInf));
}


void testSign() {
    const unsigned char ucm=0xff, ucz=0, ucp=27;
    const unsigned short usm=0xffff, usz=0, usp=2342;
    const unsigned int   uim=0xffffffff, uiz=0, uip=2342344;
    const unsigned long ulm=(unsigned long)-23423L, ulz=0, ulp=234234UL;
    const unsigned long long ullm=(unsigned long long)-234234234LL, ullz=0, ullp=234234234ULL;

    SimTK_TEST(sign(ucm)==1 && sign(ucz)==0 && sign(ucp)==1);
    SimTK_TEST(sign(usm)==1 && sign(usz)==0 && sign(usp)==1);
    SimTK_TEST(sign(uim)==1 && sign(uiz)==0 && sign(uip)==1);
    SimTK_TEST(sign(ulm)==1 && sign(ulz)==0 && sign(ulp)==1);
    SimTK_TEST(sign(ullm)==1 && sign(ullz)==0 && sign(ullp)==1);

    // Note that sign(char) doesn't exist.

    const signed char cm=-23, cz=0, cp=99;
    const short sm=-1234, sz=0, sp=23423;
    const int im=-2342343, iz=0, ip=29472383;
    const long lm=-43488, lz=0, lp=3454545;
    const long long llm=-2342342343433LL, llz=0, llp=874578478478574LL;

    SimTK_TEST(sign(cm)==-1 && sign(cz)==0 && sign(cp)==1);
    SimTK_TEST(sign(sm)==-1 && sign(sz)==0 && sign(sp)==1);
    SimTK_TEST(sign(im)==-1 && sign(iz)==0 && sign(ip)==1);
    SimTK_TEST(sign(lm)==-1 && sign(lz)==0 && sign(lp)==1);
    SimTK_TEST(sign(llm)==-1 && sign(llz)==0 && sign(llp)==1);

    const float fm=-12398.34f, fz=0, fp=4354.331f;
    const double dm=-234234.454, dz=0, dp=345345.2342;
    float mfz=-fz; double mdz=-dz;// -0

    SimTK_TEST(sign(fm)==-1 && sign(fz)==0 && sign(fp)==1);
    SimTK_TEST(sign(dm)==-1 && sign(dz)==0 && sign(dp)==1);
    SimTK_TEST(sign(mfz)==0 && sign(mdz)==0); // doesn't matter if it's -0

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

    SimTK_TEST(sign(nfm)==1 && sign(nfz)==0 && sign(nfp)==-1);
    SimTK_TEST(sign(ndm)==1 && sign(ndz)==0 && sign(ndp)==-1);
    SimTK_TEST(sign(nmfz)==0 && sign(nmdz)==0); // doesn't matter if it's -0

    const float fltInf = NTraits<float>::getInfinity();
    const double dblInf = NTraits<double>::getInfinity();
    const float mfltInf = -fltInf;
    const double mdblInf = -dblInf;
    const negator<float>& nfltInf = reinterpret_cast<const negator<float>&>(fltInf);
    const negator<double>& ndblInf = reinterpret_cast<const negator<double>&>(dblInf);

    SimTK_TEST(sign(fltInf)==1 && sign(dblInf)==1);
    SimTK_TEST(sign(mfltInf)==-1 && sign(mdblInf)==-1);
    SimTK_TEST(sign(nfltInf)==-1 && sign(ndblInf)==-1);
}

void testSquareAndCube() {
    const float fval = -23.33f;
    const double dval = -234443.441;
    const negator<float>& nfval = reinterpret_cast<const negator<float>&>(fval);
    const negator<double>& ndval = reinterpret_cast<const negator<double>&>(dval);

    // Basic test.
    SimTK_TEST_EQ(square(fval), fval*fval);
    SimTK_TEST_EQ(square(dval), dval*dval);
    SimTK_TEST_EQ(cube(fval), fval*fval*fval);
    SimTK_TEST_EQ(cube(dval), dval*dval*dval);

    // Test scalar negators.
    SimTK_TEST_EQ(square(nfval), nfval*nfval);
    SimTK_TEST_EQ(square(nfval), fval*fval);
    SimTK_TEST_EQ(square(ndval), ndval*ndval);
    SimTK_TEST_EQ(square(ndval), dval*dval);
    SimTK_TEST_EQ(cube(nfval), nfval*nfval*nfval);
    SimTK_TEST_EQ(cube(nfval), -fval*fval*fval);
    SimTK_TEST_EQ(cube(ndval), ndval*ndval*ndval);
    SimTK_TEST_EQ(cube(ndval), -dval*dval*dval);

    // Create complex and conjugate values.

    std::complex<float> fc(-234.343f, 45345e7f);
    std::complex<double> dc(-234.343, 45345e7);
    conjugate<float> fcj(-19.1e3f, -454.234f);
    conjugate<double> dcj(-19.1e3, -454.234);

    // Manual conjugates
    std::complex<float>  fcmj(fcj.real(), fcj.imag());
    std::complex<double> dcmj(dcj.real(), dcj.imag());
    SimTK_TEST(fcj == fcmj);    // sign change only; should be exact
    SimTK_TEST(dcj == dcmj);
    SimTK_TEST_EQ(fcj*fcj, fcmj*fcmj);
    SimTK_TEST_EQ(dcj*dcj, dcmj*dcmj);
    SimTK_TEST_EQ(fcj*fcj*fcj, fcmj*fcmj*fcmj);
    SimTK_TEST_EQ(dcj*dcj*dcj, dcmj*dcmj*dcmj);

    // Negators of complex an conjugate.
    negator<std::complex<float> >&   nfc  = reinterpret_cast<negator<std::complex<float> >&> (fc);
    negator<std::complex<double> >&  ndc  = reinterpret_cast<negator<std::complex<double> >&>(dc);
    negator<conjugate<float> >&      nfcj = reinterpret_cast<negator<conjugate<float> >&>    (fcj);
    negator<conjugate<double> >&     ndcj = reinterpret_cast<negator<conjugate<double> >&>   (dcj);

    // Change of sign should be exact.
    SimTK_TEST(nfc == -fc);
    SimTK_TEST(ndc == -dc);
    SimTK_TEST(nfcj == -fcj);
    SimTK_TEST(ndcj == -dcj);


    // Basic complex and conjugate tests.
    SimTK_TEST_EQ(square(fc), fc*fc);
    SimTK_TEST_EQ(cube(fc), fc*fc*fc);
    SimTK_TEST_EQ(square(dc), dc*dc);
    SimTK_TEST_EQ(cube(dc), dc*dc*dc);
    SimTK_TEST_EQ(square(fcj), fcj*fcj);
    SimTK_TEST_EQ(cube(fcj), fcj*fcj*fcj);
    SimTK_TEST_EQ(square(dcj), dcj*dcj);
    SimTK_TEST_EQ(cube(dcj), dcj*dcj*dcj);

    // Tests involving negators of complex and conjugate.
    SimTK_TEST_EQ(square(nfc), nfc*nfc); 
    SimTK_TEST_EQ(square(nfc), fc*fc);
    SimTK_TEST_EQ(square(ndc), ndc*ndc);
    SimTK_TEST_EQ(square(ndc), dc*dc);

    SimTK_TEST_EQ(cube(nfc), nfc*nfc*nfc); 
    SimTK_TEST_EQ(cube(nfc), -fc*fc*fc);
    SimTK_TEST_EQ(cube(ndc), ndc*ndc*ndc);
    SimTK_TEST_EQ(cube(ndc), -dc*dc*dc);

    SimTK_TEST_EQ(square(nfcj), nfcj*nfcj); 
    SimTK_TEST_EQ(square(nfcj), fcj*fcj);
    SimTK_TEST_EQ(square(ndcj), ndcj*ndcj);
    SimTK_TEST_EQ(square(ndcj), dcj*dcj);

    SimTK_TEST_EQ(cube(nfcj), nfcj*nfcj*nfcj); 
    SimTK_TEST_EQ(cube(nfcj), -fcj*fcj*fcj);
    SimTK_TEST_EQ(cube(ndcj), ndcj*ndcj*ndcj);
    SimTK_TEST_EQ(cube(ndcj), -dcj*dcj*dcj);
}

void testIsNumericallyEqual() {
    const float  f=1.234f, fn=1.234f+1e-5f, fe=1.234f+1e-9f;
    const double d=1.234,  dn=1.234 +1e-12, de=1.234 +1e-15;
    const negator<float>& nf=negator<float>::recast(f);
    const negator<float>& nfn=negator<float>::recast(fn);
    const negator<float>& nfe=negator<float>::recast(fe);

    SimTK_TEST(isNumericallyEqual(f,f))
    SimTK_TEST(isNumericallyEqual(f,fe));
    SimTK_TEST(!isNumericallyEqual(f,fn));
    SimTK_TEST(isNumericallyEqual(f,fn,1e-4f));
    SimTK_TEST(!isNumericallyEqual(f,fn,1e-6f));

    SimTK_TEST(CNT<float>::isNumericallyEqual(f,f));
    SimTK_TEST(CNT<float>::isNumericallyEqual(f,fe));
    SimTK_TEST(!CNT<float>::isNumericallyEqual(f,fn));
    SimTK_TEST(CNT<float>::isNumericallyEqual(f,fn,1e-4f));
    SimTK_TEST(!CNT<float>::isNumericallyEqual(f,fn,1e-6f));

    SimTK_TEST(nf.isNumericallyEqual(nf));
    SimTK_TEST(nf.isNumericallyEqual(-f));
    SimTK_TEST(!nf.isNumericallyEqual(f));

    SimTK_TEST(isNumericallyEqual(1000*f,1234));
    SimTK_TEST(isNumericallyEqual(1234,1000*f));
    SimTK_TEST(isNumericallyEqual(1000*fe,1234));
    SimTK_TEST(isNumericallyEqual(1234,1000*fe));
    SimTK_TEST(!isNumericallyEqual(1000*fn,1234));
    SimTK_TEST(!isNumericallyEqual(1234,1000*fn));

    SimTK_TEST(isNumericallyEqual(d,d));
    SimTK_TEST(isNumericallyEqual(d,de));
    SimTK_TEST(!isNumericallyEqual(d,dn));
    SimTK_TEST(isNumericallyEqual(1000*d,1234));
    SimTK_TEST(isNumericallyEqual(1234,1000*d));
    SimTK_TEST(isNumericallyEqual(1000*de,1234));
    SimTK_TEST(isNumericallyEqual(1234,1000*de));
    SimTK_TEST(!isNumericallyEqual(1000*dn,1234));
    SimTK_TEST(!isNumericallyEqual(1234,1000*dn));

    // Mixed should use float tolerance
    SimTK_TEST(isNumericallyEqual(fe,de));
    SimTK_TEST(!isNumericallyEqual((double)fe,de));
}

void testClamp() {
    const int i4=4;
    const double d325=3.25;
    const float fn325=-3.25;

    // int
    SimTK_TEST(clamp(4,i4,4)==4);
    SimTK_TEST(clamp(0,i4,9)==4);
    SimTK_TEST(clamp(5,i4,9)==5);
    SimTK_TEST(clamp(-7,i4,-5)==-5);

    // double
    SimTK_TEST(clamp(3.25,d325,3.25)==3.25);
    SimTK_TEST(clamp(0.,d325,9.)==3.25);
    SimTK_TEST(clamp(5.,d325,9.)==5);
    SimTK_TEST(clamp(-7.,d325,-5.)==-5);

    // float
    SimTK_TEST(clamp(-3.25f,fn325,-3.25f)==-3.25);
    SimTK_TEST(clamp(-9.f,fn325,0.f)==-3.25f);
    SimTK_TEST(clamp(-9.f,fn325,-5.f)==-5);
    SimTK_TEST(clamp(5.f,fn325,7.f)==5);

    // Test methods that take integer bounds.
    SimTK_TEST(clamp(0,d325,9)==3.25);
    SimTK_TEST(clamp(5,d325,9)==5);
    SimTK_TEST(clamp(-7,d325,-5)==-5);

    SimTK_TEST(clamp(-9,fn325,0)==-3.25);
    SimTK_TEST(clamp(-9,fn325,-5)==-5);
    SimTK_TEST(clamp(5,fn325,7)==5);

    SimTK_TEST(clamp(0.,d325,9)==3.25);
    SimTK_TEST(clamp(5.,d325,9)==5);
    SimTK_TEST(clamp(-7.,d325,-5)==-5);

    SimTK_TEST(clamp(-9.f,fn325,0)==-3.25);
    SimTK_TEST(clamp(-9.f,fn325,-5)==-5);
    SimTK_TEST(clamp(5.f,fn325,7)==5);

    SimTK_TEST(clamp(0,d325,9.)==3.25);
    SimTK_TEST(clamp(5,d325,9.)==5);
    SimTK_TEST(clamp(-7,d325,-5.)==-5);

    SimTK_TEST(clamp(-9,fn325,0.f)==-3.25);
    SimTK_TEST(clamp(-9,fn325,-5.f)==-5);
    SimTK_TEST(clamp(5,fn325,7.f)==5);

    int i; double d; float f;
    i=i4; 
    SimTK_TEST(clampInPlace(-2,i,3)==3 && i==3);
    d=d325;
    SimTK_TEST(clampInPlace(-2.,d,3.)==3 && d==3);
    f=fn325;
    SimTK_TEST(clampInPlace(-2,f,3)==-2 && f==-2);

    // Do a test for each of the less-common supported types.
    char c='j'; unsigned char uc=3; signed char sc=-2;
    SimTK_TEST(clamp('a',c,'e')=='e');
    SimTK_TEST(clamp('a',c,'z')=='j');
    SimTK_TEST(clamp((unsigned char)4,uc,(unsigned char)5)==4);
    SimTK_TEST(clamp((signed char)-7,sc,(signed char)-1)==-2);

    short s=-32000; unsigned short us=17; unsigned ui=4023456789U;
    SimTK_TEST(clamp((short)-29000,s,(short)400)==-29000);
    SimTK_TEST(clamp((unsigned short)4,us,(unsigned short)15)==15);
    SimTK_TEST(clamp(100000000U,ui,4010000000U)==4010000000U);

    long l=-234234L; unsigned long ul=293493849UL; 
    long long ll=-123456789123LL; unsigned long long ull=123456789123ULL;
    SimTK_TEST(clamp(-1000000L,l,-200000L)==-234234);
    SimTK_TEST(clamp(1000000UL,ul,4000000000UL)==293493849);
    SimTK_TEST(clamp(-100000000000LL,ll,27LL)==-100000000000LL);
    SimTK_TEST(clamp(-1000000000000LL,ll,27LL)==-123456789123LL);
}

void testStep() {
        // double
    SimTK_TEST(stepUp(0.)==0 && stepUp(.5)==.5 && stepUp(1.)==1);
    SimTK_TEST(0 < stepUp(.3) && stepUp(.3) < .5);
    SimTK_TEST(.5 < stepUp(.7) && stepUp(.7) < 1);
    SimTK_TEST(stepDown(0.)==1 && stepDown(.5)==.5 && stepDown(1.)==0);
    SimTK_TEST(.5 < stepDown(.3) && stepDown(.3) < 1);
    SimTK_TEST(0 < stepDown(.7) && stepDown(.7) < .5);
    SimTK_TEST(dstepUp(0.)==0 && dstepUp(.5)>0 && dstepUp(1.)==0);
    SimTK_TEST(dstepDown(0.)==0 && dstepDown(.5)<0 && dstepDown(1.)==0);
    SimTK_TEST(d2stepUp(0.)==0 && d2stepUp(1.)==0);
    SimTK_TEST(d2stepDown(0.)==0 && d2stepDown(1.)==0);
        // float
    SimTK_TEST(stepUp(0.f)==0 && stepUp(.5f)==.5f && stepUp(1.f)==1);
    SimTK_TEST(0 < stepUp(.3f) && stepUp(.3f) < .5f);
    SimTK_TEST(.5f < stepUp(.7f) && stepUp(.7f) < 1);
    SimTK_TEST(stepDown(0.f)==1 && stepDown(.5f)==.5f && stepDown(1.f)==0);
    SimTK_TEST(.5f < stepDown(.3f) && stepDown(.3f) < 1);
    SimTK_TEST(0 < stepDown(.7f) && stepDown(.7f) < .5f);
    SimTK_TEST(dstepUp(0.f)==0 && dstepUp(.5f)>0 && dstepUp(1.f)==0);
    SimTK_TEST(dstepDown(0.f)==0 && dstepDown(.5f)<0 && dstepDown(1.f)==0);
    SimTK_TEST(d2stepUp(0.f)==0 && d2stepUp(1.f)==0);
    SimTK_TEST(d2stepDown(0.f)==0 && d2stepDown(1.f)==0);
        // long double
    SimTK_TEST(stepUp(0.L)==0 && stepUp(.5L)==.5L && stepUp(1.L)==1);
    SimTK_TEST(0 < stepUp(.3L) && stepUp(.3L) < .5L);
    SimTK_TEST(.5L < stepUp(.7L) && stepUp(.7L) < 1);
    SimTK_TEST(stepDown(0.L)==1 && stepDown(.5L)==.5L && stepDown(1.L)==0);
    SimTK_TEST(.5L < stepDown(.3L) && stepDown(.3L) < 1);
    SimTK_TEST(0 < stepDown(.7L) && stepDown(.7L) < .5L);
    SimTK_TEST(dstepUp(0.L)==0 && dstepUp(.5L)>0 && dstepUp(1.L)==0);
    SimTK_TEST(dstepDown(0.L)==0 && dstepDown(.5L)<0 && dstepDown(1.L)==0);
    SimTK_TEST(d2stepUp(0.L)==0 && d2stepUp(1.L)==0);
    SimTK_TEST(d2stepDown(0.L)==0 && d2stepDown(1.L)==0);

        // int is treated as a double, but only for stepUp()/stepDown()
    SimTK_TEST(stepUp(0)==0 && stepUp(1)==1);
    SimTK_TEST(stepDown(0)==1 && stepDown(1)==0);

    // Don't know anything analytic about d3 but can test with finite
    // differencing below.

    // Central difference estimates should give around 10 
    // decimal places in double, 4 in float.
    const double dupEst = (stepUp(.799+1e-6)-stepUp(.799-1e-6))/2e-6;
    const double ddnEst = (stepDown(.799+1e-6)-stepDown(.799-1e-6))/2e-6;
    const double d2upEst = (dstepUp(.723+1e-6)-dstepUp(.723-1e-6))/2e-6;
    const double d2dnEst = (dstepDown(.723+1e-6)-dstepDown(.723-1e-6))/2e-6;
    const double d3upEst = (d2stepUp(.123+1e-6)-d2stepUp(.123-1e-6))/2e-6;
    const double d3dnEst = (d2stepDown(.123+1e-6)-d2stepDown(.123-1e-6))/2e-6;
    SimTK_TEST_EQ_TOL(dstepUp(.799), dupEst, 1e-8);
    SimTK_TEST_EQ_TOL(dstepDown(.799), ddnEst, 1e-8);
    SimTK_TEST_EQ_TOL(d2stepUp(.723), d2upEst, 1e-8);
    SimTK_TEST_EQ_TOL(d2stepDown(.723), d2dnEst, 1e-8);
    SimTK_TEST_EQ_TOL(d3stepUp(.123), d3upEst, 1e-8);
    SimTK_TEST_EQ_TOL(d3stepDown(.123), d3dnEst, 1e-8);

    const float fdupEst = (stepUp(.699f+1e-3f)-stepUp(.699f-1e-3f))/2e-3f;
    const float fddnEst = (stepDown(.699f+1e-3f)-stepDown(.699f-1e-3f))/2e-3f;
    const float fd2upEst = (dstepUp(.623f+1e-3f)-dstepUp(.623f-1e-3f))/2e-3f;
    const float fd2dnEst = (dstepDown(.623f+1e-3f)-dstepDown(.623f-1e-3f))/2e-3f;
    const float fd3upEst = (d2stepUp(.211f+1e-3f)-d2stepUp(.211f-1e-3f))/2e-3f;
    const float fd3dnEst = (d2stepDown(.211f+1e-3f)-d2stepDown(.211f-1e-3f))/2e-3f;
    SimTK_TEST_EQ_TOL(dstepUp(.699f), fdupEst, 1e-3);
    SimTK_TEST_EQ_TOL(dstepDown(.699f), fddnEst, 1e-3);
    SimTK_TEST_EQ_TOL(d2stepUp(.623f), fd2upEst, 1e-3);
    SimTK_TEST_EQ_TOL(d2stepDown(.623f), fd2dnEst, 1e-3);
    SimTK_TEST_EQ_TOL(d3stepUp(.211f), fd3upEst, 1e-3);
    SimTK_TEST_EQ_TOL(d3stepDown(.211f), fd3dnEst, 1e-3);


    // y = stepAny(y0,yrange,x0,1/xrange, x)
    // y goes from -1 to 1 as x goes from 0 to 1, exact arithmetic.
    SimTK_TEST(stepAny(-1,2,0,1,0.) == -1);
    SimTK_TEST(stepAny(-1,2,0,1,.5) == 0);
    SimTK_TEST(stepAny(-1,2,0,1,1.) == 1);
    SimTK_TEST(stepAny(-1,2,0,1,0.f) == -1);
    SimTK_TEST(stepAny(-1,2,0,1,.5f) == 0);
    SimTK_TEST(stepAny(-1,2,0,1,1.f) == 1);
    SimTK_TEST(stepAny(-1,2,0,1,0.L) == -1);
    SimTK_TEST(stepAny(-1,2,0,1,.5L) == 0);
    SimTK_TEST(stepAny(-1,2,0,1,1.L) == 1);

    // y goes from -7 down to -14 as x goes from -3.1 up to +429.3.
    const double x0=-3.1, x1=429.3, y0=-7., y1=-14.;
    const double xr=(x1-x0), ooxr=1/xr, yr=(y1-y0);
    SimTK_TEST_EQ(stepAny(y0,yr,x0,ooxr,-3.1), y0);
    SimTK_TEST_EQ(stepAny(y0,yr,x0,ooxr,429.3), y1);
    SimTK_TEST_EQ(stepAny(y0,yr,x0,ooxr,x0+xr/2), y0+yr/2);

    const float fx0=-3.1f, fx1=429.3f, fy0=-7.f, fy1=-14.f;
    const float fxr=(fx1-fx0), fooxr=1/fxr, fyr=(fy1-fy0);
    SimTK_TEST_EQ(stepAny(fy0,fyr,fx0,fooxr,-3.1f),fy0);
    SimTK_TEST_EQ(stepAny(fy0,fyr,fx0,fooxr,429.3f),fy1);
    SimTK_TEST_EQ(stepAny(fy0,fyr,fx0,fooxr,fx0+fxr/2),fy0+fyr/2);

    // Check derivatives
    const double danyEst = 
        (stepAny(y0,yr,x0,ooxr,.799+1e-6)-stepAny(y0,yr,x0,ooxr,.799-1e-6))/2e-6;
    const double d2anyEst = 
        (dstepAny(yr,x0,ooxr,.723+1e-6)-dstepAny(yr,x0,ooxr,.723-1e-6))/2e-6;
    const double d3anyEst = 
        (d2stepAny(yr,x0,ooxr,.123+1e-6)-d2stepAny(yr,x0,ooxr,.123-1e-6))/2e-6;
    SimTK_TEST_EQ_TOL(dstepAny(yr,x0,ooxr,.799), danyEst, 1e-8);
    SimTK_TEST_EQ_TOL(d2stepAny(yr,x0,ooxr,.723), d2anyEst, 1e-8);
    SimTK_TEST_EQ_TOL(d3stepAny(yr,x0,ooxr,.123), d3anyEst, 1e-8);

    const float fdanyEst = 
        (stepAny(fy0,fyr,fx0,fooxr,.799f+1e-3f)
        -stepAny(fy0,fyr,fx0,fooxr,.799f-1e-3f))/2e-3f;
    const float fd2anyEst = 
        (dstepAny(fyr,fx0,fooxr,.723f+1e-3f)
        -dstepAny(fyr,fx0,fooxr,.723f-1e-3f))/2e-3f;
    const float fd3anyEst = 
        (d2stepAny(fyr,fx0,fooxr,.123f+1e-3f)
        -d2stepAny(fyr,fx0,fooxr,.123f-1e-3f))/2e-3f;
    SimTK_TEST_EQ_TOL(dstepAny(fyr,fx0,fooxr,.799f), fdanyEst, 1e-3);
    SimTK_TEST_EQ_TOL(d2stepAny(fyr,fx0,fooxr,.723f), fd2anyEst, 1e-3);
    SimTK_TEST_EQ_TOL(d3stepAny(fyr,fx0,fooxr,.123f), fd3anyEst, 1e-3);
}

void testDeadband() {
    // float
    SimTK_TEST(deadband(-3.f, 0.5f)==-3.f);
    SimTK_TEST(deadband( 3.f, 0.5f)== 3.f);
    SimTK_TEST(deadband( 0.3f, 0.5f)== 0.f);
    SimTK_TEST(deadband(-0.3f, 0.5f)== 0.f);
    SimTK_TEST(deadband( 0.5f, 0.5f)== 0.f);
    SimTK_TEST(deadband(-0.5f, 0.5f)== 0.f);
    SimTK_TEST(deadband( 0.5000005f, 0.5f)== 0.5000005f);
    SimTK_TEST(deadband(-0.5000005f, 0.5f)==-0.5000005f);
    SimTK_TEST(isNaN(deadband(fNaN, 1e-3f)));

    SimTK_TEST(deadband(-0.25f, -3.f, 0.5f)==-3.f);
    SimTK_TEST(deadband(-0.25f,  3.f, 0.5f)== 3.f);
    SimTK_TEST(deadband(-0.25f,  0.3f, 0.5f)== 0.f);
    SimTK_TEST(deadband(-0.25f, -0.275f, 0.5f)== -0.275f);
    SimTK_TEST(deadband(-0.25f, -0.2f, 0.5f)== 0.f);
    SimTK_TEST(deadband(-0.25f,  0.5f, 0.5f)== 0.f);
    SimTK_TEST(deadband(-0.25f, -0.25f, 0.5f)== 0.f);
    SimTK_TEST(deadband(-0.25f,  0.5000005f, 0.5f)== 0.5000005f);
    SimTK_TEST(deadband(-0.25f, -0.2500005f, 0.5f)==-0.2500005f);
    SimTK_TEST(isNaN(deadband(-1e-3f, fNaN, 1e-3f)));

    // double
    SimTK_TEST(deadband(-3., 0.5)==-3.);
    SimTK_TEST(deadband( 3., 0.5)== 3.);
    SimTK_TEST(deadband( 0.3, 0.5)== 0.);
    SimTK_TEST(deadband(-0.3, 0.5)== 0.);
    SimTK_TEST(deadband( 0.5, 0.5)== 0.);
    SimTK_TEST(deadband(-0.5, 0.5)== 0.);
    SimTK_TEST(deadband( 0.5000000000005, 0.5)== 0.5000000000005);
    SimTK_TEST(deadband(-0.5000000000005, 0.5)==-0.5000000000005);
    SimTK_TEST(isNaN(deadband(NaN, 1e-3)));

    SimTK_TEST(deadband(-0.25, -3., 0.5)==-3.);
    SimTK_TEST(deadband(-0.25,  3., 0.5)== 3.);
    SimTK_TEST(deadband(-0.25,  0.3, 0.5)== 0.);
    SimTK_TEST(deadband(-0.25, -0.275, 0.5)== -0.275);
    SimTK_TEST(deadband(-0.25, -0.2, 0.5)== 0.);
    SimTK_TEST(deadband(-0.25,  0.5, 0.5)== 0.);
    SimTK_TEST(deadband(-0.25, -0.25, 0.5)== 0.);
    SimTK_TEST(deadband(-0.25,  0.5000000000005, 0.5)== 0.5000000000005);
    SimTK_TEST(deadband(-0.25, -0.2500000000005, 0.5)==-0.2500000000005);
    SimTK_TEST(isNaN(deadband(-1e-3, NaN, 1e-3)));

    // long double
    const long double lNaN = NTraits<long double>::getNaN();
    SimTK_TEST(deadband(-3.L, 0.5L)==-3.L);
    SimTK_TEST(deadband( 3.L, 0.5L)== 3.L);
    SimTK_TEST(deadband( 0.3L, 0.5L)== 0.L);
    SimTK_TEST(deadband(-0.3L, 0.5L)== 0.L);
    SimTK_TEST(deadband( 0.5L, 0.5L)== 0.L);
    SimTK_TEST(deadband(-0.5L, 0.5L)== 0.L);
    SimTK_TEST(deadband( 0.500000000005L, 0.5L)
                      == 0.500000000005L);
    SimTK_TEST(deadband(-0.500000000005L, 0.5L)
                      ==-0.500000000005L);
    SimTK_TEST(isNaN(deadband(lNaN, 1e-3L)));

    SimTK_TEST(deadband(-0.25L, -3.L, 0.5L)==-3.L);
    SimTK_TEST(deadband(-0.25L,  3.L, 0.5L)== 3.L);
    SimTK_TEST(deadband(-0.25L,  0.3L, 0.5L)== 0.L);
    SimTK_TEST(deadband(-0.25L, -0.275L, 0.5L)== -0.275L);
    SimTK_TEST(deadband(-0.25L, -0.2L, 0.5L)== 0.L);
    SimTK_TEST(deadband(-0.25L,  0.5L, 0.5L)== 0.L);
    SimTK_TEST(deadband(-0.25L, -0.25L, 0.5L)== 0.L);
    SimTK_TEST(deadband(-0.25L,  0.500000000005L, 0.5L)
                              == 0.500000000005L);
    SimTK_TEST(deadband(-0.25L, -0.250000000005L, 0.5L)
                              ==-0.250000000005L);
    SimTK_TEST(isNaN(deadband(-1e-3L, lNaN, 1e-3L)));

}

void testBitScan() {
    SimTK_TEST(lowBitIndex((unsigned char)0x0)==-1);
    SimTK_TEST(lowBitIndex((unsigned char)0x1)==0);
    SimTK_TEST(lowBitIndex((unsigned char)0xc0)==6);
    SimTK_TEST(lowBitIndex((unsigned char)0x80)==7);
    SimTK_TEST(lowBitIndex((unsigned char)0xff)==0);

    SimTK_TEST(lowBitIndex((unsigned short)0x0)==-1);
    SimTK_TEST(lowBitIndex((unsigned short)0x1)==0);
    SimTK_TEST(lowBitIndex((unsigned short)0x2340u)==6);
    SimTK_TEST(lowBitIndex((unsigned short)0x8000)==15);
    SimTK_TEST(lowBitIndex((unsigned short)0xffff)==0);

    SimTK_TEST(lowBitIndex(0x0u)==-1);
    SimTK_TEST(lowBitIndex(0x1u)==0);
    SimTK_TEST(lowBitIndex(0x12340u)==6);
    SimTK_TEST(lowBitIndex(0x80000000u)==31);
    SimTK_TEST(lowBitIndex(0xffffffffu)==0);

    // Don't assume long is larger than int.
    SimTK_TEST(lowBitIndex(0x0ul)==-1);
    SimTK_TEST(lowBitIndex(0x1ul)==0);
    SimTK_TEST(lowBitIndex(0x12340ul)==6);
    SimTK_TEST(lowBitIndex(0x80000000ul)==31);
    SimTK_TEST(lowBitIndex(0xfffffffful)==0);

    SimTK_TEST(lowBitIndex(0ull)==-1);
    SimTK_TEST(lowBitIndex(1ull)==0);
    SimTK_TEST(lowBitIndex(0x12340ull)==6);
    SimTK_TEST(lowBitIndex(0x101800000000ull)==35);
    SimTK_TEST(lowBitIndex(0x8000000000000000ull)==63);
    SimTK_TEST(lowBitIndex(0xffffffffffffffffull)==0);

    SimTK_TEST(clearLowBit((unsigned char)0x0)==(unsigned char)0);
    SimTK_TEST(clearLowBit((unsigned char)0x20)==(unsigned char)0);
    SimTK_TEST(clearLowBit((unsigned char)0x22)==(unsigned char)0x20);

    SimTK_TEST(clearLowBit((unsigned short)0x0)==(unsigned short)0);
    SimTK_TEST(clearLowBit((unsigned short)0x2000)==(unsigned short)0);
    SimTK_TEST(clearLowBit((unsigned short)0x2002)==(unsigned short)0x2000);

    SimTK_TEST(clearLowBit(0x0u)==0u);
    SimTK_TEST(clearLowBit(0x20000000u)==0u);
    SimTK_TEST(clearLowBit(0x20000002u)==0x20000000u);

    SimTK_TEST(clearLowBit(0x0ul)==0ul);
    SimTK_TEST(clearLowBit(0x20000000ul)==0ul);
    SimTK_TEST(clearLowBit(0x20000002ul)==0x20000000ul);

    SimTK_TEST(clearLowBit(0x0ull)==0ull);
    SimTK_TEST(clearLowBit(0x2000000000000000ull)==0ull);
    SimTK_TEST(clearLowBit(0x2000000000000002ull)==0x2000000000000000ull);

    SimTK_TEST(isolateLowBit((unsigned char)0x0)==(unsigned char)0);
    SimTK_TEST(isolateLowBit((unsigned char)0x20)==(unsigned char)0x20);
    SimTK_TEST(isolateLowBit((unsigned char)0x22)==(unsigned char)0x02);

    SimTK_TEST(isolateLowBit((unsigned short)0x0)==(unsigned short)0);
    SimTK_TEST(isolateLowBit((unsigned short)0x2000)==(unsigned short)0x2000);
    SimTK_TEST(isolateLowBit((unsigned short)0x2002)==(unsigned short)0x0002);

    SimTK_TEST(isolateLowBit(0x0u)==0u);
    SimTK_TEST(isolateLowBit(0x20000000u)==0x20000000u);
    SimTK_TEST(isolateLowBit(0x20000002u)==0x00000002u);

    SimTK_TEST(isolateLowBit(0x0ul)==0ul);
    SimTK_TEST(isolateLowBit(0x20000000ul)==0x20000000ul);
    SimTK_TEST(isolateLowBit(0x20000002ul)==0x00000002ul);

    SimTK_TEST(isolateLowBit(0x0ull)==0ull);
    SimTK_TEST(isolateLowBit(0x2000000000000000ull)==0x2000000000000000ull);
    SimTK_TEST(isolateLowBit(0x2000000000000002ull)==0x0000000000000002ull);

    SimTK_TEST(countSetBits((unsigned char)0x0) == 0u);
    SimTK_TEST(countSetBits((unsigned short)0x0) == 0u);
    SimTK_TEST(countSetBits((unsigned)0x0) == 0u);
    SimTK_TEST(countSetBits((unsigned long)0x0) == 0u);
    SimTK_TEST(countSetBits((unsigned long long)0x0) == 0u);

    SimTK_TEST(countSetBits((unsigned char)0x1) == 1u);
    SimTK_TEST(countSetBits((unsigned short)0x1) == 1u);
    SimTK_TEST(countSetBits((unsigned)0x1) == 1u);
    SimTK_TEST(countSetBits((unsigned long)0x1) == 1u);
    SimTK_TEST(countSetBits((unsigned long long)0x1) == 1u);

    SimTK_TEST(countSetBits((unsigned char)0xc7) == 5u);
    SimTK_TEST(countSetBits((unsigned short)0xc0f7) == 9u);
    SimTK_TEST(countSetBits((unsigned)0xc0f00007) == 9u);
    SimTK_TEST(countSetBits((unsigned long)0xc000f007) == 9u);
    SimTK_TEST(countSetBits((unsigned long long)0xc0000f0000000007) == 9u);

    SimTK_TEST(isBitSet((char)0xc7, 2));
    SimTK_TEST(!isBitSet((char)0xc7, 3));
    SimTK_TEST(isBitSet((unsigned char)0xc7, 2));
    SimTK_TEST(!isBitSet((unsigned char)0xc7, 3));
    SimTK_TEST(isBitSet((short)0xf0c7, 15)&& isBitSet((short)0xf0c7, 0));
    SimTK_TEST(!isBitSet((short)0xf0c7, 11));
    SimTK_TEST(   isBitSet((unsigned short)0xf0c7, 15)
               && isBitSet((unsigned short)0xf0c7, 0));
    SimTK_TEST(!isBitSet((unsigned short)0xf0c7, 11));
    SimTK_TEST(isBitSet((int)0xf0c7f0f0, 31)&& isBitSet((int)0xf0f0f0c7, 0));
    SimTK_TEST(!isBitSet((int)0xf0c7f0f0, 27));
    SimTK_TEST(   isBitSet((unsigned)0xf0c7f0f0, 31)
               && isBitSet((unsigned)0xf0f0f0c7, 0));
    SimTK_TEST(!isBitSet((unsigned)0xf0c7f0f0, 27));
    if (sizeof(long)==sizeof(int)) { // MSVC
        SimTK_TEST(isBitSet((long)0xf0c7f0f0, 31)&& isBitSet((long)0xf0f0f0c7, 0));
        SimTK_TEST(!isBitSet((long)0xf0c7f0f0, 27));
        SimTK_TEST(   isBitSet((unsigned long)0xf0c7f0f0, 31)
                   && isBitSet((unsigned long)0xf0f0f0c7, 0));
        SimTK_TEST(!isBitSet((unsigned long)0xf0c7f0f0, 27));
    }

    SimTK_TEST(   isBitSet(0xf0f0f0f0f0c7f0f0LL, 63)
               && isBitSet(0xf0f0f0f0f0f0f0c7LL, 0));
    SimTK_TEST(!isBitSet(0xf0f0f0f0f0c7f0f0LL, 59));
    SimTK_TEST(   isBitSet(0xf0f0f0f0f0c7f0f0ULL, 63)
               && isBitSet(0xf0f0f0f0f0f0f0c7ULL, 0));
    SimTK_TEST(!isBitSet(0xf0f0f0f0f0c7f0f0ULL, 59));
    if (sizeof(long)==sizeof(long long)) { // gcc, clang
        SimTK_TEST(   isBitSet(0xf0f0f0f0f0c7f0f0L, 63)
                   && isBitSet(0xf0f0f0f0f0f0f0c7L, 0));
        SimTK_TEST(!isBitSet(0xf0f0f0f0f0c7f0f0L, 59));
        SimTK_TEST(   isBitSet(0xf0f0f0f0f0c7f0f0UL, 63)
                   && isBitSet(0xf0f0f0f0f0f0f0c7UL, 0));
        SimTK_TEST(!isBitSet(0xf0f0f0f0f0c7f0f0UL, 59));
    }

}

unsigned long long factorial(unsigned x, unsigned long long result = 1) {
  return x <= 1? result : factorial(x - 1, x * result);
}

// Compute n choose k. Returns 1 if k is zero, 0 if k > n.
unsigned long long choose(unsigned n, unsigned k) {
    if (k > n) return 0;
    return factorial(n)/(factorial(k)*factorial(n-k));
}

void testBitCombination() {

    int comb = -1;
    int count=0;
    while (nextBitCombination(14,6,comb))
        ++count;
    SimTK_TEST(count == 3003); // choose(14,6)

    comb=-1; count=0; std::set<int> combos;
    while (nextBitCombination(4,3,comb)) {
        ++count;
        combos.insert(comb);
    }
    SimTK_TEST(count == 4);
    SimTK_TEST(combos.size()==4);
    SimTK_TEST(combos.count(7)==1);  // 0111
    SimTK_TEST(combos.count(11)==1); // 1011
    SimTK_TEST(combos.count(13)==1); // 1101
    SimTK_TEST(combos.count(14)==1); // 1110
    SimTK_TEST(combos.count(5)==0);  // not 0101 (or anything else)

    long long combl = -1;
    count = 0;
    while (nextBitCombination(39,1,combl))
        ++count;
    SimTK_TEST(count == 39);

    // Start with (8,3) combo 00010110; next should be 00011001.
    combl = 0x16;
    SimTK_TEST(nextBitCombination(8,3,combl));
    SimTK_TEST(combl == 0x19);

    // Choosing k=0 should always return 1 combo, 0.
    int combz=-1;
    SimTK_TEST(nextBitCombination(5,0,combz) && combz==0);
    SimTK_TEST(!nextBitCombination(5,0,combz));
    combz=-1;
    SimTK_TEST(nextBitCombination(0,0,combz) && combz==0);
    SimTK_TEST(!nextBitCombination(0,0,combz));
}

int main() {
    SimTK_START_TEST("TestScalar");

        SimTK_SUBTEST(testIsNaN);
        SimTK_SUBTEST(testIsInf);
        SimTK_SUBTEST(testIsFinite);
        SimTK_SUBTEST(testSignBit);
        SimTK_SUBTEST(testSign);
        SimTK_SUBTEST(testSquareAndCube);
        SimTK_SUBTEST(testIsNumericallyEqual);
        SimTK_SUBTEST(testClamp);
        SimTK_SUBTEST(testStep);
        SimTK_SUBTEST(testDeadband);
        SimTK_SUBTEST(testBitScan);
        SimTK_SUBTEST(testBitCombination);

    SimTK_END_TEST();
}

