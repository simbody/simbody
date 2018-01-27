/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-17 Stanford University and the Authors.        *
 * Authors: Antoine Falisse                                                   *
 * Contributors: Michael Sherman, Chris Dembia                                *
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
#include <adolc/adolc.h> // for jacobian() ADOL-C driver

#include <iostream>
using std::cout;
using std::endl;

using namespace SimTK;

// Test derivative of simple function with ADOLC without Simbody; just to make
// sure that ADOLC is included properly
void testDerivativeADOLC() {
    double xp[1];
    xp[0] = -2.3;

    trace_on(1);
    adouble x;
    adouble y;
    x <<= xp[0];
    y = 3*pow(x,3)+cos(x)+1;
    double y0;
    y >>= y0;
    trace_off();

    double** J;
    J = myalloc(1,1);
    jacobian(1, 1, 1, xp, J);
    SimTK_TEST(J[0][0] == 9*pow(x,2)-sin(x));
    myfree(J);
}

// Various unit tests verifying that NTraits<adouble> works properly
void testNTraitsADOLC() {
    // Widest
    constexpr bool wfad =
        std::is_same<SimTK::Widest<float, adouble>::Type, adouble>::value;
    constexpr bool wadf =
        std::is_same<SimTK::Widest<adouble, float>::Type, adouble>::value;
    constexpr bool wdad =
        std::is_same<SimTK::Widest<double, adouble>::Type, adouble>::value;
    constexpr bool wadd =
        std::is_same<SimTK::Widest<adouble, double>::Type, adouble>::value;
    constexpr bool wadad =
        std::is_same<SimTK::Widest<adouble, adouble>::Type, adouble>::value;
    SimTK_TEST(wfad);
    SimTK_TEST(wadf);
    SimTK_TEST(wdad);
    SimTK_TEST(wadd);
    SimTK_TEST(wadad);
    // Narrowest
    bool nfad =
        std::is_same<SimTK::Narrowest<float, adouble>::Type, adouble>::value;
    bool nadf =
        std::is_same<SimTK::Narrowest<adouble, float>::Type, adouble>::value;
    bool ndad =
        std::is_same<SimTK::Narrowest<double, adouble>::Type, adouble>::value;
    bool nadd =
        std::is_same<SimTK::Narrowest<adouble, double>::Type, adouble>::value;
    bool nadad =
        std::is_same<SimTK::Narrowest<adouble, adouble>::Type, adouble>::value;
    SimTK_TEST(nfad);
    SimTK_TEST(nadf);
    SimTK_TEST(ndad);
    SimTK_TEST(nadd);
    SimTK_TEST(nadad);
    // isNaN, isFinite, isInf
    adouble xad = -9.45;
    adouble xNaN = SimTK::NaN;
    adouble xInf = SimTK::Infinity;
    SimTK_TEST(isNaN(xNaN));
    SimTK_TEST(!isNaN(xad));
    SimTK_TEST(isFinite(xad));
    SimTK_TEST(!isFinite(xNaN));
    SimTK_TEST(!isFinite(xInf));
    SimTK_TEST(isInf(xInf));
    SimTK_TEST(!isInf(xad));
    // isNumericallyEqual
    double xd = -9.45;
    float xf = (float)-9.45;
    adouble yad = -9;
    int yi = -9;
    std::complex<float> cf(xf,0.);
    std::complex<double> cd(xd,0.);
    SimTK::conjugate<float> cjf(xf,0);
    SimTK::conjugate<double> cjd(xd,0.);
    SimTK_TEST(isNumericallyEqual(xad,xd));
    SimTK_TEST(isNumericallyEqual(xd,xad));
    SimTK_TEST(isNumericallyEqual(xad,xad));
    SimTK_TEST(isNumericallyEqual(xad,xf));
    SimTK_TEST(isNumericallyEqual(xf,xad));
    SimTK_TEST(isNumericallyEqual(yad,yi));
    SimTK_TEST(isNumericallyEqual(yi,yad));
    SimTK_TEST(isNumericallyEqual(cd,xad));
    SimTK_TEST(isNumericallyEqual(xad,cd));
    SimTK_TEST(isNumericallyEqual(cf,xad));
    SimTK_TEST(isNumericallyEqual(xad,cf));
    SimTK_TEST(isNumericallyEqual(cjd,xad));
    SimTK_TEST(isNumericallyEqual(xad,cjd));
    SimTK_TEST(isNumericallyEqual(cjf,xad));
    SimTK_TEST(isNumericallyEqual(xad,cjf));
}

// This test should throw an exception when using value() while taping
void testExceptionTaping() {
    adouble a = 5.;
    double b = NTraits<adouble>::value(a);
    SimTK_TEST(b == 5);

    trace_on(0);
    SimTK_TEST_MUST_THROW_EXC(NTraits<adouble>::value(a),
        SimTK::Exception::ADOLCTapingNotAllowed
    );
    trace_off();
}

// Various unit tests verifying that negator<adouble> works properly
void testNegator() {
    // Test evaluation of simple function and its derivative
    double xp[1];
    xp[0] = 2;
    trace_on(2);
    adouble x;
    adouble y;
    x <<= xp[0];
    auto result = NTraits<adouble>::pow(x,3);
    y = (negator<adouble>&)result;
    double y0;
    y >>= y0;
    trace_off();
    // function evaluation
    double f[1];
    function(2, 1, 1, xp, f);
    SimTK_TEST(f[0] == -8.);
    // derivative evaluation
    double** J;
    J = myalloc(1, 1);
    jacobian(2, 1, 1, xp, J);
    SimTK_TEST(J[0][0] == -3*NTraits<adouble>::pow(x,2));
    myfree(J);
    // isNumericallyEqual
    adouble xd = 9.45;
    SimTK_TEST(isNumericallyEqual(-xd,(negator<adouble>&)xd));
    // isNaN, isFinite, isInf
    adouble xad = -9.45;
    adouble xNaN = SimTK::NaN;
    adouble xInf = SimTK::Infinity;
    SimTK_TEST(isNaN((negator<adouble>&)xNaN));
    SimTK_TEST(!isNaN((negator<adouble>&)xad));
    SimTK_TEST(isFinite((negator<adouble>&)xad));
    SimTK_TEST(!isFinite((negator<adouble>&)xNaN));
    SimTK_TEST(!isFinite((negator<adouble>&)xInf));
    SimTK_TEST(isInf((negator<adouble>&)xInf));
    SimTK_TEST(!isInf((negator<adouble>&)xad));
    // ensure consistent behavior between double and adouble
    double a = 5;
    adouble ad = 5;
    SimTK_TEST((negator<double>)a == (negator<adouble>)ad);
    SimTK_TEST((negator<double>&)a == (negator<adouble>&)ad);
}

// Various unit tests verifying that cast() works properly
void testCast() {
    // cast an adouble in a double
    adouble a = 5.;
    double b = NTraits<adouble>::cast<double>(a);
    SimTK_TEST(b == a);
    // cast an adouble in a double when taping, this should throw an exception
    trace_on(3);
    SimTK_TEST_MUST_THROW_EXC(NTraits<adouble>::cast<double>(a),
        SimTK::Exception::ADOLCTapingNotAllowed
    );
    trace_off();
    // cast an adouble in an adouble when taping
    trace_on(4);
    adouble c = NTraits<adouble>::cast<adouble>(a);
    trace_off();
    SimTK_TEST(c == a);
}

// Various unit tests verifying that operators involving a vector and an
// adouble work properly
void testVec() {
    adouble a = -2;
    adouble b = 2;
    adouble c = -1.5;
    adouble d = -2.8;
    Vec<3,adouble,1> v;
    v[0] = b;
    v[1] = c;
    v[2] = d;
    // multiplication
    Vec<3,adouble,1> vresmr = v*a;
    SimTK_TEST(vresmr[0] == b*a);
    SimTK_TEST(vresmr[1] == c*a);
    SimTK_TEST(vresmr[2] == d*a);
    Vec<3,adouble,1> vresml = a*v;
    SimTK_TEST(vresml[0] == a*b);
    SimTK_TEST(vresml[1] == a*c);
    SimTK_TEST(vresml[2] == a*d);
    // division
    Vec<3,adouble,1> vresdr = v/a;
    SimTK_TEST(vresdr[0] == b/a);
    SimTK_TEST(vresdr[1] == c/a);
    SimTK_TEST(vresdr[2] == d/a);
    // addition
    Vec<3,adouble,1> vresar = v+a;
    SimTK_TEST(vresar[0] == b+a);
    SimTK_TEST(vresar[1] == c+a);
    SimTK_TEST(vresar[2] == d+a);
    Vec<3,adouble,1> vresal = a+v;
    SimTK_TEST(vresal[0] == a+b);
    SimTK_TEST(vresal[1] == a+c);
    SimTK_TEST(vresal[2] == a+d);
    // substraction
    Vec<3,adouble,1> vressr = v-a;
    SimTK_TEST(vressr[0] == b-a);
    SimTK_TEST(vressr[1] == c-a);
    SimTK_TEST(vressr[2] == d-a);
}

// Various unit tests verifying that operators involving a matrix and an
// adouble work properly
void testMat() {
    adouble a = -2;
    adouble b = 2;
    adouble c = -1.5;
    adouble d = -2.8;
    adouble e = 1.87;
    Mat<2,2,adouble,2,1> m;
    m[0][0] = b;
    m[1][0] = c;
    m[0][1] = d;
    m[1][1] = e;
    // multiplication
    Mat<2,2,adouble,2,1> mresmr = m*a;
    SimTK_TEST(mresmr[0][0] == b*a);
    SimTK_TEST(mresmr[1][0] == c*a);
    SimTK_TEST(mresmr[0][1] == d*a);
    SimTK_TEST(mresmr[1][1] == e*a);
    Mat<2,2,adouble,2,1> mresml = a*m;
    SimTK_TEST(mresml[0][0] == a*b);
    SimTK_TEST(mresml[1][0] == a*c);
    SimTK_TEST(mresml[0][1] == a*d);
    SimTK_TEST(mresml[1][1] == a*e);
    // division
    Mat<2,2,adouble,2,1> mresdr = m/a;
    SimTK_TEST(mresdr[0][0] == b/a);
    SimTK_TEST(mresdr[1][0] == c/a);
    SimTK_TEST(mresdr[0][1] == d/a);
    SimTK_TEST(mresdr[1][1] == e/a);
    Mat<2,2,adouble,2,1> mresdl = a/m;
    Mat<2,2,adouble,2,1> minv = a*m.invert();
    SimTK_TEST(mresdl[0][0] == minv[0][0]);
    SimTK_TEST(mresdl[1][0] == minv[1][0]);
    SimTK_TEST(mresdl[0][1] == minv[0][1]);
    SimTK_TEST(mresdl[1][1] == minv[1][1]);
    // addition
    Mat<2,2,adouble,2,1> mresar = m+a;
    SimTK_TEST(mresar[0][0] == b+a);
    SimTK_TEST(mresar[1][0] == m[1][0]);
    SimTK_TEST(mresar[0][1] == m[0][1]);
    SimTK_TEST(mresar[1][1] == e+a);
    Mat<2,2,adouble,2,1> mresal = a+m;
    SimTK_TEST(mresal[0][0] == a+b);
    SimTK_TEST(mresal[1][0] == m[1][0]);
    SimTK_TEST(mresal[0][1] == m[0][1]);
    SimTK_TEST(mresal[1][1] == a+e);
    // substraction
    Mat<2,2,adouble,2,1> mressr = m-a;
    SimTK_TEST(mressr[0][0] == b-a);
    SimTK_TEST(mressr[1][0] == m[1][0]);
    SimTK_TEST(mressr[0][1] == m[0][1]);
    SimTK_TEST(mressr[1][1] == e-a);
}


int main() {
    SimTK_START_TEST("TestADOLCCommon");
        SimTK_SUBTEST(testDerivativeADOLC);
        SimTK_SUBTEST(testNTraitsADOLC);
        SimTK_SUBTEST(testExceptionTaping);
        SimTK_SUBTEST(testNegator);
        SimTK_SUBTEST(testCast);
        SimTK_SUBTEST(testVec);
        SimTK_SUBTEST(testMat);
    SimTK_END_TEST();
}
