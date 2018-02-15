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
    const short int TraceTag = 2;
    trace_on(TraceTag);
    adouble x, y; // Declare active variables
    x <<= xp[0]; // Select independent variable
    auto result = NTraits<adouble>::pow(x,3);
    // Negate the result by reinterpretation rather than computation
    y = reinterpret_cast<const negator<adouble>&>(result);
    double y0;
    y >>= y0; // Select dependent variable
    trace_off();
    // Function evaluation. function(arg1,arg2,arg3,arg4,arg5) evaluates the
    // desired function from the tape instead of executing the corresponding
    // source code. arg1 is the tape identification, arg2 the number of
    // dependent variables m, arg3 the number of independent variables n, arg4
    // the independent vector x, and arg4 the dependent vector y = f(x) with
    // f: R(n) -> R(m). Please see ADOL-C manual for more details.
    double f[1];
    function(TraceTag, 1, 1, xp, f);
    SimTK_TEST(f[0] == -8.);
    // Derivative evaluation. gradient(arg1,arg2,arg3,arg4) evaluates the
    // the desired function derivative. arg1 is the tape identification, arg2
    // the number of independent variables n (the number of dependent variables
    // m = 1), arg3 the independent vector x, and arg4 the resulting gradient
    // of f(x). Please see ADOL-C manual for more details.
    double g[1];
    gradient(TraceTag, 1, xp, g);
    SimTK_TEST(g[0] == -3*NTraits<adouble>::pow(x,2));
    // isNumericallyEqual
    adouble xd = 9.45;
    auto& nxd = reinterpret_cast<const negator<adouble>&>(xd);
    SimTK_TEST(isNumericallyEqual(-xd,nxd));
    // isNaN, isFinite, isInf
    adouble xad = -9.45;
    adouble xNaN = SimTK::NaN;
    adouble xInf = SimTK::Infinity;
    auto& nxad = reinterpret_cast<const negator<adouble>&>(xad);
    auto& nxNaN = reinterpret_cast<const negator<adouble>&>(xNaN);
    auto& nxInf = reinterpret_cast<const negator<adouble>&>(xInf);
    SimTK_TEST(isNaN(nxNaN));
    SimTK_TEST(!isNaN(nxad));
    SimTK_TEST(isFinite(nxad));
    SimTK_TEST(!isFinite(nxNaN));
    SimTK_TEST(!isFinite(nxInf));
    SimTK_TEST(isInf(nxInf));
    SimTK_TEST(!isInf(nxad));
    // Ensure consistent behavior between double and adouble. The values should
    // be unchanged by conversion to the negated type and negated by
    // reinterpretation as the negated type.
    double a = 5;
    adouble ad = 5;
    SimTK_TEST(static_cast<negator<double>>(a) == a);
    SimTK_TEST(static_cast<negator<adouble>>(ad) == ad);
    SimTK_TEST(reinterpret_cast<negator<double>&>(a) == -a);
    SimTK_TEST(reinterpret_cast<negator<adouble>&>(ad) == -ad);
}

// Various unit tests verifying that cast() works properly
void testCast() {
    // cast an adouble to a double
    adouble a = 5.;
    double b = NTraits<adouble>::cast<double>(a);
    SimTK_TEST(b == a);
    // cast an adouble to a double when taping, this should throw an exception
    trace_on(3);
    SimTK_TEST_MUST_THROW_EXC(NTraits<adouble>::cast<double>(a),
        SimTK::Exception::ADOLCTapingNotAllowed
    );
    trace_off();
    // cast an adouble to an adouble when taping
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

// Various unit tests verifying that functions defined in Scalar work properly
// with adouble
void testScalar() {
    double a = 2;
    double b = -2;
    double c = 0;
    adouble ad = 2;
    adouble bd = -2;
    adouble cd = 0;
    // Negate the result by reinterpretation rather than computation
    auto& na = reinterpret_cast<const negator<double>&>(a);
    auto& nb = reinterpret_cast<const negator<double>&>(b);
    auto& nc = reinterpret_cast<const negator<double>&>(c);
    auto& nad = reinterpret_cast<const negator<adouble>&>(ad);
    auto& nbd = reinterpret_cast<const negator<adouble>&>(bd);
    auto& ncd = reinterpret_cast<const negator<adouble>&>(cd);
    // The following tests ensure a consistent behavior between double and
    // adouble
    // signBit()
    SimTK_TEST(signBit(a) == signBit(ad));
    SimTK_TEST(signBit(b) == signBit(bd));
    SimTK_TEST(signBit(c) == signBit(cd));
    SimTK_TEST(signBit(na) == signBit(nad));
    SimTK_TEST(signBit(nb) == signBit(nbd));
    SimTK_TEST(signBit(nc) == signBit(ncd));
    // sign()
    SimTK_TEST(sign(a) == sign(ad));
    SimTK_TEST(sign(b) == sign(bd));
    SimTK_TEST(sign(c) == sign(cd));
    SimTK_TEST(sign(na) == sign(nad));
    SimTK_TEST(sign(nb) == sign(nbd));
    SimTK_TEST(sign(nc) == sign(ncd));
    // square()
    SimTK_TEST(square(a) == square(ad));
    SimTK_TEST(square(na) == square(nad));
    // cube()
    SimTK_TEST(cube(a) == cube(ad));
    SimTK_TEST(cube(na) == cube(nad));
    // Ensure that recast(), used in cube(const negator<adouble>& x), is
    // differentiable. Function and derivative evaluations should be the same
    // when negating variables by reinterpretaion or by computation.
    double xp[1];
    xp[0] = 2;
    adouble x, y[2];
    const short int TraceTag = 5;
    trace_on(TraceTag);
    x <<= xp[0]; // select independent variables
    y[0] = cube(-x); // negate variable by computation
    // negate variable by reintepretation
    y[1] = cube(reinterpret_cast<const negator<adouble>&>(x));
    double y0[2];
    y[0] >>= y0[0]; // select dependent variables
    y[1] >>= y0[1]; // select dependent variables
    trace_off();
    // Function evaluation
    double f[2];
    function(TraceTag, 2, 1, xp, f);
    SimTK_TEST(f[0] == cube(-xp[0]));
    SimTK_TEST(f[1] == cube(-xp[0]));
    // Derivative evaluation. jacobian(arg1,arg2,arg3,arg4,arg5) evaluates the
    // the desired function derivatives. arg1 is the tape identification, arg2
    // the number of dependent variables m, arg3 the number of independent
    // variables n, arg4 the independent vector x, and arg5 the resulting
    // jacobian of f(x). Please see ADOL-C manual for more details.
    double** J = myalloc(2,1); // initialization and allocation of J, a column
    // of row pointers that will contain the derivatives. J has dimensions
    // m x n where m is the number of dependent variables and n the number of
    // independent variables.
    jacobian(TraceTag, 2, 1, xp, J);
    SimTK_TEST(J[0][0] == -3*square(xp[0]));
    SimTK_TEST(J[1][0] == -3*square(xp[0]));
    myfree(J); // clean memory
    // The following tests ensure a consistent behavior between double and
    // adouble
    int ai = 2;
    int bi = -2;
    // clampInPlace()
    {   double h = 4;
        adouble hd = 4;
        SimTK_TEST(clampInPlace(b,h,a)==clampInPlace(b,hd,a) && hd==a); }
    {   double h = 4;
        adouble hd = 4;
        SimTK_TEST(clampInPlace(b,h,a)==clampInPlace(bd,hd,ad) && hd==a); }
    {   double h = 4;
        adouble hd = 4;
        SimTK_TEST(clampInPlace(bi,h,ai)==clampInPlace(bi,hd,ai) && hd==ai); }
    {   double h = 4;
        adouble hd = 4;
        SimTK_TEST(clampInPlace(bi,h,a)==clampInPlace(bi,hd,a) && hd==a); }
    {   double h = 4;
        adouble hd = 4;
        SimTK_TEST(clampInPlace(b,h,ai)==clampInPlace(b,hd,ai) && hd==ai); }
    {   double nh = -4;
        adouble nhd = -4;
        SimTK_TEST(clampInPlace(b,nh,a)==clampInPlace(b,nhd,a) && nhd==b); }
    {   double nh = -4;
        adouble nhd = -4;
        SimTK_TEST(clampInPlace(b,nh,a)==clampInPlace(bd,nhd,ad) && nhd==b); }
    {   double nh = -4;
        adouble nhd = -4;
        SimTK_TEST(clampInPlace(bi,nh,ai)==clampInPlace(bi,nhd,ai)&&nhd==bi); }
    {   double nh = -4;
        adouble nhd = -4;
        SimTK_TEST(clampInPlace(bi,nh,a)==clampInPlace(bi,nhd,a) && nhd==b); }
    {   double nh = -4;
        adouble nhd = -4;
        SimTK_TEST(clampInPlace(b,nh,ai)==clampInPlace(b,nhd,ai) && nhd==bi); }
    {   double h = 4;
        adouble hd = 4;
        auto nh = reinterpret_cast<const negator<double>&>(h);
        auto nhd = reinterpret_cast<const negator<adouble>&>(hd);
        SimTK_TEST(clampInPlace(b,nh,a)==clampInPlace(b,nhd,a) && nhd==b); }
    {   double h = 4;
        adouble hd = 4;
        auto nh = reinterpret_cast<const negator<double>&>(h);
        auto nhd = reinterpret_cast<const negator<adouble>&>(hd);
        SimTK_TEST(clampInPlace(b,nh,a)==clampInPlace(bd,nhd,ad) && nhd==b); }
    {   double h = 4;
        adouble hd = 4;
        auto nha = reinterpret_cast<const negator<double>&>(h);
        double hb = nha;
        nha = clampInPlace(b,hb,a);
        auto nhb = reinterpret_cast<const negator<double>&>(h);
        auto nhd = reinterpret_cast<const negator<adouble>&>(hd);
        SimTK_TEST(clampInPlace(b,nhb,a)==nha &&
            clampInPlace(b,nhb,a)==clampInPlace(b,nhd,a) && nhd==b); }
    {   double h = 4;
        adouble hd = 4;
        auto nha = reinterpret_cast<const negator<double>&>(h);
        double hb = nha;
        nha = clampInPlace(b,hb,a);
        auto nhb = reinterpret_cast<const negator<double>&>(h);
        auto nhd = reinterpret_cast<const negator<adouble>&>(hd);
        SimTK_TEST(clampInPlace(b,nhb,a)==nha &&
            clampInPlace(b,nhb,a)==clampInPlace(bd,nhd,ad) && nhd==b); }
    // clamp()
    {   double h = 4;
        adouble hd = 4;
        SimTK_TEST(clamp(b,h,a) == clamp(b,hd,a) && hd==4); }
    {   double h = 4;
        adouble hd = 4;
        SimTK_TEST(clamp(b,h,a) == clamp(bd,hd,ad) && hd==4); }
    {   double h = 4;
        adouble hd = 4;
        SimTK_TEST(clamp(bi,h,ai) == clamp(bi,hd,ai) && hd==4); }
    {   double h = 4;
        adouble hd = 4;
        SimTK_TEST(clamp(bi,h,a) == clamp(bi,hd,a) && hd==4); }
    {   double h = 4;
        adouble hd = 4;
        SimTK_TEST(clamp(b,h,ai) == clamp(b,hd,ai) && hd==4); }
    {   double nh = -4;
        adouble nhd = -4;
        SimTK_TEST(clamp(b,nh,a) == clamp(b,nhd,a) && nhd==-4); }
    {   double nh = -4;
        adouble nhd = -4;
        SimTK_TEST(clamp(b,nh,a) == clamp(bd,nhd,ad) && nhd==-4); }
    {   double nh = -4;
        adouble nhd = -4;
        SimTK_TEST(clamp(bi,nh,ai) == clamp(bi,nhd,ai) && nhd==-4); }
    {   double nh = -4;
        adouble nhd = -4;
        SimTK_TEST(clamp(bi,nh,a) == clamp(bi,nhd,a) && nhd==-4); }
    {   double nh = -4;
        adouble nhd = -4;
        SimTK_TEST(clamp(b,nh,ai) == clamp(b,nhd,ai) && nhd==-4); }
    {   double h = 4;
        adouble hd = 4;
        auto nh = reinterpret_cast<const negator<double>&>(h);
        auto nhd = reinterpret_cast<const negator<adouble>&>(hd);
        SimTK_TEST(clamp(b,nh,a)==clamp(b,nhd,a) && nhd == -4); }
    {   double h = 4;
        adouble hd = 4;
        auto nh = reinterpret_cast<const negator<double>&>(h);
        auto nhd = reinterpret_cast<const negator<adouble>&>(hd);
        SimTK_TEST(clamp(b,nh,a)==clamp(bd,nhd,ad) && nhd==-4); }
    // stepUp(), stepDown(), dstepAny(), d2stepUp(), d2stepDown(), d2stepAny(),
    // d3stepUp(), d3stepDown(), d3stepAny()
    double d = 0.2;
    adouble dd = 0.2;
    SimTK_TEST(stepUp(d) == stepUp(dd));
    SimTK_TEST(stepDown(d) == stepDown(dd));
    double e = -1;
    double g = 1;
    adouble ed = -1;
    adouble gd = 1;
    SimTK_TEST(stepAny(e,a,c,g,d) == stepAny(ed,ad,cd,gd,dd));
    SimTK_TEST(dstepUp(d) == dstepUp(dd));
    SimTK_TEST(dstepDown(d) == dstepDown(dd));
    SimTK_TEST(dstepAny(a,c,g,d) == dstepAny(ad,cd,gd,dd));
    SimTK_TEST(d2stepUp(d) == d2stepUp(dd));
    SimTK_TEST(d2stepDown(d) == d2stepDown(dd));
    SimTK_TEST(d2stepAny(a,c,g,d) == d2stepAny(ad,cd,gd,dd));
    SimTK_TEST(d3stepUp(d) == d3stepUp(dd));
    SimTK_TEST(d3stepDown(d) == d3stepDown(dd));
    SimTK_TEST(d3stepAny(a,c,g,d) == d3stepAny(ad,cd,gd,dd));
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
        SimTK_SUBTEST(testScalar);
    SimTK_END_TEST();
}
