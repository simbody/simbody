/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-7 Stanford University and the Authors.         *
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

//#define SimTK_DEFAULT_PRECISION 1
//#define SimTK_DEFAULT_PRECISION 2
//#define SimTK_DEFAULT_PRECISION 4

#include "SimTKcommon/basics.h"
#include "SimTKcommon/internal/Scalar.h"

#include <iostream>
#include <iomanip>
#include <limits>
#include <complex>
#include <cstdio>

using std::cout;
using std::endl;
using std::setprecision;
using std::complex;

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS((cond), "Assertion failed");}

using namespace SimTK;

// use this to keep the compiler from whining about divides by zero
static Real getRealZero();

int main() {
  try
  { const Real zero = 0., one = 1., two = 2.;
    const Complex oneTwo(1.,2.);
    const Complex threeFour(3.,4.); 
    const complex<float> fcinf = CNT< complex<float> >::getInfinity();

        // nan tests
    const Real nan = CNT<Real>::getNaN();
    const float fnan = NTraits<float>::getNaN();
    const double dnan = NTraits<double>::getNaN();
    const long double lnan = NTraits<long double>::getNaN();
    const std::complex<float> cfnan(fnan, fnan);
    const std::complex<double> cdnan(3, dnan);
    const std::complex<long double> clnan(lnan, 0L);
    const conjugate<float> jfnan(fnan, 0.09f);
    const conjugate<double> jdnan(3, dnan);
    const conjugate<long double> jlnan(lnan, lnan);    
    const negator<Real>& nzero = reinterpret_cast<const negator<Real>&>(zero);
    const negator<Real>& ntwo = reinterpret_cast<const negator<Real>&>(two);
    const negator<float>& nfnan = reinterpret_cast<const negator<float>&>(fnan);
    const negator< std::complex<float> >& ncfnan = reinterpret_cast<const negator<std::complex<float> >&>(cfnan);
    const negator< conjugate<long double> >& njlnan = reinterpret_cast<const negator<conjugate<long double> >&>(jlnan);

    ASSERT(isNaN(nan));
    ASSERT(isNaN(-nan));
    ASSERT(isNaN(NaN)); // SimTK::NaN
    ASSERT(isNaN(-NaN)); // SimTK::NaN
    ASSERT(!isNaN(one));
    ASSERT(!isNaN(Infinity));
    ASSERT(isNaN(0./getRealZero()));
    ASSERT(!isNaN(1./getRealZero())); // Infinity

    ASSERT(isNaN(fnan)); ASSERT(isNaN(dnan)); ASSERT(isNaN(lnan)); // float,double,long double
    ASSERT(isNaN(cfnan)); ASSERT(isNaN(cdnan)); ASSERT(isNaN(clnan)); // complex<float,double,long double>
    ASSERT(!isNaN(fcinf)); // complex infinity, not NaN
    ASSERT(isNaN(jfnan)); ASSERT(isNaN(jdnan)); ASSERT(isNaN(jlnan)); // conjugate<float,double,long double>

    // Check negator behavior
    ASSERT(nzero == zero); ASSERT(-nzero == zero);
    ASSERT(ntwo == -two);  ASSERT(-ntwo == two);
    ASSERT(isNaN(nfnan));  ASSERT(isNaN(-nfnan));
    ASSERT(!isNaN(nzero)); ASSERT(!isNaN(-ntwo));
    ASSERT(isNaN(ncfnan)); ASSERT(isNaN(-ncfnan));
    ASSERT(isNaN(njlnan)); ASSERT(isNaN(-njlnan));

    
    cout << "one=" << one << " two=" << two << endl;
    cout << "oneTwo=" << oneTwo << " threeFour=" << threeFour << endl; 
    
    cout << "negator(one)=" << negator<Real>(one) << endl; 
    cout << "reinterp<negator>(one)=" << reinterpret_cast<const negator<Real>&>(one) << endl; 

    cout << "fcinf=" << fcinf << endl;
    cout << "nan=" << nan << endl;

    cout << "negator<Real>::inf=" << CNT<negator<Real> >::getInfinity() << endl;
    cout << "negator<Real>::nan=" << CNT<negator<Real> >::getNaN() << endl;
    //cout << "conj(one)=" << SimTK::conj(one) << " conj(oneTwo)=" << SimTK::conj(oneTwo) << endl;
    
    const conjugate<Real> conj34(threeFour);
    cout << "conj34=" << conj34 << endl;
    cout << "complex(conj34)=" << 
        std::complex<Real>(conj34.real(),conj34.imag()) << endl;
    cout << "-conj34=" << -conj34 << endl;

    const conjugate<Real>& reconj34 = 
        reinterpret_cast<const conjugate<Real>&>(threeFour);
    cout << "reconj34=" << reconj34 << endl;

    const negator<conjugate<Real> > negconj34(conj34);
    cout << "negconj34=" << negconj34 << endl;
    cout << "negconj34.real()=" << negconj34.real() << endl;
    cout << "negconj34.imag()=" << negconj34.imag() << endl;
    const conjugate<Real> crn = (conjugate<Real>)negconj34;
    cout << "Conj(negconj34)=" << crn << endl;

    typedef negator<conjugate<Real> > NCR;
    typedef NCR::TWithoutNegator NCRWN;
    cout << "NCR is a " << typeid(NCR).name() << endl;
    cout << "NCRWN is a " << typeid(NCRWN).name() << endl;

    const NCR& negreconj34 =
        reinterpret_cast<const negator<conjugate<Real> >&>(reconj34);
    cout << "negreconj34=" << negreconj34 << endl;
    cout << " ... is a " << typeid(negreconj34).name() << endl;
    cout << "negreconj34.normalize()=" << negreconj34.normalize() << endl;
    cout << " ... is a " << typeid(negreconj34.normalize()).name() << endl;
    cout << "(noneg)negreconj34=" 
         << CNT<NCR>::castAwayNegatorIfAny(negreconj34) << endl;
    cout << " ... is a "
         << typeid(CNT<NCR>::castAwayNegatorIfAny(negreconj34)).name() << endl;
    const NCRWN& nnn = CNT<NCR>::castAwayNegatorIfAny(negreconj34);
     cout << "(noneg)negreconj34.normalize()=" 
         << CNT<NCRWN>::normalize(nnn) << endl;
   
    const negator<conjugate<Real> >& nc_threeFour 
        = reinterpret_cast<const negator<conjugate<Real> >&>(threeFour);
    cout << "nc_threeFour=" << nc_threeFour << " conj(.)=" 
        << CNT<negator<conjugate<Real> > >::transpose(nc_threeFour) << endl;

    cout << "NC<C> nan=" << CNT<negator<conjugate<Real> > >::getNaN() << endl;
    cout << "NC<C> inf=" << CNT<negator<conjugate<Real> > >::getInfinity() << endl;

    cout << "negator<complex<float>>*long double=" <<
        typeid( negator< complex<float> >::Result<long double>::Mul ).name() << endl;
    negator< complex<long double> > nlc = 
        negator< complex<float> >::Result<long double>::Mul(complex<long double>(1,2));
    cout << "nlc=" << nlc << endl;

    cout << "NegConjugate<double>*float=" <<
        typeid( negator<conjugate<double> >::Result<float>::Mul ).name() << endl;
    negator<conjugate<double> > ncdc = 
        negator< conjugate<double> >::Result<float>::Mul(complex<double>(9,10));
    cout << "ncdc=" << ncdc << endl;

    cout << "NegConjugate<float>*complex<double>=" <<
        typeid( negator<conjugate<float> >::Result< complex<double> >::Mul ).name() << endl;
    negator< complex<double> > ndc = 
        negator<conjugate<float> >::Result< complex<double> >::Mul(complex<double>(.1,.2));
    cout << "ndc=" << ndc << endl;

    negator<Complex> x(Complex(Real(7.1),Real(1.7))), y;
    y = x; y *= Real(2);
    cout << "x=" << x << "(y=x)*=2. =" << y << endl;
    cout << "x*2.=" << x*2. << endl;
    cout << "x*y=" << x*y << endl;
    cout << "x+y=" << x+y << endl;
    cout << "x+(-y)=" << x+(-y) << endl;
    cout << "x+2=" << x+Complex(2,0) << endl;
    Complex zz = operator+(x,Complex(2,0));
    cout << "zz=op+(x,2)=" << zz << endl;

    cout << "negator<Complex> x=" << x << endl;
    cout << "square(x)=" << square(x) << " x*x=" << x*x << " diff=" << square(x)-x*x <<endl;
    cout << "cube(x)=" << cube(x) << " x*x*x=" << x*x*x << " diff=" << cube(x)-x*x*x <<endl;

    Real pp=27, nn=-14, zzz=0;
    cout << "Real: sign(27)=" << sign(pp) << " sign(-14)=" << sign(nn) << " sign(0)=" << sign(zzz) << endl;
    cout << "negator<Real>: sign(27)=" 
         <<  sign(negator<Real>::recast(pp)) << " sign(-14)=" << sign(negator<Real>::recast(nn)) 
         << " sign(0)=" << sign(negator<Real>::recast(zzz)) << endl;


    // Check mixed-mode complex & conjugate operators
    complex<float> cff;
    complex<double> dff;
    complex<long double> lff;
    cff * 3.; 3.*cff; cff /3.; 3./cff; cff + 3.;3.+cff;cff-3.;3.-cff;
    dff * 3.f; 3.f*dff; dff / 3.f; 3.f/dff;dff + 3.f;3.f+dff;dff-3.f;3.f-dff;
    lff * 3.; 3.*lff; lff/3.; 3./lff;lff + 3.;3.+lff;lff-3.;3.-lff;
    cff*3;lff/3;3/dff;3+cff;
    conjugate<float> ccf;
    conjugate<double> dcf;
    conjugate<long double> lcf;
    ccf * 3.; 3.*ccf; ccf /3.; 3./ccf; ccf + 3.;3.+ccf;ccf-3.;3.-ccf;
    dcf * 3.f; 3.f*dcf; dcf / 3.f; 3.f/dcf;dcf + 3.f;3.f+dcf;dcf-3.f;3.f-dcf;
    lcf * 3.; 3.*lcf; lcf/3.; 3./lcf;lcf + 3.;3.+lcf;lcf-3.;3.-lcf;

    // Constants in various precisions
#define STRZ_(X) #X
#define STRZ(X) STRZ_(X)


    printf("\nCONSTANTS IN DEFAULT REAL PRECISION\n");
    printf("NumDigits=%d, LosslessNumDigits=%d\n", NumDigitsReal, LosslessNumDigitsReal);
    cout << "e^(i*pi)+1=" << std::pow(E, I*Pi)+1 << endl;
    cout << "NaN=" << setprecision(LosslessNumDigitsReal) << NaN << endl;
    cout << "Infinity=" << setprecision(LosslessNumDigitsReal) << Infinity << endl;

    cout << "Eps=" << setprecision(LosslessNumDigitsReal) << Eps << endl;
    cout << "SqrtEps=" << setprecision(LosslessNumDigitsReal) << SqrtEps << endl;
    cout << "TinyReal=" << setprecision(LosslessNumDigitsReal) << TinyReal << endl;
    cout << "SignificantReal=" << setprecision(LosslessNumDigitsReal) << SignificantReal << endl;
    cout << "LeastPositiveReal=" << setprecision(LosslessNumDigitsReal) << LeastPositiveReal << endl;
    cout << "MostPositiveReal=" << setprecision(LosslessNumDigitsReal) << MostPositiveReal << endl;
    cout << "LeastNegativeReal=" << setprecision(LosslessNumDigitsReal) << LeastNegativeReal << endl;
    cout << "MostNegativeReal=" << setprecision(LosslessNumDigitsReal) << MostNegativeReal << endl;
    cout << "Zero=" << setprecision(LosslessNumDigitsReal) << Zero << endl;
    cout << "One=" << setprecision(LosslessNumDigitsReal) << One << endl;
    cout << "MinusOne=" << setprecision(LosslessNumDigitsReal) << MinusOne << endl;
    cout << "Two=" << setprecision(LosslessNumDigitsReal) << Two << endl;
    cout << "Three=" << setprecision(LosslessNumDigitsReal) << Three << endl;
    cout << "OneHalf=" << setprecision(LosslessNumDigitsReal) << OneHalf << endl;
    cout << "OneThird=" << setprecision(LosslessNumDigitsReal) << OneThird << endl;
    cout << "OneFourth=" << setprecision(LosslessNumDigitsReal) << OneFourth << endl;
    cout << "OneFifth=" << setprecision(LosslessNumDigitsReal) << OneFifth << endl;
    cout << "OneSixth=" << setprecision(LosslessNumDigitsReal) << OneSixth << endl;
    cout << "OneSeventh=" << setprecision(LosslessNumDigitsReal) << OneSeventh << endl;
    cout << "OneEighth=" << setprecision(LosslessNumDigitsReal) << OneEighth << endl;
    cout << "OneNinth=" << setprecision(LosslessNumDigitsReal) << OneNinth << endl;
    cout << "Pi=" << setprecision(LosslessNumDigitsReal) << Pi << endl;
    cout << "OneOverPi=" << setprecision(LosslessNumDigitsReal) << OneOverPi << endl;
    cout << "E=" << setprecision(LosslessNumDigitsReal) << E << endl;
    cout << "Log2E=" << setprecision(LosslessNumDigitsReal) << Log2E << endl;
    cout << "Log10E=" << setprecision(LosslessNumDigitsReal) << Log10E << endl;
    cout << "Sqrt2=" << setprecision(LosslessNumDigitsReal) << Sqrt2 << endl;
    cout << "OneOverSqrt2=" << setprecision(LosslessNumDigitsReal) << OneOverSqrt2 << endl;
    cout << "Sqrt3=" << setprecision(LosslessNumDigitsReal) << Sqrt3 << endl;
    cout << "OneOverSqrt3=" << setprecision(LosslessNumDigitsReal) << OneOverSqrt3 << endl;
    cout << "CubeRoot2=" << setprecision(LosslessNumDigitsReal) << CubeRoot2 << endl;
    cout << "CubeRoot3=" << setprecision(LosslessNumDigitsReal) << CubeRoot3 << endl;
    cout << "Ln2=" << setprecision(LosslessNumDigitsReal) << Ln2 << endl;
    cout << "Ln10=" << setprecision(LosslessNumDigitsReal) << Ln10 << endl;
    cout << "I=" << setprecision(LosslessNumDigitsReal) << I << endl;

    printf("\nSOME CONSTANTS IN VARIOUS PRECISIONS\n");

    printf("PI=%s\n", STRZ(SimTK_PI));
    cout << "f=" << setprecision(NTraits<float>::getNumDigits()+2) << NTraits<float>::getPi()
         << " d=" << setprecision(NTraits<double>::getNumDigits()+2) << NTraits<double>::getPi()
         << " ld=" << setprecision(NTraits<long double>::getNumDigits()+2) << NTraits<long double>::getPi() << endl;
    
    std::printf("1/sqrt(2)=%.18Lg\n", 1/SimTK_SQRT2);
    cout << "f=" << setprecision(NTraits<float>::getNumDigits()+2) << NTraits<float>::getOneOverSqrt2()
         << " d=" << setprecision(NTraits<double>::getNumDigits()+2) << NTraits<double>::getOneOverSqrt2()
         << " ld=" << setprecision(NTraits<long double>::getNumDigits()+2) << NTraits<long double>::getOneOverSqrt2() << endl;

    printf("Eps f=%.16Lg d=%.16Lg ld=%.16Lg\n",
        (long double)NTraits<float>::getEps(), 
        (long double)NTraits<double>::getEps(), 
        NTraits<long double>::getEps());

    printf("SqrtEps f=%.16Lg d=%.16Lg ld=%.16Lg\n",
        (long double)NTraits<float>::getSqrtEps(), 
        (long double)NTraits<double>::getSqrtEps(), 
        NTraits<long double>::getSqrtEps());

    printf("Significant f=%.16Lg d=%.16Lg ld=%.16Lg\n",
        (long double)NTraits<float>::getSignificant(), 
        (long double)NTraits<double>::getSignificant(), 
        NTraits<long double>::getSignificant());

    printf("Tiny f=%.16Lg d=%.16Lg ld=%.16Lg\n",
        (long double)NTraits<float>::getTiny(), 
        (long double)NTraits<double>::getTiny(), 
        NTraits<long double>::getTiny());

  } catch(const std::exception& e) {
      std::cout << "exception: " << e.what() << std::endl;
      return 1;
  }

    return 0; // success
}

// Try to prevent a smart optimizer from noticing this zero.
static Real getRealZero() {
    return std::sin(Real(0));
}
