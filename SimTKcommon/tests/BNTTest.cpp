/* -------------------------------------------------------------------------- *
 *                       SimTK Simbody: SimTKcommon                           *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-12 Stanford University and the Authors.        *
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

//#define SimTK_DEFAULT_PRECISION 1
//#define SimTK_DEFAULT_PRECISION 2
//#define SimTK_DEFAULT_PRECISION 4

#include "SimTKcommon.h"

#include <iostream>
#include <iomanip>
#include <limits>
#include <complex>
#include <cstdio>
#include <cmath>

using std::cout;
using std::endl;
using std::setprecision;
using std::complex;

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS((cond), "Assertion failed");}
// Assert that two real-valued calculations should produce the same result to within a few 
// machine roundoff errors.
#define ASSERT_EPS(x,y) {ASSERT(std::abs((x)-(y))<=(4*Eps));}
// Assert that two complex-valued calculations should produce the same result (in both real and
// imaginary components) to within a few machine roundoff errors.
#define ASSERT_EPSX(x,y) {ASSERT_EPS((x).real(),(y).real()); \
                          ASSERT_EPS((x).imag(),(y).imag());}

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
    const std::complex<float> cfnan(fnan, fnan);
    const std::complex<double> cdnan(3, dnan);
    const conjugate<float> jfnan(fnan, 0.09f);
    const conjugate<double> jdnan(3, dnan);
    const negator<Real>& nzero = reinterpret_cast<const negator<Real>&>(zero);
    const negator<Real>& ntwo = reinterpret_cast<const negator<Real>&>(two);
    const negator<float>& nfnan = reinterpret_cast<const negator<float>&>(fnan);
    const negator< std::complex<float> >& ncfnan = reinterpret_cast<const negator<std::complex<float> >&>(cfnan);

    writeUnformatted(std::cout, fcinf);
    Array_< negator<Complex> > arrc;
    arrc.push_back(oneTwo); 
    arrc.push_back(threeFour);
    arrc.push_back(fcinf);
    writeUnformatted(cout, arrc); cout << endl;
    writeUnformatted(cout, Vec3(1,2,3)); cout << endl;
    Vector vxxx(Vec3(4,NaN,-3));
    writeUnformatted(cout, vxxx); cout << endl;
    cout << vxxx << "\n";
    writeUnformatted(cout, Mat34( 1, 2, 3, 4,
                              5, NaN, 7, 8,
                              Infinity, 10, -Infinity, 12 ));
    cout << endl;
    writeUnformatted(cout, SymMat33( 1, 
                                 2, 3,
                                 NaN, 5, 6 )); cout << endl;

    Matrix mxxx(Mat34( 1, 2, 3, 4,
                       5, NaN, 7, 8,
                       Infinity, 10, -Infinity, 12 ));
    writeUnformatted(cout, mxxx); cout << endl;
    cout << mxxx;

    double inval;
    std::stringstream ss("1.3  nan 4  6 -3 -inf");
    while (readUnformatted(ss, inval))
        cout << "'" << String(inval) << "'\n";

    std::stringstream ss2("1 2 3 nan 4 inf");
    readUnformatted(ss2, arrc);
    cout << "arrc=" << arrc << endl;
    writeUnformatted(cout, arrc);

    ss2.clear(); ss2.seekg(0, std::ios::beg);
    Vec6 myv6;
    readUnformatted(ss2, myv6);
    writeUnformatted(cout << "myv6=", myv6);
    ss2.clear(); ss2.seekg(0, std::ios::beg);
    Matrix mym23(2,3);
    fillUnformatted(ss2, mym23);
    writeUnformatted(cout << "\nmym23=", mym23);
    ss2.clear(); ss2.seekg(0, std::ios::beg);
    Array_<float> axx(2), ayy(4);
    ArrayView_<float> avxx(axx), avyy(ayy);
    readUnformatted(ss2, avxx);
    readUnformatted(ss2, avyy);
    writeUnformatted(cout << "u axx=",axx); cout<<endl;
    writeUnformatted(cout << "u ayy=",ayy); cout<<endl;
    writeFormatted(cout << "f axx=",axx); cout<<endl;
    writeFormatted(cout << "f ayy=",ayy); cout<<endl;

    ASSERT(isNaN(nan));
    ASSERT(isNaN(-nan));
    ASSERT(isNaN(NaN)); // SimTK::NaN
    ASSERT(isNaN(-NaN)); // SimTK::NaN
    ASSERT(!isNaN(one));
    ASSERT(!isNaN(Infinity));
    ASSERT(isNaN(0./getRealZero()));
    ASSERT(!isNaN(1./getRealZero())); // Infinity

    ASSERT(isNaN(fnan)); ASSERT(isNaN(dnan)); // float,double
    ASSERT(isNaN(cfnan)); ASSERT(isNaN(cdnan)); // complex<float,double>
    ASSERT(!isNaN(fcinf)); // complex infinity, not NaN
    ASSERT(isNaN(jfnan)); ASSERT(isNaN(jdnan)); // conjugate<float,double>

    // Check negator behavior
    ASSERT(nzero == zero); ASSERT(-nzero == zero);
    ASSERT(ntwo == -two);  ASSERT(-ntwo == two);
    ASSERT(isNaN(nfnan));  ASSERT(isNaN(-nfnan));
    ASSERT(!isNaN(nzero)); ASSERT(!isNaN(-ntwo));
    ASSERT(isNaN(ncfnan)); ASSERT(isNaN(-ncfnan));

    
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
    cout << "x-y=" << x-y << endl;
    cout << "x+(-y)=" << x+(-y) << endl;

    // In gcc 4.1.2, if you remove this output line then the
    // corresponding ASSERT below it will fail! 
    cout << "-(-x)+y=" << -(-x)+y << endl;
    ASSERT_EPSX(x+y, -(-x)+y);

    ASSERT_EPSX(x+y, -((-x)+(-y)));
    ASSERT_EPSX(x+y, x-(-y));
    ASSERT_EPSX(x-y, x+(-y));
    ASSERT_EPSX(x-y, -(y-x));
    ASSERT_EPSX(x-y, -(-x+y));
    ASSERT_EPSX(-(x+y), (-x)-y);
    ASSERT_EPSX(-(x+y), -x-y);

    ASSERT_EPSX(x*y, (-x)*(-y));
    ASSERT_EPSX(x*y, -x*-y);
    ASSERT_EPSX(-(x*y), -x*y);
    ASSERT_EPSX(-(x*y), x*-y);

    ASSERT_EPSX(x/y, (-x)/(-y));
    ASSERT_EPSX(x/y, -x/-y);
    ASSERT_EPSX(-(x/y), -x/y);
    ASSERT_EPSX(-(x/y), x/-y);

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
    cff * 3.; 3.*cff; cff /3.; 3./cff; cff + 3.;3.+cff;cff-3.;3.-cff;
    dff * 3.f; 3.f*dff; dff / 3.f; 3.f/dff;dff + 3.f;3.f+dff;dff-3.f;3.f-dff;
    cff*3;3/dff;3+cff;
    conjugate<float> ccf;
    conjugate<double> dcf;
    ccf * 3.; 3.*ccf; ccf /3.; 3./ccf; ccf + 3.;3.+ccf;ccf-3.;3.-ccf;
    dcf * 3.f; 3.f*dcf; dcf / 3.f; 3.f/dcf;dcf + 3.f;3.f+dcf;dcf-3.f;3.f-dcf;

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
         << " d=" << setprecision(NTraits<double>::getNumDigits()+2) << NTraits<double>::getPi() << endl;
    
    std::printf("1/sqrt(2)=%.18Lg\n", 1/SimTK_SQRT2);
    cout << "f=" << setprecision(NTraits<float>::getNumDigits()+2) << NTraits<float>::getOneOverSqrt2()
         << " d=" << setprecision(NTraits<double>::getNumDigits()+2) << NTraits<double>::getOneOverSqrt2() << endl;

    printf("Eps f=%.16g d=%.16g\n",
        (double)NTraits<float>::getEps(), 
        (double)NTraits<double>::getEps());

    printf("SqrtEps f=%.16g d=%.16g\n",
        (double)NTraits<float>::getSqrtEps(), 
        (double)NTraits<double>::getSqrtEps());

    printf("Significant f=%.16g d=%.16g\n",
        (double)NTraits<float>::getSignificant(), 
        (double)NTraits<double>::getSignificant());

    printf("Tiny f=%.16g d=%.16g\n",
        (double)NTraits<float>::getTiny(), 
        (double)NTraits<double>::getTiny());

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
