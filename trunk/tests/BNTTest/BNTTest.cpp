/* Portions copyright (c) 2005-6 Stanford University and Michael Sherman.
 * Contributors:
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "SimTKcommon/basics.h"
#include "SimTKcommon/internal/Scalar.h"

#include <iostream>
#include <iomanip>
#include <limits>
#include <complex>
using std::cout;
using std::endl;
using std::setprecision;
using std::complex;

using namespace SimTK;

int main()
{
    const Real one = 1., two = 2.;
    const Complex oneTwo(1.,2.);
    const Complex threeFour(3.,4.); 
    const complex<float> fcinf = CNT< complex<float> >::getInfinity();
    const Real nan = CNT<Real>::getNaN();
    
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

    negator<Complex> x(Complex(7.1,1.7)), y;
    y = x; y *= 2.;
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
    printf("PI=%s\n", STRZ(SimTK_PI));
    cout << "f=" << setprecision(std::numeric_limits<float>::digits10+2) << NTraits<float>::Pi
         << " d=" << setprecision(std::numeric_limits<double>::digits10+2) << NTraits<double>::Pi
         << " ld=" << setprecision(std::numeric_limits<long double>::digits10+2) << NTraits<long double>::Pi << endl;
    
    printf("1/sqrt(2)=%.18lg\n", 1/SimTK_SQRT2);
    cout << "f=" << setprecision(std::numeric_limits<float>::digits10+2) << NTraits<float>::OneOverSqrt2
         << " d=" << setprecision(std::numeric_limits<double>::digits10+2) << NTraits<double>::OneOverSqrt2
         << " ld=" << setprecision(std::numeric_limits<long double>::digits10+2) << NTraits<long double>::OneOverSqrt2 << endl;

    printf("Eps f=%.16lg d=%.16lg ld=%.16lg\n",
        (long double)NTraits<float>::Eps, 
        (long double)NTraits<double>::Eps, 
        NTraits<long double>::Eps);

    printf("Eps_13 f=%.16lg d=%.16lg ld=%.16lg\n",
        (long double)NTraits<float>::Eps_13, 
        (long double)NTraits<double>::Eps_13, 
        NTraits<long double>::Eps_13);

    printf("Eps_78 f=%.16lg d=%.16lg ld=%.16lg\n",
        (long double)NTraits<float>::Eps_78, 
        (long double)NTraits<double>::Eps_78, 
        NTraits<long double>::Eps_78);

    printf("Tiny f=%.16lg d=%.16lg ld=%.16lg\n",
        (long double)NTraits<float>::Tiny, 
        (long double)NTraits<double>::Tiny, 
        NTraits<long double>::Tiny);

    return 0; // success
}
