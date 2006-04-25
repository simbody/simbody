/* Copyright (c) 2005 Stanford University and Michael Sherman.
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

#include "SimTKcommon.h"
#include "simmatrix/internal/Scalar.h"

#include <iostream>
using std::cout;
using std::endl;

using namespace SimTK;

int main()
{
    const Real one = 1., two = 2.;
    const Complex oneTwo(1.,2.);
    const Complex threeFour(3.,4.); 
    const FComplex fcinf = CNT<FComplex>::getInfinity();
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

    const negator<conjugate<Real> >& negreconj34 =
        reinterpret_cast<const negator<conjugate<Real> >&>(reconj34);
    cout << "negreconj34=" << negreconj34 << endl;
    
    const negator<conjugate<Real> >& nc_threeFour 
        = reinterpret_cast<const negator<conjugate<Real> >&>(threeFour);
    cout << "nc_threeFour=" << nc_threeFour << " conj(.)=" 
        << CNT<negator<conjugate<Real> > >::transpose(nc_threeFour) << endl;

    cout << "NC<C> nan=" << CNT<negator<conjugate<Real> > >::getNaN() << endl;
    cout << "NC<C> inf=" << CNT<negator<conjugate<Real> > >::getInfinity() << endl;

    cout << "negator<FComplex>*long double=" <<
        typeid( negator<FComplex>::Result<long double>::Mul ).name() << endl;
    negator<LComplex> nlc = negator<FComplex>::Result<long double>::Mul(LComplex(1,2));
    cout << "nlc=" << nlc << endl;

    cout << "NegConjugate<DReal>*float=" <<
        typeid( negator<conjugate<DReal> >::Result<float>::Mul ).name() << endl;
    negator<conjugate<DReal> > ncdc = 
        negator<conjugate<DReal> >::Result<float>::Mul(DComplex(9,10));
    cout << "ncdc=" << ncdc << endl;

    cout << "NegConjugate<FReal>*DComplex=" <<
        typeid( negator<conjugate<FReal> >::Result<DComplex>::Mul ).name() << endl;
    negator<DComplex> ndc = 
        negator<conjugate<FReal> >::Result<DComplex>::Mul(DComplex(.1,.2));
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

    return 0; // success
}
