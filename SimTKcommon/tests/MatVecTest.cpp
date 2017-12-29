/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
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

#include "SimTKcommon/SmallMatrix.h"
#include "SimTKcommon/Testing.h"

#include <cstdio>
#include <iostream>
#include <iomanip>
using std::cout;
using std::endl;
using std::setw;

#include <complex>
using std::complex;
using std::sin;
using std::cos;

using namespace SimTK;

static void dummy();

//typedef Vec<3,Complex> CVec3;
//static CVec3 f(CVec3 v) {
 //   return CVec3();
//}

static Complex f(Complex x) {
    return std::sin(x);
}

// TODO: move lots of tests here and check the answers
void testNegator() {
    Vec3 vvv(1,2,3); Row3 rrr(1,2,3);
    SimTK_TEST(vvv-vvv == Vec3(0));     // these are exact results
    SimTK_TEST(-vvv-vvv == (-2*vvv));
    SimTK_TEST(rrr-rrr == Row3(0));
    SimTK_TEST(-rrr-rrr == (-2*rrr));
}

void testElementwiseOps() {
    Vec3 vvv(1,2,3); Vec3 www(7,9,2);
    Row3 rrr(1,2,3); Row3 sss(5,4,10);

    Vec3 mv=vvv.elementwiseMultiply(www);
    Vec3 dv=vvv.elementwiseDivide(www);
    SimTK_TEST_EQ(mv, Vec3(7,18,6));
    SimTK_TEST_EQ(dv, Vec3(Real(1)/7, Real(2)/9, Real(3)/2));

    Row3 mr=rrr.elementwiseMultiply(sss);
    Row3 dr=rrr.elementwiseDivide(sss);
    SimTK_TEST_EQ(mr, Row3(5,8,30));
    SimTK_TEST_EQ(dr, Row3(Real(1)/5, Real(2)/4, Real(3)/10));

    Mat22 mmm(1, 2,
              3, 4);
    Mat22 nnn(7, 9,
              2, 3);
    SymMat22 yyy(1,
                 2, 3);
    SymMat22 zzz(5,
                 4, 10);

    Mat22 mm=mmm.elementwiseMultiply(nnn);
    Mat22 dm=mmm.elementwiseDivide(nnn);
    SymMat22 my=yyy.elementwiseMultiply(zzz);
    SymMat22 dy=yyy.elementwiseDivide(zzz);

    SimTK_TEST_EQ(mm, Mat22(7,18,6,12));
    SimTK_TEST_EQ(dm, Mat22(Real(1)/7,Real(2)/9,Real(3)/2,Real(4)/3));
    SimTK_TEST_EQ(my, SymMat22(5,8,30));
    SimTK_TEST_EQ(dy, SymMat22(Real(1)/5,Real(2)/4,Real(3)/10));
}

void testSums() {
    Mat22 m(1,2,
            3,4);
    SimTK_TEST_EQ(m.colSum(), Row2(4,6));
    SimTK_TEST_EQ(m.rowSum(), Vec2(3,7));
    SimTK_TEST(m.sum() == m.colSum()); // should be exact

    SymMat22 y(3, /*4*/
               4, 5);
    SimTK_TEST_EQ(y.colSum(), Row2(7,9)); // same for real, sym
    SimTK_TEST_EQ(y.rowSum(), Vec2(7,9));
    SimTK_TEST(y.sum() == y.colSum()); // should be exact

    Mat22 sm(y); // create fully populated symmetric matrix
    SimTK_TEST_EQ(sm.rowSum(), y.rowSum());
    SimTK_TEST_EQ(sm.colSum(), y.colSum());

    Mat<2,2,Complex> mc(1+2*I, 3+4*I,
                        5+6*I, 7+8*I);
    typedef Row<2,Complex> CRow2;
    typedef Vec<2,Complex> CVec2;
    SimTK_TEST_EQ(mc.colSum(), CRow2(6+8*I, 10+12*I));
    SimTK_TEST_EQ(mc.rowSum(), CVec2(4+6*I, 12+14*I));
    SimTK_TEST(mc.sum() == mc.colSum()); // should be exact

    // Row sum and col sum for Hermitian are conjugates; not the same.
    SymMat<2,Complex> yc( 1, /*3+6*I*/
                         3-6*I, 4);
    SimTK_TEST_EQ(yc.colSum(), CRow2(4-6*I, 7+6*I));
    SimTK_TEST_EQ(yc.rowSum(), CVec2(4+6*I, 7-6*I));
    SimTK_TEST(yc.sum() == yc.colSum()); // should be exact

    Mat<2,2,Complex> smc(yc); // create fully populated symmetric matrix
    SimTK_TEST_EQ(smc.rowSum(), yc.rowSum());
    SimTK_TEST_EQ(smc.colSum(), yc.colSum());
}

void testMiscellaneous()
{
    cout << std::setprecision(16);
    cout << "f(.3)=" << f(Complex(0.3)) << endl;
    Real h = 1e-20;
    cout << "f(.3 + i*h)/h=" << f(Complex(0.3,h)) / h << endl;

    cout << CNT< Mat<2,3, Vec<2> > >::getNaN() << endl;

    Mat<2,3, Vec<2, Mat<2,2,Complex> > > isThisNaN;
    cout << "isThisNan? " << isThisNaN << endl;

    const Complex mdc[] = {
        Complex(1.,2.),  Complex(3.,4.),   Complex(5.,6.),   Complex(7.,8.),
        Complex(9.,10.), Complex(10.,11.), Complex(.1,.26),  Complex(.3,.45),
        Complex(.5,.64), Complex(.7,.83),  Complex(.9,.102), Complex(.10,.111)   
    }; 

    cout << "*** TEST COMPLEX DOT PRODUCT ***" << endl;
    Vec<3,Complex> vdot(mdc), wdot(&mdc[3]);
    Row<3,Complex> rdot(mdc), sdot(&mdc[3]);
    cout << "v=" << vdot << " w=" << wdot << endl;
    cout << "r=" << rdot << " s=" << sdot << endl;
    cout << "v.normalize()=" << vdot.normalize() << endl;
    cout << "r.normalize()=" << rdot.normalize() << endl;


    cout << "--- dot() global function:dot(v,w), rw, vs, rs should be the same" << endl;
    cout << "vw=" << dot(vdot,wdot) << " rw" << dot(rdot,wdot) 
         << " vs" << dot(vdot,sdot) << " rs" << dot(rdot,sdot) << endl;
    cout << "--- dot operator* requires row*col meaning Hermitian transpose with sign changes" << endl;
    cout << "vw=" << ~vdot*wdot << " rw" << rdot*wdot
         << " vs" << ~vdot*~sdot << " rs" << rdot*~sdot << endl;

    cout << endl << "*** TEST COMPLEX OUTER PRODUCT ***" << endl;
    cout << "--- outer() global function:dot(v,w), rw, vs, rs should be the same" << endl;
    cout << "vw=" << outer(vdot,wdot) << " rw" << outer(rdot,wdot) 
         << " vs" << outer(vdot,sdot) << " rs" << outer(rdot,sdot) << endl;
    cout << "--- outer operator* requires col*row meaning Hermitian transpose with sign changes" << endl;
    cout << "vw=" << vdot*~wdot << " rw" << ~rdot*~wdot
         << " vs" << vdot*sdot << " rs" << ~rdot*sdot << endl;

    cout << "*** TEST COMPLEX CROSS PRODUCT ***" << endl;
    cout << "--- cross() global function:dot(v,w), rw, vs, rs should be the same" << endl;
    cout << "vw=" << cross(vdot,wdot) << " rw" << cross(rdot,wdot) 
         << " vs" << cross(vdot,sdot) << " rs" << cross(rdot,sdot) << endl;
    cout << "--- cross operator% involves NO sign changes, but returns row if either arg is a row" << endl;
    cout << "vw=" << vdot%wdot << " rw" << rdot%wdot
         << " vs" << vdot%sdot << " rs" << rdot%sdot << endl;

    cout << "*** TEST crossMat() ***" << endl;
    Mat<3,3,Complex> vcross(crossMat(vdot));
    Mat<3,3,Complex> rcross(crossMat(rdot));
    cout << "--- crossMat 3d should be same whether made from row or vec" << endl;
    cout << "vdot%wdot=" << vdot%wdot << endl;
    cout << "crossMat(v)=" << vcross << "crossMat(r)=" << rcross;
    cout << "crossMat(v)*w=" << crossMat(vdot)*wdot << " vcross*w=" << vcross*wdot << endl;


    Vec<2,Complex> vdot2 = vdot.getSubVec<2>(0);
    Vec<2,Complex> wdot2 = wdot.getSubVec<2>(0);
    Row<2,Complex> rdot2 = rdot.getSubRow<2>(0);
    Row<2,Complex> vcross2(crossMat(vdot2));
    Row<2,Complex> rcross2(crossMat(rdot2));

    cout << "--- crossMat 2d should be same whether made from row or vec" << endl;
    cout << "vdot2, wdot2=" << vdot2 << ", " << wdot2 << " vdot2%wdot2=" << vdot2%wdot2 << endl;
    cout << "crossMat(v2)=" << vcross2 << "crossMat(r2)=" << rcross2;
    cout << "crossMat(v2)*w2=" << crossMat(vdot2)*wdot2 << " vcross2*w2=" << vcross2*wdot2 << endl;

    cout << "*********\n";



    Mat<2,5,float> m25f( 1, 2, 3, 4, 5,
                         6, 7, 8, 9, 10 );
    cout << "Mat<2,5,float>=" << m25f;
    cout << "Mat<2,5,float>.normalize()=" << m25f.normalize();
    cout << "Mat<2,5,float>.sqrt()=" << m25f.sqrt();

    const Mat<1,5,Vec<2,float> >& m15v2f = 
        *reinterpret_cast<const Mat<1,5,Vec<2,float> >*>(&m25f);
    cout << "  m25f@" << &m25f << " m15v2f@" << &m15v2f << endl;
    cout << "Mat<1,5,Vec<2,float> >=" << m15v2f;;
    cout << "Mat<1,5,Vec<2,float> >.normalize()=" << m15v2f.normalize();

    const Real twoXthree[] = { 1, 2, 3,
                               4, 5, 6 };
    const Real threeXone[] = { .1, .001, .00001 };
    const Real sym33[] = { 1,
                           2, 3,
                           4, -5, 6 };
    SymMat<3> sm3(sym33);
    cout << "SymMat<3> sm3=" << sm3;
    Mat<3,3> m33sm3;
    for (int i=0; i<3; ++i) 
        for (int j=0; j<=i; ++j)
            m33sm3(i,j) = m33sm3(j,i) = sm3(i,j);
    cout << "Mat33(sm3)=" << m33sm3;
    cout << "sm3*3=" << sm3*3;
    cout << "sm3+100=" << sm3+100;
    cout << "m33sm3+100=" << m33sm3+100;
    cout << "sm3.normalize()=" << sm3.normalize();
    cout << "m33sm3.normalize()=" << m33sm3.normalize();
    cout << "sm3+=100:" << (sm3+=100.);
    cout << "m33sm3+=100:" << (m33sm3+=100.);

    Mat<3,3,Complex> whole(mdc);
    SymMat<3,Complex,9> sym = SymMat<3,Complex,9>().setFromLower(whole);
    cout << "whole=" << whole << endl;
    cout << "sym  =" << sym << "(pos~)sym  =" << sym.positionalTranspose() << endl;

    cout << "whole.real()=" << whole.real();
    cout << "whole.imag()=" << whole.imag();
    cout << "sym.real()=" << sym.real();
    cout << "sym.imag()=" << sym.imag();



    Mat<3,4,Complex>  mdcp(mdc);  cout << "*** Data looks like this: " << mdcp;
    SymMat<4,negator<Complex> > symp(reinterpret_cast<const negator<conjugate<double> >*>(mdc));
    cout << "    4x4 Sym<Neg<cmplx>> from (negator<conj>)pointer to data gives this:" << symp;
    cout << "    sym.real()=" << symp.real();
    cout << "    sym.imag()=" << symp.imag();
    cout << "   ~sym.imag()=" << ~symp.imag();
    cout << "pos~(sym.imag())=" << symp.imag().positionalTranspose();
    cout << "(pos~sym).imag()=" << symp.positionalTranspose().imag();
    cout << "   -sym.imag()=" << -symp.imag();

    symp(2,1).real() = 99.;
    cout << "after sym(2,1).real=99, sym=" << symp;

    symp.updPositionalTranspose().imag()(3,1)=123.;
    cout << "after (pos~sym).imag()(3,1)=123, (pos~sym).imag()=" << symp.positionalTranspose().imag();
    cout << "    ... sym=" << symp;

    Mat<2,3, SymMat<3,Complex> > weird(Row<3,SymMat<3,Complex> >( sym, -sym, sym ),
                                       Row<3,SymMat<3,Complex> >( sym, sym, sym ));
    cout << "weird=" << weird;
    weird *= 2.;
    cout << " weird*=2: " << weird;
    cout << " weird(1)=" << weird(1) << endl;
    cout << " weird(0,1)=" << weird(0,1) << " [0][1]=" << weird[0][1] << endl;

    cout << " typename(weird)=" << typeid(weird).name() << endl;
    cout << " typename(weird.real)=" << typeid(weird.real()).name() << endl;
    cout << " typename(weird.imag)=" << typeid(weird.imag()).name() << endl;

    cout << " weird.real()=" << weird.real();
    cout << " weird.imag()=" << weird.imag();

    cout << "sizeof(sym<3,cplx>)=" << sizeof(sym) << " sizeof(mat<2,3,sym>=" << sizeof(weird) << endl;

    Mat<2,3> m23(twoXthree);
    Mat<3,1> m31(threeXone);
    cout << "m23=" << m23 << endl;
    cout << "m31=" << m31 << endl;
    cout << "m23*-m31=" << m23*-m31 << endl;
    cout << "~ ~m31 * ~-m23=" << ~((~m31)*(~-m23)) << endl;

    Mat<2,3,Complex> c23(m23);
    Mat<3,1,Complex> c31(m31);
    cout << "c23=" << c23 << endl;
    cout << "c31=" << c31 << endl;
    cout << "c23*c31=" << c23*-c31 << endl;
    cout << "  ~c31 * ~-c23=" << (~c31)*(~-c23) << endl;
    cout << "~ ~-c31 * ~c23=" << ~((~-c31)*(~c23)) << endl;


    Mat<3,4> m34;
    Mat<3,4,Complex> cm34;


    cm34 = mdc;
    m34 = cm34.real();

    cout << "Mat<3,4,Complex> cm34=" << cm34 << endl;
    cout << "cm34.diag()=" << cm34.diag() << endl;

    cout << "cm34 + cm34=" << cm34+cm34 << endl; //INTERNAL COMPILER ERROR IN Release MODE
    cout << "~cm34 * 1000=" << ~cm34 * 1000. << endl;

    cout << "m34=" << m34 << endl;
    m34 =19.123;
    cout << "after m34=19.123, m34=" << m34 << endl;
 
    const double ddd[] = { 11, 12, 13, 14, 15, 16 }; 
    const complex<float> ccc[] = {  complex<float>(1.,2.),  
                                    complex<float>(3.,4.),
                                    complex<float>(5.,6.),
                                    complex<float>(7.,8.) };
    Vec<2,complex<float>,1> cv2(ccc);
    cout << "cv2 from array=" << cv2 << endl;
    cv2 = Vec<2,complex<float> >(complex<float>(1.,2.), complex<float>(3.,4.));
    cout << "cv2 after assignment=" << cv2 << endl;

    cout << "cv2.real()=" << cv2.real() << " cv2.imag()=" << cv2.imag() << endl;

    Vec<2,negator<complex<float> >,1>& negCv2 = (Vec<2,negator<complex<float> >,1>&)cv2;
    Vec<2,conjugate<float>,1>& conjCv2 = (Vec<2,conjugate<float>,1>&)cv2;
    Vec<2,negator<conjugate<float> >,1>& negConjCv2 = (Vec<2,negator<conjugate<float> >,1>&)cv2;

    

    Vec<2,complex<float> > testMe = cv2;
    cout << "testMe=cv2 (init)=" << testMe << endl;
    testMe = cv2;
    cout << "testMe=cv2 (assign)=" << testMe << endl;


    cout << "(cv2+cv2)/complex<float>(1000,0):" << (cv2 + cv2) / complex<float>(1000,0) << endl; 
    cout << "(cv2+cv2)/1000.f:" << (cv2 + cv2) / 1000.f << endl;
    cout << "(cv2+cv2)/1000.:" << (cv2 + cv2) / 1000. << endl;
    cout << "(cv2+cv2)/1000:" << (cv2 + cv2) / 1000 << endl;

    cout << "negCv2=" << negCv2 << endl;
    cout << "conjCv2=" << conjCv2 << endl;
    cout << "negConjCv2=" << negConjCv2 << endl;
    cout << "cv2+negCv2=" << cv2+negCv2 << endl;

    negConjCv2 = complex<float>(8,9);
    cout << "AFTER negConjCv2 = (8,9):" << endl;
    cout << "  cv2=" << cv2 << endl;
    cout << "  negCv2=" << negCv2 << endl;
    cout << "  conjCv2=" << conjCv2 << endl;
    cout << "  negConjCv2=" << negConjCv2 << endl;

    cout << "cv2:  " << cv2 << endl;
    cout << "cv2T: " << cv2.transpose() << endl; 
    cout << "-cv2: " << -cv2 << endl;
    cout << "~cv2: " << ~cv2 << endl;
    cout << "-~cv2: " << -(~cv2) << endl;
    cout << "~-cv2: " << ~(-cv2) << endl; 
    cout << "~-cv2*10000: " << (~(-cv2))*10000.f << endl;  
        
   (~cv2)[1]=complex<float>(101.1f,202.3f);
    cout << "after ~cv2[1]=(101.1f,202.3f), cv2= " << cv2 << endl;    
    (-(~cv2))[1]=complex<float>(11.1f,22.3f);
    cout << "after -~cv2[1]=(11.1f,22.3f), cv2= " << cv2 << endl; 
        
    Vec<3> dv3(ddd), ddv3(ddd+3);
    dv3[2] = 1000;
    cout << "dv3=" << dv3 << " ddv3=" << ddv3 << endl;
    cout << "100(ddv3-dv3)/1000=" << 100.* (ddv3 - dv3) / 1000. << endl; 

    Vec<3> xxx(dv3); cout << "copy of dv3 xxx=" << xxx << endl;
    Vec<3> yyy(*ddd);cout << "copy of *ddd yyy=" << yyy << endl;
    
    cout << "dv3.norm()=" << dv3.norm() << endl;
    cout << "cv2=" << cv2 << " cv2.norm()=" << cv2.norm() << endl; 
       
    const Vec<2> v2c[] = {Vec<2>(ddd),Vec<2>(ddd+1)};
    Vec<2, Vec<2> > vflt(v2c);
    cout << "vflt 2xvec2=" << vflt << endl;
    cout << "10.*vflt=" << 10.*vflt << endl;
    cout << "vflt*10.=" << vflt*10. << endl;

    int ivals[] = {0x10, 0x20, 0x30, 0x40};
    Vec<4> iv(ivals);
    cout << "iv=" << iv << endl;
    
    Vec<2, Vec<2> > v22;
    v22 = Vec<2>(&ivals[2]);
    cout << "v22=" << v22 << endl;


    // Test dot product
    {
    double d[] = {1,2,3,4,5,6,7,8};


    Vec<2> v1(&d[0]), v2(&d[2]);
    Row<2> r1(&d[4]), r2(&d[6]);
    Vec<2>::TNeg& nv1 = (Vec<2>::TNeg&)v1;

    negator<double> nd(100); cout << endl << "nd=" << nd << endl;
    cout << "nv1=" << nv1 << endl;
    cout << "nv1*nd=" << nv1*nd << endl;
    cout << "nd*nv1=" << nd*nv1 << endl;
    cout << "nv1/nd=" << nv1/nd << endl << endl;


    cout << "v1,v2=" << v1 << v2 << endl;
    cout << "r1,r2=" << r1 << r2 << endl;
    cout << "dot r1*v1 =" << dot(r1,v1) << endl;
    cout << "dot r1*nv1=" << dot(r1,nv1) << endl;
    cout << "r1*v1 =" << r1*v1 << endl;
    cout << "r1*nv1=" << r1*nv1 << endl;
    
    // outer product
    cout << " outer v1*r1=" << v1*r1 << endl;
    cout << " outer nv1*r1=" << nv1*r1 << endl;

    // cross product (2d)
    cout << "cross(v1,v2)=" << cross(v1,v2) << endl;
    cout << "v1 % v2=" << v1 % v2 << endl;
    cout << "cross(r1,v2)=" << cross(r1,v2) << endl;
    cout << "r1 % v2=" << r1 % v2 << endl;
    cout << "cross(v1,r2)=" << cross(v1,r2) << endl;
    cout << "v1 % r2=" << v1 % r2 << endl;
    cout << "cross(r1,r2)=" << cross(r1,r2) << endl;
    cout << "r1 % r2=" << r1 % r2 << endl;

    // do the cross products with 3d routines
    Vec3 v13(v1[0],v1[1],0), v23(v2[0],v2[1],0);
    Row3 r13(r1[0],r1[1],0), r23(r2[0],r2[1],0);
    cout << "cross(v13,v23)=" << cross(v13,v23) << endl;
    cout << "v13 % v23=" << v13 % v23 << endl;
    cout << "cross(r13,v23)=" << cross(r13,v23) << endl;
    cout << "r13 % v23=" << r13 % v23 << endl;
    cout << "cross(v13,r23)=" << cross(v13,r23) << endl;
    cout << "v13 % r23=" << v13 % r23 << endl;
    cout << "cross(r13,r23)=" << cross(r13,r23) << endl;
    cout << "r13 % r23=" << r13 % r23 << endl;

    v13[2]=7;
    cout << "v13=" << v13 << " 2*v13+.001=" << 2*v13+.001 << endl;
    cout << "v13 % (2*v13)=" << v13 % (2*v13) << endl;
    cout << "v13 % (2*v13+.001)=" << v13 % (2*v13+.001) << endl;

    cout << endl;

    // test constructors
    Mat<2,3> mvcols( v1, ~r1, v2 );
    Mat<3,2> mvrows( ~v1, 
                      r1,
                      r2 ); 

    cout << "mvcols=" << mvcols << endl;
    cout << "mvrows=" << mvrows << endl;   

    Vec<3,float> v2f(39.f, 40.f, 50.L);
    cout << "v2f=" << v2f << endl;

    cout << "v2f.drop1(0)=" << v2f.drop1(0) << endl;
    cout << "v2f.drop1(1)=" << v2f.drop1(1) << endl;
    cout << "v2f.drop1(2)=" << v2f.drop1(2) << endl;

    cout << "v2f.append1(3.3f)=" << v2f.append1(3.3f) << endl;
    cout << "v2f.append1((short)1)="   << v2f.append1((short)1)   << endl;
    cout << "v2f.drop1(1).append1((unsigned short)0x10)=" << v2f.drop1(1).append1((unsigned short)0x10) << endl;
    cout << "(~v2f).drop1(1).append1((unsigned short)0x10)=" << (~v2f).drop1(1).append1((unsigned short)0x10) << endl;

    cout << "v2f.insert1(0, 23.f)=" << v2f.insert1(0, 23.f) << endl;
    cout << "v2f.insert1(1, 23.f)=" << v2f.insert1(1, 23.f) << endl;
    cout << "v2f.insert1(2, 23.f)=" << v2f.insert1(2, 23.f) << endl;
    cout << "v2f.insert1(3, 23.f)=" << v2f.insert1(3, 23.f) << endl;
    cout << "v2f.insert1(2, 23.f).drop1(2)=" << v2f.insert1(2, 23.f).drop1(2) << endl;

    cout << "(~v2f).insert1(2, 23.f).drop1(2)=" << (~v2f).insert1(2, 23.f).drop1(2) << endl;

    Mat33 m33( Row3(1,     2,     3),
               Row3(4,     5,     6),
               Row3(.003f, 9.62L, 41.1) );
    cout << "m33=" << m33 << endl;  
    cout << "v13=" << v13 << endl;
    cout << "m33*v13=" << m33*v13 << endl;
    cout << "~m33*v13=" << (~m33)*v13 << endl;
    cout << "~v13*~m33=" << ~v13*(~m33) << endl;
    }

    Mat<4,3> ident43;
    ident43 = 1;
    cout << "ident43=" << ident43 << endl;

    Mat<4,3> negid43 = -ident43; // requires implicit conversion from Mat<negator<Real>>
    cout << "negid43=" << negid43 << endl;

    // Absolute value
    const Vec<3> vorig(1.,-2.,-3.);
    cout << "vorig=" << vorig << " vorig.abs()=" << vorig.abs() << endl;
    const Vec<2, Vec<3> > v2orig(vorig,-vorig);
    cout << "v2orig=" << v2orig << " v2orig.abs()=" << v2orig.abs() << endl;

    Vec<3> nvorig = -vorig;
    Vec<2, Vec<3> > nv2orig = -v2orig;

    Mat<3,2> morig(vorig,-vorig);
    cout << "morig=" << morig << " morig.abs()=" << morig.abs() << endl;

    //Mat<2,2, Vec<3> > m2orig(v2orig, (Vec<2, Vec<3> >)-v2orig);
    //cout << "m2orig=" << m2orig << " m2orig.abs()=" << m2orig.abs() << endl;

    cout << "vorig.getSubVec<2>(1)=" << vorig.getSubVec<2>(1) << endl;

    negid43.updSubMat<2,2>(2,1) = -27.;
    cout << "after negid43.updSubMat<2,2>(2,1) = -27., negid43=" << negid43;

    cout << "negid43[2].getSubRow<2>(1)=" << negid43[2].getSubRow<2>(1) << endl;

    cout << "CHECK DIAGONAL LENGTH FOR RECTANGULAR MATRICES" << endl;
    Mat<3,2, Row3, 1, 2> H;
    cout << "H[" << H.nrow() << "," << H.ncol() << "]" << endl;
    cout << "H.diag()[" << H.diag().nrow() << "," << H.diag().ncol() << "]" << endl;
    Mat<3,2, Row3, 1, 2>::TransposeType Ht;
    cout << "Ht[" << Ht.nrow() << "," << Ht.ncol() << "]" << endl;
    cout << "Ht.diag()[" << Ht.diag().nrow() << "," << Ht.diag().ncol() << "]" << endl;
}

void testMatInverse() {
    Matrix m = Test::randMatrix(20,20);
    Matrix mi = m.invert();
    Matrix id(20,20); id=1; // identity
    SimTK_TEST_EQ_SIZE(m*mi, id, 20);
}


int main() {
    SimTK_START_TEST("MatVecTest");

        SimTK_SUBTEST(testSums);
        SimTK_SUBTEST(testNegator);
        SimTK_SUBTEST(testElementwiseOps);
        SimTK_SUBTEST(testMiscellaneous);
        SimTK_SUBTEST(testMatInverse);

    SimTK_END_TEST();
}



