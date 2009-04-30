/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-9 Stanford University and the Authors.         *
 * Authors: Peter Eastman                                                     *
 * Contributors: Michael Sherman                                              *
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
#include "SimTKcommon/Testing.h"

#include <iostream>

using std::cout;
using std::endl;
using namespace SimTK;
using namespace std;

// A matrix operation of size N can be expected to achieve an
// accuracy of about N*tol where tol is the expected accuracy
// of a scalar operation.
static double NTol(int N) {return N*Test::defTol<Real>();}

template <int N>
void testInverse() {
    Mat<N,N> identity(1);
    Mat<N,N> mat = Test::randMat<N,N>();
    SimTK_TEST_NUMEQ_TOL(mat*mat.invert(), identity, NTol(N));
    mat = ~mat;
    SimTK_TEST_NUMEQ_TOL(mat*mat.invert(), identity, NTol(N));
    SimTK_TEST_NUMEQ_TOL((-mat)*(-mat).invert(), identity, NTol(N));
}

void testDotProducts() {
    Vec3 v1(1, 2, 3);
    Vec3 v2(-1, -2, -3);
    Row3 r1(0.1, 0.2, 0.3);
    Row3 r2(-0.1, -0.2, -0.3);
    SimTK_TEST_NUMEQ(dot(v1, v2), -14.0);
    SimTK_TEST_NUMEQ(dot(r1, r2), -0.14);
    SimTK_TEST_NUMEQ(dot(v1, r2), -1.4);
    SimTK_TEST_NUMEQ(dot(r1, v2), -1.4);
    SimTK_TEST_NUMEQ(r1*v2, -1.4);
    SpatialVec sv(Vec3(1, 2, 3), Vec3(4, 5, 6));
    SpatialRow sr(Row3(1, 2, 3), Row3(4, 5, 6));
    SimTK_TEST_NUMEQ(sr*sv, 91.0);
}

void testCrossProducts() {
    Vec3 w = Test::randVec3();
    Vec3 v = Test::randVec3();

    SimTK_TEST_NUMEQ( w%v, Vec3(w[1]*v[2]-w[2]*v[1],
                                w[2]*v[0]-w[0]*v[2],
                                w[0]*v[1]-w[1]*v[0]) );

    SimTK_TEST_NUMEQ(  w %  v,  cross(  w,  v));
    SimTK_TEST_NUMEQ( ~w % ~v, -cross( ~v, ~w));

    SimTK_TEST_NUMEQ( crossMat(w)*v, w % v );
    SimTK_TEST_NUMEQ( crossMat(~w), crossMat(w) );
    SimTK_TEST_NUMEQ( crossMat(-w), ~crossMat(w) );
    SimTK_TEST_NUMEQ( crossMatSq(w)*v, -w % (w%v) );
}

// Individually test 2x2, 3x3, and 4x4 because the
// smaller sizes may have specialized inline operators.
void testSymMat() {
    // 2x2
    const Vec3   a = Test::randVec3();
    const Vec<2> v = Test::randVec<2>();
    SymMat<2> sm( a[0],
                  a[1], a[2] );
    Mat<2,2>   m( a[0], a[1],
                  a[1], a[2] );

    SimTK_TEST_NUMEQ( (Mat<2,2>(sm)), m );
    SimTK_TEST_NUMEQ( sm, SymMat<2>(m) );

    SimTK_TEST_NUMEQ( sm*v, m*v );
    SimTK_TEST_NUMEQ( ~v*sm, ~v*m );

    // 3x3
    const Vec<6> a3 = Test::randVec<6>();
    const Vec<3> v3 = Test::randVec<3>();
    SymMat<3> sm3( a3[0],
                   a3[1], a3[2],
                   a3[3], a3[4], a3[5]);
    Mat<3,3>   m3( a3[0], a3[1], a3[3],
                   a3[1], a3[2], a3[4],
                   a3[3], a3[4], a3[5]);

    SimTK_TEST_NUMEQ( (Mat<3,3>(sm3)), m3 );
    SimTK_TEST_NUMEQ( sm3, SymMat<3>(m3) );

    SimTK_TEST_NUMEQ( sm3*v3, m3*v3 );
    SimTK_TEST_NUMEQ( ~v3*sm3, ~v3*m3 );

    // 4x4 (hopefully the general case)
    const Vec<10> a4 = Test::randVec<10>();
    const Vec<4>  v4 = Test::randVec<4>();
    SymMat<4> sm4( a4[0],
                   a4[1], a4[2],
                   a4[3], a4[4], a4[5],
                   a4[6], a4[7], a4[8], a4[9]);
    Mat<4,4>   m4( a4[0], a4[1], a4[3], a4[6],
                   a4[1], a4[2], a4[4], a4[7],
                   a4[3], a4[4], a4[5], a4[8],
                   a4[6], a4[7], a4[8], a4[9]);

    SimTK_TEST_NUMEQ( (Mat<4,4>(sm4)), m4 );
    SimTK_TEST_NUMEQ( sm4, SymMat<4>(m4) );

    SimTK_TEST_NUMEQ( sm4*v4, m4*v4 );
    SimTK_TEST_NUMEQ( ~v4*sm4, ~v4*m4 );


    // Complex is tricky for symmetric (really Hermitian) matrices because
    // the diagonals must be real and the corresponding off-diagonals are
    // complex conjugate pairs, NOT the same value even though the off
    // diagonal data is only stored once.
    const Vec<3,Complex> ac( Test::randComplex(), Test::randComplex(), Test::randComplex() );
    const Vec<2,Complex> vc( Test::randComplex(), Test::randComplex() );
    SymMat<2,Complex> smc( ac[0],
                           ac[1], ac[2] );
    Mat<2,2,Complex>   mc( ac[0].real(), std::conj(ac[1]),
                           ac[1],         ac[2].real() );

    // This constructor has to figure out how to generate a conjugate 
    // element for the upper right in the full Mat.
    Mat<2,2,Complex> sm2mc(smc);
    SimTK_TEST_NUMEQ( sm2mc, mc );
    SimTK_TEST_NUMEQ( smc, (SymMat<2,Complex>(mc)) );

    SimTK_TEST_NUMEQ( smc*vc, mc*vc );
    SimTK_TEST_NUMEQ( ~vc*smc, ~vc*mc );

    SimTK_TEST_NUMEQ( ~smc*vc, ~mc*vc );
    SimTK_TEST_NUMEQ( ~smc*vc, smc*vc );    
}

int main() {
    SimTK_START_TEST("TestSmallMatrix");
        SimTK_SUBTEST(testSymMat);
        SimTK_SUBTEST(testInverse<1>);
        SimTK_SUBTEST(testInverse<2>);
        SimTK_SUBTEST(testInverse<3>);
        SimTK_SUBTEST(testInverse<5>);
        SimTK_SUBTEST(testInverse<10>);
        SimTK_SUBTEST(testDotProducts);
        SimTK_SUBTEST(testCrossProducts);
    SimTK_END_TEST();
}
