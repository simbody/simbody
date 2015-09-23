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

#include "SimTKcommon.h"
#include "SimTKcommon/Testing.h"

using namespace SimTK;

#include <iostream>
#include <iomanip>
using std::cout;
using std::endl;
using std::setw;

void testPhiMatrix() {
    const Vec3 p = Test::randVec3();
    PhiMatrix phi(p);
    SpatialVec v1(Test::randSpatialVec());
    SpatialMat m1(Test::randSpatialMat());

    SimTK_TEST_EQ( phi*v1, SpatialVec(v1[0] + p%v1[1], v1[1]));
    SimTK_TEST_EQ(~phi*v1, SpatialVec(v1[0], v1[1] - p%v1[0]));

    SimTK_TEST_EQ(phi*v1, phi.toSpatialMat()*v1);
    SimTK_TEST_EQ(phi*m1, phi.toSpatialMat()*m1);

    SimTK_TEST_EQ(~phi*v1, (~phi).toSpatialMat()*v1);
    SimTK_TEST_EQ(m1*~phi, m1*(~phi).toSpatialMat());
}

// TODO: this isn't a real regression test but it does catch
// compilation problems.
void testMiscSpatialAlgebra() {
    SpatialVec v; SpatialRow r; SpatialMat m;
    Mat33 m33(1,2,3,
              4,5,6,
              7,8,9);
    Mat22 m22(10,20,
              30,40);
    cout << endl << "TEST: uninitialized should be NaN in Debug, random crap in Release" << endl;
    cout << "rawVec=" << v << " rawRow=" << r << " rawMat=" << m;

    cout << "SpatialMat::NRows=" << SpatialMat::NRows
         << " SpatialMat::ArgDepth=" << CNT<SpatialMat>::ArgDepth
         << " CNT<SpatialMat>::ArgDepth=" << SpatialMat::ArgDepth
         << " SpatialRow::ArgDepth=" << SpatialRow::ArgDepth
         << " SpatialVec::ArgDepth=" << SpatialVec::ArgDepth
         << " Mat33::ArgDepth=" << Mat33::ArgDepth
         << " CNT<Real>::ArgDepth=" << CNT<Real>::ArgDepth
         << endl;

    cout << endl << "TEST: set to 1; element wise for col/row; diagonal for mat" << endl;
    v=1;r=1;m=1;
    cout << "v=1:" << v << " r=1:" << r << " m=1:" << m;

    SpatialVec sv( Vec3(1,2,-3), Vec3(.1,-.2,.3) );
    Vec6& sv6 = Vec6::updAs(&sv[0][0]);
    cout << "sv (" << &sv << ")=" << sv << endl;
    cout << "sv6(" << &sv6 << ")=" << sv6 << endl;
    cout << "sv.normalize() =" << sv.normalize() << endl;
    cout << "sv6.normalize()=" << sv6.normalize() << endl;

    SpatialVec::TNeg& nsv = sv.updNegate();
    Vec6::TNeg& nsv6 = Vec6::TNeg::updAs(&nsv[0][0]);
    cout << "nsv (" << &nsv << ")=" << nsv << endl;
    cout << "nsv6(" << &nsv6 << ")=" << nsv6 << endl;
    cout << "nsv.normalize() =" << nsv.normalize() << endl;
    cout << "nsv6.normalize()=" << nsv6.normalize() << endl;

    sv = sv.normalize();
    cout << "after sv=sv.normalize(), sv=" << sv << endl;
    cout << "now nsv+1=" << nsv+1 << ", (nsv+1).normalize()=" << (nsv+1).normalize() << endl;
    nsv = (nsv+1).normalize();
    cout << "after nsv=(nsv+1).normalize(), nsv=" << nsv << endl;

    CNT<Vec3>::Result<double>::Mul xxx;

    cout << endl << "TEST: v+v=" << v+v;
    cout << endl << "TEST: v*3.=" << v*3.;
    //cout << endl << "TEST: 3./v=" << 3./v;


    cout << endl << "TEST: 3.f*v=" << 3.f*v;
    cout << endl << "TEST: 3*v*2=" << 3*v*2;
    cout << endl << "TEST: v/negator<Real>(3.)=" << v/negator<Real>(3.);
    cout << endl << "TEST: v/complex<Real>(2.,3.)=" << v/complex<Real>(2.,3.);
    cout << endl << "TEST: v/conjugate<Real>(2.,3.)=" << v/conjugate<Real>(2.,3.);
    const conjugate<Real> cr23(2.,3.);
    cout << endl << "TEST: v*negator<conjugate<Real>(2.,3.)>="
         << v*negator<conjugate<Real> >::recast(cr23);
    cout << endl << "TEST: v/negator<conjugate<Real>(2.,3.)>="
         << v/negator<conjugate<Real> >::recast(cr23);
    //cout << endl << "TEST: v-v=" << v-v;
    //cout << endl << "TEST: r+r=" << r+r;
    //cout << endl << "TEST: r-r=" << r-r;
    //cout << endl << "TEST: m+m=" << m+m;
    //cout << endl << "TEST: m-m=" << m-m;

    // CONFORMING MULTIPLY

    cout << endl << "TEST: multiply by identity matrix" << endl;
    cout << "m*v:" << m*v << endl;
    cout << " r*m:" << r*m << endl;
    cout << " m*m:" << m*m;
    cout << " m*m22:" << m*m22;

    cout << endl << "TEST: scalar multiply, left and right" << endl;
    cout << "v*9=" << v*9 << " 2*v=" << 2*v << endl;
    cout << "r*9=" << r*9 << " 2*r=" << 2*r << endl;
    cout << "m*9=" << m*9 << " 2*m=" << 2*m << endl;

    cout << endl << "TEST: dot product" << endl;
    cout << "r*v=" << r*v << " ~v*~r=" << ~v*~r << endl;

    cout << endl << "TEST: outer product" << endl;
    cout << "v*r=" << v*r << endl;
    cout << " ~r*~v=" << ~r*~v << endl;

    // NONCONFORMING MULTIPLY

    cout << endl << "TEST: m*m33=" << m * m33 << "   ~m*m33=" << ~m*m33;
    cout << endl << "TEST: m33*m=" << m33 * m;
    cout << endl << "TEST: r*m33=" << r * m33 << "   ~v*m33=" << ~v*m33;
    cout << endl << "TEST: m33*v=" << m33 * v << endl << endl;

    Row<2, Mat33> m1233(Mat33(2), Mat33(3));
    Vec3 v310(10);
    SpatialRow sr = ~v310 * m1233;
    cout << "v310=" << v310 << " m1233=" << m1233;
    cout << "~v310*m1233=" << sr << endl;

    cout << endl << "TEST: element multiply by non-scalar" << endl;
    Row<2> rnew = r * Vec3(1,2,3);
    cout << "r*vec(1,2,3)=" << rnew << endl;

    // SYMMAT
    SymMat<3> sy3 = SymMat<3>().setFromLower(Mat<3,3>(2, 99, 99,
                                                      3,  4, 99,
                                                      5,  6,  7));
    SymMat<3> sy3d(-5);
    cout << "sy3=" << sy3 << "  sy3d=" << sy3d;
    cout << "sy3+sy3d=" << sy3+sy3d;

    // TRICKY CASE: 1x1 matrix, 1-row, 1-vec
    Mat<1, 2, Row3> m11(Row3(1,2,3), Row3(4,5,6));
    Row<1>   r1(-.1);
    Vec<2, Real>   v1(10,20);
    cout << "r1=" << r1 << " v1=" << v1 << " m11=" << m11;
    cout << "m11*v1=" << m11*v1 << endl;
    cout << "r1*m11=" << r1*m11 << endl;;
    //cout << "r2*v1=" << r1*v1 << endl;
    //cout << "v1*r1=" << v1*r1 << endl;;


    MassProperties mprops(23, Vec3(1,2,3), UnitInertia::brick(.1,.2,.3));
    cout << "MassProperties: " << mprops;
    cout << "MassProperties.toSpatialMat: " << mprops.toSpatialMat();
    cout << "MassProperties.toMat66: " << mprops.toMat66();
}


int main()
{
    SimTK_START_TEST("SpatialAlgebraTest");
        SimTK_SUBTEST(testPhiMatrix);
        SimTK_SUBTEST(testMiscSpatialAlgebra);
    SimTK_END_TEST();
}
