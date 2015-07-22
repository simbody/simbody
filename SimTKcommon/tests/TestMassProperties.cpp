/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2009-12 Stanford University and the Authors.        *
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
using std::cout;
using std::endl;

using namespace SimTK;

// We use cross product and cross product matrices extensively in shifting
// mass properties around.

void testCrossProduct() {
    // These are exactly representable in base 2.
    Vec3 w(1.25, 3, -2.5), v(-2.75, 2.125, 5);
    Vec<3,float> wf(1.25, 3, -2.5), vf(-2.75, 2.125, 5);
    Vec3 wxv = w % v, vxw = v % w;
    Vec<3,float> wxvf = wf % vf, vxwf = vf % wf;
    SimTK_TEST(wxv == Vec3(20.3125, 0.625, 10.90625));
    SimTK_TEST(wxvf == (Vec<3,float>(20.3125, 0.625, 10.90625)));
    SimTK_TEST(vxw == -wxv);
    SimTK_TEST(vxwf == -wxvf);

    // Cross product involving one or more rows returns a row.
    SimTK_TEST((~w)%v == ~wxv);
    SimTK_TEST(w % ~v == ~wxv);
    SimTK_TEST((~w) % (~v) == ~wxv);

    SimTK_TEST( crossMat(w) == Mat33( 0,   2.5,  3,
                                     -2.5, 0,   -1.25,
                                     -3,   1.25, 0) );

    SimTK_TEST( crossMatSq(w) == ~crossMat(w)*crossMat(w) );

    const Mat33 full33 = Test::randMat33();
    const SymMat33 sym33 = Test::randSymMat33();
    const Mat33 fsym33(sym33);
    SimTK_TEST(fsym33 == sym33);

    SimTK_TEST_EQ(v % full33, Mat33( v%full33(0), v%full33(1), v%full33(2) ));
    SimTK_TEST_EQ(v % full33, crossMat(v)*full33);
    SimTK_TEST_EQ(v % sym33,  Mat33( v%fsym33(0), v%fsym33(1), v%fsym33(2) ));
    SimTK_TEST_EQ(v % sym33, crossMat(v)*fsym33);

    // m % v == m * crossMat(v) == -~(v % ~m)
    SimTK_TEST_EQ(full33 % v, full33 * crossMat(v));
    SimTK_TEST_EQ(full33 % v, -~(v % ~full33));
    SimTK_TEST_EQ(full33 % v, ~(-v % ~full33));
    SimTK_TEST_EQ(sym33 % v, fsym33 * crossMat(v));
    SimTK_TEST_EQ(sym33 % v, -~(v % sym33));
    SimTK_TEST_EQ(sym33 % v, ~(-v % sym33));

    // gratuitous det() test
    SimTK_TEST_EQ(det(sym33), det(fsym33));

    // 2D

    // These are exactly representable in base 2.
    Vec2 w2(1.25, 3), v2(-2.75, 2.125);
    Vec<2,float> w2f(1.25, 3), v2f(-2.75, 2.125);
    Real wxv2 = w2 % v2, vxw2 = v2 % w2;
    float wxv2f = w2f % v2f, vxw2f = v2f % w2f;
    // 2D cross product is the z component of the 3D cross product where
    // the z coordinates of the two vectors are zero (i.e., they lie in the x-y plane).
    SimTK_TEST(wxv2 ==  (Vec3(w2[0],w2[1],0)           % Vec3(v2[0],v2[1],0))[2]);
    SimTK_TEST(wxv2f == (Vec<3,float>(w2f[0],w2f[1],0) % Vec<3,float>(v2f[0],v2f[1],0))[2]);
    SimTK_TEST(vxw2 == -wxv2);
    SimTK_TEST(vxw2f == -wxv2f);

    // Cross product involving one or more rows returns same result in 2d.
    SimTK_TEST((~w2)%v2 == wxv2);
    SimTK_TEST(w2 % ~v2 == wxv2);
    SimTK_TEST((~w2) % (~v2) == wxv2);

    SimTK_TEST( crossMat(w2) == Row2( -w2[1], w2[0] ) );
    SimTK_TEST( crossMat(w2)*v2 == w2 % v2 );

    const Mat22 full22 = Test::randMat<2,2>();
    const SymMat22 sym22 = Test::randSymMat<2>();
    const Mat22 fsym22(sym22);
    SimTK_TEST(fsym22 == sym22);

    /* These don't exist in 2D
    SimTK_TEST_EQ(v2 % full22, Row2( v2%full22(0), v2%full22(1) ));
    SimTK_TEST_EQ(v2 % full22, crossMat(v2)*full22);
    SimTK_TEST_EQ(v2 % sym22,  Row2( v2%fsym22(0), v2%fsym22(1) ));
    SimTK_TEST_EQ(v2 % sym22, crossMat(v2)*fsym22);

    // m % v == m * crossMat(v) == -~(v % ~m)
    SimTK_TEST_EQ(full22 % v2, full22 * crossMat(v2));
    SimTK_TEST_EQ(full22 % v2, -~(v2 % ~full22));
    SimTK_TEST_EQ(full22 % v2, ~(-v2 % ~full22));
    SimTK_TEST_EQ(sym22 % v2, fsym22 * crossMat(v2));
    SimTK_TEST_EQ(sym22 % v2, -~(v2 % sym22));
    SimTK_TEST_EQ(sym22 % v2, ~(-v2 % sym22));
    */

    // gratuitous det() test
    SimTK_TEST_EQ(det(sym22), det(fsym22));
}

void testInertia() {
    const Real mass = std::abs(Test::randReal());
    const UnitInertia_<Real> G( Vec3(1,2,2.5),       // moments
                                Vec3(0.1,0.2,0.3) ); // products
    const Real Gtrace = 1+2+2.5;

    SimTK_TEST(G.trace() == Gtrace); // should be exact because .5 is power of 2

    const Inertia_<Real>  I = mass*G;
    SymMat33 sI = I.asSymMat33();
    Mat33    mI = I.toMat33();

    SimTK_TEST_EQ( sI, mass*SymMat33(1,
                                     0.1, 2,
                                     0.2, 0.3, 2.5) );

    SimTK_TEST(mI.isExactlySymmetric());
    SimTK_TEST(mI == Mat33(sI));
    SimTK_TEST(sI == SymMat33(mI));

    SimTK_TEST_EQ( I.trace(), mass*Gtrace );

    // Test Inertia*scalar
    const Real s = std::abs(Test::randReal()) + 1;
    SimTK_TEST_EQ( (I*s).toMat33(), I.toMat33()*s );
    SimTK_TEST_EQ( (s*I).toMat33(), I.toMat33()*s );
    SimTK_TEST_EQ( (G*s).toMat33(), G.toMat33()*s );
    SimTK_TEST_EQ( (s*G).toMat33(), G.toMat33()*s );
    SimTK_TEST_EQ( (I*s).asSymMat33(), I.asSymMat33()*s );
    SimTK_TEST_EQ( (s*I).asSymMat33(), I.asSymMat33()*s );
    SimTK_TEST_EQ( (G*s).asSymMat33(), G.asSymMat33()*s );
    SimTK_TEST_EQ( (s*G).asSymMat33(), G.asSymMat33()*s );

    // Test Inertia*vec (I*w)
    const Vec3 w = Test::randVec3();
    SimTK_TEST_EQ( I*w, I.toMat33()*w );
    SimTK_TEST_EQ( G*w, G.toMat33()*w );
    SimTK_TEST_EQ( I*w, I.asSymMat33()*w );
    SimTK_TEST_EQ( G*w, G.asSymMat33()*w );

    // Test inertia rotation

    const Rotation R = Test::randRotation();
    Inertia_<Real> mIR = Inertia_<Real>(~R*I.toMat33()*R);

    SimTK_TEST_EQ(mIR.asSymMat33(), I.reexpress(R).asSymMat33());
    SimTK_TEST_EQ(I.asSymMat33(),   mIR.reexpress(~R).asSymMat33());

    Inertia_<Real> J=I;
    J.reexpressInPlace(R);
    SimTK_TEST_EQ(J.asSymMat33(), mIR.asSymMat33());

    // Test inertia shifting

    // Calculate the inertia of a point mass with the same mass
    // we used above.
    const Vec3     pLoc = Test::randVec3();
    const SymMat33 psG = crossMatSq(pLoc); // unit inertia
    const SymMat33 psI = mass*psG; // inertia
    const UnitInertia_<Real> pG(psG);
    const Inertia_<Real>     pI(psI);

    // Assuming I and G are central, shifting them to pLoc should be
    // the same as adding the point inertias above.
    SimTK_TEST_EQ( G.shiftFromCentroid(pLoc), G + pG );
    SimTK_TEST_EQ( I.shiftFromMassCenter(pLoc, mass), I + pI );

    // Now try in place shifts and shifting back.
    UnitInertia_<Real> Gshft(G); Gshft.shiftFromCentroidInPlace(pLoc);
    Inertia_<Real>     Ishft(I); Ishft.shiftFromMassCenterInPlace(pLoc, mass);
    SimTK_TEST_EQ(Gshft, G+pG);
    SimTK_TEST_EQ(Ishft, I+pI);

    SimTK_TEST_EQ(Gshft.shiftToCentroid(pLoc), G);
    SimTK_TEST_EQ(Ishft.shiftToMassCenter(pLoc, mass), I);

    Gshft.shiftToCentroidInPlace(pLoc);
    Ishft.shiftToMassCenterInPlace(pLoc, mass);
    SimTK_TEST_EQ(Gshft, G);
    SimTK_TEST_EQ(Ishft, I);

    // Check that we catch bad inertias in debug mode
#ifndef NDEBUG
    Inertia bad1;
    SimTK_TEST_MUST_THROW(bad1 = Inertia(1,2,NaN));
    SimTK_TEST_MUST_THROW(bad1 = Inertia(1,2,-.000003)); // negative diag
    SimTK_TEST_MUST_THROW(bad1 = Inertia(5,1,2)); // triangle inequality violated

    // A little slop should be allowed for the triangle inequality.
    Real tooMuchSlop = 1e-3;
    Real okSlop = SignificantReal;
    SimTK_TEST_MUST_THROW(bad1 = Inertia(1, 2+tooMuchSlop, 1));
    bad1 = Inertia(1, 2+okSlop, 1);


#endif

}

// Calculate the lower half of vx*F where vx is the cross product matrix
// of v and F is a full 3x3 matrix. This result would normally be a full
// 3x3 but for the uses below we know we're only going to need the diagonal
// and lower triangle so we can save some flops by working this out by hand.
// The method is templatized so that it will work on a transposed matrix
// as efficiently as an untransposed one. (18 flops)
template <class P, int CS, int RS>
static inline SymMat<3,P>
halfCross(const Vec<3,P>& v, const Mat<3,3,P,CS,RS>& F) {
    return SymMat<3,P>
      ( v[1]*F(2,0)-v[2]*F(1,0),
        v[2]*F(0,0)-v[0]*F(2,0), v[2]*F(0,1)-v[0]*F(2,1),
        v[0]*F(1,0)-v[1]*F(0,0), v[0]*F(1,1)-v[1]*F(0,1), v[0]*F(1,2)-v[1]*F(0,2) );
}

// Calculate the lower half of G*vx where G is a full 3x3 matrix and vx
// is the cross product matrix of v. See comment above for details.
// (18 flops)
template <class P, int CS, int RS>
static inline SymMat<3,P>
halfCross(const Mat<3,3,P,CS,RS>& G, const Vec<3,P>& v) {
    return SymMat<3,P>
      ( v[2]*G(0,1)-v[1]*G(0,2),
        v[2]*G(1,1)-v[1]*G(1,2), v[0]*G(1,2)-v[2]*G(1,0),
        v[2]*G(2,1)-v[1]*G(2,2), v[0]*G(2,2)-v[2]*G(2,0), v[1]*G(2,0)-v[0]*G(2,1) );
}

// This method computes the lower half of the difference vx*F-G*vx using
// the same methods as above, but done together in order to pull out the
// common v terms. This is 33 flops, down from 42 if you call the two
// methods above and add them.
template <class P, int CS1, int RS1, int CS2, int RS2>
static inline SymMat<3,P>
halfCrossDiff(const Vec<3,P>& v, const Mat<3,3,P,CS1,RS1>& F, const Mat<3,3,P,CS2,RS2>& G) {
    return SymMat<3,P>
      ( v[1]*(F(2,0)+G(0,2)) - v[2]*(F(1,0)+G(0,1)),
        v[2]*(F(0,0)-G(1,1)) - v[0]*F(2,0) + v[1]*G(1,2),
                v[2]*(F(0,1)+G(1,0)) - v[0]*(F(2,1)+G(1,2)),
        v[0]*F(1,0) - v[2]*G(2,1) - v[1]*(F(0,0)-G(2,2)),
                v[0]*(F(1,1)-G(2,2)) - v[1]*F(0,1) + v[2]*G(2,0),
                        v[0]*(F(1,2)+G(2,1)) - v[1]*(F(0,2)+G(2,0)) );
}

void testHalfCross() {
    const Vec3 v = Test::randVec3();
    const Mat33 F = Test::randMat33();
    const Mat33 G = Test::randMat33();
    const Mat33 vx = crossMat(v);
    const Mat33 vxF = v % F;
    const Mat33 Gxv = G % v;

    SimTK_TEST_EQ(vxF, vx * F);
    SimTK_TEST_EQ(Gxv, G * vx);

    const SymMat33 hvxF( vxF(0,0),
                         vxF(1,0), vxF(1,1),
                         vxF(2,0), vxF(2,1), vxF(2,2) );
    SimTK_TEST_EQ(halfCross(v, F), hvxF);

    const SymMat33 hGxv( Gxv(0,0),
                         Gxv(1,0), Gxv(1,1),
                         Gxv(2,0), Gxv(2,1), Gxv(2,2) );
    SimTK_TEST_EQ(halfCross(G, v), hGxv);

    const SymMat33 hdiff = hvxF-hGxv;
    SimTK_TEST_EQ(halfCrossDiff(v, F, G), hdiff);

}

void testSpatialInertia() {
    const Real     mass = 1.125;
    const Vec3     com(.1, .2, .25);
    const UnitInertia gyration(1.8, 1.9, 2.1, .01, .03, .02);

    SpatialInertia si(mass,com,gyration);

    const SpatialMat msi = si.toSpatialMat();
    SimTK_TEST( msi(0,0) == mass*gyration.toMat33() );
    SimTK_TEST( msi(0,1) == mass*crossMat(com) );
    SimTK_TEST( msi(1,0) == mass*~crossMat(com) );
    SimTK_TEST( msi(1,1) == mass*Mat33(1) );

    // Using phi shifts by the negative of the shift vector.
    const Vec3 shiftVec(1,2,3);
    const SpatialMat shiftMat(Mat33(1), crossMat(-shiftVec),
                              Mat33(0),    Mat33(1));
    const SpatialMat msiShiftedManually = shiftMat*msi*~shiftMat;

    const PhiMatrix phi(-shiftVec);
    const SpatialMat msiShiftedByPhi = phi*msi*~phi;

    SimTK_TEST_EQ(msiShiftedByPhi, msiShiftedManually);

    // The SpatialInertia shift() method behaves correctly; the Articulated
    // one currently has the reverse sense, see below.
    const SpatialInertia shiftSi = si.shift(shiftVec);

    SimTK_TEST_EQ(shiftSi.toSpatialMat(), msiShiftedManually);

    SimTK_TEST_EQ(shiftSi.shift(-shiftVec).toSpatialMat(),
                  si.toSpatialMat());

    si.shiftInPlace(shiftVec);
    SimTK_TEST_EQ(si.toSpatialMat(), msiShiftedManually);

    SimTK_TEST_EQ(si.shiftInPlace(-shiftVec).toSpatialMat(), msi);
}

void testArticulatedInertia() {
    const SymMat33 mass = Test::randSymMat33();
    const Mat33    massMoment = Test::randMat33();
    const SymMat33 inertia = crossMatSq(Test::randVec3());

    ArticulatedInertia abi(mass,massMoment,inertia);

    const SpatialMat mabi = abi.toSpatialMat();
    SimTK_TEST( mabi(0,0) == inertia );
    SimTK_TEST( mabi(0,1) == massMoment );
    SimTK_TEST( mabi(1,0) == ~massMoment );
    SimTK_TEST( mabi(1,1) == mass );

    // Using phi shifts by the negative of the shift vector.
    const Vec3 shiftVec(1,2,3);
    const SpatialMat shiftMat(Mat33(1), crossMat(-shiftVec),
                              Mat33(0),    Mat33(1));
    const SpatialMat mabiShiftedManually = shiftMat*mabi*~shiftMat;

    const PhiMatrix phi(-shiftVec);
    const SpatialMat mabiShiftedByPhi = phi*mabi*~phi;

    SimTK_TEST_EQ(mabiShiftedByPhi, mabiShiftedManually);

    // Unfortunately the abi shift() method is also defined to shift by the
    // negative of the shift vector; that's an API bug.
    const Vec3 negShiftVec(-shiftVec);
    const ArticulatedInertia shiftAbi = abi.shift(negShiftVec);

    SimTK_TEST_EQ(shiftAbi.toSpatialMat(), mabiShiftedManually);

    SimTK_TEST_EQ(shiftAbi.shift(-negShiftVec).toSpatialMat(),
                  abi.toSpatialMat());

    abi.shiftInPlace(negShiftVec);
    SimTK_TEST_EQ(abi.toSpatialMat(), mabiShiftedManually);

    SimTK_TEST_EQ(abi.shiftInPlace(-negShiftVec).toSpatialMat(), mabi);
}

void testManualABIShift(const ArticulatedInertia& abi, const Array_<Vec3>& shifts, Real& out) {
    SpatialMat mabi = abi.toSpatialMat();

    SpatialMat shiftMat(Mat33(1));
    for (unsigned i=0; i<shifts.size(); ++i) {
        shiftMat(0,1) = crossMat(shifts[i]);
        mabi = shiftMat * mabi * ~shiftMat;
    }
    out = mabi(1,1)(2,2);
}

void testPhiABIShift(const ArticulatedInertia& abi, const Array_<Vec3>& shifts, Real& out) {
    SpatialMat mabi = abi.toSpatialMat();

    for (unsigned i=0; i<shifts.size(); ++i) {
        const PhiMatrix phi(shifts[i]);
        mabi = phi * mabi * ~phi;
    }
    out = mabi(1,1)(2,2);
}

void testFastABIShift(const ArticulatedInertia& abi_in, const Array_<Vec3>& shifts, Real& out) {
    ArticulatedInertia abi(abi_in);
    for (unsigned i=0; i<shifts.size(); ++i)
        abi = abi.shift(shifts[i]);
    out = abi.getMass()(2,2);
}

Array_<Vec3> shifts;
#ifdef NDEBUG
    #define NSHIFTS 2000000
#else
    #define NSHIFTS 100000
#endif

int main() {
    SimTK_START_TEST("TestMassProperties");

        SimTK_SUBTEST(testCrossProduct);
        SimTK_SUBTEST(testInertia);
        SimTK_SUBTEST(testHalfCross);
        SimTK_SUBTEST(testSpatialInertia);
        SimTK_SUBTEST(testArticulatedInertia);

        // Speed tests.
        shifts.resize(NSHIFTS);
        for (int i=0; i<NSHIFTS; ++i) shifts[i] = Test::randVec3();
        const SymMat33 mass = Test::randSymMat33();
        const Mat33    massMoment = Test::randMat33();
        const SymMat33 inertia = crossMatSq(Test::randVec3());
        const ArticulatedInertia abi(mass,massMoment,inertia);
        Real out1, out2, out3;
        SimTK_SUBTEST3(testManualABIShift, abi, shifts, out1);
        SimTK_SUBTEST3(testPhiABIShift, abi, shifts, out2);
        SimTK_SUBTEST3(testFastABIShift, abi, shifts, out3);
        SimTK_TEST_EQ(out1, out2); SimTK_TEST_EQ(out2, out3);

    SimTK_END_TEST();
}

