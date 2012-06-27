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
 * Contributors: Paul Mitiguy                                                 *
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

/**@file
 * Adhoc tests for rotation-related classes.
 */

#include "SimTKcommon.h"

#include <string>
#include <typeinfo>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
using std::cout;
using std::endl;

using namespace SimTK;

static Rotation rotate1(int i, Real a);

float  feps = std::numeric_limits<float>::epsilon();
double deps = std::numeric_limits<double>::epsilon();

// What is the smallest a for which acos(cos(a))=a to machine precision?
double  ferr(float angle)  {return double((double)std::acos(std::cos(angle))-(double)angle)/(std::fabs((double)angle)+deps);}
double derr(double angle) {return (std::acos(std::cos(angle))-angle)/(std::fabs(angle)+deps);}
double err(double a, double aest) {return std::fabs(aest-a)/std::max(std::fabs(a), deps);}
void cosTest() {

    float sfeps = std::sqrt(feps);
    double sdeps = std::sqrt(deps);

    float csf = std::cos(sfeps);
    double csd = std::cos(sdeps);

    cout << "feps=" << feps << " sfeps=" << sfeps << " cos(sf)-1" << csf-1 << endl;
    cout << "deps=" << deps << " sdeps=" << sdeps << " cos(sd)-1" << csd-1 << endl;

    double pi2 = std::acos(0.);
    double maxs=0, maxc=0, maxt=0;
    for (int i=-10000001; i<=10000001; ++i) {
        double a = pi2 - 1.03e-9*double(i)*pi2 + 0.237e-9;
        double s = std::sin(a), c = std::cos(a);
        maxs = std::max(maxs, err(a, std::asin(s)));
        maxc = std::max(maxc, err(a, std::acos(c)));
        maxt = std::max(maxt, err(a, std::atan2(s,c)));
    }
    printf("near pi/2: esin=%g ecos=%g etan2=%g\n", maxs, maxc, maxt);

    double maxs2=0, maxc2=0, maxt2=0;
    for (int i=-10000001; i<=10000001; ++i) {
        double a = 1.03e-9*double(i)*pi2 + 0.237e-9;
        double s = std::sin(a), c = std::cos(a);
        maxs2 = std::max(maxs2, err(a, std::asin(s)));
        maxc2 = std::max(maxc2, err(a, std::acos(c)));
        maxt2 = std::max(maxt2, err(a, std::atan2(s,c)));
    }
    printf("near zero: esin=%g ecos=%g etan2=%g\n", maxs2, maxc2, maxt2);

}

void quatTest() {
    const Real pi2 = std::acos(Real(0));
    Vec4 avOrig(-.1-1e-4, 7e1,-.2,.1);
    Quaternion q1,q2;
    q1.setQuaternionFromAngleAxis(avOrig);
    cout << "q1=" << q1 << endl;
    Vec4 av1 = q1.convertQuaternionToAngleAxis();
    cout << std::setprecision(18);
    cout << "avOrig=" << avOrig << " av1=" << av1 << endl;
    q2.setQuaternionFromAngleAxis(av1);
    cout << "q2-q1=" << q2-q1 << endl;
    Vec4 av2 = q2.convertQuaternionToAngleAxis();
    cout << "av2-av1=" << av2-av1 << endl;

    Rotation r1(q1);
    Vec4 av3 = r1.convertRotationToAngleAxis();
    cout << "av3=" << av3 << " av3-av1=" << av3-av1 << endl;
    Quaternion q3;q3.setQuaternionFromAngleAxis(av3);
    Rotation r2(q3);
    cout << "norm(r2*~r1)-sqrt(3)=" << (r2*~r1).norm()-std::sqrt(3.) << endl;

}

void orthoTest(String msg, const Rotation& R) {
    cout << msg << endl;
    cout << "cols=" << R(0).norm()-1 << ", " << R(1).norm()-1 << ", " << R(2).norm()-1 << endl;
    cout << "rows=" << R[0].norm()-1 << ", " << R[1].norm()-1 << ", " << R[2].norm()-1 << endl;
    cout << "perp=" << dot(R(0),R(1)) << ", " << dot(R(1),R(2)) << ", " << dot(R(0),R(2)) << endl;
}

void f(const CoordinateAxis& ax) { printf("which = %d\n", (int)ax); }
void f(CoordinateAxis::XCoordinateAxis) { printf("X\n"); }
void f(CoordinateAxis::YCoordinateAxis) { printf("Y\n"); }
void f(CoordinateAxis::ZCoordinateAxis) { printf("Z\n"); }

int a[]={9,10,11};

int main() {
    quatTest();

    cout << "-----------------\n";
    cout << "X,Y,Z=" << XAxis << "," <<YAxis << "," <<ZAxis << endl;
    printf("%d %d %d\n", (int)XAxis, (int)YAxis, (int)ZAxis);

    f(XAxis); f(YAxis); f(ZAxis);
    printf("index: %d %d %d\n", a[XAxis], a[YAxis], a[ZAxis]);

try {
    UnitVec3 u;
    cout << "should be NaN: u()=" << u << endl;
    cout << "typename(u)=" << typeid(u).name() << " ~u=" << typeid(~u).name() << endl;

    u = UnitVec3(1,2,3);
    cout << "u=UnitVec(1,2,3): " << u << " transposed: " << ~u << endl;
    cout << "u[2]=" << u[2] << " u(2)=" << u(2) << endl;
    cout << "(~u)[2]=" << (~u)[2] << " (~u)(2)=" << (~u)(2) << endl;

    cout << "u.perp()=" << u.perp() << " ~u*u.perp()=" << ~u*u.perp() << endl;
    cout << "(~u).perp()=" << (~u).perp() << " (~u).perp()*u=" << (~u).perp()*u << endl;

    Vec3 vv(u);
    Vec3 w = u;
    cout << "Vec3 vv(u)=" << vv << "  Vec3 w=u, w=" << w << endl;

    // these lines shouldn't compile
    //u = Vec3(1,2,3);
    //u += UnitVec3(1,2,3);
    //u[2] = 5.;
    //u(2) = 5.;
    //u=0.;
    // end of lines which shouldn't compile


    Rotation R;
    cout << "should be identity: R()=" << R;
    cout << "typename(R)=" << typeid(R).name() << " ~R=" << typeid(~R).name() << endl;

    Transform X;
    cout << "should be identity: X()=" << X;
    cout << "typename(X)=" << typeid(X).name() << " ~X=" << typeid(~X).name() << endl;

	R = Rotation( u, ZAxis );

    cout << "Rotation with u as z axis (norm=" << R.norm() << "): " << R; 
    cout << "~R: " << ~R;
    cout << "typename(R[1])=" << typeid(R[1]).name() << " R(1)=" << typeid(R(1)).name() << endl;

    X = Transform(R, Vec3(-1, 2, 20));
    cout << "Transform X=" << X;
    cout << "   X.R=" << X.R() << "  X.p=" << X.p() << "  X.pInv=" << X.pInv() << endl;
    cout << "   X.x=" << X.x() << "  y=" << X.y() << "  z=" << X.z() << endl;
    cout << "Transform ~X=" << ~X;
    cout << "   (~X).R=" << (~X).R() << "  (~X).p=" << (~X).p() << "  (~X).pInv=" << (~X).pInv() << endl;
    cout << "   (~X).x=" << (~X).x() << "  y=" << (~X).y() << "  z=" << (~X).z() << endl;    
    cout << "Transform X asMat34():" << X.asMat34();
    cout << "Transform X toMat34():" << X.toMat34();
    cout << "Transform X toMat44():" << X.toMat44();
    cout << "Transform ~X toMat34():" << (~X).toMat34();
    cout << "Transform ~X toMat44():" << (~X).toMat44();
    cout << "should be identity: X*~X=" << X*~X;
    cout << "should be identity: ~X*X=" << ~X*X;

    Vec3 v(1,2,3);
    cout << "v=" << v << endl;
    cout << "X.R()*v=" << X.R()*v << endl;
    cout << "X*v=" << X*v << endl;
    cout << "X*[v,0]=" << X*v.append1(0) << endl;
    cout << "X*[v,1]=" << X*v.append1(1) << endl;
    cout << "~X*(X*v)=" << ~X*(X*v) << endl;
    cout << "~X*(X*[v,0])=" << ~X*(X*v.append1(0)) << endl;
    cout << "~X*(X*[v,1])=" << ~X*(X*v.append1(1)) << endl;

    Rotation rx( 0.1,  XAxis );
    Rotation ry( 0.17, YAxis );
    Rotation rz( 0.31, ZAxis );
    Rotation rxoldy( SpaceRotationSequence, 0.1, XAxis, 0.17, YAxis );
    Rotation rxyz = rz*ry*rx; // space fixed
    // reverse order to use body fixed like space fixed
    Rotation r123( BodyRotationSequence, 0.31, ZAxis, 0.17, YAxis, 0.1, XAxis );

    cout << "Rotation r123:" << r123;
    cout << "Rotation r123.asMat33():" << r123.asMat33();
    cout << "Rotation r123.toMat33():" << r123.toMat33();

    cout << "InvRotation ~r123:" << (~r123);
    cout << "InvRotation ~r123.asMat33():" << (~r123).asMat33();
    cout << "InvRotation ~r123.toMat33():" << (~r123).toMat33();

    cout << "ry(.17)*rx(.1)=" << ry*rx;
    cout << "rxoldy(.1,.17)=" << rxoldy;
    cout << "norm(rxoldy-ry*rx)=" << (rxoldy-ry*rx).norm() << endl;
    cout << "rxyz=" << rxyz;
    cout << "r123=" << r123;

    // Test inverse assignments, copy construct
    Rotation invr123(~r123);
    Rotation invr123eq;
    invr123eq = ~r123;

    cout << "Check inverse assignment/copy." << endl;
    cout << "norm invr123*r123-identity=" << (invr123*r123-Mat33(1)).norm() << endl;
    cout << "norm invr123eq*r123-identity=" << (invr123eq*r123-Mat33(1)).norm() << endl;

    Rotation bodyXY( BodyRotationSequence, 0.03, XAxis, 0.11, YAxis, 0.0, ZAxis );
	Rotation aboutYthenX( SpaceRotationSequence, 0.11, YAxis, 0.03, XAxis );

    cout << "bodyXY(.03,.11)=" << bodyXY;
    cout << "aboutYthenoldX(.11,.03)=" << aboutYthenX;

    Rotation bodyZX( BodyRotationSequence, 0.07, ZAxis, 0.0, YAxis, -0.22, XAxis );
    Rotation aboutZThenNewX( BodyRotationSequence, 0.07, ZAxis, -0.22, XAxis );
    Rotation aboutZ( 0.07, ZAxis );
    Rotation aboutXaboutZ = Rotation(-0.22, aboutZ*Vec3(1,0,0)) * aboutZ;
    cout << "bodyZX(.07,-.22)=" << bodyZX;
    cout << "aboutZThenNewX(.07,-.22)=" << aboutZThenNewX;
    cout << "aboutXaboutZ(.07,-.22)=" << aboutXaboutZ;
    cout << "  diff norm=" << (aboutXaboutZ-aboutZThenNewX).norm() << endl;

    // Check all space-fixed sequences vs composed single rotations.
    cout << "SXY err=" << ( Rotation( SpaceRotationSequence, 0.13, XAxis, -0.29, YAxis ) - Rotation( -0.29, YAxis )*Rotation( 0.13, XAxis) ).norm() << endl;
    cout << "SXZ err=" << ( Rotation( SpaceRotationSequence, 0.13, XAxis, -0.29, ZAxis ) - Rotation( -0.29, ZAxis )*Rotation( 0.13, XAxis) ).norm() << endl;
    cout << "SZX err=" << ( Rotation( SpaceRotationSequence, 0.13, ZAxis, -0.29, XAxis ) - Rotation( -0.29, XAxis )*Rotation( 0.13, ZAxis) ).norm() << endl;
    cout << "SYZ err=" << ( Rotation( SpaceRotationSequence, 0.13, YAxis, -0.29, ZAxis ) - Rotation( -0.29, ZAxis )*Rotation( 0.13, YAxis) ).norm() << endl;
    cout << "SZY err=" << ( Rotation( SpaceRotationSequence, 0.13, ZAxis, -0.29, YAxis ) - Rotation( -0.29, YAxis )*Rotation( 0.13, ZAxis) ).norm() << endl;

    // Check all body-fixed sequences vs. transformed & composed single rotations.
    const Real a1=-.23, a2=1.09;
    const Rotation ax( a1, XAxis ); const Vec3 x(1,0,0);
    const Rotation ay( a1, YAxis ); const Vec3 y(0,1,0);
    const Rotation az( a1, ZAxis ); const Vec3 z(0,0,1);

    cout << "BXY err=" << ( Rotation( BodyRotationSequence, a1, XAxis, a2, YAxis ) - Rotation(a2,ax*y)*ax ).norm() << endl;
    cout << "BYX err=" << ( Rotation( BodyRotationSequence, a1, YAxis, a2, XAxis ) - Rotation(a2,ay*x)*ay ).norm() << endl;
    cout << "BXZ err=" << ( Rotation( BodyRotationSequence, a1, XAxis, a2, ZAxis ) - Rotation(a2,ax*z)*ax ).norm() << endl;
    cout << "BZX err=" << ( Rotation( BodyRotationSequence, a1, ZAxis, a2, XAxis ) - Rotation(a2,az*x)*az ).norm() << endl;
    cout << "BYZ err=" << ( Rotation( BodyRotationSequence, a1, YAxis, a2, ZAxis ) - Rotation(a2,ay*z)*ay ).norm() << endl;
    cout << "BZY err=" << ( Rotation( BodyRotationSequence, a1, ZAxis, a2, YAxis ) - Rotation(a2,az*y)*az ).norm() << endl;

    cout << "rotate1(0,.03)="  << ( rotate1(0,.03) - Rotation(.03,XAxis) ).norm() << endl;
    cout << "rotate1(1,.03)="  << ( rotate1(1,.03) - Rotation(.03,YAxis) ).norm() << endl;
    cout << "rotate1(2,.03)="  << ( rotate1(2,.03) - Rotation(.03,ZAxis) ).norm() << endl;

    cout << "-1 % 3=" << -1 % 3 << endl;

    const Transform X_01( Rotation( .03, Vec3(1,.1,.7) ), Vec3(1,2,3) );
    const Real masses[] = {1,2,3,4,5};
    const Real mtot = Vector(5,masses).sum();
    const Vec3 stations000[] = {Vec3(.1,.2,.3), Vec3(-2,-9,1), Vec3(.01,.02,.05),
                                Vec3(1,-1,1), Vec3(0,0,0)};
    const Vector_<Vec3> s000(5, stations000);
    const Vec3 com000 = (~Vector(5,masses) * s000)/mtot;
    Vector_<Vec3> s123(s000); 
    for (int i=0;i<5;++i) s123[i] = ~X_01*s123[i];
    const Vec3 com123 = (~Vector(5,masses) * s123)/mtot;
    Inertia I000(0), I123(0);
    for (int i=0; i<5; ++i) I000 += Inertia(s000[i],masses[i]);
    for (int i=0; i<5; ++i) I123 += Inertia(s123[i],masses[i]);
    cout << "mtot=" << mtot << endl;
    cout << "com000=" << com000 << endl;
    cout << "com123=" << com123 << endl;
    cout << "I000=" << I000;
    cout << "I123=" << I123;
    
    MassProperties mp(mtot, com000, I000);
    cout << mp;
    Real scale = std::sqrt(I123.toMat33().diag().normSqr()/3); cout << "inertia scale=" << scale << endl;
    cout << "norm(I000->123-I123)/rms(I123)=" 
        << (mp.calcTransformedInertia(X_01).toMat33()-I123.toMat33()).norm()/scale
        << endl;


    Rotation R_AB( BodyRotationSequence,   0.31, ZAxis, 0.17, YAxis, 0.1, XAxis ); 
    Rotation R_BC( BodyRotationSequence, -123.3, ZAxis, 41.1, YAxis, 14,  XAxis );
		
    Rotation Rtmp;
    cout << "R_AB*R_BC=" << R_AB*R_BC;
    Rtmp = R_AB;
    cout << "R_AB*=R_BC=" << (Rtmp*=R_BC);
    cout << "R_AB*~R_BC=" << R_AB*~R_BC;
    cout << "R_AB/R_BC=" << R_AB/R_BC;
    Rtmp=R_AB;
    cout << "R_AB/=R_BC=" << (Rtmp/=R_BC);

    cout << "COMPARE ROTATIONS" << endl;
    Rotation R_GB, R_GX; Real backOut;
    R_GB.setRotationFromAngleAboutNonUnitVector( 0.17, Vec3(1,2,3) );
    backOut = R_GB.convertRotationToAngleAxis()[0];
    cout << " in=0.17 radians, out=" << backOut << " err=" << std::abs(backOut-0.17) << endl;

    R_GB.setRotationFromAngleAboutNonUnitVector( 0.17+1e-13, Vec3(1,2,3) );
    R_GX.setRotationFromAngleAboutNonUnitVector( 0.17,       Vec3(1,2,3) );
    cout << " 0.17+1e-13:0.17 isSameToPrecision? " << R_GB.isSameRotationToWithinAngleOfMachinePrecision(R_GX)
         << " isSameToAngle(1e-12)? " << R_GB.isSameRotationToWithinAngle(R_GX, 1e-12) << endl;
    R_GB.setRotationFromAngleAboutNonUnitVector( 0.17+1e-15, Vec3(1,2,3) );
    R_GX.setRotationFromAngleAboutNonUnitVector( 0.17,       Vec3(1,2,3) );
    cout << " 0.17+1e-15:0.17 isSameToPrecision? " << R_GB.isSameRotationToWithinAngleOfMachinePrecision(R_GX)
         << " isSameToAngle(1e-18)? " << R_GB.isSameRotationToWithinAngle(R_GX, 1e-18) << endl;

    const Real pi2 = Pi/2;
    const Real pi2x = -pi2 + 1e-8;
    cout << "pi2x=pi2-" << pi2-pi2x << " sin(pi2x)-1=" << std::sin(pi2x)-1 << endl;
    const Vec3 vin(-3, pi2x, 0.1);
    Rotation b123( BodyRotationSequence, vin[0], XAxis, vin[1], YAxis, vin[2], ZAxis );  
    Mat33 m123=b123; m123[0][0] += 1e-14; m123[1][2] += 1e-14;
    //b123 = Rotation::trustMe(m123);
    //cout << "bad  b123*~b123 angle=" << (b123*~b123).convertToAngleAxis()[0] << endl;
    //b123 = Rotation(m123);
    //cout << "good b123*~b123.norm()=" << (b123*~b123).convertToAngleAxis()[0] << endl;
    b123 *= R_GX;
    b123 *= ~R_GX;
    cout << "b123=" << b123;
    Rotation cleanb123(b123.asMat33());
   // b123=cleanb123;

    cout << "cos(pi2x)=" << cos(pi2x) << "cos(pi2)=" << cos(pi2) << endl;
    cout << "atan2(02/00)-pi2=" << atan2(b123[0][2],b123[0][0])-pi2 << endl;
    Vec3 v123 = b123.convertThreeAxesRotationToThreeAngles( BodyRotationSequence, XAxis, YAxis, ZAxis );
    Vec4 aax = b123.convertRotationToAngleAxis();
    cout << "vin=" << vin << "\nvout=" << v123 << endl;
    Rotation b123x( BodyRotationSequence, v123[0], XAxis, v123[1], YAxis, v123[2], ZAxis );
    Vec4 aax2 = (~b123*b123x).convertRotationToAngleAxis();
    cout << " aax2=" << aax2 << endl;

    Rotation chrisRot( UnitVec3( Vec3(1,1,0) ), XAxis, Vec3(0,1,0), YAxis );
    cout << "chrisRot+y=" << chrisRot;
    chrisRot = Rotation( UnitVec3( Vec3(1,1,0) ), XAxis, Vec3(0,-1,0), YAxis );
    cout << "chrisRot-y=" << chrisRot;
    chrisRot = Rotation( UnitVec3( Vec3(1,1,0) ), XAxis, Vec3(0,0,1), YAxis );
    cout << "chrisRot+z=" << chrisRot;
    return 0;
}
catch(const Exception::Base& e) {
    std::cout << e.getMessage() << std::endl;
}
}

static Rotation rotate1(int i, Real a) {
    assert(0 <= i && i < 3);
    int j=(i+1)%3, k=(i+2)%3;
    Real s=std::sin(a), c=std::cos(a);
    Mat33 m;
    m(i,i)=1; m(i,j)=m(j,i)=m(i,k)=m(k,i)=0;
    m(j,j)=m(k,k)=c;
    m(k,j)=s; m(j,k)=-s;
	return Rotation(m,true);
}


// Two-angle space-fixed rotations.

//0,1
static Rotation aboutXThenOldY(const Real& xInRad, const Real& yInRad) {
    const Real s0 = std::sin(xInRad), c0 = std::cos(xInRad);
    const Real s1 = std::sin(yInRad), c1 = std::cos(yInRad);
    const Mat33 m( c1   ,  s0*s1 ,  c0*s1 ,
                    0   ,   c0   ,  -s0   ,
                  -s1   ,  s0*c1 ,  c0*c1 ); return Rotation(m,true);
}
//2,0
static Rotation aboutZThenOldX(const Real& zInRad, const Real& xInRad) {
    const Real s1 = std::sin(xInRad), c1 = std::cos(xInRad);
    const Real s0 = std::sin(zInRad), c0 = std::cos(zInRad);
    const Mat33 m(  c0   ,  -s0   ,    0   ,
                   s0*c1 ,  c0*c1 ,  -s1   ,
                   s0*s1 ,  c0*s1 ,   c1   ); return Rotation(m,true);
}
//1,2
static Rotation aboutYThenOldZ(const Real& yInRad, const Real& zInRad) {
    const Real s0 = std::sin(yInRad), c0 = std::cos(yInRad);
    const Real s1 = std::sin(zInRad), c1 = std::cos(zInRad);
    const Mat33 m( c0*c1 ,  -s1   ,  s0*c1 ,
                   c0*s1 ,   c1   ,  s0*s1 ,
                   -s0   ,    0   ,   c0   ); return Rotation(m,true);
}

//1,0
static Rotation aboutYThenOldX(const Real& yInRad, const Real& xInRad) {
    const Real s1 = std::sin(xInRad), c1 = std::cos(xInRad);
    const Real s0 = std::sin(yInRad), c0 = std::cos(yInRad);
    const Mat33 m(  c0   ,    0   ,   s0   ,
                   s0*s1 ,   c1   , -c0*s1 ,
                  -s0*c1 ,   s1   ,  c0*c1 ); return Rotation(m,true);
}
//0,2
static Rotation aboutXThenOldZ(const Real& xInRad, const Real& zInRad) {
    const Real s0 = std::sin(xInRad), c0 = std::cos(xInRad);
    const Real s1 = std::sin(zInRad), c1 = std::cos(zInRad);
    const Mat33 m(  c1   , -c0*s1 ,  s0*s1 ,
                    s1   ,  c0*c1 , -s0*c1 ,
                     0   ,   s0   ,   c0   ); return Rotation(m,true);
}

//2,1
static Rotation aboutZThenOldY(const Real& zInRad, const Real& yInRad) {
    const Real s1 = std::sin(yInRad), c1 = std::cos(yInRad);
    const Real s0 = std::sin(zInRad), c0 = std::cos(zInRad);
    const Mat33 m( c0*c1 , -s0*c1 ,   s1   ,
                    s0   ,   c0   ,    0   ,
                  -c0*s1 ,  s0*s1 ,   c1   ); return Rotation(m,true);
}
