#ifndef SimTK_SIMMATH_GEO_BICUBIC_BEZIER_PATCH_H_
#define SimTK_SIMMATH_GEO_BICUBIC_BEZIER_PATCH_H_

/* -------------------------------------------------------------------------- *
 *                        SimTK Simbody: SimTKmath                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2011 Stanford University and the Authors.           *
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

/** @file
Provides primitive operations for a single bicubic Bezier patch using either
single or double precision. **/

#include "SimTKcommon.h"
#include "simmath/internal/common.h"
#include "simmath/internal/Geo.h"

#include <cassert>
#include <cmath>
#include <algorithm>

namespace SimTK {


//==============================================================================
//                         GEO BICUBIC BEZIER PATCH
//==============================================================================
/** A primitive useful for computations involving a single bicubic Bezier
patch. Note that a bicubic Bezier spline surface would not necessarily be
composed of these, but could use the static methods here for patch 
computations. 

<h3>Theory</h3>
The primary reference for this implementation is the book "Geometric Modeling, 
3rd ed." by Michael E. Mortenson, Industrial Press 2006, chapter 8. We follow
Mortenson's notation here (with some name changes) and equation numbers are 
from the text. See CubicHermiteCurve_ and BicubicHermiteSurface_ comments for 
introductory material; here we add the Bezier description to the algebraic and 
Hermite (geometric) forms described there.

We use B for Bezier control points (rather than P), and H for Hermite 
coefficients (rather than B). We call the Bezier basis matrix Mb (same as
Mortenson) but call the Hermite basis matrix Mh (rather than Mf). 

The 16 control points are laid out like this, matching Mortenson, with the
Hermite coefficient matrix for comparison:
<pre>
        [ b11 b12 b13 b14 ]         [ h00  h01  w00  w01 ]    u=dp/du
    B = [ b21 b22 b23 b24 ]     H = [ h10  h11  w10  w11 ]    w=dp/dw
        [ b31 b32 b33 b34 ]         [ u00  u01  t00  t01 ]    t=d2p/dudw
        [ b41 b42 b43 b44 ]         [ u10  u11  t10  t11 ]      ("twist")
</pre>
We store those in a 4x4 hypermatrix; subtract one from each control point
index to get the matrix indices (sorry). 

To convert between Bezier and Hermite or algebraic forms, use
<pre>
    A = Mb B ~Mb  (although Mb is symmetric so ~Mb=Mb)
    H = (Mh^-1 Mb) B ~(Mh^-1 Mb)
    B = (Mb^-1 Mh) H ~(Mb^-1 Mh)
</pre>
where the parenthesized (very simple) matrices are given explicitly in the 
Theory section for CubicBezierCurve_. The results are:
<pre>
    [    b11         b14           3(b12-b11)             3(b14-b13)     ]
H = [    b41         b44           3(b42-b41)             3(b44-b43)     ]
    [ 3(b21-b11)  3(b24-b14)  9(b22-b21-{b12-b11})  9({b24-b14}+b13-b23) ]
    [ 3(b41-b31)  3(b44-b34)  9({b42-b41}+b31-b32)  9({b44-b34}+b33-b43) ]
</pre>
The braces in the twist terms mark common subexpressions, so the
conversion from Bezier to Hermite takes 3x28=84 flops (all entries are 
3-vectors).
<pre>
    [   h00            h00+w00/3               h01-w01/3          h01    ]
B = [h00+u00/3 {h00+u00/3}+w00/3+t00/9 {h01+u01/3}-w01/3-t01/9 h01+u01/3 ]
    [h10-u10/3 {h10-u10/3}+w10/3-t10/9 {h11-u11/3}-w11/3+t11/9 h11-u11/3 ]
    [   h10            h10+w10/3               h11-w11/3          h11    ] 
Exploiting the common subexpressions in braces, conversion from Hermite to 
Bezier takes 3x32=96 flops.

**/
template <class P>
class Geo::BicubicBezierPatch_ {
typedef P               RealP;
typedef Vec<3,RealP>    Vec3P;

public:
/** Construct an uninitialized patch; control points will be garbage. **/
BicubicBezierPatch_() {}

/** Construct a bicubic Bezier patch using the given control points B. See the
class documentation for how this matrix is defined. If instead you
have algebraic or Hermite coefficients, you can first convert them to Bezier
control points with the calcBFromA() or calcHFromA() provided here. **/
explicit BicubicBezierPatch_(const Mat<4,4,Vec3P>& controlPoints) 
:   B(controlPoints) {} 


/** Evaluate a point P(u,w) on this patch given values for the
parameters u and w in [0,1]. Values outside this range are permitted but do 
not lie on the patch. Cost is 123 flops. **/
Vec3P evalP(RealP u, RealP w) const {return evalPUsingB(B,u,w);}

/** Evaluate the tangents Pu=dP/du, Pw=dP/dw on this patch given values for the
parameters u and w in [0,1]. Values outside this range are permitted but do not
lie on the curve segment. Cost is 248 flops. **/
void evalP1(RealP u, RealP w, Vec3P& Pu, Vec3P& Pw) const 
{   return evalP1UsingB(B,u,w,Pu,Pw); }

/** Evaluate the second derivatives Puu=d2P/du2, Pww=d2P/dw2, and cross
derivative Puw=Pwu=d2P/dudw on this patch given values for the parameters u 
and w in [0,1]. Values outside this range are permitted but do not lie on the 
curve segment. Cost is 363 flops. **/
void evalP2(RealP u, RealP w, Vec3P& Puu, Vec3P& Puw, Vec3P& Pww) const 
{   evalP2UsingB(B,u,w,Puu,Puw,Pww); }

/** Evaluate the third derivatives Puuu=d3P/du3, Pwww=d3P/dw3, and cross
derivatives Puuw=Pwuu=Puwu=d3P/du2dw and Puww=Pwwu=Pwuw=d3P/dudw2 on this patch 
given values for the parameters u and w in [0,1]. Cost is 468 flops. All
higher derivatives of a cubic patch are zero. **/
void evalP3(RealP u, RealP w, Vec3P& Puuu, Vec3P& Puuw, 
                              Vec3P& Puww, Vec3P& Pwww) const 
{   evalP3UsingB(B,u,w,Puuu,Puuw,Puww,Pwww); }

/** Return a reference to the Bezier control points B that are
stored in this object. See the documentation for this class to see how the
returned matrix of control points is defined. **/
const Mat<4,4,Vec3P>& getControlPoints() const {return B;}
/** Calculate the algebraic coefficients A from the stored Bezier control
points. See the documentation for BicubicHermitePatch_ to see how the
returned matrix of coefficients is defined. Cost is 240 flops. **/
Mat<4,4,Vec3P> calcAlgebraicCoefficients() const {return calcAFromB(B);}
/** Calculate the Hermite coefficients H from the stored 
Bezier control points. See the documentation for this class to see how the
returned matrix of coefficients is defined. Cost is 84 flops. **/
Mat<4,4,Vec3P> calcHermiteCoefficients() const {return calcHFromB(B);}

/** Return the u=0 boundary curve as a Bezier curve segment. The control points
are just the first row of the patch control points B so no computation is 
required. **/
CubicBezierCurve_<P> getBoundaryCurveU0() const 
{   return CubicBezierCurve_<P>(B[0]); } // b11 b12 b13 b14
/** Return the u=1 boundary curve as a Bezier curve segment. The control points
are just the last row of the patch control points B so no computation is 
required. **/
CubicBezierCurve_<P> getBoundaryCurveU1() const 
{   return CubicBezierCurve_<P>(B[3]); } // b41 b42 b43 b44
/** Return the w=0 boundary curve as a Bezier curve segment. The control points
are just the first column of the patch control points B so no computation is 
required. **/
CubicBezierCurve_<P> getBoundaryCurveW0() const 
{   return CubicBezierCurve_<P>(B(0)); } // b11 b21 b31 b41
/** Return the w=1 boundary curve as a Bezier curve segment. The control points
are just the last column of the patch control points B so no computation is 
required. **/
CubicBezierCurve_<P> getBoundaryCurveW1() const 
{   return CubicBezierCurve_<P>(B(3)); } // b14 b24 b34 b44

/** Given a particular value u0 for patch coordinate u, create a cubic Bezier
curve segment P(w)=P(u0,w) for the isoparametric curve along the patch at fixed
u=u0. Cost is 93 flops. **/
CubicBezierCurve_<P> calcIsoCurveU(RealP u0) const 
{   return calcIsoCurveU(B, u0); }
/** Given a particular value w0 for patch coordinate w, create a cubic Bezier
curve segment P(u)=P(u,w0) for the isoparametric curve along the patch at fixed
w=w0. Cost is 93 flops. **/
CubicBezierCurve_<P> calcIsoCurveW(RealP w0) const 
{   return calcIsoCurveW(B, w0); }



/**@name                 Utility methods
These static methods work with given control points. **/
/**@{**/

/** Given Bezier control points B and values for the curve parameters 
u and w in [0..1], return the point P(u,w)=Fb(u)*B*~Fb(w) at that location, 
where Fb is a vector of Bezier basis functions. Cost is 
3x35+18=123 flops. **/
static Vec3P evalPUsingB(const Mat<4,4,Vec3P>& B, RealP u, RealP w) { 
    Row<4,P> Fbu = CubicBezierCurve_<P>::calcFb(u);     // 9 flops
    Row<4,P> Fbw = CubicBezierCurve_<P>::calcFb(w);     // 9 flops
    return Fbu * B * ~Fbw;                              // 3x35 flops
}

/** Given Bezier control points B and values for the curve parameters 
u and w in [0..1], return the tangents Pu(u,w)=dFb(u)*B*~Fb(w) and 
Pw(u,w)=Fb(u)*B*~dFb(w) at 
that location. Cost is 3x70+38=248 flops. **/
static void evalP1UsingB(const Mat<4,4,Vec3P>& B, RealP u, RealP w,
                         Vec3P& Pu, Vec3P& Pw) {
    Row<4,P> Fbu  = CubicBezierCurve_<P>::calcFb(u);     //  9 flops
    Row<4,P> Fbw  = CubicBezierCurve_<P>::calcFb(w);     //  9 flops
    Row<4,P> dFbu = CubicBezierCurve_<P>::calcDFb(u);    // 10 flops
    Row<4,P> dFbw = CubicBezierCurve_<P>::calcDFb(w);    // 10 flops
    Pu = dFbu * B * ~Fbw;                                // 3x35
    Pw = Fbu  * B * ~dFbw;                               // 3x35
}

/** Given Bezier control points B and values for the curve parameters 
u and w in [0..1], return the second derivatives Puu(u,w)=d2Fb(u)*B*~Fb(w),
Puw(u,w)=dFb(u)*B*~dFb(w) and Pww(u,w)=Fb(u)*B*~d2Fb(w) at that location. 
Cost is 3x105+48=363 flops. **/
static void evalP2UsingB(const Mat<4,4,Vec3P>& B, RealP u, RealP w,
                         Vec3P& Puu, Vec3P& Puw, Vec3P& Pww) {
    Row<4,P> Fbu   = CubicBezierCurve_<P>::calcFb(u);     //  9 flops
    Row<4,P> Fbw   = CubicBezierCurve_<P>::calcFb(w);     //  9 flops
    Row<4,P> dFbu  = CubicBezierCurve_<P>::calcDFb(u);    // 10 flops
    Row<4,P> dFbw  = CubicBezierCurve_<P>::calcDFb(w);    // 10 flops
    Row<4,P> d2Fbu = CubicBezierCurve_<P>::calcD2Fb(u);   //  5 flops
    Row<4,P> d2Fbw = CubicBezierCurve_<P>::calcD2Fb(w);   //  5 flops
    Puu = d2Fbu * B * ~Fbw;                               // 3x35
    Puw = dFbu  * B * ~dFbw;                              // 3x35
    Pww = Fbu   * B * ~d2Fbw;                             // 3x35
}

/** Given Bezier control points B and values for the curve parameters 
u and w in [0..1], return the third derivatives Puuu(u,w)=d3Fb(u)*B*~Fb(w),
Puuw(u,w)=d2Fb(u)*B*~dFb(w), Puww(u,w)=dFb(u)*B*~d2Fb(w) and 
Pwww(u,w)=Fb(u)*B*~d3Fb(w) at that location. Cost is 3x140+48=468 flops. **/
static void evalP3UsingB(const Mat<4,4,Vec3P>& B, RealP u, RealP w,
                         Vec3P& Puuu, Vec3P& Puuw, Vec3P& Puww, Vec3P& Pwww) {
    Row<4,P> Fbu   = CubicBezierCurve_<P>::calcFb(u);     //  9 flops
    Row<4,P> Fbw   = CubicBezierCurve_<P>::calcFb(w);     //  9 flops
    Row<4,P> dFbu  = CubicBezierCurve_<P>::calcDFb(u);    // 10 flops
    Row<4,P> dFbw  = CubicBezierCurve_<P>::calcDFb(w);    // 10 flops
    Row<4,P> d2Fbu = CubicBezierCurve_<P>::calcD2Fb(u);   //  5 flops
    Row<4,P> d2Fbw = CubicBezierCurve_<P>::calcD2Fb(w);   //  5 flops
    Row<4,P> d3Fbu = CubicBezierCurve_<P>::calcD3Fb(u);   //  0 
    Row<4,P> d3Fbw = CubicBezierCurve_<P>::calcD3Fb(w);   //  0
    Puuu = d3Fbu * B * ~Fbw;                              // 3x35
    Puuw = d2Fbu * B * ~dFbw;                             // 3x35
    Puww = dFbu  * B * ~d2Fbw;                            // 3x35
    Pwww = Fbu   * B * ~d3Fbw;                            // 3x35
}

/** Given a particular value u0 for patch coordinate u, create a cubic Bezier
curve segment P(w)=P(u0,w) for the isoparametric curve along the patch at fixed
u=u0. Cost is 3x28+9=93 flops. **/
static CubicBezierCurve_<P> 
calcIsoCurveU(const Mat<4,4,Vec3P>& B, RealP u0) 
{   const Row<4,Vec3P> Bu0 = CubicBezierCurve_<P>::calcFb(u0) * B;
    return CubicBezierCurve_<P>(Bu0); }

/** Given a particular value w0 for patch coordinate w, create a cubic Bezier
curve segment P(u)=P(u,w0) for the isoparametric curve along the patch at fixed
w=w0. Cost is 3x28+9=93 flops. **/
static CubicBezierCurve_<P> 
calcIsoCurveW(const Mat<4,4,Vec3P>& B, RealP w0)
{    const Vec<4,Vec3P> Bw0 = B * ~CubicBezierCurve_<P>::calcFb(w0);
     return CubicBezierCurve_<P>(Bw0); }

/** Given the Bezier control points B, return the equivalent vector algebraic
coefficients A. All coefficients are 3-vectors. Cost is 3x80=240 flops. **/
static Mat<4,4,Vec3P> calcAFromB(const Mat<4,4,Vec3P>& B) {
    typedef const Vec3P& Coef;
    Coef b11=B(0,0), b12=B(0,1), b13=B(0,2), b14=B(0,3), 
         b21=B(1,0), b22=B(1,1), b23=B(1,2), b24=B(1,3), 
         b31=B(2,0), b32=B(2,1), b33=B(2,2), b34=B(2,3), 
         b41=B(3,0), b42=B(3,1), b43=B(3,2), b44=B(3,3); 
    // First calculate Mb*B:
    //       a   b   c   d
    //       e   f   g   h
    //       p   q   r   s
    //      b11 b12 b13 b14
    Vec3P a= b41-b11+3*(b21-b31), b= b42-b12+3*(b22-b32),   // 3x16 flops
          c= b43-b13+3*(b23-b33), d= b44-b14+3*(b24-b34);
    Vec3P e= 3*(b11+b31)-6*b21,   f= 3*(b12+b32)-6*b22,     // 3x16 flops
          g= 3*(b13+b33)-6*b23,   h= 3*(b14+b34)-6*b24;
    Vec3P p= 3*(b21-b11),         q= 3*(b22-b12),           // 3x8 flops
          r= 3*(b23-b13),         s= 3*(b24-b14);

    // Then calculate (Mb*B)*~Mb. (3x40 more flops)
    return Mat<4,4,Vec3P>
       ( d-a+3*(b-c),          3*(a+c)-6*b,       3*(b-a),        a,
         h-e+3*(f-g),          3*(e+g)-6*f,       3*(f-e),        e,
         s-p+3*(q-r),          3*(p+r)-6*q,       3*(q-p),        p,
         b14-b11+3*(b12-b13),  3*(b11+b13)-6*b12, 3*(b12-b11),   b11 );
}

/** Given the vector algebraic coefficients A, return the equivalent Bezier
control points B. All coefficients are 3-vectors. Cost is 3x72=216 flops. **/
static Mat<4,4,Vec3P> calcBFromA(const Mat<4,4,Vec3P>& A) {
    typedef const Vec3P& Coef;
    Coef a33=A(0,0), a32=A(0,1), a31=A(0,2), a30=A(0,3), 
         a23=A(1,0), a22=A(1,1), a21=A(1,2), a20=A(1,3), 
         a13=A(2,0), a12=A(2,1), a11=A(2,2), a10=A(2,3), 
         a03=A(3,0), a02=A(3,1), a01=A(3,2), a00=A(3,3); 
    // First calculate Mb^-1*A:
    //      a03 a02 a01 a00
    //       a   b   c   d
    //       e   f   g   h
    //       p   q   r   s
    Vec3P a=a13/3+a03, b=a12/3+a02, c=a11/3+a01, d=a10/3+a00;   // 3x8 flops
    Vec3P e=(a23+2*a13)/3+a03, f=(a22+2*a12)/3+a02,             // 3x16 flops
          g=(a21+2*a11)/3+a01, h=(a20+2*a10)/3+a00;
    Vec3P p=a33+a23+a13+a03, q=a32+a22+a12+a02,                 // 3x12 flops
          r=a31+a21+a11+a01, s=a30+a20+a10+a00;


    // Then calculate (Mb^-1*A)*Mb^-T (3x36 more flops)
    return Mat<4,4,Vec3P>
       ( a00,   a01/3+a00,  (a02+2*a01)/3+a00,     a03+a02+a01+a00,
          d,      c/3+d,       (b+2*c)/3+d,             a+b+c+d,
          h,      g/3+h,       (f+2*g)/3+h,             e+f+g+h,
          s,      r/3+s,       (q+2*r) /3+s,            p+q+r+s     );
}

/** Given the Bezier control points B, return the equivalent vector Hermite
coefficients H. All coefficients are 3-vectors. Cost is 3x28=84 flops. **/
static Mat<4,4,Vec3P> calcHFromB(const Mat<4,4,Vec3P>& B) {
    typedef const Vec3P& Coef;
    Coef b11=B(0,0), b12=B(0,1), b13=B(0,2), b14=B(0,3), 
         b21=B(1,0), b22=B(1,1), b23=B(1,2), b24=B(1,3), 
         b31=B(2,0), b32=B(2,1), b33=B(2,2), b34=B(2,3), 
         b41=B(3,0), b42=B(3,1), b43=B(3,2), b44=B(3,3); 

    // First calculate a few temps -- see class comments. (3x4 flops)
    Vec3P b12mb11=b12-b11, b24mb14=b24-b14, b42mb41=b42-b41, b44mb34=b44-b34;

    // Then calculate (Mh^-1 Mb) B ~(Mh^-1 Mb). (3x24 more flops)
    return Mat<4,4,Vec3P>
       (     b11,         b14,          3*b12mb11,           3*(b14-b13),
             b41,         b44,          3*b42mb41,           3*(b44-b43),
         3*(b21-b11),  3*b24mb14,  9*(b22-b21-b12mb11),  9*(b24mb14+b13-b23),
         3*(b41-b31),  3*b44mb34,  9*(b42mb41+b31-b32),  9*(b44mb34+b33-b43) );
}

/** Given the vector Hermite coefficients H, return the equivalent Bezier
control points B. All coefficients are 3-vectors. Cost is 3x32=96 flops. **/
static Mat<4,4,Vec3P> calcBFromH(const Mat<4,4,Vec3P>& H) {
    typedef const Vec3P& Coef;
    Coef h00=H(0,0), h01=H(0,1), w00=H(0,2), w01=H(0,3), 
         h10=H(1,0), h11=H(1,1), w10=H(1,2), w11=H(1,3), 
         u00=H(2,0), u01=H(2,1), t00=H(2,2), t01=H(2,3), 
         u10=H(3,0), u11=H(3,1), t10=H(3,2), t11=H(3,3); 

    // First calculate a few temps -- see class comments. (3x8 flops)
    Vec3P tmp00=h00+u00/3, tmp01=h01+u01/3, 
          tmp10=h10-u10/3, tmp11=h11-u11/3;

    // Then calculate (Mb^-1 Mh) H ~(Mb^-1 Mh). (3x24 more flops)
    return Mat<4,4,Vec3P>
       (  h00,       h00+w00/3,          h01-w01/3,       h01,
         tmp00,  tmp00+w00/3+t00/9,  tmp01-w01/3-t01/9,  tmp01,
         tmp10,  tmp10+w10/3-t10/9,  tmp11-w11/3+t11/9,  tmp11,
          h10,       h10+w10/3,          h11-w11/3,       h11  );
}
/**@}**/

private:
Mat<4,4,Vec3P> B;   // 16 Bezier control points; see above for definition
};



} // namespace SimTK

#endif // SimTK_SIMMATH_GEO_BICUBIC_BEZIER_PATCH_H_
