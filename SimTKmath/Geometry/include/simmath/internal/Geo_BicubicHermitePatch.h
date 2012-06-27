#ifndef SimTK_SIMMATH_GEO_BICUBIC_HERMITE_PATCH_H_
#define SimTK_SIMMATH_GEO_BICUBIC_HERMITE_PATCH_H_

/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2011-12 Stanford University and the Authors.        *
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

/** @file
Provides primitive operations for a single bicubic Hermite patch using either
single or double precision. **/

#include "SimTKcommon.h"
#include "simmath/internal/common.h"
#include "simmath/internal/Geo.h"

#include <cassert>
#include <cmath>
#include <algorithm>

namespace SimTK {


//==============================================================================
//                         GEO BICUBIC HERMITE PATCH
//==============================================================================
/** A primitive useful for computations involving a single bicubic Hermite
patch. Note that a bicubic Hermite spline surface would not necessarily be
composed of these, but could use the static methods here for patch 
computations. 

<h3>Theory</h3>
The primary reference for this implementation is the book "Geometric Modeling, 
3rd ed." by Michael E. Mortenson, Industrial Press 2006, chapter 7. We follow
Mortenson's notation here (with some name changes) and equation numbers are 
from the text. See CubicHermiteCurve_ comments for 
introductory material; here we add the code needed for a bicubic Hermite
surface, leaning on the cubic Hermite curve code whenever possible. Note that
a bicubic surface is bounded by cubic curves and every cross section is a
cubic curve. This class deals with a bicubic patch in algebraic or Hermite
(geometric) form. The algebraic form is stored since most operations are faster
when the algebraic coefficients are available. Methods are available to 
convert quickly between the algebraic and Hermite forms. See 
BicubicBezierPatch_ for additional code for conversions to/from Bezier form
and additional operations that are better performed on Bezier control points.

We use H for Hermite coefficients, rather than B, to avoid confusion with
Bezier coefficients. We call the Hermite basis matrix Mh (rather than Mf).
With that name change, the algebraic and Hermite basis matrices match
Mortenson's:
<pre>
        [ a33 a32 a31 a30 ]         [ h00  h01  w00  w01 ]    u=dp/du
    A = [ a23 a22 a21 a20 ]     H = [ h10  h11  w10  w11 ]    w=dp/dw
        [ a13 a12 a11 a10 ]         [ u00  u01  t00  t01 ]    t=d2p/dudw
        [ a03 a02 a01 a00 ]         [ u10  u11  t10  t11 ]      ("twist")
</pre>
We store those in a 4x4 hypermatrix with Vec3 elements. Note that the element
indices do \e not match the matrix indices; instead they have different 
meanings as defined by Mortenson. For the algebraic coefficients we have
aij = A[3-i][3-j].

The patch is parameterized by u,w each in
range [0..1]. A point P(u,w) on the patch can be calculated as
<pre>
    P(u,w) = U A ~W = sum_ij( aij u^i w^j )
           = U Mh H ~Mh ~W
</pre>
where U and W are row vectors U=[u^3 u^2 u 1], W=[w^3 w^2 w 1]. From that you
can see how to interconvert between the two forms:
<pre>
    A = Mh    H ~Mh
    H = Mh^-1 A ~Mh^-1
</pre>
The Mh matrix is presented explicitly in CubicHermiteCurve_; it has a very 
simple structure. This can be used to work out low cost interconversions:

**/
template <class P>
class Geo::BicubicHermitePatch_ {
typedef P               RealP;
typedef Vec<3,RealP>    Vec3P;

public:
/** Construct an uninitialized patch; control points will be garbage. **/
BicubicHermitePatch_() {}
/** Construct a bicubic Hermite patch using the given geometry matrix B. **/
explicit BicubicHermitePatch_(const Mat<4,4,Vec3P>& A) : A(A) {} 
/** Return a reference to the algebraic coefficients A that are
stored in this object. See the documentation for this class to see how the
returned matrix of coefficients is defined. **/
const Mat<4,4,Vec3P>& getAlgebraicCoefficients() const {return A;}
/** Calculate the Hermite coefficients H from the stored 
algebraic coefficients. See the documentation for this class to see how the
returned matrix of coefficients is defined. Cost is 168 flops. **/
Mat<4,4,Vec3P> calcHermiteCoefficients() const {return calcHFromA(A);}

/** Evaluate a point P(u,w) on this patch given values for the
parameters u and w in [0,1]. Values outside this range are permitted but do 
not lie on the patch. Cost is 94 flops. **/
Vec3P evalP(RealP u, RealP w) const {return evalPUsingA(A,u,w);}

/** Evaluate the tangents Pu=dP/du, Pw=dP/dw on this patch given values for the
parameters u and w in [0,1]. Values outside this range are permitted but do not
lie on the curve segment. Cost is 148 flops. **/
void evalP1(RealP u, RealP w, Vec3P& Pu, Vec3P& Pw) const 
{   return evalP1UsingA(A,u,w,Pu,Pw); }

/** Evaluate the second derivatives Puu=d2P/du2, Pww=d2P/dw2, and cross
derivative Puw=Pwu=d2P/dudw on this patch given values for the parameters u 
and w in [0,1]. Values outside this range are permitted but do not lie on the 
curve segment. Cost is 172 flops. **/
void evalP2(RealP u, RealP w, Vec3P& Puu, Vec3P& Puw, Vec3P& Pww) const 
{   evalP2UsingA(A,u,w,Puu,Puw,Pww); }

/** Evaluate the third derivatives Puuu=d3P/du3, Pwww=d3P/dw3, and cross
derivatives Puuw=Pwuu=Puwu=d3P/du2dw and Puww=Pwwu=Pwuw=d3P/dudw2 on this patch 
given values for the parameters u and w in [0,1]. Cost is 142 flops. All
higher derivatives of a cubic patch are zero. **/
void evalP3(RealP u, RealP w, Vec3P& Puuu, Vec3P& Puuw, 
                              Vec3P& Puww, Vec3P& Pwww) const 
{   evalP3UsingA(A,u,w,Puuu,Puuw,Puww,Pwww); }

/**@name                 Utility methods
These static methods work with given coefficients. **/
/**@{**/

/** Given the vector Hermite coefficients H, return the algebraic
coefficients A. All coefficients are 3-vectors. Cost is 3x70=210 flops. **/
static Mat<4,4,Vec3P> calcAFromH(const Mat<4,4,Vec3P>& H) {
    typedef const Vec3P& Coef;
    Coef h00=H(0,0), h01=H(0,1), w00=H(0,2), w01=H(0,3), 
         h10=H(1,0), h11=H(1,1), w10=H(1,2), w11=H(1,3), 
         u00=H(2,0), u01=H(2,1), t00=H(2,2), t01=H(2,3), 
         u10=H(3,0), u11=H(3,1), t10=H(3,2), t11=H(3,3); 
    // First calculate Mh*H:
    //       a   b   c   d
    //       p   q   r   s
    //      u00 u01 t00 t01
    //      h00 h01 w00 w01
    Vec3P a=2*(h00-h10)+u00+u10,   b=2*(h01-h11)+u01+u11,
          c=2*(w00-w10)+t00+t10,   d=2*(w01-w11)+t01+t11;
    Vec3P p=3*(h10-h00)-2*u00-u10, q=3*(h11-h01)-2*u01-u11,
          r=3*(w10-w00)-2*t00-t10, s=3*(w11-w01)-2*t01-t11;

    // Then calculate (Mh*H)*~Mh.
    Vec3P bma=b-a, qmp=q-p;
    return Mat<4,4,Vec3P>
       (     c+d-2*bma,            3*bma-2*c-d,        c,   a,
             r+s-2*qmp,            3*qmp-2*r-s,        r,   p,
        2*(u00-u01)+t00+t01,  3*(u01-u00)-2*t00-t01,  t00, u00,
        2*(h00-h01)+w00+w01,  3*(h01-h00)-2*w00-w01,  w00, h00  );
}

/** Given the vector algebraic coefficients A, return the Hermite coefficients
H. All coefficients are 3-vectors. Cost is 3x56=168 flops. **/
static Mat<4,4,Vec3P> calcHFromA(const Mat<4,4,Vec3P>& A) {
    typedef const Vec3P& Coef;
    Coef a33=A(0,0), a32=A(0,1), a31=A(0,2), a30=A(0,3), 
         a23=A(1,0), a22=A(1,1), a21=A(1,2), a20=A(1,3), 
         a13=A(2,0), a12=A(2,1), a11=A(2,2), a10=A(2,3), 
         a03=A(3,0), a02=A(3,1), a01=A(3,2), a00=A(3,3); 
    // First calculate Mh^-1*A:
    //      a03 a02 a01 a00
    //       a   b   c   d
    //      a13 a12 a11 a10
    //       p   q   r   s
    Vec3P a=a33+a23+a13+a03, b=a32+a22+a12+a02, // 3x12 flops
          c=a31+a21+a11+a01, d=a30+a20+a10+a00;
    Vec3P p=3*a33+2*a23+a13, q=3*a32+2*a22+a12, // 3x16 flops
          r=3*a31+2*a21+a11, s=3*a30+2*a20+a10;

    // Then calculate (Mh^-1*A)*Mh^-T (3x28 flops)
    return Mat<4,4,Vec3P>
       ( a00,  a03+a02+a01+a00,  a01,  3*a03+2*a02+a01,
          d,       a+b+c+d,       c,      3*a+2*b+c,
         a10,  a13+a12+a11+a10,  a11,  3*a13+2*a12+a11,
          s,       p+q+r+s,       r,      3*p+2*q+r     );
}

/** Given vector algebraic coefficients A and values for the curve parameters 
u and w in [0..1], return the point P(u,w)=U*A*~W at that location. Cost is 
3x30+4=94 flops. **/
static Vec3P evalPUsingA(const Mat<4,4,Vec3P>& A, RealP u, RealP w) { 
    typedef const Vec3P& Coef;
    Coef a33=A(0,0), a32=A(0,1), a31=A(0,2), a30=A(0,3), 
         a23=A(1,0), a22=A(1,1), a21=A(1,2), a20=A(1,3), 
         a13=A(2,0), a12=A(2,1), a11=A(2,2), a10=A(2,3), 
         a03=A(3,0), a02=A(3,1), a01=A(3,2), a00=A(3,3); 

    const RealP u2 = u*u, u3 = u*u2, w2 = w*w, w3 = w*w2;
    Vec3P p =   u3*(a33*w3 + a32*w2 + a31*w + a30)
              + u2*(a23*w3 + a22*w2 + a21*w + a20)
              + u *(a13*w3 + a12*w2 + a11*w + a10)
              +    (a03*w3 + a02*w2 + a01*w + a00);
    return p;
}

/** Given vector algebraic coefficients A and values for the curve parameters 
u and w in [0..1], return the tangents Pu(u,w)=dU*A*~W and Pw(u,w)=U*A*~dW at 
that location. Cost is 3x48+4=148 flops. **/
static void evalP1UsingA(const Mat<4,4,Vec3P>& A, RealP u, RealP w,
                         Vec3P& Pu, Vec3P& Pw) {
    typedef const Vec3P& Coef;
    Coef a33=A(0,0), a32=A(0,1), a31=A(0,2), a30=A(0,3), 
         a23=A(1,0), a22=A(1,1), a21=A(1,2), a20=A(1,3), 
         a13=A(2,0), a12=A(2,1), a11=A(2,2), a10=A(2,3), 
         a03=A(3,0), a02=A(3,1), a01=A(3,2); 

    const RealP u2 = u*u, u3 = u*u2, w2 = w*w, w3 = w*w2;
    Pu =   3*u2*(a33*w3 + a32*w2 + a31*w + a30)
         + 2* u*(a23*w3 + a22*w2 + a21*w + a20)
         +      (a13*w3 + a12*w2 + a11*w + a10);
    Pw =   3*w2*(u3*a33 + u2*a23 + u*a13 + a03)
         + 2* w*(u3*a32 + u2*a22 + u*a12 + a02)
         +      (u3*a31 + u2*a21 + u*a11 + a01);
}

/** Given vector algebraic coefficients A and values for the curve parameters 
u and w in [0..1], return the second derivatives Puu(u,w)=d2U*A*~W,
Puw(u,w)=dU*A*~dW and Pww(u,w)=U*A*~d2W at that location. 
Cost is 3x56+4=172 flops. **/
static void evalP2UsingA(const Mat<4,4,Vec3P>& A, RealP u, RealP w,
                         Vec3P& Puu, Vec3P& Puw, Vec3P& Pww) {
    typedef const Vec3P& Coef;
    Coef a33=A(0,0), a32=A(0,1), a31=A(0,2), a30=A(0,3), 
         a23=A(1,0), a22=A(1,1), a21=A(1,2), a20=A(1,3), 
         a13=A(2,0), a12=A(2,1), a11=A(2,2),
         a03=A(3,0), a02=A(3,1); 

    const RealP u2 = u*u, u3 = u*u2, w2 = w*w, w3 = w*w2;
    Puu =   6*u*(a33*w3 + a32*w2 + a31*w + a30)
          + 2  *(a23*w3 + a22*w2 + a21*w + a20);
    Pww =   6*w*(u3*a33 + u2*a23 + u*a13 + a03)
          + 2  *(u3*a32 + u2*a22 + u*a12 + a02);
    Puw =   3*u2*(3*a33*w2 + 2*a32*w + a31)
          + 2*u *(3*a23*w2 + 2*a22*w + a21)
          +      (3*a13*w2 + 2*a12*w + a11);
}

/** Given vector algebraic coefficients A and values for the curve parameters 
u and w in [0..1], return the third derivatives Puuu(u,w)=d3U*A*~W,
Puuw(u,w)=d2U*A*~dW, Puww(u,w)=dU*A*~d2W and Pwww(u,w)=U*A*~d3W at that 
location. Cost is 3x46+4=142 flops. **/
static void evalP3UsingA(const Mat<4,4,Vec3P>& A, RealP u, RealP w,
                         Vec3P& Puuu, Vec3P& Puuw, Vec3P& Puww, Vec3P& Pwww) {
    typedef const Vec3P& Coef;
    Coef a33=A(0,0), a32=A(0,1), a31=A(0,2), a30=A(0,3), 
         a23=A(1,0), a22=A(1,1), a21=A(1,2), a20=A(1,3), 
         a13=A(2,0), a12=A(2,1), 
         a03=A(3,0), a02=A(3,1); 

    const RealP u2 = u*u, u3 = u*u2, w2 = w*w, w3 = w*w2;
    Puuu = 6*(a33*w3 + a32*w2 + a31*w + a30);
    Pwww = 6*(u3*a33 + u2*a23 + u*a13 + a03);
    Puuw = 6*u*(3*a33*w2 + 2*a32*w + a31) 
           +    6*a23*w2 + 4*a22*w + 2*a21;
    Puww = 6*w*(3*u2*a33 + 2*u*a23 + a13) 
           +    6*u2*a32 + 4*u*a22 + 2*a12;
}
/**@}**/

private:
Mat<4,4,Vec3P> A;   // algebraic coefficients
};



} // namespace SimTK

#endif // SimTK_SIMMATH_GEO_BICUBIC_HERMITE_PATCH_H_
