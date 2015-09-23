#ifndef SimTK_SIMMATH_GEO_CUBIC_HERMITE_CURVE_H_
#define SimTK_SIMMATH_GEO_CUBIC_HERMITE_CURVE_H_

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
Provides primitive operations for a single cubic Hermite curve using either
single or double precision. **/

#include "SimTKcommon.h"
#include "simmath/internal/common.h"
#include "simmath/internal/Geo.h"

#include <cassert>
#include <cmath>
#include <algorithm>

namespace SimTK {


//==============================================================================
//                         GEO CUBIC HERMITE CURVE
//==============================================================================
/** A primitive useful for computations involving a single cubic Hermite
curve segment in algebraic or geometric (Hermite) form. Objects of this class
contain the algebraic coefficients because most operations are more efficient
in that form, but methods are provided for easy conversion to or from Hermite
form and for working directly with the Hermite form.

Note that a cubic Hermite spline (made up of multiple segments) would not
necessarily be composed of these because they can be constructed more compactly
with shared end points. However, the primitive and inline methods here can be
used for fast curve segment computations.

<h3>Theory</h3>
The primary reference for this implementation is the book "Geometric Modeling,
3rd ed." by Michael E. Mortenson, Industrial Press 2006, chapter 3. We follow
Mortenson's notation here (with a few exceptions) and equation numbers are
from the text. We're using h's for the Hermite coefficients rather than
Mortenson's b's to avoid confusion with Bezier control points B, so that we
can use A for algebraic, H for Hermite, and B for Bezier coefficient matrices.

The curve is parameterized by a scalar u in [0..1], such that points on the
curve, and their derivatives with respect to u are given by <pre>
    P(u)   = a3 u^3 + a2 u^2 +   a1 u +   a0                               (3.2)
    Pu(u)  =        3 a3 u^2 + 2 a2 u +   a1
    Puu(u) =                   6 a3 u + 2 a2
    Puuu(u) =                           6 a3
</pre> where Pu=dP/du, Puu=d2P/du2, Puuu=d3P/du3. Note that all higher
derivatives are zero for a cubic. The 3-vectors ai are the algebraic
coefficients, and that is the algebraic form of this cubic curve. The Hermite,
or geometric, form is parameterized by position and tangent at the end points,
given by <pre>
    h0=P(0)=a0, h1=P(1)=a3+a2+a1+a0, hu0=Pu(0)=a1, hu1=Pu(1)=3*a3+2*a2+a1
</pre>
These define the Hermite coefficients h0, h1, hu0, hu1. In this form the
curve's points are <pre>
    P(u)    =    F1(u) h0 +    F2(u) h1 +    F3(u) hu0 +    F4(u) hu1      (3.6)
    Pu(u)   =   F1u(u) h0 +   F2u(u) h1 +   F3u(u) hu0 +   F4u(u) hu1
    Puu(u)  =  F1uu(u) h0 +  F2uu(u) h1 +  F3uu(u) hu0 +  F4uu(u) hu1
    Puuu(u) = F1uuu(u) h0 + F2uuu(u) h1 + F3uuu(u) hu0 + F4uuu(u) hu1
</pre> where the Fi's are the (scalar) Hermite basis functions given by
<pre>
    F1(u) =  2 u^3 - 3 u^2     + 1
    F2(u) = -2 u^3 + 3 u^2                                                 (3.4)
    F3(u) =    u^3 - 2 u^2 + u
    F4(u) =    u^3 -   u^2

    F1u(u) =  6 u^2 - 6 u       F1uu(u) =  12 u - 6   F1uuu(u) =  12       (3.7)
    F2u(u) = -6 u^2 + 6 u       F2uu(u) = -12 u + 6   F2uuu(u) = -12       (3.8)
    F3u(u) =  3 u^2 - 4 u + 1   F3uu(u) =   6 u - 4   F3uuu(u) = 6
    F4u(u) =  3 u^2 - 2 u       F4uu(u) =   6 u - 2   F4uuu(u) = 6
</pre>
In matrix notation, let Fh=[F1 F2 F3 F4], and U=[u^3 u^2 u 1]. Then
    Fh = U Mh
where Mh, the Hermite basis transformation matrix, and its inverse are:
<pre>
         [ 2 -2  1  1 ]             [ 0  0  0  1 ]
    Mh = [-3  3 -2 -1 ]   inv(Mh) = [ 1  1  1  1 ]                        (3.18)
         [ 0  0  1  0 ]             [ 0  0  1  0 ]                        (3.23)
         [ 1  0  0  0 ]             [ 3  2  1  0 ]
</pre>
(Mortenson calls this matrix "Mf".) Now we can write the algebraic and Hermite
forms in matrix notation. Let A=~[a3 a2 a1 a0], H=~[h0 h1 hu0 hu1]. We have
<pre>
    P(u) = U A
         = U Mh H
       A = Mh H                                                        (3.20-22)
       H = inv(Mh) A
</pre>
where the last two equations show how to convert between the algebraic and
geometric forms. Note that while U, Fh, and Mh are ordinary matrices, A and H
are hypermatrices since their elements are 3-vectors.

Because of the sparsity of the matrices and the many common subexpressions
above, it saves a considerable amount of computation to work out the necessary
products by hand, and this implementation does that. For example, to find the
algebraic coefficients A given the Hermite coefficients H the matrix-vector
multiply Mh*H would take 3x28=84 flops, while the hand-worked version is:
<pre>
       [ a3 ]   [ 2 (h0 - h1) +   hu0 + hu1 ]
   A = [ a2 ] = [-3 (h0 - h1) - 2 hu0 - hu1 ]
       [ a1 ]   [           hu0             ]
       [ a0 ]   [           h0              ]
</pre> which instead takes 3x8=24 flops, 3.5X faster.

@see CubicBezierCurve_, BicubicHermitePatch_, BicubicBezierPatch_
**/
template <class P>
class Geo::CubicHermiteCurve_ {
typedef P               RealP;
typedef Vec<3,RealP>    Vec3P;

public:
/** Construct an uninitialized curve; coefficients will be garbage. **/
CubicHermiteCurve_() {}
/** Construct a cubic Hermite curve using the given algebraic coefficients
A=[a3 a2 a1 a0], such that points on the curve are P(u)=sum_i(ai*u^i). If you
have Hermite coefficients H, convert them to algebraic using static method
calcAFromH(). **/
explicit CubicHermiteCurve_(const Vec<4,Vec3P>& A) : A(A) {}
/** Return a reference to the algebraic coefficients A=[a3 a2 a1 a0] that are
stored in this object. **/
const Vec<4,Vec3P>& getAlgebraicCoefficients() const {return A;}
/** Calculate the Hermite coefficients H=[h0 h1 hu0 hu1] from the stored
algebraic coefficients. Here h0=P(0), h1=P(1), hu0=(dP/du)(0), hu1=(dP/du)(1)
where P(u) is the curve evaluation function. Cost is 21 flops. **/
Vec<4,Vec3P> calcGeometricCoefficients() const {return calcHFromA(A);}
/** Evaluate a point P(u) on this curve given a value for parameter u in [0,1].
Values outside this range are permitted but do not lie on the curve segment.
Cost is 20 flops. **/
Vec3P evalP(RealP u) const {return evalPUsingA(A,u);}
/** Evaluate the tangent Pu=dP/du on this curve given a value for parameter u
in [0,1]. Values outside this range are permitted but do not lie on the curve
segment. Cost is 15 flops. **/
Vec3P evalPu(RealP u) const {return evalPuUsingA(A,u);}
/** Evaluate the second derivative Puu=d2P/du2 on this curve given a value for
parameter u in [0,1]. Values outside this range are permitted but do not lie on
the curve segment. Cost is 10 flops. **/
Vec3P evalPuu(RealP u) const {return evalPuuUsingA(A,u);}
/** Evaluate the third derivative Puuu=d3P/du3 on this curve. Parameter u is
ignored here since the 3rd derivative of a cubic curve is a constant. Cost is
3 flops. **/
Vec3P evalPuuu(RealP u) const {return evalPuuuUsingA(A,u);}

/**@name                      Utility methods
These static methods provide operations useful for working with cubic Hermite
curves. **/
/**@{**/
/** Return the row vector U=[u^3 u^2 u 1]. Cost is 2 flops. **/
static Row<4,P> calcU(RealP u) {
    const RealP u2 = u*u;
    return Row<4,P>(u*u2, u2, u, 1);
}

/** Calculate the Hermite basis functions Fh=[F1..F4] for a given value of the
parameter u. This is an optimized calculation of U*Mh, taking 10 flops. **/
static Row<4,P> calcFh(RealP u) {
    const RealP u2 = u*u, u3 = u*u2;
    const RealP t = 3*u2 - 2*u3;
    return Row<4,P>(1-t, t, u3-2*u2+u, u3-u2);
}

/** Calculate first derivatives Fhu=[F1u..F4u] of the Hermite basis functions
for a given value of the parameter u. Cost is 10 flops. **/
static Row<4,P> calcFhu(RealP u) {
    const RealP u2 = u*u, u23=3*u2;
    const RealP dt = 6*(u-u2);
    return Row<4,P>(-dt, dt, u23-4*u+1, u23-2*u);
}

/** Calculate second derivatives Fhuu=[F1uu..F4uu] of the Hermite basis
functions for a given value of the parameter u. Cost is 6 flops. **/
static Row<4,P> calcFhuu(RealP u) {
    const RealP u6  = 6*u;
    const RealP ddt = 6 - 12*u;
    return Row<4,P>(-ddt, ddt, u6-4, u6-2);
}

/** Calculate third derivatives Fhuuu=[F1uuu..F4uuu] of the Hermite basis
functions for a given value of the parameter u. For a cubic curve the third
derivative is a constants so the cost is 0 flops. **/
static Row<4,P> calcFhuuu(RealP u) {
    return Row<4,P>(12, -12, 6, 6);
}

/** Given the Hermite coefficients H=~[h0 h1 hu0 hu1], return the algebraic
coefficients A=~[a3 a2 a1 a0]. All coefficients are 3-vectors. Cost is 24
flops. **/
static Vec<4,Vec3P> calcAFromH(const Vec<4,Vec3P>& H) {
    const Vec3P& h0= H[0]; const Vec3P& h1= H[1];    // aliases for beauty
    const Vec3P& hu0=H[2]; const Vec3P& hu1=H[3];
    const Vec3P h01 = h0 - h1;
    return Vec<4,Vec3P>( 2*h01 +   hu0 + hu1,
                        -3*h01 - 2*hu0 - hu1,
                         hu0,
                         h0 );
}

/** Given the algebraic coefficients A=~[a3 a2 a1 a0], return the Hermite
coefficients H=~[h0 h1 hu0 hu1]. All coefficients are 3-vectors. Cost is 21
flops. **/
static Vec<4,Vec3P> calcHFromA(const Vec<4,Vec3P>& A) {
    const Vec3P& a3=A[0]; const Vec3P& a2=A[1];    // aliases for beauty
    const Vec3P& a1=A[2]; const Vec3P& a0=A[3];
    return Vec<4,Vec3P>(a0, a3+a2+a1+a0, a1, 3*a3+2*a2+a1);
}

/** Given algebraic coefficients A and a value for the curve parameter u, return
the point P(u) at that location. Cost is 20 flops. **/
static Vec3P evalPUsingA(const Vec<4,Vec3P>& A, RealP u) {
    const RealP u2 = u*u, u3 = u*u2;
    return u3*A[0] + u2*A[1] + u*A[2] + A[3]; // (vectors)
}
/** Given algebraic coefficients A and a value for the curve parameter u, return
the first derivative Pu(u)=dP/du at that location. Cost is 15 flops. **/
static Vec3P evalPuUsingA(const Vec<4,Vec3P>& A, RealP u) {
    return (3*u*u)*A[0] + (2*u)*A[1] + A[2];
}
/** Given algebraic coefficients A and a value for the curve parameter u, return
the second derivative Puu(u)=d2P/du2 at that location. Cost is 10 flops. **/
static Vec3P evalPuuUsingA(const Vec<4,Vec3P>& A, RealP u) {
    return (6*u)*A[0] + 2*A[1];
}
/** Given algebraic coefficients A and a value for the curve parameter u, return
the third derivative Puuu(u)=d3P/du3 at that location. The parameter u is
ignored since the third derivative of a cubic is just a constant. Cost is
3 flops. **/
static Vec3P evalPuuuUsingA(const Vec<4,Vec3P>& A, RealP u) {
    return 6*A[0];
}

/** Given Hermite coefficients H and a value for the curve parameter u, return
the point P(u) at that location. Cost is 31 flops. Note that if you need to
do this for the same curve more than twice, it is cheaper to convert to
algebraic form using calcAFromH() (24 flops) and then evaluate using A
(20 flops). **/
static Vec3P evalPUsingH(const Vec<4,Vec3P>& H, RealP u) {
    return calcFh(u)*H; // 10 + 3*7 = 31 flops
}
/** Given Hermite coefficients H and a value for the curve parameter u, return
the first derivative Pu(u)=dP/du at that location. Cost is 31 flops. Note that
if you need to do this for the same curve more than once, it is cheaper to
convert to algebraic form using calcAFromH() (24 flops) and then evaluate using
A (15 flops). **/
static Vec3P evalPuUsingH(const Vec<4,Vec3P>& H, RealP u) {
    return calcFhu(u)*H; // 10 + 3*7 = 31 flops
}
/** Given Hermite coefficients H and a value for the curve parameter u, return
the second derivative Puu(u)=d2P/du2 at that location. Cost is 27 flops. Note
that if you need to do this for the same curve more than once, it is cheaper to
convert to algebraic form using calcAFromH() (24 flops) and then evaluate using
A (10 flops). **/
static Vec3P evalPuuUsingH(const Vec<4,Vec3P>& H, RealP u) {
    return calcFhuu(u)*H; // 6 + 3*7 = 27 flops
}
/** Given Hermite coefficients H and a value for the curve parameter u, return
the third derivative Puuu(u)=d3P/du3 at that location. Cost is 21 flops. Note
that if you need to do this for the same curve more than once, it is cheaper to
convert to algebraic form using calcAFromH() (24 flops) and then evaluate using
A (3 flops). **/
static Vec3P evalPuuuUsingH(const Vec<4,Vec3P>& H, RealP u) {
    return calcFhuuu(u)*H; // 0 + 3*7 = 21 flops
}

/** Obtain the Hermite basis matrix Mh explicitly. This is mostly useful for
testing since specialized routines can save a lot of CPU time over working
directly in matrix form. This is a constant matrix so there is no computation
cost. **/
static Mat<4,4,P> getMh() {
    return Mat<4,4,P>( 2, -2,  1,  1,
                      -3,  3, -2, -1,
                       0,  0,  1,  0,
                       1,  0,  0,  0 );
}

/** Form the product of the Hermite basis matrix Mh and a 4-vector, exploiting
the structure of Mh (which is not symmetric). Cost is 8 flops. **/
static Vec<4,P> multiplyByMh(const Vec<4,P>& v) {
    const RealP v0=v[0], v1=v[1], v2=v[2], v3=v[3];
    const RealP v01 = v0-v1;
    return Vec<4,P>(2*v01+v2+v3, -3*v01-2*v2-v3, v2, v0);
}

/** Obtain the inverse inv(Mh) of the Hermite basis matrix explicitly. This is
mostly useful for testing since specialized routines can save a lot of CPU time
over working directly in matrix form. This is a constant matrix so there is no
computation cost. **/
static Mat<4,4,P> getMhInv() {
    return Mat<4,4,P>( 0,  0,  0,  1,
                       1,  1,  1,  1,
                       0,  0,  1,  0,
                       3,  2,  1,  0 );
}

/** Form the product of the inverse Hermite basis matrix inv(Mh) and a
4-vector, exploiting the structure of inv(Mh) (which is not symmetric).
Cost is 7 flops. **/
static Vec<4,P> multiplyByMhInv(const Vec<4,P>& v) {
    const RealP v0=v[0], v1=v[1], v2=v[2], v3=v[3];
    return Vec<4,P>(v3, v0+v1+v2+v3, v2, 3*v0+2*v1+v2);
}
/**@}**/

private:
Vec<4,Vec3P> A;
};



} // namespace SimTK

#endif // SimTK_SIMMATH_GEO_CUBIC_HERMITE_CURVE_H_
