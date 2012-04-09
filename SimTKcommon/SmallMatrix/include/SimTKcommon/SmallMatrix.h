#ifndef SimTK_SIMMATRIX_SMALLMATRIX_H_
#define SimTK_SIMMATRIX_SMALLMATRIX_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmatrix(tm)                       *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-7 Stanford University and the Authors.         *
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

/**@file
 * This file is the user-includeable header to be included in
 * user programs to provide fixed-length Vec & Mat classes. The
 * included internal headers declare templates needed for zero-overhead handling of small
 * vectors and matrices. The idea is to convey all needed information,
 * including the size, at compile time through the templatized types.
 * These classes have well-defined minimal storage layouts and can
 * be easily interconverted to/from arrays of scalar types, with zero 
 * conversion overhead.
 *
 * There are three generic types required: a column vector Vec, a row 
 * vector (a.k.a. "covector") Row, and a matrix Mat. Real and Complex
 * elements at single, double, and quad precision are supported.
 * Almost all operations are inline -- this package is intended to
 * provide performance as good as one could achieve by special-purpose C code.
 *
 * When the element type is a "basic numerical type" then the resulting
 * Vec, Row or Mat is also a basic numerical type, so these can
 * be built up with structured elements so that Vec< 2,Vec<3> > for
 * example is a 2-element vector whose elements are Vec<3>'s.
 *
 * These classes are designed so that many operations that look like
 * computations are free. These include: negation, conjugation, extraction
 * of real and imaginary parts from complex numerical types, and
 * hermitian or positional transpose. These are performed
 * by type casting rather than computation. That is, these operations
 * can be viewed as a change in perspective rather than an actual
 * computation.
 *
 * To achieve zero overhead, separate types are needed to convey
 * special matrix structure, although appropriate interconversions
 * are defined. So symmetric matrices are SymMat's and diagonal
 * matrices are DiagMat's. Note that symmetric matrices are more
 * properly called Hermitian in that the reflected elements are
 * conjugates of one another.
 *
 * Supported operations
 *
 * Class member functions for Vec<M>:
 *   Unary: -, +, ~ (Hermitian transpose)
 *   Assignment: copy assignment of identical object
 *               elementwise assignment from any M-row Vec or Mat
 *               scalar assignment to each element
 *   TODO
 *
 * Global operations involving only Vecs or Vec and scalar,
 * yielding a vector result. (Row has the same set.)
 *   v+s, s+v, v+v
 *   v-s, s-v, v-v
 *   v*s, s*v  (can't multiply a vector by another vector)
 *   v/s
 *   v==s, s==v, v==v (same for !=)
 *
 * Global operations mixing Row, Vec, Mat
 *   s=r*v (dot product)
 *   m=v*r (outer product)
 *   v=m*v, r=r*m
 *   m=m*m (with compatible dimensions)
 * 
 */

#include "SimTKcommon/Scalar.h"
#include "SimTKcommon/TemplatizedLapack.h"


#include "SimTKcommon/internal/ResultType.h"
#include "SimTKcommon/internal/Vec.h"
#include "SimTKcommon/internal/Row.h"
#include "SimTKcommon/internal/Mat.h"
#include "SimTKcommon/internal/SymMat.h"
#include "SimTKcommon/internal/SmallMatrixMixed.h"

// Friendly abbreviations.
namespace SimTK {

typedef Vec<1> Vec1; // just a scalar
typedef Vec<2> Vec2;
typedef Vec<3> Vec3;
typedef Vec<4> Vec4;
typedef Vec<5> Vec5;
typedef Vec<6> Vec6;
typedef Vec<7> Vec7;
typedef Vec<8> Vec8;
typedef Vec<9> Vec9;

typedef Vec<1,float> fVec1; // just a scalar
typedef Vec<2,float> fVec2;
typedef Vec<3,float> fVec3;
typedef Vec<4,float> fVec4;
typedef Vec<5,float> fVec5;
typedef Vec<6,float> fVec6;
typedef Vec<7,float> fVec7;
typedef Vec<8,float> fVec8;
typedef Vec<9,float> fVec9;

typedef Row<1> Row1; // just a scalar
typedef Row<2> Row2;
typedef Row<3> Row3;
typedef Row<4> Row4;
typedef Row<5> Row5;
typedef Row<6> Row6;
typedef Row<7> Row7;
typedef Row<8> Row8;
typedef Row<9> Row9;

typedef Row<1,float> fRow1; // just a scalar
typedef Row<2,float> fRow2;
typedef Row<3,float> fRow3;
typedef Row<4,float> fRow4;
typedef Row<5,float> fRow5;
typedef Row<6,float> fRow6;
typedef Row<7,float> fRow7;
typedef Row<8,float> fRow8;
typedef Row<9,float> fRow9;

typedef SymMat<1> SymMat11; // just a scalar
typedef SymMat<2> SymMat22;
typedef SymMat<3> SymMat33;
typedef SymMat<4> SymMat44;
typedef SymMat<5> SymMat55;
typedef SymMat<6> SymMat66;
typedef SymMat<7> SymMat77;
typedef SymMat<8> SymMat88;
typedef SymMat<9> SymMat99;

typedef SymMat<1,float> fSymMat11; // just a scalar
typedef SymMat<2,float> fSymMat22;
typedef SymMat<3,float> fSymMat33;
typedef SymMat<4,float> fSymMat44;
typedef SymMat<5,float> fSymMat55;
typedef SymMat<6,float> fSymMat66;
typedef SymMat<7,float> fSymMat77;
typedef SymMat<8,float> fSymMat88;
typedef SymMat<9,float> fSymMat99;

typedef Mat<1,1> Mat11; // This is just a scalar
typedef Mat<1,2> Mat12; // The rest here are just single Rows
typedef Mat<1,3> Mat13;
typedef Mat<1,4> Mat14;
typedef Mat<1,5> Mat15;
typedef Mat<1,6> Mat16;
typedef Mat<1,7> Mat17;
typedef Mat<1,8> Mat18;
typedef Mat<1,9> Mat19;

typedef Mat<2,1> Mat21; // Mats with 2 rows
typedef Mat<2,2> Mat22;
typedef Mat<2,3> Mat23;
typedef Mat<2,4> Mat24;
typedef Mat<2,5> Mat25;
typedef Mat<2,6> Mat26;
typedef Mat<2,7> Mat27;
typedef Mat<2,8> Mat28;
typedef Mat<2,9> Mat29;

typedef Mat<3,1> Mat31; // Mats with 3 rows
typedef Mat<3,2> Mat32;
typedef Mat<3,3> Mat33;
typedef Mat<3,4> Mat34;
typedef Mat<3,5> Mat35;
typedef Mat<3,6> Mat36;
typedef Mat<3,7> Mat37;
typedef Mat<3,8> Mat38;
typedef Mat<3,9> Mat39;

typedef Mat<4,1> Mat41; // Mats with 4 rows
typedef Mat<4,2> Mat42;
typedef Mat<4,3> Mat43;
typedef Mat<4,4> Mat44;
typedef Mat<4,5> Mat45;
typedef Mat<4,6> Mat46;
typedef Mat<4,7> Mat47;
typedef Mat<4,8> Mat48;
typedef Mat<4,9> Mat49;

typedef Mat<5,1> Mat51; // Mats with 5 rows
typedef Mat<5,2> Mat52;
typedef Mat<5,3> Mat53;
typedef Mat<5,4> Mat54;
typedef Mat<5,5> Mat55;
typedef Mat<5,6> Mat56;
typedef Mat<5,7> Mat57;
typedef Mat<5,8> Mat58;
typedef Mat<5,9> Mat59;

typedef Mat<6,1> Mat61; // Mats with 6 rows
typedef Mat<6,2> Mat62;
typedef Mat<6,3> Mat63;
typedef Mat<6,4> Mat64;
typedef Mat<6,5> Mat65;
typedef Mat<6,6> Mat66;
typedef Mat<6,7> Mat67;
typedef Mat<6,8> Mat68;
typedef Mat<6,9> Mat69;

typedef Mat<7,1> Mat71; // Mats with 7 rows
typedef Mat<7,2> Mat72;
typedef Mat<7,3> Mat73;
typedef Mat<7,4> Mat74;
typedef Mat<7,5> Mat75;
typedef Mat<7,6> Mat76;
typedef Mat<7,7> Mat77;
typedef Mat<7,8> Mat78;
typedef Mat<7,9> Mat79;

typedef Mat<8,1> Mat81; // Mats with 8 rows
typedef Mat<8,2> Mat82;
typedef Mat<8,3> Mat83;
typedef Mat<8,4> Mat84;
typedef Mat<8,5> Mat85;
typedef Mat<8,6> Mat86;
typedef Mat<8,7> Mat87;
typedef Mat<8,8> Mat88;
typedef Mat<8,9> Mat89;

typedef Mat<9,1> Mat91; // Mats with 9 rows
typedef Mat<9,2> Mat92;
typedef Mat<9,3> Mat93;
typedef Mat<9,4> Mat94;
typedef Mat<9,5> Mat95;
typedef Mat<9,6> Mat96;
typedef Mat<9,7> Mat97;
typedef Mat<9,8> Mat98;
typedef Mat<9,9> Mat99;

// Just doing some of the popular ones for now.
typedef Mat<1,1,float> fMat11;
typedef Mat<2,2,float> fMat22;
typedef Mat<3,3,float> fMat33;
typedef Mat<3,4,float> fMat34;
typedef Mat<4,3,float> fMat43;
typedef Mat<4,4,float> fMat44;
typedef Mat<5,5,float> fMat55;
typedef Mat<6,6,float> fMat66;
typedef Mat<7,7,float> fMat77;
typedef Mat<8,8,float> fMat88;
typedef Mat<9,9,float> fMat99;


} //namespace SimTK


#endif //SimTK_SIMMATRIX_SMALLMATRIX_H_
