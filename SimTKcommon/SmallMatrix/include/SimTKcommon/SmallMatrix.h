#ifndef SimTK_SIMMATRIX_SMALLMATRIX_H_
#define SimTK_SIMMATRIX_SMALLMATRIX_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-13 Stanford University and the Authors.        *
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
/** @defgroup MatVecTypedefs     Predefined typedefs
@ingroup MatVecUtilities

These typedefs provide convenient synonyms for common matrix and vector types.
Note that the typedef name may be used interchangeably with the fully
templatized names; they represent exactly the same type.

The simplest names are for types whose elements are of the compile-time default
precision type @ref SimTK::Real "Real" which is typically \c double but
can be changed to \c float at compile time. There are also names prefixed with
a lowercase "f" that always use float regardless of the default precision.

Note that there are more template parameters than are specified here;
these typedefs are using default values for them. The missing parameters
specify the spacing between elements; these typedefs always refer to types for
which the elements are packed in memory. See Vec, Row, Mat, SymMat for more
information.

For fixed size vectors and matrices, only the 2-, 3-, and 4-element sizes
are commonly used. However we provide typedefs for sizes up to 9 in case they
are needed. For larger sizes, just use the explicit templatized forms.
**/
/**@{**/
/** This is the most common 2D vector type: a column of 2 Real values stored
consecutively in memory (packed). **/
typedef Vec<2> Vec2;
/** This is the most common 3D vector type: a column of 3 Real values stored
consecutively in memory (packed). **/
typedef Vec<3> Vec3;
/** This is the most common 4D vector type: a column of 4 Real values stored
consecutively in memory (packed). **/
typedef Vec<4> Vec4;


/** This is the most common 2x2 matrix type: two packed columns of 2 Real
values each. The columns have type \c Vec2 but rows have a stride of 2 so the
row type is \c Row<2,Real,2>, \e not \c Row2. **/
typedef Mat<2, 2> Mat22;
/** This is the most common 3x3 matrix type: three packed columns of 3 Real
values each. The columns have type \c Vec3 but rows have a stride of 3 so the
row type is \c Row<3,Real,3>, \e not \c Row3. **/
typedef Mat<3, 3> Mat33;
/** This is the most common 4x4 matrix type: four packed columns of 4 Real
values each. The columns have type \c Vec4 but rows have a stride of 4 so the
row type is \c Row<4,Real,4>, \e not \c Row4. **/
typedef Mat<4, 4> Mat44;

/** A compact, 2x2 Real symmetric matrix; only 3 elements are stored. **/
typedef SymMat<2> SymMat22;
/** A compact, 3x3 Real symmetric matrix; only 6 elements are stored. **/
typedef SymMat<3> SymMat33;
/** A compact, 2x2 Real symmetric matrix; only 10 elements are stored. **/
typedef SymMat<4> SymMat44;

/** Packed, 2-element row of Real values. This is the type of a transposed Vec2
and usually does not appear explicitly in user programs. **/
typedef Row<2> Row2;
/** Packed, 3-element row of Real values. This is the type of a transposed Vec3
and usually does not appear explicitly in user programs. **/
typedef Row<3> Row3;
/** Packed, 4-element row of Real values. This is the type of a transposed Vec4
and usually does not appear explicitly in user programs. **/
typedef Row<4> Row4;
/**@}**/

/** @defgroup UncommonMatVecTypedefs    Less commonly-used typedefs
@ingroup MatVecTypedefs **/
/**@{**/

// Less-popular Vec typedefs.
typedef Vec<1> Vec1; ///< A vector of just one Real element (not too useful).
typedef Vec<5> Vec5; ///< Packed, 5-element vector of Real values.
typedef Vec<6> Vec6; ///< Packed, 6-element vector of Real values.
typedef Vec<7> Vec7; ///< Packed, 7-element vector of Real values.
typedef Vec<8> Vec8; ///< Packed, 8-element vector of Real values.
typedef Vec<9> Vec9; ///< Packed, 9-element vector of Real values.

// Less-popular Mat typedefs.
typedef Mat<1,1> Mat11; ///< 1x1 Real matrix, that is, a scalar.
typedef Mat<1,2> Mat12; ///< 1x2 Real row matrix.
typedef Mat<1,3> Mat13; ///< 1x3 Real row matrix.
typedef Mat<1,4> Mat14; ///< 1x4 Real row matrix.
typedef Mat<1,5> Mat15; ///< 1x5 Real row matrix.
typedef Mat<1,6> Mat16; ///< 1x6 Real row matrix.
typedef Mat<1,7> Mat17; ///< 1x7 Real row matrix.
typedef Mat<1,8> Mat18; ///< 1x8 Real row matrix.
typedef Mat<1,9> Mat19; ///< 1x9 Real row matrix.

typedef Mat<2,1> Mat21; ///< 2x1 Real column matrix.
typedef Mat<2,3> Mat23; ///< 2x3 Real matrix, packed by columns.
typedef Mat<2,4> Mat24; ///< 2x4 Real matrix, packed by columns.
typedef Mat<2,5> Mat25; ///< 2x5 Real matrix, packed by columns.
typedef Mat<2,6> Mat26; ///< 2x6 Real matrix, packed by columns.
typedef Mat<2,7> Mat27; ///< 2x7 Real matrix, packed by columns.
typedef Mat<2,8> Mat28; ///< 2x8 Real matrix, packed by columns.
typedef Mat<2,9> Mat29; ///< 2x9 Real matrix, packed by columns.

typedef Mat<3,1> Mat31; ///< 3x1 Real column matrix.
typedef Mat<3,2> Mat32; ///< 3x2 Real matrix, packed by columns.
typedef Mat<3,4> Mat34; ///< 3x4 Real matrix, packed by columns.
typedef Mat<3,5> Mat35; ///< 3x5 Real matrix, packed by columns.
typedef Mat<3,6> Mat36; ///< 3x6 Real matrix, packed by columns.
typedef Mat<3,7> Mat37; ///< 3x7 Real matrix, packed by columns.
typedef Mat<3,8> Mat38; ///< 3x8 Real matrix, packed by columns.
typedef Mat<3,9> Mat39; ///< 3x9 Real matrix, packed by columns.

typedef Mat<4,1> Mat41; ///< 4x1 Real column matrix.
typedef Mat<4,2> Mat42; ///< 4x2 Real matrix, packed by columns.
typedef Mat<4,3> Mat43; ///< 4x3 Real matrix, packed by columns.
typedef Mat<4,5> Mat45; ///< 4x5 Real matrix, packed by columns.
typedef Mat<4,6> Mat46; ///< 4x6 Real matrix, packed by columns.
typedef Mat<4,7> Mat47; ///< 4x7 Real matrix, packed by columns.
typedef Mat<4,8> Mat48; ///< 4x8 Real matrix, packed by columns.
typedef Mat<4,9> Mat49; ///< 4x9 Real matrix, packed by columns.

typedef Mat<5,1> Mat51; ///< 5x1 Real column matrix.
typedef Mat<5,2> Mat52; ///< 5x2 Real matrix, packed by columns.
typedef Mat<5,3> Mat53; ///< 5x3 Real matrix, packed by columns.
typedef Mat<5,4> Mat54; ///< 5x4 Real matrix, packed by columns.
typedef Mat<5,5> Mat55; ///< 5x5 Real matrix, packed by columns.
typedef Mat<5,6> Mat56; ///< 5x6 Real matrix, packed by columns.
typedef Mat<5,7> Mat57; ///< 5x7 Real matrix, packed by columns.
typedef Mat<5,8> Mat58; ///< 5x8 Real matrix, packed by columns.
typedef Mat<5,9> Mat59; ///< 5x9 Real matrix, packed by columns.

typedef Mat<6,1> Mat61; ///< 6x1 Real column matrix.
typedef Mat<6,2> Mat62; ///< 6x2 Real matrix, packed by columns.
typedef Mat<6,3> Mat63; ///< 6x3 Real matrix, packed by columns.
typedef Mat<6,4> Mat64; ///< 6x4 Real matrix, packed by columns.
typedef Mat<6,5> Mat65; ///< 6x5 Real matrix, packed by columns.
typedef Mat<6,6> Mat66; ///< 6x6 Real matrix, packed by columns.
typedef Mat<6,7> Mat67; ///< 6x7 Real matrix, packed by columns.
typedef Mat<6,8> Mat68; ///< 6x8 Real matrix, packed by columns.
typedef Mat<6,9> Mat69; ///< 6x9 Real matrix, packed by columns.

typedef Mat<7,1> Mat71; ///< 7x1 Real column matrix.
typedef Mat<7,2> Mat72; ///< 7x2 Real matrix, packed by columns.
typedef Mat<7,3> Mat73; ///< 7x3 Real matrix, packed by columns.
typedef Mat<7,4> Mat74; ///< 7x4 Real matrix, packed by columns.
typedef Mat<7,5> Mat75; ///< 7x5 Real matrix, packed by columns.
typedef Mat<7,6> Mat76; ///< 7x6 Real matrix, packed by columns.
typedef Mat<7,7> Mat77; ///< 7x7 Real matrix, packed by columns.
typedef Mat<7,8> Mat78; ///< 7x8 Real matrix, packed by columns.
typedef Mat<7,9> Mat79; ///< 7x9 Real matrix, packed by columns.

typedef Mat<8,1> Mat81; ///< 8x1 Real column matrix.
typedef Mat<8,2> Mat82; ///< 8x2 Real matrix, packed by columns.
typedef Mat<8,3> Mat83; ///< 8x3 Real matrix, packed by columns.
typedef Mat<8,4> Mat84; ///< 8x4 Real matrix, packed by columns.
typedef Mat<8,5> Mat85; ///< 8x5 Real matrix, packed by columns.
typedef Mat<8,6> Mat86; ///< 8x6 Real matrix, packed by columns.
typedef Mat<8,7> Mat87; ///< 8x7 Real matrix, packed by columns.
typedef Mat<8,8> Mat88; ///< 8x8 Real matrix, packed by columns.
typedef Mat<8,9> Mat89; ///< 8x9 Real matrix, packed by columns.

typedef Mat<9,1> Mat91; ///< 9x1 Real column matrix.
typedef Mat<9,2> Mat92; ///< 9x2 Real matrix, packed by columns.
typedef Mat<9,3> Mat93; ///< 9x3 Real matrix, packed by columns.
typedef Mat<9,4> Mat94; ///< 9x4 Real matrix, packed by columns.
typedef Mat<9,5> Mat95; ///< 9x5 Real matrix, packed by columns.
typedef Mat<9,6> Mat96; ///< 9x6 Real matrix, packed by columns.
typedef Mat<9,7> Mat97; ///< 9x7 Real matrix, packed by columns.
typedef Mat<9,8> Mat98; ///< 9x8 Real matrix, packed by columns.
typedef Mat<9,9> Mat99; ///< 9x9 Real matrix, packed by columns.

// Less-popular SymMat typedefs.
typedef SymMat<1> SymMat11; ///< 1x1 Real symmetric matrix, that is, a scalar.
typedef SymMat<5> SymMat55; ///< 5x5 compact Real symmetric matrix.
typedef SymMat<6> SymMat66; ///< 6x6 compact Real symmetric matrix.
typedef SymMat<7> SymMat77; ///< 7x7 compact Real symmetric matrix.
typedef SymMat<8> SymMat88; ///< 8x8 compact Real symmetric matrix.
typedef SymMat<9> SymMat99; ///< 9x9 compact Real symmetric matrix.

// Less-popular Row typedefs.
typedef Row<1> Row1; ///< A row vector of one Real element (not too useful).
/** Packed, 5-element row of Real values. This is the type of a transposed Vec5. **/
typedef Row<5> Row5;
/** Packed, 6-element row of Real values. This is the type of a transposed Vec6. **/
typedef Row<6> Row6;
/** Packed, 7-element row of Real values. This is the type of a transposed Vec7. **/
typedef Row<7> Row7;
/** Packed, 8-element row of Real values. This is the type of a transposed Vec8. **/
typedef Row<8> Row8;
/** Packed, 9-element row of Real values. This is the type of a transposed Vec9. **/
typedef Row<9> Row9;

// float-precision Vecs.
typedef Vec<1, float> fVec1; ///< A vector of one float element (not too useful).
typedef Vec<2, float> fVec2; ///< Packed, 2-element vector of \c float values.
typedef Vec<3, float> fVec3; ///< Packed, 3-element vector of \c float values.
typedef Vec<4, float> fVec4; ///< Packed, 4-element vector of \c float values.
typedef Vec<5, float> fVec5; ///< Packed, 5-element vector of \c float values.
typedef Vec<6, float> fVec6; ///< Packed, 6-element vector of \c float values.
typedef Vec<7, float> fVec7; ///< Packed, 7-element vector of \c float values.
typedef Vec<8, float> fVec8; ///< Packed, 8-element vector of \c float values.
typedef Vec<9, float> fVec9; ///< Packed, 9-element vector of \c float values.

// Just doing some of the popular ones for now.
typedef Mat<1,1,float> fMat11; ///< 1x1 \c float matrix, that is, a scalar.
typedef Mat<2,2,float> fMat22; ///< 2x2 \c float matrix, packed by columns.
typedef Mat<3,3,float> fMat33; ///< 3x3 \c float matrix, packed by columns.
typedef Mat<3,4,float> fMat34; ///< 3x4 \c float matrix, packed by columns.
typedef Mat<4,3,float> fMat43; ///< 4x3 \c float matrix, packed by columns.
typedef Mat<4,4,float> fMat44; ///< 4x4 \c float matrix, packed by columns.
typedef Mat<5,5,float> fMat55; ///< 5x5 \c float matrix, packed by columns.
typedef Mat<6,6,float> fMat66; ///< 6x6 \c float matrix, packed by columns.
typedef Mat<7,7,float> fMat77; ///< 7x7 \c float matrix, packed by columns.
typedef Mat<8,8,float> fMat88; ///< 8x8 \c float matrix, packed by columns.
typedef Mat<9,9,float> fMat99; ///< 9x9 \c float matrix, packed by columns.

/** A 1x1 \c float symmetric matrix, that is, a scalar. **/
typedef SymMat<1, float> fSymMat11;
typedef SymMat<2, float> fSymMat22; ///< 2x2 compact \c float symmetric matrix.
typedef SymMat<3, float> fSymMat33; ///< 3x3 compact \c float symmetric matrix.
typedef SymMat<4, float> fSymMat44; ///< 4x4 compact \c float symmetric matrix.
typedef SymMat<5, float> fSymMat55; ///< 5x5 compact \c float symmetric matrix.
typedef SymMat<6, float> fSymMat66; ///< 6x6 compact \c float symmetric matrix.
typedef SymMat<7, float> fSymMat77; ///< 7x7 compact \c float symmetric matrix.
typedef SymMat<8, float> fSymMat88; ///< 8x8 compact \c float symmetric matrix.
typedef SymMat<9, float> fSymMat99; ///< 9x9 compact \c float symmetric matrix.

// float-precision Rows.
typedef Row<1, float> fRow1; ///< A row vector of one float element (not too useful).
typedef Row<2, float> fRow2; ///< Packed, 2-element row vector of \c float values.
typedef Row<3, float> fRow3; ///< Packed, 3-element row vector of \c float values.
typedef Row<4, float> fRow4; ///< Packed, 4-element row vector of \c float values.
typedef Row<5, float> fRow5; ///< Packed, 5-element row vector of \c float values.
typedef Row<6, float> fRow6; ///< Packed, 6-element row vector of \c float values.
typedef Row<7, float> fRow7; ///< Packed, 7-element row vector of \c float values.
typedef Row<8, float> fRow8; ///< Packed, 8-element row vector of \c float values.
typedef Row<9, float> fRow9; ///< Packed, 9-element row vector of \c float values.
/**@}**/

} //namespace SimTK


#endif //SimTK_SIMMATRIX_SMALLMATRIX_H_
