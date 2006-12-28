#ifndef SimTK_SIMMATRIX_SMALLMATRIX_H_
#define SimTK_SIMMATRIX_SMALLMATRIX_H_

/* Portions copyright (c) 2005-6 Stanford University and Michael Sherman.
 * Contributors:
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

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

#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/internal/Scalar.h"

// Forward declarations.
namespace SimTK {

template <int M, class E=Real, int STRIDE=1>              class Vec;
template <int N, class E=Real, int STRIDE=1>              class Row; 
template <int M, int N, class E=Real, int CS=M, int RS=1> class Mat; // col & row spacing
template <int M, class E=Real, int RS=1>                  class SymMat;

} // namespace SimTK

#include "SimTKcommon/internal/ResultType.h"
#include "SimTKcommon/internal/Vec.h"
#include "SimTKcommon/internal/Row.h"
#include "SimTKcommon/internal/Mat.h"
#include "SimTKcommon/internal/SymMat.h"
#include "SimTKcommon/internal/SmallMatrixMixed.h"

// Friendly abbreviations.
namespace SimTK {

typedef Vec<2> Vec2;
typedef Vec<3> Vec3;
typedef Vec<4> Vec4;
typedef Vec<5> Vec5;
typedef Vec<6> Vec6;
typedef Vec<7> Vec7;

typedef Row<2> Row2;
typedef Row<3> Row3;
typedef Row<4> Row4;
typedef Row<5> Row5;
typedef Row<6> Row6;
typedef Row<7> Row7;

typedef Mat<2,2> Mat22;
typedef Mat<3,3> Mat33;
typedef Mat<3,4> Mat34;
typedef Mat<4,3> Mat43;
typedef Mat<4,4> Mat44;
typedef Mat<5,5> Mat55;
typedef Mat<6,6> Mat66;

typedef SymMat<2> SymMat22;
typedef SymMat<3> SymMat33;
typedef SymMat<4> SymMat44;
typedef SymMat<5> SymMat55;
typedef SymMat<6> SymMat66;

} //namespace SimTK


#endif //SimTK_SIMMATRIX_SMALLMATRIX_H_
