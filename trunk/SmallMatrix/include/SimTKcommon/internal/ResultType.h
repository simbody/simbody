#ifndef SimTK_SIMMATRIX_RESULT_TYPE_H_
#define SimTK_SIMMATRIX_RESULT_TYPE_H_

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
 * This file defines a template class with specializations that helps
 * determine the result type when arithmetic operators are applied
 * to Composite Numerical Types.
 */

#include "SimTKcommon/internal/common.h"

namespace SimTK {


/**
 * The classes here can be thought of as compile-time type calculators
 * which can compute the appropriate return type when the binary
 * arithmetic operators +,-,*,/ are applied to Composite Numerical
 * Type (CNT) arguments. We define these template classes:
 *      AddCNTs       MulCNTs
 *      SubCNTs       DvdCNTs
 * and specialize them for each distinct combination of arguments.
 * Each has a typedef member "Type" which gives the result type,
 * and a static method "perform()" which delegates the operation to
 * the appropriate method of one of the arguments.

 * The general template class must be sufficiently parameterized that ALL 
 * combinations of CNT arguments can be viewed as special cases.
 * That permits us to specialize properly here, and then individual
 * CNTs base their Result classes on this one, without the need
 * for specialization (which doesn't work well for internal
 * classes on some compilers).
 *
 * The argument attributes that we must be able to specialize on
 * are: shape, type (Row,Vec,Mat,SymMat) and "depth", meaning whether
 * an argument is a scalar (depth 0), composite with scalar elements 
 * (depth 1), or composite with composite elements (depth 2). We do
 * not implement the most general form at all; all meaningful
 * combinations of CNT shapes must be specialized. We distinguish
 * conforming arguments from non-conforming. The idea is to determine
 * whether we are doing a matrix/matrix operation, or a matrix/element
 * operation, and then in the latter case which of the two arguments is
 * to be treated as the matrix.
 *
 * Here are the rules:
 *    (0) if both arguments are scalars results typing will be 
 *        handled elsewhere, so we don't worry about that case here
 *    (1) if one argument is a scalar, we are doing a matrix/element
 *        operation (and of course the other one is the matrix) unless
 *        the scalar is the left operand of a divide
 *    (2) if neither argument is scalar, and the dimensions are
 *        conforming for the particular operation being performed,
 *        then we are doing a matrix/matrix operation
 *    (3) (no scalar, not conforming) If one argument has scalar
 *        elements (depth 1) and the other has composites (depth 2),
 *        treat the depth 1 argument as though it were just a scalar
 *        and perform a matrix/element operation.
 *    (4) (no scalar, not conforming, both or neither have scalar
 *        elements) This is an error and won't compile.
 *
 * Conforming matrices
 * -------------------
 *
 * Whether two CNTs are conforming depends on the operation being
 * performed. For add and subtract, they are conforming if and only
 * if both CNTs have the same outer shape (number of rows and columns)
 * and both are rows, both are vectors, or both are matrices. (Note that
 * the Mat and SymMat CNTs are both matrices so can be mixed.) For
 * conforming multiplication the number of columns of the LHS must
 * match the number of rows on the RHS. For conforming division (that is,
 * right multiplication by the inverse) the number of columns of the LHS
 * must match the number of columns on the RHS. Here are the cases, with
 * m=Mat, sy=SymMat, v=Vec, r=Row, and s=scalar:
 *
 *    Add, Sub: conforming means LHS_ncol==RHS_ncol && LHS_nrow==RHS_nrow
 *                  m=m+m m=sy+m m=m+sy sy=sy+sy 
 *                  v=v+v  r=r+r
 *         Mul: conforming means LHS_ncol==RHS_nrow
 *                  m=m*m  m=sy*m m=m*sy m=sy*sy m=v*r (outer product)
 *                  v=m*v  v=sy*v
 *                  r=r*m  r=r*sy
 *                  s=r*v (dot product)
 *         Dvd: conforming means LHS_ncol==RHS_ncol
 *                  m=m/m  m=sy/m m=m/sy m=sy/sy 
 *                  v=m/r  v=sy/r
 *                  r=r/m  r=r/sy
 *                  r=s/v  v=s/r m=s/m (inversion)
 * Note that we don't care about the element types when determining
 * conformation; we'll delegate to the elements later and let them
 * deal with their own problems recursively.
 *
 * We only allow non-conforming arguments when (a) one argument is
 * a scalar and the other composite, or (b) one composite argument
 * has scalar elements and the other has composite elements. In case
 * (b), we treat the "less composite" argument as though it were a
 * scalar, and apply it to each of the elements of the "more composite"
 * argument. The cases are (e is the less composite 'element'):
 *      m= m op e   m=e op  m
 *     sy=sy op e  sy=e op sy
 *      v= v op e   v=e op  v
 *      r= r op e   r=e op  r
 *     
 * We assume that results types for scalars are dealt with 
 * elsewhere, and that we can access ResultType for elements
 * using CNT<LHSType>::Result<RHSType>::op where 'op' is Add,
 * Sub, Mul, Dvd. Note that elements of a CNT can be arbitrarily
 * complicated CNTs themselves.
 */

// The template arguments are repeated for LHS and RHS with these meanings:
//   NRow, NCol: Dimensions of the argument. NRow=1 for rows, NCol=1 for vecs,
//               and NRow=NCol=1 for scalars.
//     ArgDepth: Complexity of the argument as described above. 0=scalar,
//               1=composite with scalar elements, 2=composite of composites
//       L or R: The complete type of the LHS or RHS argument,
//               e.g. Row<5,Vec3,7>
//       CS, RS: Column spacing and row spacing (1 each for scalars). These
//               parameters are here because they are a "don't care" item;
//               the return types are always packed.

template <int LNRow, int LNCol, int LArgDepth, class L, int LCS, int LRS,
          int RNRow, int RNCol, int RArgDepth, class R, int RCS, int RRS>
class AddCNTs {
    // MUST BE SPECIALIZED
public:
    typedef void Type;
    static void perform(const L&,const R&) {assert(false);}
};
template <int LNRow, int LNCol, int LArgDepth, class L, int LCS, int LRS,
          int RNRow, int RNCol, int RArgDepth, class R, int RCS, int RRS>
class SubCNTs {
    // MUST BE SPECIALIZED
public:
    typedef void Type;
    static void perform(const L&,const R&) {assert(false);}
};
template <int LNRow, int LNCol, int LArgDepth, class L, int LCS, int LRS,
          int RNRow, int RNCol, int RArgDepth, class R, int RCS, int RRS>
class MulCNTs {
    // MUST BE SPECIALIZED
public:
    typedef void Type;
    static void perform(const L&,const R&) {assert(false);}
};
template <int LNRow, int LNCol, int LArgDepth, class L, int LCS, int LRS,
          int RNRow, int RNCol, int RArgDepth, class R, int RCS, int RRS>
class DvdCNTs {
    // MUST BE SPECIALIZED
public:
    typedef void Type;
    static void perform(const L&,const R&) {assert(false);}
};

///////////////////////////////////////////////////////////////////////////
// RULE 1: if an argument is a scalar, perform matrix/element operation, //
//         except for scalar on left of divide which is really just an   //
//         inversion of the RHS, followed by a scalar multiply.          //
///////////////////////////////////////////////////////////////////////////

// CNT+scalar, scalar+CNT
template <int LNRow, int LNCol, class L, int LCS, int LRS, class R>
class AddCNTs<LNRow,LNCol,SCALAR_COMPOSITE_DEPTH,L,LCS,LRS,1,1,SCALAR_DEPTH,R,1,1> {
    typedef typename CNT<L>::TElement LE;
    typedef typename CNT<LE>::template Result<R>::Add EAdd;
public:
    typedef typename CNT<L>::template Substitute<EAdd>::Type Type;
    static Type perform(const L& l, const R& r) {return l.scalarAdd(r);}
};
template <int LNRow, int LNCol, class L, int LCS, int LRS, class R>
class AddCNTs<LNRow,LNCol,COMPOSITE_COMPOSITE_DEPTH,L,LCS,LRS,1,1,SCALAR_DEPTH,R,1,1> {
    typedef typename CNT<L>::TElement LE;
    typedef typename CNT<LE>::template Result<R>::Add EAdd;
public:
    typedef typename CNT<L>::template Substitute<EAdd>::Type Type;
    static Type perform(const L& l, const R& r) {return l.scalarAdd(r);}
};
template <int LNRow, int LNCol, class L, int LCS, int LRS, class R>
class AddCNTs<LNRow,LNCol,COMPOSITE_3_DEPTH,L,LCS,LRS,1,1,SCALAR_DEPTH,R,1,1> {
    typedef typename CNT<L>::TElement LE;
    typedef typename CNT<LE>::template Result<R>::Add EAdd;
public:
    typedef typename CNT<L>::template Substitute<EAdd>::Type Type;
    static Type perform(const L& l, const R& r) {return l.scalarAdd(r);}
};

// scalar+CNT
template <class L, int RNRow, int RNCol, class R, int RCS, int RRS>
class AddCNTs<1,1,SCALAR_DEPTH,L,1,1, RNRow,RNCol,SCALAR_COMPOSITE_DEPTH,R,RCS,RRS> {
    typedef typename CNT<R>::TElement RE;
    typedef typename CNT<RE>::template Result<L>::Add EAdd;
public:
    typedef typename CNT<R>::template Substitute<EAdd>::Type Type;
    static Type perform(const L& l, const R& r) {return r.scalarAdd(l);}
};
template <class L, int RNRow, int RNCol, class R, int RCS, int RRS>
class AddCNTs<1,1,SCALAR_DEPTH,L,1,1, RNRow,RNCol,COMPOSITE_COMPOSITE_DEPTH,R,RCS,RRS> {
    typedef typename CNT<R>::TElement RE;
    typedef typename CNT<RE>::template Result<L>::Add EAdd;
public:
    typedef typename CNT<R>::template Substitute<EAdd>::Type Type;
    static Type perform(const L& l, const R& r) {return r.scalarAdd(l);}
};
template <class L, int RNRow, int RNCol, class R, int RCS, int RRS>
class AddCNTs<1,1,SCALAR_DEPTH,L,1,1, RNRow,RNCol,COMPOSITE_3_DEPTH,R,RCS,RRS> {
    typedef typename CNT<R>::TElement RE;
    typedef typename CNT<RE>::template Result<L>::Add EAdd;
public:
    typedef typename CNT<R>::template Substitute<EAdd>::Type Type;
    static Type perform(const L& l, const R& r) {return r.scalarAdd(l);}
};

// CNT-scalar
template <int LNRow, int LNCol, class L, int LCS, int LRS, class R>
class SubCNTs<LNRow,LNCol,SCALAR_COMPOSITE_DEPTH,L,LCS,LRS,1,1,SCALAR_DEPTH,R,1,1> {
    typedef typename CNT<L>::TElement LE;
    typedef typename CNT<LE>::template Result<R>::Sub ESub;
public:
    typedef typename CNT<L>::template Substitute<ESub>::Type Type;
    static Type perform(const L& l, const R& r) {return l.scalarSubtract(r);}
};
template <int LNRow, int LNCol, class L, int LCS, int LRS, class R>
class SubCNTs<LNRow,LNCol,COMPOSITE_COMPOSITE_DEPTH,L,LCS,LRS,1,1,SCALAR_DEPTH,R,1,1> {
    typedef typename CNT<L>::TElement LE;
    typedef typename CNT<LE>::template Result<R>::Sub ESub;
public:
    typedef typename CNT<L>::template Substitute<ESub>::Type Type;
    static Type perform(const L& l, const R& r) {return l.scalarSubtract(r);}
};
template <int LNRow, int LNCol, class L, int LCS, int LRS, class R>
class SubCNTs<LNRow,LNCol,COMPOSITE_3_DEPTH,L,LCS,LRS,1,1,SCALAR_DEPTH,R,1,1> {
    typedef typename CNT<L>::TElement LE;
    typedef typename CNT<LE>::template Result<R>::Sub ESub;
public:
    typedef typename CNT<L>::template Substitute<ESub>::Type Type;
    static Type perform(const L& l, const R& r) {return l.scalarSubtract(r);}
};
// scalar-CNT
template <class L, int RNRow, int RNCol, class R, int RCS, int RRS>
class SubCNTs<1,1,SCALAR_DEPTH,L,1,1, RNRow,RNCol,SCALAR_COMPOSITE_DEPTH,R,RCS,RRS> {
    typedef typename CNT<R>::TElement RE;
    typedef typename CNT<L>::template Result<RE>::Sub ESub;
public:
    typedef typename CNT<R>::template Substitute<ESub>::Type Type;
    static Type perform(const L& l, const R& r) {return r.scalarSubFromLeft(l);}
};
template <class L, int RNRow, int RNCol, class R, int RCS, int RRS>
class SubCNTs<1,1,SCALAR_DEPTH,L,1,1, RNRow,RNCol,COMPOSITE_COMPOSITE_DEPTH,R,RCS,RRS> {
    typedef typename CNT<R>::TElement RE;
    typedef typename CNT<L>::template Result<RE>::Sub ESub;
public:
    typedef typename CNT<R>::template Substitute<ESub>::Type Type;
    static Type perform(const L& l, const R& r) {return r.scalarSubFromLeft(l);}
};
template <class L, int RNRow, int RNCol, class R, int RCS, int RRS>
class SubCNTs<1,1,SCALAR_DEPTH,L,1,1, RNRow,RNCol,COMPOSITE_3_DEPTH,R,RCS,RRS> {
    typedef typename CNT<R>::TElement RE;
    typedef typename CNT<L>::template Result<RE>::Sub ESub;
public:
    typedef typename CNT<R>::template Substitute<ESub>::Type Type;
    static Type perform(const L& l, const R& r) {return r.scalarSubFromLeft(l);}
};

// CNT*scalar
template <int LNRow, int LNCol, class L, int LCS, int LRS, class R>
class MulCNTs<LNRow,LNCol,SCALAR_COMPOSITE_DEPTH,L,LCS,LRS,1,1,SCALAR_DEPTH,R,1,1> {
    typedef typename CNT<L>::TElement LE;
    typedef typename CNT<LE>::template Result<R>::Mul EMul;
public:
    typedef typename CNT<L>::template Substitute<EMul>::Type Type;
    static Type perform(const L& l, const R& r) {return l.scalarMultiply(r);}
};
template <int LNRow, int LNCol, class L, int LCS, int LRS, class R>
class MulCNTs<LNRow,LNCol,COMPOSITE_COMPOSITE_DEPTH,L,LCS,LRS,1,1,SCALAR_DEPTH,R,1,1> {
    typedef typename CNT<L>::TElement LE;
    typedef typename CNT<LE>::template Result<R>::Mul EMul;
public:
    typedef typename CNT<L>::template Substitute<EMul>::Type Type;
    static Type perform(const L& l, const R& r) {return l.scalarMultiply(r);}
};
template <int LNRow, int LNCol, class L, int LCS, int LRS, class R>
class MulCNTs<LNRow,LNCol,COMPOSITE_3_DEPTH,L,LCS,LRS,1,1,SCALAR_DEPTH,R,1,1> {
    typedef typename CNT<L>::TElement LE;
    typedef typename CNT<LE>::template Result<R>::Mul EMul;
public:
    typedef typename CNT<L>::template Substitute<EMul>::Type Type;
    static Type perform(const L& l, const R& r) {return l.scalarMultiply(r);}
};
// scalar*CNT
template <class L, int RNRow, int RNCol, class R, int RCS, int RRS>
class MulCNTs<1,1,SCALAR_DEPTH,L,1,1, RNRow,RNCol,SCALAR_COMPOSITE_DEPTH,R,RCS,RRS> {
    typedef typename CNT<R>::TElement RE;
    typedef typename CNT<L>::template Result<RE>::Mul EMul;
public:
    typedef typename CNT<R>::template Substitute<EMul>::Type Type;
    static Type perform(const L& l, const R& r) {return r.scalarMultiplyFromLeft(l);}
};
template <class L, int RNRow, int RNCol, class R, int RCS, int RRS>
class MulCNTs<1,1,SCALAR_DEPTH,L,1,1, RNRow,RNCol,COMPOSITE_COMPOSITE_DEPTH,R,RCS,RRS> {
    typedef typename CNT<R>::TElement RE;
    typedef typename CNT<L>::template Result<RE>::Mul EMul;
public:
    typedef typename CNT<R>::template Substitute<EMul>::Type Type;
    static Type perform(const L& l, const R& r) {return r.scalarMultiplyFromLeft(l);}
};
template <class L, int RNRow, int RNCol, class R, int RCS, int RRS>
class MulCNTs<1,1,SCALAR_DEPTH,L,1,1, RNRow,RNCol,COMPOSITE_3_DEPTH,R,RCS,RRS> {
    typedef typename CNT<R>::TElement RE;
    typedef typename CNT<L>::template Result<RE>::Mul EMul;
public:
    typedef typename CNT<R>::template Substitute<EMul>::Type Type;
    static Type perform(const L& l, const R& r) {return r.scalarMultiplyFromLeft(l);}
};

// CNT/scalar
template <int LNRow, int LNCol, class L, int LCS, int LRS, class R>
class DvdCNTs<LNRow,LNCol,SCALAR_COMPOSITE_DEPTH,L,LCS,LRS,1,1,SCALAR_DEPTH,R,1,1> {
    typedef typename CNT<L>::TElement LE;
    typedef typename CNT<LE>::template Result<R>::Dvd EDvd;
public:
    typedef typename CNT<L>::template Substitute<EDvd>::Type Type;
    static Type perform(const L& l, const R& r) {return l.scalarDivide(r);}
};
template <int LNRow, int LNCol, class L, int LCS, int LRS, class R>
class DvdCNTs<LNRow,LNCol,COMPOSITE_COMPOSITE_DEPTH,L,LCS,LRS,1,1,SCALAR_DEPTH,R,1,1> {
    typedef typename CNT<L>::TElement LE;
    typedef typename CNT<LE>::template Result<R>::Dvd EDvd;
public:
    typedef typename CNT<L>::template Substitute<EDvd>::Type Type;
    static Type perform(const L& l, const R& r) {return l.scalarDivide(r);}
};
template <int LNRow, int LNCol, class L, int LCS, int LRS, class R>
class DvdCNTs<LNRow,LNCol,COMPOSITE_3_DEPTH,L,LCS,LRS,1,1,SCALAR_DEPTH,R,1,1> {
    typedef typename CNT<L>::TElement LE;
    typedef typename CNT<LE>::template Result<R>::Dvd EDvd;
public:
    typedef typename CNT<L>::template Substitute<EDvd>::Type Type;
    static Type perform(const L& l, const R& r) {return l.scalarDivide(r);}
};

// scalar/CNT is really scalar * inverse(CNT) so we're combining the
// inverse with a scalar multiply as handled in MulCNTs above.
template <class L, int RNRow, int RNCol, class R, int RCS, int RRS>
class DvdCNTs<1,1,SCALAR_DEPTH,L,1,1, RNRow,RNCol,SCALAR_COMPOSITE_DEPTH,R,RCS,RRS> {
    typedef typename CNT<R>::TInverse RInv;
public:
    typedef typename CNT<RInv>::template Result<L>::Mul Type;
    static Type perform(const L& l, const R& r) {return l*r.invert();}
};
template <class L, int RNRow, int RNCol, class R, int RCS, int RRS>
class DvdCNTs<1,1,SCALAR_DEPTH,L,1,1, RNRow,RNCol,COMPOSITE_COMPOSITE_DEPTH,R,RCS,RRS> {
    typedef typename CNT<R>::TInverse RInv;
public:
    typedef typename CNT<RInv>::template Result<L>::Mul Type;
    static Type perform(const L& l, const R& r) {return l*r.invert();}
};
template <class L, int RNRow, int RNCol, class R, int RCS, int RRS>
class DvdCNTs<1,1,SCALAR_DEPTH,L,1,1, RNRow,RNCol,COMPOSITE_3_DEPTH,R,RCS,RRS> {
    typedef typename CNT<R>::TInverse RInv;
public:
    typedef typename CNT<RInv>::template Result<L>::Mul Type;
    static Type perform(const L& l, const R& r) {return l*r.invert();}
};

    // CONFORMING Add and Subtract

// conforming: vec=vec+vec
template <int NRow, int LArgDepth, class LE, int LCS, int LRS,
                    int RArgDepth, class RE, int RCS, int RRS>
class AddCNTs<NRow,1,LArgDepth,Vec<NRow,LE,LRS>,LCS,LRS,
              NRow,1,RArgDepth,Vec<NRow,RE,RRS>,RCS,RRS>
{
    typedef Vec<NRow,LE,LRS> L;
    typedef Vec<NRow,RE,RRS> R;
    typedef typename CNT<LE>::template Result<RE>::Add EAdd;
public:
    typedef Vec<NRow,EAdd> Type;
    static Type perform(const L& l, const R& r) {return l.conformingAdd(r);}
};
// conforming: vec=vec-vec
template <int NRow, int LArgDepth, class LE, int LCS, int LRS,
                    int RArgDepth, class RE, int RCS, int RRS>
class SubCNTs<NRow,1,LArgDepth,Vec<NRow,LE,LRS>,LCS,LRS,
              NRow,1,RArgDepth,Vec<NRow,RE,RRS>,RCS,RRS>
{
    typedef Vec<NRow,LE,LRS> L;
    typedef Vec<NRow,RE,RRS> R;
    typedef typename CNT<LE>::template Result<RE>::Sub ESub;
public:
    typedef Vec<NRow,ESub> Type;
    static Type perform(const L& l, const R& r) {return l.conformingSubtract(r);}
};
// conforming: row=row+row
template <int NCol, int LArgDepth, class LE, int LCS, int LRS,
                    int RArgDepth, class RE, int RCS, int RRS>
class AddCNTs<1,NCol,LArgDepth,Row<NCol,LE,LCS>,LCS,LRS,
              1,NCol,RArgDepth,Row<NCol,RE,RCS>,RCS,RRS>
{
    typedef Row<NCol,LE,LCS> L;
    typedef Row<NCol,RE,RCS> R;
    typedef typename CNT<LE>::template Result<RE>::Add EAdd;
public:
    typedef Row<NCol,EAdd> Type;
    static Type perform(const L& l, const R& r) {return l.conformingAdd(r);}
};
// conforming: row=row-row
template <int NCol, int LArgDepth, class LE, int LCS, int LRS,
                    int RArgDepth, class RE, int RCS, int RRS>
class SubCNTs<1,NCol,LArgDepth,Row<NCol,LE,LCS>,LCS,LRS,
              1,NCol,RArgDepth,Row<NCol,RE,RCS>,RCS,RRS>
{
    typedef Row<NCol,LE,LCS> L;
    typedef Row<NCol,RE,RCS> R;
    typedef typename CNT<LE>::template Result<RE>::Sub ESub;
public:
    typedef Row<NCol,ESub> Type;
    static Type perform(const L& l, const R& r) {return l.conformingSubtract(r);}
};
// conforming: sym=sym+sym
template <int Dim,  int LArgDepth, class LE, int LCS, int LRS,
                    int RArgDepth, class RE, int RCS, int RRS>
class AddCNTs<Dim,Dim,LArgDepth,SymMat<Dim,LE,LRS>,LCS,LRS,
              Dim,Dim,RArgDepth,SymMat<Dim,RE,RRS>,RCS,RRS>
{
    typedef SymMat<Dim,LE,LRS> L;
    typedef SymMat<Dim,RE,RRS> R;
    typedef typename CNT<LE>::template Result<RE>::Add EAdd;
public:
    typedef SymMat<Dim,EAdd> Type;
    static Type perform(const L& l, const R& r) {return l.conformingAdd(r);}
};
// conforming: sym=sym-sym
template <int Dim,  int LArgDepth, class LE, int LCS, int LRS,
                    int RArgDepth, class RE, int RCS, int RRS>
class SubCNTs<Dim,Dim,LArgDepth,SymMat<Dim,LE,LRS>,LCS,LRS,
              Dim,Dim,RArgDepth,SymMat<Dim,RE,RRS>,RCS,RRS>
{
    typedef SymMat<Dim,LE,LRS> L;
    typedef SymMat<Dim,RE,RRS> R;
    typedef typename CNT<LE>::template Result<RE>::Sub ESub;
public:
    typedef SymMat<Dim,ESub> Type;
    static Type perform(const L& l, const R& r) {return l.conformingSubtract(r);}
};
// conforming: mat=mat+mat, mat=mat+sym, mat=sym+mat
template <int NRow, int NCol, int LArgDepth, class LE, int LCS, int LRS,
                              int RArgDepth, class RE, int RCS, int RRS>
class AddCNTs<NRow,NCol,LArgDepth,Mat<NRow,NCol,LE,LCS,LRS>,LCS,LRS,
              NRow,NCol,RArgDepth,Mat<NRow,NCol,RE,RCS,RRS>,RCS,RRS>
{
    typedef Mat<NRow,NCol,LE,LCS,LRS> L;
    typedef Mat<NRow,NCol,RE,RCS,RRS> R;
    typedef typename CNT<LE>::template Result<RE>::Add EAdd;
public:
    typedef Mat<NRow,NCol,EAdd> Type;
    static Type perform(const L& l, const R& r) {return l.conformingAdd(r);}
};
template <int Dim, int LArgDepth, class LE, int LCS, int LRS,
                   int RArgDepth, class RE, int RCS, int RRS>
class AddCNTs<Dim,Dim,LArgDepth,Mat<Dim,Dim,LE,LCS,LRS>,LCS,LRS,
              Dim,Dim,RArgDepth,SymMat<Dim,RE,RRS>,     RCS,RRS>
{
    typedef Mat<Dim,Dim,LE,LCS,LRS> L;
    typedef SymMat<Dim,RE,RRS>      R;
    typedef typename CNT<LE>::template Result<RE>::Add EAdd;
public:
    typedef Mat<Dim,Dim,EAdd> Type;
    // Let the SymMat deal with this
    static Type perform(const L& l, const R& r) {return r.conformingAdd(l);}
};
template <int Dim, int LArgDepth, class LE, int LCS, int LRS,
                   int RArgDepth, class RE, int RCS, int RRS>
class AddCNTs<Dim,Dim,LArgDepth,SymMat<Dim,LE,LRS>,     LCS,LRS,
              Dim,Dim,RArgDepth,Mat<Dim,Dim,RE,RCS,RRS>,RCS,RRS>
{
    typedef SymMat<Dim,LE,LRS>      L;
    typedef Mat<Dim,Dim,RE,RCS,RRS> R;
    typedef typename CNT<LE>::template Result<RE>::Add EAdd;
public:
    typedef Mat<Dim,Dim,EAdd> Type;
    // Let the SymMat deal with this
    static Type perform(const L& l, const R& r) {return l.conformingAdd(r);}
};

// conforming mat=mat-mat, mat=mat-sym, mat=sym-mat
template <int NRow, int NCol, int LArgDepth, class LE, int LCS, int LRS,
                              int RArgDepth, class RE, int RCS, int RRS>
class SubCNTs<NRow,NCol,LArgDepth,Mat<NRow,NCol,LE,LCS,LRS>,LCS,LRS,
              NRow,NCol,RArgDepth,Mat<NRow,NCol,RE,RCS,RRS>,RCS,RRS>
{
    typedef Mat<NRow,NCol,LE,LCS,LRS> L;
    typedef Mat<NRow,NCol,RE,RCS,RRS> R;
    typedef typename CNT<LE>::template Result<RE>::Sub ESub;
public:
    typedef Mat<NRow,NCol,ESub> Type;
    static Type perform(const L& l, const R& r) {return l.conformingSubtract(r);}
};
template <int Dim, int LArgDepth, class LE, int LCS, int LRS,
                   int RArgDepth, class RE, int RCS, int RRS>
class SubCNTs<Dim,Dim,LArgDepth,Mat<Dim,Dim,LE,LCS,LRS>,LCS,LRS,
              Dim,Dim,RArgDepth,SymMat<Dim,RE,RRS>,     RCS,RRS>
{
    typedef Mat<Dim,Dim,LE,LCS,LRS> L;
    typedef SymMat<Dim,RE,RRS>      R;
    typedef typename CNT<LE>::template Result<RE>::Sub ESub;
public:
    typedef Mat<Dim,Dim,ESub> Type;
    // Let the SymMat deal with this
    static Type perform(const L& l, const R& r) {return r.conformingSubtractFromLeft(l);}
};
template <int Dim, int LArgDepth, class LE, int LCS, int LRS,
                   int RArgDepth, class RE, int RCS, int RRS>
class SubCNTs<Dim,Dim,LArgDepth,SymMat<Dim,LE,LRS>,     LCS,LRS,
              Dim,Dim,RArgDepth,Mat<Dim,Dim,RE,RCS,RRS>,RCS,RRS>
{
    typedef SymMat<Dim,LE,LRS>      L;
    typedef Mat<Dim,Dim,RE,RCS,RRS> R;
    typedef typename CNT<LE>::template Result<RE>::Sub ESub;
public:
    typedef Mat<Dim,Dim,ESub> Type;
    // Let the SymMat deal with this
    static Type perform(const L& l, const R& r) {return l.conformingSubtract(r);}
};

    // NON-CONFORMING ADD/SUBTRACT NOT ALLOWED (except for scalars)

    // CONFORMING MULTIPLY

// conforming: element=row*vec (dot product)
template <int Dim, int LArgDepth, class LE, int LCS, int LRS,
                   int RArgDepth, class RE, int RCS, int RRS>
class MulCNTs<1,Dim,LArgDepth,Row<Dim,LE,LCS>,LCS,LRS,
              Dim,1,RArgDepth,Vec<Dim,RE,RRS>,RCS,RRS>
{
    typedef Row<Dim,LE,LCS> L;
    typedef Vec<Dim,RE,RRS> R;
    typedef typename CNT<LE>::template Result<RE>::Mul EMul;
public:
    typedef EMul Type;
    static Type perform(const L& l, const R& r) {return l.conformingMultiply(r);}
};

// conforming: mat=vec*row (outer product)
template <int Dim, int LArgDepth, class LE, int LCS, int LRS,
                   int RArgDepth, class RE, int RCS, int RRS>
class MulCNTs<Dim,1,LArgDepth,Vec<Dim,LE,LRS>,LCS,LRS,
              1,Dim,RArgDepth,Row<Dim,RE,RCS>,RCS,RRS>
{
    typedef Vec<Dim,LE,LRS> L;
    typedef Row<Dim,RE,RCS> R;
    typedef typename CNT<LE>::template Result<RE>::Mul EMul;
public:
    typedef Mat<Dim,Dim,EMul> Type;
    static Type perform(const L& l, const R& r) {return l.conformingMultiply(r);}
};

// conforming: row=row*mat
template <int LNCol, int LArgDepth, class LE, int LCS, int LRS,
          int RNCol, int RArgDepth, class RE, int RCS, int RRS>
class MulCNTs<1,    LNCol,LArgDepth,Row<LNCol,LE,LCS>,          LCS,LRS,
              LNCol,RNCol,RArgDepth,Mat<LNCol,RNCol,RE,RCS,RRS>,RCS,RRS>
{
    typedef Row<LNCol,LE,LCS>           L;
    typedef Mat<LNCol,RNCol,RE,RCS,RRS> R;
    typedef typename CNT<LE>::template Result<RE>::Mul EMul;
public:
    typedef Row<RNCol,EMul> Type;
    // Give the job to the matrix rather than the row.
    static Type perform(const L& l, const R& r) {return r.conformingMultiplyFromLeft(l);}
};
// conforming: row=row*symmat
template <int LNCol, int LArgDepth, class LE, int LCS, int LRS,
                     int RArgDepth, class RE, int RCS, int RRS>
class MulCNTs<1,    LNCol,LArgDepth,Row<LNCol,LE,LCS>,   LCS,LRS,
              LNCol,LNCol,RArgDepth,SymMat<LNCol,RE,RRS>,RCS,RRS>
{
    typedef Row<LNCol,LE,LCS>    L;
    typedef SymMat<LNCol,RE,RRS> R;
    typedef typename CNT<LE>::template Result<RE>::Mul EMul;
public:
    typedef Row<LNCol,EMul> Type;
    // Give the job to the matrix rather than the row.
    static Type perform(const L& l, const R& r) {return r.conformingMultiplyFromLeft(l);}
};

// conforming: vec=mat*vec
template <int LNRow, int LNCol, int LArgDepth, class LE, int LCS, int LRS,
                                int RArgDepth, class RE, int RCS, int RRS>
class MulCNTs<LNRow,LNCol,LArgDepth,Mat<LNRow,LNCol,LE,LCS,LRS>,LCS,LRS,
              LNCol,1,    RArgDepth,Vec<LNCol,RE,RRS>,          RCS,RRS>
{
    typedef Mat<LNRow,LNCol,LE,LCS,LRS> L;
    typedef Vec<LNCol,RE,RRS>           R;
    typedef typename CNT<LE>::template Result<RE>::Mul EMul;
public:
    typedef Vec<LNRow,EMul> Type;
    // We want the matrix to deal with this rather than the vec.
    static Type perform(const L& l, const R& r) {return l.conformingMultiply(r);}
};

// conforming: vec=sym*vec
template <int Dim, int LArgDepth, class LE, int LCS, int LRS,
                   int RArgDepth, class RE, int RCS, int RRS>
class MulCNTs<Dim,Dim,LArgDepth,SymMat<Dim,LE,LRS>,LCS,LRS,
              Dim,1,  RArgDepth,Vec<Dim,RE,RRS>,   RCS,RRS>
{
    typedef SymMat<Dim,LE,LRS> L;
    typedef Vec<Dim,RE,RRS>    R;
    typedef typename CNT<LE>::template Result<RE>::Mul EMul;
public:
    typedef Vec<Dim,EMul> Type;
    // We want the matrix to deal with this rather than the vec.
    static Type perform(const L& l, const R& r) {return l.conformingMultiply(r);}
};

// conforming: mat=mat*mat
template <int LNRow, int LNCol, int LArgDepth, class LE, int LCS, int LRS,
                     int RNCol, int RArgDepth, class RE, int RCS, int RRS>
class MulCNTs<LNRow,LNCol,LArgDepth,Mat<LNRow,LNCol,LE,LCS,LRS>,LCS,LRS,
              LNCol,RNCol,RArgDepth,Mat<LNCol,RNCol,RE,RCS,RRS>,RCS,RRS>
{
    typedef Mat<LNRow,LNCol,LE,LCS,LRS> L;
    typedef Mat<LNCol,RNCol,RE,RCS,RRS> R;
    typedef typename CNT<LE>::template Result<RE>::Mul EMul;
public:
    typedef Mat<LNRow,RNCol,EMul> Type;
    static Type perform(const L& l, const R& r) {return l.conformingMultiply(r);}
};
// conforming: mat=mat*sym
template <int LNRow, int LNCol, int LArgDepth, class LE, int LCS, int LRS,
                                int RArgDepth, class RE, int RCS, int RRS>
class MulCNTs<LNRow,LNCol,LArgDepth,Mat<LNRow,LNCol,LE,LCS,LRS>,LCS,LRS,
              LNCol,LNCol,RArgDepth,SymMat<LNCol,RE,RRS>,       RCS,RRS>
{
    typedef Mat<LNRow,LNCol,LE,LCS,LRS> L;
    typedef SymMat<LNCol,RE,RRS>        R;
    typedef typename CNT<LE>::template Result<RE>::Mul EMul;
public:
    typedef Mat<LNRow,LNCol,EMul> Type;
    // Let SymMat deal with this.
    static Type perform(const L& l, const R& r) {return r.conformingMultiplyFromLeft(l);}
};
// conforming: mat=sym*mat
template <int LDim,  int LArgDepth, class LE, int LCS, int LRS,
          int RNCol, int RArgDepth, class RE, int RCS, int RRS>
class MulCNTs<LDim,LDim, LArgDepth,SymMat<LDim,LE,LRS>,       LCS,LRS,
              LDim,RNCol,RArgDepth,Mat<LDim,RNCol,RE,RCS,RRS>,RCS,RRS>
{
    typedef SymMat<LDim,LE,LRS>        L;
    typedef Mat<LDim,RNCol,RE,RCS,RRS> R;
    typedef typename CNT<LE>::template Result<RE>::Mul EMul;
public:
    typedef Mat<LDim,RNCol,EMul> Type;
    // Let SymMat deal with this.
    static Type perform(const L& l, const R& r) {return l.conformingMultiply(r);}
};
// conforming: mat=sym*sym
template <int Dim, int LArgDepth, class LE, int LCS, int LRS,
                   int RArgDepth, class RE, int RCS, int RRS>
class MulCNTs<Dim,Dim,LArgDepth,SymMat<Dim,LE,LRS>,LCS,LRS,
              Dim,Dim,RArgDepth,SymMat<Dim,RE,RRS>,RCS,RRS>
{
    typedef SymMat<Dim,LE,LRS> L;
    typedef SymMat<Dim,RE,RRS> R;
    typedef typename CNT<LE>::template Result<RE>::Mul EMul;
public:
    typedef Mat<Dim,Dim,EMul> Type;
    static Type perform(const L& l, const R& r) {return l.conformingMultiply(r);}
};

    // TODO: CONFORMING DIVIDE

    // NON-CONFORMING MULTIPLY
template <int LNRow, int LNCol, int LArgDepth, class L, int LCS, int LRS,
          int RNRow, int RNCol, int RArgDepth, class R, int RCS, int RRS>
class MulCNTsNonConforming {
    // MUST BE SPECIALIZED
public:
    typedef void Type;
    static void perform(const L&,const R&) {assert(false);}
};

// nonconforming: left = left(2)*right(1) (like scalar multiply)
template <int LNRow, int LNCol, class L, int LCS, int LRS,
          int RNRow, int RNCol, class R, int RCS, int RRS>
class MulCNTsNonConforming<LNRow,LNCol,COMPOSITE_COMPOSITE_DEPTH,L,LCS,LRS,
                           RNRow,RNCol,SCALAR_COMPOSITE_DEPTH,   R,RCS,RRS>
{
    typedef typename CNT<L>::TElement LE;
    typedef typename CNT<R>::TElement RE;
    typedef typename CNT<LE>::template Result<R>::Mul EMul;
public:
    typedef typename CNT<L>::template Substitute<EMul>::Type Type;
    static Type perform(const L& l, const R& r) {return l.scalarMultiply(r);}
};
// nonconforming: left = left(3)*right(1) (like scalar multiply)
template <int LNRow, int LNCol, class L, int LCS, int LRS,
          int RNRow, int RNCol, class R, int RCS, int RRS>
class MulCNTsNonConforming<LNRow,LNCol,COMPOSITE_3_DEPTH,        L,LCS,LRS,
                           RNRow,RNCol,SCALAR_COMPOSITE_DEPTH,   R,RCS,RRS>
{
    typedef typename CNT<L>::TElement LE;
    typedef typename CNT<R>::TElement RE;
    typedef typename CNT<LE>::template Result<R>::Mul EMul;
public:
    typedef typename CNT<L>::template Substitute<EMul>::Type Type;
    static Type perform(const L& l, const R& r) {return l.scalarMultiply(r);}
};
// nonconforming: left = left(3)*right(2) (like scalar multiply)
template <int LNRow, int LNCol, class L, int LCS, int LRS,
          int RNRow, int RNCol, class R, int RCS, int RRS>
class MulCNTsNonConforming<LNRow,LNCol,COMPOSITE_3_DEPTH,         L,LCS,LRS,
                           RNRow,RNCol,COMPOSITE_COMPOSITE_DEPTH, R,RCS,RRS>
{
    typedef typename CNT<L>::TElement LE;
    typedef typename CNT<R>::TElement RE;
    typedef typename CNT<LE>::template Result<R>::Mul EMul;
public:
    typedef typename CNT<L>::template Substitute<EMul>::Type Type;
    static Type perform(const L& l, const R& r) {return l.scalarMultiply(r);}
};

// nonconforming: right = left(1)*right(2) (like scalar multiply)
template <int LNRow, int LNCol, class L, int LCS, int LRS,
          int RNRow, int RNCol, class R, int RCS, int RRS>
class MulCNTsNonConforming<LNRow,LNCol,SCALAR_COMPOSITE_DEPTH,   L,LCS,LRS,
                           RNRow,RNCol,COMPOSITE_COMPOSITE_DEPTH,R,RCS,RRS>
{
    typedef typename CNT<L>::TElement LE;
    typedef typename CNT<R>::TElement RE;
    typedef typename CNT<L>::template Result<RE>::Mul EMul;
public:
    typedef typename CNT<R>::template Substitute<EMul>::Type Type;
    static Type perform(const L& l, const R& r) {return r.scalarMultiplyFromLeft(l);}
};
// nonconforming: right = left(1)*right(3) (like scalar multiply)
template <int LNRow, int LNCol, class L, int LCS, int LRS,
          int RNRow, int RNCol, class R, int RCS, int RRS>
class MulCNTsNonConforming<LNRow,LNCol,SCALAR_COMPOSITE_DEPTH,   L,LCS,LRS,
                           RNRow,RNCol,COMPOSITE_3_DEPTH,R,RCS,RRS>
{
    typedef typename CNT<L>::TElement LE;
    typedef typename CNT<R>::TElement RE;
    typedef typename CNT<L>::template Result<RE>::Mul EMul;
public:
    typedef typename CNT<R>::template Substitute<EMul>::Type Type;
    static Type perform(const L& l, const R& r) {return r.scalarMultiplyFromLeft(l);}
};
// nonconforming: right = left(2)*right(3) (like scalar multiply)
template <int LNRow, int LNCol, class L, int LCS, int LRS,
          int RNRow, int RNCol, class R, int RCS, int RRS>
class MulCNTsNonConforming<LNRow,LNCol,COMPOSITE_COMPOSITE_DEPTH,   L,LCS,LRS,
                           RNRow,RNCol,COMPOSITE_3_DEPTH,R,RCS,RRS>
{
    typedef typename CNT<L>::TElement LE;
    typedef typename CNT<R>::TElement RE;
    typedef typename CNT<L>::template Result<RE>::Mul EMul;
public:
    typedef typename CNT<R>::template Substitute<EMul>::Type Type;
    static Type perform(const L& l, const R& r) {return r.scalarMultiplyFromLeft(l);}
};

    // TODO: NON-CONFORMING DIVIDE


} //namespace SimTK


#endif // SimTK_SIMMATRIX_RESULT_TYPE_H_
