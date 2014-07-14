#ifndef SimTK_SIMMATRIX_BIGMATRIX_H_
#define SimTK_SIMMATRIX_BIGMATRIX_H_

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

/** @file
 * This file defines the client side of the SimTK::Matrix classes, which
 * hold medium to large, variable-sized matrices whose elements are packed 
 * SimTK "Composite Numerical Types" (CNTs). Unlike CNTs, the implemention here 
 * is opaque, and almost all properties are captured in the implementation at 
 * run time rather than in the type at compile time. 
 *
 * Every Matrix consists logically of three pieces: 
 *  - the matrix handle 
 *  - the matrix helper
 *  - and the matrix data. 
 *
 * They are organized like this:
 * <pre>
 *      ------------            ------------
 *     |  Handle<E> | -------> |            |
 *      ------------  <------- | Helper<S>  |
 *                             |            |
 *                             |            |          --------~ ~--
 *                             |            | ------> | Data<S> ... |
 *                              ------------           --------~ ~--
 * </pre>
 * The handle is the object actually appearing in SimTK API user programs.
 * It always consists of just a single pointer, pointing to a library-side
 * "helper" object whose implementation is opaque. The handle is templatized
 * by the user's element type, which may be any packed composite numerical
 * type, including scalar types like \c float and \c complex<double>, but also
 * composite types such as \c Vec3 or \c Mat<2,2,Mat<3,3>>. A Matrix handle
 * owns the helper to which it points and must destruct the helper when
 * the handle's destructor is called.
 *
 * The helper, on the other hand, is parameterized only by the underlying scalar 
 * type. There are exactly 12 SimTK scalar types, so all can be instantiated on 
 * the library side leaving the implementation opaque and thus flexible from
 * release to release without compromising binary compatibility. (The scalar 
 * types are: the four C++ standard types float and double, 
 * complex<float>, and complex<double>; the SimTK numbers conjugate<float> and 
 * conjugate<double>; and negator<> types templatized by any of the six
 * numeric types.) The helper contains several kinds of information:
 *  - the underlying scalar type S (as its template parameter)
 *  - the number of scalars in the handle's logical element type E
 *  - whether this is an owner matrix, or just a view
 *  - the handle "commitment"; defining the range of matrix characteristics
 *      to which that handle may refer
 *  - the actual characteristics of the matrix currently represented by
 *      the helper
 *  - a virtual function table full of methods which are aware of the
 *      logical structure of the Matrix and the physical structure of 
 *      the data to support operations such as element indexing
 *  - a pointer to the underlying data, which may be shared with other 
 *      helpers
 *
 * The data itself consists only of scalars
 * S of the same type as the helper's template argument, but different 
 * helpers can look at the same data differently. For examples, when the
 * elements are composite consisting of k scalars, the helper will provide a 
 * view of the data in which its scalars are interpreted in groups of k.
 * Many other reinterpretations of the data are possible and useful, such
 * as a real-valued helper viewing only the real or imaginary part of
 * complex data, or a helper which views the data order as though it were
 * transposed.
 *
 * At most \e one matrix helper owns the matrix data and is responsible
 * for deleting that data when no longer needed. That is called an "owner"
 * helper and its associated handle is an owner handle. Normally the owner 
 * is the handle (and helper) that allocated the data, and
 * in most cases an owner can resize the data at will. Many other handles
 * may reference the same data; those non-owner handles are called "views".
 * Every view may present a different picture of the underlying data. The
 * default view is "whole" meaning that all the elements of the data are 
 * visible, and appear in their normal order. A "transpose" view also shows 
 * all the elements but the matrix dimensions and indices are reversed. 
 * Other common views are "block" to select a sub-block of a matrix, and 
 * "diagonal" which shows only the diagonal of a matrix (as a vector). 
 *
 * NOTE: Destruction of an owner destructs the data it owns
 * *regardless* of the presence of other views into that data! I.e., these
 * are not reference counted. TODO: should we change that?
 * 
 * In some cases there may be no owner helper for a particular piece of 
 * matrix data. That occurs when pre-existing memory, such as a Fortran
 * array, is used to construct a Matrix. In that case all the helpers are
 * views and the data will persist after the destruction of the last
 * referencing helper.
 *                 
 * A Matrix that is the owner of its data will be resized whenever
 * necessary, unless you take active steps to prevent that. For example, if
 * you declare a Vector, the number of rows can resize but the number of
 * columns will be locked at 1. A RowVector does the reverse. You can also
 * explicitly lock the number of rows and/or columns of a matrix to prevent
 * unwanted resizes.
 *
 * Here are the classes and short descriptions:
 * <pre>
 *   MatrixHelper<S>  interface to the opaque implementation, templatized
 *                      by scalar type only
 *   MatrixBase<CNT>  fully templatized client, contains a MatrixHelper
 * </pre>
 *
 * The rest are dataless classes all of which can be interconverted just
 * by recasting. Every one of these classes has a default conversion to
 * type Matrix_<same element type>, so users can write functions that expect
 * a Matrix argument and pass it a Vector, RowVectorView, or whatever.
 *
 * <pre>
 *   VectorBase<CNT>    these are derived from MatrixBase and add no new data,
 *   RowVectorBase<CNT> but change some of the operators and other methods to
 *                        be appropriate for 1d data.
 *
 *   Matrix_<CNT>      2d owner class     (a MatrixBase<CNT>)
 *   Vector_<CNT>      column owner class (a VectorBase<CNT>)
 *   RowVector_<CNT>   row owner class    (a RowVectorBase<CNT>)
 * </pre>
 *
 * Views are exactly the same as the corresponding owner class, but with
 * shallow construction and assignment semantics.
 *
 * <pre>
 *   MatrixView_<CNT>, VectorView_<CNT>, RowVectorView_<CNT>
 * </pre>
 *
 */

#include "SimTKcommon/Scalar.h"
#include "SimTKcommon/SmallMatrix.h"

#include "SimTKcommon/internal/MatrixHelper.h"
#include "SimTKcommon/internal/MatrixCharacteristics.h"

#include <iostream>
#include <cassert>
#include <complex>
#include <cstddef>
#include <limits>

namespace SimTK {
    template <class ELT>    class MatrixBase;
    template <class ELT>    class VectorBase;
    template <class ELT>    class RowVectorBase;

    template <class ELT = Real> class MatrixView_;
    template <class ELT = Real> class Matrix_;

    template <class ELT = Real> class VectorView_;
    template <class ELT = Real> class Vector_;

    template <class ELT = Real> class RowVectorView_;
    template <class ELT = Real> class RowVector_;

    template <class ELT, class VECTOR_CLASS> class VectorIterator;
}

#include "SimTKcommon/internal/MatrixBase.h"
#include "SimTKcommon/internal/VectorBase.h"
#include "SimTKcommon/internal/RowVectorBase.h"

#include "SimTKcommon/internal/MatrixView_.h"
#include "SimTKcommon/internal/Matrix_.h"

#include "SimTKcommon/internal/VectorView_.h"
#include "SimTKcommon/internal/Vector_.h"

#include "SimTKcommon/internal/RowVectorView_.h"
#include "SimTKcommon/internal/RowVector_.h"

#include "SimTKcommon/internal/VectorIterator.h"


namespace SimTK {

//  ------------------------ MatrixBase definitions ----------------------------

template <class ELT> inline MatrixView_<ELT> 
MatrixBase<ELT>::block(int i, int j, int m, int n) const { 
    SimTK_INDEXCHECK(i,nrow()+1,"MatrixBase::block()");
    SimTK_INDEXCHECK(j,ncol()+1,"MatrixBase::block()");
    SimTK_SIZECHECK(i+m,nrow(),"MatrixBase::block()");
    SimTK_SIZECHECK(j+n,ncol(),"MatrixBase::block()");

    MatrixHelper<Scalar> h(MatrixCommitment(),helper,i,j,m,n);    
    return MatrixView_<ELT>(h.stealRep()); 
}
    
template <class ELT> inline MatrixView_<ELT>
MatrixBase<ELT>::updBlock(int i, int j, int m, int n) { 
    SimTK_INDEXCHECK(i,nrow()+1,"MatrixBase::updBlock()");
    SimTK_INDEXCHECK(j,ncol()+1,"MatrixBase::updBlock()");
    SimTK_SIZECHECK(i+m,nrow(),"MatrixBase::updBlock()");
    SimTK_SIZECHECK(j+n,ncol(),"MatrixBase::updBlock()");

    MatrixHelper<Scalar> h(MatrixCommitment(),helper,i,j,m,n);        
    return MatrixView_<ELT>(h.stealRep()); 
}

template <class E> inline MatrixView_<typename CNT<E>::THerm>
MatrixBase<E>::transpose() const { 
    MatrixHelper<typename CNT<Scalar>::THerm> 
        h(MatrixCommitment(),
          helper, typename MatrixHelper<typename CNT<Scalar>::THerm>::TransposeView());
    return MatrixView_<typename CNT<E>::THerm>(h.stealRep()); 
}
    
template <class E> inline MatrixView_<typename CNT<E>::THerm>
MatrixBase<E>::updTranspose() {     
    MatrixHelper<typename CNT<Scalar>::THerm> 
        h(MatrixCommitment(),
          helper, typename MatrixHelper<typename CNT<Scalar>::THerm>::TransposeView());
    return MatrixView_<typename CNT<E>::THerm>(h.stealRep()); 
}

template <class E> inline VectorView_<E>
MatrixBase<E>::diag() const { 
    MatrixHelper<Scalar> h(MatrixCommitment::Vector(),
                           helper, typename MatrixHelper<Scalar>::DiagonalView());
    return VectorView_<E>(h.stealRep()); 
}
    
template <class E> inline VectorView_<E>
MatrixBase<E>::updDiag() {     
    MatrixHelper<Scalar> h(MatrixCommitment::Vector(),
                           helper, typename MatrixHelper<Scalar>::DiagonalView());
    return VectorView_<E>(h.stealRep());
}

template <class ELT> inline VectorView_<ELT> 
MatrixBase<ELT>::col(int j) const { 
    SimTK_INDEXCHECK(j,ncol(),"MatrixBase::col()");

    MatrixHelper<Scalar> h(MatrixCommitment::Vector(),
                           helper,0,j,nrow(),1);    
    return VectorView_<ELT>(h.stealRep()); 
}
    
template <class ELT> inline VectorView_<ELT>
MatrixBase<ELT>::updCol(int j) {
    SimTK_INDEXCHECK(j,ncol(),"MatrixBase::updCol()");

    MatrixHelper<Scalar> h(MatrixCommitment::Vector(),
                           helper,0,j,nrow(),1);        
    return VectorView_<ELT>(h.stealRep()); 
}

template <class ELT> inline RowVectorView_<ELT> 
MatrixBase<ELT>::row(int i) const { 
    SimTK_INDEXCHECK(i,nrow(),"MatrixBase::row()");

    MatrixHelper<Scalar> h(MatrixCommitment::RowVector(),
                           helper,i,0,1,ncol());    
    return RowVectorView_<ELT>(h.stealRep()); 
}
    
template <class ELT> inline RowVectorView_<ELT>
MatrixBase<ELT>::updRow(int i) { 
    SimTK_INDEXCHECK(i,nrow(),"MatrixBase::updRow()");

    MatrixHelper<Scalar> h(MatrixCommitment::RowVector(),
                           helper,i,0,1,ncol());        
    return RowVectorView_<ELT>(h.stealRep()); 
}

// M = diag(v) * M; v must have nrow() elements.
// That is, M[i] *= v[i].
template <class ELT> template <class EE> inline MatrixBase<ELT>& 
MatrixBase<ELT>::rowScaleInPlace(const VectorBase<EE>& v) {
	assert(v.nrow() == nrow());
	for (int i=0; i < nrow(); ++i)
		(*this)[i] *= v[i];
	return *this;
}

template <class ELT> template <class EE> inline void
MatrixBase<ELT>::rowScale(const VectorBase<EE>& v, typename MatrixBase<ELT>::template EltResult<EE>::Mul& out) const {
	assert(v.nrow() == nrow());
	out.resize(nrow(), ncol());
    for (int j=0; j<ncol(); ++j)
	    for (int i=0; i<nrow(); ++i)
		   out(i,j) = (*this)(i,j) * v[i];
}

// M = M * diag(v); v must have ncol() elements
// That is, M(i) *= v[i]
template <class ELT> template <class EE>  inline MatrixBase<ELT>& 
MatrixBase<ELT>::colScaleInPlace(const VectorBase<EE>& v) {
    assert(v.nrow() == ncol());
	for (int j=0; j < ncol(); ++j)
		(*this)(j) *= v[j];
	return *this;
}

template <class ELT> template <class EE> inline void
MatrixBase<ELT>::colScale(const VectorBase<EE>& v, typename MatrixBase<ELT>::template EltResult<EE>::Mul& out) const {
	assert(v.nrow() == ncol());
	out.resize(nrow(), ncol());
    for (int j=0; j<ncol(); ++j)
	    for (int i=0; i<nrow(); ++i)
		   out(i,j) = (*this)(i,j) * v[j];
}


// M(i,j) *= r[i]*c[j]; r must have nrow() elements; c must have ncol() elements
template <class ELT> template <class ER, class EC> inline MatrixBase<ELT>& 
MatrixBase<ELT>::rowAndColScaleInPlace(const VectorBase<ER>& r, const VectorBase<EC>& c) {
	assert(r.nrow()==nrow() && c.nrow()==ncol());
    for (int j=0; j<ncol(); ++j)
	    for (int i=0; i<nrow(); ++i)
			(*this)(i,j) *= (r[i]*c[j]);
    return *this;
}

template <class ELT> template <class ER, class EC> inline void
MatrixBase<ELT>::rowAndColScale(
    const VectorBase<ER>& r, 
    const VectorBase<EC>& c,
    typename EltResult<typename VectorBase<ER>::template EltResult<EC>::Mul>::Mul& 
                          out) const
{
	assert(r.nrow()==nrow() && c.nrow()==ncol());
    out.resize(nrow(), ncol());
    for (int j=0; j<ncol(); ++j)
	    for (int i=0; i<nrow(); ++i)
			out(i,j) = (*this)(i,j) * (r[i]*c[j]);
}

// M(i,j) = s
template <class ELT> template <class S> inline MatrixBase<ELT>& 
MatrixBase<ELT>::elementwiseAssign(const S& s) {
    for (int j=0; j<ncol(); ++j)
	for (int i=0; i<nrow(); ++i)
	    (*this)(i,j) = s;
    return *this;
}

// Set M(i,j) = M(i,j)^-1.
template <class ELT> inline MatrixBase<ELT>& 
MatrixBase<ELT>::elementwiseInvertInPlace() {
    const int nr=nrow(), nc=ncol();
    for (int j=0; j<nc; ++j)
        for (int i=0; i<nr; ++i) {
            ELT& e = updElt(i,j);
            e = CNT<ELT>::invert(e);
        }
    return *this;
}

template <class ELT> inline void
MatrixBase<ELT>::elementwiseInvert(MatrixBase<typename CNT<E>::TInvert>& out) const {
    const int nr=nrow(), nc=ncol();
    out.resize(nr,nc);
    for (int j=0; j<nc; ++j)
        for (int i=0; i<nr; ++i)
            out(i,j) = CNT<ELT>::invert((*this)(i,j));
}

// M(i,j) += s
template <class ELT> template <class S> inline MatrixBase<ELT>& 
MatrixBase<ELT>::elementwiseAddScalarInPlace(const S& s) {
    for (int j=0; j<ncol(); ++j)
	    for (int i=0; i<nrow(); ++i)
			(*this)(i,j) += s;
    return *this;
}

template <class ELT> template <class S> inline void 
MatrixBase<ELT>::elementwiseAddScalar(
    const S& s,
    typename MatrixBase<ELT>::template EltResult<S>::Add& out) const
{
    const int nr=nrow(), nc=ncol();
    out.resize(nr,nc);
    for (int j=0; j<nc; ++j)
	    for (int i=0; i<nr; ++i)
			out(i,j) = (*this)(i,j) + s;
}

// M(i,j) -= s
template <class ELT> template <class S> inline MatrixBase<ELT>& 
MatrixBase<ELT>::elementwiseSubtractScalarInPlace(const S& s) {
    for (int j=0; j<ncol(); ++j)
	    for (int i=0; i<nrow(); ++i)
			(*this)(i,j) -= s;
    return *this;
}

template <class ELT> template <class S> inline void 
MatrixBase<ELT>::elementwiseSubtractScalar(
    const S& s,
    typename MatrixBase<ELT>::template EltResult<S>::Sub& out) const
{
    const int nr=nrow(), nc=ncol();
    out.resize(nr,nc);
    for (int j=0; j<nc; ++j)
	    for (int i=0; i<nr; ++i)
			out(i,j) = (*this)(i,j) - s;
}

// M(i,j) = s - M(i,j)
template <class ELT> template <class S> inline MatrixBase<ELT>& 
MatrixBase<ELT>::elementwiseSubtractFromScalarInPlace(const S& s) {
    const int nr=nrow(), nc=ncol();
    for (int j=0; j<nc; ++j)
        for (int i=0; i<nr; ++i) {
            ELT& e = updElt(i,j);
			e = s - e;
        }
    return *this;
}

template <class ELT> template <class S> inline void 
MatrixBase<ELT>::elementwiseSubtractFromScalar(
    const S& s,
    typename MatrixBase<S>::template EltResult<ELT>::Sub& out) const
{
    const int nr=nrow(), nc=ncol();
    out.resize(nr,nc);
    for (int j=0; j<nc; ++j)
	    for (int i=0; i<nr; ++i)
			out(i,j) = s - (*this)(i,j);
}

// M(i,j) *= R(i,j); R must have same dimensions as this
template <class ELT> template <class EE> inline MatrixBase<ELT>& 
MatrixBase<ELT>::elementwiseMultiplyInPlace(const MatrixBase<EE>& r) {
    const int nr=nrow(), nc=ncol();
	assert(r.nrow()==nr && r.ncol()==nc);
    for (int j=0; j<nc; ++j)
	    for (int i=0; i<nr; ++i)
			(*this)(i,j) *= r(i,j);
    return *this;
}

template <class ELT> template <class EE> inline void 
MatrixBase<ELT>::elementwiseMultiply(
    const MatrixBase<EE>& r, 
    typename MatrixBase<ELT>::template EltResult<EE>::Mul& out) const
{
    const int nr=nrow(), nc=ncol();
	assert(r.nrow()==nr && r.ncol()==nc);
	out.resize(nr,nc);
    for (int j=0; j<nc; ++j)
	    for (int i=0; i<nr; ++i)
			out(i,j) = (*this)(i,j) * r(i,j);
}

// M(i,j) = R(i,j) * M(i,j); R must have same dimensions as this
template <class ELT> template <class EE> inline MatrixBase<ELT>& 
MatrixBase<ELT>::elementwiseMultiplyFromLeftInPlace(const MatrixBase<EE>& r) {
    const int nr=nrow(), nc=ncol();
	assert(r.nrow()==nr && r.ncol()==nc);
    for (int j=0; j<nc; ++j)
        for (int i=0; i<nr; ++i) {
            ELT& e = updElt(i,j);
			e = r(i,j) * e;
        }
    return *this;
}

template <class ELT> template <class EE> inline void 
MatrixBase<ELT>::elementwiseMultiplyFromLeft(
    const MatrixBase<EE>& r, 
    typename MatrixBase<EE>::template EltResult<ELT>::Mul& out) const
{
    const int nr=nrow(), nc=ncol();
	assert(r.nrow()==nr && r.ncol()==nc);
	out.resize(nr,nc);
    for (int j=0; j<nc; ++j)
	    for (int i=0; i<nr; ++i)
			out(i,j) =  r(i,j) * (*this)(i,j);
}

// M(i,j) /= R(i,j); R must have same dimensions as this
template <class ELT> template <class EE> inline MatrixBase<ELT>& 
MatrixBase<ELT>::elementwiseDivideInPlace(const MatrixBase<EE>& r) {
    const int nr=nrow(), nc=ncol();
	assert(r.nrow()==nr && r.ncol()==nc);
    for (int j=0; j<nc; ++j)
	    for (int i=0; i<nr; ++i)
			(*this)(i,j) /= r(i,j);
    return *this;
}

template <class ELT> template <class EE> inline void 
MatrixBase<ELT>::elementwiseDivide(
    const MatrixBase<EE>& r,
    typename MatrixBase<ELT>::template EltResult<EE>::Dvd& out) const
{
    const int nr=nrow(), nc=ncol();
	assert(r.nrow()==nr && r.ncol()==nc);
	out.resize(nr,nc);
    for (int j=0; j<nc; ++j)
	    for (int i=0; i<nr; ++i)
			out(i,j) = (*this)(i,j) / r(i,j);
}
// M(i,j) = R(i,j) / M(i,j); R must have same dimensions as this
template <class ELT> template <class EE> inline MatrixBase<ELT>& 
MatrixBase<ELT>::elementwiseDivideFromLeftInPlace(const MatrixBase<EE>& r) {
    const int nr=nrow(), nc=ncol();
	assert(r.nrow()==nr && r.ncol()==nc);
    for (int j=0; j<nc; ++j)
        for (int i=0; i<nr; ++i) {
            ELT& e = updElt(i,j);
			e = r(i,j) / e;
        }
    return *this;
}

template <class ELT> template <class EE> inline void 
MatrixBase<ELT>::elementwiseDivideFromLeft(
    const MatrixBase<EE>& r,
    typename MatrixBase<EE>::template EltResult<ELT>::Dvd& out) const
{
    const int nr=nrow(), nc=ncol();
	assert(r.nrow()==nr && r.ncol()==nc);
	out.resize(nr,nc);
    for (int j=0; j<nc; ++j)
	    for (int i=0; i<nr; ++i)
			out(i,j) = r(i,j) / (*this)(i,j);
}

/*
template <class ELT> inline MatrixView_< typename CNT<ELT>::TReal > 
MatrixBase<ELT>::real() const { 
    if (!CNT<ELT>::IsComplex) { // known at compile time
        return MatrixView_< typename CNT<ELT>::TReal >( // this is just ELT
            MatrixHelper(helper,0,0,nrow(),ncol()));    // a view of the whole matrix
    }
    // Elements are complex -- helper uses underlying precision (real) type.
    MatrixHelper<Precision> h(helper,typename MatrixHelper<Precision>::RealView);    
    return MatrixView_< typename CNT<ELT>::TReal >(h); 
}
*/


//  ----------------------------------------------------------------------------
/// @name Global operators involving Matrix objects
/// These operators take MatrixBase arguments and produce Matrix_ results.
/// @{

// + and - allow mixed element types, but will fail to compile if the elements aren't
// compatible. At run time these will fail if the dimensions are incompatible.
template <class E1, class E2>
Matrix_<typename CNT<E1>::template Result<E2>::Add>
operator+(const MatrixBase<E1>& l, const MatrixBase<E2>& r) {
    return Matrix_<typename CNT<E1>::template Result<E2>::Add>(l) += r;
}

template <class E>
Matrix_<E> operator+(const MatrixBase<E>& l, const typename CNT<E>::T& r) {
    return Matrix_<E>(l) += r;
}

template <class E>
Matrix_<E> operator+(const typename CNT<E>::T& l, const MatrixBase<E>& r) {
    return Matrix_<E>(r) += l;
}

template <class E1, class E2>
Matrix_<typename CNT<E1>::template Result<E2>::Sub>
operator-(const MatrixBase<E1>& l, const MatrixBase<E2>& r) {
    return Matrix_<typename CNT<E1>::template Result<E2>::Sub>(l) -= r;
}

template <class E>
Matrix_<E> operator-(const MatrixBase<E>& l, const typename CNT<E>::T& r) {
    return Matrix_<E>(l) -= r;
}

template <class E>
Matrix_<E> operator-(const typename CNT<E>::T& l, const MatrixBase<E>& r) {
    Matrix_<E> temp(r.nrow(), r.ncol());
    temp = l;
    return (temp -= r);
}

// Scalar multiply and divide. You might wish the scalar could be
// a templatized type "E2", but that would create horrible ambiguities since
// E2 would match not only scalar types but everything else including
// matrices.
template <class E> Matrix_<E>
operator*(const MatrixBase<E>& l, const typename CNT<E>::StdNumber& r) 
  { return Matrix_<E>(l)*=r; }

template <class E> Matrix_<E>
operator*(const typename CNT<E>::StdNumber& l, const MatrixBase<E>& r) 
  { return Matrix_<E>(r)*=l; }

template <class E> Matrix_<E>
operator/(const MatrixBase<E>& l, const typename CNT<E>::StdNumber& r) 
  { return Matrix_<E>(l)/=r; }

// Handle ints explicitly.
template <class E> Matrix_<E>
operator*(const MatrixBase<E>& l, int r) 
  { return Matrix_<E>(l)*= typename CNT<E>::StdNumber(r); }

template <class E> Matrix_<E>
operator*(int l, const MatrixBase<E>& r) 
  { return Matrix_<E>(r)*= typename CNT<E>::StdNumber(l); }

template <class E> Matrix_<E>
operator/(const MatrixBase<E>& l, int r) 
  { return Matrix_<E>(l)/= typename CNT<E>::StdNumber(r); }

/// @}

/// @name Global operators involving Vector objects
/// These operators take VectorBase arguments and produce Vector_ results.
/// @{

template <class E1, class E2>
Vector_<typename CNT<E1>::template Result<E2>::Add>
operator+(const VectorBase<E1>& l, const VectorBase<E2>& r) {
    return Vector_<typename CNT<E1>::template Result<E2>::Add>(l) += r;
}
template <class E>
Vector_<E> operator+(const VectorBase<E>& l, const typename CNT<E>::T& r) {
    return Vector_<E>(l) += r;
}
template <class E>
Vector_<E> operator+(const typename CNT<E>::T& l, const VectorBase<E>& r) {
    return Vector_<E>(r) += l;
}
template <class E1, class E2>
Vector_<typename CNT<E1>::template Result<E2>::Sub>
operator-(const VectorBase<E1>& l, const VectorBase<E2>& r) {
    return Vector_<typename CNT<E1>::template Result<E2>::Sub>(l) -= r;
}
template <class E>
Vector_<E> operator-(const VectorBase<E>& l, const typename CNT<E>::T& r) {
    return Vector_<E>(l) -= r;
}
template <class E>
Vector_<E> operator-(const typename CNT<E>::T& l, const VectorBase<E>& r) {
    Vector_<E> temp(r.size());
    temp = l;
    return (temp -= r);
}

// Scalar multiply and divide.

template <class E> Vector_<E>
operator*(const VectorBase<E>& l, const typename CNT<E>::StdNumber& r) 
  { return Vector_<E>(l)*=r; }

template <class E> Vector_<E>
operator*(const typename CNT<E>::StdNumber& l, const VectorBase<E>& r) 
  { return Vector_<E>(r)*=l; }

template <class E> Vector_<E>
operator/(const VectorBase<E>& l, const typename CNT<E>::StdNumber& r) 
  { return Vector_<E>(l)/=r; }

// Handle ints explicitly
template <class E> Vector_<E>
operator*(const VectorBase<E>& l, int r) 
  { return Vector_<E>(l)*= typename CNT<E>::StdNumber(r); }

template <class E> Vector_<E>
operator*(int l, const VectorBase<E>& r) 
  { return Vector_<E>(r)*= typename CNT<E>::StdNumber(l); }

template <class E> Vector_<E>
operator/(const VectorBase<E>& l, int r) 
  { return Vector_<E>(l)/= typename CNT<E>::StdNumber(r); }

// These are fancier "scalars"; whether they are allowed depends on
// whether the element type and the CNT are compatible.

// Vector * Vec
template <class E1, int M, class E2, int S> 
Vector_<typename CNT<E1>::template Result< Vec<M,E2,S> >::Mul>
operator*(const VectorBase<E1>& v, const Vec<M,E2,S>& s) {
    Vector_<typename CNT<E1>::template Result< Vec<M,E2,S> >::Mul> res(v.nrow());
    for (int i=0; i < v.nrow(); ++i)
        res[i] = v[i]*s; 
    return res;
}

// Vec * Vector
template <class E1, int M, class E2, int S> 
Vector_<typename Vec<M,E2,S>::template Result<E1>::Mul>
operator*(const Vec<M,E2,S>& s, const VectorBase<E1>& v) {
    Vector_<typename Vec<M,E2,S>::template Result<E1>::Mul> res(v.nrow());
    for (int i=0; i < v.nrow(); ++i)
        res[i] = s*v[i]; 
    return res;
}

// Vector * Row
template <class E1, int N, class E2, int S> 
Vector_<typename CNT<E1>::template Result< Row<N,E2,S> >::Mul>
operator*(const VectorBase<E1>& v, const Row<N,E2,S>& s) {
    Vector_<typename CNT<E1>::template Result< Row<N,E2,S> >::Mul> res(v.nrow());
    for (int i=0; i < v.nrow(); ++i)
        res[i] = v[i]*s; 
    return res;
}

// Row * Vector
template <class E1, int N, class E2, int S> 
Vector_<typename Row<N,E2,S>::template Result<E1>::Mul>
operator*(const Row<N,E2,S>& s, const VectorBase<E1>& v) {
    Vector_<typename Row<N,E2,S>::template Result<E1>::Mul> res(v.nrow());
    for (int i=0; i < v.nrow(); ++i)
        res[i] = s*v[i]; 
    return res;
}

// Vector * Mat
template <class E1, int M, int N, class E2, int S1, int S2> 
Vector_<typename CNT<E1>::template Result< Mat<M,N,E2,S1,S2> >::Mul>
operator*(const VectorBase<E1>& v, const Mat<M,N,E2,S1,S2>& s) {
    Vector_<typename CNT<E1>::template Result< Mat<M,N,E2,S1,S2> >::Mul> res(v.nrow());
    for (int i=0; i < v.nrow(); ++i)
        res[i] = v[i]*s; 
    return res;
}

// Mat * Vector
template <class E1, int M, int N, class E2, int S1, int S2> 
Vector_<typename Mat<M,N,E2,S1,S2>::template Result<E1>::Mul>
operator*(const Mat<M,N,E2,S1,S2>& s, const VectorBase<E1>& v) {
    Vector_<typename Mat<M,N,E2,S1,S2>::template Result<E1>::Mul> res(v.nrow());
    for (int i=0; i < v.nrow(); ++i)
        res[i] = s*v[i]; 
    return res;
}

// Vector * SymMat
template <class E1, int M, class E2, int S> 
Vector_<typename CNT<E1>::template Result< SymMat<M,E2,S> >::Mul>
operator*(const VectorBase<E1>& v, const SymMat<M,E2,S>& s) {
    Vector_<typename CNT<E1>::template Result< SymMat<M,E2,S> >::Mul> res(v.nrow());
    for (int i=0; i < v.nrow(); ++i)
        res[i] = v[i]*s; 
    return res;
}

// SymMat * Vector
template <class E1, int M, class E2, int S> 
Vector_<typename SymMat<M,E2,S>::template Result<E1>::Mul>
operator*(const SymMat<M,E2,S>& s, const VectorBase<E1>& v) {
    Vector_<typename SymMat<M,E2,S>::template Result<E1>::Mul> res(v.nrow());
    for (int i=0; i < v.nrow(); ++i)
        res[i] = s*v[i]; 
    return res;
}

/// @}

/// @name Global operators involving RowVector objects
/// These operators take RowVectorBase arguments and produce RowVector_ results.
/// @{

template <class E1, class E2>
RowVector_<typename CNT<E1>::template Result<E2>::Add>
operator+(const RowVectorBase<E1>& l, const RowVectorBase<E2>& r) {
    return RowVector_<typename CNT<E1>::template Result<E2>::Add>(l) += r;
}
template <class E>
RowVector_<E> operator+(const RowVectorBase<E>& l, const typename CNT<E>::T& r) {
    return RowVector_<E>(l) += r;
}
template <class E>
RowVector_<E> operator+(const typename CNT<E>::T& l, const RowVectorBase<E>& r) {
    return RowVector_<E>(r) += l;
}
template <class E1, class E2>
RowVector_<typename CNT<E1>::template Result<E2>::Sub>
operator-(const RowVectorBase<E1>& l, const RowVectorBase<E2>& r) {
    return RowVector_<typename CNT<E1>::template Result<E2>::Sub>(l) -= r;
}
template <class E>
RowVector_<E> operator-(const RowVectorBase<E>& l, const typename CNT<E>::T& r) {
    return RowVector_<E>(l) -= r;
}
template <class E>
RowVector_<E> operator-(const typename CNT<E>::T& l, const RowVectorBase<E>& r) {
    RowVector_<E> temp(r.size());
    temp = l;
    return (temp -= r);
}

// Scalar multiply and divide 

template <class E> RowVector_<E>
operator*(const RowVectorBase<E>& l, const typename CNT<E>::StdNumber& r) 
  { return RowVector_<E>(l)*=r; }

template <class E> RowVector_<E>
operator*(const typename CNT<E>::StdNumber& l, const RowVectorBase<E>& r) 
  { return RowVector_<E>(r)*=l; }

template <class E> RowVector_<E>
operator/(const RowVectorBase<E>& l, const typename CNT<E>::StdNumber& r) 
  { return RowVector_<E>(l)/=r; }

// Handle ints explicitly.
template <class E> RowVector_<E>
operator*(const RowVectorBase<E>& l, int r) 
  { return RowVector_<E>(l)*= typename CNT<E>::StdNumber(r); }

template <class E> RowVector_<E>
operator*(int l, const RowVectorBase<E>& r) 
  { return RowVector_<E>(r)*= typename CNT<E>::StdNumber(l); }

template <class E> RowVector_<E>
operator/(const RowVectorBase<E>& l, int r) 
  { return RowVector_<E>(l)/= typename CNT<E>::StdNumber(r); }


// These are fancier "scalars"; whether they are allowed depends on
// whether the element type and the CNT are compatible.

// RowVector * Vec
template <class E1, int M, class E2, int S> 
RowVector_<typename CNT<E1>::template Result< Vec<M,E2,S> >::Mul>
operator*(const RowVectorBase<E1>& v, const Vec<M,E2,S>& s) {
    RowVector_<typename CNT<E1>::template Result< Vec<M,E2,S> >::Mul> res(v.ncol());
    for (int i=0; i < v.ncol(); ++i)
        res[i] = v[i]*s; 
    return res;
}

// Vec * RowVector
template <class E1, int M, class E2, int S> 
RowVector_<typename Vec<M,E2,S>::template Result<E1>::Mul>
operator*(const Vec<M,E2,S>& s, const RowVectorBase<E1>& v) {
    RowVector_<typename Vec<M,E2,S>::template Result<E1>::Mul> res(v.ncol());
    for (int i=0; i < v.ncol(); ++i)
        res[i] = s*v[i]; 
    return res;
}

// RowVector * Row
template <class E1, int N, class E2, int S> 
RowVector_<typename CNT<E1>::template Result< Row<N,E2,S> >::Mul>
operator*(const RowVectorBase<E1>& v, const Row<N,E2,S>& s) {
    RowVector_<typename CNT<E1>::template Result< Row<N,E2,S> >::Mul> res(v.ncol());
    for (int i=0; i < v.ncol(); ++i)
        res[i] = v[i]*s; 
    return res;
}

// Row * RowVector
template <class E1, int N, class E2, int S> 
RowVector_<typename Row<N,E2,S>::template Result<E1>::Mul>
operator*(const Row<N,E2,S>& s, const RowVectorBase<E1>& v) {
    RowVector_<typename Row<N,E2,S>::template Result<E1>::Mul> res(v.ncol());
    for (int i=0; i < v.ncol(); ++i)
        res[i] = s*v[i]; 
    return res;
}

// RowVector * Mat
template <class E1, int M, int N, class E2, int S1, int S2> 
RowVector_<typename CNT<E1>::template Result< Mat<M,N,E2,S1,S2> >::Mul>
operator*(const RowVectorBase<E1>& v, const Mat<M,N,E2,S1,S2>& s) {
    RowVector_<typename CNT<E1>::template Result< Mat<M,N,E2,S1,S2> >::Mul> res(v.ncol());
    for (int i=0; i < v.ncol(); ++i)
        res[i] = v[i]*s; 
    return res;
}

// Mat * RowVector
template <class E1, int M, int N, class E2, int S1, int S2> 
RowVector_<typename Mat<M,N,E2,S1,S2>::template Result<E1>::Mul>
operator*(const Mat<M,N,E2,S1,S2>& s, const RowVectorBase<E1>& v) {
    RowVector_<typename Mat<M,N,E2,S1,S2>::template Result<E1>::Mul> res(v.ncol());
    for (int i=0; i < v.ncol(); ++i)
        res[i] = s*v[i]; 
    return res;
}

// RowVector * SymMat
template <class E1, int M, class E2, int S> 
RowVector_<typename CNT<E1>::template Result< SymMat<M,E2,S> >::Mul>
operator*(const RowVectorBase<E1>& v, const SymMat<M,E2,S>& s) {
    RowVector_<typename CNT<E1>::template Result< SymMat<M,E2,S> >::Mul> res(v.ncol());
    for (int i=0; i < v.ncol(); ++i)
        res[i] = v[i]*s; 
    return res;
}

// SymMat * RowVector
template <class E1, int M, class E2, int S> 
RowVector_<typename SymMat<M,E2,S>::template Result<E1>::Mul>
operator*(const SymMat<M,E2,S>& s, const RowVectorBase<E1>& v) {
    RowVector_<typename SymMat<M,E2,S>::template Result<E1>::Mul> res(v.ncol());
    for (int i=0; i < v.ncol(); ++i)
        res[i] = s*v[i]; 
    return res;
}

/// @}


/// @name Global operators involving mixed matrix, vector, and row vector objects
/// These operators take MatrixBase, VectorBase, and RowVectorBase arguments
/// and produce Matrix_, Vector_, and RowVector_ results.
/// @{

    // TODO: these should use LAPACK!

// Dot product
template <class E1, class E2> 
typename CNT<E1>::template Result<E2>::Mul
operator*(const RowVectorBase<E1>& r, const VectorBase<E2>& v) {
    assert(r.ncol() == v.nrow());
    typename CNT<E1>::template Result<E2>::Mul sum(0);
    for (int j=0; j < r.ncol(); ++j)
        sum += r(j) * v[j];
    return sum;
}

template <class E1, class E2> 
Vector_<typename CNT<E1>::template Result<E2>::Mul>
operator*(const MatrixBase<E1>& m, const VectorBase<E2>& v) {
    assert(m.ncol() == v.nrow());
    Vector_<typename CNT<E1>::template Result<E2>::Mul> res(m.nrow());
    for (int i=0; i< m.nrow(); ++i)
        res[i] = m[i]*v;
    return res;
}

template <class E1, class E2> 
Matrix_<typename CNT<E1>::template Result<E2>::Mul>
operator*(const MatrixBase<E1>& m1, const MatrixBase<E2>& m2) {
    assert(m1.ncol() == m2.nrow());
    Matrix_<typename CNT<E1>::template Result<E2>::Mul> 
        res(m1.nrow(),m2.ncol());

    for (int j=0; j < res.ncol(); ++j)
        for (int i=0; i < res.nrow(); ++i)
            res(i,j) = m1[i] * m2(j);

    return res;
}

/// @}

// This "private" static method is used to implement VectorView's 
// fillVectorViewFromStream() and Vector's readVectorFromStream() 
// namespace-scope static methods, which are in turn used to implement 
// VectorView's and 
// Vector's stream extraction operators ">>". This method has to be in the 
// header file so that we don't need to pass streams through the API, but it 
// is not intended for use by users and has no Doxygen presence, unlike 
// fillArrayFromStream() and readArrayFromStream() and (more commonly)
// the extraction operators.
template <class T> static inline 
std::istream& readVectorFromStreamHelper
   (std::istream& in, bool isFixedSize, Vector_<T>& out)
{
    // If already failed, bad, or eof, set failed bit and return without 
    // touching the Vector.
    if (!in.good()) {in.setstate(std::ios::failbit); return in;}

    // If the passed-in Vector isn't resizeable, then we have to treat it as
    // a fixed size VectorView regardless of the setting of the isFixedSize
    // argument.
    if (!out.isResizeable())
        isFixedSize = true; // might be overriding the argument here

    // numRequired will be ignored unless isFixedSize==true.
    const int numRequired = isFixedSize ? out.size() : 0;

    if (!isFixedSize)
        out.clear(); // We're going to replace the entire contents of the Array.

    // Skip initial whitespace. If that results in eof this may be a successful
    // read of a 0-length, unbracketed Vector. That is OK for either a
    // variable-length Vector or a fixed-length VectorView of length zero.
    std::ws(in); if (in.fail()) return in;
    if (in.eof()) {
        if (isFixedSize && numRequired != 0)
            in.setstate(std::ios_base::failbit); // zero elements not OK
        return in;
    }
    
    // Here the stream is good and the next character is non-white.
    assert(in.good());

    // Use this for raw i/o (peeks and gets).
    typename       std::iostream::int_type ch;
#ifndef NDEBUG
    const typename std::iostream::int_type EOFch = 
        std::iostream::traits_type::eof();
#endif

    // First we'll look for the optional "~". If found, the brackets become
    // required.
    bool tildeFound = false;
    ch = in.peek(); if (in.fail()) return in;
    assert(ch != EOFch); // we already checked above
    if ((char)ch == '~') {
        tildeFound = true;
        in.get(); // absorb the tilde
        // Eat whitespace after the tilde to see what's next.
        if (in.good()) std::ws(in);
        // If we hit eof after the tilde we don't like the formatting.
        if (!in.good()) {in.setstate(std::ios_base::failbit); return in;}
    }

    // Here the stream is good, the next character is non-white, and we
    // might have seen a tilde.
    assert(in.good());

    // Now see if the sequence is bare or surrounded by (), or [].
    bool lookForCloser = true;
    char openBracket, closeBracket;
    ch = in.peek(); if (in.fail()) return in;
    assert(ch != EOFch); // we already checked above

    openBracket = (char)ch;
    if      (openBracket=='(') {in.get(); closeBracket = ')';}
    else if (openBracket=='[') {in.get(); closeBracket = ']';}
    else lookForCloser = false;

    // If we found a tilde, the opening bracket was mandatory. If we didn't
    // find one then we reject the formatting.
    if (tildeFound && !lookForCloser)
    {   in.setstate(std::ios_base::failbit); return in;}

    // If lookForCloser is true, then closeBracket contains the terminating
    // delimiter, otherwise we're not going to quit until eof.

    // Eat whitespace after the opening bracket to see what's next.
    if (in.good()) std::ws(in);

    // If we're at eof now it must be because the open bracket was the
    // last non-white character in the stream, which is an error.
    if (!in.good()) {
        if (in.eof()) {
            assert(lookForCloser); // or we haven't read anything that could eof
            in.setstate(std::ios::failbit);
        }
        return in;
    }

    // istream is good and next character is non-white; ready to read first
    // value or terminator.

    // We need to figure out whether the elements are space- or comma-
    // separated and then insist on consistency.
    bool commaOK = true, commaRequired = false;
    bool terminatorSeen = false;
    int nextIndex = 0;
    while (true) {
        char c;

        // Here at the top of this loop, we have already successfully read 
        // n=nextIndex values of type T. For fixed-size reads, it might be
        // the case that n==numRequired already, but we still may need to
        // look for a closing bracket before we can declare victory.
        // The stream is good() (not at eof) but it might be the case that 
        // there is nothing but white space left; we don't know yet because
        // if we have satisfied the fixed-size count and are not expecting
        // a terminator then we should quit without absorbing the trailing
        // white space.
        assert(in.good());

        // Look for closing bracket before trying to read value.
        if (lookForCloser) {
            // Eat white space to find the closing bracket.
            std::ws(in); if (!in.good()) break; // eof?
            ch = in.peek(); assert(ch != EOFch);
            if (!in.good()) break;
            c = (char)ch;
            if (c == closeBracket) {   
                in.get(); // absorb the closing bracket
                terminatorSeen = true; 
                break; 
            }
            // next char not a closing bracket; fall through
        }

        // We didn't look or didn't find a closing bracket. The istream is good 
        // but we might be looking at white space.

        // If we already got all the elements we want, break for final checks.
        if (isFixedSize && (nextIndex == numRequired))
            break; // that's a full count.

        // Look for comma before value, except the first time.
        if (commaOK && nextIndex != 0) {
            // Eat white space to find the comma.
            std::ws(in); if (!in.good()) break; // eof?
            ch = in.peek(); assert(ch != EOFch);
            if (!in.good()) break;
            c = (char)ch;
            if (c == ',') {
                in.get(); // absorb comma
                commaRequired = true; // all commas from now on
            } else { // next char not a comma
                if (commaRequired) // bad, e.g.: v1, v2, v3 v4 
                {   in.setstate(std::ios::failbit); break; }
                else commaOK = false; // saw: v1 v2 (no commas now)
            }
            if (!in.good()) break; // might be eof
        }

        // No closing bracket yet; don't have enough elements; skipped comma 
        // if any; istream is good; might be looking at white space.
        assert(in.good());

        // Now read in an element of type T.
        // The extractor T::operator>>() will ignore leading white space.
        if (!isFixedSize)
            out.resizeKeep(out.size()+1); // grow by one (default consructed)
        in >> out[nextIndex]; if (in.fail()) break;
        ++nextIndex;

        if (!in.good()) break; // might be eof
    }

    // We will get here under a number of circumstances:
    //  - the fail bit is set in the istream, or
    //  - we reached eof
    //  - we saw a closing brace
    //  - we got all the elements we wanted (for a fixed-size read)
    // Note that it is possible that we consumed everything except some
    // trailing white space (meaning we're not technically at eof), but
    // for consistency with built-in operator>>()'s we won't try to absorb
    // that trailing white space.

    if (!in.fail()) {
        if (lookForCloser && !terminatorSeen)
            in.setstate(std::ios::failbit); // missing terminator

        if (isFixedSize && nextIndex != numRequired)
            in.setstate(std::ios::failbit); // wrong number of values
    }

    return in;
}



//------------------------------------------------------------------------------
//                          RELATED GLOBAL OPERATORS
//------------------------------------------------------------------------------
// These are logically part of the Matrix_<T> class but are not actually 
// class members; that is, they are in the SimTK namespace.

/**@name             Matrix_<T> serialization and I/O
These methods are at namespace scope but are logically part of the Vector
classes. These deal with reading and writing Vectors from and to streams,
which places an additional requirement on the element type T: the element 
must support the same operation you are trying to do on the Vector as a 
whole. **/
/**@{**/

/** Specialize for VectorBase<E> to delegate to element type E, with spaces
separating the elements. 
@relates SimTK::VectorBase **/
template <class E> inline void
writeUnformatted(std::ostream& o, const VectorBase<E>& v) {
    const int sz = v.size();
    for (int i=0; i < sz; ++i) {
        if (i != 0) o << " ";
        writeUnformatted(o, v[i]);
    }
} 
/** Raw serialization of VectorView_<E>; same as VectorBase<E>.
@relates SimTK::VectorView_ **/
template <class E> inline void
writeUnformatted(std::ostream& o, const VectorView_<E>& v) 
{   writeUnformatted(o, static_cast< const VectorBase<E> >(v)); }

/** Raw serialization of Vector_<E>; same as VectorBase<E>.
@relates SimTK::Vector_ **/
template <class E> inline void
writeUnformatted(std::ostream& o, const Vector_<E>& v) 
{   writeUnformatted(o, static_cast< const VectorBase<E> >(v)); }

/** Specialize for RowVectorBase<E> to delegate to element type E, with spaces
separating the elements; raw output is same as VectorBase. 
@relates SimTK::RowVectorBase **/
template <class E> inline void
writeUnformatted(std::ostream& o, const RowVectorBase<E>& v) 
{   writeUnformatted(o, ~v); }

/** Raw serialization of RowVectorView_<E>; same as VectorView_<E>.
@relates SimTK::RowVectorView_ **/
template <class E> inline void
writeUnformatted(std::ostream& o, const RowVectorView_<E>& v) 
{   writeUnformatted(o, static_cast< const RowVectorBase<E> >(v)); }

/** Raw serialization of RowVector_<E>; same as Vector_<E>.
@relates SimTK::RowVector_ **/
template <class E> inline void
writeUnformatted(std::ostream& o, const RowVector_<E>& v) 
{   writeUnformatted(o, static_cast< const RowVectorBase<E> >(v)); }

/** Specialize for MatrixBase<E> delegating to RowVectorBase<E> with newlines
separating the rows, but no final newline.
@relates SimTK::MatrixBase **/
template <class E> inline void
writeUnformatted(std::ostream& o, const MatrixBase<E>& v) {
    const int nr = v.nrow();
    for (int i=0; i < nr; ++i) {
        if (i != 0) o << std::endl;
        writeUnformatted(o, v[i]);
    }
}
/** Raw serialization of MatrixView_<E>; same as MatrixBase<E>.
@relates SimTK::MatrixView_ **/
template <class E> inline void
writeUnformatted(std::ostream& o, const MatrixView_<E>& v) 
{   writeUnformatted(o, static_cast< const MatrixBase<E> >(v)); }

/** Raw serialization of Vector_<E>; same as VectorBase<E>.
@relates SimTK::Matrix_ **/
template <class E> inline void
writeUnformatted(std::ostream& o, const Matrix_<E>& v) 
{   writeUnformatted(o, static_cast< const MatrixBase<E> >(v)); }

/** Read fixed-size VectorView from input stream. It is an error if there
aren't enough elements. 
@relates SimTK::VectorView_ **/
template <class E> inline bool
readUnformatted(std::istream& in, VectorView_<E>& v) {
    for (int i=0; i < v.size(); ++i)
        if (!readUnformatted(in, v[i])) return false;
    return true;
}

/** Read variable-size Vector from input stream. Reads until error or eof. 
@relates SimTK::Vector_ **/
template <class E> inline bool
readUnformatted(std::istream& in, Vector_<E>& v) {
    if (!v.isResizeable())
        return readUnformatted(in, v.updAsVectorView());

    Array_<E,int> a;
    if (!readUnformatted(in,a)) return false;
    v.resize(a.size());
    for (int i=0; i<a.size(); ++i)
        v[i] = a[i];
    return true;
}

/** Read fixed-size RowVectorView from input stream. It is an error if there
aren't enough elements.  
@relates SimTK::RowVectorView_ **/
template <class E> inline bool
readUnformatted(std::istream& in, RowVectorView_<E>& v) 
{   VectorView_<E> vt(~v);
    return readUnformatted<E>(in, vt); }

/** Read variable-size RowVector from unformatted (whitespace-separated) 
input stream. Reads until error or eof. 
@relates SimTK::RowVector_ **/
template <class E> inline bool
readUnformatted(std::istream& in, RowVector_<E>& v) 
{   Vector_<E> vt(~v);
    return readUnformatted<E>(in, vt); }

/** Read fixed-size MatrixView in row order from unformatted (whitespace-
separated) input stream. Newlines in the input have no special meaning --
we'll read them as whitespace. It is an error if there aren't enough 
elements.  
@relates SimTK::MatrixView_ **/
template <class E> inline bool
readUnformatted(std::istream& in, MatrixView_<E>& v) { 
    for (int row=0; row < v.nrow(); ++row) {
        RowVectorView_<E> oneRow(v[row]);
        if (!readUnformatted<E>(in, oneRow)) return false;
    }
    return true;
}

/** Read in new values for a Matrix without changing its size, from a stream
of whitespace-separated tokens with no other formatting recognized. Newlines in 
the input have no special meaning -- we'll read them as whitespace. It is an 
error if there aren't enough elements.  
@relates SimTK::Matrix_ **/
template <class E> inline bool
fillUnformatted(std::istream& in, Matrix_<E>& v) {
    return readUnformatted<E>(in, v.updAsMatrixView());
}

/** NOT IMPLEMENTED: read variable-size Matrix recognizing newlines as end
of row; use fillUnformatted() instead.  
@relates SimTK::Matrix_ **/
template <class E> inline bool
readUnformatted(std::istream& in, Matrix_<E>& v) {
    SimTK_ASSERT_ALWAYS(!"implemented", 
        "SimTK::readUnformatted(istream, Matrix) is not implemented; try"
        " SimTK::fillUnformatted(istream, Matrix) instead.");
    return false;
}
  
/** Output a human readable representation of a Vector to an std::ostream
(like std::cout). The format is ~[ \e elements ] where \e elements is a 
space-separated list of the Vector's contents output by invoking the "<<" 
operator on the elements. This function will not compile if the element type 
does not support the "<<" operator. No newline is issued before
or after the output. @relates Vector_ **/
template <class T> inline std::ostream&
operator<<(std::ostream& o, const VectorBase<T>& v)
{   o << "~["; 
    if (v.size()) {
        o << v[0];
        for (int i=1; i < v.size(); ++i) o << " " << v[i];
    }
    return o << "]"; 
}

/** Output a human readable representation of a RowVector to an std::ostream
(like std::cout). The format is [ \e elements ] where \e elements is a 
space-separated list of the RowVector's contents output by invoking the "<<" 
operator on the elements. This function will not compile if the element type 
does not support the "<<" operator. No newline is issued before
or after the output. @relates RowVector_ **/
template <class T> inline std::ostream&
operator<<(std::ostream& o, const RowVectorBase<T>& v)
{   o << "["; 
    if (v.size()) {
        o << v[0];
        for (int i=1; i < v.size(); ++i) o << " " << v[i];
    }
    return o << "]"; 
}

/** Output a human readable representation of a Matrix to an std::ostream
(like std::cout). The format is one row per line, with each row output as
[ \e elements ] where \e elements is a 
space-separated list of the row's contents output by invoking the "<<" 
operator on the elements. This function will not compile if the element type 
does not support the "<<" operator. A newline is issued before each row and
at the end. @relates Matrix_ **/
template <class T> inline std::ostream&
operator<<(std::ostream& o, const MatrixBase<T>& m) {
    for (int i=0;i<m.nrow();++i)
        o << std::endl << m[i];
    if (m.nrow()) o << std::endl;
    return o; 
}


/** Read in a Vector_<T> from a stream, as a sequence of space-separated or
comma-separated values optionally surrounded by parentheses (), or square 
brackets [], or the "transposed" ~() or ~[]. In the case that the transpose
operator is present, the parentheses or brackets are required, otherwise they
are optional. We will continue to read elements of 
type T from the stream until we find a reason to stop, using type T's stream 
extraction operator>>() to read in each element and resizing the Vector as
necessary. If the data is bracketed, we'll read until we hit the closing 
bracket. If it is not bracketed, we'll read until we hit eof() or get an error
such as the element extractor setting the stream's fail bit due to bad 
formatting. On successful return, the stream will be positioned right after 
the final read-in element or terminating bracket, and the stream's status will 
be good() or eof(). We will not consume trailing whitespace after bracketed 
elements; that means the stream might actually be empty even if we don't 
return eof(). If you want to know whether there is anything else in the 
stream, follow this call with the STL whitespace skipper std::ws() like this:
@code
    if (readVectorFromStream(in,vec) && !in.eof()) 
        std::ws(in); // might take us to eof
    if (in.fail()) {...} // probably a formatting error
    else {
        // Here if the stream is good() then there is more to read; if the
        // stream got used up the status is guaranteed to be eof().
    }
@endcode
A compilation error will occur if you try to use this method on an Vector_<T>
for a type T for which there is no stream extraction operator>>(). 
@note If you want to fill a resizeable Vector_<T> with a fixed amount of data 
from the stream, resize() the Vector to the appropriate length and then use 
fillVectorFromStream() instead. @see fillVectorFromStream()
@relates Vector_ **/
template <class T> static inline 
std::istream& readVectorFromStream(std::istream& in, Vector_<T>& out)
{   return readVectorFromStreamHelper<T>(in, false /*variable sizez*/, out); }



/** Read in a fixed number of elements from a stream into a Vector. We expect 
to read in exactly size() elements of type T, using type T's stream extraction 
operator>>(). This will stop reading when we've read size() elements, or set 
the fail bit in the stream if we run out of elements or if any element's 
extract operator sets the fail bit. On successful return, all size() elements 
will have been set, the stream will be positioned right after the final 
read-in element or terminating bracket, and the stream's status will be good()
or eof(). We will not consume trailing whitespace after reading all the 
elements; that means the stream might actually be empty even if we don't 
return eof(). If you want to know whether there is anything else in the 
stream, follow this call with std::ws() like this:
@code
    if (fillVectorFromStream(in,vec))
        if (!in.eof()) std::ws(in); // might take us to eof
    if (in.fail()) {...} // deal with I/O or formatting error
    // Here if the stream is good() then there is more to read; if the
    // stream got used up the status is guaranteed to be eof().
@endcode
A compilation error will occur if you try to use this method on a Vector_<T>
for a type T for which there is no stream extraction operator>>().
@note If you want to read in a variable number of elements and have the 
Vector_<T> resized as needed, use readVectorFromStream() instead.
@see readVectorFromStream()
@relates Vector_ **/
template <class T> static inline 
std::istream& fillVectorFromStream(std::istream& in, Vector_<T>& out)
{   return readVectorFromStreamHelper<T>(in, true /*fixed size*/, out); }

/** Read in a fixed number of elements from a stream into an VectorView. See
fillVectorFromStream() for more information; this works the same way.
@see fillVectorFromStream()
@relates VectorView_ **/
template <class T> static inline 
std::istream& fillVectorViewFromStream(std::istream& in, VectorView_<T>& out)
{   return readVectorFromStreamHelper<T>(in, true /*fixed size*/, out); }


/** Read Vector_<T> from a stream as a sequence of space- or comma-separated
values of type T, optionally delimited by parentheses, or brackets, and
preceded by "~". The Vector_<T> may be an owner (variable size) or a view 
(fixed size n). In the case of an owner, we'll read all the elements in 
brackets or until eof if there are no brackets. In the case of a view, there 
must be exactly n elements in brackets, or if there are no brackets we'll 
consume exactly n elements and then stop. Each element is read in with its 
own operator ">>" so this won't work if no such operator is defined for 
type T. @relates Vector_ **/
template <class T> inline
std::istream& operator>>(std::istream& in, Vector_<T>& out) 
{   return readVectorFromStream<T>(in, out); }

/** Read a (fixed size n) VectorView_<T> from a stream as a sequence of space- 
or comma-separated values of type T, optionally delimited by parentheses or 
square brackets, and preceded by "~". If there are no delimiters then we will 
read size() values and then stop. Otherwise, there must be exactly size() 
values within the brackets. Each element is read in with its own 
operator ">>" so  this won't work if no such operator is defined for type T.
@relates VectorView_ **/
template <class T> inline
std::istream& operator>>(std::istream& in, VectorView_<T>& out) 
{   return fillVectorViewFromStream<T>(in, out); }
/**@}**/  // End of Matrix serialization.

// Friendly abbreviations for vectors and matrices with scalar elements.
/** @ingroup MatVecTypedefs **/
/**@{**/
/** Variable-size column vector of Real elements; abbreviation for 
Vector_<Real>. This is the most common large-vector type in the Simbody API
and in Simbody user programs. **/
typedef Vector_<Real>           Vector;

/** Variable-size 2D matrix of Real elements; abbreviation for Matrix_<Real>.
This is the most common large-matrix type in the Simbody API and in Simbody
user programs. **/
typedef Matrix_<Real>           Matrix;

/** Variable-size row vector of Real elements; abbreviation for
RowVector_<Real>. This is the type of a transposed Vector and does not 
usually appear explicitly in the Simbody API or user programs. **/
typedef RowVector_<Real>        RowVector;
/**@}**/

/** @ingroup UncommonMatVecTypedefs **/
/**@{**/
/** Variable-size column vector of Complex (std::complex<Real>) elements. **/
typedef Vector_<Complex>        ComplexVector;
/** Variable-size 2D matrix of Complex (std::complex<Real>) elements. **/
typedef Matrix_<Complex>        ComplexMatrix;
/** Variable-size row vector of Complex (std::complex<Real>) elements. **/
typedef RowVector_<Complex>     ComplexRowVector;

/** Non-owner column vector sharing Real elements. **/
typedef VectorView_<Real>       VectorView;
/** Non-owner matrix sharing Real elements. **/
typedef MatrixView_<Real>       MatrixView;
/** Non-owner row vector sharing Real elements. **/
typedef RowVectorView_<Real>    RowVectorView;

/** Non-owner column vector sharing Complex (std::complex<Real>) elements. **/
typedef VectorView_<Complex>    ComplexVectorView;
/** Non-owner matrix sharing Complex (std::complex<Real>) elements. **/
typedef MatrixView_<Complex>    ComplexMatrixView;
/** Non-owner row vector sharing Complex (std::complex<Real>) elements. **/
typedef RowVectorView_<Complex> ComplexRowVectorView;

/** Abbreviation for Vector_<float>. **/
typedef Vector_<float>          fVector;
/** Abbreviation for Vector_<double>. **/
typedef Vector_<double>         dVector;
/** Abbreviation for Vector_<std::complex<float>>. **/
typedef Vector_<fComplex>       fComplexVector;
/** Abbreviation for Vector_<std::complex<double>>. **/
typedef Vector_<dComplex>       dComplexVector;

/** Abbreviation for VectorView_<float>. **/
typedef VectorView_<float>      fVectorView;
/** Abbreviation for VectorView_<double>. **/
typedef VectorView_<double>     dVectorView;
/** Abbreviation for VectorView_<std::complex<float>>. **/
typedef VectorView_<fComplex>   fComplexVectorView;
/** Abbreviation for VectorView_<std::complex<double>>. **/
typedef VectorView_<dComplex>   dComplexVectorView;

/** Abbreviation for RowVector_<float>. **/
typedef RowVector_<float>       fRowVector;
/** Abbreviation for RowVector_<double>. **/
typedef RowVector_<double>      dRowVector;
/** Abbreviation for RowVector_<std::complex<float>>. **/
typedef RowVector_<fComplex>    fComplexRowVector;
/** Abbreviation for RowVector_<std::complex<double>>. **/
typedef RowVector_<dComplex>    dComplexRowVector;

/** Abbreviation for RowVectorView_<float>. **/
typedef RowVectorView_<float>   fRowVectorView;
/** Abbreviation for RowVectorView_<double>. **/
typedef RowVectorView_<double>  dRowVectorView;
/** Abbreviation for RowVectorView_<std::complex<float>>. **/
typedef RowVectorView_<fComplex> fComplexRowVectorView;
/** Abbreviation for RowVectorView_<std::complex<double>>. **/
typedef RowVectorView_<dComplex> dComplexRowVectorView;

/** Abbreviation for Matrix_<float>. **/
typedef Matrix_<float>          fMatrix;
/** Abbreviation for Matrix_<double>. **/
typedef Matrix_<double>         dMatrix;
/** Abbreviation for Matrix_<std::complex<float>>. **/
typedef Matrix_<fComplex>       fComplexMatrix;
/** Abbreviation for Matrix_<std::complex<double>>. **/
typedef Matrix_<dComplex>       dComplexMatrix;

/** Abbreviation for MatrixView_<float>. **/
typedef MatrixView_<float>      fMatrixView;
/** Abbreviation for MatrixView_<double>. **/
typedef MatrixView_<double>     dMatrixView;
/** Abbreviation for MatrixView_<std::complex<float>>. **/
typedef MatrixView_<fComplex>   fComplexMatrixView;
/** Abbreviation for MatrixView_<std::complex<double>>. **/
typedef MatrixView_<dComplex>   dComplexMatrixView;
/**@}**/

} //namespace SimTK

#endif //SimTK_SIMMATRIX_BIGMATRIX_H_
