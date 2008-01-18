
/* Portions copyright (c) 2007 Stanford University and Jack Middleton.
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
 *
 * Factors systems of linear algebra equations.
 */


#include <iostream> 
#include <malloc.h>
#include <math.h>
#include <complex>
#include "SimTKcommon.h"
#include "LapackInterface.h"
#include "LinearAlgebra.h"
#include "LATraits.h"
#include "FactorRep.h"
#include "WorkSpace.h"
#include "LapackConvert.h"


namespace SimTK {

   //////////////////////
   // FactorLUDefault  //
   //////////////////////
FactorLUDefault::FactorLUDefault() {
    isFactored = false;
}
FactorLURepBase* FactorLUDefault::clone() const {
    return( new FactorLUDefault(*this));
}


   ///////////////
   // FactorLU //
   ///////////////
FactorLU::~FactorLU() {
    delete rep;
}
// default constructor
FactorLU::FactorLU() {
    rep = new FactorLUDefault();
}

// copy constructor
FactorLU::FactorLU( const FactorLU& c ) {
    rep = c.rep->clone();
}
// copy assignment operator
FactorLU& FactorLU::operator=(const FactorLU& rhs) {
    rep = rhs.rep->clone();
    return *this;
}

template <typename ELT>
void FactorLU::inverse( Matrix_<ELT>& inverse ) const {
    rep->inverse( inverse );
}


template < class ELT >
FactorLU::FactorLU( const Matrix_<ELT>& m ) {
    rep = new FactorLURep<typename CNT<ELT>::StdNumber>(m);
}

template < class ELT >
void FactorLU::factor( const Matrix_<ELT>& m ) {
    delete rep;
    rep = new FactorLURep<typename CNT<ELT>::StdNumber>(m);
}

template < typename ELT >
void FactorLU::solve( const Vector_<ELT>& b, Vector_<ELT>& x ) const {
    rep->solve( b, x );
    return;
}
template < class ELT >
void FactorLU::solve(  const Matrix_<ELT>& b, Matrix_<ELT>& x ) const {
    rep->solve(  b, x );
    return;
}

template < class ELT >
void FactorLU::getL( Matrix_<ELT>& m) const {
    rep->getL( m );
    return;
}
template < class ELT >
void FactorLU::getU( Matrix_<ELT>& m) const {
    rep->getU( m );
    return;
}
template < class ELT >
void FactorLU::getD( Matrix_<ELT>& m) const {
    rep->getD( m );
    return;
}

/*    implement in future release ?
Real FactorLU::getConditionNumber() const  {
   return( rep->getConditionNumber() );
}
template < class ELT >
void FactorLU::getErrorBounds (Vector_<ELT>& err, Vector_<ELT>& berr) const {
    rep->getErrorBounds( err, berr );
    return;
}
*/

bool FactorLU::isSingular () const {
    return( rep->isSingular() );
}

int FactorLU::getSingularIndex () const {
    return( rep->getSingularIndex() );
}
   /////////////////
   // FactorLURep //
   /////////////////
template <typename T >
    template < typename ELT >
FactorLURep<T>::FactorLURep( const Matrix_<ELT>& mat ) 
      : nRow( mat.nrow() ),
        nCol( mat.ncol() ),
        mn( (mat.nrow() < mat.ncol()) ? mat.nrow() : mat.ncol() ),
        positiveDefinite( false ),
        lu( mat.nrow()*mat.ncol() ),
        structure(mat.getMatrixStructure() ),
        condition(mat.getMatrixCondition() ),
        shape(mat.getMatrixShape() ),
        sparsity(mat.getMatrixSparsity() ),
        storage(mat.getMatrixStorage() ),
        pivots(mat.ncol())              { 

	FactorLURep<T>::factor( mat );
}
template <typename T >
FactorLURep<T>::FactorLURep() 
      : nRow(0),
        nCol(0),
        mn(0),
        lu(0),
        pivots(0)             { 
        
}
template <typename T >
FactorLURep<T>::~FactorLURep() {}
 
template <typename T >
FactorLURepBase* FactorLURep<T>::clone() const {
   return( new FactorLURep<T>(*this) );
}

template < class T >
void FactorLURep<T>::inverse(  Matrix_<T>& inverse ) const {
    Matrix_<T> iden(mn,mn);
    iden.resize(mn,mn);
    iden = 1.0;
    solve( iden, inverse );
}

template < class T >
void FactorLURep<T>::solve( const Vector_<T>& b, Vector_<T> &x ) const {

    SimTK_APIARGCHECK2_ALWAYS(b.size()==nRow,"FactorLU","solve",
       "number of rows in right hand side=%d does not match number of rows in original matrix=%d \n",
        b.size(), nRow );

    x.copyAssign(b);
/*  TODO after 1.0
    if( structure == MatrixStructures::Symmetric ) {
        if( condition == MatrixConditions::PositiveDefinite ) {
            LapackInterface::potrs<T>( 'L', nCol, 1, lu.data,  &x(0) );
        } else {
            LapackInterface::sytrs<T>( 'L', nCol, 1, lu.data, pivots.data, &x(0) );
        }
    } else {
*/
        LapackInterface::getrs<T>( 'N', nCol, 1, lu.data, pivots.data, &x(0) );
//    }

    return;
}
template <typename T >
void FactorLURep<T>::solve(  const Matrix_<T>& b, Matrix_<T>& x ) const {

    SimTK_APIARGCHECK2_ALWAYS(b.nrow()==nRow,"FactorLU","solve",
       "number of rows in right hand side=%d does not match number of rows in original matrix=%d \n",
        b.nrow(), nRow );

    x.copyAssign(b);

/* TODO after 1.0
    if( structure == MatrixStructures::Symmetric ) {
        const char uplo = 'L';
        if( condition == MatrixConditions::PositiveDefinite ) {
            LapackInterface::potrs<T>( uplo, nCol, b.ncol(), lu.data,  &x(0,0) );
        } else {
            LapackInterface::sytrs<T>( uplo, nCol, b.ncol(), lu.data, pivots.data, &x(0,0) );
        }
    } else {
*/
        const char trans = 'N';
        LapackInterface::getrs<T>( trans, nCol, b.ncol(), lu.data, pivots.data, &x(0,0) );
//    }

    return;
}
template <typename T >
void FactorLURep<T>::getL( Matrix_<T>& m) const {
       int i,j;
      
       m.resize( nRow, nCol ); 

       for(i=0;i<nRow;i++) {
           for(j=0;j<i;j++) m(j,i) = 0.0;
           for(j=0;j<nCol;j++) m(j,i) = lu.data[j*nRow+i];
       }

    return;
}
template <typename T >
void FactorLURep<T>::getU( Matrix_<T>& m) const {
    int i,j;
    m.resize( nRow, nCol );
       
   for(i = 0;i<nRow;i++) {
       for(j=0;j<i+1;j++) m(j,i) = lu.data[j*nRow+i];
       for(;j<nCol;j++) m(j,i) = 0.0;
   }
   return;
}
template <typename T >
void FactorLURep<T>::getD( Matrix_<T>& m) const {
   int i,j;


   m.resize( nRow, nCol );
   for(i = 0;i<nRow;i++) {
        for(j=0;j<nCol;j++) m(j,i) = 0.0;
   }
   for(i = 0;i<nRow;i++) m(i,i) = lu.data[i*nRow+i];
       
   return;
}
template <typename T >
Real FactorLURep<T>::getConditionNumber() const {
   return 0.0;
}
template <typename T >
void FactorLURep<T>::getErrorBounds ( Vector_<T>& err, Vector_<T>& berr) const {
    return;
}
template <typename T >
bool FactorLURep<T>::isSingular () const {

    if( singularIndex) 
        return( true );
    else 
        return( false );
}
template <typename T >
int FactorLURep<T>::getSingularIndex () const {
    return( singularIndex );
}

template <class T> 
    template<typename ELT>
void FactorLURep<T>::factor(const Matrix_<ELT>&mat )  {

    elementSize = sizeof( T );
    imagOffset = CNT<ELT>::ImagOffset;  // real/complex (usefull for debugging)
   
    // initialize the matrix we pass to LAPACK
    // converts (negated,conjugated etc.) to LAPACK format 
    LapackConvert::convertMatrixToLapack( lu.data, mat );


    int lda = nRow;
    int info;

/*  TODO after 1.0
    if( structure == MatrixStructures::Symmetric ) {
        if( condition == MatrixConditions::PositiveDefinite ) {
            positiveDefinite = true; 
            if( storage == MatrixStorageFormats::Packed ) {
//                LapackInterface::pptrf<ELT>( );     
            } else if( sparsity == MatrixSparseFormats::Banded ) {
//  TODO               LapackInterface::pbtrf<ELT>( );     
            } else if( structure == MatrixStructures::TriDiagonal ) {
//    TODO             LapackInterface::pttrf<ELT>( );     
            } else {
                char uplo = 'L';
                LapackInterface::potrf<T>(uplo, nRow, lu.data, lda, info);     
                if( info > 0 ) {
                    singularIndex = info;
                    SimTK_THROW2( SimTK::Exception::NotPositiveDefinite, 
                    getSingularIndex(), "FactorLU:LapackInterface::potrf" ); 
                }
 //           }
        }  else {

            if( storage == MatrixStorageFormats::Packed ) {
//                LapackInterface::sptrf<ELT>();
            } else {

                long workSize = nCol*LapackInterface::ilaenv<T>(1, "sytrf", "U", nCol, -1, -1, -1);
                TypedWorkSpace<T>  work( workSize );
                
                LapackInterface::sytrf<T>(nRow, nCol, lu.data, lda, pivots.data, work.data, workSize, info);
                if( info > 0 ) {
                    singularIndex = info;
                    SimTK_THROW2( SimTK::Exception::SingularMatrix, 
                    getSingularIndex(), "FactorLU:LapackInterface::sytrf" ); 
                }
 //           }
        }
    } else {

        if( sparsity == MatrixSparseFormats::Banded ) {
//    TODO          LapackInterface::gbtrf<T>(nRow, nCol kl, ku, lu.data, lda, pivots.data, info);
        } else if( structure == MatrixStructures::Tridiagonal ) {
//             double *dl, *d, *du, *du2;
//    TODO          LapackInterface::gttrf<T>(nRow, nCol, dl, d, du, du2, pivots.data, info);
            if( info > 0 ) {
                singularIndex = info;
                SimTK_THROW2( SimTK::Exception::SingularMatrix, 
                getSingularIndex(), "FactorLU:LapackInterface::gttrf" ); 
            }
        } else {
*/
            LapackInterface::getrf<T>(nRow, nCol, lu.data, lda, pivots.data, info);
            if( info > 0 ) {
                singularIndex = info;
                SimTK_THROW2( SimTK::Exception::SingularMatrix, 
                getSingularIndex(), "FactorLU:LapackInterface::getrf" ); 
            }
 //       }
 //   }

}

// instantiate
template SimTK_SIMMATH_EXPORT FactorLU::FactorLU( const Matrix_<double>& m );
template SimTK_SIMMATH_EXPORT FactorLU::FactorLU( const Matrix_<float>& m );
template SimTK_SIMMATH_EXPORT FactorLU::FactorLU( const Matrix_<std::complex<float> >& m );
template SimTK_SIMMATH_EXPORT FactorLU::FactorLU( const Matrix_<std::complex<double> >& m );
template SimTK_SIMMATH_EXPORT FactorLU::FactorLU( const Matrix_<conjugate<float> >& m );
template SimTK_SIMMATH_EXPORT FactorLU::FactorLU( const Matrix_<conjugate<double> >& m );
template SimTK_SIMMATH_EXPORT FactorLU::FactorLU( const Matrix_<negator< double> >& m );
template SimTK_SIMMATH_EXPORT FactorLU::FactorLU( const Matrix_<negator< float> >& m );
template SimTK_SIMMATH_EXPORT FactorLU::FactorLU( const Matrix_<negator< std::complex<float> > >& m );
template SimTK_SIMMATH_EXPORT FactorLU::FactorLU( const Matrix_<negator< std::complex<double> > >& m );
template SimTK_SIMMATH_EXPORT FactorLU::FactorLU( const Matrix_<negator< conjugate<float> > >& m );
template SimTK_SIMMATH_EXPORT FactorLU::FactorLU( const Matrix_<negator< conjugate<double> > >& m );

template SimTK_SIMMATH_EXPORT void FactorLU::factor( const Matrix_<double>& m );
template SimTK_SIMMATH_EXPORT void FactorLU::factor( const Matrix_<float>& m );
template SimTK_SIMMATH_EXPORT void FactorLU::factor( const Matrix_<std::complex<float> >& m );
template SimTK_SIMMATH_EXPORT void FactorLU::factor( const Matrix_<std::complex<double> >& m );
template SimTK_SIMMATH_EXPORT void FactorLU::factor( const Matrix_<conjugate<float> >& m );
template SimTK_SIMMATH_EXPORT void FactorLU::factor( const Matrix_<conjugate<double> >& m );
template SimTK_SIMMATH_EXPORT void FactorLU::factor( const Matrix_<negator< double> >& m );
template SimTK_SIMMATH_EXPORT void FactorLU::factor( const Matrix_<negator< float> >& m );
template SimTK_SIMMATH_EXPORT void FactorLU::factor( const Matrix_<negator< std::complex<float> > >& m );
template SimTK_SIMMATH_EXPORT void FactorLU::factor( const Matrix_<negator< std::complex<double> > >& m );
template SimTK_SIMMATH_EXPORT void FactorLU::factor( const Matrix_<negator< conjugate<float> > >& m );
template SimTK_SIMMATH_EXPORT void FactorLU::factor( const Matrix_<negator< conjugate<double> > >& m );

template class FactorLURep<double>;
template FactorLURep<double>::FactorLURep( const Matrix_<double>& m);
template FactorLURep<double>::FactorLURep( const Matrix_<negator<double> >& m);
template void FactorLURep<double>::factor( const Matrix_<double>& m);
template void FactorLURep<double>::factor( const Matrix_<negator<double> >& m);

template class FactorLURep<float>;
template FactorLURep<float>::FactorLURep( const Matrix_<float>& m);
template FactorLURep<float>::FactorLURep( const Matrix_<negator<float> >& m);
template void FactorLURep<float>::factor( const Matrix_<float>& m);
template void FactorLURep<float>::factor( const Matrix_<negator<float> >& m);

template class FactorLURep<std::complex<double> >;
template FactorLURep<std::complex<double> >::FactorLURep( const Matrix_<std::complex<double> >& m);
template FactorLURep<std::complex<double> >::FactorLURep( const Matrix_<negator<std::complex<double> > >& m);
template FactorLURep<std::complex<double> >::FactorLURep( const Matrix_<conjugate<double> >& m);
template FactorLURep<std::complex<double> >::FactorLURep( const Matrix_<negator<conjugate<double> > >& m);
template void FactorLURep<std::complex<double> >::factor( const Matrix_<std::complex<double> >& m);
template void FactorLURep<std::complex<double> >::factor( const Matrix_<negator<std::complex<double> > >& m);
template void FactorLURep<std::complex<double> >::factor( const Matrix_<conjugate<double> >& m);
template void FactorLURep<std::complex<double> >::factor( const Matrix_<negator<conjugate<double> > >& m);

template class FactorLURep<std::complex<float> >;
template FactorLURep<std::complex<float> >::FactorLURep( const Matrix_<std::complex<float> >& m);
template FactorLURep<std::complex<float> >::FactorLURep( const Matrix_<negator<std::complex<float> > >& m);
template FactorLURep<std::complex<float> >::FactorLURep( const Matrix_<conjugate<float> >& m);
template FactorLURep<std::complex<float> >::FactorLURep( const Matrix_<negator<conjugate<float> > >& m);
template void FactorLURep<std::complex<float> >::factor( const Matrix_<std::complex<float> >& m);
template void FactorLURep<std::complex<float> >::factor( const Matrix_<negator<std::complex<float> > >& m);
template void FactorLURep<std::complex<float> >::factor( const Matrix_<conjugate<float> >& m);
template void FactorLURep<std::complex<float> >::factor( const Matrix_<negator<conjugate<float> > >& m);

template SimTK_SIMMATH_EXPORT void FactorLU::getL<float>(Matrix_<float>&) const;
template SimTK_SIMMATH_EXPORT void FactorLU::getL<double>(Matrix_<double>&) const;
template SimTK_SIMMATH_EXPORT void FactorLU::getL<std::complex<float> >(Matrix_<std::complex<float> >&) const;
template SimTK_SIMMATH_EXPORT void FactorLU::getL<std::complex<double> >(Matrix_<std::complex<double> >&) const;
template SimTK_SIMMATH_EXPORT void FactorLU::getD<float>(Matrix_<float>&) const;
template SimTK_SIMMATH_EXPORT void FactorLU::getD<double>(Matrix_<double>&) const;
template SimTK_SIMMATH_EXPORT void FactorLU::getD<std::complex<float> >(Matrix_<std::complex<float> >&) const;
template SimTK_SIMMATH_EXPORT void FactorLU::getD<std::complex<double> >(Matrix_<std::complex<double> >&) const;
template SimTK_SIMMATH_EXPORT void FactorLU::getU<float>(Matrix_<float>&) const;
template SimTK_SIMMATH_EXPORT void FactorLU::getU<double>(Matrix_<double>&) const;
template SimTK_SIMMATH_EXPORT void FactorLU::getU<std::complex<float> >(Matrix_<std::complex<float> >&) const;
template SimTK_SIMMATH_EXPORT void FactorLU::getU<std::complex<double> >(Matrix_<std::complex<double> >&) const;
template SimTK_SIMMATH_EXPORT void FactorLU::solve<float>(const Vector_<float>&, Vector_<float>&) const;
template SimTK_SIMMATH_EXPORT void FactorLU::solve<double>(const Vector_<double>&, Vector_<double>&) const;
template SimTK_SIMMATH_EXPORT void FactorLU::solve<std::complex<float> >(const Vector_<std::complex<float> >&, Vector_<std::complex<float> >&) const;
template SimTK_SIMMATH_EXPORT void FactorLU::solve<std::complex<double> >(const Vector_<std::complex<double> >&, Vector_<std::complex<double> >&) const;
template SimTK_SIMMATH_EXPORT void FactorLU::solve<float>(const Matrix_<float>&, Matrix_<float>&) const;
template SimTK_SIMMATH_EXPORT void FactorLU::solve<double>(const Matrix_<double>&, Matrix_<double>&) const;
template SimTK_SIMMATH_EXPORT void FactorLU::solve<std::complex<float> >(const Matrix_<std::complex<float> >&, Matrix_<std::complex<float> >&) const;
template SimTK_SIMMATH_EXPORT void FactorLU::solve<std::complex<double> >(const Matrix_<std::complex<double> >&, Matrix_<std::complex<double> >&) const;
template SimTK_SIMMATH_EXPORT void FactorLU::inverse<float>(Matrix_<float>&) const;
template SimTK_SIMMATH_EXPORT void FactorLU::inverse<double>(Matrix_<double>&) const;
template SimTK_SIMMATH_EXPORT void FactorLU::inverse<std::complex<float> >(Matrix_<std::complex<float> >&) const;
template SimTK_SIMMATH_EXPORT void FactorLU::inverse<std::complex<double> >(Matrix_<std::complex<double> >&) const;


} // namespace SimTK
