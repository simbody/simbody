
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
#include "FactorRep.h"
#include "WorkSpace.h"
#include "LATraits.h"


namespace SimTK {

   ///////////////
   // FactorLU //
   ///////////////
FactorLU::~FactorLU() {
    delete rep;
}




template < class ELT >
FactorLU::FactorLU( const Matrix_<ELT>& m ) {
    rep = new FactorLURep<typename CNT<ELT>::StdNumber>(m);
}


template < typename ELT >
void FactorLU::solve( const Vector_<ELT>& b, Vector_<ELT>& x ) {
    rep->solve( b, x );
    return;
}
template < class ELT >
void FactorLU::solve(  const Matrix_<ELT>& b, Matrix_<ELT>& x ) {
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
Real FactorLU::getConditionNumber() const  {
   return( rep->getConditionNumber() );
}
template < class ELT >
void FactorLU::getErrorBounds (Vector_<ELT>& err, Vector_<ELT>& berr) const {
    rep->getErrorBounds( err, berr );
    return;
}
bool FactorLU::isSingular () const {
    return( rep->isSingular() );
}

int FactorLU::getSingularIndex () const {
    return( rep->getSingularIndex() );
}
void FactorLU::display (int opts) {
    return( rep->display(opts) );
}
   /////////////////
   // FactorLURep //
   /////////////////
template <typename T >
    template < typename ELT >
FactorLURep<T>::FactorLURep( const Matrix_<ELT>& mat ) 
      : nRow( mat.nrow() ),
        nCol( mat.ncol() ),
        lu( mat.nrow()*mat.ncol() ),
        pivots(mat.ncol())             { 
        
	FactorLURep<T>::factor( mat );
}

template <typename T >
void FactorLURep<T>::display(int options) {
   if( options > 1 ) {
       printf("LU = \n");
       for(int i = 0;i<nRow;i++) {
           for(int j = 0;j<nCol;j++) printElement( i,j);
           printf("\n");
       } 

   }
   if( options > 2 ) {
       printf("pivots = ");
       for(int j = 0;j<nCol;j++) printf("%d ",pivots.data[j]);
       printf("\n");
   }
   return;
}
template <typename T >
void FactorLURep<T>::printElement( int i, int j) {

    if( imagOffset == 1 ) {
        if( elementSize == 8 ) {
//  TODO need to specialize for complex types           printf(" %f,%f", lu.data[j*nRow+i].real(), lu.data[j*nRow+i].imag() ); 
        } else { 
//            printf(" %f,%f", lu.data[j*nRow+i].real(), lu.data[j*nRow+i].imag() ); 
        }
    } else {
        if( elementSize == 4 ) {
//            printf(" %f", lu.data[j*nRow+i] ); 
        } else {
//            printf(" %f", lu.data[j*nRow+i] ); 
        }
    } 
    
}
template <typename T >
FactorLURep<T>::~FactorLURep() {}

template < class T >
    template < class ELT >
void FactorLURep<T>::solve( const Vector_<ELT>& b, Vector_<ELT> &x ) {
    x.copyAssign(b);
// TODO check  that ELT of b is same as the factored matrix (size,imageoffset)
    LapackInterface::getrs<T>( false, nCol, 1, lu.data, pivots.data, &x(0));
    return;
}
template < class T >
void FactorLURep<T>::copyElement( const int i, const int j, T* ptr ) {
    *ptr =  lu.data[j*nRow+i];
}
template < class T >
    template < class ELT >
void FactorLURep<T>::copyElement( const int i, const int j, std::complex<ELT> *ptr ) {
    *ptr =  std::complex<ELT>(lu.data[j*nRow+i].real(), lu.data[j*nRow+i].imag() ); 
}
// TODO handle cases where length to b,x and dimensions of lu are not consistant
template <typename T >
    template < class ELT >
void FactorLURep<T>::solve(  const Matrix_<ELT>& b, Matrix_<ELT>& x ) {
    x.copyAssign(b);
    LapackInterface::getrs<ELT>( false, nCol, b.ncol(), lu.data, pivots.data, &x(0,0));
    return;
}
template <typename T >
   template < class ELT >
void FactorLURep<T>::getL( Matrix_<ELT>& m) const {
       int i,j;
      
       m.resize( nRow, nCol ); 

       for(i = 0;i<nRow;i++) {
           for(j = 0;j<i;j++) m(j,i) = 0.0;
           for(;j<nCol;j++) copyElement( j, i, &m(j,i) ) ;
       }

    return;
}
template <typename T >
    template < class ELT >
void FactorLURep<T>::getU( Matrix_<ELT>& m) const {
       int i,j;
       m.resize( nRow, nCol );
       
       for(i = 0;i<nRow;i++) {
           for(j=0;j<i+1;j++) copyElement( j,i, &m(j,i) );
           for(;j<nCol;j++) m(j,i) = 0.0;
       }
    return;
}
template <typename T >
    template < class ELT >
void FactorLURep<T>::getD( Matrix_<ELT>& m) const {
       int i,j;
       m.resize( nRow, nCol );
       for(i = 0;i<nRow;i++) {
           for(j=0;j<nCol;j++) m(j,i) = 0.0;
       }
       for(i = 0;i<nRow;i++) copyElement( i, i, &m(i,i) );
       
    return;
}
template <typename T >
Real FactorLURep<T>::getConditionNumber() const {
   return 0.0;
}
template <typename T >
     template < class ELT >
void FactorLURep<T>::getErrorBounds ( Vector_<ELT>& err, Vector_<ELT>& berr) const {
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
template <typename T >
    template < class ELT> 
void FactorLURep<T>::initLU( const Matrix_<ELT>& mat ) {
    for(int i=0;i<nCol;i++) {
        for(int j=0;j<nRow;j++) lu.data[i*nRow+j] = mat(j,i);
    }
    
    return;
}
template <typename T >
    template < class ELT > 
void FactorLURep<T>::initLU( const Matrix_<negator<ELT> >& mat ) {
    for(int i=0;i<nCol;i++) {
        for(int j=0;j<nRow;j++) {
            lu.data[i*nRow+j] = -mat(j,i);
        }
    }

    return;
}
template < class T > 
    template < class ELT > 
void FactorLURep<T>::initLU( const Matrix_<std::complex<ELT> >& mat ) {
    for(int i=0;i<nCol;i++) {
        for(int j=0;j<nRow;j++) {
            lu.data[i*nRow+j] = std::complex<ELT>(mat(j,i).real(),mat(j,i).imag());
        }
    }

    return;
}
template < class T > 
    template < class ELT > 
void FactorLURep<T>::initLU( const Matrix_<negator<std::complex<ELT> > >& mat ) {
    for(int i=0;i<nCol;i++) {
        for(int j=0;j<nRow;j++) {
            lu.data[i*nRow+j] = std::complex<ELT>(-mat(j,i).real(),-mat(j,i).imag());
        }
    }

    return;
}
template < class T > 
    template < class ELT > 
void FactorLURep<T>::initLU( const Matrix_<conjugate<ELT> >& mat ) {
    for(int i=0;i<nCol;i++) {
        for(int j=0;j<nRow;j++) {
            lu.data[i*nRow+j] = std::complex<ELT>(mat(j,i).real(),-mat(j,i).imag());
        }
    }

    return;
}
template < class T > 
    template < class ELT > 
void FactorLURep<T>::initLU( const Matrix_<negator<conjugate<ELT> > >& mat ) {
    for(int i=0;i<nCol;i++) {
        for(int j=0;j<nRow;j++) {
            lu.data[i*nRow+j] = std::complex<ELT>(-mat(j,i).real(),mat(j,i).imag());
        }
    }

    return;
}

template <class T> 
    template<typename ELT>
void FactorLURep<T>::factor(const Matrix_<ELT>&mat )  {

    elementSize = sizeof( T );
    imagOffset = CNT<ELT>::ImagOffset;  // real/complex (usefull for debugging)
   
    // allocate and initialize the matrix we pass to LAPACK
    // convert (negated,conjugated etc.) to LAPACK format 
    initLU( mat ); 

    int lda = nRow;
    int info;

    MatrixStructures::Structure structure  = mat.getMatrixStructure();
    MatrixShapes::Shape shape              = mat.getMatrixShape();
    MatrixSparseFormats::Sparsity sparsity = mat.getMatrixSparsity();
    MatrixStorageFormats::Storage storage  = mat.getMatrixStorage();
    MatrixConditions::Condition condition  = mat.getMatrixCondition();

    if( structure == MatrixStructures::Symmetric ) {
        if( condition == MatrixConditions::PositiveDefinite ) {
            if( storage == MatrixStorageFormats::Packed ) {
//                LapackInterface::pptrf<ELT>( );     
            } else if( sparsity == MatrixSparseFormats::Banded ) {
//                LapackInterface::pbtrf<ELT>( );     
            } else if( structure == MatrixStructures::TriDiagonal ) {
//                LapackInterface::pttrf<ELT>( );     
            } else {
                int kl = nRow;  // TODO poperly set these
                int ku = nCol;
// TODO check  that ELT of b is same as the factored matrix (size,imageoffset)
                LapackInterface::potrf<T>(nRow,nCol,kl,ku, lu.data, lda, pivots.data, info);     
            }
        }  else {
            if( storage == MatrixStorageFormats::Packed ) {
//                LapackInterface::sptrf<ELT>();
            } else {
                long workSize = nCol*LapackInterface::ilaenv<T>(1, "sytrf", "U", nCol, -1, -1, -1);
                WorkSpace  work( workSize*sizeof(T) );
                
                LapackInterface::sytrf<T>(nRow, nCol, lu.data, lda, pivots.data, work.getData<T>(), workSize, info);
            }
        }
    } else {
        if( sparsity == MatrixSparseFormats::Banded ) {
//             LapackInterface::gbtrf<T>(nRow, nCol kl, ku, lu.data, lda, pivots.data, info);
        } else if( structure == MatrixStructures::Triangular ) {
             double *dl, *d, *du, *du2;
//             LapackInterface::gttrf<T>(nRow, nCol, dl, d, du, du2, pivots.data, info);
        } else {
             LapackInterface::getrf<T>(nRow, nCol, lu.data, lda, pivots.data, info);
        }
    }

    if( info < 0 ) {
       // TODO arg #info bad value throw
    } 
}

// instantiate
template FactorLU::FactorLU( const Matrix_<double>& m );
template FactorLU::FactorLU( const Matrix_<float>& m );
template FactorLU::FactorLU( const Matrix_<std::complex<float> >& m );
template FactorLU::FactorLU( const Matrix_<std::complex<double> >& m );
template FactorLU::FactorLU( const Matrix_<conjugate<float> >& m );
template FactorLU::FactorLU( const Matrix_<conjugate<double> >& m );
template FactorLU::FactorLU( const Matrix_<negator< double> >& m );
template FactorLU::FactorLU( const Matrix_<negator< float> >& m );
template FactorLU::FactorLU( const Matrix_<negator< std::complex<float> > >& m );
template FactorLU::FactorLU( const Matrix_<negator< std::complex<double> > >& m );
template FactorLU::FactorLU( const Matrix_<negator< conjugate<float> > >& m );
template FactorLU::FactorLU( const Matrix_<negator< conjugate<double> > >& m );

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

template void FactorLU::getL<float>(Matrix_<float>&) const;
template void FactorLU::getL<double>(Matrix_<double>&) const;
template void FactorLU::getL<std::complex<float> >(Matrix_<std::complex<float> >&) const;
template void FactorLU::getL<std::complex<double> >(Matrix_<std::complex<double> >&) const;
template void FactorLU::getD<float>(Matrix_<float>&) const;
template void FactorLU::getD<double>(Matrix_<double>&) const;
template void FactorLU::getD<std::complex<float> >(Matrix_<std::complex<float> >&) const;
template void FactorLU::getD<std::complex<double> >(Matrix_<std::complex<double> >&) const;
template void FactorLU::getU<float>(Matrix_<float>&) const;
template void FactorLU::getU<double>(Matrix_<double>&) const;
template void FactorLU::getU<std::complex<float> >(Matrix_<std::complex<float> >&) const;
template void FactorLU::getU<std::complex<double> >(Matrix_<std::complex<double> >&) const;
template void FactorLU::solve<float>(const Vector_<float>&, Vector_<float>&);
template void FactorLU::solve<double>(const Vector_<double>&, Vector_<double>&);
template void FactorLU::solve<std::complex<float> >(const Vector_<std::complex<float> >&, Vector_<std::complex<float> >&);
template void FactorLU::solve<std::complex<double> >(const Vector_<std::complex<double> >&, Vector_<std::complex<double> >&);
template void FactorLU::solve<float>(const Matrix_<float>&, Matrix_<float>&);
template void FactorLU::solve<double>(const Matrix_<double>&, Matrix_<double>&);
template void FactorLU::solve<std::complex<float> >(const Matrix_<std::complex<float> >&, Matrix_<std::complex<float> >&);
template void FactorLU::solve<std::complex<double> >(const Matrix_<std::complex<double> >&, Matrix_<std::complex<double> >&);

} // namespace SimTK
