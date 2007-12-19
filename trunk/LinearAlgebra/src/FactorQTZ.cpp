
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
#include "FactorQTZRep.h"
#include "WorkSpace.h"
#include "LATraits.h"
#include "LapackConvert.h"


namespace SimTK {

static const double ZERO = 0.0;
static const double ONE  = 1.0;
   ///////////////
   // FactorQTZ //
   ///////////////
FactorQTZ::~FactorQTZ() {
    delete rep;
}

template < class ELT >
FactorQTZ::FactorQTZ( const Matrix_<ELT>& m, typename CNT<ELT>::TReal rcond ) {
    rep = new FactorQTZRep<typename CNT<ELT>::StdNumber>(m, rcond); 
}


template < typename ELT >
void FactorQTZ::solve( const Vector_<ELT>& b, Vector_<ELT>& x ) {
    rep->solve( b, x );
    return;
}
template < class ELT >
void FactorQTZ::solve(  const Matrix_<ELT>& b, Matrix_<ELT>& x ) {
    rep->solve(  b, x );
    return;
}

   /////////////////
   // FactorQTZRep //
   /////////////////
template <typename T >
    template < typename ELT >
FactorQTZRep<T>::FactorQTZRep( const Matrix_<ELT>& mat, typename CNT<T>::TReal rc) 
      : nRow( mat.nrow() ),
        nCol( mat.ncol() ),
        rank(0),
        qtz( mat.nrow()*mat.ncol() ),
        pivots(mat.ncol()),
        mn( (mat.nrow() < mat.ncol()) ? mat.nrow() : mat.ncol() ),
        tauGEQP3(mn),
        tauORMQR(mn),
        scaleLinSys(false),
        scaleRHS(false)             { 
        
        rcond = rc;
        for(int i=0;i<mat.ncol();i++) pivots.data[i] = 0;
	FactorQTZRep<T>::factor( mat );
}

template <typename T >
FactorQTZRep<T>::~FactorQTZRep() {}

template < class T >
void FactorQTZRep<T>::solve( const Vector_<T>& b, Vector_<T> &x ) {
    Matrix_<T> m(b.nrow(), 1 );
    m.copyAssign(b);
    Matrix_<T> r(nCol, 1 );
    doSolve( m, r );
    x.copyAssign(r);
    return;
}
// TODO handle cases where length of b,x and dimensions of  are not consistant
template <typename T >
void FactorQTZRep<T>::solve(  const Matrix_<T>& b, Matrix_<T>& x ) {
    x.resize(nCol, b.ncol() );
    Matrix_<T> tb;
    tb.copyAssign(b);
    doSolve(tb, x);
}
template <typename T >
void FactorQTZRep<T>::doSolve(  Matrix_<T>& b, Matrix_<T>& x) {
    int i,j;
    int info;
    typedef typename CNT<T>::TReal RealType;
    RealType bnrm, smlnum, bignum;
    int nrhs = b.ncol();
    int n = nCol;
    int m = nRow;
    // compute size of workspace 
    // for dormqr, dormrz:  lwork = n*nb
    long lwork1 = n*LapackInterface::ilaenv<T>(1, "ormqr", "LT ", nRow, b.ncol(), -1, -1);
    long lwork2 = n*LapackInterface::ilaenv<T>(1, "ormrz", "LUNN", rank, b.ncol(), -1, -1);
    TypedWorkSpace<T> work( lwork1>lwork2 ? lwork1 : lwork2);

    // compute norm of RHS
    bnrm = (RealType)LapackInterface::lange<T>( 'M', m, nrhs, &b(0,0), b.nrow() );

    LapackInterface::getMachinePrecision<RealType>( smlnum, bignum);
 
    // compute scale for RHS
    if( bnrm > ZERO  && bnrm < smlnum ) {
        scaleRHS = true;
        rhsScaleF = smlnum;
    } else if( bnrm > bignum ) {
        scaleRHS = true;
        rhsScaleF = bignum;
    }

    if( scaleRHS ) {  // apply scale factor to RHS
        LapackInterface::lascl<T>( 'G', 0, 0, bnrm, rhsScaleF, m, nrhs, &b(0,0), b.nrow(), info ); 
    }
    // 
    LapackInterface::ormqr<T>( 'L', 'T', nRow, b.ncol(), mn, qtz.data, nRow, tauGEQP3.data,
                               &b(0,0), b.nrow(), work.data, work.size, info );
    LapackInterface::trsm<T>( 'L', 'U', 'N', 'N', rank, b.ncol(), 1.0, qtz.data, nRow, &b(0,0), b.nrow() );

    //  zero out elements of RHS for rank deficient systems
    for( j = 0; j<nrhs; j++ ) {
        for( i = rank; i<n; i++ ) {
            b(i,j) = 0;
        }
    }
   
    if( rank < nCol ) {
// TODO  fix SimTKlapack.h  LapackInterface::ormrz<T>('L', 'T', nCol, x.ncol(), rank, nCol-rank, qtz.data, nRow, 
        int l = nCol-rank;
        LapackInterface::ormrz<T>('L', 'T', nCol, b.ncol(), rank, &l, qtz.data, nRow, 
                                   tauORMQR.data, &b(0,0), b.nrow(), work.data, work.size, info );
    }

    // adjust for pivoting
    for( j = 0; j<nrhs; j++ ) {
        for( i = 0; i<n; i++ ) {
            work.data[ pivots.data[i]-1] = b(i,j);
        }
        LapackInterface::copy<T>(n, work.data, 1, &x(0,j), 1 );
    }

 
    // compensate for scaling of linear system 
    if( scaleLinSys ) { 
        LapackInterface::lascl<T>( 'g', 0, 0, anrm, linSysScaleF, nCol, x.ncol(), &x(0,0), 1, info );
    }

    // compensate for scaling of RHS 
    if( scaleRHS  ) { 
        LapackInterface::lascl<T>('g', 0, 0, bnrm, rhsScaleF, nCol, x.ncol(), &x(0,0), 1, info);
    }
    
    return;
}

template <class T> 
    template<typename ELT>
void FactorQTZRep<T>::factor(const Matrix_<ELT>&mat )  {

    // allocate and initialize the matrix we pass to LAPACK
    // converts (negated,conjugated etc.) to LAPACK format 
    LapackConvert::convertMatrixToLapack( qtz.data, mat );

    int info;
    const int smallestSingularValue = 2;
    const int largestSingularValue = 1;
    typedef typename CNT<T>::TReal  RealType;
    RealType smlnum, bignum, smin, smax;

    // compute optimal block size
    // dtzrzf: lwork = m*nb  
    long lwork1 = nRow*LapackInterface::ilaenv<T>(1, "tzrzf", " ", nCol, -1, -1, -1);

    // dgepq3: lwork = 2*N+( N+1 )*NB
    long lwork2 = 2*nCol + (nCol+1)*LapackInterface::ilaenv<T>(1, "geqp3", " ", nRow, nCol,  -1, -1);
    TypedWorkSpace<T> work( lwork1>lwork2 ? lwork1 : lwork2 );

    LapackInterface::getMachinePrecision<RealType>( smlnum, bignum);

    // scale the input system of equations
    anrm = (RealType)LapackInterface::lange<T>( 'M', nRow, nCol, qtz.data, nRow );

    if( anrm > 0 && anrm < smlnum ) {
        scaleLinSys = true;
        linSysScaleF = smlnum;
    } else if( anrm > bignum )  {
        scaleLinSys = true;
        linSysScaleF = bignum;
    } 
    if( anrm == 0  ) { // matrix all zeros
        rank = 0;
    } else {
        if ( scaleLinSys ) {
           LapackInterface::lascl<T>( 'G', 0, 0, anrm, linSysScaleF, nRow, nCol, qtz.data, nRow, info );
        }

        // compute QR factorization with column pivoting: A = Q * R
        // Q * R is returned in qtz.data
        LapackInterface::geqp3<T>(nRow, nCol, qtz.data, nRow, pivots.data, 
                                tauGEQP3.data, work.data, work.size, info );

        // compute Rank
        work.data[0] = 1;
        work.data[mn] = 1;
        smax = CNT<T>::abs( qtz.data[0] );
        smin = smax;
        if( CNT<T>::abs(qtz.data[0]) == 0 ) {
            rank = 0;
        } else {
            T s1,s2,c1,c2;
            RealType smaxpr,sminpr;

            // Determine rank using incremental condition estimate
            for( rank=1,smaxpr=0.0,sminpr=1.0; rank<mn && smaxpr*rcond < sminpr; ) {

                LapackInterface::laic1<T>( smallestSingularValue, rank, work.data, smin,
                      &qtz.data[rank*nRow], qtz.data[(rank*nRow)+rank], sminpr, s1, c1 );

                LapackInterface::laic1<T>( largestSingularValue,  rank, &work.data[mn], smax,
                       &qtz.data[rank*nRow], qtz.data[(rank*nRow)+rank], smaxpr, s2, c2 );

                if( smaxpr*rcond < sminpr ) {
                    for(int i=0;i<rank;i++) {
                         work.data[i]    *= s1;
                         work.data[i+mn] *= s2;
                    }
                    work.data[rank] = c1;
                    work.data[rank+mn] = c2;
                    smin = sminpr;
                    smax = smaxpr;
                    rank++;
                }
            } 

            // R => T * Z
            // Matrix is rank deficient so complete the orthogonalization by 
            // applying orthogonal transformations from the right to 
            // factor R from an upper trapezoidal to an upper triangular matrix
            // T is returned in qtz.data and Z is returned in tauORMQR.data
            if( rank < nCol ) {
                LapackInterface::tzrzf<T>( rank, nCol, qtz.data, nRow, tauORMQR.data, work.data, 
                                    work.size, info );
            }
        }
    }

    return;
}

// instantiate
template SimTK_SIMMATH_EXPORT FactorQTZ::FactorQTZ( const Matrix_<double>& m, double rcond );
template SimTK_SIMMATH_EXPORT FactorQTZ::FactorQTZ( const Matrix_<float>& m, float rcond );
template SimTK_SIMMATH_EXPORT FactorQTZ::FactorQTZ( const Matrix_<std::complex<float> >& m, float rcond );
template SimTK_SIMMATH_EXPORT FactorQTZ::FactorQTZ( const Matrix_<std::complex<double> >& m, double rcond );
template SimTK_SIMMATH_EXPORT FactorQTZ::FactorQTZ( const Matrix_<conjugate<float> >& m, float rcond );
template SimTK_SIMMATH_EXPORT FactorQTZ::FactorQTZ( const Matrix_<conjugate<double> >& m, double rcond );
template SimTK_SIMMATH_EXPORT FactorQTZ::FactorQTZ( const Matrix_<negator< double> >& m, double rcond );
template SimTK_SIMMATH_EXPORT FactorQTZ::FactorQTZ( const Matrix_<negator< float> >& m, float rcond );
template SimTK_SIMMATH_EXPORT FactorQTZ::FactorQTZ( const Matrix_<negator< std::complex<float> > >& m, float rcond );
template SimTK_SIMMATH_EXPORT FactorQTZ::FactorQTZ( const Matrix_<negator< std::complex<double> > >& m, double rcond );
template SimTK_SIMMATH_EXPORT FactorQTZ::FactorQTZ( const Matrix_<negator< conjugate<float> > >& m, float rcond );
template SimTK_SIMMATH_EXPORT FactorQTZ::FactorQTZ( const Matrix_<negator< conjugate<double> > >& m, double rcond );

template class FactorQTZRep<double>;
template FactorQTZRep<double>::FactorQTZRep( const Matrix_<double>& m, double rcond);
template FactorQTZRep<double>::FactorQTZRep( const Matrix_<negator<double> >& m, double rcond);
template void FactorQTZRep<double>::factor( const Matrix_<double>& m);
template void FactorQTZRep<double>::factor( const Matrix_<negator<double> >& m);

template class FactorQTZRep<float>;
template FactorQTZRep<float>::FactorQTZRep( const Matrix_<float>& m, float rcond );
template FactorQTZRep<float>::FactorQTZRep( const Matrix_<negator<float> >& m, float rcond );
template void FactorQTZRep<float>::factor( const Matrix_<float>& m);
template void FactorQTZRep<float>::factor( const Matrix_<negator<float> >& m);

template class FactorQTZRep<std::complex<double> >;
template FactorQTZRep<std::complex<double> >::FactorQTZRep( const Matrix_<std::complex<double> >& m, double rcond);
template FactorQTZRep<std::complex<double> >::FactorQTZRep( const Matrix_<negator<std::complex<double> > >& m, double rcond);
template FactorQTZRep<std::complex<double> >::FactorQTZRep( const Matrix_<conjugate<double> >& m, double rcond);
template FactorQTZRep<std::complex<double> >::FactorQTZRep( const Matrix_<negator<conjugate<double> > >& m, double rcond);
template void FactorQTZRep<std::complex<double> >::factor( const Matrix_<std::complex<double> >& m);
template void FactorQTZRep<std::complex<double> >::factor( const Matrix_<negator<std::complex<double> > >& m);
template void FactorQTZRep<std::complex<double> >::factor( const Matrix_<conjugate<double> >& m);
template void FactorQTZRep<std::complex<double> >::factor( const Matrix_<negator<conjugate<double> > >& m);

template class FactorQTZRep<std::complex<float> >;
template FactorQTZRep<std::complex<float> >::FactorQTZRep( const Matrix_<std::complex<float> >& m, float rcond);
template FactorQTZRep<std::complex<float> >::FactorQTZRep( const Matrix_<negator<std::complex<float> > >& m, float rcond);
template FactorQTZRep<std::complex<float> >::FactorQTZRep( const Matrix_<conjugate<float> >& m, float rcond);
template FactorQTZRep<std::complex<float> >::FactorQTZRep( const Matrix_<negator<conjugate<float> > >& m, float rcond);
template void FactorQTZRep<std::complex<float> >::factor( const Matrix_<std::complex<float> >& m);
template void FactorQTZRep<std::complex<float> >::factor( const Matrix_<negator<std::complex<float> > >& m);
template void FactorQTZRep<std::complex<float> >::factor( const Matrix_<conjugate<float> >& m);
template void FactorQTZRep<std::complex<float> >::factor( const Matrix_<negator<conjugate<float> > >& m);

template SimTK_SIMMATH_EXPORT void FactorQTZ::solve<float>(const Vector_<float>&, Vector_<float>&);
template SimTK_SIMMATH_EXPORT void FactorQTZ::solve<double>(const Vector_<double>&, Vector_<double>&);
template SimTK_SIMMATH_EXPORT void FactorQTZ::solve<std::complex<float> >(const Vector_<std::complex<float> >&, Vector_<std::complex<float> >&);
template SimTK_SIMMATH_EXPORT void FactorQTZ::solve<std::complex<double> >(const Vector_<std::complex<double> >&, Vector_<std::complex<double> >&);
template SimTK_SIMMATH_EXPORT void FactorQTZ::solve<float>(const Matrix_<float>&, Matrix_<float>&);
template SimTK_SIMMATH_EXPORT void FactorQTZ::solve<double>(const Matrix_<double>&, Matrix_<double>&);
template SimTK_SIMMATH_EXPORT void FactorQTZ::solve<std::complex<float> >(const Matrix_<std::complex<float> >&, Matrix_<std::complex<float> >&);
template SimTK_SIMMATH_EXPORT void FactorQTZ::solve<std::complex<double> >(const Matrix_<std::complex<double> >&, Matrix_<std::complex<double> >&);

} // namespace SimTK
