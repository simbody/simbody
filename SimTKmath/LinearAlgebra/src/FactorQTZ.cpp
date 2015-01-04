/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-12 Stanford University and the Authors.        *
 * Authors: Jack Middleton                                                    *
 * Contributors: Michael Sherman                                              *
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
 *
 * Factors systems of linear algebra equations.
 */

#include "SimTKcommon.h"

#include "simmath/internal/common.h"
#include "simmath/LinearAlgebra.h"

#include "LapackInterface.h"
#include "FactorQTZRep.h"
#include "WorkSpace.h"
#include "LATraits.h"
#include "LapackConvert.h"
#include "SimTKcommon.h"

#include <iostream> 
#include <cmath>
#include <complex>


namespace SimTK {

   //////////////////////
   // FactorQTZDefault //
   //////////////////////
FactorQTZDefault::FactorQTZDefault() {
    isFactored = false; 
}
FactorQTZRepBase* FactorQTZDefault::clone() const {
    return( new FactorQTZDefault(*this));  
}

   ///////////////
   // FactorQTZ //
   ///////////////
FactorQTZ::~FactorQTZ() {
    delete rep;
}
// default constructor
FactorQTZ::FactorQTZ() {
    rep = new FactorQTZDefault();
}
// copy constructor
FactorQTZ::FactorQTZ( const FactorQTZ& c ) {
    rep = c.rep->clone();
}
// copy assignment operator
FactorQTZ& FactorQTZ::operator=(const FactorQTZ& rhs) {
    rep = rhs.rep->clone();
    return *this;
}

template <typename ELT>
void FactorQTZ::inverse( Matrix_<ELT>& inverse ) const {
    rep->inverse( inverse );
}
template < class ELT >
void FactorQTZ::factor( const Matrix_<ELT>& m ){
    delete rep;
  
    // if user does not supply rcond set it to max(nRow,nCol)*(eps)^7/8 (similar to matlab)
    int mnmax = (m.nrow() > m.ncol()) ? m.nrow() : m.ncol();
    rep = new FactorQTZRep<typename CNT<ELT>::StdNumber>(m, mnmax*NTraits<typename CNT<ELT>::Precision>::getSignificant());
}
template < class ELT >
void FactorQTZ::factor( const Matrix_<ELT>& m, double rcond ){
    delete rep;
    rep = new FactorQTZRep<typename CNT<ELT>::StdNumber>(m, rcond );
}
template < class ELT >
void FactorQTZ::factor( const Matrix_<ELT>& m, float rcond ){
    delete rep;
    rep = new FactorQTZRep<typename CNT<ELT>::StdNumber>(m, rcond );
}
template < class ELT >
FactorQTZ::FactorQTZ( const Matrix_<ELT>& m ) {

    // if user does not supply rcond set it to max(nRow,nCol)*(eps)^7/8 (similar to matlab)
    int mnmax = (m.nrow() > m.ncol()) ? m.nrow() : m.ncol();
    rep = new FactorQTZRep<typename CNT<ELT>::StdNumber>(m, mnmax*NTraits<typename CNT<ELT>::Precision>::getSignificant()); 
}
template < class ELT >
FactorQTZ::FactorQTZ( const Matrix_<ELT>& m, double rcond ) {
    rep = new FactorQTZRep<typename CNT<ELT>::StdNumber>(m, rcond); 
}
template < class ELT >
FactorQTZ::FactorQTZ( const Matrix_<ELT>& m, float rcond ) {
    rep = new FactorQTZRep<typename CNT<ELT>::StdNumber>(m, rcond); 
}

int FactorQTZ::getRank() const {
    return(rep->rank);
}
double FactorQTZ::getRCondEstimate() const {
    return (rep->actualRCond);
}
template < typename ELT >
void FactorQTZ::solve( const Vector_<ELT>& b, Vector_<ELT>& x ) const {
    rep->solve( b, x );
    return;
}
template < class ELT >
void FactorQTZ::solve(  const Matrix_<ELT>& b, Matrix_<ELT>& x ) const {
    rep->solve(  b, x );
    return;
}

   /////////////////
   // FactorQTZRep //
   /////////////////
template <typename T >
FactorQTZRep<T>::FactorQTZRep() 
:   mn(0),
    maxmn(0),
    nRow(0),
    nCol(0),
    scaleLinSys(false),
    linSysScaleF(NTraits<typename CNT<T>::Precision>::getNaN()),
    anrm(NTraits<typename CNT<T>::Precision>::getNaN()),
    rcond(NTraits<typename CNT<T>::Precision>::getSignificant()),
    pivots(0),
    qtz(0),
    tauGEQP3(0),
    tauORMQR(0)
{ 
} 

template <typename T >
    template < typename ELT >
FactorQTZRep<T>::FactorQTZRep( const Matrix_<ELT>& mat, typename CNT<T>::TReal rc) 
:   mn( (mat.nrow() < mat.ncol()) ? mat.nrow() : mat.ncol() ),
    maxmn( (mat.nrow() > mat.ncol()) ? mat.nrow() : mat.ncol() ),
    nRow( mat.nrow() ),
    nCol( mat.ncol() ),
    scaleLinSys(false),
    linSysScaleF(NTraits<typename CNT<T>::Precision>::getNaN()),
    anrm(NTraits<typename CNT<T>::Precision>::getNaN()),
    rcond(rc),
    pivots(mat.ncol()),
    qtz( mat.nrow()*mat.ncol() ),
    tauGEQP3(mn),
    tauORMQR(mn)    
{ 
    for(int i=0; i<mat.ncol(); ++i) 
        pivots.data[i] = 0;
    FactorQTZRep<T>::factor( mat );
    isFactored = true;
}

template <typename T >
FactorQTZRepBase* FactorQTZRep<T>::clone() const {
   return( new FactorQTZRep<T>(*this) );
}

template <typename T >
FactorQTZRep<T>::~FactorQTZRep() {}

template < class T >
void FactorQTZRep<T>::solve( const Vector_<T>& b, Vector_<T> &x ) const {

    SimTK_APIARGCHECK_ALWAYS(isFactored ,"FactorQTZ","solve",
       "No matrix was passed to FactorQTZ. \n"  );

    SimTK_APIARGCHECK2_ALWAYS(b.size()==nRow,"FactorQTZ","solve",
       "number of rows in right hand side=%d does not match number of rows in original matrix=%d \n", 
        b.size(), nRow );


    Matrix_<T> m(maxmn,1);

    for(int i=0;i<b.size();i++) {
        m(i,0) = b(i);
    }
    Matrix_<T> r(nCol, 1 );
    doSolve( m, r );
    x.copyAssign(r);
}

template <typename T >
void FactorQTZRep<T>::solve(  const Matrix_<T>& b, Matrix_<T>& x ) const {
    SimTK_APIARGCHECK_ALWAYS(isFactored ,"FactorQTZ","solve",
       "No matrix was passed to FactorQTZ. \n"  );

    SimTK_APIARGCHECK2_ALWAYS(b.nrow()==nRow,"FactorQTZ","solve",
       "number of rows in right hand side=%d does not match number of rows in original matrix=%d \n", 
        b.nrow(), nRow );

    x.resize(nCol, b.ncol() );
    Matrix_<T> tb;
    tb.resize(maxmn, b.ncol() );
    for(int j=0;j<b.ncol();j++) for(int i=0;i<b.nrow();i++) tb(i,j) = b(i,j);
    doSolve(tb, x);
}

template <typename T >
void FactorQTZRep<T>::doSolve(Matrix_<T>& b, Matrix_<T>& x) const {
    int info;
    typedef typename CNT<T>::TReal RealType;
    RealType bnrm, smlnum, bignum;
    int nrhs = b.ncol();
    int n = nCol;
    int m = nRow;
    typename CNT<T>::TReal rhsScaleF; // scale factor applied to right hand side
    bool scaleRHS = false; // true if right hand side should be scaled

    if (rank == 0) return;

    // Ask the experts for their optimal workspace sizes. The size parameters
    // here must match the calls below.
    T workSz;
    LapackInterface::ormqr<T>('L', 'T', nRow, b.ncol(), mn, 0, nRow, 
                               0, 0, b.nrow(), &workSz, -1, info );
    const int lwork1 = (int)NTraits<T>::real(workSz);

    LapackInterface::ormrz<T>('L', 'T', nCol, b.ncol(), rank, nCol-rank, 
                              0, nRow, 0, 0, 
                              b.nrow(), &workSz, -1, info );
    const int lwork2 = (int)NTraits<T>::real(workSz);
    
    TypedWorkSpace<T> work(std::max(lwork1, lwork2));

    // compute norm of RHS
    bnrm = (RealType)LapackInterface::lange<T>('M', m, nrhs, &b(0,0), b.nrow());

    LapackInterface::getMachinePrecision<RealType>(smlnum, bignum);
 
    // compute scale for RHS
    if (bnrm > Zero && bnrm < smlnum) {
        scaleRHS = true;
        rhsScaleF = smlnum;
    } else if (bnrm > bignum) {
        scaleRHS = true;
        rhsScaleF = bignum;
    }


    if (scaleRHS) {   // apply scale factor to RHS
        LapackInterface::lascl<T>('G', 0, 0, bnrm, rhsScaleF, b.nrow(), nrhs, 
                                  &b(0,0), b.nrow(), info ); 
    }
    // b1 = Q'*b0
    LapackInterface::ormqr<T>('L', 'T', nRow, b.ncol(), mn, qtz.data, 
                              nRow, tauGEQP3.data, &b(0,0), b.nrow(), 
                              work.data, work.size, info );
    // b2 = T^-1*b1 = T^-1 * Q' * b0
    LapackInterface::trsm<T>('L', 'U', 'N', 'N', rank, b.ncol(), 1.0, 
                             qtz.data, nRow, &b(0,0), b.nrow() );

    //  zero out elements of RHS for rank deficient systems
    for(int j = 0; j<nrhs; ++j) {
        for(int i = rank; i<n; ++i)
            b(i,j) = 0;
    }
   
    if (rank < nCol) {
        // b3 = Z'*b2 = Z'*T^-1*Q'*b0
        LapackInterface::ormrz<T>('L', 'T', nCol, b.ncol(), rank, nCol-rank, 
                                  qtz.data, nRow, tauORMQR.data, &b(0,0), 
                                  b.nrow(), work.data, work.size, info );
    }

    // adjust for pivoting
    Array_<T> b_pivot(n);
    for(int j = 0; j<nrhs; ++j) {
        for(int i = 0; i<n; ++i)
            b_pivot[pivots.data[i]-1] = b(i,j);

        LapackInterface::copy<T>(n, b_pivot.begin(), 1, &x(0,j), 1 );
    }

    // compensate for scaling of linear system 
    if (scaleLinSys) { 
        LapackInterface::lascl<T>('g', 0, 0, anrm, linSysScaleF, nCol, x.ncol(),
                                  &x(0,0), nCol, info );
    }

    // compensate for scaling of RHS 
    if (scaleRHS) { 
        LapackInterface::lascl<T>('g', 0, 0, bnrm, rhsScaleF, nCol, x.ncol(), 
                                  &x(0,0), nCol, info);
    }
}


template < class T >
void FactorQTZRep<T>::inverse(  Matrix_<T>& inverse ) const {
    Matrix_<T> iden(mn,mn);
    inverse.resize(mn,mn);
    iden = 1.0;
    doSolve( iden, inverse );

    return;
}

template <class T> 
    template<typename ELT>
void FactorQTZRep<T>::factor(const Matrix_<ELT>&mat )  {
    SimTK_APIARGCHECK2_ALWAYS(mat.nelt() > 0,"FactorQTZ","factor",
       "Can't factor a matrix that has a zero dimension -- got %d X %d.",
       (int)mat.nrow(), (int)mat.ncol());


    // allocate and initialize the matrix we pass to LAPACK
    // converts (negated,conjugated etc.) to LAPACK format 
    LapackConvert::convertMatrixToLapack( qtz.data, mat );

    int info;
    const int smallestSingularValue = 2;
    const int largestSingularValue = 1;
    typedef typename CNT<T>::TReal  RealType;
    RealType smlnum, bignum, smin, smax;

    if (mat.nelt() == 0) return;


    // Compute optimal size for work space for dtzrzf and dgepq3. The
    // arguments here should match the calls below, although we'll use maxRank
    // rather than rank since we don't know the rank yet.
    T workSz;
    const int maxRank = std::min(nRow, nCol);
    LapackInterface::tzrzf<T>(maxRank, nCol, 0, nRow, 0, &workSz, -1, info);
    const int lwork1 = (int)NTraits<T>::real(workSz);

    LapackInterface::geqp3<T>(nRow, nCol, 0, nRow, 0, 0, &workSz, -1, info);
    const int lwork2 = (int)NTraits<T>::real(workSz);
   
    TypedWorkSpace<T> work(std::max(lwork1, lwork2));

    LapackInterface::getMachinePrecision<RealType>( smlnum, bignum);

    // scale the input system of equations
    anrm = (RealType)LapackInterface::lange<T>('M', nRow, nCol, qtz.data, nRow);

    if (anrm > 0 && anrm < smlnum) {
        scaleLinSys = true;
        linSysScaleF = smlnum;
    } else if( anrm > bignum )  {
        scaleLinSys = true;
        linSysScaleF = bignum;
    } 
    if (anrm == 0) { // matrix all zeros
        rank = 0;
    } else {
        if (scaleLinSys) {
           LapackInterface::lascl<T>('G', 0, 0, anrm, linSysScaleF, nRow, nCol,
                                     qtz.data, nRow, info );
        }

        // compute QR factorization with column pivoting: A = Q * R
        // Q * R is returned in qtz.data
        LapackInterface::geqp3<T>(nRow, nCol, qtz.data, nRow, pivots.data, 
                                  tauGEQP3.data, work.data, work.size, info );

        // compute Rank

        smax = CNT<T>::abs( qtz.data[0] );
        smin = smax;
        if( CNT<T>::abs(qtz.data[0]) == 0 ) {
            rank = 0;
            actualRCond = 0;
        } else {
            T s1,s2,c1,c2;
            RealType smaxpr,sminpr;

            // Determine rank using incremental condition estimate
            Array_<T> xSmall(mn), xLarge(mn); // temporaries
            xSmall[0] = xLarge[0] = 1;
            for (rank=1,smaxpr=0.0,sminpr=1.0; 
                 rank<mn && smaxpr*rcond < sminpr; ) 
            {
                LapackInterface::laic1<T>(smallestSingularValue, rank, 
                    xSmall.begin(), smin, &qtz.data[rank*nRow], 
                    qtz.data[(rank*nRow)+rank], sminpr, s1, c1);

                LapackInterface::laic1<T>(largestSingularValue, rank, 
                    xLarge.begin(), smax, &qtz.data[rank*nRow], 
                    qtz.data[(rank*nRow)+rank], smaxpr, s2, c2);

                if (smaxpr*rcond < sminpr) {
                    for(int i=0; i<rank; i++) {
                         xSmall[i] *= s1;
                         xLarge[i] *= s2;
                    }
                    xSmall[rank] = c1;
                    xLarge[rank] = c2;
                    smin = sminpr;
                    smax = smaxpr;
                    actualRCond = (double)(smin/smax);
                    rank++;
                }
            } 

            // R => T * Z
            // Matrix is rank deficient so complete the orthogonalization by 
            // applying orthogonal transformations from the right to 
            // factor R from an upper trapezoidal to an upper triangular matrix
            // T is returned in qtz.data and Z is returned in tauORMQR.data
            if (rank < nCol) {
                LapackInterface::tzrzf<T>(rank, nCol, qtz.data, nRow, 
                                          tauORMQR.data, work.data, 
                                          work.size, info);
            }
        }
    }
}

// instantiate
template SimTK_SIMMATH_EXPORT FactorQTZ::FactorQTZ( const Matrix_<double>& m );
template SimTK_SIMMATH_EXPORT FactorQTZ::FactorQTZ( const Matrix_<float>& m );
template SimTK_SIMMATH_EXPORT FactorQTZ::FactorQTZ( const Matrix_<std::complex<float> >& m );
template SimTK_SIMMATH_EXPORT FactorQTZ::FactorQTZ( const Matrix_<std::complex<double> >& m );
template SimTK_SIMMATH_EXPORT FactorQTZ::FactorQTZ( const Matrix_<conjugate<float> >& m );
template SimTK_SIMMATH_EXPORT FactorQTZ::FactorQTZ( const Matrix_<conjugate<double> >& m );
template SimTK_SIMMATH_EXPORT FactorQTZ::FactorQTZ( const Matrix_<negator< double> >& m );
template SimTK_SIMMATH_EXPORT FactorQTZ::FactorQTZ( const Matrix_<negator< float> >& m );
template SimTK_SIMMATH_EXPORT FactorQTZ::FactorQTZ( const Matrix_<negator< std::complex<float> > >& m );
template SimTK_SIMMATH_EXPORT FactorQTZ::FactorQTZ( const Matrix_<negator< std::complex<double> > >& m );
template SimTK_SIMMATH_EXPORT FactorQTZ::FactorQTZ( const Matrix_<negator< conjugate<float> > >& m );
template SimTK_SIMMATH_EXPORT FactorQTZ::FactorQTZ( const Matrix_<negator< conjugate<double> > >& m );

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

template SimTK_SIMMATH_EXPORT void FactorQTZ::factor( const Matrix_<double>& m );
template SimTK_SIMMATH_EXPORT void FactorQTZ::factor( const Matrix_<float>& m );
template SimTK_SIMMATH_EXPORT void FactorQTZ::factor( const Matrix_<std::complex<float> >& m );
template SimTK_SIMMATH_EXPORT void FactorQTZ::factor( const Matrix_<std::complex<double> >& m );
template SimTK_SIMMATH_EXPORT void FactorQTZ::factor( const Matrix_<conjugate<float> >& m );
template SimTK_SIMMATH_EXPORT void FactorQTZ::factor( const Matrix_<conjugate<double> >& m );
template SimTK_SIMMATH_EXPORT void FactorQTZ::factor( const Matrix_<negator< double> >& m );
template SimTK_SIMMATH_EXPORT void FactorQTZ::factor( const Matrix_<negator< float> >& m );
template SimTK_SIMMATH_EXPORT void FactorQTZ::factor( const Matrix_<negator< std::complex<float> > >& m );
template SimTK_SIMMATH_EXPORT void FactorQTZ::factor( const Matrix_<negator< std::complex<double> > >& m );
template SimTK_SIMMATH_EXPORT void FactorQTZ::factor( const Matrix_<negator< conjugate<float> > >& m );
template SimTK_SIMMATH_EXPORT void FactorQTZ::factor( const Matrix_<negator< conjugate<double> > >& m );

template SimTK_SIMMATH_EXPORT void FactorQTZ::factor( const Matrix_<double>& m, double rcond );
template SimTK_SIMMATH_EXPORT void FactorQTZ::factor( const Matrix_<float>& m, float rcond );
template SimTK_SIMMATH_EXPORT void FactorQTZ::factor( const Matrix_<std::complex<float> >& m, float rcond );
template SimTK_SIMMATH_EXPORT void FactorQTZ::factor( const Matrix_<std::complex<double> >& m, double rcond );
template SimTK_SIMMATH_EXPORT void FactorQTZ::factor( const Matrix_<conjugate<float> >& m, float rcond );
template SimTK_SIMMATH_EXPORT void FactorQTZ::factor( const Matrix_<conjugate<double> >& m, double rcond );
template SimTK_SIMMATH_EXPORT void FactorQTZ::factor( const Matrix_<negator< double> >& m, double rcond );
template SimTK_SIMMATH_EXPORT void FactorQTZ::factor( const Matrix_<negator< float> >& m, float rcond );
template SimTK_SIMMATH_EXPORT void FactorQTZ::factor( const Matrix_<negator< std::complex<float> > >& m, float rcond );
template SimTK_SIMMATH_EXPORT void FactorQTZ::factor( const Matrix_<negator< std::complex<double> > >& m, double rcond );
template SimTK_SIMMATH_EXPORT void FactorQTZ::factor( const Matrix_<negator< conjugate<float> > >& m, float rcond );
template SimTK_SIMMATH_EXPORT void FactorQTZ::factor( const Matrix_<negator< conjugate<double> > >& m, double rcond );

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

template SimTK_SIMMATH_EXPORT void FactorQTZ::solve<float>(const Vector_<float>&, Vector_<float>&) const;
template SimTK_SIMMATH_EXPORT void FactorQTZ::solve<double>(const Vector_<double>&, Vector_<double>&) const;
template SimTK_SIMMATH_EXPORT void FactorQTZ::solve<std::complex<float> >(const Vector_<std::complex<float> >&, Vector_<std::complex<float> >&) const;
template SimTK_SIMMATH_EXPORT void FactorQTZ::solve<std::complex<double> >(const Vector_<std::complex<double> >&, Vector_<std::complex<double> >&) const;
template SimTK_SIMMATH_EXPORT void FactorQTZ::solve<float>(const Matrix_<float>&, Matrix_<float>&) const;
template SimTK_SIMMATH_EXPORT void FactorQTZ::solve<double>(const Matrix_<double>&, Matrix_<double>&) const;
template SimTK_SIMMATH_EXPORT void FactorQTZ::solve<std::complex<float> >(const Matrix_<std::complex<float> >&, Matrix_<std::complex<float> >&) const;
template SimTK_SIMMATH_EXPORT void FactorQTZ::solve<std::complex<double> >(const Matrix_<std::complex<double> >&, Matrix_<std::complex<double> >&) const;
template SimTK_SIMMATH_EXPORT void FactorQTZ::inverse<float>(Matrix_<float>&) const;
template SimTK_SIMMATH_EXPORT void FactorQTZ::inverse<double>(Matrix_<double>&) const;
template SimTK_SIMMATH_EXPORT void FactorQTZ::inverse<std::complex<float> >(Matrix_<std::complex<float> >&) const;
template SimTK_SIMMATH_EXPORT void FactorQTZ::inverse<std::complex<double> >(Matrix_<std::complex<double> >&) const;

} // namespace SimTK
