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
 *
 * Solves for singular values and singular vectors 
 */


#include "SimTKcommon.h"

#include "simmath/internal/common.h"
#include "simmath/LinearAlgebra.h"

#include "LapackInterface.h"
#include "FactorSVDRep.h"
#include "WorkSpace.h"
#include "LATraits.h"
#include "LapackConvert.h"

#include <iostream> 
#include <cstdio>
#include <cmath>
#include <complex>


namespace SimTK {

   //////////////////
   // FactorSVDDefault //
   //////////////////
FactorSVDDefault::FactorSVDDefault() {
    isFactored = false;
}
FactorSVDRepBase* FactorSVDDefault::clone() const {
    return( new FactorSVDDefault(*this));
}


   ///////////
   // FactorSVD //
   ///////////
FactorSVD::~FactorSVD() {
    delete rep;
}
FactorSVD::FactorSVD() {
   rep = new FactorSVDDefault();
}
// copy constructor
FactorSVD::FactorSVD( const FactorSVD& c ) {
    rep = c.rep->clone();
}
// copy assignment operator
FactorSVD& FactorSVD::operator=(const FactorSVD& rhs) {
    rep = rhs.rep->clone();
    return *this;
}
int FactorSVD::getRank() {
     return(rep->getRank() );
}
template < typename ELT >
void FactorSVD::solve( const Vector_<ELT>& b, Vector_<ELT>& x ) {
    rep->solve( b, x );
    return;
}
template < class ELT >
void FactorSVD::solve(  const Matrix_<ELT>& b, Matrix_<ELT>& x ) {
    rep->solve(  b, x );
    return;
}

template <typename ELT>
void FactorSVD::inverse( Matrix_<ELT>& inverse ) {
    rep->inverse( inverse );
}
template < class ELT >
FactorSVD::FactorSVD( const Matrix_<ELT>& m ) {

    // if user does not supply rcond set it to max(nRow,nCol)*(eps)^7/8 (similar to matlab)
    int mnmax = (m.nrow() > m.ncol()) ? m.nrow() : m.ncol();
    rep = new FactorSVDRep<typename CNT<ELT>::StdNumber>(m, mnmax*NTraits<typename CNT<ELT>::Precision>::getSignificant()); 
}
template < class ELT >
FactorSVD::FactorSVD( const Matrix_<ELT>& m, double rcond ) {
    rep = new FactorSVDRep<typename CNT<ELT>::StdNumber>(m, rcond);
}
template < class ELT >
FactorSVD::FactorSVD( const Matrix_<ELT>& m, float rcond ) {
    rep = new FactorSVDRep<typename CNT<ELT>::StdNumber>(m, rcond);
}

template < class ELT >
void FactorSVD::factor( const Matrix_<ELT>& m ) {
    delete rep;

    // if user does not supply rcond set it to max(nRow,nCol)*(eps)^7/8 (similar to matlab)
    int mnmax = (m.nrow() > m.ncol()) ? m.nrow() : m.ncol();
    rep = new FactorSVDRep<typename CNT<ELT>::StdNumber>(m, mnmax*NTraits<typename CNT<ELT>::Precision>::getSignificant()); 
}

template < class ELT >
void FactorSVD::factor( const Matrix_<ELT>& m, double rcond ){
    delete rep;
    rep = new FactorSVDRep<typename CNT<ELT>::StdNumber>(m, rcond );
}
template < class ELT >
void FactorSVD::factor( const Matrix_<ELT>& m, float rcond ){
    delete rep;
    rep = new FactorSVDRep<typename CNT<ELT>::StdNumber>(m, rcond );
}

template <class T> 
void FactorSVD::getSingularValuesAndVectors( Vector_<typename CNT<T>::TReal >& values, Matrix_<T>& leftVectors, Matrix_<T>& rightVectors) {

    rep->getSingularValuesAndVectors( values, leftVectors, rightVectors );

    return;
}

template <class T> 
void FactorSVD::getSingularValues( Vector_<T>& values ) {
    rep->getSingularValues( values );
    return;
}
//////////////////
   // FactorSVDRep //
   //////////////////
template <typename T >        // constructor 
    template < typename ELT >
FactorSVDRep<T>::FactorSVDRep( const Matrix_<ELT>& mat, typename CNT<T>::TReal rc):
    nCol(mat.ncol()),  
    nRow(mat.nrow()),
    mn( (mat.nrow() < mat.ncol()) ? mat.nrow() : mat.ncol() ), 
    maxmn( (mat.nrow() > mat.ncol()) ? mat.nrow() : mat.ncol() ), 
    rank(0),
    singularValues(mn),
    inputMatrix(nCol*nRow),
    structure(mat.getMatrixCharacter().getStructure())  {
    
    LapackInterface::getMachineUnderflow( abstol );
    abstol *= 0.5;
     rcond = rc;

    LapackConvert::convertMatrixToLapack( inputMatrix.data, mat );
    isFactored = true;
        
}
template <typename T >
int FactorSVDRep<T>::getRank() {

    if( rank == 0 ) {   // check if SVD has been done
        Vector_<RType> v;    
        getSingularValues( v );
    }
    return( rank );
}
template < class T >
void FactorSVDRep<T>::inverse(  Matrix_<T>& inverse ) {

    Matrix_<T> iden(mn,mn);
    inverse.resize(mn,mn);
    iden = 1.0;
    solve( iden, inverse );
 
}

template <typename T >
FactorSVDRep<T>::~FactorSVDRep() {}

template <typename T >
FactorSVDRepBase* FactorSVDRep<T>::clone() const {
   return( new FactorSVDRep<T>(*this) );
}

template < class T >
void FactorSVDRep<T>::solve( const Vector_<T>& b, Vector_<T> &x ) {

    SimTK_APIARGCHECK_ALWAYS(isFactored ,"FactorSVD","solve",
       "No matrix was passed to FactorSVD. \n"  );

    SimTK_APIARGCHECK2_ALWAYS(b.size()==nRow,"FactorSVD","solve",
       "number of rows in right hand side=%d does not match number of rows in original matrix=%d \n",
        b.size(), nRow );

    if( inputMatrix.size == 0 || b.size() == 0 ) {
        x.resize(0);
    } else { 

        Matrix_<T> m(maxmn,1);

        for(int i=0;i<b.size();i++) {
            m(i,0) = b(i);
        }
        Matrix_<T> r(nCol, 1 );
        doSolve( m, r );
        x.copyAssign(r);
    }
    return;

}

template < class T >
void FactorSVDRep<T>::solve( const Matrix_<T>& b, Matrix_<T> &x ) {

    SimTK_APIARGCHECK_ALWAYS(isFactored ,"FactorSVD","solve",
       "No matrix was passed to FactorSVD. \n"  );

    SimTK_APIARGCHECK2_ALWAYS(b.nrow()==nRow,"FactorSVD","solve",
       "number of rows in right hand side=%d does not match number of rows in original matrix=%d \n",
        b.nrow(), nRow );

    if( inputMatrix.size == 0 || b.nelt() == 0 ) {
        x.resize(0,0);
    } else { 
        x.resize( nCol, b.ncol() );
        Matrix_<T> tb;
        tb.resize(maxmn, b.ncol() );
        for(int j=0;j<b.ncol();j++) for(int i=0;i<b.nrow();i++) tb(i,j) = b(i,j);
        doSolve(tb, x);
    }


}

template <typename T >
void FactorSVDRep<T>::doSolve(  Matrix_<T>& b, Matrix_<T>& x) {
    int i,j;
    int info;
    typedef typename CNT<T>::TReal RealType;

    if( b.nelt() == 0 || inputMatrix.size == 0) return;

    TypedWorkSpace<T> tempMatrix = inputMatrix;

    x.resize(nCol, b.ncol() );
    Matrix_<T> tb;
    tb.resize(maxmn, b.ncol() );
    for(j=0;j<b.ncol();j++) for(i=0;i<b.nrow();i++) tb(i,j) = b(i,j);

    LapackInterface::gelss<T>( nRow, nCol, mn, b.ncol(), tempMatrix.data, nRow, &tb(0,0), 
                      tb.nrow(), singularValues.data, rcond, rank, info  );

    if( info > 0 ) {
        SimTK_THROW2( SimTK::Exception::ConvergedFailed,
        "FactorSVD::solve", 
        "divide and conquer singular value decomposition" );
    }
    
    for(j=0;j<b.ncol();j++) for(i=0;i<nCol;i++) x(i,j) = tb(i,j);

}

template < class T >
void FactorSVDRep<T>::getSingularValuesAndVectors( Vector_<RType>& values,  Matrix_<T>& leftVectors, Matrix_<T>& rightVectors ) {
   

    if( inputMatrix.size == 0 ) {
        leftVectors.resize(0,0);
        rightVectors.resize(0,0);
        values.resize(0);
        return;
    } else {
        leftVectors.resize(nRow,nRow);
        rightVectors.resize(nCol,nCol);
        values.resize(mn);
    }

    computeSVD( true, &values(0), &leftVectors(0,0), &rightVectors(0,0) );

    return;
}
template < class T >
void FactorSVDRep<T>::getSingularValues( Vector_<RType>& values ) {

    if( inputMatrix.size == 0 ) {
        values.resize(0);
        return;
    } else {
        values.resize(mn);
    }

    computeSVD( false, &values(0), NULL, NULL );

    return;
}
template < class T >
void FactorSVDRep<T>::computeSVD( bool computeVectors, RType* values, T* leftVectors, T* rightVectors ) {

    int info;
    char jobz;


    if( computeVectors ) {
        jobz = 'A';
    } else {
        jobz = 'N';
    }

    TypedWorkSpace<T> tempMatrix = inputMatrix;
    LapackInterface::gesdd<T>(jobz, nRow,nCol,tempMatrix.data, nRow, values,
           leftVectors, nRow, rightVectors, nCol, info);

    for(int i=0, rank=0;i<mn;i++) {
        if( values[i] > rcond*values[0] ) rank++;
    } 

    if( info > 0 ) {
        SimTK_THROW2( SimTK::Exception::ConvergedFailed,
        "FactorSVD::getSingularValuesAndVectors", 
        "divide and conquer singular value decomposition" );
    }

    return;
}

// instantiate
template SimTK_SIMMATH_EXPORT FactorSVD::FactorSVD( const Matrix_<double>& m );
template SimTK_SIMMATH_EXPORT FactorSVD::FactorSVD( const Matrix_<float>& m );
template SimTK_SIMMATH_EXPORT FactorSVD::FactorSVD( const Matrix_<std::complex<float> >& m );
template SimTK_SIMMATH_EXPORT FactorSVD::FactorSVD( const Matrix_<std::complex<double> >& m );
template SimTK_SIMMATH_EXPORT FactorSVD::FactorSVD( const Matrix_<conjugate<float> >& m );
template SimTK_SIMMATH_EXPORT FactorSVD::FactorSVD( const Matrix_<conjugate<double> >& m );
template SimTK_SIMMATH_EXPORT FactorSVD::FactorSVD( const Matrix_<negator< double> >& m );
template SimTK_SIMMATH_EXPORT FactorSVD::FactorSVD( const Matrix_<negator< float> >& m );
template SimTK_SIMMATH_EXPORT FactorSVD::FactorSVD( const Matrix_<negator< std::complex<float> > >& m );
template SimTK_SIMMATH_EXPORT FactorSVD::FactorSVD( const Matrix_<negator< std::complex<double> > >& m );
template SimTK_SIMMATH_EXPORT FactorSVD::FactorSVD( const Matrix_<negator< conjugate<float> > >& m );
template SimTK_SIMMATH_EXPORT FactorSVD::FactorSVD( const Matrix_<negator< conjugate<double> > >& m );

template SimTK_SIMMATH_EXPORT FactorSVD::FactorSVD( const Matrix_<double>& m, double rcond );
template SimTK_SIMMATH_EXPORT FactorSVD::FactorSVD( const Matrix_<float>& m, float rcond );
template SimTK_SIMMATH_EXPORT FactorSVD::FactorSVD( const Matrix_<std::complex<float> >& m, float rcond );
template SimTK_SIMMATH_EXPORT FactorSVD::FactorSVD( const Matrix_<std::complex<double> >& m, double rcond );
template SimTK_SIMMATH_EXPORT FactorSVD::FactorSVD( const Matrix_<conjugate<float> >& m, float rcond );
template SimTK_SIMMATH_EXPORT FactorSVD::FactorSVD( const Matrix_<conjugate<double> >& m, double rcond );
template SimTK_SIMMATH_EXPORT FactorSVD::FactorSVD( const Matrix_<negator< double> >& m, double rcond );
template SimTK_SIMMATH_EXPORT FactorSVD::FactorSVD( const Matrix_<negator< float> >& m, float rcond );
template SimTK_SIMMATH_EXPORT FactorSVD::FactorSVD( const Matrix_<negator< std::complex<float> > >& m, float rcond );
template SimTK_SIMMATH_EXPORT FactorSVD::FactorSVD( const Matrix_<negator< std::complex<double> > >& m, double rcond );
template SimTK_SIMMATH_EXPORT FactorSVD::FactorSVD( const Matrix_<negator< conjugate<float> > >& m, float rcond );
template SimTK_SIMMATH_EXPORT FactorSVD::FactorSVD( const Matrix_<negator< conjugate<double> > >& m, double rcond );

template SimTK_SIMMATH_EXPORT void FactorSVD::getSingularValues<float >(Vector_<float>&  );
template SimTK_SIMMATH_EXPORT void FactorSVD::getSingularValues<double >(Vector_<double>&  );

template SimTK_SIMMATH_EXPORT void FactorSVD::getSingularValuesAndVectors<float >(Vector_<float>&, Matrix_<float>&, Matrix_<float>&  );
template SimTK_SIMMATH_EXPORT void FactorSVD::getSingularValuesAndVectors<double >(Vector_<double>&, Matrix_<double>&, Matrix_<double>&  );
template SimTK_SIMMATH_EXPORT void FactorSVD::getSingularValuesAndVectors<std::complex<float> >(Vector_<float>&, Matrix_<std::complex<float> >&, Matrix_<std::complex<float> >&  );
template SimTK_SIMMATH_EXPORT void FactorSVD::getSingularValuesAndVectors<std::complex<double> >(Vector_<double>&, Matrix_<std::complex<double> >&, Matrix_<std::complex<double> >&  );

template SimTK_SIMMATH_EXPORT void FactorSVD::factor( const Matrix_<double>& m );
template SimTK_SIMMATH_EXPORT void FactorSVD::factor( const Matrix_<float>& m );
template SimTK_SIMMATH_EXPORT void FactorSVD::factor( const Matrix_<std::complex<float> >& m );
template SimTK_SIMMATH_EXPORT void FactorSVD::factor( const Matrix_<std::complex<double> >& m );
template SimTK_SIMMATH_EXPORT void FactorSVD::factor( const Matrix_<conjugate<float> >& m );
template SimTK_SIMMATH_EXPORT void FactorSVD::factor( const Matrix_<conjugate<double> >& m );
template SimTK_SIMMATH_EXPORT void FactorSVD::factor( const Matrix_<negator< double> >& m );
template SimTK_SIMMATH_EXPORT void FactorSVD::factor( const Matrix_<negator< float> >& m );
template SimTK_SIMMATH_EXPORT void FactorSVD::factor( const Matrix_<negator< std::complex<float> > >& m );
template SimTK_SIMMATH_EXPORT void FactorSVD::factor( const Matrix_<negator< std::complex<double> > >& m );
template SimTK_SIMMATH_EXPORT void FactorSVD::factor( const Matrix_<negator< conjugate<float> > >& m );
template SimTK_SIMMATH_EXPORT void FactorSVD::factor( const Matrix_<negator< conjugate<double> > >& m );

template SimTK_SIMMATH_EXPORT void FactorSVD::factor( const Matrix_<double>& m, double rcond );
template SimTK_SIMMATH_EXPORT void FactorSVD::factor( const Matrix_<float>& m, float rcond );
template SimTK_SIMMATH_EXPORT void FactorSVD::factor( const Matrix_<std::complex<float> >& m, float rcond );
template SimTK_SIMMATH_EXPORT void FactorSVD::factor( const Matrix_<std::complex<double> >& m, double rcond );
template SimTK_SIMMATH_EXPORT void FactorSVD::factor( const Matrix_<conjugate<float> >& m, float rcond );
template SimTK_SIMMATH_EXPORT void FactorSVD::factor( const Matrix_<conjugate<double> >& m, double rcond );
template SimTK_SIMMATH_EXPORT void FactorSVD::factor( const Matrix_<negator< double> >& m, double rcond );
template SimTK_SIMMATH_EXPORT void FactorSVD::factor( const Matrix_<negator< float> >& m, float rcond );
template SimTK_SIMMATH_EXPORT void FactorSVD::factor( const Matrix_<negator< std::complex<float> > >& m, float rcond );
template SimTK_SIMMATH_EXPORT void FactorSVD::factor( const Matrix_<negator< std::complex<double> > >& m, double rcond );
template SimTK_SIMMATH_EXPORT void FactorSVD::factor( const Matrix_<negator< conjugate<float> > >& m, float rcond );
template SimTK_SIMMATH_EXPORT void FactorSVD::factor( const Matrix_<negator< conjugate<double> > >& m, double rcond );

template SimTK_SIMMATH_EXPORT void FactorSVD::inverse<float>(Matrix_<float>&);
template SimTK_SIMMATH_EXPORT void FactorSVD::inverse<double>(Matrix_<double>&);
template SimTK_SIMMATH_EXPORT void FactorSVD::inverse<std::complex<float> >(Matrix_<std::complex<float> >&);
template SimTK_SIMMATH_EXPORT void FactorSVD::inverse<std::complex<double> >(Matrix_<std::complex<double> >&);

template SimTK_SIMMATH_EXPORT void FactorSVD::solve<float>(const Vector_<float>&, Vector_<float>&);
template SimTK_SIMMATH_EXPORT void FactorSVD::solve<double>(const Vector_<double>&, Vector_<double>&);
template SimTK_SIMMATH_EXPORT void FactorSVD::solve<std::complex<float> >(const Vector_<std::complex<float> >&, Vector_<std::complex<float> >&);
template SimTK_SIMMATH_EXPORT void FactorSVD::solve<std::complex<double> >(const Vector_<std::complex<double> >&, Vector_<std::complex<double> >&);
template SimTK_SIMMATH_EXPORT void FactorSVD::solve<float>(const Matrix_<float>&, Matrix_<float>&);
template SimTK_SIMMATH_EXPORT void FactorSVD::solve<double>(const Matrix_<double>&, Matrix_<double>&);
template SimTK_SIMMATH_EXPORT void FactorSVD::solve<std::complex<float> >(const Matrix_<std::complex<float> >&, Matrix_<std::complex<float> >&);
template SimTK_SIMMATH_EXPORT void FactorSVD::solve<std::complex<double> >(const Matrix_<std::complex<double> >&, Matrix_<std::complex<double> >&);

template class FactorSVDRep<double>;
template FactorSVDRep<double>::FactorSVDRep( const Matrix_<double>& m, double rcond);
template FactorSVDRep<double>::FactorSVDRep( const Matrix_<negator<double> >& m, double rcond);

template class FactorSVDRep<float>;
template FactorSVDRep<float>::FactorSVDRep( const Matrix_<float>& m, float rcond );
template FactorSVDRep<float>::FactorSVDRep( const Matrix_<negator<float> >& m, float rcond );

template class FactorSVDRep<std::complex<double> >;
template FactorSVDRep<std::complex<double> >::FactorSVDRep( const Matrix_<std::complex<double> >& m, double rcond );
template FactorSVDRep<std::complex<double> >::FactorSVDRep( const Matrix_<negator<std::complex<double> > >& m, double rcond );
template FactorSVDRep<std::complex<double> >::FactorSVDRep( const Matrix_<conjugate<double> >& m, double rcond );
template FactorSVDRep<std::complex<double> >::FactorSVDRep( const Matrix_<negator<conjugate<double> > >& m, double rcond );

template class FactorSVDRep<std::complex<float> >;
template FactorSVDRep<std::complex<float> >::FactorSVDRep( const Matrix_<std::complex<float> >& m, float rcond );
template FactorSVDRep<std::complex<float> >::FactorSVDRep( const Matrix_<negator<std::complex<float> > >& m, float rcond );
template FactorSVDRep<std::complex<float> >::FactorSVDRep( const Matrix_<conjugate<float> >& m, float rcond );
template FactorSVDRep<std::complex<float> >::FactorSVDRep( const Matrix_<negator<conjugate<float> > >& m, float rcond );

} // namespace SimTK
