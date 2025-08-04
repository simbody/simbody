/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-25 Stanford University and the Authors.        *
 * Authors: Alexander Beattie                                                 *
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
 * Solves LLT (Cholesky) factorization of real symmetric positive-definite
 * matrices
 */

#include "SimTKcommon.h"
#include "LapackInterface.h"
#include "LapackConvert.h"
#include "simmath/internal/common.h"
#include "simmath/LinearAlgebra.h"
#include "LapackInterface.h"
#include "FactorLLTRep.h"
#include "WorkSpace.h"

namespace SimTK {

///////////////
// FactorLLTDefault
///////////////
FactorLLTDefault::FactorLLTDefault() { isFactored = false; }
FactorLLTRepBase* FactorLLTDefault::clone() const {
    return new FactorLLTDefault(*this);
}

///////////
// FactorLLT //
///////////
FactorLLT::~FactorLLT() {
    delete rep;
}
// default constructor
FactorLLT::FactorLLT() {
    rep = new FactorLLTDefault();
}

// copy constructor
FactorLLT::FactorLLT( const FactorLLT& c ) {
    rep = c.rep->clone();
}
// copy assignment operator
FactorLLT& FactorLLT::operator=(const FactorLLT& rhs) {
    rep = rhs.rep->clone();
    return *this;
}

template <typename ELT>
void FactorLLT::inverse( Matrix_<ELT>& inverse ) const {
    rep->inverse( inverse );
}


template < class ELT >
FactorLLT::FactorLLT( const Matrix_<ELT>& m ) {
    rep = new FactorLLTRep<typename CNT<ELT>::StdNumber>(m);
}

template < class ELT >
void FactorLLT::factor( const Matrix_<ELT>& m ) {
    delete rep;
    rep = new FactorLLTRep<typename CNT<ELT>::StdNumber>(m);
}

template < typename ELT >
void FactorLLT::solve( const Vector_<ELT>& b, Vector_<ELT>& x ) const {
    rep->solve( b, x );
    return;
}
template < class ELT >
void FactorLLT::solve(  const Matrix_<ELT>& b, Matrix_<ELT>& x ) const {
    rep->solve(  b, x );
    return;
}
template < class ELT >
void FactorLLT::getL( Matrix_<ELT>& m) const {
    rep->getL( m );
    return;
}

//////////////////
// FactorLLTRep 
//////////////////

template <typename T >
    template < typename ELT >
FactorLLTRep<T>::FactorLLTRep( const Matrix_<ELT>& mat ) 
      : nRow( mat.nrow() ),
        nCol( mat.ncol() ),
        mn( (mat.nrow() < mat.ncol()) ? mat.nrow() : mat.ncol() ),
        singularIndex(0),          
        lu( mat.nrow()*mat.ncol() )
{ 
    FactorLLTRep<T>::factor( mat );
}
template <typename T >
FactorLLTRep<T>::FactorLLTRep() 
      : nRow(0),
        nCol(0),
        mn(0),
        singularIndex(0),        
        lu(0)
{
}
template <typename T >
FactorLLTRep<T>::~FactorLLTRep() {}
 
template <typename T >
FactorLLTRepBase* FactorLLTRep<T>::clone() const {
   return( new FactorLLTRep<T>(*this) );
}
template < class T >
void FactorLLTRep<T>::inverse(  Matrix_<T>& inverse ) const {
    Matrix_<T> iden(mn,mn);
    iden.resize(mn,mn);
    iden = 1.0;
    solve( iden, inverse );
}

template < class T >
void FactorLLTRep<T>::solve( const Vector_<T>& b, Vector_<T> &x ) const {

    SimTK_APIARGCHECK2_ALWAYS(b.size()==nRow,"FactorLLT","solve",
       "number of rows in right hand side=%d does not match number of rows in original matrix=%d \n",
        b.size(), nRow );

    x.copyAssign(b);
    LapackInterface::potrs<T>( 'L', nCol, 1, lu.data, &x(0) );

    return;
}
template <typename T >
void FactorLLTRep<T>::solve(  const Matrix_<T>& b, Matrix_<T>& x ) const {

    SimTK_APIARGCHECK2_ALWAYS(b.nrow()==nRow,"FactorLLT","solve",
       "number of rows in right hand side=%d does not match number of rows in original matrix=%d \n",
        b.nrow(), nRow );

    x.copyAssign(b);
    LapackInterface::potrs<T>( 'L', nCol, b.ncol(), lu.data, &x(0,0) );
    return;
}
template <class T> 
    template<typename ELT>
void FactorLLTRep<T>::factor(const Matrix_<ELT>&mat )  {

    SimTK_APIARGCHECK2_ALWAYS(mat.nelt() > 0,"FactorLLT","factor",
       "Can't factor a matrix that has a zero dimension -- got %d X %d.",
       (int)mat.nrow(), (int)mat.ncol());
    
    // initialize the matrix we pass to LAPACK
    // converts (negated,conjugated etc.) to LAPACK format 
    LapackConvert::convertMatrixToLapack( lu.data, mat );


    int lda = nRow;
    int info;

    LapackInterface::potrf<T>('L', nCol, lu.data, lda, info);
    if( info > 0 ) 
        singularIndex = info; // matrix is singular info = i when U(i,i) is exactly zero
    else 
        singularIndex = 0;

}

template <typename T >
void FactorLLTRep<T>::getL( Matrix_<T>& m) const {
    int i,j;
      
    m.resize( nRow, nCol ); 

    for (int i = 0; i < nRow; ++i) {
        for (int j = 0; j <= i; ++j) {
            m(i,j) = lu.data[j*nRow + i]; // lower-triangle 
        }
        for (int j = i+1; j < nCol; ++j) {
            m(i,j) = 0; // upper-triangle is zero
        }
    }
}


/////////////
// Explicit instantiations
/////////////

// instantiate
template SimTK_SIMMATH_EXPORT FactorLLT::FactorLLT( const Matrix_<double>& m );
template SimTK_SIMMATH_EXPORT FactorLLT::FactorLLT( const Matrix_<float>& m );
template SimTK_SIMMATH_EXPORT FactorLLT::FactorLLT( const Matrix_<std::complex<float> >& m );
template SimTK_SIMMATH_EXPORT FactorLLT::FactorLLT( const Matrix_<std::complex<double> >& m );
template SimTK_SIMMATH_EXPORT FactorLLT::FactorLLT( const Matrix_<conjugate<float> >& m );
template SimTK_SIMMATH_EXPORT FactorLLT::FactorLLT( const Matrix_<conjugate<double> >& m );
template SimTK_SIMMATH_EXPORT FactorLLT::FactorLLT( const Matrix_<negator< double> >& m );
template SimTK_SIMMATH_EXPORT FactorLLT::FactorLLT( const Matrix_<negator< float> >& m );
template SimTK_SIMMATH_EXPORT FactorLLT::FactorLLT( const Matrix_<negator< std::complex<float> > >& m );
template SimTK_SIMMATH_EXPORT FactorLLT::FactorLLT( const Matrix_<negator< std::complex<double> > >& m );
template SimTK_SIMMATH_EXPORT FactorLLT::FactorLLT( const Matrix_<negator< conjugate<float> > >& m );
template SimTK_SIMMATH_EXPORT FactorLLT::FactorLLT( const Matrix_<negator< conjugate<double> > >& m );

template SimTK_SIMMATH_EXPORT void FactorLLT::factor( const Matrix_<double>& m );
template SimTK_SIMMATH_EXPORT void FactorLLT::factor( const Matrix_<float>& m );
template SimTK_SIMMATH_EXPORT void FactorLLT::factor( const Matrix_<std::complex<float> >& m );
template SimTK_SIMMATH_EXPORT void FactorLLT::factor( const Matrix_<std::complex<double> >& m );
template SimTK_SIMMATH_EXPORT void FactorLLT::factor( const Matrix_<conjugate<float> >& m );
template SimTK_SIMMATH_EXPORT void FactorLLT::factor( const Matrix_<conjugate<double> >& m );
template SimTK_SIMMATH_EXPORT void FactorLLT::factor( const Matrix_<negator< double> >& m );
template SimTK_SIMMATH_EXPORT void FactorLLT::factor( const Matrix_<negator< float> >& m );
template SimTK_SIMMATH_EXPORT void FactorLLT::factor( const Matrix_<negator< std::complex<float> > >& m );
template SimTK_SIMMATH_EXPORT void FactorLLT::factor( const Matrix_<negator< std::complex<double> > >& m );
template SimTK_SIMMATH_EXPORT void FactorLLT::factor( const Matrix_<negator< conjugate<float> > >& m );
template SimTK_SIMMATH_EXPORT void FactorLLT::factor( const Matrix_<negator< conjugate<double> > >& m );

template class FactorLLTRep<double>;
template FactorLLTRep<double>::FactorLLTRep( const Matrix_<double>& m);
template FactorLLTRep<double>::FactorLLTRep( const Matrix_<negator<double> >& m);
template void FactorLLTRep<double>::factor( const Matrix_<double>& m);
template void FactorLLTRep<double>::factor( const Matrix_<negator<double> >& m);

template class FactorLLTRep<float>;
template FactorLLTRep<float>::FactorLLTRep( const Matrix_<float>& m);
template FactorLLTRep<float>::FactorLLTRep( const Matrix_<negator<float> >& m);
template void FactorLLTRep<float>::factor( const Matrix_<float>& m);
template void FactorLLTRep<float>::factor( const Matrix_<negator<float> >& m);

template class FactorLLTRep<std::complex<double> >;
template FactorLLTRep<std::complex<double> >::FactorLLTRep( const Matrix_<std::complex<double> >& m);
template FactorLLTRep<std::complex<double> >::FactorLLTRep( const Matrix_<negator<std::complex<double> > >& m);
template FactorLLTRep<std::complex<double> >::FactorLLTRep( const Matrix_<conjugate<double> >& m);
template FactorLLTRep<std::complex<double> >::FactorLLTRep( const Matrix_<negator<conjugate<double> > >& m);
template void FactorLLTRep<std::complex<double> >::factor( const Matrix_<std::complex<double> >& m);
template void FactorLLTRep<std::complex<double> >::factor( const Matrix_<negator<std::complex<double> > >& m);
template void FactorLLTRep<std::complex<double> >::factor( const Matrix_<conjugate<double> >& m);
template void FactorLLTRep<std::complex<double> >::factor( const Matrix_<negator<conjugate<double> > >& m);

template class FactorLLTRep<std::complex<float> >;
template FactorLLTRep<std::complex<float> >::FactorLLTRep( const Matrix_<std::complex<float> >& m);
template FactorLLTRep<std::complex<float> >::FactorLLTRep( const Matrix_<negator<std::complex<float> > >& m);
template FactorLLTRep<std::complex<float> >::FactorLLTRep( const Matrix_<conjugate<float> >& m);
template FactorLLTRep<std::complex<float> >::FactorLLTRep( const Matrix_<negator<conjugate<float> > >& m);
template void FactorLLTRep<std::complex<float> >::factor( const Matrix_<std::complex<float> >& m);
template void FactorLLTRep<std::complex<float> >::factor( const Matrix_<negator<std::complex<float> > >& m);
template void FactorLLTRep<std::complex<float> >::factor( const Matrix_<conjugate<float> >& m);
template void FactorLLTRep<std::complex<float> >::factor( const Matrix_<negator<conjugate<float> > >& m);

template SimTK_SIMMATH_EXPORT void FactorLLT::getL<float>(Matrix_<float>&) const;
template SimTK_SIMMATH_EXPORT void FactorLLT::getL<double>(Matrix_<double>&) const;
template SimTK_SIMMATH_EXPORT void FactorLLT::getL<std::complex<float> >(Matrix_<std::complex<float> >&) const;
template SimTK_SIMMATH_EXPORT void FactorLLT::getL<std::complex<double> >(Matrix_<std::complex<double> >&) const;
template SimTK_SIMMATH_EXPORT void FactorLLT::solve<float>(const Vector_<float>&, Vector_<float>&) const;
template SimTK_SIMMATH_EXPORT void FactorLLT::solve<double>(const Vector_<double>&, Vector_<double>&) const;
template SimTK_SIMMATH_EXPORT void FactorLLT::solve<std::complex<float> >(const Vector_<std::complex<float> >&, Vector_<std::complex<float> >&) const;
template SimTK_SIMMATH_EXPORT void FactorLLT::solve<std::complex<double> >(const Vector_<std::complex<double> >&, Vector_<std::complex<double> >&) const;
template SimTK_SIMMATH_EXPORT void FactorLLT::solve<float>(const Matrix_<float>&, Matrix_<float>&) const;
template SimTK_SIMMATH_EXPORT void FactorLLT::solve<double>(const Matrix_<double>&, Matrix_<double>&) const;
template SimTK_SIMMATH_EXPORT void FactorLLT::solve<std::complex<float> >(const Matrix_<std::complex<float> >&, Matrix_<std::complex<float> >&) const;
template SimTK_SIMMATH_EXPORT void FactorLLT::solve<std::complex<double> >(const Matrix_<std::complex<double> >&, Matrix_<std::complex<double> >&) const;
template SimTK_SIMMATH_EXPORT void FactorLLT::inverse<float>(Matrix_<float>&) const;
template SimTK_SIMMATH_EXPORT void FactorLLT::inverse<double>(Matrix_<double>&) const;
template SimTK_SIMMATH_EXPORT void FactorLLT::inverse<std::complex<float> >(Matrix_<std::complex<float> >&) const;
template SimTK_SIMMATH_EXPORT void FactorLLT::inverse<std::complex<double> >(Matrix_<std::complex<double> >&) const;


} // namespace SimTK
