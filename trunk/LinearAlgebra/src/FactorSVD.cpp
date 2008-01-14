
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
 * Solves for singular values and singular vectors 
 */


#include <iostream> 
#include <cstdio>
#include <malloc.h>
#include <math.h>
#include <complex>
#include "SimTKcommon.h"
#include "LapackInterface.h"
#include "LinearAlgebra.h"
#include "FactorSVDRep.h"
#include "WorkSpace.h"
#include "LATraits.h"
#include "LapackConvert.h"


namespace SimTK {

   //////////////////
   // FactorSVDDefault //
   //////////////////
FactorSVDDefault::FactorSVDDefault() {
    haveMatrix = false;
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

template < class ELT >
FactorSVD::FactorSVD( const Matrix_<ELT>& m ) {
    rep = new FactorSVDRep<typename CNT<ELT>::StdNumber>(m); 
}
template < class ELT >
void FactorSVD::factor( const Matrix_<ELT>& m ) {
    delete rep;
    rep = new FactorSVDRep<typename CNT<ELT>::StdNumber>(m); 
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
FactorSVDRep<T>::FactorSVDRep( const Matrix_<ELT>& mat):
    n(mat.ncol()),  
    m(mat.nrow()),
    mn( (mat.nrow() < mat.ncol()) ? mat.nrow() : mat.ncol() ), 
    inputMatrix(n*n),
    structure(mat.getMatrixStructure())  {
    
    LapackInterface::getMachineUnderflow( abstol );
    abstol *= 0.5;

    LapackConvert::convertMatrixToLapack( inputMatrix.data, mat );
    haveMatrix = true;
        
}

template <typename T >
FactorSVDRep<T>::~FactorSVDRep() {}

template < class T >
void FactorSVDRep<T>::getSingularValuesAndVectors( Vector_<RType>& values,  Matrix_<T>& leftVectors, Matrix_<T>& rightVectors ) {
   

    leftVectors.resize(m,m);
    rightVectors.resize(n,n);
    values.resize(mn);

    computeSVD( true, &values(0), &leftVectors(0,0), &rightVectors(0,0) );

    return;
}
template < class T >
void FactorSVDRep<T>::getSingularValues( Vector_<RType>& values ) {

    values.resize(mn);

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

    LapackInterface::gesdd<T>(jobz, m,n,inputMatrix.data, m, values,
           leftVectors, m, rightVectors, n, info);

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

template SimTK_SIMMATH_EXPORT void FactorSVD::getSingularValues<float >(Vector_<float>&  );
template SimTK_SIMMATH_EXPORT void FactorSVD::getSingularValues<double >(Vector_<double>&  );

template SimTK_SIMMATH_EXPORT void FactorSVD::getSingularValuesAndVectors<float >(Vector_<float>&, Matrix_<float>&, Matrix_<float>&  );
template SimTK_SIMMATH_EXPORT void FactorSVD::getSingularValuesAndVectors<double >(Vector_<double>&, Matrix_<double>&, Matrix_<double>&  );
template SimTK_SIMMATH_EXPORT void FactorSVD::getSingularValuesAndVectors<std::complex<float> >(Vector_<float>&, Matrix_<std::complex<float> >&, Matrix_<std::complex<float> >&  );
template SimTK_SIMMATH_EXPORT void FactorSVD::getSingularValuesAndVectors<std::complex<double> >(Vector_<double>&, Matrix_<std::complex<double> >&, Matrix_<std::complex<double> >&  );


template class FactorSVDRep<double>;
template FactorSVDRep<double>::FactorSVDRep( const Matrix_<double>& m);
template FactorSVDRep<double>::FactorSVDRep( const Matrix_<negator<double> >& m);

template class FactorSVDRep<float>;
template FactorSVDRep<float>::FactorSVDRep( const Matrix_<float>& m );
template FactorSVDRep<float>::FactorSVDRep( const Matrix_<negator<float> >& m );

template class FactorSVDRep<std::complex<double> >;
template FactorSVDRep<std::complex<double> >::FactorSVDRep( const Matrix_<std::complex<double> >& m );
template FactorSVDRep<std::complex<double> >::FactorSVDRep( const Matrix_<negator<std::complex<double> > >& m );
template FactorSVDRep<std::complex<double> >::FactorSVDRep( const Matrix_<conjugate<double> >& m );
template FactorSVDRep<std::complex<double> >::FactorSVDRep( const Matrix_<negator<conjugate<double> > >& m );

template class FactorSVDRep<std::complex<float> >;
template FactorSVDRep<std::complex<float> >::FactorSVDRep( const Matrix_<std::complex<float> >& m );
template FactorSVDRep<std::complex<float> >::FactorSVDRep( const Matrix_<negator<std::complex<float> > >& m );
template FactorSVDRep<std::complex<float> >::FactorSVDRep( const Matrix_<conjugate<float> >& m );
template FactorSVDRep<std::complex<float> >::FactorSVDRep( const Matrix_<negator<conjugate<float> > >& m );

} // namespace SimTK
