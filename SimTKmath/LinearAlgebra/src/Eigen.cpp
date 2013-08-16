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
 * Solves for eigen values and eigen vectors 
 */



#include "SimTKcommon.h"

#include "simmath/internal/common.h"
#include "simmath/LinearAlgebra.h"

#include "LapackInterface.h"
#include "EigenRep.h"
#include "WorkSpace.h"
#include "LATraits.h"
#include "LapackConvert.h"

#include <iostream> 
#include <cstdio>
#include <cmath>
#include <complex>

namespace SimTK {

   //////////////////
   // EigenDefault //
   //////////////////
EigenDefault::EigenDefault() {
    isFactored = false;
}
EigenRepBase* EigenDefault::clone() const {
    return( new EigenDefault(*this));
}


   ///////////
   // Eigen //
   ///////////
Eigen::~Eigen() {
    delete rep;
}
Eigen::Eigen() {
   rep = new EigenDefault();
}

template < class ELT >
Eigen::Eigen( const Matrix_<ELT>& m ) {
    rep = new EigenRep<typename CNT<ELT>::StdNumber>(m); 
}
// copy constructor
Eigen::Eigen( const Eigen& c ) {
    rep = c.rep->clone();
}
// copy assignment operator
Eigen& Eigen::operator=(const Eigen& rhs) {
    rep = rhs.rep->clone();
    return *this;
}

template < class ELT >
void Eigen::factor( const Matrix_<ELT>& m ) {
    delete rep;
    rep = new EigenRep<typename CNT<ELT>::StdNumber>(m); 
}

template <class VAL, class VEC> 
void Eigen::getAllEigenValuesAndVectors( Vector_<VAL>& values, Matrix_<VEC>& vectors) {
    rep->getAllEigenValuesAndVectors( values, vectors );
    return;
}
template <class VAL> 
void Eigen::getAllEigenValues( Vector_<VAL>& values) {
    rep->getAllEigenValues( values );
    return;
}
template <class VAL, class VEC> 
void Eigen::getFewEigenValuesAndVectors( Vector_<VAL>& values, Matrix_<VEC>& vectors, int ilow, int ihi) {
    rep->getFewEigenValuesAndVectors( values, vectors, ilow, ihi );
    return;
}
template <class T> 
void Eigen::getFewEigenVectors( Matrix_<T>& vectors, int ilow, int ihi ) {
    rep->getFewEigenVectors( vectors, ilow, ihi );
    return;
}

template <class T> 
void Eigen::getFewEigenValues( Vector_<T>& values, int ilow, int ihi ) {
    rep->getFewEigenValues( values, ilow, ihi );
    return;
}


template <class VAL, class VEC> 
void Eigen::getFewEigenValuesAndVectors( Vector_<VAL>& values, Matrix_<VEC>& vectors, typename CNT<VAL>::TReal rlow, typename CNT<VAL>::TReal rhi) {
    rep->getFewEigenValuesAndVectors( values, vectors, rlow, rhi );
    return;
}
template <class T> 
void Eigen::getFewEigenVectors( Matrix_<T>& vectors, typename CNT<T>::TReal rlow, typename CNT<T>::TReal rhi) {
    rep->getFewEigenVectors( vectors, rlow, rhi );
    return;
}

template <class T> 
void Eigen::getFewEigenValues( Vector_<T>& values, typename CNT<T>::TReal rlow, typename CNT<T>::TReal rhi) {
    rep->getFewEigenValues( values, rlow, rhi );
    return;
}

   /////////////////
   // EigenRep //
   /////////////////
template <typename T >        // constructor 
    template < typename ELT >
EigenRep<T>::EigenRep( const Matrix_<ELT>& mat):
    n(mat.ncol()),  
    inputMatrix(n*n),
    structure(mat.getMatrixCharacter().getStructure()),
    range(AllValues),
    valuesFound(0),
    vectorsInMatrix(false),
    needValues(true),
    needVectors(true) {  
         
    LapackInterface::getMachineUnderflow( abstol );
    abstol *= 0.5;

    LapackConvert::convertMatrixToLapack( inputMatrix.data, mat );
    isFactored = true;
        
}
template <typename T >
EigenRepBase* EigenRep<T>::clone() const {
   return( new EigenRep<T>(*this) );
}



template <typename T >
EigenRep<T>::~EigenRep() {}

template <typename T >
void EigenRep<T>::copyValues(Vector_<float>& values) {

    values.resize(valuesFound);
    for(int j = 0;j<valuesFound;j++ ) {
        values(j) = (float)realEigenValues.data[j];
    }
    return;
}
template <typename T >
void EigenRep<T>::copyValues(Vector_<double>& values) {

    values.resize(valuesFound);
    for(int j = 0;j<valuesFound;j++ ) {
        values(j) = realEigenValues.data[j];
    }
    return;
}
template <typename T >
void EigenRep<T>::copyValues(Vector_<std::complex<float> >& values) {

    values.resize(valuesFound);
    for(int j = 0;j<valuesFound;j++ ) {
        values(j) = complexEigenValues.data[j];
    }
}
template <typename T >
void EigenRep<T>::copyValues(Vector_<std::complex<double> >& values) {

    values.resize(valuesFound);
    for(int j = 0;j<valuesFound;j++ ) {
        values(j) = complexEigenValues.data[j];
    }
}
template <>
    template <>
void EigenRep<float>::copyVectors(Matrix_<float>& vectors) {   // symmetric
    int i,j;

    vectors.resize(n,valuesFound);
    if( vectorsInMatrix ) {
         for(j=0;j<valuesFound;j++) for(i=0;i<n;i++) vectors(i,j) = inputMatrix.data[j*n+i];
    } else {
         for(j=0;j<valuesFound;j++) for(i=0;i<n;i++) vectors(i,j) = symmetricEigenVectors.data[j*n+i];
    }

    return;
}
template <>
    template <>
void EigenRep<double>::copyVectors(Matrix_<double>& vectors) { // symmetric
    int i,j;

    vectors.resize(n,valuesFound);
    if( vectorsInMatrix ) {
        for(j=0;j<valuesFound;j++) for(i=0;i<n;i++) vectors(i,j) = inputMatrix.data[j*n+i];
    } else {
        for(j=0;j<valuesFound;j++) for(i=0;i<n;i++) vectors(i,j) = symmetricEigenVectors.data[j*n+i];
    }
    return;
}
template <>
    template <>
void EigenRep<float>::copyVectors(Matrix_<std::complex<float> >& vectors) { // non-symmetric
    int i,j;

    vectors.resize(n,valuesFound);
    for(j=0;j<valuesFound;j++) for(i=0;i<n;i++) vectors(i,j) = complexEigenVectors.data[j*n+i];

    return;
}
template <>
    template <>
void EigenRep<double>::copyVectors(Matrix_<std::complex<double> >& vectors) { // non-symmetric
    int i,j;

    vectors.resize(n,valuesFound);
    for(j=0;j<valuesFound;j++) for(i=0;i<n;i++) vectors(i,j) = complexEigenVectors.data[j*n+i];

    return;
}
template <>
   template <>
void EigenRep<std::complex<float> >::copyVectors(Matrix_<std::complex<float> >& vectors) {
    int i,j;

    vectors.resize(n,valuesFound);
/*
    if( structure ==  MatrixStructure::Symmetric ) {
        if( vectorsInMatrix ) {
            for(j=0;j<valuesFound;j++) for(i=0;i<n;i++) vectors(i,j) = inputMatrix.data[j*n+i];
        } else {
            for(j=0;j<valuesFound;j++) for(i=0;i<n;i++) vectors(i,j) = symmetricEigenVectors.data[j*n+i];
        }
    } else {
*/
        for(j=0;j<valuesFound;j++) for(i=0;i<n;i++) vectors(i,j) = complexEigenVectors.data[j*n+i];
//   }
    return;
}
template <>
   template <>
void EigenRep<std::complex<double> >::copyVectors(Matrix_<std::complex<double> >& vectors) {
    int i,j;
/*
    vectors.resize(n,valuesFound);
    if( structure ==  MatrixStructure::Symmetric ) {
        if( vectorsInMatrix ) {
            for(j=0;j<valuesFound;j++) for(i=0;i<n;i++) vectors(i,j) = inputMatrix.data[j*n+i];
        } else {
            for(j=0;j<valuesFound;j++) for(i=0;i<n;i++) vectors(i,j) = symmetricEigenVectors.data[j*n+i];
        }
    } else {
*/
        for(j=0;j<valuesFound;j++) for(i=0;i<n;i++) vectors(i,j) = complexEigenVectors.data[j*n+i];
 //   }
    return;
}
template <>
   template <>
void EigenRep<std::complex<double> >::copyVectors(Matrix_<double>& vectors) { assert(false);} // should never get called

template <>
   template <>
void EigenRep<std::complex<float> >::copyVectors(Matrix_<float>& vectors) { assert(false);} // should never get called

template < class T >
void EigenRep<T>::computeValues(bool computeVectors) {

    int info;

    T *leftVectors = NULL;           // no left vectors are calculated
    char calcLeftEigenVectors  = 'N';     // don't compute left eigen vectors
    char calcRightEigenVectors; 
    if( computeVectors )  {
        calcRightEigenVectors = 'V'; // compute  right eigen vectors
    } else {
        calcRightEigenVectors = 'N'; // don't compute right eigen vectors
    } 
    int computeLwork = -1;
    //T size[1];
/*
    if( structure ==  MatrixStructure::Symmetric ) {
         char useUpper = 'U'; // matrix is stored in upper triangle
         realEigenValues.resize(n);
         if( range == AllValues ) {

             LapackInterface::syev<T>( calcRightEigenVectors, useUpper, n,
                 inputMatrix.data, n, realEigenValues.data, info );
             if( info < 0 ) {
                 SimTK_THROW2( SimTK::Exception::ConvergedFailed,
                 "LapackInterface::syev",
                 "off-diagonal elements of intermediate tridiagonal" );
             }
             vectorsInMatrix = true;
             valuesFound = n;

         } else {
             char rangeChar; 
             if( range == IndexRange ) {
                 rangeChar = 'I';
             } else {
                 rangeChar = 'V';
             }
             LapackInterface::getMachineUnderflow( abstol );
             abstol *= 0.5;
             ifail.resize((long)n);
             symmetricEigenVectors.resize(n*n);
             LapackInterface::syevx<T>( calcRightEigenVectors, rangeChar, 
             useUpper, n, inputMatrix.data, n, lowValue, hiValue, 
             lowIndex, hiIndex, abstol, valuesFound, realEigenValues.data,
             symmetricEigenVectors.data, n, ifail.data, info );

             if( info < 0 ) {
                 SimTK_THROW2( SimTK::Exception::ConvergedFailed,
                 "LapackInterface::syev",
                 "off-diagonal elements of intermediate tridiagonal" );
             }
         }

    } else {
*/
          complexEigenVectors.resize(n*n);
          complexEigenValues.resize(n);
/*
          LapackInterface::geev<T>( calcLeftEigenVectors, calcRightEigenVectors,
               n, inputMatrix.data, n, complexEigenValues.data, leftVectors, 
               n, complexEigenVectors.data, n, size, computeLwork, info);
          if( info < 0 ) {
               SimTK_THROW2( SimTK::Exception::ConvergedFailed,
               "LapackInterface::geev",
               "QR algorithm" );
          }


          int lwork = LapackInterface::getLWork( size );
*/
          int lwork = 2000;
          valuesFound = n;
          TypedWorkSpace<T> workSpace(lwork);
    
          // compute all eigen values and eigen vectors of a general matrix

          LapackInterface::geev<T>( calcLeftEigenVectors, calcRightEigenVectors,
                              n, inputMatrix.data, n, complexEigenValues.data,
                              leftVectors, n, complexEigenVectors.data, n, 
                              workSpace.data,  lwork, info);
          if( info < 0 ) {
               SimTK_THROW2( SimTK::Exception::ConvergedFailed,
               "LapackInterface::geev",
               "QR algorithm" );
          }

 //   } 
    // copy computed eigen values and eigen vectors into caller's arguements
    if( info != 0 ) {
         // TODO THROW EXCEPTION
    } 
    needVectors = false;
    needValues = false;
    return;
}
template < class T >
void EigenRep<T>::getAllEigenValuesAndVectors( Vector_<std::complex<RType> >& values,  Matrix_<std::complex<RType> >& vectors ) {

    range = AllValues;

    if( needValues ) computeValues( true );
    copyValues( values );
    copyVectors( vectors );

    return;
}
// only for symmetric real matrix
template < class T >
void EigenRep<T>::getAllEigenValuesAndVectors( Vector_<RType>& values,  Matrix_<RType>& vectors ) {

    // symmtric matrices return real eigen values,  nonsymmetric matrices return complex eigen values
    if( !CNT<T>::IsPrecision  ) {
         SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getAllEigenValues",
         "getAllEigenValuesAndVectors(Real, Real) called with for a complex matrix   \n");
    }
    if( structure.getStructure() != MatrixStructure::Symmetric ) {
         SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getAllEigenValues",
         "getAllEigenValuesAndVectors(Real, Real) called with   for a non symmetric matrix   \n");
    }

    range = AllValues;

    if( needValues ) computeValues( true );
    copyValues( values );
    copyVectors( vectors );

    return;
}
// only for symmetric complex matrix
template < class T >
void EigenRep<T>::getAllEigenValuesAndVectors( Vector_<RType>& values,  Matrix_<std::complex<RType> >& vectors ) {

    // symmtric matrices return real eigen values,  nonsymmetric matrices return complex eigen values
    if( CNT<T>::IsPrecision  ) {
         SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getAllEigenValues",
         "getAllEigenValuesAndVectors(Real, complex) called for a real matrix   \n");
    }
    if( structure.getStructure() != MatrixStructure::Symmetric ) {
         SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getAllEigenValues",
         "getAllEigenValuesAndVectors(Real, complex) called for a non symmetric matrix   \n");
    }

    range = AllValues;

    if( needValues ) computeValues( true );
    copyValues( values );
    copyVectors( vectors );

    return;
}

template < class T >
void EigenRep<T>::getAllEigenValues( Vector_<RType>& values ) {

    // symmtric matrices return real eigen values,  nonsymmetric matrices return complex eigen values
    if(  structure.getStructure() != MatrixStructure::Symmetric ) {
         SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getAllEigenValues",
         "getAllEigenValues called with value of type real for a non symmetric matrix   \n");
    }
    range = AllValues;
    if( needValues ) computeValues( false ); 
    copyValues( values );

    return; 
}
template < class T >
void EigenRep<T>::getAllEigenValues( Vector_<std::complex<RType> >& values ) {

    range = AllValues;
    if( needValues ) computeValues( false ); 
    copyValues( values );

    return; 
}
template < class T >
void EigenRep<T>::getFewEigenValuesAndVectors( Vector_<RType>& values, Matrix_<RType>& vectors,  typename CNT<T>::TReal rlow, typename CNT<T>::TReal rhi ) {

    // symmtric matrices return real eigen values,  nonsymmetric matrices return complex eigen values
    if(  structure.getStructure() != MatrixStructure::Symmetric ) {
         SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenValuesAndVectors",
         "getFewEigenValuesAndVectors(real, real ) called for a non symmetric matrix   \n");
    }
    range = ValueRange;
    lowValue = rlow;
    hiValue = rhi;
    if( needValues ) computeValues(true);

    copyValues(values);
    copyVectors( vectors );
    

    return; 
}
template < class T >
void EigenRep<T>::getFewEigenValuesAndVectors( Vector_<RType>& values, Matrix_<std::complex<RType> >& vectors,  typename CNT<T>::TReal rlow, typename CNT<T>::TReal rhi ) {

    // symmtric matrices return real eigen values,  nonsymmetric matrices return complex eigen values
    if(  structure.getStructure() != MatrixStructure::Symmetric ) {
         SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenValuesAndVectors",
         "getFewEigenValuesAndVectors(real, complex ) called for a non symmetric matrix   \n");
    }
    if(  !CNT<T>::IsPrecision ) {
         SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenValuesAndVectors",
         "getFewEigenValuesAndVectors(real, complex ) called for a real matrix   \n");
    }
    range = ValueRange;
    lowValue = rlow;
    hiValue = rhi;
    if( needValues ) computeValues(true);

    copyValues(values);
    copyVectors( vectors );
    

    return; 
}
template < class T >
void EigenRep<T>::getFewEigenValuesAndVectors( Vector_<std::complex<RType> >& values, Matrix_<std::complex<RType> >& vectors,  typename CNT<T>::TReal rlow, typename CNT<T>::TReal rhi ) {

    range = ValueRange;
    lowValue = rlow;
    hiValue = rhi;
    if( needValues ) computeValues(true);

    copyValues(values);
    copyVectors( vectors );
    

    return; 
}

template < class T >
void EigenRep<T>::getFewEigenValuesAndVectors( Vector_<RType>& values, Matrix_<RType>& vectors,  int ilow, int ihi ) {

    // symmtric matrices return real eigen values,  nonsymmetric matrices return complex eigen values
    if(  structure.getStructure() != MatrixStructure::Symmetric ) {
         SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenValuesAndVectors",
         "getFewEigenValuesAndVectors(real, complex ) called for a non symmetric matrix   \n");
    }
    if(  !CNT<T>::IsPrecision ) {
         SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenValuesAndVectors",
         "getFewEigenValuesAndVectors(real, real ) called for a complex matrix   \n");
    }
    range = IndexRange;
    lowIndex = ilow;
    hiIndex = ihi;
    if( needValues ) computeValues( true );

    copyValues(values);
    copyVectors( vectors );

    return; 
}
template < class T >
void EigenRep<T>::getFewEigenValuesAndVectors( Vector_<RType>& values, Matrix_<std::complex<RType> >& vectors,  int ilow, int ihi ) {

    // symmtric matrices return real eigen values,  nonsymmetric matrices return complex eigen values
    if(  structure.getStructure() != MatrixStructure::Symmetric ) {
         SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenValuesAndVectors",
         "getFewEigenValuesAndVectors(real, complex ) called for a non symmetric matrix   \n");
    }
    if( CNT<T>::IsPrecision ) {
         SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenValuesAndVectors",
         "getFewEigenValuesAndVectors(real, complex ) called for a real matrix   \n");
    }
    range = IndexRange;
    lowIndex = ilow;
    hiIndex = ihi;
    if( needValues ) computeValues( true );

    copyValues(values);
    copyVectors( vectors );

    return; 
}
template < class T >
void EigenRep<T>::getFewEigenValuesAndVectors( Vector_<std::complex<RType> >& values, Matrix_<std::complex<RType> >& vectors,  int ilow, int ihi ) {

    range = IndexRange;
    lowIndex = ilow;
    hiIndex = ihi;
    if( needValues ) computeValues( true );

    copyValues(values);
    copyVectors( vectors );

    return; 
}
template < class T >
void EigenRep<T>::getFewEigenValues( Vector_<RType>& values,  typename CNT<T>::TReal rlow, typename CNT<T>::TReal rhi ) {

    // symmetric matrices return real eigen values,  nonsymmetric matrices return complex eigen values
    if(  structure.getStructure() != MatrixStructure::Symmetric ) {
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenValues",
       "getFewEigenValues(real) called for a non symmetric matrix   \n");
    }
    range = ValueRange;
    lowValue = rlow;
    hiValue = rhi;
    if( needValues ) computeValues(false);

    copyValues(values);

    return; 
}
template < class T >
void EigenRep<T>::getFewEigenValues( Vector_<std::complex<RType> >& values,  typename CNT<T>::TReal rlow, typename CNT<T>::TReal rhi ) {

    range = ValueRange;
    lowValue = rlow;
    hiValue = rhi;
    if( needValues ) computeValues(false);

    copyValues(values);

    return; 
}
template < class T >
void EigenRep<T>::getFewEigenValues( Vector_<RType>& values,  int ilow, int ihi ) {

    // symmtric matrices return real eigen values,  nonsymmetric matrices return complex eigen values
    if( structure.getStructure() != MatrixStructure::Symmetric ) {
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenValues",
       "getFewEigenValues(real) called for a non symmetric matrix   \n");
    }
    range = IndexRange;
    lowIndex = ilow;
    hiIndex = ihi;
    if( needValues ) computeValues(false);

    copyValues(values);

    return; 
}
template < class T >
void EigenRep<T>::getFewEigenValues( Vector_<std::complex<RType> >& values,  int ilow, int ihi ) {

    range = IndexRange;
    lowIndex = ilow;
    hiIndex = ihi;
    if( needValues ) computeValues(false);

    copyValues(values);

    return; 
}
template < class T >
void EigenRep<T>::getFewEigenVectors( Matrix_<std::complex<RType> >& vectors,  int ilow, int ihi ) {

    range = IndexRange;
    lowIndex = ilow;
    hiIndex = ihi;
    if( needVectors ) computeValues( true );
    copyVectors( vectors );

    return; 
}
template < class T >
void EigenRep<T>::getFewEigenVectors( Matrix_<std::complex<RType> >& vectors,  typename CNT<T>::TReal rlow, typename CNT<T>::TReal rhi ) {

    range = ValueRange;
    lowValue = rlow;
    hiValue = rhi;
    if( needVectors ) computeValues( true );
    copyVectors( vectors );

    return; 
}
template < class T >
void EigenRep<T>::getFewEigenVectors( Matrix_<RType>& vectors,  int ilow, int ihi ) {

    // symmtric matrices return real eigen values,  nonsymmetric matrices return complex eigen values
    if( structure.getStructure() != MatrixStructure::Symmetric ) {
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenValues",
       "getFewEigenVectors(real) called for a non symmetric matrix   \n");
    } 
    if( !CNT<T>::IsPrecision ) {
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenValues",
       "getFewEigenVectors(real) called for a complex matrix   \n");
    } 
    range = IndexRange;
    lowIndex = ilow;
    hiIndex = ihi;
    if( needVectors ) computeValues( true );
    copyVectors( vectors );

    return; 
}
template < class T >
void EigenRep<T>::getFewEigenVectors( Matrix_<RType>& vectors,  RType rlow, RType rhi ) {

    // symmtric matrices return real eigen values,  nonsymmetric matrices return complex eigen values
    if( structure.getStructure() != MatrixStructure::Symmetric ) {
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenValues",
       "getFewEigenVectors(real) called for a non symmetric matrix   \n");
    } 
    if( !CNT<T>::IsPrecision ) {
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenValues",
       "getFewEigenVectors(real) called for a complex matrix   \n");
    } 
    range = ValueRange;
    lowValue = rlow;
    hiValue = rhi;
    if( needVectors ) computeValues( true );
    copyVectors( vectors );

    return; 
}

// instantiate
template SimTK_SIMMATH_EXPORT Eigen::Eigen( const Matrix_<double>& m );
template SimTK_SIMMATH_EXPORT Eigen::Eigen( const Matrix_<float>& m );
template SimTK_SIMMATH_EXPORT Eigen::Eigen( const Matrix_<std::complex<float> >& m );
template SimTK_SIMMATH_EXPORT Eigen::Eigen( const Matrix_<std::complex<double> >& m );
template SimTK_SIMMATH_EXPORT Eigen::Eigen( const Matrix_<conjugate<float> >& m );
template SimTK_SIMMATH_EXPORT Eigen::Eigen( const Matrix_<conjugate<double> >& m );
template SimTK_SIMMATH_EXPORT Eigen::Eigen( const Matrix_<negator< double> >& m );
template SimTK_SIMMATH_EXPORT Eigen::Eigen( const Matrix_<negator< float> >& m );
template SimTK_SIMMATH_EXPORT Eigen::Eigen( const Matrix_<negator< std::complex<float> > >& m );
template SimTK_SIMMATH_EXPORT Eigen::Eigen( const Matrix_<negator< std::complex<double> > >& m );
template SimTK_SIMMATH_EXPORT Eigen::Eigen( const Matrix_<negator< conjugate<float> > >& m );
template SimTK_SIMMATH_EXPORT Eigen::Eigen( const Matrix_<negator< conjugate<double> > >& m );

//template SimTK_SIMMATH_EXPORT void Eigen::getAllEigenValuesAndVectors<float, float >(Vector_<float>&, Matrix_<float>& );
//template SimTK_SIMMATH_EXPORT void Eigen::getAllEigenValuesAndVectors<double, double >(Vector_<double>&, Matrix_<double>& );
//template SimTK_SIMMATH_EXPORT void Eigen::getAllEigenValuesAndVectors<float, std::complex<float> >(Vector_<float>&, Matrix_<std::complex<float> >& );
//template SimTK_SIMMATH_EXPORT void Eigen::getAllEigenValuesAndVectors<double, std::complex<double> >(Vector_<double>&, Matrix_<std::complex<double> >& );
template SimTK_SIMMATH_EXPORT void Eigen::getAllEigenValuesAndVectors<std::complex<float>, std::complex<float> >(Vector_<std::complex<float> >&, Matrix_<std::complex<float> >& );
template SimTK_SIMMATH_EXPORT void Eigen::getAllEigenValuesAndVectors<std::complex<double>, std::complex<double> >(Vector_<std::complex<double> >&, Matrix_<std::complex<double> >& );
/*
template SimTK_SIMMATH_EXPORT void Eigen::getFewEigenValuesAndVectors<float, float >(Vector_<float>&, Matrix_<float>&, int ilow, int ihi );
template SimTK_SIMMATH_EXPORT void Eigen::getFewEigenValuesAndVectors<double, double >(Vector_<double>&, Matrix_<double>&, int ilow, int ihi );
template SimTK_SIMMATH_EXPORT void Eigen::getFewEigenValuesAndVectors<float, std::complex<float> >(Vector_<float>&, Matrix_<std::complex<float> >&, int ilow, int ihi );
template SimTK_SIMMATH_EXPORT void Eigen::getFewEigenValuesAndVectors<double, std::complex<double> >(Vector_<double>&, Matrix_<std::complex<double> >&, int ilow, int ihi );
template SimTK_SIMMATH_EXPORT void Eigen::getFewEigenValuesAndVectors<std::complex<float>, std::complex<float> >(Vector_<std::complex<float> >&, Matrix_<std::complex<float> >&, int ilow, int ihi );
template SimTK_SIMMATH_EXPORT void Eigen::getFewEigenValuesAndVectors<std::complex<double>, std::complex<double> >(Vector_<std::complex<double> >&, Matrix_<std::complex<double> >&, int ilow, int ihi );

template SimTK_SIMMATH_EXPORT void Eigen::getFewEigenValuesAndVectors<float, float >(Vector_<float>&, Matrix_<float>&, float rlow, float rhi );
template SimTK_SIMMATH_EXPORT void Eigen::getFewEigenValuesAndVectors<double, double >(Vector_<double>&, Matrix_<double>&, double rlow, double rhi );
template SimTK_SIMMATH_EXPORT void Eigen::getFewEigenValuesAndVectors<float, std::complex<float> >(Vector_<float>&, Matrix_<std::complex<float> >&, float rlow, float rhi );
template SimTK_SIMMATH_EXPORT void Eigen::getFewEigenValuesAndVectors<double, std::complex<double> >(Vector_<double>&, Matrix_<std::complex<double> >&, double rlow, double rhi );
template SimTK_SIMMATH_EXPORT void Eigen::getFewEigenValuesAndVectors<std::complex<float>, std::complex<float> >(Vector_<std::complex<float> >&, Matrix_<std::complex<float> >&, float rlow, float rhi );
template SimTK_SIMMATH_EXPORT void Eigen::getFewEigenValuesAndVectors<std::complex<double>, std::complex<double> >(Vector_<std::complex<double> >&, Matrix_<std::complex<double> >&, double rlow, double rhi );
*/

//template SimTK_SIMMATH_EXPORT void Eigen::getAllEigenValues<double>(Vector_<double>& );
//template SimTK_SIMMATH_EXPORT void Eigen::getAllEigenValues<float>(Vector_<float>& );
template SimTK_SIMMATH_EXPORT void Eigen::getAllEigenValues<std::complex<double> >(Vector_<std::complex<double> >& );
template SimTK_SIMMATH_EXPORT void Eigen::getAllEigenValues<std::complex<float> >(Vector_<std::complex<float> >& );
/*
template SimTK_SIMMATH_EXPORT void Eigen::getFewEigenValues<double>(Vector_<double>&, int ilow, int ihi);
template SimTK_SIMMATH_EXPORT void Eigen::getFewEigenValues<float>(Vector_<float>&, int ilow, int ihi );
template SimTK_SIMMATH_EXPORT void Eigen::getFewEigenValues<std::complex<double> >(Vector_<std::complex<double> >&, int ilow, int ihi );
template SimTK_SIMMATH_EXPORT void Eigen::getFewEigenValues<std::complex<float> >(Vector_<std::complex<float> >&, int ilow, int ihi );
template SimTK_SIMMATH_EXPORT void Eigen::getFewEigenValues<double>(Vector_<double>&, double rlow, double rhi);
template SimTK_SIMMATH_EXPORT void Eigen::getFewEigenValues<float>(Vector_<float>&, float rlow, float rhi );
template SimTK_SIMMATH_EXPORT void Eigen::getFewEigenValues<std::complex<double> >(Vector_<std::complex<double> >&, double rlow, double rhi );
template SimTK_SIMMATH_EXPORT void Eigen::getFewEigenValues<std::complex<float> >(Vector_<std::complex<float> >&, float rlow, float rhi );

template SimTK_SIMMATH_EXPORT void Eigen::getFewEigenVectors<double>(Matrix_<double>&, int ilow, int ihi);
template SimTK_SIMMATH_EXPORT void Eigen::getFewEigenVectors<float>(Matrix_<float>&, int ilow, int ihi );
template SimTK_SIMMATH_EXPORT void Eigen::getFewEigenVectors<std::complex<double> >(Matrix_<std::complex<double> >&, int ilow, int ihi );
template SimTK_SIMMATH_EXPORT void Eigen::getFewEigenVectors<std::complex<float> >(Matrix_<std::complex<float> >&, int ilow, int ihi );
template SimTK_SIMMATH_EXPORT void Eigen::getFewEigenVectors<double>(Matrix_<double>&, double rlow, double rhi);
template SimTK_SIMMATH_EXPORT void Eigen::getFewEigenVectors<float>(Matrix_<float>&, float rlow, float rhi );
template SimTK_SIMMATH_EXPORT void Eigen::getFewEigenVectors<std::complex<double> >(Matrix_<std::complex<double> >&, double rlow, double rhi );
*/

template class EigenRep<double>;
template EigenRep<double>::EigenRep( const Matrix_<double>& m);
template EigenRep<double>::EigenRep( const Matrix_<negator<double> >& m);

template class EigenRep<float>;
template EigenRep<float>::EigenRep( const Matrix_<float>& m );
template EigenRep<float>::EigenRep( const Matrix_<negator<float> >& m );

template class EigenRep<std::complex<double> >;
template EigenRep<std::complex<double> >::EigenRep( const Matrix_<std::complex<double> >& m );
template EigenRep<std::complex<double> >::EigenRep( const Matrix_<negator<std::complex<double> > >& m );
template EigenRep<std::complex<double> >::EigenRep( const Matrix_<conjugate<double> >& m );
template EigenRep<std::complex<double> >::EigenRep( const Matrix_<negator<conjugate<double> > >& m );

template class EigenRep<std::complex<float> >;
template EigenRep<std::complex<float> >::EigenRep( const Matrix_<std::complex<float> >& m );
template EigenRep<std::complex<float> >::EigenRep( const Matrix_<negator<std::complex<float> > >& m );
template EigenRep<std::complex<float> >::EigenRep( const Matrix_<conjugate<float> >& m );
template EigenRep<std::complex<float> >::EigenRep( const Matrix_<negator<conjugate<float> > >& m );

} // namespace SimTK
