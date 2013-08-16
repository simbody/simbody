
#ifndef SimTK_SIMMATH_LAPACK_INTERFACE_H_
#define SimTK_SIMMATH_LAPACK_INTERFACE_H_

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
 * These is a templatized, C++ callable interface to LAPACK and BLAS.
 * Each method must be explicitly specialized for the supported precisions.
 */


#include "SimTKcommon.h"
#include "SimTKlapack.h"

namespace SimTK {

class LapackInterface { 
   
public:

static int getLWork( float* work);
static int getLWork( double* work);
static int getLWork( std::complex<float>* work);
static int getLWork( std::complex<double>* work);

template <class T> static
void gelss( int m, int n,  int mn, int nrhs, 
           T* a, int lda, T* b,  int ldb, typename CNT<T>::TReal* s,
           typename CNT<T>::TReal rcond, int& rank, int& info);

template <class T> static
void gesdd( char jobz, int m, int n, T* a, int lda, 
           typename CNT<T>::TReal* s, T* u, int ldu,
           T* vt, int ldvt, int& info);

template <class T> static
void geev(char jobvl, char jobvr, int n, T* a, int lda, 
    std::complex<typename CNT<T>::TReal>* values, 
    T* vl, int ldvl, std::complex<typename CNT<T>::TReal>* vr, 
    int ldvr, T* work, int lwork, int& info );

template <class T> static
void syevx( char jobz, char range, char uplo, int n, T* a, int lda, 
    typename CNT<T>::TReal vl, typename CNT<T>::TReal vu, int il, int iu, 
    typename CNT<T>::TReal abstol, int& nFound, typename CNT<T>::TReal* values, 
    T* vectors, int LDVectors, int* ifail, int& info );
                  

template <class T> static
void syev( char jobz,  char uplo, int n, T* a_eigenVectors, int lda,  
    typename CNT<T>::TReal* eigenValues, int& info );


/* solve system of linear equations using the LU factorization  */
template <class T> static 
void potrs( char uplo, const int ncol, const int nrhs, const T *lu,  T *b ); 

template <class T> static 
// TODO void sytrs( char uplo, const int ncol, const int nrhs, const T *lu, const int* pivots, T *b ); 
void sytrs( char uplo, const int ncol, const int nrhs, T *lu, int* pivots, T *b ); 


template <class T> static 
void getrs( char trans, const int ncol, const int nrhs, const T *lu, const int* pivots, T *b ); 

template <class T> static 
void getrf( const int m, const int n, T *a, const int lda, int* pivots, int& info );

template <class T> static 
void gttrf( const int m, const int n, T* dl, T* d, T* du, T* du2, int* pivots, int& info );

template <class T> static 
void gbtrf( const int m, const int n, const int kl, const int ku, T* lu, const int lda, int* pivots, int& info );

template <class T> static 
void potrf( const char& uplo, const int n,  T* lu, const int lda, int& info );

template <class T> static 
void sytrf( const char& uplo, const int n, T* a,  const int lda, int* pivots, T* work, const int lwork, int& info );

template <class T> static
int ilaenv( const int& ispec,  const char* name,  const char* opts, const int& n1, const int& n2, const int& n3, const int& n4  );

template <class T> static
void getMachinePrecision( T& smallNumber, T& bigNumber );

template <class T> static
void getMachineUnderflow( T& underFlow );

template <class T> static
void tzrzf( const int& m, const int& n,  T* a, const int& lda, T* tau, T* work, const int& lwork, int& info );

template <class T> static
void geqp3( const int& m, const int& n,  T* a, const int& lda, int *pivots, T* tau, T* work, const int& lwork, int& info );

template <class T> static
void lascl( const char& type, const int& kl, const int& ku, const typename CNT<T>::TReal& cfrom, const typename CNT<T>::TReal& cto,  const int& m, const int& n, T* a, const int& lda, int& info );

template <class T> static
double lange( const char& norm, const int& m, const int& n, const T* a, const int& lda );

template <class T> static
void ormqr(const char& side, const char& trans, const int& m, const int& n, const int& k, T* a, const int& lda, T *tau, T *c__, const int& ldc, T* work, const int& lwork, int& info);

template <class T> static
void trsm(const char& side, const char& uplo, const char& transA, const char& diag, const int& m, const int& n, const T& alpha, const T* A, const int& lda, T* B, const int& ldb );
 
template <class T> static
void ormrz(const char& side, const char& trans, const int& m, const int& n, const int& k, const int& l, T* a, const int& lda, T* tau, T* c__, const int& ldc, T* work, const int& lwork, int& info);
 
template <class T> static
void copy( const int& n, const T* x, const int& incx, T* y, const int& incy);

template <class T> static
void laic1(const int& job, const int& j, const T* x, const typename CNT<T>::TReal& sest, const T* w, const T& gamma, typename CNT<T>::TReal& sestpr, T& s, T& c__ );

}; // class LapackInterface

}   // namespace SimTK

#endif // SimTK_SIMMATH_LAPACK_INTERFACE_H_
