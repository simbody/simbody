
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
 * These is a templatized, C++ callable interface to LAPACK and BLAS.
 * Each method must be explicitly specialized for the supported precisions.
 */


#include "SimTKcommon.h"
#include "SimTKlapack.h"
#include "LapackInterface.h"

namespace SimTK {


template <> void LapackInterface::getrs<double>
    ( const bool transpose, const int ncol, const int nrhs, const double *lu, const int *pivots, double *b ) {

    int info;
    char trans;
    
    if( transpose ) 
        trans = 'T';
    else
        trans = 'N';

    dgetrs_(trans, ncol, nrhs, lu, ncol, pivots, b, ncol, info, 1  );

    return;
}


template <> void LapackInterface::getrs<float>
    ( const bool transpose, const int ncol, const int nrhs, const float *lu, const int *pivots, float *b ) {

    int info;
    char trans;

    if( transpose ) 
        trans = 'T';
    else
        trans = 'N';

    sgetrs_(trans, ncol, nrhs, lu, ncol, pivots, b, ncol, info, 1  );

    return;
}

template <> void LapackInterface::getrs<complex<float> >
    ( const bool transpose, const int ncol, const int nrhs, const std::complex<float> *lu, const int *pivots, complex<float> *b ) {

    int info;
    char trans;
    
    if( transpose ) 
        trans = 'T';
    else
        trans = 'N';

    cgetrs_(trans, ncol, nrhs, lu, ncol, pivots, b, ncol, info, 1  );

    return;
}
template <> void LapackInterface::getrs<complex<double> >
    ( const bool transpose, const int ncol, const int nrhs, const complex<double> *lu, const int *pivots, complex<double> *b ) {

    int info;
    char trans;
   
    if( transpose ) 
        trans = 'T';
    else
        trans = 'N';

    zgetrs_(trans, ncol, nrhs, lu, ncol, pivots, b, ncol, info, 1  );

    return;
}
template <> inline void LapackInterface::geev<double>
   (char jobvl, char jobvr,
    int n, double a[], int lda, double wr[], double wi[],  
    double vl[], int ldvl, double vr[], int ldvr, double work[],
    int lwork, int& info )
{

    dgeev_( jobvl, jobvr, 
             n, a, lda, wr, wi, vl, ldvl, vr, ldvr, 
             work, lwork, info, 
             1, 1);
}

template <> inline void LapackInterface::geev<float>
   (char jobvl, char jobvr,
    int n, float a[], int lda, float wr[], float wi[],  
    float vl[], int ldvl, float vr[], int ldvr, float work[],
    int lwork, int& info )
{

    sgeev_( jobvl, jobvr, 
             n, a, lda, wr, wi, vl, ldvl, vr, ldvr, 
             work, lwork, info, 
             1, 1);
}
template <> 
void  LapackInterface::getrf<double>( const int m, const int n, double *lu, const int lda,  int *pivots, int& info ) {
    dgetrf_(m, n, lu, lda, pivots, info   );
   return;
}
template <> 
void  LapackInterface::getrf<float>( const int m, const int n, float *lu, const int lda,  int *pivots, int& info ) {
    sgetrf_(m, n, lu, lda, pivots, info   );
   return;
}
template <> 
void  LapackInterface::getrf<std::complex<double> >( const int m, const int n, std::complex<double> *lu, const int lda,  int *pivots, int& info ) {
    zgetrf_(m, n, lu, lda, pivots, info   );
   return;
}
template <> 
void  LapackInterface::getrf<std::complex<float> >( const int m, const int n, std::complex<float> *lu, const int lda,  int *pivots, int& info ) {
    cgetrf_(m, n, lu, lda, pivots, info   );
   return;
}
template <>
void LapackInterface::potrf<double>( const int m, const int n, const int kl, const int ku, double* lu, const int lda, int *pivots, int& info ) { assert(false); }
template <>
void LapackInterface::potrf<float>( const int m, const int n, const int kl, const int ku, float* lu, const int lda, int *pivots, int& info ) { assert(false); }
template <>
void LapackInterface::potrf<std::complex<double> >( const int m, const int n, const int kl, const int ku, std::complex<double>* lu, const int lda, int *pivots, int& info ) { assert(false); }
template <>
void LapackInterface::potrf<std::complex<float> >( const int m, const int n, const int kl, const int ku, std::complex<float>* lu, const int lda, int *pivots, int& info ) { assert(false); }

template <> 
void LapackInterface::sytrf<float>( const char m, const int n, float* a,  const int lda, int *pivots, float* work, const int lwork, int& info ){ assert(false); }
template <> 
void LapackInterface::sytrf<double>( const char m, const int n, double* a,  const int lda, int *pivots, double* work, const int lwork, int& info ) { assert(false); }
template <> 
void LapackInterface::sytrf<std::complex<float> >( const char m, const int n, std::complex<float>* a,  const int lda, int *pivots, std::complex<float>* work, const int lwork, int& info ) { assert(false); }
template <> 
static void LapackInterface::sytrf<std::complex<double> >( const char m, const int n, std::complex<double>* a,  const int lda, int *pivots, std::complex<double>* work, const int lwork, int& info ) { assert(false); }

template <>
int LapackInterface::ilaenv<double>( int ispec,  const char* name,  const char *opts, int n1, int n2, int n3, int n4 ) { 
     char d[10];
     d[0] = 'd';
     d[1] = '\0';
     return (ilaenv_( ispec, strcat( d, name), opts, n1, n2, n3, n3, 6, strlen(opts)) ); 
}
template <>
int LapackInterface::ilaenv<float>( int ispec,  const char* name,  const char *opts, int n1, int n2, int n3, int n4 ) { 
     char s[10];
     s[0] = 's';
     s[1] = '\0';
     return (ilaenv_( ispec, strcat( s, name), opts, n1, n2, n3, n3, 6, strlen(opts)) ); 
}
template <>
int LapackInterface::ilaenv<std::complex<double> >( int ispec,  const char* name,  const char *opts, int n1, int n2, int n3, int n4 ) { 
     char z[10];
     z[0] = 'd';
     z[1] = '\0';
     return (ilaenv_( ispec, strcat( z, name), opts, n1, n2, n3, n3, 6, strlen(opts)) ); 
}
template <>
int LapackInterface::ilaenv<std::complex<float> >( int ispec,  const char* name,  const char *opts, int n1, int n2, int n3, int n4 ) { 
     char c[10];
     c[0] = 'd';
     c[1] = '\0';
     return (ilaenv_( ispec, strcat( c, name), opts, n1, n2, n3, n3, 6, strlen(opts)) ); 
}



}   // namespace SimTK

