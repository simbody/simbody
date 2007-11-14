
#ifndef SimTK_LAPACK_INTERFACE_H_
#define SimTK_LAPACK_INTERFACE_H_

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

namespace SimTK {

class LapackInterface { public:
    // MEANINGLESS IF NOT SPECIALIZED

template <class P> 
    static void geev (char jobvl, char jobvr,
    int n, P a[], int lda, P wr[], P wi[],  
    P vl[], int ldvl, P vr[], int ldvr, P work[],
    int lwork, int& info ) {assert(false);}


/* solve system of linear equations using the LU factorization  computed by getrf */
template <class T>
    static void getrs( const bool transpose, const int ncol, const int nrhs, const T *lu,
            const int *pivots, T *b )   { assert(false); }

template <class T>
    static void getrf( const int m, const int n, T *a, const int lda, int *pivots, int& info ) { assert(false); }
template <class T>
    static void gttrf( const int m, const int n, T* dl, T* d, T* du, T* du2, int *pivots, int& info ) { assert(false); }
template <class T>
    static void gbtrf( const int m, const int n, const int kl, const int ku, T* lu, const int lda, int *pivots, int& info ) { assert(false); }
template <class T>
    static void potrf( const int m, const int n, const int kl, const int ku, T* lu, const int lda, int *pivots, int& info ) { assert(false); }
template <class T>
    static void sytrf( const char m, const int n, T* a,  const int lda, int *pivots, T* work, const int lwork, int& info ) { assert(false); }

template <class T>
    static int ilaenv( int ispec,  const char* name,  const char *opts, int n1, int n2, int n3, int n4  ) { assert(false); }

}; // class LapackInterface

}   // namespace SimTK

#endif // SimTK_LAPACK_INTERFACE_H_
