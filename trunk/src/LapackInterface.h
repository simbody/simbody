
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

template <class P> static void
    geev
   (char jobvl, char jobvr,
    int n, P a[], int lda, P wr[], P wi[],  
    P vl[], int ldvl, P vr[], int ldvr, P work[],
    int lwork, int& info ) {assert(false);}
};


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


}   // namespace SimTK

#endif // SimTK_LAPACK_INTERFACE_H_
