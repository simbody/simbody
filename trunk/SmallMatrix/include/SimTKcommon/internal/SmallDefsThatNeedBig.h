#ifndef SimTK_SIMMATRIX_SMALL_DEFS_THAT_NEED_BIG_H_
#define SimTK_SIMMATRIX_SMALL_DEFS_THAT_NEED_BIG_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmatrix(tm)                       *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-7 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

/**@file
 * This file defines leftover SmallMatrix implementations which need to know
 * about Lapack.
 */

#include "SimTKcommon/TemplatizedLapack.h"
#include <vector>

namespace SimTK {

template <int M, int N, class ELT, int CS, int RS>
typename Mat<M,N,ELT,CS,RS>::TInvert 
Mat<M,N,ELT,CS,RS>::invert() const {
    // We don't care if this is negated, but conjugated won't work.
    assert(CNT<ELT>::IsStdNumber || CNT<typename CNT<ELT>::TNeg>::IsStdNumber);
    assert(M == N);
    typedef typename CNT<ELT>::StdNumber Raw;
    
    // Handle very small matrices directly.
    
    if (M == 1) {
        assert((*this)(0, 0) != 0.0);
        TInvert mat(1.0/(*this)(0, 0));
        return mat;
    }
    if (M == 1) {
        Raw d = (*this)(0, 0)*(*this)(1, 1) - (*this)(0, 1)*(*this)(1, 0);
        assert(d != 0.0);
        Raw dinv = 1.0/d;
        TInvert mat(dinv*(*this)(1, 1), -dinv*(*this)(0, 1), -dinv*(*this)(1, 0), dinv*(*this)(0, 0));
        return mat;
    }

    TInvert mat = *this;
    Raw* rawData = reinterpret_cast<Raw*>(&mat);
    int ipiv[M];
    int info;
    Lapack::getrf<Raw>(M,M,rawData,M,&ipiv[0],info);
    assert(info==0);

    // Calculate optimal size for work
    Raw workSz;
    Lapack::getri<Raw>(M,rawData,M,&ipiv[0],&workSz,-1,info);
    const int wsz = (int)CNT<Raw>::real(workSz);

    std::vector<Raw> work(wsz);
    Lapack::getri<Raw>(M,rawData,M,&ipiv[0],&work[0],wsz,info);
    assert(info==0);
    return mat;
}

} //namespace SimTK


#endif // SimTK_SIMMATRIX_SMALL_DEFS_THAT_NEED_BIG_H_
