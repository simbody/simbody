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
 * about BigMatrix. This occurs for slow or complicated operations where punting to the
 * BigMatrix classes is the easiest way to get the desired functionality.
 */

namespace SimTK {

// TODO: This is really bad! Should share space instead of recopying.
template <int M, int N, class ELT, int CS, int RS>
typename Mat<M,N,ELT,CS,RS>::TInvert 
Mat<M,N,ELT,CS,RS>::invert() const {
    Matrix_< EStandard > bigm(M,N);
    for (int j=0; j<N; ++j)
        for (int i=0; i<M; ++i)
            bigm(i,j) = (*this)(i,j);
    Matrix_< EInvert > result = bigm.invert();
    assert(result.nrow() == TInvert::NRows && result.ncol() == TInvert::NCols);
    TInvert out;
    for (int j=0; j<TInvert::NCols; ++j)
        for (int i=0; i<TInvert::NRows; ++i)
            out(i,j) = result(i,j);
    return out;
}

} //namespace SimTK


#endif // SimTK_SIMMATRIX_SMALL_DEFS_THAT_NEED_BIG_H_
