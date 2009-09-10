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

#include "SimTKcommon/Scalar.h"
#include "SimTKcommon/SmallMatrix.h"

using namespace SimTK;

template class Vec<3,Real>;
template class Vec<2,Real>;
template class Vec<3,Complex>;
template class Vec<2,Complex>;
template class Vec< 3,negator<Real> >;
template class Vec< 3,conjugate<float>,2 >;

template class Mat<3,3,Real>;
template class Mat<2,2,Real>;
template class Mat<3,3,Complex>;
template class Mat<2,2,Complex>;
template class Mat<5,5,negator< std::complex<double> > >;
template class Mat< 3,3,negator<Real> >;
template class Mat< 3,3,conjugate<float>,2 >;

template class SymMat< 4, std::complex<long double>, 7>;


template Real       SimTK::det(const Mat<1,1,Real>&);
template Real       SimTK::det(const SymMat<1,Real>&);
template Complex    SimTK::det(const Mat<2,2,Complex>&);
template Real       SimTK::det(const SymMat<2,Real>&);
template Real       SimTK::det(const Mat<3,3,Real>&);
template Real       SimTK::det(const SymMat<3,Real>&);
template Real       SimTK::det(const Mat<5,5,Real>&);
template Complex    SimTK::det(const SymMat<5,Complex>&);

template Mat<1,1,Real>::TInvert      SimTK::lapackInverse(const Mat<1,1,Real>&);
template Mat<2,2,Real>::TInvert      SimTK::lapackInverse(const Mat<2,2,Real>&);
template Mat<3,3,Conjugate>::TInvert SimTK::lapackInverse(const Mat<3,3,Conjugate>&);

template Mat<1,1,Real>::TInvert     SimTK::inverse(const Mat<1,1,Real>&);
template SymMat<1,Real>::TInvert    SimTK::inverse(const SymMat<1,Real>&);
template Mat<2,2,Complex>::TInvert  SimTK::inverse(const Mat<2,2,Complex>&);
template SymMat<2,Real>::TInvert    SimTK::inverse(const SymMat<2,Real>&);
template Mat<3,3,Real>::TInvert     SimTK::inverse(const Mat<3,3,Real>&);
template SymMat<3,Real>::TInvert    SimTK::inverse(const SymMat<3,Real>&);
template Mat<5,5,Real>::TInvert     SimTK::inverse(const Mat<5,5,Real>&);
//template SymMat<5,Complex>::TInvert SimTK::inverse(const SymMat<5,Complex>&);


