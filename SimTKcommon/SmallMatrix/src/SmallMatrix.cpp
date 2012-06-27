/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
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

template class SymMat< 4, std::complex<double>, 7>;


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


