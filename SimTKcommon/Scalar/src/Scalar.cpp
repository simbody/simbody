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

#include <limits>
#include <cmath>

#include <complex>
using std::complex;

namespace SimTK {


// These constants are global external symbols exported by the library. See
// the Scalar.h header file for information.

constexpr Real NaN      = NTraits<Real>::getNaN();
constexpr Real Infinity = NTraits<Real>::getInfinity();
constexpr Real Eps      = NTraits<Real>::getEps();
constexpr Real SqrtEps        = NTraits<Real>::getSqrtEps();
constexpr Real TinyReal       = NTraits<Real>::getTiny();
constexpr Real SignificantReal = NTraits<Real>::getSignificant();
constexpr Real LeastPositiveReal = NTraits<Real>::getLeastPositive();
constexpr Real MostPositiveReal  = NTraits<Real>::getMostPositive();
constexpr Real LeastNegativeReal = NTraits<Real>::getLeastNegative();
constexpr Real MostNegativeReal  = NTraits<Real>::getMostNegative();

constexpr int NumDigitsReal = NTraits<Real>::getNumDigits();
constexpr int LosslessNumDigitsReal = NTraits<Real>::getLosslessNumDigits();

constexpr Real Zero         = NTraits<Real>::getZero();
constexpr Real One          = NTraits<Real>::getOne();
constexpr Real MinusOne     = NTraits<Real>::getMinusOne();
constexpr Real Two          = NTraits<Real>::getTwo();
constexpr Real Three        = NTraits<Real>::getThree();

constexpr Real OneHalf      = NTraits<Real>::getOneHalf();
constexpr Real OneThird     = NTraits<Real>::getOneThird();
constexpr Real OneFourth    = NTraits<Real>::getOneFourth();
constexpr Real OneFifth     = NTraits<Real>::getOneFifth();
constexpr Real OneSixth     = NTraits<Real>::getOneSixth();
constexpr Real OneSeventh   = NTraits<Real>::getOneSeventh();
constexpr Real OneEighth    = NTraits<Real>::getOneEighth();
constexpr Real OneNinth     = NTraits<Real>::getOneNinth();
constexpr Real Pi           = NTraits<Real>::getPi();
constexpr Real OneOverPi    = NTraits<Real>::getOneOverPi();
constexpr Real E            = NTraits<Real>::getE();
constexpr Real Log2E        = NTraits<Real>::getLog2E();
constexpr Real Log10E       = NTraits<Real>::getLog10E();
constexpr Real Sqrt2        = NTraits<Real>::getSqrt2();
constexpr Real OneOverSqrt2 = NTraits<Real>::getOneOverSqrt2();
constexpr Real Sqrt3        = NTraits<Real>::getSqrt3();
constexpr Real OneOverSqrt3 = NTraits<Real>::getOneOverSqrt3();
constexpr Real CubeRoot2    = NTraits<Real>::getCubeRoot2();
constexpr Real CubeRoot3    = NTraits<Real>::getCubeRoot3();
constexpr Real Ln2          = NTraits<Real>::getLn2();
constexpr Real Ln10         = NTraits<Real>::getLn10();

const Complex I = NTraits<Complex>::getI();

// These instantiations are just here to make sure everything is working. We would
// rather have these fail to compile here than in some poor user's program.
// (sherm 090827: also, the Intel compiler 11.1.038 seems to need some of these to be
// present in the library)

template class negator<float>;
template class negator<double>;

template class negator< complex<float> >;
template class negator< complex<double> >;

template class negator< conjugate<float> >;
template class negator< conjugate<double> >;

template class CNT< negator<float> >;
template class CNT< negator<double> >;

template class CNT< complex<float> >;
template class CNT< complex<double> >;

template class CNT< negator< complex<float> > >;
template class CNT< negator< complex<double> > >;

template class CNT< conjugate<float> >;
template class CNT< conjugate<double> >;

template class CNT< negator< conjugate<float> > >;
template class CNT< negator< conjugate<double> > >;


#define INSTANTIATE_ALL_LEFT(T) \
template bool isNumericallyEqual(const T&, const complex<float>&,           double tol); \
template bool isNumericallyEqual(const T&, const complex<double>&,          double tol); \
template bool isNumericallyEqual(const T&, const conjugate<float>&,         double tol); \
template bool isNumericallyEqual(const T&, const conjugate<double>&,        double tol); \
template bool isNumericallyEqual(const T&, const float&,                    double tol); \
template bool isNumericallyEqual(const T&, const double&,                   double tol); \
template bool isNumericallyEqual(const T&, int,                             double tol)

INSTANTIATE_ALL_LEFT(complex<float>);
INSTANTIATE_ALL_LEFT(complex<double>);
INSTANTIATE_ALL_LEFT(conjugate<float>);
INSTANTIATE_ALL_LEFT(conjugate<double>);

// Don't duplicate anything instantiated with the previous macro.
#define INSTANTIATE_ALL_RIGHT(T) \
template bool isNumericallyEqual(const float&,                  const T&, double tol); \
template bool isNumericallyEqual(const double&,                 const T&, double tol); \
template bool isNumericallyEqual(int,                           const T&, double tol)

INSTANTIATE_ALL_RIGHT(complex<float>);
INSTANTIATE_ALL_RIGHT(complex<double>);
INSTANTIATE_ALL_RIGHT(conjugate<float>);
INSTANTIATE_ALL_RIGHT(conjugate<double>);

}
