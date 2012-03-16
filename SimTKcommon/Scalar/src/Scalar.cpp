/* -------------------------------------------------------------------------- *
 *                      SimTK Simbody: SimTKcommon                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-12 Stanford University and the Authors.        *
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

#include <limits>
#include <cmath>

#include <complex>
using std::complex;

namespace SimTK {


// These constants are global external symbols exported by the library. See
// the Scalar.h header file for information.

SimTK_SimTKCOMMON_EXPORT const Real NaN               = NTraits<Real>::getNaN(); 
SimTK_SimTKCOMMON_EXPORT const Real Infinity          = NTraits<Real>::getInfinity();
SimTK_SimTKCOMMON_EXPORT const Real Eps               = NTraits<Real>::getEps();
SimTK_SimTKCOMMON_EXPORT const Real SqrtEps           = NTraits<Real>::getSqrtEps();
SimTK_SimTKCOMMON_EXPORT const Real TinyReal          = NTraits<Real>::getTiny(); 
SimTK_SimTKCOMMON_EXPORT const Real SignificantReal   = NTraits<Real>::getSignificant(); 
SimTK_SimTKCOMMON_EXPORT const Real LeastPositiveReal = NTraits<Real>::getLeastPositive(); 
SimTK_SimTKCOMMON_EXPORT const Real MostPositiveReal  = NTraits<Real>::getMostPositive();  
SimTK_SimTKCOMMON_EXPORT const Real LeastNegativeReal = NTraits<Real>::getLeastNegative();
SimTK_SimTKCOMMON_EXPORT const Real MostNegativeReal  = NTraits<Real>::getMostNegative();

SimTK_SimTKCOMMON_EXPORT const int NumDigitsReal = NTraits<Real>::getNumDigits(); 
SimTK_SimTKCOMMON_EXPORT const int LosslessNumDigitsReal = NTraits<Real>::getLosslessNumDigits();

SimTK_SimTKCOMMON_EXPORT const Real Zero         = NTraits<Real>::getZero();
SimTK_SimTKCOMMON_EXPORT const Real One          = NTraits<Real>::getOne(); 
SimTK_SimTKCOMMON_EXPORT const Real MinusOne     = NTraits<Real>::getMinusOne();
SimTK_SimTKCOMMON_EXPORT const Real Two          = NTraits<Real>::getTwo(); 
SimTK_SimTKCOMMON_EXPORT const Real Three        = NTraits<Real>::getThree(); 

SimTK_SimTKCOMMON_EXPORT const Real OneHalf      = NTraits<Real>::getOneHalf();   
SimTK_SimTKCOMMON_EXPORT const Real OneThird     = NTraits<Real>::getOneThird();  
SimTK_SimTKCOMMON_EXPORT const Real OneFourth    = NTraits<Real>::getOneFourth(); 
SimTK_SimTKCOMMON_EXPORT const Real OneFifth     = NTraits<Real>::getOneFifth();  
SimTK_SimTKCOMMON_EXPORT const Real OneSixth     = NTraits<Real>::getOneSixth();  
SimTK_SimTKCOMMON_EXPORT const Real OneSeventh   = NTraits<Real>::getOneSeventh();
SimTK_SimTKCOMMON_EXPORT const Real OneEighth    = NTraits<Real>::getOneEighth(); 
SimTK_SimTKCOMMON_EXPORT const Real OneNinth     = NTraits<Real>::getOneNinth();  
SimTK_SimTKCOMMON_EXPORT const Real Pi           = NTraits<Real>::getPi();        
SimTK_SimTKCOMMON_EXPORT const Real OneOverPi    = NTraits<Real>::getOneOverPi(); 
SimTK_SimTKCOMMON_EXPORT const Real E            = NTraits<Real>::getE(); 
SimTK_SimTKCOMMON_EXPORT const Real Log2E        = NTraits<Real>::getLog2E(); 
SimTK_SimTKCOMMON_EXPORT const Real Log10E       = NTraits<Real>::getLog10E();
SimTK_SimTKCOMMON_EXPORT const Real Sqrt2        = NTraits<Real>::getSqrt2();
SimTK_SimTKCOMMON_EXPORT const Real OneOverSqrt2 = NTraits<Real>::getOneOverSqrt2();
SimTK_SimTKCOMMON_EXPORT const Real Sqrt3        = NTraits<Real>::getSqrt3();
SimTK_SimTKCOMMON_EXPORT const Real OneOverSqrt3 = NTraits<Real>::getOneOverSqrt3();
SimTK_SimTKCOMMON_EXPORT const Real CubeRoot2    = NTraits<Real>::getCubeRoot2();
SimTK_SimTKCOMMON_EXPORT const Real CubeRoot3    = NTraits<Real>::getCubeRoot3();
SimTK_SimTKCOMMON_EXPORT const Real Ln2          = NTraits<Real>::getLn2();
SimTK_SimTKCOMMON_EXPORT const Real Ln10         = NTraits<Real>::getLn10();

SimTK_SimTKCOMMON_EXPORT const Complex I = NTraits<Complex>::getI();

// These instantiations are just here to make sure everything is working. We would
// rather have these fail to compile here than in some poor user's program.
// (sherm 090827: also, the Intel compiler 11.1.038 seems to need some of these to be
// present in the library)

template class negator<float>;
template class negator<double>;
template class negator<long double>;

template class negator< complex<float> >;
template class negator< complex<double> >;
template class negator< complex<long double> >;

template class negator< conjugate<float> >;
template class negator< conjugate<double> >;
template class negator< conjugate<long double> >;

template class CNT< negator<float> >;
template class CNT< negator<double> >;
template class CNT< negator<long double> >;

template class CNT< complex<float> >;
template class CNT< complex<double> >;
template class CNT< complex<long double> >;

template class CNT< negator< complex<float> > >;
template class CNT< negator< complex<double> > >;
template class CNT< negator< complex<long double> > >;

template class CNT< conjugate<float> >;
template class CNT< conjugate<double> >;
template class CNT< conjugate<long double> >;

template class CNT< negator< conjugate<float> > >;
template class CNT< negator< conjugate<double> > >;
template class CNT< negator< conjugate<long double> > >;


#define INSTANTIATE_ALL_LEFT(T) \
template bool isNumericallyEqual(const T&, const complex<float>&,           double tol); \
template bool isNumericallyEqual(const T&, const complex<double>&,          double tol); \
template bool isNumericallyEqual(const T&, const complex<long double>&,     double tol); \
template bool isNumericallyEqual(const T&, const conjugate<float>&,         double tol); \
template bool isNumericallyEqual(const T&, const conjugate<double>&,        double tol); \
template bool isNumericallyEqual(const T&, const conjugate<long double>&,   double tol); \
template bool isNumericallyEqual(const T&, const float&,                    double tol); \
template bool isNumericallyEqual(const T&, const double&,                   double tol); \
template bool isNumericallyEqual(const T&, const long double&,              double tol); \
template bool isNumericallyEqual(const T&, int,                             double tol)

INSTANTIATE_ALL_LEFT(complex<float>);
INSTANTIATE_ALL_LEFT(complex<double>);
INSTANTIATE_ALL_LEFT(complex<long double>);
INSTANTIATE_ALL_LEFT(conjugate<float>);
INSTANTIATE_ALL_LEFT(conjugate<double>);
INSTANTIATE_ALL_LEFT(conjugate<long double>);

// Don't duplicate anything instantiated with the previous macro.
#define INSTANTIATE_ALL_RIGHT(T) \
template bool isNumericallyEqual(const float&,                  const T&, double tol); \
template bool isNumericallyEqual(const double&,                 const T&, double tol); \
template bool isNumericallyEqual(const long double&,            const T&, double tol); \
template bool isNumericallyEqual(int,                           const T&, double tol)

INSTANTIATE_ALL_RIGHT(complex<float>);
INSTANTIATE_ALL_RIGHT(complex<double>);
INSTANTIATE_ALL_RIGHT(complex<long double>);
INSTANTIATE_ALL_RIGHT(conjugate<float>);
INSTANTIATE_ALL_RIGHT(conjugate<double>);
INSTANTIATE_ALL_RIGHT(conjugate<long double>);

}
