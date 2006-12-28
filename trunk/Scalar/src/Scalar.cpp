/* Portions copyright (c) 2005-6 Stanford University and Michael Sherman.
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

#include "SimTKcommon/basics.h"
#include "SimTKcommon/internal/Scalar.h"

#include <limits>
#include <cmath>

#include <complex>
using std::complex;

namespace SimTK {

// Calculate the basic NTraits constants at each precision. These have to be
// done in the right order!

#define CALC_CONSTANTS(T) \
const T NTraits<T>::Pi            = (T)(SimTK_PI);       \
const T NTraits<T>::OneOverPi     = (T)(1/SimTK_PI);     \
const T NTraits<T>::E             = (T)(SimTK_E);        \
const T NTraits<T>::Log2E         = (T)(SimTK_LOG2E);    \
const T NTraits<T>::Log10E        = (T)(SimTK_LOG10E);   \
const T NTraits<T>::Sqrt2         = (T)(SimTK_SQRT2);    \
const T NTraits<T>::OneOverSqrt2  = (T)(1/SimTK_SQRT2);  \
const T NTraits<T>::Sqrt3         = (T)(SimTK_SQRT3);    \
const T NTraits<T>::OneOverSqrt3  = (T)(1/SimTK_SQRT3);  \
const T NTraits<T>::CubeRoot2     = (T)(SimTK_CBRT2);    \
const T NTraits<T>::CubeRoot3     = (T)(SimTK_CBRT3);    \
const T NTraits<T>::Ln2           = (T)(SimTK_LN2);      \
const T NTraits<T>::Ln10          = (T)(SimTK_LN10);     \
const T NTraits<T>::OneThird      = (T)(1.L/3.L);        \
const T NTraits<T>::OneSixth      = (T)(1.L/6.L);        \
const T NTraits<T>::OneSeventh    = (T)(1.L/7.L);        \
const T NTraits<T>::OneNinth      = (T)(1.L/9.L);        \
const T NTraits<T>::NaN           = std::numeric_limits<T>::quiet_NaN();    \
const T NTraits<T>::Infinity      = std::numeric_limits<T>::infinity();     \
const T NTraits<T>::Eps           = std::numeric_limits<T>::epsilon();      \
const T NTraits<T>::Eps_12        = std::sqrt(Eps);      \
const T NTraits<T>::Eps_13        = std::pow (Eps, OneThird);   \
const T NTraits<T>::Eps_14        = std::sqrt(Eps_12);   \
const T NTraits<T>::Eps_15        = std::pow (Eps, (T)0.2L);    \
const T NTraits<T>::Eps_18        = std::sqrt(Eps_14);   \
const T NTraits<T>::Eps_23        = Eps_13*Eps_13;       \
const T NTraits<T>::Eps_25        = Eps_15*Eps_15;       \
const T NTraits<T>::Eps_34        = Eps_12*Eps_14;       \
const T NTraits<T>::Eps_35        = Eps_15*Eps_25;       \
const T NTraits<T>::Eps_38        = Eps_18*Eps_14;       \
const T NTraits<T>::Eps_45        = Eps_25*Eps_25;       \
const T NTraits<T>::Eps_58        = Eps_18*Eps_12;       \
const T NTraits<T>::Eps_78        = Eps_18*Eps_34;       \
const T NTraits<T>::Eps_54        = Eps*Eps_14;          \
const T NTraits<T>::Tiny          = Eps_54

CALC_CONSTANTS(float);
CALC_CONSTANTS(double);
CALC_CONSTANTS(long double);

// These instantiations are just here to make sure everything is working. We would
// rather have these fail to compile here than in some poor user's program.
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

template class CNT< negator<conjugate<float> > >;
template class CNT< negator<conjugate<double> > >;
template class CNT< negator<conjugate<long double> > >;

}
