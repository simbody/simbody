#ifndef SimTK_SIMMATRIX_SCALAR_H_
#define SimTK_SIMMATRIX_SCALAR_H_

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

/** @file
 * This is a user-includable header which includes everything needed
 * to make use of SimMatrix Scalar code. More commonly, this will be
 * included from within Matrix code.
 */

#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/internal/conjugate.h"
#include "SimTKcommon/internal/NTraits.h"
#include "SimTKcommon/internal/CompositeNumericalTypes.h"
#include "SimTKcommon/internal/negator.h"

#include <complex>

namespace SimTK {

// Some scalar utilities
inline float       square(const float& x)       {return x*x;}
inline double      square(const double& x)      {return x*x;}
inline long double square(const long double& x) {return x*x;}

// Negation is free for negators, so we can square them and clean
// them up at the same time at no extra cost.
inline float       square(const negator<float>& x)       {return square(-x);}
inline double      square(const negator<double>& x)      {return square(-x);}
inline long double square(const negator<long double>& x) {return square(-x);}

// It is safer to templatize using complex classes, and doesn't make
// debugging any worse since complex is already templatized.
template <class P> inline 
std::complex<P> square(const std::complex<P>& x) {
    const P re=x.real(), im=x.imag();
    return std::complex<P>(re*re-im*im, 2*re*im);
}

// We can square a conjugate and clean it up back to complex at
// the same time at zero cost (or maybe 1 negation depending
// on how the compiler handles the "-2" below).
template <class P> inline 
std::complex<P> square(const conjugate<P>& x) {
    const P re=x.real(), nim=x.negImag();
    return std::complex<P>(re*re-nim*nim, -2*re*nim);
}

template <class P> inline
std::complex<P> square(const negator< std::complex<P> >& x) {
    return square(-x); // negation is free for negators
}

// Note return type here after squaring negator<conjugate>
// is complex, not conjugate.
template <class P> inline
std::complex<P> square(const negator< conjugate<P> >& x) {
    return square(-x); // negation is free for negators
}

inline float cube(const float& x) {return x*x*x;}
inline double cube(const double& x) {return x*x*x;}
inline long double cube(const long double& x) {return x*x*x;}

// To keep this cheap we'll defer getting rid of the negator<> until
// some other operation. We cube -x and then recast that to negator<>
// on return, for a total cost of 2 flops.
inline negator<float> cube(const negator<float>& x) {
    return negator<float>::recast(cube(-x));
}
inline negator<double> cube(const negator<double>& x) {
    return negator<double>::recast(cube(-x));
}
inline negator<long double> cube(const negator<long double>& x) {
    return negator<long double>::recast(cube(-x));
}

// Cubing a complex this way is a *lot* cheaper than doing it by
// multiplication. Cost here is 8 flops vs. 22 the other way.
template <class P> inline
std::complex<P> cube(const std::complex<P>& x) {
    const P re=x.real(), im=x.imag(), rr=re*re, ii=im*im;
    return std::complex<P>(re*(rr-3*ii), im*(3*rr-ii));
}

// Cubing a negated complex allows us to cube and eliminate the
// negator at the same time for zero extra cost. Compare the 
// expressions here to the normal cube above to see the free
// sign changes in both parts. 8 flops here.
template <class P> inline
std::complex<P> cube(const negator< std::complex<P> >& x) {
    const P nre=(-x).real(), nim=(-x).imag(), // -x is free for negator
            rr=nre*nre, ii=nim*nim;
    return std::complex<P>(nre*(3*ii-rr), nim*(ii-3*rr));
}

} // namespace SimTK

#endif //SimTK_SIMMATRIX_SCALAR_H_
