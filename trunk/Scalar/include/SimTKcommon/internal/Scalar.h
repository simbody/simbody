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

// s=sign(n)
// Return int -1,0,1 according to n<0, n==0, n>0 for any integer
// or real numeric type. Sign is not defined for complex or conjugate.
// For unsigned argument, return unsigned 0 or 1.

inline unsigned int sign(unsigned char  u) {return u==0 ? 0 : 1;}
inline unsigned int sign(unsigned short u) {return u==0 ? 0 : 1;}
inline unsigned int sign(unsigned int   u) {return u==0 ? 0 : 1;}
inline unsigned int sign(unsigned long  u) {return u==0 ? 0 : 1;}

// Don't overload for plain "char" because it may be signed or unsigned
// depending on the compiler.

inline int sign(signed char i) {return i>0 ? 1 : (i<0 ? -1 : 0);}
inline int sign(short       i) {return i>0 ? 1 : (i<0 ? -1 : 0);}
inline int sign(int         i) {return i>0 ? 1 : (i<0 ? -1 : 0);}
inline int sign(long        i) {return i>0 ? 1 : (i<0 ? -1 : 0);}

inline int sign(const float&       x) {return x>0 ? 1 : (x<0 ? -1 : 0);}
inline int sign(const double&      x) {return x>0 ? 1 : (x<0 ? -1 : 0);}
inline int sign(const long double& x) {return x>0 ? 1 : (x<0 ? -1 : 0);}

inline int sign(const negator<float>&       x) {return -sign(-x);} // -x is free
inline int sign(const negator<double>&      x) {return -sign(-x);}
inline int sign(const negator<long double>& x) {return -sign(-x);}

// y=square(x)
// Return the square of the argument for any numeric type. We promise
// to evaluate x only once. We assume the result type is the same
// as the argument type; if it won't fit caller must cast argument
// to a wider type first.

inline unsigned char  square(unsigned char  u) {return u*u;}
inline unsigned short square(unsigned short u) {return u*u;}
inline unsigned int   square(unsigned int   u) {return u*u;}
inline unsigned long  square(unsigned long  u) {return u*u;}

inline char        square(char c) {return c*c;}

inline signed char square(signed char i) {return i*i;}
inline short       square(short       i) {return i*i;}
inline int         square(int         i) {return i*i;}
inline long        square(long        i) {return i*i;}

inline float       square(const float&       x) {return x*x;}
inline double      square(const double&      x) {return x*x;}
inline long double square(const long double& x) {return x*x;}

// Negation is free for negators, so we can square them and clean
// them up at the same time at no extra cost.
inline float       square(const negator<float>&       x) {return square(-x);}
inline double      square(const negator<double>&      x) {return square(-x);}
inline long double square(const negator<long double>& x) {return square(-x);}

// It is safer to templatize using complex classes, and doesn't make
// debugging any worse since complex is already templatized. 
// 5 flops vs. 6 for general complex multiply.
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

// y=cube(x)
// Return the cube of the argument for any numeric type. We promise
// to evaluate x only once.


inline unsigned char  cube(unsigned char  u) {return u*u*u;}
inline unsigned short cube(unsigned short u) {return u*u*u;}
inline unsigned int   cube(unsigned int   u) {return u*u*u;}
inline unsigned long  cube(unsigned long  u) {return u*u*u;}

inline char        cube(char c) {return c*c*c;}

inline signed char cube(signed char i) {return i*i*i;}
inline short       cube(short       i) {return i*i*i;}
inline int         cube(int         i) {return i*i*i;}
inline long        cube(long        i) {return i*i*i;}

inline float       cube(const float&       x) {return x*x*x;}
inline double      cube(const double&      x) {return x*x*x;}
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

// Cubing a complex this way is cheaper than doing it by
// multiplication. Cost here is 8 flops vs. 11 for a square
// followed by a multiply.
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
