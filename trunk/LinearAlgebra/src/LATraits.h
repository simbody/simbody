#ifndef _LA_TRAITS_H
#define _LA_TRAITS_H
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



#include <complex>
#include "SimTKcommon.h"

namespace SimTK {

template <typename T> struct LATraits {
};
template <> struct LATraits<float>{
     static const bool isNegated    = false;
     static const bool isConjugated = false;
};
template <> struct LATraits<double>{
     static const bool isNegated     = false;
     static const bool isConjugated = false;
};
template <> struct LATraits<std::complex<float> >{
     static const bool isNegated    = false;
     static const bool isConjugated = false;
};
template <> struct LATraits<std::complex<double> >{
     static const bool isNegated     = false;
     static const bool isConjugated = false;
};
template <> struct LATraits<negator<float> >{
     static const bool isNegated    = true;
     static const bool isConjugated = false;
};
template <> struct LATraits<negator<double> >{
     static const bool isNegated     = true;
     static const bool isConjugated = false;
};
template <> struct LATraits<negator<std::complex<float> > >{
     static const bool isNegated    = true;
     static const bool isConjugated = false;
};
template <> struct LATraits<negator<std::complex<double> > >{
     static const bool isNegated     = true;
     static const bool isConjugated = false;
};
template <> struct LATraits<conjugate<float> >{
     static const bool isNegated    = false;
     static const bool isConjugated = true;
};
template <> struct LATraits<conjugate<double> >{
     static const bool isNegated     = false;
     static const bool isConjugated = true;
};
template <> struct LATraits<negator<conjugate<float> > >{
     static const bool isNegated    = true;
     static const bool isConjugated = true;
};
template <> struct LATraits<negator<conjugate<double> > >{
     static const bool isNegated     = true;
     static const bool isConjugated = true;
};

} // namespace SimTK
#endif   //  _LA_TRAITS_H
