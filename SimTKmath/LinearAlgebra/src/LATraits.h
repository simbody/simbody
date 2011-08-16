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

template <typename T> struct LANT;
        
template <> struct LANT<float>{
     static const int sign =  1;
     static const bool conjugate = false;
};
template <> struct LANT<double>{
     static const int sign     = 1;
     static const bool conjugate = false;
};
template <> struct LANT<std::complex<float> >{
     static const int sign    = 1;
     static const bool conjugate = false;
};
template <> struct LANT<std::complex<double> >{
     static const int sign     = 1;
     static const bool conjgate = false;
};
template <> struct LANT<negator<float> >{
     static const int sign    = -1;
     static const bool conjugate = false;
};
template <> struct LANT<negator<double> >{
     static const int sign     = -1;
     static const bool conjugate = false;
};
template <> struct LANT<negator<std::complex<float> > >{
     static const int sign    = -1;
     static const bool conjugate = false;
};
template <> struct LANT<negator<std::complex<double> > >{
     static const int sign     = -1;
     static const bool conjugate = false;
};
template <> struct LANT<conjugate<float> >{
     static const int sign    = 1;
     static const bool conjugate = true;
};
template <> struct LANT<conjugate<double> >{
     static const int sign     = 1;
     static const bool conjugate = true;
};
template <> struct LANT<negator<conjugate<float> > >{
     static const int sign    = -1;
     static const bool conjugate = true;
};
template <> struct LANT<negator<conjugate<double> > >{
     static const int sign     = -1;
     static const bool conjugate = true;
};

} // namespace SimTK
#endif   //  _LA_TRAITS_H
