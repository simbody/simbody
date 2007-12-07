
#ifndef _LAPACK_CONVERT_H
#define _LAPACK_CONVERT_H
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


class LapackConvert {
    public:
    template <typename T, typename ELT> static void convertMatrixToLapack( T* lapackArray, const Matrix_<ELT>& mat ); 

	static void elementToLapack( float& e, const Matrix_<float>& m, int i, int j){ e =  m(i,j); }
	static void elementToLapack( double& e, const Matrix_<double>& m, int i, int j){ e = m(i,j); }
	static void elementToLapack( std::complex<float>& e,  const Matrix_<std::complex<float> >& m, int i, int j){ e = m(i,j); }
	static void elementToLapack( std::complex<double>& e, const Matrix_<std::complex<double> >& m, int i, int j){ e = m(i,j); }
	static void elementToLapack( float& e, const Matrix_<negator<float> >& m, int i, int j){ e = -m(i,j); }
    static void elementToLapack( double& e, const Matrix_<negator<double> >& m, int i, int j){ e = -m(i,j); }
    static void elementToLapack( std::complex<float>& e,  const Matrix_<negator<std::complex<float> > >& m, int i, int j){ 
        e = -m(i,j);
    }
    static void elementToLapack( std::complex<double>& e, const Matrix_<negator<std::complex<double> > >& m, int i, int j){ 
        e = -m(i,j);
    }
    static void elementToLapack( std::complex<float>& e, const Matrix_<conjugate<float> >&m, int i, int j){ e = conj(m(i,j)); }
 
	static void elementToLapack( std::complex<double>& e,  const Matrix_<conjugate<double> >&m, int i, int j){ e = conj(m(i,j)); }
   
    static void elementToLapack( std::complex<float>& e, const Matrix_<negator<conjugate<float> > >&m, int i, int j) {
         e = -conj(m(i,j));
    }
    static void elementToLapack( std::complex<double>& e, const Matrix_<negator<conjugate<double> > >&m, int i, int j) {
         e = -conj(m(i,j));
    }
};

        

} // namespace SimTK
#endif   //  _LAPACK_CONVERT_H
