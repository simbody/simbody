
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
#include "LapackConvert.h"

namespace SimTK {


template <>
void LapackConvert::convertMatrixToLapack( std::complex<double>* lapackArray,  const Matrix_<negator<conjugate<double> > >& mat ) {
    int m = mat.nrow();
    int n = mat.ncol();
    for(int i=0;i<n;i++) {
        for(int j=0;j<m;j++)  {
             lapackArray[i*m+j] = std::complex<double>(mat(j,i).real(),mat(j,i).imag());
        }
    }
    return;
}

template <>
void LapackConvert::convertMatrixToLapack( std::complex<float>* lapackArray,  const Matrix_<negator<conjugate<float> > >& mat ) {
    int m = mat.nrow();
    int n = mat.ncol();
    for(int i=0;i<n;i++) {
        for(int j=0;j<m;j++)  {
             lapackArray[i*m+j] = std::complex<float>(mat(j,i).real(),mat(j,i).imag());
        }
    }
    return;
}
template < typename T, typename ELT>
void LapackConvert::convertMatrixToLapack ( T* lapackArray,  const Matrix_<ELT>& mat ) {
    int m = mat.nrow();
    int n = mat.ncol();
    for(int i=0;i<n;i++) {
        for(int j=0;j<m;j++)  {
             lapackArray[i*m+j] = mat(j,i);
        }
    }
    return;
}


template void LapackConvert::convertMatrixToLapack( float*  lapackArray, const Matrix_<float>& mat );
template void LapackConvert::convertMatrixToLapack( double* lapackArray, const Matrix_<double>& mat );
template void LapackConvert::convertMatrixToLapack( float*  lapackArray, const Matrix_<negator<float> >& mat );
template void LapackConvert::convertMatrixToLapack( double* lapackArray, const Matrix_<negator<double> >& mat );
template void LapackConvert::convertMatrixToLapack( std::complex<float>*  lapackArray, const Matrix_<std::complex<float> >& mat );
template void LapackConvert::convertMatrixToLapack( std::complex<double>* lapackArray, const Matrix_<std::complex<double> >& mat );
template void LapackConvert::convertMatrixToLapack( std::complex<float>*  lapackArray, const Matrix_<negator<std::complex<float> > >& mat );
template void LapackConvert::convertMatrixToLapack( std::complex<double>* lapackArray, const Matrix_<negator<std::complex<double> > >& mat );
template void LapackConvert::convertMatrixToLapack( std::complex<float>*  lapackArray, const Matrix_<conjugate<float> >& mat );
template void LapackConvert::convertMatrixToLapack( std::complex<double>* lapackArray, const Matrix_<conjugate<double> >& mat );

} // namespace SimTK
