
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

template <typename T> typename CNT<T>::StdNumber LapackConvert::elementToLapack( const Matrix_<T>&m, int i, int j) { 
    return( m(i,j) ); 
}
template float LapackConvert::elementToLapack( const Matrix_<float>&, int, int);
template double LapackConvert::elementToLapack( const Matrix_<double>&, int, int);
template std::complex<float> LapackConvert::elementToLapack( const Matrix_<std::complex<float> >&, int, int);
template std::complex<double> LapackConvert::elementToLapack( const Matrix_<std::complex<double> >&, int, int);


template <typename T> typename CNT<T>::StdNumber LapackConvert::elementToLapack( const Matrix_<negator<T> >&m, int i, int j) { 
    return(-m(i,j));
}
template float LapackConvert::elementToLapack( const Matrix_<negator<float> >&, int, int);
template double LapackConvert::elementToLapack( const Matrix_<negator<double> >&, int, int);
template std::complex<float> LapackConvert::elementToLapack( const Matrix_<negator<std::complex<float> > >&, int, int);
template std::complex<double> LapackConvert::elementToLapack( const Matrix_<negator<std::complex<double> > >&, int, int);


template <> std::complex<float> LapackConvert::elementToLapack( const Matrix_<conjugate<float> >&m, int i, int j) {
         return( conj(m(i,j)) );
     }
template <> std::complex<double> LapackConvert::elementToLapack( const Matrix_<conjugate<double> >&m, int i, int j) {
         return( conj(m(i,j)) );
     }
template <> std::complex<float> LapackConvert::elementToLapack( const Matrix_<negator<conjugate<float> > >&m, int i, int j) {
         return( -conj(m(i,j)) );
     }
template <> std::complex<double> LapackConvert::elementToLapack( const Matrix_<negator<conjugate<double> > >&m, int i, int j) {
         return( -conj(m(i,j)) );
     }

template < typename T>  
void LapackConvert::convertMatrixToLapack ( typename CNT<T>::StdNumber* lapackArray, const Matrix_<T>& mat ) {
    int m = mat.nrow();
    int n = mat.ncol();
    for(int i=0;i<n;i++) {
        for(int j=0;j<m;j++)  {
            lapackArray[i*m+j] = elementToLapack( mat, j, i );
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
template void LapackConvert::convertMatrixToLapack( std::complex<float>*  lapackArray, const Matrix_<negator<conjugate<float> > >& mat );
template void LapackConvert::convertMatrixToLapack( std::complex<double>* lapackArray, const Matrix_<negator<conjugate<double> > >& mat );

} // namespace SimTK
