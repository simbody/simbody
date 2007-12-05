
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
    template <typename T> static void convertMatrixToLapack( typename CNT<T>::StdNumber* lapackArray, const Matrix_<T>& mat ); 
    template <typename T> static typename CNT<T>::StdNumber elementToLapack( const Matrix_<T>&mat, int i, int j);
    template <typename T> static typename CNT<T>::StdNumber elementToLapack( const Matrix_<negator<T> >&mat, int i, int j);
};

        

} // namespace SimTK
#endif   //  _LAPACK_CONVERT_H
