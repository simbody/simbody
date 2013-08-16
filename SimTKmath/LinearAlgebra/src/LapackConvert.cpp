/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-12 Stanford University and the Authors.        *
 * Authors: Jack Middleton                                                    *
 * Contributors:                                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

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
    for(int c=0;c<n;c++) {
        for(int r=0;r<m;r++)  {
             lapackArray[c*m+r] = mat(r,c);
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
