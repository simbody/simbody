/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmatrix(tm)                       *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-7 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "SimTKcommon/Scalar.h"
#include "SimTKcommon/SmallMatrix.h"
#include "SimTKcommon/TemplatizedLapack.h"

#include "SimTKcommon/internal/BigMatrix.h"

#include "DataDescriptor.h"

#include <cstddef>
#include <cassert>
#include <complex>
#include <vector>

namespace SimTK {

template <class S> void
DataDescriptor<S>::fillWithScalar(const typename CNT<S>::StdNumber& val) {   
    for (int j=0; j<n; ++j) {
        S* const colBegin = &data[sindx(0,j)];
        for (int i=0; i<m; ++i) colBegin[i] = val;
    } 
}

template <class S> void
DataDescriptor<S>::fillWithElement(const int sz, const S* eltp) {
    if (sz==1) { fillWithScalar(typename CNT<S>::Number(*eltp)); return; }  
    for (int j=0; j<n; ++j) {
        S* const colBegin = &data[sindx(0,j)];
        for (int i=0; i<m; i+=sz)
            copyElement(sz,colBegin+i,eltp);
    } 
}

// Here we are given a pointer to a column's worth of elements, possibly
// with regular gaps between the elements. The idea is to copy them
// into a particular column of this data. Stride==1 means the source elements
// are packed together (that is the default for estride; see declaration). 
template <class S> void
DataDescriptor<S>::copyInOneEltColumn(int sz, int j, const S* eltp, int estride) {
    if (sz==1) { copyInOneScalarColumn(j,eltp,estride); return; }
    
    assert(j < n); assert(writable);   
    S* const colBegin = &data[sindx(0,j)];
    for (int i=0; i<m; i+=sz)
        copyElement(sz,colBegin+i,eltp+(i*estride));
}

template <class S> void
DataDescriptor<S>::copyInOneScalarColumn(int j, const S* sp, int stride) {
    assert(j<n); assert(writable);
    S* const colBegin = &data[sindx(0,j)];
    for (int i=0; i<m; ++i) colBegin[i] = sp[i*stride];
}

template <class S> void
DataDescriptor<S>::copyInDataByColumn(const S* sp, int ldim) {
    if (ldim==0) ldim=m;
    assert(ldim >= m); assert(writable);
    for (int j=0; j<n; ++j) {
        S* const colBegin = &data[sindx(0,j)];
        const S* srcColBeg = sp + j*ldim;
        for (int i=0; i<m; ++i) colBegin[i] = srcColBeg[i];
    }
}
template <class S> void
DataDescriptor<S>::invertInPlace() {
    assert(m==n && writable);
    // We don't care if this is negated, but conjugated won't work.
    assert(CNT<S>::IsStdNumber ||
           CNT<typename CNT<S>::TNeg>::IsStdNumber);

    typedef typename CNT<S>::StdNumber Raw;

    Raw* rawData = reinterpret_cast<Raw*>(data);
    std::vector<int> ipiv(m);
    int info;
    Lapack::getrf<Raw>(m,m,rawData,leadingDim,&ipiv[0],info);
    assert(info==0);

    // Calculate optimal size for work
    Raw workSz;
    Lapack::getri<Raw>(m,rawData,leadingDim,&ipiv[0],&workSz,-1,info);
    const int wsz = (int)CNT<Raw>::real(workSz);

    std::vector<Raw> work(wsz);
    Lapack::getri<Raw>(m,rawData,leadingDim,&ipiv[0],&work[0],wsz,info);
    assert(info==0);
}

// Instantiations for each of the 18 Scalar types.
template class DataDescriptor< float >;
template class DataDescriptor< double >;
template class DataDescriptor< long double >;
template class DataDescriptor< std::complex<float> >;
template class DataDescriptor< std::complex<double> >;
template class DataDescriptor< std::complex<long double> >;
template class DataDescriptor< SimTK::conjugate<float> >;
template class DataDescriptor< SimTK::conjugate<double> >;
template class DataDescriptor< SimTK::conjugate<long double> >;

template class DataDescriptor< SimTK::negator< float > >;
template class DataDescriptor< SimTK::negator< double > >;
template class DataDescriptor< SimTK::negator< long double > >;
template class DataDescriptor< SimTK::negator< std::complex<float> > >;
template class DataDescriptor< SimTK::negator< std::complex<double> > >;
template class DataDescriptor< SimTK::negator< std::complex<long double> > >;
template class DataDescriptor< SimTK::negator< SimTK::conjugate<float> > >;
template class DataDescriptor< SimTK::negator< SimTK::conjugate<double> > >;
template class DataDescriptor< SimTK::negator< SimTK::conjugate<long double> > >;

} // namespace SimTK

