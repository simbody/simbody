#ifndef SimTK_SIMMATH_GCVSPL_UTIL_H_
#define SimTK_SIMMATH_GCVSPL_UTIL_H_

/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
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

#include "SimTKcommon.h"
#include "simmath/internal/common.h"

// These are global functions.
int SimTK_gcvspl_(const SimTK::Real *, const SimTK::Real *, int *, const SimTK::Real *, const SimTK::Real *, int *, int *,
                  int *, int *, SimTK::Real *, SimTK::Real *, int *, SimTK::Real *, int *);
SimTK::Real SimTK_splder_(int *, int *, int *, SimTK::Real *, const SimTK::Real *, const SimTK::Real *, int *, SimTK::Real *, int);

namespace SimTK {

/**
 * This class provides entry points for using the GCVSPL algorithm in terms of SimTK data types.
 * In most cases, you should use the SplineFitter class rather than invoking this class directly.
 * For details, see Woltring, H.J. (1986), A FORTRAN package for generalized, cross-validatory
 * spline smoothing and differentiation. Advances in Engineering Software 8(2):104-113.
 */

class SimTK_SIMMATH_EXPORT GCVSPLUtil {
public:
    static void gcvspl(const Vector& x, const Vector& y, const Vector& wx, Real wy, int m, int md, Real val, Vector& c, Vector& wk, int& ier);
    template <int K>
    static void gcvspl(const Vector& x, const Vector_<Vec<K> >& y, const Vector& wx, Vec<K> wy, int m, int md, Real val, Vector_<Vec<K> >& c, Vector& wk, int& ier);
    static Real splder(int derivOrder, int degree, Real t, const Vector& x, const Vector& coeff);
    template <int K>
    static Vec<K> splder(int derivOrder, int degree, Real t, const Vector& x, const Vector_<Vec<K> >& coeff);
};

template <int K>
void GCVSPLUtil::gcvspl(const Vector& x, const Vector_<Vec<K> >& y, const Vector& wx, Vec<K> wy, int degree, int md, Real val, Vector_<Vec<K> >& c, Vector& wk, int& ier) {
    SimTK_APIARGCHECK_ALWAYS(degree > 0 && degree%2==1, "GCVSPLUtil", "gcvspl", "degree must be positive and odd");
    SimTK_APIARGCHECK_ALWAYS(y.size() >= x.size(), "GCVSPLUtil", "gcvspl", "y is shorter than x");
    SimTK_APIARGCHECK_ALWAYS(wx.size() >= x.size(), "GCVSPLUtil", "gcvspl", "wx and x must be the same size");
    SimTK_APIARGCHECK_ALWAYS(x.hasContiguousData(), "GCVSPLUtil", "gcvspl", "x must have contiguous storage (i.e. not be a view)");
    SimTK_APIARGCHECK_ALWAYS(wk.hasContiguousData(), "GCVSPLUtil", "gcvspl", "wk must have contiguous storage (i.e. not be a view)");
    
    // Create various temporary variables.
    
    int m = (degree+1)/2;
    int n = x.size();
    int ny = y.size();
    Vector yvec(ny*K);
    int index = 0;
    for (int j = 0; j < K; ++j)
        for (int i = 0; i < ny; ++i)
            yvec[index++] = y[i][j];
    int nc = n*K;
    Vector cvec(nc);
    wk.resize(6*(m*n+1)+n);
    int k = K;
    
    // Invoke GCV.
    
    SimTK_gcvspl_(&x[0], &yvec[0], &ny, &wx[0], &wy[0], &m, &n, &k, &md, &val, &cvec[0], &n, &wk[0], &ier);
    if (ier != 0) {
        SimTK_APIARGCHECK_ALWAYS(n >= 2*m, "GCVSPLUtil", "gcvspl", "Too few data points");
        SimTK_APIARGCHECK_ALWAYS(ier != 2, "GCVSPLUtil", "gcvspl", "The values in x must be strictly increasing");
        SimTK_APIARGCHECK_ALWAYS(ier == 0, "GCVSPLUtil", "gcvspl", "GCVSPL returned an error code");
    }
    c.resize(n);
    index = 0;
    for (int j = 0; j < K; ++j)
        for (int i = 0; i < n; ++i)
            c[i][j] = cvec[index++];
}

template <int K>
Vec<K> GCVSPLUtil::splder(int derivOrder, int degree, Real t, const Vector& x, const Vector_<Vec<K> >& coeff) {
    assert(derivOrder >= 0);
    assert(t >= x[0] && t <= x[x.size()-1]);
    assert(x.size() == coeff.size());
    assert(degree > 0 && degree%2==1);
    assert(x.hasContiguousData());
    
    // Create various temporary variables.

    
    Vec<K> result;
    int m = (degree+1)/2;
    int n = x.size();
    int interval = (int) ceil(n*(t-x[0])/(x[n-1]-x[0]));

    const int MaxCheapM = 32;
    double qbuf[2*MaxCheapM];

    double *q = qbuf; // tentatively
    if (m > MaxCheapM)
        q = new double[2*m]; // use heap instead; don't forget to delete
    
    int offset = (int) (&coeff[1][0]-&coeff[0][0]);

    // Evaluate the spline one component at a time.
    
    for (int i = 0; i < K; ++i)
        result[i] = SimTK_splder_(&derivOrder, &m, &n, &t, &x[0], &coeff[0][i], &interval, q, offset);
    
    if (m > MaxCheapM)
        delete[] q;
    
    return result;
}

} // namespace SimTK

#endif // SimTK_SIMMATH_GCVSPL_UTIL_H_


