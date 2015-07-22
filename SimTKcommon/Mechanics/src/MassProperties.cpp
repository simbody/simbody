/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
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

/**@file
 * Implementations of non-inline methods of MassProperties classes.
 */

#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/internal/MassProperties.h"

#include <iostream>

namespace SimTK {
    /////////////////////////
    //       INERTIA       //
    /////////////////////////

// Instantiate so we catch bugs now.
template class Inertia_<float>;
template class Inertia_<double>;

    /////////////////////////
    //     UNIT INERTIA    //
    /////////////////////////

// Instantiate so we catch bugs now.
template class UnitInertia_<float>;
template class UnitInertia_<double>;

    /////////////////////////
    //   MASS PROPERTIES   //
    /////////////////////////

// Instantiate so we catch bugs now.
template class MassProperties_<float>;
template class MassProperties_<double>;

    /////////////////////////
    //   SPATIAL INERTIA   //
    /////////////////////////

// Instantiate so we catch bugs now.
template class SpatialInertia_<float>;
template class SpatialInertia_<double>;

    /////////////////////////
    // ARTICULATED INERTIA //
    /////////////////////////

// Calculate the lower half of vx*F where vx is the cross product matrix
// of v and F is a full 3x3 matrix. This result would normally be a full
// 3x3 but for the uses below we know we're only going to need the diagonal
// and lower triangle so we can save some flops by working this out by hand.
// The method is templatized so that it will work on a transposed matrix
// as efficiently as an untransposed one. (18 flops)
template <class P, int CS, int RS>
static inline SymMat<3,P>
halfCross(const Vec<3,P>& v, const Mat<3,3,P,CS,RS>& F) {
    return SymMat<3,P>
      ( v[1]*F(2,0)-v[2]*F(1,0),
        v[2]*F(0,0)-v[0]*F(2,0), v[2]*F(0,1)-v[0]*F(2,1),
        v[0]*F(1,0)-v[1]*F(0,0), v[0]*F(1,1)-v[1]*F(0,1), v[0]*F(1,2)-v[1]*F(0,2) );
}

// Calculate the lower half of G*vx where G is a full 3x3 matrix and vx
// is the cross product matrix of v. See comment above for details.
// (18 flops)
template <class P, int CS, int RS>
static inline SymMat<3,P>
halfCross(const Mat<3,3,P,CS,RS>& G, const Vec<3,P>& v) {
    return SymMat<3,P>
      ( v[2]*G(0,1)-v[1]*G(0,2),
        v[2]*G(1,1)-v[1]*G(1,2), v[0]*G(1,2)-v[2]*G(1,0),
        v[2]*G(2,1)-v[1]*G(2,2), v[0]*G(2,2)-v[2]*G(2,0), v[1]*G(2,0)-v[0]*G(2,1) );
}

// This method computes the lower half of the difference vx*F-G*vx using
// the same methods as above, but done together in order to pull out the
// common v terms. This is 33 flops, down from 42 if you call the two
// methods above and add them.
template <class P, int CS1, int RS1, int CS2, int RS2>
static inline SymMat<3,P>
halfCrossDiff(const Vec<3,P>& v, const Mat<3,3,P,CS1,RS1>& F, const Mat<3,3,P,CS2,RS2>& G) {
    return SymMat<3,P>
      ( v[1]*(F(2,0)+G(0,2)) - v[2]*(F(1,0)+G(0,1)),
        v[2]*(F(0,0)-G(1,1)) - v[0]*F(2,0) + v[1]*G(1,2),
                v[2]*(F(0,1)+G(1,0)) - v[0]*(F(2,1)+G(1,2)),
        v[0]*F(1,0) - v[2]*G(2,1) - v[1]*(F(0,0)-G(2,2)),
                v[0]*(F(1,1)-G(2,2)) - v[1]*F(0,1) + v[2]*G(2,0),
                        v[0]*(F(1,2)+G(2,1)) - v[1]*(F(0,2)+G(2,0)) );
}

// We're computing
//      P' =  [ J'  F' ]  =  [ 1  sx ] [ J  F ] [ 1  0 ]
//            [~F'  M  ]     [ 0  1  ] [~F  M ] [-sx 1 ]
// like this:
//      F' = F + sx*M
//      J' = J + (sx*~F - F'*sx)
// where the parenthesized quantity is symmetric although its
// individual terms are not. Unfortunately this is a shift by -s rather than s.
// Cost is 72 flops.
template <class P> ArticulatedInertia_<P>
ArticulatedInertia_<P>::shift(const Vec3P& s) const {
    const Mat33P    Fp = F + s % M; // same meaning as sx*M but faster (33 flops)
    const SymMat33P Jp = J + halfCrossDiff(s, ~F, Fp); // sx*~F - F'*sx (39 flops)
    return ArticulatedInertia_(M, Fp, Jp);
}

// Same as above but perform the shift in place. Same flop count but less copying.
template <class P> ArticulatedInertia_<P>&
ArticulatedInertia_<P>::shiftInPlace(const Vec3P& s) {
    const Mat33P Fp = F + s % M;   // same meaning as sx*M but faster (33 flops)
    J += halfCrossDiff(s, ~F, Fp); // J + (sx*~F - F'*sx) (39 flops)
    F = Fp;
    // M doesn't change
    return *this;
}

// Instantiate so we catch bugs now.
template class ArticulatedInertia_<float>;
template class ArticulatedInertia_<double>;



} // namespace SimTK

