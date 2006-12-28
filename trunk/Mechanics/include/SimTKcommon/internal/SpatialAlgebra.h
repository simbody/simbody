#ifndef SimTK_SIMMATRIX_SPATIAL_ALGEBRA_H_
#define SimTK_SIMMATRIX_SPATIAL_ALGEBRA_H_

/* Portions copyright (c) 2005-6 Stanford University and Michael Sherman.
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

/** @file
 *
 * These are declarations for special matrices and vectors of use in implementing
 * Rodriguez and Jain's Spatial Operator Algebra.
 */

#include "SimTKcommon/basics.h"
#include "SimTKcommon/internal/Scalar.h"
#include "SimTKcommon/internal/SmallMatrix.h"

#include <iostream>

namespace SimTK {

// Spatial vectors are used for (orientation,translation) quantities.
// These include
//      spatial velocity     = (angularVelocity,linearVelocity)
//      spatial acceleration = (angularAcceleration,linearAcceleration)
//      generalized forces   = (torque,force)
// Spatial configuration has to be handled differently though since
// orientation is not a vector quantity. (We use "Transform" for this concept
// which includes an orientation matrix and a translation vector.)
typedef Vec<2,   Vec3>  SpatialVec;
typedef Row<2,   Row3>  SpatialRow;
typedef Mat<2,2, Mat33> SpatialMat;

// support for efficient matrix multiplication involving the special phi
// matrix

using namespace SimTK;

class PhiMatrixTranspose;


inline static Vec3 
unitVec(const Vec3& v) { 
    Real tmpNorm = v.norm();
    if (tmpNorm < 1e-15) tmpNorm = 1e-15;

    Vec3 ret = v/tmpNorm;

    return ret;
}

class PhiMatrix {
public:
    typedef PhiMatrixTranspose TransposeType;

    PhiMatrix() { setToNaN(); }
    PhiMatrix(const Vec3& l) : l_(l) {}

    void setToZero() { l_ = 0.; }
    void setToNaN()  { l_.setToNaN(); }

    const Vec3& l() const { return l_; }
private:
    Vec3 l_;
};

class PhiMatrixTranspose {
public:
  PhiMatrixTranspose(const PhiMatrix& phi) : phi(phi) {}
  const Vec3& l() const {return phi.l();}
private:
  const PhiMatrix& phi;
};

inline PhiMatrixTranspose
transpose(const PhiMatrix& phi)
{
    PhiMatrixTranspose ret(phi);
    return ret;
}

inline PhiMatrixTranspose
operator~(const PhiMatrix& phi) {return transpose(phi);}

inline SpatialVec
operator*(const PhiMatrix&  phi,
          const SpatialVec& v)
{
    return SpatialVec(v[0] + phi.l() % v[1],
                      v[1]);
}

inline SpatialMat
operator*(const PhiMatrix&  phi,
          const SpatialMat& m)
{
    const Mat33 x = crossMat(phi.l());
    return SpatialMat( m(0,0) + x*m(1,0), m(0,1) + x*m(1,1),
                           m(1,0)       ,     m(1,1));
}

inline SpatialVec
operator*(const PhiMatrixTranspose& phiT,
          const SpatialVec&         v)
{
    return SpatialVec(v[0],
                      v[1] + v[0] % phiT.l());
}


inline SpatialMat
operator*(const SpatialMat::THerm&  m,
          const PhiMatrixTranspose& phiT)
{
    const Mat33 x = crossMat(phiT.l());
    return SpatialMat( m(0,0) - m(0,1) * x, m(0,1),
                       m(1,0) - m(1,1) * x, m(1,1) );
}

inline SpatialMat
operator*(const SpatialMat&         m,
          const PhiMatrixTranspose& phiT)
{
    const Mat33 x = crossMat(phiT.l());
    return SpatialMat( m(0,0) - m(0,1) * x, m(0,1),
                       m(1,0) - m(1,1) * x, m(1,1) );
}

inline bool
operator==(const PhiMatrix& p1, const PhiMatrix& p2)
{
    return p1.l() == p2.l();
}

inline bool
operator==(const PhiMatrixTranspose& p1, const PhiMatrixTranspose& p2)
{
    return p1.l() == p2.l();
}
} // namespace SimTK

#endif // SimTK_SIMMATRIX_SPATIAL_ALGEBRA_H_
