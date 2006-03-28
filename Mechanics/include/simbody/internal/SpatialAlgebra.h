#ifndef SIMTK_SIMBODY_SPATIAL_ALGEBRA_H_
#define SIMTK_SIMBODY_SPATIAL_ALGEBRA_H_

/** @file
 *
 * These are declarations for special matrices and vectors of use in implementing
 * Rodriguez and Jain's Spatial Operator Algebra.
 */

#include "simbody/internal/common.h"
#include "simbody/internal/Geometry.h"

#include <iostream>

namespace simtk {

// Spatial vectors are used for (orientation,translation) quantities.
// These include
//      spatial velocity     = (angularVelocity,linearVelocity)
//      spatial acceleration = (angularAcceleration,linearAcceleration)
//      generalized forces   = (torque,force)
// Spatial configuration has to be handled differently though since
// orientation is not a vector quantity. (We use "TransformMat" for this concept
// which includes an orientation matrix and a translation vector.)
typedef Vec<2,   Vec3>  SpatialVec;
typedef Row<2,   Row3>  SpatialRow;
typedef Mat<2,2, Mat33> SpatialMat;

// support for efficient matrix multiplication involving the special phi
// matrix

using namespace simtk;

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
} // namespace simtk

#endif // SIMTK_SIMBODY_SPATIAL_ALGEBRA_H_
