#ifndef __phiMatrix_hh__
#define __phiMatrix_hh__

// support for efficient matrix multiplication involving the special phi
// matrix

#include "simbody/internal/SimbodyCommon.h"
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

    PhiMatrix() { l_(0) = 1.31e30; }
    PhiMatrix(const Vec3& l) : l_(l) {}

    const Vec3& l() const { assert( l_(0) != 1.31e30 ); return l_; }
private:
    Vec3 l_;
};

class PhiMatrixTranspose {
public:
  PhiMatrixTranspose(const PhiMatrix& phi) : phi(phi) {}
  const Vec3 l() const {return phi.l();}
private:
  const PhiMatrix& phi;
};

inline PhiMatrixTranspose
transpose(const PhiMatrix& phi)
{
    PhiMatrixTranspose ret(phi);
    return ret;
}
  
inline SpatialVec
operator*(const PhiMatrix&  phi,
          const SpatialVec& v)
{
    return SpatialVec(v[0] + phi.l() % v[1],
                      v[1]);
}

inline SpatialVec
operator*(const PhiMatrixTranspose& phiT,
          const SpatialVec&         v)
{
    return SpatialVec(v[0],
                      v[1] + v[0] % phiT.l());
}

inline SpatialMat
operator*(const SpatialMat&         m,
          const PhiMatrixTranspose& phiT)
{
    return SpatialMat( m(0,0) - m(0,1) * crossMat(phiT.l()), m(0,1),
                       m(1,0) - m(1,1) * crossMat(phiT.l()), m(1,1) );
}

#endif /*  __phiMatrix_hh__ */
