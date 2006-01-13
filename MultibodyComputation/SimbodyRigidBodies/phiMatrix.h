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
  
inline Vec6
operator*(const PhiMatrix& phi,
          const Vec6&   vec)
{
    const Vec3& v1 = vec.getSubVec<3>(0);
    const Vec3& v2 = vec.getSubVec<3>(3);

    Vec6 ret;

    ret.updSubVec<3>(0) = v1 + phi.l() % v2;
    ret.updSubVec<3>(3) = v2;

    return ret;
}

inline Vec6
operator*(const PhiMatrixTranspose& phiT,
          const Vec6&            vec)
{
    const Vec3& v1 = vec.getSubVec<3>(0);
    const Vec3& v2 = vec.getSubVec<3>(3);

    Vec6 ret;

    ret.updSubVec<3>(0) = v1;
    ret.updSubVec<3>(3) = v2 + v1 % phiT.l();

    return ret;
}

inline Mat66
operator*(const Mat66&              mat,
          const PhiMatrixTranspose& phiT)
{
    typedef Mat66::SubMat<3,3>::Type SubMat33;
    const SubMat33& m11 = mat.getSubMat<3,3>(0,0);
    const SubMat33& m12 = mat.getSubMat<3,3>(0,3);
    const SubMat33& m21 = mat.getSubMat<3,3>(3,0);
    const SubMat33& m22 = mat.getSubMat<3,3>(3,3);

    Mat66 ret;
    SubMat33& rm11 = ret.updSubMat<3,3>(0,0);
    SubMat33& rm12 = ret.updSubMat<3,3>(0,3);
    SubMat33& rm21 = ret.updSubMat<3,3>(3,0);
    SubMat33& rm22 = ret.updSubMat<3,3>(3,3);

    rm11 = m11 - m12 * crossMat(phiT.l());
    rm12 = m12;
    rm21 = m21 - m22 * crossMat(phiT.l());
    rm22 = m22;

    return ret;
}

#endif /*  __phiMatrix_hh__ */
