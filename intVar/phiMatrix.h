#ifndef __phiMatrix_hh__
#define __phiMatrix_hh__

// support for efficient matrix multiplication involving the special phi
// matrix

#include <vec3.h>
#include <fixedVector.h>
#include <fixedMatrix.h>
#include <subVector.h>
#include <subMatrix.h>

typedef FixedVector<double,6>   Vec6;
typedef FixedMatrix<double,6,6> Mat6;

class PhiMatrixTranspose;

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
          const Vec6&      vec)
{
    const SubVector<const Vec6> v1(vec,0,3);
    const SubVector<const Vec6> v2(vec,3,3);

    FixedVector<double,6> ret;
    SubVector<Vec6> rv1(ret,0,3);
    SubVector<Vec6> rv2(ret,3,3);

    rv1 = (v1 + cross(phi.l() , v2.vector())).vector();
    rv2 = v2.vector();

    return ret;
}

inline Vec6
operator*(const PhiMatrixTranspose&    phiT,
          const FixedVector<double,6>& vec)
{
    const SubVector<const Vec6> v1(vec,0,3);
    const SubVector<const Vec6> v2(vec,3,3);

    Vec6 ret;
    SubVector<Vec6> rv1(ret,0,3);
    SubVector<Vec6> rv2(ret,3,3);

    rv1 = v1.vector();
    rv2 = (v2 + cross(v1.vector() , phiT.l())).vector();

    return ret;
}

inline Mat6
operator*(const Mat6&               mat,
          const PhiMatrixTranspose& phiT)
{
    SubMatrix<const Mat6> m11(mat,0,0,3,3);
    SubMatrix<const Mat6> m12(mat,0,3,3,3);
    SubMatrix<const Mat6> m21(mat,3,0,3,3);
    SubMatrix<const Mat6> m22(mat,3,3,3,3);

    Mat6 ret;
    SubMatrix<Mat6> rm11(ret,0,0,3,3);
    SubMatrix<Mat6> rm12(ret,0,3,3,3);
    SubMatrix<Mat6> rm21(ret,3,0,3,3);
    SubMatrix<Mat6> rm22(ret,3,3,3,3);

    rm11 = m11 - m12 * crossMat(phiT.l());
    rm12 = m12;
    rm21 = m21 - m22 * crossMat(phiT.l());
    rm22 = m22;

    return ret;
}

#endif /*  __phiMatrix_hh__ */
