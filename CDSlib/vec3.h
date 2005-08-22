#ifndef __vec3_hh__
#define __vec3_hh__

#include <fixedVector.h>
#include <vector.h>
#include <cdsMath.h>
#include "Mat33.h"

class Vec3 : public FixedVector<float_type,3> {
public: 
    Vec3() : FixedVector<float_type,3>() {}
    explicit Vec3(const float_type& x) : FixedVector<float_type,3>(x) {}
    explicit Vec3(const float_type* x) : FixedVector<float_type,3>(x) {}
    template <class VEC>
    Vec3(const CDS::Vector<VEC>& v) : FixedVector<float_type,3>(v) {}
    Vec3(const FixedVector<float_type,3> &v) : FixedVector<float_type,3>(v) {}

    explicit Vec3(const float_type& x,
	    const float_type& y,
	    const float_type& z) {d_[0]=x; d_[1]=y; d_[2]=z;}
    Vec3& operator=(const FixedVector<float_type,3>& v) 
    { ((FixedVector<float_type,3>&) *this) = v; return *this; }

    template<class VEC>
    Vec3& operator=(const CDS::Vector<VEC>& v)
    { ((FixedVector<float_type,3>&) *this) = v; return *this; }
    //template<class VEC>
    //Vec3& operator=(const CDSVector<VEC>& v)
    //{ FixedVectorBase<float_type,3>::operator=(v); return *this; }

    const float_type &x() const { return d_[0]; }
    const float_type &y() const { return d_[1]; }
    const float_type &z() const { return d_[2]; }
    float_type &x() { return d_[0]; }
    float_type &y() { return d_[1]; }
    float_type &z() { return d_[2]; }
};

// cross product
inline Vec3
cross(const Vec3& v1, const Vec3& v2) {
    return Vec3(v1.y()*v2.z() - v2.y()*v1.z(),
	            v1.z()*v2.x() - v2.z()*v1.x(),
	            v1.x()*v2.y() - v2.x()*v1.y());
}

inline Vec3 
unitVec(const Vec3& v) { 
    float_type tmpNorm = CDSMath::max( norm(v), 1e-15);

    if ( tmpNorm==0. )
        throw CDS::SingularError("unitVec(Vec3): singular vector");

    Vec3 ret = v;
    ret *= 1.0/tmpNorm;

    return ret;
}

//
// matrix associated with the cross product 
//  [ M s.t. M(v) z = v X z ]
//
inline static Mat33
crossMat(const Vec3& v) {
    Mat33 r(0.0);

    r(0,1) = -v(2); r(1,0) = -r(0,1);
    r(0,2) =  v(1); r(2,0) = -r(0,2);
    r(1,2) = -v(0); r(2,1) = -r(1,2);

    return r;
}

#endif /* __vec3_hh__ */
