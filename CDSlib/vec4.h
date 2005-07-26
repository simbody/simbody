
#ifndef __vec4_hh__
#define __vec4_hh__

#include <fixedVector.h>
#include <math.h>

class Vec4 : public FixedVector<double,4> {
public: 
  Vec4() : FixedVector<double,4>() {}
  explicit Vec4(const double& x) : FixedVector<double,4>(x) {}
  explicit Vec4(const double* x) : FixedVector<double,4>(x) {}
  Vec4(const FixedVector<double,4> &v) : FixedVector<double,4>(v) {}
  explicit Vec4(const double& q1,
		const double& q2,
		const double& q3,
		const double& q4) {d_[0]=q1; d_[1]=q2; d_[2]=q3; d_[3]=q4;}
  template<class VEC>
  Vec4& operator=(const VEC& v) 
    { ((FixedVector<double,4>&) *this) = v; return *this; }

};

#endif /* __vec4_hh__ */
