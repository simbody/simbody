
#ifndef __cdsComplex_hh__
#define __cdsComplex_hh__

#include <math.h>
#include <vector.h>
#include <cdsExcept.h>

#include <iostream>
using namespace std;

namespace CDS {

template <class T>
class CDSComplex {
public:
  typedef T floatType;

  T re;
  T im;
  
  CDSComplex(const T& re=0.,
	  const T& im=0.) : re(re), im(im) {}

  CDSComplex<T>& operator+=(const CDSComplex<T>& a);
  CDSComplex<T>& operator*=(const CDSComplex<T>& a);
  CDSComplex<T>& operator/=(const CDSComplex<T>& a);
  CDSComplex<T>& operator/=(const T& a);
};

template<class T>
CDSComplex<T>&
CDSComplex<T>::operator+=(const CDSComplex<T>& a)
{
 re += a.re; 
 im += a.im;
 return *this;
} /* operator += */

template<class T>
CDSComplex<T>&
CDSComplex<T>::operator*=(const CDSComplex<T>& a)
{
 CDSComplex<T> ret(re*a.re-im*a.im,
		re*a.im+im*a.re);
 (*this) = ret;
 return *this;
} /* operator *= */

template<class T>
CDSComplex<T>&
CDSComplex<T>::operator/=(const CDSComplex<T>& z)
{
 T den = fabs(z.re) + fabs(z.im);
 if (den == 0.) 
   throw CDS::exception("CDSComplex::/=: divide by zero.\n");

 T xrden = re / den;
 T xiden = im / den;
 T yrden = z.re / den;
 T yiden = z.im / den;
 T nrm   = yrden * yrden + yiden * yiden;
 re = (xrden * yrden + xiden * yiden) / nrm;
 im = (xiden * yrden - xrden * yiden) / nrm;
 return *this;
} /* operator /= */

template<class T>
CDSComplex<T>&
CDSComplex<T>::operator/=(const T& a)
{
 re /= a;
 im /= a;
 return *this;
} /* operator /= */

template<class T>
inline CDSComplex<T>
operator+(const CDSComplex<T>& a,
	  const CDSComplex<T>& b)
{
 return CDSComplex<T>(a.re+b.re,a.im+b.im);
}

template<class T>
inline CDSComplex<T>
operator-(const CDSComplex<T>& a,
	  const CDSComplex<T>& b)
{
 return CDSComplex<T>(a.re-b.re,a.im-b.im);
}

template<class T>
inline CDSComplex<T>
operator*(const CDSComplex<T>& a,
	  const CDSComplex<T>& b)
{
 return CDSComplex<T>(a.re*b.re-a.im*b.im,
		   a.re*b.im+a.im*b.re);
}

template<class T>
inline CDSComplex<T>
operator*(const T&          a,
	  const CDSComplex<T>& b)
{
 return CDSComplex<T>( a*b.re, a*b.im );
}

template<class T>
inline CDSComplex<T>
operator/(const CDSComplex<T>& a,
	  const CDSComplex<T>& b)
{
 CDSComplex<T> tmp = a;
 tmp /= b;
 return tmp;
}

template<class T>
inline T
norm(const CDSComplex<T>& a)
{
 return ::sqrt(a.re*a.re + a.im*a.im);
}

template<class T>
inline T
abs(const CDSComplex<T>& a)
{
 return norm(a);
}

template<class T>
inline ostream&
operator<<(ostream& s,
	   const CDSComplex<T>& a)
{
 s << "{ " << a.re << ' ' << a.im << " }";
 return s;
}

template<class T>
inline CDSComplex<T>
conj(const CDSComplex<T>& z)
{
 return CDSComplex<T>(z.re,-z.im);
} /* conj */

template<class T>
inline T
real(const CDSComplex<T>& z)
{
 return z.re;
}

template<class T>
inline T
imag(const CDSComplex<T>& z)
{
 return z.im;
}

template<class VECTOR>
VECTOR
conj(const Vector<VECTOR>& v)
{
 VECTOR ret=v.v;
 for (int i=0 ; i<ret.size() ; i++)
   ret(i) = conj(ret(i));

 return ret;
} /* conj(vector) */

template<class T>
inline CDSComplex<T>
exp(const CDSComplex<T>& z)
{
 T r = ::exp(z.re);
 return CDSComplex<T>(r * cos(z.im), 
		   r * sin(z.im));
} /* exp */

template<class T>
inline CDSComplex<T>
sqrt(const CDSComplex<T>& z)
{
 T r   = ::sqrt( abs(z) );
 T phi = 0.5 * atan2(z.im,z.re);
 return r*CDSComplex<T>(cos(phi),sin(phi));
} /* sqrt */

} /* namespace CDS */


#endif /* __cdsComplex_hh__ */
