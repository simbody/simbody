
#ifndef __cdsComplex_hh__
#define __cdsComplex_hh__

#include <math.h>
#include <vector.h>
#include <cdsExcept.h>

#include <iostream>
using namespace std;

namespace CDS {

template <class T>
class Complex {
public:
  typedef T floatType;

  T re;
  T im;
  
  Complex(const T& re=0.,
	  const T& im=0.) : re(re), im(im) {}

  Complex<T>& operator+=(const Complex<T>& a);
  Complex<T>& operator*=(const Complex<T>& a);
  Complex<T>& operator/=(const Complex<T>& a);
  Complex<T>& operator/=(const T& a);
};

template<class T>
Complex<T>&
Complex<T>::operator+=(const Complex<T>& a)
{
 re += a.re; 
 im += a.im;
 return *this;
} /* operator += */

template<class T>
Complex<T>&
Complex<T>::operator*=(const Complex<T>& a)
{
 Complex<T> ret(re*a.re-im*a.im,
		re*a.im+im*a.re);
 (*this) = ret;
 return *this;
} /* operator *= */

template<class T>
Complex<T>&
Complex<T>::operator/=(const Complex<T>& z)
{
 T den = fabs(z.re) + fabs(z.im);
 if (den == 0.) 
   throw CDS::exception("Complex::/=: divide by zero.\n");

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
Complex<T>&
Complex<T>::operator/=(const T& a)
{
 re /= a;
 im /= a;
 return *this;
} /* operator /= */

template<class T>
inline Complex<T>
operator+(const Complex<T>& a,
	  const Complex<T>& b)
{
 return Complex<T>(a.re+b.re,a.im+b.im);
}

template<class T>
inline Complex<T>
operator-(const Complex<T>& a,
	  const Complex<T>& b)
{
 return Complex<T>(a.re-b.re,a.im-b.im);
}

template<class T>
inline Complex<T>
operator*(const Complex<T>& a,
	  const Complex<T>& b)
{
 return Complex<T>(a.re*b.re-a.im*b.im,
		   a.re*b.im+a.im*b.re);
}

template<class T>
inline Complex<T>
operator*(const T&          a,
	  const Complex<T>& b)
{
 return Complex<T>( a*b.re, a*b.im );
}

template<class T>
inline Complex<T>
operator/(const Complex<T>& a,
	  const Complex<T>& b)
{
 Complex<T> tmp = a;
 tmp /= b;
 return tmp;
}

template<class T>
inline T
norm(const Complex<T>& a)
{
 return ::sqrt(a.re*a.re + a.im*a.im);
}

template<class T>
inline T
abs(const Complex<T>& a)
{
 return norm(a);
}

template<class T>
inline ostream&
operator<<(ostream& s,
	   const Complex<T>& a)
{
 s << "{ " << a.re << ' ' << a.im << " }";
 return s;
}

template<class T>
inline Complex<T>
conj(const Complex<T>& z)
{
 return Complex<T>(z.re,-z.im);
} /* conj */

template<class T>
inline T
real(const Complex<T>& z)
{
 return z.re;
}

template<class T>
inline T
imag(const Complex<T>& z)
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
inline Complex<T>
exp(const Complex<T>& z)
{
 T r = ::exp(z.re);
 return Complex<T>(r * cos(z.im), 
		   r * sin(z.im));
} /* exp */

template<class T>
inline Complex<T>
sqrt(const Complex<T>& z)
{
 T r   = ::sqrt( abs(z) );
 T phi = 0.5 * atan2(z.im,z.re);
 return r*Complex<T>(cos(phi),sin(phi));
} /* sqrt */

} /* namespace CDS */


#endif /* __cdsComplex_hh__ */
