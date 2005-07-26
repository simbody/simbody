
#ifndef __cdsmath_hh__
#define __cdsmath_hh__ 1

#include <sthead.h>
#include <cmath>

namespace CDSMath {

template<class T> inline
T max(T a, T b)
{
 return (a>b?a:b);
} /* max */

template<class T> inline
T min(T a, T b)
{
 return (a<b?a:b);
} /* min */

template<class T> T 
ipow(const T&  x,
	   int i)
/* return the ith integer power of x
   - works for negative integers too*/
{
 T temp = 1;
 for (int cnt=0 ; cnt< ::abs(i) ; cnt++) temp *= x;

 if (i<0 && temp != 0.0) temp = 1.0 / temp;

 return temp;
} /* pow */


template<class T> inline
T sq(const T t)
{
 return t*t;
}

template<class T> inline
int sign(T x)
{
 return ((x<0)?-1:1);
}

//
// conversion factor between radians and degrees
//
const double RAD2DEG = 180. / PI;

};

#endif /*__cdsmath_hh__*/
