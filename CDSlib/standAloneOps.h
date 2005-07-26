
#ifndef __standAloneOps_hh__
#define __standAloneOps_hh__

namespace StandAloneOps {

template<class T>
T
operator+(const T& v1,
	  const T& v2)
{
 return T(v1) += v2;
}

template<class T>
T
operator-(const T& v1,
	  const T& v2)
{
 return T(v1) -= v2;
}

//unary minus
template<class T>
T
operator-(const T& v1)
{
 return T(v1) *= -1;
}



};

#endif /* __standAloneOps_hh__ */
