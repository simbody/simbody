#ifndef __fixedVectorDefs_h__
#define __fixedVectorDefs_h__ 1

#include <iostream>
#include <iomanip>
using namespace std;

template<class T, int SIZE>
template<class VEC>
FixedVectorBase<T,SIZE>&
FixedVectorBase<T,SIZE>::operator=(const CDS::GenericVector<VEC>& v)
{ 
 if ((FixedVectorBase<T,SIZE>*)&v == this) return *this; //self-assignment
 assert(SIZE==v.size());
 for (int i=0 ; i<SIZE ; i++)
   d_[i] = v(v.offset()+i);
 return *this;
} /* subVec(range in column of matrix) */



template<class T, int SIZE>
const FixedVectorBase<T,SIZE>&
FixedVectorBase<T,SIZE>::operator/=(const T &x)
  //
  // return v / x
  //
{
 for (int i=0 ; i<SIZE ; i++)
   d_[i] /= x;

  return *this;
} /* operator/= (T) */

template<class T, int SIZE>
const FixedVectorBase<T,SIZE>&
FixedVectorBase<T,SIZE>::operator-=(const FixedVectorBase<T,SIZE>& v)
{
 for (int i=0 ; i<SIZE ; i++)
   d_[i] -= v.d_[i];

 return *this;
} /* operator-= */

template<class T, int SIZE>
void
FixedVectorBase<T,SIZE>::copyArray(const T* a)
  //
  // copy contents of array into the vector. This assumes that the array is
  // is of the correct size.
{
 for (int i=0 ; i<SIZE ; i++)
   d_[i] = a[i];
} /* copyArray */

template<class T, int SIZE, int OFF>
template<int ms1, int ms2, int mOFF1, int mOFF2>
FixedVector<T,SIZE,OFF>
FixedVector<T,SIZE,OFF>::subCol(const FixedMatrix<T,ms1,ms2,mOFF1,mOFF2>& m,
				      int col,
				const int o1,
				const int o2)
  //
  // takes CDSAmtrix as argument because the indices specifying range depend
  // on the offset.
  //
{ 
 assert( SIZE==o2-o1+1 );
 FixedVector<T,SIZE,OFF> ret;
 for (int i=o1 ; i<=o2 ; i++)
   ret(i-o1+OFF) = m(i,col);
 return ret;
} /* subCol(range in column of matrix) */

template<class T,int SIZE1,int SIZE2>
FixedVector<T,SIZE1+SIZE2>
blockVec(const FixedVector<T,SIZE1>& v1,
	 const FixedVector<T,SIZE2>& v2)
{
 FixedVector<T,SIZE1+SIZE2> ret;
// const FixedVector<T>& v1r = v1;
// const FixedVector<T>& v2r = v2;
 for (int i=0 ; i<SIZE1 ; i++)
   ret(i) = v1(i);
 for (l_int i=0 ; i<SIZE2 ; i++)
   ret(SIZE1+i) = v2(i);
 return ret;
} /* blockVec */

template<class T,int SIZE>
T
norm(const FixedVector<T,SIZE>& v)
{
 using namespace CDS;
 return sqrt( dot(v,v) ); 
} /* norm */
template double norm(const FixedVector<double,3>&);

template<class T, int SIZE, int offset>
istream&   
operator>>(istream& s, FixedVector<T,SIZE,offset>& x)
{
 FMTFLAGS_TYPE flags = s.flags();
 s.setf(ios::skipws);
 char c; 
 s>>c; if ( c!='{' ) throw typename FixedVectorBase<T,SIZE>::IOError();
 s >> x(offset); // this will cause error if the vector is of zero size!
 for (int i=offset+1 ; i<SIZE+offset ; i++) {
   s >> c; if ( c!=',' ) throw typename FixedVectorBase<T,SIZE>::IOError();
   s >> x(i);
 }
 s>>c; if ( c!='}' ) throw typename FixedVectorBase<T,SIZE>::IOError();
 s.setf(flags);
 return s;
} /* operator>> */

template<class T, int SIZE>
ostream&   
operator<<(ostream& s, const FixedVectorBase<T,SIZE>& x)
{
 int width = s.width(); // apply width to numeric fields only
 s << setw(0) << "{ ";
 if ( SIZE ) {
   s << setw(width) << x(0); 
   for (int i=1 ; i<SIZE ; i++)
     s << ", " << setw(width) << x(i) ;
 }
 s << " }";
 return s;
} /* operator<< */

//template<class T, int SIZE>
//T
//distPow2(const FixedVectorBase<T,SIZE>& v1,
//	   const FixedVectorBase<T,SIZE>& v2,
//	   const int                  p2)
//  //
//  // perform |v1-v2|^(2*p2)
//  //
//{
// T ret = 0;
// for (int i=0 ; i <SIZE ; i++)
//   ret += sq( v1.d_[i] - v2.d_[i] );
// return pow(ret,p2);
//} /* distPow */

#endif /* __fixedVectorDefs_h__ */
