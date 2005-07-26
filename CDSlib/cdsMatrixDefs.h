
#ifndef __CDSMatrixDefs_h__
#define __CDSMatrixDefs_h__ 1

#include <iostream>
#include <iomanip>
using namespace std;

//
//constructor
//

template<class T> 
CDSMatrixBase<T>::CDSMatrixBase(const int size1,
				const int size2) 
{
 rep = new CDSMatrixRep<T>(size1,size2);
}

template<class T> 
CDSMatrixBase<T>::CDSMatrixBase(const int size1,
				const int size2,
				const T&  init)
{
 rep = new CDSMatrixRep<T>(size1,size2);
 for (int i=0 ; i<size1*size2 ; i++)
   rep->data[i] = init;
}

template<class T> 
CDSMatrixBase<T>::CDSMatrixBase(const int size1,
				const int size2,
				const T*  a)
  // constructor from array: not real safe
{
 rep = new CDSMatrixRep<T>(size1,size2);
 for (int i=0 ; i<size1*size2 ; i++)
   rep->data[i] = a[i];
}

//
//copy contructor
//

template<class T>
CDSMatrixBase<T>::CDSMatrixBase(const CDSMatrixBase<T> &m)
{
 rep = m.rep;
 rep->count++;
}

//template<class T>
//template<int OFFSET1, int OFFSET2>
//CDSMatrixBase<T>::CDSMatrixBase(const CDSMatrix<T,OFFSET1,OFFSET2> &m);   
//{
// rep = m.rep;
// rep->count++;
//}
//
//generic contructor
//

template<class T>
template<class MAT>
CDSMatrixBase<T>::CDSMatrixBase(const MAT &m)
{
 rep = new CDSMatrixRep<T>(m.rows(),m.cols());
 for (int i=0 ; i<rows() ; i++)
   for (int j=0 ; j<cols() ; j++)
     updData(i,j) = m(m.offset1()+i,m.offset2()+j);
}

//
//destructor
//

template<class T> 
CDSMatrixBase<T>::~CDSMatrixBase()
{
 if (--rep->count <= 0) delete rep;
}

template<class T>
CDSMatrixBase<T>& 
CDSMatrixBase<T>::resize(const int rows,
			 const int cols)
{
 if (rows != rep->size1 ||
     cols != rep->size2   ) {
   if (--rep->count <= 0) delete rep;
   rep = new CDSMatrixRep<T>(rows,cols);
 }
 
 //FIX: should there be a copy operation???
 return *this;
} /*resize*/

template<class T> 
CDSMatrixBase<T>&
CDSMatrixBase<T>::operator=(const CDSMatrixBase<T> &m)
{
 if (&m == this)
   return *this;

 m.rep->count++; 
 if (--rep->count <= 0) delete rep;

 rep = m.rep;

 return *this;
}

template<class T>
inline void
CDSMatrixBase<T>::set(const T &x)
{
 for (int i=0 ; i<rows()*cols() ; i++)
   rep->data[i] = x;
} /* set */


template<class T>
inline void
CDSMatrixBase<T>::setDiag(const T &x)
{
 for (int i=0 ; i<((rows()<cols())?rows():cols()) ; i++)
   updData(i,i) = x;
}


template<class T> 
template<class X>
CDSMatrixBase<T>&
CDSMatrixBase<T>::scale(const X& x)
{
 for (int i=0 ; i<rows()*cols() ; i++)
   rep->data[i] *= x;
 return *this;
} /* scale */


template<class T>
const CDSMatrixBase<T>&
CDSMatrixBase<T>::operator+=(const CDSMatrixBase<T>& m)
{
 assert( m.rep->size1 == rep->size1 );
 assert( m.rep->size2 == rep->size2 );
 splitRep();
  
 for (int i=0 ; i<rep->size1*rep->size2 ; i++)
   rep->data[i] += m.rep->data[i];
 return *this;
}

template<class T>
const CDSMatrixBase<T>&
CDSMatrixBase<T>::operator-=(const CDSMatrixBase<T>& m)
{
 assert( m.rep->size1 == rep->size1 );
 assert( m.rep->size2 == rep->size2 );
 splitRep();
  
 for (int i=0 ; i<rep->size1*rep->size2 ; i++)
   rep->data[i] -= m.rep->data[i];
 return *this;
}

template<class T> 
CDSMatrix<T>
operator+(const CDSMatrix<T>& m1,
	  const CDSMatrix<T>& m2)
{
 assert( m1.cols() == m2.cols() );
 assert( m1.rows() == m2.rows() );

 CDSMatrix<T> r = m1;
 r += m2;

 return r;
} /* operator+ (matrix,matrix) */

template<class T> 
CDSMatrix<T>
operator-(const CDSMatrix<T>& m1,
	  const CDSMatrix<T>& m2)
{
 assert( m1.cols() == m2.cols() );
 assert( m1.rows() == m2.rows() );

 CDSMatrix<T> r = m1;
 r -= m2;

 return r;
} /* operator- (matrix,matrix) */

template<class T> 
CDSMatrix<T>
operator*(const CDSMatrix<T>& m1,
	  const CDSMatrix<T>& m2)
{
 assert( m1.cols() == m2.rows() );

 CDSMatrix<T> r(m1.rows(),m2.cols(),(T)0);
 for (int i=0 ; i<m1.rows() ; i++)
   for (int k=0 ; k<m2.cols() ; k++)
     for (int j=0 ; j<m1.cols() ; j++) 
       r(i,k) += m1(i,j) * m2(j,k);

 return r;
} /* operator* (matrix,matrix) */

template<class T> 
CDSVector<T>
operator*(const CDSMatrixBase<T>& m,
	      const CDSVectorBase<T>& v)
{
 assert( m.cols() == v.size() );
 CDSVector<T> r(m.rows(),(T)0);
 
 for (int i=0 ; i<m.rows() ; i++) 
   for (int j=0 ; j<m.cols() ; j++) 
     r(i) += m.getData(i,j) * v.getData(j);
 
 return r;
} /* operator* (matrix, vector) */

//template<class T>
//CDSMatrix<T>
//transpose(const CDSMatrix<T> &m)
//{
// CDSMatrix<T> r(m.cols(),m.rows());
// for (int i=0 ; i<m.rows() ; i++)
//   for (int j=0 ; j<m.cols() ; j++) 
//     r(j,i) = m(i,j);
// return r;
//} /* transpose */


//generate a 2x2 block matrix
template<class T>
CDSMatrix<T>
blockMat22(const CDSMatrix<T>& m11,
	   const CDSMatrix<T>& m12,
	   const CDSMatrix<T>& m21,
	   const CDSMatrix<T>& m22)
{
 assert( m11.rows() == m12.rows() );
 assert( m11.cols() == m21.cols() );
 assert( m12.cols() == m22.cols() );
 assert( m21.rows() == m22.rows() );
 CDSMatrix<T> ret( m11.rows()+m21.rows() , m11.cols()+m12.cols() );
 for (int i=0 ; i<m11.rows() ; i++)
   for (int j=0 ; j<m11.cols() ; j++)
     ret(i,j) = m11(i,j);
 for (l_int i=0 ; i<m12.rows() ; i++)
   for (int j=0 ; j<m12.cols() ; j++)
     ret(i,j+m11.cols()) = m12(i,j);
 for (l_int i=0 ; i<m21.rows() ; i++)
   for (int j=0 ; j<m21.cols() ; j++)
     ret(i+m11.rows(),j) = m21(i,j);
 for (l_int i=0 ; i<m22.rows() ; i++)
   for (int j=0 ; j<m22.cols() ; j++)
     ret(i+m11.rows(),j+m11.cols()) = m22(i,j);
 return ret;
} /* blockMat22 */

template<class T>
CDSMatrix<T>
blockMat12(const CDSMatrix<T>& m1,
	   const CDSMatrix<T>& m2)
{
 assert( m1.rows() == m2.rows() );
 CDSMatrix<T> ret(m1.rows(),m1.cols()+m2.cols());
 for (int i=0 ; i<m1.rows() ; i++)
   for (int j=0 ; j<m1.cols() ; j++)
     ret(i,j) = m1(i,j);
 for (l_int i=0 ; i<m2.rows() ; i++)
   for (int j=0 ; j<m2.cols() ; j++)
     ret(i,j+m1.cols()) = m2(i,j);
 return ret;
} /* blockMat12 */

template<class T>
CDSMatrix<T>
blockMat21(const CDSMatrix<T>& m1,
	   const CDSMatrix<T>& m2)
{
 assert( m1.cols() == m2.cols() );
 CDSMatrix<T> ret(m1.rows()+m2.rows(),m1.cols());
 for (int i=0 ; i<m1.rows() ; i++)
   for (int j=0 ; j<m1.cols() ; j++)
     ret(i,j) = m1(i,j);
 for (l_int i=0 ; i<m2.rows() ; i++)
   for (int j=0 ; j<m2.cols() ; j++)
     ret(m1.rows()+i,j) = m2(i,j);
 return ret;
} /* blockMat21 */

//template<class T, int s1, int s2, int o1, int o2>
//CDSMatrix<T,s1,s2,o1,o2>
//CDSMatrix<T,s1,s2,o1,o2>::transpose()
//{
// CDSMatrix<T,s1,s2,o1,o2> r;
// for (int i=0 ; i<s1 ; i++)
//   for (int j=0 ; j<s2 ; j++) 
//     r(i,j) = (*this)(j,i);
// return r;
//} /* transpose */
//
//

//template<class T>
//static void
//ludcmp(CDSMatrix<T,1,1>& a,
//       CDSVector<int,1>& indx,
//       double&           d   )
//	    //from Numerical Recipes
//{
// assert( a.rows() == a.cols()    );
// assert( a.rows() == indx.size() );
// const int    size=a.cols();
// const double TINY=1e-20;//....
// CDSVector<double,1> vv(size);
//
// d = 1;
// for (int i=1 ; i<=size ; i++) {
//   double big=0;
//   for (int j=1 ; j<=size ; j++) {
//     double temp=fabs(a(i,j)); 
//     if ( temp > big ) big = temp;
//   }
//   if (big==0.0) throw CDS_NAMESPACE(SingularError)();
//   vv(i) = 1.0/big;
// }
// for (int j=1 ; j<=size ; j++) {
//   for (int i=1 ; i<j ; i++) {
//     T sum = a(i,j);
//     for (int k=1 ; k<i ; k++) sum -= a(i,k) * a(k,j);
//     a(i,j) = sum;
//   }
//   double big=0.0;
//   int imax=0;
//   for (l_int i=j ; i<=size ; i++) {
//     T sum = a(i,j);
//     for (int k=1 ; k<j ; k++) sum -= a(i,k) * a(k,j);
//     a(i,j) = sum;
//     double dum=vv(i)*fabs(sum);
//     if ( dum >= big) {
//       big = dum;
//       imax = i;
//     }
//   }
//   if (imax==0) throw CDS_NAMESPACE(SingularError)();
//   if (j != imax) {
//     for (int k=1 ; k<=size ; k++) {
//       T dum = a(imax,k);
//       a(imax,k) = a(j,k);
//       a(j,k) = dum;
//     }
//     d *= -1;
//     vv(imax) = vv(j);
//   }
//   indx(j) = imax;
//   if ( a(j,j) == 0.0) //singular matrix?
//     a(j,j) = TINY;
//   if (j != size) {
//     T dum = 1.0/a(j,j);
//     for (int i=j+1 ; i<=size ; i++) a(i,j) *= dum;
//   }
// }
//} /* ludcmp */

//template<class T>
//static void
//lubksb(      CDSMatrix<T,1,1>& a,
//       const CDSVector<int,1>& indx,
//	     CDSVector<T,1>&   b   )
//  //
//  // solve a x = b  (for x) -- place solution in b
//  //
//{
// assert( a.rows() == a.cols() );
// assert( a.rows() == indx.size() );
// assert( b.size() == indx.size() );
// const int size = a.rows();
// int ii=0;
// for (int i=1 ; i<=size ; i++) {
//   int ip = indx(i);
//   T sum=b(ip);
//   b(ip) = b(i);
//   if (ii)
//     for (int j=ii ; j<=i-1 ; j++)
//       sum -= a(i,j)*b(j);
//   else
//     if (sum) ii=i;
//   b(i) = sum;
// }
// for (l_int i=size ; i>=1 ; i--) {
//   T sum = b(i);
//   for (int j=i+1 ; j<=size ; j++)
//     sum -= a(i,j) * b(j);
//   b(i) = sum / a(i,i);
// }
//} /* lubksb */

//template<class T>
//CDSMatrix<T>
//inverse(const CDSMatrix<T>& m)
//{
// assert( m.rows()==m.cols() );
// const int size=m.rows();
// double d;
// CDSMatrix<T,1,1> a = m;
// //( (CDSMatrix<T,size,size,1,1>&)m );
// CDSVector<int,1> indx(size);
// ludcmp(a,indx,d);
// CDSMatrix<T> r(size,size);
// CDSVector<T,1> col(size);
// for (int j=1 ; j<=size ; j++) {
//   col.set( 0.0 );
//   col(j) = 1.0;
//   lubksb(a,indx,col);
//   for (int i=1 ; i<=size ; i++)
//     r(i-1,j-1) = col(i);
// }
// return r;
//} /* inverse */

     
//
// perform generalized orthogonal transform on m: S * m * transpose(S)
//
template<class T> 
CDSMatrix<T>
orthoTransform(const CDSMatrix<T>& m,
	       const CDSMatrix<T>& S)
  // return S * m * transpose(S);
{
 assert( m.rows() == m.cols() );
 assert( m.rows() == S.cols() );

 CDSMatrix<T> dum(S.rows(),S.cols(),(T)0);
 for (int i=0 ; i<dum.rows() ; i++) 
   for (int k=0 ; k<dum.cols() ; k++) 
     for (int j=0 ; j<m.rows() ; j++) 
       dum(i,k) += S(i,j) * m(j,k);

 CDSMatrix<T> ret(S.rows(),S.rows(),(T)0);
 for (l_int i=0 ; i<ret.rows() ; i++) 
   for (int k=0 ; k<ret.cols() ; k++) 
     for (int j=0 ; j<m.rows() ; j++) 
       ret(i,k) += dum(i,j) * S(k,j);
 
 return ret;
} /* orthoTransform */

template<class T, int moff1, int moff2>
CDSVector<T>
subCol(const CDSMatrix<T,moff1,moff2>& m,
       const int col,
       const int o1,
       const int o2)
  //
  // takes CDSMatrix as argument because the indices specifying range depend
  // on the offset.
  //
{ 
 CDSVector<T> ret(o2-o1+1);
 for (int i=o1 ; i<=o2 ; i++)
   ret(i-o1) = m(i,col);
 return ret;
} /* subVec(range in column of matrix) */

template<class T>
istream&   
operator>>(istream& s, CDSMatrixBase<T>& m)
{
 FMTFLAGS_TYPE flags = s.flags();
 s.setf(ios::skipws);
 char c; 
 s>>c; if ( c!='{' ) throw typename CDSMatrixBase<T>::IOError();
 for (int i=0 ; i<m.rows() ; i++) {
   s>>c; if ( c!='{' ) throw typename CDSMatrixBase<T>::IOError();
   s >> m.updData(i,0); // this will cause error if the matrix has a zero size!
   for (int j=1 ; j<m.cols() ; j++) {
     s >> c; if ( c!=',' ) throw typename CDSMatrixBase<T>::IOError();
     s >> m.updData(i,j);
     s>>c; if ( c!='}' ) throw typename CDSMatrixBase<T>::IOError();
   }
   if (i<m.rows()-1) {
     s >> c; if ( c!=',' ) throw typename CDSMatrixBase<T>::IOError();
   }
 }
 s>>c; if ( c!='}' ) throw typename CDSMatrixBase<T>::IOError();
 s.setf(flags);
 return s;
} /* operator>> */

template<class T>
ostream&   
operator<<(ostream& s, const CDSMatrixBase<T>& m)
{
 int width = s.width(); // apply width to numeric fields only
 s << "{ ";
 for (int i=0 ; i<m.rows() ; i++) {
   s << "{ ";
   s << setw(width) << m.getData(i,0); //will fail for vectors of size 0!
   for (int j=1 ; j<m.cols() ; j++) 
     s << ", " << setw(width) << m.getData(i,j) ;
   s << " }";
   if (i<m.rows()-1) s << ", ";
 }
 s << " }";
 return s;
}

#endif /* __cdsMatrixDefs_h__ */
