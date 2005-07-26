#ifndef __fixedMatrixDefs_h__
#define __fixedMatrixDefs_h__ 1

template<class T, int s1, int s2>
void
FixedMatrixBase<T,s1,s2>::resize(int r,
				 int c) 
{ 
 if ( r!=rows() ||
      c!=cols()   ) 
   throw CDS::exception(String("FixedMatrix::resize: ") +
			"illegal resize operation attempted.");
} /* resize */
		
template<class T, int s1, int s2>
FixedMatrixBase<T,s1,s2>&
FixedMatrixBase<T,s1,s2>::operator+=(const FixedMatrixBase<T,s1,s2>& m)
{
  for (int i=0 ; i<s1*s2 ; i++) 
  d_[i] += m.d_[i];
  return *this;
} /* operator+= */

template<class T, int s1, int s2>
FixedMatrixBase<T,s1,s2>&
FixedMatrixBase<T,s1,s2>::operator-=(const FixedMatrixBase<T,s1,s2>& m)
{
  for (int i=0 ; i<s1*s2 ; i++) 
  d_[i] -= m.d_[i];
  return *this;
} /* operator-= */

template<class T, int S1, int S2>
FixedMatrixBase<T,S1,S2>&
FixedMatrixBase<T,S1,S2>::operator*=(const FixedMatrixBase<T,S1,S2>& m)
{
 FixedMatrixBase<T,S1,S2> tmp(*this);

 set(0.0);
 for (int i=0 ; i<S1 ; i++)
   for (int j=0 ; j<S2 ; j++)
     for (int k=0 ; k<S2 ; k++)
       data(i,j) += tmp.data(i,k) * m.data(k,j);

 return *this;
} /* operator*= */

template<class T, int m1, int n1, int n2>
FixedMatrix<T,m1,n1+n2>
blockMat12(const FixedMatrix<T,m1,n1>& M1,
	   const FixedMatrix<T,m1,n2>& M2)
{
 FixedMatrix<T,m1,n1+n2> ret;
 for (int i=0 ; i<m1 ; i++)
   for (int j=0 ; j<n1 ; j++)
     ret(i,j) = M1(i,j);
 for (l_int i=0 ; i<m1 ; i++)
   for (int j=0 ; j<n2 ; j++)
     ret(i,j+n1) = M2(i,j);
 return ret;
} /* blockMat12 */

template<class T, int m1, int n1, int m2>
FixedMatrix<T,m1+m2,n1>
blockMat21(const FixedMatrix<T,m1,n1>& M1,
	   const FixedMatrix<T,m2,n1>& M2)
{
 FixedMatrix<T,m1+m2,n1> ret;
 for (int i=0 ; i<m1 ; i++)
   for (int j=0 ; j<n1 ; j++)
     ret(i,j) = M1(i,j);
 for (l_int i=0 ; i<m2 ; i++)
   for (int j=0 ; j<n1 ; j++)
     ret(i+m1,j) = M2(i,j);
 return ret;
} /* blockMat21 */

//template<class T, int s1, int s2, int o1, int o2>
//FixedMatrix<T,s1,s2,o1,o2>
//FixedMatrix<T,s1,s2,o1,o2>::transpose()
//{
// FixedMatrix<T,s1,s2,o1,o2> r;
// for (int i=0 ; i<s1 ; i++)
//   for (int j=0 ; j<s2 ; j++) 
//     r(i,j) = (*this)(j,i);
// return r;
//} /* transpose */
//
//

template<class T, int size>
static void
ludcmp(FixedMatrix<T,size,size,1,1>& a,
       FixedVector<int,size,1>&      indx,
       double&                       d)
	    //from Numerical Recipes
{
 const double TINY=1e-20;//....
 FixedVector<double,size,1> vv;

 d = 1;
 for (int i=1 ; i<=size ; i++) {
   double big=0;
   for (int j=1 ; j<=size ; j++) {
     double temp=fabs(a(i,j)); 
     if ( temp > big ) big = temp;
   }
   if (big==0.0) throw CDS_NAMESPACE(SingularError)();
   vv(i) = 1.0/big;
 }
 for (int j=1 ; j<=size ; j++) {
   for (int i=1 ; i<j ; i++) {
     T sum = a(i,j);
     for (int k=1 ; k<i ; k++) sum -= a(i,k) * a(k,j);
     a(i,j) = sum;
   }
   double big=0.0;
   int imax=0;
   for (l_int i=j ; i<=size ; i++) {
     T sum = a(i,j);
     for (int k=1 ; k<j ; k++) sum -= a(i,k) * a(k,j);
     a(i,j) = sum;
     double dum=vv(i)*fabs(sum);
     if ( dum >= big) {
       big = dum;
       imax = i;
     }
   }
   if (j != imax) {
     for (int k=1 ; k<=size ; k++) {
       T dum = a(imax,k);
       a(imax,k) = a(j,k);
       a(j,k) = dum;
     }
     d *= -1;
     vv(imax) = vv(j);
   }
   indx(j) = imax;
   if ( a(j,j) == 0.0) //singular matrix?
     a(j,j) = TINY;
   if (j != size) {
     T dum = 1.0/a(j,j);
     for (int i=j+1 ; i<=size ; i++) a(i,j) *= dum;
   }
 }
} /* ludcmp */

template<class T, int size>
static void
lubksb(      FixedMatrix<T,size,size,1,1>& a,
       const FixedVector<int,size,1>&      indx,
	     FixedVector<T,size,1>&        b)
  //
  // solve a x = b  (for x) -- place solution in b
  //
{
 int ii=0;
 for (int i=1 ; i<=size ; i++) {
   int ip = indx(i);
   T sum=b(ip);
   b(ip) = b(i);
   if (ii)
     for (int j=ii ; j<=i-1 ; j++)
       sum -= a(i,j)*b(j);
   else
     if (sum) ii=i;
   b(i) = sum;
 }
 for (l_int i=size ; i>=1 ; i--) {
   T sum = b(i);
   for (int j=i+1 ; j<=size ; j++)
     sum -= a(i,j) * b(j);
   b(i) = sum / a(i,i);
 }
} /* lubksb */

//template<class T, int size>
//FixedMatrix<T,size>
//inverse(const FixedMatrix<T,size>& m)
//{
// double d;
// FixedMatrix<T,size,size,1,1> a = m;
// //( (FixedMatrix<T,size,size,1,1>&)m );
// FixedVector<int,size,1> indx;
// ludcmp(a,indx,d);
// FixedMatrix<T,size> r;
// FixedVector<T,size,1> col;
// for (int j=1 ; j<=size ; j++) {
//   col.set( 0.0 );
//   col(j) = 1.0;
//   lubksb(a,indx,col);
//   for (int i=1 ; i<=size ; i++)
//     r(i-1,j-1) = col(i);
// }
// return r;
//} /* inverse */

     
template<class T, int s1, int s2, int o1, int o2>
istream&   
operator>>(istream& s, FixedMatrix<T,s1,s2,o1,o2>& x)
{
 FMTFLAGS_TYPE flags = s.flags();
 s.setf(ios::skipws);
 char c; 
 s>>c; if ( c!='{' ) throw typename FixedMatrix<T,s1,s2,o1,o2>::IOError();
 for (int i=o1 ; i<s1+o1 ; i++) {
   s>>c; if ( c!='{' ) throw typename FixedMatrix<T,s1,s2,o1,o2>::IOError();
   s >> x(i,o2); // this will cause error if the matrix has a zero size!
   for (int j=o2+1 ; j<s2+o2 ; j++) {
     s >> c; if ( c!=',' ) throw typename FixedMatrix<T,s1,s2,o1,o2>::IOError();
     s >> x(i,j);
     s>>c; if ( c!='}' ) throw typename FixedMatrix<T,s1,s2,o1,o2>::IOError();
   }
   if (i<s1+o1-1) {
     s >> c; if ( c!=',' ) throw typename FixedMatrix<T,s1,s2,o1,o2>::IOError();
   }
 }
 s>>c; if ( c!='}' ) throw typename FixedMatrix<T,s1,s2,o1,o2>::IOError();
 s.setf(flags);
 return s;
} /* operator>> */

template<class T, int s1, int s2>
ostream&   
operator<<(ostream& s, const FixedMatrixBase<T,s1,s2>& x)
{
 s << "{ ";
 for (int i=0 ; i<s1 ; i++) {
   s << "{ ";
   s << x.data(i,0); //will fail for vectors of size 0!
   for (int j=1 ; j<s2 ; j++) 
     s << ", " << x.data(i,j) ;
   s << " }";
   if (i<s1-1) s << ", ";
 }
 s << " }";
 return s;
}


#endif /* __fixedMatrixDefs_h__ */
