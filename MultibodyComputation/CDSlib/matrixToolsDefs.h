#ifndef __matrixToolsDefs_h__
#define __matrixToolsDefs_h__

#include "matrixTools.h"
#include "lapack.h"
#include "cdsAuto_arr.h"
#include "cdsExcept.h"
#include "cdsMath.h"
#include "cdsSStream.h"

using CDS::CDSComplex;


namespace MatrixTools {

template<class MATRIX>
CDSVector<typename MATRIX::ElementType>
getColumn(const MATRIX& m,
	      const int     col)
{
    CDSVector<typename MATRIX::ElementType> ret(m.nrow());
    for (int i=0 ; i<m.nrow() ; i++)
        ret(i) = m(i+m.offset1(),col);
    return ret;
}


template<class MATRIX>
CDSVector<typename MATRIX::ElementType>
getRow(const MATRIX& m,
       const int     row)
{
    CDSVector<typename MATRIX::ElementType> ret(m.ncol());
    for (int i=0 ; i<m.ncol() ; i++)
        ret(i) = m(row,i+m.offset2());
    return ret;
}


template<class MATRIX, class VECTOR> void
setColumn(      MATRIX& m,
	      const int     col,
	      const VECTOR& v)
{
    for (int i=0 ; i<m.nrow() ; i++)
        m(i+m.offset1(),col) = v(i+v.offset());
}


template<class MATRIX, class VECTOR> void
setRow(      MATRIX& m,
       const int     row,
       const VECTOR& v)
{
    for (int i=0 ; i<m.ncol() ; i++)
        m(row,i+m.offset2()) = v(i+v.offset());
}


template<class MATRIX>
typename MATRIX::TransposeType
transpose(const MATRIX& m)
{
    typename MATRIX::TransposeType r;
    r.resize(m.ncol(),m.nrow());
    for (int i=0 ; i<m.nrow() ; i++)
        for (int j=0 ; j<m.ncol() ; j++) 
            r(j,i) = m(i,j);
    return r;
}

} // namespace MatrixTools

namespace MatrixTools {
template<class MATRIX> MATRIX
callInverse(const MATRIX& matrix,
		    InverseResults<FullMatrix<typename MATRIX::ElementType> >)
{
    assert( matrix.nrow()==matrix.ncol() );

    const int size=matrix.nrow();

    CDSVector<int,1> ipiv(size);

    MATRIX ret;
    ret.resize(size,size);

    for (int i=0 ; i<size ; i++)
        for (int j=0 ; j<size ; j++)
            ret(i,j) = matrix(i,j); //row major

    int info=-100;
    callTRF(size,size,ret.pointer(),size,ipiv.pointer(),info);

    if ( info<0 )
        throw CDS::exception("inverse: bad argument to callTrf");
    if ( info>0 )
        throw CDS::SingularError("inverse");

    info = -100;

    CDSVector<typename MATRIX::ElementType> work(size*3);
    callTRI(size,ret.pointer(),size,ipiv.pointer(),
            work.pointer(),work.size(),info);
    if ( info<0 )
        throw CDS::exception("inverse: bad argument to callTri");
    if ( info>0 )
        throw CDS::SingularError();

    return ret;
}
}

namespace MatrixTools {
template<class MATRIX>
MATRIX
callInverse(const MATRIX& matrix,
		    InverseResults<SymmetricMatrix<typename MATRIX::ElementType> >)
{
    typedef typename MATRIX::ElementType T;

    const int size=matrix.nrow();

    CDSVector<int,1> ipiv(size);

    MATRIX ret = matrix;

    int info=0;
    callSPTRF( 'L', size, ret.pointer(), ipiv.pointer(),info);

    if ( info<0 )
        throw CDS::exception("inverse: bad argument to callSPTRF");
    if ( info>0 )
        throw CDS::SingularError("inverse");

    info = -100;

    CDSVector<T> work(size);
    callSPTRI( 'L', size, ret.pointer(), ipiv.pointer(), work.pointer(), info );

    if ( info<0 )
        throw CDS::exception("inverse: bad argument to callSPTRI");
    if ( info>0 )
        throw CDS::SingularError();

    return ret;
}
}

namespace MatrixTools {
template<class MATRIX> MATRIX
inverse(const MATRIX& matrix,
	    InverseResults<typename MATRIX::MatrixType> results)
{
    typedef typename MATRIX::ElementType T;

    if ( matrix.nrow()==0 ) return matrix;

    return callInverse(matrix,results);
}
}

template<class Number>           // a placeholder template function:
static void                      // must be specialized.
callSVD(const char   &JOBU, 
	    const char   &JOBVT,
	    const int    &M, 
	    const int    &N, 
	    Number        A[],
	    const int    &LDA,
	    Number        S[],
	    Number        U[],
	    const int    &LDU, 
	    Number        VT[], 
	    const int    &LDVT, 
	    Number        WORK[],
	    const int    &LWORK, 
	    int          &INFO )
{ throw CDS::exception("callSVD not defined"); }

template<class Number>           // a placeholder template function:
static void                      // must be specialized.
callTRI(const int &N,
	    Number     A[],
	    const int &LDA,
	    const int  IPIV[], 
	    Number     WORK[], 
	    const int &LWORK, 
	    int       &INFO )
{ throw CDS::exception("callTri not defined"); }

template<class Number>           // a placeholder template function:
static void                      // must be specialized.
callSPTRI(const char& UPLO,
	      const int&  size,
		  Number      A[],
		  int         IPIV[],
		  Number      WORK[], 
		  int        &INFO )
{ throw CDS::exception("callSPTRI not defined"); }

template<class Number>           // a placeholder template function:
static void                      // must be specialized.
callTRF(const int &M,
	    const int &N, 
	    Number A[],
	    const int &LDA, 
	    int IPIV[], 
	    int &INFO )
{ throw CDS::exception("callTrf not defined"); }

template<class Number>           // a placeholder template function:
static void                      // must be specialized.
callSPTRF(const char& UPLO,
	      const int& size,
		  Number A[],
		  int IPIV[], 
		  int &INFO )
{ throw CDS::exception("callSPTRF not defined"); }


using namespace MatrixTools;

// for full matrix
template<class MATRIX>
static EigenResults<typename MATRIX::MatrixType>
callEigen(const MATRIX&              matrix,
	      const char                 jobz,
		  EigenResults<FullMatrix<typename MATRIX::ElementType> > ret)
{
    typedef typename MATRIX::ElementType T;
    MATRIX trashedByLapack = matrix;
    ret.work.resize(3*matrix.ncol());

    CDSVector<T > values_r(matrix.nrow());
    CDSVector<T > values_i(matrix.nrow());
    CDSMatrix<T > vectors_r(1,matrix.nrow());
    CDSMatrix<T > vectors_l(1,matrix.nrow());

    char jobvl = 'N';
    char jobvr = 'N';
    if ( jobz == 'V' || jobz == 'R' ) {
        jobvr = 'V';
        ret.work.resize(4*matrix.nrow());
        vectors_r.resize(matrix.nrow(),matrix.ncol());
    }
    if ( jobz == 'L' ) {
        jobvl = 'V';
        ret.work.resize(4*matrix.nrow());
        vectors_l.resize(matrix.nrow(),matrix.ncol());
    }

    callEigenFull(jobvl,jobvr,
                  matrix.ncol(),trashedByLapack.pointer(),matrix.ncol(),
                  values_r.pointer(),values_i.pointer(),
                  vectors_l.pointer(),vectors_l.nrow(),
                  vectors_r.pointer(),vectors_r.nrow(),
                  ret.work.pointer(),ret.work.size(),
                  ret.info);

    if (ret.info>0) {
        OStringStream s;
        s << "eigen: " << ret.info 
          << " QR algorithm failed to converge." << ends;
        cerr << s.str() << endl;
        throw CDS::exception(s.str());
    } else if (ret.info<0) {
        OStringStream s;
        s << "eigen: argument " << -ret.info 
          << " has an illegal value." << ends;
        cerr << s.str() << endl;
        throw CDS::exception(s.str());
    }

    ret.eigenPairs.resize(matrix.ncol());
    for (int i=0 ; i<matrix.ncol() ; i++) {
        ret.eigenPairs[i].value = CDSComplex<T>( values_r(i), values_i(i) );
        ret.eigenPairs[i].vector.resize( matrix.ncol() );
    }

    if ( jobvl == 'V' )
        for (int i=0 ; i<matrix.ncol() ; i++)
            for (int j=0 ; j<matrix.ncol() ; j++)
                if ( values_i(i) != 0. ) {
                    ret.eigenPairs[i].vector(j) 
                        = CDSComplex<T>(vectors_l(j,i), vectors_l(j,i+1));
                    ret.eigenPairs[i+1].vector(j) 
                        = CDSComplex<T>(vectors_l(j,i), -vectors_l(j,i+1));
                    i++;
                } else 
                    ret.eigenPairs[i].vector(j) = CDSComplex<T>(vectors_l(j,i),0.);

    if ( jobvr == 'V' )
        for (int i=0 ; i<matrix.ncol() ; i++) {
            if ( values_i(i) != 0. ) {
                for (int j=0 ; j<matrix.ncol() ; j++) {
                    ret.eigenPairs[i].vector(j) 
                        = CDSComplex<T>(vectors_r(j,i), vectors_r(j,i+1));
                    ret.eigenPairs[i+1].vector(j) 
                        = CDSComplex<T>(vectors_r(j,i), -vectors_r(j,i+1));
                }
                i++;
            } else 
                for (int j=0 ; j<matrix.ncol() ; j++)
                    ret.eigenPairs[i].vector(j) = CDSComplex<T>(vectors_r(j,i),0.);
        }

    return ret;   
}

template<class MATRIX>
EigenResults<SymmetricMatrix<typename MATRIX::ElementType> >
callEigen(const MATRIX&              matrix,
	      const char                 jobz,
		  EigenResults<SymmetricMatrix<typename MATRIX::ElementType> > 
                                     ret)
{
    typedef typename MATRIX::ElementType T;
    MATRIX trashedByLapack = matrix;
    ret.work.resize(3*matrix.ncol());

    CDSVector<T> values(matrix.nrow());
    CDSMatrix<T> vectors(matrix.nrow(),matrix.ncol());

    callEigenSymmetric(jobz,'L',matrix.ncol(),trashedByLapack.pointer(),
                       values.pointer(),
                       vectors.pointer(),vectors.nrow(),
                       ret.work.pointer(),
                       ret.info);

    if (ret.info>0) {
        OStringStream s;
        s << "eigen: " << ret.info 
          << " off-diagonal elements on an intermediate tridiagonal form "
          << "did not converge to zero." << ends;
        cerr << s.str() << endl;
        throw CDS::exception(s.str());
    } else if (ret.info<0) {
        OStringStream s;
        s << "eigen: argument " << -ret.info 
          << " has an illegal value." << ends;
        cerr << s.str() << endl;
        throw CDS::exception(s.str());
    }

    ret.eigenPairs.resize(matrix.ncol());
    for (int i=0 ; i<matrix.ncol() ; i++) {
        ret.eigenPairs[i].value = values(i);
        ret.eigenPairs[i].vector.resize( matrix.ncol() );
        for (int j=0 ; j<matrix.ncol() ; j++)
            ret.eigenPairs[i].vector(j) = vectors(j,i);
    }
    return ret;
}


template<class Number>           // a placeholder template function:
static void                      // must be specialized.
callEigenSymmetric(const char   &JOBZ,
		           const char   &UPLO,
		           const int    &N,
			       Number        A[],
			       Number        W[],
			       Number        Z[],
		           const int    &LDZ,
			       Number        WORK[],
			       int          &INFO ) 
{ throw CDS::exception("callEigenSymmetric not defined"); }

template<class NUMBER>
inline void 
callEigenFull(const char   &JOBVL, 
	          const char   &JOBVR, 
	          const int    &N,
		      NUMBER        A[],
	          const int    &LDA, 
		      NUMBER        WR[],
		      NUMBER        WI[],
		      NUMBER        VL[],
	          const int    &LDVL,
		      NUMBER        VR[],
	          const int    &LDVR,
		      NUMBER        WORK[],
	          const int    &LWORK,
		      int          &INFO )
{ throw CDS::exception("callEigenFull not defined"); }

namespace MatrixTools {
template<class MATRIX>
MatrixTools::SVDResults<typename MATRIX::ElementType>
svd(const MATRIX&   matrix,
    const char      jobu,
    const char      jobvt,
	MatrixTools::SVDResults<typename MATRIX::ElementType> 
                    ret)
{
    typedef typename MATRIX::ElementType T;

    int m = matrix.nrow();
    int n = matrix.ncol();

    int minWorkSize = CDSMath::max(3*CDSMath::min(m,n)+
                      CDSMath::max(m,n),5*CDSMath::min(m,n));
    if (ret.work.size() < minWorkSize)
        ret.work.resize( minWorkSize );

    ret.sigma.resize( CDSMath::min(m,n) );

    switch ( jobu ) {
    case 'A' : ret.u.resize(m,m); break;
    case 'S' : ret.u.resize(m,CDSMath::min(m,n)); break;
    case 'N' : ret.u.resize(1,0); break;
    //case 'O' : ret.u.resize(1,0); break;
    default: throw CDS::exception(CDSString("callSVD: bad value for jobu argument: ") +
                                  jobu);
    }

    switch ( jobvt ) {
    case 'A' : ret.vT.resize(n,n); break;
    case 'S' : ret.vT.resize(CDSMath::min(m,n),n); break;
    case 'N' : ret.vT.resize(1,0); break;
    //case 'O' : ret.vT.resize(1,0); break;
    default: throw CDS::exception(CDSString("callSVD: bad value for jobvt argument: ") +
                                  jobvt);
    }

    ret.info=-100;

    CDSMatrix<T> trashedBySVD = matrix;

    trashedBySVD(0,0) =  trashedBySVD(0,0);  //force copy

    callSVD(jobu,jobvt,m,n,trashedBySVD.pointer(),m,
            ret.sigma.pointer(),
            ret.u.pointer(),ret.u.nrow(),ret.vT.pointer(),ret.vT.nrow(),
            ret.work.pointer(),ret.work.size(),ret.info);
    return ret;
}
}

//
// perform generalized orthogonal transform on m: S * m * transpose(S)
//
//template<class MATRIX> 
//MATRIX
//MatrixTools::orthoTransform(const MATRIX& m,
//		    const MATRIX& S)
//  // return S * m * transpose(S);
//{
// typedef typename MATRIX::ElementType T;
// 
// assert( m.nrow() == m.ncol() );
// assert( m.nrow() == S.ncol() );
//
// CDSMatrix<T> dum(S.nrow(),S.ncol(),(T)0);
// for (int i=0 ; i<dum.nrow() ; i++) 
//   for (int k=0 ; k<dum.ncol() ; k++) 
//     for (int j=0 ; j<m.nrow() ; j++) 
//       dum(i,k) += S(i,j) * m(j,k);
//
// MATRIX ret;
// ret.resize(S.nrow(),S.nrow());
// for (int i=0 ; i<ret.nrow() ; i++) 
//   for (int j=0 ; j<ret.ncol() ; j++) 
//     ret(i,j) = (T)0;
//
// for (int i=0 ; i<ret.nrow() ; i++) 
//   for (int k=0 ; k<ret.ncol() ; k++) 
//     for (int j=0 ; j<m.nrow() ; j++) 
//       ret(i,k) += dum(i,j) * S(k,j);
// 
// return ret;
//} /* orthoTransform */

template<> inline void 
callEigenSymmetric<double>
	(const char   &JOBZ,
     const char   &UPLO,
	 const int    &N,
	 double        A[],
     double        W[],
	 double        Z[],
	 const int    &LDZ,
	 double        WORK[],
	 int          &INFO )
{
    FORTRAN(dspev,DSPEV)(JOBZ FORTRAN_STRLEN_FOLLOWS_CALL(1),UPLO FORTRAN_STRLEN_FOLLOWS_CALL(1),
	        N,A,W,Z,LDZ,WORK,INFO
		    FORTRAN_STRLEN_ATEND_CALL(1) FORTRAN_STRLEN_ATEND_CALL(1));
}

template<> inline void 
callEigenFull<double>
             (const char   &JOBVL, 
		      const char   &JOBVR, 
		      const int    &N,
			  double        A[],
		      const int    &LDA, 
			  double        WR[],
			  double        WI[],
			  double        VL[],
		      const int    &LDVL,
			  double        VR[],
		      const int    &LDVR,
			  double        WORK[],
		      const int    &LWORK,
			  int          &INFO )
{
    FORTRAN(dgeev,DGEEV)(JOBVL FORTRAN_STRLEN_FOLLOWS_CALL(1), JOBVR FORTRAN_STRLEN_FOLLOWS_CALL(1), 
	        N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO 
		    FORTRAN_STRLEN_ATEND_CALL(1) FORTRAN_STRLEN_ATEND_CALL(1) );
}
//
//template<>
//inline void 
//callEigen<GenericMatrix<double> >(const char   &JOBZ,
//			   const char   &UPLO,
//			   const int    &N,
//				 double A[],
//				 double W[],
//				 double Z[],
//			   const int    &LDZ,
//				 double WORK[],
//				 int    &INFO )
//{
// FORTRAN(dgeev)(JOBZ,UPLO,N,A,W,Z,LDZ,WORK,INFO);
//}

template<> inline void 
callSVD<double>(const char   &JOBU, 
		        const char   &JOBVT,
		        const int    &M, 
		        const int    &N, 
		        double        A[],
		        const int    &LDA,
		        double        S[],
		        double        U[],
		        const int    &LDU, 
		        double        VT[], 
		        const int    &LDVT, 
		        double        WORK[],
		        const int    &LWORK, 
		        int          &INFO )
{
    FORTRAN(dgesvd,DGESVD)
	    (JOBU  FORTRAN_STRLEN_FOLLOWS_CALL(1),
	     JOBVT FORTRAN_STRLEN_FOLLOWS_CALL(1),
	     M,N,A,LDA,S,U,LDU,VT,LDVT,WORK,LWORK,INFO
	     FORTRAN_STRLEN_ATEND_CALL(1) 
	     FORTRAN_STRLEN_ATEND_CALL(1));
}

template<> inline void 
callTRI<double>(const int   &N,
		        double       A[],
		        const int   &LDA,
		        const int    IPIV[], 
		        double       WORK[], 
		        const int   &LWORK, 
		        int         &INFO )
{
    FORTRAN(dgetri,DGETRI)(N,A,LDA,IPIV, WORK, LWORK, INFO);
}
template<> inline void 
callTRF<double>(const int   &M,
	            const int   &N, 
	            double       A[],
	            const int   &LDA, 
	            int          IPIV[], 
	            int         &INFO )
{
    FORTRAN(dgetrf,DGETRF)(M,N, A,LDA, IPIV, INFO);
}

template<> inline void
callSPTRI<double>(const char    &UPLO,
		          const int     &size,
			      double         A[],
			      int            IPIV[],
			      double         WORK[], 
			      int           &INFO )
{
    FORTRAN(dsptri,DSPTRI)
	    (UPLO FORTRAN_STRLEN_FOLLOWS_CALL(1),
         size,A,IPIV,WORK,INFO 
	     FORTRAN_STRLEN_ATEND_CALL(1));
}

template<> inline void
callSPTRF<double>(const char &UPLO,
		          const int  &size,
			      double      A[],
			      int         IPIV[], 
			      int        &INFO )
{
    FORTRAN(dsptrf,DSPTRF)
	    (UPLO FORTRAN_STRLEN_FOLLOWS_CALL(1),
        size,A,IPIV,INFO
	    FORTRAN_STRLEN_ATEND_CALL(1));
}


//template<class Matrix1, class Matrix2>
//CDSMatrix<typename Matrix1::ElementType>
//operator*(const Matrix1& m1,
//	    const Matrix2& m2)
//{
// assert( m1.ncol() == m2.nrow() );
//
// typedef typename Matrix1::ElementType T;
//
// CDSMatrix<T> ret(m1.nrow(),m2.ncol(),(T)0);
// for (int i=0 ; i<m1.nrow() ; i++)
//   for (int k=0 ; k<m2.ncol() ; k++)
//     for (int j=0 ; j<m1.ncol() ; j++) 
//	 r(i,k) += m1(i,j) * m2(j,k);
//
// return r;
//} /* operator*(matrix1,matrix2) */

//
// eigenvalue analysis
//

namespace MatrixTools {
template<class MATRIX>
EigenResults<typename MATRIX::MatrixType>
eigen(const MATRIX&              matrix,
      const char                 jobz,
	  EigenResults<typename MATRIX::MatrixType> 
                                 ret)
{
    typedef typename MATRIX::ElementType T;

    if (matrix.nrow()==0) return ret;

    return callEigen(matrix,jobz,ret);
}
}

#endif // __matrixToolsDefs_h__
