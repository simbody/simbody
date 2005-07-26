
#ifndef __lapack_hh__
#define __lapack_hh__

/*c headers for lapack*/
#include <sthead.h>
//#include <Complex.h>


#ifdef __cplusplus
extern "C" {
extern void STDCALL FORTRAN(dgetri,DGETRI)(const int &N,
		           double A[],
		     const int &LDA,
		     const int IPIV[], 
		           double WORK[], 
		     const int &LWORK, 
		           int &INFO );
extern void STDCALL FORTRAN(dgetrf,DGETRF)(const int &M,
		     const int &N, 
		           double A[],
		     const int &LDA, 
		           int IPIV[], 
		           int &INFO );

void STDCALL FORTRAN(dgeev,DGEEV)
		   (const char   &JOBVL FORTRAN_STRLEN_FOLLOWS_DECL, 
		    const char   &JOBVR FORTRAN_STRLEN_FOLLOWS_DECL, 
		    const int    &N,
		          double A[],
		    const int    &LDA, 
		          double WR[],
		          double WI[],
		          double VL[],
		    const int    &LDVL,
		          double VR[],
		    const int    &LDVR,
		          double WORK[],
		    const int    &LWORK,
		          int    &INFO  
	        FORTRAN_STRLEN_ATEND_DECL FORTRAN_STRLEN_ATEND_DECL);

void STDCALL FORTRAN(dsyev,DSYEV)
           (const char   &JOBZ FORTRAN_STRLEN_FOLLOWS_DECL, 
		    const char   &UPLO FORTRAN_STRLEN_FOLLOWS_DECL, 
		    const int    &N,
		          double A[],
		    const int    &LDA, 
		          double W[],
		          double WORK[],
		    const int    &LWORK,
		          int    &INFO  
	        FORTRAN_STRLEN_ATEND_DECL FORTRAN_STRLEN_ATEND_DECL);

void STDCALL FORTRAN(dspev,DSPEV)
           (const char   &JOBZ FORTRAN_STRLEN_FOLLOWS_DECL, 
		    const char   &UPLO FORTRAN_STRLEN_FOLLOWS_DECL, 
		    const int    &N,
		          double A[],
		          double W[],
		          double Z[],
		    const int    &LDZ,
		          double WORK[],
		          int    &INFO  
	        FORTRAN_STRLEN_ATEND_DECL FORTRAN_STRLEN_ATEND_DECL);

void STDCALL FORTRAN(dsptri,DSPTRI)
            (const char& UPLO FORTRAN_STRLEN_FOLLOWS_DECL,
		     const int& size,
                   double A[],
                   int    IPIV[],
                   double WORK[], 
                   int &INFO  
	         FORTRAN_STRLEN_ATEND_DECL);

void STDCALL FORTRAN(dsptrf,DSPTRF)
            (const char& UPLO FORTRAN_STRLEN_FOLLOWS_DECL,
		     const int& size,
                   double A[],
                   int IPIV[], 
                   int &INFO  
	         FORTRAN_STRLEN_ATEND_DECL);

void STDCALL FORTRAN(dsyevx,DSYEVX)(const char &JOBZ,
		     const char &RANGE,
		     const char &UPLO,
		     const int  &N,
		     double A[],
		     const int &LDA,
		     const double &VL,
		     const double &VU,
		     const int &IL,
		     const int &IU,
		     const double &ABSTOL,
		     int &M,
		     double W[],
		     double Z[],
		     const int &LDZ,
		     double WORK[],
		     const int &LWORK,
		     int IWORK[],
		     int IFAIL[],
		     int &INFO );
void STDCALL FORTRAN(dgelss,DGELSS)(int &M,
	     const int &N,
	     const int &NRHS,
	     double A[],
	     const int &LDA,
	     double B[],
	     const int &LDB,
	     double S[],
	     const double &RCOND,
	     int &RANK,
	     double WORK[],
	     const int &LWORK,
	     int &INFO );
void STDCALL FORTRAN(dgesv,DGESV)(int *n,
            int *nrhs,
            double A[],
            int *lda,
            int ipiv[],
            double B[],
            int *ldb,
            int *info);
void STDCALL FORTRAN(dgesvd,DGESVD)
            (const char   &JOBU  FORTRAN_STRLEN_FOLLOWS_DECL, 
		     const char   &JOBVT FORTRAN_STRLEN_FOLLOWS_DECL,
		     const int    &M, 
		     const int    &N, 
		           double  A[],
		     const int    &LDA,
		           double  S[],
		           double  U[],
		     const int    &LDU, 
		           double  VT[], 
		     const int    &LDVT, 
		           double WORK[],
		     const int    &LWORK, 
		           int    &INFO
		     FORTRAN_STRLEN_ATEND_DECL
			 FORTRAN_STRLEN_ATEND_DECL);

//  void FORTRAN(zgeev)(const char    &JOBVL,
//		      const char    &JOBVR,
//		      const int     &N,
//		            Complex  A[],
//		      const int     &LDA,
//		            Complex  W[],
//		            Complex  VL[],
//		      const int     &LDVL,
//		            Complex  VR[], 
//		      const int     &LDVR,
//		            Complex  WORK[],
//		      const int     &LWORK,
//		            double   RWORK[],
//		            int     &INFO );
//void FORTRAN(zgesv)(const int &n,
//	    const int &nrhs,
//	    Complex A[],
//	    const int &lda,
//	    int ipiv[],
//	    Complex B[],
//	    const int &ldb,
//	    int &info);
//void FORTRAN(zgesvd)(char *JOBU,
//	     char *JOBVT,
//	     int *M,
//	     int *N,
//	     Complex A[],
//	     int *LDA,
//	     double S[],
//	     Complex U[],
//	     int *LDU,
//	     Complex VT[],
//	     int *LDVT,
//	     Complex WORK[],
//	     int *LWORK,
//	     double RWORK[],
//	     int *INFO );
//void FORTRAN(zgelss)(int *M,
//	     int *N,
//	     int *NRHS,
//	     Complex A[],
//	     int *LDA,
//	     Complex B[],
//	     int *LDB,
//	     double S[],
//	     double *RCOND,
//	     int *RANK,
//	     Complex WORK[],
//	     int *LWORK,
//	     double RWORK[],
//	     int *INFO );
//void FORTRAN(zgelsx)(int *M,
//	     int *N,
//	     int *NRHS,
//	     Complex A[],
//	     int *LDA,
//	     Complex B[],
//	     int *LDB,
//	     int JPVT[],
//	     double *RCOND,
//	     int *RANK,
//	     Complex WORK[],
//	     double RWORK[],
//	     int *INFO );
}
#endif

#endif /* __lapack_hh__ */
