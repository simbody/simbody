#ifndef netlib_h_
#define netlib_h_

/*
//:
// \file
// \brief Header file for all exported netlib (fortran-originating) functions
// \author Peter Vanroose, KULeuven
// \date   March 2002
*/

#ifdef __cplusplus
#include "vcl_complex.h"
typedef vcl_complex<float> cmplx;
typedef vcl_complex<double> dcmplx;
typedef int logical;
extern "C" {
#else
#define cmplx complex
#define dcmplx doublecomplex
#endif

#include "vxl_config.h"
#ifndef sqrtf // for VS8
#if VXL_C_MATH_HAS_SQRTF
float sqrtf(float);
#else
# define sqrtf(f) ((float)sqrt((double)(f)))
#endif
#endif

  char * F77_aloc(int Len, const char *whence);
  void   sig_die(const char *s, int kill);
  void   exit_(int *rc);
  void   s_cat(char *lp, const char *rpp[], long rnp[], long *np, long ll);
  int    s_cmp(const char *a0, const char *b0, long la, long lb);
  void   s_copy(char *a, const char *b, long la, long lb);
  double f__cabs(double, double);

  double pow_dd(const double *x, const double *y);
  double pow_di(const double *ap, const int *bp);
  int    pow_ii(const int *ap, const int *bp);
  float  pow_ri(const float *ap, const int *bp);
  float  c_abs(const cmplx *z);
  double z_abs(const dcmplx *z);
  void   r_cnjg(cmplx *r, const cmplx *z);
  void   d_cnjg(dcmplx *r, const dcmplx *z);
  void   c_div(cmplx *c, const cmplx *a, const cmplx *b);
  void   z_div(dcmplx *c, const dcmplx *a, const dcmplx *b);
  int    i_dnnt(const double *x);
  float  r_imag(const cmplx *z);
  double d_imag(const dcmplx *a);
  double d_lg10(const double *x);
  float  r_sign(const float *a, const float *b);
  double d_sign(const double *a, const double *b);
  void   z_sqrt(dcmplx *ret_value, const dcmplx *z);
  double urand_(int *iy);
  void   xerbla_(const char *srname, int *info);

  /*: Computes singular values and vectors of an mxn matrix (double version) */
  void dsvdc_(double *x, const int* ldx, /*!< (IN) matrix, m rows, n columns, stored row-wise */
              const int* m, const int* n,
              double *singular_values, /*!< (OUT) in descending order of magnitude */
              double *errors, /*!< (OUT) superdiagonal of u^T*x*v (normally 0) */
              double *u, const int* ldu, /*!< (OUT) left singular vectors */
              double *v, const int* ldv, /*!< (OUT) right singular vectors */
              double *work,  /*!< (IN/OUT) scratch work area of length m */
              const int* job, /*!< (IN) 2-digit number. First digit refers to u; 0 = do not compute, 1 = all m; 2 = only min(m,n) */
              int *info); /*!< (OUT) singular values [info] and up are correct */

  /*: Computes singular values and vectors of an mxn matrix (float version) */
  void ssvdc_(float* x, const int* ldx, const int* m, const int* n,
              float* s, float* e, float* u, const int* ldu, float* v, const int* ldv,
              float* work, const int* job, int* info);

  /*: Computes singular values and vectors of an mxn matrix (double_complex version) */
  void zsvdc_(dcmplx* x, const int* ldx, const int* m, const int* n,
              dcmplx* s, dcmplx* e, dcmplx* u, const int* ldu, dcmplx* v, const int* ldv,
              dcmplx* work, const int* job, int* info);

  /*: Computes singular values and vectors of an mxn matrix (float_complex version) */
  void csvdc_(cmplx* x, const int* ldx, const int* m, const int* n,
              cmplx* s, cmplx* e, cmplx* u, const int* ldu, cmplx* v, const int* ldv,
              cmplx* work, const int* job, int* info);

  void sggsvd_(char *jobu, char *jobv, char *jobq, int *m, int *n, int *p, int *k, int *l, float *a, int *lda,
               float *b, int *ldb, float *alpha, float *beta, float *u, int * ldu, float *v, int *ldv, float *q, int *ldq,
               float *work, int *iwork, int *info);

  void dggsvd_(char *jobu, char *jobv, char *jobq, int *m, int *n, int *p, int *k, int *l,
               double *a, int *lda, double *b, int *ldb, double *alpha, double *beta,
               double *u, int *ldu, double *v, int *ldv, double *q, int *ldq, double *work, int *iwork, int *info);

  /*: Finds eigenvalues and eigenvectors of a general matrix */
  void rg_(int* nm, int* n, double* a, double* wr, double* wi, int* matz, double* z, int* iv1, double* fv1, int* ierr);

  /*: Computes eigenvalues and eigenvectors of nxn complex general matrix */
  void zgeev_(const char* jobvl, const char* jobvr, const int* n, dcmplx* a, const int* lda, dcmplx* w,
              dcmplx* vl, const int* ldvl, dcmplx* vr, const int* ldvr,
              dcmplx* work, int* lwork, double* rwork, int* info);

  /*: Computes eigenvalues and eigenvectors of a real symmetric matrix */
  void rs_(const int* nm, /*!< (IN) leading dimension of matrices */
           const int* n, /*!< (IN) order of the square matrix a */
           double *a_matrix, /*!< (IN) real symmetric matrix */
           double *eigenvalues,    /*!< (OUT) eigenvalues in ascending order */
           const int* want_eigenvectors, /*!< (IN) set nonzero if eigenvectors wanted */
           double *eigenvectors, /*!< (OUT) eigenvectors */
           double *workspace_1_size_n, double *workspace_2_size_n, /*!< scratch */
           int* output_error_code); /*!< (OUT) normal completion code is 0 */

  /*: Computes eigenvalues and eigenvectors of a real symmetric generalized eigenproblem  ax = lambda bx.  */
  void rsg_(const int* nm, /*!< (IN) leading dimension of matrices */
            const int* n, /*!< (IN) order of the square matrices a and b */
            double *a_matrix, /*!< (IN) real symmetric matrix */
            double *b_matrix, /*!< (IN) positive definite real symm matrix */
            double *eigenvalues,    /*!< (OUT) eigenvalues in ascending order */
            const int* want_eigenvectors, /*!< (IN) set nonzero if eigenvectors wanted */
            double *eigenvectors, /*!< (OUT) eigenvectors */
            double *workspace_1_size_n, double *workspace_2_size_n, /*!< scratch */
            int* output_error_code); /*!< (OUT) normal completion code is 0 */

  /*: Computes QR factorisation of an n x p double matrix */
  void dqrdc_(double *x, const int* ldx, /*!< (IN/OUT) matrix, n rows, p columns, stored row-wise */
              const int* n, const int* p,
              double* qraux, /*!< (OUT) further info necessary to recover R part from x */
              int *jpvt,  /*!< (IN/OUT) length p; selection of pivot columns: */
                          /*   ==0 ==> any; >0 ==> initial column; <0 ==> final */
              double *work,  /*!< (IN/OUT) scratch work area of length p */
              const int* job); /*!< (IN) if == 0, no pivoting is done */

  /*: Computes QR factorisation of an n x p float matrix */
  void sqrdc_(float* x, const int* ldx, const int* n, const int* p,
              float* qraux, int* jpvt, float* work, const int* job);

  /*: Computes QR factorisation of an n x p double_complex matrix */
  void zqrdc_(dcmplx* x, const int* ldx, const int* n, const int* p,
              dcmplx* qraux, int* jpvt, dcmplx* work, const int* job);

  /*: Computes QR factorisation of an n x p float_complex matrix */
  void cqrdc_(cmplx* x, const int* ldx, const int* n, const int* p,
              cmplx* qraux, int* jpvt, cmplx* work, const int* job);

  /*: Computes coord transf etc from QR factorisation of double matrix */
  void dqrsl_(const double* x, const int* ldx, /*!< (IN) output of dqrdc_, n x k matrix */
              const int* n, const int* k, /*!< (IN) k <= min(n,p) with n,p from dqrdc_ */
              const double* qraux, /*!< (IN) qraux output of dqrdc_ */
              const double* y, /*!< (IN) n-vector to operate on */
              double* qy,  /*!< (OUT) q*y */
              double* qty, /*!< (OUT) q^T*y (conjugate transpose if complex) */
              double* b,   /*!< (OUT) solution b of min norm_2(y - x*b) */
              double* rsd, /*!< (OUT) least squares residual y - x*b = proj of y on orth complement of columns(x) */
              double* xb,  /*!< (OUT) least squares approx of x*b = proj of y on columns(x) */
              const int* job, /*!< (IN) decimal acbde: a:compute qy; c:qty; b:qty+b; d:qty+rsd; e:qty+xb */
              int* info); /*!< non-zero if r is singular and b is set. */

  /*: Computes coord transf etc from QR factorisation of float matrix */
  void sqrsl_(const float* x, const int* ldx, const int* n, const int* k,
              const float* qraux, const float* y,
              float* qy, float* qty, float* b, float* rsd, float* xb,
              const int* job, int* info);

  /*: Computes coord transf etc from QR factorisation of double_complex matrix */
  void zqrsl_(const dcmplx* x, const int* ldx, const int* n, const int* k,
              const dcmplx* qraux, const dcmplx* y,
              dcmplx* qy, dcmplx* qty, dcmplx* b, dcmplx* rsd, dcmplx* xb,
              const int* job, int* info);

  /*: Computes coord transf etc from QR factorisation of float_complex matrix */
  void cqrsl_(const cmplx* x, const int* ldx, const int* n, const int* k,
              const cmplx* qraux, const cmplx* y,
              cmplx* qy, cmplx* qty, cmplx* b, cmplx* rsd, cmplx* xb,
              int* job, int* info);

  /*: Minimizes a function using the conjugate gradient method */
  void cg_(double* x, /*!< (IN/OUT) minimizer, length n; input = starting guess */
           double* e, /*!< (OUT) max-norm of gradient */
           int* it,   /*!< (OUT) number of iterations performed */
           double* step, /*!< (IN/OUT) step size along search direction */
           const double* tolerance_on_e, const int* max_iterations,
           const int* n, /*!< (IN) number of unknowns */
           const int* m, /*!< (IN) # iterations before calc new seach direction */
           double (*cost_function)(double*),
           void (*gradient_func)(double*,double*),
           void (*both)(double*,double*,double*),
           void (*preconditioning_func)(double*,double*),
           double *work);

  /*: Computes the float cumulative distribution function value for the chi-squared distribution */
  void chscdf_(const float* x, /*!< (IN) value where the cumulative distribution must be evaluated */
               const int* nu,  /*!< (IN) # degrees of freedom */
               float* cdf);    /*!< (OUT) the function value */

  /*: Computes the double cumulative distribution function value for the chi-squared distribution */
  void dchscdf_(double* x, int* nu, double* cdf);

  /*: Self-sorting in-place generalized prime factor (complex) double fft */
  void dgpfa_(double* a, /*!< (IN/OUT) Real part of input/output vectors */
              double* b, /*!< (IN/OUT) Imaginary part of input/output vectors */
              const double* trigs, /*!< (IN) output of dsetgfpa_ (twiddle factors) */
              const int* inc, /*!< (IN) increment within each data vector (normally 1) */
              const int* jump, /*!< (IN) increment between data vectors */
              const int* n, /*!< (IN) length of the transforms; should only have 2,3,5 as prime factors */
              const int* lot, /*!< (IN) number of transforms */
              const int* isign, /*!< (IN) forward transform: +1; backward: -1 */
              const int* npqr, /*!< (IN) 3-array with the number of factors of 2,3,5 */
              int* info); /*!< (OUT) 0 if no problems */

  /*: Self-sorting in-place generalized prime factor (complex) float fft */
  void gpfa_(float* a, float* b, const float* trigs, const int* inc, const int* jump,
             const int* n, const int* lot, const int* isign, const int* nj, int* info);

  /*: Set-up routine for dgpfa_ */
  void dsetgpfa_(double* trigs, const int* n, int* ires, int* info);

  /*: Set-up routine for gpfa_ */
  void setgpfa_(float* trigs, const int* n, int* ires, int* info);

  void gpfa2f_(float *a, float *b, const float *trigs, const int *inc, const int *jump,
               const int *n, const int *mm, const int *lot, const int *isign);
  void gpfa3f_(float *a, float *b, const float *trigs, const int *inc, const int *jump,
               const int *n, const int *mm, const int *lot, const int *isign);
  void gpfa5f_(float *a, float *b, const float *trigs, const int *inc, const int *jump,
               const int *n, const int *mm, const int *lot, const int *isign);
  void dgpfa2f_(double *a, double *b, const double *trigs, const int *inc, const int *jump,
                const int *n, const int *mm, const int *lot, const int *isign);
  void dgpfa3f_(double *a, double *b, const double *trigs, const int *inc, const int *jump,
                const int *n, const int *mm, const int *lot, const int *isign);
  void dgpfa5f_(double *a, double *b, const double *trigs, const int *inc, const int *jump,
                const int *n, const int *mm, const int *lot, const int *isign);

  /*: Finds a few eigenvalues and eigenvectors at either end of the spectrum of a large sparse symmetric matrix.  */
  void dnlaso_(void (*op)(const int* n,const int* m, const double* p, double* q),
               void (*iovect)(const int* n,const int* m, double* q, const int* j, const int* k),
               const int* n, const int* nval, const int* nfig, int* nperm,
               const int* nmval, double* val,
               const int* nmvec, double* vec,
               const int* nblock, const int* maxop, const int* maxj, double* work, int* ind, int* ierr);
  void snlaso_(void (*op)(const int* n,const int* m, const float* p, float* q),
               void (*iovect)(const int* n,const int* m, float* q, const int* j, const int* k),
               const int* n, const int* nval, const int* nfig, int* nperm,
               const int* nmval, float* val,
               const int* nmvec, float* vec,
               const int* nblock, const int* maxop, const int* maxj, float* work, int* ind, int* ierr);

  /*: Factors a symmetric positive definite matrix and estimates the condition of the matrix */
  void dpoco_(double* a, int* lda, int* n, double* rcond, double* z, int* info);

  /*: Computes determinant and inverse of a certain symmetric positive definite matrix using dpoco_, dposa_ or dqrdc_ output */
  void dpodi_(double* a, const int* lda, const int* n, double* det, const int* job);

  /*: Factors a double precision symmetric positive definite matrix */
  void dpofa_(double* a, int* lda, int* n, int* info);

  /*: Solves the symmetric positive definite system a * x = b dpoco or dpofa output */
  void dposl_(const double* a, const int* lda, const int* n, double* b);

  /*: Solves the unconstrained minimization problem min F(x1..xN) */
  void lbfgs_(int* n, int* m, double* x, double* f, double* g,
              logical * diagco, double* diag, int* iprint,
              double* eps, double* xtol, double* w, int* iflag);

  /*: Minimizes the sum of the squares of m nonlin functions in n variables */
  void lmder1_(void (*fcn)(int*,int*,double*,double*,double*,int*,int*),
               int* m, int* n, double* x,
               double* fvec, double* fjac, int* ldfjac,
               double* tol, int* info, int* ipvt, double* wa, int* lwa);
  void lmder_(void (*fcn)(int*,int*,double*,double*,double*,int*,int*),
              int *m, int *n, double *x, double *fvec, double *fjac, int *ldfjac, double *ftol, double *xtol, double *gtol,
              int *maxfev, double *diag, int *mode, double *factor, int *nprint, int *info,
              int *nfev, int *njev, int *ipvt, double *qtf, double *wa1, double *wa2, double *wa3, double *wa4);

  /*: Minimizes the sum of the squares of m nonlin functions in n variables */
  void lmdif_(void (*fcn)(int*,int*,double*,double*,int*),
              int* m, int* n, double* x, double* fvec, double* ftol,
              double* xtol, double* gtol, int* maxfev, double* epsfcn, double* diag,
              int* mode, double* factor, int* nprint, int* info, int* nfev,
              double* fjac, int* ldfjac, int* ipvt, double* qtf,
              double* wa1, double* wa2, double* wa3, double* wa4, double* errors);

  /*: Solves  A*x = b */
  void lsqr_(int* m, int* n,
             void (*aprod)(int*,int*,int*,double*,double*,int*,int*,int*,double*),
             double* damp, int* leniw, int* lenrw,
             int* iw, double* rw, double* u, double* v, double* w, double* x, double* se,
             double* atol, double* btol, double* conlim, int* itnlim, int* nout, int* istop,
             int* itn, double* anorm, double* acond, double* rnorm, double* arnorm, double* xnorm);

  /*: Finds the zeros of a real polynomial */
  void rpoly_(double* op, int* degree, double* zeror, double* zeroi, logical* fail);

  void saxpy_(const int *n, const float *sa, const float *sx, const int *incx, float *sy, const int *incy);
  void daxpy_(const int *n, const double *da, const double *dx, const int *incx, double *dy, const int *incy);
  void caxpy_(const int *n, const cmplx *ca, const cmplx *cx, const int *incx, cmplx *cy, const int *incy);
  void zaxpy_(const int *n, const dcmplx *za, const dcmplx *zx, const int *incx, dcmplx *zy, const int *incy);

  void scopy_(const int *n, const float *sx, const int *incx, float *sy, const int *incy);
  void dcopy_(const int *n, const double *dx, const int *incx, double *dy, const int *incy);
  void ccopy_(const int *n, const cmplx *cx, const int *incx, cmplx *cy, const int *incy);
  void zcopy_(const int *n, const dcmplx *zx, const int *incx, dcmplx *zy, const int *incy);

 float sdot_(const int *n, const float *sx, const int *incx, const float *sy, const int *incy);
double ddot_(const int *n, const double *dx, const int *incx, const double *dy, const int *incy);
  void cdotc_(cmplx *ret_val, const int *n, const cmplx *cx, const int *incx, const cmplx *cy, const int *incy);
  void zdotc_(dcmplx *ret_val, const int *n, const dcmplx *zx, const int *incx, const dcmplx *zy, const int *incy);
  void zdotu_(dcmplx *ret_val, const int *n, const dcmplx *zx, const int *incx, const dcmplx *zy, const int *incy);

  void sscal_(const int *n, const float *sa, float *sx, const int *incx);
  void dscal_(const int *n, const double *da, double *dx, const int *incx);
  void cscal_(const int *n, const cmplx *ca, cmplx *cx, const int *incx);
  void zscal_(const int *n, const dcmplx *za, dcmplx *zx, const int *incx);
  void zdscal_(const int *n, const double *da, dcmplx *zx, const int *incx);

  void sswap_(const int *n, float *sx, const int *incx, float *sy, const int *incy);
  void dswap_(const int *n, double *dx, const int *incx, double *dy, const int *incy);
  void cswap_(const int *n, cmplx *cx, const int *incx, cmplx *cy, const int *incy);
  void zswap_(const int *n, dcmplx *zx, const int *incx, dcmplx *zy, const int *incy);

  void dgecon_(char *norm, int *n, double *a, int *lda, double *anorm, double *rcond, double *work, int *iwork, int *info);

  void dgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k, double *alpha,
              double *a, const int *lda, double *b, const int *ldb, double *beta, double *c, const int *ldc);
  void zgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k, dcmplx *alpha,
              dcmplx *a, const int *lda, dcmplx *b, const int *ldb, dcmplx *beta, dcmplx *c, const int *ldc);

  void sgemv_(const char *trans, const int *m, const int *n, float *alpha,
              float *a, const int *lda, float *x, const int *incx, float *beta, float *y, const int *incy);
  void dgemv_(const char *trans, const int *m, const int *n, double *alpha,
              double *a, const int *lda, double *x, const int *incx, double *beta, double *y, const int *incy);
  void zgemv_(const char *trans, const int *m, const int *n, dcmplx *alpha,
              dcmplx *a, const int *lda, dcmplx *x, const int *incx, dcmplx *beta, dcmplx *y, const int *incy);

  void sgeqpf_(int *m, int *n, float *a, int *lda, int *jpvt, float *tau, float *work, int *info);
  void dgeqpf_(int *m, int *n, double *a, int *lda, int *jpvt, double *tau, double *work, int *info);
  void sgeqr2_(int *m, int *n, float *a, int *lda, float *tau, float *work, int *info);
  void dgeqr2_(int *m, int *n, double *a, int *lda, double *tau, double *work, int *info);
  void dgeqrf_(int *m, int *n, double *a, int *lda, double *tau, double *work, int *lwork, int *info);

  void sger_(const int *m, const int *n, float *alpha, float *x, const int *incx, float *y, const int *incy,
             float *a, const int *lda);
  void dger_(const int *m, const int *n, double *alpha, double *x, const int *incx, double *y, const int *incy,
             double *a, const int *lda);
  void zgerc_(const int *m, const int *n, dcmplx *alpha, dcmplx *x, const int *incx, dcmplx *y, const int *incy,
              dcmplx *a, const int *lda);
  void sgerq2_(const int *m, const int *n, float *a, const int *lda, float *tau, float *work, int *info);
  void dgerq2_(const int *m, const int *n, double *a, const int *lda, double *tau, double *work, int *info);

  void dgesc2_(int *n, double *a, int *lda, double *rhs, int *ipiv, int *jpiv, double *scale);
  void dgetc2_(int *n, double *a, int *lda, int *ipiv, int *jpiv, int *info);

  void dggbak_(const char *job, const char *side, const int *n, int *ilo, int *ihi, double *lscale, double *rscale,
               const int *m, double *v, const int *ldv, int *info);
  void zgebak_(const char *job, const char *side, const int *n, int *ilo, int *ihi, double *scale,
               const int *m, dcmplx *v, const int *ldv, int *info);
  void dggbal_(const char *job, const int *n, double *a, const int *lda, double *b, const int *ldb, int *ilo, int *ihi,
               double *lscale, double *rscale, double *work, int *info);
  void zgebal_(const char *job, const int *n, dcmplx *a, const int *lda, int *ilo, int *ihi, double *scale, int *info);

  void dgges_(const char *jobvsl, const char *jobvsr, const char *sort, logical (*delctg)(double*,double*,double*),
              int *n, double *a, int *lda, double *b, int *ldb, int *sdim, double *alphar, double *alphai, double *beta,
              double *vsl, int *ldvsl, double *vsr, int *ldvsr, double *work, int *lwork, logical *bwork, int *info);

  void dgghrd_(const char *compq, const char *compz, const int *n, int *ilo, int *ihi, double *a, const int *lda,
               double *b, const int *ldb, double *q, const int *ldq, double *z, const int *ldz, int *info);
  void zgehrd_(const int *n, int *ilo, int *ihi, dcmplx *a, const int *lda, dcmplx *tau, dcmplx *work, int *lwork, int *info);
  void zgehd2_(const int *n, const int *ilo, const int *ihi, dcmplx *a, const int *lda, dcmplx *tau, dcmplx *work, int *info);

  void sggsvp_(char *jobu, char *jobv, char *jobq, int *m, int *p, int *n, float *a, int *lda, float *b, int *ldb,
               float *tola, float *tolb, int *k, int *l, float *u, int *ldu, float *v, int *ldv, float *q, int *ldq,
               int *iwork, float *tau, float *work, int *info);
  void dggsvp_(char *jobu, char *jobv, char *jobq, int *m, int *p, int *n, double *a, int *lda, double *b, int *ldb,
               double *tola, double *tolb, int *k, int *l, double *u, int *ldu, double *v, int *ldv, double *q, int *ldq,
               int *iwork, double *tau, double *work, int *info);

  void dhgeqz_(const char *job, const char *compq, const char *compz, int *n, int *ilo, int *ihi,
               double *a, int *lda, double *b, int *ldb, double *alphar, double *alphai, double *beta,
               double *q, int *ldq, double *z, int *ldz, double *work, int *lwork, int *info);

  void dlabad_(double *small, double *large);
  void dlacon_(int *n, double *v, double *x, int *isgn, double *est, int *kase);

  void dlacpy_(const char *uplo, const int *m, const int *n, double *a, const int *lda, double *b, const int *ldb);
  void slacpy_(const char *uplo, const int *m, const int *n, float *a, const int *lda, float *b, const int *ldb);
  void zlacpy_(const char *uplo, const int *m, const int *n, dcmplx *a, const int *lda, dcmplx *b, const int *ldb);

  void dladiv_(const double *a, const double *b, const double *c, const double *d, double *p, double *q);
  void zladiv_(dcmplx *ret_val, const dcmplx *x, const dcmplx *y);

  void dlag2_(double *a, int *lda, double *b, int *ldb, double *safmin, double *scale1, double *scale2,
              double *wr1, double *wr2, double *wi);
  void slags2_(logical *upper, float *a1, float *a2, float *a3, float *b1, float *b2, float *b3,
               float *csu, float *snu, float *csv, float * snv, float *csq, float *snq);
  void dlags2_(logical *upper, double *a1, double *a2, double *a3, double *b1, double *b2, double *b3,
               double *csu, double *snu, double *csv, double *snv, double *csq, double *snq);
  void dlagv2_(double *a, int *lda, double *b, int *ldb, double *alphar, double *alphai, double *beta,
               double *csl, double *snl, double *csr, double *snr);

 float slamch_(const char *cmach);
double dlamch_(const char *cmach);

 float slange_(const char *norm, const int *m, const int *n, float *a, const int *lda, float *work);
double dlange_(const char *norm, const int *m, const int *n, double *a, const int *lda, double *work);
double zlange_(const char *norm, const int *m, const int *n, dcmplx *a, const int *lda, double *work);

double dlanhs_(const char *norm, const int *n, double *a, const int *lda, double *work);
double zlanhs_(const char *norm, const int *n, dcmplx *a, const int *lda, double *work);

  void slapll_(int *n, float *x, int *incx, float *y, int *incy, float *ssmin);
  void dlapll_(int *n, double *x, int *incx, double *y, int *incy, double *ssmin);

  void slapmt_(logical *forwrd, int *m, int *n, float *x, int *ldx, int *k);
  void dlapmt_(logical *forwrd, int *m, int *n, double *x, int *ldx, int *k);

 float slapy2_(const float *x, const float *y);
double dlapy2_(const double *x, const double *y);
double dlapy3_(const double *x, const double *y, const double *z);

  void slarf_(const char *side, const int *m, const int *n, float *v, const int *incv, const float *tau,
              float *c, const int *ldc, float *work);
  void dlarf_(const char *side, const int *m, const int *n, double *v, const int *incv, const double *tau,
              double *c, const int *ldc, double *work);
  void zlarf_(const char *side, const int *m, const int *n, dcmplx *v, const int *incv, const dcmplx *tau,
              dcmplx *c, const int *ldc, dcmplx *work);

  void dlarfb_(const char *side, const char *trans, const char *direct, const char *storev,
               const int *m, const int *n, const int *k, double *v, const int *ldv, double *t, const int *ldt,
               double *c, const int *ldc, double *work, const int *ldwork);
  void zlarfb_(const char *side, const char *trans, const char *direct, const char *storev,
               const int *m, const int *n, const int *k, dcmplx *v, const int *ldv,
               dcmplx *t, const int *ldt, dcmplx *c, const int *ldc, dcmplx *work, const int *ldwork);

  void slarfg_(const int *n, float *alpha, float *x, const int *incx, float *tau);
  void dlarfg_(const int *n, double *alpha, double *x, const int *incx, double *tau);
  void zlarfg_(const int *n, dcmplx *alpha, dcmplx *x, const int *incx, dcmplx *tau);

  void dlarft_(const char *direct, const char *storev, const int *n, const int *k,
               double *v, const int *ldv, const double *tau, double *t, const int *ldt);
  void zlarft_(const char *direct, const char *storev, const int *n, const int *k, dcmplx *v, const int *ldv,
               const dcmplx *tau, dcmplx *t, int *ldt);
  void zlarfx_(const char *side, const int *m, const int *n, dcmplx *v, dcmplx *tau, dcmplx *c, const int *ldc, dcmplx *work);

  void slartg_(float *f, float *g, float *cs, float *sn, float *r);
  void dlartg_(double *f, double *g, double *cs, double *sn, double *r);

  void slas2_(float *f, float *g, float *h, float *ssmin, float * ssmax);
  void dlas2_(double *f, double *g, double *h, double *ssmin, double *ssmax);

  void dlascl_(const char *type, const int *kl, const int *ku, double *cfrom, double *cto,
               const int *m, const int *n, double *a, const int *lda, int *info);
  void zlascl_(const char *type, const int *kl, const int *ku, double *cfrom, double *cto,
               const int *m, const int *n, dcmplx *a, const int *lda, int *info);

  void slaset_(const char *uplo, const int *m, const int *n, float *alpha, float *beta, float *a, const int *lda);
  void dlaset_(const char *uplo, const int *m, const int *n, double *alpha, double *beta, double *a, const int *lda);
  void zlaset_(const char *uplo, const int *m, const int *n, dcmplx *alpha, dcmplx *beta, dcmplx *a, const int *lda);

  void slassq_(const int *n, const float *x, const int *incx, float *scale, float *sumsq);
  void dlassq_(const int *n, const double *x, const int *incx, double *scale, double *sumsq);
  void zlassq_(const int *n, const dcmplx *x, const int *incx, double *scale, double *sumsq);

  void slasv2_(float *f, float *g, float *h, float *ssmin, float *ssmax, float *snr, float *csr, float *snl, float *csl);
  void dlasv2_(double *f, double *g, double *h, double *ssmin, double *ssmax, double *snr, double * csr, double *snl, double *csl);

  void dlaswp_(int *n, double *a, int *lda, int *k1, int *k2, int *ipiv, int *incx);

  void dlatdf_(int *ijob, int *n, double *z, int *ldz, double *rhs, double *rdsum, double *rdscal, int *ipiv, int *jpiv);

 float snrm2_(const int *n, const float *x, const int *incx);
double dnrm2_(const int* n, const double* x, const int* incx);
 float scnrm2_(const int *n, const cmplx *x, const int *incx);
double dznrm2_(const int *n, const dcmplx *x, const int *incx);
double enorm_(const int *n, const double *x);

  void sorg2r_(int *m, int *n, int *k, float *a, int *lda, float *tau, float *work, int *info);
  void dorg2r_(int *m, int *n, int *k, double *a, int *lda, double *tau, double *work, int *info);
  void dorgqr_(int *m, int *n, int *k, double *a, int *lda, double *tau, double *work, int *lwork, int *info);
  void dorgr2_(int *m, int *n, int *k, double *a, int *lda, double *tau, double *work, int *info);

  void sorm2r_(const char *side, const char *trans, const int *m, const int *n, const int *k,
               float *a, const int *lda, const float *tau, float *c, const int *ldc, float *work, int *info);
  void dorm2r_(const char* side, const char* trans, const int *m, const int *n, const int *k,
               double *a, const int *lda, const double *tau, double *c, const int*ldc, double *work, int *info);
  void dormqr_(const char *side, const char *trans, const int *m, const int *n, const int *k,
               double *a, const int *lda, double *tau, double *c, const int *ldc, double *work, int *lwork, int *info);
  void sormr2_(char *side, char *trans, int *m, int *n, int*k,float*a,int*lda,float*tau,float*c,int*ldc,float*work,int*info);
  void dormr2_(char*side, char*trans, int*m, int*n, int*k, double*a, int*lda, double*tau, double*c, int*ldc, double*work, int*info);

  void lmpar_(int *n, double *r, int *ldr, int *ipvt, double *diag, double *qtb, double *delta, double *par,
              double *x, double *sdiag, double *wa1, double *wa2);
double dpmpar_(const int *i);

  void srot_(const int *n, float *sx, const int *incx, float *sy, const int *incy, const float *c, const float *s);
  void drot_(const int *n, double *dx, const int *incx, double *dy, const int *incy, const double *c, const double *s);
  void csrot_(const int *n, cmplx *cx, const int *incx, cmplx *cy, const int *incy, const float *c, const float *s);
  void zdrot_(const int *n, dcmplx *zx, const int *incx, dcmplx *zy, const int *incy, const double *c, const double *s);
  void srotg_(float *sa, float *sb, float *c, float *s);
  void drotg_(double *da, double *db, double *c, double *s);

  void drscl_(int *n, double *sa, double *sx, int *incx);

  void dtgex2_(logical *wantq, logical *wantz, int *n, double *a, int *lda, double *b, int *ldb, double *q, int *ldq,
               double *z, int *ldz, int *j1, int *n1, int *n2, double *work, int *lwork, int *info);
  void dtgexc_(logical *wantq, logical *wantz, int *n, double *a, int *lda, double *b, int *ldb, double *q, int *ldq,
               double *z, int *ldz, int *ifst, int *ilst, double *work, int *lwork, int *info);
  void dtgsen_(int *ijob, logical *wantq, logical *wantz, logical *select, int *n, double *a, int *lda, double *b, int *ldb,
               double *alphar, double *alphai, double *beta, double *q, int *ldq, double *z, int *ldz, int *m,
               double *pl, double *pr, double *dif, double *work, int *lwork, int *iwork, int *liwork, int *info);
  void stgsja_(char *jobu, char *jobv, char *jobq, int *m, int *p, int *n, int *k, int *l,
               float *a, int *lda, float *b, int *ldb, float *tola, float *tolb, float *alpha, float * beta,
               float *u, int *ldu, float *v, int *ldv, float *q, int * ldq, float *work, int *ncycle, int *info);
  void dtgsja_(char *jobu, char *jobv, char *jobq, int *m, int *p, int *n, int *k, int *l,
               double *a, int *lda, double *b, int *ldb, double *tola, double *tolb, double *alpha, double *beta,
               double *u, int *ldu, double *v, int *ldv, double *q, int * ldq, double *work, int *ncycle, int *info);
  void dtgsy2_(char *trans, int *ijob, int *m, int *n, double *a, int *lda, double *b, int *ldb,
               double *c, int *ldc, double *d, int *ldd, double *e, int *lde, double *f, int *ldf,
               double *scale, double *rdsum, double *rdscal, int *iwork, int *pq, int *info);
  void dtgsyl_(char *trans, int *ijob, int *m, int *n,
               double *a, int *lda, double *b, int *ldb, double *c, int *ldc, double *d, int *ldd, double *e, int *lde,
               double *f, int *ldf, double *scale, double *dif, double *work, int *lwork, int *iwork, int *info);

  void trans_(float *a, const int *m, const int *n, const int *mn, int *move, int *iwrk, int *iok);
  void dtrans_(double *a, const int *m, const int *n, const int *mn, int *move, int *iwrk, int *iok);

  void dtrmm_(const char *side, const char *uplo, const char *transa, const char *diag, const int *m, const int *n,
              double *alpha, double *a, const int *lda, double *b, const int *ldb);
  void ztrmm_(const char *side, const char *uplo, const char *transa, const char *diag, const int *m, const int *n,
              dcmplx *alpha, dcmplx *a, const int *lda, dcmplx *b, const int *ldb);
  void dtrmv_(const char *uplo, const char *trans, const char *diag, const int *n,
              double *a, const int *lda, double *x, const int *incx);
  void ztrmv_(const char *uplo, const char *trans, const char *diag, const int *n,
              dcmplx *a, const int *lda, dcmplx *x, const int *incx);

  void dtrsv_(const char *uplo, const char *trans, const char *diag, const int *n,
              const double *a, const int *lda, double *x, const int *incx);
  void ztrsv_(const char *uplo, const char *trans, const char *diag, const int *n,
              const dcmplx *a, const int *lda, dcmplx *x, const int *incx);
  void dlatrs_(const char *uplo, const char *trans, const char *diag, const char *normin, const int *n,
               const double *a, const int *lda, double *x, double *scale, double *cnorm, int *info);
  void zlatrs_(const char *uplo, const char *trans, const char *diag, const char *normin, const int *n,
               const dcmplx *a, const int *lda, dcmplx *x, double *scale, double *cnorm, int *info);

 float sasum_(const int *n, const float *sx, const int *incx);
double dasum_(const int *n, const double *dx, const int *incx);
double dzasum_(const int *n, const dcmplx *x, const int *incx);

  void fdjac2_(void (*fcn)(int*,int*,double*,double*,int*),
               int *m, int *n, double *x, double *fvec, double *fjac, int *ldfjac, int *iflag, double *epsfcn, double *wa);

   int isamax_(const int *n, const float *sx, const int *incx);
   int idamax_(const int *n, const double *dx, const int *incx);
   int izamax_(const int *n, const dcmplx *zx, const int *incx);
   int izmax1_(const int *n, const dcmplx *cx, const int *incx);

   int ilaenv_(const int *ispec, const char *name, const char *opts, const int *n1, const int *n2, const int *n3, const int *n4);
logical lsame_(const char *ca, const char *cb);
double pythag_(const double *a, const double *b);

  void qrfac_(int *m, int *n, double *a, int *lda, logical *pivot, int *ipvt, int *lipvt,
              double *rdiag, double *acnorm, double *wa);
  void qrsolv_(const int *n, double *r, const int *ldr, const int *ipvt, const double *diag,
               const double *qtb, double *x, double *sdiag, double *wa);

  void tql1_(const int *n, double *d, double *e, int *ierr);
  void tql2_(const int *nm, const int *n, double *d, double *e, double *z, int *ierr);
  void tred1_(const int *nm, const int *n, double *a, double *d, double *e, double *e2);
  void tred2_(const int *nm, const int *n, const double *a, double *d, double *e, double *z);

  void zhseqr_(const char *job, const char *compz, const int *n, int *ilo, int *ihi, dcmplx *h, const int *ldh,
               dcmplx *w, dcmplx *z, const int *ldz, dcmplx *work, int *lwork, int *info);
  void zlacgv_(const int *n, dcmplx *x, const int *incx);
  void zlahqr_(const logical *wantt, const logical *wantz, const int *n, const int *ilo, const int *ihi,
               dcmplx *h, const int *ldh, dcmplx *w, int *iloz, int *ihiz, dcmplx *z, const int *ldz, int *info);
  void zlahrd_(const int *n, const int *k, const int *nb, dcmplx *a, const int *lda, dcmplx *tau,
               dcmplx *t, const int *ldt, dcmplx *y, const int *ldy);
  void ztrevc_(const char *side, const char *howmny, logical *select, const int *n, dcmplx *t, const int *ldt,
               dcmplx *vl, const int *ldvl, dcmplx *vr, const int *ldvr, const int *mm, int *m,
               dcmplx *work, double *rwork, int *info);

  void zung2r_(const int *m, const int *n, const int *k, dcmplx *a, const int *lda, const dcmplx *tau,
               dcmplx *work, int *info);
  void zungqr_(const int *m, const int *n, const int *k, dcmplx *a, const int *lda, const dcmplx *tau,
               dcmplx *work, const int *lwork, int *info);
  void zunghr_(const int *n, int *ilo, int *ihi, dcmplx *a, const int *lda, const dcmplx *tau,
               dcmplx *work, const int *lwork, int *info);

  /* ITPACK functions from dsrc2c.f.  */
  int jcg_(int *nn, int *ia, int *ja, double *a, double *rhs, double *u,
           int *iwksp, int *nw, double *wksp, int *iparm, double *rparm, int *ierr);
  int jsi_(int *nn, int *ia, int *ja, double *a, double *rhs, double *u,
           int *iwksp, int *nw, double *wksp, int *iparm, double *rparm, int *ierr);
  int sor_(int *nn, int *ia, int *ja, double *a, double *rhs, double *u,
           int *iwksp, int *nw, double *wksp, int *iparm, double *rparm, int *ierr);
  int ssorcg_(int *nn, int *ia, int *ja, double *a, double *rhs, double *u,
              int *iwksp, int *nw, double *wksp, int *iparm, double *rparm, int *ierr);
  int ssorsi_(int *nn, int *ia, int *ja, double *a, double *rhs, double *u,
              int *iwksp, int *nw, double *wksp, int *iparm, double *rparm, int *ierr);
  int rscg_(int *nn, int *ia, int *ja, double *a, double *rhs, double *u,
            int *iwksp, int *nw, double *wksp, int *iparm, double *rparm, int *ierr);
  int rssi_(int *nn, int *ia, int *ja, double *a, double *rhs, double *u,
            int *iwksp, int *nw, double *wksp, int *iparm, double *rparm, int *ierr);
  int itjcg_(int *nn, int *ia, int *ja, double *a, double *u, double *u1,
             double *d__, double *d1, double *dtwd, double *tri);
  int itjsi_(int *nn, int *ia, int *ja, double *a, double *rhs, double *u, double *u1, double *d__, int *icnt);
  int itsor_(int *nn, int *ia, int *ja, double *a, double *rhs, double *u, double *wk);
  int itsrcg_(int *nn, int *ia, int *ja, double *a, double *rhs, double *u, double *u1,
              double *c__, double *c1, double *d__, double *dl, double *wk, double *tri);
  int itsrsi_(int *nn, int *ia, int *ja, double *a, double *rhs, double *u, double *u1,
              double *c__, double *d__, double *ctwd, double *wk);
  int itrscg_(int *n, int *nnb, int *ia, int *ja, double *a, double *ub, double *ub1,
              double *db, double *db1, double *wb, double *tri);
  int itrssi_(int *n, int *nnb, int *ia, int *ja, double *a, double *rhs, double *ub, double *ub1, double *db);
  int bisrch_(int *n, int *k, int *l);
  double cheby_(double *qa, double *qt, double *rrr, int *ip, double *cme, double *sme);
  int chgcon_(double *tri, double *gamold, double *rhoold, int *ibmth);
  int chgsi_(double *dtnrm, int *ibmth);
  logical chgsme_(double *oldnrm, int *icnt);
  int itpackdaxpy_(int *n, double *da, double *dx, int *incx, double *dy, int *incy);
  int itpackdcopy_(int *n, double *dx, int *incx, double *dy, int *incy);
  double itpackddot_(int *n, double *dx, int *incx, double *dy, int *incy);
  double determ_(int *n, double *tri, double *xlmda);
  int dfault_(int *iparm, double *rparm);
  int echall_(int *nn, int *ia, int *ja, double *a, double *rhs, int *iparm, double *rparm, int *icall);
  int echout_(int *iparm, double *rparm, int *imthd);
  double eigvns_(int *n, double *tri, double *d__, double *e2, int *ier);
  double eigvss_(int *n, double *tri, double *start, double *zeta, int *itmax, int *ier);
  int eqrt1s_(double *d__, double *e2, int *nn, int *m, int *isw, int *ierr);
  int ipstr_(double *omega);
  int iterm_(int *nn, double *a, double *u, double *wk, int *imthdd);
  int ivfill_(int *n, int *iv, int *ival);
  int omeg_(double *dnrm, int *iflag);
  logical omgchg_(int *ndummy);
  logical omgstr_(int *ndummy);
  int parcon_(double *dtnrm, double *c1, double *c2, double *c3, double *c4, double *gamold, double *rhotmp, int *ibmth);
  int parsi_(double *c1, double *c2, double *c3, int *ibmth);
  double pbeta_(int *nn, int *ia, int *ja, double *a, double *v, double *w1, double *w2);
  int pbsor_(int *nn, int *ia, int *ja, double *a, double *u, double *rhs);
  int permat_(int *nn, int *ia, int *ja, double *a, int *p, int *newia, int *isym, int *level, int *nout, int *ierr);
  int perror_(int *nn, int *ia, int *ja, double *a, double *rhs, double *u, double *w, double *digtt1, double *digtt2, int *idgtts);
  int pervec_(int *n, double *v, int *p);
  int pfsor_(int *nn, int *ia, int *ja, double *a, double *u, double *rhs);
  int pfsor1_(int *nn, int *ia, int *ja, double *a, double *u, double *rhs);
  int pjac_(int *nn, int *ia, int *ja, double *a, double *u, double *rhs);
  int pmult_(int *nn, int *ia, int *ja, double *a, double *u, double *w);
  int prbndx_(int *nn, int *nblack, int *ia, int *ja, int *p, int *ip, int *level, int *nout, int *ier);
  int prsblk_(int *nnb, int *nnr, int *ia, int *ja, double *a, double *ur, double *vb);
  int prsred_(int *nnb, int *nnr, int *ia, int *ja, double *a, double *ub, double *vr);
  int pssor1_(int *nn, int *ia, int *ja, double *a, double *u, double *rhs, double *fr, double *br);
  int pstop_(int *n, double *u, double *dnrm, double *ccon, int *iflag, logical *q1);
  double pvtbv_(int *n, int *ia, int *ja, double *a, double *v);
  int qsort_(int *nn, int *key, double *data, int *error);
  int sbagn_(int *n, int *nz, int *ia, int *ja, double *a, int *iwork, int *levell, int *noutt, int *ierr);
  int sbelm_(int *nn, int *ia, int *ja, double *a, double *rhs, int *iw, double *rw,
             double *tol, int *isym, int *level, int *nout, int *ier);
  int sbend_(int *n, int *nz, int *ia, int *ja, double *a, int *iwork);
  int sbini_(int *n, int *nz, int *ia, int *ja, double *a, int *iwork);
  int sbsij_(int *n, int *nz, int *ia, int *ja, double *a, int *iwork, int *ii, int *jj,
             double *vall, int *mode, int *levell, int *noutt, int *ierr);
  int scal_(int *nn, int *ia, int *ja, double *a, double *rhs, double *u, double *d__, int *level, int *nout, int *ier);
  int sum3_(int *n, double *c1, double *x1, double *c2, double *x2, double *c3, double *x3);
  double tau_(int *ii);
  double timer_(float *timdmy);
  logical tstchg_(int *ibmth);
  int unscal_(int *n, int *ia, int *ja, double *a, double *rhs, double *u, double *d__);
  int vevmw_(int *n, double *v, double *w);
  int vevpw_(int *n, double *v, double *w);
  int vfill_(int *n, double *v, double *val);
  int vout_(int *n, double *v, int *iswt, int *noutt);
  int wevmw_(int *n, double *v, double *w);
  int zbrent_(int *n, double *tri, double *eps, int *nsig, double *aa, double *bb, int *maxfnn, int *ier);

  int trapru_(double f(double*), double *a, double *b, int *m, double *trule);
  int simpru_(double f(double*), double *a, double *b, int *m, double *srule);
  /*: computes integral; input: f=integrand a,b=endpoints tol=tolerance; output: errbdd=error_estimation m=substates */
  int adaptquad_(double f(double*), double*a, double*b, double*tol, double*srmat, double*integral, double*errbdd, int*m, int*state);

#ifdef __cplusplus
}
#endif

#endif /* netlib_h_ */
