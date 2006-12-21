
#include "IpLapackSolverInterface.hpp"
#include "SimTKlapack.h"


namespace Ipopt
{
#ifdef IP_DEBUG
  static const Index dbg_verbosity = 0;
#endif

int *ipiv;
double *afact;
  LapackSolverInterface::LapackSolverInterface()
      :
      n(0),
      nz(0),
      a(NULL),
      irn_(NULL),
      jcn_(NULL),
      negevals_(-1)
  {
    DBG_START_METH("LapackSolverInterface::LapackSolverInterface()", dbg_verbosity);
  }


  LapackSolverInterface::~LapackSolverInterface()
  {
    DBG_START_METH("LapackSolverInterface::~LapackSolverInterface()", dbg_verbosity);
//    dmumps_c(&mumps_data); /* Terminate instance */

    delete [] a;
    delete [] irn_;
    delete [] jcn_;
  }

  void LapackSolverInterface::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {}

  bool LapackSolverInterface::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    return true;
  }

  ESymSolverStatus LapackSolverInterface::MultiSolve(bool new_matrix, const Index* ia, const Index* ja,
      Index nrhs, double* rhs_vals, bool check_NegEVals,
      Index numberOfNegEVals)
  {
    int i;
    double *atmp;
    DBG_START_METH("LapackSolverInterface::MultiSolve", dbg_verbosity);
    DBG_ASSERT(!check_NegEVals || ProvidesInertia());
    //DBG_ASSERT(initialized_);


    ESymSolverStatus retval = SYMSOLVER_SUCCESS;

    // check if a factorization has to be done
    // perform the factorization

    char fact;
    char uplo = 'L';
    const char& rfact = fact;
    const char& ruplo = uplo;
    double *x,*berr,*ferr,*work,rcond;
    int *iwork,lwork;
    int &rlwork = lwork;
    const int &ndim = n;
    const int &rnrhs = nrhs;
    int info;
    int &rinfo = info;
    

    atmp = new double[n*n];
    for(i=0;i<n*n;i++) atmp[i] = a[i];
    if (new_matrix) {
      retval = Factorization(ia, ja, check_NegEVals, numberOfNegEVals);
      if (retval == SYMSOLVER_SUCCESS)  {
         isFactored = 1;
      } else {
        DBG_PRINT((1, "FACTORIZATION FAILED!\n"));
         printf( "MultiSolve initial FACTORIZATION FAILED! retval = %d\n",retval);
         isFactored = 0;
      }

    }
      if (isFactored)  {
        retval =  Solve(ia, ja, nrhs, rhs_vals);
      } else {
         const int &ndim = n;
         const int &rnrhs = nrhs;
         double rcond = -1.0;
         double *s,*work;
         int ispec = 1;
         int info;
         int &rinfo = info;
         int *iwork,rank,nlvl,smlsiz,lwork;
         int &rlwork = lwork;
         const char *name = "DGELSD";
         const char *opts = "";
         s = new double[n];
         smlsiz = ilaenv_( &ispec, name, opts, &n, &n, &n, &n, 6, 0);
         if( smlsiz < 0 ) {
             printf("ilaenv arg# %d illegal value \n",smlsiz );
             return retval;
         }
//         nlvl = int(log2(n)/(smlsiz+1)) + 1;
         nlvl = int((log10((double)n)/log10(2.))/(smlsiz+1)) + 1;
         iwork = new int[3*n*nlvl + 11*n];
         lwork = 12*n + 2*n*smlsiz + 8*n*nlvl + n*nrhs + (smlsiz+1)*(smlsiz+1);
         work = new double[lwork];
         dgelsd_( ndim, ndim, rnrhs, atmp, ndim, rhs_vals, ndim, s, &rcond, &rank, work, 
                  rlwork, iwork, rinfo );
         
         delete [] work;
         delete [] s;
         delete [] iwork;
         if( info > 0 ) {
            printf( "dgelsd %d elements failed to converge to zero \n",info );
         } else if( info < 0 ) {
            printf( "dgelsd illegal arg #%d \n",info );
         } else {
            retval = SYMSOLVER_SUCCESS;
         }
      }
      delete [] atmp;
      return retval;  

  }


  double* LapackSolverInterface::GetValuesArrayPtr()
  {
    return a;
  }

  /** Initialize the local copy of the positions of the nonzero
      elements */
  ESymSolverStatus LapackSolverInterface::InitializeStructure(Index dim, Index nonzeros,
      const Index* ia, const Index* ja)
  {
    ESymSolverStatus retval = SYMSOLVER_SUCCESS;
    DBG_START_METH("LapackSolverInterface::InitializeStructure", dbg_verbosity);
    if (a) {
      delete [] a;
    }
    n = dim;
    nz = nonzeros;
    delete [] a;
    delete [] irn_;
    delete [] jcn_;
    a = new double[dim*dim];
    irn_ = new int[nz];
    jcn_ = new int[nz];
    for (Index i=0; i<nz; i++) {
      irn_[i] = ia[i];
      jcn_[i] = ja[i];
    }

//    dmumps_c(&mumps_data);


    return retval;
  }


  ESymSolverStatus LapackSolverInterface::Factorization(const Index* ia, const Index* ja,
      bool check_NegEVals, Index numberOfNegEVals)
  {
    DBG_START_METH("LapackSolverInterface::Factorization", dbg_verbosity);
    ESymSolverStatus retval = SYMSOLVER_SUCCESS;
    int info,lwork;
    int &rinfo = info;
    int &rlwork = lwork;
    double *w, *work,*atmp;
    const int &ndim = n;
    char jobzc = 'N';
    char uploc = 'L';
    const char &jobz = jobzc;
    const char &uplo = uploc;
    int *tmp_ipiv; 
    int i,j;

    
    ipiv = new int[n];
    tmp_ipiv = (int *)malloc( sizeof(int)*n );
/* printf("LapackSolverInterface::Factorization a=\n");
    for(i=0;i<n;i++) {
       for(j=0;j<n;j++) {
           printf(" %f" , a[i*n+j]);
       }
       printf(" \n");
    }
*/
    /* compute negative eigenvalues */

    negevals_ = 0;
    w = new double[n];
    lwork = 3*n;   // TODO get optimial value 
    work = new double[lwork];
/* create a tempory copy  of the A matrix becuase dsyev over writes it */
    atmp = new double[n*n];
    for(i=0;i<n*n;i++) atmp[i] = a[i];
    dsyev_(jobz, uplo, ndim, atmp, ndim, w, work, rlwork, rinfo,  1, 1);
//    delete [] atmp;
    if( rinfo != 0 ) {
         printf("dsyev failed info = %d\n",rinfo  );
         return(SYMSOLVER_FATAL_ERROR);
    }
    for(i=0;i<ndim;i++){
   //     printf(" eigenvlaue #%d = %f \n",i,w[i] );
        if( w[i] < 0.0 ) negevals_++;
    }
    delete [] w;
    delete [] work;
    if (check_NegEVals && (numberOfNegEVals!=negevals_)) {
printf("Factorization SYMSOLVER_WRONG_INERTIA numberOfNegEVals=%d negevals_=%d\n",numberOfNegEVals, negevals_);
       return SYMSOLVER_WRONG_INERTIA;
    }


    dgetrf_( ndim, ndim, a, ndim, ipiv, rinfo); 
    delete [] atmp;
/*
         lwork = ndim *4; // TODO compute proper blocks size
         work = new double[lwork];
         afact = new double[n*n];
         for(i=0;i<n*n;i++) afact[i] = a[i];

         dsytrf_(uplo, ndim, afact, ndim, ipiv, work, rlwork, rinfo);

         delete [] work;
*/
/*
printf("\n\nLapackSolverInterface::Factorization factored a=\n");
    for(i=0;i<n;i++) {
       for(j=0;j<n;j++) {
          printf(" %f" , a[i*n+j]);
       }
       printf(" \n");
    }
*/
       if( info > 0  ) {
         retval = SYMSOLVER_SINGULAR;
         printf(" LapackSolverInterface::Factorization dgetrf failed info=%d\n",info);
       }

    return retval;
  }

  ESymSolverStatus LapackSolverInterface::Solve(const Index* ia, const Index* ja, Index nrhs, double *b)
  {
    DBG_START_METH("LapackSolverInterface::Solve", dbg_verbosity);
    ESymSolverStatus retval = SYMSOLVER_SUCCESS;
    int info,cond,lwork,*iwork;
    int &rinfo = info;
    char transpose = 'N';
    char fact = 'F';
    char uplo = 'L';
    const char& trans = transpose;
    const char& rfact = fact;
    const char& ruplo = uplo;
    const int &ndim = n;
    const int &rnrhs = nrhs;
    const int &rlwork = lwork;
    double *ferr,*work;
    double *berr;
    int i,j;
    double *x,rcond,*tmp_b;
/*
//printf("\n\nLapackSolverInterface::Solve nrhs =%d dim=%d trans=%c\n",nrhs,ndim,trans);
printf("input rhs =");
    for(i=0;i<nrhs;i++) {
       for(j=0;j<n;j++) {
           printf(" %f" , b[i*n+j]);
       }
//       printf(" \n");
    }
printf("\n\nLapackSolverInterface::Solve a =\n");

    for(i=0;i<n;i++) {
       for(j=0;j<n;j++) {
           printf(" %f" , a[i*n+j]);
       }
       printf(" \n");
    }


printf("\n\n** Solve ipiv  =%x \n",ipiv);
for(i=0;i<ndim;i++) printf("%d ",*(ipiv+i)); printf("\n");
*/
    tmp_b = new double[n*nrhs];
    for(i=0;i<n*nrhs;i++) tmp_b[i] = b[i];
//    printf("ipiv="); for(i=0;i<n;i++)printf("%d ",*(ipiv+i)); printf("\n");
    dgetrs_( trans, ndim, rnrhs, a, ndim, ipiv, b, ndim, rinfo, 1); 
/*
printf("  output rhs =");
    for(i=0;i<nrhs;i++) {
       for(j=0;j<n;j++) {
           printf(" %f" , b[i*n+j]);
       }
    }
       printf(" \n");
*/

//    delete [] ipiv;
/*
    x = new double[n*nrhs];
    berr = new double[nrhs];
    ferr = new double[nrhs];
    lwork = 32*n;
    work = new double[lwork];
    iwork = new int[n];
    dsysvx_( rfact, ruplo, ndim, rnrhs, a, ndim, afact, ndim, ipiv, 
             b, ndim, x, ndim, &rcond, ferr, berr, work, rlwork, iwork, rinfo, 1, 1); 

    for(i=0;i<n*nrhs;i++) b[i] = x[i];
    
    delete [] iwork;
    delete [] x;
    delete [] work;
    delete [] berr;
    delete [] ferr;
*/
/*
printf("\n\nLapackSolverInterface::Solve solve rhs =\n");
    for(i=0;i<nrhs;i++) {
       for(j=0;j<n;j++) {
           printf(" %f" , b[i*n+j]);

       }
       printf(" \n");
    }
*/
       if( rinfo != 0 ) {
          if( rinfo > 0 && rinfo <= n  ) {
             retval = SYMSOLVER_SINGULAR;
             printf(" LapackSolverInterface::Solve dgetrs_ singular info = %d \n,rinfo");
           } else {
              retval = SYMSOLVER_FATAL_ERROR;
              printf(" LapackSolverInterface::Solve dgetrs_ failed info = %d \n,rinfo");
           }
       } else {
          retval = SYMSOLVER_SUCCESS;
       }


    return retval;
  }

  Index LapackSolverInterface::NumberOfNegEVals() const
  {
    DBG_START_METH("LapackSolverInterface::NumberOfNegEVals", dbg_verbosity);
//printf("LapackSolverInterface::NumberOfNegEVals = %d\n",negevals_);
    DBG_ASSERT(negevals >= 0);
    return negevals_;
  }

  bool LapackSolverInterface::IncreaseQuality()
  {
    return false;
  }

}//end Ipopt namespace



