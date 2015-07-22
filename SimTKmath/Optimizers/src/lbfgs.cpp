/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-12 Stanford University and the Authors.        *
 * Authors: Jack Middleton                                                    *
 * Contributors:                                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

#include "simmath/internal/common.h"


#include "LBFGSOptimizer.h"
#include "../src/Simmath_f2c.h"

#include <iostream>
#include <cmath>

#define NUMBER_OF_CORRECTIONS 5

#if SimTK_DEFAULT_PRECISION==1 // float
#define DAXPY   saxpy_
#define DDOT    sdot_
#else // double
#define DAXPY   daxpy_
#define DDOT    ddot_
#endif

using SimTK::Real;

struct lb3_1_ {
/*
C    GTOL is a DOUBLE PRECISION variable with default value 0.9, which
C        controls the accuracy of the line search routine MCSRCH. If the
C        function and gradient evaluations are inexpensive with respect
C        to the cost of the iteration (which is sometimes the case when
C        solving very large problems) it may be advantageous to set GTOL
C        to a small value. A typical small value is 0.1.  Restriction:
C        GTOL should be greater than 1.D-04.
C
C    STPMIN and STPMAX are non-negative DOUBLE PRECISION variables which
C        specify lower and upper bounds for the step in the line search.
C        Their default values are 1.D-20 and 1.D+20, respectively. These
C        values need not be modified unless the exponents are too large
C        for the machine being used, or unless the problem is extremely
C        badly scaled (in which case the exponents should be increased).
*/
  int mp, lp; /* Fortran i/o stuff.  Unused here. */
  Real gtol, stpmin, stpmax;
  Real stpawf; /* line search default step length, added by awf */
};

/*#define lb3_1 (*(struct lb3_1_ *) &lb3_)*/
#define lb3_1 lb3_


using std::cout;
using std::endl;

/* Initialized data */
const struct lb3_1_ lb3_1 = {0, 0, Real(.9), Real(1e-20), Real(1e20), Real(1)};

/* Table of constant values */
static const integer c__1 = 1;

#include "SimTKlapack.h"


static void mcstep_(Real *stx, Real *fx, Real *dx, Real *sty, Real *fy, Real *dy,
                    Real *stp, Real *fp, Real *dp, bool *brackt,
                    Real *stpmin, Real *stpmax, integer *info);
void lb1_(int *iprint, int *iter, int *nfun, Real *gnorm, int *n, int *m,
          Real *x, Real *f, Real *g, Real *stp, bool *finish);
void lbptf_(char* msg);
void lbp1d_(char* msg, int* i);


/*    ----------------------------------------------------------------------*/
/*     This file contains the LBFGS algorithm and supporting routines */

/*     **************** */
/*     LBFGS SUBROUTINE */
/*     **************** */



void SimTK::LBFGSOptimizer::lbfgs_
   (int n, int m, SimTK::Real *x, SimTK::Real *f,
    int *iprint, SimTK::Real *eps, SimTK::Real *xtol)
{

    /* System generated locals */
    Real d__1;

    /* Local variables */
    Real beta;
    integer inmc;
    integer iscn, nfev, iycn, iter;
    Real ftol;
    integer nfun, ispt, iypt;
    integer bound;
    Real gnorm;
    integer point;
    integer cp;
    Real sq, yr, ys;
    Real yy;
    integer maxfev;
    integer npt;
    Real stp, stp1;

/*        LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION */
/*                          JORGE NOCEDAL */
/*                        *** July 1990 *** */

/*     This subroutine solves the unconstrained minimization problem */

/*                      min F(x),    x= (x1,x2,...,xN), */

/*      using the limited memory BFGS method. The routine is especially */
/*      effective on problems involving a large number of variables. In */
/*      a typical iteration of this method an approximation Hk to the */
/*      inverse of the Hessian is obtained by applying M BFGS updates to */
/*      a diagonal matrix Hk0, using information from the previous M steps.  */
/*      The user specifies the number M, which determines the amount of */
/*      storage required by the routine. The user may also provide the */
/*      diagonal matrices Hk0 if not satisfied with the default choice. */
/*      The algorithm is described in "On the limited memory BFGS method */
/*      for large scale optimization", by D. Liu and J. Nocedal, */
/*      Mathematical Programming B 45 (1989) 503-528. */

/*      The user is required to calculate the function value F and its */
/*      gradient G. In order to allow the user complete control over */
/*      these computations, reverse  communication is used. The routine */
/*      must be called repeatedly under the control of the parameter */
/*      IFLAG. */

/*      The steplength is determined at each iteration by means of the */
/*      line search routine MCVSRCH, which is a slight modification of */
/*      the routine CSRCH written by More' and Thuente. */

/*      The calling statement is */

/*          CALL LBFGS(N,M,X,F,G,DIAGCO,DIAG,IPRINT,EPS,XTOL,W,IFLAG) */

/*      where */

/*     N       is an INTEGER variable that must be set by the user to the */
/*             number of variables. It is not altered by the routine. */
/*             Restriction: N>0. */

/*     M       is an INTEGER variable that must be set by the user to */
/*             the number of corrections used in the BFGS update. It */
/*             is not altered by the routine. Values of M less than 3 are */
/*             not recommended; large values of M will result in excessive */
/*             computing time. 3<= M <=7 is recommended. Restriction: M>0.  */

/*     X       is a DOUBLE PRECISION array of length N. On initial entry */
/*             it must be set by the user to the values of the initial */
/*             estimate of the solution vector. On exit with IFLAG=0, it */
/*             contains the values of the variables at the best point */
/*             found (usually a solution). */

/*     F       is a DOUBLE PRECISION variable. Before initial entry and on */
/*             a re-entry with IFLAG=1, it must be set by the user to */
/*             contain the value of the function F at the point X. */

/*     IPRINT  is an INTEGER array of length two which must be set by the */
/*             user. */

/*             IPRINT(1) specifies the frequency of the output: */
/*                IPRINT(1) < 0 : no output is generated, */
/*                IPRINT(1) = 0 : output only at first and last iteration, */
/*                IPRINT(1) > 0 : output every IPRINT(1) iterations. */

/*             IPRINT(2) specifies the type of output generated: */
/*                IPRINT(2) = 0 : iteration count, number of function */
/*                                evaluations, function value, norm of the */
/*                                gradient, and steplength, */
/*                IPRINT(2) = 1 : same as IPRINT(2)=0, plus vector of */
/*                                variables and  gradient vector at the */
/*                                initial point, */
/*                IPRINT(2) = 2 : same as IPRINT(2)=1, plus vector of */
/*                                variables, */
/*                IPRINT(2) = 3 : same as IPRINT(2)=2, plus gradient vector.*/


/*    EPS     is a positive DOUBLE PRECISION variable that must be set by */
/*            the user, and determines the accuracy with which the solution*/
/*            is to be found. The subroutine terminates when */

    //sherm 100303: this criteria makes no sense to me. It looks like the X
    // is on the wrong side.
/*                         ||G|| < EPS max(1,||X||), */

/*            where ||.|| denotes the Euclidean norm. */

/*    XTOL    is a  positive DOUBLE PRECISION variable that must be set by */
/*            the user to an estimate of the machine precision (e.g. */
/*            10**(-16) on a SUN station 3/60). The line search routine will*/
/*            terminate if the relative width of the interval of uncertainty*/
/*            is less than XTOL. */


/*    ON THE DRIVER: */

/*    The program that calls LBFGS must contain the declaration: */

/*                       EXTERNAL LB2 */

/*    LB2 is a BLOCK DATA that defines the default values of several */
/*    parameters described in the COMMON section. */

/*    COMMON: */

/*     The subroutine contains one common area, which the user may wish to */
/*    reference: */

/* awf added stpawf */

/*    MP  is an INTEGER variable with default value 6. It is used as the */
/*        unit number for the printing of the monitoring information */
/*        controlled by IPRINT. */

/*    LP  is an INTEGER variable with default value 6. It is used as the */
/*        unit number for the printing of error messages. This printing */
/*        may be suppressed by setting LP to a non-positive value. */

/*    GTOL is a DOUBLE PRECISION variable with default value 0.9, which */
/*        controls the accuracy of the line search routine MCSRCH. If the */
/*        function and gradient evaluations are inexpensive with respect */
/*        to the cost of the iteration (which is sometimes the case when */
/*        solving very large problems) it may be advantageous to set GTOL */
/*        to a small value. A typical small value is 0.1.  Restriction: */
/*        GTOL should be greater than 1.D-04. */

/*    STPMIN and STPMAX are non-negative DOUBLE PRECISION variables which */
/*        specify lower and uper bounds for the step in the line search.  */
/*        Their default values are 1.D-20 and 1.D+20, respectively. These */
/*        values need not be modified unless the exponents are too large */
/*        for the machine being used, or unless the problem is extremely */
/*        badly scaled (in which case the exponents should be increased).  */


/*  MACHINE DEPENDENCIES */

/*        The only variables that are machine-dependent are XTOL, */
/*        STPMIN and STPMAX. */


/*  GENERAL INFORMATION */

/*    Other routines called directly:  DAXPY, DDOT, LB1, MCSRCH */

/*    Input/Output  :  No input; diagnostic messages on unit MP and */
/*                     error messages on unit LP. */


/*    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

/*      BLOCK DATA LB3 */

/*     INITIALIZE */
/*     ---------- */

    char buf[256];
    Real *diag, *gradient, *w;
    int info;
    bool converged = false;
    ispt = n + (m << 1);
    iypt = ispt + n * m;
    ftol = Real(1e-4);
    maxfev = 20;

    iter = 0;
    if (n <= 0 || m <= 0) {
       SimTK_THROW1(SimTK::Exception::OptimizerFailed , "IMPROPER INPUT PARAMETERS N OR M ARE NOT POSITIVE");
    }
    diag =     new Real[n];
    w =        new Real[n*(2*m+1) + 2*m];
    gradient = new Real[n];
    nfun = 1;
    point = 0;

    /* compute initial function and gradient values */
    objectiveFuncWrapper( n, x, true, f, this);
    gradientFuncWrapper( n,  x, false, gradient, this);

    // sherm 100303:
    // Scale the absolute gradient df/dx by f and x to get a relative
    // gradient (%chg f)/(%chg x). Then if 100% change of an x doesn't
    // produce at least a change eps*f we can say that x has converged.
    // This is the convergence criteria used for BFGS in Numerical Recipes
    // in C++ 2nd ed. p433 and makes more sense to me than the one that
    // came with this LBFGS implementation.
    // This is a scaled infinity-norm.

    Real fscale = 1 / std::max(Real(0.1), std::abs(*f));

    // Make a quick attempt to quit if we're already converged. No need
    // to calculate the whole norm here; we'll stop as soon as we find
    // a non-converged element.
    converged = true;
    for (int i=0; i<n; ++i) {
        const Real xscale = std::max(Real(1), std::abs(x[i]));
        if( std::abs(gradient[i])*fscale*xscale > *eps ) {
           converged=false;
           break;
        }
    }
    if( converged ) {
        delete [] diag;
        delete [] gradient;
        delete [] w;
        return;   // check if starting at minimum
    }

    // This is the unscaled 2-norm of the gradient.
    gnorm = std::sqrt(DDOT(n, gradient, c__1, gradient, c__1));

    stp1 = Real(1) / gnorm;

    for (int i = 0; i < n; ++i) {
        diag[i] = 1.;
    }

/*     THE WORK VECTOR W IS DIVIDED AS FOLLOWS: */
/*     --------------------------------------- */
/*     THE FIRST N LOCATIONS ARE USED TO STORE THE GRADIENT AND */
/*         OTHER TEMPORARY INFORMATION. */
/*     LOCATIONS (N+1)...(N+M) STORE THE SCALARS RHO. */
/*     LOCATIONS (N+M+1)...(N+2M) STORE THE NUMBERS ALPHA USED */
/*         IN THE FORMULA THAT COMPUTES H*G. */
/*     LOCATIONS (N+2M+1)...(N+2M+NM) STORE THE LAST M SEARCH */
/*         STEPS. */
/*     LOCATIONS (N+2M+NM+1)...(N+2M+2NM) STORE THE LAST M */
/*         GRADIENT DIFFERENCES. */

/*     THE SEARCH STEPS AND GRADIENT DIFFERENCES ARE STORED IN A */
/*     CIRCULAR ORDER CONTROLLED BY THE PARAMETER POINT. */

    for (int i = 0; i < n; ++i) {
        w[ispt + i] = -gradient[i] * diag[i];
    }

/*     PARAMETERS FOR LINE SEARCH ROUTINE */


    if (iprint[0] > 0) {
        lb1_(iprint, &iter, &nfun, &gnorm, &n, &m, x, f, gradient, &stp, &converged);
    }

/*    -------------------- */
/*     MAIN ITERATION LOOP */
/*    -------------------- */

    while( !converged ) {
       ++iter;
       info = 0;
       bound = iter - 1;
       if (iter != 1) {
          if (iter > m) {
              bound = m;
          }

          ys = DDOT(n, &w[iypt + npt], c__1, &w[ispt + npt], c__1);
          yy = DDOT(n, &w[iypt + npt], c__1, &w[iypt + npt], c__1);
          for (int i = 0; i < n; ++i) {
               diag[i] = ys / yy;
          }

/*     COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980, */
/*     "Updating quasi-Newton matrices with limited storage", */
/*     Mathematics of Computation, Vol.24, No.151, pp. 773-782. */
/*     --------------------------------------------------------- */

          cp = point;
          if (point == 0) {
           cp = m;
          }
          w[n + cp-1] = Real(1) / ys;
          for (int i = 0; i < n; ++i) {
              w[i] = -gradient[i];
          }
          cp = point;
          for (int i = 0; i < bound; ++i) {
              --cp;
              if (cp == -1) {
                  cp = m - 1;
              }
              sq = DDOT(n, &w[ispt + cp * n], c__1, w, c__1);
              inmc = n + m + cp;
              iycn = iypt + cp * n;
              w[inmc] = w[n + cp] * sq;
              d__1 = -w[inmc];
              DAXPY(n, d__1, &w[iycn], c__1, w, c__1);
          }

          for (int i = 0; i < n; ++i) {
              w[i] *= diag[i];
          }

          for (int i = 0; i < bound; ++i) {
              yr = DDOT(n, &w[iypt + cp * n], c__1, w, c__1);
              beta = w[n + cp] * yr;
              inmc = n + m + cp;
              beta = w[inmc] - beta;
              iscn = ispt + cp * n;
              DAXPY(n, beta, &w[iscn], c__1, w, c__1);
              ++cp;
              if (cp == m) {
                  cp = 0;
              }
          }

/*        STORE THE NEW SEARCH DIRECTION */
/*        ------------------------------ */

          for (int i = 0; i < n; ++i) {
             w[ispt + point * n + i] = w[i];
          }

/*        OBTAIN THE ONE-DIMENSIONAL MINIMIZER OF THE FUNCTION */
/*        BY USING THE LINE SEARCH ROUTINE MCSRCH */
/*        ---------------------------------------------------- */
      }
       nfev = 0;
/* awf changed initial step from ONE to be parametrized. */
       stp = lb3_1.stpawf;
       if (iter == 1) {
           stp = stp1;
       }
       for (int i = 0; i < n; ++i) {
           w[i] = gradient[i];
       }

       do {
          SimTK::LBFGSOptimizer::mcsrch_(&n, x, f, gradient,
                                         &w[ispt + point * n],
                                         &stp, &ftol, xtol,
                                         &maxfev, &info, &nfev, diag);
          if (info == -1) {
              objectiveFuncWrapper( n, x, true, f, this);
              gradientFuncWrapper( n,  x, false, gradient, this);
          } else if (info != 1) {
              delete [] diag;
              delete [] gradient;
              delete [] w;
              if (lb3_1.lp > 0) {
                 SimTK_THROW1(SimTK::Exception::OptimizerFailed ,
                 "LBFGS LINE SEARCH FAILED POSSIBLE CAUSES: FUNCTION OR GRADIENT ARE INCORRECT OR INCORRECT TOLERANCES");
              } else if( info == 0) {
                 SimTK_THROW1(SimTK::Exception::OptimizerFailed,
                 "LBFGS ERROR: IMPROPER INPUT PARAMETERS");
             } else if( info == 2) {
                 SimTK_THROW1(SimTK::Exception::OptimizerFailed,
                 "LBFGS ERROR: RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY IS AT MOST XTOL");
             } else if( info == 3) {
                 SimTK_THROW1(SimTK::Exception::OptimizerFailed,
                 "LBFGS ERROR: MORE THAN 20 FUNCTION EVALUATIONS WERE REQUIRED AT THE PRESENT ITERATION");
             } else if( info == 4) {
                 SimTK_THROW1(SimTK::Exception::OptimizerFailed,
                 "LBFGS ERROR: THE STEP IS TOO SMALL");
             } else if( info == 5) {
                 SimTK_THROW1(SimTK::Exception::OptimizerFailed,
                 "LBFGS ERROR: THE STEP IS TOO LARGE");
             } else if( info == 6){
                 SimTK_THROW1(SimTK::Exception::OptimizerFailed,
                 "LBFGS ERROR: ROUNDING ERRORS PREVENT FURTHER PROGRESS.\n THERE MAY NOT BE A STEP WHICH SATISFIES THE SUFFICIENT DECREASE AND CURVATURE\n CONDITIONS. TOLERANCES MAY BE TOO SMALL.");
             } else if( info == 7){
                 SimTK_THROW1(SimTK::Exception::OptimizerFailed , "Error in input parameters to MCSRCH");
             } else if( info == 8){
                 SimTK_THROW1(SimTK::Exception::OptimizerFailed , "THE SEARCH DIRECTION IS NOT A DESCENT DIRECTION");
             } else if( info == 9){
                 SimTK_THROW1(SimTK::Exception::OptimizerFailed , "Error in input parameters to mcstep_");
             } else {
                 sprintf(buf, "LBFGS ERROR: info = %d \n",info );
                 SimTK_THROW1(SimTK::Exception::OptimizerFailed, SimTK::String(buf) );
             }
          }
       } while( info == -1 );

       nfun += nfev;

/*     COMPUTE THE NEW STEP AND GRADIENT CHANGE */
/*     ----------------------------------------- */

       npt = point * n;
       for (int i = 0; i < n; ++i) {
           w[ispt + npt + i] *= stp;
           w[iypt + npt + i] = gradient[i] - w[i];
       }
       ++point;
       if (point == m) {
           point = 0;
       }

/*     TERMINATION TEST */
/*     ---------------- */

        //gnorm = std::sqrt(DDOT(n, gradient, c__1, gradient, c__1));
        //xnorm = std::sqrt(DDOT(n, x, c__1, x, c__1));
        //   xnorm = std::max(1.,xnorm);
        //   if (gnorm / xnorm <= *eps) {
        //       converged = true;
        //   }


        // sherm 100303: use scaled infinity norm instead
        fscale = 1 / std::max(Real(0.1), std::abs(*f));
        gnorm = 0;
        for (int i=0; i<n; ++i) {
            const Real xscale = std::max(Real(1), std::abs(x[i]));
            gnorm = std::max(gnorm, std::abs(gradient[i])*fscale*xscale);
        }
        converged = (gnorm <= *eps);

        if (iprint[0] > 0)
            lb1_(iprint, &iter, &nfun, &gnorm, &n, &m,
                 x, f, gradient, &stp, &converged);

    }  // end while loop
    delete [] diag;
    delete [] gradient;
    delete [] w;

/*     ------------------------------------------------------------ */
/*     END OF MAIN ITERATION LOOP.  */
/*     ------------------------------------------------------------ */

} /*  LAST LINE OF SUBROUTINE lbfgs_ */


/*      SUBROUTINE LB1(IPRINT,ITER,NFUN,GNORM,N,M,X,F,G,STP,FINISH) */
/*      ** moved to c file */

/*   ---------------------------------------------------------- */

/*   These routines removed for insertion into TargetJr netlib */

/* awf       subroutine daxpy(n,da,dx,incx,dy,incy) */
/* awf c */
/* awf c     constant times a vector plus a vector. */
/* awf c     uses unrolled loops for increments equal to one. */
/* awf c     jack dongarra, linpack, 3/11/78. */
/* awf c */
/* awf       Real precision dx(1),dy(1),da */
/* awf       integer i,incx,incy,ix,iy,m,mp1,n */
/* awf c */
/* awf       if(n.le.0)return */
/* awf       if (da .eq. 0.0d0) return */
/* awf       if(incx.eq.1.and.incy.eq.1)go to 20 */
/* awf c */
/* awf c        code for unequal increments or equal increments */
/* awf c          not equal to 1 */
/* awf c */
/* awf       ix = 1 */
/* awf       iy = 1 */
/* awf       if(incx.lt.0)ix = (-n+1)*incx + 1 */
/* awf       if(incy.lt.0)iy = (-n+1)*incy + 1 */
/* awf       do 10 i = 1,n */
/* awf         dy(iy) = dy(iy) + da*dx(ix) */
/* awf         ix = ix + incx */
/* awf         iy = iy + incy */
/* awf    10 continue */
/* awf       return */
/* awf c */
/* awf c        code for both increments equal to 1 */
/* awf c */
/* awf c */
/* awf c        clean-up loop */
/* awf c */
/* awf    20 m = mod(n,4) */
/* awf       if( m .eq. 0 ) go to 40 */
/* awf       do 30 i = 1,m */
/* awf         dy(i) = dy(i) + da*dx(i) */
/* awf    30 continue */
/* awf       if( n .lt. 4 ) return */
/* awf    40 mp1 = m + 1 */
/* awf       do 50 i = mp1,n,4 */
/* awf         dy(i) = dy(i) + da*dx(i) */
/* awf         dy(i + 1) = dy(i + 1) + da*dx(i + 1) */
/* awf         dy(i + 2) = dy(i + 2) + da*dx(i + 2) */
/* awf         dy(i + 3) = dy(i + 3) + da*dx(i + 3) */
/* awf    50 continue */
/* awf       return */
/* awf       end */
/* awf C */
/* awf C */
/* awf C   ---------------------------------------------------------- */
/* awf C */
/* awf       Real precision function ddot(n,dx,incx,dy,incy) */
/* awf c */
/* awf c     forms the dot product of two vectors. */
/* awf c     uses unrolled loops for increments equal to one. */
/* awf c     jack dongarra, linpack, 3/11/78. */
/* awf c */
/* awf       Real precision dx(1),dy(1),dtemp */
/* awf       integer i,incx,incy,ix,iy,m,mp1,n */
/* awf c */
/* awf       ddot = 0.0d0 */
/* awf       dtemp = 0.0d0 */
/* awf       if(n.le.0)return */
/* awf       if(incx.eq.1.and.incy.eq.1)go to 20 */
/* awf c */
/* awf c        code for unequal increments or equal increments */
/* awf c          not equal to 1 */
/* awf c */
/* awf       ix = 1 */
/* awf       iy = 1 */
/* awf       if(incx.lt.0)ix = (-n+1)*incx + 1 */
/* awf       if(incy.lt.0)iy = (-n+1)*incy + 1 */
/* awf       do 10 i = 1,n */
/* awf         dtemp = dtemp + dx(ix)*dy(iy) */
/* awf         ix = ix + incx */
/* awf         iy = iy + incy */
/* awf    10 continue */
/* awf       ddot = dtemp */
/* awf       return */
/* awf c */
/* awf c        code for both increments equal to 1 */
/* awf c */
/* awf c */
/* awf c        clean-up loop */
/* awf c */
/* awf    20 m = mod(n,5) */
/* awf       if( m .eq. 0 ) go to 40 */
/* awf       do 30 i = 1,m */
/* awf         dtemp = dtemp + dx(i)*dy(i) */
/* awf    30 continue */
/* awf       if( n .lt. 5 ) go to 60 */
/* awf    40 mp1 = m + 1 */
/* awf       do 50 i = mp1,n,5 */
/* awf         dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) + */
/* awf      *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4) */
/* awf    50 continue */
/* awf    60 ddot = dtemp */
/* awf       return */
/* awf       end */
/*    ------------------------------------------------------------------ */

/*     ************************** */
/*     LINE SEARCH ROUTINE MCSRCH */
/*     ************************** */

/* Subroutine */
void SimTK::LBFGSOptimizer::mcsrch_
   (integer *n, Real *x, Real *f, Real *g, Real *s, Real *stp,
    Real *ftol, Real *xtol, integer *maxfev,
    integer *info, integer *nfev, Real *wa)
{
    /* Initialized data */

    const Real xtrapf = 4.;

    /* Local variables */
    Real dgxm, dgym;
    integer j, infoc;
    Real finit, width, stmin, stmax;
    bool stage1;
    Real width1, ftest1, dg, fm, fx, fy;
    bool brackt;
    Real dginit, dgtest;
    Real dgm, dgx, dgy, fxm, fym, stx, sty;


/*                     SUBROUTINE MCSRCH */

/*     A slight modification of the subroutine CSRCH of More' and Thuente. */
/*     The changes are to allow reverse communication, and do not affect */
/*     the performance of the routine. */

/*     THE PURPOSE OF MCSRCH IS TO FIND A STEP WHICH SATISFIES */
/*     A SUFFICIENT DECREASE CONDITION AND A CURVATURE CONDITION. */

/*     AT EACH STAGE THE SUBROUTINE UPDATES AN INTERVAL OF */
/*     UNCERTAINTY WITH ENDPOINTS STX AND STY. THE INTERVAL OF */
/*     UNCERTAINTY IS INITIALLY CHOSEN SO THAT IT CONTAINS A */
/*     MINIMIZER OF THE MODIFIED FUNCTION */

/*          F(X+STP*S) - F(X) - FTOL*STP*(GRADF(X)'S). */

/*     IF A STEP IS OBTAINED FOR WHICH THE MODIFIED FUNCTION */
/*     HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE DERIVATIVE, */
/*     THEN THE INTERVAL OF UNCERTAINTY IS CHOSEN SO THAT IT */
/*     CONTAINS A MINIMIZER OF F(X+STP*S). */

/*     THE ALGORITHM IS DESIGNED TO FIND A STEP WHICH SATISFIES */
/*     THE SUFFICIENT DECREASE CONDITION */

/*           F(X+STP*S) .LE. F(X) + FTOL*STP*(GRADF(X)'S), */

/*     AND THE CURVATURE CONDITION */

/*           ABS(GRADF(X+STP*S)'S)) .LE. GTOL*ABS(GRADF(X)'S). */

/*     IF FTOL IS LESS THAN GTOL AND IF, FOR EXAMPLE, THE FUNCTION */
/*     IS BOUNDED BELOW, THEN THERE IS ALWAYS A STEP WHICH SATISFIES */
/*     BOTH CONDITIONS. IF NO STEP CAN BE FOUND WHICH SATISFIES BOTH */
/*     CONDITIONS, THEN THE ALGORITHM USUALLY STOPS WHEN ROUNDING */
/*     ERRORS PREVENT FURTHER PROGRESS. IN THIS CASE STP ONLY */
/*     SATISFIES THE SUFFICIENT DECREASE CONDITION. */

/*     THE SUBROUTINE STATEMENT IS */

/*        SUBROUTINE MCSRCH(N,X,F,G,S,STP,FTOL,XTOL, MAXFEV,INFO,NFEV,WA) */
/*     WHERE */

/*       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER */
/*         OF VARIABLES. */

/*       X IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE */
/*         BASE POINT FOR THE LINE SEARCH. ON OUTPUT IT CONTAINS */
/*         X + STP*S. */

/*       F IS A VARIABLE. ON INPUT IT MUST CONTAIN THE VALUE OF F */
/*         AT X. ON OUTPUT IT CONTAINS THE VALUE OF F AT X + STP*S. */

/*       G IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE */
/*         GRADIENT OF F AT X. ON OUTPUT IT CONTAINS THE GRADIENT */
/*         OF F AT X + STP*S. */

/*       S IS AN INPUT ARRAY OF LENGTH N WHICH SPECIFIES THE */
/*         SEARCH DIRECTION. */

/*       STP IS A NONNEGATIVE VARIABLE. ON INPUT STP CONTAINS AN */
/*         INITIAL ESTIMATE OF A SATISFACTORY STEP. ON OUTPUT */
/*         STP CONTAINS THE FINAL ESTIMATE. */

/*       FTOL AND GTOL ARE NONNEGATIVE INPUT VARIABLES. (In this reverse */
/*         communication implementation GTOL is defined in a COMMON */
/*         statement.) TERMINATION OCCURS WHEN THE SUFFICIENT DECREASE */
/*         CONDITION AND THE DIRECTIONAL DERIVATIVE CONDITION ARE */
/*         SATISFIED. */

/*       XTOL IS A NONNEGATIVE INPUT VARIABLE. TERMINATION OCCURS */
/*         WHEN THE RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY */
/*         IS AT MOST XTOL. */

/*       STPMIN AND STPMAX ARE NONNEGATIVE INPUT VARIABLES WHICH */
/*         SPECIFY LOWER AND UPPER BOUNDS FOR THE STEP. (In this reverse */
/*         communication implementatin they are defined in a COMMON */
/*         statement). */

/*       MAXFEV IS A POSITIVE INTEGER INPUT VARIABLE. TERMINATION */
/*         OCCURS WHEN THE NUMBER OF CALLS TO FCN IS AT LEAST */
/*         MAXFEV BY THE END OF AN ITERATION. */

/*       INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS: */

/*         INFO = 0  IMPROPER INPUT PARAMETERS. */

/*        INFO =-1  A RETURN IS MADE TO COMPUTE THE FUNCTION AND GRADIENT.  */
/*       NFEV IS AN INTEGER OUTPUT VARIABLE SET TO THE NUMBER OF */
/*         CALLS TO FCN. */

/*       WA IS A WORK ARRAY OF LENGTH N. */

/*     SUBPROGRAMS CALLED */

/*       MCSTEP */

/*       FORTRAN-SUPPLIED...ABS,MAX,MIN */

/*     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983 */
/*     JORGE J. MORE', DAVID J. THUENTE */

    infoc = 1;

/*     CHECK THE INPUT PARAMETERS FOR ERRORS. */

    if (*n <= 0 || *stp <= 0. || *ftol < 0. || lb3_1.gtol < 0. || *xtol < 0. ||
        lb3_1.stpmin < 0. || lb3_1.stpmax < lb3_1.stpmin || *maxfev <= 0) {
        *info = 7;
        return;
    }

/*     COMPUTE THE INITIAL GRADIENT IN THE SEARCH DIRECTION */
/*     AND CHECK THAT S IS A DESCENT DIRECTION. */

    dginit = 0.;
    for (j = 0; j < *n; ++j) {
        dginit += g[j] * s[j];
    }
    if (dginit >= 0.) {
        *info = 8;
        return;
    }

/*     INITIALIZE LOCAL VARIABLES. */

    brackt = FALSE_;
    stage1 = TRUE_;
    *nfev = 0;
    finit = *f;
    dgtest = *ftol * dginit;
    width = lb3_1.stpmax - lb3_1.stpmin;
    width1 = 2*width;
    for (j = 0; j < *n; ++j) {
        wa[j] = x[j];
    }

/*     THE VARIABLES STX, FX, DGX CONTAIN THE VALUES OF THE STEP, */
/*     FUNCTION, AND DIRECTIONAL DERIVATIVE AT THE BEST STEP. */
/*     THE VARIABLES STY, FY, DGY CONTAIN THE VALUE OF THE STEP, */
/*     FUNCTION, AND DERIVATIVE AT THE OTHER ENDPOINT OF */
/*     THE INTERVAL OF UNCERTAINTY. */
/*     THE VARIABLES STP, F, DG CONTAIN THE VALUES OF THE STEP, */
/*     FUNCTION, AND DERIVATIVE AT THE CURRENT STEP. */

    stx = 0.;
    fx = finit;
    dgx = dginit;
    sty = 0.;
    fy = finit;
    dgy = dginit;

/*     START OF ITERATION. */

L30:

/*        SET THE MINIMUM AND MAXIMUM STEPS TO CORRESPOND */
/*        TO THE PRESENT INTERVAL OF UNCERTAINTY. */

    if (brackt) {
        stmin = std::min(stx,sty);
        stmax = std::max(stx,sty);
    } else {
        stmin = stx;
        stmax = *stp + xtrapf * (*stp - stx);
    }

/*        FORCE THE STEP TO BE WITHIN THE BOUNDS STPMAX AND STPMIN. */

    *stp = std::max(*stp,lb3_1.stpmin);
    *stp = std::min(*stp,lb3_1.stpmax);

/*        IF AN UNUSUAL TERMINATION IS TO OCCUR THEN LET */
/*        STP BE THE LOWEST POINT OBTAINED SO FAR. */

    if ( (brackt && (*stp <= stmin || *stp >= stmax) ) || *nfev >= *maxfev - 1
        || infoc == 0 || ( brackt && stmax - stmin <= *xtol * stmax) ) {
        *stp = stx;
    }

/*        EVALUATE THE FUNCTION AND GRADIENT AT STP */
/*        AND COMPUTE THE DIRECTIONAL DERIVATIVE. */

    for (j = 0; j < *n; ++j) {
        x[j] = wa[j] + *stp * s[j];
    }

    objectiveFuncWrapper( *n, x, true, f, this);
    gradientFuncWrapper( *n,  x, false, g, this);

    *info = 0;
    ++(*nfev);
    dg = 0.;
    for (j = 0; j < *n; ++j) {
        dg += g[j] * s[j];
    }
    ftest1 = finit + *stp * dgtest;

/*        TEST FOR CONVERGENCE. */

    if ( (brackt && (*stp <= stmin || *stp >= stmax) ) || infoc == 0) {
        *info = 6;
    }
    if (*stp == lb3_1.stpmax && *f <= ftest1 && dg <= dgtest) {
        *info = 5;
    }
    if (*stp == lb3_1.stpmin && (*f > ftest1 || dg >= dgtest)) {
        *info = 4;
    }
    if (*nfev >= *maxfev) {
        *info = 3;
    }
    if (brackt && stmax - stmin <= *xtol * stmax) {
        *info = 2;
    }
    if (*f <= ftest1 && fabs(dg) <= lb3_1.gtol * (-dginit)) {
        *info = 1;
    }

/*        CHECK FOR TERMINATION. */

    if (*info != 0) {
        // Moved exception handling one level deeper Sept 2009 cmb
        if (*info != 1) {
            // Error condition - might want to breakpoint here
            return;
        }
        return;
    }

/*        IN THE FIRST STAGE WE SEEK A STEP FOR WHICH THE MODIFIED */
/*        FUNCTION HAS A NONPOSITIVE VALUE AND NONNEGATIVE DERIVATIVE. */

    if (stage1 && *f <= ftest1 && dg >= std::min(*ftol,lb3_1.gtol) * dginit) {
        stage1 = FALSE_;
    }

/*        A MODIFIED FUNCTION IS USED TO PREDICT THE STEP ONLY IF */
/*        WE HAVE NOT OBTAINED A STEP FOR WHICH THE MODIFIED */
/*        FUNCTION HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE */
/*        DERIVATIVE, AND IF A LOWER FUNCTION VALUE HAS BEEN */
/*        OBTAINED BUT THE DECREASE IS NOT SUFFICIENT. */

    if (stage1 && *f <= fx && *f > ftest1) {

/*           DEFINE THE MODIFIED FUNCTION AND DERIVATIVE VALUES. */

        fm = *f - *stp * dgtest;
        fxm = fx - stx * dgtest;
        fym = fy - sty * dgtest;
        dgm = dg - dgtest;
        dgxm = dgx - dgtest;
        dgym = dgy - dgtest;

/*           CALL CSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY */
/*           AND TO COMPUTE THE NEW STEP. */

        mcstep_(&stx, &fxm, &dgxm, &sty, &fym, &dgym, stp, &fm, &dgm, &brackt, &stmin, &stmax, &infoc);

/*           RESET THE FUNCTION AND GRADIENT VALUES FOR F. */

        fx = fxm + stx * dgtest;
        fy = fym + sty * dgtest;
        dgx = dgxm + dgtest;
        dgy = dgym + dgtest;
    } else {

/*           CALL MCSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY */
/*           AND TO COMPUTE THE NEW STEP. */

        mcstep_(&stx, &fx, &dgx, &sty, &fy, &dgy, stp, f, &dg, &brackt, &stmin, &stmax, &infoc);
    }

/*        FORCE A SUFFICIENT DECREASE IN THE SIZE OF THE */
/*        INTERVAL OF UNCERTAINTY. */

    if (brackt) {
        if (std::abs(sty - stx) >= .66 * width1) {
            *stp = stx + (sty - stx)/2;
        }
        width1 = width;
        width = std::abs(sty - stx);
    }

/*        END OF ITERATION. */

    goto L30;

/*     LAST LINE OF SUBROUTINE MCSRCH. */

} /* mcsrch_ */

/* Subroutine */
static void mcstep_(Real *stx, Real *fx, Real *dx, Real *sty, Real *fy, Real *dy,
                    Real *stp, Real *fp, Real *dp, bool *brackt,
                    Real *stpmin, Real *stpmax, integer *info)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    Real sgnd, stpc, stpf, stpq, p, q, gamma, r, s, theta;
    bool bound;


/*     SUBROUTINE MCSTEP */

/*     THE PURPOSE OF MCSTEP IS TO COMPUTE A SAFEGUARDED STEP FOR */
/*     A LINESEARCH AND TO UPDATE AN INTERVAL OF UNCERTAINTY FOR */
/*     A MINIMIZER OF THE FUNCTION. */

/*     THE PARAMETER STX CONTAINS THE STEP WITH THE LEAST FUNCTION */
/*     VALUE. THE PARAMETER STP CONTAINS THE CURRENT STEP. IT IS */
/*     ASSUMED THAT THE DERIVATIVE AT STX IS NEGATIVE IN THE */
/*     DIRECTION OF THE STEP. IF BRACKT IS SET TRUE THEN A */
/*     MINIMIZER HAS BEEN BRACKETED IN AN INTERVAL OF UNCERTAINTY */
/*     WITH ENDPOINTS STX AND STY. */

/*     THE SUBROUTINE STATEMENT IS */

/*       SUBROUTINE MCSTEP(STX,FX,DX,STY,FY,DY,STP,FP,DP,BRACKT, */
/*                        STPMIN,STPMAX,INFO) */

/*     WHERE */

/*       STX, FX, AND DX ARE VARIABLES WHICH SPECIFY THE STEP, */
/*         THE FUNCTION, AND THE DERIVATIVE AT THE BEST STEP OBTAINED */
/*         SO FAR. THE DERIVATIVE MUST BE NEGATIVE IN THE DIRECTION */
/*         OF THE STEP, THAT IS, DX AND STP-STX MUST HAVE OPPOSITE */
/*         SIGNS. ON OUTPUT THESE PARAMETERS ARE UPDATED APPROPRIATELY. */

/*       STY, FY, AND DY ARE VARIABLES WHICH SPECIFY THE STEP, */
/*         THE FUNCTION, AND THE DERIVATIVE AT THE OTHER ENDPOINT OF */
/*         THE INTERVAL OF UNCERTAINTY. ON OUTPUT THESE PARAMETERS ARE */
/*         UPDATED APPROPRIATELY. */

/*       STP, FP, AND DP ARE VARIABLES WHICH SPECIFY THE STEP, */
/*         THE FUNCTION, AND THE DERIVATIVE AT THE CURRENT STEP. */
/*         IF BRACKT IS SET TRUE THEN ON INPUT STP MUST BE */
/*         BETWEEN STX AND STY. ON OUTPUT STP IS SET TO THE NEW STEP. */

/*       BRACKT IS A LOGICAL VARIABLE WHICH SPECIFIES IF A MINIMIZER */
/*         HAS BEEN BRACKETED. IF THE MINIMIZER HAS NOT BEEN BRACKETED */
/*         THEN ON INPUT BRACKT MUST BE SET FALSE. IF THE MINIMIZER */
/*         IS BRACKETED THEN ON OUTPUT BRACKT IS SET TRUE. */

/*       STPMIN AND STPMAX ARE INPUT VARIABLES WHICH SPECIFY LOWER */
/*         AND UPPER BOUNDS FOR THE STEP. */

/*       INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS: */
/*         IF INFO = 1,2,3,4,5, THEN THE STEP HAS BEEN COMPUTED */
/*         ACCORDING TO ONE OF THE FIVE CASES BELOW. OTHERWISE */
/*         INFO = 0, AND THIS INDICATES IMPROPER INPUT PARAMETERS. */

/*     SUBPROGRAMS CALLED */

/*       FORTRAN-SUPPLIED ... ABS,MAX,MIN,SQRT */

/*     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983 */
/*     JORGE J. MORE', DAVID J. THUENTE */

    *info = 0;

/*     CHECK THE INPUT PARAMETERS FOR ERRORS. */

    if ( ( *brackt && ( *stp <= std::min(*stx,*sty) || *stp >= std::max(*stx,*sty) ) )
        || *dx * (*stp - *stx) >= 0. || *stpmax < *stpmin) {
        *info = 9;
        return;
    }

/*     DETERMINE IF THE DERIVATIVES HAVE OPPOSITE SIGN. */

    sgnd = *dp * (*dx / std::abs(*dx));

/*     FIRST CASE. A HIGHER FUNCTION VALUE. */
/*     THE MINIMUM IS BRACKETED. IF THE CUBIC STEP IS CLOSER */
/*     TO STX THAN THE QUADRATIC STEP, THE CUBIC STEP IS TAKEN, */
/*     ELSE THE AVERAGE OF THE CUBIC AND QUADRATIC STEPS IS TAKEN. */

    if (*fp > *fx) {
        *info = 1;
        bound = TRUE_;
        theta = (*fx - *fp) * 3 / (*stp - *stx) + *dx + *dp;
        s = std::max(std::max(std::abs(theta),std::abs(*dx)),std::abs(*dp));
        d__1 = theta / s;
        gamma = s * std::sqrt(d__1 * d__1 - *dx / s * (*dp / s));
        if (*stp < *stx) {
            gamma = -gamma;
        }
        p = gamma - *dx + theta;
        q = gamma - *dx + gamma + *dp;
        r = p / q;
        stpc = *stx + r * (*stp - *stx);
        stpq = *stx + *dx / ((*fx - *fp) / (*stp - *stx) + *dx) / 2 * (*stp - *stx);
        if (std::abs(stpc - *stx) < std::abs(stpq - *stx)) {
            stpf = stpc;
        } else {
            stpf = stpc + (stpq - stpc) / 2;
        }
        *brackt = TRUE_;

/*     SECOND CASE. A LOWER FUNCTION VALUE AND DERIVATIVES OF */
/*     OPPOSITE SIGN. THE MINIMUM IS BRACKETED. IF THE CUBIC */
/*     STEP IS CLOSER TO STX THAN THE QUADRATIC (SECANT) STEP, */
/*     THE CUBIC STEP IS TAKEN, ELSE THE QUADRATIC STEP IS TAKEN. */

    } else if (sgnd < 0.) {
        *info = 2;
        bound = FALSE_;
        theta = (*fx - *fp) * 3 / (*stp - *stx) + *dx + *dp;
        s = std::max(std::max(std::abs(theta),std::abs(*dx)),std::abs(*dp));
        d__1 = theta / s;
        gamma = s * std::sqrt(d__1 * d__1 - *dx / s * (*dp / s));
        if (*stp > *stx) {
            gamma = -gamma;
        }
        p = gamma - *dp + theta;
        q = gamma - *dp + gamma + *dx;
        r = p / q;
        stpc = *stp + r * (*stx - *stp);
        stpq = *stp + *dp / (*dp - *dx) * (*stx - *stp);
        if (std::abs(stpc - *stp) > std::abs(stpq - *stp)) {
            stpf = stpc;
        } else {
            stpf = stpq;
        }
        *brackt = TRUE_;

/*     THIRD CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE */
/*     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DECREASES. */
/*     THE CUBIC STEP IS ONLY USED IF THE CUBIC TENDS TO INFINITY */
/*     IN THE DIRECTION OF THE STEP OR IF THE MINIMUM OF THE CUBIC */
/*     IS BEYOND STP. OTHERWISE THE CUBIC STEP IS DEFINED TO BE */
/*     EITHER STPMIN OR STPMAX. THE QUADRATIC (SECANT) STEP IS ALSO */
/*     COMPUTED AND IF THE MINIMUM IS BRACKETED THEN THE THE STEP */
/*     CLOSEST TO STX IS TAKEN, ELSE THE STEP FARTHEST AWAY IS TAKEN. */

    } else if (std::abs(*dp) < std::abs(*dx)) {
        *info = 3;
        bound = TRUE_;
        theta = (*fx - *fp) * 3 / (*stp - *stx) + *dx + *dp;
        s = std::max(std::max(std::abs(theta),std::abs(*dx)),std::abs(*dp));

/*        THE CASE GAMMA = 0 ONLY ARISES IF THE CUBIC DOES NOT TEND */
/*        TO INFINITY IN THE DIRECTION OF THE STEP. */

        d__1 = theta / s;
        d__1 = d__1 * d__1 - *dx / s * (*dp / s);
        gamma = s * std::sqrt((std::max(Real(0),d__1)));
        if (*stp > *stx) {
            gamma = -gamma;
        }
        p = gamma - *dp + theta;
        q = gamma + (*dx - *dp) + gamma;
        r = p / q;
        if (r < 0. && gamma != 0.) {
            stpc = *stp + r * (*stx - *stp);
        } else if (*stp > *stx) {
            stpc = *stpmax;
        } else {
            stpc = *stpmin;
        }
        stpq = *stp + *dp / (*dp - *dx) * (*stx - *stp);
        if (*brackt) {
            if (std::abs(*stp - stpc) < std::abs(*stp - stpq)) {
                stpf = stpc;
            } else {
                stpf = stpq;
            }
        } else {
            if (std::abs(*stp - stpc) > std::abs(*stp - stpq)) {
                stpf = stpc;
            } else {
                stpf = stpq;
            }
        }

/*     FOURTH CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE */
/*     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DOES */
/*     NOT DECREASE. IF THE MINIMUM IS NOT BRACKETED, THE STEP */
/*     IS EITHER STPMIN OR STPMAX, ELSE THE CUBIC STEP IS TAKEN. */

    } else {
        *info = 4;
        bound = FALSE_;
        if (*brackt) {
            theta = (*fp - *fy) * 3 / (*sty - *stp) + *dy + *dp;
            s = std::max(std::max(std::abs(theta),std::abs(*dy)),std::abs(*dp));
            d__1 = theta / s;
            gamma = s * std::sqrt(d__1 * d__1 - *dy / s * (*dp / s));
            if (*stp > *sty) {
                gamma = -gamma;
            }
            p = gamma - *dp + theta;
            q = gamma - *dp + gamma + *dy;
            r = p / q;
            stpc = *stp + r * (*sty - *stp);
            stpf = stpc;
        } else if (*stp > *stx) {
            stpf = *stpmax;
        } else {
            stpf = *stpmin;
        }
    }

/*     UPDATE THE INTERVAL OF UNCERTAINTY. THIS UPDATE DOES NOT */
/*     DEPEND ON THE NEW STEP OR THE CASE ANALYSIS ABOVE. */

    if (*fp > *fx) {
        *sty = *stp;
        *fy = *fp;
        *dy = *dp;
    } else {
        if (sgnd < 0.) {
            *sty = *stx;
            *fy = *fx;
            *dy = *dx;
        }
        *stx = *stp;
        *fx = *fp;
        *dx = *dp;
    }

/*     COMPUTE THE NEW STEP AND SAFEGUARD IT. */

    stpf = std::min(*stpmax,stpf);
    stpf = std::max(*stpmin,stpf);
    *stp = stpf;
    if (*brackt && bound) {
        if (*sty > *stx) {
            d__1 = *stx + (*sty - *stx) * .66f;
            *stp = std::min(d__1,*stp);
        } else {
            d__1 = *stx + (*sty - *stx) * .66f;
            *stp = std::max(d__1,*stp);
        }
    }
} /* mcstep_ */



void lbptf_(const char* msg)
{
  printf(msg, 0); // dummy argument avoids gcc warning
}

void lbp1d_(const char* msg, int* i)
{
  printf(msg, *i);
}

void lbp1f_(const char* msg, Real* i)
{
  printf(msg, *i);
}

static void write50(Real* v, int n)
{
  int cols = 15;
  Real vmax = 0;
  int i;
  Real vmaxscale;
  for (i = 0; i < n; ++i)
    if (std::abs(v[i]) > vmax)
      vmax = v[i];
  vmaxscale = std::log(std::abs(vmax)) / std::log(Real(10));
  vmaxscale = std::pow(Real(10), ceil(vmaxscale) - 1);
  if (vmaxscale != 1.0)
    printf("  %e x\n", vmaxscale);

  for (i = 0; i < n; ++i) {
    if (i > 0 && i%cols == 0)
      printf("\n");
    printf(" %10.5f", v[i] / vmaxscale);
  }
  printf("\n");
}

/*C
//C     -------------------------------------------------------------
//C     THIS ROUTINE PRINTS MONITORING INFORMATION. THE FREQUENCY AND
//C     AMOUNT OF OUTPUT ARE CONTROLLED BY IPRINT.
//C     -------------------------------------------------------------
*/
void lb1_( int *iprint, int *iter, int *nfun, Real *gnorm, int *n,
           int *m, Real *x, Real *f, Real *g, Real *stp, bool *finish) /* bool*/
{
  (void)m;
  --iprint;
/* C*/
/*IF (ITER.EQ.0)THEN*/
  if (*iter == 0) {
/*  30   FORMAT(' F= ',1PD10.3,'   GNORM= ',1PD10.3)*/
/*       WRITE(MP,30)F,GNORM*/
    printf(" F = %g, GNORM = %g\n", *f, *gnorm);
/*       IF (IPRINT(2).GE.1)THEN*/
    if (iprint[2] >= 1) {
/*  40   FORMAT(' VECTOR X= ')*/
/*       WRITE(MP,40)*/
      printf(" VECTOR X=\n");
/*       WRITE(MP,50) (X(I),I=1,N)*/
      write50(x, *n);
/*  60   FORMAT(' GRADIENT VECTOR G= ')*/
/*       WRITE(MP,60)*/
      printf(" GRADIENT VECTOR G=\n");
/*       WRITE(MP,50) (G(I),I=1,N)*/
      write50(g, *n);
/*       ENDIF*/
    }
/*  10   FORMAT('*************************************************')*/
    printf("*************************************************\n");
/*  70   FORMAT(/'   I   NFN',4X,'FUNC',8X,'GNORM',7X,'STEPLENGTH'/)*/
/*       WRITE(MP,70)*/
    printf("   I   NFN    FUNC        GNORM       STEPLENGTH\n");
/*ELSE*/
  } else {
/*  IF ((IPRINT(1).EQ.0).AND.(ITER.NE.1.AND..NOT.FINISH))RETURN*/
    if ((iprint[1]==0) && (*iter != 1 && !*finish))
      return;
/*  IF (IPRINT(1).NE.0)THEN*/
    if (iprint[1] != 0) {
/*    IF(MOD(ITER-1,IPRINT(1)).EQ.0.OR.FINISH)THEN*/
      if ((*iter - 1)%iprint[1] == 0 || *finish) {
/*  70  FORMAT(/'   I   NFN',4X,'FUNC',8X,'GNORM',7X,'STEPLENGTH'/)*/
/*      IF(IPRINT(2).GT.1.AND.ITER.GT.1) WRITE(MP,70)*/
        if (iprint[2] > 1 && *iter > 1)
          printf("   I   NFN    FUNC        GNORM       STEPLENGTH\n");
/*  80  FORMAT(2(I4,1X),3X,3(1PD10.3,2X))*/
/*      WRITE(MP,80)ITER,NFUN,F,GNORM,STP*/
        printf("%4d %4d    %10.3f  %10.3f  %10.3f\n", *iter, *nfun, *f, *gnorm, *stp);
      }
/*    ELSE*/
      else {
/*      RETURN*/
        return;
/*    ENDIF*/
      }
    }
/*  ELSE*/
    else {

/*  70   FORMAT(/'   I   NFN',4X,'FUNC',8X,'GNORM',7X,'STEPLENGTH'/)*/
/*    IF( IPRINT(2).GT.1.AND.FINISH) WRITE(MP,70)*/
      if (iprint[2] > 1 && *finish)
        printf("   I   NFN    FUNC        GNORM       STEPLENGTH\n");

/*  80   FORMAT(2(I4,1X),3X,3(1PD10.3,2X))*/
/*    WRITE(MP,80)ITER,NFUN,F,GNORM,STP*/
      printf("%4d %4d    %10.3f  %10.3f  %10.3f\n", *iter, *nfun, *f, *gnorm, *stp);
/*  ENDIF*/
    }

/*  IF (IPRINT(2).EQ.2.OR.IPRINT(2).EQ.3)THEN*/
    if (iprint[2] == 2 || iprint[2] == 3) {
/*    IF (FINISH)THEN*/
      if (*finish)
/*  90  FORMAT(' FINAL POINT X= ')*/
/*      WRITE(MP,90)*/
        printf(" FINAL POINT X=\n");
/*    ELSE*/
      else
/*  40  FORMAT(' VECTOR X= ')*/
/*      WRITE(MP,40)*/
        printf(" VECTOR X=\n");
/*    ENDIF*/

/*  50   FORMAT(6(2X,1PD10.3))*/
/*     WRITE(MP,50)(X(I),I=1,N)*/
      write50(x, *n);
/*     IF (IPRINT(2).EQ.3)THEN*/
      if (iprint[2] == 3) {
/*  60   FORMAT(' GRADIENT VECTOR G= ')*/
/*       WRITE(MP,60)*/
        printf(" GRADIENT VECTOR G=\n");
/*  50   FORMAT(6(2X,1PD10.3))*/
/*       WRITE(MP,50)(G(I),I=1,N)*/
        write50(g, *n);
/*     ENDIF*/
      }
/*  ENDIF*/
    }
/*  100  FORMAT(/' THE MINIMIZATION TERMINATED WITHOUT DETECTING ERRORS.',*/
/* .       /' IFLAG = 0')*/
/*  IF (FINISH) WRITE(MP,100)*/
    if (*finish)
      printf(" THE MINIMIZATION TERMINATED WITHOUT DETECTING ERRORS.\n");
  }
/*  ENDIF*/
/* C*/
/*  RETURN*/
/*  END*/
}
