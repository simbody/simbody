/*
 * -----------------------------------------------------------------
 * Problem 1: Van der Pol oscillator
 *   xdotdot - 3*(1 - x^2)*xdot + x = 0, x(0) = 2, xdot(0) = 0.
 * This second-order ODE is converted to a first-order system by
 * defining y0 = x and y1 = xdot.
 *
 * Problem 2: ydot = A * y, where A is a banded lower triangular
 * matrix derived from 2-D advection PDE.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <cpodes/cpodes.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>


#define ODE CP_EXPL

/* Shared Problem Constants */

#define ATOL RCONST(1.0e-6)   /* 1.0e-6 */
#define RTOL RCONST(0.0)      /* 0.0 */

#define ZERO   RCONST(0.0)
#define ONE    RCONST(1.0)
#define TWO    RCONST(2.0)
#define THIRTY RCONST(30.0)

/* Problem #1 Constants */

#define P1_NEQ        2
#define P1_ETA        RCONST(3.0)
#define P1_NOUT       4
#define P1_T0         RCONST(0.0)
#define P1_T1         RCONST(1.39283880203)
#define P1_DTOUT      RCONST(2.214773875)

/* Problem #2 Constants */

#define P2_MESHX      5
#define P2_MESHY      5
#define P2_NEQ        P2_MESHX*P2_MESHY
#define P2_ALPH1      RCONST(1.0)
#define P2_ALPH2      RCONST(1.0)
#define P2_NOUT       5
#define P2_ML         5
#define P2_MU         0
#define P2_T0         RCONST(0.0)
#define P2_T1         RCONST(0.01)
#define P2_TOUT_MULT  RCONST(10.0)

/* Private Helper Functions */
static void Problem1(void);
static void Problem2(void);
static realtype MaxError(N_Vector y, realtype t);
static void PrintFinalStats(void *cpode_mem);

/* Functions Called by the Solver */
static int res1(realtype t, N_Vector y, N_Vector yp, N_Vector res, void *f_data);
static int f1(realtype t, N_Vector y, N_Vector ydot, void *f_data);
static int res2(realtype t, N_Vector y, N_Vector yp, N_Vector res, void *f_data);
static int f2(realtype t, N_Vector y, N_Vector ydot, void *f_data);

int main(void)
{

  printf("Van der Pol\n\n");
  Problem1();
  printf("\n\ny' = A*y\n\n");
  Problem2();

  return(0);
}

static void Problem1(void)
{
  void *fct;
  void *cpode_mem;
  N_Vector y, yp;
  realtype reltol=RTOL, abstol=ATOL, t, tout, hu;
  int flag, iout, qu;

  y = NULL;
  yp = NULL;
  cpode_mem = NULL;

  y = N_VNew_Serial(P1_NEQ);
  NV_Ith_S(y,0) = TWO;
  NV_Ith_S(y,1) = ZERO;

  yp = N_VNew_Serial(P1_NEQ);
  
  if (ODE == CP_EXPL) {
    fct = (void *)f1;
  } else {
    fct = (void *)res1;
    f1(P1_T0, y, yp, NULL);
  }

  cpode_mem = CPodeCreate(ODE, CP_ADAMS, CP_FUNCTIONAL);
  /*  flag = CPodeSetInitStep(cpode_mem, 4.0e-9);*/
  flag = CPodeInit(cpode_mem, fct, NULL, P1_T0, y, yp, CP_SS, reltol, &abstol);

  printf("\n     t           x              xdot         qu     hu \n");
  for(iout=1, tout=P1_T1; iout <= P1_NOUT; iout++, tout += P1_DTOUT) {
    flag = CPode(cpode_mem, tout, &t, y, yp, CP_NORMAL);
    if (flag != CP_SUCCESS)  break;
    flag = CPodeGetLastOrder(cpode_mem, &qu);
    flag = CPodeGetLastStep(cpode_mem, &hu);    
    printf("%10.5f    %12.5le   %12.5le   %2d    %6.4le\n", t, NV_Ith_S(y,0), NV_Ith_S(y,1), qu, hu);
  }
  
  PrintFinalStats(cpode_mem);

  CPodeFree(&cpode_mem);
  N_VDestroy_Serial(y);
  N_VDestroy_Serial(yp);

  return;
}

static int res1(realtype t, N_Vector y, N_Vector yp, N_Vector res, void *f_data)
{
  f1(t,y,res,f_data);
  N_VLinearSum(1.0, yp, -1.0, res, res);
  return(0);
}


static int f1(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
  realtype y0, y1;
  
  y0 = NV_Ith_S(y,0);
  y1 = NV_Ith_S(y,1);

  NV_Ith_S(ydot,0) = y1;
  NV_Ith_S(ydot,1) = (ONE - SQR(y0))* P1_ETA * y1 - y0;

  return(0);
} 





static void Problem2(void)
{
  void *fct;
  void *cpode_mem;
  N_Vector y, yp;
  realtype reltol=RTOL, abstol=ATOL, t, tout, erm, hu;
  int flag, iout, qu;

  y = NULL;
  yp = NULL;
  cpode_mem = NULL;

  y = N_VNew_Serial(P2_NEQ);
  N_VConst(ZERO, y);
  NV_Ith_S(y,0) = ONE;

  yp = N_VNew_Serial(P2_NEQ);

  if (ODE == CP_EXPL) {
    fct = (void *)f2;
  } else {
    fct = (void *)res2;
    f2(P2_T0, y, yp, NULL);
  }

  cpode_mem = CPodeCreate(ODE, CP_ADAMS, CP_FUNCTIONAL);      
  /*  flag = CPodeSetInitStep(cpode_mem, 2.0e-9);*/
  flag = CPodeInit(cpode_mem, fct, NULL, P2_T0, y, yp, CP_SS, reltol, &abstol);
  
  printf("\n      t        max.err      qu     hu \n");
  for(iout=1, tout=P2_T1; iout <= P2_NOUT; iout++, tout*=P2_TOUT_MULT) {
    flag = CPode(cpode_mem, tout, &t, y, yp, CP_NORMAL);
    if (flag != CP_SUCCESS) break;
    erm = MaxError(y, t);
    flag = CPodeGetLastOrder(cpode_mem, &qu);
    flag = CPodeGetLastStep(cpode_mem, &hu);
    printf("%10.3f  %12.4le   %2d   %12.4le\n", t, erm, qu, hu);
  }
  
  PrintFinalStats(cpode_mem);
  
  CPodeFree(&cpode_mem);
  N_VDestroy_Serial(y);
  N_VDestroy_Serial(yp);

  return;
}

static int res2(realtype t, N_Vector y, N_Vector yp, N_Vector res, void *f_data)
{
  f2(t,y,res,f_data);
  N_VLinearSum(1.0, yp, -1.0, res, res);
  return(0);
}

static int f2(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
  long int i, j, k;
  realtype d, *ydata, *dydata;
  
  ydata = NV_DATA_S(y);
  dydata = NV_DATA_S(ydot);

  /*
     Excluding boundaries, 

     ydot    = f    = -2 y    + alpha1 * y      + alpha2 * y
         i,j    i,j       i,j             i-1,j             i,j-1
  */

  for (j=0; j < P2_MESHY; j++) {
    for (i=0; i < P2_MESHX; i++) {
      k = i + j * P2_MESHX;
      d = -TWO*ydata[k];
      if (i != 0) d += P2_ALPH1 * ydata[k-1];
      if (j != 0) d += P2_ALPH2 * ydata[k-P2_MESHX];
      dydata[k] = d;
    }
  }

  return(0);
}

static realtype MaxError(N_Vector y, realtype t)
{
  long int i, j, k;
  realtype *ydata, er, ex=ZERO, yt, maxError=ZERO, ifact_inv, jfact_inv=ONE;
  
  if (t == ZERO) return(ZERO);

  ydata = NV_DATA_S(y);
  if (t <= THIRTY) ex = EXP(-TWO*t); 
  
  for (j = 0; j < P2_MESHY; j++) {
    ifact_inv = ONE;
    for (i = 0; i < P2_MESHX; i++) {
      k = i + j * P2_MESHX;
      yt = RPowerI(t,i+j) * ex * ifact_inv * jfact_inv;
      er = ABS(ydata[k] - yt);
      if (er > maxError) maxError = er;
      ifact_inv /= (i+1);
    }
    jfact_inv /= (j+1);
  }
  return(maxError);
}



static void PrintFinalStats(void *cpode_mem)
{
  long int nst, nfe, nni, ncfn, netf;
  realtype h0u;
  int flag;
  
  flag = CPodeGetActualInitStep(cpode_mem, &h0u);
  flag = CPodeGetNumSteps(cpode_mem, &nst);
  flag = CPodeGetNumFctEvals(cpode_mem, &nfe);
  flag = CPodeGetNumErrTestFails(cpode_mem, &netf);
  flag = CPodeGetNumNonlinSolvIters(cpode_mem, &nni);
  flag = CPodeGetNumNonlinSolvConvFails(cpode_mem, &ncfn);

  printf("\n Final statistics:\n\n");
  printf(" Number of steps                          = %4ld \n", nst);
  printf(" Number of f-s                            = %4ld \n", nfe);
  printf(" Number of nonlinear iterations           = %4ld \n", nni);
  printf(" Number of nonlinear convergence failures = %4ld \n", ncfn);
  printf(" Number of error test failures            = %4ld \n", netf);
  printf(" Initial step size                        = %g \n\n", h0u);

}

