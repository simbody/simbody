/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2007/10/25 20:03:26 $
 * -----------------------------------------------------------------
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Van der Pol oscillator:
 *   xdotdot - 3*(1 - x^2)*xdot + x = 0, x(0) = 2, xdot(0) = 0.
 * This second-order ODE is converted to a first-order system by
 * defining y0 = x and y1 = xdot.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <cpodes/cpodes.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>


#define ODE CP_EXPL

/* Problem Constants */

#define ATOL RCONST(1.0e-6)   /* 1.0e-6 */
#define RTOL RCONST(0.0)      /* 0.0 */

#define ZERO   RCONST(0.0)
#define ONE    RCONST(1.0)
#define TWO    RCONST(2.0)
#define THIRTY RCONST(30.0)

#define P1_NEQ        2
#define P1_ETA        RCONST(3.0)
#define P1_NOUT       4
#define P1_T0         RCONST(0.0)
#define P1_T1         RCONST(1.39283880203)
#define P1_DTOUT      RCONST(2.214773875)

/* Private Helper Functions */
static void PrintFinalStats(void *cpode_mem);

/* Functions Called by the Solver */
static int res(realtype t, N_Vector y, N_Vector yp, N_Vector res, void *f_data);
static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data);

int main()
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
    fct = (void *)f;
  } else {
    fct = (void *)res;
    f(P1_T0, y, yp, NULL);
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

  return 0;
}

static int res(realtype t, N_Vector y, N_Vector yp, N_Vector res, void *f_data)
{
  f(t,y,res,f_data);
  N_VLinearSum(1.0, yp, -1.0, res, res);
  return(0);
}


static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
  realtype y0, y1;
  
  y0 = NV_Ith_S(y,0);
  y1 = NV_Ith_S(y,1);

  NV_Ith_S(ydot,0) = y1;
  NV_Ith_S(ydot,1) = (ONE - SQR(y0))* P1_ETA * y1 - y0;

  return(0);
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

