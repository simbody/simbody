#include <stdio.h>

#include <cpodes/cpodes.h>
#include <cpodes/cpodes_dense.h>
#include <nvector/nvector_serial.h>

#define Ith(v,i)    NV_Ith_S(v,i-1)
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1)

/* Problem Constants */

#define ODE   CP_EXPL

#define NEQ   3
#define Y1    RCONST(1.0)
#define Y2    RCONST(0.0)
#define Y3    RCONST(0.0)

#define RTOL  RCONST(1.0e-4)   /* 1e-4 */
#define ATOL1 RCONST(1.0e-8)   /* 1e-8 */ 
#define ATOL2 RCONST(1.0e-14)  /* 1e-14 */
#define ATOL3 RCONST(1.0e-6)   /* 1e-6 */

#define T0    RCONST(0.0)
#define T1    RCONST(0.4)
#define TMULT RCONST(10.0)
#define NOUT  12

#define ONE   RCONST(1.0)
#define ZERO  RCONST(0.0)

/* Functions Called by the Solver */
static int res(realtype t, N_Vector y, N_Vector yp, N_Vector res, void *f_data);
static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data);
static int g(realtype t, N_Vector y, N_Vector yp, realtype *gout, void *g_data);
static int jacI(int N, realtype t, realtype gm,
                N_Vector y, N_Vector yp, N_Vector r, 
                DlsMat J, void *jac_data,
                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int jacE(int N, realtype t,
                N_Vector y, N_Vector yp,
                DlsMat J, void *jac_data,
                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static void PrintFinalStats(void *cpode_mem);

/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

int main()
{
  void *fct, *jac;
  void *cpode_mem;
  N_Vector yy, yp, abstol;
  realtype reltol, t, tout;
  int flag, flagr, iout;
  int rootsfound[2];

  /* Create serial vectors of length NEQ for I.C. and abstol */
  yy = N_VNew_Serial(NEQ);
  yp = N_VNew_Serial(NEQ);
  abstol = N_VNew_Serial(NEQ); 

  /* Initialize y */
  Ith(yy,1) = Y1;
  Ith(yy,2) = Y2;
  Ith(yy,3) = Y3;

  /* Set tolerances */
  reltol = RTOL;
  Ith(abstol,1) = ATOL1;
  Ith(abstol,2) = ATOL2;
  Ith(abstol,3) = ATOL3;


  if (ODE == CP_EXPL) {
    fct = (void *)f;
    jac = (void *)jacE;
  } else {
    f(T0, yy, yp, NULL);
    fct = (void *)res;
    jac = (void *)jacI;
  }

  /* Initialize solver */
  cpode_mem = CPodeCreate(ODE, CP_BDF, CP_NEWTON);  
  flag = CPodeInit(cpode_mem, fct, NULL, T0, yy, yp, CP_SV, reltol, abstol);

  /* Set initial step size */
  /*
  {
    realtype h0;
    h0 = 8.5e-14;
    flag = CPodeSetInitStep(cpode_mem, h0);
  }
  */

  /* Call CPodeRootInit to specify the root function g with 2 components */
  flag = CPodeRootInit(cpode_mem, 2, g, NULL);

  /* Call CPDense to specify the CPDENSE dense linear solver */
  flag = CPDense(cpode_mem, NEQ);

  /* Set the Jacobian routine (comment out for internal DQ) */
  flag = CPDlsSetJacFn(cpode_mem, jac, NULL);

  /* In loop, call CPode, print results, and test for error.
     Break out of loop when NOUT preset output times have been reached.  */
  iout = 0;  tout = T1;
  while(1) {
    flag = CPode(cpode_mem, tout, &t, yy, yp, CP_NORMAL);
    printf("At t = %0.4le      y =%14.6le  %14.6le  %14.6le\n", 
           t, Ith(yy,1), Ith(yy,2), Ith(yy,3));

    if (flag < 0) break;

    if (flag == CP_ROOT_RETURN) {
      flagr = CPodeGetRootInfo(cpode_mem, rootsfound);
      printf("    rootsfound[] = %3d %3d\n", rootsfound[0],rootsfound[1]);
    }

    if (flag == CP_SUCCESS) {
      iout++;
      tout *= TMULT;
    }

    if (iout == NOUT) break;
  }

  /* Print some final statistics */
  PrintFinalStats(cpode_mem);

  /* Free memory */
  N_VDestroy_Serial(yy);
  N_VDestroy_Serial(yp);
  N_VDestroy_Serial(abstol);
  CPodeFree(&cpode_mem);

  return(0);
}


/*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */

/*
 * residual function: F(y',y) = y' - f(y)
 */

static int res(realtype t, N_Vector y, N_Vector yp, N_Vector res, void *f_data)
{
  realtype y1, y2, y3, yp1, yp2, yp3;

  y1  = Ith(y,1);  y2  = Ith(y,2);  y3  = Ith(y,3);
  yp1 = Ith(yp,1); yp2 = Ith(yp,2); yp3 = Ith(yp,3);

  Ith(res,1) = yp1 - (RCONST(-0.04)*y1 + RCONST(1.0e4)*y2*y3);
  Ith(res,2) = yp2 + (RCONST(-0.04)*y1 + RCONST(1.0e4)*y2*y3 + RCONST(3.0e7)*y2*y2);
  Ith(res,3) = yp3 - RCONST(3.0e7)*y2*y2;

  return(0);
}

/*
 * f routine. Compute function f(t,y). 
 */

static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
  realtype y1, y2, y3, yd1, yd3;

  y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);

  yd1 = Ith(ydot,1) = RCONST(-0.04)*y1 + RCONST(1.0e4)*y2*y3;
  yd3 = Ith(ydot,3) = RCONST(3.0e7)*y2*y2;
        Ith(ydot,2) = -yd1 - yd3;

  return(0);
}

/*
 * Jacobian routine. Compute J(t,y) = dF/dy' + gm * dF/dy = I - gm*df/dy.
 */
static int jacI(int N, realtype t, realtype gm,
                N_Vector y, N_Vector yp, N_Vector r, 
                DlsMat J, void *jac_data,
                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{

  realtype y1, y2, y3;

  y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);

  IJth(J,1,1) = ONE - gm*RCONST(-0.04);
  IJth(J,1,2) = -gm*RCONST(1.0e4)*y3;
  IJth(J,1,3) = -gm*RCONST(1.0e4)*y2;

  IJth(J,2,1) = -gm*RCONST(0.04); 
  IJth(J,2,2) = ONE - gm * (RCONST(-1.0e4)*y3 - RCONST(6.0e7)*y2);
  IJth(J,2,3) = -gm*RCONST(-1.0e4)*y2;

  IJth(J,3,1) = ZERO;
  IJth(J,3,2) = -gm*RCONST(6.0e7)*y2;
  IJth(J,3,3) = ONE;

  return(0);
}

/*
 * Jacobian routine. Compute J(t,y) = df/dy. *
 */

static int jacE(int N, realtype t,
                N_Vector y, N_Vector fy, 
                DlsMat J, void *jac_data,
                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype y1, y2, y3;

  y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);

  IJth(J,1,1) = RCONST(-0.04);
  IJth(J,1,2) = RCONST(1.0e4)*y3;
  IJth(J,1,3) = RCONST(1.0e4)*y2;
  IJth(J,2,1) = RCONST(0.04); 
  IJth(J,2,2) = RCONST(-1.0e4)*y3-RCONST(6.0e7)*y2;
  IJth(J,2,3) = RCONST(-1.0e4)*y2;
  IJth(J,3,1) = ZERO;
  IJth(J,3,2) = RCONST(6.0e7)*y2;
  IJth(J,3,3) = ZERO;

  return(0);
}


/*
 * g routine. Compute functions g_i(t,y) for i = 0,1. 
 */

static int g(realtype t, N_Vector y, N_Vector yp, realtype *gout, void *g_data)
{
  realtype y1, y3;

  y1 = Ith(y,1); y3 = Ith(y,3);
  gout[0] = y1 - RCONST(0.0001);
  gout[1] = y3 - RCONST(0.01);

  return(0);
}

/* 
 * Get and print some final statistics
 */

static void PrintFinalStats(void *cpode_mem)
{
  realtype h0u;
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  int flag;

  flag = CPodeGetActualInitStep(cpode_mem, &h0u);
  flag = CPodeGetNumSteps(cpode_mem, &nst);
  flag = CPodeGetNumFctEvals(cpode_mem, &nfe);
  flag = CPodeGetNumLinSolvSetups(cpode_mem, &nsetups);
  flag = CPodeGetNumErrTestFails(cpode_mem, &netf);
  flag = CPodeGetNumNonlinSolvIters(cpode_mem, &nni);
  flag = CPodeGetNumNonlinSolvConvFails(cpode_mem, &ncfn);

  flag = CPDlsGetNumJacEvals(cpode_mem, &nje);
  flag = CPDlsGetNumFctEvals(cpode_mem, &nfeLS);

  flag = CPodeGetNumGEvals(cpode_mem, &nge);

  printf("\nFinal Statistics:\n");
  printf("h0u = %g\n",h0u);
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
     nst, nfe, nsetups, nfeLS, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n",
     nni, ncfn, netf, nge);
}

