#include <stdio.h>
#include <math.h>

#include <cpodes/cpodes.h>
#include <cpodes/cpodes_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>

#define Ith(v,i)    NV_Ith_S(v,i-1)

#define TOL RCONST(1.0e-5)
#define TOL_REF RCONST(1.0e-14)

static int fref(realtype t, N_Vector yy, N_Vector fy, void *f_data);

static int f(realtype t, N_Vector yy, N_Vector fy, void *f_data);
static int proj(realtype t, N_Vector yy, N_Vector corr, 
                realtype epsProj, N_Vector err, void *pdata);

void GetSol(void *cpode_mem, N_Vector yy0, realtype tol, 
            realtype tout, booleantype proj, N_Vector yref);

void RefSol(realtype tout, N_Vector yref);

/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

int main()
{
  void *cpode_mem;
  N_Vector yref, yy0;
  realtype tol, tout;
  int i, flag;


  tout = 30.0;
 
  /* Get reference solution */
  yref = N_VNew_Serial(4);
  RefSol(tout, yref);

  /* Initialize solver */
  tol = TOL;
  cpode_mem = CPodeCreate(CP_EXPL, CP_BDF, CP_NEWTON);  
  yy0 = N_VNew_Serial(4);
  Ith(yy0,1) = 1.0;  /* x */
  Ith(yy0,2) = 0.0;  /* y */
  Ith(yy0,3) = 0.0;  /* xd */
  Ith(yy0,4) = 0.0;  /* yd */
  flag = CPodeInit(cpode_mem, (void *)f, NULL, 0.0, yy0, NULL, CP_SS, tol, &tol);
  flag = CPodeSetMaxNumSteps(cpode_mem, 50000);
  flag = CPodeSetStopTime(cpode_mem, tout);
  flag = CPodeProjDefine(cpode_mem, proj, NULL);
  flag = CPDense(cpode_mem, 4);

  for (i=0;i<5;i++) {

    printf("\n\n%.2e\n", tol);
    GetSol(cpode_mem, yy0, tol, tout, TRUE, yref);
    GetSol(cpode_mem, yy0, tol, tout, FALSE, yref);
    tol /= 10.0;
  }

  N_VDestroy_Serial(yref);
  CPodeFree(&cpode_mem);

  return(0);
}

void GetSol(void *cpode_mem, N_Vector yy0, realtype tol, 
            realtype tout, booleantype proj, N_Vector yref)
{
  N_Vector yy, yp;
  realtype t, x, y, xd, yd, g;
  int flag;
  long int nst, nfe, nsetups, nje, nfeLS, ncfn, netf;

  if (proj) {
    printf(" YES   ");
    CPodeSetProjFrequency(cpode_mem, 1);
  } else {
    CPodeSetProjFrequency(cpode_mem, 0);
    printf(" NO    ");
  }

  yy = N_VNew_Serial(4);
  yp = N_VNew_Serial(4);

  flag = CPodeReInit(cpode_mem, (void *)f, NULL, 0.0, yy0, NULL, CP_SS, tol, &tol);

  flag = CPode(cpode_mem, tout, &t, yy, yp, CP_NORMAL_TSTOP);

  x  = Ith(yy,1);
  y  = Ith(yy,2);
  g = ABS(x*x + y*y - 1.0);

  N_VLinearSum(1.0, yy, -1.0, yref, yy);

  N_VAbs(yy, yy);

  x  = Ith(yy,1);
  y  = Ith(yy,2);  
  xd = Ith(yy,3);
  yd = Ith(yy,4);


  printf("%9.2e  %9.2e  %9.2e  %9.2e  |  %9.2e  |",  
         Ith(yy,1),Ith(yy,2),Ith(yy,3),Ith(yy,4),g);

  CPodeGetNumSteps(cpode_mem, &nst);
  CPodeGetNumFctEvals(cpode_mem, &nfe);
  CPodeGetNumLinSolvSetups(cpode_mem, &nsetups);
  CPodeGetNumErrTestFails(cpode_mem, &netf);
  CPodeGetNumNonlinSolvConvFails(cpode_mem, &ncfn);

  CPDlsGetNumJacEvals(cpode_mem, &nje);
  CPDlsGetNumFctEvals(cpode_mem, &nfeLS);


  printf(" %6ld   %6ld+%-4ld  %4ld (%3ld)  |  %3ld  %3ld\n",
         nst, nfe, nfeLS, nsetups, nje, ncfn, netf);

  N_VDestroy_Serial(yy);
  N_VDestroy_Serial(yp);

  return;
}


void RefSol(realtype tout, N_Vector yref)
{
  void *cpode_mem;
  N_Vector yy, yp;
  realtype tol, t, th, thd;
  int flag;
  

  yy = N_VNew_Serial(2);
  yp = N_VNew_Serial(2);
  Ith(yy,1) = 0.0;  /* theta */
  Ith(yy,2) = 0.0;  /* thetad */
  tol = TOL_REF;

  cpode_mem = CPodeCreate(CP_EXPL, CP_BDF, CP_NEWTON);  
  flag = CPodeSetMaxNumSteps(cpode_mem, 100000);
  flag = CPodeInit(cpode_mem, (void *)fref, NULL, 0.0, yy, yp, CP_SS, tol, &tol);
  flag = CPDense(cpode_mem, 2);

  flag = CPodeSetStopTime(cpode_mem, tout);
  flag = CPode(cpode_mem, tout, &t, yy, yp, CP_NORMAL_TSTOP);
  th  = Ith(yy,1);
  thd = Ith(yy,2);
  Ith(yref,1) = cos(th);
  Ith(yref,2) = sin(th);
  Ith(yref,3) = -thd*sin(th);
  Ith(yref,4) =  thd*cos(th);
  
  N_VDestroy_Serial(yy);
  N_VDestroy_Serial(yp);
  CPodeFree(&cpode_mem);

  return;
}

static int fref(realtype t, N_Vector yy, N_Vector fy, void *f_data)
{
  realtype th, thd, g;

  g = 13.7503716373294544;

  th  = Ith(yy,1);
  thd  = Ith(yy,2);

  Ith(fy,1) = thd;
  Ith(fy,2) = -g*cos(th);

  return(0);
}



static int f(realtype t, N_Vector yy, N_Vector fy, void *f_data)
{
  realtype x, y, xd, yd, g, tmp;

  g = 13.7503716373294544;

  x  = Ith(yy,1);
  y  = Ith(yy,2);
  xd = Ith(yy,3);
  yd = Ith(yy,4);
 
  tmp = xd*xd + yd*yd - g*y;

  Ith(fy,1) = xd;
  Ith(fy,2) = yd;
  Ith(fy,3) = -x*tmp;
  Ith(fy,4) = -y*tmp - g;

  return(0);
}

static int proj(realtype t, N_Vector yy, N_Vector corr,
                realtype epsProj, N_Vector err, void *pdata)
{
  realtype x, y, xd, yd;
  realtype x_new, y_new, xd_new, yd_new;
  realtype e1, e2, e3, e4;
  realtype e1_new, e2_new, e3_new, e4_new;
  realtype R;

  /* Extract current solution */

  x  = Ith(yy,1);
  y  = Ith(yy,2);
  xd = Ith(yy,3);
  yd = Ith(yy,4);

  /* Project onto manifold */

  R = sqrt(x*x+y*y);
  
  x_new = x/R;
  y_new = y/R;

  xd_new =   xd*y_new*y_new - yd*x_new*y_new;
  yd_new = - xd*x_new*y_new + yd*x_new*x_new;

  /* Return corrections */

  Ith(corr,1) = x_new  - x;
  Ith(corr,2) = y_new  - y;
  Ith(corr,3) = xd_new - xd;
  Ith(corr,4) = yd_new - yd;

  /*      +-            -+
   *      |  y*y    -x*y |
   *  P = |              |
   *      | -x*y     x*x |
   *      +-            -+
   */
  
  /* Return err <-  P * err */

  e1 = Ith(err,1);
  e2 = Ith(err,2);
  e3 = Ith(err,3);
  e4 = Ith(err,4);

  e1_new =  y_new*y_new * e1 - x_new*y_new * e2;
  e2_new = -x_new*y_new * e1 + x_new*x_new * e2;

  e3_new =  y_new*y_new * e3 - x_new*y_new * e4;
  e4_new = -x_new*y_new * e3 + x_new*x_new * e4;

  Ith(err,1) = e1_new;
  Ith(err,2) = e2_new;
  Ith(err,3) = e3_new;
  Ith(err,4) = e4_new;

  return(0);
}
