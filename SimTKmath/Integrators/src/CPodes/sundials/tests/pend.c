#include <stdio.h>
#include <math.h>

#include <cpodes/cpodes.h>
#include <cpodes/cpodes_dense.h>
#include <nvector/nvector_serial.h>

#define Ith(v,i)    NV_Ith_S(v,i-1)
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1)

/* Problem Constants */

#define RTOL  RCONST(1.0e-8)
#define ATOL  RCONST(1.0e-8)

/* Functions Called by the Solver */
static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data);
static int proj(realtype t, N_Vector yy, N_Vector corr,
                realtype epsProj, N_Vector err, void *pdata);
static int cfun(realtype t, N_Vector yy, N_Vector cout, void *c_data);
static void PrintFinalStats(void *cpode_mem);

/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

int main()
{
  void *cpode_mem;
  N_Vector yy, yp, ctols;
  realtype reltol, abstol, t, tout, Tout;
  realtype x, y, xd, yd, g;
  int iout, Nout, flag;


  printf("double  %d\n",(int)sizeof(double));
  printf("long double  %d\n",(int)sizeof(long double));
  printf("float  %d\n",(int)sizeof(float));
  printf("realtype  %d\n",(int)sizeof(realtype));


  yy = N_VNew_Serial(4);
  yp = N_VNew_Serial(4);

  /* Initialize y */
  Ith(yy,1) = 1.0;  /* x */
  Ith(yy,2) = 0.0;  /* y */
  Ith(yy,3) = 0.0;  /* xd */
  Ith(yy,4) = 0.0;  /* yd */

  /* Set tolerances */
  reltol = RTOL;
  abstol = ATOL;

  Nout = 30;
  Tout = Nout*1.0;

  /* Initialize solver */
  cpode_mem = CPodeCreate(CP_EXPL, CP_BDF, CP_NEWTON);
  flag = CPodeInit(cpode_mem, f, NULL, 0.0, yy, yp, CP_SS, reltol, &abstol);
  flag = CPodeSetMaxNumSteps(cpode_mem, 5000);
  flag = CPodeSetStopTime(cpode_mem, Tout);
  flag = CPDense(cpode_mem, 4);

  /* USER-PROVIDED PROJECTION FUNCTION */
  /*
  flag = CPodeProjDefine(cpode_mem, proj, NULL);
  */

  /* INTERNAL PROJECTION FUNCTION */

  ctols = N_VNew_Serial(2);
  Ith(ctols,1) = 1.0e-8;
  Ith(ctols,2) = 1.0e-8;
  flag = CPodeProjInit(cpode_mem, CP_PROJ_L2NORM, CP_CNSTR_NONLIN, cfun, NULL, ctols);
  flag = CPodeSetProjTestCnstr(cpode_mem, TRUE);
  flag = CPDenseProj(cpode_mem, 2, 4, CPDIRECT_LU);

  flag = CPodeSetProjUpdateErrEst(cpode_mem, FALSE);

  /* DISABLE PROJECTION */
  /*
  CPodeSetProjFrequency(cpode_mem, 0);
  */

  /* INTEGRATE TO FINAL TIME */
  /*
  flag = CPode(cpode_mem, Tout, &t, yy, yp, CP_NORMAL_TSTOP);
  x  = Ith(yy,1);
  y  = Ith(yy,2);
  xd = Ith(yy,3);
  yd = Ith(yy,4);
  g = x*x + y*y - 1.0;
  Ith(yy,1) = x - 1.0;
  Ith(yy,2) = y - 0.0;
  Ith(yy,3) = xd - 0.0;
  Ith(yy,4) = yd - 0.0;
  printf("%14.10e  %14.10e  %14.10e  %14.10e  |  %14.10e\n",
         Ith(yy,1),Ith(yy,2),Ith(yy,3),Ith(yy,4),g);
  */

  /* INTEGRATE THROUGH A SEQUENCE OF TIMES */
  t = 0.0;
  for(iout=1; iout<=Nout; iout++) {
    tout = iout*1.0;
    flag = CPode(cpode_mem, tout, &t, yy, yp, CP_NORMAL_TSTOP);
    if (flag < 0) break;

    x  = Ith(yy,1);
    y  = Ith(yy,2);
    xd = Ith(yy,3);
    yd = Ith(yy,4);
    g = x*x + y*y - 1.0;
    printf(" -------------- %lf  %14.10lf  %14.10lf  %14.10lf  %14.10lf    %14.10lf\n",  t, x,y,xd,yd,g);
  }


  PrintFinalStats(cpode_mem);

  N_VDestroy_Serial(yy);
  N_VDestroy_Serial(yp);
  N_VDestroy_Serial(ctols);
  CPodeFree(&cpode_mem);

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

  x  = Ith(yy,1);
  y  = Ith(yy,2);
  xd = Ith(yy,3);
  yd = Ith(yy,4);

  R = sqrt(x*x+y*y);

  x_new = x/R;
  y_new = y/R;

  xd_new =   xd*y_new*y_new - yd*x_new*y_new;
  yd_new = - xd*x_new*y_new + yd*x_new*x_new;

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

  //  printf(":: %g  %g  %g  %g\n",e1,e2,e3,e4);
  //  printf(":: %g  %g  %g  %g\n",e1_new,e2_new,e3_new,e4_new);

  return(0);
}

static int cfun(realtype t, N_Vector yy, N_Vector cout, void *c_data)
{
  realtype x, y, xd, yd;

  x  = Ith(yy,1);
  y  = Ith(yy,2);
  xd = Ith(yy,3);
  yd = Ith(yy,4);

  Ith(cout,1) = x*x + y*y - 1.0;
  Ith(cout,2) = x*xd + y*yd;

  return(0);
}


static void PrintFinalStats(void *cpode_mem)
{
  realtype h0u;
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf;
  long int nproj, nce, nsetupsP, nprf;
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

  flag = CPodeGetProjStats(cpode_mem, &nproj, &nce, &nsetupsP, &nprf);

  printf("\nFinal Statistics:\n");
  printf("h0u = %g\n",h0u);
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld\n",
     nst, nfe, nsetups);
  printf("nfeLS = %-6ld nje = %ld\n",
     nfeLS, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld \n",
     nni, ncfn, netf);
  printf("nproj = %-6ld nce = %-6ld nsetupsP = %-6ld nprf = %-6ld\n",
         nproj, nce, nsetupsP, nprf);
}

