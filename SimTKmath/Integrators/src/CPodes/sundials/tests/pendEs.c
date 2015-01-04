#include <stdio.h>
#include <math.h>

#include <cpodes/cpodes.h>
#include <cpodes/cpodes_dense.h>
#include <nvector/nvector_serial.h>

/* Problem Constants */
#define g     RCONST(13.7503716373294544)
#define RTOL  RCONST(1.0e-10)
#define ATOL  RCONST(1.0e-10)

/* Functions Called by the Solver */
static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data);
static int cfun(realtype t, N_Vector yy, N_Vector cout, void *c_data);
static void PrintFinalStats(void *cpode_mem, booleantype proj);

/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

int main()
{
  void *cpode_mem;
  N_Vector yy, yp, ctols;
  realtype reltol, abstol, t, /*tout,*/ Tout;
  realtype half_pi, a, ad, E;
  int /*iout,*/ Nout, flag;

  half_pi = 2.0*atan(1.0);
  
  yy = N_VNew_Serial(2);
  yp = N_VNew_Serial(2);

  /* Initialize y */
  NV_Ith_S(yy,0) = half_pi;
  NV_Ith_S(yy,1) = 0.0;

  /* Set tolerances */
  reltol = RTOL;
  abstol = ATOL;

  Nout = 1;
  Tout = Nout*1.0;

  /* CREATE CPODES OBJECT */
  cpode_mem = CPodeCreate(CP_EXPL, CP_BDF, CP_NEWTON);  
  flag = CPodeSetMaxNumSteps(cpode_mem, 50000);
  flag = CPodeSetStopTime(cpode_mem, Tout);

  /* INITIALIZE SOLVER */
  flag = CPodeInit(cpode_mem, f, NULL, 0.0, yy, yp, CP_SS, reltol, &abstol);
  flag = CPDense(cpode_mem, 2);

  /* INTEGRATE IN ONE STEP MODE */
  t = 0.0;
  while(t<Tout) {
    flag = CPode(cpode_mem, Tout, &t, yy, yp, CP_ONE_STEP_TSTOP);
    if (flag < 0) break;
    a  = NV_Ith_S(yy,0);
    ad = NV_Ith_S(yy,1);
    E  = ad*ad - 2*g*cos(a);
    printf("%le  %14.10le  %14.10le     %14.10le\n",  t, a, ad, E);
  }

  /* INTEGRATE THROUGH A SEQUENCE OF TIMES */
  /*
  t = 0.0;
  for(iout=1; iout<=Nout; iout++) {
    tout = iout*1.0;
    flag = CPode(cpode_mem, tout, &t, yy, yp, CP_NORMAL_TSTOP);
    if (flag < 0) break;

    a  = NV_Ith_S(yy,0);
    ad = NV_Ith_S(yy,1);
    E  = ad*ad - 2*g*cos(a);
    printf("%le  %14.10le  %14.10le     %14.10le\n",  t, a, ad, E);
  }
  */

  PrintFinalStats(cpode_mem, FALSE);

  /* RE-INITIALIZE SOLVER */  
  NV_Ith_S(yy,0) = half_pi;
  NV_Ith_S(yy,1) = 0.0;

  printf("REINITIALIZE SOLVER\n\n");
  a  = NV_Ith_S(yy,0);
  ad = NV_Ith_S(yy,1);
  E  = ad*ad - 2*g*cos(a);
  printf("%14.10le  %14.10le     %14.10le\n\n",  a, ad, E);


  flag = CPodeReInit(cpode_mem, f, NULL, 0.0, yy, yp, CP_SS, reltol, &abstol);
  flag = CPDense(cpode_mem, 2);

  /* INTERNAL PROJECTION FUNCTION */  
  ctols = N_VNew_Serial(1);
  NV_Ith_S(ctols,0) = 1.0e-2;
  flag = CPodeProjInit(cpode_mem, CP_PROJ_L2NORM, CP_CNSTR_NONLIN, cfun, NULL, ctols);
  flag = CPodeSetProjNonlinConvCoef(cpode_mem, 1.0e-3);
  //  flag = CPodeSetProjTestCnstr(cpode_mem, TRUE);
  flag = CPDenseProj(cpode_mem, 1, 2, CPDIRECT_LU);

  /* INTEGRATE IN ONE STEP MODE */
  t = 0.0;
  while(t<Tout) {
    flag = CPode(cpode_mem, Tout, &t, yy, yp, CP_ONE_STEP_TSTOP);
    if (flag < 0) break;
    a  = NV_Ith_S(yy,0);
    ad = NV_Ith_S(yy,1);
    E  = ad*ad - 2*g*cos(a);
    printf("%le  %14.10le  %14.10le     %14.10le\n",  t, a, ad, E);
  }

 
  /* INTEGRATE THROUGH A SEQUENCE OF TIMES */
  /*
  t = 0.0;
  for(iout=1; iout<=Nout; iout++) {
    tout = iout*1.0;
    flag = CPode(cpode_mem, tout, &t, yy, yp, CP_NORMAL_TSTOP);
    if (flag < 0) break;

    a  = NV_Ith_S(yy,0);
    ad = NV_Ith_S(yy,1);
    E  = ad*ad - 2*g*cos(a);
    printf(" -------------- %le  %14.10le  %14.10le     %14.10le\n",  t, a, ad, E);
  }
  */

  PrintFinalStats(cpode_mem, TRUE);

  N_VDestroy_Serial(yy);
  N_VDestroy_Serial(yp);
  N_VDestroy_Serial(ctols);
  CPodeFree(&cpode_mem);

  return(0);
}


static int f(realtype t, N_Vector yy, N_Vector fy, void *f_data)
{
  realtype a, ad;

  a  = NV_Ith_S(yy,0);
  ad = NV_Ith_S(yy,1);

  NV_Ith_S(fy,0) = ad;
  NV_Ith_S(fy,1) = -g*sin(a);

  return(0);
}

static int cfun(realtype t, N_Vector yy, N_Vector cout, void *c_data)
{
  realtype a, ad;

  a  = NV_Ith_S(yy,0);
  ad = NV_Ith_S(yy,1);

  NV_Ith_S(cout,0) = ad*ad - 2*g*cos(a);

  return(0);
}

static void PrintFinalStats(void *cpode_mem, booleantype proj)
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

  printf("\nFinal Statistics:\n");
  printf("h0u = %g\n",h0u);
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld\n",
     nst, nfe, nsetups);
  printf("nfeLS = %-6ld nje = %ld\n",
     nfeLS, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld \n",
     nni, ncfn, netf);

  if (proj) {
    flag = CPodeGetProjStats(cpode_mem, &nproj, &nce, &nsetupsP, &nprf);
    printf("nproj = %-6ld nce = %-6ld nsetupsP = %-6ld nprf = %-6ld\n",
           nproj, nce, nsetupsP, nprf);
  }
}

