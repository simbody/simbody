#include <stdio.h>
#include <cpodes/cpodes.h>
#include <cpodes/cpodes_dense.h>
#include <nvector/nvector_serial.h>

#define Ith(v,i) NV_Ith_S(v,i-1)

#define RTOL  RCONST(1.0e-8)
#define ATOL  RCONST(1.0e-8)

static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data);

static int g(realtype t, N_Vector y, N_Vector yp, realtype *gout, void *g_data);

static void PrintFinalStats(void *cpode_mem);


/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

int main()
{
  realtype reltol, abstol, t, t0, tout, delt;
  N_Vector yy, yp;
  void *cpode_mem;
  int flag;
  int rdir[3], iroots[3];
  realtype gout[3];

  t0 = 0.0;
  tout = 10.0;
  delt = 0.01;

  yy = N_VNew_Serial(1);
  yp = N_VNew_Serial(1);
  Ith(yy,1) = 0.0;
  Ith(yp,1) = 0.0;

  reltol = RTOL;
  abstol = ATOL;

  cpode_mem = CPodeCreate(CP_EXPL, CP_BDF, CP_NEWTON);
  flag = CPodeInit(cpode_mem, (void *)f, NULL, t0, yy, NULL, CP_SS, reltol, &abstol);
  flag = CPDense(cpode_mem, 1);
  flag = CPodeSetStopTime(cpode_mem, tout);
  flag = CPodeRootInit(cpode_mem, 3, g, NULL);

  rdir[0] = -1;
  rdir[1] = 0;
  rdir[2] = 0;
  flag = CPodeSetRootDirection(cpode_mem, rdir);

  t = t0;
  g(t, yy, yp, gout, NULL);
  printf("%le %le   %le %le %le   0 0 0\n", t, Ith(yy,1), gout[0], gout[1], gout[2]);

  while(t<tout) {

    flag = CPode(cpode_mem, tout, &t, yy, yp, CP_ONE_STEP_TSTOP);
    if (flag < 0)  break;

    g(t, yy, yp, gout, NULL);
    printf("%le %le   %le %le %le   ", t, Ith(yy,1), gout[0], gout[1], gout[2]);
    if (flag == CP_ROOT_RETURN) {
      CPodeGetRootInfo(cpode_mem, iroots);
      printf("              %d %d %d\n", iroots[0], iroots[1], iroots[2]);
    } else {
      printf("0 0 0\n");
    }
      
  }

  //  PrintFinalStats(cpode_mem);

  /* Clean-up */

  N_VDestroy_Serial(yy);
  N_VDestroy_Serial(yp);
  CPodeFree(&cpode_mem);

  return(0);
}


/*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */

static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
  Ith(ydot,1) = t*cos(t) + sin(t);
  Ith(ydot,1) = cos(t);
  return(0);
}

static int g(realtype t, N_Vector y, N_Vector yp, realtype *gout, void *g_data)
{
  gout[0] = Ith(y,1);

  if (Ith(y,1)+0.41 >= 0) gout[1] = +0.6;
  else                   gout[1] = -0.2;

  if (Ith(y,1)-0.4 >= 0)      gout[2] = -0.4;
  else if (Ith(y,1)+0.4 <= 0) gout[2] = +0.4;
  else                        gout[2] = 0.0;

  return(0);
}


static void PrintFinalStats(void *cpode_mem)
{
  realtype h0u;
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
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

  flag = CPodeGetNumGEvals(cpode_mem, &nge);

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
  printf("nge = %ld\n", nge);
}

