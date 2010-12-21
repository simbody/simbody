#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <cpodes/cpodes.h>
#include <cpodes/cpodes_band.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>

/* Problem Constants */

#define XMAX  RCONST(2.0)    /* domain boundaries         */
#define YMAX  RCONST(1.0)
#define MX    10             /* mesh dimensions           */
#define MY    5
#define NEQ   MX*MY          /* number of equations       */
#define ATOL  RCONST(1.0e-5) /* scalar absolute tolerance */
#define T0    RCONST(0.0)    /* initial time              */
#define T1    RCONST(0.1)    /* first output time         */
#define DTOUT RCONST(0.1)    /* output time increment     */
#define NOUT  10             /* number of output times    */

#define ZERO RCONST(0.0)
#define HALF RCONST(0.5)
#define ONE  RCONST(1.0)
#define TWO  RCONST(2.0)
#define FIVE RCONST(5.0)

#define IJth(vdata,i,j) (vdata[(j-1) + (i-1)*MY])

typedef struct {
  realtype dx, dy, hdcoef, hacoef, vdcoef;
} *UserData;

/* Private Helper Functions */
static void SetIC(N_Vector u, UserData data);
static void PrintHeader(realtype reltol, realtype abstol, realtype umax);
static void PrintOutput(realtype t, realtype umax, long int nst);
static void PrintFinalStats(void *cvode_mem);

/* Functions Called by the Solver */
static int f(realtype t, N_Vector u, N_Vector udot, void *f_data);
static int Jac(int N, int mu, int ml, 
               realtype t, N_Vector u, N_Vector fu, 
               DlsMat J, void *jac_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

int main(void)
{
  realtype dx, dy, reltol, abstol, t, tout, umax;
  N_Vector u, up;
  UserData data;
  void *cvode_mem;
  int iout, flag;
  long int nst;

  u = NULL;
  data = NULL;
  cvode_mem = NULL;

  u  = N_VNew_Serial(NEQ);
  up = N_VNew_Serial(NEQ);

  reltol = ZERO;
  abstol = ATOL;

  data = (UserData) malloc(sizeof *data);
  dx = data->dx = XMAX/(MX+1);
  dy = data->dy = YMAX/(MY+1);
  data->hdcoef = ONE/(dx*dx);
  data->hacoef = HALF/(TWO*dx);
  data->vdcoef = ONE/(dy*dy);

  SetIC(u, data);

  cvode_mem = CPodeCreate(CP_EXPL, CP_BDF, CP_NEWTON);
  flag = CPodeInit(cvode_mem, (void *)f, data, T0, u, NULL, CP_SS, reltol, &abstol);

  flag = CPBand(cvode_mem, NEQ, MY, MY);
  flag = CPDlsSetJacFn(cvode_mem, (void *)Jac, data);

  /* In loop over output points: call CPode, print results, test for errors */
  umax = N_VMaxNorm(u);
  PrintHeader(reltol, abstol, umax);
  for(iout=1, tout=T1; iout <= NOUT; iout++, tout += DTOUT) {
    flag = CPode(cvode_mem, tout, &t, u, up, CP_NORMAL);
    umax = N_VMaxNorm(u);
    flag = CPodeGetNumSteps(cvode_mem, &nst);
    PrintOutput(t, umax, nst);
  }

  PrintFinalStats(cvode_mem);

  N_VDestroy_Serial(u);
  CPodeFree(&cvode_mem);
  free(data);

  return(0);
}

/* f routine. Compute f(t,u). */
static int f(realtype t, N_Vector u,N_Vector udot, void *f_data)
{
  realtype uij, udn, uup, ult, urt, hordc, horac, verdc, hdiff, hadv, vdiff;
  realtype *udata, *dudata;
  int i, j;
  UserData data;

  udata = NV_DATA_S(u);
  dudata = NV_DATA_S(udot);

  /* Extract needed constants from data */

  data = (UserData) f_data;
  hordc = data->hdcoef;
  horac = data->hacoef;
  verdc = data->vdcoef;

  /* Loop over all grid points. */

  for (j=1; j <= MY; j++) {

    for (i=1; i <= MX; i++) {

      /* Extract u at x_i, y_j and four neighboring points */

      uij = IJth(udata, i, j);
      udn = (j == 1)  ? ZERO : IJth(udata, i, j-1);
      uup = (j == MY) ? ZERO : IJth(udata, i, j+1);
      ult = (i == 1)  ? ZERO : IJth(udata, i-1, j);
      urt = (i == MX) ? ZERO : IJth(udata, i+1, j);

      /* Set diffusion and advection terms and load into udot */

      hdiff = hordc*(ult - TWO*uij + urt);
      hadv = horac*(urt - ult);
      vdiff = verdc*(uup - TWO*uij + udn);
      IJth(dudata, i, j) = hdiff + hadv + vdiff;
    }
  }

  return(0);
}

/* Jacobian routine. Compute J(t,u). */
static int Jac(int N, int mu, int ml, 
               realtype t, N_Vector u, N_Vector fu, 
               DlsMat J, void *jac_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  long int i, j, k;
  realtype *kthCol, hordc, horac, verdc;
  UserData data;
  
  /*
    The components of f = udot that depend on u(i,j) are
    f(i,j), f(i-1,j), f(i+1,j), f(i,j-1), f(i,j+1), with
      df(i,j)/du(i,j) = -2 (1/dx^2 + 1/dy^2)
      df(i-1,j)/du(i,j) = 1/dx^2 + .25/dx  (if i > 1)
      df(i+1,j)/du(i,j) = 1/dx^2 - .25/dx  (if i < MX)
      df(i,j-1)/du(i,j) = 1/dy^2           (if j > 1)
      df(i,j+1)/du(i,j) = 1/dy^2           (if j < MY)
  */

  data = (UserData) jac_data;
  hordc = data->hdcoef;
  horac = data->hacoef;
  verdc = data->vdcoef;
  
  for (j=1; j <= MY; j++) {
    for (i=1; i <= MX; i++) {
      k = j-1 + (i-1)*MY;
      kthCol = BAND_COL(J,k);

      /* set the kth column of J */

      BAND_COL_ELEM(kthCol,k,k) = -TWO*(verdc+hordc);
      if (i != 1)  BAND_COL_ELEM(kthCol,k-MY,k) = hordc + horac;
      if (i != MX) BAND_COL_ELEM(kthCol,k+MY,k) = hordc - horac;
      if (j != 1)  BAND_COL_ELEM(kthCol,k-1,k)  = verdc;
      if (j != MY) BAND_COL_ELEM(kthCol,k+1,k)  = verdc;
    }
  }

  return(0);
}


/* Set initial conditions in u vector */
static void SetIC(N_Vector u, UserData data)
{
  int i, j;
  realtype x, y, dx, dy;
  realtype *udata;

  /* Extract needed constants from data */

  dx = data->dx;
  dy = data->dy;

  /* Set pointer to data array in vector u. */

  udata = NV_DATA_S(u);

  /* Load initial profile into u vector */
  
  for (j=1; j <= MY; j++) {
    y = j*dy;
    for (i=1; i <= MX; i++) {
      x = i*dx;
      IJth(udata,i,j) = x*(XMAX - x)*y*(YMAX - y)*EXP(FIVE*x*y);
    }
  }  
}

/* Print first lines of output (problem description) */
static void PrintHeader(realtype reltol, realtype abstol, realtype umax)
{
  printf("\n2-D Advection-Diffusion Equation\n");
  printf("Mesh dimensions = %d X %d\n", MX, MY);
  printf("Total system size = %d\n", NEQ);
  printf("Tolerance parameters: reltol = %lg   abstol = %lg\n\n", reltol, abstol);
  printf("At t = %lg      max.norm(u) =%14.6le \n", T0, umax);

  return;
}

/* Print current value */
static void PrintOutput(realtype t, realtype umax, long int nst)
{
  printf("At t = %4.2f   max.norm(u) =%14.6le   nst = %4ld\n", t, umax, nst);
  return;
}

/* Get and print some final statistics */

static void PrintFinalStats(void *cvode_mem)
{
  int flag;
  long int nst, nfe, nsetups, netf, nni, ncfn, nje, nfeLS;

  flag = CPodeGetNumSteps(cvode_mem, &nst);
  flag = CPodeGetNumFctEvals(cvode_mem, &nfe);
  flag = CPodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  flag = CPodeGetNumErrTestFails(cvode_mem, &netf);
  flag = CPodeGetNumNonlinSolvIters(cvode_mem, &nni);
  flag = CPodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);

  flag = CPDlsGetNumJacEvals(cvode_mem, &nje);
  flag = CPDlsGetNumFctEvals(cvode_mem, &nfeLS);

  printf("\nFinal Statistics:\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
	 nst, nfe, nsetups, nfeLS, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %ld\n \n",
	 nni, ncfn, netf);

  return;
}

