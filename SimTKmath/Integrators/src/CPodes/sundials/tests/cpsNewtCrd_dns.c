/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2007/10/25 20:03:25 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Example problem for the event detection in CPODES: a 2-pendulum
 * Newton's craddle. Two identical pendulums of length L and
 * mass M are suspended at points (+R,0) and (-R,0), where R is the
 * radius of the balls. 
 * We consider elastic impact with a coefficient of restitution C=1.
 *
 * Case 1: The first pendulum has an initial horizontal
 *         velocity V0, while the second one is at rest.
 * Case 2: Both pendulums are initially at rest. The first one is
 *         then moved from rest after a short interval, by applying
 *         an external force, F = m*frc(t), where
 *         frc(t) = 3(1-cos(5t)) on t =[0,2*pi/5] and 0 otherwise.
 *
 * Each pendulum is modeled with 2 coordinates (x and y) and one
 * constraint (x^2 + y^2 = L^2), resulting in a 1st order system
 * with 8 diferential equations and 4 constraints.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <cpodes/cpodes.h>
#include <cpodes/cpodes_dense.h>
#include <nvector/nvector_serial.h>

/* User data structure */

typedef struct {
  realtype R;
  realtype L;
  realtype m;
  realtype g;
} *PbData;

/* Functions Called by the Solver */

static int ffun1(realtype t, N_Vector y, N_Vector ydot, void *f_data);
static int ffun2(realtype t, N_Vector y, N_Vector ydot, void *f_data);

static int gfun(realtype t, N_Vector y, N_Vector yp, realtype *gout, void *g_data);

static int cfun(realtype t, N_Vector yy, N_Vector cout, void *c_data);

static void contact(N_Vector yy, PbData data);

static void PrintFinalStats(void *cpode_mem);

/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

int main()
{
  PbData data;
  realtype R, L, m, g, V0;
  void *cpode_mem;
  N_Vector yy, yp, ctols;
  realtype reltol, abstol;
  realtype t0, tf, t;
  realtype x1, y1, x2, y2;
  realtype vx1, vy1, vx2, vy2;
  int flag, Neq, Nc;
  int rdir[1], iroots[1];

  FILE *fout;

  /* --------------------------------
   * INITIALIZATIONS
   * -------------------------------- */

  R = 0.1;  /* ball radius */
  L = 1.0;  /* pendulum length */
  m = 1.0;  /* pendulum mass */
  g = 9.8;  /* gravitational acc. */

  V0 = 3.0; /* initial velocity for pendulum 1 */

  /* Set-up user data structure */

  data = (PbData)malloc(sizeof *data);
  data->R = R;
  data->L = L;
  data->m = m;
  data->g = g;

  /* Problem dimensions */

  Neq = 2*2*2;
  Nc  = 2*2;

  /* Solution vectors */

  yy = N_VNew_Serial(Neq);
  yp = N_VNew_Serial(Neq);

  /* Integration limits */

  t0 = 0.0;
  tf = 6.0;

  /* Integration and projection tolerances */

  reltol = 1.0e-8;
  abstol = 1.0e-8;

  ctols = N_VNew_Serial(Nc);
  N_VConst(1.0e-8, ctols);

  /* Direction of monitored events 
   * (only zero-crossing with decreasing even function) */

  rdir[0] = -1;

  /* --------------------------------
   * CASE 1
   * -------------------------------- */

  fout = fopen("newton1.out","w");

  /* Initial conditions */

  NV_Ith_S(yy,0) = R;       NV_Ith_S(yy,4) = -R;
  NV_Ith_S(yy,1) = -L;      NV_Ith_S(yy,5) = -L;
  NV_Ith_S(yy,2) = V0;      NV_Ith_S(yy,6) = 0.0;
  NV_Ith_S(yy,3) = 0.0;     NV_Ith_S(yy,7) = 0.0;

  /* Initialize solver */

  cpode_mem = CPodeCreate(CP_EXPL, CP_BDF, CP_NEWTON);  
  flag = CPodeInit(cpode_mem, ffun1, data, t0, yy, yp, CP_SS, reltol, &abstol);
  flag = CPDense(cpode_mem, Neq);

  flag = CPodeRootInit(cpode_mem, 1, gfun, data);
  flag = CPodeSetRootDirection(cpode_mem, rdir);

  /* Set-up the internal projection */
  
  flag = CPodeProjInit(cpode_mem, CP_PROJ_L2NORM, CP_CNSTR_NONLIN, cfun, data, ctols);
  flag = CPodeSetProjTestCnstr(cpode_mem, TRUE);
  flag = CPDenseProj(cpode_mem, Nc, Neq, CPDIRECT_LU);

  /* Integrate in ONE_STEP mode, while monitoring events */

  t = t0;
  while(t<tf) {

    flag = CPode(cpode_mem, tf, &t, yy, yp, CP_ONE_STEP);

    if (flag < 0) break;

    x1  = NV_Ith_S(yy,0);   x2  = NV_Ith_S(yy,4);
    y1  = NV_Ith_S(yy,1);   y2  = NV_Ith_S(yy,5);
    vx1 = NV_Ith_S(yy,2);   vx2 = NV_Ith_S(yy,6);
    vy1 = NV_Ith_S(yy,3);   vy2 = NV_Ith_S(yy,7);

    fprintf(fout, "%lf    %14.10lf  %14.10lf   %14.10lf  %14.10lf   %14.10lf  %14.10lf   %14.10lf  %14.10lf",
            t, x1, y1, x2, y2, vx1, vy1, vx2, vy2);

    if (flag == CP_ROOT_RETURN) {

      CPodeGetRootInfo(cpode_mem, iroots);
      fprintf(fout, " %d\n", iroots[0]);

      /* Note: the test iroots[0]<0 is really needed ONLY if not using rdir */

      if (iroots[0] < 0) {
        /* Update velocities in yy */
        contact(yy, data);
        /* reinitialize CPODES solver */
        flag = CPodeReInit(cpode_mem, ffun1, data, t, yy, yp, CP_SS, reltol, &abstol);
      }

    } else {

      fprintf(fout, " 0\n");

    }

  }

  PrintFinalStats(cpode_mem);

  CPodeFree(&cpode_mem);
    
  fclose(fout);

  /* --------------------------------
   * CASE 2
   * -------------------------------- */

  fout = fopen("newton2.out","w");

  /* Initial conditions */

  NV_Ith_S(yy,0) = R;       NV_Ith_S(yy,4) = -R;
  NV_Ith_S(yy,1) = -L;      NV_Ith_S(yy,5) = -L;
  NV_Ith_S(yy,2) = 0.0;     NV_Ith_S(yy,6) = 0.0;
  NV_Ith_S(yy,3) = 0.0;     NV_Ith_S(yy,7) = 0.0;

  /* Initialize solver */

  cpode_mem = CPodeCreate(CP_EXPL, CP_BDF, CP_NEWTON);  
  flag = CPodeInit(cpode_mem, ffun2, data, t0, yy, yp, CP_SS, reltol, &abstol);
  flag = CPDense(cpode_mem, Neq);

  flag = CPodeRootInit(cpode_mem, 1, gfun, data);
  flag = CPodeSetRootDirection(cpode_mem, rdir);

  /* Set-up the internal projection */
  
  flag = CPodeProjInit(cpode_mem, CP_PROJ_L2NORM, CP_CNSTR_NONLIN, cfun, data, ctols);
  flag = CPodeSetProjTestCnstr(cpode_mem, TRUE);
  flag = CPDenseProj(cpode_mem, Nc, Neq, CPDIRECT_LU);

  /* Integrate in ONE_STEP mode, while monitoring events */

  t = t0;
  while(t<tf) {

    flag = CPode(cpode_mem, tf, &t, yy, yp, CP_ONE_STEP);

    if (flag < 0) break;

    x1  = NV_Ith_S(yy,0);   x2 = NV_Ith_S(yy,4);
    y1  = NV_Ith_S(yy,1);   y2 = NV_Ith_S(yy,5);
    vx1 = NV_Ith_S(yy,2);   vx2 = NV_Ith_S(yy,6);
    vy1 = NV_Ith_S(yy,3);   vy2 = NV_Ith_S(yy,7);

    fprintf(fout, "%lf    %14.10lf  %14.10lf   %14.10lf  %14.10lf   %14.10lf  %14.10lf   %14.10lf  %14.10lf",
            t, x1, y1, x2, y2, vx1, vy1, vx2, vy2);

    if (flag == CP_ROOT_RETURN) {

      CPodeGetRootInfo(cpode_mem, iroots);
      fprintf(fout, " %d\n", iroots[0]);

      /* Note: the test iroots[0]<0 is really needed ONLY if not using rdir */

      if (iroots[0] < 0) {
        /* Update velocities in yy */
        contact(yy, data);
        /* reinitialize CPODES solver */
        flag = CPodeReInit(cpode_mem, ffun2, data, t, yy, yp, CP_SS, reltol, &abstol);
      }

    } else {

      fprintf(fout, " 0\n");

    }

  }

  PrintFinalStats(cpode_mem);

  CPodeFree(&cpode_mem);
    
  fclose(fout);

  /* --------------------------------
   * CLEAN-UP
   * -------------------------------- */

  free(data);
  N_VDestroy_Serial(yy);
  N_VDestroy_Serial(yp);
  N_VDestroy_Serial(ctols);

  return(0);
}


static int ffun1(realtype t, N_Vector yy, N_Vector fy, void *f_data)
{
  PbData data;
  realtype x1, y1, x2, y2;
  realtype vx1, vy1, vx2, vy2;
  realtype ax1, ay1, ax2, ay2;
  realtype lam1, lam2;
  realtype R, L, m, g;

  data = (PbData) f_data;
  R = data->R;
  L = data->L;
  m = data->m;
  g = data->g;

  x1  = NV_Ith_S(yy,0);   x2  = NV_Ith_S(yy,4);
  y1  = NV_Ith_S(yy,1);   y2  = NV_Ith_S(yy,5);
  vx1 = NV_Ith_S(yy,2);   vx2 = NV_Ith_S(yy,6);
  vy1 = NV_Ith_S(yy,3);   vy2 = NV_Ith_S(yy,7);
  
  lam1 = m/(2*L*L)*(vx1*vx1 + vy1*vy1 - g*y1);
  lam2 = m/(2*L*L)*(vx2*vx2 + vy2*vy2 - g*y2);
 
  ax1 = -2.0*(x1-R)*lam1/m;
  ay1 = -2.0*y1*lam1/m - g;

  ax2 = -2.0*(x2+R)*lam2/m;
  ay2 = -2.0*y2*lam2/m - g;
  
  NV_Ith_S(fy,0) = vx1;    NV_Ith_S(fy, 4) = vx2;
  NV_Ith_S(fy,1) = vy1;    NV_Ith_S(fy, 5) = vy2;
  NV_Ith_S(fy,2) = ax1;    NV_Ith_S(fy, 6) = ax2;
  NV_Ith_S(fy,3) = ay1;    NV_Ith_S(fy, 7) = ay2;

  return(0);
}

static int ffun2(realtype t, N_Vector yy, N_Vector fy, void *f_data)
{
  PbData data;
  realtype x1, y1, x2, y2;
  realtype vx1, vy1, vx2, vy2;
  realtype ax1, ay1, ax2, ay2;
  realtype lam1, lam2;
  realtype R, L, m, g;
  realtype pi, frc;


  pi = 4.0*atan(1.0);
  if ( (t <= 2*pi/5.0) ) {
    frc = 3.0 * ( 1.0 - cos(5.0*t) );
  } else {
    frc = 0.0;
  }

  data = (PbData) f_data;
  R = data->R;
  L = data->L;
  m = data->m;
  g = data->g;

  x1  = NV_Ith_S(yy,0);   x2  = NV_Ith_S(yy,4);
  y1  = NV_Ith_S(yy,1);   y2  = NV_Ith_S(yy,5);
  vx1 = NV_Ith_S(yy,2);   vx2 = NV_Ith_S(yy,6);
  vy1 = NV_Ith_S(yy,3);   vy2 = NV_Ith_S(yy,7);
  
  lam1 = m/(2*L*L)*(vx1*vx1 + vy1*vy1 + frc*(x1-R) - g*y1);
  lam2 = m/(2*L*L)*(vx2*vx2 + vy2*vy2 - g*y2);
 
  ax1 = -2.0*(x1-R)*lam1/m + frc;
  ay1 = -2.0*y1*lam1/m - g;

  ax2 = -2.0*(x2+R)*lam2/m;
  ay2 = -2.0*y2*lam2/m - g;
  
  NV_Ith_S(fy,0) = vx1;    NV_Ith_S(fy, 4) = vx2;
  NV_Ith_S(fy,1) = vy1;    NV_Ith_S(fy, 5) = vy2;
  NV_Ith_S(fy,2) = ax1;    NV_Ith_S(fy, 6) = ax2;
  NV_Ith_S(fy,3) = ay1;    NV_Ith_S(fy, 7) = ay2;

  return(0);
}

static int cfun(realtype t, N_Vector yy, N_Vector c, void *c_data)
{
  PbData data;
  realtype x1, y1, x2, y2;
  realtype vx1, vy1, vx2, vy2;
  realtype R, L;

  data = (PbData) c_data;
  R = data->R;
  L = data->L;

  x1  = NV_Ith_S(yy,0);   x2  = NV_Ith_S(yy,4);
  y1  = NV_Ith_S(yy,1);   y2  = NV_Ith_S(yy,5);
  vx1 = NV_Ith_S(yy,2);   vx2 = NV_Ith_S(yy,6);
  vy1 = NV_Ith_S(yy,3);   vy2 = NV_Ith_S(yy,7);

  NV_Ith_S(c,0) = (x1-R)*(x1-R) + y1*y1 - L*L;
  NV_Ith_S(c,1) = (x1-R)*vx1 + y1*vy1;

  NV_Ith_S(c,2) = (x2+R)*(x2+R) + y2*y2 - L*L;
  NV_Ith_S(c,3) = (x2+R)*vx2 + y2*vy2;

  return(0);
}

static int gfun(realtype t, N_Vector yy, N_Vector yp, realtype *gval, void *g_data)
{
  PbData data;
  realtype x1, y1, x2, y2;
  realtype R, L, m, g;

  data = (PbData) g_data;
  R = data->R;
  L = data->L;
  m = data->m;
  g = data->g;

  x1  = NV_Ith_S(yy,0);   x2  = NV_Ith_S(yy,4);
  y1  = NV_Ith_S(yy,1);   y2  = NV_Ith_S(yy,5);

  gval[0] = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) - 4.0*R*R;

  return(0);
}


static void contact(N_Vector yy, PbData data)
{
  realtype x1, y1, x2, y2;
  realtype vx1, vy1, vx2, vy2;
  realtype vt1, vn1, vn1_, vt2, vn2, vn2_;
  realtype alpha, ca, sa;

  x1  = NV_Ith_S(yy,0);   x2  = NV_Ith_S(yy,4);
  y1  = NV_Ith_S(yy,1);   y2  = NV_Ith_S(yy,5);
  vx1 = NV_Ith_S(yy,2);   vx2 = NV_Ith_S(yy,6);
  vy1 = NV_Ith_S(yy,3);   vy2 = NV_Ith_S(yy,7);

  /* Angle of contact line */

  alpha = atan2(y2-y1, x2-x1);
  ca = cos(alpha);
  sa = sin(alpha);

  /* Normal and tangential velocities before impact
   * (rotate velocity vectors by +alpha) */

  vn1 =  ca*vx1 + sa*vy1;
  vt1 = -sa*vx1 + ca*vy1;

  vn2 =  ca*vx2 + sa*vy2;
  vt2 = -sa*vx2 + ca*vy2;

  /* New normal velocities (M1=M2 and COR=1.0) */

  vn1_ = vn2;
  vn2_ = vn1;

  vn1 = vn1_;
  vn2 = vn2_;

  /* Velocities after impact (rotate back by -alpha) */

  vx1 = ca*vn1 - sa*vt1;
  vy1 = sa*vn1 + ca*vt1;

  vx2 = ca*vn2 - sa*vt2;
  vy2 = sa*vn2 + ca*vt2;
  
  NV_Ith_S(yy,2) = vx1;   NV_Ith_S(yy,6) = vx2;
  NV_Ith_S(yy,3) = vy1;   NV_Ith_S(yy,7) = vy2;

  return;
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

