/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
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

/**@file
 * This is a duplicate of Radu's pendLr.c example converted to C++ and using
 * the SimTKcpodes interface instead of the original C one.
 */

//#define SimTK_USE_STATIC_LIBRARIES

#include "simmath/internal/SimTKcpodes.h"

// Just so we can get the version number:
#include "SimTKlapack.h"

#include <cstdio>
#include <iostream>

using SimTK::Real;
using SimTK::Vector;
using SimTK::Matrix;
using SimTK::CPodes;
using std::printf;
using std::cout;
using std::endl;

extern "C" {
    void SimTK_version_SimTKlapack(int* major, int* minor, int* build);
    void SimTK_about_SimTKlapack(const char* key, int maxlen, char* value);
}


static const Real RTOL = (Real)1e-8;
static const Real ATOL = (Real)1e-8;
static const Real CTOL = (Real)1e-8;

class PendLrSystem : public SimTK::CPodesSystem {
public:
    // Override default implementations of these virtual functions.
    int explicitODE(Real t, const Vector& y, Vector& fout) const override;
    int constraint(Real t, const Vector& y, Vector& cout) const override;
};

// Functions Called by the Solver
//static int f(Real t, N_Vector y, N_Vector ydot, void* f_data);
//static int cfun(Real t, N_Vector yy, N_Vector cout, void* c_data);

// Local static functions.
static void PrintFinalStats(CPodes&);

static void showSimTKAboutInfo() {

    int major,minor,build;
    char out[100];
    const char* keylist[] = { "version", "library", "type", "debug", "authors", "copyright", "svn_revision", 0 };

#ifdef TEST_LAPACK_VERSION
    SimTK_version_SimTKlapack(&major,&minor,&build);
    std::printf("==> SimTKlapack library version: %d.%d.%d\n", major, minor, build);
    std::printf("    SimTK_about_SimTKlapack():\n");
    for (const char** p = keylist; *p; ++p) {
        SimTK_about_SimTKlapack(*p, 100, out);
        std::printf("      about(%s)='%s'\n", *p, out);
    }
#endif

    SimTK_version_SimTKcommon(&major,&minor,&build);
    std::printf("==> SimTKcommon library version: %d.%d.%d\n", major, minor, build);
    std::printf("    SimTK_about_SimTKcommon():\n");
    for (const char** p = keylist; *p; ++p) {
        SimTK_about_SimTKcommon(*p, 100, out);
        std::printf("      about(%s)='%s'\n", *p, out);
    }

}

int main() {
  showSimTKAboutInfo();

try {
  Vector yy, yp, ctols;
  Real reltol, abstol, t, tout, Tout;
  Real x, y, xd, yd, g;
  int iout, Nout, flag;

  yy.resize(4);
  yp.resize(4);

  /* Initialize y */
  yy[0] = 1.0;  /* x */
  yy[1] = 0.0;  /* y */
  yy[2] = 0.0;  /* xd */
  yy[3] = 0.0;  /* yd */

  /* Set tolerances */
  reltol = RTOL;
  abstol = ATOL;

  Nout = 30;
  Tout = Nout*1.0;

  /* Initialize solver */
  PendLrSystem sys;
  CPodes cpode(CPodes::ExplicitODE, CPodes::BDF, CPodes::Newton);

  flag = cpode.init(sys, 0.0, yy, yp, CPodes::ScalarScalar, reltol, &abstol);
  flag = cpode.setMaxNumSteps(50000);
  flag = cpode.setStopTime(Tout);
  flag = cpode.lapackDense(4);

  /* INTERNAL PROJECTION FUNCTION */
  ctols.resize(/*3*/2);
  ctols[0] = CTOL;
  ctols[1] = CTOL;
  //ctols[2] = CTOL;
  flag = cpode.projInit(CPodes::L2Norm, CPodes::Nonlinear, ctols);
  flag = cpode.setProjTestCnstr(true);
  flag = cpode.lapackDenseProj(/*3*/2, 4, CPodes::ProjectWithQRPivot);

  /* INTEGRATE THROUGH A SEQUENCE OF TIMES */
  t = 0.0;
  for(iout=1; iout<=Nout; iout++) {
    tout = iout*1.0;
    flag = cpode.step(tout, &t, yy, yp, CPodes::NormalTstop); // CPode() in C
    if (flag < 0) break;

    x  = yy[0];
    y  = yy[1];
    xd = yy[2];
    yd = yy[3];
    g = (x*x + y*y - 1.0)/2;
    printf(" -------------- %lf  %14.10lf  %14.10lf  %14.10lf  %14.10lf    %14.10lf\n",  t, x,y,xd,yd,g);
  }

  PrintFinalStats(cpode);

  return 0;
}
catch (const std::exception& e) {
  std::cout << e.what() << std::endl;
}
  return 1;
}

static int pendODE(Real t, const Vector& yy, Vector& fy)
{
  Real x, y, xd, yd, g, tmp;

  g = 13.7503716373294544;

  x  = yy[0];
  y  = yy[1];
  xd = yy[2];
  yd = yy[3];



  // Radu's version:
  //tmp = xd*xd + yd*yd - g*y;
  tmp = (xd*xd + yd*yd - g*y)/(x*x + y*y);

  fy[0] = xd;
  fy[1] = yd;
  fy[2] = - x*tmp;
  fy[3] = - y*tmp - g;


  return(0);
}

int PendLrSystem::explicitODE(Real t, const Vector& yy,
                              Vector& fy) const
{
    return pendODE(t, yy, fy);
}

int PendLrSystem::constraint(Real t, const Vector& yy,
                             Vector& cout) const
{
  Real x, y, xd, yd;

  x  = yy[0];
  y  = yy[1];
  xd = yy[2];
  yd = yy[3];

  cout[0] = (x*x + y*y - 1.0)/2;
  cout[1] = x*xd + y*yd;
  //cout[2] = cout[0] + cout[1];

  return(0);
}


static void PrintFinalStats(CPodes& cpode)
{
  Real h0u;
  int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  int nproj, nce, nsetupsP, nprf;
  int flag;

  flag = cpode.getActualInitStep(&h0u);
  flag = cpode.getNumSteps(&nst);
  flag = cpode.getNumFctEvals(&nfe);
  flag = cpode.getNumLinSolvSetups(&nsetups);
  flag = cpode.getNumErrTestFails(&netf);
  flag = cpode.getNumNonlinSolvIters(&nni);
  flag = cpode.getNumNonlinSolvConvFails(&ncfn);

  flag = cpode.dlsGetNumJacEvals(&nje);
  flag = cpode.dlsGetNumFctEvals(&nfeLS);

  flag = cpode.getProjStats(&nproj, &nce, &nsetupsP, &nprf);

  flag = cpode.getNumGEvals(&nge);

  printf("\nFinal Statistics:\n");
  printf("h0u = %g\n",h0u);
  printf("nst = %-6d nfe  = %-6d nsetups = %-6d\n",
     nst, nfe, nsetups);
  printf("nfeLS = %-6d nje = %d\n",
     nfeLS, nje);
  printf("nni = %-6d ncfn = %-6d netf = %-6d \n",
     nni, ncfn, netf);
  printf("nproj = %-6d nce = %-6d nsetupsP = %-6d nprf = %-6d\n",
         nproj, nce, nsetupsP, nprf);
  printf("nge = %d\n", nge);

}

