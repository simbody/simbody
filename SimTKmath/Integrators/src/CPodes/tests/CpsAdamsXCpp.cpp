/* -------------------------------------------------------------------------- *
 *                         SimTK Simbody: SimTKmath                           *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-12 Stanford University and the Authors.        *
 * Authors: Jack Middleton                                                    *
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
 * This is a duplicate of Radu's cpsadamsx.c example converted to C++ and using
 * the SimTKcpodes interface instead of the original C one.
 */

//#define SimTK_USE_STATIC_LIBRARIES

#include "simmath/internal/SimTKcpodes.h"

#include <cmath>
#include <cstdio>

using SimTK::Real;
using SimTK::Vector;
using SimTK::CPodes;
using std::printf;


class Problem1System : public SimTK::CPodesSystem {
public:
    // Override default implementations of these virtual functions.
    int explicitODE(Real t, const Vector& y, Vector& fout) const;
    int implicitODE(Real t, const Vector& y, const Vector& yp, 
                    Vector& fout) const;
};


class Problem2System : public SimTK::CPodesSystem {
public:
    // Override default implementations of these virtual functions.
    int explicitODE(Real t, const Vector& y, Vector& fout) const;
    int implicitODE(Real t, const Vector& y, const Vector& yp, 
                    Vector& fout) const;
};

static const CPodes::ODEType odeType = CPodes::ExplicitODE;

/* Shared Problem Constants */

#define ATOL Real(1.0e-6)   /* 1.0e-6 */
#define RTOL Real(0.0)      /* 0.0 */

#define ZERO   Real(0.0)
#define ONE    Real(1.0)
#define TWO    Real(2.0)
#define THIRTY Real(30.0)

/* Problem #1 Constants */

#define P1_NEQ        2
#define P1_ETA        Real(3.0)
#define P1_NOUT       4
#define P1_T0         Real(0.0)
#define P1_T1         Real(1.39283880203)
#define P1_DTOUT      Real(2.214773875)

/* Problem #2 Constants */

#define P2_MESHX      5
#define P2_MESHY      5
#define P2_NEQ        P2_MESHX*P2_MESHY
#define P2_ALPH1      Real(1.0)
#define P2_ALPH2      Real(1.0)
#define P2_NOUT       5
#define P2_ML         5
#define P2_MU         0
#define P2_T0         Real(0.0)
#define P2_T1         Real(0.01)
#define P2_TOUT_MULT  Real(10.0)

/* Private Helper Functions */
static void Problem1(void);
static void Problem2(void);
static Real MaxError(const Vector& y, Real t);
static void PrintFinalStats(CPodes&);

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
  Vector y(P1_NEQ), yp(P1_NEQ);
  Real reltol=RTOL, abstol=ATOL, t, tout, hu;
  int flag, iout, qu;

  y[0] = TWO;
  y[1] = ZERO;

  Problem1System sys1;

  // Use the explicit function to initialize derivatives for the implicit one.
  if (odeType == CPodes::ImplicitODE)
    sys1.explicitODE(P1_T0, y, yp);

  CPodes cpode(odeType, CPodes::Adams, CPodes::Functional);
  /*  flag = cpode.setInitStep(4.0e-9);*/
  flag = cpode.init(sys1, P1_T0, y, yp, CPodes::ScalarScalar, reltol, &abstol);

  printf("\n     t           x              xdot         qu     hu \n");
  for(iout=1, tout=P1_T1; iout <= P1_NOUT; iout++, tout += P1_DTOUT) {
    flag = cpode.step(tout, &t, y, yp, CPodes::Normal);
    if (flag != CPodes::Success)  break;
    flag = cpode.getLastOrder(&qu);
    flag = cpode.getLastStep(&hu);    
    printf("%10.5f    %12.5le   %12.5le   %2d    %6.4le\n", t, y[0], y[1], qu, hu);
  }
  
  PrintFinalStats(cpode);

  return;
}

int Problem1System::implicitODE(Real t, const Vector& y, const Vector& yp, 
                                Vector& res) const
{
  const int ret = Problem1System::explicitODE(t,y,res);
  res = yp - res;
  return ret;
}


int Problem1System::explicitODE(Real t, const Vector& y, Vector& ydot) const
{
  Real y0, y1;
  
  y0 = y[0];
  y1 = y[1];

  ydot[0] = y1;
  ydot[1] = (ONE - y0*y0)* P1_ETA * y1 - y0;

  return(0);
} 



static void Problem2(void)
{
  Vector y(P2_NEQ), yp(P2_NEQ);
  Real reltol=RTOL, abstol=ATOL, t, tout, erm, hu;
  int flag, iout, qu;

  y = 0;
  y[0] = ONE;

  Problem2System sys2;

  // Use the explicit function to initialize derivatives for the implicit one.
  if (odeType == CPodes::ImplicitODE)
    sys2.explicitODE(P2_T0, y, yp);

  CPodes cpode(odeType, CPodes::Adams, CPodes::Functional);      
  /*  flag = cpode.setInitStep(2.0e-9);*/
  flag = cpode.init(sys2, P2_T0, y, yp, CPodes::ScalarScalar, reltol, &abstol);
  
  printf("\n      t        max.err      qu     hu \n");
  for(iout=1, tout=P2_T1; iout <= P2_NOUT; iout++, tout*=P2_TOUT_MULT) {
    flag = cpode.step(tout, &t, y, yp, CPodes::Normal);
    if (flag != CPodes::Success) break;
    erm = MaxError(y, t);
    flag = cpode.getLastOrder(&qu);
    flag = cpode.getLastStep(&hu);
    printf("%10.3f  %12.4le   %2d   %12.4le\n", t, erm, qu, hu);
  }
  
  PrintFinalStats(cpode);

  return;
}

int Problem2System::implicitODE(Real t, const Vector& y, const Vector& yp, 
                                Vector& res) const
{
  const int ret = Problem2System::explicitODE(t, y, res);
  res = yp - res;
  return ret;
}

int Problem2System::explicitODE(Real t, const Vector& y, Vector& ydot) const
{
  long int i, j, k;
  Real d;

  /*
     Excluding boundaries, 

     ydot    = f    = -2 y    + alpha1 * y      + alpha2 * y
         i,j    i,j       i,j             i-1,j             i,j-1
  */

  for (j=0; j < P2_MESHY; j++) {
    for (i=0; i < P2_MESHX; i++) {
      k = i + j * P2_MESHX;
      d = -TWO*y[k];
      if (i != 0) d += P2_ALPH1 * y[k-1];
      if (j != 0) d += P2_ALPH2 * y[k-P2_MESHX];
      ydot[k] = d;
    }
  }

  return(0);
}

static Real MaxError(const Vector& y, Real t)
{
  Real er, ex=ZERO, yt, maxError=ZERO, ifact_inv, jfact_inv=ONE;
  
  if (t == ZERO) return(ZERO);

  if (t <= THIRTY) ex = std::exp(-TWO*t); 
  
  for (int j = 0; j < P2_MESHY; j++) {
    ifact_inv = ONE;
    for (int i = 0; i < P2_MESHX; i++) {
      const int k = i + j * P2_MESHX;
      yt = std::pow(t,i+j) * ex * ifact_inv * jfact_inv;
      er = std::abs(y[k] - yt);
      if (er > maxError) maxError = er;
      ifact_inv /= (i+1);
    }
    jfact_inv /= (j+1);
  }
  return(maxError);
}



static void PrintFinalStats(CPodes& cpode)
{
  int nst, nfe, nni, ncfn, netf;
  Real h0u;
  int flag;
  
  flag = cpode.getActualInitStep(&h0u);
  flag = cpode.getNumSteps(&nst);
  flag = cpode.getNumFctEvals(&nfe);
  flag = cpode.getNumErrTestFails(&netf);
  flag = cpode.getNumNonlinSolvIters(&nni);
  flag = cpode.getNumNonlinSolvConvFails(&ncfn);

  printf("\n Final statistics:\n\n");
  printf(" Number of steps                          = %4d \n", nst);
  printf(" Number of f-s                            = %4d \n", nfe);
  printf(" Number of nonlinear iterations           = %4d \n", nni);
  printf(" Number of nonlinear convergence failures = %4d \n", ncfn);
  printf(" Number of error test failures            = %4d \n", netf);
  printf(" Initial step size                        = %g \n\n", h0u);

}

