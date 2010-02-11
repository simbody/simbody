/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmath(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-10 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

/** @file
 * This is the private (library side) implementation of the 
 * RungeKuttaMersonIntegrator and RungeKuttaMersonIntegratorRep class which 
 * is a concrete class implementing the abstract IntegratorRep.
 */

#include "SimTKcommon.h"
#include "simmath/Integrator.h"
#include "simmath/RungeKuttaMersonIntegrator.h"

#include "IntegratorRep.h"
#include "RungeKuttaMersonIntegratorRep.h"

#include <exception>
#include <limits>

using namespace SimTK;

//------------------------------------------------------------------------------
//                     RUNGE KUTTA MERSON INTEGRATOR
//------------------------------------------------------------------------------

RungeKuttaMersonIntegrator::RungeKuttaMersonIntegrator(const System& sys) 
{
    rep = new RungeKuttaMersonIntegratorRep(this, sys);
}

RungeKuttaMersonIntegrator::~RungeKuttaMersonIntegrator() {
    delete rep;
}



//------------------------------------------------------------------------------
//                   RUNGE KUTTA MERSON INTEGRATOR REP
//------------------------------------------------------------------------------

RungeKuttaMersonIntegratorRep::RungeKuttaMersonIntegratorRep
   (Integrator* handle, const System& sys) 
:   AbstractIntegratorRep(handle, sys, 4, 4, "RungeKuttaMerson",  true) {
}

// For a discussion of the Runge-Kutta-Merson method, see Hairer,
// Norsett & Wanner, Solving ODEs I, 2nd rev. ed. pp. 166-8, and table 4.1
// on page 167. This is the Butcher diagram:
//
//        0|
//      1/3|  1/3
//      1/3|  1/6  1/6
//      1/2|  1/8   0   3/8
//        1|  1/2   0  -3/2   2
//       --|----------------------------
//        1|  1/6   0    0   2/3  1/6       propagated 4th order solution
//       --|----------------------------
//        1|  1/10  0   3/10 2/5  1/5  0    embedded 3rd order solution
//
// This is a 5-stage, first-same-as-last (FSAL) 4th order method which gives 
// us an embedded 3rd order method as well, so we can extract a 4th-order 
// error estimate for the 3rd-order result, which error estimate can then be 
// used for step size control, since it will behave as h^4. We then propagate 
// the 4th order result (whose error is unknown), which Hairer calls "local 
// extrapolation". We call the initial state (t0,y0) and want (t0+h,y1). We 
// are given the initial derivative f0=f(t0,y0), which most likely is left 
// over from an evaluation at the end of the last step.
// 
// We will call the derivatives at stage f1,f2,f3,f4 but these are done with 
// only two temporaries fa and fb. (What we're calling "f" Hairer calls "k".)
bool RungeKuttaMersonIntegratorRep::attemptAStep(Real t0, Real t1, 
                  const Vector& q0, const Vector& qdot0, const Vector& qdotdot0, 
                  const Vector& u0, const Vector& udot0, const Vector& z0, 
                  const Vector& zdot0, Vector& y1err, int& errOrder, int& numIterations)
{
    assert(t1 > t0);

    statsStepsAttempted++;
    errOrder = 4;
    const Vector& y0 = getPreviousY();
    const Vector& f0 = getPreviousYDot();
    if (ytmp[0].size() != y0.size())
        for (int i=0; i<NTemps; ++i)
            ytmp[i].resize(y0.size());
    Vector& ysave = ytmp[0]; // rename temps
    Vector& fa    = ytmp[1];
    Vector& fb    = ytmp[2];

    const Real h = t1-t0;

  try
  { setAdvancedStateAndRealizeDerivatives(t0+h/3, y0 + (h/3)*f0);
    fa = getAdvancedState().getYDot(); // fa=f1

    setAdvancedStateAndRealizeDerivatives(t0+h/3, y0 + (h/6)*(f0+fa)); // f0+f1
    fa = getAdvancedState().getYDot(); // fa=f2

    setAdvancedStateAndRealizeDerivatives(t0+h/2, y0 + (h/8)*(f0 + 3*fa)); // f0+3f2
    fb = getAdvancedState().getYDot(); // fb=f3

    // We'll need this for error estimation.
    ysave = y0 + (h/2)*(f0 - 3*fa + 4*fb); // f0-3f2+4f3
    setAdvancedStateAndRealizeDerivatives(t1, ysave);
    fa = getAdvancedState().getYDot(); // fa=f4

    // Final value. This is the 4th order accurate estimate for y1=y(t0+h)+O(h^5):
    // y1 = y0 + (h/6)*(f0 + 4 f3 + f4). Evaluate through kinematics only.
    setAdvancedStateAndRealizeKinematics(t1, y0 + (h/6)*(f0 + 4*fb + fa));
    // YErr is valid now
  } catch (...) { 
    return false; 
  }

    // This is an embedded 3rd-order estimate y1hat=y(t0+h)+O(h^4). (Apparently
    // Merson thought it was 5th order, but that is only true if
    // the function is linear w/constant coefficients; not bloody likely!)
    //     y1hat = y0 + (h/10)*(f0 + 3 f2 + 4 f3 + 2 f4)
    //
    // We don't actually have any need for y1hat, just its 4th-order
    // error estimate y1hat-y1=(1/5)(y1-ysave) (easily verified from the above).

    const Vector& y1 = getAdvancedState().getY();
    for (int i=0; i<y1.size(); ++i)
        y1err[i] = 0.2*std::abs(y1[i]-ysave[i]);

    if (userProjectEveryStep != 1) {
        const Real constraintError = 
            IntegratorRep::calcWeightedRMSNorm(getAdvancedState().getYErr(),
                                               getDynamicSystemOneOverTolerances());
        if (constraintError <= consTol)
            return true; // no need to project
    }

    // Project back to manifold and reduce error estimate appropriately. This
    // requires only kinematic evaluations, so doesn't count as a stage!
    try {projectStateAndErrorEstimate(updAdvancedState(), y1err);} 
    catch (...) {return false;}

    // Don't do final evaluation yet because we won't need it if
    // we fail the error test.
    numIterations = 1;
    return true;
}

