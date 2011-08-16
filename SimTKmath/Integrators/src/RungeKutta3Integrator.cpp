/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmath(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
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
 * RungeKutta3Integrator and RungeKutta3IntegratorRep classes.
 */

#include "SimTKcommon.h"
#include "simmath/Integrator.h"
#include "simmath/RungeKutta3Integrator.h"

#include "IntegratorRep.h"
#include "RungeKutta3IntegratorRep.h"

#include <exception>
#include <limits>

using namespace SimTK;

//------------------------------------------------------------------------------
//                        RUNGE KUTTA 3 INTEGRATOR
//------------------------------------------------------------------------------

RungeKutta3Integrator::RungeKutta3Integrator(const System& sys) 
{
    rep = new RungeKutta3IntegratorRep(this, sys);
}

RungeKutta3Integrator::~RungeKutta3Integrator() {
    delete rep;
}


//------------------------------------------------------------------------------
//                      RUNGE KUTTA 3 INTEGRATOR REP
//------------------------------------------------------------------------------


RungeKutta3IntegratorRep::RungeKutta3IntegratorRep
   (Integrator* handle, const System& sys) 
:   AbstractIntegratorRep(handle, sys, 3, 3, "RungeKutta3",  true) {
}

// For a discussion of this Runge-Kutta 3(2) method, see J.C. Butcher, "The 
// Numerical Analysis of Ordinary Differential Equations", John Wiley & Sons,
// 1987, page 325. The embedded error estimate was derived using the method
// mentioned in Hairer, Norsett & Wanner, Solving ODEs I, 2nd rev. ed. on
// page 166. This is the Butcher diagram:
//
//           0|
//         1/2|  1/2
//           1|  -1    2
//          --|-------------------
//           1|  1/6  2/3  1/6       3rd order propagated solution
//          --|-------------------
//           1|   0    1    0    0   2nd order midpoint for error estimate
//
// This is a 3-stage, first-same-as-last (FSAL) 3rd order method which
// gives us an embedded 2nd order method as well, so we can extract
// a 3rd-order error estimate for the 2nd-order result, which error
// estimate can then be used for step size control, since it will
// behave as h^3. We then propagate the 3rd order result (whose error
// is unknown), which Hairer calls "local extrapolation".
// We call the initial state (t0,y0) and want (t0+h,y1). We are
// given the initial derivative f0=f(t0,y0), which most likely
// is left over from an evaluation at the end of the last step.

bool RungeKutta3IntegratorRep::attemptODEStep
   (Real t0, Real t1, 
    const Vector& q0, const Vector& qdot0, const Vector& qdotdot0, 
    const Vector& u0, const Vector& udot0, const Vector& z0, 
    const Vector& zdot0, Vector& y1err, int& errOrder, int& numIterations)
{
    assert(t1 > t0);

    statsStepsAttempted++;
    errOrder = 3;
    const Vector& y0 = getPreviousY();
    const Vector& f0 = getPreviousYDot();
    if (ytmp[0].size() != y0.size())
        for (int i=0; i<NTemps; ++i)
            ytmp[i].resize(y0.size());
    Vector& f1    = ytmp[0]; // rename temps
    Vector& f2    = ytmp[1];

    const Real h = t1-t0;

    setAdvancedStateAndRealizeDerivatives(t0+h/2, y0 + (h/2)*f0);
    f1 = getAdvancedState().getYDot();

    setAdvancedStateAndRealizeDerivatives(t1,     y0 + h*(2*f1-f0));
    f2 = getAdvancedState().getYDot();

    // Final value. This is the 3rd order accurate estimate for 
    // y1=y(t0+h)+O(h^4): y1 = y0 + (h/6)*(f0 + 4 f1 + f2). 
    // Evaluate through kinematics only; it is a waste of a stage to 
    // evaluate derivatives here since the caller will muck with this before
    // the end of the step.
    setAdvancedStateAndRealizeKinematics(t1,      y0 + (h/6)*(f0 + 4*f1 + f2));
    // YErr is valid now

    // This is an embedded 2nd-order estimate y1hat=y(t1)+O(h^3), with
    // y1hat = y0 + h*f1 (explicit midpoint method). We just need the
    // error y1-y1hat.

    const Vector& y1 = getAdvancedState().getY();
    for (int i=0; i<y1.size(); ++i)
        y1err[i] = std::abs(y1[i]-(y0[i] + h*f1[i]));

    return true;
}
