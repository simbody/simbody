/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
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

/** @file
 * This is the private (library side) implementation of the 
 * RungeKutta2Integrator and RungeKutta2IntegratorRep classes.
 */

#include "SimTKcommon.h"
#include "simmath/Integrator.h"
#include "simmath/RungeKutta2Integrator.h"

#include "IntegratorRep.h"
#include "RungeKutta2IntegratorRep.h"

#include <exception>
#include <limits>

using namespace SimTK;

//------------------------------------------------------------------------------
//                        RUNGE KUTTA 2 INTEGRATOR
//------------------------------------------------------------------------------

RungeKutta2Integrator::RungeKutta2Integrator(const System& sys) 
{
    rep = new RungeKutta2IntegratorRep(this, sys);
}


//------------------------------------------------------------------------------
//                      RUNGE KUTTA 2 INTEGRATOR REP
//------------------------------------------------------------------------------


RungeKutta2IntegratorRep::RungeKutta2IntegratorRep
   (Integrator* handle, const System& sys) 
:   AbstractIntegratorRep(handle, sys, 2, 2, "RungeKutta2",  true) {
}

// This is the explicit trapezoid rule, a Runge-Kutta 2(1) method. Here is
// the Butcher diagram:
//
//       0|
//       1|   1
//      --|--------------
//       1|  1/2  1/2       2nd order propagated solution
//      --|--------------
//        |   0    1    0   1st order Euler for error estimate
//
//
// This is a 2-stage, first-same-as-last (FSAL) 2nd order method which
// gives us an embedded 1st order method as well, so we can extract
// a 2nd-order error estimate for the 1st-order result, which error
// estimate can then be used for step size control, since it will
// behave as h^2. We then propagate the 2nd order result (whose error
// is unknown), which Hairer calls "local extrapolation".
// We call the initial state (t0,y0) and want (t0+h,y1). We are
// given the initial derivative f0=f(t0,y0), which most likely
// is left over from an evaluation at the end of the last step.

bool RungeKutta2IntegratorRep::attemptODEStep
   (Real t1, Vector& y1err, int& errOrder, int& numIterations)
{
    const Real t0 = getPreviousTime();
    assert(t1 > t0);

    statsStepsAttempted++;
    errOrder = 2;
    const Vector& y0 = getPreviousY();
    const Vector& f0 = getPreviousYDot();
    if (ytmp[0].size() != y0.size())
        for (int i=0; i<NTemps; ++i)
            ytmp[i].resize(y0.size());
    Vector& f1    = ytmp[0]; // rename temp

    const Real h = t1-t0;

    // First stage f1 = f(t1, y0+h*f0)
    setAdvancedStateAndRealizeDerivatives(t1, y0 + h*f0);
    f1 = getAdvancedState().getYDot();

    // Final value. This is the 2nd order accurate estimate for 
    // y1=y(t0+h)+O(h^3): y1 = y0 + (h/2)*(f0 + f1). 
    // Evaluate through kinematics only; it is a waste of a stage to 
    // evaluate derivatives here since the caller will muck with this before
    // the end of the step.
    setAdvancedStateAndRealizeKinematics(t1, y0 + (h/2)*(f0 + f1));
    // YErr is valid now

    // This is an embedded 1st-order estimate y1hat=y(t1)+O(h^2), with
    // y1hat = y0 + h*f1 (note that that is the the end-of-step
    // derivative). We just need the error y1-y1hat.

    const Vector& y1 = getAdvancedState().getY();
    for (int i=0; i<y1.size(); ++i)
        y1err[i] = std::abs(y1[i]-(y0[i] + h*f1[i]));

    return true;
}
