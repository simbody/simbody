/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-12 Stanford University and the Authors.        *
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
   (Real t1, Vector& y1err, int& errOrder, int& numIterations)
{
    const Real t0 = getPreviousTime();
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
