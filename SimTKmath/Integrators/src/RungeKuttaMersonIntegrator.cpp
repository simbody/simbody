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
bool RungeKuttaMersonIntegratorRep::attemptODEStep
   (Real t1, Vector& y1err, int& errOrder, int& numIterations)
{
    const Real t0 = getPreviousTime();
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

    setAdvancedStateAndRealizeDerivatives(t0+h/3, y0 + (h/3)*f0);
    fa = getAdvancedState().getYDot(); // fa=f1

    setAdvancedStateAndRealizeDerivatives(t0+h/3, y0 + (h/6)*(f0+fa)); // f0+f1
    fa = getAdvancedState().getYDot(); // fa=f2

    setAdvancedStateAndRealizeDerivatives(t0+h/2, y0 + (h/8)*(f0 + 3*fa)); // f0+3f2
    fb = getAdvancedState().getYDot(); // fb=f3

    // We'll need this for error estimation.
    ysave = y0 + (h/2)*(f0 - 3*fa + 4*fb); // f0-3f2+4f3
    setAdvancedStateAndRealizeDerivatives(t1, ysave);
    fa = getAdvancedState().getYDot(); // fa=f4

    // Final value. This is the 4th order accurate estimate for 
    // y1=y(t0+h)+O(h^5): y1 = y0 + (h/6)*(f0 + 4 f3 + f4). 
    // Evaluate through kinematics only; it is a waste of a stage to 
    // evaluate derivatives here since the caller will muck with this before
    // the end of the step.
    setAdvancedStateAndRealizeKinematics(t1, y0 + (h/6)*(f0 + 4*fb + fa));
    // YErr is valid now

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

    return true;
}

