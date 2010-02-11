/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmath(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-10 Stanford University and the Authors.        *
 * Authors: Peter Eastman, Michael Sherman                                    *
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
 * RungeKuttaFeldbergIntegrator and RungeKuttaFeldbergIntegratorRep classes.
 */

#include "SimTKcommon.h"
#include "simmath/Integrator.h"
#include "simmath/RungeKuttaFeldbergIntegrator.h"

#include "IntegratorRep.h"
#include "RungeKuttaFeldbergIntegratorRep.h"

#include <exception>
#include <limits>

using namespace SimTK;

//------------------------------------------------------------------------------
//                     RUNGE KUTTA FELDBERG INTEGRATOR
//------------------------------------------------------------------------------

RungeKuttaFeldbergIntegrator::RungeKuttaFeldbergIntegrator(const System& sys) 
{
    rep = new RungeKuttaFeldbergIntegratorRep(this, sys);
}

RungeKuttaFeldbergIntegrator::~RungeKuttaFeldbergIntegrator() {
    delete rep;
}



//------------------------------------------------------------------------------
//                   RUNGE KUTTA FELDBERG INTEGRATOR REP
//------------------------------------------------------------------------------

RungeKuttaFeldbergIntegratorRep::RungeKuttaFeldbergIntegratorRep(Integrator* handle, const System& sys) : AbstractIntegratorRep(handle, sys, 5, 5, "RungeKuttaFeldberg",  true) {
}

bool RungeKuttaFeldbergIntegratorRep::attemptAStep
   (Real t0, Real t1, 
    const Vector& q0, const Vector& qdot0, const Vector& qdotdot0, 
    const Vector& u0, const Vector& udot0, const Vector& z0, 
    const Vector& zdot0, Vector& y1err, int& errOrder, int& numIterations)
{
    const double C21	=  1.0/4.0;
    const double C22	=  1.0/4.0;

    const double C31	=  3.0/8.0;
    const double C32	=  3.0/32.0;
    const double C33	=  9.0/32.0;

    const double C41	=  12.0/13.0;
    const double C42	=  1932.0/2197.0;
    const double C43	= -7200.0/2197.0;
    const double C44	=  7296.0/2197.0;

    const double C51	=  1.0;
    const double C52	=  439.0/216.0;
    const double C53	= -8.0;
    const double C54	=  3680.0/513.0;
    const double C55	= -845.0/4104.0;

    const double C61	=  1.0/2.0;
    const double C62	= -8.0/27.0;
    const double C63	=  2.0;
    const double C64	= -3544.0/2565.0;
    const double C65	=  1859.0/4104.0;
    const double C66	= -11.0/40.0;

    const double CY1	=  25.0/216.0;;
    const double CY2	=  1408.0/2565.0;
    const double CY3	=  2197.0/4104.0;
    const double CY4	= -1.0/5.0;

    const double CE1	=  16.0/135.0-CY1;
    const double CE2	=  6656.0/12825.0-CY2;
    const double CE3	=  28561.0/56430.0-CY3;;
    const double CE4	=  -9.0/50.0-CY4;
    const double CE5	=  2.0/55.0;

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

    try {
        // Calculate the intermediate states.
        
        setAdvancedStateAndRealizeDerivatives(t0 + h*C21, y0 + h*C22*f0);
        ytmp[0] = getAdvancedState().getYDot();

        setAdvancedStateAndRealizeDerivatives(t0 + h*C31, y0 + h*C32*f0 + h*C33*ytmp[0]);
        ytmp[1] = getAdvancedState().getYDot();

        setAdvancedStateAndRealizeDerivatives(t0 + h*C41, y0 + h*C42*f0 + h*C43*ytmp[0] + h*C44*ytmp[1]);
        ytmp[2] = getAdvancedState().getYDot();

        setAdvancedStateAndRealizeDerivatives(t0 + h*C51, y0 + h*C52*f0 + h*C53*ytmp[0] + h*C54*ytmp[1] + h*C55*ytmp[2]);
        ytmp[3] = getAdvancedState().getYDot();

        setAdvancedStateAndRealizeDerivatives(t0 + h*C61, y0 + h*C62*f0 + h*C63*ytmp[0] + h*C64*ytmp[1] + h*C65*ytmp[2] + h*C66*ytmp[3]);
        ytmp[4] = getAdvancedState().getYDot();
        
        // Calculate the final state.

        setAdvancedStateAndRealizeDerivatives(t1, y0 + h*CY1*f0 + h*CY2*ytmp[1] + h*CY3*ytmp[2] + h*CY4*ytmp[3]);
        
        // Calculate the error estimate.
        
        y1err = h*CE1*f0 + h*CE2*ytmp[1] + h*CE3*ytmp[2] + h*CE4*ytmp[3] + h*CE5*ytmp[4];
    } catch (...) { 
        return false; 
    }

    if (userProjectEveryStep != 1) {
        const Real constraintError = 
            IntegratorRep::calcWeightedRMSNorm(getAdvancedState().getYErr(), getDynamicSystemOneOverTolerances());
        if (constraintError <= consTol)
            return true; // no need to project
    }

    // Project back to manifold and reduce error estimate appropriately. This
    // requires only kinematic evaluations, so doesn't count as a stage!
    try {
        projectStateAndErrorEstimate(updAdvancedState(), y1err);
    } 
    catch (...) {
        return false;
    }

    numIterations = 1;
    return true;
}

