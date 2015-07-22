/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman, Michael Sherman                                    *
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


//------------------------------------------------------------------------------
//                   RUNGE KUTTA FELDBERG INTEGRATOR REP
//------------------------------------------------------------------------------

RungeKuttaFeldbergIntegratorRep::RungeKuttaFeldbergIntegratorRep
   (Integrator* handle, const System& sys)
:   AbstractIntegratorRep(handle, sys, 5, 5, "RungeKuttaFeldberg",  true) {}

bool RungeKuttaFeldbergIntegratorRep::attemptODEStep
   (Real t1, Vector& y1err, int& errOrder, int& numIterations)
{
    const Real C21    = Real( 1.0/4.0);
    const Real C22    = Real( 1.0/4.0);

    const Real C31    = Real( 3.0/8.0);
    const Real C32    = Real( 3.0/32.0);
    const Real C33    = Real( 9.0/32.0);

    const Real C41    = Real( 12.0/13.0);
    const Real C42    = Real( 1932.0/2197.0);
    const Real C43    = Real(-7200.0/2197.0);
    const Real C44    = Real( 7296.0/2197.0);

    const Real C51    = Real( 1.0);
    const Real C52    = Real( 439.0/216.0);
    const Real C53    = Real(-8.0);
    const Real C54    = Real( 3680.0/513.0);
    const Real C55    = Real(-845.0/4104.0);

    const Real C61    = Real( 1.0/2.0);
    const Real C62    = Real(-8.0/27.0);
    const Real C63    = Real( 2.0);
    const Real C64    = Real(-3544.0/2565.0);
    const Real C65    = Real( 1859.0/4104.0);
    const Real C66    = Real(-11.0/40.0);

    const double dCY1 =  25.0/216.0;
    const double dCY2 =  1408.0/2565.0;
    const double dCY3 =  2197.0/4104.0;
    const double dCY4 = -1.0/5.0;

    const Real CY1    = Real(dCY1);
    const Real CY2    = Real(dCY2);
    const Real CY3    = Real(dCY3);
    const Real CY4    = Real(dCY4);

    const Real CE1    = Real( 16.0/135.0-dCY1);
    const Real CE2    = Real( 6656.0/12825.0-dCY2);
    const Real CE3    = Real( 28561.0/56430.0-dCY3);
    const Real CE4    = Real(-9.0/50.0-dCY4);
    const Real CE5    = Real( 2.0/55.0);

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

    // Calculate the intermediate states.

    setAdvancedStateAndRealizeDerivatives(t0 + h*C21,
        y0 + h*C22*f0);
    ytmp[0] = getAdvancedState().getYDot();

    setAdvancedStateAndRealizeDerivatives(t0 + h*C31,
        y0 + h*C32*f0 + h*C33*ytmp[0]);
    ytmp[1] = getAdvancedState().getYDot();

    setAdvancedStateAndRealizeDerivatives(t0 + h*C41,
        y0 + h*C42*f0 + h*C43*ytmp[0] + h*C44*ytmp[1]);
    ytmp[2] = getAdvancedState().getYDot();

    setAdvancedStateAndRealizeDerivatives(t0 + h*C51,
        y0 + h*C52*f0 + h*C53*ytmp[0] + h*C54*ytmp[1] + h*C55*ytmp[2]);
    ytmp[3] = getAdvancedState().getYDot();

    setAdvancedStateAndRealizeDerivatives(t0 + h*C61,
        y0 + h*C62*f0 + h*C63*ytmp[0] + h*C64*ytmp[1] + h*C65*ytmp[2]
           + h*C66*ytmp[3]);
    ytmp[4] = getAdvancedState().getYDot();

    // Calculate the final state but don't evaluate the derivatives. That
    // would be a wasted stage since the caller will muck with the state before
    // the end of the step.
    setAdvancedStateAndRealizeKinematics(t1,
        y0 + h*CY1*f0 + h*CY2*ytmp[1] + h*CY3*ytmp[2] + h*CY4*ytmp[3]);
    // YErr is valid now, but not YDot.

    // Calculate the error estimate.
    y1err = h*CE1*f0 + h*CE2*ytmp[1] + h*CE3*ytmp[2] + h*CE4*ytmp[3]
                     + h*CE5*ytmp[4];

    return true;
}

