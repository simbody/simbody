/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
 * Contributors: Michael Sherman                                              *
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

#include "simmath/ExplicitEulerIntegrator.h"

#include "IntegratorRep.h"
#include "ExplicitEulerIntegratorRep.h"

using namespace SimTK;


//------------------------------------------------------------------------------
//                        EXPLICIT EULER INTEGRATOR
//------------------------------------------------------------------------------

ExplicitEulerIntegrator::ExplicitEulerIntegrator(const System& sys) {
    rep = new ExplicitEulerIntegratorRep(this, sys);
}

ExplicitEulerIntegrator::ExplicitEulerIntegrator(const System& sys, Real stepSize) {
    rep = new ExplicitEulerIntegratorRep(this, sys);
    setInitialStepSize(stepSize);
}

ExplicitEulerIntegrator::~ExplicitEulerIntegrator() {
    delete rep;
}



//------------------------------------------------------------------------------
//                      EXPLICIT EULER INTEGRATOR REP
//------------------------------------------------------------------------------
ExplicitEulerIntegratorRep::ExplicitEulerIntegratorRep
   (Integrator* handle, const System& sys) 
:   AbstractIntegratorRep(handle, sys, 1, 1, "ExplicitEuler", true) 
{
}



//==============================================================================
//                         CREATE INTERPOLATED STATE
//==============================================================================
// Create an interpolated state at time t, which is between tPrev and tCurrent.
// If we haven't yet delivered an interpolated state in this interval, we have
// to initialize its discrete part from the advanced state.
void ExplicitEulerIntegratorRep::createInterpolatedState(Real t) {
    const System& system   = getSystem();
    const State&  advanced = getAdvancedState();
    State&        interp   = updInterpolatedState();
    interp = advanced; // pick up discrete stuff.
    const Real weight1 = (getAdvancedTime()-t) /
                         (getAdvancedTime()-getPreviousTime());
    const Real weight2 = 1.0-weight1;
    interp.updY() = weight1*getPreviousY()+weight2*getAdvancedState().getY();
    interp.updTime() = t;

    if (userProjectInterpolatedStates == 0) {
        system.realize(interp, Stage::Time);
        system.prescribeQ(interp);
        system.realize(interp, Stage::Position);
        system.prescribeU(interp);
        system.realize(interp, Stage::Velocity);
        return;
    }

    // We may need to project onto constraint manifold. Allow project()
    // to throw an exception if it fails since there is no way to recover here.
    realizeAndProjectKinematicsWithThrow(interp, ProjectOptions::LocalOnly);
}


//==============================================================================
//                  BACK UP ADVANCED STATE BY INTERPOLATION
//==============================================================================
// Interpolate the advanced state back to an earlier part of the interval,
// forgetting about the rest of the interval. This is necessary, for
// example after we have localized an event trigger to an interval tLow:tHigh
// where tHigh < tAdvanced.
void ExplicitEulerIntegratorRep::backUpAdvancedStateByInterpolation(Real t) {
    const System& system   = getSystem();
    State& advanced = updAdvancedState();

    assert(getPreviousTime() <= t && t <= advanced.getTime());
    advanced.updY() = getPreviousY() + (t-getPreviousTime())*getPreviousYDot();
    advanced.updTime() = t;

    // Ignore any user request not to project interpolated states here -- this
    // is the actual advanced state which will be propagated through the
    // rest of the trajectory so we can't allow it not to satisfy the 
    // constraints!
    // But it is OK if it just *barely* satisfies the constraints so we
    // won't get carried away if the user isn't being finicky about it.

    // Project position constraints if warranted. Allow project()
    // to throw an exception if it fails since there is no way to recover here.

    realizeAndProjectKinematicsWithThrow(advanced, ProjectOptions::LocalOnly);
}



//==============================================================================
//                            ATTEMPT DAE STEP
//==============================================================================
// Note that ExplicitEuler overrides the entire DAE step because it can't use
// the default ODE-then-DAE structure. Instead the constraint projections are
// interwoven here. The reason is that we need the end-of-step derivative
// value in order to be able to get an error estimate, and those derivatives
// must not be calculated until projection is done or we would do more than
// one function evaluation per step.
bool ExplicitEulerIntegratorRep::attemptDAEStep
   (Real t1, Vector& yErrEst, int& errOrder, int& numIterations)
{
    const System& system   = getSystem();
    State& advanced = updAdvancedState();

    statsStepsAttempted++;
    const Real h = t1 - getPreviousTime();

    // Take the step.
    advanced.updTime() = t1;
    advanced.updY()    = getPreviousY() + h*getPreviousYDot();
    yErrEst = advanced.getY(); // save unprojected Y for error estimate

    system.realize(advanced, Stage::Time);
    system.prescribeQ(advanced);
    system.realize(advanced, Stage::Position);

    // Consider position constraint projection. Note that we have to do this
    // projection prior to calculating prescribed u's since the prescription
    // can depend on q's. Prevent project() from throwing an exception since
    // failure here may be recoverable.
    Vector dummy; // no error estimate to project
    bool anyChanges;
    if (!localProjectQAndQErrEstNoThrow(advanced, dummy, anyChanges))
        return false; // convergence failure for this step

    // q's satisfy the position constraint manifold. Now work on u's.

    system.prescribeU(advanced);
    system.realize(advanced, Stage::Velocity);

    // Now try velocity constraint projection. Nothing will happen if
    // velocity constraints are already satisfied unless user has set the
    // ForceProjection option.

    if (!localProjectUAndUErrEstNoThrow(advanced, dummy, anyChanges))
        return false; // convergence failure for this step

    // Now calculate derivatives at the end of this interval/start of next
    // interval. TODO: prevent this from throwing
    try {realizeStateDerivatives(advanced);}
    catch (...) {return false;} // evaluation failed
    
    // Calculate a reference state with the explicit trapezoidal rule, and use
    // it to estimate error.
    //TODO: this is an odd mix of the unprojected Y and the projected YDot;
    //probably not right!
    yErrEst -= getPreviousY() + (h/2)*(getPreviousYDot()+advanced.getYDot());
    errOrder = 2;
    numIterations = 1;
    return true;
}


