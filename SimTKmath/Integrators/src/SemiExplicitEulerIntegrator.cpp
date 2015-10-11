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

#include "simmath/SemiExplicitEulerIntegrator.h"

#include "IntegratorRep.h"
#include "SemiExplicitEulerIntegratorRep.h"

using namespace SimTK;

//==============================================================================
//                     SEMI EXPLICIT EULER INTEGRATOR
//==============================================================================

SemiExplicitEulerIntegrator::SemiExplicitEulerIntegrator
   (const System& sys, Real stepSize) {
    rep = new SemiExplicitEulerIntegratorRep(this, sys);
    setFixedStepSize(stepSize);
}


//==============================================================================
//                   SEMI EXPLICIT EULER INTEGRATOR REP
//==============================================================================
SemiExplicitEulerIntegratorRep::SemiExplicitEulerIntegratorRep
   (Integrator* handle, const System& sys)
:   AbstractIntegratorRep(handle, sys, 1, 1, "SemiExplicitEuler",  false)
{
}

//==============================================================================
//                         CREATE INTERPOLATED STATE
//==============================================================================
// Create an interpolated state at time t, which is between tPrev and tCurrent.
// If we haven't yet delivered an interpolated state in this interval, we have
// to initialize its discrete part from the advanced state.
void SemiExplicitEulerIntegratorRep::createInterpolatedState(Real t) {
    const System& system   = getSystem();
    const State&  advanced = getAdvancedState();
    State&        interp   = updInterpolatedState();
    interp = advanced; // pick up discrete stuff.
    const Real weight1 = (getAdvancedTime()-t) /
                         (getAdvancedTime()-getPreviousTime());
    const Real weight2 = 1-weight1;
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
void SemiExplicitEulerIntegratorRep::
backUpAdvancedStateByInterpolation(Real t) {
    const System& system   = getSystem();
    State& advanced = updAdvancedState();
    const Real t0 = getPreviousTime(), t1 = advanced.getTime();
    const Vector& y0 = getPreviousY();
    const Vector& y1 = advanced.getY();

    assert(t0 <= t && t <= t1);
    assert(t1 > t0);

    const Real weight0 = (t1-t) / (t1-t0);
    const Real weight1 = 1-weight0;
    advanced.updY() = weight0*y0 + weight1*y1;
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
// Note that SemiExplicitEuler overrides the entire DAE step because it can't
// use the default ODE-then-DAE structure. Instead the constraint projections
// are interwoven here.
bool SemiExplicitEulerIntegratorRep::attemptDAEStep
   (Real t1, Vector& yErrEst, int& errOrder, int& numIterations)
{
    const System& system   = getSystem();
    State& advanced = updAdvancedState();
    Vector dummyErrEst; // for when we don't want the error estimate projected

    statsStepsAttempted++;

    const Real    t0        = getPreviousTime();       // nicer names
    const Vector& q0        = getPreviousQ();
    const Vector& u0        = getPreviousU();
    const Vector& z0        = getPreviousZ();
    const Vector& qdot0     = getPreviousQDot();
    const Vector& udot0     = getPreviousUDot();
    const Vector& zdot0     = getPreviousZDot();
    const Vector& qdotdot0  = getPreviousQDotDot();

    const Real h = t1-t0;

    // Advance the first order variables.
    // TODO: this part should be implicit in u and z to make this symplectic
    // Euler.
    advanced.updZ() = getPreviousZ() + h * getPreviousZDot();
    advanced.updU() = getPreviousU() + h * getPreviousUDot();

    // Note that changing time does not invalidate position kinematics.
    advanced.updTime() = t1;
    system.realize(advanced, Stage::Position); // old q, new t
    system.prescribeU(advanced);

    // Update qdot_t1 = N(q_t0)*u_t1 from now-advanced u.
    const int nq = advanced.getNQ();
    Vector qdot_t1(nq);
    system.multiplyByN(advanced, advanced.getU(), qdot_t1);
    advanced.updQ() = getPreviousQ() + h * qdot_t1;
    system.prescribeQ(advanced);
    system.realize(advanced, Stage::Position); // new q, new t
    system.prescribeU(advanced); // update prescribed u in case q-dependent

    // Consider position constraint projection. Note that we have to do this
    // projection prior to calculating prescribed u's since the prescription
    // can depend on q's. Prevent project() from throwing an exception since
    // failure here may be recoverable.
    bool anyChanges;
    if (!localProjectQAndQErrEstNoThrow(advanced, dummyErrEst, anyChanges))
        return false; // convergence failure for this step

    // q's satisfy the position constraint manifold. Now work on u's.

    system.prescribeU(advanced);
    system.realize(advanced, Stage::Velocity);

    // Now try velocity constraint projection. Nothing will happen if
    // velocity constraints are already satisfied unless user has set the
    // ForceProjection option.

    if (!localProjectUAndUErrEstNoThrow(advanced, dummyErrEst, anyChanges))
        return false; // convergence failure for this step

    numIterations = 1;
    return true;
}

