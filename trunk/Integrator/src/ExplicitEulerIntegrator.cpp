/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmath(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007-10 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
 * Contributors: Michael Sherman                                              *
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

// Create an interpolated state at time t, which is between tPrev and tCurrent.
// If we haven't yet delivered an interpolated state in this interval, we have
// to initialize its discrete part from the advanced state.
void ExplicitEulerIntegratorRep::createInterpolatedState(Real t) {
    const State& advanced = getAdvancedState();
    State&       interp   = updInterpolatedState();
    interp = advanced; // pick up discrete stuff.
    const Real weight1 = (getAdvancedTime()-t) /
                         (getAdvancedTime()-getPreviousTime());
    const Real weight2 = 1.0-weight1;
    interp.updY() = weight1*getPreviousY()+weight2*getAdvancedState().getY();
    interp.updTime() = t;
    getSystem().realize(interp, Stage::Velocity); // cheap  
    if (userProjectInterpolatedStates == 0) // default is to project interpolated states if they need it
        return; // leave 'em in "as is" condition
    if (userProjectEveryStep != 1) {
        const Real constraintError =  IntegratorRep::calcWeightedRMSNorm(interp.getYErr(), getDynamicSystemOneOverTolerances());
        if (constraintError <= consTol)
            return; // no need to project
    }

    // No error estimate to project here; just pass an empty Vector. Local
    // projection only is allowed; we should be close to the manifold.
    projectStateAndErrorEstimate(interp, Vector());
}

// Note that ExplicitEuler overrides the entire DAE step because it can't use
// the default ODE-then-DAE structure. Instead the constraint projections are
// interwoven here. The reason is that we need the end-of-step derivative
// value in order to be able to get an error estimate, and those derivatives
// must not be calculated until projection is done or we would do more than
// one function evaluation per step.
bool ExplicitEulerIntegratorRep::attemptDAEStep
   (Real t0, Real t1, 
    const Vector& q0, const Vector& qdot0, const Vector& qdotdot0, 
    const Vector& u0, const Vector& udot0, 
    const Vector& z0, const Vector& zdot0, 
    Vector& yErrEst, int& errOrder, int& numIterations)
{
    statsStepsAttempted++;
    const Real h = t1 - t0;
    State& advanced = updAdvancedState();

    // Take the step.
    advanced.updTime() = t1;
    advanced.updY()    = getPreviousY() + h*getPreviousYDot();

    yErrEst = advanced.getY(); // save unprojected Y

    getSystem().realize(advanced, Stage::Velocity);
    const Real consErrAfterODE = 
        calcWeightedRMSNorm(advanced.getYErr(),
                            getDynamicSystemOneOverTolerances());

    if (   userProjectEveryStep == 1 
        || consErrAfterODE > getConstraintToleranceInUse()) 
    {
        try {
            projectStateAndErrorEstimate(advanced, Vector());
        } catch (...) {
            return false; // projection failed
        }
    }

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

// Interpolate the advanced state back to an earlier part of the interval,
// forgetting about the rest of the interval. This is necessary, for
// example after we have localized an event trigger to an interval tLow:tHigh
// where tHigh < tAdvanced.
void ExplicitEulerIntegratorRep::backUpAdvancedStateByInterpolation(Real t) {
    State& advanced = updAdvancedState();
    Vector yinterp;

    assert(getPreviousTime() <= t && t <= advanced.getTime());
    advanced.updY() = getPreviousY() + (t-getPreviousTime())*getPreviousYDot();
    advanced.updTime() = t;
    getSystem().realize(advanced, Stage::Velocity); // cheap 

    // Ignore any user request not to project interpolated states here -- this
    // is the actual advanced state which will be propagated through the
    // rest of the trajectory so we can't allow it not to satisfy the 
    // constraints!

    // But it is OK if it just *barely* satisfies the constraints so we
    // won't get carried away if the user isn't being finicky about it.
    if (userProjectEveryStep != 1) {
        const Real constraintError = 
            IntegratorRep::calcWeightedRMSNorm(advanced.getYErr(),
                                               getDynamicSystemOneOverTolerances());
        if (constraintError <= consTol)
            return; // no need to project
    }

    // No error estimate to project here; just pass an empty Vector. Local
    // projection only is allowed; we should be close to the manifold.
    projectStateAndErrorEstimate(advanced, Vector());
}


