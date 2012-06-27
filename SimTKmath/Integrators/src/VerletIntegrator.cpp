/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-12 Stanford University and the Authors.        *
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

#include "simmath/VerletIntegrator.h"

#include "IntegratorRep.h"
#include "VerletIntegratorRep.h"

using namespace SimTK;

//==============================================================================
//                            VERLET INTEGRATOR
//==============================================================================

VerletIntegrator::VerletIntegrator(const System& sys) {
    rep = new VerletIntegratorRep(this, sys);
}

VerletIntegrator::VerletIntegrator(const System& sys, Real stepSize) {
    rep = new VerletIntegratorRep(this, sys);
    setFixedStepSize(stepSize);
}

VerletIntegrator::~VerletIntegrator() {
    delete rep;
}



//==============================================================================
//                          VERLET INTEGRATOR REP
//==============================================================================
VerletIntegratorRep::VerletIntegratorRep(Integrator* handle, const System& sys)
:   AbstractIntegratorRep(handle, sys, 2, 3, "Verlet",  true) 
{
}



//==============================================================================
//                            ATTEMPT DAE STEP
//==============================================================================
// Note that Verlet overrides the entire DAE step because it can't use the
// default ODE-then-DAE structure. Instead the constraint projections are
// interwoven here.
bool VerletIntegratorRep::attemptDAEStep
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

    const int nq = advanced.getNQ();
    const int nu = advanced.getNU();
    const int nz = advanced.getNZ();

    // We will catch any exceptions thrown by realize() or project() and simply
    // treat that as a failure to take a step due to the step size being too 
    // big. The idea is that the caller should reduce the step size and try 
    // again, giving up only when the step size goes below the allowed minimum.

  try
  {
    numIterations = 0;

    VectorView qErrEst  = yErrEst(    0, nq);       // all 3rd order estimates
    VectorView uErrEst  = yErrEst(   nq, nu);
    VectorView zErrEst  = yErrEst(nq+nu, nz);
    VectorView uzErrEst = yErrEst(   nq, nu+nz);    // all 2nd order estimates
    
    // Calculate the new positions q (3rd order) and initial (1st order) 
    // estimate for the velocities u and auxiliary variables z.
    
    // These are final values (the q's will get projected, though).
    advanced.updTime() = t1;
    advanced.updQ()    = q0 + h*qdot0 + (h*h/2)*qdotdot0;

    // Now make an initial estimate of first-order variable u and z.
    const Vector u1_est = u0 + h*udot0;
    const Vector z1_est = z0 + h*zdot0;

    advanced.updU() = u1_est; // u's and z's will change in advanced below
    advanced.updZ() = z1_est;

    system.realize(advanced, Stage::Time);
    system.prescribeQ(advanced);
    system.realize(advanced, Stage::Position);

    // Consider position constraint projection. (See AbstractIntegratorRep
    // for how we decide not to project.)
    const Real projectionLimit = 
        std::max(2*getConstraintToleranceInUse(), 
                    std::sqrt(getConstraintToleranceInUse()));

    bool anyChangesQ;
    if (!localProjectQAndQErrEstNoThrow(advanced, dummyErrEst, anyChangesQ, 
                                        projectionLimit))
        return false; // convergence failure for this step

    // q is now at its final integrated, prescribed, and projected value.
    // u and z still need refinement.

    system.prescribeU(advanced);
    system.realize(advanced, Stage::Velocity);

    // No u projection yet.

    // Get new values for the derivatives.
    realizeStateDerivatives(advanced);
    
    // We're going to integrate the u's and z's with the 2nd order implicit
    // trapezoid rule: u(t+h) = u(t) + h*(f(u(t))+f(u(t+h)))/2. Unfortunately 
    // this is an implicit method so we have to iterate to refine u(t+h) until
    // this equation is acceptably satisfied. We're using functional iteration 
    // here which has a very limited radius of convergence.
    
    const Real tol = std::min(1e-4, 0.1*getAccuracyInUse());
    Vector usave(nu), zsave(nz); // temporaries
    bool converged = false;
    Real prevChange = Infinity; // use this to quit early
    for (int i = 0; !converged && i < 10; ++i) {
        ++numIterations;
        // At this point we know that the advanced state has been realized
        // through the Acceleration level, so its uDot and zDot reflect
        // the u and z state values it contains.
        usave = advanced.getU();
        zsave = advanced.getZ();

        // Get these references now -- as soon as we change u or z they
        // will be invalid.
        const Vector& udot1 = advanced.getUDot();
        const Vector& zdot1 = advanced.getZDot();
        
        // Refine u and z estimates.
        advanced.setU(u0 + (h/2)*(udot0 + udot1));
        advanced.setZ(z0 + (h/2)*(zdot0 + zdot1));

        // Fix prescribed u's which may have been changed here.
        system.prescribeU(advanced);
        system.realize(advanced, Stage::Velocity);

        // No projection yet.

        // Calculate fresh derivatives UDot and ZDot.
        realizeStateDerivatives(advanced);

        // Calculate convergence as the ratio of the norm of the last delta to
        // the norm of the values prior to the last change. We're using the 
        // 2-norm but this ratio would be the same if we used the RMS norm. 
        // TinyReal is there to keep us out of trouble if we started at zero.
        
        const Real convergenceU = (advanced.getU()-usave).norm()
                                  / (usave.norm()+TinyReal);
        const Real convergenceZ = (advanced.getZ()-zsave).norm()
                                  / (zsave.norm()+TinyReal);
        const Real change = std::max(convergenceU,convergenceZ);
        converged = (change <= tol);
        if (i > 1 && (change > prevChange))
            break; // we're headed the wrong way after two iterations -- give up
        prevChange = change;
    }

    // Now that we have achieved 2nd order estimates of u and z, we can use 
    // them to calculate a 3rd order error estimate for q and 2nd order error 
    // estimates for u and z. Note that we have already realized the state with
    // the new values, so QDot reflects the new u's.

    qErrEst = q0 + (h/2)*(qdot0+advanced.getQDot()) // implicit trapezoid rule integral
              - advanced.getQ();                    // Verlet integral

    uErrEst = u1_est                    // explicit Euler integral
              - advanced.getU();        // implicit trapezoid rule integral
    zErrEst = z1_est - advanced.getZ(); // ditto for z's

    // TODO: because we're only projecting velocities here, we aren't going to 
    // get our position errors reduced here, which is a shame. Should be able 
    // to do this even though we had to project q's earlier, because the 
    // projection matrix should still be around.
    bool anyChangesU;
    if (!localProjectUAndUErrEstNoThrow(advanced, yErrEst, anyChangesU,
                                        projectionLimit))
        return false; // convergence failure for this step

    // Two different integrators were used to estimate errors: trapezoidal for 
    // Q, and explicit Euler for U and Z.  This means that the U and Z errors
    // are of a different order than the Q errors.  We therefore multiply them 
    // by h so everything will be of the same order.
    
    // TODO: (sherm) I don't think this is valid. Although it does fix the
    // order, it also changes the absolute errors (typically by reducing
    // them since h is probably < 1) which will affect whether the caller
    // decides to accept the step. Instead, a different error order should
    // be used when one of these is driving the step size.

    uErrEst *= h; zErrEst *= h; // everything is 3rd order in h now
    errOrder = 3;

    //errOrder = qErrRMS > uzErrRMS ? 3 : 2;

    return converged;
  }
  catch (std::exception ex) {
    return false;
  }
}

