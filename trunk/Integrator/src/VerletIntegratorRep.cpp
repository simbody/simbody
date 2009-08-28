/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmath(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007-2008 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
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

#include "SimTKcommon.h"
#include "VerletIntegratorRep.h"

using namespace SimTK;

VerletIntegratorRep::VerletIntegratorRep(Integrator* handle, const System& sys) : AbstractIntegratorRep(handle, sys, 2, 3, "Verlet",  true) {
}

bool VerletIntegratorRep::attemptAStep(Real t0, Real t1, 
                                       const Vector& q0, const Vector& qdot0, const Vector& qdotdot0, 
                                       const Vector& u0, const Vector& udot0, 
                                       const Vector& z0, const Vector& zdot0, 
                                       Vector& yErrEst, int& errOrder, int& numIterations)
{
    // We will catch any exceptions thrown by realize() or project() and simply treat that
    // as a failure to take a step due to the step size being too big. The idea is that the
    // caller should reduce the step size and try again, giving up only when the step size
    // goes below the allowed minimum.

  try
  {
    statsStepsAttempted++;
    errOrder = 3;
    numIterations = 0;
    const Real h = t1-t0;
    State& advanced = updAdvancedState();

    const int nq = advanced.getNQ();
    const int nu = advanced.getNU();
    const int nz = advanced.getNZ();

    VectorView qErrEst = yErrEst(    0, nq);
    VectorView uErrEst = yErrEst(   nq, nu);
    VectorView zErrEst = yErrEst(nq+nu, nz);

    Vector dummyErrEst; // use this when we don't want the error estimate projected
    
    // Calculate the new positions q (3rd order) and initial (1st order) 
    // estimate for the velocities u and auxiliary variables z.
    
    // These are final values (the q's will get projected, though).
    advanced.updTime() = t1;
    advanced.updQ()    = q0 + qdot0*h + 0.5*qdotdot0*h*h;

    const Vector u1_est = u0 + udot0*h;
    const Vector z1_est = z0 + zdot0*h;

    advanced.updU()    = u1_est; // these will change
    advanced.updZ()    = z1_est;

    // Here we'll project the q's to their final constrained values, and project
    // our estimated u's also. TODO: does this mess up the error estimation?

    getSystem().realize(advanced, Stage::Velocity);
    if (userProjectEveryStep == 1 || IntegratorRep::calcWeightedRMSNorm(advanced.getYErr(), getDynamicSystemOneOverTolerances()) > consTol)
        projectStateAndErrorEstimate(advanced, dummyErrEst);

    // Get new values for the derivatives.
    realizeStateDerivatives(advanced);
    
    // We're going to integrate the u's and z's with the 2nd order implicit
    // trapezoid rule: u(t+h) = u(t) + h*(f(u(t))+f(u(t+h)))/2. Unfortunately this is an
    // implicit method so we have to iterate to refine u(t+h) until this equation
    // is acceptably satisfied. We're using functional iteration here which has a
    // very limited radius of convergence.
    
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
        advanced.setU(u0 + 0.5*h*(udot0 + udot1));
        advanced.setZ(z0 + 0.5*h*(zdot0 + zdot1));

        // Project only the velocities; we're done with the q's for good.
        getSystem().realize(advanced, Stage::Velocity);
        if (userProjectEveryStep == 1 || IntegratorRep::calcWeightedRMSNorm(advanced.getYErr(), getDynamicSystemOneOverTolerances()) > consTol)
            projectStateAndErrorEstimate(advanced, dummyErrEst, System::ProjectOptions::VelocityOnly);

        // Calculate fresh derivatives UDot and ZDot.
        realizeStateDerivatives(advanced);

        // Calculate convergence as the ratio of the norm of the last delta to
        // the norm of the values prior to the last change. We're using the 2-norm
        // but this ratio would be the same if we used the RMS norm. TinyReal is
        // there to keep us out of trouble if we started at zero.
        
        const Real convergenceU = (advanced.getU()-usave).norm()/(usave.norm()+TinyReal);
        const Real convergenceZ = (advanced.getZ()-zsave).norm()/(zsave.norm()+TinyReal);
        const Real change = std::max(convergenceU,convergenceZ);
        converged = (change <= tol);
        if (i > 1 && (change > prevChange))
            break; // we're headed in the wrong direction after two iterations -- give up
        prevChange = change;
    }

    // Now that we have achieved 2nd order estimates of u and z, we can use them to calculate a 3rd order
    // error estimate for q and 2nd order error estimates for u and z. Note that we have already
    // realized the state with the new values, so QDot reflects the new u's.

    qErrEst = q0 + 0.5*h*(qdot0+advanced.getQDot()) // implicit trapezoid rule integral
              - advanced.getQ();                    // Verlet integral

    uErrEst = u1_est                    // explicit Euler integral
              - advanced.getU();        // implicit trapezoid rule integral
    zErrEst = z1_est - advanced.getZ(); // ditto for z's

    // TODO: because we're only projecting velocities here, we aren't going to get our
    // position errors reduced here, which is a shame. Should be able to do this even
    // though we had to project q's earlier, because the projection matrix should still
    // be around.
    if (userProjectEveryStep == 1 || IntegratorRep::calcWeightedRMSNorm(advanced.getYErr(), getDynamicSystemOneOverTolerances()) > consTol)
        projectStateAndErrorEstimate(advanced, yErrEst, System::ProjectOptions::VelocityOnly);
    
    // Two different integrators were used to estimate errors: trapezoidal for Q, and
    // explicit Euler for U and Z.  This means that the U and Z errors are of a different
    // order than the Q errors.  We therefore multiply them by h so everything will be
    // of the same order.
    
    uErrEst *= h; zErrEst *= h; // everything is 3rd order in h now
    return converged;
  }
  catch (std::exception ex) {
    return false;
  }
}
