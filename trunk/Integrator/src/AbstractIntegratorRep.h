#ifndef SimTK_SIMMATH_ABSTRACT_INTEGRATOR_REP_H_
#define SimTK_SIMMATH_ABSTRACT_INTEGRATOR_REP_H_

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

#include "simmath/internal/common.h"
#include "simmath/Integrator.h"

#include "IntegratorRep.h"

namespace SimTK {

/**
 * This class implements most of the generic functionality needed for an Integrator,
 * leaving only the actual integration method to be implemented by the subclass.  This
 * is the parent class of several different integrators.
 */

class AbstractIntegratorRep : public IntegratorRep {
public:
    AbstractIntegratorRep(Integrator* handle, const System& sys, int minOrder, int maxOrder, const std::string& methodName, bool hasErrorControl);
    void methodInitialize(const State&);
    Integrator::SuccessfulStepStatus stepTo(Real reportTime, Real scheduledEventTime);
    Real getActualInitialStepSizeTaken() const;
    Real getPreviousStepSizeTaken() const;
    Real getPredictedNextStepSize() const;
    long getNumStepsAttempted() const;
    long getNumStepsTaken() const;
    long getNumErrorTestFailures() const;
    long getNumConvergenceTestFailures() const;
    long getNumConvergentIterations() const;
    long getNumDivergentIterations() const;
    long getNumIterations() const;
    void resetMethodStatistics();
    const char* getMethodName() const;
    int getMethodMinOrder() const;
    int getMethodMaxOrder() const;
    bool methodHasErrorControl() const;

protected:
    /**
     * Given initial values for all the continuous variables y=(q,u,z) and their derivatives
     * (not necessarily what's in advancedState currently), take a trial step of size 
     * h=(t1-t0), optimistically storing the result in advancedState.
     * Also estimate the absolute error in each element of y, and store them in yErrEst.
     * Returns true if the step converged (always true for non-iterative methods), false otherwise.
     * The number of internal iterations just for this step is return in numIterations, which
     * should always be 1 for non-iterative methods.
     */
    virtual bool attemptAStep(Real t0, Real t1, 
                      const Vector& q0, const Vector& qdot0, const Vector& qdotdot0, 
                      const Vector& u0, const Vector& udot0, const Vector& z0, 
                      const Vector& zdot0, Vector& yErrEst, int& errOrder, int& numIterations) = 0;
    /**
     * Evaluate the error that occurred in the step we just attempted, and select a new
     * step size accordingly.  The default implementation should work well for most integrators.
     * 
     * @param err     the error estimate from the step that was just attempted
     * @param hWasArtificiallyLimited   tells whether the step size was artificially
     * reduced due to a scheduled event time.  If this is true, we will never attempt
     * to increase the step size.
     * @return true if the step should be accepted, false if it should be rejected and retried
     * with a smaller step size
     */
    virtual bool adjustStepSize(Real err, int errOrder, bool hWasArtificiallyLimited);
    /**
     * Create an interpolated state at time t, which is between the previous and advanced
     * times.  The default implementation uses third order Hermite spline interpolation.
     */
    virtual void createInterpolatedState(Real t);
    /**
     * Interpolate the advanced state back to an earlier part of the interval,
     * forgetting about the rest of the interval. This is necessary, for example,
     * after we have localized an event trigger to an interval tLow:tHigh where
     * tHigh < tAdvanced.  The default implementation uses third order Hermite
     * spline interpolation.
     */
    virtual void backUpAdvancedStateByInterpolation(Real t);
    long statsStepsTaken, statsStepsAttempted, statsErrorTestFailures, statsConvergenceTestFailures;

    // Iterative methods should count iterations and then classify them as 
    // iterations that led to successful convergence and those that didn't.
    long statsConvergentIterations, statsDivergentIterations;
private:
    bool takeOneStep(Real tMax, Real tReport);
    bool initialized, hasErrorControl;
    Real currentStepSize, lastStepSize, actualInitialStepSizeTaken;
    int minOrder, maxOrder;
    std::string methodName;
};

} // namespace SimTK

#endif // SimTK_SIMMATH_ABSTRACT_INTEGRATOR_REP_H_
