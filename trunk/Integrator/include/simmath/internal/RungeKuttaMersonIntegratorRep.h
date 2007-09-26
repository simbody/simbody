#ifndef SimTK_SIMMATH_RUNGE_KUTTA_MERSON_INTEGRATOR_REP_H_
#define SimTK_SIMMATH_RUNGE_KUTTA_MERSON_INTEGRATOR_REP_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmath(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-7 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
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
 * Declares the RungeKuttaMerson concrete class which implements
 * an IntegratorRep.
 */

#include "SimTKcommon.h"
#include "simmath/Integrator.h"

#include "simmath/IntegratorRep.h"

#include <exception>
#include <limits>

using namespace SimTK;

class SimTK_SIMMATH_EXPORT RungeKuttaMersonIntegratorRep : public IntegratorRep {
public:
    RungeKuttaMersonIntegratorRep(Integrator* handle, const System& sys)
        : IntegratorRep(handle, sys)
    {
    }

    // virtuals

    Real getActualInitialStepSizeTaken() const {
        assert(initialized);
        return actualInitialStepSizeTaken;
    }
    Real getPreviousStepSizeTaken() const {
        assert(initialized);
        return previousSuccessfulStepSizeTaken;   
    }
    Real getPredictedNextStepSize() const {
        assert(initialized);
        return predictedNextStep;
    }

    long getNStepsAttempted() const {
        assert(initialized);
        return statsStepsAttempted;
    }
    long getNStepsTaken() const {
        assert(initialized);
        return statsStepsTaken;
    } 

    long getNErrorTestFailures() const {
        assert(initialized);
        return statsErrorTestFailures;
    } 

    void resetMethodStatistics() {
        statsStepsAttempted = 0;
        statsStepsTaken = 0;
        statsStepSizeChanges = 0;
        statsErrorTestFailures = 0;
    }

    // This is called from the parent class's initialize().
    void methodInitialize(const State&);

    // Throws an exception if anything goes wrong, otherwise takes a step
    // and returns status to indicate what happened.
    Integrator::SuccessfulStepStatus stepTo(Real reportTime, Real advanceLimit);

    const char* getMethodName() const;
    int getMethodMinOrder() const;
    int getMethodMaxOrder() const;
    bool methodHasErrorControl() const;

private:
    // Take a Runge-Kutta-Merson step from (t0,y0) to (t0+h,y1), including
    // constraint projection and error estimation. We are passed in an
    // already evaluated y0'=f0=f(t0,y0) from realizing (t0,y0) earlier.
    // Returns true only if all evaluations are successful, we make it to y1,
    // and are able to project to within tolerance. In that case 
    // y1err[i] is a 4th-order estimate of the absolute error in y1[i].
    // y1 itself is a 4th-order accurate result, meaning it has
    // a 5th-order error, likely much smaller than the estimate,
    // but we can't estimate its error directly because this is a 4(3)
    // method (4th order with embedded 3rd order).
    //
    // The currentState's time & continuous variables are trashed by
    // this routine. If successful, currentState=(t0+h,y1) at the end.
    bool takeAnRK4MStep(Real t0, const Vector& y0, const Vector& f0,
                        Real t1, Vector& y1err);

    void backUpAdvancedStateByInterpolation(Real t);
    void createInterpolatedState(Real t);
    
    enum AttemptedStepResult {
        SuccessWithoutEvent,
        SuccessWithEvent,
        ErrorTestFailure,  // need to reduce step size
        EvaluationFailure  // couldn't realize or project
    };

    AttemptedStepResult 
    attemptAStep(Real t0, const Vector& y0, const Vector& yd0,
                 Real t1, Real tReport, 
                 Real& scalarErrorEstimate, Real& tLow, Real& tHigh);

    // Adjust step size after a successful (errEst>=0) or failed (errEst<0) step. 
    void adjustStepSize(const Real& hTaken, bool hWasArtificiallyLimited, 
                        const Real& errEst, Real& hNeededForAccuracy);

private:    
    void initializeIntegrationParameters() {
        // must be realized to Stage::Instance
        initializeStepSizes();
        resetMethodStatistics();
        previousSuccessfulStepSizeTaken = NaN;
        predictedNextStep = initStepSize;

        // finalTime is infinity if unspecified
        if (userFinalTime != -1.) finalTime = userFinalTime;
        else finalTime = Infinity;

        // max inernal steps per stepTo() call is unlimited if unspecified;
        // we use "0" to mean unlimited.
        if (userInternalStepLimit != -1) 
            maxNumInternalSteps = userInternalStepLimit; // 0 means infinity
        else maxNumInternalSteps = 0;

        // Booleans

        // If the DynamicSystem says it needs control every step in order to 
        // simulate properly, we will still return even if the user explicitly
        // tells the integrator not to do so.
        if (userReturnEveryInternalStep != -1)
            returnEveryInternalStep = (   (userReturnEveryInternalStep==1) 
                                       || getDynamicSystemHasTimeAdvancedEvents());

        if (userProjectEveryStep != -1)
            projectEveryStep = (userProjectEveryStep==1);
        if (userAllowInterpolation != -1)
            allowInterpolation = (userAllowInterpolation==1);
        if (userProjectInterpolatedStates != -1)
            projectInterpolatedStates = (userProjectInterpolatedStates==1);
    }

    void initializeStepSizes() {
        const Real MinStep = std::pow(Eps, Real(0.75)); // e.g., 1e-12 in double
        const Real dsTimescale = getDynamicSystemTimescale();
        if (userMaxStepSize != -1.) { // got max
            maxStepSize = userMaxStepSize;
            if (userInitStepSize != -1.) { // got max & init
                initStepSize = userInitStepSize;
                if (userMinStepSize != -1.)
                    minStepSize = userMinStepSize; // got min,max,init
                else minStepSize = std::min(MinStep, initStepSize);
            } else { // got max, not init
                if (userMinStepSize != -1.) {
                    minStepSize = userMinStepSize; // max & min, not init
                    initStepSize = std::max(minStepSize, dsTimescale/10.);
                } else {
                    // got max only
                    initStepSize = std::min(dsTimescale/10., maxStepSize);
                    minStepSize  = std::min(MinStep, initStepSize);
                }
            }
        } else { // didn't get max
            maxStepSize = Infinity;
            if (userInitStepSize != -1.) { // got init, not max
                initStepSize = userInitStepSize;
                if (userMinStepSize != -1.) minStepSize = userMinStepSize;
                else minStepSize = std::min(MinStep, initStepSize);
            } else { // didn't get init or max
                if (userMinStepSize != -1.) { // got only min
                    minStepSize = userMinStepSize;
                    initStepSize = std::max(dsTimescale/10., minStepSize);
                } else { // didn't get anything
                    initStepSize = dsTimescale/10.;
                    minStepSize  = std::min(MinStep, initStepSize);
                }
            }
        }
    }

    static const int NTemps = 4;

    void reconstructForNewModel() {
        const int ny = getSystem().getDefaultState().getNY();
        const int nc = getSystem().getDefaultState().getNYErr();
        const int ne = getSystem().getDefaultState().getNEvents();

        initStepSize=minStepSize=maxStepSize
            =finalTime=predictedNextStep
            =previousSuccessfulStepSizeTaken
            =actualInitialStepSizeTaken=NaN;
    
        resetMethodStatistics();

        maxNumInternalSteps = -1;
        returnEveryInternalStep = getDynamicSystemHasTimeAdvancedEvents();
        projectEveryStep = false;
        allowInterpolation = true;
        projectInterpolatedStates = true;
        errEst.resize(ny);

        for (int i=0; i<NTemps; ++i)
            ytmp[i].resize(ny);
        haveTakenAStep = false;
        initialized = false;
    }

    Real   initStepSize, minStepSize, maxStepSize;
    Real   finalTime;
    int    maxNumInternalSteps;
    bool   returnEveryInternalStep, projectEveryStep;
    bool   allowInterpolation, projectInterpolatedStates;
    Vector errEst;
    Vector ytmp[NTemps];
    Real   previousSuccessfulStepSizeTaken, predictedNextStep;
    Real   actualInitialStepSizeTaken;

    bool haveTakenAStep;
    bool initialized;

    // Method stats (caution: not allowed to affect computation)
    long statsStepsAttempted;
    long statsStepsTaken;
    long statsStepSizeChanges;
    long statsErrorTestFailures;

};

#endif // SimTK_SIMMATH_RUNGE_KUTTA_MERSON_INTEGRATOR_REP_H_


