#ifndef SimTK_SIMMATH_CPODES_INTEGRATOR_REP_H_
#define SimTK_SIMMATH_CPODES_INTEGRATOR_REP_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmath(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
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
#include "SimTKcpodes.h"

#include "simmath/internal/common.h"
#include "simmath/Integrator.h"

#include "IntegratorRep.h"

namespace SimTK {

class CPodesIntegratorRep : public IntegratorRep {
public:
    CPodesIntegratorRep(Integrator* handle, const System& sys, CPodes::LinearMultistepMethod method);
    CPodesIntegratorRep(Integrator* handle, const System& sys, CPodes::LinearMultistepMethod method, CPodes::NonlinearSystemIterationType iterationType);
    ~CPodesIntegratorRep();
    void methodInitialize(const State&);
    void methodReinitialize(Stage stage, bool shouldTerminate);
    Integrator::SuccessfulStepStatus stepTo(Real reportTime, Real scheduledEventTime);
    Real getActualInitialStepSizeTaken() const;
    Real getPreviousStepSizeTaken() const;
    Real getPredictedNextStepSize() const;
    int getNumStepsAttempted() const;
    int getNumStepsTaken() const;
    int getNumErrorTestFailures() const;
    int getNumConvergenceTestFailures() const;
    int getNumConvergentIterations() const 
       {SimTK_ASSERT_ALWAYS(false, "CPodesIntegratorRep::getNumConvergentIterations(): not implemented");}
    int getNumDivergentIterations() const
       {SimTK_ASSERT_ALWAYS(false, "CPodesIntegratorRep::getNumDivergentIterations(): not implemented");}
    int getNumIterations() const;
    void resetMethodStatistics();
    void createInterpolatedState(Real t);
    void initializeIntegrationParameters();
    void reconstructForNewModel();
    const char* getMethodName() const;
    int getMethodMinOrder() const;
    int getMethodMaxOrder() const;
    bool methodHasErrorControl() const;
    void setUseCPodesProjection();
    void setOrderLimit(int order);
    class CPodesSystemImpl;
    friend class CPodesSystemImpl;
private:
    CPodes* cpodes;
    CPodesSystemImpl* cps;
    bool initialized, useCpodesProjection;
    int statsStepsTaken, statsErrorTestFailures, statsConvergenceTestFailures;
    int statsIterations;
    int pendingReturnCode;
    Real previousStartTime, previousTimeReturned;
    Vector savedY;
    CPodes::LinearMultistepMethod method;
    void init(CPodes::LinearMultistepMethod method, CPodes::NonlinearSystemIterationType iterationType);
};

} // namespace SimTK

#endif // SimTK_SIMMATH_CPODES_INTEGRATOR_REP_H_
