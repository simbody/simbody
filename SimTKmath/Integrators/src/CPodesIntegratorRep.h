#ifndef SimTK_SIMMATH_CPODES_INTEGRATOR_REP_H_
#define SimTK_SIMMATH_CPODES_INTEGRATOR_REP_H_

/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
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

#include "SimTKcommon.h"
#include "simmath/internal/common.h"
#include "simmath/Integrator.h"
#include "simmath/internal/SimTKcpodes.h"

#include "IntegratorRep.h"


namespace SimTK {

class CPodesIntegratorRep : public IntegratorRep {
public:
    CPodesIntegratorRep(Integrator* handle, const System& sys, CPodes::LinearMultistepMethod method);
    CPodesIntegratorRep(Integrator* handle, const System& sys, CPodes::LinearMultistepMethod method, CPodes::NonlinearSystemIterationType iterationType);
    ~CPodesIntegratorRep();
    void methodInitialize(const State&) override;
    void methodReinitialize(Stage stage, bool shouldTerminate) override;
    Integrator::SuccessfulStepStatus stepTo(Real reportTime, Real scheduledEventTime) override;
    Real getActualInitialStepSizeTaken() const override;
    Real getPreviousStepSizeTaken() const override;
    Real getPredictedNextStepSize() const override;
    int getNumStepsAttempted() const override;
    int getNumStepsTaken() const override;
    int getNumErrorTestFailures() const override;
    int getNumConvergenceTestFailures() const override;
    int getNumConvergentIterations() const override
       {SimTK_ASSERT_ALWAYS(false, "CPodesIntegratorRep::getNumConvergentIterations(): not implemented");}
    int getNumDivergentIterations() const override
       {SimTK_ASSERT_ALWAYS(false, "CPodesIntegratorRep::getNumDivergentIterations(): not implemented");}
    int getNumIterations() const override;
    void resetMethodStatistics() override;
    void createInterpolatedState(Real t);
    void initializeIntegrationParameters();
    void reconstructForNewModel();
    const char* getMethodName() const override;
    int getMethodMinOrder() const override;
    int getMethodMaxOrder() const override;
    bool methodHasErrorControl() const override;
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
