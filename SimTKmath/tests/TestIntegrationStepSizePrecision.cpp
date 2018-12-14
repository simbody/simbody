/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-15 Stanford University and the Authors.        *
 * Authors: Artin Farahani                                                    *
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

/* Some Simbody integrators use checkStepSizePrecision() to throw an exception 
 * when their step size falls below the Real type's precision. The following 
 * is a test case for this functionality.
 */

#include "SimTKcommon.h"
#include "SimTKmath.h"
#include "simmath/TimeStepper.h"
#include "PendulumSystem.h"
#include "simmath/RungeKuttaMersonIntegrator.h"
#include "../Integrators/src/IntegratorRep.h"
using namespace SimTK;


void testIntegrator (Integrator& integ, PendulumSystem& sys, Real accuracy) {
    const Real t0 = 0.0;
    const Real qi[] = {1.0, 0.0};
    const Real ui[] = {0.0, 0.0};
    const Vector q0(2, qi);
    const Vector u0(2, ui);

    sys.setDefaultMass(10);
    sys.setDefaultTimeAndState(t0, q0, u0);
    integ.setAccuracy(accuracy);
    // This loose accuracy is needed for step size to reach the precision limit.
    integ.setConstraintTolerance(0.8);
    
    TimeStepper ts(sys);
    ts.setIntegrator(integ);
    ts.initialize(sys.getDefaultState());
    ts.stepTo(2.0);
}


// In this subtest, a very tight accuracy is intentionally chosen to cause the
// step size to reach its precision limit.  Therefore, we expect 
// checkStepSizePrecision() to catch this and throw an exception.
int testNoMinStepSize(PendulumSystem& sys, 
                      RungeKuttaMersonIntegrator& integ) {
    try {
        testIntegrator(integ, sys, 1e-35);
        return 1;
    }
    catch (Integrator::StepFailed& e) {
        std::printf("Success: As expected: %s\n", e.what());
        return 0;
    }
}


// In this subtest, the user specifies a min step size that is _smaller_ than
// the precision limit.  Therefore, the limit should still be reached and
// caught by checkStepSizePrecision().
int testSmallerMinStepSize(PendulumSystem& sys,
                           RungeKuttaMersonIntegrator& integ) {
    try {
        integ.setMinimumStepSize(1e-16);
        testIntegrator(integ, sys, 1e-35);
        return 1;
    }
    catch (Integrator::StepFailed& e) {
        std::printf("Success: As expected: %s\n", e.what());
        return 0;
    }
}


// In this subtest, the user specifies a min step size that is _larger_ than
// the precision limit.  Therefore, the limit should not be reached and not
// caught by checkStepSizePrecision().  We don't expect an exception to be 
// thrown here.
int testLargerMinStepSize(PendulumSystem& sys,
                          RungeKuttaMersonIntegrator& integ) {
    integ.setMinimumStepSize(1e-14);
    testIntegrator(integ, sys, 1e-35);
    return 0;
}


int main() {
    SimTK_START_TEST("TestIntegrationStepSizePrecision");

    PendulumSystem sys;
    sys.realizeTopology();
    RungeKuttaMersonIntegrator integ(sys);

    SimTK_SUBTEST2(testNoMinStepSize, sys, integ);
    SimTK_SUBTEST2(testSmallerMinStepSize, sys, integ);
    SimTK_SUBTEST2(testLargerMinStepSize, sys, integ);

    SimTK_END_TEST();
    return 0;
}

