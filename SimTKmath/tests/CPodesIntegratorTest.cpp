/* -------------------------------------------------------------------------- *
 *                          Simbody(tm): SimTKmath                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman, Peter Eastman                                    *
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

#include "IntegratorTestFramework.h"
#include "simmath/CPodesIntegrator.h"

int main () {
  try {
    PendulumSystem sys;
    sys.addEventHandler(new ZeroVelocityHandler(sys));
    sys.addEventHandler(PeriodicHandler::handler = new PeriodicHandler());
    sys.addEventHandler(new ZeroPositionHandler(sys));
    sys.addEventReporter(PeriodicReporter::reporter = new PeriodicReporter(sys));
    sys.addEventReporter(new OnceOnlyEventReporter());
    sys.addEventReporter(new DiscontinuousReporter());
    sys.realizeTopology();

    // Test with various intervals for the event handler and event reporter, ones that are either
    // large or small compared to the expected internal step size of the integrator.

    for (int i = 0; i < 4; ++i) {
        PeriodicHandler::handler->setEventInterval(i == 0 || i == 1 ? 0.01 : 2.0);
        PeriodicReporter::reporter->setEventInterval(i == 0 || i == 2 ? 0.015 : 1.5);
        
        // Test the BDF integrator in both normal and single step modes.
        
        CPodesIntegrator bdfInteg(sys, CPodes::BDF);
        testIntegrator(bdfInteg, sys);
        bdfInteg.setReturnEveryInternalStep(true);
        testIntegrator(bdfInteg, sys);
        
        // Test the Adams integrator in both normal and single step modes.
        
        CPodesIntegrator adamsInteg(sys, CPodes::Adams);
        testIntegrator(adamsInteg, sys);
        adamsInteg.setReturnEveryInternalStep(true);
        testIntegrator(adamsInteg, sys);
        
        // Calling setUseCPodesProjection() after the integrator has been initialized should
        // produce an exception.
        
        try {
            adamsInteg.setUseCPodesProjection();
            assert(false);
        }
        catch (...) {
        }
        
        // Try having CPODES do the projection instead of the System.
        
        CPodesIntegrator projInteg(sys, CPodes::BDF);
        projInteg.setUseCPodesProjection();
        testIntegrator(projInteg, sys);
    }
    cout << "Done" << endl;
    return 0;
  }
  catch (std::exception& e) {
    std::printf("FAILED: %s\n", e.what());
    return 1;
  }
}
