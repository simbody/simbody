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
#include "simmath/ExplicitEulerIntegrator.h"

int main () {
  try {
    PendulumSystem sys;
    sys.adoptEventHandler(new ZeroVelocityHandler(sys));
    sys.adoptEventHandler(PeriodicHandler::handler = new PeriodicHandler());
    sys.adoptEventHandler(new ZeroPositionHandler(sys));
    sys.adoptEventReporter(PeriodicReporter::reporter = new PeriodicReporter(sys));
    sys.adoptEventReporter(new OnceOnlyEventReporter());
    sys.adoptEventReporter(new DiscontinuousReporter());
    sys.realizeTopology();
    PeriodicHandler::handler->setEventInterval(0.01);
    PeriodicReporter::reporter->setEventInterval(0.015);
    
    // Test the integrator in both normal and single step modes.
    
    ExplicitEulerIntegrator integ(sys);
    testIntegrator(integ, sys, 1e-7);
    integ.setReturnEveryInternalStep(true);
    testIntegrator(integ, sys, 1e-7);
    cout << "Done" << endl;
    return 0;
  }
  catch (std::exception& e) {
    std::printf("FAILED: %s\n", e.what());
    return 1;
  }
}
