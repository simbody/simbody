/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmath(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-7 Stanford University and the Authors.         *
 * Authors: Michael Sherman, Peter Eastman                                    *
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
