/* -------------------------------------------------------------------------- *
 *                      Simbody(tm) Example: Long Pendulum                    *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-13 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
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

/*                      Simbody ExampleLongPendulum
This example shows how to build a linked chain of bodies programmatically,
simulate it, and produce a simple animation while it is simulating. */

#include "Simbody.h"
#include <iostream>

using namespace SimTK;

int main() {
  try {    
    // Create the system.
    
    MultibodySystem system; system.setUseUniformBackground(true);
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid pendulumBody(MassProperties(1.0, Vec3(0), Inertia(1)));
    pendulumBody.addDecoration(Transform(), DecorativeSphere(0.1));

    MobilizedBody lastBody = matter.Ground();
    for (int i = 0; i < 10; ++i) {
        MobilizedBody::Ball pendulum(lastBody,     Transform(Vec3(0)), 
                                     pendulumBody, Transform(Vec3(0, 1, 0)));
        lastBody = pendulum;
    }

    Visualizer viz(system);
    system.addEventReporter(new Visualizer::Reporter(viz, 1./30));
    
    // Initialize the system and state.
    
    system.realizeTopology();
    State state = system.getDefaultState();
    Random::Gaussian random;
    for (int i = 0; i < state.getNQ(); ++i)
        state.updQ()[i] = random.getValue();
    
    // Simulate it.

    RungeKuttaMersonIntegrator integ(system);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(10.0);

  } catch(const std::exception& e) {
    std::cout << "EXCEPTION: " << e.what() << std::endl;
    return 1;
  } catch (...) {
      std::cout << "UNKNOWN EXCEPTION\n";
      return 1;
  }
    return 0;
}
