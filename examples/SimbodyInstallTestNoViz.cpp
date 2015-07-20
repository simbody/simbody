/* -------------------------------------------------------------------------- *
 *             Simbody(tm) Example: Simbody Install Test No Viz               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2011-12 Stanford University and the Authors.        *
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

#include "Simbody.h"

#include <cstdio>
#include <exception>

using namespace SimTK;

/* This is the same as the example named Pendulum except that
it does not use the Visualizer. Instead it prints out the kinetic 
energy and total energy (the latter should be conserved). This is 
useful for checking that installation is OK on systems where the 
Visualizer wasn't built or isn't working. */


class EnergyReporter : public PeriodicEventReporter {
public:
    EnergyReporter(const MultibodySystem& system, Real interval) 
    :   PeriodicEventReporter(interval), system(system) {}

    void handleEvent(const State& s) const {
        std::cout << "t=" << s.getTime() 
             << "\tke=" << system.calcKineticEnergy(s) 
             << "\tE=" << system.calcEnergy(s)
             << std::endl;
    }
private:
    const MultibodySystem& system;
};

int main() {
  try
  { // Create the system.
    
    MultibodySystem         system;
    SimbodyMatterSubsystem  matter(system);
    GeneralForceSubsystem   forces(system);
    Force::UniformGravity   gravity(forces, matter, Vec3(0, -9.8, 0));

    Body::Rigid pendulumBody(MassProperties(1.0, Vec3(0), Inertia(1)));
    pendulumBody.addDecoration(Transform(), DecorativeSphere(0.1));
    MobilizedBody::Pin pendulum(matter.Ground(), Transform(Vec3(0)), 
                                pendulumBody,    Transform(Vec3(0, 1, 0)));

    system.adoptEventReporter(new EnergyReporter(system,1.));
   
    // Initialize the system and state.
    
    system.realizeTopology();
    State state = system.getDefaultState();
    pendulum.setOneU(state, 0, 1.0);

    
    // Simulate it.

    RungeKuttaMersonIntegrator integ(system);
    integ.setAccuracy(1e-6); // ask for *lots* of accuracy here (default is 1e-3)
    TimeStepper ts(integ);
    ts.initialize(state);
    ts.stepTo(100.0);

  } catch (const std::exception& e) {
    std::printf("EXCEPTION THROWN: %s\n", e.what());
    exit(1);

  } catch (...) {
    std::printf("UNKNOWN EXCEPTION THROWN\n");
    exit(1);
  }

    return 0;
}
