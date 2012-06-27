/* -------------------------------------------------------------------------- *
 *                     Simbody(tm) Example: Event Handler                     *
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

/* This example demonstrates the use of an event handler to implement 
discontinuous changes to the State during time stepping. */

#include "Simbody.h"

using namespace SimTK;

class BounceHandler : public TriggeredEventHandler {
public:
    BounceHandler() : TriggeredEventHandler(Stage::Position) {
        getTriggerInfo().setTriggerOnRisingSignTransition(false);
    }
    Real getValue(const State& state) const {
        return state.getQ()[0];
    }
    // Note: in general a discontinuous velocity change should be followed by
    // an impulse-momentum analysis to ensure that momentum is conserved. We're
    // not doing that here.
    void handleEvent(State& state, Real accuracy, bool& shouldTerminate) const {
        state.updU()[0] *= -1;
    }
};

int main() {
    
    // Create the system.
    
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid pendulumBody(MassProperties(1.0, Vec3(0), Inertia(1)));
    pendulumBody.addDecoration(Transform(), DecorativeSphere(0.1));
    MobilizedBody::Pin pendulum(matter.updGround(), Transform(Vec3(0)), pendulumBody, Transform(Vec3(0, 1, 0)));
   
    Visualizer viz(system);
    system.addEventReporter(new Visualizer::Reporter(viz, 1./30));
    
    system.addEventHandler(new BounceHandler());
    
    // Initialize the system and state.
    
    system.realizeTopology();
    State state = system.getDefaultState();
    pendulum.setOneU(state, 0, 1.0);
    
    // Simulate it.

    RungeKuttaMersonIntegrator integ(system);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(100.0);
}
