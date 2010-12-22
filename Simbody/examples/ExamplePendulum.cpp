/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): Pendulum Example                       *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007-10 Stanford University and the Authors.        *
 * Authors: Ajay Seth                                                         *
 * Contributors: Michael Sherman                                              *
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

/* A one-body pendulum. Doesn't get much easier than this! */

#include "SimTKsimbody.h"

using namespace SimTK;

int main() {
    
    // Create the system, with subsystems for the bodies and some forces.
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);

    // Add gravity as a force element.
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));

    // Create the body and some artwork for it.
    Body::Rigid pendulumBody(MassProperties(1.0, Vec3(0), Inertia(1)));
    pendulumBody.addDecoration(Transform(), DecorativeSphere(0.1).setColor(Red));

    // Add an instance of the body to the multibody system by connecting
    // it to Ground via a pin mobilizer.
    MobilizedBody::Pin pendulum(matter.updGround(), Transform(Vec3(0)), 
                                pendulumBody, Transform(Vec3(0, 1, 0)));

    // Visualize with default options; ask for a report every 1/30 of a second
    // to match the Visualizer's default 30 frames per second rate.
    Visualizer viz(system);
    system.updDefaultSubsystem().addEventReporter
                                (new Visualizer::Reporter(viz, 1./30));
    
    // Initialize the system and state.
    
    system.realizeTopology();
    State state = system.getDefaultState();
    pendulum.setOneU(state, 0, 1.0); // initial velocity 1 rad/sec
    
    // Simulate it.

    RungeKuttaMersonIntegrator integ(system);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(100.0);
}
