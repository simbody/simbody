/* -------------------------------------------------------------------------- *
 *                        Simbody(tm) Example: Gears                          *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
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

using namespace SimTK;

int main() {

    // Create the system.

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);

    // Create bodies.
    Body::Rigid gearBody(MassProperties(1.0, Vec3(0), Inertia(1)));
    gearBody.addDecoration(Transform(Rotation(0.5*Pi, XAxis)),
        DecorativeCylinder(1.0, 0.1));
    Body::Rigid rodBody(MassProperties(1.0, Vec3(0), Inertia(1)));
    rodBody.addDecoration(Transform(Vec3(0, 1, 0)), DecorativeCylinder(0.05, 1.0));

    // Create instances of the bodies that are connected into the multibody
    // system via Mobilizers. Note that we use the gear body twice.

    MobilizedBody::Pin gear1(matter.updGround(), Transform(Vec3(1, 0, 0)),
                             gearBody, Transform());
    MobilizedBody::Pin gear2(matter.updGround(), Transform(Vec3(-1, 0, 0)),
                             gearBody, Transform());
    MobilizedBody::Pin rod(gear2, Transform(Vec3(0, 0.8, 0.1)), rodBody, Transform());

    // Add constraints.
    Constraint::ConstantSpeed(gear1, 2*Pi); // i.e., 1 rotation per second
    Constraint::NoSlip1D(matter.updGround(), Vec3(0), UnitVec3(0, 1, 0), gear1, gear2);

    // We want the rod end point traveling along a line. We'll draw part of the
    // line to make it clear.
    Constraint::PointOnLine(matter.updGround(), UnitVec3(0, 1, 0), Vec3(0, 0, 0.1),
                            rod, Vec3(0, 2, 0));

    matter.updGround().addBodyDecoration(Vec3(0,0,.1),
        DecorativeLine(Vec3(0), 3*UnitVec3(0,1,0))
        .setColor(Red));

    // Visualize the system, reporting an output frame every 1/30 of a simulated
    // second. The Visualizer's default frame rate is 30fps, and it will slow the
    // simulation down to keep to that speed, so we'll get exactly real time
    // this way. We don't want the default ground and sky background here.
    Visualizer viz(system);
    viz.setBackgroundType(Visualizer::SolidColor); // default is white
    system.addEventReporter(new Visualizer::Reporter(viz, 1./30));

    // Initialize the system and state.

    system.realizeTopology();
    State state = system.getDefaultState();

    // Simulate it.

    RungeKuttaMersonIntegrator integ(system);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(1000.0);
}
