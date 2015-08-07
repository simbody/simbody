/* -------------------------------------------------------------------------- *
 *             Simbody(tm) Example: SimplePlanarMechanism                     *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2015 Stanford University and the Authors.           *
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

// This very simple example builds a 3-body planar mechanism that does 
// nothing but rock back and forth for 10 seconds. Note that Simbody always
// works in 3D; this mechanism is planar because of the alignment of its 
// joints not because it uses any special 2D features. The mechanism looks
// like this:
//                              @
//                     @--------+--------@
//     Y               |        |         \
//     |               |        |          \
//     |               |        |           \
//     /-----X         *        *            *
//    /
//    Z
//
// It consists of a central T-shaped body pinned to ground, and
// two pendulum bodies pinned to either side of the T. The @'s above represent
// pin joints rotating about the Z axes. Each body's mass is concentrated into
// point masses shown by *'s above. Gravity is in the -Y direction.
//
// We'll add rigid joint stops to the two arms as a very simple test of
// unilateral contact constraints.

namespace {
const Real StopCoefRest = .1;
}

int main() {
  try { // catch errors if any

    // Create the system, with subsystems for the bodies and some forces.
    MultibodySystem system; 
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::Gravity gravity(forces, matter, -YAxis, 9.8);

    // Describe a body with a point mass at (0, -3, 0) and draw a sphere there.
    Real mass = 3; Vec3 pos(0,-3,0);
    Body::Rigid bodyInfo(MassProperties(mass, pos, UnitInertia::pointMassAt(pos)));
    bodyInfo.addDecoration(pos, DecorativeSphere(.2).setOpacity(.5));

    // Create the tree of mobilized bodies, reusing the above body description.
    MobilizedBody::Pin bodyT  (matter.Ground(), Vec3(0), bodyInfo, Vec3(0));
    MobilizedBody::Pin leftArm(bodyT, Vec3(-2, 0, 0),    bodyInfo, Vec3(0));
    MobilizedBody::Pin rtArm  (bodyT, Vec3(2, 0, 0),     bodyInfo, Vec3(0));

    matter.adoptUnilateralContact(
        new HardStopLower(leftArm, MobilizerQIndex(0), -0.1, StopCoefRest));
    matter.adoptUnilateralContact(
        new HardStopUpper(leftArm, MobilizerQIndex(0), .15, StopCoefRest));

    matter.adoptUnilateralContact(
        new HardStopLower(rtArm, MobilizerQIndex(0), -0.2, StopCoefRest));
    matter.adoptUnilateralContact(
        new HardStopUpper(rtArm, MobilizerQIndex(0), .25, StopCoefRest));

    // Ask for visualization every 1/30 second.
    system.setUseUniformBackground(true); // turn off floor 
    Visualizer viz(system);
    viz.setShowSimTime(true);
    system.adoptEventReporter(new Visualizer::Reporter(viz, 1./30));
    
    // Initialize the system and state.    
    State state = system.realizeTopology();
    bodyT.setRate(state, 2);

    const double SimTime = 10;

    // Simulate with acceleration-level time stepper.
    RungeKuttaMersonIntegrator integ(system);
    TimeStepper ts(integ);
    ts.initialize(state);
    ts.stepTo(SimTime);

    // Simulate with velocity-level time stepper.
    SemiExplicitEulerTimeStepper sxe(system);
    sxe.initialize(state);
    const Real h = .001;
    const int DrawEvery = 33; // draw every nth step ~= 33ms
    int nSteps = 0;
    while (true) {
        const State& sxeState = sxe.getState();
        const double t = sxe.getTime();
        const bool isLastStep = t + h/2 > SimTime;

        if((nSteps%DrawEvery)==0 || isLastStep) {
            viz.report(sxeState);
            printf("%7.4f %9.3g %9.3g\n", t,
                   leftArm.getAngle(sxeState),
                   rtArm.getAngle(sxeState));
        }

        if (isLastStep)
            break;

        sxe.stepTo(sxeState.getTime() + h);
        ++nSteps;
    }

  } catch (const std::exception& e) {
    std::cout << "ERROR: " << e.what() << std::endl;
    return 1;
  }

    return 0;
}
