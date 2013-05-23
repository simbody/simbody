/* -------------------------------------------------------------------------- *
 *             Simbody(tm) Example: MotorWithSpeedControl                     *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
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

// This is an example that shows several things:
// 1. Integrating with a small maximum step size, using just an integrator
//    without a time stepper.
// 2. Implementing a fixed-rate motor using a Motion (prescribed motion).
// 3. Changing motor speed externally from user input.
// 4. Creating sliders and menus in the Visualizer.
//
// Note: we are allowing the user to violently change motor speed, but we are
// not performing a momentum balance so the change is non-physical modeled
// this way. A momentum balance would allow us to calculate an impulse to be
// applied to the whole mechanism that would make the change momentum conserving
// as though the change had been implemented with a very large force over a very
// short duration.

const Real InitialMotorRate   = 1; // rad/s
const Real InitialDissipation = 0.5; // damping in left joint stop

// Ids for the two sliders.
const int SliderIdMotorSpeed = 1, SliderIdDissipation = 2;

// Ids for things on the Run Menu.
const int MenuIdRun = 1;
static const int ResetItem=1, QuitItem=2;

const Real StopStiffness = 10000; // not changeable here
int main() {
    try { // catch errors if any

    // Create the system, with subsystems for the bodies and some forces.
    MultibodySystem system; 
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::Gravity gravity(forces, matter, -YAxis, 9.8);
    system.setUseUniformBackground(true); // request no ground & sky    

    // Describe a body with a point mass at (0, -3, 0) and draw a sphere there.
    Real mass = 3; Vec3 pos(0,-3,0);
    Body::Rigid bodyInfo(MassProperties(mass, pos, UnitInertia::pointMassAt(pos)));
    bodyInfo.addDecoration(pos, DecorativeSphere(.2).setOpacity(.5));

    // Create the tree of mobilized bodies, reusing the above body description.
    MobilizedBody::Pin bodyT  (matter.Ground(), Vec3(0), bodyInfo, Vec3(0));
    MobilizedBody::Pin leftArm(bodyT, Vec3(-2, 0, 0),    bodyInfo, Vec3(0));
    MobilizedBody::Pin rtArm  (bodyT, Vec3(2, 0, 0),     bodyInfo, Vec3(0,-1,0));

    // Add some damping.
    Force::MobilityLinearDamper damper(forces, leftArm, MobilizerUIndex(0), 0.1);

    // Use built-in Steady Motion as a low-budget motor model.
    Motion::Steady motor(rtArm, InitialMotorRate);

    // Add a joint stop to the left arm restricting it to q in [0,Pi/5].
    Force::MobilityLinearStop stop(forces, leftArm, MobilizerQIndex(0), 
        StopStiffness,
        InitialDissipation,
        -Pi/8,   // lower stop
         Pi/8);  // upper stop

    Visualizer viz(system);

    // Add sliders.
    viz.addSlider("Motor speed", SliderIdMotorSpeed, -10, 10, InitialMotorRate);
    viz.addSlider("Dissipation", SliderIdDissipation, 0, 10, InitialDissipation);

    // Add Run menu.
    Array_<std::pair<String,int> > runMenuItems;
    runMenuItems.push_back(std::make_pair("Reset", ResetItem));
    runMenuItems.push_back(std::make_pair("Quit", QuitItem));
    viz.addMenu("Run", MenuIdRun, runMenuItems);

    Visualizer::InputSilo* userInput = new Visualizer::InputSilo();
    viz.addInputListener(userInput);
    
    // Initialize the system and state.    
    State initState = system.realizeTopology();

    // Simulate forever with a small max step size. Check for user input
    // in between steps. Note: an alternate way to do this is to let the
    // integrator take whatever steps it wants but use a TimeStepper to 
    // manage a periodic event handler to poll for user input. Here we're 
    // treating completion of a step as an event.
    const Real MaxStepSize = 0.01; // 10ms
    const int  DrawEveryN = 3;     // 3 steps = 30ms
    RungeKuttaMersonIntegrator integ(system);
    integ.initialize(initState);
    // Don't permit interpolation because we want the state returned after
    // a step to be modifiable.
    integ.setAllowInterpolation(false);

    int stepsSinceViz = DrawEveryN-1;
    while (true) {
        if (++stepsSinceViz % DrawEveryN == 0) {
            viz.report(integ.getState());
            stepsSinceViz = 0;
        }

        // Advance time by MaxStepSize (might take multiple steps to get there).
        integ.stepBy(MaxStepSize);

        // Now poll for user input.
        int whichSlider, whichMenu, whichItem; Real newValue;

        // Did a slider move?
        if (userInput->takeSliderMove(whichSlider, newValue)) {
            State& state = integ.updAdvancedState();
            switch(whichSlider) {
            case SliderIdMotorSpeed:
                // TODO: momentum balance?
                motor.setRate(state, newValue);
                system.realize(state, Stage::Position);
                system.prescribeU(state);
                break;
            case SliderIdDissipation:
                stop.setMaterialProperties(state, StopStiffness, newValue);
                system.realize(state, Stage::Position);
                break;
            }
        }

        // Was there a menu pick?
        if (userInput->takeMenuPick(whichMenu, whichItem)) {
            if (whichItem == QuitItem) 
                break; // done

            // If Reset, stop the motor and restore default dissipation. 
            // Tell visualizer to update the sliders to match.
            // Zero out all the q's and u's.
            if (whichItem == ResetItem) {
                State& state = integ.updAdvancedState();
                motor.setRate(state, 0);
                viz.setSliderValue(SliderIdMotorSpeed, 0);

                stop.setMaterialProperties(state, StopStiffness, InitialDissipation);
                viz.setSliderValue(SliderIdDissipation, InitialDissipation);

                state.updQ() = 0; // all positions to zero
                state.updU() = 0; // all velocities to zero
                system.realize(state, Stage::Position);
                system.prescribeU(state);
            }
        }
    }

    } catch (const std::exception& e) {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
