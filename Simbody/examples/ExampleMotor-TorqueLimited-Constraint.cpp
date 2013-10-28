/* -------------------------------------------------------------------------- *
 *  Simbody(tm) Example: Torque-limited motor using speed control Constraint  *
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
#include <iostream>
#include <algorithm>
using std::cout; using std::endl;

using namespace SimTK;

/*
What this example demonstrates
------------------------------
1. Implementing a rate-controlled but torque-limited motor using a pair of 
   model elements: a speed-controlling Constraint element when operating 
   below the torque limit, and a constant-torque Force element when operating
   at the limit (meaning we won't be able to maintain desired speed).
2. Changing desired motor speed and torque limit externally from user input.
3. Performing a momentum balance analysis when enabling the constraint, so 
   that discontinuous speed changes are made consistent with Newton's laws.
4. Integrating with a small maximum step size, using just an Integrator
   without a TimeStepper to handle events.
5. Manually checking for events between integration steps, and handling 
   those events by making state changes before resuming integration.
6. Creating sliders, menus, and screen text in the Visualizer.

Compare this with ExampleMotor-TorqueLimited-Controller to see the same 
system implemented with simpler logic, using a Custom force element with a 
controller rather than a Constraint.

The advantages of using an intermittent constraint to implement this motor are: 
(1) fast execution time and (2) perfect speed tracking while the constraint is 
enabled. The primary disadvantage is coding the logic needed for making the 
switch between speed controlled and torque controlled operation and back. Here 
we are simplifying that by looking for transitions only at discrete times rather
than having the integrator isolate the event occurrence precisely. As a 
consequence, we have to deal with potentially large discrete state changes.

Strategy used here
------------------
Assume initially that we can control speed and use the constraint. Then prior
to each small integration step:
1. If speed control is active: Monitor the required torque; if it exceeds the
   limit then turn off the constraint and turn on the constant-torque element, 
   applying the torque in the same direction that the excessive constraint 
   torque would have been applied (opposing the impending speed error).
2. If torque control is active: Monitor the speed error. When the applied 
   torque is in the same direction as the error (meaning it is making the error
   worse), the speed must have gone from too slow to too fast or vice versa. Try
   matching the speed and enabling the constraint. If the required torque is 
   within range, enable speed control. Otherwise, correct the direction of the 
   applied torque to oppose the speed error and remain in torque control.
3. Periodically update the Visualizer display and poll for user input. A user
   change to the desired speed always changes us to torque control until the
   speed change has been achieved.

Alternative strategies
----------------------
1. Because we're checking for events only between steps, our choice of step size
here limits how accurately we can isolate the transition events between speed- 
and torque-control modes, requiring a small maximum step size for good event
isolation. An alternative strategy is to let the integrator take whatever steps 
it wants but use a TimeStepper and event witness functions to automatically
isolate events when they occur. That would allow much larger steps and better
event isolation.
2. A simpler implementation is to use only a force element to implement the 
motor, with a controller to generate the appropriate output torque to track a
desired speed, with a limit set. This method is illustrated in the companion
Simbody example ExampleMotor-TorqueLimited-Controller. This eliminates the need 
for any constraint-switching logic so makes for a simpler flow of control. 
However, depending on the control gains it may fail to track the desired speed 
well, or may make the problem stiff causing the integrator to require smaller 
time steps.
*/

namespace {     // Use anonymous namespace to keep global symbols private.
// Motor parameters.
const Real MaxMotorSpeed      = 10;  // rad/s
const Real MaxTorqueLimit     = 500; // N-m
const Real InitialMotorSpeed  = 1;
const Real InitialTorqueLimit = 100;

// Joint stop material parameters.
const Real StopStiffness = 10000; // stiffness in left joint stop
const Real StopDissipation = 0.5; // dissipation rate

// Integration step size, display update, user polling rates.
const Real MaxStepSize    = Real(1/240.); //  4 1/6 ms (240 Hz)
const int  DrawEveryN     = 8;            // 33 1/3 ms frame update (30 Hz)
const int  PollUserEveryN = 16;           // 66 2/3 ms user response lag (15 Hz)


//==============================================================================
//                               MY MECHANISM
//==============================================================================
// This class builds the multibody system and supports the various operations
// we'll need to turn the constraint on and off.
class MyMechanism {
public:
    MyMechanism(); // Uses default destructor.

    const MultibodySystem& getSystem() const {return m_system;}
    const State& getDefaultState()     const {return m_system.getDefaultState();}

    Real getDesiredSpeed(const State& state) const 
    {   return m_speedController.getSpeed(state); }
    Real getTorqueLimit(const State& state) const
    {   return std::abs(m_torqueController.getForce(state)); }
    Real getActualSpeed(const State& state) const 
    {   return m_rtArm.getOneU(state, MobilizerUIndex(0)); }
    Real getSpeedError(const State& state) const
    {   return getActualSpeed(state) - getDesiredSpeed(state); }

    // Returns the actual torque being applied, whether from the constraint or
    // the constant-torque element.
    Real getMotorTorque(const State& state) const {
        if (!isSpeedControlEnabled(state))
            return m_torqueController.getForce(state);
        m_system.realize(state, Stage::Acceleration);   // for multipliers
        return -m_speedController.getMultiplier(state); // watch sign convention
    }

    // Return true if speed control in effect; otherwise it's torque control.
    bool isSpeedControlEnabled(const State& state) const 
    {   return !m_speedController.isDisabled(state); }

    // Switch from speed control to torque control. This leaves velocity 
    // unchanged but causes a reduction in motor torque. It always succeeds.
    void switchToTorqueControl(State& state) const;

    // Turn on speed control if possible, meaning that the acceleration can be
    // made zero with torque below the limit. If the actual and desired 
    // speeds don't match exactly, then we must make an impulsive speed change
    // that needs to be momentum balanced.
    void tryToSwitchToSpeedControl(State& state) const;

    // Change the desired speed. This will turn speed control off if it is on,
    // until the speeds can be made to match using torque control.
    void changeDesiredSpeed(State& state, Real newSpeed) const;

    // Change the maximum torque limit. If we're already running torque limited
    // this will cause an immediate change to the torque being applied.
    void changeTorqueLimit(State& state, Real newTorqueLimit) const;

    // Update the Visualizer, including the sliders.
    void draw(const State& state) const;

    // Returns true when it is time to quit.
    bool processUserInput(State& state) const;

private:
    void constructSystem();
    void setUpVisualizer();

    MultibodySystem                 m_system;
    SimbodyMatterSubsystem          m_matter;
    GeneralForceSubsystem           m_forces;
    Visualizer                      m_viz;
    Visualizer::InputSilo*          m_userInput; // just a ref; not owned here

    MobilizedBody::Pin              m_bodyT, m_leftArm, m_rtArm;
    Constraint::ConstantSpeed       m_speedController;
    Force::MobilityConstantForce    m_torqueController;
};

// Write interesting integrator info to stdout.
void dumpIntegratorStats(const Integrator& integ);
}


//==============================================================================
//                                  MAIN
//==============================================================================
// Simulate forever with a small max step size.
int main() {
    try { // catch errors if any

    // Create the Simbody MultibodySystem and Visualizer.
    MyMechanism mech;

    // We're forcing very small step sizes so should use very low order
    // integration and loose accuracy. We must prevent interpolation so that the
    // state returned after each step is the integrator's "advanced state", 
    // which is modifiable, in case we need to make a state change.
    SemiExplicitEuler2Integrator integ(mech.getSystem());
    integ.setAccuracy(Real(1e-1)); // 10%
    integ.setAllowInterpolation(false);

    integ.initialize(mech.getDefaultState());
    unsigned stepNum = 0;
    while (true) {
        // Get access to State being advanced by the integrator. Interpolation 
        // must be off so that we're modifying the actual trajectory.
        State& state = integ.updAdvancedState();

        // Output a frame to the Visualizer if it is time.
        if (stepNum % DrawEveryN == 0)
            mech.draw(state);

        // Check for user input periodically. 
        if (stepNum % PollUserEveryN == 0 && mech.processUserInput(state))
            break; // stop if user picked Quit from the Run menu

        // Check for speed/torque control transition events and handle.
        const Real trqNow = mech.getMotorTorque(state);
        if (mech.isSpeedControlEnabled(state)) {
           if (std::abs(trqNow) > mech.getTorqueLimit(state)) {
                printf("%d: SWITCH TO TORQUE CONTROL cuz trqNow=%g\n", 
                       stepNum, trqNow);
                mech.switchToTorqueControl(state);
           } 
        } else { // Currently limiting torque.
            // If the torque is now in the same direction as the error, try
            // to switch back to speed control. If that doesn't work we'll at
            // least reverse the applied maximum torque.
            const Real errNow = mech.getSpeedError(state);
            if (errNow * trqNow > 0) {
                printf("%d: TRY SPEED CONTROL cuz errNow=%g, trqNow=%g\n",
                       stepNum, errNow, trqNow);
                mech.tryToSwitchToSpeedControl(state);
            }
        }

        // Advance time by MaxStepSize. Might take multiple internal steps to 
        // get there, depending on required accuracy.
        integ.stepBy(MaxStepSize);
        ++stepNum;
    }

    // DONE.
    dumpIntegratorStats(integ);

    } catch (const std::exception& e) {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}



//==============================================================================
//                        MY MECHANISM IMPLEMENTATION
//==============================================================================

//------------------------------- CONSTRUCTOR ----------------------------------
MyMechanism::MyMechanism() 
:   m_system(), m_matter(m_system), m_forces(m_system), m_viz(m_system), 
    m_userInput(0)
{
    constructSystem();
    setUpVisualizer();
}

//----------------------------- CONSTRUCT SYSTEM -------------------------------
void MyMechanism::constructSystem() {
    Force::Gravity(m_forces, m_matter, -YAxis, Real(9.80665));

    // Describe a body with a point mass at (0, -3, 0) and draw a sphere there.
    Real mass = 3; Vec3 pos(0,-3,0);
    Body::Rigid bodyInfo(MassProperties(mass, pos, UnitInertia::pointMassAt(pos)));
    bodyInfo.addDecoration(pos, DecorativeSphere(Real(.2)).setOpacity(.5));

    // Create the tree of mobilized bodies, reusing the above body description.
    m_bodyT   = MobilizedBody::Pin(m_matter.Ground(),Vec3(0), bodyInfo,Vec3(0));
    m_leftArm = MobilizedBody::Pin(m_bodyT,Vec3(-2,0,0), bodyInfo,Vec3(0));
    m_rtArm   = MobilizedBody::Pin(m_bodyT,Vec3(2,0,0),  bodyInfo,Vec3(0,-1,0));

    // Add some damping.
    Force::MobilityLinearDamper(m_forces, m_bodyT,   MobilizerUIndex(0), 10);
    Force::MobilityLinearDamper(m_forces, m_leftArm, MobilizerUIndex(0), 30);
    Force::MobilityLinearDamper(m_forces, m_rtArm,   MobilizerUIndex(0), 10);

    // Add a joint stop to the left arm restricting it to q in [0,Pi/5].
    Force::MobilityLinearStop(m_forces, m_leftArm, MobilizerQIndex(0), 
        StopStiffness, StopDissipation,
        -Pi/8,   // lower stop
         Pi/8);  // upper stop

    // Use built-in ConstantSpeed constraint as a low-budget motor model.
    m_speedController = Constraint::ConstantSpeed(m_rtArm, InitialMotorSpeed);

    // This is used when we're at the maximum torque.
    m_torqueController = Force::MobilityConstantForce(m_forces, m_rtArm, 
                                                      InitialTorqueLimit);
    m_torqueController.setDisabledByDefault(true);

    // We're done with the System; finalize it.
    m_system.realizeTopology();
}

//------------------------- SWITCH TO TORQUE CONTROL ---------------------------
// The switch to torque control causes a sudden drop in torque but that doesn't
// require any sudden change to velocity.
void MyMechanism::switchToTorqueControl(State& state) const {
    assert(isSpeedControlEnabled(state));

    const Real oldTrq = getMotorTorque(state);
    const Real newTrq = sign(oldTrq)*getTorqueLimit(state);
    printf("  switchToTorqueControl(): change torque from %g -> %g\n", 
        oldTrq, newTrq);

    m_speedController.disable(state);
    m_torqueController.enable(state);
    m_torqueController.setForce(state, newTrq);
    m_system.realize(state, Stage::Velocity);
}

//---------------------- TRY TO SWITCH TO SPEED CONTROL ------------------------
// We want to switch to the constraint to enforce the desired speed.
// If the actual speed and desired speed don't match, this will first require
// an impulsive change to the speed. We'll perform a momentum balance analysis
// here to calculate system wide velocity changes that will satisfy f=ma.
//   Solve GM\~G lambda = deltaV to calculate constraint impulse lambda.
//   Then f = ~G lambda is the corresponding generalized impulse.
//   And deltaU = M\f is the resulting change to the generalized speeds.
//
// But after all that, we might find that the torque required is more than
// we allow. In that case we have to put the speed back where we found it
// and switch back to torque control, but with the torque applied in the same
// direction as the excessive constraint torque would have been.
void MyMechanism::tryToSwitchToSpeedControl(State& state) const {
    assert(!isSpeedControlEnabled(state));

    const Real curSpeed = getActualSpeed(state);
    const Real desSpeed = getDesiredSpeed(state);
    printf("  tryToSwitchToSpeedControl(): torque now is %g, speed err=%g\n", 
        getMotorTorque(state), curSpeed-desSpeed);


    // Save the current velocities in case we have to put them back.
    const Vector prevU = state.getU();
    m_torqueController.disable(state);
    m_speedController.enable(state); // Tentatively enable the constraint.

    // Dynamics operators require this stage.
    m_system.realize(state, Stage::Velocity);

    // Momentum balance analysis: see comment above.
    Vector deltaV(1, desSpeed-curSpeed);
    Vector allDeltaV, lambda, f, deltaU;
    m_speedController.setMyPartInConstraintSpaceVector(state, deltaV, allDeltaV);
    m_matter.solveForConstraintImpulses(state, allDeltaV, lambda);
    m_matter.multiplyByGTranspose(state, lambda, f);
    m_matter.multiplyByMInv(state,f,deltaU);
    state.updU() += deltaU;
    // This final projection shouldn't be necessary, but just in case ...
    m_system.projectU(state); // leaves state at Stage::Velocity

    // Now that the speed constraint is satisfied, try keeping the acceleration
    // zero using the constraint and check how much torque that requires.
    m_system.realize(state, Stage::Acceleration);
    const Real requiredTorque = getMotorTorque(state);
    if (std::abs(requiredTorque) > getTorqueLimit(state)) {
        m_speedController.disable(state);
        m_torqueController.enable(state);
        const Real newMaxTrq = sign(requiredTorque)*getTorqueLimit(state);
        m_torqueController.setForce(state, newMaxTrq);
        printf("  ... NO switch to speed control cuz torque would be %g;"
               " applying %g instead.\n",
            requiredTorque, newMaxTrq);
        state.updU() = prevU; // restore previous speed
        m_system.realize(state, Stage::Velocity);
        return;
    }

    printf("  ... switched to speed control with trq=%g and impulsive speed change:"
        " %g -> %g\n",  requiredTorque, curSpeed, desSpeed);
    cout << "  ... momentum balance deltaU was" << deltaU << endl;
}

//--------------------------- CHANGE DESIRED SPEED -----------------------------
// The user has specified a new desired motor speed. This always requires a
// switch to torque control if we're in speed control now, because an 
// instantaneous speed change would require an infinite torque. The torque 
// direction depends on whether the new desired speed is larger or smaller than 
// the current actual speed.
void MyMechanism::changeDesiredSpeed(State& state, Real newDesiredSpeed) const {
    const Real actualSpeed = getActualSpeed(state);
    const Real changeDirection = (Real)sign(newDesiredSpeed - actualSpeed);
    if (changeDirection == 0)
        return; // nothing to do

    printf("Desired speed change: %g -> %g\n", actualSpeed, newDesiredSpeed);
    m_speedController.setSpeed(state, newDesiredSpeed);
    if (isSpeedControlEnabled(state)) {
        printf("...requires switch to torque control\n");
        m_speedController.disable(state);
        m_torqueController.enable(state);
    }

    // Torque control is on; set torque to move in direction of new speed.
    const Real newTrq = changeDirection*getTorqueLimit(state);
    m_torqueController.setForce(state, newTrq);
    printf("...applying max torque %g to change speed\n", newTrq);

    m_system.realize(state, Stage::Velocity);
}

//--------------------------- CHANGE TORQUE LIMIT ------------------------------
// The user has specified a new value for the magnitude of the maximum torque.
// The sign will be inherited from the sign of the torque currently in use. No
// change to the control state is done here; that will occur in the main loop
// if needed.
void MyMechanism::changeTorqueLimit(State& state, Real newTorqueLimit) const {
    const Real oldTorqueLimit = getTorqueLimit(state);
    printf("Torque limit change: %g -> %g\n", oldTorqueLimit, newTorqueLimit);

    // Put torque in same direction as current torque (might not be using the
    // torque controller though).
    const Real trq = getMotorTorque(state);
    Real dir = (Real)sign(trq); if (dir==0) dir=1;
    m_torqueController.setForce(state, dir*newTorqueLimit);
}

//---------------------------- PROCESS USER INPUT ------------------------------
namespace {
// Constants for the user interaction widgets.
// Ids for the sliders.
const int SliderIdMotorSpeed = 1, SliderIdTorqueLimit = 2, 
          SliderIdTach = 3, SliderIdTorque = 4; // these two are used for output
// Ids for things on the Run Menu.
const int MenuIdRun = 1;
const int ResetItem=1, QuitItem=2;
}

// Return true if it is time to quit.
bool MyMechanism::processUserInput(State& state) const {
    int whichSlider, whichMenu, whichItem; Real newValue;

    // Did a slider move?
    if (m_userInput->takeSliderMove(whichSlider, newValue)) {
        switch(whichSlider) {
        case SliderIdMotorSpeed:
            // This will momentum balance if necessary.
            changeDesiredSpeed(state, newValue);
            break;
        case SliderIdTorqueLimit:
            changeTorqueLimit(state, newValue);
            m_viz.setSliderRange(SliderIdTorque, -newValue, newValue); 
            break;
        }
    }

    // Was there a menu pick?
    if (m_userInput->takeMenuPick(whichMenu, whichItem)) {
        if (whichItem == QuitItem) 
            return true; // done

        // If Reset, stop the motor and zero out all the q's and u's. 
        // Tell visualizer to update the sliders to match.
        if (whichItem == ResetItem) {
            // Don't momentum balance here!
            m_speedController.setSpeed(state, 0);
            m_viz.setSliderValue(SliderIdMotorSpeed, 0)
                 .setSliderValue(SliderIdTach, 0);

            m_torqueController.setForce(state, InitialTorqueLimit);
            m_viz.setSliderValue(SliderIdTorqueLimit, InitialTorqueLimit)
                 .setSliderRange(SliderIdTorque, -InitialTorqueLimit, 
                                                  InitialTorqueLimit) 
                 .setSliderValue(SliderIdTorque, 0);

            state.updQ() = 0; // all positions to zero
            state.updU() = 0; // all velocities to zero
            m_system.projectQ(state);
            m_system.projectU(state);
        }
    }
    
    return false; // keep going
}

//----------------------------- CLASS SHOW STUFF -------------------------------
// We'll provide the Visualizer with an object of this class to generate any
// per-frame text and geometry we want to show on screen.
namespace {
class ShowStuff : public DecorationGenerator {
public:
    ShowStuff(const MyMechanism& mech) : m_mech(mech) {}

    void generateDecorations(const State&                state, 
                             Array_<DecorativeGeometry>& geometry) OVERRIDE_11
    {
        DecorativeText msg;
        msg.setIsScreenText(true);
        if (m_mech.isSpeedControlEnabled(state))
            msg.setText("OK");
        else
            msg.setText("TORQUE LIMITED err=" 
                         + String(m_mech.getSpeedError(state), "%.2g"));
        geometry.push_back(msg);
    }
private:
    const MyMechanism& m_mech;
};
}

//----------------------------- SET UP VISUALIZER ------------------------------
void MyMechanism::setUpVisualizer() {
    m_viz.setShutdownWhenDestructed(true) // make sure display window dies
         .setBackgroundType(Visualizer::SolidColor); // turn off Ground & Sky

    // Add sliders.
    m_viz.addSlider("Motor speed", SliderIdMotorSpeed, 
                    -MaxMotorSpeed, MaxMotorSpeed, InitialMotorSpeed)
         .addSlider("Torque limit", SliderIdTorqueLimit, 
                    0, MaxTorqueLimit, InitialTorqueLimit)
         .addSlider("Tach",   SliderIdTach,   
                    -MaxMotorSpeed,  MaxMotorSpeed,  0)
         .addSlider("Torque", SliderIdTorque, 
                    -InitialTorqueLimit, InitialTorqueLimit, 0);

    // Add Run menu.
    Array_<std::pair<String,int> > runMenuItems;
    runMenuItems.push_back(std::make_pair("Reset", ResetItem));
    runMenuItems.push_back(std::make_pair("Quit", QuitItem));
    m_viz.addMenu("Run", MenuIdRun, runMenuItems);

    // Add per-frame text and geometry.
    m_viz.addDecorationGenerator(new ShowStuff(*this));

    // Add an input listener so the user can talk to us.
    m_userInput = new Visualizer::InputSilo();
    m_viz.addInputListener(m_userInput);
}

//------------------------------------ DRAW ------------------------------------
// Update the Visualizer, including the output sliders.
void MyMechanism::draw(const State& state) const {
    m_viz.report(state);
    m_viz.setSliderValue(SliderIdTach, getActualSpeed(state))
         .setSliderValue(SliderIdTorque, getMotorTorque(state));
}



//==============================================================================
//                        DUMP INTEGRATOR STATS
//==============================================================================
namespace {
void dumpIntegratorStats(const Integrator& integ) {
    const int evals = integ.getNumRealizations();
    std::cout << "\nDone -- simulated " << integ.getTime() << "s with " 
            << integ.getNumStepsTaken() << " steps, avg step=" 
        << (1000*integ.getTime())/integ.getNumStepsTaken() << "ms " 
        << (1000*integ.getTime())/evals << "ms/eval\n";

    printf("Used Integrator %s at accuracy %g:\n", 
        integ.getMethodName(), integ.getAccuracyInUse());
    printf("# STEPS/ATTEMPTS = %d/%d\n",  integ.getNumStepsTaken(), 
                                          integ.getNumStepsAttempted());
    printf("# ERR TEST FAILS = %d\n",     integ.getNumErrorTestFailures());
    printf("# REALIZE/PROJECT = %d/%d\n", integ.getNumRealizations(), 
                                          integ.getNumProjections());
}
}

