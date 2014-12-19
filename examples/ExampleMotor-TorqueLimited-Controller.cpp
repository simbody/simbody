/* -------------------------------------------------------------------------- *
 *       Simbody(tm) Example: Torque-limited motor using PI controller        *
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
1. Implementing a rate-controlled but torque-limited motor using a Custom force 
   element that uses PI (Proportional-Integral control) to track speed when
   possible without exceeding the torque limit.
2. Changing desired motor speed and torque limit externally from user input.
3. Integrating using just an Integrator without a TimeStepper with a modest 
   maximum step size that allows for output of Visualizer frames and polling for
   user input between steps.
4. Creating sliders, menus, and screen text in the Visualizer.

Compare this with ExampleMotor-TorqueLimited-Constraint to see the same 
system implemented more efficiently, but with more complex logic, using
an intermittent constraint.

Strategy used here
------------------
Provide a Custom force element that implements the torque-limited motor. It 
allocates a state variable used for integral control so that a steady-state
torque can be generated. There is also proportional control that generates a
torque proportional to the tracking error. These are summed and then capped 
by the user-controlled torque limit. Two gains must be selected, with some 
tradeoff of performance for tracking effectiveness.

The integration consists of repeated stepping by the maximum step size, with
periodic output of Visualizer frames and polling for user input that may 
result in changes to the desired speed and torque limit. Nothing else need be
managed in the integration loop since the motor force element is always on.

Alternative strategies
----------------------
1. Because we are generating screen updates and polling for user input only 
between steps, we have to restrict the maximum step size accordingly. An 
alternative strategy is to let the integrator take whatever steps it wants but 
use a TimeStepper with a periodic event reporter to generate interpolated
display updates when needed, and a periodic event handler to poll for user input 
and make state changes as a result. That would allow larger steps, depending on 
how much user response lag is tolerable (usually 100ms is acceptable).
2. A more difficult approach is to use an intermittent constraint to control
speed rather than doing it with a speed-tracking force element as is done here.
That has the advantage of perfect speed tracking when the torque is below the
limit, requires no arbitrary control gains, and can provide the best performance
since a constraint-driven simulation is non-stiff. See the companion Simbody
example ExampleMotor-TorqueLimited-Constraint for an implementation of this
alternative.
*/

namespace {     // Use anonymous namespace to keep global symbols private.
// Motor parameters.
const Real MaxMotorSpeed      = 10;  // rad/s
const Real MaxTorqueLimit     = 500; // N-m
const Real InitialMotorSpeed  = 1;
const Real InitialTorqueLimit = 100;

// Joint stop material parameters.
const Real StopStiffness   = 10000; // stiffness in left joint stop
const Real StopDissipation = 0.5;   // dissipation rate 

// Integration step size, display update, user polling rates.
const Real MaxStepSize    = Real(1/30.); // 33 1/3 ms (30 Hz)
const int  DrawEveryN     = 1;           // 33 1/3 ms frame update (30 Hz)
const int  PollUserEveryN = 2;           // 66 2/3 ms user response lag (15 Hz)

// Gains for the PI controller.
const Real ProportionalGain = 1000;
const Real IntegralGain     = 50000;

//==============================================================================
//                         MY TORQUE LIMITED MOTOR
//==============================================================================
// This Force element implements a torque-limited motor that runs at a user-
// selected speed unless that would require too much torque.
class MyTorqueLimitedMotor : public Force::Custom::Implementation {
public:
    MyTorqueLimitedMotor(const MobilizedBody& mobod, MobilizerUIndex whichU)
    :   m_matter(mobod.getMatterSubsystem()), m_mobod(mobod), m_whichU(whichU)
    {}

    void setDesiredSpeed(State& state, Real speed) const {
        Real& u_desired = Value<Real>::updDowncast(
            m_matter.updDiscreteVariable(state, m_desiredUIx));
        u_desired = speed;
    }

    Real getDesiredSpeed(const State& state) const {
        const Real u_desired = Value<Real>::downcast(
            m_matter.getDiscreteVariable(state, m_desiredUIx));
        return u_desired;
    }

    void setTorqueLimit(State& state, Real speed) const {
        Real& torqueLimit = Value<Real>::updDowncast(
            m_matter.updDiscreteVariable(state, m_torqueLimitIx));
        torqueLimit = speed;
    }

    Real getTorqueLimit(const State& state) const {
        const Real torqueLimit = Value<Real>::downcast(
            m_matter.getDiscreteVariable(state, m_torqueLimitIx));
        return torqueLimit;
    }

    Real getActualSpeed(const State& state) const {
        const Real u_actual = m_mobod.getOneU(state, m_whichU);
        return u_actual;
    }

    Real getSpeedError(const State& state) const {
        return getActualSpeed(state) - getDesiredSpeed(state);
    }

    // Combine integrated torque and proportional torque, and clamp to limit.
    Real getTorque(const State& state) const {
        const Real integTrq  = m_matter.getZ(state)[m_trqIx];
        const Real trqLimit  = getTorqueLimit(state);
        const Real proportionalTrq = -ProportionalGain*getSpeedError(state);
        const Real totalTrq = clamp(-trqLimit, integTrq+proportionalTrq, 
                                     trqLimit);
        return totalTrq;
    }

    // Set the integrated torque back to zero.
    void resetIntegratedTorque(State& state) const
    {   m_matter.updZ(state)[m_trqIx] = 0; }


    // Force::Custom::Implementation virtuals:

    void calcForce(const State&         state, 
                   Vector_<SpatialVec>& bodyForces, 
                   Vector_<Vec3>&       particleForces, 
                   Vector&              mobilityForces) const override
    {
        m_mobod.applyOneMobilityForce(state, m_whichU, getTorque(state), 
                                      mobilityForces);
    }

    Real calcPotentialEnergy(const State& state) const override {
        return 0;
    }

    // Allocate state variables: 
    //   Discrete variables for desired speed and torque limit;
    //   a continuous variable for the integral torque (a "z").
    void realizeTopology(State& state) const override {
        m_desiredUIx = m_matter.allocateDiscreteVariable
           (state, Stage::Acceleration, new Value<Real>(InitialMotorSpeed));
        m_torqueLimitIx = m_matter.allocateDiscreteVariable
           (state, Stage::Acceleration, new Value<Real>(InitialTorqueLimit));
        m_trqIx = m_matter.allocateZ(state, Vector(1,Real(0)));
    }

    // Calculate the derivative for the integral control state.
    void realizeAcceleration(const State& state) const override {
        const Real integTrq = m_matter.getZ(state)[m_trqIx];
        const Real trqLimit  = getTorqueLimit(state);
        const Real abstrq=std::abs(integTrq), trqsign = (Real)sign(integTrq);

        Real trqDot = -IntegralGain*getSpeedError(state);
        // Don't ask for more torque if we're already past the limit
        if (abstrq >= trqLimit && sign(trqDot)*trqsign >= 0)
            trqDot = 0;

        m_matter.updZDot(state)[m_trqIx] = trqDot; 
    }

private:
    const SimbodyMatterSubsystem&   m_matter;
    MobilizedBody                   m_mobod;
    MobilizerUIndex                 m_whichU;

    // Topology cache; written only once by realizeTopology().
    mutable DiscreteVariableIndex   m_desiredUIx;
    mutable DiscreteVariableIndex   m_torqueLimitIx;
    mutable ZIndex                  m_trqIx;
};



//==============================================================================
//                               MY MECHANISM
//==============================================================================
// This class builds the multibody system and supports the various operations
// we'll need to turn the constraint on and off.
class MyMechanism {
public:
    MyMechanism(); // Uses default destructor.

    const MultibodySystem& getSystem() const {return m_system;}
    const State& getDefaultState() const {return m_system.getDefaultState();}

    Real getDesiredSpeed(const State& state) const 
    {   return m_controller->getDesiredSpeed(state); }
    Real getTorqueLimit(const State& state) const
    {   return m_controller->getTorqueLimit(state); }
    Real getActualSpeed(const State& state) const 
    {   return m_controller->getActualSpeed(state); }
    Real getSpeedError(const State& state) const
    {   return m_controller->getSpeedError(state); }

    // Returns the torque currently being applied.
    Real getMotorTorque(const State& state) const
    {   return m_controller->getTorque(state); }

    // Change the desired speed.
    void changeDesiredSpeed(State& state, Real newSpeed) const;

    // Change the maximum torque limit.
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
    Visualizer::InputSilo*          m_userInput;  // just a ref; not owned here

    MobilizedBody::Pin              m_bodyT, m_leftArm, m_rtArm;
    MyTorqueLimitedMotor*           m_controller; // just a ref; not owned here
};

// Write interesting integrator info to stdout.
void dumpIntegratorStats(const Integrator& integ);
}


//==============================================================================
//                                  MAIN
//==============================================================================
// Simulate forever with a maximum step size small enough to provide adequate
// polling for user input and screen updating.
int main() {
    try { // catch errors if any

    // Create the Simbody MultibodySystem and Visualizer.
    MyMechanism mech;

    // We want real time performance here and don't need high accuracy, so we'll
    // use semi explicit Euler. However, we hope to take fairly large steps so
    // need error control to avoid instability. We must also prevent 
    // interpolation so that the state returned after each step is the 
    // integrator's "advanced state", which is modifiable, so that we can update
    // it in response to user input.
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

        // Advance time by MaxStepSize. Might take multiple internal steps to 
        // get there, depending on difficulty and required accuracy.
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
    m_userInput(0), m_controller(0)
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
    m_leftArm = MobilizedBody::Pin(m_bodyT,Vec3(-2, 0, 0), bodyInfo,Vec3(0));
    m_rtArm   = MobilizedBody::Pin(m_bodyT,Vec3(2, 0, 0),  bodyInfo,Vec3(0,-1,0));

    // Add some damping.
    Force::MobilityLinearDamper(m_forces, m_bodyT, MobilizerUIndex(0), 10);
    Force::MobilityLinearDamper(m_forces, m_leftArm, MobilizerUIndex(0), 30);
    Force::MobilityLinearDamper(m_forces, m_rtArm, MobilizerUIndex(0), 10);

    // Add a joint stop to the left arm restricting it to q in [0,Pi/5].
    Force::MobilityLinearStop(m_forces, m_leftArm, MobilizerQIndex(0), 
        StopStiffness,
        StopDissipation,
        -Pi/8,   // lower stop
         Pi/8);  // upper stop

    // Use custom controller to track speed up to torque limit.
    m_controller = new MyTorqueLimitedMotor(m_rtArm, MobilizerUIndex(0));
    Force::Custom(m_forces, m_controller); // takes over ownership

    // We're done with the System; finalize it.
    m_system.realizeTopology();
}

//--------------------------- CHANGE DESIRED SPEED -----------------------------
// The user has specified a new desired motor speed.
void MyMechanism::changeDesiredSpeed(State& state, Real newDesiredSpeed) const {
    const Real actualSpeed = getActualSpeed(state);
    printf("Desired speed change: %g -> %g\n", actualSpeed, newDesiredSpeed);
    m_controller->setDesiredSpeed(state, newDesiredSpeed);
}


//--------------------------- CHANGE TORQUE LIMIT ------------------------------
// The user has specified a new value for the magnitude of the maximum torque.
// The direction will oppose the speed error; we just need magnitude here.
void MyMechanism::changeTorqueLimit(State& state, Real newTorqueLimit) const {
    const Real oldTorqueLimit = getTorqueLimit(state);
    printf("Torque limit change: %g -> %g\n", oldTorqueLimit, newTorqueLimit);
    m_controller->setTorqueLimit(state, newTorqueLimit);
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
            m_controller->resetIntegratedTorque(state);
            m_controller->setDesiredSpeed(state, 0);
            m_viz.setSliderValue(SliderIdMotorSpeed, 0)
                 .setSliderValue(SliderIdTach, 0);

            m_controller->setTorqueLimit(state, InitialTorqueLimit);
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
                             Array_<DecorativeGeometry>& geometry) override
    {
        DecorativeText msg;
        msg.setIsScreenText(true);
        if (std::abs(m_mech.getMotorTorque(state)) < m_mech.getTorqueLimit(state))
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
// Update the Visualizer, including the sliders.
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
