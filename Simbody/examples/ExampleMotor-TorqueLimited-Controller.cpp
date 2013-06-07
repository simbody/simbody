/* -------------------------------------------------------------------------- *
 *        Simbody(tm) Example: Torque Limited Motor With Speed Controller     *
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

// This is an example that shows several things:
// 1. Implementing a rate-controlled but torque-limited motor using
//    a Custom force element that uses PI (Proportional-Integral control).
// 2. Changing desired motor speed and torque limit externally from user input.
// 3. Integrating with a small maximum step size, using just an Integrator
//    without a TimeStepper to handle events.
// 4. Creating sliders, menus, and screen text in the Visualizer.
//
// Compare this with ExampleMotor-TorqueLimited-Constraint to see the same 
// system implemented more efficiently, but with more complex logic, using
// an intermittent constraint.

// Use anonymous namespace to keep global symbols private.
namespace {

// Motor parameters.
const Real MaxMotorSpeed      = 10; // rad/s
const Real MaxTorqueLimit     = 500; // N-m
const Real InitialMotorSpeed  = 1;
const Real InitialTorqueLimit = 100;

// Joint stop material parameters.
const Real StopStiffness = 10000; // stiffness in left joint stop
const Real StopDissipation = 0.5; // dissipation rate 

// Constants for the user interaction widgets.
// Ids for the sliders.
const int SliderIdMotorSpeed = 1, SliderIdTorqueLimit = 2, 
          SliderIdTach = 3, SliderIdTorque = 4; // these two are used for output
// Ids for things on the Run Menu.
const int MenuIdRun = 1;
const int ResetItem=1, QuitItem=2;

//==============================================================================
//                         MY TORQUE LIMITED MOTOR
//==============================================================================
// This Force element implements a torque-limited motor that runs at a user-
// selected speed unless that would require too much torque.

const Real IntegralGain = 50000, ProportionalGain = 1000;
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


    // Force::Custom::Implementation virtuals:

    void calcForce(const State&         state, 
                   Vector_<SpatialVec>& bodyForces, 
                   Vector_<Vec3>&       particleForces, 
                   Vector&              mobilityForces) const OVERRIDE_11
    {
        m_mobod.applyOneMobilityForce(state, m_whichU, getTorque(state), 
                                      mobilityForces);
    }

    Real calcPotentialEnergy(const State& state) const OVERRIDE_11 {
        return 0;
    }

    // Allocate state variables: 
    //   Discrete variables for desired speed and torque limit;
    //   a continuous variable for the integral torque.
    void realizeTopology(State& state) const OVERRIDE_11 {
        m_desiredUIx = m_matter.allocateDiscreteVariable
           (state, Stage::Acceleration, new Value<Real>(InitialMotorSpeed));
        m_torqueLimitIx = m_matter.allocateDiscreteVariable
           (state, Stage::Acceleration, new Value<Real>(InitialTorqueLimit));
        m_trqIx = m_matter.allocateZ(state, Vector(1,Real(0)));
    }

    void realizeAcceleration(const State& state) const OVERRIDE_11 {
        const Real integTrq = m_matter.getZ(state)[m_trqIx];
        const Real trqLimit  = getTorqueLimit(state);
        const Real abstrq=std::abs(integTrq), trqsign = sign(integTrq);

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

    // Topology cache
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
    MyMechanism() 
    :   system(), matter(system), forces(system), viz(system), userInput(0),
        controller(0)
    {
        constructSystem();
        setUpVisualizer();
    }

    const MultibodySystem& getSystem() const {return system;}

    const State& getDefaultState() const
    {   return system.getDefaultState(); }

    Real getDesiredSpeed(const State& state) const 
    {   return controller->getDesiredSpeed(state); }
    Real getTorqueLimit(const State& state) const
    {   return controller->getTorqueLimit(state); }
    Real getActualSpeed(const State& state) const 
    {   return controller->getActualSpeed(state); }
    Real getSpeedError(const State& state) const
    {   return controller->getSpeedError(state); }

    Real getMotorTorque(const State& state) const
    {   return controller->getTorque(state); }

    // Change the desired speed.
    void changeSpeed(State& state, Real newSpeed) const;

    // Change the maximum torque limit.
    void changeTorqueLimit(State& state, Real newTorqueLimit) const;

    // Update the Visualizer, including the sliders.
    void draw(const State& state) const;

    // Returns true when it is time to quit.
    bool processUserInput(State& state) const;

private:
    void constructSystem();
    void setUpVisualizer();

    MultibodySystem                 system;
    SimbodyMatterSubsystem          matter;
    GeneralForceSubsystem           forces;
    Visualizer                      viz;
    Visualizer::InputSilo*          userInput; // reference only

    MobilizedBody                   bodyT, leftArm, rtArm;

    Force::Gravity                  gravity;
    MyTorqueLimitedMotor*           controller; // reference only
};

}


//==============================================================================
//                                  MAIN
//==============================================================================
int main() {
    try { // catch errors if any

    // Create the Simbody MultibodySystem and Visualizer.
    MyMechanism mech;

    // Set initState to the System's default value.    
    State initState = mech.getDefaultState();

    // Simulate forever with a step size as big as we can take while updating
    // the screen fast enough. Check for user input in between steps. Note: an 
    // alternate way to do this is to let the integrator take whatever steps it 
    // wants, by using a TimeStepper to manage two events: a 
    // PeriodicEventReporter to draw frames however fast we want using 
    // interpolated values, and a periodic event handler to poll for user input
    // less frequently (probably 100ms would be often enough). Here we treat
    // step completion as an event and update screen & user input between steps.
    const Real MaxStepSize = 0.030;  // 30ms
    const int  DrawEveryN  = 1;      // 1 steps = 30ms
    //RungeKuttaMersonIntegrator integ(mech.getSystem());
    //RungeKutta3Integrator integ(mech.getSystem());
    SemiExplicitEuler2Integrator integ(mech.getSystem());
    //SemiExplicitEulerIntegrator integ(mech.getSystem(), .001);

    integ.setAccuracy(1e-1);
    //integ.setAccuracy(1e-3);

    // Don't permit interpolation because we want the state returned after
    // a step to be modifiable.
    integ.setAllowInterpolation(false);

    integ.initialize(initState);
    State& state = integ.updAdvancedState();

    int stepsSinceViz = DrawEveryN-1;
    while (true) {
        if (++stepsSinceViz % DrawEveryN == 0) {
            mech.draw(state);
            stepsSinceViz = 0;
        }

        // Advance time by MaxStepSize (might take multiple steps to get there).
        // Note: interpolation must be disabled if we're going to modify the
        // state at the end of each step.
        integ.stepBy(MaxStepSize);

        // Check for user input before resuming simulation.
        if (mech.processUserInput(state))
            break;
    }

    // DONE. Dump out final stats.

    const int evals = integ.getNumRealizations();
    std::cout << "Done -- simulated " << integ.getTime() << "s with " 
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

    } catch (const std::exception& e) {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}



//==============================================================================
//                        MY MECHANISM IMPLEMENTATION
//==============================================================================

//----------------------------- CONSTRUCT SYSTEM -------------------------------
void MyMechanism::constructSystem() {
    gravity = Force::Gravity(forces, matter, -YAxis, 9.80665);

    // Describe a body with a point mass at (0, -3, 0) and draw a sphere there.
    Real mass = 3; Vec3 pos(0,-3,0);
    Body::Rigid bodyInfo(MassProperties(mass, pos, UnitInertia::pointMassAt(pos)));
    bodyInfo.addDecoration(pos, DecorativeSphere(.2).setOpacity(.5));

    // Create the tree of mobilized bodies, reusing the above body description.
    bodyT   = MobilizedBody::Pin(matter.Ground(), Vec3(0), bodyInfo, Vec3(0));
    leftArm = MobilizedBody::Pin(bodyT, Vec3(-2, 0, 0),    bodyInfo, Vec3(0));
    rtArm   = MobilizedBody::Pin(bodyT, Vec3(2, 0, 0),     bodyInfo, Vec3(0,-1,0));

    // Add some damping.
    Force::MobilityLinearDamper damper1(forces, bodyT, MobilizerUIndex(0), 10);
    Force::MobilityLinearDamper damper2(forces, leftArm, MobilizerUIndex(0), 30);
    Force::MobilityLinearDamper damper3(forces, rtArm, MobilizerUIndex(0), 10);

    // Add a joint stop to the left arm restricting it to q in [0,Pi/5].
    Force::MobilityLinearStop stop(forces, leftArm, MobilizerQIndex(0), 
        StopStiffness,
        StopDissipation,
        -Pi/8,   // lower stop
         Pi/8);  // upper stop

    // Use custom controller to track speed up to torque limit.
    controller = new MyTorqueLimitedMotor(rtArm, MobilizerUIndex(0));
    Force::Custom(forces, controller); // takes over ownership

    // We're done with the System; finalize it.
    system.realizeTopology();
}

//------------------------------- CHANGE SPEED ---------------------------------
// Solve GM\~G lambda = deltaV to calculate constraint impulse lambda.
// Then f = ~G lambda is the corresponding generalized impulse.
// And deltaU = M\f is the resulting change to the generalized speeds.
void MyMechanism::changeSpeed(State& state, Real newSpeed) const {
    const Real oldSpeed = getActualSpeed(state);
    printf("Desired speed change: %g -> %g\n", oldSpeed, newSpeed);
    controller->setDesiredSpeed(state, newSpeed);
}


//--------------------------- CHANGE TORQUE LIMIT ------------------------------
void MyMechanism::changeTorqueLimit(State& state, Real newTorqueLimit) const {
    const Real oldTorqueLimit = getTorqueLimit(state);
    printf("Torque limit change: %g -> %g\n", oldTorqueLimit, newTorqueLimit);
    controller->setTorqueLimit(state, newTorqueLimit);
}


//---------------------------- PROCESS USER INPUT ------------------------------
// Return true if it is time to quit.
bool MyMechanism::processUserInput(State& state) const {
    int whichSlider, whichMenu, whichItem; Real newValue;

    // Did a slider move?
    if (userInput->takeSliderMove(whichSlider, newValue)) {
        switch(whichSlider) {
        case SliderIdMotorSpeed:
            // This will momentum balance if necessary.
            changeSpeed(state, newValue);
            break;
        case SliderIdTorqueLimit:
            changeTorqueLimit(state, newValue);
            viz.setSliderRange(SliderIdTorque, -newValue, newValue); 
            break;
        }
    }

    // Was there a menu pick?
    if (userInput->takeMenuPick(whichMenu, whichItem)) {
        if (whichItem == QuitItem) 
            return true; // done

        // If Reset, stop the motor and zero out all the q's and u's. 
        // Tell visualizer to update the sliders to match.
        if (whichItem == ResetItem) {
            // Don't momentum balance here!
            controller->setDesiredSpeed(state, 0);
            viz.setSliderValue(SliderIdMotorSpeed, 0);
            viz.setSliderValue(SliderIdTach, 0);

            controller->setTorqueLimit(state, InitialTorqueLimit);
            viz.setSliderValue(SliderIdTorqueLimit, InitialTorqueLimit);
            viz.setSliderRange(SliderIdTorque, -InitialTorqueLimit, 
                                                InitialTorqueLimit); 
            viz.setSliderValue(SliderIdTorque, 0);

            state.updQ() = 0; // all positions to zero
            state.updU() = 0; // all velocities to zero
            system.projectQ(state);
            system.projectU(state);
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
    viz.setShutdownWhenDestructed(true) // make sure display window dies
       .setBackgroundType(Visualizer::SolidColor); // turn off Ground & Sky
    
    // Add sliders.
    viz.addSlider("Motor speed", SliderIdMotorSpeed, 
                  -MaxMotorSpeed, MaxMotorSpeed, InitialMotorSpeed);
    viz.addSlider("Torque limit", SliderIdTorqueLimit, 
                  0, MaxTorqueLimit, InitialTorqueLimit);
    viz.addSlider("Tach",   SliderIdTach,   
                  -MaxMotorSpeed,  MaxMotorSpeed,  0);
    viz.addSlider("Torque", SliderIdTorque, 
                  -InitialTorqueLimit, InitialTorqueLimit, 0);

    // Add Run menu.
    Array_<std::pair<String,int> > runMenuItems;
    runMenuItems.push_back(std::make_pair("Reset", ResetItem));
    runMenuItems.push_back(std::make_pair("Quit", QuitItem));
    viz.addMenu("Run", MenuIdRun, runMenuItems);

    // Add per-frame text and geometry.
    viz.addDecorationGenerator(new ShowStuff(*this));

    // Add an input listener so the user can talk to us.
    userInput = new Visualizer::InputSilo();
    viz.addInputListener(userInput);
}

//------------------------------------ DRAW ------------------------------------
// Update the Visualizer, including the sliders.
void MyMechanism::draw(const State& state) const {
    viz.report(state);
    viz.setSliderValue(SliderIdTach, getActualSpeed(state));
    viz.setSliderValue(SliderIdTorque, getMotorTorque(state));
}
