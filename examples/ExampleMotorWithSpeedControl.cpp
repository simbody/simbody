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
#include <iostream>
#include <algorithm>
using std::cout; using std::endl;

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
const int SliderIdMotorSpeed = 1, SliderIdDissipation = 2, 
    SliderIdTach = 3, SliderIdTorque = 4;

// Ids for things on the Run Menu.
const int MenuIdRun = 1;
static const int ResetItem=1, QuitItem=2;

#define USE_TORQUE_LIMITED_MOTOR
const Real MaxTorque = 100; // N-m
const Real MaxTorqueRate = .5; // s / MaxTorque change
const Real TorqueGain = 50000, DampingGain = 1000;
const Real TorqueDecay = 100;

#ifdef USE_TORQUE_LIMITED_MOTOR
// This Force element implements a torque-limited motor that runs at a user-
// selected speed unless that would require too much torque.
//     trqDot = gain*(u_desired - u_actual)
//     if (trq > MaxTorque) trqDot = min(0, trqDot)
//     else if (trq < -MaxTorque) trqDot = max(0, trqDot)
//

class MyTorqueLimitedMotor : public Force::Custom::Implementation {
public:
    MyTorqueLimitedMotor
       (const MobilizedBody& mobod, MobilizerUIndex whichU, 
        Real gain, Real torqueLimit)
    :   m_matter(mobod.getMatterSubsystem()), m_mobod(mobod), m_whichU(whichU), 
        m_torqueGain(gain), m_torqueLimit(torqueLimit) 
    {   assert(gain >= 0 && torqueLimit >= 0); }

    void setDesiredRate(State& state, Real speed) const {
        Real& u_desired = Value<Real>::updDowncast(
            m_matter.updDiscreteVariable(state, m_desiredUIx));
        u_desired = speed;
    }
    // synonym to match name used by Motion
    void setRate(State& state, Real speed) const {setDesiredRate(state,speed);}
    // synonym to match name used by ConstantSpeed constraint
    void setSpeed(State& state, Real speed) const {setDesiredRate(state,speed);}

    Real getDesiredRate(const State& state) const {
        const Real u_desired = Value<Real>::downcast(
            m_matter.getDiscreteVariable(state, m_desiredUIx));
        return u_desired;
    }

    Real getActualRate(const State& state) const {
        const Real u_actual = m_mobod.getOneU(state, m_whichU);
        return u_actual;
    }

    Real getRateError(const State& state) const {
        return getActualRate(state) - getDesiredRate(state);
    }

    // Combine integrated torque and proportional torque, and clamp to limit.
    Real getTorque(const State& state) const {
        const Real integTrq  = m_matter.getZ(state)[m_trqIx];
        const Real proportionalTrq = -DampingGain*getRateError(state);
        const Real totalTrq = clamp(-m_torqueLimit, integTrq+proportionalTrq, 
                                     m_torqueLimit);
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

    void realizeTopology(State& state) const OVERRIDE_11 {
        m_desiredUIx = m_matter.allocateDiscreteVariable
           (state, Stage::Acceleration, new Value<Real>(0));
        m_trqIx = m_matter.allocateZ(state, Vector(1,Real(0)));
    }

    void realizeAcceleration(const State& state) const OVERRIDE_11 {
        const Real integTrq = m_matter.getZ(state)[m_trqIx];
        const Real abstrq=std::abs(integTrq), trqsign = sign(integTrq);

        Real trqDot = -m_torqueGain*getRateError(state);
        // Don't ask for more torque if we're already past the limit
        if (abstrq >= m_torqueLimit && sign(trqDot)*trqsign >= 0)
            trqDot = 0;

        //const char* delay="WantLess: ";
        //if (sign(trqDot)*trqsign >= 0) {
        //    // want more torque
        //    if (abstrq > m_torqueLimit) {
        //        trqDot = -trqsign*TorqueDecay*(abstrq-m_torqueLimit);
        //        delay =   "TooBig:   ";
        //    } else if (abstrq > 0.9*m_torqueLimit) {
        //        trqDot *= (m_torqueLimit-abstrq)/(0.9*m_torqueLimit);
        //        delay =   "Limited:  ";
        //    } else
        //        delay =   "WantMore: ";
        //}

        //printf("%suerr=%g (actual=%g, des=%g) trq=%g trqdot=%g\n", 
        //    delay, u_actual-u_desired,
        //    u_actual, u_desired, trq, trqDot);

        Real rateLimit = m_torqueLimit/MaxTorqueRate; // full-scale change
        //clampInPlace(-rateLimit, trqDot, rateLimit);

        m_matter.updZDot(state)[m_trqIx] = trqDot; 
    }

private:
    const SimbodyMatterSubsystem& m_matter;
    MobilizedBody   m_mobod;
    MobilizerUIndex m_whichU;
    Real            m_torqueGain;
    Real            m_torqueLimit;

    // Topology cache
    mutable DiscreteVariableIndex m_desiredUIx;
    mutable ZIndex                m_trqIx;
};
#endif

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
    Force::MobilityLinearDamper damper1(forces, bodyT, MobilizerUIndex(0), 10);
    Force::MobilityLinearDamper damper2(forces, leftArm, MobilizerUIndex(0), 30);
    Force::MobilityLinearDamper damper3(forces, rtArm, MobilizerUIndex(0), 10);

#ifdef USE_TORQUE_LIMITED_MOTOR
    MyTorqueLimitedMotor* motorp = 
        new MyTorqueLimitedMotor(rtArm, MobilizerUIndex(0), TorqueGain, MaxTorque);
    const MyTorqueLimitedMotor& motor = *motorp;
    Force::Custom(forces, motorp); // takes over ownership
#else
    // Use built-in Steady Motion as a low-budget motor model.
    //Motion::Steady motor(rtArm, InitialMotorRate);
    // Use built-in ConstantSpeed constraint as a low-budget motor model.
    Constraint::ConstantSpeed motor(rtArm, InitialMotorRate);
#endif

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
    viz.addSlider("Tach", SliderIdTach, -20, 20, 0);
    viz.addSlider("Torque", SliderIdTorque, -MaxTorque, MaxTorque, 0);

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
    const Real MaxStepSize = 0.01*3; // 10ms
    const int  DrawEveryN = 3/3;     // 3 steps = 30ms
    //RungeKuttaMersonIntegrator integ(system);
    //RungeKutta2Integrator integ(system);
    SemiExplicitEuler2Integrator integ(system);
    //SemiExplicitEulerIntegrator integ(system, .001);

    integ.setAccuracy(1e-1);
    //integ.setAccuracy(1e-3);

    // Don't permit interpolation because we want the state returned after
    // a step to be modifiable.
    integ.setAllowInterpolation(false);

    integ.initialize(initState);

    int stepsSinceViz = DrawEveryN-1;
    while (true) {
        if (++stepsSinceViz % DrawEveryN == 0) {
            const State& s = integ.getState();
            viz.report(s);
            const Real uActual = rtArm.getOneU(s, MobilizerUIndex(0));
            viz.setSliderValue(SliderIdTach, uActual);
#ifdef USE_TORQUE_LIMITED_MOTOR
            viz.setSliderValue(SliderIdTorque, motor.getTorque(s));
#else
            system.realize(s); // taus are acceleration stage
            //viz.setSliderValue(SliderIdTorque, 
            //                   rtArm.getOneTau(s, MobilizerUIndex(0)));
            viz.setSliderValue(SliderIdTorque, motor.getMultiplier(s));
#endif

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
                //motor.setRate(state, newValue);
                motor.setSpeed(state, newValue);
                system.realize(state, Stage::Position);
                system.prescribeU(state);
                system.realize(state, Stage::Velocity);
                system.projectU(state);
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
                //motor.setRate(state, 0);
                motor.setSpeed(state, 0);
                viz.setSliderValue(SliderIdMotorSpeed, 0);

                stop.setMaterialProperties(state, StopStiffness, InitialDissipation);
                viz.setSliderValue(SliderIdDissipation, InitialDissipation);

                state.updQ() = 0; // all positions to zero
                state.updU() = 0; // all velocities to zero
                system.realize(state, Stage::Position);
                system.prescribeU(state);
                system.realize(state, Stage::Velocity);
                system.projectU(state);
            }
        }

    }
    const int evals = integ.getNumRealizations();
    std::cout << "Done -- simulated " << integ.getTime() << "s with " 
            << integ.getNumStepsTaken() << " steps, avg step=" 
        << (1000*integ.getTime())/integ.getNumStepsTaken() << "ms " 
        << (1000*integ.getTime())/evals << "ms/eval\n";

    printf("Used Integrator %s at accuracy %g:\n", 
        integ.getMethodName(), integ.getAccuracyInUse());
    printf("# STEPS/ATTEMPTS = %d/%d\n", integ.getNumStepsTaken(), integ.getNumStepsAttempted());
    printf("# ERR TEST FAILS = %d\n", integ.getNumErrorTestFailures());
    printf("# REALIZE/PROJECT = %d/%d\n", integ.getNumRealizations(), integ.getNumProjections());

    } catch (const std::exception& e) {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
