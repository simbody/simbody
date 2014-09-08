/* -------------------------------------------------------------------------- *
 *                  Simbody(tm) Example: Operational Space                    *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-14 Stanford University and the Authors.        *
 * Authors: Chris Dembia, Jack Wang                                           *
 * Contributors: Michael Sherman                                              *
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

/**
 * This example shows how one could write an operational space controller in
 * Simbody. The model we are controlling is a human upper body model.
 * Operational space controllers are model-based, and we show how to access the
 * necessary quantities from the SimbodyMatterSubsystem. We assume the system
 * has a motor at each of its degrees of freedom.
 *
 * The task the controller will achieve has two components:
 * 1. One arm reaches for a target point
 * 2. All bodies are subject to gravity compensation (to counteract the effect
 * of gravity).
 * 
 *      TODO inertial forces vector, in TaskSpace (1)
 *      TODO cache computations (2)
 * TODO make computations efficient (3)
 * TODO document TaskSpace
 * TODO write missing method in Simbody to return a Vector.
 * TODO put nullspace subtraction elsewhere.
 * TODO missing (convenience) operators.
 * TODO separate gravity into its own quantity.
 * TODO how to avoid crashing with singularities.
 */

#include "Humanoid.h"
#include <Simbody.h>

using namespace SimTK;
using namespace std;

#define DEG(rad) convertRadiansToDegrees(rad)
#define RAD(deg) convertDegreesToRadians(deg)

//==============================================================================
// ReachingAndGravityCompensation
//==============================================================================
class ReachingAndGravityCompensation :
    public SimTK::Force::Custom::Implementation
{
public:

    /// The left hand reaches for a target location, and the controller
    /// compensates for gravity.
    /// 
    /// @param[in] stationLocationInHand The point whose location we want
    ///                                  to control.
    /// @param[in] proportionalGain Units of N-m/rad
    /// @param[in] derivativeGain Units of N-m-s/rad
    ReachingAndGravityCompensation(const Humanoid& system,
            Vec3 stationLocationInHand=Vec3(0.06, -0.04, 0),
            double proportionalGain=100, double derivativeGain=20) :
        m_system(system), m_matter(system.getMatterSubsystem()),
        m_tspace1(m_matter, system.getGravity()),
        m_tspace2(m_matter, system.getGravity()),
        m_stationLocationInHand(stationLocationInHand),
        m_proportionalGain(proportionalGain),
        m_derivativeGain(derivativeGain),
        m_desiredLeftPosInGround(Vec3(0.4, 1.5, -0.1)),
        m_desiredRightPosInGround(Vec3(0.4, 1.8, +0.2))
    {
        m_tspace1.addStationTask(m_system.getBody(Humanoid::hand_l),
                         m_stationLocationInHand);
        m_tspace2.addStationTask(m_system.getBody(Humanoid::hand_r),
                         m_stationLocationInHand);
    }

    void calcForce(const SimTK::State&                state,
                   SimTK::Vector_<SimTK::SpatialVec>& bodyForces,
                   SimTK::Vector_<SimTK::Vec3>&       particleForces,
                   SimTK::Vector&                     mobilityForces) const
                   OVERRIDE_11;

    Real calcPotentialEnergy(const SimTK::State& state) const OVERRIDE_11
    { return 0; }

    void calcDecorativeGeometryAndAppend(const State & state, Stage stage,
            Array_<DecorativeGeometry>& geometry) const OVERRIDE_11
    {
        geometry.push_back(DecorativeSphere(0.02)
                .setTransform(getTarget())
                .setColor(m_targetColor));

        const MobilizedBody& leftHand = m_system.getBody(Humanoid::hand_l);
        Vec3 leftPosInGround = leftHand.findStationLocationInGround(state,
                m_stationLocationInHand);
        geometry.push_back(DecorativeSphere(0.02)
                .setTransform(leftPosInGround)
                .setColor(Vec3(0, 1, 0)));

        const MobilizedBody& rightHand = m_system.getBody(Humanoid::hand_r);
        Vec3 rightPosInGround = rightHand.findStationLocationInGround(state,
                m_stationLocationInHand);
        geometry.push_back(DecorativeSphere(0.02)
                .setTransform(rightPosInGround)
                .setColor(Vec3(0, 1, 0)));
    }

    void realizeTopology(State& state) const OVERRIDE_11
    {
        m_tspace1.realizeTopology(state);
        m_tspace2.realizeTopology(state);
    }

    const Vec3& getTarget() const { return m_desiredLeftPosInGround; }
    Vec3& updTarget() { return m_desiredLeftPosInGround; }
    void setTarget(Vec3 pos) { m_desiredLeftPosInGround = pos; }

private:

    const Humanoid& m_system;
    const SimbodyMatterSubsystem& m_matter;
    TaskSpace m_tspace1;
    TaskSpace m_tspace2;
    const Vec3 m_stationLocationInHand;
    const double m_proportionalGain;
    const double m_derivativeGain;
    Vec3 m_desiredLeftPosInGround;
    Vec3 m_desiredRightPosInGround;
    static const unsigned int m_numTasks = 3;
    static const double m_dampingGain;
    static const Vec3 m_targetColor;
};

const Vec3 ReachingAndGravityCompensation::m_targetColor = Vec3(1, 0, 0);
const double ReachingAndGravityCompensation::m_dampingGain = 1;

/*
void ReachingAndGravityCompensation::calcForce(
               const State&                state,
               Vector_<SimTK::SpatialVec>& bodyForces,
               Vector_<SimTK::Vec3>&       particleForces,
               Vector&                     mobilityForces) const
{
    const int nu = state.getNU();
    const int m = m_numTasks;

    // Compute control law in task space (F*).
    // ---------------------------------------
    // These are approximations, assuming the user is not moving the target too
    // quickly.
    Vec3 desiredVelInGround(0);
    Vec3 desiredAccInGround(0);

    // Get info about the actual location, etc. of the left hand.
    Vec3 posInGround;
    Vec3 velInGround;
    const MobilizedBody& leftHand = m_system.getBody(Humanoid::hand_l);
    leftHand.findStationLocationAndVelocityInGround(state,
            m_stationLocationInLeftHand, posInGround, velInGround);

    // Units of acceleration.
    Vector_<Real> Fstar(desiredAccInGround +
        m_derivativeGain * (desiredVelInGround - velInGround) +
        m_proportionalGain * (m_desiredPosInGround - posInGround));

    // Task jacobian.
    // --------------
    TaskSpace ts(m_matter);
    ts.addTask(leftHand, m_stationLocationInLeftHand);

    const TaskSpace::Jacobian& J = ts.jacobian();
    const TaskSpace::Inertia& Lambda = ts.inertia();
    // TODO const TaskSpace::CoriolisForce& mu = ts.coriolisForce();
    const TaskSpace::GravityForce& p = ts.gravityForce();
    const TaskSpace::NullspaceProjection& N = ts.nullspaceProjection();
    const TaskSpace::JacobianTranspose Jt = J.transpose();
    const TaskSpace::NullspaceProjectionTranspose = N.transpose();

    // Gravity (g).
    // ------------
    // Get the gravity vector, for gravity compensation.
    Vector g;
    m_matter.multiplyBySystemJacobianTranspose(state,
            m_system.getGravity().getBodyForces(state),
            g);

    // Compute task-space force that achieves the task-space control.
    // F = Lambda F*  + mu + p
    Vector F = Lambda * Fstar + mu + p;
    // Vector F = ts.calcTaskSpaceForceFromTaskSpaceAcceleration(Fstar);
    // Vector F = ts.calcForwardDynamics(Fstar);

    // Combine the reaching task with the gravity compensation.
    mobilityForces = J.transpose() * F  + N.transpose() * (-g);
}
*/

void ReachingAndGravityCompensation::calcForce(
               const State&                state,
               Vector_<SimTK::SpatialVec>& bodyForces,
               Vector_<SimTK::Vec3>&       particleForces,
               Vector&                     mobilityForces) const
{
    // Shorthands.
    // -----------
    const TaskSpace& p1 = m_tspace1;
    const TaskSpace& p2 = m_tspace2;
    const State& s = state;

    const int nu = state.getNU();
    const int m = m_numTasks;

    const Real& kd = m_derivativeGain;
    const Real& kp = m_proportionalGain;
    const Vec3& x1_des = m_desiredLeftPosInGround;
    const Vec3& x2_des = m_desiredRightPosInGround;

    p1.setState(state);
    p2.setState(state);

    const TaskSpace::Jacobian& J = p1.getJacobian(s);
    J.value2();

    // Compute control law in task space (F*).
    // ---------------------------------------
    Vec3 xd_des(0);
    Vec3 xdd_des(0);

    // Get info about the actual location, etc. of the left hand.
    Vec3 x1, x1d;
    const MobilizedBody& leftHand = m_system.getBody(Humanoid::hand_l);
    leftHand.findStationLocationAndVelocityInGround(state,
            m_stationLocationInHand, x1, x1d);
    Vec3 x2, x2d;
    const MobilizedBody& rightHand = m_system.getBody(Humanoid::hand_r);
    rightHand.findStationLocationAndVelocityInGround(state,
            m_stationLocationInHand, x2, x2d);

    // Units of acceleration.
    Vec3 Fstar1 = xdd_des + kd * (xd_des - x1d) + kp * (x1_des - x1);
    Vector Fstar2(xdd_des + kd * (xd_des - x2d) + kp * (x2_des - x2));

    // Compute task-space force that achieves the task-space control.
    // F = Lambda Fstar + p
    Vector F1 = p1.Lambda() * Fstar1 + p1.mu() + p1.p();
    Vector F2 = p2.calcInverseDynamics(Fstar2);

    // Combine the reaching task with the gravity compensation and nullspace
    // damping.
    // Gamma = J1T F1 + N1T J1T F2 + N1T N2T (g - k u)
    const Vector& u = state.getU();
    const double& k = m_dampingGain;
    mobilityForces = p1.JT() * F1 + p1.NT() * (
                        p2.JT() * F2 + p2.NT() * (
                            p1.g() - k * u)
                        );

//    TODO what would this look like with states as arguments everywhere?
//    // Compute task-space force that achieves the task-space control.
//    // F = Lambda Fstar + p
//    Vector F1 = p1.Lambda(s) * Fstar1 + /* TODO p1.mu() + */ p1.p(s);
//    Vector F2 = p2.calcInverseDynamics(Fstar2);
//
//    // Combine the reaching task with the gravity compensation and nullspace
//    // damping.
//    // Gamma = J1T F1 + N1T J1T F2 + N1T N2T (g - k u)
//    const Vector& u = state.getU();
//    const double& k = m_dampingGain;
//    mobilityForces = p1.JT(s) * F1 +
//                     p1.NT(s) * (p2.JT(s) * F2) +
//                     p1.NT(s) * (p2.NT(s) * (p1.g(s) - k * u));
}

//==============================================================================
// UserInputHandler
//==============================================================================
/// This is a periodic event handler that interrupts the simulation on a
/// regular basis to poll the InputSilo for user input.
class UserInputHandler : public PeriodicEventHandler {
public:
    UserInputHandler(Visualizer::InputSilo& silo,
            ReachingAndGravityCompensation& controller, Real interval)
        :
        PeriodicEventHandler(interval), m_silo(silo), m_controller(controller),
        m_increment(0.05) {}
    void handleEvent(State& state, Real accuracy,
            bool& shouldTerminate) const OVERRIDE_11
    {
        while (m_silo.isAnyUserInput()) {
            unsigned key, modifiers;
            while (m_silo.takeKeyHit(key, modifiers)) {
                if (key == Visualizer::InputListener::KeyEsc) {
                    shouldTerminate = true;
                    m_silo.clear();
                    return;
                }
                else if (key == Visualizer::InputListener::KeyRightArrow) {
                    // z coordinate goes in and out of the screen.
                    m_controller.updTarget()[ZAxis] -= m_increment;
                    m_silo.clear();
                    return;
                }
                else if (key == Visualizer::InputListener::KeyLeftArrow) {
                    m_controller.updTarget()[ZAxis] += m_increment;
                    m_silo.clear();
                    return;
                }
                else if (key == Visualizer::InputListener::KeyUpArrow) {
                    m_controller.updTarget()[YAxis] += m_increment;
                    m_silo.clear();
                    return;
                }
                else if (key == Visualizer::InputListener::KeyDownArrow) {
                    m_controller.updTarget()[YAxis] -= m_increment;
                    m_silo.clear();
                    return;
                }
                else if (key == Visualizer::InputListener::KeyPageUp) {
                    m_controller.updTarget()[XAxis] += m_increment;
                    m_silo.clear();
                    return;
                }
                else if (key == Visualizer::InputListener::KeyPageDown) {
                    m_controller.updTarget()[XAxis] -= m_increment;
                    m_silo.clear();
                    return;
                }
            }
        }
    }
private:
    Visualizer::InputSilo& m_silo;
    ReachingAndGravityCompensation& m_controller;
    const Real m_increment;
};

//==============================================================================
// OutputReporter
//==============================================================================
/// This is a periodic event handler that interrupts the simulation on a
/// regular basis so that we can display information.
class OutputReporter : public PeriodicEventReporter {
public:
    OutputReporter(const Humanoid& system, Real interval) :
        PeriodicEventReporter(interval), m_system(system) {}
    void handleEvent(const State& state) const OVERRIDE_11 {}
private:
    const Humanoid& m_system;
};

//==============================================================================
// main function
//==============================================================================
int main(int argc, char **argv)
{
    try {

        // Set some options.
        const double duration = Infinity; // seconds.
        const double realTimeScale = 1.0; // unitless.
        const double frameRate = 30.0; // frames per second.

        // Create the system.
        Humanoid system;

        // Add the controller.
        ReachingAndGravityCompensation* controller =
            new ReachingAndGravityCompensation(system);
        // Force::Custom takes ownership over controller.
        Force::Custom control(system.updForceSubsystem(), controller);

        // Set up visualizer and event handlers.
        SimTK::Visualizer viz(system);
        viz.setMode(Visualizer::RealTime);
        viz.setRealTimeScale(realTimeScale);
        viz.setShowFrameRate(true);
        viz.setShowSimTime(true);
        SimTK::Visualizer::InputSilo* userInput = new Visualizer::InputSilo();
        viz.addInputListener(userInput);
        system.addEventHandler(
                new UserInputHandler(*userInput, *controller, 0.1));
        system.addEventReporter(new OutputReporter(system, .01));
        system.addEventReporter(
                new Visualizer::Reporter(viz, realTimeScale / frameRate));

        // Display message on the screen about how to start simulation.
        DecorativeText help("Any input to start; ESC to quit. "
                "Move target with arrow keys, PageUp and PageDown.");
        help.setIsScreenText(true);
        viz.addDecoration(MobilizedBodyIndex(0), Vec3(0), help);

        // Initialize the system and other related classes.
        State s;
        system.initialize(s);
        SemiExplicitEuler2Integrator integ(system);
        integ.setAccuracy(0.001);
        TimeStepper ts(system, integ);
        ts.initialize(s);
        viz.report(ts.getState());

        /* TODO
           printf("Act: u\n");
           for (int i=0; i < NumActuators; ++i) {
           printf("%2d: %d\n", i, int(system.getUIndex(Humanoid::Coordinate(i))));
           }
           */

        userInput->waitForAnyUserInput();
        userInput->clear();

        const double startCPU  = cpuTime(), startTime = realTime();

        // Simulate.
        ts.stepTo(duration);

        std::cout << "CPU time: " << cpuTime() - startCPU << " seconds. "
                  << "Real time: " << realTime() - startTime << " seconds."
                  << std::endl;

    } catch (const std::exception& e) {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
