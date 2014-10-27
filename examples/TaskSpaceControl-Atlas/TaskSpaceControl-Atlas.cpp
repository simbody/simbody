/* -------------------------------------------------------------------------- *
 *               Simbody(tm) Example: Atlas Task Space Control                *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
 * Authors: Michael Sherman, Chris Dembia                                     *
 * Contributors: John Hsu                                                     *
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

/* This example shows a task space (a.k.a. operational space) controller for
the Gazebo model of the Boston Dynamics Atlas humanoid robot used in the 
DARPA Robotics Challenge. Task space controllers are model-based, meaning that
the controller itself contains a model of the system being controlled. Since
we don't have a real robot handy, that means there will be two versions of
the Atlas system here: one that we're simulating that is taking the role of the
"real" robot, and one contained in the task space controller that we'll call
the "model" robot. In real life the internal model won't perfectly match the 
real one; we'll fake that here by introducing some sensor noise which you can
control with sliders in the user interface.

We assume the system has a direct-drive motor at each of its degrees of freedom.

The task the controller will achieve has two components:
1. An arm reaches for a target point that can be moved with arrow keys.
2. All links are subject to gravity compensation (to counteract the effect
of gravity).

You can also optionally sense the end effector position on the real robot
and have that sent to the controller so that it doesn't have to depend 
entirely on the behavior of the model robot when the real robot's sensors are
noisy. Try cranking up the noise, which causes poor tracking, and then hit
"e" to enable the end effector sensing which improves things dramatically.

For more information about operational space control, see:

Khatib, Oussama, et al. "Robotics-based synthesis of human motion."
Journal of physiology-Paris 103.3 (2009): 211-219.
*/

#include "Simbody.h"
#include "Atlas.h"

#include <iostream>

using namespace SimTK;
using namespace std;

static const int QNoise=1, UNoise=2; // sliders in the UI

//==============================================================================
//                             TASKS MEASURE
//==============================================================================
// This Measure is added to the modelRobot and is used to manage the tasks
// it is supposed to achieve and to return as its value the control torques
// that should be applied to the realRobot (that is, the simulated one).
// This should only be instantiated with T=Vector.
template <class T>
class TasksMeasure : public Measure_<T> {
public:
    SimTK_MEASURE_HANDLE_PREAMBLE(TasksMeasure, Measure_<T>);

    TasksMeasure(Atlas& modelRobot) 
    :   Measure_<T>(modelRobot.updForceSubsystem(), 
                    new Implementation(modelRobot),
                    AbstractMeasure::SetHandle()) {}


    const Vec3& getTarget() const { return getImpl().m_desiredTaskPosInGround; }
    Vec3& updTarget() { return updImpl().m_desiredTaskPosInGround; }
    void setTarget(Vec3 pos) { updImpl().m_desiredTaskPosInGround = pos; }

    void toggleGravityComp() {
        updImpl().m_compensateForGravity = !isGravityCompensationOn();}
    void toggleTask() {updImpl().m_controlTask = !getImpl().m_controlTask;}
    void toggleEndEffectorSensing() 
    {   updImpl().m_endEffectorSensing = !getImpl().m_endEffectorSensing;}

    bool isGravityCompensationOn() const 
    {   return getImpl().m_compensateForGravity; }
    bool isEndEffectorSensingOn() const 
    {   return getImpl().m_endEffectorSensing; }
    bool isTaskPointFollowingOn() const
    {   return getImpl().m_controlTask; }
    const Vec3& getTaskPointInEndEffector() const 
    {   return getImpl().m_taskPointInEndEffector; }

    SimTK_MEASURE_HANDLE_POSTSCRIPT(TasksMeasure, Measure_<T>);
};


template <class T>
class TasksMeasure<T>::Implementation : public Measure_<T>::Implementation {
public:
    Implementation(const Atlas& modelRobot,
                   Real proportionalGain=225, double derivativeGain=30) 
                   //Real proportionalGain=100, double derivativeGain=20) 
    :   Measure_<T>::Implementation(T(), 1),
        m_modelRobot(modelRobot),
        m_tspace1(m_modelRobot.getMatterSubsystem(), m_modelRobot.getGravity()),
        m_proportionalGain(proportionalGain),
        m_derivativeGain(derivativeGain),
        m_jointPositionGain(100),
        m_jointDampingGain(20),
        m_compensateForGravity(true),
        m_controlTask(false),
        m_endEffectorSensing(false),
        m_desiredTaskPosInGround(Vec3(0.4, 0.1, 1)) // Z is up
    {       
        //TODO: should have end effector body
        m_tspace1.addStationTask(m_modelRobot.getEndEffectorBody(),
                                 m_modelRobot.getEndEffectorStation());
    }


    // Default copy constructor, destructor, copy assignment are fine.

    // Implementations of virtual methods.
    Implementation* cloneVirtual() const OVERRIDE_11
    {   return new Implementation(*this); }
    int getNumTimeDerivativesVirtual() const OVERRIDE_11 {return 0;}
    Stage getDependsOnStageVirtual(int order) const OVERRIDE_11
    {   return Stage::Velocity; }

    // This is the task space controller. It returns joint torques tau as the
    // value of the enclosing Measure.
    void calcCachedValueVirtual(const State& s, int derivOrder, T& tau) const
        OVERRIDE_11;

    // TaskSpace objects require some State resources; this call is the time
    // for doing that so forward on to the TaskSpace.
    void realizeMeasureTopologyVirtual(State& modelState) const OVERRIDE_11 {
        m_tspace1.realizeTopology(modelState);
    }
private:
friend class TasksMeasure<T>;

    const Atlas&    m_modelRobot;
    TaskSpace       m_tspace1;

    const Real      m_proportionalGain;     // task space
    const Real      m_derivativeGain;
    const Real      m_jointPositionGain;    // joint space
    const Real      m_jointDampingGain;

    bool            m_compensateForGravity;
    bool            m_controlTask;
    bool            m_endEffectorSensing;
    Vec3            m_desiredTaskPosInGround;
};


//==============================================================================
//                    REACHING AND GRAVITY COMPENSATION
//==============================================================================
// This is a task-space controller that tries to move the end effector to 
// a particular target location, and applies gravity compensation and some
// joint damping as lower-priority tasks.
//
// The controller has its own internal Atlas model which in general does not
// match the "real" Atlas perfectly. Each time it is asked to
// generate control torques it reads the sensors on the real Atlas, updates
// its internal model to match. It then generates torques that would be right
// for the internal model, but returns them to be applied to the real Atlas.
class ReachingAndGravityCompensation : public Force::Custom::Implementation {
public:
    // The arm reaches for a target location, and the controller
    // compensates for gravity.
    // 
    // @param[in] taskPointInEndEffector The point whose location we want
    //                                   to control.
    // @param[in] proportionalGain Units of N-m/rad
    // @param[in] derivativeGain Units of N-m-s/rad
    ReachingAndGravityCompensation(const Atlas& realRobot) 
    :   m_modelRobot(), m_modelTasks(m_modelRobot), m_realRobot(realRobot),
        m_targetColor(Red)
    {
        m_modelRobot.initialize(m_modelState);
    }

    const Vec3& getTarget() const {return m_modelTasks.getTarget();}
    Vec3& updTarget() {return m_modelTasks.updTarget();}
    void setTarget(Vec3 pos) {m_modelTasks.setTarget(pos);}

    void toggleGravityComp() {m_modelTasks.toggleGravityComp();}
    void toggleTask() {m_modelTasks.toggleTask();}
    void toggleEndEffectorSensing() {m_modelTasks.toggleEndEffectorSensing();}

    bool isGravityCompensationOn() const 
    {   return m_modelTasks.isGravityCompensationOn(); }

    // This method calculates the needed control torques and adds them into
    // the given mobilityForces Vector which will be applied to the real Atlas.
    // The supplied State is from the real Atlas and will be used to read its
    // sensors.
    void calcForce(const SimTK::State&                realState,
                   SimTK::Vector_<SimTK::SpatialVec>& bodyForces,
                   SimTK::Vector_<SimTK::Vec3>&       particleForces,
                   SimTK::Vector&                     mobilityForces) const
                   OVERRIDE_11;

    // This controller does not contribute potential energy to the system.
    Real calcPotentialEnergy(const SimTK::State& state) const OVERRIDE_11
    { return 0; }

    // Add some useful text and graphics that changes due to user input.
    void calcDecorativeGeometryAndAppend(const State & state, Stage stage,
            Array_<DecorativeGeometry>& geometry) const OVERRIDE_11;

private:
    Atlas                m_modelRobot;   // The controller's internal model.
    TasksMeasure<Vector> m_modelTasks;
    mutable State        m_modelState;   // Temporary: State for the model robot.
    const Atlas&         m_realRobot;    // The "real" robot being controlled.
    const Vec3           m_targetColor;
};


//==============================================================================
//                           USER INPUT HANDLER
//==============================================================================
/// This is a periodic event handler that interrupts the simulation on a
/// regular basis to poll the InputSilo for user input.
class UserInputHandler : public PeriodicEventHandler {
public:
    UserInputHandler(Visualizer::InputSilo&             silo,
                     Atlas&                              realRobot,
                     ReachingAndGravityCompensation&    controller, 
                     Real                               interval)
    :   PeriodicEventHandler(interval), m_silo(silo), m_realRobot(realRobot),
        m_controller(controller), m_increment(0.05) {}

    void handleEvent(State& realState, Real accuracy,
                     bool& shouldTerminate) const OVERRIDE_11;
private:
    Visualizer::InputSilo&          m_silo;
    Atlas&                          m_realRobot;
    ReachingAndGravityCompensation& m_controller;
    const Real                      m_increment;
};


//==============================================================================
//                                  MAIN
//==============================================================================
int main(int argc, char **argv) {
  try {
    cout << "Working dir=" << Pathname::getCurrentWorkingDirectory() << endl;
    Matrix().dump("Ignore this"); // force load for debugging

    // Set some options.
    const double duration = Infinity; // seconds.

    // Create the "real" robot (the one that is being simulated).
    Atlas realRobot;

    // Add the controller.
    ReachingAndGravityCompensation* controller =
        new ReachingAndGravityCompensation(realRobot);
    // Force::Custom takes ownership over controller.
    Force::Custom control(realRobot.updForceSubsystem(), controller);

    // Set up visualizer and event handlers.
    Visualizer viz(realRobot);
    viz.setShowFrameRate(true);
    viz.setShowSimTime(true);

    viz.addSlider("Rate sensor noise", UNoise, 0, 1, 0); 
    viz.addSlider("Angle sensor noise", QNoise, 0, 1, 0); 

    Visualizer::InputSilo* userInput = new Visualizer::InputSilo();
    viz.addInputListener(userInput);
    realRobot.addEventHandler(
            new UserInputHandler(*userInput, realRobot, *controller, 0.05));
    realRobot.addEventReporter(
            new Visualizer::Reporter(viz, 1./30));

    // Display message on the screen about how to start simulation.
    DecorativeText help("Any input to start; ESC to quit.");
    help.setIsScreenText(true);
    viz.addDecoration(MobilizedBodyIndex(0), Vec3(0), help);
    help.setText("Move target: Arrows, PageUp/Down");
    viz.addDecoration(MobilizedBodyIndex(0), Vec3(0), help);

    // Initialize the real robot and other related classes.
    State s;
    realRobot.initialize(s);
    printf("ndofs=%d\n", s.getNU());

    //RungeKuttaMersonIntegrator integ(realRobot);
    SemiExplicitEuler2Integrator integ(realRobot);
    integ.setAccuracy(0.001);
    TimeStepper ts(realRobot, integ);
    ts.initialize(s);
    viz.report(ts.getState());

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



//------------------------------------------------------------------------------
//                TASKS MEASURE :: CALC CACHED VALUE VIRTUAL
//------------------------------------------------------------------------------
// Given a modelState that has been updated from the real robot's sensors, 
// generate control torques as the TasksMeasure's value. This is the only part
// of the code that is actually doing task space operations.

template <class T>
void TasksMeasure<T>::Implementation::calcCachedValueVirtual
   (const State& modelState, int derivOrder, T& tau) const
{
    SimTK_ASSERT1_ALWAYS(derivOrder==0,
        "TasksMeasure::Implementation::calcCachedValueVirtual():"
        " derivOrder %d seen but only 0 allowed.", derivOrder);

    // Shorthands.
    // -----------
    const TaskSpace& p1 = m_tspace1;

    const int nq = modelState.getNQ();
    const int nu = modelState.getNU();
    tau.resize(nu);

    const Real& kd = m_derivativeGain;
    const Real& kp = m_proportionalGain;
    const Vec3& x1_des = m_desiredTaskPosInGround;

    // Abbreviate model state for convenience.
    const State& ms = modelState;

    // Compute control law in task space (F*).
    // ---------------------------------------
    Vec3 xd_des(0);
    Vec3 xdd_des(0);

    // Get info about the actual location, etc. of the model robot.
    Vec3 x1, x1d;
    p1.findStationLocationAndVelocityInGround(ms,
            TaskSpace::StationTaskIndex(0),
            m_modelRobot.getEndEffectorStation(), x1, x1d);

    if (m_endEffectorSensing)
        x1 = m_modelRobot.getSampledEndEffectorPos(ms);


    // Units of acceleration.
    Vec3 Fstar1 = xdd_des + kd * (xd_des - x1d) + kp * (x1_des - x1);

    // Compute task-space force that achieves the task-space control.
    // F = Lambda Fstar + p
    Vector F1 = p1.Lambda(ms) * Fstar1 + p1.mu(ms) + p1.p(ms);
    //Vector F2 = p2.calcInverseDynamics(ms, Fstar2);

    // Combine the reaching task with the gravity compensation and nullspace
    // damping.
    // Gamma = J1T F1 + N1T J1T F2 + N1T N2T (g - c u)
    const Vector& q = ms.getQ();
    const Vector& u = ms.getU();
    const Real k = m_jointPositionGain;
    const Real c = m_jointDampingGain;
    Vector Mu(nu), Mq(nu);
    m_modelRobot.getMatterSubsystem().multiplyByM(ms, u, Mu);
    m_modelRobot.getMatterSubsystem().multiplyByM(ms, q, Mq);

    tau.setToZero();
    //tau +=   p1.JT(s) * F1 
    //              + p1.NT(s) * (  p2.JT(s) * F2 
    //                            + p2.NT(s) * (p1.g(s) - c * u));


    if (m_controlTask) {
        tau += p1.JT(ms) * F1;
        if (m_compensateForGravity)
            tau += p1.NT(ms) * p1.g(ms);
        tau -= p1.NT(ms) * (k*Mq + c*Mu); // damping
    } else if (m_compensateForGravity) {
        tau += p1.g(ms) - (k*Mq + c*Mu);
    } else 
        tau -= k*Mq + c*Mu;

    // Cut tau back to within effort limits.
    // TODO: can't use these limits with one-foot support!
    //const Vector& effortLimits = m_modelRobot.getEffortLimits();
    //for (int i=0; i < nu; ++i) {
    //    const Real oldtau = tau[i], effort = effortLimits[i];
    //    if (std::abs(oldtau) <= effort) continue;
    //    tau[i] = clamp(-effort, oldtau, effort);
    //    printf("Limit tau[%d]: was %g now %g\n", i, oldtau, tau[i]);
    //}

}


//------------------------------------------------------------------------------
//           REACHING AND GRAVITY COMPENSATION :: CALC FORCE
//------------------------------------------------------------------------------
// Given sensor readings from the real robot, generate control torques for it.
// We'll pass on those sensor readings to the task space controller for it to
// use to update its internal modelRobot.
void ReachingAndGravityCompensation::calcForce(
               const State&         realState,
               Vector_<SpatialVec>& bodyForces,
               Vector_<Vec3>&       particleForces,
               Vector&              mobilityForces) const
{
    // Sense the real robot and use readings to update model robot.
    // ------------------------------------------------------------
    const int nq = realState.getNQ(), nu = realState.getNU();
    const Vector& sensedQ = m_realRobot.getSampledAngles(realState);
    const Vector& sensedU = m_realRobot.getSampledRates(realState);
    for (int i=0; i < nq; ++i)
        m_modelRobot.setJointAngle(m_modelState, QIndex(i), sensedQ[i]);
    for (int i=0; i < nu; ++i)
        m_modelRobot.setJointRate(m_modelState, UIndex(i), sensedU[i]);

    // Optional: if real robot end effector location can be sensed, it can
    // be used to improve accuracy. Otherwise, estimate the end effector
    // location using the model robot.
    const Vec3 sensedEEPos = 
        m_realRobot.getSampledEndEffectorPos(realState);
    m_modelRobot.setSampledEndEffectorPos(m_modelState, sensedEEPos);

    m_modelRobot.realize(m_modelState, Stage::Velocity);

    const Vector& tau = m_modelTasks.getValue(m_modelState);

    mobilityForces += tau;
}


//------------------------------------------------------------------------------
//       REACHING AND GRAVITY COMPENSATION :: CALC DECORATIVE GEOMETRY
//------------------------------------------------------------------------------
void ReachingAndGravityCompensation::
calcDecorativeGeometryAndAppend(const State & state, Stage stage,
                                Array_<DecorativeGeometry>& geometry) const
{
    if (stage != Stage::Position) return;

    const Vec3 targetPos = m_modelTasks.getTarget();
    geometry.push_back(DecorativeSphere(0.02)
            .setTransform(targetPos)
            .setColor(m_targetColor));
    geometry.push_back(DecorativeText("Target: " +
        String(targetPos[0])+","+String(targetPos[1])+","+String(targetPos[2]))
        .setIsScreenText(true));

    const MobilizedBody& ee = m_realRobot.getEndEffectorBody();
    Vec3 taskPosInGround = ee.findStationLocationInGround(state,
                                        m_realRobot.getEndEffectorStation());
    geometry.push_back(DecorativePoint(taskPosInGround)
            .setColor(Green).setLineThickness(3));

    geometry.push_back(DecorativeText(String("TOGGLES: [t]Task point ")
        + (m_modelTasks.isTaskPointFollowingOn() ? "ON" : "OFF")
        + "...[g]Gravity comp "
        + (m_modelTasks.isGravityCompensationOn() ? "ON" : "OFF")
        + "...[e]End effector sensor "
        + (m_modelTasks.isEndEffectorSensingOn() ? "ON" : "OFF")
        )
        .setIsScreenText(true));
}


//------------------------------------------------------------------------------
//                USER INPUT HANDLER :: HANDLE EVENT
//------------------------------------------------------------------------------
void UserInputHandler::handleEvent(State& realState, Real accuracy,
                                   bool& shouldTerminate) const
{
    while (m_silo.isAnyUserInput()) {

        int whichSlider; Real sliderValue;
        while (m_silo.takeSliderMove(whichSlider, sliderValue)) {
            if (whichSlider == QNoise) {
                m_realRobot.setAngleNoise(realState, sliderValue);
                continue;
            }
            if (whichSlider == UNoise) {
                m_realRobot.setRateNoise(realState, sliderValue);
                continue;
            }
        }

        unsigned key, modifiers;
        while (m_silo.takeKeyHit(key, modifiers)) {
            if (key == Visualizer::InputListener::KeyEsc) {
                shouldTerminate = true;
                m_silo.clear();
                continue;
            }
            if (key == 'g') {
                m_controller.toggleGravityComp();
                continue;
            }
            if (key == 't') {
                m_controller.toggleTask();
                continue;
            }            
            if (key == 'e') {
                m_controller.toggleEndEffectorSensing();
                continue;
            }
            else if (key == Visualizer::InputListener::KeyRightArrow) {
                // x coordinate goes in and out of the screen.
                m_controller.updTarget()[XAxis] -= m_increment;
                continue;
            }
            else if (key == Visualizer::InputListener::KeyLeftArrow) {
                m_controller.updTarget()[XAxis] += m_increment;
                continue;
            }
            else if (key == Visualizer::InputListener::KeyUpArrow) {
                m_controller.updTarget()[ZAxis] += m_increment;
                continue;
            }
            else if (key == Visualizer::InputListener::KeyDownArrow) {
                m_controller.updTarget()[ZAxis] -= m_increment;
                continue;
            }
            else if (key == Visualizer::InputListener::KeyPageUp) {
                m_controller.updTarget()[YAxis] -= m_increment;
                continue;
            }
            else if (key == Visualizer::InputListener::KeyPageDown) {
                m_controller.updTarget()[YAxis] += m_increment;
                continue;
            }
        }
    }
}