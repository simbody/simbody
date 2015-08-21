/* -------------------------------------------------------------------------- *
 *    Simbody(tm): Gazebo Basic Controller Response                           *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2015 Stanford University and the Authors.           *
 * Authors: Michael Sherman, John Hsu                                         *
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

/* This test is drawn from the Open Source Robotics Foundation Gazebo physics
regression test "PhysicsTest::TrikeyWheelResponse". From the Gazebo test:

    The trikey model has three wheels oriented in different directions.
    it caught a corner case in ODE where the inertia was being
    truncated unequally based on cartesian orientation of the link.
    Hence this test is added to ensure we get the same dynamics behavior
    regardless the orientation of the inertia matrices in world frame.

There is a controller on each wheel ("PID" but actually has only a P gain)
that is instructed to rotate each wheel to an arbitrary target position.

The test just requires that all three wheels behave identically; it doesn't
check whether they also behave correctly!
*/

#include "Simbody.h"

#include <cassert>
#include <iostream>
using std::cout; using std::endl;

using namespace SimTK;

//#define USE_VISUALIZER

// Otherwise use a discrete controller that is spring-like. Energy is 
// conserved only if you use the spring because the controller doesn't
// store any.
//#define USE_SPRING


// Control gain
const Real Kp1 = 10; // proportional to angle error

// Target angle
const Real Target1 = 1.3; // radians

const Real TriKeyRadius = 0.27;
const Real TriKeyHalfLength = 0.65/2;
const Real TriKeyMass = 32.7;
const Vec3 TriKeyCOM(-0.003048,0.00254,0.41415);
const Inertia TriKeyCentralInertia(1.747,1.747,1.192); // Central
const Inertia TriKeyInertia(TriKeyCentralInertia
                            .shiftFromMassCenter(-TriKeyCOM, TriKeyMass));
const Vec3 TriKeyColor = Gray;

const Real WheelRadius = 0.101;
const Real WheelHalfLength = 0.025/2;
const Real WheelMass = 0.66725; 
const Vec3 WheelCOM(0,0,0);
const Inertia WheelInertia(0.00160418,0.00279814,0.00160339,   //xx,yy,zz
                           -3.598e-06,1.6199e-05,-4.1656e-05); //xy,xz,yz
const Vec3 WheelColors[] = {Red, Green, Blue};

const Real Accuracy       = 1e-3;
const Real MaxStepSize    = Real(1/1000.); // 1 ms (1000 Hz)
const int  DrawEveryN     = 5;
const Real SimTime        = 1.5;
const int  NSteps         = // make this a whole number of viz frames
    DrawEveryN*(int(SimTime/MaxStepSize/DrawEveryN+0.5));

// Use this class to hold references into the Simbody system.
struct MyMultibodySystem {
    MyMultibodySystem(); // see below

    MultibodySystem                 m_system;
    SimbodyMatterSubsystem          m_matter;
    GeneralForceSubsystem           m_forces;
    Force::DiscreteForces           m_discrete;
    #ifdef USE_VISUALIZER
    Visualizer                      m_viz;
    #endif

    MobilizedBody::Pin              m_trikey_base;
    MobilizedBody::Pin              m_wheel[3];
};

// Execute the simulation with a given integrator and accuracy, and verify that
// it produces the correct answers.
static void runOnce(const MyMultibodySystem& mbs, Integrator& integ, 
                    Real accuracy);


//==============================================================================
//                                   MAIN
//==============================================================================
int main() {
    SimTK_START_TEST("TrikeyWheelResponse");
        // Create the system.   
        MyMultibodySystem mbs;

        printf("SemiExplicitEuler2 @ %g accuracy\n", Accuracy);
        SemiExplicitEuler2Integrator sexpeul2(mbs.m_system);
        SimTK_SUBTEST3(runOnce, mbs, sexpeul2,Accuracy);

        printf("RungeKuttaMerson @ %g accuracy\n",Accuracy/1000);
        RungeKuttaMersonIntegrator rkm(mbs.m_system);
        SimTK_SUBTEST3(runOnce,mbs,rkm,Accuracy/1000);

    SimTK_END_TEST();
}


//==============================================================================
//                           MY MULTIBODY SYSTEM
//==============================================================================
// Construct the multibody system. The dampers are built in here but the springs
// are applied during execution.
MyMultibodySystem::MyMultibodySystem()
:   m_system(), m_matter(m_system), 
    m_forces(m_system), m_discrete(m_forces,m_matter)
    #ifdef USE_VISUALIZER
    , m_viz(m_system)
    #endif
{
    #ifdef USE_VISUALIZER
    m_viz.setSystemUpDirection(ZAxis);
    m_viz.setShowFrameNumber(true);
    m_viz.setShowSimTime(true);
    #endif

    Force::Gravity(m_forces, m_matter, -ZAxis, 9.81);

    // Cylinder is along Z in Gazebo, Y in Simbody
    Rotation YtoZ(Pi/2,XAxis);
    DecorativeCylinder drawTriKey(TriKeyRadius, TriKeyHalfLength);
    drawTriKey.setOpacity(0.5).setTransform(YtoZ);

    DecorativeCylinder drawWheel(WheelRadius,WheelHalfLength);
    drawWheel.setOpacity(0.5).setTransform(YtoZ);

    Body::Rigid triKeyInfo(MassProperties(TriKeyMass, Vec3(TriKeyCOM), TriKeyInertia));
    Body::Rigid wheelInfo(MassProperties(WheelMass, Vec3(WheelCOM), WheelInertia));
  
    MobilizedBody& Ground = m_matter.updGround(); // Nicer name for Ground.

    m_trikey_base = MobilizedBody::Pin(Ground,Vec3(0,0,.5), 
                                       triKeyInfo, Vec3(0));
    m_trikey_base.addBodyDecoration(Vec3(0,0,.426), 
                                    drawTriKey.setColor(TriKeyColor));
    m_trikey_base.addBodyDecoration(Vec3(0), DecorativeFrame(0.5).setColor(White));

    const Rotation ZtoMinusY(Pi/2, XAxis);
    for (int i=0; i < 3; ++i) {
        Real offs[] = {0., 1e-1, -1e-1};
        Real angle = Real(i)*2*Pi/3 + offs[i]; // 0, 120, 240
        Rotation aboutZ(angle, ZAxis);
        Transform X_IF(aboutZ*ZtoMinusY, aboutZ*Vec3(0,-.24,.1)+Vec3(offs[i],0,0));
        Transform X_OM(YtoZ, Vec3(0));
        std::cout << std::setprecision(15) << "  X_IF=" << X_IF.p() << X_IF.R();
        std::cout << std::setprecision(15) << "  X_OM=" << X_OM.p() << X_OM.R();

        m_wheel[i] = MobilizedBody::Pin(m_trikey_base, X_IF,
                                        wheelInfo, X_OM);

        m_wheel[i].addBodyDecoration(Transform(Rotation(Pi/2,XAxis),
                                               Vec3(0,-0.05,2.44839e-13)),
                                     drawWheel.setColor(WheelColors[i]));
        m_wheel[i].addBodyDecoration(Transform(),
                      DecorativeFrame(.2).setColor(Orange));

        #ifdef USE_SPRING
        // Using a spring is more accurate than the discrete spring-like controller.
        Force::MobilityLinearSpring
           (m_forces,m_wheel[i], MobilizerQIndex(0), Kp1, Target1);
        #endif
    }


    m_system.realizeTopology();

}

//==============================================================================
//                                 RUN ONCE
//==============================================================================

// Write interesting integrator info to stdout.
static void dumpIntegratorStats(const Integrator& integ);

Real calcTotalEnergy(const MyMultibodySystem& mbs, const State& state) {
    mbs.m_system.realize(state, Stage::Dynamics);
    // Calculate potential energy in controller.
    Real controllerPE = 0;
    #ifndef USE_SPRING
        for(int i=0; i < 3; ++i) {
            const Real aerr = mbs.m_wheel[i].getAngle(state)-Target1;
            controllerPE += Kp1*square(aerr)/2; // 1/2 k x^2
        }
    #endif
    return mbs.m_system.calcEnergy(state) + controllerPE;
}

// Run the system until it settles down, then check the answers.
void runOnce(const MyMultibodySystem& mbs, Integrator& integ, Real accuracy) 
{
    integ.setAllowInterpolation(false);
    integ.setAccuracy(accuracy);
    integ.initialize(mbs.m_system.getDefaultState());

    printf("Test with order %d integator %s, Accuracy=%g, MaxStepSize=%g\n",
        integ.getMethodMinOrder(), integ.getMethodName(), 
        integ.getAccuracyInUse(), MaxStepSize); 

    #ifdef USE_VISUALIZER
    mbs.m_viz.report(integ.getState());
    printf("Hit ENTER to simulate:");
    getchar();
    #endif

    //mbs.m_trikey_base.setOneU(integ.updAdvancedState(), MobilizerUIndex(0), 1);
    //mbs.m_trikey_base.lock(integ.updAdvancedState());

    printf("\n--------------------------------------------------------\n");
    #ifdef USE_SPRING
        printf("Using spring: energy should be conserved.\n");
    #else
        printf("Using controller: don't expect energy to be conserved.\n");
    #endif
    printf("--------------------------------------------------------\n");

    printf("Energy = %.15g\n", calcTotalEnergy(mbs, integ.getState()));

    unsigned stepNum = 0;
    while (true) {
        // Get access to State being advanced by the integrator. Interpolation 
        // must be off so that we're modifying the actual trajectory.
        State& state = integ.updAdvancedState();

        #ifdef USE_VISUALIZER
        // Output a frame to the Visualizer if it is time.
        if (stepNum % DrawEveryN == 0)
            mbs.m_viz.report(state);
        #endif

        if (stepNum++ == NSteps)
            break;

        #ifndef USE_SPRING
            // Apply controller forces.
            for (int i=0; i < 3; ++i) {
                const Real aerr = mbs.m_wheel[i].getAngle(state)-Target1;
                mbs.m_discrete.setOneMobilityForce(state, mbs.m_wheel[i],
                                                   MobilizerUIndex(0), -Kp1*aerr);
            }
        #endif
        //printf("Energy = %.15g\n",calcTotalEnergy(mbs,integ.getState()));

        // Advance time by MaxStepSize. Might take multiple internal steps to 
        // get there, depending on difficulty and required accuracy.
        const Real tNext = stepNum * MaxStepSize;
        do {integ.stepTo(tNext,tNext);} while (integ.getTime() < tNext);
    }

    const State& state = integ.getAdvancedState();
    mbs.m_system.realize(state);

    printf("Energy = %.15g\n",calcTotalEnergy(mbs,state));

    printf("TriKey angle=%.15g rate=%.15g\n",
           mbs.m_trikey_base.getAngle(state),mbs.m_trikey_base.getRate(state));

    for (int i=0; i < 3; ++i) {
        const Real angle = mbs.m_wheel[i].getAngle(state);
        const Real rate = mbs.m_wheel[i].getRate(state);
        printf("%.15g %.15g\n", angle, rate);
    }
    printf("\n");

    // These should be very similar since they are all treated the same.
    // They might not be right, but they should match!
    SimTK_TEST_EQ_TOL(mbs.m_wheel[0].getAngle(state),
                      mbs.m_wheel[1].getAngle(state), 1e-12);
    SimTK_TEST_EQ_TOL(mbs.m_wheel[0].getAngle(state),
                      mbs.m_wheel[2].getAngle(state), 1e-12);
    SimTK_TEST_EQ_TOL(mbs.m_wheel[0].getRate(state),
                      mbs.m_wheel[1].getRate(state), 1e-12);
    SimTK_TEST_EQ_TOL(mbs.m_wheel[0].getRate(state),
                      mbs.m_wheel[2].getRate(state), 1e-12);


    dumpIntegratorStats(integ);
}



//==============================================================================
//                        DUMP INTEGRATOR STATS
//==============================================================================
static void dumpIntegratorStats(const Integrator& integ) {
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
