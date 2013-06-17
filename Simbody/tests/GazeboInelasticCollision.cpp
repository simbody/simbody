/* -------------------------------------------------------------------------- *
 *                   Simbody(tm): Gazebo Inelastic Collision                  *
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

/* This test is drawn from the Open Source Robotics Foundation Gazebo physics
regression test "CollisionTest". In the original test there is a cube and 
a sphere, both of unit mass and same half-dimension 0.5. Initially they are 1 
unit apart in x, with the block to the left of the sphere. The block is 
accelerated to the right with a discrete applied force of 1000 applied for 1ms
for a total impulse of 1, which should produce a velocity of +1 that is steady
until the block hits the stationary sphere. The collision is supposed to be 
fully inelastic (coefficient of restitution=0), so there should be no rebound 
and the two objects should move off to the right together at half the impact 
velocity.

The test runs with fixed 1ms steps and calculates the expected result by 
integrating manually. The position and velocity changes during the first step 
with the 1000 unit force active are
   x(.001) = 1/2 a t^2 = 1/2 1e3 1e-6 = .0005 length units
   v(.001) = a t       = 1e3 1e-3     = 1 velocity unit
After that first step, all subsequent steps prior to the collision have
   deltaX = v t       = 1 1e-3       = .001 length units
   deltaV = 0 
We would like the collision to occur exactly at 1s, that is, between steps 1000
and 1001. The position x(1)=.0005 + 999*.001=.9995. So we'll place the sphere
initially so that the block surface and sphere surface are separated by .9995.

Note that while any integrator can get the velocity correct in one step, no 
first order integrator will be able to calculate x(.001) correctly in a single 
step.
  Explicit Euler:           v(.001)=v(0) + h*a(0)=1
                            x(.001)=x(0) + h*v(0)=0                 <--

  Semi-explicit Euler:      v(.001)=v(0) + h*a(0)=1
                            x(.001)=x(0) + h*v(.001)=.001.          <--

  Semi-explicit Euler 2:    v'(.0005)=v(0)+h/2 a(0)=.5
                            x'(.0005)=x(0)+h/2 v'(.0005)=.00025
                            v(.001)=v'(.0005)+h/2 a'(.0005)=1
                            x(.001)=x'(.0005)+h/2 v(.001)=.00075    <--
                       
But any 2nd order integrator would give the correct result, for example 
Explicit trapezoid rule:    v'(.001)=v(0)+h*a(0)=1
                            x(.001)=x(0)+h*(v(0)+v'(.001))/2 =.0005 <--

Explicit midpoint rule:     v'(.0005)=v(0)+h/2*a(0)=1/2
                            x(.001)=x(0)+h*v'(.0005)=.0005          <--
*/

#include "Simbody.h"

#include <cassert>
#include <iostream>
using std::cout; using std::endl;

using namespace SimTK;

//#define USE_VISUALIZER

// Write interesting integrator info to stdout.
static void dumpIntegratorStats(const Integrator& integ);

const Real Radius = 0.5;
const Real Mass = 1;
// Define an extremely stiff, lossy material.
const Real Stiffness = 1e8;
const Real Dissipation = 1000;

const Real MaxStepSize    = Real(1/1000.); // 1 ms (1000 Hz)
const int  DrawEveryN     = 33;            // 33 ms frame update (30.3 Hz)
const Real SimTime        = 2;
const int  NSteps         = // make this a whole number of viz frames
    DrawEveryN*(int(SimTime/MaxStepSize/DrawEveryN+0.5));

const Real TotalImpulse   = 1; // Applied only on the first step
const Real StepForce      = TotalImpulse/MaxStepSize;

const Real IntegAccuracy = 1e-3;
const Real CheckAccuracy = 1e-3;

struct MyMultibodySystem {
    MyMultibodySystem()
    :   m_system(), m_matter(m_system), 
        m_tracker(m_system), m_contact(m_system,m_tracker),
        m_forces(m_system), m_discrete(m_forces,m_matter)
        #ifdef USE_VISUALIZER
        , m_viz(m_system)
        #endif
    {
        ContactMaterial lossyMaterial(Stiffness,
                                      Dissipation,
                                      0,  // mu_static
                                      0,  // mu_dynamic
                                      0); // mu_viscous

        // no gravity
        Body::Rigid sphereBody
           (MassProperties(Mass, Vec3(0), UnitInertia::sphere(Radius)));
        sphereBody.addDecoration(Transform(), DecorativeSphere(Radius));
        sphereBody.addContactSurface(Transform(),
            ContactSurface(ContactGeometry::Sphere(Radius),
                           lossyMaterial));

        // TODO: using inscribed sphere as contact shape within the block.
        Body::Rigid cubeBody
           (MassProperties(Mass, Vec3(0), UnitInertia::sphere(Radius)));
        cubeBody.addDecoration(Transform(), DecorativeBrick(Vec3(Radius)));
        cubeBody.addContactSurface(Transform(),
            ContactSurface(ContactGeometry::Sphere(Radius), // TODO!
                           lossyMaterial));

        MobilizedBody Ground = m_matter.Ground();

        m_cube   = MobilizedBody::Slider(Ground,   Transform(Vec3(0,2,0)), 
                                         cubeBody, Transform(Vec3(0)));
        m_sphere = MobilizedBody::Slider(Ground,   Transform(Vec3(2-.0005,2,0)), 
                                         sphereBody, Transform(Vec3(0)));

        m_system.realizeTopology();

        #ifdef USE_VISUALIZER
        m_viz.setShowFrameNumber(true);
        m_viz.setShowSimTime(true);
        #endif
    }

    MultibodySystem                 m_system;
    SimbodyMatterSubsystem          m_matter;
    ContactTrackerSubsystem         m_tracker;
    CompliantContactSubsystem       m_contact;
    GeneralForceSubsystem           m_forces;
    Force::DiscreteForces           m_discrete;
    #ifdef USE_VISUALIZER
    Visualizer                      m_viz;
    #endif

    MobilizedBody                   m_cube;
    MobilizedBody                   m_sphere;
};

void runOnce(const MyMultibodySystem& mbs, Integrator& integ) 
{
    integ.setAllowInterpolation(false);
    integ.setAccuracy(IntegAccuracy);
    integ.initialize(mbs.m_system.getDefaultState());

    printf("Test with order %d integator %s, Accuracy=%g, MaxStepSize=%g\n",
        integ.getMethodMinOrder(), integ.getMethodName(), 
        integ.getAccuracyInUse(), MaxStepSize); 

    // These variables are the manually calculated values for the cube's
    // x coordinate and x velocity in Ground. We'll calculate these assuming
    // a perfect inelastic collision occurring at t=1 and then compare with
    // the approximate solution produced by the integrator.
    Real x=0, v=0;

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

        if (stepNum == NSteps)
            break;

        ++stepNum;

        mbs.m_discrete.clearAllBodyForces(state);
        if (stepNum == 1)
            mbs.m_discrete.setOneBodyForce(state, mbs.m_cube, 
                SpatialVec(Vec3(0), Vec3(StepForce,0,0)));

        // Advance time by MaxStepSize. Might take multiple internal steps to 
        // get there, depending on difficulty and required accuracy.
        const Real tNext = stepNum * MaxStepSize;
        do {integ.stepTo(tNext,tNext);} while (integ.getTime() < tNext);

        // From Gazebo test code in physics.cc: 
        // integrate here to see when the collision should happen
        if (stepNum == 1) {
            const Real a = StepForce/Mass;      // a is acceleration
            v += a * MaxStepSize;               // dv = a t
            x += a * square(MaxStepSize)/2;     // dx = 1/2 a t^2
            if (integ.getMethodMinOrder() == 1) {
                // 1st order method can't get this right
                const Real prevPos = mbs.m_cube.getOneQ(state,MobilizerQIndex(0));
                mbs.m_cube.setOneQ(state, MobilizerQIndex(0), x);
                printf("-------- FIRST ORDER INTEGRATION METHOD --------\n");
                printf("Cube position repaired: was=%g, now=%g\n", prevPos, x);
                printf("------------------------------------------------\n");
                mbs.m_system.realize(state, Stage::Velocity);
            }
        } else {
            const Real impulse = StepForce*MaxStepSize;
            if (stepNum > 1000)
                v = impulse / (2*Mass);  //inelastic col. w/equal mass
            x += v * MaxStepSize;
        }

        SimTK_TEST_EQ_TOL(mbs.m_cube.getBodyOriginLocation(state)[0], x, CheckAccuracy);
        SimTK_TEST_EQ_TOL(mbs.m_cube.getBodyOriginVelocity(state)[0], v, CheckAccuracy);

        //mbs.m_system.realize(state);
        //printf("after step %d t=%g (h=%g): px=%g,%g vx=%g,%g ax=%g,%g\n",
        //    stepNum, state.getTime(),integ.getPreviousStepSizeTaken(),
        //    mbs.m_cube.getBodyOriginLocation(state)[0],
        //    mbs.m_sphere.getBodyOriginLocation(state)[0],
        //    mbs.m_cube.getBodyOriginVelocity(state)[0],
        //    mbs.m_sphere.getBodyOriginVelocity(state)[0],
        //    mbs.m_cube.getBodyOriginAcceleration(state)[0],
        //    mbs.m_sphere.getBodyOriginAcceleration(state)[0]);
        //printf("  x=%g v=%g\n", x, v);
    }

    dumpIntegratorStats(integ);
}

int main() {
    SimTK_START_TEST("GazeboInelasticCollision");
        // Create the system.   
        MyMultibodySystem mbs;

        RungeKuttaMersonIntegrator rkm(mbs.m_system);
        RungeKutta3Integrator rk3(mbs.m_system);
        RungeKutta2Integrator rk2(mbs.m_system);
        SemiExplicitEuler2Integrator sexpeul2(mbs.m_system);
        ExplicitEulerIntegrator expeul(mbs.m_system);
        RungeKuttaFeldbergIntegrator rkf(mbs.m_system);
        VerletIntegrator verlet(mbs.m_system);

        SimTK_SUBTEST2(runOnce, mbs, rkm);
        SimTK_SUBTEST2(runOnce, mbs, rk3);
        SimTK_SUBTEST2(runOnce, mbs, rk2);
        SimTK_SUBTEST2(runOnce, mbs, sexpeul2);
        SimTK_SUBTEST2(runOnce, mbs, expeul);

        //sherm 130617: these aren't accurate enough
        //SimTK_SUBTEST2(runOnce, mbs, rkf);
        //SimTK_SUBTEST2(runOnce, mbs, verlet);

    SimTK_END_TEST();
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
