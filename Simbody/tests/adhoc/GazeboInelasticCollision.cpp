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
accelerated to the right with a discrete applied force of 1 until it hits the 
stationary sphere. The collision is supposed to be fully inelastic (coefficient
of restitution=0), so there should be no rebound and the two objects should 
move off to the right together at half the impact velocity.

The test runs with fixed 1ms steps and calculates the expected result by 
integrating manually.

*/

#include "Simbody.h"

#include <cassert>
#include <iostream>
using std::cout; using std::endl;

using namespace SimTK;

// Write interesting integrator info to stdout.
static void dumpIntegratorStats(const Integrator& integ);

int main() {
  try {    
    // Create the system.   
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    ContactTrackerSubsystem tracker(system);
    CompliantContactSubsystem contact(system, tracker);
    GeneralForceSubsystem forces(system);
    Force::DiscreteForces discrete(forces, matter);

    // Define an extremely stiff, lossy material.
    const Real Stiffness = 1e8;
    const Real Dissipation = 1000;
    ContactMaterial lossyMaterial(Stiffness,
                                  Dissipation,
                                  0,  // mu_static
                                  0,  // mu_dynamic
                                  0); // mu_viscous

    // no gravity
    const Real Radius = 0.5;
    const Real Mass = 1;
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

    MobilizedBody Ground = matter.Ground();

    MobilizedBody::Slider cube(Ground,   Transform(Vec3(0,2,0)), 
                               cubeBody, Transform(Vec3(0)));
    MobilizedBody::Slider sphere(Ground,   Transform(Vec3(2-.0005,2,0)), 
                                 sphereBody, Transform(Vec3(0)));

    system.realizeTopology();

    Visualizer viz(system);
    viz.setShowFrameNumber(true);
    viz.setShowSimTime(true);

    // Initialize the system and state.
    
    State state = system.getDefaultState();

    //cube.setRate(state, 1);

    const Real MaxStepSize    = Real(1/1000.); // 1 ms (1000 Hz)
    const int  DrawEveryN     = 33;            // 33 ms frame update (30.3 Hz)
    const Real SimTime        = 2;
    const int  NSteps         = 33*(int(SimTime/MaxStepSize/DrawEveryN+0.5));

    const Real TotalImpulse   = 1; // Applied only on the first step
    const Real StepForce      = TotalImpulse/MaxStepSize;

    RungeKutta3Integrator integ(system);
    //RungeKutta2Integrator integ(system);
    integ.setAccuracy(1e-4);
    //SemiExplicitEulerIntegrator integ(system, MaxStepSize);
    //SemiExplicitEuler2Integrator integ(system);
    //integ.setAccuracy(Real(1e-1)); // 10%
    integ.setAllowInterpolation(false);
    integ.initialize(state);

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

        // Output a frame to the Visualizer if it is time.
        if (stepNum % DrawEveryN == 0)
            viz.report(state);

        if (stepNum == NSteps)
            break;

        ++stepNum;

        discrete.clearAllBodyForces(state);
        if (stepNum == 1)
            discrete.setOneBodyForce(state, cube, 
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
        } else {
            const Real impulse = StepForce*MaxStepSize;
            if (stepNum > 1000)
                v = impulse / (2*Mass);  //inelastic col. w/equal mass
            x += v * MaxStepSize;
        }

        system.realize(state);
        printf("after step %d t=%g (h=%g): px=%g,%g vx=%g,%g ax=%g,%g\n",
            stepNum, state.getTime(),integ.getPreviousStepSizeTaken(),
            cube.getBodyOriginLocation(state)[0],
            sphere.getBodyOriginLocation(state)[0],
            cube.getBodyOriginVelocity(state)[0],
            sphere.getBodyOriginVelocity(state)[0],
            cube.getBodyOriginAcceleration(state)[0],
            sphere.getBodyOriginAcceleration(state)[0]);
        printf("  x=%g v=%g\n", x, v);
    }

    dumpIntegratorStats(integ);

  } catch (const std::exception& e) {
    cout << "EXCEPTION: " << e.what() << "\n";
  }
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
