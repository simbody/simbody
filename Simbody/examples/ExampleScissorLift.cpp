/* -------------------------------------------------------------------------- *
 *             Simbody(tm) Example: SimplePlanarMechanism                     *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
 * Authors: Michael Sherman                                                   *
 * Contributors: Kevin He (Roblox)                                            *
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

#include <cstdio>
using std::cout; using std::endl;

// This builds a scissor lift mechanism as suggested by Kevin He at Roblox.
// You can choose the number of scissor levels.
// This is a planar mechanism. It is rather a worst case for an internal
// coordinate multibody system since it has lots of tree degrees of freedom,
// all but one of which is removed by constraints. Execution time will grow
// as O(n^3) once the number of constraints is large enough so that factoring
// the constraint matrix dominates the execution time.
//
// The mechanism is built to operate in the X-Y plane, with Y vertical and
// X to the right. Joints are all pins in the Z direction.

static const int NumLevels = 10;
static const bool PrescribeRotor = true;

int main() {
    try { // catch errors if any

    // Create the system, with subsystems for the bodies and some forces.
    MultibodySystem system; 
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::Gravity(forces, matter, -YAxis, 9.8);

    // Turn off automatically-generated geometry so we just see what we
    // draw here.
    matter.setShowDefaultGeometry(false);

    // Height off the ground, and half width of the scissor mechanism's base.
    Real Height = 2, BaseHalfWidth = .35;

    // Definition of the rotor body, used so the simulation won't be boring.
    Real rotorMass = .5; 
    Vec3 rotorSz(.25,.05,.025);  // half dimensions of rotor body
    Body::Rigid rotorInfo(MassProperties(rotorMass, Vec3(0),
                                         UnitInertia::brick(rotorSz)));
    rotorInfo.addDecoration(Vec3(0), DecorativeBrick(rotorSz).setColor(Red));

    // Describe a long thin rectangular body, with the long direction in Y.
    Real linkMass = 1; 
    Vec3 linkSz(.1,1,.025); // half dimensions of link body
    Body::Rigid linkInfo(MassProperties(linkMass, Vec3(0), 
                                        UnitInertia::brick(linkSz)));
    linkInfo.addDecoration(Vec3(0), DecorativeBrick(linkSz).setColor(Green));

    // Attach the rotor to Ground off to the right and push into -z a little
    // so it is offset from the link.
    MobilizedBody::Pin rotor(matter.Ground(), Vec3(BaseHalfWidth + 2*rotorSz[0], 
                                                   Height, 0),
                             rotorInfo,       Vec3(rotorSz[0],0, rotorSz[2]));

    // Can let the rotor flop or prescribe it to go at a constant velocity.
    if (PrescribeRotor) 
        Constraint::ConstantSpeed(rotor, 1);

    // Create the two trees of mobilized bodies, reusing the above link 
    // description.
    MobilizedBody::Pin right1
       (rotor,           Vec3(-rotorSz[0], 0, rotorSz[2]), 
        linkInfo,        Vec3(0, -linkSz[1], -linkSz[2]));
    MobilizedBody::Pin left1 
       (matter.Ground(), Vec3(-BaseHalfWidth,Height,2*linkSz[2]), 
        linkInfo,        Vec3(0, -linkSz[1], -linkSz[2]));
    right1.setDefaultAngle(Pi/8);
    left1.setDefaultAngle(-Pi/8);

    MobilizedBody::Pin lastRight = right1, lastLeft = left1;
    Real sign = -1; // alternate initial angles
    for (int i=1; i <= NumLevels; ++i) {
        // Add cross connections between the tree ends using two 1-dof 
        // constraints (that's enough since this is planar).
        Constraint::PointInPlane(lastLeft, YAxis, 0.,  // the plane
                                 lastRight, Vec3(0));  // the point
        Constraint::PointInPlane(lastRight, YAxis, 0., // the plane
                                 lastLeft,  Vec3(0));  // the point
        if (i==NumLevels) break;

        // Add generic link pair.
        lastRight = MobilizedBody::Pin(lastRight, Vec3(0, linkSz[1], 0),
                                       linkInfo, Vec3(0, -linkSz[1], 0));
        lastLeft = MobilizedBody::Pin(lastLeft, Vec3(0, linkSz[1], 0),
                                      linkInfo, Vec3(0, -linkSz[1], 0));

        // Set default initial conditions in a zig-zag pattern.
        lastRight.setDefaultAngle(sign*2*Pi/8);
        lastLeft.setDefaultAngle(-sign*2*Pi/8);
        sign = -sign;
    }

    // Ask for visualization every 1/30 second.
    Visualizer viz(system);
    system.addEventReporter(new Visualizer::Reporter(viz, 1./30));
    
    // Initialize the system and state.    
    State state = system.realizeTopology();

    cout << "Scissors before assembly ... ENTER to assemble.\n";
    viz.report(state);
    getchar();

    Assembler(system).assemble(state);

    cout << "Scissors after assembly ... ENTER to run.\n";
    viz.report(state);
    getchar();

    // Run a simulation. There are a variety of integration settings you
    // can play with here. If you put a cap on the max step size, you should
    // use a low order integrator to avoid wasting cycles.

    const Real Accuracy = 0.01; // i.e., 1%. Default is 0.1%.

    const Real MaxStepSize = Infinity;
    //RungeKuttaMersonIntegrator integ(system); // 4th order 

    //const Real MaxStepSize = 1./30;
    RungeKutta3Integrator integ(system);        // 3rd order

    //const Real MaxStepSize = .001;
    //ExplicitEulerIntegrator integ(system);    // 1st order

    integ.setMaximumStepSize(MaxStepSize);
    integ.setAccuracy(Accuracy);

    // Maintain 1mm tolerance even at very loose integration accuracy.
    integ.setConstraintTolerance(std::min(.001, Accuracy/10));

    TimeStepper ts(system, integ);
    ts.initialize(state);

    const double startReal = realTime();
    const double startCPU = cpuTime();
    ts.stepTo(10); // Run simulation for 10s.

    // Finished simulating; dump out some stats. On Windows CPUtime is not
    // reliable if very little time is spent in this thread.
    const double timeInSec = realTime()-startReal;
    const double cpuInSec = cpuTime()-startCPU;
    const int evals = integ.getNumRealizations();
    cout << "Done -- took " << integ.getNumStepsTaken() << " steps in " <<
        timeInSec << "s for " << integ.getTime() << "s sim (avg step=" 
        << (1000*integ.getTime())/integ.getNumStepsTaken() << "ms) " 
        << (1000*integ.getTime())/evals << "ms/eval\n";
    cout << "CPUtime (not reliable when visualizing) " << cpuInSec << endl;

    printf("Used Integrator %s at accuracy %g:\n", 
        integ.getMethodName(), integ.getAccuracyInUse());
    printf("# STEPS/ATTEMPTS = %d/%d\n", integ.getNumStepsTaken(), 
        integ.getNumStepsAttempted());
    printf("# ERR TEST FAILS = %d\n", integ.getNumErrorTestFailures());
    printf("# REALIZE/PROJECT = %d/%d\n", integ.getNumRealizations(), 
        integ.getNumProjections());

    } catch (const std::exception& e) {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
