/* -------------------------------------------------------------------------- *
 *             Simbody(tm) Example: ClosedTopologyMechanism                   *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2012 Stanford University and the Authors.           *
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
using namespace SimTK;

/* This example shows two different ways to create a closed-topology mechanism:
by cutting a joint and by cutting a body. Generally it is easier to cut a body 
because then all joints are treated the same way, and the only constraint 
needed is a Weld constraint to glue the broken body back together.

The mechanism here is a spatial crank/rocker with three bodies. We're using
a ground frame orientation in which Y is up, X is right, and Z is pointing 
out of the screen at the viewer. The crank is a disk we'll call "rotor", with 
its axis connected by a Pin mobilizer (1 dof revolute internal coordinate 
joint) to ground along the X axis. The rocker is long thin rectangle pinned
at the top to ground along the Z axis, and initially hanging down in the -Y 
direction. Finally there is a "linker" body connected to a point on the 
surface of the rotor by a Ball mobilizer (3 dof spherical internal coordinate 
joint). The linker will connect to the bottom of the rocker by another ball
joint, which gives the mechanism a closed topology. One option is to implement
this final joint directly using a Ball constraint (3 constraint equations). 
Another is to split the linker into two overlapping halves, linkerA and linkerB
with half the mass properties each, connect linkerB to the rocker with a Ball 
mobilizer, and then glue the two halves back together with a Weld constraint 
(6 constraint equations). We'll do it both ways here and then perform a dynamic
simulation (applying a torque to the rotor) and output the final state to 
show that the mechanisms are equivalent.

Here the choice of which approach is better is easy -- the direct approach
gives us a 5 dof tree system and 3 constraints; the second gives an
8 dof tree system and 6 constraints which is much bigger. And the Ball 
constraint is a Simbody built-in and no harder to work with than a Ball 
mobilizer. But in general the split-body approach can be very appealing --
for example, if the loop were closed with a pin joint (a very simple 
mobilizer) the corresponding awkward constraint would have 5 equations, isn't 
a built-in, would have no convenient axis along which to apply forces or 
measure rotations, and special treatment would be required to handle multiple
rotations if they were to be tracked.

(In case you're counting, yes this mechanism has two net degrees of freedom
because the two ball joints permit the linker to rotate about its long axis.)
*/

int main() {
    try { // catch errors if any

    // --------------------------------------------------------------
    // Create the system, with subsystems for bodies and some forces.
    // --------------------------------------------------------------
    MultibodySystem system; 
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::Gravity(forces, matter, -YAxis, 9.8);

    // The rotor will be centered at (0,0,0).
    const Vec3 RockerLocation = Vec3(1.75,1.5,0); // Where to pin rocker.
    const Real TorqueOnRotor = 4;                 // How hard to crank.
    const Vec3 Mech2Offset = Vec3(5,0,0);         // Shift for 2nd mechanism.

    // Some useful rotation matrices we'll need below.
    const Rotation YtoX(-Pi/2, ZAxis); // Reorient to put Y where X was.
    const Rotation ZtoX( Pi/2, YAxis); // Reorient to put Z where X was.

    // --------------------------------------------------------------
    // Define body information. These are just convenient collections
    // of related information useful for constructing the Mobilized
    // Bodies that actually comprise the multibody system.
    // --------------------------------------------------------------  
    // Define body "rotor".
    Real rotorMass = 10, rotorRadius = 1, rotorHalfThickness = .1; 
    Body::Rigid rotorInfo(MassProperties(rotorMass, Vec3(0), 
        UnitInertia::cylinderAlongX(rotorRadius, rotorHalfThickness)));
    rotorInfo.addDecoration(YtoX, 
        DecorativeCylinder(rotorRadius, rotorHalfThickness)); // along Y

    // Define body "rocker".
    Real rockerMass = 3; Vec3 rockerHalfDims(.1, 1, .1);
    Body::Rigid rockerInfo(MassProperties(rockerMass, Vec3(0), 
        UnitInertia::brick(rockerHalfDims)));
    rockerInfo.addDecoration(Vec3(0), 
        DecorativeBrick(rockerHalfDims).setColor(Red));

    // Define a full "linker" for mechanism 1.
    Real linkerMass = 0.5; Vec3 linkerHalfDims(.6, .05, .05);
    Body::Rigid linkerInfo(MassProperties(linkerMass, Vec3(0), 
        UnitInertia::brick(linkerHalfDims)));
    linkerInfo.addDecoration(Vec3(0), 
        DecorativeBrick(linkerHalfDims).setColor(Blue));

    // Define a half "linker" for mechanism 2 (we'll use it twice). The only
    // change is to use half the mass (since UnitInertia gets scaled by mass).
    Body::Rigid halfLinkerInfo(MassProperties(linkerMass/2, Vec3(0),
        UnitInertia::brick(linkerHalfDims)));
    halfLinkerInfo.addDecoration(Vec3(0), 
        DecorativeBrick(linkerHalfDims).setColor(Blue).setOpacity(0.5));

    // --------------------------------------------------------------
    //                 MECHANISM 1 (Ball constraint)
    // --------------------------------------------------------------
    // Note that the Pin mobilizer is defined to rotate about local Z so we
    // need to create local frames with Z pointing along the bodies' X.
    MobilizedBody::Pin rotor1(matter.Ground(),  ZtoX, 
                              rotorInfo,        ZtoX);
    MobilizedBody::Pin rocker1(matter.Ground(), RockerLocation,
                               rockerInfo,      Vec3(0,rockerHalfDims[1],0));
    MobilizedBody::Ball linker1
       (rotor1,     Vec3(rotorHalfThickness,0,.8*rotorRadius),
        linkerInfo, Vec3(-linkerHalfDims[0],0,0));

    // Add a ball constraint instead of a ball mobilizer to connect
    // linker to rocker.
    Constraint::Ball(linker1, Vec3(linkerHalfDims[0],0,0),
                     rocker1, Vec3(0,-rockerHalfDims[1],0));

    // Apply a constant torque about the rotor's Pin axis.
    Force::MobilityConstantForce(forces, rotor1, 0, TorqueOnRotor);

    // --------------------------------------------------------------
    //           MECHANISM 2 (Split linker + Weld constraint)
    // --------------------------------------------------------------  
    MobilizedBody::Pin rotor2(matter.Ground(),  Transform(ZtoX,Mech2Offset), 
                              rotorInfo,        ZtoX);
    MobilizedBody::Pin rocker2(matter.Ground(), RockerLocation+Mech2Offset,
                               rockerInfo,      Vec3(0,rockerHalfDims[1],0));
    // First half-linker connects to the rotor just as above.
    MobilizedBody::Ball linker2a
       (rotor2,         Vec3(rotorHalfThickness,0,.8*rotorRadius),
        halfLinkerInfo, Vec3(-linkerHalfDims[0],0,0));
    // Second half-linker connects to the rocker at the other end.
    MobilizedBody::Ball linker2b
       (rocker2,        Vec3(0,-rockerHalfDims[1],0),
        halfLinkerInfo, Vec3(linkerHalfDims[0],0,0));

    // Now add a weld constraint to glue the half-linkers together. We're
    // taking the default which places the Weld at the body origins.
    Constraint::Weld(linker2a, linker2b);

    // Apply the same constant torque about the rotor's Pin axis.
    Force::MobilityConstantForce(forces, rotor2, 0, TorqueOnRotor);


    // --------------------------------------------------------------
    // SET UP VISUALIZATION, RUN SIMULATION
    // --------------------------------------------------------------  
    // Ask for visualization every 1/30 second.
    system.setUseUniformBackground(true); // turn off floor
    Visualizer viz(system);
    system.adoptEventReporter(new Visualizer::Reporter(viz, 1./30));
    
    // Initialize the system and obtain the default state.    
    State state = system.realizeTopology();

    viz.report(state); // draw frame
    printf("Not yet assembled ... (hit ENTER)\n");
    getchar();

    // Assemble both systems into identical configurations to a tight
    // tolerance. (Note that the linker is free to rotate about its long
    // axis so might not come out exactly the same but it doesn't matter.)
    Assembler assembler(system);
    assembler.lockMobilizer(rotor1);
    assembler.lockMobilizer(rotor2);  
    assembler.setAccuracy(1e-10);
    assembler.assemble(state);

    viz.report(state);
    printf("Assembled ... (hit ENTER)\n");
    getchar();

    // Simulate for 10 seconds.
    RungeKuttaMersonIntegrator integ(system);
    integ.setAccuracy(1e-5); // default is 1e-3
    TimeStepper ts(integ);
    ts.initialize(state);
    ts.stepTo(10);
    state = ts.getState(); // retrieve final state

    // --------------------------------------------------------------
    // OUTPUT FINAL ANGLES AND RATES - should match to accuracy.
    // --------------------------------------------------------------  
    printf("Final angles: rotor1=%g, rocker1=%g\n", 
        rotor1.getAngle(state), rocker1.getAngle(state));
    printf("              roter2=%g, rocker2=%g\n", 
        rotor2.getAngle(state), rocker2.getAngle(state));

    printf("Final rates: rotor1=%g, rocker1=%g\n", 
        rotor1.getRate(state), rocker1.getRate(state));
    printf("             roter2=%g, rocker2=%g\n", 
        rotor2.getRate(state), rocker2.getRate(state));

    } catch (const std::exception& e) {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
