/* -------------------------------------------------------------------------- *
 *                            Simbody(tm) Example                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-10 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "SimTKsimbody.h"

using namespace SimTK;

int main() {
    
    // Create the system.
    
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);

    // Create bodies.
    Body::Rigid gearBody(MassProperties(1.0, Vec3(0), Inertia(1)));
    gearBody.addDecoration(Transform(Rotation(0.5*Pi, XAxis)), 
        DecorativeCylinder(1.0, 0.1));
    Body::Rigid rodBody(MassProperties(1.0, Vec3(0), Inertia(1)));
    rodBody.addDecoration(Transform(Vec3(0, 1, 0)), DecorativeCylinder(0.05, 1.0));

    // Create instances of the bodies that are connected into the multibody
    // system via Mobilizers. Note that we use the gear body twice.

    MobilizedBody::Pin gear1(matter.updGround(), Transform(Vec3(1, 0, 0)), 
                             gearBody, Transform());
    MobilizedBody::Pin gear2(matter.updGround(), Transform(Vec3(-1, 0, 0)), 
                             gearBody, Transform());
    MobilizedBody::Pin rod(gear2, Transform(Vec3(0, 0.8, 0.1)), rodBody, Transform());

    // Add constraints.
    Constraint::ConstantSpeed(gear1, 2*Pi); // i.e., 1 rotation per second
    Constraint::NoSlip1D(matter.updGround(), Vec3(0), UnitVec3(0, 1, 0), gear1, gear2);
    
    // We want the rod end point traveling along a line. We'll draw part of the
    // line to make it clear.
    Constraint::PointOnLine(matter.updGround(), UnitVec3(0, 1, 0), Vec3(0, 0, 0.1), 
                            rod, Vec3(0, 2, 0));
    matter.updGround().addBodyDecoration(Vec3(0), // transform ignored for line
        DecorativeLine(Vec3(0,0,.1), Vec3(0,0,.1)+3*UnitVec3(0,1,0))
        .setColor(Red));
   
    // Visualize the system, reporting an output frame every 1/30 of a simulated
    // second. The Visualizer's default frame rate is 30fps, and it will slow the
    // simulation down to keep to that speed, so we'll get exactly real time 
    // this way. We don't want the default ground and sky background here.
    Visualizer viz(system);
    viz.setBackgroundType(Visualizer::SolidColor); // default is white
    system.updDefaultSubsystem().addEventReporter(new Visualizer::Reporter(viz, 1./30));
    
    // Initialize the system and state.
    
    system.realizeTopology();
    State state = system.getDefaultState();
    
    // Simulate it.

    RungeKuttaMersonIntegrator integ(system);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(1000.0);
}
