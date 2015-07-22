/* -------------------------------------------------------------------------- *
 *             Simbody(tm) - Rattleback using ellipsoid mobilizer             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2009-13 Stanford University and the Authors.        *
 * Authors: Ajay Seth                                                         *
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


/////// TODO TODO TODO ////////
// This test is currently broken because it is using the wrong kind of
// "no slip" constraint. It should be using an ellipsoid-on-halfplane
// rolling constraint but we don't have one yet.
// Thomas Uchida says he'll write one ...
/////// TODO TODO TODO ////////

#define USE_BAD_CONSTRAINTS  // won't conserve energy if these are enabled

/* This program attempts to implement a rattleback (which has an ellipsoid
shape) using Simbody's Ellipsoid mobilizer (3 rotational dofs) rather than
contact constraints. A massless base plate is used to provide two sliding
dofs along the ground plane, for a total of 5 dofs (can't lift off). For
rolling, there must be 2 dofs removed by "no slip" constraint equations.
By construction the contact point between the ellipsoid and ground occurs at
a fixed point in the base plate frame; you can see that in the animation. */

#include "Simbody.h"
#include <iostream>

using namespace SimTK;

// A periodic reporter to dump out the energy and momentum.
// The total energy should be conserved throughout the simulation.
class EnergyReporter : public PeriodicEventReporter {
public:
    EnergyReporter(const MultibodySystem& system,
                   const MobilizedBody& body, Real interval)
    :   PeriodicEventReporter(interval), system(system), body(body) {}

    void handleEvent(const State& state) const {
        system.realize(state, Stage::Dynamics);
        Real energy = system.calcEnergy(state);
        SpatialVec momentum = system.getMatterSubsystem()
                                    .calcSystemMomentumAboutGroundOrigin(state);
        std::cout << state.getTime() << " \tEnergy = " << energy
                  << "  \tAngMom: " << momentum[0].norm()
                  << "  \tLinMom: " << momentum[1].norm() << std::endl;
    }
private:
    const MultibodySystem& system;
    const MobilizedBody&   body;
};


//==============================================================================
//                                  MAIN
//==============================================================================
int main() {
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    matter.setShowDefaultGeometry(false);
    GeneralForceSubsystem forces(system);
    Force::Gravity gravity(forces, matter, -YAxis, 9.80665);
    Visualizer viz(system);

    MobilizedBody& Ground = matter.updGround(); // short name for Ground

    // Calculate the mass properties for a half-ellipsoid. a1,b1,c1 are the
    // radii (semi-axis lengths) of the full ellipsoid in x,y,z resp.
    // Body origin is still center of full ellipsoid (0,0,0) but the COM moves
    // lower in y which affects the xx and zz inertias. Note that Simbody
    // requires the inertias to be given about the body origin, *not* COM.
    // Unit inertia about the body origin is the same as for a full ellipsoid,
    // but will be weighted by half as much mass.
    // TODO: is this right?
    Real m1 = 1.0, a1 = 0.25, b1 = 0.083333333333333, c1 = 0.083333333333333;
    Real comShiftY = (3./8.)*b1; // Because it's a half ellipsoid.
    Body::Rigid halfEllipsoid(MassProperties(m1, Vec3(0, -comShiftY, 0),
                                      UnitInertia::ellipsoid(Vec3(a1,b1,c1))));
    // Add some artwork -- don't have a half-ellipsoid unfortunately.
    halfEllipsoid.addDecoration(Transform(),
        DecorativeEllipsoid(Vec3(a1,b1,c1)).setColor(Red).setResolution(10));

    // Now define a rectangular solid that we'll weld to the rattleback to
    // give it asymmetrical mass properties.
    Real m2 = 2.0, a2 = 2.0*a1, b2 = 0.02, c2 = 0.05;
    const Vec3 barHalfDims = Vec3(a2,b2,c2)/2;
    Body::Rigid barBody(MassProperties(m2, Vec3(0),
                                       UnitInertia::brick(barHalfDims)));
    barBody.addDecoration(Transform(), DecorativeBrick(barHalfDims)
                                                .setColor(Blue).setOpacity(1.));


    // Create a massless x-z base to provide the two slipping dofs.
    MobilizedBody::Slider xdir(Ground, Transform(),
                               Body::Massless(), Transform());
    const Rotation x2z(-Pi/2, YAxis); // rotate so +x moves to +z
    MobilizedBody::Slider base(xdir, x2z, Body::Massless(), x2z);
    base.addBodyDecoration(Transform(), DecorativeBrick(Vec3(0.25, 0.001, 0.25))
                                            .setColor(Orange).setOpacity(0.50));


    // Use a reverse mobilizer so that the contact point remains in a fixed
    // location of the base body.
    MobilizedBody::Ellipsoid rattle(base,          Rotation(Pi/2, XAxis),
                                    halfEllipsoid, Rotation(Pi/2, XAxis),
                                    Vec3(a1,b1,c1), // ellipsoid half-radii
                                    MobilizedBody::Reverse);

    // Weld the bar to the ellipsoid at a 45 degree angle to produce lopsided
    // inertia properties.
    MobilizedBody::Weld bar(rattle, Transform(Rotation(45*Pi/180, YAxis),
                                              Vec3(0,-b2/2.1,0)),
                            barBody, Transform());

    // Finally, the rattle cannot just slide on the surface of the ground, it
    // must roll.
    #ifdef USE_BAD_CONSTRAINTS
    // TODO: (sherm 20130620) these are the wrong constraints because they
    // ignore the acceleration term caused by the contact point moving on the
    // ellipsoid's surface. The correct constraint has to cognizant of the
    // ellipsoid geometry at the contact point. Use of these constraints fails
    // to conserve energy.

    viz.addDecoration(Ground, Vec3(0),
        DecorativeText("TODO: BROKEN -- USING INVALID NOSLIP CONSTRAINTS")
        .setIsScreenText(true));
    Constraint::NoSlip1D contactPointXdir(base, Vec3(0), UnitVec3(1,0,0),
                                          matter.updGround(), rattle);
    Constraint::NoSlip1D contactPointZdir(base, Vec3(0), UnitVec3(0,0,1),
                                          matter.updGround(), rattle);
    #endif

    // Draw a cute green box to rattle around in.
    Ground.addBodyDecoration(Vec3(0,  0*1e-5, 0), // floor
        DecorativeBrick(Vec3(.5, .00001, .5)).setColor(Green).setOpacity(.1));
    Ground.addBodyDecoration(Vec3(0.5, 0.25, 0),  // right wall
        DecorativeBrick(Vec3(1e-5, .25, .5)).setColor(Green).setOpacity(.25));
    Ground.addBodyDecoration(Vec3(-.5, .25, 0),   // left
        DecorativeBrick(Vec3(1e-5, .25, .5)).setColor(Green).setOpacity(.25));
    Ground.addBodyDecoration(Vec3(0, .25, -.5),   // back
        DecorativeBrick(Vec3(.5, .25, 1e-5)).setColor(Green).setOpacity(.25));
    Ground.addBodyDecoration(Vec3(0, .25, .5),    // front
        DecorativeBrick(Vec3(.5, .25, 1e-5)).setColor(Green).setOpacity(.1));

    // Output a visualization frame every 1/30 of a second, and output
    // energy information every 1/4 second.
    system.addEventReporter(new Visualizer::Reporter(viz, 1./30));
    system.addEventReporter(new EnergyReporter(system, rattle, 1./4));

    // We're done building the system. Create it and obtain a copy of the
    // default state.
    State state = system.realizeTopology();

    // Start this off at an angle so it will do something.
    // Caution -- this joint is reversed.
    rattle.setQToFitRotation(state, ~Rotation(10*Pi/180., YAxis));
    //rattle.setUToFitAngularVelocity(state, Vec3(0, -0.5*Pi, 1.0*Pi));

    // Set up simulation.
    RungeKuttaMersonIntegrator integ(system);
    integ.setAccuracy(1e-5);
    TimeStepper ts(system, integ);
    ts.initialize(state);

    // Simulate.
    ts.stepTo(50);
}
