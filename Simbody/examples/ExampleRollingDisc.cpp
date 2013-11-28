/* -------------------------------------------------------------------------- *
 *                 Simbody(tm) Example: Rolling Disc                          *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
 * Authors: Chris Dembia                                                      *
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

/*
A rolling disc has 3 degrees of freedom. We give the system 5 degrees of
freedom through mobilizers. Then, we remove two of these with no-slip/rolling
non-holonomic constraints. Although the disc has mass-properties of a 3D
cylinder (and is visualized as such), its geometry (for the purpose of contact)
is of a zero-height disc.  The disc rolls in the x direction of the yaw
frame/body.
*/

/// We want the camera to follow the disc.
class FollowDisc : public Visualizer::FrameController
{
public:
    void generateControls(const Visualizer & viz, const State & state,
            Array_<DecorativeGeometry> & geometry)
    {
        Vector q = state.getQ();
        std::cout << q << std::endl;
        // I figured the proper indices by guess and check.
        double x = q[1];
        double z = q[3];
        viz.setCameraTransform(Transform(Vec3(0, 5, x + 5.0)));
        viz.pointCameraAt(Vec3(x, 0, z), Vec3(0, 1, 0));
    }
};

int main() {

    // TODO Explain the degrees of freedom (transforms for mobilized bodies).

    // Parameters.
    // -----------
    // The disc is about the size of a coin.
    double radius = 1.0;
    double height = 0.01 * radius;
    double mass = 0.005;

    // We'll use this a lot.
    Vec3 rvec(0, 0, radius);

    // Create the system.
    // ------------------
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);

    // Add gravity as a force element.
    Force::UniformGravity gravity(forces, matter, Vec3(0, Real(-9.8), 0));

    // Create bodies.
    // --------------
    // This massless body is like a skidding car wheel that can turn/yaw (but
    // not fall over).
    Body::Rigid yawBody(MassProperties(0.0, Vec3(0), Inertia(0)));
    DecorativeBrick yawDec(Vec3(radius, 0.5 * height, radius));
    yawDec.setColor(Vec3(1.0, 0, 0));
    yawDec.setOpacity(0.5);
    yawBody.addDecoration(Transform(rvec), yawDec);

    // This massless body can fall over, but slips on the ground.
    Body::Rigid leanBody(MassProperties(0.0, Vec3(0), Inertia(0)));
    DecorativeBrick leanDec(Vec3(0.9 * radius, 0.6 * height, 0.9 * radius));
    leanDec.setColor(Vec3(0, 1.0, 0));
    leanDec.setOpacity(0.5);
    leanBody.addDecoration(Transform(rvec), leanDec);

    // Here's the disc, finally.
    Body::Rigid discBody(MassProperties(mass, Vec3(0),
                Inertia::cylinderAlongY(radius, height)));
    discBody.addDecoration(Transform(), DecorativeCylinder(radius, height));

    // Create mobilized bodies.
    // ------------------------
    MobilizedBody::Planar yawMB(
        matter.updGround(), Transform(Rotation(-0.5*Pi, XAxis)),
        yawBody, Transform());
    MobilizedBody::Pin leanMB(
        yawMB, Transform(Rotation(0.5*Pi, YAxis)),
        leanBody, Transform(Rotation(0.5*Pi, YAxis)));
    MobilizedBody::Pin discMB(
        leanMB, Transform(Rotation(0.5*Pi, XAxis), rvec),
        discBody, Transform(Rotation(0.5*Pi, XAxis)));

    // Constrain the disc to roll.
    // ---------------------------
    // Along the direction of travel.
    // The contact point is fixed in the leanMB frame.
    Constraint::NoSlip1D(leanMB, Vec3(0), UnitVec3(1, 0, 0),
            matter.updGround(), discMB);
    // Perpendicular to the direction of travel.
    Constraint::NoSlip1D(leanMB, Vec3(0), UnitVec3(0, 1, 0),
            matter.updGround(), discMB);

    // We want the disc to slow down so that it falls over.
    Force::MobilityLinearDamper(forces, discMB, MobilizerUIndex(0), 0.2);

    // Visualize.
    // ----------
    Visualizer viz(system);
    viz.addFrameController(new FollowDisc());
    // Report at the framerate (real-time).
    system.addEventReporter(new Visualizer::Reporter(viz, 1./30));

    // Initialize the system and state.
    // --------------------------------
    State state = system.realizeTopology();

    // We want the disc to start off spinning, in the +x direction.
    discMB.setRate(state, -2 * Pi);
    // So the disc falls over.
    leanMB.setAngle(state, 0.1*Pi);
    discMB.lock(state, Motion::Velocity);
    system.project(state);
    discMB.unlock(state);

    // Simulate the disc.
    // ------------------
    RungeKuttaMersonIntegrator integ(system);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(60.0);
};
