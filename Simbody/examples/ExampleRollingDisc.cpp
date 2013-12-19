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
is of a zero-height (or thickness) disc.  The disc rolls in the x 
direction of the yaw frame/body.
*/

/** Causes the camera to point at and follow a point on a body. This might
 be useful if your system is translating substantially, and would
 otherwise leave the camera's field of view. This class was mostly taken
 from the BodyWatcher class in the TimsBox.cpp test. **/
class BodyFollower : public Visualizer::FrameController {
public:
    /**
     @param bodyB The body to follow.
     @param stationPinB A station (point) P on the body to follow,
        expressed in B.
     @param offset Position of the camera from P, expressed in ground.
     @param upDirection A direction which should point upward as seen by
        the camera.
     **/
    explicit BodyFollower(const MobilizedBody& body,
                          const Vec3&          stationPinB = Vec3(0, 0, 0),
                          const Vec3&          offset = Vec3(1, 1, 1),
                          const Vec3&          upDirection = Vec3(0, 1, 0))
    :   m_body(body), m_stationPinB(stationPinB), m_offset(offset),
        m_upDirection(upDirection) {}
    
    void generateControls(const Visualizer&             viz,
                          const State&                  state,
                          Array_< DecorativeGeometry >& geometry) OVERRIDE_11
    {
        const Vec3 Bp = m_body.findStationLocationInGround(state,
                                                           m_stationPinB);
        const Vec3 camera_position = Bp + m_offset;
        const Rotation camera_R(UnitVec3(m_offset), ZAxis, m_upDirection,
                                viz.getSystemUpDirection().getAxis());
        viz.setCameraTransform(Transform(camera_R, camera_position));
        std::cout << "Total energy: " << m_body.getOneU(state, 0) << std::endl;
    }
    
private:
    const MobilizedBody m_body;
    const Vec3 m_stationPinB;
    const Vec3 m_offset;
    const Vec3 m_upDirection;
};

// Display the time and the system's total energy on the screen.
class TimeAndEnergyDisplayer : public DecorationGenerator {
public:
    
    TimeAndEnergyDisplayer(const MultibodySystem& system)
    : DecorationGenerator(), m_system(system) {}
    
    void generateDecorations(const State&                state,
                             Array_<DecorativeGeometry>& geometry) OVERRIDE_11
    {
        DecorativeText msgT;
        msgT.setIsScreenText(true);
        msgT.setText("Simulation time elapsed: " +
                    String(state.getTime(), "%.3g") + " s");
        geometry.push_back(msgT);
        
        DecorativeText msgE;
        msgE.setIsScreenText(true);
        msgE.setText("Total energy: " +
                    String(m_system.calcEnergy(state), "%.4g") + " J");
        geometry.push_back(msgE);
    }
    
private:
    const MultibodySystem& m_system;
};


int main() {
    
    // Ground x is the direction in which the risc rolls initially.
    // Ground y is perpendicular to the ground.
    // Ground z is the cross product between x and y, and points to the "right"
    // of the disc, if we're sitting behind the disc as it rolls. In the disc's
    // reference frame, z is the normal to the disc.

    // Parameters.
    // -----------
    // The disc is about the size of a coin.
    double radius = 0.5;
    double height = 0.01 * radius;
    double mass = 0.005;

    // We'll use this a lot.
    Vec3 rvec(0, radius, 0);

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
    // The origin of the yawBody is on the ground. Initi
    Body::Rigid yawBody(MassProperties(0.0, Vec3(0), Inertia(0)));
    DecorativeBrick yawDec(Vec3(radius, radius, 0.5 * height));
    yawDec.setColor(Vec3(1.0, 0, 0));
    yawDec.setOpacity(0.5);
    yawBody.addDecoration(Transform(rvec), yawDec);

    // This massless body can fall over, but slips on the ground.
    // The origin of the leanBody is on the ground.
    Body::Rigid leanBody(MassProperties(0.0, Vec3(0), Inertia(0)));
    // This decoration is thin in the direction normal to the disc.
    DecorativeBrick leanDec(Vec3(0.9 * radius, 0.9 * radius, 0.6 * height));
    leanDec.setColor(Vec3(0, 1.0, 0));
    leanDec.setOpacity(0.5);
    leanBody.addDecoration(Transform(rvec), leanDec);

    // Here's the disc, finally.
    // The origin of the discBody is at its centroid.
    Body::Rigid discBody(MassProperties(mass, Vec3(0),
                Inertia::cylinderAlongZ(radius, height)));
    // By default, the axis of DecorativeCylinder is its y axis. We want
    // the axis to be the z axis. The sign could be + or - but we chose - for
    // consistency with the definition of yawMB.
    discBody.addDecoration(Transform(Rotation(-0.5*Pi, XAxis)),
                           DecorativeCylinder(radius, height));

    // Create mobilized bodies.
    // ------------------------
    // Planar MB gives us x and y translation, and rotation about z. We want,
    // in terms of the ground frame, x and z translation, and rotation about y.
    // Using this same transformation for both bodies puts us back at the
    // original frame, which makes the specification of subsequent frames
    // easier.
    MobilizedBody::Planar yawMB(
        matter.updGround(), Transform(Rotation(-0.5*Pi, XAxis)),
        yawBody, Transform(Rotation(-0.5*Pi, XAxis)));
    
    // Pin MB gives us rotation about z. We want rotation about yawMB's x.
    MobilizedBody::Pin leanMB(
        yawMB, Transform(Rotation(0.5*Pi, YAxis)),
        leanBody, Transform(Rotation(0.5*Pi, YAxis)));
    
    // This DOF is the rolling. Positive values for this coordinate indicate
    // rolling forward (in the +x direction of the leanMB frame).
    // We want rotation about leanMB's -z. If we did not do these rotations,
    // positive values of this coordinate would yield movement in the -x
    // direction.
    MobilizedBody::Pin discMB(
        leanMB, Transform(Rotation(Pi, YAxis), rvec),
        discBody, Transform(Rotation(Pi, YAxis)));

    // Constrain the disc to roll.
    // ---------------------------
    // Along the direction of travel.
    // The contact point is fixed in the leanMB frame,
    // at the origin of the leanMB frame. In the leanMB frame, the direction of
    // no-slip is the forward (+x) direction. The two bodies in no-slip contact
    // are the ground and the disc.
    Constraint::NoSlip1D(leanMB, Vec3(0), UnitVec3(1, 0, 0),
            matter.updGround(), discMB);
    
    // Perpendicular to the direction of travel.
    Constraint::NoSlip1D(leanMB, Vec3(0), UnitVec3(0, 0, 1),
            matter.updGround(), discMB);

    // We want the disc to slow down so that it falls over eventually.
    // TODO We can uncomment the following line when energy is conserved otherwise.
    // TODO Force::MobilityLinearDamper(forces, discMB, MobilizerUIndex(0), 0.02);

    // Visualize.
    // ----------
    Visualizer viz(system);
    viz.addFrameController(new BodyFollower(discMB, Vec3(0),
                                            Vec3(-5, 1, 2)));
    viz.addDecorationGenerator(new TimeAndEnergyDisplayer(system));
    
    // Report at the framerate (real-time).
    system.addEventReporter(new Visualizer::Reporter(viz, 1./30));

    // Initialize the system and state.
    // --------------------------------
    State state = system.realizeTopology();

    // We want the disc to start off spinning, in the +x direction.
    discMB.setRate(state, 2*Pi);
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
