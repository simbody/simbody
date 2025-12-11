/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
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
using namespace std;

void testCalculationMethods() {
    
    // Create a system with two bodies.
    
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    MobilizedBody::Free b1(matter.Ground(), body);
    MobilizedBody::Free b2(matter.Ground(), body);
    
    // Set all the state variables to random values.

    system.realizeTopology();
    State state = system.getDefaultState();
    Random::Gaussian random;

    for (int i = 0; i < state.getNY(); ++i)
        state.updY()[i] = random.getValue();

    system.realize(state, Stage::Acceleration);
   
    // Test the low level methods for transforming points and vectors.
    
    const Vec3 point(0.5, 1, -1.5);
    SimTK_TEST_EQ(b1.findStationLocationInGround(state, Vec3(0)), b1.getBodyOriginLocation(state));
    SimTK_TEST_EQ(b1.findStationAtGroundPoint(state, b1.findStationLocationInGround(state, point)), point);
    SimTK_TEST_EQ(b2.findStationAtGroundPoint(state, b1.findStationLocationInGround(state, point)), b1.findStationLocationInAnotherBody(state, point, b2));
    SimTK_TEST_EQ(b2.findStationAtGroundPoint(state, b1.findStationLocationInGround(state, Vec3(0))).norm(), (b1.getBodyOriginLocation(state)-b2.getBodyOriginLocation(state)).norm());
    SimTK_TEST_EQ(b2.findMassCenterLocationInGround(state), b2.findStationLocationInGround(state, b2.getBodyMassCenterStation(state)));
    SimTK_TEST_EQ(b1.expressVectorInGroundFrame(state, Vec3(0)), Vec3(0));
    SimTK_TEST_EQ(b1.expressVectorInGroundFrame(state, point), b1.getBodyRotation(state)*point);
    SimTK_TEST_EQ(b1.expressGroundVectorInBodyFrame(state, b1.expressVectorInGroundFrame(state, point)), point);
    SimTK_TEST_EQ(b2.expressGroundVectorInBodyFrame(state, b1.expressVectorInGroundFrame(state, point)), b1.expressVectorInAnotherBodyFrame(state, point, b2));
    
    // Test the routines for mapping locations, velocities, and accelerations.
    
    Vec3 r, v, a;
    b1.findStationLocationVelocityAndAccelerationInGround(state, point, r, v, a);
    SimTK_TEST_EQ(v, b1.findStationVelocityInGround(state, point));
    SimTK_TEST_EQ(a, b1.findStationAccelerationInGround(state, point));
    {
        Vec3 r2, v2;
        b1.findStationLocationAndVelocityInGround(state, point, r2, v2);
        SimTK_TEST_EQ(r, r2);
        SimTK_TEST_EQ(v, v2);
    }
    SimTK_TEST_EQ(b1.findStationVelocityInGround(state, Vec3(0)), b1.getBodyOriginVelocity(state));
    SimTK_TEST_EQ(b1.findStationAccelerationInGround(state, Vec3(0)), b1.getBodyOriginAcceleration(state));
    SimTK_TEST_EQ(b1.findStationVelocityInGround(state, point), b1.findStationVelocityInAnotherBody(state, point, matter.Ground()));
}

void testFittingMethods() {

    // Define an empty system.
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));

    // Define zero and arbitrary transforms and velocities.
    Transform X_FM_zero(Vec3(0));
    SpatialVec V_FM_zero(Vec3(0), Vec3(0));
    Transform X_FM(Rotation(-1.9, Vec3(-3, 2, 4)), Vec3(-0.33, 0.66, -0.99));
    SpatialVec V_FM(Vec3(-3, .4, -.7), Vec3(1, .2, .3));

    // MobilizedBody::Free
    {
        MobilizedBody::Free mobod(matter.updGround(), body);
        system.realizeTopology();
        State state = system.getDefaultState();

        // Zero transform and velocity.
        mobod.setQToFitTransform(state, X_FM_zero);
        SimTK_TEST_EQ(mobod.getQ(state), Vec7(1, 0, 0, 0, 0, 0, 0));
        mobod.setUToFitVelocity(state, V_FM_zero);
        SimTK_TEST_EQ(mobod.getU(state), Vec6(0));

        // Arbitrary transform and velocity.
        mobod.setQToFitTransform(state, X_FM);
        Quaternion quat = X_FM.R().convertRotationToQuaternion();
        Vec3 trans = X_FM.p();
        SimTK_TEST_EQ(mobod.getQ(state),
            Vec7(quat[0], quat[1], quat[2], quat[3],
                 trans[0], trans[1], trans[2]));
        mobod.setUToFitVelocity(state, V_FM);
        SimTK_TEST_EQ(mobod.getU(state),
            Vec6(V_FM[0][0], V_FM[0][1], V_FM[0][2],
                 V_FM[1][0], V_FM[1][1], V_FM[1][2]));
    }

    // MobilizedBody::CantileverFreeBeam
    {
        Real length = 1.23;
        MobilizedBody::CantileverFreeBeam mobod(matter.updGround(),
            body, length);
        system.realizeTopology();
        State state = system.getDefaultState();

        // Zero transform and velocity.
        mobod.setQToFitTransform(state, X_FM_zero);
        SimTK_TEST_EQ(mobod.getQ(state), Vec3(0));
        mobod.setUToFitVelocity(state, V_FM_zero);
        SimTK_TEST_EQ(mobod.getU(state), Vec3(0));

        // Check that setting a Z-translation slightly larger or slightly
        // smaller than the beam length produces q = 0.
        mobod.setQToFitTranslation(state, Vec3(0, 0, 0.99 * length));
        SimTK_TEST_EQ(mobod.getQ(state), Vec3(0));
        mobod.setQToFitTranslation(state, Vec3(0, 0, 1.01 * length));
        SimTK_TEST_EQ(mobod.getQ(state), Vec3(0));
        // While we're at it, check that setting the Z-translation to exactly
        // the beam length also produces q = 0.
        mobod.setQToFitTranslation(state, Vec3(0, 0, length));
        SimTK_TEST_EQ(mobod.getQ(state), Vec3(0));

        // Reset the state to the zero transform and velocity.
        mobod.setQToFitTransform(state, X_FM_zero);
        mobod.setUToFitVelocity(state, V_FM_zero);

        // Arbitrary rotation.
        mobod.setQToFitRotation(state, X_FM.R());
        SimTK_TEST_EQ(mobod.getQ(state),
            X_FM.R().convertRotationToBodyFixedXYZ());

        // Arbitrary translation.
        Vec3 q = mobod.getQ(state);
        mobod.setQToFitTranslation(state, X_FM.p());
        // Based on RBNodeCantileverFreeBeam::setQToFitTranslationImpl().
        const UnitVec3 e(X_FM.p());
        const Real latitude = std::atan2(-e[1], e[2]);
        const Real longitude = std::atan2(e[0], e[2]);
        Real spin = q[2]; // q2 value *prior* to fitting the translation
        const Rotation R_FM =
            Rotation(SpaceRotationSequence, latitude, XAxis, longitude, YAxis) *
            Rotation(spin, ZAxis);
        SimTK_TEST_EQ(mobod.getQ(state), R_FM.convertRotationToBodyFixedXYZ());

        // Arbitrary transform.
        mobod.setQToFitTransform(state, X_FM);
        // Should be the same as the translation fitting above.
        SimTK_TEST_EQ(mobod.getQ(state), R_FM.convertRotationToBodyFixedXYZ());

        // Arbitrary angular velocity.
        q = mobod.getQ(state); // update q to the current value in the state
        mobod.setUToFitAngularVelocity(state, V_FM[0]);
        // Based on RBNodeCantileverFreeBeam::setUToFitAngularVelocityImpl().
        const Vec2 cosxy(std::cos(q[0]), std::cos(q[1]));
        const Vec2 sinxy(std::sin(q[0]), std::sin(q[1]));
        const Real oocosy = 1 / cosxy[1];
        Vec3 qdot = Rotation::convertAngVelInParentToBodyXYZDot(
            cosxy, sinxy, oocosy, V_FM[0]);
        SimTK_TEST_EQ(mobod.getU(state), qdot);

        // Arbitrary linear velocity.
        Vec3 u = mobod.getU(state); // get the current u, prior to fitting
        mobod.setUToFitLinearVelocity(state, V_FM[1]);

        // Based on RBNodeCantileverFreeBeam::setUToFitLinearVelocityImpl().
        Matrix m(3, 2, 0.0);
        m(0, 1) = (2.0 / 3.0) * length;
        m(1, 0) = -m(0, 1);
        m(2, 0) = -(8.0 / 15.0) * length * q[0];
        m(2, 1) = -(8.0 / 15.0) * length * q[1];

        // Solve for qdot0 and qdot1 via least squares.
        FactorQTZ qtz(m);
        Vector b(3, 0.0);
        b[0] = V_FM[1][0];
        b[1] = V_FM[1][1];
        b[2] = V_FM[1][2];
        Vector x(2, 0.0);
        qtz.solve(b, x);

        // Should match the qdot0 and qdot1 we solved for above, plus the
        // current value of u2.
        qdot = Vec3(x[0], x[1], u[2]);
        SimTK_TEST_EQ(mobod.getU(state), qdot);

        // Arbitrary velocity.
        mobod.setUToFitVelocity(state, V_FM);
        // Should be the same as the linear velocity fitting above.
        SimTK_TEST_EQ(mobod.getU(state), qdot);
    }
}

void testWeld() {
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -1, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    
    // Create two pendulums, each with two welded bodies.  One uses a Weld MobilizedBody,
    // and the other uses a Weld constraint.
    
    Transform inboard(Vec3(0.1, 0.5, -1));
    Transform outboard(Vec3(0.2, -0.2, 0));
    MobilizedBody::Ball p1(matter.updGround(), Vec3(0), body, Vec3(0, 1, 0));
    MobilizedBody::Ball p2(matter.updGround(), Vec3(0), body, Vec3(0, 1, 0));
    MobilizedBody::Weld c1(p1, inboard, body, outboard);
    MobilizedBody::Free c2(p2, inboard, body, outboard);
    Constraint::Weld constraint(p2, inboard, c2, outboard);

    // It is not a general test unless the Weld mobilizer has children!
    MobilizedBody::Pin wchild1(c1, inboard, body, outboard);
    MobilizedBody::Pin wchild2(c2, inboard, body, outboard);
    Force::MobilityLinearSpring(forces, wchild1, 0, 1000, 0);
    Force::MobilityLinearSpring(forces, wchild2, 0, 1000, 0);

    State state = system.realizeTopology();
    p1.setU(state, Vec3(1, 2, 3));
    p2.setU(state, Vec3(1, 2, 3));
    system.realize(state, Stage::Velocity);
    system.project(state, 1e-10);

    SimTK_TEST_EQ(c1.getBodyTransform(state), c2.getBodyTransform(state));
    SimTK_TEST_EQ(c1.getBodyVelocity(state), c2.getBodyVelocity(state));
    
    // Simulate it and see if both pendulums behave identically.
    
    RungeKuttaMersonIntegrator integ(system);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(5);
    system.realize(integ.getState(), Stage::Velocity);
    SimTK_TEST_EQ_TOL(c1.getBodyTransform(integ.getState()), 
                      c2.getBodyTransform(integ.getState()), 1e-10);
    SimTK_TEST_EQ_TOL(c1.getBodyVelocity(integ.getState()), 
                      c2.getBodyVelocity(integ.getState()), 1e-10);
}


// Create two pendulums with gimbal joints, one where the gimbals are faked
// with a series of pin joints that should provide identical generalized
// coordinates and speeds (we're expecting the new gimbal definition in which
// the generalized speeds are the Euler angle derivatives).

void testGimbal() {
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.81, 0));
    Body::Rigid lumpy(MassProperties(3.1, Vec3(.1, .2, .3), 
                      UnitInertia(1.2,1.1,1.3,.01,.02,.03)));
    Body::Rigid massless(MassProperties(0,Vec3(0),UnitInertia(0)));
    
    Vec3 inboard(0.1, 0.5, -1);
    Vec3 outboard(0.2, -0.2, 0);
    MobilizedBody::Gimbal p1(matter.updGround(), inboard, lumpy, outboard);

    Rotation axisX(Pi/2, YAxis); // rotate z into x
    Rotation axisY(-Pi/2, XAxis); // rotate z into y
    Rotation axisZ; // leave z where it is
    MobilizedBody::Pin dummy1(matter.updGround(), Transform(axisX,inboard),
        massless, Transform(axisX, Vec3(0)));
    MobilizedBody::Pin dummy2(dummy1, Transform(axisY, Vec3(0)),
        massless, Transform(axisY, Vec3(0)));
    MobilizedBody::Pin p2(dummy2, Transform(axisZ, Vec3(0)),
        lumpy, Transform(axisZ, outboard));

    MobilizedBody::Weld c1(p1, inboard, lumpy, outboard);
    MobilizedBody::Weld c2(p2, inboard, lumpy, outboard);

    c1.addBodyDecoration(Transform(), DecorativeBrick(.1*Vec3(1,2,3))
        .setColor(Cyan).setOpacity(0.2));
    c2.addBodyDecoration(Transform(), DecorativeBrick(.1*Vec3(1,2,3))
        .setColor(Red).setRepresentation(DecorativeGeometry::DrawWireframe));

    //Visualizer viz(system); viz.setBackgroundType(Visualizer::SolidColor);
    //system.addEventReporter(new Visualizer::Reporter(viz, 1./30));

    State state = system.realizeTopology();
    p1.setQ(state, Vec3(.1,.2,.3));
    dummy1.setAngle(state, .1);
    dummy2.setAngle(state, .2);
    p2.setAngle(state, .3);

    p1.setU(state, Vec3(1, 2, 3));
    dummy1.setRate(state, 1);
    dummy2.setRate(state, 2);
    p2.setRate(state, 3);

    system.realize(state, Stage::Velocity);
    system.project(state, 1e-10);

    SimTK_TEST_EQ(c1.getBodyTransform(state), c2.getBodyTransform(state));
    SimTK_TEST_EQ(c1.getBodyVelocity(state), c2.getBodyVelocity(state));
    
    // Simulate it and see if both pendulums behave identically.
    
    RungeKuttaMersonIntegrator integ(system);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(5);
    system.realize(integ.getState(), Stage::Velocity);
    SimTK_TEST_EQ_TOL(c1.getBodyTransform(integ.getState()), 
                  c2.getBodyTransform(integ.getState()), 1e-10);
    SimTK_TEST_EQ_TOL(c1.getBodyVelocity(integ.getState()), 
                      c2.getBodyVelocity(integ.getState()), 1e-10);
}


// Create two pendulums with bushing joints, one where the bushings are faked
// with a Cartesian (3 sliders) and a series of pin joints that should provide 
// identical generalized coordinates and speeds, except that the translational
// coordinates are the last three q's for the bushing but are interpreted as
// translating first, then rotating.

void testBushing() {
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.81, 0));
    Body::Rigid lumpy(MassProperties(3.1, Vec3(.1, .2, .3), 
                      UnitInertia(1.2,1.1,1.3,.01,.02,.03)));
    Body::Rigid massless(MassProperties(0,Vec3(0),UnitInertia(0)));
    
    Vec3 inboard(0.1, 0.5, -1);
    Vec3 outboard(0.2, -0.2, 0);
    MobilizedBody::Bushing p1(matter.updGround(), inboard, lumpy, outboard);

    Rotation axisX(Pi/2, YAxis); // rotate z into x
    Rotation axisY(-Pi/2, XAxis); // rotate z into y
    Rotation axisZ; // leave z where it is
    MobilizedBody::Translation dummy0(matter.updGround(), inboard,
        massless, Transform());
    MobilizedBody::Pin dummy1(dummy0, Transform(axisX,Vec3(0)),
        massless, Transform(axisX, Vec3(0)));
    MobilizedBody::Pin dummy2(dummy1, Transform(axisY, Vec3(0)),
        massless, Transform(axisY, Vec3(0)));
    MobilizedBody::Pin p2(dummy2, Transform(axisZ, Vec3(0)),
        lumpy, Transform(axisZ, outboard));

    MobilizedBody::Weld c1(p1, inboard, lumpy, outboard);
    MobilizedBody::Weld c2(p2, inboard, lumpy, outboard);

    c1.addBodyDecoration(Transform(), DecorativeBrick(.1*Vec3(1,2,3))
        .setColor(Cyan).setOpacity(0.2));
    c2.addBodyDecoration(Transform(), DecorativeBrick(.1*Vec3(1,2,3))
        .setColor(Red).setRepresentation(DecorativeGeometry::DrawWireframe));

    //Visualizer viz(system); viz.setBackgroundType(Visualizer::SolidColor);
    //system.addEventReporter(new Visualizer::Reporter(viz, 1./30));

    State state = system.realizeTopology();
    p1.setQ(state, Vec6(.1,.2,.3,1,2,3));
    dummy0.setTranslation(state, Vec3(1,2,3));
    dummy1.setAngle(state, .1);
    dummy2.setAngle(state, .2);
    p2.setAngle(state, .3);

    p1.setU(state, Vec6(1, 2, 3, -.1, -.2, -.3));
    dummy0.setVelocity(state, Vec3(-.1, -.2, -.3));
    dummy1.setRate(state, 1);
    dummy2.setRate(state, 2);
    p2.setRate(state, 3);

    system.realize(state, Stage::Velocity);
    system.project(state, 1e-10);

    SimTK_TEST_EQ(c1.getBodyTransform(state), c2.getBodyTransform(state));
    SimTK_TEST_EQ(c1.getBodyVelocity(state), c2.getBodyVelocity(state));
    
    // Simulate it and see if both pendulums behave identically.
    
    RungeKuttaMersonIntegrator integ(system);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(5);
    system.realize(integ.getState(), Stage::Velocity);
    SimTK_TEST_EQ_TOL(c1.getBodyTransform(integ.getState()), 
                  c2.getBodyTransform(integ.getState()), 1e-10);
    SimTK_TEST_EQ_TOL(c1.getBodyVelocity(integ.getState()), 
                      c2.getBodyVelocity(integ.getState()), 1e-10);
}

// A helper function to replicate the translation-rotational coupling of
// MobilizedBody::CantileverFreeBeam associated with X and Y beam deflections.
class BeamDeflectionFunction : public Function {
public:
    BeamDeflectionFunction(const Real& beamLength, bool negate) {
        m_deflectionCoefficient = (2.0 / 3.0) * beamLength;
        m_sign = negate ? -1.0 : 1.0;
    }

    Real calcValue(const Vector& x) const override {
        return x[0] + m_sign * m_deflectionCoefficient * x[1];
    }

    Real calcDerivative(const Array_<int>& derivComponents,
            const Vector& x) const override {
        if (derivComponents.size() == 1) {
            if (derivComponents[0] == 0) {
                return 1;
            } else if (derivComponents[0] == 1) {
                return m_sign * m_deflectionCoefficient;
            }
        }
        return 0;

    }

    int getArgumentSize() const override {
        return 2;
    }

    int getMaxDerivativeOrder() const override {
        return 2;
    }

private:
    Real m_deflectionCoefficient;
    Real m_sign;
};

// A helper function to replicate the translation-rotational coupling of
// MobilizedBody::CantileverFreeBeam associated with Z beam displacement.
class BeamDisplacementFunction : public Function {
public:
    BeamDisplacementFunction(const Real& beamLength) :
            m_beamLength(beamLength) {
        m_displacementCoefficient = (4.0 / 15.0) * m_beamLength;
    }

    Real calcValue(const Vector& x) const override {
        return x[0] - m_beamLength +
            m_displacementCoefficient * (x[1]*x[1] + x[2]*x[2]);
    }

    Real calcDerivative(const Array_<int>& derivComponents,
            const Vector& x) const override {
        if (derivComponents.size() == 1) {
            if (derivComponents[0] == 0) {
                return 1;
            } else if (derivComponents[0] == 1) {
                return m_displacementCoefficient * 2 * x[1];
            } else if (derivComponents[0] == 2) {
                return m_displacementCoefficient * 2 * x[2];
            }
        } else if (derivComponents.size() == 2) {
            if (derivComponents[0] == 1 && derivComponents[1] == 1) {
                return 2 * m_displacementCoefficient;
            } else if (derivComponents[0] == 2 && derivComponents[1] == 2) {
                return 2 * m_displacementCoefficient;
            } else {
                return 0;
            }
        }
        return 0;
    }

    int getArgumentSize() const override {
        return 3;
    }

    int getMaxDerivativeOrder() const override {
        return 2;
    }

private:
    Real m_beamLength;
    Real m_displacementCoefficient;
};

// Create two pendulums with cantilever-free beam joints, one where the
// cantilever-free beam is faked with a MobilizedBody::Gimbal and a
// MobilizedBody::Translation where the translational coordinates are coupled
// to the rotational coordinates via a Constraint::CoordinateCoupler using the
// same beam equations relationships as the MobilizedBody::CantileverFreeBeam.
// Both pendulums should provide identical generalized coordinates and speeds.
void testCantileverFreeBeam() {

    // Define the system.
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.81, 0));

    // Mass properties and mobilizer frames.
    Body::Rigid lumpy(MassProperties(3.1, Vec3(.1, .2, .3),
                      UnitInertia(1.2, 1.1, 1.3, .01, .02, .03)));
    Body::Rigid massless(MassProperties(0,Vec3(0),UnitInertia(0)));
    Transform inboard(Vec3(0.1, 0.5, -1));
    Transform outboard(Vec3(0.2, -0.2, 0));

    // Undeflected beam length.
    Real length = 1.23;

    // The first pendulum with the true cantilever-free beam mobilizer.
    MobilizedBody::CantileverFreeBeam pendulum1(matter.updGround(), inboard,
        lumpy, outboard, length);

    // The second pendulum with the faked cantilever-free beam mobilizer.
    MobilizedBody::Translation dummy(matter.updGround(), inboard,
        massless, Transform());
    MobilizedBody::Gimbal pendulum2(dummy, Transform(),
        lumpy, outboard);

    // Coupling of X translation with Y rotation.
    std::vector<MobilizedBodyIndex> mobilizers(2);
    std::vector<MobilizerQIndex> coordinates(2);
    mobilizers[0]  = dummy.getMobilizedBodyIndex();
    mobilizers[1]  = pendulum2.getMobilizedBodyIndex();
    coordinates[0] = MobilizerQIndex(0); // X translation
    coordinates[1] = MobilizerQIndex(1); // Y rotation
    Constraint::CoordinateCoupler coupler1(matter,
        new BeamDeflectionFunction(length, true), mobilizers, coordinates);

    // Coupling of Y translation with X rotation.
    mobilizers[0]  = dummy.getMobilizedBodyIndex();
    mobilizers[1]  = pendulum2.getMobilizedBodyIndex();
    coordinates[0] = MobilizerQIndex(1); // Y translation
    coordinates[1] = MobilizerQIndex(0); // X rotation
    Constraint::CoordinateCoupler coupler2(matter,
        new BeamDeflectionFunction(length, false), mobilizers, coordinates);

    // Coupling of Z translation with X and Y rotations.
    std::vector<MobilizedBodyIndex> mobilizers3(3);
    std::vector<MobilizerQIndex> coordinates3(3);
    mobilizers3[0]  = dummy.getMobilizedBodyIndex();
    mobilizers3[1]  = pendulum2.getMobilizedBodyIndex();
    mobilizers3[2]  = pendulum2.getMobilizedBodyIndex();
    coordinates3[0] = MobilizerQIndex(2); // Z translation
    coordinates3[1] = MobilizerQIndex(0); // X rotation
    coordinates3[2] = MobilizerQIndex(1); // Y rotation
    Constraint::CoordinateCoupler coupler3(matter,
        new BeamDisplacementFunction(length), mobilizers3, coordinates3);

    // Test bodies.
    MobilizedBody::Weld c1(pendulum1, inboard, lumpy, outboard);
    MobilizedBody::Weld c2(pendulum2, inboard, lumpy, outboard);
    c1.addBodyDecoration(Transform(), DecorativeBrick(.1*Vec3(1,2,3))
        .setColor(Cyan).setOpacity(0.2));
    c2.addBodyDecoration(Transform(), DecorativeBrick(.1*Vec3(1,2,3))
        .setColor(Red).setRepresentation(DecorativeGeometry::DrawWireframe));

    // Visualizer viz(system);
    // viz.setBackgroundType(Visualizer::SolidColor);
    // system.addEventReporter(new Visualizer::Reporter(viz, 1./30));

    // Construct the system, prescribe the motion of the second pendulum,
    // realize the system to the velocity stage, and project coordinates.
    State state = system.realizeTopology();
    pendulum2.setQ(state, Vec3(0.1, 0.2, 0.3));
    pendulum2.setU(state, Vec3(-0.1, -0.2, -0.3));
    system.realize(state, Stage::Velocity);
    system.project(state, 1e-10);

    // Now that constraints are satisfied for the second pendulum, set the
    // coordinates of the first pendulum to match the coordinates of the second
    // pendulum. This will guarantee that the kinematics of both pendulums are
    // identical.
    pendulum1.setQ(state, pendulum2.getQ(state));
    pendulum1.setU(state, pendulum2.getU(state));
    system.realize(state, Stage::Velocity);
    // Project again, just to be sure.
    system.project(state, 1e-10);

    // Check that the state is identical.
    SimTK_TEST_EQ(pendulum1.getQ(state), pendulum2.getQ(state));
    SimTK_TEST_EQ(pendulum1.getU(state), pendulum2.getU(state));

    // Check that the body transforms and velocities are identical.
    SimTK_TEST_EQ_TOL(c1.getBodyTransform(state),
        c2.getBodyTransform(state), 1e-10);
    SimTK_TEST_EQ_TOL(c1.getBodyVelocity(state),
        c2.getBodyVelocity(state), 1e-10);

    // Simulate it and see if both pendulums behave identically.
    RungeKuttaMersonIntegrator integ(system);
    integ.setAccuracy(1e-10);
    integ.setConstraintTolerance(1e-10);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(5);
    system.realize(integ.getState(), Stage::Velocity);
    SimTK_TEST_EQ_TOL(c1.getBodyTransform(integ.getState()),
                      c2.getBodyTransform(integ.getState()), 1e-10);
    SimTK_TEST_EQ_TOL(c1.getBodyVelocity(integ.getState()),
                      c2.getBodyVelocity(integ.getState()), 1e-10);
}



int main() {
    SimTK_START_TEST("TestMobilizedBody");
        SimTK_SUBTEST(testCalculationMethods);
        SimTK_SUBTEST(testFittingMethods);
        SimTK_SUBTEST(testWeld);
        SimTK_SUBTEST(testGimbal);
        SimTK_SUBTEST(testBushing);
        SimTK_SUBTEST(testCantileverFreeBeam);
    SimTK_END_TEST();
}
