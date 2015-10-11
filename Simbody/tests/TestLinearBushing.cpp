/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2009-12 Stanford University and the Authors.        *
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

/**@file
 * Test the Force::LinearBushing force element.
 */

#include "SimTKsimbody.h"
#include "SimTKcommon/Testing.h"

#include <cstdio>
#include <exception>
#include <iostream>
using std::cout; using std::endl;

//#define VISUALIZE
#ifdef VISUALIZE
    #define WAIT_FOR_INPUT(str) \
        do {printf(str); getchar();} while(false)
    #define REPORT(state) \
        do {viz.report(state);viz.zoomCameraToShowAllGeometry();} while (false)
#else
    #define WAIT_FOR_INPUT(str)
    #define REPORT(state)
#endif


using namespace SimTK;

void testParameterSetting() {
    MultibodySystem         system;
    SimbodyMatterSubsystem  matter(system);
    GeneralForceSubsystem   forces(system);

    // This is the frame on body1.
    const Transform X_B1F(Test::randRotation(), Test::randVec3());
    // This is the frame on body 2.
    const Transform X_B2M(Test::randRotation(), Test::randVec3());
    // Material properties.
    const Vec6 k = 10*Vec6(1,2,1,1,5,1);
    const Vec6 c =  1.3*Vec6(1,.1,1,.11,1,11);

    const Real Mass = 1.234;
    Body::Rigid aBody(MassProperties(Mass, Vec3(.1,.2,.3),
                          Mass*Inertia(1,1.1,1.2,0.01,0.02,0.03)));

    MobilizedBody::Free body1(matter.Ground(), Transform(),
                              aBody,           Transform());
    MobilizedBody::Free body2(matter.Ground(), Transform(),
                              aBody,           Transform());

    Force::LinearBushing bushing
        (forces, body1, Vec3(1,2,3),
                 body2, Vec3(-2,-3,-4),
         Vec6(5), Vec6(7));
    SimTK_TEST_EQ(bushing.getDefaultFrameOnBody1(),Transform(Vec3(1,2,3)));
    SimTK_TEST_EQ(bushing.getDefaultFrameOnBody2(),Transform(Vec3(-2,-3,-4)));
    SimTK_TEST_EQ(bushing.getDefaultStiffness(),Vec6(5,5,5,5,5,5));
    SimTK_TEST_EQ(bushing.getDefaultDamping(),Vec6(7,7,7,7,7,7));

    bushing.setDefaultFrameOnBody1(Vec3(-1,-2,-3))
           .setDefaultFrameOnBody2(Vec3(2,3,4))
           .setDefaultStiffness(Vec6(9))
           .setDefaultDamping(Vec6(11));
    SimTK_TEST_EQ(bushing.getDefaultFrameOnBody1(),Transform(Vec3(-1,-2,-3)));
    SimTK_TEST_EQ(bushing.getDefaultFrameOnBody2(),Transform(Vec3(2,3,4)));
    SimTK_TEST_EQ(bushing.getDefaultStiffness(),Vec6(9));
    SimTK_TEST_EQ(bushing.getDefaultDamping(),Vec6(11));

    system.realizeTopology();
    State state = system.getDefaultState();
    system.realizeModel(state);

    // See if defaults transfer to default state properly.
    SimTK_TEST_EQ(bushing.getFrameOnBody1(state),Transform(Vec3(-1,-2,-3)));
    SimTK_TEST_EQ(bushing.getFrameOnBody2(state),Transform(Vec3(2,3,4)));
    SimTK_TEST_EQ(bushing.getStiffness(state),Vec6(9));
    SimTK_TEST_EQ(bushing.getDamping(state),Vec6(11));

    bushing.setFrameOnBody1(state, X_B1F)
           .setFrameOnBody2(state, X_B2M)
           .setStiffness(state, k)
           .setDamping(state, c);

    SimTK_TEST_EQ(bushing.getFrameOnBody1(state),X_B1F);
    SimTK_TEST_EQ(bushing.getFrameOnBody2(state),X_B2M);
    SimTK_TEST_EQ(bushing.getStiffness(state),k);
    SimTK_TEST_EQ(bushing.getDamping(state),c);


}

// Here we're going to build a chain like this:
//
//   Ground --> body1 ==> body2
//
// where body1 is Free, and body2 is connected to body1 by
// a series of mobilizers designed to have the same kinematics
// as a LinearBushing element connected between them.
//
// We'll verify that the LinearBushing calculates kinematics
// that exactly matches the composite joint. Then we'll run
// for a while and see that energy+dissipation is conserved.
void testKinematicsAndEnergyConservation() {
    // Build system.
    MultibodySystem         system;
    SimbodyMatterSubsystem  matter(system);
    GeneralForceSubsystem   forces(system);

    // This is the frame on body1.
    const Transform X_B1F(Test::randRotation(), Test::randVec3());
    // This is the frame on body 2.
    const Transform X_B2M(Test::randRotation(), Test::randVec3());
    // Material properties.
    const Vec6 k = 10*Vec6(1,1,1,1,1,1);
    const Vec6 c =  1*Vec6(1,1,1,1,1,1);

    // Use ugly mass properties to make sure we test all the terms.
    const Real Mass = 1.234;
    const Vec3 HalfShape = Vec3(1,.5,.25)/2;

    Body::Rigid aBody(MassProperties(Mass, Vec3(.1,.2,.3),
                          Mass*Inertia(1,1.1,1.2,0.01,0.02,0.03)));
    aBody.addDecoration(Transform(), DecorativeEllipsoid(HalfShape)
                                            .setOpacity(0.25)
                                            .setColor(Blue));
    aBody.addDecoration(X_B1F,
                DecorativeFrame(0.5).setColor(Red));
    aBody.addDecoration(X_B2M,
                DecorativeFrame(0.5).setColor(Green));

    // Ground attachement frame.
    const Transform X_GA(
        Test::randRotation(), Test::randVec3());

    MobilizedBody::Free body1(matter.Ground(), X_GA,
                              aBody,           X_B2M);

    // This is to keep the system from flying away.
    Force::TwoPointLinearSpring(forces,matter.Ground(),Vec3(0),body1,Vec3(0),
                                10,0);

    Body::Rigid massless(MassProperties(0, Vec3(0), Inertia(0)));
    const Rotation ZtoX(Pi/2, YAxis);
    const Rotation ZtoY(-Pi/2, XAxis);
    MobilizedBody::Cartesian dummy1(body1, X_B1F,
                                    massless, Transform());
    MobilizedBody::Pin dummy2(dummy1,   ZtoX,
                              massless, ZtoX);
    MobilizedBody::Pin dummy3(dummy2,   ZtoY,
                              massless, ZtoY);
    MobilizedBody::Pin body2 (dummy3,   Transform(),    // about Z
                              aBody, X_B2M);

    // Set the actual parameters in the State below.
    Force::LinearBushing bushing
        (forces, body1, body2, Vec6(0), Vec6(0));

#ifdef VISUALIZE
    Visualizer viz(system);
    viz.setBackgroundType(Visualizer::SolidColor);
    system.addEventReporter(new Visualizer::Reporter(viz, 0.01));
#endif

    // Initialize the system and state.

    system.realizeTopology();
    State state = system.getDefaultState();

    bushing.setFrameOnBody1(state, X_B1F)
           .setFrameOnBody2(state, X_B2M)
           .setStiffness(state, k)
           .setDamping(state, c);

    REPORT(state);
    WAIT_FOR_INPUT("\nDefault state -- hit ENTER\n");

    state.updQ() = Test::randVector(state.getNQ());
    state.updU() = Test::randVector(state.getNU());
    state.updU()(3,3) = 0; // kill translational velocity

    const Real Accuracy = 1e-6;
    RungeKuttaMersonIntegrator integ(system);
    integ.setAccuracy(Accuracy);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    State istate = integ.getState();
    system.realize(istate, Stage::Velocity);
    Vector xyz  = istate.getQ()(7,3);    Vector ang  = istate.getQ()(10,3);
    Vector xyzd = istate.getQDot()(7,3); Vector angd = istate.getQDot()(10,3);
    Vec6 mq  = Vec6(ang[0],ang[1],ang[2],xyz[0],xyz[1],xyz[2]);
    Vec6 mqd = Vec6(angd[0],angd[1],angd[2],xyzd[0],xyzd[1],xyzd[2]);

    SimTK_TEST_EQ(bushing.getQ(istate), mq);
    SimTK_TEST_EQ(bushing.getQDot(istate), mqd);

    const Real initialEnergy = system.calcEnergy(istate);

    SimTK_TEST( bushing.getDissipatedEnergy(istate) == 0 );

    REPORT(integ.getState());

    cout << "t=" << integ.getTime()
         << "\nE=" << initialEnergy
         << "\nmobilizer q=" << mq
         << "\nbushing   q=" << bushing.getQ(istate)
         << "\nmobilizer qd=" << mqd
         << "\nbushing   qd=" << bushing.getQDot(istate)
         << endl;
    WAIT_FOR_INPUT("After initialize -- hit ENTER\n");

    // Simulate it.
    ts.stepTo(5.0);
    istate = integ.getState();
    system.realize(istate, Stage::Velocity);
    xyz  = istate.getQ()(7,3);    ang  = istate.getQ()(10,3);
    xyzd = istate.getQDot()(7,3); angd = istate.getQDot()(10,3);
    mq  = Vec6(ang[0],ang[1],ang[2],xyz[0],xyz[1],xyz[2]);
    mqd = Vec6(angd[0],angd[1],angd[2],xyzd[0],xyzd[1],xyzd[2]);

    SimTK_TEST_EQ(bushing.getQ(istate), mq);
    SimTK_TEST_EQ(bushing.getQDot(istate), mqd);

    // This should account for all the energy.
    const Real finalEnergy = system.calcEnergy(istate)
                                + bushing.getDissipatedEnergy(istate);

    SimTK_TEST_EQ_TOL(initialEnergy, finalEnergy, 10*Accuracy);

    // Let's find everything and see if the bushing agrees.
    const Transform& X_GB1 = body1.getBodyTransform(istate);
    const Transform& X_GB2 = body2.getBodyTransform(istate);
    const Transform X_GF = X_GB1*X_B1F;
    const Transform X_GM = X_GB2*X_B2M;
    SimTK_TEST_EQ(X_GF, bushing.getX_GF(istate));
    SimTK_TEST_EQ(X_GM, bushing.getX_GM(istate));
    SimTK_TEST_EQ(~X_GF*X_GM, bushing.getX_FM(istate));

    // Calculate velocities and ask the bushing for same.
    const SpatialVec& V_GB1 = body1.getBodyVelocity(istate);
    const SpatialVec& V_GB2 = body2.getBodyVelocity(istate);
    const SpatialVec
        V_GF(V_GB1[0], body1.findStationVelocityInGround(istate,X_B1F.p()));
    const SpatialVec
        V_GM(V_GB2[0], body2.findStationVelocityInGround(istate,X_B2M.p()));

    SimTK_TEST_EQ(V_GF, bushing.getV_GF(istate));
    SimTK_TEST_EQ(V_GM, bushing.getV_GM(istate));

    const SpatialVec
        V_FM(~X_B1F.R()*body2.findBodyAngularVelocityInAnotherBody(istate,body1),
             ~X_B1F.R()*body2.findStationVelocityInAnotherBody(istate,X_B2M.p(),body1));

    SimTK_TEST_EQ(V_FM, bushing.getV_FM(istate));

    cout << "t=" << integ.getTime()
         << "\nE=" << system.calcEnergy(istate)
         << "\nE-XW=" << finalEnergy << " final-init=" << finalEnergy-initialEnergy
         << "\nmobilizer q=" << mq
         << "\nbushing   q=" << bushing.getQ(istate)
         << "\nmobilizer qd=" << mqd
         << "\nbushing   qd=" << bushing.getQDot(istate)
         << endl;

    WAIT_FOR_INPUT("After simulation -- hit ENTER\n");
}


// Here we're going to build a chain like this:
//
//   Ground --> body1 ==> body2
//
// where body1 is Free, and body2 is connected to body1 by
// a Bushing mobilizer, which should have the same kinematics
// as a LinearBushing element connected between them.
//
// We'll verify that the LinearBushing calculates kinematics
// that exactly matches the Bushing mobilizer. Then we'll run
// for a while and see that energy+dissipation is conserved.
void testKinematicsAndEnergyConservationUsingBushingMobilizer() {
    // Build system.
    MultibodySystem         system;
    SimbodyMatterSubsystem  matter(system);
    GeneralForceSubsystem   forces(system);

    // This is the frame on body1.
    const Transform X_B1F(Test::randRotation(), Test::randVec3());
    // This is the frame on body 2.
    const Transform X_B2M(Test::randRotation(), Test::randVec3());
    // Material properties.
    const Vec6 k = 10*Vec6(1,1,1,1,1,1);
    const Vec6 c =  1*Vec6(1,1,1,1,1,1);

    // Use ugly mass properties to make sure we test all the terms.
    const Real Mass = 1.234;
    const Vec3 HalfShape = Vec3(1,.5,.25)/2;

    Body::Rigid aBody(MassProperties(Mass, Vec3(.1,.2,.3),
                          Mass*Inertia(1,1.1,1.2,0.01,0.02,0.03)));
    aBody.addDecoration(Transform(), DecorativeEllipsoid(HalfShape)
                                            .setOpacity(0.25)
                                            .setColor(Blue));
    aBody.addDecoration(X_B1F,
                DecorativeFrame(0.5).setColor(Red));
    aBody.addDecoration(X_B2M,
                DecorativeFrame(0.5).setColor(Green));

    // Ground attachement frame.
    const Transform X_GA(
        Test::randRotation(), Test::randVec3());

    MobilizedBody::Free body1(matter.Ground(), X_GA,
                              aBody,           X_B2M);

    // This is to keep the system from flying away.
    Force::TwoPointLinearSpring(forces,matter.Ground(),Vec3(0),body1,Vec3(0),
                                10,0);

    MobilizedBody::Bushing body2(body1, X_B1F,
                                 aBody, X_B2M);

    // Set the actual parameters in the State below.
    Force::LinearBushing bushing
        (forces, body1, body2, Vec6(0), Vec6(0));

#ifdef VISUALIZE
    Visualizer viz(system);
    viz.setBackgroundType(Visualizer::SolidColor);
    system.addEventReporter(new Visualizer::Reporter(viz, 0.01));
#endif

    // Initialize the system and state.

    system.realizeTopology();
    State state = system.getDefaultState();

    bushing.setFrameOnBody1(state, X_B1F)
           .setFrameOnBody2(state, X_B2M)
           .setStiffness(state, k)
           .setDamping(state, c);

    REPORT(state);
    WAIT_FOR_INPUT("\nDefault state -- hit ENTER\n");

    state.updQ() = Test::randVector(state.getNQ());
    state.updU() = Test::randVector(state.getNU());
    state.updU()(3,3) = 0; // kill translational velocity

    const Real Accuracy = 1e-6;
    RungeKuttaMersonIntegrator integ(system);
    integ.setAccuracy(Accuracy);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    State istate = integ.getState();
    system.realize(istate, Stage::Velocity);

    // Get the q's and qdot's from the Bushing mobilizer.
    Vec6 mq  = body2.getQ(istate);
    Vec6 mqd = body2.getQDot(istate);

    SimTK_TEST_EQ(bushing.getQ(istate), mq);
    SimTK_TEST_EQ(bushing.getQDot(istate), mqd);

    const Real initialEnergy = system.calcEnergy(istate);

    SimTK_TEST( bushing.getDissipatedEnergy(istate) == 0 );

    REPORT(integ.getState());

    cout << "t=" << integ.getTime()
         << "\nE=" << initialEnergy
         << "\nmobilizer q=" << mq
         << "\nbushing   q=" << bushing.getQ(istate)
         << "\nmobilizer qd=" << mqd
         << "\nbushing   qd=" << bushing.getQDot(istate)
         << endl;
    WAIT_FOR_INPUT("After initialize -- hit ENTER\n");

    // Simulate it.
    ts.stepTo(5.0);
    istate = integ.getState();
    system.realize(istate, Stage::Velocity);
    mq  = body2.getQ(istate);
    mqd = body2.getQDot(istate);

    SimTK_TEST_EQ(bushing.getQ(istate), mq);
    SimTK_TEST_EQ(bushing.getQDot(istate), mqd);

    // This should account for all the energy.
    const Real finalEnergy = system.calcEnergy(istate)
                                + bushing.getDissipatedEnergy(istate);

    SimTK_TEST_EQ_TOL(initialEnergy, finalEnergy, 10*Accuracy);

    // Let's find everything and see if the bushing agrees.
    const Transform& X_GB1 = body1.getBodyTransform(istate);
    const Transform& X_GB2 = body2.getBodyTransform(istate);
    const Transform X_GF = X_GB1*X_B1F;
    const Transform X_GM = X_GB2*X_B2M;
    SimTK_TEST_EQ(X_GF, bushing.getX_GF(istate));
    SimTK_TEST_EQ(X_GM, bushing.getX_GM(istate));
    SimTK_TEST_EQ(~X_GF*X_GM, bushing.getX_FM(istate));

    // Calculate velocities and ask the bushing for same.
    const SpatialVec& V_GB1 = body1.getBodyVelocity(istate);
    const SpatialVec& V_GB2 = body2.getBodyVelocity(istate);
    const SpatialVec
        V_GF(V_GB1[0], body1.findStationVelocityInGround(istate,X_B1F.p()));
    const SpatialVec
        V_GM(V_GB2[0], body2.findStationVelocityInGround(istate,X_B2M.p()));

    SimTK_TEST_EQ(V_GF, bushing.getV_GF(istate));
    SimTK_TEST_EQ(V_GM, bushing.getV_GM(istate));

    const SpatialVec
        V_FM(~X_B1F.R()*body2.findBodyAngularVelocityInAnotherBody(istate,body1),
             ~X_B1F.R()*body2.findStationVelocityInAnotherBody(istate,X_B2M.p(),body1));

    SimTK_TEST_EQ(V_FM, bushing.getV_FM(istate));

    cout << "t=" << integ.getTime()
         << "\nE=" << system.calcEnergy(istate)
         << "\nE-XW=" << finalEnergy << " final-init=" << finalEnergy-initialEnergy
         << "\nmobilizer q=" << mq
         << "\nbushing   q=" << bushing.getQ(istate)
         << "\nmobilizer qd=" << mqd
         << "\nbushing   qd=" << bushing.getQDot(istate)
         << endl;

    WAIT_FOR_INPUT("After simulation -- hit ENTER\n");
}


// Here we're going to build a system containing two parallel multibody
// trees, each consisting of a single body connected to ground. For the
// first system the body is on a Free mobilizer and a LinearBushing
// connects the body and Ground. In the second, a series of mobilizers
// exactly mimics the kinematics, and a set of mobility springs and
// dampers are used to mimic the forces that should be produced by the
// bushing. We'll then evaluate in some arbitrary configuration and
// make sure the mobilizer reaction forces match the corresponding
// LinearBushing force.
void testForces() {
    // Create the system.
    MultibodySystem         system;
    SimbodyMatterSubsystem  matter(system);
    GeneralForceSubsystem   forces(system);
    Force::UniformGravity   gravity(forces, matter, 0*Vec3(0, -9.8, 0));

    const Real Mass = 1;
    const Vec3 HalfShape = Vec3(1,.5,.25)/2;
    const Transform BodyAttach(Rotation(), Vec3(HalfShape[0],0,0));
    Body::Rigid brickBody(MassProperties(Mass, Vec3(.1,.2,.3),
                                Mass*Inertia(1,1.1,1.2,0.01,0.02,0.03)));
    //Body::Rigid brickBody(MassProperties(Mass, Vec3(0),
    //                        Mass*UnitInertia::ellipsoid(HalfShape)));
    brickBody.addDecoration(Transform(), DecorativeEllipsoid(HalfShape)
                                            .setOpacity(0.25)
                                            .setColor(Blue));
    brickBody.addDecoration(BodyAttach,
                DecorativeFrame(0.5).setColor(Red));

    MobilizedBody::Free brick1(matter.Ground(), Transform(),
                               brickBody,       BodyAttach);

    Body::Rigid massless(MassProperties(0, Vec3(0), Inertia(0)));
    const Rotation ZtoX(Pi/2, YAxis);
    const Rotation ZtoY(-Pi/2, XAxis);

    // This sequence is used to match the kinematics of a "forward" direction
    // bushing where Ground is body 1.
    //MobilizedBody::Cartesian dummy1(matter.Ground(), Vec3(1,1,1)+Vec3(2,0,0),
    //                                massless, Transform());
    //MobilizedBody::Pin dummy2(dummy1,   ZtoX,
    //                          massless, ZtoX);
    //MobilizedBody::Pin dummy3(dummy2,   ZtoY,
    //                          massless, ZtoY);
    //MobilizedBody::Pin brick2(dummy3,   Transform(),    // about Z
    //                          brickBody, BodyAttach);

    // This sequence is used to match the kinematics of a "reverse" direction
    // bushing where Ground is body 2 (this is trickier).
    MobilizedBody::Pin dummy1(matter.Ground(), Vec3(1,1,1)+Vec3(2,0,0),
                              massless, Transform(), MobilizedBody::Reverse);
    MobilizedBody::Pin dummy2(dummy1,   ZtoY,
                              massless, ZtoY, MobilizedBody::Reverse);
    MobilizedBody::Pin dummy3(dummy2,   ZtoX,
                              massless, ZtoX, MobilizedBody::Reverse);
    MobilizedBody::Cartesian brick2(dummy3,   Transform(),
                                    brickBody, BodyAttach, MobilizedBody::Reverse);

    const Vec6 k = Test::randVec<6>().abs(); // must be > 0
    const Vec6 c =  Test::randVec<6>().abs(); // "
    Transform GroundAttach(Rotation(), Vec3(1,1,1));
    matter.Ground().updBody().addDecoration(GroundAttach,
                DecorativeFrame(0.5).setColor(Green));

    const Vec6 k1 = k;
    const Vec6 c1 = c;
    // This is the Forward version
    //Force::LinearBushing bushing
    //    (forces, matter.Ground(), GroundAttach,
    //     brick1, BodyAttach, k1, c1);

    // This is the reverse version -- the moving body is the first body
    // for the bushing; Ground is second. This is a harder test because
    // the F frame is moving.
    Force::LinearBushing bushing
        (forces, brick1, BodyAttach,
         matter.Ground(), GroundAttach,
         k1, c1);

    Force::LinearBushing bushing2
        (forces, brick1, brick2, k, c);

    Transform GroundAttach2(Rotation(), Vec3(1,1,1) + Vec3(2,0,0));
    matter.Ground().updBody().addDecoration(GroundAttach2,
                DecorativeFrame(0.5).setColor(Green));
    const Vec6 k2 = k;
    const Vec6 c2 = c;

    // This is what you would do for a forward-direction set of
    // mobilizers to match the forward bushing.
    //Force::MobilityLinearSpring kqx(forces, dummy2, 0, k2[0], 0);
    //Force::MobilityLinearDamper cqx(forces, dummy2, 0, c2[0]);
    //Force::MobilityLinearSpring kqy(forces, dummy3, 0, k2[1], 0);
    //Force::MobilityLinearDamper cqy(forces, dummy3, 0, c2[1]);
    //Force::MobilityLinearSpring kqz(forces, brick2, 0, k2[2], 0);
    //Force::MobilityLinearDamper cqz(forces, brick2, 0, c2[2]);
    //Force::MobilityLinearSpring kpx(forces, dummy1, 0, k2[3], 0);
    //Force::MobilityLinearDamper cpx(forces, dummy1, 0, c2[3]);
    //Force::MobilityLinearSpring kpy(forces, dummy1, 1, k2[4], 0);
    //Force::MobilityLinearDamper cpy(forces, dummy1, 1, c2[4]);
    //Force::MobilityLinearSpring kpz(forces, dummy1, 2, k2[5], 0);
    //Force::MobilityLinearDamper cpz(forces, dummy1, 2, c2[5]);

    // This is the set of force elements that mimics the bushing hooked
    // up in the reverse direction.
    Force::MobilityLinearSpring kqx(forces, dummy3, 0, k2[0], 0);
    Force::MobilityLinearDamper cqx(forces, dummy3, 0, c2[0]);
    Force::MobilityLinearSpring kqy(forces, dummy2, 0, k2[1], 0);
    Force::MobilityLinearDamper cqy(forces, dummy2, 0, c2[1]);
    Force::MobilityLinearSpring kqz(forces, dummy1, 0, k2[2], 0);
    Force::MobilityLinearDamper cqz(forces, dummy1, 0, c2[2]);
    Force::MobilityLinearSpring kpx(forces, brick2, 0, k2[3], 0);
    Force::MobilityLinearDamper cpx(forces, brick2, 0, c2[3]);
    Force::MobilityLinearSpring kpy(forces, brick2, 1, k2[4], 0);
    Force::MobilityLinearDamper cpy(forces, brick2, 1, c2[4]);
    Force::MobilityLinearSpring kpz(forces, brick2, 2, k2[5], 0);
    Force::MobilityLinearDamper cpz(forces, brick2, 2, c2[5]);

#ifdef VISUALIZE
    Visualizer viz(system);
    viz.setBackgroundType(Visualizer::SolidColor);
    system.addEventReporter(new Visualizer::Reporter(viz, 0.01));
#endif


    // Initialize the system and state.

    system.realizeTopology();
    State state = system.getDefaultState();

    REPORT(state);

    Rotation RR = Test::randRotation();
    brick1.setQToFitTransform(state, RR);

    Vec3 qRRinv;
    qRRinv = Rotation(~RR).convertRotationToBodyFixedXYZ();

    brick2.setQToFitTransform(state, ~RR*Vec3(-1,-1,-1));
    dummy3.setOneQ(state, 0, qRRinv[0]);
    dummy2.setOneQ(state, 0, qRRinv[1]);
    dummy1.setOneQ(state, 0, qRRinv[2]);


    REPORT(state);
    WAIT_FOR_INPUT("testForces() trial pose -- hit ENTER\n");

    system.realize(state,Stage::Acceleration);

    //cout << "\nf=" << bushing.getF(state) << endl;
    //cout << "F_GM=" << bushing.getF_GM(state) << endl;
    //cout << "F_GF=" << bushing.getF_GF(state) << endl;

    Vector_<SpatialVec> reactions;
    matter.calcMobilizerReactionForces(state, reactions);
    //cout << "Reaction force on brick2="
     //    << reactions[brick2.getMobilizedBodyIndex()] << endl;

    SimTK_TEST_EQ(bushing.getF_GF(state),
                  reactions[brick2.getMobilizedBodyIndex()]);
}

// This test is identical to testForces() except we use a reverse Bushing
// mobilizer rather than a sequence of Translation and Pins.
//
// Here we're going to build a system containing two parallel multibody
// trees, each consisting of a single body connected to ground. For the
// first system the body is on a Free mobilizer and a LinearBushing
// connects the body and Ground, although with Ground serving as "body2" for
// the bushing (i.e., the M frame is on Ground). In the second, a reversed
// Bushing mobilizer exactly mimics the kinematics, and a set of mobility
// springs and dampers are used to mimic the forces that should be produced by
// the bushing. We'll then evaluate in some arbitrary configuration and
// make sure the mobilizer reaction forces match the corresponding
// LinearBushing force.
void testForcesUsingReverseBushingMobilizer() {
    // Create the system.
    MultibodySystem         system;
    SimbodyMatterSubsystem  matter(system);
    GeneralForceSubsystem   forces(system);
    Force::UniformGravity   gravity(forces, matter, 0*Vec3(0, -9.8, 0));

    matter.Ground().addBodyDecoration(Transform(),
        DecorativeFrame(.25).setColor(Purple));

    const Real Mass = 1;
    const Vec3 HalfShape = Vec3(1,.5,.25)/2;
    const Transform BodyAttach(Rotation(), Vec3(HalfShape[0],0,0));
    Body::Rigid brickBody(MassProperties(Mass, Vec3(.1,.2,.3),
                                Mass*Inertia(1,1.1,1.2,0.01,0.02,0.03)));
    brickBody.addDecoration(Transform(), DecorativeEllipsoid(HalfShape)
                                            .setOpacity(0.25)
                                            .setColor(Blue));
    brickBody.addDecoration(BodyAttach,
                DecorativeFrame(0.5).setColor(Red));

    MobilizedBody::Free brick1(matter.Ground(), Transform(),
                               brickBody,       BodyAttach);
    brick1.addBodyDecoration(Transform(), DecorativeText("1").setScale(.25));

    // This mobilizer is used to match the kinematics of a "reverse" direction
    // bushing where Ground is body 2.
    MobilizedBody::Bushing brick2(matter.Ground(), Vec3(1,1,1)+Vec3(2,0,0),
                                  brickBody, BodyAttach, MobilizedBody::Reverse);
    brick2.addBodyDecoration(Transform(), DecorativeText("2").setScale(.25));

    const Vec6 k = Test::randVec<6>().abs(); // must be > 0
    const Vec6 c =  Test::randVec<6>().abs(); // "
    Transform GroundAttach(Rotation(), Vec3(1,1,1));
    matter.Ground().updBody().addDecoration(GroundAttach,
                DecorativeFrame(0.5).setColor(Green));

    const Vec6 k1 = k;
    const Vec6 c1 = c;

    // This is reversed -- the moving body is the first body
    // for the bushing; Ground is second. This is a hard test because
    // the F frame is moving.
    Force::LinearBushing bushing
        (forces, brick1, BodyAttach,    // .5,0,0
         matter.Ground(), GroundAttach, // 1,1,1
         k1, c1);

    Force::LinearBushing bushing2
        (forces, brick1, brick2, k, c);

    Transform GroundAttach2(Rotation(), Vec3(1,1,1) + Vec3(2,0,0));
    matter.Ground().updBody().addDecoration(GroundAttach2,
                DecorativeFrame(0.5).setColor(Green));
    const Vec6 k2 = k;
    const Vec6 c2 = c;

    // This is the set of force elements that mimics the bushing hooked
    // up in the reverse direction.
    Force::MobilityLinearSpring kqx(forces, brick2, 0, k2[0], 0);
    Force::MobilityLinearDamper cqx(forces, brick2, 0, c2[0]);
    Force::MobilityLinearSpring kqy(forces, brick2, 1, k2[1], 0);
    Force::MobilityLinearDamper cqy(forces, brick2, 1, c2[1]);
    Force::MobilityLinearSpring kqz(forces, brick2, 2, k2[2], 0);
    Force::MobilityLinearDamper cqz(forces, brick2, 2, c2[2]);
    Force::MobilityLinearSpring kpx(forces, brick2, 3, k2[3], 0);
    Force::MobilityLinearDamper cpx(forces, brick2, 3, c2[3]);
    Force::MobilityLinearSpring kpy(forces, brick2, 4, k2[4], 0);
    Force::MobilityLinearDamper cpy(forces, brick2, 4, c2[4]);
    Force::MobilityLinearSpring kpz(forces, brick2, 5, k2[5], 0);
    Force::MobilityLinearDamper cpz(forces, brick2, 5, c2[5]);

#ifdef VISUALIZE
    Visualizer viz(system);
    viz.setBackgroundType(Visualizer::SolidColor);
    system.addEventReporter(new Visualizer::Reporter(viz, 0.01));
#endif


    // Initialize the system and state.

    system.realizeTopology();
    State state = system.getDefaultState();

    REPORT(state);

    Rotation RR = Test::randRotation();
    brick1.setQToFitTransform(state, RR);
    system.realize(state, Stage::Position);

    cout << "\nbrick1 .q=" << brick1.getQ(state) << endl;
    cout <<   "bushing.q=" << bushing.getQ(state) << endl;

    brick2.setQToFitTransform(state, Transform(RR, -Vec3(1,1,1)));

    cout << "brick2 .q=" << brick2.getQ(state) << endl;


    REPORT(state);
    WAIT_FOR_INPUT("testForcesUsingBushing() trial pose -- hit ENTER\n");

    system.realize(state,Stage::Acceleration);

    //cout << "\nf=" << bushing.getF(state) << endl;
    //cout << "F_GM=" << bushing.getF_GM(state) << endl;
    cout << "F_GF=" << bushing.getF_GF(state) << endl;

    Vector_<SpatialVec> reactions;
    matter.calcMobilizerReactionForces(state, reactions);
    cout << "Reaction force on brick2="
         << reactions[brick2.getMobilizedBodyIndex()] << endl;

    SimTK_TEST_EQ(bushing.getF_GF(state),
                  reactions[brick2.getMobilizedBodyIndex()]);
}

int main() {
    SimTK_START_TEST("TestLinearBushing");
        SimTK_SUBTEST(testParameterSetting);
        SimTK_SUBTEST(testKinematicsAndEnergyConservation);
        SimTK_SUBTEST(testKinematicsAndEnergyConservationUsingBushingMobilizer);
        SimTK_SUBTEST(testForces);
        SimTK_SUBTEST(testForcesUsingReverseBushingMobilizer);
    SimTK_END_TEST();
}

