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

/* Test the Force::Gravity element.
*/

#include "Simbody.h"

#include <cassert>
#include <iostream>
using std::cout; using std::endl;

using namespace SimTK;

//#define USE_VISUALIZER

const Vec3 Cube(.5,.5,.5); // half-dimensions of cube
const Real Mass1=100, Mass2=5, Mass3=1;
const Vec3 Centroid(.5,0,.5);
const Vec3 COM1=Centroid, COM2=Centroid, COM3=Centroid+Vec3(0,.5,0);
const UnitInertia Central(Vec3(.1), Vec3(.05));
// Simbody requires inertias to be expressed about body origin rather than COM.
const Inertia Inertia1=Mass1*Central.shiftFromCentroid(-COM1);
const Inertia Inertia2=Mass2*Central.shiftFromCentroid(-COM2);
const Inertia Inertia3=Mass3*Central.shiftFromCentroid(-COM3); // weird

Body::Rigid body1Info(MassProperties(Mass1, COM1, Inertia1));
Body::Rigid body2Info(MassProperties(Mass2, COM2, Inertia2));
Body::Rigid body3Info(MassProperties(Mass3, COM3, Inertia3));

// Make sure the constructors and default setters and getters work.
static void testConstruction() {
    MultibodySystem         mbs;
    SimbodyMatterSubsystem  matter(mbs);
    GeneralForceSubsystem   forces(mbs);

    Force::Gravity gravity1(forces, matter, -ZAxis, 50);
    SimTK_TEST(gravity1.getDefaultMagnitude()==50);
    SimTK_TEST(gravity1.getDefaultDownDirection()==UnitVec3(0,0,-1));
    SimTK_TEST(gravity1.getDefaultZeroHeight()==0);
    SimTK_TEST(gravity1.getDefaultGravityVector()==Vec3(0,0,-50));
    SimTK_TEST(gravity1.getDefaultBodyIsExcluded(MobodIndex(0)));
    SimTK_TEST(!gravity1.getDefaultBodyIsExcluded(MobodIndex(1)));

    const Vec3 grav2(1,2,3);
    Force::Gravity gravity2(forces, matter, grav2);
    SimTK_TEST_EQ(gravity2.getDefaultMagnitude(),grav2.norm());
    SimTK_TEST_EQ(gravity2.getDefaultDownDirection(),UnitVec3(grav2));
    SimTK_TEST(gravity2.getDefaultZeroHeight()==0);
    SimTK_TEST(gravity2.getDefaultGravityVector()==grav2);

    mbs.setUpDirection(XAxis); const Real mag = 16.75;
    Force::Gravity gravity3(forces, matter, mag);

    SimTK_TEST(gravity3.getDefaultMagnitude()==mag);
    SimTK_TEST(gravity3.getDefaultDownDirection()==UnitVec3(-XAxis));
    SimTK_TEST(gravity3.getDefaultZeroHeight()==0);
    SimTK_TEST(gravity3.getDefaultGravityVector()==Vec3(-mag,0,0));

    // Using the vector constructor with a zero vector should pluck the
    // direction out of the System.
    Force::Gravity gravity4(forces, matter, Vec3(0));
    SimTK_TEST(gravity4.getDefaultMagnitude()==0);
    SimTK_TEST(gravity4.getDefaultDownDirection()==UnitVec3(-XAxis));
    SimTK_TEST(gravity4.getDefaultZeroHeight()==0);
    SimTK_TEST(gravity4.getDefaultGravityVector()==Vec3(0));

    // Make sure Ground can't be included.
    gravity1.setDefaultBodyIsExcluded(MobodIndex(0), false);
    SimTK_TEST(gravity1.getDefaultBodyIsExcluded(MobodIndex(0)));

    // Make sure we can exclude a body and that it doesn't leak over into
    // another body. Note that for defaults we don't yet know how may bodies
    // there might be so we can use arbitrary body numbers.
    SimTK_TEST(!gravity1.getDefaultBodyIsExcluded(MobodIndex(13)));
    gravity1.setDefaultBodyIsExcluded(MobodIndex(13), true);
    SimTK_TEST(gravity1.getDefaultBodyIsExcluded(MobodIndex(13)));
    SimTK_TEST(!gravity1.getDefaultBodyIsExcluded(MobodIndex(5)));

    gravity1.setDefaultGravityVector(Vec3(19,20,21));
    SimTK_TEST_EQ(gravity1.getDefaultGravityVector(), Vec3(19,20,21));

    gravity1.setDefaultMagnitude(3).setDefaultDownDirection(grav2);
    SimTK_TEST(gravity1.getDefaultMagnitude()==3);
    SimTK_TEST_EQ(gravity1.getDefaultDownDirection(), UnitVec3(grav2));
    SimTK_TEST_EQ(gravity1.getDefaultGravityVector(), Real(3)*UnitVec3(grav2));

    gravity1.setDefaultZeroHeight(5.25);
    SimTK_TEST(gravity1.getDefaultZeroHeight() == 5.25);
}

// Make sure we can override default values in the State.
void testParameters() {
    MultibodySystem         mbs;
    SimbodyMatterSubsystem  matter(mbs);
    GeneralForceSubsystem   forces(mbs);

    Force::Gravity gravity(forces, matter, -ZAxis, 49, -5);
    SimTK_TEST(gravity.getDefaultZeroHeight()==-5);

    MobilizedBody::Weld mobod1(matter.Ground(),Vec3(0), body1Info, Vec3(0));
    MobilizedBody::Pin mobod2(mobod1, Vec3(0), body2Info, Vec3(0));
    MobilizedBody::Pin mobod3(mobod2, Vec3(0), body3Info,Vec3(0));
    MobilizedBody::Pin mobod4(mobod3, Vec3(0), body3Info,Vec3(0));
    MobilizedBody::Free mobod5(matter.Ground(), Vec3(0), body3Info,Vec3(0));

    State s = mbs.realizeTopology();

    // Make sure defaults made it to the state.
    SimTK_TEST(gravity.getBodyIsExcluded(s,MobodIndex(0)));
    for (MobodIndex i(1); i < matter.getNumBodies(); ++i)
        SimTK_TEST(!gravity.getBodyIsExcluded(s,i));
    SimTK_TEST(gravity.getGravityVector(s)==Vec3(0,0,-49));
    SimTK_TEST(gravity.getDownDirection(s)==UnitVec3(0,0,-1));
    SimTK_TEST(gravity.getMagnitude(s)==49);
    SimTK_TEST(gravity.getZeroHeight(s)==-5);

    // Now make some changes in the state.

    // Shouldn't be able to include ground.
    gravity.setBodyIsExcluded(s, MobodIndex(0), false);
    SimTK_TEST(gravity.getBodyIsExcluded(s,MobodIndex(0)));

    gravity.setBodyIsExcluded(s, mobod3, true);
    SimTK_TEST(gravity.getBodyIsExcluded(s,mobod3));
    SimTK_TEST(!gravity.getBodyIsExcluded(s,mobod2));
    SimTK_TEST(!gravity.getBodyIsExcluded(s,mobod4));

    // That shouldn't have changed the default.
    SimTK_TEST(!gravity.getDefaultBodyIsExcluded(mobod3));

    // When making changes in the state we are restricted to bodies that
    // actually exist.
    SimTK_TEST_MUST_THROW(gravity.setBodyIsExcluded(s, MobodIndex(13), true));

    gravity.setGravityVector(s, Vec3(1,2,3));
    SimTK_TEST_EQ(gravity.getGravityVector(s), Vec3(1,2,3));
    SimTK_TEST_EQ(gravity.getDownDirection(s), UnitVec3(Vec3(1,2,3)));
    SimTK_TEST_EQ(gravity.getMagnitude(s), Vec3(1,2,3).norm());

    gravity.setDownDirection(s, UnitVec3(9,10,11));
    SimTK_TEST_EQ(gravity.getDownDirection(s), UnitVec3(9,10,11));
    SimTK_TEST_EQ(gravity.getMagnitude(s), Vec3(1,2,3).norm()); // no change

    SimTK_TEST_MUST_THROW(gravity.setMagnitude(s, -5)); // must be >= 0
    gravity.setMagnitude(s, 5);
    SimTK_TEST_EQ(gravity.getDownDirection(s), UnitVec3(9,10,11)); // no change
    SimTK_TEST(gravity.getMagnitude(s)==5);

    // Changing gravity vector to zero should leave direction unchanged.
    gravity.setGravityVector(s, Vec3(0));
    SimTK_TEST(gravity.getGravityVector(s)==Vec3(0));
    SimTK_TEST_EQ(gravity.getDownDirection(s), UnitVec3(9,10,11));
    SimTK_TEST(gravity.getMagnitude(s)==0);

    gravity.setZeroHeight(s, 1.25);
    SimTK_TEST(gravity.getZeroHeight(s)==1.25);

}

static void testForces() {
    MultibodySystem         mbs;
    SimbodyMatterSubsystem  matter(mbs);
    GeneralForceSubsystem   forces(mbs);
    Force::DiscreteForces   discrete(forces, matter);
    #ifdef USE_VISUALIZER
    Visualizer              viz(mbs);
    viz.setSystemUpDirection(ZAxis);
    viz.setCameraTransform(
        Transform(Rotation(BodyRotationSequence,Pi/2,XAxis,Pi/8,YAxis),
        Vec3(5,-8,2)));
    #endif

    Force::Gravity gravity(forces, matter, -ZAxis, 50, 17);

    matter.Ground().addBodyDecoration(Vec3(0,0,.05),
        DecorativeFrame(2).setColor(Green));

    MobilizedBody::Weld mobod1(matter.Ground(),Vec3(0),
                               body1Info, Vec3(0));

    const Rotation ZtoY(-Pi/2, XAxis);
    MobilizedBody::Pin mobod2(mobod1, Transform(ZtoY,2*Centroid),
                              body2Info, Transform(ZtoY,Vec3(0)));
    MobilizedBody::Pin mobod3(mobod2, Transform(ZtoY,2*Centroid),
                              body3Info, Transform(ZtoY,Vec3(0)));

    State state = mbs.realizeTopology();
    state.updQ() = Test::randVector(state.getNQ());
    state.updU() = Test::randVector(state.getNU());

    mbs.realize(state, Stage::Dynamics);
    #ifdef USE_VISUALIZER
    viz.report(state);
    #endif


    const Vector_<SpatialVec>& bodyForces =
        mbs.getRigidBodyForces(state, Stage::Dynamics);

    // Mobility forces should all be zero.
    const Vector& mobilityForces =
        mbs.getMobilityForces(state, Stage::Dynamics);
    for (int i=0; i < state.getNU(); ++i)
        SimTK_TEST(mobilityForces[i] == 0);

    // Calculate body forces and torques and verify that they are correct.
    // (These are about the body origin, not the COM.)
    SimTK_TEST(bodyForces[0] == SpatialVec(Vec3(0))); // Ground
    Real g = gravity.getMagnitude(state);
    UnitVec3 d = gravity.getDownDirection(state);
    Real h0 = gravity.getZeroHeight(state);
    Real pe = 0;
    for (MobodIndex i(1); i < matter.getNumBodies(); ++i) {
        const Mobod& mobod = matter.getMobilizedBody(i);
        const Real m = mobod.getBodyMass(state);
        const Vec3 p_BC = mobod.getBodyMassCenterStation(state);
        const Vec3 p_BC_G = mobod.expressVectorInGroundFrame(state,p_BC);
        const Vec3 F = m*g*d;
        const Vec3 M = p_BC_G % F;
        SimTK_TEST_EQ(bodyForces[i], SpatialVec(M,F));

        SimTK_TEST_EQ(gravity.getBodyForce(state, i), SpatialVec(M,F));

        const Real h =
            -dot(mobod.findStationLocationInGround(state,p_BC),d)-h0;
        pe += m*g*h;
    }

    SimTK_TEST_EQ(gravity.getPotentialEnergy(state), pe);
    SimTK_TEST_EQ(mbs.calcPotentialEnergy(state), pe);

    // Change zero height and verify that the potential energy changes.
    h0 = 10.3;
    gravity.setZeroHeight(state, h0);
    mbs.realize(state, Stage::Dynamics);
    pe = 0;
    for (MobodIndex i(1); i < matter.getNumBodies(); ++i) {
        const Mobod& mobod = matter.getMobilizedBody(i);
        const Real m = mobod.getBodyMass(state);
        const Vec3 p_BC = mobod.getBodyMassCenterStation(state);
        const Real h =
            -dot(mobod.findStationLocationInGround(state,p_BC),d)-h0;
        pe += m*g*h;
    }
    SimTK_TEST_EQ(gravity.getPotentialEnergy(state), pe);
    SimTK_TEST_EQ(mbs.calcPotentialEnergy(state), pe);

    // Turn off a body and make sure it doesn't see gravity after that, and
    // that the other bodies are unchanged.
    gravity.setBodyIsExcluded(state, mobod2, true);
    mbs.realize(state, Stage::Dynamics);

    SimTK_TEST(bodyForces[0] == SpatialVec(Vec3(0))); // Ground
    pe = 0;
    for (MobodIndex i(1); i < matter.getNumBodies(); ++i) {
        const Mobod& mobod = matter.getMobilizedBody(i);
        const Real m = mobod.getBodyMass(state);
        const Vec3 p_BC = mobod.getBodyMassCenterStation(state);
        const Vec3 p_BC_G = mobod.expressVectorInGroundFrame(state,p_BC);
        const Vec3 F = m*g*d;
        const Vec3 M = p_BC_G % F;

        if (i != mobod2.getMobilizedBodyIndex()) {
            SimTK_TEST_EQ(bodyForces[i], SpatialVec(M,F));
            SimTK_TEST_EQ(gravity.getBodyForce(state, i), SpatialVec(M,F));
            const Real h =
                -dot(mobod.findStationLocationInGround(state,p_BC),d)-h0;
            pe += m*g*h;
        } else {
            SimTK_TEST(bodyForces[i]==SpatialVec(Vec3(0)));
            SimTK_TEST(gravity.getBodyForce(state, i)==SpatialVec(Vec3(0)));
        }
    }
    SimTK_TEST_EQ(gravity.getPotentialEnergy(state), pe);
    SimTK_TEST_EQ(mbs.calcPotentialEnergy(state), pe);

    // Test caching.
    const long long nevals1=gravity.getNumEvaluations();
    SimTK_TEST(gravity.isForceCacheValid(state));
    // This should not require re-evaluation.
    SimTK_TEST(gravity.getBodyForces(state)[0] == SpatialVec(Vec3(0)));
    SimTK_TEST(gravity.getNumEvaluations()==nevals1);

    gravity.invalidateForceCache(state);
    // Force re-evaluation.
    SimTK_TEST(gravity.getBodyForces(state)[0] == SpatialVec(Vec3(0)));
    const long long nevals2=gravity.getNumEvaluations();
    SimTK_TEST(nevals2==nevals1+1);

    state.invalidateAllCacheAtOrAbove(Stage::Velocity); // shouldn't invalidate
    SimTK_TEST(gravity.getBodyForces(state)[0] == SpatialVec(Vec3(0)));
    SimTK_TEST(gravity.getNumEvaluations()==nevals2);
    mbs.realize(state, Stage::Dynamics); // shouldn't reevaluate
    SimTK_TEST(gravity.getNumEvaluations()==nevals2);

    // Bring mobod2 back in. This should only invalidate Dynamics stage, but
    // should nevertheless force recomputation of gravity.
    gravity.setBodyIsExcluded(state, mobod2, false);
    SimTK_TEST(state.getSystemStage() == Stage(Stage::Dynamics-1));
    SimTK_TEST(gravity.getBodyForces(state)[0] == SpatialVec(Vec3(0)));
    const long long nevals3=gravity.getNumEvaluations();
    SimTK_TEST(nevals3==nevals2+1);

    // Make sure that setting gravity to zero works properly -- it is a
    // special case since the zeroes are precalculated. This should not
    // require a gravity evaluation.
    gravity.setMagnitude(state, 0);
    SimTK_TEST(!gravity.isForceCacheValid(state));
    for (int i=0; i < matter.getNumBodies(); ++i)
        SimTK_TEST(gravity.getBodyForces(state)[i] == SpatialVec(Vec3(0)));
    SimTK_TEST(gravity.getNumEvaluations()==nevals3);

    // Setting to non-zero should invalidate, and then require just a single
    // evaluation to respond to multiple calls.
    gravity.setMagnitude(state, 9.8);
    SimTK_TEST(!gravity.isForceCacheValid(state));
    for (int i=1; i < matter.getNumBodies(); ++i)
        SimTK_TEST(gravity.getBodyForces(state)[i] != SpatialVec(Vec3(0)));
    const long long nevals4=gravity.getNumEvaluations();
    SimTK_TEST(nevals4==nevals3+1);

    // Turn off velocities for gravity compensation test to eliminate Coriolis
    // foces (shouldn't invalidate gravity forces).
    state.updU() = 0;
    mbs.realize(state, Stage::Acceleration);
    const Vector_<SpatialVec>& gfrc = gravity.getBodyForces(state);
    Vector f;
    matter.multiplyBySystemJacobianTranspose(state, gfrc, f);
    discrete.setAllMobilityForces(state, -f);
    mbs.realize(state, Stage::Acceleration);
    SimTK_TEST_EQ(state.getUDot(), Vector(state.getNU(), Real(0)));
    SimTK_TEST(gravity.getNumEvaluations()==nevals4); // all for free?

    // Sneaking a zero in by vector should behave just like setting the
    // magnitude to zero.
    gravity.setGravityVector(state, Vec3(0));
    SimTK_TEST(!gravity.isForceCacheValid(state));
    for (int i=0; i < matter.getNumBodies(); ++i)
        SimTK_TEST(gravity.getBodyForces(state)[i] == SpatialVec(Vec3(0)));
    SimTK_TEST(gravity.getNumEvaluations()==nevals4); // no eval needed
}


//==============================================================================
//                                   MAIN
//==============================================================================
int main() {
    // Add artwork.
    DecorativeBrick drawCube(Cube); drawCube.setOpacity(0.5).setColor(Gray);
    body1Info.addDecoration(Centroid, drawCube);
    body2Info.addDecoration(Centroid, drawCube);
    body3Info.addDecoration(Centroid, drawCube);

    SimTK_START_TEST("TestGravity");
        SimTK_SUBTEST(testConstruction);
        SimTK_SUBTEST(testParameters);
        SimTK_SUBTEST(testForces);
    SimTK_END_TEST();
}
