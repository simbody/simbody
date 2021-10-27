/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2021 the Authors.                                   *
 * Authors: Frank C. Anderson                                                 *
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

#include "SimTKsimbody.h"

#include <iostream>
using std::cout;
using std::endl;

using namespace SimTK;


//=============================================================================
// MAIN
//=============================================================================
int main() {
    SimTK_START_TEST("TestExponentialSpring");

    SimTK_SUBTEST(testInitialization);
    SimTK_SUBTEST(testBouncingBlockNoDampingNoFriction);
    SimTK_SUBTEST(testSlidingBlockNoDampingNoFriction);
    SimTK_SUBTEST(testSpinningBlockNoDampingNoFriction);
    SimTK_SUBTEST(testTumblingBlockNoDampingNoFriction);
    SimTK_SUBTEST(testSlidingBlockNoDampingWithFriction);
    SimTK_SUBTEST(testSlidingBlockWithDampingWithFriction);
    SimTK_SUBTEST(testHopping3DPendulumWithDampingWithFriction);
    SimTK_SUBTEST(testMultipleContactPlanes);

    SimTK_END_TEST();
}


//=============================================================================
// SUB TESTS
//=============================================================================

//_____________________________________________________________________________
// Test things that happen before integration.
// 1. Setting and getting parameters (ExponentialSpringParameters)
// 2. Construction with default and non-default parameters
// 3. Realization, getting and setting μₛ and μₖ, resetting spring zeros
void testInitialization() {

    //----------------------------------
    // 1. Setting and getting parameters
    //----------------------------------
    // Test Equality Operator
    ExponentialSpringParameters params, paramsDef;
    SimTK_TEST(params == paramsDef);

    // Get the default parameters
    Real d0Def, d1Def, d2Def;
    paramsDef.getShapeParameters(d0Def, d1Def, d2Def);
    Real kvNormDef = paramsDef.getNormalViscosity();
    Real kpDef = paramsDef.getElasticity();
    Real kvDef = paramsDef.getViscosity();
    Real tauDef = paramsDef.getSlidingTimeConstant();
    Real vSettleDef = paramsDef.getSettleVelocity();
 
    // Make up non-default parameters
    Real delta = 0.1;
    Real d0 = d0Def + delta;
    Real d1 = d1Def + delta;
    Real d2 = d2Def + delta;
    Real kvNorm = kvNormDef + delta;
    Real kp = kpDef + delta;
    Real kv = kvDef + delta;
    Real tau = tauDef + delta;
    Real vSettle = vSettleDef + delta;

    // Test Equality Operator again by setting non-default parameters on
    // the non-default params object.
    params.setShapeParameters(d0, d1, d2);
    SimTK_TEST(params != paramsDef);
    params.setNormalViscosity(kvNorm);
    SimTK_TEST(params != paramsDef);
    params.setElasticity(kp);
    SimTK_TEST(params != paramsDef);
    params.setViscosity(kv);
    SimTK_TEST(params != paramsDef);
    params.setSlidingTimeConstant(tau);
    SimTK_TEST(params != paramsDef);
    params.setSettleVelocity(vSettle);
    SimTK_TEST(params != paramsDef);
    
    // Test the set methods by checking the individual member variables.
    Real tol = 1.0e-12 * delta;
    params.getShapeParameters(d0, d1, d2);
    SimTK_TEST_EQ(d0 - d0Def, delta);
    SimTK_TEST_EQ(d1 - d1Def, delta);
    SimTK_TEST_EQ_TOL(d2 - d2Def, delta, tol);
    SimTK_TEST_EQ(params.getNormalViscosity() - kvNormDef, delta);
    SimTK_TEST_EQ_TOL(params.getElasticity() - kpDef, delta, tol);
    SimTK_TEST_EQ(params.getViscosity() - kvDef, delta);
    SimTK_TEST_EQ(params.getSlidingTimeConstant() - tauDef, delta);
    SimTK_TEST_EQ(params.getSettleVelocity() - vSettleDef, delta);

    // Test setting viscosity for critical damping,
    // with a default mass (1.0 kg)
    Real kp1 = 100000.0;
    Real kv1 = 2.0 * sqrt(kp1);
    params.setElasticityAndViscosityForCriticalDamping(kp1);
    SimTK_TEST_EQ(params.getViscosity(), kv1);
    // with a non-default mass
    Real mass = 20.0;
    kv1 = 2.0 * sqrt(kp1 * mass);
    params.setElasticityAndViscosityForCriticalDamping(kp1, mass);
    SimTK_TEST_EQ(params.getViscosity(), kv1);

    // Return the non-default params to the original non-default values
    params.setElasticity(kp);
    params.setViscosity(kv);

    //------------------------------------------------------------
    // 2. Test construction with default an non-default parameters
    //------------------------------------------------------------
    // Construct the system and subsystems.
    MultibodySystem system; system.setUseUniformBackground(true);
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0., -9.8, 0.));
    Body::Rigid bodyProps(MassProperties(10.0, Vec3(0), Inertia(1)));

    // Create a rotation and translation for the contact plane.
    // Origin is shifted +1 on the x axis.
    // Plane is titled by +45 deg about the z axis.
    Real angle = 0.25 * SimTK::Pi;
    Rotation planeTilt(angle, ZAxis);
    Vec3 planeOrigin(1., 0., 0.);

    // Construct several exponential springs.
    Transform plane(planeTilt, planeOrigin);
    MobilizedBody::Free body(matter.Ground(), bodyProps);
    Vec3 station(0.1, -0.1, 0.1);
    Real mus = 0.7, muk = 0.5;
    ExponentialSpringForce
        sprDef0(system, plane, body, station, mus, muk);
    ExponentialSpringForce
        sprDef1(system, plane, body, station, mus, muk, paramsDef);
    ExponentialSpringForce
        spr(system, plane, body, station, mus, muk, params);

    // Test that the non-default parameter values were set properly.
    SimTK_TEST(sprDef0.getParameters() == sprDef1.getParameters());
    SimTK_TEST(spr.getParameters() != sprDef0.getParameters());
    const ExponentialSpringParameters& p = spr.getParameters();
    p.getShapeParameters(d0, d1, d2);
    SimTK_TEST_EQ(d0 - d0Def, delta);
    SimTK_TEST_EQ(d1 - d1Def, delta);
    SimTK_TEST_EQ_TOL(d2 - d2Def, delta, tol);
    SimTK_TEST_EQ(p.getNormalViscosity() - kvNormDef, delta);
    SimTK_TEST_EQ_TOL(p.getElasticity() - kpDef, delta, tol);
    SimTK_TEST_EQ(p.getViscosity() - kvDef, delta);
    SimTK_TEST_EQ(p.getSlidingTimeConstant() - tauDef, delta);
    SimTK_TEST_EQ(p.getSettleVelocity() - vSettleDef, delta);

    // Test getting contact plane, body, and station
    SimTK_TEST(spr.getContactPlane() == plane);
    SimTK_TEST(spr.getBody().isSameMobilizedBody(body));
    SimTK_TEST(spr.getStation() == station);

    //------------------------------------------------------------------------
    // 3. Realization, getting and setting μₛ and μₖ, resetting spring zeros.
    //------------------------------------------------------------------------
    // Realize through Stage::Model (construct the state)
    system.realizeTopology();
    State state = system.getDefaultState();
    system.realizeModel(state);

    // Test getting μₛ and μₖ
    SimTK_TEST(spr.getMuStatic(state) == mus);
    SimTK_TEST(spr.getMuKinetic(state) == muk);

    // Test setting μₛ and μₖ
    // μₖ ≤ μₛ
    Real musNew = mus + 1.0;
    Real mukNew = muk + 1.0;
    spr.setMuStatic(state, musNew);
    spr.setMuKinetic(state, mukNew);
    SimTK_TEST(spr.getMuStatic(state) == musNew);
    SimTK_TEST(spr.getMuKinetic(state) == mukNew);
    // μₖ ≥ μₛ (μₛ should be set to μₖ to maintain μₖ ≤ μₛ)
    mukNew = musNew + 1.0;
    spr.setMuKinetic(state, mukNew);
    SimTK_TEST(spr.getMuStatic(state) == mukNew);
    SimTK_TEST(spr.getMuKinetic(state) == mukNew);
    // μₛ ≤ μₖ (μₖ should be set to μₛ to maintain μₖ ≤ μₛ)
    musNew = musNew - 1.5;
    spr.setMuStatic(state, musNew);
    SimTK_TEST(spr.getMuStatic(state) == musNew);
    SimTK_TEST(spr.getMuKinetic(state) == musNew);
    // μₖ < 0.0 (μₖ should be set to 0.0; μₛ should not change)
    mukNew = -0.5;
    musNew = spr.getMuStatic(state);
    spr.setMuKinetic(state, mukNew);
    SimTK_TEST(spr.getMuStatic(state) == musNew);
    SimTK_TEST(spr.getMuKinetic(state) == 0.0);
    // μₛ < 0.0 (μₛ and μₖ should both be set to 0.0)
    musNew = -0.5;
    mukNew = 0.5;
    spr.setMuKinetic(state, mukNew);
    spr.setMuStatic(state, musNew);
    SimTK_TEST(spr.getMuStatic(state) == 0.0);
    SimTK_TEST(spr.getMuKinetic(state) == 0.0);

    // Test resetting the spring zero
    // Before the reset, the spring zero is the origin of the contact plane,
    // which may be displaced from the Ground origin.
    system.realize(state, Stage::Position);
    Vec3 p0(0.0, 0.0, 0.0);  // expressed in the frame of contact plane
    Vec3 p0_G = plane.shiftFrameStationToBase(p0);
    SimTK_TEST(p0 == spr.getSpringZeroPosition(state, false));
    SimTK_TEST(p0_G == spr.getSpringZeroPosition(state));
    // Reset
    // Now the spring zero should coincide with the projection of the
    // spring station onto the contact plane, which should result in
    // the elastic component of the frictional force being zero.
    spr.resetSpringZero(state);
    // Expected positions after the reset
    // express spring station in the Ground frame
    Vec3 s_G = body.findStationLocationInGround(state, station);
    // express station in the contact plane frame
    Vec3 s = plane.shiftBaseStationToFrame(s_G);
    // project onto the contact plane by setting the y-component to 0.0
    Vec3 p0After = s;  p0After[1] = 0.0;
    // express in Ground frame after the projection
    Vec3 p0After_G = plane.shiftFrameStationToBase(p0After);
    SimTK_TEST(p0After == spr.getSpringZeroPosition(state, false));
    SimTK_TEST(p0After_G == spr.getSpringZeroPosition(state));

    // Test that the elastic component of the friction force is 0.0
    system.realize(state, Stage::Dynamics);
    Vec3 fElastic = spr.getFrictionForceElasticPart(state);
    SimTK_TEST_EQ(fElastic,Vec3(0.0, 0.0, 0.0),1);
}

//_____________________________________________________________________________
// Test, in a simple case, that energy is conserved and that the strain energy
// calculation (i.e., the potential energy calculation) is correct.
// A 10.0 kg body, mobilized by a vertical Slider, is dropped from 1.0 meters.
// One ExponentialSpringForce, with no damping, acts on the body.
// The body should bounce back to the height from which it was dropped (within
// the accuracy of the integrator).
// In addition, when the body is at its lowest point (i.e., when its kinetic
// energy is zero), its total potential energy (gravitational + spring strain
// energy) should equal its initial total potential energy.
void testBouncingBlockNoDampingNoFriction() {

}

//_____________________________________________________________________________
void testSlidingBlockNoDampingNoFriction() {

}

//_____________________________________________________________________________
void testSpinningBlockNoDampingNoFriction() {

}

//_____________________________________________________________________________
void testTumblingBlockNoDampingNoFriction() {

}

//_____________________________________________________________________________
void testSlidingBlockNoDampingWithFriction() {

}

//_____________________________________________________________________________
void testSlidingBlockWithDampingWithFriction() {

}

//_____________________________________________________________________________
// Test that the spring force calculations are correct or at least internally
// consistent. The getters for the data cache are checked in this process.
// A 10 kg body, mobilized by a Free joint, is tossed horizontally from a
// height of 1.0 meter. One ExponentialSpringForce is placed at a Station that
// is displaced by 0.2 meters from the body's center of mass. The simulation
// is long enough for the spring zero to stop sliding. The spring data is
// sampled and checked every 0.1 seconds.
void testHopping3DPendulumWithDampingWithFriction() {

}


//_____________________________________________________________________________
void testMultipleContactPlanes() {

}



//=============================================================================
// Event Reporters and the Vizualizer
//=============================================================================

//_____________________________________________________________________________
// This class impmlements an event reported that stores the System State at
// regular time intervals.
class PeriodicStateStorage : public PeriodicEventReporter {
public:
    PeriodicStateStorage(const MultibodySystem& system, Real reportInterval)
        : PeriodicEventReporter(reportInterval) {
        storage = new Array_<State>;
    }
    ~PeriodicStateStorage() {
        delete storage;
    }
    void handleEvent(const State& state) const override {
        storage->push_back(state);
    }
    void clearStateStorage() {
        storage->clear();
    }
    const Array_<State>* getStateStorage() {
        return storage;
    }
private:
    Array_<State> *storage;
};




//=============================================================================
/** This class impmlements an event reported that identifies the highest
position reached by a body.  If energy is conserved, which should
be the case if there are no damping components of the contact force,
then the height should be constant across multiple bounces. */
class MaxHeightReporter : public TriggeredEventReporter {
public:
    MaxHeightReporter(const MultibodySystem& system, const MobilizedBody& body)
        : TriggeredEventReporter(Stage::Velocity), system(system), body(body) {
        getTriggerInfo().setTriggerOnRisingSignTransition(false);
    }
    Real getValue(const State& state) const {
        Vec3 vel = body.getBodyOriginVelocity(state);
        return vel[1];
    }
    void handleEvent(const State& state) const {
        system.realize(state, Stage::Velocity);
        Vec3 pos = body.getBodyOriginLocation(state);
        Vec3 vel = body.getBodyOriginVelocity(state);
        cout << endl;
        cout << state.getTime() << "\tp = " << pos << "\tv = " << vel << endl;
        cout << endl;
    }
private:
    const MultibodySystem& system;
    const MobilizedBody& body;
};


