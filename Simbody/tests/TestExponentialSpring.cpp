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

// Structure for specifying simulation options
struct SimulationOptions {
    int condition;  // which set of initial conditions
    bool damping;   // whether there is damping
    bool friction;  // whether there is friction
    bool viz;       // whether to create a visualization window
    Real tilt;      // tilt of primary contact plane about z-axis
};

// Declarations
void testInitialization();
void testBlockVerticalBounceNoFricNoDamp();
void simulateBlock(const SimulationOptions &options);
void checkSpringCalculations(MultibodySystem& system, Real acc,
    ExponentialSpringForce& spr, const Array_<State>* stateArray);
void checkConservationOfEnergy(MultibodySystem& system, Real acc,
    ExponentialSpringForce& spr, const Array_<State>* stateArray);


//=============================================================================
// Event Reporters and the Vizualizer
//=============================================================================
//_____________________________________________________________________________
// This class impmlements an event reported that records the System State at
// a regular time interval.
class PeriodicStateRecorder : public PeriodicEventReporter {
public:
    PeriodicStateRecorder(Real reportInterval)
        : PeriodicEventReporter(reportInterval) {
        storage = new Array_<State>;
    }
    ~PeriodicStateRecorder() {
        delete storage;
    }
    void handleEvent(const State& state) const override {
        storage->push_back(state);
    }
    void clearStateArray() {
        storage->clear();
    }
    const Array_<State>* getStateArray() {
        return storage;
    }
private:
    Array_<State>* storage;
};
//_____________________________________________________________________________
// This class implements an event reported that identifies the minimum and
// maximum heights reached by a body and records the State at those heights.
// This reporter is included in the testing to ensure that key events
// during contact events are not missed.  It is posible that the periodic
// reporter (see directly above) could step over such events.
class MinMaxHeightStateRecorder : public TriggeredEventReporter {
public:
    MinMaxHeightStateRecorder(const MultibodySystem& system,
            const MobilizedBody& body) :
            TriggeredEventReporter(Stage::Velocity),
            system(system), body(body) {
        storage = new Array_<State>;
        getTriggerInfo().setTriggerOnRisingSignTransition(true);
    }
    ~MinMaxHeightStateRecorder() {
        delete storage;
    }
    Real getValue(const State& state) const {
        Vec3 vel = body.getBodyOriginVelocity(state);
        return vel[1];
    }
    void handleEvent(const State& state) const {
        system.realize(state, Stage::Velocity);
        Vec3 pos = body.getBodyOriginLocation(state);
        Vec3 vel = body.getBodyOriginVelocity(state);
        //cout << state.getTime() << "\tp = " << pos << "\tv = " << vel << endl;
        storage->push_back(state);
    }
    void clearStateArray() {
        storage->clear();
    }
    const Array_<State>* getStateArray() {
        return storage;
    }

private:
    const MultibodySystem& system;
    const MobilizedBody& body;
    Array_<State>* storage;
};


//=============================================================================
// MAIN
//=============================================================================
int main() {
    SimTK_START_TEST("TestExponentialSpring");

    SimTK_SUBTEST(testInitialization);
    SimTK_SUBTEST(testBlockVerticalBounceNoFricNoDamp);


    SimTK_END_TEST();
}


//=============================================================================
// SUB TESTS
//=============================================================================
//_____________________________________________________________________________
// Test the spring force calculations for a block that bounces up and down
// on the floor without damping or friction.
void testBlockVerticalBounceNoFricNoDamp() {
    // Set Options
    SimulationOptions options;
    options.condition = 1;
    options.friction = false;
    options.damping = false;
    options.viz = false;
    options.tilt = 0.0;

    // Run the simulation
    simulateBlock(options);
}
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
    SimTK_TEST_EQ(fElastic,Vec3(0., 0., 0.));
};


//=============================================================================
// Routines called by Subtests
//=============================================================================
//_____________________________________________________________________________
// Simulate a block bouncing and sliding on the floor. The block has a mass of
// 10.0 kg and is mobilized by a Free joint. An ExponentialSpringForce acts at
// each corner of the block. A second contact plane representing a wall is
// also added to the simulation.
//
// Simulation options specify which set of initial conditions are set, whether
// damping and friction are part of the simulation, and the tilt of the floor
// plane.  
void simulateBlock(const SimulationOptions& options) {
    // Construct the system and basic subsystems.
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0., -9.8, 0.));

    // Construct the mobilized body.
    // The body will slide along the y-axis (up and down)
    // Rotate the parent frame (F) and body frame (M) so that the x axis up.
    Body::Rigid bodyProps(MassProperties(10.0, Vec3(0), Inertia(1)));
    MobilizedBody::Free body(matter.Ground(), bodyProps);

    // Add exponential springs to the block.
    Transform contactPlane(Rotation(0.0, ZAxis), Vec3(0.));  // No change.
    Vec3 station(0., 0., 0.);  // Spring acts at center of mass.
    Real mus = 0.0, muk = 0.0;  // No friction.
    ExponentialSpringParameters params;
    params.setNormalViscosity(0.0);  // No damping in the normal direction.
    params.setViscosity(0.0);  // No damping in the friction spring.
    ExponentialSpringForce
        spr(system, contactPlane, body, station, mus, muk, params);

    // Add reporters and visualization
    PeriodicStateRecorder* periodicRecorder = new PeriodicStateRecorder(0.1);
    system.addEventReporter(periodicRecorder);
    MinMaxHeightStateRecorder* minmaxRecorder =
        new MinMaxHeightStateRecorder(system, body);
    system.addEventReporter(minmaxRecorder);

    // Realize through Stage::Model (construct the State)
    system.realizeTopology();
    State state = system.getDefaultState();
    system.realizeModel(state);

    // Set the initial states
    body.setQToFitTranslation(state, Vec3(0.0, 1.0, 0.0));
    body.setU(state, Vec6(0., 0., 0., 0., 0., 0.));
    spr.resetSpringZero(state);

    // Simulate
    RungeKuttaMersonIntegrator integ(system);
    Real acc = 1.0e-5;
    integ.setAccuracy(acc);
    integ.setMaximumStepSize(0.01);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(5.0);

    // Get the time histories of states
    const Array_<State>* periodicArray = periodicRecorder->getStateArray();
    //cout << "Number of stored state arrays = " << periodicArray->size() << endl;
    const Array_<State>* minmaxArray = minmaxRecorder->getStateArray();
    //cout << "Number of stored state arrays = " << minmaxArray->size() << endl;

    // Check the spring force calculations
    checkSpringCalculations(system, acc, spr, periodicArray);
    checkSpringCalculations(system, acc, spr, minmaxArray);

    // Check energy conservation
    if((options.damping == false) && (options.friction == false))
        checkConservationOfEnergy(system, acc, spr, periodicArray);
    if((options.damping == false) && (options.friction == false))
        checkConservationOfEnergy(system, acc, spr, minmaxArray);
}
//_____________________________________________________________________________
// Check that spring forces are internally consistent for a simulation.
// system       : system being simulated
// acc          : integrator accuracy
// spr          : exponential spring
// stateArray   : time history of states
void checkSpringCalculations(MultibodySystem& system, Real acc,
    ExponentialSpringForce& spr, const Array_<State>* stateArray) {
    // Check for empty state array
    if(stateArray->size() == 0) return;

    // Check that getting data before Stage::Dynamics has been realized
    // throws an exception.
    const State& state = (*stateArray)[0];
    SimTK_TEST_MUST_THROW(spr.getNormalForce(state));

    // Loop through the recorded states
    for(unsigned int i = 0; i < stateArray->size(); ++i) {

        const State& state = (*stateArray)[i];

        // Realize through Stage::Dynamics
        system.realize(state, Stage::Dynamics);

        // Check the normal force when expressed in the contact plane
        Vec3 fyElas = spr.getNormalForceElasticPart(state, false);
        Vec3 fyDamp = spr.getNormalForceDampingPart(state, false);
        Vec3 fy = spr.getNormalForce(state, false);
        // x components should be 0.0
        SimTK_TEST(fyElas[0] == 0.0);
        SimTK_TEST(fyDamp[0] == 0.0);
        SimTK_TEST(fy[0] == 0.0);
        // z components should be 0.0
        SimTK_TEST(fyElas[2] == 0.0);
        SimTK_TEST(fyDamp[2] == 0.0);
        SimTK_TEST(fy[2] == 0.0);
        // sum elastic and damping
        Vec3 fyCalc = fyElas - fyDamp;
        // bounds on fy are enforced in realizeSubsystemDynamicsImpl()
        if(fyCalc[1] < 0.0) fyCalc[1] = 0.0;
        if(fyCalc[1] > 1000000.0) fyCalc[1] = 1000000.0;
        //cout << "fyElas= " << fyElas << "  fyDamp= " << fyDamp << "  fy= " << fy << endl;
        SimTK_TEST_EQ(fyCalc, fy);

        // Check normal force when expressed in Ground
        Vec3 fyElas_G = spr.getNormalForceElasticPart(state);
        Vec3 fyDamp_G = spr.getNormalForceDampingPart(state);
        Vec3 fy_G = spr.getNormalForce(state);
        Vec3 fyCalc_G = fyElas_G + fyDamp_G;
        // bounds on fy are enforced in realizeSubsystemDynamicsImpl()
        if(fyCalc_G[1] < 0.0) fyCalc_G[1] = 0.0;
        if(fyCalc_G[1] > 1000000.0) fyCalc_G[1] = 1000000.0;
        //cout << "fyElas= " << fyElas << "  fyDamp= " << fyDamp << "  fy= " << fy << endl;
        SimTK_TEST_EQ(fyCalc_G, fy);

        // Check magnitude of fy is same in Ground and ContacPlane
        SimTK_TEST_EQ(fy_G.norm(), fy.norm());

        // Check 0.0 ≤ sliding state ≤ 1.0
        Real sliding = spr.getSliding(state);
        //cout << "sliding= " << sliding << endl;
        SimTK_TEST(sliding < 1.0 + acc);
        SimTK_TEST(sliding > 0.0 - acc);

        // Check μₖ ≤ μ ≤ μₛ
        Real mus = spr.getMuStatic(state);
        Real muk = spr.getMuKinetic(state);
        Real mu = spr.getMu(state);
        SimTK_TEST(muk <= mu);
        SimTK_TEST(mu <= mus);

        // Check friction limit ≤ μ*fy
        Vec3 mu_fy = mu * fy_G;
        Real fricLimit = spr.getFrictionForceLimit(state);
        SimTK_TEST(fricLimit <= mu_fy.norm());

        // Check friction force expressed in the contact plane
        Vec3 fricElas = spr.getFrictionForceElasticPart(state, false);
        Vec3 fricDamp = spr.getFrictionForceDampingPart(state, false);
        Vec3 fric = spr.getFrictionForce(state, false);
        // y component should be 0.0
        SimTK_TEST(fricElas[1] == 0.0);
        SimTK_TEST(fricDamp[1] == 0.0);
        SimTK_TEST(fric[1] == 0.0);
        // sum elastic and damping
        Vec3 fricCalc = fricElas + fricDamp;
        //cout << "fricElas= " << fricElas << "  fricDamp= " << fricDamp << "  fric= " << fric << endl;
        SimTK_TEST_EQ(fricCalc, fric);

        // Check friction force expressed in Ground
        Vec3 fricElas_G = spr.getFrictionForceElasticPart(state);
        Vec3 fricDamp_G = spr.getFrictionForceDampingPart(state);
        Vec3 fric_G = spr.getFrictionForce(state);
        Vec3 fricCalc_G = fricElas_G + fricDamp_G;
        SimTK_TEST_EQ(fricCalc_G, fric_G);

        // Check magnitude of fric is same in Ground and Contact Plane
        SimTK_TEST_EQ(fric_G.norm(), fric.norm());

        // Check |friction| ≤ friction limit
        SimTK_TEST(fric.norm() <= fricLimit);

        // Check total force in Contact Plane
        Vec3 f = spr.getForce(state, false);
        Vec3 fCalc = fyCalc + fricCalc;
        SimTK_TEST_EQ(fCalc, f);

        // Check total force in Ground
        Vec3 f_G = spr.getForce(state);
        Vec3 fCalc_G = fyCalc_G + fricCalc_G;
        SimTK_TEST_EQ(fCalc_G, f_G);

        // Check magnitude of total force is same in Ground and Contact Plane
        SimTK_TEST_EQ(f_G.norm(), f.norm());

        // Check that the point at which the spring force is applied is equal
        // to the spring station when expressed in Ground.
        Vec3 p_G = spr.getStationPosition(state);
        Vec3 station = spr.getStation();
        Vec3 station_G =
            spr.getBody().findStationLocationInGround(state, station);
        SimTK_TEST_EQ(station_G, p_G);

        // Check that the spring zero lies in the Contact Plane
        Vec3 p0 = spr.getSpringZeroPosition(state, false);
        SimTK_TEST(p0[1] == 0.0);
    }
}

//_____________________________________________________________________________
// Check that spring forces are internally consistent for a simulation.
// system       : system being simulated
// acc          : integrator accuracy
// spr          : exponential spring
// stateArray   : time history of states
void checkConservationOfEnergy(MultibodySystem& system, Real acc,
    ExponentialSpringForce& spr, const Array_<State>* stateArray) {
    // Check for empty state array
    if(stateArray->size() <= 1) return;

    // Calculate the initial energy of the system
    const State& state0 = (*stateArray)[0];
    system.realize(state0, Stage::Dynamics);
    Real energy0 = system.calcEnergy(state0);

    // Loop through the recorded states
    for(unsigned int i = 0; i < stateArray->size(); ++i) {

        const State& state = (*stateArray)[i];

        // Realize through Stage::Dynamics
        system.realize(state, Stage::Dynamics);

        // Check the energy
        Real energy = system.calcEnergy(state);
        //cout << "energy= " << energy << "  energy0= " << energy0 << endl;
        SimTK_TEST_EQ_TOL(energy, energy0, acc*energy0);
    }
}


