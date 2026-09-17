/* --------------------------------------------------------------------------*
 *                       Simbody(tm): SimTKcommon                            *
 * --------------------------------------------------------------------------*
 * This is part of the SimTK biosimulation toolkit originating from          *
 * Simbios, the NIH National Center for Physics-Based Simulation of          *
 * Biological Structures at Stanford, funded under the NIH Roadmap for       *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody. *
 *                                                                           *
 * Portions copyright (c) 2021-22 the Authors.                               *
 * Authors: Frank C. Anderson                                                *
 * Contributors:                                                             *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * Unless required by applicable law or agreed to in writing, software       *
 * distributed under the License is distributed on an "AS IS" BASIS,         *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  *
 * See the License for the specific language governing permissions and       *
 * limitations under the License.                                            *
 * ------------------------------------------------------------------------- */

#include "SimTKsimbody.h"
#include <iostream>
#include <cstdlib>
#include <memory>
using std::cout;
using std::endl;
using std::unique_ptr;

using namespace SimTK;

// Item for the visualization menu.
static const int RunMenuId = 3, HelpMenuId = 7;
static const int GoItem = 1, ReplayItem = 2, QuitItem = 3;

// Global flag for turning off visualization.
// This flag overrides locally set behavior when it is set to true.
// BEFORE COMMITTING changes, make sure to set VizOff = true;
static const bool VizOff = true;

// Structure for specifying simulation options
struct SimulationOptions {
    int condition;  // which set of initial conditions
    bool damping;   // whether there is damping
    bool friction;  // whether there is friction
    bool viz;       // whether to create a visualization window
    Real tilt;      // tilt of primary contact plane about z-axis
    Real tf;        // final time of the simulation
};

// Declarations
void testInitialization();
void testBlockVerticalBounceNoDampNoFric();
void testBlockVerticalBounceNoDampNoFric();
void testBlockVerticalBounceNoDampWithFric();
void testBlockVerticalBounceWithDampNoFric();
void testBlockVerticalBounceWithDampWithFric();

void testBlockSlideNoDampNoFric();
void testBlockSlideWithDampWithFric();

void testBlockSpinNoDampNoFric();
void testBlockSpinWithDampWithFric();
void testBlockSpinWithDampWithFricTilt15();

void testBlockSpinAndSlideNoDampNoFric();
void testBlockSpinAndSlideNoDampWithFric();
void testBlockSpinAndSlideWithDampNoFric();
void testBlockSpinAndSlideWithDampWithFric();

void testBlockSpinLikeTopNoDampNoFric();
void testBlockSpinLikeTopWithDampWithFric();

void testBlockTumbleNoDampNoFric();
void testBlockTumbleWithDampWithFric();

// Lower level routines called by the sub-tests
void simulateBlock(const SimulationOptions &options);
void checkSpringCalculations(const SimulationOptions& options,
    MultibodySystem& system, Real acc, ExponentialSpringForce& spr,
    const Array_<State>* stateArray);
void checkConservationOfEnergy(MultibodySystem& system, Real acc,
    ExponentialSpringForce& spr, const Array_<State>* stateArray);


//=============================================================================
// Reporters
//=============================================================================
//_____________________________________________________________________________
// This class implements an event reporter that records the System State at
// a regular time interval.
class PeriodicStateRecorder : public PeriodicEventReporter {
public:
    PeriodicStateRecorder(Real reportInterval)
        : PeriodicEventReporter(reportInterval) {
        storage = unique_ptr<Array_<State>>(new Array_<State>());
    }
    void handleEvent(const State& state) const override {
        storage->emplace_back(state);
    }
    void clearStateArray() {
        storage->clear();
    }
    const Array_<State>* getStateArray() {
        return storage.get();
    }
private:
    unique_ptr<Array_<State>> storage;
};
//_____________________________________________________________________________
// This class implements an event reporter that identifies the minimum and
// maximum heights reached by a body and records the State at those heights.
// This reporter is included in the testing to ensure that key events
// during contact events are not missed.  It is possible that the periodic
// reporter (see directly above) could step over such events.
class MinMaxHeightStateRecorder : public TriggeredEventReporter {
public:
    MinMaxHeightStateRecorder(const MultibodySystem& system,
            const MobilizedBody& body) :
            TriggeredEventReporter(Stage::Velocity),
            system(system), body(body) {
        storage = unique_ptr<Array_<State>>(new Array_<State>());
        getTriggerInfo().setTriggerOnRisingSignTransition(true);
    }
    ~MinMaxHeightStateRecorder() {
    }
    Real getValue(const State& state) const {
        Vec3 vel = body.getBodyOriginVelocity(state);
        return vel[1];
    }
    void handleEvent(const State& state) const {
        system.realize(state, Stage::Velocity);
        Vec3 pos = body.getBodyOriginLocation(state);
        Vec3 vel = body.getBodyOriginVelocity(state);
        storage->emplace_back(state);
    }
    void clearStateArray() {
        storage->clear();
    }
    const Array_<State>* getStateArray() {
        return storage.get();
    }

private:
    const MultibodySystem& system;
    const MobilizedBody& body;
    unique_ptr<Array_<State>> storage;
};


//=============================================================================
// MAIN
//=============================================================================
int main() {
    SimTK_START_TEST("TestExponentialSpring");

    SimTK_SUBTEST(testInitialization);

    SimTK_SUBTEST(testBlockVerticalBounceNoDampNoFric);
    SimTK_SUBTEST(testBlockVerticalBounceNoDampWithFric);
    SimTK_SUBTEST(testBlockVerticalBounceWithDampNoFric);
    SimTK_SUBTEST(testBlockVerticalBounceWithDampWithFric);

    SimTK_SUBTEST(testBlockSlideNoDampNoFric);
    SimTK_SUBTEST(testBlockSlideWithDampWithFric);

    SimTK_SUBTEST(testBlockSpinNoDampNoFric);
    SimTK_SUBTEST(testBlockSpinWithDampWithFric);
    SimTK_SUBTEST(testBlockSpinWithDampWithFricTilt15);

    SimTK_SUBTEST(testBlockSpinAndSlideNoDampNoFric);
    SimTK_SUBTEST(testBlockSpinAndSlideNoDampWithFric);
    SimTK_SUBTEST(testBlockSpinAndSlideWithDampNoFric);
    SimTK_SUBTEST(testBlockSpinAndSlideWithDampWithFric);

    SimTK_SUBTEST(testBlockSpinLikeTopNoDampNoFric);
    SimTK_SUBTEST(testBlockSpinLikeTopWithDampWithFric);

    SimTK_SUBTEST(testBlockTumbleNoDampNoFric);
    SimTK_SUBTEST(testBlockTumbleWithDampWithFric);

    SimTK_END_TEST();
}


//=============================================================================
// SUB TESTS
//=============================================================================
//_____________________________________________________________________________
// Drop, no damp, no fric
void testBlockVerticalBounceNoDampNoFric() {
    // Set Options
    SimulationOptions options;
    options.condition = 1;
    options.damping = false;
    options.friction = false;
    options.viz = true;
    options.tilt = 0.0;
    options.tf = 2.4;

    // Check global viz flag
    if(VizOff) options.viz = false;

    // Run the simulation
    simulateBlock(options);
}
//_____________________________________________________________________________
// Drop, no damp, with fric
void testBlockVerticalBounceNoDampWithFric() {
    // Set Options
    SimulationOptions options;
    options.condition = 1;
    options.damping = false;
    options.friction = true;
    options.viz = true;
    options.tilt = 0.0;
    options.tf = 2.4;

    // Check global viz flag
    if(VizOff) options.viz = false;

    // Run the simulation
    simulateBlock(options);
}
//_____________________________________________________________________________
// Drop, with damp, no fric
void testBlockVerticalBounceWithDampNoFric() {
    // Set Options
    SimulationOptions options;
    options.condition = 1;
    options.damping = true;
    options.friction = false;
    options.viz = true;
    options.tilt = 0.0;
    options.tf = 3.0;

    // Check global viz flag
    if(VizOff) options.viz = false;

    // Run the simulation
    simulateBlock(options);
}
//_____________________________________________________________________________
// Drop, with damp, with fric.
void testBlockVerticalBounceWithDampWithFric() {
    // Set Options
    SimulationOptions options;
    options.condition = 1;
    options.damping = true;
    options.friction = true;
    options.viz = true;
    options.tilt = 0.0;
    options.tf = 3.0;

    // Check global viz flag
    if(VizOff) options.viz = false;

    // Run the simulation
    simulateBlock(options);
}

//_____________________________________________________________________________
// Slide, no damp, no fric
void testBlockSlideNoDampNoFric() {
    // Set Options
    SimulationOptions options;
    options.condition = 2;
    options.damping = false;
    options.friction = false;
    options.viz = true;
    options.tilt = 0.0;
    options.tf = 1.0;

    // Check global viz flag
    if(VizOff) options.viz = false;

    // Run the simulation
    simulateBlock(options);
}
//_____________________________________________________________________________
// Slide, with damp, with fric
void testBlockSlideWithDampWithFric() {
    // Set Options
    SimulationOptions options;
    options.condition = 2;
    options.damping = true;
    options.friction = true;
    options.viz = true;
    options.tilt = 0.0;
    options.tf = 3.0;

    // Check global viz flag
    if(VizOff) options.viz = false;

    // Run the simulation
    simulateBlock(options);
}

//_____________________________________________________________________________
// Spin, no damp, no fric
void testBlockSpinNoDampNoFric() {
    // Set Options
    SimulationOptions options;
    options.condition = 3;
    options.damping = false;
    options.friction = false;
    options.viz = true;
    options.tilt = 0.0;
    options.tf = 2.0;

    // Check global viz flag
    if(VizOff) options.viz = false;

    // Run the simulation
    simulateBlock(options);
}
//_____________________________________________________________________________
// Spin, with damp, with fric
void testBlockSpinWithDampWithFric() {
    // Set Options
    SimulationOptions options;
    options.condition = 3;
    options.damping = true;
    options.friction = true;
    options.viz = true;
    options.tilt = 0.0;
    options.tf = 3.0;

    // Check global viz flag
    if(VizOff) options.viz = false;

    // Run the simulation
    simulateBlock(options);
}
//_____________________________________________________________________________
// Spin, with damp, with fric, tilt = -15 deg about Z
void testBlockSpinWithDampWithFricTilt15() {
    // Set Options
    SimulationOptions options;
    options.condition = 3;
    options.damping = true;
    options.friction = true;
    options.viz = true;
    options.tilt = -15.0;
    options.tf = 8.0;

    // Check global viz flag
    if(VizOff) options.viz = false;

    // Run the simulation
    simulateBlock(options);
}

//_____________________________________________________________________________
// Spin and slide, no damp, no fric
void testBlockSpinAndSlideNoDampNoFric() {
    // Set Options
    SimulationOptions options;
    options.condition = 4;
    options.damping = false;
    options.friction = false;
    options.viz = true;
    options.tilt = 0.0;
    options.tf = 0.8;

    // Check global viz flag
    if(VizOff) options.viz = false;

    // Run the simulation
    simulateBlock(options);
}
//_____________________________________________________________________________
// Spin and slide, no damp, with fric
void testBlockSpinAndSlideNoDampWithFric() {
    // Set Options
    SimulationOptions options;
    options.condition = 4;
    options.damping = false;
    options.friction = true;
    options.viz = true;
    options.tilt = 0.0;
    options.tf = 1.6;

    // Check global viz flag
    if(VizOff) options.viz = false;

    // Run the simulation
    simulateBlock(options);
}
//_____________________________________________________________________________
// Spin and slide, with damp, no fric
void testBlockSpinAndSlideWithDampNoFric() {
    // Set Options
    SimulationOptions options;
    options.condition = 4;
    options.damping = true;
    options.friction = false;
    options.viz = true;
    options.tilt = 0.0;
    options.tf = 1.6;

    // Check global viz flag
    if(VizOff) options.viz = false;

    // Run the simulation
    simulateBlock(options);
}
//_____________________________________________________________________________
// Spin and slide, with damp, with fric
void testBlockSpinAndSlideWithDampWithFric() {
    // Set Options
    SimulationOptions options;
    options.condition = 4;
    options.damping = true;
    options.friction = true;
    options.viz = true;
    options.tilt = 0.0;
    options.tf = 5.5;

    // Check global viz flag
    if(VizOff) options.viz = false;

    // Run the simulation
    simulateBlock(options);
}

//_____________________________________________________________________________
// Spin like top, no damp, no fric
void testBlockSpinLikeTopNoDampNoFric() {
    // Set Options
    SimulationOptions options;
    options.condition = 5;
    options.damping = false;
    options.friction = false;
    options.viz = true;
    options.tilt = 0.0;
    options.tf = 2.0;

    // Check global viz flag
    if(VizOff) options.viz = false;

    // Run the simulation
    simulateBlock(options);
}
//_____________________________________________________________________________
// Spin like top, with damp, with fric
void testBlockSpinLikeTopWithDampWithFric() {
    // Set Options
    SimulationOptions options;
    options.condition = 5;
    options.damping = true;
    options.friction = true;
    options.viz = true;
    options.tilt = 0.0;
    options.tf = 2.0;

    // Check global viz flag
    if(VizOff) options.viz = false;

    // Run the simulation
    simulateBlock(options);
}

//_____________________________________________________________________________
// Tumble, no damp, no fric
void testBlockTumbleNoDampNoFric() {
    // Set Options
    SimulationOptions options;
    options.condition = 6;
    options.damping = false;
    options.friction = false;
    options.viz = true;
    options.tilt = 0.0;
    options.tf = 2.2;

    // Check global viz flag
    if(VizOff) options.viz = false;

    // Run the simulation
    simulateBlock(options);
}
//_____________________________________________________________________________
// Tumble, with damp, with fric
void testBlockTumbleWithDampWithFric() {
    // Set Options
    SimulationOptions options;
    options.condition = 6;
    options.damping = true;
    options.friction = true;
    options.viz = true;
    options.tilt = 0.0;
    options.tf = 10.0;

    // Check global viz flag
    if(VizOff) options.viz = false;

    // Run the simulation
    simulateBlock(options);
}


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
    system.setUseUniformBackground(true);
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0., -9.8, 0.));

    // Define the corners of the block
    Real hs = 0.1;
    Vec3 corner[8];
    corner[0] = Vec3(hs, -hs, hs);
    corner[1] = Vec3(hs, -hs, -hs);
    corner[2] = Vec3(-hs, -hs, -hs);
    corner[3] = Vec3(-hs, -hs, hs);
    corner[4] = Vec3(hs, hs, hs);
    corner[5] = Vec3(hs, hs, -hs);
    corner[6] = Vec3(-hs, hs, -hs);
    corner[7] = Vec3(-hs, hs, hs);

    // Construct the mobilized body.
    Body::Rigid bodyProps(MassProperties(10.0, Vec3(0), Inertia(1)));
    bodyProps.addDecoration(Transform(),
        DecorativeBrick(Vec3(hs)).setColor(Blue));
    MobilizedBody::Free body(matter.Ground(), bodyProps);

    // Parameters
    // Damping
    ExponentialSpringParameters params;
    if(!options.damping) {
        params.setNormalViscosity(0.0);
        params.setFrictionViscosity(0.0);
    }
    // Friction
    Real muk = 0.25, mus = 0.5;
    if(!options.friction) muk = mus = 0.0;
    params.setInitialMuStatic(mus);
    params.setInitialMuKinetic(muk);

    // Floor Plane
    // First point the z-axis up, then rotate about the ground z
    Real angle1 = convertDegreesToRadians(-90.0);
    Real angle2 = convertDegreesToRadians(options.tilt);
    Rotation R;
    R.setRotationFromTwoAnglesTwoAxes(SpaceRotationSequence,
        angle1, XAxis, angle2, ZAxis);
    Transform floorPlane(R, Vec3(0.));
    // Create an exponential spring at each corner of the block.
    unique_ptr<ExponentialSpringForce> sprFloor[8];
    int i;
    for(i = 0; i < 8; ++i) {
        sprFloor[i] = unique_ptr<ExponentialSpringForce>(
            new ExponentialSpringForce(forces, floorPlane, body, corner[i],
                params));
    }
    Rotation RDecor(angle2, ZAxis);
    Transform shiftedFloorPlane(RDecor, Vec3(0.0, -0.1, 0.0));
    matter.Ground().updBody().addDecoration(shiftedFloorPlane,
        DecorativeBrick(Vec3(4, 0.1, 4)).setColor(Gray).setOpacity(1.0));

    // Wall Plane
    // 1.5 m along the +x axis, but angled to face the Ground +z axis.
    // The 30 deg angle is so that motion will be generated that is
    // not perfectly aligned with the coordinate axes.
    Real angle = convertDegreesToRadians(-60.0);
    Rotation r;
    r.setRotationFromAngleAboutY(angle);
    Transform wallPlane(r, Vec3(1.5, 0., 0.));
    // Create an exponential spring to each corner of the block.
    unique_ptr<ExponentialSpringForce> sprWall[8];
    for(i = 0; i < 8; ++i) {
        sprWall[i] = unique_ptr<ExponentialSpringForce>(
            new ExponentialSpringForce(forces, wallPlane, body, corner[i],
                params));
    }
    angle1 = convertDegreesToRadians(90.0);
    angle2 = convertDegreesToRadians(30.0);
    // Add decoration
    Rotation rDecor;
    rDecor.setRotationFromTwoAnglesTwoAxes(SpaceRotationSequence,
        angle1, ZAxis, angle2, YAxis);
    Transform shiftedWallPlane(rDecor, Vec3(1.6, 0., 0.));
    matter.Ground().updBody().addDecoration(shiftedWallPlane,
        DecorativeBrick(Vec3(4, 0.1, 4)).setColor(Vec3(0.9,0.9,0.9)).
        setOpacity(0.8));

    // Periodic Reporter
    PeriodicStateRecorder* periodicRecorder =
        new PeriodicStateRecorder(1.0/60.0);
    system.addEventReporter(periodicRecorder);

    // Min Max Height Reporter
    MinMaxHeightStateRecorder* minmaxRecorder = new
        MinMaxHeightStateRecorder(system, body);
    system.addEventReporter(minmaxRecorder);

    // Add visualization
    // The visualization is intended to verify that the simulations
    // are reasonable.  Visualization should be turned off when unit
    // tests are run, so before committing changes to GitHub.
    unique_ptr<Visualizer> viz;
    unique_ptr<Visualizer::Reporter> vizReporter;
    Visualizer::InputSilo* silo;
    if(options.viz) {
        viz = unique_ptr<Visualizer>(new Visualizer(system));
        viz->setShowShadows(true);
        viz->setShutdownWhenDestructed(true);

        // ----- Run Menu ----
        silo = new Visualizer::InputSilo();
        viz->addInputListener(silo);
        Array_<std::pair<String, int> > runMenuItems;
        runMenuItems.push_back(std::make_pair("Go", GoItem));
        runMenuItems.push_back(std::make_pair("Replay", ReplayItem));
        runMenuItems.push_back(std::make_pair("Quit", QuitItem));
        viz->addMenu("Run", RunMenuId, runMenuItems);

        system.addEventReporter(new Visualizer::Reporter(*viz, 1.0 / 60.0));
    }

    // Realize through Stage::Model (construct the State)
    system.realizeTopology();
    State state = system.getDefaultState();
    system.realizeModel(state);

    // Set the initial conditions for the block
    R.setRotationFromAngleAboutUnitVector(0.0, XAxis);  // No rotation
    Vec3 x(0.5, 0.11, 0.0);     // Hair above the floor
    Vec3 w(0.);                 // No angular velocity
    Vec6 u(0.);                 // All speeds = 0.0
    switch(options.condition) {
    case 1: // Dropped
        x = Vec3(0.5, 1.0, 0.0);
        break;
    case 2: // Sliding
        x = Vec3(-0.5, 0.11, 0.0);
        u = Vec6(0, 0, 0, 4.0, 0, 0);
        break;
    case 3: // Spinning
        w = Vec3(0.0, 8.0, 0.0);
        break;
    case 4: // Sliding and Spinning
        u = Vec6(0, 0, 0, 2.0, 0, 0);
        w = Vec3(0.0, 12.0, 0.0);
        break;
    case 5: // Spinning top
        R.setRotationFromAngleAboutNonUnitVector(
            convertDegreesToRadians(54.74), Vec3(1, 0, 1));
        x = Vec3(0.5, 0.20, 0.0);
        w = Vec3(0.0, 12.0, 0.0);
        break;
    case 6: // Tumbling
        x = Vec3(-1.0, 1.0, 0.5);
        u = Vec6(0, 0, 0, 4.0, 0, 0);
        w = Vec3(-1.0, 0.0, -6.0);
        break;
    case 7: // Corner Shot
        x = Vec3(0.0, 2.5, 1.0);
        u = Vec6(0, 0, 0, 2.0, 0, -2.0);
        w = Vec3(-2.0, 0.0, 0.0);
        break;
    }
    body.setQToFitRotation(state, R);
    body.setQToFitTranslation(state, x);
    body.setU(state, u);
    body.setUToFitAngularVelocity(state, w);

    // Reset the spring zeros
    for(i = 0; i < 8; ++i) {
        sprFloor[i]->resetAnchorPoint(state);
        sprWall[i]->resetAnchorPoint(state);
    }

    // Simulate
    RungeKuttaMersonIntegrator integ(system);
    Real acc = 1.0e-5;
    integ.setAccuracy(acc);
    integ.setMaximumStepSize(0.01);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(options.tf);

    // Get the recorded state arrays
    const Array_<State>* periodicArray = periodicRecorder->getStateArray();
    const Array_<State>* minmaxArray = minmaxRecorder->getStateArray();

    // Test that getting data before Stage::Dynamics has been realized
    // throws an exception.
    const State& pstate = (*periodicArray)[0];
    SimTK_TEST_MUST_THROW(sprFloor[0]->getNormalForce(pstate));

    // Run low level tests for each of the springs
    cout << " ";
    for(i = 0; i < 8; ++i) {

        cout << " f"<< i;
        // periodic floor plane
        checkSpringCalculations(options,
            system, acc, *sprFloor[i], periodicArray);
        if((options.damping == false) && (options.friction == false))
            checkConservationOfEnergy(system, acc, *sprFloor[i], minmaxArray);
        // minmax floor plane
        checkSpringCalculations(options,
            system, acc, *sprFloor[i], periodicArray);
        if((options.damping == false) && (options.friction == false))
            checkConservationOfEnergy(system, acc, *sprFloor[i], minmaxArray);

        cout << " w" << i;
        // periodic wall plane
        checkSpringCalculations(options,
            system, acc, *sprWall[i], periodicArray);
        if((options.damping == false) && (options.friction == false))
            checkConservationOfEnergy(system, acc, *sprFloor[i], minmaxArray);
        // minmax wall plane
        checkSpringCalculations(options,
            system, acc, *sprWall[i], periodicArray);
        if((options.damping == false) && (options.friction == false))
            checkConservationOfEnergy(system, acc, *sprFloor[i], minmaxArray);
    }
    cout << endl;

    // Replay Loop
    if(options.viz) {
        silo->clear(); // forget earlier input
        while(true) {
            cout << "Choose Replay to see that again ...\n";

            int menuId, item;
            silo->waitForMenuPick(menuId, item);

            if(item == QuitItem)
                break;
            if(item != ReplayItem) {
                cout << "\aHuh? Try again.\n";
                continue;
            }

            for(double i = 0; i < (int)periodicArray->size(); i++) {
                viz->report((*periodicArray)[(int)i]);
            }
        }
    }
}
//_____________________________________________________________________________
// Check that spring forces are internally consistent for a simulation.
// system       : system being simulated
// acc          : integrator accuracy
// spr          : exponential spring
// stateArray   : time history of states
void checkSpringCalculations(const SimulationOptions& options,
    MultibodySystem& system, Real acc, ExponentialSpringForce& spr,
    const Array_<State>* stateArray)
{
    // Check for empty state array
    if(stateArray->size() == 0) return;

    // Loop through the recorded states
    for(unsigned int i = 0; i < stateArray->size(); ++i) {

        const State& state = (*stateArray)[i];

        // Realize through Stage::Dynamics
        system.realize(state, Stage::Dynamics);

        // Set a tolerance for comparisons
        Real tol = 1.0e-12;

        // Check the normal force when expressed in the contact plane
        Vec3 fzElas = spr.getNormalForceElasticPart(state, false);
        Vec3 fzDamp = spr.getNormalForceDampingPart(state, false);
        Vec3 fz = spr.getNormalForce(state, false);
        // x components should be 0.0
        SimTK_TEST(fzElas[0] == 0.0);
        SimTK_TEST(fzDamp[0] == 0.0);
        SimTK_TEST(fz[0] == 0.0);
        // z components should be 0.0
        SimTK_TEST(fzElas[1] == 0.0);
        SimTK_TEST(fzDamp[1] == 0.0);
        SimTK_TEST(fz[1] == 0.0);
        // sum elastic and damping
        Vec3 fzCalc = fzElas + fzDamp;
        SimTK_TEST_EQ(fzCalc, fz);

        // Check normal force when expressed in Ground
        Vec3 fzElas_G = spr.getNormalForceElasticPart(state);
        Vec3 fzDamp_G = spr.getNormalForceDampingPart(state);
        Vec3 fz_G = spr.getNormalForce(state);
        Vec3 fzCalc_G = fzElas_G + fzDamp_G;
        SimTK_TEST_EQ(fzCalc_G, fz_G);

        // Check magnitude of fz is same in Ground and ContacPlane
        SimTK_TEST_EQ(fz_G.norm(), fz.norm());

        // Check 0.0 ≤ sliding state ≤ 1.0
        Real sliding = spr.getSliding(state);
        SimTK_TEST(sliding <= 1.0);
        SimTK_TEST(sliding >= 0.0);

        // Check μₖ ≤ μ ≤ μₛ
        Real mus = spr.getMuStatic(state);
        Real muk = spr.getMuKinetic(state);
        Real mu = spr.getMu(state);
        SimTK_TEST(muk <= (mu+acc));
        SimTK_TEST(mu <= (mus+acc));
        if(options.friction) SimTK_TEST(mu > 0.0);

        // Check friction limit ≤ μ*fy
        Vec3 mufz = mu * fz_G;
        Real fricLimit = spr.getFrictionForceLimit(state);
        SimTK_TEST_EQ_TOL(fricLimit, mufz.norm(), tol);

        // Check friction force expressed in the contact plane
        Vec3 fricElas = spr.getFrictionForceElasticPart(state, false);
        Vec3 fricDamp = spr.getFrictionForceDampingPart(state, false);
        Vec3 fric = spr.getFrictionForce(state, false);
        // z component should be 0.0
        SimTK_TEST(fricElas[2] == 0.0);
        SimTK_TEST(fricDamp[2] == 0.0);
        SimTK_TEST(fric[2] == 0.0);
        // sum elastic and damping
        Vec3 fricCalc = fricElas + fricDamp;
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
        SimTK_TEST(fric.norm() <= (fricLimit+tol));

        // Check total force when expressed in Contact Plane
        Vec3 f = spr.getForce(state, false);
        Vec3 fCalc = fzCalc + fricCalc;
        SimTK_TEST_EQ(fCalc, f);

        // Check total force when expressed in Ground
        Vec3 f_G = spr.getForce(state);
        Vec3 fCalc_G = fzCalc_G + fricCalc_G;
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
        Vec3 p0 = spr.getAnchorPointPosition(state, false);
        SimTK_TEST(p0[2] == 0.0);
    }
}

//_____________________________________________________________________________
// Check that energy is reasonably conserved when there is no damping
// and no friction, and the only force acting on the system is Gravity.
// Accounting for dissipative energy losses due to damping and friction
// would require the damping and friction forces to be integrated over time.
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
        Real tol = state.getTime() * acc * energy0;
        Real energy = system.calcEnergy(state);
        SimTK_TEST_EQ_TOL(energy, energy0, tol);
    }
}


//_____________________________________________________________________________
// Test things that happen before integration.
// 1. Setting and getting parameters (ExponentialSpringParameters)
// 2. Construction with default and non-default parameters
// 3. Realization, getting and setting μₛ and μₖ, resetting spring zeros
// This routine is at the very end of the file bcause it is rarely looked at.
void testInitialization() {

    //----------------------------------
    // 1. Setting and getting parameters
    //----------------------------------
    // Test Equality Operator
    ExponentialSpringParameters params, paramsDef;

    // Get the default parameters
    Real d0Def, d1Def, d2Def;
    paramsDef.getShapeParameters(d0Def, d1Def, d2Def);
    Real czDef = paramsDef.getNormalViscosity();
    Real maxFzDef = paramsDef.getMaxNormalForce();
    Real kxyDef = paramsDef.getFrictionElasticity();
    Real cxyDef = paramsDef.getFrictionViscosity();
    Real vSettleDef = paramsDef.getSettleVelocity();
    Real initMusDef = paramsDef.getInitialMuStatic();
    Real initMukDef = paramsDef.getInitialMuKinetic();

    // Make up non-default parameters
    Real delta = 0.1;
    Real d0 = d0Def + delta;
    Real d1 = d1Def + delta;
    Real d2 = d2Def + delta;
    Real cz = czDef + delta;
    Real maxFz = maxFzDef + delta;
    Real kxy = kxyDef + delta;
    Real cxy = cxyDef + delta;
    Real vSettle = vSettleDef + delta;
    Real initMus = initMusDef + delta;
    Real initMuk = initMukDef + delta;

    // Test the Equality Operator
    SimTK_TEST(params==paramsDef);
    params.setShapeParameters(d0, d1, d2);
    SimTK_TEST(!(params == paramsDef));
    params.setNormalViscosity(cz);
    SimTK_TEST(!(params == paramsDef));
    params.setMaxNormalForce(maxFz);
    SimTK_TEST(!(params == paramsDef));
    params.setFrictionElasticity(kxy);
    SimTK_TEST(!(params == paramsDef));
    params.setFrictionViscosity(cxy);
    SimTK_TEST(!(params == paramsDef));
    SimTK_TEST(!(params == paramsDef));
    params.setSettleVelocity(vSettle);
    SimTK_TEST(!(params == paramsDef));
    SimTK_TEST(!(params == paramsDef));
    params.setInitialMuStatic(initMus);
    SimTK_TEST(!(params == paramsDef));
    params.setInitialMuKinetic(initMuk);
    SimTK_TEST(!(params == paramsDef));
    // Now return to the default values, one member variable at a time.
    // This may seem redundant, but different combinations of the member
    // variables have non-default values.
    params.setShapeParameters(d0Def, d1Def, d2Def);
    SimTK_TEST(!(params == paramsDef));
    params.setNormalViscosity(czDef);
    SimTK_TEST(!(params == paramsDef));
    params.setMaxNormalForce(maxFzDef);
    SimTK_TEST(!(params == paramsDef));
    params.setFrictionElasticity(kxyDef);
    SimTK_TEST(!(params == paramsDef));
    params.setFrictionViscosity(cxyDef);
    SimTK_TEST(!(params == paramsDef));
    SimTK_TEST(!(params == paramsDef));
    params.setSettleVelocity(vSettleDef);
    SimTK_TEST(!(params == paramsDef));
    SimTK_TEST(!(params == paramsDef));
    params.setInitialMuStatic(initMusDef);
    SimTK_TEST(!(params == paramsDef));
    params.setInitialMuKinetic(initMukDef);
    SimTK_TEST(params == paramsDef);

    // Test the Inequality Operator
    params.setShapeParameters(d0, d1, d2);
    SimTK_TEST(params != paramsDef);
    params.setNormalViscosity(cz);
    SimTK_TEST(params != paramsDef);
    params.setMaxNormalForce(maxFz);
    SimTK_TEST(params != paramsDef);
    params.setFrictionElasticity(kxy);
    SimTK_TEST(params != paramsDef);
    params.setFrictionViscosity(cxy);
    SimTK_TEST(params != paramsDef);
    SimTK_TEST(params != paramsDef);
    params.setSettleVelocity(vSettle);
    SimTK_TEST(params != paramsDef);
    SimTK_TEST(params != paramsDef);
    params.setInitialMuStatic(initMus);
    SimTK_TEST(params != paramsDef);
    params.setInitialMuKinetic(initMuk);
    SimTK_TEST(params != paramsDef);

    // Test the set methods by checking the individual member variables.
    Real tol = 1.0e-10 * delta;
    params.getShapeParameters(d0, d1, d2);
    SimTK_TEST_EQ(d0 - d0Def, delta);
    SimTK_TEST_EQ(d1 - d1Def, delta);
    SimTK_TEST_EQ_TOL(d2 - d2Def, delta, tol);
    SimTK_TEST_EQ(params.getNormalViscosity() - czDef, delta);
    SimTK_TEST_EQ_TOL(params.getMaxNormalForce() - maxFzDef, delta, tol);
    SimTK_TEST_EQ_TOL(params.getFrictionElasticity() - kxyDef, delta, tol);
    SimTK_TEST_EQ_TOL(params.getFrictionViscosity() - cxyDef, delta, tol);
    SimTK_TEST_EQ(params.getSettleVelocity() - vSettleDef, delta);
    SimTK_TEST_EQ(params.getInitialMuStatic() - initMusDef, delta);
    SimTK_TEST_EQ(params.getInitialMuKinetic() - initMukDef, delta);

    // Test setting viscosity for critical damping,
    // with a default mass (1.0 kg)
    Real kp1 = 100000.0;
    Real kv1 = 2.0 * sqrt(kp1);
    params.setElasticityAndViscosityForCriticalDamping(kp1);
    SimTK_TEST_EQ(params.getFrictionViscosity(), kv1);
    // with a non-default mass
    Real mass = 20.0;
    kv1 = 2.0 * sqrt(kp1 * mass);
    params.setElasticityAndViscosityForCriticalDamping(kp1, mass);
    SimTK_TEST_EQ(params.getFrictionViscosity(), kv1);

    // Return the non-default params to the original non-default values
    params.setFrictionElasticity(kxy);
    params.setFrictionViscosity(cxy);

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
    ExponentialSpringForce
        sprDef0(forces, plane, body, station);
    ExponentialSpringForce
        sprDef1(forces, plane, body, station, paramsDef);
    ExponentialSpringForce
        spr(forces, plane, body, station, params);

    // Test that the non-default parameter values were set properly.
    SimTK_TEST(sprDef0.getParameters() == sprDef1.getParameters());
    SimTK_TEST(spr.getParameters() != sprDef0.getParameters());
    const ExponentialSpringParameters& p = spr.getParameters();
    p.getShapeParameters(d0, d1, d2);
    SimTK_TEST_EQ(d0 - d0Def, delta);
    SimTK_TEST_EQ(d1 - d1Def, delta);
    SimTK_TEST_EQ_TOL(d2 - d2Def, delta, tol);
    SimTK_TEST_EQ(p.getNormalViscosity() - czDef, delta);
    SimTK_TEST_EQ_TOL(p.getMaxNormalForce() - maxFzDef, delta, tol);
    SimTK_TEST_EQ_TOL(p.getFrictionElasticity() - kxyDef, delta, tol);
    SimTK_TEST_EQ_TOL(p.getFrictionViscosity() - cxyDef, delta, tol);
    SimTK_TEST_EQ(p.getSettleVelocity() - vSettleDef, delta);
    SimTK_TEST_EQ(params.getInitialMuStatic() - initMusDef, delta);
    SimTK_TEST_EQ(params.getInitialMuKinetic() - initMukDef, delta);

    // Test getting contact plane, body, and station
    SimTK_TEST(spr.getContactPlaneTransform() == plane);
    SimTK_TEST(spr.getBody().isSameMobilizedBody(body));
    SimTK_TEST(spr.getStation() == station);

    //------------------------------------------------------------------------
    // 3. Realization, getting and setting μₛ and μₖ, resetting spring zeros.
    //------------------------------------------------------------------------
    // Realize through Stage::Model (construct the state)
    system.realizeTopology();
    State state = system.getDefaultState();
    system.realizeModel(state);

    // Test getting μₛ and μₖ via the ExponentialSpringForce API
    Real mus = params.getInitialMuStatic();
    SimTK_TEST(spr.getMuStatic(state) == mus);
    Real muk = params.getInitialMuKinetic();
    SimTK_TEST(spr.getMuKinetic(state) == muk);

    // Test getting μₛ and μₖ via the Subsystem API using indices obtained via
    // getMuStaticStateIndex() and getMuKineticStateIndex().
    const Subsystem& subsys = spr.getForceSubsystem();
    Real musByIndex = SimTK::Value<double>::downcast(
        subsys.getDiscreteVariable(state, spr.getMuStaticStateIndex()) );
    SimTK_TEST(musByIndex == mus);
    Real mukByIndex = SimTK::Value<double>::downcast(
        subsys.getDiscreteVariable(state, spr.getMuKineticStateIndex()) );
    SimTK_TEST(mukByIndex == muk);

    // Test setting μₛ and μₖ via the ExponentialSpringForce API
    Real musNew = mus + 1.0;
    spr.setMuStatic(state, musNew);
    SimTK_TEST(spr.getMuStatic(state) == musNew);
    Real mukNew = muk + 1.0;
    spr.setMuKinetic(state, mukNew);
    SimTK_TEST(spr.getMuKinetic(state) == mukNew);

    // Test that the constraint μₖ ≤ μₛ is enforced
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

    // Set the initial coordinates of the body
    Rotation R;
    R.setRotationFromAngleAboutUnitVector(0.0, XAxis);  // No rotation
    Vec3 x(0.0, 0.2, 0.0);  // 20 cm above the floor
    Vec3 w(0.);             // No angular velocity
    Vec6 u(0.);             // All speeds = 0.0
    body.setQToFitRotation(state, R);
    body.setQToFitTranslation(state, x);
    body.setU(state, u);
    body.setUToFitAngularVelocity(state, w);

    // REALIZATION-RELATED EXCEPTIONS
    // At Stage::Model - All data cache accessors must throw.
    system.realizeModel(state);
    SimTK_TEST_MUST_THROW(spr.getNormalForceElasticPart(state));
    SimTK_TEST_MUST_THROW(spr.getNormalForceDampingPart(state));
    SimTK_TEST_MUST_THROW(spr.getNormalForce(state));
    SimTK_TEST_MUST_THROW(spr.getMu(state));
    SimTK_TEST_MUST_THROW(spr.getFrictionForceLimit(state));
    SimTK_TEST_MUST_THROW(spr.getFrictionForceElasticPart(state));
    SimTK_TEST_MUST_THROW(spr.getFrictionForceDampingPart(state));
    SimTK_TEST_MUST_THROW(spr.getFrictionForce(state));
    SimTK_TEST_MUST_THROW(spr.getForce(state));
    SimTK_TEST_MUST_THROW(spr.getAnchorPointPosition(state));
    SimTK_TEST_MUST_THROW(spr.getStationVelocity(state));
    SimTK_TEST_MUST_THROW(spr.getStationPosition(state));
    // At Stage::Position - all but getStationPosition() must throw
    system.realize(state, Stage::Position);
    SimTK_TEST_MUST_THROW(spr.getNormalForceElasticPart(state));
    SimTK_TEST_MUST_THROW(spr.getNormalForceDampingPart(state));
    SimTK_TEST_MUST_THROW(spr.getNormalForce(state));
    SimTK_TEST_MUST_THROW(spr.getMu(state));
    SimTK_TEST_MUST_THROW(spr.getFrictionForceLimit(state));
    SimTK_TEST_MUST_THROW(spr.getFrictionForceElasticPart(state));
    SimTK_TEST_MUST_THROW(spr.getFrictionForceDampingPart(state));
    SimTK_TEST_MUST_THROW(spr.getFrictionForce(state));
    SimTK_TEST_MUST_THROW(spr.getForce(state));
    SimTK_TEST_MUST_THROW(spr.getAnchorPointPosition(state));
    SimTK_TEST_MUST_THROW(spr.getStationVelocity(state));
    spr.getStationPosition(state);
    // At Stage::Velocity - all but getStationPosition() Velocity() must throw
    system.realize(state, Stage::Velocity);
    SimTK_TEST_MUST_THROW(spr.getNormalForceElasticPart(state));
    SimTK_TEST_MUST_THROW(spr.getNormalForceDampingPart(state));
    SimTK_TEST_MUST_THROW(spr.getNormalForce(state));
    SimTK_TEST_MUST_THROW(spr.getMu(state));
    SimTK_TEST_MUST_THROW(spr.getFrictionForceLimit(state));
    SimTK_TEST_MUST_THROW(spr.getFrictionForceElasticPart(state));
    SimTK_TEST_MUST_THROW(spr.getFrictionForceDampingPart(state));
    SimTK_TEST_MUST_THROW(spr.getFrictionForce(state));
    SimTK_TEST_MUST_THROW(spr.getForce(state));
    SimTK_TEST_MUST_THROW(spr.getAnchorPointPosition(state));
    spr.getStationVelocity(state);
    spr.getStationPosition(state);
    // At Stage::Dynamics - None should throw an exception.
    system.realize(state, Stage::Dynamics);
    spr.getNormalForceElasticPart(state);
    spr.getNormalForceDampingPart(state);
    spr.getNormalForce(state);
    spr.getMu(state);
    spr.getFrictionForceLimit(state);
    spr.getFrictionForceElasticPart(state);
    spr.getFrictionForceDampingPart(state);
    spr.getFrictionForce(state);
    spr.getForce(state);
    spr.getAnchorPointPosition(state);
    spr.getStationVelocity(state);
    spr.getStationPosition(state);
    // At Stage::Acceleration - None should throw an exception.
    system.realize(state, Stage::Acceleration);
    spr.getNormalForceElasticPart(state);
    spr.getNormalForceDampingPart(state);
    spr.getNormalForce(state);
    spr.getMu(state);
    spr.getFrictionForceLimit(state);
    spr.getFrictionForceElasticPart(state);
    spr.getFrictionForceDampingPart(state);
    spr.getFrictionForce(state);
    spr.getForce(state);
    spr.getAnchorPointPosition(state);
    spr.getStationVelocity(state);
    spr.getStationPosition(state);

    // Test resetting the spring zero
    // After the reset the spring zero should coincide with the projection of
    // the spring station onto the contact plane, which should result in the
    // elastic component of the frictional force being zero.
    // Expected positions after the reset
    // express spring station in the Ground frame
    Vec3 s_G = body.findStationLocationInGround(state, station);
    // express station in the contact plane frame
    Vec3 s = plane.shiftBaseStationToFrame(s_G);
    // project onto the contact plane by setting the z-component to 0.0
    Vec3 p0After = s;  p0After[2] = 0.0;
    // express in Ground frame after the projection
    Vec3 p0After_G = plane.shiftFrameStationToBase(p0After);
    // Now reset
    spr.resetAnchorPoint(state);
    system.realize(state, Stage::Dynamics);
    SimTK_TEST(p0After == spr.getAnchorPointPosition(state, false));
    SimTK_TEST(p0After_G == spr.getAnchorPointPosition(state));
    // Test that the elastic component of the friction force is 0.0
    Vec3 fElastic = spr.getFrictionForceElasticPart(state);
    SimTK_TEST_EQ(fElastic, Vec3(0., 0., 0.));

    // Test getting the Sliding state via the Subsystem API using indices
    // obtained via getSlidingStateIndex().
    Real sliding = spr.getSliding(state);
    Real slidingByIndex = SimTK::Value<double>::downcast(
        subsys.getDiscreteVariable(state, spr.getSlidingStateIndex()) );
    SimTK_TEST(slidingByIndex == sliding);

    // Test getting the Elastic Anchor Point state via the Subsystem API using
    // indices obtained via getAnchorPOintStateIndex().
    Vec3 anchor = spr.getAnchorPointPosition(state, false);
    Vec3 anchorByIndex = SimTK::Value<Vec3>::downcast(
        subsys.getDiscreteVariable(state, spr.getAnchorPointStateIndex()) );
    SimTK_TEST(anchorByIndex == anchor);
};

