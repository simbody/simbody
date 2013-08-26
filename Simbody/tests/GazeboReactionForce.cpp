/* -------------------------------------------------------------------------- *
 *                     Simbody(tm): Gazebo Reaction Force                     *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
 * Authors: Michael Sherman, John Hsu                                                   *
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

/* This test is drawn from the Open Source Robotics Foundation Gazebo physics
regression test "Joint_TEST::ForceTorque". 

It is a small sphere pinned at its center to ground, rotating about X, then
an outboard brick pinned to the sphere and rotating about mutual Z. Oblique
gravity cause the brick to fall, hit the ground and come to rest. Then we check
against John Hsu's hand coded reaction force values.
*/

#include "Simbody.h"

#include <cassert>
#include <iostream>
using std::cout; using std::endl;

using namespace SimTK;

//#define USE_VISUALIZER

const Real Radius=.1; // for sphere attached to link1
const Vec3 Box(.05,.1,.2); // half-dimensions of box
const Real Mass1=10, Mass2=10;
const Vec3 COM1(0,0,.5), COM2(0,0,.5);
const Inertia Central(Vec3(.1,.1,.1)); // supposedly mass weighted already
// Simbody requires inertias to be expressed about body origin rather than COM.
const Inertia Inertia1=Central.shiftFromMassCenter(-COM1, Mass1);
const Inertia Inertia2=Central.shiftFromMassCenter(-COM2, Mass2);

// Define a stiff, lossy material.
const Real Stiffness = 1e8;
const Real Dissipation = 10;
const Real Mu_s = 0.15, Mu_d = 0.1, Mu_v = 0;
const ContactMaterial lossyMaterial(Stiffness,
                                    Dissipation,
                                    Mu_s,
                                    Mu_d,
                                    Mu_v);

const Real MaxStepSize    = Real(1/1000.); // 1 ms (1000 Hz)
const int  DrawEveryN     = 33;            // 33 ms frame update (30.3 Hz)
const Real SimTime        = 2;
const int  NSteps         = // make this a whole number of viz frames
    DrawEveryN*(int(SimTime/MaxStepSize/DrawEveryN+0.5));

// Use this class to hold references into the Simbody system.
struct MyMultibodySystem {
    MyMultibodySystem(); // see below

    MultibodySystem                 m_system;
    SimbodyMatterSubsystem          m_matter;
    ContactTrackerSubsystem         m_tracker;
    CompliantContactSubsystem       m_contact;
    GeneralForceSubsystem           m_forces;
    Force::Gravity                  m_gravity;

    #ifdef USE_VISUALIZER
    Visualizer                      m_viz;
    #endif

    MobilizedBody::Pin              m_link1, m_link2;
};

// Execute the simulation with a given integrator and accuracy, and verify that
// it produces the correct answers.
static void runOnce(const MyMultibodySystem& mbs, Integrator& integ, 
                    Real accuracy);


//==============================================================================
//                                   MAIN
//==============================================================================
int main() {
    SimTK_START_TEST("GazeboReactionForce");
        // Create the system.   
        MyMultibodySystem mbs;

        printf("LOW ACCURACY\n");
        SemiExplicitEuler2Integrator sexpeul2(mbs.m_system);
        SimTK_SUBTEST3(runOnce, mbs, sexpeul2, 1e-2);

        printf("\nHIGH ACCURACY\n");
        RungeKuttaMersonIntegrator rkm(mbs.m_system);
        SimTK_SUBTEST3(runOnce, mbs, rkm, 1e-6);

    SimTK_END_TEST();
}


//==============================================================================
//                           MY MULTIBODY SYSTEM
//==============================================================================
// Construct the multibody system.
MyMultibodySystem::MyMultibodySystem()
:   m_system(), m_matter(m_system), 
    m_tracker(m_system), m_contact(m_system,m_tracker),
    m_forces(m_system), m_gravity(m_forces,m_matter,Vec3(-30,10,-50))
    #ifdef USE_VISUALIZER
    , m_viz(m_system)
    #endif
{
    #ifdef USE_VISUALIZER
    m_viz.setSystemUpDirection(ZAxis);
    m_viz.setCameraTransform(
        Transform(Rotation(BodyRotationSequence,Pi/2,XAxis,Pi/8,YAxis),
        Vec3(5,-8,2)));
    m_viz.setShowFrameNumber(true);
    m_viz.setShowSimTime(true);
    #endif

    Body::Rigid link1Info(MassProperties(Mass1, COM1, Inertia1));
    link1Info.addDecoration(COM1, DecorativeSphere(Radius).setOpacity(.3));

    Body::Rigid link2Info(MassProperties(Mass2, COM2, Inertia2));
    link2Info.addDecoration(COM2, DecorativeBrick(Box).setOpacity(.3));

    ContactGeometry::TriangleMesh cubeMesh
        (PolygonalMesh::createBrickMesh(Box, 3));
    DecorativeMesh showMesh(cubeMesh.createPolygonalMesh());
    showMesh.setRepresentation(DecorativeGeometry::DrawWireframe);
    link2Info.addContactSurface(COM2,
        ContactSurface(cubeMesh, lossyMaterial, 1));
    link2Info.addDecoration(COM2, showMesh);

    MobilizedBody& Ground = m_matter.updGround(); // Nicer name for Ground.
    Ground.addBodyDecoration(Vec3(0,0,.01),
        DecorativeFrame(2).setColor(Green));

    // Add the Ground contact geometry. Contact half space has -XAxis normal
    // (right hand wall) so we have to rotate.
    const Rotation NegXToZ(Pi/2, YAxis);
    Ground.updBody().addContactSurface(Transform(NegXToZ,Vec3(0)),
        ContactSurface(ContactGeometry::HalfSpace(),lossyMaterial));

    const Rotation ZtoX(Pi/2, YAxis);
    m_link1 = MobilizedBody::Pin(Ground,Transform(ZtoX,COM1), 
                                 link1Info, Transform(ZtoX,COM1));

    m_link2 = MobilizedBody::Pin(m_link1, Vec3(0,0,1.5), 
                                 link2Info, Vec3(0));

    Force::MobilityLinearStop
       (m_forces, m_link1, MobilizerQIndex(0),1000000.,1.,-Pi/2,Pi/2);
    Force::MobilityLinearStop
       (m_forces, m_link2, MobilizerQIndex(0), 1000000.,1.,-Infinity,Infinity);

    m_system.realizeTopology();

}

//==============================================================================
//                                 RUN ONCE
//==============================================================================

// Reaction force information: results are the joint reaction, at the F frame 
// on the parent and M frame on the child, expressed in the parent or child 
// frame, resp. Note that Gazebo's GetForceTorque() method uses the negation of
// the joint reaction, and Gazebo's results are ordered (force,torque) rather 
// than (torque,force) as in a Simbody SpatialVec.
struct ReactionPair {
    SpatialVec reactionOnParentInParent;
    SpatialVec reactionOnChildInChild;
};
static ReactionPair getReactionPair(const State&         state,
                                    const MobilizedBody& mobod); 

// Write interesting integrator info to stdout.
static void dumpIntegratorStats(const Integrator& integ);

// Run the system until it settles down, then check the answers.
void runOnce(const MyMultibodySystem& mbs, Integrator& integ, Real accuracy) 
{
    integ.setAllowInterpolation(false);
    integ.setAccuracy(accuracy);
    integ.initialize(mbs.m_system.getDefaultState());

    printf("Test with order %d integator %s, Accuracy=%g, MaxStepSize=%g\n",
        integ.getMethodMinOrder(), integ.getMethodName(), 
        integ.getAccuracyInUse(), MaxStepSize); 

    #ifdef USE_VISUALIZER
    mbs.m_viz.report(integ.getState());
    printf("Hit ENTER to simulate:");
    getchar();
    #endif

    unsigned stepNum = 0;
    while (true) {
        // Get access to State being advanced by the integrator. Interpolation 
        // must be off so that we're modifying the actual trajectory.
        State& state = integ.updAdvancedState();

        #ifdef USE_VISUALIZER
        // Output a frame to the Visualizer if it is time.
        if (stepNum % DrawEveryN == 0)
            mbs.m_viz.report(state);
        #endif

        if (stepNum++ == NSteps)
            break;

        // Advance time by MaxStepSize. Might take multiple internal steps to 
        // get there, depending on difficulty and required accuracy.
        const Real tNext = stepNum * MaxStepSize;
        do {integ.stepTo(tNext,tNext);} while (integ.getTime() < tNext);
    }

    const State& state = integ.getAdvancedState();
    mbs.m_system.realize(state);

    ReactionPair reaction1 = getReactionPair(state, mbs.m_link1);
    ReactionPair reaction2 = getReactionPair(state, mbs.m_link2);
    cout << "t=" << state.getTime() << endl;
    cout << "  Reaction 1p=" << reaction1.reactionOnParentInParent << "\n"; 
    cout << "  Reaction 1c=" << reaction1.reactionOnChildInChild << "\n"; 
    cout << "  Reaction 2p=" << reaction2.reactionOnParentInParent << "\n"; 
    cout << "  Reaction 2c=" << reaction2.reactionOnChildInChild << "\n"; 

    // Check the answers. Note (torque,force) ordering.
    SimTK_TEST_EQ_TOL(reaction1.reactionOnParentInParent,
                      SpatialVec(Vec3(-750, 0, 450), Vec3(-600,200,-1000)), 0.5);
    SimTK_TEST_EQ_TOL(reaction1.reactionOnChildInChild,
                      SpatialVec(Vec3(750,450,0), Vec3(600,-1000,-200)), 0.5);
    SimTK_TEST_EQ_TOL(reaction2.reactionOnParentInParent,
                      SpatialVec(Vec3(-250,-150,0), Vec3(-300,500,100)), 0.5);
    SimTK_TEST_EQ_TOL(reaction2.reactionOnChildInChild,
                      SpatialVec(Vec3(250,150,0), Vec3(300,-500,-100)), 0.5);

    dumpIntegratorStats(integ);
}



//==============================================================================
//                            GET REACTION PAIR
//==============================================================================
static ReactionPair getReactionPair(const State&         state,
                                    const MobilizedBody& mobod) 
{
    SpatialVec p = mobod.findMobilizerReactionOnParentAtFInGround(state);
    SpatialVec c = mobod.findMobilizerReactionOnBodyAtMInGround(state);
    const Rotation& R_GC = mobod.getBodyRotation(state);
    const Rotation& R_GP = mobod.getParentMobilizedBody().getBodyRotation(state);
    ReactionPair pair;
    pair.reactionOnChildInChild   = ~R_GC*c; // from Ground to Child
    pair.reactionOnParentInParent = ~R_GP*p; // from Ground to Parent
    return pair;
}


//==============================================================================
//                        DUMP INTEGRATOR STATS
//==============================================================================
static void dumpIntegratorStats(const Integrator& integ) {
    const int evals = integ.getNumRealizations();
    std::cout << "\nDone -- simulated " << integ.getTime() << "s with " 
            << integ.getNumStepsTaken() << " steps, avg step=" 
        << (1000*integ.getTime())/integ.getNumStepsTaken() << "ms " 
        << (1000*integ.getTime())/evals << "ms/eval\n";

    printf("Used Integrator %s at accuracy %g:\n", 
        integ.getMethodName(), integ.getAccuracyInUse());
    printf("# STEPS/ATTEMPTS = %d/%d\n",  integ.getNumStepsTaken(), 
                                          integ.getNumStepsAttempted());
    printf("# ERR TEST FAILS = %d\n",     integ.getNumErrorTestFailures());
    printf("# REALIZE/PROJECT = %d/%d\n", integ.getNumRealizations(), 
                                          integ.getNumProjections());
}
