/* -------------------------------------------------------------------------- *
 *          Simbody(tm): Gazebo Reaction Force With Applied Force             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
 * Authors: Michael Sherman, John Hsu                                         *
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
regression test "Joint_TEST::GetForceTorqueWithAppliedForce". 

It is a stack of three cubes hinged together at their edges. The bottom block
is heavy and rests on the ground, the other two are light and have their 
positions maintained by a pair of PD controllers like this:

              / 
     link3  /   \         * = pin joint
           /     \
           \  1  /
             \  /                   z
        ------*  45 degrees         ^                     g = 0 0 -50
  link2 |     |                     |   y                       |
        |  6  |                     |  /                        |
   -----*------                     | /                         v 
  |     | 0 degrees                 ---------> x
  | 100 |
  ------- link1
  contact
   GROUND

All the cube edges are of length 1. Masses are 100,6,1 as shown. The top block 
has a COM that is offset into the +y direction by 0.5. 

Expected reaction force results:
    pin1 on inboard:  tx=-25, ty= 175, fz=-300
    pin1 on outboard: tx= 25, ty=-175, fx=1, fz= 300, 
    pin2 on inboard:  tx=-25, fz=-50
    pin2 in outboard: tx=25*pi/4, tz=-25*pi/4, fx=fz=50*pi/4 

The reason for fx=1 in the second line is that with a gain of only 50000, the
first pin joint must be at an angle of 0.0035 radians to generate a torque of 
175. That tips link2's frame in which (0,0,300)_G is reexpressed.
*/

#include "Simbody.h"

#include <cassert>
#include <iostream>
using std::cout; using std::endl;

using namespace SimTK;

//#define USE_VISUALIZER


// Control gains
const Real Kp1 = 50000; // link1-2 joint stiffness
const Real Kp2 = 10000; // link2-3 joint stiffness
const Real Cd1 = 30;    // link1-2 joint damping
const Real Cd2 = 30;    // link2-3 joint damping

// Target angles
const Real Target1 = 0;
const Real Target2 = -Pi/4;

const Vec3 Cube(.5,.5,.5); // half-dimensions of cube
const Real Mass1=100, Mass2=5, Mass3=1;
const Vec3 Centroid(.5,0,.5);
const Vec3 COM1=Centroid, COM2=Centroid, COM3=Centroid+Vec3(0,.5,0);
const UnitInertia Central(Vec3(.1), Vec3(.05));
// Simbody requires inertias to be expressed about body origin rather than COM.
const Inertia Inertia1=Mass1*Central.shiftFromCentroid(-COM1);
const Inertia Inertia2=Mass2*Central.shiftFromCentroid(-COM2);
const Inertia Inertia3=Mass3*Central.shiftFromCentroid(-COM3); // weird

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
const Real SimTime        = 5;
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
    Force::DiscreteForces           m_discrete;
    #ifdef USE_VISUALIZER
    Visualizer                      m_viz;
    #endif

    MobilizedBody                   m_link1;
    MobilizedBody::Pin              m_link2, m_link3;
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
// Construct the multibody system. The dampers are built in here but the springs
// are applied during execution.
MyMultibodySystem::MyMultibodySystem()
:   m_system(), m_matter(m_system), 
    m_tracker(m_system), m_contact(m_system,m_tracker),
    m_forces(m_system), m_discrete(m_forces,m_matter)
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

    Force::Gravity(m_forces, m_matter, -ZAxis, 50);

    DecorativeBrick drawCube(Cube); drawCube.setOpacity(0.5).setColor(Gray);

    Body::Rigid link1Info(MassProperties(Mass1, COM1, Inertia1));
    ContactGeometry::TriangleMesh cubeMesh
        (PolygonalMesh::createBrickMesh(Cube, 3));
    DecorativeMesh showMesh(cubeMesh.createPolygonalMesh());
    showMesh.setRepresentation(DecorativeGeometry::DrawWireframe);
    link1Info.addDecoration(Centroid, drawCube);
    link1Info.addContactSurface(Centroid,
        ContactSurface(cubeMesh, lossyMaterial, 1));
    link1Info.addDecoration(Centroid, showMesh);

    Body::Rigid link2Info(MassProperties(Mass2, COM2, Inertia2));
    link2Info.addDecoration(Centroid, drawCube);

    Body::Rigid link3Info(MassProperties(Mass3, COM3, Inertia3));
    link3Info.addDecoration(Centroid, drawCube);    

    MobilizedBody& Ground = m_matter.updGround(); // Nicer name for Ground.
    Ground.addBodyDecoration(Vec3(0,0,.05),
        DecorativeFrame(2).setColor(Green));

    // Add the Ground contact geometry. Contact half space has -XAxis normal
    // (right hand wall) so we have to rotate.
    const Rotation NegXToZ(Pi/2, YAxis);
    Ground.updBody().addContactSurface(Transform(NegXToZ,Vec3(0)),
        ContactSurface(ContactGeometry::HalfSpace(),lossyMaterial));
    m_link1 = MobilizedBody::Free(Ground,Vec3(0), 
                                    link1Info, Vec3(0));

    // Use this instead of the free joint to remove contact.
    //m_link1 = MobilizedBody::Weld(Ground,Vec3(0), 
    //                              link1Info, Vec3(0));

    const Rotation ZtoY(-Pi/2, XAxis);
    m_link2 = MobilizedBody::Pin(m_link1, Transform(ZtoY,2*Centroid), 
                                    link2Info, Transform(ZtoY,Vec3(0)));
    m_link3 = MobilizedBody::Pin(m_link2, Transform(ZtoY,2*Centroid), 
                                    link3Info, Transform(ZtoY,Vec3(0)));

    // It is more stable to build the springs into the mechanism like
    // this rather than apply them discretely.
    //Force::MobilityLinearSpring
    //   (m_forces, m_link2, MobilizerQIndex(0), Kp1, Target1);
    //Force::MobilityLinearSpring
    //   (m_forces, m_link3, MobilizerQIndex(0), Kp2, Target2);
    Force::MobilityLinearDamper
        (m_forces, m_link2, MobilizerUIndex(0), Cd1);
    Force::MobilityLinearDamper
        (m_forces, m_link3, MobilizerUIndex(0), Cd2);

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

        // Apply discrete spring forces.
        const Real a1err = mbs.m_link2.getAngle(state)-Target1;
        const Real a2err = mbs.m_link3.getAngle(state)-Target2;
        mbs.m_discrete.setOneMobilityForce(state, mbs.m_link2, 
            MobilizerUIndex(0), -Kp1*a1err);
        mbs.m_discrete.setOneMobilityForce(state, mbs.m_link3, 
            MobilizerUIndex(0), -Kp2*a2err);

        // Advance time by MaxStepSize. Might take multiple internal steps to 
        // get there, depending on difficulty and required accuracy.
        const Real tNext = stepNum * MaxStepSize;
        do {integ.stepTo(tNext,tNext);} while (integ.getTime() < tNext);
    }

    const State& state = integ.getAdvancedState();
    mbs.m_system.realize(state);
    Rotation R_G1 = mbs.m_link1.getBodyRotation(state);
    Vec3 a = R_G1.convertRotationToBodyFixedXYZ();
    Vec3 w = mbs.m_link1.getBodyAngularVelocity(state);
    ReactionPair reaction2 = getReactionPair(state, mbs.m_link2);
    ReactionPair reaction3 = getReactionPair(state, mbs.m_link3);
    cout << "t=" << state.getTime() << endl;
    cout << "  joint1 a=" << a << " w=" << w << endl;
    cout << "  joint2 qerr=" << mbs.m_link2.getAngle(state)-Target1
            << " u=" << mbs.m_link2.getRate(state) << endl;
    cout << "  joint3 qerr=" << mbs.m_link3.getAngle(state)-Target2
            << " u=" << mbs.m_link3.getRate(state) << endl;
    cout << "  Reaction 2p=" << reaction2.reactionOnParentInParent << "\n"; 
    cout << "  Reaction 2c=" << reaction2.reactionOnChildInChild << "\n"; 
    cout << "  Reaction 3p=" << reaction3.reactionOnParentInParent << "\n"; 
    cout << "  Reaction 3c=" << reaction3.reactionOnChildInChild << "\n"; 

    // Check the answers. Note (torque,force) ordering.
    SimTK_TEST_EQ_TOL(reaction2.reactionOnParentInParent,
                      SpatialVec(Vec3(-25, 175,0), Vec3(0,0,-300)), 0.5);
    SimTK_TEST_EQ_TOL(reaction2.reactionOnChildInChild,
                      SpatialVec(Vec3( 25,-175,0), Vec3(-1,0, 300)), 0.5);
    SimTK_TEST_EQ_TOL(reaction3.reactionOnParentInParent,
                      SpatialVec(Vec3(-25,0,0), Vec3(0,0,-50)), 0.5);
    SimTK_TEST_EQ_TOL(reaction3.reactionOnChildInChild,
                      (Pi/4)*SpatialVec(Vec3(25,0,-25), Vec3(50,0,50)), 0.5);

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
