/* -------------------------------------------------------------------------- *
 *                  Simbody(tm) Example: Operational Space                    *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-14 Stanford University and the Authors.        *
 * Authors: Chris Dembia, Jack Wang                                           *
 * Contributors: Michael Sherman                                              *
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

/**
 * This example shows how one could write an operational space controller in
 * Simbody. The model we are controlling is a human upper body model.
 * Operational space controllers are model-based, and we show how to access the
 * necessary quantities from the SimbodyMatterSubsystem. We assume the system
 * has a motor at each of its degrees of freedom.
 *
 * The task the controller will achieve has two components:
 * 1. One arm reaches for a target point
 * 2. All bodies are subject to gravity compensation (to counteract the effect
 * of gravity).
 */

#include <Simbody.h>

using namespace std;
using namespace SimTK;

#define DEG(rad) convertRadiansToDegrees(rad)
#define RAD(deg) convertDegreesToRadians(deg)

//==============================================================================
// UpperBody system
//==============================================================================

/// This system is copied from the Biped system.
/// x is directed anteriorly.
/// y is up.
/// z is to the UpperBody's right.
class UpperBody : public MultibodySystem {
public:
    UpperBody();

    const SimbodyMatterSubsystem& getMatterSubsystem() const {return m_matter;}
    SimbodyMatterSubsystem& updMatterSubsystem() {return m_matter;}

    const GeneralForceSubsystem& getForceSubsystem() const {return m_forces;}
    GeneralForceSubsystem& updForceSubsystem() {return m_forces;}

    const Force::Gravity& getGravity() const {
        return static_cast<const Force::Gravity&>(
                m_forces.getForce(m_gravityIndex));
    }

    /// Realizes the topology and model, and uses the resulting initial state
    /// to perform further internal initialization.
    void initialize(State& state) {
        state = realizeTopology();
        realizeModel(state);

        // Initialization specific to UpperBody.
        fillInCoordinateMap(state);
    }

    /// Populates a map that allow users of UpperBody to access coordinates (Q)
    /// and speeds (U) with descriptive Enums (e.g., neck_extension) so the
    /// user (e.g., SIMBICON) needn't know how State indices correspond to
    /// coordinates.
    void fillInCoordinateMap(const State& s);

    /// The names of the rigid bodies that make up the system.
    /// Used for indexing into data structures.
    enum Segment {
        trunk, head,
        upperarm_r, lowerarm_r, hand_r, upperarm_l, lowerarm_l, hand_l,
    };

    /// The names of the degrees of freedom of the system.
    /// Used for indexing into data structures.
    enum Coordinate {
        neck_extension, neck_bending, neck_rotation,
        shoulder_r_flexion, shoulder_r_adduction, shoulder_r_rotation,
        elbow_r_flexion, elbow_r_rotation,
        shoulder_l_flexion, shoulder_l_adduction, shoulder_l_rotation,
        elbow_l_flexion, elbow_l_rotation
    };

    /// Index, in the State vector, of this coordinate's value (Q).
    QIndex getQIndex(Coordinate coord) const {
        assert(!m_coordinates.empty());
        return m_coordinates.at(coord).first;
    }

    /// Coordinate coord's value (Q).
    Real getQ(const State& s, Coordinate coord) const {
        return s.getQ()[getQIndex(coord)];
    }

    /// Index, in the State vector, of this coordinate's speed (U).
    UIndex getUIndex(Coordinate coord) const {
        assert(!m_coordinates.empty());
        return m_coordinates.at(coord).second;
    }

    /// Coordinate coord's speed (U).
    Real getU(const State& s, Coordinate coord) const {
        return s.getU()[getUIndex(coord)];
    }

    /// Pull out from mobForces the mobility force corresponding to the
    /// provided coord.
    Real getForce(const Coordinate coord, const Vector& mobForces)
    {
        return mobForces[getUIndex(coord)];
    }

    /// Add in a generalized force into the proper place in `mobForces`,
    /// as determined by `coord`.
    void addInForce(const Coordinate coord, const Real force,
            Vector& mobForces) const {
        mobForces[getUIndex(coord)] += force;
    }

    const MobilizedBody& getBody(Segment seg) const {
        return m_bodies.at(seg);
    }

private:

    /// For convenience in building the model.
    /// Lower bound and upper bound should be in degrees.
    void addMobilityLinearStop(const MobilizedBody& mobod,
            unsigned int idx,
            Real lowerBound, Real upperBound,
            Real stopStiffness=DEG(1000), Real stopDissipation=1.0)
    {
        Force::MobilityLinearStop(m_forces, mobod, MobilizerQIndex(idx),
                stopStiffness, stopDissipation,
                RAD(lowerBound), RAD(upperBound));
    }

    /// The indices, in the State vector, of a coordinate value (Q) and speed
    /// (U) in a specific MobilizedBody. Used in fillInCoordinateMap().
    static std::pair<QIndex, UIndex> getQandUIndices(
            const SimTK::State& s, const SimTK::MobilizedBody& mobod, int which)
    {
        SimTK::QIndex q0 = mobod.getFirstQIndex(s);
        SimTK::UIndex u0 = mobod.getFirstUIndex(s);
        return std::make_pair(QIndex(q0 + which), UIndex(u0 + which));
    }

    // Subsystems.
    SimbodyMatterSubsystem       m_matter;
    GeneralForceSubsystem        m_forces;

    // Index of the Gravity Force in m_forces.
    ForceIndex                   m_gravityIndex;

    /// The MobilizedBody's that make up the UpperBody.
    std::map<Segment, MobilizedBody> m_bodies;

    /// The indices of each coordinate's Q and U in the State vector.
    std::map<Coordinate, std::pair<QIndex, UIndex> > m_coordinates;
};

UpperBody::UpperBody()
    : m_matter(*this), m_forces(*this)
{

    m_matter.setShowDefaultGeometry(false);

    //--------------------------------------------------------------------------
    //                          Constants, etc.
    //--------------------------------------------------------------------------

    // Miscellaneous.
    // --------------
    // Original OpenSim ellipsoid_center.vtp half dimensions.
    const Vec3 ectr(.03, .12, .03);
    // Original OpenSim sphere.vtp half dimension.
    const Real rad = .5;
    // Original OpenSim block.vtp half dimensions.
    const Vec3 blk(.05,.05,.05);

    // Some convenient rotation matrices for relabeling axes.
    const Rotation Rzmxmy(BodyRotationSequence,  Pi/2, XAxis,  Pi/2, ZAxis);
    const Rotation Rzxy(BodyRotationSequence,   -Pi/2, XAxis, -Pi/2, ZAxis);
    const Rotation Ryzx(BodyRotationSequence,    Pi/2, YAxis,  Pi/2, ZAxis);

    DecorativeSphere originMarker(0.04*rad);
    originMarker.setColor(Vec3(.1,.1,.1));


    //--------------------------------------------------------------------------
    //                          Gravity
    //--------------------------------------------------------------------------

    Force::Gravity gravity(m_forces, m_matter, -YAxis, 9.8066);
    m_gravityIndex = gravity.getForceIndex();

    //--------------------------------------------------------------------------
    //                          Body information
    //--------------------------------------------------------------------------
    // Trunk.
    // ------
    const Real trunkMass = 19.733716299458440;
    const Vec3 trunkCOM(0,0,0);
    Body trunkInfo(MassProperties(trunkMass, trunkCOM,
                                  Inertia(0.355666710204554,
                                          0.224533650416368,
                                          0.281526683481324,
                                          0.046737850487895,
                                          0.000655032101243,
                                         -0.000572934744554)));

    trunkInfo.addDecoration(Transform(Rotation(.01, ZAxis),
                                      Vec3(-0.025,-0.06,0)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(1.6,1.25,1.6)))
            .setColor(White).setOpacity(1));
    trunkInfo.addDecoration(Transform(Rotation(-.3, ZAxis),
                                      Vec3(-0.01,0.16,0)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(0.85,0.55,0.85)))
            .setColor(White).setOpacity(1));
    trunkInfo.addDecoration(Transform(Ryzx,Vec3(-0.025,0.09,0)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(1,1.43,1)))
            .setColor(White).setOpacity(1));
    trunkInfo.addDecoration(Vec3(0),originMarker);

    // Head.
    // -----
    const Real headMass = 4.340524803328146;
    const Vec3 headCOM(0.041488177946752, 0.085892319249215, 0);
    Body headInfo(MassProperties(headMass, headCOM,
                                 Inertia( 0.020576741740382,
                                          0.016008984554380,
                                          0.025555859085964,
                                         -0.004554656543977,
                                         -0.000103058383929,
                                          0.000186029116753)
                                 .shiftFromMassCenter(-headCOM,headMass)));
    headInfo.addDecoration(Vec3(0,.11,0),
        DecorativeEllipsoid(Vec3(0.2,0.22,0.2)/2.)
            .setColor(White).setOpacity(1));
    headInfo.addDecoration(Vec3(0),originMarker);

    // Upper arm.
    // ----------
    // COM z has opposite sign left to right.
    const Real upperarmMass = 2.070989783095760;
    const Vec3 upperarm_rCOM(0.003289136233947,-0.078058926824158,0.065606556342984);
    Body upperarm_rInfo(MassProperties(upperarmMass, upperarm_rCOM,
                      Inertia(0.014082101994221,
                              0.003604979502976,
                              0.013769380880709)
                      .shiftFromMassCenter(-upperarm_rCOM, upperarmMass)));
    upperarm_rInfo.addDecoration(Transform(
                    Rotation(BodyRotationSequence,-0.47,XAxis,0.13,ZAxis),
                    Vec3(0.017,-0.12,0.06)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(1,1.2,1)))
            .setColor(White).setOpacity(1));
    upperarm_rInfo.addDecoration(Vec3(0),originMarker);

    const Vec3 upperarm_lCOM(0.003289136233947,-0.078058926824158,-0.065606556342984);
    Body upperarm_lInfo(MassProperties(upperarmMass, upperarm_lCOM,
                      Inertia(0.014082101994221,
                              0.003604979502976,
                              0.013769380880709)
                      .shiftFromMassCenter(-upperarm_lCOM, upperarmMass)));
    upperarm_lInfo.addDecoration(Transform(
                    Rotation(BodyRotationSequence,0.47,XAxis,0.13,ZAxis),
                    Vec3(0.017,-0.12,-0.06)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(1,1.2,1)))
            .setColor(White).setOpacity(1));
    upperarm_lInfo.addDecoration(Vec3(0),originMarker);

    // Lower arm.
    // ----------
    // Some signs differ left to right.
    const Real lowerarmMass = 1.106702647320712;
    const Vec3 lowerarm_rCOM(0.031656703591848,-0.089369993258598,0.017231110378866);
    Body lowerarm_rInfo(MassProperties(lowerarmMass, lowerarm_rCOM,
                            Inertia( 0.003846276658463,
                                     0.001704523106360,
                                     0.004819186789386,
                                     0.001553953681336,
                                    -0.000083971410109,
                                     0.000083971410109)
                            .shiftFromMassCenter(-lowerarm_rCOM,lowerarmMass)));
    lowerarm_rInfo.addDecoration(Transform(
                    Rotation(BodyRotationSequence,-0.05,XAxis,0.45,ZAxis),
                    Vec3(0.053,-0.111,0.01)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(0.95,1.15,0.95)))
            .setColor(White).setOpacity(1));
    lowerarm_rInfo.addDecoration(Vec3(0),originMarker);

    const Vec3 lowerarm_lCOM(0.031656703591848,-0.089369993258598,-0.017231110378866);
    Body lowerarm_lInfo(MassProperties(lowerarmMass, lowerarm_lCOM,
                            Inertia( 0.003846276658463,
                                     0.001704523106360,
                                     0.004819186789386,
                                     0.001553953681336,
                                     0.000083971410109,
                                    -0.000083971410109)
                            .shiftFromMassCenter(-lowerarm_lCOM,lowerarmMass)));
    lowerarm_lInfo.addDecoration(Transform(
                    Rotation(BodyRotationSequence,0.05,XAxis,0.45,ZAxis),
                    Vec3(0.053,-0.111,-0.01)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(0.95,1.15,0.95)))
            .setColor(White).setOpacity(1));
    lowerarm_lInfo.addDecoration(Vec3(0),originMarker);

    // Hand.
    // -----
    // Some signs differ left to right.
    const Real handMass = 0.340742469583528;
    const Vec3 hand_rCOM(0.031681557027587,-0.041582042351409,-0.008872831097566);
    Body hand_rInfo(MassProperties(handMass, hand_rCOM,
                            Inertia( 0.000294382529694,
                                     0.000262531305170,
                                     0.000326233754218,
                                     0.000061772071805,
                                     0.000054050562829,
                                    -0.000036677167634)
                            .shiftFromMassCenter(-hand_rCOM,handMass)));
    hand_rInfo.addDecoration(Transform(
                    Rotation(0.8,ZAxis),
                    Vec3(0.025,-0.025,0)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(0.7,0.3,0.7)))
            .setColor(White).setOpacity(1));
    hand_rInfo.addDecoration(Vec3(0),originMarker);

    const Vec3 hand_lCOM(0.031681557027587,-0.041582042351409, 0.008872831097566);
    Body hand_lInfo(MassProperties(handMass, hand_lCOM,
                            Inertia( 0.000294382529694,
                                     0.000262531305170,
                                     0.000326233754218,
                                     0.000061772071805,
                                    -0.000054050562829,
                                     0.000036677167634)
                            .shiftFromMassCenter(-hand_lCOM,handMass)));
    hand_lInfo.addDecoration(Transform(
                    Rotation(0.8,ZAxis),
                    Vec3(0.025,-0.025,0)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(0.7,0.3,0.7)))
            .setColor(White).setOpacity(1));
    hand_lInfo.addDecoration(Vec3(0),originMarker);

    //--------------------------------------------------------------------------
    //                         Mobilized Bodies
    //--------------------------------------------------------------------------
    // Trunk.
    // ------
    m_bodies[trunk] = MobilizedBody::Weld(
        m_matter.updGround(), Vec3(0, 1.5, 0),
        trunkInfo,          Vec3(0));

    // Neck angles are: q0=extension (about z), q1=bending (x), q2=rotation (y).
    m_bodies[head] = MobilizedBody::Gimbal(
        m_bodies[trunk], Transform(Rzxy, Vec3(0.010143822053248,0.222711680750785,0)),
        headInfo, Rzxy);
    addMobilityLinearStop(m_bodies[head], 0, -80, 50); // extension
    addMobilityLinearStop(m_bodies[head], 1, -60, 50); // bending
    addMobilityLinearStop(m_bodies[head], 2, -80, 80); // rotation

    // Right arm.
    //-----------
    // Shoulder angles are q0=flexion, q1=adduction, q2=rotation
    m_bodies[upperarm_r] = MobilizedBody::Gimbal(
        m_bodies[trunk], Transform(Rzxy,
                Vec3(-0.023921136233947,0.079313926824158,0.164710443657016)),
        upperarm_rInfo, Rzxy);
    addMobilityLinearStop(m_bodies[upperarm_r], 0, -80, 160); // flexion
    addMobilityLinearStop(m_bodies[upperarm_r], 1, -45, 45); // adduction
    addMobilityLinearStop(m_bodies[upperarm_r], 2, -20, 20); // rotation

    // Elbow angles are q0=flexion, q1=rotation
    m_bodies[lowerarm_r] = MobilizedBody::Universal(
        m_bodies[upperarm_r], Transform(Rotation(-Pi/2,YAxis),
                Vec3(0.033488432642100,-0.239093933565560,0.117718445964118)),
        lowerarm_rInfo, Rotation(-Pi/2,YAxis));
    addMobilityLinearStop(m_bodies[lowerarm_r], 0, 0, 120); // flexion
    addMobilityLinearStop(m_bodies[lowerarm_r], 1, -90, 40); // rotation

    m_bodies[hand_r] = MobilizedBody::Weld(
        m_bodies[lowerarm_r], Vec3(0.110610146564261,-0.232157950907188,0.014613941476432),
        hand_rInfo, Vec3(0));

    // Left arm.
    //----------
    m_bodies[upperarm_l] = MobilizedBody::Gimbal(
        m_bodies[trunk], Transform(Rzmxmy,
                Vec3(-0.023921136233947,0.079313926824158,-0.164710443657016)),
        upperarm_lInfo, Rzmxmy);
    addMobilityLinearStop(m_bodies[upperarm_l], 0, -80, 160); // flexion
    addMobilityLinearStop(m_bodies[upperarm_l], 1, -45, 45); // adduction
    addMobilityLinearStop(m_bodies[upperarm_l], 2, -20, 20); // rotation

    m_bodies[lowerarm_l] = MobilizedBody::Universal(
        m_bodies[upperarm_l], Transform(
                Rotation(BodyRotationSequence, Pi/2,YAxis, Pi,ZAxis),
                Vec3(0.033488432642100,-0.239093933565560,-0.117718445964118)),
        lowerarm_lInfo, Rotation(BodyRotationSequence, Pi/2,YAxis, Pi,ZAxis));
    addMobilityLinearStop(m_bodies[lowerarm_l], 0, 0, 120); // flexion
    addMobilityLinearStop(m_bodies[lowerarm_l], 1, -90, 40); // rotation

    m_bodies[hand_l] = MobilizedBody::Weld(
        m_bodies[lowerarm_l], Vec3(0.110610146564261,-0.232157950907188,-0.014613941476432),
        hand_lInfo, Vec3(0));
}

void UpperBody::fillInCoordinateMap(const State& s)
{
    // So all lines below fit within 80 columns:
    std::map<Coordinate, std::pair<QIndex, UIndex> >& coords = m_coordinates;

    // Upper body.
    coords[neck_extension] = getQandUIndices(s, m_bodies[head], 0);
    coords[neck_bending] = getQandUIndices(s, m_bodies[head], 1);
    coords[neck_rotation] = getQandUIndices(s, m_bodies[head], 2);

    // right.
    coords[shoulder_r_flexion]   = getQandUIndices(s, m_bodies[upperarm_r], 0);
    coords[shoulder_r_adduction] = getQandUIndices(s, m_bodies[upperarm_r], 1);
    coords[shoulder_r_rotation]  = getQandUIndices(s, m_bodies[upperarm_r], 2);

    coords[elbow_r_flexion]  = getQandUIndices(s, m_bodies[lowerarm_r], 0);
    coords[elbow_r_rotation] = getQandUIndices(s, m_bodies[lowerarm_r], 1);

    // left.
    coords[shoulder_l_flexion]   = getQandUIndices(s, m_bodies[upperarm_l], 0);
    coords[shoulder_l_adduction] = getQandUIndices(s, m_bodies[upperarm_l], 1);
    coords[shoulder_l_rotation]  = getQandUIndices(s, m_bodies[upperarm_l], 2);

    coords[elbow_l_flexion]  = getQandUIndices(s, m_bodies[lowerarm_l], 0);
    coords[elbow_l_rotation] = getQandUIndices(s, m_bodies[lowerarm_l], 1);
}

//==============================================================================
// ReachingAndGravityCompensation
//==============================================================================
class ReachingAndGravityCompensation :
    public SimTK::Force::Custom::Implementation
{
public:

    /// The left hand reaches for a target location, and the controller
    /// compensates for gravity.
    ReachingAndGravityCompensation(const UpperBody& system) :
        m_system(system), m_matter(system.getMatterSubsystem()) {}

    void calcForce(const SimTK::State&                state,
                   SimTK::Vector_<SimTK::SpatialVec>& bodyForces,
                   SimTK::Vector_<SimTK::Vec3>&       particleForces,
                   SimTK::Vector&                     mobilityForces) const
                   OVERRIDE_11;

    Real calcPotentialEnergy(const SimTK::State& state) const OVERRIDE_11
    { return 0; }

private:

    const UpperBody& m_system;
    const SimbodyMatterSubsystem& m_matter;

};

void ReachingAndGravityCompensation::calcForce(
               const SimTK::State&                state,
               SimTK::Vector_<SimTK::SpatialVec>& bodyForces,
               SimTK::Vector_<SimTK::Vec3>&       particleForces,
               SimTK::Vector&                     mobilityForces) const
{
    m_matter.multiplyBySystemJacobianTranspose(state,
            -m_system.getGravity().getBodyForces(state),
            mobilityForces);
}

//==============================================================================
// UserInputHandler
//==============================================================================
/// This is a periodic event handler that interrupts the simulation on a
/// regular basis to poll the InputSilo for user input.
class UserInputHandler : public PeriodicEventHandler {
public:
    UserInputHandler(Visualizer::InputSilo& silo, Real interval) :
        PeriodicEventHandler(interval), m_silo(silo) {}
    void handleEvent(State& state, Real accuracy,
            bool& shouldTerminate) const OVERRIDE_11
    {
        while (m_silo.isAnyUserInput()) {
            unsigned key, modifiers;
            while (m_silo.takeKeyHit(key, modifiers)) {
                if (key == Visualizer::InputListener::KeyEsc) {
                    shouldTerminate = true;
                    m_silo.clear();
                    return;
                }
            }
        }
    }
private:
    Visualizer::InputSilo& m_silo;
};

//==============================================================================
// OutputReporter
//==============================================================================
/// This is a periodic event handler that interrupts the simulation on a
/// regular basis so that we can display information.
class OutputReporter : public PeriodicEventReporter {
public:
    OutputReporter(const UpperBody& system, Real interval) :
        PeriodicEventReporter(interval), m_system(system) {}
    void handleEvent(const State& state) const OVERRIDE_11 {}
private:
    const UpperBody& m_system;
};

//==============================================================================
// main function
//==============================================================================
int main(int argc, char **argv)
{
    try {

        // Set some options.
        const double duration = 1.5; // seconds.
        const double realTimeScale = 1.0; // unitless.
        const double frameRate = 30.0; // frames per second.

        // Create the system.
        UpperBody system;
        system.updMatterSubsystem().setShowDefaultGeometry(false);

        // Set up visualizer and event handlers.
        SimTK::Visualizer viz(system);
        viz.setMode(Visualizer::RealTime);
        viz.setRealTimeScale(realTimeScale);
        viz.setShowFrameRate(true);
        viz.setShowSimTime(true);
        SimTK::Visualizer::InputSilo* userInput = new Visualizer::InputSilo();
        viz.addInputListener(userInput);
        system.addEventHandler(new UserInputHandler(*userInput, 0.1));
        system.addEventReporter(new OutputReporter(system, .01));
        system.addEventReporter(
                new Visualizer::Reporter(viz, realTimeScale / frameRate));

        // Display message on the screen about how to start simulation.
        DecorativeText help("Any input to start; ESC to quit");
        help.setIsScreenText(true);
        viz.addDecoration(MobilizedBodyIndex(0), Vec3(0), help);

        // Add the controller.
        ReachingAndGravityCompensation* controller =
            new ReachingAndGravityCompensation(system);
        // Force::Custom takes ownership over controller.
        Force::Custom control(system.updForceSubsystem(), controller);

        // Initialize the system and other related classes.
        State s;
        system.initialize(s);
        SemiExplicitEuler2Integrator integ(system);
        integ.setAccuracy(0.001);
        TimeStepper ts(system, integ);
        ts.initialize(s);
        viz.report(ts.getState());

        /* TODO
           printf("Act: u\n");
           for (int i=0; i < NumActuators; ++i) {
           printf("%2d: %d\n", i, int(system.getUIndex(Biped::Coordinate(i))));
           }
           */

        printf("Hit ENTER to simulate ... (ESC to quit)\n");
        userInput->waitForAnyUserInput();
        userInput->clear();

        const double startCPU  = cpuTime(), startTime = realTime();

        // Simulate.
        ts.stepTo(duration);

        std::cout << "CPU time: " << cpuTime() - startCPU << " seconds. "
                  << "Real time: " << realTime() - startTime << " seconds."
                  << std::endl;

    } catch (const std::exception& e) {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
