/* -------------------------------------------------------------------------- *
 *                 Simbody(tm) Example: SIMBICON                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-12 Stanford University and the Authors.        *
 * Authors: Chris Dembia, Jack Wang                                           *
 * Contributors: Michael Sherman, Tim Dorn                                    *
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

/*

SIMBICON stands for Simple Biped Controller. It is described at
https://www.cs.ubc.ca/~van/papers/Simbicon.htm.


Brief description of the controller
===================================

The controller uses a "finite state machine" / "pose control graph": the biped
is always in some state, and all states have a target pose. The number of
states is typically 4.  The joints are torque-actuated, and the torques are
given by proportional-derivative (PD) control laws that minimize the error of a
given joint angle from that in the target pose:

    \tau = k_p (\theta_{des} - \theta) - k_d \dot{\theta}

There are some exceptions to this general control scheme. Namely, the target
angles for the torso and "swing-leg femur" are specified relative to the global
frame. These two angles are used to balance of the biped.

The parameters of the controller are used to define the motion that the biped
executes. The motions the controller can achieve consist of walking, running,
and variants on both (e.g, 'high-step walk'). The paramters are:

for each state:
    * state dwell duration
    * position balance feedback coefficient
    * position balance feedback coefficient, lateral (for 3D gait)
    * velocity balance feedback coefficient
    * velocity balance feedback coefficient, lateral (for 3D gait)
    * torso target angle
    * swing-hip target angle
    * swing-knee target angle
    * swing-ankle target angle
    * stance-knee target angle
    * stance-ankle target angle

See [1], particularly Table 1 of [1], for more information.


How to navigate this file
=========================
To see how X happens, look in method Y:

X                                        Y
---------------------------------------  ---------------------------------------
construction of the model                Biped::Biped()
control laws                             SIMBICON::addInPDControl()
                                         SIMBICON::addInBalanceControl()
logic to change state in state machine   SIMBICON::updateSIMBICONState()


Notes
=====

I (Chris Dembia) wrote this code pretty much by copying code that Michael
Sherman gave me. The code that Michael Sherman gave me was originally written
by Jack Wang and Tim Dorn.

Yin, KangKang, Kevin Loken, and Michiel van de Panne. "SIMBICON: Simple biped
locomotion control." ACM Transactions on Graphics (TOG). Vol. 26. No. 3. ACM,
2007.

*/

// Interval for updating the SIMBICON state (not same as SimTK State) of the
// biped.
#define SIMBICON_STATE_UPDATE_STEPSIZE 0.0005

#include "Simbody.h"

using namespace SimTK;

const Real D2R = Pi / 180.0;

//==============================================================================
// BIPED
//==============================================================================
/// The MultibodySystem that we will control with the SIMBICON controller
/// (though the system is independent of the controller, and coulc conceivably
/// be controlled by a different controller).
class Biped : public MultibodySystem {
public:
    Biped();

    const SimbodyMatterSubsystem& getMatterSubsystem() const {return m_matter;}
    SimbodyMatterSubsystem& updMatterSubsystem() {return m_matter;}

    const GeneralForceSubsystem& getForceSubsystem() const {return m_forces;}
    GeneralForceSubsystem& updForceSubsystem() {return m_forces;}

    /// Realizes the topology and model, and uses the resulting initial state
    /// to perform further internal initialization.
    void initialize(State& state) {
        state = realizeTopology();
        realizeModel(state);

        // Initialization specific to Biped.
        fillInCoordinateMap(state);
    }

    /// Populates a map that allow users of Biped to access coordinates (Q) and
    /// speeds (U) with descriptive Enums (e.g., neck_extension) so the user
    /// (e.g., SIMBICON) needn't know how State indices correspond to
    /// coordinates.
    void fillInCoordinateMap(const State& s);

    void setTrunkOriginPosition(State& s, Vec3 posInGround) const {
        m_bodies.at(trunk).setQToFitTranslation(s, posInGround);
    }
    void setTrunkOriginVelocity(State& s, Vec3 velocityInGround) const {
        m_bodies.at(trunk).setUToFitLinearVelocity(s, velocityInGround);
    }

    /// param[out] fLeft vertical contact force on left foot.
    /// param[out] fRight vertical contact force on right foot.
    void findContactForces(const State& s, Real& fLeft, Real& fRight) const;

    /// param[out] left left foot is in contact with the ground.
    /// param[out] right right foot is in contact with the ground.
    void findContactStatus(const State& s, bool& left, bool& right) const;

    enum Segment {
        trunk, // called 'torso' in paper.
        head,
        pelvis,
        upperarm_r,
        lowerarm_r,
        hand_r,
        upperarm_l,
        lowerarm_l,
        hand_l,
        thigh_r,
        shank_r,
        foot_r,
        toes_r,
        thigh_l,
        shank_l,
        foot_l,
        toes_l
    };

    enum Coordinate {
        neck_extension = 0,
        neck_bending = 1,
        neck_rotation = 2,
        back_tilt = 3,
        back_list = 4,
        back_rotation = 5,
        shoulder_r_flexion = 6,
        shoulder_r_adduction = 7,
        shoulder_r_rotation = 8,
        elbow_r_flexion = 9,
        elbow_r_rotation = 10,
        shoulder_l_flexion = 11,
        shoulder_l_adduction = 12,
        shoulder_l_rotation = 13,
        elbow_l_flexion = 14,
        elbow_l_rotation = 15,
        hip_r_adduction = 16,
        hip_r_flexion = 17,
        hip_r_rotation = 18,
        knee_r_extension = 19,
        ankle_r_inversion = 20,
        ankle_r_dorsiflexion = 21,
        mtp_r_dorsiflexion = 22,
        hip_l_adduction = 23,
        hip_l_flexion = 24,
        hip_l_rotation = 25,
        knee_l_extension = 26,
        ankle_l_inversion = 27,
        ankle_l_dorsiflexion = 28,
        mtp_l_dorsiflexion = 29,
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

    /// Add in a generalized force into the proper place in `mobForces`,
    /// as determined by `coord`.
    void addInForce(const Coordinate coord, const Real force,
            Vector& mobForces) const {
        Real forceAccountingForSymmetry = force;
        if (coord == hip_r_adduction) {
            // Negate so that we actually apply an adduction torque (not an
            // abduction torque).
            forceAccountingForSymmetry *= -1;
        }
        mobForces[getUIndex(coord)] += forceAccountingForSymmetry;
    }

    UnitVec3 getNormalToSagittalPlane(const State& s) const
    {
        Vec3 bodyZInGround = Vec3(m_bodies.at(pelvis).getBodyRotation(s).z());
        // Project onto a horizontal plane.
        bodyZInGround[YAxis] = 0;
        // Make into a unit vector.
        return UnitVec3(bodyZInGround);
    }

    UnitVec3 getNormalToCoronalPlane(const State& s) const
    {
        Vec3 bodyXInGround = Vec3(m_bodies.at(pelvis).getBodyRotation(s).x());
        // Project onto a horizontal plane.
        bodyXInGround[YAxis] = 0;
        // Make into a unit vector.
        return UnitVec3(bodyXInGround);
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
            Real stopStiffness=1000 / D2R, Real stopDissipation=1.0)
    {
        Force::MobilityLinearStop(m_forces, mobod, MobilizerQIndex(idx),
                stopStiffness, stopDissipation,
                lowerBound * D2R, upperBound * D2R);
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

    bool isRightFoot(const MobilizedBody& mobod) const {
        return     mobod.isSameMobilizedBody(m_bodies.at(foot_r))
                || mobod.isSameMobilizedBody(m_bodies.at(toes_r));
    }
    bool isLeftFoot(const MobilizedBody& mobod) const {
        return     mobod.isSameMobilizedBody(m_bodies.at(foot_l))
                || mobod.isSameMobilizedBody(m_bodies.at(toes_l));
    }

    // Subsystems.
    SimbodyMatterSubsystem       m_matter;
    GeneralForceSubsystem        m_forces;
    ContactTrackerSubsystem      m_tracker;
    CompliantContactSubsystem    m_contact;

    /// The MobilizedBody's that make up the Biped.
    std::map<Segment, MobilizedBody> m_bodies;

    /// The indices of each coordinate's Q and U in the State vector.
    std::map<Coordinate, std::pair<QIndex, UIndex> > m_coordinates;
};

//==============================================================================
// SIMBICON
//==============================================================================

double clamp(double x, double max) {
    assert(max > 0);
    if (x > max) return max;
    if (x < -max) return -max;
    return x;
}

/// The actual controller that specifies the torques to apply to the Biped.
class SIMBICON : public Force::Custom::Implementation {
public:

    /// The default arguments define the motion that the controller will try to
    /// execute. Only symmetrical 4-state motions are permitted. The first
    /// element of each Vec2 is for states 0 and 2. The second element is for
    /// states 1 and 3.
    ///
    ///   * deltaT: state dwell duration
    ///   * cd: position balance feedback coefficient
    ///   * cdLat: position balance feedback coefficient, lateral (for 3D gait)
    ///   * cv: velocity balance feedback coefficient
    ///   * cvLat: elocity balance feedback coefficient, lateral (for 3D gait)
    ///   * tor: torso target angle
    ///   * swh: swing-hip target angle
    ///   * swk: swing-knee target angle
    ///   * swa: swing-ankle target angle
    ///   * stk: stance-knee target angle
    ///   * sta: stance-ankle target angle
    ///
    /// See the SIMBICON paper for information about these parameters.
    ///
    SIMBICON(Biped& biped,
            Real minSIMBICONStateDuration=0.1,
            Vec2 deltaT=Vec2(0.30, NaN),
            Vec2 cd=Vec2(0.5, 0.5),
            Vec2 cdLat=Vec2(0.5, 0.5),
            Vec2 cv=Vec2(0.2, 0.2),
            Vec2 cvLat=Vec2(0.2, 0.2),
            Vec2 tor=Vec2(0.0, 0.0),
            Vec2 swh=Vec2(0.5, -0.1),
            Vec2 swk=Vec2(-1.1, -0.05),
            Vec2 swa=Vec2(0.6, 0.15),
            Vec2 stk=Vec2(-0.05, -0.1),
            Vec2 sta=Vec2(0.0, 0.0)
            );

    /// The possible states of the finite state machine.
    enum SIMBICONState {
        UNKNOWN = -1, // initial state at the beginning of a simulation.
        STATE0 = 0,
        STATE1 = 1,
        STATE2 = 2,
        STATE3 = 3
    };

    /// Allocates some state that is necessary for this controller.
    void realizeTopology(State& s) const
    {
        // NOTE: the next two assignments discard the const qualifier of this
        // method. There is a very specific reason for this! Sherm explained
        // this to me as follows:
        //
        // 1. The state has a cache. Even when the state is const, it's okay to
        // write to the cache. This is because the cache is something that is
        // computed from a (const) state. So long as the state is const, we are
        // able to regenerate the cache appropriately.
        //
        // 2. We can think of the model-building process as having a state and
        // a cache. The state is the things that make up the model, such as the
        // MobilizedBody's. The cache in this case is quantities computed from
        // this model-buidling state, such as the number of mobilized bodies.
        //
        // 3. Therefore, quantities computed about the model-building state do
        // not need to be const.
        //
        // 4. The indicies in the state of my two AutoUpdateDiscreteVariable's
        // are part of the model-building cache, and so I can remove the const
        // qualifier.
        const_cast<SIMBICON*>(this)->m_simbiconStateIndex =
            m_forces.allocateAutoUpdateDiscreteVariable(s,
                    Stage::Acceleration, new Value<SIMBICONState>(UNKNOWN),
                    Stage::Dynamics);

        const_cast<SIMBICON*>(this)->m_simbiconStateCacheIndex =
            m_forces.getDiscreteVarUpdateIndex(s, m_simbiconStateIndex);
    }

    const SIMBICONState getSIMBICONState(const State& s) const
    {
        const AbstractValue& m = m_forces.getDiscreteVariable(
                s, m_simbiconStateIndex);
        return Value<SIMBICONState>::downcast(m).get();
    }
    Real getSIMBICONStateStartTime(const State& s) const
    {
        return m_forces.getDiscreteVarLastUpdateTime(s,
                m_simbiconStateIndex);
    }

    /// The control consists of a few components. This method calls other
    /// methods that take care of these individual components.
    ///
    void calcForce(const State&         s,
                   Vector_<SpatialVec>& bodyForces,
                   Vector_<Vec3>&       particleForces,
                   Vector&              mobilityForces) const
                   OVERRIDE_11
    {
        updateSIMBICONState(s);
        addInPDControl(s, mobilityForces);
        addInBalanceControl(s, mobilityForces);
    }

    Real calcPotentialEnergy(const State& state) const OVERRIDE_11
    {
        return 0;
    }

    /// Used to identify gains for different parts of the body.
    enum GainGroup {
        neck,
        back,
        hip_flexion_adduction,
        hip_rotation,
        knee,
        arm_flexion_adduction,
        arm_rotation,
        ankle_flexion,
        ankle_inversion,
        toe
    };

private:

    /// Determine if the SIMBICON state of the state machine has changed.
    void updateSIMBICONState(const State& s) const;

    /// Compute torques defined by proportional-derivative (PD) control laws
    /// that track the target pose for the current SIMBICON state. We have
    /// such a control law for most (TODO) joints in the biped.
    void addInPDControl(const State& s, Vector& mobForces) const;

    /// Compute the torques used to balance the biped. We compute internal
    /// torques that achieve the tracking of target angles that are specified
    /// relative to the global frame.
    void addInBalanceControl(const State& s, Vector& mobForces) const;

    /// Update the cache variable.
    void setSIMBICONState(const State& s, SIMBICONState simbiconState) const
    {

        m_forces.updDiscreteVarUpdateValue(s, m_simbiconStateIndex) =
            Value<SIMBICONState>(simbiconState);
        m_forces.markDiscreteVarUpdateValueRealized(s, m_simbiconStateIndex);
    }

    /// Apply generalized force to a specific coord using a
    /// proportional-derivative control law that tracks thetad.
    /// The force is added into the appropriate place in `mobForces.`
    /// For convenience / inspection, the force is returned as well.
    Real coordPDControl(const State& s,
            const Biped::Coordinate coord, const GainGroup gainGroup,
            const Real thetad, Vector& mobForces) const
    {
        // Prepare quantities.
        Real kp = m_proportionalGains.at(gainGroup);
        Real kd = m_derivativeGains.at(gainGroup);
        Real q = m_biped.getQ(s, coord);
        Real u = m_biped.getU(s, coord);

        // PD control law:
        Real force = clamp(kp * (thetad - q) - kd * u, kp);

        // Update the proper entry in mobForces.
        m_biped.addInForce(coord, force, mobForces);

        return force;
    }

    Biped& m_biped;
    const GeneralForceSubsystem& m_forces;

    /// SIMBICON parameters.
    const Vec2 m_deltaT;
    const Vec2 m_cd;
    const Vec2 m_cdLat;
    const Vec2 m_cv;
    const Vec2 m_cvLat;
    const Vec2 m_tor;
    const Vec2 m_swh;
    const Vec2 m_swk;
    const Vec2 m_swa;
    const Vec2 m_stk;
    const Vec2 m_sta;

    std::map<GainGroup, Real> m_proportionalGains;
    std::map<GainGroup, Real> m_derivativeGains;

    const Real m_minSIMBICONStateDuration;

    // The SIMBICON controller itself has some state.
    // ----------------------------------------------
    // Written during the realizeTopology step.

    /// The SIMBICON state of the controller for the finite state machine
    /// (e.g., 0, 1, 2, or 3). We read the value of this variable using this
    //index.
    DiscreteVariableIndex m_simbiconStateIndex;

    // We set the value of the SIMBICON state using this index.
    CacheEntryIndex m_simbiconStateCacheIndex;


};

//==============================================================================
// OUTPUT REPORTER
//==============================================================================
/// This is a periodic event handler that we use to print information to the
/// visualizer during the simulation.
class OutputReporter : public DecorationGenerator {
public:
    OutputReporter(const Biped& biped, const SIMBICON* simbicon)
    :   m_biped(biped), m_simbicon(simbicon) {}

    void generateDecorations(const State&                state,
                             Array_<DecorativeGeometry>& geometry) OVERRIDE_11
    {
        DecorativeText text;
        text.setIsScreenText(true);
        std::ostringstream oss;

        bool lContact;
        bool rContact;
        m_biped.findContactStatus(state, lContact, rContact);

        // Directed distally.
        UnitVec3 trunkAxialDir = m_biped.getBody(Biped::trunk).getBodyRotation(state).y();
        // UnitVec3(YAxis) gives a vector perpendicular to ground, directed up.
        Real globalTrunkAngle = acos(dot(UnitVec3(YAxis), trunkAxialDir)) / D2R;

        // TODO

        oss << "SIMBICON state " << m_simbicon->getSIMBICONState(state) << 
            " LContact: " << lContact <<
            " RContact: " << rContact <<
            " global trunk angle: " << globalTrunkAngle
            << std::endl;

        text.setText(oss.str());
        geometry.push_back(text);
    }
private:
    const Biped& m_biped;
    const SIMBICON* m_simbicon;
};

//==============================================================================
// MAIN FUNCTION
//==============================================================================
int main(int argc, char **argv)
{
    // Parameters.
    // -----------

    // Set this < 1 for slow motion, > 1 for speedup.
    // 1.0 indicates real time.
    const Real realTimeFactor = 1.0;

    const Real framesPerSecond = 30.0;

    // Create the system.
    Biped biped;

    // Add controller to the force system. See Doxygen for Force::Custom.
    // The name of this Force is irrelevant.
    SIMBICON* simbicon = new SIMBICON(biped);
    Force::Custom simbicon1(biped.updForceSubsystem(), simbicon);

    // Start up the visualizer.
    Visualizer viz(biped);
    viz.setShowSimTime(true);
    viz.setShowFrameRate(true);

    // Update the visualization at a framerate of 1/interval.
    biped.addEventReporter(
            new Visualizer::Reporter(viz, realTimeFactor / framesPerSecond));

    // Initialize the system (this is our own method).
    State state;
    biped.initialize(state);
    biped.realize(state, Stage::Instance);

    // Set the initial conditions for the biped.
    biped.setTrunkOriginPosition(state, Vec3(0, 1.5, 0));
    // Give the biped some initial forward velocity.
    biped.setTrunkOriginVelocity(state, Vec3(1, 0, 0));

    // Show the biped.
    printf("Initial state show. Press ENTER to simulate.\n");
    viz.report(state);

    // Periodically report state-related quantities to the visualizer.
    viz.addDecorationGenerator(new OutputReporter(biped, simbicon));

    getchar();

    // Set up the numerical integration.
    SemiExplicitEuler2Integrator integ(biped);
    integ.setAccuracy(0.1);
    TimeStepper ts(biped, integ);
    ts.initialize(state);

    ts.stepTo(Infinity);

	return 0;
}


//==============================================================================
// BIPED
//==============================================================================
Biped::Biped()
    : m_matter(*this), m_forces(*this),
      m_tracker(*this), m_contact(*this, m_tracker)
{

    m_matter.setShowDefaultGeometry(false);

    //--------------------------------------------------------------------------
    //                          Constants, etc.
    //--------------------------------------------------------------------------

    // Properties of joint stops.
    const Real stopStiffness = 1000 / D2R; // N-m/deg -> N-m/radian
    const Real stopDissipation = 1;

    // Contact parameters.
    // -------------------
    // slide -> stick velocity.
    const Real transitionVelocity = .03;

    // Friction coefficients.
    const Real mu_s = 5;
    const Real mu_d = 5;
    const Real mu_v = 0;

    // Contact spheres attached to feet.
    const Real contactSphereRadius = 1.5 * .01;

    // Contact material: rubber.
    const Real rubber_density = 1100.;  // kg/m^3
    const Real rubber_young   = 0.01e9; // pascals (N/m)
    const Real rubber_poisson = 0.5;    // ratio
    const Real rubber_planestrain =
        ContactMaterial::calcPlaneStrainStiffness(rubber_young,rubber_poisson);
    const Real rubber_dissipation = /*0.005*/1;

    const ContactMaterial rubber(rubber_planestrain,rubber_dissipation,
                                   mu_s,mu_d,mu_v);

    // Contact material: concrete.
    const Real concrete_density = 2300.;  // kg/m^3
    const Real concrete_young   = 25e9;  // pascals (N/m)
    const Real concrete_poisson = 0.15;    // ratio
    const Real concrete_planestrain =
        ContactMaterial::calcPlaneStrainStiffness(concrete_young,concrete_poisson);
    const Real concrete_dissipation = 0.005;

    const ContactMaterial concrete(concrete_planestrain,concrete_dissipation,
                                   mu_s,mu_d,mu_v);

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
    //                          Gravity and ground contact
    //--------------------------------------------------------------------------

    // Gravity.
    Force::Gravity(m_forces, m_matter, -YAxis, 9.8066);

    // Contact with the ground.
    m_matter.updGround().updBody().addContactSurface(
            Transform(Rotation(-0.5 * Pi, ZAxis), Vec3(0)),
            ContactSurface(ContactGeometry::HalfSpace(), concrete));


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

    // Pelvis.
    // -------
    const Real pelvisMass = 13.924855817213411;
    const Vec3 pelvisCOM(0.036907589663647, -0.142772863495411, 0);
    Body pelvisInfo(MassProperties(pelvisMass, pelvisCOM,
                                   Inertia( 0.172382614643800,
                                            0.137961114411544,
                                            0.128551359933154,
                                           -0.010239461806632,
                                            0.001027963710884,
                                           -0.003286514395970)
                                   .shiftFromMassCenter(-pelvisCOM,pelvisMass)));
    pelvisInfo.addDecoration(Transform(Ryzx,Vec3(0.02,-0.1,0.0)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(3.3,1.5,2.4)))
            .setColor(White).setOpacity(1));
    pelvisInfo.addDecoration(Vec3(0),originMarker);

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

    // Thigh.
    // ------
    const Real thighMass = 8.082407914884000;
    const Vec3 thigh_rCOM(0,-0.178920728716523,0.001605747837523);
    Body thigh_rInfo(MassProperties(thighMass, thigh_rCOM,
                            Inertia( 0.116351777130644,
                                     0.030499980412887,
                                     0.122695077900275 )
                            .shiftFromMassCenter(-thigh_rCOM,thighMass)));
    thigh_rInfo.addDecoration(Transform(
                    Rotation(BodyRotationSequence, -.01,XAxis, -.012,ZAxis),
                    Vec3(-0.002,-0.205,0.003)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(1.3,1.9,1.3)))
            .setColor(White).setOpacity(1));
    thigh_rInfo.addDecoration(Vec3(0),originMarker);

    const Vec3 thigh_lCOM(0,-0.178920728716523,-0.001605747837523);
    Body thigh_lInfo(MassProperties(thighMass, thigh_lCOM,
                            Inertia( 0.116351777130644,
                                     0.030499980412887,
                                     0.122695077900275 )
                            .shiftFromMassCenter(-thigh_lCOM,thighMass)));
    thigh_lInfo.addDecoration(Transform(
                    Rotation(BodyRotationSequence, .01,XAxis, -.012,ZAxis),
                    Vec3(-0.002,-0.205,-0.003)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(1.3,1.9,1.3)))
            .setColor(White).setOpacity(1));
    thigh_lInfo.addDecoration(Vec3(0),originMarker);

    // Shank.
    // ------
    const Real shankMass = 3.222323418816000;
    const Vec3 shank_rCOM(0,-0.182765070363067,0.005552190835500);
    Body shank_rInfo(MassProperties(shankMass, shank_rCOM,
                            Inertia( 0.043804477493817,
                                     0.004432595936874,
                                     0.044412873014564 )
                            .shiftFromMassCenter(-shank_rCOM,shankMass)));
    shank_rInfo.addDecoration(Transform(
                    Rotation(0.030,XAxis),
                    Vec3(0.0,-0.21,-0.005)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(1.27,1.8,1.27)))
            .setColor(White).setOpacity(1));
    shank_rInfo.addDecoration(Vec3(0),originMarker);

    const Vec3 shank_lCOM(0,-0.182765070363067,-0.005552190835500);
    Body shank_lInfo(MassProperties(shankMass, shank_lCOM,
                            Inertia( 0.043804477493817,
                                     0.004432595936874,
                                     0.044412873014564 )
                            .shiftFromMassCenter(-shank_lCOM,shankMass)));
    shank_lInfo.addDecoration(Transform(
                    Rotation(-0.030,XAxis),
                    Vec3(0.0,-0.21,0.005)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(1.27,1.8,1.27)))
            .setColor(White).setOpacity(1));
    shank_lInfo.addDecoration(Vec3(0),originMarker);

    // Foot.
    // -----
    const Real footMass = 1.172905458264000;
    const Vec3 foot_rCOM(0.035606945567853,-0.051617802456029,-0.000574057583573);
    const Inertia foot_Ic(0.001313654113256, // central inertia
                          0.003659465029784,
                          0.003847129903106);
    const Real Ifac = 1; // for playing with foot inertia; should be 1
    Body foot_rInfo(MassProperties(footMass, foot_rCOM,
                      (Ifac*foot_Ic).shiftFromMassCenter(-foot_rCOM,footMass)));
    foot_rInfo.addDecoration(Vec3(0.052,-0.043,0.0),
        DecorativeBrick(blk.elementwiseMultiply(Vec3(1.65,0.6,0.8)))
            .setColor(White).setOpacity(1));
    foot_rInfo.addDecoration(Vec3(0),originMarker);

    const Vec3 foot_lCOM(0.035606945567853,-0.051617802456029,0.000574057583573);
    Body foot_lInfo(MassProperties(footMass, foot_lCOM,
                      (Ifac*foot_Ic).shiftFromMassCenter(-foot_lCOM,footMass)));
    foot_lInfo.addDecoration(Vec3(0.052,-0.043,0.0),
        DecorativeBrick(blk.elementwiseMultiply(Vec3(1.65,0.6,0.8)))
            .setColor(White).setOpacity(1));
    foot_lInfo.addDecoration(Vec3(0),originMarker);

    // Toes.
    // -----
    const Real toesMass = 0.20;
    // This inertia looks artificially large.
    const Inertia toes_Ic(0.087132432150,0.0174264864299,0.0871324321496);
    const Vec3 toes_rCOM(0.023716003435794,-0.001184334594970,-0.002484544347746);
    Body toes_rInfo(MassProperties(toesMass, toes_rCOM,
                            toes_Ic.shiftFromMassCenter(-toes_rCOM,toesMass)));
    toes_rInfo.addDecoration(Vec3(0.03,0.014,0.0),
        DecorativeBrick(blk.elementwiseMultiply(Vec3(0.6,0.3,0.8)))
            .setColor(White).setOpacity(1));
    toes_rInfo.addDecoration(Vec3(0),originMarker);

    const Vec3 toes_lCOM(0.023716003435794,-0.001184334594970,0.002484544347746);
    Body toes_lInfo(MassProperties(toesMass, toes_lCOM,
                            toes_Ic.shiftFromMassCenter(-toes_lCOM,toesMass)));
    toes_lInfo.addDecoration(Vec3(0.03,0.014,0.0),
        DecorativeBrick(blk.elementwiseMultiply(Vec3(0.6,0.3,0.8)))
            .setColor(White).setOpacity(1));
    toes_lInfo.addDecoration(Vec3(0),originMarker);

    //--------------------------------------------------------------------------
    //                         Contact at feet and toes
    //--------------------------------------------------------------------------
    const Vec2 rlatmed(.04,-.04);
    const Vec2 heelball(-.02,.125);

    ContactSurface contactBall(ContactGeometry::Sphere(contactSphereRadius),
                               rubber);

    // Use this clique for contact surfaces on the humanoid that you don't want
    // to have bump into one another. They will still contact Ground.
    const ContactCliqueId clique = ContactSurface::createNewContactClique();

    contactBall.joinClique(clique);
    DecorativeSphere contactArt(contactSphereRadius);
    contactArt.setColor(Magenta);

    const Real footsphereht = -.07;
    for (int x=0; x <= 1; ++x)
        for (int z=0; z <= 1; ++z) {
            const Vec3 rctr(heelball[x], footsphereht,  rlatmed[z]);
            const Vec3 lctr(heelball[x], footsphereht, -rlatmed[z]);
            foot_rInfo.addContactSurface(rctr, contactBall);
            foot_rInfo.addDecoration(rctr, contactArt);
            foot_lInfo.addContactSurface(lctr, contactBall);
            foot_lInfo.addDecoration(lctr, contactArt);
        }

    // Balls just at toe tips.
    const Real toetip = .06;
    for (int z=0; z <= 1; ++z) {
        const Vec3 rctr(toetip, 0,  rlatmed[z]);
        const Vec3 lctr(toetip, 0, -rlatmed[z]);
        toes_rInfo.addContactSurface(rctr, contactBall);
        toes_rInfo.addDecoration(rctr, contactArt);
        toes_lInfo.addContactSurface(lctr, contactBall);
        toes_lInfo.addDecoration(lctr, contactArt);
    }

    //--------------------------------------------------------------------------
    //                         Mobilized Bodies
    //--------------------------------------------------------------------------
    // Trunk.
    // ------
    m_bodies[trunk] = MobilizedBody::Free(
        m_matter.updGround(), Vec3(0),
        trunkInfo,          Vec3(0));

    // Neck angles are: q0=extension (about z), q1=bending (x), q2=rotation (y).
    m_bodies[head] = MobilizedBody::Gimbal(
        m_bodies[trunk], Transform(Rzxy, Vec3(0.010143822053248,0.222711680750785,0)),
        headInfo, Rzxy);
    addMobilityLinearStop(m_bodies[head], 0, -80, 50); // extension
    addMobilityLinearStop(m_bodies[head], 1, -60, 50); // bending
    addMobilityLinearStop(m_bodies[head], 2, -80, 80); // rotation

    // Back angles are: q0=tilt (about z), q1=list (x), q2=rotation (y).
    m_bodies[pelvis] = MobilizedBody::Gimbal(
        m_bodies[trunk], Transform(Rzxy, Vec3(-0.019360589663647,-0.220484136504589,0)),
        pelvisInfo, Rzxy);
    addMobilityLinearStop(m_bodies[pelvis], 0, -5, 10); // tilt
    addMobilityLinearStop(m_bodies[pelvis], 1, -5, 5); // list
    addMobilityLinearStop(m_bodies[pelvis], 2, -15, 15); // rotation

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


    // Right leg.
    //-----------
    // Hip angles are q0=flexion, q1=adduction, q2=rotation.
    m_bodies[thigh_r] = MobilizedBody::Gimbal(
        m_bodies[pelvis], Transform(Rzxy,
                Vec3(0.029343644095793,-0.180413783750097,0.117592252162477)),
        thigh_rInfo, Rzxy);
    addMobilityLinearStop(m_bodies[thigh_r], 0, -60, 165); // flexion
    addMobilityLinearStop(m_bodies[thigh_r], 1, -20, 20); // adduction
    addMobilityLinearStop(m_bodies[thigh_r], 2, -120, 20); // rotation

    // Knee angle is q0=extension
    m_bodies[shank_r] = MobilizedBody::Pin(
        m_bodies[thigh_r], Vec3(-0.005,-0.416780050422019,0.004172557002023),
        shank_rInfo, Vec3(0));
    addMobilityLinearStop(m_bodies[shank_r], 0, -165, 0); // extension

    // Ankle angles are q0=dorsiflexion, q1=inversion.
    m_bodies[foot_r] = MobilizedBody::Universal(
        m_bodies[shank_r], Transform(
            Rotation(BodyRotationSequence,Pi/2,XAxis,Pi,YAxis,Pi/2,ZAxis),
            Vec3(0,-0.420937226867266,-0.011971751580927)),
        foot_rInfo,
            Rotation(BodyRotationSequence,Pi/2,XAxis,Pi,YAxis,Pi/2,ZAxis));
    addMobilityLinearStop(m_bodies[foot_r], 0, -50, 30); // dorsiflexion
    addMobilityLinearStop(m_bodies[foot_r], 1, -2, 35); // inversion

    // Toe angle is q0=dorsiflexion
    m_bodies[toes_r] = MobilizedBody::Pin(
        m_bodies[foot_r], Vec3(0.134331942132059,-0.071956467861059,-0.000000513235827),
        toes_rInfo, Vec3(0));
    addMobilityLinearStop(m_bodies[toes_r], 0, 0, 30); // dorsiflexion

    // Left leg.
    //----------
    m_bodies[thigh_l] = MobilizedBody::Gimbal(
        m_bodies[pelvis], Transform(Rzmxmy,
                Vec3(0.029343644095793,-0.180413783750097,-0.117592252162477)),
        thigh_lInfo, Rzmxmy);
    addMobilityLinearStop(m_bodies[thigh_l], 0, -60, 165); // flexion
    addMobilityLinearStop(m_bodies[thigh_l], 1, -20, 20); // adduction
    addMobilityLinearStop(m_bodies[thigh_l], 2, -120, 20); // rotation

    m_bodies[shank_l] = MobilizedBody::Pin(
        m_bodies[thigh_l], Vec3(-0.005,-0.416780050422019,-0.004172557002023),
        shank_lInfo, Vec3(0));
    addMobilityLinearStop(m_bodies[shank_l], 0, -165, 0); // extension

    m_bodies[foot_l] = MobilizedBody::Universal(
        m_bodies[shank_l], Transform(
            Rotation(BodyRotationSequence,-Pi/2,XAxis,Pi,YAxis,-Pi/2,ZAxis),
            Vec3(0,-0.420937226867266,0.011971751580927)),
        foot_lInfo,
            Rotation(BodyRotationSequence,-Pi/2,XAxis,Pi,YAxis,-Pi/2,ZAxis));
    addMobilityLinearStop(m_bodies[foot_l], 0, -50, 30); // dorsiflexion
    addMobilityLinearStop(m_bodies[foot_l], 1, -2, 35); // inversion

    m_bodies[toes_l] = MobilizedBody::Pin(
        m_bodies[foot_l], Vec3(0.134331942132059,-0.071956467861059,0.000000513235827),
        toes_lInfo, Vec3(0));
    addMobilityLinearStop(m_bodies[toes_l], 0, 0, 30); // dorsiflexion

}

void Biped::fillInCoordinateMap(const State& s)
{
    // So all lines below fit within 80 columns:
    std::map<Coordinate, std::pair<QIndex, UIndex> >& coords = m_coordinates;

    // Upper body.
    coords[neck_extension] = getQandUIndices(s, m_bodies[head], 0);
    coords[neck_bending] = getQandUIndices(s, m_bodies[head], 1);
    coords[neck_rotation] = getQandUIndices(s, m_bodies[head], 2);

    coords[back_tilt]     = getQandUIndices(s, m_bodies[pelvis], 0);
    coords[back_list]     = getQandUIndices(s, m_bodies[pelvis], 1);
    coords[back_rotation] = getQandUIndices(s, m_bodies[pelvis], 2);

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

    // Lower body.
    // right.
    coords[hip_r_flexion]    = getQandUIndices(s, m_bodies[thigh_r], 0);
    coords[hip_r_adduction]  = getQandUIndices(s, m_bodies[thigh_r], 1);
    coords[hip_r_rotation]   = getQandUIndices(s, m_bodies[thigh_r], 2);

    coords[knee_r_extension] = getQandUIndices(s, m_bodies[shank_r], 0);

    coords[ankle_r_dorsiflexion] = getQandUIndices(s, m_bodies[foot_r], 0);
    coords[ankle_r_inversion]    = getQandUIndices(s, m_bodies[foot_r], 1);

    coords[mtp_r_dorsiflexion] = getQandUIndices(s, m_bodies[toes_r], 0);

    // left.
    coords[hip_l_flexion]    = getQandUIndices(s, m_bodies[thigh_l], 0);
    coords[hip_l_adduction]  = getQandUIndices(s, m_bodies[thigh_l], 1);
    coords[hip_l_rotation]   = getQandUIndices(s, m_bodies[thigh_l], 2);

    coords[knee_l_extension] = getQandUIndices(s, m_bodies[shank_l], 0);

    coords[ankle_l_dorsiflexion] = getQandUIndices(s, m_bodies[foot_l], 0);
    coords[ankle_l_inversion]    = getQandUIndices(s, m_bodies[foot_l], 1);

    coords[mtp_l_dorsiflexion] = getQandUIndices(s, m_bodies[toes_l], 0);
}

void Biped::findContactForces(const State& s, Real& fLeft, Real& fRight) const
    {
        const unsigned int nContacts = m_contact.getNumContactForces(s);
        const ContactSnapshot& snapshot = m_tracker.getActiveContacts(s);

        fLeft = 0;
        fRight = 0;
        for (unsigned int i = 0; i < nContacts; ++i)
        {
            const ContactForce& force = m_contact.getContactForce(s, i);
            const ContactId id = force.getContactId();
            assert(snapshot.hasContact(id));
            const Contact& contact = snapshot.getContactById(id);
            const MobilizedBody& b1 = m_tracker.getMobilizedBody(contact.getSurface1());
            const MobilizedBody& b2 = m_tracker.getMobilizedBody(contact.getSurface2());
            const bool left = isLeftFoot(b1) || isLeftFoot(b2);
            const bool right = isRightFoot(b1) || isRightFoot(b2);
            if (left)
            {
                fLeft += force.getForceOnSurface2()[1].norm();
            }
            else if (right)
            {
                fRight += force.getForceOnSurface2()[1].norm();
            }
        }
    }

void Biped::findContactStatus(const State& s, bool& left, bool& right) const
{
    Real fLeft;
    Real fRight;
    findContactForces(s, fLeft, fRight);
    left = fLeft > 0;
    right = fRight > 0;
}


//==============================================================================
// SIMBICON
//==============================================================================
Real criticallyDampedDerivativeGain(Real proportionalGain) {
    return 2.0 * std::sqrt(proportionalGain);
}

SIMBICON::SIMBICON(Biped& biped,
        Real minSIMBICONStateDuration,
        Vec2 deltaT, Vec2 cd, Vec2 cdLat, Vec2 cv,
        Vec2 cvLat, Vec2 tor, Vec2 swh, Vec2 swk, Vec2 swa, Vec2 stk, Vec2 sta)
    : m_biped(biped), m_forces(m_biped.getForceSubsystem()),
      m_minSIMBICONStateDuration(minSIMBICONStateDuration),
      m_deltaT(deltaT), m_cd(cd),
      m_cdLat(cdLat), m_cv(cv), m_cvLat(cvLat), m_tor(tor), m_swh(swh),
      m_swk(swk), m_swa(swa), m_stk(stk), m_sta(sta)

{

    // Hook up with the model.
    // -----------------------
/*
            */

    // Proportional (position) gains (kp).
    m_proportionalGains[neck] = 100;
    m_proportionalGains[back] = 300;
    m_proportionalGains[hip_flexion_adduction] = 1000;
    m_proportionalGains[hip_rotation] = 300;
    m_proportionalGains[knee] = 300;
    m_proportionalGains[arm_flexion_adduction] = 300;
    m_proportionalGains[arm_rotation] = 300;
    m_proportionalGains[ankle_flexion] = 300;
    m_proportionalGains[ankle_inversion] = 30;
    m_proportionalGains[toe] = 30;

    // Derivative (speed) gains; mostly chosen for critical damping.
    m_derivativeGains[neck] = criticallyDampedDerivativeGain(
            m_proportionalGains[neck]);
    m_derivativeGains[back] = criticallyDampedDerivativeGain(
            m_proportionalGains[back]);
    // Overdamped.
    m_derivativeGains[hip_flexion_adduction] = 100;
    m_derivativeGains[hip_rotation] = criticallyDampedDerivativeGain(
            m_proportionalGains[hip_rotation]);
    m_derivativeGains[knee] = criticallyDampedDerivativeGain(
            m_proportionalGains[knee]);
    m_derivativeGains[arm_flexion_adduction] = criticallyDampedDerivativeGain(
            m_proportionalGains[arm_flexion_adduction]);
    m_derivativeGains[arm_rotation] = criticallyDampedDerivativeGain(
            m_proportionalGains[arm_rotation]);
    m_derivativeGains[ankle_flexion] = criticallyDampedDerivativeGain(
            m_proportionalGains[ankle_flexion]);
    m_derivativeGains[ankle_inversion] = criticallyDampedDerivativeGain(
            m_proportionalGains[ankle_inversion]);
    m_derivativeGains[toe] = criticallyDampedDerivativeGain(
            m_proportionalGains[toe]);
}

void SIMBICON::updateSIMBICONState(const State& s) const
{
    bool lContact;
    bool rContact;
    m_biped.findContactStatus(s, lContact, rContact);

    const SIMBICONState simbiconState = getSIMBICONState(s);
    const Real duration = s.getTime() - getSIMBICONStateStartTime(s);

    // simbiconState stateIdx
    // ------------- --------
    // 0             0
    // 1             1
    // 2             0
    // 3             1
    const int stateIdx = simbiconState % 2;

    switch (simbiconState)
    {
        // Neither foot has contacted the ground yet.
        case UNKNOWN:

            // Entering left stance.
            if (lContact) {
                 setSIMBICONState(s, STATE0);
            }
            // Entering right stance.
            else if (rContact) {
                 setSIMBICONState(s, STATE2);
            }
            break;

        // Left stance.
        case STATE0:

            // Stay in this state for \delta t seconds.
            if (duration > m_deltaT[stateIdx]) {
                setSIMBICONState(s, STATE1);
            }
            // Already entered right stance; skip STATE1.
            else if (rContact && duration > m_minSIMBICONStateDuration) {
                setSIMBICONState(s, STATE2);
            }
            break;

        // Right foot strike.
        case STATE1:

            // Stay in this state until the right foot makes contact.
            if (rContact && duration > m_minSIMBICONStateDuration) {
                setSIMBICONState(s, STATE2);
            }
            break;

        // Right stance.
        case STATE2:

            // Stay in this state for \delta t seconds.
            if (duration > m_deltaT[stateIdx]) {
                setSIMBICONState(s, STATE3);
            }
            // Already entered left stance; skip STATE3.
            else if (lContact && duration > m_minSIMBICONStateDuration) {
                setSIMBICONState(s, STATE0);
            }

        // Left foot strike.
        case STATE3:

            // Stay in this state until the left foot makes contact.
            if (lContact && duration > m_minSIMBICONStateDuration) {
                setSIMBICONState(s, STATE0); return;
            }
            break;
    }
}

void SIMBICON::addInPDControl(const State& s, Vector& mobForces) const
{

    // For most joints, track target theta of 0.0 degrees.
    // ===================================================
    coordPDControl(s, Biped::neck_extension, neck, 0.0, mobForces);
    coordPDControl(s, Biped::neck_bending, neck, 0.0, mobForces);
    coordPDControl(s, Biped::neck_rotation, neck, 0.0, mobForces);

    coordPDControl(s, Biped::back_tilt, back, 0.0, mobForces);
    coordPDControl(s, Biped::back_list, back, 0.0, mobForces);
    coordPDControl(s, Biped::back_rotation, back, 0.0, mobForces);

    coordPDControl(s, Biped::shoulder_r_flexion, arm_flexion_adduction, 0.0, mobForces);
    coordPDControl(s, Biped::shoulder_l_flexion, arm_flexion_adduction, 0.0, mobForces);
    coordPDControl(s, Biped::shoulder_r_adduction, arm_flexion_adduction, 0.0, mobForces);
    coordPDControl(s, Biped::shoulder_l_adduction, arm_flexion_adduction, 0.0, mobForces);
    coordPDControl(s, Biped::elbow_r_flexion, arm_flexion_adduction, 0.0, mobForces);
    coordPDControl(s, Biped::elbow_l_flexion, arm_flexion_adduction, 0.0, mobForces);

    coordPDControl(s, Biped::shoulder_r_rotation, arm_rotation, 0.0, mobForces);
    coordPDControl(s, Biped::shoulder_l_rotation, arm_rotation, 0.0, mobForces);
    coordPDControl(s, Biped::elbow_r_rotation, arm_rotation, 0.0, mobForces);
    coordPDControl(s, Biped::elbow_l_rotation, arm_rotation, 0.0, mobForces);
    coordPDControl(s, Biped::hip_r_rotation, hip_rotation, 0.0, mobForces);
    coordPDControl(s, Biped::hip_l_rotation, hip_rotation, 0.0, mobForces);

    coordPDControl(s, Biped::ankle_r_inversion, ankle_inversion, 0.0, mobForces);
    coordPDControl(s, Biped::ankle_l_inversion, ankle_inversion, 0.0, mobForces);

    coordPDControl(s, Biped::mtp_r_dorsiflexion, toe, 0.0, mobForces);
    coordPDControl(s, Biped::mtp_l_dorsiflexion, toe, 0.0, mobForces);

    // Deal with limbs whose target angle is affected by the state machine.
    // ====================================================================

    // Which leg is in stance?
    // -----------------------
    Biped::Coordinate swing_hip_flexion;
    Biped::Coordinate swing_hip_adduction;
    Biped::Coordinate swing_knee_extension;
    Biped::Coordinate swing_ankle_dorsiflexion;

    Biped::Coordinate stance_hip_flexion;
    Biped::Coordinate stance_hip_adduction;
    Biped::Coordinate stance_knee_extension;
    Biped::Coordinate stance_ankle_dorsiflexion;

    const SIMBICONState simbiconState = getSIMBICONState(s);
    if (simbiconState != UNKNOWN)
    {
        if (simbiconState == STATE0 || simbiconState == STATE1)
        {
            // Left leg is in stance.
            swing_hip_flexion = Biped::hip_r_flexion;
            swing_hip_adduction = Biped::hip_r_adduction;
            swing_knee_extension = Biped::knee_r_extension;
            swing_ankle_dorsiflexion = Biped::ankle_r_dorsiflexion;

            stance_hip_flexion = Biped::hip_l_flexion;
            stance_hip_adduction = Biped::hip_l_adduction;
            stance_knee_extension = Biped::knee_l_extension;
            stance_ankle_dorsiflexion = Biped::ankle_l_dorsiflexion;
        }
        else if (simbiconState == STATE2 || simbiconState == STATE3)
        {
            // Right leg is in stance.
            swing_hip_flexion = Biped::hip_l_flexion;
            swing_hip_adduction = Biped::hip_l_adduction;
            swing_knee_extension = Biped::knee_l_extension;
            swing_ankle_dorsiflexion = Biped::ankle_l_dorsiflexion;

            stance_hip_flexion = Biped::hip_r_flexion;
            // TODO
            stance_hip_adduction = Biped::hip_r_adduction;
            stance_knee_extension = Biped::knee_r_extension;
            stance_ankle_dorsiflexion = Biped::ankle_r_dorsiflexion;
        }

        // Apply swing/stance-dependent PD control to lower limb sagittal coords
        // ---------------------------------------------------------------------
        // simbiconState stateIdx
        // ------------- --------
        // 0             0
        // 1             1
        // 2             0
        // 3             1
        const int stateIdx = simbiconState % 2;
        // Swing.
        coordPDControl(s, swing_hip_flexion, hip_flexion_adduction,
                m_swh[stateIdx], mobForces);
        // TODO
        //coordPDControl(s, swing_hip_adduction, hip_flexion_adduction,
        //        0.0, mobForces);
        coordPDControl(s, swing_knee_extension, knee,
                m_swk[stateIdx], mobForces);
        coordPDControl(s, swing_ankle_dorsiflexion, ankle_flexion,
                m_swa[stateIdx], mobForces);

        // Stance. No hip, since stance hip is used for balance.
        coordPDControl(s, stance_knee_extension, knee,
                m_stk[stateIdx], mobForces);
        coordPDControl(s, stance_ankle_dorsiflexion, ankle_flexion,
                m_sta[stateIdx], mobForces);
    }

    // TODO stance hip adduction shouldn't exist.
}

void SIMBICON::addInBalanceControl(const State& s, Vector& mobForces) const
{
    const SIMBICONState simbiconState = getSIMBICONState(s);

    if (simbiconState == UNKNOWN) return;

    // Which leg is in stance, etc.?
    // =============================
    Biped::Coordinate swing_hip_flexion;
    Biped::Coordinate swing_hip_adduction;
    Biped::Coordinate stance_hip_flexion;
    Biped::Coordinate stance_hip_adduction;
    const MobilizedBody* stanceFoot;
    const MobilizedBody* stanceThigh;
    if (simbiconState == STATE0 || simbiconState == STATE1)
    {
        // Left leg is in stance.
        swing_hip_flexion = Biped::hip_r_flexion;
        swing_hip_adduction = Biped::hip_r_adduction;

        stance_hip_flexion = Biped::hip_l_flexion;
        stance_hip_adduction = Biped::hip_r_adduction;

        stanceFoot = &m_biped.getBody(Biped::foot_l);
        stanceThigh = &m_biped.getBody(Biped::thigh_l);
    }
    else if (simbiconState == STATE2 || simbiconState == STATE3)
    {
        // Right leg is in stance.
        swing_hip_flexion = Biped::hip_l_flexion;
        swing_hip_adduction = Biped::hip_l_adduction;

        stance_hip_flexion = Biped::hip_r_flexion;
        stance_hip_adduction = Biped::hip_r_adduction;

        stanceFoot = &m_biped.getBody(Biped::foot_r);
        stanceThigh = &m_biped.getBody(Biped::thigh_r);
    }
    else throw;

    // Preliminary quantities
    // ======================
    const SimbodyMatterSubsystem& matter = m_biped.getMatterSubsystem();

    // Translation-related quantities.
    // -------------------------------

    // Whole-body mass center, expressed in ground.
    const Vec3 massCenterLoc = matter.calcSystemMassCenterLocationInGround(s);
    // This is 'v' in Yin, 2007.
    const Vec3 massCenterVel = matter.calcSystemMassCenterVelocityInGround(s);

    // Ankle, expressed in ground.
    const Transform& X = stanceFoot->getInboardFrame(s);
    Vec3 stanceAnkleLocInParent = X.p();
    Vec3 stanceAnkleLoc = stanceFoot->getParentMobilizedBody()
        .findStationLocationInGround(s, stanceAnkleLocInParent);

    // Displacement of whole-body COM from ankle, expressed in ground.
    // This is called 'd' in Yin, 2007.
    Vec3 massCenterLocFromStanceAnkle = massCenterLoc - stanceAnkleLoc;

    // Normals to the sagittal and coronal planes.
    UnitVec3 sagittalNormal = m_biped.getNormalToSagittalPlane(s);
    UnitVec3 coronalNormal = m_biped.getNormalToCoronalPlane(s);

    // Projection of 'd' and 'v' into the forward / sagittal and lateral /
    // coronal planes.  We apply balance control separately for these two
    // planes.
    Real massCenterLocFromStanceAnkleForwardMeasure =
        dot(coronalNormal, massCenterLocFromStanceAnkle);
    Real massCenterVelFromStanceAnkleForwardMeasure =
        dot(coronalNormal, massCenterVel);
    Real massCenterLocFromStanceAnkleLateralMeasure =
        dot(sagittalNormal, massCenterLocFromStanceAnkle);
    Real massCenterVelFromStanceAnkleLateralMeasure =
        dot(sagittalNormal, massCenterVel);

    // Rotation-related quantities.
    // ----------------------------
    // Hip flexion and adduction.
    // TODO UnitVec3 stanceThighAxialDir = stanceThigh->getBodyRotation(s).y();
    Vec3 upThigh = swingThigh->getBodyRotation(s).y();
    Vec3 XinSag = cross(UnitVec(YAxis), sagittalNormal);
    Vec3 projUpThigh = upThigh - dot(upThigh, sagittalNormal) * sagittalNormal;
    // TODO points distally proximally?
    Real globalStanceHipFlexionAngle = acos(dot(projUpThigh.normalize(), XinSag)) - Pi/2;
    Real globalStanceHipFlexionRate = stanceThigh->expressGroundVectorInBodyFrame(s,
            stanceThigh->getBodyAngularVelocity(s))[ZAxis];

    // Trunk extension.
    const MobilizedBody* trunk = &m_biped.getBody(Biped::trunk);
    // Directed distally.
    UnitVec3 trunkAxialDir = m_biped.getBody(Biped::trunk).getBodyRotation(s).y();
    UnitVec3 upPelvis = m_biped.getBody(Biped::pelvis).getBodyRotation(s).y();
    Vec3 projUpPelvis = upPelvis - dot(upPelvis, sagittalNormal) * sagittalNormal;
    // UnitVec3(YAxis) gives a vector perpendicular to ground, directed up.
    Real globalTrunkExtensionAngle = acos(dot(projUpPelvis.normalize(), XinSag)) * Pi/2;
    Real globalTrunkExtensionRate = trunk->expressGroundVectorInBodyFrame(s,
            trunk->getBodyAngularVelocity(s))[ZAxis];

    // simbiconState stateIdx
    // ------------- --------
    // 0             0
    // 1             1
    // 2             0
    // 3             1
    const int stateIdx = simbiconState % 2;

    // Apply PD control to swing hip.
    // ==============================
    // \theta_d = \theta_0 + c_d d + c_v v

    // Balance in the sagittal plane.
    /* TODO INCORRECT
     * j
    Real swingHipFlexionTorque =
        coordPDControl(s, swing_hip_flexion, hip_flexion_adduction,
            m_swh[stateIdx] +
            m_cd[stateIdx] * massCenterLocFromStanceAnkleForwardMeasure +
            m_cv[stateIdx] * massCenterVelFromStanceAnkleForwardMeasure,
            mobForces);
    // Balance in the coronal plane.
    // TODO sign.
    Real swingHipAdductionTorque =
        coordPDControl(s, swing_hip_adduction, hip_flexion_adduction,
            0.0 +
            m_cdLat[stateIdx] * massCenterLocFromStanceAnkleLateralMeasure +
            m_cvLat[stateIdx] * massCenterVelFromStanceAnkleLateralMeasure,
            mobForces);
    Real coordPDControl(const State& s,
            const Biped::Coordinate coord, const GainGroup gainGroup,
            const Real thetad, Vector& mobForces) const
    {
        // Prepare quantities.
        Real kp = m_proportionalGains.at(gainGroup);
        Real kd = m_derivativeGains.at(gainGroup);
        Real q = m_biped.getQ(s, coord);
        Real u = m_biped.getU(s, coord);

        // PD control law:
        Real force = clamp(kp * (thetad - q) - kd * u, kp);

        // Update the proper entry in mobForces.
        m_biped.addInForce(coord, force, mobForces);

        return force;
    }
    */
    // Balance in the sagittal plane.
    Real swingHipFlexionTorque;
    Real stanceHipFlexionTorque;
    {
        // Create scope to scrap temporary variables that we use
        // only to make the code clearer.
        Real kp = m_proportionalGains.at(hip_flexion_adduction);
        Real kd = m_derivativeGains.at(hip_flexion_adduction);
        Real thetad = m_swh[stateIdx] +
            m_cd[stateIdx] * massCenterLocFromStanceAnkleForwardMeasure +
            m_cv[stateIdx] * massCenterLocFromStanceAnkleForwardMeasure;
        Real theta = globalHipFlexionAngle;
        Real thetadot =  globalHipFlexionRate;

        swingHipFlexionTorque =
            clamp(kp * (thetad - theta) + kd * thetadot, kp);

        // We use 'trunk' instead of 'torso', which is used in the paper.
        Real trunkPDControl = kp * (0.0 - globalHipFlexionAngle) - kd * globalHipFlexionRate;
        stanceHipFlexionTorque = clamp(trunkPDControl - swingHipFlexionTorque, kp);

    }
    m_biped.addInForce(swing_hip_flexion, swingHipFlexionTorque, mobForces);
    m_biped.addInForce(stance_hip_flexion, stanceHipFlexionTorque, mobForces);

}


/* TODO
   shouldTerminate = false;
   bool lContact;
   bool rContact;
   m_biped.checkContact(s, rContact, lContact);

   if (m_simbicon.getState() == NaN)
   {
   if (rContact) m_simbicon.setState(
   }
   */



// TODO Make SIMBICONSTATEHANDLER not periodic, but triggered upon collisions?
// TODO make biped const within Simbicon?
// TODO hide all the visualization initializtion code.
// TODO make the body of main() really small.
// TODO summary of important pieces of code.
// TODO clean up model building
// TODO initial simbicon state is -1?
