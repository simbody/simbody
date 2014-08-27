#ifndef SimTK_SIMBODY_HUMANOID_H_
#define SimTK_SIMBODY_HUMANOID_H_

/* -------------------------------------------------------------------------- *
 *                        Simbody(tm) Example: Humanoid                       *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-14 Stanford University and the Authors.        *
 * Authors: Chris Dembia, Jack Wang                                           *
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

#include <Simbody.h>

using namespace std;
using namespace SimTK;

#define DEG(rad) convertRadiansToDegrees(rad)
#define RAD(deg) convertDegreesToRadians(deg)

/// The MultibodySystem that we will control with the SIMBICON controller
/// This system is independent of the controller, and could conceivably
/// be controlled by a different controller.
/// x is directed anteriorly (forward, the direction of walking).
/// y is up.
/// z is to the Humanoid's right.
class Humanoid : public MultibodySystem {
public:
    Humanoid();

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

        // Initialization specific to Humanoid.
        fillInCoordinateMap(state);
    }

    /// Populates a map that allow users of Humanoid to access coordinates (Q)
    /// and  speeds (U) with descriptive Enums (e.g., neck_extension) so the
    /// user (e.g., SIMBICON) needn't know how State indices correspond to 
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

    /// The names of the rigid bodies that make up the biped.
    /// Used for indexing into data structures.
    enum Segment {
        trunk,
        head,
        pelvis,
        upperarm_r, lowerarm_r, hand_r, upperarm_l, lowerarm_l, hand_l,
        thigh_r, shank_r, foot_r, toes_r,
        thigh_l, shank_l, foot_l, toes_l
    };

    /// The names of the degrees of freedom of the biped.
    /// Used for indexing into data structures.
    enum Coordinate {
        neck_extension, neck_bending, neck_rotation,
        back_tilt, back_list, back_rotation,
        shoulder_r_flexion, shoulder_r_adduction, shoulder_r_rotation,
        elbow_r_flexion, elbow_r_rotation,
        shoulder_l_flexion, shoulder_l_adduction, shoulder_l_rotation,
        elbow_l_flexion, elbow_l_rotation,
        hip_r_adduction, hip_r_flexion, hip_r_rotation,
        knee_r_extension,
        ankle_r_inversion, ankle_r_dorsiflexion,
        mtp_r_dorsiflexion,
        hip_l_adduction, hip_l_flexion, hip_l_rotation,
        knee_l_extension,
        ankle_l_inversion, ankle_l_dorsiflexion,
        mtp_l_dorsiflexion,
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

    /// The pelvis' z direction (which points to the biped's right) projected
    /// into the plane of the ground, and normalized into a unit vector.
    UnitVec3 getNormalToSagittalPlane(const State& s) const
    {
        Vec3 bodyZInGround = Vec3(m_bodies.at(pelvis).getBodyRotation(s).z());
        // Project onto a horizontal plane.
        bodyZInGround[YAxis] = 0;
        // Make into a unit vector.
        return UnitVec3(bodyZInGround);
    }

    /// The pelvis' x direction (which points in the direction of walking)
    /// projected into the plane of the ground, and normalized into a
    /// unit vector.
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

    // Index of the Gravity Force in m_forces.
    ForceIndex                   m_gravityIndex;

    /// The MobilizedBody's that make up the Humanoid.
    std::map<Segment, MobilizedBody> m_bodies;

    /// The indices of each coordinate's Q and U in the State vector.
    std::map<Coordinate, std::pair<QIndex, UIndex> > m_coordinates;

    #ifdef RIGID_CONTACT
        std::vector<PointPlaneContact*> m_rightContacts;
        std::vector<PointPlaneContact*> m_leftContacts;
    #endif
};

Humanoid::Humanoid()
    : m_matter(*this), m_forces(*this),
      m_tracker(*this), m_contact(*this, m_tracker)
{

    m_matter.setShowDefaultGeometry(false);

    //--------------------------------------------------------------------------
    //                          Constants, etc.
    //--------------------------------------------------------------------------
    // Contact parameters are elsewhere.

    // Properties of joint stops.
    const Real stopStiffness = DEG(1000); // N-m/deg -> N-m/radian
    const Real stopDissipation = 1;

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

    // Gravity.
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
    //                  Contact at feet and toes: Compliant
    //--------------------------------------------------------------------------
#ifndef RIGID_CONTACT

    // Parameters.
    // -----------
    // slide -> stick velocity.
    // TODO unused. const Real transitionVelocity = .03; 

    // Contact spheres attached to feet.
    const Real contactSphereRadius = 1.5 * .01;

    // Friction coefficients.
    const Real mu_s = 5;
    const Real mu_d = 5;
    const Real mu_v = 0;

    // Ground.
    // -------
    
    // Contact material: concrete.
    const Real concrete_density = 2300.;  // kg/m^3
    const Real concrete_young   = 25e9;  // pascals (N/m)
    const Real concrete_poisson = 0.15;    // ratio
    const Real concrete_planestrain =
        ContactMaterial::calcPlaneStrainStiffness(concrete_young,
                concrete_poisson);
    const Real concrete_dissipation = 0.005;
    
    const ContactMaterial concrete(concrete_planestrain,concrete_dissipation,
                                   mu_s,mu_d,mu_v);

    m_matter.updGround().updBody().addContactSurface(
            Transform(Rotation(-0.5 * Pi, ZAxis), Vec3(0)),
            ContactSurface(ContactGeometry::HalfSpace(), concrete));

    // Feet.
    // -----

    // Contact material: rubber.
    const Real rubber_density = 1100.;  // kg/m^3
    const Real rubber_young   = 0.01e9; // pascals (N/m)
    const Real rubber_poisson = 0.5;    // ratio
    const Real rubber_planestrain =
        ContactMaterial::calcPlaneStrainStiffness(rubber_young,rubber_poisson);
    const Real rubber_dissipation = /*0.005*/1;
    
    const ContactMaterial rubber(rubber_planestrain,rubber_dissipation,
                                   mu_s,mu_d,mu_v);

    ContactSurface contactBall(ContactGeometry::Sphere(contactSphereRadius),
                               rubber);
    
    // Use this clique for contact surfaces on the humanoid that you don't want
    // to have bump into one another. They will still contact Ground.
    const ContactCliqueId clique = ContactSurface::createNewContactClique();
    contactBall.joinClique(clique);

    DecorativeSphere contactArt(contactSphereRadius);
    contactArt.setColor(Magenta);

    const Vec2 rlatmed(.04,-.04);
    const Vec2 heelball(-.02,.125);

    // Feet.
    const Real footsphereht = -.07;
    for (int x = 0; x <= 1; ++x)
    {
        for (int z = 0; z <= 1; ++z)
        {
            const Vec3 rctr(heelball[x], footsphereht,  rlatmed[z]);
            const Vec3 lctr(heelball[x], footsphereht, -rlatmed[z]);

            foot_rInfo.addContactSurface(rctr, contactBall);
            foot_rInfo.addDecoration(rctr, contactArt);
            foot_lInfo.addContactSurface(lctr, contactBall);
            foot_lInfo.addDecoration(lctr, contactArt);
        }
    }

    // Balls just at toe tips.
    const Real toetip = .06;
    for (int z = 0; z <= 1; ++z)
    {
        const Vec3 rctr(toetip, 0,  rlatmed[z]);
        const Vec3 lctr(toetip, 0, -rlatmed[z]);

        toes_rInfo.addContactSurface(rctr, contactBall);
        toes_rInfo.addDecoration(rctr, contactArt);
        toes_lInfo.addContactSurface(lctr, contactBall);
        toes_lInfo.addDecoration(lctr, contactArt);
    }
#endif

    //--------------------------------------------------------------------------
    //                         Mobilized Bodies
    //--------------------------------------------------------------------------
    // Trunk.
    // ------
    m_bodies[trunk] = MobilizedBody::Free(
        m_matter.updGround(), Vec3(0, 1.5, 0),
        trunkInfo,          Vec3(0));

    // Neck angles are: q0=extension (about z), q1=bending (x), q2=rotation (y).
    m_bodies[head] = MobilizedBody::Ball(
        m_bodies[trunk], Transform(Rzxy, Vec3(0.010143822053248,0.222711680750785,0)),
        headInfo, Rzxy);
    addMobilityLinearStop(m_bodies[head], 0, -80, 50); // extension
    addMobilityLinearStop(m_bodies[head], 1, -60, 50); // bending
    addMobilityLinearStop(m_bodies[head], 2, -80, 80); // rotation

    // Back angles are: q0=tilt (about z), q1=list (x), q2=rotation (y).
    m_bodies[pelvis] = MobilizedBody::Ball(
        m_bodies[trunk], Transform(Rzxy, Vec3(-0.019360589663647,-0.220484136504589,0)),
        pelvisInfo, Rzxy);
    addMobilityLinearStop(m_bodies[pelvis], 0, -5, 10); // tilt
    addMobilityLinearStop(m_bodies[pelvis], 1, -5, 5); // list
    addMobilityLinearStop(m_bodies[pelvis], 2, -15, 15); // rotation

    // Right arm.
    //-----------
    // Shoulder angles are q0=flexion, q1=adduction, q2=rotation
    m_bodies[upperarm_r] = MobilizedBody::Ball(
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
    m_bodies[upperarm_l] = MobilizedBody::Ball(
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
    m_bodies[thigh_r] = MobilizedBody::Ball(
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
    m_bodies[thigh_l] = MobilizedBody::Ball(
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

    //--------------------------------------------------------------------------
    //                    Contact at feet and toes: Rigid
    //--------------------------------------------------------------------------

#ifdef RIGID_CONTACT
    // Friction coefficients.
    const Real CoefRest = 0.0;
    const Real mu_s = 0.8;
    const Real mu_d = 0.5;
    const Real mu_v = 0;

    // Where the points are located.
    const Vec2 rlatmed(.04,-.04);
    const Vec2 heelball(-.02,.125);

    // Feet.
    const Real footsphereht = -.07;
    for (int x = 0; x <= 1; ++x)
    {
        for (int z = 0; z <= 1; ++z)
        {
            const Vec3 rctr(heelball[x], footsphereht,  rlatmed[z]);
            const Vec3 lctr(heelball[x], footsphereht, -rlatmed[z]);

            PointPlaneContact* ppcr = new PointPlaneContact(
                    m_matter.updGround(), YAxis, 0.,
                    m_bodies[foot_r], rctr, CoefRest, mu_s, mu_d, mu_v);
            PointPlaneContact* ppcl = new PointPlaneContact(
                    m_matter.updGround(), YAxis, 0.,
                    m_bodies[foot_l], lctr, CoefRest, mu_s, mu_d, mu_v);

            m_matter.adoptUnilateralContact(ppcr);
            m_matter.adoptUnilateralContact(ppcl);

            m_rightContacts.push_back(ppcr);
            m_leftContacts.push_back(ppcl);
        }
    }

    // Balls just at toe tips.
    const Real toetip = .06;
    for (int z = 0; z <= 1; ++z)
    {
        const Vec3 rctr(toetip, 0,  rlatmed[z]);
        const Vec3 lctr(toetip, 0, -rlatmed[z]);

        PointPlaneContact * ppcr = new PointPlaneContact(
                m_matter.updGround(), YAxis, 0.,
                m_bodies[toes_r], rctr, CoefRest, mu_s, mu_d, mu_v);
        PointPlaneContact * ppcl = new PointPlaneContact(
                m_matter.updGround(), YAxis, 0.,
                m_bodies[toes_l], lctr, CoefRest, mu_s, mu_d, mu_v);

        m_matter.adoptUnilateralContact(ppcr);
        m_matter.adoptUnilateralContact(ppcl);

        m_rightContacts.push_back(ppcr);
        m_leftContacts.push_back(ppcl);
    }
#endif
}

void Humanoid::fillInCoordinateMap(const State& s)
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

void Humanoid::findContactForces(const State& s, Real& fLeft, Real& fRight)
    const
{
    fLeft = 0;
    fRight = 0;

    #ifdef RIGID_CONTACT

        for (unsigned int iRight = 0; iRight < m_rightContacts.size(); ++iRight)
        {
            const UnilateralContact* contact = m_rightContacts[iRight];
            fRight += s.getMultipliers()[contact->getContactMultiplierIndex(s)];
        }
        for (unsigned int iLeft = 0; iLeft < m_leftContacts.size(); ++iLeft)
        {
            const UnilateralContact* contact = m_leftContacts[iLeft];
            fLeft += s.getMultipliers()[contact->getContactMultiplierIndex(s)];
        }
        
    #else

        const unsigned int nContacts = m_contact.getNumContactForces(s);
        const ContactSnapshot& snapshot = m_tracker.getActiveContacts(s);
    
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

    #endif
}

void Humanoid::findContactStatus(const State& s, bool& left, bool& right) const
{
    #ifdef RIGID_CONTACT
        left = false;
        right = false;
        for (unsigned int iRight = 0; iRight < m_rightContacts.size(); ++iRight)
        {
            right = right || m_rightContacts[iRight]->isEnabled(s);
        }
        for (unsigned int iLeft = 0; iLeft < m_leftContacts.size(); ++iLeft)
        {
            left = left || m_leftContacts[iLeft]->isEnabled(s);
        }
    #else
        Real fLeft;
        Real fRight;
        findContactForces(s, fLeft, fRight);
        left = fLeft > 0;
        right = fRight > 0;
    #endif
}

#endif // SimTK_SIMBODY_BIPEDSYSTEM_H_
