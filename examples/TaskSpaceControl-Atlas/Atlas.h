#ifndef SimTK_SIMBODY_EXAMPLE_ATLAS_H_
#define SimTK_SIMBODY_EXAMPLE_ATLAS_H_

/* -------------------------------------------------------------------------- *
 *             Simbody(tm) Example: Boston Dynamics Atlas robot               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
 * Authors: Michael Sherman                                                   *
 * Contributors: Jack Wang, Chris Dembia, John Hsu                            *
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

#include "Simbody.h"
#include "URDFReader.h"

/* This is a Simbody model of the Boston Dynamics Atlas humanoid robot.

The up direction is Z.      
*/

class Atlas : public SimTK::MultibodySystem {
public:
    Atlas();

    enum Link {
        Ground      = 0,    // These are MobilizedBody indices.
        Base        = 1,
        Shoulder    = 2,
        UpperArm    = 3,
        Forearm     = 4,
        Wrist1      = 5,
        Wrist2      = 6,
        Wrist3      = 7,
        EndEffector = 8
    };
    static const int NumLinks = EndEffector + 1; // include ground

    enum Coords {
        PanCoord    = 0,    // These are the q and u coordinate indices
        LiftCoord   = 1,
        ElbowCoord  = 2,
        Wrist1Coord = 3,
        Wrist2Coord = 4,
        Wrist3Coord = 5
    };
    //static const int NumCoords = Wrist3Coord + 1;

    void setAngleNoise(SimTK::State& state, SimTK::Real qNoise) const 
    {   m_qNoise.setValue(state, qNoise); }

    void setRateNoise(SimTK::State& state, SimTK::Real uNoise) const 
    {   m_uNoise.setValue(state, uNoise); }

    void setSampledAngles(SimTK::State& state, const SimTK::Vector& angles) const {
        m_sampledAngles.setValue(state, angles);
    }

    void setSampledRates(SimTK::State& state, const SimTK::Vector& rates) const {
        m_sampledRates.setValue(state, rates);
    }

    void setSampledEndEffectorPos(SimTK::State& state, 
                                  const SimTK::Vec3& pos) const 
    {   m_sampledEndEffectorPos.setValue(state, pos); }

    SimTK::Real getAngleNoise(const SimTK::State& s) const 
    {   return m_qNoise.getValue(s); }
    SimTK::Real getRateNoise(const SimTK::State& s) const 
    {   return m_uNoise.getValue(s); }
    const SimTK::Vector& getSampledAngles(const SimTK::State& s) const
    {   return m_sampledAngles.getValue(s); }
    const SimTK::Vector& getSampledRates(const SimTK::State& s) const
    {   return m_sampledRates.getValue(s); }
    const SimTK::Vec3& getSampledEndEffectorPos(const SimTK::State& s) const
    {   return m_sampledEndEffectorPos.getValue(s); }

    // Return the Ground frame location of the body origin point of the 
    // EndEffector link.
    SimTK::Vec3 getActualEndEffectorPosition(const SimTK::State& s) const 
    {   return getBody(m_endEffectorLinkName)
                .findStationLocationInGround(s, m_endEffectorStation); }

    // Set a particular joint angle (in radians) in the given State.
    void setJointAngle
       (SimTK::State& s, SimTK::QIndex which, SimTK::Real angle) const 
    {   s.updQ()[which] = angle; }

    // Set a particular joint angular rate (in radians/s) in the given State.
    void setJointRate
       (SimTK::State& s, SimTK::UIndex which, SimTK::Real rate) const 
    {   s.updU()[which] = rate; }

    // Given a set of proposed actuator torques tau, make sure each is in
    // range given the specifications for each actuator.
    static void clampToLimits(SimTK::Vector& tau) {
        SimTK::clampInPlace(-330, tau[PanCoord],   330);
        SimTK::clampInPlace(-330, tau[LiftCoord],  330);
        SimTK::clampInPlace(-150, tau[ElbowCoord], 150);
        SimTK::clampInPlace( -54, tau[Wrist1Coord], 54);
        SimTK::clampInPlace( -54, tau[Wrist2Coord], 54);
        SimTK::clampInPlace( -54, tau[Wrist3Coord], 54);
    }

    const SimTK::SimbodyMatterSubsystem& getMatterSubsystem() const {return m_matter;}
    SimTK::SimbodyMatterSubsystem& updMatterSubsystem() {return m_matter;}

    const SimTK::GeneralForceSubsystem& getForceSubsystem() const {return m_forces;}
    SimTK::GeneralForceSubsystem& updForceSubsystem() {return m_forces;}

    const SimTK::Force::Gravity& getGravity() const {return m_gravity;}

    /// Realizes the topology and model, and uses the resulting initial state
    /// to perform further internal initialization.
    void initialize(SimTK::State& state) {
        state = realizeTopology();
        realizeModel(state);
        m_sampledAngles.setValue(state, SimTK::Vector(state.getNQ(), SimTK::Zero));
        m_sampledRates.setValue(state, SimTK::Vector(state.getNU(), SimTK::Zero));
    }

    const SimTK::MobilizedBody& getBody(const std::string& name) const {
        return m_urdfRobot.links.getLink(name).masterMobod;
    }

    const SimTK::MobilizedBody& getEndEffectorBody() const {
        return m_urdfRobot.links.getLink(m_endEffectorLinkName).masterMobod;
    }

    const Vec3& getEndEffectorStation() const {
        return m_endEffectorStation;
    }

private:
    SimTK::SimbodyMatterSubsystem               m_matter;
    SimTK::GeneralForceSubsystem                m_forces;
    SimTK::ContactTrackerSubsystem              m_tracker;
    SimTK::CompliantContactSubsystem            m_contact;

    SimTK::Force::Gravity                       m_gravity;
    //SimTK::MobilizedBody                        m_bodies[NumLinks];
    SimTK::Measure_<SimTK::Vector>::Variable    m_sampledAngles;
    SimTK::Measure_<SimTK::Vector>::Variable    m_sampledRates;
    SimTK::Measure_<SimTK::Vec3>::Variable      m_sampledEndEffectorPos;
    SimTK::Measure_<SimTK::Real>::Variable      m_qNoise;
    SimTK::Measure_<SimTK::Real>::Variable      m_uNoise;

    SimTK::MultibodyGraphMaker                  m_mbGraph;
    URDFRobot                                   m_urdfRobot;
    std::string                                 m_endEffectorLinkName;
    Vec3                                        m_endEffectorStation;
};



#endif // SimTK_SIMBODY_EXAMPLE_ATLAS_H_
