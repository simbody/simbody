/* -------------------------------------------------------------------------- *
 *                 Simbody(tm) Example: SIMBICON                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-12 Stanford University and the Authors.        *
 * Authors: Jack Wang, Tim Dorn                                               *
 * Contributors: Chris Dembia                                                 *
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

/* SIMBICON stands for Simple Biped Controller. It is described at
https://www.cs.ubc.ca/~van/papers/Simbicon.htm.

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
    * velocity balance feedback coefficient
    * torso target angle
    * swing-hip target angle
    * swing-knee target angle
    * swing-ankle target angle
    * stance-knee target angle
    * stance-ankle target angle

See [1], particularly Table 1 of [1], for more information.

I (Chris Dembia) wrote this code pretty much by copying code that Michael
Sherman gave me. The code that Michael Sherman gave me was originally written
by Jack Wang and Tim Dorn.

Yin, KangKang, Kevin Loken, and Michiel van de Panne. "SIMBICON: Simple biped
locomotion control." ACM Transactions on Graphics (TOG). Vol. 26. No. 3. ACM,
2007.

*/

#include "Simbody.h"

using namespace SimTK;

//==============================================================================
// BIPED
//==============================================================================
/// The MultibodySystem that we will control with the SIMBICON controller
/// (though the system is independent of the controller, and coulc conceivably
/// be controlled by a different controller).
class Biped : public MultibodySystem {
public:
    Biped() 
        : m_matter(*this), m_forces(*this) {}

    const SimbodyMatterSubsystem& getMatterSubsystem() const {return m_matter;}
    SimbodyMatterSubsystem& updMatterSubsystem() {return m_matter;}

    const GeneralForceSubsystem& getForceSubsystem() const {return m_forces;}
    GeneralForceSubsystem& updForceSubsystem() {return m_forces;}

    SimbodyMatterSubsystem       m_matter;
    GeneralForceSubsystem        m_forces;
};

//==============================================================================
// SIMBICON
//==============================================================================
/// The actual controller that specifies the torques to apply to the Biped.
class SIMBICON : public Force::Custom::Implementation {
public:
    SIMBICON(const Biped& biped)
        : m_biped(biped) {}
    void calcForce(const State&         state, 
                   Vector_<SpatialVec>& bodyForces, 
                   Vector_<Vec3>&       particleForces, 
                   Vector&              mobilityForces) const 
                   OVERRIDE_11
    {
        // TODO
    }

    Real calcPotentialEnergy(const State& state) const OVERRIDE_11
    {
        return 0;
    }
private:
    const Biped& m_biped;
};

//==============================================================================
// OUTPUT REPORTER
//==============================================================================
/// This is a periodic event handler that we use to print information to the
/// screen during the simulation.
class OutputReporter : public PeriodicEventReporter {
public:
    OutputReporter(const Biped& biped, Real interval)
        : PeriodicEventReporter(interval), m_biped(biped) {}

    void handleEvent(const State& state) const OVERRIDE_11
    {
        // TODO
    }
private:
    const Biped& m_biped;
};

//==============================================================================
// MAIN FUNCTION
//==============================================================================
int main(int argc, char **argv)
{
    // Create system.
    Biped biped;
    biped.addEventReporter(new OutputReporter(biped, 0.01));

    // Add controller to the force system. See Doxygen for Force::Custom.
    // The name of this Force is irrelevant.
    Force::Custom simbicon1(biped.updForceSubsystem(), new SIMBICON(biped));
	return 0; 
}












