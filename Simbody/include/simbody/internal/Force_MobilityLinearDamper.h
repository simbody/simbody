#ifndef SimTK_SIMBODY_FORCE_MOBILITY_LINEAR_DAMPER_H_
#define SimTK_SIMBODY_FORCE_MOBILITY_LINEAR_DAMPER_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-13 Stanford University and the Authors.        *
 * Authors: Peter Eastman, Michael Sherman                                    *
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

#include "SimTKcommon.h"
#include "simbody/internal/Force.h"

/** @file
This contains the user-visible API ("handle" class) for the SimTK::Force 
subclass Force::MobilityLinearDamper and is logically part of Force.h. The file
assumes that Force.h will have included all necessary declarations. **/

namespace SimTK {

/**
 * A linear damper on a mobility coordinate. The
 * damping constant c is provided, with the generated force
 * being -c*u where u is the mobility's generalize speed.
 * This is meaningful on any mobility, since all our
 * generalized speeds have physical meaning. This is not
 * a potential force and hence does not contribute to
 * potential energy.
 */

class SimTK_SIMBODY_EXPORT Force::MobilityLinearDamper : public Force {
public:
    /**
     * Create a %MobilityLinearDamper.
     * 
     * @param forces     the subsystem to which this force should be added
     * @param body       the body to which the force should be applied
     * @param coordinate the index of the coordinate in the body's u vector to which the force should be applied
     * @param damping    the damping constant
     */
    MobilityLinearDamper(GeneralForceSubsystem& forces, const MobilizedBody& body, int coordinate, Real damping);
    
    /** Default constructor creates an empty handle. **/
    MobilityLinearDamper() {}

    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(MobilityLinearDamper, MobilityLinearDamperImpl, Force);
};

} // namespace SimTK

#endif // SimTK_SIMBODY_FORCE_MOBILITY_LINEAR_DAMPER_H_
