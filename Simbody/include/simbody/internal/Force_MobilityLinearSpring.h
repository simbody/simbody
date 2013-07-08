#ifndef SimTK_SIMBODY_FORCE_MOBILITY_LINEAR_SPRING_H_
#define SimTK_SIMBODY_FORCE_MOBILITY_LINEAR_SPRING_H_

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
subclass Force::MobilityLinearSpring and is logically part of Force.h. The file
assumes that Force.h will have included all necessary declarations. **/

namespace SimTK {

/**
 * A linear spring along or around a mobility coordinate. The
 * stiffness k is provided, along with an arbitrary "zero" 
 * coordinate value q0 at which the spring generates no force.
 * The generated force is k*(q-q0), and potential energy is 
 * pe = 1/2 k (q-q0)^2.
 * This is not meaningful unless the mobility coordinate is such that qdot=u 
 * for that coordinate.  In particular, do not use this on a coordinate
 * which is part of a quaternion.
 */

class SimTK_SIMBODY_EXPORT Force::MobilityLinearSpring : public Force {
public:
    /**
     * Create a %MobilityLinearSpring.
     * 
     * @param forces     the subsystem to which this force should be added
     * @param body       the body to which the force should be applied
     * @param coordinate the index of the coordinate in the body's u vector to which the force should be applied
     * @param k          the spring constant
     * @param q0         the value of the coordinate at which the force is 0
     */
    MobilityLinearSpring(GeneralForceSubsystem& forces, const MobilizedBody& body, int coordinate, Real k, Real q0);
    
    /** Default constructor creates an empty handle. **/
    MobilityLinearSpring() {}

    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(MobilityLinearSpring, MobilityLinearSpringImpl, Force);
};

} // namespace SimTK

#endif // SimTK_SIMBODY_FORCE_MOBILITY_LINEAR_SPRING_H_
