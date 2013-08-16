#ifndef SimTK_SIMBODY_LOCAL_ENERGY_MINIMIZER_H_
#define SimTK_SIMBODY_LOCAL_ENERGY_MINIMIZER_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
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
#include "simbody/internal/common.h"
#include "simbody/internal/MultibodySystem.h"

namespace SimTK {

/**
 * This class performs local potential energy minimization of a MultibodySystem.
 * Only positions (generalized coordinates q) are changed; velocities are ignored.
 */

class SimTK_SIMBODY_EXPORT LocalEnergyMinimizer {
public:
    /**
     * Find the local potential energy minimum of a MultibodySystem.
     * 
     * @param system       the system whose energy should be minimized
     * @param state        on entry, this should contain the starting state from which to begin the search.
     *                     On exit, it contains a (usually nearby) state with revised generalized coordinate
     *                     values q at which the energy is a local minimum.
     * @param tolerance    a tolerance for the energy minimization.  The search ends when no component of the
     *                     energy gradient is larger than this.
     */
    static void minimizeEnergy(const MultibodySystem& system, State& state, Real tolerance);
private:
    class OptimizerFunction;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_LOCAL_ENERGY_MINIMIZER_H_
