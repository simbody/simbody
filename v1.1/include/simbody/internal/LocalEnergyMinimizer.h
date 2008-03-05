#ifndef SimTK_SIMBODY_LOCAL_ENERGY_MINIMIZER_H_
#define SimTK_SIMBODY_LOCAL_ENERGY_MINIMIZER_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/MultibodySystem.h"

namespace SimTK {

/**
 * This class performs local energy minimization of a MultibodySystem.
 */

class SimTK_SIMBODY_EXPORT LocalEnergyMinimizer {
public:
    /**
     * Find the local energy minimum of a MultibodySystem.
     * 
     * @param system       the system whose energy should be minimized
     * @param state        on entry, this should contain the starting state from which to begin the search.
     *                     On exit, it contains a (usually nearby) state at which the energy is a local minimum.
     * @param tolerance    a tolerance for the energy minimization.  The search ends when no component of the
     *                     energy gradient is larger than this.
     */
    static void minimizeEnergy(const MultibodySystem& system, State& state, Real tolerance);
private:
    class OptimizerFunction;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_LOCAL_ENERGY_MINIMIZER_H_
