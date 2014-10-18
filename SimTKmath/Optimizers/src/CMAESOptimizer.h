#ifndef _SimTK_CMAES_OPTIMIZER_H_
#define _SimTK_CMAES_OPTIMIZER_H_

/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-14 Stanford University and the Authors.        *
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

#include "simmath/internal/common.h"
#include "simmath/internal/OptimizerRep.h"
#include "c-cmaes/cmaes_interface.h"

namespace SimTK {

class CMAESOptimizer: public Optimizer::OptimizerRep {
public:
    CMAESOptimizer(const OptimizerSystem& sys);
    OptimizerRep* clone() const OVERRIDE_11;
    Real optimize(SimTK::Vector& results) OVERRIDE_11;
    OptimizerAlgorithm getAlgorithm() const OVERRIDE_11 { return CMAES; }

private:
    // Bridges Simbody's interface for options with c-cmaes.
    void processSettingsBeforeCMAESInit(cmaes_t& evo);
    void processSettingsAfterCMAESInit(cmaes_t& evo);

    // The seed for cmaes, which we provide to cmaes_init.
    // In cmaes this is an unsigned int. But, Simbody users are able to provide
    // an int. So, we check to see if the int they provide is nonnegative.
    // If we made this variable an unsigned int, we might let a user error
    // silently pass.
    int m_seed;

};

} // namespace SimTK
#endif // _SimTK_CMAES_OPTIMIZER_H_
