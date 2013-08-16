#ifndef SimTK_SIMMATH_RUNGE_KUTTA_FELDBERG_INTEGRATOR_H_
#define SimTK_SIMMATH_RUNGE_KUTTA_FELDBERG_INTEGRATOR_H_

/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman, Peter Eastman                                    *
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
#include "simmath/internal/common.h"
#include "simmath/Integrator.h"

namespace SimTK {

/**
 * This is an Integrator based on the Runge-Kutta-Feldberg algorithm.  It is an
 * error controlled, fifth order explicit integrator.
 */

class RungeKuttaFeldbergIntegratorRep;

class SimTK_SIMMATH_EXPORT RungeKuttaFeldbergIntegrator : public Integrator {
public:
    explicit RungeKuttaFeldbergIntegrator(const System& sys);
};

} // namespace SimTK

#endif // SimTK_SIMMATH_RUNGE_KUTTA_FELDBERG_INTEGRATOR_H_


