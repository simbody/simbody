#ifndef SimTK_SIMMATH_RUNGE_KUTTA_3_INTEGRATOR_H_
#define SimTK_SIMMATH_RUNGE_KUTTA_3_INTEGRATOR_H_

/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
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
class RungeKutta3IntegratorRep;

/**
 * This is a 3rd order Runge-Kutta Integrator using coefficents from J.C. Butcher's
 * book "The Numerical Analysis of Ordinary Differential Equations", John Wiley & Sons,
 * 1987, page 325.  It is an error controlled, third order, three stage explicit integrator
 * with an embedded 2nd order error estimate.
 */
class SimTK_SIMMATH_EXPORT RungeKutta3Integrator : public Integrator {
public:
    explicit RungeKutta3Integrator(const System& sys);
};

} // namespace SimTK

#endif // SimTK_SIMMATH_RUNGE_KUTTA_3_INTEGRATOR_H_


