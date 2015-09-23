#ifndef SimTK_SIMMATH_VERLET_INTEGRATOR_H_
#define SimTK_SIMMATH_VERLET_INTEGRATOR_H_

/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
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
#include "simmath/internal/common.h"
#include "simmath/Integrator.h"

namespace SimTK {

class VerletIntegratorRep;

/**
 * This is an Integrator based on the velocity Verlet algorithm.  It is a third order, semi-explicit integrator.
 * Velocity independent forces only need to be evaluated once per time step, but forces which depend on velocity
 * must be evaluated multiple times.  This makes it a very efficient algorithm for systems where most of the time
 * is spent evaluating forces that depend only on position.
 *
 * Although this is a third order integrator, the velocities reported at each time step are only accurate to lower
 * order.  If any forces depend on velocity, this may lead to a reduction in the overall accuracy of integration,
 * since the forces being integrated are less accurate than the method used to integrate them.  Whether this
 * actually happens in a particular case depends on the magnitude of the velocity dependent forces, how sensitive
 * they are to errors in velocity, and how the system is affected by those forces.
 *
 * When this integrator is used with fixed size time steps, it is symplectic.  This means that it is extremely
 * good at conserving energy, and will usually produce much less variation in energy than most other integrators
 * would at the same step size.  This makes it a good choice for problems where accurate energy conservation over
 * long time periods is important.  To use it in this way, use the constructor which takes a step size, or call
 * setFixedStepSize().
 *
 * Alternatively, it may be used in variable step mode, in which case the step size is selected based on the
 * accuracy specified by calling setAccuracy().  In this case it is no longer symplectic, so the energy will
 * fluctuate more over time.  Because the step size is adjusted based on the local error in each step, however,
 * the trajectory will generally be more accurate in variable step size mode than would be obtained with the
 * same number of fixed size steps spanning the same amount of time.
 *
 * Another possible strategy is to set only a maximum step size but not a minimum step size.  It should be chosen
 * such that most time steps will use the fixed maximum step size, but smaller steps can be used when necessary to
 * preserve accuracy.  This leads to fairly good long term energy conservation, while still maintaining reasonable
 * accuracy when unusually large forces transiently occur.
 */

class SimTK_SIMMATH_EXPORT VerletIntegrator : public Integrator {
public:
    /**
     * Create a VerletIntegrator for integrating a System with variable size steps.
     */
    explicit VerletIntegrator(const System& sys);
    /**
     * Create a VerletIntegrator for integrating a System with fixed size steps.
     */
    VerletIntegrator(const System& sys, Real stepSize);
};

} // namespace SimTK

#endif // SimTK_SIMMATH_VERLET_INTEGRATOR_H_


