#ifndef SimTK_SIMMATH_CPODES_INTEGRATOR_H_
#define SimTK_SIMMATH_CPODES_INTEGRATOR_H_

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
#include "simmath/internal/SimTKcpodes.h"

namespace SimTK {

class CPodesIntegratorRep;

/**
 * This is an Integrator based on the CPODES library.  It is an error controlled, variable order
 * implicit integrator.  It provides a good combination of accuracy, stability, and speed, and
 * is a good choice for integrating stiff problems.
 * 
 * When creating a CPodesIntegrator, you can specify various options for how to perform
 * the implicit integration: the linear multistep method to use (Adams or BDF), and the nonlinear system
 * iteration type (Newton iteration or functional iteration).  For stiff problems, the recommended choices
 * are BDF with Newton iteration.  For non-stiff problems, using Adams and/or functional iteration may
 * provide better performance.  Note that Adams is <i>never</i> recommended for systems that include
 * constraints.
 */

class SimTK_SIMMATH_EXPORT CPodesIntegrator : public Integrator {
public:
    /**
     * Create a CPodesIntegrator for integrating a System.
     */
    explicit CPodesIntegrator(const System& sys, CPodes::LinearMultistepMethod method=CPodes::BDF);
    /**
     * Create a CPodesIntegrator for integrating a System.  The nonlinear system iteration type is chosen automatically
     * based on the linear multistep method: Newton iteration for BDF (the default), and functional iteration for Adams.
     */
    CPodesIntegrator(const System& sys, CPodes::LinearMultistepMethod method, CPodes::NonlinearSystemIterationType iterationType);
    /**
     * CPODES provides its own mechanism for projecting the system onto the constraint manifold.  By default,
     * CPodesIntegrator uses the System's project() method for doing projection, which is usually more
     * efficient.  Invoking this method tells it to use the CPODES mechanism instead.
     * 
     * This method must be invoked before the integrator is initialized.  Invoking it after initialization
     * will produce an exception.
     */
    void setUseCPodesProjection();
    /**
     * Restrict the integrator to lower orders than it is otherwise capable of (up to 12 for Adams, 5 for BDF).  This method
     * may only be used to decrease the maximum order permitted, never to increase it.  Once you specify an order limit, calling it
     * again with a larger value will fail.
     */
    void setOrderLimit(int order);
};

} // namespace SimTK

#endif // SimTK_SIMMATH_CPODES_INTEGRATOR_H_
