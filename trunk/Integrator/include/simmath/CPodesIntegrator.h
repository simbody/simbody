#ifndef SimTK_SIMMATH_CPODES_INTEGRATOR_H_
#define SimTK_SIMMATH_CPODES_INTEGRATOR_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmath(tm)                         *
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
#include "SimTKcpodes.h"

#include "simmath/internal/common.h"
#include "simmath/Integrator.h"

namespace SimTK {

/**
 * This is an integrator based on the CPODES library.  It is an error controlled,
 * variable order implicit integrator.
 * 
 * When creating a CPodesIntegrator, you can specify various options for how to perform
 * the implicit integration: the linear multistep method to use, and the nonlinear system
 * iteration type.  For stiff problems, the recommended choices are BDF with Newton iteration.
 * For non-stiff problems, the recommended choices are Adams with functional iteration.
 */

class CPodesIntegratorRep;

class SimTK_SIMMATH_EXPORT CPodesIntegrator : public Integrator {
public:
    /**
     * Create a CPodesIntegrator for integrating a System.
     */
    CPodesIntegrator(const System& sys, CPodes::LinearMultistepMethod method=CPodes::BDF);
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
     * Restrict the integrator to a lower order than it is otherwise capable of (12 for Adams, 5 for BDF).  This method
     * may only be used to decrease the maximum order permitted, never to increase it.  Once you specify an order limit, calling it
     * again with a larger value will fail.
     */
    void setOrderLimit(int order);
};

} // namespace SimTK

#endif // SimTK_SIMMATH_CPODES_INTEGRATOR_H_
